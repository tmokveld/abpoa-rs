//! AbPOA aligner
//!
//! Module tries to keep the unsafe FFI contained while exposing methods for
//! running alignments, updating graphs, and exporting consensus/MSA/GFA outputs

use crate::encode::{
    Alphabet, decode_aa, decode_aa_code, decode_dna, decode_dna_code, encode_aa, encode_dna,
};
use crate::graph::Graph;
pub use crate::output::{
    dot::{EdgeLabel, EdgePenWidth, PogDotOptions, RankDir},
    pog::PogWriteOptions,
};
use crate::params::{NodeId, OutputMode, Parameters, SentinelNode};
use crate::result::{EncodedMsaResult, MsaResult};
use crate::{Error, Result, sys};
use libc;
use std::{
    env,
    ffi::CString,
    fs::{self, OpenOptions},
    io::Write,
    marker::PhantomData,
    os::raw::c_char,
    path::{Path, PathBuf},
    process, ptr,
    ptr::NonNull,
    rc::Rc,
    slice,
    time::{SystemTime, UNIX_EPOCH},
};

#[cfg(unix)]
use std::os::unix::io::AsRawFd;

/// Node ids (exclusive) delimiting a subgraph to align against
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubgraphRange {
    pub beg: NodeId,
    pub end: NodeId,
}

impl SubgraphRange {
    fn as_raw(self) -> (i32, i32) {
        (self.beg.0, self.end.0)
    }
}

/// CIGAR operations emitted by abPOA
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match,
    Insertion,
    Deletion,
    Mismatch,
    SoftClip,
    HardClip,
}

impl CigarOp {
    fn from_raw(value: u64) -> Self {
        match value & 0xf {
            v if v == sys::ABPOA_CMATCH as u64 => CigarOp::Match,
            v if v == sys::ABPOA_CINS as u64 => CigarOp::Insertion,
            v if v == sys::ABPOA_CDEL as u64 => CigarOp::Deletion,
            v if v == sys::ABPOA_CDIFF as u64 => CigarOp::Mismatch,
            v if v == sys::ABPOA_CSOFT_CLIP as u64 => CigarOp::SoftClip,
            v if v == sys::ABPOA_CHARD_CLIP as u64 => CigarOp::HardClip,
            _ => unreachable!("abPOA should only emit the defined cigar opcodes"),
        }
    }

    /// Single-letter CIGAR code for this operation
    pub fn as_char(self) -> char {
        match self {
            CigarOp::Match => 'M',
            CigarOp::Insertion => 'I',
            CigarOp::Deletion => 'D',
            CigarOp::Mismatch => 'X',
            CigarOp::SoftClip => 'S',
            CigarOp::HardClip => 'H',
        }
    }
}

/// Decoded representation of a single abPOA graph CIGAR entry
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphCigarOp {
    /// Graph node aligned to a query position (match or mismatch)
    Aligned {
        op: CigarOp,
        node_id: NodeId,
        query_index: usize,
    },
    /// Graph node skipped from the query
    Deletion { node_id: NodeId, len: u32 },
    /// Query-only run (insertion or clipping) ending at `end` with length `len`
    QueryRun { op: CigarOp, end: usize, len: u32 },
}

impl GraphCigarOp {
    /// Length in bases that this operation consumes on the query
    pub fn query_len(&self) -> usize {
        match self {
            GraphCigarOp::Aligned { .. } => 1,
            GraphCigarOp::Deletion { .. } => 0,
            GraphCigarOp::QueryRun { len, .. } => (*len).max(1) as usize,
        }
    }

    /// CIGAR opcode for this entry
    pub fn op(&self) -> CigarOp {
        match self {
            GraphCigarOp::Aligned { op, .. } => *op,
            GraphCigarOp::Deletion { .. } => CigarOp::Deletion,
            GraphCigarOp::QueryRun { op, .. } => *op,
        }
    }

    /// Query span covered by this operation, if any
    pub fn query_range(&self) -> Option<std::ops::RangeInclusive<usize>> {
        match self {
            GraphCigarOp::Aligned { query_index, .. } => Some(*query_index..=*query_index),
            GraphCigarOp::Deletion { .. } => None,
            GraphCigarOp::QueryRun { end, len, .. } => {
                let len = (*len).max(1) as usize;
                let start = end.saturating_add(1).saturating_sub(len);
                Some(start..=*end)
            }
        }
    }
}

/// Iterator over decoded graph CIGAR entries
pub struct GraphCigar<'a> {
    entries: std::slice::Iter<'a, u64>,
}

impl<'a> GraphCigar<'a> {
    fn new(entries: &'a [u64]) -> Self {
        Self {
            entries: entries.iter(),
        }
    }
}

impl<'a> Iterator for GraphCigar<'a> {
    type Item = GraphCigarOp;

    fn next(&mut self) -> Option<Self::Item> {
        let raw = *self.entries.next()?;
        let op = CigarOp::from_raw(raw);
        let len_or_query = ((raw >> 4) & 0x3fff_ffff) as u32;
        let high_bits = raw >> 34;

        Some(match op {
            CigarOp::Match | CigarOp::Mismatch => GraphCigarOp::Aligned {
                op,
                node_id: NodeId(high_bits as i32),
                query_index: len_or_query as usize,
            },
            CigarOp::Deletion => GraphCigarOp::Deletion {
                node_id: NodeId(high_bits as i32),
                len: len_or_query.max(1),
            },
            CigarOp::Insertion | CigarOp::SoftClip | CigarOp::HardClip => GraphCigarOp::QueryRun {
                op,
                end: high_bits as usize,
                len: len_or_query.max(1),
            },
        })
    }
}

/// Low-level alignment result that owns the C-allocated graph cigar
pub struct RawAlignment {
    res: sys::abpoa_res_t,
    // Not Send/Sync: shares upstream global buffers with other abPOA calls
    _not_send_sync: PhantomData<Rc<()>>,
}

impl RawAlignment {
    fn new() -> Self {
        Self {
            res: sys::abpoa_res_t {
                n_cigar: 0,
                m_cigar: 0,
                graph_cigar: ptr::null_mut(),
                node_s: 0,
                node_e: 0,
                query_s: 0,
                query_e: 0,
                n_aln_bases: 0,
                n_matched_bases: 0,
                best_score: 0,
            },
            _not_send_sync: PhantomData,
        }
    }

    fn as_mut_ptr(&mut self) -> *mut sys::abpoa_res_t {
        &mut self.res
    }

    pub(crate) fn as_raw(&self) -> sys::abpoa_res_t {
        self.res
    }

    /// Number of cigar operations produced by the alignment
    pub fn cigar_len(&self) -> i32 {
        self.res.n_cigar
    }

    /// Alignment score reported by abPOA
    pub fn best_score(&self) -> i32 {
        self.res.best_score
    }

    /// Number of query bases that participated in the alignment
    pub fn aligned_bases(&self) -> i32 {
        self.res.n_aln_bases
    }

    /// Number of exact matches reported by abPOA
    pub fn matched_bases(&self) -> i32 {
        self.res.n_matched_bases
    }

    /// Iterate over the decoded graph CIGAR for this alignment
    pub fn cigar(&self) -> GraphCigar<'_> {
        let entries = if self.res.n_cigar <= 0 || self.res.graph_cigar.is_null() {
            &[]
        } else {
            // Safety: `graph_cigar` is allocated and populated by abPOA when `n_cigar > 0`
            unsafe { slice::from_raw_parts(self.res.graph_cigar, self.res.n_cigar as usize) }
        };
        GraphCigar::new(entries)
    }

    /// Decode the graph CIGAR into a human-readable three-line alignment
    ///
    /// `query` must be encoded with the same alphabet used to configure the aligner; a
    /// `Graph` view can be borrowed via [`Aligner::graph`] for use here
    pub fn format_alignment(
        &self,
        graph: &Graph<'_>,
        query: &[u8],
        alphabet: Alphabet,
    ) -> Result<String> {
        let mut graph_row = String::new();
        let mut markers = String::new();
        let mut query_row = String::new();

        for entry in self.cigar() {
            match entry {
                GraphCigarOp::Aligned {
                    op,
                    node_id,
                    query_index,
                } => {
                    let node = graph.node(node_id)?;
                    let graph_base = decode_base(alphabet, node.base);
                    let query_base = *query.get(query_index).ok_or(Error::InvalidInput(
                        "graph cigar referenced a query index past sequence length".into(),
                    ))?;
                    graph_row.push(graph_base);
                    query_row.push(decode_base(alphabet, query_base));
                    markers.push(if op == CigarOp::Match { '|' } else { '.' });
                }
                GraphCigarOp::Deletion { node_id, len } => {
                    let node = graph.node(node_id)?;
                    let graph_base = decode_base(alphabet, node.base);
                    for _ in 0..len.max(1) {
                        graph_row.push(graph_base);
                        markers.push(' ');
                        query_row.push('-');
                    }
                }
                GraphCigarOp::QueryRun { end, len, .. } => {
                    let len = len.max(1) as usize;
                    let start = end.saturating_add(1).saturating_sub(len);
                    for offset in 0..len {
                        let idx = start + offset;
                        let query_base = *query.get(idx).ok_or(Error::InvalidInput(
                            "graph cigar referenced a query index past sequence length".into(),
                        ))?;
                        graph_row.push('-');
                        markers.push(' ');
                        query_row.push(decode_base(alphabet, query_base));
                    }
                }
            }
        }

        Ok(format!("{}\n{}\n{}", graph_row, markers, query_row))
    }
}

impl Drop for RawAlignment {
    fn drop(&mut self) {
        // Safety: `graph_cigar` is allocated by abPOA when `n_cigar > 0`; free with libc to
        // mirror the C allocation and avoid leaking between calls
        if self.res.n_cigar > 0 && !self.res.graph_cigar.is_null() {
            unsafe { libc::free(self.res.graph_cigar as *mut libc::c_void) };
            self.res.graph_cigar = ptr::null_mut();
            self.res.n_cigar = 0;
        }
    }
}

/// Wrapper around `abpoa_t`
///
/// `Aligner` is stateful: methods that generate consensus/MSA/GFA outputs may
/// update the underlying `Parameters` (for example `outputs` or quality-weighted
/// consensus mode). These changes persist for subsequent operations on the same
/// aligner. If you need a clean configuration, construct a new `Aligner` with
/// fresh `Parameters`
pub struct Aligner {
    raw: NonNull<sys::abpoa_t>,
    params: Parameters,
    // Not Send/Sync: abPOA uses global mutable tables and allocators without locking
    _not_send_sync: PhantomData<Rc<()>>,
}

/// Collection of sequences plus optional metadata used as input to MSA calls
pub struct SequenceBatch<'a> {
    sequences: &'a [&'a [u8]],
    names: Option<&'a [&'a str]>,
    quality_weights: Option<&'a [&'a [i32]]>,
}

impl<'a> SequenceBatch<'a> {
    /// Construct a batch from raw sequence slices
    pub fn from_sequences(sequences: &'a [&'a [u8]]) -> Self {
        Self {
            sequences,
            names: None,
            quality_weights: None,
        }
    }

    /// Attach per-sequence names (must match the sequence count)
    pub fn with_names(mut self, names: &'a [&'a str]) -> Self {
        self.names = Some(names);
        self
    }

    /// Attach per-base quality weights (must match sequence count and lengths)
    pub fn with_quality_weights(mut self, quality_weights: &'a [&'a [i32]]) -> Self {
        self.quality_weights = Some(quality_weights);
        self
    }

    pub(crate) fn sequences(&self) -> &'a [&'a [u8]] {
        self.sequences
    }

    pub(crate) fn names(&self) -> Option<&'a [&'a str]> {
        self.names
    }

    pub(crate) fn quality_weights(&self) -> Option<&'a [&'a [i32]]> {
        self.quality_weights
    }
}

impl Aligner {
    /// Allocate a new aligner handle with default parameters
    pub fn new() -> Result<Self> {
        Self::with_params(Parameters::new()?)
    }

    /// Allocate a new aligner handle using the provided parameters
    pub fn with_params(params: Parameters) -> Result<Self> {
        // Safety: calls into the C API to allocate an aligner; abPOA defines the allocation layout
        // and requires `abpoa_free` to release it
        let raw = unsafe { sys::abpoa_init() };
        let raw = NonNull::new(raw).ok_or(Error::NullPointer("abpoa_init returned null"))?;

        Ok(Self {
            raw,
            params,
            _not_send_sync: PhantomData,
        })
    }

    /// Restore an aligner from a serialized graph (FASTA/GFA) in one call (uses default Parameters,
    /// but disables MSA output when `preserve_read_ids` is `false`).
    pub fn from_graph_file<P: AsRef<Path>>(path: P, preserve_read_ids: bool) -> Result<Self> {
        let mut params = Parameters::configure()?;
        if !preserve_read_ids {
            params.set_outputs(OutputMode::CONSENSUS);
        }
        params.set_use_read_ids(preserve_read_ids);
        params.set_incremental_graph_file(path)?;
        let mut aligner = Self::with_params(params)?;
        aligner.restore_graph()?;
        Ok(aligner)
    }

    pub(crate) fn as_mut_ptr(&mut self) -> *mut sys::abpoa_t {
        self.raw.as_ptr()
    }

    pub(crate) fn as_ptr(&self) -> *const sys::abpoa_t {
        self.raw.as_ptr()
    }

    /// Borrow the underlying parameters immutably
    pub fn params(&self) -> &Parameters {
        &self.params
    }

    #[cfg(test)]
    /// Borrow the underlying parameters mutably for test-only reconfiguration
    pub(crate) fn params_mut(&mut self) -> &mut Parameters {
        &mut self.params
    }

    /// Borrow the underlying graph for read-only inspection
    pub fn graph(&self) -> Result<Graph<'_>> {
        // Safety: the graph pointer is owned by this aligner and remains live for `'self`
        let abg = unsafe { (*self.as_ptr()).abg };
        let abs = unsafe { (*self.as_ptr()).abs };
        Graph::new(abg, abs)
    }

    /// Manually add a node with an already-encoded base to the underlying graph
    pub fn add_node(&mut self, base: u8) -> Result<NodeId> {
        let graph = self.graph_mut()?;
        // Safety: the graph pointer is owned by this aligner and live for the duration of the call
        let id = unsafe { sys::abpoa_add_graph_node(graph, base) };
        Self::invalidate_graph(graph);
        Ok(NodeId(id))
    }

    /// Insert or update a directed edge between two nodes
    pub fn add_edge(
        &mut self,
        from: NodeId,
        to: NodeId,
        weight: i32,
        check_existing: bool,
    ) -> Result<()> {
        if weight <= 0 {
            return Err(Error::InvalidInput("edge weight must be positive".into()));
        }
        let graph = self.graph_mut()?;
        Self::validate_node_id(graph, from)?;
        Self::validate_node_id(graph, to)?;
        if from == to {
            return Err(Error::InvalidInput(
                "edge endpoints must differ to avoid self-cycles".into(),
            ));
        }

        // Safety: graph and node ids are valid
        unsafe {
            sys::abpoa_add_graph_edge(
                graph,
                from.0,
                to.0,
                check_existing as i32,
                weight,
                0,
                0,
                0,
                0,
                0,
            );
        }
        Self::invalidate_graph(graph);
        Ok(())
    }

    /// Recompute topological indices for all nodes using the existing buffers
    pub fn refresh_node_indices(&mut self) -> Result<()> {
        let graph = self.graph_mut()?;
        Self::ensure_index_capacity(graph)?;

        // Safety: index buffers are allocated and sized to at least `node_n`
        let src = SentinelNode::Source.as_raw();
        let sink = SentinelNode::Sink.as_raw();
        unsafe { sys::abpoa_BFS_set_node_index(graph, src, sink) };
        Ok(())
    }

    /// Recompute remaining path distances from each node to the sink
    pub fn refresh_node_remaining(&mut self) -> Result<()> {
        let graph = self.graph_mut()?;
        Self::ensure_remain_capacity(graph)?;

        // Safety: remain buffers are allocated and sized to at least `node_n`
        let src = SentinelNode::Source.as_raw();
        let sink = SentinelNode::Sink.as_raw();
        unsafe { sys::abpoa_BFS_set_node_remain(graph, src, sink) };
        Ok(())
    }

    /// Ensure the graph is topologically sorted and dependent buffers are refreshed
    pub fn ensure_topological(&mut self) -> Result<()> {
        let graph = self.graph_mut()?;
        // Safety: abPOA allocates and fills topology buffers when sorting
        unsafe { sys::abpoa_topological_sort(graph, self.params.as_mut_ptr()) };
        Ok(())
    }

    /// Reset the underlying graph and DP matrices to prepare for incremental alignment
    pub fn reset(&mut self, ref_len: usize) -> Result<()> {
        let ref_len = to_i32(ref_len, "reference length exceeds i32")?;
        // Safety: `raw` and `params` are live pointers; abPOA expects `abpoa_reset` before
        // running incremental alignment
        unsafe { sys::abpoa_reset(self.as_mut_ptr(), self.params.as_mut_ptr(), ref_len) };
        Ok(())
    }

    /// Restore a previously serialized graph using the path stored in `Parameters::set_incremental_graph_file`
    pub fn restore_graph(&mut self) -> Result<()> {
        // Safety: abPOA requires `incr_fn` to be set to a valid path; return an error if the
        // caller forgot to configure it
        let raw_params = unsafe { self.params.as_mut_ptr().as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?;
        if raw_params.incr_fn.is_null() {
            return Err(Error::InvalidInput(
                "set an incremental graph path before calling restore_graph".into(),
            ));
        }

        let restored =
            unsafe { sys::abpoa_restore_graph(self.as_mut_ptr(), self.params.as_mut_ptr()) };
        if restored.is_null() {
            return Err(Error::NullPointer("failed to restore graph from file"));
        }
        Ok(())
    }

    /// Compute the exclusive node ids that bound an included region
    pub fn subgraph_nodes(
        &mut self,
        include_begin: NodeId,
        include_end: NodeId,
    ) -> Result<SubgraphRange> {
        let mut exc_beg = 0;
        let mut exc_end = 0;

        // Safety: aligner and parameters are live; abpoa will write the two exclusive node ids
        unsafe {
            sys::abpoa_subgraph_nodes(
                self.as_mut_ptr(),
                self.params.as_mut_ptr(),
                include_begin.0,
                include_end.0,
                &mut exc_beg,
                &mut exc_end,
            )
        };

        Ok(SubgraphRange {
            beg: NodeId(exc_beg),
            end: NodeId(exc_end),
        })
    }

    /// Align a pre-encoded sequence to the current graph and return the raw alignment result
    pub fn align_sequence_raw(&mut self, encoded_seq: &[u8]) -> Result<RawAlignment> {
        if encoded_seq.is_empty() {
            return Err(Error::InvalidInput("cannot align an empty sequence".into()));
        }
        let qlen = to_i32(encoded_seq.len(), "sequence length exceeds i32")?;
        let mut res = RawAlignment::new();

        // abPOA returns -1 when the graph is still empty; allow that case so the first sequence
        // can be added without a pre-existing graph
        let status = unsafe {
            sys::abpoa_align_sequence_to_graph(
                self.as_mut_ptr(),
                self.params.as_mut_ptr(),
                encoded_seq.as_ptr() as *mut u8,
                qlen,
                res.as_mut_ptr(),
            )
        };
        if status != 0 {
            if status == -1 && self.graph_is_empty() {
                return Ok(res);
            }
            return Err(Error::Abpoa {
                func: "abpoa_align_sequence_to_graph",
                code: status,
            });
        }

        Ok(res)
    }

    /// Align a sequence to a specific subgraph and return the raw alignment result
    pub fn align_sequence_to_subgraph(
        &mut self,
        range: SubgraphRange,
        encoded_seq: &[u8],
    ) -> Result<RawAlignment> {
        if encoded_seq.is_empty() {
            return Err(Error::InvalidInput("cannot align an empty sequence".into()));
        }
        let qlen = to_i32(encoded_seq.len(), "sequence length exceeds i32")?;
        let mut res = RawAlignment::new();
        let (beg, end) = range.as_raw();

        let status = unsafe {
            sys::abpoa_align_sequence_to_subgraph(
                self.as_mut_ptr(),
                self.params.as_mut_ptr(),
                beg,
                end,
                encoded_seq.as_ptr() as *mut u8,
                qlen,
                res.as_mut_ptr(),
            )
        };
        if status != 0 {
            if status == -1 && self.graph_is_empty() {
                return Ok(res);
            }
            return Err(Error::Abpoa {
                func: "abpoa_align_sequence_to_subgraph",
                code: status,
            });
        }

        Ok(res)
    }

    fn add_alignment_inner(
        &mut self,
        encoded_seq: &[u8],
        weights: Option<&[i32]>,
        res: &RawAlignment,
        read_id: i32,
        total_reads: i32,
    ) -> Result<()> {
        if encoded_seq.is_empty() {
            return Err(Error::InvalidInput("cannot add an empty sequence".into()));
        }
        if read_id < 0 || total_reads <= 0 {
            return Err(Error::InvalidInput(
                "read_id and total_reads must be non-negative, total_reads > 0".into(),
            ));
        }
        if read_id >= total_reads {
            return Err(Error::InvalidInput(
                format!("read_id {read_id} out of range 0..{total_reads}").into(),
            ));
        }
        if let Some(w) = weights
            && w.len() != encoded_seq.len()
        {
            return Err(Error::InvalidInput(
                "quality weights must match sequence length".into(),
            ));
        }
        let qlen = to_i32(encoded_seq.len(), "sequence length exceeds i32")?;
        self.ensure_sequence_count(total_reads)?;

        if weights.is_some() {
            // Enable quality-weighted consensus for this graph
            self.params.set_use_quality(true);
        }

        let w_ptr = weights
            .map(|w| w.as_ptr() as *mut i32)
            .unwrap_or(ptr::null_mut());

        let status = unsafe {
            sys::abpoa_add_graph_alignment(
                self.as_mut_ptr(),
                self.params.as_mut_ptr(),
                encoded_seq.as_ptr() as *mut u8,
                w_ptr,
                qlen,
                ptr::null_mut(),
                res.as_raw(),
                read_id,
                total_reads,
                1,
            )
        };
        if status != 0 {
            return Err(Error::Abpoa {
                func: "abpoa_add_graph_alignment",
                code: status,
            });
        }

        Ok(())
    }

    /// Add an aligned sequence to the graph using the provided alignment
    pub fn add_alignment(
        &mut self,
        encoded_seq: &[u8],
        res: &RawAlignment,
        read_id: i32,
        total_reads: i32,
    ) -> Result<()> {
        self.add_alignment_inner(encoded_seq, None, res, read_id, total_reads)
    }

    /// Add an aligned sequence to the graph using per-base quality weights
    fn add_alignment_with_quality(
        &mut self,
        encoded_seq: &[u8],
        weights: &[i32],
        res: &RawAlignment,
        read_id: i32,
        total_reads: i32,
    ) -> Result<()> {
        self.add_alignment_inner(encoded_seq, Some(weights), res, read_id, total_reads)
    }

    fn align_and_add_batch(
        &mut self,
        encoded: &[Vec<u8>],
        prepared_weights: Option<&PreparedWeights>,
        read_id_offset: i32,
        total_reads: i32,
    ) -> Result<()> {
        let raw_params = unsafe { self.params.as_mut_ptr().as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?;
        let alphabet = self.alphabet();
        let amb_strand = raw_params.amb_strand() != 0 && alphabet == Alphabet::Dna;
        let max_mat = raw_params.max_mat;

        let abs_ptr = unsafe { (*self.as_mut_ptr()).abs };
        if abs_ptr.is_null() {
            return Err(Error::NullPointer("abpoa sequence container was null"));
        }
        let is_rc_ptr = if amb_strand {
            let ptr = unsafe { (*abs_ptr).is_rc };
            if ptr.is_null() {
                return Err(Error::NullPointer(
                    "abpoa reverse complement flags buffer was null",
                ));
            }
            Some(ptr)
        } else {
            None
        };

        for (idx, seq) in encoded.iter().enumerate() {
            let read_id = read_id_offset
                .checked_add(to_i32(idx, "too many sequences for abpoa")?)
                .ok_or(Error::InvalidInput("too many sequences for abpoa".into()))?;

            if let Some(is_rc_ptr) = is_rc_ptr {
                // Safety: `is_rc_ptr` is allocated to `m_seq` entries by abPOA and `read_id`
                // is within the `total_reads` bound ensured by `ensure_sequence_count`
                unsafe {
                    *is_rc_ptr.add(read_id as usize) = 0;
                }
            }

            let mut alignment = self.align_sequence_raw(seq)?;
            let mut seq_to_add: &[u8] = seq;
            let mut rc_seq_storage: Option<Vec<u8>> = None;
            let mut rc_weights: Option<Vec<i32>> = None;
            let mut is_rc = false;

            if amb_strand && !self.graph_is_empty() {
                let graph_nodes = unsafe {
                    let abg = (*self.as_ptr()).abg;
                    if abg.is_null() {
                        0
                    } else {
                        (*abg).node_n.max(0) as usize
                    }
                };
                let min_len = seq.len().min(graph_nodes.saturating_sub(2));
                let threshold = (min_len as f64) * (max_mat as f64) * 0.3333_f64;
                if (alignment.best_score() as f64) < threshold {
                    let rc_seq: Vec<u8> = seq
                        .iter()
                        .rev()
                        .map(|&b| if b < 4 { 3 - b } else { 4 })
                        .collect();
                    let rc_alignment = self.align_sequence_raw(&rc_seq)?;
                    if rc_alignment.best_score() > alignment.best_score() {
                        alignment = rc_alignment;
                        rc_seq_storage = Some(rc_seq);
                        is_rc = true;
                        if let Some(weights) = prepared_weights {
                            let row = &weights.rows()[idx];
                            rc_weights = Some(row.iter().rev().copied().collect());
                        }
                    }
                }
            }

            if is_rc {
                seq_to_add = rc_seq_storage.as_ref().ok_or(Error::InvalidInput(
                    "reverse complement sequence missing".into(),
                ))?;
            }

            if let Some(weights) = prepared_weights {
                let row = &weights.rows()[idx];
                if is_rc {
                    let rc_row = rc_weights.as_ref().ok_or(Error::InvalidInput(
                        "reverse complement weights missing".into(),
                    ))?;
                    self.add_alignment_with_quality(
                        seq_to_add,
                        rc_row,
                        &alignment,
                        read_id,
                        total_reads,
                    )?;
                } else {
                    self.add_alignment_with_quality(
                        seq_to_add,
                        row,
                        &alignment,
                        read_id,
                        total_reads,
                    )?;
                }
            } else {
                self.add_alignment(seq_to_add, &alignment, read_id, total_reads)?;
            }

            if is_rc && let Some(is_rc_ptr) = is_rc_ptr {
                // Safety: `is_rc_ptr` is allocated to `m_seq` entries by abPOA and `read_id`
                // is within the `total_reads` bound ensured by `ensure_sequence_count`
                unsafe {
                    *is_rc_ptr.add(read_id as usize) = 1;
                }
            }
        }

        Ok(())
    }

    /// Add an aligned sequence to the graph within a subgraph window
    pub fn add_subgraph_alignment(
        &mut self,
        range: SubgraphRange,
        encoded_seq: &[u8],
        res: &RawAlignment,
        read_id: i32,
        total_reads: i32,
        include_both_ends: bool,
    ) -> Result<()> {
        if encoded_seq.is_empty() {
            return Err(Error::InvalidInput("cannot add an empty sequence".into()));
        }
        if read_id < 0 || total_reads <= 0 {
            return Err(Error::InvalidInput(
                "read_id and total_reads must be non-negative, total_reads > 0".into(),
            ));
        }
        if read_id >= total_reads {
            return Err(Error::InvalidInput(
                format!("read_id {read_id} out of range 0..{total_reads}").into(),
            ));
        }
        let qlen = to_i32(encoded_seq.len(), "sequence length exceeds i32")?;
        self.ensure_sequence_count(total_reads)?;
        let (beg, end) = range.as_raw();

        let status = unsafe {
            sys::abpoa_add_subgraph_alignment(
                self.as_mut_ptr(),
                self.params.as_mut_ptr(),
                beg,
                end,
                encoded_seq.as_ptr() as *mut u8,
                ptr::null_mut(),
                qlen,
                ptr::null_mut(),
                res.as_raw(),
                read_id,
                total_reads,
                include_both_ends as i32,
            )
        };
        if status != 0 {
            return Err(Error::Abpoa {
                func: "abpoa_add_subgraph_alignment",
                code: status,
            });
        }

        Ok(())
    }

    /// Build a graph incrementally from a batch of sequences without returning a result yet
    ///
    /// This aligns each sequence against the current graph in input order. Minimizer seeding,
    /// guide-tree partitioning, and progressive POA settings are only used by [`msa`] /
    /// [`msa_encoded`] and do not affect this incremental API
    pub fn msa_in_place(&mut self, batch: SequenceBatch<'_>) -> Result<()> {
        let seqs = batch.sequences();
        if seqs.is_empty() {
            return Ok(());
        }

        self.reset_cached_outputs()?;
        self.params
            .set_use_quality(batch.quality_weights().is_some());
        let alphabet = self.alphabet();
        let encoded = encode_sequences(seqs, alphabet);
        let seq_lens: Vec<i32> = encoded
            .iter()
            .map(|seq| to_i32(seq.len(), "sequence length exceeds i32"))
            .collect::<Result<_>>()?;
        let max_len = encoded.iter().map(Vec::len).max().unwrap_or(0);
        self.reset(max_len)?;

        let total_reads = to_i32(encoded.len(), "too many sequences for abpoa")?;
        self.store_batch_in_abs(&batch, 0, total_reads)?;
        let prepared_weights = prepare_quality_weights(batch.quality_weights(), &seq_lens)?;
        self.align_and_add_batch(&encoded, prepared_weights.as_ref(), 0, total_reads)?;

        Ok(())
    }

    /// Add more sequences to the existing graph, maintaining read ids
    ///
    /// Like [`msa_in_place`], this uses direct sequence-to-graph alignment and ignores minimizer
    /// seeding or guide-tree parameters
    pub fn add_sequences(&mut self, batch: SequenceBatch<'_>) -> Result<()> {
        let new_seqs = batch.sequences();
        if new_seqs.is_empty() {
            return Ok(());
        }

        let alphabet = self.alphabet();
        let encoded = encode_sequences(new_seqs, alphabet);
        let seq_lens: Vec<i32> = encoded
            .iter()
            .map(|seq| to_i32(seq.len(), "sequence length exceeds i32"))
            .collect::<Result<_>>()?;
        let current = self.sequence_count()?;
        let added = to_i32(encoded.len(), "too many sequences for abpoa")?;
        let total_reads = current
            .checked_add(added)
            .ok_or(Error::InvalidInput("too many sequences for abpoa".into()))?;

        self.store_batch_in_abs(&batch, current, total_reads)?;
        let prepared_weights = prepare_quality_weights(batch.quality_weights(), &seq_lens)?;
        self.align_and_add_batch(&encoded, prepared_weights.as_ref(), current, total_reads)?;

        Ok(())
    }

    /// Generate consensus and/or MSA output from the current graph
    pub fn finalize_msa(&mut self, outputs: OutputMode) -> Result<MsaResult> {
        self.finalize_msa_inner(outputs, |abc, alphabet| unsafe {
            MsaResult::from_raw(abc, alphabet)
        })
    }

    /// Generate consensus and/or MSA output from the current graph in the raw abPOA alphabet
    pub fn finalize_msa_encoded(&mut self, outputs: OutputMode) -> Result<EncodedMsaResult> {
        self.finalize_msa_inner(outputs, |abc, _| unsafe { EncodedMsaResult::from_raw(abc) })
    }

    fn finalize_msa_inner<T>(
        &mut self,
        outputs: OutputMode,
        convert: impl FnOnce(*const sys::abpoa_cons_t, Alphabet) -> T,
    ) -> Result<T> {
        crate::runtime::ensure_output_tables();
        self.reset_cached_outputs()?;
        if self.graph_is_empty() {
            let alphabet = self.alphabet();
            return Ok(convert(ptr::null(), alphabet));
        }

        if outputs.is_empty() {
            return Err(Error::InvalidInput(
                "enable consensus and/or msa output before finalizing".into(),
            ));
        }

        let alphabet = self.alphabet();
        self.params.set_outputs_for_call(outputs);

        if outputs.contains(OutputMode::MSA) {
            // Safety: aligner and parameters are live; abPOA will populate `abc` with MSA/consensus
            unsafe { sys::abpoa_generate_rc_msa(self.as_mut_ptr(), self.params.as_mut_ptr()) };
        } else {
            // Safety: aligner and parameters are live; abPOA will populate `abc` with consensus
            unsafe { sys::abpoa_generate_consensus(self.as_mut_ptr(), self.params.as_mut_ptr()) };
        }

        let abc = unsafe { (*self.as_ptr()).abc };
        if abc.is_null() {
            return Err(Error::NullPointer(
                "abpoa returned a null consensus pointer",
            ));
        }

        Ok(convert(abc, alphabet))
    }

    /// Write consensus sequences in FASTA form to a Rust `Write`
    pub fn write_consensus_fasta(&mut self, writer: &mut impl Write) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        crate::runtime::ensure_output_tables();
        self.reset_cached_outputs()?;
        let alphabet = self.alphabet();
        self.params.set_outputs_for_call(OutputMode::CONSENSUS);

        unsafe { sys::abpoa_generate_consensus(self.as_mut_ptr(), self.params.as_mut_ptr()) };
        let abc = unsafe { (*self.as_ptr()).abc };
        if abc.is_null() {
            return Err(Error::NullPointer(
                "abpoa returned a null consensus pointer",
            ));
        }

        let result = unsafe { MsaResult::from_raw(abc, alphabet) };
        let batch_index = unsafe { self.params.as_mut_ptr().as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?
            .batch_index;

        (|| -> Result<()> {
            for (idx, cluster) in result.clusters.iter().enumerate() {
                write!(writer, ">Consensus_sequence")?;
                if batch_index > 0 {
                    write!(writer, "_{}", batch_index)?;
                }
                if result.clusters.len() > 1 {
                    write!(writer, "_{} ", idx + 1)?;
                    for (read_idx, read_id) in cluster.read_ids.iter().enumerate() {
                        if read_idx > 0 {
                            write!(writer, ",")?;
                        }
                        write!(writer, "{}", read_id)?;
                    }
                }
                writeln!(writer)?;
                writeln!(writer, "{}", cluster.consensus)?;
            }
            Ok(())
        })()
    }

    /// Write RC-MSA output (optionally including consensus rows) to a Rust `Write`
    pub fn write_msa_fasta(&mut self, writer: &mut impl Write) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        crate::runtime::ensure_output_tables();
        self.reset_cached_outputs()?;
        self.params
            .set_outputs_for_call(OutputMode::CONSENSUS | OutputMode::MSA);
        let alphabet = self.alphabet();
        let decode_row: fn(&[u8]) -> String = match alphabet {
            Alphabet::Dna => decode_dna,
            Alphabet::AminoAcid => decode_aa,
        };

        unsafe { sys::abpoa_generate_consensus(self.as_mut_ptr(), self.params.as_mut_ptr()) };
        unsafe { sys::abpoa_generate_rc_msa(self.as_mut_ptr(), self.params.as_mut_ptr()) };
        let abc_ptr = unsafe { (*self.as_ptr()).abc };
        let abc = unsafe { abc_ptr.as_ref() }.ok_or(Error::NullPointer(
            "abpoa returned a null consensus pointer",
        ))?;
        let abs_ptr = unsafe { (*self.as_ptr()).abs };
        let abs = unsafe { abs_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa sequence container was null"))?;
        let n_seq = abc.n_seq.max(0) as usize;
        let msa_len = abc.msa_len.max(0) as usize;

        if abc.msa_len <= 0 {
            Ok(())
        } else if abc.msa_base.is_null() {
            Err(Error::NullPointer("abpoa returned a null msa buffer"))
        } else {
            (|| -> Result<()> {
                for idx in 0..n_seq {
                    let name = format_msa_name(abs, idx)?;
                    let row_ptr = unsafe { *abc.msa_base.add(idx) };
                    if row_ptr.is_null() {
                        return Err(Error::NullPointer("abpoa returned a null msa row"));
                    }
                    let row = unsafe { slice::from_raw_parts(row_ptr, msa_len) };
                    writeln!(writer, ">{}", name)?;
                    writeln!(writer, "{}", decode_row(row))?;
                }

                if abc.n_cons > 0 {
                    for cons_idx in 0..abc.n_cons.max(0) as usize {
                        let row_ptr = unsafe { *abc.msa_base.add(n_seq + cons_idx) };
                        if row_ptr.is_null() {
                            return Err(Error::NullPointer(
                                "abpoa returned a null consensus msa row",
                            ));
                        }
                        let ids = consensus_read_ids(abc, cons_idx);
                        let row = unsafe { slice::from_raw_parts(row_ptr, msa_len) };
                        write!(writer, ">Consensus_sequence")?;
                        if abc.n_cons > 1 {
                            write!(writer, "_{} ", cons_idx + 1)?;
                            for (idx, read_id) in ids.iter().enumerate() {
                                if idx > 0 {
                                    write!(writer, ",")?;
                                }
                                write!(writer, "{}", read_id)?;
                            }
                        }
                        writeln!(writer)?;
                        writeln!(writer, "{}", decode_row(row))?;
                    }
                }
                Ok(())
            })()
        }
    }

    /// Write GFA output to the given path
    pub fn write_gfa_to_path<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        crate::runtime::ensure_output_tables();
        self.reset_cached_outputs()?;
        let abs_ptr = unsafe { (*self.as_ptr()).abs };
        let has_consensus = has_consensus_sequence(abs_ptr)?;
        let desired_outputs = if has_consensus {
            // Avoid regenerating consensus if the graph already carries one from a restore
            OutputMode::MSA
        } else {
            OutputMode::CONSENSUS | OutputMode::MSA
        };
        self.params.set_outputs_for_call(desired_outputs);

        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(path.as_ref())?;
        // Safety: `dup` creates an owned fd for `fdopen`/`fclose` without double-closing `file`
        let dup_fd = unsafe { libc::dup(file.as_raw_fd()) };
        if dup_fd < 0 {
            return Err(std::io::Error::last_os_error().into());
        }

        // Safety: `dup_fd` is a valid fd opened for read/write; `fdopen` takes ownership
        let fp = unsafe { libc::fdopen(dup_fd, c"w+".as_ptr() as *const c_char) };

        if fp.is_null() {
            unsafe { libc::close(dup_fd) };
            return Err(Error::NullPointer("failed to open FILE for GFA output"));
        }

        let write_result = {
            unsafe {
                sys::abpoa_generate_gfa(
                    self.as_mut_ptr(),
                    self.params.as_mut_ptr(),
                    fp as *mut sys::FILE,
                );
                libc::fflush(fp);
            }
            Ok(())
        };

        unsafe { libc::fclose(fp) };
        // abpoa_generate_gfa may allocate consensus buffers when out_cons is enabled
        unsafe { sys::abpoa_clean_msa_cons(self.as_mut_ptr()) };

        write_result
    }

    /// Write GFA output to a Rust `Write`
    ///
    /// This is a convenience adapter around [`write_gfa_to_path`] and uses a temporary file
    pub fn write_gfa(&mut self, writer: &mut impl Write) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        let path = temp_gfa_path()?;
        let write_result = self.write_gfa_to_path(&path);
        let read_result = match write_result {
            Ok(()) => {
                let contents = fs::read(&path)?;
                writer.write_all(&contents)?;
                Ok(())
            }
            Err(err) => Err(err),
        };

        let _ = fs::remove_file(&path);
        read_result
    }

    /// Write a GraphViz view of the current partial order graph (POG) to `path`
    ///
    /// Supported extensions:
    /// - `.dot`: write GraphViz DOT
    /// - `.png` / `.pdf`: write an image using GraphViz `dot`
    ///
    /// When writing `.png`/`.pdf`, a DOT file is written next to `path` at `path + ".dot"`
    /// (mirroring upstream)
    ///
    /// `options` can be:
    /// - [`PogDotOptions`] (or `&PogDotOptions`) to use the Rust DOT emitter
    /// - [`PogWriteOptions::Upstream`] to call upstream `abpoa_dump_pog` (uses `system()` and may
    ///   `exit()` on failure)
    pub fn write_pog_to_path<P: AsRef<Path>>(
        &mut self,
        path: P,
        options: impl Into<PogWriteOptions>,
    ) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        let path = path.as_ref();
        let ext = path
            .extension()
            .and_then(|ext| ext.to_str())
            .ok_or(Error::InvalidInput(
                "graph dump path must end with .dot, .png, or .pdf".into(),
            ))?;

        match options.into() {
            PogWriteOptions::Rust(dot_options) => {
                self.ensure_topological()?;
                let graph = self.graph()?;

                if ext == "dot" {
                    let mut file = OpenOptions::new()
                        .write(true)
                        .create(true)
                        .truncate(true)
                        .open(path)?;
                    return crate::output::dot::write_pog_dot(
                        &graph,
                        self.alphabet(),
                        &mut file,
                        &dot_options,
                    );
                }

                if ext != "png" && ext != "pdf" {
                    return Err(Error::InvalidInput(
                        "graph dump path must end with .dot, .png, or .pdf".into(),
                    ));
                }

                let mut dot_path = std::ffi::OsString::from(path.as_os_str());
                dot_path.push(".dot");
                let dot_path = PathBuf::from(dot_path);

                {
                    let mut file = OpenOptions::new()
                        .write(true)
                        .create(true)
                        .truncate(true)
                        .open(&dot_path)?;
                    crate::output::dot::write_pog_dot(
                        &graph,
                        self.alphabet(),
                        &mut file,
                        &dot_options,
                    )?;
                }

                let status = process::Command::new("dot")
                    .arg(format!("-T{ext}"))
                    .arg("-o")
                    .arg(path)
                    .arg(&dot_path)
                    .status()
                    .map_err(|err| {
                        std::io::Error::other(format!("failed to execute `dot`: {err}"))
                    })?;

                if !status.success() {
                    return Err(std::io::Error::other(format!(
                        "`dot` exited with status {status}"
                    ))
                    .into());
                }

                Ok(())
            }
            PogWriteOptions::Upstream => {
                if ext != "png" && ext != "pdf" {
                    return Err(Error::InvalidInput(
                        "upstream graph dump requires a .png or .pdf path".into(),
                    ));
                }

                let path_str = path.to_str().ok_or(Error::InvalidInput(
                    "graph dump path must be valid UTF-8".into(),
                ))?;
                let c_path = CString::new(path_str).map_err(|_| {
                    Error::InvalidInput("graph dump path cannot contain null bytes".into())
                })?;

                let raw = unsafe { self.params.as_mut_ptr().as_mut() }
                    .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?;
                let previous = raw.out_pog;
                // Safety: `c_path` is a valid C string; `strdup` allocates an owned copy
                raw.out_pog = unsafe { libc::strdup(c_path.as_ptr()) };
                if raw.out_pog.is_null() {
                    raw.out_pog = previous;
                    return Err(Error::NullPointer("failed to store graph dump path"));
                }

                // Safety: `self`/`self.params` are valid for the duration of the call
                unsafe { sys::abpoa_dump_pog(self.as_mut_ptr(), self.params.as_mut_ptr()) };

                unsafe {
                    libc::free(raw.out_pog as *mut libc::c_void);
                    raw.out_pog = previous;
                }

                Ok(())
            }
        }
    }

    /// Run abPOA's one-shot MSA on a set of sequences and return the consensus/MSA output
    ///
    /// This calls into abPOA `abpoa_msa`, so minimizer seeding, guide-tree
    /// partitioning, and progressive POA parameters will be used if enabled
    pub fn msa(&mut self, batch: SequenceBatch<'_>, outputs: OutputMode) -> Result<MsaResult> {
        self.msa_one_shot_inner(batch, outputs, |abc, alphabet| unsafe {
            MsaResult::from_raw(abc, alphabet)
        })
    }

    /// Run abPOA's one-shot MSA on a set of sequences and return encoded consensus/MSA output
    pub fn msa_encoded(
        &mut self,
        batch: SequenceBatch<'_>,
        outputs: OutputMode,
    ) -> Result<EncodedMsaResult> {
        self.msa_one_shot_inner(batch, outputs, |abc, _| unsafe {
            EncodedMsaResult::from_raw(abc)
        })
    }

    fn msa_one_shot_inner<T>(
        &mut self,
        batch: SequenceBatch<'_>,
        outputs: OutputMode,
        convert: impl FnOnce(*const sys::abpoa_cons_t, Alphabet) -> T,
    ) -> Result<T> {
        let seqs = batch.sequences();
        if seqs.is_empty() {
            let alphabet = self.alphabet();
            return Ok(convert(ptr::null(), alphabet));
        }
        if seqs.iter().any(|seq| seq.is_empty()) {
            return Err(Error::InvalidInput("cannot align an empty sequence".into()));
        }
        if outputs.is_empty() {
            return Err(Error::InvalidInput(
                "enable consensus and/or msa output to collect results".into(),
            ));
        }

        crate::runtime::ensure_output_tables();
        self.reset_cached_outputs()?;
        self.params
            .set_use_quality(batch.quality_weights().is_some());
        self.params.set_outputs_for_call(outputs);

        let alphabet = self.alphabet();
        let encoded = encode_sequences(seqs, alphabet);
        let mut seq_lens: Vec<i32> = encoded
            .iter()
            .map(|seq| to_i32(seq.len(), "sequence length exceeds i32"))
            .collect::<Result<_>>()?;

        let max_len = encoded.iter().map(Vec::len).max().unwrap_or(0);
        self.reset(max_len)?;

        let n_seq = to_i32(encoded.len(), "too many sequences for abpoa")?;
        let mut seq_ptrs: Vec<*mut u8> =
            encoded.iter().map(|seq| seq.as_ptr() as *mut u8).collect();

        let prepared_weights = prepare_quality_weights(batch.quality_weights(), &seq_lens)?;
        let mut qual_ptrs: Vec<*mut i32> = prepared_weights
            .as_ref()
            .map(|w| {
                w.rows()
                    .iter()
                    .map(|row| row.as_ptr() as *mut i32)
                    .collect()
            })
            .unwrap_or_default();
        let qual_ptr = if qual_ptrs.is_empty() {
            ptr::null_mut()
        } else {
            qual_ptrs.as_mut_ptr()
        };

        let tmp_fp = TempCFile::new()?;
        // Safety: all pointers passed are valid for the duration of the call and match the
        // configured alphabet, abPOA will not keep them after returning
        let status = unsafe {
            sys::abpoa_msa(
                self.as_mut_ptr(),
                self.params.as_mut_ptr(),
                n_seq,
                ptr::null_mut(),
                seq_lens.as_mut_ptr(),
                seq_ptrs.as_mut_ptr(),
                qual_ptr,
                tmp_fp.as_ptr() as *mut sys::FILE,
            )
        };
        if status != 0 {
            return Err(Error::Abpoa {
                func: "abpoa_msa",
                code: status,
            });
        }

        // Store the original sequences and optional names on the aligner for downstream graph
        // inspection helpers
        self.store_batch_in_abs(&batch, 0, n_seq)?;

        let abc = unsafe { (*self.as_ptr()).abc };
        if abc.is_null() {
            return Err(Error::NullPointer(
                "abpoa returned a null consensus pointer",
            ));
        }

        Ok(convert(abc, alphabet))
    }

    fn graph_is_empty(&self) -> bool {
        // Safety: `abg` is owned by the aligner and valid for the lifetime of `self`
        unsafe {
            let abg = (*self.as_ptr()).abg;
            abg.is_null() || (*abg).node_n <= 2
        }
    }

    fn alphabet(&self) -> Alphabet {
        // Parameters default to DNA, only DNA and amino acids are supported
        self.params.get_alphabet().unwrap()
    }

    fn ensure_sequence_count(&mut self, total_reads: i32) -> Result<()> {
        // Safety: `abs` is allocated alongside the aligner; update the total number of reads the
        // graph accounts for so consensus/MSA generation sizes buffers correctly
        let abs_ptr = unsafe { (*self.as_mut_ptr()).abs };
        if abs_ptr.is_null() {
            return Err(Error::NullPointer("abpoa sequence container was null"));
        }
        // Safety: `abs_ptr` is owned by this aligner. abPOA reallocates its internal arrays
        // (seq/name/is_rc/etc.) based on `n_seq`
        let abs = unsafe { abs_ptr.as_mut() }
            .ok_or(Error::NullPointer("abpoa sequence container was null"))?;
        if abs.n_seq < total_reads {
            abs.n_seq = total_reads;
            unsafe {
                sys::abpoa_realloc_seq(abs_ptr);
            }
        }
        Ok(())
    }

    fn store_batch_in_abs(
        &mut self,
        batch: &SequenceBatch<'_>,
        read_id_offset: i32,
        total_reads: i32,
    ) -> Result<()> {
        let seqs = batch.sequences();
        let names = batch.names();
        if let Some(names) = names
            && names.len() != seqs.len()
        {
            return Err(Error::InvalidInput(
                "names length must match sequence count".into(),
            ));
        }

        let abs_ptr = unsafe { (*self.as_mut_ptr()).abs };
        if abs_ptr.is_null() {
            return Err(Error::NullPointer("abpoa sequence container was null"));
        }

        // Safety: `abs_ptr` is owned by this aligner and stays valid for the lifetime of `self`
        // We update `n_seq` and let abPOA grow its internal arrays to fit
        let abs = unsafe { abs_ptr.as_mut() }
            .ok_or(Error::NullPointer("abpoa sequence container was null"))?;
        abs.n_seq = total_reads;
        unsafe {
            sys::abpoa_realloc_seq(abs_ptr);
        }

        for (idx, seq) in seqs.iter().enumerate() {
            let read_id = read_id_offset
                .checked_add(to_i32(idx, "too many sequences for abpoa")?)
                .ok_or(Error::InvalidInput("too many sequences for abpoa".into()))?;
            let read_idx = usize::try_from(read_id)
                .map_err(|_| Error::InvalidInput("read id cannot be negative".into()))?;

            let seq_len = to_i32(seq.len(), "sequence length exceeds i32")?;
            // Safety: `abs.seq` points to an array sized to at least `n_seq` entries above
            unsafe {
                sys::abpoa_cpy_str(abs.seq.add(read_idx), seq.as_ptr() as *mut i8, seq_len);
            }

            if let Some(names) = names {
                let name_bytes = names[idx].as_bytes();
                let name_len = to_i32(name_bytes.len(), "name length exceeds i32")?;
                // Safety: `abs.name` points to an array sized to at least `n_seq` entries above
                unsafe {
                    sys::abpoa_cpy_str(
                        abs.name.add(read_idx),
                        name_bytes.as_ptr() as *mut i8,
                        name_len,
                    );
                }
            }
        }

        Ok(())
    }

    fn graph_mut(&mut self) -> Result<&mut sys::abpoa_graph_t> {
        // Safety: `abg` is allocated alongside the aligner and lives for `'self`
        let abg = unsafe { (*self.as_mut_ptr()).abg };
        let mut abg =
            NonNull::new(abg).ok_or(Error::NullPointer("abpoa graph pointer was null"))?;
        // Safety: the graph pointer is valid for the duration of the borrow
        Ok(unsafe { abg.as_mut() })
    }

    fn validate_node_id(graph: &sys::abpoa_graph_t, id: NodeId) -> Result<()> {
        if id.0 < 0 {
            return Err(Error::InvalidInput("node id cannot be negative".into()));
        }
        if id.0 >= graph.node_n {
            return Err(Error::InvalidInput(
                "node id out of bounds for current graph".into(),
            ));
        }
        Ok(())
    }

    fn ensure_index_capacity(graph: &sys::abpoa_graph_t) -> Result<()> {
        if graph.index_to_node_id.is_null() || graph.node_id_to_index.is_null() {
            return Err(Error::NullPointer(
                "graph indices not allocated; call ensure_topological first",
            ));
        }
        if graph.index_rank_m < graph.node_n {
            return Err(Error::InvalidInput(
                "graph indices too small; call ensure_topological to rebuild them".into(),
            ));
        }
        Ok(())
    }

    fn ensure_remain_capacity(graph: &sys::abpoa_graph_t) -> Result<()> {
        Self::ensure_index_capacity(graph)?;
        if graph.node_id_to_max_remain.is_null() {
            return Err(Error::NullPointer(
                "remain buffer not allocated; call ensure_topological first",
            ));
        }
        Ok(())
    }

    fn invalidate_graph(graph: &mut sys::abpoa_graph_t) {
        graph.set_is_topological_sorted(0);
        graph.set_is_called_cons(0);
        graph.set_is_set_msa_rank(0);
    }

    fn sequence_count(&self) -> Result<i32> {
        // Safety: `abs` is owned alongside the aligner and valid for the lifetime of `self`
        let abs = unsafe { (*self.as_ptr()).abs };
        if abs.is_null() {
            return Err(Error::NullPointer("abpoa sequence container was null"));
        }
        // Safety: we only read a single integer field
        let count = unsafe { (*abs).n_seq };
        Ok(count.max(0))
    }

    fn reset_cached_outputs(&mut self) -> Result<()> {
        unsafe { sys::abpoa_clean_msa_cons(self.as_mut_ptr()) };
        let graph = self.graph_mut()?;
        graph.set_is_called_cons(0);
        graph.set_is_set_msa_rank(0);
        Ok(())
    }
}

impl Drop for Aligner {
    fn drop(&mut self) {
        // Safety: `raw` came from `abpoa_init` and is exclusively owned here; abPOA expects a
        // matching `abpoa_free` for each init
        unsafe { sys::abpoa_free(self.raw.as_ptr()) }
    }
}

struct TempCFile {
    fp: *mut libc::FILE,
}

impl TempCFile {
    fn new() -> Result<Self> {
        // Safety: `tmpfile` returns an owned FILE* or null on error
        let fp = unsafe { libc::tmpfile() };
        if fp.is_null() {
            return Err(Error::NullPointer(
                "failed to open temporary FILE for abPOA output",
            ));
        }
        Ok(Self { fp })
    }

    fn as_ptr(&self) -> *mut libc::FILE {
        self.fp
    }
}

impl Drop for TempCFile {
    fn drop(&mut self) {
        if !self.fp.is_null() {
            unsafe {
                libc::fclose(self.fp);
            }
        }
    }
}

struct PreparedWeights {
    _storage: Vec<Vec<i32>>,
}

impl PreparedWeights {
    fn rows(&self) -> &[Vec<i32>] {
        &self._storage
    }
}

fn encode_sequences(seqs: &[&[u8]], alphabet: Alphabet) -> Vec<Vec<u8>> {
    match alphabet {
        Alphabet::Dna => seqs.iter().map(|seq| encode_dna(seq)).collect(),
        Alphabet::AminoAcid => seqs.iter().map(|seq| encode_aa(seq)).collect(),
    }
}

fn prepare_quality_weights(
    weights: Option<&[&[i32]]>,
    lengths: &[i32],
) -> Result<Option<PreparedWeights>> {
    let Some(weights) = weights else {
        return Ok(None);
    };
    if weights.len() != lengths.len() {
        return Err(Error::InvalidInput(
            "quality weights length must match sequence count".into(),
        ));
    }

    let mut storage = Vec::with_capacity(weights.len());
    for (row, &len) in weights.iter().zip(lengths) {
        let expected = len.max(0) as usize;
        if row.len() != expected {
            return Err(Error::InvalidInput(
                "quality weights must match each sequence length".into(),
            ));
        }
        storage.push(row.to_vec());
    }
    Ok(Some(PreparedWeights { _storage: storage }))
}

fn has_consensus_sequence(abs: *const sys::abpoa_seq_t) -> Result<bool> {
    let Some(abs) = (unsafe { abs.as_ref() }) else {
        return Err(Error::NullPointer("abpoa sequence container was null"));
    };
    if abs.name.is_null() {
        return Ok(false);
    }
    let count = abs.n_seq.max(0) as usize;
    for idx in 0..count {
        let name = unsafe { abs.name.add(idx).as_ref() }
            .ok_or(Error::NullPointer("abpoa sequence name pointer was null"))?;
        if name.l > 0 && !name.s.is_null() {
            let bytes = unsafe { slice::from_raw_parts(name.s as *const u8, name.l as usize) };
            if bytes.starts_with(b"Consensus_sequence") {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn decode_base(alphabet: Alphabet, code: u8) -> char {
    match alphabet {
        Alphabet::Dna => decode_dna_code(code),
        Alphabet::AminoAcid => decode_aa_code(code),
    }
}

fn to_i32(value: usize, context: &'static str) -> Result<i32> {
    i32::try_from(value).map_err(|_| Error::InvalidInput(context.into()))
}

fn sequence_name(abs: &sys::abpoa_seq_t, idx: usize) -> Result<String> {
    if abs.name.is_null() {
        return Ok(format!("Seq_{}", idx + 1));
    }
    let kstr = unsafe { abs.name.add(idx).as_ref() }
        .ok_or(Error::NullPointer("abpoa sequence name pointer was null"))?;
    if kstr.l > 0 && !kstr.s.is_null() {
        // Safety: `s` points to a buffer at least `l` bytes long when `l > 0`
        let bytes = unsafe { slice::from_raw_parts(kstr.s as *const u8, kstr.l as usize) };
        Ok(String::from_utf8_lossy(bytes).into_owned())
    } else {
        Ok(format!("Seq_{}", idx + 1))
    }
}

fn sequence_is_rc(abs: &sys::abpoa_seq_t, idx: usize) -> Result<bool> {
    if abs.is_rc.is_null() {
        return Ok(false);
    }
    // Safety: `is_rc` is sized to `n_seq` by abPOA; `idx` is derived from that bound
    let flag = unsafe { *abs.is_rc.add(idx) };
    Ok(flag != 0)
}

fn format_msa_name(abs: &sys::abpoa_seq_t, idx: usize) -> Result<String> {
    let mut name = sequence_name(abs, idx)?;
    if sequence_is_rc(abs, idx)? {
        name.push_str("_reverse_complement");
    }
    Ok(name)
}

fn consensus_read_ids(abc: &sys::abpoa_cons_t, idx: usize) -> Vec<i32> {
    if abc.clu_n_seq.is_null() || abc.clu_read_ids.is_null() {
        return Vec::new();
    }
    // Safety: read count and pointers are filled by abPOA when consensus is generated
    let count = unsafe { *abc.clu_n_seq.add(idx) }.max(0) as usize;
    if count == 0 {
        return Vec::new();
    }
    let ids_ptr = unsafe { *abc.clu_read_ids.add(idx) };
    if ids_ptr.is_null() {
        return Vec::new();
    }
    unsafe { slice::from_raw_parts(ids_ptr, count) }.to_vec()
}

fn temp_gfa_path() -> std::io::Result<PathBuf> {
    let mut path = env::temp_dir();
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    path.push(format!("abpoa_gfa_{}_{}.gfa", process::id(), timestamp));
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aligner_init_and_free() {
        let mut aligner = Aligner::new().unwrap();
        assert!(!aligner.as_mut_ptr().is_null());
    }

    #[test]
    fn from_graph_file_restores_graph_and_toggle_read_ids() {
        let mut aligner = Aligner::new().unwrap();
        aligner
            .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
            .unwrap();

        let path = std::env::temp_dir().join(format!(
            "abpoa_restore_{}_{}.gfa",
            process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        aligner.write_gfa_to_path(&path).unwrap();

        let mut restored = Aligner::from_graph_file(&path, false).unwrap();
        let graph = restored.graph().unwrap();
        assert!(
            !graph.is_empty(),
            "restored graph should include aligned nodes"
        );
        let raw_params = unsafe { restored.params_mut().as_mut_ptr().as_ref() }.unwrap();
        assert!(
            raw_params.use_read_ids() == 0,
            "preserve_read_ids flag should propagate to parameters"
        );

        std::fs::remove_file(&path).unwrap_or_default();
    }

    #[test]
    fn raw_alignment_exposes_cigar_and_counts() {
        let mut aligner = Aligner::new().unwrap();
        {
            let params = aligner.params_mut();
            params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
            params.set_alphabet(Alphabet::Dna).unwrap();
        }

        let ref_seq = encode_dna(b"ACGT");
        let query = encode_dna(b"AGGT");
        let max_len = ref_seq.len().max(query.len());
        aligner.reset(max_len).unwrap();

        let first = aligner.align_sequence_raw(&ref_seq).unwrap();
        aligner.add_alignment(&ref_seq, &first, 0, 2).unwrap();

        let second = aligner.align_sequence_raw(&query).unwrap();

        assert_eq!(second.aligned_bases(), query.len() as i32);
        assert!(second.matched_bases() < second.aligned_bases());

        let ops: Vec<_> = second.cigar().collect();
        assert_eq!(ops.len() as i32, second.cigar_len());
        assert!(
            ops.iter()
                .any(|op| matches!(op, GraphCigarOp::Aligned { .. })),
            "graph cigar should contain aligned nodes"
        );
    }

    #[test]
    fn dna_msa_variants() {
        let sequences = [b"ACGT".as_ref(), b"ACGT".as_ref(), b"ACGG".as_ref()];
        let mut aligner = Aligner::new().unwrap();
        let decoded = aligner
            .msa(
                SequenceBatch::from_sequences(&sequences),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(decoded.msa.len(), sequences.len());
        assert_eq!(decoded.msa[0].len(), 4);
        assert!(!decoded.clusters.is_empty());
        assert_eq!(decoded.clusters[0].consensus, "ACGT");

        let encoded = aligner
            .msa_encoded(
                SequenceBatch::from_sequences(&sequences),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(encoded.msa.len(), decoded.msa.len());
        for (enc_row, dec_row) in encoded.msa.iter().zip(decoded.msa.iter()) {
            assert_eq!(decode_dna(enc_row), *dec_row);
        }

        assert_eq!(encoded.clusters.len(), decoded.clusters.len());
        for (enc_cluster, dec_cluster) in encoded.clusters.iter().zip(decoded.clusters.iter()) {
            assert_eq!(enc_cluster.read_ids, dec_cluster.read_ids);
            assert_eq!(decode_dna(&enc_cluster.consensus), dec_cluster.consensus);
            assert_eq!(enc_cluster.coverage, dec_cluster.coverage);
        }

        let sequences_with_gap = [b"ACGT".as_ref(), b"AGT".as_ref()];
        let mut aligner = Aligner::new().unwrap();
        let gap_result = aligner
            .msa(
                SequenceBatch::from_sequences(&sequences_with_gap),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(gap_result.msa.len(), sequences_with_gap.len());
        assert!(
            gap_result
                .msa
                .iter()
                .all(|row| row.len() == gap_result.msa[0].len()),
            "msa rows should be padded to the same length"
        );
        assert!(
            gap_result.msa.iter().any(|row| row.contains('-')),
            "msa should surface gap sentinels as '-'"
        );

        let sequences_quality = [b"ACGT".as_ref(), b"ACGT".as_ref()];
        let qual_a = [30, 30, 30, 30];
        let qual_b = [10, 10, 10, 10];
        let qualities: [&[i32]; 2] = [qual_a.as_slice(), qual_b.as_slice()];
        let mut aligner = Aligner::new().unwrap();

        let qual_result = aligner
            .msa(
                SequenceBatch::from_sequences(&sequences_quality).with_quality_weights(&qualities),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(qual_result.msa.len(), sequences_quality.len());
        assert!(!qual_result.clusters.is_empty());
    }

    #[test]
    fn consensus_writer_matches_alignment() {
        let mut aligner = Aligner::new().unwrap();
        aligner
            .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGT"]))
            .unwrap();

        let mut buffer = Vec::new();
        aligner.write_consensus_fasta(&mut buffer).unwrap();
        let output = String::from_utf8(buffer).unwrap();
        assert!(
            output.contains(">Consensus_sequence"),
            "consensus header should be present"
        );
        assert!(
            output.contains("ACGT"),
            "consensus sequence should match aligned reads"
        );
    }

    #[test]
    fn msa_writer_includes_names_and_reverse_complements() {
        let mut params = Parameters::configure().unwrap();
        params.set_ambiguous_strand(true);
        let mut aligner = Aligner::with_params(params).unwrap();
        let sequences = [b"ACGTTGC".as_ref(), b"GCAACGT".as_ref()]; // read2 is the RC of read1
        let names = ["read1", "read2"];

        aligner
            .msa(
                SequenceBatch::from_sequences(&sequences).with_names(&names),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        let mut buffer = Vec::new();
        aligner.write_msa_fasta(&mut buffer).unwrap();
        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains(">read1"), "msa should include read1 header");
        assert!(output.contains(">read2"), "msa should include read2 header");
        assert!(
            output.contains(">Consensus_sequence"),
            "msa writer should include consensus rows"
        );
        assert!(
            output.contains("_reverse_complement"),
            "msa output should flag reverse complements when ambiguous strand is enabled"
        );
    }

    #[test]
    fn gfa_writer_generates_header() {
        let mut aligner = Aligner::new().unwrap();
        aligner
            .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
            .unwrap();

        let mut buffer = Vec::new();
        aligner.write_gfa(&mut buffer).unwrap();
        let output = String::from_utf8(buffer).unwrap();
        assert!(
            output.contains("H\tVN:Z:1.0"),
            "GFA header should be present"
        );
        assert!(
            output.contains("\nS\t"),
            "GFA output should include at least one segment"
        );
    }

    #[test]
    fn graph_alignment_round_trip() {
        let mut aligner = Aligner::new().unwrap();
        {
            let params = aligner.params_mut();
            params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
            params.set_alphabet(Alphabet::Dna).unwrap();
        }

        let seqs = [encode_dna(b"ACGT"), encode_dna(b"ACGG")];
        let max_len = seqs.iter().map(Vec::len).max().unwrap();
        aligner.reset(max_len).unwrap();

        let first = aligner.align_sequence_raw(&seqs[0]).unwrap();
        aligner.add_alignment(&seqs[0], &first, 0, 2).unwrap();

        let second = aligner.align_sequence_raw(&seqs[1]).unwrap();
        aligner.add_alignment(&seqs[1], &second, 1, 2).unwrap();

        let params_ptr = {
            let params = aligner.params_mut();
            params.as_mut_ptr()
        };
        unsafe { sys::abpoa_generate_consensus(aligner.as_mut_ptr(), params_ptr) };
        let abc = unsafe { (*aligner.as_ptr()).abc };
        let result = unsafe { MsaResult::from_raw(abc, Alphabet::Dna) };
        unsafe { sys::abpoa_clean_msa_cons(aligner.as_mut_ptr()) };

        let mut expected_aligner = Aligner::new().unwrap();
        let expected = expected_aligner
            .msa(
                SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(result.clusters.len(), 1);
        assert_eq!(expected.clusters.len(), 1);
        assert_eq!(
            result.clusters[0].consensus, expected.clusters[0].consensus,
            "incremental consensus should match one-shot run"
        );
    }

    #[test]
    fn incremental_api_matches_one_shot() {
        let mut aligner = Aligner::new().unwrap();
        aligner
            .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
            .unwrap();
        let incremental = aligner
            .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
            .unwrap();

        let mut expected_aligner = Aligner::new().unwrap();
        let expected = expected_aligner
            .msa(
                SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(incremental.clusters.len(), expected.clusters.len());
        assert_eq!(
            incremental.clusters[0].consensus, expected.clusters[0].consensus,
            "incremental consensus should match one-shot run"
        );
        assert_eq!(incremental.msa.len(), expected.msa.len());
    }

    #[test]
    fn finalize_msa_recomputes_after_cleanup() {
        let mut aligner = Aligner::new().unwrap();
        aligner
            .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGT"]))
            .unwrap();

        let first = aligner
            .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
            .unwrap();
        assert!(!first.clusters.is_empty(), "initial consensus should exist");

        let second = aligner
            .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
            .unwrap();
        assert_eq!(second.clusters.len(), first.clusters.len());
        assert_eq!(second.clusters[0].consensus, first.clusters[0].consensus);
        assert_eq!(second.msa.len(), first.msa.len());
    }

    #[test]
    fn adding_sequences_updates_graph() {
        let mut aligner = Aligner::new().unwrap();
        aligner
            .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
            .unwrap();
        aligner
            .add_sequences(SequenceBatch::from_sequences(&[b"ACGA"]))
            .unwrap();

        let incremental = aligner
            .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
            .unwrap();

        let mut expected_aligner = Aligner::new().unwrap();
        let expected = expected_aligner
            .msa(
                SequenceBatch::from_sequences(&[b"ACGT", b"ACGG", b"ACGA"]),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(incremental.msa.len(), expected.msa.len());
        assert_eq!(
            incremental.clusters[0].consensus, expected.clusters[0].consensus,
            "incremental consensus should match one-shot run"
        );
    }

    #[test]
    fn subgraph_alignment_matches_one_shot() {
        let mut aligner = Aligner::new().unwrap();
        {
            let params = aligner.params_mut();
            params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
            params.set_alphabet(Alphabet::Dna).unwrap();
        }

        let seqs = [
            encode_dna(b"ACGT"),
            encode_dna(b"ACGG"),
            encode_dna(b"ACGA"),
        ];
        let max_len = seqs.iter().map(Vec::len).max().unwrap();
        aligner.reset(max_len).unwrap();

        let total_reads = to_i32(seqs.len(), "too many sequences for abpoa").unwrap();
        let whole_graph = SubgraphRange {
            beg: SentinelNode::Source.as_node_id(),
            end: SentinelNode::Sink.as_node_id(),
        };

        let first = aligner
            .align_sequence_to_subgraph(whole_graph, &seqs[0])
            .unwrap();
        aligner
            .add_subgraph_alignment(whole_graph, &seqs[0], &first, 0, total_reads, false)
            .unwrap();

        let last_node = to_i32(seqs[0].len() + 1, "sequence length exceeds i32").unwrap();
        let range = aligner
            .subgraph_nodes(NodeId(2), NodeId(last_node))
            .unwrap();
        assert_eq!(range.beg, SentinelNode::Source.as_node_id());
        assert_eq!(range.end, SentinelNode::Sink.as_node_id());

        let second = aligner.align_sequence_to_subgraph(range, &seqs[1]).unwrap();
        aligner
            .add_subgraph_alignment(range, &seqs[1], &second, 1, total_reads, false)
            .unwrap();

        let third = aligner.align_sequence_to_subgraph(range, &seqs[2]).unwrap();
        aligner
            .add_subgraph_alignment(range, &seqs[2], &third, 2, total_reads, false)
            .unwrap();

        let incremental = aligner
            .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
            .unwrap();

        let mut expected_aligner = Aligner::new().unwrap();
        let expected = expected_aligner
            .msa(
                SequenceBatch::from_sequences(&[b"ACGT", b"ACGG", b"ACGA"]),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(
            incremental.clusters[0].consensus,
            expected.clusters[0].consensus
        );
        assert_eq!(incremental.msa.len(), expected.msa.len());
        assert_eq!(incremental.msa[0].len(), expected.msa[0].len());
    }

    #[test]
    fn format_alignment_renders_indels() {
        fn check_case<F>(reference: &[u8], query: &[u8], expected_pretty: &str, expects_op: F)
        where
            F: Fn(&GraphCigarOp) -> bool,
        {
            let mut aligner = Aligner::new().unwrap();
            {
                let params = aligner.params_mut();
                params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
                params.set_alphabet(Alphabet::Dna).unwrap();
            }

            let reference = encode_dna(reference);
            let query = encode_dna(query);
            let max_len = reference.len().max(query.len());
            aligner.reset(max_len).unwrap();

            let base_alignment = aligner.align_sequence_raw(&reference).unwrap();
            aligner
                .add_alignment(&reference, &base_alignment, 0, 2)
                .unwrap();

            let query_alignment = aligner.align_sequence_raw(&query).unwrap();
            let ops: Vec<_> = query_alignment.cigar().collect();
            assert!(
                ops.iter().any(expects_op),
                "graph cigar should include expected indel"
            );

            let graph = aligner.graph().unwrap();
            let pretty = query_alignment
                .format_alignment(&graph, &query, Alphabet::Dna)
                .unwrap();
            assert_eq!(pretty, expected_pretty);
        }

        check_case(b"ACGT", b"AACGT", "-ACGT\n ||||\nAACGT", |op| {
            matches!(
                op,
                GraphCigarOp::QueryRun {
                    op: CigarOp::Insertion,
                    len,
                    ..
                } if *len >= 1
            )
        });

        check_case(
            b"ACGT",
            b"ACG",
            "ACGT\n||| \nACG-",
            |op| matches!(op, GraphCigarOp::Deletion { len, .. } if *len >= 1),
        );
    }

    #[test]
    fn msa_encoded_matches_decoded_for_amino_acids() {
        let mut aligner = Aligner::new().unwrap();
        {
            let params = aligner.params_mut();
            params.set_alphabet(Alphabet::AminoAcid).unwrap();
        }

        let sequences = [b"ACDE".as_ref(), b"ACDF".as_ref()];

        let decoded = aligner
            .msa(
                SequenceBatch::from_sequences(&sequences),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();
        assert_eq!(decoded.msa.len(), sequences.len());
        assert!(
            !decoded.clusters.is_empty(),
            "consensus cluster should be returned for amino acid input"
        );
        assert!(
            decoded.clusters[0].consensus.starts_with("ACD"),
            "consensus should decode amino acid alphabet"
        );

        let encoded = aligner
            .msa_encoded(
                SequenceBatch::from_sequences(&sequences),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();

        assert_eq!(encoded.msa.len(), decoded.msa.len());
        for (enc_row, dec_row) in encoded.msa.iter().zip(decoded.msa.iter()) {
            assert_eq!(decode_aa(enc_row), *dec_row);
        }

        assert_eq!(encoded.clusters.len(), decoded.clusters.len());
        for (enc_cluster, dec_cluster) in encoded.clusters.iter().zip(decoded.clusters.iter()) {
            assert_eq!(enc_cluster.read_ids, dec_cluster.read_ids);
            assert_eq!(decode_aa(&enc_cluster.consensus), dec_cluster.consensus);
            assert_eq!(enc_cluster.coverage, dec_cluster.coverage);
        }
    }

    #[test]
    fn manual_graph_editing_and_topology_sorting() {
        let mut aligner = Aligner::new().unwrap();
        {
            let params = aligner.params_mut();
            params.set_alphabet(Alphabet::Dna).unwrap();
        }

        let base = encode_dna(b"A")[0];
        let src = SentinelNode::Source.as_node_id();
        let sink = SentinelNode::Sink.as_node_id();
        let mid = aligner.add_node(base).unwrap();
        aligner.add_edge(src, mid, 1, false).unwrap();
        aligner.add_edge(mid, sink, 1, false).unwrap();

        aligner.ensure_topological().unwrap();
        let graph = aligner.graph().unwrap();
        assert_eq!(graph.node_count(), 3);

        let src_idx = graph.node_index(src).unwrap();
        let mid_idx = graph.node_index(mid).unwrap();
        let sink_idx = graph.node_index(sink).unwrap();
        assert!(src_idx < mid_idx && mid_idx < sink_idx);

        let weight = graph.node_weight(src).unwrap();
        assert_eq!(weight, 1);
    }

    #[test]
    fn refresh_helpers_require_allocated_buffers() {
        let mut aligner = Aligner::new().unwrap();
        let src = SentinelNode::Source.as_node_id();
        let sink = SentinelNode::Sink.as_node_id();
        let mid = aligner.add_node(encode_dna(b"C")[0]).unwrap();
        aligner.add_edge(src, mid, 1, false).unwrap();
        aligner.add_edge(mid, sink, 1, false).unwrap();

        assert!(
            matches!(
                aligner.refresh_node_indices(),
                Err(Error::NullPointer(_) | Error::InvalidInput(_))
            ),
            "refresh should fail until topology buffers exist"
        );
        assert!(
            matches!(
                aligner.refresh_node_remaining(),
                Err(Error::NullPointer(_) | Error::InvalidInput(_))
            ),
            "refresh should fail until topology buffers exist"
        );

        aligner.ensure_topological().unwrap();

        aligner.refresh_node_indices().unwrap();
        aligner.refresh_node_remaining().unwrap();
    }
}
