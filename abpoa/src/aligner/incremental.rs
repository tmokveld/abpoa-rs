use super::{
    encode_sequences, to_i32, validate_quality_weights, Aligner, RawAlignment, SequenceBatch,
};
use crate::encode::Alphabet;
use crate::graph::Graph;
use crate::params::{NodeId, OutputMode, Parameters, SentinelNode};
use crate::result::{EncodedMsaResult, EncodedMsaView, MsaResult};
use crate::{sys, Error, Result};
use std::{marker::PhantomData, path::Path, ptr, ptr::NonNull};

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
            graph_tracks_read_ids: true,
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

    /// Borrow the underlying parameters mutably.
    ///
    /// Note: changing settings like `use_read_ids` does not retroactively update an existing
    /// graph. If you build a graph without read ids, MSA/GFA/multi-consensus outputs are not
    /// available unless you rebuild the graph with read ids enabled.
    pub fn params_mut(&mut self) -> &mut Parameters {
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
        // `abpoa_reset` clears the graph; read-id tracking starts fresh.
        self.graph_tracks_read_ids = true;
        Ok(())
    }

    /// Restore a previously serialized graph using the path stored in `Parameters::set_incremental_graph_file`
    pub fn restore_graph(&mut self) -> Result<()> {
        let empty_before = self.graph_is_empty();
        // Safety: abPOA requires `incr_fn` to be set to a valid path; return an error if the
        // caller forgot to configure it
        let params_ptr = self.params.as_mut_ptr();
        let raw_params = unsafe { params_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?;
        if raw_params.incr_fn.is_null() {
            return Err(Error::InvalidInput(
                "set an incremental graph path before calling restore_graph".into(),
            ));
        }

        let restored = unsafe { sys::abpoa_restore_graph(self.as_mut_ptr(), params_ptr) };
        if restored.is_null() {
            return Err(Error::NullPointer("failed to restore graph from file"));
        }
        let uses_read_ids = raw_params.use_read_ids() != 0;
        if empty_before {
            self.graph_tracks_read_ids = uses_read_ids;
        } else if !uses_read_ids {
            self.graph_tracks_read_ids = false;
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
        if let Some(w) = weights {
            if w.len() != encoded_seq.len() {
                return Err(Error::InvalidInput(
                    "quality weights must match sequence length".into(),
                ));
            }
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

        let params_ptr = self.params.as_mut_ptr();
        let status = unsafe {
            sys::abpoa_add_graph_alignment(
                self.as_mut_ptr(),
                params_ptr,
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
        let uses_read_ids = unsafe { params_ptr.as_ref() }
            .map(|raw| raw.use_read_ids() != 0)
            .unwrap_or(false);
        if !uses_read_ids {
            self.graph_tracks_read_ids = false;
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
        quality_weights: Option<&[&[i32]]>,
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
                        if let Some(weights) = quality_weights {
                            let row = weights[idx];
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

            if let Some(weights) = quality_weights {
                let row = weights[idx];
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

            if is_rc {
                if let Some(is_rc_ptr) = is_rc_ptr {
                    // Safety: `is_rc_ptr` is allocated to `m_seq` entries by abPOA and `read_id`
                    // is within the `total_reads` bound ensured by `ensure_sequence_count`
                    unsafe {
                        *is_rc_ptr.add(read_id as usize) = 1;
                    }
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

        let params_ptr = self.params.as_mut_ptr();
        let status = unsafe {
            sys::abpoa_add_subgraph_alignment(
                self.as_mut_ptr(),
                params_ptr,
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
        let uses_read_ids = unsafe { params_ptr.as_ref() }
            .map(|raw| raw.use_read_ids() != 0)
            .unwrap_or(false);
        if !uses_read_ids {
            self.graph_tracks_read_ids = false;
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
        let quality_weights = validate_quality_weights(batch.quality_weights(), &seq_lens)?;
        self.align_and_add_batch(&encoded, quality_weights, 0, total_reads)?;

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
        let quality_weights = validate_quality_weights(batch.quality_weights(), &seq_lens)?;
        self.align_and_add_batch(&encoded, quality_weights, current, total_reads)?;

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

    /// Generate consensus and/or MSA output from the current graph as zero-copy encoded views
    pub fn finalize_msa_view_encoded(&mut self, outputs: OutputMode) -> Result<EncodedMsaView<'_>> {
        self.finalize_msa_inner(outputs, EncodedMsaView::new)
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
        if outputs.contains(OutputMode::MSA) && !self.graph_tracks_read_ids {
            return Err(Error::InvalidInput(
                "cannot generate MSA from a graph built without read ids; rebuild the graph with read ids enabled before adding sequences"
                    .into(),
            ));
        }
        if outputs.contains(OutputMode::CONSENSUS)
            && self.params.max_consensus() > 1
            && !self.graph_tracks_read_ids
        {
            return Err(Error::InvalidInput(
                "cannot generate multiple consensus sequences from a graph built without read ids; rebuild the graph with read ids enabled before adding sequences"
                    .into(),
            ));
        }

        let alphabet = self.alphabet();
        self.params.set_outputs_for_call(outputs);

        let needs_msa_rank = outputs.contains(OutputMode::MSA)
            || (outputs.contains(OutputMode::CONSENSUS) && self.consensus_needs_msa_rank()?);
        if needs_msa_rank {
            self.ensure_msa_rank_buffer()?;
        }

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
}
