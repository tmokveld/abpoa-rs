use crate::encode::{Alphabet, decode_aa_code, decode_dna_code};
use crate::graph::Graph;
use crate::params::NodeId;
use crate::{Error, Result, sys};
use libc;
use std::{marker::PhantomData, ptr, rc::Rc, slice};

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

impl Iterator for GraphCigar<'_> {
    type Item = GraphCigarOp;

    fn next(&mut self) -> Option<Self::Item> {
        // See: abpoa-sys/abPOA/src/abpoa_align.h
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
    pub(super) fn new() -> Self {
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

    pub(super) fn as_mut_ptr(&mut self) -> *mut sys::abpoa_res_t {
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
    /// `Graph` view can be borrowed via [`crate::Aligner::graph`] for use here
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

fn decode_base(alphabet: Alphabet, code: u8) -> char {
    match alphabet {
        Alphabet::Dna => decode_dna_code(code),
        Alphabet::AminoAcid => decode_aa_code(code),
    }
}
