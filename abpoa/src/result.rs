//! Structures representing consensus clusters and MSA output

use crate::encode::{Alphabet, decode_aa, decode_dna};
use crate::params::NodeId;
use crate::sys;
use std::{marker::PhantomData, rc::Rc, slice};

/// Consensus cluster containing a consensus string and per-base metadata
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cluster {
    pub read_ids: Vec<i32>,
    pub consensus: String,
    /// Graph node id corresponding to each consensus base (same length as `consensus`)
    pub node_ids: Vec<NodeId>,
    pub coverage: Vec<i32>,
    /// FASTQ-encoded quality bytes (Phred+33) for each consensus base (same length as `consensus`)
    ///
    /// If per-base input qualities are not provided, abPOA derives these scores from the
    /// per-base consensus coverage and the number of reads in the cluster (see upstream
    /// `abpoa_cons_phred_score`).
    pub phred: Vec<u8>,
}

/// Output of a one-shot MSA run: per-sequence alignments plus consensus clusters
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MsaResult {
    pub msa: Vec<String>,
    pub clusters: Vec<Cluster>,
}

/// Consensus cluster carrying encoded bases instead of decoded strings.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EncodedCluster {
    pub read_ids: Vec<i32>,
    pub consensus: Vec<u8>,
    /// Graph node id corresponding to each consensus base (same length as `consensus`)
    pub node_ids: Vec<NodeId>,
    pub coverage: Vec<i32>,
    /// FASTQ-encoded quality bytes (Phred+33) for each consensus base (same length as `consensus`)
    ///
    /// If per-base input qualities are not provided, abPOA derives these scores from the
    /// per-base consensus coverage and the number of reads in the cluster (see upstream
    /// `abpoa_cons_phred_score`).
    pub phred: Vec<u8>,
}

/// MSA and consensus output in the raw abPOA integer alphabet.
///
/// The `msa` rows and `consensus` sequences are copies of the underlying C
/// buffers (0..m-1 codes for the configured alphabet) and can be decoded
/// using [`crate::encode::decode_dna`] or [`crate::encode::decode_aa`] if
/// needed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EncodedMsaResult {
    pub msa: Vec<Vec<u8>>,
    pub clusters: Vec<EncodedCluster>,
}

/// Zero-copy borrowed view into abPOA's encoded MSA/consensus buffers.
///
/// This view borrows the underlying output buffers owned by an [`crate::Aligner`]. It is
/// invalidated by any subsequent call that regenerates or clears MSA/consensus output (such as
/// `finalize_msa*`, `msa*`, `reset`, or graph restoration). Rust's borrow checker prevents using
/// the aligner mutably while a view exists.
pub struct EncodedMsaView<'a> {
    abc: *const sys::abpoa_cons_t,
    alphabet: Alphabet,
    _owner: PhantomData<&'a sys::abpoa_t>,
    // Not Send/Sync: matches `Aligner`'s thread-safety guarantees.
    _not_send_sync: PhantomData<Rc<()>>,
}

impl<'a> EncodedMsaView<'a> {
    pub(crate) fn new(abc: *const sys::abpoa_cons_t, alphabet: Alphabet) -> Self {
        Self {
            abc,
            alphabet,
            _owner: PhantomData,
            _not_send_sync: PhantomData,
        }
    }

    fn abc(&self) -> Option<&sys::abpoa_cons_t> {
        // Safety: `abc` is either null (empty output) or points to the aligner-owned
        // consensus/MSA buffers, which remain valid for the lifetime of this view.
        unsafe { self.abc.as_ref() }
    }

    /// Alphabet used to interpret the encoded bases.
    pub fn alphabet(&self) -> Alphabet {
        self.alphabet
    }

    /// Number of input sequences included in the MSA output.
    pub fn sequence_count(&self) -> usize {
        self.abc().map(|abc| abc.n_seq.max(0) as usize).unwrap_or(0)
    }

    /// Number of consensus clusters generated.
    pub fn cluster_count(&self) -> usize {
        self.abc()
            .map(|abc| abc.n_cons.max(0) as usize)
            .unwrap_or(0)
    }

    /// Length (columns) of the MSA output.
    pub fn msa_len(&self) -> usize {
        self.abc()
            .map(|abc| abc.msa_len.max(0) as usize)
            .unwrap_or(0)
    }

    /// Borrow a single MSA row for an input sequence by index.
    pub fn msa_row(&self, index: usize) -> Option<&[u8]> {
        let abc = self.abc()?;
        if abc.msa_len <= 0 || abc.msa_base.is_null() {
            return None;
        }
        let n_seq = abc.n_seq.max(0) as usize;
        if index >= n_seq {
            return None;
        }
        let row_ptr = unsafe { *abc.msa_base.add(index) };
        if row_ptr.is_null() {
            return None;
        }
        let msa_len = abc.msa_len.max(0) as usize;
        // Safety: `abpoa_allocate_rc_msa` allocates `msa_len` bytes for each row when
        // `msa_len > 0`. The buffers are owned by the aligner and live for the view lifetime.
        Some(unsafe { slice::from_raw_parts(row_ptr, msa_len) })
    }

    /// Iterate over MSA rows for input sequences (excludes consensus rows).
    pub fn msa_rows(&self) -> EncodedMsaRows<'_> {
        EncodedMsaRows {
            view: self,
            next: 0,
            total: self.sequence_count(),
        }
    }

    /// Borrow a consensus row embedded in the RC-MSA output by cluster index.
    ///
    /// This is only present when MSA and consensus are generated together.
    pub fn consensus_msa_row(&self, cluster_index: usize) -> Option<&[u8]> {
        let abc = self.abc()?;
        if abc.msa_len <= 0 || abc.msa_base.is_null() {
            return None;
        }
        let n_seq = abc.n_seq.max(0) as usize;
        let n_cons = abc.n_cons.max(0) as usize;
        if cluster_index >= n_cons {
            return None;
        }
        let row_ptr = unsafe { *abc.msa_base.add(n_seq + cluster_index) };
        if row_ptr.is_null() {
            return None;
        }
        let msa_len = abc.msa_len.max(0) as usize;
        // Safety: consensus MSA rows are allocated alongside the sequence MSA rows when
        // `abpoa_allocate_rc_msa` is called with `n_cons > 0`.
        Some(unsafe { slice::from_raw_parts(row_ptr, msa_len) })
    }

    /// Borrow a consensus cluster by index.
    pub fn cluster(&self, index: usize) -> Option<EncodedClusterView<'_>> {
        if index < self.cluster_count() {
            Some(EncodedClusterView { view: self, index })
        } else {
            None
        }
    }

    /// Iterate over consensus clusters.
    pub fn clusters(&self) -> EncodedClusters<'_> {
        EncodedClusters {
            view: self,
            next: 0,
            total: self.cluster_count(),
        }
    }
}

/// Iterator over encoded MSA rows.
pub struct EncodedMsaRows<'a> {
    view: &'a EncodedMsaView<'a>,
    next: usize,
    total: usize,
}

impl<'a> Iterator for EncodedMsaRows<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.next < self.total {
            let idx = self.next;
            self.next += 1;
            self.view.msa_row(idx)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.total.saturating_sub(self.next);
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for EncodedMsaRows<'a> {
    fn len(&self) -> usize {
        self.total.saturating_sub(self.next)
    }
}

/// Borrowed view of a single encoded consensus cluster.
pub struct EncodedClusterView<'a> {
    view: &'a EncodedMsaView<'a>,
    index: usize,
}

impl<'a> EncodedClusterView<'a> {
    fn abc(&self) -> Option<&sys::abpoa_cons_t> {
        self.view.abc()
    }

    /// Read ids for this cluster (may be empty when unavailable)
    pub fn read_ids(&self) -> &[i32] {
        let Some(abc) = self.abc() else {
            return &[];
        };
        if abc.clu_n_seq.is_null() || abc.clu_read_ids.is_null() {
            return &[];
        }
        let count = unsafe { *abc.clu_n_seq.add(self.index) }.max(0) as usize;
        if count == 0 {
            return &[];
        }
        let ids_ptr = unsafe { *abc.clu_read_ids.add(self.index) };
        if ids_ptr.is_null() {
            return &[];
        }
        // Safety: `abpoa_allocate_cons` allocates `n_seq` entries for each cluster; abPOA fills
        // `clu_n_seq` for the actual count.
        unsafe { slice::from_raw_parts(ids_ptr, count) }
    }

    /// Length of the consensus sequence for this cluster
    pub fn consensus_len(&self) -> usize {
        let Some(abc) = self.abc() else {
            return 0;
        };
        if abc.cons_len.is_null() {
            return 0;
        }
        unsafe { *abc.cons_len.add(self.index) }.max(0) as usize
    }

    /// Encoded consensus sequence for this cluster
    pub fn consensus(&self) -> &[u8] {
        let Some(abc) = self.abc() else {
            return &[];
        };
        let len = self.consensus_len();
        if len == 0 || abc.cons_base.is_null() {
            return &[];
        }
        let ptr = unsafe { *abc.cons_base.add(self.index) };
        if ptr.is_null() {
            return &[];
        }
        // Safety: `abpoa_allocate_cons` allocates `n_node` bytes per cluster consensus; abPOA
        // sets `cons_len[index]` bases, and the buffers live for the view lifetime.
        unsafe { slice::from_raw_parts(ptr, len) }
    }

    /// Graph node ids corresponding to each consensus base (raw `i32` ids)
    pub fn node_ids_raw(&self) -> &[i32] {
        let Some(abc) = self.abc() else {
            return &[];
        };
        let len = self.consensus_len();
        if len == 0 || abc.cons_node_ids.is_null() {
            return &[];
        }
        let ptr = unsafe { *abc.cons_node_ids.add(self.index) };
        if ptr.is_null() {
            return &[];
        }
        // Safety: `abpoa_allocate_cons` allocates `n_node` ids per cluster consensus.
        unsafe { slice::from_raw_parts(ptr, len) }
    }

    /// Coverage per consensus base
    pub fn coverage(&self) -> &[i32] {
        let Some(abc) = self.abc() else {
            return &[];
        };
        let len = self.consensus_len();
        if len == 0 || abc.cons_cov.is_null() {
            return &[];
        }
        let ptr = unsafe { *abc.cons_cov.add(self.index) };
        if ptr.is_null() {
            return &[];
        }
        // Safety: coverage arrays are allocated by `abpoa_allocate_cons` for each cluster
        unsafe { slice::from_raw_parts(ptr, len) }
    }

    /// Raw per-base Phred scores (ASCII codes, Phred+33) for consensus FASTQ output
    pub fn phred_scores_raw(&self) -> &[i32] {
        let Some(abc) = self.abc() else {
            return &[];
        };
        let len = self.consensus_len();
        if len == 0 || abc.cons_phred_score.is_null() {
            return &[];
        }
        let ptr = unsafe { *abc.cons_phred_score.add(self.index) };
        if ptr.is_null() {
            return &[];
        }
        // Safety: phred arrays are allocated by `abpoa_allocate_cons` for each cluster
        unsafe { slice::from_raw_parts(ptr, len) }
    }
}

/// Iterator over consensus clusters in an [`EncodedMsaView`].
pub struct EncodedClusters<'a> {
    view: &'a EncodedMsaView<'a>,
    next: usize,
    total: usize,
}

impl<'a> Iterator for EncodedClusters<'a> {
    type Item = EncodedClusterView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next < self.total {
            let idx = self.next;
            self.next += 1;
            Some(EncodedClusterView {
                view: self.view,
                index: idx,
            })
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.total.saturating_sub(self.next);
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for EncodedClusters<'a> {
    fn len(&self) -> usize {
        self.total.saturating_sub(self.next)
    }
}

/// Placeholder for a single alignment row (reserved for future APIs)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub sequence: String,
}

impl MsaResult {
    /// # Safety
    /// `abc` must point to a valid `abpoa_cons_t` populated by abPOA; the memory must outlive
    /// this call and remain readable for the duration of the copy into Rust-owned buffers
    pub(crate) unsafe fn from_raw(abc: *const sys::abpoa_cons_t, alphabet: Alphabet) -> Self {
        let decode_row: fn(&[u8]) -> String = match alphabet {
            Alphabet::Dna => decode_dna,
            Alphabet::AminoAcid => decode_aa,
        };

        let parsed = unsafe { parse_raw(abc) };
        let msa = parsed.msa.iter().map(|row| decode_row(row)).collect();
        let clusters = parsed
            .clusters
            .into_iter()
            .map(|cluster| Cluster {
                read_ids: cluster.read_ids,
                consensus: decode_row(&cluster.consensus),
                node_ids: cluster.node_ids,
                coverage: cluster.coverage,
                phred: cluster.phred,
            })
            .collect();

        Self { msa, clusters }
    }
}

impl EncodedMsaResult {
    /// # Safety
    /// `abc` must point to a valid `abpoa_cons_t` populated by abPOA; the memory must outlive
    /// this call and remain readable for the duration of the copy into Rust-owned buffers
    pub(crate) unsafe fn from_raw(abc: *const sys::abpoa_cons_t) -> Self {
        let parsed = unsafe { parse_raw(abc) };
        let msa = parsed.msa;
        let clusters = parsed
            .clusters
            .into_iter()
            .map(|cluster| EncodedCluster {
                read_ids: cluster.read_ids,
                consensus: cluster.consensus,
                node_ids: cluster.node_ids,
                coverage: cluster.coverage,
                phred: cluster.phred,
            })
            .collect();

        Self { msa, clusters }
    }
}

struct RawCluster {
    read_ids: Vec<i32>,
    consensus: Vec<u8>,
    node_ids: Vec<NodeId>,
    coverage: Vec<i32>,
    phred: Vec<u8>,
}

struct RawParsedResult {
    msa: Vec<Vec<u8>>,
    clusters: Vec<RawCluster>,
}

/// Parse MSA and consensus output from an abPOA consensus struct.
///
/// # Safety
/// `abc` must point to a valid `abpoa_cons_t` populated by abPOA. All internal pointers must be
/// readable for the duration of the copies performed here.
unsafe fn parse_raw(abc: *const sys::abpoa_cons_t) -> RawParsedResult {
    let Some(abc) = (unsafe { abc.as_ref() }) else {
        return RawParsedResult {
            msa: Vec::new(),
            clusters: Vec::new(),
        };
    };

    let n_seq = abc.n_seq.max(0) as usize;
    let msa_len = abc.msa_len.max(0) as usize;
    let msa = if msa_len == 0 || abc.msa_base.is_null() {
        Vec::new()
    } else {
        (0..n_seq)
            .map(|i| {
                let row_ptr = unsafe { *abc.msa_base.add(i) };
                if row_ptr.is_null() {
                    Vec::new()
                } else {
                    unsafe { slice::from_raw_parts(row_ptr, msa_len) }.to_vec()
                }
            })
            .collect()
    };

    let n_cons = abc.n_cons.max(0) as usize;
    let clusters = if n_cons == 0 || abc.cons_len.is_null() || abc.cons_base.is_null() {
        Vec::new()
    } else {
        (0..n_cons)
            .map(|cons_idx| {
                let len = unsafe { *abc.cons_len.add(cons_idx) }.max(0) as usize;
                let consensus_ptr = unsafe { *abc.cons_base.add(cons_idx) };
                let consensus = if len == 0 || consensus_ptr.is_null() {
                    Vec::new()
                } else {
                    unsafe { slice::from_raw_parts(consensus_ptr, len) }.to_vec()
                };

                let coverage = if abc.cons_cov.is_null() {
                    Vec::new()
                } else {
                    let cov_ptr = unsafe { *abc.cons_cov.add(cons_idx) };
                    if len == 0 || cov_ptr.is_null() {
                        Vec::new()
                    } else {
                        unsafe { slice::from_raw_parts(cov_ptr, len) }.to_vec()
                    }
                };

                let node_ids = if abc.cons_node_ids.is_null() {
                    Vec::new()
                } else {
                    let ids_ptr = unsafe { *abc.cons_node_ids.add(cons_idx) };
                    if len == 0 || ids_ptr.is_null() {
                        Vec::new()
                    } else {
                        unsafe { slice::from_raw_parts(ids_ptr, len) }
                            .iter()
                            .copied()
                            .map(NodeId)
                            .collect()
                    }
                };

                let phred = if abc.cons_phred_score.is_null() {
                    Vec::new()
                } else {
                    let phred_ptr = unsafe { *abc.cons_phred_score.add(cons_idx) };
                    if len == 0 || phred_ptr.is_null() {
                        Vec::new()
                    } else {
                        unsafe { slice::from_raw_parts(phred_ptr, len) }
                            .iter()
                            .copied()
                            .map(|score| u8::try_from(score).unwrap_or(b'!'))
                            .collect()
                    }
                };

                let read_ids = if abc.clu_n_seq.is_null() || abc.clu_read_ids.is_null() {
                    Vec::new()
                } else {
                    let count = unsafe { *abc.clu_n_seq.add(cons_idx) }.max(0) as usize;
                    if count == 0 {
                        Vec::new()
                    } else {
                        let ids_ptr = unsafe { *abc.clu_read_ids.add(cons_idx) };
                        if ids_ptr.is_null() {
                            Vec::new()
                        } else {
                            unsafe { slice::from_raw_parts(ids_ptr, count) }.to_vec()
                        }
                    }
                };

                RawCluster {
                    read_ids,
                    consensus,
                    node_ids,
                    coverage,
                    phred,
                }
            })
            .collect()
    };

    RawParsedResult { msa, clusters }
}
