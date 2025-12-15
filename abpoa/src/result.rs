//! Structures representing consensus clusters and MSA output

use crate::encode::{Alphabet, decode_aa, decode_dna};
use crate::params::NodeId;
use crate::sys;
use std::slice;

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
