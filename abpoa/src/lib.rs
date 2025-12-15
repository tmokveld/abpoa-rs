//! Safe Rust wrapper around the abPOA partial order aligner

//! ## Platform support
//!
//! `abpoa` currently supports Unix-like targets only

#[cfg(not(unix))]
compile_error!(
    "`abpoa` currently supports Unix-like targets only; upstream abPOA depends on POSIX APIs."
);

#[cfg(unix)]
use std::slice;

#[cfg(unix)]
pub use abpoa_sys as sys;

#[cfg(unix)]
pub mod aligner;
#[cfg(unix)]
pub mod encode;
#[cfg(unix)]
pub mod graph;
#[cfg(unix)]
pub mod output;
#[cfg(unix)]
pub mod params;
#[cfg(unix)]
pub mod result;

#[cfg(unix)]
pub use crate::aligner::{
    Aligner, CigarOp, GraphCigar, GraphCigarOp, RawAlignment, SequenceBatch, SubgraphRange,
};
#[cfg(unix)]
pub use crate::encode::Alphabet;
#[cfg(unix)]
pub use crate::graph::{Graph, NodeRef, NodeView, Sequence, SequenceIter, SequenceStr, Sequences};
#[cfg(unix)]
pub use crate::output::dot::{EdgeLabel, EdgePenWidth, PogDotOptions, RankDir};
#[cfg(unix)]
pub use crate::output::pog::PogWriteOptions;
#[cfg(unix)]
pub use crate::params::{
    AlignMode, ConsensusAlgorithm, GapMode, GapPenalty, NodeId, OutputMode, Parameters, Scoring,
    SentinelNode, Verbosity,
};
#[cfg(unix)]
pub use crate::result::{Alignment, Cluster, EncodedCluster, EncodedMsaResult, MsaResult};

#[cfg(unix)]
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("abPOA call {func} failed with code {code}")]
    Abpoa { func: &'static str, code: i32 },

    #[error("invalid input: {0}")]
    InvalidInput(std::borrow::Cow<'static, str>),

    #[error("null pointer: {0}")]
    NullPointer(&'static str),

    #[error(transparent)]
    Io(#[from] std::io::Error),
}

#[cfg(unix)]
pub type Result<T> = std::result::Result<T, Error>;

/// Shared view of abPOA's precomputed integer log2 table for 16-bit values
#[cfg(unix)]
pub fn log_table_65536() -> &'static [u8] {
    // Safety: `AB_LOG_TABLE_65536` points to static storage sized to 2^16 entries and lives
    // for the duration of the process
    unsafe {
        let ptr = sys::AB_LOG_TABLE_65536.as_ptr();
        slice::from_raw_parts(ptr, 65_536)
    }
}

/// Shared view of abPOA's precomputed bit-count table for 16-bit values
#[cfg(unix)]
pub fn bit_table_16() -> &'static [u8] {
    // Safety: `AB_BIT_TABLE_16` points to static storage sized to 2^16 entries and lives
    // for the duration of the process
    unsafe {
        let ptr = sys::AB_BIT_TABLE_16.as_ptr();
        slice::from_raw_parts(ptr, 65_536)
    }
}
