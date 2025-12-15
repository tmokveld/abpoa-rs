//! Helpers for exporting the partial order graph (POG) to files.

use crate::output::dot::PogDotOptions;

/// How to export a partial order graph (POG) when writing to a file path.
#[derive(Debug, Clone, PartialEq)]
pub enum PogWriteOptions {
    /// Use the Rust DOT emitter and invoke GraphViz `dot` directly (no shell).
    ///
    /// When writing `png`/`pdf`, a DOT file is written next to the output at `path + ".dot"`
    /// (mirroring upstream).
    Rust(PogDotOptions),

    /// Use upstream `abpoa_dump_pog`.
    ///
    /// Note: upstream calls `system("dot ...")` and will `exit()` on failure, so prefer
    /// [`PogWriteOptions::Rust`].
    Upstream,
}

impl Default for PogWriteOptions {
    fn default() -> Self {
        Self::Rust(PogDotOptions::default())
    }
}

impl From<PogDotOptions> for PogWriteOptions {
    fn from(value: PogDotOptions) -> Self {
        Self::Rust(value)
    }
}

impl From<&PogDotOptions> for PogWriteOptions {
    fn from(value: &PogDotOptions) -> Self {
        Self::Rust(value.clone())
    }
}
