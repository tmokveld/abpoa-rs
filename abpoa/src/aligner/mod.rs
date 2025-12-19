//! AbPOA aligner
//!
//! Module tries to keep the unsafe FFI contained while exposing methods for
//! running alignments, updating graphs, and exporting consensus/MSA/GFA outputs.

mod core;
mod incremental;
mod io;
mod one_shot;
mod raw_alignment;

use crate::encode::{Alphabet, encode_aa, encode_dna};
pub use crate::output::{
    dot::{EdgeLabel, EdgePenWidth, PogDotOptions, RankDir},
    pog::PogWriteOptions,
};
use crate::params::{NodeId, Parameters};
use crate::{Error, Result, sys};
use libc;
use std::{marker::PhantomData, ptr::NonNull, rc::Rc};

pub use incremental::SubgraphRange;
pub use raw_alignment::{CigarOp, GraphCigar, GraphCigarOp, RawAlignment};

/// Wrapper around `abpoa_t`
///
/// `Aligner` is stateful: methods that generate consensus/MSA/GFA outputs may
/// update the underlying `Parameters` (for example enabling quality-weighted scoring based on
/// per-base weights). These changes persist for subsequent operations on the same aligner.
///
/// Output generation is controlled by [`Parameters::set_outputs`] (defaults to generating both
/// consensus and MSA output).
pub struct Aligner {
    raw: NonNull<sys::abpoa_t>,
    params: Parameters,
    graph_tracks_read_ids: bool,
    read_id_count: i32,
    // Not Send/Sync: abPOA uses global mutable tables and allocators without locking.
    _not_send_sync: PhantomData<Rc<()>>,
}

/// Collection of sequences plus optional metadata used as input to MSA calls.
///
/// Construction validates sequence lengths, names, and quality weights.
#[derive(Debug)]
pub struct SequenceBatch<'a> {
    sequences: &'a [&'a [u8]],
    names: Option<&'a [&'a str]>,
    quality_weights: Option<&'a [&'a [i32]]>,
}

impl<'a> SequenceBatch<'a> {
    /// Construct a batch from raw sequence slices.
    pub fn from_sequences(sequences: &'a [&'a [u8]]) -> Result<Self> {
        validate_sequences(sequences)?;
        Ok(Self {
            sequences,
            names: None,
            quality_weights: None,
        })
    }

    /// Attach per-sequence names (must match the sequence count).
    pub fn with_names(mut self, names: &'a [&'a str]) -> Result<Self> {
        validate_names(self.sequences, names)?;
        self.names = Some(names);
        Ok(self)
    }

    /// Attach per-base quality weights (must match sequence count and lengths).
    pub fn with_quality_weights(mut self, quality_weights: &'a [&'a [i32]]) -> Result<Self> {
        validate_quality_weights(self.sequences, quality_weights)?;
        self.quality_weights = Some(quality_weights);
        Ok(self)
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
    fn graph_is_empty(&self) -> bool {
        // Safety: `abg` is owned by the aligner and valid for the lifetime of `self`.
        unsafe {
            let abg = (*self.as_ptr()).abg;
            abg.is_null() || (*abg).node_n <= 2
        }
    }

    fn alphabet(&self) -> Alphabet {
        // Parameters default to DNA, only DNA and amino acids are supported.
        self.params.get_alphabet().unwrap()
    }

    fn bump_read_id_count(&mut self, read_id: i32) -> Result<i32> {
        if read_id < 0 {
            return Err(Error::InvalidInput("read_id must be non-negative".into()));
        }
        let required = read_id.checked_add(1).ok_or(Error::InvalidInput(
            "read_id exceeds supported range".into(),
        ))?;
        if required > self.read_id_count {
            self.read_id_count = required;
        }
        Ok(self.read_id_count)
    }

    fn set_read_id_count(&mut self, total_reads: i32) -> Result<()> {
        if total_reads < 0 {
            return Err(Error::InvalidInput(
                "total_reads must be non-negative".into(),
            ));
        }
        if total_reads < self.read_id_count {
            return Err(Error::InvalidInput(
                format!(
                    "total_reads {total_reads} cannot be smaller than existing sequence count {}",
                    self.read_id_count
                )
                .into(),
            ));
        }
        self.read_id_count = total_reads;
        Ok(())
    }

    fn ensure_sequence_count(&mut self, total_reads: i32) -> Result<()> {
        // Safety: `abs` is allocated alongside the aligner; update the total number of reads the
        // graph accounts for so consensus/MSA generation sizes buffers correctly.
        let abs_ptr = unsafe { (*self.as_mut_ptr()).abs };
        if abs_ptr.is_null() {
            return Err(Error::NullPointer("abpoa sequence container was null"));
        }
        // Safety: `abs_ptr` is owned by this aligner. abPOA reallocates its internal arrays
        // (seq/name/is_rc/etc.) based on `n_seq`.
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
        let required = read_id_offset
            .checked_add(to_i32(seqs.len(), "too many sequences for abpoa")?)
            .ok_or(Error::InvalidInput("too many sequences for abpoa".into()))?;
        if total_reads < required {
            return Err(Error::InvalidInput(
                format!(
                    "total_reads {total_reads} cannot be smaller than required {required}"
                )
                .into(),
            ));
        }
        self.set_read_id_count(total_reads)?;

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
        // Safety: `abg` is allocated alongside the aligner and lives for `'self`.
        let abg = unsafe { (*self.as_mut_ptr()).abg };
        let mut abg =
            NonNull::new(abg).ok_or(Error::NullPointer("abpoa graph pointer was null"))?;
        // Safety: the graph pointer is valid for the duration of the borrow.
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

    fn consensus_needs_msa_rank(&mut self) -> Result<bool> {
        let params_ptr = self.params.as_mut_ptr()?;
        let raw = unsafe { params_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?;
        Ok(raw.max_n_cons > 1 || raw.cons_algrm == sys::ABPOA_MF as i32)
    }

    fn ensure_msa_rank_buffer(&mut self) -> Result<()> {
        let graph = self.graph_mut()?;
        let required = graph
            .index_rank_m
            .max(graph.node_m)
            .max(graph.node_n)
            .max(0) as usize;
        if required == 0 {
            return Ok(());
        }
        let bytes = required
            .checked_mul(std::mem::size_of::<i32>())
            .ok_or(Error::InvalidInput(
                "msa rank buffer size exceeds addressable memory".into(),
            ))?;

        let current = graph.node_id_to_msa_rank as *mut libc::c_void;
        let new_ptr = if current.is_null() {
            // Safety: allocate a buffer sized to at least the graph capacity. abPOA will free it
            // with `free()` as part of `abpoa_free_graph`.
            (unsafe { libc::calloc(required, std::mem::size_of::<i32>()) }) as *mut i32
        } else {
            // Safety: `node_id_to_msa_rank` is allocated with the system allocator (abPOA uses
            // malloc/realloc); resizing preserves a valid allocation for upstream to free.
            (unsafe { libc::realloc(current, bytes) }) as *mut i32
        };
        if new_ptr.is_null() {
            return Err(Error::NullPointer("failed to allocate msa rank buffer"));
        }
        graph.node_id_to_msa_rank = new_ptr;
        Ok(())
    }

    fn invalidate_graph(graph: &mut sys::abpoa_graph_t) {
        graph.set_is_topological_sorted(0);
        graph.set_is_called_cons(0);
        graph.set_is_set_msa_rank(0);
    }

    fn sequence_count(&mut self) -> Result<i32> {
        // Safety: `abs` is owned alongside the aligner and valid for the lifetime of `self`.
        let abs = unsafe { (*self.as_ptr()).abs };
        if abs.is_null() {
            return Err(Error::NullPointer("abpoa sequence container was null"));
        }
        // Safety: we only read a single integer field.
        let count = unsafe { (*abs).n_seq };
        let count = count.max(0);
        self.read_id_count = count;
        Ok(count)
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
        // matching `abpoa_free` for each init.
        unsafe { sys::abpoa_free(self.raw.as_ptr()) }
    }
}

fn encode_sequences(seqs: &[&[u8]], alphabet: Alphabet) -> Vec<Vec<u8>> {
    match alphabet {
        Alphabet::Dna => seqs.iter().map(|seq| encode_dna(seq)).collect(),
        Alphabet::AminoAcid => seqs.iter().map(|seq| encode_aa(seq)).collect(),
    }
}

fn validate_sequences(sequences: &[&[u8]]) -> Result<()> {
    to_i32(sequences.len(), "too many sequences for abpoa")?;
    for seq in sequences {
        if seq.is_empty() {
            return Err(Error::InvalidInput("cannot align an empty sequence".into()));
        }
        to_i32(seq.len(), "sequence length exceeds i32")?;
    }
    Ok(())
}

fn validate_names(sequences: &[&[u8]], names: &[&str]) -> Result<()> {
    if names.len() != sequences.len() {
        return Err(Error::InvalidInput(
            "names length must match sequence count".into(),
        ));
    }
    Ok(())
}

fn validate_quality_weights(sequences: &[&[u8]], weights: &[&[i32]]) -> Result<()> {
    if weights.len() != sequences.len() {
        return Err(Error::InvalidInput(
            "quality weights length must match sequence count".into(),
        ));
    }

    for (row, seq) in weights.iter().zip(sequences) {
        if row.len() != seq.len() {
            return Err(Error::InvalidInput(
                "quality weights must match each sequence length".into(),
            ));
        }
        if row.iter().any(|weight| *weight < 0) {
            return Err(Error::InvalidInput(
                "quality weights must be non-negative".into(),
            ));
        }
    }
    Ok(())
}

fn to_i32(value: usize, context: &'static str) -> Result<i32> {
    i32::try_from(value).map_err(|_| Error::InvalidInput(context.into()))
}

#[cfg(test)]
mod tests;
