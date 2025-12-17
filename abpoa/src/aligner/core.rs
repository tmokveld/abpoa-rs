use super::{Aligner, to_i32};
use crate::graph::Graph;
use crate::params::Parameters;
use crate::{Error, Result, sys};
use std::{marker::PhantomData, ptr::NonNull};

impl Aligner {
    /// Allocate a new aligner handle with default parameters.
    pub fn new() -> Result<Self> {
        Self::with_params(Parameters::new()?)
    }

    /// Allocate a new aligner handle using the provided parameters.
    pub fn with_params(params: Parameters) -> Result<Self> {
        // Safety: calls into the C API to allocate an aligner; abPOA defines the allocation layout
        // and requires `abpoa_free` to release it.
        let raw = unsafe { sys::abpoa_init() };
        let raw = NonNull::new(raw).ok_or(Error::NullPointer("abpoa_init returned null"))?;

        Ok(Self {
            raw,
            params,
            graph_tracks_read_ids: true,
            _not_send_sync: PhantomData,
        })
    }

    pub(crate) fn as_mut_ptr(&mut self) -> *mut sys::abpoa_t {
        self.raw.as_ptr()
    }

    pub(crate) fn as_ptr(&self) -> *const sys::abpoa_t {
        self.raw.as_ptr()
    }

    /// Borrow the underlying parameters immutably.
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

    /// Borrow the underlying graph for read-only inspection.
    pub fn graph(&self) -> Result<Graph<'_>> {
        // Safety: the graph pointer is owned by this aligner and remains live for `'self`
        let abg = unsafe { (*self.as_ptr()).abg };
        let abs = unsafe { (*self.as_ptr()).abs };
        Graph::new(abg, abs)
    }

    /// Reset the underlying graph and DP matrices.
    pub fn reset(&mut self, ref_len: usize) -> Result<()> {
        let ref_len = to_i32(ref_len, "reference length exceeds i32")?;
        let params_ptr = self.params.as_mut_ptr()?;
        // Safety: `raw` and `params` are live pointers; abPOA expects `abpoa_reset` before
        // running incremental alignment.
        unsafe { sys::abpoa_reset(self.as_mut_ptr(), params_ptr, ref_len) };
        // `abpoa_reset` clears the graph; read-id tracking starts fresh.
        self.graph_tracks_read_ids = true;
        Ok(())
    }
}
