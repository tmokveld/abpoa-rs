//! Parameter ownership and configuration access to `abpoa_para_t`

use libc;
use std::{
    ffi::CString,
    marker::PhantomData,
    mem,
    path::Path,
    ptr::{self, NonNull},
    rc::Rc,
};

use crate::{Error, Result, encode::Alphabet, sys};

/// Wrapper around `abpoa_para_t`
pub struct Parameters {
    raw: NonNull<sys::abpoa_para_t>,
    custom_matrix: Option<CustomMatrix>,
    alphabet: Alphabet,
    preserve_read_ids: bool,
    dirty: bool,
    // Not Send/Sync: abPOA mutates global decoder tables during parameter setup
    _not_send_sync: PhantomData<Rc<()>>,
}

#[derive(Clone)]
struct CustomMatrix {
    m: i32,
    values: Vec<i32>,
}

bitflags::bitflags! {
    /// Outputs to request when generating MSA/consensus results
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct OutputMode: u32 {
        const CONSENSUS = 1 << 0;
        const MSA       = 1 << 1;
    }
}

impl Default for OutputMode {
    fn default() -> Self {
        Self::CONSENSUS | Self::MSA
    }
}

impl Parameters {
    /// Allocate a new parameter object with the default abPOA configuration
    pub fn new() -> Result<Self> {
        let mut params = Self::allocate()?;
        params.finalize();
        Ok(params)
    }

    /// Start configuring a parameter object in-place
    ///
    /// This is equivalent to `Parameters::new()` but does not eagerly run finalization,
    /// changes are finalized automatically on the next FFI call
    pub fn configure() -> Result<Self> {
        Self::allocate()
    }

    /// Allocate a new parameter object
    fn allocate() -> Result<Self> {
        // Safety: calls into the C API to allocate a parameter object; abPOA owns the allocation
        // layout and expects `abpoa_free_para` to release it
        let raw = unsafe { sys::abpoa_init_para() };
        let raw = NonNull::new(raw).ok_or(Error::NullPointer("abpoa_init_para returned null"))?;

        let mut params = Self {
            raw,
            custom_matrix: None,
            alphabet: Alphabet::Dna,
            preserve_read_ids: true,
            dirty: true,
            _not_send_sync: PhantomData,
        };
        params.set_alphabet(Alphabet::Dna)?;
        params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
        params.set_use_read_ids(true);
        Ok(params)
    }

    pub(crate) fn as_mut_ptr(&mut self) -> *mut sys::abpoa_para_t {
        self.finalize();
        self.raw.as_ptr()
    }

    fn mark_dirty(&mut self) {
        self.dirty = true;
    }

    fn finalize(&mut self) {
        if !self.dirty {
            return;
        }
        let needs_output_tables = unsafe {
            // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
            let raw = self.raw.as_ref();
            raw.out_msa() != 0
                || raw.out_gfa() != 0
                || raw.max_n_cons > 1
                || raw.cons_algrm == sys::ABPOA_MF as i32
        };

        if needs_output_tables {
            crate::runtime::ensure_output_tables();
        }

        crate::runtime::with_abpoa_global_write_lock(|| unsafe {
            // Safety: we hold the crate-global lock to avoid concurrent mutation of abPOA's
            // global tables from `abpoa_post_set_para`.
            let raw = self.raw.as_mut();

            let out_msa = raw.out_msa();
            let out_gfa = raw.out_gfa();
            let max_n_cons = raw.max_n_cons;
            let cons_algrm = raw.cons_algrm;

            if needs_output_tables {
                // Avoid rewriting global output lookup tables in `abpoa_post_set_para`.
                raw.set_out_msa(0);
                raw.set_out_gfa(0);
                raw.max_n_cons = 1;
                raw.cons_algrm = sys::ABPOA_HB as i32;
            }

            sys::abpoa_post_set_para(self.raw.as_ptr());

            if needs_output_tables {
                raw.set_out_msa(out_msa);
                raw.set_out_gfa(out_gfa);
                raw.max_n_cons = max_n_cons;
                raw.cons_algrm = cons_algrm;
            }

            raw.set_use_read_ids((self.preserve_read_ids || needs_output_tables) as u8);
        });
        self.apply_custom_matrix();
        self.dirty = false;
    }

    /// Configure the alphabet used for encoding sequences
    pub fn set_alphabet(&mut self, alphabet: Alphabet) -> Result<()> {
        let m = alphabet_size(alphabet);

        let current_m = unsafe { self.raw.as_ref().m };
        if current_m != m {
            self.resize_matrix(m)?;
        }
        if let Some(custom) = &self.custom_matrix {
            if custom.m != m {
                self.custom_matrix = None;
            }
        }

        self.alphabet = alphabet;
        self.mark_dirty();
        Ok(())
    }

    fn resize_matrix(&mut self, m: i32) -> Result<()> {
        if m <= 0 {
            return Err(Error::InvalidInput("alphabet size must be positive".into()));
        }
        let new_size = (m as usize)
            .checked_mul(m as usize)
            .and_then(|v| v.checked_mul(mem::size_of::<i32>()))
            .ok_or(Error::InvalidInput("score matrix size overflow".into()))?;

        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            let raw = self.raw.as_mut();
            let new_ptr = libc::realloc(raw.mat as *mut libc::c_void, new_size) as *mut i32;
            if new_ptr.is_null() {
                return Err(Error::NullPointer(
                    "failed to resize scoring matrix for alphabet",
                ));
            }
            raw.mat = new_ptr;
            raw.m = m;
        }

        Ok(())
    }

    /// Current alphabet, if it matches a supported encoding
    pub fn get_alphabet(&self) -> Option<Alphabet> {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        let raw = unsafe { self.raw.as_ref() };
        if raw.m == alphabet_size(self.alphabet) {
            Some(self.alphabet)
        } else {
            None
        }
    }

    /// Set whether quality-weighted scoring is enabled.
    pub fn set_use_quality(&mut self, enabled: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_use_qv(enabled as u8);
        }
        self
    }

    pub(crate) fn set_outputs_for_call(&mut self, outputs: OutputMode) {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            let raw = self.raw.as_mut();
            raw.set_out_cons(outputs.contains(OutputMode::CONSENSUS) as u8);
            raw.set_out_msa(outputs.contains(OutputMode::MSA) as u8);
        }
        let needs_read_ids = unsafe {
            // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
            let raw = self.raw.as_ref();
            raw.out_msa() != 0
                || raw.out_gfa() != 0
                || raw.max_n_cons > 1
                || raw.cons_algrm == sys::ABPOA_MF as i32
        };
        let use_read_ids = self.preserve_read_ids || needs_read_ids;
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_use_read_ids(use_read_ids as u8);
        }
    }

    pub fn set_use_read_ids(&mut self, enabled: bool) -> &mut Self {
        self.preserve_read_ids = enabled;
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_use_read_ids(enabled as u8);
        }
        self.mark_dirty();
        self
    }

    /// Configure minimizer seeding parameters and enable seeding
    ///
    /// - `k`: minimizer k-mer size
    /// - `w`: minimizer window size (one minimizer picked per `w` consecutive k-mers)
    /// - `min_w`: minimum POA window length when partitioning by minimizer anchors
    pub fn set_minimizer_seeding(&mut self, k: i32, w: i32, min_w: i32) -> Result<&mut Self> {
        if k <= 0 || w <= 0 || min_w < 0 {
            return Err(Error::InvalidInput(
                "seeding parameters must be positive (min_w may be zero)".into(),
            ));
        }

        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            let raw = self.raw.as_mut();
            raw.k = k;
            raw.w = w;
            raw.min_w = min_w;
            raw.set_disable_seeding(0);
        }

        self.mark_dirty();
        Ok(self)
    }

    /// Disable or enable minimizer seeding
    pub fn set_disable_seeding(&mut self, disable: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_disable_seeding(disable as u8);
        }
        self.mark_dirty();
        self
    }

    /// Enable progressive partial-order alignment
    pub fn set_progressive_poa(&mut self, enabled: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_progressive_poa(enabled as u8);
        }
        self.mark_dirty();
        self
    }

    /// Adjust path scoring to factor in node coverage
    pub fn set_inc_path_score(&mut self, enabled: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().inc_path_score = enabled as i32;
        }
        self.mark_dirty();
        self
    }

    /// Enable subgraph-aware alignment and consensus counting
    pub fn set_sub_alignment(&mut self, enabled: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_sub_aln(enabled as u8);
        }
        self.mark_dirty();
        self
    }

    /// Sort input reads by length before alignment
    ///
    /// Deprecated: upstream abPOA only applies this when reading sequences from a file
    /// (`abpoa_msa1`), which this crate does not expose. It has no effect on
    /// [`crate::Aligner::msa`], [`crate::Aligner::msa_in_place`], or [`crate::Aligner::add_sequences`]
    ///
    /// Sort your input sequences yourself before calling alignment methods
    #[deprecated(
        since = "0.1.0",
        note = "No effect in the Rust APIs; upstream only applies this in its file-reading path (abpoa_msa1), which is not exposed. Sort inputs yourself."
    )]
    pub fn set_sort_input_seq(&mut self, enabled: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().sort_input_seq = enabled as i32;
        }
        self.mark_dirty();
        self
    }

    /// Allow abPOA to try the reverse complement when the forward strand scores poorly
    pub fn set_ambiguous_strand(&mut self, enabled: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_amb_strand(enabled as u8);
        }
        self.mark_dirty();
        self
    }

    /// When ambiguous, place gaps on the right-most position (wfa2-like)
    pub fn set_gap_on_right(&mut self, enabled: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_put_gap_on_right(enabled as u8);
        }
        self.mark_dirty();
        self
    }

    /// Push gap placement toward the end of the alignment
    pub fn set_gap_at_end(&mut self, enabled: bool) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().set_put_gap_at_end(enabled as u8);
        }
        self.mark_dirty();
        self
    }

    /// Set the alignment mode (global/local/extend)
    pub fn set_align_mode(&mut self, mode: AlignMode) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().align_mode = mode.as_raw();
        }
        self.mark_dirty();
        self
    }

    /// Set the scoring scheme for the alignment
    pub fn set_scoring_scheme(&mut self, scoring: Scoring) -> Result<&mut Self> {
        scoring.validate()?;

        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            let raw = self.raw.as_mut();
            raw.match_ = scoring.match_score;
            raw.mismatch = scoring.mismatch;
            raw.gap_mode = scoring.gaps.mode().as_raw();
            match scoring.gaps {
                GapPenalty::Linear { gap } => {
                    raw.gap_open1 = 0;
                    raw.gap_ext1 = gap;
                    raw.gap_open2 = 0;
                    raw.gap_ext2 = 0;
                }
                GapPenalty::Affine {
                    gap_open,
                    gap_extend,
                } => {
                    raw.gap_open1 = gap_open;
                    raw.gap_ext1 = gap_extend;
                    raw.gap_open2 = 0;
                    raw.gap_ext2 = 0;
                }
                GapPenalty::Convex {
                    gap_open,
                    gap_extend,
                    second_gap_open,
                    second_gap_extend,
                } => {
                    raw.gap_open1 = gap_open;
                    raw.gap_ext1 = gap_extend;
                    raw.gap_open2 = second_gap_open;
                    raw.gap_ext2 = second_gap_extend;
                }
            }
        }
        self.mark_dirty();
        Ok(self)
    }

    /// Configure adaptive band parameters
    pub fn set_band(&mut self, extra_b: i32, extra_f: f32) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            let raw = self.raw.as_mut();
            raw.wb = extra_b;
            raw.wf = extra_f;
        }
        self.mark_dirty();
        self
    }

    /// Configure the Z-drop score in extension alignment.
    ///
    /// `None` disables Z-drop (maps to `-1`, the upstream default).
    pub fn set_zdrop(&mut self, zdrop: Option<i32>) -> Result<&mut Self> {
        let zdrop = match zdrop {
            None => -1,
            Some(value) if value > 0 => value,
            Some(_) => {
                return Err(Error::InvalidInput(
                    "zdrop must be positive when set".into(),
                ));
            }
        };

        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().zdrop = zdrop;
        }
        self.mark_dirty();
        Ok(self)
    }

    /// Configure the end bonus in extension alignment.
    ///
    /// `None` disables the end bonus (maps to `-1`, the upstream default).
    pub fn set_end_bonus(&mut self, end_bonus: Option<i32>) -> Result<&mut Self> {
        let end_bonus = match end_bonus {
            None => -1,
            Some(value) if value > 0 => value,
            Some(_) => {
                return Err(Error::InvalidInput(
                    "end_bonus must be positive when set".into(),
                ));
            }
        };

        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().end_bonus = end_bonus;
        }
        self.mark_dirty();
        Ok(self)
    }

    /// Configure which outputs abPOA should generate when asked to finalize
    pub fn set_outputs(&mut self, outputs: OutputMode) -> OutputMode {
        let previous = self.outputs();
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            let raw = self.raw.as_mut();
            raw.set_out_cons(outputs.contains(OutputMode::CONSENSUS) as u8);
            raw.set_out_msa(outputs.contains(OutputMode::MSA) as u8);
        }
        self.mark_dirty();
        previous
    }

    pub fn outputs(&self) -> OutputMode {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        let raw = unsafe { self.raw.as_ref() };
        let mut output = OutputMode::empty();
        output.set(OutputMode::CONSENSUS, raw.out_cons() != 0);
        output.set(OutputMode::MSA, raw.out_msa() != 0);
        output
    }

    /// Maximum number of consensus sequences configured for this parameter set
    pub fn max_consensus(&self) -> i32 {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe { self.raw.as_ref() }.max_n_cons
    }

    /// Set the maximum number of consensus sequences to emit
    pub fn set_max_consensus(&mut self, max_n_cons: i32) -> Result<&mut Self> {
        if max_n_cons < 1 {
            return Err(Error::InvalidInput(
                "max_n_cons must be positive to request consensus output".into(),
            ));
        }
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().max_n_cons = max_n_cons;
        }
        self.mark_dirty();
        Ok(self)
    }

    /// Configure the minimum cluster frequency threshold
    pub fn set_min_cluster_freq(&mut self, min_freq: f64) -> Result<&mut Self> {
        if !(0.0..=1.0).contains(&min_freq) {
            return Err(Error::InvalidInput(
                "min_freq must be between 0.0 and 1.0 (inclusive)".into(),
            ));
        }
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().min_freq = min_freq;
        }
        self.mark_dirty();
        Ok(self)
    }

    /// Configure consensus calling strategy
    ///
    /// - `HeaviestBundle` walks the graph using edge weights (quality scores if
    ///   `set_use_quality(true)` was set)
    /// - `MostFrequent` performs per-column majority voting, respecting
    ///   `set_sub_alignment(true)` to count spanning reads
    ///
    /// `max_n_cons` and `min_freq` control clustering when you want multiple consensus
    /// sequences and apply to both algorithms
    pub fn set_consensus(
        &mut self,
        alg: ConsensusAlgorithm,
        max_n_cons: i32,
        min_freq: f64,
    ) -> Result<&mut Self> {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            let raw = self.raw.as_mut();
            raw.cons_algrm = alg.as_raw();
        }
        self.set_max_consensus(max_n_cons)?;
        self.set_min_cluster_freq(min_freq)?;
        self.mark_dirty();
        Ok(self)
    }

    /// Control verbosity level
    pub fn set_verbosity(&mut self, level: Verbosity) -> &mut Self {
        // Safety: `raw` is uniquely owned and points to a live `abpoa_para_t`
        unsafe {
            self.raw.as_mut().verbose = level.as_raw();
        }
        self.mark_dirty();
        self
    }

    fn apply_custom_matrix(&mut self) {
        if let Some(custom) = &self.custom_matrix {
            debug_assert_eq!(
                custom.values.len(),
                (custom.m as usize) * (custom.m as usize)
            );
            // Safety: the matrix buffer was resized when the custom matrix was recorded
            unsafe {
                let raw = self.raw.as_mut();
                ptr::copy_nonoverlapping(custom.values.as_ptr(), raw.mat, custom.values.len());
                let (max_mat, min_mis) = matrix_extents(&custom.values);
                raw.max_mat = max_mat;
                raw.min_mis = min_mis;
            }
        }
    }

    /// Load a scoring matrix from a file and use it instead of match/mismatch scores
    pub fn set_score_matrix_file<P: AsRef<Path>>(&mut self, path: P) -> Result<&mut Self> {
        let path = path.as_ref().to_str().ok_or(Error::InvalidInput(
            "score matrix path must be valid UTF-8".into(),
        ))?;
        let c_path = CString::new(path).map_err(|_| {
            Error::InvalidInput("score matrix path cannot contain null bytes".into())
        })?;

        // Safety: `raw` is uniquely owned; replace any previous path to avoid leaking C-allocated
        // storage
        unsafe {
            let raw = self.raw.as_mut();
            if !raw.mat_fn.is_null() {
                libc::free(raw.mat_fn as *mut libc::c_void);
                raw.mat_fn = ptr::null_mut();
            }
            raw.mat_fn = libc::strdup(c_path.as_ptr());
            if raw.mat_fn.is_null() {
                return Err(Error::NullPointer(
                    "failed to store score matrix path for abPOA",
                ));
            }
            raw.use_score_matrix = 1;
        }
        self.custom_matrix = None;
        self.mark_dirty();

        Ok(self)
    }

    /// Provide an in-memory scoring matrix, overrides the matrix generated in `finalize`
    pub fn set_custom_matrix(&mut self, m: i32, matrix: &[i32]) -> Result<&mut Self> {
        if m <= 0 {
            return Err(Error::InvalidInput(
                "custom matrix dimension must be positive".into(),
            ));
        }
        let matrix_alphabet = match m {
            5 => Alphabet::Dna,
            27 => Alphabet::AminoAcid,
            _ => {
                return Err(Error::InvalidInput(
                    "custom matrix dimension must be 5 (DNA) or 27 (amino acids)".into(),
                ));
            }
        };
        if matrix_alphabet != self.alphabet {
            return Err(Error::InvalidInput(
                "custom matrix dimension must match configured alphabet; call set_alphabet first"
                    .into(),
            ));
        }
        let expected = (m as usize)
            .checked_mul(m as usize)
            .ok_or(Error::InvalidInput(
                "custom matrix dimensions overflow native size".into(),
            ))?;
        if matrix.len() != expected {
            return Err(Error::InvalidInput(
                "custom matrix length must be m*m entries".into(),
            ));
        }

        self.resize_matrix(m)?;
        self.custom_matrix = Some(CustomMatrix {
            m,
            values: matrix.to_vec(),
        });
        self.apply_custom_matrix();

        // Safety: `raw` is uniquely owned; ensure we do not try to read a matrix from disk
        unsafe {
            let raw = self.raw.as_mut();
            raw.use_score_matrix = 0;
            if !raw.mat_fn.is_null() {
                libc::free(raw.mat_fn as *mut libc::c_void);
                raw.mat_fn = ptr::null_mut();
            }
        }

        self.mark_dirty();
        Ok(self)
    }

    /// Provide a path to an existing graph for incremental alignment/restoration
    pub fn set_incremental_graph_file<P: AsRef<Path>>(&mut self, path: P) -> Result<&mut Self> {
        let path = path.as_ref().to_str().ok_or(Error::InvalidInput(
            "incremental graph path must be valid UTF-8".into(),
        ))?;
        let c_path = CString::new(path).map_err(|_| {
            Error::InvalidInput("incremental graph path cannot contain null bytes".into())
        })?;

        // Safety: `raw` is uniquely owned; replace any previous path to avoid leaking C-allocated
        // storage
        unsafe {
            let raw = self.raw.as_mut();
            if !raw.incr_fn.is_null() {
                libc::free(raw.incr_fn as *mut libc::c_void);
                raw.incr_fn = std::ptr::null_mut();
            }
            raw.incr_fn = libc::strdup(c_path.as_ptr());
            if raw.incr_fn.is_null() {
                return Err(Error::NullPointer(
                    "failed to store incremental graph path for abPOA",
                ));
            }
        }

        self.mark_dirty();
        Ok(self)
    }

    /// Clear any configured incremental graph path
    pub fn clear_incremental_graph(&mut self) -> &mut Self {
        // Safety: `raw` is uniquely owned; freeing here avoids later double frees because we set
        // the pointer back to null
        unsafe {
            let raw = self.raw.as_mut();
            if !raw.incr_fn.is_null() {
                libc::free(raw.incr_fn as *mut libc::c_void);
                raw.incr_fn = std::ptr::null_mut();
            }
        }
        self.mark_dirty();
        self
    }
}

impl Drop for Parameters {
    fn drop(&mut self) {
        // Safety: `raw` came from `abpoa_init_para` and is exclusively owned here; abPOA expects a
        // matching `abpoa_free_para` for each init
        unsafe { sys::abpoa_free_para(self.raw.as_ptr()) }
    }
}

fn matrix_extents(entries: &[i32]) -> (i32, i32) {
    if entries.is_empty() {
        return (0, 0);
    }

    let mut max_mat = i32::MIN;
    let mut min_mis = 0;
    for &value in entries {
        if value > max_mat {
            max_mat = value;
        }
        let neg = value.saturating_neg();
        if neg > min_mis {
            min_mis = neg;
        }
    }

    (max_mat, min_mis)
}

fn alphabet_size(alphabet: Alphabet) -> i32 {
    match alphabet {
        Alphabet::Dna => 5,
        Alphabet::AminoAcid => 27,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verbosity {
    None,
    Info,
    Debug,
    LongDebug,
}

impl Verbosity {
    pub fn as_raw(self) -> i32 {
        match self {
            Verbosity::None => sys::ABPOA_NONE_VERBOSE as i32,
            Verbosity::Info => sys::ABPOA_INFO_VERBOSE as i32,
            Verbosity::Debug => sys::ABPOA_DEBUG_VERBOSE as i32,
            Verbosity::LongDebug => sys::ABPOA_LONG_DEBUG_VERBOSE as i32,
        }
    }

    pub fn from_raw(value: i32) -> Option<Self> {
        match value {
            v if v == sys::ABPOA_NONE_VERBOSE as i32 => Some(Verbosity::None),
            v if v == sys::ABPOA_INFO_VERBOSE as i32 => Some(Verbosity::Info),
            v if v == sys::ABPOA_DEBUG_VERBOSE as i32 => Some(Verbosity::Debug),
            v if v == sys::ABPOA_LONG_DEBUG_VERBOSE as i32 => Some(Verbosity::LongDebug),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignMode {
    Global,
    Local,
    Extend,
}

impl AlignMode {
    pub fn as_raw(self) -> i32 {
        match self {
            AlignMode::Global => sys::ABPOA_GLOBAL_MODE as i32,
            AlignMode::Local => sys::ABPOA_LOCAL_MODE as i32,
            AlignMode::Extend => sys::ABPOA_EXTEND_MODE as i32,
        }
    }

    pub fn from_raw(value: i32) -> Option<Self> {
        match value {
            v if v == sys::ABPOA_GLOBAL_MODE as i32 => Some(AlignMode::Global),
            v if v == sys::ABPOA_LOCAL_MODE as i32 => Some(AlignMode::Local),
            v if v == sys::ABPOA_EXTEND_MODE as i32 => Some(AlignMode::Extend),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapMode {
    Linear,
    Affine,
    Convex,
}

impl GapMode {
    pub fn as_raw(self) -> i32 {
        match self {
            GapMode::Linear => sys::ABPOA_LINEAR_GAP as i32,
            GapMode::Affine => sys::ABPOA_AFFINE_GAP as i32,
            GapMode::Convex => sys::ABPOA_CONVEX_GAP as i32,
        }
    }

    pub fn from_raw(value: i32) -> Option<Self> {
        match value {
            v if v == sys::ABPOA_LINEAR_GAP as i32 => Some(GapMode::Linear),
            v if v == sys::ABPOA_AFFINE_GAP as i32 => Some(GapMode::Affine),
            v if v == sys::ABPOA_CONVEX_GAP as i32 => Some(GapMode::Convex),
            _ => None,
        }
    }
}

/// Gap penalty configuration to pair with substitution scores
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapPenalty {
    Linear {
        gap: i32,
    },
    Affine {
        gap_open: i32,
        gap_extend: i32,
    },
    Convex {
        gap_open: i32,
        gap_extend: i32,
        second_gap_open: i32,
        second_gap_extend: i32,
    },
}

impl GapPenalty {
    fn mode(self) -> GapMode {
        match self {
            GapPenalty::Linear { .. } => GapMode::Linear,
            GapPenalty::Affine { .. } => GapMode::Affine,
            GapPenalty::Convex { .. } => GapMode::Convex,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch: i32,
    pub gaps: GapPenalty,
}

impl Scoring {
    /// Linear gap model
    pub fn linear(match_score: i32, mismatch: i32, gap: i32) -> Self {
        Self {
            match_score,
            mismatch,
            gaps: GapPenalty::Linear { gap },
        }
    }

    /// Affine gap model
    pub fn affine(match_score: i32, mismatch: i32, gap_open: i32, gap_extend: i32) -> Self {
        Self {
            match_score,
            mismatch,
            gaps: GapPenalty::Affine {
                gap_open,
                gap_extend,
            },
        }
    }

    /// Convex gap model
    pub fn convex(
        match_score: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
        second_gap_open: i32,
        second_gap_extend: i32,
    ) -> Self {
        Self {
            match_score,
            mismatch,
            gaps: GapPenalty::Convex {
                gap_open,
                gap_extend,
                second_gap_open,
                second_gap_extend,
            },
        }
    }

    fn validate(self) -> Result<()> {
        if self.match_score < 0 {
            return Err(Error::InvalidInput("match score cannot be negative".into()));
        }
        if self.mismatch < 0 {
            return Err(Error::InvalidInput(
                "mismatch score cannot be negative".into(),
            ));
        }
        match self.gaps {
            GapPenalty::Linear { gap } => {
                if gap <= 0 {
                    return Err(Error::InvalidInput("linear gap must be positive".into()));
                }
            }
            GapPenalty::Affine {
                gap_open,
                gap_extend,
            } => {
                if gap_open <= 0 {
                    return Err(Error::InvalidInput(
                        "affine gap open must be positive".into(),
                    ));
                }
                if gap_extend <= 0 {
                    return Err(Error::InvalidInput(
                        "affine gap extension must be positive".into(),
                    ));
                }
            }
            GapPenalty::Convex {
                gap_open,
                gap_extend,
                second_gap_open,
                second_gap_extend,
            } => {
                if gap_open <= 0 || second_gap_open <= 0 {
                    return Err(Error::InvalidInput(
                        "convex gap opens must be positive".into(),
                    ));
                }
                if gap_extend <= 0 || second_gap_extend <= 0 {
                    return Err(Error::InvalidInput(
                        "convex gap extensions must be positive".into(),
                    ));
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsensusAlgorithm {
    HeaviestBundle,
    MostFrequent,
}

impl ConsensusAlgorithm {
    pub fn as_raw(self) -> i32 {
        match self {
            ConsensusAlgorithm::HeaviestBundle => sys::ABPOA_HB as i32,
            ConsensusAlgorithm::MostFrequent => sys::ABPOA_MF as i32,
        }
    }

    pub fn from_raw(value: i32) -> Option<Self> {
        match value {
            v if v == sys::ABPOA_HB as i32 => Some(ConsensusAlgorithm::HeaviestBundle),
            v if v == sys::ABPOA_MF as i32 => Some(ConsensusAlgorithm::MostFrequent),
            _ => None,
        }
    }
}

/// Newtype for graph node identifiers in abPOA
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(pub i32);

/// Sentinel node identifiers reserved by abPOA (source and sink)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SentinelNode {
    Source,
    Sink,
}

impl SentinelNode {
    pub const fn as_raw(self) -> i32 {
        match self {
            SentinelNode::Source => sys::ABPOA_SRC_NODE_ID as i32,
            SentinelNode::Sink => sys::ABPOA_SINK_NODE_ID as i32,
        }
    }

    pub const fn as_node_id(self) -> NodeId {
        NodeId(self.as_raw())
    }

    pub const fn from_raw(value: i32) -> Option<Self> {
        match value {
            v if v == sys::ABPOA_SRC_NODE_ID as i32 => Some(SentinelNode::Source),
            v if v == sys::ABPOA_SINK_NODE_ID as i32 => Some(SentinelNode::Sink),
            _ => None,
        }
    }
}

impl From<SentinelNode> for NodeId {
    fn from(node: SentinelNode) -> Self {
        node.as_node_id()
    }
}

impl From<SentinelNode> for i32 {
    fn from(node: SentinelNode) -> Self {
        node.as_raw()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;
    use std::ffi::CStr;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn parameters_init_and_free() {
        let mut params = Parameters::new().unwrap();
        assert!(!params.as_mut_ptr().is_null());
    }

    #[test]
    fn enum_mappings_roundtrip() {
        assert_eq!(
            AlignMode::from_raw(AlignMode::Global.as_raw()),
            Some(AlignMode::Global)
        );
        assert_eq!(
            GapMode::from_raw(GapMode::Affine.as_raw()),
            Some(GapMode::Affine)
        );
        assert_eq!(
            ConsensusAlgorithm::from_raw(ConsensusAlgorithm::MostFrequent.as_raw()),
            Some(ConsensusAlgorithm::MostFrequent)
        );
        assert_eq!(
            Verbosity::from_raw(Verbosity::Debug.as_raw()),
            Some(Verbosity::Debug)
        );
        assert_eq!(
            SentinelNode::from_raw(SentinelNode::Source.as_raw()),
            Some(SentinelNode::Source)
        );
        assert_eq!(
            SentinelNode::from_raw(SentinelNode::Sink.as_raw()),
            Some(SentinelNode::Sink)
        );
        assert_eq!(
            NodeId::from(SentinelNode::Source).0,
            SentinelNode::Source.as_raw()
        );
        assert!(SentinelNode::from_raw(123).is_none());
        assert!(AlignMode::from_raw(123).is_none());
    }

    #[test]
    fn alphabet_configuration_and_matrix_resize() {
        let mut params = Parameters::new().unwrap();
        assert_eq!(params.get_alphabet(), Some(Alphabet::Dna));

        params.set_alphabet(Alphabet::AminoAcid).unwrap();
        assert_eq!(params.get_alphabet(), Some(Alphabet::AminoAcid));

        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.m, 27);
        assert!(!raw.mat.is_null(), "score matrix is allocated");
    }

    #[test]
    fn quality_and_heuristic_flags_toggle() {
        let mut params = Parameters::new().unwrap();
        params.set_use_quality(true);
        assert_eq!(unsafe { params.raw.as_ref() }.use_qv(), 1);

        params.set_use_quality(false);
        assert_eq!(unsafe { params.raw.as_ref() }.use_qv(), 0);

        #[allow(deprecated)]
        params
            .set_progressive_poa(true)
            .set_inc_path_score(true)
            .set_sub_alignment(true)
            .set_sort_input_seq(true)
            .set_gap_on_right(true)
            .set_gap_at_end(true);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.progressive_poa(), 1);
        assert_eq!(raw.inc_path_score, 1);
        assert_eq!(raw.sub_aln(), 1);
        assert_eq!(raw.sort_input_seq, 1);
        assert_eq!(raw.put_gap_on_right(), 1);
        assert_eq!(raw.put_gap_at_end(), 1);

        params.set_gap_on_right(false).set_gap_at_end(false);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.put_gap_on_right(), 0);
        assert_eq!(raw.put_gap_at_end(), 0);
    }

    #[test]
    fn seeding_controls_and_disable_flag() {
        let mut params = Parameters::new().unwrap();
        params.set_minimizer_seeding(11, 5, 100).unwrap();
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.k, 11);
        assert_eq!(raw.w, 5);
        assert_eq!(raw.min_w, 100);
        assert_eq!(raw.disable_seeding(), 0);

        assert!(params.set_minimizer_seeding(0, 5, 1).is_err());
        params.set_disable_seeding(true);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.disable_seeding(), 1);
    }

    #[test]
    fn zdrop_and_end_bonus_are_configurable() {
        let mut params = Parameters::new().unwrap();
        assert!(!params.dirty);

        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.zdrop, -1);
        assert_eq!(raw.end_bonus, -1);

        params.set_zdrop(Some(10)).unwrap();
        params.set_end_bonus(Some(5)).unwrap();
        assert!(params.dirty);
        let _ = params.as_mut_ptr();
        assert!(!params.dirty);

        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.zdrop, 10);
        assert_eq!(raw.end_bonus, 5);

        params.set_zdrop(None).unwrap();
        params.set_end_bonus(None).unwrap();
        let _ = params.as_mut_ptr();

        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.zdrop, -1);
        assert_eq!(raw.end_bonus, -1);

        assert!(params.set_zdrop(Some(0)).is_err());
        assert!(params.set_end_bonus(Some(0)).is_err());
    }

    #[test]
    fn scoring_models_and_substitution_scores() {
        let mut params = Parameters::new().unwrap();
        params.set_scoring_scheme(Scoring::linear(3, 5, 2)).unwrap();
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.match_, 3);
        assert_eq!(raw.mismatch, 5);
        assert_eq!(raw.gap_mode, GapMode::Linear.as_raw());
        assert_eq!(raw.gap_open1, 0);
        assert_eq!(raw.gap_ext1, 2);
        assert_eq!(raw.gap_open2, 0);
        assert_eq!(raw.gap_ext2, 0);

        params
            .set_scoring_scheme(Scoring::affine(4, 6, 7, 8))
            .unwrap();
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.match_, 4);
        assert_eq!(raw.mismatch, 6);
        assert_eq!(raw.gap_mode, GapMode::Affine.as_raw());
        assert_eq!(raw.gap_open1, 7);
        assert_eq!(raw.gap_ext1, 8);
        assert_eq!(raw.gap_open2, 0);
        assert_eq!(raw.gap_ext2, 0);

        params
            .set_scoring_scheme(Scoring::convex(5, 7, 9, 10, 11, 12))
            .unwrap();
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.match_, 5);
        assert_eq!(raw.mismatch, 7);
        assert_eq!(raw.gap_mode, GapMode::Convex.as_raw());
        assert_eq!(raw.gap_open1, 9);
        assert_eq!(raw.gap_ext1, 10);
        assert_eq!(raw.gap_open2, 11);
        assert_eq!(raw.gap_ext2, 12);

        params
            .set_scoring_scheme(Scoring::linear(11, 13, 2))
            .unwrap();
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.match_, 11);
        assert_eq!(raw.mismatch, 13);
    }

    #[test]
    fn scoring_scheme_validates_inputs() {
        let mut params = Parameters::new().unwrap();
        assert!(
            params
                .set_scoring_scheme(Scoring::linear(-1, 2, 3))
                .is_err()
        );
        assert!(params.set_scoring_scheme(Scoring::linear(1, 1, 0)).is_err());
        assert!(
            params
                .set_scoring_scheme(Scoring::affine(1, 1, 0, 2))
                .is_err()
        );
        assert!(
            params
                .set_scoring_scheme(Scoring::convex(1, 1, 2, 0, 3, 1))
                .is_err()
        );
    }

    #[test]
    fn outputs_and_finalize_flags() {
        let mut params = Parameters::new().unwrap();
        assert_eq!(params.outputs(), OutputMode::CONSENSUS | OutputMode::MSA);

        let previous = params.set_outputs(OutputMode::MSA);
        assert_eq!(previous, OutputMode::CONSENSUS | OutputMode::MSA);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.out_cons(), 0);
        assert_eq!(raw.out_msa(), 1);
        assert_eq!(params.outputs(), OutputMode::MSA);

        params
            .set_align_mode(AlignMode::Local)
            .set_scoring_scheme(Scoring::affine(5, 7, 4, 2))
            .unwrap()
            .set_band(6, 0.5);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.align_mode, AlignMode::Local.as_raw());
        assert_eq!(raw.wb, 6);
        assert!((raw.wf - 0.5).abs() < f32::EPSILON);

        let _ = params.as_mut_ptr();
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.gap_mode, GapMode::Affine.as_raw());
        assert_eq!(raw.use_read_ids(), 1);
        assert_eq!(raw.wb, -1, "local mode disables banding");

        params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.out_cons(), 1);
        assert_eq!(raw.out_msa(), 1);
        assert_eq!(params.outputs(), OutputMode::CONSENSUS | OutputMode::MSA);
    }

    #[test]
    fn use_read_ids_forced_for_msa_output() {
        let mut params = Parameters::configure().unwrap();
        params.set_use_read_ids(false);
        let _ = params.as_mut_ptr();

        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.out_msa(), 1);
        assert_eq!(raw.use_read_ids(), 1);
    }

    #[test]
    fn set_outputs_for_call_updates_use_read_ids() {
        let mut params = Parameters::new().unwrap();
        params.set_outputs(OutputMode::CONSENSUS);
        params.set_use_read_ids(false);
        let _ = params.as_mut_ptr();

        assert!(!params.dirty);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.out_msa(), 0);
        assert_eq!(raw.use_read_ids(), 0);

        params.set_outputs_for_call(OutputMode::MSA);

        assert!(!params.dirty);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.out_msa(), 1);
        assert_eq!(raw.use_read_ids(), 1);
    }

    #[test]
    fn consensus_and_verbosity_limits() {
        let mut params = Parameters::new().unwrap();
        params.set_max_consensus(2).unwrap();
        params.set_min_cluster_freq(0.4).unwrap();
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.max_n_cons, 2);
        assert!((raw.min_freq - 0.4).abs() < f64::EPSILON);
        assert!(params.set_max_consensus(0).is_err());
        assert!(params.set_min_cluster_freq(1.5).is_err());

        params
            .set_consensus(ConsensusAlgorithm::HeaviestBundle, 1, 0.2)
            .unwrap();
        params.set_verbosity(Verbosity::Debug);
        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.cons_algrm, ConsensusAlgorithm::HeaviestBundle.as_raw());
        assert_eq!(raw.max_n_cons, 1);
        assert!((raw.min_freq - 0.2).abs() < f64::EPSILON);
        assert_eq!(raw.verbose, Verbosity::Debug.as_raw());
    }

    #[test]
    fn score_matrix_and_custom_matrix_are_applied() {
        let mut params = Parameters::new().unwrap();
        let mut path = env::temp_dir();
        let suffix = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        path.push(format!("abpoa-matrix-{suffix}.txt"));
        let file_contents = "\
# simple 5x5 matrix
 A C G T N
A 4 -3 -3 -3 0
C -3 8 -3 -3 0
G -3 -3 4 -3 0
T -3 -3 -3 4 0
N 0 0 0 0 0
";
        fs::write(&path, file_contents).unwrap();

        params.set_score_matrix_file(&path).unwrap();
        let _ = params.as_mut_ptr();

        let raw = unsafe { params.raw.as_ref() };
        assert_eq!(raw.use_score_matrix, 1);
        let stored_path = unsafe { CStr::from_ptr(raw.mat_fn) }
            .to_string_lossy()
            .into_owned();
        assert!(
            stored_path.contains("abpoa-matrix-"),
            "stored path should include configured file name"
        );
        let matrix = unsafe { std::slice::from_raw_parts(raw.mat, 25) };
        assert_eq!(matrix[0], 4);
        assert_eq!(matrix[1], -3);
        assert_eq!(raw.max_mat, 8);
        assert_eq!(raw.min_mis, 3);

        let custom = vec![
            6, -1, -1, -1, 0, //
            -1, 6, -1, -1, 0, //
            -1, -1, 8, -1, 0, //
            -1, -1, -1, 6, 0, //
            0, 0, 0, 0, 0,
        ];
        params.set_custom_matrix(5, &custom).unwrap();
        let _ = params.as_mut_ptr();
        let raw = unsafe { params.raw.as_ref() };
        let restored = unsafe { std::slice::from_raw_parts(raw.mat, custom.len()) };
        assert_eq!(restored, custom.as_slice());
        assert_eq!(raw.use_score_matrix, 0);
        assert_eq!(raw.max_mat, 8);
        assert_eq!(raw.min_mis, 1);

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn incremental_graph_path_lifecycle() {
        let mut params = Parameters::new().unwrap();
        params.set_incremental_graph_file("graph.gfa").unwrap();
        let raw = unsafe { params.raw.as_ref() };
        assert!(!raw.incr_fn.is_null());

        params.clear_incremental_graph();
        let raw = unsafe { params.raw.as_ref() };
        assert!(raw.incr_fn.is_null());
    }
}
