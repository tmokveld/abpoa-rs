#![allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
#![allow(clippy::all)]

//! Raw FFI bindings to the vendored abPOA C library.
//!
//! Most users should prefer the safe wrapper crate `abpoa`

pub mod sys {
    #![allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
    #![allow(clippy::all)]

    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

unsafe extern "C" {
    /// Reallocate the internal sequence buffers for `abpoa_seq_t`
    pub fn abpoa_realloc_seq(abs: *mut sys::abpoa_seq_t) -> *mut sys::abpoa_seq_t;

    /// Copy a byte slice into an `abpoa_str_t` (null-terminating it)
    pub fn abpoa_cpy_str(str_: *mut sys::abpoa_str_t, s: *mut libc::c_char, l: libc::c_int);

    /// Initialize the global 16-bit integer log2 lookup table (`ab_LogTable65536`)
    pub fn set_65536_table();

    /// Initialize the global 16-bit bit-count lookup table (`ab_bit_table16`)
    pub fn set_bit_table16();
}

/// Pointer wrapper for accessing abPOA's internal lookup tables
pub struct ExternPtr<T>(*const T);

unsafe impl<T> Sync for ExternPtr<T> {}

impl<T> ExternPtr<T> {
    /// Raw pointer to the underlying C storage
    pub fn as_ptr(&self) -> *const T {
        self.0
    }
}

/// Lookup table mapping DNA base bytes to the 0..=4 alphabet used by abPOA
///
/// This matches the upstream `ab_nt4_table`: `A/a => 0`, `C/c => 1`, `G/g => 2`, `T/t/U/u => 3`,
/// and everything else maps to `4` (`N`)
pub static AB_NT4_TABLE: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_nt4_table"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Lookup table mapping bytes to canonical DNA bases
///
/// This matches the upstream `ab_nt256_table`: indices `0..=4` map to `A/C/G/T/N`, `5` maps to
/// `'-'` (gap), `27` maps to `'-'` (alternate gap sentinel), and ASCII bytes like `b'A'`/`b'a'`
/// are normalized to their uppercase base
pub static AB_NT256_TABLE: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_nt256_table"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Lookup table mapping amino-acid bytes to the 0..=26 alphabet used by abPOA
///
/// Unknown residues map to `26` (same as `'*'`)
pub static AB_AA26_TABLE: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_aa26_table"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Lookup table mapping bytes to canonical amino-acid residue letters
///
/// This matches the upstream `ab_aa256_table`: indices `0..=26` map to residue letters (with
/// `26` => `'*'`) and `27` maps to `'-'` (gap)
pub static AB_AA256_TABLE: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_aa256_table"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Pre-computed integer log2 for `[0, 2^16)`
///
/// The underlying table is populated by calling [`set_65536_table`]
pub static AB_LOG_TABLE_65536: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_LogTable65536"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Pre-computed popcount (bit-count) table for 16-bit values
///
/// The underlying table is populated by calling [`set_bit_table16`].
pub static AB_BIT_TABLE_16: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_bit_table16"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

pub use sys::*;
