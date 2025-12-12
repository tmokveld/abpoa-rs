#![allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
#![allow(clippy::all)]

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

/// Transforms AaCcGgTtNn to 0-4
pub static AB_NT4_TABLE: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_nt4_table"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Transforms 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
pub static AB_NT256_TABLE: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_nt256_table"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Transforms amino acids to 0-25
pub static AB_AA26_TABLE: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_aa26_table"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Transforms 0-25 to amino acids.
pub static AB_AA256_TABLE: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_aa256_table"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Pre-computed ilog2 for [0, 2^16)
pub static AB_LOG_TABLE_65536: ExternPtr<u8> = {
    unsafe extern "C" {
        #[link_name = "ab_LogTable65536"]
        static inner: u8;
    }
    // Safety: `inner` is a statically allocated table in the abPOA C library
    // and is valid throughout the process
    ExternPtr(unsafe { &inner as *const u8 })
};

/// Pre-computed bit flag table for 16-bit values
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
