use std::sync::{Mutex, Once};

use crate::sys;

static ABPOA_GLOBAL_WRITE_LOCK: Mutex<()> = Mutex::new(());
static OUTPUT_TABLES_INIT: Once = Once::new();

pub(crate) fn with_abpoa_global_write_lock<T>(f: impl FnOnce() -> T) -> T {
    let guard = ABPOA_GLOBAL_WRITE_LOCK
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    let result = f();
    drop(guard);
    result
}

pub(crate) fn ensure_output_tables() {
    OUTPUT_TABLES_INIT.call_once(|| {
        with_abpoa_global_write_lock(|| unsafe {
            // Safety: these C functions write to global lookup tables in the abPOA library.
            // They are deterministic and intended to be run once during process initialization.
            sys::set_65536_table();
            sys::set_bit_table16();
        });
    });
}

