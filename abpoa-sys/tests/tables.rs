use abpoa_ffi_sys::{AB_AA26_TABLE, AB_AA256_TABLE};

#[test]
fn aa_tables_have_expected_values() {
    unsafe {
        let aa26_ptr = AB_AA26_TABLE.as_ptr();
        let aa256_ptr = AB_AA256_TABLE.as_ptr();

        // ab_aa26_table maps ASCII amino-acid chars -> 0..=25 (or 26 for unknown).
        assert_eq!(*aa26_ptr.add(b'A' as usize), 0);
        assert_eq!(*aa26_ptr.add(b'C' as usize), 1);
        assert_eq!(*aa26_ptr.add(b'G' as usize), 2);
        assert_eq!(*aa26_ptr.add(b'T' as usize), 3);
        assert_eq!(*aa26_ptr.add(b'N' as usize), 4);

        // ab_aa256_table maps 0..=25 -> ASCII amino-acid chars.
        assert_eq!(*aa256_ptr.add(0), b'A');
        assert_eq!(*aa256_ptr.add(1), b'C');
        assert_eq!(*aa256_ptr.add(2), b'G');
        assert_eq!(*aa256_ptr.add(3), b'T');
        assert_eq!(*aa256_ptr.add(4), b'N');

        // The two tables should not alias each other.
        assert_ne!(aa26_ptr, aa256_ptr);
    }
}
