use abpoa_sys::sys;

#[test]
fn init_and_free_handles() {
    unsafe {
        let params = sys::abpoa_init_para();
        assert!(!params.is_null(), "abpoa_init_para returned null");
        sys::abpoa_free_para(params);

        let handle = sys::abpoa_init();
        assert!(!handle.is_null(), "abpoa_init returned null");
        sys::abpoa_free(handle);
    }
}
