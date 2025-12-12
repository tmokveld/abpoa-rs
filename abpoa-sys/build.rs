use std::{
    env,
    path::{Path, PathBuf},
};

use glob::glob;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let target_triple = env::var("TARGET").unwrap_or_default();
    let is_x86 = matches!(target_arch.as_str(), "x86_64" | "x86" | "i686");
    let is_linux = target_os == "linux";

    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let abpoa_dir = manifest_dir.join("abPOA");

    let include_dir = abpoa_dir.join("include");
    let src_dir = abpoa_dir.join("src");

    let mut include_paths = vec![include_dir.clone()];

    // Rerun when any vendored source/header files change
    let glob_patterns = [
        format!("{}/src/*.c", abpoa_dir.display()),
        format!("{}/src/*.h", abpoa_dir.display()),
        format!("{}/include/*.h", abpoa_dir.display()),
        format!("{}/include/simde/**/*.h", abpoa_dir.display()),
    ];
    for pattern in glob_patterns {
        for path in glob(&pattern)
            .expect("Failed to read glob pattern")
            .flatten()
        {
            println!("cargo:rerun-if-changed={}", path.display());
        }
    }

    println!(
        "cargo:rerun-if-changed={}",
        include_dir.join("abpoa.h").display()
    );
    println!("cargo:rerun-if-changed={}", src_dir.display());

    let mut base_build = cc::Build::new();
    base_build.warnings(false);
    base_build.include(&include_dir);
    base_build.include(&src_dir);
    if is_linux {
        base_build.define("_GNU_SOURCE", None);
    }
    base_build.flag_if_supported("-std=c99");
    base_build.flag_if_supported("-Wno-unused-function");
    base_build.flag_if_supported("-Wno-misleading-indentation");
    base_build.flag_if_supported("-Wno-unused-parameter");
    base_build.flag_if_supported("-Wno-sign-compare");

    // On non-x86 targets we must use SIMDe
    let mut clang_args = vec![
        "-Wno-unused-function".into(),
        "-Wno-misleading-indentation".into(),
        "-Wno-unused-parameter".into(),
        "-Wno-sign-compare".into(),
    ];
    if is_linux {
        clang_args.push("-D_GNU_SOURCE".into());
    }
    let use_simde = !is_x86;
    if use_simde {
        base_build.define("USE_SIMDE", None);
        base_build.define("SIMDE_ENABLE_NATIVE_ALIASES", None);
        clang_args.push("-DUSE_SIMDE".into());
        clang_args.push("-DSIMDE_ENABLE_NATIVE_ALIASES".into());
        // abPOA expects __AVX2__ even on ARM to pick the wider SIMD kernels when using SIMDe
        if matches!(target_arch.as_str(), "aarch64" | "arm" | "arm64") {
            clang_args.push("-D__AVX2__".into());
            if target_os == "macos" {
                clang_args.push("-mcpu=apple-m1".into());
            } else {
                clang_args.push("-march=armv8-a+simd".into());
            }
        }
    } else {
        // Let x86/x86_64 use native intrinsics via dispatch variants below
    }

    let base_c_files = [
        "abpoa_align.c",
        "abpoa_graph.c",
        "abpoa_output.c",
        "abpoa_plot.c",
        "abpoa_seed.c",
        "abpoa_seq.c",
        "abpoa_simd.c",
        "kalloc.c",
        "kstring.c",
        "utils.c",
    ];

    for file in base_c_files {
        base_build.file(src_dir.join(file));
    }

    let use_dispatch = is_x86;

    base_build.compile("abpoa");

    // If not on x86, build the single SIMD object (SIMDe path) with arch-specific flags
    if !use_dispatch {
        let mut align_build = cc::Build::new();
        align_build.warnings(false);
        align_build.include(&include_dir);
        align_build.include(&src_dir);
        if is_linux {
            align_build.define("_GNU_SOURCE", None);
        }
        align_build.flag_if_supported("-std=c99");
        align_build.flag_if_supported("-Wno-unused-function");
        align_build.flag_if_supported("-Wno-misleading-indentation");
        align_build.flag_if_supported("-Wno-unused-parameter");
        align_build.flag_if_supported("-Wno-sign-compare");
        if use_simde {
            align_build.define("USE_SIMDE", None);
            align_build.define("SIMDE_ENABLE_NATIVE_ALIASES", None);
            // abPOA expects __AVX2__ even on ARM to pick the wider SIMD kernels when using SIMDe
            if matches!(target_arch.as_str(), "aarch64" | "arm" | "arm64") {
                align_build.define("__AVX2__", None);
                if target_os == "macos" {
                    align_build.flag_if_supported("-mcpu=apple-m1");
                } else {
                    align_build.flag_if_supported("-march=armv8-a+simd");
                }
            }
        }
        align_build.file(src_dir.join("abpoa_align_simd.c"));
        align_build.compile("abpoa_align_simd");
    } else {
        // For x86/x86_64, also build the runtime-dispatch variants and dispatcher to match upstream
        build_dispatch_variants(&include_dir, &src_dir, is_linux);
    }

    // abPOA depends on zlib and libm; pthreads is part of libSystem on macOS
    if target_os != "windows" {
        println!("cargo:rustc-link-lib=z");
        println!("cargo:rustc-link-lib=m");
        println!("cargo:rustc-link-lib=pthread");
    } else {
        println!("cargo:rustc-link-lib=z");
    }

    include_paths.push(include_dir);
    let header = find_header(&include_paths).unwrap_or_else(|| PathBuf::from("abpoa.h"));

    generate_bindings(&header, &include_paths, &clang_args, &target_triple);
}

fn find_header(include_paths: &[PathBuf]) -> Option<PathBuf> {
    include_paths
        .iter()
        .map(|p| p.join("abpoa.h"))
        .find(|p| p.exists())
}

fn build_dispatch_variants(include_dir: &PathBuf, src_dir: &PathBuf, is_linux: bool) {
    // Dispatcher (CPUID-based)
    let mut dispatch = cc::Build::new();
    dispatch
        .warnings(false)
        .include(include_dir)
        .include(src_dir)
        .define("ABPOA_SIMD_DISPATCH", None)
        .flag_if_supported("-std=c99")
        .flag_if_supported("-Wno-unused-function")
        .flag_if_supported("-Wno-misleading-indentation")
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-sign-compare")
        .file(src_dir.join("abpoa_dispatch_simd.c"));
    if is_linux {
        dispatch.define("_GNU_SOURCE", None);
    }
    dispatch.compile("abpoa_dispatch");

    // Helper to compile a single SIMD variant of abpoa_align_simd.c with the given flags.
    let compile_variant = |name: &str, extra_flags: &[&str]| {
        let mut b = cc::Build::new();
        b.warnings(false)
            .include(include_dir)
            .include(src_dir)
            .define("ABPOA_SIMD_DISPATCH", None)
            .flag_if_supported("-std=c99")
            .flag_if_supported("-Wno-unused-function")
            .flag_if_supported("-Wno-misleading-indentation")
            .flag_if_supported("-Wno-unused-parameter")
            .flag_if_supported("-Wno-sign-compare");
        if is_linux {
            b.define("_GNU_SOURCE", None);
        }
        for f in extra_flags {
            b.flag_if_supported(f);
        }
        b.file(src_dir.join("abpoa_align_simd.c"));
        b.compile(name);
    };

    // SSE2 variant (clear __SSE4_1__ to force SSE2, mirrors upstream Makefile)
    compile_variant("abpoa_align_simd_sse2", &["-msse2", "-U__SSE4_1__"]);
    // SSE4.1 variant
    compile_variant("abpoa_align_simd_sse41", &["-msse4.1"]);
    // AVX2 variant
    compile_variant("abpoa_align_simd_avx2", &["-mavx2"]);
    // AVX512BW variant
    compile_variant("abpoa_align_simd_avx512bw", &["-mavx512bw"]);
}

fn generate_bindings(
    header: &Path,
    include_paths: &[PathBuf],
    clang_args: &[String],
    target_triple: &str,
) {
    let mut builder = bindgen::Builder::default()
        .header(header.to_string_lossy())
        .allowlist_function("abpoa_.*")
        .allowlist_type("abpoa_.*")
        .allowlist_var("ABPOA_.*")
        .clang_arg("-DABPOA_RUST_BINDGEN=1")
        .clang_arg(format!("--target={target_triple}"));

    for dir in include_paths {
        builder = builder.clang_arg(format!("-I{}", dir.display()));
    }
    for arg in clang_args {
        builder = builder.clang_arg(arg);
    }

    let bindings = builder
        .generate()
        .expect("Unable to generate abPOA bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap()).join("bindings.rs");
    bindings
        .write_to_file(out_path)
        .expect("Could not write bindings!");
}
