//! Sets `cfg(haricot_debug_fixture)` when `tests/fixtures/debug_reads.pod5` exists: the end-to-end tests
//! (`tests/e2e.rs`) run on it, and are reported as ignored without it (the file is not distributed).

fn main() {
    println!("cargo::rustc-check-cfg=cfg(haricot_debug_fixture)");
    println!("cargo::rerun-if-changed=tests/fixtures");
    if std::path::Path::new("tests/fixtures/debug_reads.pod5").is_file() {
        println!("cargo::rustc-cfg=haricot_debug_fixture");
    }
}
