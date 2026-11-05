# vendor/ — patched third-party crates

## arrow-ipc 53.4.1 (Apache-2.0, https://github.com/apache/arrow-rs)

`arrow-ipc/` is the crates.io release 53.4.1 of the Apache Arrow IPC crate, wired in through
`[patch.crates-io]` in `../Cargo.toml` (the version `pod5lib` pulls in). Two changes:

- `src/reader.rs::read_dictionary_impl`: reads delta dictionary batches (below);
- `src/lib.rs`: `#![allow(mismatched_lifetime_syntaxes)]`, a lint added to Rust after this release.

**Why.** The Rust pod5 reader (`pod5lib`, on arrow 53) failed with
`Arrow error: Invalid argument error: delta dictionary batches not supported` on pod5s written by the
Python pod5 writer (`pod5` 0.3.39) when a dictionary column (run info, end reason, pore type) gains a
value after the first record batch was flushed, e.g. reads of several runs written in an order where a run
first appears after the first batch (a permuted read order, `checks/pod5_permute.py`).

**The patch.** Upstream 53.4.1 returns that error for any dictionary batch with `isDelta`. The patched
reader reads the delta batch like any other and appends its values to the dictionary already held for that
id (the Arrow columnar format: a delta dictionary batch extends the existing dictionary; record batches
index into the extended one). Non-delta dictionary batches take the unchanged path, so files without deltas
read exactly as before. A delta for an id with no dictionary yet, or of another value type, is an error.
Later arrow-rs releases support deltas natively; moving to one means moving `pod5lib`'s arrow major version.
