# RawHash2 ↔ Rust hash-compatibility fixtures

These files pin the Rust e-mer hashing (`src/emer/{quantize,sketch}.rs`) to
RawHash2's reference implementation. The Rust side is exercised by
`tests/hash_compat.rs`, which runs on every `cargo test`; the C reference is run
**once** to bake the golden hashes, so the test itself needs no system deps
(libhdf5/pod5/etc.) or even the RawHash2 sources.

## Files

| file | committed? | what |
|------|-----------|------|
| `dump_hashes.c` | yes | tiny C harness: generates the signal and dumps RawHash2 hashes |
| `regenerate.sh` | yes | builds the harness against `extern/RawHash2/src` and regenerates everything below |
| `signal.txt` | yes | 2000 deterministic pseudo-random "normalized event" values (one f32 per line) |
| `*.tsv` | yes | golden `pos<TAB>hash` dumps, one per parameter configuration |

`signal.txt` is written with `%.9g`, which uniquely round-trips an `f32`, so the
C generator, the C dumper, and Rust's `f32::from_str` all parse identical
values — the comparison can't drift on float formatting.

## What is compared

All configs use the **dense** sketch (`w = 0`, RawHash2's `ri_sketch_reg`), so
seeds come out in event order and the test compares the `(pos, hash)` sequence
position-by-position. `hash` is the bare 64-bit value with RawHash2's packing
stripped (`x >> RI_HASH_SHIFT`); `pos` is the first-event index of the e-mer.
This exercises `dynamic_quantize`, the sig-diff collapse (`diff`), e-mer
packing, and `hash64` end to end.

Minimizer selection (`w > 0`) is intentionally **not** compared: `sketch.rs`
documents a deliberate deviation from RawHash2's duplicate-handling rotation.

The configurations (see `CONFIGS` in `tests/hash_compat.rs` and the `emit` lines
in `regenerate.sh`) cover the R10.4.1 default (`e=8, diff=0`), the sig-diff
filter disabled (`diff=-1`), and a couple of other `e` values.

## Regenerating

Only needed if you change `signal.txt`, the configs, or RawHash2's sketch code:

```bash
tests/rawhash_compat/regenerate.sh   # needs a C compiler + RawHash2 sources in extern/RawHash2
```

It is deterministic — same seed, same signal, same goldens.
