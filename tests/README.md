# Tests

| tests | what they check | where |
|---|---|---|
| unit tests | each module of the library and of the programs, on small constructed inputs | `#[cfg(test)] mod tests` in `src/` and `src/bin/` |
| `e2e.rs` | two real reads: their events and seeds, the overlap that `haricot assemble`'s detection finds between them, and the two POD5 readers | `tests/` |
| `hash_compat.rs` | the e-mer hashes, bit for bit against RawHash2's | `tests/`; goldens in `rawhash_compat/` |
| `eval/` | the evaluation scripts and `checks/cmp_run.sh`, on small synthetic inputs | `tests/eval/` |

```bash
cargo test --release --locked                  # unit tests, e2e.rs, hash_compat.rs
python -m unittest discover -s tests/eval      # in the haricot environment (envs/haricot.yml)
```

`e2e.rs` runs on `tests/fixtures/debug_reads.pod5`, two overlapping reads of *E. coli* that are not distributed
([`fixtures/README.md`](fixtures/README.md) gives their source). Without the file, `cargo test` reports its three tests
as ignored.

RawHash2's own sketching code wrote the goldens once ([`rawhash_compat/README.md`](rawhash_compat/README.md)), so
`hash_compat.rs` needs neither RawHash2 nor its dependencies.

The assemblies themselves are checked by the scripts in `checks/`: invariance to the order and the names of the reads,
and byte identity of two runs ([Checks](../reproduce/README.md#checks)).
