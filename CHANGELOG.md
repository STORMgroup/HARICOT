# Changelog

## 1.0.1 (2026-11-13)

Code and documentation cleanup. Assemblies are unchanged: every output file was byte-identical to the manuscript's
runs on the 54 of its 56 test sets that were rerun, with the published forests.

- **One crate.** The pipeline programs are binaries of the `haricot` crate (`src/bin/haricot-*.rs`); the code they
  share is `src/pipeline/`. The binary names are unchanged.
- **CLI.** `haricot` has one subcommand, `assemble`. The subcommands `baseline`, `events`, `sketch`, `overlap` and
  `info` are removed, with the library code only they reached (the bin-DP overlapper, the OLC string graph, syncmer
  and strobemer seeding, the other quantizers). None is part of the manuscript's method.
- **Keep rule.** The stage-2 keep rule is named as in the manuscript: `haricot-layout --keep-level depth-adaptive`
  (module `pipeline::keep_cut`). The earlier names (`policy7`, `--s2-drop`, `haricot-forest --s2-dir`,
  `pipeline::lo`) are still accepted.
- **Reproduction.** `reproduce/` holds one script per step, from the public raw reads to the results table; its
  README replaces `docs/REPRODUCE.md`. `reproduce/expected.sha256` holds the SHA-256 of each test run's
  `unitigs.gfa`.
- **Documentation.** READMEs for `tests/` and `evaluation/`; each program's documentation names its step of
  `haricot assemble` and its algorithm; help texts and script headers describe what the code does.
- **Packaging.** `rust-version` is 1.96, the pinned toolchain. `vendor/arrow-ipc` carries Apache Arrow's
  `LICENSE.txt` and `NOTICE.txt`. Dependencies no code uses are removed (`hdf5-metno` and the `fast5` feature,
  `thiserror`, `log`, `env_logger`).

## 1.0.0 (2026-11-05)

The release of the manuscript's runs (tag `v1.0.0`; its `Cargo.toml` still said 0.1.0).
