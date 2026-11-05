//! # HARICOT
//!
//! **Genome assembly directly from raw nanopore signal**: the wtdbg2
//! fuzzy-Bruijn-graph approach carried over from nucleotide k-mers to raw-signal
//! "e-mers". It avoids basecalling, and the all-vs-all overlap cost of OLC-style raw
//! assembly (cf. Rawsamble).
//!
//! ## Module map
//!
//! | stage | module |
//! |-------|--------|
//! | raw signal I/O (pod5; fast5 behind a feature) | [`signal::io`] |
//! | normalization | [`signal::normalize`] |
//! | event segmentation (dual t-test) | [`signal::events`] |
//! | quantization | [`emer::quantize`] |
//! | e-mer sketching / seed selection | [`emer::sketch`], [`mod@emer::select`] |
//! | read overlap (KBM analog): seed index, Hough overlapper or banded bin-DP | [`overlap`] |
//! | fuzzy Bruijn graph: data model + construction | [`graph::fbg`], [`graph::build`] |
//! | graph cleaning (edge-cov filter, transitive reduction, tips, bubbles) | [`graph::clean`] |
//! | layout → unitigs/contigs (wtdbg2's fragment layout), GFA | [`layout::unitig`], [`layout::frg`], [`layout::gfa`] |
//! | read order: the tie rank by signal content, the canonical pair frame | [`readorder`] |
//! | the full assembler, `haricot assemble` | [`mod@assemble`] |
//!
//! Signal-processing defaults follow the RawHash2 **R10.4.1** preset (see
//! [`params`]). Ported code cites the C reference implementations by `file:line`, as
//! `extern/RawHash2/...` and `extern/wtdbg2/...` (clone them there to follow a citation;
//! they are not distributed), or by file name alone (`kbm.h`, `wtdbg.c`, ...).

/// Process-wide allocator.
///
/// Declared in the library rather than in `main.rs` deliberately: the pipeline programs
/// (`pipeline/src/bin/`) are separate binaries, and a `#[global_allocator]` in a linked
/// library applies to the whole program, so this one declaration covers every program
/// instead of each having to remember to opt in.
///
/// It is a pure performance choice with no effect on results. The overlap stage is
/// 91–99% of wall time and its dominant phase (`dp_align`, 70.6% on E. coli fwd-30x)
/// allocates a short-lived scratch `Vec<f64>` for every slope evaluation of every read
/// pair — order 1.2e9 allocations, from all 32 threads at once. glibc's per-arena
/// locking turns that into the scaling limit of a phase that is otherwise
/// embarrassingly parallel.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

pub mod assemble;
pub mod cli;
pub mod emer;
pub mod graph;
pub mod layout;
pub mod overlap;
pub mod params;
pub mod pipeline;
pub mod readorder;
pub mod signal;

pub use params::Params;
