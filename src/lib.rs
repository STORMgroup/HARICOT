//! # HARICOT
//!
//! **Genome assembly directly from raw nanopore signals**, without basecalling: wtdbg2's fuzzy de
//! Bruijn graph built over raw-signal *e-mers* (runs of quantized signal events) instead of
//! nucleotide k-mers. The output is a layout (unitigs and contigs as paths of read segments), with
//! no consensus sequence.
//!
//! `haricot assemble` ([`mod@assemble`]) runs the pipeline as a chain of programs (`src/bin/`):
//! candidate overlap detection by a Hough transform over shared e-mer seeds (`haricot-detect`), the
//! pair features and the stage-1 forest (`haricot-pair-features`, `haricot-forest`), the graph
//! features and the stage-2 forest (`haricot-graph-features`, `haricot-forest`), the
//! depth-adaptive keep cut, the graph build, cleaning, layout and end bridges (`haricot-layout
//! layout`), the contig-end rescue and its forest (`haricot-rescue`, `haricot-pair-features`,
//! `haricot-forest`), and the final build (`haricot-layout final`). The
//! manuscript's Supplementary Algorithms S1-S9 give the pseudocode; module docs cite them.
//!
//! ## Module map
//!
//! | stage | module |
//! |-------|--------|
//! | raw signal I/O (pod5) | [`signal::io`] |
//! | normalization | [`mod@signal::normalize`] |
//! | event segmentation (dual t-test) | [`signal::events`] |
//! | quantization | [`emer::quantize`] |
//! | e-mer sketching (minimizers) | [`mod@emer::sketch`] |
//! | read overlap (KBM analog): seed index, pre-filter, Hough overlapper | [`overlap`] |
//! | fuzzy de Bruijn graph: data model + construction | [`graph::fbg`], [`graph::build`] |
//! | graph cleaning (edge filters, transitive reduction, tips, bubbles) | [`mod@graph::clean`] |
//! | layout → unitigs/contigs (wtdbg2's fragment layout), GFA | [`layout::unitig`], [`layout::frg`], [`layout::gfa`] |
//! | read order: the tie rank by signal content, the canonical pair frame | [`readorder`] |
//! | parameters (the R10.4.1 preset) | [`params`] |
//! | reading and sketching a pod5, the candidate overlaps | [`pipeline`] |
//! | the pipeline programs' candidate cache | [`pipeline::cache`] |
//! | the pipeline programs' strict command line | [`pipeline::cliargs`] |
//! | `haricot-detect`'s flags, graph build and layout outputs | [`pipeline::fbg_build`] |
//! | the seed-occurrence table (stage-1 features, contig-end rescue) | [`pipeline::occtable`] |
//! | the keep cut: drop lists, depth estimate, the depth-adaptive keep rule | [`pipeline::keep_cut`] |
//! | the fill: each candidate's path straightened | [`pipeline::fill`] |
//! | the end bridges | [`pipeline::endbridge`] |
//! | the contig-end rescue's admission | [`pipeline::admit`] |
//! | the full assembler, `haricot assemble` | [`mod@assemble`] |
//!
//! Signal-processing defaults follow the RawHash2 **R10.4.1** preset (see
//! [`params`]). Ported code cites the C reference implementations by `file:line`, as
//! `extern/RawHash2/...` and `extern/wtdbg2/...` (clone them there to follow a citation;
//! they are not distributed), or by file name alone (`kbm.h`, `wtdbg.c`, ...).

/// Process-wide allocator.
///
/// Declared in the library rather than in `main.rs` deliberately: the pipeline programs
/// (`src/bin/`) are separate binaries, and a `#[global_allocator]` in a linked
/// library applies to the whole program, so this one declaration covers every program
/// instead of each having to remember to opt in.
///
/// It is a pure performance choice with no effect on results. The overlap stage is
/// 91–99% of wall time and its dominant phase (the per-pair alignment, 70.6% on a 30x E. coli set)
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
