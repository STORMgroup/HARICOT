//! Command-line interface.
//!
//! - `assemble` — the full assembler (reads -> layout), as run for the manuscript ([`crate::assemble`])

use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(
    name = "haricot",
    version,
    about = "HARICOT: genome assembly from raw nanopore signal (R10.4.1, forward-strand reads)",
    after_help = "Assemble a POD5 file:\n  haricot assemble reads.pod5 -o out\n\nRun `haricot assemble --help` for its options."
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Assemble a POD5 file: the full pipeline, as run for the manuscript.
    ///
    /// Writes OUT/unitigs.gfa (the layout: unitigs as paths of read segments), OUT/contigs.gfa, the path and node
    /// tables, OUT/threads.tsv (read placements) and OUT/time.tsv (wall time and peak memory per step); step logs
    /// in OUT/work/.
    #[command(after_help = "Examples:\n  haricot assemble reads.pod5 -o asm\n  \
        haricot assemble reads.pod5 -o asm -t 16 --models path/to/forests\n\n\
        Input: R10.4.1 reads, 5 or 4 kHz. Strands are not merged: reads of the two strands assemble separately.\n\
        The run log states whether the forests read are the published ones (SHA-256).")]
    Assemble {
        /// Input reads (.pod5).
        reads: PathBuf,
        /// Output directory (must not exist, or be empty, or hold an earlier run with --force).
        #[arg(short, long)]
        out: PathBuf,
        /// Directory with the forests (stage1.*, stage2.*, rescue.*).
        #[arg(short, long, default_value = crate::assemble::DEFAULT_MODELS)]
        models: PathBuf,
        /// Threads per step [default: all available cores].
        #[arg(short, long)]
        threads: Option<usize>,
        /// Keep the intermediate files (candidates, features, forest scores) in OUT/work.
        #[arg(long)]
        keep_work: bool,
        /// Replace OUT if it holds an earlier run.
        #[arg(long)]
        force: bool,
    },
}

/// Parse args and dispatch. Called by `main`.
pub fn run() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Command::Assemble {
            reads,
            out,
            models,
            threads,
            keep_work,
            force,
        } => crate::assemble::run(&crate::assemble::Options {
            reads,
            out,
            models,
            threads,
            keep_work,
            force,
        }),
    }
}
