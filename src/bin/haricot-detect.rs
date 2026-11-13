//! `haricot-detect`: candidate detection, the `detect` step of `haricot assemble` (Algorithm S2,
//! part 1).
//!
//! From a pod5: segment each read's signal into events, sketch the events into e-mer seeds, and
//! find the candidate overlaps by the Hough transform over shared seeds (`haricot::overlap`);
//! `--save-overlaps` writes them as a candidate cache (`pipeline::cache`) for the later steps, and
//! `<outdir>` gets the read tables (read_names.tsv, read_content.tsv).
//!
//! Usage:
//!   haricot-detect <reads.pod5> --production [--outdir DIR] [--save-overlaps FILE] --no-graph
//! Unknown `--` flags are an error. `--production` stands for the detection and graph flags of
//! `haricot assemble` (`fbg_build::DETECT_FLAGS`, `fbg_build::GRAPH_FLAGS`); its first step is
//! `haricot-detect <reads.pod5> --production --save-overlaps FILE --no-graph`. Without it the defaults
//! are `Params::r10_4_1()`'s. `--no-graph` is accepted for that command line: the program never
//! builds a graph.
//!
//! The flags live in `src/pipeline/fbg_build.rs`, shared with `haricot-layout`.

use std::fs;
use std::time::Instant;

use anyhow::{Context, Result};

use haricot::{overlap, pipeline};

use haricot::pipeline::{cache, fbg_build};

use fbg_build::{Opts, prof};

fn main() -> Result<()> {
    let t_total = Instant::now();
    let Opts {
        input,
        outdir,
        params,
        save_overlaps,
        ..
    } = fbg_build::parse(std::env::args().skip(1))?;
    fs::create_dir_all(&outdir)?;

    // Overlaps depend only on the sketch and overlap parameters, not on the graph's, so the
    // later steps read them from the cache. `counts`/`names` follow sketch (= read) order, so a
    // node's anchor read index resolves to its name. `content`: every read's signal content
    // hash, the key of the tie rank (`readorder`).
    let (counts, names, overlaps, content) = {
        let input = input.as_ref().context("need <reads.pod5>")?;
        let t = Instant::now();
        // Read and sketch in ONE streamed pass. Loading every read first makes the whole
        // signal and the whole seed store live at the same time, which is the
        // largest single term in peak RSS (32.65 of 49.40 GB on a D. melanogaster set,
        // 7.97 of 11.30 GB on a 30x E. coli set). Chunking bounds the signal's share
        // to a constant. 4 GiB of signal in flight: faster than 1 GiB chunks on every
        // dataset measured, because sketching a chunk whose signal is still cache-resident
        // more than repays the lost read-ahead once chunks are large.
        let sketches = {
            let s = pipeline::sketch_pod5_streaming(input, None, 4 << 30, &params)
                .with_context(|| format!("reading {}", input.display()))?;
            let n: usize = s.iter().map(|x| x.seeds.len()).sum();
            prof(
                "read_sketch",
                t,
                &format!("reads={} seeds={n} (streamed)", s.len()),
            );
            s
        };
        // Every later step keys reads by id, and the canonical tie-breaks rank them by it: a
        // repeated id is an input error.
        haricot::readorder::ensure_distinct_ids(sketches.iter().map(|s| s.name.as_str()))
            .with_context(|| format!("reading {}", input.display()))?;
        // ... and so is a repeated SIGNAL: the tie rank orders reads by their content hash,
        // which must be a total order (an identical signal is a duplicated read).
        let hashes: Vec<u64> = sketches.iter().map(|s| s.content).collect();
        if let Some((x, y)) = haricot::readorder::content_duplicate(&hashes) {
            anyhow::bail!(
                "reading {}: reads {} and {} have the same signal content hash {:016x} (an \
                     identical signal, i.e. a duplicated read, or a 64-bit collision)",
                input.display(),
                sketches[x].name,
                sketches[y].name,
                hashes[x]
            );
        }
        // Written BEFORE the overlap stage on purpose: a process that reads the overlap
        // stage's stderr on a pipe can map read numbers to names while detection runs
        // (written after it, the reader would wait on the file while the overlap stage
        // blocks on a full pipe).
        {
            use std::io::Write as _;
            let path = outdir.join("read_names.tsv");
            let f =
                fs::File::create(&path).with_context(|| format!("writing {}", path.display()))?;
            let mut w = std::io::BufWriter::new(f);
            for (i, sk) in sketches.iter().enumerate() {
                writeln!(w, "{i}\t{}", sk.name)?;
            }
        }
        let t = Instant::now();
        let overlaps = pipeline::overlaps(&sketches, &params);
        prof(
            "overlap",
            t,
            &format!(
                "overlaps={} kernel={}",
                overlaps.len(),
                overlap::kernel_name()
            ),
        );
        // Coordinate space (merged events), NOT raw n_events: bins/anchors/K-bin
        // windows are all indexed in the retained stream. See
        // haricot::emer::retained_quantized.
        let counts: Vec<u32> = sketches.iter().map(|s| s.n_coord_events).collect();
        let names: Vec<String> = sketches.iter().map(|s| s.name.clone()).collect();
        let content: Vec<u64> = sketches.iter().map(|s| s.content).collect();
        (counts, names, overlaps, Some(content))
    };
    fbg_build::write_read_tables(&outdir, &names, content.as_deref())?;

    if let Some(path) = &save_overlaps {
        let t = Instant::now();
        let bytes = cache::encode(&counts, &names, &overlaps, content.as_deref())
            .context("encoding overlap cache")?;
        fs::write(path, bytes).with_context(|| format!("writing {}", path.display()))?;
        prof("save_overlaps", t, &format!("overlaps={}", overlaps.len()));
        eprintln!(
            "saved {} overlaps ({} reads) to {}",
            overlaps.len(),
            counts.len(),
            path.display()
        );
    }
    fbg_build::overlap_stats(&overlaps, &params);
    prof("TOTAL", t_total, "end-to-end (--no-graph)");
    Ok(())
}
