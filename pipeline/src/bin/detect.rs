//! Candidate detection, and the graph build and layout from a candidate set.
//!
//! From a pod5: sketch the reads, find the candidate overlaps, `--save-overlaps` them as a cache, and
//! (without `--no-graph`) build the graph. From a cache (`--load-overlaps`): build the graph and layout
//! only. The build writes `<outdir>/fbg.gfa` and `<outdir>/fbg_nodes.tsv` (node, anchor_read, ev_beg,
//! ev_end, cov, anchor_coord_events), and with `--unitig` / `--dump-paths` the unitig and contig files.
//!
//! Usage:
//!   haricot-detect <reads.pod5> [--outdir DIR] [overlap flags] [--save-overlaps FILE] [--no-graph]
//!   haricot-detect --load-overlaps FILE [--outdir DIR] [graph flags] [--clean] [--unitig]
//!       [--dump-paths] [--frg-layout]
//! Unknown `--` flags are an error. `--production` stands for the detection and graph flags of
//! `haricot assemble` (`fbg_build::DETECT_FLAGS`, `fbg_build::GRAPH_FLAGS`); its first step is
//! `haricot-detect <reads.pod5> --production --save-overlaps FILE --no-graph`. Without it the defaults
//! are `Params::r10_4_1()`'s.
//!
//! The flags, the graph build and the layout outputs live in `src/fbg_build.rs`, shared with
//! `haricot-layout`.

use std::fs;
use std::time::Instant;

use anyhow::{Context, Result};

use haricot::{overlap, pipeline, signal};

#[path = "../cache.rs"]
mod cache;
#[path = "../fbg_build.rs"]
mod fbg_build;

use fbg_build::{Opts, prof};

fn main() -> Result<()> {
    let t_total = Instant::now();
    let Opts {
        input,
        outdir,
        params,
        layout: lay,
        no_graph,
        ovl_min_z,
        save_overlaps,
        load_overlaps,
        max_reads,
        canonical_order,
        stream_sketch,
        stream_chunk_bytes,
    } = fbg_build::parse(std::env::args().skip(1))?;
    fs::create_dir_all(&outdir)?;

    // Overlaps depend only on the sketch/overlap params, NOT on the graph/coord
    // params — so cache them once and reuse across build-parameter iterations,
    // which is the expensive (all-vs-all) stage. `counts`/`names` mirror sketch
    // (= read) order so node provenance (anchor_read index) still resolves.
    // `content`: every read's signal content hash, the key of the tie rank (`readorder`); `None` for a
    // cache without them.
    let (counts, names, mut overlaps, content) = if let Some(path) = &load_overlaps {
        let bytes = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        let c = cache::decode(&bytes).context("decoding overlap cache")?;
        haricot::readorder::ensure_distinct_ids(c.names.iter().map(String::as_str))
            .with_context(|| format!("reading {}", path.display()))?;
        (c.counts, c.names, c.ovl, c.content)
    } else {
        let input = input
            .as_ref()
            .context("need <reads.pod5> (or --load-overlaps CACHE)")?;
        let t = Instant::now();
        // Read and sketch in ONE streamed pass unless --no-stream-sketch selects
        // the two-pass form. Loading every read first makes the whole
        // signal and the whole seed store live at the same time, which is the
        // largest single term in peak RSS (32.65 of D. melanogaster fwd's 49.40 GB,
        // 7.97 of E. coli fwd-30x's 11.30 GB). Chunking bounds the signal's share
        // to a constant. The sketches are identical either way — the pod5 iterator
        // is consumed once in order, so reads, order and `read_id`s all match.
        let mut sketches = if stream_sketch {
            let s = pipeline::sketch_pod5_streaming(input, max_reads, stream_chunk_bytes, &params)
                .with_context(|| format!("reading {}", input.display()))?;
            let n: usize = s.iter().map(|x| x.seeds.len()).sum();
            prof(
                "read_sketch",
                t,
                &format!("reads={} seeds={n} (streamed)", s.len()),
            );
            s
        } else {
            let reads = signal::read_signals_limit(input, max_reads)
                .with_context(|| format!("reading {}", input.display()))?;
            prof("read_signals", t, &format!("reads={}", reads.len()));
            let t = Instant::now();
            let s = pipeline::sketch_reads(&reads, &params);
            let n: usize = s.iter().map(|x| x.seeds.len()).sum();
            // The raw signal is dead once the sketch exists; holding it through the
            // overlap stage was pure peak-RSS tax.
            drop(reads);
            prof("sketch", t, &format!("reads={} seeds={n}", s.len()));
            s
        };
        // Every later step keys reads by id, and the canonical tie-breaks rank them by it: a
        // repeated id is an input error under any numbering, not only --canonical-order's.
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
        if canonical_order {
            haricot::readorder::canonicalize(&mut sketches)?;
            eprintln!(
                "canonical read order: {} reads numbered by read id",
                sketches.len()
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

    // Offline re-thresholding. Applied before the save so a filtered subset can be
    // materialised, and before the graph so one permissive pass serves every operating
    // point. `NaN` (synthesized overlaps) fails `<` and therefore survives.
    if let Some(min_z) = ovl_min_z {
        let before = overlaps.len();
        overlaps.retain(|o| o.z >= min_z || o.z.is_nan());
        eprintln!(
            "--ovl-min-z {min_z}: kept {} of {before} overlaps ({:.1}%)",
            overlaps.len(),
            100.0 * overlaps.len() as f64 / before.max(1) as f64
        );
    }

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
    if no_graph {
        prof("TOTAL", t_total, "end-to-end (--no-graph)");
        return Ok(());
    }
    let reads = fbg_build::Reads {
        counts: &counts,
        names: &names,
        content: content.as_deref(),
    };
    fbg_build::build(&outdir, &params.graph, &lay, &reads, &overlaps, None)?;
    prof("TOTAL", t_total, "end-to-end");
    Ok(())
}
