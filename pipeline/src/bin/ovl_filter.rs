//! Drop overlaps from a cache by READ or by PAIR, and report what went; dump a cache as TSV.
//!
//! `haricot assemble` applies the forests' drop lists this way (`--drop-pairs`, in
//! `haricot-layout`). As a tool it also answers whether layout damage is a property of
//! particular READS (a chimera poisons every overlap it is in) or of particular PAIRS.
//!
//! Usage:
//!   haricot-ovl-filter <in.bin> [--out <out.bin>] [--drop-reads F] [--drop-pairs F] [--min-z F]
//!                      [--max-overhang-bins N] [--scores F --min-score X] [--dump F]
//!                      [--dump-canonical F] [--depth]
//!                      [--read-model [--read-model-tsv F] [--agree REF.bin]]
//!
//! `--drop-reads` takes one read NAME per line and removes every overlap touching one.
//! `--drop-pairs` takes `name_a<TAB>name_b` per line, in either order. `--dump` writes the
//! surviving overlaps as a TSV keyed by read NAME, including the aligned spans and each
//! read's event count -- everything an OVERHANG needs, which `Overlap` alone cannot give
//! because it does not carry read lengths. `--dump-canonical` writes them as a function of the
//! read set alone (`canonical_dump`), for the exactness checks across read orders.
//!
//! `--max-overhang-bins` is wtdbg2's `is_dovetail_overlap` (`wtdbg.h:435-456`, default
//! `max_overhang = 256 bp = 1 bin`), which the overlapper does not apply. It is applied here
//! because the cache carries the aligned bin spans and `counts` carries each read's event count.
//! `--scores` applies a model's per-pair score from a TSV (`a b ... score`).
//!
//! `--read-model` applies the evidence-based read model (`haricot::overlap::read_model`,
//! wtdbg2's `clip_read_core`) to the
//! overlaps the other filters kept: each read is cut at the junctions its partners show (no
//! partner across, partners stopping on one side and starting on the other) and keeps only its
//! best-supported segment's overlaps. `--read-model-tsv` writes one row per cut read;
//! `--agree REF.bin` counts how many of the dropped pairs a reference cache (e.g. one built from
//! the truth) also lacks.
//!
//! `--depth` prints the label-free depth estimate of the surviving overlaps (`%.1f`): per read
//! bin, 1 + the overlaps covering it; each read's median over its bins; the median over reads
//! (a median of an even count is the mean of the two middle values). `haricot assemble`'s policy7
//! reads it (`haricot-layout`, W/depth.txt).
//!
//! The filters, the dump, the depth estimate and the read model live in `src/lo.rs`, shared with
//! `haricot-layout`.
use anyhow::{Context, Result};
use haricot::overlap::Overlap;

#[path = "../cache.rs"]
mod cache;
#[path = "../lo.rs"]
mod lo;

fn opt(a: &[String], n: &str) -> Option<String> {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .cloned()
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let in_path = args.get(1).context(
        "usage: haricot-ovl-filter <in.bin> [--out OUT.bin] [--dump F] [--dump-canonical F] \
         [--depth] [filters]",
    )?;
    let out_path = opt(&args, "--out");
    let dump_path = opt(&args, "--dump");
    let canon_path = opt(&args, "--dump-canonical");
    let depth = args.iter().any(|a| a == "--depth");
    anyhow::ensure!(
        out_path.is_some() || dump_path.is_some() || canon_path.is_some() || depth,
        "--out, --dump, --dump-canonical or --depth is required"
    );
    // `content` (the reads' signal content hashes, when the cache has them) is passed through to
    // `--out` and keys `--dump-canonical`'s frames.
    let cache::Cache {
        counts,
        names,
        ovl,
        content,
    } = cache::read(in_path)?;
    let idx = lo::name_index(&names);
    let mut f = lo::Filters::default();

    if let Some(p) = opt(&args, "--drop-reads") {
        let txt = std::fs::read_to_string(&p).with_context(|| format!("reading {p}"))?;
        let (reads, hit, miss) = lo::parse_drop_reads(&txt, &idx);
        f.drop_reads = reads;
        eprintln!("--drop-reads {p}: {hit} reads matched, {miss} not in the cache");
    }

    if let Some(p) = opt(&args, "--drop-pairs") {
        let txt = std::fs::read_to_string(&p).with_context(|| format!("reading {p}"))?;
        let (pairs, hit, miss) = lo::parse_drop_pairs(&txt, &idx);
        f.drop_pairs = pairs;
        eprintln!("--drop-pairs {p}: {hit} pairs matched, {miss} not in the cache");
    }

    f.min_z = opt(&args, "--min-z")
        .and_then(|v| v.parse().ok())
        .unwrap_or(f32::NEG_INFINITY);
    f.max_oh = opt(&args, "--max-overhang-bins")
        .and_then(|v| v.parse().ok())
        .unwrap_or(-1);
    f.min_score = opt(&args, "--min-score")
        .and_then(|v| v.parse().ok())
        .unwrap_or(f64::NEG_INFINITY);
    if let Some(p) = opt(&args, "--scores") {
        let txt = std::fs::read_to_string(&p).with_context(|| format!("reading {p}"))?;
        let mut hdr = true;
        for line in txt.lines() {
            let c: Vec<&str> = line.split('\t').collect();
            if hdr {
                hdr = false;
                if c.first() == Some(&"a") {
                    continue;
                }
            }
            let (Some(&i), Some(&j)) = (idx.get(c[0]), idx.get(c[1])) else {
                continue;
            };
            if let Some(s) = c.last().and_then(|v| v.parse().ok()) {
                f.scores.insert(if i < j { (i, j) } else { (j, i) }, s);
            }
        }
        eprintln!(
            "--scores {p}: {} pairs scored, cut {}",
            f.scores.len(),
            f.min_score
        );
    }

    let before = ovl.len();
    let mut t = lo::Tally::default();
    let out: Vec<Overlap> = ovl
        .into_iter()
        .filter(|o| f.keep(o, &counts, &mut t))
        .collect();
    eprintln!("{}", t.line(before, out.len()));
    let out = if args.iter().any(|a| a == "--read-model") {
        let tsv = opt(&args, "--read-model-tsv");
        let agree = opt(&args, "--agree");
        lo::read_model_apply(out, &counts, &names, tsv.as_deref(), agree.as_deref())?
    } else {
        out
    };
    if let Some(p) = dump_path {
        std::fs::write(&p, lo::dump_tsv(&out, &counts, &names))
            .with_context(|| format!("writing {p}"))?;
        eprintln!("dumped {} overlaps to {p}", out.len());
    }
    if let Some(p) = canon_path {
        std::fs::write(
            &p,
            canonical_dump(&out, &counts, &names, content.as_deref())?,
        )
        .with_context(|| format!("writing {p}"))?;
        eprintln!(
            "dumped {} overlaps in their canonical frames to {p}",
            out.len()
        );
    }
    if depth {
        let refs: Vec<&Overlap> = out.iter().collect();
        println!("{}", lo::py_fixed1(lo::depth_estimate(&refs, &counts)));
    }
    if let Some(p) = out_path {
        cache::write(&p, &counts, &names, &out, content.as_deref())?;
    }
    Ok(())
}

/// `--dump-canonical`: the overlaps as a function of the read set alone, for the exactness checks
/// across read orders (`checks/sym_compare.py`). Every overlap in its canonical pair frame
/// (`readorder::in_canonical_frame`: the longer read as the query `q`, the other `t`; equal lengths
/// by the tie rank, the reads' content-hash order, as the pipeline; a cache without content hashes,
/// by read id), keyed by the two read names, rows sorted by them; every field (`slope` to 17
/// significant digits, `z` exact) and the path cells (`q:t`, in the frame's bins).
fn canonical_dump(
    ovl: &[Overlap],
    counts: &[u32],
    names: &[String],
    content: Option<&[u64]>,
) -> Result<String> {
    let rank = match content {
        Some(c) => haricot::readorder::rank_by_content(c)?,
        None => haricot::readorder::rank_by_name(names),
    };
    let mut rows: Vec<(&str, &str, String)> = ovl
        .iter()
        .map(|o| {
            let c = haricot::readorder::in_canonical_frame(o, counts, &rank);
            let path: Vec<String> = c.path.iter().map(|(q, t)| format!("{q}:{t}")).collect();
            let (q, t) = (names[c.a as usize].as_str(), names[c.b as usize].as_str());
            let row = format!(
                "{q}\t{t}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.16e}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                c.z,
                c.n_anchors,
                c.matched,
                c.cnt,
                c.aligned_bins,
                c.gap,
                c.bin_events,
                c.slope,
                c.a_bin_span.0,
                c.a_bin_span.1,
                c.b_bin_span.0,
                c.b_bin_span.1,
                counts[c.a as usize],
                counts[c.b as usize],
                path.join(",")
            );
            (q, t, row)
        })
        .collect();
    rows.sort_unstable_by(|x, y| (x.0, x.1).cmp(&(y.0, y.1)));
    let mut s = String::from(
        "q\tt\tz\tn_anchors\tmatched\tcnt\taligned_bins\tgap\tbin_events\tslope\
         \tq_lo\tq_hi\tt_lo\tt_hi\tq_ev\tt_ev\tpath\n",
    );
    for (_, _, r) in rows {
        s.push_str(&r);
    }
    Ok(s)
}
