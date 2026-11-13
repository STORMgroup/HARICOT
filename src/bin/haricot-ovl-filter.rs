//! `haricot-ovl-filter`: dump a candidate cache as a function of the read set alone, for the
//! exactness checks across read orders and read names (`checks/invariance.sh`). Not a step of
//! `haricot assemble`.
//!
//! Usage:
//!   haricot-ovl-filter <in.bin> --dump-canonical F
//!
//! `--dump-canonical` writes the cache's overlaps, each in its canonical pair frame and keyed by
//! read names, to F (`canonical_dump`).
use anyhow::{Context, Result};
use haricot::overlap::Overlap;

use haricot::pipeline::cache;

/// The value after flag `n` in `a`.
fn opt(a: &[String], n: &str) -> Option<String> {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .cloned()
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let usage = "usage: haricot-ovl-filter <in.bin> --dump-canonical F";
    let in_path = args.get(1).context(usage)?;
    let p = opt(&args, "--dump-canonical").context(usage)?;
    // `content` (the reads' signal content hashes) keys `--dump-canonical`'s frames.
    let cache::Cache {
        counts,
        names,
        ovl,
        content,
    } = cache::read(in_path)?;
    std::fs::write(
        &p,
        canonical_dump(&ovl, &counts, &names, content.as_deref())?,
    )
    .with_context(|| format!("writing {p}"))?;
    eprintln!(
        "dumped {} overlaps in their canonical frames to {p}",
        ovl.len()
    );
    Ok(())
}

/// `--dump-canonical`: the overlaps as a function of the read set alone, for the exactness checks
/// across read orders (`checks/sym_compare.py`). Every overlap in its canonical pair frame
/// (`readorder::in_canonical_frame`: the longer read as the query `q`, the other `t`; equal lengths
/// by the tie rank, the reads' content-hash order, as the pipeline), keyed by the two read names,
/// rows sorted by them; every field (`slope` to 17 significant digits, `z` exact) and the path cells
/// (`q:t`, in the frame's bins).
fn canonical_dump(
    ovl: &[Overlap],
    counts: &[u32],
    names: &[String],
    content: Option<&[u64]>,
) -> Result<String> {
    let rank = cache::tie_rank(content)?;
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
