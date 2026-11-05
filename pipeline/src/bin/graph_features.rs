//! Stage-2 features: what the ACCEPTED graph says about each accepted overlap.
//!
//! Stage 1 (`haricot-forest --stage 1`) scores each permissive candidate on its own evidence and its
//! read context. What is left after it are mostly good alignments to the wrong copy of a repeat
//! (on D. melanogaster, ~150 false candidates behind every misjoin). Those are visible only against
//! the other accepted overlaps: at a true locus, the reads that overlap A where A meets B also overlap B,
//! at consistent offsets; at a repeat copy, A's partners over that region split into two groups
//! that do not overlap each other.
//!
//! For every candidate with a stage-1 score >= `--t1` (the stage-1 accepted set), over that set:
//!   mut_a / mut_b    share of A's (B's) other accepted partners over the shared region that are
//!                    also accepted partners of B (A)
//!   cons_a / cons_b  share of those whose offsets compose: off(A->C) + off(C->B) ~ off(A->B)
//!   nbr_a / nbr_b    how many such partners, over the dataset median
//!   nb_score         mean stage-1 score of the accepted overlaps around the pair
//!   rank_a / rank_b  rank of this pair's stage-1 score among A's (B's) partners over the region
//!
//! Usage: haricot-graph-features CAND.bin SCORES.tsv --t1 T --out F.tsv [--tol 576] [--no-symmetric]
//!   every pair in its canonical frame (the longer read as `a`; see the code), the rows in the
//!   cache's names and order: the default (`--symmetric` is accepted); `--no-symmetric`: every pair
//!   in the cache's orientation
//!   SCORES.tsv: name_a name_b cls score (haricot-forest's or training/train_stage1.py's
//!   *.scores.tsv / *.oob.tsv)
//!   An unknown flag, a flag without its value and a value that does not parse are errors
//!   (`cliargs.rs`).
use std::collections::HashMap;
use std::io::Write;

use anyhow::{Context, Result};
use haricot::overlap::Overlap;
use rayon::prelude::*;

#[path = "../cache.rs"]
mod cache;
#[path = "../cliargs.rs"]
mod cliargs;

/// The flags that take a value, and the switches (`cliargs.rs`).
const VALUED: &[&str] = &["--t1", "--out", "--tol"];
const SWITCHES: &[&str] = &["--symmetric", "--no-symmetric"];

fn theil(pts: &[(f64, f64)]) -> Option<(f64, f64)> {
    if pts.len() < 2 {
        return None;
    }
    let mut sl = Vec::with_capacity(pts.len() * (pts.len() - 1) / 2);
    for i in 0..pts.len() {
        for j in i + 1..pts.len() {
            let dx = pts[j].0 - pts[i].0;
            if dx.abs() > 1e-9 {
                sl.push((pts[j].1 - pts[i].1) / dx);
            }
        }
    }
    if sl.is_empty() {
        return None;
    }
    sl.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let s = sl[sl.len() / 2];
    let mut ic: Vec<f64> = pts.iter().map(|p| p.1 - s * p.0).collect();
    ic.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some((s, ic[ic.len() / 2]))
}

/// One accepted overlap as seen from one of its reads.
#[derive(Clone, Copy)]
struct Half {
    partner: u32,
    lo: u32,
    hi: u32,
    /// Where this read's event 0 lands in the partner's coordinates.
    origin_in_partner: f64,
    score: f32,
}

fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let usage = "usage: haricot-graph-features CAND.bin SCORES.tsv --t1 T --out F.tsv";
    let args = cliargs::Cli::parse(&argv, VALUED, SWITCHES, usage)?;
    args.positional(2, usage)?;
    let (cand_p, sc_p) = (&args.pos[0], &args.pos[1]);
    let t1: f32 = args.opt("--t1").context("--t1")?.parse()?;
    let out_p = args.opt("--out").context("--out")?;
    let tol: f64 = args.num("--tol", 576.0)?;

    let cache::Cache {
        counts,
        names,
        ovl: cand,
        content,
    } = cache::read(cand_p)?;
    // The canonical pair frame (the default; `--no-symmetric` for the owner's): every candidate in its
    // canonical pair frame (`readorder::query_first`: the longer read as `a`), so the `_a` / `_b`
    // columns are the longer / shorter read's and no value depends on which read owns the pair; each
    // read's accepted partners are visited in tie-rank order, so the sums do not depend on the
    // numbering either. The rows keep the cache's names and order.
    let symmetric = args.symmetric()?;
    // the tie rank (the reads' signal content-hash order); used, and required, only in that frame
    let rank = if symmetric {
        cache::tie_rank(content.as_deref())?
    } else {
        Vec::new()
    };
    let cand_c: Vec<Overlap> = if symmetric {
        cand.iter()
            .map(|o| haricot::readorder::in_canonical_frame(o, &counts, &rank))
            .collect()
    } else {
        Vec::new()
    };
    let fc: &[Overlap] = if symmetric { &cand_c } else { &cand };
    let bin = cand.first().map(|o| o.bin_events).unwrap_or(192);
    let idx: HashMap<&str, u32> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect();
    let mut score: HashMap<(u32, u32), f32> = HashMap::new();
    for line in std::fs::read_to_string(sc_p)?.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if let (Some(&a), Some(&b), Ok(s)) = (idx.get(f[0]), idx.get(f[1]), f[3].parse::<f32>()) {
            score.insert((a.min(b), a.max(b)), s);
        }
    }

    // The accepted set, as halves per read.
    let lines: Vec<Option<(f64, f64)>> = fc
        .par_iter()
        .map(|o| {
            let pts: Vec<(f64, f64)> = o
                .path
                .iter()
                .map(|&(a, b)| ((a * bin + bin / 2) as f64, (b * bin + bin / 2) as f64))
                .collect();
            theil(&pts)
        })
        .collect();
    let mut halves: Vec<Vec<Half>> = vec![Vec::new(); names.len()];
    let mut accepted: HashMap<(u32, u32), usize> = HashMap::new();
    for (k, o) in fc.iter().enumerate() {
        let pk = (o.a.min(o.b), o.a.max(o.b));
        let Some(&s) = score.get(&pk) else {
            continue;
        };
        let Some((sl, c)) = lines[k] else { continue };
        if s < t1 {
            continue;
        }
        accepted.insert(pk, k);
        halves[o.a as usize].push(Half {
            partner: o.b,
            lo: o.a_bin_span.0,
            hi: o.a_bin_span.1,
            origin_in_partner: c,
            score: s,
        });
        halves[o.b as usize].push(Half {
            partner: o.a,
            lo: o.b_bin_span.0,
            hi: o.b_bin_span.1,
            origin_in_partner: if sl.abs() > 1e-6 { -c / sl } else { 0.0 },
            score: s,
        });
    }
    if symmetric {
        for h in halves.iter_mut() {
            h.sort_by_key(|x| rank[x.partner as usize]);
        }
    }
    let is_acc = |x: u32, y: u32| accepted.contains_key(&(x.min(y), x.max(y)));
    let origin = |x: u32, y: u32| -> Option<f64> {
        halves[x as usize]
            .iter()
            .find(|h| h.partner == y)
            .map(|h| h.origin_in_partner)
    };
    // Median number of region-partners, for normalisation.
    let region = |r: u32, lo: u32, hi: u32, skip: u32| -> Vec<Half> {
        #[allow(clippy::manual_div_ceil)] // half the interval's length, floored
        let need = ((hi - lo + 1) / 2).max(1);
        halves[r as usize]
            .iter()
            .copied()
            .filter(|h| h.partner != skip && h.hi.min(hi) + 1 >= h.lo.max(lo) + need)
            .collect()
    };
    // Over EVERY accepted candidate. A sample (say the first 20,000 of `accepted` in the iteration
    // order of a randomly seeded HashMap) would give a different median and
    // different `nbr_*` features on every run of an input with more accepted candidates.
    let keys: Vec<usize> = accepted.values().copied().collect();
    let mut nsz: Vec<u32> = keys
        .par_iter()
        .map(|&k| {
            let (a, lo, hi, b) = (fc[k].a, fc[k].a_bin_span.0, fc[k].a_bin_span.1, fc[k].b);
            #[allow(clippy::manual_div_ceil)] // as in `region`
            let need = ((hi - lo + 1) / 2).max(1);
            halves[a as usize]
                .iter()
                .filter(|h| h.partner != b && h.hi.min(hi) + 1 >= h.lo.max(lo) + need)
                .count() as u32
        })
        .collect();
    nsz.sort_unstable();
    let nmed = nsz.get(nsz.len() / 2).copied().unwrap_or(1).max(1) as f64;
    eprintln!(
        "{} candidates, {} accepted at t1 {t1}; median region partners {nmed}",
        cand.len(),
        accepted.len()
    );

    // Rows in candidate order. `accepted` is a std HashMap, whose iteration order is random per
    // process, and the stage-2 forest bootstraps by row position: an unordered file made every
    // retrain a different forest whatever its random_state.
    let mut order: Vec<((u32, u32), usize)> = accepted.iter().map(|(&p, &k)| (p, k)).collect();
    order.sort_unstable_by_key(|&(_, k)| k);
    let rows: Vec<String> = order
        .par_iter()
        .map(|&((ka, kb), k)| {
            // the features' frame (canonical unless --no-symmetric); the row keeps the cache's names
            let o = &fc[k];
            let (a, b) = (o.a, o.b);
            let s0 = score[&(ka, kb)];
            let off_ab = lines[k].map(|x| x.1).unwrap_or(0.0);
            let side = |x: u32, y: u32, lo: u32, hi: u32, off_xy: f64| -> (f64, f64, f64, f64, f64) {
                let n = region(x, lo, hi, y);
                if n.is_empty() {
                    return (f64::NAN, f64::NAN, 0.0, f64::NAN, 0.0);
                }
                let (mut m, mut c) = (0usize, 0usize);
                let mut sum = 0.0f64;
                let mut higher = 0usize;
                for h in &n {
                    sum += h.score as f64;
                    if h.score > s0 {
                        higher += 1;
                    }
                    if is_acc(y, h.partner) {
                        m += 1;
                        // off(x->c) + off(c->y) ~ off(x->y)
                        if let Some(cy) = origin(h.partner, y)
                            && (h.origin_in_partner + cy - off_xy).abs() <= tol
                        {
                            c += 1;
                        }
                    }
                }
                let nn = n.len() as f64;
                (m as f64 / nn, c as f64 / nn, nn / nmed, sum / nn, higher as f64 / nn)
            };
            let (ma, ca, na, sa, ra) = side(a, b, o.a_bin_span.0, o.a_bin_span.1, off_ab);
            let off_ba = origin(b, a).unwrap_or(0.0);
            let (mb, cb, nb, sb, rb) = side(b, a, o.b_bin_span.0, o.b_bin_span.1, off_ba);
            let nbs = if sa.is_nan() {
                sb
            } else if sb.is_nan() {
                sa
            } else {
                (sa + sb) / 2.0
            };
            format!(
                "{}\t{}\t{s0:.4}\t{ma:.4}\t{mb:.4}\t{ca:.4}\t{cb:.4}\t{na:.3}\t{nb:.3}\t{nbs:.4}\t{ra:.4}\t{rb:.4}",
                names[ka as usize], names[kb as usize]
            )
        })
        .collect();
    let mut w = std::io::BufWriter::new(std::fs::File::create(&out_p)?);
    writeln!(
        w,
        "name_a\tname_b\ts1\tmut_a\tmut_b\tcons_a\tcons_b\tnbr_a\tnbr_b\tnb_score\trank_a\trank_b"
    )?;
    for r in rows {
        writeln!(w, "{r}")?;
    }
    eprintln!("wrote {out_p}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The calls in the tree parse, all in the canonical frame: `haricot assemble`'s, the
    /// training script's, training/train_invariance.py's (`--symmetric`, stated).
    #[test]
    fn the_calls_in_the_tree_parse() {
        let cli = |l: &str| {
            let v: Vec<String> = l.split_whitespace().map(String::from).collect();
            cliargs::Cli::parse(&v, VALUED, SWITCHES, "usage")
        };
        for l in [
            "C S --t1 0.529462250692771 --out f.tsv",
            "C S --t1 0.529462250692771 --out f.tsv --symmetric",
        ] {
            let a = cli(l).unwrap();
            a.positional(2, "").unwrap();
            assert!(a.symmetric().unwrap());
            assert_eq!(a.num("--tol", 576.0).unwrap(), 576.0);
            assert!(a.opt("--t1").unwrap().parse::<f32>().is_ok());
        }
        assert!(
            !cli("C S --t1 0.5 --out f --no-symmetric")
                .unwrap()
                .symmetric()
                .unwrap()
        );
        assert!(cli("C S --t1 0.5 --out f --tol").is_err());
        assert!(
            cli("C S --t1 0.5 --out f --tol 5x6")
                .unwrap()
                .num("--tol", 576.0)
                .is_err()
        );
        assert!(cli("C S --t1 0.5 --out f --sym").is_err());
    }
}
