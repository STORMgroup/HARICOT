//! The contig-end rescue: a significance that can see a SHORT overlap in a LARGE anchor cloud.
//!
//! Three modes:
//!
//! * `haricot-rescue rescue` (fn `rescue`): the pool of a build's contig-end read pairs that the
//!   candidates do not hold and the per-offset scan finds significant (`haricot assemble` step 4);
//! * `haricot-rescue admit` (fn `admit`): the mutual best of the pool pairs the rescue forest
//!   accepts, added to a cache (in `haricot assemble`, `haricot-layout final` runs the same code);
//! * `haricot-rescue <reads.pod5> <pairs.tsv> ...`: the scan on listed pairs, as a tool.
//!
//! The detector's `z` (`dp::diagonal_significance_parts`) takes the densest offset window of the
//! whole cloud and compares it to the null maximum over all windows, with ONE expected count `E`
//! for the window at that peak. A dovetail between two long reads sits on a diagonal that crosses
//! only a corner of the read box: its window collects the overlap's anchors plus the few random
//! ones a short diagonal can hold, while the central diagonals, which cross the whole box, collect
//! many more random anchors. On E. coli 6x the true contig-end dovetails (1-6 kb between 16-45k-event
//! reads, 4-19k anchors) score `z` 0.4-3.6, indistinguishable from the 2,270 false end pairs, and
//! their fitted slopes (0.70-1.38) show the Hough fit never found their diagonal.
//!
//! Here every offset window gets its OWN expectation: its count is tested against the random
//! anchors a diagonal of that length can hold. Per slope in a small grid around 1:
//!   1. offsets `c = a - s*b` of every anchor, counted in windows of width `2*tol` sliding by `tol`;
//!   2. a cheap expectation per window, uniform density over the read box times the band area of
//!      that diagonal inside it (`2*tol * its length`), ranks the windows by Poisson tail;
//!   3. the best `--refine` windows get the exact marginal-preserving expectation the detector uses
//!      (`E = (1/N) #{(i,j) : |a_i - s*b_j - c| <= tol}`, `dp.rs`), and their Poisson tail
//!      (`dp::tail_nlp`, one window);
//!   4. the pair's score is the best window's `-log10 p`, minus `log10 M` for the M windows x slopes.
//!
//! The best window's geometry is reported as the overlap it implies: a dovetail a->b (a's suffix on
//! b's prefix), b->a, a containment, or an internal (chimeric) segment.
//!
//! Anchors are built as the detector builds them: seeds from the pipeline's sketch, hashes that
//! occur more than `--max-seed-occ` times in the whole input masked.
//!
//! With `--emit BASE.bin`, every listed pair (the admitted ones) is instead turned into an
//! [`Overlap`] along its best window's diagonal, clipped to the read box (a dovetail runs from one
//! read's end to the other's start, as the fill would draw it), appended to BASE's overlaps and
//! written to `--out` as a cache (`z` = NaN: not a detector significance).
//!
//! Usage: haricot-rescue <reads.pod5> <pairs.tsv> --out <scores.tsv | out.bin> [--emit BASE.bin]
//!                 [--max-seed-occ 2000] [--bin 192] [--slopes 0.93,1,1.07] [--refine 10]
//!                 [--end-margin M]
//!        haricot-rescue rescue ... (see fn `rescue`)
//!        haricot-rescue admit ... (see fn `admit`)
//! `pairs.tsv`: read indices in the first two columns (a header line is skipped).
//! In every mode an unknown flag, a flag without its value and a value that does not parse are
//! errors (`cliargs.rs`).
use anyhow::{Context, Result};
use haricot::emer;
use haricot::overlap::dp::tail_nlp;
use haricot::overlap::{Anchor, Overlap, SeedCounts};
use haricot::params::Params;
use haricot::params::Strand;
use haricot::pipeline;
use haricot::signal;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io::Write;

#[path = "../admit.rs"]
mod admit;
#[path = "../cache.rs"]
mod cache;
#[path = "../cliargs.rs"]
mod cliargs;
#[path = "../occtable.rs"]
mod occtable;

use cliargs::Cli;

/// Each mode's flags that take a value, and its switches (`cliargs.rs`).
const RESCUE_VALUED: &[&str] = &[
    "--out-dir",
    "--end-nodes",
    "--t",
    "--max-seed-occ",
    "--bin",
    "--slopes",
    "--refine",
    "--end-margin",
    "--occ",
];
const RESCUE_SWITCHES: &[&str] = &["--dovetails-only", "--symmetric", "--no-symmetric"];
const ADMIT_VALUED: &[&str] = &["--out", "--allow", "--out-pairs"];
const ADMIT_SWITCHES: &[&str] = &["--symmetric", "--no-symmetric"];
const SCORE_VALUED: &[&str] = &[
    "--out",
    "--max-seed-occ",
    "--bin",
    "--slopes",
    "--refine",
    "--end-margin",
    "--emit",
];

/// Best window of one pair at one slope: (nlp of the single window, offset centre, obs, E) and the
/// overlap shape its own anchors support.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Win {
    nlp: f64,
    c: f64,
    obs: u32,
    e: f64,
    shape: &'static str,
}

/// The shape a window's ANCHORS support: which read ends they reach, within `m` events. A real
/// overlap reaches ends (a dovetail two, a containment the contained read's both); anchors that stop
/// short of every end are a shared segment inside both reads (a repeat copy, a chimera's piece).
#[allow(clippy::too_many_arguments)]
fn anchor_shape(
    anchors: &[Anchor],
    la: f64,
    lb: f64,
    s: f64,
    c: f64,
    tol: f64,
    m: f64,
) -> &'static str {
    let (mut amin, mut amax, mut bmin, mut bmax) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for an in anchors {
        let (x, y) = (an.a_pos as f64, an.b_pos as f64);
        if (x - s * y - c).abs() <= tol {
            amin = amin.min(x);
            amax = amax.max(x);
            bmin = bmin.min(y);
            bmax = bmax.max(y);
        }
    }
    if amin > amax {
        return "none";
    }
    let (a0, a1, b0, b1) = (amin <= m, amax >= la - m, bmin <= m, bmax >= lb - m);
    if b0 && b1 {
        "b_in_a"
    } else if a0 && a1 {
        "a_in_b"
    } else if a1 && b0 {
        "dove_ab" // a's suffix on b's prefix: a comes first
    } else if a0 && b1 {
        "dove_ba"
    } else {
        "internal"
    }
}

/// Length (in b) of the diagonal `a = s*b + c` inside the box `[0, la] x [0, lb]`.
fn diag_len(la: f64, lb: f64, s: f64, c: f64) -> f64 {
    let lo = (-c / s).max(0.0);
    let hi = ((la - c) / s).min(lb);
    (hi - lo).max(0.0)
}

/// Exact marginal-preserving expectation for the window centred at `c` (as the detector's null).
fn marginal_e(anchors: &[Anchor], sb_sorted: &[f64], c: f64, tol: f64) -> f64 {
    let mut m = 0u64;
    for an in anchors {
        let center = an.a_pos as f64 - c;
        let l = sb_sorted.partition_point(|&v| v < center - tol);
        let r = sb_sorted.partition_point(|&v| v <= center + tol);
        m += (r - l) as u64;
    }
    m as f64 / anchors.len().max(1) as f64
}

/// The per-offset scan at slope `s`. Returns the best window whose anchors reach read ends (`m`
/// events; `m <= 0` accepts any window) and the number of windows scanned.
#[allow(clippy::too_many_arguments)]
fn scan(
    anchors: &[Anchor],
    la: f64,
    lb: f64,
    s: f64,
    tol: f64,
    refine: usize,
    m: f64,
) -> (Option<Win>, usize) {
    let n = anchors.len();
    if n == 0 {
        return (None, 0);
    }
    let w = 2.0 * tol;
    let c_min = -s * lb;
    let n_win = (((la - c_min) / tol).ceil() as usize).max(1);
    // counts per tol-wide cell; a window is two adjacent cells
    let mut cell = vec![0u32; n_win + 1];
    for an in anchors {
        let c = an.a_pos as f64 - s * an.b_pos as f64;
        let k = ((c - c_min) / tol).floor();
        if k >= 0.0 && (k as usize) < cell.len() {
            cell[k as usize] += 1;
        }
    }
    let rho = n as f64 / (la * lb).max(1.0);
    let mut cand: Vec<(f64, usize)> = Vec::with_capacity(n_win);
    for k in 0..n_win {
        let obs = cell[k] + cell.get(k + 1).copied().unwrap_or(0);
        if obs < 3 {
            continue;
        }
        let c = c_min + (k as f64 + 1.0) * tol;
        let e_box = rho * w * diag_len(la, lb, s, c);
        cand.push((tail_nlp(obs, e_box.max(1e-3), 1.0), k));
    }
    if cand.is_empty() {
        return (None, n_win);
    }
    cand.sort_by(|x, y| y.0.total_cmp(&x.0));
    let mut sb: Vec<f64> = anchors.iter().map(|a| s * a.b_pos as f64).collect();
    sb.sort_unstable_by(f64::total_cmp);
    let mut best: Option<Win> = None;
    for &(_, k) in cand.iter().take(refine) {
        let obs = cell[k] + cell.get(k + 1).copied().unwrap_or(0);
        let c = c_min + (k as f64 + 1.0) * tol;
        let shape = if m > 0.0 {
            anchor_shape(anchors, la, lb, s, c, tol, m)
        } else {
            "any"
        };
        if shape == "internal" || shape == "none" {
            continue;
        }
        let e = marginal_e(anchors, &sb, c, tol);
        let nlp = tail_nlp(obs, e.max(1e-3), 1.0);
        if best.is_none_or(|b| nlp > b.nlp) {
            best = Some(Win {
                nlp,
                c,
                obs,
                e,
                shape,
            });
        }
    }
    (best, n_win)
}

/// The overlap a window implies, in read coordinates, and its shape.
fn geometry(la: f64, lb: f64, s: f64, c: f64, tol: f64) -> (&'static str, f64) {
    let b_lo = (-c / s).max(0.0);
    let b_hi = ((la - c) / s).min(lb);
    let (a_lo, a_hi) = (s * b_lo + c, s * b_hi + c);
    let m = 2.0 * tol;
    let (a_start, a_end, b_start, b_end) = (a_lo <= m, a_hi >= la - m, b_lo <= m, b_hi >= lb - m);
    let shape = match (a_start, a_end, b_start, b_end) {
        (_, true, true, false) => "dove_ab", // a's suffix on b's prefix: a comes first
        (true, false, _, true) => "dove_ba",
        (_, _, true, true) => "b_in_a",
        (true, true, _, _) => "a_in_b",
        _ => "internal",
    };
    (shape, (b_hi - b_lo).max(0.0))
}

/// An overlap along `a = s*b + c`, clipped to the read box, one cell per `a` bin.
#[allow(clippy::too_many_arguments)]
fn overlap_for(
    a: u32,
    b: u32,
    la: f64,
    lb: f64,
    s: f64,
    c: f64,
    bin: u32,
    n: u32,
    obs: u32,
) -> Option<Overlap> {
    let b_lo = (-c / s).max(0.0);
    let b_hi = ((la - c) / s).min(lb);
    if b_hi <= b_lo {
        return None;
    }
    let (a_lo, a_hi) = (s * b_lo + c, s * b_hi + c);
    let bf = bin as f64;
    let (ab_lo, ab_hi) = (
        (a_lo / bf).floor() as u32,
        ((a_hi - 1.0).max(a_lo) / bf).floor() as u32,
    );
    let max_b = ((lb - 1.0) / bf).floor().max(0.0) as u32;
    let path: Vec<(u32, u32)> = (ab_lo..=ab_hi)
        .filter_map(|q| {
            // the bin centre, kept inside the overlap (the last bin's centre can lie past the read end)
            let bb = (((q as f64 + 0.5) * bf).clamp(a_lo, a_hi - 1.0) - c) / s / bf;
            (bb >= 0.0).then(|| (q, (bb.floor() as u32).min(max_b)))
        })
        .collect();
    let fold = |f: fn(&(u32, u32)) -> u32| {
        path.iter()
            .map(f)
            .fold((u32::MAX, 0u32), |(l, h), x| (l.min(x), h.max(x)))
    };
    let (asp, bsp) = (fold(|x| x.0), fold(|x| x.1));
    let aligned = (asp.1 - asp.0 + 1).min(bsp.1 - bsp.0 + 1);
    Some(Overlap {
        a,
        b,
        strand: Strand::Forward,
        n_anchors: n,
        matched: obs,
        cnt: obs,
        aligned_bins: aligned,
        bin_events: bin,
        gap: 0,
        slope: s,
        a_bin_span: asp,
        b_bin_span: bsp,
        path,
        z: f32::NAN,
    })
}

/// Detector settings shared by every mode: the seed cap and bin of the candidate set.
fn scan_params(args: &Cli) -> Result<Params> {
    let mut params = Params::r10_4_1();
    params.overlap.max_seed_occ = args.num("--max-seed-occ", 2000)?;
    params.overlap.bin_events = args.num("--bin", 192)?;
    Ok(params)
}

/// One read's `(hash, pos)` seeds, sorted.
type Seeds = Vec<(u64, u32)>;

/// The rescue's reads, in the candidate cache's order: every read's length in retained events, the
/// masked seeds of the reads sketched (empty for the others), and the pod5's signal content hashes of
/// the reads sketched, as `(read index, hash)`.
type CacheReads = (Vec<u32>, Vec<Seeds>, Vec<(usize, u64)>);

/// The masked, hash-sorted `(hash, pos)` seeds of every read and the reads' lengths, as the
/// detector sees them (hashes occurring more than the cap over the whole input dropped).
fn masked_seeds(pod5: &str, params: &Params) -> Result<(Vec<pipeline::ReadSketch>, Vec<Seeds>)> {
    let sketches =
        pipeline::sketch_pod5_streaming(std::path::Path::new(pod5), None, 4 << 30, params)?;
    let occ = SeedCounts::count(&sketches, |s| s.seeds.len(), |s, i| s.seeds[i].hash());
    let cap = params.overlap.max_seed_occ;
    let seeds: Vec<Seeds> = sketches
        .par_iter()
        .map(|s| {
            let mut v: Vec<(u64, u32)> = s
                .seeds
                .iter()
                .filter(|x| occ.get(x.hash()) <= cap)
                .map(|x| (x.hash(), x.pos))
                .collect();
            v.sort_unstable();
            v
        })
        .collect();
    eprintln!("{} reads sketched, cap {cap}", sketches.len());
    Ok((sketches, seeds))
}

/// [`masked_seeds`] + [`in_cache_order`] for the reads `want` alone (read indices in the candidate
/// cache's order, `names`), the mask from a seed-occurrence table (`occtable.rs`, haricot-pair-features
/// `--save-occ` over the same pod5) instead of every read's sketch: the pod5 is decoded once and
/// only the wanted reads are segmented and sketched. Returns every read's length from the cache
/// (`counts`, checked against the sketch for every wanted read), the masked, sorted seeds of
/// the wanted reads (empty for the others) and the wanted reads' signal content hashes
/// (`readorder::content_hash` of the same pA samples the sketch reads), by read index.
fn wanted_seeds(
    pod5: &str,
    params: &Params,
    names: &[String],
    counts: &[u32],
    want: &[u32],
    occ_p: &str,
) -> Result<CacheReads> {
    let occ = occtable::load(occ_p)?.counts;
    let cap = params.overlap.max_seed_occ;
    let at: FxHashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let want: rustc_hash::FxHashSet<usize> = want.iter().map(|&x| x as usize).collect();
    let mut seeds: Vec<Seeds> = vec![Vec::new(); names.len()];
    let mut content: Vec<(usize, u64)> = Vec::with_capacity(want.len());
    signal::for_each_adc_chunk(std::path::Path::new(pod5), None, 1 << 30, |chunk| {
        let got: Vec<(usize, u32, Seeds, u64)> = chunk
            .par_iter()
            .filter_map(|r| {
                let &i = at.get(r.id.as_str())?;
                if !want.contains(&i) {
                    return None;
                }
                // pipeline::sketch_reads_from's steps, for this read
                let pa = r.signal_pa();
                let events = signal::signal_to_events(&pa, params);
                let retained = emer::retained_quantized(&events, params);
                let sk = emer::sketch_from_retained(&retained, 0, Strand::Forward, params);
                let mut v: Seeds = sk
                    .iter()
                    .filter(|x| occ.get(x.hash()) <= cap)
                    .map(|x| (x.hash(), x.pos))
                    .collect();
                v.sort_unstable();
                let h = haricot::readorder::content_hash(&pa);
                Some((i, retained.len() as u32, v, h))
            })
            .collect();
        for (i, n, v, h) in got {
            anyhow::ensure!(
                n == counts[i],
                "{}: {n} events vs cache {}",
                names[i],
                counts[i]
            );
            seeds[i] = v;
            content.push((i, h));
        }
        Ok(())
    })?;
    let found = content.len();
    anyhow::ensure!(
        found == want.len(),
        "{found} of the {} wanted reads in the pod5",
        want.len()
    );
    eprintln!("{found} end reads segmented, cap {cap}");
    Ok((counts.to_vec(), seeds, content))
}

/// The sketches and seeds in the candidate cache's read order (`names`), matched by read id, so the
/// rescue reads any pod5 order of the cache's reads (a canonically numbered cache, or a pod5 written in
/// another order); the identity when the pod5 is in cache order. Errors unless both hold the same reads.
fn in_cache_order(
    sketches: Vec<pipeline::ReadSketch>,
    seeds: Vec<Seeds>,
    names: &[String],
) -> Result<(Vec<pipeline::ReadSketch>, Vec<Seeds>)> {
    if names.len() == sketches.len() && names.iter().zip(&sketches).all(|(n, s)| *n == s.name) {
        return Ok((sketches, seeds));
    }
    anyhow::ensure!(
        names.len() == sketches.len(),
        "the candidate cache holds {} reads, the pod5 {}",
        names.len(),
        sketches.len()
    );
    let perm: Vec<usize> = {
        let at: FxHashMap<&str, usize> = sketches
            .iter()
            .enumerate()
            .map(|(i, s)| (s.name.as_str(), i))
            .collect();
        anyhow::ensure!(at.len() == sketches.len(), "the pod5 repeats a read id");
        names
            .iter()
            .map(|n| {
                at.get(n.as_str())
                    .copied()
                    .with_context(|| format!("read {n} of the candidate cache is not in the pod5"))
            })
            .collect::<Result<_>>()?
    };
    let mut sk: Vec<Option<pipeline::ReadSketch>> = sketches.into_iter().map(Some).collect();
    let mut se: Vec<Option<Seeds>> = seeds.into_iter().map(Some).collect();
    let sketches = perm.iter().map(|&i| sk[i].take().unwrap()).collect();
    let seeds = perm.iter().map(|&i| se[i].take().unwrap()).collect();
    eprintln!("pod5 reads put in the candidate cache's order (by read id)");
    Ok((sketches, seeds))
}

/// Read `i`'s anchors (`seeds`, its masked seeds) against every read after it in `index` (hash ->
/// `(local read index, pos)`, ascending in the index): `(partner, anchor)`, grouped by partner.
/// Per partner, the multiset [`pair_anchors`] gives.
fn anchors_after(
    index: &FxHashMap<u64, Vec<(u32, u32)>>,
    seeds: &[(u64, u32)],
    i: u32,
) -> Vec<(u32, Anchor)> {
    let mut hits: Vec<(u32, Anchor)> = Vec::new();
    for &(h, pa) in seeds {
        if let Some(list) = index.get(&h) {
            let start = list.partition_point(|e| e.0 <= i);
            hits.extend(list[start..].iter().map(|&(lj, pb)| {
                (
                    lj,
                    Anchor {
                        a_pos: pa,
                        b_pos: pb,
                    },
                )
            }));
        }
    }
    hits.sort_unstable_by_key(|h| h.0);
    hits
}

/// Every anchor of the pair `(a, b)`: a merge-join of the two reads' masked seeds.
#[cfg(test)]
fn pair_anchors(va: &[(u64, u32)], vb: &[(u64, u32)]) -> Vec<Anchor> {
    let mut anchors: Vec<Anchor> = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < va.len() && j < vb.len() {
        match va[i].0.cmp(&vb[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                let h = va[i].0;
                let (i0, j0) = (i, j);
                while i < va.len() && va[i].0 == h {
                    i += 1;
                }
                while j < vb.len() && vb[j].0 == h {
                    j += 1;
                }
                for x in &va[i0..i] {
                    for y in &vb[j0..j] {
                        anchors.push(Anchor {
                            a_pos: x.1,
                            b_pos: y.1,
                        });
                    }
                }
            }
        }
    }
    anchors
}

/// The per-offset scan settings.
struct ScanCfg {
    slopes: Vec<f64>,
    refine: usize,
    tol: f64,
    margin: f64,
}

/// `--slopes` (default `0.93,1,1.07`): the slope grid, comma-separated.
fn slopes(args: &Cli) -> Result<Vec<f64>> {
    let v = args.opt("--slopes").unwrap_or_else(|| "0.93,1,1.07".into());
    v.split(',')
        .map(|x| x.parse().with_context(|| format!("--slopes {v}: {x}")))
        .collect()
}

impl ScanCfg {
    fn from_args(args: &Cli, params: &Params) -> Result<Self> {
        let tol = params.overlap.bin_events as f64;
        Ok(ScanCfg {
            slopes: slopes(args)?,
            refine: args.num("--refine", 10)?,
            tol,
            // how close a window's anchors must come to a read end to count as reaching it
            // (0 = no check)
            margin: args.num("--end-margin", 2.0 * tol)?,
        })
    }

    /// `(anchors, best (window, slope), windows scanned)` of one pair.
    fn best_of(&self, anchors: &[Anchor], la: f64, lb: f64) -> (Option<(Win, f64)>, usize) {
        let mut best: Option<(Win, f64)> = None;
        let mut m = 0usize;
        for &s in &self.slopes {
            let (w, nw) = scan(anchors, la, lb, s, self.tol, self.refine, self.margin);
            m += nw;
            if let Some(w) = w
                && best.is_none_or(|(b, _)| w.nlp > b.nlp)
            {
                best = Some((w, s));
            }
        }
        (best, m)
    }
}

/// The score column of a pair's row as the selection reads it back: the pair's best window
/// `-log10 p` minus `log10` of the windows scanned, printed to 3 decimals and parsed again.
/// The selection compares THIS value (not the unrounded one) to its floor, as the score table
/// (`pool_detail.tsv`, `--out scores.tsv`) prints it.
fn score_as_read(w: &Win, m: usize) -> f64 {
    format!("{:.3}", w.nlp - (m.max(1) as f64).log10())
        .parse()
        .unwrap()
}

/// The shape column of a pair's row (see [`anchor_shape`] / [`geometry`]).
fn row_shape(w: &Win, la: f64, lb: f64, s: f64, tol: f64) -> &'static str {
    if w.shape == "any" {
        geometry(la, lb, s, w.c, tol).0
    } else {
        w.shape
    }
}

/// The shape of an implied overlap seen from the other read (`a_in_b` <-> `b_in_a`, `dove_ab` <->
/// `dove_ba`).
fn mirror_shape(shape: &'static str) -> &'static str {
    match shape {
        "a_in_b" => "b_in_a",
        "b_in_a" => "a_in_b",
        "dove_ab" => "dove_ba",
        "dove_ba" => "dove_ab",
        x => x,
    }
}

/// One pair's best window as its row reads it: `(score, shape, window, slope)`, the shape in the
/// owner's orientation (`a`, `b` of lengths `la`, `lb`; `anchors` as `a` gathered them, `a_pos` on
/// `a`). `flip` (`--symmetric`, when the owner is not the query of the canonical pair frame,
/// `readorder::query_first`): the scan runs with the two reads swapped and the shape is mirrored
/// back, so the score and the shape do not depend on which read owns the pair.
fn scan_row(
    cfg: &ScanCfg,
    anchors: &[Anchor],
    la: f64,
    lb: f64,
    flip: bool,
) -> Option<(f64, &'static str, Win, f64)> {
    let (best, m) = if flip {
        let swapped: Vec<Anchor> = anchors
            .iter()
            .map(|x| Anchor {
                a_pos: x.b_pos,
                b_pos: x.a_pos,
            })
            .collect();
        cfg.best_of(&swapped, lb, la)
    } else {
        cfg.best_of(anchors, la, lb)
    };
    let (w, s) = best?;
    let sc = score_as_read(&w, m);
    let shape = if flip {
        mirror_shape(row_shape(&w, lb, la, s, cfg.tol))
    } else {
        row_shape(&w, la, lb, s, cfg.tol)
    };
    Some((sc, shape, w, s))
}

/// The overlap of a row's window ([`overlap_for`]) in the owner's orientation; with `flip` (see
/// [`scan_row`]) drawn in the scan's swapped frame and transposed back.
#[allow(clippy::too_many_arguments)]
fn row_overlap(
    a: u32,
    b: u32,
    la: f64,
    lb: f64,
    w: &Win,
    s: f64,
    bin: u32,
    n: u32,
    flip: bool,
) -> Option<Overlap> {
    if flip {
        overlap_for(b, a, lb, la, s, w.c, bin, n, w.obs)
            .map(|o| haricot::readorder::transpose_overlap(&o))
    } else {
        overlap_for(a, b, la, lb, s, w.c, bin, n, w.obs)
    }
}

/// A contig end: `(contig, is_right)`.
type End = (u32, bool);

/// The contig ends of a build and the reads at each: the first / last `end_nodes` nodes of every
/// contig of >= 2 mapped steps (`unitig_paths.tsv`, nodes by `(anchor_read, ev_beg)` from
/// `fbg_nodes.tsv`), and the reads threaded through them (`T` rows of the threads dump). Returns
/// every read's ends, keyed by read index.
fn read_ends(build: &str, threads: &str, end_nodes: usize) -> Result<FxHashMap<u32, Vec<End>>> {
    let rd = |p: String| std::fs::read_to_string(&p).with_context(|| format!("reading {p}"));
    let nodes = rd(format!("{build}/fbg_nodes.tsv"))?;
    let mut key2node: FxHashMap<(String, u32), u32> = FxHashMap::default();
    for l in nodes.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        key2node.insert((f[1].to_string(), f[2].parse()?), f[0][1..].parse()?);
    }
    // unitig_paths.tsv: contig (`u<N>`), step, anchor_read, ev_beg, ... -> per contig, its nodes by
    // step
    let up = rd(format!("{build}/unitig_paths.tsv"))?;
    let mut steps: Vec<(u32, u32, Option<u32>)> = Vec::new();
    for l in up.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        let node = key2node.get(&(f[2].to_string(), f[3].parse()?)).copied();
        // contig ids are `u<N>`; only their identity matters
        let contig = f[0].trim_start_matches(|c: char| !c.is_ascii_digit());
        steps.push((contig.parse()?, f[1].parse()?, node));
    }
    steps.sort_by_key(|s| (s.0, s.1));
    let mut node_end: FxHashMap<u32, Vec<End>> = FxHashMap::default();
    let mut k = 0;
    while k < steps.len() {
        let u = steps[k].0;
        let mut e = k;
        while e < steps.len() && steps[e].0 == u {
            e += 1;
        }
        let ns: Vec<u32> = steps[k..e].iter().filter_map(|s| s.2).collect();
        k = e;
        if ns.len() < 2 {
            continue;
        }
        let m = end_nodes.min(ns.len());
        for &x in &ns[..m] {
            node_end.entry(x).or_default().push((u, false));
        }
        for &x in &ns[ns.len() - m..] {
            node_end.entry(x).or_default().push((u, true));
        }
    }
    let mut ends: FxHashMap<u32, Vec<End>> = FxHashMap::default();
    for l in rd(threads.to_string())?.lines() {
        let Some(rest) = l.strip_prefix("T\t") else {
            continue;
        };
        let f: Vec<&str> = rest.split('\t').collect();
        let (rid, node): (u32, u32) = (f[0].parse()?, f[1].parse()?);
        if let Some(es) = node_end.get(&node) {
            ends.entry(rid).or_default().extend(es.iter().copied());
        }
    }
    for es in ends.values_mut() {
        es.sort_unstable();
        es.dedup();
    }
    Ok(ends)
}

/// Whether an implied overlap of `shape` fits the end pair `(ea, eb)` (ea holds the pair's lower
/// read): a at a right end and b at a left end needs a dovetail a->b, the reverse b->a, and a
/// containment fits any unless `dovetails_only`.
fn fits(ea: End, eb: End, shape: &str, dovetails_only: bool) -> bool {
    match shape {
        "a_in_b" | "b_in_a" => !dovetails_only,
        "dove_ab" => ea.1 && !eb.1,
        "dove_ba" => !ea.1 && eb.1,
        _ => false,
    }
}

/// `haricot-rescue rescue`: the contig-end rescue pool of a build, in one pass over the pod5.
///
/// Every pair of reads at ends of two DIFFERENT contigs that the candidate cache does not already
/// hold is scored by the per-offset scan; a pair joins the POOL when its score is at least
/// `log10(pairs tested) + --t` (a Bonferroni floor over this input's own universe) and its implied
/// overlap fits one of its end pairs. The universe of pairs is scored in memory and never written
/// out (on A. thaliana chr1 30x it is 32.6M pairs).
///
/// Writes to `--out-dir`: `pool.tsv` (read indices), `pool_detail.tsv` (score, shape and fitting
/// end pairs of every pool pair, for `admit`), `pool_ovl.bin` (the pool's overlaps) and
/// `cand_pool.bin` (the candidate cache plus the pool, for the stage-1 features).
///
/// Every pair is scanned in its canonical pair frame (the longer read as the query,
/// `readorder::query_first`; `scan_row`), its shape and overlap mirrored back to the owner's
/// orientation, so the pool does not depend on which read owns a pair: the default (`--symmetric`
/// is accepted); `--no-symmetric` scans every pair in the owner's frame.
///
/// `--occ OCC.bin`: the repeat mask from a seed-occurrence table (haricot-pair-features `--save-occ` over the
/// same pod5; `occtable.rs`) and only the end reads segmented, instead of sketching every read: the
/// same pool (the scan reads the end reads' seeds and lengths alone).
///
/// Usage: haricot-rescue rescue <pod5> <build_dir> <threads.tsv> <cand.bin> --out-dir D [--end-nodes 3]
///        [--t 4] [--dovetails-only] [--max-seed-occ 2000] [--bin 192] [--no-symmetric] [--occ OCC.bin]
fn rescue(argv: &[String]) -> Result<()> {
    let usage = "usage: haricot-rescue rescue <pod5> <build_dir> <threads.tsv> <cand.bin> --out-dir D \
                 [--end-nodes 3] [--t 4] [--dovetails-only] [--max-seed-occ 2000] [--bin 192] \
                 [--no-symmetric] [--occ OCC.bin]";
    let args = Cli::parse(argv, RESCUE_VALUED, RESCUE_SWITCHES, usage)?;
    args.positional(4, usage)?;
    let (pod5, build, threads, cand_p) = (&args.pos[0], &args.pos[1], &args.pos[2], &args.pos[3]);
    let out = args.opt("--out-dir").context("--out-dir is required")?;
    std::fs::create_dir_all(&out)?;
    let end_nodes: usize = args.num("--end-nodes", 3)?;
    let t: f64 = args.num("--t", 4.0)?;
    let dovetails_only = args.has("--dovetails-only");
    let params = scan_params(&args)?;
    let cfg = ScanCfg::from_args(&args, &params)?;
    // the canonical pair frame unless --no-symmetric (`scan_row`)
    let symmetric = args.symmetric()?;

    let t0 = std::time::Instant::now();
    let ends = read_ends(build, threads, end_nodes)?;
    let cache::Cache {
        counts,
        names,
        ovl: cand,
        content,
    } = cache::read(cand_p)?;
    let have: rustc_hash::FxHashSet<(u32, u32)> = cand.iter().map(|o| (o.a, o.b)).collect();

    // The universe: read pairs x < y with an end of each on different contigs, not candidates.
    let mut at: Vec<u32> = ends.keys().copied().collect();
    at.sort_unstable();
    let universe: Vec<(u32, u32)> = at
        .par_iter()
        .enumerate()
        .flat_map_iter(|(i, &x)| {
            let (ex, have, ends, at) = (&ends[&x], &have, &ends, &at);
            at[i + 1..].iter().filter_map(move |&y| {
                let ey = &ends[&y];
                let joined = ex.iter().any(|a| ey.iter().any(|b| a.0 != b.0));
                (joined && !have.contains(&(x, y))).then_some((x, y))
            })
        })
        .collect();
    let floor = (universe.len().max(1) as f64).log10() + t;
    eprintln!(
        "{} reads at {} contig ends; {} end-read pairs not among the {} candidates; floor {floor:.2}",
        at.len(),
        {
            let mut e: Vec<End> = ends.values().flatten().copied().collect();
            e.sort_unstable();
            e.dedup();
            e.len()
        },
        universe.len(),
        have.len()
    );

    let t_universe = t0.elapsed().as_secs_f64();
    // Every read's length and the masked seeds the scan reads (the end reads'): with `--occ`, the
    // mask from the occurrence table and only the end reads segmented; otherwise every read sketched.
    // `pod_content`: the pod5's signal content hashes of the reads sketched here (every read, or the
    // end reads), `(read index in the cache, hash)`, for the check against the cache's below.
    let (ev, seeds, pod_content): CacheReads = match args.opt("--occ") {
        Some(p) => wanted_seeds(pod5, &params, &names, &counts, &at, &p)?,
        None => {
            let (sketches, seeds) = masked_seeds(pod5, &params)?;
            let (sketches, seeds) = in_cache_order(sketches, seeds, &names)?;
            let ev: Vec<u32> = sketches.iter().map(|s| s.n_coord_events).collect();
            let pc: Vec<(usize, u64)> = sketches.iter().map(|s| s.content).enumerate().collect();
            (ev, seeds, pc)
        }
    };
    let t_seeds = t0.elapsed().as_secs_f64();
    let lens: Vec<f64> = ev.iter().map(|&e| e as f64).collect();
    let bin = params.overlap.bin_events;
    // the tie rank: the reads' signal content-hash order, the cache's (checked against the pod5's
    // for every read sketched here)
    let rank = if symmetric {
        if let Some(c) = &content {
            anyhow::ensure!(
                pod_content.iter().all(|&(i, h)| c[i] == h),
                "the pod5's signals are not the candidate cache's (content hashes differ)"
            );
        }
        cache::tie_rank(content.as_deref())?
    } else {
        Vec::new()
    };
    let flip_of = |a: u32, b: u32| {
        let (a, b) = (a as usize, b as usize);
        symmetric && !haricot::readorder::query_first(ev[a], ev[b], rank[a], rank[b])
    };

    // Score the universe; keep the pairs over the floor whose overlap fits an end pair.
    type PoolRow = (
        (u32, u32),
        f64,
        &'static str,
        Vec<(End, End)>,
        Option<Overlap>,
    );
    let score = |a: u32, b: u32, anchors: &[Anchor]| -> Option<PoolRow> {
        let (la, lb) = (lens[a as usize], lens[b as usize]);
        let flip = flip_of(a, b);
        let (sc, shape, w, s) = scan_row(&cfg, anchors, la, lb, flip)?;
        // A NaN score never joins the pool.
        if sc.is_nan() || sc < floor {
            return None;
        }
        let mut eps: Vec<(End, End)> = Vec::new();
        for &ea in &ends[&a] {
            for &eb in &ends[&b] {
                if ea.0 != eb.0 && fits(ea, eb, shape, dovetails_only) {
                    eps.push((ea, eb));
                }
            }
        }
        if eps.is_empty() {
            return None;
        }
        let o = row_overlap(a, b, la, lb, &w, s, bin, anchors.len() as u32, flip);
        Some(((a, b), sc, shape, eps, o))
    };
    // A pair's anchors are every (x, y) of the two reads' masked seeds with one hash, which
    // `pair_anchors` merge-joins per pair. Gathered instead per end read against an index of the
    // end reads' seeds: the same multiset per pair (and every step of the scoring is
    // order-free), but the cost follows the anchors rather than pairs x seeds -- the merge was
    // ~26k steps for each of 32.6M pairs on A. thaliana chr1 30x.
    let mut index: FxHashMap<u64, Vec<(u32, u32)>> = FxHashMap::default();
    for (li, &x) in at.iter().enumerate() {
        for &(h, p) in &seeds[x as usize] {
            index.entry(h).or_default().push((li as u32, p));
        }
    }
    let mut by_x: Vec<Vec<u32>> = vec![Vec::new(); at.len()];
    let pos: FxHashMap<u32, u32> = at.iter().enumerate().map(|(i, &x)| (x, i as u32)).collect();
    for &(x, y) in &universe {
        by_x[pos[&x] as usize].push(y);
    }
    let t_index = t0.elapsed().as_secs_f64();
    let pool: Vec<PoolRow> = (0..at.len())
        .into_par_iter()
        .flat_map_iter(|i| {
            let (x, ys) = (at[i], &by_x[i]);
            let mut rows = Vec::new();
            if ys.is_empty() {
                return rows.into_iter();
            }
            let hits = anchors_after(&index, &seeds[x as usize], i as u32);
            let mut anchors: Vec<Anchor> = Vec::new();
            for &y in ys {
                let ly = pos[&y];
                let lo = hits.partition_point(|h| h.0 < ly);
                let hi = hits.partition_point(|h| h.0 <= ly);
                anchors.clear();
                anchors.extend(hits[lo..hi].iter().map(|h| h.1));
                rows.extend(score(x, y, &anchors));
            }
            rows.into_iter()
        })
        .collect();
    eprintln!("pool: {} pairs fit their ends above the floor", pool.len());
    let t_score = t0.elapsed().as_secs_f64();
    eprintln!(
        "PHASE ends+universe {t_universe:.1} s, sketch+mask {:.1} s, end index {:.1} s, scoring {:.1} s",
        t_seeds - t_universe,
        t_index - t_seeds,
        t_score - t_index
    );

    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{out}/pool.tsv"))?);
    for ((a, b), ..) in &pool {
        writeln!(w, "{a}\t{b}")?;
    }
    w.flush()?;
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{out}/pool_detail.tsv"))?);
    writeln!(w, "#universe\t{}\tfloor\t{floor}", universe.len())?;
    writeln!(w, "a\tb\tscore\tshape\tends")?;
    for ((a, b), sc, shape, eps, _) in &pool {
        let es: Vec<String> = eps
            .iter()
            .map(|(x, y)| {
                let side = |e: &End| if e.1 { 'R' } else { 'L' };
                format!("{}{}-{}{}", x.0, side(x), y.0, side(y))
            })
            .collect();
        writeln!(w, "{a}\t{b}\t{sc}\t{shape}\t{}", es.join(";"))?;
    }
    w.flush()?;
    let pool_ovl: Vec<Overlap> = pool.iter().filter_map(|p| p.4.clone()).collect();
    std::fs::write(
        format!("{out}/pool_ovl.bin"),
        bincode::serialize(&pool_ovl)?,
    )?;
    let mut all = cand;
    all.extend(pool_ovl);
    all.sort_by_key(|o| (o.a, o.b));
    cache::write(
        &format!("{out}/cand_pool.bin"),
        &counts,
        &names,
        &all,
        content.as_deref(),
    )?;
    eprintln!("wrote {out}/pool.tsv, pool_detail.tsv, pool_ovl.bin, cand_pool.bin");
    Ok(())
}

/// `haricot-rescue admit`: the mutual best among the pool pairs the rescue forest allows
/// (`--allow`), added to a base cache, from `rescue`'s outputs and without the pod5 (the pool's
/// overlaps are the ones `--emit` would draw: the same function of the same anchors).
///
/// A pair is kept when, for one of its fitting end pairs, it is the best-scoring allowed pair at
/// BOTH ends (ties: the first in the order of the pair's two tie ranks, smaller first, so the pick
/// does not depend on the read numbering: the default, `--symmetric` accepted; with
/// `--no-symmetric`, the first in read-index order).
///
/// Usage: haricot-rescue admit <pool_detail.tsv> <pool_ovl.bin> <base.bin> --allow ACCEPTED.tsv
///        --out OUT.bin [--out-pairs ADMIT.tsv] [--no-symmetric]
fn admit(argv: &[String]) -> Result<()> {
    let usage = "usage: haricot-rescue admit <pool_detail.tsv> <pool_ovl.bin> <base.bin> --allow ACCEPTED.tsv \
                 --out OUT.bin [--out-pairs ADMIT.tsv] [--no-symmetric]";
    let args = Cli::parse(argv, ADMIT_VALUED, ADMIT_SWITCHES, usage)?;
    args.positional(3, usage)?;
    let (detail, pool_p, base_p) = (&args.pos[0], &args.pos[1], &args.pos[2]);
    let out = args.opt("--out").context("--out is required")?;
    let allow_p = args.opt("--allow").context("--allow is required")?;
    let cache::Cache {
        counts,
        names,
        mut ovl,
        content,
    } = cache::read(base_p)?;
    let symmetric = args.symmetric()?;
    let out_pairs = args.opt("--out-pairs");
    let line = admit::admit(
        detail,
        pool_p,
        &allow_p,
        &mut ovl,
        content.as_deref(),
        symmetric,
        out_pairs.as_deref(),
    )?;
    eprintln!("{line}");
    cache::write(&out, &counts, &names, &ovl, content.as_deref())?;
    eprintln!("wrote {out}");
    Ok(())
}

fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().collect();
    match argv.get(1).map(String::as_str) {
        Some("rescue") => return rescue(&argv[2..]),
        Some("admit") => return admit(&argv[2..]),
        _ => {}
    }
    let usage = "usage: haricot-rescue <pod5> <pairs.tsv> --out F [--emit BASE.bin] [--max-seed-occ 2000] \
                 [--bin 192] [--slopes 0.93,1,1.07] [--refine 10] [--end-margin M] | rescue ... | admit ...";
    let args = Cli::parse(&argv[1..], SCORE_VALUED, &[], usage)?;
    args.positional(2, usage)?;
    let (pod5, pairs_p) = (&args.pos[0], &args.pos[1]);
    let out_p = args.opt("--out").context("--out is required")?;
    let params = scan_params(&args)?;
    let slopes: Vec<f64> = slopes(&args)?;
    let refine: usize = args.num("--refine", 10)?;
    let tol = params.overlap.bin_events as f64;
    // how close a window's anchors must come to a read end to count as reaching it (0 = no check)
    let margin: f64 = args.num("--end-margin", 2.0 * tol)?;

    let text = std::fs::read_to_string(pairs_p).with_context(|| format!("reading {pairs_p}"))?;
    let pairs: Vec<(u32, u32)> = text
        .lines()
        .filter_map(|l| {
            let mut f = l.split('\t');
            match (f.next()?.parse().ok(), f.next()?.parse().ok()) {
                (Some(a), Some(b)) => Some((a, b)),
                _ => None,
            }
        })
        .collect();
    eprintln!("{} pairs; tol {tol} events; slopes {slopes:?}", pairs.len());

    let sketches =
        pipeline::sketch_pod5_streaming(std::path::Path::new(pod5), None, 4 << 30, &params)?;
    let occ = SeedCounts::count(&sketches, |s| s.seeds.len(), |s, i| s.seeds[i].hash());
    let cap = params.overlap.max_seed_occ;
    let seeds: Vec<Vec<(u64, u32)>> = sketches
        .iter()
        .map(|s| {
            let mut v: Vec<(u64, u32)> = s
                .seeds
                .iter()
                .filter(|x| occ.get(x.hash()) <= cap)
                .map(|x| (x.hash(), x.pos))
                .collect();
            v.sort_unstable();
            v
        })
        .collect();
    let lens: Vec<f64> = sketches.iter().map(|s| s.n_coord_events as f64).collect();
    eprintln!("{} reads sketched, cap {cap}", sketches.len());

    let best_of = |a: u32, b: u32| -> (usize, Option<(Win, f64)>, usize) {
        let (va, vb) = (&seeds[a as usize], &seeds[b as usize]);
        let mut anchors: Vec<Anchor> = Vec::new();
        let (mut i, mut j) = (0usize, 0usize);
        while i < va.len() && j < vb.len() {
            match va[i].0.cmp(&vb[j].0) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    let h = va[i].0;
                    let (i0, j0) = (i, j);
                    while i < va.len() && va[i].0 == h {
                        i += 1;
                    }
                    while j < vb.len() && vb[j].0 == h {
                        j += 1;
                    }
                    for x in &va[i0..i] {
                        for y in &vb[j0..j] {
                            anchors.push(Anchor {
                                a_pos: x.1,
                                b_pos: y.1,
                            });
                        }
                    }
                }
            }
        }
        let (la, lb) = (lens[a as usize], lens[b as usize]);
        let mut best: Option<(Win, f64)> = None;
        let mut m = 0usize;
        for &s in &slopes {
            let (w, nw) = scan(&anchors, la, lb, s, tol, refine, margin);
            m += nw;
            if let Some(w) = w
                && best.is_none_or(|(b, _)| w.nlp > b.nlp)
            {
                best = Some((w, s));
            }
        }
        (anchors.len(), best, m)
    };

    if let Some(base_p) = args.opt("--emit") {
        let cache::Cache {
            counts,
            names,
            mut ovl,
            content,
        } = cache::read(&base_p)?;
        anyhow::ensure!(
            names.len() == sketches.len() && names.iter().zip(&sketches).all(|(n, s)| *n == s.name),
            "the base cache's reads are not the pod5's reads in the same order"
        );
        let bin = params.overlap.bin_events;
        let add: Vec<Overlap> = pairs
            .par_iter()
            .filter_map(|&(a, b)| {
                let (n, best, _) = best_of(a, b);
                let (w, s) = best?;
                overlap_for(
                    a,
                    b,
                    lens[a as usize],
                    lens[b as usize],
                    s,
                    w.c,
                    bin,
                    n as u32,
                    w.obs,
                )
            })
            .collect();
        eprintln!("emit: {} base overlaps + {} admitted", ovl.len(), add.len());
        ovl.extend(add);
        ovl.sort_by_key(|o| (o.a, o.b));
        cache::write(&out_p, &counts, &names, &ovl, content.as_deref())?;
        eprintln!("wrote {out_p}");
        return Ok(());
    }

    let rows: Vec<String> = pairs
        .par_iter()
        .map(|&(a, b)| {
            let (la, lb) = (lens[a as usize], lens[b as usize]);
            let (n_anchors, best, m) = best_of(a, b);
            match best {
                Some((w, s)) => {
                    let (box_shape, ov) = geometry(la, lb, s, w.c, tol);
                    let shape = if w.shape == "any" { box_shape } else { w.shape };
                    format!(
                        "{a}\t{b}\t{}\t{la}\t{lb}\t{:.3}\t{:.3}\t{s}\t{:.0}\t{}\t{:.3}\t{shape}\t{ov:.0}",
                        n_anchors,
                        w.nlp - (m.max(1) as f64).log10(),
                        w.nlp,
                        w.c,
                        w.obs,
                        w.e
                    )
                }
                None => format!("{a}\t{b}\t{n_anchors}\t{la}\t{lb}\tnan\tnan\tnan\tnan\t0\tnan\tnone\t0"),
            }
        })
        .collect();
    let mut o = std::io::BufWriter::new(std::fs::File::create(&out_p)?);
    writeln!(
        o,
        "a\tb\tn\tla\tlb\tscore\tnlp_win\tslope\toffset\tobs\te\tshape\tov_events"
    )?;
    for r in rows {
        writeln!(o, "{r}")?;
    }
    eprintln!("wrote {out_p}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The calls of `haricot assemble` and of a replay parse, with their frames.
    #[test]
    fn the_calls_in_the_tree_parse() {
        let v = |l: &str| l.split_whitespace().map(String::from).collect::<Vec<_>>();
        for (l, sym) in [
            (
                "P eb thr.tsv C --out-dir q --dovetails-only --occ occ.bin",
                true,
            ),
            (
                "P eb thr.tsv C --out-dir q --t 4 --dovetails-only --symmetric --max-seed-occ 2000 \
                 --bin 192 --occ occ.bin",
                true,
            ),
            (
                "P eb thr.tsv C --out-dir q --t 4 --dovetails-only --max-seed-occ 2000 --bin 192 \
                 --no-symmetric",
                false,
            ),
        ] {
            let a = Cli::parse(&v(l), RESCUE_VALUED, RESCUE_SWITCHES, "").unwrap();
            a.positional(4, "").unwrap();
            assert_eq!(a.symmetric().unwrap(), sym, "{l}");
            assert!(a.has("--dovetails-only"));
            let p = scan_params(&a).unwrap();
            assert_eq!((p.overlap.max_seed_occ, p.overlap.bin_events), (2000, 192));
            assert_eq!(a.num("--t", 4.0).unwrap(), 4.0);
            let c = ScanCfg::from_args(&a, &p).unwrap();
            assert_eq!(
                (c.slopes, c.refine, c.margin),
                (vec![0.93, 1.0, 1.07], 10, 384.0)
            );
        }
        let l = "D O K --allow A --out q.bin --out-pairs a.tsv --no-symmetric";
        let a = Cli::parse(&v(l), ADMIT_VALUED, ADMIT_SWITCHES, "").unwrap();
        assert!(a.positional(3, "").is_ok() && !a.symmetric().unwrap());
        assert!(
            Cli::parse(
                &v("P eb t C --out-dir q --dovetail-only"),
                RESCUE_VALUED,
                RESCUE_SWITCHES,
                ""
            )
            .is_err()
        );
        let bad = Cli::parse(
            &v("P eb t C --slopes 1,x"),
            RESCUE_VALUED,
            RESCUE_SWITCHES,
            "",
        )
        .unwrap();
        assert!(slopes(&bad).is_err());
    }

    /// The per-read gather gives every pair the merge-join's anchors, as a multiset.
    #[test]
    fn gathered_anchors_are_the_merge_joins() {
        let mut h: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut rnd = move |m: u64| {
            h ^= h << 13;
            h ^= h >> 7;
            h ^= h << 17;
            h % m
        };
        let reads: Vec<Vec<(u64, u32)>> = (0..12)
            .map(|_| {
                let mut v: Vec<(u64, u32)> =
                    (0..300).map(|_| (rnd(40), rnd(5000) as u32)).collect();
                v.sort_unstable();
                v
            })
            .collect();
        let mut index: FxHashMap<u64, Vec<(u32, u32)>> = FxHashMap::default();
        for (li, r) in reads.iter().enumerate() {
            for &(hh, p) in r {
                index.entry(hh).or_default().push((li as u32, p));
            }
        }
        let key = |a: &Anchor| (a.a_pos, a.b_pos);
        for i in 0..reads.len() {
            let hits = anchors_after(&index, &reads[i], i as u32);
            for j in i + 1..reads.len() {
                let mut got: Vec<_> = hits
                    .iter()
                    .filter(|x| x.0 == j as u32)
                    .map(|x| key(&x.1))
                    .collect();
                let mut want: Vec<_> = pair_anchors(&reads[i], &reads[j]).iter().map(key).collect();
                got.sort_unstable();
                want.sort_unstable();
                assert_eq!(got, want, "pair {i} {j}");
            }
        }
    }

    /// A 60-anchor dovetail in the corner of a box that also holds 3,000 uniform random anchors:
    /// the per-offset scan finds it at the corner, as a dovetail a->b of the planted length.
    #[test]
    fn corner_dovetail_stands_out_of_a_large_random_cloud() {
        let (la, lb) = (20_000.0, 24_000.0);
        let mut h: u64 = 7;
        let mut rnd = || {
            h = h
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (h >> 33) as f64 / (1u64 << 31) as f64
        };
        let mut anchors: Vec<Anchor> = (0..3000)
            .map(|_| Anchor {
                a_pos: (rnd() * la) as u32,
                b_pos: (rnd() * lb) as u32,
            })
            .collect();
        // overlap of 1,500 events: a in [18500, 20000] <-> b in [0, 1500]
        for k in 0..60 {
            let b = k as f64 * 25.0;
            anchors.push(Anchor {
                a_pos: (18_500.0 + b + 30.0 * (rnd() - 0.5)) as u32,
                b_pos: b as u32,
            });
        }
        let (w, _) = scan(&anchors, la, lb, 1.0, 192.0, 5, 384.0);
        let w = w.unwrap();
        assert!(w.nlp > 20.0, "{w:?}");
        assert!((w.c - 18_500.0).abs() < 400.0, "{w:?}");
        assert_eq!(w.shape, "dove_ab");
        let (shape, ov) = geometry(la, lb, 1.0, w.c, 192.0);
        assert_eq!(shape, "dove_ab");
        assert!((ov - 1_500.0).abs() < 400.0, "{ov}");
    }

    /// A segment shared INSIDE both reads (a repeat copy) is dense but reaches no read end: with the
    /// end check it is not taken as the overlap; without it, it is.
    #[test]
    fn an_internal_shared_segment_is_not_an_overlap() {
        let (la, lb) = (20_000.0, 24_000.0);
        let anchors: Vec<Anchor> = (0..80)
            .map(|k| Anchor {
                a_pos: 8_000 + 25 * k,
                b_pos: 11_000 + 25 * k,
            })
            .collect();
        let (w, _) = scan(&anchors, la, lb, 1.0, 192.0, 5, 384.0);
        assert!(w.is_none(), "{w:?}");
        let (w, _) = scan(&anchors, la, lb, 1.0, 192.0, 5, 0.0);
        assert!(w.is_some());
    }

    /// The emitted overlap of a dovetail a->b runs from a's last bins to b's first ones, one cell
    /// per `a` bin, on the window's diagonal.
    #[test]
    fn emitted_dovetail_overlap_spans_the_corner() {
        let o = overlap_for(3, 9, 20_000.0, 24_000.0, 1.0, 18_500.0, 192, 3060, 60).unwrap();
        assert_eq!(o.a_bin_span, (96, 104)); // 18,500 / 192 = 96.4; 19,999 / 192 = 104.2
        assert_eq!(o.b_bin_span, (0, 7)); // b = a - 18,500 over [0, 1,500)
        assert_eq!(o.path.len(), 9);
        assert!(
            o.path
                .windows(2)
                .all(|w| w[1].0 == w[0].0 + 1 && w[1].1 >= w[0].1)
        );
        assert_eq!(o.aligned_bins, 8);
        assert!(o.z.is_nan());
        // a line that misses the box emits nothing
        assert!(overlap_for(3, 9, 20_000.0, 24_000.0, 1.0, 25_000.0, 192, 10, 5).is_none());
    }

    /// `--symmetric`: a pair gathered by either read gives the same score, the mirrored shape and the
    /// transposed overlap. Read y (20,000 events) ends on read x's (24,000) start; x is the query.
    #[test]
    fn the_rescue_scan_does_not_depend_on_the_owner() {
        let (ly, lx) = (20_000.0, 24_000.0);
        let mut h: u64 = 5;
        let mut rnd = || {
            h = h
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (h >> 33) as f64 / (1u64 << 31) as f64
        };
        // owner y: a_pos on y, b_pos on x
        let mut from_y: Vec<Anchor> = (0..2000)
            .map(|_| Anchor {
                a_pos: (rnd() * ly) as u32,
                b_pos: (rnd() * lx) as u32,
            })
            .collect();
        for k in 0..60 {
            let b = k as f64 * 25.0;
            from_y.push(Anchor {
                a_pos: (18_500.0 + b + 30.0 * (rnd() - 0.5)) as u32,
                b_pos: (b * 1.04) as u32,
            });
        }
        let from_x: Vec<Anchor> = from_y
            .iter()
            .map(|x| Anchor {
                a_pos: x.b_pos,
                b_pos: x.a_pos,
            })
            .collect();
        let mut params = Params::r10_4_1();
        params.overlap.bin_events = 192;
        let none = Cli::parse(&[], RESCUE_VALUED, RESCUE_SWITCHES, "").unwrap();
        let cfg = ScanCfg::from_args(&none, &params).unwrap();
        // numbering 1: y = 0 owns the pair and is not the query (flip); numbering 2: x = 0 owns it
        let (s1, sh1, w1, sl1) = scan_row(&cfg, &from_y, ly, lx, true).unwrap();
        let (s2, sh2, w2, sl2) = scan_row(&cfg, &from_x, lx, ly, false).unwrap();
        assert_eq!((s1.to_bits(), sl1.to_bits()), (s2.to_bits(), sl2.to_bits()));
        assert_eq!((w1.obs, w1.c.to_bits()), (w2.obs, w2.c.to_bits()));
        assert_eq!(
            (sh1, sh2),
            ("dove_ab", "dove_ba"),
            "y's suffix on x's prefix"
        );
        let o1 = row_overlap(0, 1, ly, lx, &w1, sl1, 192, 2060, true).unwrap();
        let o2 = row_overlap(0, 1, lx, ly, &w2, sl2, 192, 2060, false).unwrap();
        let t2 = haricot::readorder::transpose_overlap(&o2);
        assert_eq!(
            (o1.path, o1.a_bin_span, o1.b_bin_span),
            (t2.path, t2.a_bin_span, t2.b_bin_span)
        );
        // without the frame, the owner y scans in its own frame
        assert!(scan_row(&cfg, &from_y, ly, lx, false).is_some());
    }

    #[test]
    fn a_pure_random_cloud_scores_near_zero_after_the_window_correction() {
        let (la, lb) = (20_000.0, 24_000.0);
        let mut h: u64 = 11;
        let mut rnd = || {
            h = h
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (h >> 33) as f64 / (1u64 << 31) as f64
        };
        let anchors: Vec<Anchor> = (0..3000)
            .map(|_| Anchor {
                a_pos: (rnd() * la) as u32,
                b_pos: (rnd() * lb) as u32,
            })
            .collect();
        let (w, m) = scan(&anchors, la, lb, 1.0, 192.0, 5, 0.0);
        let score = w.unwrap().nlp - (m as f64).log10();
        assert!(score < 2.0, "{score}");
    }
}
