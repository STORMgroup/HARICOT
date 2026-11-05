//! The SIGNAL-EXACT oracle: perfect overlaps whose paths come from each event's own
//! reference position, not from a rate model.
//!
//! A rate model puts read base `e / n_events * rlen` under event `e` -- one uniform rate
//! per read. `truth/event_ref.py` places each event through the basecaller's move table (event
//! -> read base) and a CIGAR walk (read base -> reference), under two aligners, with every
//! event labelled CONF / SINGLE / LOWQ / CONFLICT / UNPLACED. Here a pair's path is built
//! bin by bin: the median reference position of read A's bin `i` is looked up among read
//! B's events, and the cell is `(i, e_B / bin)`. Local speed variation, adapter samples,
//! soft clips and indels between the two reads are therefore exact, and a bin with no
//! counterpart in B within `--tol` bp (a deletion in B, B's clipped end) leaves a hole
//! rather than an interpolated cell.
//!
//! Membership is the same-strand intersection of two placed SEGMENTS (an alignment
//! record's run of events), using only events whose status is allowed:
//!   `--status conf`     CONF only (both aligners agree, mapq >= 10)          [default]
//!   `--status lenient`  CONF + SINGLE (one aligner, mapq >= 10)
//! `--primary-only` drops supplementary segments, i.e. treats a chimera's second locus
//! as unplaced rather than as the true overlap it is in the signal.
//! `--cut-physical` drops only the supplementary segments that are MOLECULE junctions: a
//! fold-back (opposite strand to one of the read's primary segments of the same aligner,
//! overlapping it or within 2 kb) or an inter-locus chimera (> 1 Mb from every primary, or
//! on another contig). Every other split -- same-strand splits, inversions against the
//! reference -- is kept: on A. thaliana those are mostly real joins of the sample.
//!
//! `--truth-out` writes EVERY same-strand read pair whose placed segments (any status
//! but CONFLICT) intersect by >= `--truth-min` bp, with the oracle's verdict, so a
//! classifier can train on exact labels and treat the rest as UNKNOWN.
//!
//! Usage:
//!   haricot-truth-oracle IN_CACHE OUT_CACHE --evref P.evref.bin --segs P.segs.tsv
//!              [--min-ovl 2000] [--min-bins 8] [--status conf|lenient]
//!              [--primary-only | --cut-physical] [--tol 60] [--truth-out T.tsv] [--truth-min 500]
use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result, bail};
use haricot::overlap::Overlap;
use haricot::params::Strand;

const CONF: u8 = 1;
const SINGLE: u8 = 2;
const CONFLICT: u8 = 4;

fn flag<T: std::str::FromStr>(a: &[String], n: &str, d: T) -> T {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(d)
}
fn opt(a: &[String], n: &str) -> Option<String> {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .cloned()
}

struct EvRef {
    rf: Vec<i64>,
    seg: Vec<i32>,
    status: Vec<u8>,
}

fn read_evref(path: &str) -> Result<HashMap<String, EvRef>> {
    let buf = std::fs::read(path).with_context(|| format!("reading {path}"))?;
    let mut off = 0usize;
    let u32_at = |o: usize| u32::from_le_bytes(buf[o..o + 4].try_into().unwrap());
    let mut out = HashMap::new();
    while off < buf.len() {
        let nl = u32_at(off) as usize;
        off += 4;
        let name = String::from_utf8(buf[off..off + nl].to_vec())?;
        off += nl;
        let n = u32_at(off) as usize;
        off += 4;
        let rf = (0..n)
            .map(|k| i64::from_le_bytes(buf[off + 8 * k..off + 8 * k + 8].try_into().unwrap()))
            .collect();
        off += 8 * n;
        off += 4 * n; // read base: not needed here
        let seg = (0..n)
            .map(|k| i32::from_le_bytes(buf[off + 4 * k..off + 4 * k + 4].try_into().unwrap()))
            .collect();
        off += 4 * n;
        let status = buf[off..off + n].to_vec();
        off += n;
        out.insert(name, EvRef { rf, seg, status });
    }
    Ok(out)
}

/// One placed segment of one read, restricted to the events a run may use.
struct Seg {
    read: u32,
    minus: bool,
    primary: bool,
    lo: i64,
    hi: i64,
    /// `(ref, event)` sorted by ref, for the B-side lookup.
    by_ref: Vec<(i64, u32)>,
}

impl Seg {
    /// The event of this segment placed nearest to `r`, if within `tol`.
    fn lookup(&self, r: i64, tol: i64) -> Option<u32> {
        let k = self.by_ref.partition_point(|&(x, _)| x < r);
        let mut best: Option<(i64, u32)> = None;
        for j in [k.wrapping_sub(1), k] {
            if let Some(&(x, e)) = self.by_ref.get(j) {
                let d = (x - r).abs();
                if d <= tol && best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, e));
                }
            }
        }
        best.map(|(_, e)| e)
    }
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let usage = "usage: haricot-truth-oracle IN OUT --evref P.evref.bin --segs P.segs.tsv";
    let cin = args.get(1).context(usage)?;
    let cout = args.get(2).context(usage)?;
    let evp = opt(&args, "--evref").context("--evref is required")?;
    let segp = opt(&args, "--segs").context("--segs is required")?;
    let min_ovl: i64 = flag(&args, "--min-ovl", 2000);
    let min_bins: usize = flag(&args, "--min-bins", 8);
    let tol: i64 = flag(&args, "--tol", 60);
    let status_mode = opt(&args, "--status").unwrap_or_else(|| "conf".into());
    let primary_only = args.iter().any(|a| a == "--primary-only");
    let cut_physical = args.iter().any(|a| a == "--cut-physical");
    if primary_only && cut_physical {
        bail!("--primary-only and --cut-physical are alternatives");
    }
    let truth_min: i64 = flag(&args, "--truth-min", 500);
    let allowed = |s: u8| match status_mode.as_str() {
        "lenient" => s == CONF || s == SINGLE,
        _ => s == CONF,
    };
    if !matches!(status_mode.as_str(), "conf" | "lenient") {
        bail!("--status must be conf or lenient");
    }

    let bytes = std::fs::read(cin).with_context(|| format!("reading {cin}"))?;
    let (counts, names, _det): (Vec<u32>, Vec<String>, Vec<Overlap>) =
        bincode::deserialize(&bytes).context("decoding cache")?;
    let bin: u32 = flag(
        &args,
        "--bin",
        _det.first().map(|o| o.bin_events).unwrap_or(192),
    );
    drop(_det);

    let mut ev = read_evref(&evp)?;
    let (mut n_have, mut n_bad) = (0usize, 0usize);
    let mut per_read: Vec<Option<EvRef>> = Vec::with_capacity(names.len());
    for (i, nm) in names.iter().enumerate() {
        match ev.remove(nm) {
            Some(e) if e.rf.len() == counts[i] as usize => {
                n_have += 1;
                per_read.push(Some(e));
            }
            Some(e) => {
                n_bad += 1;
                if n_bad <= 5 {
                    eprintln!(
                        "  count mismatch {nm}: evref {} vs cache {}",
                        e.rf.len(),
                        counts[i]
                    );
                }
                per_read.push(None);
            }
            None => per_read.push(None),
        }
    }
    eprintln!(
        "evref: {n_have} of {} cache reads matched event-for-event, {n_bad} count mismatches",
        names.len()
    );
    if n_bad > 0 {
        bail!("event counts disagree with the cache: the span extraction is not the pipeline's");
    }

    // Segment metadata: strand and primary flag per global segment id.
    let mut meta: HashMap<i32, (bool, bool)> = HashMap::new();
    // --cut-physical: (read, aligner) -> its primary segments' (minus, ref_lo, ref_hi) on the
    // global axis, and every supplementary segment, to find the molecule junctions.
    type Prim = HashMap<(String, String), Vec<(bool, i64, i64)>>; // (read, aligner) -> (minus, lo, hi)
    let mut prim: Prim = HashMap::new();
    let mut supp: Vec<(i32, String, String, bool, i64, i64)> = Vec::new();
    for line in std::fs::read_to_string(&segp)?.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        let (sid, minus, primary) = (f[0].parse()?, f[5] == "-", f[3] == "1");
        meta.insert(sid, (minus, primary));
        if cut_physical {
            let (lo, hi): (i64, i64) = (f[9].parse()?, f[10].parse()?);
            let key = (f[1].to_string(), f[2].to_string());
            if primary {
                prim.entry(key).or_default().push((minus, lo, hi));
            } else {
                supp.push((sid, key.0, key.1, minus, lo, hi));
            }
        }
    }
    const FOLD_GAP: i64 = 2_000;
    const CHIMERA_DIST: i64 = 1_000_000;
    let mut cut: HashSet<i32> = HashSet::new();
    if cut_physical {
        let (mut n_fold, mut n_chim, mut n_kept) = (0usize, 0usize, 0usize);
        for (sid, rd, al, minus, lo, hi) in &supp {
            let Some(ps) = prim.get(&(rd.clone(), al.clone())) else {
                n_kept += 1;
                continue;
            };
            // gap between two intervals; negative when they overlap
            let gap = |plo: i64, phi: i64| (*lo - phi).max(plo - *hi);
            if ps
                .iter()
                .any(|&(pm, plo, phi)| pm != *minus && gap(plo, phi) <= FOLD_GAP)
            {
                n_fold += 1;
                cut.insert(*sid);
            } else if ps.iter().all(|&(_, plo, phi)| gap(plo, phi) > CHIMERA_DIST) {
                n_chim += 1;
                cut.insert(*sid);
            } else {
                n_kept += 1;
            }
        }
        eprintln!(
            "cut-physical: {} supplementary segments: {n_fold} fold-back and {n_chim} inter-locus segments cut, {n_kept} kept",
            supp.len()
        );
    }

    // Build the usable segments: `strict` for the oracle, `any` for the truth table.
    let build = |ok: &dyn Fn(u8) -> bool| -> Vec<Seg> {
        let mut segs = Vec::new();
        for (r, e) in per_read.iter().enumerate() {
            let Some(e) = e else { continue };
            let mut by: HashMap<i32, Vec<(i64, u32)>> = HashMap::new();
            for k in 0..e.rf.len() {
                if e.seg[k] >= 0 && e.rf[k] >= 0 && ok(e.status[k]) {
                    by.entry(e.seg[k]).or_default().push((e.rf[k], k as u32));
                }
            }
            for (sid, mut v) in by {
                let (minus, primary) = meta.get(&sid).copied().unwrap_or((false, true));
                if (primary_only && !primary) || cut.contains(&sid) {
                    continue;
                }
                v.sort_unstable();
                segs.push(Seg {
                    read: r as u32,
                    minus,
                    primary,
                    lo: v[0].0,
                    hi: v[v.len() - 1].0 + 1,
                    by_ref: v,
                });
            }
        }
        segs.sort_by_key(|s| (s.minus, s.lo));
        segs
    };
    let segs = build(&|s| allowed(s));
    let any_segs = build(&|s| s != CONFLICT && s != 0);
    eprintln!(
        "segments: {} usable ({} supplementary), {} for the truth table",
        segs.len(),
        segs.iter().filter(|s| !s.primary).count(),
        any_segs.len()
    );

    // The path of one segment pair; `a` is the lower read id.
    let path_of = |sa: &Seg, sb: &Seg| -> Vec<(u32, u32)> {
        let ea = per_read[sa.read as usize].as_ref().unwrap();
        let n_a = counts[sa.read as usize];
        let a_max = (n_a / bin).saturating_sub(1);
        let b_max = (counts[sb.read as usize] / bin).saturating_sub(1);
        let (lo, hi) = (sa.lo.max(sb.lo), sa.hi.min(sb.hi));
        let mut path = Vec::new();
        let mut buf: Vec<i64> = Vec::with_capacity(bin as usize);
        // Only A bins whose events touch the shared interval.
        let first = sa.by_ref.iter().filter(|&&(x, _)| x >= lo && x < hi);
        let (mut bmin, mut bmax) = (u32::MAX, 0u32);
        for &(_, e) in first {
            bmin = bmin.min(e / bin);
            bmax = bmax.max(e / bin);
        }
        if bmin == u32::MAX {
            return path;
        }
        let sid = ea.seg[sa.by_ref[0].1 as usize];
        for i in bmin..=bmax.min(a_max) {
            buf.clear();
            let (e0, e1) = (i * bin, ((i + 1) * bin).min(n_a));
            for k in e0..e1 {
                let k = k as usize;
                if ea.seg[k] == sid && ea.rf[k] >= 0 && allowed(ea.status[k]) {
                    buf.push(ea.rf[k]);
                }
            }
            if buf.len() < (bin / 4) as usize {
                continue;
            }
            buf.sort_unstable();
            let r = buf[buf.len() / 2];
            if r < lo || r >= hi {
                continue;
            }
            if let Some(eb) = sb.lookup(r, tol) {
                let bb = eb / bin;
                if bb <= b_max {
                    path.push((i, bb));
                }
            }
        }
        path
    };

    let mut best: HashMap<(u32, u32), Overlap> = HashMap::new();
    let mut n_pairs = 0usize;
    let mut holes_tot = 0u64;
    let mut span_tot = 0u64;
    for (p, x) in segs.iter().enumerate() {
        for y in &segs[p + 1..] {
            if y.minus != x.minus || y.lo >= x.hi {
                break;
            }
            if y.read == x.read || x.hi.min(y.hi) - x.lo.max(y.lo) < min_ovl {
                continue;
            }
            n_pairs += 1;
            let (sa, sb) = if x.read < y.read { (x, y) } else { (y, x) };
            let path = path_of(sa, sb);
            if path.len() < min_bins {
                continue;
            }
            let (a0, a1) = (path[0].0, path[path.len() - 1].0);
            let blo = path.iter().map(|c| c.1).min().unwrap();
            let bhi = path.iter().map(|c| c.1).max().unwrap();
            let n = path.len() as u32;
            let holes = (a1 - a0 + 1) - n;
            let o = Overlap {
                a: sa.read,
                b: sb.read,
                strand: Strand::Forward,
                n_anchors: n,
                matched: n * bin,
                cnt: n,
                aligned_bins: n,
                bin_events: bin,
                gap: holes,
                slope: if bhi > blo {
                    (a1 - a0) as f64 / (bhi - blo) as f64
                } else {
                    1.0
                },
                a_bin_span: (a0, a1),
                b_bin_span: (blo, bhi),
                path,
                z: f32::NAN,
            };
            let key = (sa.read, sb.read);
            if best.get(&key).is_none_or(|b| b.path.len() < o.path.len()) {
                best.insert(key, o);
            }
        }
    }
    let mut oracle: Vec<Overlap> = best.into_values().collect();
    oracle.sort_by_key(|o| (o.a, o.b));
    for o in &oracle {
        holes_tot += o.gap as u64;
        span_tot += (o.a_bin_span.1 - o.a_bin_span.0 + 1) as u64;
    }
    eprintln!(
        "oracle: {} overlaps from {n_pairs} segment pairs >= {min_ovl} bp (status={status_mode}, \
         primary_only={primary_only}); holes {:.2}% of spanned bins",
        oracle.len(),
        100.0 * holes_tot as f64 / span_tot.max(1) as f64
    );

    if let Some(tp) = opt(&args, "--truth-out") {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(&tp)?);
        writeln!(
            w,
            "a\tb\tname_a\tname_b\tstrand\tov_bp\tsupp\tin_oracle\toracle_cells\t\
             a_lo\ta_hi\tb_lo\tb_hi"
        )?;
        // (a, b) -> (oracle index, a_bin_span, b_bin_span)
        type Spans = HashMap<(u32, u32), (usize, (u32, u32), (u32, u32))>;
        let inor: Spans = oracle
            .iter()
            .map(|o| ((o.a, o.b), (o.path.len(), o.a_bin_span, o.b_bin_span)))
            .collect();
        let mut seen: HashMap<(u32, u32), i64> = HashMap::new();
        let mut rows: Vec<(u32, u32, bool, i64, bool)> = Vec::new();
        for (p, x) in any_segs.iter().enumerate() {
            for y in &any_segs[p + 1..] {
                if y.minus != x.minus || y.lo >= x.hi {
                    break;
                }
                let ov = x.hi.min(y.hi) - x.lo.max(y.lo);
                if y.read == x.read || ov < truth_min {
                    continue;
                }
                let key = (x.read.min(y.read), x.read.max(y.read));
                if seen.get(&key).is_some_and(|&v| v >= ov) {
                    continue;
                }
                seen.insert(key, ov);
                rows.push((key.0, key.1, x.minus, ov, !(x.primary && y.primary)));
            }
        }
        let mut last: HashMap<(u32, u32), usize> = HashMap::new();
        for (k, r) in rows.iter().enumerate() {
            last.insert((r.0, r.1), k);
        }
        let mut n_rows = 0usize;
        for (k, r) in rows.iter().enumerate() {
            if last[&(r.0, r.1)] != k {
                continue;
            }
            n_rows += 1;
            let (cells, asp, bsp) = inor
                .get(&(r.0, r.1))
                .copied()
                .unwrap_or((0, (0, 0), (0, 0)));
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.0,
                r.1,
                names[r.0 as usize],
                names[r.1 as usize],
                if r.2 { '-' } else { '+' },
                r.3,
                r.4 as u8,
                (cells > 0) as u8,
                cells,
                asp.0,
                asp.1,
                bsp.0,
                bsp.1
            )?;
        }
        eprintln!("truth: {n_rows} read pairs >= {truth_min} bp -> {tp}");
    }

    std::fs::write(cout, bincode::serialize(&(&counts, &names, &oracle))?)
        .with_context(|| format!("writing {cout}"))?;
    eprintln!("wrote {cout}");
    Ok(())
}
