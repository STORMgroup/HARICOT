//! The end bridges, a step of `haricot-layout layout`: dropped candidates that join two contig ends are
//! added back (`#[path]`-included like `cache.rs`).
//!
//! A contig end is the first / last `--end-nodes` nodes of a contig of at least `--min-steps` steps of the
//! build in OPDIR; its reads are every read the build threaded through those nodes (the THREADS dump's `T`
//! rows). A dropped candidate (`--drop`) is an END BRIDGE when one read sits at an end of one contig and the
//! other at an end of a different contig. Per end, the best bridge is the one with the highest stage-1 score
//! (`--scores`; a pair without one scores 0), an exact score tie going to the smaller tie key: the pair's two
//! read names in order, or with `--content-ties` (`haricot assemble`'s setting) the two reads' content
//! hashes in order (OPDIR/read_content.tsv). A bridge that is the best at both ends of one of its end pairs
//! is MUTUAL. The selections, each written as `<name>.drop.tsv` (the drop list minus the bridges it adds
//! back):
//!
//! * `q16_end_all`  every end bridge;
//! * `q16_mbb`      mutual-best over every bridge, score >= `--theta`;
//! * `q16_mbbo`     the same over ORIENTED bridges: a dovetail whose first read (the overlap starts deeper
//!   into it, and the other read runs on past it) sits at a contig's R end and the other at another's L end;
//! * `q16_mbbs`     oriented bridges that are mutual-best over every bridge;
//! * `q16_mbbc`     (`haricot assemble`'s) as `q16_mbbo`, with containments (no dovetail order) on any two
//!   ends;
//! * `q16_mbbx`     `q16_mbbc` plus every other bridge of score >= theta between two ends `q16_mbbc` joins;
//! * `q16_true_all`, `q16_true_end`: the drop list itself (no labels are read).
//!
//! Every collection is keyed rather than ordered, and the per-end pick is a strict order (score, then the
//! tie key), so no iteration order reaches an output; the drop lists are sorted by the two names. Scores
//! are parsed as `f64` (correctly rounded; the scores file holds each value's shortest round-trip digits).
//!
//! The options are [`parse_args`]'s (defaults `--end-nodes 3 --min-steps 2 --theta 0.5`); `haricot-layout`
//! passes `--only q16_mbbc` and, by default, `--content-ties`. The candidates' geometry comes from a
//! `haricot-ovl-filter --dump` TSV (DUMP) or from the overlaps themselves ([`run_on`], as `haricot-layout`
//! does): the same values (every dumped column it reads is an integer field of the overlap or a read
//! length), interned in the same order.
#![allow(dead_code)]

use anyhow::{Context, Result, bail, ensure};
use haricot::overlap::Overlap;
use rustc_hash::{FxHashMap, FxHashSet};

/// The command line.
#[derive(Debug, PartialEq)]
pub struct Args {
    pub opdir: String,
    pub threads: String,
    pub dump: String,
    pub drop: String,
    pub scores: String,
    pub out_dir: String,
    pub end_nodes: usize,
    pub min_steps: usize,
    pub theta: f64,
    pub only: Option<String>,
    pub content_ties: bool,
}

pub fn parse_args(argv: &[String]) -> Result<Args> {
    let usage = "usage: haricot-end-bridges OPDIR THREADS DUMP --drop S2.drop.tsv --scores S1.scores.tsv \
                 --out-dir D [--end-nodes 3] [--min-steps 2] [--theta 0.5] [--only NAME] [--content-ties]";
    let mut pos = Vec::new();
    let (mut drop, mut scores, mut out_dir, mut only) = (None, None, None, None);
    let (mut end_nodes, mut min_steps, mut theta, mut content_ties) =
        (3usize, 2usize, 0.5f64, false);
    let mut it = argv.iter().skip(1);
    while let Some(a) = it.next() {
        let mut val = |name: &str| -> Result<String> {
            it.next()
                .cloned()
                .with_context(|| format!("{name} needs a value; {usage}"))
        };
        match a.as_str() {
            "--drop" => drop = Some(val("--drop")?),
            "--scores" => scores = Some(val("--scores")?),
            "--out-dir" => out_dir = Some(val("--out-dir")?),
            "--only" => only = Some(val("--only")?),
            "--end-nodes" => end_nodes = val("--end-nodes")?.parse()?,
            "--min-steps" => min_steps = val("--min-steps")?.parse()?,
            "--theta" => theta = val("--theta")?.parse()?,
            "--content-ties" => content_ties = true,
            "--noext" | "--features" | "--report" | "--picks" | "--s2-scores" => {
                bail!("{a}: not supported; {usage}")
            }
            s if s.starts_with("--") => bail!("unknown option {s}; {usage}"),
            _ => pos.push(a.clone()),
        }
    }
    ensure!(pos.len() == 3, "{usage}");
    Ok(Args {
        opdir: pos[0].clone(),
        threads: pos[1].clone(),
        dump: pos[2].clone(),
        drop: drop.with_context(|| format!("--drop is required; {usage}"))?,
        scores: scores.with_context(|| format!("--scores is required; {usage}"))?,
        out_dir: out_dir.with_context(|| format!("--out-dir is required; {usage}"))?,
        end_nodes,
        min_steps,
        theta,
        only,
        content_ties,
    })
}

/// Read names, interned: a pair is two ids, compared by their names (string order).
#[derive(Default)]
pub struct Names {
    pub id: FxHashMap<String, u32>,
    pub name: Vec<String>,
}

impl Names {
    fn get(&mut self, s: &str) -> u32 {
        if let Some(&i) = self.id.get(s) {
            return i;
        }
        let i = self.name.len() as u32;
        self.id.insert(s.to_string(), i);
        self.name.push(s.to_string());
        i
    }

    /// `pkey(x, y)`: the two names in sort order.
    fn pair(&mut self, x: &str, y: &str) -> Pair {
        let (a, b) = (self.get(x), self.get(y));
        if x < y { (a, b) } else { (b, a) }
    }
}

/// A TSV with a header line, columns by name.
pub struct Tsv {
    text: String,
    cols: FxHashMap<String, usize>,
}

impl Tsv {
    fn read(path: &str) -> Result<Tsv> {
        let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
        let head = text
            .lines()
            .next()
            .with_context(|| format!("{path} is empty"))?;
        let cols = head
            .split('\t')
            .enumerate()
            .map(|(i, c)| (c.to_string(), i))
            .collect();
        Ok(Tsv { text, cols })
    }

    fn col(&self, name: &str, path: &str) -> Result<usize> {
        self.cols
            .get(name)
            .copied()
            .with_context(|| format!("{path}: no column {name}"))
    }

    /// The data rows, split on tabs.
    fn rows(&self) -> impl Iterator<Item = Vec<&str>> {
        self.text.lines().skip(1).map(|l| l.split('\t').collect())
    }
}

pub fn field<'a>(row: &[&'a str], i: usize, path: &str) -> Result<&'a str> {
    row.get(i)
        .copied()
        .with_context(|| format!("{path}: a row without column {i}"))
}

pub fn int(s: &str, path: &str) -> Result<i64> {
    s.parse()
        .with_context(|| format!("{path}: {s:?} is not an integer"))
}

/// A contig end: (contig, R side).
pub type End = (u32, bool);
/// A read pair: two interned names, the smaller name first (`pkey`).
pub type Pair = (u32, u32);
/// Bridges and their end pairs.
pub type Bridges = FxHashMap<Pair, Vec<(End, End)>>;
/// Per end, the best bridge's score and pair.
pub type Best = FxHashMap<End, (f64, Pair)>;
/// A unitig path row: step, file order, node (`None` when the build has no such node).
pub type PathRow = (i64, usize, Option<i64>);
/// A selection: its name and add-back set (the bridges it adds back).
pub type Arm = (&'static str, FxHashSet<Pair>);

/// `"{:,}"`: an integer with thousands separators.
pub fn commas(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

/// `%g`: six significant digits, trailing zeros dropped, exponent form outside 1e-4 .. 1e6.
pub fn py_g(x: f64) -> String {
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    if !x.is_finite() {
        return if x.is_nan() {
            "nan".into()
        } else if x > 0.0 {
            "inf".into()
        } else {
            "-inf".into()
        };
    }
    let e = format!("{x:.5e}");
    let (mant, exp) = e.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let trim = |s: &str| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if (-4..6).contains(&exp) {
        trim(&format!("{:.*}", (5 - exp).max(0) as usize, x))
    } else {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", trim(mant), exp.abs())
    }
}

/// One report line: the number of pairs a selection holds.
pub fn describe(name: &str, n: usize) -> String {
    format!("  {name:<12} {:>7} pairs", commas(n))
}

/// The selection, from the parsed inputs; returns the report lines and each selection's add-back set.
pub struct Selection {
    pub lines: Vec<String>,
    pub arms: Vec<Arm>,
    pub dropped: FxHashSet<Pair>,
}

/// Candidates held in memory: the overlaps, their reads' lengths and names (what `haricot-ovl-filter --dump`
/// writes of them).
#[derive(Clone, Copy)]
pub struct Mem<'a> {
    pub ovl: &'a [Overlap],
    pub counts: &'a [u32],
    pub names: &'a [String],
}

/// The candidates' geometry: every candidate pair, which read comes first in the genome (none for
/// a containment or a crossing), and the containments.
#[derive(Default)]
struct Geometry {
    cand: FxHashSet<Pair>,
    first: FxHashMap<Pair, u32>,
    contained: FxHashSet<Pair>,
}

impl Geometry {
    /// One candidate: its two read names and `[a_lo, a_hi, b_lo, b_hi, a_ev, b_ev, bin_events]`.
    fn add(&mut self, nm: &mut Names, x: &str, y: &str, v: [i64; 7]) {
        let [alo, ahi, blo, bhi, aev, bev, be] = v;
        let p = nm.pair(x, y);
        self.cand.insert(p);
        let (na, nb) = (aev.div_euclid(be) + 1, bev.div_euclid(be) + 1);
        let (ta, tb) = (na - 1 - ahi, nb - 1 - bhi);
        if alo >= blo && ta <= tb && (alo > blo || ta < tb) {
            self.first.insert(p, nm.get(x));
        } else if blo >= alo && tb <= ta && (blo > alo || tb < ta) {
            self.first.insert(p, nm.get(y));
        } else if (alo <= blo && ta <= tb) || (blo <= alo && tb <= ta) {
            self.contained.insert(p);
        }
    }
}

/// The selection over the candidates of `a.dump`.
pub fn run(a: &Args, nm: &mut Names) -> Result<Selection> {
    run_on(a, nm, None)
}

/// [`run`] with the candidates from memory (`Some`) instead of `a.dump`.
pub fn run_on(a: &Args, nm: &mut Names, mem: Option<Mem>) -> Result<Selection> {
    // read_names.tsv: index -> name
    let rn = format!("{}/read_names.tsv", a.opdir);
    let mut names: FxHashMap<i64, u32> = FxHashMap::default();
    for line in std::fs::read_to_string(&rn)
        .with_context(|| format!("reading {rn}"))?
        .lines()
    {
        let (i, n) = line
            .split_once('\t')
            .with_context(|| format!("{rn}: {line:?}"))?;
        ensure!(!n.contains('\t'), "{rn}: {line:?}");
        let id = nm.get(n);
        names.insert(int(i, &rn)?, id);
    }
    // the tie order: content hashes (OPDIR/read_content.tsv) or the names
    let mut content: FxHashMap<u32, u64> = FxHashMap::default();
    if a.content_ties {
        let p = format!("{}/read_content.tsv", a.opdir);
        ensure!(
            std::path::Path::new(&p).exists(),
            "--content-ties: no {p} (a cache without content hashes: re-run the detection)"
        );
        for line in std::fs::read_to_string(&p)?.lines().skip(1) {
            let (h, n) = line
                .split_once('\t')
                .with_context(|| format!("{p}: {line:?}"))?;
            ensure!(!n.contains('\t'), "{p}: {line:?}");
            let h = u64::from_str_radix(h, 16).with_context(|| format!("{p}: hash {h:?}"))?;
            let id = nm.get(n);
            content.insert(id, h);
        }
    }

    // contig ends: node -> {(contig, side)}
    let fp = format!("{}/fbg_nodes.tsv", a.opdir);
    let fbg = Tsv::read(&fp)?;
    let (cn, cr, cb) = (
        fbg.col("node", &fp)?,
        fbg.col("anchor_read", &fp)?,
        fbg.col("ev_beg", &fp)?,
    );
    let mut key2node: FxHashMap<(String, i64), i64> = FxHashMap::default();
    for row in fbg.rows() {
        let node = field(&row, cn, &fp)?;
        let n = int(node.get(1..).unwrap_or(""), &fp)?;
        let key = (
            field(&row, cr, &fp)?.to_string(),
            int(field(&row, cb, &fp)?, &fp)?,
        );
        key2node.insert(key, n);
    }
    let up = format!("{}/unitig_paths.tsv", a.opdir);
    let ut = Tsv::read(&up)?;
    let (uc, us, ur, ub) = (
        ut.col("contig", &up)?,
        ut.col("step", &up)?,
        ut.col("anchor_read", &up)?,
        ut.col("ev_beg", &up)?,
    );
    // rows by contig, then by step (a stable sort: equal steps keep the file order)
    let mut contig_id: FxHashMap<String, u32> = FxHashMap::default();
    let mut by_contig: Vec<Vec<PathRow>> = Vec::new();
    for (k, row) in ut.rows().enumerate() {
        let c = field(&row, uc, &up)?;
        let next = contig_id.len() as u32;
        let ci = *contig_id.entry(c.to_string()).or_insert(next);
        if ci as usize == by_contig.len() {
            by_contig.push(Vec::new());
        }
        let step = int(field(&row, us, &up)?, &up)?;
        let key = (
            field(&row, ur, &up)?.to_string(),
            int(field(&row, ub, &up)?, &up)?,
        );
        by_contig[ci as usize].push((step, k, key2node.get(&key).copied()));
    }
    let mut node_end: FxHashMap<i64, FxHashSet<End>> = FxHashMap::default();
    let mut n_ends = 0usize;
    for (ci, rows) in by_contig.iter_mut().enumerate() {
        rows.sort_by_key(|&(step, k, _)| (step, k));
        let ns: Vec<i64> = rows.iter().filter_map(|r| r.2).collect();
        if ns.len() < a.min_steps {
            continue;
        }
        let e = a.end_nodes.min(ns.len());
        for &x in &ns[..e] {
            node_end.entry(x).or_default().insert((ci as u32, false));
        }
        // the last e nodes (all of ns when e is 0)
        for &x in &ns[if e == 0 { 0 } else { ns.len() - e }..] {
            node_end.entry(x).or_default().insert((ci as u32, true));
        }
        n_ends += 2;
    }
    // reads at an end: read -> {(contig, side)}
    let mut read_end: FxHashMap<u32, FxHashSet<End>> = FxHashMap::default();
    for line in std::fs::read_to_string(&a.threads)
        .with_context(|| format!("reading {}", a.threads))?
        .lines()
    {
        if !line.starts_with("T\t") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        ensure!(f.len() == 5, "{}: a T row of {} fields", a.threads, f.len());
        if let Some(es) = node_end.get(&int(f[2], &a.threads)?)
            && !es.is_empty()
        {
            let rid = int(f[1], &a.threads)?;
            let r = *names
                .get(&rid)
                .with_context(|| format!("read {rid} is not in {rn}"))?;
            read_end.entry(r).or_default().extend(es.iter().copied());
        }
    }

    // the dropped candidates, the stage-1 scores, the candidates' geometry
    let mut dropped: FxHashSet<Pair> = FxHashSet::default();
    let dtext = std::fs::read_to_string(&a.drop).with_context(|| format!("reading {}", a.drop))?;
    ensure!(
        !dtext.trim().is_empty(),
        "{}: no columns to parse (an empty drop list)",
        a.drop
    );
    for line in dtext.lines().filter(|l| !l.is_empty()) {
        let (x, y) = line
            .split_once('\t')
            .with_context(|| format!("{}: {line:?}", a.drop))?;
        ensure!(!y.contains('\t'), "{}: {line:?}", a.drop);
        dropped.insert(nm.pair(x, y));
    }
    let st = Tsv::read(&a.scores)?;
    let (sa, sb, ss) = (
        st.col("name_a", &a.scores)?,
        st.col("name_b", &a.scores)?,
        st.col("score", &a.scores)?,
    );
    let mut s1: FxHashMap<Pair, f64> = FxHashMap::default();
    for row in st.rows() {
        let p = nm.pair(field(&row, sa, &a.scores)?, field(&row, sb, &a.scores)?);
        let v = field(&row, ss, &a.scores)?;
        let v: f64 = v
            .parse()
            .with_context(|| format!("{}: score {v:?}", a.scores))?;
        s1.insert(p, v);
    }
    let score = |p: &Pair| s1.get(p).copied().unwrap_or(0.0);
    let mut geo = Geometry::default();
    match mem {
        None => {
            let dt = Tsv::read(&a.dump)?;
            let dc: Vec<usize> = [
                "a",
                "b",
                "a_lo",
                "a_hi",
                "b_lo",
                "b_hi",
                "a_ev",
                "b_ev",
                "bin_events",
            ]
            .iter()
            .map(|c| dt.col(c, &a.dump))
            .collect::<Result<_>>()?;
            for row in dt.rows() {
                let (x, y) = (field(&row, dc[0], &a.dump)?, field(&row, dc[1], &a.dump)?);
                let v: Vec<i64> = dc[2..]
                    .iter()
                    .map(|&i| int(field(&row, i, &a.dump)?, &a.dump))
                    .collect::<Result<_>>()?;
                geo.add(nm, x, y, [v[0], v[1], v[2], v[3], v[4], v[5], v[6]]);
            }
        }
        // the dump's rows: one per overlap, in cache order (`haricot-ovl-filter --dump`)
        Some(m) => {
            for o in m.ovl {
                let v = [
                    o.a_bin_span.0 as i64,
                    o.a_bin_span.1 as i64,
                    o.b_bin_span.0 as i64,
                    o.b_bin_span.1 as i64,
                    m.counts[o.a as usize] as i64,
                    m.counts[o.b as usize] as i64,
                    o.bin_events as i64,
                ];
                geo.add(nm, &m.names[o.a as usize], &m.names[o.b as usize], v);
            }
        }
    }
    let Geometry {
        cand,
        first,
        contained,
    } = geo;

    // end bridges, oriented ones, and the containment-aware set
    let mut bridges: Bridges = FxHashMap::default();
    for &p in &dropped {
        let (Some(ex), Some(ey)) = (read_end.get(&p.0), read_end.get(&p.1)) else {
            continue;
        };
        let combos: Vec<(End, End)> = ex
            .iter()
            .flat_map(|&e1| ey.iter().map(move |&e2| (e1, e2)))
            .filter(|(e1, e2)| e1.0 != e2.0)
            .collect();
        if !combos.is_empty() {
            bridges.insert(p, combos);
        }
    }
    let mut obridges: Bridges = FxHashMap::default();
    for (&p, combos) in &bridges {
        let Some(&f0) = first.get(&p) else { continue };
        let oc: Vec<(End, End)> = combos
            .iter()
            .map(|&(e1, e2)| if p.0 == f0 { (e1, e2) } else { (e2, e1) })
            .filter(|(r, l)| r.1 && !l.1)
            .collect();
        if !oc.is_empty() {
            obridges.insert(p, oc);
        }
    }
    let mut cbridges = obridges.clone();
    for (&p, combos) in &bridges {
        if contained.contains(&p) {
            cbridges.insert(p, combos.clone());
        }
    }

    // per end, the best bridge: the highest score, an exact tie to the smaller tie key
    let tie_key = |p: Pair| -> Result<(u64, u64)> {
        let get = |r: u32| {
            content
                .get(&r)
                .copied()
                .with_context(|| format!("read {} has no content hash", nm.name[r as usize]))
        };
        let (x, y) = (get(p.0)?, get(p.1)?);
        Ok((x.min(y), x.max(y)))
    };
    let names_of = |p: Pair| (&nm.name[p.0 as usize], &nm.name[p.1 as usize]);
    let name_lt = |p: Pair, q: Pair| names_of(p) < names_of(q);
    let smaller = |p: Pair, q: Pair| -> Result<bool> {
        Ok(if a.content_ties {
            tie_key(p)? < tie_key(q)?
        } else {
            name_lt(p, q)
        })
    };
    let best_at = |br: &Bridges| -> Result<Best> {
        let mut best: Best = FxHashMap::default();
        for (&p, combos) in br {
            let v = score(&p);
            for &(e1, e2) in combos {
                for e in [e1, e2] {
                    let take = match best.get(&e) {
                        None => true,
                        Some(&(bv, bp)) => v > bv || (v == bv && smaller(p, bp)?),
                    };
                    if take {
                        best.insert(e, (v, p));
                    }
                }
            }
        }
        Ok(best)
    };
    let mutual_best = |br: &Bridges, best: &Best| {
        br.iter()
            .filter(|(p, combos)| {
                combos
                    .iter()
                    .any(|(e1, e2)| best[e1].1 == **p && best[e2].1 == **p)
            })
            .map(|(&p, _)| (p, score(&p)))
            .collect::<FxHashMap<Pair, f64>>()
    };
    let ball = best_at(&bridges)?;
    let mb = mutual_best(&bridges, &ball);
    let mbo = mutual_best(&obridges, &best_at(&obridges)?);
    let mbs = mutual_best(&obridges, &ball);
    let bc = best_at(&cbridges)?;
    let mbc = mutual_best(&cbridges, &bc);
    let at_least = |m: &FxHashMap<Pair, f64>, th: f64| -> FxHashSet<Pair> {
        m.iter()
            .filter(|(_, v)| **v >= th)
            .map(|(&p, _)| p)
            .collect()
    };
    let (mbb, mbbo, mbbs, mbbc) = (
        at_least(&mb, a.theta),
        at_least(&mbo, a.theta),
        at_least(&mbs, a.theta),
        at_least(&mbc, a.theta),
    );
    // the bundle: every other bridge (score >= theta) between two ends q16_mbbc joins
    let unordered = |(e1, e2): (End, End)| if e1 <= e2 { (e1, e2) } else { (e2, e1) };
    let mut joined: FxHashSet<(End, End)> = FxHashSet::default();
    for p in &mbbc {
        for &c in &cbridges[p] {
            if bc[&c.0].1 == *p && bc[&c.1].1 == *p {
                joined.insert(unordered(c));
            }
        }
    }
    let mut mbbx = mbbc.clone();
    for (&p, combos) in &cbridges {
        if score(&p) >= a.theta && combos.iter().any(|&c| joined.contains(&unordered(c))) {
            mbbx.insert(p);
        }
    }

    // the report: the pairs each selection holds
    let end_all: FxHashSet<Pair> = bridges.keys().copied().collect();
    let oriented: FxHashSet<Pair> = obridges.keys().copied().collect();
    let mut lines = vec![format!(
        "{} contig ends, {} reads at an end; {} dropped of {} candidates",
        commas(n_ends),
        commas(read_end.len()),
        commas(dropped.len()),
        commas(cand.len())
    )];
    lines.push(describe("end_all", end_all.len()));
    lines.push(describe("oriented", oriented.len()));
    for th in [0.0, 0.3, 0.5, 0.7, 0.9] {
        let g = py_g(th);
        lines.push(describe(&format!("mbb>={g}"), at_least(&mb, th).len()));
        lines.push(describe(&format!("mbbo>={g}"), at_least(&mbo, th).len()));
        lines.push(describe(&format!("mbbs>={g}"), at_least(&mbs, th).len()));
        lines.push(describe(&format!("mbbc>={g}"), at_least(&mbc, th).len()));
    }
    lines.push(describe(&format!("mbbx>={}", py_g(a.theta)), mbbx.len()));
    let empty = FxHashSet::default();
    let arms = vec![
        ("q16_true_all", empty.clone()),
        ("q16_true_end", empty),
        ("q16_end_all", end_all),
        ("q16_mbb", mbb),
        ("q16_mbbo", mbbo),
        ("q16_mbbs", mbbs),
        ("q16_mbbc", mbbc),
        ("q16_mbbx", mbbx),
    ];
    Ok(Selection {
        lines,
        arms,
        dropped,
    })
}

/// A selection's drop list: the dropped pairs not in `add`, sorted by the two names, `name_a<TAB>name_b`
/// per line.
pub fn drop_list(sel: &Selection, nm: &Names, add: &FxHashSet<Pair>) -> String {
    let mut keep: Vec<_> = sel
        .dropped
        .iter()
        .filter(|p| !add.contains(*p))
        .map(|p| {
            (
                nm.name[p.0 as usize].as_str(),
                nm.name[p.1 as usize].as_str(),
            )
        })
        .collect();
    keep.sort_unstable();
    let mut s = String::new();
    for (x, y) in keep {
        s.push_str(x);
        s.push('\t');
        s.push_str(y);
        s.push('\n');
    }
    s
}

/// Write the selections' drop lists (`only`: that one alone) to `out_dir/<name>.drop.tsv`; returns each
/// selection's name and list.
pub fn write_arms(
    sel: &Selection,
    nm: &Names,
    only: Option<&str>,
    out_dir: &str,
) -> Result<Vec<(&'static str, String)>> {
    let arms: Vec<&Arm> = match only {
        Some(o) => {
            let v: Vec<_> = sel.arms.iter().filter(|(n, _)| *n == o).collect();
            ensure!(!v.is_empty(), "--only {o}: no such selection");
            v
        }
        None => sel.arms.iter().collect(),
    };
    let mut out = Vec::new();
    for (name, add) in arms {
        let text = drop_list(sel, nm, add);
        let path = format!("{out_dir}/{name}.drop.tsv");
        std::fs::write(&path, &text).with_context(|| format!("writing {path}"))?;
        out.push((*name, text));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    /// A candidate as the test's dump rows describe it (`a_lo 5 a_hi 9 b_lo 0 b_hi 4`, bin 10).
    fn ovl(a: u32, b: u32) -> Overlap {
        Overlap {
            a,
            b,
            strand: haricot::params::Strand::Forward,
            n_anchors: 1,
            matched: 1,
            cnt: 1,
            aligned_bins: 1,
            bin_events: 10,
            gap: 0,
            slope: 1.0,
            a_bin_span: (5, 9),
            b_bin_span: (0, 4),
            path: Vec::new(),
            z: 5.0,
        }
    }

    #[test]
    fn report_formats_are_pythons() {
        assert_eq!(commas(0), "0");
        assert_eq!(commas(999), "999");
        assert_eq!(commas(1000), "1,000");
        assert_eq!(commas(12_075), "12,075");
        assert_eq!(commas(1_234_567), "1,234,567");
        for (x, want) in [
            (0.0, "0"),
            (0.3, "0.3"),
            (0.5, "0.5"),
            (0.9, "0.9"),
            (1.0, "1"),
            (0.55, "0.55"),
            (1e-5, "1e-05"),
            (123456.0, "123456"),
            (1234567.0, "1.23457e+06"),
            (0.0001, "0.0001"),
        ] {
            assert_eq!(py_g(x), want, "{x}");
        }
        assert_eq!(describe("end_all", 16), "  end_all           16 pairs");
        assert_eq!(describe("mbbx>=0.5", 1234), "  mbbx>=0.5      1,234 pairs");
    }

    #[test]
    fn arguments_are_the_production_options() {
        let a = parse_args(&s(&[
            "haricot-end-bridges",
            "op",
            "thr",
            "dump",
            "--drop",
            "d",
            "--scores",
            "s",
            "--out-dir",
            "o",
            "--only",
            "q16_mbbc",
            "--content-ties",
        ]))
        .unwrap();
        assert_eq!((a.end_nodes, a.min_steps, a.theta), (3, 2, 0.5));
        assert!(a.content_ties && a.only.as_deref() == Some("q16_mbbc"));
        assert!(
            parse_args(&s(&[
                "haricot-end-bridges",
                "op",
                "thr",
                "dump",
                "--drop",
                "d",
                "--scores",
                "s"
            ]))
            .is_err()
        );
        let research = s(&[
            "haricot-end-bridges",
            "op",
            "thr",
            "dump",
            "--drop",
            "d",
            "--scores",
            "s",
            "--out-dir",
            "o",
            "--noext",
            "x",
        ]);
        assert!(parse_args(&research).is_err());
    }

    /// Builds an OPDIR + inputs where two bridges tie on score at one end: the pick follows the names
    /// without `--content-ties` and the content hashes with it.
    #[test]
    fn score_ties_follow_names_or_content() {
        let dir = std::env::temp_dir().join(format!("eb_select_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let w = |name: &str, text: &str| std::fs::write(dir.join(name), text).unwrap();
        // contigs u0 (nodes 0-3) and u1 (nodes 4-7); reads: a at u0's R end, b and c at u1's L end
        w("read_names.tsv", "0\tra\n1\trb\n2\trc\n3\trd\n");
        let mut nodes =
            String::from("node\tanchor_read\tev_beg\tev_end\tcov\tanchor_coord_events\n");
        let mut paths = String::from(
            "contig\tstep\tanchor_read\tev_beg\tev_end\tdir\tadvance_events\tcoord_events\n",
        );
        for k in 0..8 {
            nodes += &format!("n{k}\tanchor\t{}\t0\t1\t1\n", k * 100);
            let c = if k < 4 { "u0" } else { "u1" };
            paths += &format!("{c}\t{}\tanchor\t{}\t0\t+\t0\t0\n", k % 4, k * 100);
        }
        w("fbg_nodes.tsv", &nodes);
        w("unitig_paths.tsv", &paths);
        // T rows: read node ...: a on node 3 (u0's R end), b and c on node 4 (u1's L end), d on node 1
        w(
            "thr.tsv",
            "C\tx\tx\tx\tx\tx\tx\nT\t0\t3\t0\t1\nT\t1\t4\t0\t1\nT\t2\t4\t0\t1\nT\t3\t1\t0\t1\n",
        );
        w("drop.tsv", "ra\trb\nra\trc\nrd\trb\n");
        // equal scores for (ra, rb) and (ra, rc): a tie at u0's R end and at u1's L end
        w(
            "scores.tsv",
            "name_a\tname_b\tcls\tscore\nra\trb\t?\t0.8\nra\trc\t?\t0.8\nrd\trb\t?\t0.1\n",
        );
        // dovetails: the first read is the one the overlap starts deeper into
        let dump = "a\tb\tz\tn_anchors\tmatched\tcnt\taligned_bins\tgap\tbin_events\tslope\ta_lo\ta_hi\tb_lo\tb_hi\ta_ev\tb_ev\n";
        let row =
            |x: &str, y: &str| format!("{x}\t{y}\t5\t1\t1\t1\t1\t0\t10\t1\t5\t9\t0\t4\t99\t99\n");
        w(
            "dump.tsv",
            &(dump.to_string() + &row("ra", "rb") + &row("ra", "rc") + &row("rd", "rb")),
        );
        let d = dir.to_str().unwrap();
        let run_with = |content: Option<&str>| -> String {
            if let Some(c) = content {
                w("read_content.tsv", c);
            }
            let mut argv = s(&[
                "haricot-end-bridges",
                d,
                &format!("{d}/thr.tsv"),
                &format!("{d}/dump.tsv"),
            ]);
            argv.extend(s(&[
                "--drop",
                &format!("{d}/drop.tsv"),
                "--scores",
                &format!("{d}/scores.tsv"),
            ]));
            argv.extend(s(&["--out-dir", d, "--only", "q16_mbbc"]));
            if content.is_some() {
                argv.push("--content-ties".into());
            }
            let a = parse_args(&argv).unwrap();
            let mut nm = Names::default();
            let sel = run(&a, &mut nm).unwrap();
            // the same candidates from memory: the same report and the same drop lists
            let names: Vec<String> = ["ra", "rb", "rc", "rd"].map(String::from).to_vec();
            let ovls = [ovl(0, 1), ovl(0, 2), ovl(3, 1)];
            let m = Mem {
                ovl: &ovls,
                counts: &[99; 4],
                names: &names,
            };
            let mut nm2 = Names::default();
            let sel2 = run_on(&a, &mut nm2, Some(m)).unwrap();
            assert_eq!(sel.lines, sel2.lines);
            assert_eq!(sel.arms.len(), sel2.arms.len());
            for ((n1, add1), (n2, add2)) in sel.arms.iter().zip(&sel2.arms) {
                assert_eq!(n1, n2);
                assert_eq!(drop_list(&sel, &nm, add1), drop_list(&sel2, &nm2, add2));
            }
            let (_, add) = sel.arms.iter().find(|(n, _)| *n == "q16_mbbc").unwrap();
            let mut v: Vec<String> = add
                .iter()
                .map(|p| format!("{}-{}", nm.name[p.0 as usize], nm.name[p.1 as usize]))
                .collect();
            v.sort();
            v.join(",")
        };
        // names: ("ra", "rb") < ("ra", "rc")
        assert_eq!(run_with(None), "ra-rb");
        // content: rc's hash below rb's makes (ra, rc) the smaller key
        assert_eq!(
            run_with(Some(
                "content\tname\n0000000000000010\tra\n00000000000000f0\trb\n0000000000000020\trc\n0000000000000030\trd\n"
            )),
            "ra-rc"
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
