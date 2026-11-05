//! The contig-end rescue's admission (`haricot-rescue admit`) on a base overlap set, shared with
//! `haricot-layout` (`#[path]`-included like `cache.rs`,
//! which the including program must also declare at its root). The method is documented at
//! `haricot-rescue`'s `admit`.
#![allow(dead_code)]

use std::io::Write;

use anyhow::Result;
use haricot::overlap::Overlap;
use rustc_hash::FxHashMap;

/// The mutual best among the pool pairs `allow_p` allows (`pool_detail.tsv` at `detail`, the
/// pool's overlaps at `pool_p`), added to `ovl` (the base set, then sorted by pair); with
/// `symmetric`, ties by the reads' tie rank (`content`, the base cache's content hashes).
/// `out_pairs`: the admitted pairs, one per line. Returns the report line (`admit: ...`).
pub fn admit(
    detail: &str,
    pool_p: &str,
    allow_p: &str,
    ovl: &mut Vec<Overlap>,
    content: Option<&[u64]>,
    symmetric: bool,
    out_pairs: Option<&str>,
) -> Result<String> {
    let mut allow: rustc_hash::FxHashSet<(u32, u32)> = Default::default();
    for l in std::fs::read_to_string(allow_p)?.lines() {
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.len() >= 2 {
            let (x, y): (u32, u32) = (f[0].parse()?, f[1].parse()?);
            allow.insert((x.min(y), x.max(y)));
        }
    }
    type Row = ((u32, u32), f64, Vec<(String, String)>);
    let mut rows: Vec<Row> = Vec::new();
    for l in std::fs::read_to_string(detail)?.lines() {
        if l.starts_with('#') || l.starts_with("a\t") {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        let (x, y): (u32, u32) = (f[0].parse()?, f[1].parse()?);
        if !allow.contains(&(x.min(y), x.max(y))) {
            continue;
        }
        let eps = f[4]
            .split(';')
            .filter_map(|ep| ep.split_once('-'))
            .map(|(p, q)| (p.to_string(), q.to_string()))
            .collect();
        rows.push(((x, y), f[2].parse()?, eps));
    }
    if symmetric {
        // ties among equal scores by the tie rank (the reads' signal content-hash order)
        let rank = crate::cache::tie_rank(content)?;
        let key = |&(x, y): &(u32, u32)| {
            let (p, q) = (rank[x as usize], rank[y as usize]);
            (p.min(q), p.max(q))
        };
        rows.sort_by_key(|r| key(&r.0));
    }
    let mut best: FxHashMap<&str, (f64, (u32, u32))> = FxHashMap::default();
    for (p, sc, eps) in &rows {
        for (e1, e2) in eps {
            for e in [e1.as_str(), e2.as_str()] {
                if best.get(e).is_none_or(|b| *sc > b.0) {
                    best.insert(e, (*sc, *p));
                }
            }
        }
    }
    let mut adm: Vec<(u32, u32)> = rows
        .iter()
        .filter(|(p, _, eps)| {
            eps.iter()
                .any(|(e1, e2)| best[e1.as_str()].1 == *p && best[e2.as_str()].1 == *p)
        })
        .map(|r| r.0)
        .collect();
    adm.sort_unstable();
    if let Some(p) = out_pairs {
        let mut w = std::io::BufWriter::new(std::fs::File::create(p)?);
        for (x, y) in &adm {
            writeln!(w, "{x}\t{y}")?;
        }
    }
    let pool: Vec<Overlap> = bincode::deserialize(&std::fs::read(pool_p)?)?;
    let by_pair: FxHashMap<(u32, u32), &Overlap> = pool.iter().map(|o| ((o.a, o.b), o)).collect();
    let add: Vec<Overlap> = adm
        .iter()
        .filter_map(|p| by_pair.get(p).map(|o| (*o).clone()))
        .collect();
    let line = format!(
        "admit: {} allowed pool pairs, {} admitted; {} base overlaps + {}",
        rows.len(),
        adm.len(),
        ovl.len(),
        add.len()
    );
    ovl.extend(add);
    ovl.sort_by_key(|o| (o.a, o.b));
    Ok(line)
}
