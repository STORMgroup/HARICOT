//! wtdbg2's SECOND layout level: fragments linked by the reads that run between them.
//!
//! wtdbg2 does not lay out the node graph's unitigs (`wtdbg.c:1137-1188`). After cleaning it
//! deletes every branching node (`mask_all_branching_nodes_graph`, `wtdbg-graph.h:3105`), so the
//! rest of the graph falls into unbranched FRAGMENTS (`gen_unitigs_graph`, `unitigs2frgs_graph`,
//! `:3543`). Links between fragments come from the READS: a read's placements, in read order,
//! name the fragments it runs through (`gen_lnks_graph`, `:3760`; `scan_rd_lnk_core`, `:3591`). A
//! link carried by fewer than [`FrgParams::weak`] reads is weak (`:3850`) and is cut where its
//! fragment end has other links (`cut_weak_branches_frg_graph`, `:3080`). Contigs are the linear
//! unique paths of what is left (`gen_contigs_graph`, `:3262`).
//!
//! Why it matters: a branching node is where a unitig stops -- a braid, a collapsed repeat, the
//! junction a chimeric read welds -- and the node graph only knows pairwise adjacency there. A read
//! that runs from one fragment into the next is molecule-level evidence of the continuation, and
//! reads are 20-40 kb where branch regions are a few nodes. Measured with a prototype from
//! correct overlaps of A. thaliana chr4 30x: 5.2 -> 14.6 Mb; on 33 bacterial sets: 60.8 -> 73.3%
//! of the 1 kb oracle.
//!
//! Needs the reads' node placements after cleaning: [`FuzzyGraph::read_threads`], recorded by the
//! build when `GraphParams::frg_layout` (or `read_threading`) is on. Placements on nodes that
//! cleaning removed or merged away are skipped. Forward-strand only, like the rest of the pipeline.
//!
//! Simplifications of wtdbg2's pass list (`wtdbg.c:1145-1186`): its link offsets, binary-link cut,
//! boomerang, bubble and repeat-detach passes on fragments are not ported; a transitive reduction
//! (drop x->z when x->y->z) and a short-fragment tip trim stand in for its Myers reduction and tip
//! trim. The masked branching nodes stay as one-node contigs (wtdbg2 deletes them), so the result
//! still partitions the live nodes and [`crate::layout::unitig_links`] applies unchanged.

use rustc_hash::{FxHashMap, FxHashSet};

use crate::graph::FuzzyGraph;
use crate::graph::ends::{End, EndAdj};
use crate::graph::fbg::Dir;
use crate::layout::unitig::Unitig;

/// Which fragment pairs one read links.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinkMode {
    /// wtdbg2's rule (`scan_rd_lnk_core`): one link per read, from the first fragment it visits to
    /// the last. Intermediate adjacencies come from other reads; the transitive reduction removes
    /// the long link again when they exist.
    FirstLast,
    /// Every consecutive pair of fragments the read visits.
    Consecutive,
}

/// Parameters of the fragment layout.
#[derive(Debug, Clone, Copy)]
pub struct FrgParams {
    /// A link carried by fewer reads is weak, and cut where its fragment end has other links.
    /// wtdbg2 uses its node-coverage minimum (`max_node_cov_sg = node_cov`, `wtdbg.c:966`): 3 by
    /// default, 2 when a genome size puts the depth at <= 40x (`wtdbg.c:923-932`).
    pub weak: u32,
    pub links: LinkMode,
    /// A fragment of at most this many nodes, linked on one side only, into a fragment end that
    /// has other links, is a tip and loses its link.
    pub tip_nodes: usize,
}

impl Default for FrgParams {
    fn default() -> Self {
        FrgParams {
            weak: 3,
            links: LinkMode::FirstLast,
            tip_nodes: 2,
        }
    }
}

/// What the fragment layout did, for the log.
#[derive(Debug, Default, Clone, Copy)]
pub struct FrgStats {
    pub masked: usize,
    pub fragments: usize,
    pub links: usize,
    pub weak_cut: usize,
    pub transitive: usize,
    pub tips: usize,
    pub contigs: usize,
    /// Contigs made of two or more fragments.
    pub joined: usize,
}

/// A read counts as leaving a fragment at its tail (entering at its head) if its last (first)
/// placement there is within this many nodes of the end. Node formation can leave a read's region
/// uncommitted at a fragment's end node (Omega conflict, node drop), and one node of slack keeps
/// such a read's link; a read that leaves from the MIDDLE of a fragment does not assert an
/// end-to-end continuation.
const END_SLACK: usize = 1;

/// One fragment visit of a read: (fragment, first index visited, last index visited).
type Visit = (usize, usize, usize);

/// Lay the cleaned graph out through fragments linked by reads. Returns contigs as [`Unitig`]s
/// (fragments concatenated in link order; each masked branching node a one-node contig), so every
/// live node is in exactly one. Falls back to plain unitigs if the graph has a reversed edge or no
/// read placements were recorded.
pub fn frg_contigs(graph: &FuzzyGraph, p: &FrgParams) -> (Vec<Unitig>, FrgStats) {
    let mut st = FrgStats::default();
    let reversed = graph.edges().iter().any(|e| {
        graph.is_live(e.from)
            && graph.is_live(e.to)
            && (e.from_dir != Dir::Fwd || e.to_dir != Dir::Fwd)
    });
    if reversed || graph.read_threads().is_empty() {
        log::warn!("frg layout: reversed edges or no read placements; plain unitigs instead");
        return (crate::layout::gen_unitigs(graph), st);
    }
    let ends = EndAdj::new(graph);
    let live: Vec<u32> = (0..graph.n_nodes() as u32)
        .filter(|&n| graph.is_live(n))
        .collect();
    // 1. Mask every branching node.
    let masked: FxHashSet<u32> = live
        .iter()
        .copied()
        .filter(|&n| ends.deg(n, End::L) > 1 || ends.deg(n, End::R) > 1)
        .collect();
    st.masked = masked.len();
    // 2. Unbranched fragments over the rest (forward strand: enter at L, leave at R).
    let next = |n: u32| -> Option<(u32, usize)> {
        match ends.at(n, End::R) {
            [(m, End::L, ei)] if !masked.contains(m) => Some((*m, *ei)),
            _ => None,
        }
    };
    let has_prev: FxHashSet<u32> = live
        .iter()
        .filter(|n| !masked.contains(n))
        .filter_map(|&n| next(n).map(|(m, _)| m))
        .collect();
    let mut frgs: Vec<Vec<u32>> = Vec::new();
    let mut spans: Vec<i64> = Vec::new();
    let mut frg_of: FxHashMap<u32, (usize, usize)> = FxHashMap::default();
    let len_of = |n: u32| graph.nodes()[n as usize].len_events as i64;
    let mut walk =
        |start: u32, frgs: &mut Vec<Vec<u32>>, frg_of: &mut FxHashMap<u32, (usize, usize)>| {
            let (fi, mut chain, mut span) = (frgs.len(), Vec::new(), 0i64);
            let mut x = Some((start, None::<usize>));
            while let Some((n, ei)) = x {
                if frg_of.contains_key(&n) {
                    break;
                }
                frg_of.insert(n, (fi, chain.len()));
                chain.push(n);
                span += len_of(n) + ei.map_or(0, |e| graph.edges()[e].offset as i64);
                x = next(n).map(|(m, e)| (m, Some(e)));
            }
            frgs.push(chain);
            spans.push(span);
        };
    for &n in &live {
        if !masked.contains(&n) && !has_prev.contains(&n) {
            walk(n, &mut frgs, &mut frg_of);
        }
    }
    for &n in &live {
        if !masked.contains(&n) && !frg_of.contains_key(&n) {
            walk(n, &mut frgs, &mut frg_of); // an unbranched cycle
        }
    }
    let nf = frgs.len();
    st.fragments = nf;
    // 3. Links from the reads.
    let mut cov: FxHashMap<(usize, usize), u32> = FxHashMap::default();
    let mut rids: Vec<&u32> = graph.read_threads().keys().collect();
    rids.sort_unstable();
    for rid in rids {
        let mut places = graph.read_threads()[rid].clone();
        places.sort_by_key(|t| t.beg);
        let mut seq: Vec<Visit> = Vec::new(); // in read order
        for t in &places {
            let Some(&(f, i)) = frg_of.get(&t.node) else {
                continue;
            };
            match seq.last_mut() {
                Some(last) if last.0 == f => last.2 = i,
                _ => seq.push((f, i, i)),
            }
        }
        let pairs: Vec<(Visit, Visit)> = match p.links {
            LinkMode::FirstLast if seq.len() >= 2 => vec![(seq[0], seq[seq.len() - 1])],
            LinkMode::FirstLast => Vec::new(),
            LinkMode::Consecutive => seq.windows(2).map(|w| (w[0], w[1])).collect(),
        };
        for ((fa, _, a_hi), (fb, b_lo, _)) in pairs {
            if fa != fb && a_hi + 1 + END_SLACK >= frgs[fa].len() && b_lo <= END_SLACK {
                *cov.entry((fa, fb)).or_insert(0) += 1;
            }
        }
    }
    st.links = cov.len();
    let mut out_l: Vec<FxHashMap<usize, u32>> = vec![FxHashMap::default(); nf];
    let mut in_l: Vec<FxHashMap<usize, u32>> = vec![FxHashMap::default(); nf];
    for (&(x, y), &c) in &cov {
        out_l[x].insert(y, c);
        in_l[y].insert(x, c);
    }
    let mut keys: Vec<(usize, usize)> = cov.keys().copied().collect();
    keys.sort_unstable();
    // 4. Weak links, cut where the fragment end has other links (decided on the full link set,
    //    then applied, as wtdbg2 collects its cuts first).
    let weak: Vec<(usize, usize)> = keys
        .iter()
        .copied()
        .filter(|&(x, y)| cov[&(x, y)] < p.weak && (out_l[x].len() > 1 || in_l[y].len() > 1))
        .collect();
    for &(x, y) in &weak {
        out_l[x].remove(&y);
        in_l[y].remove(&x);
    }
    st.weak_cut = weak.len();
    // 5. Transitive reduction: x->z goes when x->y->z exists.
    let mut tr: Vec<(usize, usize)> = Vec::new();
    for x in 0..nf {
        for &y in out_l[x].keys() {
            for &z in out_l[y].keys() {
                if z != y && z != x && out_l[x].contains_key(&z) {
                    tr.push((x, z));
                }
            }
        }
    }
    tr.sort_unstable();
    tr.dedup();
    for &(x, z) in &tr {
        out_l[x].remove(&z);
        in_l[z].remove(&x);
    }
    st.transitive = tr.len();
    // 6. Tips: a short fragment linked on one side only, into an end with other links.
    for f in 0..nf {
        if frgs[f].len() > p.tip_nodes {
            continue;
        }
        if in_l[f].is_empty() && out_l[f].len() == 1 {
            let y = *out_l[f].keys().next().expect("len 1");
            if in_l[y].len() > 1 {
                in_l[y].remove(&f);
                out_l[f].clear();
                st.tips += 1;
            }
        } else if out_l[f].is_empty() && in_l[f].len() == 1 {
            let x = *in_l[f].keys().next().expect("len 1");
            if out_l[x].len() > 1 {
                out_l[x].remove(&f);
                in_l[f].clear();
                st.tips += 1;
            }
        }
    }
    // 7. Contigs: linear unique paths.
    let uniq_next: Vec<Option<usize>> = (0..nf)
        .map(|x| match out_l[x].keys().collect::<Vec<_>>()[..] {
            [&y] if in_l[y].len() == 1 => Some(y),
            _ => None,
        })
        .collect();
    let has_p: FxHashSet<usize> = uniq_next.iter().flatten().copied().collect();
    let mut seen = vec![false; nf];
    let mut paths: Vec<Vec<usize>> = Vec::new();
    let chain = |f0: usize, seen: &mut Vec<bool>, paths: &mut Vec<Vec<usize>>| {
        let (mut path, mut x) = (Vec::new(), Some(f0));
        while let Some(f) = x {
            if seen[f] {
                break;
            }
            seen[f] = true;
            path.push(f);
            x = uniq_next[f];
        }
        paths.push(path);
    };
    for f in 0..nf {
        if !has_p.contains(&f) {
            chain(f, &mut seen, &mut paths);
        }
    }
    for f in 0..nf {
        if !seen[f] {
            chain(f, &mut seen, &mut paths);
        }
    }
    st.joined = paths.iter().filter(|p| p.len() > 1).count();
    let mut out: Vec<Unitig> = paths
        .iter()
        .map(|path| Unitig {
            nodes: path
                .iter()
                .flat_map(|&f| frgs[f].iter().map(|&n| (n, Dir::Fwd)))
                .collect(),
            len_events: path.iter().map(|&f| spans[f]).sum::<i64>().max(1) as u32,
        })
        .collect();
    let mut masked_sorted: Vec<u32> = masked.into_iter().collect();
    masked_sorted.sort_unstable();
    out.extend(masked_sorted.into_iter().map(|n| Unitig {
        nodes: vec![(n, Dir::Fwd)],
        len_events: len_of(n).max(1) as u32,
    }));
    st.contigs = paths.len();
    (out, st)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::fbg::Edge;
    use crate::graph::threading::{ReadPlacement, ReadThreads};

    fn e(from: u32, to: u32) -> Edge {
        Edge {
            from,
            from_dir: Dir::Fwd,
            to,
            to_dir: Dir::Fwd,
            cov: 3,
            support: 3.0,
            offset: 0,
        }
    }

    /// A collapsed repeat x: a1-a2-a3 -> x -> b1-b2 and d1-d2 -> x -> c1-c2. Reads run a -> x -> b
    /// (three) and d -> x -> c (three); one read a -> x -> c is weak. Contigs: a+b, d+c, and x alone.
    fn repeat_graph(extra_weak: bool) -> (FuzzyGraph, [u32; 10]) {
        let mut g = FuzzyGraph::new();
        let n: Vec<u32> = (0..10).map(|_| g.add_node(3, 80)).collect();
        let [a1, a2, a3, x, b1, b2, d1, d2, c1, c2] = n[..] else {
            unreachable!()
        };
        for (u, v) in [
            (a1, a2),
            (a2, a3),
            (a3, x),
            (x, b1),
            (b1, b2),
            (d1, d2),
            (d2, x),
            (x, c1),
            (c1, c2),
        ] {
            g.add_edge(e(u, v));
        }
        let mut rt = ReadThreads::default();
        let path = |ns: &[u32]| -> Vec<ReadPlacement> {
            ns.iter()
                .enumerate()
                .map(|(i, &node)| ReadPlacement {
                    node,
                    beg: 4 * i as u32,
                    end: 4 * i as u32 + 4,
                })
                .collect()
        };
        for r in 0..3 {
            rt.insert(r, path(&[a2, a3, x, b1, b2]));
            rt.insert(10 + r, path(&[d1, d2, x, c1]));
        }
        if extra_weak {
            rt.insert(20, path(&[a3, x, c1, c2]));
        }
        g.set_read_threads(rt);
        (g, [a1, a2, a3, x, b1, b2, d1, d2, c1, c2])
    }

    fn contig_sets(us: &[Unitig]) -> Vec<Vec<u32>> {
        let mut v: Vec<Vec<u32>> = us
            .iter()
            .map(|u| u.nodes.iter().map(|x| x.0).collect())
            .collect();
        v.sort();
        v
    }

    #[test]
    fn reads_link_fragments_through_a_masked_repeat() {
        for mode in [LinkMode::FirstLast, LinkMode::Consecutive] {
            for weak in [false, true] {
                let (g, [a1, a2, a3, x, b1, b2, d1, d2, c1, c2]) = repeat_graph(weak);
                let p = FrgParams {
                    links: mode,
                    ..FrgParams::default()
                };
                let (us, st) = frg_contigs(&g, &p);
                assert_eq!(st.masked, 1);
                assert_eq!(
                    contig_sets(&us),
                    vec![vec![a1, a2, a3, b1, b2], vec![x], vec![d1, d2, c1, c2]],
                    "{mode:?} weak={weak}"
                );
                assert_eq!(st.weak_cut, usize::from(weak));
                // Partition: every live node exactly once.
                let mut all: Vec<u32> = us
                    .iter()
                    .flat_map(|u| u.nodes.iter().map(|x| x.0))
                    .collect();
                all.sort_unstable();
                assert_eq!(all, (0..10).collect::<Vec<u32>>());
            }
        }
    }

    /// A read that leaves a fragment from its middle does not link it; without any placements
    /// the layout falls back to plain unitigs.
    #[test]
    fn mid_fragment_exits_do_not_link_and_no_placements_fall_back() {
        let (mut g, [a1, _a2, _a3, x, b1, ..]) = repeat_graph(false);
        let mut rt = ReadThreads::default();
        for r in 0..3 {
            // leaves a1-a2-a3 at a1, two nodes before its tail (one is the slack), then x, b1
            rt.insert(
                r,
                vec![
                    ReadPlacement {
                        node: a1,
                        beg: 0,
                        end: 4,
                    },
                    ReadPlacement {
                        node: x,
                        beg: 8,
                        end: 12,
                    },
                    ReadPlacement {
                        node: b1,
                        beg: 12,
                        end: 16,
                    },
                ],
            );
        }
        g.set_read_threads(rt);
        let (us, st) = frg_contigs(&g, &FrgParams::default());
        assert_eq!(st.joined, 0);
        assert_eq!(us.len(), 5, "4 fragments + the masked node");
        g.set_read_threads(ReadThreads::default());
        let (us, st) = frg_contigs(&g, &FrgParams::default());
        assert_eq!(st.fragments, 0, "fell back");
        assert_eq!(us.len(), crate::layout::gen_unitigs(&g).len());
    }
}
