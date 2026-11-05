//! Layout: collapse the cleaned graph into unitigs.
//!
//! Port of wtdbg2's `gen_unitigs_graph` (`wtdbg-graph.h:3128`): walk maximal
//! non-branching paths to form unitigs.
//!
//! **This is orientation-aware, and must stay that way.** A node has two physical
//! ends; an edge attaches to one of them ([`End`]). Two nodes merge into a chain
//! only when the joining edge is the *sole* edge on **both** of the ends it
//! touches. Undirected degree 2 is NOT sufficient: a node whose two edges land on
//! the *same* end is a convergence, and merging through it invents a path that
//! does not exist (see [`crate::graph::ends`]).
//!
//! Two invariants make that class of bug impossible, and are asserted in tests:
//!
//! 1. the unitigs **partition** the live nodes (every node in exactly one);
//! 2. every live edge is either *internal* to a unitig or emitted as exactly one
//!    [`UnitigLink`].
//!
//! Together they force the unitig graph to be a true contraction of the FBG, so
//! component and cycle counts are preserved. Verify changes against
//! `gfatools asm -u`: it must agree on component and cycle counts.

use rustc_hash::{FxHashMap, FxHashSet};

use crate::graph::FuzzyGraph;
use crate::graph::ends::{End, EndAdj, entry_end, exit_end, from_end, to_end};
use crate::graph::fbg::Dir;
use crate::graph::threading;

/// A maximal non-branching path through the graph.
#[derive(Debug, Clone, Default)]
pub struct Unitig {
    /// Ordered `(node, orientation)` along the path. `Fwd` = entered at the node's
    /// L end and left by its R end.
    pub nodes: Vec<(u32, Dir)>,
    /// Span in events: `sum(node lengths) + sum(internal edge offsets)`. The
    /// offsets are negative for overlaps, so the shared bins between consecutive
    /// nodes are counted once — matching `gfatools asm -u` (`LN:i:672` for two
    /// 384-event nodes overlapping by 96).
    pub len_events: u32,
}

/// A link between two unitig ends, induced by the FBG edge that joins them.
#[derive(Debug, Clone)]
pub struct UnitigLink {
    pub from: usize,
    /// `Fwd` = the link leaves `from`'s right end (GFA `+`).
    pub from_dir: Dir,
    pub to: usize,
    /// `Fwd` = the link enters `to`'s left end (GFA `+`).
    pub to_dir: Dir,
    pub overlap_events: u32,
}

/// Collapse the live graph into unitigs (maximal non-branching paths).
///
/// The unitigs partition the live nodes: a junction node is **not** duplicated into
/// each incident unitig — duplicating it counts a degree-`d` junction `d` times and
/// inflates the reported lengths.
pub fn gen_unitigs(graph: &FuzzyGraph) -> Vec<Unitig> {
    gen_unitigs_threaded(graph, 1)
}

/// [`gen_unitigs`], with the read-corroboration floor for threaded repeat traversal
/// (`GraphParams::thread_min_reads`). `1` reproduces `gen_unitigs` exactly.
pub fn gen_unitigs_threaded(graph: &FuzzyGraph, thread_min_reads: u32) -> Vec<Unitig> {
    let ends = EndAdj::new(graph);
    // Read-threaded repeat traversal. Empty unless `GraphParams::read_threading` was on,
    // in which case this whole mechanism is inert and the walk is exactly as before.
    //
    // KNOWN LIMITATION: a repeat node threaded by k disjoint paths can still only belong
    // to ONE unitig, so k-1 of them break there anyway. Resolving that properly means
    // UNZIPPING the node into k copies at build time; this extends whichever path
    // reaches it first (deterministic — nodes are scanned in id order — but arbitrary).
    // It is still strictly better than stopping all k, which is what happens without it.
    let tr = graph.transitions();
    let tadj = threading::ThreadAdj::new(graph);
    let len_of = |n: u32| graph.nodes()[n as usize].len_events as i64;
    let offset_of = |ei: usize| graph.edges()[ei].offset as i64;
    let next = |n: u32, d: Dir| ends.next(n, d);
    let prev = |n: u32, d: Dir| ends.prev(n, d);

    let mut visited: FxHashSet<u32> = FxHashSet::default();
    let mut unitigs = Vec::new();
    for s in 0..graph.n_nodes() as u32 {
        if !graph.is_live(s) || visited.contains(&s) {
            continue;
        }
        // Rewind to the chain's start. `seen` stops a circular chain from looping
        // forever; it then starts at an arbitrary node of the cycle.
        let mut start = (s, Dir::Fwd);
        let mut seen: FxHashSet<u32> = FxHashSet::default();
        seen.insert(s);
        while let Some((m, dm, _)) = prev(start.0, start.1) {
            if !seen.insert(m) {
                break;
            }
            start = (m, dm);
        }
        // Walk forward, accumulating span as sum(len) + sum(offset).
        let mut nodes = vec![start];
        visited.insert(start.0);
        let mut span = len_of(start.0);
        let mut cur = start;
        let mut prev_node: Option<u32> = None;
        loop {
            let step = match next(cur.0, cur.1) {
                Some(s) => Some(s),
                // The chain hit a branch. A branch is a locus several places in the
                // genome map onto, and which continuation is right is a property of the
                // READS passing through, not of the node or of either edge. If exactly
                // one continuation is read-threaded — and no other entry claims it —
                // carry on through instead of stopping. Inert unless threading was
                // recorded, so the default layout is unchanged.
                None if !tr.is_empty() => prev_node
                    // (a) leave a branching node by the read-threaded continuation
                    .and_then(|p| threading::threaded_step(&tadj, tr, p, cur.0, thread_min_reads))
                    // (b) or ENTER a junction, but only when we can already see that
                    // its exit will be threaded. Without that look-ahead the walk
                    // could claim a repeat node it then cannot leave, blocking a path
                    // that could have gone through.
                    .or_else(|| {
                        let (m, _dm, _ei) = ends.step_out(cur.0, cur.1)?;
                        threading::threaded_step(&tadj, tr, cur.0, m, thread_min_reads).map(|_| m)
                    })
                    .and_then(|w| tadj.edge(cur.0, w).map(|ei| (w, Dir::Fwd, ei))),
                None => None,
            };
            let Some((m, dm, ei)) = step else { break };
            if visited.contains(&m) {
                break; // closed a cycle, or a repeat another path already claimed
            }
            visited.insert(m);
            span += offset_of(ei) + len_of(m);
            nodes.push((m, dm));
            prev_node = Some(cur.0);
            cur = (m, dm);
        }
        unitigs.push(Unitig {
            nodes,
            len_events: span.max(1) as u32,
        });
    }
    unitigs
}

/// `node -> index of the unitig containing it`. Total, given the partition invariant.
fn unitig_of(unitigs: &[Unitig]) -> FxHashMap<u32, usize> {
    let mut m = FxHashMap::default();
    for (ui, u) in unitigs.iter().enumerate() {
        for &(n, _) in &u.nodes {
            m.insert(n, ui);
        }
    }
    m
}

/// A unitig's two free ends as `(node, end)`.
fn free_ends(u: &Unitig) -> ((u32, End), (u32, End)) {
    let (n1, d1) = u.nodes[0];
    let (nk, dk) = *u.nodes.last().unwrap();
    ((n1, entry_end(d1)), (nk, exit_end(dk)))
}

/// Build the unitig graph: every live FBG edge not internal to a unitig becomes one
/// link between two unitig ends.
///
/// Note a "two unitigs link if they share an endpoint node" rule does not work here:
/// it relies on junction duplication, and links all `k` unitigs at a junction
/// pairwise — a `k`-clique where the truth is a star, inventing cycles.
pub fn unitig_links(graph: &FuzzyGraph, unitigs: &[Unitig]) -> Vec<UnitigLink> {
    let uof = unitig_of(unitigs);
    let fe: Vec<((u32, End), (u32, End))> = unitigs.iter().map(free_ends).collect();
    let mut links = Vec::new();
    for e in graph.edges() {
        if !graph.is_live(e.from) || !graph.is_live(e.to) || e.from == e.to {
            continue;
        }
        let (ua, ub) = (uof[&e.from], uof[&e.to]);
        let (a_end, b_end) = ((e.from, from_end(e.from_dir)), (e.to, to_end(e.to_dir)));
        // An edge on an interior end is internal to its unitig: skip it. Only the
        // two free ends can carry a link.
        let from_dir = if a_end == fe[ua].1 {
            Dir::Fwd // leaves ua's right end
        } else if a_end == fe[ua].0 {
            Dir::Rev // leaves ua's left end
        } else {
            continue;
        };
        let to_dir = if b_end == fe[ub].0 {
            Dir::Fwd // enters ub's left end
        } else if b_end == fe[ub].1 {
            Dir::Rev // enters ub's right end
        } else {
            continue;
        };
        links.push(UnitigLink {
            from: ua,
            from_dir,
            to: ub,
            to_dir,
            overlap_events: if e.offset < 0 { (-e.offset) as u32 } else { 0 },
        });
    }
    links
}

/// A contig: an ordered chain of unitig indices laid out by best-overlap.
#[derive(Debug, Clone, Default)]
pub struct Contig {
    /// Indices into the `gen_unitigs` result, in layout order.
    pub unitigs: Vec<usize>,
    /// Total length in events (sum of unitig lengths).
    pub len_events: u32,
}

/// Lay the unitig graph out into linear contigs by greedy **mutual-best**
/// extension: at each junction node (where several unitig ends meet) every unitig
/// end proposes its highest-coverage neighbour, and two ends are joined only when
/// the choice is *mutual*.
///
/// DEVIATION from wtdbg2: its second-level layout (`gen_contigs_graph`,
/// `wtdbg-graph.h:3262`) traces linear-unique paths over a *fragment* graph whose
/// nodes are whole unitigs; that fragment graph is not built here, so the chaining
/// rule is mutual-best-by-coverage over the unitig graph directly.
///
/// The effect is to linearise a bushy unitig graph: the dominant tiling of each
/// locus chains into one long contig, while lower-coverage side branches (including
/// the phase-shifted double-tiling artefacts left behind by transitive reduction /
/// bubble popping) fall out as short separate contigs.
pub fn gen_contigs(graph: &FuzzyGraph) -> Vec<Contig> {
    gen_contigs_from(graph, &gen_unitigs(graph))
}

/// [`gen_contigs`] over a unitig set the caller already has.
///
/// A `Contig` stores unitig INDICES, so it is only meaningful against the exact slice it
/// was built from. Recomputing the unitigs inside here would give indices into a different
/// list whenever `gen_unitigs_threaded` returns a different partition (an out-of-bounds
/// panic). Taking the slice makes that coupling impossible to get wrong.
pub fn gen_contigs_from(graph: &FuzzyGraph, unitigs: &[Unitig]) -> Vec<Contig> {
    let nu = unitigs.len();
    if nu == 0 {
        return Vec::new();
    }
    // Mean node coverage of each unitig (the "weight" for best-overlap).
    let weight: Vec<f64> = unitigs
        .iter()
        .map(|u| {
            let s: u64 = u
                .nodes
                .iter()
                .map(|&(n, _)| graph.nodes()[n as usize].cov as u64)
                .sum();
            s as f64 / u.nodes.len().max(1) as f64
        })
        .collect();

    // Unitig ends: 0 = left (head), 1 = right (tail). Adjacency comes from
    // `unitig_links`, not from a shared endpoint *node*: junctions are not
    // duplicated into every incident unitig, so there is no shared node to key on.
    let mut at_end: FxHashMap<(usize, u8), Vec<(usize, u8)>> = FxHashMap::default();
    // Link overlap, keyed both ways, so a chain can subtract the bins it shares.
    type EndPair = ((usize, u8), (usize, u8));
    let mut ov: FxHashMap<EndPair, u32> = FxHashMap::default();
    for l in unitig_links(graph, unitigs) {
        let a = (l.from, if l.from_dir == Dir::Fwd { 1u8 } else { 0u8 });
        let b = (l.to, if l.to_dir == Dir::Fwd { 0u8 } else { 1u8 });
        at_end.entry(a).or_default().push(b);
        at_end.entry(b).or_default().push(a);
        ov.insert((a, b), l.overlap_events);
        ov.insert((b, a), l.overlap_events);
    }

    // best[(ui,end)] = the highest-weight unitig end joined to it (excluding self).
    // None if the end is free (a dead end).
    let mut best: FxHashMap<(usize, u8), (usize, u8)> = FxHashMap::default();
    for ui in 0..nu {
        for end in 0u8..2 {
            let mut pick: Option<((usize, u8), f64)> = None;
            if let Some(cands) = at_end.get(&(ui, end)) {
                for &(uj, ej) in cands {
                    if uj == ui {
                        continue;
                    }
                    let w = weight[uj];
                    if pick.map(|(_, pw)| w > pw).unwrap_or(true) {
                        pick = Some(((uj, ej), w));
                    }
                }
            }
            if let Some((cand, _)) = pick {
                best.insert((ui, end), cand);
            }
        }
    }

    // Mutual-best links: join (ui,ei)-(uj,ej) only when each is the other's best.
    let mutual = |a: (usize, u8)| -> Option<(usize, u8)> {
        best.get(&a).copied().filter(|&b| best.get(&b) == Some(&a))
    };

    // Walk chains of mutually-linked unitig ends into contigs.
    let mut used = vec![false; nu];
    let mut contigs = Vec::new();
    for start in 0..nu {
        if used[start] {
            continue;
        }
        // Walk backward from `start`'s head to find a chain end, then walk forward.
        // Simpler: collect the whole chain `start` belongs to by following mutual
        // links from both ends.
        let mut chain = std::collections::VecDeque::new();
        chain.push_back(start);
        used[start] = true;
        // Overlap consumed by the links joining this chain: consecutive unitigs
        // share `overlap_events` bins, so summing their lengths double-counts them.
        let mut ov_sum: i64 = 0;
        // extend from the tail end (end 1) forward
        let mut cur = (start, 1u8);
        while let Some((uj, ej)) = mutual(cur) {
            if used[uj] {
                break;
            }
            ov_sum += ov.get(&(cur, (uj, ej))).copied().unwrap_or(0) as i64;
            used[uj] = true;
            chain.push_back(uj);
            cur = (uj, 1 - ej); // leave by the opposite end
        }
        // extend from the head end (end 0) backward
        let mut cur = (start, 0u8);
        while let Some((uj, ej)) = mutual(cur) {
            if used[uj] {
                break;
            }
            ov_sum += ov.get(&(cur, (uj, ej))).copied().unwrap_or(0) as i64;
            used[uj] = true;
            chain.push_front(uj);
            cur = (uj, 1 - ej);
        }
        let unitigs_in: Vec<usize> = chain.into_iter().collect();
        // Subtract the shared bins, exactly as `Unitig::len_events` does for the
        // nodes inside a unitig. Summing unitig lengths raw inflates every contig by
        // ~(n_unitigs - 1) * overlap — at bin96 that is 96 events (~118 bp) per join,
        // and it lands straight in contigs.gfa's `LN:i:` and contig N50.
        let raw: i64 = unitigs_in
            .iter()
            .map(|&u| unitigs[u].len_events as i64)
            .sum();
        let len_events = (raw - ov_sum).max(1) as u32;
        contigs.push(Contig {
            unitigs: unitigs_in,
            len_events,
        });
    }
    contigs
}

#[cfg(test)]
mod threading_tests {
    use super::*;
    use crate::graph::fbg::Edge;
    use crate::graph::threading::Transitions;

    fn e(from: u32, to: u32) -> Edge {
        Edge {
            from,
            from_dir: Dir::Fwd,
            to,
            to_dir: Dir::Fwd,
            cov: 2,
            support: 2.0,
            offset: 0,
        }
    }

    /// REGRESSION: a `Contig` stores unitig INDICES, so `gen_contigs` must never
    /// recompute a different unitig set than the caller holds. Threaded traversal makes
    /// the two partitions differ, and an index into the wrong one runs out of bounds.
    #[test]
    fn contigs_index_into_the_unitigs_they_were_built_from() {
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(4, 80);
        }
        for (a, b) in [(0, 2), (1, 2), (2, 3), (2, 4)] {
            g.add_edge(e(a, b));
        }
        let mut tr = Transitions::default();
        tr.insert((0, 2, 3), 1);
        tr.insert((1, 2, 4), 1);
        g.set_transitions(tr);
        let unitigs = gen_unitigs_threaded(&g, 1);
        let contigs = gen_contigs_from(&g, &unitigs);
        for c in &contigs {
            for &ui in &c.unitigs {
                assert!(
                    ui < unitigs.len(),
                    "contig indexes unitig {ui} of {}",
                    unitigs.len()
                );
            }
        }
        // And the unthreaded pair must still agree with each other.
        let plain = gen_unitigs(&g);
        for c in gen_contigs(&g) {
            for ui in c.unitigs {
                assert!(ui < plain.len());
            }
        }
    }

    /// A repeat node entered and left two ways. Without threading BOTH chains stop at
    /// the repeat, so the layout is four short unitigs. With reads threading the node,
    /// one chain carries straight through — which is the entire point: a longer CORRECT
    /// contig, not a longer one.
    #[test]
    fn threading_carries_a_unitig_through_a_repeat() {
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(4, 80);
        }
        for (a, b) in [(0, 2), (1, 2), (2, 3), (2, 4)] {
            g.add_edge(e(a, b));
        }
        let plain = gen_unitigs(&g);
        let longest_plain = plain.iter().map(|u| u.nodes.len()).max().unwrap();
        assert_eq!(longest_plain, 1, "every chain stops at the branch");

        let mut tr = Transitions::default();
        tr.insert((0, 2, 3), 1);
        tr.insert((1, 2, 4), 1);
        g.set_transitions(tr);
        let threaded = gen_unitigs(&g);
        let longest = threaded.iter().map(|u| u.nodes.len()).max().unwrap();
        assert!(
            longest >= 3,
            "one path should run 0 -> 2 -> 3 (or 1 -> 2 -> 4), got {longest}"
        );
        // The repeat belongs to exactly one unitig; the other path still breaks. That is
        // the documented limitation — unzipping is what would fix it.
        let uses_2: usize = threaded
            .iter()
            .filter(|u| u.nodes.iter().any(|&(n, _)| n == 2))
            .count();
        assert_eq!(uses_2, 1, "node 2 is claimed once");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::fbg::{Dir, Edge, FuzzyGraph};

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

    /// A linear chain is one unitig covering every node.
    #[test]
    fn linear_chain_is_one_unitig() {
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(3, 72);
        }
        for i in 0..4 {
            g.add_edge(e(i, i + 1));
        }
        let u = gen_unitigs(&g);
        assert_eq!(u.len(), 1);
        assert_eq!(u[0].nodes.len(), 5);
        assert_eq!(u[0].len_events, 5 * 72);
    }

    /// A Y junction yields three unitigs. Note *which* three: node 1's L end holds
    /// only the edge from 0, so 0 and 1 merge; its R end forks to 2 and 3, which
    /// stay separate. The junction node 1 belongs to exactly one unitig, not to each
    /// incident one.
    #[test]
    fn y_junction_splits_unitigs() {
        let mut g = FuzzyGraph::new();
        for _ in 0..4 {
            g.add_node(3, 72);
        }
        // 0-1 (deg1 ends), 1 is junction with 2 and 3.
        g.add_edge(e(0, 1));
        g.add_edge(e(1, 2));
        g.add_edge(e(1, 3));
        let u = gen_unitigs(&g);
        assert_eq!(u.len(), 3, "three arms off the junction");
        let mut sizes: Vec<usize> = u.iter().map(|x| x.nodes.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![1, 1, 2], "0+1 merge; 2 and 3 stand alone");
        assert_eq!(total_nodes(&u), 4, "junction node is not duplicated");
    }

    fn total_nodes(u: &[Unitig]) -> usize {
        u.iter().map(|x| x.nodes.len()).sum()
    }

    /// REGRESSION: A -> B <- C. Both edges enter B's LEFT end, so B has undirected
    /// degree 2 but is a convergence, not a chain, and must not be merged through.
    /// An undirected walk (`deg(n) == 2` ⇒ interior) folds A-B-C into ONE unitig,
    /// inventing a path from A to C that the graph does not contain — which shows up
    /// as the unitig graph *gaining* components and cycles, something a contraction
    /// cannot do.
    #[test]
    fn same_end_convergence_is_not_merged() {
        let mut g = FuzzyGraph::new();
        for _ in 0..3 {
            g.add_node(3, 72);
        }
        g.add_edge(e(0, 1)); // 0's R -> 1's L
        g.add_edge(e(2, 1)); // 2's R -> 1's L  (same end of 1)
        let u = gen_unitigs(&g);
        assert_eq!(
            u.len(),
            3,
            "A->B<-C is a convergence on B's L end: no merge is legal (undirected \
             code merged it into 1)"
        );
        assert_eq!(total_nodes(&u), 3);
    }

    /// Contraction invariants: unitigs partition the live nodes, and every live edge
    /// is either internal to a unitig or becomes exactly one link. Together they
    /// force component and cycle counts to be preserved.
    #[test]
    fn unitigs_partition_nodes_and_conserve_edges() {
        let mut g = FuzzyGraph::new();
        for _ in 0..8 {
            g.add_node(3, 72);
        }
        // chain 0-1-2, junction at 2 forking to 3 and 5; 3-4 and 5-6 chains that
        // rejoin at 7 — a bubble, so the graph has a cycle to preserve.
        for (a, b) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (2, 5),
            (3, 4),
            (5, 6),
            (4, 7),
            (6, 7),
        ] {
            g.add_edge(e(a, b));
        }
        let u = gen_unitigs(&g);
        let links = unitig_links(&g, &u);

        // 1. partition: every live node in exactly one unitig.
        let mut seen: Vec<u32> = u
            .iter()
            .flat_map(|x| x.nodes.iter().map(|&(n, _)| n))
            .collect();
        let n_before = seen.len();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), n_before, "a node appears in two unitigs");
        assert_eq!(
            seen.len(),
            g.n_live_nodes(),
            "unitigs must cover every node"
        );

        // 2. edge conservation: internal + linked == every live edge. A linear
        // unitig of k nodes internalises k-1 edges.
        let internal: usize = u.iter().map(|x| x.nodes.len() - 1).sum();
        assert_eq!(
            internal + links.len(),
            g.n_edges(),
            "every edge must be internal to a unitig or emitted as one link"
        );
    }

    /// Length counts the shared bins once: two 72-event nodes overlapping by 24
    /// span 120, not 144 — summing node lengths raw over-reports. `gfatools asm -u`
    /// agrees (LN:i:672 for 2x384 at 96 overlap).
    #[test]
    fn unitig_length_subtracts_overlap() {
        let mut g = FuzzyGraph::new();
        g.add_node(3, 72);
        g.add_node(3, 72);
        g.add_edge(Edge {
            from: 0,
            from_dir: Dir::Fwd,
            to: 1,
            to_dir: Dir::Fwd,
            cov: 3,
            support: 3.0,
            offset: -24, // 24 events of overlap
        });
        let u = gen_unitigs(&g);
        assert_eq!(u.len(), 1);
        assert_eq!(u[0].len_events, 72 - 24 + 72, "overlap counted once");
    }

    /// A single linear chain lays out as exactly one contig spanning it.
    #[test]
    fn linear_chain_is_one_contig() {
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(3, 72);
        }
        for i in 0..4 {
            g.add_edge(e(i, i + 1));
        }
        let c = gen_contigs(&g);
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].len_events, 5 * 72);
    }

    /// At a junction, best-overlap merges the backbone with its higher-coverage
    /// continuation into one contig and leaves the low-coverage side branch as a
    /// separate short contig.
    #[test]
    fn best_overlap_picks_heavier_continuation() {
        let mut g = FuzzyGraph::new();
        // Backbone 0-1-2-3 (cov 9), strong continuation 3-4-5 (cov 9), weak side 3-6 (cov 2).
        for _ in 0..4 {
            g.add_node(9, 72); // 0,1,2,3
        }
        g.add_node(9, 72); // 4
        g.add_node(9, 72); // 5
        g.add_node(2, 72); // 6 (weak)
        for i in 0..3 {
            g.add_edge(e(i, i + 1)); // 0-1-2-3
        }
        g.add_edge(e(3, 4));
        g.add_edge(e(4, 5));
        g.add_edge(e(3, 6));
        let c = gen_contigs(&g);
        // One long contig (backbone + strong continuation) + one short side contig.
        let lens: Vec<u32> = c.iter().map(|x| x.len_events).collect();
        let longest = *lens.iter().max().unwrap();
        assert!(
            longest >= 6 * 72,
            "backbone should merge with the strong arm, got {lens:?}"
        );
        assert!(
            c.len() >= 2,
            "the weak side branch is a separate contig, got {lens:?}"
        );
    }

    /// REGRESSION: a contig must subtract the overlap of the links joining its
    /// unitigs, exactly as `Unitig::len_events` does for the nodes inside a unitig.
    /// Summing unitig lengths raw inflated every contig by ~(n_unitigs-1)*overlap,
    /// straight into contigs.gfa's LN:i: and contig N50.
    ///
    /// The other tests all use `e()`, which hardcodes `offset: 0` — so no overlap
    /// exists to lose and they cannot catch this.
    ///
    /// The fixture needs BOTH a junction (so more than one unitig exists, and a
    /// *link* joins them) and nonzero offsets (so there is an overlap to lose). An
    /// unbranched chain is a single unitig with no links, where `ov_sum` is
    /// trivially 0 and the assertion passes on `gen_unitigs`' span math alone.
    ///
    /// Node 3's R end branches (to 4 and to the weak 6), splitting the graph into
    /// A=[0,1,2,3], B=[4,5], C=[6]. Mutual-best joins A-B over the 3→4 link.
    #[test]
    fn contig_length_subtracts_link_overlap() {
        fn eo(from: u32, to: u32, offset: i32, cov: u32) -> Edge {
            Edge {
                from,
                from_dir: Dir::Fwd,
                to,
                to_dir: Dir::Fwd,
                cov,
                support: cov as f32,
                offset,
            }
        }
        const OV: i32 = -24; // every join shares 24 events
        let mut g = FuzzyGraph::new();
        for _ in 0..6 {
            g.add_node(9, 72); // 0..5: the strong backbone + continuation
        }
        g.add_node(2, 72); // 6: weak side branch, loses best-overlap to [4,5]

        g.add_edge(eo(0, 1, OV, 9)); // ┐
        g.add_edge(eo(1, 2, OV, 9)); // ├ inside unitig A
        g.add_edge(eo(2, 3, OV, 9)); // ┘
        g.add_edge(eo(3, 4, OV, 9)); // the LINK joining unitig A to unitig B
        g.add_edge(eo(4, 5, OV, 9)); // inside unitig B
        g.add_edge(eo(3, 6, OV, 2)); // weak side branch -> its own contig

        let u = gen_unitigs(&g);
        assert_eq!(
            u.len(),
            3,
            "the junction at 3's R end must split A=[0,1,2,3], B=[4,5], C=[6]"
        );

        let c = gen_contigs(&g);
        let mut lens: Vec<u32> = c.iter().map(|x| x.len_events).collect();
        lens.sort_unstable();
        // A = 4*72 - 3*24 = 216, B = 2*72 - 1*24 = 120, C = 72.
        // A and B are joined by a link overlapping 24, so the contig spans
        // 216 + 120 - 24 = 312. Summing raw would give 336.
        assert_eq!(
            lens,
            vec![72, 312],
            "the joining link's 24-event overlap must be subtracted once"
        );
    }
}
