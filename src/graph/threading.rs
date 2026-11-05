//! Read threading: path-level evidence the FBG builds and then discards.
//!
//! `build_graph_annotated` threads each read through the committed regions it belongs to
//! — a per-read *path* through the node graph — and then reduces every such path to its
//! consecutive pairs, keeping only an edge coverage count. The longer context is thrown
//! away at exactly the point where it is the only thing that can resolve a repeat.
//!
//! A branching node is a locus that several places in the genome map onto. Which
//! continuation is correct is not a property of the node or of either edge — both edges
//! are real — it is a property of the *reads that pass through*. `min_edge_cov`,
//! `min_edge_support` and `max_node_cov` are all blind to it by construction, which is
//! why a relative edge gate reaches only 1.4-1.8x selectivity on misjoining edges
//! (measured on the development sets).
//!
//! A **transition** is an ordered node triple `(u, v, w)` that at least one read threads
//! in that order. It is threshold-free: a read either threads it or it does not.
//!
//! Recorded only when [`GraphParams::read_threading`] is on, so the default graph and
//! every existing layout are unchanged.

use rustc_hash::FxHashMap;

use crate::graph::ends::{End, EndAdj};
use crate::graph::fbg::FuzzyGraph;

/// Where one read sits in the node graph: the node, and the read's own bin range in it.
///
/// This is the placement the build stage computes and then reduces to edges. Two reads
/// that share a node have a relative offset that the node fixes — `beg` is in each read's
/// OWN bin coordinates — which is a geometric constraint on any overlap between them,
/// independent of whether that overlap was detected.
#[derive(Debug, Clone, Copy)]
pub struct ReadPlacement {
    pub node: u32,
    pub beg: u32,
    pub end: u32,
}

/// Read id → every node it was threaded through, in read order.
pub type ReadThreads = FxHashMap<u32, Vec<ReadPlacement>>;

/// Ordered node triples, each mapped to **how many reads** thread it.
///
/// The count is corroboration, not a score: a triple threaded by one read rests on one
/// molecule, which may be chimeric or misplaced. Requiring more is
/// [`GraphParams::thread_min_reads`].
pub type Transitions = FxHashMap<(u32, u32, u32), u32>;

/// How much of the branching in a graph read threading could actually resolve.
#[derive(Debug, Default, Clone, Copy)]
pub struct BranchStats {
    /// `(in-edge, node)` pairs where the node's far end branches.
    pub entries: usize,
    /// Entries with exactly one read-supported continuation — resolvable.
    pub resolvable: usize,
    /// Entries with more than one — a genuine ambiguity no threading can settle.
    pub ambiguous: usize,
    /// Entries with none — the incoming path is not continued by any read, so the
    /// correct action is to STOP rather than to pick an edge.
    pub dead: usize,
    /// Branching ends seen (deg > 1), for context.
    pub branch_ends: usize,
}

/// Count, for every way of entering a branching node, whether the reads passing through
/// single out one continuation.
///
/// Edges are Fwd/Fwd throughout while the overlapper is forward-strand only
/// (`overlap::mod`), so triples are keyed on plain node ids in read order. This will
/// need orientation once cross-strand edges exist.
pub fn branch_stats(graph: &FuzzyGraph, tr: &Transitions, min_reads: u32) -> BranchStats {
    let ends = EndAdj::new(graph);
    // node -> the nodes preceding / following it, from the live edge set.
    let mut out: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    let mut inc: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for e in graph.edges() {
        if !graph.is_live(e.from) || !graph.is_live(e.to) || e.from == e.to {
            continue;
        }
        out.entry(e.from).or_default().push(e.to);
        inc.entry(e.to).or_default().push(e.from);
    }
    let mut s = BranchStats::default();
    for n in 0..graph.n_nodes() as u32 {
        if !graph.is_live(n) {
            continue;
        }
        if ends.deg(n, End::R) > 1 {
            s.branch_ends += 1;
        }
        let Some(succs) = out.get(&n) else { continue };
        if succs.len() < 2 {
            continue; // not a branch: nothing to resolve
        }
        let empty = Vec::new();
        for &u in inc.get(&n).unwrap_or(&empty) {
            s.entries += 1;
            let k = succs
                .iter()
                .filter(|&&w| tr.get(&(u, n, w)).is_some_and(|&c| c >= min_reads))
                .count();
            match k {
                0 => s.dead += 1,
                1 => s.resolvable += 1,
                _ => s.ambiguous += 1,
            }
        }
    }
    s
}

/// Successor / predecessor lists over the live edge set, built once so a walk does not
/// rescan every edge per step.
pub struct ThreadAdj {
    succ: FxHashMap<u32, Vec<u32>>,
    pred: FxHashMap<u32, Vec<u32>>,
    edge_of: FxHashMap<(u32, u32), usize>,
    /// False when any live edge is reversed; threading then disables itself rather than
    /// mis-orienting.
    pub fwd_only: bool,
}

impl ThreadAdj {
    pub fn new(graph: &FuzzyGraph) -> Self {
        use crate::graph::fbg::Dir;
        let mut a = ThreadAdj {
            succ: FxHashMap::default(),
            pred: FxHashMap::default(),
            edge_of: FxHashMap::default(),
            fwd_only: true,
        };
        for (ei, e) in graph.edges().iter().enumerate() {
            if !graph.is_live(e.from) || !graph.is_live(e.to) || e.from == e.to {
                continue;
            }
            if e.from_dir != Dir::Fwd || e.to_dir != Dir::Fwd {
                a.fwd_only = false;
            }
            a.succ.entry(e.from).or_default().push(e.to);
            a.pred.entry(e.to).or_default().push(e.from);
            a.edge_of.entry((e.from, e.to)).or_insert(ei);
        }
        a
    }

    /// Index of the edge `from -> to`, if it is live.
    pub fn edge(&self, from: u32, to: u32) -> Option<usize> {
        self.edge_of.get(&(from, to)).copied()
    }
}

/// The unique read-threaded continuation from `u -> v`, if there is exactly one and it
/// is unique in BOTH directions.
///
/// Extending `u -> v` to `w` requires:
///   1. some read threads `(u, v, w)`;
///   2. `w` is the only successor of `v` that any read continues to from `u`;
///   3. `u` is the only predecessor of `v` that any read arrives from before `w`.
///
/// Condition 3 is what stops two different paths from both claiming `w` — without it the
/// winner would be decided by iteration order. Nothing here is a threshold: a read either
/// threads a triple or it does not.
///
/// Forward-strand only, like the rest of the pipeline: every edge is `Fwd/Fwd`
/// (`overlap::mod` pairs same-strand seed hits only), so triples are plain node ids in
/// read order. Returns `None` when any live edge is reversed, so this silently disables
/// itself rather than mis-orienting the day cross-strand edges land.
pub fn threaded_step(
    adj: &ThreadAdj,
    tr: &Transitions,
    u: u32,
    v: u32,
    min_reads: u32,
) -> Option<u32> {
    if !adj.fwd_only {
        return None;
    }
    let ok = |a: u32, b: u32, c: u32| tr.get(&(a, b, c)).is_some_and(|&n| n >= min_reads);
    let empty = Vec::new();
    let succ = adj.succ.get(&v).unwrap_or(&empty);
    let mut cands = succ.iter().copied().filter(|&w| ok(u, v, w));
    let w = cands.next()?;
    if cands.next().is_some() {
        return None; // ambiguous: the reads do not single out a continuation
    }
    // Mutual uniqueness: no other predecessor of `v` threads into `w` either.
    if adj
        .pred
        .get(&v)
        .unwrap_or(&empty)
        .iter()
        .any(|&p| p != u && ok(p, v, w))
    {
        return None;
    }
    Some(w)
}

/// What [`unzip`] did.
#[derive(Debug, Default, Clone, Copy)]
pub struct UnzipStats {
    /// Repeat nodes split.
    pub nodes: usize,
    /// Copies created — one per read-threaded path through a split node.
    pub copies: usize,
    /// Nodes fully consumed: every edge moved onto a copy, so the original is deleted.
    pub consumed: usize,
}

/// Split each repeat node into one copy per read-threaded path through it.
///
/// [`threaded_step`] lets a walk carry ONE chain through a branch; the other *k*−1 paths
/// still break there, because a node can belong to only one unitig. Unzipping is the
/// principled form: a node that *k* disjoint read-threaded paths pass through is a locus
/// with *k* genomic copies collapsed onto it, so it becomes *k* nodes and every path gets
/// its own. Just as important, entries whose continuation the reads do NOT single out are
/// left on the original node instead of being claimed by whichever walk arrives first.
///
/// The pairing rule is exactly [`threaded_step`]'s — `u` matches `w` when some read
/// threads `(u, v, w)` at least `min_reads` times, `w` is the only such successor of `u`,
/// and `u` the only such predecessor of `w`. Only nodes with **two or more** matched
/// pairs are split: one pair is what `threaded_step` already handles, and splitting it
/// would duplicate sequence for nothing.
///
/// Each copy carries the coverage of the reads that thread it, so a *k*-copy repeat ends
/// up as *k* nodes at roughly `cov / k` rather than one node at `cov` — which is what
/// makes the duplicated sequence honest rather than inflationary.
///
/// Node ids are append-only, so existing ids and the `annos` they index stay valid.
/// Forward-strand only, like the rest of the pipeline: returns immediately if any live
/// edge is reversed.
pub fn unzip(
    graph: &mut FuzzyGraph,
    annos: &mut Vec<crate::graph::NodeAnno>,
    min_reads: u32,
) -> UnzipStats {
    use crate::graph::fbg::{Dir, Edge};

    let mut st = UnzipStats::default();
    let tr = graph.transitions().clone();
    if tr.is_empty() {
        return st;
    }
    let adj = ThreadAdj::new(graph);
    if !adj.fwd_only {
        return st;
    }
    let ok = |a: u32, b: u32, c: u32| tr.get(&(a, b, c)).is_some_and(|&n| n >= min_reads);

    // Decide everything against the adjacency above, then apply in one pass, so the
    // graph is not mutated while it is still being read.
    let mut drop: FxHashMap<(u32, u32), ()> = FxHashMap::default();
    let mut add: Vec<Edge> = Vec::new();
    let mut consume: Vec<u32> = Vec::new();
    let empty = Vec::new();

    for v in 0..graph.n_nodes() as u32 {
        if !graph.is_live(v) {
            continue;
        }
        let preds = adj.pred.get(&v).unwrap_or(&empty);
        let succs = adj.succ.get(&v).unwrap_or(&empty);
        if preds.len() < 2 && succs.len() < 2 {
            continue;
        }
        let mut matched: Vec<(u32, u32, u32)> = Vec::new();
        for &u in preds {
            let mut it = succs.iter().copied().filter(|&w| ok(u, v, w));
            let Some(w) = it.next() else { continue };
            if it.next().is_some() {
                continue; // this entry has several supported continuations
            }
            if preds.iter().any(|&p| p != u && ok(p, v, w)) {
                continue; // another entry claims the same continuation
            }
            matched.push((u, w, tr.get(&(u, v, w)).copied().unwrap_or(0)));
        }
        if matched.len() < 2 {
            continue;
        }
        st.nodes += 1;
        let len = graph.nodes()[v as usize].len_events;
        let anno = annos[v as usize];
        for &(u, w, n) in &matched {
            let off_in = adj.edge(u, v).map(|i| graph.edges()[i].offset).unwrap_or(0);
            let off_out = adj.edge(v, w).map(|i| graph.edges()[i].offset).unwrap_or(0);
            let copy = graph.add_node(n, len);
            annos.push(anno);
            graph.set_node_support(copy, n as f32);
            st.copies += 1;
            drop.insert((u, v), ());
            drop.insert((v, w), ());
            add.push(Edge {
                from: u,
                from_dir: Dir::Fwd,
                to: copy,
                to_dir: Dir::Fwd,
                cov: n,
                support: n as f32,
                offset: off_in,
            });
            add.push(Edge {
                from: copy,
                from_dir: Dir::Fwd,
                to: w,
                to_dir: Dir::Fwd,
                cov: n,
                support: n as f32,
                offset: off_out,
            });
        }
        // Every incident edge routed onto a copy ⇒ nothing is left of the original.
        if matched.len() == preds.len() && matched.len() == succs.len() {
            consume.push(v);
            st.consumed += 1;
        }
    }

    if st.nodes == 0 {
        return st;
    }
    graph
        .edges_mut()
        .retain(|e| !drop.contains_key(&(e.from, e.to)));
    for e in add {
        graph.add_edge(e);
    }
    for v in consume {
        graph.nodes_mut()[v as usize].del = true;
    }
    st
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::fbg::{Dir, Edge};

    fn edge(from: u32, to: u32) -> Edge {
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

    /// The case the whole module exists for: a node entered from two sides and leaving
    /// to two sides. Edge counts cannot say which pairing is right — both edges are
    /// real. The reads can.
    #[test]
    fn threading_resolves_a_repeat_a_count_cannot() {
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(4, 80);
        }
        // 0 -> 2 -> 3 and 1 -> 2 -> 4, sharing the repeat node 2.
        for (a, b) in [(0, 2), (1, 2), (2, 3), (2, 4)] {
            g.add_edge(edge(a, b));
        }
        let mut tr: Transitions = Transitions::default();
        tr.insert((0, 2, 3), 2);
        tr.insert((1, 2, 4), 2);
        let s = branch_stats(&g, &tr, 1);
        assert_eq!(s.entries, 2, "two ways in");
        assert_eq!(s.resolvable, 2, "each entry has one supported continuation");
        assert_eq!(s.ambiguous, 0);
        assert_eq!(s.dead, 0);

        // With no read spanning the repeat, the same topology is unresolvable — and the
        // right answer is to stop, not to guess.
        let s = branch_stats(&g, &Transitions::default(), 1);
        assert_eq!(s.resolvable, 0);
        assert_eq!(s.dead, 2);
    }

    /// Unzipping is what `threaded_step` cannot do: with a repeat entered and left two
    /// ways, the walk can only carry ONE chain through, so the other still breaks.
    /// Splitting the node lets both run, and the coverage goes with them.
    #[test]
    fn unzip_splits_a_repeat_so_every_threaded_path_runs() {
        use crate::graph::NodeAnno;
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(4, 80);
        }
        // 0 -> 2 -> 3 and 1 -> 2 -> 4, sharing the repeat node 2 (cov 4 = 2 + 2).
        for (a, b) in [(0, 2), (1, 2), (2, 3), (2, 4)] {
            g.add_edge(edge(a, b));
        }
        let mut tr = Transitions::default();
        tr.insert((0, 2, 3), 2);
        tr.insert((1, 2, 4), 2);
        g.set_transitions(tr);
        let mut annos: Vec<NodeAnno> = (0..5)
            .map(|_| NodeAnno {
                anchor_read: 0,
                window: 0,
                ev_beg: 0,
                ev_end: 80,
            })
            .collect();

        let st = unzip(&mut g, &mut annos, 1);
        assert_eq!(st.nodes, 1, "one repeat split");
        assert_eq!(st.copies, 2, "one copy per threaded path");
        assert_eq!(st.consumed, 1, "every edge moved, so the original is gone");
        assert!(!g.is_live(2), "the collapsed node is retired");
        assert_eq!(annos.len(), g.n_nodes(), "annos stay indexable by node id");

        // Both paths are now simple chains, so the plain walk runs them end to end.
        let unitigs = gen_unitigs_for_test(&g);
        let long = unitigs.iter().filter(|u| u.len() >= 3).count();
        assert_eq!(long, 2, "BOTH paths run through, not one: {unitigs:?}");
    }

    /// Node ids are append-only so `annos` stays valid, and a copy carries only the
    /// coverage of the reads that thread it — a 2-copy repeat becomes two cov-2 nodes,
    /// not two cov-4 ones, which is what stops duplication from inflating the layout.
    #[test]
    fn unzip_splits_coverage_between_the_copies() {
        use crate::graph::NodeAnno;
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(4, 80);
        }
        for (a, b) in [(0, 2), (1, 2), (2, 3), (2, 4)] {
            g.add_edge(edge(a, b));
        }
        let mut tr = Transitions::default();
        tr.insert((0, 2, 3), 3);
        tr.insert((1, 2, 4), 1);
        g.set_transitions(tr);
        let mut annos: Vec<NodeAnno> = (0..5)
            .map(|_| NodeAnno {
                anchor_read: 0,
                window: 0,
                ev_beg: 0,
                ev_end: 80,
            })
            .collect();
        unzip(&mut g, &mut annos, 1);
        let covs: Vec<u32> = (5..g.n_nodes() as u32)
            .filter(|&n| g.is_live(n))
            .map(|n| g.nodes()[n as usize].cov)
            .collect();
        assert_eq!(covs, vec![3, 1], "each copy carries its own reads");
    }

    /// Only a node the reads actually thread twice may be split. An ambiguous repeat —
    /// both entries supporting both exits — must be left alone, because guessing there
    /// is exactly how a misjoin is made.
    #[test]
    fn unzip_refuses_an_ambiguous_repeat() {
        use crate::graph::NodeAnno;
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(4, 80);
        }
        for (a, b) in [(0, 2), (1, 2), (2, 3), (2, 4)] {
            g.add_edge(edge(a, b));
        }
        let mut tr = Transitions::default();
        for t in [(0, 2, 3), (0, 2, 4), (1, 2, 3), (1, 2, 4)] {
            tr.insert(t, 5);
        }
        g.set_transitions(tr);
        let mut annos: Vec<NodeAnno> = (0..5)
            .map(|_| NodeAnno {
                anchor_read: 0,
                window: 0,
                ev_beg: 0,
                ev_end: 80,
            })
            .collect();
        let st = unzip(&mut g, &mut annos, 1);
        assert_eq!(st.nodes, 0, "ambiguity is not resolved by splitting");
        assert_eq!(g.n_nodes(), 5, "no copies made");
    }

    /// Node paths of the live graph, as node-id lists — a stand-in for `gen_unitigs`,
    /// which lives in `layout` and would be a circular dependency here.
    fn gen_unitigs_for_test(g: &FuzzyGraph) -> Vec<Vec<u32>> {
        let adj = ThreadAdj::new(g);
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::new();
        for n in 0..g.n_nodes() as u32 {
            if !g.is_live(n) || seen.contains(&n) {
                continue;
            }
            // only start at a node with no live predecessor
            if adj.pred.get(&n).is_some_and(|p| !p.is_empty()) {
                continue;
            }
            let mut path = vec![n];
            seen.insert(n);
            let mut cur = n;
            while let Some(succ) = adj.succ.get(&cur) {
                if succ.len() != 1 || seen.contains(&succ[0]) {
                    break;
                }
                cur = succ[0];
                seen.insert(cur);
                path.push(cur);
            }
            out.push(path);
        }
        out
    }

    /// `threaded_step` must refuse whenever the reads do not single out ONE
    /// continuation, in either direction — otherwise the winner is iteration order.
    #[test]
    fn threaded_step_requires_mutual_uniqueness() {
        let mut g = FuzzyGraph::new();
        for _ in 0..5 {
            g.add_node(4, 80);
        }
        for (a, b) in [(0, 2), (1, 2), (2, 3), (2, 4)] {
            g.add_edge(edge(a, b));
        }
        let mut tr = Transitions::default();
        tr.insert((0, 2, 3), 2);
        tr.insert((1, 2, 4), 2);
        let adj = ThreadAdj::new(&g);
        assert_eq!(threaded_step(&adj, &tr, 0, 2, 1), Some(3));
        assert_eq!(threaded_step(&adj, &tr, 1, 2, 1), Some(4));

        // Both entries thread into 3: ambiguous from 3's side, so neither may take it.
        let mut tr2 = Transitions::default();
        tr2.insert((0, 2, 3), 1);
        tr2.insert((1, 2, 3), 1);
        assert_eq!(
            threaded_step(&adj, &tr2, 0, 2, 1),
            None,
            "another predecessor claims 3"
        );

        // One entry threading to two successors is ambiguous from its own side.
        let mut tr3 = Transitions::default();
        tr3.insert((0, 2, 3), 1);
        tr3.insert((0, 2, 4), 1);
        assert_eq!(
            threaded_step(&adj, &tr3, 0, 2, 1),
            None,
            "two continuations"
        );

        // No triple at all: stop, do not guess.
        assert_eq!(threaded_step(&adj, &Transitions::default(), 0, 2, 1), None);

        // One read is not corroboration: `min_reads = 2` refuses a singly-threaded
        // triple, which is the whole point of counting rather than flagging.
        let mut tr4 = Transitions::default();
        tr4.insert((0, 2, 3), 1);
        assert_eq!(threaded_step(&adj, &tr4, 0, 2, 1), Some(3));
        assert_eq!(threaded_step(&adj, &tr4, 0, 2, 2), None);
    }
}
