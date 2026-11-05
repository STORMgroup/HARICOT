//! Graph simplification passes (de-tangling).
//!
//! Ports of wtdbg2's cleaning pipeline (`extern/wtdbg2/wtdbg-graph.h`), in wtdbg2's order:
//!   1. low-coverage edge filter
//!   2. transitive reduction (`myers_transitive_reduction_graph`, line 1877)
//!   3. tip trimming                (`trim_tips_graph`,            line 2381)
//!   4. bubble popping              (`pop_bubbles_graph`,          line 2763): the
//!      per-end two-arm popper [`pop_bubbles`], or with `branched_bubbles` wtdbg2's
//!      branched popper [`pop_branched_bubbles`] (`haricot assemble` uses 2: its unsafe
//!      mode too, alternated with tip trimming to convergence)
//!
//! [`clean`] also runs [`merge_parallel_nodes`] each round, which has no wtdbg2
//! counterpart — see its doc comment — and the relative edge gate
//! [`filter_weak_edges`] (`min_edge_support`, off by default).
//!
//! **These passes are orientation-aware** ([`crate::graph::ends::EndAdj`]). A node
//! has two ends and every traversal here — the transitive-reduction path search,
//! the tip walk, the bubble branch trace — enters a node by one end and leaves by
//! the *other*.
//!
//! An undirected view is lossy here: undirected degree 2 also matches a
//! **convergence** (both edges on one end, nothing on the other), and walking
//! through one traces a path that does not exist. See `ends.rs` for the argument.
//!
//! [`n_components`] and the `deg(n) == 1` dead-end test still use the undirected
//! [`FuzzyGraph::live_adjacency`], which is correct: connectivity and "how many
//! edges does this node have at all" are genuinely undirected notions.
//!
//! Still simpler than wtdbg2's bidirected string graph: we have no reverse
//! complement, since the overlapper is forward-strand only. Node ids stay stable —
//! removal is the soft-delete [`FuzzyGraph::remove_node`] tombstone.

use std::collections::VecDeque;

use rustc_hash::{FxHashMap, FxHashSet};

use crate::graph::ends::{End, EndAdj, entry_end, exit_end, from_end, to_end};
use crate::graph::fbg::{Dir, FuzzyGraph};
use crate::params::GraphParams;

/// Drop edges below `min_edge_cov`.
pub fn filter_low_coverage_edges(graph: &mut FuzzyGraph, params: &GraphParams) -> usize {
    let before = graph.n_edges();
    graph.edges_mut().retain(|e| e.cov >= params.min_edge_cov);
    before - graph.n_edges()
}

/// Drop edges below `min_edge_support · min(cov(from), cov(to))` — the *relative*
/// counterpart of [`filter_low_coverage_edges`]. No-op when `min_edge_support` is
/// `0.0` (the default), so the cleaned graph is unchanged until a run opts in.
///
/// Judging an edge against the coverage of the nodes it joins is what makes the
/// test depth-aware: the same `cov = 2` is strong support between two 4x nodes and
/// almost none between two 40x nodes. Edges touching an already-deleted node are
/// left alone — `prune_dangling_edges` owns those.
pub fn filter_weak_edges(graph: &mut FuzzyGraph, params: &GraphParams) -> usize {
    if params.min_edge_support <= 0.0 {
        return 0;
    }
    let f = params.min_edge_support as f64;
    // `support` is the read count (see `Node::support`); only `threading::unzip` sets it
    // otherwise, to its copies' share.
    let sup: Vec<f64> = graph.nodes().iter().map(|n| n.support as f64).collect();
    let del: Vec<bool> = graph.nodes().iter().map(|n| n.del).collect();
    let before = graph.n_edges();
    graph.edges_mut().retain(|e| {
        let (u, v) = (e.from as usize, e.to as usize);
        if del[u] || del[v] {
            return true;
        }
        let local = sup[u].min(sup[v]);
        (e.support as f64) >= f * local
    });
    before - graph.n_edges()
}

/// Genomic step length implied by an edge: how far the `to` node sits past the
/// `from` node, in events. `offset` is the inter-segment gap (negative = overlap),
/// so the start-to-start step is `node_events + offset`, floored at 1.
#[inline]
fn step_len(node_events: u32, offset: i32) -> i64 {
    (node_events as i64 + offset as i64).max(1)
}

/// Transitive reduction (Myers, with offset fuzz).
///
/// Removes an edge `u→w` when an alternative path `u→…→w` of comparable genomic
/// length already exists (so the direct edge is a redundant chord — typically a
/// "skip" edge bridging a gap that other reads fill with a real chain). Follows
/// wtdbg2's intent (`myers_transitive_reduction_graph`, `wtdbg-graph.h:1877`): an
/// edge implied by a longer two-/multi-hop path is reducible within an offset
/// fuzz. The search is length-budgeted and capped at `MAX_HOPS`, because a single
/// missing node can make the implied path several hops long.
///
/// `fuzz` tolerance: the source node's length (one node of slack) plus half the
/// edge's length, to absorb median-offset noise that grows with distance.
///
/// Node lengths are read from the graph rather than from a params field: a node's
/// length is `node_bins * bin_events` (see [`GraphParams::node_bins`]), and
/// `GraphParams` alone cannot know `bin_events`. Taking it from the node also keeps
/// the step length per-edge.
///
/// **Orientation-aware** ([`EndAdj`], see [`has_alt_path`]): the alternate path
/// must leave `u` by the end the chord leaves and arrive at `w` by the end the
/// chord arrives, traversing every intermediate node end-to-end. An undirected
/// search can justify deleting a real edge with a path that turns around inside a
/// node.
///
/// `step` is measured from the node being *left* — using `edges[ei].from`'s length
/// unconditionally is wrong for an edge traversed backwards.
pub fn transitive_reduction(graph: &mut FuzzyGraph, params: &GraphParams) -> usize {
    const MAX_HOPS: usize = 8;
    let _ = params;
    // The `removed` mask is consulted by the search, so each decision is made
    // against the graph as already reduced. That keeps the pass connectivity-
    // preserving: an edge is dropped only while an alternate path is still live, so
    // the last link between two nodes survives. Deciding all edges against the
    // unreduced graph instead can drop two mutually-supporting chords at once and
    // disconnect it.
    let ends = EndAdj::new(graph);
    let edges = graph.edges().to_vec();
    let node_len: Vec<u32> = graph.nodes().iter().map(|n| n.len_events).collect();
    // How far the far node sits past `leaving`, in events.
    let step = |ei: usize, leaving: u32| step_len(node_len[leaving as usize], edges[ei].offset);
    let elen = |ei: usize| step(ei, edges[ei].from);
    let mut removed: Vec<bool> = vec![false; edges.len()];

    // Process longest edges first: long chords are the redundant skip-edges; keeping
    // the short chain edges is what we want.
    let mut order: Vec<usize> = (0..edges.len())
        .filter(|&ei| {
            !graph.nodes()[edges[ei].from as usize].del && !graph.nodes()[edges[ei].to as usize].del
        })
        .collect();
    order.sort_by_key(|&ei| std::cmp::Reverse(elen(ei)));

    let mut n_removed = 0usize;
    for ei in order {
        let (u, w) = (edges[ei].from, edges[ei].to);
        let (ue, we) = (from_end(edges[ei].from_dir), to_end(edges[ei].to_dir));
        let target_len = elen(ei);
        let fuzz = node_len[u as usize] as i64 + target_len / 2;
        if has_alt_path(
            &ends, &step, &removed, u, ue, w, we, ei, target_len, fuzz, MAX_HOPS,
        ) {
            removed[ei] = true;
            n_removed += 1;
        }
    }

    let mut idx = 0usize;
    graph.edges_mut().retain(|_| {
        let keep = !removed[idx];
        idx += 1;
        keep
    });
    n_removed
}

/// Bounded Dijkstra: is there an alternative path from `u`'s `ue` end to `w`'s `we`
/// end (≥ 2 hops, excluding edge `skip_ei`) whose shortest length lands in
/// `[target - fuzz, target + fuzz]`?
///
/// Using shortest distances (with a visited map) keeps the search polynomial — the
/// naive all-paths DFS blows up in dense regions. The bridging chords we target
/// have an alternate fill-chain spanning the *same* genomic distance, so their
/// shortest alternate ≈ `target`; a genuine chain has no comparable alternate. We
/// require the alt ≥ `target − fuzz` so a far-shorter coincidental path (a real
/// short link) does not trigger removal of a longer real edge.
///
/// **Orientation-aware.** The search state is `(node, entry_end)`: a path arriving
/// at a node by one end must leave by the **other**, and the search starts only from
/// the edges on the end the chord actually leaves. Walking node-to-node instead would
/// enter and leave a node by the *same* end — a path that does not physically exist —
/// and delete a real edge as a "redundant chord" on the strength of it.
#[allow(clippy::too_many_arguments)]
fn has_alt_path<F: Fn(usize, u32) -> i64>(
    ends: &EndAdj,
    step: &F,
    removed: &[bool],
    u: u32,
    ue: End,
    w: u32,
    we: End,
    skip_ei: usize,
    target: i64,
    fuzz: i64,
    max_hops: usize,
) -> bool {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let budget = target + fuzz;
    // (Reverse(dist), node, entry_end, hops) min-heap.
    let mut heap: BinaryHeap<Reverse<(i64, u32, End, usize)>> = BinaryHeap::new();
    let mut best: FxHashMap<(u32, End), i64> = FxHashMap::default();
    // Leave u by the same end the chord does.
    for &(m, em, ei) in ends.at(u, ue) {
        if ei == skip_ei || removed[ei] || m == u {
            continue;
        }
        let d = step(ei, u);
        if d <= budget {
            heap.push(Reverse((d, m, em, 1)));
        }
    }
    while let Some(Reverse((dist, node, entry, hops))) = heap.pop() {
        if let Some(&b) = best.get(&(node, entry)) {
            if b <= dist {
                continue;
            }
        }
        best.insert((node, entry), dist);
        // Arriving at w by the chord's own end is what makes the chord redundant.
        if node == w && entry == we && hops >= 2 && dist >= target - fuzz {
            return true;
        }
        if hops >= max_hops {
            continue;
        }
        // A real path leaves by the opposite end from the one it arrived on.
        for &(k, ek, ei) in ends.at(node, entry.flip()) {
            if ei == skip_ei || removed[ei] || k == node {
                continue;
            }
            let nd = dist + step(ei, node);
            if nd <= budget && best.get(&(k, ek)).map(|&b| nd < b).unwrap_or(true) {
                heap.push(Reverse((nd, k, ek, hops + 1)));
            }
        }
    }
    false
}

/// Merge parallel duplicate nodes: group live nodes by their exact set of graph
/// neighbours and fold each group into one representative. Coverage is summed into
/// the representative; edges are redirected and de-duplicated. Because every edge of
/// a folded node is redirected onto the representative, the merge preserves
/// connectivity. Returns the number of nodes merged away.
///
/// Nodes sharing a neighbour set are parallel arms between identical flanks — the
/// signature expected of two tilings of one locus keyed on different anchor reads
/// (the per-anchor build lets each read pick its own window grid), which
/// transitive reduction and simple bubble popping leave behind as a bubble web.
/// Two *adjacent* nodes cannot be grouped: each lists the other while self is
/// excluded from the key, so their neighbour sets differ.
pub fn merge_parallel_nodes(graph: &mut FuzzyGraph) -> usize {
    let mut merged_total = 0usize;
    loop {
        let adj = graph.live_adjacency();
        // distinct-neighbour degree per node (used to tell bubble *arms* from the
        // bubble *endpoints*: both pairs share a neighbour set, but an arm's flanks
        // are junctions deg>=3 while an endpoint's "flanks" are the deg-2 arms).
        let distinct_deg = |n: u32| -> usize {
            let mut v: Vec<u32> = adj[n as usize]
                .iter()
                .map(|&(m, _)| m)
                .filter(|&m| m != n)
                .collect();
            v.sort_unstable();
            v.dedup();
            v.len()
        };
        // neighbour set (sorted, self excluded) -> live nodes with that exact set
        let mut groups: FxHashMap<Vec<u32>, Vec<u32>> = FxHashMap::default();
        for n in 0..graph.n_nodes() as u32 {
            if !graph.is_live(n) {
                continue;
            }
            let mut nbrs: Vec<u32> = adj[n as usize]
                .iter()
                .map(|&(m, _)| m)
                .filter(|&m| m != n)
                .collect();
            nbrs.sort_unstable();
            nbrs.dedup();
            if nbrs.len() < 2 {
                continue; // need >=2 shared flanks to be a confident parallel
            }
            groups.entry(nbrs).or_default().push(n);
        }
        // remap[node] = representative it folds into (default: itself)
        let mut remap: Vec<u32> = (0..graph.n_nodes() as u32).collect();
        let mut any = false;
        for (key, members) in groups {
            if members.len() < 2 {
                continue;
            }
            // Safety guard: only merge arms whose shared flanks are junctions
            // (deg>=3). In an isolated bubble the source/sink pair also shares a
            // neighbour set, but its "flanks" are the deg-2 arms, so this guard
            // skips it rather than folding the two endpoints into one node.
            if !key.iter().all(|&f| distinct_deg(f) >= 3) {
                continue;
            }
            let rep = *members
                .iter()
                .max_by_key(|&&m| graph.nodes()[m as usize].cov)
                .unwrap();
            let mut covsum = 0u64;
            for &m in &members {
                if m != rep {
                    remap[m as usize] = rep;
                    covsum += graph.nodes()[m as usize].cov as u64;
                    graph.remove_node(m);
                    merged_total += 1;
                    any = true;
                }
            }
            let rep_cov = &mut graph.nodes_mut()[rep as usize].cov;
            *rep_cov = (*rep_cov as u64 + covsum).min(u32::MAX as u64) as u32;
        }
        if !any {
            break;
        }
        // Redirect edges through the remap, then drop self-loops and merge parallels.
        for e in graph.edges_mut() {
            e.from = remap[e.from as usize];
            e.to = remap[e.to as usize];
        }
        dedup_edges(graph);
    }
    merged_total
}

/// Collapse parallel edges (same *oriented* endpoints, attaching to the same
/// *ends*) into one with summed coverage and the median offset, and drop
/// self-loops. Used after node merges.
///
/// Keyed on the full `(from, from_dir, to, to_dir)`: `from`/`to` carries the
/// read-order orientation the build stage recorded, so A->B and B->A are distinct
/// edges and must not be folded together (folding them would also median offsets of
/// opposite sign against each other). The **dirs** are part of the key, and are
/// preserved onto the surviving edge, because they are what says which physical END
/// each edge attaches to — the field [`EndAdj`] reads.
///
/// Keying on `(from, to)` alone and re-emitting as `Fwd`/`Fwd` would be a no-op only
/// while `build.rs` emits Fwd/Fwd exclusively (forward-strand overlapper); it would
/// silently rewrite every end assignment the moment cross-strand edges land, which is
/// precisely when it would be hardest to notice.
fn dedup_edges(graph: &mut FuzzyGraph) {
    use crate::graph::fbg::{Dir, Edge};
    let key_of = |d: Dir| match d {
        Dir::Fwd => 0u8,
        Dir::Rev => 1u8,
    };
    let dir_of = |k: u8| if k == 0 { Dir::Fwd } else { Dir::Rev };
    // (from, from_dir, to, to_dir) -> (cov, support, offsets)
    type Acc = FxHashMap<(u32, u8, u32, u8), (u32, f32, Vec<i32>)>;
    let mut acc: Acc = FxHashMap::default();
    for e in graph.edges() {
        if e.from == e.to {
            continue; // self-loop from a merge
        }
        let ent = acc
            .entry((e.from, key_of(e.from_dir), e.to, key_of(e.to_dir)))
            .or_insert((0, 0.0, Vec::new()));
        ent.0 += e.cov;
        // Support is additive for the same reason `cov` is: merging two edges pools
        // the reads asserting them.
        ent.1 += e.support;
        ent.2.push(e.offset);
    }
    let edges = graph.edges_mut();
    edges.clear();
    for ((from, fd, to, td), (cov, support, mut offs)) in acc {
        offs.sort_unstable();
        edges.push(Edge {
            from,
            from_dir: dir_of(fd),
            to,
            to_dir: dir_of(td),
            cov,
            support,
            offset: offs[offs.len() / 2],
        });
    }
}

/// Trim short dead-end branches (tips) up to `tip_steps` nodes.
///
/// Port of `trim_tips_graph` (`wtdbg-graph.h:2381`): a tip is a short chain that
/// dead-ends and hangs off a junction. We walk inward from each dead end along the
/// chain; if a junction is reached within `tip_steps`, the spur nodes (up to but
/// excluding the junction) are removed. A whole short *isolated* component with no
/// continuation is left intact — whether it terminates in a linear dead end or in a
/// convergence (two arms meeting with nothing beyond). Iterated to convergence
/// (trimming exposes new tips).
///
/// **Orientation-aware** ([`EndAdj`]): the walk steps with `step_out`, and only
/// continues when the chain is unbranched on *both* ends it passes through.
/// Undirected degree — "d == 2 ⇒ step to the neighbour I did not come from" — would
/// walk straight through a *convergence* (two edges on one end, nothing on the
/// other), traversing a path that does not exist; the nodes it collected would then
/// be deleted as a "tip". A dead end is still an undirected notion (one edge total),
/// so the `start` test uses plain degree.
pub fn trim_tips(graph: &mut FuzzyGraph, params: &GraphParams) -> usize {
    let mut removed = 0usize;
    loop {
        let adj = graph.live_adjacency();
        let ends = EndAdj::new(graph);
        let deg = |n: u32| adj[n as usize].len();
        let mut to_remove: FxHashSet<u32> = FxHashSet::default();

        for start in 0..graph.n_nodes() as u32 {
            if !graph.is_live(start) || deg(start) != 1 {
                continue;
            }
            // Orient `start` so the walk heads *inward*: its single edge sits on
            // one end, and we leave by that end.
            let dir = if ends.deg(start, End::R) == 1 {
                Dir::Fwd
            } else {
                Dir::Rev
            };
            let mut path = vec![start];
            let mut cur = (start, dir);
            let mut hit_junction = false;
            while path.len() <= params.tip_steps as usize {
                // The walk only ever stands on a node whose exit end carries exactly
                // one edge (`start` is oriented that way; interiors have one per
                // end), so `step_out` always names the next node `m`. Use it rather
                // than `next`, which *also* demands m's entry end be unbranched and
                // so cannot name the node it stopped in front of — and we must
                // inspect that node to tell a junction from a terminus.
                let Some((m, dm, _)) = ends.step_out(cur.0, cur.1) else {
                    break; // no continuation at all ⇒ nothing to trim toward
                };
                if m == start {
                    break; // tiny cycle guard
                }
                // NOTHING BEYOND `m` ⇒ the component terminates here, so the whole
                // thing is an isolated contig, not a spur off a backbone — leave it.
                // wtdbg2 agrees (`trim_tip_core_graph`, wtdbg-graph.h:2312 returns 0
                // when the trace ends at a dead end).
                //
                // This must be judged by m's FAR end alone. `deg(L) + deg(R) <= 1`
                // catches only a *linear* dead end and misses a **convergence
                // terminus** (degL=2, degR=0: two arms meeting with nothing past
                // them), which sums to 2, and `is_chain_interior` is false for a
                // terminus exactly as for a junction — so testing either one first
                // erases the component.
                if ends.deg(m, exit_end(dm)) == 0 {
                    break; // hit_junction stays false ⇒ nothing removed
                }
                if !ends.is_chain_interior(m) {
                    hit_junction = true;
                    break; // junction — stop, and do not include it in the tip
                }
                path.push(m);
                cur = (m, dm);
            }
            if hit_junction && path.len() <= params.tip_steps as usize {
                for n in path {
                    to_remove.insert(n);
                }
            }
        }

        if to_remove.is_empty() {
            break;
        }
        for n in &to_remove {
            graph.remove_node(*n);
        }
        graph.prune_dangling_edges();
        removed += to_remove.len();
    }
    removed
}

/// Pop bubbles: collapse short alternative paths between the same two endpoints.
///
/// Port of `pop_bubbles_graph` (`wtdbg-graph.h:2763`), simplified to the common
/// **simple bubble**: a source `s` that branches, from which two internally
/// unbranched paths reconverge at the same node `t` within `bubble_steps`. The
/// lower-coverage side's interior nodes are removed (the surviving side keeps the
/// stronger consensus). Iterated to convergence.
///
/// **Orientation-aware** ([`EndAdj`]): a branch is traced with `step_out`, and a
/// node is interior only if it is a true chain interior (one edge on *each* end).
/// Undirected degree would walk through a convergence — two edges on one end — and
/// so trace paths that do not exist, then delete their interiors as the weaker arm
/// of a "bubble". A source must branch on a *single* end (`deg(end) >= 2`);
/// undirected degree ≥ 3 both over-counts (a node with one edge per end has nothing
/// to branch) and misses a pure convergence (deg 2 on one end, which IS a branch
/// point).
pub fn pop_bubbles(graph: &mut FuzzyGraph, params: &GraphParams) -> usize {
    let mut removed = 0usize;
    loop {
        let ends = EndAdj::new(graph);
        let mut popped_any = false;
        let mut to_remove: FxHashSet<u32> = FxHashSet::default();

        'src: for s in 0..graph.n_nodes() as u32 {
            if !graph.is_live(s) {
                continue;
            }
            for send in [End::L, End::R] {
                if ends.deg(s, send) < 2 {
                    continue; // s does not branch on this end
                }
                // Trace each branch leaving THIS end of s to its first
                // non-chain-interior node. Record exit -> (branch_interior,
                // cov_sum); two branches sharing an exit form a bubble.
                //
                // `exits` MUST be scoped per end. A bubble is two arms leaving the
                // SAME end of s and reconverging; arms leaving opposite ends are
                // just the chain continuing both ways, and if they happen to meet
                // (a cycle through s) sharing one map calls it a bubble and deletes
                // an arm's interior, breaking the cycle.
                //
                // It is keyed on `(exit node, entry end)`, not the bare node id, for
                // the mirror-image reason: two arms off the same end of s can arrive
                // at OPPOSITE ends of the same node t, which is a cycle through t,
                // not a bubble. A bubble's arms must re-enter t by the *same* end.
                let mut exits: FxHashMap<(u32, End), (Vec<u32>, u64)> = FxHashMap::default();
                for &(first, fem, _ei) in ends.at(s, send) {
                    let mut interior: Vec<u32> = Vec::new();
                    let mut cov_sum = 0u64;
                    let mut steps = 0u32;
                    let mut cur = (first, if fem == End::L { Dir::Fwd } else { Dir::Rev });
                    loop {
                        if to_remove.contains(&cur.0) {
                            break;
                        }
                        cov_sum += graph.nodes()[cur.0 as usize].cov as u64;
                        if !ends.is_chain_interior(cur.0) || cur.0 == s {
                            // exit node (junction or dead end) — close this branch.
                            if cur.0 != s && steps >= 1 {
                                // We arrive at `cur.0` oriented `cur.1`, so the end
                                // we re-enter by is `entry_end(cur.1)`.
                                let exit_key = (cur.0, entry_end(cur.1));
                                if let Some(prev_branch) = exits.get(&exit_key) {
                                    // Two branches from s meet at the same exit:
                                    // a bubble. Remove the weaker interior.
                                    let (ref pa, ca) = *prev_branch;
                                    let (this_int, cb) = (interior.clone(), cov_sum);
                                    let weaker = if ca <= cb { pa.clone() } else { this_int };
                                    if !weaker.is_empty() {
                                        for n in weaker {
                                            to_remove.insert(n);
                                        }
                                        popped_any = true;
                                        continue 'src;
                                    }
                                } else {
                                    exits.insert(exit_key, (interior.clone(), cov_sum));
                                }
                            }
                            break;
                        }
                        interior.push(cur.0);
                        steps += 1;
                        if steps > params.bubble_steps {
                            break;
                        }
                        match ends.step_out(cur.0, cur.1) {
                            Some((m, dm, _)) => cur = (m, dm),
                            None => break,
                        }
                    }
                }
            }
        }

        if to_remove.is_empty() || !popped_any {
            break;
        }
        for n in &to_remove {
            graph.remove_node(*n);
        }
        graph.prune_dangling_edges();
        removed += to_remove.len();
    }
    removed
}

/// One search entry of [`branched_bubble_at`]: wtdbg2's `bt_t` (`wtdbg-graph.h:2498`), an
/// arrival at `node` by `edge`, `step` edges from the source.
#[derive(Debug, Clone, Copy)]
struct Bt {
    node: u32,
    dir: Dir,
    /// The edge this arrival came by; `None` for the source (and the index-0 sentinel).
    edge: Option<usize>,
    step: u32,
    /// The entry this one was expanded from; 0 = none.
    parent: usize,
    /// Sum of edge coverage along the arrival path (heap tie-break only).
    score: i64,
    /// Unsafe mode: reached from the source along a simple chain.
    ind: bool,
    /// The node's other entry this arrival met; 0 = first arrival.
    ending: usize,
}

/// wtdbg2's cap on one search's entries (`WT_MAX_BTIDX`, 16-bit back-pointers).
const MAX_BT: usize = 0xFFFF;

/// Search one bubble from source `s`, left in direction `sdir`: a port of
/// `pop_bubble_core_graph` (`wtdbg-graph.h:2665`). Returns the edges to cut and the nodes
/// of the region, or `None` if the region is not a bubble within `max_step` edges.
///
/// Breadth-first by step from `s`. A node is CLOSED once every edge on its entry end has
/// been traversed from inside the search, and only closed nodes are expanded; so if any
/// node is also entered from outside, it never closes and nothing is popped. When the
/// frontier is a single closed node `t` and no node is open, every path from `s` that
/// does not dead-end reaches `t`: the region is a bubble of any shape. Of its paths the
/// one kept is the LONGEST in steps (at each node, the latest-closing arrival): a skip
/// edge over a node is the shorter path, and the fuller tiling carries the sequence. The
/// other arms are cut where they rejoin the kept path, back to where they left it.
/// Arms that dead-end inside are left attached (tips; trimming handles them).
///
/// Unsafe mode (`safe = false`) adds wtdbg2's early exit: two arrivals at one node both
/// on simple chains from `s` cut one chain at once, the shorter unless tied.
///
/// DEVIATION: an arm is cut when its NODE's kept arrival is on the kept path; wtdbg2
/// tests the one entry the arrival met (`ending`), which misses a third arrival at a node
/// and leaves that arm for a later pass.
fn branched_bubble_at(
    graph: &FuzzyGraph,
    ends: &EndAdj,
    s: u32,
    sdir: Dir,
    max_step: u32,
    safe: bool,
    touched: &FxHashSet<u32>,
) -> Option<(Vec<usize>, Vec<u32>)> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    if ends.deg(s, exit_end(sdir)) < 2 {
        return None;
    }
    let entering = |em: End| {
        if em == entry_end(Dir::Fwd) {
            Dir::Fwd
        } else {
            Dir::Rev
        }
    };
    let sentinel = Bt {
        node: s,
        dir: sdir,
        edge: None,
        step: 0,
        parent: 0,
        score: 0,
        ind: false,
        ending: 0,
    };
    let mut bts = vec![
        sentinel,
        Bt {
            ind: !safe,
            ..sentinel
        },
    ];
    // node -> (its kept arrival, entry edges not yet traversed, single entry edge)
    let mut state: FxHashMap<u32, (usize, usize, bool)> = FxHashMap::default();
    state.insert(s, (1, 0, true));
    // wtdbg2's heap order: fewest steps first, then the highest score.
    let mut heap = BinaryHeap::new();
    heap.push(Reverse((0u32, 0i64, 1usize)));
    let mut unclosed = 0usize;
    while let Some(Reverse((_, _, b))) = heap.pop() {
        if bts.len() >= MAX_BT || bts[b].step >= max_step || touched.contains(&bts[b].node) {
            return None;
        }
        if bts[b].ind && !state[&bts[b].node].2 {
            bts[b].ind = false;
        }
        let cur = bts[b];
        let lst = bts.len();
        for &(m, em, ei) in ends.at(cur.node, exit_end(cur.dir)) {
            if m == s {
                return None; // back to the source: a cycle
            }
            bts.push(Bt {
                node: m,
                dir: entering(em),
                edge: Some(ei),
                step: cur.step + 1,
                parent: b,
                score: cur.score + graph.edges()[ei].cov as i64,
                ind: false,
                ending: 0,
            });
        }
        if cur.ind && (cur.parent == 0 || lst + 1 == bts.len()) {
            for x in &mut bts[lst..] {
                x.ind = true;
            }
        }
        for i in lst..bts.len() {
            let n = bts[i].node;
            if touched.contains(&n) {
                return None;
            }
            match state.get(&n).copied() {
                None => {
                    let unv = ends.deg(n, entry_end(bts[i].dir));
                    state.insert(n, (i, unv, unv == 1));
                    unclosed += 1;
                }
                Some((f, _, _)) => {
                    bts[i].ending = f;
                    if bts[i].dir != bts[f].dir {
                        return None; // re-entered by the other end: a cycle
                    }
                    if !safe {
                        let (ti, fi) = (bts[i].ind, bts[f].ind);
                        let pair = match (ti, fi) {
                            (true, true) if bts[i].step == bts[f].step => Some((i, f)),
                            (true, true) => Some((f, i)),
                            (true, false) => Some((i, f)),
                            (false, true) => Some((f, i)),
                            _ => None,
                        };
                        if let Some((b1, b2)) = pair {
                            return Some(cut_redundant_chain(&bts, b1, b2));
                        }
                    }
                }
            }
            let st = state.get_mut(&n).expect("inserted above");
            st.1 = st.1.saturating_sub(1);
            if st.1 == 0 {
                let f = st.0;
                if bts[i].step > bts[f].step {
                    bts[f].ending = i;
                    st.0 = i;
                    bts[i].ending = 0;
                }
                let k = st.0;
                if ends.deg(n, exit_end(bts[k].dir)) > 0 {
                    heap.push(Reverse((bts[k].step, -bts[k].score, k)));
                }
                unclosed -= 1;
            }
        }
        if heap.len() == 1 && unclosed == 0 {
            let Reverse((_, _, t)) = *heap.peek().expect("len 1");
            return Some(cut_bubble_arms(&bts, &state, t));
        }
    }
    None
}

/// Unsafe-mode early exit (`safe_cut_redundant_edges_graph`, `wtdbg-graph.h:2618`): cut
/// `b1`'s arrival path back to where it meets `b2`'s.
fn cut_redundant_chain(bts: &[Bt], mut b1: usize, mut b2: usize) -> (Vec<usize>, Vec<u32>) {
    let mut cut = Vec::new();
    while b1 != b2 {
        if bts[b1].step >= bts[b2].step {
            cut.extend(bts[b1].edge);
            b1 = bts[b1].parent;
        } else {
            b2 = bts[b2].parent;
        }
    }
    (cut, bts[1..].iter().map(|b| b.node).collect())
}

/// Keep the arrival path of `t`; cut every other arm where it rejoins a kept node, back
/// to where it left the kept path (`pop_bubble2_backtrace_graph`, `wtdbg-graph.h:2543`).
fn cut_bubble_arms(
    bts: &[Bt],
    state: &FxHashMap<u32, (usize, usize, bool)>,
    t: usize,
) -> (Vec<usize>, Vec<u32>) {
    // 1 = kept path, 2 = undecided, 0 = cut
    let mut keep = vec![2u8; bts.len()];
    let mut idx = t;
    while idx != 0 {
        keep[idx] = 1;
        idx = bts[idx].parent;
    }
    for i in 1..bts.len() {
        if keep[i] == 1 || keep[state[&bts[i].node].0] != 1 {
            continue;
        }
        let mut idx = i;
        while idx != 0 && keep[idx] == 2 {
            keep[idx] = 0;
            idx = bts[idx].parent;
        }
    }
    let cut = (1..bts.len())
        .filter(|&i| keep[i] == 0)
        .filter_map(|i| bts[i].edge)
        .collect();
    (cut, bts[1..].iter().map(|b| b.node).collect())
}

/// Pop every bubble wtdbg2's popper finds ([`branched_bubble_at`] from each node end),
/// iterated to convergence. Regions popped in one pass do not overlap (a search touching
/// an already-popped region is abandoned until the next pass). After the cuts, the
/// popped arms left without a connection to anything outside them are deleted -- wtdbg2
/// deletes them as isolated nodes -- and arms still hanging off the kept path stay for
/// tip trimming. Returns the number of bubbles popped.
pub fn pop_branched_bubbles(graph: &mut FuzzyGraph, params: &GraphParams, safe: bool) -> usize {
    let mut popped = 0usize;
    loop {
        let ends = EndAdj::new(graph);
        let mut touched: FxHashSet<u32> = FxHashSet::default();
        let mut cut: FxHashSet<usize> = FxHashSet::default();
        let mut regions: Vec<Vec<u32>> = Vec::new();
        for s in 0..graph.n_nodes() as u32 {
            for d in [Dir::Fwd, Dir::Rev] {
                if !graph.is_live(s) || touched.contains(&s) {
                    break;
                }
                let Some((edges, nodes)) =
                    branched_bubble_at(graph, &ends, s, d, params.bubble_steps, safe, &touched)
                else {
                    continue;
                };
                if edges.is_empty() {
                    continue;
                }
                cut.extend(edges);
                touched.insert(s);
                touched.extend(nodes.iter().copied());
                regions.push(nodes);
                popped += 1;
            }
        }
        if cut.is_empty() {
            break;
        }
        let mut i = 0usize;
        graph.edges_mut().retain(|_| {
            let keep = !cut.contains(&i);
            i += 1;
            keep
        });
        // Delete a region's arms that no longer reach outside the region's cut-off part.
        let adj = graph.live_adjacency();
        for nodes in &regions {
            let region: FxHashSet<u32> = nodes.iter().copied().collect();
            let mut done: FxHashSet<u32> = FxHashSet::default();
            for &n in nodes {
                if done.contains(&n) {
                    continue;
                }
                let mut comp = vec![n];
                let mut q = VecDeque::from([n]);
                done.insert(n);
                let mut escapes = false;
                while let Some(x) = q.pop_front() {
                    for &(y, _) in &adj[x as usize] {
                        if !region.contains(&y) {
                            escapes = true;
                        } else if done.insert(y) {
                            comp.push(y);
                            q.push_back(y);
                        }
                    }
                }
                if !escapes {
                    for x in comp {
                        graph.remove_node(x);
                    }
                }
            }
        }
        graph.prune_dangling_edges();
    }
    popped
}

/// Run the full cleaning pipeline in wtdbg2's order. Returns
/// `(edges_filtered, edges_reduced, tip_nodes, bubble_nodes)`.
pub fn clean(graph: &mut FuzzyGraph, params: &GraphParams) -> (usize, usize, usize, usize) {
    let filtered = filter_low_coverage_edges(graph, params);
    // wtdbg2 alternates reduction with tip/bubble removal: each removal can expose a
    // new chord, so the passes below are iterated (4 rounds). Every pass is
    // connectivity-preserving
    // (transitive reduction only drops edges with a live alternate path; tip/bubble
    // removal only drops dead-end spurs and degree-2 bubble interiors).
    let mut reduced = 0;
    let mut tips = 0;
    let mut bubbles = 0;
    let mut merged = 0;
    let mut weak = 0;
    for _ in 0..4 {
        // Collapse parallel double-tilings first so reduction/tip/bubble work on a
        // de-duplicated graph, then again each round as new parallels are exposed.
        merged += merge_parallel_nodes(graph);
        // Relative edge gate ahead of reduction: it changes which chords exist, and
        // it is re-run each round because tip/bubble removal changes which nodes are
        // live (and `merge_parallel_nodes` rewrites the representative's coverage).
        weak += filter_weak_edges(graph, params);
        reduced += transitive_reduction(graph, params);
        tips += trim_tips(graph, params);
        bubbles += match params.branched_bubbles {
            0 => pop_bubbles(graph, params),
            _ => pop_branched_bubbles(graph, params, true),
        };
    }
    // wtdbg2 turns its popper's unsafe mode on once the safe passes find nothing more,
    // alternating it with tip trimming to convergence (`wtdbg.c:1080-1110`).
    if params.branched_bubbles >= 2 {
        loop {
            let t = trim_tips(graph, params);
            let b = pop_branched_bubbles(graph, params, false);
            tips += t;
            bubbles += b;
            if t + b == 0 {
                break;
            }
        }
    }
    graph.prune_dangling_edges();
    log::info!(
        "clean: merge_parallel removed {merged} nodes, transitive_reduction removed {reduced} \
         edges, edge-cov filter removed {filtered}, relative edge-support filter removed {weak} \
         edges, trim_tips removed {tips} nodes, pop_bubbles removed {bubbles} (nodes; bubbles if branched_bubbles={bb})",
        bb = params.branched_bubbles
    );
    (filtered, reduced, tips, bubbles)
}

/// Number of connected components over live nodes (helper for tests/diagnostics).
pub fn n_components(graph: &FuzzyGraph) -> usize {
    let adj = graph.live_adjacency();
    let mut seen: FxHashSet<u32> = FxHashSet::default();
    let mut comps = 0usize;
    for s in 0..graph.n_nodes() as u32 {
        if !graph.is_live(s) || seen.contains(&s) {
            continue;
        }
        comps += 1;
        let mut q = VecDeque::new();
        q.push_back(s);
        seen.insert(s);
        while let Some(n) = q.pop_front() {
            for &(nbr, _) in &adj[n as usize] {
                if seen.insert(nbr) {
                    q.push_back(nbr);
                }
            }
        }
    }
    comps
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::fbg::{Dir, Edge, FuzzyGraph};
    use crate::params::Params;

    fn e(from: u32, to: u32, cov: u32, offset: i32) -> Edge {
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

    /// A braid p -> s -> {a, b} -> {c, d} -> t -> q (every a/b links to every c/d): no arm
    /// is a simple chain, so the simple-bubble port cannot pop it; wtdbg2's popper keeps one
    /// path through, and the arms cut off from everything are deleted.
    #[test]
    fn branched_bubble_braid_pops_to_one_chain() {
        let p = Params::r10_4_1().graph;
        let braid = || {
            let mut g = FuzzyGraph::new();
            let n: Vec<u32> = (0..8).map(|_| g.add_node(3, 72)).collect();
            let [pp, s, a, b, c, d, t, q] = n[..] else {
                unreachable!()
            };
            for (x, y) in [
                (pp, s),
                (s, a),
                (s, b),
                (a, c),
                (b, c),
                (a, d),
                (b, d),
                (c, t),
                (d, t),
                (t, q),
            ] {
                g.add_edge(e(x, y, 2, 0));
            }
            (g, n)
        };
        let (mut g, _) = braid();
        assert_eq!(pop_bubbles(&mut g, &p), 0, "the simple port sees no bubble");
        let (mut g, n) = braid();
        assert_eq!(pop_branched_bubbles(&mut g, &p, true), 1);
        let live: Vec<u32> = n.iter().copied().filter(|&x| g.is_live(x)).collect();
        assert_eq!(
            live,
            vec![n[0], n[1], n[2], n[4], n[6], n[7]],
            "p s a c t q"
        );
        assert_eq!(crate::layout::gen_unitigs(&g).len(), 1);
    }

    /// A node of the region also entered from outside never closes: no bubble. Here `a` is
    /// entered from x (blocks the search from s) and left to z (blocks the search back from
    /// t). Without them, a skip edge over one node is the shorter path: the tiling is kept.
    #[test]
    fn branched_bubble_needs_a_closed_region_and_keeps_the_longer_path() {
        let p = Params::r10_4_1().graph;
        let mut g = FuzzyGraph::new();
        let n: Vec<u32> = (0..8).map(|_| g.add_node(3, 72)).collect();
        let [pp, s, a, t, q, x, z, zz] = n[..] else {
            unreachable!()
        };
        for (u, v) in [
            (pp, s),
            (s, a),
            (a, t),
            (s, t),
            (t, q),
            (x, a),
            (a, z),
            (z, zz),
        ] {
            g.add_edge(e(u, v, 2, 0));
        }
        assert_eq!(
            pop_branched_bubbles(&mut g, &p, true),
            0,
            "a is entered from x"
        );
        for r in [x, z, zz] {
            g.remove_node(r);
        }
        g.prune_dangling_edges();
        assert_eq!(pop_branched_bubbles(&mut g, &p, true), 1);
        assert!(g.is_live(a));
        assert!(
            g.edges().iter().all(|y| !(y.from == s && y.to == t)),
            "skip edge cut"
        );
        assert_eq!(g.n_edges(), 4);
    }

    /// Unsafe mode: two simple chains from s meeting at t, where t is also entered from y
    /// and s also left to z (so neither search direction closes a region), still lose the
    /// shorter chain.
    #[test]
    fn branched_bubble_unsafe_cuts_the_shorter_simple_chain() {
        let p = Params::r10_4_1().graph;
        let mut g = FuzzyGraph::new();
        let n: Vec<u32> = (0..8).map(|_| g.add_node(3, 72)).collect();
        let [pp, s, a, b, c, t, y, z] = n[..] else {
            unreachable!()
        };
        for (u, v) in [
            (pp, s),
            (s, a),
            (a, t),
            (s, b),
            (b, c),
            (c, t),
            (y, t),
            (s, z),
        ] {
            g.add_edge(e(u, v, 2, 0));
        }
        assert_eq!(pop_branched_bubbles(&mut g, &p, true), 0);
        assert_eq!(pop_branched_bubbles(&mut g, &p, false), 1);
        assert!(
            !g.is_live(a),
            "the shorter chain s -> a -> t is cut off and deleted"
        );
        assert!(g.is_live(b) && g.is_live(c) && g.is_live(y));
    }

    /// Weighting is what makes the relative gate usable at low depth. Two nodes at
    /// cov 4 with one weak edge: by COUNT the edge is 1 of 4 and survives `f = 0.20`;
    /// by STRENGTH it is one z=5 assertion against nodes built from z=50 ones, and it
    /// dies. That distinction is the whole point -- at the dmel arm's mean node
    /// coverage of 3.80 a count has almost nothing to work with.
    #[test]
    fn weighted_support_kills_an_edge_a_count_keeps() {
        let mut p = Params::r10_4_1().graph;
        p.min_edge_support = 0.2;
        let mut g = FuzzyGraph::new();
        let a = g.add_node(4, 80);
        let b = g.add_node(4, 80);
        let mut e = e(a, b, 1, 0);

        // Unweighted: support == cov, 1 >= 0.2 * 4, kept.
        g.add_edge(e.clone());
        assert_eq!(filter_weak_edges(&mut g, &p), 0, "count keeps it");

        // Weighted: the nodes were built from strong overlaps, the edge from one weak
        // one. Same topology, same threshold, opposite verdict.
        g.set_node_support(a, 200.0);
        g.set_node_support(b, 200.0);
        e.support = 5.0;
        let edges = g.edges_mut();
        edges.clear();
        edges.push(e);
        assert_eq!(filter_weak_edges(&mut g, &p), 1, "strength drops it");
    }

    /// The relative gate is depth-aware where the absolute one is not: the *same*
    /// `cov = 3` edge survives between two 8x nodes and dies between two 40x nodes.
    /// That is the whole point, so it is what the test pins.
    #[test]
    fn filter_weak_edges_is_depth_relative() {
        let mut p = Params::r10_4_1().graph;
        p.min_edge_support = 0.2;

        // Shallow locus: 3 / min(8, 8) = 0.375 >= 0.2, keep.
        let mut g = FuzzyGraph::new();
        g.add_node(8, 80);
        g.add_node(8, 80);
        g.add_edge(e(0, 1, 3, 40));
        assert_eq!(
            filter_weak_edges(&mut g, &p),
            0,
            "3 of 8x is strong support"
        );

        // Deep locus, identical edge: 3 / min(40, 40) = 0.075 < 0.2, drop.
        let mut g = FuzzyGraph::new();
        g.add_node(40, 80);
        g.add_node(40, 80);
        g.add_edge(e(0, 1, 3, 40));
        assert_eq!(filter_weak_edges(&mut g, &p), 1, "3 of 40x is a spur");

        // Gated on the SHALLOWER endpoint, so a deep node cannot condemn an edge
        // into a genuinely shallow one: 3 / min(40, 8) = 0.375 >= 0.2, keep.
        let mut g = FuzzyGraph::new();
        g.add_node(40, 80);
        g.add_node(8, 80);
        g.add_edge(e(0, 1, 3, 40));
        assert_eq!(
            filter_weak_edges(&mut g, &p),
            0,
            "min(), not max() or mean()"
        );
    }

    /// Off by default: the cleaned graph must be untouched until a run opts in.
    #[test]
    fn filter_weak_edges_is_off_by_default() {
        let p = Params::r10_4_1().graph;
        assert_eq!(p.min_edge_support, 0.0);
        let mut g = FuzzyGraph::new();
        g.add_node(40, 80);
        g.add_node(40, 80);
        g.add_edge(e(0, 1, 1, 40)); // 1 of 40x — would die at any positive f
        assert_eq!(filter_weak_edges(&mut g, &p), 0);
    }

    /// `dedup_edges` must key on the *oriented* `(from, to)` pair. Keying on the
    /// unordered pair and negating `offset` to match would (a) fold the genuinely
    /// distinct 2->1 and 1->2 into one edge and (b) median opposite-signed offsets
    /// against each other. Here ids are set by hand so the high->low edge is the one
    /// that would get flipped.
    #[test]
    fn dedup_edges_keeps_opposite_orientations_distinct() {
        let mut g = FuzzyGraph::new();
        for _ in 0..3 {
            g.add_node(5, 80);
        }
        g.add_edge(e(2, 1, 3, 40)); // high -> low: the flip candidate
        g.add_edge(e(1, 2, 2, 40)); // the genuine opposite orientation
        g.add_edge(e(1, 2, 4, 60)); // a true parallel edge, must fold into the above
        dedup_edges(&mut g);

        let mut seen: Vec<(u32, u32, u32, i32)> = g
            .edges()
            .iter()
            .map(|x| (x.from, x.to, x.cov, x.offset))
            .collect();
        seen.sort_unstable();
        assert_eq!(
            seen,
            vec![(1, 2, 6, 60), (2, 1, 3, 40)],
            "2->1 and 1->2 must stay distinct; only the two 1->2 edges may fold \
             (cov 2+4=6), and no offset may be sign-flipped"
        );
    }

    #[test]
    fn low_coverage_edges_are_filtered() {
        let mut g = FuzzyGraph::new();
        let a = g.add_node(5, 80);
        let b = g.add_node(5, 80);
        g.add_edge(e(a, b, 1, 0));
        g.add_edge(e(a, b, 3, 0));

        let mut params = Params::r10_4_1().graph;
        params.min_edge_cov = 2;
        let removed = filter_low_coverage_edges(&mut g, &params);
        assert_eq!(removed, 1);
        assert_eq!(g.n_edges(), 1);
    }

    /// A chord 0—2 over the chain 0—1—2 (equal genomic lengths) is reduced.
    #[test]
    fn transitive_chord_is_reduced() {
        let mut g = FuzzyGraph::new();
        for _ in 0..3 {
            g.add_node(3, 72);
        }
        // node_events=72, offset 0 → step 72; chord 0->2 has offset 72 → step 144.
        g.add_edge(e(0, 1, 3, 0));
        g.add_edge(e(1, 2, 3, 0));
        g.add_edge(e(0, 2, 2, 72)); // skip edge: length 144 ≈ 72+72
        // Node length (72) comes from the nodes themselves, not from params.
        let p = Params::r10_4_1().graph;
        let removed = transitive_reduction(&mut g, &p);
        assert_eq!(removed, 1, "the chord 0-2 should be reduced");
        // The chain edges survive.
        assert_eq!(g.n_edges(), 2);
        assert_eq!(super::n_components(&g), 1);
    }

    /// A genuine chain (no chord) is left untouched.
    #[test]
    fn plain_chain_is_untouched() {
        let mut g = FuzzyGraph::new();
        for _ in 0..4 {
            g.add_node(3, 72);
        }
        g.add_edge(e(0, 1, 3, 0));
        g.add_edge(e(1, 2, 3, 0));
        g.add_edge(e(2, 3, 3, 0));
        let p = Params::r10_4_1().graph;
        assert_eq!(transitive_reduction(&mut g, &p), 0);
    }

    /// A short dead-end spur off a junction is trimmed; the long backbone arms
    /// (longer than `tip_steps`) survive.
    #[test]
    fn tip_is_trimmed() {
        let mut g = FuzzyGraph::new();
        for _ in 0..8 {
            g.add_node(3, 72);
        }
        // Backbone 0-1-2-3-4-5-6 (junction at 3), short spur 3-7.
        for i in 0..6 {
            g.add_edge(e(i, i + 1, 3, 0));
        }
        g.add_edge(e(3, 7, 1, 0));
        let mut p = Params::r10_4_1().graph;
        p.tip_steps = 2; // arms (len 3 to the junction) survive; spur (len 1) trimmed
        let removed = trim_tips(&mut g, &p);
        assert_eq!(removed, 1, "spur node 7 should be trimmed");
        assert!(!g.is_live(7));
        assert!(g.is_live(0) && g.is_live(6), "backbone arms survive");
        assert_eq!(super::n_components(&g), 1);
    }

    /// A simple bubble (two parallel interiors between s and t) collapses to one.
    #[test]
    fn bubble_is_popped() {
        let mut g = FuzzyGraph::new();
        // s=0, t=3, path A: 0-1-3 (cov hi), path B: 0-2-3 (cov lo).
        g.add_node(3, 72); // 0
        g.add_node(5, 72); // 1 (strong)
        g.add_node(1, 72); // 2 (weak)
        g.add_node(3, 72); // 3
        // give s and t degree >=3 so they read as junctions
        g.add_node(3, 72); // 4 extra neighbor of 0
        g.add_node(3, 72); // 5 extra neighbor of 3
        g.add_edge(e(0, 1, 3, 0));
        g.add_edge(e(1, 3, 3, 0));
        g.add_edge(e(0, 2, 1, 0));
        g.add_edge(e(2, 3, 1, 0));
        g.add_edge(e(0, 4, 3, 0));
        g.add_edge(e(3, 5, 3, 0));
        let mut p = Params::r10_4_1().graph;
        p.bubble_steps = 40;
        let removed = pop_bubbles(&mut g, &p);
        assert_eq!(removed, 1, "the weak interior node 2 should be removed");
        assert!(!g.is_live(2));
        assert!(g.is_live(1), "the strong interior survives");
    }

    /// Parallel arms of an embedded bubble merge into one node; the source/sink
    /// (whose flanks are the deg-2 arms) must not be merged.
    #[test]
    fn parallel_arms_merge_endpoints_survive() {
        let mut g = FuzzyGraph::new();
        // backbone u(0)-s(1)-{A(2),B(3)}-t(4)-v(5)
        for _ in 0..6 {
            g.add_node(4, 72);
        }
        g.add_edge(e(0, 1, 4, 0)); // u-s
        g.add_edge(e(1, 2, 4, 0)); // s-A
        g.add_edge(e(1, 3, 4, 0)); // s-B
        g.add_edge(e(2, 4, 4, 0)); // A-t
        g.add_edge(e(3, 4, 4, 0)); // B-t
        g.add_edge(e(4, 5, 4, 0)); // t-v
        let merged = merge_parallel_nodes(&mut g);
        assert_eq!(merged, 1, "arms A,B should fold into one node");
        assert!(g.is_live(1) && g.is_live(4), "source s and sink t survive");
        assert_eq!(super::n_components(&g), 1, "graph stays connected");
        // exactly one of A,B remains, as a degree-2 node between s and t.
        assert!(
            g.is_live(2) ^ g.is_live(3),
            "one arm survives, one is merged away"
        );
    }

    /// Safety: an isolated bubble (s,t have no backbone continuation, so every node
    /// is degree 2) is left untouched rather than collapsing its source/sink pair.
    #[test]
    fn isolated_bubble_is_not_merged() {
        let mut g = FuzzyGraph::new();
        for _ in 0..4 {
            g.add_node(4, 72); // s=0, A=1, B=2, t=3
        }
        g.add_edge(e(0, 1, 4, 0));
        g.add_edge(e(1, 3, 4, 0));
        g.add_edge(e(0, 2, 4, 0));
        g.add_edge(e(2, 3, 4, 0));
        let merged = merge_parallel_nodes(&mut g);
        assert_eq!(
            merged, 0,
            "no junction flanks (all deg 2) → conservatively skip"
        );
    }

    /// REGRESSION: a whole isolated linear component must survive trim_tips — it is
    /// not a tip, it is a (short) contig. The docstring has always promised this and
    /// HEAD's undirected code enforced it with an explicit `if d <= 1 { break; }`.
    /// wtdbg2 does the same (`trim_tip_core_graph`, wtdbg-graph.h:2312).
    ///
    /// The orientation-aware rewrite dropped that guard: `is_chain_interior` is
    /// `deg(L)==1 && deg(R)==1`, which is false for a DEAD END just as it is for a
    /// junction, and both collapsed into `hit_junction = true`. On 0-1-2 the walk
    /// from 0 flags a junction at dead-end 2 and removes {0,1}; the walk from 2
    /// removes {2,1}; the union erases the component. Any isolated chain of
    /// <= tip_steps+1 nodes vanished.
    ///
    /// `tip_is_trimmed` cannot catch this: its tip_steps=2 with length-3 arms means
    /// the walk exits on the loop bound before ever reaching a dead end.
    #[test]
    fn isolated_linear_component_is_not_trimmed() {
        let mut p = Params::r10_4_1().graph;
        p.tip_steps = 10; // the default
        let mut g = FuzzyGraph::new();
        for _ in 0..3 {
            g.add_node(5, 80);
        }
        g.add_edge(e(0, 1, 3, 0));
        g.add_edge(e(1, 2, 3, 0));
        let removed = trim_tips(&mut g, &p);
        assert_eq!(removed, 0, "an isolated chain is a contig, not a tip");
        assert!(
            g.is_live(0) && g.is_live(1) && g.is_live(2),
            "the whole component was erased"
        );
    }

    /// REGRESSION: an isolated component terminating in a **convergence** is not a
    /// tip either. `0R→2L` and `1R→2L` leave node 2 with degL=2, degR=0 — nothing
    /// lies beyond it, so the component is a (branched) contig, not two spurs off a
    /// backbone. A `deg(L) + deg(R) <= 1` terminus test sums to 2 here and misses
    /// it, trimming both arms and leaving a lone node.
    #[test]
    fn isolated_convergence_component_is_not_trimmed() {
        let mut p = Params::r10_4_1().graph;
        p.tip_steps = 10;
        let mut g = FuzzyGraph::new();
        for _ in 0..3 {
            g.add_node(5, 80);
        }
        g.add_edge(e(0, 2, 3, 0));
        g.add_edge(e(1, 2, 3, 0));
        let removed = trim_tips(&mut g, &p);
        assert_eq!(removed, 0, "nothing lies beyond the convergence at 2");
        assert!(
            g.is_live(0) && g.is_live(1) && g.is_live(2),
            "both arms were trimmed into node 2, erasing the component"
        );
    }

    /// A convergence that *does* have a continuation is a real junction, so short
    /// arms hanging off it are real tips. Guards the fix above against over-reach.
    ///
    /// `0R→2L`, `1R→2L` are the length-1 arms; `2R→3L→4L→5L→6L` is a backbone
    /// LONGER than `tip_steps`, so the walk inward from node 6 exits on the length
    /// bound and the backbone is not itself trimmed as a tip.
    #[test]
    fn arms_off_a_convergence_with_continuation_are_still_trimmed() {
        let mut p = Params::r10_4_1().graph;
        p.tip_steps = 2; // arms (len 1) trim; the 4-node backbone does not
        let mut g = FuzzyGraph::new();
        for _ in 0..7 {
            g.add_node(5, 80);
        }
        g.add_edge(e(0, 2, 3, 0)); // arm
        g.add_edge(e(1, 2, 3, 0)); // arm
        g.add_edge(e(2, 3, 3, 0)); // backbone beyond the convergence
        g.add_edge(e(3, 4, 3, 0));
        g.add_edge(e(4, 5, 3, 0));
        g.add_edge(e(5, 6, 3, 0));
        let removed = trim_tips(&mut g, &p);
        assert_eq!(
            removed, 2,
            "both short arms hang off a junction with a backbone"
        );
        assert!(!g.is_live(0) && !g.is_live(1), "the two spurs must trim");
        assert!(
            (2..=6).all(|n| g.is_live(n)),
            "the junction and its backbone survive"
        );
    }

    fn ed(from: u32, fd: Dir, to: u32, td: Dir, cov: u32) -> Edge {
        Edge {
            from,
            from_dir: fd,
            to,
            to_dir: td,
            cov,
            support: cov as f32,
            offset: 0,
        }
    }

    /// REGRESSION: arms leaving OPPOSITE ends of the source are not a bubble. A
    /// bubble is two arms off the SAME end that reconverge; arms off opposite ends
    /// are the chain continuing in both directions, and if they meet (a cycle
    /// through s) they must be left alone.
    ///
    /// `exits` was declared outside the `for send in [L, R]` loop, so the L-arm and
    /// the R-arm collided on their shared exit and one interior was deleted,
    /// breaking the cycle.
    ///
    /// Node 0 branches on BOTH ends (L -> 1,2 ; R -> 3,4). Interiors 1 and 3 both
    /// reach exit 5. Only a per-end `exits` keeps them apart.
    #[test]
    fn arms_off_opposite_ends_are_not_a_bubble() {
        let p = Params::r10_4_1().graph;
        let mut g = FuzzyGraph::new();
        for _ in 0..7 {
            g.add_node(5, 80);
        }
        g.add_edge(ed(1, Dir::Fwd, 0, Dir::Fwd, 9)); // 1R -> 0L
        g.add_edge(ed(2, Dir::Fwd, 0, Dir::Fwd, 3)); // 2R -> 0L   (dead-end arm)
        g.add_edge(ed(0, Dir::Fwd, 3, Dir::Fwd, 9)); // 0R -> 3L
        g.add_edge(ed(0, Dir::Fwd, 4, Dir::Fwd, 3)); // 0R -> 4L   (dead-end arm)
        g.add_edge(ed(5, Dir::Fwd, 1, Dir::Fwd, 9)); // 5R -> 1L
        g.add_edge(ed(3, Dir::Fwd, 5, Dir::Fwd, 9)); // 3R -> 5L
        g.add_edge(ed(5, Dir::Fwd, 6, Dir::Fwd, 9)); // 5R -> 6L  (makes 5 a junction)

        let removed = pop_bubbles(&mut g, &p);
        assert_eq!(
            removed, 0,
            "1 (off 0's L end) and 3 (off 0's R end) share exit 5, but they are not \
             a bubble — they are opposite directions of one chain"
        );
        assert!(
            g.is_live(1) && g.is_live(3),
            "neither interior may be popped"
        );
    }

    /// REGRESSION: the mirror image of the test above, at the *arrival* end. Both
    /// arms leave the SAME end of 0 (its L end), but reach node 5 by OPPOSITE ends
    /// — 1 lands on 5's L, 2 lands on 5's R. That is a cycle passing through 5, not
    /// a bubble: a bubble's arms must re-enter the exit by the same end. Keying
    /// `exits` on the bare node id calls it a bubble and pops a real chain interior.
    #[test]
    fn arms_arriving_at_opposite_ends_of_one_exit_are_not_a_bubble() {
        let p = Params::r10_4_1().graph;
        let mut g = FuzzyGraph::new();
        for _ in 0..7 {
            g.add_node(5, 80);
        }
        g.add_edge(ed(1, Dir::Fwd, 0, Dir::Fwd, 3)); // 1R -> 0L
        g.add_edge(ed(2, Dir::Fwd, 0, Dir::Fwd, 9)); // 2R -> 0L   (same end of 0)
        g.add_edge(ed(5, Dir::Rev, 1, Dir::Fwd, 3)); // 5L -- 1L   (arrives at 5's L)
        g.add_edge(ed(5, Dir::Fwd, 2, Dir::Fwd, 9)); // 5R -> 2L   (arrives at 5's R)
        g.add_edge(ed(6, Dir::Fwd, 5, Dir::Fwd, 9)); // 6R -> 5L

        let removed = pop_bubbles(&mut g, &p);
        assert_eq!(
            removed, 0,
            "arms off 0's L end arrive at OPPOSITE ends of 5 — a cycle, not a bubble"
        );
        assert!(
            g.is_live(1) && g.is_live(2),
            "node 1 is a legitimate chain interior and must not be popped"
        );
    }

    /// Positive control for the test above: when two arms off the same end of 0 do
    /// arrive at the SAME end of the exit, that IS a bubble and the weaker arm pops.
    /// Guards the `(node, entry end)` key against never matching.
    ///
    /// Arm strength is the sum of NODE coverage along the arm (not edge cov), so
    /// the two interiors are given different node covs.
    #[test]
    fn arms_arriving_at_the_same_end_of_one_exit_are_a_bubble() {
        let p = Params::r10_4_1().graph;
        let mut g = FuzzyGraph::new();
        g.add_node(5, 80); // 0: source junction
        g.add_node(9, 80); // 1: strong arm interior
        g.add_node(3, 80); // 2: weak arm interior
        g.add_node(5, 80); // 3: exit junction
        g.add_node(5, 80); // 4: continuation past the exit
        g.add_node(5, 80); // 5: continuation before the source
        g.add_edge(ed(0, Dir::Fwd, 1, Dir::Fwd, 9)); // 0R -> 1L
        g.add_edge(ed(0, Dir::Fwd, 2, Dir::Fwd, 3)); // 0R -> 2L   (same end of 0)
        g.add_edge(ed(1, Dir::Fwd, 3, Dir::Fwd, 9)); // 1R -> 3L   (both arrive at
        g.add_edge(ed(2, Dir::Fwd, 3, Dir::Fwd, 3)); // 2R -> 3L    3's L end)
        g.add_edge(ed(3, Dir::Fwd, 4, Dir::Fwd, 9)); // 3R -> 4L  (3 is a junction)
        g.add_edge(ed(5, Dir::Fwd, 0, Dir::Fwd, 9)); // 5R -> 0L  (0 is a junction)

        let removed = pop_bubbles(&mut g, &p);
        assert_eq!(removed, 1, "a real bubble: the weaker arm pops");
        assert!(!g.is_live(2), "arm 2 (cov 3) is weaker than arm 1 (cov 9)");
        assert!(g.is_live(1), "the stronger arm survives");
    }
}
