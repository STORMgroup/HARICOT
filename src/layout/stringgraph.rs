//! OLC-style read/string overlap graph — an alternative *layout* backend.
//!
//! Where the fuzzy-Bruijn graph ([`crate::graph`]) slices overlaps into shared
//! K-bin nodes, this backend keeps each **read** as a node and turns every
//! accepted [`Overlap`] into a directed edge, then reduces the graph the way an
//! overlap-layout-consensus assembler does:
//!
//!   1. **Classify** each overlap from its bin spans into *containment*
//!      (one read fully inside the other) vs *dovetail* (suffix-prefix) vs
//!      *internal* (partial — a likely repeat/chimera).
//!   2. **Take contained reads off the backbone** — their dovetails add no layout
//!      information and create tangles, so the dovetails are dropped and the read
//!      is re-attached as a leaf of its best container (see [`StringGraph`]).
//!   3. Keep dovetails as directed edges and apply **Myers transitive reduction**
//!      (remove `u→w` when `u→v→w` exists): the read-string-graph reduction from
//!      Myers 2005, "The fragment assembly string graph" (Bioinformatics 21).
//!      wtdbg2's bin-graph analog is `myers_transitive_reduction_graph`
//!      (`extern/wtdbg2/wtdbg-graph.h:1877`).
//!
//! Where each cluster of reads is a clique of overlaps, this is intended to collapse
//! the cluster to a single simple chain → one connected component with no branching.
//! Cross-cluster repeat-bridge overlaps are kept out by the dovetail acceptance floor
//! ([`StringGraphParams::min_dovetail_bins`]).

use crate::graph::fbg::Dir;
use crate::overlap::Overlap;

/// Tunables for the string-graph layout backend.
#[derive(Debug, Clone, Copy)]
pub struct StringGraphParams {
    /// A read end is treated as "flush" (fully aligned, no overhang) when its
    /// hang is within this many events. Larger ⇒ more overlaps read as clean
    /// dovetails/containments instead of internal.
    pub flush_tol_events: u32,
    /// Minimum aligned span (in bins) for a dovetail to be trusted as an edge.
    /// Raising this above the overlap gate's `min_aln_bins` drops the weak
    /// cross-cluster "repeat bridge" overlaps that would fuse two clusters.
    pub min_dovetail_bins: u32,
    /// Minimum aligned span (in bins) to trust a **containment**. Containment
    /// attaches a whole read into a component, so a false one is more damaging
    /// than a false dovetail — this is gated separately (and usually stricter).
    pub min_contain_bins: u32,
    /// Transitive-reduction fuzz, in events (Myers' `FUZZ`): slack on the
    /// two-hop length budget so a slightly-longer indirect path still reduces a
    /// direct edge.
    pub trans_fuzz_events: u32,
}

impl Default for StringGraphParams {
    fn default() -> Self {
        StringGraphParams {
            // ~2 DP bins of slack at the default bin_events=24.
            flush_tol_events: 48,
            // 16 bins ≈ 384 events: above the gate-floor (≈8-13 bin) cross-cluster
            // repeat overlaps, well below the within-cluster median (~50 bins), so
            // genuine edges survive while repeat bridges are dropped. Calibrated on
            // overlap_clusters_long alone (4 clusters → 4 clean components) — one
            // fixture, so treat it as a starting point rather than a tuned value.
            min_dovetail_bins: 16,
            min_contain_bins: 16,
            trans_fuzz_events: 48,
        }
    }
}

/// How an overlap relates the two reads (forward strand only, as the overlapper
/// emits).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OvlClass {
    /// `child` is fully contained in `container`.
    Contained { child: u32, container: u32 },
    /// A clean suffix-prefix dovetail `from → to`.
    Dovetail(StrEdge),
    /// Partial / off-end overlap — kept out of the layout (repeat or chimera).
    Internal,
}

/// A directed dovetail edge in the read string graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StrEdge {
    pub from: u32,
    pub to: u32,
    pub from_dir: Dir,
    pub to_dir: Dir,
    /// Overlap length, in events (for the GFA `*M` field).
    pub overlap_events: u32,
    /// Myers edge length = events the `to` read adds beyond the overlap (its far
    /// overhang). Drives transitive reduction's length budget.
    pub ext_events: u32,
}

/// The reduced read string graph.
///
/// Every read is a node. Non-contained reads form the **dovetail backbone**
/// ([`edges`](StringGraph::edges), after transitive reduction). Each contained
/// read is **attached** as a leaf to its best-aligned container
/// ([`contain_edges`](StringGraph::contain_edges)) so it stays in the right
/// component (and the cluster keeps its node count) instead of being discarded.
#[derive(Debug, Clone, Default)]
pub struct StringGraph {
    /// `contained[r]` ⇒ read `r` is contained (attached as a leaf, not on the
    /// dovetail backbone).
    pub contained: Vec<bool>,
    /// Reduced dovetail backbone edges (between non-contained reads).
    pub edges: Vec<StrEdge>,
    /// `(container, child)` attachment edges for contained reads.
    pub contain_edges: Vec<(u32, u32)>,
}

/// Midpoint (events) of read `r`'s aligned bin span `(lo, hi)` (hi inclusive).
fn span_mid_events(bin_lo: u32, bin_hi: u32, bin_events: u32) -> f64 {
    (bin_lo as f64 + (bin_hi as f64 + 1.0)) * 0.5 * bin_events as f64
}

/// Classify a single overlap given the two reads' event lengths.
///
/// The bin-DP backtraces between anchor cells and does not extend to read ends, so
/// its spans are a short *interior* aligned core and testing them directly for
/// end-flushness misclassifies real overlaps. Instead we do what an OLC overlapper
/// does: take the alignment **diagonal** `a = slope·b + c` fit from the core, then
/// **extrapolate** it to project read `b`'s whole `[0, len_b)` interval into read
/// `a`'s coordinate frame and compare the two intervals. `slope` comes from the
/// overlap ([`Overlap::slope`], 1.0 for the wtdbg2 baseline).
///
/// `len` is indexed by read id.
pub fn classify(ovl: &Overlap, len: &[u32], p: &StringGraphParams) -> OvlClass {
    if ovl.aligned_bins < p.min_dovetail_bins {
        return OvlClass::Internal;
    }
    let tol = p.flush_tol_events as f64;
    let bin = ovl.bin_events;
    let (la, lb) = (len[ovl.a as usize] as f64, len[ovl.b as usize] as f64);
    let slope = ovl.slope;

    // Diagonal through the core's midpoint: a = slope·b + c.
    let a_mid = span_mid_events(ovl.a_bin_span.0, ovl.a_bin_span.1, bin);
    let b_mid = span_mid_events(ovl.b_bin_span.0, ovl.b_bin_span.1, bin);
    let c = a_mid - slope * b_mid;

    // Read b's full extent projected into read a's coordinate frame.
    let b0 = c; // a-coord of b's start (b_pos = 0)
    let b1 = slope * lb + c; // a-coord of b's end (b_pos = lb)
    // a occupies [0, la].

    // Reject pairs whose projected intervals barely meet (extrapolation artifact
    // or off-diagonal noise): require a real overlap of the projected intervals.
    let ov_lo = b0.max(0.0);
    let ov_hi = b1.min(la);
    let overlap_events = ov_hi - ov_lo;
    if overlap_events <= tol {
        return OvlClass::Internal;
    }

    // Containment: one read's projected interval sits inside the other (±tol).
    let b_in_a = b0 >= -tol && b1 <= la + tol;
    let a_in_b = b0 <= tol && b1 >= la - tol;
    if b_in_a || a_in_b {
        // Containment attaches a whole read, so a false one is costly: gate it
        // more strictly than a dovetail. A too-weak core ⇒ treat as Internal
        // (dropped) rather than risk a cross-cluster bridge.
        if ovl.aligned_bins < p.min_contain_bins {
            return OvlClass::Internal;
        }
        return if b_in_a && a_in_b {
            if la <= lb {
                OvlClass::Contained {
                    child: ovl.a,
                    container: ovl.b,
                }
            } else {
                OvlClass::Contained {
                    child: ovl.b,
                    container: ovl.a,
                }
            }
        } else if b_in_a {
            OvlClass::Contained {
                child: ovl.b,
                container: ovl.a,
            }
        } else {
            OvlClass::Contained {
                child: ovl.a,
                container: ovl.b,
            }
        };
    }

    let ov = overlap_events as u32;
    // Dovetail a→b: b starts after a (b0 > 0) and extends past a's end (b1 > la).
    // a's suffix == b's prefix. b adds `b1 - la` events to the right.
    if b0 > tol && b1 > la + tol {
        return OvlClass::Dovetail(StrEdge {
            from: ovl.a,
            to: ovl.b,
            from_dir: Dir::Fwd,
            to_dir: Dir::Fwd,
            overlap_events: ov,
            ext_events: (b1 - la).max(0.0) as u32,
        });
    }
    // Dovetail b→a: b starts before a (b0 < 0) and ends inside a (b1 < la).
    // b's suffix == a's prefix. a adds `la - b1` events to the right.
    if b0 < -tol && b1 < la - tol {
        return OvlClass::Dovetail(StrEdge {
            from: ovl.b,
            to: ovl.a,
            from_dir: Dir::Fwd,
            to_dir: Dir::Fwd,
            overlap_events: ov,
            ext_events: (la - b1).max(0.0) as u32,
        });
    }

    OvlClass::Internal
}

/// Build the reduced read string graph from accepted overlaps and read lengths.
pub fn build_string_graph(
    n_reads: usize,
    overlaps: &[Overlap],
    len: &[u32],
    p: &StringGraphParams,
) -> StringGraph {
    let mut contained = vec![false; n_reads];
    // Best container per contained read, chosen by the *strongest* overlap
    // (`aligned_bins`). This tends to route a contained read to a same-cluster
    // container: on the cluster fixtures the spurious cross-cluster containments sit
    // near the gate floor, so they lose to the read's genuine within-cluster one.
    let mut container_of: Vec<Option<(u32, u32)>> = vec![None; n_reads]; // (container, aligned_bins)
    let mut edges: Vec<StrEdge> = Vec::new();

    // Pass 1: classify; mark contained reads and collect dovetails.
    for ovl in overlaps {
        match classify(ovl, len, p) {
            OvlClass::Contained { child, container } => {
                contained[child as usize] = true;
                let better = match container_of[child as usize] {
                    None => true,
                    Some((_, best_aln)) => ovl.aligned_bins > best_aln,
                };
                if better {
                    container_of[child as usize] = Some((container, ovl.aligned_bins));
                }
            }
            OvlClass::Dovetail(e) => edges.push(e),
            OvlClass::Internal => {}
        }
    }

    // Drop any dovetail touching a contained read — it carries no backbone value.
    edges.retain(|e| !contained[e.from as usize] && !contained[e.to as usize]);
    let edges = transitive_reduction(edges, n_reads, p.trans_fuzz_events);

    // Attach each contained read to its best-aligned container (a leaf edge).
    let mut contain_edges = Vec::new();
    for child in 0..n_reads {
        if contained[child] {
            if let Some((cont, _)) = container_of[child] {
                contain_edges.push((cont, child as u32));
            }
        }
    }
    StringGraph {
        contained,
        edges,
        contain_edges,
    }
}

/// Myers (2005) transitive reduction over the directed dovetail edges: remove
/// `u→w` whenever a two-hop `u→v→w` exists within the length budget. Faithful to
/// the marking scheme in §"Transitive Edge Reduction"; wtdbg2's bin analog is
/// `myers_transitive_reduction_graph` (`wtdbg-graph.h:1877`).
fn transitive_reduction(edges: Vec<StrEdge>, n_reads: usize, fuzz: u32) -> Vec<StrEdge> {
    // Adjacency: per node, its out-edges as (to, ext, edge_index), sorted by ext.
    let mut adj: Vec<Vec<(u32, u32, usize)>> = vec![Vec::new(); n_reads];
    for (i, e) in edges.iter().enumerate() {
        adj[e.from as usize].push((e.to, e.ext_events, i));
    }
    for v in adj.iter_mut() {
        v.sort_by_key(|&(_, ext, _)| ext);
    }

    // Mark state per node: 0 vacant, 1 inplay, 2 eliminated.
    const VACANT: u8 = 0;
    const INPLAY: u8 = 1;
    const ELIMINATED: u8 = 2;
    let mut mark = vec![VACANT; n_reads];
    let mut reduced = vec![false; edges.len()];

    for v in 0..n_reads {
        let out = &adj[v];
        if out.is_empty() {
            continue;
        }
        for &(w, _, _) in out {
            mark[w as usize] = INPLAY;
        }
        // Longest reachable single hop + fuzz bounds the two-hop budget.
        let longest = out.last().map(|&(_, ext, _)| ext).unwrap_or(0) + fuzz;
        for &(w, ext_vw, _) in out {
            if mark[w as usize] != INPLAY {
                continue;
            }
            for &(x, ext_wx, _) in &adj[w as usize] {
                if ext_vw + ext_wx > longest {
                    break; // sorted by ext ⇒ all further are longer too
                }
                if mark[x as usize] == INPLAY {
                    mark[x as usize] = ELIMINATED;
                }
            }
        }
        for &(w, _, idx) in out {
            if mark[w as usize] == ELIMINATED {
                reduced[idx] = true;
            }
            mark[w as usize] = VACANT;
        }
    }

    edges
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !reduced[*i])
        .map(|(_, e)| e)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Strand;

    /// Build a minimal [`Overlap`] for classification tests (only the fields
    /// `classify` reads are populated).
    fn ovl(a: u32, b: u32, abin: (u32, u32), bbin: (u32, u32), aligned_bins: u32) -> Overlap {
        Overlap {
            a,
            b,
            strand: Strand::Forward,
            n_anchors: 0,
            matched: 0,
            cnt: 0,
            aligned_bins,
            bin_events: 10,
            gap: 0,
            slope: 1.0,
            a_bin_span: abin,
            b_bin_span: bbin,
            path: vec![],
            // Test fixture: z is not under test. NaN = not measured.
            z: f32::NAN,
        }
    }

    fn tol5() -> StringGraphParams {
        StringGraphParams {
            flush_tol_events: 5,
            min_dovetail_bins: 8,
            min_contain_bins: 8,
            trans_fuzz_events: 48,
        }
    }

    #[test]
    fn classifies_dovetail_a_to_b() {
        // a core mid 75, b core mid 25 → diagonal c=50: b projects to a-coords
        // [50, 150). a is [0,100). a's suffix overlaps b's prefix ⇒ a→b.
        let len = vec![100u32, 100];
        let o = ovl(0, 1, (5, 9), (0, 4), 10);
        match classify(&o, &len, &tol5()) {
            OvlClass::Dovetail(e) => {
                assert_eq!((e.from, e.to), (0, 1));
                assert_eq!(e.from_dir, Dir::Fwd);
            }
            other => panic!("expected dovetail a→b, got {other:?}"),
        }
    }

    #[test]
    fn classifies_dovetail_b_to_a() {
        // Mirror: b projects to [-50, 50) ⇒ b's suffix overlaps a's prefix ⇒ b→a.
        let len = vec![100u32, 100];
        let o = ovl(0, 1, (0, 4), (5, 9), 10);
        match classify(&o, &len, &tol5()) {
            OvlClass::Dovetail(e) => assert_eq!((e.from, e.to), (1, 0)),
            other => panic!("expected dovetail b→a, got {other:?}"),
        }
    }

    #[test]
    fn classifies_containment() {
        // b (len 50) projects to a-coords [30, 80) ⊂ a's [0,200) ⇒ b contained.
        let len = vec![200u32, 50];
        let o = ovl(0, 1, (3, 7), (0, 4), 10);
        assert_eq!(
            classify(&o, &len, &tol5()),
            OvlClass::Contained {
                child: 1,
                container: 0
            }
        );
    }

    #[test]
    fn weak_overlap_is_internal() {
        let len = vec![100u32, 100];
        let p = StringGraphParams {
            min_dovetail_bins: 8,
            ..Default::default()
        };
        let o = ovl(0, 1, (5, 9), (0, 4), 4); // only 4 aligned bins
        assert_eq!(classify(&o, &len, &p), OvlClass::Internal);
    }

    /// A clique-like tiling reduces to a simple chain: reads 0..4 each dovetail to
    /// the next AND skip to further reads; transitive reduction must leave only
    /// the 4 adjacent edges.
    #[test]
    fn transitive_reduction_collapses_clique_to_chain() {
        // Reads of length 100, each starting 20 ev later → 0,20,40,60,80.
        let n = 5;
        let len = vec![100u32; n];
        let bin = 10u32;
        let mut overlaps = Vec::new();
        // For every i<j that overlap (start_j < end_i = start_i+100), make a→b.
        let start = [0u32, 20, 40, 60, 80];
        for i in 0..n {
            for j in (i + 1)..n {
                let si = start[i];
                let sj = start[j];
                if sj >= si + 100 {
                    continue; // no overlap
                }
                // overlap region in absolute ev: [sj, si+100)
                let ov_lo = sj;
                let ov_hi = si + 100;
                // read i (=a) local coords: ov - si ; read j (=b) local: ov - sj
                let a_lo = (ov_lo - si) / bin;
                let a_hi = (ov_hi - si - 1) / bin;
                let b_lo = (ov_lo - sj) / bin;
                let b_hi = (ov_hi - sj - 1) / bin;
                overlaps.push(ovl(
                    i as u32,
                    j as u32,
                    (a_lo, a_hi),
                    (b_lo, b_hi),
                    (a_hi - a_lo + 1).min(b_hi - b_lo + 1),
                ));
            }
        }
        let p = StringGraphParams {
            min_dovetail_bins: 1,
            flush_tol_events: 5,
            trans_fuzz_events: 5,
            ..Default::default()
        };
        let g = build_string_graph(n, &overlaps, &len, &p);
        assert!(g.contained.iter().all(|&c| !c), "no read is contained here");
        // Exactly the 4 adjacent edges 0→1→2→3→4 survive.
        assert_eq!(g.edges.len(), 4, "edges: {:?}", g.edges);
        let mut pairs: Vec<(u32, u32)> = g.edges.iter().map(|e| (e.from, e.to)).collect();
        pairs.sort();
        assert_eq!(pairs, vec![(0, 1), (1, 2), (2, 3), (3, 4)]);
    }
}
