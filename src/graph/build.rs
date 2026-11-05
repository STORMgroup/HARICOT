//! Fuzzy-Bruijn-graph construction from bin-DP overlaps.
//!
//! A port of wtdbg2's node/edge build (`extern/wtdbg2/wtdbg.h`), working entirely in **bins**. The
//! domain change vs wtdbg2: a bin is `bin_events` events (not 256 bp), and a node
//! is **K = `node_bins` consecutive DP bins**, so its length is the *derived*
//! `node_bins * bin_events` events.
//!
//! K is the parameter and is never re-derived from a length: `K = round(node_events
//! / bin_events)` is exactly the inference [`Grid`] exists to eliminate, because it
//! silently reinterprets a node whenever `bin_events` changes (wtdbg2 converts once,
//! `wtdbg.c:638-640`, after which `reglen` *is* K). `--node-events` survives only as
//! a deprecated alias that rounds up to a whole bin and warns.
//!
//! Stages:
//!   1. **`hit2rdregs`** (`wtdbg.h:459`): slice each overlap's bin alignment into
//!      K-bin windows on the **query** read's grid; each window projects the other
//!      read's aligned interval into a *region* keyed by `(query, window)`. A
//!      region is *closed* if its projected length leaves `[K−1, K+1]` — the ±1
//!      "fuzzy" tolerance (`wtdbg.h:570`).
//!   2. **`build_nodes`** (`wtdbg.h:1642`): regions sharing a K-bin key are a
//!      candidate node.
//!   3. **Ω conflict resolution** (`mul_update_regs_graph`, `wtdbg.h:1048`):
//!      process candidates by coverage descending (see the Stage 3 comment for why
//!      that is what wtdbg2 executes); per read drop duplicate regions and close any that overlap
//!      an already-committed K-bin by more than `regovl`; reject the node if too few
//!      survive (`min_node_cov`) or too many conflicted (`node_drop`); else commit
//!      and thread each region onto its read.
//!   4. **`build_edges`** (`wtdbg.h:2019`): connect consecutive regions on a read
//!      (read threading); the median per-read offset is the edge offset.
//!
//! As in wtdbg2, regions are anchored on **both** the query and target grids
//! (`hit2rdregs` blocks 1 & 2, `wtdbg.h:534`/`:596`), with globally-keyed node ids
//! `(read, window)`, so a node accumulates every read that aligns to that locus
//! regardless of overlap direction. Ω conflict resolution (`mul_update_regs_graph`,
//! `wtdbg.h:1048`) commits the highest-coverage candidate first: one read's tiling of
//! a locus claims the shared bins and the opposite-grid duplicate then conflicts and
//! is dropped, rather than both surviving as a phase-shifted double-tiling (see the
//! `both_grids_collapse_to_single_tiling` test).
//!
//! Every overlap upstream is built with `Strand::Forward` (`overlap::mod`,
//! `overlap::kbm`), so all edge directions here are forward.

use std::collections::HashMap;

use crate::graph::fbg::{Dir, Edge, FuzzyGraph};
use crate::overlap::Overlap;
use crate::params::GraphParams;

/// A region: an interval of a read's bins assigned to one candidate K-bin.
/// Analog of wtdbg2 `rd_reg_t` (`wtdbg.h:46`).
#[derive(Debug, Clone, Copy)]
struct RdReg {
    /// Global K-bin key `((query_read << 32) | window)` — unique per anchor+window.
    key: u64,
    /// Read this interval belongs to.
    rid: u32,
    /// Start bin in this read's own coordinates.
    beg: u32,
    /// End bin (exclusive).
    end: u32,
    /// Closed by the ±1 fuzz (not a valid K-bin projection).
    closed: bool,
}

/// A committed region threaded onto a read (analog of `reg_t`, `wtdbg.h:58`).
#[derive(Debug, Clone, Copy)]
struct Threaded {
    node: u32,
    beg: u32,
    end: u32,
}

/// Provenance of a node, emitted alongside the graph for diagnostics/coloring:
/// the anchor (query) read it was keyed on and its event range on that read.
#[derive(Debug, Clone, Copy)]
pub struct NodeAnno {
    /// Index of the anchor read this K-bin was keyed on.
    pub anchor_read: u32,
    /// K-bin window index on the anchor read.
    pub window: u32,
    /// Event range `[ev_beg, ev_end)` of the K-bin on the anchor read.
    pub ev_beg: u32,
    pub ev_end: u32,
}

#[inline]
fn mk_key(anchor: u32, window: u32) -> u64 {
    ((anchor as u64) << 32) | window as u64
}

/// The bin/K-bin grid a build runs on. Derived once per build from the params and
/// the bin size the overlaps were aligned at, so that K, the node length and
/// `regovl` cannot disagree with each other.
///
/// K comes straight from [`GraphParams::node_bins`]. Do *not* re-derive it by
/// dividing an event count by `bin_events`: that lets a node's declared length
/// differ from the ground it actually covers.
#[derive(Debug, Clone, Copy)]
struct Grid {
    /// Bin size in events (from the overlaps).
    bin_events: u32,
    /// K: bins per node.
    k: u32,
    /// Region-overlap conflict threshold, in bins.
    regovl: u32,
}

impl Grid {
    fn new(overlaps: &[Overlap], params: &GraphParams) -> Grid {
        Grid {
            bin_events: overlaps[0].bin_events.max(1),
            k: params.node_bins.max(1),
            regovl: params.regovl_bins,
        }
    }
    /// Node length in events. The single definition of a node's size.
    #[inline]
    fn node_events(&self) -> u32 {
        self.k * self.bin_events
    }
}

/// Emit K-bin candidate regions for **one** alignment direction. wtdbg2 runs this
/// for both directions of every hit (`hit2rdregs` blocks 1 & 2): `owner` is the
/// read whose window grid defines the node ids, `partner` is the other read.
/// `proj[owner_bin] = (partner_lo, partner_hi)` is the alignment projected onto the
/// owner's bins; `span` is the owner's aligned bin span. Each owner window becomes a
/// node carrying the partner's projected interval plus the owner's own interval, so
/// the same locus keyed on either read collects the same set of reads.
fn emit_grid_regions(
    regions: &mut Vec<RdReg>,
    self_emitted: &mut std::collections::HashSet<u64>,
    owner: u32,
    partner: u32,
    proj: &HashMap<u32, (u32, u32)>,
    span: (u32, u32),
    k: u32,
) {
    for w in (span.0 / k)..=(span.1 / k) {
        let (lo, hi) = (w * k, w * k + k);
        let (mut p_beg, mut p_end, mut matched) = (u32::MAX, 0u32, 0u32);
        for ob in lo..hi {
            if let Some(&(pl, ph)) = proj.get(&ob) {
                p_beg = p_beg.min(pl);
                p_end = p_end.max(ph);
                matched += 1;
            }
        }
        if matched == 0 {
            continue; // min_node_mats: need >=1 matched bin in the window
        }
        let key = mk_key(owner, w);
        // ±1 fuzz: close if the projected length leaves [K-1, K+1] (wtdbg.h:570).
        let len = p_end - p_beg + 1;
        let closed = len + 1 < k || len > k + 1;
        regions.push(RdReg {
            key,
            rid: partner,
            beg: p_beg,
            end: p_end + 1,
            closed,
        });
        // The owner read's own region for this window (emitted once per key).
        if self_emitted.insert(key) {
            regions.push(RdReg {
                key,
                rid: owner,
                beg: lo,
                end: hi,
                closed: false,
            });
        }
    }
}

/// Build the fuzzy-Bruijn graph from bin-DP overlaps. Node extent comes from the
/// overlaps alone.
pub fn build_graph(overlaps: &[Overlap], params: &GraphParams) -> FuzzyGraph {
    build_graph_annotated(overlaps, params).0
}

/// Like [`build_graph`] but also returns per-node provenance ([`NodeAnno`]),
/// indexed by node id, for diagnostics and reference-position coloring.
pub fn build_graph_annotated(
    overlaps: &[Overlap],
    params: &GraphParams,
) -> (FuzzyGraph, Vec<NodeAnno>) {
    build_graph_annotated_ranked(overlaps, params, None)
}

/// [`build_graph_annotated`] with an optional **canonical tie-break**. `read_rank[r]` is read `r`'s
/// tie rank, an order that depends neither on the reads' numbering nor on their names
/// (`haricot-detect --canonical-ties` passes `readorder::rank_by_content`, the order of the reads'
/// signal content hashes). Stage 3 commits candidates by coverage descending, and among
/// equal coverage by the K-bin key `(anchor read number, window)`: with `read_rank`, by
/// `(read_rank[anchor], window)` instead, so the commit order, the node ids (assigned in commit order)
/// and every pass after the build that iterates nodes or edges by id do not depend on how the reads
/// are numbered or named. `None`: ties by the K-bin key.
pub fn build_graph_annotated_ranked(
    overlaps: &[Overlap],
    params: &GraphParams,
    read_rank: Option<&[u32]>,
) -> (FuzzyGraph, Vec<NodeAnno>) {
    build_graph_dumping(overlaps, params, read_rank, None)
}

/// [`build_graph_annotated_ranked`] that also writes the Stage-3 read placements to `threads_dump`
/// (`None`: none), as `haricot-layout` writes OUT/threads.tsv. The file is created after Stage 2; an
/// empty overlap set writes none.
pub fn build_graph_dumping(
    overlaps: &[Overlap],
    params: &GraphParams,
    read_rank: Option<&[u32]>,
    threads_dump: Option<&std::path::Path>,
) -> (FuzzyGraph, Vec<NodeAnno>) {
    let mut g = FuzzyGraph::new();
    let mut annos: Vec<NodeAnno> = Vec::new();
    if overlaps.is_empty() {
        return (g, annos);
    }

    // All overlaps were aligned at the same bin size.
    let grid = Grid::new(overlaps, params);
    let (bin_events, k, regovl) = (grid.bin_events, grid.k, grid.regovl);

    // ---- Stage 1: hit2rdregs — slice each overlap into K-bin regions on BOTH the
    // query (a) and target (b) grids, so the same locus keyed on either read
    // collects the same reads and the redundant tiling can be resolved away. ----
    let mut regions: Vec<RdReg> = Vec::new();
    let mut self_emitted: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for o in overlaps {
        // Project the alignment both ways: a_bin -> b-range and b_bin -> a-range.
        let mut proj_ab: HashMap<u32, (u32, u32)> = HashMap::new();
        let mut proj_ba: HashMap<u32, (u32, u32)> = HashMap::new();
        for &(a_bin, b_bin) in &o.path {
            let e = proj_ab.entry(a_bin).or_insert((b_bin, b_bin));
            e.0 = e.0.min(b_bin);
            e.1 = e.1.max(b_bin);
            let e = proj_ba.entry(b_bin).or_insert((a_bin, a_bin));
            e.0 = e.0.min(a_bin);
            e.1 = e.1.max(a_bin);
        }
        if proj_ab.is_empty() {
            continue;
        }
        emit_grid_regions(
            &mut regions,
            &mut self_emitted,
            o.a,
            o.b,
            &proj_ab,
            o.a_bin_span,
            k,
        );
        emit_grid_regions(
            &mut regions,
            &mut self_emitted,
            o.b,
            o.a,
            &proj_ba,
            o.b_bin_span,
            k,
        );
    }

    // ---- Stage 2: build_nodes — group regions by K-bin key. ----
    regions.retain(|r| !r.closed);
    regions.sort_by_key(|r| (r.key, r.rid, r.beg));
    let mut candidates: Vec<(u64, Vec<RdReg>)> = Vec::new();
    for r in regions {
        match candidates.last_mut() {
            Some((k, v)) if *k == r.key => v.push(r),
            _ => candidates.push((r.key, vec![r])),
        }
    }

    // De-duplicate per read within each candidate (keep the longest region), and
    // record the raw depth = distinct supporting reads. Keep the K-bin `key` so the
    // committed node can be annotated with its anchor read + window.
    let mut cand: Vec<(u64, u32, Vec<RdReg>)> = Vec::with_capacity(candidates.len());
    for (key, regs) in candidates {
        let mut by_rid: HashMap<u32, RdReg> = HashMap::new();
        for r in regs {
            by_rid
                .entry(r.rid)
                .and_modify(|cur| {
                    if (r.end - r.beg) > (cur.end - cur.beg) {
                        *cur = r;
                    }
                })
                .or_insert(r);
        }
        let regs: Vec<RdReg> = by_rid.into_values().collect();
        let depth = regs.len() as u32;
        cand.push((key, depth, regs));
    }

    // ---- Stage 3: Ω conflict resolution. Process candidates by **coverage
    // descending** (wtdbg2 `mul_update_regs_graph`, `wtdbg.h:1048`): the
    // best-supported tiling of a locus commits first and claims its reads' bins, so
    // the opposite-grid duplicate then conflicts and is dropped. Key breaks ties
    // deterministically. ----
    //
    //
    // The order is coverage descending. wtdbg2 looks like it sorts by distance from the
    // median candidate depth (`wtdbg.h:1816`), but its `gcov` (`wtdbg.h:1809-1813`) is
    // not a median: the accumulation loop runs while `rank < nds->size` over a histogram
    // whose counts sum to exactly `nds->size`, so `gcov == max_depth + 1`, and sorting
    // ascending on `num_diff(cnt, gcov) = gcov - cnt` is sorting descending on `cnt` —
    // the commented-out `:1815`. Tie-break by key for determinism (wtdbg2 ties on
    // `rank` ascending).
    //
    // ⚠ Ordering on overlap-strength-weighted support instead reorders which tiling of a
    // locus commits first: measured, it cost ~50% of refNG50 on two bacteria (E. coli
    // 886,287 → 423,884, S. aureus D7 846,764 → 468,909) for +0.6% on a D. melanogaster
    // arm. The most NUMEROUS tiling is the one to commit.
    match read_rank {
        None => cand.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0))),
        Some(rank) => {
            // The key with the anchor's number replaced by its rank: a total order (window < 2^32).
            let ranked = |key: u64| (rank[(key >> 32) as usize], key & 0xffff_ffff);
            cand.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| ranked(a.0).cmp(&ranked(b.0))));
        }
    }

    // rid -> committed intervals tagged with the node that claimed them.
    let mut claims: HashMap<u32, Vec<(u32, u32, u32)>> = HashMap::new();
    let mut threads: HashMap<u32, Vec<Threaded>> = HashMap::new(); // rid -> committed regions

    // The Stage-3 dump (`threads_dump`): one `C` row per candidate K-bin (owner read, window, raw
    // depth, free regions, verdict) and one `T` row per committed region (read, node, bins). The
    // pipeline reads it back: the end-bridge and contig-end rescue steps (`haricot-layout`,
    // `haricot-rescue rescue`) take each read's node placements from it.
    let mut census = threads_dump.map(|p| {
        std::io::BufWriter::new(std::fs::File::create(p).expect("creating the threads dump"))
    });
    for (key, raw_depth, regs) in cand {
        let (owner, window) = ((key >> 32) as u32, (key & 0xffff_ffff) as u32);
        if raw_depth > params.max_node_cov {
            if let Some(w) = census.as_mut() {
                use std::io::Write;
                writeln!(w, "C\t{owner}\t{window}\t{raw_depth}\t-\tmaxcov\t-").ok();
            }
            continue; // high-coverage repeat (wtdbg.h:1883)
        }
        // Classify each region against already-committed K-bins: the node it conflicts
        // with (overlap > regovl ⇒ its bins are taken; `None` = free).
        let mut region_node: Vec<Option<u32>> = Vec::with_capacity(regs.len());
        for r in &regs {
            let conflict = claims.get(&r.rid).and_then(|owned| {
                owned
                    .iter()
                    .find(|&&(b, e, _)| e.min(r.end).saturating_sub(b.max(r.beg)) > regovl)
                    .map(|&(_, _, nd)| nd)
            });
            region_node.push(conflict);
        }
        let n_free = region_node.iter().filter(|x| x.is_none()).count() as u32;
        // Commit a new node from the free (non-conflicting) regions (wtdbg2: survivors =
        // free), unless too few survive or too many conflicted.
        let pass = n_free;
        if pass < params.min_node_cov
            || (pass as f64) < raw_depth as f64 * (1.0 - params.node_drop as f64)
        {
            if let Some(w) = census.as_mut() {
                use std::io::Write;
                let why = if pass < params.min_node_cov {
                    "mincov"
                } else {
                    "nodedrop"
                };
                writeln!(w, "C\t{owner}\t{window}\t{raw_depth}\t{pass}\t{why}\t-").ok();
            }
            continue;
        }
        let node = g.add_node(pass, grid.node_events());
        if let Some(w) = census.as_mut() {
            use std::io::Write;
            writeln!(w, "C\t{owner}\t{window}\t{raw_depth}\t{pass}\tnode\t{node}").ok();
        }
        let ev_beg = window * grid.node_events();
        annos.push(NodeAnno {
            anchor_read: owner,
            window,
            ev_beg,
            ev_end: ev_beg + grid.node_events(),
        });
        // Thread the free regions through the new node, claiming their bins. The regions
        // are one per read (deduplicated above), so the node's coverage is `pass`.
        for (r, rn) in regs.iter().zip(&region_node) {
            if rn.is_none() {
                claims.entry(r.rid).or_default().push((r.beg, r.end, node));
                threads.entry(r.rid).or_default().push(Threaded {
                    node,
                    beg: r.beg,
                    end: r.end,
                });
            }
        }
    }

    if let Some(w) = census.as_mut() {
        use std::io::Write;
        let mut rids: Vec<&u32> = threads.keys().collect();
        rids.sort_unstable();
        for rid in rids {
            for t in &threads[rid] {
                writeln!(w, "T\t{rid}\t{}\t{}\t{}", t.node, t.beg, t.end).ok();
            }
        }
    }
    drop(census);

    // ---- Stage 4: build_edges — read threading. Connect consecutive committed
    // regions on each read; accumulate coverage and median offset per canonical
    // edge. ----
    //
    // Deviation from a literal `max_edge_span` gate (wtdbg2's `g->max_edge_span`
    // check, `wtdbg.h:2053`): by default we connect consecutive committed regions on
    // a read regardless of the gap between them. Rationale — a read physically spans
    // its own DNA, so two committed loci adjacent in the read's own ordering are
    // contiguous; the intervening "gap" only means no node was committed there (low
    // local coverage or conflict-rejected candidates), not a discontinuity. Capping
    // the span at `max_edge_span_events=1024` (≈43 bins) instead dropped these
    // bridges; each read cluster was then observed to break into 8-14 components.
    // `max_edge_span_events == 0` selects the bridge-everything behaviour; any
    // positive value caps the span.
    // i64 so the "unlimited" sentinel can exceed any realistic gap without wrap.
    let max_span_bins: i64 = if params.max_edge_span_events == 0 {
        i64::MAX
    } else {
        (params.max_edge_span_events as f64 / bin_events as f64).round() as i64
    };
    // (from, from_dir, to, to_dir) -> offsets in events
    let mut edge_offsets: HashMap<(u32, u8, u32, u8), Vec<i32>> = HashMap::new();
    for regs in threads.values_mut() {
        regs.sort_by_key(|t| t.beg);
        for pair in regs.windows(2) {
            let (r1, r2) = (pair[0], pair[1]);
            if r1.node == r2.node {
                continue;
            }
            let gap_bins = r2.beg as i64 - r1.end as i64;
            if gap_bins > max_span_bins {
                continue; // beyond max_edge_span — leave a break
            }
            let offset = (gap_bins * bin_events as i64) as i32;
            // Orientation is *read order*: r1 precedes r2 along the read, so the
            // edge runs r1 -> r2. Sketching is forward-strand only
            // (`pipeline::sketch_reads`), so every read agrees on this order.
            // Do NOT canonicalize on node id: `from`/`to` is what the GFA writer
            // turns into the L-line orientation, and node ids come from
            // candidate-region iteration order, which is unrelated to genome
            // position. Keying on id order also negates `offset`, colliding with
            // the writer's and `clean::step_len`'s "negative = overlap".
            let key = (r1.node, Dir::Fwd as u8, r2.node, Dir::Fwd as u8);
            edge_offsets.entry(key).or_default().push(offset);
        }
    }
    // Read placements (for threading and the fragment layout) and read-threaded triples: the path
    // context the pairwise reduction below discards.
    if params.read_threading || params.frg_layout {
        let mut rt = crate::graph::threading::ReadThreads::default();
        for (&rid, regs) in threads.iter() {
            rt.insert(
                rid,
                regs.iter()
                    .map(|t| crate::graph::threading::ReadPlacement {
                        node: t.node,
                        beg: t.beg,
                        end: t.end,
                    })
                    .collect(),
            );
        }
        g.set_read_threads(rt);
    }
    if params.read_threading {
        let mut tr = crate::graph::threading::Transitions::default();
        for regs in threads.values() {
            for w in regs.windows(3) {
                if w[0].node != w[1].node && w[1].node != w[2].node {
                    *tr.entry((w[0].node, w[1].node, w[2].node)).or_insert(0) += 1;
                }
            }
        }
        g.set_transitions(tr);
    }

    // Insert in key order: the map's iteration order is randomly seeded, which made the edge
    // ids — and so the GFA's link order — differ from run to run (the content never did).
    let mut edge_offsets: Vec<_> = edge_offsets.into_iter().collect();
    edge_offsets.sort_unstable_by_key(|e| e.0);
    for ((from, _fd, to, _td), mut offs) in edge_offsets {
        let cov = offs.len() as u32;
        // `clean` filters again, so this only bites `--raw`, which skips cleaning —
        // without it `--min-edge-cov` is a no-op on the raw graph.
        if cov < params.min_edge_cov {
            continue;
        }
        offs.sort_unstable();
        let offset = offs[offs.len() / 2]; // median
        g.add_edge(Edge {
            from,
            from_dir: Dir::Fwd,
            to,
            to_dir: Dir::Fwd,
            cov,
            support: cov as f32,
            offset,
        });
    }

    (g, annos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Params, Strand};

    /// Make an overlap with a slope-1 path covering `nbins` aligned bins.
    fn overlap(a: u32, b: u32, nbins: u32, bin_events: u32) -> Overlap {
        let path: Vec<(u32, u32)> = (0..nbins).map(|i| (i, i)).collect();
        Overlap {
            a,
            b,
            strand: Strand::Forward,
            n_anchors: nbins,
            matched: nbins * bin_events,
            cnt: nbins,
            aligned_bins: nbins,
            bin_events,
            gap: 0,
            slope: 1.0,
            a_bin_span: (0, nbins - 1),
            b_bin_span: (0, nbins - 1),
            path,
            // Test fixture: z is not under test. NaN = not measured.
            z: f32::NAN,
        }
    }

    /// Walkthrough §4.5: three reads over a shared locus → a clean linear chain
    /// with every node and edge supported by 3 reads.
    fn overlap_shifted(a: u32, b: u32, a_off: u32, nbins: u32, bin_events: u32) -> Overlap {
        let path: Vec<(u32, u32)> = (0..nbins).map(|i| (a_off + i, i)).collect();
        Overlap {
            a,
            b,
            strand: Strand::Forward,
            n_anchors: nbins,
            matched: nbins * bin_events,
            cnt: nbins,
            aligned_bins: nbins,
            bin_events,
            gap: 0,
            slope: 1.0,
            a_bin_span: (a_off, a_off + nbins - 1),
            b_bin_span: (0, nbins - 1),
            path,
            // Test fixture: z is not under test. NaN = not measured.
            z: f32::NAN,
        }
    }

    #[test]
    fn three_reads_make_a_linear_chain() {
        let mut p = Params::r10_4_1().graph;
        p.node_bins = 3; // a node is 3 DP bins (72 events at bin_events 24)
        p.regovl_bins = 1;
        p.min_node_cov = 2;
        let be = 24;

        // a<b<c all cover 30 bins identically (slope 1).
        let ovls = vec![
            overlap(0, 1, 30, be),
            overlap(0, 2, 30, be),
            overlap(1, 2, 30, be),
        ];
        let g = build_graph(&ovls, &p);

        // 30 bins / K=3 = 10 windows → 10 nodes, a linear chain of 9 edges.
        assert_eq!(g.n_nodes(), 10, "expected 10 K-bin nodes");
        assert_eq!(g.n_edges(), 9, "expected a linear chain");
        assert!(
            g.nodes().iter().all(|n| n.cov == 3),
            "every node covered by 3 reads"
        );
        assert!(
            g.edges().iter().all(|e| e.cov == 3),
            "every edge covered by 3 reads"
        );
        // Linear: each node touched by at most 2 edges.
        let mut deg = vec![0u32; g.n_nodes()];
        for e in g.edges() {
            deg[e.from as usize] += 1;
            deg[e.to as usize] += 1;
        }
        assert!(deg.iter().all(|&d| d <= 2), "chain nodes have degree <= 2");
    }

    /// Orientation regression: edges must be oriented by **read order**, never by
    /// node id. `from`/`to` is what the GFA writer renders as the L-line
    /// orientation (`L x + y +` = x's right end joins y's left end), so on a linear
    /// chain every internal node must be the `to` of one edge and the `from` of the
    /// other — its two edges must use *opposite* ends.
    ///
    /// Canonicalizing on node id instead put both edges on the same end of an
    /// internal node whenever its id was the min or max of the three (p = 2/3 for
    /// random ids; 65% of degree-2 nodes on the `full_runs` assemblies), and negated
    /// `offset` for the flipped half, colliding with the "negative = overlap"
    /// contract that `gfa.rs` and `clean::step_len` both read.
    ///
    /// The fixture is a **staggered ladder**, not a uniform-coverage stack: node ids
    /// are handed out coverage-descending (Stage 3), so a chain of equal-coverage
    /// nodes gets ids in position order and canonicalizing on id is a silent no-op.
    /// The ladder's depth ramp (5,5,5,5,4,…,2) is what decorrelates id from position
    /// and exposes the fault — a uniform stack reports 0 same-end nodes even when the
    /// orientation is wrong, so it cannot catch this.
    #[test]
    fn chain_edges_are_oriented_by_read_order_not_node_id() {
        let mut p = Params::r10_4_1().graph;
        p.node_bins = 3;
        p.regovl_bins = 1;
        p.min_node_cov = 2;
        let be = 24;

        // Read i starts 6 bins after read i-1 and spans 30 bins; overlap every pair
        // that still shares ≥1 bin. Depth ramps up then down across the ladder.
        let n = 6usize;
        let mut ovls = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = ((j - i) * 6) as u32;
                if d < 30 {
                    ovls.push(overlap_shifted(i as u32, j as u32, d, 30 - d, be));
                }
            }
        }
        let g = build_graph(&ovls, &p);

        let mut indeg = vec![0u32; g.n_nodes()];
        let mut outdeg = vec![0u32; g.n_nodes()];
        for e in g.edges() {
            outdeg[e.from as usize] += 1;
            indeg[e.to as usize] += 1;
        }
        // Coverage-descending ids over a depth ramp ⇒ id order != chain order.
        assert!(
            g.nodes().windows(2).any(|w| w[0].cov != w[1].cov),
            "fixture must have non-uniform coverage or ids stay position-ordered \
             and the id canonicalization becomes an invisible no-op"
        );
        let internal: Vec<usize> = (0..g.n_nodes())
            .filter(|&i| indeg[i] + outdeg[i] == 2)
            .collect();
        assert_eq!(internal.len(), 14, "expected 14 internal chain nodes");
        for i in internal {
            assert_eq!(
                (indeg[i], outdeg[i]),
                (1, 1),
                "node {i} carries both edges on the same end (in={}, out={}) — edges \
                 got canonicalized on node id rather than read order",
                indeg[i],
                outdeg[i]
            );
        }
    }

    /// Both-grids regression: a single overlap generates candidate nodes on the
    /// query's grid and on the target's grid for the same loci. Conflict resolution
    /// must collapse them to one ~10-node linear chain, not ~20 nodes forming a
    /// phase-shifted double-tiling.
    #[test]
    fn both_grids_collapse_to_single_tiling() {
        let mut p = Params::r10_4_1().graph;
        p.node_bins = 3; // a node is 3 DP bins
        p.regovl_bins = 1;
        p.min_node_cov = 2;
        let be = 24;

        // One overlap, two reads fully aligned over 30 bins (slope 1).
        let g = build_graph(&[overlap(0, 1, 30, be)], &p);

        // 30 bins / K=3 = 10 windows → one tiling of 10 nodes, not two.
        assert_eq!(
            g.n_nodes(),
            10,
            "double-tiling not collapsed (got {} nodes)",
            g.n_nodes()
        );
        assert_eq!(g.n_edges(), 9, "expected a single linear chain");
        assert!(
            g.nodes().iter().all(|n| n.cov == 2),
            "every node covered by both reads"
        );
        // No junction: a single overlap must be linear (degree <= 2 everywhere).
        let mut deg = vec![0u32; g.n_nodes()];
        for e in g.edges() {
            deg[e.from as usize] += 1;
            deg[e.to as usize] += 1;
        }
        assert!(
            deg.iter().all(|&d| d <= 2),
            "collapsed tiling must be linear, degrees={deg:?}"
        );
    }

    #[test]
    fn no_overlaps_is_empty_graph() {
        let p = Params::r10_4_1().graph;
        let g = build_graph(&[], &p);
        assert_eq!(g.n_nodes(), 0);
        assert_eq!(g.n_edges(), 0);
    }

    /// The same overlaps under another read numbering (`perm[old] = new`): ids mapped, a pair whose
    /// lower id changes stored from the other read's side (spans swapped, cells transposed, slope
    /// inverted), sorted by `(a, b)`, as `ovl_relabel` writes them.
    fn relabelled(ovls: &[Overlap], perm: &[u32]) -> Vec<Overlap> {
        let mut out: Vec<Overlap> = ovls
            .iter()
            .map(|o| {
                let (a, b) = (perm[o.a as usize], perm[o.b as usize]);
                if a < b {
                    Overlap { a, b, ..o.clone() }
                } else {
                    let mut path: Vec<(u32, u32)> = o.path.iter().map(|&(x, y)| (y, x)).collect();
                    path.sort_unstable();
                    Overlap {
                        a: b,
                        b: a,
                        a_bin_span: o.b_bin_span,
                        b_bin_span: o.a_bin_span,
                        path,
                        slope: 1.0 / o.slope,
                        ..o.clone()
                    }
                }
            })
            .collect();
        out.sort_by_key(|o| (o.a, o.b));
        out
    }

    /// With the canonical tie-break the build does not depend on the read numbering: a staggered
    /// ladder (many equal-coverage candidates, so the tie-break decides the commit order) built under
    /// two numberings gives the same nodes, provenance and edges, node for node.
    #[test]
    fn ranked_build_is_independent_of_the_read_numbering() {
        let mut p = Params::r10_4_1().graph;
        p.node_bins = 3;
        p.regovl_bins = 1;
        p.min_node_cov = 2;
        let be = 24;
        let n = 6u32;
        let mut ovls = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                let off = 3 * (j - i);
                ovls.push(overlap_shifted(i, j, off, 36 - off, be));
            }
        }
        // read r's "name" is r: rank = identity under the original numbering
        let rank_x: Vec<u32> = (0..n).collect();
        let perm: Vec<u32> = vec![4, 0, 5, 2, 1, 3];
        let mut rank_y = vec![0u32; n as usize];
        for (old, &new) in perm.iter().enumerate() {
            rank_y[new as usize] = rank_x[old];
        }
        let (gx, ax) = build_graph_annotated_ranked(&ovls, &p, Some(&rank_x));
        let (gy, ay) = build_graph_annotated_ranked(&relabelled(&ovls, &perm), &p, Some(&rank_y));
        // the fixture has equal-coverage nodes, so the tie-break is exercised
        let mut covs: Vec<u32> = gx.nodes().iter().map(|x| x.cov).collect();
        covs.sort_unstable();
        assert!(
            covs.windows(2).any(|w| w[0] == w[1]),
            "no coverage ties: {covs:?}"
        );
        let node = |g: &FuzzyGraph| -> Vec<(u32, u32)> {
            g.nodes().iter().map(|x| (x.cov, x.len_events)).collect()
        };
        assert_eq!(node(&gx), node(&gy));
        let anno_x: Vec<(u32, u32, u32)> = ax
            .iter()
            .map(|x| (x.anchor_read, x.window, x.ev_beg))
            .collect();
        let anno_y: Vec<(u32, u32, u32)> = ay
            .iter()
            .map(|x| (rank_y[x.anchor_read as usize], x.window, x.ev_beg))
            .collect();
        assert_eq!(anno_x, anno_y);
        let edge = |g: &FuzzyGraph| -> Vec<(u32, u32, u32, i32)> {
            g.edges()
                .iter()
                .map(|e| (e.from, e.to, e.cov, e.offset))
                .collect()
        };
        assert_eq!(edge(&gx), edge(&gy));
        // and with no rank, the ranked build is the plain build
        let (g0, a0) = build_graph_annotated(&ovls, &p);
        let (g1, a1) = build_graph_annotated_ranked(&ovls, &p, None);
        assert_eq!((node(&g0), edge(&g0)), (node(&g1), edge(&g1)));
        assert_eq!(a0.len(), a1.len());
    }
}
