//! Node-**end** adjacency: the orientation-aware view of the graph.
//!
//! A node is a segment with two physical ends, L and R. An edge attaches to one
//! specific end of each node it touches, exactly as in GFA: `L a + b +` joins a's
//! **R** end to b's **L** end.
//!
//! **Why this module exists.** [`FuzzyGraph::live_adjacency`] collapses that down
//! to "n is adjacent to m", discarding which end.
//! [`build_graph_annotated`](crate::graph::build::build_graph_annotated) orients
//! every edge by read order, so the ends carry real information and the undirected
//! view is lossy in a way that silently invents structure:
//!
//! * **Undirected degree 2 does not mean "chain interior".** A node whose two
//!   edges land on the *same* end is a **convergence** — two things arriving at
//!   one side, with nothing leaving the other. Walking "in one edge and out the
//!   other" through it traverses a path that does not physically exist. Measured
//!   on one full E. coli fwd run, 175 of 14,123 degree-2 nodes (1.2%) were
//!   convergences.
//! * The same applies to any *path* search: a genuine path enters a node at one
//!   end and leaves by the **other**.
//!
//! Use [`EndAdj`] wherever a traversal, chain-walk or path search is meant to
//! follow real sequence. Genuine undirected notions — connected components, "is
//! this node isolated" — can keep using [`FuzzyGraph::live_adjacency`].

use rustc_hash::FxHashMap;

use crate::graph::FuzzyGraph;
use crate::graph::fbg::Dir;

/// Which physical end of a node an edge attaches to.
///
/// `Ord` is derived only so `(node, end)` can key an ordered container / priority
/// queue; L-before-R carries no meaning.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum End {
    L,
    R,
}

impl End {
    /// The node's other end.
    pub fn flip(self) -> End {
        match self {
            End::L => End::R,
            End::R => End::L,
        }
    }
}

/// The end of an edge's `from` node that the edge leaves by (GFA: `L a + ...`
/// leaves a's R end).
pub fn from_end(d: Dir) -> End {
    match d {
        Dir::Fwd => End::R,
        Dir::Rev => End::L,
    }
}

/// The end of an edge's `to` node that the edge arrives at (GFA: `L ... b +`
/// enters b's L end).
pub fn to_end(d: Dir) -> End {
    match d {
        Dir::Fwd => End::L,
        Dir::Rev => End::R,
    }
}

/// Traversing a node oriented `d`: the end we leave by.
pub fn exit_end(d: Dir) -> End {
    match d {
        Dir::Fwd => End::R,
        Dir::Rev => End::L,
    }
}

/// Traversing a node oriented `d`: the end we entered by.
pub fn entry_end(d: Dir) -> End {
    match d {
        Dir::Fwd => End::L,
        Dir::Rev => End::R,
    }
}

/// Per-node-end adjacency over live nodes/edges:
/// `(node, end) -> [(other_node, other_end, edge_index)]`.
///
/// Self-loops are excluded: they cannot participate in a linear chain, and
/// admitting one would make both of a node's ends look singly-linked, merging the
/// node into itself.
/// `(node, end)` -> the `(other node, its end, edge index)` of every edge at that end.
type EndMap = FxHashMap<(u32, End), Vec<(u32, End, usize)>>;

#[derive(Debug, Default, Clone)]
pub struct EndAdj {
    map: EndMap,
}

impl EndAdj {
    pub fn new(graph: &FuzzyGraph) -> Self {
        let mut map: EndMap = FxHashMap::default();
        for (ei, e) in graph.edges().iter().enumerate() {
            if !graph.is_live(e.from) || !graph.is_live(e.to) || e.from == e.to {
                continue;
            }
            let (ae, be) = (from_end(e.from_dir), to_end(e.to_dir));
            map.entry((e.from, ae)).or_default().push((e.to, be, ei));
            map.entry((e.to, be)).or_default().push((e.from, ae, ei));
        }
        EndAdj { map }
    }

    /// Edges attached to one end of a node.
    pub fn at(&self, n: u32, end: End) -> &[(u32, End, usize)] {
        self.map.get(&(n, end)).map(|v| &v[..]).unwrap_or(&[])
    }

    /// Number of edges on one end.
    pub fn deg(&self, n: u32, end: End) -> usize {
        self.at(n, end).len()
    }

    /// The unique edge continuing off `(n, end)` — present only when that end
    /// carries exactly one edge **and** the far end does too. Anything else is a
    /// branch, and terminates a chain.
    pub fn simple_link(&self, n: u32, end: End) -> Option<(u32, End, usize)> {
        let v = self.at(n, end);
        if v.len() != 1 {
            return None;
        }
        let (m, em, ei) = v[0];
        if self.deg(m, em) != 1 {
            return None;
        }
        Some((m, em, ei))
    }

    /// Step forward out of node `n` oriented `d`: the next `(node, orientation,
    /// edge)` when the chain continues without branching.
    pub fn next(&self, n: u32, d: Dir) -> Option<(u32, Dir, usize)> {
        let (m, em, ei) = self.simple_link(n, exit_end(d))?;
        // We arrive at m by `em`, so m is traversed Fwd iff we entered its L end.
        Some((m, if em == End::L { Dir::Fwd } else { Dir::Rev }, ei))
    }

    /// Step out of `(n, d)` by the end we exit, requiring only that **our** side is
    /// unbranched — the destination may be a junction.
    ///
    /// Differs from [`next`](Self::next), which additionally demands the far end be
    /// unbranched. Use this to walk *to* a branch point and learn its identity: at
    /// the node where two paths reconverge the entry end has degree ≥ 2, so `next`
    /// stops one step short and cannot name it.
    /// [`pop_bubbles`](crate::graph::clean::pop_bubbles) needs the name.
    pub fn step_out(&self, n: u32, d: Dir) -> Option<(u32, Dir, usize)> {
        let v = self.at(n, exit_end(d));
        if v.len() != 1 {
            return None;
        }
        let (m, em, ei) = v[0];
        Some((m, if em == End::L { Dir::Fwd } else { Dir::Rev }, ei))
    }

    /// Step backward into node `n` oriented `d`.
    pub fn prev(&self, n: u32, d: Dir) -> Option<(u32, Dir, usize)> {
        let (m, em, ei) = self.simple_link(n, entry_end(d))?;
        // m hands off to us from `em`, i.e. m *exits* by em ⇒ Fwd iff that is its R end.
        Some((m, if em == End::R { Dir::Fwd } else { Dir::Rev }, ei))
    }

    /// Is `n` a true chain interior — exactly one edge on **each** end?
    ///
    /// This is the orientation-aware replacement for `live_adjacency`'s
    /// `deg(n) == 2`, which also accepts a convergence (both edges on one end).
    pub fn is_chain_interior(&self, n: u32) -> bool {
        self.deg(n, End::L) == 1 && self.deg(n, End::R) == 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::fbg::{Edge, FuzzyGraph};

    fn e(from: u32, fd: Dir, to: u32, td: Dir) -> Edge {
        Edge {
            from,
            from_dir: fd,
            to,
            to_dir: td,
            cov: 3,
            support: 3.0,
            offset: 0,
        }
    }

    /// The distinction the whole module exists for: A->B->C makes B a chain
    /// interior; A->B<-C gives B the same *undirected* degree 2 but both edges
    /// land on B's L end, so it is a convergence and must not be walked through.
    #[test]
    fn convergence_is_not_a_chain_interior() {
        let mut chain = FuzzyGraph::new();
        for _ in 0..3 {
            chain.add_node(3, 72);
        }
        chain.add_edge(e(0, Dir::Fwd, 1, Dir::Fwd)); // 0R -> 1L
        chain.add_edge(e(1, Dir::Fwd, 2, Dir::Fwd)); // 1R -> 2L
        let a = EndAdj::new(&chain);
        assert!(
            a.is_chain_interior(1),
            "A->B->C: B is a real chain interior"
        );
        assert_eq!((a.deg(1, End::L), a.deg(1, End::R)), (1, 1));

        let mut conv = FuzzyGraph::new();
        for _ in 0..3 {
            conv.add_node(3, 72);
        }
        conv.add_edge(e(0, Dir::Fwd, 1, Dir::Fwd)); // 0R -> 1L
        conv.add_edge(e(2, Dir::Fwd, 1, Dir::Fwd)); // 2R -> 1L  (same end of 1)
        let a = EndAdj::new(&conv);
        assert_eq!(
            conv.live_adjacency()[1].len(),
            2,
            "undirected degree is 2 either way — which is why it misleads"
        );
        assert!(
            !a.is_chain_interior(1),
            "A->B<-C: both edges on B's L end, so B is a convergence"
        );
        assert_eq!((a.deg(1, End::L), a.deg(1, End::R)), (2, 0));
        assert!(a.next(1, Dir::Fwd).is_none(), "nothing leaves B's R end");
    }

    /// `next`/`prev` are inverses along a chain, and orientation propagates.
    #[test]
    fn next_and_prev_traverse_a_chain() {
        let mut g = FuzzyGraph::new();
        for _ in 0..3 {
            g.add_node(3, 72);
        }
        g.add_edge(e(0, Dir::Fwd, 1, Dir::Fwd));
        g.add_edge(e(1, Dir::Fwd, 2, Dir::Fwd));
        let a = EndAdj::new(&g);
        let (n, d, _) = a.next(0, Dir::Fwd).expect("0 -> 1");
        assert_eq!((n, d), (1, Dir::Fwd));
        let (n, d, _) = a.next(n, d).expect("1 -> 2");
        assert_eq!((n, d), (2, Dir::Fwd));
        assert!(a.next(n, d).is_none(), "2's R end is free");
        let (n, d, _) = a.prev(2, Dir::Fwd).expect("2 <- 1");
        assert_eq!((n, d), (1, Dir::Fwd));
    }

    /// A node reached by a `-` link is entered at its R end and traversed Rev.
    #[test]
    fn reverse_link_flips_orientation() {
        let mut g = FuzzyGraph::new();
        for _ in 0..2 {
            g.add_node(3, 72);
        }
        // 0's R end joins 1's R end: 1 is entered at R, so it is traversed Rev.
        g.add_edge(e(0, Dir::Fwd, 1, Dir::Rev));
        let a = EndAdj::new(&g);
        let (n, d, _) = a.next(0, Dir::Fwd).expect("0 -> 1");
        assert_eq!((n, d), (1, Dir::Rev));
        assert_eq!((a.deg(1, End::R), a.deg(1, End::L)), (1, 0));
    }
}
