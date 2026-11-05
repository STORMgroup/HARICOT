//! Fuzzy Bruijn Graph (FBG) data model.
//!
//! Adapted from wtdbg2 (`extern/wtdbg2/wtdbg.h`, `wtdbg-graph.h`). Unlike a
//! classical de Bruijn graph of fixed-length k-mers, a node here is a
//! coverage-aggregated read **segment** of ~`node_events` events; edges record
//! the observed offset between consecutive segments along reads. That tolerance
//! for imprecise segment boundaries is what wtdbg2 relies on for noisy long reads.
//!
//! The data model and basic graph operations live here; construction from
//! overlaps is in [`crate::graph::build`].

/// Orientation of a node as traversed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dir {
    Fwd,
    Rev,
}

impl Dir {
    /// `+`/`-` as used in GFA.
    pub fn sign(self) -> char {
        match self {
            Dir::Fwd => '+',
            Dir::Rev => '-',
        }
    }
}

/// A graph node: a coverage-aggregated read segment.
#[derive(Debug, Clone)]
pub struct Node {
    /// Index of this node within the graph.
    pub id: u32,
    /// Number of reads supporting this node.
    pub cov: u32,
    /// Support as the relative edge gate (`clean::filter_weak_edges`) reads it: `cov as
    /// f32` as built; `threading::unzip` gives each copy of a split node its share.
    pub support: f32,
    /// Nominal length of the segment, in events.
    pub len_events: u32,
    /// Tombstone: set by cleaning passes that "remove" a node. Node ids must stay
    /// stable across cleaning (provenance / GFA names reference them), so removal
    /// is a soft delete — the slot survives but is skipped by [`FuzzyGraph::live_nodes`],
    /// the GFA writer, and traversal.
    pub del: bool,
}

/// A directed edge between two oriented nodes.
#[derive(Debug, Clone)]
pub struct Edge {
    pub from: u32,
    pub from_dir: Dir,
    pub to: u32,
    pub to_dir: Dir,
    /// Number of reads supporting this edge.
    pub cov: u32,
    /// Support; see [`Node::support`]. Written to the GFA as `es:f`.
    pub support: f32,
    /// Observed offset between segment boundaries (events); may be negative
    /// (overlap) or positive (gap).
    pub offset: i32,
}

/// The fuzzy Bruijn graph.
#[derive(Debug, Default, Clone)]
pub struct FuzzyGraph {
    nodes: Vec<Node>,
    edges: Vec<Edge>,
    /// Read-threaded node triples; empty unless `GraphParams::read_threading` is on.
    /// See [`crate::graph::threading`].
    transitions: crate::graph::threading::Transitions,
    /// Where each read sits in the graph; recorded under the same flag.
    read_threads: crate::graph::threading::ReadThreads,
}

impl FuzzyGraph {
    pub fn new() -> Self {
        FuzzyGraph::default()
    }

    /// Add a node, returning its id.
    pub fn add_node(&mut self, cov: u32, len_events: u32) -> u32 {
        let id = self.nodes.len() as u32;
        self.nodes.push(Node {
            id,
            cov,
            support: cov as f32,
            len_events,
            del: false,
        });
        id
    }

    /// Read-threaded node triples (empty when threading was not recorded).
    pub fn transitions(&self) -> &crate::graph::threading::Transitions {
        &self.transitions
    }

    /// Install the read-threaded triples.
    pub fn set_transitions(&mut self, tr: crate::graph::threading::Transitions) {
        self.transitions = tr;
    }

    /// Where each read was threaded (empty when threading was not recorded).
    pub fn read_threads(&self) -> &crate::graph::threading::ReadThreads {
        &self.read_threads
    }

    /// Install the per-read placements.
    pub fn set_read_threads(&mut self, rt: crate::graph::threading::ReadThreads) {
        self.read_threads = rt;
    }

    /// Override a node's weighted support (see [`Node::support`]).
    pub fn set_node_support(&mut self, id: u32, support: f32) {
        if let Some(n) = self.nodes.get_mut(id as usize) {
            n.support = support;
        }
    }

    /// Add a directed edge.
    pub fn add_edge(&mut self, edge: Edge) {
        self.edges.push(edge);
    }

    /// Total node slots (including tombstoned ones — ids index into this).
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Total edges (including any incident on tombstoned nodes — call
    /// [`prune_dangling_edges`](Self::prune_dangling_edges) to drop those).
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Number of non-tombstoned nodes.
    pub fn n_live_nodes(&self) -> usize {
        self.nodes.iter().filter(|n| !n.del).count()
    }

    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// Is node `id` present (not tombstoned)?
    pub fn is_live(&self, id: u32) -> bool {
        self.nodes.get(id as usize).map(|n| !n.del).unwrap_or(false)
    }

    /// Iterate non-tombstoned nodes.
    pub fn live_nodes(&self) -> impl Iterator<Item = &Node> {
        self.nodes.iter().filter(|n| !n.del)
    }

    /// Tombstone a node (soft delete; keeps ids stable). Incident edges remain in
    /// the vec until [`prune_dangling_edges`](Self::prune_dangling_edges) runs but
    /// are ignored by [`live_adjacency`](Self::live_adjacency) and the GFA writer.
    pub fn remove_node(&mut self, id: u32) {
        if let Some(n) = self.nodes.get_mut(id as usize) {
            n.del = true;
        }
    }

    /// Drop edges that touch a tombstoned node. Returns the number removed.
    pub fn prune_dangling_edges(&mut self) -> usize {
        let before = self.edges.len();
        let nodes = &self.nodes;
        self.edges
            .retain(|e| !nodes[e.from as usize].del && !nodes[e.to as usize].del);
        before - self.edges.len()
    }

    /// Mutable access to edges (used by graph-cleaning passes).
    pub(crate) fn edges_mut(&mut self) -> &mut Vec<Edge> {
        &mut self.edges
    }

    /// Mutable access to nodes (used by the node-merge cleaning pass to fold
    /// coverage into a surviving representative).
    pub(crate) fn nodes_mut(&mut self) -> &mut Vec<Node> {
        &mut self.nodes
    }

    /// Build an undirected adjacency list over **live** nodes/edges. Each entry is
    /// `(neighbor, edge_index)` so callers can read cov/offset from `edges()`.
    ///
    /// This discards which **end** of a node each edge attaches to, so it answers
    /// only genuinely undirected questions — connectivity, "how many edges does
    /// this node have at all". Traversals that must follow real sequence use
    /// [`EndAdj`](crate::graph::ends::EndAdj) instead.
    pub fn live_adjacency(&self) -> Vec<Vec<(u32, usize)>> {
        let mut adj = vec![Vec::new(); self.nodes.len()];
        for (ei, e) in self.edges.iter().enumerate() {
            if self.nodes[e.from as usize].del || self.nodes[e.to as usize].del {
                continue;
            }
            adj[e.from as usize].push((e.to, ei));
            adj[e.to as usize].push((e.from, ei));
        }
        adj
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_and_inspect_graph() {
        let mut g = FuzzyGraph::new();
        let a = g.add_node(3, 80);
        let b = g.add_node(2, 80);
        g.add_edge(Edge {
            from: a,
            from_dir: Dir::Fwd,
            to: b,
            to_dir: Dir::Fwd,
            cov: 2,
            support: 2.0,
            offset: -20,
        });
        assert_eq!(g.n_nodes(), 2);
        assert_eq!(g.n_edges(), 1);
        assert_eq!(g.nodes()[0].id, 0);
        assert_eq!(g.edges()[0].offset, -20);
    }
}
