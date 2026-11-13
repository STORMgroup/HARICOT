//! Read placements: path-level evidence the FBG builds alongside its edges.
//!
//! `build_graph_annotated` threads each read through the committed regions it belongs to
//! — a per-read *path* through the node graph — and then reduces every such path to its
//! consecutive pairs, keeping only an edge coverage count. The placements are kept for the
//! fragment layout ([`crate::layout::frg_contigs`]), which links fragments by the reads that
//! pass through them.

use rustc_hash::FxHashMap;

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
