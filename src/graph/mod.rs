//! Assembly graph: the fuzzy Bruijn graph and its simplification passes.

pub mod build;
pub mod clean;
pub mod ends;
pub mod fbg;
pub mod threading;

pub use build::{
    NodeAnno, build_graph, build_graph_annotated, build_graph_annotated_ranked, build_graph_dumping,
};
pub use clean::clean;
pub use ends::{End, EndAdj};
pub use fbg::{Dir, Edge, FuzzyGraph, Node};
