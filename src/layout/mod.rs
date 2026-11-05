//! Layout & output: unitigs/contigs and GFA serialization.

pub mod frg;
pub mod gfa;
pub mod stringgraph;
pub mod unitig;

pub use frg::{FrgParams, FrgStats, LinkMode, frg_contigs};
pub use gfa::{
    to_gfa_string, write_gfa, write_header, write_link, write_segment, write_string_graph_gfa,
};
pub use stringgraph::{StringGraph, StringGraphParams, build_string_graph};
pub use unitig::{
    Contig, Unitig, UnitigLink, gen_contigs, gen_contigs_from, gen_unitigs, gen_unitigs_threaded,
    unitig_links,
};
