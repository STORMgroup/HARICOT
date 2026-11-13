//! Layout & output: unitigs/contigs and GFA serialization.

pub mod frg;
pub mod gfa;
pub mod unitig;

pub use frg::{FrgParams, FrgStats, frg_contigs};
pub use gfa::{write_gfa, write_header, write_link, write_segment};
pub use unitig::{Contig, Unitig, UnitigLink, gen_contigs_from, gen_unitigs, unitig_links};
