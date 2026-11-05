//! E-mer domain: quantization of events and sketching into hashable seeds.
//!
//! This is the bridge between raw-signal processing ([`crate::signal`]) and the
//! assembly graph ([`crate::overlap`], [`crate::graph`]): events become
//! quantized e-mers, then seeds.

pub mod quantize;
pub mod select;
pub mod sketch;

pub use quantize::{dynamic_quantize, normal_quantile, normal_quantile_edges};
pub use select::{SelectCfg, SketchMethod, select};
pub use sketch::{QEvent, Seed, hash64, retained_quantized, sketch, sketch_from_retained};
