//! GFA v1 output for the assembly graph.
//!
//! Segments carry no nucleotide sequence (this is a raw-signal assembler), so
//! the sequence field is `*` and the segment length is reported in *events*
//! via the standard `LN:i` tag plus a `dp:i` depth tag. Edge overlaps are
//! derived from the node offset. Visualize with Bandage / gfatools.

use std::io::{self, Write};

use crate::graph::FuzzyGraph;
use crate::layout::stringgraph::StringGraph;

/// Segment-name prefix so node ids are valid GFA identifiers.
const SEG_PREFIX: &str = "n";

/// Write the GFA v1 header line (`H\tVN:Z:1.0`).
pub fn write_header<W: Write>(out: &mut W) -> io::Result<()> {
    writeln!(out, "H\tVN:Z:1.0")
}

/// Write a segment line: `S <name> * LN:i:<len_events>`, plus a `dp:i:<depth>`
/// tag when `depth` is given. Sequence is always `*` (raw-signal assembler).
pub fn write_segment<W: Write>(
    out: &mut W,
    name: &str,
    len_events: u32,
    depth: Option<u32>,
) -> io::Result<()> {
    match depth {
        Some(d) => writeln!(out, "S\t{name}\t*\tLN:i:{len_events}\tdp:i:{d}"),
        None => writeln!(out, "S\t{name}\t*\tLN:i:{len_events}"),
    }
}

/// Write a link line: `L <from> <fd> <to> <td> <overlap>M`, plus an `ec:i:<cov>`
/// tag when `edge_cov` is given and an `es:f:<support>` tag when `edge_sup` is.
/// `fd`/`td` are orientation signs (`+`/`-`).
#[allow(clippy::too_many_arguments)] // one GFA L-line: its fields, and two optional tags
pub fn write_link<W: Write>(
    out: &mut W,
    from: &str,
    from_dir: char,
    to: &str,
    to_dir: char,
    overlap: u32,
    edge_cov: Option<u32>,
    edge_sup: Option<f32>,
) -> io::Result<()> {
    let sup = match edge_sup {
        Some(v) => format!("\tes:f:{v:.4}"),
        None => String::new(),
    };
    match edge_cov {
        Some(c) => writeln!(
            out,
            "L\t{from}\t{from_dir}\t{to}\t{to_dir}\t{overlap}M\tec:i:{c}{sup}"
        ),
        None => writeln!(
            out,
            "L\t{from}\t{from_dir}\t{to}\t{to_dir}\t{overlap}M{sup}"
        ),
    }
}

/// Write `graph` as GFA v1 to `out`.
pub fn write_gfa<W: Write>(graph: &FuzzyGraph, out: &mut W) -> io::Result<()> {
    write_header(out)?;
    for node in graph.nodes() {
        if node.del {
            continue; // tombstoned by a cleaning pass
        }
        write_segment(
            out,
            &format!("{SEG_PREFIX}{}", node.id),
            node.len_events,
            Some(node.cov),
        )?;
    }
    for e in graph.edges() {
        if graph.nodes()[e.from as usize].del || graph.nodes()[e.to as usize].del {
            continue; // edge incident on a tombstoned node
        }
        // Negative offset = overlap; clamp gaps to 0M.
        let overlap = if e.offset < 0 { (-e.offset) as u32 } else { 0 };
        // ec:i: edge coverage (supporting reads) — handy for Bandage / diagnostics.
        write_link(
            out,
            &format!("{SEG_PREFIX}{}", e.from),
            e.from_dir.sign(),
            &format!("{SEG_PREFIX}{}", e.to),
            e.to_dir.sign(),
            overlap,
            Some(e.cov),
            // es:f: weighted edge support — equals `ec` unless evidence weighting is on.
            Some(e.support),
        )?;
    }
    Ok(())
}

/// Convenience: render the graph to a `String`.
pub fn to_gfa_string(graph: &FuzzyGraph) -> String {
    let mut buf = Vec::new();
    // Writing to a Vec<u8> is infallible.
    write_gfa(graph, &mut buf).expect("writing GFA to an in-memory buffer cannot fail");
    String::from_utf8(buf).expect("GFA output is valid UTF-8")
}

/// Write an OLC read [`StringGraph`] as GFA v1: one segment per read (length in
/// events), one link per reduced dovetail edge, plus containment attachments.
pub fn write_string_graph_gfa<W: Write>(
    graph: &StringGraph,
    names: &[String],
    len_events: &[u32],
    out: &mut W,
) -> io::Result<()> {
    write_header(out)?;
    for (rid, name) in names.iter().enumerate() {
        write_segment(out, name, len_events[rid], None)?;
    }
    for e in &graph.edges {
        write_link(
            out,
            &names[e.from as usize],
            e.from_dir.sign(),
            &names[e.to as usize],
            e.to_dir.sign(),
            e.overlap_events,
            None,
            None,
        )?;
    }
    // Containment attachments: container(+) -> child(+), overlap = child length.
    for &(container, child) in &graph.contain_edges {
        write_link(
            out,
            &names[container as usize],
            '+',
            &names[child as usize],
            '+',
            len_events[child as usize],
            None,
            None,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::fbg::{Dir, Edge, FuzzyGraph};

    #[test]
    fn emits_header_segments_and_links() {
        let mut g = FuzzyGraph::new();
        let a = g.add_node(3, 80);
        let b = g.add_node(2, 80);
        g.add_edge(Edge {
            from: a,
            from_dir: Dir::Fwd,
            to: b,
            to_dir: Dir::Rev,
            cov: 2,
            support: 2.0,
            offset: -20,
        });
        let gfa = to_gfa_string(&g);
        assert!(gfa.starts_with("H\tVN:Z:1.0\n"));
        assert!(gfa.contains("S\tn0\t*\tLN:i:80\tdp:i:3"));
        assert!(gfa.contains("L\tn0\t+\tn1\t-\t20M"));
    }

    #[test]
    fn empty_graph_is_just_header() {
        let g = FuzzyGraph::new();
        assert_eq!(to_gfa_string(&g), "H\tVN:Z:1.0\n");
    }
}
