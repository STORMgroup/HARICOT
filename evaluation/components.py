#!/usr/bin/env python3
"""Connected components of an assembly graph: every segment (no length floor) and every link.

Graphs: a GFA (S and L lines: HARICOT's unitigs.gfa, Flye's assembly_graph.gfa, Rawsamble's asm.gfa) or wtdbg2's
fragment graph (asm.ctg.dot or .dot.gz: `F<n> [label=...` nodes, `F<a> -> F<b>` edges).

Usage: components.py GRAPH.{gfa,dot,dot.gz} [GRAPH ...]   -> one line per graph: path<TAB>segments<TAB>components
"""
import gzip
import re
import sys


def read_gfa(path):
    segs, links = set(), set()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                segs.add(f[1])
            elif f[0] == "L":
                links.add((f[1], f[3]))
    return segs, links


def read_dot(path):
    segs, links = set(), set()
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            m = re.match(r"\s*(F\d+) \[label=", line)
            if m:
                segs.add(m.group(1))
                continue
            m = re.match(r"\s*(F\d+) -> (F\d+) \[label=\"([+-])([+-])", line)
            if m:
                links.add((m.group(1), m.group(2)))
    return segs, links


def components(segs, links):
    """(segments, components): the segments with every link endpoint, grouped by the links (union-find)."""
    nodes = set(segs) | {x for l in links for x in l}
    parent = {u: u for u in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in links:
        parent[find(a)] = find(b)
    return len(nodes), len({find(u) for u in nodes})


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for path in sys.argv[1:]:
        dot = path.endswith(".dot") or path.endswith(".dot.gz")
        segs, links = read_dot(path) if dot else read_gfa(path)
        n, c = components(segs, links)
        print(f"{path}\t{n}\t{c}")


if __name__ == "__main__":
    main()
