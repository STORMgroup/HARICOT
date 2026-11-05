#!/usr/bin/env python3
"""EVALUATION ONLY: wtdbg2's own contig layout written WITHOUT its consensus (the consensus control).

`asm.ctg.lay.gz` (wtdbg2 -o asm) lists per contig its edges in order -- `E <offset> <node1> <dir1> <node2>
<dir2>`, offset in contig bp -- and under each edge the read segments that span it -- `S <read> <dir> <beg>
<len> <seq>`, the sequence already in the contig's orientation (what wtpoa-cns builds its consensus from).
Here each edge contributes ONE read segment, the first listed, cut to the distance to the next edge's
offset; the last edge contributes its whole segment: one read per step cut at the layout's offset, no consensus,
as anchor_seq.py writes our layouts (wtdbg2 lists several reads per edge, not one anchor per node, so the
first listed is taken; its offsets sit on wtdbg2's 1 kb grid). Scoring this and wtdbg2's consensus (asm.fa) with
identical QUAST settings measures what QUAST charges for the missing consensus on a layout we did not build.

Usage: wtdbg2_lay.py asm.ctg.lay.gz OUT.fa [--min-len 1000] [--pick first|longest]
"""
import argparse
import gzip
import sys


def contigs(path):
    name, edges, cur = None, [], None
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    if cur:
                        edges.append(cur)
                    yield name, edges
                name, edges, cur = line[1:].split()[0], [], None
            elif line.startswith("E\t"):
                if cur:
                    edges.append(cur)
                cur = [int(line.split("\t")[1]), []]
            elif line.startswith("S\t") and cur is not None:
                f = line.rstrip("\n").split("\t")
                cur[1].append(f[5] if len(f) > 5 else "")
    if name is not None:
        if cur:
            edges.append(cur)
        yield name, edges


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lay")
    ap.add_argument("out")
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--pick", choices=("first", "longest"), default="first")
    a = ap.parse_args()
    n = w = tot = empty = 0
    with open(a.out, "w") as o:
        for name, edges in contigs(a.lay):
            n += 1
            parts = []
            for k, (off, segs) in enumerate(edges):
                segs = [s for s in segs if s]
                if not segs:
                    empty += 1
                    continue
                s = segs[0] if a.pick == "first" else max(segs, key=len)
                if k + 1 < len(edges):
                    s = s[:max(edges[k + 1][0] - off, 0)]
                parts.append(s)
            seq = "".join(parts)
            if len(seq) < a.min_len:
                continue
            o.write(f">{name}\n")
            for i in range(0, len(seq), 100):
                o.write(seq[i:i + 100] + "\n")
            w += 1
            tot += len(seq)
    print(f"wtdbg2_lay: {n} contigs, {w} written ({tot / 1e6:.2f} Mb); edges without a read segment {empty}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
