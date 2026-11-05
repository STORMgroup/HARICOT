#!/usr/bin/env python3
"""EVALUATION ONLY: the sequence a layout asserts, read off its anchor reads (no consensus).

The layout (OUTDIR/unitig_paths.tsv, one row per path step) gives for each step the node's anchor read (the read
that defines the node), the node's retained events on it [ev_beg, ev_end), the direction the path crosses the node,
and advance_events: the node's length plus the edge's offset to the next node, i.e. where the next node starts, in
this node's anchor events (pipeline/src/fbg_build.rs, the unitig_paths.tsv writer). Each step writes its anchor
read's bases over exactly that stretch, and the unitig is the steps' pieces in path order:
  dir +   events [ev_beg, ev_beg + advance)
  dir -   events [ev_end - advance, ev_end), reverse-complemented
  the last step of a unitig: its whole node, [ev_beg, ev_end)
An event range [s, t) is the read's bases [B(s), B(t)), B(e) = the base of the first event at or after e that has
one (BASES: the basecaller's move table; an event with no base, adapter or past the last base, covers no bases),
the read's length if none. Nothing else: no choice of read, no check, no tolerance, no reference. Where a stretch
runs past its anchor read's end, the read's bases stop there and the rest is written as N, its length the events
past the end x the read's bases per event (counted). A step whose anchor read has no basecall or no base map is an
error: the read is to be basecalled, not stood in for.

Inputs: OUTDIR/unitig_paths.tsv; BASES (an evref.bin's base field: truth/truth.sh);
READS (FASTQ, gzip or plain: the basecalls). Unitigs shorter than --min-len (QUAST's --min-contig) are not written.
--pieces writes each step's piece (N runs are not pieces).

Usage: anchor_seq.py OUTDIR BASES READS OUT.fa [--pieces P.tsv] [--min-len 1000]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.evio import read_evbase, read_fastq, read_paths  # noqa: E402

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def base_at(bs, e, n_bases):
    """B(e): the base of the first event at or after e that has one; the read's length if none."""
    for x in bs[max(e, 0):]:
        if x >= 0:
            return int(x)
    return n_bases


def step_piece(bs, seq, beg, end, dirn, advance, last):
    """(event lo, event hi, base lo, base hi, N count) of one step on its anchor read."""
    n, L = len(bs), len(seq)
    if dirn == "+":
        lo, hi = beg, end if last else beg + advance
    else:
        lo, hi = (beg if last else end - advance), end
    over = max(hi - n, 0) + max(-lo, 0)
    lo, hi = max(lo, 0), min(hi, n)
    b_lo = base_at(bs, lo, L)
    b_hi = base_at(bs, hi, L) if hi < n else L
    return lo, hi, b_lo, b_hi, round(over * L / n) if over else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("bases")
    ap.add_argument("reads")
    ap.add_argument("out")
    ap.add_argument("--pieces", help="write each step's piece")
    ap.add_argument("--min-len", type=int, default=1000)
    a = ap.parse_args()

    paths = read_paths(a.outdir)
    need = {s[0] for st in paths.values() for s in st}
    bases = read_evbase(a.bases, need)
    seqs = read_fastq(a.reads, need)
    missing = sorted(need - (bases.keys() & seqs.keys()))
    if missing:
        sys.exit(f"anchor_seq: {len(missing)} anchor reads without a basecall or a base map (basecall them): "
                 f"{' '.join(missing[:20])}{' ...' if len(missing) > 20 else ''}")

    c = dict(unitigs=0, written=0, bp=0, steps=0, minus=0, over_steps=0, n_bp=0)  # steps: of the written unitigs
    pieces = []
    with open(a.out, "w") as o:
        for u, st in paths.items():
            c["unitigs"] += 1
            parts, rows, n_ns = [], [], []
            for i, (nm, beg, end, dirn, adv, _) in enumerate(st):
                lo, hi, b_lo, b_hi, n_n = step_piece(bases[nm], seqs[nm], beg, end, dirn, adv, i == len(st) - 1)
                s = seqs[nm][b_lo:b_hi]
                parts.append((s if dirn == "+" else s[::-1].translate(COMP)) + "N" * n_n)
                rows.append((u, nm, lo, hi, b_lo, b_hi, dirn, "anchor"))
                n_ns.append(n_n)
            seq = "".join(parts)
            if len(seq) < a.min_len:
                continue
            o.write(f">{u}\n")
            for i in range(0, len(seq), 100):
                o.write(seq[i:i + 100] + "\n")
            c["written"] += 1
            c["bp"] += len(seq)
            c["steps"] += len(st)
            c["minus"] += sum(s[3] == "-" for s in st)
            c["over_steps"] += sum(x > 0 for x in n_ns)
            c["n_bp"] += sum(n_ns)
            pieces += rows
    if a.pieces:
        with open(a.pieces, "w") as o:
            o.write("unitig\tread\tev_lo\tev_hi\tbase_lo\tbase_hi\tdir\tkind\n")
            for r in pieces:
                o.write("\t".join(map(str, r)) + "\n")
    print(f"anchor_seq: {c['unitigs']} unitigs, {c['written']} written ({c['bp'] / 1e6:.4f} Mb, >= {a.min_len} bp); "
          f"their {c['steps']:,} steps ({c['minus']:,} crossed '-'), {c['over_steps']:,} running past their anchor read "
          f"({c['n_bp']:,} bp written as N)", file=sys.stderr)


if __name__ == "__main__":
    main()
