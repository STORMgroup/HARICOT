#!/usr/bin/env python3
"""Exact event -> read base -> reference placement, from dorado move tables and TWO aligners.

A uniform-rate model places a retained event at `e / n_events * rlen` read bases: one
rate per read. That ignores translocation-speed variation along the read, the
adapter samples before dorado's trim point, and the fact that segmentation density
tracks the sequence. Here each event is placed by the basecaller's own move table:

    event span (raw samples, haricot-truth-events)  ->  mid sample s
    block = (s - sp - ts) // stride                 (dorado `mv`, `ts`, `sp`)
    read base = (# emitted bases up to that block) - 1
    reference = CIGAR walk of the alignment covering that base

MEMBERSHIP comes from two aligners with different seeding, minimap2 and winnowmap
(repeat-aware), each keeping primary AND supplementary records, so a chimera's second
locus is placed rather than silently clipped. Per event:

    1 CONF      both aligners place it within --tol bp, and either record has mapq >= 10
    2 SINGLE    only one aligner places it, at mapq >= 10
    3 LOWQ      placed, but no record reaching mapq 10 (an unresolved repeat)
    4 CONFLICT  both place it, somewhere different
    0 UNPLACED  adapter, clip, reverse-strand fold-back of a forward read, or unmapped

Outputs (prefix OUT):
    OUT.evref.bin    per read: u32 name_len, name, u32 n, n x i64 ref, n x i32 base,
                     n x i32 seg, n x u8 status   (ref = global coord, -1 unplaced)
    OUT.segs.tsv     one row per placed alignment segment
    OUT.anchored.tsv one row per read: its primary placement (name, read0, rlen, aln_start,
                     aln_end, clips, mapq, nm, strand, status fractions, segment count),
                     from the primary minimap2 record (winnowmap if none)
    OUT.summary.txt  validation: moves vs pod5 sample counts, and how wrong the uniform
                     rate model is, per read

Global coordinates: contigs in header order, each offset by the previous ones' length
plus --gap.

Usage: event_ref.py SPANS.bin MOVES.bam MM2.bam WM.bam --out PREFIX [--frag-map FRAG.tsv] [--tol 100]
                    [--min-mapq 10] [--gap 10000000] [--fwd-only]
"""
import argparse
import struct
import sys
from collections import defaultdict

import numpy as np
import pysam

UNPLACED, CONF, SINGLE, LOWQ, CONFLICT = 0, 1, 2, 3, 4

# CIGAR op codes (pysam): M I D N S H P = X
Q_ADV = np.array([1, 1, 0, 0, 1, 1, 0, 1, 1], dtype=bool)   # consumes query (H: full read)
R_ADV = np.array([1, 0, 1, 1, 0, 0, 0, 1, 1], dtype=bool)
ALN = np.array([1, 0, 0, 0, 0, 0, 0, 1, 1], dtype=bool)     # base aligned to a ref base
INS = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0], dtype=bool)


def read_spans(path):
    buf = open(path, "rb").read()
    off, out = 0, {}
    while off < len(buf):
        nl = int.from_bytes(buf[off:off + 4], "little")
        off += 4
        name = buf[off:off + nl].decode()
        off += nl
        ns, n = struct.unpack_from("<II", buf, off)
        off += 8
        st = np.frombuffer(buf, dtype="<u4", count=n, offset=off)
        off += 4 * n
        en = np.frombuffer(buf, dtype="<u4", count=n, offset=off)
        off += 4 * n
        out[name] = (ns, st, en)
    return out


def contig_offsets(bam, gap):
    off, cur = {}, 0
    for sq in bam.header.to_dict().get("SQ", []):
        off[sq["SN"]] = cur
        cur += sq["LN"] + gap
    return off


def base_to_ref(rec, off):
    """(read-orientation base -> global ref) for one record; -1 where not aligned.

    Insertion bases take the reference position of the next aligned base, so the map is
    monotone along the alignment. Length = the FULL read, hard clips included."""
    ops = np.array(rec.cigartuples, dtype=np.int64)
    code, n = ops[:, 0], ops[:, 1]
    qadv = np.where(Q_ADV[code], n, 0)
    radv = np.where(R_ADV[code], n, 0)
    ref_at_op = rec.reference_start + np.cumsum(radv) - radv
    tot = int(qadv.sum())
    opi = np.repeat(np.arange(len(code)), qadv)
    within = np.arange(tot) - np.repeat(np.cumsum(qadv) - qadv, qadv)
    c = code[opi]
    ref = np.where(ALN[c], ref_at_op[opi] + within, np.where(INS[c], ref_at_op[opi], -1))
    ref = np.where(ref >= 0, ref + off[rec.reference_name], -1)
    if rec.is_reverse:
        ref = ref[::-1]
    return ref


def load_alignments(path, want, off_of):
    """name -> list of (is_primary, mapq, strand, contig, base->ref array)."""
    recs = defaultdict(list)
    with pysam.AlignmentFile(path, check_sq=False) as bam:
        off = off_of(bam)
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.query_name not in want:
                continue
            recs[r.query_name].append(
                (not r.is_supplementary, r.mapping_quality, "-" if r.is_reverse else "+",
                 r.reference_name, base_to_ref(r, off), r))
    return recs


def combine(recs, L):
    """Per-base (ref, record index) with primary first, then by mapq."""
    ref = np.full(L, -1, dtype=np.int64)
    rid = np.full(L, -1, dtype=np.int32)
    order = sorted(range(len(recs)), key=lambda i: (not recs[i][0], -recs[i][1]))
    for i in order:
        m = recs[i][4]
        if len(m) != L:
            continue  # length mismatch: a different basecall; reported in the summary
        take = (m >= 0) & (ref < 0)
        ref[take] = m[take]
        rid[take] = i
    return ref, rid


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spans")
    ap.add_argument("moves")
    ap.add_argument("mm2")
    ap.add_argument("wm")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=int, default=100, help="aligner agreement, bp")
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--gap", type=int, default=10_000_000)
    ap.add_argument("--frag-map", help="fragment.py DST.frag.tsv: place each fragment "
                    "through its PARENT's move table and alignments")
    ap.add_argument("--fwd-only", action="store_true",
                    help="treat reverse-strand records as unplaced (fold-backs of a "
                         "forward read set); leave off when a sample-matched reference "
                         "may place whole reads on the minus strand")
    a = ap.parse_args()

    spans = read_spans(a.spans)
    alias = {}  # fragment -> (parent, sample offset, base_lo, base_hi)
    if a.frag_map:
        for line in open(a.frag_map):
            f = line.rstrip("\n").split("\t")
            if f[0] == "frag":
                continue
            alias[f[0]] = (f[1], int(f[3]), int(f[5]), int(f[6]))
    want = {alias[n][0] if n in alias else n for n in spans}
    print(f"spans: {len(spans)} reads, {sum(len(v[1]) for v in spans.values()):,} events",
          file=sys.stderr)

    moves = {}
    with pysam.AlignmentFile(a.moves, check_sq=False) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.query_name not in want:
                continue
            if not r.has_tag("mv"):
                continue
            mv = np.frombuffer(r.get_tag("mv"), dtype=np.int8)
            moves[r.query_name] = (
                int(mv[0]), mv[1:].astype(np.int32), int(r.get_tag("ts")),
                int(r.get_tag("ns")), int(r.get_tag("sp")) if r.has_tag("sp") else 0,
                r.query_length)
    print(f"moves: {len(moves)} reads", file=sys.stderr)

    offs = {}

    def off_of(bam):
        o = contig_offsets(bam, a.gap)
        offs.update(o)
        return o

    mm = load_alignments(a.mm2, want, off_of)
    wm = load_alignments(a.wm, want, off_of)
    print(f"alignments: mm2 {len(mm)} reads, wm {len(wm)} reads", file=sys.stderr)

    fbin = open(a.out + ".evref.bin", "wb")
    fseg = open(a.out + ".segs.tsv", "w")
    fseg.write("seg\tread\taligner\tprimary\tcontig\tstrand\tmapq\tev_lo\tev_hi\t"
               "ref_lo\tref_hi\tn_conf\n")
    fanc = open(a.out + ".anchored.tsv", "w")
    fanc.write("name\tread0\trlen\taln_start\taln_end\tclip5\tclip3\tmapq\tnm\tstrand\t"
               "frac_conf\tfrac_placed\tn_segs\n")

    seg_id = 0
    tally = np.zeros(5, dtype=np.int64)
    ns_mismatch, len_mismatch, no_moves = 0, 0, 0
    uni_err = []          # |move-table base - uniform-rate base| per event, sampled
    per_read_p90 = []
    rng = np.random.default_rng(0)
    for name, (n_samp, st, en) in spans.items():
        n = len(st)
        ref = np.full(n, -1, dtype=np.int64)
        base = np.full(n, -1, dtype=np.int32)
        seg = np.full(n, -1, dtype=np.int32)
        status = np.zeros(n, dtype=np.uint8)
        parent, s0 = (alias[name][0], alias[name][1]) if name in alias else (name, 0)
        mvr = moves.get(parent)
        if mvr is None or n == 0:
            no_moves += mvr is None
        else:
            stride, flags, ts, ns, sp, L = mvr
            if ns != n_samp and name not in alias:
                ns_mismatch += 1
            cum = np.cumsum(flags)
            mid = (st.astype(np.int64) + en.astype(np.int64)) // 2 + s0
            blk = (mid - sp - ts) // stride
            ok = (blk >= 0) & (blk < len(cum))
            b = np.full(n, -1, dtype=np.int64)
            b[ok] = cum[blk[ok]] - 1
            b[(b < 0) | (b >= L)] = -1
            base[:] = b
            # How wrong the uniform `e / n * rlen` model is.
            good = b >= 0
            if good.sum() > 10 and name not in alias:
                uni = np.arange(n) / n * L
                d = np.abs(uni[good] - b[good])
                per_read_p90.append(np.percentile(d, 90))
                uni_err.append(rng.choice(d, size=min(200, len(d)), replace=False))

            ra, rb = mm.get(parent, []), wm.get(parent, [])
            if a.fwd_only:
                ra = [x for x in ra if x[2] == "+"]
                rb = [x for x in rb if x[2] == "+"]
            if any(len(x[4]) != L for x in ra + rb):
                len_mismatch += 1
            ref_a, id_a = combine(ra, L)
            ref_b, id_b = combine(rb, L)
            bb = np.where(good, b, 0)
            ea = np.where(good, ref_a[bb], -1)
            eb = np.where(good, ref_b[bb], -1)
            ia = np.where(good, id_a[bb], -1)
            ib = np.where(good, id_b[bb], -1)
            q_a = np.array([x[1] for x in ra] + [0])[ia]
            q_b = np.array([x[1] for x in rb] + [0])[ib]
            s_a = np.array([x[2] for x in ra] + ["."])[ia]
            s_b = np.array([x[2] for x in rb] + ["."])[ib]
            both = (ea >= 0) & (eb >= 0)
            agree = both & (np.abs(ea - eb) <= a.tol) & (s_a == s_b)
            hiq = np.maximum(np.where(ea >= 0, q_a, 0), np.where(eb >= 0, q_b, 0)) >= a.min_mapq
            only_a = (ea >= 0) & (eb < 0)
            only_b = (eb >= 0) & (ea < 0)
            status[agree & hiq] = CONF
            status[(only_a & (q_a >= a.min_mapq)) | (only_b & (q_b >= a.min_mapq))] = SINGLE
            status[(agree | only_a | only_b) & ~hiq] = LOWQ
            status[both & ~agree] = CONFLICT
            use_a = ea >= 0
            ref[:] = np.where(use_a, ea, eb)
            ref[status == CONFLICT] = -1
            ref[status == UNPLACED] = -1
            # Segments: runs of one (aligner, record) over placed events.
            key = np.where(use_a, ia, 10_000 + ib)
            pidx = np.flatnonzero(ref >= 0)
            if len(pidx):
                kp = key[pidx]
                cuts = np.flatnonzero(np.diff(kp)) + 1
                for lo, hi in zip(np.r_[0, cuts], np.r_[cuts, len(pidx)]):
                    ev = pidx[lo:hi]
                    k = int(kp[lo])
                    rec = ra[k] if k < 10_000 else rb[k - 10_000]
                    seg[ev] = seg_id
                    rr = ref[ev]
                    fseg.write(f"{seg_id}\t{name}\t{'mm2' if k < 10_000 else 'wm'}\t"
                               f"{int(rec[0])}\t{rec[3]}\t{rec[2]}\t{rec[1]}\t{ev[0]}\t"
                               f"{ev[-1] + 1}\t{int(rr.min())}\t{int(rr.max()) + 1}\t"
                               f"{int((status[ev] == CONF).sum())}\n")
                    seg_id += 1

            # anchored.tsv from the primary minimap2 record (winnowmap if none).
            prim = [x for x in ra if x[0]] or [x for x in rb if x[0]]
            if name in alias:
                pl = ref[ref >= 0]
                if len(pl):
                    fl = alias[name][3] - alias[name][2]
                    q = prim[0][1] if prim else 0
                    sd = prim[0][2] if prim else "+"
                    fanc.write(f"{name}\t{int(pl.min())}\t{fl}\t{int(pl.min())}\t"
                               f"{int(pl.max()) + 1}\t0\t0\t{q}\t-1\t{sd}\t"
                               f"{(status == CONF).mean():.4f}\t{(ref >= 0).mean():.4f}\t"
                               f"{len(set(seg[seg >= 0].tolist()))}\n")
            elif prim:
                p = prim[0][5]
                o = offs[p.reference_name]
                ops = p.cigartuples
                c5 = ops[0][1] if ops[0][0] in (4, 5) else 0
                c3 = ops[-1][1] if ops[-1][0] in (4, 5) else 0
                if p.is_reverse:
                    c5, c3 = c3, c5
                s0, e0 = p.reference_start + o, p.reference_end + o
                nm = p.get_tag("NM") if p.has_tag("NM") else -1
                fanc.write(f"{name}\t{s0 - c5 if not p.is_reverse else e0 + c5}\t{L}\t{s0}\t"
                           f"{e0}\t{c5}\t{c3}\t{p.mapping_quality}\t{nm}\t"
                           f"{'-' if p.is_reverse else '+'}\t"
                           f"{(status == CONF).mean():.4f}\t{(ref >= 0).mean():.4f}\t"
                           f"{len(set(seg[seg >= 0].tolist()))}\n")
        tally += np.bincount(status, minlength=5)
        fbin.write(struct.pack("<I", len(name)))
        fbin.write(name.encode())
        fbin.write(struct.pack("<I", n))
        fbin.write(ref.astype("<i8").tobytes())
        fbin.write(base.astype("<i4").tobytes())
        fbin.write(seg.astype("<i4").tobytes())
        fbin.write(status.tobytes())
    fbin.close()
    fseg.close()
    fanc.close()

    tot = tally.sum()
    lines = [
        f"reads {len(spans)}  events {tot:,}  no_moves {no_moves}  "
        f"ns!=pod5 {ns_mismatch}  cigar_len!=seq {len_mismatch}",
        "event status: " + "  ".join(
            f"{k} {v / max(tot, 1):.2%}" for k, v in
            zip(("unplaced", "conf", "single", "lowq", "conflict"), tally)),
        f"segments {seg_id}",
    ]
    if uni_err:
        d = np.concatenate(uni_err)
        pr = np.array(per_read_p90)
        lines.append(
            "uniform-rate model error vs move table (bases): "
            f"p50 {np.percentile(d, 50):.0f}  p90 {np.percentile(d, 90):.0f}  "
            f"p99 {np.percentile(d, 99):.0f}  max {d.max():.0f}; "
            f"reads whose p90 error > 1 bin (267 bp): {(pr > 267).mean():.2%}")
    open(a.out + ".summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines), file=sys.stderr)


if __name__ == "__main__":
    main()
