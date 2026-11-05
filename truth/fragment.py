#!/usr/bin/env python3
"""Cut reads into contiguous signal fragments of a target BASE length, at exact bases.

A cut every `target_bp * 12.5` samples, with the basecall sliced in proportion, would be
wrong twice over: the rate varies along a read, and some sets are sampled at 4 kHz, not
5 kHz. Here every cut is placed at the sample where base `k * L` begins according to the
parent's dorado move table, so a fragment's signal and its basecall slice describe exactly
the same bases.

A fragment keeps the parent's calibration and run metadata; only the signal window and
the read id differ: `uuid5(NAMESPACE_OID, "<parent>:L<L>:<k>")`. Fragment 0 keeps the
leading adapter samples and the last fragment keeps any signal beyond the basecalled
record, exactly as a shorter real read would.
A trailing fragment with fewer than `--min-frac * L` bases is dropped.

Writes:
    DST.pod5         the fragments
    DST.frag.tsv     frag_id, parent, k, sample_lo, sample_hi, base_lo, base_hi
    DST.fq.gz        fragment basecalls (read orientation)

event_ref.py --frag-map DST.frag.tsv then places each fragment's events through its
parent's move table and alignments: the truth is inherited, not re-derived.

Usage: fragment.py SRC.pod5 MOVES.bam DST --target-bp 4000 [--min-frac 0.6]
"""
import argparse
import gzip
import sys
import uuid

import numpy as np
import pod5
import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("moves")
    ap.add_argument("dst")
    ap.add_argument("--target-bp", type=int, required=True)
    ap.add_argument("--min-frac", type=float, default=0.6)
    a = ap.parse_args()
    L = a.target_bp

    mv = {}
    with pysam.AlignmentFile(a.moves, check_sq=False) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or not r.has_tag("mv"):
                continue
            m = np.frombuffer(r.get_tag("mv"), dtype=np.int8)
            mv[r.query_name] = (int(m[0]), m[1:].astype(np.int32), int(r.get_tag("ts")),
                                int(r.get_tag("sp")) if r.has_tag("sp") else 0,
                                r.get_forward_sequence(),
                                pysam.qualities_to_qualitystring(r.get_forward_qualities()))
    print(f"moves for {len(mv)} parents", file=sys.stderr)

    n_in = n_out = n_skip = 0
    bases_in = bases_out = 0
    with pod5.Reader(a.src) as rd, pod5.Writer(a.dst + ".pod5") as wr, \
            open(a.dst + ".frag.tsv", "w") as fm, gzip.open(a.dst + ".fq.gz", "wt") as fq:
        fm.write("frag\tparent\tk\tsample_lo\tsample_hi\tbase_lo\tbase_hi\n")
        for rec in rd.reads():
            pid = str(rec.read_id)
            n_in += 1
            if pid not in mv:
                n_skip += 1
                continue
            stride, flags, ts, sp, seq, qual = mv[pid]
            sig = rec.signal
            nb = len(seq)
            bases_in += nb
            # Sample at which each emitted base begins.
            emit = np.flatnonzero(flags == 1)
            base_sample = sp + ts + emit.astype(np.int64) * stride
            if nb < a.min_frac * L:
                continue
            b0, k = 0, 0
            while b0 < nb:
                b1 = min(b0 + L, nb)
                if nb - b1 < a.min_frac * L:
                    b1 = nb  # the remainder is too short to stand alone: absorb it
                s0 = 0 if k == 0 else int(base_sample[b0])
                s1 = len(sig) if b1 >= nb else int(base_sample[b1])
                fid = str(uuid.uuid5(uuid.NAMESPACE_OID, f"{pid}:L{L}:{k}"))
                wr.add_read(pod5.Read(
                    read_id=uuid.UUID(fid), pore=rec.pore, calibration=rec.calibration,
                    read_number=rec.read_number, start_sample=rec.start_sample + s0,
                    median_before=rec.median_before, end_reason=rec.end_reason,
                    run_info=rec.run_info, signal=sig[s0:s1]))
                fm.write(f"{fid}\t{pid}\t{k}\t{s0}\t{s1}\t{b0}\t{b1}\n")
                fq.write(f"@{fid}\n{seq[b0:b1]}\n+\n{qual[b0:b1]}\n")
                n_out += 1
                bases_out += b1 - b0
                b0, k = b1, k + 1
    print(f"{a.src} -> {a.dst}.pod5: {n_in} reads ({n_skip} without moves) -> {n_out} "
          f"fragments; {bases_out / max(bases_in, 1):.1%} of bases kept; mean fragment "
          f"{bases_out / max(n_out, 1):,.0f} bp", file=sys.stderr)


if __name__ == "__main__":
    main()
