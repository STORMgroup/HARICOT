"""evaluation/anchor_seq.py on a synthetic genome G: reads are exact pieces of G with 2 bases per retained event, so a
layout whose advances are right must come out as G itself (forward) or its reverse complement (the same layout
crossed '-'). Expected values are worked out in the comments."""
import os
import random
import struct
import tempfile
import unittest

from evaltest import run, tsv_rows, write

random.seed(1)
G = "".join(random.choice("ACGT") for _ in range(1200))
COMP = str.maketrans("ACGT", "TGCA")
HDR = "contig\tstep\tanchor_read\tev_beg\tev_end\tdir\tadvance_events\tcoord_events\n"


def rc(s):
    return s[::-1].translate(COMP)


# read -> (its bases, its retained events' bases): A = G[0:600], B = G[200:1000], 2 bases per event; A2 = A whose
# first three events have no base (adapter)
READS = {"A": (G[0:600], [2 * e for e in range(300)]),
         "B": (G[200:1000], [2 * e for e in range(400)]),
         "A2": (G[0:600], [-1, -1, -1] + [2 * e for e in range(3, 300)])}


def case(d, steps, reads=READS, *extra):
    """steps: (unitig, read, ev_beg, ev_end, dir, advance) rows -> {unitig: sequence}, pieces rows, the run."""
    write(f"{d}/run/unitig_paths.tsv", HDR + "".join(
        f"{u}\t{i}\t{r}\t{b}\t{e}\t{dr}\t{adv}\t0\n" for i, (u, r, b, e, dr, adv) in enumerate(steps)))
    with open(f"{d}/ev.bin", "wb") as fh:  # evbase.bin: u32 name_len, name, u32 n, n x i32 base
        for nm, (_, bs) in reads.items():
            fh.write(struct.pack("<I", len(nm)) + nm.encode() + struct.pack(f"<I{len(bs)}i", len(bs), *bs))
    write(f"{d}/reads.fq", "".join(f"@{nm} x\n{s}\n+\n{'I' * len(s)}\n" for nm, (s, _) in reads.items()))
    r = run("evaluation/anchor_seq.py", f"{d}/run", f"{d}/ev.bin", f"{d}/reads.fq", f"{d}/out.fa", "--pieces",
            f"{d}/p.tsv", *extra, check=False)
    out, cur = {}, None
    if os.path.exists(f"{d}/out.fa"):
        with open(f"{d}/out.fa") as fh:
            for line in fh:
                if line.startswith(">"):
                    cur = line[1:].strip()
                    out[cur] = ""
                else:
                    out[cur] += line.strip()
    return out, (tsv_rows(f"{d}/p.tsv") if os.path.exists(f"{d}/p.tsv") else []), r


class AnchorSeq(unittest.TestCase):
    def test_forward_path_is_the_genome(self):
        # A [0,100) advance 150: A events [0,150) = A bases [0,300) = G[0:300]
        # B [50,150) advance 100: B events [50,150) = B bases [100,300) = G[300:500]
        # B [150,250) last: the whole node, B bases [300,500) = G[500:700]
        with tempfile.TemporaryDirectory() as d:
            out, p, _ = case(d, [("u1", "A", 0, 100, "+", 150), ("u1", "B", 50, 150, "+", 100),
                                 ("u1", "B", 150, 250, "+", 100)], READS, "--min-len", "0")
            self.assertEqual(out["u1"], G[0:700])
            self.assertEqual([(x["read"], x["ev_lo"], x["ev_hi"], x["base_lo"], x["base_hi"]) for x in p],
                             [("A", "0", "150", "0", "300"), ("B", "50", "150", "100", "300"),
                              ("B", "150", "250", "300", "500")])

    def test_reverse_path_is_its_reverse_complement(self):
        # the same nodes crossed '-': each step's advance is measured from its node's right end
        # B [150,250) advance 100: events [150,250) = G[500:700], reversed
        # B [50,150) advance 150: events [0,150) = B bases [0,300) = G[200:500], reversed
        # A [0,100) last: G[0:200], reversed
        with tempfile.TemporaryDirectory() as d:
            out, _, _ = case(d, [("u2", "B", 150, 250, "-", 100), ("u2", "B", 50, 150, "-", 150),
                                 ("u2", "A", 0, 100, "-", 100)], READS, "--min-len", "0")
            self.assertEqual(out["u2"], rc(G[0:700]))

    def test_stretch_past_the_read_end_is_n(self):
        # A [200,300) advance 150: events [200,350), A has 300: A bases [400,600) = G[400:600], then 50 events past
        # the end x A's 600 / 300 bases per event = 100 N; B [250,350) last: B bases [500,700) = G[700:900]
        with tempfile.TemporaryDirectory() as d:
            out, _, r = case(d, [("u3", "A", 200, 300, "+", 150), ("u3", "B", 250, 350, "+", 100)], READS,
                             "--min-len", "0")
            self.assertEqual(out["u3"], G[400:600] + "N" * 100 + G[700:900])
            self.assertIn("1 running past their anchor read (100 bp written as N)", r.stderr)

    def test_events_without_a_base_cover_none(self):
        # A2's events 0-2 have no base: B(0) = the base of event 3 = 6, so the node [0,100) is G[6:200]
        with tempfile.TemporaryDirectory() as d:
            out, _, _ = case(d, [("u4", "A2", 0, 100, "+", 100)], READS, "--min-len", "0")
            self.assertEqual(out["u4"], G[6:200])

    def test_missing_read_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            out, _, r = case(d, [("u5", "A", 0, 100, "+", 100), ("u5", "C", 0, 100, "+", 100)], READS,
                             "--min-len", "0")
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("1 anchor reads without a basecall or a base map", r.stderr)
            self.assertIn("C", r.stderr)

    def test_short_unitigs_are_not_written(self):
        # u1 is 700 bp: below the default --min-len 1000
        with tempfile.TemporaryDirectory() as d:
            out, p, _ = case(d, [("u1", "A", 0, 100, "+", 150), ("u1", "B", 50, 150, "+", 100),
                                 ("u1", "B", 150, 250, "+", 100)])
            self.assertEqual(out, {})
            self.assertEqual(p, [])


if __name__ == "__main__":
    unittest.main()
