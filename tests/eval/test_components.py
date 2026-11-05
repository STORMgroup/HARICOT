"""evaluation/components.py: connected components over every segment and link of a GFA or a wtdbg2 .dot graph."""
import gzip
import os
import tempfile
import unittest

from evaltest import run, write


class Components(unittest.TestCase):
    def test_gfa(self):
        with tempfile.TemporaryDirectory() as d:
            # a-b-c joined (one link written in both orientations), d alone, e only named by a link to f
            write(f"{d}/g.gfa", "H\tVN:Z:1.0\nS\ta\t*\tLN:i:5\nS\tb\tACGT\nS\tc\t*\nS\td\t*\n"
                                "L\ta\t+\tb\t+\t0M\nL\tb\t-\ta\t-\t0M\nL\tb\t+\tc\t-\t0M\nL\te\t+\tf\t+\t0M\n")
            out = run("evaluation/components.py", f"{d}/g.gfa").stdout.split("\t")
            self.assertEqual([int(x) for x in out[1:]], [6, 3])  # a b c d e f; {a,b,c} {d} {e,f}

    def test_dot_gz(self):
        with tempfile.TemporaryDirectory() as d:
            text = ('digraph {\nF1 [label="F1"]\nF2 [label="F2"]\nF3 [label="F3"]\n'
                    'F1 -> F2 [label="++ 3"]\nF2 -> F1 [label="-- 3"]\n}\n')
            with gzip.open(f"{d}/asm.ctg.dot.gz", "wt") as fh:
                fh.write(text)
            out = run("evaluation/components.py", f"{d}/asm.ctg.dot.gz").stdout.split("\t")
            self.assertEqual([int(x) for x in out[1:]], [3, 2])

    def test_empty_graph(self):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/g.gfa", "H\tVN:Z:1.0\n")
            self.assertEqual(run("evaluation/components.py", f"{d}/g.gfa").stdout.rstrip("\n").split("\t")[1:],
                             ["0", "0"])
            self.assertTrue(os.path.exists(f"{d}/g.gfa"))


if __name__ == "__main__":
    unittest.main()
