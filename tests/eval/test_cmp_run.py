"""checks/cmp_run.sh REF NEW [NREADS]: byte identity of a recorded run, what it skips (logs, timings, scorer
outputs, the invariance check's canonical dumps), the content-hash appendix of a candidate cache, missing and added files,
usage errors."""
import os
import shutil
import struct
import subprocess
import tempfile
import unittest

from evaltest import ROOT, write

CMP = os.path.join(ROOT, "checks", "cmp_run.sh")


def cmp_run(ref, new, *nreads):
    r = subprocess.run(["bash", CMP, ref, new, *map(str, nreads)], capture_output=True, text=True)
    return r.returncode, r.stdout


class CmpRun(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.mkdtemp()
        self.ref, self.new = f"{self.t}/ref", f"{self.t}/new"
        for name, text in (("unitigs.gfa", "S\tu0\t*\tLN:i:5\n"), ("work/unitig_paths.tsv", "contig\tstep\n"),
                           ("work/detect.log", "a\n"), ("time.tsv", "detect\t1\t2\t3\t4\t0\n"),
                           ("metrics.tsv", "label\nk\n"), ("work/canon/ovl.tsv", "x\n")):
            write(f"{self.ref}/{name}", text)
        with open(f"{self.ref}/work/ovl.bin", "wb") as fh:
            fh.write(b"CANDIDATES")
        shutil.copytree(self.ref, self.new)

    def tearDown(self):
        shutil.rmtree(self.t)

    def test_identical_and_skipped(self):
        # logs, time.tsv, metrics.tsv and work/canon/ differ: not compared
        for name in ("work/detect.log", "time.tsv", "metrics.tsv", "work/canon/ovl.tsv"):
            write(f"{self.new}/{name}", "other\n")
        rc, out = cmp_run(self.ref, self.new)
        self.assertEqual(rc, 0, out)
        self.assertIn("compared 3 files: 3 same, 0 same + content-hash appendix, 0 differ, 0 missing", out)
        self.assertIn("added (in NEW only; logs not listed): none", out)

    def test_differs_missing_added(self):
        write(f"{self.new}/unitigs.gfa", "S\tu0\t*\tLN:i:6\n")
        os.remove(f"{self.new}/work/unitig_paths.tsv")
        write(f"{self.new}/extra.tsv", "x\n")
        rc, out = cmp_run(self.ref, self.new)
        self.assertEqual(rc, 1)
        self.assertIn("  DIFFERS: unitigs.gfa", out)
        self.assertIn("  missing: work/unitig_paths.tsv", out)
        self.assertIn("added (in NEW only; logs not listed): extra.tsv", out)
        self.assertIn("compared 2 files: 1 same, 0 same + content-hash appendix, 1 differ, 1 missing", out)

    def test_content_hash_appendix(self):
        # NEW's cache = REF's bytes + a bincode Vec<u64> of 2 content hashes (8-byte length + 8 bytes each)
        with open(f"{self.new}/work/ovl.bin", "ab") as fh:
            fh.write(struct.pack("<QQQ", 2, 11, 22))
        rc, out = cmp_run(self.ref, self.new, 2)
        self.assertEqual(rc, 0, out)
        self.assertIn("compared 3 files: 2 same, 1 same + content-hash appendix, 0 differ, 0 missing", out)
        rc, out = cmp_run(self.ref, self.new, 3)  # another read count: the appendix does not fit
        self.assertEqual(rc, 1)
        self.assertIn("  DIFFERS: work/ovl.bin", out)
        rc, out = cmp_run(self.ref, self.new)  # without NREADS a longer cache differs
        self.assertEqual(rc, 1)

    def test_usage_errors(self):
        r = subprocess.run(["bash", CMP, "-h"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 2)
        self.assertIn("cmp_run.sh REF NEW [NREADS]", r.stdout)
        r = subprocess.run(["bash", CMP, f"{self.t}/absent", self.new], capture_output=True, text=True)
        self.assertEqual(r.returncode, 2)  # REF is not a directory
        os.makedirs(f"{self.t}/empty")
        rc, out = cmp_run(f"{self.t}/empty", self.new)
        self.assertEqual(rc, 1)  # nothing compared is not a pass
        self.assertIn("compared 0 files", out)


if __name__ == "__main__":
    unittest.main()
