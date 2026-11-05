"""evaluation/table.py: QUAST fields, time totals (a signal-killed step counts as failed), components and Rawsamble's
protocol, per set and tool."""
import tempfile
import unittest

from evaltest import run, tsv_rows, write

REPORT = ("Assembly\t{l}\n# contigs\t3\nLargest contig\t900000\nNG50\t850000\nNGA50\t800000\n# misassemblies\t2\n"
          "Genome fraction (%)\t98.765\n")


class Table(unittest.TestCase):
    def test_rows(self):
        with tempfile.TemporaryDirectory() as d:
            r = f"{d}/k"
            write(f"{r}/ours/time.tsv", "detect\t10.50\t2097152\t100.00\t2.00\t0\nfinal\t1.00\t1048576\t1.00\t0.50\t0\n")
            write(f"{r}/flye.time.tsv", "flye\t20.40\t3145728\t30.00\t1.00\t0\n")
            write(f"{r}/wtdbg2.time.tsv", "wtdbg2\t5.00\t1048576\t9.00\t1.00\t0\nwtpoa_cns\t1.00\t524288\t2.00\t0.00\t0\n")
            # a miniasm crash: GNU time's signal line, then a row whose exit reads 0
            write(f"{r}/rawsamble.time.tsv", "rh_index\t3.00\t1048576\t5.00\t1.00\t0\nCommand terminated by signal 6\n"
                                             "miniasm\t1.00\t524288\t1.00\t0.00\t0\nminiasm_top25\t1.00\t524288\t1.00\t0.00\t0\n")
            write(f"{r}/rawsamble/miniasm_protocol.txt", "DEVIATION: miniasm needed the density filter TOPN=25\n")
            for label in ("ours_seq", "flye", "wtdbg2", "rawsamble"):
                write(f"{r}/score/quast_{label}/report.tsv", REPORT.format(l=label))
            write(f"{r}/score/components.tsv", "tool\tlevel\tgraph\tsegments\tcomponents\n"
                                               "ours\tlayout\tours/unitigs.gfa\t9\t2\nflye\tgraph\tflye/g.gfa\t4\t1\n")
            run("evaluation/table.py", d, "k", "--out", f"{d}/t.tsv")
            rows = {x["tool"]: x for x in tsv_rows(f"{d}/t.tsv")}
            self.assertEqual(list(rows), ["ours_seq", "flye", "wtdbg2", "wtd_lay", "rawsamble"])
            o = rows["ours_seq"]
            self.assertEqual((o["contigs"], o["largest"], o["NGA50"], o["genome_frac"]), ("3", "900000", "800000", "98.765"))
            self.assertEqual((o["wall_s"], o["cpu_s"], o["rss_gb"], o["failed_steps"], o["components"]),
                             ("12", "104", "2.00", "0", "2"))
            self.assertEqual(rows["flye"]["components"], "1")
            self.assertEqual((rows["wtd_lay"]["contigs"], rows["wtd_lay"]["wall_s"]), ("", "6"))
            rs = rows["rawsamble"]
            self.assertEqual((rs["failed_steps"], rs["rawsamble_protocol"], rs["components"]), ("1", "DEVIATION", ""))


if __name__ == "__main__":
    unittest.main()
