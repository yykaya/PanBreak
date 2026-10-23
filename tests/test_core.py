import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from panbreak.annotate import ltr_length, period, tsd
from panbreak.bubbles import Graph, Reference, walk
from panbreak.common import canonical, rc
from panbreak.genotype import call_gt
from panbreak.kmers import junction_kmers
from panbreak.popgen import bh, fisher_greater, hudson, wilson

OPTS = {"ploidy": 2, "error": 0.02, "min_alt": 2, "min_ref": 2, "min_depth": 2, "min_gq": 10}


def toy_gfa(path):
    nodes = {"1": "ACGTACGTAC" * 3, "2": "TTGACCATGG" * 12, "3": "GGCATTACCA" * 3, "4": "CATCATCAGG" * 15, "5": "AACCGGTTAA" * 3}
    with open(path, "w") as fh:
        fh.write("H\tVN:Z:1.0\n")
        for n, s in nodes.items():
            fh.write(f"S\t{n}\t{s}\n")
        fh.write("P\tREF#1#REF_x\t1+,2+,3+,5+\t*\n")
        fh.write("P\tDEL#1#DEL_x\t1+,3+,5+\t*\n")
        fh.write("P\tINS#1#INS_x\t1+,2+,3+,4+,5+\t*\n")
        fh.write("P\tDUP#1#DUP_x\t1+,2+,3+,2+,3+,5+\t*\n")
        fh.write("P\tINS2#1#INS2_x\t1+,2+,3+,4+,5+\t*\n")
    return nodes


class TestSequence(unittest.TestCase):
    def test_rc_canonical(self):
        self.assertEqual(rc("AACG"), "CGTT")
        self.assertEqual(canonical("TTTT"), "AAAA")

    def test_junction_kmers_span_breakpoint(self):
        w = "A" * 30 + "C" * 30
        ks = junction_kmers(w, 30, 31, 8)
        self.assertEqual(len(ks), 16)
        for km, off in ks:
            self.assertTrue(8 <= off <= 31 - 8)
            self.assertEqual(km[:off], "A" * off)

    def test_period_and_tsd(self):
        self.assertEqual(period("AT" * 30), 2)
        self.assertIsNone(period("ACGTTGCAAGCTTACGGATCCA"))
        left, te, right = "GGCATTACCTGATTA", "TG" + "CATGCATTTACGG" * 10 + "CA", "CCGTAGGTAC"
        ctx = left + te + "GATTA" + right
        self.assertEqual(tsd(ctx, len(left), len(left) + len(te)), "GATTA")

    def test_ltr_length(self):
        import random
        rng = random.Random(1)
        ltr = "".join(rng.choice("ACGT") for _ in range(250))
        body = "".join(rng.choice("ACGT") for _ in range(2000))
        self.assertTrue(240 <= ltr_length(ltr + body + ltr) <= 260)


class TestGraph(unittest.TestCase):
    def test_walk_events(self):
        with tempfile.TemporaryDirectory() as d:
            toy_gfa(os.path.join(d, "t.gfa"))
            g = Graph(os.path.join(d, "t.gfa"))
            ref = Reference(g, "REF")
            ev = {p.split("#")[0]: walk(g, ref, p, 50)[1] for p in g.paths}
            self.assertEqual(ev["REF"], [])
            self.assertEqual([(e["type"], e["a"], e["b"]) for e in ev["DEL"]], [("DEL", 30, 150)])
            self.assertEqual([(e["type"], e["alt_len"]) for e in ev["INS"]], [("INS", 150)])
            self.assertEqual([e["type"] for e in ev["DUP"]], ["DUP"])
            from panbreak.bubbles import cluster_events, kmer_set
            evs = []
            for p in g.paths:
                s, es, span, bm = walk(g, ref, p, 50)
                for e in es:
                    e["kmers"] = kmer_set(s[e["i0"]:e["i1"]])
                    evs.append(e)
            cl = cluster_events(evs, 50, 0.8, 0.3)
            self.assertEqual(sorted((c["type"], len(c["members"])) for c in cl), [("DEL", 1), ("DUP", 1), ("INS", 2)])


class TestModel(unittest.TestCase):
    def test_genotypes(self):
        self.assertEqual(call_gt(0, 20, "biallelic", 1 / 3, OPTS)["GT"], "0/0")
        self.assertEqual(call_gt(20, 0, "biallelic", 1 / 3, OPTS)["GT"], "1/1")
        self.assertEqual(call_gt(7, 13, "biallelic", 1 / 3, OPTS)["GT"], "0/1")
        self.assertEqual(call_gt(1, 0, "biallelic", 1 / 3, OPTS)["call"], "weak")
        self.assertEqual(call_gt(0, 1, "biallelic", 1 / 3, OPTS)["call"], "nodata")
        p = call_gt(8, 16, "presence", 0.5, OPTS)
        self.assertEqual(p["carrier"], 1)
        self.assertAlmostEqual(p["copy_estimate"], 2.0)

    def test_haploid(self):
        o = dict(OPTS, ploidy=1)
        self.assertEqual(call_gt(15, 0, "biallelic", 0.5, o)["GT"], "1")


class TestPopgen(unittest.TestCase):
    def test_stats(self):
        lo, hi = wilson(5, 10)
        self.assertTrue(lo < 0.5 < hi)
        self.assertLess(fisher_greater(9, 1, 1, 9), 0.01)
        self.assertEqual(bh([0.01, 0.04, 0.5]), [0.03, 0.06, 0.5])
        num, den = hudson(0.9, 40, 0.1, 40)
        self.assertGreater(num / den, 0.6)


if __name__ == "__main__":
    unittest.main()
