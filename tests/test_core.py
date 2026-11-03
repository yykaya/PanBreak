import csv
import io
import os
import random
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

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


def random_locus(d):
    rng = random.Random(7)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    nodes = {"1": seq(300), "2": seq(400), "3": seq(300), "4": seq(350), "5": seq(300)}
    paths = {"REF": "1 2 3 5", "DEL": "1 3 5", "INS": "1 2 3 4 5", "INS2": "1 2 3 4 5"}
    gfa = os.path.join(d, "locus.gfa")
    with open(gfa, "w") as fh:
        fh.write("H\tVN:Z:1.0\n")
        for n, s in nodes.items():
            fh.write(f"S\t{n}\t{s}\n")
        for g, p in paths.items():
            fh.write(f"P\t{g}#1#{g}_x\t{','.join(x + '+' for x in p.split())}\t*\n")
    return gfa, {g: "".join(nodes[x] for x in p.split()) for g, p in paths.items()}


def write_reads(path, sample, haplotypes, length=100, step=4):
    import pysam
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": 10000}], "RG": [{"ID": sample, "SM": sample}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        i = 0
        for s in haplotypes:
            for st in range(0, len(s) - length, step):
                r = pysam.AlignedSegment(out.header)
                r.query_name = f"r{i}"
                r.query_sequence = s[st:st + length] if i % 2 else rc(s[st:st + length])
                r.flag = 4
                r.query_qualities = pysam.qualitystring_to_array("I" * length)
                r.set_tag("RG", sample)
                out.write(r)
                i += 1


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


class TestGFAFormats(unittest.TestCase):
    def nodes(self):
        rng = random.Random(3)
        seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        return {"1": seq(300), "2": seq(400), "3": seq(300), "4": seq(350), "5": seq(300)}

    def events(self, gfa, ref):
        g = Graph(gfa)
        r = Reference(g, ref)
        return r, {p: [(e["type"], e["a"], e["b"], e["alt_len"]) for e in walk(g, r, p, 50)[1]] for p in g.paths}

    def test_walk_lines(self):
        from panbreak.gfa import normalise
        n = self.nodes()
        with tempfile.TemporaryDirectory() as d:
            src, out = os.path.join(d, "mc.gfa"), os.path.join(d, "locus.gfa")
            with open(src, "w") as fh:
                for k, s in n.items():
                    fh.write(f"S\t{k}\t{s}\n")
                fh.write("W\tREF\t0\tchr1\t5000\t6300\t>1>2>3>5\n")
                fh.write("W\tDEL\t1\tctg1\t0\t900\t>1>3>5\n")
                fh.write("W\tINS\t2\tctg2\t0\t1650\t<5<4<3<2<1\n")
            normalise(src, out)
            r, ev = self.events(out, "REF")
            self.assertEqual(r.path, "REF#0#chr1:5000-6300")
            self.assertEqual(r.offset, 5000)
            self.assertEqual(ev["DEL#1#ctg1"], [("DEL", 5300, 5700, 0)])
            self.assertEqual(ev["INS#2#ctg2"], [("INS", 6000, 6000, 350)])

    def test_rgfa_gaf(self):
        from panbreak.gfa import normalise
        n = self.nodes()
        sn = {"1": ("REF#1#REF_x", 0, 0), "2": ("REF#1#REF_x", 300, 0), "3": ("REF#1#REF_x", 700, 0),
              "5": ("REF#1#REF_x", 1000, 0), "4": ("INS#1#INS_x", 300, 1)}
        with tempfile.TemporaryDirectory() as d:
            src, gaf, out = os.path.join(d, "g.rgfa"), os.path.join(d, "a.gaf"), os.path.join(d, "locus.gfa")
            with open(src, "w") as fh:
                for k, s in n.items():
                    name, so, sr = sn[k]
                    fh.write(f"S\ts{k}\t{s}\tLN:i:{len(s)}\tSN:Z:{name}\tSO:i:{so}\tSR:i:{sr}\n")
                for a, b in [(1, 2), (2, 3), (3, 5), (1, 3), (3, 4), (4, 5)]:
                    fh.write(f"L\ts{a}\t+\ts{b}\t+\t0M\n")
            with open(gaf, "w") as fh:
                fh.write("DEL#1#DEL_x\t900\t0\t900\t+\t>s1>s3>s5\t900\t0\t900\t900\t900\t60\ttp:A:P\n")
                fh.write("INS#1#INS_x\t1650\t0\t1650\t+\t>REF#1#REF_x:0-1000>INS#1#INS_x:300-650>REF#1#REF_x:1000-1300\t1650\t0\t1650\t1650\t1650\t60\ttp:A:P\n")
                fh.write("INS#1#INS_x\t1650\t0\t300\t+\t>s1\t300\t0\t300\t300\t300\t5\ttp:A:S\n")
            normalise(src, out, [gaf])
            r, ev = self.events(out, "REF")
            self.assertEqual(r.path, "REF#1#REF_x")
            self.assertEqual(sorted(ev), ["DEL#1#DEL_x", "INS#1#INS_x", "REF#1#REF_x"])
            self.assertEqual(ev["DEL#1#DEL_x"], [("DEL", 300, 700, 0)])
            self.assertEqual(ev["INS#1#INS_x"], [("INS", 1000, 1000, 350)])
            with open(out) as fh:
                body = fh.read()
            self.assertNotIn("\ts1\t", body)
            self.assertIn("L\t3\t+\t5\t+\t0M", body)


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


class TestEndToEnd(unittest.TestCase):
    def test_bubbles_kmers_genotype(self):
        from panbreak.cli import main
        with tempfile.TemporaryDirectory() as d:
            gfa, hap = random_locus(d)
            bams = os.path.join(d, "bams")
            os.makedirs(bams)
            for sample, haps in {"s_ref": ["REF", "REF"], "s_del": ["DEL", "DEL"], "s_het": ["REF", "DEL"], "s_ins": ["INS", "REF"]}.items():
                write_reads(os.path.join(bams, sample + ".bam"), sample, [hap[h] for h in haps])
            main(["bubbles", "-g", gfa, "-r", "REF", "-o", os.path.join(d, "b")])
            main(["kmers", "-b", os.path.join(d, "b"), "-t", "1"])
            main(["genotype", "-k", os.path.join(d, "b"), "-B", bams, "-t", "1", "-o", os.path.join(d, "g")])
            with open(os.path.join(d, "b", "bubbles.tsv")) as fh:
                types = {r["bubble"]: r["type"] for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(sorted(types.values()), ["DEL", "INS"])
            with open(os.path.join(d, "g", "genotypes.tsv")) as fh:
                gt = {(r["sample"], types[r["bubble"]]): (r["GT"], r["FT"]) for r in csv.DictReader(fh, delimiter="\t")}
            expect = {("s_ref", "DEL"): "0/0", ("s_del", "DEL"): "1/1", ("s_het", "DEL"): "0/1", ("s_ins", "DEL"): "0/0",
                      ("s_ref", "INS"): "0/0", ("s_del", "INS"): "0/0", ("s_het", "INS"): "0/0", ("s_ins", "INS"): "0/1"}
            self.assertEqual({k: v[0] for k, v in gt.items()}, expect)
            self.assertTrue(all(v[1] == "PASS" for v in gt.values()))
            self.assertTrue(os.path.exists(os.path.join(d, "g", "genotypes.vcf.gz")))

    @unittest.skipUnless(shutil.which("minigraph") and shutil.which("odgi"), "minigraph and odgi not installed")
    def test_minigraph_builder(self):
        import pysam
        from panbreak.cli import main
        rng = random.Random(11)
        seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        a, x, b, ins, c = seq(5000), seq(1200), seq(4000), seq(900), seq(5000)
        genomes = {"R": a + x + b + c, "D": a + b + c, "I": a + x + b + ins + c, "J": a + x + b + ins + c}
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "locus.fa")
            with open(fa, "w") as fh:
                for g, s in genomes.items():
                    fh.write(f">{g}#1#{g}_x\n{s}\n")
            pysam.tabix_compress(fa, fa + ".gz")
            pysam.faidx(fa + ".gz")
            main(["graph", "-i", fa + ".gz", "-r", "R", "-b", "minigraph", "-t", "1", "-o", os.path.join(d, "g")])
            main(["bubbles", "-g", os.path.join(d, "g", "locus.gfa"), "-r", "R", "-o", os.path.join(d, "b")])
            with open(os.path.join(d, "b", "bubbles.tsv")) as fh:
                rows = [(r["type"], r["carriers"]) for r in csv.DictReader(fh, delimiter="\t")]
            self.assertEqual(sorted(rows), [("DEL", "D"), ("INS", "I,J")])

    def test_check(self):
        from panbreak.cli import main
        buf = io.StringIO()
        with redirect_stdout(buf):
            main(["check"])
        self.assertIn("pysam", buf.getvalue())
        with self.assertRaises(SystemExit):
            with redirect_stdout(io.StringIO()):
                main(["check", "--require", "no-such-tool"])


if __name__ == "__main__":
    unittest.main()
