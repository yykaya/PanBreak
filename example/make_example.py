import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from panbreak.common import rc, write_fasta, write_tsv

BASES = "ACGT"


def rand_seq(rng, n, gc=0.36):
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices(BASES, weights=w, k=n))


def mutate(rng, seq, rate):
    s = list(seq)
    for i in rng.sample(range(len(s)), int(len(s) * rate)):
        s[i] = rng.choice([b for b in BASES if b != s[i]])
    return "".join(s)


def make_te_library(rng):
    ltr = rand_seq(rng, 300)
    copia = "TG" + ltr[2:-2] + "CA" + rand_seq(rng, 4100, 0.42) + "TG" + ltr[2:-2] + "CA"
    helitron = "TC" + rand_seq(rng, 2094, 0.33) + "CTAG"
    tir = rand_seq(rng, 14)
    hat = tir + rand_seq(rng, 972, 0.40) + rc(tir)
    gypsy_ltr = rand_seq(rng, 420)
    gypsy = gypsy_ltr + rand_seq(rng, 4200) + gypsy_ltr
    line = rand_seq(rng, 3100, 0.40) + "A" * 25
    return {"Copia_1#LTR/Copia": copia, "Helitron_1#DNA/Helitron": helitron, "hAT_1#DNA/hAT": hat,
            "Gypsy_2#LTR/Gypsy": gypsy, "L1_3#LINE/L1": line}


GENES = {
    "geneW": {"strand": "+", "exons": [(31001, 31600), (32201, 33800)], "cds": (31101, 33600)},
    "geneX": {"strand": "-", "exons": [(44001, 44450), (44701, 44800), (45201, 45260), (45501, 45560), (45801, 45900), (46201, 46300), (50201, 50500)],
              "cds": (44741, 50350)},
    "geneY": {"strand": "+", "exons": [(56001, 56800), (57201, 58500)], "cds": (56101, 58300)},
}


def variants(te):
    return {
        "TE_DEL": ("del", 36001, 37000),
        "UTR_DEL": ("del", 44011, 44455),
        "STR_EXP": ("ins_unit", 46720, 46601, 2),
        "DEL22": ("del", 47001, 49200),
        "DEL15": ("del", 47301, 48800),
        "LTR_INS": ("ins_te", 48500, "Copia_1#LTR/Copia", 5),
        "DUP": ("dup", 45601, 50100, 2),
        "HEL_INS": ("ins_te", 50650, "Helitron_1#DNA/Helitron", 0),
    }


HAPLOTYPES = {
    "ATG1": [], "ATG2": ["DEL22"], "ATG3": ["DEL22", "UTR_DEL", "TE_DEL"], "ATG4": ["DEL15"],
    "ATG5": ["LTR_INS"], "ATG6": ["HEL_INS", "TE_DEL"], "ATG7": ["UTR_DEL", "STR_EXP"], "ATG8": ["DUP"],
}
REVERSED = {"ATG4", "ATG7"}
POPULATIONS = {
    "North": {"ATG1": 0.45, "ATG5": 0.15, "ATG6": 0.20, "ATG7": 0.20},
    "Central": {"ATG2": 0.45, "ATG3": 0.20, "ATG1": 0.20, "ATG4": 0.15},
    "South": {"ATG7": 0.35, "ATG8": 0.30, "ATG1": 0.20, "ATG3": 0.15},
}


def build_genome(rng, anc, te_copies, names):
    seq = mutate(rng, anc, 0.0015)
    vs = variants(te_copies)
    edits = []
    for n in names:
        v = vs[n]
        if v[0] == "del":
            edits.append((v[1], "del", v[2]))
        elif v[0] == "ins_te":
            edits.append((v[1], "ins", te_copies[v[2]], v[3]))
        elif v[0] == "ins_unit":
            edits.append((v[1], "unit", v[2], v[3]))
        elif v[0] == "dup":
            edits.append((v[2], "dup", v[1], v[3]))
    shifts = []
    for e in sorted(edits, key=lambda e: -e[0]):
        p = e[0]
        if e[1] == "del":
            seq = seq[:p - 1] + seq[e[2]:]
            shifts.append((p, -(e[2] - p + 1), e[2]))
        elif e[1] == "ins":
            t = seq[p - e[3]:p] if e[3] else ""
            ins = e[2] + t
            seq = seq[:p] + ins + seq[p:]
            shifts.append((p + 1, len(ins), p))
        elif e[1] == "unit":
            unit = seq[e[2] - 1:p]
            seq = seq[:p] + unit * e[3] + seq[p:]
            shifts.append((p + 1, len(unit) * e[3], p))
        elif e[1] == "dup":
            unit = seq[e[2] - 1:p]
            seq = seq[:p] + unit * e[3] + seq[p:]
            shifts.append((p + 1, len(unit) * e[3], p))

    def lift(x):
        y = x
        for p, d, end in shifts:
            if d < 0 and p <= x <= end:
                return None
            if x >= p:
                y += d if d > 0 else d
        return y
    return seq, lift


def gff_lines(contig, lift, length, reverse, offset):
    out = []
    for g, m in GENES.items():
        def conv(a, b):
            la, lb = lift(a), lift(b)
            if la is None or lb is None:
                return None
            la, lb = la + offset, lb + offset
            return (length - lb + 1, length - la + 1) if reverse else (la, lb)
        st = m["strand"] if not reverse else {"+": "-", "-": "+"}[m["strand"]]
        gs = conv(m["exons"][0][0], m["exons"][-1][1])
        if gs is None:
            continue
        out.append(f"{contig}\tsim\tgene\t{gs[0]}\t{gs[1]}\t.\t{st}\t.\tID={g};Name={g}")
        out.append(f"{contig}\tsim\tmRNA\t{gs[0]}\t{gs[1]}\t.\t{st}\t.\tID={g}.1;Parent={g}")
        for a, b in m["exons"]:
            ex = conv(a, b)
            if ex:
                out.append(f"{contig}\tsim\texon\t{ex[0]}\t{ex[1]}\t.\t{st}\t.\tParent={g}.1")
        cs, ce = m["cds"]
        for a, b in m["exons"]:
            a2, b2 = max(a, cs), min(b, ce)
            if a2 <= b2:
                c = conv(a2, b2)
                if c:
                    out.append(f"{contig}\tsim\tCDS\t{c[0]}\t{c[1]}\t.\t{st}\t0\tParent={g}.1")
    return sorted(out, key=lambda l: int(l.split("\t")[3]))


def simulate_reads(rng, haps, cov, r1, r2, read_len=150, frag_mu=380, frag_sd=60, err=0.002):
    comp = str.maketrans("ACGT", "TGCA")
    with open(r1, "w") as f1, open(r2, "w") as f2:
        i = 0
        for h in haps:
            n = int(cov * len(h) / (2 * read_len) / len(haps))
            for _ in range(n):
                fl = max(2 * read_len, min(800, int(rng.gauss(frag_mu, frag_sd))))
                st = rng.randrange(0, len(h) - fl)
                frag = h[st:st + fl]
                if rng.random() < 0.5:
                    frag = frag.translate(comp)[::-1]
                a, b = frag[:read_len], frag.translate(comp)[::-1][:read_len]
                a, b = [mutate(rng, x, err) if rng.random() < 0.5 else x for x in (a, b)]
                q = "F" * read_len
                f1.write(f"@r{i}/1\n{a}\n+\n{q}\n")
                f2.write(f"@r{i}/2\n{b}\n+\n{q}\n")
                i += 1


def main():
    ap = argparse.ArgumentParser(description="simulate the ATG1..ATG8 example data set")
    ap.add_argument("-o", "--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))
    ap.add_argument("--per-pop", type=int, default=40)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("-t", "--threads", type=int, default=8)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = a.out
    for sub in ("assemblies", "bams"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    lib = make_te_library(rng)
    with open(os.path.join(out, "te_library.fa"), "w") as fh:
        for k, v in lib.items():
            write_fasta(fh, k, v)
    copies = {k: mutate(rng, v, 0.02) for k, v in lib.items()}
    anc = rand_seq(rng, 100000)
    hat = mutate(rng, lib["hAT_1#DNA/hAT"], 0.03)
    anc = anc[:36000] + hat + anc[36000 + len(hat):]
    genomes, rows = {}, []
    for g, names in HAPLOTYPES.items():
        seq, lift = build_genome(rng, anc, copies, names)
        offset = 0
        contig = "chr1"
        if g == "ATG8":
            pre = rand_seq(rng, 7000)
            seq, offset, contig = pre + seq, 7000, "ctg_8"
        rev = g in REVERSED
        if rev:
            seq = rc(seq)
        genomes[g] = seq if not rev else rc(seq)
        fa = os.path.join(out, "assemblies", f"{g}.fa")
        with open(fa, "w") as fh:
            write_fasta(fh, contig, seq)
        with open(os.path.join(out, "assemblies", f"{g}.gff3"), "w") as fh:
            fh.write("##gff-version 3\n")
            fh.write("\n".join(gff_lines(contig, lift, len(seq), rev, offset)) + "\n")
        rows.append({"genome": g, "fasta": os.path.join("data", "assemblies", f"{g}.fa"), "gff": os.path.join("data", "assemblies", f"{g}.gff3"),
                     "variants": ",".join(names) or "none", "orientation": "-" if rev else "+"})
    write_tsv(os.path.join(out, "assemblies.tsv"), rows, ["genome", "fasta", "gff"])
    write_tsv(os.path.join(out, "assembly_truth.tsv"), rows, ["genome", "variants", "orientation"])
    ref = os.path.join(out, "assemblies", "ATG1.fa")
    minimap2, samtools = shutil.which("minimap2"), shutil.which("samtools")
    subprocess.run([samtools, "faidx", ref], check=True)
    samples = []
    n = 0
    for pop, freq in POPULATIONS.items():
        for i in range(a.per_pop):
            n += 1
            s = f"S{n:03d}"
            h1 = rng.choices(list(freq), weights=list(freq.values()))[0]
            h2 = rng.choices(list(freq), weights=list(freq.values()))[0] if rng.random() < 0.08 else h1
            cov = round(rng.uniform(8, 25), 1) if rng.random() > 0.04 else 1.5
            samples.append({"sample": s, "population": pop, "haplotype1": h1, "haplotype2": h2, "coverage": cov})
    with tempfile.TemporaryDirectory() as tmp:
        for k, sm in enumerate(samples, 1):
            s = sm["sample"]
            haps = [mutate(rng, genomes[sm["haplotype1"]], 0.0005), mutate(rng, genomes[sm["haplotype2"]], 0.0005)]
            r1, r2 = os.path.join(tmp, "r1.fq"), os.path.join(tmp, "r2.fq")
            simulate_reads(rng, haps, sm["coverage"], r1, r2)
            full = os.path.join(tmp, f"{s}.full.bam")
            mm = subprocess.Popen([minimap2, "-ax", "sr", "-t", str(a.threads), "-R", f"@RG\\tID:{s}\\tSM:{s}", ref, r1, r2],
                                  stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            subprocess.run([samtools, "sort", "-@", "2", "-o", full, "-"], stdin=mm.stdout, check=True, stderr=subprocess.DEVNULL)
            mm.wait()
            subprocess.run([samtools, "index", full], check=True)
            bam = os.path.join(out, "bams", f"{s}.bam")
            subprocess.run([samtools, "view", "-b", "-o", bam, full, "chr1:20000-75000"], check=True)
            subprocess.run([samtools, "index", bam], check=True)
            if k % 20 == 0:
                print(f"{k}/{len(samples)} samples", file=sys.stderr, flush=True)
    write_tsv(os.path.join(out, "samples.tsv"), samples, ["sample", "population", "haplotype1", "haplotype2", "coverage"])
    print(f"example data written to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
