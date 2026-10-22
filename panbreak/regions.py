import os
import re
import tempfile

from .common import bgzip_index, die, log, mkdir, opener, rc, read_fasta, read_tsv, run, tool, write_fasta, write_tsv

GENE_TYPES = {"gene", "mRNA", "transcript", "ncRNA_gene", "pseudogene", "lnc_RNA"}


def gff_attrs(text):
    out = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def gene_names(a):
    names = set()
    for key in ("ID", "Name", "gene_id", "gene", "locus_tag", "Parent"):
        for v in a.get(key, "").split(","):
            if v:
                v = re.sub(r"^(gene|transcript|mRNA):", "", v)
                names.update({v, re.sub(r"_\d+$", "", v), v.split(".")[0]})
    return names


def find_gene(gff, gene):
    hits = []
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in GENE_TYPES:
                continue
            a = gff_attrs(f[8])
            if gene in gene_names(a):
                rank = (f[2] != "gene", a.get("extra_copy_number", "0") != "0", int(f[3]))
                hits.append((rank, f[0], int(f[3]), int(f[4]), f[6]))
    if not hits:
        return None
    hits.sort()
    return hits[0][1:]


class Genome:
    def __init__(self, path):
        self.path = path
        self.fa = None
        self.cache = {}
        try:
            import pysam
            self.fa = pysam.FastaFile(path)
        except Exception:
            self.fa = None

    def length(self, contig):
        if self.fa is not None:
            return self.fa.get_reference_length(contig)
        return len(self.seq(contig))

    def seq(self, contig):
        if contig not in self.cache:
            for name, s in read_fasta(self.path):
                if name == contig:
                    self.cache[contig] = s
                    break
            else:
                die(f"contig {contig} not in {self.path}")
        return self.cache[contig]

    def fetch(self, contig, start, end):
        if self.fa is not None:
            return self.fa.fetch(contig, start - 1, end).upper()
        return self.seq(contig)[start - 1:end]


def locate_by_alignment(genome_fa, query_fa, threads, preset):
    minimap2 = tool("minimap2")
    with tempfile.NamedTemporaryFile("w+", suffix=".paf") as paf:
        run([minimap2, "-x", preset, "-c", "--secondary=no", "-t", threads, genome_fa, query_fa], stdout=paf)
        paf.seek(0)
        groups = {}
        qlen = None
        for line in paf:
            f = line.split("\t")
            qlen = int(f[1])
            key = (f[5], f[4])
            g = groups.setdefault(key, {"match": 0, "ts": [], "te": [], "qs": [], "qe": []})
            g["match"] += int(f[9])
            g["ts"].append(int(f[7]))
            g["te"].append(int(f[8]))
            g["qs"].append(int(f[2]))
            g["qe"].append(int(f[3]))
    if not groups:
        return None
    (contig, strand), g = max(groups.items(), key=lambda kv: kv[1]["match"])
    if g["match"] < 0.3 * qlen:
        return None
    ts, te = min(g["ts"]), max(g["te"])
    left, right = min(g["qs"]), qlen - max(g["qe"])
    if strand == "+":
        return contig, max(1, ts + 1 - left), te + right, "+"
    return contig, max(1, ts + 1 - right), te + left, "-"


def extract(args):
    out = mkdir(args.out)
    mkdir(os.path.join(out, "regions"))
    table = read_tsv(args.assemblies)
    if not table or "genome" not in table[0] or "fasta" not in table[0]:
        die("assemblies table needs columns genome, fasta (and gff for --gene)")
    names = [r["genome"] for r in table]
    if args.ref not in names:
        die(f"reference {args.ref} not in assemblies table")
    table.sort(key=lambda r: r["genome"] != args.ref)
    query_fa = None
    ref_strand = "+"
    if args.gene:
        ref_row = table[0]
        if not ref_row.get("gff"):
            die("--gene needs a gff column for the reference genome")
        hit = find_gene(ref_row["gff"], args.gene)
        if not hit:
            die(f"{args.gene} not found in {ref_row['gff']}")
        ref_strand = hit[3]
        if any(not r.get("gff") for r in table):
            gref = Genome(ref_row["fasta"])
            query_fa = os.path.join(out, "regions", "reference_locus.fa")
            with open(query_fa, "w") as fh:
                write_fasta(fh, "locus", gref.fetch(hit[0], hit[1], hit[2]))
    elif args.ref_locus:
        from .common import parse_region
        contig, s, e = parse_region(args.ref_locus)
        gref = Genome(table[0]["fasta"])
        query_fa = os.path.join(out, "regions", "reference_locus.fa")
        with open(query_fa, "w") as fh:
            write_fasta(fh, "locus", gref.fetch(contig, s, e))
    else:
        die("give --gene (with gff column) or --ref-locus contig:start-end")
    rows = []
    combined = os.path.join(out, "locus.fa")
    with open(combined, "w") as all_fh:
        for r in table:
            g = r["genome"]
            hit, method = None, None
            if args.gene and r.get("gff"):
                hit = find_gene(r["gff"], args.gene)
                method = "annotation"
            if hit is None and query_fa:
                hit = locate_by_alignment(r["fasta"], query_fa, args.threads, args.preset)
                if hit is not None:
                    hit = hit[:3] + (ref_strand if hit[3] == "+" else {"+": "-", "-": "+"}[ref_strand],)
                method = "alignment"
            if hit is None:
                log(f"{g}: locus not found, skipped")
                continue
            contig, gs, ge, strand = hit
            genome = Genome(r["fasta"])
            clen = genome.length(contig)
            flip = strand != ref_strand
            up, down = (args.flank, args.flank) if args.up is None else (args.up, args.down)
            if strand == "+":
                lo, hi = gs - up, ge + down
            else:
                lo, hi = gs - down, ge + up
            lo, hi = max(1, lo), min(clen, hi)
            seq = genome.fetch(contig, lo, hi)
            if flip:
                seq = rc(seq)
                rgs, rge = hi - ge + 1, hi - gs + 1
            else:
                rgs, rge = gs - lo + 1, ge - lo + 1
            path = f"{g}#1#{g}_{args.locus}"
            with open(os.path.join(out, "regions", f"{g}.fa"), "w") as fh:
                write_fasta(fh, path, seq)
            write_fasta(all_fh, path, seq)
            rows.append({"genome": g, "path": path, "contig": contig, "start": lo, "end": hi,
                         "orientation": "-" if flip else "+", "gene_start": rgs, "gene_end": rge,
                         "gene_strand": ref_strand, "length": len(seq), "method": method})
            log(f"{g}: {contig}:{lo}-{hi} ({'reverse' if flip else 'forward'}, {method})")
    write_tsv(os.path.join(out, "regions.tsv"), rows,
              ["genome", "path", "contig", "start", "end", "orientation", "gene_start", "gene_end", "gene_strand", "length", "method"])
    bgzip_index(combined, combined + ".gz")
    os.remove(combined)
    log(f"{len(rows)} regions written to {out}/locus.fa.gz")
