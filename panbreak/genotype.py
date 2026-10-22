import glob
import math
import os
from collections import defaultdict
from multiprocessing import Pool

from .common import die, log, mkdir, parse_region, rc, read_fasta, read_tsv, write_tsv

STATE = {}


def bam_list(items):
    out = []
    for it in items:
        if os.path.isdir(it):
            out += sorted(glob.glob(os.path.join(it, "*.bam")) + glob.glob(os.path.join(it, "*.cram")))
        elif it.endswith((".bam", ".cram")):
            out.append(it)
        else:
            out += [l.strip().split("\t")[-1] for l in open(it) if l.strip() and not l.startswith("#")]
    if not out:
        die("no BAM/CRAM files given")
    return out


def sample_name(bam, path, mode):
    if mode == "file":
        return os.path.basename(path).split(".")[0]
    rg = bam.header.to_dict().get("RG", [])
    sm = {r.get("SM") for r in rg if r.get("SM")}
    return sm.pop() if len(sm) == 1 else os.path.basename(path).split(".")[0]


def init(kmap, k, windows, opts):
    STATE.update({"kmap": kmap, "k": k, "windows": windows, "opts": opts})


def scan(path):
    import pysam
    kmap, k, windows, o = STATE["kmap"], STATE["k"], STATE["windows"], STATE["opts"]
    try:
        bam = pysam.AlignmentFile(path, reference_filename=o["cram_ref"]) if path.endswith(".cram") else pysam.AlignmentFile(path)
    except Exception as e:
        return {"path": path, "error": str(e)}
    sample = sample_name(bam, path, o["sample_from"])
    if windows:
        its = []
        for c, s, e in windows:
            try:
                its.append(bam.fetch(c, max(0, s), e))
            except ValueError:
                pass
        reads = (r for it in its for r in it)
    else:
        reads = bam.fetch(until_eof=True) if bam.has_index() else bam.fetch(until_eof=True)
    seen = set()
    counts = defaultdict(lambda: [0, 0])
    jx = defaultdict(lambda: defaultdict(int))
    n = 0
    for r in reads:
        if r.is_secondary or r.is_supplementary or r.is_qcfail or (r.is_duplicate and not o["keep_dups"]):
            continue
        if r.mapping_quality < o["min_mapq"] and not r.is_unmapped:
            continue
        key = (r.query_name, r.is_read1)
        if key in seen:
            continue
        seen.add(key)
        s = r.query_sequence
        if not s:
            continue
        n += 1
        hits = set()
        for i in range(len(s) - k + 1):
            h = kmap.get(s[i:i + k])
            if h:
                hits.update(h)
        if not hits:
            continue
        alt_b = {b for b, a, j in hits if a == 1}
        for b, a, j in hits:
            if a == 1:
                jx[b][f"a{j}"] += 1
            elif b not in alt_b:
                jx[b][f"r{j}"] += 1
        for b in alt_b:
            counts[b][1] += 1
        for b in {b for b, a, j in hits if a == 0} - alt_b:
            counts[b][0] += 1
    bam.close()
    return {"path": path, "sample": sample, "reads": n, "counts": {b: tuple(v) for b, v in counts.items()},
            "jx": {b: dict(v) for b, v in jx.items()}}


def call_gt(a, r, model, h, o):
    p1 = o["ploidy"] == 1
    miss = "." if p1 else "./."
    if model == "presence":
        cn = round(1 + 2 * a / r, 2) if r else None
        if a >= o["min_alt"]:
            return {"GT": "1" if p1 else "./1", "GQ": 99, "PL": ".", "FT": "PASS", "call": "carrier", "dosage": None, "carrier": 1, "copy_estimate": cn}
        if a == 0 and r >= o["min_ref"]:
            gq = min(99, int(round(10 * (r / 2) / math.log(10))))
            return {"GT": "0" if p1 else "0/0", "GQ": gq, "PL": ".", "FT": "PASS" if gq >= o["min_gq"] else "LowGQ",
                    "call": "hom_ref" if gq >= o["min_gq"] else "nodata", "dosage": 0 if gq >= o["min_gq"] else None,
                    "carrier": 0 if gq >= o["min_gq"] else None, "copy_estimate": cn}
        if a > 0:
            return {"GT": miss, "GQ": 0, "PL": ".", "FT": "LowAltSupport", "call": "weak", "dosage": None, "carrier": None, "copy_estimate": cn}
        return {"GT": miss, "GQ": 0, "PL": ".", "FT": "NoData", "call": "nodata", "dosage": None, "carrier": None, "copy_estimate": None}
    e = o["error"]
    if p1:
        gts, ps, dose = ["0", "1"], [e, 1 - e], [0, 1]
    else:
        h = min(max(h, 0.05), 0.95)
        gts, ps, dose = ["0/0", "0/1", "1/1"], [e, h, 1 - e], [0, 1, 2]
    if a + r < o["min_depth"]:
        return {"GT": miss, "GQ": 0, "PL": ".", "FT": "LowAltSupport" if a else "NoData", "call": "weak" if a else "nodata",
                "dosage": None, "carrier": None, "copy_estimate": None}
    ll = [a * math.log(p) + r * math.log(1 - p) for p in ps]
    m = max(ll)
    pl = [int(round(-10 * (x - m) / math.log(10))) for x in ll]
    i = pl.index(0)
    gq = min(99, sorted(pl)[1])
    out = {"GT": gts[i], "GQ": gq, "PL": ",".join(map(str, pl)), "dosage": dose[i], "carrier": int(dose[i] > 0), "copy_estimate": None}
    if dose[i] > 0 and a < o["min_alt"]:
        out.update({"FT": "LowAltSupport", "call": "weak", "dosage": None, "carrier": None})
    elif dose[i] == 0 and r < o["min_ref"]:
        out.update({"FT": "LowRefSupport", "call": "nodata", "dosage": None, "carrier": None})
    elif gq < o["min_gq"]:
        out.update({"FT": "LowGQ", "call": "lowqual", "dosage": None, "carrier": None})
    else:
        out.update({"FT": "PASS", "call": {0: "hom_ref", 1: "het" if not p1 else "hom_alt", 2: "hom_alt"}[dose[i]]})
    return out


def resolve_windows(junctions, bubbles, fasta, region, pad):
    import pysam
    fa = pysam.FastaFile(fasta)
    contigs = [region[0]] if region else list(fa.references)
    seqs = {}
    for c in contigs:
        s, e = (region[1] - 1, region[2]) if region else (0, fa.get_reference_length(c))
        seqs[c] = (s, fa.fetch(c, s, e).upper())
    wins, missing = [], []
    for b in bubbles:
        found = []
        for j in junctions.get(b, []):
            w = j["window"]
            for c, (s0, sq) in seqs.items():
                for probe, shift in ((w, int(j["breakpoint"])), (rc(w), len(w) - int(j["breakpoint"]))):
                    i = sq.find(probe)
                    if i != -1 and sq.find(probe, i + 1) == -1:
                        found.append((c, s0 + i + shift))
        if not found:
            missing.append(b)
            continue
        c = found[0][0]
        ps = [p for cc, p in found if cc == c]
        wins.append((b, c, min(ps) - pad, max(ps) + pad))
    return wins, missing


def genotype(args):
    out = mkdir(args.out)
    d = args.kmers
    km = read_tsv(os.path.join(d, "junction_kmers.tsv"))
    summ = {r["bubble"]: r for r in read_tsv(os.path.join(d, "junction_summary.tsv"))}
    k = len(km[0]["kmer"]) if km else 31
    bubbles = [b for b, s in summ.items() if s["genotypable"] == "yes"]
    bi = {b: i for i, b in enumerate(bubbles)}
    kmap = defaultdict(set)
    for r in km:
        if r["bubble"] not in bi:
            continue
        v = (bi[r["bubble"]], 1 if r["allele"] == "alt" else 0, int(r["junction"]))
        kmap[r["kmer"]].add(v)
        kmap[rc(r["kmer"])].add(v)
    kmap = {x: tuple(v) for x, v in kmap.items()}
    windows = None
    if args.bam_ref:
        jdir = args.bubbles or d
        junc = defaultdict(list)
        for j in read_tsv(os.path.join(jdir, "bubble_junctions.tsv")):
            if j["allele"] == "ref":
                junc[j["bubble"]].append(j)
        region = parse_region(args.bam_region) if args.bam_region else None
        wins, missing = resolve_windows(junc, bubbles, args.bam_ref, region, args.window)
        if missing:
            log(f"no unique reference anchor for {', '.join(missing)}: all reads are scanned")
            windows = None
        else:
            windows = [(c, s, e) for _, c, s, e in wins]
            write_tsv(os.path.join(out, "read_windows.tsv"), [{"bubble": b, "contig": c, "start": s, "end": e} for b, c, s, e in wins])
            log(f"reads fetched from {len(windows)} windows on {args.bam_ref}")
    paths = bam_list(args.bams)
    log(f"{len(paths)} alignment files, {len(bubbles)} genotypable bubbles, {len(kmap) // 2} junction k-mers, {args.threads} workers")
    opts = {"min_mapq": args.min_mapq, "keep_dups": args.keep_duplicates, "cram_ref": args.cram_ref, "sample_from": args.sample_from}
    rows, qc = [], []
    gopt = {"ploidy": args.ploidy, "error": args.error, "min_alt": args.min_alt, "min_ref": args.min_ref,
            "min_depth": args.min_depth, "min_gq": args.min_gq}
    with Pool(args.threads, initializer=init, initargs=(kmap, k, windows, opts)) as pool:
        for i, res in enumerate(pool.imap_unordered(scan, paths, chunksize=4), 1):
            if "error" in res:
                log(f"skipped {res['path']}: {res['error']}")
                continue
            called = 0
            for b in bubbles:
                ref_n, alt_n = res["counts"].get(bi[b], (0, 0))
                s = summ[b]
                h = float(s["het_alt_fraction"]) if s["het_alt_fraction"] not in ("NA", "") else 0.5
                g = call_gt(alt_n, ref_n, s["model"], h, gopt)
                jd = res["jx"].get(bi[b], {})
                called += g["FT"] == "PASS"
                rows.append({"sample": res["sample"], "bubble": b, "model": s["model"], "alt_reads": alt_n, "ref_reads": ref_n,
                             "alt_fraction": round(alt_n / (alt_n + ref_n), 3) if alt_n + ref_n else "NA",
                             "alt_junction_reads": ",".join(str(jd.get(f"a{j}", 0)) for j in range(int(s["alt_junctions"]) or 1)),
                             "ref_junction_reads": ",".join(str(jd.get(f"r{j}", 0)) for j in range(max(1, int(s["ref_junctions"])))),
                             **g})
            qc.append({"sample": res["sample"], "file": res["path"], "reads_scanned": res["reads"], "bubbles": len(bubbles),
                       "called": called, "call_rate": round(called / len(bubbles), 3) if bubbles else "NA",
                       "median_junction_reads": sorted(sum(res["counts"].get(bi[b], (0, 0))) for b in bubbles)[len(bubbles) // 2] if bubbles else 0})
            if i % 100 == 0:
                log(f"{i}/{len(paths)} files")
    write_outputs(rows, qc, out, args.bubbles or d, args.ploidy)


GT_FIELDS = ["sample", "bubble", "model", "alt_reads", "ref_reads", "alt_fraction", "alt_junction_reads", "ref_junction_reads",
             "GT", "GQ", "PL", "FT", "call", "dosage", "carrier", "copy_estimate"]


def write_outputs(rows, qc, out, bdir, ploidy):
    rows.sort(key=lambda r: (r["sample"], r["bubble"]))
    write_tsv(os.path.join(out, "genotypes.tsv"), rows, GT_FIELDS)
    if qc:
        write_tsv(os.path.join(out, "sample_qc.tsv"), sorted(qc, key=lambda r: r["sample"]))
    samples = sorted({r["sample"] for r in rows})
    bubbles = sorted({r["bubble"] for r in rows})
    by = {(r["sample"], r["bubble"]): r for r in rows}
    for name, key in (("dosage_matrix.tsv", "dosage"), ("carrier_matrix.tsv", "carrier")):
        mat = []
        for s in samples:
            m = {"sample": s}
            for b in bubbles:
                v = by.get((s, b), {}).get(key)
                m[b] = "NA" if v in (None, "", "NA") else v
            mat.append(m)
        write_tsv(os.path.join(out, name), mat, ["sample"] + bubbles)
    write_vcf(os.path.join(out, "genotypes.vcf"), rows, samples, bubbles, bdir)
    summary = defaultdict(lambda: defaultdict(int))
    for r in rows:
        summary[r["bubble"]][r["call"]] += 1
    for b in bubbles:
        log(f"{b}: " + "  ".join(f"{k} {v}" for k, v in sorted(summary[b].items())))
    log(f"genotypes for {len(samples)} samples written to {out}")


def write_vcf(path, rows, samples, bubbles, bdir):
    import pysam
    btab = os.path.join(bdir, "bubbles.annotated.tsv")
    binfo = {r["bubble"]: r for r in read_tsv(btab if os.path.exists(btab) else os.path.join(bdir, "bubbles.tsv"))}
    (rname, rseq), = list(read_fasta(os.path.join(bdir, "reference.fa")))
    contig, r0, r1 = parse_region(rname)
    by = {(r["sample"], r["bubble"]): r for r in rows}
    alt_sym = {"DEL": "<DEL>", "INS": "<INS>", "DUP": "<DUP:TANDEM>", "INV": "<INV>", "REPL": "<CPX>"}
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.3\n##source=panbreak\n")
        fh.write(f"##contig=<ID={contig},length={r1}>\n")
        for a, desc in (("DEL", "Deletion"), ("INS", "Insertion"), ("DUP:TANDEM", "Tandem duplication"), ("INV", "Inversion"),
                        ("CPX", "Replacement of reference sequence by non-reference sequence")):
            fh.write(f'##ALT=<ID={a},Description="{desc}">\n')
        for i, n, t, desc in (("SVTYPE", 1, "String", "Bubble type"), ("END", 1, "Integer", "End position"),
                              ("SVLEN", 1, "Integer", "Reference length (DEL, DUP, INV, CPX) or inserted length (INS)"),
                              ("ALTLEN", 1, "Integer", "Length of the non-reference allele in the representative assembly"),
                              ("SVCLASS", 1, "String", "Annotated class"), ("TE", 1, "String", "Best transposable element hit"),
                              ("CONTEXT", 1, "String", "Gene context"), ("ASM_CARRIERS", 1, "Integer", "Assemblies carrying the bubble"),
                              ("MODEL", 1, "String", "Genotype model (biallelic or presence)")):
            fh.write(f'##INFO=<ID={i},Number={n},Type={t},Description="{desc}">\n')
        for i, n, t, desc in (("GT", 1, "String", "Genotype"), ("GQ", 1, "Integer", "Genotype quality"),
                              ("AD", "R", "Integer", "Reads with reference and alternative junction k-mers"),
                              ("DP", 1, "Integer", "Informative reads"), ("PL", "G", "Integer", "Phred-scaled genotype likelihoods"),
                              ("FT", 1, "String", "Sample filter"), ("CN", 1, "Float", "Copy number estimate (presence model)")):
            fh.write(f'##FORMAT=<ID={i},Number={n},Type={t},Description="{desc}">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for b in sorted(bubbles, key=lambda x: int(binfo[x]["pos"])):
            bi = binfo[b]
            pos = max(1, int(bi["pos"]))
            refb = rseq[pos - r0] if 0 <= pos - r0 < len(rseq) else "N"
            svlen = int(bi["alt_len"]) if bi["type"] == "INS" else int(bi["ref_len"])
            model = next((by[(s, b)]["model"] for s in samples if (s, b) in by), "biallelic")
            info = f"SVTYPE={bi['type']};END={bi['end']};SVLEN={svlen};ALTLEN={bi['alt_len']};ASM_CARRIERS={bi['n_carriers']};MODEL={model}"
            for key, col in (("SVCLASS", "sv_class"), ("TE", "te_family"), ("CONTEXT", "context")):
                v = bi.get(col)
                if v and v not in (".", "NA"):
                    info += f";{key}=" + v.replace(" ", "_").replace(";", ",").replace("=", ":")
            cols = []
            for s in samples:
                r = by.get((s, b))
                if r is None:
                    cols.append("./.")
                    continue
                cn = r["copy_estimate"] if r["copy_estimate"] not in (None, "NA") else "."
                cols.append(f"{r['GT']}:{r['GQ']}:{r['ref_reads']},{r['alt_reads']}:{int(r['ref_reads']) + int(r['alt_reads'])}:{r['PL']}:{r['FT']}:{cn}")
            fh.write(f"{contig}\t{pos}\t{b}\t{refb}\t{alt_sym.get(bi['type'], '<SV>')}\t.\tPASS\t{info}\tGT:GQ:AD:DP:PL:FT:CN\t" + "\t".join(cols) + "\n")
    pysam.tabix_compress(path, path + ".gz", force=True)
    pysam.tabix_index(path + ".gz", preset="vcf", force=True)
    os.remove(path)


def merge(args):
    rows, qc = [], []
    for d in args.inputs:
        rows += read_tsv(os.path.join(d, "genotypes.tsv"))
        q = os.path.join(d, "sample_qc.tsv")
        if os.path.exists(q):
            qc += read_tsv(q)
    for r in rows:
        for key in ("dosage", "carrier"):
            r[key] = None if r[key] in ("NA", "") else int(r[key])
        r["copy_estimate"] = None if r["copy_estimate"] in ("NA", "") else r["copy_estimate"]
    write_outputs(rows, qc, mkdir(args.out), args.bubbles, args.ploidy)
