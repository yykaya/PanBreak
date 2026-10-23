import os
import shutil
import subprocess
import tempfile
from collections import defaultdict

import numpy as np

from .common import canonical, log, opener, parse_region, rc, read_fasta, read_tsv, write_fasta, write_tsv
from .regions import gff_attrs


def size_str(n):
    n = int(n)
    return f"{n / 1000:.1f} kb" if n >= 1000 else f"{n} bp"


def identity(x, y):
    n = min(len(x), len(y))
    if n == 0:
        return 0.0
    a = np.frombuffer(x[:n].encode(), dtype=np.uint8)
    b = np.frombuffer(y[:n].encode(), dtype=np.uint8)
    return float((a == b).mean())


def load_models(gff, region_row, contig, offset, length):
    if not gff:
        return []
    lo, hi = offset + 1, offset + length
    feats = []
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            s, e, st = int(f[3]), int(f[4]), f[6]
            if f[0] == contig:
                pass
            elif region_row and f[0] == region_row["contig"]:
                r0, r1 = int(region_row["start"]), int(region_row["end"])
                if e < r0 or s > r1:
                    continue
                if region_row["orientation"] == "+":
                    s, e = s - r0 + 1, e - r0 + 1
                else:
                    s, e, st = r1 - e + 1, r1 - s + 1, {"+": "-", "-": "+"}.get(st, st)
            else:
                continue
            if e < lo or s > hi:
                continue
            feats.append((f[2], s, e, st, gff_attrs(f[8])))
    genes, tx_gene, parts = {}, {}, defaultdict(list)
    for t, s, e, st, a in feats:
        if t in ("gene", "ncRNA_gene", "pseudogene", "transposable_element_gene"):
            gid = a.get("ID") or a.get("gene_id") or a.get("Name")
            genes[gid] = {"id": a.get("Name") or a.get("gene_id") or gid, "start": s, "end": e, "strand": st, "tx": None}
    for t, s, e, st, a in feats:
        if t in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "tRNA", "miRNA"):
            g = (a.get("Parent") or a.get("gene_id") or "").split(",")[0]
            tids = [x for x in (a.get("ID"), a.get("transcript_id")) if x]
            if g in genes and genes[g]["tx"] is None and tids:
                genes[g]["tx"] = tids[0]
                for x in tids:
                    tx_gene[x] = g
    for t, s, e, st, a in feats:
        if t in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
            keys = [a["transcript_id"]] if a.get("transcript_id") else a.get("Parent", "").split(",")
            for p in keys:
                if p in tx_gene:
                    parts[tx_gene[p]].append((t, s, e))
    models = []
    for gid, g in genes.items():
        exons = sorted({(s, e) for t, s, e in parts[gid] if t == "exon"})
        cds = sorted({(s, e) for t, s, e in parts[gid] if t == "CDS"})
        if not exons:
            exons = cds or [(g["start"], g["end"])]
        rev = g["strand"] == "-"
        order = sorted(exons, key=lambda x: -x[0] if rev else x[0])
        introns = []
        srt = sorted(exons)
        for (s1, e1), (s2, e2) in zip(srt, srt[1:]):
            if s2 > e1 + 1:
                introns.append((e1 + 1, s2 - 1))
        introns = sorted(introns, key=lambda x: -x[0] if rev else x[0])
        utr5, utr3 = [], []
        explicit = [(t, s, e) for t, s, e in parts[gid] if t.endswith("UTR")]
        if explicit:
            utr5 = [(s, e) for t, s, e in explicit if t.startswith("five")]
            utr3 = [(s, e) for t, s, e in explicit if t.startswith("three")]
        elif cds:
            c0, c1 = cds[0][0], cds[-1][1]
            left = [(s, min(e, c0 - 1)) for s, e in exons if s < c0]
            right = [(max(s, c1 + 1), e) for s, e in exons if e > c1]
            utr5, utr3 = (right, left) if rev else (left, right)
        models.append({"id": g["id"], "start": g["start"], "end": g["end"], "strand": g["strand"],
                       "exons": order, "introns": introns, "cds": cds, "utr5": utr5, "utr3": utr3,
                       "tss": g["end"] if rev else g["start"], "tes": g["start"] if rev else g["end"]})
    return sorted(models, key=lambda m: m["start"])


def context(models, s, e, promoter, downstream):
    hits, cats = [], []
    rank = {"coding": 0, "UTR": 1, "exonic": 2, "intronic": 3, "promoter": 4, "downstream": 5}
    best_tss = None
    for m in models:
        if m["strand"] == "+":
            tss_d = s - m["tss"] if s >= m["tss"] else e - m["tss"] if e < m["tss"] else 0
        else:
            tss_d = m["tss"] - e if e <= m["tss"] else m["tss"] - s if s > m["tss"] else 0
        if best_tss is None or abs(tss_d) < abs(best_tss[1]):
            best_tss = (m["id"], tss_d)
        lab = []
        ov = lambda iv: [i for i, (a, b) in enumerate(iv, 1) if a <= e and s <= b]
        ex, it = ov(m["exons"]), ov(m["introns"])
        c = None
        if any(a <= e and s <= b for a, b in m["cds"]):
            c = "coding"
            lab.append("CDS")
            if any(a <= e and s <= b for a, b in m["utr5"] + m["utr3"]):
                lab.append("5' UTR" if any(a <= e and s <= b for a, b in m["utr5"]) else "3' UTR")
        elif any(a <= e and s <= b for a, b in m["utr5"] + m["utr3"]):
            c = "UTR"
            lab.append("5' UTR" if any(a <= e and s <= b for a, b in m["utr5"]) else "3' UTR")
        elif ex:
            c = "exonic"
        elif it:
            c = "intronic"
        if ex:
            lab.append(("exon " if len(ex) == 1 else "exons ") + (f"{ex[0]}" if len(ex) == 1 else f"{ex[0]}-{ex[-1]}"))
        if it:
            lab.append(("intron " if len(it) == 1 else "introns ") + (f"{it[0]}" if len(it) == 1 else f"{it[0]}-{it[-1]}"))
        if c is None:
            up = (m["tss"] - e) if m["strand"] == "+" else (s - m["tss"])
            dn = (s - m["tes"]) if m["strand"] == "+" else (m["tes"] - e)
            if 0 <= up <= promoter:
                c = "promoter"
                lab.append(f"{up} bp upstream of TSS")
            elif 0 <= dn <= downstream:
                c = "downstream"
                lab.append(f"{dn} bp downstream of TES")
        if c:
            cats.append(c)
            hits.append(f"{m['id']}: {', '.join(dict.fromkeys(lab))}" if lab else m["id"])
    cat = min(cats, key=lambda x: rank[x]) if cats else "intergenic"
    return cat, "; ".join(hits) if hits else "intergenic", best_tss


def te_name_class(name):
    if "#" in name:
        n, c = name.split("#", 1)
        return n, c
    return name, "Unknown"


def blast_te(seqs, te_lib, blastn, threads, tmp):
    q = os.path.join(tmp, "q.fa")
    with open(q, "w") as fh:
        for k, s in seqs.items():
            if len(s) >= 30:
                write_fasta(fh, k, s)
    mk = shutil.which("makeblastdb") or os.path.join(os.path.dirname(blastn), "makeblastdb")
    db = os.path.join(tmp, "te")
    subprocess.run([mk, "-in", te_lib, "-dbtype", "nucl", "-out", db], check=True, stdout=subprocess.DEVNULL)
    res = subprocess.run([blastn, "-task", "blastn", "-query", q, "-db", db, "-evalue", "1e-10", "-num_threads", str(threads),
                          "-max_target_seqs", "100", "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"],
                         check=True, capture_output=True, text=True).stdout
    hits = defaultdict(lambda: defaultdict(list))
    for line in res.splitlines():
        f = line.split("\t")
        hits[f[0]][f[1]].append((float(f[2]), int(f[3]), int(f[4]), int(f[5]), min(int(f[6]), int(f[7])), max(int(f[6]), int(f[7])),
                                 float(f[9]), int(f[10]), int(f[11])))
    out = {}
    for qid, subj in hits.items():
        def cov(iv, n):
            iv = sorted(iv)
            tot, cs, ce = 0, None, None
            for a, b in iv:
                if cs is None or a > ce + 1:
                    if cs is not None:
                        tot += ce - cs + 1
                    cs, ce = a, b
                else:
                    ce = max(ce, b)
            if cs is not None:
                tot += ce - cs + 1
            return tot / n
        scored = []
        for sid, hs in subj.items():
            bits = sum(h[6] for h in hs)
            qn, sn = hs[0][7], hs[0][8]
            pid = sum(h[0] * h[1] for h in hs) / sum(h[1] for h in hs)
            scored.append((bits, sid, cov([(h[2], h[3]) for h in hs], qn), cov([(h[4], h[5]) for h in hs], sn), pid,
                           min(h[2] for h in hs), max(h[3] for h in hs)))
        scored.sort(reverse=True)
        bits, sid, qc, sc, pid, qs, qe = scored[0]
        qn = next(iter(subj.values()))[0][7]
        allq = cov([(h[2], h[3]) for hs in subj.values() for h in hs], qn)
        out[qid] = {"te": sid, "te_cov": qc, "te_any_cov": allq, "te_subject_cov": sc, "te_identity": pid, "q_start": qs, "q_end": qe}
    return out


def kmer_te(seqs, te_lib, k=21):
    lib = []
    for name, s in read_fasta(te_lib):
        lib.append((name, {canonical(s[i:i + k]) for i in range(len(s) - k + 1)}, len(s)))
    out = {}
    for qid, s in seqs.items():
        if len(s) < 50:
            continue
        ks = [canonical(s[i:i + k]) for i in range(len(s) - k + 1)]
        best = None
        for name, kset, n in lib:
            hit = [i for i, x in enumerate(ks) if x in kset]
            if not hit:
                continue
            covered = np.zeros(len(s), dtype=bool)
            for i in hit:
                covered[i:i + k] = True
            qc = covered.mean()
            sc = min(1.0, len({ks[i] for i in hit}) / max(1, n - k + 1))
            if best is None or qc > best["te_cov"]:
                idx = np.flatnonzero(covered)
                best = {"te": name, "te_cov": float(qc), "te_any_cov": float(qc), "te_subject_cov": sc, "te_identity": None,
                        "q_start": int(idx[0]) + 1, "q_end": int(idx[-1]) + 1}
        if best:
            out[qid] = best
    return out


def period(seq, max_period=600):
    if len(seq) < 8:
        return None
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    for p in range(1, min(max_period, len(seq) // 2) + 1):
        if (a[p:] == a[:-p]).mean() >= 0.9:
            return p
    return None


def ltr_length(seq, w=50):
    n = len(seq)
    if n < 600:
        return 0
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    best = 0
    for o in range(0, 60, 6):
        seed = seq[o:o + 15]
        j = seq.find(seed, n // 2)
        while j != -1:
            d = j - o
            m = (a[:n - d] == a[d:]).astype(float)
            roll = np.convolve(m, np.ones(w) / w, "valid")
            good = np.flatnonzero(roll >= 0.8)
            if len(good) and good[0] <= 30:
                stop = np.flatnonzero(roll[good[0]:] < 0.8)
                L = (stop[0] if len(stop) else len(roll) - good[0]) + w - 1
                if 80 <= L <= d:
                    best = max(best, int(L))
            j = seq.find(seed, j + 1)
        if best:
            break
    return best


def tir_length(seq):
    if len(seq) < 200:
        return 0
    r = rc(seq)
    best = 0
    for L in range(10, min(300, len(seq) // 4)):
        if identity(seq[:L], r[:L]) >= 0.85:
            best = L
        elif L > best + 20:
            break
    return best if best >= 10 else 0


def tsd(context, te_s, te_e, search=25, slack=8, min_len=5):
    left = context[max(0, te_s - search):te_s + slack]
    right = context[max(0, te_e - slack):te_e + search]
    best = ""
    for i in range(len(left)):
        for j in range(len(right)):
            k = 0
            while i + k < len(left) and j + k < len(right) and left[i + k] == right[j + k]:
                k += 1
            if k > len(best) and len(left) - (i + k) <= slack and j <= slack:
                best = left[i:i + k]
    return best if len(best) >= min_len else ""


def kmer_positions(ref, k=21):
    idx = defaultdict(list)
    for i in range(len(ref) - k + 1):
        idx[ref[i:i + k]].append(i)
    return idx


def dup_hits(seq, idx, k=21, exclude=None):
    pos, n = [], 0
    for i in range(0, len(seq) - k + 1):
        km = seq[i:i + k]
        n += 1
        for p in idx.get(km, []) + idx.get(rc(km), []):
            if exclude is None or not (exclude[0] - k < p < exclude[1]):
                pos.append(p)
                break
    return (len(pos) / n if n else 0.0), (float(np.median(pos)) if pos else None)


def classify(seq, kind, te, site, ref_idx, ref_off, exclude, args):
    out = {"class": "unclassified", "te_family": ".", "te_class": ".", "te_cov": 0.0, "te_identity": None,
           "full_length": ".", "ltr_len": 0, "tir_len": 0, "dup_frac": 0.0, "dup_where": ".", "motif": "."}
    if not seq:
        return out
    if te:
        n, c = te_name_class(te["te"])
        out.update({"te_family": n, "te_class": c, "te_cov": te["te_cov"], "te_identity": te["te_identity"],
                    "full_length": "yes" if te["te_subject_cov"] >= 0.8 else "no"})
    frac, med = dup_hits(seq, ref_idx, exclude=None if exclude is None else (exclude[0] - ref_off, exclude[1] - ref_off))
    out["dup_frac"] = frac
    if med is not None and frac >= 0.5:
        out["dup_where"] = "adjacent" if abs(med + ref_off - site) <= len(seq) + 100 else f"{int(med + ref_off)}"
    p = period(seq)
    if p and p <= 6:
        out["motif"] = seq[:p]
    elif p:
        out["motif"] = f"{p}-bp unit"
    if out["dup_where"] == "." and not (p and p < len(seq) // 2):
        out["ltr_len"] = ltr_length(seq)
        out["tir_len"] = tir_length(seq)
    if te and te["te_cov"] >= args.te_min_cov:
        out["class"] = f"TE {out['te_class']}"
    elif p and p <= 6:
        out["class"] = f"microsatellite ({seq[:p]})n"
    elif out["dup_where"] == "adjacent" and p and p < len(seq):
        out["class"] = f"tandem repeat expansion ({len(seq) // p} x {p} bp)"
    elif out["dup_where"] == "adjacent":
        out["class"] = "tandem duplication"
    elif frac >= 0.8:
        out["class"] = "dispersed duplication"
    elif te and te["te_any_cov"] >= 0.2:
        out["class"] = f"TE-derived ({out['te_class']}, partial)"
    return out


def annotate(args):
    d = args.bubbles
    rows = read_tsv(os.path.join(d, "bubbles.tsv"))
    seqs = dict(read_fasta(os.path.join(d, "bubble_sequences.fa")))
    alt = {k.split("|")[0]: v for k, v in seqs.items() if k.split("|")[1] == "alt"}
    refs = {k.split("|")[0]: v for k, v in seqs.items() if k.split("|")[1] == "ref"}
    (rname, rseq), = list(read_fasta(os.path.join(d, "reference.fa")))
    contig, r0, r1 = parse_region(rname)
    roff = r0 - 1
    reg = None
    if args.regions:
        regs = read_tsv(args.regions)
        reg = next((r for r in regs if r["path"].split("#")[-1] == contig or r["genome"] == args.ref), None)
    models = load_models(args.gff, reg, contig, roff, len(rseq))
    log(f"{len(models)} gene models on the reference path")
    query = {}
    for r in rows:
        b = r["bubble"]
        if r["type"] in ("INS", "REPL", "INV"):
            query[b + "|alt"] = alt.get(b, "")
        if r["type"] in ("DEL", "REPL", "INV", "DUP"):
            query[b + "|ref"] = refs.get(b, "")
    te = {}
    if args.te_lib:
        blastn = args.blastn or os.environ.get("BLASTN") or shutil.which("blastn")
        if blastn and args.te_method != "kmer":
            with tempfile.TemporaryDirectory() as tmp:
                te = blast_te(query, args.te_lib, blastn, args.threads, tmp)
            log(f"TE search: blastn, {len(te)} of {len(query)} allele sequences hit the library")
        else:
            te = kmer_te(query, args.te_lib)
            log(f"TE search: 21-mer containment, {len(te)} of {len(query)} allele sequences hit the library")
    idx = kmer_positions(rseq)
    altwin = defaultdict(dict)
    jf = os.path.join(d, "bubble_junctions.tsv")
    if os.path.exists(jf):
        for j in read_tsv(jf):
            if j["allele"] == "alt":
                altwin[j["bubble"]][int(j["junction"])] = j
    for r in rows:
        b, t = r["bubble"], r["type"]
        a, e = int(r["pos"]), int(r["end"])
        s1, e1 = (a, a + 1) if t == "INS" else (a + 1, e)
        cat, detail, tss = context(models, s1, e1, args.promoter, args.downstream)
        r.update({"context": cat, "context_detail": detail,
                  "nearest_gene": tss[0] if tss else ".", "tss_distance": tss[1] if tss else "NA"})
        ac = classify(alt.get(b, ""), "alt", te.get(b + "|alt"), a, idx, roff, None, args) if t in ("INS", "REPL", "INV") else None
        rcl = classify(refs.get(b, ""), "ref", te.get(b + "|ref"), a, idx, roff, (a, e), args) if t in ("DEL", "REPL", "INV", "DUP") else None
        main = ac if t in ("INS", "REPL", "INV") else rcl
        r["alt_class"] = ac["class"] if ac else "."
        r["ref_class"] = rcl["class"] if rcl else "."
        for k in ("te_family", "te_class", "te_cov", "te_identity", "full_length", "ltr_len", "tir_len", "dup_frac", "dup_where", "motif"):
            r[k] = main[k]
        r["tsd"] = "."
        hit = te.get(b + "|alt")
        if t == "INS" and hit and main["class"].startswith("TE"):
            w = altwin.get(b, {})
            left = w[0]["window"][:int(w[0]["breakpoint"])] if 0 in w else rseq[max(0, a - roff - 40):a - roff]
            right = w[1]["window"][int(w[1]["breakpoint"]):] if 1 in w else rseq[a - roff:a - roff + 40]
            ctx = left + alt.get(b, "") + right
            r["tsd"] = tsd(ctx, len(left) + hit["q_start"] - 1, len(left) + hit["q_end"]) or "."
        if t == "DUP":
            cp = r.get("copies") or "2"
            r["sv_class"] = f"tandem duplication ({float(cp):g} copies)"
        elif t == "DEL":
            r["sv_class"] = "deletion" if rcl["class"] == "unclassified" else f"deletion of {rcl['class']}"
        elif t == "INS":
            if ac["class"] == "unclassified":
                r["sv_class"] = "insertion"
            elif ac["class"].startswith(("tandem repeat expansion", "microsatellite")):
                r["sv_class"] = ac["class"]
            else:
                r["sv_class"] = ac["class"] + " insertion"
        elif t == "INV":
            r["sv_class"] = "inversion"
        else:
            r["sv_class"] = "replacement" + ("" if ac["class"] == "unclassified" else f" by {ac['class']}")
        size = int(r["alt_len"]) if t == "INS" else int(r["ref_len"])
        where = detail.split(": ", 1)[1] if ": " in detail else detail
        r["label"] = f"{size_str(size)} {r['sv_class']} ({where})"
    for r in rows:
        r["linked"] = "."
        if r["type"] == "INS" and r["dup_where"] == "adjacent":
            for o in rows:
                if o["type"] == "DUP" and o["carriers"] == r["carriers"] and int(o["pos"]) - 50 <= int(r["pos"]) <= int(o["end"]) + 50:
                    r["linked"] = o["bubble"]
                    o["linked"] = r["bubble"]
                    r["label"] = r["label"].replace("tandem duplication insertion", f"tandem copy of {o['bubble']}")
    fields = list(rows[0].keys()) if rows else []
    out = args.out or os.path.join(d, "bubbles.annotated.tsv")
    write_tsv(out, rows, fields)
    for r in rows:
        log(f"{r['bubble']} {r['type']:<4} {r['n_carriers']:>3} carriers  {r['label']}")
    log(f"written {out}")
