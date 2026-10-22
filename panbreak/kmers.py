import os
from collections import defaultdict
from multiprocessing import Pool

from .common import canonical, log, mkdir, rc, read_fasta, read_tsv, write_fasta, write_tsv


def junction_kmers(window, bp, k, min_side):
    out = []
    for st in range(max(0, bp - k + min_side), min(len(window) - k, bp - min_side) + 1):
        km = window[st:st + k]
        if "N" not in km:
            out.append((km, bp - st))
    return out


def genome_sequences(paths_fa, extra):
    seqs = defaultdict(list)
    for name, s in read_fasta(paths_fa):
        seqs[name.split("#")[0]].append(s)
    for f in extra or []:
        for name, s in read_fasta(f):
            seqs[name.split("#")[0]].append(s)
    return {g: "N".join(v) for g, v in seqs.items()}


def count_background(args):
    path, kmers, k = args
    want = set(kmers) | {rc(x) for x in kmers}
    counts = defaultdict(int)
    for _, s in read_fasta(path):
        for i in range(len(s) - k + 1):
            x = s[i:i + k]
            if x in want:
                counts[canonical(x)] += 1
    return path, counts


def build(args):
    d = args.bubbles
    out = mkdir(args.out or d)
    btab = os.path.join(d, "bubbles.annotated.tsv")
    bubbles = {r["bubble"]: r for r in read_tsv(btab if os.path.exists(btab) else os.path.join(d, "bubbles.tsv"))}
    matrix = {r["genome"]: r for r in read_tsv(os.path.join(d, "assembly_matrix.tsv"))}
    genomes = genome_sequences(os.path.join(d, "genome_paths.fa"), args.fasta)
    rows, seen = [], set()
    for j in read_tsv(os.path.join(d, "bubble_junctions.tsv")):
        for km, off in junction_kmers(j["window"], int(j["breakpoint"]), args.k, args.min_side):
            key = (j["bubble"], j["allele"], j["junction"], km)
            if key in seen:
                continue
            seen.add(key)
            rows.append({"bubble": j["bubble"], "type": j["type"], "allele": j["allele"], "junction": int(j["junction"]),
                         "ref_pos": j["ref_pos"], "source": j["source"], "kmer": km, "canonical": canonical(km), "offset": off})
    log(f"{len(rows)} candidate junction {args.k}-mers from {len(bubbles)} bubbles")
    present = {}
    uniq = {r["kmer"] for r in rows}
    for g, s in genomes.items():
        present[g] = {km for km in uniq if km in s or rc(km) in s}
    bg = defaultdict(int)
    if args.background:
        with Pool(min(args.threads, len(args.background))) as pool:
            for path, counts in pool.imap_unordered(count_background, [(p, list(uniq), args.k) for p in args.background]):
                for km, n in counts.items():
                    bg[km] = max(bg[km], n)
                log(f"background {os.path.basename(path)} scanned")
    presence = {b for b, br in bubbles.items() if br["type"] == "DUP" or (br["type"] == "INS" and br.get("dup_where") == "adjacent")}
    kept, dropped = [], defaultdict(int)
    for r in rows:
        b = r["bubble"]
        status = {g: matrix[g].get(b, ".") for g in matrix}
        car = [g for g, v in status.items() if v == "1" and g in present]
        non = [g for g, v in status.items() if v in ("0", "x") and g in present]
        refs = [g for g, v in status.items() if v == "0" and g in present]
        km = r["kmer"]
        r["n_carrier_asm"] = sum(km in present[g] for g in car)
        r["n_noncarrier_asm"] = sum(km in present[g] for g in non)
        if r["allele"] == "alt":
            if r["n_noncarrier_asm"]:
                dropped["alt k-mer in a non-carrier assembly"] += 1
                continue
            if not r["n_carrier_asm"]:
                dropped["alt k-mer in no carrier assembly"] += 1
                continue
        else:
            if b not in presence and r["n_carrier_asm"]:
                dropped["ref k-mer in a carrier assembly"] += 1
                continue
            if refs and not any(km in present[g] for g in refs):
                dropped["ref k-mer in no reference-allele assembly"] += 1
                continue
        if args.background and bg.get(r["canonical"], 0) > args.max_copies:
            dropped["multi-copy in background genome"] += 1
            continue
        kept.append(r)
    for k, v in dropped.items():
        log(f"dropped {v} k-mers: {k}")
    by = defaultdict(list)
    for r in kept:
        by[r["bubble"]].append(r)
    summ, asm_rows = [], {g: {"genome": g} for g in sorted(matrix)}
    for b, br in bubbles.items():
        ks = by.get(b, [])
        alt = {r["kmer"] for r in ks if r["allele"] == "alt"}
        ref = {r["kmer"] for r in ks if r["allele"] == "ref"}
        ja = len({r["junction"] for r in ks if r["allele"] == "alt"})
        jr = len({r["junction"] for r in ks if r["allele"] == "ref"})
        model = "presence" if b in presence else "biallelic"
        tp = fn = fp = tn = 0
        for g in sorted(matrix):
            st = matrix[g].get(b, ".")
            if g not in present:
                asm_rows[g][b] = "."
                continue
            ha, hr = len(alt & present[g]), len(ref & present[g])
            call = "1" if ha else ("0" if hr else ".")
            asm_rows[g][b] = call
            if st == "1":
                tp += call == "1"
                fn += call != "1"
            elif st in ("0", "x"):
                fp += call == "1"
                tn += call != "1"
        reason = []
        if not ja:
            reason.append("no specific alt k-mer")
        if not jr:
            reason.append("no specific ref k-mer" if model == "biallelic" else "no local depth k-mer")
        summ.append({"bubble": b, "type": br["type"], "model": model, "alt_junctions": ja, "ref_junctions": jr,
                     "alt_kmers": len(alt), "ref_kmers": len(ref),
                     "het_alt_fraction": round(ja / (ja + jr), 3) if model == "biallelic" and ja + jr else "NA",
                     "asm_carriers": int(br["n_carriers"]), "asm_recall": round(tp / (tp + fn), 3) if tp + fn else "NA",
                     "asm_false_carriers": fp, "genotypable": "no" if reason else "yes", "reason": "; ".join(reason) or "."})
    write_tsv(os.path.join(out, "junction_kmers.tsv"), kept,
              ["bubble", "type", "allele", "junction", "ref_pos", "source", "offset", "kmer", "canonical", "n_carrier_asm", "n_noncarrier_asm"])
    with open(os.path.join(out, "junction_kmers.fa"), "w") as fh:
        for i, r in enumerate(kept):
            write_fasta(fh, f"{r['bubble']}|{r['allele']}|j{r['junction']}|{i}", r["kmer"])
    write_tsv(os.path.join(out, "junction_summary.tsv"), summ)
    write_tsv(os.path.join(out, "assembly_kmer_genotypes.tsv"), [asm_rows[g] for g in sorted(asm_rows)], ["genome"] + list(bubbles))
    ok = sum(s["genotypable"] == "yes" for s in summ)
    log(f"{len(kept)} specific k-mers kept; {ok}/{len(summ)} bubbles genotypable")
    for s in summ:
        log(f"{s['bubble']} {s['type']:<4} alt {s['alt_kmers']:>3} ({s['alt_junctions']} jx)  ref {s['ref_kmers']:>3} ({s['ref_junctions']} jx)  "
            f"assembly recall {s['asm_recall']}  false carriers {s['asm_false_carriers']}  {s['genotypable']}")
