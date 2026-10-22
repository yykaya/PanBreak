import math
import os
from collections import Counter, defaultdict
from itertools import combinations

from .common import log, mkdir, read_tsv, write_tsv


def wilson(k, n, z=1.96):
    if n == 0:
        return None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def lchoose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_greater(a, b, c, d):
    n1, n2, m = a + b, c + d, a + c
    lo, hi = max(0, m - n2), min(n1, m)
    denom = lchoose(n1 + n2, m)
    p = sum(math.exp(lchoose(n1, x) + lchoose(n2, m - x) - denom) for x in range(a, hi + 1))
    return min(1.0, p)


def bh(ps):
    idx = sorted(range(len(ps)), key=lambda i: ps[i])
    q, prev = [0.0] * len(ps), 1.0
    for rank, i in reversed(list(enumerate(idx, 1))):
        prev = min(prev, ps[i] * len(ps) / rank)
        q[i] = prev
    return q


def hudson(p1, n1, p2, n2):
    if n1 < 2 or n2 < 2:
        return None, None
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num, den


def stats(args):
    out = mkdir(args.out)
    gts = read_tsv(args.genotypes)
    meta = {r["sample"]: r for r in read_tsv(args.samples)} if args.samples else {}
    col = args.group_col
    pop = {s: (meta[s].get(col) or "NA") if s in meta else "unassigned" for s in {r["sample"] for r in gts}}
    if args.drop_unassigned:
        pop = {s: p for s, p in pop.items() if p != "unassigned"}
    btab = None
    if args.bubbles:
        f = os.path.join(args.bubbles, "bubbles.annotated.tsv")
        btab = {r["bubble"]: r for r in read_tsv(f if os.path.exists(f) else os.path.join(args.bubbles, "bubbles.tsv"))}
    pops = sorted(set(pop.values()))
    by = defaultdict(list)
    for r in gts:
        if r["sample"] in pop:
            by[r["bubble"]].append(r)
    carried = defaultdict(set)
    for r in gts:
        if r["FT"] == "PASS" and r["carrier"] == "1":
            carried[r["sample"]].add(r["bubble"])
    for r in gts:
        r["other"] = False
        if r["FT"] != "PASS" and btab and not args.no_overlap_rescue:
            ov = set(btab.get(r["bubble"], {}).get("overlaps", ".").split(","))
            if ov & carried[r["sample"]]:
                r["other"] = True
    ploidy = args.ploidy
    rows, enrich, fst_rows = [], [], []
    for b in sorted(by):
        rs = by[b]
        model = rs[0]["model"]
        groups = {"ALL": rs}
        for p in pops:
            groups[p] = [r for r in rs if pop[r["sample"]] == p]
        stat = {}
        for g, grs in groups.items():
            n = len(grs)
            called = [r for r in grs if r["FT"] == "PASS"]
            other = sum(r["other"] for r in grs)
            nc = len(called) + other
            car = sum(r["carrier"] == "1" for r in called)
            hom = sum(r["call"] == "hom_alt" for r in called)
            het = sum(r["call"] == "het" for r in called)
            if model == "presence" or ploidy == 1:
                ac, an = car, nc
            else:
                ac, an = 2 * hom + het, 2 * nc
            af = ac / an if an else None
            lo, hi = wilson(ac, an)
            exp_het = 2 * af * (1 - af) if af is not None else None
            obs_het = het / nc if nc and model != "presence" and ploidy == 2 else None
            hom_ref_n = sum(r["call"] == "hom_ref" for r in called) + other
            fis = 1 - obs_het / exp_het if exp_het and obs_het is not None else None
            stat[g] = {"n": nc, "ac": ac, "an": an, "af": af, "car": car}
            rows.append({"bubble": b, "population": g, "samples": n, "called": nc, "other_allele": other, "call_rate": round(nc / n, 3) if n else None,
                         "hom_ref": hom_ref_n, "het": het, "hom_alt": hom,
                         "carriers": car, "carrier_freq": round(car / nc, 4) if nc else None,
                         "allele_count": ac, "allele_number": an, "AF": round(af, 4) if af is not None else None,
                         "AF_lo95": round(lo, 4) if lo is not None else None, "AF_hi95": round(hi, 4) if hi is not None else None,
                         "exp_het": round(exp_het, 4) if exp_het is not None else None,
                         "obs_het": round(obs_het, 4) if obs_het is not None else None, "F_IS": round(fis, 3) if fis is not None else None})
        tot = stat["ALL"]
        for p in pops:
            a, n = stat[p]["car"], stat[p]["n"]
            c, m = tot["car"] - a, tot["n"] - n
            if n == 0:
                continue
            orr = ((a + 0.5) * (m - c + 0.5)) / ((n - a + 0.5) * (c + 0.5))
            enrich.append({"bubble": b, "population": p, "carriers": a, "called": n, "carriers_elsewhere": c, "called_elsewhere": m,
                           "odds_ratio": round(orr, 3), "p_enrichment": fisher_greater(a, n - a, c, m - c)})
        for p1, p2 in combinations(pops, 2):
            s1, s2 = stat[p1], stat[p2]
            if s1["af"] is None or s2["af"] is None:
                continue
            num, den = hudson(s1["af"], s1["an"], s2["af"], s2["an"])
            if num is None:
                continue
            fst_rows.append({"bubble": b, "pop1": p1, "pop2": p2, "AF1": round(s1["af"], 4), "AF2": round(s2["af"], 4),
                             "fst_num": num, "fst_den": den, "Fst": round(num / den, 4) if den > 0 else None})
    qs = bh([e["p_enrichment"] for e in enrich]) if enrich else []
    for e, q in zip(enrich, qs):
        e["q_enrichment"] = q
    mean_fst = []
    for p1, p2 in combinations(pops, 2):
        fr = [f for f in fst_rows if f["pop1"] == p1 and f["pop2"] == p2]
        num, den = sum(f["fst_num"] for f in fr), sum(f["fst_den"] for f in fr)
        if den > 0:
            mean_fst.append({"bubble": "ALL", "pop1": p1, "pop2": p2, "Fst": round(num / den, 4), "bubbles": len(fr)})
    if btab:
        for r in rows:
            bi = btab.get(r["bubble"], {})
            r["type"] = bi.get("type", "")
            r["label"] = bi.get("label", "")
    write_tsv(os.path.join(out, "population_frequencies.tsv"), rows)
    write_tsv(os.path.join(out, "population_enrichment.tsv"), sorted(enrich, key=lambda e: e["p_enrichment"]))
    write_tsv(os.path.join(out, "pairwise_fst.tsv"), [{k: v for k, v in f.items() if not k.startswith("fst_")} for f in fst_rows] + mean_fst,
              ["bubble", "pop1", "pop2", "AF1", "AF2", "Fst", "bubbles"])
    wide = []
    for b in sorted(by):
        w = {"bubble": b}
        if btab:
            w["label"] = btab.get(b, {}).get("label", "")
        for r in rows:
            if r["bubble"] == b:
                w[r["population"]] = r["AF"] if r["AF"] is not None else "NA"
        wide.append(w)
    write_tsv(os.path.join(out, "af_matrix.tsv"), wide, ["bubble"] + (["label"] if btab else []) + ["ALL"] + pops)
    haplotypes(gts, pop, pops, out, args, btab or {})
    log(f"population statistics for {len(by)} bubbles in {len(pops)} groups written to {out}")


def haplotypes(gts, pop, pops, out, args, btab):
    by = defaultdict(dict)
    for r in gts:
        if r["sample"] in pop:
            by[r["sample"]][r["bubble"]] = r
    bubbles = sorted({r["bubble"] for r in gts})
    asm = {}
    if args.bubbles and os.path.exists(os.path.join(args.bubbles, "assembly_matrix.tsv")):
        for r in read_tsv(os.path.join(args.bubbles, "assembly_matrix.tsv")):
            if all(r.get(b) in ("0", "1", "x") for b in bubbles):
                key = "+".join(b for b in bubbles if r[b] == "1") or "REF"
                asm.setdefault(key, []).append(r["genome"])
    srows = []
    for s, d in sorted(by.items()):
        parts, complete = [], True
        carried = {b for b, r in d.items() if r["FT"] == "PASS" and r["carrier"] == "1"}
        for b in bubbles:
            r = d.get(b)
            if r is None or r["FT"] != "PASS":
                ov = set(btab.get(b, {}).get("overlaps", ".").split(","))
                if not ov & carried:
                    complete = False
                continue
            if r["carrier"] == "1":
                parts.append(b + ("(het)" if r["call"] == "het" else ""))
        hap = "+".join(parts) or "REF"
        srows.append({"sample": s, "population": pop[s], "haplotype": hap if complete else hap + "?", "complete": "yes" if complete else "no",
                      "matching_assemblies": ",".join(asm.get(hap, [])) or "."})
    write_tsv(os.path.join(out, "sample_haplotypes.tsv"), srows)
    cnt = Counter((r["haplotype"], r["population"]) for r in srows if r["complete"] == "yes")
    tot = Counter(r["population"] for r in srows if r["complete"] == "yes")
    haps = sorted({h for h, _ in cnt}, key=lambda h: -sum(v for (hh, _), v in cnt.items() if hh == h))
    hrows = []
    for h in haps:
        row = {"haplotype": h, "matching_assemblies": ",".join(asm.get(h, [])) or ".", "total": sum(v for (hh, _), v in cnt.items() if hh == h)}
        for p in pops:
            row[p] = cnt[(h, p)]
            row[p + "_freq"] = round(cnt[(h, p)] / tot[p], 4) if tot[p] else None
        hrows.append(row)
    write_tsv(os.path.join(out, "haplotype_frequencies.tsv"), hrows,
              ["haplotype", "matching_assemblies", "total"] + [x for p in pops for x in (p, p + "_freq")])
