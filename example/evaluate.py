import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from panbreak.common import read_tsv, write_tsv


def main():
    ap = argparse.ArgumentParser(description="compare panbreak genotypes with the simulated haplotypes")
    ap.add_argument("--samples", default="data/samples.tsv")
    ap.add_argument("--matrix", default="results/03_bubbles/assembly_matrix.tsv")
    ap.add_argument("--bubbles", default="results/03_bubbles/bubbles.annotated.tsv")
    ap.add_argument("--genotypes", default="results/05_genotypes/genotypes.tsv")
    ap.add_argument("-o", "--out", default="results/evaluation.tsv")
    a = ap.parse_args()
    mat = {r["genome"]: r for r in read_tsv(a.matrix)}
    samples = {r["sample"]: r for r in read_tsv(a.samples)}
    binfo = {r["bubble"]: r for r in read_tsv(a.bubbles)}
    st = defaultdict(lambda: defaultdict(int))
    for g in read_tsv(a.genotypes):
        s, b = samples[g["sample"]], g["bubble"]
        hs = [mat[s["haplotype1"]][b], mat[s["haplotype2"]][b]]
        truth = sum(h == "1" for h in hs)
        other = any(h == "x" for h in hs)
        x = st[b]
        x["samples"] += 1
        x["true_carriers"] += truth > 0
        if g["FT"] != "PASS":
            x["no_call"] += 1
            x["no_call_other_allele"] += other and truth == 0
            continue
        x["called"] += 1
        car = g["carrier"] == "1"
        if g["model"] == "presence":
            x["concordant"] += car == (truth > 0)
        else:
            x["concordant"] += int(g["dosage"]) == truth
        x["carriers_found"] += car and truth > 0
        x["false_carriers"] += car and truth == 0
        x["missed_carriers"] += (not car) and truth > 0
    rows = []
    for b in sorted(st, key=lambda b: int(binfo[b]["pos"])):
        x = st[b]
        rows.append({"bubble": b, "type": binfo[b]["type"], "label": binfo[b]["label"], **x,
                     "call_rate": round(x["called"] / x["samples"], 3),
                     "call_rate_excl_other_allele": round(x["called"] / max(1, x["samples"] - x["no_call_other_allele"]), 3),
                     "concordance": round(x["concordant"] / x["called"], 4) if x["called"] else "NA",
                     "sensitivity": round(x["carriers_found"] / max(1, x["carriers_found"] + x["missed_carriers"]), 4)})
    tot = defaultdict(int)
    for x in st.values():
        for k, v in x.items():
            tot[k] += v
    rows.append({"bubble": "ALL", "type": "", "label": "", **tot, "call_rate": round(tot["called"] / tot["samples"], 3),
                 "call_rate_excl_other_allele": round(tot["called"] / max(1, tot["samples"] - tot["no_call_other_allele"]), 3),
                 "concordance": round(tot["concordant"] / tot["called"], 4),
                 "sensitivity": round(tot["carriers_found"] / max(1, tot["carriers_found"] + tot["missed_carriers"]), 4)})
    fields = ["bubble", "type", "label", "samples", "true_carriers", "called", "no_call", "no_call_other_allele", "call_rate",
              "call_rate_excl_other_allele", "concordant", "concordance", "carriers_found", "missed_carriers", "false_carriers", "sensitivity"]
    write_tsv(a.out, rows, fields)
    for r in rows:
        print(f"{r['bubble']:<5} {r['type']:<4} called {r['called']:>3}/{r['samples']}  concordance {r['concordance']}  "
              f"sensitivity {r['sensitivity']}  false carriers {r['false_carriers']}  no-call (other allele) {r['no_call_other_allele']}")


if __name__ == "__main__":
    main()
