import os
import re
import statistics
from collections import defaultdict

from .common import canonical, die, log, mkdir, opener, parse_path_name, rc, read_tsv, write_fasta, write_tsv


class Graph:
    def __init__(self, gfa):
        self.seq = {}
        self.paths = {}
        with opener(gfa) as fh:
            for line in fh:
                if line.startswith("S\t"):
                    f = line.rstrip("\n").split("\t")
                    self.seq[f[1]] = f[2].upper()
                elif line.startswith("P\t"):
                    f = line.rstrip("\n").split("\t")
                    self.paths[f[1]] = [(s[:-1], s[-1]) for s in f[2].split(",") if s]
                elif line.startswith("W\t"):
                    f = line.rstrip("\n").split("\t")
                    name = f"{f[1]}#{f[2]}#{f[3]}:{f[4]}-{f[5]}"
                    self.paths[name] = [(n, "+" if d == ">" else "-") for d, n in re.findall(r"([<>])([^<>]+)", f[6])]
        if not self.paths:
            die(f"no P or W lines in {gfa}")
        self.info = {p: parse_path_name(p) for p in self.paths}

    def node(self, n, o):
        s = self.seq[n]
        return s if o == "+" else rc(s)

    def path_seq(self, p):
        return "".join(self.node(n, o) for n, o in self.paths[p])

    def genomes(self):
        out = defaultdict(list)
        for p, i in self.info.items():
            out[i["sample"]].append(p)
        return out


class Reference:
    def __init__(self, graph, ref):
        cands = [p for p, i in graph.info.items() if i["sample"] == ref]
        if not cands:
            die(f"no path for reference {ref} in the graph")
        self.path = max(cands, key=lambda p: sum(len(graph.seq[n]) for n, _ in graph.paths[p]))
        self.contig = graph.info[self.path]["contig"]
        self.offset = graph.info[self.path]["offset"]
        self.pos = defaultdict(list)
        self.orient = {}
        x = self.offset
        for n, o in graph.paths[self.path]:
            self.pos[n].append(x)
            self.orient.setdefault(n, o)
            x += len(graph.seq[n])
        self.seq = graph.path_seq(self.path)
        self.end = x

    def at(self, a, b):
        a, b = max(a, self.offset), min(b, self.end)
        return self.seq[a - self.offset:b - self.offset]


def walk(graph, ref, path, min_sv):
    seq, events, bmap = [], [], {}
    prev, pending, flipped, pstart, cur = None, 0, 0, 0, 0
    first = None
    for n, o in graph.paths[path]:
        s = graph.node(n, o)
        cand = ref.pos.get(n)
        if cand is not None and prev is not None and o != ref.orient[n]:
            if pending == 0:
                pstart = cur
            pending += len(s)
            flipped += len(s)
            seq.append(s)
            cur += len(s)
            continue
        if cand is None:
            if prev is not None:
                if pending == 0:
                    pstart = cur
                pending += len(s)
            seq.append(s)
            cur += len(s)
            continue
        x = cand[0] if prev is None or len(cand) == 1 else min(cand, key=lambda v: abs(v - prev))
        if prev is None:
            first = x
        else:
            d = x - prev
            ev = None
            if d > min_sv and pending >= min_sv:
                ev = ("INV" if flipped >= 0.5 * pending else "REPL", prev, x, pending, pstart, cur)
            elif d > min_sv:
                ev = ("DEL", prev, x, 0, cur, cur)
            elif pending >= min_sv:
                ev = ("INS", prev, prev, pending, pstart, cur)
            if ev:
                events.append(dict(zip(("type", "a", "b", "alt_len", "i0", "i1"), ev)))
            if d < -min_sv:
                events.append({"type": "DUP", "a": x, "b": prev, "alt_len": prev - x, "i0": cur, "i1": cur})
        bmap.setdefault(x, cur)
        bmap.setdefault(x + len(s), cur + len(s))
        seq.append(s)
        cur += len(s)
        prev, pending, flipped = x + len(s), 0, 0
    span = (first, prev) if first is not None else None
    return "".join(seq), events, span, bmap


def kmer_set(seq, k=15):
    return {canonical(seq[i:i + k]) for i in range(0, max(0, len(seq) - k + 1))}


def jaccard(x, y):
    if not x or not y:
        return 0.0
    return len(x & y) / len(x | y)


def size_of(e):
    if e["type"] == "INS":
        return e["alt_len"]
    if e["type"] in ("REPL", "INV"):
        return max(e["b"] - e["a"], e["alt_len"])
    return e["b"] - e["a"]


def cluster_events(events, tol, min_ratio, min_jaccard):
    clusters = []
    for e in sorted(events, key=lambda e: (e["type"], e["a"], e["b"])):
        hit = None
        for c in clusters:
            if c["type"] != e["type"] or abs(c["a"] - e["a"]) > tol or abs(c["b"] - e["b"]) > tol:
                continue
            s1, s2 = size_of(c), size_of(e)
            if min(s1, s2) / max(s1, s2, 1) < min_ratio:
                continue
            if e["type"] in ("INS", "REPL", "INV") and jaccard(c["kmers"], e["kmers"]) < min_jaccard:
                continue
            hit = c
            break
        if hit is None:
            clusters.append({"type": e["type"], "a": e["a"], "b": e["b"], "members": [e], "kmers": e["kmers"]})
        else:
            hit["members"].append(e)
            hit["a"] = int(statistics.median(m["a"] for m in hit["members"]))
            hit["b"] = int(statistics.median(m["b"] for m in hit["members"]))
    return clusters


def overlaps(a1, b1, a2, b2, pad=0):
    return a1 - pad <= b2 and a2 - pad <= b1


def call(args):
    out = mkdir(args.out)
    graph = Graph(args.gfa)
    ref = Reference(graph, args.ref)
    log(f"graph: {len(graph.seq)} nodes, {len(graph.paths)} paths; reference path {ref.path} ({len(ref.seq)} bp, offset {ref.offset})")
    walks, events, bmaps = {}, [], {}
    spans = defaultdict(list)
    for p in graph.paths:
        g = graph.info[p]["sample"]
        s, evs, span, bm = walk(graph, ref, p, args.min_sv)
        walks[p] = s
        bmaps[p] = bm
        if span:
            spans[g].append(span)
        for e in evs:
            e.update({"genome": g, "path": p})
            e["alt_seq"] = s[e["i0"]:e["i1"]]
            e["kmers"] = kmer_set(e["alt_seq"]) if e["type"] in ("INS", "REPL", "INV") else set()
            events.append(e)
    if not events:
        log("no structural differences >= min_sv relative to the reference")
    clusters = cluster_events(events, args.tolerance, args.min_size_ratio, args.min_jaccard)
    clusters.sort(key=lambda c: (c["a"], c["b"], c["type"]))
    genomes = sorted(graph.genomes())
    regions = {r["genome"]: r for r in read_tsv(args.regions)} if args.regions else {}
    rref = regions.get(args.ref)
    rows, mrows, jrows = [], {g: {"genome": g} for g in genomes}, []
    seqfh = open(os.path.join(out, "bubble_sequences.fa"), "w")
    F = args.window_flank
    for i, c in enumerate(clusters, 1):
        bid = f"{args.prefix}{i:03d}"
        c["id"] = bid
        carriers = sorted({m["genome"] for m in c["members"]})
        rep = min(c["members"], key=lambda m: (abs(m["a"] - c["a"]) + abs(m["b"] - c["b"]), m["genome"] == args.ref, m["genome"]))
        c["rep"] = rep
        status = {}
        for g in genomes:
            if g in carriers:
                status[g] = "1"
            elif not any(s[0] <= c["a"] - args.cover_pad and s[1] >= c["b"] + args.cover_pad for s in spans.get(g, [])):
                status[g] = "."
            elif any(e["genome"] == g and (e["type"] == "DUP") == (c["type"] == "DUP") and overlaps(c["a"], c["b"], e["a"], e["b"], 20) for e in events):
                status[g] = "x"
            else:
                status[g] = "0"
            mrows[g][bid] = status[g]
        c["status"] = status
        alt_lens = [m["alt_len"] for m in c["members"]]
        row = {"bubble": bid, "type": c["type"], "pos": c["a"], "end": c["b"] if c["type"] != "INS" else c["a"],
               "ref_len": c["b"] - c["a"], "alt_len": int(statistics.median(alt_lens)),
               "n_carriers": len(carriers), "n_reference": sum(v == "0" for v in status.values()),
               "n_other_allele": sum(v == "x" for v in status.values()), "n_missing": sum(v == "." for v in status.values()),
               "carriers": ",".join(carriers), "representative": rep["genome"]}
        if c["type"] == "DUP":
            row["copies"] = statistics.median([1 + sum(1 for e in events if e["genome"] == g and e["type"] == "DUP" and overlaps(c["a"], c["b"], e["a"], e["b"])) for g in carriers])
        if rref:
            s0, e0, ori = int(rref["start"]), int(rref["end"]), rref["orientation"]
            if ori == "+":
                row.update({"chrom": rref["contig"], "chrom_start": s0 + c["a"], "chrom_end": s0 + max(c["a"], c["b"] - 1)})
            else:
                row.update({"chrom": rref["contig"], "chrom_start": e0 - max(c["a"], c["b"] - 1), "chrom_end": e0 - c["a"]})
        rows.append(row)
        write_fasta(seqfh, f"{bid}|ref|{args.ref}|{ref.contig}:{c['a'] + 1}-{c['b']}", ref.at(c["a"], c["b"]) if c["b"] > c["a"] else "")
        write_fasta(seqfh, f"{bid}|alt|{rep['genome']}", rep["alt_seq"])
        seen = set()
        for m in sorted(c["members"], key=lambda m: m is not rep):
            cs = walks[m["path"]]
            alt_j = [m["i0"], m["i1"]] if c["type"] in ("INS", "REPL", "INV") else [m["i0"]]
            for k, j in enumerate(alt_j):
                if j - F >= 0 and j + F <= len(cs) and ("alt", k, cs[j - F:j + F]) not in seen:
                    seen.add(("alt", k, cs[j - F:j + F]))
                    jrows.append({"bubble": bid, "type": c["type"], "allele": "alt", "junction": k, "source": m["genome"],
                                  "ref_pos": c["a"] if k == 0 else c["b"], "window": cs[j - F:j + F], "breakpoint": F})
        ref_j = [c["a"]] if c["type"] == "INS" else [c["a"], c["b"]]
        for k, j in enumerate(ref_j):
            w = ref.at(j - F, j + F)
            if len(w) == 2 * F:
                seen.add(("ref", k, w))
                jrows.append({"bubble": bid, "type": c["type"], "allele": "ref", "junction": k, "source": args.ref,
                              "ref_pos": j, "window": w, "breakpoint": F})
            for p in graph.paths:
                g = graph.info[p]["sample"]
                if status.get(g) != "0" or g == args.ref or j not in bmaps[p]:
                    continue
                i, cs = bmaps[p][j], walks[p]
                w = cs[i - F:i + F] if i - F >= 0 and i + F <= len(cs) else ""
                if len(w) == 2 * F and ("ref", k, w) not in seen:
                    seen.add(("ref", k, w))
                    jrows.append({"bubble": bid, "type": c["type"], "allele": "ref", "junction": k, "source": g,
                                  "ref_pos": j, "window": w, "breakpoint": F})
    seqfh.close()
    for r in rows:
        r["overlaps"] = ",".join(o["bubble"] for o in rows if o is not r and overlaps(r["pos"], r["end"], o["pos"], o["end"])) or "."
    fields = ["bubble", "type", "pos", "end", "ref_len", "alt_len", "copies", "n_carriers", "n_reference", "n_other_allele",
              "n_missing", "carriers", "representative", "overlaps"] + (["chrom", "chrom_start", "chrom_end"] if rref else [])
    write_tsv(os.path.join(out, "bubbles.tsv"), rows, fields)
    write_tsv(os.path.join(out, "assembly_matrix.tsv"), [mrows[g] for g in genomes], ["genome"] + [r["bubble"] for r in rows])
    write_tsv(os.path.join(out, "bubble_junctions.tsv"), jrows, ["bubble", "type", "allele", "junction", "source", "ref_pos", "breakpoint", "window"])
    write_tsv(os.path.join(out, "path_events.tsv"),
              [{"genome": e["genome"], "path": e["path"], "type": e["type"], "pos": e["a"], "end": e["b"], "ref_len": e["b"] - e["a"],
                "alt_len": e["alt_len"], "bubble": next(c["id"] for c in clusters if e in c["members"])} for e in events],
              ["genome", "path", "type", "pos", "end", "ref_len", "alt_len", "bubble"])
    with open(os.path.join(out, "reference.fa"), "w") as fh:
        write_fasta(fh, f"{ref.contig}:{ref.offset + 1}-{ref.end}", ref.seq)
    with open(os.path.join(out, "genome_paths.fa"), "w") as fh:
        for g, ps in graph.genomes().items():
            for p in ps:
                write_fasta(fh, p, walks[p])
    counts = defaultdict(int)
    for r in rows:
        counts[r["type"]] += 1
    log(f"{len(rows)} bubbles ({', '.join(f'{k} {v}' for k, v in sorted(counts.items()))}) from {len(events)} path events in {len(genomes)} genomes")
