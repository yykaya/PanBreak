import os
import re
from bisect import bisect_right
from collections import defaultdict

from .common import die, log, opener, parse_path_name, run, tool

STEP = re.compile(r"([<>])([^<>]+)")


def tags(fields):
    out = {}
    for t in fields:
        k, _, v = t.partition(":")
        out[k] = v.partition(":")[2]
    return out


class GFA:
    def __init__(self, path):
        self.order, self.seq, self.tag = [], {}, {}
        self.links, self.paths = set(), {}
        walks = defaultdict(list)
        with opener(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == "S":
                    self.order.append(f[1])
                    self.seq[f[1]] = f[2].upper()
                    self.tag[f[1]] = tags(f[3:])
                elif f[0] == "L":
                    self.links.add((f[1], f[2], f[3], f[4]))
                elif f[0] == "P":
                    self.paths[f[1]] = [(s[:-1], s[-1]) for s in f[2].split(",") if s]
                elif f[0] == "W":
                    steps = [(n, "+" if d == ">" else "-") for d, n in STEP.findall(f[6])]
                    start = int(f[4]) if f[4] != "*" else 0
                    end = int(f[5]) if f[5] != "*" else start + sum(len(self.seq.get(n, "")) for n, _ in steps)
                    walks[(f[1], f[2], f[3])].append((start, end, steps))
        for (s, h, c), ws in walks.items():
            for start, end, steps in ws:
                name = f"{s}#{h}#{c}" if len(ws) == 1 and start == 0 else f"{s}#{h}#{c}:{start}-{end}"
                self.paths[name] = steps
        self.index = defaultdict(list)
        for n in self.order:
            t = self.tag[n]
            if "SN" in t and "SO" in t:
                self.index[t["SN"]].append((int(t["SO"]), n))
        for v in self.index.values():
            v.sort()

    def stable(self, sn, a, b, strand):
        segs = self.index.get(sn)
        if not segs:
            die(f"GAF refers to {sn}, which is not a stable sequence (SN tag) of the graph")
        starts = [s for s, _ in segs]
        i = max(0, bisect_right(starts, a) - 1)
        out = []
        while i < len(segs) and segs[i][0] < b:
            s, n = segs[i]
            if s + len(self.seq[n]) > a:
                out.append((n, "+"))
            i += 1
        return out if strand == "+" else [(n, "-" if o == "+" else "+") for n, o in reversed(out)]

    def gaf_steps(self, f):
        path = f[5]
        if path[0] not in "<>":
            return self.stable(path, int(f[7]), int(f[8]), f[4])
        steps = []
        for d, tok in STEP.findall(path):
            if tok in self.seq:
                steps.append((tok, "+" if d == ">" else "-"))
            else:
                m = re.match(r"^(.+):(\d+)-(\d+)$", tok)
                if not m:
                    die(f"unknown segment {tok} in GAF path")
                steps.extend(self.stable(m.group(1), int(m.group(2)), int(m.group(3)), "+" if d == ">" else "-"))
        if f[4] == "-":
            steps = [(n, "-" if o == "+" else "+") for n, o in reversed(steps)]
        return steps

    def add_gaf(self, path, min_mapq=0):
        alns = defaultdict(list)
        with opener(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12 or f[5] == "*":
                    continue
                t = tags(f[12:])
                if t.get("tp", "P") != "P" or int(f[11]) < min_mapq:
                    continue
                alns[f[0]].append((int(f[2]), int(f[3]), f))
        n = 0
        for q, al in alns.items():
            al.sort()
            for qs, qe, f in al:
                name = q if len(al) == 1 else f"{q}:{qs}-{qe}"
                if name in self.paths:
                    continue
                self.paths[name] = self.gaf_steps(f)
                n += 1
        return n

    def add_rank0(self):
        n = 0
        for sn, segs in self.index.items():
            if all(self.tag[s].get("SR", "0") == "0" for _, s in segs) and sn not in self.paths:
                self.paths[sn] = [(s, "+") for _, s in segs]
                n += 1
        return n

    def write(self, path, keep=None):
        names = [p for p in self.paths if keep is None or parse_path_name(p)["sample"] in keep]
        used = {n for p in names for n, _ in self.paths[p]} if names else set(self.order)
        ids = {n: n for n in self.order}
        if not all(n.isdigit() for n in self.order):
            ids = {n: str(i) for i, n in enumerate(self.order, 1)}
        links = set(self.links)
        for p in names:
            st = self.paths[p]
            for (a, ao), (b, bo) in zip(st, st[1:]):
                if (a, ao, b, bo) not in links and (b, flip(bo), a, flip(ao)) not in links:
                    links.add((a, ao, b, bo))
        with open(path, "w") as fh:
            fh.write("H\tVN:Z:1.0\n")
            for n in self.order:
                if n in used:
                    fh.write(f"S\t{ids[n]}\t{self.seq[n]}\n")
            for a, ao, b, bo in sorted(links):
                if a in used and b in used:
                    fh.write(f"L\t{ids[a]}\t{ao}\t{ids[b]}\t{bo}\t0M\n")
            for p in names:
                fh.write(f"P\t{p}\t{','.join(ids[n] + o for n, o in self.paths[p])}\t*\n")
        return len(names)


def flip(o):
    return "-" if o == "+" else "+"


def normalise(gfa, out, gafs=(), min_mapq=0, keep=None):
    g = GFA(gfa)
    if not g.seq:
        die(f"no segments in {gfa}")
    if g.index:
        r0 = g.add_rank0()
        if r0:
            log(f"{r0} reference path(s) rebuilt from rGFA rank-0 segments")
    for gaf in gafs:
        k = g.add_gaf(gaf, min_mapq)
        log(f"{os.path.basename(gaf)}: {k} path(s) from primary alignments")
    missing = [p for p, st in g.paths.items() if any(n not in g.seq for n, _ in st)]
    for p in missing:
        log(f"skipped path {p}: it uses segments that are not in the graph")
        del g.paths[p]
    if not g.paths:
        die(f"no paths in {gfa}; minigraph graphs need --gaf alignments of the assemblies (minigraph -cxasm --vc)")
    n = g.write(out, keep)
    log(f"{len(g.seq)} segments, {n} paths written to {out}")
    return out


def convert(args):
    keep = set(args.samples) if args.samples else None
    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)
    if not args.region:
        normalise(args.gfa, args.out, args.gaf or [], args.min_mapq, keep)
        return
    odgi = args.odgi or tool("odgi")
    stem = os.path.splitext(os.path.abspath(args.out))[0]
    full = normalise(args.gfa, stem + ".full.gfa", args.gaf or [], args.min_mapq, keep)
    run([odgi, "build", "-g", full, "-o", stem + ".full.og", "-O", "-t", args.threads])
    run([odgi, "extract", "-i", stem + ".full.og", "-r", args.region, "-c", 0, "-d", args.max_distance,
         "-O", "-t", args.threads, "-o", stem + ".raw.og"])
    run([odgi, "sort", "-i", stem + ".raw.og", "-p", "Ygs", "-O", "-t", args.threads, "-o", stem + ".og"])
    with open(args.out, "w") as fh:
        run([odgi, "view", "-i", stem + ".og", "-g"], stdout=fh)
    for f in (full, stem + ".full.og", stem + ".raw.og"):
        if not args.keep_intermediate and os.path.exists(f):
            os.remove(f)
    log(f"locus graph for {args.region} written to {args.out}")
