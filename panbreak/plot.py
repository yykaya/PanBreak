import os
from collections import defaultdict

import matplotlib
import numpy as np

from .annotate import load_models
from .common import log, mkdir, opener, parse_region, read_fasta, read_tsv

CLASS_COL = {"deletion": "#B2182B", "deletion_te": "#D6604D", "LTR": "#8C510A", "Helitron": "#C9A227", "DNA": "#A08C5B",
             "TE": "#BF812D", "tandem": "#35978F", "dup": "#2166AC", "micro": "#80CDC1", "inv": "#762A83", "other": "#7F7F7F"}
CALL_COL = {"hom_ref": "#E3E3E3", "nodata": "#FFFFFF", "weak": "#F6E8C3", "lowqual": "#F6E8C3", "other": "#C9C3E6"}


def setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
                         "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "svg.fonttype": "none",
                         "pdf.fonttype": 42})
    return plt


def bubble_colour(b):
    t, c = b["type"], (b.get("sv_class") or "") + " " + (b.get("alt_class") or "") + " " + (b.get("ref_class") or "")
    if t == "DUP":
        return CLASS_COL["dup"]
    if t == "INV":
        return CLASS_COL["inv"]
    if t == "DEL":
        return CLASS_COL["deletion_te"] if "TE" in c else CLASS_COL["deletion"]
    if "LTR" in c:
        return CLASS_COL["LTR"]
    if "Helitron" in c:
        return CLASS_COL["Helitron"]
    if "TE DNA" in c or "/DNA" in c or "TIR" in c:
        return CLASS_COL["DNA"]
    if "TE" in c:
        return CLASS_COL["TE"]
    if "microsatellite" in c:
        return CLASS_COL["micro"]
    if "tandem" in c or "duplication" in c:
        return CLASS_COL["tandem"]
    return CLASS_COL["other"]


def save(fig, out, name, fmts):
    for f in fmts:
        fig.savefig(os.path.join(out, f"{name}.{f}"), dpi=400, bbox_inches="tight")
    log(f"{name}: {', '.join(fmts)}")


def load_bubbles(d):
    f = os.path.join(d, "bubbles.annotated.tsv")
    return read_tsv(f if os.path.exists(f) else os.path.join(d, "bubbles.tsv"))


def locus(args, plt):
    from matplotlib.patches import Arc, FancyArrowPatch, Rectangle
    bubbles = load_bubbles(args.bubbles)
    (rname, rseq), = list(read_fasta(os.path.join(args.bubbles, "reference.fa")))
    contig, r0, r1 = parse_region(rname)
    reg = None
    if args.regions:
        reg = next((r for r in read_tsv(args.regions) if r["genome"] == args.ref), None)
    models = load_models(args.gff, reg, contig, r0 - 1, len(rseq)) if args.gff else []
    ngen = len(read_tsv(os.path.join(args.bubbles, "assembly_matrix.tsv")))
    af = {}
    if args.popstats:
        for r in read_tsv(os.path.join(args.popstats, "population_frequencies.tsv")):
            if r["population"] == "ALL" and r["AF"] not in ("NA", ""):
                af[r["bubble"]] = float(r["AF"])
    fig, ax = plt.subplots(figsize=(args.width / 25.4, args.height / 25.4))
    span = r1 - r0
    kb = lambda x: x / 1000
    ax.plot([kb(r0), kb(r1)], [0, 0], color="#BDBDBD", lw=2.2, solid_capstyle="butt", zorder=1)
    for m in models:
        ax.plot([kb(m["start"]), kb(m["end"])], [0, 0], color="#525252", lw=0.9, zorder=2)
        for s, e in m["exons"]:
            ax.add_patch(Rectangle((kb(s), -0.07), kb(e - s + 1), 0.14, color="#252525", lw=0, zorder=3))
        x0, x1 = (m["start"], m["end"]) if m["strand"] == "+" else (m["end"], m["start"])
        ax.annotate("", xy=(kb(x0 + (x1 - x0) * 0.12), 0.2), xytext=(kb(x0), 0.2),
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#252525", mutation_scale=5), zorder=3)
        lx = m["start"] - 0.006 * span if m["strand"] == "+" else m["end"] + 0.006 * span
        ax.text(kb(lx), 0.08, m["id"], ha="right" if m["strand"] == "+" else "left", va="bottom", style="italic", fontsize=6.5)
    labels = []
    for b in bubbles:
        col = bubble_colour(b)
        s, e, t = int(b["pos"]), int(b["end"]), b["type"]
        n = int(b["n_carriers"])
        lw = 0.6 + 2.2 * np.sqrt(n / max(ngen, 1))
        if t in ("DEL", "REPL", "INV"):
            h = 0.35 + 0.6 * min(1, (e - s) / (0.15 * span))
            ax.add_patch(Arc((kb((s + e) / 2), 0), kb(e - s), 2 * h, theta1=180, theta2=360, color=col, lw=lw,
                             ls="--" if t == "INV" else "-", zorder=4))
            labels.append((kb((s + e) / 2), -h, -1, b, col))
        if t in ("INS", "REPL"):
            size = int(b["alt_len"])
            r = 0.004 * span * (1 + np.log10(max(size, 10)) / 2)
            hgt = 0.25 + 0.22 * np.log10(max(size, 10))
            th = np.linspace(-np.pi / 2, 3 * np.pi / 2, 100)
            ax.plot(kb(s + r * np.cos(th)), hgt / 2 + hgt / 2 * np.sin(th) + 0.02, color=col, lw=lw, zorder=4)
            labels.append((kb(s), hgt + 0.02, 1, b, col))
        if t == "DUP":
            h = 1.05
            ax.add_patch(Arc((kb((s + e) / 2), 0.05), kb(e - s), 2 * h, theta1=0, theta2=180, color=col, lw=lw, zorder=4))
            ax.add_patch(FancyArrowPatch((kb(s + 0.03 * (e - s)), 0.3), (kb(s + 0.002 * (e - s)), 0.08), arrowstyle="-|>",
                                         mutation_scale=6, color=col, lw=0.8, zorder=5))
            labels.append((kb((s + e) / 2), h + 0.05, 1, b, col))
    for side in (1, -1):
        ls = sorted([l for l in labels if l[2] == side], key=lambda l: l[0])
        for i, (x, y, _, b, col) in enumerate(ls):
            lx = kb(r0) + kb(span) * (i + 0.5) / len(ls)
            ly = side * (1.55 if i % 2 == 0 else 2.05)
            txt = f"{b['bubble']}  {b.get('label', b['type'])}\n{b['n_carriers']}/{ngen} assemblies"
            if b["bubble"] in af:
                txt += f", AF {af[b['bubble']]:.2f}"
            ax.plot([x, lx], [y, ly - side * 0.17], color="#9E9E9E", lw=0.3, zorder=2)
            ax.text(lx, ly, txt, ha="center", va="bottom" if side > 0 else "top", color=col, fontsize=5.6, linespacing=1.05)
    ax.set_xlim(kb(r0) - kb(span) * 0.01, kb(r1) + kb(span) * 0.01)
    ax.set_ylim(-2.6, 2.6)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.set_yticks([])
    ax.set_xlabel(f"{contig} (kb, reference path {args.ref})")
    save(fig, args.out, "locus_bubbles", args.formats)
    plt.close(fig)


POP_PAL = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"]
STATE_ORDER = ["alt", "het", "other", "ref", "nc"]
STATE_LAB = {"alt": "alternative allele", "het": "heterozygous", "other": "overlapping allele", "ref": "reference allele", "nc": "no call"}


def sample_states(gts, binfo):
    by = defaultdict(dict)
    carried = defaultdict(set)
    for r in gts:
        by[r["sample"]][r["bubble"]] = r
        if r["FT"] == "PASS" and r["carrier"] == "1":
            carried[r["sample"]].add(r["bubble"])
    states = {}
    for s, d in by.items():
        st = {}
        for b, r in d.items():
            if r["FT"] == "PASS":
                st[b] = "het" if r["call"] == "het" else ("alt" if r["carrier"] == "1" else "ref")
            elif set(binfo[b].get("overlaps", ".").split(",")) & carried[s]:
                st[b] = "other"
            else:
                st[b] = "nc"
        states[s] = st
    return states


def ribbon(x0, x1, y0, y1, h, n=24):
    t = np.linspace(0, 1, n)
    s = t * t * (3 - 2 * t)
    xs = x0 + (x1 - x0) * t
    top = y0 + (y1 - y0) * s
    return np.column_stack([np.concatenate([xs, xs[::-1]]), np.concatenate([top, (top + h)[::-1]])])


def flow(args, plt):
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Patch, Rectangle
    import matplotlib.colors as mc
    gts = read_tsv(os.path.join(args.genotypes, "genotypes.tsv"))
    binfo = {b["bubble"]: b for b in load_bubbles(args.bubbles)}
    bubbles = sorted({r["bubble"] for r in gts}, key=lambda b: (int(binfo[b]["pos"]), b))
    meta = {r["sample"]: r for r in read_tsv(args.samples)} if args.samples else {}
    states = sample_states(gts, binfo)
    samples = sorted(states)
    pop = {s: meta.get(s, {}).get(args.group_col, "all") for s in samples}
    pops = sorted(set(pop.values()))
    pal = POP_PAL if len(pops) <= len(POP_PAL) else [mc.to_hex(matplotlib.colormaps["tab20"](i % 20)) for i in range(len(pops))]
    pcol = {p: pal[i % len(pal)] for i, p in enumerate(pops)}
    asm_hap = defaultdict(list)
    mfile = os.path.join(args.bubbles, "assembly_matrix.tsv")
    if os.path.exists(mfile):
        for r in read_tsv(mfile):
            if all(r.get(b) in ("0", "1", "x") for b in bubbles):
                asm_hap["+".join(b for b in bubbles if r[b] == "1") or "REF"].append(r["genome"])
    hap = {s: "+".join(b for b in bubbles if states[s].get(b) in ("alt", "het")) or "REF" for s in samples}
    counts = defaultdict(int)
    for s in samples:
        counts[hap[s]] += 1
    top = sorted(counts, key=lambda h: -counts[h])[:args.max_haplotypes]
    for s in samples:
        if hap[s] not in top or counts[hap[s]] < 2:
            hap[s] = "other haplotypes"
    hcount = defaultdict(int)
    for s in samples:
        hcount[hap[s]] += 1
    haps = sorted(hcount, key=lambda h: (h == "other haplotypes", -hcount[h], h))
    match = {h: asm_hap.get(h, []) for h in haps}
    rank = {h: i for i, h in enumerate(haps)}
    key = {s: (rank[hap[s]], pops.index(pop[s]), s) for s in samples}
    cols = [("pop", None)] + [("bub", b) for b in bubbles] + [("hap", None)]
    unit = 1.0
    pad = max(1.5, 0.25 * len(samples) / max(8, len(pops), min(len(haps), 20)))
    w = 0.16
    ypos, nodes = [], []
    for ci, (kind, b) in enumerate(cols):
        if kind == "pop":
            groups = [(p, [s for s in samples if pop[s] == p]) for p in pops]
        elif kind == "hap":
            groups = [(h, [s for s in samples if hap[s] == h]) for h in haps]
        else:
            groups = [(st, [s for s in samples if states[s].get(b, "nc") == st]) for st in STATE_ORDER]
        y, pos = 0.0, {}
        for g, ss in groups:
            if not ss:
                continue
            ss.sort(key=lambda s: key[s])
            nodes.append((ci, kind, b, g, y, len(ss) * unit))
            for k, s in enumerate(ss):
                pos[s] = y + k * unit
            y += len(ss) * unit + pad
        ypos.append(pos)
    height = max(n[4] + n[5] for n in nodes)
    polys, colors = [], []
    for ci in range(len(cols) - 1):
        for s in samples:
            polys.append(ribbon(ci + w / 2, ci + 1 - w / 2, ypos[ci][s], ypos[ci + 1][s], unit))
            colors.append(mc.to_rgba(pcol[pop[s]], 0.38))
    fig, ax = plt.subplots(figsize=(args.width / 25.4, max(70, min(150, 40 + 0.4 * len(haps) * 6)) / 25.4))
    ax.add_collection(PolyCollection(polys, facecolors=colors, edgecolors="none", zorder=1))
    last_lab = [-1e9]
    for ci, kind, b, g, y, h in nodes:
        if kind == "pop":
            fc, ec = pcol[g], "none"
        elif kind == "hap":
            fc, ec = "#4D4D4D", "none"
        elif g == "alt":
            fc, ec = bubble_colour(binfo[b]), "none"
        elif g == "het":
            fc, ec = mc.to_hex(0.45 * np.array(mc.to_rgb(bubble_colour(binfo[b]))) + 0.55), "none"
        elif g == "other":
            fc, ec = CALL_COL["other"], "none"
        elif g == "ref":
            fc, ec = "#BDBDBD", "none"
        else:
            fc, ec = "white", "#9E9E9E"
        ax.add_patch(Rectangle((ci - w / 2, y), w, h, fc=fc, ec=ec, lw=0.4, zorder=3))
        if kind == "pop":
            if h >= 0.012 * height:
                ax.text(ci - w / 2 - 0.06, y + h / 2, f"{g} ({int(h / unit)})", ha="right", va="center",
                        fontsize=6.0 if len(pops) <= 8 else 4.8, color=pcol[g])
        elif kind == "hap":
            m = match.get(g, [])
            lab = g if g == "other haplotypes" else (f"{m[0]}-like" + (f" +{len(m) - 1}" if len(m) > 1 else "") if m else "no assembly match")
            detail = "" if g in ("other haplotypes", "REF") else g
            if h >= 0.008 * height and y + h / 2 - last_lab[0] >= 0.022 * height:
                last_lab[0] = y + h / 2
                ax.text(ci + w / 2 + 0.06, y + h / 2, f"{lab} ({int(h / unit)})" + (f"  {detail}" if detail else ""),
                        ha="left", va="center", fontsize=5.4 if len(haps) <= 12 else 4.6, color="#252525")
        elif g in ("alt", "het") and h >= 3 * unit:
            ax.text(ci, y + h / 2, str(int(h / unit)), ha="center", va="center", fontsize=4.8, color="white", rotation=90, zorder=4)
    for ci, (kind, b) in enumerate(cols):
        if kind == "bub":
            bi = binfo[b]
            size = int(bi["alt_len"]) if bi["type"] == "INS" else int(bi["ref_len"])
            sz = f"{size / 1000:.1f} kb" if size >= 1000 else f"{size} bp"
            ax.text(ci, -pad * 0.8, f"{b}\n{bi['type']}\n{sz}", ha="center", va="bottom", fontsize=5.6, color=bubble_colour(bi), linespacing=1.0)
        else:
            ax.text(ci + (w / 2 if kind == "pop" else 0.1), -pad * 0.8, "population" if kind == "pop" else "locus haplotype\n(assembly match)",
                    ha="right" if kind == "pop" else "left", va="bottom", fontsize=5.8, fontweight="bold")
    ax.set_xlim(-1.0, len(cols) + 1.4)
    ax.set_ylim(height + pad, -pad * 5.5)
    ax.axis("off")
    handles = [Patch(facecolor=pcol[p], label=p.replace("_", " "), alpha=0.7) for p in pops]
    handles += [Patch(facecolor="#7F7F7F", label="alternative (class colour)"), Patch(facecolor="#C8C8C8", label="heterozygous (light)"),
                Patch(facecolor=CALL_COL["other"], label="overlapping allele"), Patch(facecolor="#BDBDBD", label="reference"),
                Patch(facecolor="white", edgecolor="#9E9E9E", label="no call")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=min(8, len(handles)), frameon=False, fontsize=5.6 if len(pops) <= 8 else 5.0,
              handlelength=1.2, columnspacing=1.0)
    save(fig, args.out, "haplotype_flow", args.formats)
    plt.close(fig)


def graph_layout(args, plt):
    from bisect import bisect_right
    from matplotlib.collections import LineCollection
    import matplotlib.colors as mc
    from .bubbles import Graph, Reference
    g = Graph(args.gfa)
    ref = Reference(g, args.ref)
    ids = list(g.seq)
    if all(i.isdigit() for i in ids):
        ids.sort(key=int)
    x0, x = {}, 0
    for n in ids:
        x0[n] = x
        x += len(g.seq[n])
    total = x
    L = lambda n: len(g.seq[n])
    depth = defaultdict(int)
    users = defaultdict(set)
    for p, steps in g.paths.items():
        gen = g.info[p]["sample"]
        for (n, o), (m, q) in zip(steps, steps[1:]):
            users[frozenset((n, m))].add(gen)
        for n, o in steps:
            depth[n] += 1
    starts, xs = [], []
    pos = ref.offset
    for n, o in g.paths[ref.path]:
        starts.append(pos)
        xs.append((x0[n], o, L(n)))
        pos += L(n)

    def gx(p):
        i = max(0, bisect_right(starts, p) - 1)
        nx, o, ln = xs[i]
        d = min(max(p - starts[i], 0), ln)
        return nx + (d if o == "+" else ln - d)
    bubbles = load_bubbles(args.bubbles)
    binfo = []
    for b in bubbles:
        a, e = gx(int(b["pos"])), gx(int(b["end"]))
        lo, hi = min(a, e), max(a, e)
        if b["type"] in ("INS", "REPL"):
            hi = max(hi, lo + int(b["alt_len"]))
        binfo.append((lo - 200, hi + 200, set(b["carriers"].split(",")), b))
    edges = []
    with opener(args.gfa) as fh:
        for line in fh:
            if line.startswith("L\t"):
                f = line.split("\t")
                a, ao, c, co = f[1], f[2], f[3], f[4]
                xa = x0[a] + L(a) if ao == "+" else x0[a]
                xc = x0[c] if co == "+" else x0[c] + L(c)
                if abs(xc - xa) > 1:
                    edges.append((xa, xc, users.get(frozenset((a, c)), set())))
    names = sorted(g.paths, key=lambda p: (g.info[p]["sample"] != args.ref, p))
    top_mm, row_mm = 66, (2.6 if len(names) <= 12 else 1.5)
    path_mm = max(18, min(170, len(names) * row_mm))
    fig = plt.figure(figsize=(args.width / 25.4, (top_mm + path_mm + 12) / 25.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[top_mm, path_mm], hspace=0.04, left=0.09, right=0.98, top=0.95, bottom=0.1 * 80 / (top_mm + path_mm + 12))
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    kb = 1000.0
    th = np.linspace(0, np.pi, 50)
    segs, cols, lws, apex = [], [], [], {}
    lmax = np.log10(max(total, 10))
    for xa, xc, us in sorted(edges, key=lambda e: -abs(e[1] - e[0])):
        span = abs(xc - xa)
        mid, r = (xa + xc) / 2, span / 2
        fwd = xc > xa
        hit = None
        if span >= 30 and us and args.ref not in us:
            near = lambda lo, hi: max(0, lo - max(xa, xc), min(xa, xc) - hi)
            cand = [(len(cs), near(lo, hi), b) for lo, hi, cs, b in binfo if us <= cs and near(lo, hi) <= 2000]
            hit = min(cand, key=lambda t: (t[0], t[1]))[2] if cand else None
        h = (0.35 + 1.9 * np.log10(span + 1) / lmax) if fwd else (0.25 + 0.75 * np.log10(span + 1) / lmax)
        base, sgn = (0.32, 1) if fwd else (-0.32, -1)
        segs.append(np.column_stack([(mid - r * np.cos(th)) / kb, base + sgn * h * np.sin(th)]))
        if hit is not None:
            cols.append(bubble_colour(hit))
            lws.append(1.2)
            if fwd and (hit["bubble"] not in apex or apex[hit["bubble"]][1] < base + h):
                apex[hit["bubble"]] = (mid / kb, base + h, hit)
        else:
            cols.append("#B0B0B0")
            lws.append(0.35 if span < 30 else 0.6)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws, zorder=2))
    cmap = matplotlib.colormaps["viridis"]
    npath = len(g.paths)
    ax.add_collection(LineCollection([[(x0[n] / kb, 0), ((x0[n] + L(n)) / kb, 0)] for n in ids],
                                     colors=[cmap(min(1.0, depth[n] / max(npath, 1))) for n in ids], linewidths=10, zorder=3))
    visitors = defaultdict(set)
    for p, steps in g.paths.items():
        for n, o in steps:
            visitors[n].add(g.info[p]["sample"])
    alt_nodes = defaultdict(list)
    for n in ids:
        vs = visitors[n]
        if not vs or args.ref in vs:
            continue
        cand = [(len(cs), b) for lo, hi, cs, b in binfo if vs <= cs and b["type"] in ("INS", "REPL", "INV", "DUP")]
        if cand and L(n) >= 5:
            alt_nodes[min(cand, key=lambda t: t[0])[1]["bubble"]].append(n)
    for b in bubbles:
        ns = alt_nodes.get(b["bubble"])
        if not ns:
            continue
        col = bubble_colour(b)
        ax.add_collection(LineCollection([[(x0[n] / kb, 0.2), ((x0[n] + L(n)) / kb, 0.2)] for n in ns], colors=col, linewidths=3.2, zorder=4))
        xm = (min(x0[n] for n in ns) + max(x0[n] + L(n) for n in ns)) / 2 / kb
        if b["bubble"] not in apex:
            apex[b["bubble"]] = (xm, 0.32, b)
    placed = []
    for xm, top, b in sorted(apex.values(), key=lambda t: t[0]):
        y = top + 0.12
        while any(abs(xm - px) < total / kb * 0.045 and abs(y - py) < 0.22 for px, py in placed):
            y += 0.24
        placed.append((xm, y))
        ax.text(xm, y, b["bubble"], ha="center", va="bottom", fontsize=5.8, color=bubble_colour(b), fontweight="bold")
    models = []
    if args.gff:
        reg = next((r for r in read_tsv(args.regions) if r["genome"] == args.ref), None) if args.regions else None
        models = load_models(args.gff, reg, ref.contig, ref.offset, len(ref.seq))
    gy = -1.55
    for i, m in enumerate(sorted(models, key=lambda m: gx(m["start"]))):
        a, e = gx(m["start"]), gx(m["end"])
        ax.plot([min(a, e) / kb, max(a, e) / kb], [gy, gy], color="#525252", lw=0.7, zorder=3)
        for s0, t0 in m["exons"]:
            ea, eb = gx(s0), gx(t0)
            ax.plot([min(ea, eb) / kb, max(ea, eb) / kb], [gy, gy], color="#252525", lw=3.2, solid_capstyle="butt", zorder=4)
        ax.text(min(a, e) / kb, gy - 0.2 - (0.32 if i % 2 else 0), m["id"], fontsize=5.4 if len(models) <= 6 else 4.6, style="italic", va="top")
    rowfont = 5.6 if len(names) <= 20 else max(3.0, 5.6 * 20 / len(names) ** 0.95)
    for k, p in enumerate(names):
        visits = defaultdict(int)
        for n, o in g.paths[p]:
            visits[n] += 1
        ax2.add_collection(LineCollection([[(x0[n] / kb, -k), ((x0[n] + L(n)) / kb, -k)] for n in visits],
                                          colors=["#111111" if visits[n] > 1 else POP_PAL[k % len(POP_PAL)] for n in visits],
                                          linewidths=max(1.0, min(3.0, 110 / max(len(names), 1))), zorder=3))
        ax2.text(-0.01 * total / kb, -k, g.info[p]["sample"], ha="right", va="center", fontsize=rowfont)
    lx = -0.01 * total / kb
    for y, t, bold in ((1.2, "edges\nforward", False), (0, "nodes", True), (-0.75, "edges\nbackward", False), (gy, "genes", False)):
        ax.text(lx, y, t, ha="right", va="center", fontsize=5.8 if bold else 5.4, fontweight="bold" if bold else "normal")
    ax2.text(-0.075 * total / kb, 1.2, "paths", ha="left", va="bottom", fontsize=5.8, fontweight="bold")
    ax.set_xlim(-0.075 * total / kb, total / kb * 1.01)
    ax.set_ylim(gy - 1.0, max([p[1] for p in placed] + [2.6]) + 0.35)
    ax2.set_ylim(-len(names) + 0.4, 1.6)
    for a_ in (ax, ax2):
        a_.set_yticks([])
        for sp in ("left", "right", "top"):
            a_.spines[sp].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax2.set_xlabel("sorted graph coordinate (kb; odgi sort -p Ygs)")
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=mc.Normalize(0, npath))
    cax = fig.add_axes([0.80, 0.965, 0.15, 0.008])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("paths per node", fontsize=5.5, labelpad=1)
    cb.ax.tick_params(labelsize=5, length=1.5)
    cb.outline.set_linewidth(0.3)
    ax.text(0.0, 1.04, "coloured bar on nodes: sequence private to a bubble's carriers; black path segments: node visited more than once",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.2, color="#555555")
    save(fig, args.out, "graph_layout", args.formats)
    plt.close(fig)


def af_dots(args, plt):
    rows = read_tsv(os.path.join(args.popstats, "population_frequencies.tsv"))
    binfo = {b["bubble"]: b for b in load_bubbles(args.bubbles)}
    bubbles = sorted({r["bubble"] for r in rows}, key=lambda b: int(binfo[b]["pos"]))
    pops = ["ALL"] + sorted({r["population"] for r in rows} - {"ALL"})
    d = {(r["bubble"], r["population"]): r for r in rows}
    fig, ax = plt.subplots(figsize=((40 + 14 * len(pops)) / 25.4, (20 + 7 * len(bubbles)) / 25.4))
    for i, b in enumerate(bubbles):
        for j, p in enumerate(pops):
            r = d.get((b, p))
            if not r or r["AF"] in ("NA", ""):
                continue
            a = float(r["AF"])
            ax.scatter(j, i, s=8 + 260 * a, color=bubble_colour(binfo[b]), alpha=0.9 if a > 0 else 0.25, lw=0, zorder=3)
            ax.text(j + 0.32, i, f"{a:.2f}", va="center", fontsize=5.5, color="#333333")
    ax.set_yticks(range(len(bubbles)))
    ax.set_yticklabels([f"{b}  {binfo[b].get('label', binfo[b]['type'])}" for b in bubbles], fontsize=6)
    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels([f"{p}\n(n={d.get((bubbles[0], p), {}).get('called', '')})" for p in pops], fontsize=6.5)
    ax.set_xlim(-0.6, len(pops) - 0.2)
    ax.set_ylim(len(bubbles) - 0.4, -0.6)
    ax.grid(axis="y", color="#EEEEEE", lw=0.5, zorder=0)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    ax.set_title("Allele frequency by population", fontsize=7, loc="left")
    save(fig, args.out, "population_af", args.formats)
    plt.close(fig)


def support(args, plt):
    gts = read_tsv(os.path.join(args.genotypes, "genotypes.tsv"))
    binfo = {b["bubble"]: b for b in load_bubbles(args.bubbles)}
    bubbles = sorted({r["bubble"] for r in gts}, key=lambda b: int(binfo[b]["pos"]))
    ncol = min(5, len(bubbles))
    nrow = (len(bubbles) + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * 34 / 25.4, nrow * 32 / 25.4), squeeze=False)
    cmap = {"hom_ref": "#9E9E9E", "het": "#F4A582", "hom_alt": "#B2182B", "carrier": "#2166AC", "weak": "#E7B800",
            "nodata": "#000000", "lowqual": "#E7B800"}
    for ax, b in zip(axs.flat, bubbles):
        rs = [r for r in gts if r["bubble"] == b]
        for c in ("hom_ref", "nodata", "weak", "lowqual", "het", "hom_alt", "carrier"):
            xs = [np.sqrt(int(r["ref_reads"])) for r in rs if r["call"] == c]
            ys = [np.sqrt(int(r["alt_reads"])) for r in rs if r["call"] == c]
            if xs:
                ax.scatter(np.array(xs) + np.random.uniform(-0.08, 0.08, len(xs)), np.array(ys) + np.random.uniform(-0.08, 0.08, len(ys)),
                           s=5, color=cmap[c], lw=0, alpha=0.8, label=c.replace("_", " "))
        ax.set_title(f"{b} {binfo[b]['type']}", fontsize=6.5, color=bubble_colour(binfo[b]))
        ax.set_xlabel("√ ref reads", fontsize=6)
        ax.set_ylabel("√ alt reads", fontsize=6)
        ax.tick_params(labelsize=5.5)
        for sp in ("right", "top"):
            ax.spines[sp].set_visible(False)
    for ax in list(axs.flat)[len(bubbles):]:
        ax.axis("off")
    hs, ls = [], []
    for ax in axs.flat:
        h, l = ax.get_legend_handles_labels()
        for hh, ll in zip(h, l):
            if ll not in ls:
                hs.append(hh)
                ls.append(ll)
    fig.legend(hs, ls, loc="lower center", ncol=len(ls), frameon=False, fontsize=6, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout()
    save(fig, args.out, "read_support", args.formats)
    plt.close(fig)


def plots(args):
    plt = setup()
    mkdir(args.out)
    np.random.seed(1)
    locus(args, plt)
    if args.gfa:
        graph_layout(args, plt)
    if args.genotypes:
        flow(args, plt)
        support(args, plt)
    if args.popstats:
        af_dots(args, plt)
