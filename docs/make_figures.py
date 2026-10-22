import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, ".."))
from panbreak.common import read_tsv

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
                     "svg.fonttype": "none", "pdf.fonttype": 42})
OUT = os.path.join(here, "figures")
RES = os.path.join(here, "..", "example", "results")
GEN = ["#E15759", "#4E79A7", "#B6992D", "#59A14F", "#9C755F", "#76B7B2", "#F28E2B", "#B07AA1"]
DEL, TE, HEL, DUP, TAN, REFC = "#B2182B", "#8C510A", "#C9A227", "#2166AC", "#35978F", "#9E9E9E"
PRESENT, ABSENT = "#3E9B4F", "#C9BFB3"
MONO = {"family": "monospace", "fontsize": 5.6}


def box(ax, x, y, w, h, title=None, fc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fc, ec="#3A3A3A", lw=0.8))
    if title:
        ax.text(x + 1.2, y + h - 1.4, title, ha="left", va="top", fontsize=8, fontweight="bold", style="italic")


def arrow(ax, x0, y0, x1, y1, text=None, side="right", col="#555555"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="simple,head_length=6,head_width=7,tail_width=2.4",
                                 color=col, lw=0, mutation_scale=1))
    if text:
        ax.text((x0 + x1) / 2 + (1.2 if side == "right" else -1.2), (y0 + y1) / 2, text, ha="left" if side == "right" else "right",
                va="center", fontsize=6.5, linespacing=1.15)


def kmer_ticks(ax, x, y, bp_x, n=8, k=6.0, dy=0.55, col=PRESENT, lw=1.0):
    for i in range(n):
        s = bp_x - k + 1.0 + i * (k - 2.0) / (n - 1)
        ax.plot([s, s + k], [y - i * dy, y - i * dy], color=col, lw=lw, solid_capstyle="butt")


def run_numbers():
    st = [l.split("\t") for l in open(os.path.join(RES, "02_graph", "locus_stats.tsv")) if not l.startswith("#")][-1]
    js = read_tsv(os.path.join(RES, "04_kmers", "junction_summary.tsv"))
    kept = len(read_tsv(os.path.join(RES, "04_kmers", "junction_kmers.tsv")))
    ok = sum(r["genotypable"] == "yes" for r in js)
    fp = sum(int(r["asm_false_carriers"]) for r in js)
    return {"nodes": int(st[1]), "kept": kept, "bubbles": len(js), "ok": ok, "fp": fp}


def figure_workflow():
    N = run_numbers()
    fig = plt.figure(figsize=(180 / 25.4, 250 / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 139)
    ax.axis("off")

    box(ax, 2, 117, 40, 20, "Genome assemblies ATG1 ... ATGn")
    for i in range(8):
        col, row = i % 2, i // 2
        x, y = 4 + col * 19, 131 - row * 3.3
        ax.text(x, y + 1.0, f"ATG{i + 1}", fontsize=6, fontweight="bold", color=GEN[i])
        ax.add_patch(Rectangle((x, y - 0.35), 16, 0.7, color=GEN[i], lw=0))
        lx = x + 6.5 + (0.6 if i in (3, 6) else 0)
        ax.add_patch(Rectangle((lx, y - 0.9), 3.2, 1.8, fc="none", ec="black", lw=0.6))
    ax.text(22, 118.0, "annotation (e.g. liftoff GFF3) or minimap2 anchors", ha="center", fontsize=5.8, color="#555555")
    arrow(ax, 43, 127, 54, 127)
    ax.text(48.5, 129.5, "panbreak regions", ha="center", fontsize=6.5)
    ax.text(48.5, 124.6, "gene ± flank\nreference orientation", ha="center", va="top", fontsize=5.6, color="#555555")
    box(ax, 56, 117, 42, 20, "Locus sequences (PanSN)")
    for i in range(8):
        y = 132.3 - i * 1.9
        L = 30 + (6 if i == 7 else 0) + (3 if i == 4 else 0) - (2 if i in (1, 2) else 0)
        ax.plot([58, 58 + L * 0.72], [y, y], color=GEN[i], lw=1.6, solid_capstyle="butt")
        ax.text(59 + L * 0.72, y, f"ATG{i + 1}#1#ATG{i + 1}_x", fontsize=5.2, va="center", family="monospace")
    arrow(ax, 77, 116.5, 77, 111, "panbreak graph  (pggb -p 90 -s 5000; odgi sort, view, viz)", "left")

    box(ax, 2, 80, 96, 30, "1  Local pangenome graph")
    nodes = np.linspace(6, 58, 14)
    y0 = 92.5
    for i, x in enumerate(nodes):
        ax.add_patch(Rectangle((x - 1.1, y0 - 0.8), 2.2, 1.6, fc="#D9D9D9", ec="#555555", lw=0.5, zorder=3))
        ax.text(x, y0, str(i + 1), ha="center", va="center", fontsize=4.8, zorder=4)
        if i:
            ax.plot([nodes[i - 1] + 1.1, x - 1.1], [y0, y0], color="#555555", lw=0.6, zorder=2)
    ax.add_patch(Arc(((nodes[3] + nodes[6]) / 2, y0 - 0.8), nodes[6] - nodes[3], 7, theta1=180, theta2=360, color=DEL, lw=1.4))
    ax.text(nodes[6] + 0.5, y0 - 3.9, "skip = deletion", ha="left", fontsize=5.8, color=DEL)
    xi = (nodes[8] + nodes[9]) / 2
    ax.add_patch(Rectangle((xi - 2.2, y0 + 3.6), 4.4, 1.6, fc=TE, ec="#333333", lw=0.5, zorder=3))
    ax.plot([nodes[8] + 1.1, xi - 2.2], [y0 + 0.8, y0 + 4.4], color=TE, lw=1.0)
    ax.plot([xi + 2.2, nodes[9] - 1.1], [y0 + 4.4, y0 + 0.8], color=TE, lw=1.0)
    ax.text(xi - 2.6, y0 + 4.4, "non-reference node\n= insertion", ha="right", va="center", fontsize=5.8, color=TE)
    ax.add_patch(Arc(((nodes[10] + nodes[13]) / 2, y0 + 0.8), nodes[13] - nodes[10], 9, theta1=0, theta2=180, color=DUP, lw=1.4))
    ax.add_patch(FancyArrowPatch((nodes[10] + 0.9, y0 + 2.6), (nodes[10] + 0.2, y0 + 1.0), arrowstyle="-|>", mutation_scale=6, color=DUP, lw=1.0))
    ax.text((nodes[10] + nodes[13]) / 2 + 2, y0 + 6.0, "back edge = tandem copy", ha="center", fontsize=5.8, color=DUP)
    for i, (g, path) in enumerate((("ATG1", "reference path"), ("ATG2", "skips 5-6"), ("ATG5", "visits insertion node"), ("ATG8", "re-walks 11-14"))):
        y = 86.6 - i * 1.55
        ax.plot([6, 48], [y, y], color=GEN[int(g[-1]) - 1], lw=1.1)
        ax.text(49, y, f"{g}: {path}", fontsize=5.4, va="center", color=GEN[int(g[-1]) - 1])
    gfa = ["H  VN:Z:1.0", "S  1  ATTCAAAACAACTCAGTCAC...", "S  2  G", "S  3  C", "L  1 + 2 + 0M",
           "L  1 + 3 + 0M", "P  ATG1#1#ATG1_x  8+,10+,11+,12+,...", "P  ATG2#1#ATG2_x  8+,10+,11+,13+,...",
           "P  ATG8#1#ATG8_x  ...,20+,21+,23+,..."]
    ax.add_patch(Rectangle((70, 82), 27, 23, fc="#F7F7F7", ec="#BDBDBD", lw=0.5))
    ax.text(71, 104, f"locus.gfa (example, {N['nodes']:,} nodes)", fontsize=5.8, fontweight="bold", va="top")
    for i, l in enumerate(gfa):
        ax.text(71, 101.6 - i * 2.1, l, va="top", **MONO)
    arrow(ax, 50, 79.5, 50, 75, "panbreak bubbles + annotate", "right")

    box(ax, 2, 49, 96, 25, "2  Bubble calling and graph annotation")
    gx0, gx1, gy = 6, 56, 66
    ax.plot([gx0, gx1], [gy, gy], color="#BDBDBD", lw=3, solid_capstyle="butt")
    for s, e in ((14, 16), (20, 20.6), (23, 23.6), (26, 26.6), (40, 42)):
        ax.add_patch(Rectangle((s, gy - 0.9), e - s, 1.8, color="#252525", lw=0))
    ax.plot([16, 40], [gy, gy], color="#525252", lw=1.0)
    ax.text(28, gy + 1.6, "geneX", style="italic", fontsize=6.5, ha="center")
    ax.text(41, gy - 2.4, "TSS", fontsize=5.4, ha="center", color="#555555")
    ax.add_patch(Arc((32, gy - 0.9), 9, 6, theta1=180, theta2=360, color=DEL, lw=1.6))
    ax.add_patch(Arc((17.0, gy - 0.9), 2.4, 3, theta1=180, theta2=360, color=DEL, lw=1.4))
    th = np.linspace(-np.pi / 2, 3 * np.pi / 2, 60)
    ax.plot(36 + 1.2 * np.cos(th), gy + 2.4 + 2.2 * np.sin(th), color=TE, lw=1.5)
    ax.plot(44 + 0.9 * np.cos(th), gy + 2.0 + 1.8 * np.sin(th), color=HEL, lw=1.5)
    ax.add_patch(Arc((32, gy + 0.9), 22, 9, theta1=0, theta2=180, color=DUP, lw=1.3))
    ax.text(31, 54.0, "walk each path along the reference path: skip > min_sv = DEL, non-reference run = INS,\n"
            "both = REPL, reversed run = INV, backward jump = DUP; cluster events across genomes\n"
            "(breakpoints ± 50 bp, size ratio ≥ 0.8, 15-mer Jaccard ≥ 0.3)", ha="center", va="center", fontsize=5.8, color="#333333")
    rows = [("B005", "DEL", "2.2 kb", "intron 1", "–", DEL), ("B002", "DEL", "445 bp", "3' UTR, exon 7", "–", DEL),
            ("B007", "INS", "4.7 kb", "intron 1", "LTR/Copia, LTRs, TSD", TE),
            ("B009", "INS", "2.1 kb", "promoter (−150)", "Helitron, full length", HEL),
            ("B004", "DUP", "3.5 kb", "intron 1", "tandem copy (+B008)", DUP), ("B003", "INS", "240 bp", "intron 1", "tandem repeat 2 x 120 bp", TAN)]
    hx = [59, 64.5, 69, 74, 84.5]
    for x, h in zip(hx, ("bubble", "type", "size", "context", "sequence class")):
        ax.text(x, 70.5, h, fontsize=5.6, fontweight="bold")
    for i, r in enumerate(rows):
        y = 68.2 - i * 2.4
        for j, (x, v) in enumerate(zip(hx, r[:5])):
            ax.text(x, y, v, fontsize=5.4, color=r[5] if j < 2 else "#222222")
    ax.text(77.5, 51.0, "GFF3 context · BLASTN vs TE library · TSD/LTR/TIR · self-similarity", ha="center", fontsize=5.4, color="#555555")
    arrow(ax, 50, 48.5, 50, 44, "panbreak kmers  (k = 31, ≥ 8 bp on each side, assembly specificity filter)", "right")

    box(ax, 2, 22, 96, 21, "3  Allele-specific junction k-mers")
    y = 36.2
    ax.text(5, y + 1.4, "ref", fontsize=6, color=REFC, fontweight="bold")
    ax.plot([8, 18], [y, y], color=REFC, lw=2.2, solid_capstyle="butt")
    ax.add_patch(Rectangle((18, y - 1.1), 10, 2.2, fc="#F2D4D4", ec=DEL, lw=0.6))
    ax.text(23, y, "deleted", ha="center", va="center", fontsize=5.2, color=DEL)
    ax.plot([28, 38], [y, y], color=REFC, lw=2.2, solid_capstyle="butt")
    kmer_ticks(ax, 0, y - 1.7, 18, n=5, k=5, dy=0.5, col=REFC)
    kmer_ticks(ax, 0, y - 1.7, 28, n=5, k=5, dy=0.5, col=REFC)
    y2 = 29.4
    ax.text(5, y2 + 1.6, "alt", fontsize=6, color=DEL, fontweight="bold")
    ax.plot([8, 18], [y2, y2], color=REFC, lw=2.2, solid_capstyle="butt")
    ax.plot([18, 28], [y2, y2], color=REFC, lw=2.2, solid_capstyle="butt")
    ax.plot([18, 18], [y2 - 1.2, y2 + 1.2], color=DEL, lw=1.0)
    kmer_ticks(ax, 0, y2 - 1.7, 18, n=5, k=5, dy=0.5, col=DEL)
    ax.text(30, y2, "1 alt junction\n2 ref junctions", fontsize=5.4, va="center")
    ax.text(23, 24.4, "all 31-mers spanning a breakpoint, from every\nassembly with that allele (pan-junction windows)",
            ha="center", va="center", fontsize=5.4, color="#555555")
    mx, my = 50, 39
    ncol = 14
    rng = np.random.default_rng(3)
    car = [1, 2]
    ax.text(mx + 5.8, my + 1.6, "alt k-mers", fontsize=5.6, ha="center", color=DEL)
    ax.text(mx + 17, my + 1.6, "ref k-mers", fontsize=5.6, ha="center", color="#555555")
    for gi in range(8):
        yy = my - gi * 1.9
        ax.text(mx - 1, yy, f"ATG{gi + 1}" + (" *" if gi in car else ""), ha="right", va="center", fontsize=5.4, color=GEN[gi])
        for c in range(ncol):
            alt = c < 7
            present = (gi in car) == alt
            if alt and gi not in car and c == 5:
                present = True
            ax.scatter(mx + 1 + c * 1.6, yy, s=9, color=PRESENT if present else ABSENT, lw=0, zorder=3)
        ax.plot([mx, mx + 22.5], [yy, yy], color="#DDDDDD", lw=0.4, zorder=1)
    ax.add_patch(Rectangle((mx + 1 + 5 * 1.6 - 0.7, my - 7 * 1.9 - 0.9), 1.6, 8 * 1.9 + 0.1, fc="none", ec="black", lw=0.6, ls="--"))
    ax.text(mx + 1 + 5 * 1.6, my - 8 * 1.9 + 0.2, "dropped", ha="center", fontsize=5.2)
    ax.text(mx + 25, 39.5, "alt: present in carriers (*),\n  absent from all others\nref: absent from carriers\n"
            "optional genome-wide\n  copy-number filter", fontsize=5.5, va="top")
    ax.text(mx + 25, 30.2, f"example: {N['kept']} specific k-mers\n{N['ok']}/{N['bubbles']} bubbles genotypable\nrecall 1.0, {N['fp']} false carriers",
            fontsize=5.5, va="top", color=PRESENT)
    ax.scatter([mx + 26], [24.0], s=9, color=PRESENT)
    ax.text(mx + 27, 24.0, "present", fontsize=5.4, va="center")
    ax.scatter([mx + 34], [24.0], s=9, color=ABSENT)
    ax.text(mx + 35, 24.0, "absent", fontsize=5.4, va="center")
    arrow(ax, 50, 21.5, 50, 17.5, "panbreak genotype  →  panbreak popstats", "right")

    box(ax, 2, 1, 96, 16, "4  Short-read genotyping and population analysis")
    jx = 16
    ax.plot([5, 28], [9, 9], color=REFC, lw=2.0)
    ax.plot([jx, jx], [5.5, 13], color="black", lw=0.6, ls=":")
    reads = [(9, 12.5, DEL), (11, 11.6, DEL), (13.5, 10.7, DEL), (8, 7.4, REFC), (12, 6.5, REFC), (19, 12.5, "#DDDDDD")]
    for x, yy, c in reads:
        ax.plot([x, x + 7], [yy, yy], color=c, lw=1.6, solid_capstyle="butt")
    ax.text(5, 3.8, "subset BAM / CRAM (reads around the locus,\nsoft-clipped and unmapped mates included)", fontsize=5.4, va="center")
    ax.text(30, 12.8, "a = reads with an alt junction k-mer\nr = reads with a ref junction k-mer", fontsize=5.6, va="top")
    ax.text(30, 7.8, "L(g) = a·ln p_g + r·ln(1 − p_g)\np = {ε, J_a/(J_a+J_r), 1 − ε}\nGT, GQ, PL, AD, FT  →  VCF",
            fontsize=5.6, va="top", family="DejaVu Sans")
    pops = ["North", "Central", "South"]
    af = {"B005": [0.0, 0.68, 0.17], "B007": [0.19, 0.0, 0.0], "B003": [0.21, 0.0, 0.46]}
    cols = {"B005": DEL, "B007": TE, "B003": TAN}
    bx = 56
    for j, p in enumerate(pops):
        ax.text(bx + 4 + j * 11, 3.0, p, ha="center", fontsize=5.6)
        for k, (b, v) in enumerate(af.items()):
            h = v[j] * 9
            ax.add_patch(Rectangle((bx + 1 + j * 11 + k * 2.2, 4.2), 2.0, h, color=cols[b], lw=0))
    for k, b in enumerate(af):
        ax.text(bx + 34, 12.5 - k * 2.0, b, fontsize=5.6, color=cols[b], va="center")
    ax.text(bx + 16, 15.2, "AF, carrier frequency, Wilson CI, Fisher enrichment,\nHudson Fst, locus haplotypes ↔ assemblies",
            ha="center", va="top", fontsize=5.4, color="#555555")
    for f in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(OUT, f"fig1_workflow.{f}"), dpi=300)
    plt.close(fig)


def hap(ax, x0, x1, y, segs, lw=2.4):
    for s, e, c in segs:
        ax.plot([s, e], [y, y], color=c, lw=lw, solid_capstyle="butt")


def figure_junctions():
    fig, axs = plt.subplots(1, 3, figsize=(180 / 25.4, 62 / 25.4))
    for ax in axs:
        ax.set_xlim(0, 60)
        ax.set_ylim(0, 40)
        ax.axis("off")
    ax = axs[0]
    ax.set_title("Deletion (DEL)", fontsize=8, color=DEL, loc="left")
    hap(ax, 2, 58, 33, [(4, 22, REFC), (22, 38, "#F2D4D4"), (38, 56, REFC)], lw=3)
    ax.text(2, 33, "ref", ha="right", va="center", fontsize=6)
    for bp in (22, 38):
        kmer_ticks(ax, 0, 30.5, bp, n=6, k=9, dy=0.7, col="#7F7F7F")
    hap(ax, 2, 58, 17, [(4, 22, REFC), (22, 40, REFC)], lw=3)
    ax.plot([22, 22], [15.5, 18.5], color=DEL, lw=1.2)
    ax.text(2, 17, "alt", ha="right", va="center", fontsize=6, color=DEL)
    kmer_ticks(ax, 0, 14.5, 22, n=6, k=9, dy=0.7, col=DEL)
    ax.text(30, 4.5, "J_a = 1, J_r = 2  →  het expectation p = 1/3\nbiallelic model", ha="center", fontsize=6)
    ax = axs[1]
    ax.set_title("TE insertion (INS)", fontsize=8, color=TE, loc="left")
    hap(ax, 2, 58, 33, [(4, 30, REFC), (30, 56, REFC)], lw=3)
    ax.add_patch(Rectangle((28.6, 32.2), 1.4, 1.6, color="#F4A582", lw=0))
    kmer_ticks(ax, 0, 30.5, 30, n=6, k=9, dy=0.7, col="#7F7F7F")
    hap(ax, 2, 58, 17, [(2, 14, REFC), (14, 46, TE), (46, 58, REFC)], lw=3)
    for x in (14, 44.6):
        ax.add_patch(Rectangle((x - 1.4 if x == 14 else x, 16.2), 1.4, 1.6, color="#F4A582", lw=0))
    for x in (14.2, 41.0):
        ax.add_patch(Polygon([[x, 18.6], [x + 4.8, 18.6], [x + 4.8, 19.6]], color="#5A3505", lw=0))
    ax.text(30, 20.6, "LTR > ... < LTR", ha="center", fontsize=5.4, color=TE)
    for bp in (14, 46):
        kmer_ticks(ax, 0, 14.5, bp, n=6, k=9, dy=0.7, col=TE)
    ax.text(12.6, 21.4, "TSD", fontsize=5.2, color="#D6604D", ha="center")
    ax.text(30, 4.5, "host|TE and TE|host junctions\nJ_a = 2, J_r = 1  →  p = 2/3", ha="center", fontsize=6)
    ax = axs[2]
    ax.set_title("Tandem duplication (DUP)", fontsize=8, color=DUP, loc="left")
    hap(ax, 2, 58, 33, [(4, 20, REFC), (20, 40, "#C6DBEF"), (40, 56, REFC)], lw=3)
    for bp in (20, 40):
        kmer_ticks(ax, 0, 30.5, bp, n=6, k=9, dy=0.7, col="#7F7F7F")
    hap(ax, 2, 58, 17, [(2, 10, REFC), (10, 28, "#C6DBEF"), (28, 46, "#C6DBEF"), (46, 58, REFC)], lw=3)
    ax.plot([28, 28], [15.5, 18.5], color=DUP, lw=1.2)
    kmer_ticks(ax, 0, 14.5, 28, n=6, k=9, dy=0.7, col=DUP)
    ax.text(19, 19.3, "unit", ha="center", fontsize=5.4, color=DUP)
    ax.text(37, 19.3, "unit", ha="center", fontsize=5.4, color=DUP)
    ax.text(30, 4.5, "unit end|unit start junction (alt)\nref junctions persist → depth;  CN ≈ 1 + 2a/r\npresence model",
            ha="center", fontsize=6)
    fig.tight_layout(w_pad=0.5)
    for f in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(OUT, f"fig2_junction_kmers.{f}"), dpi=300)
    plt.close(fig)


def figure_model():
    gts = read_tsv(os.path.join(RES, "05_genotypes", "genotypes.tsv"))
    fig, axs = plt.subplots(1, 3, figsize=(180 / 25.4, 58 / 25.4))
    cmap = {"hom_ref": "#9E9E9E", "het": "#F4A582", "hom_alt": DEL, "carrier": DUP, "weak": "#E7B800", "lowqual": "#E7B800", "nodata": "black"}
    for ax, b, h, title in ((axs[0], "B005", 1 / 3, "B005  2.2-kb deletion (p_het = 1/3)"), (axs[1], "B007", 2 / 3, "B007  LTR/Copia insertion (p_het = 2/3)")):
        rs = [r for r in gts if r["bubble"] == b]
        mx = float(np.percentile([max(int(r["alt_reads"]), int(r["ref_reads"])) for r in rs], 95)) + 5
        for frac, ls in ((0.02, ":"), (h, "--"), (0.98, ":")):
            xs = np.linspace(0, mx, 50)
            ax.plot(xs, xs * frac / (1 - frac), color="#BBBBBB", lw=0.7, ls=ls)
        for c in ("nodata", "lowqual", "hom_ref", "het", "hom_alt"):
            pts = [(int(r["ref_reads"]), int(r["alt_reads"])) for r in rs if r["call"] == c]
            if pts:
                p = np.array(pts, dtype=float) + np.random.default_rng(1).uniform(-0.25, 0.25, (len(pts), 2))
                ax.scatter(p[:, 0], p[:, 1], s=6, color=cmap[c], lw=0, label=c.replace("_", " "), alpha=0.85)
        ax.set_xlim(-0.8, mx)
        ax.set_ylim(-0.8, mx)
        ax.set_xlabel("reads with ref junction k-mers (r)", fontsize=6)
        ax.set_ylabel("reads with alt junction k-mers (a)", fontsize=6)
        ax.set_title(title, fontsize=6.5, loc="left")
        ax.tick_params(labelsize=5.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axs[1].legend(frameon=False, fontsize=5.5, loc="upper right")
    ax = axs[2]
    a, r, e = 6, 11, 0.02
    ps = {"0/0": e, "0/1": 1 / 3, "1/1": 1 - e}
    ll = {g: a * np.log(p) + r * np.log(1 - p) for g, p in ps.items()}
    m = max(ll.values())
    pl = {g: -10 * (v - m) / np.log(10) for g, v in ll.items()}
    ax.bar(list(pl), list(pl.values()), color=["#9E9E9E", "#F4A582", DEL], width=0.6)
    for i, (g, v) in enumerate(pl.items()):
        ax.text(i, v + 2, f"{abs(v):.0f}", ha="center", fontsize=6)
    ax.set_ylabel("PL (phred-scaled, 0 = best)", fontsize=6)
    ax.set_title(f"example: a = {a}, r = {r}  →  GT 0/1, GQ = {int(sorted(pl.values())[1])}", fontsize=6.5, loc="left")
    ax.tick_params(labelsize=6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    for f in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(OUT, f"fig3_genotype_model.{f}"), dpi=300)
    plt.close(fig)


def figure_gfa():
    fig = plt.figure(figsize=(180 / 25.4, 88 / 25.4))
    ax = fig.add_axes([0, 0, 0.62, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.axis("off")
    nodes = {"1": (6, 40, "ACGTTAG"), "2": (24, 46, "C"), "3": (24, 34, "T"), "4": (40, 40, "GATTACA"),
             "5": (58, 54, "TG..TE..CA"), "6": (76, 40, "CCATGA")}
    wd = {"1": 12, "2": 5, "3": 5, "4": 12, "5": 14, "6": 12}
    cols = {"5": TE}
    for n, (x, y, sq) in nodes.items():
        ax.add_patch(FancyBboxPatch((x, y - 2.2), wd[n], 4.4, boxstyle="round,pad=0,rounding_size=0.8",
                                    fc=cols.get(n, "#E6E6E6"), ec="#333333", lw=0.7, zorder=4))
        ax.text(x + wd[n] / 2, y, sq, ha="center", va="center", fontsize=5.4, family="monospace", zorder=5,
                color="white" if n in cols else "black")
        ax.text(x + wd[n] / 2, y + 3.4, f"S {n}", ha="center", fontsize=5.6, fontweight="bold", zorder=5)

    def edge(a, b, rad=0.0, col="#333333", ls="-"):
        xa, ya = nodes[a][0] + wd[a], nodes[a][1]
        xb, yb = nodes[b][0], nodes[b][1]
        ax.add_patch(FancyArrowPatch((xa, ya), (xb, yb), connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>", mutation_scale=7,
                                     color=col, lw=0.9, ls=ls, zorder=3, shrinkA=0, shrinkB=0))
    for a, b in (("1", "2"), ("1", "3"), ("2", "4"), ("3", "4"), ("4", "6"), ("4", "5"), ("5", "6")):
        edge(a, b, col=TE if "5" in (a, b) else "#333333")
    edge("2", "6", rad=-0.45, col=DEL)
    xa, ya = nodes["6"][0] + wd["6"], nodes["6"][1]
    ax.add_patch(FancyArrowPatch((xa, ya - 1.5), (nodes["4"][0], ya - 1.5), connectionstyle="arc3,rad=-0.55", arrowstyle="-|>",
                                 mutation_scale=7, color=DUP, lw=1.0, zorder=3))
    walks = [("ATG1", ["1", "2", "4", "6"], GEN[0]), ("ATG2", ["1", "3", "4", "6"], GEN[1]), ("ATG3", ["1", "2", "6"], GEN[2]),
             ("ATG5", ["1", "2", "4", "5", "6"], GEN[4]), ("ATG8", ["1", "2", "4", "6", "4", "6"], GEN[7])]
    for k, (gname, w, c) in enumerate(walks):
        y = 15.5 - k * 2.6
        ax.text(2, y, gname, fontsize=5.8, color=c, va="center", fontweight="bold")
        xx = 11
        for n in w:
            ax.add_patch(Rectangle((xx, y - 0.8), 6.5, 1.6, fc=cols.get(n, c), ec="none", alpha=0.9))
            ax.text(xx + 3.25, y, n, ha="center", va="center", fontsize=5, color="white")
            xx += 7.2
        ax.text(xx + 1, y, {"ATG1": "reference path", "ATG2": "SNP bubble (node 3)", "ATG3": "deletion: edge 2 -> 6 skips node 4",
                            "ATG5": "insertion: node 5 only in ATG5", "ATG8": "duplication: back edge 6 -> 4"}[gname],
                fontsize=5.4, va="center", color="#333333")
    ax.text(2, 59, "Nodes (S), edges (L) and paths (P)", fontsize=8, fontweight="bold", style="italic")
    ax.text(37, 51.2, "deletion edge", fontsize=5.6, color=DEL, ha="center")
    ax.text(83, 26.5, "back edge", fontsize=5.6, color=DUP, ha="center")
    ax.text(65, 59, "insertion node", fontsize=5.6, color=TE, ha="center")
    ax.text(26.5, 52.5, "SNP bubble", fontsize=5.6, color="#555555", ha="center")
    bx = fig.add_axes([0.63, 0.04, 0.36, 0.92])
    bx.axis("off")
    lines = [("H", "VN:Z:1.0", "#555555"), ("S", "1  ACGTTAG", "black"), ("S", "2  C", "black"), ("S", "3  T", "black"),
             ("S", "4  GATTACA", "black"), ("S", "5  TG...CA", TE), ("S", "6  CCATGA", "black"), ("L", "1 + 2 + 0M", "black"),
             ("L", "1 + 3 + 0M", "black"), ("L", "2 + 4 + 0M", "black"), ("L", "3 + 4 + 0M", "black"), ("L", "4 + 6 + 0M", "black"),
             ("L", "4 + 5 + 0M", TE), ("L", "5 + 6 + 0M", TE), ("L", "2 + 6 + 0M", DEL), ("L", "6 + 4 + 0M", DUP),
             ("P", "ATG1#1#ATG1_x  1+,2+,4+,6+", GEN[0]), ("P", "ATG2#1#ATG2_x  1+,3+,4+,6+", GEN[1]),
             ("P", "ATG3#1#ATG3_x  1+,2+,6+", GEN[2]), ("P", "ATG5#1#ATG5_x  1+,2+,4+,5+,6+", GEN[4]),
             ("P", "ATG8#1#ATG8_x  1+,2+,4+,6+,4+,6+", GEN[7])]
    bx.add_patch(Rectangle((0, 0), 1, 1, transform=bx.transAxes, fc="#F7F7F7", ec="#BDBDBD", lw=0.5))
    bx.text(0.04, 0.965, "locus.gfa", transform=bx.transAxes, fontsize=6.5, fontweight="bold", va="top")
    for i, (t, rest, c) in enumerate(lines):
        y = 0.9 - i * 0.041
        bx.text(0.04, y, t, transform=bx.transAxes, va="top", fontweight="bold", color="#333333", **MONO)
        bx.text(0.11, y, rest, transform=bx.transAxes, va="top", color=c, **MONO)
    for f in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(OUT, f"fig0_gfa_model.{f}"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    figure_gfa()
    figure_workflow()
    figure_junctions()
    figure_model()
