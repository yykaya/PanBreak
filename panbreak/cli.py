import argparse
import os
import sys

from . import __version__


def kv(text):
    out = {}
    for item in text or []:
        k, v = item.split("=", 1)
        out[k] = int(v)
    return out


def parser():
    p = argparse.ArgumentParser(prog="panbreak", description="Local pangenome graph construction, bubble annotation and "
                                "junction k-mer genotyping of structural variants from short-read alignments.")
    p.add_argument("--version", action="version", version=f"panbreak {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True, metavar="command")

    r = sub.add_parser("regions", help="extract the locus (+/- flank) from every assembly in reference orientation")
    r.add_argument("-a", "--assemblies", required=True, help="TSV with columns genome, fasta[, gff]")
    r.add_argument("-r", "--ref", required=True, help="reference genome name (must be in the table)")
    r.add_argument("-g", "--gene", help="gene ID/Name to locate through each genome's annotation (e.g. liftoff GFF3)")
    r.add_argument("-L", "--ref-locus", help="contig:start-end of the locus in the reference; located in other genomes with minimap2")
    r.add_argument("-l", "--locus", default="locus", help="locus name used in sequence names <genome>#1#<genome>_<locus> [locus]")
    r.add_argument("-f", "--flank", type=int, default=20000, help="bp added on both sides [20000]")
    r.add_argument("--up", type=int, help="bp upstream of the gene (overrides --flank)")
    r.add_argument("--down", type=int, help="bp downstream of the gene (overrides --flank)")
    r.add_argument("--preset", default="asm20", help="minimap2 preset for alignment-based location [asm20]")
    r.add_argument("-t", "--threads", type=int, default=4)
    r.add_argument("-o", "--out", required=True)

    g = sub.add_parser("graph", help="build the local graph with pggb and sort, view and visualise it with odgi")
    g.add_argument("-i", "--fasta", required=True, help="bgzipped, indexed PanSN FASTA (locus.fa.gz from 'regions')")
    g.add_argument("-r", "--ref", required=True, help="reference genome (PanSN sample name)")
    g.add_argument("-e", "--engine", default="docker", choices=["docker", "singularity", "apptainer", "native"])
    g.add_argument("--image", default=None, help="pggb container image (digest-pinned default)")
    g.add_argument("-p", "--identity", type=int, default=90, help="pggb -p [90]")
    g.add_argument("-s", "--segment", type=int, default=5000, help="pggb -s [5000]")
    g.add_argument("-n", "--haplotypes", type=int, default=0, help="pggb -n [number of sequences]")
    g.add_argument("--pggb-args", default="", help="extra pggb arguments, quoted")
    g.add_argument("--regions", help="regions.tsv from 'regions' (needed for --trim-up/--trim-down)")
    g.add_argument("--trim-up", type=int, help="keep this many bp upstream of the gene in every path")
    g.add_argument("--trim-down", type=int, help="keep this many bp downstream of the gene in every path")
    g.add_argument("--extend", nargs="*", default=[], metavar="GENOME=BP", help="extra downstream bp for single genomes (e.g. tandem copies beyond the gene)")
    g.add_argument("--odgi", help="odgi binary (default: from the container, or PATH with --engine native)")
    g.add_argument("--layout", action="store_true", help="also compute a 2D layout and drawing")
    g.add_argument("--skip-pggb", action="store_true", help="reuse an existing pggb/ output")
    g.add_argument("--viz-width", type=int, default=1500)
    g.add_argument("--viz-height", type=int, default=500)
    g.add_argument("-t", "--threads", type=int, default=8)
    g.add_argument("-o", "--out", required=True)

    b = sub.add_parser("bubbles", help="walk every path against the reference path and cluster SV bubbles")
    b.add_argument("-g", "--gfa", required=True, help="locus GFA (graph/locus.gfa)")
    b.add_argument("-r", "--ref", required=True)
    b.add_argument("--regions", help="regions.tsv, adds reference chromosome coordinates")
    b.add_argument("-m", "--min-sv", type=int, default=50, help="minimum skipped / inserted length [50]")
    b.add_argument("--tolerance", type=int, default=50, help="breakpoint tolerance for clustering events across genomes [50]")
    b.add_argument("--min-size-ratio", type=float, default=0.8, help="minimum size ratio within a bubble [0.8]")
    b.add_argument("--min-jaccard", type=float, default=0.3, help="minimum 15-mer Jaccard of inserted sequences within a bubble [0.3]")
    b.add_argument("--cover-pad", type=int, default=200, help="a path must span the bubble +/- this to be genotyped [200]")
    b.add_argument("--window-flank", type=int, default=60, help="junction window half-width written for k-mer design [60]")
    b.add_argument("--prefix", default="B", help="bubble ID prefix [B]")
    b.add_argument("-o", "--out", required=True)

    a = sub.add_parser("annotate", help="annotate bubbles: gene context, TE family, TSD/LTR/TIR, duplication, tandem repeat")
    a.add_argument("-b", "--bubbles", required=True, help="output folder of 'bubbles'")
    a.add_argument("-r", "--ref", required=True)
    a.add_argument("--gff", help="GFF3 of the reference genome (chromosome or locus coordinates)")
    a.add_argument("--regions", help="regions.tsv to convert chromosome GFF coordinates to the locus")
    a.add_argument("--te-lib", help="TE library FASTA (names like family#Class/Superfamily)")
    a.add_argument("--te-method", default="auto", choices=["auto", "blast", "kmer"])
    a.add_argument("--blastn", help="blastn binary")
    a.add_argument("--te-min-cov", type=float, default=0.5, help="TE coverage of the allele sequence to call it a TE [0.5]")
    a.add_argument("--promoter", type=int, default=1000, help="bp upstream of a TSS counted as promoter [1000]")
    a.add_argument("--downstream", type=int, default=500, help="bp downstream of a TES counted as downstream [500]")
    a.add_argument("-t", "--threads", type=int, default=4)
    a.add_argument("-o", "--out", help="output TSV [bubbles/bubbles.annotated.tsv]")

    k = sub.add_parser("kmers", help="design allele-specific junction k-mers and filter them against all assemblies")
    k.add_argument("-b", "--bubbles", required=True)
    k.add_argument("-k", type=int, default=31, help="k-mer size [31]")
    k.add_argument("--min-side", type=int, default=8, help="minimum bases on each side of the breakpoint [8]")
    k.add_argument("--fasta", nargs="*", help="extra assembly sequences (PanSN FASTA, e.g. locus.fa.gz) for the specificity filter")
    k.add_argument("--background", nargs="*", help="whole-genome FASTAs: drop k-mers with more copies than --max-copies")
    k.add_argument("--max-copies", type=int, default=1)
    k.add_argument("-t", "--threads", type=int, default=4)
    k.add_argument("-o", "--out", help="output folder [bubbles folder]")

    t = sub.add_parser("genotype", help="count junction k-mers in reads of BAM/CRAM files and call genotypes")
    t.add_argument("-k", "--kmers", required=True, help="folder with junction_kmers.tsv and junction_summary.tsv")
    t.add_argument("-B", "--bams", required=True, nargs="+", help="BAM/CRAM files, folders, or a list file")
    t.add_argument("-b", "--bubbles", help="bubbles folder (VCF annotation, read windows) [same as --kmers]")
    t.add_argument("--bam-ref", help="FASTA the reads were aligned to; reads are then fetched only around the bubbles")
    t.add_argument("--bam-region", help="contig:start-end of the locus in --bam-ref (speeds up anchoring)")
    t.add_argument("--window", type=int, default=300, help="bp around the anchored breakpoints [300]")
    t.add_argument("--cram-ref", help="reference for CRAM decoding")
    t.add_argument("--sample-from", default="header", choices=["header", "file"], help="sample name from @RG SM or the file name")
    t.add_argument("--ploidy", type=int, default=2, choices=[1, 2])
    t.add_argument("--error", type=float, default=0.02, help="per-read allele error rate in the likelihood [0.02]")
    t.add_argument("--min-alt", type=int, default=2, help="alt reads needed for an alt genotype [2]")
    t.add_argument("--min-ref", type=int, default=2, help="ref reads needed for a hom-ref genotype [2]")
    t.add_argument("--min-depth", type=int, default=2, help="informative reads needed for any call [2]")
    t.add_argument("--min-gq", type=int, default=10, help="minimum GQ for PASS [10]")
    t.add_argument("--min-mapq", type=int, default=0, help="skip mapped reads below this MAPQ [0]")
    t.add_argument("--keep-duplicates", action="store_true")
    t.add_argument("-t", "--threads", type=int, default=4)
    t.add_argument("-o", "--out", required=True)

    mg = sub.add_parser("merge", help="merge genotype folders from chunked runs")
    mg.add_argument("inputs", nargs="+")
    mg.add_argument("-b", "--bubbles", required=True)
    mg.add_argument("--ploidy", type=int, default=2, choices=[1, 2])
    mg.add_argument("-o", "--out", required=True)

    ps = sub.add_parser("popstats", help="allele frequencies, enrichment, Fst and locus haplotypes per population")
    ps.add_argument("-G", "--genotypes", required=True, help="genotypes.tsv")
    ps.add_argument("-s", "--samples", help="TSV with a sample column and a grouping column")
    ps.add_argument("-c", "--group-col", default="population")
    ps.add_argument("-b", "--bubbles", help="bubbles folder (labels, overlaps, assembly haplotypes)")
    ps.add_argument("--ploidy", type=int, default=2, choices=[1, 2])
    ps.add_argument("--drop-unassigned", action="store_true")
    ps.add_argument("--no-overlap-rescue", action="store_true", help="do not count no-calls of samples carrying an overlapping bubble as non-carriers")
    ps.add_argument("-o", "--out", required=True)

    pl = sub.add_parser("plot", help="locus bubbles, graph layout, haplotype flow, read support and population frequencies")
    pl.add_argument("-b", "--bubbles", required=True)
    pl.add_argument("-r", "--ref", required=True)
    pl.add_argument("--gff")
    pl.add_argument("--gfa", help="locus GFA for the graph layout plot")
    pl.add_argument("--regions")
    pl.add_argument("--genotypes", help="genotype folder")
    pl.add_argument("--popstats", help="popstats folder")
    pl.add_argument("-s", "--samples")
    pl.add_argument("-c", "--group-col", default="population")
    pl.add_argument("--width", type=float, default=180, help="locus plot width in mm [180]")
    pl.add_argument("--height", type=float, default=80, help="locus plot height in mm [80]")
    pl.add_argument("--max-haplotypes", type=int, default=12, help="haplotypes shown in the flow plot; the rest are pooled [12]")
    pl.add_argument("--formats", nargs="+", default=["png", "pdf"])
    pl.add_argument("-o", "--out", required=True)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.cmd == "regions":
        from .regions import extract
        extract(args)
    elif args.cmd == "graph":
        from .graph import DEFAULT_IMAGE, build
        args.image = args.image or DEFAULT_IMAGE
        args.extend = kv(args.extend)
        build(args)
    elif args.cmd == "bubbles":
        from .bubbles import call
        call(args)
    elif args.cmd == "annotate":
        from .annotate import annotate
        annotate(args)
    elif args.cmd == "kmers":
        from .kmers import build
        build(args)
    elif args.cmd == "genotype":
        from .genotype import genotype
        genotype(args)
    elif args.cmd == "merge":
        from .genotype import merge
        merge(args)
    elif args.cmd == "popstats":
        from .popgen import stats
        stats(args)
    elif args.cmd == "plot":
        from .plot import plots
        plots(args)


if __name__ == "__main__":
    main(sys.argv[1:])
