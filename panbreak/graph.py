import glob
import os
import shlex
import shutil
from collections import defaultdict

from .common import die, log, mkdir, parse_path_name, read_fasta, read_tsv, run, tool, write_fasta
from .gfa import normalise

DEFAULT_IMAGE = "ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474"
CACTUS_IMAGE = "quay.io/comparative-genomics-toolkit/cactus@sha256:fb4f1bab95d9c3b003d1036989f977e8cde2ffcf90a27691b992bd7d2d1a253d"


def launcher(engine, image, workdir):
    if engine == "native":
        return lambda t: [tool(t)]
    if engine == "docker":
        docker = tool("docker")
        return lambda t: [docker, "run", "--rm", "--platform", "linux/amd64", "-v", f"{workdir}:/data", "-w", "/data", image, t]
    if engine in ("singularity", "apptainer"):
        exe = tool(engine)
        return lambda t: [exe, "exec", "--bind", f"{workdir}:/data", "--pwd", "/data", f"docker://{image}", t]
    die(f"unknown engine {engine}")


def split_fasta(fasta, out, ref, contig_names=False):
    d = os.path.join(out, "fa")
    shutil.rmtree(d, ignore_errors=True)
    mkdir(d)
    files = {}
    for name, seq in read_fasta(fasta):
        g = parse_path_name(name)["sample"]
        files.setdefault(g, os.path.join(d, f"{g}.fa"))
        with open(files[g], "a") as fh:
            write_fasta(fh, parse_path_name(name)["contig"] if contig_names else name, seq)
    order = [ref] + sorted(g for g in files if g != ref)
    return order, files


def run_pggb(args, out, pggb, n):
    if not args.skip_build:
        shutil.rmtree(os.path.join(out, "pggb"), ignore_errors=True)
        with open(os.path.join(out, "pggb.log"), "w") as lg:
            run(pggb("pggb") + ["-i", "input.fa.gz", "-o", "pggb", "-n", n, "-t", args.threads,
                                "-p", args.identity, "-s", args.segment, "-V", args.ref] + shlex.split(args.pggb_args),
                cwd=out, stdout=lg, stderr=lg)
    ogs = sorted(glob.glob(os.path.join(out, "pggb", "*.final.og")))
    if not ogs:
        die("pggb produced no *.final.og, see pggb.log")
    return os.path.relpath(ogs[0], out)


def run_minigraph(args, out, odgi):
    d = os.path.join(out, "minigraph")
    if not args.skip_build:
        mg = tool("minigraph")
        shutil.rmtree(d, ignore_errors=True)
        mkdir(os.path.join(d, "gaf"))
        order, files = split_fasta(args.fasta, d, args.ref)
        opts = shlex.split(args.minigraph_args)
        with open(os.path.join(d, "locus.rgfa"), "w") as fh, open(os.path.join(out, "minigraph.log"), "w") as lg:
            run([mg, "-cxggs", "-t", args.threads] + opts + [files[g] for g in order], stdout=fh, stderr=lg)
            for g in order:
                with open(os.path.join(d, "gaf", f"{g}.gaf"), "w") as gh:
                    run([mg, "-cxasm", "--vc", "-t", args.threads] + opts + [os.path.join(d, "locus.rgfa"), files[g]], stdout=gh, stderr=lg)
    if not os.path.exists(os.path.join(d, "locus.rgfa")):
        die("no minigraph/locus.rgfa, see minigraph.log")
    normalise(os.path.join(d, "locus.rgfa"), os.path.join(d, "locus.gfa"), sorted(glob.glob(os.path.join(d, "gaf", "*.gaf"))))
    run(odgi("odgi") + ["build", "-g", "minigraph/locus.gfa", "-o", "minigraph/locus.og", "-O", "-t", args.threads], cwd=out)
    return "minigraph/locus.og"


def run_cactus(args, out, odgi):
    d = os.path.join(out, "cactus")
    if not args.skip_build:
        shutil.rmtree(d, ignore_errors=True)
        mkdir(d)
        order, files = split_fasta(args.fasta, d, args.ref, contig_names=True)
        with open(os.path.join(d, "seqfile.txt"), "w") as fh:
            for g in order:
                fh.write(f"{g if g == args.ref else g + '.1'}\tcactus/fa/{g}.fa\n")
        cactus = launcher(args.engine, args.cactus_image, out)
        with open(os.path.join(out, "cactus.log"), "w") as lg:
            run(cactus("cactus-pangenome") + ["cactus/js", "cactus/seqfile.txt", "--outDir", "cactus/out", "--outName", "locus",
                                              "--reference", args.ref, "--gfa", "full", "--maxCores", args.threads] + shlex.split(args.cactus_args),
                cwd=out, stdout=lg, stderr=lg)
    full = sorted(glob.glob(os.path.join(d, "out", "locus.full.gfa*")))
    if not full:
        die("no cactus/out/locus.full.gfa.gz, see cactus.log")
    normalise(full[0], os.path.join(d, "locus.gfa"))
    run(odgi("odgi") + ["build", "-g", "cactus/locus.gfa", "-o", "cactus/locus.og", "-O", "-t", args.threads], cwd=out)
    return "cactus/locus.og"


def build(args):
    out = os.path.abspath(mkdir(args.out))
    pggb = launcher(args.engine, args.image, out) if args.builder == "pggb" or not (args.odgi or shutil.which("odgi") or os.environ.get("ODGI")) else None
    if args.odgi:
        odgi = lambda t: [args.odgi]
    elif args.engine == "native" or os.environ.get("ODGI") or (args.builder != "pggb" and shutil.which("odgi")):
        odgi = lambda t: [tool("odgi")]
    else:
        odgi = pggb
    fa = os.path.join(out, "input.fa.gz")
    if not os.path.exists(args.fasta + ".fai"):
        die(f"{args.fasta} must be bgzipped and indexed (panbreak regions does this)")
    shutil.copyfile(args.fasta, fa)
    shutil.copyfile(args.fasta + ".fai", fa + ".fai")
    if os.path.exists(args.fasta + ".gzi"):
        shutil.copyfile(args.fasta + ".gzi", fa + ".gzi")
    paths = [l.split("\t")[0] for l in open(fa + ".fai")]
    if not any(p.split("#")[0] == args.ref for p in paths):
        die(f"reference {args.ref} not among the sequences of {args.fasta}")
    n = args.haplotypes or len(paths)
    if args.builder == "minigraph":
        og = run_minigraph(args, out, odgi)
    elif args.builder == "cactus":
        og = run_cactus(args, out, odgi)
    else:
        og = run_pggb(args, out, pggb, n)
    src = og
    if args.regions and (args.trim_up is not None or args.trim_down is not None):
        up = args.trim_up if args.trim_up is not None else 0
        down = args.trim_down if args.trim_down is not None else 0
        with open(os.path.join(out, "locus_paths.txt"), "w") as fh:
            run(odgi("odgi") + ["paths", "-i", og, "-L"], cwd=out, stdout=fh)
        names = defaultdict(list)
        for p in open(os.path.join(out, "locus_paths.txt")).read().split():
            names[parse_path_name(p)["sample"]].append(p)
        with open(os.path.join(out, "locus_paths.bed"), "w") as bed:
            for r in read_tsv(args.regions):
                gs, ge, ln = int(r["gene_start"]), int(r["gene_end"]), int(r["length"])
                lo, hi = (gs - up, ge + down) if r["gene_strand"] == "+" else (gs - down, ge + up)
                extra = args.extend.get(r["genome"], 0)
                if r["gene_strand"] == "+":
                    hi += extra
                else:
                    lo -= extra
                for p in names.get(r["genome"], [r["path"]]):
                    bed.write(f"{p}\t{max(0, lo - 1)}\t{min(ln, hi)}\n")
        run(odgi("odgi") + ["extract", "-i", og, "-b", "locus_paths.bed", "-c", "0", "-O", "-t", args.threads, "-o", "locus.raw.og"], cwd=out)
        src = "locus.raw.og"
    run(odgi("odgi") + ["sort", "-i", src, "-p", "Ygs", "-O", "-t", args.threads, "-o", "locus.og"], cwd=out)
    with open(os.path.join(out, "locus.gfa"), "w") as fh:
        run(odgi("odgi") + ["view", "-i", "locus.og", "-g"], cwd=out, stdout=fh)
    with open(os.path.join(out, "locus_stats.tsv"), "w") as fh:
        run(odgi("odgi") + ["stats", "-i", "locus.og", "-S"], cwd=out, stdout=fh)
    with open(os.path.join(out, "locus_paths.txt"), "w") as fh:
        run(odgi("odgi") + ["paths", "-i", "locus.og", "-L"], cwd=out, stdout=fh)
    run(odgi("odgi") + ["viz", "-i", "locus.og", "-o", "locus_1D.png", "-x", args.viz_width, "-y", args.viz_height, "-a", 10, "-t", args.threads], cwd=out)
    run(odgi("odgi") + ["viz", "-i", "locus.og", "-o", "locus_1D_depth.png", "-x", args.viz_width, "-y", args.viz_height, "-a", 10, "-m", "-t", args.threads], cwd=out)
    run(odgi("odgi") + ["viz", "-i", "locus.og", "-o", "locus_1D_orientation.png", "-x", args.viz_width, "-y", args.viz_height, "-a", 10, "-z", "-t", args.threads], cwd=out)
    if args.layout:
        run(odgi("odgi") + ["layout", "-i", "locus.og", "-o", "locus.lay", "-T", "locus_layout.tsv", "-t", args.threads], cwd=out)
        run(odgi("odgi") + ["draw", "-i", "locus.og", "-c", "locus.lay", "-p", "locus_2D.png", "-H", 1200, "-C", "-w", 30], cwd=out)
    stats = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "locus_stats.tsv")) if not l.startswith("#")]
    if stats:
        log("locus graph (length, nodes, edges, paths):", " ".join(stats[-1]))
    log(f"GFA written to {out}/locus.gfa")
