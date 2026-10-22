import glob
import os
import shlex
import shutil

from .common import die, log, mkdir, read_tsv, run, tool

DEFAULT_IMAGE = "ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474"


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


def build(args):
    out = os.path.abspath(mkdir(args.out))
    pggb = launcher(args.engine, args.image, out)
    if args.odgi:
        odgi = lambda t: [args.odgi]
    elif args.engine == "native" or os.environ.get("ODGI"):
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
    if not args.skip_pggb:
        shutil.rmtree(os.path.join(out, "pggb"), ignore_errors=True)
        with open(os.path.join(out, "pggb.log"), "w") as lg:
            run(pggb("pggb") + ["-i", "input.fa.gz", "-o", "pggb", "-n", n, "-t", args.threads,
                                "-p", args.identity, "-s", args.segment, "-V", args.ref] + shlex.split(args.pggb_args),
                cwd=out, stdout=lg, stderr=lg)
    ogs = sorted(glob.glob(os.path.join(out, "pggb", "*.final.og")))
    if not ogs:
        die("pggb produced no *.final.og, see pggb.log")
    og = os.path.relpath(ogs[0], out)
    src = og
    if args.regions and (args.trim_up is not None or args.trim_down is not None):
        up = args.trim_up if args.trim_up is not None else 0
        down = args.trim_down if args.trim_down is not None else 0
        with open(os.path.join(out, "locus_paths.bed"), "w") as bed:
            for r in read_tsv(args.regions):
                gs, ge, ln = int(r["gene_start"]), int(r["gene_end"]), int(r["length"])
                lo, hi = (gs - up, ge + down) if r["gene_strand"] == "+" else (gs - down, ge + up)
                extra = args.extend.get(r["genome"], 0)
                if r["gene_strand"] == "+":
                    hi += extra
                else:
                    lo -= extra
                bed.write(f"{r['path']}\t{max(0, lo - 1)}\t{min(ln, hi)}\n")
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
