import importlib
import os
import platform
import shutil
import subprocess
import sys

from . import __version__

PACKAGES = [("pysam", "pysam"), ("numpy", "numpy"), ("matplotlib", "matplotlib")]

TOOLS = [
    ("minimap2", "MINIMAP2", ["--version"], "regions (-L, alignment mode)"),
    ("odgi", "ODGI", ["version"], "graph (--engine native, --odgi)"),
    ("pggb", "PGGB", ["--version"], "graph (--engine native)"),
    ("docker", "DOCKER", ["--version"], "graph (--engine docker)"),
    ("singularity", "SINGULARITY", ["--version"], "graph (--engine singularity)"),
    ("apptainer", "APPTAINER", ["--version"], "graph (--engine apptainer)"),
    ("blastn", "BLASTN", ["-version"], "annotate (optional, else 21-mer search)"),
    ("makeblastdb", "MAKEBLASTDB", ["-version"], "annotate (with blastn)"),
    ("samtools", "SAMTOOLS", ["--version"], "workflow/subset_bams.sh"),
]


def first_line(cmd):
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    for line in res.stdout.splitlines():
        line = line.strip()
        if line:
            return line[:60]
    return ""


def find(name, env):
    path = os.environ.get(env)
    if path:
        return path if os.path.exists(path) or shutil.which(path) else None
    return shutil.which(name)


def check(args):
    rows, missing = [], []
    for mod, name in PACKAGES:
        try:
            m = importlib.import_module(mod)
            rows.append((name, "ok", getattr(m, "__version__", "?"), "all", os.path.dirname(m.__file__)))
        except ImportError:
            rows.append((name, "missing", "", "all", ""))
            missing.append(name)
    for name, env, vflag, used in TOOLS:
        path = find(name, env)
        if path:
            rows.append((name, "ok", first_line([path] + vflag), used, path))
        else:
            rows.append((name, "not found", "", used, ""))
            if name in (args.require or []):
                missing.append(name)
    print(f"panbreak {__version__}  python {platform.python_version()}  {platform.system()} {platform.machine()}")
    print()
    head = ("name", "status", "version", "used by", "path")
    wid = [max(len(str(r[i])) for r in rows + [head]) for i in range(4)]
    for r in [head] + rows:
        print("  ".join(str(r[i]).ljust(wid[i]) for i in range(4)) + "  " + r[4])
    engines = [n for n in ("docker", "singularity", "apptainer", "pggb") if find(n, n.upper())]
    print()
    print("graph engines available:", ", ".join("native" if e == "pggb" else e for e in engines) or "none")
    unknown = [n for n in args.require or [] if n not in {t[0] for t in TOOLS}]
    if unknown:
        print("unknown tool(s) in --require:", ", ".join(unknown), file=sys.stderr)
        raise SystemExit(2)
    if missing:
        print("missing:", ", ".join(missing), file=sys.stderr)
        raise SystemExit(1)
