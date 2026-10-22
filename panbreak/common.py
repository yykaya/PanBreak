import csv
import gzip
import os
import re
import shutil
import subprocess
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
PANSN = re.compile(r"^(?P<sample>[^#]+)#(?P<hap>[^#]+)#(?P<contig>.+?)(?::(?P<start>\d+)-(?P<end>\d+))?$")


def rc(seq):
    return seq.translate(COMP)[::-1]


def canonical(kmer):
    r = rc(kmer)
    return kmer if kmer <= r else r


def opener(path, mode="rt"):
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path):
    name, buf = None, []
    with opener(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf).upper()


def write_fasta(fh, name, seq, width=80):
    fh.write(">" + name + "\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def read_tsv(path):
    with opener(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t"))


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with opener(path, "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k, "")) for k in fields})


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        if v != v:
            return "NA"
        return f"{v:.4g}"
    if isinstance(v, (list, tuple, set)):
        return ",".join(map(str, v)) if v else "."
    return v


def parse_path_name(name):
    m = PANSN.match(name)
    if not m:
        base, _, rng = name.partition(":")
        start = int(rng.split("-")[0]) if rng and "-" in rng else 0
        return {"sample": base, "hap": "1", "contig": base, "offset": start}
    return {"sample": m.group("sample"), "hap": m.group("hap"), "contig": m.group("contig"),
            "offset": int(m.group("start")) if m.group("start") else 0}


def parse_region(text):
    m = re.match(r"^(.+):([\d,]+)-([\d,]+)$", text)
    if not m:
        raise SystemExit(f"cannot parse region {text!r}, expected contig:start-end")
    return m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))


def log(*args):
    print("[panbreak]", *args, file=sys.stderr, flush=True)


def die(msg):
    log("error:", msg)
    raise SystemExit(1)


def tool(name, env=None):
    path = os.environ.get(env or name.upper().replace("-", "_"))
    if path:
        return path
    found = shutil.which(name)
    if not found:
        die(f"{name} not found in PATH (or set ${env or name.upper()})")
    return found


def run(cmd, cwd=None, stdout=None, stderr=None):
    log("$", " ".join(map(str, cmd)))
    res = subprocess.run(list(map(str, cmd)), cwd=cwd, stdout=stdout, stderr=stderr)
    if res.returncode != 0:
        die(f"command failed ({res.returncode}): {' '.join(map(str, cmd))}")
    return res


def mkdir(path):
    os.makedirs(path, exist_ok=True)
    return path


def bgzip_index(src, dst):
    import pysam
    pysam.tabix_compress(src, dst, force=True)
    pysam.faidx(dst)
    return dst
