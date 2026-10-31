from __future__ import annotations
import gzip, hashlib, json, os, re, shutil, subprocess, sys, tempfile, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

def load_config():
    values = {}
    path = Path(os.environ.get("METAG_CONFIG", ROOT / "f" / "metag.conf"))
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        values[k.strip()] = v.strip()
    return values

CFG = load_config()

DIR_REF = Path(CFG.get("DIR_REF", "/mnt/f/annot/metag/raw"))
DIR_META = Path(CFG.get("DIR_META", "/mnt/f/annot/metag/0.metadata"))
DIR_LIB = Path(CFG.get("DIR_LIB", "/mnt/f/annot/metag/lib"))
DIR_RAW = Path(CFG.get("DIR_RAW", "/mnt/f/metag/raw"))
DIR_SPLIT = Path(CFG.get("DIR_SPLIT", "/mnt/f/metag/split"))
DIR_OUT = Path(CFG.get("DIR_OUT", "/mnt/d/analysis/metag"))
THREADS = int(CFG.get("THREADS", "16"))
PART_COMPRESSED_GIB = float(CFG.get("PART_COMPRESSED_GIB", "3"))
SOURMASH_K = int(CFG.get("SOURMASH_K", "31"))
SOURMASH_SCALED = int(CFG.get("SOURMASH_SCALED", "1000"))
SOURMASH_BUILD_WORKERS = int(CFG.get("SOURMASH_BUILD_WORKERS", "16"))

REF_GROUPS = ("viruses", "microbes", "plants", "invertebrates", "vertebrates")
FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
FASTQ_RE = re.compile(r"(.+?)(?:_R?([12])(?:_001)?|\.([12]))?\.(?:fastq|fq)(?:\.gz)?$", re.I)

def run(cmd, *, cwd=None, env=None, stdout=None, stderr=None, check=True):
    cmd = [str(x) for x in cmd]
    print("+", " ".join(shlex_quote(x) for x in cmd), flush=True)
    return subprocess.run(cmd, cwd=cwd, env=env, stdout=stdout, stderr=stderr, check=check)

def shlex_quote(s):
    import shlex
    return shlex.quote(str(s))

def require(program):
    path = shutil.which(program)
    if not path:
        raise RuntimeError(f"Required program not found: {program}")
    return path

def sha256(path: Path, block=8 << 20):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(block), b""):
            h.update(b)
    return h.hexdigest()

def atomic_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.partial")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    tmp.replace(path)

def fastq_records(path: Path):
    raw = gzip.open(path, "rb") if str(path).endswith(".gz") else path.open("rb")
    with raw as f:
        n = 0
        while True:
            h = f.readline()
            if not h:
                break
            s, p, q = f.readline(), f.readline(), f.readline()
            n += 1
            if not s or not p or not q or not h.startswith(b"@") or not p.startswith(b"+"):
                raise RuntimeError(f"Invalid/truncated FASTQ {path} record {n}")
            if len(s.rstrip(b"\r\n")) != len(q.rstrip(b"\r\n")):
                raise RuntimeError(f"Sequence/quality length mismatch {path} record {n}")
            yield (h, s, p, q)

def read_id(header: bytes) -> str:
    x = header.decode(errors="replace").strip().split()[0]
    if x.startswith("@"):
        x = x[1:]
    if x.endswith("/1") or x.endswith("/2"):
        x = x[:-2]
    return x

def paired_id(h1: bytes, h2: bytes) -> str:
    a, b = read_id(h1), read_id(h2)
    if a == b:
        return a
    # SRA style SRRxxx.N.1 / .2
    ma = re.match(r"^([SED]RR\d+\.\d+)\.1$", a)
    mb = re.match(r"^([SED]RR\d+\.\d+)\.2$", b)
    if ma and mb and ma.group(1) == mb.group(1):
        return ma.group(1)
    raise RuntimeError(f"Paired FASTQ IDs differ: {a!r} vs {b!r}")

def reference_files(*, sizes=None, states=None, incomplete=None, progress=None):
    """Scan once with directory-entry metadata, avoiding repeated WSL stat calls."""
    files = []
    for group in REF_GROUPS:
        folder = DIR_REF / group
        if not folder.is_dir():
            continue
        count = 0
        if progress:
            progress(group, count, "start")
        folders = [folder]
        while folders:
            current = folders.pop()
            with os.scandir(current) as entries:
                for entry in entries:
                    if entry.is_dir(follow_symlinks=False):
                        folders.append(Path(entry.path))
                        continue
                    name = entry.name.lower()
                    if incomplete is not None and name.endswith(".aria2") and name[:-6].endswith(FASTA_SUFFIXES):
                        incomplete.append(Path(entry.path))
                    if not name.endswith(FASTA_SUFFIXES) or not entry.is_file():
                        continue
                    path = Path(entry.path)
                    files.append(path)
                    if sizes is not None or states is not None:
                        info = entry.stat()
                        if sizes is not None:
                            sizes[path] = info.st_size
                        if states is not None:
                            states[str(path.relative_to(DIR_REF))] = [
                                info.st_size, info.st_mtime_ns, info.st_ctime_ns,
                                info.st_dev, info.st_ino,
                            ]
                    count += 1
                    if progress and count % 5000 == 0:
                        progress(group, count, "scan")
        if progress:
            progress(group, count, "done")
    return sorted(files)

def assert_reference_policy():
    if CFG.get("KEEP_REFERENCE_FASTA", "TRUE").upper() != "TRUE":
        raise RuntimeError("KEEP_REFERENCE_FASTA must remain TRUE in this clean pipeline")

