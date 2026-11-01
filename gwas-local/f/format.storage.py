#!/usr/bin/env python3
"""Stage GWAS references and build reusable lead coordinate indexes."""

import argparse
from contextlib import closing, contextmanager
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile


@contextmanager
def locked(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def identity(path):
    info = path.stat()
    return str(path.resolve()), info.st_size, info.st_mtime_ns


def check_space(root, size):
    free = shutil.disk_usage(root).free
    if free < size + 1024**3:
        raise RuntimeError(f"Insufficient local cache space: need {size + 1024**3:,} bytes, have {free:,}: {root}")


def stage(root, sources):
    # Keep the requested basenames: a MAGMA prefix may consist of symlinks to
    # differently named source files, but its local .bed/.bim/.fam must match.
    requested = [Path(p).absolute() for p in sources]
    sources = list(dict.fromkeys(requested))
    if len({p.name for p in sources}) != len(sources):
        raise ValueError("Reference group contains duplicate file names")
    before = [identity(p) for p in sources]
    cache_identity = before
    if any(p.name != Path(info[0]).name for p, info in zip(sources, before)):
        cache_identity = ([p.name for p in sources], before)
    key = hashlib.sha256(json.dumps(cache_identity).encode()).hexdigest()[:24]
    root.mkdir(parents=True, exist_ok=True)
    target = root / key
    with locked(root / (key + ".lock")):
        ready = target / ".complete"
        try:
            valid = json.loads(ready.read_text()) == [list(info) for info in before] and all(
                (target / p.name).is_file() and (target / p.name).stat().st_size == info[1]
                for p, info in zip(sources, before)
            )
        except (OSError, ValueError):
            valid = False
        if not valid:
            check_space(root, sum(info[1] for info in before))
            work = Path(tempfile.mkdtemp(prefix=key + ".tmp.", dir=root))
            try:
                for source in sources:
                    print(f"Cache reference: {source}", file=sys.stderr, flush=True)
                    shutil.copy2(source, work / source.name)
                if [identity(p) for p in sources] != before:
                    raise RuntimeError("Reference changed while being copied; retry after its writer finishes")
                (work / ".complete").write_text(json.dumps(before) + "\n")
                if target.exists():
                    # Only a damaged/incomplete cache generation, never source data.
                    shutil.rmtree(target)
                work.rename(target)
            finally:
                if work.exists():
                    shutil.rmtree(work)
    for source in requested:
        print(target / source.name)


@contextmanager
def reference_text(source):
    """Stream compressed references without loading a chromosome into memory."""
    if source.suffix == ".zst":
        with subprocess.Popen(["zstd", "-q", "-dc", "--", str(source)],
                              stdout=subprocess.PIPE, text=True) as proc:
            try:
                yield proc.stdout
            finally:
                proc.stdout.close()
            if proc.wait():
                raise RuntimeError(f"Could not decompress reference: {source}")
    elif source.suffix in (".gz", ".bgz"):
        with gzip.open(source, "rt") as handle:
            yield handle
    else:
        with source.open() as handle:
            yield handle


def norm_chrom(value):
    value = value.upper()
    if value.startswith("CHR"):
        value = value[3:]
    value = {"X": "23", "Y": "24", "M": "25", "MT": "25"}.get(value, value)
    return value.lstrip("0") or "0"


def lead_index(root, source, kind, sort_memory):
    """Build one compact coordinate index and ambiguity list per reference version.

    The original PVAR/BIM is scanned only on a cache miss. Genotype files and
    their variant order are untouched; these sorted views are for ID matching.
    """
    before = identity(source)
    meta = {"version": 3, "kind": kind, "source": list(before)}
    key = hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest()[:24]
    root = root / "index-v3" / kind
    root.mkdir(parents=True, exist_ok=True)
    target = root / key
    with locked(root / (key + ".lock")):
        ready = target / ".complete"
        try:
            saved = json.loads(ready.read_text())
            valid = saved["identity"] == meta and all(
                (target / name).stat().st_size == size
                for name, size in saved["files"].items()
            ) and {"reference.bgz", "reference.bgz.csi"} <= saved["files"].keys()
            if kind == "pvar":
                valid = valid and "ambiguous.snp" in saved["files"]
            if valid:
                valid = subprocess.run(["tabix", "-l", str(target / "reference.bgz")],
                                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0
        except (OSError, ValueError, KeyError, TypeError):
            valid = False
        if not valid:
            check_space(root, before[1] * 2)
            print(f"Build lead {kind} coordinate index: {source}", file=sys.stderr, flush=True)
            work = Path(tempfile.mkdtemp(prefix=key + ".tmp.", dir=root))
            try:
                compact = work / "rows.tsv"
                all_ids, multi_ids = work / "ids", work / "multi"
                columns = None
                rows = 0
                with reference_text(source) as lines, compact.open("w") as out, \
                        all_ids.open("w") as ids, multi_ids.open("w") as multi:
                    for line in lines:
                        if kind == "pvar" and line.startswith("##"):
                            continue
                        if kind == "pvar" and line.startswith("#CHROM"):
                            header = line.upper().split()
                            columns = [header.index(name) for name in ("#CHROM", "POS", "ID", "REF", "ALT")]
                            continue
                        if not line.strip():
                            continue
                        if kind == "pvar":
                            if columns is None:
                                raise ValueError(f"Invalid .pvar header: {source}")
                            # Discard INFO/sample payload without splitting every field.
                            fields = line.split(None, max(columns) + 1)
                            ch, pos, snp, ref, alt = (fields[i] for i in columns)
                            if snp not in ("", "."):
                                ids.write(snp + "\n")
                                if "," in alt:
                                    multi.write(snp + "\n")
                            row = [norm_chrom(ch), pos, snp, ref, alt]
                        else:
                            fields = line.split()
                            if len(fields) < 6:
                                continue
                            row = fields[:6]
                            ch, pos = row[0], row[3]
                            row[0] = norm_chrom(ch)
                        if not pos.isascii() or not pos.isdigit() or int(pos) <= 0:
                            continue
                        row[1 if kind == "pvar" else 3] = str(int(pos))
                        out.write("\t".join(row) + "\n")
                        rows += 1
                if kind == "pvar" and columns is None:
                    raise ValueError(f"Invalid .pvar header: {source}")
                env = dict(os.environ, LC_ALL="C")
                sort_args = ["sort", "--parallel=1", "-T", str(work), "-S", sort_memory]
                sorted_rows = work / "sorted.tsv"
                with sorted_rows.open("wb") as out:
                    if kind == "pvar":
                        out.write(b"#CHROM\tPOS\tID\tREF\tALT\n")
                        out.flush()
                    pos_col = "2" if kind == "pvar" else "4"
                    subprocess.run(sort_args + ["-t", "\t", "-k1,1V", f"-k{pos_col},{pos_col}n", str(compact)],
                                   stdout=out, env=env, check=True)
                cached = work / "reference.bgz"
                with cached.open("wb") as out:
                    subprocess.run(["bgzip", "-@", "1", "-c", str(sorted_rows)], stdout=out, check=True)
                # CSI also supports coordinates beyond the TBI 512-Mbp limit.
                subprocess.run(["tabix", "-C", "-s", "1", "-b", pos_col, "-e", pos_col,
                                "-S", "1" if kind == "pvar" else "0", str(cached)], check=True)
                if kind == "pvar":
                    sorted_ids, duplicates = work / "ids.sorted", work / "duplicates"
                    with sorted_ids.open("wb") as out:
                        subprocess.run(sort_args + [str(all_ids)], stdout=out, env=env, check=True)
                    with duplicates.open("wb") as out:
                        subprocess.run(["uniq", "-d", str(sorted_ids)], stdout=out, env=env, check=True)
                    with (work / "ambiguous.snp").open("wb") as out:
                        subprocess.run(sort_args + ["-u", str(multi_ids), str(duplicates)],
                                       stdout=out, env=env, check=True)
                if identity(source) != before:
                    raise RuntimeError("Reference changed while indexing; retry after its writer finishes")
                subprocess.run(["tabix", "-l", str(cached)], stdout=subprocess.DEVNULL, check=True)
                keep = {"reference.bgz", "reference.bgz.csi"}
                if kind == "pvar":
                    keep.add("ambiguous.snp")
                for path in work.iterdir():
                    if path.name not in keep:
                        path.unlink()
                (work / ".complete").write_text(json.dumps({"identity": meta, "rows": rows,
                    "files": {name: (work / name).stat().st_size for name in keep}}) + "\n")
                if target.exists():
                    shutil.rmtree(target)
                work.rename(target)
                print(f"Lead {kind} index ready: {rows:,} rows -> {target}", file=sys.stderr, flush=True)
            finally:
                if work.exists():
                    shutil.rmtree(work)
    print(target / "reference.bgz")


def seed_rsid(source, target):
    """Import old project caches once; an existing local generation is authoritative."""
    if not source.is_dir() or source.resolve() == target.resolve():
        return
    target.mkdir(parents=True, exist_ok=True)
    with locked(target / ".seed.lock"):
        for old in sorted(source.glob("v1.GRCh*/*")):
            if not old.is_dir():
                continue
            new = target / old.relative_to(source)
            if new.exists():
                continue
            new.parent.mkdir(parents=True, exist_ok=True)
            # Cooperate with the mapper and reference builder before copying.
            with locked(old / "reference.lock"), locked(old / "positions.lock"):
                files = [old / name for name in (
                    "dbsnp.tsv.gz", "dbsnp.tsv.gz.tbi", "dbsnp.tsv.gz.csi",
                    "positions.sqlite3", "positions.sqlite3-journal", "positions.sqlite3-wal",
                ) if (old / name).is_file()]
                if not files:
                    continue
                check_space(new.parent, sum(p.stat().st_size for p in files))
                work = Path(tempfile.mkdtemp(prefix=new.name + ".tmp.", dir=new.parent))
                try:
                    print(f"Import rsID cache: {old} -> {new}", file=sys.stderr, flush=True)
                    for path in files:
                        shutil.copy2(path, work / path.name)
                    database = work / "positions.sqlite3"
                    if database.exists():
                        # Recover an interrupted transaction on the COPY. Preserve
                        # the source database and its rollback journal together.
                        with closing(sqlite3.connect(database)) as conn:
                            checks = conn.execute("PRAGMA quick_check").fetchall()
                            if checks != [("ok",)]:
                                raise RuntimeError(f"Invalid rsID cache: {old}: {checks[:5]}")
                            conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
                    (work / ".seeded").write_text(str(old) + "\n")
                    work.rename(new)
                finally:
                    if work.exists():
                        shutil.rmtree(work)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    copy = commands.add_parser("stage")
    copy.add_argument("--root", type=Path, required=True)
    copy.add_argument("files", nargs="+")
    index = commands.add_parser("lead-index")
    index.add_argument("--root", type=Path, required=True)
    index.add_argument("--reference", type=Path, required=True)
    index.add_argument("--kind", choices=("pvar", "bim"), required=True)
    index.add_argument("--sort-memory", default="512M")
    seed = commands.add_parser("seed-rsid")
    seed.add_argument("--source", type=Path, required=True)
    seed.add_argument("--target", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "stage":
        stage(args.root, args.files)
    elif args.command == "lead-index":
        lead_index(args.root, args.reference, args.kind, args.sort_memory)
    else:
        seed_rsid(args.source, args.target)


if __name__ == "__main__":
    main()
