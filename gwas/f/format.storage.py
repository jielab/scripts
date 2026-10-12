#!/usr/bin/env python3
"""Stage reusable GWAS inputs on a local filesystem without changing originals."""

import argparse
from contextlib import closing, contextmanager
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
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
        valid = ready.is_file() and all(
            (target / p.name).is_file() and (target / p.name).stat().st_size == info[1]
            for p, info in zip(sources, before)
        )
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
    seed = commands.add_parser("seed-rsid")
    seed.add_argument("--source", type=Path, required=True)
    seed.add_argument("--target", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "stage":
        stage(args.root, args.files)
    else:
        seed_rsid(args.source, args.target)


if __name__ == "__main__":
    main()
