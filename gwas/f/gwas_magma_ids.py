#!/usr/bin/env python3
"""Map coordinate GWAS IDs to MAGMA reference rsIDs in the declared build.

The reusable SQLite cache contains only rsIDs present in the MAGMA BIM.
Coordinate matching requires both alleles; ambiguous matches are discarded.
"""
import argparse
from contextlib import closing
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import sys
import shutil
import tempfile
import time


def chromosome(value):
    value = value.removeprefix("chr").upper()
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(value) or int(value)


def signature(*paths):
    return hashlib.sha256(json.dumps([
        (str(Path(p).resolve()), Path(p).stat().st_size, Path(p).stat().st_mtime_ns)
        for p in paths
    ]).encode()).hexdigest()[:24]


def reference_cache(bim, dbsnp, cache_dir):
    target = Path(cache_dir) / (signature(bim, dbsnp) + ".sqlite3")
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(str(target) + ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if target.exists():
            return target
        print("Build MAGMA rsID cache (one dbSNP scan): " + str(target), file=sys.stderr, flush=True)
        ids = set()
        with open(bim, buffering=8 * 1024 * 1024) as source:
            for line in source:
                snp = line.split()[1]
                if snp.startswith("rs") and snp[2:].isdigit():
                    ids.add(int(snp[2:]))
        if not ids:
            raise ValueError("Coordinate mapping requires rsIDs in the MAGMA BIM")
        temp = Path(str(target) + f".tmp.{os.getpid()}")
        workspace = tempfile.TemporaryDirectory(prefix="magma-rsid-build-", dir="/tmp")
        local_db = Path(workspace.name) / "reference.sqlite3"
        try:
            # A Connection context manager commits but does not close the file.
            # Close before publishing: renaming an open database on WSL/DrvFS
            # can leave the destination unavailable to subsequent readers.
            with closing(sqlite3.connect(local_db)) as conn:
                conn.execute("PRAGMA journal_mode=OFF")
                conn.execute("PRAGMA synchronous=OFF")
                conn.execute("CREATE TABLE variants (chr INTEGER, pos INTEGER, a TEXT, b TEXT, rs INTEGER)")
                batch = []
                started = time.monotonic()
                with open(dbsnp, "rb", buffering=8 * 1024 * 1024) as compressed, gzip.open(compressed, "rt") as source:
                    for row_number, line in enumerate(source, 1):
                        if row_number % 5000000 == 0:
                            print(f"MAGMA cache: scanned {row_number:,} dbSNP rows in {time.monotonic()-started:.0f}s", file=sys.stderr, flush=True)
                        ch, pos, snp, ref, alt = line.rstrip().split("\t")[:5]
                        if not snp.startswith("rs") or not snp[2:].isdigit() or int(snp[2:]) not in ids:
                            continue
                        try:
                            ch = chromosome(ch)
                        except ValueError:
                            continue
                        if not 1 <= ch <= 23:
                            continue
                        for allele in alt.split(","):
                            a, b = sorted((ref.upper(), allele.upper()))
                            batch.append((ch, int(pos), a, b, int(snp[2:])))
                        if len(batch) >= 10000:
                            conn.executemany("INSERT INTO variants VALUES (?,?,?,?,?)", batch)
                            batch.clear()
                    conn.executemany("INSERT INTO variants VALUES (?,?,?,?,?)", batch)
                conn.execute("CREATE INDEX coordinates ON variants(chr,pos)")
                conn.commit()
            shutil.copyfile(local_db, temp)
            temp.replace(target)
            print("MAGMA rsID cache ready: " + str(target), file=sys.stderr, flush=True)
        finally:
            temp.unlink(missing_ok=True)
            workspace.cleanup()
    return target


def map_rows(rows, output, snploc, audit, cache):
    counts = dict(rsid=0, mapped=0, unmatched=0, ambiguous=0)
    conn = sqlite3.connect(Path(cache).resolve().as_uri() + "?mode=ro&immutable=1", uri=True) if cache else None
    if conn:
        # reference_cache publishes a closed database atomically and never edits
        # it. Immutable readers avoid per-variant locks/stat calls on DrvFS.
        # mmap shares cached pages between workers without a private DB copy.
        conn.execute("PRAGMA mmap_size=2147483648")
    previous = None
    matches = {}
    try:
        with open(rows) as source, open(output, "w") as out, open(snploc, "w") as loc:
            for line in source:
                snp, p, n, ch, pos, ea, nea = line.rstrip("\n").split("\t")
                ch, pos = chromosome(ch), int(pos)
                if snp.startswith("rs") and snp[2:].isdigit():
                    counts["rsid"] += 1
                else:
                    if conn is None:
                        raise ValueError("Non-rsID input requires a dbSNP mapping cache")
                    key = (ch, pos)
                    if key != previous:
                        matches = {}
                        for a, b, rs in conn.execute("SELECT a,b,rs FROM variants WHERE chr=? AND pos=?", key):
                            matches.setdefault((a, b), set()).add(rs)
                        previous = key
                    candidates = matches.get(tuple(sorted((ea.upper(), nea.upper()))), set())
                    if len(candidates) != 1:
                        counts["ambiguous" if candidates else "unmatched"] += 1
                        continue
                    snp = "rs" + str(next(iter(candidates)))
                    counts["mapped"] += 1
                out.write(f"{snp}\t{p}\t{n}\n")
                loc.write(f"{snp}\t{ch}\t{pos}\n")
    finally:
        if conn:
            conn.close()
    with open(audit, "w") as out:
        out.write("status\trows\n")
        for key, value in counts.items():
            out.write(f"{key}\t{value}\n")
    print("MAGMA IDs: " + json.dumps(counts), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("rows", "output", "snploc", "audit", "bim", "dbsnp", "cache"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    with open(args.rows) as source:
        need_mapping = any(not (s := line.split("\t", 1)[0]).startswith("rs") or not s[2:].isdigit() for line in source)
    cache = reference_cache(args.bim, args.dbsnp, args.cache) if need_mapping else None
    map_rows(args.rows, args.output, args.snploc, args.audit, cache)


if __name__ == "__main__":
    main()
