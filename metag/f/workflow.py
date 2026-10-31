#!/usr/bin/env python3
"""Command orchestration: automatic startup checks, then the named workload."""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
import signal
from pathlib import Path

from common import ROOT, DIR_REF, DIR_LIB, DIR_OUT, DIR_RAW, CFG, PART_COMPRESSED_GIB, assert_reference_policy, require, run
import gpu_check
import ref_prepare
import check_ref
from build_limits import memory_gib, in_workflow_scope, run_workflow_limited


def check_libraries():
    from library import completed_library, file_state
    from metacache_query import database_parts
    saved = completed_library(DIR_LIB/"sourmash")
    if saved.get("output") != file_state(DIR_LIB/"sourmash/references.zip"):
        raise RuntimeError("Sourmash library differs from its completion receipt; run ./1.library.sh")
    parts = database_parts()
    for name in ("nodes.dmp", "names.dmp", "merged.dmp"):
        if not (DIR_LIB/"metacache/taxonomy"/name).is_file():
            raise RuntimeError(f"Library taxonomy missing: {name}; run ./1.library.sh")
    print(f"PASS: completed sourmash library and {len(parts)} MetaCache partitions")


def check_software(command):
    assert_reference_policy()
    print(f"Reference FASTA: {DIR_REF}\nPermanent libraries: {DIR_LIB}\nAnalysis/work files: {DIR_OUT}")
    print("Reference FASTA policy: KEEP (source files are never deleted)")
    if command == "build-kraken":
        programs, plugins = ("kraken2-build", "pigz", "perl"), ()
    else:
        programs = ["sourmash"]
        plugins = ("singlesketch", "fastgather")
        if command == "build-sourmash":
            plugins = ()
        elif command in ("classify", "split"):
            plugins = ("singlesketch", "fastgather")
        elif command == "build-metacache":
            programs, plugins = [], ()
    for program in programs:
        require(program)
    for plugin in plugins:
        proc = subprocess.run(["sourmash", "scripts", plugin, "--help"],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if proc.returncode:
            raise RuntimeError(f"Missing sourmash branchwater command {plugin}; run ./install.sh")
    if command not in ("build-sourmash", "build-kraken"):
        mc = Path(os.environ.get("METAG_TOOLS", Path.home()/".cache/metag/tools"))/"metacache"
        for name in ("metacache",):
            if not os.access(mc/name, os.X_OK):
                raise RuntimeError(f"Missing MetaCache program {mc/name}; run ./install.sh")
        print(f"MetaCache CPU: {mc/'metacache'}")
        if os.access(mc/"metacache_gpu", os.X_OK):
            print(f"MetaCache GPU: {mc/'metacache_gpu'}")
        elif command in ("build-db", "build-metacache"):
            raise RuntimeError("MetaCache GPU builder required; library builds do not fall back to CPU")
        else:
            print("MetaCache GPU: unavailable; queries can use CPU")
    print("PASS: required software")


def positive_gib(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("partition size must be a positive finite number")
    return number


def nonnegative_int(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("database size must be >= 0")
    return number


def worker_count(value):
    number=int(value)
    if not 1<=number<=32:
        raise argparse.ArgumentTypeError('workers must be between 1 and 32')
    return number


def parser():
    ap = argparse.ArgumentParser(prog=os.environ.get("METAG_PROG", "./metag.sh"), description="Library building and sample preparation")
    sp = ap.add_subparsers(dest="command", required=True)
    p = sp.add_parser("check-ref", help="Check local references and cache the validated inventory/mapping")
    p.add_argument("--force", action="store_true", help="rerun the complete reference check")
    for command, aliases, description in (
        ("build-db", ["prepare"], "Build sourmash and partitioned MetaCache databases"),
        ("build-sourmash", [], "Build the reference FracMinHash sketch database"),
        ("build-metacache", [], "Build the partitioned per-read classification database"),
    ):
        p = sp.add_parser(command, aliases=aliases, description=description, help=description)
        if command != "build-metacache":
            p.add_argument("--workers",type=worker_count,help="parallel CPU processes for sourmash (default: f/metag.conf)")
        if command != "build-sourmash":
            p.add_argument("--prepare-workers",type=int,choices=range(1,5),default=2,
                           help="parallel input preparation processes [%(default)s]")
            p.add_argument("--memory-limit-gb", type=memory_gib, default=16,
                           help="entire build workflow RAM cap, 1-16 GiB; swap disabled [%(default)s]")
            p.add_argument("--part-gib", type=positive_gib, default=PART_COMPRESSED_GIB,
                           help="compressed reference GiB per MetaCache partition [%(default)s]")
    p = sp.add_parser("classify", help="Classify and split one sample")
    p.add_argument("--r1", type=Path, required=True)
    p.add_argument("--r2", type=Path)
    p.add_argument("--name", required=True)
    p.add_argument("--out-split", type=Path)
    p = sp.add_parser("split", help="Classify and split datasets")
    p.add_argument("--dat-list", default="")
    p = sp.add_parser("build-kraken", help="Build the optional Kraken2 comparator database")
    p.add_argument("--max-db-gb", type=nonnegative_int, default=72)
    return ap


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in ("check_ref", "ref-check"):
        argv[0] = "check-ref"  # Accept the earlier spellings without advertising them.
    if argv and argv[0] in ("split", "classify"):
        from stream_pipeline import cli as split_cli
        return split_cli(argv)
    # Retain old diagnostic commands for existing scripts, outside the main usage.
    if argv and argv[0] in ("check", "gpu-check"):
        if len(argv) != 1:
            raise RuntimeError(f"{argv[0]} does not take arguments")
        if argv[0] == "check":
            check_software("check")
        elif argv[0] == "gpu-check":
            gpu_check.main()
        return 0
    args = parser().parse_args(argv)
    command = "build-db" if args.command == "prepare" else args.command
    if command == "split":
        # Public batch entry points must not silently lose singleton/duplicate
        # inputs or map distinct subfolders to the same sample output name.
        # Keep this orchestration-only check outside classification fingerprints.
        from sample_discovery import validate_split_datasets
        validate_split_datasets(DIR_RAW, args.dat_list)
    if command == "check-ref":
        return check_ref.main(force=args.force)
    if args.command == "prepare":
        print("prepare -> build-db (reference preparation + BOTH database builds)")
    total = 5 if command == "build-db" else 4

    def stage(number, title):
        print(f"\n[{number}/{total}] {title}", flush=True)

    stage(1, "Check software and reference retention policy")
    check_software(command)
    stage(2, "Check GPU / CUDA / available VRAM")
    gpu_check.main()
    if command in ("build-db", "build-sourmash", "build-metacache"):
        stage(3, "Reuse check-ref result, or prepare references if needed")
        check_ref.main()
    elif command in ("classify", "split"):
        stage(3, "Check completed permanent libraries")
        check_libraries()
    else:
        stage(3, "Check local reference files")
        ref_prepare.check_references()

    if command in ("build-db", "build-sourmash"):
        import sourmash_build
        stage(4, "BUILD DATABASE: sourmash reference FracMinHash sketches")
        if args.workers is None:
            sourmash_build.main()
        else:
            sourmash_build.main(args.workers)
    if command in ("build-db", "build-metacache"):
        import metacache_build
        stage(total, "BUILD DATABASE: partitioned MetaCache per-read classification index")
        metacache_build.main(args.part_gib, memory_gb=args.memory_limit_gb, prepare_workers=args.prepare_workers)
    elif command in ("classify", "split"):
        stage(4, "Classify reads and write species bins")
        run([sys.executable, ROOT/"f/pipeline.py", command, *argv[1:]])
    elif command == "build-kraken":
        stage(4, "BUILD DATABASE: optional Kraken2 comparator")
        run(["bash", ROOT/"f/kraken_build.sh", args.max_db_gb])
    print(f"\nPASS: {command} complete", flush=True)
    return 0


def cli(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in ('build-db','prepare','build-metacache','build-sourmash'):
        args=parser().parse_args(argv)
        limit=getattr(args,'memory_limit_gb',16)
        if not in_workflow_scope(limit):
            return run_workflow_limited([sys.executable, str(Path(__file__).resolve()), *argv],
                                        memory_gb=limit).returncode
        def stopped(signum, frame):
            raise KeyboardInterrupt
        signal.signal(signal.SIGTERM, stopped)
    return main(argv)


if __name__ == "__main__":
    try:
        raise SystemExit(cli())
    except KeyboardInterrupt:
        print("Stopped by request; completed batches are retained. Rerun ./1.library.sh to resume.",file=sys.stderr)
        raise SystemExit(130)
    except (RuntimeError, OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)

