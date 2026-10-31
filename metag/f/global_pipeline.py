#!/usr/bin/env python3
"""Separate global-index workflow, intentionally a different classifier."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
from types import SimpleNamespace

from common import CFG, DIR_LIB, DIR_OUT, THREADS, assert_reference_policy
from global_limits import ensure_memory_scope

KRAKEN_VERSION = "2.17.2"
KRAKEN_COMMIT = "8c190b1b668825935dbf6dee5f969227dc8269bb"
# The official 2.17.2 tag retains 2.17.1 in its VERSION file.  A version
# string alone therefore cannot identify the source tested for this route.
KRAKEN_REPORTED_VERSIONS = {KRAKEN_VERSION, "2.17.1"}
KRAKEN_REQUIRED_FILES = {
    "kraken2", "kraken2-build", "kraken2-inspect", "classify",
    "estimate_capacity", "build_db", "dump_table", "k2mask",
    "kraken2lib.pm", "build_kraken2_db.sh", "scan_fasta_file.pl",
    "make_seqid2taxid_map.pl", "mask_low_complexity.sh",
}
KRAKEN_REQUIRED_EXECUTABLES = {
    "kraken2", "kraken2-build", "kraken2-inspect", "classify",
    "estimate_capacity", "build_db", "dump_table", "k2mask",
}


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _reported_version(text):
    match = re.search(r"Kraken version\s+(\S+)", str(text), re.I)
    return match.group(1) if match else str(text).strip()


def _installation_receipt(executable):
    """Verify the complete pinned installation before executing its wrapper."""
    folder = executable.parent
    path = folder.parent / "install.json"
    try:
        receipt = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Missing or invalid pinned Kraken installation receipt: {path}. "
                           "Run ./3.global.sh install for the verified installation.") from exc
    if (not isinstance(receipt, dict) or receipt.get("source_tag") != KRAKEN_VERSION
            or receipt.get("source_commit") != KRAKEN_COMMIT):
        raise RuntimeError(f"Kraken installation source does not match the verified tag/commit: {path}")
    files = receipt.get("files")
    if not isinstance(files, dict) or not KRAKEN_REQUIRED_FILES.issubset(files):
        raise RuntimeError(f"Incomplete Kraken installation file manifest: {path}")
    for name, expected in files.items():
        if (not isinstance(name, str) or Path(name).name != name or name in (".", "..")
                or not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected)):
            raise RuntimeError(f"Invalid Kraken installation file manifest: {path}")
        candidate = folder / name
        if not candidate.is_file() or _file_sha256(candidate) != expected:
            raise RuntimeError(f"Kraken installation hash mismatch: {candidate}. "
                               "Run ./3.global.sh install to restore the pinned installation.")
    actual_files = {entry.name for entry in folder.iterdir() if entry.is_file()}
    if actual_files != set(files):
        raise RuntimeError(f"Kraken installation files differ from the receipt: {path}")
    declared = _reported_version(receipt.get("reported_version", ""))
    if declared not in KRAKEN_REPORTED_VERSIONS:
        raise RuntimeError(f"Unsupported reported version in pinned Kraken installation receipt: {path}")
    return receipt, path, declared


def positive_int(value):
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def nonnegative_float(value):
    try:
        number = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("must be finite and nonnegative")
    return number


def positive_float(value):
    number = nonnegative_float(value)
    if number == 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def confidence_value(value):
    number = nonnegative_float(value)
    if number > 1:
        raise argparse.ArgumentTypeError("confidence must lie between 0 and 1")
    return number


def parser():
    ap = argparse.ArgumentParser(prog="./3.global.sh", description=(
        "Global Kraken 2: one main classification pass per sample. This changes the "
        "classifier; minimizer sampling may change sensitivity. With no module, runs all: "
        "install -> plan -> build -> benchmark (first selected sample) -> split (all selected samples). "
        "Run 'install' to install the pinned tool."))
    sp = ap.add_subparsers(dest="command", required=True)
    for name, description in (
        ("all", "Run install, plan, build, first-sample benchmark, then split all selected samples (default)"),
        ("plan", "Read-only resource/reference/sample plan; no building or querying"),
        ("build", "Build/resume one global index from validated local references"),
        ("benchmark", "Classify a bounded prefix of one sample; do not publish species bins"),
        ("split", "Classify selected samples and split original FASTQ records by taxonomy"),
    ):
        sub = sp.add_parser(name, help=description)
        sub.add_argument("--db", type=Path, default=DIR_LIB / "kraken2_global/db")
        sub.add_argument("--threads", type=positive_int, default=THREADS)
        sub.add_argument("--kraken-dir", type=Path, help="Complete Kraken 2 2.17.2 installation directory")
        sub.add_argument("--min-free-gib", type=nonnegative_float,
                         default=float(CFG.get("MIN_FREE_GIB", "20")))
        sub.add_argument("--memory-reserve-gib", type=nonnegative_float, default=8,
                         help="Additional internal available-memory reserve; default 8 GiB")
        sub.add_argument("--memory-limit-gib", type=positive_float, default=40,
                         help="Hard cgroup memory limit for this job and all children; default 40 GiB; no swap")
        if name in ("all", "plan", "build"):
            sub.add_argument("--db-gib", type=positive_float, default=24,
                             help="Hash table budget in GiB; may downsample minimizers; default 24")
            sub.add_argument("--chunk-mib", type=positive_int, default=512,
                             help="Derived FASTA chunk size in MiB of bases; default 512")
            sub.add_argument("--no-masking", action="store_true",
                             help="Explicitly disable low-complexity masking (changes the method)")
        if name != "build":
            sub.add_argument("--dat-list", default="", help="Comma-separated datasets under DIR_RAW")
            sub.add_argument("--sample", default="", help="Unique dataset/run or run ID")
            sub.add_argument("--output-suffix", default="kraken2",
                             help="Separate output version suffix; default kraken2")
        if name in ("all", "benchmark", "split"):
            sub.add_argument("--confidence", type=confidence_value, default=0,
                             help="Confidence threshold; speed-oriented default 0")
            sub.add_argument("--quick", action=argparse.BooleanOptionalAction, default=True,
                             help="Stop native search once minimum hit groups are found; default enabled")
            sub.add_argument("--read-reports", action=argparse.BooleanOptionalAction, default=False,
                             help="Write per-read assignment tables; default off to save CPU and I/O")
            sub.add_argument("--minimum-hit-groups", type=positive_int, default=2)
            sub.add_argument("--minimum-base-quality", type=int, choices=range(0, 94), default=0,
                             metavar="0..93")
            sub.add_argument("--work-root", type=Path, default=DIR_OUT / "runs_global")
            sub.add_argument("--keep-work", action="store_true",
                             help="Keep compact taxids after complete output/benchmark")
        if name in ("all", "benchmark"):
            sub.add_argument("--benchmark-fragments", type=positive_int, default=1_000_000,
                             help="Fragment/pair prefix; default 1,000,000")
    return ap


def configure_tools(args, *, required):
    tool_root = Path(os.environ.get("METAG_TOOLS", str(Path.home()/".cache/metag/tools")))
    selected = args.kraken_dir or os.environ.get("METAG_KRAKEN2_DIR") or None
    pinned = tool_root / "kraken2" / KRAKEN_VERSION / "bin"
    if selected is not None:
        selected = Path(selected).expanduser().resolve()
        if not selected.is_dir():
            raise RuntimeError(f"Kraken installation directory does not exist: {selected}")
    elif pinned.is_dir():
        selected = pinned
    # An explicitly selected installation must never silently fall back to an
    # unrelated wrapper already in PATH.  A plan performs no tool execution.
    if selected is None:
        executable = shutil.which("kraken2")
    else:
        candidate = selected / "kraken2"
        executable = str(candidate) if candidate.is_file() and os.access(candidate, os.X_OK) else None
    if not required:
        return {"executable": executable, "preferred_pinned_directory": str(pinned),
                "required_source_tag": KRAKEN_VERSION, "required_source_commit": KRAKEN_COMMIT,
                "installation_receipt_required": True}
    if not executable:
        raise RuntimeError(f"Kraken 2 is missing{(' in '+str(selected)) if selected else ''}. "
                           "Run ./3.global.sh install, then repeat this command.")
    executable = Path(executable).resolve()
    receipt, receipt_path, declared = _installation_receipt(executable)
    probe = subprocess.run([str(executable), "--version"], capture_output=True, text=True, timeout=15)
    match = re.search(r"Kraken version\s+(\S+)", probe.stdout + probe.stderr, re.I)
    if probe.returncode or not match or match.group(1) != declared:
        raise RuntimeError(f"The pinned Kraken installation receipt reports {declared}; found "
                           f"{(match.group(1) if match else 'unknown version')} at {executable}. "
                           "Run ./3.global.sh install for the separate pinned installation.")
    for name in sorted(KRAKEN_REQUIRED_EXECUTABLES):
        path = executable.parent / name
        if not path.is_file() or not os.access(path, os.X_OK):
            raise RuntimeError(f"Incomplete Kraken 2 installation: executable {path} is missing")
    # Use matching wrappers/natives together, even when the initial PATH entry
    # was a symlink to one wrapper.
    os.environ["PATH"] = str(executable.parent) + os.pathsep + os.environ.get("PATH", "")
    args.kraken_installation = {
        "source_tag": receipt["source_tag"], "source_commit": receipt["source_commit"],
        "reported_version": match.group(1), "version_text": (probe.stdout + probe.stderr).strip(),
        "version_label_mismatch": match.group(1) != receipt["source_tag"],
        "executable": str(executable), "install_receipt": str(receipt_path),
        "install_receipt_sha256": _file_sha256(receipt_path),
    }
    return executable


def discover_samples(args):
    from stream_pipeline import discover
    return discover(SimpleNamespace(command="split", dat_list=args.dat_list,
                                    sample=args.sample, output_suffix=args.output_suffix))


def plan(args):
    from global_build import plan_build
    from global_resources import memory_status
    report = {"workflow": "global Kraken 2", "classifier_is_changed": True,
              "tools": configure_tools(args, required=False), "memory": memory_status()}
    report["resource_limits"] = {"memory_limit_gib": args.memory_limit_gib,
                                 "swap_limit_gib": 0, "outside_headroom_gib": 8}
    errors = []
    try:
        report["build"] = plan_build(args.db_gib, args.min_free_gib, args.memory_reserve_gib, args.db,
                                     chunk_bases=args.chunk_mib * (1 << 20))
        report["build"]["parameters"].update(masking=not args.no_masking,
                                             chunk_bases=args.chunk_mib * (1 << 20))
    except (RuntimeError, OSError, ValueError) as exc:
        report["build_error"] = str(exc)
        errors.append("reference/build plan")
    try:
        samples = discover_samples(args)
        report["samples"] = {"count": len(samples), "datasets": dict(Counter(s["dataset"] for s in samples)),
            "input_file_bytes": sum(p.stat().st_size for s in samples for p in (s["r1"], s["r2"]) if p),
            "fragment_counts": "not counted; compressed bytes cannot establish total fragment count",
            "output_suffix": args.output_suffix,
            "first_sample": {k: str(v) if isinstance(v, Path) else v for k, v in samples[0].items()} if samples else None}
    except (RuntimeError, OSError, ValueError) as exc:
        report["sample_error"] = str(exc)
        errors.append("sample discovery")
    report["plan_complete"] = not errors
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 2 if errors else 0


def parse_args(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or (argv[0].startswith("-") and argv[0] not in ("-h", "--help")):
        argv.insert(0, "all")
    return parser().parse_args(argv)


def run_all(args):
    from global_build import build_database, database_snapshot
    from global_query import benchmark, run_sample

    print("[1/5 install] Verify/install pinned Kraken 2", flush=True)
    if args.kraken_dir or os.environ.get("METAG_KRAKEN2_DIR"):
        # An explicit installation is validated, never silently replaced.
        configure_tools(args, required=True)
    else:
        subprocess.run(["bash", str(Path(__file__).with_name("install_global_kraken.sh")),
                        "--threads", str(args.threads)], check=True)
    print("[2/5 plan] Reference, samples and resources", flush=True)
    status = plan(args)
    if status:
        return status
    samples = discover_samples(args)
    if not samples:
        raise RuntimeError("No samples selected")
    executable = configure_tools(args, required=True)
    print("[3/5 build] Build/resume global database", flush=True)
    args.chunk_bases = args.chunk_mib * (1 << 20)
    build_database(args)
    snapshot = database_snapshot(args.db)
    print(f"[4/5 benchmark] {samples[0]['name']}; first {args.benchmark_fragments:,} fragments", flush=True)
    benchmark(samples[0], snapshot, executable, args)
    print(f"[5/5 split] {len(samples)} selected samples", flush=True)
    for index, sample in enumerate(samples, 1):
        print(f"[Sample {index}/{len(samples)}] {sample['name']}", flush=True)
        run_sample(sample, snapshot, executable, args)
    print("PASS all global stages", flush=True)
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parse_args(argv)
    assert_reference_policy()
    if args.command != "plan":
        status = ensure_memory_scope(args, argv)
        if status is not None:
            return status
    if args.command == "all":
        return run_all(args)
    if args.command == "plan":
        return plan(args)
    executable = configure_tools(args, required=True)
    if args.command == "build":
        from global_build import build_database
        args.chunk_bases = args.chunk_mib * (1 << 20)
        result = build_database(args)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    from global_build import database_snapshot
    from global_query import benchmark, run_sample
    snapshot = database_snapshot(args.db)
    samples = discover_samples(args)
    if not samples:
        raise RuntimeError("No samples selected")
    if args.command == "benchmark" and len(samples) != 1:
        raise RuntimeError("Benchmark requires exactly one sample; pass --sample dataset/run")
    print(f"Global database: {snapshot['path']}\nClassifier: Kraken 2 source tag {KRAKEN_VERSION} "
          f"(reports {args.kraken_installation['reported_version']}; CPU)\n"
          f"Selected samples: {len(samples)}; one main query per sample; output suffix={args.output_suffix}", flush=True)
    if args.command == "benchmark":
        benchmark(samples[0], snapshot, executable, args)
    else:
        for index, sample in enumerate(samples, 1):
            print(f"[Sample {index}/{len(samples)}] {sample['name']}", flush=True)
            run_sample(sample, snapshot, executable, args)
    return 0


def _terminate(signum, frame):
    raise KeyboardInterrupt(f"Received signal {signum}")


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, _terminate)
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Stopped. Completed database chunks/query receipts are retained for retry.", file=sys.stderr)
        sys.exit(130)
    except (RuntimeError, OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
