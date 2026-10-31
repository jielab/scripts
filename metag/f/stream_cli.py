"""Shared public options for the bounded, full-reference split workflow."""
import argparse
import os
from pathlib import Path
from common import CFG, DIR_OUT


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def nonnegative_gib(value):
    import math
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("must be a finite nonnegative GiB value")
    return number


def batch_value(value):
    number = positive_int(value)
    if number < 512:
        raise argparse.ArgumentTypeError("batch size must be >=512 (MetaCache batches are constrained by query windows)")
    return number


def add_options(parser):
    parser.add_argument("--sample", default="", help="select a discovered run ID or dataset/run (must be unique)")
    parser.add_argument("--plan", action="store_true", help="show real inputs, outputs and completed databases, without querying")
    parser.add_argument("--benchmark", action="store_true", help="test one auto-discovered sample using a bounded prefix; never publish species bins")
    parser.add_argument("--benchmark-fragments", type=positive_int, default=100000)
    parser.add_argument("--benchmark-parts", type=positive_int, default=3)
    parser.add_argument("--benchmark-batches", default="512,4096,8192")
    parser.add_argument("--batch-size", type=batch_value, default=batch_value(CFG.get("QUERY_BATCH_SIZE", "4096")))
    parser.add_argument("--backend", choices=("auto", "gpu", "cpu"), default="gpu",
                        help="query backend (default: gpu; no automatic CPU fallback)")
    parser.add_argument("--query-cache", choices=("auto", "none", "fasta"), default=CFG.get("QUERY_CACHE", "auto"))
    parser.add_argument("--cache-gib", type=nonnegative_gib, default=float(CFG.get("QUERY_CACHE_GIB", "100")))
    parser.add_argument("--min-free-gib", type=nonnegative_gib, default=float(CFG.get("MIN_FREE_GIB", "20")))
    parser.add_argument("--work-root", type=Path, default=DIR_OUT / "runs_v2")
    parser.add_argument("--output-suffix", default="", help="optional suffix for new output versions; existing outputs are never overwritten")
    parser.add_argument("--screen", action="store_true", help="also produce the optional sourmash sample report; full-reference query scope is unchanged")
    parser.add_argument("--keep-work", action="store_true", help="retain binary state/query cache after a completed split (default: clean completed scratch)")


def parser():
    ap = argparse.ArgumentParser(prog=os.environ.get("METAG_PROG", "metag"), description="Full-reference species split with compact transactional merging")
    sp = ap.add_subparsers(dest="command", required=True)
    batch = sp.add_parser("split", help="discover FASTQs automatically")
    batch.add_argument("--dat-list", default="",
                       help="comma-separated datasets (default: all datasets under DIR_RAW)")
    add_options(batch)
    single = sp.add_parser("classify", help="explicit FASTQ input")
    single.add_argument("--r1", type=Path, required=True)
    single.add_argument("--r2", type=Path)
    single.add_argument("--name", required=True)
    single.add_argument("--out-split", type=Path)
    add_options(single)
    return ap
