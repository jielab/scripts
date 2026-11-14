#!/usr/bin/env python3
"""Diagnose why metagenomic fragments remain unclassified / above species rank.

The report combines the exact existing Kraken split summary/counts with a bounded
streaming QC of the ORIGINAL single-end or paired FASTQs joined by ordinal to
``read-assignments.tsv.gz``.  It does not reclassify reads and does not require
loading the multi-terabyte Kraken database.

This is deliberately diagnostic rather than a generic FastQC replacement.  The
same read metrics are accumulated separately for Kraken ``unclassified``,
``higher_rank`` and ``assigned_species`` fragments, allowing poor quality / Ns /
low complexity to be tested as explanations for low taxonomic recall.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
from storage_layout import metadata_dir
from qc_paths import configured_split_root

from fastq_stream import fragments, read_id, new_digest, update_digest
from result_output import managed_outputs, prepare as prepare_result, finish as finish_result, replace_bool
from common import atomic_json, require_compute, load_config
from assignment_io import find_assignment_source, open_assignment_text

ADAPTERS = (b"AGATCGGAAGAGC", b"CTGTCTCTTATACACATCT")


def _entropy(seq: bytes) -> float:
    if not seq:
        return 0.0
    counts = Counter(seq.upper())
    n = len(seq)
    value = 0.0
    for count in counts.values():
        p = count / n
        value -= p * math.log2(p)
    return value


def _unique_kmer_fraction(seq: bytes, k: int = 15) -> float:
    seq = seq.upper()
    total = len(seq) - k + 1
    if total <= 0:
        return 1.0
    return len({seq[i:i+k] for i in range(total)}) / total


def _max_homopolymer(seq: bytes) -> int:
    seq = seq.upper()
    best = run = 0
    previous = None
    for base in seq:
        if base == previous:
            run += 1
        else:
            previous, run = base, 1
        best = max(best, run)
    return best


def _measure_fragment(pair):
    """Calculate sequence metrics once; reuse them across overlapping strata."""
    measured = []
    for rec in pair:
        if rec is None:
            continue
        seq = rec[1].rstrip(b"\r\n").upper()
        qual = rec[3].rstrip(b"\r\n")
        scores = [q - 33 for q in qual]
        if any(q < 0 or q > 93 for q in scores):
            raise RuntimeError("FASTQ quality contains a non-Phred+33 character")
        entropy = _entropy(seq)
        unique15 = _unique_kmer_fraction(seq)
        measured.append((
            len(seq), seq.count(b"G") + seq.count(b"C"), seq.count(b"N"),
            sum(q >= 20 for q in scores), sum(q >= 30 for q in scores), sum(scores),
            entropy, unique15, int(entropy < 1.2 or unique15 < 0.35),
            int(any(adapter in seq for adapter in ADAPTERS)),
            int(_max_homopolymer(seq) >= 15),
        ))
    return measured


class Metrics:
    __slots__ = ("fragments", "reads", "bases", "gc", "n", "q20", "q30", "qsum",
                 "low_complexity_reads", "adapter_reads", "long_homopolymer_reads",
                 "lengths", "entropy_sum", "unique15_sum")

    def __init__(self):
        self.fragments = self.reads = self.bases = 0
        self.gc = self.n = self.q20 = self.q30 = self.qsum = 0
        self.low_complexity_reads = self.adapter_reads = self.long_homopolymer_reads = 0
        self.lengths = Counter()
        self.entropy_sum = self.unique15_sum = 0.0

    def add_fragment(self, pair):
        self.add_measurements(_measure_fragment(pair))

    def add_measurements(self, measured):
        self.fragments += 1
        for (length, gc, n, q20, q30, qsum, entropy, unique15,
             low_complexity, adapter, homopolymer) in measured:
            self.reads += 1
            self.bases += length
            self.lengths[length] += 1
            self.gc += gc
            self.n += n
            self.qsum += qsum
            self.q20 += q20
            self.q30 += q30
            self.entropy_sum += entropy
            self.unique15_sum += unique15
            self.low_complexity_reads += low_complexity
            self.adapter_reads += adapter
            self.long_homopolymer_reads += homopolymer

    def report(self):
        return {
            "fragments": self.fragments,
            "reads": self.reads,
            "bases": self.bases,
            "mean_read_length": self.bases / self.reads if self.reads else None,
            "length_distribution": dict(sorted(self.lengths.items())),
            "gc_fraction": self.gc / self.bases if self.bases else None,
            "n_fraction": self.n / self.bases if self.bases else None,
            "mean_phred": self.qsum / self.bases if self.bases else None,
            "q20_fraction": self.q20 / self.bases if self.bases else None,
            "q30_fraction": self.q30 / self.bases if self.bases else None,
            "mean_sequence_entropy_bits": self.entropy_sum / self.reads if self.reads else None,
            "mean_unique_15mer_fraction": self.unique15_sum / self.reads if self.reads else None,
            "low_complexity_read_fraction": self.low_complexity_reads / self.reads if self.reads else None,
            "adapter_signature_read_fraction": self.adapter_reads / self.reads if self.reads else None,
            "homopolymer_ge15_read_fraction": self.long_homopolymer_reads / self.reads if self.reads else None,
        }


def _load_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Cannot read JSON {path}: {exc}") from exc


def _read_counts(path):
    totals = Counter()
    species = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"assigned_taxid", "assigned_rank", "name", "group", "status", "fragments"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise RuntimeError(f"counts.tsv missing {required}: {path}")
        for row in reader:
            n = int(row["fragments"])
            if n < 0:
                raise RuntimeError("counts.tsv contains a negative fragment count")
            totals[("status", row["status"])] += n
            totals[("group", row["group"])] += n
            totals[("rank", row["assigned_rank"])] += n
            if row["status"] == "assigned_species":
                species.append({"taxid": int(row["assigned_taxid"]), "name": row["name"],
                                "group": row["group"], "fragments": n})
    species.sort(key=lambda x: (-x["fragments"], x["taxid"]))
    return totals, species


def _exact_summary(summary, counts_path):
    totals, species = _read_counts(counts_path)
    total = int(summary["fragments"])
    status = {k[1]: v for k, v in totals.items() if k[0] == "status"}
    if total <= 0 or sum(status.values()) != total:
        raise RuntimeError("counts.tsv fragments do not match summary.json")
    groups = {k[1]: v for k, v in totals.items() if k[0] == "group"}
    ranks = {k[1]: v for k, v in totals.items() if k[0] == "rank"}
    thresholds = {}
    for threshold in (1, 2, 5, 10, 20, 50, 100, 500, 1000):
        kept = [row for row in species if row["fragments"] >= threshold]
        thresholds[str(threshold)] = {
            "species_taxa": len(kept),
            "fragments": sum(row["fragments"] for row in kept),
            "fraction_of_species_level_fragments": (sum(row["fragments"] for row in kept) /
                                                    max(1, status.get("assigned_species", 0))),
        }
    return {
        "fragments": total,
        "reads": int(summary.get("reads", total * (2 if summary.get("paired", True) else 1))),
        "total_bases": int(summary.get("total_bases", 0)),
        "status_fragments": status,
        "status_fraction": {k: v / total for k, v in status.items()},
        "group_fragments": groups,
        "rank_fragments": ranks,
        "species_taxa": len(species),
        "species_materialization_thresholds": thresholds,
        "top_species": species[:50],
    }


def diagnose(r1, r2, assignments, *, limit=2_000_000,
             expected_fragments=None, expected_digest=None):
    """Join an ordinal prefix, or validate all original bytes for ``limit=0``.

    A prefix is not a random sample and cannot establish that the unseen input
    or its gzip trailers are intact. Full scans require the saved input digest
    so changed sequences/qualities cannot pass just because read IDs match.
    """
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")
    if expected_fragments is not None and (isinstance(expected_fragments, bool) or
            not isinstance(expected_fragments, int) or expected_fragments < 1):
        raise ValueError("expected_fragments must be a positive integer")
    if not limit and not expected_digest:
        raise ValueError("Full-input diagnosis requires summary.json input_digest")
    metrics = defaultdict(Metrics)
    scanned = 0
    digest = new_digest() if not limit else None
    with open_assignment_text(Path(assignments)) as assignment_handle:
        reader = csv.DictReader(assignment_handle, delimiter="\t")
        required = {"ordinal", "r1_read_id", "r2_read_id", "assigned_taxid", "group", "status"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise RuntimeError(f"Assignments missing columns {required}")
        status_col = "kraken_status" if "kraken_status" in reader.fieldnames else "status"
        group_col = "kraken_group" if "kraken_group" in reader.fieldnames else "group"
        stream = fragments(r1, r2, limit)
        try:
            for scanned, pair in enumerate(stream, 1):
                if expected_fragments is not None and scanned > expected_fragments:
                    raise RuntimeError("FASTQ contains more fragments than the split summary")
                row = next(reader, None)
                if row is None:
                    raise RuntimeError(f"Assignments ended before FASTQ at ordinal {scanned}")
                if int(row["ordinal"]) != scanned:
                    raise RuntimeError(f"Assignment ordinal mismatch: {row['ordinal']} != {scanned}")
                a, b = pair
                rid1 = read_id(a[0]).decode("utf-8", "backslashreplace")
                rid2 = read_id(b[0]).decode("utf-8", "backslashreplace") if b else ""
                if row["r1_read_id"] != rid1 or row["r2_read_id"] != rid2:
                    raise RuntimeError(f"FASTQ/assignment read ID mismatch at ordinal {scanned}")
                if digest is not None:
                    update_digest(digest, a, b)
                status, group = row[status_col], row[group_col]
                measured = _measure_fragment((a, b))
                metrics["all"].add_measurements(measured)
                metrics["status:" + status].add_measurements(measured)
                metrics["group:" + group].add_measurements(measured)
                if scanned % 250_000 == 0:
                    print(f"QC joined {scanned:,} fragments", flush=True)
        finally:
            stream.close()
        if not limit and next(reader, None) is not None:
            raise RuntimeError("Assignments contain more rows than the full FASTQ input")
    if not scanned:
        raise RuntimeError("No FASTQ fragments scanned")
    if expected_fragments is not None:
        expected_scan = min(limit, expected_fragments) if limit else expected_fragments
        if scanned != expected_scan:
            raise RuntimeError(f"FASTQ fragment count differs from the split summary: {scanned} != {expected_scan}")
    if digest is not None and digest.hexdigest() != expected_digest:
        raise RuntimeError("FASTQ input digest differs from the split summary; original record bytes changed")
    return {key: value.report() for key, value in sorted(metrics.items())}, scanned


def _sample_parts(value):
    """Restrict sample selectors to one dataset/sample below configured roots."""
    parts = str(value).split("/")
    if len(parts) != 2 or any(x in ("", ".", "..") for x in parts) or "\\" in str(value):
        raise ValueError("--sample must be DATASET/SAMPLE without absolute or parent paths")
    return parts


def _validate_layout(summary, r2):
    total = int(summary["fragments"])
    paired = summary.get("paired")
    if paired is not None and not isinstance(paired, bool):
        raise ValueError("summary paired must be a JSON boolean")
    if paired is None:
        reads = int(summary.get("reads", total * (2 if r2 else 1)))
        if reads not in (total, total * 2):
            raise ValueError("summary reads is inconsistent with single/paired fragment counts")
        paired = reads == total * 2
    if paired != bool(r2):
        raise ValueError("Input FASTQ single/paired mode differs from the split summary")
    if int(summary.get("reads", total * (2 if paired else 1))) != total * (2 if paired else 1):
        raise ValueError("summary reads and fragments disagree with its single/paired mode")
    return paired


@managed_outputs
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replace", type=replace_bool, default=True)
    parser.add_argument("--config")
    parser.add_argument("--root")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--sample", metavar="DATASET/SAMPLE", help="Committed sample below the configured output_dir")
    source.add_argument("--output-dir", help="Existing sample output; infers summary/counts/assignments")
    parser.add_argument("--summary")
    parser.add_argument("--counts")
    parser.add_argument("--assignments", help="read-assignments.tsv.gz or species.tar")
    parser.add_argument("--r1")
    parser.add_argument("--r2")
    parser.add_argument("--limit", type=int, default=2_000_000,
                        help="QC prefix fragments; 0 scans all [2,000,000]")
    parser.add_argument("--out", help="Report directory [work_dir/qc/diagnose/DATASET/SAMPLE]")
    parser.add_argument("--allow-local", action="store_true")
    args = parser.parse_args(argv)
    require_compute(args.allow_local)
    if args.limit < 0:
        parser.error("--limit must be >=0")
    cfg = load_config(args.config, args.root)
    if args.sample:
        try:
            args.output_dir = str(configured_split_root(cfg).joinpath(*_sample_parts(args.sample)))
        except ValueError as exc:
            parser.error(str(exc))
    output = Path(args.output_dir).resolve() if args.output_dir else None
    if output is not None:
        if not (metadata_dir(output) / "complete.json").is_file():
            parser.error("--output-dir must contain a committed split (qc/complete.json)")
        original = _load_json(metadata_dir(output) / "provenance.json")["context"]["sample"]
        if args.r1 is None and args.r2 is None:
            args.r1, args.r2 = original["r1"], original.get("r2")
    if not args.r1:
        parser.error("provide --sample/--output-dir, or --r1 FASTQ [--r2 FASTQ]")
    summary_path = Path(args.summary) if args.summary else (metadata_dir(output) / "summary.json" if output else None)
    counts_path = Path(args.counts) if args.counts else (metadata_dir(output) / "counts.tsv" if output else None)
    if summary_path is None or counts_path is None:
        parser.error("provide --output-dir or both --summary and --counts")
    assignments = Path(args.assignments) if args.assignments else (find_assignment_source(output) if output else None)
    if assignments is None:
        parser.error("provide --assignments or --output-dir")
    summary = _load_json(summary_path)
    if args.sample and summary.get("sample_id") not in (None, args.sample):
        parser.error("Requested --sample differs from sample_id in the split summary")
    try:
        paired = _validate_layout(summary, args.r2)
    except ValueError as exc:
        parser.error(str(exc))
    if summary.get("paired") is None:
        summary["paired"] = paired
    if not args.limit:
        if summary.get("validated_scope") == "prefix" or summary.get("requested_limit", 0):
            parser.error("--limit 0 requires a split of the entire original input, not a prefix split")
        if not summary.get("input_digest"):
            parser.error("--limit 0 requires input_digest in the saved split summary")
    exact = _exact_summary(summary, counts_path)
    r1 = str(Path(args.r1).expanduser().resolve(strict=True))
    r2 = str(Path(args.r2).expanduser().resolve(strict=True)) if args.r2 else None
    if r2 and Path(r1).samefile(r2):
        parser.error("Paired inputs must be different FASTQ files")
    if args.out:
        out = Path(args.out).expanduser().resolve()
    else:
        try:
            parts = _sample_parts(summary.get("sample_id", args.sample or ""))
        except ValueError:
            parser.error("Provide --out when the split summary has no DATASET/SAMPLE identifier")
        out = Path(cfg["work_dir"]).joinpath("qc", "diagnose", *parts).resolve()
    if not prepare_result(out, replace=args.replace, cfg=cfg, protected=[r1, r2, assignments, summary_path, counts_path]):
        return
    qc, scanned = diagnose(r1, r2, assignments, limit=args.limit,
                           expected_fragments=exact["fragments"],
                           expected_digest=summary.get("input_digest"))
    report = {
        "schema": "metag-sample-diagnostic-v1",
        "sample_id": summary.get("sample_id"),
        "existing_result": exact,
        "qc_prefix_fragments": scanned,
        "qc_scope": {
            "selection": "entire_input" if not args.limit else "first_n_fragments",
            "requested_limit": args.limit,
            "fraction_of_split_fragments": scanned / exact["fragments"],
            "input_digest_validated": not bool(args.limit),
            "sampling_note": "A prefix is not a random or representative sample of the full library; unseen input bytes are not validated" if args.limit else
                             "All original FASTQ record bytes and fragment counts match the saved split summary",
        },
        "qc": qc,
        "interpretation": {
            "quality_test": "compare status:unclassified with status:assigned_species for Q30/N/complexity differences",
            "classification_test": "Similar measured QC does not establish the cause of unclassified reads; inspect reference coverage, classifier thresholds, conserved sequence and orthogonal alignments",
            "species_test": "higher_rank means taxonomic evidence exists but is not species-specific; do not count it as no-hit",
        },
    }
    out.mkdir(parents=True, exist_ok=True)
    atomic_json(out / "sample-diagnostic.json", report)
    rows = []
    for key, values in qc.items():
        rows.append({"stratum": key, **{k: v for k, v in values.items() if k != "length_distribution"}})
    fields = list(rows[0])
    with (out / "sample-qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    finish_result(out)
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
