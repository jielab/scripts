#!/usr/bin/env python3
"""Independent human-read validation with BWA-MEM2.

Use --scope candidates to check the existing human split FASTQs without the
original raw reads; use --scope full to screen all original reads and compare
against committed assignments. Every selected record is validated, with mates
kept together. Candidate-only support cannot measure missed human reads.

A human-only GRCh38 reference establishes *human compatibility*. For a stricter
species-specific check, align to a competitive human+chimp+gorilla reference and
set --contig-regex to the prefix used for human contigs (for example '^hs_').
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import json
import mmap
from pathlib import Path
from storage_layout import metadata_dir
from qc_paths import configured_split_root

from fastq_stream import fragments, new_digest, update_digest, read_id
from result_output import managed_outputs, prepare as prepare_result, finish as finish_result, replace_bool
from common import atomic_json, load_config, require_compute
from qc_human_screen import _pair_positive, host_context, screen_human_fastqs
from assignment_io import open_assignment_text, find_assignment_source
from data_readme import write_validation_readme


def count_input(r1, r2):
    digest = new_digest()
    count = 0
    for a, b in fragments(r1, r2):
        count += 1
        update_digest(digest, a, b)
        if count % 1_000_000 == 0:
            print(f"FASTQ preflight {count:,} fragments", flush=True)
    if not count:
        raise RuntimeError("FASTQ pair is empty")
    return count, digest.hexdigest()


def count_assignments(path):
    count = 0
    with open_assignment_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "ordinal" not in reader.fieldnames:
            raise RuntimeError("Assignments lack ordinal column")
        for count, row in enumerate(reader, 1):
            if int(row["ordinal"]) != count:
                raise RuntimeError(f"Assignments ordinal mismatch at {count}")
    if not count:
        raise RuntimeError("Assignment table is empty")
    return count


def compare_assignments(path, calls_path, expected, pair_policy, r1=None, r2=None):
    table = Counter()
    mate_table = Counter()
    kraken_human = kraken_human_aligned = alignment_positive = alignment_positive_kraken_human = 0
    aligned_unclassified = aligned_higher_rank = aligned_other_species = 0
    stream = fragments(r1, r2) if r1 else None
    if Path(calls_path).stat().st_size != expected:
        raise RuntimeError("Human call file length differs from expected fragments")
    with Path(calls_path).open("rb") as raw:
        calls = mmap.mmap(raw.fileno(), expected, access=mmap.ACCESS_READ)
        try:
            with open_assignment_text(path) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {"ordinal", "assigned_taxid", "group", "status"}
                if not reader.fieldnames or not required <= set(reader.fieldnames):
                    raise RuntimeError(f"Assignments missing required columns {required}: {path}")
                count = 0
                for row in reader:
                    count += 1
                    ordinal = int(row["ordinal"])
                    if ordinal != count or ordinal > expected:
                        raise RuntimeError(f"Assignments ordinal mismatch at row {count}: {ordinal}")
                    if stream is not None:
                        pair = next(stream, None)
                        records = () if pair is None else ((pair[0],) if pair[1] is None else pair)
                        if r2 is None and row.get("r2_read_id"):
                            raise RuntimeError("Single-end FASTQ cannot be compared with paired assignments")
                        if pair is None or any(row.get(col) != read_id(rec[0]).decode("utf-8", "backslashreplace")
                                               for col, rec in zip(("r1_read_id", "r2_read_id"), records)):
                            raise RuntimeError(f"FASTQ/assignment read ID mismatch at ordinal {ordinal}")
                    bits = calls[ordinal - 1]
                    positive = _pair_positive(bits, pair_policy)
                    is_human = row.get("kraken_group", row["group"]) == "human"
                    status = row.get("kraken_status", row["status"])
                    label = "kraken_human" if is_human else status
                    table[("aligned" if positive else "not_aligned", label)] += 1
                    evidence = ("both_mates" if bits & 3 == 3 else
                                "one_mate" if bits & 3 else "neither_mate") if r2 else (
                                "single_read_pass" if bits & 1 else "single_read_fail")
                    mate_table[(label, row.get("kraken_group", row["group"]), evidence)] += 1
                    kraken_human += int(is_human)
                    kraken_human_aligned += int(is_human and positive)
                    alignment_positive += int(positive)
                    alignment_positive_kraken_human += int(positive and is_human)
                    aligned_unclassified += int(positive and status == "unclassified")
                    aligned_higher_rank += int(positive and status == "higher_rank")
                    aligned_other_species += int(positive and status == "assigned_species" and not is_human)
                if count != expected:
                    raise RuntimeError(f"Assignments rows differ from FASTQ fragments: {count} != {expected}")
                if stream is not None and next(stream, None) is not None:
                    raise RuntimeError("FASTQ contains more fragments than assignments")
        finally:
            calls.close()
            if stream is not None:
                stream.close()
    return {
        "fragments": expected,
        "kraken_human_fragments": kraken_human,
        "kraken_human_alignment_supported": kraken_human_aligned,
        "kraken_human_alignment_support_fraction": kraken_human_aligned / kraken_human if kraken_human else None,
        "alignment_positive_fragments": alignment_positive,
        "alignment_positive_also_kraken_human": alignment_positive_kraken_human,
        "kraken_human_fraction_of_alignment_positive": alignment_positive_kraken_human / alignment_positive if alignment_positive else None,
        "alignment_positive_kraken_unclassified": aligned_unclassified,
        "alignment_positive_kraken_higher_rank": aligned_higher_rank,
        "alignment_positive_kraken_other_species": aligned_other_species,
        "crosstab": [{"alignment": key[0], "kraken": key[1], "fragments": value}
                     for key, value in sorted(table.items())],
        "mate_evidence": [{"kraken": key[0], "kraken_group": key[1], "alignment_evidence": key[2], "fragments": value}
                          for key, value in sorted(mate_table.items())],
        "interpretation": "agreement against an independent alignment screen, not a gold-standard sensitivity/PPV estimate"
    }


def _safe_sample_id(value):
    path = Path(value)
    if not value or path.is_absolute() or any(part in {".", ".."} for part in value.split("/")):
        raise ValueError("Sample ID must be a relative DATASET/SAMPLE path without . or ..")
    return path


def candidate_sources(output):
    """Discover flat human bins; above-species mammal calls are not human calls."""
    output = Path(output)
    root = output / "human"
    sources, recognised = [], set()
    if not root.is_dir():
        return sources
    for left in sorted(root.glob("*.R1.fastq.gz")):
        stem = left.name.removesuffix(".R1.fastq.gz")
        right = left.with_name(stem + ".R2.fastq.gz")
        if not right.is_file():
            raise RuntimeError(f"Human bin is missing R2: {left}")
        sources.append({"r1": str(left.resolve(strict=True)), "r2": str(right.resolve(strict=True)),
                        "bin": "human/" + stem})
        recognised.update((left, right))
    for single in sorted(root.glob("*.SE.fastq.gz")):
        sources.append({"r1": str(single.resolve(strict=True)), "r2": None,
                        "bin": "human/" + single.name.removesuffix(".SE.fastq.gz")})
        recognised.add(single)
    other = [path for path in root.rglob("*") if path.is_file() and
             path.name.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")) and path not in recognised]
    if other:
        raise RuntimeError("Unrecognised/orphan human FASTQ; supply --target-fastq explicitly: " + str(other[0]))
    if len({bool(source["r2"]) for source in sources}) > 1:
        raise RuntimeError("Human bins mix single-end and paired-end inputs; validate them separately")
    return sorted(sources, key=lambda row: row["r1"])


def count_sources(sources):
    digest = new_digest()
    total = 0
    for source in sources:
        count = 0
        for a, b in fragments(source["r1"], source.get("r2")):
            count += 1
            total += 1
            update_digest(digest, a, b)
        source["fragments"] = count
    return total, digest.hexdigest()


def write_candidate_evidence(sources, calls_path, expected, pair_policy, destination, evidence_database):
    """Join bounded, ordinal-sorted alignment evidence to original candidate IDs."""
    import sqlite3
    per_bin = []
    evidence = sqlite3.connect(str(evidence_database))
    evidence.row_factory = sqlite3.Row
    alignments = iter(evidence.execute("SELECT * FROM alignments ORDER BY ordinal, mate"))
    attributes = ("flag", "contig", "position", "cigar", "mapq", "aligned_fraction", "identity", "passed")
    try:
        with Path(calls_path).open("rb") as calls, Path(destination).open("w", newline="", encoding="utf-8") as out:
            fields = (["ordinal", "bin", "r1_read_id", "r2_read_id"] +
                      [f"r{mate}_{key}" for mate in (1, 2) for key in attributes] + ["pair_pass", "evidence"])
            writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            ordinal = 0
            for source in sources:
                counts = Counter()
                for a, b in fragments(source["r1"], source.get("r2")):
                    raw = calls.read(1)
                    if not raw:
                        raise RuntimeError("Candidate call file ended before input FASTQ")
                    bits = raw[0]
                    ordinal += 1
                    label = ("both_mates" if bits & 3 == 3 else "one_mate" if bits & 3 else "neither_mate") if b else (
                             "single_read_pass" if bits & 1 else "single_read_fail")
                    counts["fragments"] += 1
                    counts["alignment_supported"] += int(_pair_positive(bits, pair_policy))
                    counts[label] += 1
                    row = {"ordinal": ordinal, "bin": source.get("bin", "explicit"),
                           "r1_read_id": read_id(a[0]).decode("utf-8", "backslashreplace"),
                           "r2_read_id": read_id(b[0]).decode("utf-8", "backslashreplace") if b else "",
                           "pair_pass": int(_pair_positive(bits, pair_policy)), "evidence": label}
                    for mate in (1, 2) if b else (1,):
                        alignment = next(alignments, None)
                        if alignment is None or alignment["ordinal"] != ordinal or alignment["mate"] != mate:
                            raise RuntimeError("Alignment evidence does not match candidate fragment/mate order")
                        for key in attributes:
                            row[f"r{mate}_{key}"] = alignment[key]
                    writer.writerow(row)
                per_bin.append({"bin": source.get("bin", "explicit"), **dict(counts)})
            if ordinal != expected or calls.read(1) or next(alignments, None) is not None:
                raise RuntimeError("Candidate input/call/evidence counts differ")
    finally:
        evidence.close()
    Path(evidence_database).unlink()
    return per_bin


def _reserve_output(explicit, cfg, sample_id, scope):
    return (Path(explicit).expanduser().absolute() if explicit else
            Path(cfg["work_dir"]) / "qc" / "human" / _safe_sample_id(sample_id) / scope)


def _summary_rows(report):
    screen = report.get("screen") or {}
    total = int(report["input"]["fragments"])
    rows = [("input_fragments", total, total),
            ("human_compatible_fragments", int(screen.get("host_positive_fragments", 0)), total),
            ("any_mate_pass", int(screen.get("any_mate_positive", 0)), total),
            ("neither_mate_pass", int(screen.get("neither_mate_positive", total)), total)]
    if report["input"]["paired"]:
        rows.extend([("both_mates_pass", int(screen.get("both_mates_positive", 0)), total),
                     ("exactly_one_mate_pass", int(screen.get("exactly_one_mate_positive", 0)), total),
                     ("both_mates_concordant_pass", int(screen.get("both_mates_concordant_positive", 0)), total)])
    agreement = report.get("kraken_agreement")
    if agreement:
        rows.extend([
            ("kraken_human_fragments", agreement["kraken_human_fragments"], total),
            ("kraken_human_alignment_supported", agreement["kraken_human_alignment_supported"], agreement["kraken_human_fragments"]),
            ("aligned_kraken_unclassified", agreement["alignment_positive_kraken_unclassified"], total),
            ("aligned_kraken_higher_rank", agreement["alignment_positive_kraken_higher_rank"], total),
            ("aligned_kraken_other_species", agreement["alignment_positive_kraken_other_species"], total)])
    return [{"metric": key, "fragments": count, "denominator_fragments": denom,
             "fraction": count / denom if denom else ""} for key, count, denom in rows]


@managed_outputs
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replace", type=replace_bool, default=True)
    parser.add_argument("--config")
    parser.add_argument("--root")
    parser.add_argument("--sample", help="DATASET/SAMPLE under the configured split output root")
    parser.add_argument("--scope", choices=("candidates", "full"), default="full",
                        help="candidates: existing human bins; full: all original reads plus Kraken comparison [full]")
    parser.add_argument("--output-dir", help="Split sample directory; --scope candidates does not need original raw FASTQs")
    parser.add_argument("--target-fastq", nargs="+", metavar="FASTQ",
                        help="Explicit input: one single-end FASTQ or R1 R2; process every selected record")
    parser.add_argument("--reference", help="Indexed FASTA; default: built by 1.build-lib.sh human")
    parser.add_argument("--bwa-mem2", help="Default: pinned copy installed by 1.build-lib.sh human")
    parser.add_argument("--threads", type=int, help="Aligner threads [min(24, configured classify_cores)]")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-aligned-fraction", type=float, default=0.80)
    parser.add_argument("--min-identity", type=float, default=0.95)
    parser.add_argument("--pair-policy", choices=("any_mate", "both_mates"), default="any_mate")
    parser.add_argument("--contig-regex", help="Human contig regex for a competitive reference, e.g. '^hs_'")
    parser.add_argument("--assignments", help="Existing read-assignments.tsv.gz or species.tar; full scope only")
    parser.add_argument("--summary", help="Existing summary.json; candidate scope uses only sample-level denominator")
    parser.add_argument("--out", help="Validation directory; default is work_dir/qc/human/DATASET/SAMPLE/SCOPE")
    parser.add_argument("--allow-local", action="store_true")
    args = parser.parse_args(argv)
    if args.target_fastq and len(args.target_fastq) not in (1, 2):
        parser.error("--target-fastq accepts one single-end file or two files in R1 R2 order")
    if args.scope == "candidates" and args.assignments:
        parser.error("Candidates use a subset of split reads; --assignments requires --scope full")
    require_compute(args.allow_local)
    cfg = load_config(args.config, args.root)
    args.threads = args.threads if args.threads is not None else min(24, int(cfg.get("classify_cores", 24)))
    if args.threads < 1:
        parser.error("--threads must be positive")
    # This independent QC job has its own scheduler allocation; the classification
    # job's thread setting is not a ceiling on an explicitly requested QC job.
    cfg["classify_cores"] = max(args.threads, int(cfg.get("classify_cores", args.threads)))
    if args.sample:
        if args.output_dir:
            parser.error("Use either --sample or --output-dir")
        args.output_dir = str(configured_split_root(cfg) / _safe_sample_id(args.sample))
    output = Path(args.output_dir).resolve(strict=True) if args.output_dir else None
    saved = None
    if output is not None:
        default_summary = metadata_dir(output) / "summary.json"
        if args.summary is None and default_summary.is_file():
            args.summary = str(default_summary)
    if args.summary:
        saved = json.loads(Path(args.summary).read_text(encoding="utf-8"))
        if args.sample and saved.get("sample_id") not in (None, args.sample):
            parser.error(f"Requested sample {args.sample} differs from split summary {saved['sample_id']}")
    sample_id = args.sample or (saved or {}).get("sample_id") or (
        str(Path(output.parent.name) / output.name) if output else "explicit-input")
    sample_id = str(_safe_sample_id(sample_id))
    sources = []
    if args.target_fastq:
        sources = [{"r1": str(Path(args.target_fastq[0]).resolve(strict=True)),
                    "r2": str(Path(args.target_fastq[1]).resolve(strict=True)) if len(args.target_fastq) == 2 else None,
                    "bin": "explicit"}]
    elif args.scope == "candidates":
        if output is None:
            parser.error("--scope candidates requires --sample, --output-dir or --target-fastq")
        sources = candidate_sources(output)
    elif output is not None:
        if not (metadata_dir(output) / "complete.json").is_file():
            parser.error("Full scope requires a committed split (qc/complete.json)")
        original = json.loads((metadata_dir(output) / "provenance.json").read_text())["context"]["sample"]
        sources = [{"r1": str(Path(original["r1"]).resolve(strict=True)),
                    "r2": str(Path(original["r2"]).resolve(strict=True)) if original.get("r2") else None,
                    "bin": "original_input"}]
    else:
        parser.error("Provide --target-fastq FASTQ [R2], --sample or --output-dir")
    if args.scope == "full" and output is not None:
        if not (metadata_dir(output) / "complete.json").is_file():
            parser.error("Full scope comparison requires a committed split (qc/complete.json)")
        if saved is None:
            parser.error("Committed split is missing summary.json; full comparison cannot verify the saved input identity")
        args.assignments = args.assignments or str(find_assignment_source(output))
    paired = bool(sources[0].get("r2")) if sources else bool((saved or {}).get("paired", True))
    if args.pair_policy == "both_mates" and not paired:
        parser.error("both_mates requires paired FASTQs")
    if any(source.get("r2") and Path(source["r1"]).samefile(source["r2"]) for source in sources):
        parser.error("Paired inputs must be different FASTQ files")
    out = _reserve_output(args.out, cfg, sample_id, args.scope)
    protected = [output, args.reference, args.assignments, args.summary]
    protected += [source.get(mate) for source in sources for mate in ("r1", "r2")]
    if not prepare_result(out, replace=args.replace, cfg=cfg, protected=protected):
        return
    args.reference = args.reference or str(Path(cfg["human_reference_dir"]) / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
    args.bwa_mem2 = args.bwa_mem2 or str(Path(cfg["tools_dir"]) / "bwa-mem2/2.2.1/bwa-mem2-launch")
    cfg.update(human_host_enabled=True,
               human_host_reference=str(Path(args.reference).expanduser().absolute()),
               human_host_executable=args.bwa_mem2,
               human_host_threads=args.threads,
               human_host_min_mapq=args.min_mapq,
               human_host_min_aligned_fraction=args.min_aligned_fraction,
               human_host_min_identity=args.min_identity,
               human_host_pair_policy=args.pair_policy,
               human_host_rescue_policy="unresolved_vertebrate",
               human_host_contig_regex=args.contig_regex)
    context = host_context(cfg) if args.scope == "full" else None
    digest = None
    if args.scope == "candidates":
        expected, digest = count_sources(sources)
        known_human = (saved or {}).get("final_human_fragments", (saved or {}).get("groups", {}).get("human"))
        if not args.target_fastq and known_human is not None and expected != int(known_human):
            raise RuntimeError(f"Human candidate bins contain {expected} fragments but split summary records {known_human}; check incomplete/archived output")
        if not sources and known_human is None:
            raise RuntimeError("No readable human FASTQ bins were found; provide retained candidate FASTQs with --target-fastq")
    elif saved:
        if "paired" in saved and bool(saved["paired"]) != paired:
            parser.error("Target FASTQ single/paired mode differs from the split summary")
        expected = int(saved["fragments"])
        digest = saved.get("input_digest")
        if not isinstance(digest, str) or len(digest) != 64:
            parser.error("Full comparison requires the saved input_digest; a count alone does not establish input identity")
        if saved.get("requested_limit", 0) != 0 or saved.get("validated_scope", "entire_input") != "entire_input":
            parser.error("Full comparison requires an entire-input split summary, not a prefix pilot")
    elif args.assignments:
        expected = count_assignments(args.assignments)
    else:
        expected, digest = count_sources(sources)
    if expected < 0 or (args.scope == "full" and expected == 0):
        parser.error("Full input must contain at least one fragment")
    context = context or (host_context(cfg) if expected else {"meaning": "no human candidates; alignment was not run"})
    screen = None
    if expected:
        sample = {"sample_id": sample_id, "sample": sample_id, **sources[0]}
        screen = screen_human_fastqs(cfg, sample, out, expected, expected_digest=digest,
                                    input_sources=sources,
                                    evidence_database=out / ".human-alignment.sqlite" if args.scope == "candidates" else None)
    report = {
        "schema": "metag-human-validation-v1" if args.scope == "full" else "metag-human-validation-v2",
        "status": "complete" if expected else "no_candidates",
        "sample_id": sample_id, "scope": args.scope, "screen": screen,
        "input": {"r1": sources[0]["r1"] if len(sources) == 1 else None,
                  "r2": sources[0].get("r2") if len(sources) == 1 else None,
                  "sources": sources, "paired": paired, "fragments": expected,
                  "scope": "all_records_in_selected_human_bins" if args.scope == "candidates" else "all_records",
                  "selection": "explicit --target-fastq" if args.target_fastq else
                               "existing human split bins" if args.scope == "candidates" else "split provenance",
                  "raw_reads_retained_in_validation_output": False},
        "sample_fragments": int(saved["fragments"]) if saved else None,
        "reference_interpretation": context["meaning"],
        "interpretation": ("Candidate-only alignment support; missed human-compatible reads and sample-wide human fraction cannot be estimated from this subset."
                           if args.scope == "candidates" else
                           "Full-input human-compatible alignment screen; human-only reference is not a species-specific ground truth."),
        "output_directory": str(out),
    }
    if args.scope == "candidates" and expected:
        report["candidate_bins"] = write_candidate_evidence(sources, screen["calls"], expected, args.pair_policy,
                                                            out / "human-read-evidence.tsv", out / ".human-alignment.sqlite")
    if args.assignments:
        report["kraken_agreement"] = compare_assignments(args.assignments, screen["calls"], expected, args.pair_policy,
                                                        sources[0]["r1"], sources[0].get("r2"))
        with (out / "human-validation-crosstab.tsv").open("w", newline="", encoding="utf-8") as handle:
            rows = report["kraken_agreement"]["crosstab"]
            writer = csv.DictWriter(handle, fieldnames=("alignment", "kraken", "fragments"),
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    with (out / "human-validation.tsv").open("w", newline="", encoding="utf-8") as handle:
        rows = _summary_rows(report)
        writer = csv.DictWriter(handle, fieldnames=("metric", "fragments", "denominator_fragments", "fraction"),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    write_validation_readme(out, sources[0]["r1"] if sources else None, sources[0].get("r2") if sources else None,
                            args.reference, args.assignments, args.summary, scope=args.scope, sources=sources)
    atomic_json(out / "human-validation.json", report)
    finish_result(out)
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
