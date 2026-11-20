#!/usr/bin/env python3
"""Matched, seeded threshold sensitivity on a random sample of complete raw reads.

The unchanged full UNION is queried three times, sequentially. This does not
change production labels, build an index, materialize species FASTQs, or claim
that newly classified reads are correct. One full raw scan validates the saved
split digest and preserves original paired records in the sampled input.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import gzip
import itertools
import json
from pathlib import Path
import random
import struct
import subprocess
import sys
import time

from common import atomic_json, engine_snapshot, load_config, require_compute, utcnow
from fastq_stream import fragments, new_digest, read_id, update_digest
from ordinal_query import input_states, query_once, validate_snapshot
from project_cache import load_taxonomy
from qc_paths import configured_split_root
from result_output import managed_outputs, prepare, finish, replace_bool
from storage_layout import metadata_dir

SCHEMA = "metag-threshold-diagnostic-v1"
PROFILES = (("c10g2", 0.1, 2), ("c0g2", 0.0, 2), ("c0g1", 0.0, 1))
STATUSES = ("assigned_species", "higher_rank", "unclassified")


def sample_parts(value):
    parts = value.split("/")
    if len(parts) != 2 or any(not x or x in (".", "..") or "\\" in x for x in parts):
        raise ValueError("--sample must be DATASET/SAMPLE")
    return parts


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def extract_sample(sample, expected, expected_digest, requested, seed, out):
    """Uniform sampling without replacement over complete fragment ordinals."""
    positive(expected, "Saved fragment count")
    positive(requested, "--fragments")
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise ValueError("--seed must be a nonnegative integer")
    if not expected_digest:
        raise ValueError("The completed split must record an input_digest")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    selected = sorted(random.Random(seed).sample(range(expected), min(expected, requested)))
    wanted = iter(selected)
    next_selected = next(wanted, None)
    source_states = input_states(sample)
    paired = sample.get("r2") is not None
    paths = [out / "sample.R1.fastq"] + ([out / "sample.R2.fastq"] if paired else [])
    full_digest, selected_digest = new_digest(), new_digest()
    scanned = kept = bases = 0
    started = time.monotonic()
    with ExitStack() as stack:
        handles = [stack.enter_context(path.open("wb", buffering=1 << 20)) for path in paths]
        ordinals = stack.enter_context((out / "sample-ordinals.tsv").open("w", newline="", encoding="utf-8"))
        writer = csv.writer(ordinals, delimiter="\t", lineterminator="\n")
        writer.writerow(("sampled_ordinal", "original_ordinal", "r1_read_id", "r2_read_id"))
        stream = fragments(sample["r1"], sample.get("r2"))
        stack.callback(stream.close)
        for index, (a, b) in enumerate(stream):
            scanned += 1
            if scanned > expected:
                raise RuntimeError("Original input has more fragments than the completed split")
            update_digest(full_digest, a, b)
            if index == next_selected:
                update_digest(selected_digest, a, b)
                for handle, record in zip(handles, (a, b) if paired else (a,)):
                    handle.writelines(record)
                    bases += len(record[1].rstrip(b"\r\n"))
                kept += 1
                writer.writerow((kept, index + 1, read_id(a[0]).decode("utf-8", "backslashreplace"),
                                 read_id(b[0]).decode("utf-8", "backslashreplace") if b else ""))
                next_selected = next(wanted, None)
            if scanned % 1_000_000 == 0:
                print(f"Random-sample scan {scanned:,}/{expected:,}; selected {kept:,}", flush=True)
    if scanned != expected or full_digest.hexdigest() != expected_digest:
        raise RuntimeError("Original FASTQ bytes/order/count differ from the completed split")
    if input_states(sample) != source_states or kept != len(selected) or next_selected is not None:
        raise RuntimeError("Original inputs changed or random sample is incomplete")
    return {"requested_fragments": requested, "fragments": kept, "original_fragments": expected,
            "seed": seed, "selection": "simple random fragment ordinals without replacement",
            "scope": "whole original input; not a prefix", "paired": paired,
            "original_input_digest": full_digest.hexdigest(), "input_digest": selected_digest.hexdigest(),
            "full_input_digest_verified": True, "total_bases": bases,
            "fraction_sampled": kept / expected, "scan_seconds": time.monotonic() - started,
            "r1": str(paths[0]), "r2": str(paths[1]) if paired else None,
            "original_ordinals": [x + 1 for x in selected]}


def worker(path):
    job = json.loads(Path(path).read_text())
    if job.get("schema") != SCHEMA:
        raise ValueError("Invalid threshold worker descriptor")
    require_compute(job.get("allow_local", False))
    cfg = job["config"]
    if cfg.get("database_layout") != "partitioned" or not cfg.get("union_memory_mapping"):
        raise ValueError("Threshold diagnosis requires the full mapped UNION")
    validate_snapshot(job["snapshot"])
    profile = tuple(job["profile"])
    if profile not in PROFILES:
        raise ValueError("Unknown threshold profile")
    cfg = dict(cfg, confidence=profile[1], minimum_hit_groups=profile[2],
               classifier_engine="kraken", classify_workers=1, union_memory_mapping=True)
    out = Path(job["out"])
    out.mkdir(parents=True, exist_ok=True)
    with load_taxonomy(cfg, job["snapshot"], [job["sample"]["dataset"]]) as taxonomy:
        result = query_once(cfg, job["sample"], job["snapshot"], out, taxonomy=taxonomy)
        taxonomy.validate_current()
    if (result["fragments"] != job["expected_fragments"] or
            result["input_digest"] != job["expected_digest"]):
        raise RuntimeError("Threshold worker did not classify the exact matched random sample")
    atomic_json(out / "query.json", result)


def run_profiles(cfg, snapshot, sample, sampled, out, allow_local=False):
    """Wait for every worker before launching the next; retain failed evidence."""
    records = []
    for name, confidence, groups in PROFILES:
        folder = Path(out) / name
        folder.mkdir()
        payload = {"schema": SCHEMA, "config": dict(cfg, union_memory_mapping=True),
                   "snapshot": snapshot, "sample": sample, "profile": [name, confidence, groups],
                   "expected_fragments": sampled["fragments"], "expected_digest": sampled["input_digest"],
                   "out": str(folder), "allow_local": allow_local}
        job = folder / "worker.json"
        atomic_json(job, payload)
        record = {"profile": name, "confidence": confidence, "minimum_hit_groups": groups,
                  "status": "failed", "query": None, "error": ""}
        try:
            with (folder / "worker.log").open("w") as log:
                result = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(job)],
                                        stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError(f"Worker exit {result.returncode}; see {folder / 'worker.log'}")
            query = json.loads((folder / "query.json").read_text())
            if query["fragments"] != sampled["fragments"] or query["input_digest"] != sampled["input_digest"]:
                raise RuntimeError("Worker result differs from the shared random input")
            record.update(status="complete", query=query)
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
        records.append(record)
    return records


def compact_taxids(path, expected):
    path = Path(path)
    if path.stat().st_size != expected * 8:
        raise RuntimeError(f"Incomplete compact taxids: {path}")
    with path.open("rb") as handle:
        return [x[0] for x in struct.iter_unpack("<Q", handle.read())]


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize(records, sampled, taxonomy, out):
    """Conserve the matched sample in summaries and every pairwise transition."""
    summaries, taxids, details = [], {}, {}
    count = sampled["fragments"]
    for record in records:
        row = {k:record[k] for k in ("profile", "confidence", "minimum_hit_groups", "status", "error")}
        if record["status"] == "complete":
            query = record["query"]
            try:
                ids = compact_taxids(query["taxids"], count)
                infos = [taxonomy.describe(tid) for tid in ids]
                status = Counter(x["status"] for x in infos)
                if set(status) - set(STATUSES) or sum(status.values()) != count:
                    raise RuntimeError("Invalid normalized status or count in threshold query")
                taxids[record["profile"]], details[record["profile"]] = ids, infos
                row.update(fragments=count, **{s:status[s] for s in STATUSES},
                    human_fragments=sum(x["group"] == "human" for x in infos),
                    plant_fragments=sum(x["group"] == "plant" for x in infos),
                    plant_species_fragments=sum(x["group"] == "plant" and x["status"] == "assigned_species" for x in infos),
                    plant_higher_rank_fragments=sum(x["group"] == "plant" and x["status"] == "higher_rank" for x in infos),
                    query_seconds=query["seconds"], index_load_seconds=query.get("batch_index_load_seconds"),
                    memory_mode="mmap" if query.get("memory_mapping") else "ram")
            except Exception as exc:
                record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
                row.update(status="failed", error=record["error"])
        summaries.append(row)
    fields = ["profile", "confidence", "minimum_hit_groups", "status", "fragments", *STATUSES,
              "human_fragments", "plant_fragments", "plant_species_fragments", "plant_higher_rank_fragments",
              "query_seconds", "index_load_seconds", "memory_mode", "error"]
    write_tsv(Path(out) / "threshold-summary.tsv", summaries, fields)
    transitions, targets, agreements = [], [], []
    for left, right in itertools.combinations(details, 2):
        status_counts = Counter((a["status"], b["status"]) for a,b in zip(details[left], details[right]))
        for (a,b), n in sorted(status_counts.items()):
            transitions.append(dict(left_profile=left, right_profile=right, left_status=a,
                                    right_status=b, fragments=n, fraction=n/count))
        for target in ("human", "plant"):
            target_counts = Counter((a["group"] == target, b["group"] == target)
                                    for a,b in zip(details[left], details[right]))
            for (a,b), n in sorted(target_counts.items()):
                targets.append(dict(left_profile=left, right_profile=right, target=target,
                                    left_target=int(a), right_target=int(b), fragments=n, fraction=n/count))
        same = sum(a == b for a,b in zip(taxids[left], taxids[right]))
        both_species = [(a,b) for a,b in zip(details[left], details[right])
                        if a["species_taxid"] and b["species_taxid"]]
        species_same = sum(a["species_taxid"] == b["species_taxid"] for a,b in both_species)
        agreements.append(dict(left_profile=left, right_profile=right, fragments=count,
                               same_taxid=same, changed_taxid=count-same, both_species=len(both_species),
                               same_species=species_same, different_species=len(both_species)-species_same))
    write_tsv(Path(out) / "threshold-transitions.tsv", transitions,
              ["left_profile", "right_profile", "left_status", "right_status", "fragments", "fraction"])
    write_tsv(Path(out) / "threshold-targets.tsv", targets,
              ["left_profile", "right_profile", "target", "left_target", "right_target", "fragments", "fraction"])
    write_tsv(Path(out) / "threshold-agreement.tsv", agreements,
              ["left_profile", "right_profile", "fragments", "same_taxid", "changed_taxid", "both_species", "same_species", "different_species"])
    with gzip.open(Path(out) / "threshold-taxids.tsv.gz", "wt", encoding="utf-8", newline="", compresslevel=3) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sampled_ordinal", "original_ordinal", *taxids])
        for i, original in enumerate(sampled["original_ordinals"]):
            writer.writerow([i+1, original, *[values[i] for values in taxids.values()]])
    return summaries


@managed_outputs
def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "--worker":
        if len(argv) != 2:
            raise ValueError("--worker needs exactly one descriptor")
        return worker(argv[1])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out")
    parser.add_argument("--fragments", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20261114)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--allow-local", action="store_true")
    parser.add_argument("--replace", type=replace_bool, default=True)
    parser.add_argument("--keep-sample", action="store_true")
    args = parser.parse_args(argv)
    require_compute(args.allow_local)
    parts = sample_parts(args.sample)
    positive(args.fragments, "--fragments")
    cfg = load_config(args.config)
    if cfg.get("database_layout") != "partitioned":
        raise ValueError("Threshold diagnosis uses the existing partitioned full UNION configuration")
    cfg = dict(cfg, classifier_engine="kraken", classify_workers=1, union_memory_mapping=True,
               min_free_gib=min(float(cfg.get("min_free_gib", 1)), 1.0))
    if args.threads is not None:
        cfg["classify_threads"] = positive(args.threads, "--threads")
    source = configured_split_root(cfg).joinpath(*parts)
    metadata = metadata_dir(source)
    receipt = json.loads((metadata / "complete.json").read_text())
    summary = json.loads((metadata / "summary.json").read_text())
    provenance = json.loads((metadata / "provenance.json").read_text())
    if receipt.get("status") != "complete" or summary.get("sample_id") != args.sample:
        raise ValueError("Requested sample must have a matching completed split")
    if summary.get("validated_scope") != "entire_input" or summary.get("requested_limit", 0):
        raise ValueError("Random whole-input diagnosis needs a completed full-input split")
    sample = dict(provenance["context"]["sample"])
    if sample.get("sample_id") != args.sample or bool(sample.get("r2")) != bool(summary.get("paired")):
        raise ValueError("Sample identity/paired layout differs between split metadata")
    if sample.get("r2") and Path(sample["r1"]).samefile(sample["r2"]):
        raise ValueError("Paired inputs must be different files")
    # Hold base-quality filtering constant while varying only confidence and
    # hit groups. Existing production thresholds are recorded for comparison.
    cfg["minimum_base_quality"] = int(provenance["context"].get("minimum_base_quality", cfg.get("minimum_base_quality", 0)))
    snapshot = engine_snapshot(cfg)  # ordinary cheap library check, never --full
    original_database = provenance["context"].get("database_identity")
    if not original_database or snapshot["identity"] != original_database:
        raise ValueError("Current UNION differs from the completed split; threshold-only comparison needs the same library")
    out = Path(args.out).expanduser().absolute() if args.out else Path(cfg["work_dir"]).joinpath("qc", "threshold", *parts).absolute()
    request = dict(sample_id=args.sample, seed=args.seed, requested_fragments=args.fragments,
                   original_fragments=summary["fragments"], original_input_digest=summary["input_digest"],
                   database_identity=snapshot["identity"], profiles=[list(p) for p in PROFILES],
                   minimum_base_quality=cfg["minimum_base_quality"])
    prior = out / "threshold-diagnostic.json"
    if not args.replace and prior.is_file() and json.loads(prior.read_text()).get("request") != request:
        raise ValueError("Existing threshold diagnostic uses different inputs/settings; use --replace TRUE or another --out")
    if not prepare(out, replace=args.replace, cfg=cfg, protected=[source, sample["r1"], sample.get("r2")]):
        return
    sampled = extract_sample(sample, int(summary["fragments"]), summary["input_digest"], args.fragments, args.seed, out)
    selected = dict(sample, r1=sampled["r1"], r2=sampled["r2"])
    records = run_profiles(cfg, snapshot, selected, sampled, out, args.allow_local)
    with load_taxonomy(cfg, snapshot, [parts[0]]) as taxonomy:
        profiles = summarize(records, sampled, taxonomy, out)
        taxonomy.validate_current()
    validate_snapshot(snapshot)
    successful = all(record["status"] == "complete" for record in records)
    report = {"schema": SCHEMA, "created_utc": utcnow(), "request": request,
              "status": "complete" if successful else "incomplete", "source_split": str(source),
              "original_production_parameters": {k:provenance["context"].get(k) for k in
                    ("confidence", "minimum_hit_groups", "minimum_base_quality")},
              "sampling": {k:v for k,v in sampled.items() if k != "original_ordinals"},
              "profiles": profiles, "workers": records, "sample_retained": args.keep_sample or not successful,
              "interpretation": [
                  "All profiles query the same random paired fragments against the same complete UNION.",
                  "c10g2 is the configured production threshold convention; c0g2 and c0g1 diagnose confidence/hit-group sensitivity.",
                  "A lower threshold causing more calls does not establish correct species, human or plant identification.",
                  "Persistently unclassified reads can reflect absent references, quality, divergence or uninformative sequence; this test cannot separate those causes.",
                  "Rare human/plant candidates may be absent from this random sample; validate all saved candidates independently.",
                  "Mapping avoids a private full-index preload but does not impose a hard RSS cap; inspect each query.log and native metrics."]}
    if successful and not args.keep_sample:
        for key in ("r1", "r2"):
            if sampled.get(key):
                Path(sampled[key]).unlink()
    atomic_json(out / "threshold-diagnostic.json", report)
    if not successful:
        raise RuntimeError("Threshold diagnosis incomplete; other profiles were still attempted. See threshold-diagnostic.json and retained worker logs")
    finish(out)
    print(f"Threshold sensitivity complete: {out / 'threshold-summary.tsv'}", flush=True)


if __name__ == "__main__":
    main()
