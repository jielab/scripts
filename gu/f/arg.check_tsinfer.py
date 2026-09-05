#!/usr/bin/env python3
"""Fast structural validation for a completed tsinfer ARG and its sidecars."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import tskit

from csv_compat import enable_wide_csv_fields


enable_wide_csv_fields()


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: ARG post-build check failed: {message}")


def resolved(value: str | Path) -> str:
    return str(Path(value).resolve())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", required=True, type=Path)
    parser.add_argument("--sample-map", required=True, type=Path)
    parser.add_argument("--qc-json", required=True, type=Path)
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--chr", required=True)
    args = parser.parse_args()

    try:
        ts = tskit.load(str(args.trees))
    except Exception as exc:
        fail(f"unreadable tree sequence {args.trees}: {exc}")
    if ts.num_trees < 1 or ts.num_sites < 1 or ts.num_samples < 1:
        fail(
            f"empty tree sequence: trees={ts.num_trees} sites={ts.num_sites} "
            f"samples={ts.num_samples}"
        )
    if ts.sequence_length <= 0:
        fail(f"invalid sequence length: {ts.sequence_length}")

    try:
        qc = json.loads(args.qc_json.read_text(encoding="utf-8"))
    except Exception as exc:
        fail(f"unreadable QC JSON {args.qc_json}: {exc}")
    expected_paths = {
        "input_vcf": resolved(args.input_vcf),
        "output_trees": resolved(args.trees),
        "sample_map": resolved(args.sample_map),
    }
    for key, expected in expected_paths.items():
        if resolved(qc.get(key, "")) != expected:
            fail(f"QC {key} does not match {expected}")
    expected_counts = {
        "sample_nodes": ts.num_samples,
        "sites": ts.num_sites,
        "mutations": ts.num_mutations,
        "nodes": ts.num_nodes,
        "edges": ts.num_edges,
        "trees": ts.num_trees,
    }
    for key, expected in expected_counts.items():
        if qc.get(key) != expected:
            fail(f"QC {key}={qc.get(key)!r}, tree sequence has {expected}")
    if float(qc.get("sequence_length", -1)) != float(ts.sequence_length):
        fail("QC sequence_length does not match tree sequence")
    if qc.get("ancestral_state_source") != "INFO/AA (VCF-Zarr variant_AA)":
        fail("QC does not record INFO/AA as the ancestral-state source")

    try:
        with args.sample_map.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["tree_node_id", "sample", "haplotype"]:
                fail(f"unexpected sample-map header: {reader.fieldnames}")
            rows = list(reader)
    except Exception as exc:
        fail(f"unreadable sample map {args.sample_map}: {exc}")
    if len(rows) != ts.num_samples:
        fail(f"sample-map rows={len(rows)}, expected {ts.num_samples}")
    try:
        mapped_nodes = [int(row["tree_node_id"]) for row in rows]
        if any(not row["sample"] or int(row["haplotype"]) < 1 for row in rows):
            fail("sample map contains an empty sample or invalid haplotype")
    except (TypeError, ValueError) as exc:
        fail(f"invalid sample-map value: {exc}")
    if len(mapped_nodes) != len(set(mapped_nodes)):
        fail("sample map contains duplicate tree node IDs")
    if set(mapped_nodes) != set(map(int, ts.samples())):
        fail("sample-map node IDs do not equal tree-sequence sample nodes")
    if args.chr == 'X' and 'sample_ploidies' not in qc:
        fail('X QC lacks biological sample ploidies; rebuild with the updated X workflow')
    if 'sample_ploidies' in qc:
        observed = {}
        for row in rows:
            observed.setdefault(row['sample'], []).append(int(row['haplotype']))
        expected = {name: list(range(1, int(ploidy) + 1)) for name, ploidy in qc['sample_ploidies'].items()}
        if {name: sorted(haps) for name, haps in observed.items()} != expected:
            fail('Sample map differs from biological sample ploidies')

    print(
        f"[ARG CHECK] PASS chr={args.chr} trees={ts.num_trees} "
        f"sites={ts.num_sites} samples={ts.num_samples} sequence_length={ts.sequence_length:g}"
    )


if __name__ == "__main__":
    main()
