#!/usr/bin/env python3
"""Export a native ARG tree and sample map to GU TRACE's input contract."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import tskit


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, path)


def read_map(path: Path, method: str) -> list[tuple[int, str, int]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        rows = []
        if method == "needle":
            required = {"eid", "hap", "sample_node"}
            if not required <= fields:
                raise SystemExit(
                    f"Needle sample map requires {sorted(required)}: {path}"
                )
            for row in reader:
                rows.append(
                    (int(row["sample_node"]), str(row["eid"]), int(row["hap"]) + 1)
                )
        else:
            required = {"tree_node_id", "sample", "haplotype"}
            if not required <= fields:
                raise SystemExit(
                    f"tsinfer sample map requires {sorted(required)}: {path}"
                )
            for row in reader:
                rows.append(
                    (
                        int(row["tree_node_id"]),
                        str(row["sample"]),
                        int(row["haplotype"]),
                    )
                )
    return rows


def validate(ts: tskit.TreeSequence, rows: list[tuple[int, str, int]], chrom: str) -> dict:
    if min(ts.num_trees, ts.num_samples, ts.num_sites) < 1:
        raise SystemExit(f"Empty tree sequence for chr{chrom}")
    nodes = list(map(int, ts.samples()))
    mapped = [row[0] for row in rows]
    if len(rows) != len(nodes) or set(mapped) != set(nodes) or len(set(mapped)) != len(mapped):
        raise SystemExit(
            f"TRACE sample map does not bijectively cover chr{chrom} sample nodes: "
            f"rows={len(rows)} samples={len(nodes)}"
        )
    positions = [site.position for site in ts.sites()]
    span = (max(positions) - min(positions)) / max(1.0, float(ts.sequence_length))
    if any(hap < 1 for _, _, hap in rows):
        raise SystemExit(f"TRACE haplotype numbers must be one-based for chr{chrom}")
    return {
        "chromosome": str(chrom),
        "trees": ts.num_trees,
        "nodes": ts.num_nodes,
        "samples": ts.num_samples,
        "sites": ts.num_sites,
        "mutations": ts.num_mutations,
        "sequence_length": ts.sequence_length,
        "site_span_fraction": span,
    }


def write_tree(source: Path, dest: Path, ts: tskit.TreeSequence, method: str) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f".{dest.name}.tmp.{os.getpid()}")
    tmp.unlink(missing_ok=True)
    offset = float(ts.metadata.get('offset', 0)) if method == 'needle' and isinstance(ts.metadata, dict) else 0
    if offset:
        # Native ARG-Needle coordinates are relative to the first input SNP.
        # GU loci and TRACE maps are chromosome coordinates, so translate once.
        tables = ts.dump_tables()
        tables.sequence_length += offset
        tables.edges.left += offset
        tables.edges.right += offset
        tables.sites.position += offset
        if len(tables.migrations):
            tables.migrations.left += offset
            tables.migrations.right += offset
        tables.metadata_schema = tskit.MetadataSchema.permissive_json()
        tables.metadata = {**ts.metadata, 'offset': 0, 'source_arg_offset': offset,
                           'coordinate_system': 'chromosome_bp'}
        tables.time_units = 'generations'
        tables.tree_sequence().dump(tmp)
        os.replace(tmp, dest)
        return 'chromosome_coordinate_copy'
    if ts.time_units == "generations":
        os.symlink(os.path.relpath(source.resolve(), dest.parent.resolve()), tmp)
        os.replace(tmp, dest)
        return "relative_symlink"
    if method != "needle" or ts.time_units not in ("unknown", ""):
        raise SystemExit(
            f"TRACE requires node times in generations; got {ts.time_units!r}: {source}"
        )
    tables = ts.dump_tables()
    tables.time_units = "generations"
    tables.provenances.add_row(
        record=json.dumps(
            {
                "schema_version": "1.0.0",
                "software": {"name": "arg.trace", "version": "1"},
                "parameters": {
                    "source_method": "needle",
                    "operation": "declare_ARG_Needle_normalized_node_times_as_generations",
                },
            },
            separators=(",", ":"),
        )
    )
    tables.tree_sequence().dump(tmp)
    os.replace(tmp, dest)
    return "converted_copy"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True, choices=("needle", "tsinfer"))
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--sample-map", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--chr", required=True)
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    for path in (args.tree, args.sample_map):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing ARG input: {path}")
    out_tree = args.out_dir / f"chr{args.chr}.trees"
    out_map = args.out_dir / f"chr{args.chr}.sample_map.tsv"
    out_qc = args.out_dir / f"chr{args.chr}.arg_qc.json"
    out_meta = args.out_dir / f"chr{args.chr}.trace.meta.json"
    source_signature = {
        "tree": [str(args.tree.resolve()), args.tree.stat().st_size, args.tree.stat().st_mtime_ns],
        "sample_map": [
            str(args.sample_map.resolve()),
            args.sample_map.stat().st_size,
            args.sample_map.stat().st_mtime_ns,
        ],
        "method": args.method,
        "format": "trace",
    }
    outputs = (out_tree, out_map, out_qc, out_meta)
    if args.check:
        for path in outputs:
            if not path.exists() or path.stat().st_size == 0:
                raise SystemExit(f"Missing TRACE-format ARG output: {path}")
        old = json.loads(out_meta.read_text())
        if old.get("source") != source_signature:
            raise SystemExit(
                f"TRACE-format ARG source provenance differs for chr{args.chr}; "
                "rerun build with --format trace"
            )
        out_ts = tskit.load(out_tree)
        rows = read_map(out_map, "tsinfer")
        stats = validate(out_ts, rows, args.chr)
        if out_ts.time_units != "generations":
            raise SystemExit(f"TRACE tree time units are not generations: {out_tree}")
        print(
            f"TRACE format CHECK PASS method={args.method} chr{args.chr} "
            f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
        )
        return

    if not args.replace and all(path.exists() and path.stat().st_size > 0 for path in outputs):
        try:
            old = json.loads(out_meta.read_text())
            if old.get("source") == source_signature:
                out_ts = tskit.load(out_tree)
                rows = read_map(out_map, "tsinfer")
                stats = validate(out_ts, rows, args.chr)
                if out_ts.time_units == "generations":
                    print(
                        f"TRACE format SKIP method={args.method} chr{args.chr} "
                        f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
                    )
                    return
        except (OSError, ValueError, json.JSONDecodeError, tskit.FileFormatError):
            pass

    ts = tskit.load(args.tree)
    rows = read_map(args.sample_map, args.method)
    stats = validate(ts, rows, args.chr)
    storage = write_tree(args.tree, out_tree, ts, args.method)
    map_text = "tree_node_id\tsample\thaplotype\n" + "".join(
        f"{node}\t{sample}\t{hap}\n" for node, sample, hap in rows
    )
    atomic_text(out_map, map_text)
    out_ts = tskit.load(out_tree)
    stats = validate(out_ts, read_map(out_map, "tsinfer"), args.chr)
    if out_ts.time_units != "generations":
        raise SystemExit(f"TRACE tree time units are not generations: {out_tree}")
    qc = {
        **stats,
        "method": args.method,
        "format": "trace",
        "time_units": out_ts.time_units,
        "storage": storage,
    }
    atomic_text(out_qc, json.dumps(qc, indent=2, sort_keys=True) + "\n")
    atomic_text(
        out_meta,
        json.dumps({"source": source_signature, "output": qc}, indent=2, sort_keys=True)
        + "\n",
    )
    print(
        f"TRACE format PASS method={args.method} chr{args.chr} storage={storage} "
        f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
    )


if __name__ == "__main__":
    main()
