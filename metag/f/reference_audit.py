#!/usr/bin/env python3
"""Audit the partitioned metagenomics reference library, with a human-specific report.

This is intentionally metadata/receipt-only: it does not reread terabytes of FASTA
or hash files. It proves which frozen reference paths/taxids were actually planned
and how many prepared sequence bases were assigned to them.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path

from common import atomic_json, load_config
from storage_layout import library_root
from taxonomy import Taxonomy


def _read_json(path):
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected JSON object: {path}")
    return value


def audit(cfg):
    base = Path(library_root(cfg)).resolve()
    plan = _read_json(base / "plan.json")
    tasks = plan.get("tasks") or plan.get("context", {}).get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise RuntimeError(f"No partitioned reference tasks in {base / 'plan.json'}")
    taxonomy = Taxonomy(Path(cfg["metadata_dir"]) / "taxonomy")

    frozen_map = {}
    map_path = Path(cfg["metadata_dir"]) / "reference_taxids.tsv"
    with map_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"path", "taxid"} <= set(reader.fieldnames):
            raise RuntimeError(f"Invalid reference map: {map_path}")
        for row in reader:
            if row["path"] in frozen_map:
                raise RuntimeError(f"Duplicate path in reference map: {row['path']}")
            frozen_map[row["path"]] = int(row["taxid"])

    records = []
    seen_paths = set()
    for task in tasks:
        task_id = int(task["task_id"])
        manifest_path = base / task["manifest"]["path"]
        manifest = _read_json(manifest_path)
        rows = manifest.get("files")
        if not isinstance(rows, list) or not rows:
            raise RuntimeError(f"Empty/invalid manifest: {manifest_path}")
        receipt_path = base / "prepared" / f"task{task_id:06d}" / "receipt.json"
        receipt = _read_json(receipt_path)
        sources = receipt.get("context", {}).get("stats", {}).get("sources")
        if not isinstance(sources, list) or len(sources) != len(rows):
            raise RuntimeError(f"Prepared source counters do not match manifest: {receipt_path}")
        for row, source in zip(rows, sources):
            relative = row["path"]
            if source.get("path") != relative:
                raise RuntimeError(f"Prepared source order/path mismatch: {receipt_path}: {relative}")
            taxid = int(row["taxid"])
            if frozen_map.get(relative) != taxid:
                raise RuntimeError(f"Frozen metadata taxid differs from task manifest: {relative}")
            if relative in seen_paths:
                raise RuntimeError(f"Reference path appears in multiple tasks: {relative}")
            seen_paths.add(relative)
            info = taxonomy.describe(taxid)
            records.append({
                "path": relative,
                "accession": row.get("accession") or "",
                "source_taxid": taxid,
                "resolved_taxid": info["taxid"],
                "species_taxid": info["species_taxid"],
                "name": info["name"],
                "group": info["group"],
                "category": row.get("category") or "",
                "records": int(source["records"]),
                "bases": int(source["bases"]),
                "compressed_bytes": int(row.get("state", {}).get("size", 0)),
                "task_id": task_id,
            })
    if seen_paths != set(frozen_map):
        missing = sorted(set(frozen_map) - seen_paths)
        extra = sorted(seen_paths - set(frozen_map))
        raise RuntimeError(f"Plan/reference map coverage differs: missing={missing[:5]} extra={extra[:5]}")

    by_group = defaultdict(lambda: {"references": 0, "records": 0, "bases": 0, "compressed_bytes": 0})
    by_category = defaultdict(lambda: {"references": 0, "records": 0, "bases": 0, "compressed_bytes": 0})
    by_taxid = defaultdict(lambda: {"references": 0, "records": 0, "bases": 0, "compressed_bytes": 0})
    for row in records:
        for table, key in ((by_group, row["group"]), (by_category, row["category"]), (by_taxid, row["source_taxid"])):
            item = table[key]
            item["references"] += 1
            item["records"] += row["records"]
            item["bases"] += row["bases"]
            item["compressed_bytes"] += row["compressed_bytes"]

    humans = [row for row in records if row["group"] == "human"]
    exact_9606 = [row for row in records if row["source_taxid"] == 9606]
    summary = {
        "schema": "metag-reference-audit-v1",
        "library": str(base),
        "reference_map": str(map_path),
        "tasks": len(tasks),
        "references": len(records),
        "records": sum(row["records"] for row in records),
        "bases": sum(row["bases"] for row in records),
        "groups": dict(sorted(by_group.items())),
        "categories": dict(sorted(by_category.items())),
        "human": {
            "present": bool(humans),
            "reference_files": len(humans),
            "exact_taxid_9606_files": len(exact_9606),
            "records": sum(row["records"] for row in humans),
            "bases": sum(row["bases"] for row in humans),
            "compressed_bytes": sum(row["compressed_bytes"] for row in humans),
            "taxids": dict(Counter(str(row["source_taxid"]) for row in humans)),
            "meaning": "reference rows whose frozen NCBI lineage contains Homo sapiens taxid 9606"
        },
    }
    return summary, records, humans, by_taxid


def enrich_assembly_metadata(records, metadata_dir):
    """Join frozen NCBI assembly-summary quality fields by accession.

    The existing library need not be rebuilt; this only reads the small metadata
    summaries already used during reference preparation.  It allows downstream
    representative-genome selection to prefer reference/representative and
    complete/chromosome assemblies instead of simply preferring the largest file.
    """
    wanted = {row["accession"] for row in records if row.get("accession")}
    if not wanted:
        return records, 0
    folder = Path(metadata_dir) / "summaries"
    if not folder.is_dir():
        return records, 0
    found = {}
    fields_keep = ("refseq_category", "species_taxid", "organism_name", "version_status",
                   "assembly_level", "release_type", "genome_rep", "excluded_from_refseq")
    for path in sorted(folder.glob("*assembly_summary*")):
        if not path.is_file() or path.name.endswith((".aria2", ".partial", ".download")):
            continue
        opener = gzip.open if path.name.endswith(".gz") else open
        header = None
        try:
            with opener(path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("#"):
                        values = line.lstrip("#").strip().split("\t")
                        if values and values[0] == "assembly_accession":
                            header = values
                        continue
                    if not line.strip() or header is None:
                        continue
                    cols = line.rstrip("\r\n").split("\t")
                    if len(cols) != len(header):
                        continue
                    acc = cols[0]
                    if acc not in wanted:
                        continue
                    row = dict(zip(header, cols))
                    meta = {key: row.get(key, "") for key in fields_keep}
                    # Prefer a current record over a historical duplicate.
                    previous = found.get(acc)
                    if previous is None or (previous.get("version_status") != "latest" and meta.get("version_status") == "latest"):
                        found[acc] = meta
        except (OSError, EOFError, gzip.BadGzipFile):
            continue
    for row in records:
        row.update(found.get(row.get("accession"), {}))
    return records, len(found)


def choose_species_representatives(records):
    """Choose one existing raw assembly per NCBI species for downstream profiling.

    This does not change the Kraken library.  When frozen NCBI assembly-summary
    metadata are available, reference/representative and complete/chromosome
    assemblies are preferred before namespace and size.  The choice is deterministic.
    """
    grouped = defaultdict(list)
    for row in records:
        species = int(row.get("species_taxid") or 0)
        if species:
            grouped[species].append(row)
    result = []
    for species, rows in grouped.items():
        def key(row):
            accession = str(row.get("accession") or "")
            namespace = 0 if accession.startswith("GCF_") else (1 if accession.startswith("GCA_") else 2)
            category = str(row.get("refseq_category") or "").lower()
            category_score = 0 if category == "reference genome" else (1 if category == "representative genome" else 2)
            level = str(row.get("assembly_level") or "").lower()
            level_score = {"complete genome": 0, "chromosome": 1, "scaffold": 2, "contig": 3}.get(level, 4)
            latest_score = 0 if str(row.get("version_status") or "").lower() == "latest" else 1
            excluded_score = 0 if not str(row.get("excluded_from_refseq") or "").strip() else 1
            return (latest_score, category_score, level_score, excluded_score, namespace,
                    -int(row.get("bases", 0)), -int(row.get("records", 0)), str(row.get("path", "")))
        chosen = min(rows, key=key)
        result.append({
            "species_taxid": species, "name": chosen["name"], "group": chosen["group"],
            "path": chosen["path"], "accession": chosen["accession"],
            "source_taxid": chosen["source_taxid"], "bases": chosen["bases"],
            "records": chosen["records"], "compressed_bytes": chosen["compressed_bytes"],
            "candidate_reference_files": len(rows),
            "refseq_category": chosen.get("refseq_category", ""),
            "assembly_level": chosen.get("assembly_level", ""),
            "version_status": chosen.get("version_status", ""),
            "excluded_from_refseq": chosen.get("excluded_from_refseq", ""),
            "selection": "latest_refseq_category_assembly_level_GCF_then_bases",
        })
    result.sort(key=lambda row: (row["group"], row["name"], row["species_taxid"]))
    return result


def _write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config")
    parser.add_argument("--root")
    parser.add_argument("--out")
    parser.add_argument("--require-human", action="store_true")
    parser.add_argument("--skip-assembly-metadata", action="store_true",
                        help="Do not scan frozen NCBI assembly summaries for representative-quality fields")
    args = parser.parse_args(argv)
    cfg = load_config(args.config, args.root)
    summary, records, humans, by_taxid = audit(cfg)
    if not args.skip_assembly_metadata:
        records, matched_assemblies = enrich_assembly_metadata(records, cfg["metadata_dir"])
        summary["assembly_metadata_matches"] = matched_assemblies
    representatives = choose_species_representatives(records)
    summary["species_representatives"] = len(representatives)
    out = Path(args.out) if args.out else Path(cfg["work_dir"]) / "qc" / "reference-audit"
    out.mkdir(parents=True, exist_ok=True)
    atomic_json(out / "reference-audit.json", summary)
    fields = ["path", "accession", "source_taxid", "resolved_taxid", "species_taxid", "name",
              "group", "category", "records", "bases", "compressed_bytes", "task_id"]
    _write_tsv(out / "references.tsv", records, fields)
    _write_tsv(out / "human-references.tsv", humans, fields)
    _write_tsv(out / "species-representatives.tsv", representatives,
               ["species_taxid", "name", "group", "path", "accession", "source_taxid",
                "bases", "records", "compressed_bytes", "candidate_reference_files",
                "refseq_category", "assembly_level", "version_status", "excluded_from_refseq", "selection"])
    tax_rows = []
    taxonomy = Taxonomy(Path(cfg["metadata_dir"]) / "taxonomy")
    for taxid, values in sorted(by_taxid.items(), key=lambda kv: (-kv[1]["bases"], kv[0])):
        info = taxonomy.describe(taxid)
        tax_rows.append({"source_taxid": taxid, "name": info["name"], "group": info["group"], **values})
    _write_tsv(out / "taxa.tsv", tax_rows,
               ["source_taxid", "name", "group", "references", "records", "bases", "compressed_bytes"])
    print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))
    print(f"AUDIT_DIR\t{out}")
    if args.require_human and not summary["human"]["present"]:
        raise SystemExit("ERROR: no frozen human reference lineage was found")


if __name__ == "__main__":
    main()
