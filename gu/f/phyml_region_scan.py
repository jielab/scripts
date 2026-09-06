#!/usr/bin/env python3
"""Anchor-independent regional archaic-haplotype inventory for GU PhyML.

This is deliberately the first interpretation layer after ``core.py`` has
constructed modern haplotypes for the complete BED core.  It does not use BED
column 4, ancestral-state filters, YRI/AFR frequencies, or a risk allele.
Consequently, a missing or incorrect anchor can never erase a regional
archaic-like match.

The module writes three non-destructive products:

* ``region_matches.unfiltered.tsv.gz``: every modern H type x archaic reference;
* ``region_match_segments.unfiltered.tsv.gz``: every local high-similarity tract;
* ``region_candidate_families.unfiltered.tsv.gz``: overlap-clustered tract families.

These are similarity candidates, not introgression calls.  Derived-state,
outgroup-frequency, recurrence, anchor linkage and candidate-focused tree
filters are applied later by ``evidence.py``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

from comm import enable_wide_csv_fields, read_tsv_rows as read_tsv, write_tsv_rows as write_tsv


enable_wide_csv_fields()

BASES = {"A", "C", "G", "T"}


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default


def as_float(value: object, default: float | None = None) -> float | None:
    try:
        x = float(str(value).strip())
        return x if math.isfinite(x) else default
    except (TypeError, ValueError):
        return default


def fmt(value: float | None) -> str:
    return "NA" if value is None or not math.isfinite(value) else f"{value:.12g}"


def truth(value: object) -> bool:
    return str(value).strip().upper() in {"1", "TRUE", "T", "YES", "Y", "PASS"}


def safe_slug(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._-")
    return text or "unknown"


def lineage_name(value: str) -> str:
    text = str(value).strip()
    if re.search(r"neander|altai|vindija|chagyr", text, re.I):
        return "Neanderthal"
    if re.search(r"denis", text, re.I):
        return "Denisovan"
    return text or "Archaic"


def write_text_if_changed(path: Path, text: str) -> bool:
    old = path.read_text(errors="replace") if path.is_file() else None
    if old == text:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return True


def remove_phyml_derivatives(phy: Path) -> None:
    for suffix in (
        "_phyml_tree.txt", "_phyml_stats.txt", "_phyml_boot_trees.txt",
        "_phyml_boot_stats.txt", "_phyml_tree.png", ".phyml.log",
        ".phyml.run.status.tsv",
    ):
        Path(f"{phy}{suffix}").unlink(missing_ok=True)


def parse_copies(value: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for token in str(value).split(";"):
        token = token.strip()
        if not token:
            continue
        sample, hap = token.rsplit(":", 1) if ":" in token else (token, "")
        if sample:
            out.append((sample, hap))
    return out


def locus_dirs(out: Path) -> list[Path]:
    root = out / "loci"
    if (root / "sites.tsv").is_file():
        return [root]
    if not root.is_dir():
        return []
    return sorted(path for path in root.iterdir() if path.is_dir() and (path / "sites.tsv").is_file())


def final_rows_by_locus(out: Path, filename: str) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(out / "final" / filename):
        grouped[row.get("locus_id", "")].append(row)
    return grouped


def infer_locus_id(locus_dir: Path, meta_map: dict[str, list[dict[str, str]]]) -> str:
    if locus_dir.name != "loci":
        return locus_dir.name
    ids = [key for key in meta_map if key]
    if len(ids) != 1:
        fail(f"flat locus directory requires exactly one locus summary; found {len(ids)}")
    return ids[0]


def validate_sequences(locus_id: str, n_sites: int, named: dict[str, str]) -> None:
    bad = {key: len(value) for key, value in named.items() if len(value) != n_sites}
    if bad:
        detail = ", ".join(f"{key}={value}" for key, value in sorted(bad.items()))
        fail(f"{locus_id}: sequence/site length mismatch: sites={n_sites}; {detail}")


def full_similarity(query: str, reference: str) -> tuple[int, int, float | None]:
    compared = matched = 0
    for q, r in zip(query, reference):
        if q not in BASES or r not in BASES:
            continue
        compared += 1
        matched += int(q == r)
    return compared, matched, matched / compared if compared else None


def split_callable_blocks(
    positions: Sequence[int], query: str, reference: str, max_gap_bp: int,
) -> list[list[tuple[int, int, int]]]:
    """Return blocks of (source_index, position, is_match) for callable sites."""
    blocks: list[list[tuple[int, int, int]]] = []
    current: list[tuple[int, int, int]] = []
    for i, (position, q, r) in enumerate(zip(positions, query, reference)):
        if q not in BASES or r not in BASES:
            continue
        if current and position - current[-1][1] > max_gap_bp:
            blocks.append(current)
            current = []
        current.append((i, position, int(q == r)))
    if current:
        blocks.append(current)
    return blocks


def scan_local_segments(
    positions: Sequence[int],
    query: str,
    reference: str,
    min_sites: int,
    seed_prop: float,
    max_gap_bp: int,
) -> list[dict[str, object]]:
    """Find local high-identity segments without anchor or population filters.

    Fixed-size qualifying windows seed a segment. Overlapping/adjacent seeds are
    merged, then greedily expanded while the complete segment retains the same
    minimum identity. This is intentionally permissive and descriptive; every
    later biological filter is stored separately.
    """
    out: list[dict[str, object]] = []
    for block in split_callable_blocks(positions, query, reference, max_gap_bp):
        n = len(block)
        if n < min_sites:
            continue
        prefix = [0]
        for _, _, match in block:
            prefix.append(prefix[-1] + match)
        seeds: list[tuple[int, int]] = []
        for left in range(0, n - min_sites + 1):
            right = left + min_sites - 1
            matches = prefix[right + 1] - prefix[left]
            if matches / min_sites >= seed_prop:
                seeds.append((left, right))
        if not seeds:
            total_matches = prefix[-1]
            if total_matches / n >= seed_prop:
                seeds = [(0, n - 1)]
            else:
                continue
        merged: list[list[int]] = []
        for left, right in seeds:
            if not merged or left > merged[-1][1] + 1:
                merged.append([left, right])
            else:
                merged[-1][1] = max(merged[-1][1], right)
        for left, right in merged:
            # Expand into immediately adjacent callable sites only while the
            # whole candidate retains the seed identity threshold.
            while left > 0:
                m = prefix[right + 1] - prefix[left]
                candidate_m = m + block[left - 1][2]
                candidate_n = right - left + 2
                if candidate_m / candidate_n < seed_prop:
                    break
                left -= 1
            while right + 1 < n:
                m = prefix[right + 1] - prefix[left]
                candidate_m = m + block[right + 1][2]
                candidate_n = right - left + 2
                if candidate_m / candidate_n < seed_prop:
                    break
                right += 1
            matches = prefix[right + 1] - prefix[left]
            callable_n = right - left + 1
            if callable_n < min_sites or matches / callable_n < seed_prop:
                continue
            source_indexes = [block[i][0] for i in range(left, right + 1)]
            start_pos = block[left][1]
            end_pos = block[right][1]
            out.append({
                "start": start_pos - 1,
                "end": end_pos,
                "span_bp": end_pos - start_pos + 1,
                "n_callable": callable_n,
                "n_match": matches,
                "n_mismatch": callable_n - matches,
                "prop_match": matches / callable_n,
                "source_indexes": source_indexes,
            })
    # Collapse exact duplicates that can arise when a complete block is also a
    # merged fixed-window segment.
    unique: dict[tuple[int, int, int, int], dict[str, object]] = {}
    for row in out:
        key = (as_int(row["start"]), as_int(row["end"]), as_int(row["n_callable"]), as_int(row["n_match"]))
        unique[key] = row
    return sorted(unique.values(), key=lambda row: (as_int(row["start"]), as_int(row["end"])))


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap <= 0:
        return 0.0
    return overlap / max(1, min(a_end - a_start, b_end - b_start))


def cluster_families(
    locus_id: str,
    segment_rows: list[dict[str, object]],
    min_overlap: float,
) -> tuple[list[dict[str, object]], dict[str, str]]:
    """Cluster overlapping local matches within each archaic lineage."""
    families: list[dict[str, object]] = []
    segment_to_family: dict[str, str] = {}
    serial = 0
    by_lineage: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in segment_rows:
        by_lineage[str(row["lineage"])].append(row)
    for lineage in sorted(by_lineage):
        rows = sorted(by_lineage[lineage], key=lambda r: (as_int(r["start"]), as_int(r["end"]), str(r["segment_id"])))
        uf = UnionFind(len(rows))
        for i, left in enumerate(rows):
            left_end = as_int(left["end"])
            for j in range(i + 1, len(rows)):
                right = rows[j]
                if as_int(right["start"]) >= left_end:
                    break
                if overlap_fraction(
                    as_int(left["start"]), left_end,
                    as_int(right["start"]), as_int(right["end"]),
                ) >= min_overlap:
                    uf.union(i, j)
        groups: dict[int, list[dict[str, object]]] = defaultdict(list)
        for i, row in enumerate(rows):
            groups[uf.find(i)].append(row)
        ordered = sorted(groups.values(), key=lambda g: (min(as_int(x["start"]) for x in g), max(as_int(x["end"]) for x in g)))
        for group in ordered:
            serial += 1
            family_id = f"F{serial:04d}_{safe_slug(lineage)}"
            hap_n: dict[str, int] = {}
            refs = set()
            for row in group:
                hap_n[str(row["hap_id"])] = max(hap_n.get(str(row["hap_id"]), 0), as_int(row["hap_n"]))
                refs.add(str(row["archaic"]))
                segment_to_family[str(row["segment_id"])] = family_id
            starts = [as_int(row["start"]) for row in group]
            ends = [as_int(row["end"]) for row in group]
            props = [as_float(row["prop_match"], 0.0) or 0.0 for row in group]
            families.append({
                "locus_id": locus_id,
                "family_id": family_id,
                "lineage": lineage,
                "start": min(starts),
                "end": max(ends),
                "span_bp": max(ends) - min(starts),
                "n_segments": len(group),
                "n_haplotype_types": len(hap_n),
                "n_haplotype_copies": sum(hap_n.values()),
                "n_archaic_references": len(refs),
                "archaic_references": ",".join(sorted(refs)),
                "haplotype_ids": ",".join(sorted(hap_n)),
                "segment_ids": ",".join(sorted(str(row["segment_id"]) for row in group)),
                "max_prop_match": fmt(max(props) if props else None),
                "median_prop_match": fmt(median(props) if props else None),
                "filter_state": "UNFILTERED",
                "interpretation": "anchor_independent_archaic_similarity_family_not_introgression_call",
            })
    return families, segment_to_family


def make_region_common_tree(
    locus_dir: Path,
    locus_id: str,
    sites: list[dict[str, str]],
    haps: list[dict[str, str]],
    archaic: list[dict[str, str]],
    ancestral: str,
    families: list[dict[str, object]],
    min_copies: int,
    min_sites: int,
    lineage_filter: str = "",
) -> dict[str, object]:
    """Build an anchor-independent common-haplotype tree.

    ``lineage_filter`` creates an old-analysis-compatible lineage-specific tree
    in addition to the broad all-archaic tree. Neither tree uses the anchor.
    """
    family_by_hap: dict[str, list[str]] = defaultdict(list)
    for family in families:
        for hap_id in str(family.get("haplotype_ids", "")).split(","):
            if hap_id:
                family_by_hap[hap_id].append(str(family["family_id"]))
    total_haplotype_copies = sum(as_int(row.get("n"), len(parse_copies(row.get("copies", "")))) for row in haps)
    common = [row for row in haps if as_int(row.get("n")) >= min_copies]
    records: list[dict[str, object]] = []
    for row in sorted(common, key=lambda r: (-as_int(r.get("n")), r.get("hap_id", ""))):
        records.append({
            "label": row.get("hap_id", ""),
            "seq": row.get("seq", "").upper(),
            "role": "modern_common",
            "lineage": "Modern",
            "n": row.get("n", ""),
            "total_haplotype_copies": total_haplotype_copies,
            "hap_frequency": fmt(as_int(row.get("n")) / total_haplotype_copies if total_haplotype_copies else None),
            "copies": row.get("copies", ""),
            "family_ids": ",".join(sorted(family_by_hap.get(row.get("hap_id", ""), []))),
        })
    for row in archaic:
        row_lineage = lineage_name(row.get("lineage") or row.get("archaic", ""))
        if lineage_filter and row_lineage != lineage_filter:
            continue
        records.append({
            "label": row.get("archaic", "Archaic"),
            "seq": row.get("seq", "").upper(),
            "role": "archaic",
            "lineage": row_lineage,
            "n": 1,
            "total_haplotype_copies": total_haplotype_copies,
            "hap_frequency": "",
            "copies": "",
            "family_ids": "",
        })
    if ancestral and sum(base in BASES for base in ancestral) >= min_sites:
        records.append({
            "label": "Ancestral",
            "seq": ancestral,
            "role": "outgroup",
            "lineage": "Ancestral",
            "n": 1,
            "total_haplotype_copies": total_haplotype_copies,
            "hap_frequency": "",
            "copies": "",
            "family_ids": "",
        })

    suffix = f".{safe_slug(lineage_filter)}" if lineage_filter else ""
    tree_scope = "region_common_unfiltered_lineage" if lineage_filter else "region_common_unfiltered_all_archaic"
    phy = locus_dir / f"haplotypes.region_common{suffix}.phy"
    meta = Path(f"{phy}.meta.tsv")
    manifest = locus_dir / f"haplotypes.region_common{suffix}.manifest.tsv"
    site_file = locus_dir / f"haplotypes.region_common{suffix}.sites.tsv"
    if len(common) < 2:
        for path in (phy, meta, manifest, site_file):
            path.unlink(missing_ok=True)
        remove_phyml_derivatives(phy)
        return {
            "locus_id": locus_id,
            "tree_scope": tree_scope,
            "lineage_filter": lineage_filter or "ALL",
            "input_status": "too_few_common_haplotype_types",
            "common_min_copies": min_copies,
            "n_common_haplotype_types": len(common),
            "n_tree_sites": 0,
            "n_tree_tips": 0,
            "phy_file": "",
        }

    keep_indexes: list[int] = []
    archaic_sequences = [str(row["seq"]) for row in records if row["role"] == "archaic"]
    for i in range(len(sites)):
        if not any(seq[i] in BASES for seq in archaic_sequences):
            continue
        called = {str(row["seq"])[i] for row in records if str(row["seq"])[i] in BASES}
        if len(called) >= 2:
            keep_indexes.append(i)
    if len(keep_indexes) < min_sites:
        for path in (phy, meta, manifest, site_file):
            path.unlink(missing_ok=True)
        remove_phyml_derivatives(phy)
        return {
            "locus_id": locus_id,
            "tree_scope": tree_scope,
            "lineage_filter": lineage_filter or "ALL",
            "input_status": "too_few_variable_archaic_callable_sites",
            "common_min_copies": min_copies,
            "n_common_haplotype_types": len(common),
            "n_tree_sites": len(keep_indexes),
            "n_tree_tips": 0,
            "phy_file": "",
        }

    phy_lines = [f"{len(records)} {len(keep_indexes)}"]
    manifest_rows: list[dict[str, object]] = []
    for index, row in enumerate(records, 1):
        short = f"R{index:09d}"
        sequence = "".join(str(row["seq"])[i] for i in keep_indexes)
        phy_lines.append(f"{short:<10} {sequence}")
        manifest_rows.append({
            "phy_label": short,
            "label": row["label"],
            "role": row["role"],
            "lineage": row["lineage"],
            "n": row["n"],
            "total_haplotype_copies": row["total_haplotype_copies"],
            "hap_frequency": row["hap_frequency"],
            "copies": row["copies"],
            "family_ids": row["family_ids"],
            "tree_scope": tree_scope,
            "lineage_filter": lineage_filter or "ALL",
            "common_min_copies": min_copies,
            "n_tree_sites": len(keep_indexes),
        })
    changed = write_text_if_changed(phy, "\n".join(phy_lines) + "\n")
    if changed:
        remove_phyml_derivatives(phy)
    write_tsv(meta, ["phy_label", "label"], manifest_rows)
    write_tsv(
        manifest,
        ["phy_label", "label", "role", "lineage", "n", "total_haplotype_copies", "hap_frequency", "copies", "family_ids", "tree_scope", "lineage_filter", "common_min_copies", "n_tree_sites"],
        manifest_rows,
    )
    write_tsv(
        site_file,
        ["tree_index", "source_index", "chr", "pos", "id", "ref", "alt"],
        (
            {
                "tree_index": tree_i + 1,
                "source_index": source_i + 1,
                "chr": sites[source_i].get("chr", ""),
                "pos": sites[source_i].get("pos", ""),
                "id": sites[source_i].get("id", ""),
                "ref": sites[source_i].get("ref", ""),
                "alt": sites[source_i].get("alt", ""),
            }
            for tree_i, source_i in enumerate(keep_indexes)
        ),
    )
    return {
        "locus_id": locus_id,
        "tree_scope": tree_scope,
        "lineage_filter": lineage_filter or "ALL",
        "input_status": "ready",
        "common_min_copies": min_copies,
        "n_common_haplotype_types": len(common),
        "n_tree_sites": len(keep_indexes),
        "n_tree_tips": len(records),
        "phy_file": str(phy.resolve()),
    }


def prepare(args: argparse.Namespace) -> None:
    out: Path = args.out
    meta_map = final_rows_by_locus(out, "loci.tsv")
    dirs = locus_dirs(out)
    if not dirs and not meta_map:
        fail(f"no PhyML locus summary or artifacts under {out}")

    all_haplotype_rows: list[dict[str, object]] = []
    all_archaic_rows: list[dict[str, object]] = []
    all_match_rows: list[dict[str, object]] = []
    all_segment_rows: list[dict[str, object]] = []
    all_similarity_candidate_rows: list[dict[str, object]] = []
    all_family_rows: list[dict[str, object]] = []
    all_locus_rows: list[dict[str, object]] = []
    all_tree_inputs: list[dict[str, object]] = []

    represented_loci = {infer_locus_id(path, meta_map) for path in dirs}
    for missing_locus_id, meta_rows in meta_map.items():
        if not missing_locus_id or missing_locus_id in represented_loci:
            continue
        meta = meta_rows[0] if meta_rows else {}
        raw_status = str(meta.get("status", "") or "missing_locus_artifacts")
        all_locus_rows.append({
            "locus_id": missing_locus_id,
            "chr": meta.get("chr", ""),
            "core_start": meta.get("core_start", ""),
            "core_end": meta.get("core_end", ""),
            "x_class": meta.get("x_class", ""),
            "n_sites": as_int(meta.get("n_sites")),
            "n_haplotype_types": as_int(meta.get("n_haplotypes")),
            "n_haplotype_copies": 0,
            "haplotype_denominator": "unavailable",
            "n_archaic_references": 0,
            "n_pairwise_matches": 0,
            "n_full_similarity_pass": 0,
            "n_similarity_candidates": 0,
            "n_local_segments": 0,
            "n_unfiltered_families": 0,
            "stage1_status": f"not_evaluable:{raw_status}",
            "anchor_used": 0,
            "ancestral_filter_used": 0,
            "population_filter_used": 0,
            "interpretation": "raw_sequence_artifact_unavailable;not_a_negative_introgression_result",
        })

    for locus_dir in dirs:
        locus_id = infer_locus_id(locus_dir, meta_map)
        meta = meta_map.get(locus_id, [{}])[0]
        sites = read_tsv(locus_dir / "sites.tsv")
        haps = read_tsv(locus_dir / "haplotypes.tsv")
        archaic = read_tsv(locus_dir / "archaic.tsv")
        ancestral_rows = read_tsv(locus_dir / "ancestral.tsv")
        if not sites or not haps or not archaic:
            all_locus_rows.append({
                "locus_id": locus_id,
                "chr": meta.get("chr", ""),
                "core_start": meta.get("core_start", ""),
                "core_end": meta.get("core_end", ""),
                "x_class": meta.get("x_class", ""),
                "n_sites": len(sites),
                "n_haplotype_types": len(haps),
                "n_haplotype_copies": sum(as_int(row.get("n")) for row in haps),
                "haplotype_denominator": "unavailable_or_partial",
                "n_archaic_references": len(archaic),
                "n_pairwise_matches": 0,
                "n_full_similarity_pass": 0,
                "n_similarity_candidates": 0,
                "n_local_segments": 0,
                "n_unfiltered_families": 0,
                "stage1_status": "missing_sites_haplotypes_or_archaic",
                "anchor_used": 0,
                "ancestral_filter_used": 0,
                "population_filter_used": 0,
                "interpretation": "not_evaluable;missing_raw_artifact_not_a_negative_introgression_result",
            })
            continue
        positions = [as_int(row.get("pos"), -1) for row in sites]
        if any(pos < 0 for pos in positions):
            fail(f"{locus_id}: malformed sites.tsv positions")
        ancestral = ancestral_rows[0].get("seq", "").upper() if ancestral_rows else "N" * len(sites)
        named = {"Ancestral": ancestral}
        named.update({row.get("archaic", "Archaic"): row.get("seq", "").upper() for row in archaic})
        named.update({row.get("hap_id", "H"): row.get("seq", "").upper() for row in haps})
        validate_sequences(locus_id, len(sites), named)
        total_haplotype_copies = sum(
            as_int(row.get("n"), len(parse_copies(row.get("copies", "")))) for row in haps
        )

        # Freeze an explicit Stage-1 snapshot before any diagnostic, population,
        # recurrence, or anchor-based annotation is applied.  final/haplotypes.tsv
        # remains backward-compatible and is later enriched with evidence fields;
        # this table is the immutable regional discovery inventory.
        for hap in haps:
            all_haplotype_rows.append({
                "locus_id": locus_id,
                "chr": meta.get("chr", sites[0].get("chr", "")),
                "core_start": meta.get("core_start", ""),
                "core_end": meta.get("core_end", ""),
                "hap_id": hap.get("hap_id", ""),
                "hap_n": as_int(hap.get("n"), len(parse_copies(hap.get("copies", "")))),
                "total_haplotype_copies": total_haplotype_copies,
                "hap_frequency": fmt(
                    as_int(hap.get("n"), len(parse_copies(hap.get("copies", "")))) / total_haplotype_copies
                    if total_haplotype_copies else None
                ),
                "copies": hap.get("copies", ""),
                "seq": hap.get("seq", "").upper(),
                "anchor_used": 0,
                "ancestral_filter_used": 0,
                "population_filter_used": 0,
                "recurrence_filter_used": 0,
                "filter_state": "UNFILTERED",
                "interpretation": "complete_anchor_independent_modern_haplotype_inventory",
            })
        for ref in archaic:
            all_archaic_rows.append({
                "locus_id": locus_id,
                "chr": meta.get("chr", sites[0].get("chr", "")),
                "core_start": meta.get("core_start", ""),
                "core_end": meta.get("core_end", ""),
                "archaic": ref.get("archaic", "Archaic"),
                "lineage": lineage_name(ref.get("lineage") or ref.get("archaic", "")),
                "seq": ref.get("seq", "").upper(),
                "anchor_used": 0,
                "ancestral_filter_used": 0,
                "population_filter_used": 0,
                "filter_state": "UNFILTERED",
                "interpretation": "regional_archaic_reference_projection",
            })

        locus_segments: list[dict[str, object]] = []
        segment_serial = 0
        for hap in haps:
            hap_id = hap.get("hap_id", "")
            query = hap.get("seq", "").upper()
            hap_n = as_int(hap.get("n"), len(parse_copies(hap.get("copies", ""))))
            pair_rows: list[dict[str, object]] = []
            for ref in archaic:
                ref_name = ref.get("archaic", "Archaic")
                ref_seq = ref.get("seq", "").upper()
                lineage = lineage_name(ref.get("lineage") or ref_name)
                compared, matched, prop = full_similarity(query, ref_seq)
                pair = {
                    "locus_id": locus_id,
                    "chr": meta.get("chr", sites[0].get("chr", "")),
                    "core_start": meta.get("core_start", ""),
                    "core_end": meta.get("core_end", ""),
                    "hap_id": hap_id,
                    "hap_n": hap_n,
                    "total_haplotype_copies": total_haplotype_copies,
                    "hap_frequency": fmt(hap_n / total_haplotype_copies if total_haplotype_copies else None),
                    "copies": hap.get("copies", ""),
                    "archaic": ref_name,
                    "lineage": lineage,
                    "n_compared": compared,
                    "n_match": matched,
                    "n_mismatch": compared - matched,
                    "prop_match": fmt(prop),
                    "full_similarity_pass": int(
                        compared >= args.full_min_callable_sites
                        and prop is not None
                        and prop >= args.full_min_prop
                    ),
                    "rank_within_haplotype": "",
                    "anchor_used": 0,
                    "ancestral_filter_used": 0,
                    "population_filter_used": 0,
                    "interpretation": "unfiltered_full_region_acgt_identity_not_introgression_call",
                }
                pair_rows.append(pair)
                local = scan_local_segments(
                    positions,
                    query,
                    ref_seq,
                    args.scan_min_callable_sites,
                    args.scan_seed_prop,
                    args.scan_max_gap_bp,
                )
                for local_index, segment in enumerate(local, 1):
                    segment_serial += 1
                    segment_id = f"S{segment_serial:06d}"
                    source_indexes = list(segment.pop("source_indexes"))
                    row = {
                        "locus_id": locus_id,
                        "segment_id": segment_id,
                        "family_id": "",
                        "chr": meta.get("chr", sites[0].get("chr", "")),
                        "hap_id": hap_id,
                        "hap_n": hap_n,
                        "total_haplotype_copies": total_haplotype_copies,
                        "hap_frequency": fmt(hap_n / total_haplotype_copies if total_haplotype_copies else None),
                        "copies": hap.get("copies", ""),
                        "archaic": ref_name,
                        "lineage": lineage,
                        "local_rank_within_pair": local_index,
                        **segment,
                        "first_source_site_index": source_indexes[0] + 1 if source_indexes else "",
                        "last_source_site_index": source_indexes[-1] + 1 if source_indexes else "",
                        "source_site_indexes": ",".join(str(i + 1) for i in source_indexes),
                        "scan_min_callable_sites": args.scan_min_callable_sites,
                        "scan_seed_prop": args.scan_seed_prop,
                        "scan_max_gap_bp": args.scan_max_gap_bp,
                        "anchor_used": 0,
                        "ancestral_filter_used": 0,
                        "population_filter_used": 0,
                        "filter_state": "UNFILTERED",
                        "interpretation": "anchor_independent_local_archaic_similarity_not_introgression_call",
                    }
                    locus_segments.append(row)
            pair_rows.sort(
                key=lambda row: (
                    -(as_float(row["prop_match"], -1.0) or -1.0),
                    -as_int(row["n_match"]),
                    str(row["archaic"]),
                )
            )
            for rank, pair in enumerate(pair_rows, 1):
                pair["rank_within_haplotype"] = rank
                all_match_rows.append(pair)

        families, seg_family = cluster_families(locus_id, locus_segments, args.family_min_overlap)
        for row in locus_segments:
            row["family_id"] = seg_family.get(str(row["segment_id"]), "")
        all_segment_rows.extend(locus_segments)
        all_family_rows.extend(families)

        # A compact Stage-1 table of every similarity candidate.  This table is
        # still anchor/ancestral/population agnostic and deliberately retains
        # singleton modern haplotypes.  The exhaustive all-pairs matrix remains
        # in region_matches.unfiltered.tsv.gz.
        segments_by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
        for segment in locus_segments:
            segments_by_pair[(str(segment.get("hap_id", "")), str(segment.get("archaic", "")))].append(segment)
        for pair in [row for row in all_match_rows if row.get("locus_id") == locus_id]:
            local = segments_by_pair.get((str(pair.get("hap_id", "")), str(pair.get("archaic", ""))), [])
            full_pass = as_int(pair.get("full_similarity_pass")) == 1
            if not full_pass and not local:
                continue
            family_ids = sorted({str(row.get("family_id", "")) for row in local if row.get("family_id")})
            best_local = max(local, key=lambda row: (
                as_float(row.get("prop_match"), -1.0) or -1.0,
                as_int(row.get("n_callable")),
                as_int(row.get("span_bp")),
            ), default=None)
            reasons = []
            if full_pass:
                reasons.append("full_region_similarity")
            if local:
                reasons.append("local_similarity_segment")
            all_similarity_candidate_rows.append({
                **pair,
                "n_local_segments": len(local),
                "family_ids": ",".join(family_ids),
                "best_local_start": best_local.get("start", "") if best_local else "",
                "best_local_end": best_local.get("end", "") if best_local else "",
                "best_local_span_bp": best_local.get("span_bp", "") if best_local else "",
                "best_local_n_callable": best_local.get("n_callable", "") if best_local else "",
                "best_local_n_match": best_local.get("n_match", "") if best_local else "",
                "best_local_prop_match": best_local.get("prop_match", "") if best_local else "",
                "similarity_candidate_reason": ";".join(reasons),
                "filter_state": "UNFILTERED",
                "interpretation": "anchor_independent_similarity_candidate_not_introgression_call",
            })
        common_tree_lineages = [""] + sorted({
            lineage_name(row.get("lineage") or row.get("archaic", ""))
            for row in archaic
            if lineage_name(row.get("lineage") or row.get("archaic", ""))
        })
        for common_lineage in common_tree_lineages:
            tree_input = make_region_common_tree(
                locus_dir,
                locus_id,
                sites,
                haps,
                archaic,
                ancestral,
                families,
                args.common_tree_min_copies,
                args.common_tree_min_sites,
                common_lineage,
            )
            all_tree_inputs.append(tree_input)
        desired_common_phy = {
            Path(str(row.get("phy_file")))
            for row in all_tree_inputs
            if row.get("locus_id") == locus_id and row.get("phy_file")
        }
        for stale_phy in locus_dir.glob("haplotypes.region_common*.phy"):
            if stale_phy.resolve() in {path.resolve() for path in desired_common_phy}:
                continue
            remove_phyml_derivatives(stale_phy)
            stale_phy.unlink(missing_ok=True)
            Path(f"{stale_phy}.meta.tsv").unlink(missing_ok=True)
            stale_phy.with_name(stale_phy.name.replace(".phy", ".manifest.tsv")).unlink(missing_ok=True)
            stale_phy.with_name(stale_phy.name.replace(".phy", ".sites.tsv")).unlink(missing_ok=True)

        all_locus_rows.append({
            "locus_id": locus_id,
            "chr": meta.get("chr", sites[0].get("chr", "")),
            "core_start": meta.get("core_start", ""),
            "core_end": meta.get("core_end", ""),
            "x_class": meta.get("x_class", ""),
            "n_sites": len(sites),
            "n_haplotype_types": len(haps),
            "n_haplotype_copies": total_haplotype_copies,
            "haplotype_denominator": (
                "male_nonpar_X_chromosome_copies" if meta.get("x_class") == "male_nonpar"
                else "diploid_or_PAR_chromosome_copies"
            ),
            "n_archaic_references": len(archaic),
            "n_pairwise_matches": len(haps) * len(archaic),
            "n_full_similarity_pass": sum(
                as_int(row["full_similarity_pass"])
                for row in all_match_rows if row["locus_id"] == locus_id
            ),
            "n_similarity_candidates": sum(1 for row in all_similarity_candidate_rows if row.get("locus_id") == locus_id),
            "n_local_segments": len(locus_segments),
            "n_unfiltered_families": len(families),
            "stage1_status": "complete",
            "anchor_used": 0,
            "ancestral_filter_used": 0,
            "population_filter_used": 0,
            "interpretation": "complete_anchor_independent_inventory;not_an_introgression_call",
        })

    final = out / "final"
    write_tsv(
        final / "region_haplotypes.unfiltered.tsv.gz",
        [
            "locus_id", "chr", "core_start", "core_end", "hap_id", "hap_n",
            "total_haplotype_copies", "hap_frequency", "copies", "seq",
            "anchor_used", "ancestral_filter_used", "population_filter_used",
            "recurrence_filter_used", "filter_state", "interpretation",
        ],
        all_haplotype_rows,
    )
    write_tsv(
        final / "region_archaic_references.unfiltered.tsv.gz",
        [
            "locus_id", "chr", "core_start", "core_end", "archaic", "lineage", "seq",
            "anchor_used", "ancestral_filter_used", "population_filter_used",
            "filter_state", "interpretation",
        ],
        all_archaic_rows,
    )
    write_tsv(
        final / "region_matches.unfiltered.tsv.gz",
        [
            "locus_id", "chr", "core_start", "core_end", "hap_id", "hap_n",
            "total_haplotype_copies", "hap_frequency", "copies",
            "archaic", "lineage", "n_compared", "n_match", "n_mismatch", "prop_match",
            "full_similarity_pass", "rank_within_haplotype", "anchor_used",
            "ancestral_filter_used", "population_filter_used", "interpretation",
        ],
        all_match_rows,
    )
    write_tsv(
        final / "region_haplotype_candidates.unfiltered.tsv.gz",
        [
            "locus_id", "chr", "core_start", "core_end", "hap_id", "hap_n",
            "total_haplotype_copies", "hap_frequency", "copies",
            "archaic", "lineage", "n_compared", "n_match", "n_mismatch", "prop_match",
            "full_similarity_pass", "rank_within_haplotype", "n_local_segments", "family_ids",
            "best_local_start", "best_local_end", "best_local_span_bp", "best_local_n_callable",
            "best_local_n_match", "best_local_prop_match", "similarity_candidate_reason",
            "anchor_used", "ancestral_filter_used", "population_filter_used", "filter_state",
            "interpretation",
        ],
        all_similarity_candidate_rows,
    )
    write_tsv(
        final / "region_match_segments.unfiltered.tsv.gz",
        [
            "locus_id", "segment_id", "family_id", "chr", "hap_id", "hap_n",
            "total_haplotype_copies", "hap_frequency", "copies",
            "archaic", "lineage", "local_rank_within_pair", "start", "end", "span_bp",
            "n_callable", "n_match", "n_mismatch", "prop_match", "first_source_site_index",
            "last_source_site_index", "source_site_indexes", "scan_min_callable_sites",
            "scan_seed_prop", "scan_max_gap_bp", "anchor_used", "ancestral_filter_used",
            "population_filter_used", "filter_state", "interpretation",
        ],
        all_segment_rows,
    )
    write_tsv(
        final / "region_candidate_families.unfiltered.tsv.gz",
        [
            "locus_id", "family_id", "lineage", "start", "end", "span_bp", "n_segments",
            "n_haplotype_types", "n_haplotype_copies", "n_archaic_references",
            "archaic_references", "haplotype_ids", "segment_ids", "max_prop_match",
            "median_prop_match", "filter_state", "interpretation",
        ],
        all_family_rows,
    )
    write_tsv(
        final / "region_scan_loci.unfiltered.tsv.gz",
        [
            "locus_id", "chr", "core_start", "core_end", "x_class", "n_sites", "n_haplotype_types",
            "n_haplotype_copies", "haplotype_denominator", "n_archaic_references", "n_pairwise_matches",
            "n_full_similarity_pass", "n_similarity_candidates", "n_local_segments", "n_unfiltered_families",
            "stage1_status", "anchor_used", "ancestral_filter_used", "population_filter_used",
            "interpretation",
        ],
        all_locus_rows,
    )
    write_tsv(
        final / "region_common_tree_inputs.tsv",
        [
            "locus_id", "tree_scope", "lineage_filter", "input_status", "common_min_copies",
            "n_common_haplotype_types", "n_tree_sites", "n_tree_tips", "phy_file",
        ],
        all_tree_inputs,
    )
    parameters = {
        "schema_version": 4,
        "stage": "anchor_independent_regional_discovery",
        "anchor_is_used": False,
        "ancestral_filter_is_used": False,
        "population_filter_is_used": False,
        "recurrence_filter_is_used_for_inventory": False,
        "outputs": {
            "all_modern_haplotypes": "region_haplotypes.unfiltered.tsv.gz",
            "all_archaic_references": "region_archaic_references.unfiltered.tsv.gz",
            "all_haplotype_reference_pairs": "region_matches.unfiltered.tsv.gz",
            "permissive_similarity_candidates": "region_haplotype_candidates.unfiltered.tsv.gz",
            "local_similarity_segments": "region_match_segments.unfiltered.tsv.gz",
            "overlap_clustered_families": "region_candidate_families.unfiltered.tsv.gz",
            "common_haplotype_tree_inputs": "region_common_tree_inputs.tsv",
        },
        "thresholds": {
            "scan_min_callable_sites": args.scan_min_callable_sites,
            "scan_seed_prop": args.scan_seed_prop,
            "scan_max_gap_bp": args.scan_max_gap_bp,
            "full_min_callable_sites": args.full_min_callable_sites,
            "full_min_prop": args.full_min_prop,
            "family_min_overlap": args.family_min_overlap,
            "common_tree_min_copies": args.common_tree_min_copies,
            "common_tree_min_sites": args.common_tree_min_sites,
            "region_tree_bootstrap_min": args.min_tree_bootstrap,
            "region_tree_min_modern_types": args.min_tree_modern_types,
        },
        "warning": "Stage-1 similarity is an exhaustive discovery inventory, not an introgression call.",
    }
    write_text_if_changed(
        final / "region_scan_parameters.json",
        json.dumps(parameters, indent=2, sort_keys=True) + "\n",
    )
    print(
        f"PHYML REGION SCAN: loci={len(all_locus_rows)} haplotypes={len(all_haplotype_rows)} "
        f"all_pairs={len(all_match_rows)} "
        f"similarity_candidates={len(all_similarity_candidate_rows)} "
        f"local_segments={len(all_segment_rows)} unfiltered_families={len(all_family_rows)} "
        f"output={final}",
        flush=True,
    )


@dataclass
class NewickNode:
    label: str = ""
    children: list["NewickNode"] = field(default_factory=list)
    support: float | None = None
    node_id: str = ""


def parse_newick(text: str) -> NewickNode:
    source = "".join(text.split())
    index = 0
    serial = 0

    def token() -> str:
        nonlocal index
        start = index
        while index < len(source) and source[index] not in ":,();":
            index += 1
        return source[start:index]

    def branch_length() -> None:
        nonlocal index
        if index < len(source) and source[index] == ":":
            index += 1
            while index < len(source) and source[index] not in ",();":
                index += 1

    def subtree() -> NewickNode:
        nonlocal index, serial
        if index >= len(source):
            raise ValueError("unexpected end of Newick")
        if source[index] != "(":
            label = token()
            if not label:
                raise ValueError("empty Newick tip")
            node = NewickNode(label=label)
            branch_length()
            return node
        index += 1
        children = [subtree()]
        while index < len(source) and source[index] == ",":
            index += 1
            children.append(subtree())
        if index >= len(source) or source[index] != ")":
            raise ValueError("unterminated Newick clade")
        index += 1
        label = token()
        support = None
        try:
            support = float(label) if label else None
        except ValueError:
            pass
        serial += 1
        node = NewickNode(label=label, children=children, support=support, node_id=f"N{serial}")
        branch_length()
        return node

    root = subtree()
    if index < len(source) and source[index] == ";":
        index += 1
    if index != len(source):
        raise ValueError(f"unexpected Newick content at offset {index}")
    return root


def relabel_newick(text: str, labels: dict[str, str]) -> str:
    for short in sorted(labels, key=len, reverse=True):
        text = re.sub(rf"(?<=[(,]){re.escape(short)}(?=[:),;])", labels[short], text)
    return "".join(text.split())


def bootstrap_values(newick: str) -> list[float]:
    return [
        float(value)
        for value in re.findall(r"\)(-?[0-9]+(?:\.[0-9]+)?)(?=[:),;])", newick)
        if float(value) >= 0
    ]


def enumerate_archaic_modern_clades(
    newick: str,
    manifest: list[dict[str, str]],
    min_bootstrap: float,
    min_modern_tips: int,
) -> list[dict[str, object]]:
    """Enumerate every single-lineage archaic + modern edge in the common tree."""
    if not newick:
        return []
    by_label = {row.get("label", ""): row for row in manifest}
    root = parse_newick(newick)
    descendants: dict[str, set[str]] = {}
    internal: list[NewickNode] = []

    def visit(node: NewickNode) -> set[str]:
        if not node.children:
            return {node.label}
        tips = set().union(*(visit(child) for child in node.children))
        descendants[node.node_id] = tips
        internal.append(node)
        return tips

    all_tips = visit(root)
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, ...]] = set()
    for node in internal:
        if node is root:
            continue
        down = descendants[node.node_id]
        for side_name, side in (("descendant", down), ("complement", all_tips - down)):
            if "Ancestral" in side:
                continue
            archaic = {tip for tip in side if by_label.get(tip, {}).get("role") == "archaic"}
            modern = {tip for tip in side if by_label.get(tip, {}).get("role") == "modern_common"}
            if not archaic or not modern:
                continue
            lineages = {by_label[tip].get("lineage", "") for tip in archaic}
            if len(lineages) != 1:
                continue
            lineage = next(iter(lineages))
            key = tuple(sorted(side))
            if key in seen:
                continue
            seen.add(key)
            support = node.support if node.support is not None and node.support >= 0 else None
            modern_copies = sum(as_int(by_label[tip].get("n")) for tip in modern)
            family_ids = sorted({
                family
                for tip in modern
                for family in by_label[tip].get("family_ids", "").split(",")
                if family
            })
            rows.append({
                "node": node.node_id,
                "side": side_name,
                "lineage": lineage,
                "bootstrap": fmt(support),
                "n_tips": len(side),
                "n_modern_haplotype_types": len(modern),
                "n_modern_haplotype_copies": modern_copies,
                "n_archaic_tips": len(archaic),
                "modern_haplotypes": ",".join(sorted(modern)),
                "archaic_references": ",".join(sorted(archaic)),
                "family_ids": ",".join(family_ids),
                "clade_pass": int(
                    support is not None
                    and support >= min_bootstrap
                    and len(modern) >= min_modern_tips
                ),
                "rule": f"unfiltered_region_common_tree;bootstrap>={min_bootstrap:g};modern_types>={min_modern_tips}",
                "interpretation": "tree_supported_archaic_affinity_not_outgroup_or_anchor_filtered",
            })
    rows.sort(
        key=lambda row: (
            -as_int(row["clade_pass"]),
            -(as_float(row["bootstrap"], -1.0) or -1.0),
            -as_int(row["n_modern_haplotype_copies"]),
            str(row["lineage"]),
            str(row["node"]),
        )
    )
    return rows


def summarize(args: argparse.Namespace) -> None:
    """Summarize every anchor-independent common-haplotype tree.

    The broad all-archaic tree and each lineage-specific tree are retained as
    separate rows.  These Stage-1 topologies are descriptive regional evidence;
    they do not use the named anchor, ancestral-derived-state filter, or
    YRI/AFR frequency filter and therefore are never promoted directly to an
    introgression call.
    """
    out: Path = args.out
    meta_map = final_rows_by_locus(out, "loci.tsv")
    tree_rows: list[dict[str, object]] = []
    clade_rows: list[dict[str, object]] = []
    family_rows = read_tsv(out / "final" / "region_candidate_families.unfiltered.tsv.gz")
    # family_id -> (pass, bootstrap, lineage_filter, node)
    family_best: dict[str, tuple[int, float, str, str]] = {}

    for locus_dir in locus_dirs(out):
        locus_id = infer_locus_id(locus_dir, meta_map)
        manifest_paths = sorted(locus_dir.glob("haplotypes.region_common*.manifest.tsv"))
        if not manifest_paths:
            # Preserve an explicit row when preparation could not build any
            # common-haplotype input for this locus.
            tree_rows.append({
                "locus_id": locus_id,
                "tree_scope": "region_common_unfiltered",
                "lineage_filter": "ALL",
                "tree_status": "input_not_ready",
                "phy_file": "",
                "tree_file": "",
                "stats_file": "",
                "plot_file": "",
                "n_tree_tips": 0,
                "n_bootstrap_nodes": 0,
                "bootstrap_min": "",
                "bootstrap_median": "",
                "bootstrap_max": "",
                "n_archaic_modern_clades": 0,
                "n_passing_archaic_modern_clades": 0,
                "best_clade_lineage": "",
                "best_clade_bootstrap": "",
                "best_clade_modern_haplotype_types": "",
                "best_clade_modern_haplotype_copies": "",
                "interpretation": "anchor_independent_common_haplotype_tree_input_not_ready",
                "tree_newick": "",
            })
            continue

        for manifest_path in manifest_paths:
            stem = manifest_path.name.removesuffix(".manifest.tsv")
            phy = locus_dir / f"{stem}.phy"
            meta_path = Path(f"{phy}.meta.tsv")
            tree = Path(f"{phy}_phyml_tree.txt")
            stats = Path(f"{phy}_phyml_stats.txt")
            manifest = read_tsv(manifest_path)
            lineage_filter = manifest[0].get("lineage_filter", "ALL") if manifest else "ALL"
            tree_scope = manifest[0].get("tree_scope", "region_common_unfiltered") if manifest else "region_common_unfiltered"

            if phy.is_file() and tree.is_file() and tree.stat().st_size and meta_path.is_file():
                label_map = {
                    row["phy_label"]: row["label"]
                    for row in read_tsv(meta_path)
                    if row.get("phy_label") and row.get("label")
                }
                newick = relabel_newick(tree.read_text(errors="replace"), label_map)
                supports = bootstrap_values(newick)
                status = "complete"
                clades = enumerate_archaic_modern_clades(
                    newick,
                    manifest,
                    args.min_tree_bootstrap,
                    args.min_tree_modern_types,
                )
            else:
                newick = ""
                supports = []
                status = args.default_tree_status if phy.is_file() else "input_not_ready"
                clades = []

            for rank, clade in enumerate(clades, 1):
                full = {
                    "locus_id": locus_id,
                    "tree_scope": tree_scope,
                    "lineage_filter": lineage_filter,
                    "clade_rank": rank,
                    **clade,
                }
                clade_rows.append(full)
                for family in str(clade.get("family_ids", "")).split(","):
                    if not family:
                        continue
                    score = (
                        as_int(clade.get("clade_pass")),
                        as_float(clade.get("bootstrap"), -1.0) or -1.0,
                        lineage_filter,
                        str(clade.get("node", "")),
                    )
                    if score > family_best.get(family, (0, -1.0, "", "")):
                        family_best[family] = score

            tree_rows.append({
                "locus_id": locus_id,
                "tree_scope": tree_scope,
                "lineage_filter": lineage_filter,
                "tree_status": status,
                "phy_file": str(phy.resolve()) if phy.is_file() else "",
                "tree_file": str(tree.resolve()) if tree.is_file() else "",
                "stats_file": str(stats.resolve()) if stats.is_file() else "",
                "plot_file": str(tree.with_suffix(".png").resolve()) if tree.with_suffix(".png").is_file() else "",
                "n_tree_tips": len(manifest),
                "n_bootstrap_nodes": len(supports),
                "bootstrap_min": min(supports) if supports else "",
                "bootstrap_median": median(supports) if supports else "",
                "bootstrap_max": max(supports) if supports else "",
                "n_archaic_modern_clades": len(clades),
                "n_passing_archaic_modern_clades": sum(as_int(row["clade_pass"]) for row in clades),
                "best_clade_lineage": clades[0].get("lineage", "") if clades else "",
                "best_clade_bootstrap": clades[0].get("bootstrap", "") if clades else "",
                "best_clade_modern_haplotype_types": clades[0].get("n_modern_haplotype_types", "") if clades else "",
                "best_clade_modern_haplotype_copies": clades[0].get("n_modern_haplotype_copies", "") if clades else "",
                "interpretation": "anchor_independent_common_haplotype_tree;not_authoritative_introgression_call",
                "tree_newick": newick,
            })

    updated_families: list[dict[str, object]] = []
    for family in family_rows:
        row: dict[str, object] = dict(family)
        best = family_best.get(family.get("family_id", ""))
        row["region_common_tree_linked"] = int(best is not None)
        row["region_common_tree_clade_pass"] = best[0] if best else 0
        row["region_common_tree_best_bootstrap"] = fmt(best[1]) if best else ""
        row["region_common_tree_lineage_filter"] = best[2] if best else ""
        row["region_common_tree_node"] = best[3] if best else ""
        updated_families.append(row)

    final = out / "final"
    write_tsv(
        final / "region_common_trees.tsv",
        [
            "locus_id", "tree_scope", "lineage_filter", "tree_status", "phy_file", "tree_file", "stats_file", "plot_file",
            "n_tree_tips", "n_bootstrap_nodes", "bootstrap_min", "bootstrap_median",
            "bootstrap_max", "n_archaic_modern_clades", "n_passing_archaic_modern_clades",
            "best_clade_lineage", "best_clade_bootstrap", "best_clade_modern_haplotype_types",
            "best_clade_modern_haplotype_copies", "interpretation", "tree_newick",
        ],
        tree_rows,
    )
    write_tsv(
        final / "region_tree_clades.unfiltered.tsv.gz",
        [
            "locus_id", "tree_scope", "lineage_filter", "clade_rank", "node", "side", "lineage", "bootstrap",
            "n_tips", "n_modern_haplotype_types", "n_modern_haplotype_copies",
            "n_archaic_tips", "modern_haplotypes", "archaic_references", "family_ids",
            "clade_pass", "rule", "interpretation",
        ],
        clade_rows,
    )
    if family_rows:
        fields = list(family_rows[0].keys())
        for extra in (
            "region_common_tree_linked", "region_common_tree_clade_pass",
            "region_common_tree_best_bootstrap", "region_common_tree_lineage_filter",
            "region_common_tree_node",
        ):
            if extra not in fields:
                fields.append(extra)
        write_tsv(final / "region_candidate_families.unfiltered.tsv.gz", fields, updated_families)
    print(
        f"PHYML REGION TREE SUMMARY: trees={len(tree_rows)} clades={len(clade_rows)} "
        f"passing={sum(as_int(row['clade_pass']) for row in clade_rows)} "
        f"output={final / 'region_common_trees.tsv'}",
        flush=True,
    )


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--scan-min-callable-sites", type=int, default=12)
    parser.add_argument("--scan-seed-prop", type=float, default=0.90)
    parser.add_argument("--scan-max-gap-bp", type=int, default=50_000)
    parser.add_argument("--full-min-callable-sites", type=int, default=10)
    parser.add_argument("--full-min-prop", type=float, default=0.80)
    parser.add_argument("--family-min-overlap", type=float, default=0.50)
    parser.add_argument("--common-tree-min-copies", type=int, default=11)
    parser.add_argument("--common-tree-min-sites", type=int, default=10)
    parser.add_argument("--min-tree-bootstrap", type=float, default=70.0)
    parser.add_argument("--min-tree-modern-types", type=int, default=2)
    parser.add_argument("--default-tree-status", choices=("not_requested", "not_run"), default="not_run")


def validate_args(args: argparse.Namespace) -> None:
    for name in ("scan_seed_prop", "full_min_prop", "family_min_overlap"):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            fail(f"--{name.replace('_', '-')} must be between 0 and 1")
    for name in (
        "scan_min_callable_sites", "scan_max_gap_bp", "full_min_callable_sites",
        "common_tree_min_copies", "common_tree_min_sites", "min_tree_modern_types",
    ):
        if getattr(args, name) < 1:
            fail(f"--{name.replace('_', '-')} must be positive")
    if not 0 <= args.min_tree_bootstrap <= 100:
        fail("--min-tree-bootstrap must be between 0 and 100")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("prepare", "summarize", "all"):
        child = sub.add_parser(action)
        add_args(child)
    args = parser.parse_args()
    validate_args(args)
    if args.action in {"prepare", "all"}:
        prepare(args)
    if args.action in {"summarize", "all"}:
        summarize(args)


if __name__ == "__main__":
    main()
