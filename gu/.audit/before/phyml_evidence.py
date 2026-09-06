#!/usr/bin/env python3
"""Evidence-aware post-processing for the GU locus/PhyML workflow.

GU now uses an explicitly layered interpretation:

0. ``region_scan`` (written by region_scan.py): every regional modern
   haplotype x archaic-reference match, independent of BED column 4;
1. ``diagnostic_candidate``: a regional match carrying a contiguous set of
   archaic-derived alleles that is rare in a selected outgroup (YRI by default);
2. ``anchor_link``: optional zoom-in at the exact BED column-4 marker; and
3. ``candidate_tree_support``: bootstrap support for a pre-specified candidate
   group clustering with the same archaic lineage in a candidate-focused tree.

The anchor is never a gate for the regional candidate. Missing or incorrect
column 4 changes only the anchor-link layer. Raw sequence identity and the
anchor-independent inventory are retained but are never interpreted alone as
an introgression call.
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

# 1000 Genomes population -> superpopulation.  This is used only when the
# supplied sample table contains population but not superpopulation.
POP_TO_SUPERPOP = {
    "ACB": "AFR", "ASW": "AFR", "ESN": "AFR", "GWD": "AFR",
    "LWK": "AFR", "MSL": "AFR", "YRI": "AFR",
    "CLM": "AMR", "MXL": "AMR", "PEL": "AMR", "PUR": "AMR",
    "CDX": "EAS", "CHB": "EAS", "CHS": "EAS", "JPT": "EAS", "KHV": "EAS",
    "CEU": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR", "TSI": "EUR",
    "BEB": "SAS", "GIH": "SAS", "ITU": "SAS", "PJL": "SAS", "STU": "SAS",
}

SAMPLE_ALIASES = ("sample", "iid", "id", "sample_id", "individual", "individual_id")
SEX_ALIASES = ("sex", "gender")
POP_ALIASES = ("population", "pop", "population_code", "pop_code")
SUPERPOP_ALIASES = (
    "superpopulation", "super_population", "superpop", "super_pop",
    "superpopulation_code", "super_pop_code", "ancestry",
)

EVIDENCE_HAP_FIELDS = [
    "sequence_similarity_pass",
    "sequence_similarity_interpretation",
    "unfiltered_region_family_ids",
    "unfiltered_region_segment_count",
    "unfiltered_region_best_segment_prop",
    "regional_stage1_inventory_present",
    "diagnostic_lineage",
    "lineage_reference_count",
    "n_lineage_diagnostic_sites",
    "n_diagnostic_compared",
    "n_diagnostic_match",
    "diagnostic_match_prop",
    "n_contiguous_diagnostic_match",
    "candidate_start",
    "candidate_end",
    "candidate_bp",
    "anchor_exact",
    "anchor_diagnostic",
    "anchor_match",
    "named_anchor_found",
    "named_anchor_id",
    "named_anchor_pos",
    "named_anchor_variant_type",
    "named_anchor_allele_index_set",
    "named_anchor_allele_set",
    "named_anchor_mixed",
    "named_anchor_lineage_consensus_allele_index",
    "named_anchor_matches_lineage",
    "named_anchor_called_copies",
    "named_anchor_lineage_matching_copies",
    "named_anchor_lineage_matching_fraction",
    "named_anchor_support_status",
    "anchor_linked_candidate_copies",
    "population_copy_counts",
    "superpopulation_copy_counts",
    "outgroup_label",
    "outgroup_copy_count",
    "outgroup_exact_haplotype_frequency",
    "outgroup_candidate_pattern_copies",
    "outgroup_candidate_pattern_frequency",
    "outgroup_haplotype_frequency",
    "nonoutgroup_copy_count",
    "nonoutgroup_candidate_pattern_copies",
    "diagnostic_candidate_pass",
    "diagnostic_candidate_reason",
]

EVIDENCE_LOCUS_FIELDS = [
    "sequence_similarity_only",
    "sequence_similarity_interpretation",
    "legacy_tree_interpretation",
    "n_unfiltered_region_families",
    "n_unfiltered_region_segments",
    "regional_stage1_status",
    "evidence_lineage",
    "evidence_lineage_reference_count",
    "named_anchor_found",
    "named_anchor_id",
    "named_anchor_pos",
    "named_anchor_variant_type",
    "named_anchor_lineage_consensus_allele_index",
    "named_anchor_candidate_match",
    "named_anchor_support_status",
    "n_anchor_linked_candidate_haplotypes",
    "n_anchor_linked_candidate_copies",
    "anchor_linked_introgression_call",
    "anchor_linked_introgression_confidence",
    "introgression_scope",
    "outgroup_label",
    "n_lineage_diagnostic_sites",
    "n_candidate_haplotypes",
    "n_candidate_copies",
    "candidate_start",
    "candidate_end",
    "candidate_bp",
    "candidate_tree_status",
    "candidate_tree_bootstrap",
    "candidate_tree_pass",
    "candidate_tree_strong",
    "multi_evidence_strong",
    "candidate_tree_candidate_purity",
    "candidate_tree_candidate_sensitivity",
    "candidate_tree_candidates_in_clade",
    "candidate_tree_archaic_in_clade",
    "introgression_call",
    "introgression_confidence",
    "independent_segment_caller_required",
    "evidence_reason",
]


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def norm_header(value: str) -> str:
    value = str(value).strip().lower().lstrip("#")
    return re.sub(r"[^a-z0-9]+", "_", value).strip("_")


def safe_slug(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._-")
    return value or "unknown"


def truth(value: object) -> bool:
    return str(value).strip().upper() in {"1", "TRUE", "T", "YES", "Y", "PASS"}


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default


def as_float(value: object, default: float | None = None) -> float | None:
    try:
        out = float(str(value).strip())
        return out if math.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def format_float(value: float | None) -> str:
    return "NA" if value is None or not math.isfinite(value) else f"{value:.12g}"


def write_text_if_changed(path: Path, text: str) -> bool:
    """Atomically write text and return True only when content changed."""
    old = path.read_text(errors="replace") if path.is_file() else None
    if old == text:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{Path('/proc/self').resolve().name}")
    tmp.write_text(text)
    tmp.replace(path)
    return True


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def remove_phyml_derivatives(phy: Path) -> None:
    candidates = [
        Path(f"{phy}_phyml_tree.txt"),
        Path(f"{phy}_phyml_stats.txt"),
        Path(f"{phy}_phyml_boot_trees.txt"),
        Path(f"{phy}_phyml_boot_stats.txt"),
        Path(f"{phy}.phyml.log"),
        Path(f"{phy}_phyml_tree.png"),
    ]
    for path in candidates:
        path.unlink(missing_ok=True)


@dataclass(frozen=True)
class SampleInfo:
    sample: str
    sex: str = "unknown"
    population: str = ""
    superpopulation: str = ""


@dataclass
class Panel:
    samples: dict[str, SampleInfo]
    requested_outgroup: tuple[str, ...]
    fallback_outgroup: tuple[str, ...]
    active_outgroup: tuple[str, ...] = ()
    active_label: str = ""
    fallback_used: bool = False

    @staticmethod
    def _matches(info: SampleInfo, groups: Sequence[str]) -> bool:
        pop = info.population.upper()
        superpop = info.superpopulation.upper()
        return any(group.upper() in {pop, superpop} for group in groups)

    def choose_outgroup(self) -> None:
        requested_n = sum(self._matches(info, self.requested_outgroup) for info in self.samples.values())
        if requested_n:
            self.active_outgroup = self.requested_outgroup
            self.active_label = ",".join(self.requested_outgroup)
            return
        fallback_n = sum(self._matches(info, self.fallback_outgroup) for info in self.samples.values())
        if fallback_n:
            self.active_outgroup = self.fallback_outgroup
            self.active_label = ",".join(self.fallback_outgroup)
            self.fallback_used = True
            return
        self.active_outgroup = ()
        self.active_label = "UNAVAILABLE"

    def group(self, sample: str) -> str:
        info = self.samples.get(sample)
        if info is None or not (info.population or info.superpopulation):
            return "UNKNOWN"
        if self.active_outgroup and self._matches(info, self.active_outgroup):
            return "OUTGROUP"
        return "NONOUTGROUP"


def split_groups(value: str) -> tuple[str, ...]:
    groups = tuple(x.upper() for x in re.split(r"[,; ]+", value.strip()) if x)
    return groups


def read_panel(path: Path | None, outgroup: str, fallback: str) -> Panel:
    requested = split_groups(outgroup)
    fallback_groups = split_groups(fallback) if fallback.upper() not in {"", "NONE", "NA"} else ()
    panel = Panel({}, requested, fallback_groups)
    if path is None or not path.is_file() or path.stat().st_size == 0:
        panel.choose_outgroup()
        return panel

    lines = [line.rstrip("\n\r") for line in path.open() if line.strip()]
    if not lines:
        panel.choose_outgroup()
        return panel
    delimiter = "\t" if "\t" in lines[0] else None
    header = lines[0].split(delimiter) if delimiter else re.split(r"\s+", lines[0].strip())
    normalized = [norm_header(x) for x in header]

    def index_of(aliases: Sequence[str]) -> int | None:
        for alias in aliases:
            alias = norm_header(alias)
            if alias in normalized:
                return normalized.index(alias)
        return None

    sample_i = index_of(SAMPLE_ALIASES)
    if sample_i is None:
        fail(f"sample table {path} has no sample/IID column")
    sex_i = index_of(SEX_ALIASES)
    pop_i = index_of(POP_ALIASES)
    super_i = index_of(SUPERPOP_ALIASES)

    for line_no, line in enumerate(lines[1:], 2):
        fields = line.split(delimiter) if delimiter else re.split(r"\s+", line.strip())
        if sample_i >= len(fields):
            continue
        sample = fields[sample_i].strip()
        if not sample:
            continue
        sex = fields[sex_i].strip().lower() if sex_i is not None and sex_i < len(fields) else "unknown"
        pop = fields[pop_i].strip().upper() if pop_i is not None and pop_i < len(fields) else ""
        superpop = fields[super_i].strip().upper() if super_i is not None and super_i < len(fields) else ""
        if not superpop and pop:
            superpop = POP_TO_SUPERPOP.get(pop, "")
        panel.samples[sample] = SampleInfo(sample, sex, pop, superpop)
    panel.choose_outgroup()
    return panel


def parse_copies(value: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for token in str(value).split(";"):
        token = token.strip()
        if not token:
            continue
        if ":" in token:
            sample, hap = token.rsplit(":", 1)
        else:
            sample, hap = token, ""
        if sample:
            out.append((sample, hap))
    return out


def lineage_name(value: str) -> str:
    text = str(value).strip()
    if re.search(r"neander|altai|vindija|chagyr", text, re.I):
        return "Neanderthal"
    if re.search(r"denis", text, re.I):
        return "Denisovan"
    return text or "Archaic"


def locus_dirs(out: Path) -> list[Path]:
    root = out / "loci"
    if (root / "sites.tsv").is_file():
        return [root]
    return sorted(path for path in root.iterdir() if path.is_dir() and (path / "sites.tsv").is_file()) if root.is_dir() else []


def final_rows_by_locus(out: Path, filename: str) -> dict[str, list[dict[str, str]]]:
    rows = read_tsv(out / "final" / filename)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row.get("locus_id", "")].append(row)
    return grouped


def infer_locus_id(locus_dir: Path, locus_meta: dict[str, list[dict[str, str]]]) -> str:
    if locus_dir.name != "loci":
        return locus_dir.name
    nonempty = [key for key in locus_meta if key]
    if len(nonempty) == 1:
        return nonempty[0]
    selected = read_tsv(locus_dir / "selected_region.tsv")
    ids = {row.get("locus_id", "") for row in selected if row.get("locus_id", "")}
    if len(ids) == 1:
        return next(iter(ids))
    fail(f"cannot infer flat locus_id under {locus_dir}")


def validate_sequences(locus_id: str, n_sites: int, named: dict[str, str]) -> None:
    bad = {name: len(seq) for name, seq in named.items() if len(seq) != n_sites}
    if bad:
        detail = ", ".join(f"{name}={length}" for name, length in sorted(bad.items()))
        fail(f"{locus_id}: sequence/site length mismatch (sites={n_sites}; {detail})")


def lineage_consensus(
    archaic_rows: list[dict[str, str]],
    lineage: str,
    n_sites: int,
    min_refs: int,
) -> tuple[str, list[int], list[list[str]]]:
    refs = [row for row in archaic_rows if lineage_name(row.get("lineage") or row.get("archaic", "")) == lineage]
    sequences = [row.get("seq", "").upper() for row in refs]
    consensus: list[str] = []
    agreeing_counts: list[int] = []
    agreeing_refs: list[list[str]] = []
    for i in range(n_sites):
        called = [(row.get("archaic", ""), seq[i]) for row, seq in zip(refs, sequences) if seq[i] in BASES]
        counts = Counter(base for _, base in called)
        if not counts:
            consensus.append("N")
            agreeing_counts.append(0)
            agreeing_refs.append([])
            continue
        top_n = max(counts.values())
        tops = sorted(base for base, count in counts.items() if count == top_n)
        if len(tops) != 1 or top_n < min_refs:
            consensus.append("N")
            agreeing_counts.append(top_n)
            agreeing_refs.append([])
            continue
        base = tops[0]
        consensus.append(base)
        agreeing_counts.append(top_n)
        agreeing_refs.append([name for name, observed in called if observed == base])
    return "".join(consensus), agreeing_counts, agreeing_refs


def best_contiguous_cluster(positions: list[int], max_gap: int) -> list[int]:
    if not positions:
        return []
    positions = sorted(set(positions))
    clusters: list[list[int]] = [[positions[0]]]
    for pos in positions[1:]:
        if pos - clusters[-1][-1] > max_gap:
            clusters.append([pos])
        else:
            clusters[-1].append(pos)
    return max(clusters, key=lambda cluster: (len(cluster), cluster[-1] - cluster[0], -cluster[0]))


def hamming_similarity(seq: str, ref: str, indexes: Sequence[int]) -> tuple[int, int, float | None]:
    compared = matched = 0
    for i in indexes:
        if i >= len(seq) or i >= len(ref):
            continue
        a, b = seq[i], ref[i]
        if a not in BASES or b not in BASES:
            continue
        compared += 1
        matched += int(a == b)
    return compared, matched, matched / compared if compared else None


def metadata_value(meta: dict[str, str], key: str, default: str = "") -> str:
    return str(meta.get(key, default))


def merge_fields(path: Path, key_fields: Sequence[str], additions: dict[tuple[str, ...], dict[str, object]], fields: Sequence[str]) -> None:
    rows = read_tsv(path)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    for field_name in fields:
        if field_name not in fieldnames:
            fieldnames.append(field_name)
    for row in rows:
        key = tuple(row.get(field, "") for field in key_fields)
        extra = additions.get(key, {})
        for field_name in fields:
            if field_name in extra:
                row[field_name] = extra[field_name]
            elif field_name not in row:
                row[field_name] = ""
    write_tsv(path, fieldnames, rows)


def prepare(args: argparse.Namespace) -> None:
    out: Path = args.out
    panel = read_panel(args.sample_file, args.outgroup, args.outgroup_fallback)
    locus_meta_map = final_rows_by_locus(out, "loci.tsv")
    final_hap_map = final_rows_by_locus(out, "haplotypes.tsv")
    anchor_summary_map = final_rows_by_locus(out, "anchors.tsv")
    anchor_hap_map = final_rows_by_locus(out, "anchor_haplotypes.tsv")
    anchor_copy_map = final_rows_by_locus(out, "anchor_copies.tsv")
    region_segment_rows = read_tsv(out / "final" / "region_match_segments.unfiltered.tsv.gz")
    region_family_rows = read_tsv(out / "final" / "region_candidate_families.unfiltered.tsv.gz")
    region_locus_rows = final_rows_by_locus(out, "region_scan_loci.unfiltered.tsv.gz")
    region_segments_by_key: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for region_row in region_segment_rows:
        region_segments_by_key[(
            region_row.get("locus_id", ""),
            region_row.get("hap_id", ""),
            lineage_name(region_row.get("lineage", "")),
        )].append(region_row)
    region_families_by_locus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for region_row in region_family_rows:
        region_families_by_locus[region_row.get("locus_id", "")].append(region_row)
    dirs = locus_dirs(out)
    if not dirs and not locus_meta_map:
        fail(f"no PhyML locus summary or artifacts found under {out}")

    evidence_site_rows: list[dict[str, object]] = []
    evidence_hap_rows: list[dict[str, object]] = []
    evidence_locus_rows: list[dict[str, object]] = []
    legacy_hap_additions: dict[tuple[str, ...], dict[str, object]] = {}
    legacy_locus_additions: dict[tuple[str, ...], dict[str, object]] = {}
    desired_phy_files: set[Path] = set()

    represented_loci = {infer_locus_id(path, locus_meta_map) for path in dirs}
    for missing_locus_id, meta_rows in locus_meta_map.items():
        if not missing_locus_id or missing_locus_id in represented_loci:
            continue
        meta = meta_rows[0]
        status = meta.get("status", "") or "missing_locus_artifacts"
        row = {
            "locus_id": missing_locus_id,
            "genome_build": meta.get("genome_build", ""),
            "chr": meta.get("chr", ""),
            "core_start": meta.get("core_start", ""),
            "core_end": meta.get("core_end", ""),
            "name": meta.get("name", missing_locus_id),
            "sequence_similarity_only": int(truth(meta.get("direct_match_pass", 0))),
            "sequence_similarity_interpretation": "legacy_direct_match_is_raw_acgt_identity_only",
            "legacy_tree_interpretation": "all_unique_haplotype_tree_is_exploratory_not_authoritative",
            "n_unfiltered_region_families": len(region_families_by_locus.get(missing_locus_id, [])),
            "n_unfiltered_region_segments": sum(1 for x in region_segment_rows if x.get("locus_id") == missing_locus_id),
            "regional_stage1_status": (region_locus_rows.get(missing_locus_id, [{}])[0].get("stage1_status", "missing")),
            "evidence_lineage": "",
            "evidence_lineage_reference_count": 0,
            "named_anchor_found": 0,
            "named_anchor_id": meta.get("name", missing_locus_id),
            "named_anchor_pos": "",
            "named_anchor_variant_type": "",
            "named_anchor_lineage_consensus_allele_index": "",
            "named_anchor_candidate_match": "",
            "named_anchor_support_status": "anchor_not_evaluable",
            "n_anchor_linked_candidate_haplotypes": 0,
            "n_anchor_linked_candidate_copies": 0,
            "anchor_linked_introgression_call": "not_evaluable",
            "anchor_linked_introgression_confidence": "NA",
            "introgression_scope": "not_evaluable",
            "outgroup_label": panel.active_label,
            "n_lineage_diagnostic_sites": 0,
            "n_candidate_haplotypes": 0,
            "n_candidate_copies": 0,
            "candidate_start": "",
            "candidate_end": "",
            "candidate_bp": "",
            "candidate_tree_status": "not_run",
            "candidate_tree_bootstrap": "",
            "candidate_tree_pass": 0,
            "candidate_tree_strong": 0,
            "multi_evidence_strong": 0,
            "candidate_tree_candidate_purity": "",
            "candidate_tree_candidate_sensitivity": "",
            "candidate_tree_candidates_in_clade": "",
            "candidate_tree_archaic_in_clade": "",
            "introgression_call": "not_evaluable",
            "introgression_confidence": "NA",
            "independent_segment_caller_required": 1,
            "evidence_reason": f"no_sequence_artifact:{status}",
            "outgroup_fallback_used": int(panel.fallback_used),
            "n_unknown_sample_copies": 0,
            "unknown_samples_first10": "",
        }
        evidence_locus_rows.append(row)
        legacy_locus_additions[(missing_locus_id,)] = {key: row.get(key, "") for key in EVIDENCE_LOCUS_FIELDS}

    for locus_dir in dirs:
        locus_id = infer_locus_id(locus_dir, locus_meta_map)
        meta_rows = locus_meta_map.get(locus_id, [])
        meta = meta_rows[0] if meta_rows else {"locus_id": locus_id, "name": locus_id}
        anchor_summary_rows = anchor_summary_map.get(locus_id, [])
        anchor_summary = anchor_summary_rows[0] if anchor_summary_rows else {}
        anchor_haps = {row.get("hap_id", ""): row for row in anchor_hap_map.get(locus_id, [])}
        anchor_state_counts: dict[str, Counter[str]] = defaultdict(Counter)
        for anchor_copy in anchor_copy_map.get(locus_id, []):
            hap_id = anchor_copy.get("hap_id", "")
            allele_index = anchor_copy.get("allele_index", "")
            if hap_id and allele_index not in {"", "NA"}:
                anchor_state_counts[hap_id][str(allele_index)] += 1
        sites = read_tsv(locus_dir / "sites.tsv")
        archaic_rows = read_tsv(locus_dir / "archaic.tsv")
        ancestral_rows = read_tsv(locus_dir / "ancestral.tsv")
        hap_rows = read_tsv(locus_dir / "haplotypes.tsv")
        if not sites or not archaic_rows or not hap_rows:
            evidence_locus_rows.append({
                "locus_id": locus_id,
                "evidence_lineage": "",
                "outgroup_label": panel.active_label,
                "n_lineage_diagnostic_sites": 0,
                "n_candidate_haplotypes": 0,
                "n_candidate_copies": 0,
                "introgression_call": "not_evaluable",
                "introgression_confidence": "NA",
                "evidence_reason": "missing_sites_archaic_or_haplotype_artifact",
            })
            continue

        n_sites = len(sites)
        ancestral = (ancestral_rows[0].get("seq", "") if ancestral_rows else "N" * n_sites).upper()
        archaic_named = {row.get("archaic", f"archaic_{i}"): row.get("seq", "").upper() for i, row in enumerate(archaic_rows)}
        hap_named = {row.get("hap_id", f"H{i}"): row.get("seq", "").upper() for i, row in enumerate(hap_rows)}
        validate_sequences(locus_id, n_sites, {"Ancestral": ancestral, **archaic_named, **hap_named})

        final_haps = {row.get("hap_id", ""): row for row in final_hap_map.get(locus_id, [])}
        copy_info: dict[str, list[tuple[str, str]]] = {}
        hap_group_counts: dict[str, Counter[str]] = {}
        hap_population_counts: dict[str, Counter[str]] = {}
        hap_superpopulation_counts: dict[str, Counter[str]] = {}
        total_group_counts: Counter[str] = Counter()
        unknown_samples: set[str] = set()
        for row in hap_rows:
            hap_id = row.get("hap_id", "")
            copies = parse_copies(row.get("copies", ""))
            copy_info[hap_id] = copies
            groups: Counter[str] = Counter()
            pop_counts: Counter[str] = Counter()
            super_counts: Counter[str] = Counter()
            for sample, _ in copies:
                group = panel.group(sample)
                groups[group] += 1
                total_group_counts[group] += 1
                info = panel.samples.get(sample)
                if info is not None:
                    pop_counts[info.population or "UNKNOWN"] += 1
                    super_counts[info.superpopulation or "UNKNOWN"] += 1
                else:
                    pop_counts["UNKNOWN"] += 1
                    super_counts["UNKNOWN"] += 1
                if group == "UNKNOWN":
                    unknown_samples.add(sample)
            # If copies are absent, retain the total count but population
            # evidence remains unavailable for those unassigned copies.
            declared_n = as_int(row.get("n"), len(copies))
            if declared_n > len(copies):
                groups["UNKNOWN"] += declared_n - len(copies)
                total_group_counts["UNKNOWN"] += declared_n - len(copies)
            hap_group_counts[hap_id] = groups
            hap_population_counts[hap_id] = pop_counts
            hap_superpopulation_counts[hap_id] = super_counts

        # Per-site modern allele counts by outgroup status.
        site_counts: list[dict[str, Counter[str]]] = []
        for i in range(n_sites):
            counts = {
                "OUTGROUP": Counter(),
                "NONOUTGROUP": Counter(),
                "UNKNOWN": Counter(),
                "ALL": Counter(),
            }
            for row in hap_rows:
                hap_id = row.get("hap_id", "")
                seq = row.get("seq", "").upper()
                base = seq[i]
                if base not in BASES:
                    continue
                for group, n in hap_group_counts[hap_id].items():
                    counts[group][base] += n
                    counts["ALL"][base] += n
            site_counts.append(counts)

        lineages = sorted({lineage_name(row.get("lineage") or row.get("archaic", "")) for row in archaic_rows})
        lineage_payload: dict[str, dict[str, object]] = {}
        for lineage in lineages:
            refs = [row for row in archaic_rows if lineage_name(row.get("lineage") or row.get("archaic", "")) == lineage]
            min_refs = args.min_neanderthal_refs if lineage == "Neanderthal" else args.min_other_lineage_refs
            consensus, agreeing_counts, agreeing_refs = lineage_consensus(archaic_rows, lineage, n_sites, min_refs)
            diagnostic_indexes: list[int] = []
            for i, (site, archaic_base, ancestral_base) in enumerate(zip(sites, consensus, ancestral)):
                out_counts = site_counts[i]["OUTGROUP"]
                nonout_counts = site_counts[i]["NONOUTGROUP"]
                out_called = sum(out_counts.values())
                nonout_called = sum(nonout_counts.values())
                out_freq = out_counts[archaic_base] / out_called if archaic_base in BASES and out_called else None
                nonout_freq = nonout_counts[archaic_base] / nonout_called if archaic_base in BASES and nonout_called else None
                ancestral_diff = archaic_base in BASES and ancestral_base in BASES and archaic_base != ancestral_base
                outgroup_covered = bool(panel.active_outgroup) and out_called >= args.min_outgroup_copies
                diagnostic_pass = bool(
                    ancestral_diff
                    and outgroup_covered
                    and out_freq is not None
                    and out_freq <= args.max_outgroup_allele_freq
                    and nonout_counts[archaic_base] >= 1
                )
                if diagnostic_pass:
                    diagnostic_indexes.append(i)
                evidence_site_rows.append({
                    "locus_id": locus_id,
                    "chr": site.get("chr", meta.get("chr", "")),
                    "pos": site.get("pos", ""),
                    "id": site.get("id", ""),
                    "ref": site.get("ref", ""),
                    "alt": site.get("alt", ""),
                    "lineage": lineage,
                    "lineage_consensus": archaic_base,
                    "ancestral": ancestral_base,
                    "lineage_reference_count": len(refs),
                    "lineage_agreeing_reference_count": agreeing_counts[i],
                    "lineage_agreeing_references": ",".join(agreeing_refs[i]),
                    "outgroup_label": panel.active_label,
                    "outgroup_called_copies": out_called,
                    "outgroup_archaic_allele_copies": out_counts[archaic_base] if archaic_base in BASES else 0,
                    "outgroup_archaic_allele_frequency": format_float(out_freq),
                    "nonoutgroup_called_copies": nonout_called,
                    "nonoutgroup_archaic_allele_copies": nonout_counts[archaic_base] if archaic_base in BASES else 0,
                    "nonoutgroup_archaic_allele_frequency": format_float(nonout_freq),
                    "archaic_differs_from_ancestral": int(ancestral_diff),
                    "diagnostic_site_pass": int(diagnostic_pass),
                    "diagnostic_site_reason": (
                        "pass" if diagnostic_pass
                        else "outgroup_metadata_unavailable" if not panel.active_outgroup
                        else "ancestral_unavailable_or_same" if not ancestral_diff
                        else "outgroup_coverage_low" if not outgroup_covered
                        else "outgroup_frequency_high" if out_freq is not None and out_freq > args.max_outgroup_allele_freq
                        else "archaic_allele_absent_nonoutgroup"
                    ),
                })
            lineage_payload[lineage] = {
                "refs": refs,
                "consensus": consensus,
                "diagnostic_indexes": diagnostic_indexes,
                "reference_count": len(refs),
            }

        # Exact named anchor, if it survived SNP/MAC filtering.  We do not treat
        # the nearest eligible SNP as the named anchor for evidence purposes.
        anchor_names = {locus_id, meta.get("name", "")}
        anchor_indexes = [i for i, site in enumerate(sites) if site.get("id", "") in anchor_names]
        exact_anchor_i = anchor_indexes[0] if anchor_indexes else None

        # Precompute which diagnostic indexes each H type matches.  Candidate
        # motif frequencies are then cached by index set, avoiding repeated
        # character-by-character H x H scans at polymorphic loci.
        lineage_match_sets: dict[str, dict[str, frozenset[int]]] = {}
        pattern_count_cache: dict[str, dict[frozenset[int], tuple[int, int]]] = defaultdict(dict)
        for lineage, payload in lineage_payload.items():
            consensus = str(payload["consensus"])
            diagnostic_indexes = list(payload["diagnostic_indexes"])
            lineage_match_sets[lineage] = {
                row.get("hap_id", ""): frozenset(
                    i for i in diagnostic_indexes
                    if row.get("seq", "").upper()[i] in BASES
                    and row.get("seq", "").upper()[i] == consensus[i]
                )
                for row in hap_rows
            }

        scored_by_lineage: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in hap_rows:
            hap_id = row.get("hap_id", "")
            seq = row.get("seq", "").upper()
            declared_n = as_int(row.get("n"), len(copy_info.get(hap_id, [])))
            groups = hap_group_counts[hap_id]
            out_n = groups["OUTGROUP"]
            nonout_n = groups["NONOUTGROUP"]
            out_total = total_group_counts["OUTGROUP"]
            out_exact_hap_freq = out_n / out_total if out_total else None
            raw = final_haps.get(hap_id, {})
            raw_similarity_pass = truth(raw.get("direct_match_pass", 0))
            named_anchor = anchor_haps.get(hap_id, {})
            population_copy_counts = ",".join(
                f"{key}:{value}" for key, value in sorted(hap_population_counts[hap_id].items())
            )
            superpopulation_copy_counts = ",".join(
                f"{key}:{value}" for key, value in sorted(hap_superpopulation_counts[hap_id].items())
            )

            lineage_scores: list[dict[str, object]] = []
            for lineage, payload in lineage_payload.items():
                consensus = str(payload["consensus"])
                diagnostic_indexes = list(payload["diagnostic_indexes"])
                regional_segments = region_segments_by_key.get((locus_id, hap_id, lineage), [])
                regional_family_ids = sorted({row.get("family_id", "") for row in regional_segments if row.get("family_id", "")})
                regional_best_prop = max(
                    (as_float(row.get("prop_match"), 0.0) or 0.0 for row in regional_segments),
                    default=None,
                )
                compared, matched, prop = hamming_similarity(seq, consensus, diagnostic_indexes)
                own_match_set = lineage_match_sets[lineage].get(hap_id, frozenset())
                matched_positions = [as_int(sites[i].get("pos")) for i in own_match_set]
                cluster = best_contiguous_cluster(matched_positions, args.max_diagnostic_gap_bp)
                cluster_n = len(cluster)
                cluster_set = set(cluster)
                cluster_indexes = frozenset(
                    i for i in own_match_set if as_int(sites[i].get("pos")) in cluster_set
                )
                candidate_start = cluster[0] - 1 if cluster else ""
                candidate_end = cluster[-1] if cluster else ""
                candidate_bp = cluster[-1] - cluster[0] + 1 if cluster else ""
                # Frequency of the archaic-derived diagnostic motif, not merely
                # the exact full-length H identifier.  The latter can be absent
                # in YRI even when the same diagnostic motif is common there.
                pattern_out_n = pattern_nonout_n = 0
                if cluster_indexes:
                    cached = pattern_count_cache[lineage].get(cluster_indexes)
                    if cached is None:
                        pattern_out_n = sum(
                            hap_group_counts[other_id]["OUTGROUP"]
                            for other_id, other_matches in lineage_match_sets[lineage].items()
                            if cluster_indexes.issubset(other_matches)
                        )
                        pattern_nonout_n = sum(
                            hap_group_counts[other_id]["NONOUTGROUP"]
                            for other_id, other_matches in lineage_match_sets[lineage].items()
                            if cluster_indexes.issubset(other_matches)
                        )
                        pattern_count_cache[lineage][cluster_indexes] = (pattern_out_n, pattern_nonout_n)
                    else:
                        pattern_out_n, pattern_nonout_n = cached
                out_pattern_freq = pattern_out_n / out_total if cluster_indexes and out_total else None
                anchor_diagnostic = int(exact_anchor_i is not None and exact_anchor_i in diagnostic_indexes)
                anchor_match = (
                    int(seq[exact_anchor_i] == consensus[exact_anchor_i])
                    if exact_anchor_i is not None and seq[exact_anchor_i] in BASES and consensus[exact_anchor_i] in BASES
                    else ""
                )
                named_consensus_field = (
                    "neanderthal_consensus_allele_index"
                    if lineage == "Neanderthal"
                    else "denisovan_consensus_allele_index"
                    if lineage == "Denisovan"
                    else ""
                )
                named_consensus_index = (
                    anchor_summary.get(named_consensus_field, "") if named_consensus_field else ""
                )
                named_index_set = named_anchor.get("allele_index_set", "")
                named_mixed = as_int(named_anchor.get("mixed_anchor_within_snp_haplotype", 0))
                named_counts = anchor_state_counts.get(hap_id, Counter())
                named_called_copies = sum(named_counts.values())
                named_matching_copies = named_counts.get(str(named_consensus_index), 0) if named_consensus_index else 0
                named_matching_fraction = (
                    named_matching_copies / named_called_copies if named_called_copies and named_consensus_index else None
                )
                named_match = (
                    int(named_called_copies > 0 and named_matching_copies == named_called_copies)
                    if named_consensus_index else ""
                )
                if as_int(anchor_summary.get("anchor_provided", 1)) == 0:
                    named_support_status = "named_anchor_not_provided"
                elif not as_int(anchor_summary.get("exact_anchor_found", 0)):
                    named_support_status = "named_anchor_not_found"
                elif not named_consensus_index:
                    named_support_status = "lineage_anchor_uninformative_or_representation_mismatch"
                elif named_called_copies == 0:
                    named_support_status = "named_anchor_not_phase_assigned"
                elif named_matching_copies == named_called_copies:
                    named_support_status = "all_copies_match_lineage_anchor"
                elif named_matching_copies > 0:
                    named_support_status = "some_copies_match_lineage_anchor"
                else:
                    named_support_status = "no_copy_matches_lineage_anchor"
                reasons: list[str] = []
                if not panel.active_outgroup:
                    reasons.append("outgroup_metadata_unavailable")
                if len(diagnostic_indexes) < args.min_diagnostic_sites:
                    reasons.append("too_few_lineage_diagnostic_sites")
                if compared < args.min_diagnostic_sites:
                    reasons.append("too_few_diagnostic_sites_compared")
                if matched < args.min_diagnostic_sites:
                    reasons.append("too_few_diagnostic_matches")
                if cluster_n < args.min_diagnostic_sites:
                    reasons.append("no_contiguous_diagnostic_tract")
                if prop is None or prop < args.min_diagnostic_match_prop:
                    reasons.append("diagnostic_match_proportion_low")
                if declared_n < args.min_candidate_copies:
                    reasons.append("candidate_copy_count_low")
                if nonout_n < args.min_candidate_copies:
                    reasons.append("nonoutgroup_copy_count_low")
                if out_pattern_freq is None:
                    reasons.append("outgroup_candidate_pattern_frequency_unavailable")
                elif out_pattern_freq > args.max_outgroup_haplotype_freq:
                    reasons.append("candidate_diagnostic_pattern_present_in_outgroup")
                candidate_pass = not reasons
                anchor_linked_candidate_copies = named_matching_copies if candidate_pass else 0
                score = {
                    "locus_id": locus_id,
                    "hap_id": hap_id,
                    "n": declared_n,
                    "copies": row.get("copies", ""),
                    "raw_best_archaic": raw.get("best_archaic", ""),
                    "raw_best_lineage": raw.get("best_lineage", ""),
                    "raw_n_compared": raw.get("n_compared", ""),
                    "raw_n_match": raw.get("n_match", ""),
                    "raw_prop_match": raw.get("prop_match", ""),
                    "sequence_similarity_pass": int(raw_similarity_pass),
                    "sequence_similarity_interpretation": "raw_acgt_identity_not_an_introgression_call",
                    "unfiltered_region_family_ids": ",".join(regional_family_ids),
                    "unfiltered_region_segment_count": len(regional_segments),
                    "unfiltered_region_best_segment_prop": format_float(regional_best_prop),
                    "regional_stage1_inventory_present": int(locus_id in region_locus_rows or bool(regional_segments)),
                    "diagnostic_lineage": lineage,
                    "lineage_reference_count": payload["reference_count"],
                    "n_lineage_diagnostic_sites": len(diagnostic_indexes),
                    "n_diagnostic_compared": compared,
                    "n_diagnostic_match": matched,
                    "diagnostic_match_prop": format_float(prop),
                    "n_contiguous_diagnostic_match": cluster_n,
                    "candidate_start": candidate_start,
                    "candidate_end": candidate_end,
                    "candidate_bp": candidate_bp,
                    "anchor_exact": int(exact_anchor_i is not None),
                    "anchor_diagnostic": anchor_diagnostic,
                    "anchor_match": anchor_match,
                    "named_anchor_found": as_int(anchor_summary.get("exact_anchor_found", 0)),
                    "named_anchor_id": anchor_summary.get("anchor_id", meta.get("name", locus_id)),
                    "named_anchor_pos": anchor_summary.get("pos", ""),
                    "named_anchor_variant_type": anchor_summary.get("variant_type", ""),
                    "named_anchor_allele_index_set": named_index_set,
                    "named_anchor_allele_set": named_anchor.get("allele_set", ""),
                    "named_anchor_mixed": named_mixed,
                    "named_anchor_lineage_consensus_allele_index": named_consensus_index,
                    "named_anchor_matches_lineage": named_match,
                    "named_anchor_called_copies": named_called_copies,
                    "named_anchor_lineage_matching_copies": named_matching_copies,
                    "named_anchor_lineage_matching_fraction": format_float(named_matching_fraction),
                    "named_anchor_support_status": named_support_status,
                    "anchor_linked_candidate_copies": anchor_linked_candidate_copies,
                    "population_copy_counts": population_copy_counts,
                    "superpopulation_copy_counts": superpopulation_copy_counts,
                    "outgroup_label": panel.active_label,
                    "outgroup_copy_count": out_n,
                    "outgroup_exact_haplotype_frequency": format_float(out_exact_hap_freq),
                    "outgroup_candidate_pattern_copies": pattern_out_n,
                    "outgroup_candidate_pattern_frequency": format_float(out_pattern_freq),
                    # Backward-compatible alias; now intentionally refers to
                    # the diagnostic-pattern frequency used for the decision.
                    "outgroup_haplotype_frequency": format_float(out_pattern_freq),
                    "nonoutgroup_copy_count": nonout_n,
                    "nonoutgroup_candidate_pattern_copies": pattern_nonout_n,
                    "diagnostic_candidate_pass": int(candidate_pass),
                    "diagnostic_candidate_reason": "pass" if candidate_pass else ";".join(reasons),
                    "seq": seq,
                }
                lineage_scores.append(score)
                scored_by_lineage[lineage].append(score)

            # One authoritative lineage per modern haplotype.  Retain all
            # lineage rows in evidence_haplotypes.tsv, but annotate the legacy
            # table with the best candidate-aware lineage only.
            def score_key(item: dict[str, object]) -> tuple:
                prop = as_float(item["diagnostic_match_prop"], -1.0) or -1.0
                return (
                    as_int(item["diagnostic_candidate_pass"]),
                    as_int(item["n_contiguous_diagnostic_match"]),
                    prop,
                    as_int(item["n_diagnostic_match"]),
                    as_int(item["lineage_reference_count"]),
                    str(item["diagnostic_lineage"]),
                )

            best_score = max(lineage_scores, key=score_key) if lineage_scores else None
            for item in lineage_scores:
                evidence_hap_rows.append(item)
            if best_score:
                legacy_hap_additions[(locus_id, hap_id)] = {key: best_score.get(key, "") for key in EVIDENCE_HAP_FIELDS}

        locus_candidates: list[dict[str, object]] = []
        for lineage, scores in scored_by_lineage.items():
            candidates = [row for row in scores if as_int(row["diagnostic_candidate_pass"]) == 1]
            if not candidates:
                continue
            candidate_types = len(candidates)
            candidate_copies = sum(as_int(row["n"]) for row in candidates)
            anchor_linked_scores = [row for row in candidates if as_int(row.get("anchor_linked_candidate_copies")) > 0]
            anchor_linked_copies = sum(as_int(row.get("anchor_linked_candidate_copies")) for row in candidates)
            starts = [as_int(row["candidate_start"], -1) for row in candidates if str(row["candidate_start"]) not in {"", "NA"}]
            ends = [as_int(row["candidate_end"], -1) for row in candidates if str(row["candidate_end"]) not in {"", "NA"}]
            locus_candidates.append({
                "lineage": lineage,
                "scores": candidates,
                "candidate_types": candidate_types,
                "candidate_copies": candidate_copies,
                # Number of candidate copies for which the exact named marker
                # could be assigned to a phased chromosome.  Keep this
                # distinct from copies that carry the lineage-consensus anchor
                # allele: an indel can be present in the VCF yet unassignable
                # to a SNP-defined H type when phase/representation differs.
                "anchor_called_copies": sum(as_int(row.get("named_anchor_called_copies")) for row in candidates),
                "anchor_linked_types": len(anchor_linked_scores),
                "anchor_linked_copies": anchor_linked_copies,
                "anchor_mixed_types": sum(as_int(row.get("named_anchor_mixed")) for row in candidates),
                "start": min(starts) if starts else "",
                "end": max(ends) if ends else "",
                "bp": max(ends) - min(starts) if starts and ends else "",
                "n_diagnostic_sites": max(as_int(row["n_lineage_diagnostic_sites"]) for row in candidates),
                "lineage_reference_count": as_int(lineage_payload[lineage]["reference_count"]),
            })

        best_locus_candidate = max(
            locus_candidates,
            key=lambda item: (
                as_int(item["candidate_copies"]),
                as_int(item.get("anchor_linked_copies")),
                as_int(item["candidate_types"]),
                as_int(item["n_diagnostic_sites"]),
                str(item["lineage"]),
            ),
            default=None,
        )

        if not panel.active_outgroup:
            prep_call = "not_evaluable"
            prep_confidence = "NA"
            prep_reason = "outgroup_population_metadata_unavailable"
        elif not any(lineage_payload[lineage]["diagnostic_indexes"] for lineage in lineage_payload):
            prep_call = "not_evaluable"
            prep_confidence = "NA"
            prep_reason = "no_callable_archaic_derived_sites_rare_in_outgroup"
        elif best_locus_candidate is None:
            prep_call = "not_supported"
            prep_confidence = "none"
            prep_reason = "no_haplotype_passed_diagnostic_sequence_criteria"
        else:
            prep_call = "candidate_sequence_only"
            prep_confidence = "low"
            prep_reason = "diagnostic_candidate_defined;candidate_tree_pending"

        if as_int(anchor_summary.get("anchor_provided", 1)) == 0:
            locus_anchor_status = "named_anchor_not_provided"
        elif not as_int(anchor_summary.get("exact_anchor_found", 0)):
            locus_anchor_status = "named_anchor_not_found"
        elif not best_locus_candidate:
            locus_anchor_status = "no_regional_candidate"
        elif not str(
            anchor_summary.get(
                "neanderthal_consensus_allele_index"
                if best_locus_candidate["lineage"] == "Neanderthal"
                else "denisovan_consensus_allele_index",
                "",
            )
        ).strip():
            locus_anchor_status = "lineage_anchor_uninformative_or_representation_mismatch"
        elif as_int(best_locus_candidate.get("anchor_called_copies")) == 0:
            locus_anchor_status = "named_anchor_not_phase_assigned"
        elif as_int(best_locus_candidate.get("anchor_linked_copies")) > 0:
            locus_anchor_status = "candidate_copies_match_lineage_anchor"
        elif as_int(best_locus_candidate.get("anchor_mixed_types")) > 0:
            locus_anchor_status = "candidate_haplotype_mixed_at_named_anchor"
        else:
            locus_anchor_status = "candidate_does_not_match_lineage_anchor"

        raw_direct = truth(meta.get("direct_match_pass", 0))
        locus_row: dict[str, object] = {
            "locus_id": locus_id,
            "genome_build": meta.get("genome_build", ""),
            "chr": meta.get("chr", sites[0].get("chr", "")),
            "core_start": meta.get("core_start", ""),
            "core_end": meta.get("core_end", ""),
            "name": meta.get("name", locus_id),
            "sequence_similarity_only": int(raw_direct),
            "sequence_similarity_interpretation": "legacy_direct_match_is_raw_acgt_identity_only",
            "legacy_tree_interpretation": "all_unique_haplotype_tree_is_exploratory_not_authoritative",
            "n_unfiltered_region_families": len(region_families_by_locus.get(locus_id, [])),
            "n_unfiltered_region_segments": sum(1 for x in region_segment_rows if x.get("locus_id") == locus_id),
            "regional_stage1_status": (region_locus_rows.get(locus_id, [{}])[0].get("stage1_status", "missing")),
            "evidence_lineage": best_locus_candidate["lineage"] if best_locus_candidate else "",
            "evidence_lineage_reference_count": best_locus_candidate["lineage_reference_count"] if best_locus_candidate else 0,
            "named_anchor_found": as_int(anchor_summary.get("exact_anchor_found", 0)),
            "named_anchor_id": anchor_summary.get("anchor_id", meta.get("name", locus_id)),
            "named_anchor_pos": anchor_summary.get("pos", ""),
            "named_anchor_variant_type": anchor_summary.get("variant_type", ""),
            "named_anchor_lineage_consensus_allele_index": (
                anchor_summary.get(
                    "neanderthal_consensus_allele_index"
                    if best_locus_candidate and best_locus_candidate["lineage"] == "Neanderthal"
                    else "denisovan_consensus_allele_index",
                    "",
                ) if best_locus_candidate else ""
            ),
            "named_anchor_candidate_match": (
                int(as_int(best_locus_candidate.get("anchor_linked_copies")) > 0) if best_locus_candidate else ""
            ),
            "named_anchor_support_status": locus_anchor_status,
            "n_anchor_linked_candidate_haplotypes": best_locus_candidate.get("anchor_linked_types", 0) if best_locus_candidate else 0,
            "n_anchor_linked_candidate_copies": best_locus_candidate.get("anchor_linked_copies", 0) if best_locus_candidate else 0,
            "anchor_linked_introgression_call": (
                "candidate_sequence_only" if best_locus_candidate and as_int(best_locus_candidate.get("anchor_linked_copies")) > 0
                else "not_evaluable" if locus_anchor_status in {
                    "named_anchor_not_provided",
                    "named_anchor_not_found",
                    "named_anchor_not_phase_assigned",
                    "lineage_anchor_uninformative_or_representation_mismatch",
                }
                else "not_supported"
            ),
            "anchor_linked_introgression_confidence": (
                "low" if best_locus_candidate and as_int(best_locus_candidate.get("anchor_linked_copies")) > 0
                else "NA" if locus_anchor_status in {
                    "named_anchor_not_provided",
                    "named_anchor_not_found",
                    "named_anchor_not_phase_assigned",
                    "lineage_anchor_uninformative_or_representation_mismatch",
                }
                else "none"
            ),
            "introgression_scope": (
                "named_anchor_linked_region" if best_locus_candidate and as_int(best_locus_candidate.get("anchor_linked_copies")) > 0
                else "regional_candidate_anchor_uninformative" if best_locus_candidate and locus_anchor_status in {
                    "named_anchor_not_provided",
                    "named_anchor_not_found",
                    "named_anchor_not_phase_assigned",
                    "lineage_anchor_uninformative_or_representation_mismatch",
                }
                else "regional_candidate_not_linked_to_named_anchor" if best_locus_candidate
                else "no_candidate"
            ),
            "outgroup_label": panel.active_label,
            "outgroup_fallback_used": int(panel.fallback_used),
            "n_lineage_diagnostic_sites": best_locus_candidate["n_diagnostic_sites"] if best_locus_candidate else max(
                (len(payload["diagnostic_indexes"]) for payload in lineage_payload.values()), default=0
            ),
            "n_candidate_haplotypes": best_locus_candidate["candidate_types"] if best_locus_candidate else 0,
            "n_candidate_copies": best_locus_candidate["candidate_copies"] if best_locus_candidate else 0,
            "candidate_start": best_locus_candidate["start"] if best_locus_candidate else "",
            "candidate_end": best_locus_candidate["end"] if best_locus_candidate else "",
            "candidate_bp": best_locus_candidate["bp"] if best_locus_candidate else "",
            "candidate_tree_status": "pending",
            "candidate_tree_bootstrap": "",
            "candidate_tree_pass": 0,
            "candidate_tree_strong": 0,
            "multi_evidence_strong": 0,
            "candidate_tree_candidate_purity": "",
            "candidate_tree_candidate_sensitivity": "",
            "candidate_tree_candidates_in_clade": "",
            "candidate_tree_archaic_in_clade": "",
            "introgression_call": prep_call,
            "introgression_confidence": prep_confidence,
            "independent_segment_caller_required": 1,
            "evidence_reason": prep_reason,
            "n_unknown_sample_copies": total_group_counts["UNKNOWN"],
            "unknown_samples_first10": ",".join(sorted(unknown_samples)[:10]),
        }
        evidence_locus_rows.append(locus_row)
        legacy_locus_additions[(locus_id,)] = {key: locus_row.get(key, "") for key in EVIDENCE_LOCUS_FIELDS}

        # Build one candidate-focused tree per lineage.  It contains only
        # pre-specified candidate haplotypes, hard/common modern controls, all
        # archaic references, and one ancestral outgroup.
        for candidate_payload in locus_candidates:
            lineage = str(candidate_payload["lineage"])
            candidates = sorted(
                candidate_payload["scores"],
                key=lambda row: (
                    -as_int(row["n"]),
                    -as_int(row["n_contiguous_diagnostic_match"]),
                    -(as_float(row["diagnostic_match_prop"], 0.0) or 0.0),
                    str(row["hap_id"]),
                ),
            )[: args.max_candidate_haplotype_types]
            candidate_ids = {str(row["hap_id"]) for row in candidates}
            consensus = str(lineage_payload[lineage]["consensus"])
            diag_indexes = list(lineage_payload[lineage]["diagnostic_indexes"])
            other_scores = [row for row in scored_by_lineage[lineage] if str(row["hap_id"]) not in candidate_ids]
            copy_only_reasons = {"candidate_copy_count_low", "nonoutgroup_copy_count_low"}
            supporting_context = [
                row for row in other_scores
                if set(filter(None, str(row.get("diagnostic_candidate_reason", "")).split(";")))
                and set(filter(None, str(row.get("diagnostic_candidate_reason", "")).split(";"))).issubset(copy_only_reasons)
            ]
            supporting_context = sorted(
                supporting_context,
                key=lambda row: (
                    -as_int(row["n_diagnostic_match"]),
                    -(as_float(row["diagnostic_match_prop"], 0.0) or 0.0),
                    -as_int(row["n_contiguous_diagnostic_match"]),
                    str(row["hap_id"]),
                ),
            )[: args.max_supporting_contexts]
            supporting_ids = {str(row["hap_id"]) for row in supporting_context}
            control_scores = [row for row in other_scores if str(row["hap_id"]) not in supporting_ids]
            recurrent_controls = [row for row in control_scores if as_int(row["n"]) >= args.min_control_copies]
            singleton_controls = [row for row in control_scores if as_int(row["n"]) < args.min_control_copies]
            # Match the classic locus analysis by preferring haplotypes observed
            # at least twice; singletons are used only if recurrent controls are
            # insufficient. Hard controls have high diagnostic affinity but
            # fail at least one non-circular criterion.
            control_pool = recurrent_controls + singleton_controls
            hard = sorted(
                control_pool,
                key=lambda row: (
                    -as_int(row["n_diagnostic_match"]),
                    -(as_float(row["diagnostic_match_prop"], 0.0) or 0.0),
                    -as_int(row["n"]),
                    str(row["hap_id"]),
                ),
            )
            hard_n = args.evidence_controls // 2
            controls = hard[:hard_n]
            selected_control_ids = {str(row["hap_id"]) for row in controls}
            common = sorted(
                [row for row in control_pool if str(row["hap_id"]) not in selected_control_ids],
                key=lambda row: (
                    int(as_int(row["n"]) < args.min_control_copies),
                    -as_int(row["n"]),
                    str(row["hap_id"]),
                ),
            )
            controls += common[: max(0, args.evidence_controls - len(controls))]

            source_haps = {row.get("hap_id", ""): row for row in hap_rows}
            tree_records: list[dict[str, object]] = []
            for row in candidates:
                hap_id = str(row["hap_id"])
                source = source_haps[hap_id]
                tree_records.append({
                    "label": hap_id,
                    "seq": source.get("seq", "").upper(),
                    "role": "candidate",
                    "lineage": lineage,
                    "n": source.get("n", ""),
                })
            for row in supporting_context:
                hap_id = str(row["hap_id"])
                source = source_haps[hap_id]
                tree_records.append({
                    "label": hap_id,
                    "seq": source.get("seq", "").upper(),
                    "role": "candidate_context",
                    "lineage": lineage,
                    "n": source.get("n", ""),
                })
            for row in controls:
                hap_id = str(row["hap_id"])
                source = source_haps[hap_id]
                tree_records.append({
                    "label": hap_id,
                    "seq": source.get("seq", "").upper(),
                    "role": "control",
                    "lineage": "modern",
                    "n": source.get("n", ""),
                })
            for row in archaic_rows:
                tree_records.append({
                    "label": row.get("archaic", "Archaic"),
                    "seq": row.get("seq", "").upper(),
                    "role": "archaic",
                    "lineage": lineage_name(row.get("lineage") or row.get("archaic", "")),
                    "n": 1,
                })
            if ancestral and sum(base in BASES for base in ancestral) >= args.min_tree_sites:
                tree_records.append({
                    "label": "Ancestral",
                    "seq": ancestral,
                    "role": "outgroup",
                    "lineage": "Ancestral",
                    "n": 1,
                })

            # Build the tree over the observed candidate tract (optionally
            # with a small flank), not the entire user-supplied locus.  This
            # prevents a short ABO/X candidate from being diluted by unrelated
            # variation elsewhere in a broad interval.  Candidate definition
            # used only diagnostic alleles; topology uses every variable site
            # available inside this tract, so the tree is not a re-count of the
            # diagnostic score.
            expected_archaic_sequences = [
                str(row["seq"]) for row in tree_records
                if row["role"] == "archaic" and row["lineage"] == lineage
            ]
            raw_tree_start = as_int(candidate_payload.get("start"), 0)
            raw_tree_end = as_int(candidate_payload.get("end"), 0)
            tree_start = max(0, raw_tree_start - args.tree_flank_bp)
            tree_end = raw_tree_end + args.tree_flank_bp
            keep_indexes: list[int] = []
            for i in range(n_sites):
                pos = as_int(sites[i].get("pos"), -1)
                if raw_tree_end > raw_tree_start and not (tree_start < pos <= tree_end):
                    continue
                if not any(seq[i] in BASES for seq in expected_archaic_sequences):
                    continue
                called = {str(row["seq"])[i] for row in tree_records if str(row["seq"])[i] in BASES}
                if len(called) >= 2:
                    keep_indexes.append(i)
            if len(keep_indexes) < args.min_tree_sites:
                continue

            slug = safe_slug(lineage)
            phy = locus_dir / f"haplotypes.evidence.{slug}.phy"
            meta_path = Path(f"{phy}.meta.tsv")
            manifest_path = locus_dir / f"haplotypes.evidence.{slug}.manifest.tsv"
            sites_path = locus_dir / f"haplotypes.evidence.{slug}.sites.tsv"
            desired_phy_files.add(phy)
            labels: list[dict[str, object]] = []
            phy_lines = [f"{len(tree_records)} {len(keep_indexes)}"]
            for index, row in enumerate(tree_records, 1):
                short = f"E{index:09d}"
                sequence = "".join(str(row["seq"])[i] for i in keep_indexes)
                phy_lines.append(f"{short:<10} {sequence}")
                labels.append({
                    "phy_label": short,
                    "label": row["label"],
                    "role": row["role"],
                    "lineage": row["lineage"],
                    "expected_lineage": lineage,
                    "n": row["n"],
                    "tree_region_start": tree_start,
                    "tree_region_end": tree_end,
                    "tree_region_bp": max(0, tree_end - tree_start),
                    "n_tree_sites": len(keep_indexes),
                })
            phy_text = "\n".join(phy_lines) + "\n"
            changed = write_text_if_changed(phy, phy_text)
            if changed:
                remove_phyml_derivatives(phy)
            write_tsv(meta_path, ["phy_label", "label"], labels)
            write_tsv(
                manifest_path,
                [
                    "phy_label", "label", "role", "lineage", "expected_lineage", "n",
                    "tree_region_start", "tree_region_end", "tree_region_bp", "n_tree_sites",
                ],
                labels,
            )
            write_tsv(
                sites_path,
                ["tree_index", "source_index", "chr", "pos", "id", "ref", "alt", "diagnostic_for_expected_lineage"],
                (
                    {
                        "tree_index": tree_i + 1,
                        "source_index": source_i + 1,
                        "chr": sites[source_i].get("chr", meta.get("chr", "")),
                        "pos": sites[source_i].get("pos", ""),
                        "id": sites[source_i].get("id", ""),
                        "ref": sites[source_i].get("ref", ""),
                        "alt": sites[source_i].get("alt", ""),
                        "diagnostic_for_expected_lineage": int(source_i in diag_indexes),
                    }
                    for tree_i, source_i in enumerate(keep_indexes)
                ),
            )

    # Remove stale evidence PHYLIP inputs for lineages no longer selected.
    for phy in (args.out / "loci").glob("**/haplotypes.evidence.*.phy"):
        if phy not in desired_phy_files:
            remove_phyml_derivatives(phy)
            phy.unlink(missing_ok=True)
            Path(f"{phy}.meta.tsv").unlink(missing_ok=True)
            phy.with_name(phy.name.replace(".phy", ".manifest.tsv")).unlink(missing_ok=True)
            phy.with_name(phy.name.replace(".phy", ".sites.tsv")).unlink(missing_ok=True)

    site_fields = [
        "locus_id", "chr", "pos", "id", "ref", "alt", "lineage",
        "lineage_consensus", "ancestral", "lineage_reference_count",
        "lineage_agreeing_reference_count", "lineage_agreeing_references",
        "outgroup_label", "outgroup_called_copies", "outgroup_archaic_allele_copies",
        "outgroup_archaic_allele_frequency", "nonoutgroup_called_copies",
        "nonoutgroup_archaic_allele_copies", "nonoutgroup_archaic_allele_frequency",
        "archaic_differs_from_ancestral", "diagnostic_site_pass", "diagnostic_site_reason",
    ]
    hap_fields = [
        "locus_id", "hap_id", "n", "copies", "raw_best_archaic", "raw_best_lineage",
        "raw_n_compared", "raw_n_match", "raw_prop_match", *EVIDENCE_HAP_FIELDS,
        "seq",
    ]
    locus_fields = [
        "locus_id", "genome_build", "chr", "core_start", "core_end", "name",
        *EVIDENCE_LOCUS_FIELDS,
        "outgroup_fallback_used", "n_unknown_sample_copies", "unknown_samples_first10",
    ]
    final = out / "final"
    write_tsv(final / "evidence_sites.tsv", site_fields, evidence_site_rows)
    write_tsv(final / "evidence_haplotypes.tsv", hap_fields, evidence_hap_rows)
    write_tsv(final / "evidence_loci.tsv", locus_fields, evidence_locus_rows)

    # Explicit stage boundaries make it impossible to confuse the exhaustive
    # regional inventory with later filtered subsets.  evidence_haplotypes.tsv
    # is the complete fate table; the following files are convenient pass-only
    # views and never replace or delete Stage-1 rows.
    diagnostic_sites = [row for row in evidence_site_rows if as_int(row.get("diagnostic_site_pass")) == 1]
    diagnostic_candidates = [row for row in evidence_hap_rows if as_int(row.get("diagnostic_candidate_pass")) == 1]
    anchor_zoom_candidates = [
        row for row in diagnostic_candidates
        if as_int(row.get("named_anchor_lineage_matching_copies")) > 0
    ]
    write_tsv(final / "diagnostic_sites.filtered.tsv", site_fields, diagnostic_sites)
    write_tsv(final / "diagnostic_candidates.filtered.tsv", hap_fields, diagnostic_candidates)
    write_tsv(final / "anchor_zoom_candidates.tsv", hap_fields, anchor_zoom_candidates)

    params = {
        "schema_version": 3,
        "method": "region_first_then_diagnostic_filter_anchor_optional_tree_confirm",
        "regional_stage1_files": {
            "all_haplotypes": "region_haplotypes.unfiltered.tsv.gz",
            "all_haplotype_reference_pairs": "region_matches.unfiltered.tsv.gz",
            "local_segments": "region_match_segments.unfiltered.tsv.gz",
            "families": "region_candidate_families.unfiltered.tsv.gz",
        },
        "filtered_stage_files": {
            "all_site_fates": "evidence_sites.tsv",
            "passing_diagnostic_sites": "diagnostic_sites.filtered.tsv",
            "all_haplotype_lineage_fates": "evidence_haplotypes.tsv",
            "passing_diagnostic_candidates": "diagnostic_candidates.filtered.tsv",
            "anchor_linked_zoom": "anchor_zoom_candidates.tsv",
        },
        "anchor_is_gate": False,
        "outgroup_requested": args.outgroup,
        "outgroup_fallback": args.outgroup_fallback,
        "outgroup_active": panel.active_label,
        "outgroup_fallback_used": panel.fallback_used,
        "sample_file": str(args.sample_file.resolve()) if args.sample_file and args.sample_file.exists() else "",
        "thresholds": {
            "max_outgroup_allele_freq": args.max_outgroup_allele_freq,
            "max_outgroup_haplotype_freq": args.max_outgroup_haplotype_freq,
            "min_outgroup_copies": args.min_outgroup_copies,
            "min_diagnostic_sites": args.min_diagnostic_sites,
            "min_diagnostic_match_prop": args.min_diagnostic_match_prop,
            "min_candidate_copies": args.min_candidate_copies,
            "max_diagnostic_gap_bp": args.max_diagnostic_gap_bp,
            "min_neanderthal_refs": args.min_neanderthal_refs,
            "min_other_lineage_refs": args.min_other_lineage_refs,
            "evidence_controls": args.evidence_controls,
            "max_supporting_contexts": args.max_supporting_contexts,
            "min_control_copies": args.min_control_copies,
            "max_candidate_haplotype_types": args.max_candidate_haplotype_types,
            "min_tree_sites": args.min_tree_sites,
            "tree_flank_bp": args.tree_flank_bp,
            "min_tree_bootstrap": args.min_tree_bootstrap,
            "strong_tree_bootstrap": args.strong_tree_bootstrap,
            "min_candidate_purity": args.min_candidate_purity,
            "min_candidate_sensitivity": args.min_candidate_sensitivity,
            "min_tree_candidate_types": args.min_tree_candidate_types,
            "strong_min_diagnostic_sites": args.strong_min_diagnostic_sites,
            "strong_min_candidate_copies": args.strong_min_candidate_copies,
            "strong_min_archaic_tips": args.strong_min_archaic_tips,
        },
    }
    write_text_if_changed(final / "evidence_parameters.json", json.dumps(params, indent=2, sort_keys=True) + "\n")

    merge_fields(out / "final" / "haplotypes.tsv", ("locus_id", "hap_id"), legacy_hap_additions, EVIDENCE_HAP_FIELDS)
    merge_fields(out / "final" / "loci.tsv", ("locus_id",), legacy_locus_additions, EVIDENCE_LOCUS_FIELDS)
    print(
        f"PHYML EVIDENCE PREPARE: loci={len(evidence_locus_rows)} "
        f"diagnostic_sites={sum(as_int(row['diagnostic_site_pass']) for row in evidence_site_rows)} "
        f"candidate_haplotypes={sum(as_int(row['diagnostic_candidate_pass']) for row in evidence_hap_rows)} "
        f"outgroup={panel.active_label} output={final}",
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
    return [float(x) for x in re.findall(r"\)(-?[0-9]+(?:\.[0-9]+)?)(?=[:),;])", newick) if float(x) >= 0]


def evidence_tree_result(
    newick: str,
    manifest: list[dict[str, str]],
    min_bootstrap: float,
    strong_bootstrap: float,
    min_purity: float,
    min_sensitivity: float,
    min_modern_tips: int,
) -> dict[str, object]:
    expected_lineages = {row.get("expected_lineage", "") for row in manifest if row.get("expected_lineage", "")}
    expected_lineage = next(iter(expected_lineages)) if len(expected_lineages) == 1 else ""
    by_label = {row.get("label", ""): row for row in manifest}
    candidates_all = {
        label for label, row in by_label.items()
        if row.get("role") == "candidate" and row.get("lineage") == expected_lineage
    }
    blank = {
        "expected_lineage": expected_lineage,
        "tree_status": "complete" if newick else "not_run",
        "tree_has_ancestral_outgroup": 0,
        "candidate_clade_n_tips": 0,
        "n_candidate_tips_total": len(candidates_all),
        "candidate_clade_node": "",
        "candidate_clade_side": "",
        "candidate_clade_bootstrap": "",
        "candidate_clade_pass": 0,
        "candidate_clade_strong": 0,
        "candidate_purity": "",
        "candidate_sensitivity": "",
        "candidate_f1": "",
        "n_candidate_tips_in_clade": 0,
        "n_control_tips_in_clade": 0,
        "n_candidate_context_tips_in_clade": 0,
        "n_expected_archaic_tips_in_clade": 0,
        "n_other_archaic_tips_in_clade": 0,
        "candidate_tips_in_clade": "",
        "control_tips_in_clade": "",
        "candidate_context_tips_in_clade": "",
        "expected_archaic_tips_in_clade": "",
        "tree_call_reason": "tree_missing" if not newick else "no_candidate_archaic_edge",
    }
    if not newick or not expected_lineage:
        return blank
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
    blank["tree_has_ancestral_outgroup"] = int("Ancestral" in all_tips)
    candidates: list[tuple] = []
    for node in internal:
        if node is root:
            continue
        down = descendants[node.node_id]
        for side_name, side in (("descendant", down), ("complement", all_tips - down)):
            if "Ancestral" in side:
                continue
            cand = {tip for tip in side if tip in candidates_all}
            controls = {
                tip for tip in side
                if by_label.get(tip, {}).get("role") in {"control", "candidate"} and tip not in cand
            }
            context = {
                tip for tip in side
                if by_label.get(tip, {}).get("role") == "candidate_context"
            }
            expected_archaic = {
                tip for tip in side
                if by_label.get(tip, {}).get("role") == "archaic"
                and by_label.get(tip, {}).get("lineage") == expected_lineage
            }
            other_archaic = {
                tip for tip in side
                if by_label.get(tip, {}).get("role") == "archaic"
                and by_label.get(tip, {}).get("lineage") != expected_lineage
            }
            if len(cand) < min_modern_tips or not expected_archaic or other_archaic:
                continue
            purity = len(cand) / (len(cand) + len(controls)) if cand or controls else 0.0
            sensitivity = len(cand) / len(candidates_all) if candidates_all else 0.0
            f1 = 2 * purity * sensitivity / (purity + sensitivity) if purity + sensitivity else 0.0
            support = node.support if node.support is not None and node.support >= 0 else None
            passed = bool(
                support is not None
                and support >= min_bootstrap
                and purity >= min_purity
                and sensitivity >= min_sensitivity
            )
            strong = bool(
                passed
                and support is not None
                and support >= strong_bootstrap
                and purity >= max(min_purity, 0.90)
                and sensitivity >= max(min_sensitivity, 0.75)
            )
            candidates.append((
                int(passed), int(strong), support if support is not None else -1.0,
                f1, sensitivity, purity, len(cand), len(expected_archaic),
                -len(controls), len(context), side_name, node, cand, controls, context, expected_archaic, other_archaic,
            ))
    if not candidates:
        return blank
    selected = max(candidates)
    passed_i, strong_i, support, f1, sensitivity, purity, _, _, _, _, side_name, node, cand, controls, context, expected_archaic, other_archaic = selected
    return {
        **blank,
        "candidate_clade_node": node.node_id,
        "candidate_clade_side": side_name,
        "candidate_clade_n_tips": len(cand) + len(controls) + len(context) + len(expected_archaic) + len(other_archaic),
        "candidate_clade_bootstrap": "" if support < 0 else support,
        "candidate_clade_pass": passed_i,
        "candidate_clade_strong": strong_i,
        "candidate_purity": format_float(purity),
        "candidate_sensitivity": format_float(sensitivity),
        "candidate_f1": format_float(f1),
        "n_candidate_tips_in_clade": len(cand),
        "n_control_tips_in_clade": len(controls),
        "n_candidate_context_tips_in_clade": len(context),
        "n_expected_archaic_tips_in_clade": len(expected_archaic),
        "n_other_archaic_tips_in_clade": len(other_archaic),
        "candidate_tips_in_clade": ",".join(sorted(cand)),
        "control_tips_in_clade": ",".join(sorted(controls)),
        "candidate_context_tips_in_clade": ",".join(sorted(context)),
        "expected_archaic_tips_in_clade": ",".join(sorted(expected_archaic)),
        "tree_call_reason": (
            "strong_candidate_lineage_clade" if strong_i
            else "supported_candidate_lineage_clade" if passed_i
            else "candidate_edge_below_bootstrap_purity_or_sensitivity"
        ),
    }


def summarize(args: argparse.Namespace) -> None:
    out: Path = args.out
    evidence_loci = read_tsv(out / "final" / "evidence_loci.tsv")
    if not evidence_loci:
        fail(f"missing evidence preparation output: {out / 'final' / 'evidence_loci.tsv'}")
    tree_rows: list[dict[str, object]] = []
    trees_by_locus: dict[str, list[dict[str, object]]] = defaultdict(list)

    for locus_dir in locus_dirs(out):
        locus_meta = final_rows_by_locus(out, "loci.tsv")
        locus_id = infer_locus_id(locus_dir, locus_meta)
        for manifest_path in sorted(locus_dir.glob("haplotypes.evidence.*.manifest.tsv")):
            stem = manifest_path.name.replace(".manifest.tsv", "")
            phy = locus_dir / f"{stem}.phy"
            tree = Path(f"{phy}_phyml_tree.txt")
            stats = Path(f"{phy}_phyml_stats.txt")
            meta = Path(f"{phy}.meta.tsv")
            manifest = read_tsv(manifest_path)
            expected = next((row.get("expected_lineage", "") for row in manifest if row.get("expected_lineage", "")), "")
            if tree.is_file() and tree.stat().st_size and meta.is_file():
                label_map = {row["phy_label"]: row["label"] for row in read_tsv(meta) if row.get("phy_label") and row.get("label")}
                newick = relabel_newick(tree.read_text(errors="replace"), label_map)
                result = evidence_tree_result(
                    newick,
                    manifest,
                    args.min_tree_bootstrap,
                    args.strong_tree_bootstrap,
                    args.min_candidate_purity,
                    args.min_candidate_sensitivity,
                    args.min_tree_candidate_types,
                )
                supports = bootstrap_values(newick)
                status = "complete"
            else:
                newick = ""
                result = evidence_tree_result(
                    "",
                    manifest,
                    args.min_tree_bootstrap,
                    args.strong_tree_bootstrap,
                    args.min_candidate_purity,
                    args.min_candidate_sensitivity,
                    args.min_tree_candidate_types,
                )
                supports = []
                status = args.default_tree_status
                result["tree_status"] = status
                result["tree_call_reason"] = "tree_not_requested" if status == "not_requested" else "tree_missing"
            row = {
                "locus_id": locus_id,
                "expected_lineage": expected,
                "tree_status": status,
                "tree_file": str(tree.resolve()) if tree.is_file() else "",
                "stats_file": str(stats.resolve()) if stats.is_file() else "",
                "plot_file": str(tree.with_suffix(".png").resolve()) if tree.with_suffix(".png").is_file() else "",
                "phy_file": str(phy.resolve()) if phy.is_file() else "",
                "tree_region_start": manifest[0].get("tree_region_start", "") if manifest else "",
                "tree_region_end": manifest[0].get("tree_region_end", "") if manifest else "",
                "tree_region_bp": manifest[0].get("tree_region_bp", "") if manifest else "",
                "n_tree_sites": manifest[0].get("n_tree_sites", "") if manifest else "",
                "n_tree_tips": len(manifest),
                "n_bootstrap_nodes": len(supports),
                "bootstrap_min": min(supports) if supports else "",
                "bootstrap_median": median(supports) if supports else "",
                "bootstrap_max": max(supports) if supports else "",
                **result,
                "tree_newick": newick,
            }
            tree_rows.append(row)
            trees_by_locus[locus_id].append(row)

    updated_loci: list[dict[str, object]] = []
    legacy_locus_additions: dict[tuple[str, ...], dict[str, object]] = {}
    legacy_tree_additions: dict[tuple[str, ...], dict[str, object]] = {}
    for source in evidence_loci:
        row: dict[str, object] = dict(source)
        locus_id = source.get("locus_id", "")
        sequence_candidates = as_int(source.get("n_candidate_haplotypes"))
        sequence_copies = as_int(source.get("n_candidate_copies"))
        lineage = source.get("evidence_lineage", "")
        matching_trees = [tree for tree in trees_by_locus.get(locus_id, []) if tree.get("expected_lineage") == lineage]
        tree = max(
            matching_trees,
            key=lambda item: (
                as_int(item.get("candidate_clade_pass")),
                as_int(item.get("candidate_clade_strong")),
                as_float(item.get("candidate_clade_bootstrap"), -1.0) or -1.0,
                as_float(item.get("candidate_f1"), -1.0) or -1.0,
            ),
            default=None,
        )

        prep_call = source.get("introgression_call", "")
        tree_bootstrap = as_float(tree.get("candidate_clade_bootstrap"), -1.0) if tree else -1.0
        tree_purity = as_float(tree.get("candidate_purity"), 0.0) if tree else 0.0
        tree_sensitivity = as_float(tree.get("candidate_sensitivity"), 0.0) if tree else 0.0
        anchor_linked = as_int(source.get("n_anchor_linked_candidate_copies")) > 0
        anchor_evaluable = str(source.get("named_anchor_support_status", "")) not in {
            "named_anchor_not_provided",
            "named_anchor_not_found",
            "named_anchor_not_phase_assigned",
            "lineage_anchor_uninformative_or_representation_mismatch",
        }
        # ``strong_candidate`` is multi-component by construction.  A high
        # bootstrap alone cannot upgrade a three-site or singleton signal.
        # Conversely, a known locus such as chr3 is not mechanically downgraded
        # merely because one candidate-focused edge is 70--89 when diagnostic,
        # outgroup, carrier, lineage and topology evidence all converge.
        multi_component_strong = bool(
            tree
            and as_int(tree.get("candidate_clade_pass")) == 1
            and (tree_bootstrap or -1.0) >= args.min_tree_bootstrap
            and as_int(source.get("n_lineage_diagnostic_sites")) >= args.strong_min_diagnostic_sites
            and sequence_copies >= args.strong_min_candidate_copies
            and as_int(source.get("evidence_lineage_reference_count")) >= 2
            and as_int(tree.get("tree_has_ancestral_outgroup")) == 1
            and as_int(tree.get("n_expected_archaic_tips_in_clade")) >= args.strong_min_archaic_tips
            and (tree_purity or 0.0) >= max(args.min_candidate_purity, 0.90)
            and (tree_sensitivity or 0.0) >= max(args.min_candidate_sensitivity, 0.75)
        )
        if prep_call == "not_evaluable":
            call = "not_evaluable"
            confidence = "NA"
            reason = source.get("evidence_reason", "not_evaluable")
        elif sequence_candidates == 0:
            call = "not_supported"
            confidence = "none"
            reason = "no_haplotype_passed_diagnostic_sequence_criteria"
        elif tree is None:
            call = "candidate_sequence_only"
            confidence = "low"
            reason = "candidate_tree_not_constructed"
        elif tree.get("tree_status") != "complete":
            call = "candidate_sequence_only"
            confidence = "low"
            reason = str(tree.get("tree_call_reason", "candidate_tree_not_complete"))
        elif multi_component_strong:
            call = "strong_candidate"
            confidence = "high"
            reason = "multiple_archaic_derived_outgroup_filtered_diagnostics_and_lineage_concordant_tree"
        elif as_int(tree.get("candidate_clade_pass")) == 1:
            call = "supported_candidate"
            confidence = "moderate"
            reason = "diagnostic_outgroup_filtered_tract_and_lineage_concordant_tree"
        elif sequence_candidates < args.min_tree_candidate_types:
            call = "candidate_sequence_only"
            confidence = "low"
            reason = "insufficient_candidate_haplotype_types_for_tree_confirmation"
        else:
            call = "candidate_sequence_only"
            confidence = "low"
            reason = str(tree.get("tree_call_reason", "candidate_tree_not_supported"))

        if call in {"strong_candidate", "supported_candidate"} and anchor_linked:
            anchor_call = call
            anchor_confidence = confidence
            introgression_scope = "named_anchor_linked_region"
        elif call in {"strong_candidate", "supported_candidate", "candidate_sequence_only"} and not anchor_evaluable:
            anchor_call = "not_evaluable"
            anchor_confidence = "NA"
            introgression_scope = "regional_candidate_anchor_uninformative"
        elif call in {"strong_candidate", "supported_candidate", "candidate_sequence_only"}:
            anchor_call = "not_supported"
            anchor_confidence = "none"
            introgression_scope = "regional_candidate_not_linked_to_named_anchor"
        else:
            anchor_call = "not_evaluable" if not anchor_evaluable else "not_supported"
            anchor_confidence = "NA" if not anchor_evaluable else "none"
            introgression_scope = "not_evaluable" if call == "not_evaluable" else "no_candidate"

        row.update({
            "candidate_tree_status": tree.get("tree_status", "not_constructed") if tree else "not_constructed",
            "candidate_tree_bootstrap": tree.get("candidate_clade_bootstrap", "") if tree else "",
            "candidate_tree_pass": tree.get("candidate_clade_pass", 0) if tree else 0,
            "candidate_tree_strong": tree.get("candidate_clade_strong", 0) if tree else 0,
            "multi_evidence_strong": int(multi_component_strong),
            "anchor_linked_introgression_call": anchor_call,
            "anchor_linked_introgression_confidence": anchor_confidence,
            "introgression_scope": introgression_scope,
            "candidate_tree_candidate_purity": tree.get("candidate_purity", "") if tree else "",
            "candidate_tree_candidate_sensitivity": tree.get("candidate_sensitivity", "") if tree else "",
            "candidate_tree_candidates_in_clade": tree.get("candidate_tips_in_clade", "") if tree else "",
            "candidate_tree_archaic_in_clade": tree.get("expected_archaic_tips_in_clade", "") if tree else "",
            "introgression_call": call,
            "introgression_confidence": confidence,
            "independent_segment_caller_required": 1,
            "evidence_reason": reason,
        })
        updated_loci.append(row)
        legacy_locus_additions[(locus_id,)] = {key: row.get(key, "") for key in EVIDENCE_LOCUS_FIELDS}
        authoritative_tree_update: dict[str, object] = {
            "tree_scope": "candidate_focused_authoritative" if tree and tree.get("tree_status") == "complete" else "exploratory_tree_with_authoritative_no_call",
            "tree_interpretation": "candidate_first_outgroup_filtered_tree_confirmation",
            "authoritative_introgression_call": call,
            "authoritative_evidence_lineage": lineage,
            "authoritative_candidate_tree_status": row.get("candidate_tree_status", ""),
            "authoritative_candidate_tree_bootstrap": row.get("candidate_tree_bootstrap", ""),
            "authoritative_candidate_tree_pass": row.get("candidate_tree_pass", 0),
            "candidate_lineage": lineage,
            "tree_has_ancestral_outgroup": tree.get("tree_has_ancestral_outgroup", 0) if tree else 0,
            "candidate_clade_pass": row.get("candidate_tree_pass", 0),
            "candidate_clade_bootstrap": row.get("candidate_tree_bootstrap", ""),
            "candidate_clade_node": tree.get("candidate_clade_node", "") if tree else "",
            "candidate_clade_side": tree.get("candidate_clade_side", "") if tree else "",
            "candidate_clade_n_tips": tree.get("candidate_clade_n_tips", 0) if tree else 0,
            "candidate_clade_modern_tips": tree.get("n_candidate_tips_in_clade", 0) if tree else 0,
            "candidate_clade_archaic_tips": tree.get("n_expected_archaic_tips_in_clade", 0) if tree else 0,
            "candidate_clade_specificity": tree.get("candidate_purity", "") if tree else "",
            "candidate_clade_tips": (
                ",".join(filter(None, [
                    str(tree.get("candidate_tips_in_clade", "")),
                    str(tree.get("expected_archaic_tips_in_clade", "")),
                ])) if tree else ""
            ),
            "candidate_clade_rule": "diagnostic_tract_then_same_lineage_tree;bootstrap_purity_sensitivity",
        }
        if tree and tree.get("tree_status") == "complete":
            authoritative_tree_update.update({
                "tree_status": "complete",
                "tree_file": tree.get("tree_file", ""),
                "stats_file": tree.get("stats_file", ""),
                "plot_file": tree.get("plot_file", ""),
                "n_bootstrap_nodes": tree.get("n_bootstrap_nodes", ""),
                "bootstrap_min": tree.get("bootstrap_min", ""),
                "bootstrap_median": tree.get("bootstrap_median", ""),
                "bootstrap_max": tree.get("bootstrap_max", ""),
                "tree_newick": tree.get("tree_newick", ""),
            })
        legacy_tree_additions[(locus_id,)] = authoritative_tree_update

    tree_fields = [
        "locus_id", "expected_lineage", "tree_status", "tree_file", "stats_file", "plot_file", "phy_file",
        "tree_region_start", "tree_region_end", "tree_region_bp", "n_tree_sites", "n_tree_tips",
        "n_bootstrap_nodes", "bootstrap_min", "bootstrap_median", "bootstrap_max",
        "tree_has_ancestral_outgroup", "n_candidate_tips_total", "candidate_clade_node", "candidate_clade_side",
        "candidate_clade_n_tips",
        "candidate_clade_bootstrap", "candidate_clade_pass", "candidate_clade_strong",
        "candidate_purity", "candidate_sensitivity", "candidate_f1",
        "n_candidate_tips_in_clade", "n_control_tips_in_clade",
        "n_candidate_context_tips_in_clade",
        "n_expected_archaic_tips_in_clade", "n_other_archaic_tips_in_clade",
        "candidate_tips_in_clade", "control_tips_in_clade",
        "candidate_context_tips_in_clade",
        "expected_archaic_tips_in_clade", "tree_call_reason", "tree_newick",
    ]
    locus_fields = list(evidence_loci[0].keys())
    for field_name in EVIDENCE_LOCUS_FIELDS:
        if field_name not in locus_fields:
            locus_fields.append(field_name)
    write_tsv(out / "final" / "evidence_trees.tsv", tree_fields, tree_rows)
    write_tsv(out / "final" / "evidence_loci.tsv", locus_fields, updated_loci)
    merge_fields(out / "final" / "loci.tsv", ("locus_id",), legacy_locus_additions, EVIDENCE_LOCUS_FIELDS)
    primary_trees = out / "final" / "trees.tsv"
    exploratory_trees = out / "final" / "exploratory_trees.tsv"
    existing_primary_rows = read_tsv(primary_trees)
    # Preserve the legacy all-haplotype tree once.  On a repeated summarize
    # call, final/trees.tsv is already the candidate-focused authoritative
    # table; copying it over exploratory_trees.tsv would destroy the original
    # exploratory topology and make cache/reclassification behavior unstable.
    already_authoritative = bool(existing_primary_rows) and all(
        str(row.get("tree_interpretation", "")) == "candidate_first_outgroup_filtered_tree_confirmation"
        for row in existing_primary_rows
    )
    if existing_primary_rows and not already_authoritative:
        write_tsv(exploratory_trees, list(existing_primary_rows[0].keys()), existing_primary_rows)
    merge_fields(
        primary_trees,
        ("locus_id",),
        legacy_tree_additions,
        (
            "tree_scope", "tree_interpretation", "authoritative_introgression_call",
            "authoritative_evidence_lineage", "authoritative_candidate_tree_status",
            "authoritative_candidate_tree_bootstrap", "authoritative_candidate_tree_pass",
            "candidate_lineage", "tree_has_ancestral_outgroup", "candidate_clade_pass",
            "candidate_clade_bootstrap", "candidate_clade_node", "candidate_clade_side",
            "candidate_clade_n_tips", "candidate_clade_modern_tips",
            "candidate_clade_archaic_tips", "candidate_clade_specificity",
            "candidate_clade_tips", "candidate_clade_rule",
            "tree_status", "tree_file", "stats_file", "plot_file",
            "n_bootstrap_nodes", "bootstrap_min", "bootstrap_median",
            "bootstrap_max", "tree_newick",
        ),
    )
    counts = Counter(str(row.get("introgression_call", "")) for row in updated_loci)
    print(
        "PHYML EVIDENCE SUMMARY: "
        + " ".join(f"{key}={value}" for key, value in sorted(counts.items()))
        + f" output={out / 'final' / 'evidence_loci.tsv'}",
        flush=True,
    )


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--sample-file", type=Path)
    parser.add_argument("--outgroup", default="YRI")
    parser.add_argument("--outgroup-fallback", default="AFR")
    parser.add_argument("--max-outgroup-allele-freq", type=float, default=0.05)
    parser.add_argument("--max-outgroup-haplotype-freq", type=float, default=0.05)
    parser.add_argument("--min-outgroup-copies", type=int, default=20)
    parser.add_argument("--min-diagnostic-sites", type=int, default=3)
    parser.add_argument("--min-diagnostic-match-prop", type=float, default=0.80)
    parser.add_argument("--min-candidate-copies", type=int, default=2)
    parser.add_argument("--max-diagnostic-gap-bp", type=int, default=50_000)
    parser.add_argument("--min-neanderthal-refs", type=int, default=2)
    parser.add_argument("--min-other-lineage-refs", type=int, default=1)
    parser.add_argument("--evidence-controls", type=int, default=24)
    parser.add_argument("--max-supporting-contexts", type=int, default=24)
    parser.add_argument("--min-control-copies", type=int, default=2)
    parser.add_argument("--max-candidate-haplotype-types", type=int, default=1000)
    parser.add_argument("--min-tree-sites", type=int, default=10)
    parser.add_argument("--tree-flank-bp", type=int, default=5_000)
    parser.add_argument("--min-tree-bootstrap", type=float, default=70.0)
    parser.add_argument("--strong-tree-bootstrap", type=float, default=90.0)
    parser.add_argument("--strong-min-diagnostic-sites", type=int, default=10)
    parser.add_argument("--strong-min-candidate-copies", type=int, default=5)
    parser.add_argument("--strong-min-archaic-tips", type=int, default=2)
    parser.add_argument("--min-candidate-purity", type=float, default=0.80)
    parser.add_argument("--min-candidate-sensitivity", type=float, default=0.50)
    parser.add_argument("--min-tree-candidate-types", type=int, default=1)
    parser.add_argument("--default-tree-status", choices=("not_requested", "not_run"), default="not_run")


def validate_args(args: argparse.Namespace) -> None:
    for name in ("max_outgroup_allele_freq", "max_outgroup_haplotype_freq", "min_diagnostic_match_prop", "min_candidate_purity", "min_candidate_sensitivity"):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            fail(f"--{name.replace('_', '-')} must be between 0 and 1")
    for name in (
        "min_outgroup_copies", "min_diagnostic_sites", "min_candidate_copies",
        "max_diagnostic_gap_bp", "min_neanderthal_refs", "min_other_lineage_refs",
        "max_candidate_haplotype_types", "min_tree_sites", "min_tree_candidate_types",
        "strong_min_diagnostic_sites", "strong_min_candidate_copies", "strong_min_archaic_tips", "min_control_copies",
    ):
        if getattr(args, name) < 1:
            fail(f"--{name.replace('_', '-')} must be positive")
    if args.evidence_controls < 0:
        fail("--evidence-controls cannot be negative")
    if args.max_supporting_contexts < 0:
        fail("--max-supporting-contexts cannot be negative")
    if args.tree_flank_bp < 0:
        fail("--tree-flank-bp cannot be negative")
    if not 0 <= args.min_tree_bootstrap <= 100 or not 0 <= args.strong_tree_bootstrap <= 100:
        fail("tree bootstrap thresholds must be between 0 and 100")
    if args.strong_tree_bootstrap < args.min_tree_bootstrap:
        fail("--strong-tree-bootstrap must be >= --min-tree-bootstrap")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("prepare", "summarize", "all"):
        child = sub.add_parser(action)
        add_common_arguments(child)
    args = parser.parse_args()
    validate_args(args)
    if args.action in {"prepare", "all"}:
        prepare(args)
    if args.action in {"summarize", "all"}:
        summarize(args)


if __name__ == "__main__":
    main()
