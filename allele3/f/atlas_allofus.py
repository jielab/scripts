#!/usr/bin/env python3
"""Join completed All of Us AF-only queries to exact local Atlas substitutions.

Inputs must use the same target-file content and GRCh38 coordinates. Preserve
population frequencies and availability separately from AVI scores; do not
reconstruct REF frequencies, counts, denominators or functional distances.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import json
import math
from pathlib import Path
import re

from analysis_integrity import InputPath, identity, transactional, validate_manifest
from allofus_shards import VERSION as ALLOFUS_VERSION, AF_COLUMNS, TARGET_FIELDS
from atlas_scores import (VERSION as ATLAS_VERSION, BASES, MISSING, VARIANT_FIELDS,
                          chromosome, chromosome_key, genome_build, table, writer)

VERSION = "atlas-allofus-join-1"
FIELDS = ["site_id", "chrom", "pos", "build", "ref", "alt", "allele_role", "target_alts", "observed_from_targets",
          "avi_raw", "avi_phred", "atlas_status", "atlas_resource_id", "atlas_source_version",
          "allofus_availability", "allofus_match_status", "allofus_requested_alt", "allofus_observed_in_target",
          "allofus_additional_reported_alt", "allofus_source_ref", "allofus_evidence", *AF_COLUMNS, "gvs_max_af"]
LIMITS = [
    "This is an exact chrom/pos/REF/ALT join of completed derivatives from the same target-file content. It is not an allele discovery scan or a source-cohort independence test.",
    "Eight All of Us frequency columns and gvs_max_af retain their source meanings. gvs_max_af is a cross-population maximum indicator, not a ninth population. No denominators, counts, REF frequencies or joint allele probabilities are reconstructed.",
    "Unreported theoretical ALTs have NA frequencies, not zero. A reported AF=0 is preserved and does not demonstrate AC=0. Missing metrics, missing rows and uncallable alleles are not extinction.",
    "observed_from_targets, allofus_requested_alt and allofus_observed_in_target describe target-list membership. allofus_additional_reported_alt describes an additional summary row, not independently established genotype-level tri-allelism.",
    "The structural REF row has NA AF and NA AVI scores. REF is not assumed ancestral or most frequent. AVI raw/PHRED are reference-relative predictions; no ALT-to-ALT distance, pleiotropy count or selection statistic is inferred.",
    "Frequency provenance and prediction provenance remain separate. This join does not combine WGS counts with gnomAD denominators, ancient observations or imputed genotypes.",
]


def is_missing(value):
    return str(value).strip().upper() in MISSING


def number(value, label, frequency=False, nonnegative=False):
    if is_missing(value):
        return "NA"
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a numeric value or NA") from exc
    if not math.isfinite(parsed) or (frequency and not 0 <= parsed <= 1) or (nonnegative and parsed < 0):
        raise ValueError(f"Invalid {label}: {value}")
    return value


def bit(value, label):
    if value not in {"0", "1"}:
        raise ValueError(f"{label} must be an explicit 0/1")
    return int(value)


def row_key(row):
    if genome_build(row["build"]) != "GRCh38":
        raise ValueError("All join rows must declare GRCh38")
    if not re.fullmatch(r"[0-9]+", row["pos"]) or int(row["pos"]) < 1:
        raise ValueError("All join rows need positive 1-based positions")
    if len(row["ref"]) != 1 or row["ref"] not in BASES or not row["site_id"] or is_missing(row["site_id"]):
        raise ValueError("Join rows require a canonical REF and a site_id")
    return chromosome(row["chrom"]), int(row["pos"])


def require_inventory(manifest, filenames):
    if not set(filenames) <= set(manifest["output_inventory"]):
        raise ValueError("Required join tables are absent from the validated source output inventory")


def validate_sources(allofus_dir, atlas_dir):
    aou = validate_manifest(allofus_dir, require_evidence=True)
    atlas = validate_manifest(atlas_dir)
    if (aou.get("schema") != ALLOFUS_VERSION or aou.get("status") != "complete_descriptive_unverified"
            or aou.get("evidence") != "frequency_only" or "targets" not in aou or "target_count" not in aou):
        raise ValueError("--allofus-dir must be a completed allofus_shards.query AF-only result")
    if atlas.get("schema") != ATLAS_VERSION or atlas.get("status") not in {"complete", "complete_with_exclusions", "not_available"}:
        raise ValueError("--atlas-dir must be a completed atlas_scores result, including an explicit no-score audit if unavailable")
    if genome_build(aou.get("build")) != "GRCh38" or genome_build(atlas.get("genome_build")) != "GRCh38" or atlas.get("coordinate_system") != "1-based":
        raise ValueError("Both completed sources must use GRCh38, with 1-based Atlas coordinates")
    require_inventory(aou, ["target_alleles.tsv.gz", "target_sites.tsv"])
    require_inventory(atlas, ["atlas_variants.tsv.gz", "target_availability.tsv"])
    target_a, target_b = aou["targets"], atlas.get("source_identity", {}).get("targets", {})
    if not re.fullmatch(r"[a-fA-F0-9]{64}", str(target_a.get("sha256", ""))) or not re.fullmatch(r"[a-fA-F0-9]{64}", str(target_b.get("sha256", ""))):
        raise ValueError("Both source manifests must retain full target-file SHA256 identities")
    if any(target_a.get(k) != target_b.get(k) for k in ("sha256", "size")) or target_a.get("size", 0) <= 0:
        raise ValueError("Target-file content identities differ between All of Us and Atlas; do not join unrelated target requests")
    return aou, atlas


def read_atlas(directory, manifest):
    sites, site_ids = {}, {}
    avi_resources = [r for r in manifest.get("resources", []) if r.get("kind") == "avi" and r.get("status") == "ready"]
    if len(avi_resources) != 1:
        raise ValueError("Atlas manifest must identify its one active AVI resource")
    source = avi_resources[0]
    for _, row in table(Path(directory) / "atlas_variants.tsv.gz", VARIANT_FIELDS):
        key = row_key(row)
        sid, ref, allele = row["site_id"], row["ref"], row["allele"]
        if row["resource_id"] != source["resource_id"] or row["source_version"] != source["source_version"]:
            raise ValueError("Atlas row provenance conflicts with its resource manifest")
        alts = row["target_alts"].split(",")
        if not alts or len(alts) != len(set(alts)) or any(len(v) != 1 or v not in BASES or v == ref for v in alts):
            raise ValueError("Invalid target ALT membership in Atlas rows")
        if sid in site_ids and site_ids[sid] != key:
            raise ValueError("Atlas site_id refers to multiple coordinates")
        site_ids[sid] = key
        site = sites.setdefault(key, dict(site_id=sid, ref=ref, alts=set(alts), rows={}))
        if site["site_id"] != sid or site["ref"] != ref or site["alts"] != set(alts):
            raise ValueError("Atlas rows have inconsistent site identity, REF or target ALTs")
        if len(allele) != 1 or allele not in BASES or allele in site["rows"]:
            raise ValueError("Atlas requires exactly one row for each of the four bases")
        reference = allele == ref
        if row["allele_role"] != ("REF" if reference else "ALT"):
            raise ValueError("Atlas REF/ALT role does not match its exact nucleotide")
        if reference:
            if not all(is_missing(row[k]) for k in ("avi_raw", "avi_phred", "observed_from_targets")) or row["status"] != "reference_label_no_substitution":
                raise ValueError("Atlas REF must retain NA metrics and target-observation flag")
        elif bit(row["observed_from_targets"], "observed_from_targets") != int(allele in site["alts"]):
            raise ValueError("Atlas target membership flags conflict with target_alts")
        row["avi_raw"] = number(row["avi_raw"], "avi_raw")
        row["avi_phred"] = number(row["avi_phred"], "avi_phred", nonnegative=True)
        if not reference:
            numeric = sum(row[k] != "NA" for k in ("avi_raw", "avi_phred"))
            expected = {"available"} if numeric == 2 else {"available_partial"} if numeric else {"score_missing", "not_scored", "not_indexed_contig"}
            if row["status"] not in expected:
                raise ValueError("Atlas availability status conflicts with stored scores")
        site["rows"][allele] = row
    if not sites or len(sites) != manifest.get("counts", {}).get("target_sites") or any(set(s["rows"]) != set(BASES) for s in sites.values()):
        raise ValueError("Atlas is missing target sites or required theoretical-base rows")
    return sites


def read_allofus(directory, manifest, sites):
    rows, coordinates, requested_seen = {}, set(), {}
    for _, row in table(Path(directory) / "target_alleles.tsv.gz", TARGET_FIELDS):
        key = row_key(row)
        if key not in sites:
            raise ValueError("All of Us and Atlas target coordinate sets differ; region-filtered requests need matching target subsets")
        site = sites[key]
        if row["site_id"] != site["site_id"] or row["ref"] != site["ref"]:
            raise ValueError("All of Us / Atlas site_id or REF mismatch at the exact coordinate")
        alt = row["alt"]
        if len(alt) != 1 or alt not in BASES or alt == row["ref"] or (key, alt) in rows:
            raise ValueError("All of Us query requires unique canonical nonreference substitutions")
        if (not is_missing(row["source_ref"]) and row["source_ref"] != row["ref"]) or row["match_status"] == "reference_mismatch":
            raise ValueError("All of Us source REF conflicts with the harmonized target REF")
        if row["match_status"] not in {"matched_frequency_only", "matched_metrics_missing", "outside_prepared_windows", "no_reported_snv_rows_in_window", "not_in_summary"}:
            raise ValueError("Unrecognized All of Us query availability status")
        requested = bit(row["requested_alt"], "requested_alt")
        observed = bit(row["observed_in_target"], "observed_in_target")
        additional = bit(row["additional_reported_alt"], "additional_reported_alt")
        if requested != observed or requested != int(alt in site["alts"]) or requested and additional:
            raise ValueError("All of Us / Atlas target-list membership flags disagree")
        matched = row["match_status"] in {"matched_frequency_only", "matched_metrics_missing"}
        if additional != int(not requested and matched):
            raise ValueError("All of Us additional ALT flag conflicts with its availability status")
        for col in [*AF_COLUMNS, "gvs_max_af"]:
            row[col] = number(row[col], col, frequency=True)
        known_af = any(row[col] != "NA" for col in AF_COLUMNS)
        if (matched and is_missing(row["source_ref"])) or (not matched and any(row[col] != "NA" for col in [*AF_COLUMNS, "gvs_max_af"])):
            raise ValueError("All of Us AF values or source REF conflict with row availability")
        if (row["match_status"] == "matched_frequency_only" and not known_af) or (row["match_status"] == "matched_metrics_missing" and known_af):
            raise ValueError("All of Us metric availability is inconsistent")
        expected_evidence = "frequency_only" if known_af else "metrics_missing" if matched else "unavailable"
        if row["evidence"] != expected_evidence:
            raise ValueError("All of Us row evidence differs from its AF-only contract")
        rows[(key, alt)] = row
        coordinates.add(key)
        if requested:
            requested_seen.setdefault(key, set()).add(alt)
    if coordinates != set(sites) or len(coordinates) != manifest["target_count"]:
        raise ValueError("All of Us and Atlas completed target-site sets differ")
    if any(requested_seen.get(key, set()) != site["alts"] for key, site in sites.items()):
        raise ValueError("All of Us query omitted a target ALT instead of preserving an unavailable row")
    return rows


def protected(raw):
    return [InputPath(Path(raw[arg]) / filename) for arg, filenames in (
        ("allofus_dir", ["manifest.json", "target_alleles.tsv.gz", "target_sites.tsv"]),
        ("atlas_dir", ["manifest.json", "atlas_variants.tsv.gz", "target_availability.tsv"])) for filename in filenames]


@transactional("manifest.json", input_dirs=("allofus_dir", "atlas_dir"), input_prefixes=(), resolve_inputs=protected)
def run(allofus_dir, atlas_dir, out):
    aou, atlas = validate_sources(allofus_dir, atlas_dir)
    sites = read_atlas(atlas_dir, atlas)
    aou_rows = read_allofus(allofus_dir, aou, sites)
    counts = Counter(target_sites=len(sites))
    with gzip.open(Path(out) / "allofus_atlas.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        out_rows = writer(handle, FIELDS)
        for key in sorted(sites, key=lambda v: (chromosome_key(v[0]), v[1])):
            site = sites[key]
            for allele in [site["ref"], *[v for v in BASES if v != site["ref"]]]:
                ar = site["rows"][allele]
                aour = aou_rows.get((key, allele))
                reference = allele == site["ref"]
                record = dict(site_id=site["site_id"], chrom=key[0], pos=key[1], build="GRCh38", ref=site["ref"], alt=allele,
                              allele_role=ar["allele_role"], target_alts=ar["target_alts"], observed_from_targets=ar["observed_from_targets"],
                              avi_raw=ar["avi_raw"], avi_phred=ar["avi_phred"], atlas_status=ar["status"],
                              atlas_resource_id=ar["resource_id"], atlas_source_version=ar["source_version"])
                for source in ("match_status", "requested_alt", "observed_in_target", "additional_reported_alt", "source_ref", "evidence"):
                    record["allofus_" + source] = aour[source] if aour else "NA"
                for col in [*AF_COLUMNS, "gvs_max_af"]:
                    record[col] = aour[col] if aour else "NA"
                availability = "reference_AF_unavailable" if reference else "no_exact_alt_row_in_query" if aour is None else aour["match_status"]
                record["allofus_availability"] = availability
                out_rows.writerow(record)
                counts["rows"] += 1
                counts["allofus_status_" + availability] += 1
                if not reference:
                    af_available = any(record[col] != "NA" for col in AF_COLUMNS)
                    avi_available = ar["status"] in {"available", "available_partial"}
                    counts["substitutions_with_AF"] += af_available
                    counts["substitutions_with_AVI"] += avi_available
                    counts["substitutions_with_AF_and_AVI"] += af_available and avi_available
                    counts["substitutions_with_complete_metrics"] += ar["status"] == "available" and all(record[col] != "NA" for col in [*AF_COLUMNS, "gvs_max_af"])
    # Revalidate complete inventories, not only the subset of files read in this join.
    validate_sources(allofus_dir, atlas_dir)
    print(f"[build-check] atlas_allofus=GRCh38 same_target_SHA256=checked sites={len(sites)} exact_site_REF_ALT=checked", flush=True)
    print(f"[dataset] allofus=frequency_only populations=8 max_AF=separate AVI=reference_relative jointly_available_substitutions={counts['substitutions_with_AF_and_AVI']}", flush=True)
    return dict(status="complete" if counts["substitutions_with_complete_metrics"] == 3 * len(sites) else "complete_with_exclusions",
                schema=VERSION, build="GRCh38", coordinate_system="1-based", counts=dict(counts), limits=LIMITS,
                target_identity=aou["targets"], target_identity_check="identical full SHA256 and size; exact site/REF/ALT and membership cross-check",
                source_manifests=dict(allofus=identity(Path(allofus_dir) / "manifest.json", full=True),
                                      atlas=identity(Path(atlas_dir) / "manifest.json", full=True)),
                allofus_provenance=dict(source_identity=aou.get("source_identity"), plan_id=aou.get("plan_id"), evidence="frequency_only"),
                atlas_provenance=dict(source_identity=atlas.get("source_identity"), resources=atlas.get("resources")),
                denominator_policy="not available in AF-only source; never reconstructed, pooled or transferred")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--allofus-dir", required=True, help="Completed allofus_shards query directory")
    p.add_argument("--atlas-dir", required=True, help="Completed atlas_scores directory for exactly the same target-file content")
    p.add_argument("--out", required=True)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        result = run(args.allofus_dir, args.atlas_dir, args.out)
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        raise SystemExit(f"atlas-allofus: {exc}")
    print(json.dumps(dict(status=result["status"], counts=result["counts"]), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
