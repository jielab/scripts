#!/usr/bin/env python3
"""HERV allele catalogs and seeded validation panels from completed site caches.

This is a conditional candidate design: both cases and controls already have
at least two observed ALT nucleotides with AC>=5. It does not estimate HERV
enrichment or test differences in the allele frequencies used for matching.
"""
from __future__ import annotations
from analysis_integrity import transactional

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
import os
from pathlib import Path
import random
import re

import scan
from allele_common import IntervalMask, af_bin, an_bin, canonical_chrom, chrom_key, open_text


FAMILIES = ("ERV1", "ERVK")
AUTOSOMES = {str(i) for i in range(1, 23)}
FIELDS = ["chrom", "pos", "ref", "alts", "site_id", "family", "group", "match_id",
          "overlap_family", "case_region", "candidate_families", "context", "canonical_context",
          "partition", "comparison_partition", "body_bin", "spatial_overlap", "family_partitions",
          "primary_alt", "primary_af", "primary_an", "secondary_alt", "secondary_af",
          "primary_af_bin", "an_bin", "secondary_af_bin", "n_alt_ac5", "n_alt_ac10",
          "ac10_supported", "observed_alts_all", "alt_ac_all", "alt_an_all",
          "ref_presence_status", "strict_filter_complete", "missing_strict_filters",
          "match_stratum", "sampling_probability", "sampling_weight"]
MASK_OPTIONS = {"lcr": "architecture_lcr_bed", "segdup": "architecture_segdup_bed",
                "low_mappability": "architecture_low_mappability_bed"}


def metadata(path):
    p = Path(path).expanduser().resolve()
    s = p.stat()
    return {"path": str(p), "size": s.st_size, "mtime_ns": s.st_mtime_ns}


def write_json(path, value):
    p = Path(path); tmp = p.with_name(p.name + ".partial")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(tmp, p)


def write_tsv(path, fields, rows):
    p = Path(path); tmp = p.with_name(p.name + ".partial")
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, p)


class Reservoir:
    """Uniform bounded sample; sample size is independent of stream length."""
    def __init__(self, capacity, seed):
        self.capacity = capacity; self.rng = random.Random(seed)
        self.n = 0; self.rows = []

    def add(self, row):
        self.n += 1
        if self.capacity is None or len(self.rows) < self.capacity:
            self.rows.append(row)
        else:
            j = self.rng.randrange(self.n)
            if j < self.capacity:
                self.rows[j] = row


def stratified_quotas(counts, maximum):
    """Proportional largest-remainder quotas; one per chromosome when feasible."""
    counts = {k: n for k, n in counts.items() if n}
    budget = min(maximum, sum(counts.values()))
    if not budget:
        return {k: 0 for k in counts}
    keys = sorted(counts, key=chrom_key)
    base = int(budget >= len(keys))
    quotas = {k: base for k in keys}
    remaining = budget - sum(quotas.values())
    available = {k: counts[k] - quotas[k] for k in keys}
    total = sum(available.values())
    raw = {k: remaining * available[k] / total if total else 0 for k in keys}
    for k in keys:
        quotas[k] += math.floor(raw[k])
    left = budget - sum(quotas.values())
    order = sorted(keys, key=lambda k: (-(raw[k] - math.floor(raw[k])), chrom_key(k)))
    for k in order[:left]:
        quotas[k] += 1
    return quotas


def load_sources(modern_dir):
    root = Path(modern_dir).expanduser().resolve()
    run_path = root / "run.json"
    run = scan.load_run_metadata(root)
    if run.get("status") not in {"scanned", "complete"}:
        raise ValueError("Use a completed/scanned modern run.json")
    options = run.get("options", run.get("identity", {}).get("options", {}))
    if options.get("build", "GRCh38") not in {"38", "GRCh38", "hg38"}:
        raise ValueError("Candidates currently require GRCh38 modern shards and annotation")
    shards = []
    for raw in run.get("shards", []):
        d = Path(raw).expanduser()
        if not d.is_absolute(): d = root / d
        d = d.resolve()
        receipt = json.loads((d / "COMPLETE.json").read_text())
        if not run.get("signature") or receipt.get("signature") != run["signature"]:
            raise ValueError(f"Shard completion signature differs from run.json: {d}")
        site = d / "multiallelic_sites.tsv.gz"
        saved = receipt.get("outputs", {}).get(site.name)
        if not site.is_file() or (saved and saved.get("size") != site.stat().st_size):
            raise ValueError(f"Missing or size-changed site cache: {site}")
        info = json.loads((d / "scan.json").read_text())
        if not info.get("complete"):
            raise ValueError(f"Incomplete scan.json: {d}")
        scope = info.get("scope", [])
        if not scope:
            raise ValueError(f"Missing chromosome scope in {d}/scan.json")
        shards.append({"directory": d, "site": site, "flags": set(info.get("fields_read", [])) & {"lcr", "segdup"},
                       "scope": scope, "metadata": metadata(site), "scan_metadata": metadata(d / "scan.json")})
    if not shards or len({s["directory"] for s in shards}) != len(shards):
        raise ValueError("No shards or duplicate shard directories in run.json")
    shards.sort(key=lambda s: (min(chrom_key(c) for c in s["scope"]), str(s["directory"])))
    return run, options, shards


def split_annotation(source, options, cache):
    """Reuse verified-by-metadata chromosome cache read-only, otherwise split once."""
    source_meta = metadata(source)
    def valid_files(directory, manifest):
        for i in range(1, 23):
            name = f"chr{i}.tsv.gz"; p = directory / name
            saved = manifest.get("files", {}).get(name, {})
            if not p.is_file() or not saved: return False
            actual = metadata(p)
            if any(saved.get(k) != actual[k] for k in ("size", "mtime_ns")): return False
        return True
    previous = options.get("architecture_annotation_cache")
    if previous:
        marker = Path(previous) / "manifest.json"
        if marker.is_file():
            manifest = json.loads(marker.read_text())
            recorded = manifest.get("source", {})
            same = all(recorded.get(k) == v for k, v in source_meta.items())
            if same and valid_files(Path(previous), manifest):
                return Path(previous), "original_read_only"
    cache.mkdir(parents=True, exist_ok=True)
    marker = cache / "manifest.json"
    if marker.is_file():
        old = json.loads(marker.read_text())
        if old.get("source") == source_meta and valid_files(cache, old):
            return cache, "derived_chromosome_cache"
    handles = {}
    try:
        for i in range(1, 23): handles[str(i)] = gzip.open(cache / f"chr{i}.tsv.gz.partial", "wt", compresslevel=2)
        with open_text(source) as f:
            for ln, line in enumerate(f, 1):
                if not line.strip() or line.startswith("#"): continue
                # The final joined-RMSK field may be empty; preserve trailing tabs.
                a = line.rstrip("\r\n").split("\t") if "\t" in line else line.split()
                if len(a) != 14:
                    raise ValueError(f"{source}:{ln}: expected 14 joined RepeatMasker columns, got {len(a)}")
                chrom = canonical_chrom(a[1])
                if chrom in handles: handles[chrom].write(line if line.endswith("\n") else line + "\n")
    finally:
        for h in handles.values(): h.close()
    for i in range(1, 23): os.replace(cache / f"chr{i}.tsv.gz.partial", cache / f"chr{i}.tsv.gz")
    write_json(marker, {"source": source_meta, "method": "autosomal split; metadata consistency, no full-file checksums",
                        "files": {f"chr{i}.tsv.gz": metadata(cache / f"chr{i}.tsv.gz") for i in range(1, 23)}})
    return cache, "derived_chromosome_cache"


def load_masks(options):
    # Compact merged interval arrays are reused between both passes.
    return {name: IntervalMask(options[key], [str(i) for i in range(1, 23)])
            for name, key in MASK_OPTIONS.items() if options.get(key)}


def parse_site(row):
    """Preserve observed ALT identities; never infer an observed REF allele."""
    chrom = canonical_chrom(row["chrom"]); pos = int(row["pos"])
    ref = row["ref"].upper()
    alts = row["alts_by_af"].split(",")
    ac = [int(x) for x in row["alt_ac"].split(",")]
    an = [int(x) for x in row["alt_an"].split(",")]
    if len({len(alts), len(ac), len(an)}) != 1 or len(set(alts)) != len(alts):
        raise ValueError(f"Misaligned/duplicated ALT arrays at {chrom}:{pos}")
    if pos < 1 or ref not in "ACGT" or len(ref) != 1 or any(a not in {"A", "C", "G", "T"} or a == ref for a in alts):
        raise ValueError(f"Invalid observed SNV identity at {chrom}:{pos}")
    if any(n <= 0 or c < 1 or c > n for c, n in zip(ac, an)):
        raise ValueError(f"Invalid ALT AC/AN at {chrom}:{pos}")
    n5 = sum(c >= 5 for c in ac); n10 = sum(c >= 10 for c in ac)
    if "n_alt_ac5" in row and int(row["n_alt_ac5"]) != n5:
        raise ValueError(f"Cached AC5 count disagrees with observed ALTs at {chrom}:{pos}")
    if "n_alt_ac10" in row and int(row["n_alt_ac10"]) != n10:
        raise ValueError(f"Cached AC10 count disagrees with observed ALTs at {chrom}:{pos}")
    if n5 < 2: return None
    kept = sorted([(a, c, n, c / n) for a, c, n in zip(alts, ac, an) if c >= 5], key=lambda x: (-x[3], x[0]))
    context = row.get("context", "NA")
    canonical = row.get("canonical_context") or scan.canonical_context(context, ref)[0]
    first, second = kept[:2]
    return {"chrom": chrom, "pos": pos, "ref": ref, "alts": ",".join(x[0] for x in kept),
            "site_id": f"{chrom}:{pos}:{ref}", "context": context, "canonical_context": canonical,
            "primary_alt": first[0], "primary_af": first[3], "primary_an": first[2],
            "secondary_alt": second[0], "secondary_af": second[3],
            "primary_af_bin": af_bin(first[3]), "an_bin": an_bin(first[2]), "secondary_af_bin": af_bin(second[3]),
            "n_alt_ac5": n5, "n_alt_ac10": n10, "ac10_supported": int(n10 >= 2),
            "observed_alts_all": ",".join(alts), "alt_ac_all": ",".join(map(str, ac)),
            "alt_an_all": ",".join(map(str, an)), "ref_presence_status": "not_inferable_from_sites_summary"}


class ChromAnnotation:
    def __init__(self, cache, chrom, *, end_bp, body_bins):
        from herv_modern import AlignedBodies
        self.elements = scan.parse_joined_rmsk(cache / f"chr{chrom}.tsv.gz", chromosome=chrom)
        self.body = scan.ERVOverlapSweep(self.elements)
        self.spatial = AlignedBodies(self.elements, FAMILIES, end_bp=end_bp, bins=body_bins)
        self.last_spatial = {}
        self.boundaries = {fam: scan.BoundarySweep([e for e in self.elements if e.rep_class == "LTR" and e.family == fam])
                           for fam in FAMILIES}

    def classify(self, chrom, pos):
        overlap = self.body.locate(chrom, pos)
        self.last_spatial = self.spatial.assign(overlap, pos)
        body_families = {e.family for e in overlap if e.rep_class == "LTR" and e.family in FAMILIES}
        near200, near2000 = set(), set()
        for fam, sweep in self.boundaries.items():
            for edge in ("5p", "3p"):
                near = sweep.nearest(chrom, pos, edge, max_distance=2000)
                if near:
                    near2000.add(fam)
                    if near["distance"] < 200: near200.add(fam)
        return body_families, near200, bool(near2000), any(e.rep_class == "LTR" for e in overlap), any(
            e.rep_class in {"Simple_repeat", "Low_complexity"} for e in overlap)


def eligible_rows(shards, annotation_cache, masks, counts, *, end_bp, body_bins):
    last = None; current = None; annotation = None
    for shard in shards:
        print(f"[herv-candidates] read {shard['site']}", flush=True)
        with gzip.open(shard["site"], "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for raw in reader:
                chrom = canonical_chrom(raw["chrom"]); key = (chrom_key(chrom), int(raw["pos"]))
                if last is not None and key <= last:
                    raise ValueError("Site caches must be globally chromosome/position sorted and have no duplicate positions")
                last = key; counts["rows_read"] += 1
                if chrom not in AUTOSOMES: continue
                if "n_alt_ac5" in raw and int(raw["n_alt_ac5"]) < 2: continue
                row = parse_site(raw)
                if row is None: continue
                counts["ac5_eligible"] += 1
                flags = set(re.split("[,;|]", raw.get("region_flags", ".")))
                missing = sorted(set(MASK_OPTIONS) - set(masks) - shard["flags"])
                # Match the original ArchitectureAccumulator: an explicit BED
                # defines that mask; INFO is a fallback only when no BED exists.
                excluded = any(mask.contains(chrom, row["pos"]) for mask in masks.values())
                excluded |= any(name in flags for name in ("lcr", "segdup") if name not in masks)
                if excluded:
                    counts["excluded_difficult_regions"] += 1; continue
                if current != chrom:
                    current = chrom; annotation = ChromAnnotation(annotation_cache, chrom, end_bp=end_bp, body_bins=body_bins)
                body, near, near_control, any_ltr, simple = annotation.classify(chrom, row["pos"])
                if simple:
                    counts["excluded_simple_low_complexity"] += 1; continue
                row.update(strict_filter_complete=int(not missing), missing_strict_filters=",".join(missing) or ".")
                row["_spatial_parts"] = annotation.last_spatial
                row["family_partitions"] = ";".join(f"{family}:{value[0]}" for family,value in sorted(annotation.last_spatial.items())) or "."
                row["match_stratum"] = "|".join(str(row[k]) for k in ("chrom", "canonical_context", "primary_af_bin", "an_bin", "secondary_af_bin"))
                yield row, body, near, near_control, any_ltr


def match_key(row):
    # A per-shard missing annotation must not silently weaken just one match arm.
    return row["match_stratum"], row["missing_strict_filters"]


@transactional("candidate_manifest.json")
def run_candidates(args):
    if args.max_sites < 0 or args.controls_per_case < 0:
        raise ValueError("--max-sites must be nonnegative (0 means all); --controls-per-case must be nonnegative")
    end_bp, body_bins = args.end_bp, args.body_bins
    if end_bp < 1 or body_bins < 2:
        raise ValueError("--end-bp must be positive and --body-bins at least two")
    spatial_definition = {"schema": "aligned-record-partition-1", "end_bp": end_bp, "body_bins": body_bins}
    catalog_all = bool(getattr(args, "all_alleles_catalog", False))
    root = Path(args.modern_dir).expanduser().resolve(); out = Path(args.out).expanduser().resolve()
    run, options, shards = load_sources(root)
    if out == root or any(out == s["directory"] or out.is_relative_to(s["directory"]) for s in shards):
        raise ValueError("Candidate output must not overwrite the original modern run/shards")
    identity = {"schema": 3, "spatial_definition": spatial_definition,
                "code": {name: metadata(Path(__file__).with_name(name))
                         for name in ("herv_candidates.py", "herv_modern.py", "scan.py", "scan_store.py", "0.common.py")},
                "modern_run": metadata(root / "run.json"), "source_signature": run.get("signature"),
                "sites": [s["metadata"] for s in shards], "scan_metadata": [s["scan_metadata"] for s in shards],
                "erv_joined": metadata(args.erv_joined),
                "masks": {k: metadata(options[v]) for k, v in MASK_OPTIONS.items() if options.get(v)},
                "max_sites_per_family": args.max_sites, "all_alleles_catalog": catalog_all,
                "controls_per_case": args.controls_per_case, "seed": args.seed,
                "method": "body_or_abs_boundary_lt200; controls_noLTR_and_abs_boundary_ge2000; AC5-supported_ALT_AF_matching_v1"}
    receipt = out / "candidate_manifest.json"
    if receipt.is_file():
        old = json.loads(receipt.read_text())
        if old.get("status") == "complete" and old.get("identity") == identity:
            if all((out / p).is_file() and metadata(out / p) == fp for p, fp in old.get("outputs", {}).items()) and old.get("outputs"):
                print(f"[herv-candidates] reuse {out}", flush=True)
                return old
    out.mkdir(parents=True, exist_ok=True)
    annotation_cache, cache_mode = split_annotation(args.erv_joined, options, out / ".cache" / "joined")
    masks = load_masks(options)
    pools = {}; pass1 = Counter()
    for row, body, near, _, _ in eligible_rows(shards, annotation_cache, masks, pass1, end_bp=end_bp, body_bins=body_bins):
        families = body | near
        if not families: continue
        # A multiply annotated locus belongs to one seeded family for selection;
        # original memberships are retained, preventing duplicated target loci.
        family = random.Random(f"{args.seed}:family:{row['site_id']}").choice(sorted(families))
        part, body_bin, spatial_overlap = row["_spatial_parts"].get(family, ("boundary_only", -1, "outside_aligned_body"))
        row.update(family=family, group="case", overlap_family=",".join(sorted(body)) or "boundary_only",
                   candidate_families=",".join(sorted(families)), case_region=("body_and_boundary" if family in body & near else "body" if family in body else "boundary"),
                   partition=part, comparison_partition=part, body_bin=body_bin, spatial_overlap=spatial_overlap)
        key = family, row["chrom"]
        if key not in pools: pools[key] = Reservoir(args.max_sites or None, f"{args.seed}:case:{family}:{row['chrom']}")
        pools[key].add(row)
    cases = []; sampling = []
    for family in FAMILIES:
        counts = {c: p.n for (f, c), p in pools.items() if f == family}
        quotas = stratified_quotas(counts, args.max_sites or sum(counts.values()))
        for chrom in sorted(counts, key=chrom_key):
            pool = pools[family, chrom]; quota = quotas[chrom]
            random.Random(f"{args.seed}:quota:{family}:{chrom}").shuffle(pool.rows)
            probability = quota / pool.n
            for row in pool.rows[:quota]:
                row.update(sampling_probability=probability, sampling_weight=1 / probability)
                cases.append(row)
            sampling.append({"family": family, "chrom": chrom, "eligible_cases": pool.n, "selected_cases": quota,
                             "sampling_probability": probability, "sampling_weight": 1 / probability if probability else "NA"})
    cases.sort(key=lambda r: (r["family"], chrom_key(r["chrom"]), r["pos"]))
    demands = defaultdict(list)
    for i, row in enumerate(cases, 1):
        row["match_id"] = f"HERV{i:06d}"
        demands[match_key(row)].append(row)
    control_pools = {k: Reservoir(len(v) * args.controls_per_case, f"{args.seed}:control:{k}")
                     for k, v in demands.items() if args.controls_per_case}
    pass2 = Counter()
    if control_pools:
        for row, _, _, near_control, any_ltr in eligible_rows(shards, annotation_cache, masks, pass2, end_bp=end_bp, body_bins=body_bins):
            if any_ltr or near_control: continue
            key = match_key(row)
            if key in control_pools: control_pools[key].add(row)
    targets = list(cases); audit = []; used_controls = set()
    for key, group in sorted(demands.items()):
        pool = control_pools.get(key)
        controls = list(pool.rows) if pool else []
        random.Random(f"{args.seed}:allocate:{key}").shuffle(controls)
        matched = Counter()
        for i, control in enumerate(controls):
            case = group[i % len(group)]
            if control["site_id"] in used_controls:
                raise ValueError("Duplicate control locus encountered")
            used_controls.add(control["site_id"]); matched[case["match_id"]] += 1
            control.update(family=case["family"], group="control", match_id=case["match_id"], overlap_family="nonHERV",
                           candidate_families=".", case_region="nonLTR_distant_control", sampling_probability="NA", sampling_weight="NA",
                           partition="nonLTR_control", comparison_partition=case["partition"], body_bin=-1, spatial_overlap="nonLTR_control")
            targets.append(control)
        for case in group:
            n = matched[case["match_id"]]
            audit.append({"match_id": case["match_id"], "family": case["family"], "site_id": case["site_id"],
                          "match_stratum": case["match_stratum"], "missing_strict_filters": case["missing_strict_filters"],
                          "requested_controls": args.controls_per_case, "matched_controls": n,
                          "eligible_controls_in_stratum": pool.n if pool else 0, "cases_in_stratum": len(group),
                          "status": "complete" if n == args.controls_per_case else "insufficient_exact_controls"})
    if (any(metadata(s["site"]) != s["metadata"] or metadata(s["directory"] / "scan.json") != s["scan_metadata"] for s in shards)
            or metadata(root / "run.json") != identity["modern_run"]
            or metadata(args.erv_joined) != identity["erv_joined"]
            or any(metadata(options[MASK_OPTIONS[k]]) != saved for k, saved in identity["masks"].items())):
        raise ValueError("Inputs changed during candidate selection; rerun with stable inputs")
    targets.sort(key=lambda r: (r["match_id"], r["group"] != "case", chrom_key(r["chrom"]), r["pos"]))
    write_tsv(out / "targets.tsv", FIELDS, targets)
    write_tsv(out / "match_audit.tsv", ["match_id", "family", "site_id", "match_stratum", "missing_strict_filters",
              "requested_controls", "matched_controls", "eligible_controls_in_stratum", "cases_in_stratum", "status"], audit)
    write_tsv(out / "sampling.tsv", ["family", "chrom", "eligible_cases", "selected_cases", "sampling_probability", "sampling_weight"], sampling)
    output_names = ["targets.tsv", "match_audit.tsv", "sampling.tsv"]
    if catalog_all:
        write_allele_catalogs(out, targets)
        output_names += ["site_catalog.tsv.gz", "allele_catalog.tsv.gz"]
    result = {"status": "complete", "identity": identity, "spatial_definition": spatial_definition,
              "case_boundary_bp": 200, "control_boundary_exclusion_bp": 2000,
              "annotation_cache_mode": cache_mode,
              "cases_by_family": dict(Counter(r["family"] for r in cases)), "n_cases": len(cases), "n_controls": len(used_controls),
              "n_unmatched_cases": sum(r["matched_controls"] == 0 for r in audit), "n_incomplete_matches": sum(r["status"] != "complete" for r in audit),
              "pass1": dict(pass1), "pass2": dict(pass2),
              "all_selected_strict_complete": bool(targets) and all(r["strict_filter_complete"] for r in targets),
              "coordinate_build": "GRCh38", "site_coordinates": "1-based", "boundaries": "observed joined repeat edges, not TSS; abs(distance)<200 bp",
              "case_selection": ("ALL AC5 count-supported multi-ALT HERV cases; no site cap; no ancient statistic or window peak used" if args.max_sites == 0 else
                                 "AC5 count-supported multi-ALT; per-family chromosome-stratified reservoir; no ancient statistic or window peak used"),
              "all_observed_alt_catalog": catalog_all,
              "catalog_scope": "All selected modern loci with all observed ALT bases (including additional ALT AC<5, labelled separately). REF is a label, not an inferred observed count. All cached eligible HERV loci only when max_sites=0; not all single-ALT genomic positions.",
              "overlapping_family_assignment": "one seeded family per locus; memberships retained in candidate_families",
              "controls": "unique loci globally; both arms AC5 multi-ALT; no aligned LTR body and >=2000 bp from ERV1/ERVK observed edges",
              "matching": "chrom, canonical trinucleotide, primary ALT AF bin, primary ALT AN bin, secondary ALT AF bin; same missing-filter signature",
              "af_definition": "highest and second-highest AC/AN among AC>=5 observed ALT alleles; no inferred REF AF",
              "estimand_warning": "Candidate follow-up conditional on contemporary multi-ALT occupancy and AF/AN matching. Not an HERV enrichment estimate, an allele-age/retention test, or an AF difference test.",
              "input_validation": "completion signatures, file sizes and modification metadata; no added full-file checksums",
              "outputs": {p: metadata(out / p) for p in output_names}}
    write_json(receipt, result)
    print(f"[herv-candidates] {len(cases)} cases, {len(used_controls)} unique controls; {out / 'targets.tsv'}", flush=True)
    return result


def write_allele_catalogs(out, targets):
    """Preserve every observed modern base and count for the complete time panel.

    A locus is selected because at least two ALTs meet AC5. Extra lower-count
    bases remain explicitly lower-count, and reference counts are never inferred
    by subtracting unrelated per-ALT denominators.
    """
    site_fields = [*FIELDS, "all_alt_ac", "all_alt_an", "catalog_scope"]
    allele_fields = ["allele_id", "site_id", "chrom", "pos", "ref", "allele", "role",
                     "modern_ac", "modern_an", "modern_af", "ac5_supported", "ac10_supported",
                     "family", "group", "match_id", "partition", "comparison_partition",
                     "candidate_families", "family_partitions", "strict_filter_complete", "missing_strict_filters"]
    site_path, allele_path = Path(out)/"site_catalog.tsv.gz", Path(out)/"allele_catalog.tsv.gz"
    temporary = [p.with_name(p.name+".partial") for p in (site_path, allele_path)]
    with gzip.open(temporary[0], "wt", newline="", compresslevel=3) as sf, gzip.open(temporary[1], "wt", newline="", compresslevel=3) as af:
        sw = csv.DictWriter(sf, site_fields, delimiter="\t", extrasaction="ignore")
        aw = csv.DictWriter(af, allele_fields, delimiter="\t", extrasaction="ignore")
        sw.writeheader(); aw.writeheader()
        for row in targets:
            alts = row["observed_alts_all"].split(",")
            ac, an = row["alt_ac_all"].split(","), row["alt_an_all"].split(",")
            if not len(alts) == len(ac) == len(an):
                raise ValueError(f"Catalog allele/count lengths disagree: {row['site_id']}")
            sw.writerow(dict(row, alts=",".join(alts), all_alt_ac=",".join(ac), all_alt_an=",".join(an),
                             catalog_scope="modern_AC5_multiALT_loci_all_observed_bases"))
            common = {k:row.get(k, ".") for k in allele_fields}
            aw.writerow(common | dict(allele_id=row["site_id"]+":"+row["ref"], allele=row["ref"], role="reference_label",
                        modern_ac="NA", modern_an="NA", modern_af="NA", ac5_supported="NA", ac10_supported="NA"))
            for base, count, denominator in zip(alts, ac, an):
                count, denominator = int(count), int(denominator)
                role = "primary" if base == row["primary_alt"] else "secondary" if base == row["secondary_alt"] else "additional_supported" if count >= 5 else "additional_low_count"
                aw.writerow(common | dict(allele_id=row["site_id"]+":"+base, allele=base, role=role,
                            modern_ac=count, modern_an=denominator, modern_af=count/denominator,
                            ac5_supported=int(count >= 5), ac10_supported=int(count >= 10)))
    for temp, destination in zip(temporary, (site_path, allele_path)):
        os.replace(temp, destination)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--modern-dir", required=True, help="Completed modern run.json plus original shards")
    p.add_argument("--erv-joined", required=True, help="GRCh38 14-column joined RepeatMasker table")
    p.add_argument("--out", required=True, help="Separate derived candidate output directory")
    p.add_argument("--max-sites", type=int, default=500, help="Maximum CASE sites PER FAMILY; 0 retains ALL qualified HERV cases")
    p.add_argument("--all-alleles-catalog", action="store_true", help="Also write full-base site_catalog.tsv.gz and allele_catalog.tsv.gz for ancient trajectories")
    p.add_argument("--controls-per-case", type=int, default=4)
    p.add_argument("--end-bp", type=int, default=200)
    p.add_argument("--body-bins", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261027)
    args = p.parse_args(argv)
    try:
        run_candidates(args)
    except (ValueError, OSError, KeyError) as e:
        p.exit(2, f"herv-candidates: {e}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
