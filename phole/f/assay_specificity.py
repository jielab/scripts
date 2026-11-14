#!/usr/bin/env python3
"""Screen assay identity hypotheses from reported pQTL loci, without diagnosing them.

A significant cis signal in one assay and its absence from another assay's
selected summary is not a statistical test of discordance. This module keeps
that distinction explicit, preserves individual probes, and exports external
cis variants for a separate full-GWAS lookup. Alternative targets are hypotheses
supported by a shared genetic fingerprint, not identifications of captured
proteins. No effect direction is inferred from these summaries.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path

DATASETS = ("ukb", "ckb", "decode", "scot")
PLATFORMS = {"ukb": "Olink", "ckb": "Olink", "decode": "SomaScan", "scot": "MS"}
CODING = {"protein_coding", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
          "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"}
SCHEMA = "phole-assay-specificity-v3"
GWAS_THRESHOLD = Decimal("5e-8")
BAD_SIGNAL_FLAGS = {"identical_alleles_in_id", "invalid_effect_or_se", "reported_nonsignificant"}
BAD_ASSAY_FLAGS = {"file_dataset_build_conflict", "target_or_build_unresolved"}
NO_CIS_STATES = {"no_gwas_cis_certified", "no_cis_reported_uncertified", "cis_search_censored", "cis_unknown"}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, fields, rows):
    path = Path(path)
    tmp = Path(str(path) + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "." if row.get(key) is None else row.get(key, ".") for key in fields})
    tmp.replace(path)


def write_json(path, value):
    path = Path(path)
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def pvalue(value):
    try:
        p = Decimal(str(value))
        if not p.is_finite() or not 0 <= p <= 1:
            raise ValueError()
        return p
    except (InvalidOperation, ValueError):
        raise ValueError(f"Invalid P value: {value}") from None


def split_flags(value):
    if isinstance(value, list):
        return set(value)
    return {part for part in str(value or "").split(";") if part and part != "."}


def chrom(value):
    c = re.sub(r"^chr", "", str(value), flags=re.I).upper()
    return {"23": "X", "24": "Y", "25": "MT", "M": "MT"}.get(c, c)


def chrom_key(value):
    return (0, int(value)) if value.isdigit() else (1, {"X": "01", "Y": "02", "MT": "03"}.get(value, value))


def identity_verified(protein):
    if "identity_verified" in protein:
        return protein["identity_verified"] is True
    coordinate = protein.get("coordinate", {})
    if "identity_verified" in coordinate:
        return coordinate["identity_verified"] is True
    target = protein.get("target")
    return bool(target and target.get("status") == "single_gene" and target.get("biotype") in CODING)


def coordinate_status(protein):
    coordinate = protein.get("coordinate", {})
    return coordinate.get("coordinate_status", coordinate.get("status", "reference_only" if identity_verified(protein) else "unresolved"))


def same_coding_gene(left, right):
    """Different assay labels/proteoforms of one coding gene are not wrong mothers."""
    def identities(protein, field):
        return {str(t[field]).split(".")[0] for t in protein.get("intervals", []) if t.get(field) not in {None, "", "."}}
    left_ids, right_ids = identities(left, "gene_id"), identities(right, "gene_id")
    if left_ids and right_ids:
        return bool(left_ids & right_ids)
    return bool(identities(left, "gene") & identities(right, "gene"))


def summary_valid(assay):
    flags = split_flags(assay.get("flags"))
    return (assay.get("status") == "complete" and not (flags & BAD_ASSAY_FLAGS)
            and not any(flag.startswith("metadata_") for flag in flags))


def assay_state(protein, dataset, assay, threshold):
    if not identity_verified(protein):
        return "coordinate_or_identity_unresolved"
    if assay.get("status") != "complete":
        return "input_" + assay.get("status", "unavailable")
    if not summary_valid(assay):
        return "summary_invalid"
    lead = assay.get("best_cis")
    if lead:
        return "cis_strict" if pvalue(lead["p"]) <= threshold else "cis_gwas_only"
    if assay.get("cis_search_complete") is True:
        return "no_gwas_cis_certified"
    flags = split_flags(assay.get("flags"))
    if "cis_search_censored_by_legacy_export" in flags:
        return "cis_search_censored"
    if assay.get("cis") == 0 or "legacy_cis_window_not_fully_certified" in flags:
        return "no_cis_reported_uncertified"
    return "cis_unknown"


def normalize_signals(rows, proteins):
    """Retain per-assay observations; a merged JSON lead cannot stand for all probes."""
    index = {}
    for protein in proteins:
        for dataset in DATASETS:
            for assay in protein["data"][dataset].get("assays", []):
                key = (protein["prot"], dataset, assay["alias"])
                if key in index:
                    raise ValueError(f"Duplicate assay identity in cross-check: {key}")
                index[key] = assay
    result = []
    seen = set()
    for row in rows:
        key = (row["prot"], row["dataset"], row["assay"])
        if key not in index:
            raise ValueError(f"signals.tsv contains an assay absent from the cross-check payload: {key}")
        signal = dict(prot=key[0], dataset=key[1], assay=key[2], chr=chrom(row["CHR"]),
                      pos=int(row["POS"]), snp=row["SNP"], p=str(row["P"]), kind=row["CLASS"],
                      build=str(row.get("build", ".")).removeprefix("GRCh"), flags=split_flags(row.get("flags")))
        signal["p_decimal"] = pvalue(signal["p"])
        if signal["pos"] < 1 or signal["kind"] not in {"cis", "trans", "unknown"}:
            raise ValueError("Invalid signal position/class in cross-check bundle")
        identity = (*key, signal["chr"], signal["pos"], signal["snp"], signal["kind"])
        # An exporter can repeat a raw row; its multiplicity is not independent evidence.
        if identity in seen:
            previous = next(s for s in result if (s["prot"], s["dataset"], s["assay"], s["chr"], s["pos"], s["snp"], s["kind"]) == identity)
            if previous["p_decimal"] != signal["p_decimal"] or previous["build"] != signal["build"]:
                raise ValueError(f"Conflicting duplicate signal: {identity}")
            previous["flags"].update(signal["flags"])
            previous["clean"] = previous["assay_valid"] and not (previous["flags"] & BAD_SIGNAL_FLAGS) and previous["build"] in {"37", "38"}
            continue
        seen.add(identity)
        signal["assay_valid"] = summary_valid(index[key])
        signal["clean"] = signal["assay_valid"] and not (signal["flags"] & BAD_SIGNAL_FLAGS) and signal["build"] in {"37", "38"}
        result.append(signal)
    by_assay = defaultdict(list)
    for signal in result:
        by_assay[(signal["prot"], signal["dataset"], signal["assay"])].append(signal)
    for key, assay in index.items():
        rows_for_assay = by_assay.get(key, [])
        cis = [s for s in rows_for_assay if s["kind"] == "cis" and s["p_decimal"] <= GWAS_THRESHOLD]
        lead = assay.get("best_cis")
        if bool(cis) != bool(lead):
            raise ValueError(f"signals.tsv and payload cis calls disagree: {key}")
        if cis and not any(s["snp"] == lead["snp"] and s["chr"] == str(lead["chr"]) and s["pos"] == int(lead["pos"])
                           and s["p_decimal"] == pvalue(lead["p"]) for s in cis):
            raise ValueError(f"signals.tsv and payload cis lead disagree: {key}")
    return result, by_assay, index


def make_hotspot_index(signals):
    exact = defaultdict(set)
    regional = defaultdict(list)
    for signal in signals:
        if not signal["assay_valid"] or signal["kind"] not in {"cis", "trans"} or signal["p_decimal"] > GWAS_THRESHOLD:
            continue
        exact[(signal["dataset"], signal["build"], signal["chr"], signal["pos"], signal["snp"])].add(signal["prot"])
        if signal["kind"] == "trans":
            regional[(signal["dataset"], signal["build"], signal["chr"])].append((signal["pos"], signal["prot"]))
    regional = {key: sorted(set(value)) for key, value in regional.items()}
    positions = {key: [v[0] for v in values] for key, values in regional.items()}

    def counts(signal, halfwidth=500000):
        key = (signal["dataset"], signal["build"], signal["chr"])
        values = regional.get(key, [])
        points = positions.get(key, [])
        left = bisect.bisect_left(points, signal["pos"] - halfwidth)
        right = bisect.bisect_right(points, signal["pos"] + halfwidth)
        names = {name for _, name in values[left:right]}
        return len(exact.get((*key, signal["pos"], signal["snp"]), set())), len(names)

    return counts


def variant_match(left, right):
    if left["build"] != right["build"] or left["chr"] != right["chr"]:
        return None
    if left["pos"] == right["pos"]:
        if left["snp"] == right["snp"] and re.fullmatch(r"rs\d+", left["snp"]):
            return "same_rsid_and_position"
        return "same_position_unverified_alleles"
    return "regional_proximity"


def profile_loci(signals, span):
    """Position clusters have bounded total span, not chaining/LD independence."""
    ordered = sorted(signals, key=lambda s: (s["build"], chrom_key(s["chr"]), s["pos"], s["snp"]))
    loci = []
    for signal in ordered:
        if (not loci or signal["build"] != loci[-1]["build"] or signal["chr"] != loci[-1]["chr"]
                or signal["pos"] - loci[-1]["start"] > span):
            loci.append(dict(build=signal["build"], chr=signal["chr"], start=signal["pos"], end=signal["pos"], signals=[signal]))
        else:
            loci[-1]["end"] = signal["pos"]
            loci[-1]["signals"].append(signal)
    return loci


def match_profiles(left, right, span, hotspot_counts):
    """Maximum-cardinality ordered one-to-one matching; then exact rsIDs/distance.

    Distinct loci within each profile are position ordered. With the bounded
    interval compatibility below, an optimal non-crossing alignment avoids
    repeatedly counting the same hotspot for several summary rows.
    """
    a, b = profile_loci(left, span), profile_loci(right, span)
    cache = {}

    def pair(i, j):
        key = (i, j)
        if key in cache:
            return cache[key]
        x, y = a[i], b[j]
        if x["build"] != y["build"] or x["chr"] != y["chr"] or max(x["end"], y["end"]) - min(x["start"], y["start"]) > span:
            cache[key] = None
            return None
        options = [(variant_match(s, t) == "same_rsid_and_position", abs(s["pos"] - t["pos"]), s, t)
                   for s in x["signals"] for t in y["signals"]]
        exact, distance, s, t = min(options, key=lambda v: (-v[0], v[1], v[2]["snp"], v[3]["snp"]))
        ecount, rcount = hotspot_counts(s)
        cache[key] = dict(left=s, right=t, exact=exact, distance=distance, nonhotspot=ecount <= 10 and rcount <= 10)
        return cache[key]

    # Scores are (cardinality, exact matches, negative distance). The fourth
    # object is the explicit matched-pair list for audit and hotspot annotation.
    table = [[((0, 0, 0), []) for _ in range(len(b) + 1)] for _ in range(len(a) + 1)]
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            options = [table[i - 1][j], table[i][j - 1]]
            hit = pair(i - 1, j - 1)
            if hit is not None:
                score, pairs = table[i - 1][j - 1]
                options.append(((score[0] + 1, score[1] + int(hit["exact"]), score[2] - hit["distance"]), pairs + [hit]))
            table[i][j] = max(options, key=lambda option: option[0])
    score, pairs = table[-1][-1]
    denominator = len(a) + len(b) - score[0]
    return dict(left_n=len(a), right_n=len(b), shared=score[0], exact=score[1],
                nonhotspot=sum(hit["nonhotspot"] for hit in pairs),
                jaccard=score[0] / denominator if denominator else 0.0,
                regions=";".join(f'{hit["left"]["chr"]}:{hit["left"]["pos"]}~{hit["right"]["pos"]}' for hit in pairs) or ".",
                matches=[dict(BUILD=hit["left"]["build"], CHR=hit["left"]["chr"], PRIMARY_POS=hit["left"]["pos"],
                              PRIMARY_SNP=hit["left"]["snp"], PRIMARY_P=hit["left"]["p"], PRIMARY_CLASS=hit["left"]["kind"],
                              REFERENCE_POS=hit["right"]["pos"], REFERENCE_SNP=hit["right"]["snp"], REFERENCE_P=hit["right"]["p"],
                              REFERENCE_CLASS=hit["right"]["kind"], EXACT_RSID=int(hit["exact"]), DISTANCE_BP=hit["distance"],
                              NONHOTSPOT=int(hit["nonhotspot"])) for hit in pairs])


def analyse(payload, signal_rows, primary="all", match_span=1000000):
    if primary not in {*DATASETS, "all"} or match_span < 0:
        raise ValueError("Invalid primary dataset or match span")
    proteins = payload["proteins"]
    lookup = {protein["prot"]: protein for protein in proteins}
    if len(lookup) != len(proteins):
        raise ValueError("Duplicate canonical protein in cross-check")
    thresholds = {d: pvalue(payload["summary"]["strong_thresholds"][d]) for d in DATASETS}
    if any(not 0 < threshold <= GWAS_THRESHOLD for threshold in thresholds.values()):
        raise ValueError("Invalid strict threshold in cross-check")
    signals, by_assay, assay_index = normalize_signals(signal_rows, proteins)
    count_hotspot = make_hotspot_index(signals)
    selected = DATASETS if primary == "all" else (primary,)
    strict_by_protein = defaultdict(list)
    clean_profiles = defaultdict(list)
    for signal in signals:
        if signal["p_decimal"] <= thresholds[signal["dataset"]] and signal["assay_valid"]:
            if signal["kind"] == "cis":
                strict_by_protein[(signal["prot"], signal["dataset"])].append(signal)
            if signal["clean"] and signal["kind"] in {"cis", "trans"}:
                clean_profiles[(signal["prot"], signal["dataset"], signal["assay"])].append(signal)
    discordance, triage, queries, excluded_queries = [], [], [], []
    candidate_assays = set()
    profile_cache = {}

    def profile_compare(left_key, right_key):
        key = (left_key, right_key)
        if key not in profile_cache:
            profile_cache[key] = match_profiles(clean_profiles.get(left_key, []), clean_profiles.get(right_key, []), match_span, count_hotspot)
        return profile_cache[key]

    for protein in proteins:
        name = protein["prot"]
        states = {d: [assay_state(protein, d, a, thresholds[d]) for a in protein["data"][d].get("assays", [])] for d in DATASETS}
        external_strict = {d: strict_by_protein.get((name, d), []) for d in DATASETS}
        for dataset in selected:
            for assay in protein["data"][dataset].get("assays", []):
                state = assay_state(protein, dataset, assay, thresholds[dataset])
                if state == "cis_strict":
                    continue
                key = (name, dataset, assay["alias"])
                all_trans = [s for s in by_assay.get(key, []) if s["kind"] == "trans" and s["p_decimal"] <= GWAS_THRESHOLD]
                clean_trans = [s for s in clean_profiles.get(key, []) if s["kind"] == "trans"]
                counts = [count_hotspot(s) for s in all_trans]
                exact_only_generic = bool(counts) and all(x > 10 for x, _ in counts)
                region_only_generic = bool(counts) and all(y > 10 for _, y in counts)
                other_strong = [d for d in DATASETS if d != dataset and PLATFORMS[d] != PLATFORMS[dataset] and external_strict[d]]
                clean_other = [d for d in other_strong if any(s["clean"] for s in external_strict[d])]
                primary_cis_assays = sum(s in {"cis_strict", "cis_gwas_only"} for s in states[dataset])
                same_platform = [d for d in DATASETS if d != dataset and PLATFORMS[d] == PLATFORMS[dataset] and external_strict[d]]
                same_platform_gap = [d for d in DATASETS if d != dataset and PLATFORMS[d] == PLATFORMS[dataset] and states[d]
                                     and all(s in NO_CIS_STATES for s in states[d])]
                context = ("other_same_platform_cohort_has_cis" if same_platform else
                           "SomaScan_and_MS_support_target" if "decode" in clean_other and "scot" in clean_other else
                           "MS_supports_target" if "scot" in clean_other else
                           "different_binding_platform_supports_target" if clean_other else
                           "reference_signal_QC_review" if other_strong else "no_strong_orthogonal_cis_observed")
                base = dict(PROT=name, PRIMARY=dataset, PRIMARY_ASSAY=assay["alias"], PRIMARY_PLATFORM=PLATFORMS[dataset],
                            PRIMARY_STATE=state, PRIMARY_CIS_P=(assay.get("best_cis") or {}).get("p"),
                            PRIMARY_CIS_SEARCH_COMPLETE=int(assay.get("cis_search_complete") is True),
                            PRIMARY_ASSAYS_EXPECTED=len(states[dataset]), PRIMARY_ASSAYS_WITH_CIS=primary_cis_assays,
                            PRIMARY_STRICT_TRANS_LOCI=len(profile_loci(clean_trans, match_span)), PRIMARY_REPORTED_TRANS_ROWS=len(all_trans),
                            ALL_TRANS_EXACT_NON_SPECIFIC=int(exact_only_generic), ALL_TRANS_REGIONAL_HOTSPOTS=int(region_only_generic),
                            STRONG_ORTHOGONAL_CIS_COHORTS=";".join(other_strong) or ".", CLEAN_ORTHOGONAL_CIS_COHORTS=";".join(clean_other) or ".",
                            SAME_PLATFORM_CIS_COHORTS=";".join(same_platform) or ".", SAME_PLATFORM_UNREPORTED_CIS_COHORTS=";".join(same_platform_gap) or ".",
                            COORDINATE_STATUS=coordinate_status(protein), TARGET_IDENTITY_VERIFIED=int(identity_verified(protein)), CONTEXT=context,
                            EVIDENCE_SCOPE="reported_summary_screen", WRONG_PROTEIN_CONFIRMED=0,
                            NONREPLICATION_TESTED=0, FLAGS=";".join(sorted(split_flags(assay.get("flags")))) or ".")
                triage.append(base)
                for reference in other_strong:
                    donor_signals = external_strict[reference]
                    clean_donors = [s for s in donor_signals if s["clean"]]
                    best = min(clean_donors or donor_signals, key=lambda s: (s["p_decimal"], s["assay"], s["pos"], s["snp"]))
                    raw_best = min(donor_signals, key=lambda s: (s["p_decimal"], s["assay"], s["pos"], s["snp"]))
                    if state not in NO_CIS_STATES | {"cis_gwas_only"}:
                        screen = "unassessable_primary"
                    elif not clean_donors:
                        screen = "reference_QC_review"
                    elif state == "cis_gwas_only":
                        screen = "cis_strength_difference"
                    elif state == "no_gwas_cis_certified":
                        screen = "certified_threshold_discordance"
                    else:
                        screen = "reported_cis_discordance"
                    row = {**base, "REFERENCE": reference, "REFERENCE_PLATFORM": PLATFORMS[reference], "REFERENCE_CIS_P": best["p"],
                           "REFERENCE_CIS_SNP": best["snp"], "REFERENCE_CIS_CHR": best["chr"], "REFERENCE_CIS_POS": best["pos"],
                           "REFERENCE_STRICT_ASSAYS": ";".join(sorted({s["assay"] for s in donor_signals})),
                           "REFERENCE_ASSAYS_EXPECTED": len(states[reference]), "REFERENCE_CLEAN_STRICT_ASSAYS": len({s["assay"] for s in clean_donors}),
                           "REFERENCE_RAW_BEST_CIS_P": raw_best["p"], "REFERENCE_RAW_BEST_CIS_SNP": raw_best["snp"],
                           "REFERENCE_RAW_BEST_CIS_FLAGS": ";".join(sorted(raw_best["flags"])) or ".",
                           "REFERENCE_FLAGS": ";".join(sorted({flag for s in donor_signals for flag in s["flags"]})) or ".", "SCREEN_CLASS": screen}
                    discordance.append(row)
                    if state in NO_CIS_STATES and clean_donors:
                        candidate_assays.add(key)
                    if state not in NO_CIS_STATES | {"cis_gwas_only"}:
                        continue
                    for signal in donor_signals:
                        identity = (name, dataset, assay["alias"], reference, signal["assay"], signal["build"], signal["chr"], signal["pos"], signal["snp"])
                        query = dict(QUERY_ID="SCIS_" + hashlib.sha256(json.dumps(identity, separators=(",", ":")).encode()).hexdigest()[:20],
                                     PROT=name, PRIMARY=dataset, PRIMARY_DATASET=dataset, PRIMARY_ASSAY=assay["alias"],
                                     REFERENCE=reference, REFERENCE_DATASET=reference, REFERENCE_ASSAY=signal["assay"],
                                     BUILD=signal["build"], CHR=signal["chr"], POS=signal["pos"], SNP=signal["snp"],
                                     REFERENCE_P=signal["p"], REFERENCE_THRESHOLD=str(thresholds[reference]), QUERY_TYPE="own_cis_external",
                                     COORDINATE_STATUS=coordinate_status(protein), PRIMARY_STATE=state,
                                     REFERENCE_FLAGS=";".join(sorted(signal["flags"])) or ".", QUERY_ELIGIBLE=int(signal["clean"]),
                                     EXCLUSION_REASON="." if signal["clean"] else "reference_signal_or_summary_QC")
                        (queries if signal["clean"] else excluded_queries).append(query)

    # Potential alternative targets must be individually mapped coding genes
    # with a strong cis association on a different measurement platform.
    intervals_by_chr = defaultdict(list)
    flank = int(payload["config"]["cis_flank"])
    for protein in proteins:
        if not identity_verified(protein):
            continue
        for interval in protein.get("intervals", []):
            if interval.get("biotype") not in CODING:
                continue
            intervals_by_chr[(str(interval.get("build", "38")), interval["chr"])].append(
                (max(1, int(interval["start"]) - flank), int(interval["end"]) + flank, protein["prot"]))
    alternatives = []
    for key in sorted(candidate_assays):
        name, dataset, assay = key
        protein = lookup[name]
        clean_trans = [s for s in clean_profiles.get(key, []) if s["kind"] == "trans"]
        for signal in clean_trans:
            other_names = sorted({target for start, end, target in intervals_by_chr.get((signal["build"], signal["chr"]), [])
                                  if start <= signal["pos"] <= end and target != name})
            for target in other_names:
                other_protein = lookup[target]
                if same_coding_gene(protein, other_protein):
                    continue
                for reference in DATASETS:
                    if PLATFORMS[reference] == PLATFORMS[dataset]:
                        continue
                    cis_hits = [s for s in strict_by_protein.get((target, reference), []) if s["clean"]
                                and s["chr"] == signal["chr"] and s["build"] == signal["build"] and abs(s["pos"] - signal["pos"]) <= match_span]
                    # Each alternative reagent is a separate hypothesis. Retain
                    # its best matching cis anchor, without combining probes.
                    by_probe = defaultdict(list)
                    for hit in cis_hits:
                        by_probe[hit["assay"]].append(hit)
                    for donor_assay, hits in sorted(by_probe.items()):
                        hit = min(hits, key=lambda s: (variant_match(signal, s) != "same_rsid_and_position",
                                                      s["pos"] != signal["pos"], abs(s["pos"] - signal["pos"]), s["p_decimal"], s["snp"]))
                        relation = variant_match(signal, hit)
                        exact_count, regional_count = count_hotspot(signal)
                        alternative_profile = profile_compare(key, (target, reference, donor_assay))
                        own_profiles = [(a["alias"], profile_compare(key, (name, reference, a["alias"])))
                                        for a in protein["data"][reference].get("assays", []) if summary_valid(a)]
                        # The strongest matching intended-target probe is a
                        # conservative comparator against an alternative-target claim.
                        own_assay, own_profile = max(own_profiles, key=lambda value: (value[1]["shared"], value[1]["jaccard"], value[1]["exact"], value[0])) if own_profiles else (
                            ".", dict(right_n=0, shared=0, exact=0, jaccard=0.0))
                        own_status = "not_measured_or_unusable" if not own_profiles else "no_clean_strict_loci" if not own_profile["right_n"] else "observed"
                        other_target = other_protein.get("target") or {}
                        my_target = protein.get("target") or {}
                        nearby = (other_target.get("chr") == my_target.get("chr") and other_target.get("chr") is not None
                                  and abs((other_target.get("start", 0) + other_target.get("end", 0)) / 2
                                          - (my_target.get("start", 0) + my_target.get("end", 0)) / 2) <= match_span)
                        flags = []
                        if exact_count > 10:
                            flags.append("exact_variant_pleiotropic")
                        if regional_count > 10:
                            flags.append("regional_hotspot")
                        if nearby:
                            flags.append("neighboring_coding_targets")
                        same_target_trans = []
                        for d in DATASETS:
                            if d == dataset or PLATFORMS[d] == PLATFORMS[dataset]:
                                continue
                            if any(t["kind"] == "trans" and t["chr"] == signal["chr"] and t["build"] == signal["build"]
                                   and abs(t["pos"] - signal["pos"]) <= match_span
                                   for a in protein["data"][d].get("assays", []) for t in clean_profiles.get((name, d, a["alias"]), [])):
                                same_target_trans.append(d)
                        if same_target_trans:
                            flags.append("intended_target_trans_also_reproduced")
                        target_ref = other_protein.get("coordinate", {})
                        shared_anchor = target_ref.get("shared_anchor_targets", "")
                        if shared_anchor and shared_anchor != ".":
                            flags.append("alternative_coordinate_shared_with_other_targets")
                        profile_more = (own_status == "observed" and alternative_profile["shared"] > own_profile["shared"]
                                        and alternative_profile["jaccard"] > own_profile["jaccard"])
                        alternatives.append(dict(
                            HYPOTHESIS_ID="ALT_" + hashlib.sha256(json.dumps((name, dataset, assay, target, reference, donor_assay), separators=(",", ":")).encode()).hexdigest()[:20],
                            PROT=name, PRIMARY=dataset, PRIMARY_ASSAY=assay, ALTERNATIVE_PROT=target, REFERENCE=reference, REFERENCE_ASSAY=donor_assay,
                            BUILD=signal["build"], CHR=signal["chr"], PRIMARY_POS=signal["pos"], PRIMARY_SNP=signal["snp"], PRIMARY_P=signal["p"],
                            ALTERNATIVE_CIS_POS=hit["pos"], ALTERNATIVE_CIS_SNP=hit["snp"], ALTERNATIVE_CIS_P=hit["p"],
                            LEAD_MATCH=relation, LEAD_DISTANCE_BP=abs(signal["pos"] - hit["pos"]),
                            PRIMARY_EXACT_VARIANT_TARGETS=exact_count, PRIMARY_REGION_TRANS_TARGETS=regional_count,
                            PROFILE_PRIMARY_LOCI=alternative_profile["left_n"], PROFILE_ALTERNATIVE_LOCI=alternative_profile["right_n"],
                            PROFILE_SHARED_REGIONS=alternative_profile["shared"], PROFILE_SHARED_EXACT_RSID=alternative_profile["exact"],
                            PROFILE_SHARED_NONHOTSPOT=alternative_profile["nonhotspot"], PROFILE_JACCARD=alternative_profile["jaccard"],
                            PROFILE_MATCHED_REGIONS=alternative_profile["regions"], OWN_REFERENCE_ASSAY=own_assay, OWN_PROFILE_STATUS=own_status,
                            OWN_PROFILE_LOCI=own_profile["right_n"], OWN_PROFILE_SHARED_REGIONS=own_profile["shared"],
                            OWN_PROFILE_SHARED_EXACT_RSID=own_profile["exact"], OWN_PROFILE_JACCARD=own_profile["jaccard"],
                            PROFILE_MORE_SIMILAR_TO_ALTERNATIVE=int(profile_more), INTENDED_TARGET_TRANS_REFERENCE_COHORTS=";".join(same_target_trans) or ".",
                            ALTERNATIVE_MS_CIS=int(any(s["clean"] for s in strict_by_protein.get((target, "scot"), []))),
                            FLAGS=";".join(flags) or ".", EVIDENCE_SCOPE="selected_locus_fingerprint", DIRECTION_TESTED=0, COLOCALIZATION_TESTED=0, WRONG_PROTEIN_CONFIRMED=0))
    fingerprint_matches = []
    fingerprint_seen = set()
    for row in alternatives:
        if row["HYPOTHESIS_ID"] in fingerprint_seen:
            continue
        fingerprint_seen.add(row["HYPOTHESIS_ID"])
        key = (row["PROT"], row["PRIMARY"], row["PRIMARY_ASSAY"])
        for comparison, target, assay in (("alternative", row["ALTERNATIVE_PROT"], row["REFERENCE_ASSAY"]),
                                           ("intended", row["PROT"], row["OWN_REFERENCE_ASSAY"])):
            if assay == ".":
                continue
            profile = profile_compare(key, (target, row["REFERENCE"], assay))
            for match in profile["matches"]:
                fingerprint_matches.append(dict(HYPOTHESIS_ID=row["HYPOTHESIS_ID"], PROT=row["PROT"], PRIMARY=row["PRIMARY"],
                                                PRIMARY_ASSAY=row["PRIMARY_ASSAY"], COMPARISON=comparison, REFERENCE_PROT=target,
                                                REFERENCE=row["REFERENCE"], REFERENCE_ASSAY=assay, **match))
    ambiguity = defaultdict(set)
    exact_ambiguity = defaultdict(set)
    for row in alternatives:
        key = (row["PROT"], row["PRIMARY"], row["PRIMARY_ASSAY"], row["CHR"], row["PRIMARY_POS"], row["PRIMARY_SNP"])
        ambiguity[key].add(row["ALTERNATIVE_PROT"])
        if row["LEAD_MATCH"] == "same_rsid_and_position":
            exact_ambiguity[key].add(row["ALTERNATIVE_PROT"])
    for row in alternatives:
        key = (row["PROT"], row["PRIMARY"], row["PRIMARY_ASSAY"], row["CHR"], row["PRIMARY_POS"], row["PRIMARY_SNP"])
        row["ALTERNATIVE_TARGETS_AT_LOCUS"] = len(ambiguity[key])
        row["EXACT_ALTERNATIVE_TARGETS_AT_VARIANT"] = len(exact_ambiguity[key])
        specific = (row["LEAD_MATCH"] == "same_rsid_and_position" and row["PRIMARY_EXACT_VARIANT_TARGETS"] <= 10
                    and row["EXACT_ALTERNATIVE_TARGETS_AT_VARIANT"] == 1 and row["INTENDED_TARGET_TRANS_REFERENCE_COHORTS"] == "."
                    and not (split_flags(row["FLAGS"]) & {"neighboring_coding_targets", "alternative_coordinate_shared_with_other_targets"}))
        fingerprint = (specific and row["PROFILE_SHARED_REGIONS"] >= 2 and row["PROFILE_SHARED_NONHOTSPOT"] >= 1
                       and row["PROFILE_MORE_SIMILAR_TO_ALTERNATIVE"])
        regional_caution = row["PRIMARY_REGION_TRANS_TARGETS"] > 10 or row["ALTERNATIVE_TARGETS_AT_LOCUS"] > 1
        row["FOLLOW_UP_CLASS"] = ("specific_fingerprint_with_regional_caution" if fingerprint and regional_caution else
                                  "specific_anchor_and_fingerprint" if fingerprint else
                                  "specific_anchor_only" if specific else "hotspot_or_ambiguous_anchor" if row["FLAGS"] != "." or row["ALTERNATIVE_TARGETS_AT_LOCUS"] > 1 else "regional_only")
        row["INTERPRETATION"] = ("targeted_validation_priority" if specific else
                                 "shared_regulatory_pattern" if row["INTENDED_TARGET_TRANS_REFERENCE_COHORTS"] != "." else "regional_overlap_hypothesis")
    fingerprint_classes = {"specific_anchor_and_fingerprint", "specific_fingerprint_with_regional_caution"}
    order = {"specific_anchor_and_fingerprint": 0, "specific_fingerprint_with_regional_caution": 1,
             "specific_anchor_only": 2, "regional_only": 3, "hotspot_or_ambiguous_anchor": 4}
    alternatives.sort(key=lambda r: (order[r["FOLLOW_UP_CLASS"]], -r["PROFILE_SHARED_NONHOTSPOT"], -r["PROFILE_SHARED_REGIONS"],
                                     -r["PROFILE_JACCARD"], r["PROT"], r["PRIMARY"], r["PRIMARY_ASSAY"], r["ALTERNATIVE_PROT"], r["REFERENCE_ASSAY"]))
    discordance.sort(key=lambda r: (r["PRIMARY"], r["PROT"], r["PRIMARY_ASSAY"], r["REFERENCE"]))
    queries.sort(key=lambda r: (r["PRIMARY"], r["PROT"], r["PRIMARY_ASSAY"], r["REFERENCE"], r["REFERENCE_ASSAY"], chrom_key(r["CHR"]), r["POS"], r["SNP"]))
    excluded_queries.sort(key=lambda r: r["QUERY_ID"])
    canonical = []
    grouped = defaultdict(list)
    for row in discordance:
        grouped[(row["PROT"], row["PRIMARY"])].append(row)
    for (name, dataset), rows in sorted(grouped.items()):
        data = lookup[name]["data"][dataset]
        all_states = [assay_state(lookup[name], dataset, a, thresholds[dataset]) for a in data["assays"]]
        strict_refs = sorted({row["REFERENCE"] for row in rows})
        clean_refs = sorted({row["REFERENCE"] for row in rows if row["REFERENCE_CLEAN_STRICT_ASSAYS"] > 0})
        eligible = [row for row in rows if row["SCREEN_CLASS"] in {"reported_cis_discordance", "certified_threshold_discordance"}]
        candidates = [row for row in alternatives if row["PROT"] == name and row["PRIMARY"] == dataset]
        canonical.append(dict(PROT=name, PRIMARY=dataset, PRIMARY_CANONICAL_CIS=data.get("cis"), PRIMARY_STATUS=data.get("status"),
                              ALL_PRIMARY_ASSAYS_NO_REPORTED_CIS=int(bool(all_states) and all(s in NO_CIS_STATES for s in all_states)),
                              ALL_PRIMARY_CIS_SEARCHES_COMPLETE=int(bool(data["assays"]) and all(a.get("cis_search_complete") is True for a in data["assays"])),
                              PRIMARY_ASSAY_STATES=";".join(f'{a["alias"]}:{s}' for a, s in zip(data["assays"], all_states)),
                              STRONG_REFERENCE_COHORTS=";".join(strict_refs), CLEAN_STRONG_REFERENCE_COHORTS=";".join(clean_refs) or ".",
                              REFERENCE_QC_REVIEW_ONLY=int(bool(strict_refs) and not clean_refs),
                              N_DISCORDANT_PRIMARY_ASSAYS=len({row["PRIMARY_ASSAY"] for row in eligible}), MS_CIS_SUPPORT=int("scot" in clean_refs),
                              SAME_PLATFORM_CIS_COHORTS=rows[0]["SAME_PLATFORM_CIS_COHORTS"],
                              N_ALTERNATIVE_PROTEINS=len({row["ALTERNATIVE_PROT"] for row in candidates}),
                              SPECIFIC_FINGERPRINT_ALTERNATIVES=";".join(sorted({row["ALTERNATIVE_PROT"] for row in candidates if row["FOLLOW_UP_CLASS"] in fingerprint_classes})) or ".",
                              COORDINATE_STATUS=coordinate_status(lookup[name]), WRONG_PROTEIN_CONFIRMED=0))
    summary = dict(schema=SCHEMA, primary=primary, evidence_scope="reported_locus_hypothesis_screen",
                   counts=dict(proteins=len(proteins), assay_discordance_rows=len(discordance), canonical_contrasts=len(canonical),
                               broad_canonical_no_reported_cis_with_orthogonal_cis=sum(row["ALL_PRIMARY_ASSAYS_NO_REPORTED_CIS"] for row in canonical),
                               no_cis_reported_primary_assays_with_clean_orthogonal_cis=len(candidate_assays),
                               canonical_no_cis_reported_with_orthogonal_cis=sum(row["ALL_PRIMARY_ASSAYS_NO_REPORTED_CIS"] and bool(row["N_DISCORDANT_PRIMARY_ASSAYS"]) for row in canonical),
                               specificity_queries=len(queries), excluded_queries=len(excluded_queries), alternative_anchor_rows=len(alternatives),
                               specific_fingerprint_hypotheses=len({r["HYPOTHESIS_ID"] for r in alternatives if r["FOLLOW_UP_CLASS"] in fingerprint_classes}),
                               confirmed_wrong_proteins=0),
                   by_primary={d: dict(broad_no_reported_cis_with_orthogonal_cis=sum(r["PRIMARY"] == d and r["ALL_PRIMARY_ASSAYS_NO_REPORTED_CIS"] for r in canonical),
                                       no_reported_cis_with_clean_orthogonal_cis=sum(r["PRIMARY"] == d and r["ALL_PRIMARY_ASSAYS_NO_REPORTED_CIS"] and bool(r["N_DISCORDANT_PRIMARY_ASSAYS"]) for r in canonical),
                                       discordance_rows=sum(r["PRIMARY"] == d for r in discordance),
                                       queries=sum(r["PRIMARY"] == d for r in queries)) for d in selected},
                   rules=dict(primary_cis_absence="requires explicit cis_search_complete for certified threshold absence; otherwise reported/unknown",
                              strong_p_thresholds={d: str(thresholds[d]) for d in DATASETS}, profile_locus_span_bp=match_span,
                              profile_matching="one-to-one bounded-span ordered matching; cardinality then exact-rsID then distance",
                              nonspecific_exact_variant_targets_gt=10, regional_hotspot_halfwidth_bp=500000,
                              regional_hotspot_trans_targets_gt=10, reference_cohorts_counted_once=True),
                   caveats=["Significance in one assay and nonsignificance/unreported in another is not a difference test.",
                            "SomaScan is also affinity based; its cis signal may reflect a protein-altering/epitope effect.",
                            "A cis association supports target identity but does not rule out cross-reactivity or identify causal abundance regulation.",
                            "Alternative targets and Jaccard scores use selected loci with unknown discovery completeness, not full genetic correlation.",
                            "Shared lead variants and local fingerprints can reflect true trans regulation, protein complexes, LD or technical effects.",
                            "No beta orientation, colocalization, power-adjusted absence or experimental capture identity is tested here.",
                            "Specificity queries are selected by external positive findings; follow-up is conditional, not a new unbiased discovery scan."])
    return dict(summary=summary, discordance=discordance, triage=triage, canonical=canonical, alternatives=alternatives,
                fingerprint_matches=fingerprint_matches, queries=queries, excluded_queries=excluded_queries,
                exclusions=[row for row in discordance if row["SCREEN_CLASS"] in {"reference_QC_review", "unassessable_primary"}])


def load_bundle(cross_dir, protein_map, verify_source_files=False):
    cross_dir, protein_map = Path(cross_dir), Path(protein_map)
    files = {name: cross_dir / name for name in ("cross-check.data.json", "provenance.json", "summary.json", "signals.tsv", "input-files.tsv")}
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise ValueError("Incomplete cross-check bundle: " + ", ".join(missing))
    payload = json.loads(files["cross-check.data.json"].read_text())
    provenance = json.loads(files["provenance.json"].read_text())
    summary = json.loads(files["summary.json"].read_text())
    if digest(protein_map) != provenance.get("mapping_sha256"):
        raise ValueError("Protein map SHA256 does not match the cross-check input; rerun cross-check with this map")
    if payload.get("generated") != provenance.get("generated") or payload.get("summary") != summary:
        raise ValueError("Mixed or stale cross-check bundle: generated receipt/summary differs")
    if payload.get("config", {}).get("strong_thresholds") != summary.get("strong_thresholds"):
        raise ValueError("Mixed or stale cross-check strict thresholds")
    mapping = read_tsv(protein_map)
    if not mapping or {r["prot"] for r in mapping} != {p["prot"] for p in payload["proteins"]} or len(mapping) != len(payload["proteins"]):
        raise ValueError("Protein map identities/row count differ from cross-check")
    if {"CHR", "POS"}.issubset(mapping[0]) and any("coordinate" not in p for p in payload["proteins"]):
        raise ValueError("Coordinate-bearing map requires a new coordinate-aware cross-check; old outputs are not compatible")
    if summary.get("counts", {}).get("proteins") != len(payload["proteins"]):
        raise ValueError("Cross-check protein-count receipt mismatch")
    for name, expected in provenance.get("output_sha256", {}).items():
        path = cross_dir / name
        if Path(name).is_absolute() or ".." in Path(name).parts or not path.is_file() or digest(path) != expected:
            raise ValueError(f"Cross-check output checksum mismatch: {name}")
    verified = 0
    if verify_source_files:
        for row in read_tsv(files["input-files.tsv"]):
            path = Path(row["file"])
            if not path.is_absolute():
                # Existing cross-checks may have stored cwd-relative paths.
                # Resolve only the recorded cwd path or an unambiguous root path.
                candidates = [path, Path(payload["root"]) / {"ukb": "prot", "ckb": "ckb", "decode": "decode", "scot": "scot"}[row["dataset"]] / "mplot" / (row["assay"] + ".sig.txt")]
                present = {p.resolve() for p in candidates if p.is_file()}
                if len(present) != 1:
                    raise ValueError(f"Cannot uniquely locate source summary for verification: {row['dataset']}:{row['assay']}")
                path = next(iter(present))
            if not path.is_file() or digest(path) != row.get("sha256"):
                raise ValueError(f"Source summary changed after cross-check: {row['dataset']}:{row['assay']}")
            verified += 1
    hashes = {name: digest(path) for name, path in files.items()}
    hashes["protein_map"] = digest(protein_map)
    return payload, read_tsv(files["signals.tsv"]), dict(files=hashes, verified_source_summaries=verified,
                                                        source_summary_verification=bool(verify_source_files), cross_check_generated=payload["generated"])


DISCORDANCE_FIELDS = ["PROT", "PRIMARY", "PRIMARY_ASSAY", "PRIMARY_PLATFORM", "PRIMARY_STATE", "PRIMARY_CIS_P",
                     "PRIMARY_CIS_SEARCH_COMPLETE", "PRIMARY_ASSAYS_EXPECTED", "PRIMARY_ASSAYS_WITH_CIS", "PRIMARY_STRICT_TRANS_LOCI",
                     "PRIMARY_REPORTED_TRANS_ROWS", "ALL_TRANS_EXACT_NON_SPECIFIC", "ALL_TRANS_REGIONAL_HOTSPOTS", "REFERENCE", "REFERENCE_PLATFORM",
                     "REFERENCE_CIS_P", "REFERENCE_CIS_SNP", "REFERENCE_CIS_CHR", "REFERENCE_CIS_POS", "REFERENCE_STRICT_ASSAYS",
                     "REFERENCE_ASSAYS_EXPECTED", "REFERENCE_CLEAN_STRICT_ASSAYS", "REFERENCE_FLAGS", "SCREEN_CLASS", "CONTEXT",
                     "REFERENCE_RAW_BEST_CIS_P", "REFERENCE_RAW_BEST_CIS_SNP", "REFERENCE_RAW_BEST_CIS_FLAGS",
                     "STRONG_ORTHOGONAL_CIS_COHORTS", "CLEAN_ORTHOGONAL_CIS_COHORTS", "SAME_PLATFORM_CIS_COHORTS",
                     "SAME_PLATFORM_UNREPORTED_CIS_COHORTS", "COORDINATE_STATUS", "TARGET_IDENTITY_VERIFIED", "FLAGS",
                     "EVIDENCE_SCOPE", "NONREPLICATION_TESTED", "WRONG_PROTEIN_CONFIRMED"]
TRIAGE_FIELDS = [x for x in DISCORDANCE_FIELDS if not x.startswith("REFERENCE") and x != "SCREEN_CLASS"]
QUERY_FIELDS = ["QUERY_ID", "PROT", "PRIMARY", "PRIMARY_DATASET", "PRIMARY_ASSAY", "REFERENCE", "REFERENCE_DATASET", "REFERENCE_ASSAY",
                "BUILD", "CHR", "POS", "SNP", "REFERENCE_P", "REFERENCE_THRESHOLD", "QUERY_TYPE", "COORDINATE_STATUS", "PRIMARY_STATE",
                "REFERENCE_FLAGS", "QUERY_ELIGIBLE", "EXCLUSION_REASON"]
CANONICAL_FIELDS = ["PROT", "PRIMARY", "PRIMARY_CANONICAL_CIS", "PRIMARY_STATUS", "ALL_PRIMARY_ASSAYS_NO_REPORTED_CIS",
                    "ALL_PRIMARY_CIS_SEARCHES_COMPLETE", "PRIMARY_ASSAY_STATES", "STRONG_REFERENCE_COHORTS", "CLEAN_STRONG_REFERENCE_COHORTS",
                    "REFERENCE_QC_REVIEW_ONLY", "N_DISCORDANT_PRIMARY_ASSAYS", "MS_CIS_SUPPORT", "SAME_PLATFORM_CIS_COHORTS", "N_ALTERNATIVE_PROTEINS",
                    "SPECIFIC_FINGERPRINT_ALTERNATIVES", "COORDINATE_STATUS", "WRONG_PROTEIN_CONFIRMED"]
ALTERNATIVE_FIELDS = ["HYPOTHESIS_ID", "PROT", "PRIMARY", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "REFERENCE", "REFERENCE_ASSAY", "FOLLOW_UP_CLASS", "INTERPRETATION",
                      "BUILD", "CHR", "PRIMARY_POS", "PRIMARY_SNP", "PRIMARY_P", "ALTERNATIVE_CIS_POS", "ALTERNATIVE_CIS_SNP", "ALTERNATIVE_CIS_P",
                      "LEAD_MATCH", "LEAD_DISTANCE_BP", "ALTERNATIVE_TARGETS_AT_LOCUS", "EXACT_ALTERNATIVE_TARGETS_AT_VARIANT", "PRIMARY_EXACT_VARIANT_TARGETS", "PRIMARY_REGION_TRANS_TARGETS",
                      "PROFILE_PRIMARY_LOCI", "PROFILE_ALTERNATIVE_LOCI", "PROFILE_SHARED_REGIONS", "PROFILE_SHARED_EXACT_RSID",
                      "PROFILE_SHARED_NONHOTSPOT", "PROFILE_JACCARD", "PROFILE_MATCHED_REGIONS", "OWN_REFERENCE_ASSAY", "OWN_PROFILE_STATUS", "OWN_PROFILE_LOCI",
                      "OWN_PROFILE_SHARED_REGIONS", "OWN_PROFILE_SHARED_EXACT_RSID", "OWN_PROFILE_JACCARD", "PROFILE_MORE_SIMILAR_TO_ALTERNATIVE",
                      "INTENDED_TARGET_TRANS_REFERENCE_COHORTS", "ALTERNATIVE_MS_CIS", "FLAGS", "EVIDENCE_SCOPE", "DIRECTION_TESTED",
                      "COLOCALIZATION_TESTED", "WRONG_PROTEIN_CONFIRMED"]
FINGERPRINT_FIELDS = ["HYPOTHESIS_ID", "PROT", "PRIMARY", "PRIMARY_ASSAY", "COMPARISON", "REFERENCE_PROT", "REFERENCE", "REFERENCE_ASSAY",
                      "BUILD", "CHR", "PRIMARY_POS", "PRIMARY_SNP", "PRIMARY_P", "PRIMARY_CLASS", "REFERENCE_POS", "REFERENCE_SNP", "REFERENCE_P",
                      "REFERENCE_CLASS", "EXACT_RSID", "DISTANCE_BP", "NONHOTSPOT"]


def run(cross_dir, protein_map, out, primary="all", match_span=1000000, verify_source_files=False):
    payload, signals, receipt = load_bundle(cross_dir, protein_map, verify_source_files)
    result = analyse(payload, signals, primary, match_span)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for name, fields, key in (("assay-discordance.tsv", DISCORDANCE_FIELDS, "discordance"),
                              ("assay-triage.tsv", TRIAGE_FIELDS, "triage"),
                              ("protein-discordance.tsv", CANONICAL_FIELDS, "canonical"),
                              ("alternative-targets.tsv", ALTERNATIVE_FIELDS, "alternatives"),
                              ("fingerprint-matches.tsv", FINGERPRINT_FIELDS, "fingerprint_matches"),
                              ("candidate-exclusions.tsv", DISCORDANCE_FIELDS, "exclusions"),
                              ("specificity-queries.tsv", QUERY_FIELDS, "queries"),
                              ("specificity-queries.excluded.tsv", QUERY_FIELDS, "excluded_queries")):
        write_tsv(out / name, fields, result[key])
    result["summary"].update(generated=datetime.now(timezone.utc).isoformat(timespec="seconds"), input_provenance=receipt)
    result["summary"]["output_sha256"] = {path.name: digest(path) for path in sorted(out.glob("*.tsv"))
                                         if path.name in {"assay-discordance.tsv", "assay-triage.tsv", "protein-discordance.tsv", "alternative-targets.tsv", "fingerprint-matches.tsv", "candidate-exclusions.tsv", "specificity-queries.tsv", "specificity-queries.excluded.tsv"}}
    result["summary"]["code_sha256"] = digest(__file__)
    write_json(out / "assay-specificity.summary.json", result["summary"])
    write_json(out / "assay-specificity.data.json", result)
    return result["summary"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cross-dir", type=Path, required=True)
    parser.add_argument("--protein-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--primary", choices=[*DATASETS, "all"], default="all")
    parser.add_argument("--match-span", type=int, default=1000000, help="Maximum bounded locus span/profile match (bp); not an LD block")
    parser.add_argument("--verify-source-files", action="store_true", help="Verify all .sig.txt hashes against cross-check input receipts")
    args = parser.parse_args()
    if args.match_span < 0:
        parser.error("match-span must be nonnegative")
    print(json.dumps(run(args.cross_dir, args.protein_map, args.out, args.primary, args.match_span, args.verify_source_files), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
