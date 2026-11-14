#!/usr/bin/env python3
"""Build exhaustive, probe-preserving assay-target fingerprint dossiers.

This is a selected-locus hypothesis screen, not a test of captured identity.
The nominating cis anchor is removed before describing additional support;
physical windows, 5 Mb groups and chromosomes are not independent LD units.
No directions, colocalization probabilities or off-target probabilities are
inferred from selected summary statistics.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import assay_specificity as spec

SCHEMA = "phole-fingerprint-dossiers-v1"
DATASETS = spec.DATASETS
PLATFORMS = spec.PLATFORMS
TIERS = ("T1_distinctive_heldout", "T2_specific_full_profile", "T3_ambiguous_multilocus",
         "T4_specific_anchor_only", "T5_regional_or_hotspot_context")
PRIOR_CLASSES = {"specific_anchor_and_fingerprint", "specific_fingerprint_with_regional_caution"}
MODES = ("FULL", "REGION_HELDOUT", "BLOCK5MB_HELDOUT", "CHROMOSOME_HELDOUT")
BLOCK = 5_000_000
METRICS = ("PRIMARY_LOCI", "ALTERNATIVE_LOCI", "SHARED_REGIONS", "EXACT_RSID", "JACCARD",
           "NONHOTSPOT_PRIMARY", "NONHOTSPOT_BOTH", "NONHOTSPOT_BOTH_COMMON_UNIVERSE",
           "INTENDED_MAX_SHARED", "INTENDED_MAX_JACCARD", "INTENDED_BEST_ASSAY", "INTENDED_PROFILE_STATUS",
           "PREFERRED_TO_ALL_INTENDED")


def uid(prefix, values):
    return prefix + hashlib.sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()[:20]


def flagset(value):
    return spec.split_flags(value)


def joined(values):
    return ";".join(sorted(set(str(v) for v in values))) or "."


def complete_hypothesis_id(row):
    return uid("ALT_", (row["PROT"], row["PRIMARY"], row["PRIMARY_ASSAY"], row["ALTERNATIVE_PROT"],
                         row["REFERENCE"], row["REFERENCE_ASSAY"]))


def load_specificity(cross_dir, protein_map, specificity_dir, primary):
    payload, rows, cross_receipt = spec.load_bundle(cross_dir, protein_map)
    directory = Path(specificity_dir)
    summary_file = directory / "assay-specificity.summary.json"
    data_file = directory / "assay-specificity.data.json"
    if not summary_file.is_file() or not data_file.is_file():
        raise ValueError("Incomplete specificity bundle; run cross-check specificity first")
    summary = json.loads(summary_file.read_text())
    if summary.get("schema") != spec.SCHEMA:
        raise ValueError("Unsupported specificity schema")
    if summary.get("primary") not in {"all", primary} or (primary == "all" and summary.get("primary") != "all"):
        raise ValueError("Specificity bundle does not cover the requested primary datasets")
    for name, expected in cross_receipt["files"].items():
        if summary.get("input_provenance", {}).get("files", {}).get(name) != expected:
            raise ValueError("Specificity input hashes do not match this cross-check/map: " + name)
    if summary.get("code_sha256") != spec.digest(spec.__file__):
        raise ValueError("Specificity producer changed; rerun specificity before fingerprint dossiers")
    required = {"alternative-targets.tsv", "assay-discordance.tsv", "assay-triage.tsv", "fingerprint-matches.tsv"}
    if not required.issubset(summary.get("output_sha256", {})):
        raise ValueError("Specificity output receipts are incomplete")
    for name, expected in summary["output_sha256"].items():
        path = directory / name
        if Path(name).name != name or not path.is_file() or spec.digest(path) != expected:
            raise ValueError("Stale specificity output: " + name)
    data = json.loads(data_file.read_text())
    if data.get("summary") != summary:
        raise ValueError("Specificity data and summary disagree")
    # TSVs are the checksum-protected analytical source. The data JSON must
    # describe identical alternative rows, not merely carry the same receipt.
    alternatives = spec.read_tsv(directory / "alternative-targets.tsv")
    rendered = [{k: "." if r.get(k) is None else str(r.get(k, ".")) for k in spec.ALTERNATIVE_FIELDS}
                for r in data.get("alternatives", [])]
    if rendered != alternatives:
        raise ValueError("Specificity JSON alternative rows differ from protected TSV")
    inputs = {k: v for k, v in cross_receipt["files"].items()}
    inputs["specificity/assay-specificity.data.json"] = spec.digest(data_file)
    inputs["specificity/assay-specificity.summary.json"] = spec.digest(summary_file)
    inputs.update({"specificity/" + k: v for k, v in summary["output_sha256"].items()})
    return payload, rows, alternatives, summary, inputs


def normalize_context(payload, signal_rows):
    proteins = payload["proteins"]
    lookup = {p["prot"]: p for p in proteins}
    if len(lookup) != len(proteins):
        raise ValueError("Duplicate canonical protein")
    thresholds = {d: spec.pvalue(payload["summary"]["strong_thresholds"][d]) for d in DATASETS}
    if any(not 0 < value <= spec.GWAS_THRESHOLD for value in thresholds.values()):
        raise ValueError("Invalid cross-check strict threshold")
    signals, by_assay, assay_index = spec.normalize_signals(signal_rows, proteins)
    profiles = defaultdict(list)
    strict_cis = defaultdict(list)
    for s in signals:
        if s["clean"] and s["kind"] in {"cis", "trans"} and s["p_decimal"] <= thresholds[s["dataset"]]:
            profiles[(s["prot"], s["dataset"], s["assay"])].append(s)
            if s["kind"] == "cis":
                strict_cis[(s["prot"], s["dataset"])].append(s)
    return lookup, thresholds, signals, by_assay, assay_index, profiles, strict_cis


def enumerate_extra_seeds(payload, context, selected, span, cis_states, same_cohort=False):
    """Enumerate declared extra populations without relabelling primary cis.

    The weak-cis extension uses reported GWAS cis below 5e-8 but above the
    strict threshold. A separate same-cohort screen uses the same no/weak-cis
    population and allows an alternative cis from the primary cohort. In both
    cases clean intended cis on an orthogonal platform remains required.
    Authoritative v3 inputs/results are never mutated.
    """
    lookup, thresholds, signals, by_assay, assay_index, profiles, strict_cis = context
    flank = int(payload["config"]["cis_flank"])
    intervals = defaultdict(list)
    for p in lookup.values():
        if not spec.identity_verified(p):
            continue
        for t in p.get("intervals", []):
            if t.get("biotype") in spec.CODING:
                intervals[(str(t.get("build", "38")), spec.chrom(t["chr"]))].append(
                    (max(1, int(t["start"]) - flank), int(t["end"]) + flank, p["prot"]))
    seeds = []
    eligible_assays = set()
    for key, assay in sorted(assay_index.items()):
        name, dataset, alias = key
        p = lookup[name]
        if dataset not in selected or spec.assay_state(p, dataset, assay, thresholds[dataset]) not in cis_states:
            continue
        donors = [d for d in DATASETS if PLATFORMS[d] != PLATFORMS[dataset] and strict_cis.get((name, d))]
        if not donors:
            continue
        eligible_assays.add(key)
        for s in profiles.get(key, []):
            if s["kind"] != "trans":
                continue
            targets = sorted({target for start, end, target in intervals.get((s["build"], s["chr"]), [])
                              if start <= s["pos"] <= end and target != name})
            for target in targets:
                if spec.same_coding_gene(p, lookup[target]):
                    continue
                for reference in ((dataset,) if same_cohort else DATASETS):
                    if not same_cohort and PLATFORMS[reference] == PLATFORMS[dataset]:
                        continue
                    hits = defaultdict(list)
                    for t in strict_cis.get((target, reference), []):
                        if t["chr"] == s["chr"] and t["build"] == s["build"] and abs(t["pos"] - s["pos"]) <= span:
                            hits[t["assay"]].append(t)
                    for reference_assay, values in sorted(hits.items()):
                        t = min(values, key=lambda t: (spec.variant_match(s, t) != "same_rsid_and_position",
                                                       s["pos"] != t["pos"], abs(s["pos"] - t["pos"]), t["p_decimal"], t["snp"]))
                        row = dict(PROT=name, PRIMARY=dataset, PRIMARY_ASSAY=alias, ALTERNATIVE_PROT=target,
                                   REFERENCE=reference, REFERENCE_ASSAY=reference_assay, BUILD=s["build"], CHR=s["chr"],
                                   PRIMARY_POS=s["pos"], PRIMARY_SNP=s["snp"], PRIMARY_P=s["p"],
                                   ALTERNATIVE_CIS_POS=t["pos"], ALTERNATIVE_CIS_SNP=t["snp"], ALTERNATIVE_CIS_P=t["p"],
                                   LEAD_MATCH=spec.variant_match(s, t), LEAD_DISTANCE_BP=abs(s["pos"] - t["pos"]),
                                   FOLLOW_UP_CLASS="same_cohort_context" if same_cohort else "weak_cis_extension", FLAGS=".",
                                   SCREEN_STRATUM="weak_cis" if spec.assay_state(p, dataset, assay, thresholds[dataset]) == "cis_gwas_only" else "no_reported_cis")
                        row["HYPOTHESIS_ID"] = complete_hypothesis_id(row)
                        seeds.append(row)
    return seeds, eligible_assays


def merged_windows(windows):
    result = []
    for build, chrom, start, end in sorted(set(windows), key=lambda x: (x[0], spec.chrom_key(x[1]), x[2], x[3])):
        if result and result[-1][:2] == (build, chrom) and start <= result[-1][3] + 1:
            result[-1] = (build, chrom, result[-1][2], max(result[-1][3], end))
        else:
            result.append((build, chrom, start, end))
    return result


def holdout_windows(protein, anchors, flank, span):
    windows = [(str(t.get("build", "38")), spec.chrom(t["chr"]), max(1, int(t["start"]) - flank), int(t["end"]) + flank)
               for t in protein.get("intervals", [])]
    for a in anchors:
        low = min(int(a["PRIMARY_POS"]), int(a["ALTERNATIVE_CIS_POS"]))
        high = max(int(a["PRIMARY_POS"]), int(a["ALTERNATIVE_CIS_POS"]))
        windows.append((a["BUILD"], spec.chrom(a["CHR"]), max(1, low - span), high + span))
    region = merged_windows(windows)
    blocks = list(region)
    for a in anchors:
        for pos in (int(a["PRIMARY_POS"]), int(a["ALTERNATIVE_CIS_POS"])):
            start = ((pos - 1) // BLOCK) * BLOCK + 1
            blocks.append((a["BUILD"], spec.chrom(a["CHR"]), start, start + BLOCK - 1))
    return region, merged_windows(blocks), {(b, c) for b, c, _, _ in region}


def windows_text(windows):
    return ";".join(f"{b}:{c}:{start}-{end}" for b, c, start, end in windows) or "."


def filtered_profile(profile, span, windows=(), chromosomes=()):
    """Drop complete original clusters intersecting the exclusion, before rematching."""
    keep = []
    for locus in spec.profile_loci(profile, span):
        if (locus["build"], locus["chr"]) in chromosomes:
            continue
        if any(b == locus["build"] and c == locus["chr"] and start <= locus["end"] and end >= locus["start"]
               for b, c, start, end in windows):
            continue
        keep.extend(locus["signals"])
    return keep


def annotated_match(left, right, span, counts, left_dataset, right_dataset, common_counts, denominator):
    result = spec.match_profiles(left, right, span, counts)
    matched = []
    for hit in result["matches"]:
        l = dict(dataset=left_dataset, build=hit["BUILD"], chr=hit["CHR"], pos=hit["PRIMARY_POS"], snp=hit["PRIMARY_SNP"])
        r = dict(dataset=right_dataset, build=hit["BUILD"], chr=hit["CHR"], pos=hit["REFERENCE_POS"], snp=hit["REFERENCE_SNP"])
        le, lr = counts(l); re, rr = counts(r)
        cle, clr = common_counts(l); cre, crr = common_counts(r)
        matched.append({**hit, "PRIMARY_EXACT_TARGETS": le, "PRIMARY_REGION_TRANS_TARGETS": lr,
                        "REFERENCE_EXACT_TARGETS": re, "REFERENCE_REGION_TRANS_TARGETS": rr,
                        "NONHOTSPOT_PRIMARY": int(le <= 10 and lr <= 10),
                        "NONHOTSPOT_REFERENCE": int(re <= 10 and rr <= 10),
                        "NONHOTSPOT_BOTH": int(max(le, lr, re, rr) <= 10),
                        "COMMON_UNIVERSE_TARGETS": denominator,
                        "COMMON_PRIMARY_EXACT_TARGETS": cle, "COMMON_PRIMARY_REGION_TRANS_TARGETS": clr,
                        "COMMON_REFERENCE_EXACT_TARGETS": cre, "COMMON_REFERENCE_REGION_TRANS_TARGETS": crr,
                        "COMMON_PRIMARY_EXACT_FRACTION": cle / denominator if denominator else None,
                        "COMMON_REFERENCE_EXACT_FRACTION": cre / denominator if denominator else None,
                        "COMMON_PRIMARY_REGION_FRACTION": clr / denominator if denominator else None,
                        "COMMON_REFERENCE_REGION_FRACTION": crr / denominator if denominator else None,
                        "NONHOTSPOT_BOTH_COMMON_UNIVERSE": int(max(cle, clr, cre, crr) <= 10) if denominator else None})
    result["matches"] = matched
    result["nonhotspot_primary"] = sum(r["NONHOTSPOT_PRIMARY"] for r in matched)
    result["nonhotspot_both"] = sum(r["NONHOTSPOT_BOTH"] for r in matched)
    result["nonhotspot_both_common"] = sum(r["NONHOTSPOT_BOTH_COMMON_UNIVERSE"] == 1 for r in matched) if denominator else None
    return result


def intended_comparison(alternative, intended, has_full_profile):
    """Preference must beat every usable intended probe in both descriptive metrics.

    No measured full strict profile means unknown, not a zero comparator. An
    observed profile entirely removed by the declared holdout can have zero
    retained loci, but that conditioned comparison remains descriptive.
    """
    if not intended:
        return dict(max_shared=None, max_jaccard=None, best_assay=".", status="not_measured_or_unusable", preferred=0)
    best = max(intended, key=lambda item: (item[1]["shared"], item[1]["jaccard"], item[1]["exact"], item[0]))
    status = "observed" if has_full_profile else "no_clean_strict_loci"
    max_shared = max(value["shared"] for _, value in intended)
    max_jaccard = max(value["jaccard"] for _, value in intended)
    preferred = bool(has_full_profile and alternative["shared"] > max_shared and alternative["jaccard"] > max_jaccard)
    return dict(max_shared=max_shared, max_jaccard=max_jaccard, best_assay=best[0], status=status, preferred=int(preferred))


def apply_metrics(row, prefix, comparison, intended):
    row.update({prefix + "_PRIMARY_LOCI": comparison["left_n"], prefix + "_ALTERNATIVE_LOCI": comparison["right_n"],
                prefix + "_SHARED_REGIONS": comparison["shared"], prefix + "_EXACT_RSID": comparison["exact"],
                prefix + "_JACCARD": comparison["jaccard"], prefix + "_NONHOTSPOT_PRIMARY": comparison["nonhotspot_primary"],
                prefix + "_NONHOTSPOT_BOTH": comparison["nonhotspot_both"],
                prefix + "_NONHOTSPOT_BOTH_COMMON_UNIVERSE": comparison["nonhotspot_both_common"],
                prefix + "_INTENDED_MAX_SHARED": intended["max_shared"], prefix + "_INTENDED_MAX_JACCARD": intended["max_jaccard"],
                prefix + "_INTENDED_BEST_ASSAY": intended["best_assay"], prefix + "_INTENDED_PROFILE_STATUS": intended["status"],
                prefix + "_PREFERRED_TO_ALL_INTENDED": intended["preferred"]})


def analyse(payload, signal_rows, seed_rows, primary="all", span=1_000_000, include_weak_cis=True):
    if primary not in {*DATASETS, "all"} or span < 0:
        raise ValueError("Invalid primary or match span")
    selected = DATASETS if primary == "all" else (primary,)
    context = normalize_context(payload, signal_rows)
    lookup, thresholds, signals, by_assay, assay_index, profiles, strict_cis = context
    counts = spec.make_hotspot_index(signals)
    common_cache = {}

    def common(left, right):
        key = tuple(sorted((left, right)))
        if key not in common_cache:
            names = {name for name, p in lookup.items() if all(any(spec.summary_valid(a) for a in p["data"][d].get("assays", [])) for d in key)}
            common_cache[key] = (spec.make_hotspot_index([s for s in signals if s["prot"] in names]), len(names))
        return common_cache[key]

    main_seeds = [dict(r, SCREEN_STRATUM="no_reported_cis") for r in seed_rows if r["PRIMARY"] in selected]
    weak_seeds, weak_assays = enumerate_extra_seeds(payload, context, selected, span, {"cis_gwas_only"}) if include_weak_cis else ([], set())
    seeds = main_seeds + [dict(r, SCREEN_STRATUM="weak_cis") for r in weak_seeds]
    seed_keys = set(); deduplicated = []
    for r in seeds:
        if r["HYPOTHESIS_ID"] != complete_hypothesis_id(r):
            raise ValueError("Alternative hypothesis ID does not match its assay identities")
        for role, dataset, protein, alias, pos, snp, p, kind in (
            ("primary", r["PRIMARY"], r["PROT"], r["PRIMARY_ASSAY"], r["PRIMARY_POS"], r["PRIMARY_SNP"], r["PRIMARY_P"], "trans"),
            ("alternative", r["REFERENCE"], r["ALTERNATIVE_PROT"], r["REFERENCE_ASSAY"], r["ALTERNATIVE_CIS_POS"], r["ALTERNATIVE_CIS_SNP"], r["ALTERNATIVE_CIS_P"], "cis")):
            if not any(s["build"] == str(r["BUILD"]) and s["chr"] == spec.chrom(r["CHR"]) and s["pos"] == int(pos)
                       and s["snp"] == snp and s["p_decimal"] == spec.pvalue(p) and s["kind"] == kind
                       for s in profiles.get((protein, dataset, alias), [])):
                raise ValueError(f"Alternative {role} anchor is not a clean strict signal: {r['HYPOTHESIS_ID']}")
        key = (r["HYPOTHESIS_ID"], str(r["BUILD"]), spec.chrom(r["CHR"]), int(r["PRIMARY_POS"]), r["PRIMARY_SNP"],
               int(r["ALTERNATIVE_CIS_POS"]), r["ALTERNATIVE_CIS_SNP"])
        if key in seed_keys:
            continue
        seed_keys.add(key); deduplicated.append(r)
    seeds = deduplicated
    anchor_alternatives = defaultdict(set); anchor_exact_alternatives = defaultdict(set)
    for r in seeds:
        key = (r["PRIMARY"], r["PROT"], r["PRIMARY_ASSAY"], str(r["BUILD"]), spec.chrom(r["CHR"]), int(r["PRIMARY_POS"]), r["PRIMARY_SNP"])
        anchor_alternatives[key].add(r["ALTERNATIVE_PROT"])
        if r["LEAD_MATCH"] == "same_rsid_and_position":
            anchor_exact_alternatives[key].add(r["ALTERNATIVE_PROT"])
    groups = defaultdict(list)
    for r in seeds:
        groups[r["HYPOTHESIS_ID"]].append(r)
    hypotheses = []; anchors_out = []; matches_out = []; probe_out = []; exclusions = []
    hypothesis_windows = {}; hypothesis_modes = {}; hypothesis_keys = {}
    flank = int(payload["config"]["cis_flank"])
    for hid, anchors in sorted(groups.items()):
        seed = anchors[0]
        name, dataset, alias, target, reference, donor = (seed[k] for k in ("PROT", "PRIMARY", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "REFERENCE", "REFERENCE_ASSAY"))
        p, q = lookup[name], lookup[target]
        primary_key = (name, dataset, alias); alternate_key = (target, reference, donor)
        a = assay_index[primary_key]
        state = spec.assay_state(p, dataset, a, thresholds[dataset])
        expected_state = spec.NO_CIS_STATES if seed["SCREEN_STRATUM"] == "no_reported_cis" else {"cis_gwas_only"}
        if state not in expected_state or not spec.identity_verified(q) or spec.same_coding_gene(p, q):
            raise ValueError("Candidate is outside its declared primary cis/target identity scope")
        pair_id = uid("PAIR_", (dataset, name, target))
        own_keys = [(name, reference, x["alias"]) for x in p["data"][reference].get("assays", []) if spec.summary_valid(x)]
        own_cis = strict_cis.get((name, reference), [])
        best_cis = min(own_cis, key=lambda s: (s["p_decimal"], s["assay"], s["snp"])) if own_cis else None
        common_counts, denominator = common(dataset, reference)
        specific_anchors = []; all_cautions = set()
        for anchor in anchors:
            key = (dataset, name, alias, str(anchor["BUILD"]), spec.chrom(anchor["CHR"]), int(anchor["PRIMARY_POS"]), anchor["PRIMARY_SNP"])
            ls = dict(dataset=dataset, build=key[3], chr=key[4], pos=key[5], snp=key[6])
            rs = dict(dataset=reference, build=key[3], chr=key[4], pos=int(anchor["ALTERNATIVE_CIS_POS"]), snp=anchor["ALTERNATIVE_CIS_SNP"])
            le, lr = counts(ls); re, rr = counts(rs)
            my_target, other_target = p.get("target") or {}, q.get("target") or {}
            nearby = (my_target.get("chr") is not None and my_target.get("chr") == other_target.get("chr") and
                      abs((my_target.get("start", 0) + my_target.get("end", 0)) / 2 - (other_target.get("start", 0) + other_target.get("end", 0)) / 2) <= span)
            shared_coordinate = q.get("coordinate", {}).get("shared_anchor_targets", ".") not in (None, "", ".")
            reproduced = [d for d in DATASETS if PLATFORMS[d] != PLATFORMS[dataset]
                          and any(s["kind"] == "trans" and s["build"] == key[3] and s["chr"] == key[4] and abs(s["pos"] - key[5]) <= span
                                  for x in p["data"][d].get("assays", []) for s in profiles.get((name, d, x["alias"]), []))]
            specific = bool(anchor["LEAD_MATCH"] == "same_rsid_and_position" and le <= 10 and
                            len(anchor_exact_alternatives[key]) == 1 and not nearby and not shared_coordinate and not reproduced)
            caution = set()
            if le > 10: caution.add("primary_exact_variant_hotspot")
            if lr > 10: caution.add("primary_regional_hotspot")
            if re > 10: caution.add("reference_exact_variant_hotspot")
            if rr > 10: caution.add("reference_regional_hotspot")
            if len(anchor_alternatives[key]) > 1: caution.add("multiple_regional_alternative_targets")
            if len(anchor_exact_alternatives[key]) > 1: caution.add("multiple_exact_alternative_targets")
            if nearby: caution.add("neighboring_coding_targets")
            if shared_coordinate: caution.add("alternative_coordinate_shared")
            if reproduced: caution.add("intended_target_anchor_trans_also_reported")
            all_cautions.update(caution)
            item = dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, PROT=name, PRIMARY=dataset, PRIMARY_ASSAY=alias,
                        ALTERNATIVE_PROT=target, REFERENCE=reference, REFERENCE_ASSAY=donor,
                        SCREEN_STRATUM=seed["SCREEN_STRATUM"], BUILD=key[3], CHR=key[4], PRIMARY_POS=key[5],
                        PRIMARY_SNP=key[6], PRIMARY_P=str(anchor["PRIMARY_P"]), ALTERNATIVE_CIS_POS=int(anchor["ALTERNATIVE_CIS_POS"]),
                        ALTERNATIVE_CIS_SNP=anchor["ALTERNATIVE_CIS_SNP"], ALTERNATIVE_CIS_P=str(anchor["ALTERNATIVE_CIS_P"]),
                        LEAD_MATCH=anchor["LEAD_MATCH"], LEAD_DISTANCE_BP=abs(key[5] - int(anchor["ALTERNATIVE_CIS_POS"])),
                        PRIMARY_EXACT_TARGETS=le, PRIMARY_REGION_TRANS_TARGETS=lr,
                        REFERENCE_EXACT_TARGETS=re, REFERENCE_REGION_TRANS_TARGETS=rr,
                        ALTERNATIVE_TARGETS_AT_LOCUS=len(anchor_alternatives[key]), EXACT_ALTERNATIVE_TARGETS_AT_VARIANT=len(anchor_exact_alternatives[key]),
                        INTENDED_TARGET_TRANS_COHORTS=joined(reproduced), ANCHOR_SPECIFIC=int(specific), CAUTIONS=joined(caution))
            anchors_out.append(item); specific_anchors.append(item)
        best_anchor = min(specific_anchors, key=lambda x: (-x["ANCHOR_SPECIFIC"], x["LEAD_MATCH"] != "same_rsid_and_position",
                                                         x["LEAD_DISTANCE_BP"], x["PRIMARY_EXACT_TARGETS"], spec.pvalue(x["PRIMARY_P"]),
                                                         spec.chrom_key(x["CHR"]), x["PRIMARY_POS"], x["PRIMARY_SNP"]))
        region, block, chromosomes = holdout_windows(q, anchors, flank, span)
        hypothesis_windows[hid] = (region, block, chromosomes)
        hypothesis_keys[hid] = (primary_key, alternate_key, own_keys)
        row = dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, PROT=name, PRIMARY=dataset, PRIMARY_PLATFORM=PLATFORMS[dataset], PRIMARY_ASSAY=alias,
                   ALTERNATIVE_PROT=target, REFERENCE=reference, REFERENCE_PLATFORM=PLATFORMS[reference], REFERENCE_ASSAY=donor,
                   SCREEN_STRATUM=seed["SCREEN_STRATUM"], PRIMARY_STATE=state, PRIMARY_CIS_P=(a.get("best_cis") or {}).get("p"),
                   PRIMARY_CIS_SEARCH_COMPLETE=int(a.get("cis_search_complete") is True),
                   PRIMARY_ASSAYS_EXPECTED=len(p["data"][dataset].get("assays", [])),
                   PRIMARY_ASSAYS_WITH_STRICT_CIS=len({s["assay"] for s in strict_cis.get((name, dataset), [])}),
                   REFERENCE_INTENDED_ASSAYS_EXPECTED=len(p["data"][reference].get("assays", [])), REFERENCE_INTENDED_ASSAYS_USABLE=len(own_keys),
                   REFERENCE_INTENDED_ASSAYS_WITH_STRICT_CIS=len({s["assay"] for s in own_cis}), SAME_REFERENCE_INTENDED_STRICT_CIS=int(bool(own_cis)),
                   REFERENCE_INTENDED_CIS_P=best_cis["p"] if best_cis else None,
                   REFERENCE_INTENDED_CIS_SNP=best_cis["snp"] if best_cis else None,
                   REFERENCE_INTENDED_CIS_ASSAY=best_cis["assay"] if best_cis else None,
                   INTENDED_STRICT_CIS_COHORTS=joined(d for d in DATASETS if strict_cis.get((name, d))),
                   ALTERNATIVE_STRICT_CIS_COHORTS=joined(d for d in DATASETS if strict_cis.get((target, d))),
                   N_ANCHOR_NOMINATIONS=len(anchors), ANCHOR_SPECIFIC=int(any(x["ANCHOR_SPECIFIC"] for x in specific_anchors)),
                   EXACT_ANCHOR_MULTILOCUS=0,
                   ANCHOR_SNP=best_anchor["PRIMARY_SNP"], ANCHOR_CHR=best_anchor["CHR"], ANCHOR_POS=best_anchor["PRIMARY_POS"],
                   ANCHOR_PRIMARY_P=best_anchor["PRIMARY_P"], ANCHOR_REFERENCE_SNP=best_anchor["ALTERNATIVE_CIS_SNP"],
                   ANCHOR_REFERENCE_P=best_anchor["ALTERNATIVE_CIS_P"], ANCHOR_LEAD_MATCH=best_anchor["LEAD_MATCH"],
                   ANCHOR_PRIMARY_EXACT_TARGETS=best_anchor["PRIMARY_EXACT_TARGETS"], ANCHOR_REFERENCE_EXACT_TARGETS=best_anchor["REFERENCE_EXACT_TARGETS"],
                   ANCHOR_PRIMARY_REGION_TARGETS=best_anchor["PRIMARY_REGION_TRANS_TARGETS"], ANCHOR_REFERENCE_REGION_TARGETS=best_anchor["REFERENCE_REGION_TRANS_TARGETS"],
                   COMMON_UNIVERSE_TARGETS=denominator, REGION_EXCLUDED=windows_text(region), BLOCK5MB_EXCLUDED=windows_text(block),
                   CHROMOSOMES_EXCLUDED=joined(f"{b}:{c}" for b, c in chromosomes),
                   PRIOR_PRIORITY_CASE=int(any(x.get("FOLLOW_UP_CLASS") in PRIOR_CLASSES for x in anchors)),
                   DIRECTION_TESTED=0, COLOCALIZATION_TESTED=0, WRONG_PROTEIN_CONFIRMED=0)
        modes = {}
        for mode in MODES:
            def restricted(key):
                if mode == "FULL": return profiles.get(key, [])
                if mode == "CHROMOSOME_HELDOUT": return filtered_profile(profiles.get(key, []), span, chromosomes=chromosomes)
                return filtered_profile(profiles.get(key, []), span, windows=region if mode == "REGION_HELDOUT" else block)
            left = restricted(primary_key)
            result = annotated_match(left, restricted(alternate_key), span, counts, dataset, reference, common_counts, denominator)
            intended = [(key[2], annotated_match(left, restricted(key), span, counts, dataset, reference, common_counts, denominator)) for key in own_keys]
            own = intended_comparison(result, intended, any(profiles.get(key) for key in own_keys))
            apply_metrics(row, mode, result, own); modes[mode] = (result, intended)
            for role, target_name, reference_alias, comparison in [("alternative", target, donor, result)] + [("intended", name, name_alias, value) for name_alias, value in intended]:
                probe_out.append(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, MODE=mode, COMPARISON=role, PROT=name, PRIMARY=dataset,
                                      PRIMARY_ASSAY=alias, REFERENCE_PROT=target_name, REFERENCE=reference, REFERENCE_ASSAY=reference_alias,
                                      PRIMARY_LOCI=comparison["left_n"], REFERENCE_LOCI=comparison["right_n"], SHARED_REGIONS=comparison["shared"],
                                      EXACT_RSID=comparison["exact"], JACCARD=comparison["jaccard"], NONHOTSPOT_PRIMARY=comparison["nonhotspot_primary"],
                                      NONHOTSPOT_BOTH=comparison["nonhotspot_both"], NONHOTSPOT_BOTH_COMMON_UNIVERSE=comparison["nonhotspot_both_common"]))
                matches_out.extend(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, MODE=mode, COMPARISON=role, PROT=name, PRIMARY=dataset,
                                        PRIMARY_ASSAY=alias, REFERENCE_PROT=target_name, REFERENCE=reference, REFERENCE_ASSAY=reference_alias, **hit)
                                   for hit in comparison["matches"])
        hypothesis_modes[hid] = modes
        full_specific = bool(row["ANCHOR_SPECIFIC"] and row["SAME_REFERENCE_INTENDED_STRICT_CIS"] and row["FULL_SHARED_REGIONS"] >= 2 and row["FULL_NONHOTSPOT_PRIMARY"] >= 1 and row["FULL_PREFERRED_TO_ALL_INTENDED"])
        distinct = bool(full_specific and row["REGION_HELDOUT_SHARED_REGIONS"] >= 1 and row["REGION_HELDOUT_NONHOTSPOT_BOTH"] >= 1 and row["REGION_HELDOUT_PREFERRED_TO_ALL_INTENDED"])
        row["EXACT_ANCHOR_MULTILOCUS"] = int(row["FULL_SHARED_REGIONS"] >= 2 and any(x["LEAD_MATCH"] == "same_rsid_and_position" for x in specific_anchors))
        row["TIER"] = (TIERS[0] if distinct else TIERS[1] if full_specific else TIERS[2] if row["FULL_SHARED_REGIONS"] >= 2 and row["FULL_PREFERRED_TO_ALL_INTENDED"]
                       else TIERS[3] if row["ANCHOR_SPECIFIC"] else TIERS[4])
        reasons = set()
        if not row["ANCHOR_SPECIFIC"]: reasons.add("anchor_not_specific")
        if not row["SAME_REFERENCE_INTENDED_STRICT_CIS"]: reasons.add("same_reference_intended_strict_cis_unavailable")
        if row["FULL_SHARED_REGIONS"] < 2: reasons.add("fewer_than_two_full_profile_regions")
        if not row["FULL_NONHOTSPOT_PRIMARY"]: reasons.add("no_primary_nonhotspot_full_profile_region")
        if not row["FULL_PREFERRED_TO_ALL_INTENDED"]: reasons.add("full_profile_not_better_than_all_intended_probes")
        if not row["REGION_HELDOUT_SHARED_REGIONS"]: reasons.add("no_matching_region_after_anchor_exclusion")
        if not row["REGION_HELDOUT_NONHOTSPOT_BOTH"]: reasons.add("no_both_platform_nonhotspot_secondary_support")
        if not row["REGION_HELDOUT_PREFERRED_TO_ALL_INTENDED"]: reasons.add("anchor_heldout_profile_not_better_than_all_intended_probes")
        if not row["PRIMARY_CIS_SEARCH_COMPLETE"]: all_cautions.add("primary_cis_absence_or_strength_not_full_window_certified")
        if row["SCREEN_STRATUM"] == "weak_cis": all_cautions.add("weak_cis_extension_not_absent_cis")
        if not row["CHROMOSOME_HELDOUT_PREFERRED_TO_ALL_INTENDED"]: all_cautions.add("preference_not_retained_after_anchor_chromosome_removal")
        row["STRICT_TIER_NOT_MET_REASONS"] = joined(reasons)
        row["CAUTIONS"] = joined(all_cautions)
        row["INTERPRETATION"] = "selected_locus_assay_target_hypothesis"
        hypotheses.append(row)
        if reasons:
            exclusions.append(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, PROT=name, PRIMARY=dataset, PRIMARY_ASSAY=alias,
                                   ALTERNATIVE_PROT=target, REFERENCE=reference, REFERENCE_ASSAY=donor,
                                   SCREEN_STRATUM=row["SCREEN_STRATUM"], TIER=row["TIER"], REASONS=joined(reasons)))
    hypotheses.sort(key=hypothesis_order)
    pair_groups = defaultdict(list)
    for row in hypotheses: pair_groups[row["PAIR_ID"]].append(row)
    pairs = []
    for pair_id, members in pair_groups.items():
        best = min(members, key=hypothesis_order)
        # The unprefixed profile metrics belong to the explicitly identified
        # representative hypothesis. Counts below aggregate only identities.
        pairs.append({**best, "REPRESENTATIVE_HYPOTHESIS_ID": best["HYPOTHESIS_ID"],
                      "N_ASSAY_HYPOTHESES": len(members), "N_PRIMARY_ASSAYS": len({x["PRIMARY_ASSAY"] for x in members}),
                      "N_REFERENCE_COHORTS": len({x["REFERENCE"] for x in members}),
                      "N_REFERENCE_ASSAYS": len({(x["REFERENCE"], x["REFERENCE_ASSAY"]) for x in members}),
                      "REFERENCE_COHORTS": joined(x["REFERENCE"] for x in members),
                      "SCREEN_STRATA": joined(x["SCREEN_STRATUM"] for x in members),
                      "PRIOR_PRIORITY_CASE": int(any(x["PRIOR_PRIORITY_CASE"] for x in members)),
                      "PROFILE_METRIC_SCOPE": "representative_assay_hypothesis_not_pooled_probes"})
    pairs.sort(key=hypothesis_order)
    # Every exact-anchor multilocus case is retained in the casebook even if
    # nonspecificity or an intended comparator prevents a priority tier.
    exact_multi = {r["HYPOTHESIS_ID"] for r in hypotheses if r["FULL_SHARED_REGIONS"] >= 2 and any(a["LEAD_MATCH"] == "same_rsid_and_position" for a in groups[r["HYPOTHESIS_ID"]])}
    display_ids = {r["HYPOTHESIS_ID"] for r in hypotheses if r["TIER"] != TIERS[4] or r["PRIOR_PRIORITY_CASE"]} | exact_multi
    display_pairs = {r["PAIR_ID"] for r in hypotheses if r["HYPOTHESIS_ID"] in display_ids}
    context_rows, context_matches, sensitivity, profile_rows = build_context(
        [r for r in hypotheses if r["HYPOTHESIS_ID"] in display_ids], lookup, thresholds, by_assay, assay_index,
        profiles, counts, common, hypothesis_windows, span)
    same_cohort, same_matches, same_profiles, same_eligible = same_cohort_scan(payload, context, selected, span, include_weak_cis, counts, common)
    counts_summary = dict(same_cohort_eligible_assays=same_eligible, same_cohort_eligible_hypotheses=len(same_cohort),
                          same_cohort_similarity_candidates=sum(r["HIGH_SIMILARITY"] for r in same_cohort),
                          same_cohort_similarity_pairs=len({r["PAIR_ID"] for r in same_cohort if r["HIGH_SIMILARITY"]}),
                          source_no_cis_anchor_rows=len(main_seeds), weak_cis_anchor_rows=len(weak_seeds), weak_cis_eligible_assays=len(weak_assays),
                          anchor_rows=len(seeds), assay_hypotheses=len(hypotheses), directed_primary_protein_pairs=len(pairs),
                          directed_canonical_protein_pairs=len({(r["PROT"], r["ALTERNATIVE_PROT"]) for r in pairs}),
                          prior_priority_hypotheses=sum(r["PRIOR_PRIORITY_CASE"] for r in hypotheses), exact_anchor_multilocus_hypotheses=len(exact_multi),
                          casebook_hypotheses=len(display_ids), casebook_pairs=len(display_pairs), confirmed_wrong_proteins=0,
                          by_tier={tier: sum(r["TIER"] == tier for r in hypotheses) for tier in TIERS},
                          pair_by_tier={tier: sum(r["TIER"] == tier for r in pairs) for tier in TIERS},
                          by_stratum=dict(Counter(r["SCREEN_STRATUM"] for r in hypotheses)),
                          by_primary={d: dict(assay_hypotheses=sum(r["PRIMARY"] == d for r in hypotheses),
                                              directed_pairs=sum(r["PRIMARY"] == d for r in pairs),
                                              prior_priority=sum(r["PRIMARY"] == d and r["PRIOR_PRIORITY_CASE"] for r in hypotheses)) for d in selected})
    summary = dict(schema=SCHEMA, primary=primary, counts=counts_summary, cross_check_generated=payload.get("generated"),
                   rules=dict(main_universe="v3 no-reported-cis assay hypotheses; source alternatives preserved",
                              weak_cis_extension=bool(include_weak_cis), same_reference_intended_strict_cis_required_for_T1_T2=True,
                              match_span_bp=span, full_match="ordered one-to-one bounded-span regions; exact rsID/position prioritized",
                              region_holdout="complete clusters intersecting alternative target cis windows plus every nominating anchor pair +/- match span",
                              block_holdout="region holdout plus 1-based 5 Mb blocks containing nominating primary/reference anchor positions",
                              chromosome_holdout="every chromosome appearing in region holdout; PAR intervals retained",
                              intended_comparator="every usable same-name probe; alternative must exceed both maximal shared-count and maximal Jaccard",
                              hotspot="more than 10 canonical proteins at exact variant OR trans targets within +/-500 kb; separately in both cohorts",
                              common_universe="canonical proteins with at least one usable summary in both compared cohorts; descriptive sensitivity only, never a tier rescue",
                              independent_ld_loci=False, same_platform_context="descriptive assay consistency; never counted as external replication",
                              sensitivity_spans_bp=sorted({100000, 250000, span}),
                              same_cohort_screen="separate complete eligible-population scan: >=4 matched regions and Jaccard>=0.6; descriptive only; does not change cross-platform tiers"),
                   tier_definitions={TIERS[0]: "specific full-profile case, with alternative preference and at least one both-cohort nonhotspot region after anchor exclusion",
                                     TIERS[1]: "exact specific anchor, same-reference intended strict cis and multi-region profile with primary nonhotspot support preferred to every intended probe; T1 holdout condition not met",
                                     TIERS[2]: "multi-region alternative preference but anchor/reference qualification unresolved",
                                     TIERS[3]: "specific exact anchor, insufficient multi-region preference",
                                     TIERS[4]: "regional, hotspot, comparator or target ambiguity; retained as exhaustive context"},
                   interpretation="Hypothesis-generating descriptive screen. No observed wrong-protein confirmations, effect directions, colocalization or capture identities.",
                   casebook_scope="All T1-T4 hypotheses, every prior priority case, and every exact-anchor full-profile multi-region hypothesis; complete pair/assay TSVs remain exhaustive. Separate same_cohort_screen contains every fixed-threshold high-similarity case.")
    return dict(summary=summary, pairs=pairs, hypotheses=hypotheses, anchors=anchors_out, matched_loci=matches_out,
                probe_comparisons=probe_out, exclusions=exclusions, same_platform_context=context_rows, context_matches=context_matches,
                sensitivity=sensitivity, profile_signals=profile_rows, display_ids=display_ids,
                same_cohort_eligible=same_cohort, same_cohort_screen=[r for r in same_cohort if r["HIGH_SIMILARITY"]],
                same_cohort_screen_matches=same_matches, same_cohort_profile_signals=same_profiles)


def hypothesis_order(row):
    return (TIERS.index(row["TIER"]), -int(row.get("PRIOR_PRIORITY_CASE", 0)),
            row["SCREEN_STRATUM"] != "no_reported_cis", -row["REGION_HELDOUT_NONHOTSPOT_BOTH"],
            -row["REGION_HELDOUT_SHARED_REGIONS"], -row["FULL_SHARED_REGIONS"], -row["FULL_JACCARD"],
            row["PRIMARY"], row["PROT"], row["ALTERNATIVE_PROT"], row["REFERENCE"], row["PRIMARY_ASSAY"], row["REFERENCE_ASSAY"])


def build_context(hypotheses, lookup, thresholds, by_assay, assay_index, profiles, counts, common, windows, span):
    contexts = []; matched = []; sensitivity = []; profile_rows = []
    for h in hypotheses:
        hid, pair_id = h["HYPOTHESIS_ID"], h["PAIR_ID"]
        name, dataset, alias, target, reference, donor = (h[k] for k in ("PROT", "PRIMARY", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "REFERENCE", "REFERENCE_ASSAY"))
        primary_key = (name, dataset, alias)
        def keys(protein, cohort):
            return [(protein, cohort, a["alias"]) for a in lookup[protein]["data"][cohort].get("assays", []) if spec.summary_valid(a)]
        own_reference = keys(name, reference)
        alt_primary = keys(target, dataset)
        alt_reference = keys(target, reference)
        comparisons = []
        # Source assay is fixed. Other assays remain separate, never pooled.
        for right in own_reference: comparisons.append(("primary_intended__reference_intended", primary_key, right))
        for right in alt_reference: comparisons.append(("primary_intended__reference_alternative", primary_key, right))
        for left in alt_primary:
            for right in own_reference: comparisons.append(("primary_alternative__reference_intended", left, right))
            for right in alt_reference: comparisons.append(("primary_alternative__reference_alternative", left, right))
        for left in alt_primary: comparisons.append(("same_primary_cohort_alternative", primary_key, left))
        for left in own_reference:
            for right in alt_reference: comparisons.append(("same_reference_cohort_targets", left, right))
        for role, left, right in comparisons:
            cc, denominator = common(left[1], right[1])
            result = annotated_match(profiles.get(left, []), profiles.get(right, []), span, counts, left[1], right[1], cc, denominator)
            cid = uid("CTX_", (hid, role, left, right))
            contexts.append(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, CONTEXT_ID=cid, COMPARISON=role,
                                 LEFT_PROT=left[0], LEFT_DATASET=left[1], LEFT_ASSAY=left[2], RIGHT_PROT=right[0], RIGHT_DATASET=right[1], RIGHT_ASSAY=right[2],
                                 LEFT_LOCI=result["left_n"], RIGHT_LOCI=result["right_n"], SHARED_REGIONS=result["shared"], EXACT_RSID=result["exact"], JACCARD=result["jaccard"],
                                 NONHOTSPOT_BOTH=result["nonhotspot_both"], COMMON_UNIVERSE_TARGETS=denominator,
                                 INDEPENDENT_VALIDATION=0, INTERPRETATION="descriptive_target_by_cohort_consistency"))
            matched.extend(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, CONTEXT_ID=cid, COMPARISON=role,
                                LEFT_PROT=left[0], LEFT_DATASET=left[1], LEFT_ASSAY=left[2], RIGHT_PROT=right[0], RIGHT_DATASET=right[1], RIGHT_ASSAY=right[2], **hit)
                           for hit in result["matches"])
        roles = [("primary_intended", primary_key)] + [("reference_intended", k) for k in own_reference] + [("primary_alternative", k) for k in alt_primary] + [("reference_alternative", k) for k in alt_reference]
        for role, key in roles:
            strict = {(s["chr"], s["pos"], s["snp"], s["kind"]) for s in profiles.get(key, [])}
            for s in sorted(by_assay.get(key, []), key=lambda s: (spec.chrom_key(s["chr"]), s["pos"], s["snp"])):
                profile_rows.append(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, ROLE=role, PROT=key[0], DATASET=key[1], ASSAY=key[2],
                                         BUILD=s["build"], CHR=s["chr"], POS=s["pos"], SNP=s["snp"], P=s["p"], CLASS=s["kind"],
                                         IN_STRICT_PROFILE=int((s["chr"], s["pos"], s["snp"], s["kind"]) in strict),
                                         QC_USABLE=int(s["clean"]), FLAGS=joined(s["flags"])))
        cc, denominator = common(dataset, reference)
        alternate_key = (target, reference, donor)
        region, block, chromosomes = windows[hid]
        # Frozen exclusions from the nomination span prevent tuning the
        # nomination/heldout universe while varying matching tolerance.
        for tolerance in sorted({100000, 250000, span}):
            for mode in MODES:
                def restricted(key):
                    if mode == "FULL": return profiles.get(key, [])
                    if mode == "CHROMOSOME_HELDOUT": return filtered_profile(profiles.get(key, []), span, chromosomes=chromosomes)
                    return filtered_profile(profiles.get(key, []), span, windows=region if mode == "REGION_HELDOUT" else block)
                left = restricted(primary_key)
                result = annotated_match(left, restricted(alternate_key), tolerance, counts, dataset, reference, cc, denominator)
                own = [(key[2], annotated_match(left, restricted(key), tolerance, counts, dataset, reference, cc, denominator)) for key in own_reference]
                better = intended_comparison(result, own, any(profiles.get(key) for key in own_reference))
                row = dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, PROT=name, PRIMARY=dataset, ALTERNATIVE_PROT=target, REFERENCE=reference,
                           MATCH_SPAN_BP=tolerance, MODE=mode, FROZEN_EXCLUSION_SPAN_BP=span)
                apply_metrics(row, "PROFILE", result, better); sensitivity.append(row)
    return contexts, matched, sensitivity, profile_rows


def same_cohort_scan(payload, context, selected, span, include_weak_cis, counts, common):
    """Independent descriptive screen; same-cohort matches never raise a tier.

    Candidate eligibility is the no/weak-primary-cis, clean orthogonal own-cis
    population. Alternatives require their own clean strict cis in that same
    cohort. The fixed >=4 regions and Jaccard>=0.6 rule is not a p-value test.
    """
    lookup, thresholds, signals, by_assay, assay_index, profiles, strict_cis = context
    states = spec.NO_CIS_STATES | ({"cis_gwas_only"} if include_weak_cis else set())
    seeds, eligible = enumerate_extra_seeds(payload, context, selected, span, states, same_cohort=True)
    grouped = defaultdict(list)
    for r in seeds: grouped[r["HYPOTHESIS_ID"]].append(r)
    all_rows = []; all_matches = []; profile_rows = []
    flank = int(payload["config"]["cis_flank"])
    for hid, anchor_rows in sorted(grouped.items()):
        seed = min(anchor_rows, key=lambda x: (x["LEAD_MATCH"] != "same_rsid_and_position", x["LEAD_DISTANCE_BP"], spec.pvalue(x["PRIMARY_P"]), x["PRIMARY_SNP"]))
        name, dataset, alias, target, donor = (seed[k] for k in ("PROT", "PRIMARY", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "REFERENCE_ASSAY"))
        left_key, right_key = (name, dataset, alias), (target, dataset, donor)
        cc, denominator = common(dataset, dataset)
        result = annotated_match(profiles.get(left_key, []), profiles.get(right_key, []), span, counts, dataset, dataset, cc, denominator)
        primary_assay = assay_index[left_key]
        reference_results = {"intended": [], "alternative": []}
        for role, protein in (("intended", name), ("alternative", target)):
            for cohort in DATASETS:
                if PLATFORMS[cohort] == PLATFORMS[dataset]: continue
                paired, count = common(dataset, cohort)
                for assay in lookup[protein]["data"][cohort].get("assays", []):
                    key = (protein, cohort, assay["alias"])
                    if not spec.summary_valid(assay) or not profiles.get(key): continue
                    comparison = annotated_match(profiles.get(left_key, []), profiles[key], span, counts, dataset, cohort, paired, count)
                    reference_results[role].append((cohort, assay["alias"], comparison))
        best = {role: max(values, key=lambda v: (v[2]["shared"], v[2]["jaccard"], v[2]["exact"], v[0], v[1])) if values else None
                for role, values in reference_results.items()}
        own_cis = [sig for cohort in DATASETS if PLATFORMS[cohort] != PLATFORMS[dataset] for sig in strict_cis.get((name, cohort), [])]
        own_best = min(own_cis, key=lambda x: (x["p_decimal"], x["dataset"], x["assay"], x["snp"])) if own_cis else None
        region, block, chromosomes = holdout_windows(lookup[target], anchor_rows, flank, span)
        held = annotated_match(filtered_profile(profiles.get(left_key, []), span, windows=region),
                               filtered_profile(profiles.get(right_key, []), span, windows=region), span, counts, dataset, dataset, cc, denominator)
        pair_id = uid("PAIR_", (dataset, name, target))
        row = dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, PROT=name, PRIMARY=dataset, PRIMARY_PLATFORM=PLATFORMS[dataset],
                   PRIMARY_ASSAY=alias, ALTERNATIVE_PROT=target, ALTERNATIVE_ASSAY=donor, SCREEN_STRATUM=seed["SCREEN_STRATUM"],
                   PRIMARY_STATE=spec.assay_state(lookup[name], dataset, primary_assay, thresholds[dataset]),
                   PRIMARY_CIS_P=(primary_assay.get("best_cis") or {}).get("p"), PRIMARY_CIS_SEARCH_COMPLETE=int(primary_assay.get("cis_search_complete") is True),
                   ANCHOR_CHR=seed["CHR"], ANCHOR_POS=seed["PRIMARY_POS"], ANCHOR_SNP=seed["PRIMARY_SNP"], ANCHOR_PRIMARY_P=seed["PRIMARY_P"],
                   ALTERNATIVE_CIS_SNP=seed["ALTERNATIVE_CIS_SNP"], ALTERNATIVE_CIS_P=seed["ALTERNATIVE_CIS_P"], ANCHOR_MATCH=seed["LEAD_MATCH"],
                   N_ANCHOR_NOMINATIONS=len(anchor_rows), PRIMARY_LOCI=result["left_n"], ALTERNATIVE_LOCI=result["right_n"],
                   SHARED_REGIONS=result["shared"], EXACT_RSID=result["exact"], JACCARD=result["jaccard"], NONHOTSPOT_BOTH=result["nonhotspot_both"],
                   REGION_HELDOUT_SHARED_REGIONS=held["shared"], REGION_HELDOUT_EXACT_RSID=held["exact"], REGION_HELDOUT_JACCARD=held["jaccard"],
                   REGION_HELDOUT_NONHOTSPOT_BOTH=held["nonhotspot_both"], REGION_EXCLUDED=windows_text(region),
                   ORTHOGONAL_INTENDED_CIS_COHORTS=joined(x["dataset"] for x in own_cis),
                   ORTHOGONAL_INTENDED_BEST_CIS_P=own_best["p"] if own_best else None,
                   ORTHOGONAL_INTENDED_BEST_CIS_COHORT=own_best["dataset"] if own_best else None,
                   HIGH_SIMILARITY=int(result["shared"] >= 4 and result["jaccard"] >= .6),
                   INDEPENDENT_VALIDATION=0, WRONG_PROTEIN_CONFIRMED=0,
                   INTERPRETATION="same_cohort_descriptive_similarity_not_independent_confirmation")
        for role in ("intended", "alternative"):
            value = best[role]; prefix = "BEST_ORTHOGONAL_" + role.upper()
            row.update({prefix + "_COHORT": value[0] if value else None, prefix + "_ASSAY": value[1] if value else None,
                        prefix + "_SHARED_REGIONS": value[2]["shared"] if value else None,
                        prefix + "_EXACT_RSID": value[2]["exact"] if value else None,
                        prefix + "_JACCARD": value[2]["jaccard"] if value else None})
        all_rows.append(row)
        all_matches.extend(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, PROT=name, PRIMARY=dataset, PRIMARY_ASSAY=alias,
                                ALTERNATIVE_PROT=target, ALTERNATIVE_ASSAY=donor, MODE="FULL", **m) for m in result["matches"])
        if row["HIGH_SIMILARITY"]:
            all_matches.extend(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, PROT=name, PRIMARY=dataset, PRIMARY_ASSAY=alias,
                                    ALTERNATIVE_PROT=target, ALTERNATIVE_ASSAY=donor, MODE="REGION_HELDOUT", **m) for m in held["matches"])
            case_keys = [("primary_intended", left_key), ("primary_alternative", right_key)]
            for role, protein in (("external_intended", name), ("external_alternative", target)):
                for cohort in DATASETS:
                    if cohort == dataset: continue
                    for a in lookup[protein]["data"][cohort].get("assays", []):
                        if spec.summary_valid(a): case_keys.append((role, (protein, cohort, a["alias"])))
            for role, key in case_keys:
                for sig in sorted(by_assay.get(key, []), key=lambda s: (spec.chrom_key(s["chr"]), s["pos"], s["snp"])):
                    profile_rows.append(dict(PAIR_ID=pair_id, HYPOTHESIS_ID=hid, ROLE=role, PROT=key[0], DATASET=key[1], ASSAY=key[2],
                                             BUILD=sig["build"], CHR=sig["chr"], POS=sig["pos"], SNP=sig["snp"], P=sig["p"], CLASS=sig["kind"],
                                             IN_STRICT_PROFILE=int(sig in profiles.get(key, [])), QC_USABLE=int(sig["clean"]), FLAGS=joined(sig["flags"])))
    all_rows.sort(key=lambda r: (-r["HIGH_SIMILARITY"], -r["JACCARD"], -r["SHARED_REGIONS"], r["PRIMARY"], r["PROT"], r["ALTERNATIVE_PROT"], r["PRIMARY_ASSAY"], r["ALTERNATIVE_ASSAY"]))
    return all_rows, all_matches, profile_rows, len(eligible)


BASE_HYPOTHESIS_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "PROT", "PRIMARY", "PRIMARY_PLATFORM", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "REFERENCE", "REFERENCE_PLATFORM", "REFERENCE_ASSAY", "SCREEN_STRATUM", "TIER", "PRIOR_PRIORITY_CASE", "PRIMARY_STATE", "PRIMARY_CIS_P", "PRIMARY_CIS_SEARCH_COMPLETE", "PRIMARY_ASSAYS_EXPECTED", "PRIMARY_ASSAYS_WITH_STRICT_CIS", "REFERENCE_INTENDED_ASSAYS_EXPECTED", "REFERENCE_INTENDED_ASSAYS_USABLE", "REFERENCE_INTENDED_ASSAYS_WITH_STRICT_CIS", "SAME_REFERENCE_INTENDED_STRICT_CIS", "REFERENCE_INTENDED_CIS_P", "REFERENCE_INTENDED_CIS_SNP", "REFERENCE_INTENDED_CIS_ASSAY", "INTENDED_STRICT_CIS_COHORTS", "ALTERNATIVE_STRICT_CIS_COHORTS", "N_ANCHOR_NOMINATIONS", "ANCHOR_SPECIFIC", "EXACT_ANCHOR_MULTILOCUS", "ANCHOR_SNP", "ANCHOR_CHR", "ANCHOR_POS", "ANCHOR_PRIMARY_P", "ANCHOR_REFERENCE_SNP", "ANCHOR_REFERENCE_P", "ANCHOR_LEAD_MATCH", "ANCHOR_PRIMARY_EXACT_TARGETS", "ANCHOR_REFERENCE_EXACT_TARGETS", "ANCHOR_PRIMARY_REGION_TARGETS", "ANCHOR_REFERENCE_REGION_TARGETS", "COMMON_UNIVERSE_TARGETS", "REGION_EXCLUDED", "BLOCK5MB_EXCLUDED", "CHROMOSOMES_EXCLUDED"]
HYPOTHESIS_FIELDS = BASE_HYPOTHESIS_FIELDS + [m + "_" + f for m in MODES for f in METRICS] + ["STRICT_TIER_NOT_MET_REASONS", "CAUTIONS", "INTERPRETATION", "DIRECTION_TESTED", "COLOCALIZATION_TESTED", "WRONG_PROTEIN_CONFIRMED"]
PAIR_FIELDS = HYPOTHESIS_FIELDS + ["REPRESENTATIVE_HYPOTHESIS_ID", "N_ASSAY_HYPOTHESES", "N_PRIMARY_ASSAYS", "N_REFERENCE_COHORTS", "N_REFERENCE_ASSAYS", "REFERENCE_COHORTS", "SCREEN_STRATA", "PROFILE_METRIC_SCOPE"]
MATCH_EXTRA_FIELDS = ["PRIMARY_EXACT_TARGETS", "PRIMARY_REGION_TRANS_TARGETS", "REFERENCE_EXACT_TARGETS", "REFERENCE_REGION_TRANS_TARGETS", "NONHOTSPOT_PRIMARY", "NONHOTSPOT_REFERENCE", "NONHOTSPOT_BOTH", "COMMON_UNIVERSE_TARGETS", "COMMON_PRIMARY_EXACT_TARGETS", "COMMON_PRIMARY_REGION_TRANS_TARGETS", "COMMON_REFERENCE_EXACT_TARGETS", "COMMON_REFERENCE_REGION_TRANS_TARGETS", "COMMON_PRIMARY_EXACT_FRACTION", "COMMON_REFERENCE_EXACT_FRACTION", "COMMON_PRIMARY_REGION_FRACTION", "COMMON_REFERENCE_REGION_FRACTION", "NONHOTSPOT_BOTH_COMMON_UNIVERSE"]
MATCH_FIELDS = ["PAIR_ID", "MODE"] + spec.FINGERPRINT_FIELDS + MATCH_EXTRA_FIELDS
PROBE_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "MODE", "COMPARISON", "PROT", "PRIMARY", "PRIMARY_ASSAY", "REFERENCE_PROT", "REFERENCE", "REFERENCE_ASSAY", "PRIMARY_LOCI", "REFERENCE_LOCI", "SHARED_REGIONS", "EXACT_RSID", "JACCARD", "NONHOTSPOT_PRIMARY", "NONHOTSPOT_BOTH", "NONHOTSPOT_BOTH_COMMON_UNIVERSE"]
ANCHOR_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "PROT", "PRIMARY", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "REFERENCE", "REFERENCE_ASSAY", "SCREEN_STRATUM", "BUILD", "CHR", "PRIMARY_POS", "PRIMARY_SNP", "PRIMARY_P", "ALTERNATIVE_CIS_POS", "ALTERNATIVE_CIS_SNP", "ALTERNATIVE_CIS_P", "LEAD_MATCH", "LEAD_DISTANCE_BP", "PRIMARY_EXACT_TARGETS", "PRIMARY_REGION_TRANS_TARGETS", "REFERENCE_EXACT_TARGETS", "REFERENCE_REGION_TRANS_TARGETS", "ALTERNATIVE_TARGETS_AT_LOCUS", "EXACT_ALTERNATIVE_TARGETS_AT_VARIANT", "INTENDED_TARGET_TRANS_COHORTS", "ANCHOR_SPECIFIC", "CAUTIONS"]
EXCLUSION_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "PROT", "PRIMARY", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "REFERENCE", "REFERENCE_ASSAY", "SCREEN_STRATUM", "TIER", "REASONS"]
CONTEXT_ID_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "CONTEXT_ID", "COMPARISON", "LEFT_PROT", "LEFT_DATASET", "LEFT_ASSAY", "RIGHT_PROT", "RIGHT_DATASET", "RIGHT_ASSAY"]
CONTEXT_FIELDS = CONTEXT_ID_FIELDS + ["LEFT_LOCI", "RIGHT_LOCI", "SHARED_REGIONS", "EXACT_RSID", "JACCARD", "NONHOTSPOT_BOTH", "COMMON_UNIVERSE_TARGETS", "INDEPENDENT_VALIDATION", "INTERPRETATION"]
CONTEXT_MATCH_FIELDS = CONTEXT_ID_FIELDS + ["BUILD", "CHR", "PRIMARY_POS", "PRIMARY_SNP", "PRIMARY_P", "PRIMARY_CLASS", "REFERENCE_POS", "REFERENCE_SNP", "REFERENCE_P", "REFERENCE_CLASS", "EXACT_RSID", "DISTANCE_BP", "NONHOTSPOT"] + MATCH_EXTRA_FIELDS
PROFILE_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "ROLE", "PROT", "DATASET", "ASSAY", "BUILD", "CHR", "POS", "SNP", "P", "CLASS", "IN_STRICT_PROFILE", "QC_USABLE", "FLAGS"]
SENSITIVITY_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "PROT", "PRIMARY", "ALTERNATIVE_PROT", "REFERENCE", "MATCH_SPAN_BP", "MODE", "FROZEN_EXCLUSION_SPAN_BP"] + ["PROFILE_" + m for m in METRICS]

SAME_COHORT_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "PROT", "PRIMARY", "PRIMARY_PLATFORM", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "ALTERNATIVE_ASSAY", "SCREEN_STRATUM", "PRIMARY_STATE", "PRIMARY_CIS_P", "PRIMARY_CIS_SEARCH_COMPLETE", "ANCHOR_CHR", "ANCHOR_POS", "ANCHOR_SNP", "ANCHOR_PRIMARY_P", "ALTERNATIVE_CIS_SNP", "ALTERNATIVE_CIS_P", "ANCHOR_MATCH", "N_ANCHOR_NOMINATIONS", "PRIMARY_LOCI", "ALTERNATIVE_LOCI", "SHARED_REGIONS", "EXACT_RSID", "JACCARD", "NONHOTSPOT_BOTH", "REGION_HELDOUT_SHARED_REGIONS", "REGION_HELDOUT_EXACT_RSID", "REGION_HELDOUT_JACCARD", "REGION_HELDOUT_NONHOTSPOT_BOTH", "REGION_EXCLUDED", "ORTHOGONAL_INTENDED_CIS_COHORTS", "ORTHOGONAL_INTENDED_BEST_CIS_P", "ORTHOGONAL_INTENDED_BEST_CIS_COHORT"] + ["BEST_ORTHOGONAL_" + role + "_" + metric for role in ("INTENDED", "ALTERNATIVE") for metric in ("COHORT", "ASSAY", "SHARED_REGIONS", "EXACT_RSID", "JACCARD")] + ["HIGH_SIMILARITY", "INDEPENDENT_VALIDATION", "WRONG_PROTEIN_CONFIRMED", "INTERPRETATION"]
SAME_COHORT_MATCH_FIELDS = ["PAIR_ID", "HYPOTHESIS_ID", "PROT", "PRIMARY", "PRIMARY_ASSAY", "ALTERNATIVE_PROT", "ALTERNATIVE_ASSAY", "MODE"] + ["BUILD", "CHR", "PRIMARY_POS", "PRIMARY_SNP", "PRIMARY_P", "PRIMARY_CLASS", "REFERENCE_POS", "REFERENCE_SNP", "REFERENCE_P", "REFERENCE_CLASS", "EXACT_RSID", "DISTANCE_BP", "NONHOTSPOT"] + MATCH_EXTRA_FIELDS


def run(cross_dir, protein_map, specificity_dir, out, primary="all", include_weak_cis=True):
    payload, rows, seeds, source_summary, inputs = load_specificity(cross_dir, protein_map, specificity_dir, primary)
    span = int(source_summary["rules"]["profile_locus_span_bp"])
    result = analyse(payload, rows, seeds, primary, span, include_weak_cis)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    tables = [("pair-summary.tsv", PAIR_FIELDS, "pairs"), ("assay-evidence.tsv", HYPOTHESIS_FIELDS, "hypotheses"),
              ("anchor-evidence.tsv", ANCHOR_FIELDS, "anchors"), ("matched-loci.tsv", MATCH_FIELDS, "matched_loci"),
              ("probe-comparisons.tsv", PROBE_FIELDS, "probe_comparisons"), ("exclusions.tsv", EXCLUSION_FIELDS, "exclusions"),
              ("same-platform-context.tsv", CONTEXT_FIELDS, "same_platform_context"), ("context-matches.tsv", CONTEXT_MATCH_FIELDS, "context_matches"),
              ("profile-signals.tsv", PROFILE_FIELDS, "profile_signals"), ("match-span-sensitivity.tsv", SENSITIVITY_FIELDS, "sensitivity"),
              ("same-cohort-eligible.tsv", SAME_COHORT_FIELDS, "same_cohort_eligible"),
              ("same-cohort-screen.tsv", SAME_COHORT_FIELDS, "same_cohort_screen"),
              ("same-cohort-screen-matches.tsv", SAME_COHORT_MATCH_FIELDS, "same_cohort_screen_matches"),
              ("same-cohort-profile-signals.tsv", PROFILE_FIELDS, "same_cohort_profile_signals")]
    for filename, fields, key in tables:
        spec.write_tsv(out / filename, fields, result[key])
    summary = result["summary"]
    summary.update(generated=datetime.now(timezone.utc).isoformat(timespec="seconds"), input_sha256=inputs,
                   code_sha256={"fingerprint_dossiers.py": spec.digest(__file__), "assay_specificity.py": spec.digest(spec.__file__)})
    # Avoid a circular JSON hash: data.summary equals the final summary with
    # output_sha256 removed. The receipt hashes every TSV and the data JSON.
    data = {"summary": dict(summary), "data_display_scope": summary["casebook_scope"], "pairs": result["pairs"]}
    for key in ("hypotheses", "anchors", "matched_loci", "probe_comparisons", "same_platform_context", "context_matches", "profile_signals", "sensitivity"):
        data[key] = [r for r in result[key] if r["HYPOTHESIS_ID"] in result["display_ids"]]
    data["same_cohort_screen"] = result["same_cohort_screen"]
    same_ids = {r["HYPOTHESIS_ID"] for r in result["same_cohort_screen"]}
    data["same_cohort_screen_matches"] = [r for r in result["same_cohort_screen_matches"] if r["HYPOTHESIS_ID"] in same_ids]
    data["same_cohort_profile_signals"] = result["same_cohort_profile_signals"]
    spec.write_json(out / "fingerprint-dossiers.data.json", data)
    summary["output_sha256"] = {filename: spec.digest(out / filename) for filename, _, _ in tables}
    summary["output_sha256"]["fingerprint-dossiers.data.json"] = spec.digest(out / "fingerprint-dossiers.data.json")
    spec.write_json(out / "fingerprint-dossiers.summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("cross-dir", "protein-map", "specificity-dir", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--primary", choices=[*DATASETS, "all"], default="all")
    parser.add_argument("--no-weak-cis", action="store_true", help="Restrict to the original no-reported-cis universe")
    args = parser.parse_args()
    print(json.dumps(run(args.cross_dir, args.protein_map, args.specificity_dir, args.out, args.primary, not args.no_weak_cis), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
