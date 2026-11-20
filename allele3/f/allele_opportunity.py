#!/usr/bin/env python3
"""Audit three possible nonreference substitutions against reported AC/AN.

This is a bounded, exact join of existing target and Atlas derivatives. It does
not reread the 8.8-billion-row Atlas, reconstruct REF counts from split summaries,
or treat unreported nucleotides as extinct. All summaries are conditional on the
input sampling frame. AVI comparisons are descriptive because AVI was trained
using gnomAD frequency-derived labels. SHAP attributions are never converted to
functional breadth, signs of biological regulation, or allele distances.
"""
from __future__ import annotations

# Allele Opportunity
# ======================================================================
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import statistics

from analysis_integrity import identity, transactional

VERSION = "allele-opportunity-1.0"
BASES = "ACGT"
MISSING = {"", ".", "NA", "NAN", "NONE", "NULL"}
ATTRIBUTIONS = (
    "MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS", "MAX_ABS_DNASE",
    "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CAGE", "MAX_ABS_PROCAP",
    "MAX_ABS_RNA_SEQ", "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE",
    "CACTUS_241_WAY", "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST",
    "PHASTCONS_470_WAY", "IS_INSERTION", "IS_DELETION",
)
LIMITS = [
    "Targets define the sampling frame. A sample selected for multi-ALT status cannot estimate the genome-wide probability of acquiring another allele.",
    "Count support is evaluated independently for each reported ALT. REF counts and shared four-state frequencies cannot be reconstructed from split AC/AN, even when AN happens to agree.",
    "not_in_source_rows is missing observation evidence, not AC=0, callability, failure to arise, or extinction. A reported zero count remains a distinct state.",
    "Nonreference AF ranks are current ranks, not mutation order. Secondary ALT is not automatically the third observed nucleotide: REF presence still requires joint counts.",
    "AVI uses gnomAD frequency-derived training labels. Its relation to gnomAD rarity/count support is descriptive and potentially circular, including AVI SHAP attributions.",
    "The splicing score is a reference-relative aggregate with incomplete coverage. Missing values remain NA, not zero. It is not a signed tissue-specific effect or allele-to-allele distance.",
    "AVI feature attributions explain the impact model; their signs do not identify up/down regulation and their number is not phenotypic pleiotropy.",
    "Matching, repeated controls, ancestry, LD, mutation supply, and selection are not modeled by these descriptive summaries. No inferential P values are emitted.",
    "Input coordinate/REF consistency is checked; this derivative join does not independently check the genome FASTA, index integrity, or cohort independence.",
]


def number(value, name, *, integer=False, lower=None, upper=None, allow_missing=False):
    if str(value).strip().upper() in MISSING:
        if allow_missing:
            return None
        raise ValueError(f"Missing {name}")
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric {name}: {value!r}") from exc
    if not math.isfinite(x) or (integer and x != int(x)) or (lower is not None and x < lower) or (upper is not None and x > upper):
        raise ValueError(f"Invalid {name}: {value!r}")
    return int(x) if integer else x


def table(path, required):
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or set(required) - set(fields):
            raise ValueError(f"{path}: missing/duplicate columns; required {sorted(required)}")
        for i, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"{path}:{i}: malformed tab-separated row")
            yield row


def chromosome(value):
    c = str(value).strip().removeprefix("chr").upper()
    c = "MT" if c == "M" else c
    if c not in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
        raise ValueError(f"Invalid chromosome {value!r}")
    return c


def key(row):
    ref = row["ref"].upper()
    if ref not in BASES or len(ref) != 1:
        raise ValueError("Only exact single-nucleotide REF values are accepted")
    if row.get("build", "GRCh38").lower() not in {"grch38", "hg38", "38"}:
        raise ValueError("Every input must be GRCh38; no relabeling or implicit liftOver")
    return chromosome(row["chrom"]), number(row["pos"], "POS", integer=True, lower=1), ref


def load_targets(path, max_targets):
    required = {"site_id", "chrom", "pos", "ref", "alts", "observed_alts_all", "alt_ac_all", "alt_an_all"}
    targets, ids, positions = {}, {}, set()
    for r in table(path, required):
        k = key(r)
        sid = r["site_id"]
        if sid.strip().upper() in MISSING or sid in ids or k[:2] in positions:
            raise ValueError("Target IDs and genomic positions must be unique")
        aa = r["observed_alts_all"].upper().split(",")
        ac, an = r["alt_ac_all"].split(","), r["alt_an_all"].split(",")
        if not 1 <= len(aa) <= 3 or len(aa) != len(set(aa)) or any(len(a) != 1 or a not in BASES or a == k[2] for a in aa):
            raise ValueError(f"{sid}: invalid observed_alts_all")
        if len(ac) != len(aa) or len(an) != len(aa):
            raise ValueError(f"{sid}: ALT/AC/AN vectors have different lengths")
        counts = {}
        for a, c, n in zip(aa, ac, an):
            c = number(c, f"{sid}/{a}/AC", integer=True, lower=0)
            n = number(n, f"{sid}/{a}/AN", integer=True, lower=1)
            if c > n:
                raise ValueError(f"{sid}: AC exceeds AN")
            counts[a] = dict(ac=c, an=n, af=c / n)
        listed = r["alts"].upper().split(",")
        if not listed or len(listed) != len(set(listed)) or not set(listed) <= set(aa):
            raise ValueError(f"{sid}: target alts must be a unique subset of reported source ALTs")
        r["alts"] = ",".join(listed)
        ranks = sorted(counts, key=lambda a: (-counts[a]["af"], a))
        # AF ranks are allowed with differing AN, but that limitation stays explicit.
        r.update(_key=k, _counts=counts, _rank=ranks)
        targets[k] = r
        ids[sid] = k
        positions.add(k[:2])
        if len(targets) > max_targets:
            raise ValueError(f"Target cap {max_targets} exceeded; shard targets rather than loading a whole genome")
    if not targets:
        raise ValueError("No target rows")
    return targets, ids


def load_variants(path, targets, ids):
    required = {"site_id", "chrom", "pos", "build", "ref", "allele", "allele_role", "avi_raw", "avi_phred", "status", "resource_id", "source_version"}
    result, versions = {}, set()
    for r in table(path, required):
        k = key(r)
        if k not in targets or ids.get(r["site_id"]) != k:
            raise ValueError("Atlas and target identities do not match exactly; use derivatives for the same targets")
        a = r["allele"]
        if len(a) != 1 or a not in BASES or (k, a) in result:
            raise ValueError("Atlas requires unique exact A/C/G/T rows")
        if r["allele_role"] != ("REF" if a == k[2] else "ALT"):
            raise ValueError("Atlas allele_role disagrees with REF")
        raw = number(r["avi_raw"], "AVI raw", allow_missing=True)
        phred = number(r["avi_phred"], "AVI PHRED", lower=0, allow_missing=True)
        if a == k[2] and (raw is not None or phred is not None):
            raise ValueError("REF is a sequence label, not a substitution with an AVI score")
        if a != k[2]:
            n = int(raw is not None) + int(phred is not None)
            expected = {2: {"available"}, 1: {"available_partial"}, 0: {"not_scored", "not_indexed_contig", "score_missing"}}[n]
            if r["status"] not in expected:
                raise ValueError("AVI availability status contradicts numeric data")
        r.update(_raw=raw, _phred=phred)
        versions.add((r["resource_id"], r["source_version"]))
        result[k, a] = r
    if len(versions) != 1:
        raise ValueError("AVI contrasts require a single resource and source_version")
    if len(result) != 4 * len(targets):
        raise ValueError("Atlas must contain all four base rows for every target, including explicit REF/missing rows")
    return result


def load_metrics(path, targets, ids):
    result, seen, attribution_values, splicing_versions = {}, set(), [], set()
    if not path:
        return result, attribution_values
    required = {"site_id", "chrom", "pos", "build", "ref", "allele", "resource_id", "source_version", "metric", "value", "status"}
    for r in table(path, required):
        k = key(r)
        if k not in targets or ids.get(r["site_id"]) != k:
            raise ValueError("Atlas metrics and targets do not match exactly")
        a, metric = r["allele"], r["metric"]
        if a not in BASES or len(a) != 1:
            raise ValueError("Invalid metric allele")
        item = k, a, r["resource_id"], metric
        if item in seen:
            raise ValueError("Duplicate resource/metric/allele row")
        seen.add(item)
        value = number(r["value"], metric, allow_missing=True)
        if a == k[2] and value is not None:
            raise ValueError("REF substitution metric must be missing")
        if value is not None and r["status"] != "available":
            raise ValueError("Metric status contradicts numeric value")
        if value is None and r["status"] == "available":
            raise ValueError("Available metric is missing")
        if metric == "alphagenome_splicing":
            if (k, a) in result:
                raise ValueError("More than one splicing resource; choose a single comparable version")
            if value is not None and value < 0:
                raise ValueError("The aggregate splicing impact score must be nonnegative")
            semantic = r.get("metric_semantics", "")
            if semantic and ("shap" in semantic.lower() or "attribution" in semantic.lower()):
                raise ValueError("Splicing score cannot be supplied as an AVI attribution")
            result[k, a] = dict(value=value, status=r["status"], source_version=r["source_version"])
            splicing_versions.add((r["resource_id"], r["source_version"]))
        elif r["resource_id"] == "avi_attributions":
            if metric not in ATTRIBUTIONS:
                raise ValueError(f"Unknown AVI attribution {metric}; review its source semantics explicitly")
            attribution_values.append(dict(site_id=r["site_id"], chrom=k[0], pos=k[1], ref=k[2], allele=a,
                metric=metric, value=value, status=r["status"], source_version=r["source_version"],
                semantics="signed_AVI_feature_attribution_not_molecular_direction"))
    if len(splicing_versions) > 1:
        raise ValueError("Splicing contrasts require a single resource and source_version")
    return result, attribution_values


def substitution_class(ref, alt):
    if ref == alt:
        return "reference_label"
    return "transition" if {ref, alt} in ({"A", "G"}, {"C", "T"}) else "transversion"


def reported_status(count, threshold):
    if count is None:
        return "not_in_source_rows"
    if count["ac"] == 0:
        return "reported_zero_count"
    return "count_supported" if count["ac"] >= threshold else "reported_below_threshold"


def write_tsv(path, records, fields=None):
    records = list(records)
    fields = fields or list(dict.fromkeys(k for r in records for k in r)) or ["status"]
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in records:
            w.writerow({k: "NA" if r.get(k) is None else r.get(k, "NA") for k in fields})


def comparisons(target, variants, metrics, threshold):
    k, counts = target["_key"], target["_counts"]
    supported = [a for a in target["_rank"] if counts[a]["ac"] >= threshold]
    common = dict(site_id=target["site_id"], chrom=k[0], pos=k[1], ref=k[2],
                  family=target.get("family", "all"), group=target.get("group", "unstratified"),
                  partition=target.get("partition", "unassigned"), threshold=threshold,
                  n_supported_nonreference=len(supported), ref_presence="not_reconstructed",
                  role_interpretation="current_nonreference_AF_rank_not_mutation_order")
    if len(supported) < 2:
        return [dict(common, contrast="no_supported_nonreference_pair", status="not_estimable")]
    first, second = supported[:2]
    pairs = [("second_minus_first", second, first)]
    pairs += [("unreported_minus_second", a, second) for a in BASES if a != k[2] and a not in counts]
    pairs += [("below_threshold_minus_second", a, second) for a in target["_rank"] if a not in supported]
    if len(supported) == 3:
        pairs += [("third_nonreference_minus_second", supported[2], second)]
    rr = []
    for name, left, right in pairs:
        l, r = variants[k, left], variants[k, right]
        ls, rs = metrics.get((k, left), {}).get("value"), metrics.get((k, right), {}).get("value")
        lp, rp = l["_phred"], r["_phred"]
        rr.append(dict(common, contrast=name, status="descriptive_selected_targets", left_allele=left,
            right_allele=right, left_AC=counts.get(left, {}).get("ac"), right_AC=counts[right]["ac"],
            left_AN=counts.get(left, {}).get("an"), right_AN=counts[right]["an"],
            left_source_status=reported_status(counts.get(left), threshold), right_source_status="count_supported",
            primary_secondary_AF_tie=int(counts[first]["af"] == counts[second]["af"]),
            source_AN_equal=int(len({c["an"] for c in counts.values()}) == 1),
            left_substitution=substitution_class(k[2], left), right_substitution=substitution_class(k[2], right),
            left_avi_phred=lp, right_avi_phred=rp, avi_phred_difference=lp-rp if lp is not None and rp is not None else None,
            left_splicing=ls, right_splicing=rs, splicing_difference=ls-rs if ls is not None and rs is not None else None,
            AVI_use="descriptive_frequency_label_circularity", difference_meaning="score_difference_not_allelic_functional_distance"))
    return rr


def summarize(contrasts):
    groups = defaultdict(list)
    for r in contrasts:
        groups[r["family"], r["group"], r["threshold"], r["contrast"]].append(r)
    rows = []
    for (fam, group, threshold, contrast), rr in sorted(groups.items()):
        out = dict(family=fam, group=group, threshold=threshold, contrast=contrast,
                   n_sites=len(rr), n_AF_ties=sum(r.get("primary_secondary_AF_tie", 0) for r in rr),
                   estimand="unweighted_descriptive_summary_conditional_on_selected_targets")
        for score in ("avi_phred", "splicing"):
            pairs = [r for r in rr if r.get(score+"_difference") is not None]
            dd = [r[score+"_difference"] for r in pairs]
            out.update({score+"_paired_n": len(dd), score+"_missing_pair_n": len(rr)-len(dd),
                        score+"_median_left": statistics.median(r["left_"+score] for r in pairs) if pairs else None,
                        score+"_median_right": statistics.median(r["right_"+score] for r in pairs) if pairs else None,
                        score+"_median_within_site_difference": statistics.median(dd) if dd else None,
                        score+"_median_absolute_within_site_difference": statistics.median(abs(d) for d in dd) if dd else None,
                        score+"_left_greater_fraction": sum(d > 0 for d in dd)/len(dd) if dd else None,
                        score+"_equal_fraction": sum(d == 0 for d in dd)/len(dd) if dd else None})
        rows.append(out)
    return rows


@transactional("manifest.json", input_files=("targets", "atlas_variants", "atlas_metrics"))
def run(args):
    if args.build != "GRCh38" or args.max_targets < 1:
        raise ValueError("Require GRCh38 and a positive bounded max-targets")
    thresholds = tuple(sorted(set(int(x) for x in args.thresholds.split(","))))
    if not thresholds or min(thresholds) < 1:
        raise ValueError("Thresholds must be positive integer allele counts")
    targets, ids = load_targets(args.targets, args.max_targets)
    variants = load_variants(args.atlas_variants, targets, ids)
    metrics, attribution_rows = load_metrics(args.atlas_metrics, targets, ids)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    opportunities, contrasts, audit = [], [], Counter()
    for k, t in targets.items():
        counts = t["_counts"]
        audit["target_sites"] += 1
        audit[f"sites_with_{len(counts)}_reported_nonreference_bases"] += 1
        audit["sites_with_differing_source_AN"] += len({c["an"] for c in counts.values()}) > 1
        for a in BASES:
            v, s = variants[k, a], metrics.get((k, a), {})
            row = dict(site_id=t["site_id"], chrom=k[0], pos=k[1], build="GRCh38", ref=k[2], allele=a,
                family=t.get("family", "all"), group=t.get("group", "unstratified"), partition=t.get("partition", "unassigned"),
                match_id=t.get("match_id", "NA"), allele_role="REF" if a == k[2] else "ALT",
                target_membership=int(a in t["alts"].split(",")), source_row_reported=int(a in counts) if a != k[2] else None,
                AC=counts.get(a, {}).get("ac"), AN=counts.get(a, {}).get("an"), AF=counts.get(a, {}).get("af"),
                nonreference_AF_rank=t["_rank"].index(a)+1 if a in counts else None,
                nonreference_AF_rank_tied=int(sum(c["af"] == counts[a]["af"] for c in counts.values()) > 1) if a in counts else None,
                AN_common_across_reported_rows=int(len({c["an"] for c in counts.values()}) == 1),
                ref_presence="not_reconstructed_from_split_summary", substitution=substitution_class(k[2], a),
                context=t.get("context", "NA"), avi_raw=v["_raw"], avi_phred=v["_phred"], avi_status=v["status"],
                alphagenome_splicing=s.get("value"), splicing_status=s.get("status", "not_in_derivative"),
                atlas_source_version=v["source_version"], observation_source=args.source_cohort,
                source_sampling=args.source_sampling, AVI_use="descriptive_frequency_label_circularity")
            for threshold in thresholds:
                row[f"AC{threshold}_status"] = "reference_unmeasured" if a == k[2] else reported_status(counts.get(a), threshold)
                if a != k[2]:
                    audit[f"AC{threshold}_"+row[f"AC{threshold}_status"]] += 1
            if a != k[2]:
                audit["theoretical_nonreference_substitutions"] += 1
                audit["AVI_scored_substitutions"] += v["_phred"] is not None
                audit["splicing_scored_substitutions"] += s.get("value") is not None
            opportunities.append(row)
        for threshold in thresholds:
            contrasts.extend(comparisons(t, variants, metrics, threshold))
    write_tsv(out / "nucleotide_opportunities.tsv.gz", opportunities)
    write_tsv(out / "within_site_score_contrasts.tsv.gz", contrasts)
    write_tsv(out / "contrast_summary.tsv", summarize(contrasts))
    write_tsv(out / "observation_audit.tsv", [dict(metric=k, value=v) for k, v in sorted(audit.items())])
    write_tsv(out / "avi_attributions.tsv.gz", attribution_rows,
              ["site_id", "chrom", "pos", "ref", "allele", "metric", "value", "status", "source_version", "semantics"])
    manifest = dict(status="complete", schema=VERSION, build="GRCh38", counts=dict(audit), thresholds=thresholds,
                    source_sampling=args.source_sampling, source_cohort=args.source_cohort,
                    sources={n: identity(getattr(args, n), full=True) for n in ("targets", "atlas_variants", "atlas_metrics") if getattr(args, n)},
                    rank_semantics="current_nonreference_AF_order_not_birth_order", limits=LIMITS,
                    independent_replication_established=False, genome_wide_birth_probability_estimated=False,
                    pleiotropy_measured=False, selection_test_performed=False)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--targets", required=True, help="Candidates or another count-bearing panel with observed_alts_all,alt_ac_all,alt_an_all")
    p.add_argument("--atlas-variants", required=True, help="Existing atlas_variants.tsv.gz, four nucleotide rows/site")
    p.add_argument("--atlas-metrics", help="Existing atlas_metrics.tsv.gz; missing optional splicing/attributions remain missing")
    p.add_argument("--build", required=True, choices=["GRCh38"])
    p.add_argument("--source-cohort", required=True, help="Explicit observation-cohort/version label, not an independence assertion")
    p.add_argument("--source-sampling", required=True, choices=["selected_multialt_candidates", "prespecified_panel", "all_eligible_source_sites"])
    p.add_argument("--thresholds", default="5,10")
    p.add_argument("--max-targets", type=int, default=100000)
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    result = run(args)
    print(json.dumps({"status": result["status"], "target_sites": result["counts"]["target_sites"], "out": args.out}))
    return 0


# Opportunity Model
# ======================================================================
_description_opportunity_model = 'Within-site association between molecular predictors and count-supported ALTs.\n\nRequires all three possible nonreference bases with jointly measured WGS counts\nand shared four-base callability. The exact conditional likelihood conditions on\nthe number of supported ALTs, not on their identities. It is an association\nmodel, not a test of mutation birth, fitness or balancing selection.\n\nThe coefficient and cluster_robust_se columns use log-odds units per within-site\npredictor SD; ci95_low and ci95_high are exponentiated odds-ratio limits.\n'

import itertools
import re

import numpy as np
from scipy import optimize, special, stats


MODEL_VERSION = "opportunity-conditional-logit-1.0"
MODEL_BASES = set("ACGT")
MODEL_MISSING = {"", ".", "NA", "NAN", "NONE", "NULL"}
MODEL_DEFAULT_PREDICTORS = "molecular_breadth,effect_magnitude,log_mutation_rate"
MODEL_REQUIRED = {"site_id", "chrom", "pos", "ref", "alt", "build", "population", "cohort",
            "AC", "AN", "ref_AC", "callable_all_four", "ld_block", "mutation_rate"}
MODEL_COEFFICIENT_FIELDS = ["predictor", "within_site_sd", "coefficient", "odds_ratio_per_within_site_sd",
                      "cluster_robust_se", "ci95_low", "ci95_high", "p_value", "q_value_within_model",
                      "inference_status"]
MODEL_LIMITATIONS = [
    "The outcome is AC>=min_ac for each of the three possible nonreference nucleotides, conditional on the site/cohort/population's number of supported ALTs. This is not a minor-allele rank or mutation-order outcome.",
    "Joint WGS counts, four-base callability and raw/experimental molecular provenance are explicit caller declarations plus arithmetic checks. A TSV alone cannot verify the underlying sequencing or rule out renamed AVI/SHAP values.",
    "Molecular breadth must be prespecified using comparable molecular tracks, effect thresholds/calibration and a treatment of correlated traits. A count of nonzero AVI feature attributions does not meet this input contract.",
    "REF is a reference-sequence label, not necessarily the ancestral or major allele. All four counts must sum to one shared AN; no REF count or unobserved ALT count is imputed.",
    "Predictors are centered within each three-ALT set. Scaling uses all complete callable sets before inspecting their k=0/1/2/3 support pattern; between-site-only predictors are unidentifiable.",
    "The mutation-rate coefficient is estimated, not fixed as an offset. Its supplied sequence context and direction must match the hypothesis; a REF-relative rate is not automatically an ancestral mutation rate.",
    "The model assumes the stated conditional-logit choice probabilities. Count thresholds, mutation-rate measurement error, shared demography, ascertainment and residual technical variation may affect the association.",
    "The supplied LD-block IDs must define one consistent genomic partition across populations and cohorts. Scores are clustered by chromosome and block, including repeated sites across populations/cohorts. Cohort overlap and long-range dependence still require study-specific assessment.",
    "Wald inference uses a cluster sandwich with CR1 finite-sample correction and t reference with G-1 degrees of freedom, only for >=20 blocks and an identifiable, converged, nonseparated fit. The coefficient and cluster_robust_se columns are on the log-odds scale; ci95_low and ci95_high are exponentiated odds-ratio confidence limits per within-site predictor SD. BH adjustment covers this model's predictors only, not other thresholds or analyses.",
    "Prespecified target sets support target-conditional association only; they do not identify genome-wide rates. Ascertainment on multi-ALT status is refused for this primary model. Selection on unknown or related features can still bias a declared sampling frame.",
    "Association of molecular breadth with count support cannot establish pleiotropy-mediated survival, balancing selection, causal fitness effects or mutation birth times. Independent-cohort replication and additional population-genetic evidence are required.",
]


def model_predictor_names(value):
    names = [v.strip() for v in value.split(",")]
    if not names or any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", n) for n in names) or len(set(names)) != len(names):
        raise ValueError("--predictors must list distinct valid column names")
    if len(names) > 20:
        raise ValueError("At most 20 prespecified predictors are supported")
    for name in names:
        lower, collapsed = name.lower(), re.sub(r"[^a-z0-9]", "", name.lower())
        if any(v in collapsed for v in ("shap", "phred", "featureimportance")) or re.search(r"(^|_)avi($|_)", lower) or lower.startswith("avi"):
            raise ValueError("AVI/SHAP/PHRED/feature_importance predictors are forbidden: " + name)
        if lower in {v.lower() for v in MODEL_REQUIRED - {"mutation_rate"}}:
            raise ValueError("Count/outcome/identity columns cannot be predictors: " + name)
    return names


def model_number(value, label, integer=False):
    if value.strip().upper() in MODEL_MISSING:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: expected a finite number or explicit NA") from exc
    if not math.isfinite(value) or integer and (value < 0 or value > 2 ** 53 or value != int(value)):
        raise ValueError(f"{label}: invalid numeric value")
    return int(value) if integer else value


def model_read_dataset(path, predictors, minimum=5, maximum=200000):
    """Validate complete joint count sets; incomplete/missing sets are audited."""
    groups, id_locations, location_labels = {}, {}, {}
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    needed = MODEL_REQUIRED | (set(predictors) - {"log_mutation_rate"})
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)) or not needed <= set(reader.fieldnames):
            raise ValueError("Input requires unique TSV headers including: " + ",".join(sorted(needed)))
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Line {line}: malformed TSV row")
            row = {k: v.strip() for k, v in row.items()}
            if row["build"].lower() not in {"38", "grch38", "hg38"}:
                raise ValueError(f"Line {line}: every row must declare GRCh38")
            chrom = re.sub(r"^chr", "", row["chrom"], flags=re.I).upper()
            if chrom == "M":
                chrom = "MT"
            if chrom not in {*(str(i) for i in range(1, 23)), "X", "Y", "MT"}:
                raise ValueError(f"Line {line}: noncanonical chromosome")
            pos = model_number(row["pos"], f"Line {line} POS", integer=True)
            ref, alt = row["ref"].upper(), row["alt"].upper()
            if pos is None or pos < 1 or ref not in MODEL_BASES or alt not in MODEL_BASES or ref == alt:
                raise ValueError(f"Line {line}: need a positive position and distinct canonical REF/ALT bases")
            for name in ("site_id", "population", "cohort", "ld_block"):
                if row[name].upper() in MODEL_MISSING:
                    raise ValueError(f"Line {line}: missing {name}")
            if row["callable_all_four"] not in {"0", "1"}:
                raise ValueError(f"Line {line}: callable_all_four must be explicit 0/1")
            sid, block = row["site_id"], row["ld_block"]
            location, label = (chrom, pos), (sid, ref, block)
            if sid in id_locations and id_locations[sid] != (chrom, pos, ref):
                raise ValueError(f"site_id {sid} identifies contradictory coordinates/REF")
            if location in location_labels and location_labels[location] != label:
                raise ValueError(f"{chrom}:{pos} has inconsistent site_id/REF/LD block across cohorts/populations")
            id_locations[sid], location_labels[location] = (chrom, pos, ref), label
            key = (sid, row["population"], row["cohort"])
            if key not in groups:
                if len(groups) >= maximum:
                    raise ValueError(f"--max-sites={maximum} exceeded; supply a bounded prespecified input")
                groups[key] = dict(site_id=sid, chrom=chrom, pos=pos, ref=ref, population=key[1], cohort=key[2], ld_block=block, rows={})
            if alt in groups[key]["rows"]:
                raise ValueError(f"Duplicate ALT {alt} in joint count set {key}")
            ac, an, ref_ac = (model_number(row[n], f"Line {line} {n}", integer=True) for n in ("AC", "AN", "ref_AC"))
            rate = model_number(row["mutation_rate"], f"Line {line} mutation_rate")
            if rate is not None and rate <= 0:
                raise ValueError(f"Line {line}: mutation_rate must be positive")
            values = []
            for name in predictors:
                value = ((None if rate is None else math.log(rate)) if name == "log_mutation_rate"
                         else rate if name == "mutation_rate" else model_number(row[name], f"Line {line} {name}"))
                values.append(value)
            for name, value in zip(predictors, values):
                if name in {"molecular_breadth", "effect_magnitude"} and value is not None and value < 0:
                    raise ValueError(f"Line {line}: {name} must be nonnegative")
            groups[key]["rows"][alt] = (ac, an, ref_ac, int(row["callable_all_four"]), values, rate)
    audit, designs, outcomes, clusters = [], [], [], []
    for group in groups.values():
        rows = [group["rows"][a] for a in sorted(group["rows"])]
        record = {k: v for k, v in group.items() if k != "rows"}
        record.update(n_alts=len(rows), k_supported="NA", eligibility="excluded", reason="")
        for index, name in ((1, "AN"), (2, "ref_AC"), (3, "callable_all_four")):
            if len({r[index] for r in rows if r[index] is not None}) > 1:
                raise ValueError(f"Joint set {group['site_id']}: inconsistent {name}")
        reason = "incomplete_three_alt_set" if set(group["rows"]) != MODEL_BASES - {group["ref"]} else ""
        if not reason and any(v is None for r in rows for v in r[:3]):
            reason = "missing_joint_counts"
        if not reason and (rows[0][1] < 1 or sum(r[0] for r in rows) + rows[0][2] != rows[0][1]):
            raise ValueError(f"Joint set {group['site_id']}: sum(three ALT AC)+ref_AC must equal positive shared AN")
        if not reason and not rows[0][3]:
            reason = "not_callable_all_four"
        if not reason and any(r[5] is None or any(v is None for v in r[4]) for r in rows):
            reason = "missing_predictor_or_mutation_rate"
        if reason:
            record["reason"] = reason
        else:
            y = [int(r[0] >= minimum) for r in rows]
            designs.append([r[4] for r in rows]); outcomes.append(y)
            clusters.append((group["chrom"], group["ld_block"]))
            k = sum(y)
            record.update(k_supported=k, eligibility="informative" if k in (1, 2) else "scaling_only",
                          reason="" if k in (1, 2) else "k0_no_supported_alt" if k == 0 else "k3_all_alts_supported")
        audit.append(record)
    x = np.asarray(designs, dtype=float).reshape((-1, 3, len(predictors)))
    y = np.asarray(outcomes, dtype=int).reshape((-1, 3))
    return x, y, clusters, audit


def model_choice_design(x, y):
    k = y.sum(axis=1)
    if not np.all(np.isin(y, (0, 1))) or not np.all(np.isin(k, (1, 2))):
        raise ValueError("Conditional likelihood requires binary y with k=1 or k=2")
    # For k=2 the omitted allele has probability softmax(-eta); no subset approximation.
    sign = np.where(k == 1, 1., -1.)
    chosen = np.where(k == 1, y.argmax(axis=1), y.argmin(axis=1))
    return x * sign[:, None, None], chosen


def model_objective(beta, z, chosen):
    eta = z @ beta
    log_normalizer = special.logsumexp(eta, axis=1)
    probs = np.exp(eta - log_normalizer[:, None])
    mean = np.einsum("ni,nip->np", probs, z)
    score = z[np.arange(len(z)), chosen] - mean
    value = float(np.sum(log_normalizer - eta[np.arange(len(z)), chosen]))
    hessian = np.einsum("ni,nip,niq->pq", probs, z, z) - mean.T @ mean
    return value, -score.sum(axis=0), hessian, score


def model_separation_check(z, chosen):
    differences = z[np.arange(len(z)), chosen, None, :] - z
    mask = np.ones(z.shape[:2], dtype=bool); mask[np.arange(len(z)), chosen] = False
    contrasts = differences[mask]
    # A bounded direction with all margins >=0 and a positive average proves
    # complete or quasi separation. Zero is always feasible in this LP.
    lp = optimize.linprog(-contrasts.mean(axis=0), A_ub=-contrasts, b_ub=np.zeros(len(contrasts)),
                          bounds=[(-1, 1)] * z.shape[2], method="highs")
    if not lp.success:
        return "separation_check_failed", dict(message=lp.message)
    margins = contrasts @ lp.x
    separated = margins.min() >= -1e-7 and margins.mean() > 1e-7
    reason = ("complete_separation" if margins.min() > 1e-7 else "quasi_separation") if separated else "none"
    return reason, dict(min_margin=float(margins.min()), mean_margin=float(margins.mean()), tolerance=1e-7)


def model_bh(values):
    values = np.asarray(values); order = np.argsort(values)
    adjusted = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty(len(values)); result[order] = np.clip(adjusted, 0, 1)
    return result


def model_fit(x, y, clusters, predictors):
    n, p = len(x), len(predictors)
    rows = [dict.fromkeys(MODEL_COEFFICIENT_FIELDS, "NA") | dict(predictor=name) for name in predictors]
    diagnostics = dict(complete_callable_sets=n, predictors=p, inference_available=False)
    def unavailable(reason, estimated=False):
        diagnostics.update(reason=reason, status="complete_with_exclusions" if estimated else "not_estimable")
        for row in rows: row["inference_status"] = reason
        return rows, diagnostics
    if not n:
        return unavailable("no_complete_callable_sets")
    with np.errstate(over="ignore", invalid="ignore"):
        centered = x - x.mean(axis=1, keepdims=True)
        sd = np.sqrt(np.mean(centered ** 2, axis=(0, 1)))
    for row, value in zip(rows, sd): row["within_site_sd"] = float(value) if math.isfinite(value) else "NA"
    if np.any(~np.isfinite(sd)) or np.any(sd <= 1e-12):
        return unavailable("zero_or_invalid_within_site_predictor_variance")
    eligible = np.isin(y.sum(axis=1), (1, 2))
    diagnostics.update(informative_sets=int(eligible.sum()), scaling_sets=n,
                       k0_sets=int((y.sum(axis=1) == 0).sum()), k3_sets=int((y.sum(axis=1) == 3).sum()))
    if eligible.sum() <= p:
        return unavailable("insufficient_informative_sets")
    z, chosen = model_choice_design(centered[eligible] / sd[None, None, :], y[eligible])
    block_ids = list(itertools.compress(clusters, eligible))
    block_map = {key: i for i, key in enumerate(sorted(set(block_ids)))}
    cluster_index = np.array([block_map[key] for key in block_ids])
    groups, n = len(block_map), len(z)
    diagnostics.update(genomic_blocks=groups, minimum_inference_blocks=20)
    separated, check = model_separation_check(z, chosen)
    diagnostics["separation_check"] = check
    if separated != "none":
        return unavailable(separated)
    def loss(beta):
        value, gradient, _, _ = model_objective(beta, z, chosen)
        return value / n, gradient / n
    result = optimize.minimize(loss, np.zeros(p), jac=True, method="BFGS", options=dict(maxiter=1000, gtol=1e-8))
    value, gradient, hessian, scores = model_objective(result.x, z, chosen)
    diagnostics.update(optimizer_success=bool(result.success), optimizer_message=str(result.message),
                       optimizer_iterations=int(result.nit), mean_gradient_max=float(np.max(np.abs(gradient / n))),
                       negative_log_likelihood=float(value))
    if not result.success or not np.isfinite(result.x).all() or np.max(np.abs(gradient / n)) > 1e-6:
        return unavailable("optimizer_not_converged")
    eigen = np.linalg.eigvalsh(hessian)
    diagnostics["hessian_eigenvalues"] = eigen.tolist()
    if not np.isfinite(eigen).all() or eigen[0] <= max(1e-10, eigen[-1] * 1e-10):
        return unavailable("near_singular_hessian")
    for row, beta in zip(rows, result.x):
        row.update(coefficient=float(beta), odds_ratio_per_within_site_sd=float(np.exp(beta)) if beta < 700 else "Inf")
    if groups < 20 or groups <= p:
        return unavailable("insufficient_genomic_blocks_for_inference", estimated=True)
    cluster_scores = np.zeros((groups, p)); np.add.at(cluster_scores, cluster_index, scores)
    inverse = np.linalg.inv(hessian)
    correction = groups / (groups - 1) * (n - 1) / (n - p)
    covariance = correction * inverse @ (cluster_scores.T @ cluster_scores) @ inverse
    diagonal = np.diag(covariance)
    if not np.isfinite(diagonal).all() or np.any(diagonal <= 0):
        return unavailable("invalid_cluster_sandwich_variance", estimated=True)
    se = np.sqrt(diagonal); quantile = stats.t.ppf(.975, groups - 1)
    pvalues = 2 * stats.t.sf(np.abs(result.x / se), groups - 1)
    for row, beta, error, pv, qv in zip(rows, result.x, se, pvalues, model_bh(pvalues)):
        low, high = beta - quantile * error, beta + quantile * error
        row.update(cluster_robust_se=float(error), ci95_low=float(np.exp(low)) if low > -745 else 0.,
                   ci95_high=float(np.exp(high)) if high < 700 else "Inf", p_value=float(pv),
                   q_value_within_model=float(qv), inference_status="cluster_robust_association")
    diagnostics.update(status="complete", reason="cluster_robust_association", inference_available=True,
                       sandwich_correction=float(correction), t_degrees_freedom=groups - 1,
                       hessian_condition_number=float(eigen[-1] / eigen[0]))
    return rows, diagnostics


def model_write_table(path, rows, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


@transactional("manifest.json", input_files=("table",), input_prefixes=())
def model_run(args):
    if args.sampling_frame == "ascertained_multialt":
        raise ValueError("The primary model refuses ascertained_multialt inputs; provide all callable or outcome-independent prespecified sites")
    if args.count_evidence != "joint_wgs_genotypes" or args.predictor_evidence not in {"raw_molecular_predictions", "experimental_molecular_effects"}:
        raise ValueError("Explicit joint-WGS count and raw/experimental predictor evidence declarations are required")
    if args.max_sites < 1 or args.min_ac not in (5, 10) or args.build not in {"GRCh38", "38", "hg38"}:
        raise ValueError("Require GRCh38, --min-ac 5 or 10, and positive --max-sites")
    predictors = model_predictor_names(args.predictors)
    x, y, clusters, audit = model_read_dataset(args.table, predictors, args.min_ac, args.max_sites)
    coefficients, diagnostics = model_fit(x, y, clusters, predictors)
    out = Path(args.out)
    model_write_table(out / "coefficients.tsv", coefficients, MODEL_COEFFICIENT_FIELDS)
    audit_fields = ["site_id", "chrom", "pos", "ref", "population", "cohort", "ld_block", "n_alts", "k_supported", "eligibility", "reason"]
    model_write_table(out / "site_audit.tsv.gz", audit, audit_fields)
    counts = dict(Counter(r["reason"] or "informative" for r in audit))
    (out / "diagnostics.json").write_text(json.dumps(diagnostics, indent=2, allow_nan=False) + "\n")
    return dict(status=diagnostics["status"], schema=MODEL_VERSION, build="GRCh38", diagnostics=diagnostics,
                exclusions_and_informative_sets=counts, source=identity(args.table, full=True),
                parameters={k: getattr(args, k) for k in ("min_ac", "predictors", "sampling_frame", "count_evidence", "predictor_evidence", "max_sites")},
                scaling="Population SD of within-set centered predictors over all complete callable three-ALT sets before conditioning on k; one SD per coefficient",
                primary_scope="genome_wide_callable_frame" if args.sampling_frame == "all_callable_sites" else "prespecified_target_conditional_association",
                source_evidence_verification="caller declarations and numeric/schema consistency; underlying WGS/predictor provenance not independently verified",
                limitations=MODEL_LIMITATIONS)


def model_parser():
    p = argparse.ArgumentParser(description=_description_opportunity_model)
    p.add_argument("--table", required=True, help="Three rows/site/population/cohort, one for each possible nonreference base, including true AC=0")
    p.add_argument("--build", required=True, choices=("GRCh38", "38", "hg38"))
    p.add_argument("--min-ac", type=int, default=5, choices=(5, 10))
    p.add_argument("--predictors", default=MODEL_DEFAULT_PREDICTORS)
    p.add_argument("--sampling-frame", required=True, choices=("all_callable_sites", "prespecified_sites", "ascertained_multialt"))
    p.add_argument("--count-evidence", required=True, choices=("joint_wgs_genotypes",))
    p.add_argument("--predictor-evidence", required=True, choices=("raw_molecular_predictions", "experimental_molecular_effects"))
    p.add_argument("--max-sites", type=int, default=200000, help="Maximum joint site/population/cohort sets in this prespecified analysis")
    p.add_argument("--out", required=True)
    return p


def model_main(argv=None):
    try:
        result = model_run(model_parser().parse_args(argv))
    except (ValueError, OSError, RuntimeError) as exc:
        raise SystemExit("opportunity-model: " + str(exc))
    print(json.dumps(dict(status=result["status"], diagnostics=result["diagnostics"]), allow_nan=False))
    return 0 if result["status"] != "not_estimable" else 2


def cli(argv=None):
    """Dispatch auxiliary commands owned by this module."""
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    routes = {'model': model_main}
    entry = routes[argv.pop(0)] if argv and argv[0] in routes else main
    return entry(argv)



if __name__ == "__main__":
    raise SystemExit(cli())

