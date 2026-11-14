#!/usr/bin/env python3
"""Continuous multi-ALT density scores from completed observed-variant counters.

This is an SDS-inspired density analysis, NOT Singleton Density Score (SDS),
a recent-selection test, or a population count of three/four observed bases.
Legacy strata establish supported nonreference ALT counts but not REF presence.
No raw VCF, AlphaGenome prediction catalogue, or selected hotspot list is read.

Two compact-counter passes fit a joint category model and score all windows.
The target window uses leave-one-chromosome-out (LOCO) expectations. Its
standardization reference windows use a stricter double exclusion: both the
target chromosome and the reference window's chromosome are removed from the
background. Thus target-chromosome counts cannot affect its reference mean/SD.
Empirical Z is a descriptive standardized score; no normal-tail P is assigned.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import html
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

from analysis_integrity import identity, transactional
import genome_windows as gw

VERSION = "allele-density-1.0"
CHANNELS = ("multi_alt_ge2", "three_alt_eq3")
CHANNEL_LABELS = {
    "multi_alt_ge2": "Sites with at least two AC-supported nonreference ALT bases",
    "three_alt_eq3": "Sites with three AC-supported nonreference ALT bases",
}
REFERENCES = [
    {"title": "Field et al. (2016), Detection of human adaptation during the past 2000 years",
     "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5182071/"},
    {"title": "Original Singleton Density Score implementation",
     "url": "https://github.com/yairf/SDS"},
]
LIMITATIONS = [
    "This is a newly defined multi-ALT density score, not the 2016 Singleton Density Score. It does not infer terminal branch lengths, recent allele-frequency change, or a selection time scale.",
    "The denominator is threshold-specific observed passing variant loci with at least one nonreference ALT AC>=threshold; it is not all genomic bases, callable bases, people, or minor-allele carriers.",
    "raw_loci_per_mb uses the recorded window span in bases; it is a descriptive span-denominator density and is not corrected for callable-base coverage.",
    "multi_alt_ge2 includes sites with two or three supported nonreference ALTs. three_alt_eq3 counts three supported nonreference ALTs. Without jointly observed REF counts, neither channel proves respectively three or four observed population alleles.",
    "The legacy counters do not necessarily establish a common genotype sample set across decomposed ALT records; joint AC/AN/genotype validation is still needed for true all-base claims.",
    "Exact context/change/primary-AF/AN matching controls recorded strata only. Recurrent mutation, CpG effects, ancestry, sample size, mapping, callability, QC and variant ascertainment can still affect scores.",
    "The joint Dirichlet prior is explicit shrinkage. Aggregate binomial deviance and marginal predictive variance are working descriptions; linked sites are not independent Bernoulli observations.",
    "Double-LOCO reference fitting excludes the target chromosome entirely, but reference windows can still share training data and long-range genomic effects. Empirical Z is not a validated standard-normal statistic; no P value, q value, significance or selection claim is produced.",
    "Calibration includes every eligible nonoverlapping reference window, including zero observed target counts and enriched windows. Its distribution can contain biological and technical heterogeneity; it is not an independently established neutral null.",
    "All sizes and AC tiers must be declared before interpretation. AC5 and AC10 are sensitivity analyses of the same source, not independent replication.",
    "Only aligned unions of complete native counter bins are possible. Counts cannot be subdivided, interpolated or retrospectively masked into finer windows. Raw source/window counts remain unchanged.",
    "A prespecified exclusion BED removes each overlapping WHOLE native bin from background training and vetoes overlapping target/reference windows. No count is prorated; this is a conservative whole-bin exclusion, not a reconstructed base-level strict-mask scan.",
    "The declared source scope defines the population of genomic windows in this analysis. Explicit caller-supplied counters do not independently certify whole-genome coverage or freedom from prior region selection.",
    "Higher or lower density does not by itself establish tolerance, pathogenicity, health benefit, positive selection or balancing selection. Those hypotheses require independently justified functional, ancient and phenotype analyses.",
]
FIELDS = [
    "window_id", "chrom", "start0", "end0", "build", "window_bp", "span_bp", "step_bp",
    "threshold", "is_primary", "channel", "denominator", "raw_n_eligible", "raw_observed",
    "raw_observed_rate", "raw_loci_per_mb", "raw_exactly_two_alt_sites", "raw_three_alt_sites",
    "supported_n_eligible", "supported_observed", "expected", "observed_rate", "expected_rate",
    "oe", "full_scope_oe", "log2_oe_pseudocount", "unsupported_n_eligible", "unsupported_observed",
    "background_support_fraction", "spatial_support_fraction", "source_bins_present", "source_bins_required",
    "strata_count", "training_min_sites", "training_min_chromosomes", "predictive_variance",
    "mds_rawdeviance", "mds_pearson", "mds_z", "score_eligibility", "standardization_status",
    "reference_windows", "reference_chromosomes", "reference_deviance_mean", "reference_deviance_sd",
    "is_nonoverlapping_reference_grid", "mask_overlap_bp", "observed_exactly_three_allele_sites",
    "observed_exactly_four_allele_sites", "reference_presence_status", "expectation_model",
]


class JointRates:
    """Joint [one ALT, exactly two ALT, three ALT] model with LOCO vectors."""

    def __init__(self, minimum):
        self.minimum = minimum
        self.by_stratum = defaultdict(lambda: np.zeros((22, 3), dtype=np.int64))
        self._sealed = False

    def observe(self, row):
        if self._sealed:
            raise ValueError("Background cannot change after prediction begins")
        chrom_index = int(gw.chromosome(row["chrom"])) - 1
        self.by_stratum[gw.counter_key(row)][chrom_index] += (
            row["n_sites"], row["n_multi"], row["n_three"])

    @lru_cache(maxsize=8192)
    def predict(self, reference_chrom, key):
        """Return 22 backgrounds, one for each target chromosome.

        At target==reference this is ordinary LOCO. Else exclude BOTH target
        and reference chromosomes. Vectorization avoids 22 source-file passes.
        At least two remaining chromosomes with eligible loci support a stratum.
        """
        self._sealed = True
        counts = self.by_stratum.get(key)
        if counts is None:
            return None
        reference = int(gw.chromosome(reference_chrom)) - 1
        totals = counts.sum(axis=0)
        training = totals[np.newaxis, :] - counts[reference] - counts
        training[reference] += counts[reference]
        present = counts[:, 0] > 0
        training_chromosomes = present.sum() - int(present[reference]) - present.astype(int)
        training_chromosomes[reference] += int(present[reference])
        training_n = training[:, 0].astype(float)
        valid = (training_n >= self.minimum) & (training_chromosomes >= 2)
        # Dirichlet(1/2,1/2,1/2) over [N-M, M-Q, Q], preserving pQ<=pM.
        concentration = training_n + 1.5
        probabilities = np.column_stack((training[:, 1] + 1., training[:, 2] + .5)) / concentration[:, np.newaxis]
        return probabilities, concentration, valid, training_n, training_chromosomes


def signed_binomial_deviance(observed, expected, denominator):
    """Aggregate-mean binomial signed deviance; not the full stratum likelihood."""
    o, e, n = np.broadcast_arrays(np.asarray(observed, dtype=float), np.asarray(expected, dtype=float),
                                  np.asarray(denominator, dtype=float))
    result = np.full(o.shape, np.nan, dtype=float)
    valid = (n > 0) & (e > 0) & (e < n) & (o >= 0) & (o <= n)
    if np.any(valid):
        oo, ee, nn = o[valid], e[valid], n[valid]
        terms = np.zeros(oo.shape)
        positive = oo > 0
        terms[positive] += oo[positive] * np.log(oo[positive] / ee[positive])
        complement = oo < nn
        terms[complement] += ((nn[complement] - oo[complement]) *
                              np.log((nn[complement] - oo[complement]) / (nn[complement] - ee[complement])))
        result[valid] = np.sign(oo-ee) * np.sqrt(np.maximum(0., 2.*terms))
    return result


class Moments:
    """Stable streaming mean and sample variance, with chromosome provenance."""

    def __init__(self):
        self.n = 0
        self.mean = 0.
        self.m2 = 0.
        self.chromosomes = set()

    def observe(self, value, chrom):
        if not math.isfinite(value):
            raise ValueError("Nonfinite reference deviance")
        self.n += 1
        delta = value - self.mean
        self.mean += delta / self.n
        self.m2 += delta * (value-self.mean)
        self.chromosomes.add(chrom)

    @property
    def sd(self):
        return math.sqrt(max(0., self.m2/(self.n-1))) if self.n > 1 else None


def summarize_counts(counts, chrom, rates):
    """Scores and support for one window under ordinary/double LOCO fits."""
    n = np.zeros(22, dtype=np.int64)
    observed = np.zeros((22, 2), dtype=np.int64)
    expected = np.zeros((22, 2))
    variance = np.zeros((22, 2))
    min_training_n = np.full(22, np.inf)
    min_training_chromosomes = np.full(22, np.inf)
    raw = np.zeros(3, dtype=np.int64)
    for key, values in counts.items():
        nn, mm, qq = values
        raw += values
        prediction = rates.predict(chrom, key)
        if prediction is None or nn == 0:
            continue
        probability, concentration, valid, training_n, training_chromosomes = prediction
        n[valid] += nn
        observed[valid] += [mm, qq]
        expected[valid] += nn * probability[valid]
        # Marginal Dirichlet-multinomial/Beta-binomial variance. Counts for the
        # SAME stratum are combined across native bins before the n*(n-1) term.
        variance[valid] += (nn * probability[valid] * (1.-probability[valid]) *
                            (1. + (nn-1.)/(concentration[valid, np.newaxis]+1.)))
        min_training_n[valid] = np.minimum(min_training_n[valid], training_n[valid])
        min_training_chromosomes[valid] = np.minimum(min_training_chromosomes[valid], training_chromosomes[valid])
    deviance = signed_binomial_deviance(observed, expected, n[:, np.newaxis])
    pearson = np.full((22, 2), np.nan)
    np.divide(observed-expected, np.sqrt(variance), out=pearson, where=variance > 0)
    return dict(raw=raw, n=n, observed=observed, expected=expected, variance=variance,
                deviance=deviance, pearson=pearson, min_training_n=min_training_n,
                min_training_chromosomes=min_training_chromosomes)


def eligibility(window, values, target, channel, args):
    if window["source_bins_present"] != window["source_bins_required"]:
        return "incomplete_source_bin_coverage"
    if window["end0"]-window["start0"] != window["window_bp"]:
        return "truncated_scope_edge"
    if window["mask_overlap_bp"]:
        return "overlaps_exclusion_bed"
    if values["n"][target] < args.min_sites:
        return "too_few_supported_eligible_sites"
    raw_n = int(values["raw"][0])
    if not raw_n or values["n"][target]/raw_n < args.min_support:
        return "insufficient_exact_background_support"
    if values["expected"][target, channel] < args.min_expected:
        return "too_little_expected_support"
    if not math.isfinite(float(values["deviance"][target, channel])) or not math.isfinite(float(values["pearson"][target, channel])):
        return "nonfinite_working_score"
    return "eligible"


def finite(value):
    return float(value) if math.isfinite(float(value)) else None


def window_rows(window, counts, rates, args, reference_moments, source_chromosomes):
    values = summarize_counts(counts, window["chrom"], rates)
    source_chrom = gw.chromosome(window["chrom"])
    own = int(source_chrom)-1
    raw_n, raw_multi, raw_three = map(int, values["raw"])
    grid = window["start0"] % window["window_bp"] == 0
    # A target's reference residual cannot use its own chromosome anywhere:
    # predict() excludes target+reference from the fitted rate for this slot.
    if grid:
        for target_chrom in source_chromosomes:
            if target_chrom == source_chrom:
                continue
            target = int(target_chrom)-1
            for channel, name in enumerate(CHANNELS):
                if eligibility(window, values, target, channel, args) == "eligible":
                    reference_moments[(target_chrom, window["window_bp"], window["threshold"], name)].observe(
                        float(values["deviance"][target, channel]), source_chrom)
    result = []
    supported_n = int(values["n"][own])
    for channel, name in enumerate(CHANNELS):
        raw_observed = (raw_multi, raw_three)[channel]
        observed = int(values["observed"][own, channel])
        expected = float(values["expected"][own, channel]) if supported_n else None
        supported = supported_n/raw_n if raw_n else None
        status = eligibility(window, values, own, channel, args)
        result.append(dict(
            window_id=window["window_id"], chrom=window["chrom"], start0=window["start0"], end0=window["end0"],
            build="GRCh38", window_bp=window["window_bp"], span_bp=window["end0"]-window["start0"], step_bp=args.step_bp,
            threshold=window["threshold"], is_primary=int(window["threshold"] == args.primary_threshold),
            channel=name, denominator=gw.DENOMINATOR, raw_n_eligible=raw_n, raw_observed=raw_observed,
            raw_observed_rate=raw_observed/raw_n if raw_n else None,
            raw_loci_per_mb=raw_observed*1e6/(window["end0"]-window["start0"]),
            raw_exactly_two_alt_sites=raw_multi-raw_three, raw_three_alt_sites=raw_three,
            supported_n_eligible=supported_n, supported_observed=observed,
            expected=expected, observed_rate=observed/supported_n if supported_n else None,
            expected_rate=expected/supported_n if supported_n else None,
            oe=observed/expected if expected else None,
            full_scope_oe=raw_observed/expected if expected and supported_n == raw_n and window["source_bins_present"] == window["source_bins_required"] else None,
            log2_oe_pseudocount=math.log2((observed+.5)/(expected+.5)) if expected is not None else None,
            unsupported_n_eligible=raw_n-supported_n, unsupported_observed=raw_observed-observed,
            background_support_fraction=supported,
            spatial_support_fraction=window["source_bins_present"]/window["source_bins_required"],
            source_bins_present=window["source_bins_present"], source_bins_required=window["source_bins_required"],
            strata_count=len(counts), training_min_sites=finite(values["min_training_n"][own]),
            training_min_chromosomes=finite(values["min_training_chromosomes"][own]),
            predictive_variance=float(values["variance"][own, channel]) if supported_n else None,
            mds_rawdeviance=float(values["deviance"][own, channel]) if status == "eligible" else None,
            mds_pearson=float(values["pearson"][own, channel]) if status == "eligible" else None,
            mds_z=None, score_eligibility=status, is_nonoverlapping_reference_grid=int(grid),
            mask_overlap_bp=window["mask_overlap_bp"],
            observed_exactly_three_allele_sites=None, observed_exactly_four_allele_sites=None,
            reference_presence_status="unknown_from_legacy_strata", expectation_model="joint_Dirichlet_half_exact_stratum_LOCO",
        ))
    return result


def collect_windows(files, native, scopes, args, masks):
    rates = JointRates(args.min_training_sites)
    audit = Counter(training_excluded_native_bins=0, training_excluded_strata_rows=0, training_used_strata_rows=0)
    for threshold in args.threshold_values:
        for label in ("n_sites", "n_multi", "n_three"):
            audit[f"AC{threshold}:training_excluded_{label}"] = 0
    clock = time.monotonic()
    mask_index = gw.IntervalOverlapIndex(masks)
    last_excluded_bin = None
    for row in gw.read_counters(files, native, args.threshold_values, scopes):
        audit["source_strata_rows"] += 1
        for label in ("n_sites", "n_multi", "n_three"):
            audit[f"AC{row['threshold']}:{label}"] += row[label]
        if gw.overlap_bp(mask_index, row["chrom"], row["start0"], row["end0"]):
            current_bin = (row["chrom"], row["start0"], row["end0"])
            if current_bin != last_excluded_bin:
                audit["training_excluded_native_bins"] += 1
                last_excluded_bin = current_bin
            audit["training_excluded_strata_rows"] += 1
            for label in ("n_sites", "n_multi", "n_three"):
                audit[f"AC{row['threshold']}:training_excluded_{label}"] += row[label]
        else:
            rates.observe(row)
            audit["training_used_strata_rows"] += 1
        if time.monotonic()-clock > 50:
            print(f"[allele-density] background strata={audit['source_strata_rows']:,}", flush=True)
            clock = time.monotonic()
    print(f"[allele-density] background complete; compact strata={audit['source_strata_rows']:,}", flush=True)
    groups = iter(itertools.groupby(gw.read_counters(files, native, args.threshold_values, scopes),
                                   key=lambda row: (row["chrom"], row["start0"], row["end0"])))
    pending = next(groups, None)
    rows = []
    moments = defaultdict(Moments)
    chromosomes = sorted(scopes, key=int)
    physical_windows = 0
    for chrom in chromosomes:
        for scope_start, scope_end in scopes[chrom]:
            if scope_start % native:
                raise ValueError("Scan scopes must start on native-bin boundaries; source counts cannot be clipped")
            active = []

            def flush(end):
                nonlocal active
                retained = []
                for window, counts in active:
                    if window["end0"] <= end:
                        rows.extend(window_rows(window, counts, rates, args, moments, chromosomes))
                    else:
                        retained.append((window, counts))
                active = retained

            for lo in range(scope_start, scope_end, native):
                flush(lo)
                hi = min(lo+native, scope_end)
                if lo % args.step_bp == 0:
                    for width in args.size_values:
                        end = min(lo+width, scope_end)
                        for threshold in args.threshold_values:
                            window = dict(window_id=f"mds_chr{chrom}_{lo}_{end}_w{width}_ac{threshold}",
                                          chrom="chr"+chrom, start0=lo, end0=end, window_bp=width, threshold=threshold,
                                          source_bins_present=0, source_bins_required=math.ceil((end-lo)/native),
                                          mask_overlap_bp=gw.overlap_bp(mask_index, chrom, lo, end))
                            active.append((window, {}))
                            physical_windows += 1
                            if physical_windows > args.max_windows:
                                raise ValueError("--max-windows exceeded; predeclare fewer widths or increase the resource limit")
                for window, _ in active:
                    window["source_bins_present"] += 1
                if pending is not None and pending[0] == (chrom, lo, hi):
                    tile_totals = defaultdict(lambda: [0, 0, 0])
                    for row in pending[1]:
                        values = [row[name] for name in ("n_sites", "n_multi", "n_three")]
                        tile_totals[row["threshold"]] = [a+b for a, b in zip(tile_totals[row["threshold"]], values)]
                        key = gw.counter_key(row)
                        for window, counts in active:
                            if window["threshold"] == row["threshold"]:
                                existing = counts.setdefault(key, [0, 0, 0])
                                for j in range(3):
                                    existing[j] += values[j]
                    previous = None
                    for threshold in args.threshold_values:
                        values = tile_totals[threshold]
                        if values[0] > hi-lo:
                            raise ValueError(f"chr{chrom}:{lo}-{hi}: total eligible loci exceed source-bin bases")
                        if previous is not None and any(a > b for a, b in zip(values, previous)):
                            raise ValueError(f"chr{chrom}:{lo}-{hi}: counts increase at stricter AC threshold")
                        previous = values
                    pending = next(groups, None)
                elif pending is not None and (int(pending[0][0]), pending[0][1]) < (int(chrom), lo):
                    raise ValueError("Source-bin order disagrees with declared scope")
                if time.monotonic()-clock > 50:
                    print(f"[allele-density] scored channel rows={len(rows):,}; at=chr{chrom}:{hi}", flush=True)
                    clock = time.monotonic()
            flush(scope_end)
    if pending is not None:
        raise ValueError("Unconsumed source bins outside declared scope")
    audit["physical_window_tiers"] = physical_windows
    audit["density_channel_rows"] = len(rows)
    rates.predict.cache_clear()
    return rows, moments, audit


def standardize(rows, moments, args):
    calibrations = {}
    for row in rows:
        key = (gw.chromosome(row["chrom"]), row["window_bp"], row["threshold"], row["channel"])
        if key not in calibrations:
            reference = moments.get(key, Moments())
            if reference.n < args.min_calibration_windows or len(reference.chromosomes) < 2:
                status = "insufficient_double_LOCO_reference_windows"
            elif (reference.sd is None or not math.isfinite(reference.sd) or
                  reference.sd <= 1e-12*max(1., abs(reference.mean))):
                status = "zero_or_negligible_reference_variance"
            else:
                status = "empirically_standardized_not_a_normal_test"
            calibrations[key] = dict(
                chrom=row["chrom"], window_bp=row["window_bp"], threshold=row["threshold"], channel=row["channel"],
                status=status, reference_windows=reference.n, reference_chromosomes=len(reference.chromosomes),
                reference_chromosome_list=",".join("chr"+c for c in sorted(reference.chromosomes, key=int)),
                minimum_reference_windows=args.min_calibration_windows,
                reference_deviance_mean=reference.mean if reference.n else None, reference_deviance_sd=reference.sd,
                background_exclusion="target_chromosome_and_reference_window_chromosome",
                eligibility="N_E_support_nonoverlap_scope_and_external_mask_only_no_observed_target_count_gate",
                inference="descriptive_empirical_Z_no_P_or_q_value",
            )
        record = calibrations[key]
        for field in ("reference_windows", "reference_chromosomes", "reference_deviance_mean", "reference_deviance_sd"):
            row[field] = record[field]
        row["standardization_status"] = (record["status"] if row["score_eligibility"] == "eligible"
                                         else "ineligible:"+row["score_eligibility"])
        if row["score_eligibility"] == "eligible" and record["status"] == "empirically_standardized_not_a_normal_test":
            row["mds_z"] = (row["mds_rawdeviance"]-record["reference_deviance_mean"])/record["reference_deviance_sd"]
    return list(calibrations.values())


def validate_args(args):
    if args.build not in {"GRCh38", "hg38", "38"}:
        raise ValueError("Explicit GRCh38 build declaration is required")
    args.size_values = gw.ints(args.window_sizes, "--window-sizes")
    args.threshold_values = gw.ints(args.thresholds, "--thresholds")
    if not set(args.threshold_values) <= {5, 10} or args.primary_threshold not in args.threshold_values:
        raise ValueError("Use recorded AC5/AC10 thresholds and include the primary threshold")
    if min(args.tile_bp, args.step_bp, args.min_sites, args.min_training_sites, args.max_windows) < 1:
        raise ValueError("Widths, count/support minimums and resource limits must be positive")
    if not math.isfinite(args.min_expected) or args.min_expected <= 0 or not 0 < args.min_support <= 1:
        raise ValueError("Require finite --min-expected>0 and 0<--min-support<=1")
    if args.min_calibration_windows < 20:
        raise ValueError("Require at least 20 unselected reference windows; default is 50")
    if args.tiles and not args.scope:
        raise ValueError("--tiles requires the complete original --scope; absent bins are otherwise unknown")
    if not args.tiles and args.scope:
        raise ValueError("--scope applies only to explicit --tiles; store uses validated recorded scopes")


def load_source(args):
    sources = {}
    if args.store:
        import scan
        record = scan.load_run(args.store)
        options = record["options"]
        if not set(args.threshold_values) <= set(map(int, options.get("thresholds", []))):
            raise ValueError("Requested AC tiers were not included in the completed scan")
        build = str(options.get("architecture_build", options.get("build", "38"))).lower()
        if build not in {"38", "hg38", "grch38"}:
            raise ValueError("Completed source scan does not declare GRCh38")
        files = [Path(shard)/"strata.tsv.gz" for shard in record["shards"]]
        scopes = [Path(shard)/"scan_scope.tsv" for shard in record["shards"]]
        native = gw.integer(options.get("window_size", args.tile_bp), "source window_size", 1)
        metadata = dict(original_scan_signature=record["signature"], source_options=options,
                        ref_observation="not_recorded_by_legacy_strata", source_validation="completed_scan_store")
        sources["scan_manifest"] = identity(Path(args.store)/"run.json", full=True)
    else:
        files = [Path(path).resolve() for path in args.tiles]
        scopes = [Path(args.scope).resolve()]
        native = args.tile_bp
        metadata = dict(original_scan_signature="caller_supplied_counters_not_independently_verified",
                        ref_observation="not_recorded_by_legacy_strata", source_validation="explicit_schema_and_scope_only")
    if not files or len(files) > 128 or len(set(files)) != len(files):
        raise ValueError("Require 1..128 unique compact counter files")
    for path in (*files, *scopes):
        if not path.is_file():
            raise ValueError(f"Completed counter/scope file missing: {path}; no raw VCF rescan is started")
    gw.validate_grid(native, args)
    scope_intervals = gw.intervals(scopes, allow_overlap=False)
    if len(scope_intervals) < 3:
        raise ValueError("At least three autosomal source chromosomes are needed for LOCO scores; Z needs at least four with support")
    sources["counters"] = [identity(path) for path in files]
    sources["scopes"] = [identity(path, full=True) for path in scopes]
    return files, native, scope_intervals, metadata, sources


def plot_scores(out, rows):
    """Static scientific figures: signed density scores, with no P threshold."""
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    widths = sorted({row["window_bp"] for row in rows})
    thresholds = sorted({row["threshold"] for row in rows})
    chromosomes = sorted({gw.chromosome(row["chrom"]) for row in rows}, key=int)
    spans = {chrom: max(row["end0"] for row in rows if gw.chromosome(row["chrom"]) == chrom) for chrom in chromosomes}
    offsets, running = {}, 0
    for chrom in chromosomes:
        offsets[chrom] = running
        running += spans[chrom]
    figures = []
    palette = ("#237e8a", "#7465a0")
    with plt.rc_context({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
                         "svg.fonttype": "none", "figure.facecolor": "white", "savefig.facecolor": "white"}):
        for width_index, width in enumerate(widths):
            figure, axes = plt.subplots(2, len(thresholds), figsize=(max(8, 6.4*len(thresholds)), 7),
                                        sharex=True, squeeze=False, constrained_layout=True)
            for channel_index, channel in enumerate(CHANNELS):
                for threshold_index, threshold in enumerate(thresholds):
                    ax = axes[channel_index, threshold_index]
                    panel = [row for row in rows if row["window_bp"] == width and row["threshold"] == threshold and row["channel"] == channel]
                    estimable = [row for row in panel if row["mds_z"] is not None]
                    for chrom_index, chrom in enumerate(chromosomes):
                        points = [row for row in estimable if gw.chromosome(row["chrom"]) == chrom]
                        ax.scatter([offsets[chrom]+(row["start0"]+row["end0"])/2 for row in points],
                                   [row["mds_z"] for row in points], s=10, alpha=.75,
                                   color=palette[chrom_index % 2], linewidths=0, rasterized=len(points) > 10000)
                    ax.axhline(0, color="#566673", linewidth=.7, linestyle="--", zorder=0)
                    ax.margins(y=.16)
                    ax.set_title(("At least two supported ALT bases" if channel_index == 0 else "Three supported ALT bases") + f" · AC ≥ {threshold}")
                    ax.set_ylabel("Empirical density Z\n(descriptive)")
                    ax.text(.99, .97, f"Z available: {len(estimable):,} / {len(panel):,}", transform=ax.transAxes,
                            ha="right", va="top", fontsize=9, color="#536273")
                    if not estimable:
                        ax.text(.5, .5, "Insufficient reference/support for Z", transform=ax.transAxes,
                                ha="center", va="center", fontsize=10, color="#536273")
                    ax.spines[["top", "right"]].set_visible(False)
                    ax.set_xticks([offsets[c]+spans[c]/2 for c in chromosomes], chromosomes)
                    if channel_index == 1:
                        ax.set_xlabel("Autosomal chromosome")
            figure.suptitle(f"Multi-ALT density · {width:,}-bp windows\nDouble-LOCO standardization; positive Z is above the reference mean", fontsize=14)
            stem = "density_scores" if width_index == 0 else f"density_scores_w{width}"
            figure.savefig(out/(stem+".png"), dpi=180)
            figure.savefig(out/(stem+".svg"))
            plt.close(figure)
            figures.append(dict(window_bp=width, png=stem+".png", svg=stem+".svg",
                                omitted_unavailable_z_rows=sum(row["window_bp"] == width and row["mds_z"] is None for row in rows)))
    return figures


def write_report(path, rows, audit, figures):
    """A small self-contained review page, with every window retained in the TSV."""
    sections = []
    for channel in CHANNELS:
        eligible = [row for row in rows if row["channel"] == channel and row["mds_z"] is not None]
        examples = sorted(eligible, key=lambda row: abs(row["mds_z"]), reverse=True)[:20]
        table_rows = []
        for row in examples:
            table_rows.append("<tr>"+"".join("<td>"+html.escape(str(value))+"</td>" for value in (
                f"{row['chrom']}:{row['start0']}-{row['end0']}", row["threshold"], row["supported_observed"],
                f"{row['expected']:.3f}", f"{row['oe']:.3f}", f"{row['mds_rawdeviance']:.3f}",
                f"{row['mds_z']:.3f}"))+"</tr>")
        sections.append(f"<h2>{html.escape(CHANNEL_LABELS[channel])}</h2><p>{len(eligible):,} window/tier rows have an empirical Z. "
                        "The following 20 largest absolute scores are descriptive review examples, not significant loci.</p>"
                        "<table><thead><tr><th>Window (BED0)</th><th>AC</th><th>Observed</th><th>Expected</th><th>O/E</th>"
                        "<th>Signed deviance</th><th>Empirical Z</th></tr></thead><tbody>"+"".join(table_rows)+"</tbody></table>")
    citation_links = "".join(f'<li><a href="{html.escape(reference["url"])}">{html.escape(reference["title"])}</a></li>' for reference in REFERENCES)
    limitations = "".join("<li>"+html.escape(value)+"</li>" for value in LIMITATIONS)
    figure_html = "".join(f'<figure><img src="{html.escape(item["png"])}" style="max-width:100%;height:auto" '
                          f'alt="Signed multi-ALT density scores for {item["window_bp"]}-base windows">'
                          f'<figcaption>{item["window_bp"]:,}-bp windows; {item["omitted_unavailable_z_rows"]:,} rows with unavailable Z are not plotted. '
                          f'<a href="{html.escape(item["svg"])}">Download SVG</a>. The zero line is a reference center, not a significance threshold.</figcaption></figure>'
                          for item in figures)
    path.write_text("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Multi-ALT density scores</title>"
                    "<style>body{font:16px system-ui,sans-serif;max-width:1180px;margin:36px auto;padding:0 24px;line-height:1.6;color:#172b3a}"
                    "h1,h2{color:#123f55}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:7px 10px;border-bottom:1px solid #d7e1e6;text-align:right}"
                    "th:first-child,td:first-child{text-align:left}li{margin-bottom:8px}.note{padding:18px;background:#edf4f7;border-left:4px solid #25728b}</style>"
                    "<h1>Multi-ALT density scores</h1><p class=\"note\">Positive scores mean an excess of count-supported multi-ALT loci relative to the declared background. "
                    "They do not demonstrate selection, tolerance or lower disease risk. REF presence is unknown in these counters.</p>"
                    f"<p>Compact input strata: {audit['source_strata_rows']:,}. Window/tier/channel rows: {len(rows):,}. "
                    "Complete results: <a href=\"../results/density_windows.tsv.gz\">density_windows.tsv.gz</a>. "
                    "Reference fitting: <a href=\"../qc/calibration.tsv\">calibration.tsv</a>. "
                    "Provenance: <a href=\"../qc/manifest.json\">manifest.json</a>.</p>"
                    "<h2>Method</h2><p>Exact context, mutation class, primary-allele frequency bin, allele-number bin and AC tier define background strata. "
                    "A joint Dirichlet(½,½,½) model describes one, exactly two and three supported ALT bases. The scored chromosome is held out. "
                    "Each standardization reference also holds out its own chromosome, so the target chromosome contributes no data to its reference fit. "
                    "Reference windows use the same scale and AC tier, do not overlap within a scale, and are included without an observed-target-count or enrichment gate. "
                    "Prespecified exclusion regions remove overlapping whole native bins from background training and veto overlapping scored/reference windows; raw counts are retained. "
                    "Empirical Z centers and scales signed aggregate binomial deviance; no normal-tail P value is calculated.</p>"
                    +figure_html+"".join(sections)+"<h2>Interpretation and limits</h2><ul>"+limitations+"</ul><h2>SDS inspiration</h2>"
                    "<p>The original SDS infers relative terminal branch lengths from per-individual nearest-singleton distances around biallelic SNPs. "
                    "This module uses a different input and statistic, borrowing only the idea of density plus conditional standardization.</p><ul>"
                    +citation_links+"</ul></html>\n", encoding="utf-8")


@transactional("qc/manifest.json", input_dirs=("store",), input_files=("scope",), input_prefixes=(), resolve_inputs=gw.protected_sources)
def run(args):
    validate_args(args)
    files, native, scopes, source_metadata, sources = load_source(args)
    masks = gw.intervals(args.exclude_bed, allow_overlap=True) if args.exclude_bed else {}
    sources["exclude_bed"] = [identity(path, full=True) for path in args.exclude_bed]
    rows, moments, audit = collect_windows(files, native, scopes, args, masks)
    if sources["counters"] != [identity(path) for path in files]:
        raise ValueError("Compact counters changed between the two passes")
    rows.sort(key=lambda row: (int(gw.chromosome(row["chrom"])), row["start0"], row["window_bp"], row["threshold"], row["channel"]))
    calibrations = standardize(rows, moments, args)
    out = Path(args.out)
    for name in ("results", "reports", "qc"):
        (out / name).mkdir(parents=True, exist_ok=True)
    gw.write_table(out/"results/density_windows.tsv.gz", rows, FIELDS)
    gw.write_table(out/"qc/calibration.tsv", calibrations, list(calibrations[0]) if calibrations else
                   ["chrom", "window_bp", "threshold", "channel", "status"])
    audit["raw_score_rows"] = sum(row["mds_rawdeviance"] is not None for row in rows)
    audit["empirical_z_rows"] = sum(row["mds_z"] is not None for row in rows)
    audit["unavailable_true_three_four_allele_rows"] = len(rows)
    for row in rows:
        audit["eligibility:"+row["score_eligibility"]] += 1
    gw.write_table(out/"qc/audit.tsv", [dict(metric=key, value=value) for key, value in sorted(audit.items())], ["metric", "value"])
    figures = plot_scores(out / "reports", rows)
    write_report(out/"reports/density.html", rows, audit, figures)
    status = "complete" if audit["empirical_z_rows"] else "complete_with_exclusions"
    return dict(
        status=status, schema=VERSION, mode="observed_multi_ALT_density_from_compact_counters", build="GRCh38", source_bin_bp=native,
        counts=dict(audit), source_metadata=source_metadata, source_identity=sources,
        parameters={name: getattr(args, name) for name in ("window_sizes", "step_bp", "thresholds", "primary_threshold", "min_sites",
                    "min_expected", "min_support", "min_training_sites", "min_calibration_windows", "max_windows")},
        channels=CHANNEL_LABELS, denominator=gw.DENOMINATOR,
        prior="Dirichlet(0.5,0.5,0.5) over [N-M,M-Q,Q] where M=at least two supported ALT bases and Q=three supported ALT bases",
        expectation="For held-out chromosome C: p_multi=(M_train+1)/(N_train+1.5), p_three=(Q_train+0.5)/(N_train+1.5); E=sum N_window,stratum*p_stratum",
        predictive_variance="sum_s N_s*p_s*(1-p_s)*(1+(N_s-1)/(N_training,s+2.5)); working marginal model, not LD calibration",
        raw_score="sign(O-E)*sqrt(2*[O*log(O/E)+(N-O)*log((N-O)/(N-E))]); zero log terms defined as zero",
        standardized_score="(mds_rawdeviance - double_LOCO_reference_mean)/double_LOCO_reference_sample_SD",
        calibration="All eligible unselected same-scale/same-tier nonoverlapping other-chromosome windows; each reference background excludes target AND reference chromosomes; >=2 remaining training chromosomes per stratum and >=2 reference chromosomes",
        inference="descriptive_scores_only_no_P_values_no_q_values_no_2x_nomination_gate",
        true_observed_allele_counts="unavailable_from_legacy_strata_REF_presence_and_joint_sample_set_not_established",
        mask_policy="prespecified_exclusion_BED_omits_overlapping_whole_native_bins_from_training_and_vetoes_overlapping_target_reference_windows; raw_window_counts_unchanged_no_prorating",
        external_true_allele_contract="A future joint-base adapter must establish common sample set, site-wide AN, passing A/C/G/T AC, REF AC and threshold before exact observed-three/four-base statistics can replace these proxies; this module does not infer them",
        source_filter_change="none; all declared source scope is scored without hotspot filtering; no raw VCF, Atlas catalogue or individual genotypes are reread",
        raw_density="raw_loci_per_mb=raw_observed/(end0-start0)*1e6; span_denominator_not_callable_corrected",
        figures=figures, references=REFERENCES, limitations=LIMITATIONS,
    )


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--store", "--scan-dir", help="Validated completed scan root containing compact strata caches")
    source.add_argument("--tiles", nargs="+", help="Explicit completed strata.tsv(.gz) counter files")
    p.add_argument("--scope", help="Original complete BED0/scan_scope.tsv for --tiles")
    p.add_argument("--build", default="GRCh38", choices=("GRCh38", "hg38", "38"))
    p.add_argument("--tile-bp", type=int, default=1000000)
    p.add_argument("--window-sizes", default="1000000")
    p.add_argument("--step-bp", type=int, default=1000000)
    p.add_argument("--thresholds", default="5,10")
    p.add_argument("--primary-threshold", type=int, default=5, choices=(5, 10))
    p.add_argument("--min-sites", type=int, default=500)
    p.add_argument("--min-expected", type=float, default=1.)
    p.add_argument("--min-support", type=float, default=1.)
    p.add_argument("--min-training-sites", type=int, default=100)
    p.add_argument("--min-calibration-windows", type=int, default=50)
    p.add_argument("--exclude-bed", action="append", default=[], help="Omit touching native bins from background and veto touching score/reference windows; raw counts retained, never prorated")
    p.add_argument("--max-windows", type=int, default=250000)
    p.add_argument("--out", required=True)
    return p


def main(argv=None):
    try:
        result = run(parser().parse_args(argv))
    except (ValueError, OSError, RuntimeError) as error:
        raise SystemExit("allele-density: "+str(error))
    print(json.dumps(dict(status=result["status"], mode=result["mode"], counts=result["counts"]), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
