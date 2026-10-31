#!/usr/bin/env python3
"""Compare HERV sequence and LD partitions on the same physical windows.

This module measures an association between an annotation and existing LD maps.
It does not estimate recombination from genotypes or identify a causal driver.
Use maps made with comparable samples/SNPs/algorithms for biological contrasts.

Example:
  python f/herv_ld.py --ld-block-dir /mnt/f/LD/block \
    --erv-joined /mnt/f/annot/allele3/repeats/hg38/rmskJoinedCurrent.txt.gz \
    --out /path/to/results/ld --chrom-sizes /path/to/hg38.fa.fai

Optional covariates: a tab-delimited table with chrom,start,end (BED0), exactly
matching the physical grid, and gc,recombination,mappability,snp_density.
Missing rows/values remain missing; they are never silently replaced with zero.
"""
from __future__ import annotations
from analysis_integrity import transactional

import argparse
from bisect import bisect_left, bisect_right
from collections import defaultdict
import csv
import gzip
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scan import parse_joined_rmsk, read_blocks  # Existing BED0/joined parser.
from allele_common import canonical_chrom

VERSION = "herv_ld_1.0"
DEFAULT_ANCESTRIES = ("AFR", "EAS", "EUR", "SAS")
DEFAULT_FAMILIES = ("ERV1", "ERVK", "ERVL", "ERVL-MaLR")
REQUIRED_COVARIATES = ("gc", "recombination", "mappability", "snp_density")
METRICS = ("boundary_density_per_mb", "log2_bp_weighted_median_block_bp")


def merge_intervals(intervals):
    """Union BED0 half-open intervals, including touching intervals."""
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError(f"Invalid interval {start}-{end}")
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def intersect_intervals(left, right):
    """Intersection of two sorted, nonoverlapping interval lists."""
    result, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        if start < end:
            result.append((start, end))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return result


class IntervalIndex:
    def __init__(self, intervals):
        self.intervals = merge_intervals(intervals)
        self.starts = [v[0] for v in self.intervals]
        self.ends = [v[1] for v in self.intervals]

    def clip(self, start, end):
        lo = bisect_right(self.ends, start)
        hi = bisect_left(self.starts, end)
        return [(max(a, start), min(b, end)) for a, b in self.intervals[lo:hi]]

    def bp(self, start, end):
        return sum(b-a for a, b in self.clip(start, end))

    def spans_boundary(self, position):
        # Both bases flanking an LD junction must belong to common coverage.
        i = bisect_right(self.starts, position)-1
        return i >= 0 and self.starts[i] < position < self.ends[i]


def weighted_median(values, weights):
    values, weights = np.asarray(values), np.asarray(weights)
    if not len(values) or not np.any(weights > 0):
        return float("nan")
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[order])
    return float(values[order[np.searchsorted(cumulative, cumulative[-1]/2)]])


def read_sizes(path):
    result = {}
    with open(path) as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            row = line.split()
            chrom = canonical_chrom(row[0])
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            length = int(row[1])
            if length <= 0 or chrom in result:
                raise ValueError(f"Invalid/duplicate chromosome size: {line.strip()}")
            result[chrom] = length
    if not result:
        raise ValueError("Chromosome sizes contain no autosomes")
    return result


def build_profiles(block_sets, elements, *, sizes=None, window_bp=1_000_000,
                   min_common_coverage=.95, min_window_fraction=.5,
                   families=DEFAULT_FAMILIES):
    """One row per physical window and ancestry; use identical common bases.

    A boundary is a shared end/start of two adjacent blocks, not an edge next to
    an unobserved gap. Block lengths refer to original blocks; their weights are
    the number of common-covered bases in this window, not clipped lengths.
    """
    if window_bp < 1 or not 0 < min_common_coverage <= 1 or not 0 < min_window_fraction <= 1:
        raise ValueError("Invalid physical-window settings")
    autosomes = {str(i) for i in range(1, 23)}
    ancestry = list(block_sets)
    by_anc, max_ends = {}, defaultdict(int)
    for anc, blocks in block_sets.items():
        by_chrom = defaultdict(list)
        for block in blocks:
            if block.chrom not in autosomes:
                continue
            by_chrom[block.chrom].append(block)
            max_ends[block.chrom] = max(max_ends[block.chrom], block.end)
        by_anc[anc] = by_chrom
    lengths = dict(sizes or max_ends)
    for chrom, end in max_ends.items():
        if chrom not in lengths or end > lengths[chrom]:
            raise ValueError(f"LD intervals exceed chromosome sizes on chr{chrom}")
    if not lengths:
        raise ValueError("No autosomal LD intervals")

    repeats = defaultdict(list)
    for element in elements:
        if element.rep_class == "LTR" and element.family in families and element.chrom in lengths:
            for start, end in element.segments:
                if end > lengths[element.chrom] and sizes is not None:
                    raise ValueError(f"Repeat annotation exceeds chromosome size on chr{element.chrom}")
                repeats[element.chrom, element.family].append((start, end))
    repeat_index = {key: IntervalIndex(value) for key, value in repeats.items()}
    empty = IntervalIndex([])
    rows, windows = [], []
    for chrom in sorted(lengths, key=int):
        coverage, block_index, junctions = {}, {}, {}
        for anc in ancestry:
            blocks = sorted(by_anc[anc].get(chrom, []), key=lambda b: b.start)
            coverage[anc] = IntervalIndex((b.start, b.end) for b in blocks)
            block_index[anc] = (blocks, [b.start for b in blocks], [b.end for b in blocks])
            junctions[anc] = [right.start for left, right in zip(blocks, blocks[1:])
                              if left.end == right.start]
        common = coverage[ancestry[0]].intervals
        for anc in ancestry[1:]:
            common = intersect_intervals(common, coverage[anc].intervals)
        common_index = IntervalIndex(common)
        for start in range(0, lengths[chrom], window_bp):
            end = min(start+window_bp, lengths[chrom])
            shared = common_index.clip(start, end)
            common_bp = sum(b-a for a, b in shared)
            eligible = (common_bp/(end-start) >= min_common_coverage
                        and (end-start)/window_bp >= min_window_fraction)
            base = dict(chrom=chrom, start=start, end=end, window_bp=end-start,
                        common_bp=common_bp, common_coverage_fraction=common_bp/(end-start),
                        eligible=int(eligible))
            for family in families:
                index = repeat_index.get((chrom, family), empty)
                bp = sum(index.bp(a, b) for a, b in shared)
                base[family+"_aligned_bp"] = bp
                base[family+"_fraction"] = bp/common_bp if common_bp else float("nan")
            windows.append(base)
            for anc in ancestry:
                blocks, starts, ends = block_index[anc]
                overlap_weights, original_lengths = [], []
                # A block may cross several common-coverage fragments: count it
                # once and sum its overlap before taking either median.
                i, j = bisect_right(ends, start), bisect_left(starts, end)
                for block in blocks[i:j]:
                    overlap = sum(max(0, min(b, block.end)-max(a, block.start)) for a, b in shared)
                    if overlap:
                        original_lengths.append(block.end-block.start)
                        overlap_weights.append(overlap)
                bounds = junctions[anc]
                lo, hi = bisect_left(bounds, start), bisect_left(bounds, end)
                n_boundary = sum(common_index.spans_boundary(p) for p in bounds[lo:hi])
                median = weighted_median(original_lengths, overlap_weights)
                rows.append(base | dict(
                    ancestry=anc, covered_bp=coverage[anc].bp(start, end),
                    coverage_fraction=coverage[anc].bp(start, end)/(end-start),
                    n_blocks_intersecting_common=len(original_lengths),
                    median_intersecting_block_bp=float(np.median(original_lengths)) if original_lengths else float("nan"),
                    bp_weighted_median_block_bp=median,
                    log2_bp_weighted_median_block_bp=math.log2(median) if median > 0 else float("nan"),
                    internal_boundary_count=n_boundary,
                    boundary_density_per_mb=n_boundary*1_000_000/common_bp if common_bp else float("nan")))
    return rows, windows, lengths


def load_covariates(path, windows, columns=REQUIRED_COVARIATES):
    """Join covariates to exact, canonical BED0 physical-window keys."""
    if not path:
        return {}, [], {"provided": False, "matched_windows": 0}
    opener = gzip.open if str(path).endswith(".gz") else open
    target = {(w["chrom"], w["start"], w["end"]) for w in windows}
    values, ignored = {}, 0
    with opener(path, "rt") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"chrom", "start", "end"}.issubset(reader.fieldnames or []):
            raise ValueError("Covariate TSV requires chrom,start,end (BED0 physical grid)")
        present = [name for name in columns if name in reader.fieldnames]
        if not present:
            raise ValueError("Covariate TSV has none of the requested covariate columns")
        seen = set()
        for row in reader:
            key = (canonical_chrom(row["chrom"]), int(row["start"]), int(row["end"]))
            if key in seen:
                raise ValueError(f"Duplicate covariate window: {key}")
            seen.add(key)
            if key not in target:
                ignored += 1
                continue
            result = {}
            for name in present:
                raw = row[name].strip()
                result[name] = float("nan") if raw in {"", ".", "NA", "NaN", "nan"} else float(raw)
                if math.isinf(result[name]):
                    raise ValueError(f"Infinite covariate {name} at {key}")
            values[key] = result
    if not values:
        raise ValueError("No covariate rows match the physical window grid; check build/window size")
    return values, present, dict(provided=True, matched_windows=len(values),
                                missing_windows=len(target)-len(values), ignored_rows=ignored)


def _center_by_chromosome(matrix, weights, chromosomes):
    centered = np.asarray(matrix, float).copy()
    for chrom in np.unique(chromosomes):
        use = chromosomes == chrom
        centered[use] -= np.average(centered[use], weights=weights[use], axis=0)
    return centered


def fit_association(x, y, chromosomes, weights, covariates=None, covariate_names=(), *,
                    bootstrap=1000, seed=20261027, min_chromosomes=8, min_windows=30):
    """Weighted within-chromosome slope, resampling whole chromosomes.

    x is HERV percentage (0--100), so beta is per percentage point. Chromosome
    fixed effects are removed by weighted centering. Paired population-difference
    outcomes use this same routine and the same physical windows.
    """
    x, y, weights = map(lambda v: np.asarray(v, float), (x, y, weights))
    chromosomes = np.asarray(chromosomes, str)
    covariates = np.empty((len(x), 0)) if covariates is None else np.asarray(covariates, float)
    if covariates.ndim != 2 or len(covariates) != len(x):
        raise ValueError("Covariates must have one row per window")
    use = (np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0)
           & np.isfinite(covariates).all(axis=1))
    x, y, chromosomes, weights, covariates = [v[use] for v in (x, y, chromosomes, weights, covariates)]
    groups = np.unique(chromosomes)
    result = dict(n_windows=len(x), n_chromosomes=len(groups), beta_per_percentage_point=None,
                  ci_low=None, ci_high=None, bootstrap_p=None, fdr=None,
                  bootstrap_successful=0, test_status="insufficient_windows_or_chromosomes",
                  covariates_used="", covariates_absorbed_or_redundant="")
    if len(x) < min_windows or len(groups) < min_chromosomes:
        return result
    within = _center_by_chromosome(np.column_stack([x, y, covariates]), weights, chromosomes)
    x, yc, cc = within[:, 0], within[:, 1], within[:, 2:]
    if np.sum(weights*yc*yc) <= 1e-20:
        result.update(beta_per_percentage_point=None, test_status="no_within_chromosome_outcome_variation")
        return result
    if np.sum(weights*x*x) <= 1e-15:
        result["test_status"] = "no_within_chromosome_exposure_variation"
        return result
    kept, dropped = [], []
    cdesign = np.empty((len(x), 0))
    # Choose identifiable covariates without considering the exposure. A
    # confounder collinear with HERV must not be dropped to rescue its estimate.
    for i, name in enumerate(covariate_names):
        scale = math.sqrt(np.average(cc[:, i]**2, weights=weights))
        if scale <= 1e-10:
            dropped.append(name)
            continue
        column = cc[:, i]/scale
        trial = np.column_stack([cdesign, column])
        if np.linalg.matrix_rank(trial*np.sqrt(weights)[:, None]) > cdesign.shape[1]:
            kept.append(name)
            cdesign = trial
        else:
            dropped.append(name)
    result["covariates_used"] = ",".join(kept)
    result["covariates_absorbed_or_redundant"] = ",".join(dropped)
    design = np.column_stack([x, cdesign])
    weighted = design*np.sqrt(weights)[:, None]
    if np.linalg.matrix_rank(weighted) != design.shape[1]:
        result["test_status"] = "exposure_not_identifiable_after_adjustment"
        return result
    xx, xy = [], []
    for chrom in groups:
        use = chromosomes == chrom
        xx.append(design[use].T @ (design[use]*weights[use, None]))
        xy.append(design[use].T @ (yc[use]*weights[use]))
    xx, xy = np.asarray(xx), np.asarray(xy)
    beta = float(np.linalg.solve(xx.sum(axis=0), xy.sum(axis=0))[0])
    result.update(beta_per_percentage_point=beta, test_status="estimated_no_bootstrap")
    if bootstrap == 0:
        return result
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(bootstrap):
        counts = rng.multinomial(len(groups), np.full(len(groups), 1/len(groups)))
        gram = np.tensordot(counts, xx, axes=1)
        if np.linalg.matrix_rank(gram) != len(gram):
            continue
        estimate = np.linalg.solve(gram, np.tensordot(counts, xy, axes=1))[0]
        if np.isfinite(estimate):
            estimates.append(estimate)
    result["bootstrap_successful"] = len(estimates)
    if len(estimates) < max(20, .8*bootstrap):
        result["test_status"] = "bootstrap_unstable"
        return result
    estimates = np.asarray(estimates)
    lower, upper = np.quantile(estimates, [.025, .975])
    # Centered two-sided bootstrap test; descriptive association, not a causal
    # or randomized-test P value. BH correction spans all reported slope tests.
    pvalue = (1+np.sum(np.abs(estimates-beta) >= abs(beta)))/(len(estimates)+1)
    result.update(ci_low=float(lower), ci_high=float(upper), bootstrap_p=float(pvalue), test_status="complete")
    return result


def _mean_difference(y, chromosomes, weights, bootstrap, seed):
    groups = np.unique(chromosomes)
    sums = np.array([np.sum(y[chromosomes == c]*weights[chromosomes == c]) for c in groups])
    denominators = np.array([np.sum(weights[chromosomes == c]) for c in groups])
    result = dict(mean_difference=float(sums.sum()/denominators.sum()), mean_ci_low=None, mean_ci_high=None)
    if bootstrap and len(groups) >= 3:
        counts = np.random.default_rng(seed).multinomial(len(groups), np.full(len(groups), 1/len(groups)), bootstrap)
        values = counts @ sums/(counts @ denominators)
        result.update(mean_ci_low=float(np.quantile(values, .025)), mean_ci_high=float(np.quantile(values, .975)))
    return result


def analyze(profiles, families, ancestries, covariate_values, covariate_names, *,
            bootstrap=1000, seed=20261027, min_chromosomes=8, min_windows=30):
    index = defaultdict(dict)
    for row in profiles:
        if row["eligible"]:
            index[row["chrom"], row["start"], row["end"]][row["ancestry"]] = row
    keys = sorted(index, key=lambda k: (int(k[0]), k[1]))
    # Every slope and every paired contrast shares a complete-case window set.
    keys = [key for key in keys if set(index[key]) == set(ancestries)
            and all(np.isfinite(index[key][anc][metric]) for anc in ancestries for metric in METRICS)
            and all(np.isfinite(covariate_values.get(key, {}).get(name, np.nan)) for name in covariate_names)]
    if set(REQUIRED_COVARIATES).issubset(covariate_names):
        level = "adjusted_observational_association"
    elif covariate_names:
        level = "descriptive_partial_adjustment"
    else:
        level = "descriptive_without_sequence_recombination_covariates"
    base = [index[key][ancestries[0]] for key in keys]
    chromosomes = np.asarray([key[0] for key in keys])
    weights = np.asarray([row["common_bp"]/1_000_000 for row in base])
    covariates = np.array([[covariate_values[key][name] for name in covariate_names] for key in keys]) if covariate_names else np.empty((len(keys), 0))
    if covariate_names and not keys:
        covariates = np.empty((0, len(covariate_names)))
    association, contrasts = [], []
    for family in families:
        x = np.asarray([row[family+"_fraction"]*100 for row in base])
        for metric in METRICS:
            outcomes = {anc: np.asarray([index[key][anc][metric] for key in keys]) for anc in ancestries}
            kwargs = dict(bootstrap=bootstrap, seed=seed, min_chromosomes=min_chromosomes, min_windows=min_windows)
            for anc in ancestries:
                fit = fit_association(x, outcomes[anc], chromosomes, weights, covariates, covariate_names, **kwargs)
                association.append(dict(family=family, ancestry=anc, metric=metric, evidence_level=level) | fit)
            for first, second in itertools.combinations(ancestries, 2):
                delta = outcomes[first]-outcomes[second]
                fit = fit_association(x, delta, chromosomes, weights, covariates, covariate_names, **kwargs)
                means = _mean_difference(delta, chromosomes, weights, bootstrap, seed) if keys else dict(mean_difference=None, mean_ci_low=None, mean_ci_high=None)
                contrasts.append(dict(family=family, ancestry_a=first, ancestry_b=second, metric=metric,
                                      contrast=f"{first}-minus-{second}", evidence_level=level) | means | fit)
    tests = [row for row in association+contrasts if row["bootstrap_p"] is not None]
    ordered = sorted(tests, key=lambda row: row["bootstrap_p"])
    previous = 1.
    for i in range(len(ordered)-1, -1, -1):
        previous = min(previous, ordered[i]["bootstrap_p"]*len(ordered)/(i+1))
        ordered[i]["fdr"] = previous
    return association, contrasts, dict(eligible_physical_windows=len(index), complete_case_windows=len(keys),
                                       dropped_incomplete_windows=len(index)-len(keys), evidence_level=level,
                                       slope_tests=len(association)+len(contrasts), tested_slopes=len(tests),
                                       fdr_scope="all families, both outcomes, all ancestry-specific and paired-difference slope tests")


def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No output rows for {path}")
    path = Path(path)
    temporary = path.with_name(path.name+".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(temporary, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if value is None or isinstance(value, (float, np.floating)) and not np.isfinite(value)
                             else value for key, value in row.items()})
    temporary.replace(path)


def plot_results(profiles, association, contrasts, families, ancestries, sizes, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 180, "pdf.fonttype": 42})
    eligible = [row for row in profiles if row["eligible"]]
    if not eligible:
        return []
    offsets, total = {}, 0
    for chrom in sorted(sizes, key=int):
        offsets[chrom] = total
        total += sizes[chrom]
    palette = ("#9C4A75", "#177D8D", "#D18E24", "#547749", "#6966A6", "#626262")
    colors = {anc: palette[i % len(palette)] for i, anc in enumerate(ancestries)}
    files = []
    def save(fig, name):
        for extension in ("png", "pdf"):
            file = out/(name+"."+extension)
            fig.savefig(file, bbox_inches="tight")
            files.append(file.name)
        plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True)
    first = [r for r in eligible if r["ancestry"] == ancestries[0]]
    x = [(offsets[r["chrom"]]+(r["start"]+r["end"])/2)/1e6 for r in first]
    for family in families:
        axes[0].plot(x, [r[family+"_fraction"]*100 for r in first], lw=.7, alpha=.8, label=family)
    axes[0].set_ylabel("Aligned HERV sequence (%)\nof common-covered bases")
    axes[0].legend(ncol=len(families), frameon=False, loc="upper right")
    for anc in ancestries:
        group = [r for r in eligible if r["ancestry"] == anc]
        gx = [(offsets[r["chrom"]]+(r["start"]+r["end"])/2)/1e6 for r in group]
        axes[1].scatter(gx, [r["boundary_density_per_mb"] for r in group], s=3, alpha=.6, label=anc, color=colors[anc])
        axes[2].scatter(gx, [r["bp_weighted_median_block_bp"]/1e3 for r in group], s=3, alpha=.6, label=anc, color=colors[anc])
    axes[1].set_ylabel("Internal LD boundaries\nper common-covered Mb")
    axes[1].legend(ncol=len(ancestries), frameon=False)
    axes[2].set_ylabel("Base-weighted median\nLD block length (kb)")
    axes[2].set_yscale("log")
    ticks = [(offsets[c]+sizes[c]/2)/1e6 for c in sorted(sizes, key=int)]
    axes[2].set_xticks(ticks, ["chr"+c for c in sorted(sizes, key=int)], rotation=45)
    axes[2].set_xlabel("Identical physical windows and common LD-map coverage")
    for ax in axes:
        for offset in offsets.values():
            ax.axvline(offset/1e6, color="#CCCCCC", lw=.4, zorder=0)
    fig.suptitle("HERV sequence and LD-map architecture | descriptive association")
    fig.tight_layout()
    save(fig, "ld_landscape")

    fig, axes = plt.subplots(1, 2, figsize=(12, max(4, len(families)*1.2)), sharey=True)
    for ax, metric in zip(axes, METRICS):
        for ai, anc in enumerate(ancestries):
            group = {r["family"]: r for r in association if r["ancestry"] == anc and r["metric"] == metric}
            for fi, family in enumerate(families):
                row = group[family]
                beta = row["beta_per_percentage_point"]
                if beta is None:
                    continue
                yy = fi+(ai-(len(ancestries)-1)/2)*.13
                ax.plot(beta, yy, "o", ms=4, color=colors[anc], label=anc if fi == 0 else None)
                if row["ci_low"] is not None:
                    ax.plot([row["ci_low"], row["ci_high"]], [yy, yy], color=colors[anc], lw=1.2)
        ax.axvline(0, color="#555555", lw=.8)
        ax.set_xlabel("LD outcome change per +1 percentage point HERV\n(chromosome fixed effects; chromosome bootstrap 95% CI)")
        ax.set_title("Boundary density per Mb" if metric == METRICS[0] else "log2 block length")
        ax.set_yticks(range(len(families)), families)
        ax.grid(axis="x", color="#EEEEEE")
    axes[0].legend(frameon=False, ncol=2)
    fig.suptitle("Within-chromosome association; no causal interpretation")
    fig.tight_layout()
    save(fig, "ld_association")

    pairs = [f"{a}-minus-{b}" for a, b in itertools.combinations(ancestries, 2)]
    fig, axes = plt.subplots(1, 2, figsize=(12, max(4, len(pairs)*.5)))
    for ax, metric in zip(axes, METRICS):
        group = {(r["contrast"], r["family"]): r for r in contrasts if r["metric"] == metric}
        matrix = np.array([[group[pair, family]["beta_per_percentage_point"]
                            if group[pair, family]["beta_per_percentage_point"] is not None else np.nan
                            for family in families] for pair in pairs], dtype=float)
        finite = np.isfinite(matrix)
        maximum = max(float(np.max(np.abs(matrix[finite]))) if finite.any() else 0, 1e-8)
        im = ax.imshow(np.ma.masked_invalid(matrix), vmin=-maximum, vmax=maximum, cmap="RdBu_r", aspect="auto")
        ax.set_xticks(range(len(families)), families, rotation=30)
        ax.set_yticks(range(len(pairs)), [p.replace("-minus-", " − ") for p in pairs])
        ax.set_title("Boundary density contrast" if metric == METRICS[0] else "log2 block-length contrast")
        for i, pair in enumerate(pairs):
            for j, family in enumerate(families):
                row = group[pair, family]
                if row["fdr"] is not None and row["fdr"] < .05:
                    ax.text(j, i, "•", ha="center", va="center", color="black")
        fig.colorbar(im, ax=ax, fraction=.04, label="Paired slope difference per +1 pp HERV")
    fig.suptitle("Population contrasts at the same physical windows | dot: BH FDR < 0.05")
    fig.tight_layout()
    save(fig, "ld_ancestry_contrasts")
    return files


def file_info(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ld-block-dir", default="/mnt/f/LD/block")
    p.add_argument("--erv-joined", default="/mnt/f/annot/allele3/repeats/hg38/rmskJoinedCurrent.txt.gz")
    p.add_argument("--out", required=True)
    p.add_argument("--build", choices=("37", "38"), default="38")
    p.add_argument("--annotation-build", choices=("37", "38"), default="38", help="Declared build of joined annotation; must match --build")
    p.add_argument("--ancestries", default=",".join(DEFAULT_ANCESTRIES))
    p.add_argument("--families", default=",".join(DEFAULT_FAMILIES))
    p.add_argument("--chrom-sizes", help="Two-column chromosome sizes or FASTA .fai; inferred from LD envelope if omitted")
    p.add_argument("--window-bp", type=int, default=1_000_000)
    p.add_argument("--min-common-coverage", type=float, default=.95)
    p.add_argument("--min-window-fraction", type=float, default=.5)
    p.add_argument("--covariates", help="BED0 exact-grid TSV: chrom,start,end plus numerical covariates")
    p.add_argument("--covariate-columns", default=",".join(REQUIRED_COVARIATES))
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261027)
    p.add_argument("--min-chromosomes", type=int, default=8)
    p.add_argument("--min-windows", type=int, default=30)
    p.add_argument("--no-plots", action="store_true")
    return p


@transactional("manifest.json")
def run(args):
    if args.build != args.annotation_build:
        raise ValueError("LD BED and RepeatMasker builds must match; this module does not liftover")
    if args.bootstrap < 0 or args.min_chromosomes < 3 or args.min_windows < 3:
        raise ValueError("Require bootstrap >= 0, min-chromosomes >= 3, min-windows >= 3")
    ancestries = tuple(v.strip() for v in args.ancestries.split(",") if v.strip())
    families = tuple(v.strip() for v in args.families.split(",") if v.strip())
    if len(ancestries) < 2 or len(set(ancestries)) != len(ancestries) or not families or len(set(families)) != len(families):
        raise ValueError("Require >=2 unique ancestries and unique nonempty families")
    paths = {anc: Path(args.ld_block_dir)/f"{anc}.{args.build}.bed" for anc in ancestries}
    sets, block_stats = {}, {}
    for anc, path in paths.items():
        sets[anc], block_stats[anc] = read_blocks(path, anc, args.build)
        if not sets[anc]:
            raise ValueError(f"Empty LD partition: {path}")
    elements = parse_joined_rmsk(args.erv_joined)
    sizes = read_sizes(args.chrom_sizes) if args.chrom_sizes else None
    profiles, windows, lengths = build_profiles(sets, elements, sizes=sizes, window_bp=args.window_bp,
        min_common_coverage=args.min_common_coverage, min_window_fraction=args.min_window_fraction, families=families)
    columns = tuple(v.strip() for v in args.covariate_columns.split(",") if v.strip())
    covariate_values, covariate_names, covariate_qc = load_covariates(args.covariates, windows, columns)
    for row in profiles:
        values = covariate_values.get((row["chrom"], row["start"], row["end"]), {})
        row.update({"cov_"+name: values.get(name, float("nan")) for name in covariate_names})
    association, contrasts, analysis_qc = analyze(profiles, families, ancestries, covariate_values, covariate_names,
        bootstrap=args.bootstrap, seed=args.seed, min_chromosomes=args.min_chromosomes, min_windows=args.min_windows)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out/"window_profiles.tsv.gz", profiles)
    write_tsv(out/"association.tsv", association)
    write_tsv(out/"ancestry_contrasts.tsv", contrasts)
    figures = [] if args.no_plots else plot_results(profiles, association, contrasts, families, ancestries, lengths, out)
    manifest = dict(version=VERSION, status="complete", build=args.build, annotation_build=args.annotation_build,
        inputs={"ld_partitions": {anc: file_info(path) for anc, path in paths.items()},
                "erv_joined": file_info(args.erv_joined),
                "chrom_sizes": file_info(args.chrom_sizes) if args.chrom_sizes else None,
                "covariates": file_info(args.covariates) if args.covariates else None},
        settings=vars(args), block_counts=block_stats,
        coordinate_convention="BED0 half-open; autosomes only",
        genome_extent="chromosome sizes" if sizes else "largest LD interval end per chromosome; not a callable genome mask",
        common_coverage="intersection of interval unions of ALL requested ancestries, identically used for every exposure and outcome",
        boundary_definition="junctions of touching blocks with both flanking bases in common coverage; gap edges are excluded",
        repeat_definition="union of aligned joined RepeatMasker segments per family; consensus tails and nested insertion gaps excluded",
        block_length_definition="length of original block, weighted by common bases overlapped in the physical window",
        model="WLS within chromosome (chromosome fixed effects); per +1 percentage point HERV aligned fraction; complete cases shared across all tests",
        uncertainty="whole-chromosome bootstrap, percentile 95% CI; centered two-sided bootstrap P; global BH for slopes",
        covariates=covariate_qc | {"included": covariate_names, "missing_key_covariates": sorted(set(REQUIRED_COVARIATES)-set(covariate_names))},
        analysis=analysis_qc,
        interpretation="Association between repeat annotation and supplied LD partitions; no genotypic LD recalculation, recombination causality, or viral driving-force inference",
        map_comparability="Sample size, SNP ascertainment, LD method and map construction must be comparable; this module cannot establish that from BED intervals",
        tables=["window_profiles.tsv.gz", "association.tsv", "ancestry_contrasts.tsv"], figures=figures)
    temporary = out/"manifest.json.tmp"
    temporary.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
    temporary.replace(out/"manifest.json")
    return manifest


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        result = run(args)
    except (ValueError, OSError) as error:
        raise SystemExit(f"herv_ld: {error}") from error
    print(json.dumps({"status": result["status"], "out": str(Path(args.out).resolve()), **result["analysis"]}))


if __name__ == "__main__":
    main()
