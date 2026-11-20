#!/usr/bin/env python3
"""Annotation-blind genomic windows for count-supported multi-ALT discovery.

The completed scanner's strata.tsv.gz is a sufficient-count cache at a fixed
native resolution. It contains a THRESHOLD-SPECIFIC denominator: a site enters
N only if at least one observed passing nonreference allele has AC>=threshold.
O counts sites with at least two such alleles. REF presence is not established.
This differs from the fixed base-eligible denominator of the HERV architecture
counters. No HERV annotation is used to construct or rank these windows.

--store verifies and reuses completed strata caches; --tiles + --scope accepts
the same explicit counter schema. Both read compact counters twice, never raw
VCF or per-site arrays. Only aligned unions of whole source bins are possible.
--summary screens/aggregates an existing window_summary.tsv, retaining its
historical expected values and fallback support. It cannot reconstruct finer
bins, stricter base filters, joint four-base counts or an independent cohort.

One-sided negative-binomial (or Poisson-limit) probabilities are WORKING-MODEL
screening probabilities. Dispersion is estimated using all eligible,
nonoverlapping reference windows on OTHER chromosomes, before any observed/OE
gate. BY correction covers all estimable sizes and AC tiers in this invocation.
This does not prove LD-calibrated frequentist FDR, biological tolerance,
positive/balancing selection, mutation birth or extinction. Nominations require
independent quality review and replication. O/E>=2 is a frozen effect-size gate,
not a promise that any region will pass.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import itertools
import json
import math
from pathlib import Path
import re
import time

import numpy as np
from scipy import stats

from analysis_integrity import InputPath, identity, transactional

VERSION = "genome-windows-1.1"
AMBIGUOUS_TO_N = str.maketrans("RYSWKMBDHV", "NNNNNNNNNN")
AUTOSOMES = {str(i) for i in range(1, 23)}
MISSING = {"", ".", "NA", "NAN", "NONE", "NULL"}
COUNTER_FIELDS = {"chrom", "start0", "end0", "threshold", "context", "primary_change",
                  "primary_af_bin", "an_bin", "n_sites", "n_multi", "n_three"}
WINDOW_FIELDS = [
    "window_id", "chrom", "start0", "end0", "build", "window_bp", "step_bp", "threshold", "is_primary",
    "raw_n_eligible", "raw_observed_multi", "raw_three_alt_sites", "n_eligible", "observed_multi",
    "expected_multi", "oe", "full_scope_oe", "n_unsupported", "unsupported_observed_multi",
    "background_support_fraction", "spatial_support_fraction", "source_bins_present", "source_bins_required",
    "strata_count", "training_exact_sites", "training_context_change_sites", "training_threshold_sites",
    "background_rate_variance", "expected_status", "expected_model", "source_mode", "denominator",
    "mask_overlap_bp", "mask_overlap_fraction", "mask_status", "calibration_windows", "calibration_chromosomes",
    "nb2_alpha", "working_nb2_alpha", "ci95_status", "ci95_low", "ci95_high", "p_value", "q_value", "fdr_method",
    "inference_status", "effect_gate_pass", "count_gate_pass", "strict_2x_supported", "nominated",
    "nomination_reason",
]
DENOMINATOR = "observed_passing_SNV_sites_with_at_least_one_nonREF_ALT_AC_ge_threshold"
LIMITATIONS = [
    "IUPAC ambiguous flanks are pooled as N without inferring a nucleotide. All source counts are retained; identical normalized strata within a source bin are summed. The central reference remains unambiguous C/T.",
    "N is threshold-specific, unlike the fixed base-eligible denominator in the HERV architecture analysis. O counts >=2 supported nonreference alleles; the reference label does not establish an observed third allele in a common sample set.",
    "Discovery is blind to HERV labels. The underlying source still has variant ascertainment, ancestry, sample-size, mapping, calling and QC effects; it is not a mutation-rate or genome-wide callable-base estimator.",
    "Unobserved source bins in summary mode are not zero. Coarser windows require complete spatial support. Counts are never split, interpolated or prorated into finer windows or around masks.",
    "Supported O and N are reduced together when exact LOCO support is missing. Raw O/N and unsupported counts remain visible. Full-scope O/E remains NA if any source stratum is unsupported.",
    "The NB2 mean/dispersion model is a screening approximation, not the exact joint distribution of linked sites. Dispersion uses unselected nonoverlapping windows from other chromosomes; expected values of those reference windows retain their own fitted backgrounds and may share training data.",
    "Conditional-model confidence intervals treat fitted expected values and residual dispersion as given. They do not propagate all background-fit, LD, callability, multiple-comparison or winner's-curse uncertainty.",
    "Each reported finite confidence interval must reproduce both target tail probabilities. Intervals outside finite floating-point representation or failing this check are NA and cannot support nomination; their valid model-tail probabilities remain in the declared BY family.",
    "BY controls the declared multiple-testing family only when its component P values are valid. The working-model probabilities here have not been independently calibrated for every genomic region; q values are not a guarantee of biological false-discovery rate.",
    "AC10 uses the same underlying cohort and is sensitivity evidence, not independent replication. Testing several sizes/AC tiers is included in one BY family; changing parameters after viewing results requires reporting the additional search.",
    "A whole-window exclusion BED is a nomination veto with overlap audit. It does not recreate base-level strict-mask counters, change source counts, or certify all remaining windows as technically reliable.",
    "Merged loci are unions of overlapping nominated windows; the union is not a newly tested interval. Lead-window estimates and minimum member q are descriptive labels; counts/expectations/ratios are never added across overlapping windows.",
    "Higher occupancy does not establish tolerance, pathogenicity or positive/balancing selection. Functional and ancient analyses must use independently justified inputs, matched controls, ancestry/geography/callability checks and replication.",
]


def chromosome(value):
    c = re.sub(r"^chr", "", str(value).strip(), flags=re.I).upper()
    return "MT" if c == "M" else c


def integer(value, label, minimum=0):
    if not re.fullmatch(r"[0-9]+", str(value).strip()):
        raise ValueError(f"{label}: expected an integer, got {value!r}")
    result = int(value)
    if result < minimum:
        raise ValueError(f"{label}: must be >= {minimum}")
    return result


def number(value, label, missing=True):
    if str(value).strip().upper() in MISSING:
        if missing:
            return None
        raise ValueError(f"{label}: missing value")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: expected a finite number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label}: nonfinite value")
    return result


def ints(value, label):
    values = [integer(v, label, 1) for v in value.split(",")]
    if not values or len(set(values)) != len(values):
        raise ValueError(f"{label}: supply distinct positive comma-separated integers")
    return sorted(values)


def table(path, required):
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)) or not required <= set(reader.fieldnames):
            raise ValueError(f"{path}: requires unique TSV headers including {sorted(required)}")
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"{path}:{line}: malformed TSV record")
            yield {k: v.strip() for k, v in row.items()}


def intervals(paths, *, allow_overlap):
    grouped = defaultdict(list)
    for path in paths:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line, text in enumerate(handle, 1):
                if not text.strip() or text.startswith(("#", "track ", "browser ")):
                    continue
                cells = text.rstrip("\r\n").split("\t")
                if cells[:3] in (["chrom", "start0", "end0"], ["chrom", "start", "end"]):
                    continue
                if len(cells) < 3:
                    raise ValueError(f"{path}:{line}: need BED0 chromosome/start/end")
                c = chromosome(cells[0])
                lo, hi = integer(cells[1], "BED start"), integer(cells[2], "BED end", 1)
                if hi <= lo:
                    raise ValueError(f"{path}:{line}: empty/reversed interval")
                if c in AUTOSOMES:
                    grouped[c].append((lo, hi))
    merged = {}
    for c, values in grouped.items():
        result = []
        for lo, hi in sorted(values):
            if result and lo < result[-1][1] and not allow_overlap:
                raise ValueError(f"Overlapping source scan scopes on chr{c}; duplicated shards cannot be added")
            if result and lo <= result[-1][1]:
                result[-1] = (result[-1][0], max(result[-1][1], hi))
            else:
                result.append((lo, hi))
        merged[c] = result
    return merged


class IntervalOverlapIndex:
    """Build chromosome start arrays once, without changing scope dictionaries."""
    def __init__(self, intervals_by_chrom):
        self.intervals = intervals_by_chrom
        self.starts = {c: tuple(lo for lo, _ in values) for c, values in intervals_by_chrom.items()}


def overlap_bp(index, c, start, end):
    c = chromosome(c)
    values = index.intervals.get(c, ()) if isinstance(index, IntervalOverlapIndex) else index.get(c, ())
    if not values:
        return 0
    if isinstance(index, IntervalOverlapIndex):
        i = max(0, bisect_right(index.starts[c], start) - 1)
    else:
        # Retain compatibility with plain interval mappings without an O(M)
        # start-list allocation. Repeated scan queries use the index above.
        i = max(0, bisect_right(values, start, key=lambda v: v[0]) - 1)
    total = 0
    while i < len(values) and values[i][0] < end:
        lo, hi = values[i]
        total += max(0, min(hi, end) - max(lo, start))
        i += 1
    return total


def counter_key(row):
    return (row["threshold"], row["context"], row["primary_change"], row["primary_af_bin"], row["an_bin"])


def read_counter_file(path, native, thresholds, scopes):
    last_tile, keys = None, set()
    for r in table(path, COUNTER_FIELDS):
        c = chromosome(r["chrom"])
        if c not in AUTOSOMES:
            continue
        t = integer(r["threshold"], "threshold", 1)
        if t not in thresholds:
            continue
        lo, hi = integer(r["start0"], "start0"), integer(r["end0"], "end0", 1)
        if lo % native or hi <= lo or hi > lo + native:
            raise ValueError(f"{path}: source interval {c}:{lo}-{hi} does not match --tile-bp={native}; bins cannot be subdivided")
        if not any(a <= lo and hi <= b and hi == min(lo+native, b) for a, b in scopes.get(c, ())):
            raise ValueError(f"{path}: source bin {c}:{lo}-{hi} is outside its recorded scan scope")
        tile = (int(c), lo, hi)
        if last_tile is not None and tile < last_tile:
            raise ValueError(f"{path}: source bins must be chromosome/position sorted")
        if tile != last_tile:
            keys.clear()
        last_tile = tile
        n, o, third = (integer(r[k], k) for k in ("n_sites", "n_multi", "n_three"))
        if not 0 <= third <= o <= n <= hi-lo:
            raise ValueError(f"{path}: require 0<=n_three<=n_multi<=n_sites<=bin span")
        if not re.fullmatch(r"[ACGTRYSWKMBDHVN]{3}", r["context"]) or r["context"][1] not in "CT":
            raise ValueError(f"{path}: context must be an IUPAC trinucleotide with an unambiguous C/T centre")
        if not re.fullmatch(r"[CT]>[ACGT]", r["primary_change"]) or r["primary_change"][0] != r["context"][1] or r["primary_change"][0] == r["primary_change"][-1]:
            raise ValueError(f"{path}: primary_change disagrees with canonical context")
        if any(r[k].upper() in MISSING for k in ("primary_af_bin", "an_bin")):
            raise ValueError(f"{path}: missing AF/AN stratum labels")
        r.update(chrom=c, start0=lo, end0=hi, threshold=t, n_sites=n, n_multi=o, n_three=third)
        key = counter_key(r)
        if key in keys:
            raise ValueError(f"{path}: duplicate stratum in bin {c}:{lo}-{hi}: {key}")
        keys.add(key)
        # Historical scanners retained IUPAC flanks without complementing them.
        # Their exact ambiguity cannot be recovered from counts; retain every
        # count and pool ambiguous flanks as N, including existing N strata.
        r["context"] = r["context"].translate(AMBIGUOUS_TO_N)
        yield r


def read_counters(files, native, thresholds, scopes):
    """Merge already sorted files without loading their rows into memory."""
    import heapq
    key = lambda r: (int(r["chrom"]), r["start0"], r["end0"], *counter_key(r))
    # Each input's strata are also sorted by their key in the existing scanner.
    # A local sort within one native bin makes alternate header-preserving
    # exporters legal while keeping memory bounded by one bin per source file.
    def sorted_file(path):
        stream = read_counter_file(path, native, thresholds, scopes)
        for _, group in itertools.groupby(stream, key=lambda r: (r["chrom"], r["start0"], r["end0"])):
            pooled = {}
            for row in group:
                k = counter_key(row)
                if k in pooled:
                    for field in ("n_sites", "n_multi", "n_three"):
                        pooled[k][field] += row[field]
                else:
                    pooled[k] = dict(row)
            yield from sorted(pooled.values(), key=key)
    previous = None
    for row in heapq.merge(*(sorted_file(p) for p in files), key=key):
        current = key(row)
        if current == previous:
            raise ValueError(f"Duplicated source stratum across files: {current}")
        previous = current
        yield row


class Rates:
    def __init__(self, minimum):
        self.minimum = minimum
        self.total = defaultdict(lambda: [0, 0])
        self.by_chrom = defaultdict(lambda: [0, 0])
        self.support = defaultdict(set)
        self.cache = {}

    def observe(self, row):
        k, c = counter_key(row), row["chrom"]
        for dest in (self.total[k], self.by_chrom[(c, k)]):
            dest[0] += row["n_sites"]
            dest[1] += row["n_multi"]
        if row["n_sites"]:
            self.support[k].add(c)

    def predict(self, chrom, key):
        ck = (chrom, key)
        if ck not in self.cache:
            n, m = (self.total[key][i]-self.by_chrom[ck][i] for i in (0, 1))
            if n < self.minimum or len(self.support[key]-{chrom}) < 2:
                self.cache[ck] = None
            else:
                p = (m+.5)/(n+1.)
                self.cache[ck] = (p, p*(1-p)/(n+2.), n)
        return self.cache[ck]


def empty_window(c, lo, hi, width, step, threshold, primary, mode, bins_required):
    return dict(window_id=f"gw_chr{c}_{lo}_{hi}_w{width}_ac{threshold}", chrom="chr"+c,
                start0=lo, end0=hi, build="GRCh38", window_bp=width, step_bp=step, threshold=threshold,
                is_primary=int(threshold==primary), raw_n_eligible=0, raw_observed_multi=0, raw_three_alt_sites=0,
                n_eligible=0, observed_multi=0, expected_multi=0., n_unsupported=0, unsupported_observed_multi=0,
                source_bins_present=0, source_bins_required=bins_required, strata_count=0,
                training_exact_sites=0, training_context_change_sites=0, training_threshold_sites=0,
                background_rate_variance=0. if mode=="native_counters" else None,
                expected_model="exact_context_change_AF_AN_LOCO" if mode=="native_counters" else "legacy_hierarchical_LOCO_expected_preserved",
                source_mode=mode, denominator=DENOMINATOR,
                mask_overlap_bp=0, mask_overlap_fraction=0., mask_status="not_assessed")


def finalize_window(row):
    n, e = row["raw_n_eligible"], row["expected_multi"]
    row["n_unsupported"] = n-row["n_eligible"]
    row["unsupported_observed_multi"] = row["raw_observed_multi"]-row["observed_multi"]
    row["background_support_fraction"] = row["n_eligible"]/n if n else None
    row["spatial_support_fraction"] = row["source_bins_present"]/row["source_bins_required"] if row["source_bins_required"] else None
    row["expected_status"] = "no_support" if not row["n_eligible"] else "full" if row["n_unsupported"]==0 else "partial"
    if not row["n_eligible"]:
        row["expected_multi"] = e = None
    row["oe"] = row["observed_multi"]/e if e else None
    row["full_scope_oe"] = row["raw_observed_multi"]/e if e and not row["n_unsupported"] and row["spatial_support_fraction"]==1 else None
    return row


def from_counters(files, native, scopes, args):
    rates = Rates(args.min_training_sites)
    audit = Counter()
    clock = time.monotonic()
    for r in read_counters(files, native, args.threshold_values, scopes):
        rates.observe(r)
        audit["source_strata_rows"] += 1
        for name, field_name in (("eligible_sites", "n_sites"), ("multi_sites", "n_multi"), ("three_ALT_sites", "n_three")):
            audit[f"AC{r['threshold']}:{name}"] += r[field_name]
        if time.monotonic()-clock > 60:
            print(f"[genome-windows] background: {audit['source_strata_rows']:,} compact strata", flush=True)
            clock = time.monotonic()
    print(f"[genome-windows] background complete: {audit['source_strata_rows']:,} compact strata", flush=True)
    groups = iter(itertools.groupby(read_counters(files, native, args.threshold_values, scopes),
                                   key=lambda r: (r["chrom"], r["start0"], r["end0"])))
    pending = next(groups, None)
    result = []
    for c in sorted(scopes, key=int):
        for scope_start, scope_end in scopes[c]:
            if scope_start % native:
                raise ValueError("Scan scopes must start on a native-bin boundary; interval aggregates cannot be clipped")
            active = []
            def flush(end):
                nonlocal active
                remaining = []
                for row, counts in active:
                    if row["end0"] <= end:
                        row["strata_count"] = len(counts)
                        row["background_rate_variance"] = math.fsum(n*n*v for n, v in counts.values() if v is not None)
                        result.append(finalize_window(row))
                    else:
                        remaining.append((row, counts))
                active = remaining
            for lo in range(scope_start, scope_end, native):
                flush(lo)
                hi = min(lo+native, scope_end)
                if lo % args.step_bp == 0:
                    for width in args.size_values:
                        end = min(lo+width, scope_end)
                        # Partial terminal windows have complete source coverage
                        # but no equal-size calibration; inference marks them.
                        for threshold in args.threshold_values:
                            active.append((empty_window(c, lo, end, width, args.step_bp, threshold,
                                                        args.primary_threshold, "native_counters", math.ceil((end-lo)/native)), {}))
                            if len(result)+len(active)>args.max_windows:
                                raise ValueError("--max-windows exceeded; predeclare fewer sizes or a larger resource limit")
                for row, _ in active:
                    row["source_bins_present"] += 1
                if pending is not None and pending[0] == (c, lo, hi):
                    tile_counts = defaultdict(lambda: [0, 0, 0])
                    for r in pending[1]:
                        for i, name in enumerate(("n_sites", "n_multi", "n_three")):
                            tile_counts[r["threshold"]][i] += r[name]
                        key = counter_key(r)
                        prediction = rates.predict(c, key)
                        for row, counts in active:
                            if row["threshold"] != r["threshold"]:
                                continue
                            row["raw_n_eligible"] += r["n_sites"]
                            row["raw_observed_multi"] += r["n_multi"]
                            row["raw_three_alt_sites"] += r["n_three"]
                            existing = counts.setdefault(key, [0, prediction[1] if prediction else None])
                            existing[0] += r["n_sites"]
                            if prediction is not None:
                                row["n_eligible"] += r["n_sites"]
                                row["observed_multi"] += r["n_multi"]
                                row["expected_multi"] += r["n_sites"]*prediction[0]
                                row["training_exact_sites"] += r["n_sites"]
                    previous_counts = None
                    for threshold in args.threshold_values:
                        counts = tile_counts[threshold]
                        if counts[0] > hi-lo:
                            raise ValueError(f"Source bin chr{c}:{lo}-{hi} has more eligible sites than bases")
                        if previous_counts is not None and any(a>b for a,b in zip(counts, previous_counts)):
                            raise ValueError(f"Source bin chr{c}:{lo}-{hi}: counts increase with a stricter AC threshold")
                        previous_counts = counts
                    pending = next(groups, None)
                elif pending is not None and (int(pending[0][0]), pending[0][1]) < (int(c), lo):
                    raise ValueError("Source bin sequence does not match the declared scope")
            flush(scope_end)
    if pending is not None:
        raise ValueError("Unconsumed source bins outside declared scan scope")
    return result, audit


def summary_native(path, args):
    required = {"chrom", "start0", "end0", "threshold", "n_sites", "n_multi", "n_three", "expected_status", "expected_multi",
                "scored_sites", "scored_observed_multi", "expected_coverage", "unscored_sites",
                "training_exact_sites", "training_context_change_sites", "training_threshold_sites"}
    bins, widths, audit = {}, set(), Counter()
    for raw in table(path, required):
        c = chromosome(raw["chrom"])
        t = integer(raw["threshold"], "summary threshold", 1)
        if c not in AUTOSOMES or t not in args.threshold_values:
            continue
        lo, hi = integer(raw["start0"], "start0"), integer(raw["end0"], "end0", 1)
        if hi<=lo:
            raise ValueError("Legacy summary has a nonpositive bin width")
        widths.add(hi-lo)
        key = c, lo, t
        if key in bins:
            raise ValueError(f"Duplicate summary bin: {key}")
        n, o, third, sn, so, un = (integer(raw[k], k) for k in
                                  ("n_sites", "n_multi", "n_three", "scored_sites", "scored_observed_multi", "unscored_sites"))
        if not (0<=third<=o<=n<=hi-lo and sn+un==n and 0<=so<=sn and so<=o and o-so<=un):
            raise ValueError("Legacy summary count/support arithmetic is inconsistent")
        e = number(raw["expected_multi"], "expected_multi")
        coverage = number(raw["expected_coverage"], "expected_coverage")
        status = raw["expected_status"]
        wanted = "complete" if sn and un==0 else "partial" if sn else "unavailable"
        if status != wanted or n and (coverage is None or not math.isclose(coverage, sn/n, rel_tol=1e-8, abs_tol=1e-12)):
            raise ValueError("Legacy summary expected status/coverage disagrees with counts")
        if sn and (e is None or not 0<=e<=sn) or not sn and e not in (None, 0):
            raise ValueError("Legacy expected value must describe exactly the supported subset")
        for name, observed in (("scored_observed_over_expected", so), ("observed_over_expected", o)):
            if name in raw:
                ratio = number(raw[name], name)
                if name=="observed_over_expected" and un and ratio is not None:
                    raise ValueError("Partial legacy support must not have a full-scope O/E")
                if ratio is not None and (not e or not math.isclose(ratio, observed/e, rel_tol=1e-8, abs_tol=1e-12)):
                    raise ValueError("Legacy O/E disagrees with source counts/expected")
        row = dict(c=c, start0=lo, end0=hi, threshold=t, raw_n_eligible=n, raw_observed_multi=o,
                   raw_three_alt_sites=third, n_eligible=sn, observed_multi=so, expected_multi=e or 0.)
        for name in ("training_exact_sites", "training_context_change_sites", "training_threshold_sites"):
            row[name] = integer(raw.get(name, "0"), name)
        if all(k in raw for k in ("training_exact_sites", "training_context_change_sites", "training_threshold_sites")):
            if sum(row[k] for k in ("training_exact_sites", "training_context_change_sites", "training_threshold_sites"))!=sn:
                raise ValueError("Legacy background support levels do not sum to scored_sites")
        bins[key] = row
        audit["source_summary_rows"] += 1
        audit[f"AC{t}:eligible_sites"] += n
        audit[f"AC{t}:multi_sites"] += o
        audit[f"AC{t}:three_ALT_sites"] += third
    if not bins or len(widths)!=1:
        raise ValueError("Legacy summary requires one explicit native grid width and nonempty selected autosomal rows")
    native = widths.pop()
    if any(lo % native for _, lo, _ in bins):
        raise ValueError("Legacy bins are not aligned to their native width")
    return bins, native, audit


def from_summary(path, args):
    bins, native, audit = summary_native(path, args)
    validate_grid(native, args)
    result = []
    chromosomes = sorted({k[0] for k in bins}, key=int)
    fields = ("raw_n_eligible", "raw_observed_multi", "raw_three_alt_sites", "n_eligible", "observed_multi", "expected_multi",
              "training_exact_sites", "training_context_change_sites", "training_threshold_sites")
    for c in chromosomes:
        starts = [k[1] for k in bins if k[0]==c]
        stop = max(starts)+native
        for lo in range(min(starts)//args.step_bp*args.step_bp, stop, args.step_bp):
            for width in args.size_values:
                hi = lo+width
                for threshold in args.threshold_values:
                    row = empty_window(c, lo, hi, width, args.step_bp, threshold, args.primary_threshold,
                                       "legacy_summary", width//native)
                    row["strata_count"] = None
                    for start in range(lo, hi, native):
                        source = bins.get((c, start, threshold))
                        if source is not None:
                            row["source_bins_present"] += 1
                            for field_name in fields:
                                row[field_name] += source[field_name]
                    if row["source_bins_present"]:
                        result.append(finalize_window(row))
                        if len(result)>args.max_windows:
                            raise ValueError("--max-windows exceeded")
    return result, native, audit


def validate_grid(native, args):
    if native<1 or args.step_bp % native or any(w % native for w in args.size_values):
        raise ValueError(f"Window sizes and step must be integer multiples of source bin width {native}; finer resolution requires genuine finer counters")
    if any(w<args.step_bp for w in args.size_values):
        raise ValueError("Step cannot exceed any window size; otherwise the declared genome scan has avoidable gaps")
    if max(args.size_values)//args.step_bp>100:
        raise ValueError("At most 100 overlapping starts per largest window are supported")


def model_eligibility(row, args):
    if row["spatial_support_fraction"] is None or row["spatial_support_fraction"] < 1:
        return "incomplete_source_bin_coverage"
    if row["end0"]-row["start0"] != row["window_bp"]:
        return "truncated_scope_edge_no_equal_size_calibration"
    if row["n_eligible"] < args.min_sites:
        return "too_few_eligible_sites"
    if row["expected_multi"] is None or row["expected_multi"] < args.min_expected:
        return "too_little_expected_support"
    if row["background_support_fraction"] is None or row["background_support_fraction"] < args.min_support:
        return "insufficient_background_support"
    return ""


def dispersion_calibration(rows, args):
    # Eligibility is fixed by denominator/expected/support and external masks,
    # never by an observed count threshold, ratio, P value or nomination.
    moments = defaultdict(lambda: np.zeros(4))
    for row in rows:
        if model_eligibility(row, args) or row["start0"] % row["window_bp"] or row["mask_overlap_bp"]:
            continue
        o, e = row["observed_multi"], row["expected_multi"]
        key = (row["window_bp"], row["threshold"], chromosome(row["chrom"]))
        moments[key] += [1., (o-e)**2-e, e*e, o-e]
    calibrations, out = {}, []
    for width, threshold, chrom in sorted({(r["window_bp"], r["threshold"], chromosome(r["chrom"])) for r in rows}):
        selected = [v for (w, t, c), v in moments.items() if w==width and t==threshold and c!=chrom and v[0]>0]
        totals = np.sum(selected, axis=0) if selected else np.zeros(4)
        n, nchrom = int(totals[0]), len(selected)
        ok = n>=args.min_calibration_windows and nchrom>=2 and totals[2]>0
        alpha = max(0., float(totals[1]/totals[2])) if ok else None
        record = dict(chrom="chr"+chrom, window_bp=width, threshold=threshold,
                      reference_windows=n, reference_chromosomes=nchrom, nb2_alpha=alpha,
                      minimum_reference_windows=args.min_calibration_windows,
                      status="working_model_calibrated" if ok else "insufficient_reference_windows",
                      reference_mean_residual=float(totals[3]/n) if n else None,
                      method="NB2_moment_alpha_over_all_eligible_nonoverlapping_other_chromosome_windows",
                      uncertainty="fitted_working_model_not_LD_or_independent_cohort_validation")
        calibrations[(width, threshold, chrom)] = record
        out.append(record)
    return calibrations, out


def tail_interval(observed, expected, alpha):
    """Working tail plus a numerically validated conditional mean-ratio CI.

    scipy's beta PPF can return a finite endpoint near its numeric boundary
    even when the requested quantile is far below representable probabilities.
    Finiteness alone is therefore insufficient. Failed or unrepresentable
    intervals return (p, None, None); valid p values still enter the BY family.
    """
    if observed<0 or int(observed)!=observed or not math.isfinite(expected) or expected<=0 or not math.isfinite(alpha) or alpha<0:
        raise ValueError("Invalid count/expected/dispersion for a tail probability")
    if alpha<1e-12:
        p = stats.poisson.sf(observed-1, expected)
        lo = stats.gamma.ppf(.025, observed) if observed else 0.
        hi = stats.gamma.ppf(.975, observed+1)
        lower_tail = lambda mean: stats.poisson.sf(observed-1, mean)
        upper_tail = lambda mean: stats.poisson.cdf(observed, mean)
    else:
        shape = 1./alpha
        p = stats.nbinom.sf(observed-1, shape, shape/(shape+expected))
        prob_hi = stats.beta.ppf(.975, shape, observed) if observed else 1.
        prob_lo = stats.beta.ppf(.025, shape, observed+1)
        # Multiplication before division avoids reciprocal overflow when a
        # small shape can still yield a representable endpoint mean.
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            lo = shape*(1.-prob_hi)/prob_hi if observed else 0.
            hi = shape*(1.-prob_lo)/prob_lo
        lower_tail = lambda mean: stats.nbinom.sf(observed-1, shape, shape/(shape+mean))
        upper_tail = lambda mean: stats.nbinom.cdf(observed, shape, shape/(shape+mean))
    if not math.isfinite(p) or not 0<=p<=1:
        raise ValueError("Working-model tail probability failed numerically")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        ratio_lo, ratio_hi = lo/expected, hi/expected
    valid = all(math.isfinite(x) for x in (lo, hi, ratio_lo, ratio_hi)) and 0<=lo<=hi and 0<=ratio_lo<=ratio_hi
    if valid:
        # Check the actual means reconstructed from the published ratios too;
        # conversion underflow/overflow must not change the claimed interval.
        mean_lo, mean_hi = float(ratio_lo)*expected, float(ratio_hi)*expected
        valid = all(math.isfinite(x) for x in (mean_lo, mean_hi))
        if valid:
            valid = (not observed or math.isclose(lower_tail(mean_lo), .025, rel_tol=0., abs_tol=1e-8)) and math.isclose(upper_tail(mean_hi), .025, rel_tol=0., abs_tol=1e-8)
    if not valid:
        return float(p), None, None
    return float(p), float(ratio_lo), float(ratio_hi)


def by_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if not len(p):
        return []
    if np.any(~np.isfinite(p)) or np.any((p<0)|(p>1)):
        raise ValueError("Invalid P values for BY adjustment")
    n = len(p)
    order = np.argsort(p, kind="stable")
    harmonic = math.fsum(1./i for i in range(1, n+1))
    adjusted = np.minimum.accumulate((p[order]*n*harmonic/np.arange(1, n+1))[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.clip(adjusted, 0, 1)
    return result.tolist()


def infer(rows, args, masks):
    mask_index = masks if isinstance(masks, IntervalOverlapIndex) else IntervalOverlapIndex(masks)
    for row in rows:
        bp = overlap_bp(mask_index, row["chrom"], row["start0"], row["end0"])
        row.update(mask_overlap_bp=bp, mask_overlap_fraction=bp/(row["end0"]-row["start0"]),
                   mask_status="overlapping_exclusion" if bp else "no_overlap_with_supplied_masks" if args.exclude_bed else "not_assessed")
    calibrations, cal_rows = dispersion_calibration(rows, args)
    tested = []
    for row in rows:
        cal = calibrations[(row["window_bp"], row["threshold"], chromosome(row["chrom"]))]
        row.update(calibration_windows=cal["reference_windows"], calibration_chromosomes=cal["reference_chromosomes"],
                   nb2_alpha=cal["nb2_alpha"], working_nb2_alpha=None, ci95_status="not_estimable", ci95_low=None, ci95_high=None,
                   p_value=None, q_value=None, fdr_method="BY_all_estimable_sizes_and_AC_tiers",
                   effect_gate_pass=int(row["oe"] is not None and row["oe"]>=args.min_oe),
                   count_gate_pass=int(row["observed_multi"]>=args.min_observed), strict_2x_supported=0, nominated=0)
        reason = model_eligibility(row, args)
        if not reason and cal["nb2_alpha"] is None:
            reason = "insufficient_reference_windows"
        prefix = "legacy_summary_" if row["source_mode"]=="legacy_summary" else ""
        if reason:
            row.update(inference_status=prefix+reason, nomination_reason=reason)
            continue
        # Full counters expose posterior rate uncertainty. The added component
        # is conservative and may overlap residual dispersion; no raw input is
        # used to pretend exact background uncertainty is known in summary mode.
        alpha = cal["nb2_alpha"] + (row["background_rate_variance"] or 0.)/row["expected_multi"]**2
        pv, lo, hi = tail_interval(row["observed_multi"], row["expected_multi"], alpha)
        interval_valid = lo is not None and hi is not None
        row.update(p_value=pv, ci95_low=lo, ci95_high=hi, working_nb2_alpha=alpha,
                   ci95_status="validated_conditional_interval" if interval_valid else "not_representable_or_numerically_unavailable",
                   inference_status=prefix+"working_model_screen"+("" if interval_valid else "_CI_unavailable"),
                   strict_2x_supported=int(interval_valid and lo>=args.min_oe))
        tested.append(row)
    for row, q in zip(tested, by_adjust([r["p_value"] for r in tested])):
        row["q_value"] = q
        failed = []
        if not row["effect_gate_pass"]: failed.append("oe_below_gate")
        if not row["count_gate_pass"]: failed.append("observed_count_below_gate")
        if q>args.fdr: failed.append("BY_q_above_gate")
        if row["ci95_low"] is None or row["ci95_high"] is None: failed.append("conditional_CI_unavailable")
        elif row["ci95_low"]<=1: failed.append("conditional_CI_includes_one")
        if row["mask_overlap_bp"]: failed.append("whole_window_mask_veto")
        row["nominated"] = int(not failed)
        row["nomination_reason"] = ";".join(failed) if failed else "legacy_summary_candidate_requires_validation" if row["source_mode"]=="legacy_summary" else "model_based_candidate_requires_validation"
    return cal_rows, len(tested)


def merge_loci(rows):
    nominated = sorted((r for r in rows if r["nominated"]), key=lambda r: (r["threshold"], int(chromosome(r["chrom"])), r["start0"], r["end0"]))
    groups, active = [], []
    for row in nominated:
        if active and (row["threshold"]!=active[0]["threshold"] or row["chrom"]!=active[0]["chrom"] or row["start0"]>=max(v["end0"] for v in active)):
            groups.append(active); active=[]
        active.append(row)
    if active: groups.append(active)
    loci, members = [], []
    for index, values in enumerate(groups, 1):
        lead = min(values, key=lambda r: (r["q_value"], -r["oe"], r["window_bp"], r["start0"]))
        label = f"GW_AC{lead['threshold']}_{index:06d}"
        loci.append(dict(locus_id=label, chrom=lead["chrom"], start0=min(r["start0"] for r in values),
                         end0=max(r["end0"] for r in values), threshold=lead["threshold"], is_primary=lead["is_primary"],
                         member_windows=len(values), window_sizes=",".join(map(str, sorted({r["window_bp"] for r in values}))),
                         lead_window_id=lead["window_id"], lead_window_oe=lead["oe"], minimum_member_q=lead["q_value"],
                         inference="union_of_overlapping_nominated_windows_not_a_new_interval_test"))
        members.extend(dict(locus_id=label, window_id=r["window_id"], window_bp=r["window_bp"], threshold=r["threshold"],
                            is_lead=int(r is lead)) for r in values)
    return loci, members


def write_table(path, rows, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if row.get(k) is None else row.get(k, "") for k in fields})


def protected_sources(raw):
    result = [InputPath(p) for p in (raw.get("tiles") or [])+(raw.get("exclude_bed") or [])]
    if raw.get("store"):
        from scan_store import load_run_metadata
        run = load_run_metadata(raw["store"])
        root = Path(raw["store"]).expanduser().resolve()
        result.append(InputPath(root/"run.json"))
        for name in ("relocation.json", "adoption_receipt.json", "input_locations.json"):
            if (root/name).is_file():
                result.append(InputPath(root/name))
        for entry in run["shards"]:
            shard = Path(entry)
            result.append(InputPath(shard, "directory"))
            # File-level identities are also placed in the standard top-level
            # inputs receipt so the batch dispatcher can reject stale reuse
            # without reopening genomic data or relying on directory mtimes.
            result.extend(InputPath(shard/name) for name in
                          ("strata.tsv.gz", "scan_scope.tsv", "COMPLETE.json"))
            if (shard/"scan.json").is_file():
                result.append(InputPath(shard/"scan.json"))
    return result


def validate_args(args):
    if args.build not in {"GRCh38", "hg38", "38"}:
        raise ValueError("Explicit GRCh38 build declaration is required")
    args.size_values = ints(args.window_sizes, "--window-sizes")
    args.threshold_values = ints(args.thresholds, "--thresholds")
    if not set(args.threshold_values)<={5, 10} or args.primary_threshold not in args.threshold_values:
        raise ValueError("Use AC5/AC10 tiers and include --primary-threshold among --thresholds")
    if args.step_bp<1 or args.min_sites<1 or args.min_observed<1 or not math.isfinite(args.min_expected) or args.min_expected<=0 or args.min_training_sites<1:
        raise ValueError("Step, count/support minimums must be positive")
    if not 0<args.min_support<=1 or not 0<args.fdr<1 or not math.isfinite(args.min_oe) or args.min_oe<2:
        raise ValueError("Require 0<support<=1, 0<FDR<1 and a predeclared O/E gate >=2")
    if args.min_calibration_windows<20 or args.max_windows<1:
        raise ValueError("Need >=20 calibration windows (default 50) and a positive output-window limit")
    if args.tiles and not args.scope:
        raise ValueError("--tiles requires --scope with the complete source scan intervals; missing bins are otherwise unknown")
    if not args.tiles and args.scope:
        raise ValueError("--scope is only for explicit --tiles input; store scopes and legacy summary bins are retained")


@transactional("manifest.json", input_dirs=("store",), input_files=("summary", "scope"),
               input_prefixes=(), resolve_inputs=protected_sources)
def run(args):
    validate_args(args)
    out = Path(args.out)
    sources = {}
    if args.summary:
        sources["summary"] = identity(args.summary, full=True)
        windows, native, audit = from_summary(args.summary, args)
        mode = "legacy_summary_screen_only"
        source_metadata = dict(historical_expected="retained; exact/context-change/threshold fallback counts preserved",
                               original_scan_signature="not_available_from_summary_alone", strict_base_mask="not_reconstructed")
    else:
        if args.store:
            import scan
            run_record = scan.load_run(args.store)
            options = run_record["options"]
            if not set(args.threshold_values) <= set(map(int, options.get("thresholds", []))):
                raise ValueError("Requested AC tiers were not recorded by the completed source scan")
            build = str(options.get("architecture_build", options.get("build", "38"))).lower()
            if build not in {"38", "grch38", "hg38"}:
                raise ValueError("Source scan must declare GRCh38")
            files = [Path(s)/"strata.tsv.gz" for s in run_record["shards"]]
            scope_files = [Path(s)/"scan_scope.tsv" for s in run_record["shards"]]
            native = integer(options.get("window_size", args.tile_bp), "source window_size", 1)
            source_metadata = dict(original_scan_signature=run_record["signature"], source_options=options,
                                   strict_base_mask="original_scan_filtering_only; no HERV-prefilter or retrospective strict mask")
            sources["scan_manifest"] = identity(Path(args.store)/"run.json", full=True)
        else:
            files = [Path(p).resolve() for p in args.tiles]
            scope_files = [Path(args.scope).resolve()]
            native = args.tile_bp
            source_metadata = dict(original_scan_signature="caller_supplied_counters_not_independently_verified",
                                   strict_base_mask="as_encoded_by_counter_source")
        if not files or len(set(files))!=len(files) or len(files)>128:
            raise ValueError("Require 1..128 unique compact counter files; merge compatible source bins explicitly before larger fan-in")
        for path in (*files, *scope_files):
            if not path.is_file():
                raise ValueError(f"Required completed cache missing: {path}. No raw VCF rescan is started.")
        validate_grid(native, args)
        scopes = intervals(scope_files, allow_overlap=False)
        if len(scopes)<3:
            raise ValueError("At least three autosomal source chromosomes are required for exact LOCO support")
        sources["counters"] = [identity(p) for p in files]
        sources["scopes"] = [identity(p, full=True) for p in scope_files]
        windows, audit = from_counters(files, native, scopes, args)
        if sources["counters"] != [identity(p) for p in files]:
            raise ValueError("Compact input counters changed between passes")
        mode = "native_counter_scan"
    windows.sort(key=lambda r: (int(chromosome(r["chrom"])), r["start0"], r["window_bp"], r["threshold"]))
    masks = intervals(args.exclude_bed, allow_overlap=True) if args.exclude_bed else {}
    sources["exclude_bed"] = [identity(p, full=True) for p in args.exclude_bed]
    calibrations, tested = infer(windows, args, masks)
    nominations = [r for r in windows if r["nominated"]]
    loci, members = merge_loci(windows)
    write_table(out/"all_windows.tsv.gz", windows, WINDOW_FIELDS)
    write_table(out/"nominated_windows.tsv", nominations, WINDOW_FIELDS)
    write_table(out/"calibration.tsv", calibrations, list(calibrations[0]) if calibrations else
                ["chrom", "window_bp", "threshold", "status"])
    locus_fields = ["locus_id", "chrom", "start0", "end0", "threshold", "is_primary", "member_windows", "window_sizes",
                    "lead_window_id", "lead_window_oe", "minimum_member_q", "inference"]
    write_table(out/"loci.tsv", loci, locus_fields)
    write_table(out/"locus_members.tsv", members, ["locus_id", "window_id", "window_bp", "threshold", "is_lead"])
    with (out/"loci.bed").open("w") as handle:
        for r in loci:
            handle.write(f"{r['chrom']}\t{r['start0']}\t{r['end0']}\t{r['locus_id']}\n")
    audit.update(output_windows=len(windows), estimable_tests=tested, nominated_windows=len(nominations),
                 primary_nominated_windows=sum(r["is_primary"] for r in nominations), merged_loci=len(loci),
                 conditional_interval_unavailable_windows=sum(r["p_value"] is not None and r["ci95_low"] is None for r in windows))
    write_table(out/"scan_audit.tsv", [dict(metric=k, value=v) for k,v in sorted(audit.items())], ["metric", "value"])
    status = "complete" if tested else "complete_with_exclusions"
    return dict(status=status, schema=VERSION, mode=mode, build="GRCh38", source_bin_bp=native,
                source_metadata=source_metadata, source_identity=sources, counts=dict(audit),
                parameters={k: getattr(args, k) for k in ("window_sizes", "step_bp", "thresholds", "primary_threshold", "min_sites",
                    "min_observed", "min_expected", "min_support", "min_training_sites", "min_calibration_windows", "min_oe", "fdr", "max_windows")},
                testing_family="all estimable autosomal windows across every requested size and AC tier before observed/OE/mask nomination gates",
                confidence_interval="conditional equal-tail mean-ratio interval under fitted NB2/Poisson working distribution; not simultaneous or full uncertainty",
                calibration="all eligible unselected nonoverlapping reference windows; target chromosome omitted from residual pool",
                candidate_rule="OE>=min_oe, supported O>=min_observed, all denominator/background/spatial/calibration gates, BY q<=fdr, numerically validated finite conditional CI with lower>1, no exclusion-BED overlap",
                strict_2x_supported="additional flag: conditional CI lower >= configured min_oe; does not remove modelling assumptions",
                denominator=DENOMINATOR, source_filter_change="none; no raw source modification or VCF rescan",
                limitations=LIMITATIONS)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--store", help="Completed scan root; validates original manifests and compact counter shards")
    source.add_argument("--tiles", nargs="+", help="Explicit strata TSV(.gz) files with the original counter schema")
    source.add_argument("--summary", help="Existing window_summary.tsv; preserves historical expected/fallback model")
    p.add_argument("--scope", help="Complete BED0/scan_scope.tsv coverage for explicit --tiles")
    p.add_argument("--tile-bp", type=int, default=1000000, help="Native source width for --tiles; never invented by splitting counts")
    p.add_argument("--build", required=True, choices=("GRCh38", "hg38", "38"))
    p.add_argument("--window-sizes", default="1000000", help="Prespecified comma-separated widths, e.g. 1000000,2000000,5000000")
    p.add_argument("--step-bp", type=int, default=1000000)
    p.add_argument("--thresholds", default="5,10")
    p.add_argument("--primary-threshold", type=int, default=5, choices=(5,10))
    p.add_argument("--min-sites", type=int, default=500)
    p.add_argument("--min-observed", type=int, default=25)
    p.add_argument("--min-expected", type=float, default=10.)
    p.add_argument("--min-support", type=float, default=.99)
    p.add_argument("--min-training-sites", type=int, default=100)
    p.add_argument("--min-calibration-windows", type=int, default=50)
    p.add_argument("--min-oe", type=float, default=2.)
    p.add_argument("--fdr", type=float, default=.05)
    p.add_argument("--exclude-bed", action="append", default=[], help="BED0 whole-window nomination veto; repeatable; counts are never prorated")
    p.add_argument("--max-windows", type=int, default=250000)
    p.add_argument("--out", required=True)
    return p


def main(argv=None):
    try:
        result = run(parser().parse_args(argv))
    except (ValueError, OSError, RuntimeError) as exc:
        raise SystemExit("genome-windows: "+str(exc))
    print(json.dumps(dict(status=result["status"], mode=result["mode"], counts=result["counts"]), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
