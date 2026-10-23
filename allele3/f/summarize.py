#!/usr/bin/env python3
"""Stream compact occupancy strata into descriptive, mutation-context summaries.

The expected values are chromosome-held-out *conditional observation* rates.
They are neither mutation rates nor a selection test. No network calls are made.
Only Python's standard library is required; matplotlib is an optional renderer.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable


REQUIRED = {
    "chrom", "start0", "end0", "threshold", "context", "primary_change",
    "primary_af_bin", "an_bin", "n_sites", "n_multi", "n_three",
}
LEVELS = ("exact", "context_change", "threshold")
AUTOSOMES = {str(i) for i in range(1, 23)}


def canonical_chrom(value):
    value = str(value).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value.isdigit():
        return str(int(value))
    return "MT" if value in {"M", "MT"} else value


def chrom_key(chrom):
    if chrom in AUTOSOMES:
        return (0, int(chrom))
    return (1, {"X": 0, "Y": 1, "MT": 2}.get(chrom, 3), chrom)


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith((".gz", ".bgz", ".bgzip")) else open(path, encoding="utf-8")


def parse_region(region):
    """Return (canonical chromosome, zero-based start, exclusive end).

    A string uses the usual VCF/samtools one-based inclusive chr:start-end form.
    A tuple/list is explicitly zero-based half-open. A chromosome-only string
    selects that chromosome without bounding its coordinate span.
    """
    if region is None or region == "":
        return None
    if isinstance(region, (tuple, list)) and len(region) == 3:
        chrom, start, end = region
        start, end = int(start), int(end)
    elif isinstance(region, str):
        if ":" not in region:
            return canonical_chrom(region), 0, None
        chrom, coordinates = region.rsplit(":", 1)
        if "-" not in coordinates:
            raise ValueError("region must be chr:start-end (one-based inclusive)")
        start, end = coordinates.replace(",", "").split("-", 1)
        start, end = int(start) - 1, int(end)
    else:
        raise ValueError("region must be a chr:start-end string or (chrom,start0,end0)")
    if start < 0 or end <= start:
        raise ValueError("region has invalid coordinates")
    return canonical_chrom(chrom), start, end


def merge_intervals(intervals):
    out = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError(f"invalid zero-based interval: {start}, {end}")
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out


def intersect_intervals(left, right):
    """Intersection of two sorted disjoint interval collections."""
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        if max(a, c) < min(b, d):
            result.append((max(a, c), min(b, d)))
        if b <= d:
            i += 1
        else:
            j += 1
    return result


def load_bed(path):
    intervals = defaultdict(list)
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: BED requires at least 3 columns")
            chrom = canonical_chrom(fields[0])
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid BED coordinates") from exc
            intervals[chrom].append((start, end))
    return {chrom: merge_intervals(rows) for chrom, rows in intervals.items()}


def _find_strata(shard_dirs):
    files = []
    for item in shard_dirs:
        item = Path(item)
        if item.is_file():
            path = item
        else:
            candidates = [item / "strata.tsv.gz", item / "strata.tsv"]
            path = next((p for p in candidates if p.is_file()), None)
            if path is None:
                raise FileNotFoundError(f"No strata.tsv.gz or strata.tsv in {item}")
        files.append(path)
    if not files:
        raise ValueError("at least one strata file is required")
    if len({str(p.resolve()) for p in files}) != len(files):
        raise ValueError("duplicate strata input file would double-count observations")
    return files


def read_strata(files, window_size, region):
    for path in files:
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = REQUIRED - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{path}: missing columns: {','.join(sorted(missing))}")
            for line_number, row in enumerate(reader, 2):
                try:
                    for column in ("start0", "end0", "threshold", "n_sites", "n_multi", "n_three"):
                        row[column] = int(row[column])
                except (ValueError, TypeError) as exc:
                    raise ValueError(f"{path}:{line_number}: count/coordinate must be an integer") from exc
                row["chrom"] = canonical_chrom(row["chrom"])
                start, end = row["start0"], row["end0"]
                if start < 0 or end <= start or start // window_size != (end - 1) // window_size:
                    raise ValueError(f"{path}:{line_number}: interval invalid or crosses a {window_size}-bp summary window")
                if not 0 <= row["n_three"] <= row["n_multi"] <= row["n_sites"]:
                    raise ValueError(f"{path}:{line_number}: require 0 <= n_three <= n_multi <= n_sites")
                if row["threshold"] < 1:
                    raise ValueError(f"{path}:{line_number}: threshold must be positive")
                if region:
                    rchrom, rstart, rend = region
                    if row["chrom"] != rchrom or end <= rstart or (rend is not None and start >= rend):
                        continue
                    # Counts cannot be clipped after aggregation. Root scan must
                    # already enforce region; bins may extend across its edge.
                row["window_start0"] = (start // window_size) * window_size
                yield row


def _count():
    return {"n_sites": 0, "n_multi": 0, "n_three": 0}


def _add_count(target, row):
    for field in ("n_sites", "n_multi", "n_three"):
        target[field] += row[field]


def _expected_count():
    return {"expected_multi": 0.0, "scored_sites": 0, "scored_observed_multi": 0,
            "exact_sites": 0, "context_change_sites": 0, "threshold_sites": 0,
            "unscored_sites": 0, "training_sites_min": None, "training_sites_max": None,
            "unavailable_reasons": set()}


def _add_expected(target, row, prediction):
    if prediction["rate"] is None:
        target["unscored_sites"] += row["n_sites"]
        target["unavailable_reasons"].add(prediction["reason"])
        return
    target["expected_multi"] += row["n_sites"] * prediction["rate"]
    target["scored_sites"] += row["n_sites"]
    target["scored_observed_multi"] += row["n_multi"]
    target[prediction["level"] + "_sites"] += row["n_sites"]
    ntrain = prediction["n_train"]
    target["training_sites_min"] = ntrain if target["training_sites_min"] is None else min(target["training_sites_min"], ntrain)
    target["training_sites_max"] = ntrain if target["training_sites_max"] is None else max(target["training_sites_max"], ntrain)


def _keys(row):
    exact = (row["threshold"], row["context"], row["primary_change"], row["primary_af_bin"], row["an_bin"])
    return exact, exact[:3], exact[:1]


class HeldOutRates:
    """Compact sufficient counts; target chromosome is subtracted at all levels."""

    def __init__(self, min_training_sites, restricted=False):
        self.min_training_sites = min_training_sites
        self.restricted = restricted
        self.total = defaultdict(lambda: [0, 0])
        self.by_chrom = defaultdict(lambda: [0, 0])
        self.support = defaultdict(set)
        self.autosomes = set()
        self.cache = {}

    def observe(self, row):
        chrom = row["chrom"]
        if chrom not in AUTOSOMES or row["n_sites"] == 0:
            return
        self.autosomes.add(chrom)
        for level, key in zip(LEVELS, _keys(row)):
            gkey = (level, key)
            ckey = (chrom, level, key)
            self.total[gkey][0] += row["n_sites"]
            self.total[gkey][1] += row["n_multi"]
            self.by_chrom[ckey][0] += row["n_sites"]
            self.by_chrom[ckey][1] += row["n_multi"]
            self.support[gkey].add(chrom)

    def predict(self, row):
        chrom = row["chrom"]
        cache_key = (chrom, _keys(row)[0])
        if cache_key in self.cache:
            return self.cache[cache_key]
        reason = None
        if chrom not in AUTOSOMES:
            reason = "non_autosome_descriptive_only"
        elif self.restricted:
            reason = "region_restricted_no_cross_chromosome_background"
        elif len(self.autosomes) < 2:
            reason = "fewer_than_two_represented_autosomes"
        if reason is None:
            for level, key in zip(LEVELS, _keys(row)):
                total = self.total.get((level, key), (0, 0))
                own = self.by_chrom.get((chrom, level, key), (0, 0))
                ntrain, multitrain = total[0] - own[0], total[1] - own[1]
                if ntrain >= self.min_training_sites:
                    result = {"rate": (multitrain + 0.5) / (ntrain + 1), "level": level,
                              "n_train": ntrain, "n_train_multi": multitrain,
                              "n_training_autosomes": len(self.support[(level, key)] - {chrom}), "reason": ""}
                    self.cache[cache_key] = result
                    return result
            reason = "insufficient_held_out_training_sites"
        result = {"rate": None, "level": "unavailable", "n_train": None,
                  "n_train_multi": None, "n_training_autosomes": 0, "reason": reason}
        self.cache[cache_key] = result
        return result


def wilson_interval(success, total, z=1.959963984540054):
    if total == 0:
        return None, None
    p = success / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    width = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    return max(0.0, center - width), min(1.0, center + width)


def _raw_fields(counts):
    n, multi, three = counts["n_sites"], counts["n_multi"], counts["n_three"]
    low, high = wilson_interval(multi, n)
    return {**counts, "n_single_alt": n - multi, "n_two_alt": multi - three,
            "multi_fraction": multi / n if n else None, "multi_wilson_low": low, "multi_wilson_high": high,
            "three_alt_fraction": three / n if n else None}


def _expected_fields(raw, expected):
    expected = expected or _expected_count()
    scored, unscored = expected["scored_sites"], expected["unscored_sites"]
    e = expected["expected_multi"] if scored else None
    complete = scored > 0 and unscored == 0
    return {
        "expected_status": "complete" if complete else "partial" if scored else "unavailable",
        "expected_multi": e,
        "observed_over_expected": raw["n_multi"] / e if complete and e else None,
        "scored_sites": scored,
        "scored_observed_multi": expected["scored_observed_multi"],
        "scored_observed_over_expected": expected["scored_observed_multi"] / e if e else None,
        "expected_coverage": scored / raw["n_sites"] if raw["n_sites"] else None,
        "training_exact_sites": expected["exact_sites"],
        "training_context_change_sites": expected["context_change_sites"],
        "training_threshold_sites": expected["threshold_sites"],
        "unscored_sites": unscored,
        "training_sites_min": expected["training_sites_min"],
        "training_sites_max": expected["training_sites_max"],
        "expected_unavailable_reason": ";".join(sorted(expected["unavailable_reasons"])),
    }


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else row.get(key, "") for key in fieldnames})


RAW_COLUMNS = ["n_sites", "n_multi", "n_three", "n_single_alt", "n_two_alt", "multi_fraction",
               "multi_wilson_low", "multi_wilson_high", "three_alt_fraction"]
EXPECTED_COLUMNS = ["expected_status", "expected_multi", "observed_over_expected", "scored_sites",
                    "scored_observed_multi", "scored_observed_over_expected", "expected_coverage",
                    "training_exact_sites", "training_context_change_sites", "training_threshold_sites",
                    "unscored_sites", "training_sites_min", "training_sites_max", "expected_unavailable_reason"]


def _load_scan_scopes(files):
    """Optional root-written scope records prevent unscanned intervals becoming zeros."""
    intervals = defaultdict(list)
    found = 0
    for directory in dict.fromkeys(p.parent for p in files):
        path = directory / "scan_scope.tsv"
        if not path.is_file():
            continue
        found += 1
        with open(path, encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"chrom", "start0", "end0"}.issubset(reader.fieldnames or []):
                raise ValueError(f"{path}: expected chrom,start0,end0 columns")
            for row in reader:
                intervals[canonical_chrom(row["chrom"])].append((int(row["start0"]), int(row["end0"])))
    total = len(set(p.parent for p in files))
    if found not in (0, total):
        raise ValueError("scan_scope.tsv must be supplied for every shard or for none")
    return ({chrom: merge_intervals(values) for chrom, values in intervals.items()}, bool(found))


def _callable_windows(path, scopes, has_scopes, seen_chroms, window_size, region, warnings):
    bed = load_bed(path)
    allowed_chroms = set(scopes) if has_scopes else set(seen_chroms)
    if region:
        allowed_chroms = {region[0]} & (allowed_chroms or {region[0]})
    if not has_scopes:
        warnings.append("No scan_scope.tsv: callable zero windows are limited to chromosomes represented in strata (or the explicit region). Root scan must have covered each such chromosome/region completely.")
    by_window = defaultdict(int)
    by_chrom = defaultdict(int)
    for chrom in sorted(allowed_chroms, key=chrom_key):
        intervals = bed.get(chrom, [])
        if has_scopes:
            intervals = intersect_intervals(intervals, scopes.get(chrom, []))
        if region:
            _, start, end = region
            limit = end if end is not None else max((b for _, b in intervals), default=start)
            intervals = intersect_intervals(intervals, [(start, limit)]) if limit > start else []
        for start, end in intervals:
            by_chrom[chrom] += end - start
            while start < end:
                wstart = (start // window_size) * window_size
                stop = min(end, wstart + window_size)
                by_window[(chrom, wstart)] += stop - start
                start = stop
    return by_window, by_chrom


def _make_plots(out_dir, chromosome_rows, window_rows, context_rows, threshold_rows):
    try:
        os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "allele3-matplotlib"))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return [], ["matplotlib is unavailable; all machine-readable tables were produced"]
    paths = []
    if not threshold_rows:
        return paths, ["No strata rows; no figures were rendered"]
    thresholds = [r["threshold"] for r in threshold_rows]
    selected = thresholds[0]
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

    def save(fig, name):
        fig.tight_layout()
        for extension in ("png", "pdf"):
            path = out_dir / f"{name}.{extension}"
            fig.savefig(path, dpi=180, bbox_inches="tight")
            paths.append(str(path))
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    chroms = sorted({r["chrom"] for r in chromosome_rows}, key=chrom_key)
    xmap = {chrom: i for i, chrom in enumerate(chroms)}
    for threshold in thresholds:
        rows = [r for r in chromosome_rows if r["threshold"] == threshold and r["multi_fraction"] is not None]
        ax.plot([xmap[r["chrom"]] for r in rows], [100 * r["multi_fraction"] for r in rows], marker="o", markersize=3, linewidth=1, label=f"AC ≥ {threshold}")
    ax.set_xticks(range(len(chroms)), chroms)
    ax.set(xlabel="Chromosome", ylabel="Sites with ≥2 observed ALT alleles (%)", title="Multi-ALT fraction among eligible SNV sites")
    ax.legend(frameon=False)
    save(fig, "chromosome_occupancy")

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True, gridspec_kw={"height_ratios": [1.2, 1]})
    rows = [r for r in window_rows if r["threshold"] == selected]
    offsets = {}
    cursor = 0
    ticks, ticklabels = [], []
    spans = [max(r["end0"] for r in rows if r["chrom"] == c) -
             min(r["start0"] for r in rows if r["chrom"] == c)
             for c in chroms if any(r["chrom"] == c for r in rows)]
    gap = max(1, min(1_000_000, int(max(spans, default=1) * .02)))
    for index, chrom in enumerate(chroms):
        subset = [r for r in rows if r["chrom"] == chrom]
        if not subset:
            continue
        low, high = min(r["start0"] for r in subset), max(r["end0"] for r in subset)
        offsets[chrom] = cursor - low
        ticks.append(cursor + (high - low) / 2)
        ticklabels.append(chrom)
        color = "#306FA2" if index % 2 == 0 else "#67A998"
        valid = [r for r in subset if r["multi_fraction"] is not None]
        top.scatter([r["start0"] + offsets[chrom] for r in valid], [100 * r["multi_fraction"] for r in valid], s=7, color=color)
        modeled = [r for r in subset if r["observed_over_expected"] is not None]
        bottom.scatter([r["start0"] + offsets[chrom] for r in modeled], [r["observed_over_expected"] for r in modeled], s=7, color=color)
        cursor += high - low + gap
    top.set(ylabel="Multi-ALT sites (%)", title=f"Conditional genome distribution · AC ≥ {selected}")
    bottom.axhline(1, color="#777777", linewidth=0.8, linestyle="--")
    bottom.set(ylabel="Observed / held-out expected", xlabel="Chromosome")
    bottom.set_xticks(ticks, ticklabels)
    if not any(r["observed_over_expected"] is not None for r in rows):
        bottom.text(0.5, 0.55, "Cross-chromosome expectation unavailable", transform=bottom.transAxes, ha="center")
    save(fig, "window_distribution")

    context_aggregates = defaultdict(_count)
    for row in context_rows:
        if row["threshold"] == selected:
            _add_count(context_aggregates[row["context"]], row)
    contexts = sorted(context_aggregates)
    fig, ax = plt.subplots(figsize=(max(8, len(contexts) * 0.28), 4.6))
    ax.bar(contexts, [100 * context_aggregates[c]["n_multi"] / context_aggregates[c]["n_sites"] if context_aggregates[c]["n_sites"] else 0 for c in contexts], color="#306FA2")
    ax.tick_params(axis="x", labelrotation=90)
    ax.set(ylabel="Multi-ALT sites (%)", xlabel="Pyrimidine-canonical trinucleotide", title=f"Context distribution · AC ≥ {selected}")
    save(fig, "context_occupancy")
    return paths, []


def summarize(shard_dirs, out_dir, window_size=1_000_000, min_training_sites=100,
              make_plots=True, callable_bed=None, region=None, thresholds=None):
    """Summarize compact strata in two passes without loading a full dataframe.

    Root scanning must already enforce the requested region and callable BED;
    aggregated counts cannot be filtered by individual site coordinates here.
    Optional per-shard scan_scope.tsv bounds zero-window materialization.
    """
    if window_size <= 0 or min_training_sites < 1:
        raise ValueError("window_size and min_training_sites must be positive")
    if isinstance(shard_dirs, (str, Path)):
        shard_dirs = [shard_dirs]
    files = _find_strata(shard_dirs)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    parsed_region = parse_region(region)
    warnings = []
    raw_chr, raw_window = defaultdict(_count), defaultdict(_count)
    raw_context, raw_threshold = defaultdict(_count), defaultdict(_count)
    if thresholds is not None:
        for threshold in thresholds:
            threshold = int(threshold)
            if threshold < 1:
                raise ValueError("thresholds must be positive integers")
            raw_threshold[threshold]
    rates = HeldOutRates(min_training_sites, restricted=parsed_region is not None)
    seen_chroms, exact_rows = set(), defaultdict(_count)
    n_strata_rows = 0
    for row in read_strata(files, window_size, parsed_region):
        n_strata_rows += 1
        chrom, threshold = row["chrom"], row["threshold"]
        seen_chroms.add(chrom)
        _add_count(raw_chr[(chrom, threshold)], row)
        _add_count(raw_window[(chrom, row["window_start0"], threshold)], row)
        _add_count(raw_context[(threshold, row["context"], row["primary_change"])], row)
        _add_count(raw_threshold[threshold], row)
        _add_count(exact_rows[(chrom, _keys(row)[0])], row)
        rates.observe(row)
    scopes, has_scopes = _load_scan_scopes(files)
    if has_scopes:
        for chrom in scopes:
            if parsed_region and chrom != parsed_region[0]:
                continue
            for threshold in raw_threshold:
                raw_chr[(chrom, threshold)]
    callable_window, callable_chr = {}, {}
    if callable_bed is not None:
        callable_window, callable_chr = _callable_windows(callable_bed, scopes, has_scopes, seen_chroms, window_size, parsed_region, warnings)
        # Every scanned callable window appears even if no variant passed QC.
        for (chrom, start), bases in callable_window.items():
            if bases:
                for threshold in raw_threshold:
                    raw_window[(chrom, start, threshold)]
                    raw_chr[(chrom, threshold)]
        for (chrom, start, threshold), counts in raw_window.items():
            if counts["n_sites"] > 0 and callable_window.get((chrom, start), 0) == 0:
                raise ValueError(f"observed sites in {chrom}:{start} have zero callable denominator; scan and summary masks/scopes differ")
            if counts["n_sites"] > callable_window.get((chrom, start), 0):
                raise ValueError(f"more eligible SNV sites than callable bases in {chrom}:{start}; check duplicate shards and scan/summary mask agreement")
    expected_chr, expected_window, expected_threshold = defaultdict(_expected_count), defaultdict(_expected_count), defaultdict(_expected_count)
    for row in read_strata(files, window_size, parsed_region):
        prediction = rates.predict(row)
        _add_expected(expected_chr[(row["chrom"], row["threshold"])], row, prediction)
        _add_expected(expected_window[(row["chrom"], row["window_start0"], row["threshold"])], row, prediction)
        _add_expected(expected_threshold[row["threshold"]], row, prediction)

    chromosome_rows, window_rows, context_rows, threshold_rows = [], [], [], []
    for (chrom, threshold), counts in sorted(raw_chr.items(), key=lambda item: (chrom_key(item[0][0]), item[0][1])):
        bases = callable_chr.get(chrom, 0) if callable_bed is not None else None
        chromosome_rows.append({"chrom": chrom, "threshold": threshold, **_raw_fields(counts),
                                "callable_bases": bases, "multi_per_callable_mb": counts["n_multi"] * 1e6 / bases if bases else None,
                                **_expected_fields(counts, expected_chr.get((chrom, threshold)))})
    for (chrom, start, threshold), counts in sorted(raw_window.items(), key=lambda item: (chrom_key(item[0][0]), item[0][1], item[0][2])):
        bases = callable_window.get((chrom, start), 0) if callable_bed is not None else None
        end = start + window_size
        # Coordinates label the complete bin; the callable denominator is clipped
        # to scan scope/region and may be smaller than the displayed bin width.
        window_rows.append({"chrom": chrom, "start0": start, "end0": end, "threshold": threshold,
                            **_raw_fields(counts), "callable_bases": bases,
                            "multi_per_callable_mb": counts["n_multi"] * 1e6 / bases if bases else None,
                            **_expected_fields(counts, expected_window.get((chrom, start, threshold)))})
    for (threshold, context, change), counts in sorted(raw_context.items()):
        context_rows.append({"threshold": threshold, "context": context, "primary_change": change, **_raw_fields(counts)})
    for threshold, counts in sorted(raw_threshold.items()):
        threshold_rows.append({"threshold": threshold, **_raw_fields(counts), **_expected_fields(counts, expected_threshold.get(threshold))})
    expected_context_rows = []
    for (chrom, key), counts in sorted(exact_rows.items(), key=lambda item: (chrom_key(item[0][0]), item[0][1])):
        threshold, context, change, af_bin, an_bin = key
        row = {"chrom": chrom, "threshold": threshold, "context": context, "primary_change": change,
               "primary_af_bin": af_bin, "an_bin": an_bin}
        prediction = rates.predict(row)
        expected_context_rows.append({**row, "observed_sites": counts["n_sites"], "observed_multi": counts["n_multi"],
                                     "training_level": prediction["level"], "training_sites": prediction["n_train"],
                                     "training_multi": prediction["n_train_multi"], "training_autosomes": prediction["n_training_autosomes"],
                                     "held_out_rate": prediction["rate"],
                                     "expected_multi": counts["n_sites"] * prediction["rate"] if prediction["rate"] is not None else None,
                                     "unavailable_reason": prediction["reason"]})
    specs = [
        ("chromosome_summary.tsv", chromosome_rows, ["chrom", "threshold"] + RAW_COLUMNS + ["callable_bases", "multi_per_callable_mb"] + EXPECTED_COLUMNS),
        ("window_summary.tsv", window_rows, ["chrom", "start0", "end0", "threshold"] + RAW_COLUMNS + ["callable_bases", "multi_per_callable_mb"] + EXPECTED_COLUMNS),
        ("context_summary.tsv", context_rows, ["threshold", "context", "primary_change"] + RAW_COLUMNS),
        ("occupancy_thresholds.tsv", threshold_rows, ["threshold"] + RAW_COLUMNS + EXPECTED_COLUMNS),
        ("expected_context.tsv", expected_context_rows, ["chrom", "threshold", "context", "primary_change", "primary_af_bin", "an_bin", "observed_sites", "observed_multi", "training_level", "training_sites", "training_multi", "training_autosomes", "held_out_rate", "expected_multi", "unavailable_reason"]),
    ]
    paths = []
    for filename, rows, columns in specs:
        path = out_dir / filename
        _write_tsv(path, rows, columns)
        paths.append(str(path))
    if make_plots:
        try:
            plot_paths, plot_warnings = _make_plots(out_dir, chromosome_rows, window_rows, context_rows, threshold_rows)
            paths.extend(plot_paths)
            warnings.extend(plot_warnings)
        except Exception as exc:
            warnings.append(f"Optional figure rendering failed; tables remain valid: {type(exc).__name__}: {exc}")
    manifest = {
        "status": "complete", "inputs": [str(p) for p in files], "outputs": paths,
        "n_strata_rows": n_strata_rows, "n_windows": len({(r["chrom"], r["start0"]) for r in window_rows}),
        "thresholds": sorted(raw_threshold), "represented_autosomes": sorted(rates.autosomes, key=chrom_key),
        "window_size": window_size, "min_training_sites": min_training_sites,
        "region": region, "callable_bed": str(callable_bed) if callable_bed is not None else None,
        "scan_scope_supplied": has_scopes,
        "expected_method": "leave-one-autosome-out conditional observation probability; exact context/change/AF-bin/AN-bin, then context/change, then threshold; (multi+0.5)/(sites+1)",
        "expected_interpretation": "Conditional on an eligible polymorphic SNV site. Not a mutation rate, retention estimate, significance test, or evidence of balancing selection.",
        "wilson_interpretation": "Descriptive binomial intervals only; neighboring sites are not assumed independent for hypothesis testing.",
        "zero_window_policy": "Explicit callable zero windows only within supplied scan scope, or represented chromosomes/explicit region when scan scope is absent.",
        "window_coordinates": "0-based half-open complete bins. Callable bases are clipped to scan scope and region.",
        "warnings": warnings,
    }
    manifest_path = out_dir / "summarize_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    manifest["manifest_path"] = str(manifest_path)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shards", nargs="+")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--window-size", type=int, default=1_000_000)
    parser.add_argument("--min-training-sites", type=int, default=100)
    parser.add_argument("--callable-bed")
    parser.add_argument("--region")
    parser.add_argument("--thresholds", help="comma-separated AC thresholds; also defines rows for an empty scan")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    thresholds = [int(x) for x in args.thresholds.split(",")] if args.thresholds else None
    result = summarize(args.shards, args.out_dir, args.window_size, args.min_training_sites,
                       not args.no_plots, args.callable_bed, args.region, thresholds)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
