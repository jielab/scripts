#!/usr/bin/env python3
"""Cohort IBDmix support peaks from per-person union coverage (BED0 coordinates).

These are ranked candidate regions, not calibrated introgression probabilities.
IBDmix is unphased: each covered person contributes at most one copy, and the
denominator is twice the actual tested cohort, including people with zero calls.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile

import numpy as np

from results_normalize import (
    CHROM_LENGTHS, density_run_path, density_run_signature, reference_role,
    source_class, whole_chromosome_run,
)

VERSION = 1
LINEAGES = ("Neanderthal", "Denisovan")
SQL = """SELECT source_class,chr,sample_id,start,end,source,score
    FROM segments_evidence
    WHERE dataset_id=? AND genome_build=? AND method='ibdmix'
      AND reference_role='scientific_reference'
      AND source_class IN ('Neanderthal','Denisovan')
      AND sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None')
    ORDER BY source_class,chr,sample_id,start,end"""
FIELDS = (
    "rank", "segment_code", "chr", "start", "end", "length_bp", "source_class",
    "genome_build", "n_tested", "peak_start", "peak_end", "peak_fraction",
    "peak_raw_fraction", "peak_effective_copies", "mean_fraction",
    "background_fraction", "extension_threshold", "core_threshold",
    "core_windows", "n_windows", "excess_support_kb", "n_carriers",
    "n_references", "reference_names", "max_lod", "n_calls",
)


def canonical_ref(ref):
    return "Chagyr" if ref == "Chagyrskaya" else ref


def tested_rosters(con, database, dataset, build):
    """Certify full-autosome denominators from saved caller rosters only.

    Missing provenance is not inferred from positive calls or annotation panels.
    Historical controls and locus-only runs cannot certify these windows.
    """
    runs = con.execute("""SELECT chr,raw_file FROM method_runs
        WHERE dataset_id=? AND genome_build=? AND method='ibdmix'
        AND status='complete' AND evidence_eligible=1""", (dataset, build)).fetchall()
    paths = {raw: density_run_path(raw, database) for _, raw in runs}
    targets, warnings = {}, []
    for chrom, raw in runs:
        if chrom not in {str(i) for i in range(1, 23)}:
            continue
        path = paths[raw]
        if not path.is_file():
            warnings.append(f"chr{chrom}: missing run metadata; excluded")
            continue
        text = path.read_text().replace("\\t", "\t")
        if not whole_chromosome_run(text, chrom):
            warnings.append(f"chr{chrom}: locus-only/unverified genome scope; excluded")
            continue
        records = dict(line.split("\t", 1) for line in text.splitlines() if "\t" in line)
        roster = path.parent / "samples" / ("C" + chrom) / "ALL.txt"
        if not roster.is_file():
            roster = path.parent / "samples" / "ALL.txt"
        members = set(roster.read_text().split()) if roster.is_file() else set()
        if not members:
            warnings.append(f"chr{chrom}: missing/empty actual caller roster; excluded")
            continue
        background_only = records.get("background_filter") == "1" and records.get("export_denisovan") != "1"
        for ref in records.get("refs", "").split():
            ref = canonical_ref(ref)
            if reference_role(ref) != "scientific_reference":
                continue
            lineage = source_class(ref)
            if lineage not in LINEAGES or (lineage == "Denisovan" and background_only):
                continue
            targets.setdefault((lineage, chrom), {}).setdefault(ref, set()).update(members)
    return targets, paths, sorted(set(warnings))


def union_intervals(intervals):
    """Merge sorted per-person intervals across all scientific references."""
    left = right = None
    for lo, hi in intervals:
        if hi <= lo:
            continue
        if left is not None and lo <= right:
            right = max(right, hi)
        else:
            if left is not None:
                yield left, right
            left, right = lo, hi
    if left is not None:
        yield left, right


def add_interval(boundary_bp, full_windows, start, end, window_bp):
    """Accumulate exact covered bp using O(1) updates per union interval."""
    first, last = start // window_bp, (end - 1) // window_bp
    if first == last:
        boundary_bp[first] += end - start
    else:
        boundary_bp[first] += (first + 1) * window_bp - start
        boundary_bp[last] += end - last * window_bp
        full_windows[first + 1] += 1
        full_windows[last] -= 1


def window_profiles(rows, targets, lengths, window_bp):
    profiles = []
    # Include completely tested chromosomes with no positive calls in background.
    for (lineage, chrom), refs in sorted(targets.items()):
        if chrom not in lengths:
            continue
        length = lengths[chrom]
        n = len(set().union(*refs.values()))
        starts = np.arange(0, length, window_bp, dtype=np.int64)
        profiles.append(dict(source_class=lineage, chr=chrom, n_tested=n,
            starts=starts, ends=np.minimum(starts + window_bp, length),
            boundary=np.zeros(len(starts)), full=np.zeros(len(starts) + 1)))
    lookup = {(p["source_class"], p["chr"]): p for p in profiles}
    for (lineage, chrom, sample), group in itertools.groupby(rows, key=lambda r: r[:3]):
        p = lookup.get((lineage, chrom))
        if p is None:
            # groupby advances past this group's rows on the next iteration.
            continue
        refs = targets[(lineage, chrom)]
        intervals = ((max(0, r[3]), min(lengths[chrom], r[4])) for r in group
                     if sample in refs.get(canonical_ref(r[5]), set()))
        for lo, hi in union_intervals(intervals):
            add_interval(p["boundary"], p["full"], lo, hi, window_bp)
    for p in profiles:
        widths = p["ends"] - p["starts"]
        p["covered_bp"] = p.pop("boundary") + np.cumsum(p.pop("full")[:-1]) * widths
        p["raw"] = p["covered_bp"] / (2 * p["n_tested"] * widths)
    return profiles


def median_smooth(values, size=3):
    """No padding/reflection across chromosome edges; use available neighbors."""
    radius = size // 2
    return np.array([np.median(values[max(0, i-radius):i+radius+1]) for i in range(len(values))])


def true_runs(mask):
    changes = np.diff(np.r_[False, mask, False].astype(np.int8))
    return zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1))


def call_loci(profiles, dataset, build, window_bp=50_000, floor=0.01, smooth_windows=3):
    backgrounds = {}
    for lineage in LINEAGES:
        group = [p for p in profiles if p["source_class"] == lineage]
        numerator = sum(float(p["covered_bp"].sum()) for p in group)
        denominator = sum(2 * p["n_tested"] * int(p["ends"][-1]) for p in group)
        backgrounds[lineage] = numerator / denominator if denominator else None
    loci = []
    for p in profiles:
        lineage, chrom = p["source_class"], p["chr"]
        baseline = max(floor, backgrounds[lineage])
        core_threshold = 2 * baseline
        smooth = median_smooth(p["raw"], smooth_windows)
        # Raw gate prevents smoothing from manufacturing support in a valley.
        extent = (smooth >= baseline) & (p["raw"] >= baseline)
        core = (smooth >= core_threshold) & (p["raw"] >= core_threshold)
        # A short terminal bin cannot certify a full 50 kb core window.
        core &= (p["ends"] - p["starts"]) == window_bp
        for a, b in true_runs(extent):
            if not any(hi - lo >= 2 for lo, hi in true_runs(core[a:b])):
                continue
            start, end = int(p["starts"][a]), int(p["ends"][b-1])
            if end - start < 2 * window_bp:
                continue
            peak = a + int(np.argmax(smooth[a:b]))
            widths = p["ends"][a:b] - p["starts"][a:b]
            # Area is in fraction*kb (not overlapping-window counts or a p-value).
            excess = float(np.sum(np.maximum(0, smooth[a:b] - baseline) * widths) / 1000)
            identity = f"{dataset}|{build}|{lineage}|{chrom}|{start}|{end}"
            code = "IBD-" + hashlib.sha256(identity.encode()).hexdigest()[:12].upper()
            loci.append(dict(segment_code=code, chr=chrom, start=start, end=end,
                length_bp=end-start, source_class=lineage, genome_build=build,
                n_tested=p["n_tested"], peak_start=int(p["starts"][peak]),
                peak_end=int(p["ends"][peak]), peak_fraction=float(smooth[peak]),
                peak_raw_fraction=float(p["raw"][peak]),
                peak_effective_copies=float(p["raw"][peak] * 2 * p["n_tested"]),
                mean_fraction=float(np.average(p["raw"][a:b], weights=widths)),
                background_fraction=backgrounds[lineage], extension_threshold=baseline,
                core_threshold=core_threshold, core_windows=int(core[a:b].sum()),
                n_windows=int(b-a), excess_support_kb=excess))
    loci.sort(key=lambda r: (-r["excess_support_kb"], -r["peak_fraction"],
                            int(r["chr"]), r["start"], r["source_class"]))
    for rank, row in enumerate(loci, 1):
        row["rank"] = rank
    return loci, backgrounds


def annotate_loci(rows, loci, targets):
    """One streaming pass: no hundred repeated chromosome range scans."""
    index = {}
    for row in loci:
        row.update(_carriers=set(), _references=set(), max_lod=None, n_calls=0)
        index.setdefault((row["source_class"], row["chr"]), []).append(row)
    index = {key: sorted(group, key=lambda r: r["start"]) for key, group in index.items()}
    bounds = {key: ([r["start"] for r in group], [r["end"] for r in group]) for key, group in index.items()}
    for lineage, chrom, sample, start, end, source, score in rows:
        key = (lineage, chrom)
        ref = canonical_ref(source)
        if key not in index or sample not in targets[key].get(ref, set()):
            continue
        starts, ends = bounds[key]
        for i in range(bisect.bisect_right(ends, start), bisect.bisect_left(starts, end)):
            row = index[key][i]
            row["_carriers"].add(sample)
            row["_references"].add(ref)
            row["n_calls"] += 1
            if score is not None and math.isfinite(score):
                row["max_lod"] = score if row["max_lod"] is None else max(score, row["max_lod"])
    for row in loci:
        row["n_carriers"] = len(row.pop("_carriers"))
        refs = sorted(row.pop("_references"))
        row["n_references"] = len(refs)
        row["reference_names"] = ",".join(refs)


def database_stamp(database):
    stat = database.stat()
    return (str(database.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def prepare(database, cache_dir, dataset, build, window_bp=50_000, floor=0.01, smooth_windows=3):
    if window_bp <= 0 or not 0.01 <= floor <= 0.25 or smooth_windows < 3 or smooth_windows % 2 == 0:
        raise ValueError("Require positive window, floor in [0.01,0.25], and odd smoothing width >=3")
    lengths = CHROM_LENGTHS.get(build.replace("GRCh", ""))
    if lengths is None:
        raise ValueError(f"Unsupported genome build: {build}")
    before = database_stamp(database)
    with sqlite3.connect(f"file:{database.resolve()}?mode=ro", uri=True) as con:
        targets, paths, warnings = tested_rosters(con, database, dataset, build)
    provenance = density_run_signature(paths.values())
    settings = dict(window_bp=window_bp, smooth_windows=smooth_windows,
                    floor_fraction=floor, core_multiplier=2, min_core_windows=2)
    signature = dict(version=VERSION, database=before, provenance=provenance,
                     dataset=dataset, build=build, settings=settings,
                     code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:24]
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / key
    if (dest / "loci.json").is_file():
        return dest / "loci.json"
    with tempfile.TemporaryDirectory(prefix="gu-ibdmix-loci-") as work:
        # Indexed SQLite access on /mnt is expensive; stage one immutable copy.
        local = Path(work) / "input.sqlite"
        shutil.copyfile(database, local)
        if database_stamp(database) != before:
            raise RuntimeError("Database changed during snapshot; retry")
        with sqlite3.connect(f"file:{local}?mode=ro", uri=True) as con:
            profiles = window_profiles(con.execute(SQL, (dataset, build)), targets, lengths, window_bp)
            loci, backgrounds = call_loci(profiles, dataset, build, window_bp, floor, smooth_windows)
            if loci:
                annotate_loci(con.execute(SQL, (dataset, build)), loci, targets)
        if database_stamp(database) != before or density_run_signature(paths.values()) != provenance:
            raise RuntimeError("Database or caller roster changed during aggregation; retry")
        meta = dict(signature=signature, settings=settings, dataset=dataset, build=build,
                    backgrounds=backgrounds, candidate_count=len(loci), warnings=warnings,
                    tested=[dict(source_class=p["source_class"], chr=p["chr"], n_tested=p["n_tested"])
                            for p in profiles],
                    definition="Per-person union bp / (window bp * 2 * tested individuals); unphased conservative copy support, not dosage or a calibrated probability")
        stage = Path(tempfile.mkdtemp(prefix=".building-", dir=cache_dir))
        try:
            (stage / "loci.json").write_text(json.dumps(dict(metadata=meta, loci=loci), allow_nan=False))
            with (stage / "candidates.tsv").open("w") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerows(loci)
            with (stage / "windows.tsv").open("w") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("source_class", "chr", "start", "end", "n_tested", "raw_fraction", "smooth_fraction"))
                for p in profiles:
                    smooth = median_smooth(p["raw"], smooth_windows)
                    writer.writerows((p["source_class"], p["chr"], int(lo), int(hi), p["n_tested"], float(raw), float(sm))
                                     for lo, hi, raw, sm in zip(p["starts"], p["ends"], p["raw"], smooth))
            try:
                os.rename(stage, dest)
            except OSError:
                if not (dest / "loci.json").is_file():
                    raise
        finally:
            if stage.exists():
                shutil.rmtree(stage)
    return dest / "loci.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, default=Path("/tmp/gu-shiny-evidence"))
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--build", required=True)
    parser.add_argument("--window-bp", type=int, default=50_000)
    parser.add_argument("--floor", type=float, default=0.01)
    parser.add_argument("--smooth-windows", type=int, default=3)
    args = parser.parse_args()
    print(prepare(args.database, args.cache_dir, args.dataset, args.build,
                  args.window_bp, args.floor, args.smooth_windows))


if __name__ == "__main__":
    main()
