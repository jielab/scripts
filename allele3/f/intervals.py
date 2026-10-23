"""Callable BED handling and streaming conversion of gnomAD all-sites AN."""
from __future__ import annotations

from array import array
from bisect import bisect_right
import csv
import math
import os
from pathlib import Path

from core import dump_json, fingerprint, open_text


def canonical(chrom):
    value = str(chrom).removeprefix("chr")
    if value.isdigit():
        return str(int(value))
    return "MT" if value in ("M", "MT") else value


def order(chrom):
    c = canonical(chrom)
    return (0, int(c)) if c.isdigit() else (1, {"X": 23, "Y": 24, "MT": 25}.get(c, 100), c)


class IntervalMask:
    """Merged 0-based half-open intervals, compact arrays for point lookup."""

    def __init__(self, path, allowed_chroms=None):
        allowed = None if allowed_chroms is None else {canonical(c) for c in allowed_chroms}
        self.starts, self.ends = {}, {}
        intervals = {}
        with open_text(path) as f:
            for line_no, line in enumerate(f, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                cells = line.rstrip("\n").split("\t")
                if len(cells) < 3:
                    raise ValueError(f"Invalid BED at {path}:{line_no}")
                chrom = canonical(cells[0])
                try:
                    start, end = int(cells[1]), int(cells[2])
                except ValueError as e:
                    raise ValueError(f"BED must use integer 0-based half-open coordinates: {path}:{line_no}") from e
                if start < 0 or end <= start:
                    raise ValueError(f"Invalid BED interval at {path}:{line_no}")
                if allowed is None or chrom in allowed:
                    intervals.setdefault(chrom, []).append((start, end))
        for chrom, rows in intervals.items():
            rows.sort()
            merged = []
            for start, end in rows:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.starts[chrom] = array("Q", (s for s, _ in merged))
            self.ends[chrom] = array("Q", (e for _, e in merged))

    def contains(self, chrom, pos1):
        c, point = canonical(chrom), int(pos1) - 1
        starts = self.starts.get(c)
        if starts is None:
            return False
        i = bisect_right(starts, point) - 1
        return i >= 0 and point < self.ends[c][i]


def make_callable(all_sites_an, output_bed, min_an=100000, chromosomes=None):
    """Stream a scalar AN table into merged BED; never builds a per-base map.

    Supported headers: locus,AN (Hail export: chr1:123) or chrom,pos,AN.
    Input coordinates are always 1-based; output BED is 0-based half-open.
    AN is a coverage criterion, not a reconstructed joint allele denominator.
    """
    if min_an < 1:
        raise ValueError("min_an must be positive")
    source, output = Path(all_sites_an), Path(output_bed)
    if source.resolve() == output.resolve():
        raise ValueError("Input and output must be different files")
    if output.exists():
        raise FileExistsError(f"Output exists: {output}; choose a new path")
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(output.name + ".partial")
    allowed = None if not chromosomes else {canonical(c) for c in chromosomes}
    n_rows = n_kept = n_intervals = 0
    previous = None
    current = None
    try:
        with open_text(source) as inp, tmp.open("w", encoding="utf-8") as out:
            lines = (line for line in inp if line.strip() and not line.startswith("##"))
            reader = csv.DictReader(lines, delimiter="\t")
            names = {n.lstrip("#").lower(): n for n in (reader.fieldnames or [])}
            an_key = names.get("an")
            locus_key = names.get("locus")
            chrom_key = next((names[x] for x in ("chrom", "chr", "contig", "chromosome") if x in names), None)
            pos_key = next((names[x] for x in ("pos", "position") if x in names), None)
            if not an_key or not (locus_key or (chrom_key and pos_key)):
                raise ValueError("all-sites AN table requires header locus+AN or chrom+pos+AN; AN must be scalar")
            for row in reader:
                n_rows += 1
                if locus_key:
                    try:
                        raw_chrom, raw_pos = row[locus_key].rsplit(":", 1)
                    except ValueError as e:
                        raise ValueError(f"Invalid Hail locus on data row {n_rows}") from e
                else:
                    raw_chrom, raw_pos = row[chrom_key], row[pos_key]
                chrom, pos = canonical(raw_chrom), int(raw_pos)
                if pos < 1:
                    raise ValueError("all-sites positions must be 1-based positive integers")
                key = (order(chrom), pos)
                if previous is not None and key <= previous:
                    raise ValueError("all-sites AN must be naturally chromosome/position sorted, without duplicates")
                previous = key
                try:
                    an_value = float(row[an_key])
                except ValueError as e:
                    raise ValueError("AN must be a scalar number, not a Hail array; use the official TSV export") from e
                if not math.isfinite(an_value) or an_value < 0 or not an_value.is_integer():
                    raise ValueError(f"Invalid AN on row {n_rows}")
                if allowed is not None and chrom not in allowed:
                    continue
                if an_value < min_an:
                    continue
                n_kept += 1
                if current is not None and chrom == current[0] and pos - 1 == current[2]:
                    current[2] = pos
                else:
                    if current is not None:
                        out.write("\t".join(map(str, current)) + "\n")
                        n_intervals += 1
                    current = [chrom, pos - 1, pos]
            if current is not None:
                out.write("\t".join(map(str, current)) + "\n")
                n_intervals += 1
        os.replace(tmp, output)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    result = {"input": fingerprint(source), "output": str(output.resolve()),
              "min_an": min_an, "input_rows": n_rows, "callable_bases": n_kept,
              "merged_intervals": n_intervals, "coordinates": "BED 0-based half-open",
              "warning": "Overall pre-split AN defines the mask only; it does not recover joint REF/ALT counts."}
    dump_json(str(output) + ".json", result)
    return result
