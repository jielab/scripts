#!/usr/bin/env python3
"""Current PhyML report helpers for TSV publication and interval arithmetic."""
from __future__ import annotations

import csv
import gzip
import os
import re
import tempfile
from pathlib import Path
from typing import Iterable


def chrom(value):
    """Normalize the chromosome notation used by current input files."""
    text = re.sub(r"^chr", "", str(value), flags=re.I)
    return "X" if text == "23" else text


def truth(value):
    return str(value).strip().upper() in {"1", "TRUE", "T", "YES", "PASS"}


def rows(path):
    """Read a current TSV or TSV.GZ file without converting sample identities."""
    if path is None:
        return []
    path = Path(path)
    if not path.is_file():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        result = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in result):
        raise ValueError(f"Malformed TSV fields: {path}")
    return result


def atomic_text(path: Path, text: str):
	path.parent.mkdir(parents=True, exist_ok=True)
	fd, name = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
	try:
		with os.fdopen(fd, "w", encoding="utf-8", newline="") as f:
			f.write(text)
		os.replace(name, path)
	finally:
		if os.path.exists(name):
			os.unlink(name)


def write_tsv(path: Path, rows: list[dict], fields=None):
	import io

	fields = fields or list(dict.fromkeys(k for r in rows for k in r))
	if not fields:
		fields = ["status"]
	s = io.StringIO(newline="")
	w = csv.DictWriter(s, fieldnames=fields, delimiter="\t", extrasaction="raise", lineterminator="\n")
	w.writeheader()
	for row in rows:
		w.writerow({k: ("" if v is None else v) for k, v in row.items()})
	atomic_text(path, s.getvalue())


def union(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
	out = []
	for st, en in sorted(set(intervals)):
		if st < 0 or en <= st:
			raise ValueError(f"Invalid half-open interval: {st}, {en}")
		if out and st <= out[-1][1]:
			out[-1] = (out[-1][0], max(en, out[-1][1]))
		else:
			out.append((st, en))
	return out


def overlap(a, b):
	return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def coverage(candidate, intervals):
	if candidate[0] < 0 or candidate[1] <= candidate[0]:
		raise ValueError(f"Invalid candidate interval {candidate}")
	ints = list(intervals)
	clipped = [(max(candidate[0], s), min(candidate[1], e)) for s, e in ints if overlap(candidate, (s, e))]
	reduced = union(clipped)
	bp = sum(e - s for s, e in reduced)
	size = candidate[1] - candidate[0]
	return dict(
		overlap_bp=bp,
		union_fraction=bp / size,
		best_single_fraction=max((overlap(candidate, x) / size for x in ints), default=0),
		n_covered_blocks=len(reduced),
		longest_block_fraction=max(((e - s) / size for s, e in reduced), default=0),
	)
