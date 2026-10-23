"""Offline, optional annotation of same-site SNV candidates with local scores.

No prediction, download, API request, P value or selection inference is performed.
The streaming backend requires candidate shards and scores in natural
chromosome/position order; alleles at the same position may occur in any order.
The tabix backend also accepts other chromosome orders when each chromosome
forms one contiguous block with sorted positions. It streams one non-overlapping
coordinate block at a time, so it does not materialize a million-base block.

Candidate coordinates are 1-based. ``score_build`` explicitly asserts the
shared assembly of candidate and score inputs; the CLI verifies this against
the scan manifest. Incompatible input declarations are rejected. TSV score
input requires a header. Explicit ``score_format=
'avi'`` additionally accepts the official six-column, headerless AVI format:
chrom, pos, ref, alt, raw_score, phred. Official AVI input requires GRCh38.
"""

from __future__ import annotations

import csv
import gzip
import importlib
import json
import math
import os
import re
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Iterator


REQUIRED_CANDIDATE_COLUMNS = {
    "chrom", "pos", "ref", "alt", "status", "is_primary", "n_alt",
    "ac", "an", "af", "context",
}
SCORE_COLUMNS = ["local_score", "score_status", "score_reason"]
PAIR_COLUMNS = [
    "chrom", "pos", "ref", "context", "n_alt", "primary_alt",
    "observed_alt", "unobserved_alt", "primary_score", "observed_score",
    "unobserved_score", "delta_observed_minus_unobserved",
    "observed_ac", "observed_an", "observed_af",
]
ALLOWED_STATUSES = {
    "observed_pass", "not_listed_in_release", "present_failed_filter",
    "present_missing_counts", "present_below_threshold",
}
AVI_CAVEAT = (
    "AVI was trained using gnomAD v4.1 allele-frequency proxy labels. Its "
    "association with gnomAD frequency or occupancy is not independent evidence "
    "for selection. AVI, scalar score similarity and AVI feature attributions "
    "are not counts of biological functions or measures of pleiotropy. This "
    "module provides exploratory annotations and arithmetic score differences "
    "only; it performs no hypothesis test or causal inference."
)
_MISSING = {"", ".", "na", "n/a", "null", "none"}
_ALIASES = {
    "chrom": {"chrom", "chromosome", "chr"},
    "pos": {"pos", "position"},
    "ref": {"ref", "reference", "reference_allele"},
    "alt": {"alt", "alternate", "alternative", "alternate_allele"},
}


def canonical_chrom(value: str) -> str:
    """Normalize only standard chromosome aliases; never liftover coordinates."""
    value = str(value).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    if not value:
        raise ValueError("Empty chromosome name")
    if value.isdigit():
        return "chr" + str(int(value))
    if value.upper() in {"X", "Y"}:
        return "chr" + value.upper()
    if value.upper() in {"M", "MT"}:
        return "chrM"
    return "chr" + value


def _chrom_order(chrom: str) -> tuple:
    name = canonical_chrom(chrom)[3:]
    if name.isdigit():
        return (0, int(name))
    if name in {"X", "Y", "M"}:
        return (0, {"X": 1000001, "Y": 1000002, "M": 1000003}[name])
    return (1, tuple((0, int(p)) if p.isdigit() else (1, p)
                     for p in re.split(r"(\d+)", name)))


def _site_key(chrom: str, pos: int) -> tuple:
    return _chrom_order(chrom), pos


def _build(value: str) -> str:
    compact = str(value).strip().lower().replace("_", "").replace("-", "")
    if compact in {"grch38", "hg38", "38", "b38"}:
        return "GRCh38"
    if compact in {"grch37", "hg19", "37", "b37"}:
        return "GRCh37"
    return str(value).strip()


def _check_build_comment(line: str, source: str, expected_build: str) -> None:
    match = re.match(r"^#+\s*(build|assembly|genome_build|reference_build|reference)\s*[:=]\s*(.+)$", line, re.I)
    if not match:
        return
    value = match.group(2).strip().strip('"\'')
    embedded = re.search(r"GRCh(?:37|38)|hg(?:19|38)", value, re.I)
    if embedded:
        value = embedded.group(0)
    elif match.group(1).lower() == "reference":
        return  # A reference FASTA pathname without an assembly is not a declaration.
    if _build(value) != expected_build:
        raise ValueError(f"Build mismatch in {source}: declared {value!r}; candidates require {expected_build}")


@contextmanager
def _text_input(path: Path):
    opener = gzip.open if path.suffix.lower() in {".gz", ".bgz", ".bgzf"} else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        yield handle


def _header_names(line: str) -> list[str]:
    return line.lstrip("#").rstrip("\r\n").split("\t")


def _normalized_names(names: list[str]) -> list[str]:
    return [x.strip().lstrip("#").lower() for x in names]


def _coordinate_indices(names: list[str]) -> dict[str, int] | None:
    normalized = _normalized_names(names)
    indices = {}
    for key, aliases in _ALIASES.items():
        matches = [i for i, name in enumerate(normalized) if name in aliases]
        if len(matches) > 1:
            raise ValueError(f"Ambiguous score header: multiple columns for {key}")
        if not matches:
            return None
        indices[key] = matches[0]
    return indices


def _score_column_index(names: list[str], chosen: str, score_format: str) -> int:
    normalized = _normalized_names(names)
    requested = chosen.strip().lower()
    if score_format == "avi":
        requested = {"raw": "raw_score", "avi_raw": "raw_score",
                     "avi_phred": "phred"}.get(requested, requested)
    if normalized.count(requested) != 1:
        raise ValueError(f"Score column {chosen!r} must appear exactly once; available columns: {names}")
    return normalized.index(requested)


def _variant_fields(chrom: str, pos: str, ref: str, alt: str, source: str) -> tuple[str, int, str, str]:
    try:
        position = int(pos)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid 1-based position {pos!r} in {source}") from exc
    ref, alt = ref.strip().upper(), alt.strip().upper()
    if position < 1 or ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"} or ref == alt:
        raise ValueError(f"Expected a single-base substitution in {source}: {chrom}:{pos}:{ref}>{alt}")
    return canonical_chrom(chrom), position, ref, alt


class _ScoreTable:
    def __init__(self, path: Path, chosen: str, score_format: str, expected_build: str):
        self.path = path
        self.chosen = chosen
        self.score_format = score_format
        self.expected_build = expected_build
        self.skip_through_line = 0
        self.stats = Counter()
        self.names: list[str] = []
        self.indices: dict[str, int] = {}
        self.score_index = -1
        self.build_indices: list[int] = []
        self._inspect()

    def _inspect(self) -> None:
        if self.score_format == "avi":
            self.names = ["chrom", "pos", "ref", "alt", "raw_score", "phred"]
            self.indices = {name: i for i, name in enumerate(self.names[:4])}
            self.score_index = _score_column_index(self.names, self.chosen, "avi")
        with _text_input(self.path) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                if line.startswith("#"):
                    _check_build_comment(line.strip(), str(self.path), self.expected_build)
                names = _header_names(line)
                indices = _coordinate_indices(names)
                if indices is not None:
                    if self.score_format == "tsv":
                        self.names, self.indices = names, indices
                        self.score_index = _score_column_index(names, self.chosen, "tsv")
                        self.build_indices = [i for i, name in enumerate(_normalized_names(names))
                                              if name in {"build", "assembly", "genome_build", "reference_build"}]
                    else:
                        if indices != {"chrom": 0, "pos": 1, "ref": 2, "alt": 3}:
                            raise ValueError("Explicit AVI format requires chrom,pos,ref,alt in the first four columns")
                        # Official positional format is fixed, even when a header exists.
                        if len(names) < 6:
                            raise ValueError("AVI input requires six columns")
                        normalized = _normalized_names(names)
                        if normalized[4] not in {"raw", "raw_score", "rawscore", "avi_raw"} or normalized[5] not in {"phred", "avi_phred"}:
                            raise ValueError("Explicit AVI header must place raw_score and phred in columns five and six")
                    self.skip_through_line = line_number
                    return
                if line.startswith("#"):
                    continue
                if self.score_format == "tsv":
                    raise ValueError(f"Score TSV requires a header with chrom,pos,ref,alt and {self.chosen!r}: {self.path}. Use score_format='avi' only for official six-column AVI input.")
                self.skip_through_line = line_number - 1
                self.parse(line, f"{self.path}:{line_number}", count=False)
                return
        if self.score_format == "tsv":
            raise ValueError(f"No score TSV header found: {self.path}")

    def parse(self, line: str, source: str, count: bool = True):
        values = line.rstrip("\r\n").split("\t")
        highest = max(*self.indices.values(), self.score_index, *self.build_indices)
        if len(values) <= highest:
            raise ValueError(f"Truncated score row in {source}")
        chrom, pos, ref, alt = _variant_fields(*(values[self.indices[x]] for x in ("chrom", "pos", "ref", "alt")), source)
        for index in self.build_indices:
            if _build(values[index]) != self.expected_build:
                raise ValueError(f"Score row build mismatch in {source}: {values[index]!r}")
        raw = values[self.score_index].strip()
        if raw.lower() in _MISSING:
            score = None
        else:
            try:
                score = float(raw)
            except ValueError as exc:
                raise ValueError(f"Nonnumeric score {raw!r} in {source}") from exc
            if not math.isfinite(score):
                raise ValueError(f"Nonfinite score {raw!r} in {source}")
        if count:
            self.stats["score_rows_read"] += 1
        return chrom, pos, ref, alt, score

    def records(self):
        with _text_input(self.path) as handle:
            for line_number, line in enumerate(handle, 1):
                if line_number <= self.skip_through_line or not line.strip():
                    continue
                if line.startswith("#"):
                    _check_build_comment(line.strip(), str(self.path), self.expected_build)
                    continue
                yield self.parse(line, f"{self.path}:{line_number}")


def _score_sites(records: Iterable, stats: Counter) -> Iterator:
    previous = None
    current = None
    scores = {}
    for chrom, pos, ref, alt, score in records:
        key = _site_key(chrom, pos)
        if previous is not None and key < previous:
            raise ValueError("Score input is not in natural chromosome/position order (chr1,chr2,...,chr10,...,chrX). Supply sorted local scores or an indexed bgzip table with backend='tabix'.")
        if current is not None and key != current:
            yield current, scores
            scores = {}
        allele = (ref, alt)
        if allele in scores:
            if scores[allele] != score:
                raise ValueError(f"Conflicting duplicate score: {chrom}:{pos}:{ref}>{alt}")
            stats["identical_duplicate_score_rows"] += 1
        else:
            scores[allele] = score
        previous = current = key
    if current is not None:
        yield current, scores


class _Candidates:
    def __init__(self, paths: list[Path], expected_build: str, natural_order: bool = True):
        self.paths = paths
        self.expected_build = expected_build
        self.natural_order = natural_order
        self.header: list[str] | None = None

    def records(self):
        previous = None
        previous_chrom = None
        finished_chroms = set()
        for path in self.paths:
            with _text_input(path) as handle:
                header = None
                for line_number, line in enumerate(handle, 1):
                    if not line.strip():
                        continue
                    if line.startswith("#"):
                        _check_build_comment(line.strip(), str(path), self.expected_build)
                    if header is None:
                        names = [x.strip() for x in _header_names(line)]
                        if REQUIRED_CANDIDATE_COLUMNS.issubset(names):
                            if len(set(names)) != len(names):
                                raise ValueError(f"Duplicate candidate header fields in {path}")
                            if set(names) & set(SCORE_COLUMNS):
                                raise ValueError(f"Candidate input already contains reserved score columns: {path}")
                            header = names
                            if self.header is None:
                                self.header = names
                            elif self.header != names:
                                raise ValueError("Candidate shard headers differ; preserve the same column order in every shard")
                            continue
                        if line.startswith("#"):
                            continue
                        missing = sorted(REQUIRED_CANDIDATE_COLUMNS - set(names))
                        raise ValueError(f"Missing candidate columns in {path}: {missing}")
                    if line.startswith("#"):
                        continue
                    values = line.rstrip("\r\n").split("\t")
                    if len(values) != len(header):
                        raise ValueError(f"Candidate row has wrong column count: {path}:{line_number}")
                    row = dict(zip(header, values))
                    chrom, pos, ref, alt = _variant_fields(row["chrom"], row["pos"], row["ref"], row["alt"], f"{path}:{line_number}")
                    row.update(chrom=chrom, pos=str(pos), ref=ref, alt=alt)
                    if row["status"] not in ALLOWED_STATUSES:
                        raise ValueError(f"Unknown candidate status {row['status']!r} in {path}")
                    if row["is_primary"] not in {"0", "1"} or row.get("primary_tie", "0") not in {"0", "1"}:
                        raise ValueError(f"is_primary and primary_tie must be 0 or 1 in {path}")
                    if row["n_alt"] not in {"0", "1", "2", "3"}:
                        raise ValueError(f"n_alt must be 0, 1, 2, or 3 in {path}")
                    for field in ("build", "assembly", "genome_build", "reference_build"):
                        if field in row and _build(row[field]) != self.expected_build:
                            raise ValueError(f"Candidate build mismatch: {row[field]!r} in {path}")
                    key = _site_key(chrom, pos)
                    if self.natural_order and previous is not None and key < previous:
                        raise ValueError("Candidate shards must be supplied in natural chromosome/position order, without overlapping or reversed shards")
                    if not self.natural_order:
                        if chrom != previous_chrom:
                            if chrom in finished_chroms:
                                raise ValueError("Candidate input revisits a chromosome; each chromosome must form one contiguous block")
                            if previous_chrom is not None:
                                finished_chroms.add(previous_chrom)
                        elif previous is not None and key < previous:
                            raise ValueError("Candidate positions must be sorted within each chromosome")
                    previous = key
                    previous_chrom = chrom
                    yield key, row
                if header is None:
                    raise ValueError(f"No candidate header found in {path}")

    def sites(self):
        current = None
        rows = []
        for key, row in self.records():
            if current is not None and key != current:
                yield current, rows
                rows = []
            current = key
            if any(x["ref"] == row["ref"] and x["alt"] == row["alt"] for x in rows):
                raise ValueError(f"Duplicate candidate: {row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}")
            if rows and rows[0]["ref"] != row["ref"]:
                raise ValueError(f"Conflicting candidate reference bases: {row['chrom']}:{row['pos']}")
            rows.append(row)
        if current is not None:
            yield current, rows


def _next(iterator):
    return next(iterator, None)


def _matched_sites(candidate_sites, score_sites):
    current_score = _next(score_sites)
    for key, candidates in candidate_sites:
        while current_score is not None and current_score[0] < key:
            current_score = _next(score_sites)
        yield candidates, current_score[1] if current_score is not None and current_score[0] == key else {}


def _tabix_matched_sites(candidate_sites, table: _ScoreTable, tabix, block_size: int, stats: Counter):
    contigs = {}
    for name in tabix.contigs:
        canonical = canonical_chrom(name)
        if canonical in contigs:
            raise ValueError(f"Ambiguous tabix contig aliases for {canonical}: {contigs[canonical]}, {name}")
        contigs[canonical] = name
    block = None
    score_sites = iter(())
    current_score = None
    for key, candidates in candidate_sites:
        chrom, pos = candidates[0]["chrom"], int(candidates[0]["pos"])
        start = ((pos - 1) // block_size) * block_size
        next_block = (chrom, start)
        if next_block != block:
            # Exhaust the prior queried block to validate all its duplicates/order.
            for _ in score_sites:
                pass
            block = next_block
            current_score = None
            if chrom not in contigs:
                score_sites = iter(())
                stats["blocks_without_score_contig"] += 1
            else:
                stats["tabix_regions_fetched"] += 1
                lines = tabix.fetch(contigs[chrom], start, start + block_size)
                records = (table.parse(line, f"{table.path}:{chrom}:{start + 1}-{start + block_size}")
                           for line in lines if line.strip() and not line.startswith("#"))
                score_sites = iter(_score_sites(records, table.stats))
                current_score = _next(score_sites)
        while current_score is not None and current_score[0] < key:
            current_score = _next(score_sites)
        yield candidates, current_score[1] if current_score is not None and current_score[0] == key else {}
    for _ in score_sites:
        pass


def _number(value: float | None) -> str:
    return "." if value is None else format(value, ".17g")


def _pair(rows: list[dict[str, str]], values: dict[tuple[str, str], float | None], stats: Counter):
    if len(rows) != 3 or any(row["n_alt"] != "2" for row in rows):
        return None
    if any(row.get("primary_tie", "0") == "1" for row in rows):
        stats["pair_sites_excluded_primary_tie"] += 1
        return None
    primary = [r for r in rows if r["is_primary"] == "1" and r["status"] == "observed_pass"]
    observed = [r for r in rows if r["is_primary"] == "0" and r["status"] == "observed_pass"]
    absent = [r for r in rows if r["status"] == "not_listed_in_release" and r["is_primary"] == "0"]
    if len(primary) != 1 or len(observed) != 1 or len(absent) != 1:
        stats["pair_sites_excluded_status"] += 1
        return None
    p, o, u = primary[0], observed[0], absent[0]
    primary_score = values.get((p["ref"], p["alt"]))
    observed_score = values.get((o["ref"], o["alt"]))
    unobserved_score = values.get((u["ref"], u["alt"]))
    if observed_score is None or unobserved_score is None:
        stats["pair_sites_excluded_missing_score"] += 1
        return None
    difference = observed_score - unobserved_score
    if not math.isfinite(difference):
        raise ValueError(f"Score difference overflow at {o['chrom']}:{o['pos']}")
    return dict(
        chrom=o["chrom"], pos=o["pos"], ref=o["ref"], context=o["context"], n_alt="2",
        primary_alt=p["alt"], observed_alt=o["alt"], unobserved_alt=u["alt"],
        primary_score=_number(primary_score), observed_score=_number(observed_score),
        unobserved_score=_number(unobserved_score), delta_observed_minus_unobserved=_number(difference),
        observed_ac=o["ac"], observed_an=o["an"], observed_af=o["af"],
    )


def preflight_scores(
    scores_path,
    score_column: str,
    score_build: str = "GRCh38",
    backend: str = "auto",
    score_format: str = "tsv",
) -> dict[str, Any]:
    """Check local score prerequisites before expensive candidate construction.

    Only the header (or first AVI row) and optional tabix index are opened. No
    score-table scan, candidate read, region fetch or output write is performed.
    The returned runtime objects are reused by ``annotate_scores``; the CLI may
    discard them after this early check. Row values/order remain validated by
    the actual streaming annotation.
    """
    expected_build = _build(score_build)
    if expected_build not in {"GRCh38", "GRCh37"}:
        raise ValueError("score_build must explicitly identify GRCh38 or GRCh37; no liftover is performed")
    if score_format not in {"tsv", "avi"}:
        raise ValueError("score_format must be 'tsv' or 'avi'")
    if score_format == "avi" and expected_build != "GRCh38":
        raise ValueError(f"score_build={score_build!r} does not match official AVI's GRCh38 assembly; no liftover is performed")
    if backend == "streaming":
        backend = "stream"
    if backend not in {"auto", "stream", "tabix"}:
        raise ValueError("backend must be 'auto', 'stream', or 'tabix'")
    score_path = Path(scores_path).resolve()
    if not score_path.is_file():
        raise FileNotFoundError(f"Local input file not found: {score_path}")
    table = _ScoreTable(score_path, score_column, score_format, expected_build)
    has_index = any(Path(str(score_path) + ext).is_file() for ext in (".tbi", ".csi"))
    pysam = None
    if backend != "stream" and has_index:
        try:
            pysam = importlib.import_module("pysam")
        except ImportError:
            pass
    if backend == "tabix" and (pysam is None or not has_index):
        raise ValueError("backend='tabix' requires pysam and an adjacent .tbi or .csi index; use backend='stream' for sorted local input")
    chosen_backend = "tabix" if backend != "stream" and pysam is not None and has_index else "stream"
    if chosen_backend == "tabix":
        # Presence alone does not establish a usable BGZF/index pair. Opening
        # and closing checks it without fetching genomic intervals.
        index = pysam.TabixFile(str(score_path))
        index.close()
    return dict(expected_build=expected_build, score_path=score_path, table=table,
                backend_requested=backend, backend=chosen_backend, pysam=pysam)


def annotate_scores(
    candidate_files: list,
    scores_path,
    out_dir,
    score_column: str,
    score_build: str = "GRCh38",
    score_label: str = "local",
    score_version: str = "unspecified",
    backend: str = "auto",
    block_size: int = 1_000_000,
    score_format: str = "tsv",
) -> dict[str, Any]:
    """Annotate local candidates, returning counts, provenance and output paths.

    ``backend`` is auto, stream (streaming is an alias), or tabix. Auto uses
    pysam/tabix only when pysam and an adjacent .tbi/.csi index are available.
    Otherwise it performs a single bounded-memory sorted merge scan. Existing
    final outputs are replaced only after annotation succeeds. ``pairs.tsv.gz``
    is always written with a header and contains only eligible numeric contrasts;
    an empty file is a valid result. Failure removes this call's partial files.
    """
    if not isinstance(block_size, int) or isinstance(block_size, bool) or block_size < 1:
        raise ValueError("block_size must be a positive integer")
    paths = [Path(p).resolve() for p in candidate_files]
    if not paths:
        raise ValueError("candidate_files must contain at least one local TSV shard")
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Local input file not found: {path}")
    prepared = preflight_scores(scores_path, score_column, score_build, backend, score_format)
    expected_build, score_path, table = (prepared[x] for x in ("expected_build", "score_path", "table"))
    backend, chosen_backend, pysam = (prepared[x] for x in ("backend_requested", "backend", "pysam"))
    destination = Path(out_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    outputs = {
        "candidate_scores": destination / "candidate_scores.tsv.gz",
        "pairs": destination / "pairs.tsv.gz",
        "metadata": destination / "score_annotation.json",
    }
    if set(outputs.values()) & set([*paths, score_path]):
        raise ValueError("Output paths must not overwrite input candidate or score files")
    partials = {key: Path(str(path) + ".partial") for key, path in outputs.items()}
    if any(path.exists() for path in partials.values()):
        raise ValueError("A .partial score output already exists; use another output directory or remove the incomplete files before retrying")
    stats = Counter(total_candidates=0, available=0, absent=0, candidate_sites=0, pairs=0)
    tabix = None
    try:
        candidates = _Candidates(paths, expected_build, natural_order=chosen_backend != "tabix")
        candidate_sites = iter(candidates.sites())
        first = _next(candidate_sites)  # Establish and validate the common header.
        from itertools import chain
        candidate_sites = chain(() if first is None else (first,), candidate_sites)
        if chosen_backend == "tabix":
            tabix = pysam.TabixFile(str(score_path))
            matches = _tabix_matched_sites(candidate_sites, table, tabix, block_size, stats)
        else:
            matches = _matched_sites(candidate_sites, iter(_score_sites(table.records(), table.stats)))
        with gzip.open(partials["candidate_scores"], "wt", encoding="utf-8", newline="", compresslevel=6) as candidate_out, gzip.open(partials["pairs"], "wt", encoding="utf-8", newline="", compresslevel=6) as pair_out:
            annotation_writer = csv.DictWriter(candidate_out, fieldnames=candidates.header + SCORE_COLUMNS, delimiter="\t", lineterminator="\n")
            pair_writer = csv.DictWriter(pair_out, fieldnames=PAIR_COLUMNS, delimiter="\t", lineterminator="\n")
            annotation_writer.writeheader()
            pair_writer.writeheader()
            for rows, values in matches:
                stats["candidate_sites"] += 1
                for row in rows:
                    variant = (row["ref"], row["alt"])
                    score = values.get(variant)
                    reason = "matched" if score is not None else ("missing_score_value" if variant in values else "no_matching_record")
                    status = "available" if score is not None else "absent"
                    annotation_writer.writerow({**row, "local_score": _number(score), "score_status": status, "score_reason": reason})
                    stats["total_candidates"] += 1
                    stats[status] += 1
                    stats[f"reason_{reason}"] += 1
                pair = _pair(rows, values, stats)
                if pair is not None:
                    pair_writer.writerow(pair)
                    stats["pairs"] += 1
        metadata = {
            "schema_version": 1,
            "score_build": expected_build, "score_label": str(score_label),
            "score_version": str(score_version), "score_column": score_column,
            "score_format": score_format, "scores_path": str(score_path),
            "candidate_files": [str(p) for p in paths],
            "backend_requested": backend, "backend": chosen_backend,
            "candidate_order": "contiguous_chromosomes_sorted_positions" if chosen_backend == "tabix" else "natural_chromosome_position",
            "block_size": block_size if chosen_backend == "tabix" else None,
            "score_validation_scope": "queried_coordinate_blocks" if chosen_backend == "tabix" else "streamed_prefix_through_last_candidate",
            "counts": {**dict(stats), **dict(table.stats)},
            "outputs": {key: str(path) for key, path in outputs.items()},
            "warnings": [AVI_CAVEAT, "not_listed_in_release means absent from this released callset, not proven never mutated or extinct. Score distances do not identify ancestry or the allele from which a mutation arose."],
            "interpretation": "Exploratory local annotation and within-site observed-minus-not-listed arithmetic contrast; no P values or causal claims.",
            "sources": ["https://doi.org/10.64898/2026.09.16.26363192", "https://github.com/Ensembl/VEP_plugins/blob/release/116/AVI.pm"],
        }
        with open(partials["metadata"], "wt", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
        if tabix is not None:
            tabix.close()
            tabix = None
        for key in ("candidate_scores", "pairs", "metadata"):
            os.replace(partials[key], outputs[key])
        return metadata
    finally:
        if tabix is not None:
            tabix.close()
        for path in partials.values():
            if path.exists():
                path.unlink()
