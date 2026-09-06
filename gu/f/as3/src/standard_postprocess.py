"""Canonical ArchaicSeeker3 BED post-processing.

The standard output is built from unmerged raw segments in this order:

1. retain raw segments of at least 5 kb with score >= 0;
2. exact-merge retained segments using the requested merge distance;
3. write the combined result and ancestry-specific BED files.

Raw BED and SNP-detail files are never modified.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src.merge_bed_segments import merge_segments_exact_streaming


logger = logging.getLogger(__name__)

BED_COLUMNS = [
    "chr",
    "start_pos",
    "end_pos",
    "haplotype",
    "ancestry_label",
    "num_snps",
    "avg_prob",
    "archaic_snps",
    "african_snps",
    "sample_hap_id",
]

STANDARD_MIN_LENGTH_BP = 5_000
STANDARD_MIN_SCORE = 0.0

RAW_BED_NAME = "introgression.raw.bed"
RAW_SNPS_NAME = "introgression.raw.snps.gz"
COMBINED_BED_NAME = "introgression.bed"
ANCESTRY_BED_NAMES = {
    1: "introgression.denisovan.bed",
    2: "introgression.neanderthal.bed",
    3: "introgression.mosaic.bed",
}

RAW_BED_CHUNK_ROWS = 100_000


def read_raw_bed(path: Path) -> pd.DataFrame:
    """Read an AS3 10-column raw BED, including a valid empty file."""

    try:
        return pd.read_csv(path, sep="\t", header=None, names=BED_COLUMNS)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=BED_COLUMNS)


def filter_raw_segments(
    segments: pd.DataFrame,
    *,
    min_length_bp: int = STANDARD_MIN_LENGTH_BP,
    min_score: float = STANDARD_MIN_SCORE,
) -> pd.DataFrame:
    """Return raw segments passing the canonical pre-merge thresholds."""

    if min_length_bp < 0:
        raise ValueError("min_length_bp must be non-negative")

    if segments.empty:
        return segments.copy()

    lengths = segments["end_pos"] - segments["start_pos"]
    keep = (lengths >= min_length_bp) & (segments["avg_prob"] >= min_score)
    filtered = segments.loc[keep].copy()
    logger.info(
        "Standard raw filter (length >= %d bp, score >= %s): %d -> %d segments",
        min_length_bp,
        min_score,
        len(segments),
        len(filtered),
    )
    return filtered


def read_and_filter_raw_bed(
    path: Path,
    *,
    min_length_bp: int = STANDARD_MIN_LENGTH_BP,
    min_score: float = STANDARD_MIN_SCORE,
    chunk_rows: int = RAW_BED_CHUNK_ROWS,
) -> tuple[int, pd.DataFrame]:
    """Read and filter raw BED rows in bounded pandas chunks."""

    if min_length_bp < 0:
        raise ValueError("min_length_bp must be non-negative")
    if chunk_rows <= 0:
        raise ValueError("chunk_rows must be positive")

    raw_count = 0
    filtered_parts: list[pd.DataFrame] = []
    try:
        chunks = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=BED_COLUMNS,
            chunksize=chunk_rows,
        )
        for chunk in chunks:
            raw_count += len(chunk)
            lengths = chunk["end_pos"] - chunk["start_pos"]
            keep = (lengths >= min_length_bp) & (chunk["avg_prob"] >= min_score)
            if keep.any():
                filtered_parts.append(chunk.loc[keep].copy())
    except pd.errors.EmptyDataError:
        pass

    if filtered_parts:
        filtered = pd.concat(filtered_parts, ignore_index=True)
    else:
        filtered = pd.DataFrame(columns=BED_COLUMNS)

    logger.info(
        "Standard raw filter (length >= %d bp, score >= %s): %d -> %d "
        "segments (chunk_rows=%d)",
        min_length_bp,
        min_score,
        raw_count,
        len(filtered),
        chunk_rows,
    )
    return raw_count, filtered


def _write_bed(segments: pd.DataFrame, path: Path) -> None:
    ordered = segments.reindex(columns=BED_COLUMNS)
    ordered.to_csv(path, sep="\t", index=False, header=False)
    logger.info("Written %d segments to %s", len(ordered), path)


def run_standard_postprocess(
    output_dir: Path,
    *,
    merge_distance: int,
    min_length_bp: int = STANDARD_MIN_LENGTH_BP,
    min_score: float = STANDARD_MIN_SCORE,
) -> dict[str, int | str]:
    """Create canonical combined and ancestry-specific BED outputs.

    The returned counts are suitable for logs and regression tests. The raw
    files are read-only inputs and remain available for custom reprocessing.
    """

    if merge_distance < 0:
        raise ValueError("merge_distance must be non-negative")

    output_dir = Path(output_dir)
    raw_bed_path = output_dir / RAW_BED_NAME
    raw_snps_path = output_dir / RAW_SNPS_NAME
    combined_bed_path = output_dir / COMBINED_BED_NAME

    if not raw_bed_path.is_file():
        raise FileNotFoundError(f"Raw AS3 BED not found: {raw_bed_path}")

    raw_count, filtered_segments = read_and_filter_raw_bed(
        raw_bed_path,
        min_length_bp=min_length_bp,
        min_score=min_score,
    )

    if filtered_segments.empty or merge_distance == 0:
        merged_segments = filtered_segments.copy()
        if not merged_segments.empty:
            merged_segments = merged_segments.sort_values(
                ["sample_hap_id", "chr", "start_pos"]
            )
    else:
        if not raw_snps_path.is_file():
            raise FileNotFoundError(
                f"Raw AS3 SNP details required for exact merge: {raw_snps_path}"
            )
        merged_segments = merge_segments_exact_streaming(
            filtered_segments,
            str(raw_snps_path),
            merge_distance,
            min_snps_per_segment=2,
            mosaic_minority_threshold=0.20,
        )

    if merged_segments.empty:
        merged_segments = pd.DataFrame(columns=BED_COLUMNS)

    _write_bed(merged_segments, combined_bed_path)

    ancestry_counts: dict[int, int] = {}
    for label, filename in ANCESTRY_BED_NAMES.items():
        ancestry_segments = merged_segments.loc[
            merged_segments["ancestry_label"] == label
        ].copy()
        _write_bed(ancestry_segments, output_dir / filename)
        ancestry_counts[label] = len(ancestry_segments)

    return {
        "raw_segments": raw_count,
        "filtered_segments": len(filtered_segments),
        "merged_segments": len(merged_segments),
        "denisovan_segments": ancestry_counts[1],
        "neanderthal_segments": ancestry_counts[2],
        "mosaic_segments": ancestry_counts[3],
        "merge_distance": merge_distance,
        "combined_bed": str(combined_bed_path),
    }
