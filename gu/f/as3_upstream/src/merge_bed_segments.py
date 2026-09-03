#!/usr/bin/env python3
"""
BED Segment Merger (No Genetic Map Version)

This standalone script takes a raw BED file (from ArchaicSeeker3 with
merge_distance=0) and merges nearby segments.

This version is designed for BED files WITHOUT genetic map (NO cM columns).

Usage:
    python merge_bed_segments.py --input raw.bed --output merged.bed --merge-distance 10000

Features:
    - Exact mode merges by sample/haplotype regardless of the pre-merge label
    - Uses 10-column BED format (NO cM columns)
    - Recalculates statistics (num_snps, avg_prob) after merging
    - Supports exact merge mode using SNP details file
"""

import argparse
import csv
import pandas as pd
import numpy as np
import sys
import logging
import gzip
import os

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# EXACT MERGE: SNP-level detail functions for exact merge replication
# ============================================================================

def load_snp_details(snp_details_path):
    """
    Load SNP details file (gzipped TSV)

    Args:
        snp_details_path: Path to .raw.snps.gz file

    Returns:
        dict: {(chr, sample_hap_id, start_pos, end_pos): snp_data_dict}
    """
    if not os.path.exists(snp_details_path):
        raise FileNotFoundError(f"SNP details file not found: {snp_details_path}")

    logger.info(f"Loading SNP details from {snp_details_path}")

    with gzip.open(snp_details_path, 'rt', encoding='utf-8') as f:
        df = pd.read_csv(f, sep='\t')

    # Build lookup dictionary using composite key
    snp_details_dict = {}
    for _, row in df.iterrows():
        key = (row['chr'], row['sample_hap_id'], row['start_pos'], row['end_pos'])
        snp_details_dict[key] = {
            'positions': np.array([int(x) for x in row['snp_positions'].split(',')]),
            'states': np.array([int(x) for x in row['snp_states'].split(',')]),
            'prob_label1': np.array([float(x) for x in row['snp_prob_label1'].split(',')]),
            'prob_label2': np.array([float(x) for x in row['snp_prob_label2'].split(',')]),
        }

    logger.info(f"Loaded SNP details for {len(snp_details_dict)} segments")
    return snp_details_dict


def _segment_key(chr_value, sample_hap_id, start_pos, end_pos):
    """Return the normalized composite key shared by BED and SNP-detail rows."""
    return (str(chr_value), str(sample_hap_id), int(start_pos), int(end_pos))


class _ExactMergeAccumulator:
    """Memory-light state for one final exact-merge group."""

    __slots__ = (
        "chr",
        "haplotype",
        "sample_hap_id",
        "min_pos",
        "max_pos",
        "n1",
        "n2",
        "prob1_sum",
        "prob2_sum",
        "detail_rows",
    )

    def __init__(self, chr_value, haplotype, sample_hap_id):
        self.chr = chr_value
        self.haplotype = haplotype
        self.sample_hap_id = sample_hap_id
        self.min_pos = None
        self.max_pos = None
        self.n1 = 0
        self.n2 = 0
        self.prob1_sum = 0.0
        self.prob2_sum = 0.0
        self.detail_rows = 0


def _build_exact_merge_plan(df, merge_distance):
    """Plan merge groups without loading any SNP-detail arrays."""
    key_to_group = {}
    accumulators = []

    for _, sample_segments in df.groupby("sample_hap_id"):
        ordered = sample_segments.sort_values("start_pos")
        group_index = None
        previous_end = None

        for row in ordered.itertuples(index=False):
            start_pos = int(row.start_pos)
            end_pos = int(row.end_pos)
            if group_index is None or start_pos - previous_end > merge_distance:
                group_index = len(accumulators)
                accumulators.append(
                    _ExactMergeAccumulator(row.chr, row.haplotype, row.sample_hap_id)
                )

            key = _segment_key(row.chr, row.sample_hap_id, start_pos, end_pos)
            if key in key_to_group:
                raise ValueError(
                    "Duplicate filtered AS3 segment key cannot be exact-merged: "
                    f"{key}"
                )
            key_to_group[key] = group_index
            previous_end = end_pos

    logger.info(
        "Prepared streaming exact-merge plan: %d retained segments -> %d groups",
        len(key_to_group),
        len(accumulators),
    )
    return key_to_group, accumulators


def _parse_detail_array(row, field, dtype, key):
    value = row.get(field)
    if value is None:
        raise ValueError(f"SNP-detail column {field!r} is missing for segment {key}")
    if value == "":
        return np.array([], dtype=dtype)

    values = np.fromstring(value, sep=",", dtype=dtype)
    expected = value.count(",") + 1
    if len(values) != expected:
        raise ValueError(
            f"Malformed {field} for segment {key}: expected {expected} values, "
            f"parsed {len(values)}"
        )
    return values


def _accumulate_detail_row(accumulator, row, key):
    positions = _parse_detail_array(row, "snp_positions", np.int64, key)
    states = _parse_detail_array(row, "snp_states", np.int8, key)
    prob1 = _parse_detail_array(row, "snp_prob_label1", np.float64, key)
    prob2 = _parse_detail_array(row, "snp_prob_label2", np.float64, key)

    lengths = {len(positions), len(states), len(prob1), len(prob2)}
    if len(lengths) != 1 or not positions.size:
        raise ValueError(
            f"Inconsistent or empty SNP-detail arrays for segment {key}: "
            f"positions={len(positions)} states={len(states)} "
            f"prob1={len(prob1)} prob2={len(prob2)}"
        )

    row_min = int(positions.min())
    row_max = int(positions.max())
    accumulator.min_pos = (
        row_min if accumulator.min_pos is None else min(accumulator.min_pos, row_min)
    )
    accumulator.max_pos = (
        row_max if accumulator.max_pos is None else max(accumulator.max_pos, row_max)
    )

    label1 = states == 1
    label2 = states == 2
    n1 = int(label1.sum())
    n2 = int(label2.sum())
    accumulator.n1 += n1
    accumulator.n2 += n2
    if n1:
        accumulator.prob1_sum += float(prob1[label1].sum(dtype=np.float64))
    if n2:
        accumulator.prob2_sum += float(prob2[label2].sum(dtype=np.float64))
    accumulator.detail_rows += 1


def merge_segments_exact_streaming(
    df,
    snp_details_path,
    merge_distance,
    min_snps_per_segment=2,
    mosaic_minority_threshold=0.20,
    progress_rows=250_000,
):
    """Exact-merge segments while streaming the gzip SNP-detail table.

    The legacy implementation materialized the full gzip table as a pandas
    DataFrame and then created four NumPy arrays for every raw segment. Large
    chromosomes therefore required tens of GiB and spent hours reclaiming
    cgroup memory. This implementation plans merge membership from the already
    filtered BED, scans the gzip once, and reduces each needed detail row into
    numeric group totals immediately.
    """
    if df.empty:
        return pd.DataFrame()
    if progress_rows <= 0:
        raise ValueError("progress_rows must be positive")
    if not os.path.exists(snp_details_path):
        raise FileNotFoundError(f"SNP details file not found: {snp_details_path}")

    key_to_group, accumulators = _build_exact_merge_plan(df, merge_distance)
    required_count = len(key_to_group)
    matched_count = 0
    scanned_count = 0
    required_columns = {
        "chr",
        "sample_hap_id",
        "start_pos",
        "end_pos",
        "snp_positions",
        "snp_states",
        "snp_prob_label1",
        "snp_prob_label2",
    }

    logger.info(
        "Streaming SNP details from %s for %d retained segments",
        snp_details_path,
        required_count,
    )
    csv_limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(csv_limit)
            break
        except OverflowError:
            csv_limit //= 10
    with gzip.open(snp_details_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = required_columns.difference(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(
                "SNP details file is missing columns: "
                + ", ".join(sorted(missing_columns))
            )

        for scanned_count, row in enumerate(reader, start=1):
            key = _segment_key(
                row["chr"],
                row["sample_hap_id"],
                row["start_pos"],
                row["end_pos"],
            )
            group_index = key_to_group.pop(key, None)
            if group_index is not None:
                _accumulate_detail_row(accumulators[group_index], row, key)
                matched_count += 1

            if scanned_count % progress_rows == 0:
                logger.info(
                    "Streaming SNP details: scanned=%d matched=%d remaining=%d",
                    scanned_count,
                    matched_count,
                    len(key_to_group),
                )

    missing_count = len(key_to_group)
    logger.info(
        "Finished streaming SNP details: scanned=%d matched=%d required=%d missing=%d",
        scanned_count,
        matched_count,
        required_count,
        missing_count,
    )
    if missing_count:
        examples = list(key_to_group)[:10]
        logger.warning(
            "SNP details were missing for %d retained segments; first keys: %s",
            missing_count,
            examples,
        )
    key_to_group.clear()

    records = []
    for accumulator in accumulators:
        total_snps = accumulator.n1 + accumulator.n2
        if accumulator.detail_rows == 0 or total_snps < min_snps_per_segment:
            continue

        if accumulator.n1 and accumulator.n2:
            minority_ratio = min(accumulator.n1, accumulator.n2) / total_snps
            if minority_ratio >= mosaic_minority_threshold:
                label = 3
            else:
                label = 1 if accumulator.n1 >= accumulator.n2 else 2
        elif accumulator.n1:
            label = 1
        elif accumulator.n2:
            label = 2
        else:
            continue

        if label == 1:
            avg_prob = accumulator.prob1_sum / accumulator.n1
        elif label == 2:
            avg_prob = accumulator.prob2_sum / accumulator.n2
        else:
            avg_prob = (
                accumulator.prob1_sum + accumulator.prob2_sum
            ) / total_snps

        records.append(
            (
                accumulator.chr,
                accumulator.min_pos,
                accumulator.max_pos,
                accumulator.haplotype,
                label,
                total_snps,
                avg_prob,
                accumulator.n1,
                accumulator.n2,
                accumulator.sample_hap_id,
            )
        )

    if not records:
        return pd.DataFrame()

    result_df = pd.DataFrame.from_records(
        records,
        columns=[
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
        ],
    )
    result_df = result_df.sort_values(["sample_hap_id", "chr", "start_pos"])
    logger.info("Streaming exact merge: %d -> %d segments", len(df), len(result_df))
    return result_df


def merge_segments_exact(df, snp_details_dict, merge_distance,
                         min_snps_per_segment=2,
                         mosaic_minority_threshold=0.20):
    """
    Exact merge replication using SNP-level details

    This function 100% replicates the main program's merge logic:
    1. Merges adjacent segments regardless of label (not grouped by label!)
    2. Recalculates label based on combined SNP counts
    3. Recalculates avg_prob from original SNP probabilities
    4. Applies min_snps_per_segment filter after merging

    Args:
        df: Raw BED DataFrame
        snp_details_dict: SNP details lookup dictionary
        merge_distance: Maximum distance (bp) to merge segments
        min_snps_per_segment: Minimum SNPs after merging
        mosaic_minority_threshold: Threshold for Mosaic classification

    Returns:
        DataFrame: Merged segments
    """
    if df.empty:
        return df

    merged_segments = []

    # Group by sample_hap_id only (NOT by label!)
    for sample_hap_id, group in df.groupby('sample_hap_id'):
        group = group.sort_values('start_pos').copy()

        if len(group) == 0:
            continue

        # Start with first segment
        current_seg_rows = [group.iloc[0]]

        for i in range(1, len(group)):
            next_seg = group.iloc[i]
            current_end = current_seg_rows[-1]['end_pos']
            gap = next_seg['start_pos'] - current_end

            if gap <= merge_distance:
                # Merge: add to current accumulated segments
                current_seg_rows.append(next_seg)
            else:
                # Don't merge: process accumulated segments
                merged_seg = _process_merged_segment_exact(
                    current_seg_rows, snp_details_dict,
                    min_snps_per_segment, mosaic_minority_threshold
                )
                if merged_seg is not None:
                    merged_segments.append(merged_seg)

                # Start new accumulation
                current_seg_rows = [next_seg]

        # Process last accumulated segments
        merged_seg = _process_merged_segment_exact(
            current_seg_rows, snp_details_dict,
            min_snps_per_segment, mosaic_minority_threshold
        )
        if merged_seg is not None:
            merged_segments.append(merged_seg)

    if not merged_segments:
        return pd.DataFrame()

    result_df = pd.DataFrame(merged_segments)
    result_df = result_df.sort_values(['sample_hap_id', 'chr', 'start_pos'])

    logger.info(f"Exact merge: {len(df)} -> {len(result_df)} segments")
    return result_df


def _process_merged_segment_exact(seg_rows, snp_details_dict,
                                   min_snps_per_segment,
                                   mosaic_minority_threshold):
    """
    Process merged segment - exact replication of main program logic

    Args:
        seg_rows: List of segment rows (pd.Series) to be merged
        snp_details_dict: SNP details lookup
        min_snps_per_segment: Min SNPs threshold
        mosaic_minority_threshold: Mosaic threshold

    Returns:
        dict or None: Merged segment data
    """
    # Collect all SNPs from all segments being merged
    all_positions = []
    all_states = []
    all_prob1 = []
    all_prob2 = []

    for seg_row in seg_rows:
        # Lookup SNP details using composite key
        key = (seg_row['chr'], seg_row['sample_hap_id'],
               seg_row['start_pos'], seg_row['end_pos'])

        if key not in snp_details_dict:
            logger.warning(f"SNP details not found for segment: {key}")
            continue

        snp_data = snp_details_dict[key]
        all_positions.extend(snp_data['positions'].tolist())
        all_states.extend(snp_data['states'].tolist())
        all_prob1.extend(snp_data['prob_label1'].tolist())
        all_prob2.extend(snp_data['prob_label2'].tolist())

    if len(all_positions) == 0:
        return None

    # Convert to numpy arrays
    all_positions = np.array(all_positions)
    all_states = np.array(all_states)
    all_prob1 = np.array(all_prob1)
    all_prob2 = np.array(all_prob2)

    # Count SNPs by state
    n1 = np.sum(all_states == 1)
    n2 = np.sum(all_states == 2)
    total_snps = n1 + n2

    # Apply min_snps filter (EXACTLY as in main program)
    if total_snps < min_snps_per_segment:
        return None

    # Re-judge label (EXACTLY as in main program)
    label = 0
    if n1 > 0 and n2 > 0:
        minority_ratio = min(n1, n2) / total_snps
        if minority_ratio >= mosaic_minority_threshold:
            label = 3  # Mosaic
        else:
            label = 1 if n1 >= n2 else 2
    elif n1 > 0:
        label = 1
    elif n2 > 0:
        label = 2

    if label == 0:
        return None

    # Recalculate avg_prob (EXACTLY as in main program)
    prob_values = []
    if label == 1:
        # Only use prob_label1 for SNPs with state==1
        prob_values = all_prob1[all_states == 1]
    elif label == 2:
        # Only use prob_label2 for SNPs with state==2
        prob_values = all_prob2[all_states == 2]
    elif label == 3:
        # Use both: prob_label1 for state==1, prob_label2 for state==2
        prob_values = np.concatenate([
            all_prob1[all_states == 1],
            all_prob2[all_states == 2]
        ])

    avg_prob = np.mean(prob_values) if len(prob_values) > 0 else np.nan

    # Build merged segment (10 columns, NO cM)
    first_seg = seg_rows[0]
    last_seg = seg_rows[-1]

    return {
        'chr': first_seg['chr'],
        'start_pos': int(np.min(all_positions)),
        'end_pos': int(np.max(all_positions)),
        'haplotype': first_seg['haplotype'],
        'ancestry_label': label,
        'num_snps': total_snps,
        'avg_prob': avg_prob,
        'archaic_snps': n1,
        'african_snps': n2,
        'sample_hap_id': first_seg['sample_hap_id']
    }


# ============================================================================
# Approximate merge functions (used when SNP details file is not available)
# ============================================================================

def read_bed_file(filepath):
    """
    Read a 10-column BED file without header (NO cM columns)

    Args:
        filepath: Path to BED file

    Returns:
        DataFrame with named columns
    """
    column_names = [
        'chr', 'start_pos', 'end_pos',
        'haplotype', 'ancestry_label', 'num_snps', 'avg_prob',
        'archaic_snps', 'african_snps', 'sample_hap_id'
    ]

    try:
        df = pd.read_csv(filepath, sep='\t', header=None, names=column_names)
        logger.info(f"Read {len(df)} segments from {filepath}")
        return df
    except Exception as e:
        logger.error(f"Error reading BED file: {e}")
        sys.exit(1)


def filter_segments(df, min_snps=None, min_prob=None, max_label=None):
    """
    Filter segments by various criteria

    Args:
        df: Input DataFrame
        min_snps: Minimum number of SNPs (optional)
        min_prob: Minimum average probability (optional)
        max_label: Maximum ancestry label to keep (optional, e.g., 2 to exclude Mosaic)

    Returns:
        Filtered DataFrame
    """
    original_count = len(df)

    if min_snps is not None:
        df = df[df['num_snps'] >= min_snps]
        logger.info(f"Filtered by min_snps={min_snps}: {original_count} -> {len(df)}")
        original_count = len(df)

    if min_prob is not None:
        df = df[df['avg_prob'] >= min_prob]
        logger.info(f"Filtered by min_prob={min_prob}: {original_count} -> {len(df)}")
        original_count = len(df)

    if max_label is not None:
        df = df[df['ancestry_label'] <= max_label]
        logger.info(f"Filtered by max_label={max_label}: {original_count} -> {len(df)}")

    return df


def merge_segments(df, merge_distance, recalculate_stats=True):
    """
    Merge nearby segments of the same ancestry label (approximate method)

    Args:
        df: Input DataFrame with segments
        merge_distance: Maximum distance (bp) to merge segments
        recalculate_stats: Whether to recalculate num_snps and avg_prob

    Returns:
        DataFrame with merged segments
    """
    if df.empty:
        return df

    # Group by sample_hap_id and ancestry_label
    merged_segments = []

    for (sample_hap_id, label), group in df.groupby(['sample_hap_id', 'ancestry_label']):
        # Sort by start position
        group = group.sort_values('start_pos').copy()

        if len(group) == 0:
            continue

        # Initialize first segment
        current_seg = group.iloc[0].to_dict()

        for i in range(1, len(group)):
            next_seg = group.iloc[i].to_dict()
            gap = next_seg['start_pos'] - current_seg['end_pos']

            if gap <= merge_distance:
                # Merge: extend current segment
                current_seg['end_pos'] = max(current_seg['end_pos'], next_seg['end_pos'])

                if recalculate_stats:
                    # Sum SNP counts
                    current_seg['num_snps'] += next_seg['num_snps']
                    current_seg['archaic_snps'] += next_seg['archaic_snps']
                    current_seg['african_snps'] += next_seg['african_snps']

                    # Weighted average of probabilities
                    total_snps = current_seg['num_snps']
                    current_seg['avg_prob'] = (
                        current_seg['avg_prob'] * (total_snps - next_seg['num_snps']) +
                        next_seg['avg_prob'] * next_seg['num_snps']
                    ) / total_snps if total_snps > 0 else current_seg['avg_prob']
            else:
                # Gap too large, save current and start new
                merged_segments.append(current_seg)
                current_seg = next_seg.copy()

        # Add the last segment
        merged_segments.append(current_seg)

    if not merged_segments:
        return pd.DataFrame()

    # Convert back to DataFrame
    merged_df = pd.DataFrame(merged_segments)

    # Ensure correct column order (10 columns, NO cM)
    column_order = [
        'chr', 'start_pos', 'end_pos',
        'haplotype', 'ancestry_label', 'num_snps', 'avg_prob',
        'archaic_snps', 'african_snps', 'sample_hap_id'
    ]
    merged_df = merged_df[column_order]

    # Sort by sample_hap_id, chr, start_pos
    merged_df = merged_df.sort_values(['sample_hap_id', 'chr', 'start_pos'])

    logger.info(f"Merged segments: {len(df)} -> {len(merged_df)} (merge_distance={merge_distance})")

    return merged_df


def write_bed_file(df, filepath):
    """
    Write BED file without header

    Args:
        df: DataFrame to write
        filepath: Output file path
    """
    try:
        df.to_csv(filepath, sep='\t', index=False, header=False)
        logger.info(f"Written {len(df)} segments to {filepath}")
    except Exception as e:
        logger.error(f"Error writing BED file: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Merge BED segments from ArchaicSeeker3 output (No Genetic Map version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic merge with 10kb distance
  python merge_bed_segments.py -i raw.bed -o merged.bed -d 10000

  # Filter before merging
  python merge_bed_segments.py -i raw.bed -o merged.bed -d 10000 --min-snps 5 --min-prob 0.8

  # Exclude Mosaic segments (label=3)
  python merge_bed_segments.py -i raw.bed -o merged.bed -d 10000 --max-label 2

  # No merge (just filter and reformat)
  python merge_bed_segments.py -i raw.bed -o filtered.bed -d 0 --min-snps 10
        '''
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input BED file (raw segments)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output BED file (merged segments)')
    parser.add_argument('-d', '--merge-distance', type=int, default=10000,
                        help='Maximum distance (bp) to merge segments (default: 10000)')

    # Filtering options
    parser.add_argument('--min-snps', type=int, default=None,
                        help='Filter: minimum number of SNPs per segment')
    parser.add_argument('--min-prob', type=float, default=None,
                        help='Filter: minimum average probability per segment')
    parser.add_argument('--max-label', type=int, default=None,
                        help='Filter: maximum ancestry label (e.g., 2 to exclude Mosaic=3)')

    # Advanced options
    parser.add_argument('--no-recalculate', action='store_true',
                        help='Do not recalculate statistics after merging')
    parser.add_argument('--verbose', action='store_true',
                        help='Verbose output')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("=" * 60)
    logger.info("BED Segment Merger (No Genetic Map Version)")
    logger.info("=" * 60)

    # Read input
    df = read_bed_file(args.input)

    # Check for SNP details file for exact merge mode
    snp_details_path = args.input.replace('.raw.bed', '.raw.snps.gz')
    use_exact_merge = False

    if os.path.exists(snp_details_path):
        logger.info("=" * 60)
        logger.info("EXACT MERGE MODE ENABLED")
        logger.info(f"Found SNP details file: {snp_details_path}")
        logger.info("Will use 100% exact merge replication")
        logger.info("=" * 60)
        use_exact_merge = True
    else:
        logger.info("=" * 60)
        logger.info("APPROXIMATE MERGE MODE")
        logger.info(f"SNP details file not found: {snp_details_path}")
        logger.info("Using approximate merge (results may differ slightly from main program)")
        logger.info("=" * 60)

    # Filter if requested (only in approximate mode)
    if not use_exact_merge:
        if args.min_snps or args.min_prob or args.max_label:
            logger.info("Applying filters...")
            df = filter_segments(df, args.min_snps, args.min_prob, args.max_label)
    else:
        if args.min_snps or args.min_prob or args.max_label:
            logger.warning("Filtering options are ignored in exact merge mode")
            logger.warning("Exact merge applies fixed filters as in main program")

    # Merge
    if args.merge_distance > 0:
        logger.info(f"Merging segments (merge_distance={args.merge_distance})...")
        if use_exact_merge:
            # EXACT MERGE: Use exact replication
            df = merge_segments_exact_streaming(
                df, snp_details_path, args.merge_distance,
                min_snps_per_segment=2,
                mosaic_minority_threshold=0.20
            )
        else:
            # Approximate merge
            df = merge_segments(df, args.merge_distance, recalculate_stats=not args.no_recalculate)
    else:
        logger.info("Merge distance is 0, skipping merge step")

    # Write output
    write_bed_file(df, args.output)

    logger.info("=" * 60)
    logger.info("Done!")
    logger.info("=" * 60)

    # Summary
    logger.info(f"Input:  {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Total segments: {len(df)}")


if __name__ == '__main__':
    main()
