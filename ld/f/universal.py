"""Consensus cuts from McManus et al., Methods S1, pp. e6-e7.

Input intervals use 1-based, closed coordinates at retained variants.  Gaps
between ancestry blocks are intentional (Methods S1, p. e1); they are never
filled before computing the consensus.

The paper calls s_b a ``start point'' and e_b an ``endpoint''.  This independent
implementation interprets the candidate endpoint domain as the right ends
{e_b}.  An empty C_p has no maximum, so it generates no candidate.  Both choices
are explicit because the restricted author code cannot resolve the terminology
or its degenerate/singleton cases.  The audit also computes the alternative
start-and-end, nonempty-C_p interpretation, without choosing by target counts.

Zero-width cuts are valid.  A zero-width cut at p belongs to the block on its
left; a positive-width cut [s,e] has the left core ending at s and the right
core beginning at e.  Its interior is left unassigned, as in the paper's
``regions intervening the optimal cuts''.  This is not a gap-filled partition.
"""

from bisect import bisect_left
from pathlib import Path
from typing import Mapping, Sequence

METHOD_ID = "mcmanus2023-right-endpoints-nonempty-cp-v1"


def read_bed(path, chrom):
    """Read BED into strict, nonoverlapping 1-based closed support intervals."""
    rows = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith(("#", "track ", "browser ")):
            continue
        fields = line.split()
        if len(fields) < 3:
            raise ValueError(f"Malformed BED row {number}: {path}")
        c, start, end = fields[:3]
        start, end = int(start), int(end)
        if c.removeprefix("chr") != str(chrom).removeprefix("chr") or start < 0 or end <= start:
            raise ValueError(f"Invalid BED row {number}: {path}")
        if rows and start < rows[-1][1]:
            raise ValueError(f"Unsorted/overlapping BED row {number}: {path}")
        rows.append((start + 1, end))
    if not rows:
        raise ValueError(f"Empty BED: {path}")
    return rows


def _validate_blocks(pop_blocks):
    if not pop_blocks:
        raise ValueError("At least one ancestry partition is required")
    for population, rows in pop_blocks.items():
        if not rows:
            raise ValueError(f"Empty ancestry partition: {population}")
        previous_end = None
        for start, end in rows:
            if not isinstance(start, int) or not isinstance(end, int) or start < 1 or end < start:
                raise ValueError(f"Invalid closed block for {population}: {(start, end)}")
            if previous_end is not None and start <= previous_end:
                raise ValueError(f"Unsorted/overlapping blocks for {population}")
            previous_end = end


def eligible_candidates(blocks, beta=2000, endpoint_domain="right"):
    """Enumerate the narrowest candidate for every supported endpoint.

    A prefix maximum is sufficient: if it is less than p, every previous block
    has expired and C_p is empty.  Otherwise that maximum is exactly max(C_p).
    Strict s_b < p is essential; starts at p must not enter C_p.
    """
    if not isinstance(beta, int) or beta < 0:
        raise ValueError("beta must be a nonnegative integer")
    if endpoint_domain not in ("right", "both"):
        raise ValueError("endpoint_domain must be 'right' or 'both'")
    blocks = sorted(blocks)
    if not blocks:
        raise ValueError("No input blocks")
    points = sorted({e for _, e in blocks} if endpoint_domain == "right"
                    else {p for block in blocks for p in block})
    candidates, empty = [], []
    index, max_end = 0, None
    for point in points:
        while index < len(blocks) and blocks[index][0] < point:
            end = blocks[index][1]
            max_end = end if max_end is None else max(max_end, end)
            index += 1
        if max_end is None or max_end < point:
            empty.append(point)
            continue
        if max_end - point <= beta:
            candidates.append((point, max_end))
    return candidates, empty


def optimize_cuts(candidates):
    """Maximize disjoint cut count, then minimize total cut width exactly.

    Closed cuts must not touch.  Stable sorting and retaining the earlier
    solution on a complete objective tie make the result deterministic.
    """
    candidates = sorted(set(candidates), key=lambda cut: (cut[1], cut[0]))
    ends = [end for _, end in candidates]
    scores, parents, take = [(0, 0)], [], []
    for index, (start, end) in enumerate(candidates):
        if end < start:
            raise ValueError("A cut cannot have negative width")
        previous = bisect_left(ends, start, 0, index)
        proposed = (scores[previous][0] + 1, scores[previous][1] - (end - start))
        chosen = proposed > scores[-1]
        scores.append(proposed if chosen else scores[-1])
        parents.append(previous)
        take.append(chosen)
    selected, index = [], len(candidates)
    while index:
        if take[index - 1]:
            selected.append(candidates[index - 1])
            index = parents[index - 1]
        else:
            index -= 1
    return selected[::-1]


def select_cuts(blocks, beta=2000):
    """Public default selection; excludes the left terminal marker."""
    candidates, _ = eligible_candidates(blocks, beta, "right")
    return optimize_cuts(candidates)


def cuts_to_bed(cuts, start_1based, end_1based):
    """Export nonoverlapping BED cores without inventing extra boundaries."""
    if start_1based < 1 or end_1based < start_1based:
        raise ValueError("Invalid terminal range")
    rows = []
    current = start_1based
    previous_end = None
    for start, end in cuts:
        if not start_1based <= start <= end <= end_1based:
            raise ValueError("Cut outside retained coordinate range")
        if previous_end is not None and start <= previous_end:
            raise ValueError("Closed cuts overlap or touch")
        if start >= current:
            rows.append((current - 1, start))
        current = end + 1 if start == end else end
        previous_end = end
    if current <= end_1based:
        rows.append((current - 1, end_1based))
    if not rows:
        raise ValueError("Consensus produced no nonempty blocks")
    for index, (start, end) in enumerate(rows):
        if end <= start or (index and start < rows[index - 1][1]):
            raise ValueError("Consensus BED is empty or overlapping")
    return rows


def compute_consensus(pop_blocks: Mapping[str, Sequence[tuple[int, int]]], beta=2000):
    """Return BED cores, selected cuts, terminal metadata and audit diagnostics."""
    _validate_blocks(pop_blocks)
    blocks = [block for rows in pop_blocks.values() for block in rows]
    lo = min(start for start, _ in blocks)
    hi = max(end for _, end in blocks)
    candidates, empty = eligible_candidates(blocks, beta, "right")
    cuts = optimize_cuts(candidates)
    rows = cuts_to_bed(cuts, lo, hi)
    alternative, alternative_empty = eligible_candidates(blocks, beta, "both")
    alternative_cuts = optimize_cuts(alternative)
    alternative_rows = cuts_to_bed(alternative_cuts, lo, hi)
    for start, end in cuts:
        if end - start > beta:
            raise ValueError("Cut exceeds beta")
        if any(left < start and right > end for left, right in blocks):
            raise ValueError("A consensus cut bisects an input block")
        if not any(left < start <= right for left, right in blocks):
            raise ValueError("A consensus cut has an empty spanning set")
    diagnostics = {
        "method_id": METHOD_ID,
        "endpoint_domain": "right_ends_e_b",
        "empty_spanning_set_policy": "skip_no_maximum_defined",
        "zero_width_cut_assignment": "left_block",
        "positive_width_cut_policy": "retain_endpoints_exclude_interior",
        "beta_bp": beta,
        "input_counts": {pop: len(partition) for pop, partition in pop_blocks.items()},
        "input_singleton_blocks": sum(start == end for start, end in blocks),
        "candidate_count": len(candidates),
        "empty_right_endpoints_skipped": len(empty),
        "empty_right_endpoint_positions": empty,
        "selected_cut_count": len(cuts),
        "internal_zero_width_cuts": sum(start == end and lo < start < hi for start, end in cuts),
        "empty_spanning_selected_cuts": 0,
        "total_cut_width_bp": sum(end - start for start, end in cuts),
        "excluded_cut_interior_bp": sum(max(0, end - start - 1) for start, end in cuts),
        "universal_blocks": len(rows),
        "short_blocks_le_beta": sum(end - start <= beta for start, end in rows),
        "minimum_block_bp": min(end - start for start, end in rows),
        "start_1based": lo,
        "end_1based": hi,
        "both_endpoints_nonempty_cp_blocks": len(alternative_rows),
        "both_endpoints_empty_candidates_skipped": len(alternative_empty),
        "endpoint_interpretations_identical": cuts == alternative_cuts,
        "structurally_valid": True,
        "author_coordinate_identity_verified": False,
    }
    return {"bed": rows, "cuts": cuts, "start_1based": lo, "end_1based": hi,
            "diagnostics": diagnostics}
