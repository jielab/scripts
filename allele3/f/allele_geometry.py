#!/usr/bin/env python3
"""Within-site geometry of three reference-relative scalar AVI predictions.

This module is called after the existing chromosome/position/REF/ALT join.
It neither reads the Atlas nor stores sites. AF ranks label observed ALTs in
one population, never mutation age or rank among all four DNA bases.
"""
from __future__ import annotations

import itertools
import math

BASES = frozenset("ACGT")
SCALES = ("RAW", "PHRED")
PERMUTATION_INTERPRETATION = "descriptive_rank_exchangeability_not_site_discovery"
GEOMETRY_FIELDS = [
    "chrom", "pos", "ref", "population", "scale", "geometry_status",
    "count_mode", "n_observed_alt", "n_qualified_alt", "n_ranked_alt",
    "ranked_alts", "ranked_afs", "ranked_scores", "prediction_alts",
    "prediction_scores", "score_range", "distance_top_middle",
    "distance_top_least", "distance_middle_least", "weak_pair_distance_rank",
    "normalized_weak_pair_contrast", "normalized_least_alt_contrast",
    "weak_pair_permutation_p", "least_alt_permutation_p", "permutation_assignments",
    "permutation_interpretation", "unreported_alt", "unreported_evidence",
    "normalized_observed_unreported_contrast",
]


def finite(value):
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def packed(values):
    return ",".join("NA" if value is None else format(value, ".12g") if isinstance(value, (float, int)) else str(value) for value in values)


def le(left, right):
    return left <= right or math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-14)


def three_scores(values):
    """Return distances and two contrasts; caller requires a positive range."""
    top, middle, least = values
    span = max(values) - min(values)
    d12, d13, d23 = abs(top-middle), abs(top-least), abs(middle-least)
    weak = (d23-(d12+d13)/2) / span
    least_contrast = (d23-d13) / span
    distances = (d12/span, d13/span, d23/span)
    selected = d23/span
    lower = sum(d < selected and not math.isclose(d, selected, rel_tol=1e-12, abs_tol=1e-14) for d in distances)
    equal = sum(math.isclose(d, selected, rel_tol=1e-12, abs_tol=1e-14) for d in distances)
    return d12, d13, d23, weak, least_contrast, lower+(equal+1)/2


def analyze_site(row, predicted, qualified_indices):
    """Return ``(geometry_rows, contrasts)`` for one exact matched SNV site.

    ``predicted`` maps ALT base -> {"RAW": score, "PHRED": score}; REF has no
    scalar score. ``qualified_indices`` is the main analysis's threshold-qualified
    ALT index list. Each contrast has the existing site's numerator/denominator
    shape, so callers can reuse their genomic-block bootstrap and FDR family.

    Six score-to-rank permutations are descriptive exchangeability diagnostics.
    Mutation classes and ascertainment can break that exchangeability; these
    values are not per-site discoveries or a mutation/selection null model.
    """
    ref = row["ref"]
    alts = row["alts"]
    if ref not in BASES or len(alts) != len(set(alts)) or any(b not in BASES-{ref} for b in alts):
        raise ValueError("Allele geometry requires distinct non-REF A/C/G/T bases")
    if any(b not in BASES-{ref} for b in predicted):
        raise ValueError("AVI geometry cannot assign a prediction to REF or a non-SNV ALT")
    indices = list(qualified_indices)
    if len(indices) != len(set(indices)) or any(i < 0 or i >= len(alts) for i in indices):
        raise ValueError("Invalid qualified ALT indices")
    afs = [finite(x) for x in row["alt_afs"]]
    source = row["alt_afs"] if row["count_mode"] == "AF_only" else row["alt_acs"]
    evidence = [finite(x) for x in source]
    observed = {i for i, value in enumerate(evidence) if value is not None and value > 0}
    qualified = [i for i in indices if i in observed]
    ranked = sorted((i for i in qualified if afs[i] is not None and afs[i] > 0), key=lambda i: (-afs[i], alts[i]))
    ties = len({afs[i] for i in ranked}) != len(ranked)
    possible = sorted(BASES-{ref})
    unreported = [b for b in possible if b not in {alts[i] for i in observed}]
    third_base = unreported[0] if len(observed) == 2 and len(unreported) == 1 else None
    third_evidence = None
    if third_base is not None:
        if third_base not in alts:
            third_evidence = "not_reported"
        else:
            value = evidence[alts.index(third_base)]
            third_evidence = "reported_missing" if value is None else "reported_zero" if value == 0 else "invalid_evidence"
    output, contrasts = [], []
    for scale in SCALES:
        scores = {b: finite(predicted.get(b, {}).get(scale)) for b in possible}
        complete = all(scores[b] is not None for b in possible)
        span = max(scores.values())-min(scores.values()) if complete else None
        result = dict.fromkeys(GEOMETRY_FIELDS)
        result.update(chrom=row["chrom"], pos=row["pos"], ref=ref, population=row.get("population", "ALL"),
            scale=scale, count_mode=row["count_mode"], n_observed_alt=len(observed), n_qualified_alt=len(qualified),
            n_ranked_alt=len(ranked), ranked_alts=packed([alts[i] for i in ranked]), ranked_afs=packed([afs[i] for i in ranked]),
            ranked_scores=packed([scores[alts[i]] for i in ranked]), prediction_alts=packed(possible),
            prediction_scores=packed([scores[b] for b in possible]), score_range=span,
            unreported_alt=third_base, unreported_evidence=third_evidence,
            permutation_assignments=0, permutation_interpretation=PERMUTATION_INTERPRETATION)
        if not complete:
            result["geometry_status"] = "incomplete_three_ALT_predictions"
        elif span == 0:
            result["geometry_status"] = "equal_three_ALT_predictions"
        elif len(ranked) < 2:
            result["geometry_status"] = "insufficient_qualified_positive_AF"
        elif ties:
            result["geometry_status"] = "AF_rank_tie"
        elif len(observed) == len(qualified) == len(ranked) == 3:
            result["geometry_status"] = "three_observed_ranked_ALTs"
            values = [scores[alts[i]] for i in ranked]
            d12,d13,d23,weak,least,pair_rank = three_scores(values)
            permutations = [three_scores(permutation) for permutation in itertools.permutations(values)]
            # Keep all six assignments even when score values tie. Every mapping
            # of the three nucleotide identities is a separate rank assignment.
            result.update(distance_top_middle=d12, distance_top_least=d13, distance_middle_least=d23,
                weak_pair_distance_rank=pair_rank, normalized_weak_pair_contrast=weak,
                normalized_least_alt_contrast=least,
                weak_pair_permutation_p=sum(le(item[3],weak) for item in permutations)/6,
                least_alt_permutation_p=sum(le(item[4],least) for item in permutations)/6,
                permutation_assignments=6)
            for metric,value in [("normalized_weak_pair_distance_minus_other_pairs",weak),
                                 ("normalized_least_ALT_distance_to_middle_minus_top",least)]:
                contrasts.append(dict(scale=scale, metric=metric, numerator=value, denominator=1.))
        elif len(observed) == len(qualified) == len(ranked) == 2 and third_base is not None:
            result["geometry_status"] = "two_observed_ranked_ALTs"
            top, extra = [scores[alts[i]] for i in ranked]
            contrast = (abs(extra-top)-abs(scores[third_base]-top))/span
            result["normalized_observed_unreported_contrast"] = contrast
            contrasts.append(dict(scale=scale, metric="normalized_extra_ALT_distance_to_top_minus_unreported",
                numerator=contrast, denominator=1.))
        else:
            result["geometry_status"] = "observed_ALT_below_threshold_or_missing_AF"
        output.append(result)
    return output, contrasts
