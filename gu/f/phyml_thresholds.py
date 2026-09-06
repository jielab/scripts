"""Audited threshold compatibility with Bald/archaic introgression/locus.R.

Only numerical gates change. Region discovery, derived-state diagnostics and
candidate-focused tree confirmation remain the GU workflow.
"""
from __future__ import annotations
import math

PROFILES = {
    "strict": dict(min_diagnostic_sites=3, min_diagnostic_match_prop=0.8,
                   min_candidate_copies=2, max_diagnostic_gap_bp=50000,
                   min_candidate_purity=0.8, min_candidate_sensitivity=0.5,
                   max_ils_probability=1.0),
    "legacy_compatible": dict(min_diagnostic_sites=1, min_diagnostic_match_prop=0.0,
                              min_candidate_copies=11, max_diagnostic_gap_bp=0,
                              min_candidate_purity=0.0, min_candidate_sensitivity=0.0,
                              max_ils_probability=0.1),
}


def apply_profile(args):
    for key, value in PROFILES[args.threshold_profile].items():
        if getattr(args, key, None) is None:
            setattr(args, key, value)


def ils_probability(length_bp, recomb_cm_mb=0.53, split_years=550000,
                    archaic_age_years=50000, generation_years=29):
    """Gamma(shape=2) survival; same model as old ils_p, without 1-CDF cancellation.

    This is a model-based screen, not a calibrated introgression probability.
    The caller supplies the observed diagnostic span, never the flanked region.
    """
    x = max(0, length_bp) * recomb_cm_mb * 1e-8 * (
        (2 * split_years - archaic_age_years) / generation_years)
    return math.exp(-x) * (1 + x)
