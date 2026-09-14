"""ILS probability used by the GWAS risk-core workflow.
"""
from __future__ import annotations
import math


# 🚩 Incomplete lineage sorting probability
def ils_probability(length_bp, recomb_cm_mb=0.53, split_years=550000,
                    archaic_age_years=50000, generation_years=29):
    """Gamma(shape=2) survival; same model as old ils_p, without 1-CDF cancellation.

    This is a model-based screen, not a calibrated introgression probability.
    The caller supplies the observed diagnostic span, never the flanked region.
    """
    x = max(0, length_bp) * recomb_cm_mb * 1e-8 * (
        (2 * split_years - archaic_age_years) / generation_years)
    return math.exp(-x) * (1 + x)
