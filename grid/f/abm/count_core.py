#!/usr/bin/env python3
"""All-in-radius reference counting, with no fitted prediction/error head.

Every eligible in-radius donor contributes one unit.  No nearest-neighbour cap,
outcome-based donor pruning, query phenotype, or fitted regression enters the
prediction.  Complete sparse memberships can be streamed to local-only files.
"""
from __future__ import annotations

# Core

from dataclasses import dataclass
from statistics import NormalDist
import time

import numpy as np
import pandas as pd

from abm.support import _validate_inputs


@dataclass
class CountResult:
    table: pd.DataFrame
    device_used: str


def validate_outcomes(y, task):
    y = np.asarray(y, dtype=float)
    if y.ndim != 1 or not len(y) or not np.isfinite(y).all():
        raise ValueError("Donors require finite observed outcomes")
    if task not in {"continuous", "binary"}:
        raise ValueError("task must be continuous or binary")
    if task == "binary" and not np.isin(y, [0., 1.]).all():
        raise ValueError("Binary outcomes must be 0/1; censored endpoints need an explicit survival adapter")
    return y


def summarize_neighbors(y, groups, task, *, alpha=.05):
    """Return arithmetic sufficient statistics and descriptive uncertainty.

    The mean interval is a normal approximation using cluster sums by declared
    donor group. It is not an individual prediction interval or a conditional
    coverage guarantee. Binary outcome intervals use Wilson with group ESS as
    a conservative descriptive effective count, not an exact binomial model.
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    empty = dict(n_twins=0, n_groups=0, ess_group=0., outcome_sum=0.,
                 p_count=np.nan, neighbor_variance=np.nan,
                 neighbor_sd=np.nan, donor_q025=np.nan, donor_q975=np.nan,
                 mean_se_cluster=np.nan, mean_low=np.nan, mean_high=np.nan,
                 n_cases=0 if task == "binary" else np.nan,
                 n_controls=0 if task == "binary" else np.nan,
                 label_disagreement=np.nan)
    if not n:
        return empty
    _, inv, counts = np.unique(np.asarray(groups, str), return_inverse=True, return_counts=True)
    g = len(counts)
    ess = float(n * n / np.square(counts.astype(float)).sum())
    mean = float(y.mean())
    variance = float(np.mean((y - mean) ** 2))
    se = np.nan
    if g > 1:
        centered_sums = np.bincount(inv, weights=y - mean)
        se = float(np.sqrt(g / (g - 1) * np.square(centered_sums).sum() / n**2))
    z = NormalDist().inv_cdf(1 - alpha / 2)
    low, high = mean - z * se, mean + z * se
    if task == "binary":
        denominator = 1 + z*z/ess
        center = (mean + z*z/(2*ess)) / denominator
        width = z * np.sqrt(mean*(1-mean)/ess + z*z/(4*ess*ess)) / denominator
        low, high = center-width, center+width
    return dict(n_twins=n, n_groups=g, ess_group=ess,
                outcome_sum=float(y.sum()), p_count=mean,
                neighbor_variance=variance, neighbor_sd=float(np.sqrt(variance)),
                donor_q025=float(np.quantile(y, .025)), donor_q975=float(np.quantile(y, .975)),
                mean_se_cluster=se, mean_low=float(low), mean_high=float(high),
                n_cases=int(y.sum()) if task == "binary" else np.nan,
                n_controls=int(n-y.sum()) if task == "binary" else np.nan,
                label_disagreement=float(2*mean*(1-mean)) if task == "binary" else np.nan)


def count_predict(donor_z, query_z, *, donor_ids, query_ids, donor_groups,
                  query_groups, donor_y, task, radius, min_count=20,
                  min_group_ess=15., device="cpu", query_batch=64,
                  donor_batch=2048, trace_sink=None, groups_verified=False):
    """Freeze all memberships before donor-outcome aggregation.

    trace_sink(start, ids, indptr, indices, distances) receives *all* in-radius
    memberships. Indices refer to the supplied donor bank.  All query identities
    survive, including empty and unsupported neighborhoods.
    """
    if not np.isfinite(radius) or radius <= 0 or not np.isfinite(radius**2):
        raise ValueError("radius must be finite and positive")
    if isinstance(min_count, bool) or not isinstance(min_count, (int, np.integer)) or min_count < 1:
        raise ValueError("min_count must be a positive integer")
    if not np.isfinite(min_group_ess) or min_group_ess <= 0:
        raise ValueError("min_group_ess must be positive")
    dz, qz, di, qi, dg, qg = _validate_inputs(
        donor_z, query_z, donor_ids, query_ids, donor_groups, query_groups,
        1, query_batch, donor_batch)
    y = validate_outcomes(donor_y, task)
    if len(y) != len(dz):
        raise ValueError("Donor outcomes and donor geometry have different row counts")
    if device != "cpu":
        raise ValueError("ABM uses CPU only")
    used = "cpu"
    prior = float(y.mean())
    canonical = np.argsort(di, kind="stable")
    rows = []
    last = time.monotonic()
    for start in range(0, len(qz), query_batch):
        q = np.array(qz[start:start+query_batch], dtype=np.float64, copy=True)
        nq = len(q)
        indices = [[] for _ in range(nq)]
        distances = [[] for _ in range(nq)]
        eligible = np.zeros(nq, dtype=np.int64)
        nearest = np.full(nq, np.inf)
        qnorm = np.einsum("ij,ij->i", q, q)
        for ds in range(0, len(dz), donor_batch):
            ii = canonical[ds:ds+donor_batch]
            d = np.asarray(dz[ii], dtype=np.float64)
            dnorm = np.einsum("ij,ij->i", d, d)
            dd = np.maximum(qnorm[:, None]+dnorm[None, :]-2*(q@d.T), 0)
            if not np.isfinite(dd).all():
                raise ValueError("Distance overflow in frozen genetic representation")
            # Both devices use the identical direct expression at the boundary.
            tolerance = 64*np.finfo(float).eps*(qnorm[:, None]+dnorm[None, :]+radius**2+1)
            er, ec = np.nonzero(np.abs(dd-radius**2) <= tolerance)
            if len(er):
                delta = q[er]-d[ec]
                dd[er, ec] = np.einsum("ij,ij->i", delta, delta)
            ok = ((qi[start:start+nq, None] != di[ii][None, :]) &
                  (qg[start:start+nq, None] != dg[ii][None, :]))
            eligible += ok.sum(1)
            nearest = np.minimum(nearest, np.where(ok, dd, np.inf).min(1))
            inside = ok & (dd <= radius**2)
            for j in range(nq):
                take = np.flatnonzero(inside[j])
                if len(take):
                    indices[j].append(ii[take])
                    distances[j].append(np.sqrt(dd[j, take]))
        flat_i, flat_d, offsets = [], [], [0]
        for j in range(nq):
            ii = np.concatenate(indices[j]) if indices[j] else np.empty(0, dtype=np.int64)
            dd = np.concatenate(distances[j]) if distances[j] else np.empty(0, dtype=float)
            stats = summarize_neighbors(y[ii], dg[ii], task)
            enough = stats["n_twins"] >= min_count and stats["ess_group"] >= min_group_ess-1e-10
            reason = ("no_eligible_donors" if not eligible[j] else "no_twins_in_radius" if not len(ii)
                      else "too_few_twins" if len(ii) < min_count
                      else "low_group_ess" if not enough else "supported")
            row = dict(id=qi[start+j], group_id=qg[start+j], **stats,
                       n_eligible=int(eligible[j]), d_nearest=float(np.sqrt(nearest[j])),
                       radius=float(radius), supported=bool(enough), support_reason=reason,
                       p_abm=stats["p_count"] if enough else np.nan,
                       p_population=stats["p_count"] if enough else prior,
                       p_global_prior=prior, fallback_used=not enough,
                       prediction_rule="all_in_radius_equal_weight_mean" if enough else "global_donor_mean_fallback",
                       groups_verified=bool(groups_verified))
            rows.append(row)
            flat_i.append(ii); flat_d.append(dd); offsets.append(offsets[-1]+len(ii))
        if trace_sink is not None:
            trace_sink(start, qi[start:start+nq], np.asarray(offsets, np.int64),
                       np.concatenate(flat_i), np.concatenate(flat_d))
        if time.monotonic()-last >= 45:
            print(f"[GRID ABM] {used}: {start+nq:,}/{len(qz):,} queries; {len(dz):,} donors; all in-radius outcomes counted", flush=True)
            last = time.monotonic()
    return CountResult(pd.DataFrame(rows), used)


def freeze_policy(calibration, coverages=(.2, .4, .5, .6, .8, 1.), heterogeneity_quantile=.9):
    """Use only development query geometry and training-donor outcomes.

    These are empirical rank thresholds, not error regressions, calibrated
    individual risks, or promised test-set coverages.
    """
    if "y" in calibration:
        raise ValueError("Selection calibration must not contain query outcomes")
    if not len(calibration):
        raise ValueError("A nonempty development calibration roster is required")
    scores = np.where(calibration.supported, calibration.ess_group, -np.inf)
    finite = calibration.loc[calibration.supported, "neighbor_variance"].dropna()
    dispersion = float(finite.quantile(heterogeneity_quantile)) if len(finite) else None
    thresholds = {}
    for coverage in coverages:
        if not 0 < coverage <= 1:
            raise ValueError("Coverage must be in (0,1]")
        # Inverted-CDF rank avoids interpolation with -infinity.
        target_n = max(1, int(np.ceil(coverage*len(scores)-1e-12)))
        ix = max(0, len(scores)-target_n)
        threshold = float(np.sort(scores)[ix])
        thresholds[str(float(coverage))] = threshold if np.isfinite(threshold) else None
    p = calibration.p_population.to_numpy(float)
    edges = sorted(set(float(x) for x in np.quantile(p, np.arange(.1, 1., .1))))
    return dict(score="ess_group", score_direction="larger_is_more_support",
                threshold_role="development_calibration_without_query_y",
                thresholds=thresholds, risk_bin_edges=edges,
                dispersion_threshold=dispersion, heterogeneity_quantile=heterogeneity_quantile,
                uncertainty_statement="descriptive_neighbor_heterogeneity_not_individual_error_guarantee")


def apply_policy(predictions, policy):
    if "y" in predictions:
        raise ValueError("Predictions and selection must be frozen before attaching query outcomes")
    out = predictions.copy()
    for coverage, threshold in policy["thresholds"].items():
        suffix = f"{float(coverage):.2f}".replace(".", "p")
        out["selected_"+suffix] = out.supported & (True if threshold is None else out.ess_group >= threshold)
    cut = policy["dispersion_threshold"]
    out["evidence_state"] = np.where(~out.supported, "low_support",
        np.where(False if cut is None else out.neighbor_variance > cut,
                 "supported_heterogeneous", "supported_lower_dispersion"))
    out["risk_bin"] = np.searchsorted(policy["risk_bin_edges"], out.p_population, side="right")
    return out
