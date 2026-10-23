#!/usr/bin/env python3
"""One frozen-roster evaluator for all scientific C1 methods.

No fit, threshold tuning or winner selection occurs here. All method contrasts
use exactly the same people, horizon, censoring estimator and family draws.
Selection changes the target population; it is therefore reported separately
from model effects. Aggregate tables contain no individual identifiers.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import abm

METRICS = ("Harrell_C", "AUC_IPCW", "Brier_IPCW", "LogLoss_IPCW", "AUPRC_IPCW")
COUNT_EDGES = np.array([0., 1., 5., 10., 20., 40., 64., 128., 256., 512., 1024., np.inf])
FALLBACK_RISK_EDGES = np.array([.01, .02, .05, .1, .2, .5])
ERROR_EDGES = np.array([0., .01, .02, .05, .1, .2, .4, .6, 1.00000001])
_INVALID_IDS = {"", "nan", "none", "na", "<na>"}


def _ids(values, label):
    raw = np.asarray(values)
    if raw.ndim != 1 or pd.isna(raw).any():
        raise ValueError(f"{label} must be a complete one-dimensional identifier array")
    out = raw.astype(str)
    if any(x.strip().lower() in _INVALID_IDS for x in out):
        raise ValueError(f"{label} contains an empty identifier")
    return out


def _vector(values, n, label, *, boolean=False):
    raw = np.asarray(values)
    if raw.shape != (n,):
        raise ValueError(f"{label} must provide exactly one value per participant")
    if boolean:
        if raw.dtype.kind not in "biuf" or not np.isfinite(raw).all() or not np.isin(raw, [0, 1]).all():
            raise ValueError(f"{label} must be a complete boolean mask")
        return raw.astype(bool)
    return np.asarray(raw, float)


def _align(predictions, truth):
    if not predictions:
        raise ValueError("At least one frozen prediction method is required")
    ids = _ids(truth["eid"], "Truth IDs")
    n = len(ids)
    if not n or len(set(ids)) != n:
        raise ValueError("Truth IDs must be nonempty and unique")
    time = _vector(truth["time"], n, "Truth time")
    event = _vector(truth["event"], n, "Truth event")
    if not np.isfinite(time).all() or (time <= 0).any() or not np.isin(event, [0, 1]).all():
        raise ValueError("Invalid survival outcomes")
    aligned, groups = {}, None
    for name, pred in predictions.items():
        if not isinstance(name, str) or not name:
            raise ValueError("Method names must be nonempty strings")
        pids = _ids(pred["eid"], f"{name} IDs")
        if len(pids) != n or len(set(pids)) != n or set(pids) != set(ids):
            raise ValueError(f"{name}: all methods must have exactly the same participant ID set")
        lookup = {v: j for j, v in enumerate(pids)}
        order = np.array([lookup[v] for v in ids], int)
        current_groups = _ids(pred["group"], f"{name} families")
        if len(current_groups) != n:
            raise ValueError(f"{name}: invalid family roster")
        current_groups = current_groups[order]
        if groups is not None and not np.array_equal(groups, current_groups):
            raise ValueError(f"{name}: family assignments disagree across methods")
        groups = current_groups
        probability = _vector(pred["probability"], n, f"{name} probability")[order]
        finite = np.isfinite(probability)
        if ((probability[finite] < 0) | (probability[finite] > 1)).any():
            raise ValueError(f"{name}: probability outside [0,1]")
        qc = _vector(pred["technical_qc"], n, f"{name} technical_qc", boolean=True)[order]
        supported = (None if pred.get("supported") is None else
                     _vector(pred["supported"], n, f"{name} supported", boolean=True)[order])
        masks = {}
        for q, mask in pred.get("masks", {}).items():
            q = float(q)
            if not np.isfinite(q) or not 0 < q <= 1 or q in masks:
                raise ValueError(f"{name}: invalid or duplicate coverage target")
            masks[q] = _vector(mask, n, f"{name} mask {q:g}", boolean=True)[order]
        diagnostics = {key: _vector(value, n, f"{name} diagnostics {key}")[order]
                       for key, value in pred.get("diagnostics", {}).items()}
        optional = {key: None if pred.get(key) is None else
                    _vector(pred[key], n, f"{name} {key}")[order]
                    for key in ("predicted_error", "predictability", "residual_predictability")}
        if optional["predicted_error"] is not None:
            good = np.isfinite(optional["predicted_error"])
            if ((optional["predicted_error"][good] < 0) |
                (optional["predicted_error"][good] > 1)).any():
                raise ValueError(f"{name}: predicted squared error must lie in [0,1]")
        aligned[name] = dict(probability=probability, technical_qc=qc, supported=supported,
            masks=masks, diagnostics=diagnostics, original=pred, **optional)
    master_qc = np.logical_and.reduce([p["technical_qc"] for p in aligned.values()])
    common_finite = np.logical_and.reduce([np.isfinite(p["probability"]) for p in aligned.values()])
    return ids, groups, time, event.astype(int), aligned, master_qc, common_finite


class _MetricPlan:
    """Cached rank plans, avoiding repeated O(N log N) sorting in bootstrap."""

    def __init__(self, time, event, probability, y, w, *, harrell=True):
        self.p = np.asarray(probability, float)
        self.y, self.w = np.asarray(y, int), np.asarray(w, float)
        self.error = (self.y - self.p) ** 2
        self.loss = abm.logloss_rows(self.y, self.p)
        # Descending scores with all ties collapsed reproduces average precision
        # (step-function PR area), rather than a trapezoidal PR interpolation.
        self.order = np.argsort(-self.p, kind="stable")
        ordered = self.p[self.order]
        self.starts = np.r_[0, np.flatnonzero(ordered[1:] != ordered[:-1]) + 1] if len(ordered) else np.array([], int)
        self.harrell = abm._harrell_bootstrap_plan(time, event, self.p) if harrell else None

    def evaluate(self, multiplicity):
        m = np.asarray(multiplicity, float)
        ww = self.w * m
        total = float(ww.sum())
        out = dict.fromkeys(METRICS, np.nan)
        out["Harrell_comparable_pairs"] = 0.
        if self.harrell is not None:
            out["Harrell_C"], out["Harrell_comparable_pairs"] = self.harrell(m)
        if total <= 0:
            return out
        out["Brier_IPCW"] = float(ww @ self.error / total)
        out["LogLoss_IPCW"] = float(ww @ self.loss / total)
        case, control = ww * self.y, ww * (1 - self.y)
        pos, neg = float(case.sum()), float(control.sum())
        if pos <= 0 or neg <= 0:
            return out
        a = np.add.reduceat(case[self.order], self.starts)
        b = np.add.reduceat(control[self.order], self.starts)
        # For each descending score group: lower scores are concordant controls.
        lower = neg - np.cumsum(b)
        out["AUC_IPCW"] = float(np.sum(a * (lower + .5 * b)) / (pos * neg))
        cumulative_cases = np.cumsum(a)
        cumulative_total = cumulative_cases + np.cumsum(b)
        precision = np.divide(cumulative_cases, cumulative_total,
                              out=np.zeros_like(cumulative_cases), where=cumulative_total > 0)
        out["AUPRC_IPCW"] = float(np.sum(a * precision) / pos)
        return out


def _point_metrics(time, event, probability, km, y, w):
    # Canonical score/weight normalization remains owned by the existing core.
    row = abm.metrics(time, event, probability, km, include_harrell=False, include_uno=False)
    plan = _MetricPlan(time, event, probability, y, w)
    row.update(plan.evaluate(np.ones(len(time))))
    total = float(np.sum(w))
    observed = float(np.sum(w * y) / total) if total else np.nan
    expected = float(np.sum(w * probability) / total) if total else np.nan
    row.update(observed_risk_IPCW=observed, expected_risk_IPCW=expected,
               observed_expected_ratio=(observed / expected if expected > 0 else np.nan),
               calibration_in_the_large_difference=(observed - expected),
               AUPRC_definition="IPCW weighted average precision; tied scores grouped",
               calibration_fit_on_test=False)
    return row


def _interval(point, draws, B, minimum, family_n):
    values = np.asarray(draws, float)
    values = values[np.isfinite(values)]
    if not np.isfinite(point):
        status = "insufficient_information"
    elif B == 0:
        status = "not_requested"
    elif family_n < 2:
        status = "insufficient_independent_groups"
    elif len(values) < minimum:
        status = "insufficient_bootstrap"
    else:
        status = "evaluated"
    lo, hi = (np.quantile(values, [.025, .975]) if status == "evaluated" else (np.nan, np.nan))
    return dict(lower=float(lo), upper=float(hi), bootstrap_valid=len(values), requested_bootstrap=B,
                minimum_valid_bootstrap=minimum, inference_status=status)


def _bootstrap(plans, masks, inverse, family_n, B, seed):
    """The same global family draw is used for every arm and every subset.

    Reinitializing the seeded generator per evaluation block reproduces the exact
    draws without retaining a B x participant matrix. Subsets never resample only
    their own selected families, and group sizes are preserved by multiplicities.
    """
    out = {key: {metric: [] for metric in METRICS} for key in plans}
    if B == 0 or family_n < 2:
        return out
    rng = np.random.default_rng(seed)
    family_maps = {key: inverse[masks[key]] for key in plans}
    for _ in range(B):
        counts = np.bincount(rng.integers(0, family_n, size=family_n), minlength=family_n)
        for key, plan in plans.items():
            values = plan.evaluate(counts[family_maps[key]])
            for metric in METRICS:
                out[key][metric].append(values[metric])
    return out


def _edges_from_baseline(baseline, method_name):
    original = baseline["original"]
    edges = original.get("risk_edges")
    if edges is None:
        return FALLBACK_RISK_EDGES.copy(), "predeclared_probability_bins"
    edges = np.asarray(edges, float)
    if edges.ndim != 1 or not np.isfinite(edges).all() or (np.diff(edges) <= 0).any():
        raise ValueError("Frozen baseline risk_edges must be finite and strictly increasing")
    if ((edges < 0) | (edges > 1)).any():
        raise ValueError("Frozen baseline risk_edges must be probability boundaries")
    source = str(original.get("risk_edges_source", "provided_frozen_development_edges"))
    lower_source = source.lower()
    if lower_source.startswith(("test", "evaluation", "outer_test")) or any(
            token in lower_source for token in ("test_quantile", "evaluation_quantile")):
        raise ValueError("Risk matching boundaries must not be fitted to evaluation data")
    return edges, f"{method_name}:{source}"


def _weighted_mean(values, weights):
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    return float(np.average(values[ok], weights=weights[ok])) if ok.any() else np.nan


def _support_tables(aligned, baseline_name, common, groups, y, w, risk_labels, risk_source):
    distribution, errors, status = [], [], []
    base = aligned[baseline_name]["probability"]
    for name, pred in aligned.items():
        diag = pred["diagnostics"]
        if "eligible_match_count" not in diag:
            status.append(dict(model=name, component="matched_twin_count", status="unavailable",
                               reason="This method did not provide a frozen neighbourhood diagnostic"))
            continue
        count = diag["eligible_match_count"]
        good = common & np.isfinite(count)
        if (count[good] < 0).any() or not np.allclose(count[good], np.round(count[good])):
            raise ValueError(f"{name}: eligible_match_count must be a nonnegative integer count")
        used = diag.get("matched_count", np.full(len(common), np.nan))
        if (used[good & np.isfinite(used)] < 0).any() or (used[good & np.isfinite(used)] > count[good & np.isfinite(used)]).any():
            raise ValueError(f"{name}: used twins cannot exceed eligible in-radius twins")
        row = dict(model=name, N=int(good.sum()), common_N=int(common.sum()),
                   missing_count_N=int(np.sum(common & ~np.isfinite(count))),
                   count_definition="all eligible incident in-radius donors before K truncation",
                   count_comparison_scope="within a frozen geometry, reference population and radius",
                   historical_count_note="prevalent source support is separate and never added to incident twin count",
                   accuracy_guarantee=False)
        for key in ("eligible_match_count", "matched_count", "eligible_reference_count", "twin_fraction",
                    "family_ess", "matching_radius", "radius_ratio",
                    "prevalent_eligible_reference_count", "prevalent_matched_count", "prevalent_used_count",
                    "prevalent_family_ess", "prevalent_matching_radius", "prevalent_support_available"):
            values = np.asarray(diag.get(key, np.full(len(common), np.nan)))[good]
            values = values[np.isfinite(values)]
            row[key + "_mean"] = float(values.mean()) if len(values) else np.nan
            row[key + "_median"] = float(np.median(values)) if len(values) else np.nan
            row[key + "_min"] = float(values.min()) if len(values) else np.nan
            row[key + "_max"] = float(values.max()) if len(values) else np.nan
        row["zero_twin_fraction"] = float(np.mean(count[good] == 0)) if good.any() else np.nan
        cap_ok = good & np.isfinite(used)
        row["fraction_with_more_eligible_than_used"] = float(np.mean(count[cap_ok] > used[cap_ok])) if cap_ok.any() else np.nan
        distribution.append(row)
        p = pred["probability"]
        for lo, hi in zip(COUNT_EDGES[:-1], COUNT_EDGES[1:]):
            count_mask = good & (count >= lo) & (count < hi)
            for risk_bin in ["all"] + list(np.unique(risk_labels[common])):
                mask = count_mask if risk_bin == "all" else count_mask & (risk_labels == risk_bin)
                if not mask.any():
                    continue
                errors.append(dict(model=name, count_lower=int(lo), count_upper=(float(hi) if np.isfinite(hi) else np.inf),
                    risk_bin=str(risk_bin), risk_bin_source=risk_source, N=int(mask.sum()),
                    known_N=int(np.sum(mask & (w > 0))), cases=int(np.sum(y[mask])),
                    family_N=len(np.unique(groups[mask])),
                    mean_raw_count=float(np.mean(count[mask])),
                    Brier_IPCW=_weighted_mean((y[mask] - p[mask]) ** 2, w[mask]),
                    LogLoss_IPCW=_weighted_mean(abm.logloss_rows(y[mask], p[mask]), w[mask]),
                    baseline_Brier_IPCW=_weighted_mean((y[mask] - base[mask]) ** 2, w[mask]),
                    mean_baseline_risk=_weighted_mean(base[mask], w[mask]),
                    mean_predicted_error=(_weighted_mean(pred["predicted_error"][mask], w[mask])
                                          if pred["predicted_error"] is not None else np.nan),
                    **abm.ipcw_information(y[mask], w[mask]),
                    interpretation="descriptive error within fixed support and baseline-risk strata; no accuracy guarantee"))
        status.append(dict(model=name, component="matched_twin_count", status="evaluated",
                           reason="Raw count, used count, reference size and risk-stratified errors reported separately"))
    return distribution, errors, status


def _predictability_tables(aligned, common, y, w, risk_labels, risk_source):
    rows = []
    for name, pred in aligned.items():
        estimated = pred["predicted_error"]
        if estimated is None:
            continue
        good = common & np.isfinite(estimated)
        actual = (y - pred["probability"]) ** 2
        for risk_bin in ["all"] + list(np.unique(risk_labels[common])):
            scope = good if risk_bin == "all" else good & (risk_labels == risk_bin)
            # Includes a full-stratum row plus fixed expected squared-error bins.
            bins = [(np.nan, np.nan, scope)] + [(lo, hi, scope & (estimated >= lo) & (estimated < hi))
                      for lo, hi in zip(ERROR_EDGES[:-1], ERROR_EDGES[1:])]
            for lo, hi, mask in bins:
                if not mask.any():
                    continue
                expected = _weighted_mean(estimated[mask], w[mask])
                observed = _weighted_mean(actual[mask], w[mask])
                row = dict(model=name, risk_bin=str(risk_bin), risk_bin_source=risk_source,
                           predicted_error_lower=lo, predicted_error_upper=hi,
                           N=int(mask.sum()), known_N=int(np.sum(mask & (w > 0))),
                           cases=int(np.sum(y[mask])), mean_predicted_error=expected,
                           observed_Brier_IPCW=observed, observed_minus_expected_error=observed - expected,
                           interpretation="Expected squared-error calibration; not individual certainty",
                           test_calibration_fit=False)
                for key in ("predictability", "residual_predictability"):
                    score = pred[key]
                    row[key + "_mean"] = _weighted_mean(score[mask], w[mask]) if score is not None else np.nan
                    known = mask & (w > 0) & (np.isfinite(score) if score is not None else False)
                    usable = score is not None and known.sum() >= 3 and len(np.unique(score[known])) > 1 and len(np.unique(actual[known])) > 1
                    row[key + "_error_spearman"] = float(spearmanr(score[known], actual[known]).statistic) if usable else np.nan
                rows.append(row)
    return rows


def evaluate_methods(predictions: dict[str, dict[str, Any]], truth: dict[str, Any], km, *,
                     bootstrap=1000, seed=2026, grouping_unit="family", primary_coverage=.5):
    """Return aggregate DataFrames from frozen predictions on one shared roster.

    All methods are evaluated on every selector's exactly identical subset.
    Inference is predeclared for all-cohort model contrasts and each selector's
    own method versus the baseline on primary coverage. Other coverage points
    are descriptive. Risk-matched selection controls are outcome-blind, use fixed
    development risk boundaries, and remain a batch-level composition control.
    """
    if isinstance(bootstrap, bool) or int(bootstrap) != bootstrap or bootstrap < 0:
        raise ValueError("bootstrap must be a nonnegative integer")
    if grouping_unit not in {"family", "individual_unverified"}:
        raise ValueError("Unknown bootstrap grouping_unit")
    if not 0 < float(primary_coverage) <= 1:
        raise ValueError("primary_coverage must be in (0,1]")
    B, minimum = int(bootstrap), max(20, int(math.ceil(.8 * bootstrap)))
    ids, groups, time, event, aligned, master_qc, finite = _align(predictions, truth)
    common = master_qc & finite
    names = list(aligned)
    baseline_name = "glmnet" if "glmnet" in aligned else names[0]
    history_pair = {"abm_prevalent", "abm_prevalent_control"} <= set(names)
    y, w = km.labels_weights(time, event)
    risk_edges, risk_source = _edges_from_baseline(aligned[baseline_name], baseline_name)
    risk = aligned[baseline_name]["probability"]
    risk_labels = np.searchsorted(risk_edges, risk, side="right")
    families = np.unique(groups[common])
    family_index = {g: j for j, g in enumerate(families)}
    inverse = np.array([family_index.get(g, -1) for g in groups], int)
    family_n = len(families)
    subsets = {"all_common": common}
    requests = {"all_common": master_qc}
    subset_metadata = {"all_common": dict(selector="none", coverage_target=1., subset_kind="all_cohort")}
    controls, coverage = {}, []
    for name, pred in aligned.items():
        for q, frozen in sorted(pred["masks"].items()):
            key = f"selected::{name}::{q:g}"
            requested = frozen & master_qc
            selected = requested & finite
            subsets[key], requests[key] = selected, requested
            subset_metadata[key] = dict(selector=name, coverage_target=q, subset_kind="selected")
            control = abm.risk_matched_mask(selected, common, risk, ids, risk_edges, seed + 412)
            if control.sum() != selected.sum():
                raise AssertionError("Risk-matched control changed the selected sample size")
            controls[key] = control
            supported = pred["supported"]
            coverage.append(dict(selector=name, coverage_target=q, original_N=len(ids),
                common_N=int(common.sum()), requested_selected_N=int(requested.sum()), selected_N=int(selected.sum()),
                actual_coverage_all=float(selected.mean()),
                actual_coverage_common=float(selected.sum() / common.sum()) if common.any() else np.nan,
                cases=int(np.sum(y[selected])), controls=int(np.sum(selected & (w > 0) & (y == 0))),
                case_coverage=float(np.sum(y[selected]) / np.sum(y[common])) if np.sum(y[common]) else np.nan,
                control_coverage=float(np.sum(selected & (w > 0) & (y == 0)) / np.sum(common & (w > 0) & (y == 0)))
                    if np.sum(common & (w > 0) & (y == 0)) else np.nan,
                family_N=len(np.unique(groups[selected])), riskmatched_N=int(control.sum()),
                selected_riskmatched_overlap_N=int(np.sum(selected & control)),
                supported_fraction=(float(np.mean(supported[common])) if supported is not None and common.any() else np.nan),
                requested_roster_hash=abm.fingerprint(ids[requested]), evaluation_roster_hash=abm.fingerprint(ids[selected]),
                riskmatched_roster_hash=abm.fingerprint(ids[control]),
                risk_matching_edges_source=risk_source, threshold_source="frozen_development_gate",
                **abm.ipcw_information(y[selected], w[selected])))
    performance, lookup = [], {}
    for subset, mask in subsets.items():
        metadata = subset_metadata[subset]
        for name, pred in aligned.items():
            row = _point_metrics(time[mask], event[mask], pred["probability"][mask], km, y[mask], w[mask])
            row.update(model=name, subset=subset, **metadata, horizon=float(km.horizon),
                requested_N=int(requests[subset].sum()), evaluated_N=int(mask.sum()),
                excluded_nonfinite_N=int(np.sum(requests[subset] & ~finite)),
                master_QC_excluded_N=int(np.sum(~master_qc)), original_N=len(ids),
                family_N=len(np.unique(groups[mask])), coverage_all=float(mask.mean()),
                coverage_common=float(mask.sum() / common.sum()) if common.any() else np.nan,
                evaluation_roster_hash=abm.fingerprint(ids[mask]),
                requested_roster_hash=abm.fingerprint(ids[requests[subset]]),
                same_people_within_subset=True, bootstrap_unit=grouping_unit)
            for metric in METRICS:
                row[metric + "_lower"] = row[metric + "_upper"] = np.nan
                row[metric + "_bootstrap_valid"] = 0
                row[metric + "_inference_status"] = "not_predeclared"
            performance.append(row)
            lookup[(subset, name)] = row
    contrasts, riskmatched = [], []

    def add_contrast(name, reference, subset, metric, point, draws, *, kind="model_effect_same_people",
                     reference_subset=None, reference_mask=None):
        mask = subsets[subset]
        refmask = mask if reference_mask is None else reference_mask
        contrasts.append(dict(model=name, reference=reference, subset=subset,
            reference_subset=subset if reference_subset is None else reference_subset,
            selector=subset_metadata[subset]["selector"], coverage_target=subset_metadata[subset]["coverage_target"],
            contrast=kind, metric=metric, difference=point, difference_definition="model minus reference",
            model_N=int(mask.sum()), reference_N=int(refmask.sum()),
            same_people=bool(np.array_equal(mask, refmask)),
            model_roster_hash=abm.fingerprint(ids[mask]), reference_roster_hash=abm.fingerprint(ids[refmask]),
            bootstrap_unit=grouping_unit, paired_family_draws=True,
            uncertainty="conditional on frozen models, gates, risk controls and build censoring estimator",
            **_interval(point, draws, B, minimum, family_n)))

    # Inference blocks share the exact full-population family sequence by seed.
    # Cache only one block's score-rank plans, not every coverage point at once.
    for subset, mask in subsets.items():
        meta = subset_metadata[subset]
        primary = subset == "all_common" or np.isclose(meta["coverage_target"], primary_coverage, rtol=0, atol=1e-12)
        if not primary:
            continue
        needed = names if subset == "all_common" else list(dict.fromkeys([baseline_name, meta["selector"]]))
        history_block = history_pair and (subset == "all_common" or meta["selector"] in
                                          {"abm_prevalent", "abm_prevalent_control"})
        if history_block:
            needed = list(dict.fromkeys([*needed, "abm_prevalent", "abm_prevalent_control"]))
        plans, plan_masks = {}, {}
        for name in needed:
            key = ("selected", name)
            plans[key] = _MetricPlan(time[mask], event[mask], aligned[name]["probability"][mask], y[mask], w[mask])
            plan_masks[key] = mask
            if subset in controls:
                control = controls[subset]
                key = ("riskmatched", name)
                plans[key] = _MetricPlan(time[control], event[control], aligned[name]["probability"][control], y[control], w[control])
                plan_masks[key] = control
        draws = _bootstrap(plans, plan_masks, inverse, family_n, B, seed)
        for name in needed:
            row = lookup[(subset, name)]
            for metric in METRICS:
                interval = _interval(row[metric], draws[("selected", name)][metric], B, minimum, family_n)
                for col in ("lower", "upper", "bootstrap_valid", "inference_status"):
                    row[metric + "_" + col] = interval[col]
                if name != baseline_name:
                    delta = row[metric] - lookup[(subset, baseline_name)][metric]
                    paired = (np.asarray(draws[("selected", name)][metric]) -
                              np.asarray(draws[("selected", baseline_name)][metric]))
                    add_contrast(name, baseline_name, subset, metric, delta, paired)
                if subset in controls:
                    point_control = plans[("riskmatched", name)].evaluate(np.ones(int(controls[subset].sum())))
                    delta = row[metric] - point_control[metric]
                    paired = (np.asarray(draws[("selected", name)][metric]) -
                              np.asarray(draws[("riskmatched", name)][metric]))
                    add_contrast(name, name, subset, metric, delta, paired,
                        kind="selection_vs_riskmatched_composition", reference_subset="riskmatched::" + subset,
                        reference_mask=controls[subset])
        if history_block:
            historical, nohistory = "abm_prevalent", "abm_prevalent_control"
            for metric in METRICS:
                delta = lookup[(subset, historical)][metric] - lookup[(subset, nohistory)][metric]
                paired = (np.asarray(draws[("selected", historical)][metric]) -
                          np.asarray(draws[("selected", nohistory)][metric]))
                add_contrast(historical, nohistory, subset, metric, delta, paired,
                             kind="historical_increment_same_architecture")
    # All requested coverage points retain descriptive selected/control scores.
    for subset, control in controls.items():
        selector = subset_metadata[subset]["selector"]
        for name in list(dict.fromkeys([baseline_name, selector])):
            selected_row = lookup[(subset, name)]
            control_row = _point_metrics(time[control], event[control], aligned[name]["probability"][control], km,
                                         y[control], w[control])
            row = dict(selector=selector, model=name, coverage_target=subset_metadata[subset]["coverage_target"],
                selected_N=selected_row["N"], riskmatched_N=control_row["N"],
                selected_cases=selected_row["cases"], riskmatched_cases=control_row["cases"],
                selected_roster_hash=selected_row["evaluation_roster_hash"],
                riskmatched_roster_hash=abm.fingerprint(ids[control]),
                same_people=bool(np.array_equal(subsets[subset], control)),
                risk_matching_edges_source=risk_source,
                interpretation="same-size fixed-risk-bin composition control; membership may depend on prediction batch")
            for metric in METRICS:
                row["selected_" + metric], row["riskmatched_" + metric] = selected_row[metric], control_row[metric]
                row["difference_" + metric] = selected_row[metric] - control_row[metric]
            riskmatched.append(row)
    roster = []
    for name, pred in aligned.items():
        roster.append(dict(model=name, original_N=len(ids), model_QC_pass_N=int(pred["technical_qc"].sum()),
            model_nonfinite_N=int(np.sum(~np.isfinite(pred["probability"]))),
            master_QC_pass_N=int(master_qc.sum()), excluded_common_nonfinite_N=int(np.sum(master_qc & ~finite)),
            common_N=int(common.sum()), original_roster_hash=abm.fingerprint(ids),
            common_roster_hash=abm.fingerprint(ids[common]), shared_truth=True,
            shared_family_assignments=True, shared_censoring_estimator=True,
            horizon=float(km.horizon), grouping_unit=grouping_unit))
    distribution, errors, component_status = _support_tables(aligned, baseline_name, common, groups, y, w, risk_labels, risk_source)
    calibration = _predictability_tables(aligned, common, y, w, risk_labels, risk_source)
    status = [dict(component="unified_evaluation", status="evaluated" if common.any() else "empty_common_population",
        baseline=baseline_name, common_N=int(common.sum()), independent_groups=family_n,
        grouping_unit=grouping_unit, independent_family_support_verified=grouping_unit == "family",
        bootstrap_requested=B, minimum_valid_bootstrap=minimum,
        primary_coverage=primary_coverage, risk_matching_edges_source=risk_source,
        outcome_model_or_gate_fit_on_test=False, winner_selected_on_test=False,
        note="Fixed-horizon IPCW net risk; full-followup strict-time Harrell C; internal conditional uncertainty")]
    status.extend(component_status)
    for name, pred in aligned.items():
        error = pred["predicted_error"]
        available = common & np.isfinite(error) if error is not None else np.zeros(len(common), bool)
        status.append(dict(model=name, component="individual_expected_error",
            status="evaluated" if available.any() else "unavailable",
            common_N=int(common.sum()), finite_predicted_error_N=int(available.sum()),
            unavailable_predicted_error_N=int(np.sum(common & ~available)),
            assessment_status=str(pred["original"].get("assessment_status", "unspecified")),
            reason="Only finite frozen expected errors can be calibrated; support-only ranking is not an accuracy estimate"))
    if "abm_prevalent" in names:
        status.append(dict(model="abm_prevalent", component="historical_increment",
            status="evaluated" if history_pair else "matched_control_unavailable",
            reference="abm_prevalent_control", reason="Same architecture and incident rows isolate incremental historical information"))
    return {key: pd.DataFrame(rows) for key, rows in dict(performance=performance, contrasts=contrasts,
        coverage=coverage, roster_audit=roster, support_distribution=distribution,
        support_error_by_count=errors, risk_matched_selection=riskmatched,
        predictability_calibration=calibration, evaluation_status=status).items()}
