#!/usr/bin/env python3
"""C1 common model evaluation, selector inference and P0/P1/P2 reliability."""
from __future__ import annotations


# Common held-out model evaluation
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


# Frozen selector inference
import hashlib, numpy as np, pandas as pd
from c1_protocol import COVERAGES, PRIMARY_SELECTOR, eid_hash

ROUTES=("risk_resid_error","logloss_resid","local_stability","model_consensus","hybrid_predictability",
        "support_quality","absolute_error","rank_error","support_only","lowrisk","confidence","random")

def _stable_hash(value,salt):
    return int.from_bytes(hashlib.sha256(f"{salt}:{value}".encode()).digest()[:8],"big")

def _metrics(time,event,score,km):
    time=np.asarray(time,float);event=np.asarray(event,int);score=np.asarray(score,float)
    y,w=km.labels_weights(time,event);ok=(w>0)&np.isfinite(score)
    out=dict(N=int(ok.sum()),cases=int(y[ok].sum()),Brier=np.nan,LogLoss=np.nan,AUC=np.nan)
    if not ok.any(): return out
    p=np.clip(score[ok],1e-7,1-1e-7)
    out["Brier"]=float(np.average((y[ok]-p)**2,weights=w[ok]))
    out["LogLoss"]=float(np.average(-(y[ok]*np.log(p)+(1-y[ok])*np.log(1-p)),weights=w[ok]))
    if len(np.unique(y[ok]))==2:
        from sklearn.metrics import roc_auc_score
        out["AUC"]=float(roc_auc_score(y[ok],p,sample_weight=w[ok]))
    return out

def _risk_match(selected,qc,risk,ids,salt,bins=10):
    selected=np.asarray(selected,bool);qc=np.asarray(qc,bool);risk=np.asarray(risk,float);ids=np.asarray(ids,str)
    ok=qc&np.isfinite(risk);out=np.zeros(len(ids),bool)
    selected=selected&ok
    if not ok.any():return out
    q=np.unique(np.quantile(risk[ok],np.linspace(0,1,bins+1)))
    if len(q)<3: labels=np.zeros(len(ids),int)
    else:
        q[0],q[-1]=-np.inf,np.inf;labels=np.searchsorted(q[1:-1],risk,side="right")
    for b in np.unique(labels[ok]):
        pool=np.flatnonzero(ok&(labels==b));n=int(np.sum(selected&(labels==b)))
        order=sorted(pool,key=lambda j:_stable_hash(ids[j],salt))
        out[np.asarray(order[:n],int)]=True
    return out

def _aurc(cov,val):
    x=np.asarray(cov,float);y=np.asarray(val,float);ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<2 or not ok.all():return np.nan
    x,y=x[ok],y[ok];o=np.argsort(x);x,y=x[o],y[o]
    return float(np.trapezoid(y,x)/(x[-1]-x[0])) if x[-1]>x[0] else np.nan

def _finite_mean(values):
    values=np.asarray(values,float)
    finite=values[np.isfinite(values)]
    return float(finite.mean()) if len(finite) else np.nan

def infer_selector(prediction,time,event,groups,km,coverages=(.2,.4,.5,.6,.8,1.0),routes=ROUTES,B=1000,seed=2026,
                   grouping_unit="individual_unverified"):
    if grouping_unit not in {"family", "individual_unverified"}:
        raise ValueError("Invalid grouping unit")
    if tuple(coverages)!=COVERAGES or PRIMARY_SELECTOR not in routes:
        raise ValueError("C1 LOCKPAPER requires the locked coverage grid and primary selector")
    if int(B)!=B or B<1: raise ValueError("Bootstrap replicates must be positive")
    if pd.isna(np.asarray(groups)).any(): raise ValueError("Missing family/group identifiers")
    ids=np.asarray(prediction["eid"],str);groups=np.asarray(groups,str)
    time=np.asarray(time,float);event=np.asarray(event,int)
    qc=np.asarray(prediction["diagnostics"]["technical_qc"],bool)
    risk=np.asarray(prediction["probability"]["elasticnet_full"],float)
    if not len(ids) or len(set(ids))!=len(ids) or any(len(x)!=len(ids) for x in (time,event,groups,qc,risk)):
        raise ValueError("Invalid selector evaluation roster")
    if pd.isna(groups).any() or np.any(np.char.strip(groups)==""):
        raise ValueError("Missing family/group identifiers")
    qc=qc&np.isfinite(risk)
    masks={k:np.asarray(v,bool)&qc for k,v in prediction.get("masks",{}).items()}
    if any(f"{PRIMARY_SELECTOR}:{c:g}" not in masks for c in coverages):
        raise ValueError("Frozen primary selector must provide all six coverage masks")
    rows=[];controls={}
    for route in routes:
        for cov in coverages:
            key=f"{route}:{cov:g}"
            if key not in masks:continue
            sel=masks[key];rm=_risk_match(sel,qc,risk,ids,f"lockpaper:{route}:{cov:g}")
            controls[(route,float(cov))]=rm
            a=_metrics(np.asarray(time)[sel],np.asarray(event)[sel],risk[sel],km)
            b=_metrics(np.asarray(time)[rm],np.asarray(event)[rm],risk[rm],km)
            rows.append(dict(route=route,coverage_target=float(cov),coverage=float(sel.mean()),N=a["N"],cases=a["cases"],
              selected_N=int(sel.sum()),riskmatched_total_N=int(rm.sum()),
              selected_eid_hash=eid_hash(ids[sel]),riskmatched_eid_hash=eid_hash(ids[rm]),
              riskmatched_N=b["N"],riskmatched_cases=b["cases"],Brier=a["Brier"],riskmatched_Brier=b["Brier"],
              delta_Brier=a["Brier"]-b["Brier"],LogLoss=a["LogLoss"],riskmatched_LogLoss=b["LogLoss"],
              delta_LogLoss=a["LogLoss"]-b["LogLoss"],AUC=a["AUC"],riskmatched_AUC=b["AUC"],delta_AUC=a["AUC"]-b["AUC"]))
    curve=pd.DataFrame(rows)
    point=pd.DataFrame([dict(route=r,AURC_delta_Brier=_aurc(g.coverage_target,g.delta_Brier),
                             AURC_delta_LogLoss=_aurc(g.coverage_target,g.delta_LogLoss),
                             mean_delta_AUC=_finite_mean(g.delta_AUC))
                        for r,g in curve.groupby("route",sort=False)])
    u,inv=np.unique(groups,return_inverse=True);order=np.argsort(inv,kind="stable")
    pos=np.split(order,np.cumsum(np.bincount(inv))[:-1]);rng=np.random.default_rng(seed+9181)
    boots={r:[] for r in point.route}
    for _ in range(int(B)):
        ix=np.concatenate([pos[j] for j in rng.integers(0,len(u),len(u))])
        for route in boots:
            rr=[]
            for cov in coverages:
                key=f"{route}:{cov:g}"
                if key not in masks or (route,float(cov)) not in controls:continue
                sm=masks[key][ix];cm=controls[(route,float(cov))][ix]
                a=_metrics(np.asarray(time)[ix][sm],np.asarray(event)[ix][sm],risk[ix][sm],km)
                b=_metrics(np.asarray(time)[ix][cm],np.asarray(event)[ix][cm],risk[ix][cm],km)
                rr.append((float(cov),a["Brier"]-b["Brier"],a["LogLoss"]-b["LogLoss"],a["AUC"]-b["AUC"]))
            if len(rr)>=2:
                z=np.asarray(rr,float);boots[route].append((_aurc(z[:,0],z[:,1]),_aurc(z[:,0],z[:,2]),_finite_mean(z[:,3])))
    infer=[]
    minimum_valid=max(20,int(np.ceil(.8*B)))
    for _,p in point.iterrows():
        arr=np.asarray(boots[p.route],float)
        def ci(j):
            v=arr[:,j] if arr.ndim==2 else np.array([]);v=v[np.isfinite(v)]
            return (float(np.quantile(v,.025)),float(np.quantile(v,.975)),len(v)) if len(v)>=minimum_valid else (np.nan,np.nan,len(v))
        c0,c1,c2=ci(0),ci(1),ci(2)
        def inference_status(value,interval):
            return ("insufficient_information" if not np.isfinite(value) else
                    "evaluated" if interval[2]>=minimum_valid else "insufficient_bootstrap")
        s0=inference_status(p.AURC_delta_Brier,c0)
        s1=inference_status(p.AURC_delta_LogLoss,c1)
        s2=inference_status(p.mean_delta_AUC,c2)
        infer.append(dict(route=p.route,AURC_delta_Brier=p.AURC_delta_Brier,AURC_delta_Brier_lo=c0[0],AURC_delta_Brier_hi=c0[1],
                          AURC_delta_LogLoss=p.AURC_delta_LogLoss,AURC_delta_LogLoss_lo=c1[0],AURC_delta_LogLoss_hi=c1[1],
                          mean_delta_AUC=p.mean_delta_AUC,mean_delta_AUC_lo=c2[0],mean_delta_AUC_hi=c2[1],
                          valid_boot=min(c0[2],c1[2],c2[2]),requested_boot=int(B),minimum_valid_boot=minimum_valid,
                          Brier_bootstrap_valid=c0[2],LogLoss_bootstrap_valid=c1[2],AUC_bootstrap_valid=c2[2],
                          Brier_inference_status=s0,LogLoss_inference_status=s1,AUC_inference_status=s2,
                          status=s0 if s0!="evaluated" else s1 if s1!="evaluated" else s2,
                          bootstrap_unit=grouping_unit,
                          independent_family_support_verified=grouping_unit == "family",
                          uncertainty="conditional on frozen models, selector, risk-matched controls and censoring estimator"))
    return curve,pd.DataFrame(infer)


# P0/P1/P2 reliability evaluation
import argparse
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import os
import tempfile
import fcntl
import sys
from typing import Any

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import roc_auc_score, average_precision_score

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import abm  # noqa: E402

from c1_protocol import RELIABILITY_CONTRACT_VERSION
VERSION = RELIABILITY_CONTRACT_VERSION
PROTOCOL_PATH = HERE.parent/"config/c1_reliability_protocol.json"
DEFAULT_COVERAGES = (0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00)
PRIMARY_COVERAGE = 0.50
RISK_BINS = 50
MIN_EVENTS = 20
MIN_CONTROLS = 20


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")


def weighted_mean(x, w):
    x, w = np.asarray(x, float), np.asarray(w, float)
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not ok.any():
        return np.nan
    return float(np.sum(x[ok] * w[ok]) / np.sum(w[ok]))


def logloss_rows(y, p):
    p = np.clip(np.asarray(p, float), 1e-8, 1 - 1e-8)
    y = np.asarray(y, float)
    return -(y * np.log(p) + (1 - y) * np.log(1 - p))


def safe_auc(y, p, w):
    ok = np.isfinite(y) & np.isfinite(p) & np.isfinite(w) & (w > 0)
    if ok.sum() < 3 or len(np.unique(y[ok])) < 2:
        return np.nan
    return float(roc_auc_score(y[ok], p[ok], sample_weight=w[ok]))


def safe_auprc(y, p, w):
    ok = np.isfinite(y) & np.isfinite(p) & np.isfinite(w) & (w > 0)
    if ok.sum() < 3 or len(np.unique(y[ok])) < 2:
        return np.nan
    return float(average_precision_score(y[ok], p[ok], sample_weight=w[ok]))


def point_metrics(y, w, p, mask):
    mask = np.asarray(mask, bool) & np.isfinite(p) & (w >= 0)
    known = mask & (w > 0)
    if not known.any():
        return dict(N=int(mask.sum()), known_N=0, cases=0, Brier=np.nan, LogLoss=np.nan,
                    AUC=np.nan, AUPRC=np.nan, observed_risk=np.nan, expected_risk=np.nan)
    brier = weighted_mean((y - p) ** 2, w * known)
    ll = weighted_mean(logloss_rows(y, p), w * known)
    return dict(N=int(mask.sum()), known_N=int(known.sum()), cases=int(np.sum(y[known] == 1)),
                Brier=brier, LogLoss=ll, AUC=safe_auc(y[known], p[known], w[known]),
                AUPRC=safe_auprc(y[known], p[known], w[known]),
                observed_risk=weighted_mean(y, w * known), expected_risk=weighted_mean(p, w * known))


def as_1d(d: dict[str, Any], key: str, n: int, default=0.0):
    v = d.get(key)
    if v is None:
        return np.full(n, default, float)
    v = np.asarray(v)
    if v.ndim != 1 or len(v) != n:
        raise ValueError(f"Diagnostic {key} is not an aligned vector")
    return v.astype(float)


def align_prediction(reference_ids, pred: dict[str, Any]) -> dict[str, Any]:
    ids = np.asarray(pred["eid"]).astype(str)
    ref = np.asarray(reference_ids).astype(str)
    if len(np.unique(ids)) != len(ids) or len(np.unique(ref)) != len(ref):
        raise ValueError("Duplicate identifiers in prediction roster")
    if len(ids) == len(ref) and np.array_equal(ids, ref):
        return pred
    pos = {v: i for i, v in enumerate(ids)}
    if len(pos) != len(ids) or any(v not in pos for v in ref):
        raise ValueError("Prediction roster does not match the requested roster")
    ix = np.array([pos[v] for v in ref], int)
    out = dict(pred)
    for k in ("eid", "group", "probability", "technical_qc", "supported", "predicted_error",
              "predictability", "residual_predictability"):
        if k in out and out[k] is not None:
            v = np.asarray(out[k])
            if v.ndim == 1 and len(v) == len(ids):
                out[k] = v[ix]
    if "diagnostics" in out:
        out["diagnostics"] = {k: (np.asarray(v)[ix] if np.asarray(v).ndim == 1 and len(v) == len(ids) else v)
                              for k, v in out["diagnostics"].items()}
    return out


def completed_method_dirs(scope: Path, method: str):
    root = scope / "methods" / method
    if not root.is_dir():
        return []
    candidates = []
    for d in root.iterdir():
        if not d.is_dir():
            continue
        if (d / "COMPLETE.json").is_file() and (d / "test_predictions_LOCAL_ONLY.joblib").is_file()\
                and (d / "calibration_audit_predictions_LOCAL_ONLY.joblib").is_file():
            candidates.append(d)
    return sorted(candidates, key=lambda p: (p.stat().st_mtime_ns, p.name), reverse=True)


def locate_method(scope: Path, method: str) -> Path | None:
    # Prefer the exact method path recorded by the latest completed unified
    # evaluation.  This avoids accidentally mixing historical signatures when
    # several completed method directories coexist.
    eval_root = scope / "evaluation"
    contracts = sorted(eval_root.glob("*/evaluation_contract.json"),
                       key=lambda p: (p.stat().st_mtime_ns, str(p)), reverse=True) if eval_root.is_dir() else []
    if (scope / "C1_FINAL_INDEX.json").is_file():
        contracts.insert(0, scope / "C1_FINAL_INDEX.json")
    for contract_path in contracts:
        try:
            obj = json_read(contract_path)
            value = obj.get("methods", {}).get(method, {}).get("path")
            if value:
                d = Path(value)
                if not d.is_absolute():
                    d = (scope / d).resolve()
                if (d / "COMPLETE.json").is_file() and (d / "test_predictions_LOCAL_ONLY.joblib").is_file() and\
                   (d / "calibration_audit_predictions_LOCAL_ONLY.joblib").is_file():
                    return d
        except Exception:
            pass
    dirs = completed_method_dirs(scope, method)
    return dirs[0] if dirs else None


def locate_prepared_source(scope: Path, test_ids) -> Path:
    root = scope / "_LOCAL_ONLY" / "inputs"
    if not root.is_dir():
        raise FileNotFoundError(f"Prepared C1 input cache not found: {root}")
    target_hash = abm.fingerprint(np.asarray(test_ids).astype(str))
    matches = []
    for d in root.iterdir():
        contract_path = d / "host_contract.json"
        dev = d / "development_LOCAL_ONLY.npz"
        truth = d / "test_truth_LOCAL_ONLY.npz"
        if not (d.is_dir() and contract_path.is_file() and dev.is_file() and truth.is_file()):
            continue
        try:
            contract = json_read(contract_path)
        except Exception:
            continue
        if contract.get("incident_test_ids_hash") == target_hash:
            matches.append(d)
    if len(matches) != 1:
        raise RuntimeError(f"Expected exactly one prepared input bundle for test roster {target_hash}; found {len(matches)}")
    return matches[0]


def load_km(method_dir: Path):
    km_path = method_dir / "censoring_model_LOCAL_ONLY.joblib"
    if km_path.is_file():
        return abm.load_model(km_path)
    model_path = method_dir / "model_LOCAL_ONLY.joblib"
    if not model_path.is_file():
        raise FileNotFoundError(model_path)
    obj = abm.load_model(model_path)
    if isinstance(obj, dict) and "km" in obj:
        return obj["km"]
    if hasattr(obj, "km"):
        return obj.km
    raise ValueError(f"No shared censoring estimator found in {model_path}")


def development_truth(source: Path, calibration_ids, km):
    d = abm.read_npz(source / "development_LOCAL_ONLY.npz")
    ids = np.asarray(d["eid"]).astype(str)
    roles = np.asarray(d["role"]).astype(str)
    pos = {v: i for i, v in enumerate(ids)}
    ref = np.asarray(calibration_ids).astype(str)
    if any(v not in pos for v in ref):
        raise ValueError("Calibration-audit prediction IDs are absent from prepared development cohort")
    ix = np.array([pos[v] for v in ref], int)
    if not np.all(roles[ix] == "calibration_audit"):
        raise ValueError("Reliability heads may only use calibration_audit development outcomes")
    y_all, w_all = km.labels_weights(np.asarray(d["time"], float), np.asarray(d["event"]))
    return np.asarray(y_all)[ix].astype(float), np.asarray(w_all)[ix].astype(float)


def test_truth(source: Path, test_ids, km):
    t = abm.read_npz(source / "test_truth_LOCAL_ONLY.npz")
    ids = np.asarray(t["eid"]).astype(str)
    pos = {v: i for i, v in enumerate(ids)}
    ref = np.asarray(test_ids).astype(str)
    if any(v not in pos for v in ref):
        raise ValueError("Frozen test prediction IDs are absent from test truth")
    ix = np.array([pos[v] for v in ref], int)
    y, w = km.labels_weights(np.asarray(t["time"], float)[ix], np.asarray(t["event"])[ix])
    return np.asarray(y, float), np.asarray(w, float)


def logit(p):
    p = np.clip(np.asarray(p, float), 1e-7, 1 - 1e-7)
    return np.log(p / (1 - p))


def risk_features(p):
    p = np.clip(np.asarray(p, float), 1e-7, 1 - 1e-7)
    return np.c_[p, p * (1 - p), logit(p)]


INCIDENT_DIAGNOSTICS = ("eligible_match_count", "family_ess", "twin_fraction", "radius_ratio", "nearest_radius_ratio")
PREVALENT_DIAGNOSTICS = ("prevalent_matched_count", "prevalent_used_count", "prevalent_eligible_reference_count",
                        "prevalent_family_ess", "prevalent_matching_radius", "prevalent_support_available")
HEAD_NAMES = ("P0_risk_only", "P1_incident_support", "P2_incident_plus_prevalent")
HEAD_PAIRS = ((HEAD_NAMES[1], HEAD_NAMES[0]), (HEAD_NAMES[2], HEAD_NAMES[1]))


def require_diag(diagnostics, key, n):
    if key not in diagnostics:
        raise KeyError(f"Required reliability diagnostic missing: {key}")
    x = np.asarray(diagnostics[key],float)
    if x.shape != (n,):
        raise ValueError(f"{key}: wrong shape; expected ({n},)")
    if np.isinf(x).any():
        raise ValueError(f"{key}: infinite values")
    if np.any(x[np.isfinite(x)] < 0):
        raise ValueError(f"{key}: negative support values")
    return x


def incident_support_features(p, diagnostics):
    p = np.asarray(p,float); n = len(p)
    count,fam,frac,radius,nearest = [require_diag(diagnostics,key,n) for key in INCIDENT_DIAGNOSTICS]
    return np.c_[risk_features(p),np.log1p(count),np.log1p(fam),frac,np.minimum(radius,100),np.minimum(nearest,100)]


def prevalent_support_features(p, incident_diag, history_diag):
    base = incident_support_features(p,incident_diag)
    matched,used,ref,fam,radius,avail = [require_diag(history_diag,key,len(p)) for key in PREVALENT_DIAGNOSTICS]
    frac = np.divide(matched,np.maximum(ref,1),out=np.zeros_like(matched),where=ref>0)
    frac[~np.isfinite(ref) | ~np.isfinite(matched)] = np.nan
    return np.c_[base,np.log1p(ref),np.log1p(matched),np.log1p(used),np.log1p(fam),frac,np.minimum(radius,100),avail]


@dataclass
class ReliabilityHead:
    name: str
    min_events: int = MIN_EVENTS
    min_controls: int = MIN_CONTROLS
    seed: int = 2026
    risk_bins: int = RISK_BINS
    error_model: Any = None
    residual_model: Any = None
    risk_edges: np.ndarray | None = None
    thresholds: dict[float, float] | None = None
    rules: dict[float, Any] | None = None
    tie_seed: int = 2026
    status: str = "not_fitted"
    features_used: str = ""

    def fit(self, X, p, y, w, qc, coverages=DEFAULT_COVERAGES, *, ids):
        X, p, y, w, qc = np.asarray(X, float), np.asarray(p, float), np.asarray(y, float), np.asarray(w, float), np.asarray(qc, bool)
        ids = _ids(ids, "reliability calibration IDs")
        if len(ids) != len(p):
            raise ValueError("Reliability calibration IDs are misaligned")
        known = qc & np.isfinite(p) & np.isfinite(w) & (w > 0) & np.isin(y, [0, 1]) & np.all(np.isfinite(X), axis=1)
        cases, controls = int(np.sum(y[known] == 1)), int(np.sum(y[known] == 0))
        if known.sum() < 100 or cases < self.min_events or controls < self.min_controls:
            self.status = f"insufficient_calibration_audit_support:N={known.sum()};cases={cases};controls={controls}"
            return self
        q = np.linspace(0, 1, self.risk_bins + 1)
        self.risk_edges = np.unique(np.quantile(p[known], q))[1:-1]
        labels = np.searchsorted(self.risk_edges, p, side="right")
        loss = (y - p) ** 2
        residual = loss.copy()
        for b in np.unique(labels[known]):
            use = known & (labels == b)
            residual[use] -= weighted_mean(loss[use], w[use])
        leaf = max(15, min(100, int(known.sum() / 12)))
        common = dict(max_iter=250, max_leaf_nodes=5, learning_rate=.04,
                      min_samples_leaf=leaf, l2_regularization=10,
                      early_stopping=False, random_state=self.seed)
        self.error_model = HistGradientBoostingRegressor(**common).fit(X[known], loss[known], sample_weight=w[known] / w[known].mean())
        self.residual_model = HistGradientBoostingRegressor(**{**common, "random_state": self.seed + 17}).fit(
            X[known], residual[known], sample_weight=w[known] / w[known].mean())
        # Tree predictions contain genuine score ties. Including every score
        # equal to a quantile can turn a 50% rule into an 86% rule. Reuse the
        # established outcome-blind EID-hash boundary, learned only here.
        valid = qc & np.all(np.isfinite(X), axis=1)
        score = -self.residual_model.predict(X[valid])
        self.thresholds, self.rules = {}, {}
        for cov in coverages:
            cov = float(cov)
            rule = abm.FrozenRule().fit(score, np.ones(len(score), bool),
                                       ids[valid], cov, self.tie_seed)
            self.rules[cov] = rule
            self.thresholds[cov] = rule.threshold
        self.status = f"fitted:N={known.sum()};cases={cases};controls={controls}"
        return self

    def project(self, X, qc, coverages=DEFAULT_COVERAGES, *, ids):
        X, qc = np.asarray(X, float), np.asarray(qc, bool)
        ids = _ids(ids, "reliability projection IDs")
        if len(ids) != len(qc):
            raise ValueError("Reliability projection IDs are misaligned")
        valid = qc & np.all(np.isfinite(X), axis=1)
        n = len(qc)
        err = np.full(n, np.nan)
        score = np.full(n, np.nan)
        masks = {}
        if self.error_model is not None and valid.any():
            err[valid] = np.clip(self.error_model.predict(X[valid]), 0, 1)
            score[valid] = -self.residual_model.predict(X[valid])
            if self.rules is None:
                raise ValueError("Recompute reliability heads to install frozen tie rules; risk models are reusable")
            for cov in coverages:
                masks[float(cov)] = self.rules[float(cov)].apply(score, valid, ids)
        else:
            for cov in coverages:
                masks[float(cov)] = np.zeros(n, bool)
        return dict(expected_error=err, score=score, masks=masks, valid=valid)


def frozen_risk_edges(p_cal, qc_cal, bins=RISK_BINS):
    p = np.asarray(p_cal, float); qc = np.asarray(qc_cal, bool)
    v = p[qc & np.isfinite(p)]
    if len(v) < bins:
        bins = max(5, int(np.sqrt(max(len(v), 1))))
    if not len(v):
        return np.array([], dtype=float)
    return np.unique(np.quantile(v, np.linspace(0, 1, bins + 1)))[1:-1]


def within_risk_standardized_difference(loss, selected, p, w, risk_edges):
    selected = np.asarray(selected, bool)
    labels = np.searchsorted(risk_edges, p, side="right")
    good = np.isfinite(loss) & np.isfinite(p) & np.isfinite(w) & (w > 0)
    rows = []
    for b in np.unique(labels[good]):
        g = good & (labels == b)
        a, c = g & selected, g & ~selected
        if not a.any() or not c.any() or w[a].sum() <= 0 or w[c].sum() <= 0:
            continue
        rows.append((float(w[g].sum()), weighted_mean(loss[a], w[a]), weighted_mean(loss[c], w[c]),
                     int(a.sum()), int(c.sum())))
    if not rows:
        return dict(difference=np.nan, bins_used=0, selected_N=0, control_N=0)
    ww = np.array([r[0] for r in rows], float); ww /= ww.sum()
    diff = float(np.sum(ww * np.array([r[1] - r[2] for r in rows])))
    return dict(difference=diff, bins_used=len(rows), selected_N=sum(r[3] for r in rows), control_N=sum(r[4] for r in rows))


def family_bootstrap_indices(groups, seed, B):
    groups = np.asarray(groups).astype(str)
    uniq, inverse = np.unique(groups, return_inverse=True)
    rng = np.random.default_rng(seed)
    for _ in range(B):
        draw = rng.integers(0, len(uniq), len(uniq))
        counts = np.bincount(draw, minlength=len(uniq))
        yield counts[inverse].astype(float)


def ci(values, alpha=.05):
    v = np.asarray(values, float); v = v[np.isfinite(v)]
    if len(v) < 20:
        return np.nan, np.nan, len(v)
    return float(np.quantile(v, alpha / 2)), float(np.quantile(v, 1 - alpha / 2)), int(len(v))


def trapz_target_aurc(points, metric):
    # An absent ordinate invalidates this integral; never silently change its domain.
    if len(points) < 2 or any(not np.isfinite(r[metric]) for r in points):
        return np.nan
    pts = sorted((float(r["coverage_target"]),float(r[metric])) for r in points)
    x,y = map(np.asarray,zip(*pts))
    if not np.isfinite(x).all() or len(np.unique(x)) != len(x):
        raise ValueError("AURC requires a unique finite target coverage grid")
    return float(np.trapezoid(y,x))


def trapz_auc(points, metric):
    """Compatibility name: all primary AURCs use the target coverage axis."""
    return trapz_target_aurc(points,metric)


def inference_status(estimate, valid, requested):
    if not np.isfinite(estimate): return "insufficient_information"
    return "evaluated" if valid >= max(20,int(np.ceil(.8*requested))) else "insufficient_bootstrap"


def bootstrap_interval(values, estimate, requested):
    lo,hi,n = ci(values)
    status = inference_status(estimate,n,requested)
    if status != "evaluated": lo=hi=np.nan
    return dict(lower=lo,upper=hi,bootstrap_valid=n,inference_status=status)


def frozen_support_cutpoints(calibration_support):
    values=np.asarray(calibration_support,float)
    values=values[np.isfinite(values)]
    if not len(values): return np.array([],float)
    # Value cutpoints preserve ties; empty groups are allowed.
    return np.unique(np.quantile(values,[.25,.5,.75]))


def support_rank_within_risk(support, p, w, edges, *, support_cutpoints):
    """Optional descriptive groups, with cutpoints frozen in calibration_audit."""
    support,p,w=map(lambda v:np.asarray(v,float),(support,p,w))
    result=np.full(len(p),np.nan)
    good=np.isfinite(support)&np.isfinite(p)&np.isfinite(w)&(w>0)
    result[good]=np.searchsorted(support_cutpoints,support[good],side="right")+1
    return result


def fixed_effect_continuous_support_trend(loss, support, p, w, edges):
    support=np.asarray(support,float); p=np.asarray(p,float); w=np.asarray(w,float)
    x=np.full(len(support),np.nan); ok_support=np.isfinite(support)&(support>=0)
    x[ok_support]=np.log1p(support[ok_support]); y=np.asarray(loss,float)
    labels=np.searchsorted(edges,p,side="right")
    good=np.isfinite(x)&np.isfinite(y)&np.isfinite(p)&np.isfinite(w)&(w>0)
    xr=np.full(len(x),np.nan);yr=np.full(len(y),np.nan)
    for b in np.unique(labels[good]):
        g=good&(labels==b)
        if g.sum()<4:continue
        xr[g]=x[g]-weighted_mean(x[g],w[g]);yr[g]=y[g]-weighted_mean(y[g],w[g])
    ok=good&np.isfinite(xr)&np.isfinite(yr)
    den=np.sum(w[ok]*xr[ok]**2)
    # Avoid interpreting round-off from a constant support variable as a slope.
    if not ok.any() or np.ptp(x[ok])==0:den=0
    slope=np.sum(w[ok]*xr[ok]*yr[ok])/den if den>np.finfo(float).eps else np.nan
    return float(slope),int(ok.sum())


def risk_coverage_rows(heads, projected, y, w, p, qc):
    rows = []
    for name, proj in projected.items():
        if heads[name].error_model is None:
            continue
        for cov, mask in sorted(proj["masks"].items()):
            m = mask & qc
            row = dict(head=name, coverage_target=cov, actual_coverage=float(np.mean(m[qc])) if qc.any() else np.nan,
                       mean_expected_error=weighted_mean(proj["expected_error"][m], w[m]) if m.any() else np.nan)
            row["coverage_deviation"] = row["actual_coverage"] - cov
            row["coverage_calibration_warning"] = bool(np.isclose(cov,PRIMARY_COVERAGE) and abs(row["coverage_deviation"])>.03)
            row.update(point_metrics(y, w, p, m))
            rows.append(row)
    for name in heads:
        selected=[row for row in rows if row["head"]==name]
        if selected:
            maximum=max(abs(row["coverage_deviation"]) for row in selected)
            for row in selected:row["max_abs_coverage_deviation"]=maximum
    return rows


def bootstrap_aurc_contrasts(heads, projected, y, w, p, qc, groups, B, seed):
    point_rows=risk_coverage_rows(heads,projected,y,w,p,qc)
    by_head={h:[r for r in point_rows if r["head"]==h] for h in heads}
    points=[]
    for a,b in HEAD_PAIRS:
        if not by_head.get(a) or not by_head.get(b):continue
        if [r["coverage_target"] for r in by_head[a]] != [r["coverage_target"] for r in by_head[b]]:
            raise ValueError("Paired AURC requires identical target grids")
        for metric in ("Brier","LogLoss"):
            va,vb=trapz_target_aurc(by_head[a],metric),trapz_target_aurc(by_head[b],metric)
            rows=by_head[a]+by_head[b]
            points.append(dict(contrast=a+"_minus_"+b,model=a,reference=b,metric="AURC_"+metric,
                difference=va-vb,model_aurc=va,reference_aurc=vb,coverage_axis="predeclared_target",
                coverage_min=min(r["coverage_target"] for r in rows),coverage_max=max(r["coverage_target"] for r in rows),
                max_abs_coverage_deviation=max(abs(r["coverage_deviation"]) for r in rows),
                coverage_calibration_warning=any(r["coverage_calibration_warning"] for r in rows)))
    boots={r["contrast"]+r["metric"]:[] for r in points}
    losses={"Brier":(y-p)**2,"LogLoss":logloss_rows(y,p)}
    for mult in family_bootstrap_indices(groups,seed,B):
        wb=w*mult; curves={}
        for h,proj in projected.items():
            if heads[h].error_model is None:continue
            curves[h]=[]
            for cov,mask in sorted(proj["masks"].items()):
                use=mask&qc
                curves[h].append(dict(coverage_target=cov,**{m:weighted_mean(loss[use],wb[use]) for m,loss in losses.items()}))
        for row in points:
            metric=row["metric"].removeprefix("AURC_")
            boots[row["contrast"]+row["metric"]].append(trapz_target_aurc(curves[row["model"]],metric)-trapz_target_aurc(curves[row["reference"]],metric))
    for row in points:
        row.update(bootstrap_interval(boots[row["contrast"]+row["metric"]],row["difference"],B),
            bootstrap_unit="family",inference="conditional on frozen risk model and reliability heads")
    return point_rows,points


def standardized_contrasts(heads, projected, y, w, p, qc, groups, risk_edges, B, seed, coverage=PRIMARY_COVERAGE):
    rows = []
    w = np.where(qc, w, 0.0)
    losses = {"Brier": (y-p)**2, "LogLoss": logloss_rows(y,p)}
    keys = [k for k in ("P0_risk_only", "P1_incident_support", "P2_incident_plus_prevalent") if k in projected and heads[k].error_model is not None]
    point = {}
    for name in keys:
        sel = projected[name]["masks"][float(coverage)] & qc
        for metric, loss in losses.items():
            res = within_risk_standardized_difference(loss, sel, p, w, risk_edges)
            key=(name,metric); point[key]=res
            rows.append(dict(head=name, coverage_target=coverage, metric=metric,
                             difference_selected_minus_nonselected=res["difference"], bins_used=res["bins_used"],
                             selected_N=res["selected_N"], nonselected_N=res["control_N"],
                             selected_mean_risk=weighted_mean(p[sel],w[sel]),
                             nonselected_mean_risk=weighted_mean(p[~sel & qc],w[~sel & qc])))
    vals={k:[] for k in point}
    for mult in family_bootstrap_indices(groups, seed+101, B):
        wb=w*mult
        for name in keys:
            sel=projected[name]["masks"][float(coverage)] & qc
            for metric,loss in losses.items():
                vals[(name,metric)].append(within_risk_standardized_difference(loss,sel,p,wb,risk_edges)["difference"])
    for r in rows:
        r.update(bootstrap_interval(vals[(r["head"],r["metric"])],r["difference_selected_minus_nonselected"],B),risk_bins=len(risk_edges)+1,
                 primary_interpretation="negative means high-predictability selection has lower error at matched risk strata")
    return rows


def standardized_selected_loss(loss, selected, p, w, risk_edges):
    loss,p,w=np.asarray(loss,float),np.asarray(p,float),np.asarray(w,float)
    selected=np.asarray(selected,bool);labels=np.searchsorted(risk_edges,p,side="right")
    good=np.isfinite(loss)&np.isfinite(p)&np.isfinite(w)&(w>0)
    estimates=[];mass=[]
    for b in np.unique(labels[good]):
        g=good&(labels==b);use=g&selected
        if use.any():estimates.append(weighted_mean(loss[use],w[use]));mass.append(w[g].sum())
    if not mass:return np.nan,0
    return float(np.average(estimates,weights=mass)),len(mass)


def _paired_standardized_loss(loss, a, b, p, w, edges, common_bins=None):
    labels=np.searchsorted(edges,p,side="right")
    good=np.isfinite(loss)&np.isfinite(p)&np.isfinite(w)&(w>0)
    bins=[int(k) for k in np.unique(labels[good]) if (good&a&(labels==k)).any() and (good&b&(labels==k)).any()] if common_bins is None else common_bins
    mass=[];av=[];bv=[]
    for k in bins:
        g=good&(labels==k)
        if not g.any():continue
        if not (g&a).any() or not (g&b).any():return np.nan,np.nan,bins,np.nan
        mass.append(w[g].sum());av.append(weighted_mean(loss[g&a],w[g&a]));bv.append(weighted_mean(loss[g&b],w[g&b]))
    if not mass:return np.nan,np.nan,bins,0.
    return float(np.average(av,weights=mass)),float(np.average(bv,weights=mass)),bins,float(sum(mass)/w[good].sum())


def _risk_balance(p, a, b, w):
    good=np.isfinite(p)&np.isfinite(w)&(w>0);a=a&good;b=b&good
    ma,mb=weighted_mean(p[a],w[a]),weighted_mean(p[b],w[b])
    va=weighted_mean((p[a]-ma)**2,w[a]);vb=weighted_mean((p[b]-mb)**2,w[b])
    scale=np.sqrt((va+vb)/2)
    smd=(ma-mb)/scale if np.isfinite(scale) and scale>0 else 0. if ma==mb else np.nan
    if not a.any() or not b.any():return ma,mb,smd,np.nan
    grid=np.sort(np.unique(p[a|b]));cdf=[]
    for mask in (a,b):
        order=np.argsort(p[mask]);values=p[mask][order];weights=w[mask][order]
        cumulative=np.r_[0,np.cumsum(weights)/weights.sum()]
        cdf.append(cumulative[np.searchsorted(values,grid,side="right")])
    return ma,mb,float(smd),float(np.max(np.abs(cdf[0]-cdf[1])))


def direct_standardized_contrasts(heads, projected, y, w, p, qc, groups, risk_edges, B, seed, coverage=PRIMARY_COVERAGE):
    w=np.where(qc,w,0.);losses={"Brier":(y-p)**2,"LogLoss":logloss_rows(y,p)}
    rows=[];specs=[]
    for model,reference in HEAD_PAIRS:
        if any(k not in heads or heads[k].error_model is None for k in (model,reference)):continue
        a=projected[model]["masks"][float(coverage)]&qc;b=projected[reference]["masks"][float(coverage)]&qc
        ma,mb,smd,ks=_risk_balance(p,a,b,w)
        for metric,loss in losses.items():
            va,vb,bins,mass=_paired_standardized_loss(loss,a,b,p,w,risk_edges)
            union=int((a|b).sum());overlap=int((a&b).sum())
            rows.append(dict(contrast=model+"_minus_"+reference,model=model,reference=reference,metric=metric,
                coverage_target=coverage,model_value=va,reference_value=vb,difference=va-vb,
                model_selected_N=int(a.sum()),reference_selected_N=int(b.sum()),overlap_N=overlap,
                Jaccard=overlap/union if union else np.nan,model_mean_risk=ma,reference_mean_risk=mb,
                risk_SMD=smd,risk_KS=ks,bins_used=len(bins),reference_weight_fraction=mass,
                standardization="shared populated risk strata; common all-person IPCW weights",
                bootstrap_unit="family",coverage_calibration_warning=bool(qc.any() and max(abs(a[qc].mean()-coverage),abs(b[qc].mean()-coverage))>.03)))
            specs.append((loss,a,b,bins))
    draws=[[] for _ in rows]
    for mult in family_bootstrap_indices(groups,seed+101,B):
        for values,(loss,a,b,bins) in zip(draws,specs):
            va,vb,_,_=_paired_standardized_loss(loss,a,b,p,w*mult,risk_edges,bins)
            values.append(va-vb)
    for row,values in zip(rows,draws):row.update(bootstrap_interval(values,row["difference"],B))
    return rows


def support_trends(diagnostics, y, w, p, qc, groups, risk_edges, B, seed):
    rows=[];w=np.where(qc,w,0.)
    for key in ("eligible_match_count","family_ess"):
        try:support=require_diag(diagnostics,key,len(p))
        except (KeyError,ValueError):
            for metric in ("Brier","LogLoss"):
                rows.append(dict(support_metric=key,transformation="log1p",risk_adjustment="frozen risk-bin fixed effects",
                    error_metric=metric,slope=np.nan,lower=np.nan,upper=np.nan,bootstrap_valid=0,N_used=0,
                    inference_status="unavailable_missing_or_invalid_support_diagnostics"))
            continue
        for metric,loss in {"Brier":(y-p)**2,"LogLoss":logloss_rows(y,p)}.items():
            slope,n=fixed_effect_continuous_support_trend(loss,support,p,w,risk_edges)
            boot=[fixed_effect_continuous_support_trend(loss,support,p,w*mult,risk_edges)[0]
                  for mult in family_bootstrap_indices(groups,seed+307,B)]
            rows.append(dict(support_metric=key,transformation="log1p",risk_adjustment="frozen risk-bin fixed effects",
                error_metric=metric,slope=slope,N_used=n,bootstrap_unit="family",**bootstrap_interval(boot,slope,B)))
    return rows


def _head_features(name, glm, history):
    p=np.asarray(glm["probability"],float);diag=glm.get("diagnostics",{})
    if name==HEAD_NAMES[0]:return risk_features(p)
    if name==HEAD_NAMES[1]:return incident_support_features(p,diag)
    if history is None:raise KeyError("No prevalent projection")
    return prevalent_support_features(p,diag,history.get("diagnostics",{}))


def _unavailable_head(head, error):
    side="incident" if head.name==HEAD_NAMES[1] else "prevalent"
    reason="missing" if isinstance(error,KeyError) else "invalid"
    head.status=f"unavailable_{reason}_{side}_support_diagnostics"
    head.error_model=head.residual_model=None;head.thresholds=None;head.rules=None
    return head


def fit_heads(glm_cal, hist_cal, ycal, wcal, coverages, seed, risk_bins=RISK_BINS):
    p=np.asarray(glm_cal["probability"],float)
    qc=np.asarray(glm_cal.get("technical_qc",np.ones(len(p),bool)),bool)
    heads={}
    features=["p;p*(1-p);logit(p)",
        "P0 + log1p(eligible_match_count);log1p(family_ess);twin_fraction;radius_ratio;nearest_radius_ratio",
        "P1 + log1p(prevalent_eligible_reference_count);log1p(prevalent_matched_count);log1p(prevalent_used_count);log1p(prevalent_family_ess);matched/reference fraction;prevalent_matching_radius;prevalent_support_available"]
    for i,name in enumerate(HEAD_NAMES):
        head=ReliabilityHead(name=name,seed=seed+100*i,tie_seed=seed+91,risk_bins=risk_bins,features_used=features[i])
        try:
            X=_head_features(name,glm_cal,hist_cal)
            heads[name]=head.fit(X,p,ycal,wcal,qc,coverages,ids=glm_cal["eid"])
        except (KeyError,ValueError) as error:
            if name==HEAD_NAMES[0]:raise
            heads[name]=_unavailable_head(head,error)
    return heads


def project_heads(heads, glm_test, hist_test, coverages):
    p=np.asarray(glm_test["probability"],float)
    qc=np.asarray(glm_test.get("technical_qc",np.ones(len(p),bool)),bool)
    projected={}
    for name,head in heads.items():
        if head.error_model is None:continue
        try:projected[name]=head.project(_head_features(name,glm_test,hist_test),qc,coverages,ids=glm_test["eid"])
        except (KeyError,ValueError) as error:
            if name==HEAD_NAMES[0]:raise
            _unavailable_head(head,error)
    return projected


def primary_reliability_status(heads):
    fitted=lambda name:name in heads and heads[name].error_model is not None
    if fitted(HEAD_NAMES[0]) and fitted(HEAD_NAMES[1]):return "completed_primary_P1_vs_P0"
    if fitted(HEAD_NAMES[0]):return "insufficient_primary_P1_support"
    return "insufficient_primary_reliability_support"


def load_reliability_protocol(path=PROTOCOL_PATH):
    path=Path(path);protocol=json_read(path)
    if protocol.get("no_test_tuning") is not True or protocol.get("bootstrap_unit")!="family":
        raise ValueError("Reliability protocol must preserve family bootstrap and no-test-tuning")
    return protocol,sha_file(path)


def discovery_status(scope, protocol, host_contract):
    trait,layer=scope.parent.parent.name,scope.parent.name
    if layer!=protocol["primary_layer"]:return "secondary_layer"
    if trait in protocol["discovery_traits"]:return "discovery"
    if trait in protocol["locked_replication_traits"]:
        # Existing historical predictions do not prove their results were unseen.
        if host_contract.get("test_results_previously_inspected") is False:
            return "locked_replication"
        return "replication_candidate_inspection_unverified"
    return "exploratory"


def prevalent_support_audit(hist_dir, hist_test):
    if hist_dir is None or hist_test is None:
        return [dict(status="unavailable_no_prevalent_projection")]
    manifest=json_read(Path(hist_dir)/"MODEL_FROZEN.json")
    audit=manifest.get("source_transform",{});n=audit.get("source_used_N",0)
    # Frozen source metadata only; never derive a source sample from test outcomes.
    known_duration=audit.get("source_used_duration_known_N")
    unknown_med=audit.get("source_used_medication_unknown_N")
    families=audit.get("source_families")
    if families is None:
        obj=abm.load_model(Path(hist_dir)/"model_LOCAL_ONLY.joblib")
        transformer=obj.get("transformer") if isinstance(obj,dict) else None
        group=getattr(transformer,"source_groups_",None)
        if group is not None:families=len(np.unique(group))
    diag=hist_test.get("diagnostics",{});length=len(hist_test["eid"])
    def median(key):
        try:v=require_diag(diag,key,length);return float(np.nanmedian(v)) if np.isfinite(v).any() else np.nan
        except (KeyError,ValueError):return np.nan
    try:available=require_diag(diag,"prevalent_support_available",length)
    except (KeyError,ValueError):available=None
    return [dict(status=audit.get("status","unavailable_source_audit"),source_N=n,source_families=families,
        duration_known_pct=100*known_duration/n if n and known_duration is not None else np.nan,
        medication_known_pct=100*(n-unknown_med)/n if n and unknown_med is not None else np.nan,
        support_available_pct=100*float(np.mean(available)) if available is not None else np.nan,
        median_eligible_reference_N=median("prevalent_eligible_reference_count"),
        median_matched_N=median("prevalent_matched_count"),median_family_ESS=median("prevalent_family_ess"))]


def find_history_method(scope: Path):
    d=locate_method(scope,"abm_prevalent")
    return d


def write_tables(out, tables):
    for name, frame in tables.items():
        if any(str(c).lower() in {"eid", "person_id", "participant_id"} for c in frame.columns):
            raise AssertionError("Participant identifier in public reliability table")
        frame.to_csv(out/(name+".tsv"), sep="\t", index=False)
    with pd.ExcelWriter(out/"C1_reliability_evaluation.xlsx", engine="openpyxl") as writer:
        for name, frame in tables.items():
            frame.to_excel(writer, sheet_name=name[:31], index=False)


def evaluate_scope(scope: Path, bootstrap: int, seed: int, coverages, primary_coverage: float,
                   risk_bins: int, replace: bool, *, method_paths=None, prepared_source=None, protocol_path=PROTOCOL_PATH):
    """Evaluate the exact frozen methods from the caller, or a completed local run.

    FINAL is a stable output location. Reuse compares data and settings only;
    editing a Python source file never invalidates a risk model or these outputs.
    """
    scope = Path(scope)
    root = scope/"reliability"
    root.mkdir(parents=True, exist_ok=True)
    # This lock also protects standalone calls against an integrated C1 run.
    with (root/".evaluation.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _evaluate_scope(scope, bootstrap, seed, coverages, primary_coverage, risk_bins,
                               replace, method_paths, prepared_source, protocol_path)


def _evaluate_scope(scope, bootstrap, seed, coverages, primary_coverage, risk_bins,
                    replace, method_paths, prepared_source, protocol_path):
    if bootstrap < 1 or risk_bins < 1 or not coverages or primary_coverage not in coverages:
        raise ValueError("Positive bootstrap/risk bins and a declared primary coverage are required")
    if any(not 0 < c <= 1 for c in coverages):
        raise ValueError("Coverage must be in (0, 1]")
    protocol,protocol_hash = load_reliability_protocol(protocol_path)
    if (list(map(float,coverages)) != protocol["coverage_grid"] or primary_coverage != protocol["primary_coverage"] or risk_bins != protocol["risk_bins"]):
        raise ValueError("Reliability settings differ from the locked protocol; supply an explicit separate protocol for a new study")
    previous_contract=scope/"reliability/FINAL/reliability_contract.json"
    if previous_contract.is_file():
        previous=json_read(previous_contract)
        if (previous.get("protocol_sha256") and previous["protocol_sha256"]!=protocol_hash
                and previous.get("test_outcomes_opened")):
            raise ValueError("The reliability protocol was already evaluated here; use a separate study root instead of changing it after test access")
    glm_dir = method_paths.get("glmnet") if method_paths is not None else locate_method(scope, "glmnet")
    if glm_dir is None:
        return dict(status="skipped_missing_glmnet", scope=str(scope))
    glm_dir = Path(glm_dir)
    from c1_correlate import _verified_stored_complete
    if _verified_stored_complete(glm_dir) is None:
        raise ValueError("Frozen glmnet integrity check failed")
    glm_cal = abm.load_model(glm_dir/"calibration_audit_predictions_LOCAL_ONLY.joblib")
    glm_test = abm.load_model(glm_dir/"test_predictions_LOCAL_ONLY.joblib")
    glm_cal = align_prediction(glm_cal["eid"], glm_cal)
    glm_test = align_prediction(glm_test["eid"], glm_test)
    test_ids = np.asarray(glm_test["eid"]).astype(str)
    source = Path(prepared_source) if prepared_source is not None else locate_prepared_source(scope, test_ids)
    host_contract = json_read(source/"host_contract.json")
    if host_contract.get("incident_test_ids_hash") != abm.fingerprint(test_ids):
        raise ValueError("Prepared input test roster differs from frozen predictions")
    family_verified = host_contract.get("independent_family_support_verified",
        host_contract.get("config", {}).get("grouping_unit") == "family")
    if not family_verified or "group" not in glm_test:
        raise ValueError("Verified family components are required for reliability bootstrap")
    groups = np.asarray(glm_test["group"]).astype(str)
    if groups.shape != test_ids.shape or np.any(np.isin(groups, ["", "nan", "None"])):
        raise ValueError("Missing or misaligned family components")
    km = load_km(glm_dir)
    ycal, wcal = development_truth(source, glm_cal["eid"], km)
    hist_dir = (method_paths.get("abm_prevalent") if method_paths is not None
                else find_history_method(scope))
    hist_cal = hist_test = None
    if hist_dir is not None:
        hist_dir = Path(hist_dir)
        if _verified_stored_complete(hist_dir) is None:
            raise ValueError("Frozen prevalent model integrity check failed")
        hist_cal = align_prediction(glm_cal["eid"], abm.load_model(hist_dir/"calibration_audit_predictions_LOCAL_ONLY.joblib"))
        hist_test = align_prediction(test_ids, abm.load_model(hist_dir/"test_predictions_LOCAL_ONLY.joblib"))
    heads = fit_heads(glm_cal, hist_cal, ycal, wcal, coverages, seed, risk_bins)
    projected = project_heads(heads, glm_test, hist_test, coverages)
    ptest = np.asarray(glm_test["probability"], float)
    common = np.asarray(glm_test.get("technical_qc", np.ones(len(ptest), bool)), bool) & np.isfinite(ptest)
    risk_edges = frozen_risk_edges(np.asarray(glm_cal["probability"], float),
        np.asarray(glm_cal.get("technical_qc", np.ones(len(ycal), bool)), bool), risk_bins)
    status = [dict(component="reliability_head", head=name, status=h.status, features=h.features_used,
                   calibration_audit_N=len(ycal), calibration_audit_cases=int(np.sum((wcal>0)&(ycal==1))))
              for name, h in heads.items()]
    selection = {name: {str(cov): dict(selected_N=int(mask.sum()),
                 selected_roster_hash=abm.fingerprint(test_ids[mask]))
                 for cov, mask in proj["masks"].items()} for name, proj in projected.items()}
    contract = dict(version=VERSION, protocol=protocol, protocol_sha256=protocol_hash,
        discovery_or_replication=discovery_status(scope,protocol,host_contract), scope=str(scope), glmnet_method=str(glm_dir),
        prevalent_method=str(hist_dir) if hist_dir else None, prepared_source=str(source),
        bootstrap=bootstrap, bootstrap_unit="family", seed=seed,
        coverages=list(map(float, coverages)), primary_coverage=float(primary_coverage),
        risk_bins=int(risk_bins), risk_edges=risk_edges.tolist(),
        risk_edges_source="calibration_audit_predictions",
        selection_rule="calibration_frozen_score_and_eid_hash_tie_boundary_v1",
        heads={k:dict(status=v.status, features_used=v.features_used, thresholds=v.thresholds,
                      threshold_ties={str(c): r.tie for c, r in (v.rules or {}).items()}) for k,v in heads.items()},
        frozen_test_selections=selection, test_outcomes_opened=False)
    out = scope/"reliability"/"FINAL"
    with tempfile.TemporaryDirectory(prefix=".pending-", dir=out.parent) as tmp:
        stage = Path(tmp)/"result"; stage.mkdir()
        private = stage/"_LOCAL_ONLY"; private.mkdir(mode=0o700)
        joblib.dump(dict(heads=heads, risk_edges=risk_edges, projected=projected, contract=contract),
                    private/"reliability_heads.joblib", compress=3)
        write_json(stage/"FROZEN_BEFORE_TEST.json", contract)
        # No test-truth file is loaded or hashed until this frozen record exists.
        supported = any(h.error_model is not None for h in heads.values())
        result_status = primary_reliability_status(heads)
        data_files = {"glmnet_cal": sha_file(glm_dir/"calibration_audit_predictions_LOCAL_ONLY.joblib"),
                      "glmnet_test": sha_file(glm_dir/"test_predictions_LOCAL_ONLY.joblib"),
                      "development": sha_file(source/"development_LOCAL_ONLY.npz"),
                      "censoring_model": sha_file(glm_dir/("censoring_model_LOCAL_ONLY.joblib" if (glm_dir/"censoring_model_LOCAL_ONLY.joblib").is_file() else "model_LOCAL_ONLY.joblib"))}
        if hist_dir is not None:
            for role in ("test", "calibration_audit"):
                data_files["prevalent_"+role] = sha_file(hist_dir/(role+"_predictions_LOCAL_ONLY.joblib"))
        if supported:
            data_files["test_truth"] = sha_file(source/"test_truth_LOCAL_ONLY.npz")
        request = dict(contract=contract, data=data_files)
        # JSON round trip normalizes float keys in threshold dictionaries.
        request = json.loads(json.dumps(request, default=str))
        marker = out/"COMPLETE.json"
        if marker.is_file() and not replace:
            saved = json_read(marker)
            if saved.get("request") == request and saved.get("files") and all(
                    (out/n).is_file() and sha_file(out/n)==value for n,value in saved["files"].items()):
                print(f"[LE8] DONE C1 REUSE RELIABILITY {scope.parent.parent.name}/{scope.parent.name}",flush=True)
                return dict(status=saved["status"], reused=True, scope=str(scope), output=str(out))
        print(f"[LE8] START C1 RECOMPUTE RELIABILITY {scope.parent.parent.name}/{scope.parent.name}",flush=True)
        tables = {}
        if supported:
            # First semantic test-truth read, after models, masks and bins are saved.
            y, w = test_truth(source, test_ids, km)
            population = point_metrics(y,w,ptest,np.isfinite(ptest))
            w = np.where(common, w, 0.0)
            risk_rows, aurc = bootstrap_aurc_contrasts(heads, projected, y, w, ptest, common, groups, bootstrap, seed)
            std = standardized_contrasts(heads, projected, y, w, ptest, common, groups, risk_edges, bootstrap, seed, primary_coverage)
            increment = direct_standardized_contrasts(heads,projected,y,w,ptest,common,groups,risk_edges,bootstrap,seed,primary_coverage)
            trend = support_trends(glm_test["diagnostics"], y, w, ptest, common, groups, risk_edges, bootstrap, seed)
            tables = {"reliability_population_performance":pd.DataFrame([dict(model="glmnet",population="all frozen test participants",**population)]),
                      "risk_coverage":pd.DataFrame(risk_rows), "aurc_contrasts":pd.DataFrame(aurc),
                      "risk_standardized_selection":pd.DataFrame(std), "support_trend":pd.DataFrame(trend),
                      "reliability_increment_50":pd.DataFrame(increment)}
            for row in status:
                row.update(test_N=int(common.sum()), test_cases=int(np.sum((w>0)&(y==1))))
        tables["reliability_prevalent_support_audit"] = pd.DataFrame(prevalent_support_audit(hist_dir,hist_test))
        tables["reliability_status"] = pd.DataFrame(status).assign(primary_status=result_status)
        for name in ("risk_coverage","aurc_contrasts","reliability_increment_50","risk_standardized_selection","support_trend"):
            if name not in tables or tables[name].empty:
                tables[name]=pd.DataFrame([dict(status=result_status,inference_status="unavailable")])
        write_tables(stage, tables)
        public = dict(contract, status=result_status, test_outcomes_opened=supported,
            primary_hypothesis="P1 incident support versus P0 risk-only, with identical frozen glmnet risk",
            secondary_hypothesis="P2 prevalent support versus P1 incident support",
            primary_metrics=["IPCW Brier AURC", "IPCW log-loss AURC", "50% risk-standardized Brier/log-loss"],
            discovery_note="HTN is discovery; COPD is replication only if test performance was not previously inspected",
            inference="family bootstrap conditional on fitted risk models and reliability heads")
        write_json(stage/"reliability_contract.json", public)
        files = {str(p.relative_to(stage)):sha_file(p) for p in stage.rglob("*") if p.is_file()}
        write_json(stage/"COMPLETE.json", dict(request=request, files=files, version=VERSION, status=result_status))
        # Publish atomically, retaining the old result if staging fails.
        previous = Path(tmp)/"previous"
        if out.exists(): out.rename(previous)
        try: stage.rename(out)
        except BaseException:
            if previous.exists(): previous.rename(out)
            raise
    return dict(status=result_status, scope=str(scope), output=str(out), heads={k:v.status for k,v in heads.items()})


def parser():
    p=argparse.ArgumentParser(description=__doc__,allow_abbrev=False)
    p.add_argument("--analysis-root","--dir-analysis-root",type=Path,default=Path(os.getenv("LE8_ANALYSIS_ROOT","/mnt/d/analysis/le8")))
    p.add_argument("--traits","--Y",dest="traits",default="all",help="CSV trait names or 'all' to scan analysis-root")
    p.add_argument("--layer","--biom",dest="layer",default="prot",help="CSV prot,met")
    p.add_argument("--bootstrap",type=int,default=1000)
    p.add_argument("--seed",type=int,default=2026)
    p.add_argument("--cores",type=int,default=8)
    p.add_argument("--coverages",default=",".join(map(str,DEFAULT_COVERAGES)))
    p.add_argument("--primary-coverage",type=float,default=PRIMARY_COVERAGE)
    p.add_argument("--risk-bins",type=int,default=RISK_BINS)
    p.add_argument("--replace",nargs="?",const="TRUE",choices=["TRUE","FALSE"],default="FALSE",help="Recompute reliability aggregates only")
    p.add_argument("--reliability-protocol",type=Path,default=PROTOCOL_PATH)
    p.add_argument("--dry-run",action="store_true")
    p.add_argument("--preflight",action="store_true")
    p.add_argument("--strict",action="store_true",help="Stop on the first trait error")
    return p


def main(argv=None):
    os.umask(0o077);args=parser().parse_args(argv)
    import re
    if args.bootstrap<20 or args.risk_bins<10 or args.cores<1:raise ValueError("bootstrap>=20, risk-bins>=10 and cores>=1 required")
    coverages=tuple(float(x) for x in args.coverages.split(","))
    protocol,_=load_reliability_protocol(args.reliability_protocol)
    if list(coverages)!=protocol["coverage_grid"] or args.risk_bins!=protocol["risk_bins"] or args.primary_coverage!=protocol["primary_coverage"]:
        raise ValueError("Options differ from the frozen reliability protocol")
    layers=args.layer.split(",")
    if not set(layers)<={"prot","met"}:raise ValueError("Invalid layers")
    traits=sorted(d.name for d in args.analysis_root.iterdir() if d.is_dir() and any((d/l/"c1_correlate").is_dir() for l in layers)) if args.traits=="all" else args.traits.split(",")
    if any(not re.fullmatch(r"[A-Za-z0-9_]+",t) for t in traits):raise ValueError("Invalid traits")
    if args.dry_run or args.preflight:
        print(json.dumps(dict(action="reliability_only",traits=traits,layers=layers,risk_model_refit=False,
            protocol=protocol,bootstrap=args.bootstrap,replace=args.replace),indent=2));return 0
    summary=[]
    from threadpoolctl import threadpool_limits
    for trait in traits:
        for layer in layers:
            scope=args.analysis_root/trait/layer/"c1_correlate"
            try:
                if not scope.is_dir():result=dict(status="skipped_missing_scope",scope=str(scope))
                else:
                    with threadpool_limits(args.cores):
                        result=evaluate_scope(scope,args.bootstrap,args.seed,coverages,args.primary_coverage,
                            args.risk_bins,args.replace=="TRUE",protocol_path=args.reliability_protocol)
                summary.append(dict(result,trait=trait,layer=layer))
                print(f"[LE8 C1] RELIABILITY {trait}/{layer}: {result['status']}",flush=True)
            except Exception as error:
                summary.append(dict(trait=trait,layer=layer,status="error",error=str(error)))
                print(f"[LE8 C1] RELIABILITY ERROR {trait}/{layer}: {error}",file=sys.stderr,flush=True)
                if args.strict:raise
    from c1_report import write_reliability_cross_trait_summary
    for layer in layers:write_reliability_cross_trait_summary(args.analysis_root,layer)
    return 1 if any(r["status"]=="error" for r in summary) else 0


if __name__ == "__main__":
    import sys
    from c1_evaluation import main as reliability_main
    sys.exit(reliability_main())
