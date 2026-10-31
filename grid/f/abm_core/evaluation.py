"""Outcome-blind matching and evaluation for the two ABM objectives.

This module does not fit a predictor, selector, threshold, or matching stratum.
Those must be frozen using development data before calling it. In particular,
``selected`` and ``match_stratum`` may not be constructed from query outcomes.

Positive gain means reduced loss. Negative predictability difference means
that the selected people have lower *baseline* loss than matched controls.
The confidence intervals condition on the frozen models, selections and
matching; they do not include model-fitting or threshold-selection uncertainty.
"""
from __future__ import annotations

import hashlib
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from .losses import row_loss


def _bool_column(values: pd.Series, name: str) -> np.ndarray:
    if values.isna().any() or not values.isin([True, False, 0, 1]).all():
        raise ValueError(f"{name} must contain only nonmissing booleans/0/1")
    return values.astype(bool).to_numpy()


def _identity(frame: pd.DataFrame, group_col: str) -> tuple[np.ndarray, np.ndarray]:
    for col in ("id", group_col):
        if col not in frame or frame[col].isna().any():
            raise ValueError(f"Missing column or missing identity: {col}")
    ids = frame["id"].astype(str).to_numpy()
    groups = frame[group_col].astype(str).to_numpy()
    if len(np.unique(ids)) != len(ids) or np.any(ids == "") or np.any(groups == ""):
        raise ValueError("IDs must be unique and IDs/groups must be nonempty")
    return ids, groups


def _rank(identifier: str, seed: int, kind: str) -> bytes:
    return hashlib.sha256(f"{seed}\0{kind}\0{identifier}".encode("utf-8")).digest()


def make_matched_controls(
    frame: pd.DataFrame,
    *,
    seed: int = 20261020,
    group_col: str = "group_id",
    selected_col: str = "selected",
    match_strata_col: str = "match_stratum",
    control_eligible_col: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Choose disjoint, same-stratum controls without reading any outcome.

    A control may not belong to *any* selected family/group. Matching is 1:1
    without replacement. In an undersupplied stratum, a deterministic hash
    selects a subset of selected IDs and available controls. There is no
    outcome-dependent rematching, risk-bin expansion, or backfilling with
    selected people. ``match_stratum`` encodes frozen baseline-risk bins;
    ancestry labels do not restrict matches. Missing strata remain unmatched.
    An optional, precomputed ``control_eligible_col`` restricts the control
    pool (for example, supported & assessment_available). That eligibility
    must be frozen without query outcomes. None means all unselected queries
    are initially eligible, subject to strata and family exclusions.

    The returned positions refer to the current row order; IDs and assignments
    are invariant to permutations of that order. Pairs are individual-level
    artifacts and should remain in the project's local-only output area.
    """
    ids, groups = _identity(frame, group_col)
    if selected_col not in frame or match_strata_col not in frame:
        raise ValueError("Matching requires selected and frozen match-stratum columns")
    selected = _bool_column(frame[selected_col], selected_col)
    raw_strata = frame[match_strata_col]
    present = raw_strata.notna().to_numpy()
    strata = raw_strata.fillna("").astype(str).to_numpy()
    selected_groups = set(groups[selected])
    outside_selected_groups = np.array(
        [g not in selected_groups for g in groups], dtype=bool
    )
    if control_eligible_col is not None:
        if control_eligible_col not in frame:
            raise ValueError(f"Control eligibility column missing: {control_eligible_col}")
        control_allowed = _bool_column(frame[control_eligible_col], control_eligible_col)
    else:
        control_allowed = np.ones(len(frame), dtype=bool)
    eligible_control = (~selected) & outside_selected_groups & present & control_allowed
    records: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for stratum in sorted(set(strata[selected & present])):
        left = np.flatnonzero(selected & present & (strata == stratum))
        right = np.flatnonzero(eligible_control & (strata == stratum))
        left = sorted(left, key=lambda i: (_rank(ids[i], seed, "selected"), ids[i]))
        right = sorted(right, key=lambda i: (_rank(ids[i], seed, "control"), ids[i]))
        n = min(len(left), len(right))
        diagnostics.append({
            "match_stratum": stratum,
            "n_selected": len(left),
            "n_available_controls": len(right),
            "n_pairs": n,
            "n_unmatched_selected": len(left) - n,
        })
        for a, b in zip(left[:n], right[:n]):
            records.append({
                "selected_id": ids[a], "control_id": ids[b],
                "selected_position": int(a), "control_position": int(b),
                "match_stratum": stratum,
            })
    columns = ["selected_id", "control_id", "selected_position", "control_position", "match_stratum"]
    pairs = pd.DataFrame.from_records(records, columns=columns)
    n_selected = int(selected.sum())
    n_pairs = len(pairs)
    reason = "ok" if n_pairs == n_selected and n_selected else (
        "no_selected_queries" if n_selected == 0 else
        "no_disjoint_matched_controls" if n_pairs == 0 else
        "partial_matching_control_shortage_or_missing_stratum"
    )
    return pairs, {
        "status": reason,
        "n_selected": n_selected,
        "n_pairs": n_pairs,
        "n_unmatched_selected": n_selected - n_pairs,
        "n_selected_missing_stratum": int((selected & ~present).sum()),
        "n_excluded_unselected_same_group": int((~selected & ~outside_selected_groups & present).sum()),
        "n_excluded_unselected_control_ineligible": int(
            (~selected & outside_selected_groups & present & ~control_allowed).sum()
        ),
        "control_eligible_column": control_eligible_col,
        "control_pool": "all_unselected" if control_eligible_col is None else "caller_frozen_eligible_unselected",
        "matching_fraction_of_selected": n_pairs / n_selected if n_selected else None,
        "outcome_blind": True,
        "replacement": False,
        "selected_groups_excluded": True,
        "strata": diagnostics,
    }


def individual_loss(y: np.ndarray, prediction: np.ndarray, loss_name: str) -> np.ndarray:
    """Return proper pointwise loss, with NaN for missing labels/predictions."""
    if loss_name not in {"mse", "brier", "logloss"}:
        raise ValueError(f"Unknown loss: {loss_name}")
    # One implementation for both head training and evaluation, including the
    # 1e-7 probability clipping convention for LogLoss.
    return row_loss(y, prediction, task="continuous" if loss_name == "mse" else "binary", loss_name=loss_name)


def prediction_metrics(y: np.ndarray, prediction: np.ndarray, task: str) -> dict[str, Any]:
    """Cohort-level metrics; an individual R2/AUC is deliberately not defined."""
    y = np.asarray(y, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    ok = np.isfinite(y) & np.isfinite(prediction)
    y, p = y[ok], prediction[ok]
    n = int(ok.sum())
    if task == "continuous":
        if not n:
            return {"n": 0, "mse": None, "rmse": None, "mae": None, "prediction_r2": None}
        error = y - p
        mse = float(np.mean(error ** 2))
        denominator = float(np.sum((y - y.mean()) ** 2))
        return {
            "n": n, "mse": mse, "rmse": float(np.sqrt(mse)),
            "mae": float(np.mean(np.abs(error))),
            "prediction_r2": float(1 - np.sum(error ** 2) / denominator) if denominator > 0 else None,
        }
    if task != "binary":
        raise ValueError(f"Unknown task: {task}")
    if not n:
        return {"n": 0, "cases": 0, "prevalence": None, "brier": None,
                "logloss": None, "auc": None, "auprc": None}
    has_both = len(np.unique(y)) == 2
    return {
        "n": n, "cases": int(np.sum(y)), "prevalence": float(np.mean(y)),
        "brier": float(np.mean(individual_loss(y, p, "brier"))),
        "logloss": float(np.mean(individual_loss(y, p, "logloss"))),
        "auc": float(roc_auc_score(y, p)) if has_both else None,
        "auprc": float(average_precision_score(y, p)) if has_both else None,
    }


def _effect(
    values_a: np.ndarray, mask_a: np.ndarray,
    values_b: np.ndarray | None = None, mask_b: np.ndarray | None = None,
) -> tuple[dict[str, Any], tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None]]:
    a = mask_a & np.isfinite(values_a)
    b = None if values_b is None else mask_b & np.isfinite(values_b)
    value = None
    if a.any() and (b is None or b.any()):
        value = float(np.mean(values_a[a]) - (np.mean(values_b[b]) if b is not None else 0))
    return ({
        "estimate": value, "ci_low": None, "ci_high": None,
        "n_left": int(a.sum()), "n_right": int(b.sum()) if b is not None else None,
        "bootstrap_valid_replicates": 0,
        "status": "ok" if value is not None else "not_estimable",
    }, (values_a, a, values_b, b))


def _bootstrap_effects(
    group_codes: np.ndarray,
    effects: list[tuple[dict[str, Any], tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None]]],
    n_boot: int, seed: int,
) -> None:
    """Resample original groups once per draw; reuse draws for every contrast.

    Aggregated loss sums/counts preserve pairing between base and candidate and
    dependence among relatives. Selected/control means are recomputed inside
    each draw, including the entire gain-enrichment difference.
    """
    n_groups = int(group_codes.max() + 1) if len(group_codes) else 0
    if n_boot == 0 or n_groups < 2:
        for result, _ in effects:
            if result["estimate"] is not None:
                result["ci_status"] = "bootstrap_disabled" if n_boot == 0 else "insufficient_groups"
        return
    columns: list[np.ndarray] = []
    mapping: list[tuple[int, int, int | None, int | None]] = []
    for _, (values_a, mask_a, values_b, mask_b) in effects:
        ia = len(columns)
        columns.extend([
            np.bincount(group_codes[mask_a], weights=values_a[mask_a], minlength=n_groups),
            np.bincount(group_codes[mask_a], minlength=n_groups),
        ])
        ib = None
        if values_b is not None:
            ib = len(columns)
            columns.extend([
                np.bincount(group_codes[mask_b], weights=values_b[mask_b], minlength=n_groups),
                np.bincount(group_codes[mask_b], minlength=n_groups),
            ])
        mapping.append((ia, ia + 1, ib, None if ib is None else ib + 1))
    aggregates = np.column_stack(columns)
    draws: list[list[float]] = [[] for _ in effects]
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        # Exactly n_groups sampled with replacement; no independent row/pair bootstrap.
        counts = np.bincount(rng.integers(0, n_groups, size=n_groups), minlength=n_groups)
        totals = np.einsum("g,gk->k", counts, aggregates, optimize=False)
        for j, (a_sum, a_n, b_sum, b_n) in enumerate(mapping):
            if totals[a_n] <= 0 or (b_n is not None and totals[b_n] <= 0):
                continue
            value = totals[a_sum] / totals[a_n]
            if b_sum is not None:
                value -= totals[b_sum] / totals[b_n]
            if np.isfinite(value):
                draws[j].append(float(value))
    for (result, _), samples in zip(effects, draws):
        result["bootstrap_valid_replicates"] = len(samples)
        if result["estimate"] is None:
            result["ci_status"] = "not_estimable"
        elif len(samples) >= max(2, int(np.ceil(0.9 * n_boot))):
            result["ci_low"], result["ci_high"] = [float(x) for x in np.quantile(samples, [0.025, 0.975])]
            result["ci_status"] = "ok"
        else:
            result["ci_status"] = "too_few_valid_group_bootstrap_replicates"


def evaluate_predictions(
    frame: pd.DataFrame,
    *,
    task: str,
    loss_name: str,
    selector: str,
    candidate_id: str,
    coverage: float | None = None,
    n_boot: int = 1000,
    seed: int = 20261020,
    min_matched_pairs: int = 20,
    min_matched_groups: int = 10,
    group_col: str = "group_id",
    selected_col: str = "selected",
    match_strata_col: str = "match_stratum",
    control_eligible_col: str | None = None,
) -> dict[str, Any]:
    """Evaluate one frozen selector/coverage and one fixed candidate.

    Required columns: id, group_id, y, p_base, p_candidate, supported, selected,
    match_stratum (the last three names except ``supported`` are configurable).
    ``coverage`` is the nominal development-frozen operating point, not a
    request to select query people again. Labels may be missing. Baseline
    predictions must cover every query; unavailable candidate predictions
    cause a recorded baseline fallback without changing the frozen selection.
    """
    if task not in {"continuous", "binary"}:
        raise ValueError("task must be continuous or binary")
    if (task == "continuous" and loss_name != "mse") or (
        task == "binary" and loss_name not in {"brier", "logloss"}
    ):
        raise ValueError("Use mse for continuous outcomes, brier/logloss for binary outcomes")
    if not isinstance(n_boot, (int, np.integer)) or n_boot < 0:
        raise ValueError("n_boot must be a nonnegative integer")
    if not isinstance(min_matched_pairs, (int, np.integer)) or min_matched_pairs < 2:
        raise ValueError("min_matched_pairs must be an integer of at least 2")
    if not isinstance(min_matched_groups, (int, np.integer)) or min_matched_groups < 2:
        raise ValueError("min_matched_groups must be an integer of at least 2")
    if coverage is not None and not 0 <= coverage <= 1:
        raise ValueError("coverage must be a fraction in [0,1]")
    ids, groups = _identity(frame, group_col)
    for col in ("y", "p_base", "p_candidate", "supported", selected_col, match_strata_col):
        if col not in frame:
            raise ValueError(f"Required column missing: {col}")
    selected = _bool_column(frame[selected_col], selected_col)
    supported = _bool_column(frame["supported"], "supported")
    if np.any(selected & ~supported):
        raise ValueError("Frozen selection must imply supported=True")
    y = frame["y"].to_numpy(dtype=float)
    base = frame["p_base"].to_numpy(dtype=float)
    candidate = frame["p_candidate"].to_numpy(dtype=float)
    if not np.isfinite(base).all():
        raise ValueError("Every query requires a finite baseline prediction")
    if np.isinf(y).any() or np.isinf(candidate).any():
        raise ValueError("Use NaN for missing labels/candidate predictions, not infinity")
    if task == "binary":
        if not np.isin(y[np.isfinite(y)], [0, 1]).all():
            raise ValueError("Observed binary outcomes must be 0 or 1")
        for name, p in (("p_base", base), ("p_candidate", candidate)):
            finite = p[np.isfinite(p)]
            if np.any((finite < 0) | (finite > 1)):
                raise ValueError(f"{name} must contain probabilities in [0,1]")

    # Matching occurs before accessing losses/label availability. Do not rematch later.
    pairs, matching = make_matched_controls(
        frame, seed=seed, group_col=group_col, selected_col=selected_col,
        match_strata_col=match_strata_col, control_eligible_col=control_eligible_col,
    )
    switch = selected & np.isfinite(candidate)
    policy = np.where(switch, candidate, base)
    observed = np.isfinite(y)
    candidate_observed = observed & np.isfinite(candidate)
    base_loss = individual_loss(y, base, loss_name)
    candidate_loss = individual_loss(y, candidate, loss_name)
    policy_loss = individual_loss(y, policy, loss_name)
    gain = base_loss - candidate_loss
    policy_gain = base_loss - policy_loss
    policy_gain_vs_candidate = candidate_loss - policy_loss
    all_mask = np.ones(len(frame), dtype=bool)
    pair_left = pairs["selected_position"].to_numpy(dtype=int)
    pair_right = pairs["control_position"].to_numpy(dtype=int)
    complete_base = observed[pair_left] & observed[pair_right]
    complete_gain = candidate_observed[pair_left] & candidate_observed[pair_right]
    base_left_mask = np.zeros(len(frame), dtype=bool)
    base_right_mask = np.zeros(len(frame), dtype=bool)
    gain_left_mask = np.zeros(len(frame), dtype=bool)
    gain_right_mask = np.zeros(len(frame), dtype=bool)
    base_left_mask[pair_left[complete_base]] = True
    base_right_mask[pair_right[complete_base]] = True
    gain_left_mask[pair_left[complete_gain]] = True
    gain_right_mask[pair_right[complete_gain]] = True

    effects = {
        "always_candidate_gain_vs_base": _effect(gain, all_mask),
        "policy_gain_vs_base": _effect(policy_gain, all_mask),
        "policy_gain_vs_candidate": _effect(policy_gain_vs_candidate, candidate_observed),
        "selected_gain_vs_base": _effect(gain, selected),
        "matched_selected_base_loss": _effect(base_loss, base_left_mask),
        "matched_control_base_loss": _effect(base_loss, base_right_mask),
        "predictability_loss_difference": _effect(base_loss, base_left_mask, base_loss, base_right_mask),
        "matched_selected_gain_vs_base": _effect(gain, gain_left_mask),
        "matched_control_gain_vs_base": _effect(gain, gain_right_mask),
        "gain_enrichment": _effect(gain, gain_left_mask, gain, gain_right_mask),
    }
    _, group_codes = np.unique(groups, return_inverse=True)
    for result, (_, left_mask, _, right_mask) in effects.values():
        result["n_unique_groups_left"] = int(len(np.unique(group_codes[left_mask])))
        result["n_unique_groups_right"] = (
            int(len(np.unique(group_codes[right_mask]))) if right_mask is not None else None
        )
    _bootstrap_effects(group_codes, list(effects.values()), n_boot, seed)
    result_effects = {name: value[0] for name, value in effects.items()}
    # A tiny matched subset may have a numerical mean but not adequate evidence.
    # Retain its descriptive estimate and explicitly withhold its interval.
    for names, n_complete, left_mask, right_mask in (
        (("predictability_loss_difference", "matched_selected_base_loss", "matched_control_base_loss"), int(complete_base.sum()), base_left_mask, base_right_mask),
        (("matched_selected_gain_vs_base", "matched_control_gain_vs_base", "gain_enrichment"),
         int(complete_gain.sum()), gain_left_mask, gain_right_mask),
    ):
        enough_groups = min(
            len(np.unique(group_codes[left_mask])), len(np.unique(group_codes[right_mask]))
        ) >= min_matched_groups
        if 0 < n_complete < min_matched_pairs:
            for name in names:
                result_effects[name].update(
                    status="insufficient_matched_pairs", ci_low=None, ci_high=None,
                    ci_status="insufficient_matched_pairs",
                )
        elif n_complete and not enough_groups:
            for name in names:
                result_effects[name].update(
                    status="insufficient_matched_groups", ci_low=None, ci_high=None,
                    ci_status="insufficient_matched_groups",
                )
    n_query, n_observed = len(frame), int(observed.sum())
    n_supported = int(supported.sum())
    selected_evaluable = int((selected & observed).sum())
    # A candidate's metrics on fewer rows must not masquerade as full-cohort coverage.
    always_candidate_status = "complete" if np.isfinite(candidate).all() else "partial_prediction_coverage"
    matching["n_pairs_with_both_outcomes"] = int(complete_base.sum())
    matching["n_pairs_with_both_outcomes_and_candidates"] = int(complete_gain.sum())
    matching["min_complete_pairs_for_inference"] = int(min_matched_pairs)
    matching["min_groups_per_side_for_inference"] = int(min_matched_groups)
    matching["n_unique_base_groups_left"] = int(len(np.unique(group_codes[base_left_mask])))
    matching["n_unique_base_groups_right"] = int(len(np.unique(group_codes[base_right_mask])))
    matching["n_unique_gain_groups_left"] = int(len(np.unique(group_codes[gain_left_mask])))
    matching["n_unique_gain_groups_right"] = int(len(np.unique(group_codes[gain_right_mask])))
    if len(pairs) and not complete_gain.any():
        matching["gain_status"] = "no_complete_outcome_candidate_pairs"
    else:
        matching["gain_status"] = (
            "not_estimable_no_pairs" if not complete_gain.any() else
            "insufficient_matched_pairs" if complete_gain.sum() < min_matched_pairs else
            "insufficient_matched_groups" if min(matching["n_unique_gain_groups_left"],
                                                matching["n_unique_gain_groups_right"]) < min_matched_groups else
            "estimable"
        )
    return {
        "schema_version": "abm.evaluation.1",
        "task": task, "loss": loss_name, "selector": selector,
        "candidate_id": candidate_id,
        "effect_signs": {
            "gain": "positive means candidate/policy reduces loss versus baseline",
            "policy_gain_vs_candidate": "positive means policy reduces loss versus always-candidate on paired candidate-evaluable queries",
            "gain_enrichment": "positive means more candidate gain in selected than matched controls",
            "predictability_loss_difference": "negative means lower baseline loss in selected than matched controls",
        },
        "coverage": {
            "nominal_development_coverage": coverage,
            "n_query": n_query,
            "n_evaluable": n_observed,
            "n_supported_query": n_supported,
            "n_supported_evaluable": int((supported & observed).sum()),
            "n_selected_query": int(selected.sum()),
            "n_selected_evaluable": selected_evaluable,
            "n_switched_query": int(switch.sum()),
            "n_switched_evaluable": int((switch & observed).sum()),
            "actual_selected_coverage_query": float(selected.mean()) if n_query else None,
            "actual_selected_coverage_evaluable": selected_evaluable / n_observed if n_observed else None,
            "actual_selected_coverage_supported": int(selected.sum()) / n_supported if n_supported else None,
            "actual_switch_coverage_query": float(switch.mean()) if n_query else None,
        },
        "whole_cohort": {
            "base": prediction_metrics(y, base, task),
            "candidate": prediction_metrics(y, candidate, task),
            "candidate_coverage_status": always_candidate_status,
            "policy": prediction_metrics(y, policy, task),
            "base_on_candidate_evaluable": prediction_metrics(y[candidate_observed], base[candidate_observed], task),
            "policy_on_candidate_evaluable": prediction_metrics(y[candidate_observed], policy[candidate_observed], task),
            "always_candidate_gain_vs_base": result_effects["always_candidate_gain_vs_base"],
            "policy_gain_vs_base": result_effects["policy_gain_vs_base"],
            "policy_gain_vs_candidate": result_effects["policy_gain_vs_candidate"],
        },
        "selected": {
            "base": prediction_metrics(y[selected], base[selected], task),
            "candidate": prediction_metrics(y[selected], candidate[selected], task),
            "gain_vs_base": result_effects["selected_gain_vs_base"],
        },
        "matched": {
            "selected_base_loss": result_effects["matched_selected_base_loss"],
            "control_base_loss": result_effects["matched_control_base_loss"],
            "predictability_loss_difference": result_effects["predictability_loss_difference"],
            "selected_gain_vs_base": result_effects["matched_selected_gain_vs_base"],
            "control_gain_vs_base": result_effects["matched_control_gain_vs_base"],
            "gain_enrichment": result_effects["gain_enrichment"],
            "selected_base": prediction_metrics(y[base_left_mask], base[base_left_mask], task),
            "control_base": prediction_metrics(y[base_right_mask], base[base_right_mask], task),
            "selected_candidate": prediction_metrics(y[gain_left_mask], candidate[gain_left_mask], task),
            "control_candidate": prediction_metrics(y[gain_right_mask], candidate[gain_right_mask], task),
        },
        "matching": matching,
        "match_pairs": pairs.to_dict(orient="records"),
        "policy": {
            "rule": "candidate where frozen selected and candidate available; baseline everywhere else",
            "n_selected_without_candidate": int((selected & ~np.isfinite(candidate)).sum()),
            "matching_changes_policy": False,
        },
        "bootstrap": {
            "unit": group_col,
            "n_original_groups": int(len(np.unique(groups))),
            "n_requested_replicates": int(n_boot), "seed": int(seed),
            "method": "percentile 95% CI; original groups sampled with replacement; matched assignments fixed",
            "conditional_on": "frozen predictor, selector, support rule, thresholds and matches",
            "does_not_include": "training, threshold estimation, match selection or GWAS discovery uncertainty",
        },
    }
