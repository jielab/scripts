#!/usr/bin/env python3
"""C1 support audits, existing-result evaluation and development training-budget diagnostics."""
from __future__ import annotations


# Frozen reference-support audit
import importlib
from collections import Counter

import numpy as np
import pandas as pd


COUNT_EDGES = (0, 1, 5, 10, 20, 40, 64, 128, 256, 512, 1024, np.inf)
PRIMARY_MODELS = ("elasticnet_full", "twin_clean", "twin_integrated_clean")
_QUANTILES = (("min", 0), ("p05", .05), ("p25", .25), ("median", .5),
              ("p75", .75), ("p95", .95), ("p99", .99), ("max", 1))
_DEVELOPMENT_SCOPES = {"calibration_fit", "calibration_audit", "development_calibration"}
_SEMANTICS = {
    "eligible_match_count": "uncapped eligible donors within the frozen clean-bank radius",
    "matched_count": "legacy name: used donors after K cap; never total eligible matches",
    "family_count": "distinct families among used donors",
    "family_ess": "effective number of families among used weighted donors",
    "donor_ess": "effective number of used weighted donors",
    "eligible_reference_count": "positive-weight clean-bank donors excluding query ID/family",
    "twin_fraction": "eligible_match_count / eligible_reference_count; same geometry/radius only",
    "matching_radius": "frozen clean-bank geometric radius; comparable only within the same feature geometry",
    "radius_ratio": "farthest used donor distance / frozen radius; legacy 100 means no donor",
    "nearest_radius_ratio": "nearest used donor distance / frozen radius; legacy 100 means no donor",
    "pre_prediction_confidence": "reference support score 0-100; not probability of correctness",
    "disease_risk": "frozen elasticnet_full horizon disease probability",
    "risk_resid_error": "frozen selector score; higher predicts smaller risk-residual squared error",
}


def _vector(value, n, name, *, dtype=float, scalar=False):
    x = np.asarray(value, dtype=dtype)
    if scalar and x.ndim == 0:
        x = np.full(n, x, dtype=dtype)
    if x.shape != (n,):
        raise ValueError(f"{name} must have one value per prediction row")
    return x


def _ids(value, name, *, unique=False):
    x = np.asarray(value)
    if x.ndim != 1 or pd.isna(x).any():
        raise ValueError(f"Invalid {name}")
    x = x.astype(str)
    if any(v.strip().lower() in {"", "nan", "none", "na", "<na>"} for v in x):
        raise ValueError(f"Missing {name}")
    if unique and len(set(x)) != len(x):
        raise ValueError(f"Duplicate {name}")
    return x


def _boolean(value, n, name):
    x = _vector(value, n, name, dtype=None)
    if pd.isna(x).any() or not np.isin(x, [False, True, 0, 1]).all():
        raise ValueError(f"{name} must be a complete Boolean mask")
    return x.astype(bool)


def _number(value, default=np.nan):
    return default if value is None else value


def _clean_bank_denominator(model, ids, groups):
    """Reconstruct the same eligibility exclusions as TwinBank.predict.

    A query's family is removed first; a same-ID donor in another family is
    removed separately.  No neighbour search or access to donor outcomes is
    necessary.  Missing bank fields mean unavailable, never a guessed size.
    """
    bank = getattr(model, "banks", {}).get("clean")
    if bank is None or any(not hasattr(bank, k) for k in ("w", "ids", "groups")):
        return np.full(len(ids), np.nan), "unavailable"
    bi = _ids(bank.ids, "clean-bank IDs", unique=True)
    bg = _ids(bank.groups, "clean-bank groups")
    bw = _vector(bank.w, len(bi), "clean-bank weights")
    if len(bg) != len(bi) or np.any(bw < 0) or not np.isfinite(bw).all():
        raise ValueError("Invalid clean-bank eligibility fields")
    positive = bw > 0
    family_n = Counter(bg[positive])
    id_group = dict(zip(bi[positive], bg[positive]))
    total = int(positive.sum())
    denominator = np.asarray([
        total - family_n[g] - int(i in id_group and id_group[i] != g)
        for i, g in zip(ids, groups)
    ], dtype=float)
    return denominator, "reconstructed_from_clean_bank_positive_weights_and_ID_family_exclusions"


def _support_fields(model, prediction, ids, groups):
    n = len(ids)
    diagnostic = prediction.get("diagnostics", {})
    qc = _boolean(diagnostic["technical_qc"], n, "technical_qc")
    fields, sources = {}, {}
    for name in _SEMANTICS:
        if name in {"disease_risk", "risk_resid_error", "twin_fraction", "eligible_reference_count"}:
            continue
        fields[name] = (_vector(diagnostic[name], n, name).copy() if name in diagnostic
                        else np.full(n, np.nan))
        sources[name] = "prediction.diagnostics" if name in diagnostic else "unavailable"
    if "matching_radius" not in diagnostic:
        radius = getattr(getattr(model, "banks", {}).get("clean"), "radius", None)
        if radius is not None and np.isfinite(radius) and radius > 0:
            fields["matching_radius"] = np.full(n, float(radius))
            sources["matching_radius"] = "frozen_clean_bank_radius"
    # TwinABM._raw queries its banks only for QC-passing rows. Zeros placed in
    # other rows are placeholders, not observed zero neighbours or zero ESS.
    for name in fields:
        if name != "pre_prediction_confidence":
            fields[name][~qc] = np.nan
    for name in ("eligible_match_count", "matched_count", "family_count"):
        x = fields[name]
        known = np.isfinite(x)
        if np.isinf(x).any() or np.any(x[known] < 0) or np.any(x[known] != np.floor(x[known])):
            raise ValueError(f"{name} must contain nonnegative integer counts or NaN")
    for name in ("family_ess", "donor_ess"):
        x = fields[name]
        if np.isinf(x).any() or np.any(x[np.isfinite(x)] < 0):
            raise ValueError(f"{name} must be nonnegative or NaN")
    confidence = fields["pre_prediction_confidence"]
    if (np.isinf(confidence).any() or
            np.any((confidence[np.isfinite(confidence)] < 0) | (confidence[np.isfinite(confidence)] > 100))):
        raise ValueError("pre_prediction_confidence must be in [0, 100] or NaN")
    raw, used = fields["eligible_match_count"], fields["matched_count"]
    both = np.isfinite(raw) & np.isfinite(used)
    if np.any(used[both] > raw[both]):
        raise ValueError("Used matched_count cannot exceed uncapped eligible_match_count")
    k = float(_number(getattr(model.config, "k", None)))
    if np.isfinite(k) and (k < 1 or int(k) != k or np.any(used[np.isfinite(used)] > k)):
        raise ValueError("Invalid K or used matched_count exceeds K")
    if "eligible_reference_count" in diagnostic:
        denominator = _vector(diagnostic["eligible_reference_count"], n,
                              "eligible_reference_count", scalar=True).copy()
        denominator_source = "prediction.diagnostics"
    else:
        denominator, denominator_source = _clean_bank_denominator(model, ids, groups)
    denominator[~qc] = np.nan
    finite = np.isfinite(denominator)
    if (np.isinf(denominator).any() or np.any(denominator[finite] < 0)
            or np.any(denominator[finite] != np.floor(denominator[finite]))):
        raise ValueError("eligible_reference_count must be a nonnegative integer or NaN")
    if np.any(raw[finite & np.isfinite(raw)] > denominator[finite & np.isfinite(raw)]):
        raise ValueError("Raw match count exceeds the eligible reference denominator")
    fraction = np.full(n, np.nan)
    positive = np.isfinite(raw) & finite & (denominator > 0)
    fraction[positive] = raw[positive] / denominator[positive]
    fraction_source = ("derived_from_raw_count_and_" + denominator_source
                       if positive.any() else "unavailable")
    if "twin_fraction" in diagnostic:
        supplied = _vector(diagnostic["twin_fraction"], n, "twin_fraction").copy()
        supplied[~qc] = np.nan
        finite_fraction = np.isfinite(supplied)
        if np.isinf(supplied).any() or np.any((supplied[finite_fraction] < 0) | (supplied[finite_fraction] > 1)):
            raise ValueError("twin_fraction must be in [0, 1] or NaN")
        comparable = positive & finite_fraction
        if not np.allclose(supplied[comparable], fraction[comparable], rtol=1e-7, atol=1e-10):
            raise ValueError("twin_fraction disagrees with raw_count / eligible_reference_count")
        # A zero denominator is not an observed zero twin fraction.
        supplied[finite & (denominator == 0)] = np.nan
        fields["twin_fraction"] = supplied
        sources["twin_fraction"] = "prediction.diagnostics"
    else:
        fields["twin_fraction"], sources["twin_fraction"] = fraction, fraction_source
    fields["eligible_reference_count"] = denominator
    sources["eligible_reference_count"] = denominator_source
    # Radius ratios use a historical sentinel when no donor is used. Exclude
    # that sentinel from radius quantiles; the zero-count stratum remains.
    for name in ("radius_ratio", "nearest_radius_ratio"):
        fields[name][np.isfinite(used) & (used == 0)] = np.nan
    if "supported" in diagnostic:
        supported = _boolean(diagnostic["supported"], n, "supported").astype(float)
        support_source = "prediction.diagnostics"
    elif all(hasattr(model.config, key) for key in ("min_twins", "min_family_ess")):
        known = np.isfinite(used) & np.isfinite(fields["family_ess"])
        supported = np.full(n, np.nan)
        supported[known] = ((used[known] >= model.config.min_twins)
                            & (fields["family_ess"][known] >= model.config.min_family_ess))
        support_source = "reconstructed_from_used_count_family_ess_and_frozen_minima"
    else:
        supported, support_source = np.full(n, np.nan), "unavailable"
    supported[~qc] = np.nan
    return fields, sources, supported, support_source, k


def _summary(values, mask):
    x = np.asarray(values)[mask]
    x = x[np.isfinite(x)]
    out = {"finite_N": len(x), "missing_N": int(mask.sum()) - len(x),
           "mean": float(x.mean()) if len(x) else np.nan}
    out.update({name: float(np.quantile(x, q)) if len(x) else np.nan for name, q in _QUANTILES})
    return out


def _fraction(numerator, denominator):
    return float(numerator / denominator) if denominator else np.nan


def _metric_row(abm, model, time, event, groups, probability, mask, supported):
    out = abm.metrics(time[mask], event[mask], probability[mask], model.km)
    out.update(event_N=int(event[mask].sum()), family_N=len(set(groups[mask])),
               supported_N=int(np.sum(supported[mask] == 1)),
               support_known_N=int(np.isfinite(supported[mask]).sum()),
               horizon=float(model.km.horizon),
               case_control_pairs=int(out["cases"]) * int(out["controls"]),
               Uno_IPCW_comparable_pair_weight=out.get("comparable_pairs", np.nan),
               same_people_within_stratum=True)
    minimum_cases = int(getattr(model.config, "audit_min_cases", 30))
    minimum_controls = int(getattr(model.config, "audit_min_controls", 30))
    out["AUC_status"] = (
        "empty_subset" if not out["N"] else
        "insufficient_case_control_information" if not np.isfinite(out["AUC_IPCW"]) else
        "low_case_control_information" if out["cases"] < minimum_cases or out["controls"] < minimum_controls else
        "evaluated"
    )
    out["audit_min_cases"], out["audit_min_controls"] = minimum_cases, minimum_controls
    return out


def _reference(reference, prediction, ids, groups, grouping_unit):
    if reference is None:
        return None
    if reference.get("evaluation_scope") not in _DEVELOPMENT_SCOPES:
        raise ValueError("reference_prediction must explicitly identify a development calibration scope")
    ri = _ids(reference["eid"], "reference IDs", unique=True)
    if len(set(ri) & set(ids)):
        raise ValueError("Reference calibration and evaluation participants must be disjoint")
    if grouping_unit == "family":
        rg = _ids(reference["group"], "reference family groups")
        if len(rg) != len(ri) or len(set(rg) & set(groups)):
            raise ValueError("Reference calibration and evaluation families must be disjoint")
    ph, rh = prediction.get("frozen_model_hash"), reference.get("frozen_model_hash")
    if ph is not None and rh is not None and ph != rh:
        raise ValueError("Reference and evaluation predictions use different frozen models")
    return ri


def _ordering_table(abm, model, prediction, reference, reference_ids,
                    time, event, groups, probabilities, common, supported, fields):
    """Ordering diagnostics, not a fitted calibration model or individual guarantee."""
    if "elasticnet_full" not in probabilities:
        return pd.DataFrame([{"status": "unavailable_target_model", "target_model": "elasticnet_full"}])
    entries = [("pre_prediction_confidence", fields["pre_prediction_confidence"],
                "diagnostics", "higher_means_more_reference_support_not_correctness_probability")]
    for key in ("risk_resid_error", "absolute_error", "support_quality"):
        if key in prediction.get("selector_scores", {}):
            entries.append((key, _vector(prediction["selector_scores"][key], len(time), key),
                            "selector_scores", "higher_means_more_predicted_predictability"))
    for key in ("absolute_error", "risk_residual_squared_error", "risk_residual_logloss"):
        if key in prediction.get("estimated_error", {}):
            entries.append((key, _vector(prediction["estimated_error"][key], len(time), key),
                            "estimated_error", "lower_means_smaller_predicted_error"))
    rows = []
    for name, score, source, direction in entries:
        # The original absolute-error learner targets config.primary, whereas
        # the newer risk-residual extension targets elasticnet_full. Comparing
        # the former to full-model Brier would silently change its estimand.
        target = (getattr(model.config, "primary", None) if name == "absolute_error"
                  else "elasticnet_full")
        if target not in probabilities:
            rows.append(dict(score=name, score_source=source, status="unavailable_target_model",
                             target_model=target, interpretation=direction))
            continue
        reference_score = None
        if reference is not None and name in reference.get(source, {}):
            reference_score = _vector(reference[source][name], len(reference_ids), f"reference {name}")
            reference_qc = _boolean(reference["diagnostics"]["technical_qc"], len(reference_ids), "reference technical_qc")
            reference_score = reference_score[np.isfinite(reference_score) & reference_qc]
        if reference_score is not None and len(reference_score):
            internal = np.unique(np.quantile(reference_score, np.arange(1, 10) / 10))
            edges = np.r_[-np.inf, internal, np.inf]
            bin_source = "independent_development_calibration_quantiles"
        elif name == "pre_prediction_confidence":
            edges = np.r_[-np.inf, np.arange(10, 100, 10), np.inf]
            bin_source = "fixed_10_point_support_score_intervals"
        else:
            rows.append(dict(score=name, score_source=source, status="unavailable_development_reference",
                             target_model=target, interpretation=direction,
                             cutpoints_source="none_test_quantiles_forbidden"))
            continue
        for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            mask = common & np.isfinite(score) & (score >= lo) & (score < hi)
            row = _metric_row(abm, model, time, event, groups, probabilities[target], mask, supported)
            row.update(score=name, score_source=source, bin=index, lower=lo, upper=hi,
                       interval="[lower, upper)", target_model=target,
                       cutpoints_source=bin_source, mean_score=float(score[mask].mean()) if mask.any() else np.nan,
                       reference_scope=reference.get("evaluation_scope") if reference is not None else "unavailable",
                       interpretation=direction, evaluation_only=True,
                       absolute_error_calibration=(source == "estimated_error" and name == "absolute_error"))
            # Only predicted absolute squared error is on the Brier scale.
            # A residual error or a support score must not be called calibrated
            # probability of correctness, or compared numerically to Brier.
            if row["absolute_error_calibration"]:
                _, weights = model.km.labels_weights(time[mask], event[mask])
                valid = weights > 0
                row["IPCW_mean_predicted_squared_error"] = (
                    float(np.average(score[mask][valid], weights=weights[valid])) if valid.any() else np.nan)
                row["observed_minus_predicted_squared_error"] = (
                    row["Brier_IPCW"] - row["IPCW_mean_predicted_squared_error"])
            rows.append(row)
    return pd.DataFrame(rows)


def _count_gate_tables(abm, model, prediction, reference, reference_ids, time, event,
                       ids, groups, probabilities, common, supported, raw, bootstrap):
    """Predeclared count gates, with outcome-blind risk-composition controls."""
    if "elasticnet_full" not in probabilities:
        missing = pd.DataFrame([{"status": "unavailable_target_model"}])
        return missing, missing.copy(), missing.copy()
    selector_module = importlib.import_module("c1_evaluation")
    risk = probabilities["elasticnet_full"]
    pool = common & np.isfinite(raw)
    rows = []
    seed = int(getattr(model.config, "seed", 2026))
    for threshold in COUNT_EDGES[1:-1]:
        selected = pool & (raw >= threshold)
        risk_edges = getattr(model, "risk_edges", None)
        if risk_edges is not None:
            edges = np.asarray(risk_edges, float)
            if edges.ndim != 1 or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
                raise ValueError("Frozen model risk edges must be finite and increasing")
            matched = abm.risk_matched_mask(selected, pool, risk, ids, edges, seed + 7241 + int(threshold))
            matching_source = "frozen_development_risk_edges"
        else:
            matched = selector_module._risk_match(selected, pool, risk, ids, f"support_count:{seed}:{threshold:g}")
            matching_source = "query_risk_deciles_X_only_descriptive_control"
        a = _metric_row(abm, model, time, event, groups, risk, selected, supported)
        b = _metric_row(abm, model, time, event, groups, risk, matched, supported)
        row = dict(threshold=threshold, rule="uncapped_eligible_match_count >= threshold",
                   threshold_source="fixed_predeclared_count_grid", target_model="elasticnet_full",
                   common_pool_N=int(pool.sum()), selected_N=a["N"], riskmatched_N=b["N"],
                   coverage=_fraction(a["N"], len(ids)), coverage_denominator="all_query_rows",
                   coverage_of_common_evaluation=_fraction(a["N"], int(common.sum())),
                   coverage_of_count_observed_pool=_fraction(a["N"], int(pool.sum())),
                   selected_event_N=a["event_N"], riskmatched_event_N=b["event_N"],
                   selected_cases=a["cases"], riskmatched_cases=b["cases"],
                   selected_family_N=a["family_N"], riskmatched_family_N=b["family_N"],
                   selected_mean_risk=float(risk[selected].mean()) if selected.any() else np.nan,
                   riskmatched_mean_risk=float(risk[matched].mean()) if matched.any() else np.nan,
                   overlap_N=int((selected & matched).sum()), matching_source=matching_source,
                   count_is_disease_risk=False, status=a["status"], AUC_status=a["AUC_status"],
                   uncertainty="descriptive_point_estimates_no_CI_or_p_value",
                   control_interpretation="equal counts within risk bins; overlap allowed; no exact continuous-risk balance guarantee")
        for metric in ("Brier_IPCW", "LogLoss_IPCW", "Harrell_C", "AUC_IPCW"):
            row[metric], row["riskmatched_" + metric] = a[metric], b[metric]
            row["delta_" + metric] = a[metric] - b[metric]
        rows.append(row)
    curve = pd.DataFrame(rows)
    unavailable = {"route": "raw_count", "status": "not_requested" if not bootstrap else "unavailable_development_reference",
                   "requested_boot": bootstrap, "interpretation": "raw count alone is not the deployed error selector"}
    if not bootstrap or reference is None:
        status = pd.DataFrame([unavailable])
        return curve, status, status.copy()
    from c1_protocol import COVERAGES, PRIMARY_SELECTOR
    if ("eligible_match_count" not in reference.get("diagnostics", {}) or
            any(f"{PRIMARY_SELECTOR}:{c:g}" not in prediction.get("masks", {}) for c in COVERAGES)):
        unavailable["status"] = "unavailable_reference_count_or_frozen_primary_masks"
        status = pd.DataFrame([unavailable])
        return curve, status, status.copy()
    nr = len(reference_ids)
    ref_raw = _vector(reference["diagnostics"]["eligible_match_count"], nr, "reference raw count")
    ref_qc = _boolean(reference["diagnostics"]["technical_qc"], nr, "reference technical_qc")
    ref_eligible = ref_qc & np.isfinite(ref_raw)
    if not ref_eligible.any():
        unavailable["status"] = "unavailable_finite_reference_count"
        status = pd.DataFrame([unavailable])
        return curve, status, status.copy()
    if np.any(ref_raw[ref_eligible] < 0) or np.any(ref_raw[ref_eligible] != np.floor(ref_raw[ref_eligible])):
        raise ValueError("Reference raw count must contain nonnegative integers")
    tie_seed = seed + 7242
    ref_ties = np.array([abm.tie_hash(v, tie_seed) for v in reference_ids], np.uint64)
    query_ties = np.array([abm.tie_hash(v, tie_seed) for v in ids], np.uint64)
    eligible_idx = np.flatnonzero(ref_eligible)
    order = eligible_idx[np.lexsort((ref_ties[eligible_idx], -ref_raw[eligible_idx]))]
    masks = {key: _boolean(mask, len(ids), key).copy() & pool
             for key, mask in prediction.get("masks", {}).items()}
    thresholds = {}
    for cov in COVERAGES:
        wanted = min(len(order), int(np.ceil(cov * nr)))
        boundary = order[wanted - 1]
        threshold = float(ref_raw[boundary]) if wanted < len(order) else -np.inf
        tie = int(ref_ties[boundary])
        masks[f"raw_count:{cov:g}"] = pool & ((raw > threshold) | ((raw == threshold) & (query_ties <= tie)))
        thresholds[float(cov)] = threshold
    # The existing locked inference compares every gate using elasticnet_full,
    # joint family resampling, and the same risk-matched-control construction.
    # Use a shallow copy and replacement dictionaries; never alter predictions.
    evaluation_prediction = dict(prediction, masks=masks,
                                 diagnostics=dict(prediction["diagnostics"], technical_qc=pool))
    calibrated_curve, inference = selector_module.infer_selector(
        evaluation_prediction, time, event, groups, model.km,
        routes=(PRIMARY_SELECTOR, "raw_count"), B=bootstrap, seed=seed,
        grouping_unit=getattr(model.config, "grouping_unit", "individual_unverified"))
    calibrated_curve = calibrated_curve.loc[calibrated_curve.route == "raw_count"].copy()
    inference = inference.loc[inference.route == "raw_count"].copy()
    calibrated_curve = calibrated_curve.drop(columns=[c for c in calibrated_curve if "eid" in c.lower()])
    calibrated_curve["raw_count_threshold"] = calibrated_curve.coverage_target.map(thresholds)
    for table in (calibrated_curve, inference):
        table["threshold_source"] = "independent_development_calibration_raw_count_and_fixed_hash_ties"
        table["reference_scope"] = reference["evaluation_scope"]
        table["reference_N"] = nr
        table["reference_eligible_N"] = int(ref_eligible.sum())
        table["target_model"] = "elasticnet_full"
        table["same_common_population"] = True
        table["interpretation"] = "exploratory_count_only_gate_not_the_deployed_error_selector"
    return curve, calibrated_curve, inference


def build_support_audit(model, prediction, time, event, truth_ids, *,
                        reference_prediction=None, bootstrap=None) -> dict[str, pd.DataFrame]:
    """Evaluate frozen support on aggregate tables with one common model roster.

    Outcomes/IDs must exactly align with ``prediction``. Probability NaNs in
    *any* input model are excluded jointly from every error comparison; support
    distributions retain all query rows. The three primary models are reported
    when present. No raw or per-participant output is returned.

    ``reference_prediction`` must be an explicitly labelled, disjoint development
    calibration prediction. It supplies optional ordering/gating boundaries,
    never test outcomes. With no reference, fixed count/support-score bins are
    still evaluated; other ordering quantiles remain unavailable.

    ``bootstrap=None`` (or 0) skips the additional count-only selector bootstrap.
    A positive integer requests the existing locked family/group inference when
    both reference counts and the six frozen primary masks are available. Its
    intervals are conditional on frozen models, gates, controls and censoring;
    they do not include training or model-selection uncertainty.
    """
    abm = importlib.import_module("abm")  # dynamic: abm may call this module
    ids = _ids(prediction["eid"], "prediction IDs", unique=True)
    groups = _ids(prediction["group"], "prediction groups")
    n = len(ids)
    if len(groups) != n or not np.array_equal(ids, _ids(truth_ids, "truth IDs", unique=True)):
        raise ValueError("Truth IDs/order and groups must match the frozen prediction roster")
    time = _vector(time, n, "time")
    event = _vector(event, n, "event")
    if not np.isfinite(time).all() or np.any(time <= 0) or not np.isin(event, [0, 1]).all():
        raise ValueError("Invalid evaluation outcomes")
    event = event.astype(int)
    frozen_hash = prediction.get("frozen_model_hash")
    if frozen_hash is not None and frozen_hash != abm.model_hash(model):
        raise ValueError("Predictions do not match the supplied frozen model")
    if bootstrap is None:
        bootstrap = 0
    if isinstance(bootstrap, bool) or not np.isfinite(bootstrap) or int(bootstrap) != bootstrap or bootstrap < 0:
        raise ValueError("bootstrap must be a nonnegative integer or None")
    bootstrap = int(bootstrap)
    diagnostic = prediction.get("diagnostics", {})
    qc = _boolean(diagnostic["technical_qc"], n, "technical_qc")
    selected = _boolean(prediction["selected"], n, "selected") & qc
    probabilities = {name: _vector(values, n, f"probability {name}")
                     for name, values in prediction.get("probability", {}).items()}
    if not probabilities:
        raise ValueError("At least one frozen probability model is required")
    for name, values in probabilities.items():
        finite = np.isfinite(values)
        if np.any((values[finite] < 0) | (values[finite] > 1)):
            raise ValueError(f"{name} contains scores outside [0, 1], not absolute probabilities")
    common_finite = np.logical_and.reduce([np.isfinite(values) for values in probabilities.values()])
    common = qc & common_finite
    models = [name for name in PRIMARY_MODELS if name in probabilities] or list(probabilities)
    fields, sources, supported, support_source, k = _support_fields(model, prediction, ids, groups)
    fields["disease_risk"] = probabilities.get("elasticnet_full", np.full(n, np.nan))
    sources["disease_risk"] = "probability.elasticnet_full" if "elasticnet_full" in probabilities else "unavailable"
    fields["risk_resid_error"] = _vector(prediction.get("selector_scores", {}).get(
        "risk_resid_error", np.full(n, np.nan)), n, "risk_resid_error")
    sources["risk_resid_error"] = "selector_scores" if "risk_resid_error" in prediction.get("selector_scores", {}) else "unavailable"
    reference_ids = _reference(reference_prediction, prediction, ids, groups,
                               getattr(model.config, "grouping_unit", "individual_unverified"))
    subsets = {"all": np.ones(n, bool), "selected": selected, "rejected": ~selected}
    distribution = []
    raw, used = fields["eligible_match_count"], fields["matched_count"]
    for subset, mask in subsets.items():
        finite_used = mask & np.isfinite(used)
        finite_raw = mask & np.isfinite(raw)
        for name, values in fields.items():
            row = dict(subset=subset, diagnostic=name, N=int(mask.sum()),
                       technical_qc_N=int((mask & qc).sum()), common_evaluation_N=int((mask & common).sum()),
                       family_N=len(set(groups[mask])), source=sources[name], semantics=_SEMANTICS[name],
                       K=k, supported_N=int(np.sum(supported[mask] == 1)),
                       support_known_N=int(np.isfinite(supported[mask]).sum()), support_source=support_source,
                       used_K_saturated_N=int(np.sum(finite_used & (used == k))) if np.isfinite(k) else np.nan,
                       used_K_saturated_fraction=_fraction(int(np.sum(finite_used & (used == k))), int(finite_used.sum())) if np.isfinite(k) else np.nan,
                       raw_count_exceeds_K_N=int(np.sum(finite_raw & (raw > k))) if np.isfinite(k) else np.nan,
                       raw_count_exceeds_K_fraction=_fraction(int(np.sum(finite_raw & (raw > k))), int(finite_raw.sum())) if np.isfinite(k) else np.nan,
                       raw_zero_N=int(np.sum(finite_raw & (raw == 0))),
                       distribution_scope="all_query_rows_with_unevaluated_QC_failure_support_marked_missing",
                       fraction_comparability="only_with_same_feature_geometry_and_frozen_radius")
            row.update(_summary(values, mask))
            row["status"] = "evaluated" if row["finite_N"] else "empty_subset" if not mask.any() else "unavailable"
            distribution.append(row)
    error_rows = []
    count_masks = [(lo, hi, f"[{lo:g},{hi:g})", np.isfinite(raw) & (raw >= lo) & (raw < hi))
                   for lo, hi in zip(COUNT_EDGES[:-1], COUNT_EDGES[1:])]
    count_masks.append((np.nan, np.nan, "missing_raw_count", ~np.isfinite(raw)))
    for subset, subset_mask in subsets.items():
        for lower, upper, label, count_mask in count_masks:
            requested = subset_mask & qc & count_mask
            mask = requested & common_finite
            for name in models:
                row = _metric_row(abm, model, time, event, groups, probabilities[name], mask, supported)
                row.update(subset=subset, model=name, count_bin=label, count_lower=lower, count_upper=upper,
                           count_source="eligible_match_count_uncapped", bin_source="fixed_predeclared_count_grid",
                           requested_N=int(requested.sum()), excluded_nonfinite_N=int((requested & ~common_finite).sum()),
                           common_model_count=len(probabilities),
                           support_diagnostic_bank="clean_molecular_bank",
                           prediction_semantics="frozen_risk_including_fallback_where_used",
                           grouping_unit=getattr(model.config, "grouping_unit", "individual_unverified"),
                           diagnostic_status="missing_raw_count" if label == "missing_raw_count" else "available")
                error_rows.append(row)
    ordering = _ordering_table(abm, model, prediction, reference_prediction, reference_ids,
                               time, event, groups, probabilities, common, supported, fields)
    gate, gate_curve, gate_inference = _count_gate_tables(
        abm, model, prediction, reference_prediction, reference_ids, time, event, ids, groups,
        probabilities, common, supported, raw, bootstrap)
    status = pd.DataFrame([dict(
        query_N=n, technical_qc_N=int(qc.sum()), common_evaluation_N=int(common.sum()),
        excluded_nonfinite_N=int((qc & ~common_finite).sum()),
        common_model_count=len(probabilities), common_models=";".join(probabilities),
        reported_models=";".join(models), missing_primary_models=";".join(sorted(set(PRIMARY_MODELS) - set(probabilities))),
        raw_count_source=sources["eligible_match_count"], used_count_source=sources["matched_count"],
        denominator_source=sources["eligible_reference_count"], twin_fraction_source=sources["twin_fraction"],
        frozen_clean_bank_radius=_number(getattr(getattr(model, "banks", {}).get("clean"), "radius", None)),
        model_prediction_hash_status="verified" if frozen_hash is not None else "unavailable_in_prediction",
        reference_scope=reference_prediction.get("evaluation_scope") if reference_prediction is not None else "unavailable",
        bootstrap_requested=bootstrap, common_population_rule="technical_QC_and_finite_for_every_input_probability_model",
        same_people_within_stratum=True, test_fitting=False, contains_participant_identifiers=False,
        inference_scope="conditional_frozen_model_audit", support_is_correctness_probability=False,
        status="completed" if common.any() else "no_common_evaluation_population")])
    return {
        "support_distribution": pd.DataFrame(distribution),
        "support_error_by_count": pd.DataFrame(error_rows),
        "support_error_ordering": ordering,
        "support_count_gate": gate,
        "support_count_gate_curve": gate_curve,
        "support_count_gate_inference": gate_inference,
        "support_audit_status": status,
    }


# Existing-result audit
import argparse
from dataclasses import replace
import gzip
import json
from pathlib import Path
import re

import joblib
import numpy as np
import pandas as pd

import abm
from le8_common import sha256, write_tables

REPORT_VERSION = "2026-10-21.fair-evaluation-v1"
REQUIRED = (
    "fit/model_LOCAL_ONLY.joblib", "fit/MODEL_FROZEN.json",
    "test_predictions_LOCAL_ONLY.joblib", "test_truth_LOCAL_ONLY.npz",
    "host_contract.json",
)
INTERNAL_AGGREGATES = ("metrics", "contrasts", "prediction_finiteness", "calibration",
                       "confidence_validation", "strata_diagnostics")


def _read_json(path):
    return json.loads(Path(path).read_text())


def _aggregate_tables(tables):
    """Reject participant-level columns/cells before creating any report files."""
    checked = {}
    for name, values in tables.items():
        frame = pd.DataFrame(values)
        for column in frame.columns:
            label = str(column).lower()
            aggregate_hash = label.endswith(("_hash", "_sha256"))
            if not aggregate_hash and (label in {"id", "ids", "group", "groups", "family"} or re.search(
                    r"(?:^|_)(?:eids?|participant_ids?|sample_ids?|subject_ids?|family_ids?|group_ids?)(?:$|_)", label)):
                raise ValueError(f"Participant identifiers cannot be published: {name}/{column}")
            if frame[column].map(lambda x: isinstance(x, (dict, list, tuple, set, np.ndarray))).any():
                raise ValueError(f"Only scalar aggregate cells can be published: {name}/{column}")
        checked[name] = frame
    return checked


def _record_sources(sources, additions):
    """Never replace an earlier checksum when a component reloads a source."""
    for filename, digest in additions.items():
        key = str(Path(filename).resolve())
        if key in sources and sources[key] != digest:
            raise ValueError(f"Source changed between report components: {filename}")
        sources[key] = digest


def _publish(out, tables, metadata):
    """Only publish to a fresh directory, with a final completion receipt."""
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"Choose a new report directory; will not replace: {out}")
    tables = _aggregate_tables(tables)
    for filename, digest in metadata.get("sources", {}).items():
        if sha256(Path(filename)) != digest:
            raise ValueError(f"Source changed during reporting: {filename}")
    out.mkdir(parents=True)
    write_tables(out, tables, "Frozen fair evaluation", metadata["layer"],
                 metadata["trait"], workbook="fair_evaluation.xlsx")
    abm.json_write(out / "report.json", metadata)
    (out / "README.txt").write_text(
        "This report reuses frozen predictions; no model was fitted.\n"
        "Harrell C uses strict ordered follow-up times; IPCW AUC uses the frozen horizon and KM.\n"
        "Compare models within a table/subset with identical N and roster hash.\n"
        "All versus selected changes the evaluated population, and is not an algorithm comparison.\n"
        "Risk-matched controls have identical counts within risk bins, not exact continuous-risk balance.\n"
        "Support count is not disease risk or a probability of correctness.\n"
        "Different reference sizes, geometries and radii must be reported when comparing counts.\n"
        "Confidence intervals condition on the frozen fit, selectors and censoring model.\n"
        "Unavailable components are not validated findings. No participant IDs are exported.\n",
        encoding="utf-8",
    )
    for filename, digest in metadata.get("sources", {}).items():
        if sha256(Path(filename)) != digest:
            raise ValueError(f"Source changed during reporting: {filename}")
    files = {p.name: sha256(p) for p in out.iterdir() if p.is_file()}
    abm.json_write(out / "COMPLETE.json", dict(status="completed", report_status=metadata["status"], files=files))
    return out


def _reference_prediction(model, development):
    """Project only the allowed audit role; do not relax model overlap checks."""
    if not development.is_file():
        return None, "unavailable_development_file"
    with np.load(development, allow_pickle=False) as data:
        required = {"role", "eid", "group", "X", "feature_names"}
        if not required <= set(data.files):
            raise ValueError("Development file lacks audit projection fields")
        ix = np.asarray(data["role"], str) == "calibration_audit"
        if not ix.any():
            return None, "unavailable_calibration_audit_role"
        # time/event are deliberately never accessed by this function.
        X, ids, groups, names = data["X"][ix], data["eid"][ix], data["group"][ix], data["feature_names"]
    before = abm.model_hash(model)
    pred = model.predict(X, ids, groups, names, allow_development_audit=True)
    pred = abm.apply_predictability_extension(model, pred, X, ids, groups)
    if (pred.get("evaluation_scope") != "calibration_audit" or abm.model_hash(model) != before
            or pred.get("frozen_model_hash") != before):
        raise ValueError("Audit projection changed the frozen model or used an invalid role")
    if (not np.array_equal(np.asarray(pred["eid"], str), np.asarray(ids, str)) or
            not np.array_equal(np.asarray(pred["group"], str), np.asarray(groups, str))):
        raise ValueError("Audit projection changed the requested participant/family roster")
    return pred, "projected_frozen_calibration_audit_without_outcomes"


def _training_tables(model, contract):
    audit = getattr(model, "training_audit", {})
    roles = [dict(role=k, **v) for k, v in audit.get("roles", {}).items()]
    composition = audit.get("pruning_composition")
    if not composition:
        known_cases, retained_cases = audit.get("known_cases"), audit.get("retained_cases")
        composition = dict(status="legacy_known_case_summary_eligibility_denominator_not_reconstructed",
                           known_build_cases=known_cases, retained_cases=retained_cases,
                           eligible_reference_N=audit.get("known_donors"), retained_reference_N=audit.get("retained_donors"),
                           retention_relative_to_known_build_cases=(retained_cases / known_cases
                               if known_cases and retained_cases is not None else np.nan),
                           interpretation="Known build cases may include QC-ineligible cases; this ratio is not attributed exclusively to pruning")
    return {
        "training_roles": pd.DataFrame(roles),
        "training_pruning": pd.DataFrame(audit.get("pruning", [])),
        "training_case_retention": pd.DataFrame([composition]),
        "training_comparability": pd.DataFrame([dict(
            status="information_and_training_budget_not_proven_matched_to_C1",
            same_heldout_people_does_not_imply_same_training_budget=True,
            ABM_molecular_budget=model.config.panel_size,
            ABM_selected_panel_names=json.dumps(getattr(model, "panel_names", [])),
            ABM_selected_molecular_names=json.dumps(audit.get("molecular_panel_features", [])),
            clinical_feature_names=json.dumps(contract.get("clinical_feature_names", [])),
            reference_experiment=getattr(model.config, "reference_experiment", "legacy"),
            radius_policy=getattr(model.config, "reference_radius_policy", "arm_calibrated"),
            twin_random_control=getattr(model.config, "twin_random_control", False),
            action="Use c1_audit.py budget for ABM internal common-input build-budget curves; broad C1 comparators remain separately labelled",
        )]),
    }


def run_report(root, trait, layer, output, *, source=None, bootstrap=1000,
               count_bootstrap=0, reference=True, include_gpu=False):
    if (isinstance(bootstrap, bool) or isinstance(count_bootstrap, bool) or
            not np.isfinite(bootstrap) or not np.isfinite(count_bootstrap) or
            int(bootstrap) != bootstrap or int(count_bootstrap) != count_bootstrap or
            bootstrap < 1 or count_bootstrap < 0):
        raise ValueError("bootstrap must be a positive integer; count_bootstrap must be a nonnegative integer")
    bootstrap, count_bootstrap = int(bootstrap), int(count_bootstrap)
    root, out = Path(root), Path(output)
    scope = root / trait / layer
    twin = Path(source) if source is not None else scope / "c1_abm" / "abm_twin"
    if out.exists():
        raise FileExistsError(f"Choose a new report directory; will not replace: {out}")
    # Keep new reports out of the frozen fit tree and its completion receipts.
    frozen_roots = (twin, scope / "c1_correlate", scope / "c1_abm_GPU")
    if any(out.resolve().is_relative_to(path.resolve()) for path in frozen_roots):
        raise ValueError("Report output must be outside frozen ABM/C1/GPU source directories")
    metadata = dict(version=REPORT_VERSION, trait=trait, layer=layer,
                    bootstrap=bootstrap, count_bootstrap=count_bootstrap,
                    training_performed=False, test_predictions_recomputed=False,
                    exploratory=True, sources={},
                    co_primary_metrics="Harrell_C full follow-up strict-time; AUC_IPCW frozen horizon",
                    secondary_metrics="Brier_IPCW; LogLoss_IPCW; Uno_C_horizon",
                    same_training_information_with_C1=False,
                    censoring="reuse frozen marginal training KM, including its original tie convention")
    missing = [name for name in REQUIRED if not (twin / name).is_file()]
    if missing:
        metadata.update(status="not_available", missing=missing)
        tables = {"availability": pd.DataFrame([dict(
            component="frozen_ABM", status="not_available", missing="; ".join(missing),
            explanation="Public aggregates or a preparation workbook cannot replace private frozen predictions")])}
        return _publish(out, tables, metadata)

    sources = {str((twin / name).resolve()): sha256(twin / name) for name in REQUIRED}
    marker, contract = _read_json(twin / REQUIRED[1]), _read_json(twin / "host_contract.json")
    model = abm.load_model(twin / REQUIRED[0])
    pred = abm.load_model(twin / REQUIRED[2])  # trusted local LE8 artifact only
    truth = abm.read_npz(twin / REQUIRED[3])
    frozen_hash = abm.model_hash(model)
    if marker.get("hash") != frozen_hash or pred.get("frozen_model_hash") != frozen_hash:
        raise ValueError("Frozen model, marker and prediction hashes disagree")
    if model.config.trait != trait or model.config.layer != layer:
        raise ValueError("Requested trait/layer differs from frozen model")
    if not np.array_equal(np.asarray(pred["eid"], str), np.asarray(truth["eid"], str)):
        raise ValueError("Frozen prediction/truth IDs or their order differ")
    if len(set(map(str, pred["eid"]))) != len(pred["eid"]):
        raise ValueError("Duplicate frozen participant IDs")
    grouping = getattr(model.config, "grouping_unit", "individual_unverified")
    if contract.get("grouping_unit", grouping) != grouping:
        raise ValueError("Host and frozen model grouping provenance differ")
    print("[LE8 fair evaluation] Frozen fit/roster verified; evaluating saved predictions", flush=True)
    result = abm.evaluate(model, pred, truth["time"], truth["event"], truth["eid"],
                          bootstrap=bootstrap, full_curves=True, fit_calibration=False)
    tables = {"internal_" + name: pd.DataFrame(result[name]) for name in INTERNAL_AGGREGATES if name in result}
    tables["internal_status"] = pd.DataFrame([result["status"]])
    tables.update(_training_tables(model, contract))

    ref, ref_status = None, "disabled_by_user"
    if reference:
        development = twin / "development_LOCAL_ONLY.npz"
        if development.is_file():
            _record_sources(sources, {str(development): sha256(development)})
        print("[LE8 fair evaluation] Building development-only support reference; no refit", flush=True)
        ref, ref_status = _reference_prediction(model, development)
    tables.update(build_support_audit(model, pred, truth["time"], truth["event"], truth["eid"],
                                     reference_prediction=ref, bootstrap=count_bootstrap))
    availability = [dict(component="frozen_ABM", status="evaluated"),
                    dict(component="support_development_reference", status=ref_status)]
    benchmark_path = scope / "c1_correlate/prediction/_LOCAL_ONLY/benchmark.json.gz"
    if benchmark_path.is_file():
        from c1_report import compare
        _record_sources(sources, {str(benchmark_path): sha256(benchmark_path)})
        with gzip.open(benchmark_path, "rt") as stream:
            benchmark = json.load(stream)
        outer_hash = contract.get("outer_roster_sha256")
        if not outer_hash or benchmark.get("lockpaper_protocol", {}).get("outer_roster_sha256") != outer_hash:
            raise ValueError("C1/ABM outer-roster provenance differs or is unavailable")
        comparison = compare(benchmark, pred, truth, model.km, replace(model.config, bootstrap=bootstrap))
        tables.update({"C1_" + key: value for key, value in comparison.items()})
        availability.append(dict(component="C1_shared_test_comparison", status=str(comparison["status"].iloc[0]["status"])))
    else:
        availability.append(dict(component="C1_shared_test_comparison", status="not_available"))
    if include_gpu:
        from c1_report import load_inputs, paired_cohort, performance, primary_tables
        if source is not None and twin.resolve() != (scope / "c1_abm/abm_twin").resolve():
            raise ValueError("GPU comparison requires its original CPU source; do not mix a new experimental fit")
        if not (scope / "c1_abm_GPU/COMPLETE.json").is_file():
            availability.append(dict(component="GPU_shared_test_comparison", status="not_available"))
        else:
            gpu_root = scope / "c1_abm_GPU"
            transformer = gpu_root / "transformer" if (gpu_root / "transformer").is_dir() else gpu_root
            required_gpu = (twin / "COMPLETE.json", gpu_root / "host_contract.json",
                            transformer / "test_predictions_LOCAL_ONLY.tsv.gz", transformer / "MODEL_FROZEN.json")
            missing_gpu = [str(path.relative_to(scope)) for path in required_gpu if not path.is_file()]
            loaded = None
            if not missing_gpu:
                try:
                    loaded = load_inputs(scope)
                except FileNotFoundError as exc:
                    missing_gpu = [str(exc)]
            if missing_gpu:
                availability.append(dict(component="GPU_shared_test_comparison", status="not_available",
                                         detail="; ".join(missing_gpu)))
            else:
                gm, gp, gt, gpu, gc, _, gf, conventional = loaded
                if abm.model_hash(gm) != frozen_hash:
                    raise ValueError("GPU comparison loaded a different frozen CPU model")
                _record_sources(sources, gf)
                data = paired_cohort(gp, gt, gpu, gm.km, conventional)
                data.update(trait=trait, layer=layer, grouping_unit=grouping,
                            outer_roster_sha256=gc["outer_roster_sha256"])
                performance_rows, differences = performance(data, gm.km, bootstrap, model.config.seed)
                tables.update({"GPU_" + key: value for key, value in primary_tables(performance_rows, differences).items()})
                tables["GPU_cohort"] = pd.DataFrame([data["audit"]])
                availability.append(dict(component="GPU_shared_test_comparison", status="evaluated"))
    if abm.model_hash(model) != frozen_hash:
        raise ValueError("Reporting changed the frozen model")
    tables["availability"] = pd.DataFrame(availability)
    metadata.update(status="evaluated", frozen_model_hash=frozen_hash, grouping_unit=grouping,
                    source_scope=str(twin), sources=sources, support_reference=ref_status,
                    components=availability)
    return _publish(out, tables, metadata)


def fair_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--trait", "--Y", dest="trait", required=True)
    parser.add_argument("--layer", "--biom", dest="layer", choices=("prot", "met"), required=True)
    parser.add_argument("--source", type=Path, help="Alternative trusted abm_twin source directory")
    parser.add_argument("--output", type=Path, required=True, help="Fresh directory outside frozen source")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--count-bootstrap", type=int, default=0,
                        help="Optional count-only coverage inference; 0 leaves descriptive count tables")
    parser.add_argument("--skip-reference", action="store_true",
                        help="Skip frozen development audit projection; non-fixed ordering bins then unavailable")
    parser.add_argument("--include-gpu", action="store_true", help="Also verify/re-evaluate saved GPU predictions")
    args = parser.parse_args()
    out = run_report(args.analysis_root, args.trait, args.layer, args.output, source=args.source,
                     bootstrap=args.bootstrap, count_bootstrap=args.count_bootstrap,
                     reference=not args.skip_reference, include_gpu=args.include_gpu)
    print(f"[LE8 fair evaluation] DONE: {out}", flush=True)


# Development training-budget experiment
import argparse
from dataclasses import asdict
import gc
import hashlib
import json
import math
import os
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import abm
from c1_protocol import eid_hash

VERSION = "2026-10-21.abm-training-budget-v1"
DEVELOPMENT_FIELDS = {"X", "time", "event", "eid", "group", "role", "feature_names"}
FEATURE_FIELDS = {"X", "eid", "group", "feature_names"}
INPUT_CONTRACT_FIELDS = ("trait", "layer", "horizon", "grouping_unit", "always_include_indices")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_fractions(values):
    fractions = [float(v) for v in values.split(",")] if isinstance(values, str) else list(map(float, values))
    if (not fractions or len(set(fractions)) != len(fractions) or
            any(not np.isfinite(v) or not 0 < v <= 1 for v in fractions)):
        raise ValueError("Fractions must be distinct finite numbers in (0, 1]")
    return tuple(sorted(fractions))


def _config_values(path):
    document = json.loads(Path(path).read_text())
    if not isinstance(document, dict) or not isinstance(document.get("config", document), dict):
        raise ValueError("Configuration must be a JSON object")
    return dict(document.get("config", document)), "config" in document


def resolve_config_contract(path, source=None):
    """A prepared source owns input identity; an algorithm template cannot reset it."""
    path = Path(path)
    override, frozen_override = _config_values(path)
    values = dict(override)
    audit = dict(override_config_path=str(path), override_config_sha256=sha256_file(path),
                 override_kind="standalone_resolved_config", inherited_fields=[], template_placeholders_replaced={})
    if source is not None:
        source_path = Path(source) / "fit/MODEL_FROZEN.json"
        original, _ = _config_values(source_path)
        missing = set(INPUT_CONTRACT_FIELDS) - set(original)
        if missing:
            raise ValueError("Source frozen configuration lacks input-contract fields: " + ", ".join(sorted(missing)))
        packaged_templates = {Path(__file__).resolve().parents[1] / "config" / "abm.json"}
        template = (not frozen_override and
                    (path.resolve() in packaged_templates or
                     override.get("reference_experiment") == "balanced_reference_20261021"))
        same_source = path.resolve() == source_path.resolve()
        audit.update(source_config_path=str(source_path), source_config_sha256=sha256_file(source_path),
                     override_kind="source_frozen_configuration" if same_source else
                                   "algorithm_template" if template else "explicit_configuration")
        # These are placeholders in the distributed templates, never identities
        # of a prepared source. Non-placeholder conflicting values still fail.
        placeholders = dict(trait="cvd_cad", layer="prot", horizon=10., grouping_unit="family")
        values = {**original, **override}
        for field in INPUT_CONTRACT_FIELDS:
            inherited = field not in override
            if field == "always_include_indices":
                source_value = tuple(original[field])
                override_value = tuple(override.get(field, ()))
                if override_value and override_value != source_value:
                    raise ValueError("always_include_indices conflicts with the prepared source clinical schema")
                inherited = inherited or not override_value
                if inherited and source_value:
                    audit["clinical_schema_resolution"] = "empty_or_missing_indices_inherited_from_source"
            else:
                source_value = original[field]
                override_value = override.get(field, source_value)
                if override_value != source_value:
                    if not template or override_value != placeholders[field]:
                        raise ValueError(field + " conflicts with the prepared source input contract; --config only overrides algorithm settings")
                    inherited = True
                    audit["template_placeholders_replaced"][field] = dict(template=override_value, source=source_value)
            values[field] = source_value
            if inherited:
                audit["inherited_fields"].append(field)
        audit["algorithm_overrides"] = sorted(key for key, value in override.items()
                                               if key not in INPUT_CONTRACT_FIELDS and value != original.get(key))
    for key in ("coverages", "c_grid", "always_include_indices"):
        if key in values:
            values[key] = tuple(values[key])
    config = abm.Config(**values)
    config.validate()
    audit["resolved"] = {field: list(getattr(config, field)) if field == "always_include_indices" else getattr(config, field)
                         for field in INPUT_CONTRACT_FIELDS}
    return config, audit


def resolved_config(path, source=None):
    return resolve_config_contract(path, source)[0]


def load_development(path, *, include_features=False):
    # NpzFile is lazy: a dry run does not decompress the potentially large X.
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != DEVELOPMENT_FIELDS:
            raise ValueError("Development NPZ must contain exactly " + ", ".join(sorted(DEVELOPMENT_FIELDS)))
        out = {key: data[key] for key in data.files if include_features or key != "X"}
    ids, groups, roles = (np.asarray(out[key], str) for key in ("eid", "group", "role"))
    n = len(ids)
    invalid = {"", "nan", "none", "na", "<na>"}
    if (not n or len(set(ids)) != n or any(len(out[key]) != n for key in ("time", "event", "group", "role")) or
            any(v.strip().lower() in invalid for v in np.r_[ids, groups])):
        raise ValueError("Invalid development identifiers or lengths")
    if (set(roles) != set(abm.ROLES) or not np.isfinite(out["time"]).all() or
            np.any(out["time"] <= 0) or not np.isin(out["event"], [0, 1]).all()):
        raise ValueError("Development must contain valid survival outcomes and exactly the five original roles")
    role_by_group = {}
    for group, role in zip(groups, roles):
        if group in role_by_group and role_by_group[group] != role:
            raise ValueError("Family/group crosses development roles")
        role_by_group[group] = role
    names = np.asarray(out["feature_names"], str)
    if not len(names) or len(set(names)) != len(names):
        raise ValueError("Invalid feature schema")
    if include_features and (out["X"].shape != (n, len(names)) or np.isinf(out["X"]).any()):
        raise ValueError("Invalid development feature matrix")
    out.update(eid=ids, group=groups, role=roles, feature_names=names)
    return out


def budget_plan(development, fractions, seed, horizon):
    """Select outcome-blind hash prefixes of whole build groups, preserving rows."""
    fractions = parse_fractions(fractions)
    roles, groups, ids = (development[key] for key in ("role", "group", "eid"))
    build = roles == "build"
    families, counts = np.unique(groups[build], return_counts=True)
    order = sorted(range(len(families)), key=lambda j: (
        hashlib.sha256(f"{VERSION}:{seed}:{families[j]}".encode()).digest(), families[j]))
    ordered_groups = families[order]
    cumulative = np.cumsum(counts[order])
    role_rows = []
    for role in abm.ROLES:
        mask = roles == role
        role_rows.append(dict(role=role, N=int(mask.sum()), group_N=len(set(groups[mask])), eid_hash=eid_hash(ids[mask])))
    budgets, masks = [], {}
    for i, fraction in enumerate(fractions):
        target = int(math.ceil(fraction * int(build.sum())))
        group_count = int(np.searchsorted(cumulative, target, side="left") + 1)
        chosen = np.isin(groups, ordered_groups[:group_count]) & build
        keep = ~build | chosen
        key = f"b{i + 1:02d}_f{fraction:.12g}".replace(".", "p")
        known = (development["time"] > horizon) | ((development["event"] == 1) & (development["time"] <= horizon))
        cases = (development["event"] == 1) & (development["time"] <= horizon)
        budgets.append(dict(budget=key, target_fraction=fraction, target_build_N=target,
            offered_build_N=int(chosen.sum()), offered_build_fraction=float(chosen.sum() / build.sum()),
            offered_build_group_N=group_count, offered_build_known_N=int(np.sum(chosen & known)),
            offered_build_cases=int(np.sum(chosen & cases)), offered_build_eid_hash=eid_hash(ids[chosen]),
            extra_development_N=int((~build).sum()), total_development_N=int(keep.sum()),
            extra_fitting_role_N=int(np.sum((roles != "build") & (roles != "calibration_audit"))),
            extra_development_eid_hash=eid_hash(ids[~build]),
            **{role + "_N": int(np.sum(roles == role)) for role in abm.ROLES if role != "build"}))
        masks[key] = keep
    return dict(version=VERSION, budget_seed=int(seed), horizon=float(horizon),
        roles=role_rows, budgets=budgets,
        sampling="fixed hash order of whole build groups; shortest prefix reaching ceil(fraction * original build N)",
        group_size_rounding="actual N may exceed target N; groups are never split",
        evaluation_censoring_budget=budgets[-1]["budget"],
        unchanged_roles=list(abm.ROLES[1:]),
        scope="ABM-internal equal offered build budget and input schema; not an external C1 clinical-feature match",
        test_truth_read=False, models_fitted=False), masks


def subset_development(development, keep):
    return {key: value.copy() if key == "feature_names" else value[keep].copy()
            for key, value in development.items()}


def load_features(path, development):
    with np.load(path, allow_pickle=False) as source:
        if set(source.files) != FEATURE_FIELDS:
            raise ValueError("Projection NPZ must contain only X, eid, group, feature_names; no outcomes")
        features = {key: source[key] for key in source.files}
    X, ids, groups = abm.TwinABM.validate_features(features["X"], features["eid"], features["group"])
    if not np.array_equal(features["feature_names"].astype(str), development["feature_names"]):
        raise ValueError("Development/test feature schemas differ")
    if X.shape[1] != len(development["feature_names"]):
        raise ValueError("Test feature column count differs")
    if set(ids) & set(development["eid"]) or set(groups) & set(development["group"]):
        raise ValueError("Development/test identifiers or families overlap")
    features.update(X=X, eid=ids, group=groups)
    return features


def read_truth(path, ids):
    with np.load(path, allow_pickle=False) as source:
        if set(source.files) != {"eid", "time", "event"}:
            raise ValueError("Test truth NPZ must contain exactly eid, time, event")
        truth = {key: source[key] for key in source.files}
    if (not np.array_equal(np.asarray(truth["eid"], str), ids) or
            np.asarray(truth["time"]).shape != ids.shape or np.asarray(truth["event"]).shape != ids.shape or
            not np.isfinite(truth["time"]).all() or np.any(truth["time"] <= 0) or
            not np.isin(truth["event"], [0, 1]).all()):
        raise ValueError("Invalid test-truth roster or outcomes")
    return truth


def risk_arms(prediction, config):
    arms = tuple(sorted(name for name in prediction["probability"]
                        if name != "fallback_policy" and not name.startswith("twin_molecular_")))
    required = {"elasticnet_full", "twin_full", "twin_clean"}
    if config.lightgbm:
        required.add("lightgbm_full")
    if not required <= set(arms):
        raise ValueError("A required internal risk arm is missing")
    return arms


def training_arm_counts(model, development, arms):
    """Audit effective fitting counts; do not equate pruned and offered sizes."""
    build = development["role"] == "build"
    y, w = model.km.labels_weights(development["time"][build], development["event"][build])
    base = (w > 0) & model.prep.technical_qc(development["X"][build])
    donor = model.donor_table
    if not np.array_equal(np.asarray(donor["eid"], str), development["eid"][build]):
        raise ValueError("Donor/build roster mismatch")
    eligible = np.asarray(donor["donor_eligible"], bool)
    masks = dict(full=base, clean=base & donor["retained"], soft=base & (donor["soft_weight"] > 0),
                 random=base & donor["random_retained"])
    rows = {}
    for arm in arms:
        if arm.startswith(("elasticnet_", "lightgbm_")):
            mask = masks[arm.rsplit("_", 1)[1]]
            scope = "global model rows with positive training weight"
        elif arm.startswith("twin_"):
            suffix = arm.rsplit("_", 1)[1]
            mask = eligible if suffix == "full" else eligible & donor["retained" if suffix == "clean" else "random_retained"]
            scope = "reference-library donors"
        elif arm == "knn_unweighted_clean":
            mask, scope = eligible & donor["retained"], "reference-library donors"
        else:
            mask, scope = eligible, "eligible build rows for alternative risk route"
        rows[arm] = dict(effective_training_N=int(mask.sum()), effective_training_cases=int(y[mask].sum()),
            effective_training_controls=int(np.sum(mask & (y == 0))), training_count_scope=scope)
    return rows


def evaluate_frozen(entries, truth, features, km, out, plan):
    """All budgets and arms share one outcome-blind QC/finite evaluation mask."""
    ids = features["eid"]
    contract = plan.get("input_contract", {}).get("resolved", {})
    common_qc, common_finite = np.ones(len(ids), bool), np.ones(len(ids), bool)
    predictions = {}
    expected_arms = tuple(entries[0]["arms"])
    for entry in entries:
        path = out / entry["prediction_file"]
        if sha256_file(path) != entry["prediction_sha256"]:
            raise ValueError("Frozen prediction file changed before evaluation")
        pred = abm.load_model(path)
        if (not np.array_equal(np.asarray(pred["eid"], str), ids) or
                not np.array_equal(np.asarray(pred["group"], str), features["group"]) or
                pred["frozen_model_hash"] != entry["model_hash"] or tuple(entry["arms"]) != expected_arms):
            raise ValueError("Budget prediction roster, model hash, or risk arms differ")
        qc = np.asarray(pred["diagnostics"]["technical_qc"], bool)
        if qc.shape != ids.shape:
            raise ValueError("Invalid saved QC mask")
        common_qc &= qc
        for arm in expected_arms:
            p = np.asarray(pred["probability"][arm], float)
            if p.shape != ids.shape or np.any((p[np.isfinite(p)] < 0) | (p[np.isfinite(p)] > 1)):
                raise ValueError("Invalid saved probability")
            common_finite &= np.isfinite(p)
        predictions[entry["budget"]] = pred
    common = common_qc & common_finite
    rows = []
    for entry in entries:
        budget = next(row for row in plan["budgets"] if row["budget"] == entry["budget"])
        for arm in expected_arms:
            result = abm.metrics(truth["time"][common], truth["event"][common],
                                 predictions[entry["budget"]]["probability"][arm][common], km)
            rows.append(dict(**budget, **entry["training_counts"][arm], **result, model=arm, subset="all_common",
                trait=contract.get("trait", ""), layer=contract.get("layer", ""),
                grouping_unit=contract.get("grouping_unit", "individual_unverified"),
                requested_test_N=len(ids), evaluated_N=int(common.sum()), excluded_N=int((~common).sum()),
                excluded_QC_N=int((~common_qc).sum()), excluded_nonfinite_after_QC_N=int(np.sum(common_qc & ~common_finite)),
                requested_test_eid_hash=eid_hash(ids), comparison_eid_hash=eid_hash(ids[common]),
                same_people=True, selection_policy_used=False,
                evaluation_scope="previously_inspected_or_unverified_exploratory",
                evaluation_censoring_budget=plan["evaluation_censoring_budget"],
                evaluation_weights="same frozen largest-budget build-derived censoring KM for every budget and arm",
                IPCW_horizon=km.horizon, inference_status="descriptive_no_bootstrap",
                feature_selection="refit within each budget build; input schema and mandatory clinical columns fixed",
                selected_panel_N=entry["selected_panel_N"], selected_panel_sha256=entry["selected_panel_sha256"]))
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "training_budget_performance.tsv", sep="\t", index=False)
    return frame


def run_experiment(paths, output, config, plan, masks):
    output = Path(output)
    if output.exists():
        raise FileExistsError("Use a fresh output directory; budget experiments are never overwritten")
    development = load_development(paths["development"], include_features=True)
    current_plan, current_masks = budget_plan(development, [row["target_fraction"] for row in plan["budgets"]],
                                              plan["budget_seed"], config.horizon)
    if any(plan.get(key) != value for key, value in current_plan.items()) or any(not np.array_equal(masks[key], current_masks[key]) for key in current_masks):
        raise ValueError("Development inputs changed after budget planning")
    features = load_features(paths["test_features"], development)
    # This hashes inputs needed for fitting/projection only, never test truth.
    source_keys = ("development", "test_features", "config") + (("source_config",) if "source_config" in paths else ())
    sources = {key: dict(path=str(paths[key]), sha256=sha256_file(paths[key])) for key in source_keys}
    if "input_contract" in plan:
        audit = plan["input_contract"]
        actual = {field: list(getattr(config, field)) if field == "always_include_indices" else getattr(config, field)
                  for field in INPUT_CONTRACT_FIELDS}
        if actual != audit["resolved"]:
            raise ValueError("Resolved input contract changed after planning")
        if sources["config"]["sha256"] != audit["override_config_sha256"]:
            raise ValueError("Algorithm configuration changed after planning")
        if "source_config" in sources and sources["source_config"]["sha256"] != audit["source_config_sha256"]:
            raise ValueError("Source input contract changed after planning")
    output.mkdir(parents=True, exist_ok=False)
    frozen_plan = dict(plan, config=asdict(config), sources=sources,
        test_feature_eid_hash=eid_hash(features["eid"]),
        mandatory_clinical_features=[str(development["feature_names"][i]) for i in config.always_include_indices],
        environment={key: os.environ.get(key) for key in ("ABM_FINAL_SELECTOR", "ABM_MOE_CLUSTERS")})
    abm.json_write(output / "PLAN_FROZEN.json", frozen_plan)
    pd.DataFrame(plan["budgets"]).to_csv(output / "training_budget_plan.tsv", sep="\t", index=False)
    entries, evaluation_km = [], None
    for budget in plan["budgets"]:
        directory = output / budget["budget"]
        directory.mkdir()
        subset = subset_development(development, masks[budget["budget"]])
        data_path = directory / "development_LOCAL_ONLY.npz"
        np.savez_compressed(data_path, **subset)
        with threadpool_limits(config.cores):
            # A fresh fit sees only this budget's build rows. No full-build
            # normalizer, selected panel, coefficients or fitted model is reused.
            model = abm.fit_command(data_path, directory / "fit", config)
            model_digest = abm.model_hash(model)
            prediction = model.predict(features["X"], features["eid"], features["group"], features["feature_names"])
            prediction = abm.apply_predictability_extension(model, prediction, features["X"], features["eid"], features["group"])
        if abm.model_hash(model) != model_digest:
            raise ValueError("Projection changed the frozen fitted model")
        arms = risk_arms(prediction, config)
        if entries and arms != tuple(entries[0]["arms"]):
            raise ValueError("Risk-arm availability differs across budgets; evaluation aborted")
        prediction_path = directory / "test_predictions_LOCAL_ONLY.joblib"
        joblib.dump(prediction, prediction_path, compress=3)
        panel = list(map(str, model.panel_names))
        entry = dict(budget=budget["budget"], model_hash=model_digest, arms=list(arms),
            prediction_file=str(prediction_path.relative_to(output)), prediction_sha256=sha256_file(prediction_path),
            training_counts=training_arm_counts(model, subset, arms),
            selected_panel_N=len(panel), selected_panel_sha256=hashlib.sha256("\n".join(panel).encode()).hexdigest())
        entries.append(entry)
        abm.json_write(directory / "PREDICTIONS_FROZEN.json", entry)
        if budget["budget"] == plan["evaluation_censoring_budget"]:
            evaluation_km = model.km
        del model, prediction, subset
        gc.collect()
    # This marker precedes the first open of test truth. Failed/incomplete budget
    # runs never evaluate a subset of the requested experiments.
    abm.json_write(output / "ALL_PREDICTIONS_FROZEN.json", dict(version=VERSION, entries=entries, test_truth_read=False))
    truth = read_truth(paths["test_truth"], features["eid"])
    truth_hash = sha256_file(paths["test_truth"])
    frame = evaluate_frozen(entries, truth, features, evaluation_km, output, plan)
    abm.json_write(output / "COMPLETE.json", dict(version=VERSION, status="completed", budget_N=len(entries),
        all_predictions_frozen_before_test_truth=True, test_truth_sha256=truth_hash,
        performance_sha256=sha256_file(output / "training_budget_performance.tsv"),
        scope=plan["scope"], evaluation="all-only common QC/finite roster; fixed shared censoring model",
        test_provenance="previously_inspected_or_unverified_exploratory",
        external_validation=False, clinical_release=False, confidence_intervals="not computed; descriptive learning curve"))
    return frame


def budget_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", type=Path, help="Existing abm_twin directory containing prepared NPZ files")
    p.add_argument("--development", type=Path)
    p.add_argument("--test-features", type=Path)
    p.add_argument("--test-truth", type=Path)
    p.add_argument("--config", type=Path, help="Algorithm template with --source; otherwise an explicit resolved input config")
    p.add_argument("--fractions", default=".5,.75,1")
    p.add_argument("--budget-seed", type=int, help="Outcome-blind family sampling seed; defaults to fitting config seed")
    p.add_argument("--out", type=Path, help="Fresh experiment directory, required only with --run")
    action = p.add_mutually_exclusive_group()
    action.add_argument("--dry-run", action="store_true", help="Default: print plan without fitting or reading test truth")
    action.add_argument("--run", action="store_true", help="Explicitly perform all requested fresh fits and final evaluation")
    return p


def budget_main(argv=None):
    p = budget_parser()
    args = p.parse_args(argv)
    names = dict(development="development_LOCAL_ONLY.npz", test_features="test_X_LOCAL_ONLY.npz",
                 test_truth="test_truth_LOCAL_ONLY.npz", config="fit/MODEL_FROZEN.json")
    paths = {key: getattr(args, key) or (args.source / name if args.source else None) for key, name in names.items()}
    if paths["development"] is None or paths["config"] is None:
        p.error("Provide --source or both --development and --config")
    if args.run and (args.out is None or any(path is None for path in paths.values())):
        p.error("--run requires --out, test features and test truth paths")
    config, contract = resolve_config_contract(paths["config"], args.source)
    development = load_development(paths["development"])
    if args.source is not None:
        paths["source_config"] = args.source / "fit/MODEL_FROZEN.json"
        with np.load(args.source / "development_LOCAL_ONLY.npz", allow_pickle=False) as original:
            if not np.array_equal(np.asarray(original["feature_names"], str), development["feature_names"]):
                raise ValueError("Prepared feature schema differs from the source clinical input contract")
    if any(i < 0 or i >= len(development["feature_names"]) for i in config.always_include_indices):
        raise ValueError("Mandatory clinical column indices do not match the prepared schema")
    plan, masks = budget_plan(development, args.fractions, config.seed if args.budget_seed is None else args.budget_seed, config.horizon)
    contract["mandatory_clinical_features"] = [str(development["feature_names"][i]) for i in config.always_include_indices]
    plan["input_contract"] = contract
    if not args.run:
        print(json.dumps(plan, ensure_ascii=False, indent=2, allow_nan=False))
        return plan
    return run_experiment(paths, args.out, config, plan, masks)


def main():
    import sys
    actions = {"fair": fair_main, "budget": budget_main}
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        print("C1 internal diagnostics: c1_audit.py {fair|budget} [options]")
        return 0
    action = sys.argv.pop(1)
    if action not in actions:
        raise SystemExit("Unknown C1 audit action: " + action)
    return actions[action]()


if __name__ == "__main__":
    import sys
    sys.exit(main())
