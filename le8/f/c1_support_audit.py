"""Aggregate-only audit of frozen twin support and prediction error.

``eligible_match_count`` is the number of eligible donors inside the frozen
radius.  The legacy ``matched_count`` is the number *used*, after the K cap.
Neither is disease risk, a probability of correctness, nor an independent
family count.  Fractions normalize reference-bank size, but do not make counts
from different feature geometries or radii directly comparable.

No outcomes enter a gate or a bin boundary here.  The evaluator never fits or
changes a risk model, censoring estimator, selector, or input prediction.  All
returned objects contain aggregates only, with no participant identifiers.
"""
from __future__ import annotations

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
    selector_module = importlib.import_module("c1_selector_inference")
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
    from c1_lock_roster import COVERAGES, PRIMARY_SELECTOR
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
