"""Prediction metrics with explicit recorded-incidence and survival estimands.

The recorded task intentionally labels everyone without a recorded event in
the window as a control, matching the HealthFlux primary classification task.
The survival task instead censors at the provided follow-up/competing-death
time. Censoring death estimates net/cause-specific risk, not a competing-risk
cumulative incidence. IPCW assumes independent censoring conditional on the
chosen reference population. Its censoring KM must come from development
data, never be fitted silently on test outcomes.

References: scikit-survival cumulative_dynamic_auc and brier_score API docs.
https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.cumulative_dynamic_auc.html
https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.brier_score.html
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


def _vector(value, dtype=float):
    result = np.asarray(value, dtype=dtype)
    if result.ndim != 1:
        raise ValueError("Expected a one-dimensional array")
    return result


def _same_length(*arrays):
    if len({len(x) for x in arrays}) != 1:
        raise ValueError("Metric arrays have different lengths")


def binary_metrics(labels, scores, sample_weight=None):
    """AUROC and AP; undefined one-class results stay NaN, not zero."""
    labels, scores = _vector(labels, bool), _vector(scores)
    _same_length(labels, scores)
    weights = np.ones(len(labels)) if sample_weight is None else _vector(sample_weight)
    _same_length(labels, weights)
    valid = np.isfinite(scores) & np.isfinite(weights) & (weights > 0)
    y, p, w = labels[valid], scores[valid], weights[valid]
    cases = int(y.sum())
    controls = int((~y).sum())
    result = {"n": len(y), "cases": cases, "controls": controls,
              "auroc": float("nan"), "average_precision": float("nan")}
    if cases and controls:
        result["auroc"] = float(roc_auc_score(y, p, sample_weight=w))
        result["average_precision"] = float(average_precision_score(y, p, sample_weight=w))
    return result


@dataclass(frozen=True)
class CensoringKM:
    """Reverse Kaplan–Meier estimate of the censoring survival function.

    Disease events precede censoring at tied times. Cases use G(t−), controls
    G(horizon). No extrapolation beyond the reference follow-up is allowed.
    """
    times: np.ndarray
    survival: np.ndarray
    max_time: float
    n: int

    @classmethod
    def fit(cls, observed_time, event):
        t, e = _vector(observed_time), _vector(event, bool)
        _same_length(t, e)
        valid = np.isfinite(t) & (t > 0)
        t, e = t[valid], e[valid]
        if not len(t):
            raise ValueError("No positive finite development follow-up for censoring KM")
        order = np.argsort(t, kind="stable")
        t, e = t[order], e[order]
        times, starts, counts = np.unique(t, return_index=True, return_counts=True)
        disease_events = np.add.reduceat(e.astype(np.int64), starts)
        censored = counts - disease_events
        at_risk = len(t) - starts - disease_events
        factors = np.ones(len(times), dtype=float)
        have_censor = censored > 0
        factors[have_censor] = 1.0 - censored[have_censor] / at_risk[have_censor]
        return cls(times, np.cumprod(factors), float(times[-1]), len(t))

    def predict(self, time, *, left_limit=False):
        query = np.asarray(time, dtype=float)
        idx = np.searchsorted(self.times, query, side="left" if left_limit else "right") - 1
        result = np.ones(query.shape, dtype=float)
        nonnegative = idx >= 0
        result[nonnegative] = self.survival[idx[nonnegative]]
        result[query > self.max_time] = np.nan
        return result


def ipcw_weights(observed_time, event, horizon_years, reference, *, weight_cap=20.0):
    """Return cases, controls, weights, status for a cumulative/dynamic task.

    A censored observation at or before the horizon contributes zero weight.
    A nonpositive censoring survival is a support failure, not a weight that
    can be repaired by clipping. Positive extreme weights may be truncated;
    the reported status and counts expose that prespecified truncation.
    """
    t, e = _vector(observed_time), _vector(event, bool)
    _same_length(t, e)
    tau = float(horizon_years)
    if not np.isfinite(tau) or tau <= 0:
        raise ValueError("horizon_years must be positive")
    if weight_cap is not None and (not np.isfinite(weight_cap) or weight_cap < 1):
        raise ValueError("weight_cap must be >= 1 or None")
    cases = e & (t <= tau)
    controls = t > tau
    weights = np.zeros(len(t), dtype=float)
    info = {"status": "ok", "censor_survival_horizon": float("nan"),
            "weights_truncated": 0, "reference_n": 0 if reference is None else reference.n}
    if reference is None:
        info["status"] = "missing_development_censoring_reference"
        return cases, controls, weights, info
    if tau > reference.max_time:
        info["status"] = "horizon_outside_development_followup"
        return cases, controls, weights, info
    g_tau = float(reference.predict(np.array([tau]))[0])
    info["censor_survival_horizon"] = g_tau
    g_cases = reference.predict(t[cases], left_limit=True)
    if not np.isfinite(g_tau) or g_tau <= 0 or np.any(~np.isfinite(g_cases)) or np.any(g_cases <= 0):
        info["status"] = "insufficient_censoring_support"
        return cases, controls, weights, info
    weights[cases] = 1.0 / g_cases
    weights[controls] = 1.0 / g_tau
    if weight_cap is not None:
        info["weights_truncated"] = int(np.sum(weights > weight_cap))
        np.minimum(weights, weight_cap, out=weights)
    if not np.any(cases) or not np.any(controls):
        info["status"] = "no_cases_or_controls_at_horizon"
    return cases, controls, weights, info


def harrell_c_index(observed_time, event, score, *, return_details=False):
    """Harrell concordance in O(n log n), higher score meaning earlier event.

    Tied scores receive half credit. Equal event times are not compared.
    An event and censoring at the same time are comparable, with the event
    ordered first, matching standard survival concordance conventions.
    """
    t, e, s = _vector(observed_time), _vector(event, bool), _vector(score)
    _same_length(t, e, s)
    valid = np.isfinite(t) & (t > 0) & np.isfinite(s)
    t, e, s = t[valid], e[valid], s[valid]
    result = {"c_index": float("nan"), "comparable": 0, "concordant": 0, "tied_scores": 0}
    if not len(t) or not e.any():
        return result if return_details else result["c_index"]
    ranks = np.searchsorted(np.unique(s), s) + 1
    tree = np.zeros(int(ranks.max()) + 1, dtype=np.int64)
    total = 0

    def add(rank):
        nonlocal total
        total += 1
        while rank < len(tree):
            tree[rank] += 1
            rank += rank & -rank

    def prefix(rank):
        ans = 0
        while rank:
            ans += int(tree[rank])
            rank -= rank & -rank
        return ans

    order = np.argsort(-t, kind="stable")
    ordered_time = t[order]
    starts = np.r_[0, np.flatnonzero(ordered_time[1:] != ordered_time[:-1]) + 1]
    ends = np.r_[starts[1:], len(order)]
    for lo, hi in zip(starts, ends):
        group = order[lo:hi]
        censored = group[~e[group]]
        failed = group[e[group]]
        for j in censored:
            add(int(ranks[j]))
        for j in failed:
            rank = int(ranks[j])
            low = prefix(rank - 1)
            tied = prefix(rank) - low
            result["comparable"] += total
            result["concordant"] += low
            result["tied_scores"] += tied
        for j in failed:
            add(int(ranks[j]))
    if result["comparable"]:
        result["c_index"] = (result["concordant"] + 0.5 * result["tied_scores"]) / result["comparable"]
    return result if return_details else float(result["c_index"])


def calibration_bins(probabilities, observed_time, event, *, horizon_years,
                     reference=None, bins=10, weight_cap=20.0, edges=None):
    """Descriptive IPW calibration; returns rows rather than fitting a calibrator.

    All quantiles depend only on predictions. Earlier-censored people remain
    in the reported n but have zero IPCW weight. The weighted incidence is a
    ratio estimator, while the Brier score uses the original n denominator.
    """
    p, t, e = _vector(probabilities), _vector(observed_time), _vector(event, bool)
    _same_length(p, t, e)
    valid = np.isfinite(p) & (p >= 0) & (p <= 1) & np.isfinite(t) & (t > 0)
    p, t, e = p[valid], t[valid], e[valid]
    if not len(p):
        return []
    cases, controls, weights, info = ipcw_weights(t, e, horizon_years, reference, weight_cap=weight_cap)
    if edges is None:
        cuts = np.unique(np.quantile(p, np.linspace(0, 1, max(1, int(bins)) + 1)))
        cuts = np.r_[-np.inf, cuts[1:-1], np.inf]
    else:
        cuts = _vector(edges)
        if len(cuts) < 2 or not np.all(cuts[1:] > cuts[:-1]):
            raise ValueError("Calibration edges must be strictly increasing")
        if cuts[0] > p.min() or cuts[-1] < p.max():
            raise ValueError("Calibration edges do not cover all predictions")
    groups = np.searchsorted(cuts[1:-1], p, side="right")
    rows = []
    for b in range(len(cuts) - 1):
        choose = groups == b
        if not choose.any():
            continue
        wsum = float(weights[choose].sum())
        rows.append({"bin": b + 1, "n": int(choose.sum()), "predicted_mean": float(p[choose].mean()),
                     "predicted_min": float(p[choose].min()), "predicted_max": float(p[choose].max()),
                     "observed_cases": int(cases[choose].sum()), "observed_controls": int(controls[choose].sum()),
                     "ipcw_weight_sum": wsum,
                     "ipcw_observed": float(np.sum(weights[choose] * cases[choose]) / wsum) if wsum else float("nan"),
                     "status": info["status"]})
    return rows


def evaluate_endpoint(event_time, censor_time, scores, *, horizon_years=5.0,
                      ref_time=None, ref_event=None, probabilities=None,
                      probability_horizon_years=None, weight_cap=20.0,
                      risk_mask=None, compute_c_index=True):
    """Evaluate one endpoint using times in years after a fixed landmark.

    event_time is +inf if no disease was recorded. censor_time must be finite
    and is min(administrative end, loss, competing death) for a nonfatal end-
    point. Exclude prevalent disease (event_time <= 0) before all evaluations.
    Ref times/events describe observed survival in the TRAINING risk set.
    """
    event_t, censor_t, score = _vector(event_time), _vector(censor_time), _vector(scores)
    _same_length(event_t, censor_t, score)
    mask = (event_t > 0) & ~np.isnan(event_t) & np.isfinite(censor_t) & (censor_t > 0) & np.isfinite(score)
    if risk_mask is not None:
        risk = _vector(risk_mask, bool)
        _same_length(mask, risk)
        mask &= risk
    probabilities_ = None
    if probabilities is not None:
        probabilities_ = _vector(probabilities)
        _same_length(probabilities_, mask)
        probabilities_ = probabilities_[mask]
    event_t, censor_t, score = event_t[mask], censor_t[mask], score[mask]
    tau = float(horizon_years)
    if not np.isfinite(tau) or tau <= 0:
        raise ValueError("horizon_years must be positive")
    record_labels = event_t <= tau
    record = binary_metrics(record_labels, score)
    observed_time = np.minimum(event_t, censor_t)
    observed_event = np.isfinite(event_t) & (event_t <= censor_t)
    reference = None
    if (ref_time is None) != (ref_event is None):
        raise ValueError("Both ref_time and ref_event are required for IPCW")
    if ref_time is not None:
        try:
            reference = CensoringKM.fit(ref_time, ref_event)
        except ValueError:
            reference = None
    cases, controls, weights, info = ipcw_weights(observed_time, observed_event, tau, reference,
                                                weight_cap=weight_cap)
    result = {"horizon_years": tau, "n": record["n"], "recorded_cases": record["cases"],
              "recorded_controls": record["controls"], "recorded_auroc": record["auroc"],
              "recorded_ap": record["average_precision"], "ipcw_cases": int(cases.sum()),
              "ipcw_controls": int(controls.sum()), "ipcw_auroc": float("nan"),
              "ipcw_brier": float("nan"), "harrell_c_index": float("nan"),
              "ipcw_status": info["status"], "censor_survival_horizon": info["censor_survival_horizon"],
              "weights_truncated": info["weights_truncated"], "censor_reference_n": info["reference_n"],
              "probability_status": "not_supplied", "excluded": int((~mask).sum())}
    if compute_c_index:
        result["harrell_c_index"] = harrell_c_index(observed_time, observed_event, score)
    if info["status"] == "ok":
        result["ipcw_auroc"] = binary_metrics(cases, score, weights)["auroc"]
    if probabilities_ is not None:
        if probability_horizon_years is None or not np.isclose(float(probability_horizon_years), tau):
            result["probability_status"] = "different_or_unknown_probability_horizon"
        elif not np.all(np.isfinite(probabilities_) & (probabilities_ >= 0) & (probabilities_ <= 1)):
            result["probability_status"] = "invalid_probabilities"
        elif info["status"] in ("ok", "no_cases_or_controls_at_horizon") and len(probabilities_):
            result["probability_status"] = "ok"
            result["ipcw_brier"] = float(np.sum(weights * (cases.astype(float) - probabilities_) ** 2) / len(probabilities_))
        else:
            result["probability_status"] = "insufficient_censoring_support"
    return result


def _time_arrays(frame, endpoint=None):
    """Convert the shared prediction-table contract to survival arrays."""
    import pandas as pd
    if {"event_time", "censor_time"}.issubset(frame.columns):
        return frame["event_time"].to_numpy(float), frame["censor_time"].to_numpy(float)
    required = {"landmark_day", "event_day", "exit_day"}
    if not required.issubset(frame.columns):
        raise ValueError("Prediction data need event_time/censor_time or landmark_day/event_day/exit_day")

    def days(column):
        values = frame[column]
        if pd.api.types.is_numeric_dtype(values):
            return pd.to_numeric(values, errors="coerce").to_numpy(float)
        parsed = pd.to_datetime(values, errors="coerce", utc=True)
        return (parsed - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds().to_numpy() / 86400.0

    lm, ev, end = days("landmark_day"), days("event_day"), days("exit_day")
    ev = np.where(np.isfinite(ev), ev, np.inf)
    if "death_day" in frame and str(endpoint).lower() not in {"death", "all_cause_death", "all-cause-death"}:
        death = days("death_day")
        end = np.minimum(end, np.where(np.isfinite(death), death, np.inf))
    return (ev - lm) / 365.25, (end - lm) / 365.25


def evaluate_predictions(predictions, *, reference=None, horizons=(5.0, 10.0),
                         split="test", weight_cap=20.0, compute_c_index=True):
    """Return metric and calibration DataFrames for a long prediction table.

    Columns: eid, split, endpoint, score and times as documented by
    _time_arrays. Optional model, probability, prediction_horizon_years.
    A separate reference may omit scores but must contain training outcome
    times. Without reference, training rows of predictions supply censoring.
    """
    import pandas as pd
    frame = pd.DataFrame(predictions).copy()
    required = {"endpoint", "score", "split"}
    if not required.issubset(frame.columns):
        raise ValueError("Prediction table requires endpoint, score and split columns")
    if "model" not in frame:
        frame["model"] = "flux"
    reference_frame = frame[frame["split"] == "train"] if reference is None else pd.DataFrame(reference)
    if reference is not None and "split" in reference_frame and not np.all(reference_frame["split"].isin(["train", "development"])):
        raise ValueError("Censoring reference must contain development/training rows only")
    rows, calibration = [], []
    for (model, endpoint), group in frame[frame["split"] == split].groupby(["model", "endpoint"], sort=True):
        event_t, censor_t = _time_arrays(group, endpoint)
        ref = reference_frame[reference_frame["endpoint"] == endpoint] if "endpoint" in reference_frame else reference_frame
        if "model" in ref:
            ref = ref[ref["model"] == model]
        ref_event_t, ref_censor_t = _time_arrays(ref, endpoint) if len(ref) else (np.array([]), np.array([]))
        ref_valid = (ref_event_t > 0) & np.isfinite(ref_censor_t) & (ref_censor_t > 0)
        rt = np.minimum(ref_event_t[ref_valid], ref_censor_t[ref_valid])
        re = ref_event_t[ref_valid] <= ref_censor_t[ref_valid]
        probability = group["probability"].to_numpy(float) if "probability" in group else None
        probability_horizons = group["prediction_horizon_years"].dropna().unique() if "prediction_horizon_years" in group else []
        ph = float(probability_horizons[0]) if len(probability_horizons) == 1 else None
        concordance = None
        for tau in horizons:
            result = evaluate_endpoint(event_t, censor_t, group["score"].to_numpy(float), horizon_years=tau,
                                       ref_time=rt, ref_event=re, probabilities=probability,
                                       probability_horizon_years=ph, weight_cap=weight_cap,
                                       compute_c_index=compute_c_index and concordance is None)
            if compute_c_index:
                if concordance is None:
                    concordance = result["harrell_c_index"]
                result["harrell_c_index"] = concordance
            result.update(model=model, endpoint=endpoint, split=split)
            rows.append(result)
            if result["probability_status"] == "ok":
                eligible = (event_t > 0) & np.isfinite(censor_t) & (censor_t > 0) & np.isfinite(group["score"].to_numpy(float))
                km = CensoringKM.fit(rt, re)
                for item in calibration_bins(probability[eligible], np.minimum(event_t, censor_t)[eligible],
                                             (event_t <= censor_t)[eligible], horizon_years=tau,
                                             reference=km, weight_cap=weight_cap):
                    item.update(model=model, endpoint=endpoint, split=split, horizon_years=float(tau))
                    calibration.append(item)
    return pd.DataFrame(rows), pd.DataFrame(calibration)


def macro_average(metrics, *, group_columns=("model", "split", "horizon_years")):
    """Equal-endpoint averages with a separate eligible-endpoint count per metric."""
    import pandas as pd
    frame = pd.DataFrame(metrics)
    if frame.empty:
        return pd.DataFrame()
    numeric = [x for x in ("recorded_auroc", "recorded_ap", "ipcw_auroc", "ipcw_brier", "harrell_c_index") if x in frame]
    result = []
    group_columns = [x for x in group_columns if x in frame]
    groups = frame.groupby(group_columns, dropna=False) if group_columns else [((), frame)]
    for key, group in groups:
        key = key if isinstance(key, tuple) else (key,)
        row = dict(zip(group_columns, key))
        row["endpoints"] = len(group)
        for column in numeric:
            values = pd.to_numeric(group[column], errors="coerce").to_numpy(float)
            values = values[np.isfinite(values)]
            row[column] = float(values.mean()) if len(values) else float("nan")
            row[column + "_endpoints"] = len(values)
        result.append(row)
    return pd.DataFrame(result)
