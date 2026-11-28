#!/usr/bin/env python3
"""Optional held-out mortality prognosis for a frozen AOD landmark export.

Usage:
    python aod_validate.py --input aod.tsv.gz --Y cvd_afib --output-dir .

Requires NumPy and SciPy only when this optional module is run. It never fits
clusters, chooses a split, or includes ADuLT. All models use the same participants.
The export must contain a frozen train/test split, assigned clusters, features
measured no later than the landmark, and subsequent mortality/censoring times.
Landmarks before cohort recruitment must be excluded upstream: delayed entry is
not implemented here. IPCW uses a marginal training censoring distribution.
Only aggregate validation.json and validation.tsv are written.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import date
import gzip
import json
import math
from pathlib import Path
import sys

try:
    import numpy as np
    from scipy.optimize import minimize
except ImportError as exc:
    raise SystemExit("Optional AOD prognosis requires numpy and scipy in YGB_PYTHON.") from exc


VERSION = "1.0.0"
MISSING = {"", ".", "NA", "NaN", "nan", "None", "null", "NULL", "N/A"}
MIN_TRAIN_N = 100
MIN_TEST_N = 50
MIN_TRAIN_DEATHS = 20
MIN_TEST_DEATHS = 10
MIN_HORIZON_CASES = 10
MIN_HORIZON_CONTROLS = 10
RIDGE = 1e-4  # Fixed penalty on the event-averaged partial likelihood.


def number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else float("nan")
    except (ValueError, TypeError):
        return float("nan")


def category(value):
    value = str(value).strip() if value is not None else ""
    if value in MISSING:
        return "(missing)"
    numeric = number(value)
    return format(numeric, ".12g") if math.isfinite(numeric) else value


def read_records(path):
    path = Path(path)
    with (gzip.open if path.suffix == ".gz" else open)(
        path, "rt", encoding="utf-8-sig", newline=""
    ) as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Input must have a unique TSV header.")
        result = []
        for rowno, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed TSV row {rowno}.")
            result.append(row)
        return result, reader.fieldnames


class Preprocessor:
    """Train-only median imputation, scaling, and frozen dummy vocabularies."""

    def __init__(self, numeric=(), categorical=()):
        self.numeric = list(numeric)  # (name, transform) pairs
        self.categorical = list(categorical)
        self.specs = []
        self.names = []
        self.dropped = []

    @staticmethod
    def raw_numeric(records, name, transform):
        values = np.array([number(row.get(name)) for row in records], dtype=float)
        if transform == "log1p":
            values[values < 0] = np.nan
            values = np.log1p(values)
        elif transform != "identity":
            raise ValueError(f"Unknown transform {transform}.")
        return values

    def fit(self, records):
        self.specs, self.names, self.dropped = [], [], []
        for name, transform in self.numeric:
            values = self.raw_numeric(records, name, transform)
            finite = np.isfinite(values)
            if not finite.any():
                self.dropped.append({"feature": name, "reason": "all_missing_in_training"})
                continue
            median = float(np.median(values[finite]))
            filled = np.where(finite, values, median)
            mean, scale = float(filled.mean()), float(filled.std())
            spec = {"kind": "numeric", "name": name, "transform": transform,
                    "median": median, "mean": mean, "scale": scale,
                    "missing_indicator": bool(finite.any() and (~finite).any()),
                    "value_column": bool(scale > 1e-12)}
            self.specs.append(spec)
            if spec["value_column"]:
                self.names.append(name if transform == "identity" else f"log1p({name})")
            else:
                self.dropped.append({"feature": name, "reason": "constant_in_training"})
            if spec["missing_indicator"]:
                self.names.append(name + ":missing")
        for name in self.categorical:
            counts = Counter(category(row.get(name)) for row in records)
            if not counts or set(counts) == {"(missing)"}:
                self.dropped.append({"feature": name, "reason": "all_missing_in_training"})
                continue
            reference = sorted(counts, key=lambda key: (-counts[key], key))[0]
            levels = [reference] + sorted(set(counts) - {reference})
            self.specs.append({"kind": "categorical", "name": name, "levels": levels,
                               "reference": reference})
            self.names.extend(f"{name}={level}" for level in levels[1:])
            if len(levels) == 1:
                self.dropped.append({"feature": name, "reason": "constant_in_training"})
        return self

    def transform(self, records):
        columns = []
        for spec in self.specs:
            if spec["kind"] == "numeric":
                values = self.raw_numeric(records, spec["name"], spec["transform"])
                finite = np.isfinite(values)
                if spec["value_column"]:
                    filled = np.where(finite, values, spec["median"])
                    columns.append((filled - spec["mean"]) / spec["scale"])
                if spec["missing_indicator"]:
                    columns.append((~finite).astype(float))
            else:
                values = np.array([category(row.get(spec["name"])) for row in records])
                columns.extend((values == level).astype(float) for level in spec["levels"][1:])
        return np.column_stack(columns) if columns else np.empty((len(records), 0))

    def unseen(self, records):
        return {
            spec["name"]: sum(category(row.get(spec["name"])) not in spec["levels"] for row in records)
            for spec in self.specs if spec["kind"] == "categorical"
        }

    def summary(self):
        return {"columns": self.names, "dropped": self.dropped,
                "fit_scope": "training participants only",
                "numeric_missing": "training median; missing indicator if present in training",
                "categorical_missing": "training missing category when observed",
                "unseen_category": "all dummy columns zero; explicitly reported",
                "parameters": self.specs}


class CoxModel:
    """Cox partial likelihood with Breslow ties and a small fixed ridge penalty."""

    def __init__(self, ridge=RIDGE):
        self.ridge = ridge

    def prepare(self, x, times, events):
        order = np.argsort(times, kind="stable")
        self.x = np.asarray(x, dtype=float)[order]
        self.times = np.asarray(times, dtype=float)[order]
        self.events = np.asarray(events, dtype=int)[order]
        self.unique_times, self.starts = np.unique(self.times, return_index=True)
        all_deaths = np.add.reduceat(self.events, self.starts)
        keep = all_deaths > 0
        self.event_times = self.unique_times[keep]
        self.event_starts = self.starts[keep]
        self.deaths = all_deaths[keep].astype(float)
        self.n_events = int(self.events.sum())
        self.event_x = self.x[self.events == 1].sum(axis=0)
        return self

    def objective(self, beta):
        eta = self.x @ beta
        shift = float(np.max(eta)) if len(eta) else 0.0
        weights = np.exp(eta - shift)
        s0 = np.cumsum(weights[::-1])[::-1]
        s1 = np.cumsum((weights[:, None] * self.x)[::-1], axis=0)[::-1]
        selected = s0[self.event_starts]
        means = np.empty((len(self.event_starts), self.x.shape[1]))
        logs = np.empty(len(self.event_starts))
        stable = selected > 1e-280
        means[stable] = s1[self.event_starts[stable]] / selected[stable, None]
        logs[stable] = np.log(selected[stable]) + shift
        # Extremely separated linear predictors can underflow a late risk set.
        # Recenter only the affected risk set rather than clipping its predictors.
        for j in np.flatnonzero(~stable):
            start = self.event_starts[j]
            local_eta = eta[start:]
            local_shift = float(local_eta.max())
            local_w = np.exp(local_eta - local_shift)
            local_sum = float(local_w.sum())
            logs[j] = math.log(local_sum) + local_shift
            means[j] = (local_w @ self.x[start:]) / local_sum
        value = (float(self.deaths @ logs) - float(eta[self.events == 1].sum())) / self.n_events
        gradient = (self.deaths @ means - self.event_x) / self.n_events
        value += 0.5 * self.ridge * float(beta @ beta)
        gradient += self.ridge * beta
        return value, gradient

    def fit(self, x, times, events):
        self.prepare(x, times, events)
        if self.n_events < 1:
            raise ValueError("Cox fit requires a training death.")
        if self.x.shape[1]:
            fit = minimize(self.objective, np.zeros(self.x.shape[1]), method="L-BFGS-B",
                           jac=True, options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-7})
            self.beta = fit.x
            gradient_norm = float(np.max(np.abs(self.objective(self.beta)[1])))
            if not fit.success or not np.isfinite(self.beta).all() or gradient_norm > 1e-4:
                raise ValueError("Cox optimizer did not converge to the specified gradient tolerance.")
            self.optimization = {"converged": True, "iterations": int(fit.nit),
                                 "maximum_absolute_gradient": gradient_norm}
        else:
            self.beta = np.empty(0)
            self.optimization = {"converged": True, "iterations": 0,
                                 "maximum_absolute_gradient": 0.0}
        eta = self.x @ self.beta
        self.center = float(np.max(eta))
        weights = np.exp(eta - self.center)
        s0 = np.cumsum(weights[::-1])[::-1]
        selected = s0[self.event_starts]
        if np.any(selected <= 0):
            raise ValueError("Unstable baseline hazard; model not estimable.")
        increments = self.deaths / selected
        self.baseline_cumulative = np.cumsum(increments)
        if not np.isfinite(self.baseline_cumulative).all():
            raise ValueError("Nonfinite fitted baseline hazard; model not estimable.")
        return self

    def predict_risk(self, x, horizon):
        index = int(np.searchsorted(self.event_times, horizon, side="right") - 1)
        if index < 0:
            return np.zeros(len(x))
        hazard = float(self.baseline_cumulative[index])
        log_cumulative = math.log(hazard) + np.asarray(x) @ self.beta - self.center
        return -np.expm1(-np.exp(np.clip(log_cumulative, -745, 40)))


class KaplanMeier:
    def __init__(self, times, events):
        order = np.argsort(times, kind="stable")
        times, events = np.asarray(times)[order], np.asarray(events, dtype=int)[order]
        self.times, starts, counts = np.unique(times, return_index=True, return_counts=True)
        self.at_risk = len(times) - starts
        self.deaths = np.add.reduceat(events, starts)
        fraction = 1.0 - self.deaths / self.at_risk
        self.survival = np.cumprod(fraction)
        denominator = self.at_risk * (self.at_risk - self.deaths)
        greenwood = np.divide(self.deaths, denominator, out=np.zeros_like(fraction),
                              where=denominator > 0)
        self.greenwood = np.cumsum(greenwood)
        self.maximum_time = float(times.max())

    def value(self, query, left=False):
        query = np.asarray(query)
        indices = np.searchsorted(self.times, query, side="left" if left else "right") - 1
        return np.where(indices < 0, 1.0, self.survival[np.maximum(indices, 0)])


class Fenwick:
    def __init__(self, n):
        self.tree = np.zeros(n + 1, dtype=np.int64)

    def add(self, index):
        index += 1
        while index < len(self.tree):
            self.tree[index] += 1
            index += index & -index

    def prefix(self, count):
        value = 0
        while count:
            value += int(self.tree[count])
            count -= count & -count
        return value


def harrell_c(times, events, score, tied_tolerance=1e-8):
    """O(n log n); death/censor time ties comparable, death/death ties excluded."""
    times, events, score = np.asarray(times), np.asarray(events), np.asarray(score)
    levels = np.unique(score)
    tree = Fenwick(len(levels))
    rank = np.searchsorted(levels, score)
    order = np.argsort(-times, kind="stable")
    _, starts, sizes = np.unique(-times[order], return_index=True, return_counts=True)
    comparable = concordant = tied = added = 0
    for start, size in zip(starts, sizes):
        group = order[start:start + size]
        censored = group[events[group] == 0]
        deaths = group[events[group] == 1]
        for i in censored:
            tree.add(int(rank[i]))
            added += 1
        for i in deaths:
            lower = int(np.searchsorted(levels, score[i] - tied_tolerance, side="left"))
            upper = int(np.searchsorted(levels, score[i] + tied_tolerance, side="right"))
            less = tree.prefix(lower)
            equal = tree.prefix(upper) - less
            concordant += less
            tied += equal
            comparable += added
        for i in deaths:
            tree.add(int(rank[i]))
            added += 1
    return {"estimate": (concordant + 0.5 * tied) / comparable if comparable else None,
            "comparable_pairs": int(comparable), "concordant_pairs": int(concordant),
            "risk_tied_pairs": int(tied),
            "time_ties": "death/censor comparable; death/death excluded",
            "status": "ok" if comparable else "not_estimable"}


def horizon_metrics(train_times, train_events, test_times, test_events, risks, horizon):
    """Cumulative/dynamic IPCW AUC and unnormalized IPCW Brier mean."""
    t, e, p = np.asarray(test_times), np.asarray(test_events), np.asarray(risks)
    cases = (e == 1) & (t <= horizon)
    controls = t > horizon
    counts = {"cases_by_horizon": int(cases.sum()), "observed_beyond_horizon": int(controls.sum()),
              "censored_by_horizon": int(((e == 0) & (t <= horizon)).sum())}
    unavailable = lambda reason: {"status": "not_estimable", "reason": reason, **counts,
                                 "auc": None, "brier": None}
    censor_km = KaplanMeier(train_times, 1 - np.asarray(train_events))
    g_h = float(censor_km.value(horizon))
    if float(np.max(train_times)) < horizon:
        return unavailable("training follow-up does not support the horizon")
    if g_h <= 0.05:
        return unavailable("training censoring survival at horizon is <=0.05; unstable IPCW")
    if cases.sum() < MIN_HORIZON_CASES or controls.sum() < MIN_HORIZON_CONTROLS:
        return unavailable("fewer than 10 horizon deaths or 10 participants observed beyond horizon")
    w_case = 1.0 / censor_km.value(t[cases], left=True)
    w_control = np.full(int(controls.sum()), 1.0 / g_h)
    control_order = np.argsort(p[controls], kind="stable")
    control_scores = p[controls][control_order]
    control_weight = w_control[control_order]
    cumulative = np.concatenate(([0.0], np.cumsum(control_weight)))
    case_scores = p[cases]
    left = np.searchsorted(control_scores, case_scores - 1e-12, side="left")
    right = np.searchsorted(control_scores, case_scores + 1e-12, side="right")
    favorable = cumulative[left] + 0.5 * (cumulative[right] - cumulative[left])
    auc = float(w_case @ favorable / (w_case.sum() * w_control.sum()))
    brier = float((w_case @ (1.0 - p[cases]) ** 2 +
                   w_control @ p[controls] ** 2) / len(t))
    return {"status": "ok", **counts, "auc": auc, "brier": brier,
            "censoring_survival_at_horizon": g_h, "censoring_fit": "training sample",
            "case_weight": "1 / G(T-)", "control_weight": "1 / G(horizon)",
            "brier_denominator": int(len(t)), "weights_truncated": False}


def km_observed_risk(times, events, horizon):
    if not len(times):
        return {"status": "not_estimable", "estimate": None, "reason": "empty group"}
    km = KaplanMeier(times, events)
    survival = float(km.value(horizon))
    if km.maximum_time < horizon and survival > 0:
        return {"status": "not_estimable", "estimate": None,
                "reason": "group follow-up does not support horizon"}
    result = {"status": "ok", "estimate": 1.0 - survival, "lower95": None, "upper95": None}
    index = int(np.searchsorted(km.times, horizon, side="right") - 1)
    if 0 < survival < 1 and index >= 0:
        se = math.sqrt(float(km.greenwood[index])) / abs(math.log(survival))
        transformed = math.log(-math.log(survival))
        low_s = math.exp(-math.exp(min(700, transformed + 1.96 * se)))
        high_s = math.exp(-math.exp(max(-700, transformed - 1.96 * se)))
        result.update(lower95=1.0 - high_s, upper95=1.0 - low_s)
    elif survival == 0:
        result["interval_note"] = "survival reached zero; Greenwood interval unavailable"
    else:
        result["interval_note"] = "no horizon deaths; Greenwood interval unavailable"
    return result


def calibration(train_risk, test_risk, test_times, test_events, horizon):
    # Cut points come only from training predictions. Holdout outcomes are used
    # solely for evaluation, never for fitting or regrouping the predictions.
    cuts = np.unique(np.quantile(train_risk, [0.2, 0.4, 0.6, 0.8]))
    cuts = cuts[(cuts > np.min(train_risk)) & (cuts < np.max(train_risk))]
    labels = np.searchsorted(cuts, test_risk, side="right")
    groups = []
    for group in range(len(cuts) + 1):
        keep = labels == group
        n = int(keep.sum())
        if not n:
            continue
        groups.append({"group": group + 1, "n": n,
                       "deaths_by_horizon": int(((test_times[keep] <= horizon) &
                                                (test_events[keep] == 1)).sum()),
                       "observed_beyond_horizon": int((test_times[keep] > horizon).sum()),
                       "mean_predicted_risk": float(test_risk[keep].mean()),
                       "observed_risk": km_observed_risk(test_times[keep], test_events[keep], horizon),
                       "small_group": n < 20})
    return {"group_cutpoints": cuts.tolist(), "cutpoints_fit": "training predictions",
            "groups": groups,
            "overall": {"mean_predicted_risk": float(np.mean(test_risk)),
                        "observed_risk": km_observed_risk(test_times, test_events, horizon)},
            "method": "1 - Kaplan-Meier for all-cause mortality; log-log Greenwood intervals",
            "scope": "mortality calibration; no recurrence cumulative-incidence claims"}


def validate_records(records, trait, horizon=3.0, synthetic=False, fieldnames=None):
    if not math.isfinite(horizon) or horizon <= 0:
        raise ValueError("Horizon must be positive.")
    fields = set(fieldnames or (records[0].keys() if records else []))
    cluster_col = f"{trait}.cluster"
    required = {"split", cluster_col, "age_at_landmark", "sex", "landmark_date", "prognosis_reason",
                "death_followup_years", "death_event"}
    out = {
        "schema_version": VERSION, "trait": trait, "endpoint": "all_cause_death",
        "horizon_years": horizon, "status": "not_estimable",
        "data_kind": "synthetic fixture; not a UK Biobank result" if synthetic else
                     "participant export supplied to this run",
        "design": {
            "time_origin": "end of the feature window (landmark)",
            "population": "assigned-cluster participants alive and observed after the landmark",
            "split": "pre-existing train/test split; never created or selected here",
            "same_evaluation_sample": True,
            "cluster_fit": "must have been fitted on training feature histories without future outcomes",
            "outcome_used_for_clustering": False,
            "ADuLT_in_predictors": False,
            "recruitment": "landmark must be on/after recruitment; no delayed-entry fit",
            "availability": "features must have been ascertainable by the landmark",
            "ipcw_assumption": "independent censoring with a common marginal censoring distribution in train/test",
            "cox_ties": "Breslow",
            "ridge": RIDGE,
            "penalty_scale": "event-averaged negative partial likelihood",
            "uncertainty": "point estimates and descriptive paired differences; no delta confidence intervals",
            "eligibility_limits": {"min_train_n": MIN_TRAIN_N, "min_test_n": MIN_TEST_N,
                                   "min_train_deaths": MIN_TRAIN_DEATHS,
                                   "min_test_deaths": MIN_TEST_DEATHS,
                                   "min_training_events_per_parameter": 5},
        },
        "cohort": {"input_rows": len(records)}, "models": {}, "notes": [],
    }
    absent = sorted(required - fields)
    if absent:
        out["reason"] = "Missing required export columns: " + ", ".join(absent)
        return out
    exclusions, selected = Counter(), []
    seen_ids = set()
    for row in records:
        split = str(row.get("split", "")).strip().lower()
        cluster = category(row.get(cluster_col))
        t, event = number(row.get("death_followup_years")), number(row.get("death_event"))
        if split not in {"train", "test"}:
            exclusions["not_train_or_test"] += 1
            continue
        if cluster in {"(missing)", "0", "unassigned", "unknown"}:
            exclusions["no_assigned_cluster"] += 1
            continue
        if str(row.get("prognosis_reason")) != "eligible":
            exclusions["not_eligible_for_prognosis"] += 1
            continue
        try:
            date.fromisoformat(str(row.get("landmark_date", "")))
        except ValueError:
            exclusions["missing_or_invalid_landmark"] += 1
            continue
        if not math.isfinite(t) or t <= 0 or event not in {0.0, 1.0}:
            exclusions["no_valid_positive_mortality_followup"] += 1
            continue
        if "eid" in fields:
            identifier = category(row.get("eid"))
            if identifier != "(missing)" and identifier in seen_ids:
                out["reason"] = "Duplicate participant rows: this validator requires one landmark row per person."
                return out
            seen_ids.add(identifier)
        selected.append(row)
    train = [row for row in selected if row["split"].strip().lower() == "train"]
    test = [row for row in selected if row["split"].strip().lower() == "test"]
    tt = np.array([number(row["death_followup_years"]) for row in train])
    te = np.array([int(number(row["death_event"])) for row in train])
    vt = np.array([number(row["death_followup_years"]) for row in test])
    ve = np.array([int(number(row["death_event"])) for row in test])
    out["cohort"].update(eligible_rows=len(selected), exclusions=dict(exclusions),
                         train_n=len(train), test_n=len(test),
                         train_deaths=int(te.sum()), test_deaths=int(ve.sum()))
    if len(train) < MIN_TRAIN_N or len(test) < MIN_TEST_N:
        out["reason"] = "Need at least 100 training and 50 held-out eligible participants."
        return out
    if te.sum() < MIN_TRAIN_DEATHS or ve.sum() < MIN_TEST_DEATHS:
        out["reason"] = "Need at least 20 training and 10 held-out deaths."
        return out
    if not any(math.isfinite(number(row.get("age_at_landmark"))) for row in train):
        out["reason"] = "Age at landmark is entirely missing in training."
        return out
    if all(category(row.get("sex")) == "(missing)" for row in train):
        out["notes"].append("Sex is unavailable: these models are not adjusted for sex.")
    count_col = "n_bouts" if any(math.isfinite(number(row.get("n_bouts"))) for row in train) else "n_days"
    base_numeric = [("age_at_landmark", "identity")]
    if "disease_duration_years" in fields:
        base_numeric.append(("disease_duration_years", "identity"))
    else:
        out["notes"].append("Disease duration was not exported; time since first occurrence is not adjusted.")
    simple = [(count_col, "log1p"), ("time_since_last_days", "identity")]
    trajectory = [(name, "identity") for name in ("span_years", "active_fraction", "late_minus_early")]
    designs = {
        "M0": ("Age, sex and disease duration", base_numeric, ["sex"], []),
        "M1": ("M0 plus log count and recency", base_numeric + simple, ["sex"],
               [count_col, "time_since_last_days"]),
        "M2": ("M1 plus continuous trajectory", base_numeric + simple + trajectory, ["sex"],
               ["span_years", "active_fraction", "late_minus_early"]),
        "M3": ("M1 plus cluster category", base_numeric + simple, ["sex", cluster_col], [cluster_col]),
    }
    for name, (description, numeric, categorical, incremental) in designs.items():
        summary = {"description": description, "status": "not_estimable", "train_n": len(train),
                   "test_n": len(test), "train_deaths": int(te.sum()), "test_deaths": int(ve.sum())}
        out["models"][name] = summary
        if name in {"M1", "M2", "M3"} and any(
            not any(math.isfinite(number(row.get(column))) for row in train)
            for column in (count_col, "time_since_last_days")
        ):
            summary["reason"] = "Count or recency is entirely missing in training."
            continue
        if name == "M3" and len({category(row.get(cluster_col)) for row in train}) < 2:
            summary["reason"] = "Only one assigned training cluster; no cluster contrast."
            continue
        pre = Preprocessor(numeric, categorical).fit(train)
        x_train, x_test = pre.transform(train), pre.transform(test)
        summary["preprocessing"] = pre.summary()
        summary["unseen_test_categories"] = pre.unseen(test)
        # An unseen cluster cannot be interpreted as the reference cluster.
        if name == "M3" and summary["unseen_test_categories"].get(cluster_col, 0):
            summary["reason"] = "Held-out clusters absent from training; a frozen cluster model is required."
            continue
        if incremental and not any(
            any(column == feature or column.startswith(feature + ":") or
                column.startswith(feature + "=") or column == f"log1p({feature})"
                for column in pre.names)
            for feature in incremental
        ):
            summary["reason"] = "No variable incremental feature is available in training."
            continue
        p = x_train.shape[1]
        summary["parameters"] = p
        if te.sum() < max(MIN_TRAIN_DEATHS, 5 * p):
            summary["reason"] = "Fewer than five training deaths per fitted parameter."
            continue
        try:
            fitted = CoxModel().fit(x_train, tt, te)
        except (ValueError, FloatingPointError, OverflowError) as exc:
            summary["reason"] = str(exc)
            continue
        train_risk = fitted.predict_risk(x_train, horizon)
        test_risk = fitted.predict_risk(x_test, horizon)
        score = x_test @ fitted.beta
        summary.update(
            status="ok", optimization=fitted.optimization,
            coefficients=[{"term": term, "coefficient": float(value)}
                          for term, value in zip(pre.names, fitted.beta)],
            harrell_c=harrell_c(vt, ve, score),
            horizon=horizon_metrics(tt, te, vt, ve, test_risk, horizon),
            calibration=calibration(train_risk, test_risk, vt, ve, horizon),
        )
        if float(tt.max()) < horizon:
            summary["calibration"] = {"status": "not_estimable",
                                      "reason": "training follow-up does not support horizon"}
    reference = out["models"].get("M1", {})
    for name, model in out["models"].items():
        if name == "M1" or model.get("status") != "ok" or reference.get("status") != "ok":
            continue
        delta = {}
        for label, a, b in (
            ("harrell_c", model["harrell_c"].get("estimate"), reference["harrell_c"].get("estimate")),
            ("auc", model["horizon"].get("auc"), reference["horizon"].get("auc")),
            ("brier", model["horizon"].get("brier"), reference["horizon"].get("brier")),
        ):
            delta[label] = a - b if a is not None and b is not None else None
        model["delta_vs_M1"] = delta
    if any(model.get("status") == "ok" for model in out["models"].values()):
        out["status"] = "ok"
    else:
        out["reason"] = "No model met estimation criteria; see per-model reasons."
    return out


def write_outputs(result, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "validation.json").open("w", encoding="utf-8") as out:
        json.dump(result, out, ensure_ascii=False, indent=2, allow_nan=False)
        out.write("\n")
    fields = ["model", "status", "train_n", "test_n", "train_deaths", "test_deaths",
              "harrell_c", "horizon_years", "ipcw_auc", "ipcw_brier", "reason"]
    with (output_dir / "validation.tsv").open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for name, model in result.get("models", {}).items():
            writer.writerow({
                **{key: model.get(key, "") for key in fields if key in model},
                "model": name, "horizon_years": result["horizon_years"],
                "harrell_c": model.get("harrell_c", {}).get("estimate"),
                "ipcw_auc": model.get("horizon", {}).get("auc"),
                "ipcw_brier": model.get("horizon", {}).get("brier"),
            })


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--Y", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--horizon-years", type=float, default=3.0)
    parser.add_argument("--synthetic", action="store_true", help="Label a demonstration as synthetic.")
    args = parser.parse_args(argv)
    try:
        data, fields = read_records(args.input)
        result = validate_records(data, args.Y, args.horizon_years, args.synthetic, fields)
        write_outputs(result, args.output_dir)
    except (ValueError, OSError) as exc:
        print(f"AOD prognosis: {exc}", file=sys.stderr)
        return 2
    print(f"AOD mortality validation: {result['status']}; {args.output_dir}/validation.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
