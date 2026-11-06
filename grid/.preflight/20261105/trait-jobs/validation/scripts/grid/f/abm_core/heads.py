"""Two interpretable objectives learned from a SINGLE twin-count covariate.

Count bins need not be monotone.  The predictability head learns both absolute
base loss and a separately named risk-adjusted loss.  The benefit head learns
ABSOLUTE base-minus-candidate loss; positive means an expected improvement.
No query phenotype is accepted by predict().  Intervals concern calibration-bin
means conditional on the fitted predictors/bank, not the error of one person.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import hashlib

import numpy as np
import pandas as pd

from .losses import row_loss, assessment_loss


def boolean_values(series):
    """CSV-safe booleans: the string 'False' must never become True."""
    mapping = {True: True, False: False, 1: True, 0: False,
               "True": True, "False": False, "true": True, "false": False}
    values = []
    for value in series:
        if pd.isna(value) or value not in mapping:
            raise ValueError("Boolean values must be true/false or 1/0, with no missing values")
        values.append(mapping[value])
    return np.asarray(values, dtype=bool)


def _edges(values, bins):
    values = np.asarray(values, dtype=float)
    if bins < 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Invalid data or number of bins")
    edges = np.unique(np.quantile(values, np.linspace(0, 1, bins + 1)[1:-1]))
    return edges[(edges > values.min()) & (edges < values.max())]


def _identity_checks(frame):
    for col in ("id", "group_id"):
        if col not in frame or frame[col].isna().any():
            raise ValueError(f"Missing {col}; family-unverified IDs must still be explicit")
    if frame["id"].astype(str).duplicated().any():
        raise ValueError("Duplicate individual IDs")


def assert_disjoint(first, second):
    for frame in (first, second):
        _identity_checks(frame)
    for col in ("id", "group_id"):
        if set(first[col].astype(str)) & set(second[col].astype(str)):
            raise ValueError(f"Development roles overlap on {col}")


def _valid_training(frame, task, loss_name):
    _identity_checks(frame)
    required = {"n_twins", "supported", "y", "p_base", "p_candidate"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing training columns: {required - set(frame.columns)}")
    counts = frame["n_twins"].to_numpy(dtype=float)
    support = boolean_values(frame.supported)
    known = np.isfinite(counts)
    if np.isinf(counts).any() or np.any(support & ~known) or (counts[known] < 0).any() or (counts[known] != np.floor(counts[known])).any():
        raise ValueError("n_twins must be an uncapped nonnegative integer count")
    base = row_loss(frame.y, frame.p_base, task, loss_name)
    candidate = row_loss(frame.y, frame.p_candidate, task, loss_name)
    observed_supported = support & np.isfinite(frame.y.to_numpy(float))
    if np.any(observed_supported & (~np.isfinite(base) | ~np.isfinite(candidate))):
        raise ValueError("Both frozen predictors must cover observed supported training rows; "
                         "candidate missingness must not redefine the predictability target")
    use = observed_supported
    if not use.any():
        raise ValueError("No supported, observed, paired development records")
    return frame.loc[use].copy(), base[use], candidate[use]


def _bin_bootstrap(target, groups, prior_mean, prior_strength, n_boot, rng):
    """Resample calibration groups, preserving all observations within a group."""
    unique, inverse = np.unique(np.asarray(groups, dtype=str), return_inverse=True)
    group_n = np.bincount(inverse).astype(float)
    group_sums = np.vstack([np.bincount(inverse, weights=target[:, j])
                            for j in range(target.shape[1])]).T
    estimates = np.empty((n_boot, target.shape[1]))
    for b in range(n_boot):
        multiplicity = rng.multinomial(len(unique), np.full(len(unique), 1 / len(unique)))
        estimates[b] = (multiplicity @ group_sums + prior_strength * prior_mean) / (
            multiplicity @ group_n + prior_strength)
    return np.quantile(estimates, [0.025, 0.975], axis=0)


@dataclass
class CountHeads:
    task: str
    loss_name: str
    count_edges: list
    risk_edges: list
    risk_error_means: list
    bin_estimates: list
    bin_intervals: list
    bin_calibration_n: list
    bin_calibration_groups: list
    count_range: list
    min_calibration_groups: int
    prior_strength: float
    n_boot: int
    seed: int
    bin_calibration_events: list | None = None

    @classmethod
    def unavailable(cls, task, *, loss_name=None, seed=20261020, min_calibration_groups=20, n_boot=1000):
        """An explicit NA head when development cannot support any calibration."""
        return cls(task, assessment_loss(task, loss_name), [], [], [0.],
                   [[np.nan, np.nan, np.nan]], [[[np.nan] * 3, [np.nan] * 3]],
                   [0], [0], [0, 0], min_calibration_groups, 20., n_boot, seed)

    @classmethod
    def fit(cls, gate, calibration, *, task="continuous", loss_name=None,
            n_bins=8, risk_bins=10, prior_strength=20.0,
            min_calibration_groups=20, n_boot=1000, seed=20261020):
        assert_disjoint(gate, calibration)
        if prior_strength < 0 or min_calibration_groups < 2 or n_boot < 20:
            raise ValueError("Invalid prior strength, calibration-group minimum or bootstrap count")
        name = assessment_loss(task, loss_name)
        g, gb, gc = _valid_training(gate, task, name)
        c, cb, cc = _valid_training(calibration, task, name)
        count_edges = _edges(np.log1p(g.n_twins.to_numpy(float)), n_bins)
        risk_edges = _edges(g.p_base.to_numpy(float), risk_bins)
        rg = np.searchsorted(risk_edges, g.p_base.to_numpy(float), side="right")
        rc = np.searchsorted(risk_edges, c.p_base.to_numpy(float), side="right")
        risk_means = np.array([gb[rg == b].mean() if np.any(rg == b) else gb.mean()
                               for b in range(len(risk_edges) + 1)])
        # Column order is contractual: raw base loss, risk-adjusted base loss,
        # raw base-minus-candidate loss.  The latter is NOT risk-residualized.
        tg = np.c_[gb, gb - risk_means[rg], gb - gc]
        tc = np.c_[cb, cb - risk_means[rc], cb - cc]
        bg = np.searchsorted(count_edges, np.log1p(g.n_twins.to_numpy(float)), side="right")
        bc = np.searchsorted(count_edges, np.log1p(c.n_twins.to_numpy(float)), side="right")
        estimates, intervals, ns, ngroups, events = [], [], [], [], []
        rng = np.random.default_rng(seed)
        global_prior = tg.mean(axis=0)
        for b in range(len(count_edges) + 1):
            train = bg == b
            use = bc == b
            prior = tg[train].mean(axis=0) if train.any() else global_prior
            n = int(use.sum())
            groups = c.loc[use, "group_id"].astype(str).to_numpy()
            group_n = len(np.unique(groups))
            estimate = (prior * prior_strength + tc[use].sum(axis=0)) / max(prior_strength + n, 1)
            estimate[0] = max(estimate[0], 0.0)
            if n and group_n >= min_calibration_groups:
                interval = _bin_bootstrap(tc[use], groups, prior, prior_strength, n_boot, rng)
            else:
                interval = np.full((2, 3), np.nan)
            estimates.append(estimate.tolist())
            intervals.append(interval.tolist())
            ns.append(n)
            ngroups.append(group_n)
            events.append(int(c.loc[use, "y"].sum()) if task == "binary" else None)
        return cls(task, name, count_edges.tolist(), risk_edges.tolist(), risk_means.tolist(),
                   estimates, intervals, ns, ngroups,
                   [int(c.n_twins.min()), int(c.n_twins.max())],
                   min_calibration_groups, float(prior_strength), n_boot, seed, events)

    def predict(self, features):
        forbidden = {"y", "outcome", "truth", "loss_base", "loss_candidate", "realized_gain"}
        if forbidden & set(features.columns):
            raise ValueError("Inference features must not include query outcome or realized losses")
        _identity_checks(features)
        counts = features.n_twins.to_numpy(dtype=float)
        support = boolean_values(features.supported)
        known = np.isfinite(counts)
        if np.isinf(counts).any() or np.any(support & ~known) or (counts[known] < 0).any() or (counts[known] != np.floor(counts[known])).any():
            raise ValueError("Invalid twin counts")
        # A genotype-QC failure is unknown support, not zero genetic twins.
        # Preserve its NaN count in the output and keep it in the query roster.
        bins = np.searchsorted(self.count_edges, np.log1p(np.where(known, counts, 0)), side="right")
        estimates = np.asarray(self.bin_estimates, dtype=float)[bins]
        intervals = np.asarray(self.bin_intervals, dtype=float)[bins]
        calibrated = (np.asarray(self.bin_calibration_groups)[bins] >= self.min_calibration_groups)
        in_range = (counts >= self.count_range[0]) & (counts <= self.count_range[1])
        available = support & calibrated & in_range
        reasons = np.full(len(features), "available", dtype=object)
        reasons[~in_range] = "count_outside_calibration_range"
        reasons[~calibrated] = "insufficient_calibration_groups"
        reasons[~support] = "insufficient_genetic_support"
        reasons[~known] = "genetic_count_unavailable"
        out = features.copy().reset_index(drop=True)
        out["count_bin"] = bins
        out["assessment_available"] = available
        out["assessment_reason"] = reasons
        out["head_calibration_n"] = np.asarray(self.bin_calibration_n)[bins]
        out["head_calibration_groups"] = np.asarray(self.bin_calibration_groups)[bins]
        names = ("expected_loss_base", "risk_adjusted_error", "expected_gain")
        for j, name in enumerate(names):
            out[name] = np.where(available, estimates[:, j], np.nan)
            out[name + "_ci_low"] = np.where(available, intervals[:, 0, j], np.nan)
            out[name + "_ci_high"] = np.where(available, intervals[:, 1, j], np.nan)
        if self.task == "continuous":
            out["expected_rmse_base"] = np.sqrt(out.expected_loss_base)
        # Count-only percentile is relative to calibration support, not a
        # guaranteed individual accuracy.  Do not silently invent personal R2.
        return out

    def calibration_table(self):
        rows = []
        for b, estimate in enumerate(self.bin_estimates):
            row = {"count_bin": b, "n": self.bin_calibration_n[b],
                         "n_groups": self.bin_calibration_groups[b],
                         "events": self.bin_calibration_events[b] if self.bin_calibration_events else None,
                         "expected_loss_base": estimate[0],
                         "risk_adjusted_error": estimate[1], "expected_gain": estimate[2],
                         "assessment_available": self.bin_calibration_groups[b] >= self.min_calibration_groups,
                         "interval_scope": "calibration_bin_mean_not_individual_error"}
            for j, name in enumerate(("expected_loss_base", "risk_adjusted_error", "expected_gain")):
                row[name + "_ci_low"] = self.bin_intervals[b][0][j]
                row[name + "_ci_high"] = self.bin_intervals[b][1][j]
            rows.append(row)
        return pd.DataFrame(rows)

    def to_dict(self):
        return asdict(self)


def hash_fraction(ids, seed):
    return np.asarray([int.from_bytes(hashlib.sha256(f"{seed}:{i}".encode()).digest()[:8], "big")
                       / 2**64 for i in np.asarray(ids, dtype=str)])


@dataclass
class FrozenRule:
    """Calibration-frozen threshold plus a label-free deterministic tie rule.

    Coverage is conditional on calibrated, genetically supported eligibility.
    Actual coverage in a new dataset is measured, never forced by reranking it.
    """
    objective: str
    coverage: float
    threshold: float
    tie_cutoff: float
    seed: int
    min_gain: float = 0.0

    @classmethod
    def fit(cls, assessed_calibration, *, objective, coverage=0.5,
            min_gain=0.0, seed=20261020):
        if objective not in {"predictability", "benefit"} or not 0 < coverage <= 1:
            raise ValueError("Invalid selector objective or coverage")
        if not np.isfinite(min_gain) or min_gain < 0:
            raise ValueError("min_gain must be finite and nonnegative")
        use = boolean_values(assessed_calibration.assessment_available)
        scores = (-assessed_calibration.risk_adjusted_error.to_numpy(float) if objective == "predictability"
                  else assessed_calibration.expected_gain.to_numpy(float))
        use &= np.isfinite(scores)
        if not use.any():
            return cls(objective, coverage, np.inf, 0.0, seed, min_gain)
        ids = assessed_calibration.id.astype(str).to_numpy()[use]
        values = scores[use]
        h = hash_fraction(ids, seed)
        order = np.lexsort((h, -values))
        count = min(len(order), max(1, int(np.ceil(coverage * len(order)))))
        cut = order[count - 1]
        if coverage == 1:
            return cls(objective, coverage, -np.inf, 1.0, seed, min_gain)
        return cls(objective, coverage, float(values[cut]), float(h[cut]), seed, min_gain)

    def apply(self, assessed):
        if "y" in assessed:
            raise ValueError("Apply the frozen selector before attaching query outcomes")
        if not np.isfinite(self.min_gain) or self.min_gain < 0:
            raise ValueError("min_gain must be finite and nonnegative")
        available = boolean_values(assessed.assessment_available)
        scores = (-assessed.risk_adjusted_error.to_numpy(float) if self.objective == "predictability"
                  else assessed.expected_gain.to_numpy(float))
        h = hash_fraction(assessed.id.astype(str), self.seed)
        selected = available & np.isfinite(scores) & (
            (scores > self.threshold) | ((scores == self.threshold) & (h <= self.tie_cutoff)))
        if self.objective == "benefit":
            selected &= assessed.expected_gain.to_numpy(float) > self.min_gain
        return selected
