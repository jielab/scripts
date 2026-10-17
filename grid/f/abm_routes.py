#!/usr/bin/env python3
"""GRID ABM multiroute reference-agent implementation.

This module is deliberately additive to the current GRID twin core.  It keeps the
frozen outer test roster and the five development roles, but fixes two scientific
weaknesses of the 2026-10-14 core:

1. a matched reference bank becomes an actual local prediction route rather than
   only a confidence/support diagnostic; and
2. "predictable" and "benefits from ABM versus fallback PRS" are learned as
   different quantities.  Error/gain targets are residualized within frozen
   fallback-prediction strata so a selector cannot win merely by selecting low-risk
   or low-variance people.

No query/test outcome is accepted by ``apply_abm_routes``.  tune_gate outcomes fit
meta-models; calibration_fit freezes calibrators, percentile maps and coverage
rules.  audit/test outcomes are evaluation-only.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
import pandas as pd
from scipy.special import logit
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.neighbors import KDTree

from twin_core import FrozenRule, Geometry, hash_number, loss

VERSION = "GRID-ABM-FINAL-2026-10-17"


@dataclass(frozen=True)
class ABMRouteConfig:
    enabled: bool = True
    candidate_route: str = "ABM_twin_clean"
    policy_selector: str = "benefit_over_fallback"
    predictability_selector: str = "risk_resid_error"
    primary_coverage: float = 0.50
    coverages: tuple = (0.20, 0.40, 0.50, 0.60, 0.80, 1.00)
    risk_bins: int = 10
    prior_strength: float = 20.0
    gate_iterations: int = 180
    gate_min_leaf: int = 40
    knn_k_grid: tuple = (16, 32, 64, 128)
    moe_clusters: tuple = (4, 8, 12)
    min_cluster_n: int = 100
    min_binary_class_n: int = 10
    bootstrap: int = 1000
    publication_mode: bool = False
    require_verified_groups: bool = True
    require_gwas_nonoverlap: bool = True

    def validate(self):
        routes = {
            "ABM_twin_full", "ABM_twin_clean",
            "ABM_twin_integrated_full", "ABM_twin_integrated_clean",
            "ABM_knn_full", "ABM_knn_clean", "ABM_mixture_experts",
            "ABM_full", "ABM_clean", "ABM_soft", "ABM_random",
        }
        selectors = {
            "risk_resid_error", "logloss_resid_error", "benefit_over_fallback",
            "local_stability", "model_consensus", "hybrid_predictability",
            "hybrid_benefit",
        }
        if self.candidate_route not in routes:
            raise ValueError(f"Unknown ABM candidate_route: {self.candidate_route}")
        if self.policy_selector not in selectors or self.predictability_selector not in selectors:
            raise ValueError("Unknown ABM selector")
        if self.primary_coverage not in self.coverages:
            raise ValueError("ABM primary coverage must be predeclared")
        if any(not 0 < float(q) <= 1 for q in self.coverages) or len(set(self.coverages)) != len(self.coverages):
            raise ValueError("Invalid ABM coverages")
        if self.risk_bins < 3 or self.prior_strength < 0 or self.gate_iterations < 20:
            raise ValueError("Invalid ABM gate settings")
        if self.gate_min_leaf < 5 or self.bootstrap < 0:
            raise ValueError("Invalid ABM support/bootstrap settings")
        if not self.knn_k_grid or min(self.knn_k_grid) < 3:
            raise ValueError("Invalid KNN k grid")
        if not self.moe_clusters or min(self.moe_clusters) < 2:
            raise ValueError("Invalid mixture-of-experts cluster grid")


def load_abm_routes_config(path=None):
    values = {}
    if path and Path(path).is_file():
        values = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("coverages", "knn_k_grid", "moe_clusters"):
        if key in values:
            values[key] = tuple(values[key])
    cfg = ABMRouteConfig(**values)
    cfg.validate()
    return cfg


def _safe_json(value):
    if isinstance(value, dict):
        return {str(k): _safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return _safe_json(value.tolist())
    if isinstance(value, np.generic):
        return _safe_json(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_source_provenance(destination, root, route_config=None):
    """Bind a scientific run to source bytes and, when available, one git commit."""
    root = Path(root).resolve()
    destination = Path(destination)
    paths = [
        root / "grid.sh",
        root / "f/twin.py",
        root / "f/twin_core.py",
        root / "f/abm_routes.py",
        root / "config/twin.json",
        root / "config/abm_routes.json",
    ]
    files = {str(p.relative_to(root)): _sha256(p) for p in paths if p.is_file()}
    commit = None
    try:
        commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True, timeout=5,
        ).stdout.strip()
    except Exception:
        pass
    payload = {
        "format": "GRID-source-provenance-final",
        "abm_version": VERSION,
        "git_commit": commit,
        "files_sha256": files,
        "route_config": asdict(route_config) if route_config else None,
    }
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "SOURCE_FROZEN.json").write_text(
        json.dumps(_safe_json(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return payload


class RouteCalibrator:
    """Independent calibration for ABM secondary routes.

    Binary traits use Platt calibration. Continuous traits use a frozen affine
    map; this is needed for local/KNN/MoE routes whose raw scale need not equal
    the phenotype scale. Calibration labels come only from calibration_fit.
    """
    def fit(self, prediction, y, binary):
        p = np.asarray(prediction, float)
        y = np.asarray(y, float)
        ok = np.isfinite(p) & np.isfinite(y)
        if ok.sum() < 30:
            raise ValueError("Too few calibration rows for ABM route")
        self.binary = bool(binary)
        if self.binary:
            if len(np.unique(y[ok])) < 2:
                raise ValueError("Binary route calibration needs both classes")
            z = logit(np.clip(p[ok], 1e-7, 1 - 1e-7))[:, None]
            self.model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=3000).fit(z, y[ok])
            self.intercept = self.slope = None
        else:
            x = np.c_[np.ones(ok.sum()), p[ok]]
            coef, *_ = np.linalg.lstsq(x, y[ok], rcond=None)
            self.intercept, self.slope = map(float, coef)
            self.model = None
        return self

    def predict(self, prediction):
        p = np.asarray(prediction, float)
        if self.binary:
            valid = np.isfinite(p)
            out = np.full(p.shape, np.nan)
            if valid.any():
                z = logit(np.clip(p[valid], 1e-7, 1 - 1e-7))[:, None]
                out[valid] = self.model.predict_proba(z)[:, 1]
            return out
        return self.intercept + self.slope * p


def _calibrate_radius(reference_z, reference_ids, reference_groups, query_z, query_ids, query_groups, k, quantile):
    """Freeze a kth-neighbour caliper using tune_model only; outcomes are never used."""
    reference_z=np.asarray(reference_z,float); query_z=np.asarray(query_z,float)
    reference_ids=np.asarray(reference_ids,str); reference_groups=np.asarray(reference_groups,str)
    query_ids=np.asarray(query_ids,str); query_groups=np.asarray(query_groups,str)
    if len(reference_z)<k or len(query_z)<10:
        raise ValueError("Too few rows to calibrate ABM reference radius")
    tree=KDTree(reference_z,leaf_size=40); kth=[]; request=min(len(reference_z),int(k)+64)
    for start in range(0,len(query_z),256):
        stop=min(start+256,len(query_z)); dist,ix=tree.query(query_z[start:stop],k=request)
        for j,(ds,nn) in enumerate(zip(dist,ix)):
            good=(reference_ids[nn]!=query_ids[start+j])&(reference_groups[nn]!=query_groups[start+j])
            if good.sum()<k and request<len(reference_z):
                ds,nn=tree.query(query_z[start+j:start+j+1],k=len(reference_z));ds,nn=ds[0],nn[0]
                good=(reference_ids[nn]!=query_ids[start+j])&(reference_groups[nn]!=query_groups[start+j])
            if good.sum()>=k:kth.append(float(ds[good][k-1]))
    if len(kth)<10:raise ValueError("Too few independent tune queries for ABM reference radius")
    radius=float(np.quantile(kth,float(quantile)))
    if not np.isfinite(radius) or radius<=1e-12:raise ValueError("Collapsed ABM reference radius")
    return radius


class ResidualTwinBank:
    """Local OOF-residual borrowing in the same frozen genetic geometry.

    Build donors contribute only their out-of-fold residual. Query outcomes are
    unavailable. Same-family/self donors are excluded.  The correction is shrunk
    toward zero by effective independent-group support.
    """
    def fit(self, z, residual, ids, groups, keep, radius, config, donor_loss=None):
        keep = np.asarray(keep, bool) & np.isfinite(residual)
        if keep.sum() < max(100, config.knn_k_grid[0]):
            raise ValueError("Too few donors for residual twin bank")
        self.z = np.asarray(z, float)[keep]
        self.residual = np.asarray(residual, float)[keep]
        self.ids = np.asarray(ids, str)[keep]
        self.groups = np.asarray(groups, str)[keep]
        self.donor_loss = (
            np.asarray(donor_loss, float)[keep] if donor_loss is not None
            else np.zeros(int(keep.sum()), float)
        )
        self.radius = float(radius)
        self.k = int(max(config.knn_k_grid))
        self.prior_strength = float(config.prior_strength)
        if not np.isfinite(self.radius) or self.radius <= 0:
            raise ValueError("Residual twin bank requires a positive frozen radius")
        return self

    def predict(self, z, ids, groups, fallback, min_twins, min_group_ess, binary):
        q = np.asarray(z, float)
        ids = np.asarray(ids, str)
        groups = np.asarray(groups, str)
        fallback = np.asarray(fallback, float)
        tree = KDTree(self.z, leaf_size=40)
        corrected = fallback.copy()
        n = len(q)
        diag = {
            "eligible_match_count": np.zeros(n),
            "matched_count": np.zeros(n),
            "reference_groups": np.zeros(n),
            "effective_reference_groups": np.zeros(n),
            "donor_ess": np.zeros(n),
            "nearest_radius_ratio": np.full(n, 100.0),
            "radius_ratio": np.full(n, 100.0),
            "local_residual": np.zeros(n),
            "donor_oof_loss_final": np.full(n, 100.0),
            "borrow_fraction": np.zeros(n),
        }
        def neighbors():
            # Bound temporary radius-query memory on full UKB cohorts.
            for start in range(0, n, 256):
                near, distances = tree.query_radius(q[start:start + 256], self.radius, return_distance=True)
                yield from zip(near, distances)
        for i, (idx, ds) in enumerate(neighbors()):
            allowed = (self.ids[idx] != ids[i]) & (self.groups[idx] != groups[i])
            idx, ds = idx[allowed], ds[allowed]
            diag["eligible_match_count"][i] = len(idx)
            if not len(idx):
                continue
            order = np.lexsort((self.ids[idx], ds))[: min(self.k, len(idx))]
            idx, ds = idx[order], ds[order]
            w = np.exp(-0.5 * np.square(ds / self.radius))
            w /= max(w.sum(), 1e-12)
            _, inv = np.unique(self.groups[idx], return_inverse=True)
            group_w = np.bincount(inv, weights=w)
            gess = 1.0 / max(np.sum(group_w ** 2), 1e-12)
            dess = 1.0 / max(np.sum(w ** 2), 1e-12)
            local = float(w @ self.residual[idx])
            borrow = gess / (gess + self.prior_strength) if self.prior_strength > 0 else 1.0
            corrected[i] = fallback[i] + borrow * local
            diag["matched_count"][i] = len(idx)
            diag["reference_groups"][i] = len(group_w)
            diag["effective_reference_groups"][i] = gess
            diag["donor_ess"][i] = dess
            diag["nearest_radius_ratio"][i] = ds[0] / self.radius
            diag["radius_ratio"][i] = ds[-1] / self.radius
            diag["local_residual"][i] = local
            diag["donor_oof_loss_final"][i] = float(w @ self.donor_loss[idx])
            diag["borrow_fraction"][i] = borrow
        if binary:
            corrected = np.clip(corrected, 1e-7, 1 - 1e-7)
        diag["matching_supported_final"] = (
            (diag["matched_count"] >= min_twins)
            & (diag["effective_reference_groups"] >= min_group_ess)
        )
        return corrected, diag


class KNNResidualBank:
    """Independent unweighted-PCA geometry, used as a mechanism control."""
    def fit(self, z, residual, ids, groups, keep, dimensions, seed):
        keep = np.asarray(keep, bool) & np.isfinite(residual)
        if keep.sum() < 100:
            raise ValueError("Too few donors for KNN residual route")
        d = max(1, min(int(dimensions), z.shape[1], int(keep.sum()) - 1))
        self.pca = PCA(n_components=d, svd_solver="randomized", random_state=seed).fit(np.asarray(z)[keep])
        self.z = self.pca.transform(np.asarray(z)[keep]) / math.sqrt(d)
        self.residual = np.asarray(residual, float)[keep]
        self.ids = np.asarray(ids, str)[keep]
        self.groups = np.asarray(groups, str)[keep]
        return self

    def predict(self, z, ids, groups, fallback, k, binary):
        q = self.pca.transform(np.asarray(z)) / math.sqrt(self.pca.n_components_)
        tree = KDTree(self.z, leaf_size=40)
        request = min(len(self.z), max(int(k) + 64, int(k)))
        dist, ix = tree.query(q, k=request)
        out = np.asarray(fallback, float).copy()
        for i, (ds, nn) in enumerate(zip(dist, ix)):
            ok = (self.ids[nn] != str(ids[i])) & (self.groups[nn] != str(groups[i]))
            ds, nn = ds[ok][: int(k)], nn[ok][: int(k)]
            if len(nn) < 3:
                continue
            scale = max(float(ds[-1]), 1e-6)
            w = np.exp(-0.5 * np.square(ds / scale))
            w /= max(w.sum(), 1e-12)
            out[i] += float(w @ self.residual[nn])
        if binary:
            out = np.clip(out, 1e-7, 1 - 1e-7)
        return out


class MixtureExperts:
    """X-only genetic clusters with cluster-specific global outcome models."""
    def _fit_model(self, x, y, binary, seed):
        if binary:
            model = LogisticRegression(C=0.1, solver="lbfgs", max_iter=3000, random_state=seed)
            model.fit(x, y)
            return model
        return Ridge(alpha=100.0).fit(x, y)

    def fit(self, z_cluster, z_model, y, valid, n_clusters, seed, min_n, min_class):
        valid = np.asarray(valid, bool) & np.isfinite(y)
        if valid.sum() < max(min_n * 2, 200):
            raise ValueError("Too few rows for mixture-of-experts")
        self.binary = bool(np.isin(np.unique(np.asarray(y)[valid]), [0, 1]).all() and len(np.unique(np.asarray(y)[valid])) <= 2)
        self.kmeans = KMeans(n_clusters=int(n_clusters), n_init=10, random_state=seed).fit(np.asarray(z_cluster)[valid])
        self.global_model = self._fit_model(np.asarray(z_model)[valid], np.asarray(y)[valid], self.binary, seed)
        labels = self.kmeans.predict(np.asarray(z_cluster))
        self.models = {}
        self.fallback_clusters = []
        for c in range(int(n_clusters)):
            use = valid & (labels == c)
            yy = np.asarray(y)[use]
            usable = use.sum() >= min_n
            if self.binary:
                usable &= len(np.unique(yy)) == 2 and min(np.sum(yy == 0), np.sum(yy == 1)) >= min_class
            if not usable:
                self.models[c] = None
                self.fallback_clusters.append(int(c))
            else:
                self.models[c] = self._fit_model(np.asarray(z_model)[use], yy, self.binary, seed + c + 1)
        return self

    def predict(self, z_cluster, z_model):
        labels = self.kmeans.predict(np.asarray(z_cluster))
        p = np.empty(len(labels), float)
        for c in np.unique(labels):
            use = labels == c
            model = self.models.get(int(c)) or self.global_model
            if self.binary:
                p[use] = model.predict_proba(np.asarray(z_model)[use])[:, 1]
            else:
                p[use] = model.predict(np.asarray(z_model)[use])
        return p, labels


def _primary_row_loss(y, p, binary):
    y = np.asarray(y, float); p = np.asarray(p, float)
    return np.square(y - p)


def _logloss_rows(y, p):
    p = np.clip(np.asarray(p, float), 1e-7, 1 - 1e-7)
    y = np.asarray(y, float)
    return -(y * np.log(p) + (1 - y) * np.log1p(-p))


def _risk_residual(target, fallback, bins):
    """Subtract mean target within prediction strata; no query/test labels here."""
    target = np.asarray(target, float); fallback = np.asarray(fallback, float)
    ok = np.isfinite(target) & np.isfinite(fallback)
    if ok.sum() < 50:
        raise ValueError("Too few tune_gate rows for risk-residualized target")
    q = np.unique(np.quantile(fallback[ok], np.linspace(0, 1, int(bins) + 1)))
    if len(q) < 3:
        edges = np.array([-np.inf, np.inf])
        means = np.array([float(np.mean(target[ok]))])
        labels = np.zeros(len(target), int)
    else:
        q[0], q[-1] = -np.inf, np.inf
        edges = q
        labels = np.searchsorted(edges[1:-1], fallback, side="right")
        global_mean = float(np.mean(target[ok]))
        means = np.array([
            float(np.mean(target[ok & (labels == j)])) if np.any(ok & (labels == j)) else global_mean
            for j in range(len(edges) - 1)
        ])
    return target - means[np.clip(labels, 0, len(means) - 1)], edges, means


def _percentile(values, reference):
    x = np.asarray(values, float)
    ref = np.sort(np.asarray(reference, float)[np.isfinite(reference)])
    out = np.full(len(x), np.nan)
    if len(ref):
        ok = np.isfinite(x)
        out[ok] = np.searchsorted(ref, x[ok], side="right") / len(ref)
    return out


def _local_stability(diag):
    n = len(diag["technical_qc"])
    donor = np.asarray(diag.get("donor_oof_loss_final", diag.get("donor_oof_loss", np.full(n, 100.0))), float)
    radius = np.asarray(diag.get("radius_ratio_final", diag.get("radius_ratio", np.full(n, 100.0))), float)
    nearest = np.asarray(diag.get("nearest_radius_ratio_final", diag.get("nearest_radius_ratio", radius)), float)
    ess = np.asarray(diag.get("effective_reference_groups_final", diag.get("effective_reference_groups", np.zeros(n))), float)
    matched = np.asarray(diag.get("matched_count_final", diag.get("matched_count", np.zeros(n))), float)
    borrow = np.asarray(diag.get("borrow_fraction_final", np.zeros(n)), float)
    return (
        -donor - np.log1p(np.maximum(radius, 0)) - 0.5 * np.log1p(np.maximum(nearest, 0))
        + 0.15 * np.log1p(np.maximum(ess, 0)) + 0.05 * np.log1p(np.maximum(matched, 0))
        + 0.10 * borrow
    )


def _consensus(probabilities):
    names = [n for n in (
        "ABM_full", "ABM_clean", "ABM_soft", "ABM_random",
        "ABM_twin_full", "ABM_twin_clean",
        "ABM_twin_integrated_full", "ABM_twin_integrated_clean",
        "ABM_knn_full", "ABM_knn_clean", "ABM_mixture_experts",
        "CSx_reference", "Disco_reference",
    ) if n in probabilities]
    if len(names) < 2:
        return np.zeros(len(next(iter(probabilities.values())))), names
    mat = np.vstack([np.asarray(probabilities[n], float) for n in names]).T
    return -np.nanstd(mat, axis=1), names


def _gate_features(probabilities, diag, candidate, fallback):
    p1 = np.asarray(probabilities[candidate], float)
    p0 = np.asarray(probabilities[fallback], float)
    consensus, _ = _consensus(probabilities)
    return np.c_[
        p1, p0, p1 - p0,
        np.asarray(diag.get("pre_prediction_confidence", np.zeros(len(p1))), float) / 100.0,
        np.log1p(np.asarray(diag.get("matched_count_final", diag.get("matched_count", np.zeros(len(p1)))), float)),
        np.log1p(np.asarray(diag.get("effective_reference_groups_final", diag.get("effective_reference_groups", np.zeros(len(p1)))), float)),
        np.asarray(diag.get("nearest_radius_ratio_final", diag.get("nearest_radius_ratio", np.full(len(p1), 100.0))), float),
        np.asarray(diag.get("radius_ratio_final", diag.get("radius_ratio", np.full(len(p1), 100.0))), float),
        np.asarray(diag.get("donor_oof_loss_final", diag.get("donor_oof_loss", np.full(len(p1), 100.0))), float),
        np.asarray(diag.get("local_residual_final", np.zeros(len(p1))), float),
        np.asarray(diag.get("borrow_fraction_final", np.zeros(len(p1))), float),
        consensus,
    ]


@dataclass
class ABMExtension:
    config: ABMRouteConfig
    candidate: str
    fallback: str
    twin_full: Any
    twin_clean: Any
    twin_integrated_full: Any
    twin_integrated_clean: Any
    integrated_geometry: Any
    knn_full: Any
    knn_clean: Any
    knn_full_k: int
    knn_clean_k: int
    moe: Any
    moe_clusters: int
    calibrators: dict
    gate_models: dict
    rules: dict
    percentile_reference: dict
    risk_target_edges: dict
    risk_target_means: dict
    consensus_models: tuple
    tuning: list
    audit: dict

    def alternative_probabilities(self, model, X, ids, groups, base_probabilities):
        z = model.prep.transform(X)
        zg = model.geometry.transform(z[:, model.prep.profile_positions])
        p = dict(base_probabilities)
        fallback = np.asarray(p[self.fallback], float)
        binary = model.config.trait_type == "binary"
        raw_full, d_full = self.twin_full.predict(
            zg, ids, groups, fallback, model.config.min_twins, model.config.min_group_ess, binary
        )
        raw_clean, d_clean = self.twin_clean.predict(
            zg, ids, groups, fallback, model.config.min_twins, model.config.min_group_ess, binary
        )
        zg_integrated=self.integrated_geometry.transform(z)
        raw_integrated_full, d_integrated_full=self.twin_integrated_full.predict(
            zg_integrated,ids,groups,fallback,model.config.min_twins,model.config.min_group_ess,binary
        )
        raw_integrated_clean, d_integrated_clean=self.twin_integrated_clean.predict(
            zg_integrated,ids,groups,fallback,model.config.min_twins,model.config.min_group_ess,binary
        )
        graw_full = self.knn_full.predict(z[:, model.prep.profile_positions], ids, groups, fallback, self.knn_full_k, binary)
        graw_clean = self.knn_clean.predict(z[:, model.prep.profile_positions], ids, groups, fallback, self.knn_clean_k, binary)
        moe_raw, cluster = self.moe.predict(z[:, model.prep.profile_positions], z)
        raw = {
            "ABM_twin_full": raw_full,
            "ABM_twin_clean": raw_clean,
            "ABM_twin_integrated_full": raw_integrated_full,
            "ABM_twin_integrated_clean": raw_integrated_clean,
            "ABM_knn_full": graw_full,
            "ABM_knn_clean": graw_clean,
            "ABM_mixture_experts": moe_raw,
        }
        for name, values in raw.items():
            p[name] = self.calibrators[name].predict(values)
            if binary:
                p[name] = np.clip(p[name], 1e-7, 1 - 1e-7)
        diag = {
            "matched_count_final": d_clean["matched_count"],
            "eligible_match_count_final": d_clean["eligible_match_count"],
            "reference_groups_final": d_clean["reference_groups"],
            "effective_reference_groups_final": d_clean["effective_reference_groups"],
            "donor_ess_final": d_clean["donor_ess"],
            "nearest_radius_ratio_final": d_clean["nearest_radius_ratio"],
            "radius_ratio_final": d_clean["radius_ratio"],
            "local_residual_final": d_clean["local_residual"],
            "donor_oof_loss_final": d_clean["donor_oof_loss_final"],
            "borrow_fraction_final": d_clean["borrow_fraction"],
            "matching_supported_final": d_clean["matching_supported_final"],
            "mixture_expert_cluster_final": cluster,
            "integrated_matched_count": d_integrated_clean["matched_count"],
            "integrated_effective_reference_groups": d_integrated_clean["effective_reference_groups"],
            "integrated_radius_ratio": d_integrated_clean["radius_ratio"],
            "integrated_matching_supported": d_integrated_clean["matching_supported_final"],
        }
        return p, diag

    def scores(self, probabilities, diag):
        gx = _gate_features(probabilities, diag, self.candidate, self.fallback)
        raw = {
            "risk_resid_error": -self.gate_models["risk_resid_error"].predict(gx),
            "benefit_over_fallback": self.gate_models["benefit_over_fallback"].predict(gx),
            "local_stability": _local_stability(diag),
        }
        if "logloss_resid_error" in self.gate_models:
            raw["logloss_resid_error"] = -self.gate_models["logloss_resid_error"].predict(gx)
        raw["model_consensus"], _ = _consensus(probabilities)
        pct = {k: _percentile(v, self.percentile_reference[k]) for k, v in raw.items() if k in self.percentile_reference}
        pred_parts = [pct[k] for k in ("risk_resid_error", "local_stability", "model_consensus") if k in pct]
        if "logloss_resid_error" in pct:
            pred_parts.append(pct["logloss_resid_error"])
        raw["hybrid_predictability"] = np.nanmean(np.vstack(pred_parts), axis=0)
        benefit = pct.get("benefit_over_fallback", np.full(len(gx), np.nan))
        raw["hybrid_benefit"] = np.nanmean(
            np.vstack([benefit, raw["hybrid_predictability"], pct.get("local_stability", benefit)]), axis=0
        )
        return raw


def _route_raw_predictions(model, ext_partial, X, ids, groups):
    base, core_diag = model._raw(X, ids, groups)
    p, route_diag = ext_partial.alternative_probabilities(model, X, ids, groups, base)
    diag = dict(core_diag)
    diag.update(route_diag)
    return p, diag


def _route_loss(y, p, binary):
    return _logloss_rows(y, p) if binary else np.square(np.asarray(y) - np.asarray(p))


def _fit_hgb(x, y, known, cfg, seed):
    known = np.asarray(known, bool) & np.isfinite(y) & np.isfinite(x).all(axis=1)
    if known.sum() < max(100, 2 * cfg.gate_min_leaf):
        raise ValueError("Too few supported tune_gate rows for ABM")
    m = HistGradientBoostingRegressor(
        max_iter=cfg.gate_iterations,
        max_leaf_nodes=7,
        learning_rate=0.035,
        min_samples_leaf=max(cfg.gate_min_leaf, min(150, int(known.sum() / 12))),
        l2_regularization=15.0,
        early_stopping=False,
        random_state=seed,
    )
    m.fit(x[known], np.asarray(y)[known])
    return m


def fit_abm_routes(model, development, config=None):
    """Fit ABM using development roles only; never audit/test outcomes."""
    cfg = config or ABMRouteConfig()
    cfg.validate()
    if not cfg.enabled:
        return model
    X = np.asarray(development["X"], float)
    y = np.asarray(development["y"], float)
    ids = np.asarray(development["ids"], str)
    groups = np.asarray(development["groups"], str)
    roles = np.asarray(development["roles"], str)
    if "test" in set(roles):
        raise ValueError("ABM fit accepts development roles only")
    required = {"build", "tune_model", "tune_gate", "calibration_fit", "audit"}
    if set(roles) != required:
        raise ValueError("ABM requires the frozen five-role development design")
    b = roles == "build"; t = roles == "tune_model"; g = roles == "tune_gate"; cal = roles == "calibration_fit"
    z = model.prep.transform(X)
    zg_unweighted = z[:, model.prep.profile_positions]
    zg = model.geometry.transform(zg_unweighted)
    donor = model.donors
    if len(donor["OOF_prediction"]) != b.sum():
        raise ValueError("ABM donor OOF table is not aligned to build role")
    residual = y[b] - np.asarray(donor["OOF_prediction"], float)
    eligible = np.asarray(donor["donor_eligible"], bool)
    clean = np.asarray(donor["retained"], bool) & eligible
    donor_loss = np.asarray(donor["OOF_loss"], float)
    radius = float(model.banks["full"].radius)

    twin_full = ResidualTwinBank().fit(zg[b], residual, ids[b], groups[b], eligible, radius, cfg, donor_loss)
    twin_clean = ResidualTwinBank().fit(zg[b], residual, ids[b], groups[b], clean, radius, cfg, donor_loss)

    # LOCKPAPER analogue: primary twin geometry is GENETIC ONLY.  An integrated
    # genetic+clinical geometry is retained solely as a prespecified sensitivity.
    qc_b=model.prep.qc(X[b]); qc_t=model.prep.qc(X[t])
    integrated_geometry=Geometry().fit(z[b][qc_b],model.predictors["ABM_full"].coefficients,model.config)
    zg_integrated=integrated_geometry.transform(z)
    integrated_radius=_calibrate_radius(
        zg_integrated[b][eligible],ids[b][eligible],groups[b][eligible],
        zg_integrated[t][qc_t],ids[t][qc_t],groups[t][qc_t],model.config.k,model.config.radius_quantile
    )
    twin_integrated_full=ResidualTwinBank().fit(zg_integrated[b],residual,ids[b],groups[b],eligible,integrated_radius,cfg,donor_loss)
    twin_integrated_clean=ResidualTwinBank().fit(zg_integrated[b],residual,ids[b],groups[b],clean,integrated_radius,cfg,donor_loss)
    knn_full = KNNResidualBank().fit(zg_unweighted[b], residual, ids[b], groups[b], eligible, model.config.max_profile_dims, model.config.seed + 2001)
    knn_clean = KNNResidualBank().fit(zg_unweighted[b], residual, ids[b], groups[b], clean, model.config.max_profile_dims, model.config.seed + 2002)

    # A lightweight partial object lets us reuse one route-generation path while tuning.
    class Partial: pass
    partial = Partial()
    partial.config = cfg
    partial.fallback = model.fallback
    partial.twin_full = twin_full; partial.twin_clean = twin_clean
    partial.knn_full = knn_full; partial.knn_clean = knn_clean
    partial.calibrators = {}

    binary = model.config.trait_type == "binary"
    base_t, _ = model._raw(X[t], ids[t], groups[t])
    fallback_t = np.asarray(base_t[model.fallback], float)
    tuning = []
    kgrid = [int(k) for k in cfg.knn_k_grid if int(k) <= max(3, eligible.sum())]
    if not kgrid:
        kgrid = [min(16, int(eligible.sum()))]
    def tune_knn(bank, name):
        best = None
        for k in kgrid:
            pred = bank.predict(zg_unweighted[t], ids[t], groups[t], fallback_t, k, binary)
            score = float(np.mean(_route_loss(y[t], pred, binary)))
            tuning.append({"route": name, "parameter": k, "tune_loss": score})
            if best is None or score < best[0]: best = (score, k)
        return int(best[1])
    k_full = tune_knn(knn_full, "ABM_knn_full")
    k_clean = tune_knn(knn_clean, "ABM_knn_clean")

    moe_best = None
    qc_b = model.prep.qc(X[b])
    for q in cfg.moe_clusters:
        q = int(q)
        try:
            mm = MixtureExperts().fit(
                zg_unweighted[b], z[b], y[b], qc_b, q, model.config.seed + 2100 + q,
                cfg.min_cluster_n, cfg.min_binary_class_n,
            )
            pred, _ = mm.predict(zg_unweighted[t], z[t])
            score = float(np.mean(_route_loss(y[t], pred, binary)))
            tuning.append({"route": "ABM_mixture_experts", "parameter": q, "tune_loss": score})
            if moe_best is None or score < moe_best[0]: moe_best = (score, q, mm)
        except ValueError:
            continue
    if moe_best is None:
        q = int(cfg.moe_clusters[0])
        moe = MixtureExperts().fit(
            zg_unweighted[b], z[b], y[b], qc_b, q, model.config.seed + 2100 + q,
            max(40, cfg.min_cluster_n // 2), max(3, cfg.min_binary_class_n // 2),
        )
        moe_best = (np.nan, q, moe)
    _, moe_q, moe = moe_best

    # Temporary route bundle, first without calibrators/gates/rules.
    partial.knn_full_k = k_full; partial.knn_clean_k = k_clean; partial.moe = moe
    def raw_alt(scope):
        base, core = model._raw(X[scope], ids[scope], groups[scope])
        fallback = np.asarray(base[model.fallback], float)
        tf, df = twin_full.predict(zg[scope], ids[scope], groups[scope], fallback, model.config.min_twins, model.config.min_group_ess, binary)
        tc, dc = twin_clean.predict(zg[scope], ids[scope], groups[scope], fallback, model.config.min_twins, model.config.min_group_ess, binary)
        tif, dif = twin_integrated_full.predict(zg_integrated[scope],ids[scope],groups[scope],fallback,model.config.min_twins,model.config.min_group_ess,binary)
        tic, dic = twin_integrated_clean.predict(zg_integrated[scope],ids[scope],groups[scope],fallback,model.config.min_twins,model.config.min_group_ess,binary)
        kf = knn_full.predict(zg_unweighted[scope], ids[scope], groups[scope], fallback, k_full, binary)
        kc = knn_clean.predict(zg_unweighted[scope], ids[scope], groups[scope], fallback, k_clean, binary)
        me, cluster = moe.predict(zg_unweighted[scope], z[scope])
        raw = {"ABM_twin_full": tf, "ABM_twin_clean": tc,
               "ABM_twin_integrated_full": tif, "ABM_twin_integrated_clean": tic,
               "ABM_knn_full": kf, "ABM_knn_clean": kc, "ABM_mixture_experts": me}
        diag = dict(core)
        diag.update({
            "matched_count_final": dc["matched_count"], "eligible_match_count_final": dc["eligible_match_count"],
            "reference_groups_final": dc["reference_groups"], "effective_reference_groups_final": dc["effective_reference_groups"],
            "donor_ess_final": dc["donor_ess"], "nearest_radius_ratio_final": dc["nearest_radius_ratio"],
            "radius_ratio_final": dc["radius_ratio"], "local_residual_final": dc["local_residual"],
            "donor_oof_loss_final": dc["donor_oof_loss_final"], "borrow_fraction_final": dc["borrow_fraction"],
            "matching_supported_final": dc["matching_supported_final"], "mixture_expert_cluster_final": cluster,
            "integrated_matched_count": dic["matched_count"],
            "integrated_effective_reference_groups": dic["effective_reference_groups"],
            "integrated_radius_ratio": dic["radius_ratio"],
            "integrated_matching_supported": dic["matching_supported_final"],
        })
        return base, raw, diag

    base_cal, raw_cal, diag_cal = raw_alt(cal)
    calibrators = {name: RouteCalibrator().fit(pred, y[cal], binary) for name, pred in raw_cal.items()}
    p_cal = dict(base_cal)
    p_cal.update({name: calibrators[name].predict(pred) for name, pred in raw_cal.items()})
    if binary:
        for name in raw_cal: p_cal[name] = np.clip(p_cal[name], 1e-7, 1 - 1e-7)

    # Candidate is prespecified.  If an optional route is technically unavailable,
    # fail rather than silently pick the test-set winner.
    candidate = cfg.candidate_route
    if candidate not in p_cal:
        if candidate in base_cal:
            pass
        else:
            raise ValueError("Configured ABM candidate route was not fitted: " + candidate)
    fallback = model.fallback

    base_g, raw_g, diag_g = raw_alt(g)
    p_g = dict(base_g); p_g.update({name: calibrators[name].predict(pred) for name, pred in raw_g.items()})
    if binary:
        for name in raw_g: p_g[name] = np.clip(p_g[name], 1e-7, 1 - 1e-7)
    support_g = np.asarray(diag_g["technical_qc"], bool) & np.asarray(diag_g["matching_supported_final"], bool)
    cand_loss = _primary_row_loss(y[g], p_g[candidate], binary)
    base_loss = _primary_row_loss(y[g], p_g[fallback], binary)
    residual_loss, loss_edges, loss_means = _risk_residual(cand_loss, p_g[fallback], cfg.risk_bins)
    gain = base_loss - cand_loss
    residual_gain, gain_edges, gain_means = _risk_residual(gain, p_g[fallback], cfg.risk_bins)
    gx = _gate_features(p_g, diag_g, candidate, fallback)
    gates = {
        "risk_resid_error": _fit_hgb(gx, residual_loss, support_g, cfg, model.config.seed + 3001),
        "benefit_over_fallback": _fit_hgb(gx, residual_gain, support_g, cfg, model.config.seed + 3002),
    }
    risk_edges = {"squared_error": loss_edges, "benefit": gain_edges}
    risk_means = {"squared_error": loss_means, "benefit": gain_means}
    if binary:
        ll = _logloss_rows(y[g], p_g[candidate])
        llr, lle, llm = _risk_residual(ll, p_g[fallback], cfg.risk_bins)
        gates["logloss_resid_error"] = _fit_hgb(gx, llr, support_g, cfg, model.config.seed + 3003)
        risk_edges["logloss"] = lle; risk_means["logloss"] = llm

    # Build raw calibration scores, freeze percentile references, then coverage rules.
    gx_cal = _gate_features(p_cal, diag_cal, candidate, fallback)
    raw_scores = {
        "risk_resid_error": -gates["risk_resid_error"].predict(gx_cal),
        "benefit_over_fallback": gates["benefit_over_fallback"].predict(gx_cal),
        "local_stability": _local_stability(diag_cal),
    }
    if "logloss_resid_error" in gates:
        raw_scores["logloss_resid_error"] = -gates["logloss_resid_error"].predict(gx_cal)
    raw_scores["model_consensus"], consensus_names = _consensus(p_cal)
    eligible_cal = np.asarray(diag_cal["technical_qc"], bool) & np.asarray(diag_cal["matching_supported_final"], bool)
    refs = {name: np.sort(np.asarray(values)[eligible_cal & np.isfinite(values)]) for name, values in raw_scores.items()}
    pct = {k: _percentile(v, refs[k]) for k, v in raw_scores.items()}
    pred_parts = [pct[k] for k in ("risk_resid_error", "local_stability", "model_consensus")]
    if "logloss_resid_error" in pct: pred_parts.append(pct["logloss_resid_error"])
    raw_scores["hybrid_predictability"] = np.nanmean(np.vstack(pred_parts), axis=0)
    raw_scores["hybrid_benefit"] = np.nanmean(
        np.vstack([pct["benefit_over_fallback"], raw_scores["hybrid_predictability"], pct["local_stability"]]), axis=0
    )
    rules = {}
    for name, values in raw_scores.items():
        for q in cfg.coverages:
            rules[(name, float(q))] = FrozenRule().fit(values, eligible_cal, ids[cal], float(q), model.config.seed + 3100)

    audit = {
        "version": VERSION,
        "candidate": candidate,
        "fallback": fallback,
        "policy_selector": cfg.policy_selector,
        "predictability_selector": cfg.predictability_selector,
        "target_semantics": "loss/gain residualized within frozen fallback-prediction strata",
        "risk_only_is_not_predictability": True,
        "tune_gate_N": int(g.sum()),
        "tune_gate_supported_N": int(support_g.sum()),
        "calibration_N": int(cal.sum()),
        "test_outcomes_used": False,
        "audit_outcomes_used": False,
        "local_prediction": "OOF donor residual borrowing with independent-group ESS shrinkage",
        "primary_twin_geometry": "genetic profile only; clinical covariates excluded from distance",
        "integrated_geometry_sensitivity": "integrated genetic plus clinical standardized predictors; sensitivity only",
        "alternative_routes": ["integrated-geometry Twin sensitivity", "independent PCA KNN residual", "X-only KMeans mixture-of-experts"],
        "knn_full_k": k_full, "knn_clean_k": k_clean, "moe_clusters": int(moe_q),
    }
    ext = ABMExtension(
        cfg, candidate, fallback, twin_full, twin_clean, twin_integrated_full, twin_integrated_clean, integrated_geometry,
        knn_full, knn_clean, k_full, k_clean, moe, int(moe_q), calibrators, gates, rules, refs, risk_edges, risk_means,
        tuple(consensus_names), tuning, audit,
    )
    model.grid_abm = ext
    return model


def apply_abm_routes(model, prediction, X, ids, groups):
    """Outcome-blind projection. Keeps the old core selector for backward audit."""
    ext = getattr(model, "grid_abm", None)
    if ext is None or not ext.config.enabled:
        return prediction
    if not np.array_equal(np.asarray(ids, str), np.asarray(prediction["eid"], str)):
        raise ValueError("ABM projection roster differs from core projection")
    if not np.array_equal(np.asarray(groups, str), np.asarray(prediction["group"], str)):
        raise ValueError("ABM projection groups differ from core projection")
    base_predictions=dict(prediction["predictions"])
    if "ABM_policy" in base_predictions:
        base_predictions["ABM_support_policy"]=base_predictions.pop("ABM_policy")
    p, diag_add = ext.alternative_probabilities(model, X, ids, groups, base_predictions)
    diag = dict(prediction["diagnostics"]); diag.update(diag_add)
    scores = ext.scores(p, diag)
    eligible = np.asarray(diag["technical_qc"], bool) & np.asarray(diag["matching_supported_final"], bool)
    masks = dict(prediction.get("masks", {}))
    for name, values in scores.items():
        for q in ext.config.coverages:
            masks[f"{name}:{q:g}"] = ext.rules[(name, float(q))].apply(values, eligible, np.asarray(ids, str))
    policy = masks[f"{ext.config.policy_selector}:{ext.config.primary_coverage:g}"]
    predictable = masks[f"{ext.config.predictability_selector}:{ext.config.primary_coverage:g}"]
    candidate = np.asarray(p[ext.candidate], float)
    fallback = np.asarray(p[ext.fallback], float)
    p["ABM_candidate_route"] = candidate
    p["ABM_policy"] = np.where(np.asarray(diag["technical_qc"], bool), np.where(policy, candidate, fallback), np.nan)
    out = dict(prediction)
    out["selected_support"] = np.asarray(prediction.get("selected"),bool).copy()
    out["selected"] = predictable
    out["reason"] = np.where(
        ~np.asarray(diag["technical_qc"],bool), prediction.get("reason","technical_QC_failed"),
        np.where(~np.asarray(diag["matching_supported_final"],bool),"no_adequate_genetic_twins",
                 np.where(predictable,"predictability_selected","predictability_deferred"))
    )
    out["predictions"] = p
    out["diagnostics"] = diag
    out["masks"] = masks
    out["selector_scores_final"] = scores
    out["selected_policy"] = policy
    out["selected_predictability"] = predictable
    out["candidate_route"] = ext.candidate
    out["fallback_route"] = ext.fallback
    out["policy_selector"] = ext.config.policy_selector
    out["predictability_selector"] = ext.config.predictability_selector
    out["abm_version"] = VERSION
    out["selection_semantics_final"] = (
        "primary selector = risk_resid_error; incremental-benefit selector is secondary and used only for ABM_policy. "
        "Both targets were frozen on development data after residualization within fallback-prediction strata; audit/test outcomes never set thresholds."
    )
    return out


def _metric_row(y, p, binary):
    y = np.asarray(y, float); p = np.asarray(p, float)
    if not len(y):
        return {"N": 0, "primary_loss": np.nan, "MSE": np.nan, "prediction_R2": np.nan,
                "AUC": np.nan, "AUPRC": np.nan, "Brier": np.nan, "LogLoss": np.nan}
    if binary:
        pp = np.clip(p, 1e-7, 1 - 1e-7)
        two = len(np.unique(y)) == 2
        brier = float(np.mean(np.square(y - p)))
        return {
            "N": len(y), "primary_loss": brier,
            "AUC": float(roc_auc_score(y, p)) if two else np.nan,
            "AUPRC": float(average_precision_score(y, p)) if two else np.nan,
            "Brier": brier,
            "LogLoss": float(np.mean(_logloss_rows(y, pp))),
            "MSE": np.nan, "prediction_R2": np.nan,
        }
    mse = float(np.mean(np.square(y - p)))
    var = float(np.var(y))
    return {
        "N": len(y), "primary_loss": mse, "MSE": mse,
        "prediction_R2": 1 - mse / var if var > 0 else np.nan,
        "AUC": np.nan, "AUPRC": np.nan, "Brier": np.nan, "LogLoss": np.nan,
    }


def _composition_matched(selected, eligible, fallback, ancestry, ids, edges, seed):
    """Exact count within ancestry x frozen fallback-prediction bins; outcome-blind."""
    selected = np.asarray(selected, bool); eligible = np.asarray(eligible, bool)
    fallback = np.asarray(fallback, float); ancestry = np.asarray(ancestry, str); ids = np.asarray(ids, str)
    bins = np.searchsorted(np.asarray(edges), fallback, side="right")
    out = np.zeros(len(ids), bool)
    for pop in np.unique(ancestry):
        for b in np.unique(bins[ancestry == pop]):
            scope = (ancestry == pop) & (bins == b)
            pool = np.flatnonzero(eligible & scope)
            count = int(np.sum(selected & scope))
            order = sorted(pool, key=lambda j: hash_number(ids[j], seed))
            out[np.asarray(order[:count], int)] = True
    return out


def _aurc(frame, column):
    d = frame[np.isfinite(frame["actual_coverage"]) & np.isfinite(frame[column])].sort_values("actual_coverage")
    if len(d) < 2 or d.actual_coverage.max() <= d.actual_coverage.min(): return np.nan
    return float(np.trapz(d[column], d.actual_coverage) / (d.actual_coverage.max() - d.actual_coverage.min()))


def _bootstrap_loss_contrast(y, candidate, fallback, selected, control, groups, binary, B, seed):
    if B <= 0:
        return {"candidate_minus_fallback_lower":np.nan,"candidate_minus_fallback_upper":np.nan,"candidate_valid_boot":0,
                "selector_minus_control_lower":np.nan,"selector_minus_control_upper":np.nan,"selector_valid_boot":0}
    y = np.asarray(y, float); candidate = np.asarray(candidate, float); fallback = np.asarray(fallback, float)
    groups = np.asarray(groups, str); selected = np.asarray(selected, bool); control = np.asarray(control, bool)
    ug, inv = np.unique(groups, return_inverse=True)
    rng = np.random.default_rng(seed); a=[]; s=[]
    def weighted_loss(p, mask, gw):
        valid = np.asarray(mask, bool) & np.isfinite(y) & np.isfinite(p)
        w = gw[inv] * valid
        total = w.sum()
        if total <= 0: return np.nan
        err = np.where(valid, np.square(y - p), 0.0)
        return float(np.sum(w * err) / total)
    for _ in range(B):
        gw = np.bincount(rng.integers(0, len(ug), len(ug)), minlength=len(ug))
        # candidate-vs-fallback on the identical selected people
        a.append(weighted_loss(candidate, selected, gw) - weighted_loss(fallback, selected, gw))
        # selector composition effect on a fixed fallback model
        s.append(weighted_loss(fallback, selected, gw) - weighted_loss(fallback, control, gw))
    def summarise(x):
        x=np.asarray(x,float);x=x[np.isfinite(x)]
        if len(x)<max(20,math.ceil(.8*B)): return (np.nan,np.nan,len(x))
        q=np.quantile(x,[.025,.975]);return (float(q[0]),float(q[1]),len(x))
    al,au,an=summarise(a);sl,su,sn=summarise(s)
    return {"candidate_minus_fallback_lower":al,"candidate_minus_fallback_upper":au,"candidate_valid_boot":an,
            "selector_minus_control_lower":sl,"selector_minus_control_upper":su,"selector_valid_boot":sn}


def _bootstrap_same_subset_delta(y, first, second, mask, groups, B, seed):
    """Group-bootstrap loss(first)-loss(second) on identical frozen participants."""
    y=np.asarray(y,float);first=np.asarray(first,float);second=np.asarray(second,float)
    mask=np.asarray(mask,bool)&np.isfinite(y)&np.isfinite(first)&np.isfinite(second)
    groups=np.asarray(groups,str);ug,inv=np.unique(groups,return_inverse=True)
    def point(weights=None):
        if weights is None:
            w=mask.astype(float)
        else:
            w=weights[inv]*mask
        den=w.sum()
        if den<=0:return np.nan
        delta = np.square(y[mask]-first[mask])-np.square(y[mask]-second[mask])
        return float(np.sum(w[mask]*delta)/den)
    delta=point();vals=[];rng=np.random.default_rng(seed)
    for _ in range(int(B)):
        gw=np.bincount(rng.integers(0,len(ug),len(ug)),minlength=len(ug))
        vals.append(point(gw))
    vals=np.asarray(vals,float);vals=vals[np.isfinite(vals)]
    enough=len(vals)>=max(20,math.ceil(.8*max(1,int(B))))
    lo,hi=(np.quantile(vals,[.025,.975]) if enough else (np.nan,np.nan))
    return dict(delta_loss=delta,lower=float(lo),upper=float(hi),valid_boot=int(len(vals)),requested_boot=int(B))


def publication_status(model, config, external_validation=False):
    grouping = str(getattr(model.config, "grouping_status", "individual_unverified"))
    overlap = str(getattr(model.config, "upstream_gwas_overlap", "unknown")).strip().lower()
    group_ok = grouping not in {"individual_unverified", "unknown", "unverified"}
    overlap_ok = overlap in {"verified_nonoverlap", "verified_no_overlap", "no_overlap", "none"}
    blockers=[]
    if config.require_verified_groups and not group_ok: blockers.append("relatedness/grouping_not_verified")
    if config.require_gwas_nonoverlap and not overlap_ok: blockers.append("upstream_GWAS_target_overlap_unknown")
    if not external_validation: blockers.append("no_external_or_new_untouched_validation")
    return {
        "publication_ready_confirmatory": not blockers,
        "grouping_status": grouping,
        "grouping_verified": group_ok,
        "upstream_gwas_overlap": overlap,
        "gwas_nonoverlap_verified": overlap_ok,
        "external_validation": bool(external_validation),
        "blockers": ";".join(blockers),
        "current_claim_scope": "exploratory_internal_validation" if blockers else "confirmatory_candidate",
    }


def evaluate_abm_routes(model, prediction, truth, config=None, bootstrap=None):
    """ABM selective-prediction evaluation on frozen masks.

    Selector quality is measured on a fixed fallback model versus an exact-count
    ancestry x frozen-risk matched outcome-blind control. Candidate gain is measured
    on the same selected people.  These two questions are never conflated.
    """
    ext = getattr(model, "grid_abm", None)
    if ext is None:
        return {}
    cfg = config or ext.config
    truth = pd.DataFrame(truth).copy(); truth["eid"] = truth.eid.astype(str)
    ids = np.asarray(prediction["eid"], str)
    if truth.eid.duplicated().any() or set(truth.eid) != set(ids):
        raise ValueError("ABM truth must cover exactly the frozen query roster")
    truth = truth.set_index("eid").loc[ids]
    y = pd.to_numeric(truth.y, errors="coerce").to_numpy(float)
    known = np.isfinite(y)
    ancestry = truth.ancestry.fillna("UNASSIGNED").astype(str).to_numpy() if "ancestry" in truth else np.full(len(ids), "UNASSIGNED")
    qc = np.asarray(prediction["diagnostics"]["technical_qc"], bool)
    support = np.asarray(prediction["diagnostics"]["matching_supported_final"], bool)
    eligible = known & qc & support
    p = prediction["predictions"]
    candidate = np.asarray(p["ABM_candidate_route"], float)
    fallback = np.asarray(p[ext.fallback], float)
    policy = np.asarray(p["ABM_policy"], float)
    binary = model.config.trait_type == "binary"
    # These bins are frozen during development, including when no test labels exist.
    edges = model.risk_edges
    rows=[]
    route_map = {
        "predictability": ext.config.predictability_selector,
        "benefit": ext.config.policy_selector,
    }
    for purpose, route in route_map.items():
        for q in ext.config.coverages:
            key=f"{route}:{q:g}"
            membership = np.asarray(prediction["masks"][key], bool)
            control_membership = _composition_matched(membership, qc & support, fallback, ancestry, ids, edges, model.config.seed + 4200 + int(q*1000))
            selected = known & membership
            control = known & control_membership
            cm=_metric_row(y[selected],candidate[selected],binary); fm=_metric_row(y[selected],fallback[selected],binary); rm=_metric_row(y[control],fallback[control],binary)
            total_cases=np.sum(known & (y==1)) if binary else 0
            rows.append({
                "purpose":purpose,"selector":route,"coverage_target":q,
                "actual_coverage":float(membership.sum()/len(ids)),"selected_N":int(selected.sum()),"control_N":int(control.sum()),
                "selected_roster_N":int(membership.sum()),"control_roster_N":int(control_membership.sum()),
                "case_coverage":float(np.sum(selected&(y==1))/total_cases) if binary and total_cases else np.nan,
                "candidate_primary_loss":cm["primary_loss"],"fallback_primary_loss":fm["primary_loss"],
                "risk_ancestry_control_primary_loss":rm["primary_loss"],
                "candidate_gain_same_subset":fm["primary_loss"]-cm["primary_loss"],
                "selector_delta_loss_vs_control":fm["primary_loss"]-rm["primary_loss"],
                "candidate_R2":cm["prediction_R2"],"fallback_R2":fm["prediction_R2"],
                "candidate_AUC":cm["AUC"],"fallback_AUC":fm["AUC"],
                "candidate_Brier":cm["Brier"],"fallback_Brier":fm["Brier"],
                "candidate_LogLoss":cm["LogLoss"],"fallback_LogLoss":fm["LogLoss"],
                "control_AUC":rm["AUC"],"control_Brier":rm["Brier"],"control_LogLoss":rm["LogLoss"],
                "control_definition":"exact count within ancestry x frozen fallback-prediction bins; no outcomes used for membership",
            })
    curve=pd.DataFrame(rows)
    summaries=[]
    for purpose,g in curve.groupby("purpose",sort=False):
        summaries.append({
            "purpose":purpose,"selector":g.selector.iloc[0],
            "AURC_candidate_gain":_aurc(g,"candidate_gain_same_subset"),
            "AURC_selector_delta_loss_vs_control":_aurc(g,"selector_delta_loss_vs_control"),
            "mean_case_coverage":float(np.nanmean(g.case_coverage)) if binary else np.nan,
            "interpretation":"gain>0 means candidate beats fallback on identical selected people; selector delta<0 means selected people have lower fixed-fallback loss than a risk/ancestry matched outcome-blind control",
        })
    primary=[]
    B=cfg.bootstrap if bootstrap is None else int(bootstrap)
    for purpose,route in route_map.items():
        q=cfg.primary_coverage; membership=np.asarray(prediction["masks"][f"{route}:{q:g}"],bool)
        control_membership=_composition_matched(membership,qc&support,fallback,ancestry,ids,edges,model.config.seed+4200+int(q*1000))
        selected=known & membership; control=known & control_membership
        boot=_bootstrap_loss_contrast(y,candidate,fallback,selected,control,np.asarray(prediction["group"],str),binary,B,model.config.seed+4800)
        cm=_metric_row(y[selected],candidate[selected],binary);fm=_metric_row(y[selected],fallback[selected],binary);rm=_metric_row(y[control],fallback[control],binary)
        primary.append({
            "purpose":purpose,"selector":route,"coverage_target":q,"selected_N":int(selected.sum()),
            "candidate_minus_fallback_loss":cm["primary_loss"]-fm["primary_loss"],
            "selector_minus_control_loss":fm["primary_loss"]-rm["primary_loss"],
            **boot,
            "bootstrap_unit":str(getattr(model.config,"grouping_status","individual_unverified")),
            "inference":"conditional on frozen model, ABM routes, thresholds and current exploratory roster",
        })
    allmask=known&qc
    # LOCKPAPER-style mechanistic contrasts: same people, same endpoint, group-bootstrap.
    cleaning_specs=[
        ("genetic_twin_library_cleaning","ABM_twin_clean","ABM_twin_full"),
        ("integrated_twin_library_cleaning","ABM_twin_integrated_clean","ABM_twin_integrated_full"),
        ("genetic_vs_integrated_twin","ABM_twin_clean","ABM_twin_integrated_clean"),
        ("global_training_cleaning","ABM_clean","ABM_full"),
        ("pruning_specificity","ABM_clean","ABM_random"),
    ]
    cleaning=[]
    for j,(label,first,second) in enumerate(cleaning_specs):
        if first not in p or second not in p:continue
        row=_bootstrap_same_subset_delta(y,p[first],p[second],allmask,np.asarray(prediction["group"],str),B,model.config.seed+5200+j)
        row.update(contrast=label,model=first,reference=second,subset="all_technical_QC",metric="Brier" if binary else "MSE",direction="negative_favours_model",bootstrap_unit=str(getattr(model.config,"grouping_status","individual_unverified")))
        cleaning.append(row)
    all_candidate=_metric_row(y[allmask],candidate[allmask],binary);all_fallback=_metric_row(y[allmask],fallback[allmask],binary);all_policy=_metric_row(y[allmask],policy[allmask],binary)
    status={
        "version":VERSION,"candidate":ext.candidate,"fallback":ext.fallback,
        "policy_selector":cfg.policy_selector,"predictability_selector":cfg.predictability_selector,
        "primary_coverage":cfg.primary_coverage,"known_outcome_N":int(known.sum()),"technical_QC_N":int((known&qc).sum()),
        "matching_supported_N":int(eligible.sum()),"policy_selected_N":int(np.sum(known&prediction["selected_policy"])),
        "predictability_selected_N":int(np.sum(known&prediction["selected_predictability"])),
        "all_candidate_primary_loss":all_candidate["primary_loss"],"all_fallback_primary_loss":all_fallback["primary_loss"],
        "all_policy_primary_loss":all_policy["primary_loss"],
        "policy_gain_vs_fallback":all_fallback["primary_loss"]-all_policy["primary_loss"],
        "test_outcomes_used_for_fit_or_selection":False,
        **publication_status(model,cfg,external_validation=False),
    }
    return {
        "selector_curve":curve,
        "selector_summary":pd.DataFrame(summaries),
        "primary_contrasts":pd.DataFrame(primary),
        "route_tuning":pd.DataFrame(ext.tuning),
        "reference_cleaning":pd.DataFrame(cleaning),
        "abm_status":pd.DataFrame([status]),
    }
