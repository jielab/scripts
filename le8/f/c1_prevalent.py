#!/usr/bin/env python3
"""Development-only disease-history features for prospective ABM prediction.

This transformer accepts NO event, follow-up time or incident outcome argument.
Prevalent profiles describe a postdiagnosis source domain, not observed future
cases. They provide auxiliary molecular-similarity features; the downstream
risk model, model selection and probability calibration must use incident data.
Age/sex residualization and duration/treatment weighting do not identify causal
effects or reconstruct a person's prediagnosis molecular profile.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import warnings

import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances


@dataclass(frozen=True)
class PrevalentConfig:
    seed: int = 12345
    n_components: int = 12
    n_prototypes: int = 6
    min_prevalent: int = 20
    min_incident_build: int = 30
    max_feature_missing: float = 0.4
    max_sample_missing: float = 0.4
    age_band_years: float = 10.0
    match_sex: bool = True
    age_overlap_margin_years: float = 5.0
    duration_half_weight_years: float = 5.0
    weighting_policy: str = "duration_treatment"
    unknown_duration_weight: float = 0.5
    medication_weight: float = 0.5
    unknown_medication_weight: float = 0.75
    radius_quantile: float = 0.50
    radius_neighbour: int = 5
    max_neighbours: int = 64
    nuisance_ridge: float = 1.0
    chunk_size: int = 256

    def validate(self):
        if min(self.n_components, self.n_prototypes, self.min_prevalent,
               self.min_incident_build, self.radius_neighbour,
               self.max_neighbours, self.chunk_size) < 1:
            raise ValueError("Prevalent configuration counts must be positive")
        if not 0 < self.radius_quantile < 1:
            raise ValueError("radius_quantile must be strictly between zero and one")
        if self.weighting_policy not in {"none", "duration_treatment"}:
            raise ValueError("weighting_policy must be none or duration_treatment")
        if min(self.age_band_years, self.duration_half_weight_years) <= 0:
            raise ValueError("Age band and duration scale must be positive")
        if self.age_overlap_margin_years < 0 or self.nuisance_ridge < 0:
            raise ValueError("Negative overlap margin or nuisance penalty")
        for value in (self.max_feature_missing, self.max_sample_missing,
                      self.unknown_duration_weight, self.medication_weight,
                      self.unknown_medication_weight):
            if not 0 < value <= 1:
                raise ValueError("Missingness thresholds and source weights must be in (0, 1]")


def _vector(value, n, name, dtype=str):
    x = np.asarray(value, dtype=dtype)
    if x.ndim != 1 or len(x) != n:
        raise ValueError(f"{name} must have one entry per row")
    return x


def _metadata(metadata, n):
    metadata = {} if metadata is None else metadata
    if "treatment" in metadata and "medication_exposed" not in metadata:
        metadata = dict(metadata, medication_exposed=metadata["treatment"])
    result = {}
    for name in ("age", "duration_years", "medication_exposed"):
        result[name] = _vector(metadata.get(name, np.full(n, np.nan)), n, name, float)
        if np.isinf(result[name]).any():
            raise ValueError(f"Infinite {name}")
    sex = _vector(metadata.get("sex", np.full(n, "")), n, "sex", object)
    # Numeric labels 0 and 0.0 must agree across R/pandas encodings.
    def normalize(v):
        if v is None:
            return ""
        s = str(v).strip()
        if s.lower() in {"", "nan", "none", "<na>", "na", "unknown", "missing", "unavailable"}:
            return ""
        try:
            f = float(s)
            if np.isfinite(f):
                return str(int(f)) if f.is_integer() else str(f)
        except ValueError:
            pass
        return s
    result["sex"] = np.asarray([normalize(v) for v in sex], str)
    med = result["medication_exposed"]
    if np.any(np.isfinite(med) & ~np.isin(med, [0.0, 1.0])):
        raise ValueError("medication_exposed must be 0, 1, or missing")
    duration = result["duration_years"]
    if np.any(np.isfinite(duration) & (duration < 0)):
        raise ValueError("Prevalent duration cannot be negative")
    return result


def _take(metadata, mask):
    return {key: value[mask] for key, value in metadata.items()}


class PrevalentPrototypeTransformer:
    """Fixed-width source-domain support features, never source-domain risk.

    ``incident_roles`` must include all held-out roles supplied by the caller;
    an external forbidden-family list can cover holdouts held in a separate
    file. Every source row must have explicit ``prevalent_roles == 'build'``.
    Calling code must preserve its family roster for the full cohort, including
    the prevalent domain. Any detected cross-role family is a hard error.

    Missing age/sex prevents history matching for that row, but never removes
    the person from the downstream incident-risk population. Missing duration
    and medication receive disclosed conservative weights.
    """

    def __init__(self, config=None):
        self.config = config or PrevalentConfig()
        self.config.validate()

    @property
    def feature_names(self):
        return ["prevalent::similarity", "prevalent::near_fraction",
                "prevalent::family_ess_fraction"] + [
                    f"prevalent::route_{j + 1}" for j in range(self.config.n_prototypes)
                ]

    @property
    def feature_names_out_(self):
        return self.feature_names

    def _unavailable(self, reason):
        self.available_ = False
        self.audit_["status"] = reason
        self.audit_["source_used_N"] = 0
        return self

    def fit(self, incident_X, incident_ids, incident_groups, incident_roles,
            prevalent_X, prevalent_ids, prevalent_groups, prevalent_roles,
            incident_metadata=None, prevalent_metadata=None, feature_names=None,
            *, forbidden_groups=(), forbidden_ids=()):
        xi, xp = np.asarray(incident_X, float), np.asarray(prevalent_X, float)
        if xi.ndim != 2 or xp.ndim != 2 or xi.shape[1] != xp.shape[1]:
            raise ValueError("Incident and prevalent molecular schemas must match")
        if np.isinf(xi).any() or np.isinf(xp).any():
            raise ValueError("Infinite molecular input")
        ni, npv = len(xi), len(xp)
        ii, ig, ir = [_vector(a, ni, n) for a, n in (
            (incident_ids, "incident_ids"), (incident_groups, "incident_groups"),
            (incident_roles, "incident_roles"))]
        pi, pg, pr = [_vector(a, npv, n) for a, n in (
            (prevalent_ids, "prevalent_ids"), (prevalent_groups, "prevalent_groups"),
            (prevalent_roles, "prevalent_roles"))]
        allowed = {"build", "tune_model", "tune_gate", "calibration_fit",
                   "calibration_audit", "test"}
        if not set(ir).issubset(allowed) or np.any(pr != "build"):
            raise ValueError("All source rows must be explicitly build; unknown incident role")
        for ids, label in ((ii, "incident"), (pi, "prevalent")):
            if len(set(ids)) != len(ids) or np.any(ids == ""):
                raise ValueError(f"Duplicate or empty {label} IDs")
        if np.any(ig == "") or np.any(pg == ""):
            raise ValueError("Family groups must be nonempty")
        family_roles = {}
        for group, role in zip(ig, ir):
            if group in family_roles and family_roles[group] != role:
                raise ValueError("Incident family spans roles")
            family_roles[group] = role
        excluded_groups = set(ig[ir != "build"]) | set(map(str, forbidden_groups))
        excluded_ids = set(ii) | set(map(str, forbidden_ids))
        if set(pg) & excluded_groups or set(pi) & excluded_ids:
            raise ValueError("Prevalent source overlaps incident IDs or held-out families")
        self.input_columns_ = xi.shape[1]
        self.input_feature_names_ = None if feature_names is None else tuple(map(str, feature_names))
        if feature_names is not None and len(feature_names) != self.input_columns_:
            raise ValueError("Molecular feature-name schema length mismatch")
        mi, mp = _metadata(incident_metadata, ni), _metadata(prevalent_metadata, npv)
        build = ir == "build"
        self.audit_ = dict(
            schema="le8-prevalent-transfer-v1", config=asdict(self.config),
            status="initializing", input_build_N=int(build.sum()),
            source_offered_N=npv, source_used_N=0,
            source_roles="build_only", incidence_labels_consumed=0,
            source_can_enter_incident_denominator=False,
            family_exclusion_checked=True,
            source_unknown_duration_N=int(np.isnan(mp["duration_years"]).sum()),
            source_unknown_medication_N=int(np.isnan(mp["medication_exposed"]).sum()),
            adjustment="incident-build age/sex residualization; source duration/treatment weighting",
            source_weighting_policy=self.config.weighting_policy,
            duration_weighting_interpretation="prespecified sensitivity heuristic, not reverse-causation correction",
            causal_interpretation="none; postdiagnosis differences and survivor bias remain",
            feature_names=self.feature_names,
            source_id_sha256=hashlib.sha256("\n".join(sorted(pi)).encode()).hexdigest(),
        )
        self.available_ = False
        if build.sum() < self.config.min_incident_build:
            return self._unavailable("insufficient_incident_build")
        if npv < self.config.min_prevalent:
            return self._unavailable("insufficient_prevalent_build")
        xb = xi[build]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            vari = np.nanvar(xb, axis=0)
        self.keep_ = np.flatnonzero((np.mean(~np.isfinite(xb), axis=0) <= self.config.max_feature_missing)
                                   & np.isfinite(vari) & (vari > 1e-12))
        if len(self.keep_) < 2:
            return self._unavailable("insufficient_variable_molecules")
        build_qc = np.mean(~np.isfinite(xi[:, self.keep_]), axis=1) <= self.config.max_sample_missing
        build &= build_qc & np.isfinite(mi["age"])
        if self.config.match_sex:
            build &= mi["sex"] != ""
        if build.sum() < self.config.min_incident_build:
            return self._unavailable("insufficient_build_metadata_or_assay_support")
        xb, mb = xi[build][:, self.keep_], _take(mi, build)
        self.median_ = np.nanmedian(xb, axis=0)
        xb = np.where(np.isfinite(xb), xb, self.median_)
        self.mean_, self.scale_ = xb.mean(0), np.maximum(xb.std(0), 1e-8)
        self.age_mean_, self.age_scale_ = float(mb["age"].mean()), max(float(mb["age"].std()), 1.0)
        self.age_min_, self.age_max_ = np.quantile(mb["age"], [0.01, 0.99])
        self.age_min_ -= self.config.age_overlap_margin_years
        self.age_max_ += self.config.age_overlap_margin_years
        self.sex_levels_ = tuple(sorted(set(mb["sex"]) - {""}))
        nuisance = self._nuisance(mb)
        z = (xb - self.mean_) / self.scale_
        penalty = np.eye(nuisance.shape[1]) * self.config.nuisance_ridge
        penalty[0, 0] = 0
        self.nuisance_coef_ = np.linalg.solve(nuisance.T @ nuisance + penalty, nuisance.T @ z)
        residual = z - nuisance @ self.nuisance_coef_
        components = min(self.config.n_components, residual.shape[1], len(residual) - 1)
        self.pca_ = PCA(n_components=components, whiten=False, random_state=self.config.seed,
                        svd_solver="randomized" if components < min(residual.shape) else "full")
        zb = self.pca_.fit_transform(residual)
        self.pc_scale_ = np.maximum(zb.std(0), 1e-6)
        zb /= self.pc_scale_
        source_ok = self._metadata_ok(mp) & (
            np.mean(~np.isfinite(xp[:, self.keep_]), axis=1) <= self.config.max_sample_missing)
        self.audit_["source_metadata_or_assay_excluded_N"] = int((~source_ok).sum())
        if source_ok.sum() < self.config.min_prevalent:
            return self._unavailable("insufficient_prevalent_overlap")
        self.source_ids_, self.source_groups_ = pi[source_ok], pg[source_ok]
        self.source_metadata_ = _take(mp, source_ok)
        self.source_z_ = self._project(xp[source_ok], self.source_metadata_)
        duration = self.source_metadata_["duration_years"]
        dw = np.where(np.isfinite(duration), 1.0 / (1.0 + duration / self.config.duration_half_weight_years),
                      self.config.unknown_duration_weight)
        med = self.source_metadata_["medication_exposed"]
        mw = np.where(np.isnan(med), self.config.unknown_medication_weight,
                      np.where(med == 1, self.config.medication_weight, 1.0))
        self.source_weight_ = dw * mw
        if self.config.weighting_policy == "none":
            self.source_weight_[:] = 1.0
        # A deterministic source prototype partition; no test data or outcomes.
        nproto = min(self.config.n_prototypes, len(self.source_z_))
        self.cluster_ = MiniBatchKMeans(n_clusters=nproto, random_state=self.config.seed,
                                       n_init=3, batch_size=min(1024, len(self.source_z_)))
        self.cluster_.fit(self.source_z_, sample_weight=self.source_weight_)
        self.source_route_ = self.cluster_.labels_
        # Radius depends on outcome-blind BUILD-to-source distances, frozen before tune/test.
        kth = []
        for start in range(0, len(zb), self.config.chunk_size):
            stop = min(start + self.config.chunk_size, len(zb))
            dist = euclidean_distances(zb[start:stop], self.source_z_)
            ok = self._candidate_mask(_take(mb, slice(start, stop)), ii[build][start:stop], ig[build][start:stop])
            dist[~ok] = np.inf
            rank = min(self.config.radius_neighbour, dist.shape[1]) - 1
            values = np.partition(dist, rank, axis=1)[:, rank]
            kth.extend(values[np.isfinite(values)].tolist())
        if not kth:
            return self._unavailable("no_build_source_matching_support")
        self.radius_ = max(float(np.quantile(kth, self.config.radius_quantile)), 1e-8)
        self.available_ = True
        self.audit_.update(status="available", source_used_N=len(self.source_ids_),
                           build_transform_N=int(build.sum()), matching_radius=self.radius_,
                           source_weight_sum=float(self.source_weight_.sum()),
                           source_weight_min=float(self.source_weight_.min()),
                           source_weight_max=float(self.source_weight_.max()),
                           source_weight_median=float(np.median(self.source_weight_)),
                           source_used_duration_known_N=int(np.isfinite(duration).sum()),
                           source_used_duration_unknown_N=int(np.isnan(duration).sum()),
                           source_used_medication_exposed_N=int(np.sum(med == 1)),
                           source_used_medication_unexposed_N=int(np.sum(med == 0)),
                           source_used_medication_unknown_N=int(np.isnan(med).sum()),
                           source_used_duration_min_years=float(np.nanmin(duration)) if np.isfinite(duration).any() else None,
                           source_used_duration_max_years=float(np.nanmax(duration)) if np.isfinite(duration).any() else None,
                           source_weight_ESS=float(self.source_weight_.sum() ** 2 / np.square(self.source_weight_).sum()),
                           age_overlap_min=float(self.age_min_), age_overlap_max=float(self.age_max_),
                           molecular_features_N=len(self.keep_), representation_components=components)
        return self

    def _nuisance(self, metadata):
        age = np.where(np.isfinite(metadata["age"]), metadata["age"], self.age_mean_)
        a = (age - self.age_mean_) / self.age_scale_
        # Quadratic age and sex offsets describe domain location, not risk effects.
        columns = [np.ones(len(a)), a, a * a]
        columns += [(metadata["sex"] == s).astype(float) for s in self.sex_levels_[1:]]
        return np.column_stack(columns)

    def _metadata_ok(self, metadata):
        age = metadata["age"]
        ok = np.isfinite(age) & (age >= self.age_min_) & (age <= self.age_max_)
        if self.config.match_sex:
            ok &= np.isin(metadata["sex"], self.sex_levels_)
        return ok

    def _project(self, X, metadata):
        z = np.asarray(X, float)[:, self.keep_]
        z = (np.where(np.isfinite(z), z, self.median_) - self.mean_) / self.scale_
        z -= self._nuisance(metadata) @ self.nuisance_coef_
        return self.pca_.transform(z) / self.pc_scale_

    def _candidate_mask(self, metadata, ids=None, groups=None):
        ok = np.abs(metadata["age"][:, None] - self.source_metadata_["age"][None, :]) <= self.config.age_band_years
        if self.config.match_sex:
            ok &= metadata["sex"][:, None] == self.source_metadata_["sex"][None, :]
        if ids is not None:
            ok &= np.asarray(ids, str)[:, None] != self.source_ids_[None, :]
        if groups is not None:
            ok &= np.asarray(groups, str)[:, None] != self.source_groups_[None, :]
        return ok

    def transform_with_diagnostics(self, X, metadata=None, ids=None, groups=None, *, feature_names=None):
        if not hasattr(self, "available_"):
            raise ValueError("Fit history transformer before transform")
        x = np.asarray(X, float)
        if x.ndim != 2 or x.shape[1] != self.input_columns_ or np.isinf(x).any():
            raise ValueError("History-transform molecular schema mismatch")
        if feature_names is not None and tuple(map(str, feature_names)) != self.input_feature_names_:
            raise ValueError("History-transform molecular names/order mismatch")
        n = len(x)
        ids = None if ids is None else _vector(ids, n, "query_ids")
        groups = None if groups is None else _vector(groups, n, "query_groups")
        metadata = _metadata(metadata, n)
        result = np.zeros((n, len(self.feature_names)))
        diagnostics = {name: np.zeros(n, dtype=float) for name in (
            "prevalent_eligible_reference_count", "prevalent_matched_count",
            "prevalent_used_count", "prevalent_family_ess", "prevalent_support_available")}
        diagnostics["prevalent_matching_radius"] = np.full(n, getattr(self, "radius_", np.nan))
        if not self.available_:
            return result, diagnostics
        z = self._project(x, metadata)
        query_ok = self._metadata_ok(metadata) & (
            np.mean(~np.isfinite(x[:, self.keep_]), axis=1) <= self.config.max_sample_missing)
        for start in range(0, n, self.config.chunk_size):
            stop = min(start + self.config.chunk_size, n)
            candidate = self._candidate_mask(_take(metadata, slice(start, stop)),
                                             None if ids is None else ids[start:stop],
                                             None if groups is None else groups[start:stop])
            candidate &= query_ok[start:stop, None]
            distance = euclidean_distances(z[start:stop], self.source_z_)
            matched = candidate & (distance <= self.radius_)
            for j in range(stop - start):
                row = start + j
                eligible = np.flatnonzero(candidate[j])
                found = np.flatnonzero(matched[j])
                diagnostics["prevalent_eligible_reference_count"][row] = len(eligible)
                diagnostics["prevalent_matched_count"][row] = len(found)
                diagnostics["prevalent_support_available"][row] = bool(query_ok[row] and len(eligible))
                if not len(found):
                    continue
                # Count all matches before K truncation; one source domain only.
                chosen = found[np.argsort(distance[j, found], kind="stable")[:self.config.max_neighbours]]
                weights = self.source_weight_[chosen] * np.exp(-0.5 * np.square(distance[j, chosen] / self.radius_))
                denom = self.source_weight_[eligible].sum()
                mass = float(weights.sum())
                family, inverse = np.unique(self.source_groups_[chosen], return_inverse=True)
                family_weight = np.bincount(inverse, weights=weights)
                ess = mass * mass / float(np.square(family_weight).sum()) if mass else 0.0
                diagnostics["prevalent_used_count"][row] = len(chosen)
                diagnostics["prevalent_family_ess"][row] = ess
                result[row, :3] = [mass / denom if denom else 0.0,
                                   len(found) / len(eligible), ess / max(len(eligible), 1)]
                route_mass = np.bincount(self.source_route_[chosen], weights=weights,
                                         minlength=self.config.n_prototypes)
                result[row, 3:] = route_mass / denom if denom else 0.0
        return result, diagnostics

    def transform(self, X, metadata=None, ids=None, groups=None, *, feature_names=None):
        return self.transform_with_diagnostics(X, metadata, ids, groups, feature_names=feature_names)[0]

    def support_diagnostics(self, X, metadata=None, ids=None, groups=None, *, feature_names=None):
        return self.transform_with_diagnostics(X, metadata, ids, groups, feature_names=feature_names)[1]
