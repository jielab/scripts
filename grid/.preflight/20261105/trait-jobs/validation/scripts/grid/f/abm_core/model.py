"""Frozen ABM assessment and routing, independent of predictor and device."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .heads import CountHeads, FrozenRule, assert_disjoint, boolean_values
from .losses import assessment_loss

IDENTITY_FIELDS = ("baseline_id", "candidate_id", "bank_sha256", "geometry_sha256",
                   "predictor_bundle_sha256", "role_manifest_sha256")
COVERAGES = (0.2, 0.4, 0.5, 0.6, 0.8, 1.0)


def coverage_tag(coverage):
    return f"c{int(round(100 * coverage)):03d}"


def _jsonable(value):
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(x) for x in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        if np.isnan(value):
            return None
        if np.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(_jsonable(value), ensure_ascii=False,
                                     indent=2, allow_nan=False) + "\n", encoding="utf-8")


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check_manifest(manifest):
    missing = set(IDENTITY_FIELDS) - set(manifest)
    if missing or any(not str(manifest.get(k, "")).strip() for k in IDENTITY_FIELDS):
        raise ValueError(f"Frozen input identity fields required: {sorted(missing)}")
    for key in IDENTITY_FIELDS:
        if key.endswith("sha256"):
            v = str(manifest[key])
            if len(v) != 64 or any(c not in "0123456789abcdef" for c in v):
                raise ValueError(f"{key} must be the actual SHA256 digest of the frozen input")


@dataclass
class ABMAssessor:
    heads: CountHeads
    rules: dict
    manifest: dict
    coverages: tuple = COVERAGES
    primary_coverage: float = 0.5
    min_gain: float = 0.0

    @classmethod
    def fit(cls, gate, head_calibration, policy_calibration, *, manifest,
            task="continuous", loss_name=None, head_options=None,
            coverages=COVERAGES, primary_coverage=0.5, min_gain=0.0):
        _check_manifest(manifest)
        loss_name = assessment_loss(task, loss_name)
        if not np.isfinite(min_gain) or min_gain < 0:
            raise ValueError("min_gain must be finite and nonnegative")
        for first, second in ((gate, head_calibration), (gate, policy_calibration),
                              (head_calibration, policy_calibration)):
            assert_disjoint(first, second)
        if tuple(coverages) != COVERAGES or primary_coverage not in coverages:
            raise ValueError("The six prespecified coverages and primary 50% operating point are required")
        if primary_coverage != 0.5:
            raise ValueError("Primary coverage is frozen at 0.50 for this study version")
        options = head_options or {}
        if not (boolean_values(gate.supported) & np.isfinite(gate.y.to_numpy(float))).any() or not (
                boolean_values(head_calibration.supported) & np.isfinite(head_calibration.y.to_numpy(float))).any():
            heads = CountHeads.unavailable(task, loss_name=loss_name, **{k: options[k] for k in
                ("seed", "min_calibration_groups", "n_boot") if k in options})
        else:
            heads = CountHeads.fit(gate, head_calibration, task=task,
                                   loss_name=loss_name, **options)
        # Policy-calibration labels are deliberately not read to set coverage.
        # Their held-out policy validation is a separate reported evaluation.
        features = policy_calibration.drop(columns=["y"], errors="ignore")
        assessed = heads.predict(features)
        rules = {}
        for objective in ("predictability", "benefit"):
            for coverage in coverages:
                key = f"{objective}_{coverage_tag(coverage)}"
                rules[key] = FrozenRule.fit(assessed, objective=objective,
                                           coverage=coverage, min_gain=min_gain,
                                           seed=heads.seed)
        return cls(heads, rules, dict(manifest, loss_name=loss_name,
                   matching="frozen_risk_bin_disjoint_group_exclusion", ancestry_use="report_label_only"),
                   tuple(coverages), primary_coverage, min_gain)

    def assess(self, features, *, manifest):
        """Label-free assessment. Caller must join by IDs, never by row position."""
        _check_manifest(manifest)
        for key in IDENTITY_FIELDS:
            if manifest[key] != self.manifest[key]:
                raise ValueError(f"Frozen input changed: {key}. Recalibrate ABM for the new bank/model.")
        out = self.heads.predict(features)
        if not {"p_base", "p_candidate"}.issubset(out.columns):
            raise ValueError("Routing requires aligned frozen base and candidate predictions")
        base = out.p_base.to_numpy(float)
        candidate = out.p_candidate.to_numpy(float)
        if not np.isfinite(base).all():
            raise ValueError("The deployable baseline must cover every query")
        if self.heads.task == "binary":
            for p in (base, candidate):
                if ((p[np.isfinite(p)] < 0) | (p[np.isfinite(p)] > 1)).any():
                    raise ValueError("Binary inputs must be calibrated probabilities")
        for key, rule in self.rules.items():
            out[f"selected_{key}"] = rule.apply(out)
        primary = out[f"selected_benefit_{coverage_tag(self.primary_coverage)}"].to_numpy(bool)
        switched = primary & np.isfinite(candidate)
        out["selected"] = primary
        out["switched"] = switched
        out["p_policy"] = np.where(switched, candidate, base)
        out["policy_route"] = np.where(switched, self.manifest["candidate_id"], self.manifest["baseline_id"])
        out["policy_reason"] = np.where(switched, "positive_expected_gain_and_frozen_threshold",
                                         "baseline_by_frozen_benefit_rule")
        out.loc[primary & ~np.isfinite(candidate), "policy_reason"] = "candidate_unavailable_baseline_fallback"
        out.loc[~out.assessment_available, "policy_reason"] = "assessment_unavailable_baseline_fallback"
        out["control_eligible"] = boolean_values(out.supported) & boolean_values(out.assessment_available)
        # Matching strata may be assembled without labels on a target cohort.
        # The risk boundaries are frozen by gate, never recalculated on target.
        # Ancestry is optional reporting metadata and cannot restrict controls.
        risk_bin = np.searchsorted(self.heads.risk_edges, base, side="right")
        out["match_stratum"] = np.char.add("risk", risk_bin.astype(str))
        return out

    def dataset_assessment(self, assessed):
        """Pre-phenotype cohort summary; no inferred R2, AUC or personal guarantee."""
        available = boolean_values(assessed.assessment_available)
        supported = boolean_values(assessed.supported)
        mean_loss = float(assessed.loc[available, "expected_loss_base"].mean()) if available.any() else None
        summary = {
            "n_query": len(assessed), "n_genetically_supported": int(supported.sum()),
            "n_assessment_available": int(available.sum()),
            "assessment_coverage": float(available.mean()) if len(available) else None,
            "expected_mean_loss_assessable_subset": mean_loss,
            "full_dataset_expected_loss": mean_loss if available.all() and len(available) else None,
            "full_dataset_status": "all_assessable" if available.all() and len(available)
                                   else "not_identified_for_unassessed_queries",
            "mean_expected_gain_assessable_subset": float(assessed.loc[available, "expected_gain"].mean())
                                                    if available.any() else None,
            "n_selected_for_candidate": int(assessed.selected.sum()),
            "n_actual_candidate_switches": int(assessed.switched.sum()),
            "loss_unit": self.heads.loss_name,
            "interval_scope": "Calibration-bin means conditional on frozen predictors and bank; not personal error bounds",
            "transport_assumption": "Count-to-loss/gain calibration remains applicable to this target cohort",
            "bank_sha256": self.manifest["bank_sha256"],
            "count_distribution": {str(q): float(assessed.n_twins.quantile(q))
                                   if assessed.n_twins.notna().any() else None
                                   for q in (0, .1, .25, .5, .75, .9, 1)},
            "support_reasons": assessed.support_reason.value_counts().to_dict()
                               if "support_reason" in assessed else {},
            "groups_verified": bool(boolean_values(assessed.groups_verified).all())
                               if "groups_verified" in assessed and len(assessed) else False,
        }
        if self.heads.task == "continuous" and mean_loss is not None:
            summary["expected_rmse_assessable_subset"] = float(np.sqrt(mean_loss))
        return summary

    def save(self, path):
        write_json(path, {"schema_version": "abm.assessor.2", "heads": self.heads.to_dict(),
                          "rules": {k: vars(v) for k, v in self.rules.items()},
                          "manifest": self.manifest, "coverages": self.coverages,
                          "primary_coverage": self.primary_coverage, "min_gain": self.min_gain})

    @classmethod
    def load(cls, path):
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        if d.get("schema_version") != "abm.assessor.2":
            raise ValueError("Unsupported ABM model schema; rerun ABM for pooled risk-only matching")
        if "loss_name" not in d["heads"]:
            raise ValueError("Assessor schema requires its explicit learning loss")
        assessment_loss(d["heads"]["task"], d["heads"]["loss_name"])
        if d["manifest"].get("loss_name", d["heads"]["loss_name"]) != d["heads"]["loss_name"]:
            raise ValueError("Assessor manifest/head loss mismatch")
        rules = {}
        for key, value in d["rules"].items():
            value["threshold"] = float(value["threshold"])
            rules[key] = FrozenRule(**value)
        return cls(CountHeads(**d["heads"]), rules, d["manifest"], tuple(d["coverages"]),
                   d["primary_coverage"], d["min_gain"])
