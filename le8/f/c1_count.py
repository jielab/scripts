"""Transparent molecular twin-count risk, with no fitted regression prediction head.

Only build outcomes enter the donor bank and incidence prior. Build features learn
the standardized molecular PCA space; tune_model FEATURES freeze the radius.
Neither query truth nor tune_gate/calibration outcomes enter risk or support rules.
The existing fixed-horizon IPCW evaluator is used without changing its denominator.
"""
from __future__ import annotations

from contextlib import ExitStack
import csv
import gzip
import json
import os
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

import abm
from c1_models import _subset, finish_projection

VERSION = "LE8-transparent-twin-count-2026-11-02-v1"
DONOR_FIELDS = ["query_eid", "query_group", "donor_rank", "donor_eid", "donor_group", "distance",
                "kernel_weight", "ipcw_weight", "normalized_weight", "donor_outcome"]


def _abstain(risk, diagnostics, qc, prior, cfg):
    diagnostics["technical_qc"] = np.asarray(qc, bool).copy()
    diagnostics["supported"] &= qc
    supported = diagnostics["supported"]
    risk[~supported] = prior
    diagnostics["prior_event_rate"] = np.full(len(risk), prior)
    diagnostics["prior_strength"] = np.full(len(risk), cfg.prior_strength)
    diagnostics["effective_borrow_fraction"] = np.where(supported, diagnostics["borrow_fraction"], 0.)
    diagnostics["abstained"] = ~supported
    diagnostics["risk_for_decision"] = np.where(supported, risk, np.nan)
    diagnostics["abstention_code"] = np.select(
        [~qc, diagnostics["matched_count"] == 0, diagnostics["matched_count"] < cfg.min_twins,
         diagnostics["family_ess"] < cfg.min_family_ess], [1, 2, 3, 4], default=0)
    return risk, diagnostics


def predict_frozen_count(model_directory, X, ids, groups, feature_names, output_directory):
    """Score new individuals from a local fitted bank; there is no outcome input.

    Writes exact donor evidence and a query table below an explicitly LOCAL_ONLY
    folder. Model files are trusted local Python objects; do not load strangers'
    joblib files. All schema and family/self exclusions remain enforced.
    """
    model_directory = Path(model_directory)
    manifest = json.loads((model_directory / "MODEL_FROZEN.json").read_text())
    if manifest.get("method") != "abm_count" or manifest.get("query_regression_head") is not False:
        raise ValueError("A frozen transparent abm_count model is required")
    fitted = joblib.load(model_directory / "model_LOCAL_ONLY.joblib")
    if list(map(str, feature_names)) != fitted["feature_names"]:
        raise ValueError("Query feature names/order differ from the frozen model")
    X, ids, groups = abm.TwinABM.validate_features(X, ids, groups)
    cfg = abm.Config(**fitted["config"])
    private = Path(output_directory) / "count_projection_LOCAL_ONLY"
    private.mkdir(parents=True, mode=0o700, exist_ok=True)
    private.chmod(0o700)
    with gzip.open(private / "donors.tsv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DONOR_FIELDS, delimiter="\t")
        writer.writeheader()

        def evidence(i, row):
            for rank in range(len(row["donor_eid"])):
                writer.writerow(dict(query_eid=ids[i], query_group=groups[i], donor_rank=rank + 1,
                                     **{k: v[rank] for k, v in row.items()}))

        risk, d = fitted["bank"].predict(fitted["prep"].transform(X), ids, groups,
            np.full(len(ids), fitted["prior"]), evidence=evidence)
    risk, d = _abstain(risk, d, fitted["prep"].technical_qc(X), fitted["prior"], cfg)
    result = fitted["assessment"].project(risk, d, ids, groups)
    result["metadata"] = dict(method="abm_count", frozen_source=str(model_directory), query_outcomes_read=False)
    joblib.dump(result, private / "predictions.joblib", compress=3)
    pd.DataFrame(dict(eid=ids, group=groups, probability=risk, **d)).to_csv(
        private / "queries.tsv.gz", sep="\t", index=False)
    for path in private.iterdir():
        if path.is_file():
            path.chmod(0o600)
    return result


class GeometryAssessment:
    """Outcome-blind support ranking, explicitly not a correctness probability."""
    status = "geometry_only_no_error_regression"
    gate_N = 0
    gate_cases = 0
    definition = ("log1p(eligible_match_count) + log1p(family_ess); geometric support only; "
                  "expected error and correctness probability are not estimated")
    calibration_role = "calibration_fit features only; no outcome calibration of count risk or support"

    @staticmethod
    def score(d):
        return np.log1p(d["eligible_match_count"]) + np.log1p(d["family_ess"])

    def fit(self, config, probability, diagnostics, ids):
        score = self.score(diagnostics)
        eligible = diagnostics["technical_qc"] & diagnostics["supported"]
        self.rules = {float(q): abm.FrozenRule().fit(score, eligible, ids, q, config.seed + 1601)
                      for q in config.coverages}
        qc = diagnostics["technical_qc"] & np.isfinite(probability)
        self.risk_edges = (np.unique(np.quantile(probability[qc], np.linspace(0, 1, 11)))[1:-1]
                           if qc.any() else np.array([], float))
        return self

    def project(self, p, diagnostics, ids, groups):
        score = self.score(diagnostics)
        qc = np.asarray(diagnostics["technical_qc"], bool)
        supported = qc & diagnostics["supported"]
        missing = np.full(len(p), np.nan)
        return dict(eid=np.asarray(ids).astype(str), group=np.asarray(groups).astype(str),
            probability=np.asarray(p), technical_qc=qc, supported=supported,
            masks={q: rule.apply(score, supported, ids) for q, rule in self.rules.items()},
            predicted_error=missing.copy(), predictability=missing.copy(),
            residual_predictability=score, predictability_scale="geometry_support_not_error_probability",
            predictability_is_correctness_probability=False, diagnostics=diagnostics,
            risk_edges=self.risk_edges, risk_edges_source="calibration_fit_features_only",
            assessment_status=self.status)


def fit_count(c, out):
    """Fit a frozen donor bank; query risk is a completely reconstructible rate.

    p_local = sum(kernel * IPCW * donor_y) / sum(kernel * IPCW)
    p_shrunk = (family_ESS*p_local + prior_strength*p_build)/(family_ESS+prior_strength)
    Unsupported/QC-failed people explicitly abstain. Their all-person evaluation
    probability is p_build, and risk_for_decision is missing. A supported flag is
    evidence availability, never a guarantee that an individual outcome is known.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    private = out / "donor_evidence_LOCAL_ONLY"
    private.mkdir(mode=0o700)
    b = c.parts["build"]
    donors = b[c.w[b] > 0]
    if not len(donors):
        raise ValueError("No observed, technically valid build donors")
    columns = np.flatnonzero(~np.isin(c.prep.keep, c.cfg.always_include_indices))
    if not len(columns):
        raise ValueError("Twin-count matching requires molecular features")
    geometry = abm.FullProfileGeometry().fit(c.z[b][c.qc[b]], columns, c.cfg)
    bank = abm.TwinBank().fit(c.z[donors], c.y[donors], c.w[donors], c.ids[donors],
                            c.groups[donors], np.zeros(len(donors)), c.cfg, geometry)
    tune = c.parts["tune_model"]
    tune = tune[c.qc[tune]]
    bank.calibrate_radius(c.z[tune], c.ids[tune], c.groups[tune])
    prior = float(np.average(c.y[donors], weights=c.w[donors]))
    prior = float(np.clip(prior, 1e-7, 1 - 1e-7))
    with ExitStack() as stack:
        writers = {}
        for role in ("calibration_audit", "test"):
            handle = stack.enter_context(gzip.open(private / f"{role}_donors.tsv.gz", "wt", newline=""))
            writers[role] = csv.DictWriter(handle, fieldnames=DONOR_FIELDS, delimiter="\t")
            writers[role].writeheader()

        def evidence(i, row):
            writer = writers.get(c.roles[i])
            if writer is not None:
                for rank in range(len(row["donor_eid"])):
                    writer.writerow(dict(query_eid=c.ids[i], query_group=c.groups[i], donor_rank=rank + 1,
                                         **{k: v[rank] for k, v in row.items()}))

        risk, diagnostics = bank.predict(c.z, c.ids, c.groups, np.full(len(c.ids), prior), evidence=evidence)
    # Do not quietly use a thin neighborhood as if it supplied an ABM prediction.
    risk, diagnostics = _abstain(risk, diagnostics, c.qc, prior, c.cfg)
    cal = c.parts["calibration_fit"]
    assessment = GeometryAssessment().fit(c.cfg, risk[cal], _subset(diagnostics, cal), c.ids[cal])
    metadata = dict(method="abm_count", engine="transparent frozen twin counts", count_model_version=VERSION,
        actual_device="cpu", query_regression_head=False, probability_calibrator=False,
        geometry="build-only standardized molecular PCA; clinical columns excluded",
        donor_roles=["build"], fit_build_N=len(donors), build_weighted_incidence=prior,
        donor_case_N=int(c.y[donors].sum()), donor_control_N=int((c.y[donors] == 0).sum()),
        query_formula="(family_ess*local_event_rate + prior_strength*prior_event_rate)/(family_ess+prior_strength)",
        local_rate_formula="sum(kernel_weight*ipcw_weight*donor_outcome)/sum(kernel_weight*ipcw_weight)",
        probability_clip=[1e-7, 1 - 1e-7],
        neighbor_weighting=c.cfg.neighbor_weighting,
        kernel_formula=("1 (simple donor counts with explicit IPCW)" if c.cfg.neighbor_weighting == "uniform"
                        else "exp(-0.5*(distance/frozen_radius)^2)"),
        self_and_family_excluded=True, case_loss_pruning=False,
        unsupported_behavior="abstain; risk_for_decision=NA; all-person evaluation retains build incidence prior",
        abstention_codes={"0": "supported", "1": "technical_qc_failed", "2": "no_eligible_donors",
                          "3": "too_few_donors", "4": "insufficient_independent_family_ess"},
        evidence_directory=private.name, evidence_privacy="participant-level LOCAL_ONLY; never public report",
        interpretation="donor outcomes explain the estimate, not a causal effect or a guaranteed future outcome",
        expected_error_estimated=False)
    saved = finish_projection(c, risk, diagnostics, out,
        fitted=dict(prep=c.prep, km=c.km, bank=bank, geometry=geometry, prior=prior),
        metadata=metadata, assessment=assessment)
    for role in ("calibration_audit", "test"):
        ix = c.parts[role]
        query = pd.DataFrame(dict(eid=c.ids[ix], group=c.groups[ix], probability=risk[ix],
                                  **_subset(diagnostics, ix)))
        query.to_csv(private / f"{role}_queries.tsv.gz", sep="\t", index=False)
    for path in out.rglob("*"):
        if path.is_file() and ("LOCAL_ONLY" in str(path.relative_to(out))):
            path.chmod(0o600)
    return saved


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description="Project frozen transparent twin counts; query NPZ contains features/IDs/groups only")
    p.add_argument("--model-dir", type=Path, required=True)
    p.add_argument("--input", type=Path, required=True, help="NPZ: X, eid, group, feature_names; no outcome columns")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)
    if os.getenv("submit") == "hpc" and not os.getenv("LSB_JOBID"):
        raise RuntimeError("Frozen projection on HPC requires an LSF allocation")
    os.umask(0o077)
    query = abm.read_npz(args.input)
    if set(query) != {"X", "eid", "group", "feature_names"}:
        raise ValueError("Query NPZ must contain exactly X, eid, group, feature_names; truth is not accepted")
    result = predict_frozen_count(args.model_dir, query["X"], query["eid"], query["group"],
                                  query["feature_names"], args.out)
    print(json.dumps(dict(N=len(result["eid"]), supported_N=int(result["supported"].sum()),
        abstained_N=int((~result["supported"]).sum()), output=str(args.out / "count_projection_LOCAL_ONLY"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
