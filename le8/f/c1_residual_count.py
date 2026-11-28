"""Clinical background plus transparent molecular-neighbor residual counts.

This is a separate, named arm. It does not change the pure ``abm`` estimator.
The only outcome contribution is the IPCW mean of build donors' honest clinical
OOF residuals. Molecular proximity, family exclusions, shrinkage and every
donor's contribution remain inspectable. Unsupported queries use clinical risk.
"""
from __future__ import annotations

from contextlib import ExitStack
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances

import abm
from c1_clinical import clinical_background
from c1_count import frozen_count_bank, GeometryAssessment
from c1_models import _subset, finish_projection

VERSION = "LE8-C1-clinical-residual-count-2026-11-16-v1"
RADIUS_QUANTILES = (.50, .75, .90)
CORRECTION_STRENGTHS = (0., .25, .50, 1.)
DONOR_FIELDS = ("query_eid", "query_group", "donor_rank", "donor_eid", "donor_group",
    "distance", "ipcw_weight", "normalized_weight", "donor_outcome",
    "donor_clinical_oof_probability", "donor_residual")


class ResidualCountBank:
    def __init__(self, bank, donor_oof_probability):
        self.bank = bank
        self.oof = np.asarray(donor_oof_probability, float)
        if self.oof.shape != bank.y.shape or not np.isfinite(self.oof).all():
            raise ValueError("Every count donor needs an aligned clinical OOF probability")
        if np.any((self.oof < 0) | (self.oof > 1)):
            raise ValueError("Invalid clinical donor probability")
        if bank.config.neighbor_weighting != "uniform":
            raise ValueError("The declared residual-count arm uses uniform IPCW counts")
        self.residual = bank.y - self.oof
        self.brier = self.residual ** 2
        self.brier_excess = self.brier - self.oof * (1-self.oof)
        pp = np.clip(self.oof, 1e-7, 1-1e-7)
        entropy = -pp*np.log(pp) - (1-pp)*np.log1p(-pp)
        self.logloss_excess = abm.logloss_rows(bank.y, pp) - entropy

    def project(self, z, ids, groups, radii, *, evidence=None):
        """One molecular projection and one radius query for all declared radii.

        Output depends on the query and the frozen bank, never other queries in
        the batch. ``evidence`` is allowed only for a single already chosen radius.
        """
        bank = self.bank
        radii = tuple(float(r) for r in radii)
        if not radii or any(not np.isfinite(r) or r <= 0 for r in radii):
            raise ValueError("Positive frozen molecular radii are required")
        if evidence is not None and len(radii) != 1:
            raise ValueError("Donor evidence is written only for the frozen chosen radius")
        ids, groups = np.asarray(ids).astype(str), np.asarray(groups).astype(str)
        zz = bank.pca.transform(z) / np.sqrt(bank.pca.n_components_)
        if len(ids) != len(zz) or groups.shape != ids.shape:
            raise ValueError("Misaligned molecular query roster")
        n = len(ids)
        fields = ("eligible_match_count", "matched_count", "eligible_reference_count",
            "twin_fraction", "family_count", "family_ess", "donor_ess", "matching_radius",
            "radius_ratio", "nearest_radius_ratio", "borrow_fraction", "matched_case_count",
            "matched_control_count", "weighted_case_mass", "weighted_total_mass",
            "local_event_rate", "local_clinical_oof_probability", "local_clinical_residual",
            "local_clinical_oof_brier", "local_clinical_oof_brier_excess",
            "local_clinical_oof_logloss_excess")
        result = [{name: np.zeros(n) for name in fields} for _ in radii]
        levels, count = np.unique(bank.groups, return_counts=True)
        family_size = dict(zip(levels.tolist(), count.tolist()))
        donor_ids = dict(zip(bank.ids.tolist(), bank.groups.tolist()))
        denominator = np.asarray([len(bank.y) - family_size.get(g, 0) -
            int(i in donor_ids and donor_ids[i] != g) for i,g in zip(ids,groups)])
        for radius, d in zip(radii, result):
            d["matching_radius"].fill(radius)
            d["radius_ratio"].fill(100.)
            d["nearest_radius_ratio"].fill(100.)
            d["eligible_reference_count"] = denominator.astype(float)
        for start in range(0, n, 128):
            query = zz[start:start+128]
            if bank.index is not None:
                near, distances = bank.index.query_radius(query, r=max(radii), return_distance=True)
            else:
                dense = euclidean_distances(query, bank.Z)
                near = [np.flatnonzero(dd <= max(radii)) for dd in dense]
                distances = [dd[ix] for dd,ix in zip(dense,near)]
            for offset, (nn, ds) in enumerate(zip(near, distances)):
                i = start+offset
                good = (bank.ids[nn] != ids[i]) & (bank.groups[nn] != groups[i]) & (bank.w[nn] > 0)
                nn, ds = nn[good], ds[good]
                order = np.lexsort((bank.ids[nn], ds))
                nn, ds = nn[order], ds[order]
                for radius, d in zip(radii, result):
                    k = np.searchsorted(ds, radius, side="right")
                    if not k:
                        continue
                    idx, distance = nn[:k], ds[:k]
                    raw = bank.w[idx]
                    weight = raw/raw.sum()
                    _, inv = np.unique(bank.groups[idx], return_inverse=True)
                    mass = np.bincount(inv, weights=weight)
                    fess = 1./np.sum(mass**2)
                    values = dict(eligible_match_count=k, matched_count=k,
                        twin_fraction=k/denominator[i] if denominator[i] > 0 else 0.,
                        family_count=len(mass), family_ess=fess, donor_ess=1./np.sum(weight**2),
                        radius_ratio=distance[-1]/radius, nearest_radius_ratio=distance[0]/radius,
                        borrow_fraction=fess/(fess+bank.config.prior_strength),
                        matched_case_count=float(bank.y[idx].sum()),
                        matched_control_count=float(np.sum(bank.y[idx] == 0)),
                        weighted_case_mass=float(raw@bank.y[idx]), weighted_total_mass=float(raw.sum()),
                        local_event_rate=float(weight@bank.y[idx]),
                        local_clinical_oof_probability=float(weight@self.oof[idx]),
                        local_clinical_residual=float(weight@self.residual[idx]),
                        local_clinical_oof_brier=float(weight@self.brier[idx]),
                        local_clinical_oof_brier_excess=float(weight@self.brier_excess[idx]),
                        local_clinical_oof_logloss_excess=float(weight@self.logloss_excess[idx]))
                    for key,value in values.items():
                        d[key][i] = value
                    if evidence is not None:
                        evidence(i, dict(donor_eid=bank.ids[idx], donor_group=bank.groups[idx],
                            distance=distance, ipcw_weight=raw, normalized_weight=weight,
                            donor_outcome=bank.y[idx], donor_clinical_oof_probability=self.oof[idx],
                            donor_residual=self.residual[idx]))
        for d in result:
            d["supported"] = ((d["matched_count"] >= bank.config.min_twins) &
                (d["family_ess"] >= bank.config.min_family_ess))
        return result


def corrected_probability(baseline, diagnostics, strength, qc):
    baseline = np.asarray(baseline, float)
    supported = np.asarray(qc, bool) & np.asarray(diagnostics["supported"], bool)
    correction = float(strength) * diagnostics["borrow_fraction"] * diagnostics["local_clinical_residual"]
    out=baseline.copy()
    if float(strength)!=0.:
        out[supported]=np.clip(baseline[supported]+correction[supported],1e-7,1-1e-7)
    return out


def fit_clinical_count(c, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    private = out/"tmp/donor_evidence"; private.mkdir(parents=True, mode=0o700)
    clinical = clinical_background(c, require_oof=True)
    base, geometry, donors, _ = frozen_count_bank(c)
    if not np.array_equal(donors, clinical.donors):
        raise AssertionError("Clinical OOF and molecular count donors differ")
    bank = ResidualCountBank(base, clinical.oof_probability)
    radii = tuple(float(np.quantile(base.radius_reference_distances, q)) for q in RADIUS_QUANTILES)
    tune = c.parts["tune_model"]
    tune_diagnostics = bank.project(c.z[tune], c.ids[tune], c.groups[tune], radii)
    trials, best = [], None
    # Zero correction comes first and is exactly the calibrated clinical probability.
    # It can win: the pipeline must not manufacture a molecular increment.
    for strength in CORRECTION_STRENGTHS:
        for quantile,radius,d in zip(RADIUS_QUANTILES, radii, tune_diagnostics):
            pred = corrected_probability(clinical.probability[tune], d, strength, c.qc[tune])
            use = c.w[tune] > 0
            loss = float(np.average(abm.logloss_rows(c.y[tune][use], pred[use]), weights=c.w[tune][use]))
            trials.append(dict(radius_quantile=quantile, radius=radius, correction_strength=strength,
                tune_logloss=loss, tune_supported_fraction=float(np.mean(d["supported"] & c.qc[tune]))))
            if best is None or loss < best[0]:
                best = loss, quantile, radius, strength
    _, quantile, radius, strength = best
    with ExitStack() as stack:
        writers = {}
        for role in ("calibration_audit", "test"):
            handle = stack.enter_context(gzip.open(private/f"{role}_donors.tsv.gz", "wt", newline=""))
            writers[role] = csv.DictWriter(handle, fieldnames=DONOR_FIELDS, delimiter="\t")
            writers[role].writeheader()
        def evidence(i, row):
            writer = writers.get(c.roles[i])
            if writer is not None:
                for rank in range(len(row["donor_eid"])):
                    writer.writerow(dict(query_eid=c.ids[i], query_group=c.groups[i], donor_rank=rank+1,
                        **{k:v[rank] for k,v in row.items()}))
        query=np.flatnonzero(c.roles!="build")
        local=bank.project(c.z[query],c.ids[query],c.groups[query],(radius,),
            evidence=lambda i,row:evidence(int(query[i]),row))[0]
        d={key:np.zeros(len(c.ids),dtype=np.asarray(value).dtype) for key,value in local.items()}
        for key,value in local.items():d[key][query]=value
    d["technical_qc"] = c.qc.copy()
    d["supported"] &= c.qc
    probability = corrected_probability(clinical.probability, d, strength, c.qc)
    probability[~d["supported"]] = clinical.probability[~d["supported"]]
    d.update(clinical_background_probability=clinical.probability.copy(),
        clinical_background_raw_probability=clinical.raw_probability.copy(),
        correction_strength=np.full(len(c.ids), strength),
        effective_residual_correction=np.where(d["supported"],
            strength*d["borrow_fraction"]*d["local_clinical_residual"], 0.),
        clinical_fallback=~d["supported"],
        molecular_increment_used=d["supported"] & (strength > 0))
    cal = c.parts["calibration_fit"]
    assessment = GeometryAssessment().fit(c.cfg, probability[cal], _subset(d,cal), c.ids[cal])
    metadata = dict(method="abm_clinical", engine="clinical background plus molecular residual counts",
        actual_device="cpu", version=VERSION, clinical_columns_shared_with_comparators=True,
        clinical_columns=[c.names[c.prep.keep[j]] for j in clinical.columns],
        molecular_geometry="same build-only PCA as pure abm; clinical columns excluded",
        donor_roles=["build"], fit_build_N=len(donors), donor_case_N=int(c.y[donors].sum()),
        donor_control_N=int(np.sum(c.y[donors] == 0)), case_loss_pruning=False,
        donor_residual_definition="incident horizon_y minus calibrated clinical family-OOF probability",
        historical_prevalent_outcomes_used=False, donor_OOF_audit=clinical.oof_audit,
        final_query_formula="clip(p_clinical_calibrated + alpha*family_ESS/(family_ESS+prior_strength)*IPCW_mean(donor_y-donor_clinical_calibrated_OOF),1e-7,1-1e-7); unsupported=clinical",
        probability_calibrator="clinical baseline and each OOF fold calibrated separately on calibration_fit; no second map after residual correction",
        selected_radius_quantile=quantile, selected_radius=radius, selected_correction_strength=strength,
        zero_correction_is_valid=True, tuning=trials, tuning_role="tune_model only; evaluates exactly the deployed calibrated-scale probability",
        calibration_role="calibration_fit supplies clinical probability maps before residual-strength tuning; reused development role, not independent validation",
        development_tuning_reuse="tune_model selects clinical penalties and residual settings; tuning losses are not validation estimates",
        geometry_settings=len(radii), PCA_fits=1, tune_radius_queries=1,
        unsupported_behavior="exact clinical baseline fallback; retained in all-person evaluation",
        self_and_family_excluded=True, expected_error_estimated=False,
        evidence_directory=str(private.relative_to(out)), evidence_privacy="participant-level LOCAL_ONLY")
    saved = finish_projection(c, probability, d, out,
        fitted=dict(prep=c.prep, km=c.km, clinical=clinical, bank=bank, geometry=geometry,
            radius=radius, correction_strength=strength), metadata=metadata, assessment=assessment)
    pd.DataFrame(dict(eid=c.ids[donors], group=c.groups[donors], incident_y=c.y[donors],
        clinical_oof_probability=clinical.oof_probability, oof_fold=clinical.oof_fold,
        donor_residual=c.y[donors]-clinical.oof_probability)).to_csv(
            private/"build_residuals.tsv.gz",sep="\t",index=False)
    for role in ("calibration_audit", "test"):
        ix=c.parts[role]
        pd.DataFrame(dict(eid=c.ids[ix],group=c.groups[ix],probability=probability[ix],
            **_subset(d,ix))).to_csv(private/f"{role}_queries.tsv.gz",sep="\t",index=False)
    for path in out.rglob("*"):
        if path.is_file() and ("LOCAL_ONLY" in path.name or "tmp" in path.relative_to(out).parts):
            path.chmod(0o600)
    return saved


def predict_frozen_clinical_count(fitted, X, ids, groups, feature_names):
    """Outcome-free standalone projection; intentionally accepts no truth input."""
    if list(map(str,feature_names)) != list(map(str,fitted["feature_names"])):
        raise ValueError("Molecular/clinical feature order differs from frozen schema")
    X,ids,groups=abm.TwinABM.validate_features(X,ids,groups)
    z=fitted["prep"].transform(X);qc=fitted["prep"].technical_qc(X)
    clinical=fitted["clinical"]
    d=fitted["bank"].project(z,ids,groups,(fitted["radius"],))[0]
    d["technical_qc"]=qc;d["supported"] &= qc
    baseline=clinical.predict(z)
    p=corrected_probability(baseline,d,fitted["correction_strength"],qc)
    p[~d["supported"]]=baseline[~d["supported"]]
    return dict(eid=ids,group=groups,probability=p,technical_qc=qc,supported=d["supported"],diagnostics=d)
