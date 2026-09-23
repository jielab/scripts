"""Independent supervised prot/met -> horizon disease risk workflow."""
from pathlib import Path
import importlib.metadata
import json
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, confusion_matrix
from common import dump, log, words, fingerprints, digest
from io_data import prepare, read_table, numeric, map_metabolites
from preprocess import outcomes, MolecularPreprocessor, MetadataDesign
from pipeline import split_people, unknown_rows
from models import tune_logistic, RiskCalibrator
from survival import CensoringKM, metrics
from tf_model import select_features, train_transformer, predict, restore_predictor

VERSION = "TF-1.0.0"
PRIMARY = "tabicl_finetuned"


def manifest(a, out, model_info):
    inputs = [] if a.demo else [a.phe_file, a.omics_file]
    inputs += [v for v in [a.split_file, a.met_map if a.biom == "met" and a.met_input == "raw" and not a.demo else ""] if v]
    value = dict(version=VERSION, config=vars(a), pretrained=model_info,
        dependencies={name: importlib.metadata.version(name) for name in
                      ["torch", "tabicl", "scikit-learn", "numpy", "pandas"]},
        inputs=fingerprints(inputs, a.full_input_hash),
        code=fingerprints(list(Path(__file__).parent.glob("*.py"))+
                          [Path(__file__).parent.parent/"panome_TF.sh"], True))
    value["signature"] = digest(value)
    dump(out/"manifest.json", value)


def cohort(a, out):
    folder = out/"input"
    folder.mkdir(exist_ok=True)
    p = prepare(a, folder)
    raw = np.load(folder/"raw.npy")
    features = (folder/"features.txt").read_text().splitlines()
    p, audit = outcomes(p, a)
    eligible = p.eligible.to_numpy(bool)
    raw, p = raw[eligible], p.loc[eligible].reset_index(drop=True)
    p["split"] = split_people(p, a)
    prep = MolecularPreprocessor(a.feature_missing, words(a.residualize), words(a.categorical), a.transform)
    build = p.split.eq("build").to_numpy()
    scaled = prep.scale_transform(raw[build])
    keep = (np.mean(~np.isfinite(scaled), axis=0) < a.feature_missing) & (np.nanstd(scaled, axis=0) > 1e-8)
    if keep.sum() < 3:
        raise ValueError("Fewer than three assays pass training QC")
    missing = np.mean(~np.isfinite(prep.scale_transform(raw[:, keep])), axis=1)
    qc = missing <= a.sample_missing
    pd.DataFrame({a.id_col: p[a.id_col], "split": p.split, "missing_fraction": missing,
                  "included": qc}).to_csv(out/"sample_qc.csv", index=False)
    raw, p, missing = raw[qc], p.loc[qc].reset_index(drop=True), missing[qc]
    masks = {name: p.split.eq(name).to_numpy() for name in ["build", "tune", "calibration", "test"]}
    if min(mask.sum() for mask in masks.values()) < 30:
        raise ValueError("Each outer partition needs at least 30 people after QC")
    prep.fit(raw[masks["build"]], p.loc[masks["build"]], allowed=keep)
    x, observed = prep.transform(raw, p)
    clinical = MetadataDesign(words(a.covariates), words(a.categorical)).fit(p.loc[masks["build"]])
    c = clinical.transform(p)
    if c.shape[1] == 0:
        raise ValueError("At least one clinical covariate required for the comparator")
    km = CensoringKM().fit(p.loc[masks["build"]], a.horizon, a.min_censor_survival)
    y, w = np.zeros(len(p), int), np.zeros(len(p))
    development = ~masks["test"]
    y[development], w[development] = km.labels_weights(p.loc[development])
    for name in ["build", "tune", "calibration"]:
        known = masks[name] & (w > 0)
        if min(np.sum(y[known] == label) for label in [0, 1]) < a.min_events:
            raise ValueError(f"Too few known cases/controls in {name}; adjust horizon or --min-events for demos")
    p[[a.id_col, "split"]].to_csv(out/"split.csv", index=False)
    audit.update(after_qc=len(p), retained_features=len(prep.keep), role_counts=p.split.value_counts().to_dict(),
                 outcome_definition=f"incident disease by {a.horizon:g} years; death treated as censoring", family_split=bool(a.group_col))
    dump(out/"cohort_audit.json", audit)
    return p, x, c, observed, y, w, masks, prep, clinical, km, features


def decision_threshold(y, w, risk):
    valid = w > 0
    fpr, tpr, thresholds = roc_curve(y[valid], risk[valid], sample_weight=w[valid])
    ok = np.isfinite(thresholds) & (thresholds >= 0) & (thresholds <= 1)
    return float(thresholds[ok][np.argmax((tpr-fpr)[ok])])


def train(a, out, model_info):
    manifest(a, out, model_info)
    log("START", "TF prepare", f"Y={a.trait}; biom={a.biom}")
    p, x, c, observed, y, w, masks, prep, clinical, km, features = cohort(a, out)
    build, tune, cal, test = [masks[s] for s in ["build", "tune", "calibration", "test"]]
    selected, importance = select_features(x[build], y[build], w[build], a.max_features)
    names = np.asarray(features)[prep.keep]
    pd.DataFrame({"feature": names, "build_weighted_association": importance,
                  "selected": np.isin(np.arange(len(names)), selected)}).to_csv(out/"feature_selection.csv", index=False)
    z = x[:, selected]
    if a.with_clinical:
        z = np.c_[z, c]
    known = build & (w > 0)
    if min(np.bincount(y[known], minlength=2)) < a.folds:
        raise ValueError("Not enough known build cases/controls for context exclusion folds")
    val = np.flatnonzero(tune & (w > 0))
    if a.validation_samples and len(val) > a.validation_samples:
        from sklearn.model_selection import train_test_split
        val, _ = train_test_split(val, train_size=a.validation_samples, stratify=y[val], random_state=a.seed+53)
    p.iloc[val][[a.id_col]].to_csv(out/"early_stop_ids.csv", index=False)
    ids = p[a.id_col].to_numpy(str)
    groups = p[a.group_col].to_numpy(str) if a.group_col else None
    log("DONE", "TF prepare", f"N={len(p)}; retained={x.shape[1]}; selected={len(selected)}; train_known={known.sum()}; tune_early_stop={len(val)}")
    estimator, checkpoint, context = train_transformer(z[known], y[known], w[known], ids[known],
        None if groups is None else groups[known], z[val], y[val], w[val], a, out)
    holdout = ~build
    raw = {PRIMARY: predict(estimator, z[holdout], a.predict_batch_size)}
    estimator.model_.cpu()
    del estimator
    # Fit comparators to exactly the same selected molecular inputs.
    models, tuning = {}, []
    for name, block in [("clinical", c), ("elasticnet", z)]:
        model, rows = tune_logistic(block, y, w, build, tune, name, a, .5 if name == "elasticnet" else 0)
        models[name] = model
        tuning.extend(rows)
        raw[name] = model.predict_proba(block[holdout])[:, 1]
    pd.DataFrame(tuning).to_csv(out/"comparator_tuning.csv", index=False)
    raw["constant"] = np.full(holdout.sum(), np.average(y[build], weights=w[build]))
    calrel, tunerel, testrel = cal[holdout], tune[holdout], test[holdout]
    calibrators, thresholds, probabilities = {}, {}, {}
    for name, risk in raw.items():
        calibrator = RiskCalibrator().fit(risk[calrel], y[cal], w[cal])
        calibrators[name] = calibrator
        probabilities[name] = calibrator.predict(risk)
        raw_threshold = decision_threshold(y[tune], w[tune], risk[tunerel])
        thresholds[name] = float(calibrator.predict(np.array([raw_threshold]))[0])
    info = pd.DataFrame({a.id_col: ids[test], "missing_fraction": 1-observed[test].mean(1)})
    for name in raw:
        info[name+"_raw"] = raw[name][testrel]
        info[name] = probabilities[name][testrel]
        info[name+"_predicted_Y"] = (probabilities[name][testrel] >= thresholds[name]).astype(int)
    bundle = dict(version=VERSION, config=vars(a), primary=PRIMARY, prep=prep, clinical=clinical,
        censoring=km, features=features, selected=selected, transformer=checkpoint,
        context_x=z[known][context], context_y=y[known][context], context_ids=ids[known][context],
        context_groups=None if groups is None else groups[known][context], models=models,
        calibrators=calibrators, thresholds=thresholds, constant=float(raw["constant"][0]), pretrained=model_info)
    joblib.dump(bundle, out/"model_bundle.joblib", compress=3)
    info.to_csv(out/"test_predictions.csv", index=False)
    p.loc[test, [a.id_col, "time", "event"]+([a.group_col] if a.group_col else [])].to_csv(out/"test_outcomes.csv", index=False)
    dump(out/"thresholds.json", dict(method="Youden J on tune only; mapped through calibration fit", thresholds=thresholds))
    dump(out/"MODEL_FROZEN.json", dict(version=VERSION, model=PRIMARY, test_outcomes_used_for_training=False,
        artifact=fingerprints([out/"model_bundle.joblib"], True)))
    dump(out/"TRAIN_DONE.json", dict(version=VERSION))
    log("DONE", "TF model_frozen")


def evaluate(out):
    if not (out/"MODEL_FROZEN.json").is_file():
        raise ValueError("No frozen TF model")
    bundle = joblib.load(out/"model_bundle.joblib")
    cfg = bundle["config"]
    p = pd.read_csv(out/"test_outcomes.csv", dtype={cfg["id_col"]: str})
    pred = pd.read_csv(out/"test_predictions.csv", dtype={cfg["id_col"]: str})
    if not p[cfg["id_col"]].equals(pred[cfg["id_col"]]):
        raise ValueError("Test prediction/outcome IDs differ")
    y, w = bundle["censoring"].labels_weights(p)
    rows, calibration = [], []
    for name in bundle["calibrators"]:
        for stage, column in [("raw", name+"_raw"), ("calibrated", name)]:
            risk = pred[column].to_numpy()
            rows.append(dict(model=name, stage=stage, n=len(y), known=int((w > 0).sum()), **metrics(y, w, risk)))
        label = pred[name+"_predicted_Y"].to_numpy()
        tn, fp, fn, tp = confusion_matrix(y, label, labels=[0, 1], sample_weight=w).ravel()
        rows[-1].update(threshold=bundle["thresholds"][name], sensitivity=tp/(tp+fn) if tp+fn else None,
                        specificity=tn/(tn+fp) if tn+fp else None)
        bins = pd.qcut(pred[name], 10, duplicates="drop")
        for interval in bins.cat.categories:
            mask = bins.eq(interval).to_numpy()
            if w[mask].sum() > 0:
                calibration.append(dict(model=name, bin=str(interval), n=int(mask.sum()),
                    predicted=float(np.average(pred.loc[mask, name], weights=w[mask])),
                    observed=float(np.average(y[mask], weights=w[mask]))))
    result = pd.DataFrame(rows)
    result.to_csv(out/"test_metrics.csv", index=False)
    pd.DataFrame(calibration).to_csv(out/"test_calibration.csv", index=False)
    pred["observed_Y"] = np.where(w > 0, y, np.nan)
    pred["IPCW"] = w
    pred.to_csv(out/"test_individuals.csv", index=False)
    report = ["# Panome TF: supervised molecular disease prediction", "",
        f"Outcome: {cfg['trait']}; molecular layer: {cfg['biom']}; horizon: {cfg['horizon']} years.",
        "Synthetic data." if cfg["demo"] else "Research predictions on the supplied cohort.", "",
        "| Model | IPCW AUC | IPCW AUPRC | IPCW Brier |", "|---|---:|---:|---:|"]
    for r in result[result.stage.eq("calibrated")].itertuples():
        report.append(f"| {r.model} | {r.AUC_IPCW:.4f} | {r.AUPRC_IPCW:.4f} | {r.Brier_IPCW:.4f} |")
    report += ["", "The model, features, context, calibration and decision threshold were frozen before test scoring.",
        "TabICLv2 uses labeled build context plus a learned classifier head; it does not literally copy donor Y.",
        "Fine-tuning uses exact IPCW query loss; context sampling is only an approximate weighting mechanism.",
        "Epoch zero is eligible for selection if fine-tuning does not improve tune logloss. Inspect training_summary.json.",
        "Feature screening is univariate; excluded assays may carry interactions. Compare --max-features settings on development data only.",
        "Death is censored: these are net risks under the censoring assumptions, not competing-risk cumulative incidences.",
        "A predicted label is a thresholded risk estimate, not a diagnosis. Test rankings are descriptive, not a model-selection rule."]
    (out/"REPORT.md").write_text("\n".join(report)+"\n")
    dump(out/"DONE.json", dict(version=VERSION, n_test=len(y), results_synthetic=cfg["demo"]))
    log("DONE", "TF evaluation", f"N={len(y)}; known={int((w > 0).sum())}")


def project(a, out):
    bundle = joblib.load(out/"model_bundle.joblib")
    cfg = bundle["config"]
    idcol = cfg["id_col"]
    cols = list(dict.fromkeys([idcol]+words(cfg["covariates"])+words(cfg["residualize"])+([cfg["group_col"]] if cfg["group_col"] else [])))
    p = read_table(a.phe_file, idcol, cols, a.r_bin)
    omics = read_table(a.omics_file, idcol, r_bin=a.r_bin)
    mapping = None
    if cfg["biom"] == "met" and a.met_input == "raw":
        omics, mapping = map_metabolites(omics, cfg["met_map"], idcol, nonnegative=cfg["transform"] == "log1p")
    elif cfg["biom"] == "prot":
        omics.columns = [c if c == idcol else str(c).upper() for c in omics]
    if omics.columns.duplicated().any():
        raise ValueError("Assay names collide after normalization")
    if not set(omics[idcol]) <= set(p[idcol]):
        raise ValueError("Query IDs lack baseline metadata")
    if set(omics[idcol]) & set(bundle["context_ids"]):
        raise ValueError("Queries overlap labeled training context; supply independent samples")
    p = p.set_index(idcol).loc[omics[idcol]].reset_index()
    if cfg["group_col"]:
        group = p[cfg["group_col"]]
        if group.isna().any() or set(group.astype(str)) & set(bundle["context_groups"]):
            raise ValueError("Missing query families or families overlapping labeled context")
    absent = sorted(set(bundle["features"])-set(omics))
    omics = omics.reindex(columns=[idcol]+bundle["features"])
    x, observed = bundle["prep"].transform(numeric(omics, bundle["features"], "TF projection"), p)
    z = x[:, bundle["selected"]]
    c = bundle["clinical"].transform(p)
    if cfg["with_clinical"]:
        z = np.c_[z, c]
    a.model_path = cfg["model_path"]
    estimator = restore_predictor(bundle, a)
    raw = {PRIMARY: predict(estimator, z, a.predict_batch_size),
        "clinical": bundle["models"]["clinical"].predict_proba(c)[:, 1],
        "elasticnet": bundle["models"]["elasticnet"].predict_proba(z)[:, 1],
        "constant": np.full(len(p), bundle["constant"])}
    missing = 1-observed.mean(1)
    unknown = unknown_rows(p, [bundle["clinical"], getattr(bundle["prep"], "design", None)])
    result = pd.DataFrame({idcol: p[idcol], "missing_fraction": missing, "unknown_category": unknown,
                           "passes_input_qc": (missing <= cfg["sample_missing"]) & ~unknown})
    for name, risk in raw.items():
        result[name+"_raw"] = risk
        result[name] = bundle["calibrators"][name].predict(risk)
        result[name+"_predicted_Y"] = (result[name] >= bundle["thresholds"][name]).astype(int)
    result["released_risk"] = result[PRIMARY].where(result.passes_input_qc)
    output = Path(a.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    if mapping is not None:
        mapping.to_csv(output.with_name(output.stem+"_metabolite_mapping.csv"), index=False)
    dump(output.with_name(output.stem+"_audit.json"), dict(absent_assays=absent, outcome_columns_required=False,
         n=len(p), passes_input_qc=int(result.passes_input_qc.sum())))
    log("DONE", "TF project", f"N={len(p)}; file={output}")
