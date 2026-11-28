"""Decompose fitted Disco genetic terms into regional slopes and one nonregional PC offset."""
import json
import joblib
import numpy as np
import pandas as pd
from abm2 import RAW_COLUMNS
from core2 import POPS, NumericDesign, aligned, atomic_json, atomic_table, check_roles, read_table
from disco2 import DiscoDivas
from output_layout import project_path, trait_path


def decompose(model, frame):
    """Derivative of the interpolated genetic score with respect to the four raw discovery scores."""
    names = list(model.anchor_names_)
    weights = model.interpolation_weights(frame)
    coefficients = np.stack([model.anchor_models_[a]["score_coef"] for a in names])
    coefficients = coefficients / model.score_scale_[None, :] / model.normalizer_.scale_[:, None]
    slopes = weights @ coefficients
    pcs = frame[list(model.pc_cols_)].to_numpy(float)
    pcfit = np.column_stack((np.ones(len(frame)), pcs)) @ model.normalizer_.coef_
    centering = coefficients @ model.score_mean_
    offset = -np.sum(weights * (centering + (pcfit + model.normalizer_.mean_) / model.normalizer_.scale_), axis=1)
    gates = np.zeros((len(frame), len(POPS)))
    gates[:, [POPS.index(a) for a in names]] = weights
    return gates, slopes, offset


def fit(cfg, root, trait, out):
    frame = read_table(project_path(root, "tmp/cohort/participants.tsv.gz"))
    check_roles(frame)
    baseline = aligned(read_table(trait_path(root, trait, "baseline.predictions.tsv.gz")), frame.id)
    dev = aligned(read_table(project_path(root, f"tmp/cohort/{trait}.development.tsv.gz")), frame.id, allow_missing=True)
    upstream = aligned(read_table(trait_path(root, trait, "predictions.tsv.gz")), frame.id)
    meta = json.loads(trait_path(root, trait, "benchmark.json").read_text())
    columns = meta["geometry_columns"]  # 2.csx records these in canonical population order.
    if len(columns) != len(POPS) or len(set(columns)) != len(columns):
        raise ValueError("Expected one raw score per CSx channel")
    frame[columns] = baseline[RAW_COLUMNS].to_numpy(float)
    y = dev.y.to_numpy(float)
    train = frame.role.eq("train").to_numpy()
    if dev.loc[frame.role.eq("test"), "y"].notna().any():
        raise ValueError("Development input contains TEST labels")
    folds = baseline.oof_fold.to_numpy(int)
    if frame.loc[train].assign(fold=folds[train]).groupby("group_id").fold.nunique().max() != 1:
        raise ValueError("Families cross Disco OOF folds")
    fitted = joblib.load(trait_path(root, trait, "disco/model.joblib"))["disco"]
    if tuple(columns) != tuple(fitted.score_cols_) or set(fitted.fit_ids_) != set(frame.loc[train & np.isfinite(y), "id"].astype(str)):
        raise ValueError("Disco and CSx belong to different inputs or training rosters")
    groups = np.unique(folds[train])
    if len(groups) < 2 or min(groups) < 0:
        raise ValueError("Missing Disco family folds")
    gates, slopes, offset = decompose(fitted, frame)
    cv = cfg["benchmark"]["covariates"]
    for k in groups:
        held = train & (folds == k)
        tr = train & (folds != k) & np.isfinite(y)
        model = DiscoDivas(**fitted.get_params(deep=False))
        model.fit(frame.loc[tr], y[tr], covariates=NumericDesign(cv).fit_transform(frame.loc[tr]))
        gates[held], slopes[held], offset[held] = decompose(model, frame.loc[held])
    reconstructed = np.sum(frame[columns].to_numpy(float) * slopes, axis=1) + offset
    expected = upstream["disco.genetic_score"].to_numpy(float)
    if not np.isfinite(reconstructed).all() or not np.allclose(reconstructed, expected, rtol=3e-4, atol=3e-4):
        raise ValueError("Global scores reconstructed from regional Disco slopes disagree with 3.disco. Check the fitted model, OOF folds and four CSx score columns.")
    data = pd.DataFrame({"id": frame.id, "pc_offset": offset})
    for j, a in enumerate(POPS):
        data["gate_"+a] = gates[:, j]
        data["slope_"+a] = slopes[:, j]
    atomic_table(data, out/"weights.tsv.gz")
    atomic_json(dict(source="3.disco", trait=trait, channels=list(POPS), folds=len(groups),
                     max_abs_reconstruction_error=float(np.max(np.abs(reconstructed-expected))),
                     test_labels_used=False), out/"model.json")
