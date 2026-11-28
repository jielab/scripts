"""CALIBRATION-only convex ensemble of original and ancestry-interpolated regional predictors."""
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from core2 import aligned, atomic_json, atomic_table, check_roles, read_table
from output_layout import project_path, trait_path


def simplex_weights(predictions, outcomes):
    # Center on a frozen reference. This keeps height in cm numerically well-conditioned.
    delta = predictions - predictions[:, [0]]
    gram = delta.T @ delta / len(delta)
    cross = delta.T @ (outcomes - predictions[:, 0]) / len(delta)
    scale = max(np.max(np.abs(gram)), np.max(np.abs(cross)), 1e-12)
    gram, cross = gram/scale, cross/scale
    p = predictions.shape[1]
    fit = minimize(lambda w: (.5*w@gram@w-w@cross, gram@w-cross), np.full(p, 1/p),
                   jac=True, method="SLSQP", bounds=[(0., 1.)]*p,
                   constraints=[{"type":"eq", "fun":lambda w: w.sum()-1,
                                 "jac":lambda w: np.ones(p)}], options={"ftol":1e-12,"maxiter":300})
    if not fit.success:
        raise ValueError("Convex stack optimization failed: "+fit.message)
    w = np.clip(fit.x, 0, 1)
    return w/w.sum()


def fit(cfg, root, trait, out):
    cohort = read_table(project_path(root, "tmp/cohort/participants.tsv.gz"))
    check_roles(cohort)
    dev = aligned(read_table(project_path(root, f"tmp/cohort/{trait}.development.tsv.gz")), cohort.id, allow_missing=True)
    if dev.loc[cohort.role.eq("test"),"y"].notna().any():
        raise ValueError("Development input contains TEST labels")
    cal = cohort.role.eq("calibration").to_numpy() & np.isfinite(dev.y.to_numpy(float))
    if not cal.any():
        raise ValueError("Missing CALIBRATION outcomes for stacking")
    components = {"disco":"disco.prediction", "block.glmnet":"block.glmnet.prediction",
                  "block.glmnet.disco":"block.glmnet.disco.unblended",
                  "blocknet.disco":"blocknet.disco.unblended"}
    cols = []
    for name, column in components.items():
        path = trait_path(root, trait, "predictions.tsv.gz") if name=="disco" else root/trait/name/"predictions.tsv.gz"
        frame = read_table(path, ["id", column])
        cols.append(aligned(frame, cohort.id)[column].to_numpy(float))
    matrix = np.column_stack(cols)
    if not np.isfinite(matrix).all() or (trait=="t2dm" and ((matrix<0)|(matrix>1)).any()):
        raise ValueError("Invalid frozen predictions for stacking")
    weights = simplex_weights(matrix[cal], dev.y.to_numpy(float)[cal])
    atomic_table(pd.DataFrame({"id":cohort.id, "block.stack.prediction":matrix@weights}), out/"predictions.tsv.gz")
    atomic_json(dict(trait=trait, weights=dict(zip(components,weights.tolist())), calibration=int(cal.sum()),
                     objective="Brier" if trait=="t2dm" else "MSE", test_labels_used=False), out/"model.json")
