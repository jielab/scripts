"""Individual additive losses; individual R2/AUC are deliberately not defined."""
from __future__ import annotations

import numpy as np


def assessment_loss(task, loss_name=None):
    if task not in {"continuous", "binary"}:
        raise ValueError("task must be continuous or binary")
    name = loss_name or ("brier" if task == "binary" else "mse")
    allowed = {"mse"} if task == "continuous" else {"brier", "logloss"}
    if name not in allowed:
        raise ValueError(f"Unsupported assessment loss {name!r} for {task}; use {sorted(allowed)}")
    return name


def row_loss(y, prediction, task="continuous", loss_name=None):
    y = np.asarray(y, dtype=float)
    p = np.asarray(prediction, dtype=float)
    if y.shape != p.shape or y.ndim != 1:
        raise ValueError("Outcome and prediction must be aligned one-dimensional arrays")
    if task not in {"continuous", "binary"}:
        raise ValueError("task must be continuous or binary")
    name = assessment_loss(task, loss_name)
    if task == "continuous" and name != "mse":
        raise ValueError("Continuous primary loss is mse")
    if task == "binary":
        valid_y = y[np.isfinite(y)]
        if not np.isin(valid_y, [0, 1]).all():
            raise ValueError("Binary outcomes must be 0, 1 or missing")
        valid_p = p[np.isfinite(p)]
        if ((valid_p < 0) | (valid_p > 1)).any():
            raise ValueError("Binary predictions must be calibrated probabilities")
        if name not in {"brier", "logloss"}:
            raise ValueError("Binary loss must be brier or logloss")
    if name in {"mse", "brier"}:
        return (y - p) ** 2
    p = np.clip(p, 1e-7, 1 - 1e-7)
    return -(y * np.log(p) + (1 - y) * np.log1p(-p))
