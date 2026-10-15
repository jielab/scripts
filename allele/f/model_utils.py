"""Grouped validation shared by the two exploratory model stages."""

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict


def grouped_predictions(model, features, outcome, groups):
    """Return honest held-out predictions, or explicitly mark CV unavailable."""
    groups = np.asarray(groups)
    outcome = np.asarray(outcome)
    if len(np.unique(outcome)) < 2:
        return np.full(len(outcome), np.nan), None, "outcome has only one class"
    if len(np.unique(groups)) < 2:
        return np.full(len(outcome), np.nan), None, "fewer than two validation groups"
    cv = list(
        GroupKFold(n_splits=min(5, len(np.unique(groups)))).split(
            features, outcome, groups
        )
    )
    if any(len(np.unique(outcome[train])) < 2 for train, _ in cv):
        return np.full(len(outcome), np.nan), None, "a training fold has only one class"
    pred = cross_val_predict(model, features, outcome, cv=cv, method="predict_proba")[
        :, 1
    ]
    return pred, float(roc_auc_score(outcome, pred)), "ok"
