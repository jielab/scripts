"""Pooled Clinical/MoE comparators sharing the ABM development roles."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from .predictors import fit_linear, fit_moe, PredictionCalibration


def fit_comparators(frame, features, csx_columns, config, include_moe=True):
    roles, y = frame.role.to_numpy(str), frame.y.to_numpy(float)
    binary = config["_binary"]
    build, cal = roles == "build", roles == "prediction_calibration"
    results, models, status = {}, {}, []
    clinical = [c for c in features if not c.startswith("csx.")]
    x = frame[features].to_numpy(float)
    scale = StandardScaler().fit(x[build])
    z = scale.transform(x)
    for name in (["Clinical", "csx-moe"] if include_moe else ["Clinical"]):
        if name == "csx-moe" and build.sum() < max(
                200, 2 * int(config.get("moe_min_n", 100)),
                max(config.get("moe_clusters", [4, 8, 16]))):
            status.append(dict(method=name, status="unavailable",
                               reason="insufficient_build_rows_for_prespecified_moe",
                               n_build=int(build.sum())))
            continue
        if name == "Clinical":
            cols = [features.index(c) for c in clinical]
            m, history = fit_linear(z[:, cols], y, roles, binary, config)
            raw = m.predict_proba(z[:, cols])[:, 1] if binary else m.predict(z[:, cols])
        else:
            cols = csx_columns
            m, history = fit_moe(z, y, roles, binary, csx_columns, config)
            raw, _ = m.predict(z[:, cols], z)
        calibrator = PredictionCalibration.fit(raw[cal], y[cal], binary)
        results[name] = calibrator.predict(raw)
        models[name] = dict(model=m, scaler=scale, columns=cols, calibration=calibrator, history=history)
        status.append(dict(method=name, status="completed", fit_role="build",
                           tuning_role="tune_model", calibration_role="prediction_calibration"))
    return results, models, pd.DataFrame(status)
