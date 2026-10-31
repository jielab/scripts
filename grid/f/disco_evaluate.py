#!/usr/bin/env python3
"""Development-only outcome adapter for pinned Disco scores; final never fits."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from abm_core.pipeline import parser, configuration
from abm_core.grid_inputs import load_grid
from abm_core.predictors import fit_linear, PredictionCalibration, require_labels
from abm_core.evaluation import prediction_metrics
from kinship import analysis_guard, compatible, active
from runtime_state import file_identity
from resume import atomic_json, run_lock
from disco_publication import METHOD, read_performance


class ComparisonUnavailable(ValueError):
    pass


def evaluate(root, trait, score_path):
    root, score_path = Path(root).resolve(), Path(score_path).resolve()
    out = root / "disco" / trait
    out.mkdir(parents=True, exist_ok=True)
    with analysis_guard(root), run_lock(out):
        existing, _ = read_performance(root, trait)
        if existing is not None:
            print(f"[Disco] REUSE {trait} verified outcome evaluation")
            return
        try:
            args = parser().parse_args(["--analysis-root", str(root), "--traits", trait, "--device", "cpu"])
            config = configuration(args)
            frame, features, _, metadata = load_grid(args, config, trait)
            if not active(root) or not compatible(metadata, active(root)):
                raise ComparisonUnavailable("same-roster outcome comparison requires verified active kinship")
            from twin import read_table
            scores = read_table(score_path).copy()
            id_column = next((c for c in ("eid", "IID", "id") if c in scores), None)
            score_column = next((c for c in ("disco", "PRS") if c in scores), None)
            if not id_column or not score_column:
                raise ComparisonUnavailable("Disco scores lack participant key or pinned score column")
            scores[id_column] = scores[id_column].astype(str)
            if scores[id_column].duplicated().any() or not set(frame.id) <= set(scores[id_column]):
                raise ComparisonUnavailable("Disco requirements not met for the complete frozen cohort; no ancestry subset omitted")
            score = scores.set_index(id_column).loc[frame.id, score_column].to_numpy(float)
            if not np.isfinite(score).all():
                raise ComparisonUnavailable("Disco has undefined scores in the frozen cohort")
            # Same nuisance inputs as GRID; no outcome-based ancestry exclusions.
            side = [c for c in features if not c.startswith("csx.")]
            x = np.column_stack([score, frame[side].to_numpy(float)])
            if not np.isfinite(x).all():
                raise ComparisonUnavailable("same-roster Disco side inputs are incomplete")
            roles, y = frame.role.to_numpy(str), frame.y.to_numpy(float)
            binary = trait == "t2dm"
            for role in ("build", "tune_model", "prediction_calibration"):
                require_labels(y[roles == role], role, binary)
            scaler = StandardScaler().fit(x[roles == "build"])
            x = scaler.transform(x)
            fit_y = y.copy()
            fit_y[~np.isin(roles, ["build", "tune_model"])] = np.nan
            model, _ = fit_linear(x, fit_y, roles, binary, config.get("predictors", {}))
            raw = model.predict_proba(x)[:, 1] if binary else model.predict(x)
            cal = roles == "prediction_calibration"
            calibrator = PredictionCalibration.fit(raw[cal], y[cal], binary)
            test = roles == "test"
            metrics = prediction_metrics(y[test], calibrator.predict(raw[test]), "binary" if binary else "continuous")
            local = out / "_LOCAL_ONLY"
            local.mkdir(exist_ok=True)
            frame[["id", "group_id", "role"]].to_csv(local / "roles.tsv.gz", sep="\t", index=False)
            provenance = {k: metadata[k] for k in ("grouping_status", "groups_verified", "kinship_version", "kinship_token")}
            table = pd.DataFrame([dict(trait=trait, method=METHOD, source="3.disco.sh", stage="test", status="completed",
                training_role="build", tuning_role="tune_model", calibration_role="prediction_calibration",
                roster_status="exact_frozen_test_all_retained", **provenance, **metrics)])
            table.to_csv(out / "performance.tsv", sep="\t", index=False)
            atomic_json(out / "performance.json", dict(schema="grid.disco.performance.1", status="completed", trait=trait,
                **provenance, role_counts={k:int(v) for k,v in frame.role.value_counts().sort_index().items()},
                roster_identity=file_identity(root / "prsformer/cache/data.tsv.gz"),
                files={str(p):file_identity(p) for p in (score_path, out / "performance.tsv", local / "roles.tsv.gz",
                    root / "cache/abm/prepared.joblib")}))
            valid, reason = read_performance(root, trait)
            if valid is None:
                raise ComparisonUnavailable(reason)
            print(f"[Disco] completed same-roster {trait} outcome evaluation")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            local = out / "_LOCAL_ONLY"
            local.mkdir(exist_ok=True)
            atomic_json(local / "evaluation_failure.json", dict(error=str(exc), error_type=type(exc).__name__))
            reason = str(exc) if isinstance(exc, ComparisonUnavailable) else "Frozen cohort inputs or development fitting requirements failed validation; inspect local Disco diagnostics"
            pd.DataFrame([dict(trait=trait, method=METHOD, source="3.disco.sh", status="unavailable", reason=reason)]).to_csv(out / "performance.tsv", sep="\t", index=False)
            atomic_json(out / "performance.json", dict(schema="grid.disco.performance.1", status="unavailable", reason=reason))
            print(f"[Disco] unavailable {trait}: {reason}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis-root", required=True)
    p.add_argument("--trait", required=True, choices=["height", "ldl", "t2dm"])
    p.add_argument("--scores", required=True)
    a = p.parse_args()
    evaluate(a.analysis_root, a.trait, a.scores)


if __name__ == "__main__":
    main()
