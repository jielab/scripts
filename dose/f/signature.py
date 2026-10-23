#!/usr/bin/env python3
"""signature analysis steps. Run with dose.sh; imports do not execute analysis."""
import importlib.util
import sys
from pathlib import Path

if "dose_common" not in sys.modules:
    spec = importlib.util.spec_from_file_location("dose_common", Path(__file__).with_name("0.common.py"))
    common = importlib.util.module_from_spec(spec)
    sys.modules["dose_common"] = common
    spec.loader.exec_module(common)
else:
    common = sys.modules["dose_common"]


# ===========================================================================
# signature / protein_lasso
def run_protein_lasso():
    """Build a DOSE-related proteomic signature using Gaussian LASSO."""

    from pathlib import Path
    output_dir = common.output_dir; jobs = common.jobs
    import json
    import os
    import sys

    import gc
    import joblib
    import numpy as np
    import pandas as pd
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import Lasso, LassoCV, lasso_path
    from sklearn.model_selection import KFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    # Match Dapu-Li/Sugar_depression/s20_Lasso.py:
    # 10-fold LassoCV, six candidate penalties, and random_state=42.
    SEED = 42
    N_FOLDS = 10
    ALPHAS = np.logspace(-4, 1, 6)
    MAX_ITER = 10000

    try:
        import matplotlib.pyplot as plt
    except ImportError:  # Plots are useful but must not prevent the analysis.
        plt = None


    def project_root() -> Path:
        return output_dir()


    def normalize_eid(values) -> pd.Series:
        numeric = pd.to_numeric(values, errors="raise")
        if numeric.isna().any() or not np.all(np.equal(numeric, np.floor(numeric))):
            raise ValueError("eid must be nonmissing integer-like values")
        return numeric.astype("int64").astype(str)


    def main() -> None:
        root = project_root()
        out = root / "ukb/06signature/protein"
        out.mkdir(parents=True, exist_ok=True)
        input_file = out / "protein_LASSO_input.csv.gz"
        if not input_file.exists():
            raise FileNotFoundError(
                f"{input_file} is missing. Run ./dose.sh signature --steps prepare, which prepares the input in R first."
            )
        data = pd.read_csv(input_file, low_memory=False)
        if not {"eid", "dose"}.issubset(data.columns):
            raise ValueError("Protein LASSO input lacks eid or dose")
        data["eid"] = normalize_eid(data["eid"])
        if data["eid"].duplicated().any():
            raise ValueError("Duplicated eid in protein LASSO input")
        features = [c for c in data.columns if c not in {"eid", "dose"}]
        n_fdr_features = len(features)
        data = data.loc[np.isfinite(pd.to_numeric(data["dose"], errors="coerce"))].copy()
        x = data[features].apply(pd.to_numeric, errors="coerce")
        x = x.replace([np.inf, -np.inf], np.nan)
        usable = x.notna().sum(axis=0).ge(1000) & x.nunique(dropna=True).gt(1)
        dropped = usable.index[~usable].tolist()
        features = usable.index[usable].tolist()
        if not features:
            raise RuntimeError("No selected protein passed LASSO input QC")
        x = x[features].to_numpy(dtype=float, copy=True)
        y_raw = data["dose"].to_numpy(dtype=float, copy=True)
        eids = data["eid"].to_numpy(copy=True)
        del data
        gc.collect()
        y_mean, y_sd = float(y_raw.mean()), float(y_raw.std(ddof=1))
        if not np.isfinite(y_sd) or y_sd <= 0:
            raise RuntimeError("DOSE has zero/invalid SD")
        y = (y_raw - y_mean) / y_sd

        # Persist the costly full fit and every outer fold. Transformed full-data
        # matrices are released before OOF; only one fold owns work buffers.
        cache = common.lasso_checkpoint_cache(out)
        fit_full = cache.cache(common.fit_lasso_full)
        fit_fold = cache.cache(common.predict_lasso_fold)
        params = dict(alphas=ALPHAS, folds=N_FOLDS, max_iter=MAX_ITER, seed=SEED, log1p=False)
        print("Protein full LASSO: " + ("CACHE HIT" if fit_full.check_call_in_cache(x, y, **params) else "RUN"), flush=True)
        fitted = fit_full(x, y, **params)
        pre, lasso_cv, model = fitted["pre"], fitted["cv"], fitted["model"]
        selected, full = fitted["selected"], fitted["full"]
        path_alphas, path_coefs = fitted["path_alphas"], fitted["path_coefs"]
        alpha = float(lasso_cv.alpha_)
        oof_cv = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
        oof = np.full(y.shape[0], np.nan, dtype=float)
        for fold, (train_idx, test_idx) in enumerate(oof_cv.split(x), start=1):
            hit = fit_fold.check_call_in_cache(x, y, train_idx, test_idx, **params)
            print(f"Protein OOF {fold}/{N_FOLDS}: {'CACHE HIT' if hit else 'RUN'}", flush=True)
            prediction, fold_alpha = fit_fold(x, y, train_idx, test_idx, **params)
            oof[test_idx] = prediction
            del prediction, train_idx, test_idx
            gc.collect()
        if not np.all(np.isfinite(oof)):
            raise RuntimeError("OOF prediction contains missing or non-finite values")
        reference = common.lasso_reference_score(x, selected, pre, model)
        del x
        gc.collect()

        full_z = (full - full.mean()) / full.std(ddof=1)
        oof_z = (oof - oof.mean()) / oof.std(ddof=1)
        reference_z = (reference - reference.mean()) / reference.std(ddof=1)
        scores = pd.DataFrame({
            "eid": eids,
            "prot_score": full,
            "prot_score_z": full_z,
            "prot_score_oof": oof,
            "prot_score_oof_z": oof_z,
            "prot_score_reference": reference,
            "prot_score_reference_z": reference_z,
        })
        scores.to_csv(out / "protein_signature_scores.csv", index=False)

        imputer = pre.named_steps["imputer"]
        scaler = pre.named_steps["scaler"]
        coefficients = pd.DataFrame({
            "protein": features,
            "coefficient_standardized": model.coef_,
            "selected": np.abs(model.coef_) > 1e-12,
            "imputation_median": imputer.statistics_,
            "feature_mean_after_imputation": scaler.mean_,
            "feature_scale": scaler.scale_,
        }).sort_values(["selected", "coefficient_standardized"], ascending=[False, False])
        coefficients.to_csv(out / "protein_LASSO_coefficients.csv", index=False)
        coefficients.loc[coefficients["selected"]].to_csv(out / "protein_LASSO_selected.csv", index=False)
        (out / "protein_LASSO_selected.txt").write_text(
            "\n".join(coefficients.loc[coefficients["selected"], "protein"]) + "\n", encoding="utf-8"
        )
        mse = pd.DataFrame({
            "alpha": lasso_cv.alphas_,
            "mean_cv_mse": lasso_cv.mse_path_.mean(axis=1),
            "se_cv_mse": lasso_cv.mse_path_.std(axis=1, ddof=1) / np.sqrt(N_FOLDS),
        })
        mse.to_csv(out / "protein_LASSO_CV_MSE.csv", index=False)
        pd.DataFrame(path_coefs.T, index=path_alphas, columns=features).rename_axis("alpha").to_csv(out / "protein_LASSO_path.csv")
        joblib.dump({"preprocessor": pre, "model": model, "features": features}, out / "protein_LASSO_model.joblib")
        metrics = {
            "seed": SEED, "folds": N_FOLDS, "n_participants": int(len(eids)),
            "n_fdr_features": int(n_fdr_features),
            "n_qc_features": int(len(features)), "n_dropped_features": int(len(dropped)),
            "n_selected": int(selected.size), "alpha_min": alpha,
            "alpha_candidates": ALPHAS.tolist(),
            "alpha_selection_rule": "minimum mean 10-fold CV MSE",
            "dose_mean": y_mean, "dose_sd": y_sd,
            "full_R2": fitted["full_r2"],
            "oof_correlation_with_dose_z": float(np.corrcoef(oof, y)[0, 1]),
            "primary_score_column": "prot_score_z",
            "sensitivity_score_column": "prot_score_oof_z",
            "reference_style_score_column": "prot_score_reference_z",
            "reference_repository": "Dapu-Li/Sugar_depression/s20_Lasso.py",
        }
        (out / "protein_LASSO_summary.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")

        if plt is not None:
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.errorbar(np.log10(mse["alpha"]), mse["mean_cv_mse"], yerr=mse["se_cv_mse"], fmt="o-")
            ax.axvline(np.log10(alpha), color="red", linestyle="--", label=f"alpha={alpha:.3g}")
            ax.set(xlabel="log10(alpha)", ylabel="10-fold CV MSE", title="Proteomic LASSO CV")
            ax.legend()
            fig.tight_layout(); fig.savefig(out / "protein_LASSO_CV_MSE.png", dpi=300); plt.close(fig)
            fig, ax = plt.subplots(figsize=(7, 5))
            for j in range(path_coefs.shape[0]):
                ax.plot(-np.log10(path_alphas), path_coefs[j], alpha=0.7)
            ax.set(xlabel="-log10(alpha)", ylabel="Coefficient", title="Proteomic LASSO path")
            fig.tight_layout(); fig.savefig(out / "protein_LASSO_path.png", dpi=300); plt.close(fig)
        else:
            print("matplotlib is not installed; CSV results were saved and plots were skipped.")
        print(f"Protein LASSO complete: N={len(eids)}, input={len(features)}, selected={selected.size}, alpha={alpha:.6g}")


    try:
        main()
    except Exception as exc:
        print(f"[s5 protein LASSO] ERROR: {exc}", file=sys.stderr)
        raise


# ===========================================================================
# signature / metabolite_lasso
def run_metabolite_lasso():
    """Build a DOSE-related metabolomic signature using Gaussian LASSO."""

    from pathlib import Path
    output_dir = common.output_dir
    from contextlib import contextmanager
    import gc
    import json
    import os
    import sys
    import threading
    import time

    PROGRESS_INTERVAL_SECONDS = 30


    def log(message):
        # Flush every line so R's system2() receives it without Python buffering.
        print(f"[{time.strftime('%H:%M:%S')}] [s5 metabolite] {message}", flush=True)


    @contextmanager
    def progress_stage(label, interval=PROGRESS_INTERVAL_SECONDS):
        """Report the current stage immediately and every 30 seconds while waiting.

        A RUNNING line reports elapsed time, not a convergence percentage or ETA.
        The small monitor thread does not perform model fitting or use randomness.
        """
        started = time.perf_counter()
        cpu_started = time.process_time()
        stopped = threading.Event()
        log(f"START {label}")

        def report_elapsed():
            while not stopped.wait(interval):
                elapsed = (time.perf_counter() - started) / 60
                cpu_elapsed = (time.process_time() - cpu_started) / 60
                log(f"RUNNING {label} | elapsed={elapsed:.1f} min | CPU={cpu_elapsed:.1f} min")

        monitor = threading.Thread(target=report_elapsed, daemon=True)
        monitor.start()
        succeeded = False
        try:
            yield
            succeeded = True
        finally:
            stopped.set()
            monitor.join()
            status = "DONE" if succeeded else "STOPPED"
            log(f"{status} {label} | elapsed={(time.perf_counter() - started) / 60:.1f} min")


    # Display a message before importing the heavier numerical/plotting packages.
    with progress_stage("Load Python packages"):
        import joblib
        import numpy as np
        import pandas as pd
        from sklearn.impute import SimpleImputer
        from sklearn.base import clone
        from sklearn.linear_model import Lasso, LassoCV, lasso_path
        from sklearn.model_selection import KFold
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import FunctionTransformer, StandardScaler

        try:
            import matplotlib.pyplot as plt
        except ImportError:  # Plots are optional; numerical outputs remain available.
            plt = None

    # Use the LASSO tuning convention in Dapu-Li/Sugar_depression/s20_Lasso.py.
    SEED = 42
    N_FOLDS = 10
    ALPHAS = np.logspace(-4, 1, 6)
    # The reference used 10,000 iterations, but the present 168k-participant
    # metabolite matrix did not converge at that limit. A higher ceiling changes
    # neither the alpha grid nor the minimum-CV-error selection rule.
    MAX_ITER = 100000

    def project_root() -> Path:
        return output_dir()


    def normalize_eid(values) -> pd.Series:
        numeric = pd.to_numeric(values, errors="raise")
        if numeric.isna().any() or not np.all(np.equal(numeric, np.floor(numeric))):
            raise ValueError("eid must be nonmissing integer-like values")
        return numeric.astype("int64").astype(str)


    def checked_log1p(x):
        x = np.asarray(x, dtype=float)
        if np.any(np.isfinite(x) & (x < 0)):
            raise ValueError("Negative metabolite value found; log1p is invalid")
        return np.log1p(x)


    def main() -> None:
        run_started = time.perf_counter()
        log(f"SETTINGS | CV folds={N_FOLDS} | alphas={ALPHAS.tolist()} | max_iter={MAX_ITER} | jobs=1")
        root = project_root()
        out = root / "ukb/06signature/metabolite"
        out.mkdir(parents=True, exist_ok=True)
        input_file = out / "metabolite_LASSO_input.csv.gz"
        if not input_file.exists():
            raise FileNotFoundError(
                f"{input_file} is missing. Run ./dose.sh signature --steps prepare, which prepares the input in R first."
            )
        with progress_stage(f"Read input: {input_file.name}"):
            data = pd.read_csv(input_file, low_memory=False)
        log(f"INPUT | rows={len(data):,} | columns={len(data.columns)}")
        if not {"eid", "dose"}.issubset(data.columns):
            raise ValueError("Metabolite LASSO input lacks eid or dose")
        with progress_stage("Input QC and matrix preparation"):
            data["eid"] = normalize_eid(data["eid"])
            if data["eid"].duplicated().any():
                raise ValueError("Duplicated eid in metabolite LASSO input")
            features = [c for c in data.columns if c not in {"eid", "dose"}]
            n_fdr_features = len(features)
            data = data.loc[np.isfinite(pd.to_numeric(data["dose"], errors="coerce"))].copy()
            xdf = data[features].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
            usable = xdf.notna().sum(axis=0).ge(1000) & xdf.nunique(dropna=True).gt(1) & xdf.min(axis=0, skipna=True).ge(0)
            dropped = usable.index[~usable].tolist()
            features = usable.index[usable].tolist()
            if not features:
                raise RuntimeError("No selected metabolite passed LASSO input QC")
            x = xdf[features].to_numpy(dtype=float, copy=True)
            y_raw = data["dose"].to_numpy(dtype=float, copy=True)
            eids = data["eid"].to_numpy(copy=True)
            del data, xdf
            gc.collect()
            y_mean, y_sd = float(y_raw.mean()), float(y_raw.std(ddof=1))
            if not np.isfinite(y_sd) or y_sd <= 0:
                raise RuntimeError("DOSE has zero/invalid SD")
            y = (y_raw - y_mean) / y_sd
        log(f"QC | N={len(eids):,} | candidates={n_fdr_features} | retained={len(features)} | dropped={len(dropped)}")

        cache = common.lasso_checkpoint_cache(out)
        fit_full = cache.cache(common.fit_lasso_full)
        fit_fold = cache.cache(common.predict_lasso_fold)
        params = dict(alphas=ALPHAS, folds=N_FOLDS, max_iter=MAX_ITER, seed=SEED, log1p=True)
        hit = fit_full.check_call_in_cache(x, y, **params)
        with progress_stage(f"Full-sample LASSO: {'CACHE HIT' if hit else 'RUN'}"):
            fitted = fit_full(x, y, **params)
        pre, lasso_cv, model = fitted["pre"], fitted["cv"], fitted["model"]
        selected, full = fitted["selected"], fitted["full"]
        path_alphas, path_coefs = fitted["path_alphas"], fitted["path_coefs"]
        alpha = float(lasso_cv.alpha_)
        log(f"MODEL | best_alpha={alpha:.6g} | selected={selected.size}")
        oof_cv = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
        oof = np.full(y.shape[0], np.nan, dtype=float)
        for fold, (train_idx, test_idx) in enumerate(oof_cv.split(x), start=1):
            hit = fit_fold.check_call_in_cache(x, y, train_idx, test_idx, **params)
            with progress_stage(f"OOF {fold}/{N_FOLDS}: {'CACHE HIT' if hit else 'RUN'}"):
                prediction, fold_alpha = fit_fold(x, y, train_idx, test_idx, **params)
            oof[test_idx] = prediction
            log(f"OOF COMPLETE {fold}/{N_FOLDS} | alpha={fold_alpha:.6g}")
            del prediction, train_idx, test_idx
            gc.collect()
        if not np.all(np.isfinite(oof)):
            raise RuntimeError("OOF prediction contains missing or non-finite values")
        with progress_stage("Reference-style score calculation"):
            reference = common.lasso_reference_score(x, selected, pre, model, log1p=True)
        del x
        gc.collect()

        full_z = (full - full.mean()) / full.std(ddof=1)
        oof_z = (oof - oof.mean()) / oof.std(ddof=1)
        reference_z = (reference - reference.mean()) / reference.std(ddof=1)
        scores = pd.DataFrame({
            "eid": eids,
            "met_score": full,
            "met_score_z": full_z,
            "met_score_oof": oof,
            "met_score_oof_z": oof_z,
            "met_score_reference": reference,
            "met_score_reference_z": reference_z,
        })
        with progress_stage("Save participant scores"):
            scores.to_csv(out / "metabolite_signature_scores.csv", index=False)

        log("START Save coefficients and CV results")
        imputer = pre.named_steps["imputer"]
        scaler = pre.named_steps["scaler"]
        coefficients = pd.DataFrame({
            "metabolite": features,
            "coefficient_standardized_log1p": model.coef_,
            "selected": np.abs(model.coef_) > 1e-12,
            "imputation_median_raw": imputer.statistics_,
            "feature_mean_log1p": scaler.mean_,
            "feature_scale_log1p": scaler.scale_,
        }).sort_values(["selected", "coefficient_standardized_log1p"], ascending=[False, False])
        coefficients.to_csv(out / "metabolite_LASSO_coefficients.csv", index=False)
        coefficients.loc[coefficients["selected"]].to_csv(out / "metabolite_LASSO_selected.csv", index=False)
        (out / "metabolite_LASSO_selected.txt").write_text(
            "\n".join(coefficients.loc[coefficients["selected"], "metabolite"]) + "\n", encoding="utf-8"
        )
        mse = pd.DataFrame({
            "alpha": lasso_cv.alphas_,
            "mean_cv_mse": lasso_cv.mse_path_.mean(axis=1),
            "se_cv_mse": lasso_cv.mse_path_.std(axis=1, ddof=1) / np.sqrt(N_FOLDS),
        })
        mse.to_csv(out / "metabolite_LASSO_CV_MSE.csv", index=False)
        log("DONE Save coefficients and CV results")
        pd.DataFrame(path_coefs.T, index=path_alphas, columns=features).rename_axis("alpha").to_csv(out / "metabolite_LASSO_path.csv")
        joblib.dump({
            "imputer": imputer, "scaler": scaler, "model": model,
            "features": features, "transformation": "log1p"
        }, out / "metabolite_LASSO_model.joblib")
        metrics = {
            "seed": SEED, "folds": N_FOLDS, "n_participants": int(len(eids)),
            "max_iter": MAX_ITER, "oof_outer_jobs": 1,
            "n_fdr_features": int(n_fdr_features),
            "n_qc_features": int(len(features)), "n_dropped_features": int(len(dropped)),
            "n_selected": int(selected.size), "alpha_min": alpha,
            "alpha_candidates": ALPHAS.tolist(),
            "alpha_selection_rule": "minimum mean 10-fold CV MSE",
            "dose_mean": y_mean, "dose_sd": y_sd,
            "full_R2": fitted["full_r2"],
            "oof_correlation_with_dose_z": float(np.corrcoef(oof, y)[0, 1]),
            "primary_score_column": "met_score_z",
            "sensitivity_score_column": "met_score_oof_z",
            "reference_style_score_column": "met_score_reference_z",
            "reference_repository": "Dapu-Li/Sugar_depression/s20_Lasso.py",
            "transformation": "median imputation, log1p, StandardScaler",
        }
        (out / "metabolite_LASSO_summary.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
        log("SAVED Model and summary")
        if plt is not None:
            log("START Save plots")
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.errorbar(np.log10(mse["alpha"]), mse["mean_cv_mse"], yerr=mse["se_cv_mse"], fmt="o-")
            ax.axvline(np.log10(alpha), color="red", linestyle="--", label=f"alpha={alpha:.3g}")
            ax.set(xlabel="log10(alpha)", ylabel="10-fold CV MSE", title="Metabolomic LASSO CV")
            ax.legend(); fig.tight_layout(); fig.savefig(out / "metabolite_LASSO_CV_MSE.png", dpi=300); plt.close(fig)
            fig, ax = plt.subplots(figsize=(7, 5))
            for j in range(path_coefs.shape[0]):
                ax.plot(-np.log10(path_alphas), path_coefs[j], alpha=0.7)
            ax.set(xlabel="-log10(alpha)", ylabel="Coefficient", title="Metabolomic LASSO path")
            fig.tight_layout(); fig.savefig(out / "metabolite_LASSO_path.png", dpi=300); plt.close(fig)
            log("DONE Save plots")
        else:
            log("matplotlib is not installed; CSV results were saved and plots were skipped.")
        log(f"COMPLETE | N={len(eids):,} | input={len(features)} | selected={selected.size} | alpha={alpha:.6g} | total={(time.perf_counter() - run_started) / 60:.1f} min")


    try:
        main()
    except Exception as exc:
        print(f"[s5 metabolite LASSO] ERROR: {exc}", file=sys.stderr, flush=True)
        raise


STEPS = {
    "protein_lasso": run_protein_lasso,
    "metabolite_lasso": run_metabolite_lasso,
}

if __name__ == "__main__":
    common.run_cli(STEPS)
