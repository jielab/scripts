#!/usr/bin/env python3
"""ml analysis steps. Run with dose.sh; imports do not execute analysis."""
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
# ml / compare
def run_compare():
    """Compare ML algorithms under identical Region_code 0-9 CV."""

    import numpy as np
    import pandas as pd
    OUT = common.ml_tools().OUT; TARGET = common.ml_tools().TARGET; load_ml = common.ml_tools().load_ml; model_names = common.ml_tools().model_names; make_model = common.ml_tools().make_model; predict_score = common.ml_tools().predict_score; calc_metrics = common.ml_tools().calc_metrics; HAVE_XGB = common.ml_tools().HAVE_XGB

    df, features, folds = load_ml()
    if not HAVE_XGB:
        raise RuntimeError("XGBoost is required for the core comparison.")

    models = model_names(include_svm=True, include_catboost=True)
    print("Models to compare:", ", ".join(models))
    print("Features:", len(features))
    print("Folds:", folds)

    rows, preds = [], []
    fit_fold = common.ml_checkpoint_cache(OUT, "compare").cache(common.fit_predict_ml_fold)
    for model_name in models:
        print(f"\n=== {model_name} ===")
        for fold in folds:
            tr = df["region_code"] != fold
            te = ~tr
            X_train, X_test = df.loc[tr, features], df.loc[te, features]
            y_train = df.loc[tr, TARGET].astype(int)
            y_test = df.loc[te, TARGET].astype(int)

            model = make_model(model_name, y_train)
            score, is_prob = fit_fold(model, X_train, y_train, X_test)
            m = calc_metrics(y_test, score, is_prob)
            m.update({
                "Method": model_name, "Fold": fold,
                "n_train": int(tr.sum()), "n_test": int(te.sum()),
                "case_train": int(y_train.sum()), "case_test": int(y_test.sum()),
                "n_features": len(features),
            })
            rows.append(m)
            preds.append(pd.DataFrame({
                "eid": df.loc[te, "eid"].values,
                "region_code": fold, "y_true": y_test.values,
                "score": score, "is_probability": is_prob,
                "Method": model_name, "FeatureSet": f"All{len(features)}",
            }))
            print(f"fold {fold}: AUC={m['AUC']:.4f} AUPRC={m['AUPRC']:.4f} BalAcc={m['BalancedAccuracy_Youden']:.4f}")

    metrics = pd.DataFrame(rows)
    metric_cols = [
        "AUC", "AUPRC", "Sensitivity_Youden", "Specificity_Youden",
        "Accuracy_Youden", "Youden_Index", "BalancedAccuracy_Youden",
        "LogLoss", "Brier"
    ]
    agg = {}
    for c in metric_cols:
        agg[c + "_mean"] = (c, "mean")
        agg[c + "_sd"] = (c, "std")
    summary = metrics.groupby("Method", as_index=False).agg(**agg)

    event_rate = float(df[TARGET].mean())
    summary["EventRate_overall"] = event_rate
    summary["AUPRC_enrichment"] = summary["AUPRC_mean"] / event_rate
    summary["rank_AUC"] = summary["AUC_mean"].rank(ascending=False, method="min")
    summary["rank_AUPRC"] = summary["AUPRC_mean"].rank(ascending=False, method="min")
    summary["Discrimination_mean_rank"] = summary[["rank_AUC", "rank_AUPRC"]].mean(axis=1)
    summary = summary.sort_values(["Discrimination_mean_rank", "rank_AUC", "rank_AUPRC"]).reset_index(drop=True)
    summary["Suggested_best_discrimination"] = np.arange(len(summary)) == 0

    metrics.to_csv(OUT / "Model_Comparison_ByFold.csv", index=False)
    summary.to_csv(OUT / "Model_Comparison_Summary.csv", index=False)
    pd.concat(preds, ignore_index=True).to_csv(OUT / "Model_Comparison_Predictions.csv", index=False)

    print("\nModel comparison complete")
    print(summary[[
        "Method", "AUC_mean", "AUPRC_mean", "AUPRC_enrichment",
        "Sensitivity_Youden_mean", "Specificity_Youden_mean",
        "Youden_Index_mean", "Discrimination_mean_rank",
        "Suggested_best_discrimination"
    ]].to_string(index=False))


# ===========================================================================
# ml / importance
def run_importance():
    """
    Step 2 — Region-stratified 10-fold LightGBM feature importance.

    Outputs cross-validated:
      - TotalGain
      - Split frequency
      - native LightGBM TreeSHAP mean absolute contribution
      - Ensemble mean of normalized Gain/Split/SHAP

    Primary MYTH-style ranking is TotalGain.
    """

    import numpy as np
    import pandas as pd

    OUT = common.ml_tools().OUT; TARGET = common.ml_tools().TARGET; load_ml = common.ml_tools().load_ml; make_model = common.ml_tools().make_model; normalize_importance = common.ml_tools().normalize_importance

    df, features, folds = load_ml()

    gain_sum = np.zeros(len(features))
    split_sum = np.zeros(len(features))
    shap_sum = np.zeros(len(features))
    fold_rows = []
    fit_fold = common.ml_checkpoint_cache(OUT, "importance").cache(common.fit_importance_ml_fold)

    for fold in folds:
        print(f"fold {fold}")

        tr = df["region_code"] != fold
        te = ~tr

        X_train = df.loc[tr, features]
        X_test = df.loc[te, features]
        y_train = df.loc[tr, TARGET].astype(int)

        model = make_model("LightGBM", y_train)
        gain, split, shap_abs = fit_fold(model, X_train, y_train, X_test)

        gain_n = normalize_importance(gain)
        split_n = normalize_importance(split)
        shap_n = normalize_importance(shap_abs)

        gain_sum += gain_n
        split_sum += split_n
        shap_sum += shap_n

        fold_rows.append(pd.DataFrame({
            "Fold": fold,
            "Food": features,
            "Gain": gain_n,
            "Split": split_n,
            "SHAP": shap_n,
        }))

    imp = pd.DataFrame({
        "Food": features,
        "TotalGain_cv": gain_sum / len(folds),
        "TotalSplit_cv": split_sum / len(folds),
        "ShapValues_cv": shap_sum / len(folds),
    })

    # Compatibility with the MYTH naming, where "TotalCover" in the public script
    # is actually based on split frequency.
    imp["TotalCover_cv"] = imp["TotalSplit_cv"]
    imp["Ensemble"] = imp[
        ["TotalGain_cv", "TotalSplit_cv", "ShapValues_cv"]
    ].mean(axis=1)

    imp["Rank_TotalGain"] = imp["TotalGain_cv"].rank(
        ascending=False, method="min"
    ).astype(int)
    imp["Rank_SHAP"] = imp["ShapValues_cv"].rank(
        ascending=False, method="min"
    ).astype(int)
    imp["Rank_Ensemble"] = imp["Ensemble"].rank(
        ascending=False, method="min"
    ).astype(int)

    imp = imp.sort_values(["Rank_TotalGain", "Rank_SHAP"]).reset_index(drop=True)

    pd.concat(fold_rows, ignore_index=True).to_csv(
        OUT / "Importance_ByFold.csv", index=False
    )
    imp.to_csv(OUT / "Food_Importance.csv", index=False)

    print("\nImportance complete")
    print(
        imp[
            ["Rank_TotalGain", "Food", "TotalGain_cv", "ShapValues_cv", "Ensemble"]
        ].to_string(index=False)
    )


# ===========================================================================
# ml / select
def run_select():
    """
    Importance-ranked sequential prefix selection.

    PRIMARY stopping rule (MODERN-like):
    stop at the first newly added food for which the decrease in mean CV log loss
    is < 0.005; retain the prefix immediately before that food.

    MinLoss and OneSE are sensitivity summaries only.
    """

    import numpy as np
    import pandas as pd
    from sklearn.metrics import log_loss, roc_auc_score, average_precision_score
    OUT = common.ml_tools().OUT; TARGET = common.ml_tools().TARGET; load_ml = common.ml_tools().load_ml; make_model = common.ml_tools().make_model; predict_score = common.ml_tools().predict_score

    DELTA_THRESHOLD = 0.005

    df, all_features, folds = load_ml()
    imp = pd.read_csv(OUT / "Food_Importance.csv").sort_values("TotalGain_cv", ascending=False).reset_index(drop=True)
    ranking = imp["Food"].tolist()
    rows = []
    fit_fold = common.ml_checkpoint_cache(OUT, "select").cache(common.fit_predict_ml_fold)

    for k in range(1, len(ranking) + 1):
        use = ranking[:k]
        added = ranking[k - 1]
        fold_loss, fold_auc, fold_auprc = [], [], []
        oof = np.full(len(df), np.nan)
        print(f"k={k:02d}/{len(ranking)} | add {added}")

        for fold in folds:
            tr = df["region_code"] != fold
            te = ~tr
            X_train, X_test = df.loc[tr, use], df.loc[te, use]
            y_train = df.loc[tr, TARGET].astype(int)
            y_test = df.loc[te, TARGET].astype(int)
            model = make_model("LightGBM", y_train)
            pred, is_prob = fit_fold(model, X_train, y_train, X_test)
            if not is_prob:
                raise RuntimeError("LightGBM must return probabilities.")
            oof[np.where(te.to_numpy())[0]] = pred
            fold_loss.append(log_loss(y_test, pred, labels=[0, 1]))
            fold_auc.append(roc_auc_score(y_test, pred))
            fold_auprc.append(average_precision_score(y_test, pred))

        fold_loss = np.asarray(fold_loss)
        fold_auc = np.asarray(fold_auc)
        fold_auprc = np.asarray(fold_auprc)
        rows.append({
            "k": k, "Food_added": added, "Foods": ";".join(use),
            "loss_cv_mean": fold_loss.mean(),
            "loss_cv_sd": fold_loss.std(ddof=1),
            "loss_cv_se": fold_loss.std(ddof=1) / np.sqrt(len(folds)),
            "loss_all_oof": log_loss(df[TARGET].astype(int), oof, labels=[0, 1]),
            "auc_cv_mean": fold_auc.mean(), "auc_cv_sd": fold_auc.std(ddof=1),
            "auc_all_oof": roc_auc_score(df[TARGET].astype(int), oof),
            "auprc_cv_mean": fold_auprc.mean(), "auprc_cv_sd": fold_auprc.std(ddof=1),
            "auprc_all_oof": average_precision_score(df[TARGET].astype(int), oof),
        })

    res = pd.DataFrame(rows)
    res["Delta_loss"] = res["loss_cv_mean"].shift(1) - res["loss_cv_mean"]

    fail = res.loc[(res["k"] > 1) & (res["Delta_loss"] < DELTA_THRESHOLD)]
    if len(fail):
        stop_at_k = int(fail.iloc[0]["k"])
        k_delta = stop_at_k - 1
        stop_food = str(fail.iloc[0]["Food_added"])
    else:
        stop_at_k = np.nan
        k_delta = int(res["k"].max())
        stop_food = ""

    i_min = int(res["loss_cv_mean"].idxmin())
    k_min = int(res.loc[i_min, "k"])
    loss_min = float(res.loc[i_min, "loss_cv_mean"])
    se_min = float(res.loc[i_min, "loss_cv_se"])
    one_se_threshold = loss_min + se_min
    k_1se = int(res.loc[res["loss_cv_mean"] <= one_se_threshold, "k"].min())

    summary = pd.DataFrame([
        {"RankingMethod": "TotalGain", "criterion": "DeltaLogLoss005", "role": "PRIMARY",
         "k": k_delta, "stop_at_k": stop_at_k, "stop_food": stop_food, "threshold": DELTA_THRESHOLD},
        {"RankingMethod": "TotalGain", "criterion": "MinLoss", "role": "SENSITIVITY",
         "k": k_min, "stop_at_k": np.nan, "stop_food": "", "threshold": loss_min},
        {"RankingMethod": "TotalGain", "criterion": "OneSE", "role": "SENSITIVITY",
         "k": k_1se, "stop_at_k": np.nan, "stop_food": "", "threshold": one_se_threshold},
    ])

    res.to_csv(OUT / "SequentialSelection_TotalGain.csv", index=False)
    summary.to_csv(OUT / "Selection_Summary.csv", index=False)
    imp.head(k_delta).to_csv(OUT / "SelectedFeatures_Delta005.csv", index=False)
    imp.head(k_delta).to_csv(OUT / "Primary_Selected_Features.csv", index=False)
    imp.head(k_min).to_csv(OUT / "SelectedFeatures_MinLoss.csv", index=False)
    imp.head(k_1se).to_csv(OUT / "SelectedFeatures_1SE.csv", index=False)

    print("\nSequential selection complete")
    print(summary.to_string(index=False))
    print(f"\nPrimary Delta<0.005 selected k = {k_delta}")
    if np.isfinite(stop_at_k):
        print(f"First sub-threshold decrease occurred when adding k={stop_at_k}: {stop_food}")
    print("\nPrimary selected foods:")
    print("\n".join(imp.head(k_delta)["Food"].tolist()))


# ===========================================================================
# ml / validate
def run_validate():
    """
    MODERN-like multi-method selection validation.

    A) model-based importance: XGBoost gain, Random-Forest impurity
    B) geographic-CV mean absolute SHAP: LGB/XGB/RF
    C) RFE top-k: LGB/XGB/RF/Logistic/SVM
    D) geographic 10-fold predictive comparison across selected feature sets

    Default CORE feature sets:
      AllP, Primary_Delta005, SHAP_Consensus{k}, RFE_Consensus{k}
    Set FEATURE_SET_MODE='extended' for method-specific top-k sets.

    OPTIMIZED VERSION:
    - Importance, SHAP and All-P tree predictions share one fit per fold.
    - Completed fold results and RFE results are cached under OUT/selection_cache_v2.
      Data, estimator parameters, source code and dependency versions invalidate
      caches. Use only your own trusted cache directory (joblib uses pickle).
    - Preserves model constructors, geographic folds, RFE step=1 and 5000 SHAP
      subjects/fold. Reuse assumes deterministic model constructors/seeds.
    - DOSE_SHAP_MAX=500 is an OPTIONAL faster exploratory setting; it can change
      SHAP ranks/selected foods and therefore downstream results. Default: 5000.
    - By default does not reuse unverified s2_1 CSVs. All-P non-tree models may therefore
      require new fits on the first run. Output filenames are unchanged/overwritten.
    - DOSE_REUSE_COMPARE=1 explicitly asserts that existing s2_1 results use the
      same data/features/configuration; subject/fold/outcome checks still run.
    - Exact RF SHAP is split across tree groups in up to 4 processes (DOSE_SHAP_JOBS).
      Signed contributions are combined BEFORE absolute values. More processes
      consume more memory. Set DOSE_SHAP_JOBS=1 to disable multiprocessing.
    - DOSE_RF_SHAP_MAX overrides only the RF explanation sample cap (default 5000).
    - DOSE_THREADS optionally overrides model training threads, not tree counts.
    - SelectionValidation_Timing.csv records wall time/cache hits for each task;
      SelectionValidation_RunSettings.json records settings for this run.
    - Feature selection is still performed before predictive comparison: this is
      not nested CV and not a new independent validation of the selection process.
    """

    import os
    import time
    import hashlib
    import importlib.metadata
    import json
    from pathlib import Path
    ml = common.ml_tools()
    from joblib import Memory, Parallel, delayed, parallel_config, cpu_count
    import numpy as np
    import pandas as pd
    from sklearn.feature_selection import RFE

    OUT = common.ml_tools().OUT; TARGET = common.ml_tools().TARGET; SEED = common.ml_tools().SEED; load_ml = common.ml_tools().load_ml; make_model = common.ml_tools().make_model; make_rfe_model = common.ml_tools().make_rfe_model; model_names = common.ml_tools().model_names; predict_score = common.ml_tools().predict_score; calc_metrics = common.ml_tools().calc_metrics; normalize_importance = common.ml_tools().normalize_importance; HAVE_XGB = common.ml_tools().HAVE_XGB

    if not HAVE_XGB:
        raise RuntimeError("XGBoost is required for selection validation.")

    import xgboost as xgb
    try:
        import shap
    except Exception as e:
        raise RuntimeError("Install shap first: py -m pip install shap") from e

    SHAP_MAX_PER_FOLD = int(os.getenv("DOSE_SHAP_MAX", "5000"))
    RF_SHAP_MAX_PER_FOLD = int(os.getenv('DOSE_RF_SHAP_MAX', str(SHAP_MAX_PER_FOLD)))
    RF_SHAP_JOBS = int(os.getenv('DOSE_SHAP_JOBS', str(min(common.jobs(), cpu_count()))))
    if RF_SHAP_JOBS < 1 or RF_SHAP_MAX_PER_FOLD < 1:
        raise ValueError('DOSE_SHAP_JOBS and DOSE_RF_SHAP_MAX must be positive')
    FEATURE_SET_MODE = "core"   # "core" or "extended"

    df, features, folds = load_ml()
    fingerprint = hashlib.sha256(Path(__file__).read_bytes() + Path(ml.__file__).read_bytes())
    for package in ['numpy', 'pandas', 'scikit-learn', 'joblib', 'lightgbm', 'xgboost', 'catboost', 'shap']:
        try:
            fingerprint.update((package + importlib.metadata.version(package)).encode())
        except importlib.metadata.PackageNotFoundError:
            pass
    cache = Memory(OUT / 'selection_cache_v2' / fingerprint.hexdigest()[:16], verbose=0)
    if common.replace_results():
        cache.clear(warn=False)
    if df['eid'].duplicated().any():
        raise ValueError('Duplicate eid in ML cohort')
    if df['region_code'].isna().any() or not df[TARGET].isin([0, 1]).all():
        raise ValueError('Missing fold or non-binary outcome in ML input')
    for fold in folds:
        if df.loc[df['region_code'] == fold, TARGET].nunique() != 2:
            raise ValueError(f'Fold {fold} needs both cases and controls')
    if SHAP_MAX_PER_FOLD < 1:
        raise ValueError('DOSE_SHAP_MAX must be positive')
    y_all = df[TARGET].astype(int)
    primary = pd.read_csv(OUT / "Primary_Selected_Features.csv")
    primary_set = primary["Food"].tolist()
    if not primary_set or len(set(primary_set)) != len(primary_set) or not set(primary_set).issubset(features):
        raise ValueError('Primary_Selected_Features.csv has invalid/duplicate foods')
    k_primary = len(primary_set)
    # Match every alternative selection method to the primary LightGBM size.
    TOP_N = k_primary
    print("Primary selected k =", k_primary)
    print("Primary foods =", ", ".join(primary_set))
    print(f'Training threads={ml.THREADS}; exact RF SHAP processes={RF_SHAP_JOBS}', flush=True)
    settings = dict(SHAP_MAX_PER_FOLD=SHAP_MAX_PER_FOLD, RF_SHAP_MAX_PER_FOLD=RF_SHAP_MAX_PER_FOLD,
                    RF_SHAP_JOBS=RF_SHAP_JOBS, TrainingThreads=ml.THREADS,
                    ReuseComparison=os.getenv('DOSE_REUSE_COMPARE', '0'),
                    CohortN=len(df), Events=int(y_all.sum()), Features=features,
                    Primary=primary_set, SelectionTopN=TOP_N, SourceFingerprint=fingerprint.hexdigest())
    (OUT / 'SelectionValidation_RunSettings.json').write_text(json.dumps(settings, indent=2), encoding='utf-8')
    timings = []


    def timed_cached(function, stage, method, fold, *args):
        hit = function.check_call_in_cache(*args)
        print(f'[{stage}] {method} fold={fold}: {"CACHE HIT" if hit else "RUN"}', flush=True)
        started = time.perf_counter()
        result = function(*args)
        elapsed = time.perf_counter() - started
        timings.append(dict(Stage=stage, Method=method, Fold=str(fold), CacheHit=hit, Seconds=elapsed))
        pd.DataFrame(timings).to_csv(OUT / 'SelectionValidation_Timing.csv', index=False)
        print(f'[{stage}] completed in {elapsed:.1f}s', flush=True)
        return result

    # ------------------------------------------------------------------
    # A. XGB gain and RF impurity importance, geographic CV
    # ------------------------------------------------------------------
    importance_rows = []  # collected from the SAME fitted models as SHAP below
    all_predictions = {}

    # ------------------------------------------------------------------
    # B. Cross-validated SHAP: LGB/XGB/RF
    # ------------------------------------------------------------------
    rng = np.random.default_rng(SEED)
    shap_rows = []

    def extract_rf_shap(model, X):
        explainer = shap.TreeExplainer(model)
        sv = explainer.shap_values(X, check_additivity=False)
        if isinstance(sv, list):
            arr = np.asarray(sv[-1])
        else:
            arr = np.asarray(sv)
            if arr.ndim == 3:
                arr = arr[:, :, -1]
        if arr.ndim != 2:
            raise RuntimeError(f"Unexpected RandomForest SHAP shape: {arr.shape}")
        return arr


    def rf_shard_shap(model, X):
        # Each worker owns only its tree subset. Native thread pools stay at one.
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=1):
            return extract_rf_shap(model, X)


    def parallel_rf_shap(model, X, jobs):
        from copy import copy
        if jobs == 1:
            return extract_rf_shap(model, X)
        groups = np.array_split(np.arange(len(model.estimators_)), min(jobs, len(model.estimators_)))
        shards = []
        for group in groups:
            shard = copy(model)
            shard.estimators_ = [model.estimators_[i] for i in group]
            shard.n_estimators = len(group)
            shard.n_jobs = 1
            shards.append(shard)
        with parallel_config(backend='loky', inner_max_num_threads=1):
            parts = Parallel(n_jobs=len(shards), verbose=5)(delayed(rf_shard_shap)(shard, X) for shard in shards)
        # TreeSHAP is linear in the forest's tree-average. Combine SIGNED values
        # first; mean(abs(SHAP)) of individual trees is NOT forest SHAP.
        arr = np.zeros_like(parts[0], dtype=float)
        for group, part in zip(groups, parts):
            arr += part * (len(group) / len(model.estimators_))
        return arr


    # Cache only small results, not potentially gigabyte-sized fitted forests.
    # joblib keys include input data, estimator parameters, and function source.
    # The namespace also fingerprints this script, utils and dependency versions.
    @cache.cache
    def fit_explain_predict(model_name, model, X_train, y_train, X_test, X_eval, shap_jobs):
        started = time.perf_counter()
        model.fit(X_train, y_train)
        print(f"  fit completed in {time.perf_counter()-started:.1f}s; SHAP starts", flush=True)
        importance = None
        if model_name == 'LightGBM':
            arr = model.booster_.predict(X_eval, pred_contrib=True)[:, :-1]
        elif model_name == 'XGBoost':
            raw = model.get_booster().get_score(importance_type='gain')
            importance = normalize_importance(np.array([raw.get(f, 0.0) for f in X_train.columns]))
            arr = model.get_booster().predict(xgb.DMatrix(X_eval, feature_names=list(X_train.columns)), pred_contribs=True)[:, :-1]
        else:
            importance = normalize_importance(np.asarray(model.feature_importances_))
            print(f'  RF: {len(model.estimators_)} trees, max fitted depth={max(t.tree_.max_depth for t in model.estimators_)}, SHAP rows={len(X_eval)}, processes={shap_jobs}', flush=True)
            arr = parallel_rf_shap(model, X_eval, shap_jobs)
        vals = normalize_importance(np.abs(arr).mean(axis=0))
        return importance, vals, predict_score(model, X_test)


    @cache.cache
    def fit_predict_cached(model, X_train, y_train, X_test):
        model.fit(X_train, y_train)
        return predict_score(model, X_test)


    @cache.cache
    def rfe_cached(est, X, y, top_n):
        selector = RFE(estimator=est, n_features_to_select=top_n, step=1)
        selector.fit(X, y)
        return selector.ranking_, selector.support_

    for model_name in ["LightGBM", "XGBoost", "RandomForest"]:
        print(f"\n=== CV SHAP: {model_name} ===")
        for fold in folds:
            tr = df["region_code"] != fold
            te = ~tr
            X_train = df.loc[tr, features]
            y_train = df.loc[tr, TARGET].astype(int)
            X_test = df.loc[te, features]
            model = make_model(model_name, y_train)

            shap_limit = RF_SHAP_MAX_PER_FOLD if model_name == 'RandomForest' else SHAP_MAX_PER_FOLD
            if len(X_test) > shap_limit:
                idx = rng.choice(len(X_test), shap_limit, replace=False)
                X_eval = X_test.iloc[np.sort(idx)]
            else:
                X_eval = X_test

            print(f"  fold {fold}: fit/SHAP (or cached result)", flush=True)
            imp, vals, prediction = timed_cached(fit_explain_predict, 'Fit+SHAP', model_name, fold, model_name, model, X_train, y_train, X_test, X_eval, RF_SHAP_JOBS if model_name == 'RandomForest' else 1)
            all_predictions[(model_name, fold)] = prediction
            if imp is not None:
                importance_rows.append(pd.DataFrame({'Method': model_name, 'Fold': fold, 'Food': features, 'Importance': imp}))
            shap_rows.append(pd.DataFrame({
                "Method": model_name, "Fold": fold, "Food": features,
                "MeanAbsSHAP": vals, "n_shap": len(X_eval),
            }))
            print(f"fold {fold}: n_SHAP={len(X_eval)}")

    importance_fold = pd.concat(importance_rows, ignore_index=True)
    importance_fold.to_csv(OUT / 'CrossModel_Importance_ByFold.csv', index=False)
    for method, z in importance_fold.groupby('Method'):
        out = z.groupby('Food', as_index=False)['Importance'].mean()
        out['Rank'] = out['Importance'].rank(ascending=False, method='min').astype(int)
        out.sort_values('Rank').to_csv(OUT / f'{method}_Importance.csv', index=False)

    shap_fold = pd.concat(shap_rows, ignore_index=True)
    shap_fold.to_csv(OUT / "SHAP_ByFold.csv", index=False)
    shap_summary = shap_fold.groupby(["Method", "Food"], as_index=False).agg(
        SHAP_mean=("MeanAbsSHAP", "mean"), SHAP_sd=("MeanAbsSHAP", "std")
    )
    shap_summary["Rank"] = shap_summary.groupby("Method")["SHAP_mean"].rank(ascending=False, method="min").astype(int)
    shap_summary = shap_summary.sort_values(["Method", "Rank"])
    shap_summary.to_csv(OUT / "SHAP_Summary.csv", index=False)

    shap_wide = shap_summary.pivot(index="Food", columns="Method", values="Rank").reset_index()
    shap_rank_cols = [c for c in shap_wide.columns if c != "Food"]
    shap_wide["MeanRank"] = shap_wide[shap_rank_cols].mean(axis=1)
    shap_wide["MedianRank"] = shap_wide[shap_rank_cols].median(axis=1)
    shap_wide = shap_wide.sort_values(["MeanRank", "MedianRank", "Food"])
    shap_wide.to_csv(OUT / "SHAP_CrossModel_Ranks.csv", index=False)

    # ------------------------------------------------------------------
    # C. RFE top-k robustness
    # ------------------------------------------------------------------
    rfe_models = ["LightGBM", "XGBoost", "RandomForest", "Logistic", "SVM"]
    rfe_rows = []
    X_all = df[features]
    for model_name in rfe_models:
        print(f"\n=== RFE top-{TOP_N}: {model_name} ===")
        est = make_rfe_model(model_name, y_all)
        ranking, support = timed_cached(rfe_cached, 'RFE', model_name, 'all', est, X_all, y_all, min(TOP_N, len(features)))
        rfe_rows.append(pd.DataFrame({
            "Food": features, "Method": model_name,
            "RFE_rank": ranking.astype(int),
            "Selected": support.astype(bool),
        }))

    rfe = pd.concat(rfe_rows, ignore_index=True)
    rfe.to_csv(OUT / "RFE_Ranks.csv", index=False)
    rfe_wide = rfe.pivot(index="Food", columns="Method", values="RFE_rank").reset_index()
    rfe_rank_cols = [c for c in rfe_wide.columns if c != "Food"]
    rfe_wide["MeanRank"] = rfe_wide[rfe_rank_cols].mean(axis=1)
    rfe_wide["MedianRank"] = rfe_wide[rfe_rank_cols].median(axis=1)
    selected_freq = rfe.groupby("Food")["Selected"].sum().rename("n_RFE_selected")
    rfe_wide = rfe_wide.merge(selected_freq, on="Food", how="left")
    rfe_wide = rfe_wide.sort_values(["n_RFE_selected", "MeanRank", "MedianRank", "Food"], ascending=[False, True, True, True])
    rfe_wide.to_csv(OUT / "RFE_Consensus.csv", index=False)

    # ------------------------------------------------------------------
    # D. Feature-set definitions
    # ------------------------------------------------------------------
    def method_shap_top(method, n=TOP_N):
        z = shap_summary.loc[shap_summary["Method"] == method].sort_values("Rank")
        return z["Food"].tolist()[:min(n, len(z))]

    def top_from_csv(path, n=TOP_N):
        z = pd.read_csv(path).sort_values("Rank")
        return z["Food"].tolist()[:min(n, len(z))]

    feature_sets = {
        f"All{len(features)}": features,
        f"PrimaryDelta005_{k_primary}": primary_set,
        f"SHAP_Consensus{TOP_N}": shap_wide["Food"].tolist()[:TOP_N],
        f"RFE_Consensus{TOP_N}": rfe_wide["Food"].tolist()[:TOP_N],
    }

    if FEATURE_SET_MODE == "extended":
        feature_sets.update({
            f"LGB_SHAP{TOP_N}": method_shap_top("LightGBM"),
            f"XGB_SHAP{TOP_N}": method_shap_top("XGBoost"),
            f"RF_SHAP{TOP_N}": method_shap_top("RandomForest"),
            f"XGB_Gain{TOP_N}": top_from_csv(OUT / "XGBoost_Importance.csv"),
            f"RF_Impurity{TOP_N}": top_from_csv(OUT / "RandomForest_Importance.csv"),
        })
        for m in rfe_models:
            z = rfe.loc[rfe["Method"] == m].sort_values(["RFE_rank", "Food"])
            feature_sets[f"RFE_{m}{TOP_N}"] = z.loc[z["Selected"], "Food"].tolist()[:TOP_N]

    set_rows = []
    for set_name, fs in feature_sets.items():
        for rank, food in enumerate(fs, start=1):
            set_rows.append({"FeatureSet": set_name, "n_features": len(fs), "rank_in_set": rank, "Food": food})
    pd.DataFrame(set_rows).to_csv(OUT / "Feature_Sets_For_Validation.csv", index=False)

    # ------------------------------------------------------------------
    # E. Geographic 10-fold predictive comparison
    # ------------------------------------------------------------------
    pred_models = model_names(include_svm=True, include_catboost=True)
    print("\nPrediction models:", ", ".join(pred_models))
    print("Feature sets:", ", ".join(feature_sets.keys()))
    metric_rows, pred_rows = [], []

    # All-P tree predictions reuse B. Other models use input-keyed checkpoints.
    # Unverified CSVs from a previous cohort/configuration are NOT reused.
    all_name = f"All{len(features)}"

    # Optional reuse of the user's already completed s2_1 run. CSVs do not contain
    # feature values/model provenance, so reuse requires explicit user opt-in.
    # Verify the observable cohort/fold/label contract for EACH non-tree model.
    if os.getenv('DOSE_REUSE_COMPARE', '0') == '1':
        old_file = OUT / 'Model_Comparison_Predictions.csv'
        if old_file.exists():
            old = pd.read_csv(old_file, dtype={'eid': str})
            expected = df[['eid', 'region_code', TARGET]].copy()
            expected['eid'] = expected['eid'].astype(str)
            required = {'eid', 'region_code', 'y_true', 'score', 'is_probability', 'Method', 'FeatureSet'}
            if required.issubset(old.columns):
                for method in pred_models:
                    if method in ['LightGBM', 'XGBoost', 'RandomForest']:
                        continue
                    z = old.loc[(old['Method'] == method) & (old['FeatureSet'] == all_name)]
                    if len(z) != len(expected) or z['eid'].duplicated().any():
                        print(f'Cannot reuse {method}: missing/duplicate cohort rows', flush=True)
                        continue
                    z = expected.merge(z, on='eid', how='left', suffixes=('_current', '_old'), validate='one_to_one')
                    valid = (z['region_code_current'].eq(z['region_code_old']).all()
                             and z[TARGET].eq(z['y_true']).all() and np.isfinite(z['score']).all()
                             and z['is_probability'].isin([True, False]).all())
                    if not valid:
                        print(f'Cannot reuse {method}: cohort/fold/label/score mismatch', flush=True)
                        continue
                    print(f'Reusing {method} All-P from s2_1 by explicit opt-in; caller must ensure identical features/configuration', flush=True)
                    for fold in folds:
                        part = z.loc[z['region_code_current'] == fold]
                        if part['is_probability'].nunique() != 1:
                            raise ValueError(f'Mixed probability flags for {method}/{fold}')
                        all_predictions[(method, fold)] = (part['score'].to_numpy(), bool(part['is_probability'].iloc[0]))
            else:
                print('Comparison CSV missing required columns; computing missing models instead', flush=True)
        else:
            print('Comparison CSV absent; computing missing models instead', flush=True)

    for set_name, fs in feature_sets.items():
        for model_name in pred_models:
            print(f"\n=== prediction: {set_name} | {model_name} ===")
            for fold in folds:
                tr = df["region_code"] != fold
                te = ~tr
                X_train, X_test = df.loc[tr, fs], df.loc[te, fs]
                y_train = df.loc[tr, TARGET].astype(int)
                y_test = df.loc[te, TARGET].astype(int)
                model = make_model(model_name, y_train)
                if set_name == all_name and (model_name, fold) in all_predictions:
                    score, is_prob = all_predictions[(model_name, fold)]
                else:
                    score, is_prob = timed_cached(fit_predict_cached, f'Predict:{set_name}', model_name, fold, model, X_train, y_train, X_test)
                m = calc_metrics(y_test, score, is_prob)
                m.update({
                    "Method": model_name, "FeatureSet": set_name, "Fold": fold,
                    "n_features": len(fs), "n_test": int(te.sum()), "case_test": int(y_test.sum()),
                })
                metric_rows.append(m)
                pred_rows.append(pd.DataFrame({
                    "eid": df.loc[te, "eid"].values, "region_code": fold,
                    "y_true": y_test.values, "score": score,
                    "is_probability": is_prob, "Method": model_name, "FeatureSet": set_name,
                }))
                print(f"fold {fold}: AUC={m['AUC']:.4f} AUPRC={m['AUPRC']:.4f}")

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "SelectionValidation_ByFold.csv", index=False)
    pd.concat(pred_rows, ignore_index=True).to_csv(OUT / "SelectionValidation_Predictions.csv", index=False)

    metric_cols = ["AUC", "AUPRC", "Sensitivity_Youden", "Specificity_Youden", "Accuracy_Youden", "Youden_Index", "BalancedAccuracy_Youden"]
    agg = {}
    for c in metric_cols:
        agg[c + "_mean"] = (c, "mean")
        agg[c + "_sd"] = (c, "std")
    summary = metrics.groupby(["Method", "FeatureSet", "n_features"], as_index=False).agg(**agg)
    summary.to_csv(OUT / "SelectionValidation_Summary.csv", index=False)

    # Compact selection-method summary.
    base_imp = pd.read_csv(OUT / "Food_Importance.csv")[["Food", "Rank_TotalGain", "Rank_SHAP"]]
    xgb_imp = pd.read_csv(OUT / "XGBoost_Importance.csv")[["Food", "Rank"]].rename(columns={"Rank": "Rank_XGB_Gain"})
    rf_imp = pd.read_csv(OUT / "RandomForest_Importance.csv")[["Food", "Rank"]].rename(columns={"Rank": "Rank_RF_Impurity"})
    shap_cons = shap_wide[["Food", "MeanRank"]].rename(columns={"MeanRank": "MeanRank_SHAP"})
    rfe_cons = rfe_wide[["Food", "MeanRank", "n_RFE_selected"]].rename(columns={"MeanRank": "MeanRank_RFE"})
    tab = base_imp.merge(xgb_imp, on="Food", how="left").merge(rf_imp, on="Food", how="left").merge(shap_cons, on="Food", how="left").merge(rfe_cons, on="Food", how="left")
    tab["Primary_Delta005"] = tab["Food"].isin(primary_set)
    tab = tab.sort_values(["Primary_Delta005", "Rank_TotalGain"], ascending=[False, True])
    tab.to_csv(OUT / "Selection_Methods_Summary.csv", index=False)

    print("\nSelection validation complete")
    print("\nPrimary selected set:")
    print(", ".join(primary_set))
    print("\nPrediction summary:")
    print(summary.sort_values(["FeatureSet", "AUC_mean"], ascending=[True, False]).to_string(index=False))


# ===========================================================================
# ml / bootstrap
def run_bootstrap():
    """
    ml / bootstrap

    Post-hoc evaluation of the already generated OOF predictions.

    Adds:
    1) Bootstrap 95% CIs for Primary-set performance:
       AUC, AUPRC, Accuracy, Sensitivity/Recall, Specificity,
       Precision, F1, Youden index and balanced accuracy.
    2) Paired DeLong tests for AUC.
    3) Paired bootstrap comparisons for Delta AUC and Delta AUPRC.
    4) Confusion matrices for the four core models shown in the MYTH-like figure.

    IMPORTANT
    ---------
    - No ML model is retrained here.
    - Input is SelectionValidation_Predictions.csv.
    - Bootstrap resampling is paired whenever two predictions are compared.
    - Threshold-dependent metrics use a fixed OOF Youden cutoff estimated once
      from the complete OOF predictions for that model/feature set.
      They are descriptive; AUC/AUPRC remain the main discrimination metrics.

    Environment variable:
        DOSE_BOOT_N
    can override the default 1000 bootstrap replicates.

    OPTIMIZED VERSION:
    - Exact case-containing score-block AUC/AP calculation, including ties.
    - Vectorized midranks and cached threshold categories.
    - Same bootstrap sampling, replicate count, seeds, fixed thresholds and output
      schemas as the original; only floating-point rounding may differ.
    - DOSE_ML_OUT optionally overrides the input/output folder.
    - Completed bootstrap groups are cached under OUT/bootstrap_cache_v2; input
      arrays, seeds, replicate counts, code and dependencies determine cache keys.
      Use only a trusted local cache directory (joblib uses pickle).
    - Inference remains conditional on existing OOF scores: it does not include
      feature-selection/retraining uncertainty. Original bootstrap sign-tail
      summaries are retained for compatibility; they are not a null-resampled test.
    """

    from pathlib import Path
    output_dir = common.output_dir
    import os
    import math
    import time
    import hashlib
    import importlib.metadata
    from joblib import Memory
    import numpy as np
    import pandas as pd
    from scipy.stats import norm, rankdata
    from sklearn.metrics import (
        roc_auc_score, average_precision_score, roc_curve, confusion_matrix
    )

    SEED = 12345
    N_BOOT = int(os.getenv("DOSE_BOOT_N", "1000"))
    CORE_MODELS = ["LightGBM", "XGBoost", "RandomForest", "Logistic"]

    PROJECT = output_dir()
    OUT = Path(os.getenv('DOSE_ML_OUT', str(PROJECT / 'ukb/03ml')))
    PRED_FILE = OUT / "SelectionValidation_Predictions.csv"

    rng_master = np.random.default_rng(SEED)
    if N_BOOT < 1:
        raise ValueError('DOSE_BOOT_N must be positive')
    fingerprint = hashlib.sha256(Path(__file__).read_bytes())
    for package in ['numpy', 'pandas', 'scipy', 'scikit-learn', 'joblib']:
        fingerprint.update((package + importlib.metadata.version(package)).encode())
    cache = Memory(OUT / 'bootstrap_cache_v2' / fingerprint.hexdigest()[:16], verbose=0)
    if common.replace_results():
        cache.clear(warn=False)


    # ======================================================================
    # Utilities
    # ======================================================================

    def optimal_youden(y, score):
        fpr, tpr, thresholds = roc_curve(y, score)
        j = tpr - fpr
        i = int(np.nanargmax(j))
        return float(thresholds[i])


    def safe_div(a, b):
        return float(a / b) if b > 0 else np.nan


    def point_metrics(y, score, threshold):
        y = np.asarray(y, dtype=int)
        score = np.asarray(score, dtype=float)
        pred = (score >= threshold).astype(int)

        tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()

        sensitivity = safe_div(tp, tp + fn)
        specificity = safe_div(tn, tn + fp)
        precision = safe_div(tp, tp + fp)
        recall = sensitivity
        f1 = (
            2 * precision * recall / (precision + recall)
            if np.isfinite(precision + recall) and (precision + recall) > 0
            else np.nan
        )
        accuracy = safe_div(tp + tn, tp + tn + fp + fn)
        youden = sensitivity + specificity - 1
        balanced = np.nanmean([sensitivity, specificity])

        return {
            "AUC": roc_auc_score(y, score),
            "AUPRC": average_precision_score(y, score),
            "Accuracy": accuracy,
            "Sensitivity": sensitivity,
            "Specificity": specificity,
            "Precision": precision,
            "Recall": recall,
            "F1": f1,
            "Youden": youden,
            "BalancedAccuracy": balanced,
            "TN": int(tn),
            "FP": int(fp),
            "FN": int(fn),
            "TP": int(tp),
        }


    class ScoreCache:
        """
        Fast weighted AUC/AUPRC for bootstrap weights.

        Score order is fixed once; each bootstrap replicate only changes subject
        weights, avoiding an O(N log N) sort at every replicate.
        """

        def __init__(self, y, score):
            self.y = np.asarray(y, dtype=np.int8)
            self.score = np.asarray(score, dtype=float)
            self.n = len(self.y)

            self.order = np.argsort(-self.score, kind="mergesort")
            self.y_sorted = self.y[self.order]
            self.score_sorted = self.score[self.order]

            # End index of each unique-score block.
            self.ends = np.r_[
                np.flatnonzero(np.diff(self.score_sorted) != 0),
                self.n - 1,
            ]
            # Only score blocks containing cases contribute to AUC/AP sums.
            # Retain full blocks for tied positive/negative scores (exact, not an
            # approximation). Especially effective for rare outcomes such as IS.
            self.pos_idx = np.flatnonzero(self.y_sorted)
            group = np.searchsorted(self.ends, self.pos_idx)
            self.pos_groups, self.pos_group_index = np.unique(group, return_inverse=True)
            self.pos_ends = self.ends[self.pos_groups]
            starts = np.r_[0, self.ends[:-1] + 1]
            self.pos_starts = starts[self.pos_groups]

        def auc_ap(self, weights):
            w = np.asarray(weights, dtype=float)[self.order]

            cumulative = np.cumsum(w)
            pos_mass = np.bincount(self.pos_group_index, weights=w[self.pos_idx], minlength=len(self.pos_groups))
            tp = np.cumsum(pos_mass)
            pos = tp[-1] if len(tp) else 0.0
            neg = cumulative[-1] - pos
            if pos <= 0 or neg <= 0:
                return np.nan, np.nan

            total_at_end = cumulative[self.pos_ends]
            fp = total_at_end - tp
            total_before = np.zeros(len(self.pos_starts))
            mask = self.pos_starts > 0
            total_before[mask] = cumulative[self.pos_starts[mask] - 1]
            tied_neg = total_at_end - total_before - pos_mass
            auc = float(np.sum(pos_mass * (neg - fp + 0.5 * tied_neg)) / (pos * neg))
            precision = np.divide(tp, total_at_end, out=np.zeros_like(tp), where=total_at_end > 0)
            ap = float(np.sum(pos_mass * precision) / pos)
            return auc, ap


    def weighted_threshold_metrics(y, score, threshold, weights, categories=None):
        y = np.asarray(y, dtype=np.int8)
        score = np.asarray(score, dtype=float)
        w = np.asarray(weights, dtype=float)

        if categories is None:
            categories = 2 * y + (score >= threshold).astype(np.int8)
        tn, fp, fn, tp = np.bincount(categories, weights=w, minlength=4)

        sensitivity = safe_div(tp, tp + fn)
        specificity = safe_div(tn, tn + fp)
        precision = safe_div(tp, tp + fp)
        recall = sensitivity
        f1 = (
            2 * precision * recall / (precision + recall)
            if np.isfinite(precision + recall) and (precision + recall) > 0
            else np.nan
        )
        accuracy = safe_div(tp + tn, tp + tn + fp + fn)
        youden = sensitivity + specificity - 1
        balanced = np.nanmean([sensitivity, specificity])

        return {
            "Accuracy": accuracy,
            "Sensitivity": sensitivity,
            "Specificity": specificity,
            "Precision": precision,
            "Recall": recall,
            "F1": f1,
            "Youden": youden,
            "BalancedAccuracy": balanced,
        }


    def bootstrap_weights(n, rng):
        """
        Ordinary nonparametric bootstrap represented by integer subject weights.
        Equivalent to drawing n subject indices with replacement.
        """
        idx = rng.integers(0, n, size=n)
        return np.bincount(idx, minlength=n)


    @cache.cache
    def bootstrap_metrics(y, score, threshold, n_boot, seed):
        y = np.asarray(y, dtype=int)
        score = np.asarray(score, dtype=float)
        n = len(y)

        cache = ScoreCache(y, score)
        categories = 2 * y + (score >= threshold).astype(np.int8)
        rng = np.random.default_rng(seed)

        metrics = [
            "AUC", "AUPRC", "Accuracy", "Sensitivity", "Specificity",
            "Precision", "Recall", "F1", "Youden", "BalancedAccuracy"
        ]
        reps = {m: np.empty(n_boot, dtype=float) for m in metrics}

        for b in range(n_boot):
            w = bootstrap_weights(n, rng)
            auc, ap = cache.auc_ap(w)
            tm = weighted_threshold_metrics(y, score, threshold, w, categories)

            reps["AUC"][b] = auc
            reps["AUPRC"][b] = ap
            for m in tm:
                reps[m][b] = tm[m]

            if (b + 1) % max(100, n_boot // 10) == 0:
                print(f"  bootstrap {b + 1}/{n_boot}")

        return reps


    def summarize_bootstrap(point, reps):
        rows = []
        for metric, vals in reps.items():
            vals = vals[np.isfinite(vals)]
            rows.append({
                "Metric": metric,
                "Estimate": float(point[metric]),
                "Lower95": float(np.quantile(vals, 0.025)),
                "Upper95": float(np.quantile(vals, 0.975)),
                "BootstrapMedian": float(np.median(vals)),
                "n_boot_valid": int(len(vals)),
            })
        return pd.DataFrame(rows)


    # ======================================================================
    # Fast paired DeLong AUC test
    # ======================================================================

    def compute_midrank(x):
        return rankdata(np.asarray(x, dtype=float), method='average')


    def fast_delong(predictions_sorted_transposed, label_1_count):
        m = int(label_1_count)
        n = predictions_sorted_transposed.shape[1] - m
        k = predictions_sorted_transposed.shape[0]

        positive = predictions_sorted_transposed[:, :m]
        negative = predictions_sorted_transposed[:, m:]

        tx = np.empty((k, m))
        ty = np.empty((k, n))
        tz = np.empty((k, m + n))

        for r in range(k):
            tx[r, :] = compute_midrank(positive[r, :])
            ty[r, :] = compute_midrank(negative[r, :])
            tz[r, :] = compute_midrank(predictions_sorted_transposed[r, :])

        aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / (2.0 * n)

        v01 = (tz[:, :m] - tx) / n
        v10 = 1.0 - (tz[:, m:] - ty) / m

        sx = np.cov(v01)
        sy = np.cov(v10)
        delong_cov = sx / m + sy / n
        return aucs, delong_cov


    def delong_test(y, score_a, score_b):
        y = np.asarray(y, dtype=int)
        score_a = np.asarray(score_a, dtype=float)
        score_b = np.asarray(score_b, dtype=float)

        if not (len(y) == len(score_a) == len(score_b)):
            raise ValueError("DeLong inputs have different lengths.")

        # Positives first, as required by fast_delong.
        order = np.argsort(-y)
        y_ord = y[order]
        n_pos = int(y_ord.sum())

        preds = np.vstack([score_a[order], score_b[order]])
        aucs, cov = fast_delong(preds, n_pos)

        var = cov[0, 0] + cov[1, 1] - 2.0 * cov[0, 1]
        var = max(float(var), 0.0)

        if var == 0:
            z = np.nan
            p = np.nan
        else:
            z = float((aucs[0] - aucs[1]) / math.sqrt(var))
            p = float(2.0 * norm.sf(abs(z)))

        return float(aucs[0]), float(aucs[1]), z, p


    # ======================================================================
    # Paired bootstrap comparison
    # ======================================================================

    @cache.cache
    def paired_bootstrap_delta(y, score_a, score_b, n_boot, seed):
        y = np.asarray(y, dtype=int)
        score_a = np.asarray(score_a, dtype=float)
        score_b = np.asarray(score_b, dtype=float)
        n = len(y)

        cache_a = ScoreCache(y, score_a)
        cache_b = ScoreCache(y, score_b)
        rng = np.random.default_rng(seed)

        d_auc = np.empty(n_boot, dtype=float)
        d_ap = np.empty(n_boot, dtype=float)

        for b in range(n_boot):
            # Crucially, the SAME bootstrap weights are applied to both predictions.
            w = bootstrap_weights(n, rng)
            auc_a, ap_a = cache_a.auc_ap(w)
            auc_b, ap_b = cache_b.auc_ap(w)

            d_auc[b] = auc_a - auc_b
            d_ap[b] = ap_a - ap_b

            if (b + 1) % max(100, n_boot // 10) == 0:
                print(f"  paired bootstrap {b + 1}/{n_boot}")

        def summarise(x):
            x = x[np.isfinite(x)]
            p = 2.0 * min(
                (np.sum(x <= 0) + 1) / (len(x) + 1),
                (np.sum(x >= 0) + 1) / (len(x) + 1),
            )
            return (
                float(np.quantile(x, 0.025)),
                float(np.quantile(x, 0.975)),
                float(min(p, 1.0)),
            )

        return summarise(d_auc), summarise(d_ap)


    def aligned_pair(pred, feature_a, method_a, feature_b, method_b):
        a = pred.loc[
            (pred["FeatureSet"] == feature_a) & (pred["Method"] == method_a),
            ["eid", "y_true", "score"]
        ].rename(columns={"score": "score_a", "y_true": "y_a"})

        b = pred.loc[
            (pred["FeatureSet"] == feature_b) & (pred["Method"] == method_b),
            ["eid", "y_true", "score"]
        ].rename(columns={"score": "score_b", "y_true": "y_b"})

        z = a.merge(b, on="eid", how="inner", validate="one_to_one")
        if len(z) != len(a) or len(z) != len(b):
            raise RuntimeError(
                f"Prediction alignment failed: {feature_a}/{method_a} vs "
                f"{feature_b}/{method_b}; n={len(a)}, {len(b)}, merged={len(z)}"
            )
        if not np.array_equal(z["y_a"].to_numpy(), z["y_b"].to_numpy()):
            raise RuntimeError("Outcome labels differ after prediction alignment.")

        return (
            z["y_a"].to_numpy(dtype=int),
            z["score_a"].to_numpy(dtype=float),
            z["score_b"].to_numpy(dtype=float),
        )


    # ======================================================================
    # Main
    # ======================================================================

    if not PRED_FILE.exists():
        raise FileNotFoundError(PRED_FILE)

    required = {"eid", "y_true", "score", "Method", "FeatureSet"}
    started = time.perf_counter()
    pred = pd.read_csv(PRED_FILE, usecols=lambda c: c in required,
                       dtype={'eid': str, 'Method': 'category', 'FeatureSet': 'category'})
    print(f'Read predictions: {len(pred):,} rows in {time.perf_counter()-started:.1f}s', flush=True)
    missing = required - set(pred.columns)
    if missing:
        raise RuntimeError(f"Prediction file missing columns: {sorted(missing)}")
    if pred[list(required)].isna().any().any():
        raise ValueError('Missing values in required prediction columns; rerun selection validation')
    if not pred['y_true'].isin([0, 1]).all() or not np.isfinite(pred['score']).all():
        raise ValueError('Expected binary y_true and finite scores')
    if pred.duplicated(['FeatureSet', 'Method', 'eid']).any():
        raise ValueError('Duplicate subject within feature-set/model')
    for key, z in pred.groupby(['FeatureSet', 'Method'], observed=True):
        if z['y_true'].nunique() != 2:
            raise ValueError(f'Both classes are required: {key}')

    primary_sets = [
        x for x in pred["FeatureSet"].dropna().unique()
        if str(x).startswith("PrimaryDelta005_")
    ]
    all_sets = [
        x for x in pred["FeatureSet"].dropna().unique()
        if str(x).startswith("All") and str(x)[3:].isdigit()
    ]

    if len(primary_sets) != 1:
        raise RuntimeError(f"Expected exactly one primary feature set; found {primary_sets}")
    if len(all_sets) != 1:
        raise RuntimeError(f"Expected exactly one all-feature set; found {all_sets}")

    PRIMARY = primary_sets[0]
    ALLSET = all_sets[0]
    k_primary = len(pd.read_csv(OUT / "Primary_Selected_Features.csv"))
    if PRIMARY != f"PrimaryDelta005_{k_primary}":
        raise RuntimeError("Primary selection changed; rerun ./dose.sh ml --steps validate")
    feature_comparators = [ALLSET, f"SHAP_Consensus{k_primary}", f"RFE_Consensus{k_primary}"]
    missing_sets = set(feature_comparators) - set(pred["FeatureSet"].unique())
    if missing_sets:
        raise RuntimeError(f"Missing size-matched feature sets {missing_sets}; rerun ./dose.sh ml --steps validate")

    print("Bootstrap/paired evaluation")
    print("Primary feature set:", PRIMARY)
    print("All-candidate set:", ALLSET)
    print("Bootstrap replicates:", N_BOOT)


    # ----------------------------------------------------------------------
    # A. Bootstrap performance CIs for primary core models
    # ----------------------------------------------------------------------

    summary_rows = []
    cm_rows = []

    for j, method in enumerate(CORE_MODELS):
        z = pred.loc[
            (pred["FeatureSet"] == PRIMARY) & (pred["Method"] == method)
        ].copy()

        if z.empty:
            print(f"Skipping missing model: {method}")
            continue

        if z["eid"].duplicated().any():
            raise RuntimeError(f"Duplicated eid in {PRIMARY}/{method}")

        y = z["y_true"].to_numpy(dtype=int)
        score = z["score"].to_numpy(dtype=float)

        cutoff = optimal_youden(y, score)
        point = point_metrics(y, score, cutoff)

        print(f"\n=== bootstrap performance: {PRIMARY} | {method} ===")
        print(
            f"N={len(y)} events={y.sum()} cutoff={cutoff:.6g} "
            f"AUC={point['AUC']:.4f} AUPRC={point['AUPRC']:.4f}"
        )

        reps = bootstrap_metrics(
            y, score, cutoff, N_BOOT, SEED + 1000 + j
        )
        s = summarize_bootstrap(point, reps)
        s.insert(0, "Method", method)
        s.insert(1, "FeatureSet", PRIMARY)
        s["N"] = len(y)
        s["events"] = int(y.sum())
        s["YoudenCutoff"] = cutoff
        summary_rows.append(s)

        for true_label, pred_label, n in [
            (0, 0, point["TN"]),
            (0, 1, point["FP"]),
            (1, 0, point["FN"]),
            (1, 1, point["TP"]),
        ]:
            cm_rows.append({
                "Method": method,
                "FeatureSet": PRIMARY,
                "True": true_label,
                "Predicted": pred_label,
                "N": n,
                "YoudenCutoff": cutoff,
                "Sensitivity": point["Sensitivity"],
                "Specificity": point["Specificity"],
                "Accuracy": point["Accuracy"],
                "Precision": point["Precision"],
                "F1": point["F1"],
                "AUC": point["AUC"],
                "AUPRC": point["AUPRC"],
            })

    bootstrap_summary = pd.concat(summary_rows, ignore_index=True)
    bootstrap_summary.to_csv(
        OUT / "Bootstrap_Performance_Summary.csv", index=False
    )
    pd.DataFrame(cm_rows).to_csv(
        OUT / "Confusion_Matrix_Primary.csv", index=False
    )


    # ----------------------------------------------------------------------
    # B. Primary LightGBM vs alternative feature sets
    # ----------------------------------------------------------------------

    feature_cmp_rows = []

    for j, comparator in enumerate(feature_comparators):
        print(f"\n=== paired feature-set comparison: {PRIMARY} vs {comparator} | LightGBM ===")

        y, s_primary, s_comp = aligned_pair(
            pred,
            PRIMARY, "LightGBM",
            comparator, "LightGBM"
        )

        auc_p = roc_auc_score(y, s_primary)
        auc_c = roc_auc_score(y, s_comp)
        ap_p = average_precision_score(y, s_primary)
        ap_c = average_precision_score(y, s_comp)

        _, _, z_delong, p_delong = delong_test(
            y, s_primary, s_comp
        )

        (auc_l, auc_u, auc_p_boot), (ap_l, ap_u, ap_p_boot) = (
            paired_bootstrap_delta(
                y, s_primary, s_comp,
                N_BOOT,
                SEED + 2000 + j,
            )
        )

        feature_cmp_rows.append({
            "Method": "LightGBM",
            "Primary": PRIMARY,
            "Comparator": comparator,
            "N": len(y),
            "events": int(y.sum()),
            "AUC_Primary": auc_p,
            "AUC_Comparator": auc_c,
            "Delta_AUC": auc_p - auc_c,
            "Delta_AUC_Lower95_boot": auc_l,
            "Delta_AUC_Upper95_boot": auc_u,
            "DeLong_z": z_delong,
            "DeLong_P": p_delong,
            "Bootstrap_P_DeltaAUC": auc_p_boot,
            "AUPRC_Primary": ap_p,
            "AUPRC_Comparator": ap_c,
            "Delta_AUPRC": ap_p - ap_c,
            "Delta_AUPRC_Lower95_boot": ap_l,
            "Delta_AUPRC_Upper95_boot": ap_u,
            "Bootstrap_P_DeltaAUPRC": ap_p_boot,
            "n_boot": N_BOOT,
        })

    feature_cmp = pd.DataFrame(feature_cmp_rows)
    feature_cmp.to_csv(
        OUT / "Paired_FeatureSet_Comparison_LightGBM.csv",
        index=False
    )


    # ----------------------------------------------------------------------
    # C. Primary-set model comparisons versus LightGBM
    # ----------------------------------------------------------------------

    model_cmp_rows = []

    for j, comparator in enumerate(
        [m for m in CORE_MODELS if m != "LightGBM"]
    ):
        if comparator not in set(pred["Method"].unique()):
            continue

        print(f"\n=== paired model comparison: LightGBM vs {comparator} | {PRIMARY} ===")

        y, s_lgb, s_comp = aligned_pair(
            pred,
            PRIMARY, "LightGBM",
            PRIMARY, comparator
        )

        auc_lgb = roc_auc_score(y, s_lgb)
        auc_c = roc_auc_score(y, s_comp)
        ap_lgb = average_precision_score(y, s_lgb)
        ap_c = average_precision_score(y, s_comp)

        _, _, z_delong, p_delong = delong_test(
            y, s_lgb, s_comp
        )

        (auc_l, auc_u, auc_p_boot), (ap_l, ap_u, ap_p_boot) = (
            paired_bootstrap_delta(
                y, s_lgb, s_comp,
                N_BOOT,
                SEED + 3000 + j,
            )
        )

        model_cmp_rows.append({
            "FeatureSet": PRIMARY,
            "PrimaryModel": "LightGBM",
            "ComparatorModel": comparator,
            "N": len(y),
            "events": int(y.sum()),
            "AUC_LightGBM": auc_lgb,
            "AUC_Comparator": auc_c,
            "Delta_AUC": auc_lgb - auc_c,
            "Delta_AUC_Lower95_boot": auc_l,
            "Delta_AUC_Upper95_boot": auc_u,
            "DeLong_z": z_delong,
            "DeLong_P": p_delong,
            "Bootstrap_P_DeltaAUC": auc_p_boot,
            "AUPRC_LightGBM": ap_lgb,
            "AUPRC_Comparator": ap_c,
            "Delta_AUPRC": ap_lgb - ap_c,
            "Delta_AUPRC_Lower95_boot": ap_l,
            "Delta_AUPRC_Upper95_boot": ap_u,
            "Bootstrap_P_DeltaAUPRC": ap_p_boot,
            "n_boot": N_BOOT,
        })

    pd.DataFrame(model_cmp_rows).to_csv(
        OUT / "Paired_Model_Comparison_Primary.csv",
        index=False
    )


    # ----------------------------------------------------------------------
    # Console summary
    # ----------------------------------------------------------------------

    print("\nBootstrap evaluation complete")
    print("Outputs:")
    print(" -", OUT / "Bootstrap_Performance_Summary.csv")
    print(" -", OUT / "Confusion_Matrix_Primary.csv")
    print(" -", OUT / "Paired_FeatureSet_Comparison_LightGBM.csv")
    print(" -", OUT / "Paired_Model_Comparison_Primary.csv")

    print("\nPrimary model performance with bootstrap CIs:")
    show = bootstrap_summary.loc[
        bootstrap_summary["Metric"].isin(
            ["AUC", "AUPRC", "Accuracy", "Precision", "Recall", "F1"]
        ),
        ["Method", "Metric", "Estimate", "Lower95", "Upper95"]
    ]
    print(show.to_string(index=False))

    print("\nLightGBM paired feature-set comparisons:")
    print(feature_cmp[
        [
            "Comparator",
            "Delta_AUC",
            "Delta_AUC_Lower95_boot",
            "Delta_AUC_Upper95_boot",
            "DeLong_P",
            "Delta_AUPRC",
            "Delta_AUPRC_Lower95_boot",
            "Delta_AUPRC_Upper95_boot",
            "Bootstrap_P_DeltaAUPRC",
        ]
    ].to_string(index=False))


STEPS = {
    "compare": run_compare,
    "importance": run_importance,
    "select": run_select,
    "validate": run_validate,
    "bootstrap": run_bootstrap,
}

if __name__ == "__main__":
    common.run_cli(STEPS)
