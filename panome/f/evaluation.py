"""Evaluation starts only from frozen predictions and model metadata."""
from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from survival import metrics
from common import dump, log, VERSION


def paired_intervals(y, w, prediction, masks, groups, primary, count, seed):
    contrasts = [(primary, name) for name in ["elasticnet", "protein_elasticnet", "transformer_random_panel",
        "transformer_diversity_panel", "transformer_fullbank", "pca_same_panel", "embedding_knn_same_panel",
        "mlp_same_panel", "uniform_same_panel", "metric_same_panel", "equal_heads_same_panel", "no_pretrain_same_panel", "ssl_same_panel", "permuted_values_same_panel",
        "random_donors_same_panel", "copy1_reliable_panel"] if name in prediction]
    contrasts += [("panome_clinical", "clinical"), ("panome_clinical", "elasticnet")]
    contrasts += [("transformer_fullbank", "mlp_fullbank")] if "mlp_fullbank" in prediction else []
    rng = np.random.default_rng(seed)
    rows = []
    for subset, mask in masks.items():
        idx = np.flatnonzero(mask)
        if len(idx)<30 or len(np.unique(y[idx][w[idx]>0]))<2: continue
        units = [idx[groups[idx]==g] for g in np.unique(groups[idx])] if groups is not None else None
        names = sorted(set(v for pair in contrasts for v in pair))
        point = {name:metrics(y[idx], w[idx], prediction[name][idx]) for name in names}
        draws = {(m,b,key):[] for m,b in contrasts for key in ["AUC_IPCW", "Brier_IPCW", "LogLoss_IPCW"]}
        for _ in range(count):
            take = np.concatenate([units[j] for j in rng.integers(0,len(units),len(units))]) if units is not None else rng.choice(idx,len(idx),replace=True)
            value = {name:metrics(y[take], w[take], prediction[name][take]) for name in names}
            for (m,b,key), vals in draws.items():
                delta = value[m][key]-value[b][key]
                if np.isfinite(delta): vals.append(delta)
        for (m,b,key), vals in draws.items():
            ci = np.quantile(vals,[.025,.975]) if len(vals)>=20 else [np.nan,np.nan]
            rows.append(dict(subset=subset,model=m,reference=b,metric=key,
                delta=point[m][key]-point[b][key],lower=ci[0],upper=ci[1],replicates=len(vals),
                uncertainty="conditional_on_fitted_model; exploratory_unadjusted_intervals"))
    return pd.DataFrame(rows)


def evaluate_frozen(out):
    if not (out/"MODEL_FROZEN.json").exists(): raise ValueError("Missing frozen-model record")
    bundle = joblib.load(out/"model_bundle.joblib")
    cfg, primary = bundle["config"], bundle["primary"]
    people = pd.read_csv(out/"test_outcomes.csv", dtype={cfg["id_col"]:str})
    info = pd.read_csv(out/"test_individuals.csv", dtype={cfg["id_col"]:str})
    if not people[cfg["id_col"]].equals(info[cfg["id_col"]]): raise ValueError("Test outcome/prediction ID mismatch")
    names = json.loads((out/"prediction_columns.json").read_text())
    predictions = {name:info[name].to_numpy() for name in names}
    y,w = bundle["censoring"].labels_weights(people)
    accepted = info.supported_match.to_numpy(bool)
    masks = {"all":np.ones(len(y),bool), "supported":accepted, "rejected":~accepted}
    rows = []
    for subset, mask in masks.items():
        for name, prob in predictions.items():
            rows.append(dict(model=name, subset=subset, n=int(mask.sum()), known=int((w[mask]>0).sum()),
                cases=int(y[mask].sum()), coverage=float(mask.mean()), **metrics(y[mask],w[mask],prob[mask])))
    result = pd.DataFrame(rows); result.to_csv(out/"test_metrics.csv", index=False)
    groups = people[cfg["group_col"]].to_numpy(str) if cfg["group_col"] else None
    paired_intervals(y,w,predictions,{key:masks[key] for key in ["all","supported"]},groups,
                     primary,cfg["bootstrap"],cfg["seed"]+823).to_csv(out/"paired_contrasts.csv",index=False)
    curves = []
    for quantile, rule in bundle["coverage_rules"].items():
        mask,_ = rule.apply(info,info.missing_fraction.to_numpy(),info.unknown_category.to_numpy(bool),cfg["sample_missing"])
        for name, prob in predictions.items():
            curves.append(dict(quantile=quantile,model=name,coverage=float(mask.mean()),n=int(mask.sum()),
                               **metrics(y[mask],w[mask],prob[mask])))
    pd.DataFrame(curves).to_csv(out/"coverage_curve.csv",index=False)
    # Compare expected-error scores with actual squared error without claiming
    # that lower outcome prevalence alone proves better matching.
    errbins = pd.qcut(info.estimated_squared_error, 5, labels=False, duplicates="drop")
    error_rows = []
    for b in sorted(errbins.dropna().unique()):
        mask = errbins.eq(b).to_numpy()
        error_rows.append(dict(bin=int(b)+1,n=int(mask.sum()),expected_error=float(info.loc[mask,"estimated_squared_error"].mean()),
            observed_error=float(np.average((y[mask]-predictions[primary][mask])**2,weights=w[mask])),
            event_rate_IPCW=float(np.average(y[mask],weights=w[mask])),
            brier_gain_vs_elasticnet=float(np.mean(w[mask]*((y[mask]-predictions["elasticnet"][mask])**2-(y[mask]-predictions[primary][mask])**2)))))
    pd.DataFrame(error_rows).to_csv(out/"support_error_audit.csv",index=False)
    calibration = []
    for name,prob in predictions.items():
        bins = pd.qcut(pd.Series(prob),5,labels=False,duplicates="drop")
        for b in sorted(bins.dropna().unique()):
            mask = bins.eq(b).to_numpy()
            calibration.append(dict(model=name,bin=int(b)+1,n=int(mask.sum()),
                predicted=float(prob[mask].mean()),observed_IPCW=float(np.average(y[mask],weights=w[mask])) if w[mask].sum()>0 else np.nan))
    pd.DataFrame(calibration).to_csv(out/"test_calibration.csv",index=False)
    demographic = pd.read_csv(out/"test_demographics.csv")
    if "age" in demographic:
        demographic["age_group"] = pd.cut(demographic.age,[0,50,60,70,150],right=False).astype(str)
    subgroup = []
    for column in [c for c in ["sex","center","age_group"] if c in demographic]:
        for level,sub in demographic.groupby(column,dropna=False):
            ix = sub.index.to_numpy()
            for name in [primary,"elasticnet","panome_clinical"]:
                subgroup.append(dict(variable=column,level=level,model=name,n=len(ix),cases=int(y[ix].sum()),
                    supported_fraction=float(accepted[ix].mean()),**metrics(y[ix],w[ix],predictions[name][ix])))
    pd.DataFrame(subgroup).to_csv(out/"subgroup_metrics.csv",index=False)
    raw_rows = []
    for name in names:
        col = name+"_raw"
        if col in info:
            raw_rows.append(dict(model=name,**metrics(y,w,info[col].to_numpy())))
    pd.DataFrame(raw_rows).to_csv(out/"raw_test_metrics.csv",index=False)
    compare = result[result.subset.eq("all")].copy()
    for name,bank in bundle["banks"].items():
        ix = compare.model.eq(name)
        compare.loc[ix,"reference_count"] = len(bank.ids)
        compare.loc[ix,"reference_case_fraction"] = bank.y.mean()
    compare.to_csv(out/"approach_comparison.csv",index=False)
    write_report(out,bundle,result,info)
    figures(out,bundle,result)
    dump(out/"DONE.json",dict(version=VERSION,n_test=len(y),cases_by_horizon=int(y.sum()),
        supported_fraction=float(accepted.mean()),panel_status=bundle["readiness"]["status"],
        results_synthetic=cfg["demo"]))
    log("DONE","evaluation",f"N={len(y)}; horizon events={y.sum()}; support={accepted.mean():.1%}")


def write_report(out,bundle,result,info):
    primary = bundle["primary"]
    all_rows = result[result.subset.eq("all")].sort_values("AUC_IPCW",ascending=False)
    text = ["# Panome 5 run report", "", "SYNTHETIC demonstration; not UKB evidence." if bundle["config"]["demo"] else "Research results on the supplied input cohort.", "",
        f"Prespecified primary: `{primary}`. Reference readiness: `{bundle['readiness']['status']}`.",
        f"Readiness reasons: {bundle['readiness']['reasons']}", "",
        "| Model | AUC IPCW | Brier IPCW |", "|---|---:|---:|"]
    for row in all_rows.itertuples(): text.append(f"| {row.model} | {row.AUC_IPCW:.4f} | {row.Brier_IPCW:.5f} |")
    text += ["", "## Questions to answer before claiming a useful individual reference model", "",
        "1. Does the learned encoder + copying beat same-reference PCA, Euclidean copying, MLP, and random reference controls?",
        "2. Does a 100-person panel preserve useful signal from the full bank? Does class stratification prevent constant COPY?",
        "3. Do reference label permutation and fixed-neighbour deletion measurably change predictions? A head-only predictor is not evidence of useful outcome borrowing.",
        "4. Does masked reconstruction beat a zero build-mean profile and PCA reconstruction on exactly the same masked entries? Stable but constant predictions are not informative explanations.",
        "5. Does support improve error relative to a comparator, or mostly reject higher-risk people? Read subgroup_metrics and support_error_audit together.",
        "6. Are same-risk molecular contrasts reproducible across seeds/splits/residualization, and externally validated? Token names are not biological pathway claims.", "",
        "## Individual interpretation", "",
        "The uncalibrated estimate is exactly the sum of nonnegative reference outcome contributions and an explicit build-prior contribution. Monotone calibration follows this sum. COPY1 outputs an observed binary outcome, not a person's certain future or a probability estimate.",
        "Deleting a reference holds the retrieved neighbour set fixed and renormalizes the remaining weights. Module masking changes information availability; neither operation is a causal intervention.",
        "Reference identities alone do not establish explanatory validity. Use the recorded contribution identity, deletion sensitivity, missing-protein reconstruction and controls.", "",
        "## Validation limits", "",
        "The outer test half is not used for fitting, model selection, calibration, panel utility or support thresholds. Calibration is split into fitting and independent audit halves. Inner neural validation is only an early-stopping device; OOF reference-fit scores come from separate cross-fitting of the configured quality teachers.",
        "Report all prespecified models. Test-set ranking is exploratory: selecting a winner after seeing this report needs a new untouched cohort. Bootstrap intervals condition on the fitted model and are unadjusted for multiple contrasts; use the repeat runner for training instability.",
        "A readiness pass is an audit point-estimate screen, not proof of individual accuracy. A failed panel does not invalidate the source cohort. Death is censored, so the endpoint is fixed-horizon net risk, not a competing-risk cumulative incidence. IPCW uses a marginal build-set censoring model and assumes independent censoring.",
        ""]
    (out/"REPORT.md").write_text("\n".join(text),encoding="utf-8")


def figures(out,bundle,result):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    main = result[result.subset.eq("all")].sort_values("AUC_IPCW")
    fig,ax = plt.subplots(figsize=(10,max(6,.27*len(main))))
    colors = ["#c04a31" if s==bundle["primary"] else "#4a718e" for s in main.model]
    ax.barh(main.model,main.AUC_IPCW,color=colors)
    ax.axvline(.5,color="grey",ls="--"); ax.set_xlim(0,1); ax.set_xlabel("Full test IPCW AUC")
    fig.tight_layout(); fig.savefig(out/"Fig_model_comparison.png",dpi=180); plt.close(fig)
    if (out/"masked_reconstruction.csv").exists():
        dat = pd.read_csv(out/"masked_reconstruction.csv").groupby("method").masked_RMSE.mean().sort_values()
        fig,ax = plt.subplots(figsize=(8,4)); ax.barh(dat.index,dat.values,color="#4a718e")
        ax.set_xlabel("Mean person-level masked RMSE (lower is better)")
        fig.tight_layout(); fig.savefig(out/"Fig_masked_reconstruction.png",dpi=180); plt.close(fig)
    dat = pd.read_csv(out/"coverage_curve.csv")
    fig,ax = plt.subplots(figsize=(7,4))
    for name in [bundle["primary"],"elasticnet","panome_clinical"]:
        sub = dat[dat.model.eq(name)]
        ax.plot(sub.coverage,sub.Brier_IPCW,"o-",label=name)
    ax.set_xlabel("Supported fraction; thresholds frozen on tune")
    ax.set_ylabel("IPCW Brier on identical people"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out/"Fig_coverage.png",dpi=180); plt.close(fig)
