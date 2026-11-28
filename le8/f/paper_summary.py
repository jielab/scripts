"""Build manuscript tables and figures from saved aggregate evidence."""
import argparse
from html import escape
import os
from pathlib import Path
from shutil import copy2
import re
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dx_compare import TIMING, report
from dx_layout import DATES, READ_DATES, PRIMARY_DX, scope_for_read
from file_lock import file_lock
from paper_focus import paper_focus, focus_sections, highlights, highlight_figures, primary_sections, primary_figures

KEY = ["trait", "layer", "dx"]
MODULE_FILES = {
    "c1_correlate": ("c1.cohort.xlsx", "evaluation.xlsx"),
    "c2_cause": ("c2.mr.xlsx", "c2.genetic_leadtime.xlsx"),
    "c3_coloc": ("c3.coloc.xlsx",),
    "c4_connect": ("c4.cohort.xlsx", "c4.explain.xlsx", "c4.imaging.xlsx"),
    "c5_cellulation": ("c5.cellulation.xlsx",),
}
P1 = "P1_incident_support"
P0 = "P0_risk_only"
RAW = "raw_count_support"


def subset(df, **filters):
    if not set(filters) <= set(df):
        return df.iloc[:0]
    for col, values in filters.items():
        df = df[df[col].isin(values)]
    return df


def stack(frames):
    use = [f for f in frames if not f.empty]
    return pd.concat(use, ignore_index=True, sort=False) if use else pd.DataFrame()


def is_dated(layer, dates):
    return any((layer/"main/results"/d).is_dir() for d in dates)


def load_inputs(root, requested, layers, dates):
    from current_release import require_current
    require_current(root)
    traits = requested or sorted(p.name for p in root.iterdir() if p.is_dir()
        and re.fullmatch(r"[A-Za-z0-9_]+",p.name)
        and ((p/"aod/main/results/dx_counts.tsv").is_file()
             or any(is_dated(p/b,dates) for b in layers)))
    if not traits:
        raise ValueError("No current analysis scopes found. Run le8.sh first.")
    live = report(root,traits,layers,dates)
    inventory=[]
    for t in traits:
        for b in layers:
            for d in dates:
                scope=scope_for_read(root/t/b,d)
                for mod,files in MODULE_FILES.items():
                    if mod!="c1_correlate" and d!=PRIMARY_DX: continue
                    for name in files:
                        f=scope/mod/name
                        if f.is_file():
                            inventory.append(dict(trait=t,layer=b,dx=d,module=mod,file=str(f.relative_to(root))))
    live["availability"]["source"]=["current scope" if scope_for_read(root/t/b,d).is_dir() else "no eligible scope"
        for t,b,d in live["availability"][KEY].itertuples(index=False,name=None)]
    live["module_files"]=pd.DataFrame(inventory)
    return live


def date_auc_comparison(performance, paired):
    cols = ["trait","layer","dx","model","horizon","AUC_IPCW","evaluated_N","cases"]
    required = {"trait","layer","reference","model","metric","analysis","difference","lower","upper",
                "reference_value","model_value","test_known_N","test_cases"}
    if not set(cols) <= set(performance) or not required <= set(paired):
        return pd.DataFrame()
    rows = paired.loc[(paired.analysis=="saved_C1_common_test") & (paired.metric=="AUC_IPCW")].copy()
    if rows.empty:
        return pd.DataFrame()
    rows["reference_dx"] = rows.reference.str.split(".",n=1,regex=False).str[0]
    rows["method"] = rows.model.str.split(".",n=1,regex=False).str[1]
    rows = rows.rename(columns={"difference":"paired_delta_AUC","lower":"paired_lower","upper":"paired_upper",
                                "reference_value":"common_reference_AUC","model_value":"common_any1_AUC"})
    native = performance.loc[performance.horizon==10,cols].rename(columns={"model":"method"})
    for which in ("reference","any1"):
        target = "reference_dx" if which=="reference" else "dx"
        right = native.rename(columns={"dx":target,**{c:f"{which}_{c}" for c in ("AUC_IPCW","evaluated_N","cases")}})
        rows = rows.merge(right.drop(columns="horizon"),on=["trait","layer","method",target],how="left",validate="many_to_one")
    rows["own_cohort_delta_AUC"] = rows.any1_AUC_IPCW-rows.reference_AUC_IPCW
    return rows[["trait","layer","method","reference_dx","reference_AUC_IPCW","any1_AUC_IPCW","own_cohort_delta_AUC",
                 "reference_evaluated_N","reference_cases","any1_evaluated_N","any1_cases","test_known_N","test_cases",
                 "common_reference_AUC","common_any1_AUC","paired_delta_AUC","paired_lower","paired_upper"]]


def select_tables(data):
    performance = subset(data["performance"],subset=["all_common"])
    contrasts = subset(data["model_contrasts"],subset=["all_common"],same_people=[True],
                       contrast=["historical_increment_same_architecture","molecular_residual_increment_same_clinical_background","model_effect_same_people"])
    if not contrasts.empty:
        pairs = {("abm_prevalent","abm_prevalent_control"),("abm_clinical","clinical"),
                 ("abm","glmnet"),("abm_deep","glmnet"),("lightgbm","glmnet")}
        contrasts = contrasts.loc[[pair in pairs for pair in zip(contrasts.model,contrasts.reference)]]
        if "evaluation_roster_hash" in performance and "model_roster_hash" in contrasts:
            den = performance[[*KEY,"model","horizon","evaluated_N","cases","evaluation_roster_hash"]].drop_duplicates()
            contrasts = contrasts.merge(den,how="left",left_on=[*KEY,"model","model_roster_hash"],
                        right_on=[*KEY,"model","evaluation_roster_hash"],validate="many_to_one")
    flows = data["transitions"]
    reclass = pd.DataFrame()
    if not flows.empty:
        rows = []
        for (t,a,b),group in flows.groupby(["trait","dx_a","dx_b"],sort=False):
            for state,src in group.groupby("state_a",sort=False):
                total = src.N.sum()
                for r in src.itertuples():
                    rows.append(dict(trait=t,dx_a=a,dx_b=b,state_a=state,state_b=r.state_b,N=r.N,
                                     source_state_N=total,pct_of_source_state=100*r.N/total if total else np.nan))
        reclass = pd.DataFrame(rows)
    washout = []
    if not reclass.empty:
        for trait, rows in reclass.groupby("trait", sort=False):
            initial = subset(rows, dx_a=["icd"], dx_b=["ref"], state_a=["incident"])
            checked = subset(rows, dx_a=["ref"], dx_b=["icd1"], state_a=["prevalent"])
            if initial.empty or checked.empty:
                continue
            discordant = initial.loc[initial.state_b.eq("prevalent"), "N"].sum()
            retained = checked.loc[checked.state_b.eq("incident"), "N"].sum()
            excluded = checked.loc[checked.state_b.eq("excluded_year1"), "N"].sum()
            # ICD1 changes only ICD incident diagnoses on post-baseline days 1–365.
            if retained + excluded != discordant:
                raise ValueError(f"{trait}: incompatible ICD, REF and ICD1 transition counts")
            washout.append(dict(trait=trait, icd_incident_N=initial.N.sum(),
                                ref_prevalent_among_icd_incident=discordant,
                                excluded_by_icd1=excluded, still_incident_in_icd1=retained,
                                pct_still_incident_of_discordant=100*retained/discordant if discordant else np.nan))
    rosters = pd.DataFrame()
    cols = [*KEY,"horizon","evaluated_N","cases","evaluation_roster_hash"]
    if set(cols) <= set(performance):
        rosters = performance[cols].drop_duplicates()
        # all_common must identify one roster per layer/date/horizon across every model.
        if rosters.duplicated([*KEY,"horizon"]).any():
            raise ValueError("Conflicting all_common rosters or event counts in the same scope")
        anchor = rosters[rosters.dx.eq("icd")].drop(columns="dx").rename(columns={
            "evaluated_N":"icd_test_N","cases":"icd_cases","evaluation_roster_hash":"icd_roster_hash"})
        rosters = rosters.merge(anchor,on=["trait","layer","horizon"],how="left",validate="many_to_one")
        known = rosters.evaluation_roster_hash.notna() & rosters.icd_roster_hash.notna()
        rosters["test_membership_vs_icd"] = np.select(
            [~known,rosters.evaluation_roster_hash.eq(rosters.icd_roster_hash),rosters.evaluated_N.ne(rosters.icd_test_N)],
            ["unavailable","identical ordered IDs","different size"],default="different IDs or ordering")
    # Older portable workbooks retain the primary selector estimates and CIs
    # in wide form even when the full risk-matched table was not exported.
    summary = data["ABM_support"]
    gradients, summary_tests = [], []
    for row in summary.to_dict("records"):
        label = {k:row[k] for k in KEY}
        for col in ("twin_count_Brier_slope","family_ESS_Brier_slope","P1_minus_P0_50pct_Brier","P1_minus_P0_50pct_LogLoss"):
            ci = re.fullmatch(r"\s*\[([^,]+),\s*([^\]]+)\]\s*",str(row.get(col+"_CI","")))
            if not ci or pd.isna(row.get(col)):
                continue
            item = dict(**label,estimate=col,difference=row[col],lower=float(ci[1]),upper=float(ci[2]))
            if col.endswith("_slope"):
                gradients.append(item)
            elif row.get("primary_conclusion") not in (None,"primary_unavailable"):
                summary_tests.append(dict(**item,model=P1,reference=P0,risk_model="glmnet",metric=col.rsplit("_",1)[1],
                    risk_bins=10,primary_risk_bins=True,primary_risk_model=True,coverage_target=.5,
                    table_source="ABM_support summary; detailed selection counts unavailable"))
    tests = data["selector_tests"]
    if summary_tests:
        keys = [*KEY,"model","reference","risk_model","risk_bins","metric"]
        tests = stack([tests,pd.DataFrame(summary_tests)]).drop_duplicates(keys,keep="first")
    return {"availability":data["availability"],"cohort_counts":data["cohort_counts"],"reclassification":reclass,
            "first_year_exclusion":pd.DataFrame(washout),
            "glmnet_path":data["glmnet_path"],
            "date_auc_comparison":date_auc_comparison(performance,data["timing_paired"]),
            **{name:data[name] for name in (*TIMING,"any1_sources","any1_discordance")},
            "ham_sources":data["ham_sources"],"historical_borrowing":data["historical_borrowing"],"assay_counts":data["assay_counts"],"all_common":performance,
            "test_rosters":rosters,"paired_models":contrasts,"selector_tests":tests,"support_gradients":pd.DataFrame(gradients),"ABM_support":data["ABM_support"],
            "glmnet_neighbors":data["glmnet_neighbors"],"glmnet_neighbor_status":data["glmnet_neighbor_status"],
            "coverage":subset(data["coverage"],coverage_target=[.5]),"model_status":data["model_status"],
            "LDL_interaction":subset(data["LDL_interaction"],dx=[PRIMARY_DX,"default"]),"LDL_slopes":subset(data["LDL_slopes"],dx=[PRIMARY_DX,"default"]),"module_files":data["module_files"]}


def make_figures(tables, stage):
    names = []
    counts = tables["cohort_counts"]
    if {"incident","prevalent","excluded_year1","uncertain"} <= set(counts) and len(counts):
        traits = list(counts.trait.unique())
        flows = tables["reclassification"]
        with_flows = len(subset(flows,dx_a=["icd"],dx_b=["ref","any1","ham"],state_a=["incident"]))>0
        fig, panels = plt.subplots(len(traits),2 if with_flows else 1,
                                   figsize=(12 if with_flows else 8,3*len(traits)),squeeze=False,layout="constrained")
        for panel,trait in zip(panels,traits):
            ax = panel[0]
            df = counts[counts.trait.eq(trait)].set_index("dx")
            df = df.reindex([d for d in READ_DATES if d in df.index])
            df[["incident","prevalent","excluded_year1","uncertain"]].plot.bar(ax=ax,color=["#237d86","#b45736","#dbc17e","#80868b"],rot=0)
            ax.set(title=trait, xlabel="",ylabel="Participants")
            ax.legend(ncol=2,fontsize=8)
            if with_flows:
                right=panel[1]
                moved=subset(flows,trait=[trait],dx_a=["icd"],dx_b=["ref","any1","ham"],state_a=["incident"])
                known=moved.dx_b.unique()
                changes=moved[moved.state_b.ne("incident")]
                if len(changes):
                    graph=changes.pivot(index="dx_b",columns="state_b",values="pct_of_source_state").reindex(known).fillna(0)
                    colors={"prevalent":"#b45736","uncertain":"#80868b","control":"#237d86"}
                    graph.plot.barh(stacked=True,ax=right,color=[colors.get(c,"#dbc17e") for c in graph])
                    right.legend(fontsize=8)
                right.set(title="Reclassified among ICD incident cases",ylabel="",xlabel="% of ICD incident cases")
        names.append("Fig1_disease_dates.png")
        fig.savefig(stage/names[-1],dpi=190); plt.close(fig)
    slopes = subset(tables["support_gradients"],estimate=["twin_count_Brier_slope"]).assign(model="twin_count")
    for name, df, title, axis in (
        ("Fig2_Yang_prediction.png",subset(tables["paired_models"],model=["abm_prevalent"],
             reference=["abm_prevalent_control"],metric=["AUC_IPCW"],contrast=["historical_increment_same_architecture"]),
             "Increment from prevalent-case information","Change in IPCW AUC"),
        ("Fig3_twin_support.png",subset(tables["selector_tests"],model=[P1,RAW],metric=["Brier"],
              primary_risk_bins=[True],primary_risk_model=[True],risk_model=["glmnet"]),
             "Prediction error at 50% target coverage","Change in risk-standardized Brier score"),
        ("Fig4_neighbour_gradient.png",slopes,"Molecular neighbour count and prediction error",
             "Brier slope per unit log(1 + eligible neighbour count)"),
        ("Fig5_clinical_increment.png",subset(tables["paired_models"],model=["abm_clinical"],reference=["clinical"],
             metric=["AUC_IPCW"],contrast=["molecular_residual_increment_same_clinical_background"]),
             "Clinical risk with molecular-neighbour residual correction","Change in IPCW AUC")):
        if not {"difference","lower","upper"} <= set(df):
            continue
        values = df[["difference","lower","upper"]].apply(pd.to_numeric,errors="coerce")
        df = df.loc[np.isfinite(values).all(axis=1)].copy()
        if df.empty:
            continue
        order = [d for d in (*READ_DATES,"default") if d in set(df.dx)]
        df["dx"] = pd.Categorical(df.dx,order,ordered=True)
        df = df.sort_values(["trait","layer","dx","model"])
        scopes = list(df.groupby(["trait","layer"],sort=False))
        fig, panels = plt.subplots(len(scopes),1,figsize=(8, max(2.8,2.7*len(scopes))),squeeze=False,layout="constrained")
        for (ax,),((trait,layer),part) in zip(panels,scopes):
            for i, (method,g) in enumerate(part.groupby("model",sort=False)):
                y = g.dx.cat.codes.to_numpy() + (i*.24-.12 if name.startswith("Fig3") else 0)
                color = "#237d86" if i==0 else "#ba613a"
                ax.hlines(y,g.lower,g.upper,color=color)
                ax.scatter(g.difference,y,color=color,s=19,label={P1:"P1 minus P0",RAW:"Raw support minus risk-matched","twin_count":"Risk-adjusted error slope","abm_clinical":"Corrected minus clinical"}.get(method,"Yang minus no history"))
            ax.axvline(0,color="0.5",linewidth=.8)
            ax.set(yticks=range(len(order)),yticklabels=order,ylabel=f"{trait} / {layer}",xlabel=axis)
            ax.invert_yaxis(); ax.legend(fontsize=8,loc="best")
        fig.suptitle(title)
        fig.savefig(stage/name,dpi=190); plt.close(fig); names.append(name)
    score=subset(tables["timing_scores"],model=["omics"])
    if not score.empty:
        parts=list(score.groupby(["trait","layer"],sort=False))
        fig,panels=plt.subplots(len(parts),1,figsize=(10, max(3,3.5*len(parts))),squeeze=False,layout="constrained")
        order=["prevalent_anchor","early_with_prior","early_with_undated","early_day_1_90","early_day_91_180",
               "early_day_181_365","later_incident_anchor","uncertain_timing","control"]
        for (ax,),((trait,layer),part) in zip(panels,parts):
            part=part.set_index("profile_group").reindex(order).dropna(subset=["median"])
            y=np.arange(len(part))
            ax.hlines(y,part.q25,part.q75,color="#688d92")
            ax.scatter(part["median"],y,s=24,color="#215f66")
            ax.set(yticks=y,yticklabels=part.index,xlim=(0,1),xlabel="Held-out molecular similarity to prevalent anchors (median and IQR)",title=f"{trait} / {layer}")
            ax.invert_yaxis()
        names.append("Fig6_timing_profiles.png")
        fig.savefig(stage/names[-1],dpi=190);plt.close(fig)
    paired=subset(tables["timing_paired"],analysis=["controlled_training"],metric=["AUC_IPCW","Brier_IPCW"])
    if not paired.empty:
        fig,axes=plt.subplots(1,2,figsize=(10,max(3,len(paired)*.32)),layout="constrained")
        for ax,(metric,title) in zip(axes,(("AUC_IPCW","10-year IPCW AUC"),("Brier_IPCW","10-year IPCW Brier"))):
            part=paired[paired.metric.eq(metric)].sort_values(["trait","layer"])
            y=np.arange(len(part))
            ax.axvline(0,c="0.6",linewidth=1)
            ax.hlines(y,part.lower,part.upper,color="#237d86");ax.scatter(part.difference,y,s=22,color="#237d86")
            ax.set(yticks=y,yticklabels=part.trait+" / "+part.layer,title=title,xlabel="Filtered minus inclusive training")
            ax.invert_yaxis()
        names.append("Fig7_paired_training.png")
        fig.savefig(stage/names[-1],dpi=190);plt.close(fig)
    pairs=tables.get("focus_contrasts",pd.DataFrame())
    if not pairs.empty:
        panels=[("H1_remove_ambiguous","clean_all","clean/clinical_omics_ridge","Omit ambiguous training records"),
                ("H2_remove_low_support_train","clean_high","clean_high/twins","Omit low-support training participants")]
        fig,axes=plt.subplots(1,2,figsize=(12,5),layout="constrained")
        for ax,(analysis,test,model,title) in zip(axes,panels):
            df=pairs[(pairs.metric=="AUC_IPCW") & (pairs.analysis==analysis) & (pairs.test_set==test) & (pairs.model==model)]
            df=df.sort_values(["trait","layer"])
            y=np.arange(len(df))
            ax.hlines(y,df.lower,df.upper,color="#227a76",linewidth=2)
            ax.scatter(df.difference,y,color="#163c50",s=25,zorder=3)
            ax.axvline(0,color="0.55",lw=1,ls="--")
            ax.set_yticks(y,df.trait+" / "+df.layer);ax.invert_yaxis()
            ax.set_title(title,fontsize=11);ax.set_xlabel("Change in 10-year IPCW AUC (paired 95% CI)")
        name="Fig8_timing_and_twin_training.png"
        fig.savefig(stage/name,dpi=190);plt.close(fig);names.append(name)
    return names


def html_report(tables, figures):
    body = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>LE8 | Disease dates and molecular twins</title>',
            '<style>body{font:15px/1.55 system-ui;color:#24343b;max-width:1200px;margin:3em auto;padding:0 1.2em}.scroll{overflow:auto}th,td{padding:.3em .55em;border-bottom:1px solid #ccd5d7;text-align:right}table{border-collapse:collapse;font-variant-numeric:tabular-nums}img{max-width:100%;width:780px}h2{margin-top:2em}</style>',
            '<h1>Disease timing and the predictive value of molecular twins</h1>',
            '<p><a href="../results/paper_tables.xlsx">Download tables with denominators, intervals and source availability</a></p>',
            '<p>ICD uses the recorded ICD-10 date. REF uses the reference first occurrence. ICD1 and REF1 exclude cases diagnosed on days 1–365 after baseline. Day zero counts as prevalent. Any1 excludes cases with a first diagnosis date in days 1–365 in any eligible source, including discordant baseline histories. Other documented baseline cases are prevalent. The Any1 risk set requires observation beyond day 365. The default C1 report compares clinical and molecular neighbor increments within ICD and within Any1. Eligibility can differ across endpoints. C2–C5 use ICD once per trait and assay.</p>']
    def section(title, df, cols):
        if df.empty:return
        body.append("<h2>"+escape(title)+"</h2>")
        cols = [c for c in cols if c in df]
        body.append('<div class="scroll">'+df[cols].to_html(index=False,border=0,na_rep="",float_format=lambda x:f"{x:.5g}")+"</div>"
                    if len(df) and cols else "<p>No eligible estimates available for this selection.</p>")
    if "focus_four_groups" in tables or "focus_predictivity" in tables:
        body.append("<p>The first-year ICD-10 cohort is pooled across days 1–365. REF, GP and SRD can overlap. The workbook gives distinct counts, source combinations and the measured distribution of each analyte in focus_four_features. The similarity score is trained on ICD timing anchors. Higher values indicate greater resemblance to the prevalent anchor.</p>")
        body.append("<p>Per-person predictivity is one minus estimated squared prediction error for the clinical plus molecular ridge model. The frozen estimator uses observed errors in other development families. Compare predicted and realized Brier loss in the outer test. Negative values favor the local estimator or the higher-predictivity group in the two contrasts, respectively. Intervals condition on the fitted predictors.</p>")
    primary_sections(tables,section)
    for name in figures:
        body.append(f'<p><a href="../results/{name}"><img src="../results/{name}" alt="{escape(name)}"></a></p>')
    if len(tables.get("focus_contrasts",[])):
        body.append("<p>The paired training comparisons change training inclusion while holding the evaluation participants, ICD labels, follow-up horizon and IPCW weights fixed. Twin counts use independent reference families. The high-support threshold is frozen on a separate outcome-blind development split. Baseline molecular profile scores measure similarity to recorded timing anchors. They are not adjudicated diagnoses. Family intervals condition on the fitted models.</p>")
    focus_sections(tables,section)
    section("Endpoint counts before assay matching",tables["cohort_counts"],["trait","dx","cohort_N","control","incident","prevalent","excluded_year1","uncertain","short_followup","unavailable"])
    f = tables["reclassification"]
    if len(f):
        section("Who changed classification relative to ICD?",f[(f.dx_a=="icd") & (f.state_a!=f.state_b)],
                list(f.columns))
    section("Among ICD incident cases with a pre-baseline or baseline REF date, how many survive a first-year exclusion?",
            tables["first_year_exclusion"],list(tables["first_year_exclusion"].columns))
    section("Fate of ICD incident cases with baseline REF records",tables["any1_discordance"],
            ["trait","icd1_state","any1_state","N","source_N","pct_of_discordant"])
    section("Any-source first dates and competing histories",tables["any1_sources"],
            ["trait","any1_state","any_prior","any_undated","any_first_source","any_early_sources","N"])
    body.append("<p>The controlled ridge comparison changes only training inclusion of first-year cases. Both arms predict the same outcome through year 10 from baseline in the same participants eligible at day 365. Intervals use paired family resampling conditional on the fitted models. Saved C1 contrasts also use aligned held-out participants and Any1 labels, but other training differences between definitions can contribute.</p>")
    section("Change in 10-year AUC on each definition’s cohort and on the same eligible people with any1 labels",
            tables["date_auc_comparison"],["trait","layer","method","reference_dx","own_cohort_delta_AUC","test_known_N",
            "test_cases","common_reference_AUC","common_any1_AUC","paired_delta_AUC","paired_lower","paired_upper"])
    section("Matched held-out prediction after training exclusions",tables["timing_paired"],
            [*KEY,"analysis","reference","model","metric","reference_value","model_value","difference","lower","upper","test_known_N","test_cases"])
    section("Discrimination of recorded timing anchors",tables["timing_anchor_validation"],
            [*KEY,"model","prevalent_N","later_incident_N","heldout_anchor_AUC","prevalent_sensitivity_at_half","later_incident_specificity_at_half"])
    section("Blinded molecular similarity in the ambiguous groups",tables["timing_scores"],
            [*KEY,"model","profile_group","N","q25","median","q75","fraction_above_half"])
    section("Frozen omics profile by baseline lipid-lowering medication",tables["timing_medication"],
            [*KEY,"baseline_lipid_treatment","profile_group","N","median","q25","q75","heldout_anchor_AUC","prevalent_anchor_N","later_incident_anchor_N"])
    section("Timing analysis coverage and training counts",tables["timing_status"],
            [*KEY,"analysis","status","training_full_N","training_filtered_N","early_training_cases_removed","test_eligible_N","test_known_N"])
    section("Molecular cohorts",tables["assay_counts"],[*KEY,*[c for c in tables["assay_counts"] if c not in [*KEY,"endpoint_sha256"]]])
    body.append("<p>The history contrast uses the same incident test participants and the same risk-model architecture. The reported AUC horizon and event count refer to that scope. A positive AUC difference or negative Brier difference favors the model using prevalent-case information.</p>")
    section("Adding prevalent-case information to risk prediction",subset(tables["paired_models"],model=["abm_prevalent"],reference=["abm_prevalent_control"],contrast=["historical_increment_same_architecture"],metric=["AUC_IPCW","Brier_IPCW"]),[*KEY,"horizon","model_N","cases","metric","difference","lower","upper","inference_status"])
    section("How much did the fitted model borrow?",tables["historical_borrowing"],[*KEY,"alpha","source_offered_N","source_used_N","source_families","history_use_status","status","reason"])
    section("Performance on the common test roster",tables["all_common"],[*KEY,"model","horizon","evaluated_N","cases","controls","coverage_all","Harrell_C","AUC_IPCW","Brier_IPCW","LogLoss_IPCW"])
    section("Test roster compared with the ICD anchor",tables["test_rosters"],[*KEY,"horizon","evaluated_N","cases","icd_test_N","icd_cases","test_membership_vs_icd"])
    body.append("<p>P1 versus P0 asks whether incident twin support helps select participants with lower glmnet error beyond baseline predicted risk. Raw support versus a risk-matched selection is a secondary comparison. P2 is disabled in the current CVD selector protocol. Bootstrap intervals condition on the fitted models and selectors.</p>")
    section("Selection at 50% target coverage",subset(tables["selector_tests"],risk_bins=[10],risk_model=["glmnet"]),[*KEY,"model","reference","risk_model","metric","model_selected_N","model_cases","reference_weight_fraction","difference","lower","upper","inference_status"])
    section("Twin count and prediction error within risk strata",subset(tables["support_gradients"],estimate=["twin_count_Brier_slope"]),[*KEY,"difference","lower","upper"])
    section("Protocol-level selector results",tables["ABM_support"],[*KEY,"P1_status","P2_status","P1_minus_P0_50pct_Brier","P1_minus_P0_50pct_Brier_CI","primary_conclusion"])
    section("LDL genetic-component × residual interaction",tables["LDL_interaction"],[*KEY,"analysis","start_year","end_year","events","clustering_unit","cluster_N","beta_GxR","ci_low","ci_high","p_interaction","p_interval_holm","status"])
    section("LASSO penalty-path audit",tables["glmnet_path"],[*KEY,"path_status","solver_jerr","maxit_used","coordinate_passes","solver_seconds","lambda_count","selected_index","lambda","tune_logloss","reason"])
    section("Availability",tables["availability"],[*KEY,"source","execution","prediction","input_matches_aod","abm_deep_status"])
    body.append("</html>")
    return "\n".join(body)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis-root",type=Path,required=True)
    p.add_argument("--Y")
    p.add_argument("--biom",default="prot,met")
    p.add_argument("--dx-date",default="icd,any1")
    a=p.parse_args()
    root=a.analysis_root.resolve()
    out=root/"paper"
    for part in ("main/results","main/reports","qc","tmp","jobs","logs"):
        (out/part).mkdir(parents=True,exist_ok=True)
    with file_lock(out/"qc/summary.lock"), tempfile.TemporaryDirectory(dir=out/"tmp") as work:
        stage=Path(work)
        data=load_inputs(root,a.Y.split(",") if a.Y else None,a.biom.split(","),a.dx_date.split(","))
        tables=select_tables(data)
        traits = a.Y.split(",") if a.Y else sorted(set(tables["availability"].get("trait",[])))
        tables.update(paper_focus(root,traits,a.biom.split(","),a.dx_date.split(",")))
        if all(frame.empty for key,frame in tables.items() if key not in ("availability","focus_sources")):
            raise ValueError("No eligible scientific tables found for the requested scope; previous paper retained.")
        tables.update(highlights(tables))
        figs=primary_figures(tables,stage)+make_figures(tables,stage)+highlight_figures(tables,stage)
        with pd.ExcelWriter(stage/"paper_tables.xlsx",engine="openpyxl") as book:
            for name,frame in tables.items():
                frame.to_excel(book,sheet_name=name,index=False,freeze_panes=(1,0))
        (stage/"index.html").write_text(html_report(tables,figs),encoding="utf-8")
        tables["availability"].to_csv(stage/"availability.tsv",sep="\t",index=False)
        publish = {out/"main/results"/name:stage/name for name in ["paper_tables.xlsx",*figs]}
        # Publish the current gallery only. A scoped report leaves figures from other runs on disk.
        publish[out/"qc/availability.tsv"] = stage/"availability.tsv"
        publish[out/"main/reports/index.html"] = stage/"index.html"
        backup = stage/"previous"
        backup.mkdir()
        touched = []
        try:
            for target,source in publish.items():
                if target.is_file():
                    try:
                        os.link(target,backup/target.name)
                    except OSError:
                        copy2(target,backup/target.name)
                if source is None:
                    target.unlink(missing_ok=True)
                else:
                    os.replace(source,target)
                touched.append(target)
        except BaseException:
            for target in reversed(touched):
                previous = backup/target.name
                if previous.is_file():
                    os.replace(previous,target)
                else:
                    target.unlink(missing_ok=True)
            raise
        message=f"{len(tables['availability'])} requested scopes; {len(tables['paired_models'])} paired estimates; {len(figs)} figures."
        with (out/"logs/paper.log").open("a") as log:
            print(message,file=log)
        print(message, out/"main/reports/index.html",sep="\n")


if __name__=="__main__":
    main()
