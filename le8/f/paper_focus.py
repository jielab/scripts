"""Read the paired C1 experiments and selected aggregate native C2–C5 exports."""
from pathlib import Path
import numpy as np
from file_lock import file_lock
from zipfile import ZipFile
import json
import pandas as pd

from dx_layout import scope_for_read

EXPORTS = {
 "c2_cause":{"c2.subanalysis_status.csv":"C2_status"},
 "c3_coloc":{"c3.subanalysis_status.csv":"C3_status","c3.coloc_summary.csv":"C3_loci"},
 "c4_connect":{"c4.focus.contrasts.csv":"C4_panel_contrasts","c4.focus.proxy_accuracy.csv":"C4_proxy_accuracy",
               "c4.focus.deployed_concept_fidelity.csv":"C4_concept_fidelity"},
 "c5_cellulation":{"c5.cellulation_status.csv":"C5_status","c5.CIGMA_enrichment.csv":"C5_enrichment"}
}


def aggregate_exports(folder,targets):
    remaining=dict(targets)
    for name in list(remaining):
        path=folder/name
        if path.is_file():
            yield remaining.pop(name),pd.read_csv(path),name
    for book in sorted(folder.glob("*.xlsx")):
        if not remaining:break
        with ZipFile(book) as z:
            if "le8/manifest.json" not in z.namelist():continue
            files=json.loads(z.read("le8/manifest.json")).get("files",{})
            for name in list(remaining):
                if name in files:
                    with z.open(files[name]["part"]) as stream:
                        try: df=pd.read_csv(stream)
                        except pd.errors.EmptyDataError: df=pd.DataFrame()
                    yield remaining.pop(name),df,book.name+"/"+name


def native_comparison(base,layer,provenance):
    from timing_study import compare, native_predictions
    private=base/"tmp/c1_timing_twins"/(layer+"_risk_LOCAL_ONLY.tsv.gz")
    from c1_protocol import sha
    if not private.is_file() or json.loads(provenance.private_hashes).get(private.name)!=sha(private):return pd.DataFrame()
    hashes={"icd":provenance.icd_sha256}
    pred=native_predictions(base/layer,"icd","glmnet",hashes)
    if pred is None:return pd.DataFrame()
    d=pd.read_csv(private,sep="\t",dtype={"eid":str,"group":str})
    d["glmnet"]=d.eid.map(pred)
    out=[]
    for tag,mask in (("clean_all",d.clean),("clean_high",d.clean & d.high)):
        ok=mask & d.glmnet.between(0,1) & np.isfinite(d.twins)
        sub=d[ok]
        if len(sub) and min(sub.y.sum(),len(sub)-sub.y.sum())>=5:
            rows=compare(sub.y.to_numpy(),sub.w.to_numpy(),sub.group.to_numpy(),
                         sub.glmnet.to_numpy(),sub.twins.to_numpy(),300,"native_vs_focus",
                         "native_ICD_glmnet","clean/twins")
            out.extend(dict(row,test_set=tag,focus_known_N=int(mask.sum())) for row in rows)
    return pd.DataFrame(out)


def paper_focus(root,traits,layers,dates):
    if "icd" not in dates:return {}
    rows={}
    def add(key,df,trait,view):
        rows.setdefault(key,[]).append(df.assign(trait=trait,layer=view,dx="icd"))
    for trait in traits:
        base=Path(root)/trait
        manifest=base/"aod/qc/dx_manifest.json"
        fingerprint=json.loads(manifest.read_text()).get("endpoint_sha256") if manifest.is_file() else None
        for layer in [*layers,*(["prot_met"] if {"prot","met"}<=set(layers) else [])]:
            file=base/"main/results/c1_timing_twins"/(layer+".xlsx")
            status="absent"
            if file.is_file():
                with file_lock(base/"qc/c1_timing_twins/study.lock"), pd.ExcelFile(file) as book:
                    prov=book.parse("provenance").iloc[0]
                    status="matched" if (prov.trait==trait and prov.layer==layer and fingerprint and
                                         prov.snapshot_sha256==fingerprint) else "unverified_snapshot"
                    if status=="matched":
                        for key in book.sheet_names:add("focus_"+key,book.parse(key),trait,layer)
                        if layer in ("prot","met"):
                            paired=native_comparison(base,layer,prov)
                            add("focus_native_glmnet",paired,trait,layer)
                            add("focus_sources",pd.DataFrame([dict(source="native vs focus paired scores",
                                status="available" if len(paired) else "missing_or_insufficient_common_outcomes")]),trait,layer)
            add("focus_sources",pd.DataFrame([dict(source="paired timing/twins",status=status)]),trait,layer)
            if layer=="prot_met":continue
            scope=scope_for_read(base/layer,"icd")
            for mod,targets in EXPORTS.items():
                # The module workbooks belong to a native ICD scope validated against AOD.
                binding=scope/"qc/dx.json"
                if not binding.is_file() or not fingerprint:continue
                meta=json.loads(manifest.read_text())
                endpoint=next((r["endpoint_sha256"] for r in meta["counts"] if r["dx"]=="icd"),None)
                if json.loads(binding.read_text()).get("endpoint_sha256")!=endpoint:continue
                for key,df,source in aggregate_exports(scope/mod,targets):
                    if key=="C4_panel_contrasts" and {"stratum","landmark"}<=set(df):
                        df=df[df.stratum.eq("All") & df.landmark.eq(0)]
                    add(key,df.assign(source_file=source),trait,layer)
    return {k:pd.concat(v,ignore_index=True) for k,v in rows.items()}


def primary_sections(tables,section):
    for key,title,fields in (
        ("four_groups","1. First ICD-10 diagnosis in days 1–365: baseline molecular profiles by evidence of prior disease",
         ["group","N_all_assayed","N_outer_test","score_median","score_q25","score_q75","AUC_vs_no_prior","AUC_lower","AUC_upper"]),
        ("four_anchors","ICD-10 reference groups used to interpret molecular similarity",
         ["test_prevalent","prevalent_median","test_later","later_median","anchor_AUC","anchor_lower","anchor_upper","status"]),
        ("predictivity","2. Does expected error at baseline identify people with lower observed 10-year prediction error?",
         ["band","eligible_N","known_N","cases","mean_risk","mean_expected_Brier","observed_Brier","neighbor_coverage"]),
        ("predictivity_test","Does local support add information within similar predicted risk?",
         ["contrast","estimate","lower","upper","N","cases","families","shared_risk_bins","valid_bootstraps","status"])):
        section(title,tables.get("focus_"+key,pd.DataFrame()),["trait","layer",*fields])


def primary_figures(tables,stage):
    import matplotlib.pyplot as plt
    names=[]
    profiles=tables.get("focus_four_groups",pd.DataFrame())
    if "score_median" in profiles and profiles.score_median.notna().any():
        views=list(profiles.groupby(["trait","layer"],sort=True))
        nc=min(3,len(views));nr=(len(views)+nc-1)//nc
        fig,axes=plt.subplots(nr,nc,figsize=(4.6*nc,2.65*nr),squeeze=False,layout="constrained")
        groups=["no_prior_any_source","REF_prior","GP_prior","SRD_prior"]
        for ax,((trait,assay),frame) in zip(axes.flat,views):
            d=frame.set_index("group").reindex(groups);y=np.arange(len(groups))
            ax.hlines(y,d.score_q25,d.score_q75,color="#356e87",lw=2)
            ax.plot(d.score_median,y,"o",color="#143a50")
            ax.set_yticks(y,["No prior record","REF prior","GP prior","SRD prior"])
            ax.set_xlim(0,1);ax.invert_yaxis();ax.set_title(f"{trait} · {assay}")
            ax.set_xlabel("Prevalent-anchor similarity, median and IQR")
        for ax in list(axes.flat)[len(views):]:ax.set_visible(False)
        file="Fig12_ICD_first_year_profiles.png";fig.savefig(stage/file,dpi=180);plt.close(fig);names.append(file)
    tests=tables.get("focus_predictivity_test",pd.DataFrame())
    if "estimate" in tests and tests.estimate.notna().any():
        fields=(("local_minus_risk_only_error_prediction_MSE","Local minus risk-bin error estimator","Difference in loss-prediction MSE"),
                ("higher_minus_lower_observed_Brier_within_risk","Higher versus lower predictivity","Risk-bin standardized Brier difference"))
        n=tests.groupby("contrast").size().max()
        fig,axes=plt.subplots(1,2,figsize=(11,2.6+.30*n),squeeze=False,layout="constrained")
        for ax,(key,title,xlabel) in zip(axes.flat,fields):
            d=tests[tests.contrast.eq(key)].sort_values(["trait","layer"]); y=np.arange(len(d))
            ax.axvline(0,ls="--",lw=1,color=".6")
            ax.hlines(y,d.lower,d.upper,color="#356e87");ax.plot(d.estimate,y,"o",color="#143a50")
            ax.set_yticks(y,[f"{r.trait} · {r.layer}" for r in d.itertuples()])
            ax.invert_yaxis();ax.set_title(title);ax.set_xlabel(xlabel)
        file="Fig13_individual_predictivity.png";fig.savefig(stage/file,dpi=180);plt.close(fig);names.append(file)
    return names


def focus_sections(tables,section):
    fields=["trait","layer","test_set","analysis","reference","model","difference","lower","upper","test_known_N","test_cases"]
    pairs=tables.get("focus_contrasts",pd.DataFrame())
    if len(pairs):
        primary=pairs.metric.eq("AUC_IPCW") & (
            (pairs.analysis.eq("H1_remove_ambiguous") & pairs.test_set.eq("clean_all") & pairs.model.eq("clean/clinical_omics_ridge")) |
            (pairs.analysis.eq("H2_remove_low_support_train") & pairs.test_set.eq("clean_high") & pairs.model.eq("clean_high/twins")))
        section("Paired primary tests on identical held-out participants",pairs[primary],fields)
        incr=pairs[(pairs.metric=="AUC_IPCW") & pairs.analysis.eq("strong_model_neighbor_increment") &
                   pairs.test_set.isin(["clean_all","clean_high"])]
        section("Does local correction improve a strong clinical plus molecular ridge predictor?",incr,fields)
    for key,title,cols in (
        ("neighbor_models","Clinical and molecular model performance",["trait","layer","dx","model","evaluated_N","cases","AUC_IPCW","Brier_IPCW"]),
        ("clinical_neighbor_increment","Clinical risk with local correction",["trait","layer","dx","metric","difference","lower","upper","model_N"]),
        ("glmnet_neighbor_status","Molecular glmnet residual donors",["trait","layer","dx","status","donor_role","donor_N","donor_families","donor_cases","radius","test_supported_N"]),
        ("glmnet_neighbors","Increment above frozen molecular glmnet",["trait","layer","dx","subset","metric","reference_value","model_value","difference","lower","upper","eligible_N","test_known_N","test_cases","coverage"]),
        ("focus_icd_entry_flow","Entry into the first-year ICD-10 cohort",["trait","layer","split","stage","N"]),
        ("focus_icd_first_tests","Metabolite score in the first-ICD cohort",["trait","layer","comparison","baseline_lipid","estimator","positive_N","reference_N","families","estimate","lower","upper","status"]),
        ("focus_icd_first_scores","First-ICD score and upper-tail fraction by prior evidence",["trait","layer","period","history","baseline_lipid","N","median","q25","q75","high_N","fraction_high"]),
        ("focus_icd_any1_overlap","Which Any1 early-with-prior participants meet first-ICD entry?",["trait","layer","saved_ICD_first","ICD_any_version_prior_or_day0","first_year_sources","other_history","N"]),
        ("focus_original_scores","Metabolites in the Any1 timing groups",["trait","layer","profile_group","N","median","q25","q75"]),
        ("focus_early_median_contrast","Any1 early records by prior history",["trait","layer","stratum","prior_N","no_prior_N","prior_median","no_prior_median","median_difference","lower","upper"]),
        ("focus_early_adjusted_contrast","Adjusted mean score difference in Any1 early cases",["trait","layer","adjustment","N","families","prior_N","mean_score_difference","lower","upper"]),
        ("focus_native_glmnet","Native molecular glmnet and clean-trained twins on the same people",
          ["trait","layer","test_set","metric","reference_value","model_value","difference","lower","upper","test_known_N","test_cases"]),
        ("focus_selection","Does molecular twin support identify lower prediction error at similar clinical risk?",
          ["trait","layer","eligible_set","model","difference","lower","upper","high_known_N","low_known_N","high_cases","low_cases","overlap_strata"]),
        ("focus_timing_anchor_validation","Can baseline omics distinguish the two recorded timing anchors?",
          ["trait","layer","model","prevalent_N","later_incident_N","heldout_anchor_AUC"]),
        ("focus_timing_scores","Molecular similarity in the excluded groups",["trait","layer","model","profile_group","N","median","q25","q75"]),
        ("focus_twin_geometry","Prospective twin rule and coverage",["trait","layer","variables","reference_families","valid_reference_families","threshold","gate_coverage"]),
        ("C4_panel_contrasts","LE8 panel comparisons at equal molecular budgets",["trait","layer","model","reference","contrast","horizon","delta_AUC","delta_lo","delta_hi","valid"]),
        ("C2_status","Genetic analyses and input availability",["trait","layer","analysis","status","reason","detail","message"]),
        ("C3_status","Colocalization and input availability",["trait","layer","analysis","status","reason","detail","message"]),
        ("C5_status","Cell annotations and evidence origin",["trait","layer","analysis","status","detail"])):
        df=tables.get(key,pd.DataFrame())
        if len(df):
            if key in ("focus_timing_scores","focus_original_scores"):df=df[df.model.eq("omics")]
            section(title,df,cols)


def highlights(tables):
    """Use only paired model results already read for this report's live scopes."""
    perf=tables["all_common"]
    pairs=tables["paired_models"]
    out={}
    if not perf.empty:
        out["neighbor_models"]=perf[perf.model.isin(["clinical","abm_clinical","glmnet"])]
    if not pairs.empty:
        out["clinical_neighbor_increment"]=pairs[pairs.model.eq("abm_clinical") & pairs.reference.eq("clinical")]
    return out


def highlight_figures(tables,stage):
    import matplotlib.pyplot as plt
    names=[]
    def panels(groups,height):
        cols=min(3,len(groups))
        fig,axes=plt.subplots((len(groups)+cols-1)//cols,cols,squeeze=False,
                             figsize=(4.8*cols,height*((len(groups)+cols-1)//cols)),layout="constrained")
        for ax in list(axes.flat)[len(groups):]:
            ax.set_visible(False)
        return fig,axes.flat

    scores=tables.get("focus_original_scores",pd.DataFrame())
    if len(scores):
        views=list(scores[scores.model.eq("omics")].groupby("trait",sort=True))
        order=["later_incident_anchor","early_day_1_90","early_day_91_180","early_day_181_365","early_with_prior","prevalent_anchor"]
        labels=["Later incident","Days 1–90","Days 91–180","Days 181–365","Early with prior","Prevalent"]
        if views:
            fig,axes=panels(views,3.2)
            for ax,(trait,frame) in zip(axes,views):
                d=frame.set_index("profile_group").reindex(order)
                y=np.arange(len(d));ax.hlines(y,d.q25,d.q75,color="#356e87");ax.plot(d["median"],y,"o",color="#143a50")
                ax.set_yticks(y,[f"{label} (n={int(n):,})" if pd.notna(n) else label for label,n in zip(labels,d.N)])
                ax.set_xlim(0,1);ax.invert_yaxis();ax.set_title(str(trait));ax.set_xlabel("Prevalent-anchor similarity, median and IQR")
            file="Fig9_early_metabolites.png";fig.savefig(stage/file,dpi=190);plt.close(fig);names.append(file)

    icd=tables.get("focus_icd_first_scores",pd.DataFrame())
    if len(icd):
        chosen=icd[(icd.baseline_lipid=="all") & ((icd.period=="first_year_ICD") |
              ((icd.period=="later_ICD") & (icd.history=="no_prior_record_found")))]
        views=list(chosen.groupby("trait",sort=True))
        if views:
            fig,axes=panels(views,2+max(len(d) for _,d in views)*.32)
            for ax,(trait,d) in zip(axes,views):
                labels=[f"{r.history.replace('_',' ')} (n={r.N:,})" if r.period=="first_year_ICD"
                        else f"Later ICD, no prior (n={r.N:,})" for r in d.itertuples()]
                y=np.arange(len(d));ax.barh(y,d.fraction_high,color="#326579")
                ax.set_yticks(y,labels,fontsize=8);ax.invert_yaxis();ax.set_xlim(0,1);ax.set_title(str(trait))
                ax.set_xlabel("Fraction above later-ICD 90th percentile")
            file="Fig11_first_ICD.png";fig.savefig(stage/file,dpi=190);plt.close(fig);names.append(file)

    rows=[]
    for label,key in (("Clinical","clinical_neighbor_increment"),("Molecular glmnet","glmnet_neighbors")):
        d=tables.get(key,pd.DataFrame())
        if d.empty:continue
        d=d[d.metric.eq("AUC_IPCW") & d.subset.eq("all_common")]
        if len(d):rows.append(d.assign(background=label))
    if rows:
        contrasts=pd.concat(rows,ignore_index=True)
        views=list(contrasts.groupby(["layer","background"],sort=True))
        names_trait=sorted(contrasts.trait.unique());y=np.arange(len(names_trait))
        fig,axes=panels(views,2+.49*len(names_trait))
        for ax,((layer,background),frame) in zip(axes,views):
            ax.axvline(0,ls="--",lw=1,color=".65")
            dates=sorted(frame.dx.unique(),key=lambda x:(x!="icd",x!="any1",x))
            for n,dx in enumerate(dates):
                d=frame[frame.dx.eq(dx)].set_index("trait").reindex(names_trait)
                offset=(n-(len(dates)-1)/2)*(.65/max(1,len(dates)))
                color=plt.get_cmap("tab10")(n%10)
                ax.hlines(y+offset,d.lower,d.upper,color=color,lw=1.5)
                ax.plot(d.difference,y+offset,"o",color=color,label=dx,ms=4)
            ax.set_yticks(y,names_trait);ax.invert_yaxis();ax.set_title(f"{layer} · {background}")
            ax.legend(frameon=False);ax.set_xlabel("Paired change in 10-year IPCW AUC (95% CI)")
        file="Fig10_neighbor_increment.png";fig.savefig(stage/file,dpi=190);plt.close(fig);names.append(file)
    return names
