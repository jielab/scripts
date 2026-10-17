#!/usr/bin/env python3
"""Generate manuscript-ready C1 figures/tables from frozen C1/ABM outputs.

Outputs are descriptive/reporting only: this script never fits, selects, tunes or
changes a model. Main figures correspond exactly to the locked Results 1-3.
"""
from __future__ import annotations
import argparse, gzip, hashlib, io, json, math, zipfile
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from c1_lock_roster import grouping_audit, GROUPING_SOURCE

RELEASE="2026-10-17.c1-lockpaper-v1"

def read_table(path):
    path=Path(path)
    data = path.read_bytes() if path.is_file() else None
    if data is None:
        for book in sorted(path.parent.glob("*.xlsx")):
            with zipfile.ZipFile(book) as archive:
                if "le8/manifest.json" not in archive.namelist(): continue
                entry=json.loads(archive.read("le8/manifest.json")).get("files",{}).get(path.name)
                if not entry: continue
                data=archive.read(entry["part"])
                if hashlib.sha256(data).hexdigest()!=entry["sha256"]:
                    raise ValueError("Corrupt frozen table: " + str(path))
                break
    if data is None: return pd.DataFrame()
    if path.suffix==".gz": data=gzip.decompress(data)
    try:
        return pd.read_csv(io.BytesIO(data),sep="\t" if ".tsv" in path.name else ",")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

def savefig(fig,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"),bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"),dpi=300,bbox_inches="tight")
    plt.close(fig)

def fig1_conceptual(out):
    fig,ax=plt.subplots(figsize=(14,8));ax.set_xlim(0,14);ax.set_ylim(0,8);ax.axis("off")
    def box(x,y,w,h,title,body,rounding=.12):
        p=FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.03,rounding_size={rounding}",
                         fill=False,linewidth=1.5)
        ax.add_patch(p);ax.text(x+.18,y+h-.35,title,fontsize=13,fontweight="bold",va="top")
        ax.text(x+.18,y+h-.82,body,fontsize=10,va="top",wrap=True)
    def arr(a,b):
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=15,linewidth=1.3))
    ax.text(.2,7.65,"Figure 1. From population-average risk to individual molecular predictability",
            fontsize=16,fontweight="bold",va="top")
    box(.3,4.7,3.7,2.2,"A  Population-average prediction",
        "Conventional glmnet / LightGBM\nlearn one global mapping\nX → disease risk\n\nPopulation summary: AUC / C-index")
    box(5.0,4.7,3.7,2.2,"B  Molecular neighbourhood",
        "Each participant is embedded in an\nomics-defined neighbourhood.\nReference agents are evaluated for\nsupport, stability and reliability.")
    box(10.0,4.7,3.7,2.2,"C  Dual individual output",
        "Disease risk: 17%\nMolecular predictability: HIGH\n\nRisk answers: how high?\nPredictability answers: how reliable?")
    arr((4.05,5.8),(4.9,5.8));arr((8.75,5.8),(9.9,5.8))
    box(.3,1.0,4.0,2.2,"D  Locked validation design",
        "One outcome-blind outer test roster\nshared by conventional C1 and ABM.\nABM development is split internally into\nbuild / tune / gate / calibration roles.")
    box(5.0,1.0,4.0,2.2,"E  Primary predictability estimand",
        "At 20/40/50/60/80/100% coverage:\ncompare frozen-selector participants with\npredicted-risk–matched outcome-blind controls.\nPrimary: ΔBrier / ΔLogLoss AURC.")
    box(9.7,1.0,4.0,2.2,"F  Reference-agent specificity",
        "Twin distance: molecular-only (primary).\nClinical+omics geometry: sensitivity.\nTest whether cleaning helps local twin inference\nbut not global population regression.")
    arr((4.3,2.1),(4.9,2.1));arr((9.0,2.1),(9.6,2.1))
    ax.text(.3,.42,"Conceptual figure; values shown for illustration only. No test outcome enters individual selection.",
            fontsize=9)
    savefig(fig,Path(out)/"Fig1_conceptual_framework")

def _benchmark(perf,layer):
    d=perf.copy()
    if d.empty or not {"stage","family","model","metric","value"}<=set(d.columns):
        return pd.DataFrame()
    d=d[(d.stage=="c1_correlate")&(d.family=="C1 conventional holdout")]
    if d.empty:return d
    preferred="Pradeep-style / glmnet" if layer=="prot" else "All-metabolite / glmnet"
    if "panel" not in d: return pd.DataFrame()
    q=d[d.panel==preferred]
    rows=[]
    for metric in ["IPCW AUC","Harrell C"]:
        z=q[q.metric==metric]
        for model in ["Model 0","Biomarkers","Combined"]:
            x=z[z.model==model]
            if len(x)>1: raise ValueError(f"Ambiguous frozen benchmark: {layer}/{preferred}/{model}/{metric}")
            if len(x): rows.append(x.iloc[0])
    return pd.DataFrame(rows)

def fig2_population(root,out):
    rows=[]
    for layer in ["prot","met"]:
        perf=read_table(Path(root)/layer/"c1_abm"/"performance"/"performance.csv")
        b=_benchmark(perf,layer)
        if len(b):
            b=b.copy();b["layer"]=layer;rows.append(b)
    d=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
    if d.empty:return d
    auc=d[d.metric=="IPCW AUC"].copy()
    fig,ax=plt.subplots(figsize=(9,5.5))
    layers=["prot","met"];models=["Model 0","Biomarkers","Combined"];x=np.arange(len(layers));w=.23
    for j,m in enumerate(models):
        vals=[]
        for lay in layers:
            z=auc[(auc.layer==lay)&(auc.model==m)]
            vals.append(float(z.value.iloc[0]) if len(z) else np.nan)
        ax.bar(x+(j-1)*w,vals,width=w,label={"Model 0":"Clinical","Biomarkers":"Omics","Combined":"Clinical + omics"}[m])
    ax.set_xticks(x,["Proteomics","Metabolomics"]);ax.set_ylabel("10-year IPCW AUC")
    ax.set_title("Figure 2. Population-average prediction on the shared holdout")
    ax.legend(frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(.5,-.24))
    finite=pd.to_numeric(auc.value,errors="coerce").dropna()
    if len(finite): ax.set_ylim(min(.5,finite.min()-.03),min(1,finite.max()+.03))
    savefig(fig,Path(out)/"Fig2_population_average_prediction")
    return d

def fig3_predictability(root,out):
    allcurve=[];allsum=[];strata=[]
    for layer in ["prot","met"]:
        base=Path(root)/layer/"c1_abm"/"baseline_comparison"
        curve=read_table(base/"lockpaper_selector_curve.tsv.gz")
        summ=read_table(base/"lockpaper_selector_inference.tsv.gz")
        if len(curve):
            z=curve[curve.route=="risk_resid_error"].copy();z["layer"]=layer;allcurve.append(z)
        if len(summ):
            z=summ[summ.route=="risk_resid_error"].copy();z["layer"]=layer;allsum.append(z)
        met=read_table(Path(root)/layer/"c1_abm"/"abm_twin"/"test_report"/"metrics.tsv.gz")
        if len(met):
            z=met[(met.model=="elasticnet_full")&met.subset.isin(["selected","all","rejected"])].copy()
            z=z.drop_duplicates(["model","subset"]);z["layer"]=layer;strata.append(z)
    c=pd.concat(allcurve,ignore_index=True) if allcurve else pd.DataFrame()
    s=pd.concat(allsum,ignore_index=True) if allsum else pd.DataFrame()
    st=pd.concat(strata,ignore_index=True) if strata else pd.DataFrame()
    if c.empty:return c,s
    fig,axs=plt.subplots(1,3,figsize=(16,5))
    if len(st):
        order=["selected","all","rejected"];x=np.arange(3);w=.34
        for j,layer in enumerate(["prot","met"]):
            vals=[]
            for subset in order:
                q=st[(st.layer==layer)&(st.subset==subset)]
                vals.append(float(q.AUC_IPCW.iloc[0]) if len(q) else np.nan)
            axs[0].bar(x+(j-.5)*w,vals,width=w,label="Proteomics" if layer=="prot" else "Metabolomics")
        axs[0].set_xticks(x,["High predictability","All","Low predictability"])
        axs[0].set_ylabel("10-year IPCW AUC");axs[0].set_title("A  Predictability strata")
        axs[0].legend(frameon=False)
    for layer,g in c.groupby("layer"):
        g=g.sort_values("coverage_target");label="Proteomics" if layer=="prot" else "Metabolomics"
        axs[1].plot(100*g.coverage_target,g.delta_Brier,marker="o",label=label)
        axs[2].plot(100*g.coverage_target,g.delta_LogLoss,marker="o",label=label)
    for ax,title,ylabel in [(axs[1],"B  Error beyond risk composition","Δ IPCW Brier vs risk-matched"),
                            (axs[2],"C  Log-loss beyond risk composition","Δ IPCW log-loss vs risk-matched")]:
        ax.axhline(0,ls="--",linewidth=1);ax.set_xlabel("Population coverage (%)");ax.set_ylabel(ylabel);ax.set_title(title);ax.legend(frameon=False)
    fig.suptitle("Figure 3. Frozen individual predictability: discrimination and error",fontweight="bold")
    savefig(fig,Path(out)/"Fig3_predictability_coverage")
    return c,s


def _contrast(root,layer):
    d=read_table(Path(root)/layer/"c1_abm"/"abm_twin"/"test_report"/"contrasts.tsv.gz")
    if d.empty:return d
    keep=["twin_library_cleaning","twin_integrated_library_cleaning","all_cohort_training_effect","pruning_specificity"]
    z=d[(d.metric=="AUC_IPCW")&d.contrast.isin(keep)].copy();z["layer"]=layer
    return z

def fig4_cleaning(root,out):
    z=pd.concat([_contrast(root,x) for x in ["prot","met"]],ignore_index=True)
    if z.empty:return z
    order=["twin_library_cleaning","twin_integrated_library_cleaning","all_cohort_training_effect","pruning_specificity"]
    labels={"twin_library_cleaning":"Molecular-only twin: clean − full",
            "twin_integrated_library_cleaning":"Clinical+omics twin: clean − full",
            "all_cohort_training_effect":"Global elastic-net: clean − full",
            "pruning_specificity":"Global elastic-net: profile-pruned − random-pruned"}
    fig,ax=plt.subplots(figsize=(9,6))
    y=0;yt=[];yl=[]
    for layer in ["prot","met"]:
        for con in order:
            q=z[(z.layer==layer)&(z.contrast==con)]
            if not len(q):continue
            r=q.iloc[0]
            lo,hi,pt=float(r.lower),float(r.upper),float(r.delta)
            ax.hlines(y,lo,hi,linewidth=1.4);ax.plot(pt,y,"o")
            yt.append(y);yl.append(("Prot: " if layer=="prot" else "Met: ")+labels[con]);y+=1
        y+=.6
    ax.axvline(0,ls="--",linewidth=1);ax.set_yticks(yt,yl);ax.set_xlabel("Paired ΔIPCW AUC (95% group-bootstrap CI)")
    ax.set_title("Figure 4. Reference-agent cleaning: local and global model contrasts")
    savefig(fig,Path(out)/"Fig4_reference_agent_cleaning")
    return z

def supplement_selectors(root,out):
    rows=[]
    for layer in ["prot","met"]:
        d=read_table(Path(root)/layer/"c1_abm"/"baseline_comparison"/"lockpaper_selector_inference.tsv.gz")
        if len(d): d=d.copy();d["layer"]=layer;rows.append(d)
    z=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
    if z.empty:return z
    piv=z.pivot(index="route",columns="layer",values="AURC_delta_Brier")
    fig,ax=plt.subplots(figsize=(8,7));im=ax.imshow(piv.fillna(0).to_numpy(),aspect="auto")
    ax.set_yticks(range(len(piv)),piv.index);ax.set_xticks(range(len(piv.columns)),["Metabolomics" if x=="met" else "Proteomics" for x in piv.columns])
    for i in range(len(piv)):
        for j in range(len(piv.columns)):
            v=piv.iloc[i,j];ax.text(j,i,"NA" if not np.isfinite(v) else f"{v:.4f}",ha="center",va="center",fontsize=8)
    ax.set_title("Supplementary Figure. Selector-specific error reduction vs risk-matched controls")
    fig.colorbar(im,ax=ax,label="ΔBrier AURC (negative is better)")
    savefig(fig,Path(out)/"FigS_selector_multiroute")
    return z

def manuscript_grouping_audit(root):
    rows = []
    for layer in ("prot", "met"):
        path = Path(root) / layer / "c1_abm/abm_twin/host_contract.json"
        if path.is_file():
            # Describe the actual frozen fit, never a newly created global roster.
            audit = grouping_audit(json.loads(path.read_text()))
            rows.append(dict(layer=layer, **audit))
    return pd.DataFrame(rows)


def tables(out,bench,curve,summ,clean,selectors,grouping=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    with pd.ExcelWriter(out/"C1_LOCKPAPER_main_tables.xlsx",engine="openpyxl") as w:
        bench.to_excel(w,sheet_name="Table1_population_prediction",index=False)
        summ.to_excel(w,sheet_name="Table2_predictability",index=False)
        clean.to_excel(w,sheet_name="Table3_reference_cleaning",index=False)
        curve.to_excel(w,sheet_name="Fig3_source",index=False)
        selectors.to_excel(w,sheet_name="Supplement_selectors",index=False)
        if grouping is not None:
            grouping.to_excel(w,sheet_name="Supplement_kinship_audit",index=False)
    bench.to_csv(out/"Table1_population_prediction.tsv",sep="\t",index=False)
    summ.to_csv(out/"Table2_predictability.tsv",sep="\t",index=False)
    clean.to_csv(out/"Table3_reference_cleaning.tsv",sep="\t",index=False)

def main(argv=None):
    p=argparse.ArgumentParser();p.add_argument("--analysis-root",required=True);p.add_argument("--Y",default="cvd_cad");p.add_argument("--out")
    a=p.parse_args(argv);root=Path(a.analysis_root)/a.Y;out=Path(a.out) if a.out else root/"paper_c1"
    out.mkdir(parents=True,exist_ok=True)
    fig1_conceptual(out)
    bench=fig2_population(root,out);curve,summ=fig3_predictability(root,out);clean=fig4_cleaning(root,out);selectors=supplement_selectors(root,out)
    grouping = manuscript_grouping_audit(root)
    tables(out,bench,curve,summ,clean,selectors,grouping)
    grouped = (len(grouping) == 2 and {"grouping_source", "family_disjoint_outer_split", "family_disjoint_ABM_roles"} <= set(grouping)
               and grouping.grouping_source.eq(GROUPING_SOURCE).all()
               and grouping.family_disjoint_outer_split.eq(True).all()
               and grouping.family_disjoint_ABM_roles.eq(True).all())
    methods = (
        "Genetically related participants were grouped using the UK Biobank KING relatedness data. "
        "All participants belonging to the same kinship-connected component were assigned to the same "
        "outer data split, preventing close relatives from crossing the development-test boundary. "
        "Components were constructed from all released related pairs before intersection with the analytical cohort. "
        "Participants absent from the relatedness graph received singleton groups. No participants were removed "
        "because of relatedness. The same groups were used for ABM development-role assignment and "
        "kinship-component block bootstrap inference. These groups represent genetic connectivity, not known pedigrees."
        if grouped else "Kinship-aware Methods language is pending verified kinship-aware fits for both omic layers."
    )
    (out/"Methods_kinship.txt").write_text(methods + "\n")
    required=[f"{name}.{ext}" for name in ("Fig1_conceptual_framework","Fig2_population_average_prediction",
        "Fig3_predictability_coverage","Fig4_reference_agent_cleaning","FigS_selector_multiroute") for ext in ("pdf","png")]
    missing=[name for name in required if not (out/name).is_file()]
    complete=not missing and all(not d.empty and set(d.layer)=={"prot","met"} for d in (bench,curve,summ,clean,selectors))
    (out/"LOCKPAPER_REPORT.json").write_text(json.dumps(dict(release=RELEASE,status="complete" if complete else "partial",
        missing_outputs=missing,source="frozen aggregate outputs only; no fitting or method selection"),indent=2)+"\n")
    (out/"LOCKPAPER_RELEASE.txt").write_text(RELEASE+"\n")
    print(out)
if __name__=="__main__":main()
