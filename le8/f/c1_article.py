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
REPORT_VERSION="2026-10-21.same-people-reporting-v1"
LAYERS=("prot", "met")
LAYER_COLORS={"prot":"#1f77b4", "met":"#ff7f0e"}
CLEANING_CONTRASTS=("twin_library_cleaning", "twin_integrated_library_cleaning",
                    "training_effect_same_subset", "pruning_specificity")

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

def unavailable_figure(out, name, title, detail):
    """Replace a previous figure with an explicit unavailable state."""
    fig,ax=plt.subplots(figsize=(10,4));ax.axis("off")
    ax.set_title(title)
    ax.text(.5,.55,"Results not available",ha="center",fontsize=15)
    ax.text(.5,.35,detail,ha="center",va="center",wrap=True,fontsize=10)
    savefig(fig,Path(out)/name)

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

def _benchmark(perf,layer,subset="all_common",trait=None):
    """Extract conventional predictions evaluated on verified identical people.

    Native C1 holdout performance is deliberately not a fallback: sharing the
    outer split does not imply sharing the same complete-case participants.
    """
    if perf.empty:return pd.DataFrame()
    required={"subset","model","N","comparison_N","comparison_eid_hash",
              "same_people","status","Harrell_C","AUC_IPCW_10y"}
    if not required<=set(perf):
        raise ValueError("Same-person benchmark metadata missing: "+", ".join(sorted(required-set(perf))))
    d=perf[perf.subset.eq(subset)].copy()
    if "layer" in d:d=d[d.layer.eq(layer)]
    if trait is not None and "trait" in d:d=d[d.trait.eq(trait)]
    preferred="Pradeep-style / glmnet" if layer=="prot" else "All-metabolite / glmnet"
    models=["Model 0","Biomarkers","Combined"]
    names={"C1_"+m+"::"+preferred:m for m in models}
    d["model"]=d.model.astype(str).str.strip()
    q=d[d.model.isin(names)].copy()
    if len(q)!=3 or q.model.nunique()!=3:
        raise ValueError(f"Expected three unique benchmarks for {layer}/{subset}/{preferred}")
    if not q.same_people.map(lambda x:str(x).strip().lower() in ("true","1")).all():
        raise ValueError("Benchmark does not certify identical participants")
    if not q.status.isin(["evaluated","completed"]).all():
        raise ValueError("Benchmark evaluation is incomplete")
    n=pd.to_numeric(q.N,errors="coerce");comparison_n=pd.to_numeric(q.comparison_N,errors="coerce")
    hashes=q.comparison_eid_hash.fillna("").astype(str).str.strip()
    if (not np.isfinite(n).all() or (n<=0).any() or (n%1!=0).any()
            or n.nunique()!=1 or not n.eq(comparison_n).all()
            or hashes.nunique()!=1 or not hashes.str.fullmatch(r"[0-9a-fA-F]{64}").all()):
        raise ValueError("Benchmark participant counts or participant hashes differ")
    if "outer_roster_sha256" in q and q.outer_roster_sha256.nunique(dropna=False)!=1:
        raise ValueError("Benchmark outer rosters differ")
    rows=[]
    for _,r in q.iterrows():
        for metric,column in [("IPCW AUC","AUC_IPCW_10y"),("Harrell C","Harrell_C")]:
            value=float(r[column])
            if not np.isfinite(value) or not 0<=value<=1:
                raise ValueError(f"Invalid {column} for {layer}/{r.model}")
            row=r.to_dict()
            row.update(model=names[r.model],source_model=r.model,panel=preferred,
                       metric=metric,value=value,layer=layer,
                       stage="c1_correlate",family="C1 same-person holdout")
            rows.append(row)
    return pd.DataFrame(rows)

def _load_benchmark(root,layer,trait=None):
    candidates=[(Path(root)/layer/"c1_abm_comparison/same_subset_primary_performance.tsv.gz","all_common"),
                (Path(root)/layer/"c1_abm/baseline_comparison/same_subset_primary_performance.tsv.gz","all")]
    for path,subset in candidates:
        try:
            perf=read_table(path)
            if perf.empty:continue
            result=_benchmark(perf,layer,subset,trait)
            result["source_table"]=str(path)
            return result,"Verified identical participants within this omic layer"
        except (ValueError,KeyError) as exc:
            # A contradictory preferred table must not be hidden by a fallback.
            return pd.DataFrame(),f"Invalid same-person table {path.name}: {exc}"
    return pd.DataFrame(),"Completed same-person holdout results unavailable; training screening AUC is not validation"

def fig2_population(root,out,trait=None):
    rows=[];availability={}
    for layer in LAYERS:
        b,detail=_load_benchmark(root,layer,trait)
        availability[layer]=detail
        if len(b):rows.append(b)
    d=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
    d.attrs["availability"]=availability
    if d.empty:
        unavailable_figure(out,"Fig2_population_average_prediction","Figure 2. Population prediction",
                           "Completed validation on identical C1/ABM participants is unavailable.\nTraining screening AUCs are not validation estimates.")
        return d
    fig,axs=plt.subplots(1,2,figsize=(12,5.5))
    models=["Model 0","Biomarkers","Combined"];x=np.arange(len(LAYERS));w=.23
    labels=[]
    for layer,label in zip(LAYERS,["Proteomics","Metabolomics"]):
        q=d[d.layer.eq(layer)]
        labels.append(label+(f"\nN = {int(q.N.iloc[0]):,}" if len(q) else "\nNot available"))
    for ax,metric,ylabel in zip(axs,["Harrell C","IPCW AUC"],["Harrell C-index","10-year IPCW AUC"]):
        values=d[d.metric.eq(metric)]
        for j,m in enumerate(models):
            vals=[]
            for layer in LAYERS:
                q=values[(values.layer==layer)&(values.model==m)]
                vals.append(float(q.value.iloc[0]) if len(q) else np.nan)
            ax.bar(x+(j-1)*w,vals,width=w,label={"Model 0":"Clinical","Biomarkers":"Omics","Combined":"Clinical + omics"}[m])
        ax.set_xticks(x,labels);ax.set_ylabel(ylabel)
        finite=pd.to_numeric(values.value,errors="coerce").dropna()
        if len(finite):ax.set_ylim(min(.5,finite.min()-.03),min(1,finite.max()+.03))
    axs[0].legend(frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(1.08,-.28))
    fig.suptitle("Figure 2. Population prediction on common C1/ABM participants")
    fig.text(.5,.015,"Identical participants across models within each omic layer; omic layers have different cohorts.",ha="center",fontsize=9)
    fig.subplots_adjust(bottom=.24,wspace=.3)
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
    if c.empty:
        unavailable_figure(out,"Fig3_predictability_coverage","Figure 3. Frozen predictability",
                           "Frozen-selector holdout coverage and error estimates are unavailable.")
        return c,s
    fig,axs=plt.subplots(1,3,figsize=(16,6.7))
    if len(st):
        order=["selected","all","rejected"];x=np.arange(3);w=.34
        for j,layer in enumerate(["prot","met"]):
            vals=[]
            for subset in order:
                q=st[(st.layer==layer)&(st.subset==subset)]
                vals.append(float(q.AUC_IPCW.iloc[0]) if len(q) else np.nan)
            axs[0].bar(x+(j-.5)*w,vals,width=w,color=LAYER_COLORS[layer],label="Proteomics" if layer=="prot" else "Metabolomics")
        axs[0].set_xticks(x,["Selected","All","Rejected"])
        axs[0].set_ylabel("10-year IPCW AUC");axs[0].set_title("A  Frozen selector strata")
        axs[0].legend(frameon=False)
    for layer,g in c.groupby("layer"):
        g=g.sort_values("coverage_target");label="Proteomics" if layer=="prot" else "Metabolomics"
        axs[1].plot(100*g.coverage_target,g.delta_Brier,marker="o",color=LAYER_COLORS[layer],label=label)
        axs[2].plot(100*g.coverage_target,g.delta_LogLoss,marker="o",color=LAYER_COLORS[layer],label=label)
    for ax,title,ylabel in [(axs[1],"B  Error beyond risk composition","Δ IPCW Brier vs risk-matched"),
                            (axs[2],"C  Log-loss beyond risk composition","Δ IPCW log-loss vs risk-matched")]:
        ax.axhline(0,ls="--",linewidth=1);ax.set_xlabel("Target coverage (%)");ax.set_ylabel(ylabel);ax.set_title(title);ax.legend(frameon=False)
    fig.suptitle("Figure 3. Frozen individual predictability: discrimination and error",fontweight="bold")
    targets=sorted(pd.to_numeric(c.coverage_target,errors="coerce").dropna().unique())
    coverage_rows=[];coverage_labels=[]
    for layer,label in zip(LAYERS,["Proteomics","Metabolomics"]):
        g=c[c.layer.eq(layer)]
        if g.empty:continue
        row=[]
        for target in targets:
            q=g[np.isclose(pd.to_numeric(g.coverage_target,errors="coerce"),target)]
            if len(q)==1 and {"coverage","selected_N"}<=set(q):
                actual=float(q.coverage.iloc[0]);selected_n=float(q.selected_N.iloc[0])
                row.append(f"{100*actual:.1f}%\nN={int(selected_n):,}" if np.isfinite(actual) and np.isfinite(selected_n) else "Not available")
            else:row.append("Not available")
        coverage_rows.append(row);coverage_labels.append(label)
    ax_table=fig.add_axes([.40,.035,.57,.15]);ax_table.axis("off")
    if coverage_rows:
        table=ax_table.table(cellText=coverage_rows,rowLabels=coverage_labels,
                             colLabels=[f"Target {100*x:g}%" for x in targets],loc="center",cellLoc="center")
        table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.65)
    fig.text(.40,.225,"Actual population coverage and selected N at each target",fontsize=9)
    fig.subplots_adjust(bottom=.35,top=.86,wspace=.35)
    savefig(fig,Path(out)/"Fig3_predictability_coverage")
    return c,s


def _contrast(root,layer):
    d=read_table(Path(root)/layer/"c1_abm"/"abm_twin"/"test_report"/"contrasts.tsv.gz")
    if d.empty:return d
    required={"metric","contrast","subset","reference_subset","model","reference","delta","lower","upper"}
    if not required<=set(d):return pd.DataFrame()
    # Both local and global changes are conditional on the same selected people.
    z=d[(d.metric=="AUC_IPCW")&d.contrast.isin(CLEANING_CONTRASTS)
        &d.subset.eq("selected")&d.reference_subset.eq("selected")].copy()
    metrics=read_table(Path(root)/layer/"c1_abm/abm_twin/test_report/metrics.tsv.gz")
    if metrics.empty or not {"model","subset","N","AUC_IPCW"}<=set(metrics):return pd.DataFrame()
    rows=[]
    for _,r in z.iterrows():
        a=metrics[(metrics.model==r.model)&metrics.subset.eq("selected")]
        b=metrics[(metrics.model==r.reference)&metrics.subset.eq("selected")]
        if len(a)!=1 or len(b)!=1:continue
        a,b=a.iloc[0],b.iloc[0]
        if not np.isfinite(float(a.N)) or float(a.N)<=0 or a.N!=b.N:continue
        if any(col in metrics and a[col]!=b[col] for col in ("known_N","cases","controls","requested_N")):continue
        if not np.isclose(float(r.delta),float(a.AUC_IPCW)-float(b.AUC_IPCW),rtol=1e-7,atol=1e-10):continue
        row=r.to_dict();row.update(layer=layer,comparison_N=int(a.N),comparison_scope="selected ABM test participants")
        if "Harrell_C" in metrics:
            row["delta_Harrell_C"]=float(a.Harrell_C)-float(b.Harrell_C)
        # No Harrell interval is inferred from a stored Uno C interval.
        rows.append(row)
    z=pd.DataFrame(rows)
    if len(z) and (z.contrast.duplicated().any() or z.comparison_N.nunique()!=1):
        return pd.DataFrame()
    return z

def fig4_cleaning(root,out):
    z=pd.concat([_contrast(root,x) for x in ["prot","met"]],ignore_index=True)
    if z.empty:
        unavailable_figure(out,"Fig4_reference_agent_cleaning","Figure 4. Reference-agent cleaning",
                           "Paired local and global contrasts on the same selected participants are unavailable.")
        return z
    order=CLEANING_CONTRASTS
    labels={"twin_library_cleaning":"Molecular-only twin: clean − full",
            "twin_integrated_library_cleaning":"Clinical+omics twin: clean − full",
            "training_effect_same_subset":"Global elastic-net: clean − full",
            "pruning_specificity":"Global elastic-net: profile-pruned − random-pruned"}
    fig,ax=plt.subplots(figsize=(11,6))
    y=0;yt=[];yl=[]
    for layer in ["prot","met"]:
        for con in order:
            q=z[(z.layer==layer)&(z.contrast==con)]
            if not len(q):continue
            r=q.iloc[0]
            lo,hi,pt=float(r.lower),float(r.upper),float(r.delta)
            ax.hlines(y,lo,hi,linewidth=1.4,color=LAYER_COLORS[layer]);ax.plot(pt,y,"o",color=LAYER_COLORS[layer])
            yt.append(y);yl.append(("Prot: " if layer=="prot" else "Met: ")+labels[con]+f" (N={int(r.comparison_N):,})");y+=1
        y+=.6
    ax.axvline(0,ls="--",linewidth=1);ax.set_yticks(yt,yl);ax.set_xlabel("Paired ΔIPCW AUC (95% group-bootstrap CI)")
    ax.set_title("Figure 4. Reference-agent cleaning on the same selected participants")
    savefig(fig,Path(out)/"Fig4_reference_agent_cleaning")
    return z

def supplement_selectors(root,out):
    rows=[]
    for layer in ["prot","met"]:
        d=read_table(Path(root)/layer/"c1_abm"/"baseline_comparison"/"lockpaper_selector_inference.tsv.gz")
        if len(d): d=d.copy();d["layer"]=layer;rows.append(d)
    z=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
    if z.empty:
        unavailable_figure(out,"FigS_selector_multiroute","Supplementary Figure. Selector comparison",
                           "Held-out selector inference is unavailable.")
        return z
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


def reporting_status(root,bench,curve,summ,clean,selectors):
    """Describe evidence availability independently of intermediate filenames."""
    rows=[]
    modules=[("population_prediction",bench,6),("predictability_curve",curve,6),
             ("predictability_inference",summ,1),("reference_cleaning",clean,4),
             ("selector_sensitivity",selectors,1)]
    for layer in LAYERS:
        for module,frame,minimum in modules:
            q=frame[frame.layer.eq(layer)] if "layer" in frame else pd.DataFrame()
            status="complete" if len(q)>=minimum else ("incomplete" if len(q) else "not_available")
            detail="Frozen validation estimates available" if status=="complete" else "Completed validation estimates unavailable"
            if module=="population_prediction":
                detail=bench.attrs.get("availability",{}).get(layer,detail)
                if status!="complete" and (Path(root)/layer/"c1_correlate/prediction/prediction.results.xlsx").is_file():
                    detail+="; an intermediate prediction workbook does not establish completed validation"
            if module=="predictability_curve" and len(q):
                if not {"coverage","selected_N"}<=set(q) or q[["coverage","selected_N"]].isna().any().any():
                    status="incomplete";detail="Actual coverage or selected participant counts unavailable"
            rows.append(dict(layer=layer,module=module,status=status,rows=len(q),detail=detail))
    return pd.DataFrame(rows)

def tables(out,bench,curve,summ,clean,selectors,grouping=None,reporting=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    with pd.ExcelWriter(out/"C1_LOCKPAPER_main_tables.xlsx",engine="openpyxl") as w:
        if reporting is not None:
            reporting.to_excel(w,sheet_name="Table0_reporting_status",index=False)
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
    if reporting is not None:
        reporting.to_csv(out/"Table0_reporting_status.tsv",sep="\t",index=False)

def main(argv=None):
    p=argparse.ArgumentParser();p.add_argument("--analysis-root",required=True);p.add_argument("--Y",default="cvd_cad");p.add_argument("--out")
    a=p.parse_args(argv);root=Path(a.analysis_root)/a.Y;out=Path(a.out) if a.out else root/"paper_c1"
    out.mkdir(parents=True,exist_ok=True)
    fig1_conceptual(out)
    bench=fig2_population(root,out,trait=a.Y);curve,summ=fig3_predictability(root,out);clean=fig4_cleaning(root,out);selectors=supplement_selectors(root,out)
    grouping = manuscript_grouping_audit(root)
    reporting=reporting_status(root,bench,curve,summ,clean,selectors)
    tables(out,bench,curve,summ,clean,selectors,grouping,reporting)
    grouped = (len(grouping) == 2 and {"layer", "grouping_unit", "grouping_source", "family_disjoint_outer_split", "family_disjoint_ABM_roles"} <= set(grouping)
               and set(grouping.layer) == {"prot", "met"}
               and grouping.grouping_unit.eq("family").all()
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
    complete=grouped and not missing and reporting.status.eq("complete").all()
    (out/"LOCKPAPER_REPORT.json").write_text(json.dumps(dict(release=RELEASE,reporting_version=REPORT_VERSION,
        status="complete" if complete else "incomplete",components=reporting.to_dict("records"),
        kinship_validation="complete" if grouped else "pending",
        kinship_grouping_source=GROUPING_SOURCE if grouped else "unverified",
        missing_outputs=missing,source="frozen aggregate outputs only; no fitting or method selection"),indent=2)+"\n")
    (out/"LOCKPAPER_RELEASE.txt").write_text(RELEASE+"\n")
    print(out)
if __name__=="__main__":main()
