"""Aggregate tables and publication-oriented figures; never changes a fitted model."""
from pathlib import Path
import html
import json
import shutil
import importlib
import subprocess
import numpy as np
import pandas as pd
from core2 import atomic_table,read_table,identified

def architecture(cfg,root,trait):
    root=Path(root);folder=root/f"tmp/{trait}/regions"
    roster=identified(read_table(root/"tmp/cohort/participants.tsv.gz"));samples=read_table(folder/"ids.tsv")
    features=read_table(folder/"features.tsv");x=np.load(folder/"scores.npy",mmap_mode="r")
    ids=set(roster.loc[roster.role.eq("train"),"id"]);rows=np.flatnonzero(samples.id.isin(ids))
    mean=np.zeros(x.shape[1]);m2=mean.copy();seen=0
    for start in range(0,len(rows),256):
        a=np.asarray(x[rows[start:start+256]],float);bn=len(a);bm=a.mean(0);delta=bm-mean
        m2+=((a-bm)**2).sum(0)+delta**2*seen*bn/(seen+bn);mean+=delta*bn/(seen+bn);seen+=bn
    variance=m2/seen;key=next(c for c in ("map","population","pop") if c in features);records=[]
    for label in features[key].unique():
        ix=features[key].eq(label).to_numpy();v=variance[ix];v=v[v>1e-14]
        if not len(v):continue
        v=np.sort(v)[::-1];total=v.sum()
        records.append({"trait":trait,"map":label,"n_variable_blocks":len(v),"effective_blocks":total**2/np.square(v).sum(),
            "effective_fraction":total**2/np.square(v).sum()/len(v),"top_1pct_variance_fraction":v[:max(1,int(np.ceil(.01*len(v))))].sum()/total,
            "top_10pct_variance_fraction":v[:max(1,int(np.ceil(.10*len(v))))].sum()/total,"role":"TRAIN_X_only",
            "interpretation":"concentration of computed regional-score variance; not LD-independent explained trait variance or heritability"})
    atomic_table(pd.DataFrame(records),root/f"qc/{trait}.regional_architecture.tsv")

def figures(performance,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":220})
    traits=list(performance.trait.unique());colors={"csx.population":"#27649a","disco":"#b25d28","abm":"#287d62","prscs.best":"#8452a0","abm.deep":"#a3346c"}
    n_methods=performance.loc[performance.selection.eq("all")].groupby("trait").size().max()
    fig,axes=plt.subplots(1,len(traits),figsize=(max(5,4.6*len(traits)),max(4,.28*n_methods+1.2)),squeeze=False)
    for ax,trait in zip(axes[0],traits):
        d=performance[(performance.trait==trait)&performance.selection.eq("all")].copy()
        key="auroc" if trait=="t2dm" else "r2"
        d=d.sort_values(key);pos=np.arange(len(d))
        ax.barh(pos,d[key],color=[colors.get(m,"#8796a5") for m in d.method]);ax.set_yticks(pos,labels=d.method)
        ax.set_xlabel("AUROC" if trait=="t2dm" else "Predictive R²");ax.set_title(trait);ax.axvline(0,color="#555",lw=.6)
    fig.tight_layout();fig.savefig(out/"overall_performance.png");fig.savefig(out/"overall_performance.pdf");plt.close(fig)
    fig,axes=plt.subplots(1,len(traits),figsize=(max(5,4.2*len(traits)),3.7),squeeze=False)
    for ax,trait in zip(axes[0],traits):
        d=performance[performance.trait.eq(trait)&performance.selection.str.startswith("selected_support_")]
        for method,color in colors.items():
            r=d[d.method.eq(method)].sort_values("actual_coverage")
            if r.empty:continue
            ax.plot(r.actual_coverage,r.selected_minus_all_loss,"o-",label=method,color=color,ms=4)
            ax.fill_between(r.actual_coverage,r.selection_low,r.selection_high,color=color,alpha=.13)
        ax.axhline(0,color="#666",lw=.7);ax.set_title(trait);ax.set_xlabel("Actual selected fraction of test cohort")
        ax.set_ylabel("Selected − full-test Brier" if trait=="t2dm" else "Selected − full-test MSE")
        ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(out/"support_risk_coverage.png");fig.savefig(out/"support_risk_coverage.pdf");plt.close(fig)

def final_report(cfg,root):
    root=Path(root);out=root/"final";out.mkdir(parents=True,exist_ok=True)
    data=[pd.read_csv(out/f"{t}.performance.tsv",sep="\t") for t in cfg["traits"]]
    ancestry=[pd.read_csv(out/f"{t}.ancestry.tsv",sep="\t") for t in cfg["traits"]]
    all_perf=pd.concat(data,ignore_index=True);all_ancestry=pd.concat(ancestry,ignore_index=True)
    atomic_table(all_perf,out/"all.performance.tsv");atomic_table(all_ancestry,out/"all.ancestry.tsv")
    figures(all_perf,out)
    whole=all_perf[all_perf.selection.eq("all")]
    columns=[c for c in ["trait","method","n_evaluable","r2","auroc","auprc","primary_loss","calibration_intercept","calibration_slope","method_minus_primary_loss","paired_low","paired_high"] if c in whole]
    selected=all_perf[all_perf.selection.eq("selected_support_0p20")]
    selected_columns=[c for c in ["trait","method","actual_coverage","selected_minus_all_loss","selection_low","selection_high","risk_standardized_difference","risk_difference_low","risk_difference_high","overlap_selected_fraction","overlap_rejected_fraction","overlap_cells","risk_valid_bootstrap_replicates","risk_status"] if c in selected]
    models=[]
    for trait in cfg["traits"]:
        for label,path in [("CPU ABM",root/f"tmp/{trait}/abm/model.json"),("Disco",root/f"tmp/{trait}/disco/model.json"),
                           ("Regional glmnet / LightGBM",root/f"tmp/{trait}/block_models/model.json"),("Block deep",root/f"tmp/{trait}/deep/model.json")]:
            if path.is_file():models.append(f"<h3>{html.escape(trait)} / {label}</h3><pre>{html.escape(path.read_text())}</pre>")
    links=" ".join(f'<a href="{html.escape(t)}.report.html">{html.escape(t)} detail</a>'+
                   (f'<a href="{html.escape(t)}.twins.html">{html.escape(t)} twin support</a>' if (out/f"{t}.twins.html").is_file() else "")
                   for t in cfg["traits"])
    status=html.escape(cfg["cohort"].get("validation_status","exploratory_previously_inspected_test"))
    text=f'''<!doctype html><html><head><meta charset="utf-8"><title>GRID FINAL analysis</title>
<style>body{{font:15px system-ui;margin:32px auto;max-width:1400px;padding:0 24px;color:#18313b}}h1,h2{{color:#174b65}}p{{max-width:1000px;line-height:1.65}}table{{border-collapse:collapse;font-size:12px}}td,th{{padding:7px;border:1px solid #d7e0e5}}th{{background:#eef5f8}}img{{max-width:100%}}pre{{white-space:pre-wrap;font-size:11px;background:#f4f7f8;padding:16px}}a{{margin-right:18px}}</style></head><body>
<h1>GRID: consensus-block prediction and local support</h1>
<p>Validation status: <b>{status}</b>. Input genotypes: GRCh38 imputed PGEN. Every comparator uses the same held-out participants and covariates. Training and tuning budgets are shared. Test phenotypes enter only evaluation. Counts of matched references refer to a frozen genetic-score representation, not identical genomes.</p>
<p>{links}<a href="all.performance.tsv">All metrics TSV</a><a href="all.ancestry.tsv">Ancestry metrics TSV</a></p>
<h2>Full test cohort</h2><img src="overall_performance.png" alt="Full test performance">
{whole[columns].to_html(index=False,float_format=lambda x:f"{x:.6g}",border=0)}
<h2>Does high support transfer across predictors?</h2><p>The same calibration-frozen CPU support mask is applied to CSx, Disco, PRS-CS, ABM and every available regional comparator. Deep-embedding selectors are separately labelled exploratory. Negative selected-minus-all loss means lower mean error in the selected population. Composition-adjusted selected-versus-rejected contrasts and their overlap diagnostics answer a different question and must be read together.</p>
<img src="support_risk_coverage.png" alt="Support risk coverage curves">
{selected[selected_columns].to_html(index=False,float_format=lambda x:f"{x:.6g}",border=0)}
<p>Bootstrap intervals are pointwise, conditional on fitted models, and resample declared family groups. They do not include GWAS or model-training uncertainty. An OTH ancestry label is not a complete definition of admixture. Large family ESS is a concentration diagnostic, not a guarantee of statistical independence or individual accuracy.</p>
<h2>Frozen model choices</h2>{''.join(models)}</body></html>'''
    (out/"index.html").write_text(text)
    return out/"index.html"

def check_inputs(cfg):
    from genomics import pfile
    inputs=cfg["inputs"];errors=[]
    for c in cfg["chromosomes"]:
        try:pfile(Path(inputs["genotype_dir"])/f"chr{c}")
        except (FileNotFoundError,ValueError) as exc:errors.append(str(exc))
    source=Path(cfg["analysis_root"])/"tmp/cohort/participants.tsv.gz"
    if not Path(source).is_file():errors.append(f"Missingphenotype/cohort:{source}")
    try:
        from consensus3 import read_consensus
        _,catalog=read_consensus(inputs["consensus_dir"],expected_chromosomes=cfg["chromosomes"],build=inputs["consensus_build"])
    except (ValueError,FileNotFoundError) as exc:errors.append(str(exc));catalog=None
    for key in ("plink","rscript"):
        executable=cfg["runtime"][key]
        if not shutil.which(executable):errors.append(f"Missingexecutable:{executable}")
    for module in ("numpy", "pandas", "scipy", "sklearn", "h5py", "yaml", "threadpoolctl", "pgenlib", "matplotlib", "zstandard", "joblib") + (("lightgbm",) if cfg["block_models"].get("enabled") else ()):
        try: importlib.import_module(module)
        except (ImportError, OSError) as exc: errors.append(f"Cannot import {module}: {exc}")
    if shutil.which(cfg["runtime"]["rscript"]):
        packages=["jsonlite"]
        if cfg["block_models"].get("enabled") and cfg["block_models"].get("glmnet"): packages.append("glmnet")
        expression='p <- c('+','.join(json.dumps(p) for p in packages)+'); missing <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop(paste(missing, collapse=", "))'
        result=subprocess.run([cfg["runtime"]["rscript"],"--vanilla","-e",expression],capture_output=True,text=True)
        if result.returncode: errors.append("R dependency check: "+result.stderr.strip())
    print(json.dumps({"genotype_dir":inputs["genotype_dir"],"reference_build":"GRCh37","target_build":"GRCh38",
        "consensus_blocks":catalog["n_blocks"] if catalog else None,"LD_rebuild_requested":False,
        "errors":errors,"status":"PASS" if not errors else "INPUTS_NEEDED"},indent=2))
    if errors:raise ValueError("Inputcheckfailed; fixlistedpaths/toolsbeforefullrun")

def twin_report(cfg,root,trait):
    """Aggregate frozen neighbourhoods; never re-search or inspect query outcomes."""
    root=Path(root);work=root/f"tmp/{trait}/abm"
    d=read_table(work/"abm.predictions.tsv.gz")
    expected={"role","ancestry","n_support","support_group_ess","supported"}
    if not expected<=set(d):raise ValueError("ABM output lacks the frozen twin support fields")
    from evaluate2 import boolean_mask
    d["supported"]=boolean_mask(d.supported)
    for name in ("n_support", "support_group_ess"):
        d[name] = pd.to_numeric(d[name], errors="raise")
    rows=[]
    for (role,ancestry),g in d.groupby(["role","ancestry"],dropna=False,sort=True):
        row=dict(trait=trait,role=str(role),ancestry=str(ancestry),n=len(g),
            supported_fraction=float(g.supported.mean()),support_median=float(g.n_support.median()),
            support_p90=float(g.n_support.quantile(.9)),family_ess_median=float(g.support_group_ess.median()))
        for name in ("nearest_distance","nearest_rms_ratio","kth_rms_ratio"):
            if name in g:row[name+"_median"]=float(pd.to_numeric(g[name]).replace([np.inf,-np.inf],np.nan).median())
        rows.append(row)
    result=pd.DataFrame(rows)
    atomic_table(result,root/f"qc/{trait}.twins.tsv")
    out=root/"final";out.mkdir(exist_ok=True,parents=True)
    text=f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>{html.escape(trait)} frozen twin support</title>
    <style>body{{font:15px system-ui;max-width:1300px;margin:36px auto;padding:0 20px}}table{{border-collapse:collapse}}td,th{{border:1px solid #ddd;padding:7px;text-align:right}}</style>
    <h1>{html.escape(trait)}: frozen block-profile support</h1>
    <p>Counts include every eligible donor within the frozen radius, excluding the query family.
    They are not limited by the prediction k or the number of displayed donor examples.
    Exact search uses the declared feature representation; these are genetic-score neighbours, not biological identical twins.</p>
    <p>Nearest/RMS and kth/RMS ratios close to 1 indicate distance concentration.
    Individual donor evidence remains in the local tmp/{html.escape(trait)}/abm/twins.tsv.gz file.
    This report reads no query outcomes and performs no matching again.</p>{result.to_html(index=False,escape=True)}</html>"""
    (out/f"{trait}.twins.html").write_text(text)
    return result

