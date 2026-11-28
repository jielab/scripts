"""Paired LDL analyses on frozen predictions. Participant-level inputs stay out of the paper output."""
import numpy as np
import pandas as pd

from blocknet_twins import POLICIES,PREFIX
from core2 import aligned,atomic_table,read_table
from evaluate2 import boolean_mask,interval,paired_cluster_contrasts,risk_standardization

ABLATIONS = (
    ("direct.blocknet","direct.linear"),
    ("direct.regional_ridge","direct.linear"),
    ("direct.global_mlp","direct.linear"),
    ("direct.blocknet","direct.regional_ridge"),
    ("direct.blocknet","direct.global_mlp"),
)


def twin_audit(cohort,truth,benchmark,abm,features,folder,reps,seed,workers,trait):
    out=folder/f"{trait}.blocknet_twins.tsv"
    if features is None or "direct.blocknet.prediction" not in benchmark:
        atomic_table(pd.DataFrame([{"status":"not available"}]),out)
        return
    roster=cohort[cohort.role.eq("test")]
    frame=aligned(features,roster.id)
    baseline=aligned(benchmark,roster.id)
    p=baseline["direct.blocknet.prediction"].to_numpy(float)
    linear=baseline["direct.linear.prediction"].to_numpy(float)
    y=aligned(truth,roster.id).y.to_numpy(float)
    raw=aligned(abm,roster.id)
    policies=("all","raw_global",*POLICIES)
    masks=np.array([np.ones(len(frame),bool),boolean_mask(raw.selected_twin_count_0p20),
                    *[boolean_mask(frame[PREFIX+n]) for n in POLICIES]])
    groups=roster.group_id.to_numpy(str);pop=roster.ancestry.to_numpy(str)
    cells=np.char.add(np.char.add(pop,":"),frame.twin_risk_bin.to_numpy(str))
    rows=[]
    for anc in ["ALL",*sorted(set(pop))]:
        idx=np.ones(len(p),bool) if anc=="ALL" else pop==anc
        sel=masks[:,idx];g=groups[idx];yy=y[idx];pp=p[idx];obs=np.isfinite(yy)
        means,selections,contrasts,_=paired_cluster_contrasts(yy,np.stack([linear[idx],pp]),sel,g,reps,seed,workers)
        losses=means[1];boot=selections[:,1,:]
        gain=means[1]-means[0]
        gain_boot=contrasts[:,1,:]
        risk=policies.index("risk_only")
        evaluable_families=[len(np.unique(g[mask&obs])) for mask in sel]
        def family_interval(draws,*arms):
            if any(evaluable_families[a]<2 for a in arms):
                return np.nan,np.nan,0,"insufficient_families"
            return interval(draws,reps)
        for j,label in enumerate(policies):
            mask=sel[j];ok=mask&obs;loss=losses[j]
            var=np.var(yy[ok]) if ok.any() else np.nan
            delta=loss-losses[0]
            lo,hi,valid,ci_status=family_interval(boot[:,j],j,0)
            vsrisk=loss-losses[risk]
            rl,rh,rv,_=family_interval(boot[:,j]-boot[:,risk],j,risk)
            count=frame.loc[idx,"twin_count"].to_numpy(float)[mask]
            rec=dict(status="available" if ok.any() else "no_evaluable_rows",ancestry=anc,policy=label,
                     n_roster=len(mask),n_selected=int(mask.sum()),n_evaluable=int(ok.sum()),
                     n_evaluable_families=evaluable_families[j],coverage=mask.mean(),
                     median_latent_twins=np.nanmedian(count) if np.isfinite(count).any() else np.nan,
                     mse=loss,mae=np.mean(np.abs(yy[ok]-pp[ok])) if ok.any() else np.nan,
                     predictive_r2=1-loss/var if var>0 else np.nan,selected_minus_all_mse=delta,
                     ci_low=lo,ci_high=hi,ci_status=ci_status,valid_draws=valid,selected_minus_risk_only_mse=vsrisk,
                     risk_only_ci_low=rl,risk_only_ci_high=rh,risk_only_valid_draws=rv)
            gl,gh,_,_=family_interval(gain_boot[:,j],j)
            dl,dh,_,_=family_interval(gain_boot[:,j]-gain_boot[:,0],j,0)
            rec.update(blocknet_minus_linear_mse=gain[j],gain_ci_low=gl,gain_ci_high=gh,
                       blocknet_minus_linear_r2=-gain[j]/var if var>0 else np.nan,
                       gain_selected_minus_all_mse=gain[j]-gain[0],gain_selection_ci_low=dl,gain_selection_ci_high=dh)
            if anc=="ALL" and j and mask.any() and (~mask).any():
                rec.update(risk_standardization(yy,pp[None,:],mask,g,cells[idx],reps,seed,workers)[0])
            rows.append(rec)
    # Does a count-aware model forecast squared errors better on the SAME individuals?
    heads=np.array([frame[name+".twin_expected_loss"].to_numpy(float) for name in POLICIES[-3:]])
    valid=np.isfinite(heads).all(0)&np.isfinite(y)
    if valid.any():
        true_loss=np.square(y[valid]-p[valid])
        m,_,b,_=paired_cluster_contrasts(true_loss,heads[:,valid],[np.ones(valid.sum(),bool)],groups[valid],reps,seed,workers)
        for i,name in enumerate(POLICIES[-3:]):
            row=next(r for r in rows if r["ancestry"]=="ALL" and r["policy"]==name)
            lo,hi,n,_=interval(b[:,i,0],reps)
            row.update(error_forecast_n=int(valid.sum()),error_forecast_mse=m[i,0],
                       error_forecast_mse_minus_risk_only=m[i,0]-m[0,0],error_forecast_ci_low=lo,error_forecast_ci_high=hi)
    atomic_table(pd.DataFrame(rows),out)


def plot(tables,path):
    import matplotlib.pyplot as plt
    fig,(left,right)=plt.subplots(1,2,figsize=(11,4.2),layout="constrained")
    names={"direct.linear":"Global ridge","direct.global_mlp":"Global MLP","direct.regional_ridge":"Regional ridge","direct.blocknet":"BlockNet"}
    overall=tables["overall"]
    values=overall[(overall.trait=="ldl")&(overall.selection=="all")].set_index("method")
    order=[k for k in names if k in values.index]
    for i,k in enumerate(order):
        left.scatter(values.at[k,"r2"],i,color="#00765E" if k=="direct.blocknet" else "#295b86")
        left.annotate(f"{values.at[k,'r2']:.4f}",(values.at[k,"r2"],i),xytext=(0,8),textcoords="offset points",ha="center")
    left.set(yticks=range(len(order)),yticklabels=[names[k] for k in order],xlabel="Full TEST predictive R²",title="LDL prediction")
    left.margins(x=.22,y=.3)
    twins=tables["blocknet_twins"]
    if "policy" in twins:
        d=twins.loc[(twins.ancestry=="ALL")&(twins.policy!="all")].reset_index(drop=True)
        for i,r in d.iterrows():
            if np.isfinite([r.ci_low,r.ci_high]).all():
                right.plot([r.ci_low,r.ci_high],[i,i],color="#397988")
            right.scatter(r.selected_minus_all_mse,i,color="#00765E")
        right.set_yticks(range(len(d)),d.policy.str.replace("_"," "))
    else:
        right.text(.5,.5,"Screening estimates unavailable",transform=right.transAxes,ha="center")
    right.axvline(0,color="gray",linewidth=1)
    right.set(xlabel="Selected − all TEST MSE (lower is better)",title="Frozen LDL screening rules")
    for ax in (left,right):ax.spines[["top","right"]].set_visible(False)
    fig.savefig(path,dpi=240);plt.close(fig)


def boundary_inputs(cfg,root):
    """Use only the two named analyses in this project, on their common upstream cohort."""
    import json
    from pathlib import Path
    from driver2 import verified_optional_predictions
    from output_layout import upstream_root,trait_path
    main=upstream_root(root)/"main"
    roots={"consensus":main,"ancestry":main/"ancestry"}
    active=next((name for name,r in roots.items() if r==Path(root)),None)
    if active is None or not cfg["block_deep"]["enabled"]:
        return {}, "standard boundary comparison not requested"
    paths={};contract=None
    for layout,r in roots.items():
        meta=trait_path(r,"ldl","regions/metadata.json")
        if not meta.is_file():
            return {},"awaiting regional scores: "+layout
        md=json.loads(meta.read_text())
        if md.get("layout")!=layout or md.get("representation")!="ancestry_channels" or md.get("block_build","GRCh38")!="GRCh38":
            return {},"requires uniform and ancestry GRCh38 four-channel inputs"
        ready=verified_optional_predictions(cfg,r,"ldl",skipped={},stages={"direct_blocknet"})
        if not ready:
            return {},"awaiting completed matched BlockNet: "+layout
        info=json.loads(trait_path(r,"ldl","direct_blocknet/model.json").read_text())
        signature=(md.get("geometry_phi"),info.get("device"),info.get("side_columns"),
                   info.get("n_build"),info.get("n_tune_model"),info.get("n_prediction_calibration"))
        if contract is not None and signature!=contract:
            return {},"different CSx phi, device, covariates or development split; comparison withheld"
        contract=signature
        paths[layout]=ready["direct_blocknet"]
    return paths,"available"


def boundary_audit(cfg,root,cohort,truth,folder):
    from manuscript_stats import family_contrasts
    paths,status=boundary_inputs(cfg,root)
    output=folder/"ldl.boundaries.tsv"
    if not paths:
        atomic_table(pd.DataFrame([dict(status=status)]),output);return
    test=cohort.loc[cohort.role.eq("test")]
    methods={};baselines={}
    for layout,path in paths.items():
        d=aligned(read_table(path,columns=["id","direct.blocknet.prediction","direct.linear.prediction"]),test.id)
        methods[layout]=d["direct.blocknet.prediction"].to_numpy(float)
        baselines[layout]=d["direct.linear.prediction"].to_numpy(float)
    if not np.allclose(baselines["consensus"],baselines["ancestry"],rtol=1e-8,atol=1e-8):
        atomic_table(pd.DataFrame([dict(status="matched global ridge predictions differ; comparison withheld")]),output)
        return
    yy=aligned(truth,test.id).y.to_numpy(float)
    table=family_contrasts(yy,methods,test.group_id.to_numpy(str),test.ancestry.to_numpy(str),
                           np.ones(len(test),bool),"ldl",cfg["evaluation"]["bootstrap"],cfg["cohort"]["seed"],
                           only_pairs=[("consensus","ancestry")])
    table=table.loc[(table.selection=="all")&table.metric.isin(["MSE","predictive_R2","relative_loss_reduction_pct"])].copy()
    table["status"]="available"
    atomic_table(table,output)

