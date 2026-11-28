"""Incremental paired evaluation on a fixed held-out roster, with aggregate-only publication."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve,roc_curve
from core2 import aligned,atomic_json,atomic_table,check_roles,read_table,stamp
from evaluate2 import FrozenDiscrimination,boolean_mask,interval,metrics
from file_lock import file_lock
from output_layout import project_path,trait_path
from methods import published

REFS=("csx.population","disco")
ARTIFACTS=("performance","contrasts","curves","screening")
HERE=Path(__file__).parent

def short(name):
    return "cs."+name[6:] if name.startswith("prscs.") else name

def discover(cfg,root,t):
    sources={}
    for label,path in (("csx",trait_path(root,t,"baseline.predictions.tsv.gz")),
                       ("disco",trait_path(root,t,"predictions.tsv.gz"))):
        if path.is_file():sources[label]=path
        elif label=="csx":raise FileNotFoundError(path)
    skipped={}
    records={}
    state={}
    def ready(mark):
        if mark in state:return state[mark]
        state[mark]=None
        doc=published(mark)
        if doc:
            for f in doc.get("inputs",{}).get("files",[]):
                parent=Path(f["path"])
                current=stamp(parent) if parent.is_file() else {"path":str(parent),"missing":True}
                if current!=f or (parent.name=="complete.json" and not ready(parent)):
                    doc=None;break
        state[mark]=doc
        return doc
    for mark in sorted((root/t).glob("*/complete.json")):
        doc=ready(mark)
        if not doc:
            skipped[mark.parent.name]="incomplete or stale dependency";continue
        if doc["trait"]!=t:raise ValueError("Method receipt belongs to a different trait: "+str(mark))
        records[doc["method"]]=mark
        if doc.get("predictions") and doc["method"]!="blocknet.twins":
            sources[doc["method"]]=mark.parent/doc["predictions"]
    from methods import PARENTS,folder
    for m in PARENTS:
        if folder(root,t,m).is_dir() and m not in records:skipped.setdefault(m,"incomplete")
    return sources,records,skipped


def predcols(path,source):
    cols=pd.read_csv(path,sep="\t",nrows=0).columns
    return {short(c[:-11]):c for c in cols if c.endswith(".prediction") and
            (source!="disco" or c.startswith("disco."))}

def bootstrap(y,p,refs,groups,binary,reps,seed):
    """One family resample for each candidate and its references, including ties in ranking."""
    _,ix=np.unique(groups,return_inverse=True)
    ng=int(ix.max())+1
    counts=np.bincount(ix)
    ys=np.bincount(ix,weights=y);sq=np.bincount(ix,weights=y*y)
    diffs=np.array([np.bincount(ix,weights=(p-y)**2-(r-y)**2) for r in refs])
    ranks=[FrozenDiscrimination(y,p),*[FrozenDiscrimination(y,r) for r in refs]] if binary else []
    def draw(w):
        n=w@counts
        delta=diffs@w/n
        if not binary:
            var=w@sq/n-(w@ys/n)**2
            return np.column_stack([delta,-delta/var if var>0 else np.full(len(refs),np.nan)])
        scores=np.array([r.scores(w[ix]) for r in ranks])
        return np.column_stack([delta,scores[0]-scores[1:]])
    point=draw(np.ones(ng))
    rng=np.random.default_rng(seed)
    b=np.empty((reps,*point.shape));b[:]=np.nan
    if ng>=2:
        for i in range(reps):b[i]=draw(np.bincount(rng.integers(0,ng,ng),minlength=ng))
    return point,b


def curve_points(y,p,cal,classification,probability):
    """Calibration-set score deciles, test histograms, and full-score ROC / PR."""
    curves=[]
    def add(kind,x,z,n=None,observed=None):
        curves.append(dict(curve=kind,x=float(x),y=float(z),n=n,observed=observed))
    edges=np.unique(np.quantile(cal,np.arange(.1,1,.1)))
    bins=np.searchsorted(edges,p,side="right")
    for k in np.unique(bins):
        ok=bins==k;n=int(ok.sum());obs=float(y[ok].mean());avg=float(p[ok].mean())
        add("outcome by score bin",k+1,obs,n,avg)
        if probability or not classification:add("calibration",avg,obs,n)
    lo,hi=np.min(p),np.max(p)
    if hi<=lo:hi=lo+1e-8
    breaks=np.linspace(lo,hi,41)
    for label,mask in ([("case",y==1),("control",y==0)] if classification else [("all",np.ones(len(y),bool))]):
        ns,_=np.histogram(p[mask],bins=breaks)
        for x,n in zip((breaks[:-1]+breaks[1:])/2,ns):
            add("score distribution · "+label,x,n/max(1,mask.sum()),int(n))
    if classification and len(np.unique(y))==2:
        fpr,tpr,_=roc_curve(y,p);pr,rec,_=precision_recall_curve(y,p)
        for kind,x,z in (("ROC",fpr,tpr),("precision–recall",rec,pr)):
            for i in np.unique(np.linspace(0,len(x)-1,min(180,len(x)),dtype=int)):add(kind,x[i],z[i])
    return curves


def evaluate_method(t,name,y,p,cal,cohort,calpop,refs,cov,cfg):
    perf=[];pairs=[];curves=[]
    observed=np.isfinite(y)
    for anc in ["ALL",*sorted(cohort.ancestry.unique())]:
        mask=observed.copy()
        if anc!="ALL":mask &= cohort.ancestry.eq(anc).to_numpy()
        if not mask.any():continue
        pc=cal if anc=="ALL" else cal[calpop==anc]
        group=cohort.group_id.to_numpy()[mask];yy=y[mask];pp=p[mask]
        for endpoint in ([t,"ldl.high"] if t=="ldl" else [t]):
            bi=endpoint in ("t2dm","ldl.high")
            high=endpoint=="ldl.high"
            target=(yy>=cfg["evaluation"].get("ldl_threshold",4.9)).astype(float) if high else yy
            base=dict(trait=t,endpoint=endpoint,ancestry=anc,method=name)
            if high:
                scores=FrozenDiscrimination(target,pp).scores(np.ones(len(pp)))
                row=dict(auroc=scores[0],auprc=scores[1],n_cases=int(target.sum()),event_fraction=target.mean())
            else:row=metrics(target,pp,bi,fit_calibration=True,covariate=cov[mask] if cov is not None else None)
            perf.append(dict(**base,**{**row,"n":len(yy),"families":len(np.unique(group))}))
            others=[k for k in refs if k!=name]
            if others:
                reps=cfg["evaluation"]["bootstrap"];seed=cfg["cohort"]["seed"]
                # A high-LDL ranking score is in mmol/L. Its MSE against a 0/1 label has no interpretation.
                pt,bs=bootstrap(target,pp,[refs[k][mask] for k in others],group,bi,reps,seed)
                labels=("Brier","AUROC","AP") if bi else ("MSE","predictive_R2")
                for j,ref in enumerate(others):
                    for k,label in enumerate(labels):
                        if high and k==0:continue
                        low,up,n,status=interval(bs[:,j,k],reps)
                        pairs.append(dict(**base,reference=ref,metric=label,difference=pt[j,k],
                                          ci_low=low,ci_high=up,valid_draws=n,ci_status=status))
            if len(pc):
                curves.extend(dict(**base,**r) for r in curve_points(target,pp,pc,bi,bi and not high))
    return [pd.DataFrame(x) for x in (perf,pairs,curves)]


def read_scores(path,cols,ids):
    # Read numeric scores directly, keeping the participant IDs as strings.
    d=pd.read_csv(path,sep="\t",usecols=["id",*cols],dtype={"id":str})
    if d.id.isna().any() or d.id.duplicated().any():raise ValueError("Invalid prediction IDs: "+str(path))
    d=d.set_index("id").reindex(ids)
    if not np.isfinite(d.to_numpy(float)).all():raise ValueError("Incomplete CALIBRATION/TEST predictions: "+str(path))
    return d


def cache_result(folder,key,compute,names):
    receipt=folder/"cache.json"
    if receipt.is_file() and json.loads(receipt.read_text())==key and all((folder/(k+".tsv")).is_file() for k in names):return
    receipt.unlink(missing_ok=True)
    for k,frame in zip(names,compute(),strict=True):
        atomic_table(frame if len(frame.columns) else pd.DataFrame(columns=["trait"]),folder/(k+".tsv"))
    atomic_json(key,receipt)


def evaluate(cfg,root,t):
    paper=project_path(root,"grid.paper/main")
    cache=project_path(root,"grid.paper/logs")/"cache"
    cache=cache/t
    with file_lock(cache/"evaluate.lock"):
        src,records,skipped=discover(cfg,root,t)
        initial={k:stamp(v) for k,v in src.items()}
        inputs=[project_path(root,"tmp/cohort/participants.tsv.gz"),project_path(root,f"tmp/cohort/_private/{t}.test.tsv.gz")]
        cohort=read_table(inputs[0],["id","role","group_id","ancestry"]);check_roles(cohort)
        test=cohort[cohort.role.eq("test")].reset_index(drop=True)
        cal=cohort[cohort.role.eq("calibration")].reset_index(drop=True)
        ids=pd.Index(pd.concat([test.id,cal.id],ignore_index=True))
        y=aligned(read_table(inputs[1]),test.id).y.to_numpy(float)
        n=len(test)
        basekey=dict(inputs=[stamp(p) for p in inputs],evaluation=cfg["evaluation"],seed=cfg["cohort"]["seed"],
                     code=[stamp(p) for p in (Path(__file__),HERE/"evaluate2.py")])
        sources={};refs={};cov=None
        for owner,path in src.items():
            for name,col in predcols(path,owner).items():
                if name in sources:raise ValueError("Duplicate prediction name: "+name)
                if not re.fullmatch(r"[A-Za-z0-9_.-]+",name):raise ValueError("Invalid method identifier "+name)
                sources[name]=(owner,path,col)
        for name in (*REFS,"covariate","universal.global.auto","universal.global.fixed","universal.fixed","csx.auto","csx.meta","block.stack","direct.linear","direct.global_mlp","block.glmnet","blocknet.disco.global","block.glmnet.disco.offset"):
            if name in sources:
                _,path,col=sources[name];values=read_scores(path,[col],ids)[col].to_numpy(float)[:n]
                if name=="covariate":cov=values
                else:refs[name]=values
        if not "csx.population" in refs:raise ValueError("Missing CSx population comparator")
        # Catch a preserved Disco snapshot joined to a different CSx run.
        if "disco" in src and "csx.population.prediction" in pd.read_csv(src["disco"],sep="\t",nrows=0):
            both=read_scores(src["disco"],["csx.population.prediction"],test.id)
            if not np.allclose(both.iloc[:,0],refs["csx.population"],rtol=1e-12,atol=1e-12):raise ValueError("Disco snapshot disagrees with CSx")
        basekey["references"]=[stamp(src[k]) for k in ("csx","disco") if k in src]
        for name,(owner,path,col) in sources.items():
            comparators={r:v for r,v in refs.items() if r in REFS or
                    (name in ("universal.auto","universal.fixed") and r in
                     ("universal.global."+name.rsplit(".",1)[-1],"block.stack","csx.auto","csx.meta")) or
                    (name=="universal.auto" and r=="universal.fixed") or
                    (r.startswith("direct.") and name=="direct.blocknet") or
                    (r=="block.glmnet" and name in {"blocknet.disco","block.glmnet.disco","block.stack"}) or
                    (r=="blocknet.disco.global" and name=="blocknet.disco") or
                    (r=="block.glmnet.disco.offset" and name=="block.glmnet.disco")}
            key=dict(basekey,prediction=stamp(path),column=col,
                     matched_controls={r:stamp(sources[r][1]) for r in comparators})
            def compute():
                p=read_scores(path,[col],ids)[col].to_numpy(float)
                if t=="t2dm" and not ((p>=0)&(p<=1)).all():raise ValueError("T2DM prediction must be a probability: "+name)
                return evaluate_method(t,name,y,p[:n],p[n:],test,cal.ancestry.to_numpy(),comparators,cov,cfg)
            cache_result(cache/name,key,compute,ARTIFACTS[:3])
        extra=[]
        if "abm" in src:
            from evaluate2 import evaluate_predictions
            path=cache/"abm-screening";key=dict(basekey,abm=stamp(src["abm"]))
            def screening():
                # Only the stable comparator set participates; adding an unrelated learner costs no second ABM bootstrap.
                frame=read_table(src.get("disco",src["csx"])).rename(columns=short)
                d=read_table(src["abm"]);path.mkdir(parents=True,exist_ok=True)
                p,_,_=evaluate_predictions(cohort,read_table(inputs[1]),frame,d,path,trait=t,
                               bootstrap=cfg["evaluation"]["bootstrap"],seed=cfg["cohort"]["seed"],workers=cfg["runtime"]["threads"])
                return [p]
            cache_result(path,key,screening,["screening"]);extra.append(["screening",str(path/"screening.tsv")])
        if all(k in src for k in ("direct.linear","direct.blocknet","abm")) and "blocknet.twins" in records:
            from ldl_report import twin_audit
            path=cache/"blocknet-twins";twin=records["blocknet.twins"].parent/"predictions.tsv.gz"
            key=dict(basekey,models=[stamp(src[k]) for k in ("direct.linear","direct.blocknet","abm")]+[stamp(twin),stamp(HERE/"ldl_report.py")])
            def screen_twins():
                p=read_table(src["direct.blocknet"]).merge(read_table(src["direct.linear"]),on="id",validate="one_to_one")
                path.mkdir(parents=True,exist_ok=True)
                twin_audit(cohort,read_table(inputs[1]),p,read_table(src["abm"]),read_table(twin),path,
                           cfg["evaluation"]["bootstrap"],cfg["cohort"]["seed"],cfg["runtime"]["threads"],t)
                return [pd.read_csv(path/f"{t}.blocknet_twins.tsv",sep="\t").assign(trait=t)]
            cache_result(path,key,screen_twins,["twins"]);extra.append(["ldl_twins",str(path/"twins.tsv")])
        current,ready,_=discover(cfg,root,t)
        if src!=current or records!=ready or any(stamp(p)!=initial[k] for k,p in current.items()):
            raise RuntimeError(f"{t}: producers changed during evaluation; repeat the incremental paper merge")
        atomic_json(dict(methods=list(sources),evaluation=cfg["evaluation"],sources={k:stamp(v) for k,v in src.items()},extra=extra,
                         records={k:stamp(v) for k,v in records.items()},omitted=skipped),cache/"inventory.json")


def publish(cfg,root):
    from methods import write_index
    write_index(root)
    from contextlib import ExitStack
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paper=project_path(root,"grid.paper/main");tables=paper/"results";reports=paper/"reports"
    tables.mkdir(parents=True,exist_ok=True);reports.mkdir(parents=True,exist_ok=True)
    cache=project_path(root,"grid.paper/logs")/"cache"
    bundles={k:[] for k in (*ARTIFACTS,"ldl_twins","availability","ld_maps","cohort")}
    with file_lock(cache/"merge.lock"),ExitStack() as locks:
        traits=[p.parent.name for p in sorted(cache.glob("*/inventory.json"))]
        for t in traits:
            inv=json.loads((cache/t/"inventory.json").read_text())
            src,records,_=discover(cfg,root,t)
            # Another trait may have gained or lost methods since the last incremental report.
            if (inv["sources"]!={k:stamp(v) for k,v in src.items()}
                or inv["records"]!={k:stamp(v) for k,v in records.items()}
                or inv.get("evaluation")!=cfg["evaluation"]):
                evaluate(cfg,root,t)
            locks.enter_context(file_lock(cache/t/"evaluate.lock",exclusive=False))
            inv=json.loads((cache/t/"inventory.json").read_text())
            src,records,skip=discover(cfg,root,t)
            # Retain only methods with evaluated source data.
            current={k for source,f in src.items() for k in predcols(f,source)}
            for m in inv["methods"]:
                if m in current:
                    for k in ARTIFACTS[:3]:bundles[k].append(pd.read_csv(cache/t/m/(k+".tsv"),sep="\t"))
            for k,p in inv["extra"]:
                if ("abm" in records if k=="screening" else all(n in records for n in ("blocknet.twins","abm","direct.blocknet","direct.linear"))):
                    bundles[k].append(pd.read_csv(p,sep="\t"))
            for m in dict.fromkeys([*records,*src,*skip]):
                bundles["availability"].append(pd.DataFrame([dict(trait=t,method=m,status=skip.get(m,"complete"))]))
        from manuscript import extract_context
        from collections import defaultdict
        context=defaultdict(list)
        extract_context(root,traits,context)
        for frame in context["model_settings"]:
            for i,r in frame.iterrows():
                inv=json.loads((cache/r.trait/"inventory.json").read_text())
                frame.loc[i,"optional_methods_included"]=",".join(inv["sources"])
        for k in ("discovery","anchor_training","score_coverage","csx_weights","model_settings","ld_maps"):
            bundles[k]=context[k]
        co=read_table(project_path(root,"tmp/cohort/participants.tsv.gz"),["role","ancestry","group_id"])
        bundles["cohort"].append(co.groupby(["role","ancestry"]).group_id.agg(n="size",families="nunique").reset_index())
        combined={k:pd.concat(v,ignore_index=True).drop_duplicates() for k,v in bundles.items() if v}
        if "performance" in combined:
            h=combined["performance"].query("endpoint=='height' and ancestry=='ALL'").copy()
            if len(h):
                h=h[["method","r2","rmse","n"]].copy()
                h["r2_minus_target"]=h.r2-.75;h["target_reached"]=h.r2>.75
                combined["height_target"]=h
        bundles["height_target"]=[]
        for k in bundles:
            if k in combined:atomic_table(combined[k],tables/(k+".tsv"))
            else:(tables/(k+".tsv")).unlink(missing_ok=True)
        points=combined["performance"];pairs=combined["contrasts"];curves=combined["curves"]
        for e,df in points[points.ancestry.eq("ALL")].groupby("endpoint"):
            metric="auroc" if e in {"t2dm","ldl.high"} else "r2"
            df=df.sort_values(metric)
            fig,ax=plt.subplots(figsize=(8,max(4,.32*len(df)+1.2)),layout="constrained")
            ax.scatter(df[metric],df.method,color="#14746f")
            ax.set(xlabel="AUROC" if metric=="auroc" else "Predictive R²",title=e+" · held-out performance")
            ax.spines[["right","top"]].set_visible(False);fig.savefig(reports/(e+".png"),dpi=210);plt.close(fig)
        for e,dd in curves[curves.ancestry.eq("ALL")].groupby("endpoint"):
            kinds=["ROC","precision–recall"] if e in ("ldl.high","t2dm") else ["calibration","outcome by score bin"]
            fig,axes=plt.subplots(1,2,figsize=(10.5,4.5),layout="constrained")
            for ax,kind in zip(axes,kinds):
                view=dd.loc[dd.curve.eq(kind)&dd.method.isin(["csx.population","disco","block.glmnet","block.glmnet.disco","blocknet.disco","block.stack","direct.blocknet"])]
                for m,g in view.groupby("method"):ax.plot(g.x,g.y,label=m,linewidth=1.5)
                labels={"ROC":("False-positive rate","Sensitivity"),"precision–recall":("Recall","Precision"),
                        "calibration":("Mean prediction","Observed mean"),"outcome by score bin":("Score bin","Observed mean")}
                ax.set(xlabel=labels[kind][0],ylabel=labels[kind][1],title=kind)
                ax.spines[["right","top"]].set_visible(False)
            axes[0].legend(fontsize=7);fig.suptitle(e);fig.savefig(reports/(e+".curves.png"),dpi=210);plt.close(fig)
        # All HTML data are aggregate only; the report also works completely offline.
        keys=["trait","endpoint","ancestry","method"]
        pc=[c for c in (*keys,"r2","auroc","auprc","rmse","partial_r2_vs_covariate","n","n_cases","brier","mae") if c in points]
        export={"performance":json.loads(points[pc].to_json(orient="records",double_precision=9)),
                "contrasts":json.loads(pairs.to_json(orient="records",double_precision=9)),
                "curves":json.loads(curves.to_json(orient="records",double_precision=7))}
        if "ldl_twins" in combined:
            fields=["ancestry","policy","coverage","n_evaluable","median_latent_twins","mse","selected_minus_all_mse","ci_low","ci_high","selected_minus_risk_only_mse"]
            export["twins"]=json.loads(combined["ldl_twins"][fields].to_json(orient="records"))
        import gzip,base64
        data=json.dumps(export,ensure_ascii=True,separators=(",",":"),allow_nan=False).encode()
        html=(HERE/"methods_report.html").read_text().replace("__THRESHOLD__",str(cfg["evaluation"].get("ldl_threshold",4.9))).replace("__DATA__",base64.b64encode(gzip.compress(data)).decode())
        temp=reports/"index.partial";temp.write_text(html);temp.replace(reports/"index.html")
        (paper/"index.html").write_text('<meta charset="utf-8"><meta http-equiv="refresh" content="0; url=reports/index.html"><a href="reports/index.html">GRID report</a>')
