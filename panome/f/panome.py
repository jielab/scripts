#!/usr/bin/env python3
"""Panome: individuals, molecular neighborhoods and conditional prediction."""
from pathlib import Path
import argparse, copy, importlib.metadata, importlib.util, json, os, platform, shutil, time
from common import VERSION, STAGES, words, dump, digest, fingerprints, log, run_lock, stage_log
HOME = Path(__file__).resolve().parents[1]

def parser():
    class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass
    p = argparse.ArgumentParser(prog='./panome.sh',description=__doc__,formatter_class=HelpFormatter,
        epilog="""Usage examples（可直接复制）：
  # 检查真实 PPP 输入及环境
  ./panome.sh -Y cvd_cad --biom prot --preflight

  # PPP 主分析；默认输出 /mnt/d/analysis/panome/cvd_cad/prot/v3/
  ./panome.sh -Y cvd_cad --biom prot

  # 5000 人试运行，结果独立保存
  ./panome.sh -Y cvd_cad --biom prot --run-name pilot_v3 --max-samples 5000

  # 分别分析蛋白组和代谢组
  ./panome.sh -Y cvd_cad --biom prot,met

  # 汇总已有 v3 结果；不重新训练
  ./panome.sh final -Y cvd_cad --biom prot
""")
    p.add_argument("module",nargs="?",choices=["run","final"],default="run")
    p.add_argument("-Y","--Y","--trait",dest="trait",default="cvd_cad")
    p.add_argument("--biom",default="prot",help="prot, met, or comma-separated independent runs")
    p.add_argument("--ukb-phe",default=os.getenv("UKB_PHE","/mnt/d/data/ukb/phe"))
    p.add_argument("--phe-file",help="Default <ukb-phe>/Rdata/all.rds")
    p.add_argument("--omics-file",help="Explicit matrix, or a path containing {biom}")
    p.add_argument("--input-source",choices=["raw","cleaned"],default="raw")
    p.add_argument("--met-input",choices=["auto","raw","named"],default="auto")
    p.add_argument("--met-map",help="Default <ukb-phe>/common/met.lst")
    p.add_argument("--analysis-root",default=os.getenv("PANOME_ANALYSIS_ROOT","/mnt/d/analysis/panome"))
    p.add_argument("--run-name",default="v3")
    p.add_argument("--outcome-type",choices=["survival","quantitative"],default="survival")
    p.add_argument("--target-col",help="Quantitative outcome column; defaults to Y")
    p.add_argument("--id-col",default="eid")
    p.add_argument("--baseline-col",default="date_attend")
    p.add_argument("--diagnosis-col",help="Default fod_icd10_<Y>")
    p.add_argument("--death-col",default="date_death")
    p.add_argument("--lost-col",default="date_lost")
    p.add_argument("--end-date",default=os.getenv("DATE_FOLLOW_END","2023-04-01"))
    p.add_argument("--healthy-date-cols",default="")
    p.add_argument("--disease-evidence-col",default="")
    p.add_argument("--covariates",default="age,sex,tdi,PC1,PC2,center")
    p.add_argument("--residualize",default=None,help="Default age,sex,<biom>.plate; empty string disables")
    p.add_argument("--categorical",default=None)
    p.add_argument("--group-col",default="",help="Complete family/relatedness-component ID")
    p.add_argument("--transform",choices=["auto","none","log1p"],default="auto")
    p.add_argument("--exclude-features",default="",help="Exact assay names; no automatic GDF15/BNP removal")
    for flag,default in [("feature-missing",.2),("sample-missing",.2),("corruption",.1),
                         ("min-state-fraction",.02),("neural-lr",.001),("neural-weight-decay",.001),
                         ("primary-horizon",10),("pair-caliper",.01)]:
        p.add_argument("--"+flag,type=float,default=default)
    integers = {"latent":20,"hidden":512,"epochs":100,"patience":12,"batch-size":256,
        "neural-epochs":100,"time-bins":8,"token-modules":32,"token-dim":32,
        "attention-heads":4,"attention-layers":2,"neighbors":30,"graph-dims":10,
        "max-states":20,"stability-repeats":10,"min-expert-events":30,"pwas-top":30,
        "coxnet-alphas":20,"tree-estimators":100,"min-events":20,"bootstrap":200,
        "attribution-samples":500,"ig-steps":64,"max-pairs":100,"seed":2026,
        "cores":int(os.getenv("PANOME_THREADS","4")),"log-every":10,"demo-features":80,"max-samples":0}
    for flag,default in integers.items():
        p.add_argument("--"+flag,type=int,default=default)
    p.add_argument("--device",choices=["cpu","cuda","auto"],default="auto")
    p.add_argument("--neural",action=argparse.BooleanOptionalAction,default=True)
    p.add_argument("--tree",action=argparse.BooleanOptionalAction,default=True)
    p.add_argument("--resolutions",default="0.1,0.2,0.5,1,1.5")
    p.add_argument("--horizons",default="5,10")
    p.add_argument("--landmarks",default="0,2,5")
    p.add_argument("--match-model",default="clinical_elasticnet")
    p.add_argument("--pgs",action="store_true",help="Optional exact-matched PGS overlay")
    p.add_argument("--pgs-file",help="Else layer.pgs.rds if available, then all.rds")
    p.add_argument("--pgs-source",default="unspecified")
    p.add_argument("--pgs-overlap",choices=["none","unknown","yes"],default="unknown")
    p.add_argument("--r-bin",default=os.getenv("R_BIN","Rscript"))
    p.add_argument("--from",dest="from_stage",choices=STAGES,default=STAGES[0])
    p.add_argument("--to",dest="to_stage",choices=STAGES,default=STAGES[-1])
    p.add_argument("--steps",default="")
    for flag in ["replace","full-input-hash","preflight","dry-run","demo"]:
        p.add_argument("--"+flag,action="store_true")
    p.add_argument("--version",action="version",version=VERSION)
    return p

def resolve(a,layer):
    a = copy.deepcopy(a)
    a.biom = layer
    root = Path(a.ukb_phe)
    a.phe_file = a.phe_file or str(root/"Rdata/all.rds")
    a.met_map = a.met_map or str(root/"common/met.lst")
    a.target_col = a.target_col or a.trait
    a.diagnosis_col = a.diagnosis_col or f"fod_icd10_{a.trait}"
    if a.omics_file:
        a.omics_file = a.omics_file.replace("{biom}",layer)
    elif a.input_source=="cleaned":
        a.omics_file = str(root/f"Rdata/{layer}.rds")
    else:
        a.omics_file = str(root/("rap/raw/prot.tab.gz" if layer=="prot" else "rap/met.tab.gz"))
    if a.met_input=="auto":
        a.met_input = "raw" if layer=="met" and Path(a.omics_file).name=="met.tab.gz" else "named"
    a.residualize = a.residualize if a.residualize is not None else f"age,sex,{layer}.plate"
    a.categorical = a.categorical if a.categorical is not None else f"sex,center,{layer}.plate"
    if a.transform=="auto":
        a.transform = "log1p" if layer=="met" and not a.demo else "none"
    for key in ["resolutions","horizons","landmarks"]:
        setattr(a,key,[float(x) for x in words(getattr(a,key))])
    if a.pgs:
        candidate = root/f"Rdata/{layer}.pgs.rds"
        a.pgs_file = (a.pgs_file.replace("{biom}",layer) if a.pgs_file
                      else str(candidate if candidate.is_file() else a.phe_file))
    for value in [a.trait,a.biom,a.run_name]:
        if not value or "/" in value or "\\" in value or value in [".",".."]:
            raise ValueError("Y, biom and run-name must be simple directory names")
    for key in ["latent","hidden","epochs","patience","batch_size","neural_epochs","time_bins",
                "token_modules","token_dim","attention_heads","attention_layers","neighbors",
                "graph_dims","max_states","min_expert_events","pwas_top","coxnet_alphas",
                "tree_estimators","min_events","ig_steps","max_pairs","cores","log_every","demo_features"]:
        if getattr(a,key)<1:
            raise ValueError(f"{key} must be positive")
    if a.latent<2 or a.token_dim%a.attention_heads:
        raise ValueError("latent >= 2 and token-dim divisible by attention-heads required")
    for key in ["bootstrap","stability_repeats","attribution_samples","max_samples"]:
        if getattr(a,key)<0:
            raise ValueError(f"{key} must be nonnegative")
    for key in ["feature_missing","sample_missing","min_state_fraction"]:
        if not 0<getattr(a,key)<1:
            raise ValueError(f"{key} must be in (0,1)")
    if not 0<=a.corruption<1 or a.pair_caliper<=0 or a.neural_lr<=0 or a.neural_weight_decay<0:
        raise ValueError("Invalid corruption, caliper or optimizer parameter")
    if not a.resolutions or min(a.resolutions)<=0 or not a.horizons or min(a.horizons)<=0 or min(a.landmarks,default=0)<0:
        raise ValueError("Invalid resolution/horizon/landmark list")
    if a.primary_horizon not in a.horizons and a.outcome_type=="survival":
        raise ValueError("primary-horizon must be in horizons")
    forbidden = {a.target_col} if a.outcome_type=="quantitative" else {
        a.diagnosis_col,a.death_col,a.lost_col,a.disease_evidence_col,"event","time"}
    if forbidden.intersection(words(a.covariates)+words(a.residualize)):
        raise ValueError("Outcome/follow-up columns cannot be covariates or residualization inputs")
    return a

def preflight(a):
    packages = {"numpy":"numpy","pandas":"pandas","scipy":"scipy","scikit-learn":"sklearn",
        "scikit-survival":"sksurv","torch":"torch","igraph":"igraph","leidenalg":"leidenalg",
        "pynndescent":"pynndescent","joblib":"joblib","matplotlib":"matplotlib","openpyxl":"openpyxl",
        "lifelines":"lifelines","statsmodels":"statsmodels","pyreadr":"pyreadr"}
    missing = [name for name,mod in packages.items() if importlib.util.find_spec(mod) is None]
    if missing:
        raise RuntimeError("Missing dependencies: "+", ".join(missing))
    versions = {name:importlib.metadata.version(name) for name in packages}
    files = [] if a.demo else [a.phe_file,a.omics_file]
    if not a.demo and a.biom=="met" and a.met_input=="raw":
        files.append(a.met_map)
    if a.pgs:
        files.append(a.pgs_file)
    absent = [file for file in files if not Path(file).is_file()]
    if absent:
        raise FileNotFoundError("Missing inputs: "+", ".join(absent))
    if a.device=="cuda":
        import torch
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but unavailable")
    log("READY",a.biom,f"Python {platform.python_version()}; outcome={a.outcome_type}")
    return versions,files

def load_data(root):
    import numpy as np
    import pandas as pd
    cfg = json.loads((root/"manifest.json").read_text())["config"]
    d = root/"s2_preprocess"
    return (np.load(d/"x.npy",mmap_mode="r"),np.load(d/"observed.npy",mmap_mode="r"),
            pd.read_csv(d/"cohort.csv",dtype={cfg["id_col"]:str}),
            (d/"features.txt").read_text().splitlines())

def preprocess_stage(a,out,root):
    import numpy as np
    import pandas as pd
    import joblib
    from preprocess import outcomes,split_people,MolecularPreprocessor
    src = root/"s1_prepare"
    x = np.load(src/"raw.npy",mmap_mode="r")
    p = pd.read_csv(src/"phenotype.csv",dtype={a.id_col:str})
    p,audit = outcomes(p,a)
    keep = p.eligible.to_numpy()
    p,x = p.loc[keep].reset_index(drop=True),np.asarray(x[keep])
    p["split"] = split_people(p,a)
    train = p.split.to_numpy()=="train"
    finite = np.where(np.isfinite(x[train]),x[train],np.nan)
    retained = (np.mean(~np.isfinite(finite),axis=0)<a.feature_missing)&(np.nanstd(finite,axis=0)>1e-8)
    if retained.sum()<3:
        raise ValueError("Too few training-defined molecular features")
    row_keep = np.mean(~np.isfinite(x[:,retained]),axis=1)<a.sample_missing
    audit["sample_missing_exclusions"] = int((~row_keep).sum())
    p,x = p.loc[row_keep].reset_index(drop=True),x[row_keep]
    train = p.split.to_numpy()=="train"
    summary = {}
    for split in ["train","validation","test"]:
        sub = p[p.split==split]
        summary[split] = dict(n=len(sub))
        if len(sub)<30:
            raise ValueError(f"{split}: fewer than 30 eligible people")
        if a.outcome_type=="survival":
            summary[split]["events"] = int(sub.event.sum())
            if sub.event.sum()<a.min_events:
                raise ValueError(f"{split}: {sub.event.sum()} events; minimum {a.min_events}")
        elif sub.target.nunique()<2:
            raise ValueError(f"{split}: target has no variation")
    if min(retained.sum(),train.sum()-1)<a.latent:
        raise ValueError("latent exceeds the training matrix rank bound")
    pre = MolecularPreprocessor(a.feature_missing,words(a.residualize),
        words(a.categorical),a.transform).fit(x[train],p.loc[train],allowed=retained)
    values,mask = pre.transform(x,p)
    np.save(out/"x.npy",values); np.save(out/"observed.npy",mask)
    p.to_csv(out/"cohort.csv",index=False)
    names = np.array((src/"features.txt").read_text().splitlines())
    features = names[pre.keep]
    (out/"features.txt").write_text("\n".join(features)+"\n")
    pd.DataFrame(dict(feature=names,retained=np.isin(np.arange(len(names)),pre.keep))).to_csv(out/"feature_qc.csv",index=False)
    joblib.dump(pre,out/"preprocessor.joblib")
    audit.update(n_after_qc=len(p),features=len(features),split_summary=summary,
                 group_split=bool(a.group_col),outcome_type=a.outcome_type)
    dump(out/"cohort_audit.json",audit)
    log("COHORT",a.biom,f"N={len(p)}, features={len(features)}, splits={summary}")

def representation_stage(a,out,root):
    import numpy as np
    import joblib
    from sklearn.decomposition import PCA
    from representation import Autoencoder,MolecularTransformer,fit_reconstruction,groups_from_pca,initialize
    x,observed,p,features = load_data(root)
    tr = p.split.to_numpy()=="train"
    pca = PCA(n_components=a.latent,svd_solver="randomized",random_state=a.seed).fit(x[tr])
    np.save(out/"pca.npy",pca.transform(x).astype("float32")); joblib.dump(pca,out/"pca.joblib")
    initialize(a)
    ae,summary = fit_reconstruction(Autoencoder(x.shape[1],a.latent,a.hidden),
                                    x,observed,p.split.to_numpy(),a,out,"ae")
    np.save(out/"ae.npy",ae)
    summaries = dict(ae=summary,pca_explained_variance=float(pca.explained_variance_ratio_.sum()))
    if a.neural:
        initialize(a)
        groups = groups_from_pca(pca,a.token_modules,a.seed)
        model = MolecularTransformer(x.shape[1],groups,a.token_dim,a.attention_heads,a.attention_layers)
        embedding,stats = fit_reconstruction(model,x,observed,p.split.to_numpy(),a,out,"transformer")
        np.save(out/"transformer.npy",embedding)
        dump(out/"token_modules.json",[dict(module=j+1,features=[features[k] for k in group])
                                      for j,group in enumerate(groups)])
        summaries["transformer"] = stats
    dump(out/"summary.json",summaries)

def graph_stage(a,out,root):
    import numpy as np
    import joblib
    from graph import PersonAtlas
    _,_,p,_ = load_data(root)
    z = np.load(root/"s3_representation/ae.npy")
    tr = p.split.to_numpy()=="train"
    groups = p[a.group_col].to_numpy(str) if a.group_col else None
    atlas = PersonAtlas().fit(z[tr],p.loc[tr,a.id_col].to_numpy(str),
                              None if groups is None else groups[tr],a,out)
    projected = atlas.project(z,p[a.id_col].to_numpy(str),groups)
    np.savez_compressed(out/"graph_coordinates.npz",**projected)
    joblib.dump(atlas,out/"atlas.joblib")
    p.loc[tr,[a.id_col]].assign(discovery_state=atlas.labels+1).to_csv(out/"graph_node_ids.csv",index=False)
    dump(out/"summary.json",atlas.summary)

def predict_stage(a,out,root):
    import numpy as np
    import joblib
    from prediction import run_predictions
    from pgs import pgs_overlay
    x,mask,p,features = load_data(root)
    src = root/"s3_representation"
    z,pca = np.load(src/"ae.npy"),np.load(src/"pca.npy")
    transformer = np.load(src/"transformer.npy") if a.neural else None
    graph = dict(np.load(root/"s4_graph/graph_coordinates.npz"))
    atlas = joblib.load(root/"s4_graph/atlas.joblib")
    c = run_predictions(x,mask,p,z,pca,transformer,graph,atlas,features,a,out)
    pgs_overlay(x,mask,p,c,features,a,out)

def report_stage(a,out,root):
    import numpy as np
    import joblib
    from interpretation import interpret
    x,mask,p,features = load_data(root)
    z = np.load(root/"s3_representation/ae.npy")
    graph = dict(np.load(root/"s4_graph/graph_coordinates.npz"))
    atlas = joblib.load(root/"s4_graph/atlas.joblib")
    clinical = joblib.load(root/"s5_predict/clinical_preprocessor.joblib").transform(p)
    interpret(x,mask,p,z,graph,atlas,features,clinical,a,root,out)

def completed(directory,signature):
    path = directory/"DONE.json"
    if not path.exists():
        return False
    done = json.loads(path.read_text())
    if done["signature"]!=signature:
        raise ValueError(f"Stage signature mismatch: {directory}")
    actual = fingerprints([directory/name for name in done["outputs"]])
    if actual != done["fingerprints"]:
        raise ValueError(f"Stage outputs changed or missing: {directory}")
    return True

def run_one(a):
    root = Path(a.analysis_root).resolve()/a.trait/a.biom/a.run_name
    log("OUTPUT",a.biom,str(root))
    if a.module=="final":
        if a.replace:
            raise ValueError("final does not replace analysis data")
        if a.dry_run:
            return
        if not (root/"manifest.json").is_file():
            raise FileNotFoundError(root/"manifest.json")
        with run_lock(root):
            from final import export_final
            export_final(root)
        return
    stages = words(a.steps) if a.steps else STAGES[STAGES.index(a.from_stage):STAGES.index(a.to_stage)+1]
    if not stages or len(stages)!=len(set(stages)) or any(s not in STAGES for s in stages):
        raise ValueError("Invalid stage selection")
    stages = sorted(stages,key=STAGES.index)
    if a.dry_run:
        print(json.dumps(dict(stages=stages,phenotype=a.phe_file,omics=a.omics_file,
            met_map=a.met_map if a.biom=="met" else None,transform=a.transform,
            covariates=a.covariates,residualize=a.residualize,pgs=a.pgs_file),indent=2))
        return
    versions,files = preflight(a)
    if a.preflight:
        return
    excluded = {"module","from_stage","to_stage","steps","replace","preflight","dry_run"}
    config = {key:value for key,value in vars(a).items() if key not in excluded}
    code = fingerprints(list((HOME/"f").glob("*.py"))+list((HOME/"f").glob("*.R"))+[HOME/"panome.sh"])
    payload = dict(schema=3,version=VERSION,config=config,versions=versions,
                   inputs=fingerprints(files,a.full_input_hash),code=code)
    signature = digest(payload)
    root.mkdir(parents=True,exist_ok=True)
    from threadpoolctl import threadpool_limits
    with run_lock(root),threadpool_limits(limits=a.cores):
        manifest = root/"manifest.json"
        if manifest.exists() and not a.replace:
            old = json.loads(manifest.read_text())
            if old.get("signature")!=signature:
                raise ValueError("Configuration/code/input differs. Use a new --run-name or explicit --replace. Old caches are not imported.")
        elif any(root.glob("s[1-6]_*")) and not a.replace:
            raise ValueError("Existing stages lack a matching manifest; choose another run-name")
        if a.replace:
            if stages[0]!=STAGES[0]:
                raise ValueError("--replace requires s1_prepare")
            for stage in STAGES+["publication"]:
                directory = root/stage
                if directory.exists():
                    if directory.is_symlink():
                        raise ValueError("Refusing to replace a symlinked output directory")
                    shutil.rmtree(directory)
        dump(manifest,dict(signature=signature,**payload))
        from io_data import prepare
        funcs = [lambda a,out,root:prepare(a,out),preprocess_stage,representation_stage,
                 graph_stage,predict_stage,report_stage]
        for name in stages:
            ix = STAGES.index(name)
            for prior in STAGES[:ix]:
                if not completed(root/prior,signature):
                    raise ValueError(f"{name} needs completed {prior}")
            directory = root/name
            if completed(directory,signature):
                log("REUSE",name)
                continue
            # An interrupted stage may contain models that a retry no longer produces.
            # Never include those stale artifacts in the new completion manifest.
            if directory.is_symlink():
                raise ValueError("Refusing to rerun a symlinked stage directory")
            if directory.exists():
                shutil.rmtree(directory)
            directory.mkdir()
            start = time.monotonic()
            with stage_log(name):
                funcs[ix](a,directory,root)
            outputs = sorted(p.name for p in directory.iterdir() if p.is_file() and p.name!="DONE.json")
            dump(directory/"DONE.json",dict(signature=signature,seconds=time.monotonic()-start,
                outputs=outputs,fingerprints=fingerprints([directory/name for name in outputs])))
        if "s6_report" in stages:
            from final import export_final
            export_final(root)
    log("DONE","panome",str(root))

def main():
    original = parser().parse_args()
    layers = words(original.biom)
    if not layers or len(set(layers))!=len(layers) or any(b not in ["prot","met"] for b in layers):
        raise ValueError("biom must be prot, met or prot,met")
    if len(layers)>1 and original.omics_file and "{biom}" not in original.omics_file:
        raise ValueError("For several layers use --omics-file with {biom}, or defaults")
    for layer in layers:
        run_one(resolve(original,layer))

if __name__=="__main__":
    try:
        main()
    except Exception as exc:
        log("ERROR","panome",str(exc))
        raise
