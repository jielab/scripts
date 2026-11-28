"""Producer manifests and the GRID dependency graph. Completed scientific fits are independent."""
from pathlib import Path
import json

from core2 import atomic_json, read_table, stamp
from file_lock import file_lock
from output_layout import project_path,trait_path

RECEIPT_FORMAT = "grid-method-20261127"

PARENTS={
 "universal.prepare.auto":(), "universal.auto":("universal.prepare.auto",),
 "universal.prepare.fixed":(), "universal.fixed":("universal.prepare.fixed",),
 "prepare":(), "architecture":("prepare",),
 "block.glmnet":("prepare",), "block.lightgbm":("prepare",),
 "disco.weights":("prepare",),
 "block.glmnet.disco":("prepare","disco.weights"),"blocknet.disco":("prepare","disco.weights"),"blocknet.disco.global":("prepare","disco.weights"),
 "block.stack":("block.glmnet","block.glmnet.disco","blocknet.disco"),
 "abm":("prepare",),"abm.deep":("prepare",),"twins":("abm",),
 "direct.linear":("prepare",),"direct.global_mlp":("prepare",),
 "direct.blocknet":("prepare",),"direct.regional_ridge":("prepare",),
 "blocknet.twins":("direct.blocknet","abm"),"robustness":("abm",)}
DEEP={"blocknet.disco","blocknet.disco.global","abm.deep","direct.blocknet","direct.global_mlp"}
ALIASES={"block_models":["block.glmnet","block.lightgbm"],"abm_deep":["abm.deep"],
         "direct_blocknet":["direct.blocknet","direct.linear"],"blocknet_ablation":["direct.global_mlp","direct.regional_ridge"],
         "blocknet_twins":["blocknet.twins"]}

PUBLIC_METHODS = tuple(m for m in PARENTS if m not in {"prepare", "disco.weights"} and ".prepare." not in m)


def selection(value):
    if value is None or value.strip()=="all":
        return [m for m in PUBLIC_METHODS if m!="robustness"]
    names=[]
    for item in value.split(","):
        for name in ALIASES.get(item.strip(),[item.strip()]):
            if name not in PUBLIC_METHODS:
                raise ValueError(f"Unknown analysis method {name!r}. Choose from "+", ".join(PUBLIC_METHODS)+
                                 ". Use ./grid.paper.sh run to generate the report.")
            if name not in names:names.append(name)
    return names


def folder(root,trait,method):
    return trait_path(root,trait,"regions" if method=="prepare" else method)

def marker(root,trait,method):
    return folder(root,trait,method)/"complete.json"

def source_files(method):
    extra={
      "prepare":["genomics.py","consensus3.py","scoring_parallel.py","ancestry_regions.py"],
      "disco.weights":["disco_weights.py","disco2.py","benchmark2.py"],"architecture":["report2.py"],"twins":["report2.py"],"abm":["abm2.py","block_neighbors.py"],
      "abm.deep":["block_deep.py","block_neighbors.py"],"blocknet.twins":["blocknet_twins.py"],
      "robustness":["abm2.py","evaluate2.py"]}
    if method=="block.stack":return [Path(__file__).with_name("block_stack.py")]
    files=extra.get(method,[])
    if method.startswith("direct."):files=["direct_blocknet.py","blocknet.py","blocknet_ablation.py"]
    if method.startswith("block.") or method.startswith("blocknet.disco"):files=["block_models.py","block_glmnet.R"]
    if method.endswith(".disco") or method.startswith("blocknet.disco"):files=files+["disco_blocks.py","blocknet.py"]
    return [Path(__file__).with_name(f) for f in files]

def contract(cfg,root,t,m):
    from driver2 import dependencies
    if m.startswith("universal."):
        from universal_blocks import inputs_for
        return inputs_for(cfg,root,t,m)
    if m=="prepare":return dependencies(cfg,root,"regions",t)
    sec=["cohort","benchmark"]
    if m.startswith("direct.") or m.startswith("blocknet.disco"):sec+=["direct_blocknet"]
    if m in DEEP:sec+=["block_deep"]
    if (m.endswith(".disco") or m.startswith("blocknet.disco")) and "disco_blocks" in cfg:sec+=["disco_blocks"]
    if m.startswith("block."):sec+=["block_models"]
    if m in ("abm","twins","robustness"):sec+=["abm"]
    if m=="blocknet.twins":sec+=["ldl_twins"]
    files=[project_path(root,"tmp/cohort/participants.tsv.gz"),
           project_path(root,f"tmp/cohort/{t}.development.tsv.gz"),
           trait_path(root,t,"baseline.predictions.tsv.gz"),*source_files(m)]
    if m in {"abm","block.stack","disco.weights","robustness"} or m.endswith(".disco") or m.startswith("blocknet.disco"):files+=[trait_path(root,t,"predictions.tsv.gz")]
    if m=="disco.weights":
        files.extend(trait_path(root,t,k) for k in ("disco/model.joblib","benchmark.json"))
    files.extend(marker(root,t,p) for p in PARENTS[m])
    return dict(method=m,settings={k:cfg[k] for k in sec},files=[stamp(p) for p in files])

def published(path):
    if not path.is_file():return None
    doc=json.loads(path.read_text())
    if doc.get("format") != RECEIPT_FORMAT:
        return None
    for f,record in doc["outputs"].items():
        target=path.parent/f
        if not target.is_file() or stamp(target)!=record:return None
    return doc

def run(cfg,root,t,m,force=False):
    from driver2 import execute,run_abm,regional_inputs
    from types import SimpleNamespace
    out=folder(root,t,m);out.mkdir(parents=True,exist_ok=True)
    lock=project_path(root,f"tmp/locks/{t}.{m}.lock")
    mark=marker(root,t,m)
    from contextlib import ExitStack
    with ExitStack() as locks:
        for p in sorted(PARENTS[m]):
            locks.enter_context(file_lock(project_path(root,f"tmp/locks/{t}.{p}.lock"),exclusive=False))
            if not published(marker(root,t,p)):raise ValueError(f"Missing completed {t}/{p}")
        locks.enter_context(file_lock(lock))
        state=contract(cfg,root,t,m)
        old=published(mark)
        if old and old["inputs"]==state and not force:return
        mark.unlink(missing_ok=True)
        scorefile="predictions.tsv.gz"
        if m.startswith("universal."):
            from universal_blocks import fit, prepare
            mode=m.split(".")[-1]
            if ".prepare." in m:
                outputs=prepare(cfg,root,t,mode,out);scorefile=None
            else:
                outputs=fit(cfg,root,t,mode,out)
        elif m=="prepare":
            execute(cfg,root,"regions",t,SimpleNamespace(replace=force));scorefile=None
            outputs=["scores.npy","ids.tsv","features.tsv","metadata.json","region_counts.tsv","sum_check.tsv"]
        elif m=="architecture":
            from report2 import architecture
            architecture(cfg,root,t)
            (project_path(root,f"qc/{t}.regional_architecture.tsv")).replace(out/"architecture.tsv")
            outputs=["architecture.tsv"];scorefile=None
        elif m=="twins":
            from report2 import twin_report
            twin_report(cfg,root,t,out);outputs=["support.tsv",f"{t}.twins.html"];scorefile=None
        elif m=="robustness":
            from driver2 import evaluate
            import pandas as pd
            rows=[]
            for geometry in ("consensus","global_prs","regional","unit_block"):
                loc=out/geometry
                run_abm(cfg,root,t,geometry,loc)
                performance,_,_=evaluate(cfg,root,t,loc/"abm.predictions.tsv.gz",loc/"evaluation")
                rows.append(performance.assign(geometry=geometry))
            from core2 import atomic_table
            atomic_table(pd.concat(rows,ignore_index=True),out/"geometry.tsv");scorefile=None;outputs=["geometry.tsv"]
        elif m=="abm":
            run_abm(cfg,root,t);scorefile="abm.predictions.tsv.gz";outputs=[scorefile,"model.json"]
        elif m=="blocknet.twins":
            from blocknet_twins import fit
            fit(cfg,root,t);outputs=[scorefile,"model.json","model.joblib"]
        elif m.startswith("direct."):
            from direct_blocknet import fit
            fit(cfg,root,t,m);outputs=[scorefile,"model.json","model.joblib"]
            if m=="direct.blocknet":outputs.append("embedding.npy")
        elif m=="block.stack":
            from block_stack import fit
            fit(cfg,root,t,out);scorefile="predictions.tsv.gz";outputs=[scorefile,"model.json"]
        elif m=="disco.weights":
            from disco_weights import fit
            fit(cfg,root,t,out);scorefile=None;outputs=["weights.tsv.gz","model.json"]
        elif m.endswith(".disco") or m.startswith("blocknet.disco"):
            from disco_blocks import fit
            fit(cfg,root,t,m,out);outputs=[scorefile,"model.json","model.joblib"]
        else:
            c=read_table(project_path(root,"tmp/cohort/participants.tsv.gz"));x,_=regional_inputs(root,t,c)
            common=(c,read_table(project_path(root,f"tmp/cohort/{t}.development.tsv.gz")),
                    read_table(trait_path(root,t,"baseline.predictions.tsv.gz")),x,
                    read_table(trait_path(root,t,"regions/features.tsv")),out)
            kw=dict(trait=t,seed=cfg["cohort"]["seed"],threads=cfg["runtime"]["threads"])
            if m=="abm.deep":
                from block_deep import fit_block_deep
                fit_block_deep(*common,config=cfg["block_deep"],**kw)
            else:
                from block_models import fit_block_models
                opt=dict(cfg["block_models"],glmnet=m=="block.glmnet",lightgbm=m=="block.lightgbm")
                fit_block_models(*common,config=opt,rscript=cfg["runtime"]["rscript"],**kw)
            outputs=[scorefile,"model.json"]
        if state!=contract(cfg,root,t,m):raise RuntimeError(f"{t}/{m}: upstream inputs changed during fit")
        atomic_json(dict(format=RECEIPT_FORMAT,method=m,trait=t,predictions=scorefile,inputs=state,
                         outputs={f:stamp(out/f) for f in outputs}),mark)

def write_index(root):
    """Lightweight method inventory for the main analysis, links to published completion records."""
    import html
    import pandas as pd
    from core2 import atomic_table
    from output_layout import TRAITS
    with file_lock(project_path(root,"tmp/locks/index.lock")):
        rows=[]
        for trait in sorted(TRAITS):
            for mark in sorted((root/trait).glob("*/complete.json")):
                if not mark.is_file():
                    continue
                try:
                    doc=json.loads(mark.read_text())
                except FileNotFoundError:
                    continue
                path=mark.parent.relative_to(root)
                rows.append(dict(trait=trait,method=doc["method"],directory=str(path),
                                 predictions=str(path/doc["predictions"]) if doc.get("predictions") else ""))
        atomic_table(pd.DataFrame(rows,columns=["trait","method","directory","predictions"]),root/"results/methods.tsv")
        path=root/"reports/index.html"
        path.parent.mkdir(parents=True,exist_ok=True)
        lines=['<!doctype html><meta charset="utf-8"><title>GRID methods</title><h1>Completed GRID methods</h1><ul>']
        for row in rows:
            name=html.escape(row["trait"]+" / "+row["method"])
            link=html.escape(row["directory"]+"/complete.json",quote=True)
            lines.append(f'<li><a href="../{link}">{name}</a></li>')
        tmp=path.with_suffix(".partial")
        tmp.write_text("\n".join(lines+["</ul>"]))
        tmp.replace(path)
