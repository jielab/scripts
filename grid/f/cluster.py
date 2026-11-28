"""LSF dispatch for GRID jobs. Commands specify CPU placement; the queue supplies hard limits."""
import os
import json
from pathlib import Path
import shlex
import subprocess

from core2 import atomic_json
from file_lock import file_lock
from lsf_submit import require_submitted,submit_attempt


def active_ids(name):
    p=subprocess.run(["bjobs","-noheader","-o","jobid","-J",name],capture_output=True,text=True)
    if p.returncode and not any(s in (p.stderr+p.stdout).lower() for s in ("not found","no unfinished job","no job")):
        raise RuntimeError(p.stderr or p.stdout)
    return [int(s.strip()) for s in p.stdout.splitlines() if s.strip().isdigit()]


def submit(script,name,cores,queue,gpu,logroot,parents=(),patterns=(),*,force=False,state_dir=None,signature=None,repeat=False,optional_parents=()):
    path=Path(script)
    state=Path(state_dir) if state_dir is not None else path.parent.parent
    receipt=path.parent/(name+".submission.json")
    intent=state/(name+".active.json")
    with file_lock(state/(name+".submit.lock")):
        existing=active_ids(name)
        previous=[]
        if repeat:
            head,tail=name.rsplit(".",1)
            previous=sorted(set(existing+active_ids(head+".*."+tail)))
            if previous:
                index=2
                while active_ids(f"{head}.{index}.{tail}"):index+=1
                name=f"{head}.{index}.{tail}"
        elif existing:
            saved=json.loads(intent.read_text()) if intent.is_file() else {}
            if force or len(existing)!=1 or saved.get("job_id")!=existing[0] or saved.get("signature")!=signature:
                raise ValueError(f"{name} is active with a different or unverified request: {existing}")
            return existing[0]
        expression=['done("'+p+'")' for p in patterns]
        expression.extend(f"done({jid})" for jid in parents)
        expression.extend(f"ended({jid})" for jid in (*optional_parents, *previous))
        cmd=["bsub","-J",name,"-n",str(cores),"-R","span[hosts=1]","-oo",str(logroot/(name+".%J.out")),"-eo",str(logroot/(name+".%J.err"))]
        queue=queue if queue and queue!="auto" else os.environ.get("HPC_CPU_QUEUE")
        if queue:cmd+=["-q",queue]
        if gpu:cmd+=["-gpu","num=1"]
        if expression:cmd+=["-w"," && ".join(dict.fromkeys(expression))]
        result=submit_attempt(cmd,path=path,persist=lambda v:atomic_json(v,receipt))
        jid=require_submitted(result,receipt)
        if not repeat:atomic_json({"job_id":int(jid),"signature":signature},intent)
        print(f"[SUBMITTED] {name} {jid}",flush=True)
        return int(jid)


def wait_jobs(jobs):
    for name,jid in jobs.items():
        subprocess.run(["bwait","-w",f"ended({jid})"],check=True)
        state=subprocess.check_output(["bjobs","-a","-noheader","-o","stat",str(jid)],text=True).strip()
        if state!="DONE":raise RuntimeError(f"{name} {jid}: {state}")


def numbered(jobs,component,args,worker_python):
    from workflow_layout import new_run_dir
    run=new_run_dir(component/"jobs");logs=component/"logs"/run.name;logs.mkdir(parents=True,exist_ok=True)
    ids={};names={}
    for j in jobs:
        name=f"{args.entry}.{j.key}"
        names[j.key]=name
        resources=j.resources
        cores=resources["cores"]
        env={k:v for k,v in j.env.items() if k!="GRID_MEMORY_CAP_GB"}
        env.update(WORKFLOW_CORES=str(cores),WORKFLOW_PYTHON=worker_python)
        script=run/(name+".sh")
        script.write_text("#!/usr/bin/env bash\nset -euo pipefail\ncd "+shlex.quote(j.cwd)+"\n"+
                          "\n".join("export "+k+"="+shlex.quote(str(v)) for k,v in env.items())+
                          "\nexec "+shlex.join(j.argv)+"\n")
        script.chmod(0o700)
        parents=[k for k in j.depends if k]
        if args.plan_only:
            print(name,"after",", ".join(names[k] for k in parents) or "—");continue
        patterns=[]
        if args.entry=="csx" and parents and all(".chr" in k and k.endswith(".infer") for k in parents):
            patterns=[f"csx.{j.key.split('.')[0]}.chr*.infer"]
        ids[j.key]=submit(script,name,cores,resources["queue"],resources.get("gpus"),logs,
                          [] if patterns else [ids[k] for k in parents],patterns)
    atomic_json({"jobs":names,"submitted":ids},run/"plan.json")
    if args.wait:wait_jobs(ids)

def placement(cfg,root,t,method,cores,queue):
    """On queue 33 reserve a whole node for an estimated >96 GiB fit. Numerical parallelism remains bounded."""
    if queue in (None,"auto"):queue=os.environ.get("HPC_CPU_QUEUE","33")
    if queue!="33" or method not in {"universal.auto","universal.fixed","block.glmnet","block.glmnet.disco","block.lightgbm"}:return cores
    from output_layout import project_path,trait_path
    import pandas as pd
    qc=project_path(root,"qc/cohort.tsv")
    if not qc.is_file():return cores
    counts=pd.read_csv(qc,sep="\t")
    rows=counts.loc[counts.role.isin(["train","tune"]),"n"].sum()
    catalog=trait_path(root,t,"regions/features.tsv")
    if method.startswith("universal."):
        catalog=trait_path(root,t,f"universal.prepare.{method.split('.')[-1]}/features.tsv")
    if catalog.is_file():width=len(pd.read_csv(catalog,sep="\t",usecols=["map"]))
    else:
        from ancestry_regions import read_maps
        maps=read_maps(cfg)
        width=len(maps["shared"][0])*(1 if method.startswith("universal.") else 4) if "shared" in maps else sum(len(b) for b,_ in maps.values())
    if method=="block.glmnet.disco" and cfg["regions"]["layout"]=="consensus":width//=4
    # R's double-precision matrix, copies for fitting, and concurrent alpha paths.
    multiplier=2+4*min(cores,len(cfg["block_models"]["alpha_grid"])) if method!="block.lightgbm" else 3
    gib=rows*width*8*multiplier/2**30
    return 40 if gib>96 else cores
