#!/usr/bin/env python3
"""Plan independently selectable methods and incremental report updates."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from driver2 import configure,parser
from core2 import atomic_json
from methods import DEEP,PARENTS,contract,marker,published,selection
from output_layout import project_path
from workflow_layout import new_run_dir
from workflow_mode import resolve_mode


from cluster import active_ids,submit,wait_jobs,placement


def schedule(cfg,args):
    root=Path(cfg["analysis_root"]);traits=cfg["traits"]
    suffix=("" if root.name=="main" else "."+root.name)
    prefix="grid"+suffix
    report=args.paper
    full_run = not report and args.module == "run" and (args.method is None or args.method == selection("all"))
    names=[]
    if not report:
        if args.module=="check":
            from report2 import check_inputs
            from grid_workflow import validate_inputs
            validate_inputs(cfg,"check");return check_inputs(cfg)
        if args.module=="regions":
            names=["prepare"]
        elif args.module=="prepare":
            names=(["universal.prepare.auto","universal.prepare.fixed"] if cfg["regions"]["layout"]=="consensus" else ["prepare"])
        elif args.module=="run":
            if args.method is not None:
                names=args.method
            else:
                names=selection("all")
        else:
            names=selection(args.module)
        if cfg["regions"]["layout"]!="consensus":
            if args.method is not None and args.method!=selection("all") and any(m.startswith("universal.") for m in names):
                raise ValueError("Universal models need one shared GRCh38 BED.")
            names=[n for n in names if not n.startswith("universal.")]
        if args.without_deep or not cfg["block_deep"]["enabled"]:
            names=[n for n in names if n not in DEEP|{"blocknet.twins","block.stack"}]
    jobs={};hpc=resolve_mode(args.submit)=="hpc"
    def add(t,m,parents=(),patterns=(),force=False):
        p=("eval"+suffix if m=="evaluate" else "paper"+suffix if m=="publish" else prefix)
        name=f"{prefix}.eval.{t}" if m=="evaluate" else (f"{p}.{t}.{m}" if t else f"{p}.{m}")
        if name in jobs:return name
        gpu=m in DEEP and cfg["block_deep"]["device"]=="cuda"
        cores=min(1 if m in {"twins","publish","architecture","evaluate","block.stack","disco.weights"} else 3 if m in {"universal.auto","universal.fixed","block.glmnet","block.glmnet.disco"} else
                  args.gpu_cores if gpu else args.job_cores,args.job_cores)
        jobs[name]=dict(trait=t,method=m,parents=list(parents),patterns=list(patterns),force=force,cores=cores,gpu=gpu)
        return name
    def model(t,m):
        name=f"{prefix}.{t}.{m}"
        if name in jobs:return name
        parents=[parent for p in PARENTS[m] if (parent:=model(t,p))]
        force=args.replace and m in names
        if not parents and not force:
            receipt=published(marker(root,t,m))
            if receipt and receipt["inputs"]==contract(cfg,root,t,m):
                if not hpc or args.plan_only or not active_ids(name):
                    return None
        return add(t,m,parents,force=force)
    if not report:
        for t in traits:
            for m in names:
                model(t,m)
        if not jobs and not full_run:
            print("All requested methods are current.");return
    if report or full_run:
        evaluations = []
        for t in traits:
            producers = [name for name, job in jobs.items() if job["trait"] == t]
            evaluations.append(add(t, "evaluate", parents=producers, patterns=[f"{prefix}.{t}.*"]))
        add(None, "publish", parents=evaluations)
    taskroot=project_path(root,"jobs")
    if report:
        from output_layout import upstream_root
        taskroot=upstream_root(root)/"paper"/"jobs"/suffix.lstrip(".")
    batch=new_run_dir(taskroot)
    logroot=project_path(root,"grid.paper/logs" if report else "logs")/batch.name
    logroot.mkdir(parents=True,exist_ok=True)
    atomic_json(cfg,batch/"config.json");atomic_json(jobs,batch/"plan.json")
    # A full run records one DAG; report logs belong to paper/logs in both invocation modes.
    paper_logs = project_path(root, "grid.paper/logs") / batch.name if report else project_path(root, "grid.paper/logs") / ("grid-" + batch.name)
    if report or full_run:
        paper_logs.mkdir(parents=True, exist_ok=True)
    active={}
    for name,j in jobs.items():
        worker_logs = paper_logs if j["method"] in {"evaluate", "publish"} else logroot
        args2=[sys.executable,str(Path(__file__).resolve()),"worker",str(batch),name]
        env=os.environ.copy();threads=min(j["cores"],cfg["runtime"]["threads"])
        for k in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS","WORKFLOW_CORES"):env[k]=str(threads)
        script=batch/(name+".sh");script.write_text("#!/usr/bin/env bash\nset -euo pipefail\n"+
            "\n".join("export "+k+"="+env[k] for k in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS","WORKFLOW_CORES"))+"\nexec "+shlex.join(args2)+"\n")
        script.chmod(0o700)
        if args.plan_only:
            print(name,"after",", ".join(j["parents"] or j["patterns"]) or "—");continue
        if not hpc:
            with (worker_logs/(name+".out")).open("w") as stdout,(worker_logs/(name+".err")).open("w") as stderr:
                rc=subprocess.run([str(script)],stdout=stdout,stderr=stderr,env=env).returncode
            if rc:raise RuntimeError(f"{name} failed: {worker_logs}/{name}.err")
            continue
        # Bind dependencies to this run's job IDs, including any GRID fits already
        # active when a separate paper command is submitted.
        parents = {active[p] for p in j["parents"] if p in active}
        for pattern in j["patterns"]:
            parents.update(active_ids(pattern))
        q=(args.gpu_queue or cfg["runtime"]["gpu_queue"]) if j["gpu"] else args.queue
        slots=placement(cfg,root,j["trait"],j["method"],j["cores"],q) if not j["gpu"] else j["cores"]
        signature={k:v for k,v in cfg.items() if k not in {"runtime","traits","evaluation","worker_python"}}
        active[name]=submit(script,"grid.paper" if j["method"] == "publish" else name,slots,q,j["gpu"],worker_logs,sorted(parents),force=j["force"],
                            state_dir=project_path(root,"tmp/locks/submit"),signature=signature,
                            repeat=j["method"] in {"evaluate","publish"})
    if args.wait and hpc and not args.plan_only:wait_jobs(active)


def work(batch,name):
    from threadpoolctl import threadpool_limits
    cfg=json.loads((batch/"config.json").read_text());j=json.loads((batch/"plan.json").read_text())[name]
    root=Path(cfg["analysis_root"]);t=j["trait"];m=j["method"]
    cfg["runtime"]["threads"]=min(j["cores"],cfg["runtime"]["threads"])
    os.environ["GRID_RSCRIPT"]=cfg["runtime"]["rscript"]
    with threadpool_limits(limits=cfg["runtime"]["threads"]):
        if m in ("evaluate","publish"):
            from methods_report import evaluate
            if t:
                return evaluate(cfg,root,t)
            from paper_figures import publish_manuscript
            return publish_manuscript(cfg,root)
        from methods import run,write_index
        run(cfg,root,t,m,j["force"]);write_index(root)

def main(*,paper=False):
    cli=parser("grid",paper=paper)
    args=cli.parse_args()
    if args.method is not None and args.module!="run":
        cli.error("Use --method with ./grid.sh run")
    if args.job_cores is None and resolve_mode(args.submit)=="hpc":
        from driver2 import shared_cpu_budget
        args.job_cores=min(20,shared_cpu_budget(args.submit))
    schedule(configure(args),args)


if __name__=="__main__":
    if sys.argv[1:2]==["predict"]:
        from universal_blocks import predict_cli
        predict_cli(sys.argv[2:])
    elif sys.argv[1:2]==["worker"]:work(Path(sys.argv[2]),sys.argv[3])
    else:main()
