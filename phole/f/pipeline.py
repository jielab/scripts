#!/usr/bin/env python3
"""Run PHOLE's summary audit and signed locus-by-assay workflow.

Planning imports only the Python standard library and never scans GWAS inputs.
The environment is detected automatically. Both local and HPC execution run the
same generated worker.cmd; HPC uses the shared LSF runner. --plan-only writes only.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from contextlib import ExitStack
from pathlib import Path

F_DIR=Path(__file__).resolve().parent
PROJECT=F_DIR.parent
SCRIPTS=PROJECT.parent
sys.path.insert(0,str(SCRIPTS/"0f"))
from workflow_mode import resolve_mode
from workflow_layout import new_run_dir
DATASETS={"ukb":"prot","ckb":"ckb","decode":"decode","scot":"scot"}
ORDER=("cross-check","rebuild-summaries","rebuild-plan","fingerprint-plan","manifest","specificity-check","fingerprint-check","prepare","edge-check","fit","discover","validate","plot")
ALL=tuple(x for x in ORDER if x not in {"rebuild-plan","fingerprint-plan"})
SCHEMA="phole-current-pipeline"


def atomic_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix="."+path.name,dir=path.parent)
    try:
        with os.fdopen(fd,"w") as stream:
            json.dump(value,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write("\n")
        os.replace(name,path)
    finally:Path(name).unlink(missing_ok=True)


def stamp():return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream,delimiter="\t"))


def file_hash(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bool_text(value):
    text=str(value).upper()
    if text not in {"TRUE","FALSE"}:raise argparse.ArgumentTypeError("requires TRUE or FALSE")
    return text


def arguments(argv=None):
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Modules (comma separated):
  all         cross-check -> rebuild-summaries -> manifest -> specificity-check -> fingerprint-check -> prepare -> edge-check -> fit -> discover -> validate -> plot
  cross-check Audit loci, coordinates, assay specificity, complete fingerprint dossiers and ML readiness.
  rebuild-plan  Write four exact gwas/format.sh mplot commands; do not execute them.
  fingerprint-plan  Freeze different-target regional follow-up from existing summaries; no complete GWAS reads.
  rebuild-summaries  Regenerate summaries and source metadata from complete GWAS, then rerun cross-check.
  manifest    Discover unique complete summary paths, or check --assay-manifest.
  specificity-check  Query external cis leads in full GWAS, retaining nonsignificant and missing results separately.
  fingerprint-check  Read complete GWAS regions for candidate targets AND alternative proteins, with same-cohort and intended-target controls.
  prepare     Extract fixed, allele-harmonized locus rows from complete summaries.
  edge-check  Test exact-variant direction and conservative cohort replication.
  fit,discover,validate,plot  Run selected signed-graph stages using existing inputs.

For copy-ready Qiming commands, see the Examples section of bash phole.sh -h.

Complete primary candidate enumeration and signed summary inputs are required.
Selected .sig.txt or MAGMA statistics are not signed effect matrices.
""")
    p.add_argument("modules",nargs="?",default="all")
    p.add_argument("--submit",choices=["auto","local","hpc"],default=None,help="Execution environment [auto]; explicit overrides remain available")
    p.add_argument("--plan-only",dest="plan",action="store_true",help="Write worker/plan files only; no input discovery, analysis or submission")
    p.add_argument("--wait",nargs="?",const="TRUE",type=bool_text,default=os.environ.get("SUBMIT_WAIT","FALSE"))
    p.add_argument("--python",default=os.environ.get("PHOLE_PYTHON",os.environ.get("WORKFLOW_PYTHON",sys.executable)),help="Scientific interpreter [PHOLE_PYTHON / WORKFLOW_PYTHON]")
    p.add_argument("--queue",default=os.environ.get("PHOLE_QUEUE",os.environ.get("HPC_CPU_QUEUE","33")))
    p.add_argument("--cores",type=int,default=40)
    p.add_argument("--memory-gb",type=int,default=32)
    p.add_argument("--hours",type=float,default=120)
    p.add_argument("--dir-gwas-root",type=Path,default=Path(os.environ.get("DIR_GWAS","gwas")),help="Four-cohort full-GWAS root (UKB subfolder is prot)")
    p.add_argument("--dir-cross-root",type=Path,help="Four-cohort mplot root [completed PHOLE rebuild, otherwise --dir-gwas-root]")
    p.add_argument("--dir-out",type=Path,default=Path(os.environ.get("ANALYSIS_ROOT","."))/"phole",help="Persistent output root")
    p.add_argument("--protein-map",type=Path,help="Protein mapping [<gwas-root>/prot.names.tsv]; no bundled fallback")
    p.add_argument("--annotation",type=Path,help="Protein annotation [<gwas-root>/prot.annot.tsv]")
    p.add_argument("--sentinels",type=Path,default=Path(os.environ.get("OS_SNP_LIST",str(Path(os.environ.get("DIR_UKB",str(SCRIPTS.parent/"phe/ukb")))/"common/snp.lst"))),help="Build sentinel table [shared OS_SNP_LIST / DIR_UKB/common/snp.lst]")
    p.add_argument("--primary",choices=DATASETS,default="ukb")
    p.add_argument("--category",default="common")
    p.add_argument("--assay-manifest",type=Path,help="Explicit complete-summary manifest takes precedence over path discovery")
    p.add_argument("--candidate-loci",type=Path)
    p.add_argument("--dataset-builds",default="ukb=38,ckb=38,decode=38,scot=38",help="User-confirmed actual files [ukb=38,ckb=38,decode=38,scot=38]; override for different inputs")
    p.add_argument("--dataset-ancestries",default="ukb=EUR,ckb=EAS,decode=EUR,scot=EUR",help="Cohort ancestries [ukb=EUR,ckb=EAS,decode=EUR,scot=EUR]; pass an empty string for unknown")
    p.add_argument("--refresh-summaries",action="store_true",help="After the initial cross-check, rebuild four mplot summaries and audit their output")
    p.add_argument("--allow-scan",action="store_true",default=True,help="Allow sequential full-summary reads if tabix is unavailable (default)")
    p.add_argument("--no-allow-scan",dest="allow_scan",action="store_false",help="Require indexed extraction; disable sequential full-file reads")
    p.add_argument("--tabix",default="tabix")
    p.add_argument("--cis-flank",type=int,default=500000)
    p.add_argument("--coordinate-window",choices=["body","point"],default="body",help="Use corroborated gene bodies or supplied points before adding cis-flank")
    p.add_argument("--map-build",choices=["37","38"],default="38",help="Declared build of optional protein-map CHR/POS")
    p.add_argument("--trans-rich",type=int,default=5)
    p.add_argument("--discovery-protein-list",type=Path)
    p.add_argument("--target-protein-list",type=Path,help="Optional canonical PROTs retained as matrix targets in every cohort")
    p.add_argument("--target-selection-scope",choices=["prespecified","primary_only","external_screened","unknown"],default="unknown")
    p.add_argument("--specificity-queries",type=Path,help="External-selected own-cis query list [specificity/specificity-queries.tsv]")
    p.add_argument("--skip-specificity-window-scan",action="store_true",help="Query external sentinel variants only; do not inspect full target cis intervals")
    p.add_argument("--with-fingerprint-gwas",action="store_true",help="After cross-check, discover complete GWAS paths and run different-protein regional follow-up; does not rewrite .sig.txt")
    p.add_argument("--fingerprint-pairs",type=Path,help="Explicit frozen different-target hypotheses; otherwise select current T1/T2 candidates")
    p.add_argument("--fingerprint-selection",choices=["core-six","all-priority"],default="core-six",help="Summary-based nomination rule; core-six means current T1/T2 tier, not a forced case count")
    p.add_argument("--fingerprint-region-flank",type=int,default=250000,help="Full-GWAS window on each side of union summary leads; gene cis windows separately use --cis-flank")
    p.add_argument("--locus-bed",type=Path)
    p.add_argument("--locus-build",choices=["37","38"])
    p.add_argument("--locus-bed-coordinates",choices=["bed0","one-based"],default="bed0")
    p.add_argument("--folds",type=int,default=5)
    p.add_argument("--rank-grid",default="0,2,4,8")
    p.add_argument("--iterations",type=int,default=30)
    p.add_argument("--ridge",type=float,default=5.)
    p.add_argument("--seed",type=int,default=20261105)
    p.add_argument("--min-targets",type=int,default=10)
    p.add_argument("--validation-min-targets",type=int,help="Optional coverage threshold for external validation [--min-targets]")
    p.add_argument("--min-null",type=int,default=30)
    p.add_argument("--max-candidates",type=int,default=50)
    p.add_argument("--replace",type=bool_text,default="TRUE")
    p.add_argument("--run-config",type=Path,help=argparse.SUPPRESS)
    p.add_argument("--worker",action="store_true",help=argparse.SUPPRESS)
    a=p.parse_args(argv)
    if a.run_config:
        config=json.load(sys.stdin) if a.run_config == Path("-") else json.loads(a.run_config.read_text())
        if config.get("schema")!=SCHEMA:raise ValueError("Unrecognized PHOLE plan schema")
        required=set(vars(a)) - {"run_config", "worker"}
        saved=config.get("arguments", {})
        if not required.issubset(saved):
            raise ValueError("Incomplete current PHOLE plan; generate a fresh --plan-only request")
        a=argparse.Namespace(**saved)
        for key in ("dir_gwas_root","dir_cross_root","dir_out","protein_map","annotation","sentinels","assay_manifest","candidate_loci","discovery_protein_list","target_protein_list","specificity_queries","fingerprint_pairs","locus_bed","run_dir"):
            if getattr(a,key,None) is not None:setattr(a,key,Path(getattr(a,key)))
        a.worker=True;a.plan=False
    else:
        if a.dir_cross_root is None:
            summaries=a.dir_out.expanduser().resolve()/"summaries"
            receipt=control_path(a,"rebuild-summaries/summary.json")
            report=json.loads(receipt.read_text()) if receipt.is_file() else {}
            complete=(report.get("summary_root")==str(summaries) and set(report.get("datasets",()))==set(DATASETS)
                      and report.get("tasks",0)>0 and not report.get("failures") and not report.get("limit")
                      and all((summaries/folder/"mplot").is_dir() for folder in DATASETS.values()))
            a.dir_cross_root=summaries if complete else a.dir_gwas_root
        a.protein_map=a.protein_map or a.dir_gwas_root/"prot.names.tsv"
        a.annotation=a.annotation or a.dir_gwas_root/"prot.annot.tsv"
        for key,value in vars(a).copy().items():
            if isinstance(value,Path):setattr(a,key,value.expanduser().resolve())
    a.submit=resolve_mode(a.submit)
    requested=ALL if a.modules=="all" else tuple(x.strip() for x in a.modules.split(",") if x.strip())
    if not requested or set(requested)-set(ORDER):p.error("Unsupported module list: "+a.modules)
    if a.with_fingerprint_gwas:
        if "cross-check" not in requested:p.error("--with-fingerprint-gwas requires cross-check in the module list")
        requested=(*requested,"manifest","fingerprint-check")
    a.selected=[x for x in ORDER if x in requested]
    if a.refresh_summaries and "rebuild-summaries" in a.selected:p.error("Use rebuild-summaries or --refresh-summaries, not both")
    if a.refresh_summaries and "cross-check" not in a.selected:p.error("--refresh-summaries requires all or cross-check in the module list")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",a.queue):p.error("Invalid queue name")
    if not a.category or a.category in {".",".."} or any(x in a.category for x in "/\\\r\n\t"):p.error("Invalid category name")
    if a.cores<1 or a.memory_gb<1 or not 0<a.hours<=10000:p.error("cores/memory-gb/hours must be positive")
    if a.folds<3 or a.iterations<1 or not math.isfinite(a.ridge) or a.ridge<0 or a.min_targets<2 or a.min_null<1 or a.max_candidates<1 or a.cis_flank<0:
        p.error("Invalid statistical or window setting")
    if a.validation_min_targets is not None and a.validation_min_targets<2:p.error("validation-min-targets must be at least 2")
    if a.fingerprint_region_flank<0:p.error("fingerprint-region-flank must be nonnegative")
    try:
        ranks=[int(x) for x in a.rank_grid.split(",")]
        if any(x<0 for x in ranks) or not ranks:raise ValueError()
    except ValueError:p.error("rank-grid needs nonnegative comma-separated integers")
    return a


def control_path(a,name):
    return a.dir_out/"logs"/name


def paths(a):
    return {"cross":a.dir_out/"cross-check","matrix":a.dir_out/"matrix","fit":a.dir_out/"fit",
            "edge-check":a.dir_out/"edge-check",
            "specificity":a.dir_out/"specificity","specificity-check":a.dir_out/"qc/specificity-check",
            "fingerprints":a.dir_out/"qc/fingerprints",
            "fingerprint-check":a.dir_out/"qc/fingerprint-check",
            "fingerprint-plan":a.dir_out/"qc/fingerprint-plan",
            "readiness":a.dir_out/"readiness",
            "discover":a.dir_out/"discover","validate":a.dir_out/"validate","plot":a.dir_out/"plot",
            "manifest":a.assay_manifest or control_path(a,"assay-manifest.tsv"),
            "candidate":a.candidate_loci or a.dir_out/"cross-check/candidate-loci.tsv",
            "build":control_path(a,"build-receipt.tsv"),"summaries":a.dir_out/"summaries"}


def command(a,module,scope="full_summary_candidates",cross_root=None):
    p=paths(a);base=[a.python,str(F_DIR)]
    def py(name,*args):return [base[0],str(F_DIR/name),*map(str,args)]
    if module=="cross-check":
        cmd=py("cross_check.py","--gwas-root",cross_root or a.dir_cross_root,"--protein-map",a.protein_map,
               "--annotation",a.annotation,"--sentinels",a.sentinels,"--out",p["cross"],"--primary",a.primary,"--cis-flank",a.cis_flank,"--trans-rich",a.trans_rich,
               "--coordinate-window",a.coordinate_window,"--map-build",a.map_build)
        if a.dataset_builds:cmd += ["--builds",a.dataset_builds]
        return cmd
    if module=="rebuild-summaries":
        return py("rebuild_summaries.py","--gwas-root",a.dir_gwas_root,"--cross-dir",p["cross"],
                  "--out",control_path(a,"rebuild-summaries"),"--summary-root",p["summaries"],"--cores",a.cores,"--cis-flank",a.cis_flank,
                  "--dataset-builds",a.dataset_builds)
    if module=="manifest":
        cmd=py("prepare.py","--export-manifest",p["manifest"],"--gwas-root",a.dir_gwas_root,"--protein-map",p["cross"]/"prot.names.resolved.tsv",
               "--category",a.category,"--candidate-loci",p["candidate"],"--dataset-ancestries",a.dataset_ancestries)
        if a.dataset_builds:cmd += ["--dataset-builds",a.dataset_builds]
        return cmd
    if module=="prepare":
        cmd=py("prepare.py","--assay-manifest",p["manifest"],"--candidate-loci",p["candidate"],"--protein-annotation",p["cross"]/"coordinates-used.tsv",
               "--cross-check-dir",p["cross"],"--out",p["matrix"],"--state-dir",control_path(a,"matrix-state"),"--discovery",a.primary,"--cis-flank",a.cis_flank,
               "--tabix",a.tabix,"--replace",a.replace,"--panel-completeness",scope)
        if a.allow_scan:cmd += ["--allow-scan"]
        for opt in ("discovery_protein_list","target_protein_list","locus_bed","locus_build"):
            if getattr(a,opt):cmd += ["--"+opt.replace("_","-"),str(getattr(a,opt))]
        if a.locus_bed:cmd += ["--locus-bed-coordinates",a.locus_bed_coordinates]
        if a.target_protein_list:cmd += ["--target-selection-scope",a.target_selection_scope]
        return cmd
    if module=="specificity":
        return py("assay_specificity.py","--cross-dir",p["cross"],"--protein-map",a.protein_map,"--out",p["specificity"],"--primary","all")
    if module=="fingerprints":
        return py("fingerprint_dossiers.py","--cross-dir",p["cross"],"--protein-map",a.protein_map,
                  "--specificity-dir",p["specificity"],"--out",p["fingerprints"],"--primary","all")
    if module=="readiness":
        return py("ml_readiness.py","--cross-dir",p["cross"],"--out",p["readiness"],"--primary",a.primary)
    if module=="review-report":
        return py("review_report.py","--cross-dir",p["cross"],"--specificity-dir",p["specificity"],"--readiness-dir",p["readiness"],
                  "--fingerprints-dir",p["fingerprints"])
    if module=="specificity-check":
        cmd=py("specificity_check.py","--assay-manifest",p["manifest"],"--queries",a.specificity_queries or p["specificity"]/"specificity-queries.tsv",
               "--out",p["specificity-check"],"--state-dir",control_path(a,"specificity-state"),"--tabix",a.tabix,"--primary",a.primary)
        if a.allow_scan:cmd += ["--allow-scan"]
        if not a.skip_specificity_window_scan:
            cmd += ["--cis-annotation",str(p["cross"]/"coordinates-used.tsv"),"--scan-cis-window","--cis-flank",str(a.cis_flank)]
        return cmd
    if module in {"fingerprint-plan","fingerprint-check"}:
        cmd=py("fingerprint_check.py","--cross-dir",p["cross"],"--fingerprints-dir",p["fingerprints"],
               "--out",p[module],"--selection",a.fingerprint_selection,
               "--region-flank",a.fingerprint_region_flank,"--cis-flank",a.cis_flank,"--tabix",a.tabix)
        if a.fingerprint_pairs:cmd += ["--pairs",str(a.fingerprint_pairs)]
        if module=="fingerprint-plan":cmd += ["--plan-only"]
        else:
            cmd += ["--assay-manifest",str(p["manifest"])]
            if a.allow_scan:cmd += ["--allow-scan"]
        return cmd
    if module=="fingerprint-report":
        return py("fingerprint_report.py","--check-dir",p["fingerprint-check"],"--cross-dir",p["cross"])
    if module=="fit":return py("fit.py","--matrix-dir",p["matrix"]/a.primary,"--out",p["fit"],"--folds",a.folds,
        "--rank-grid",a.rank_grid,"--iterations",a.iterations,"--ridge",a.ridge,"--seed",a.seed,"--replace",a.replace)
    if module=="edge-check":return py("edge_check.py","--matrix-root",p["matrix"],"--out",p["edge-check"],"--primary",a.primary)
    if module=="discover":return py("discover.py","--matrix-dir",p["matrix"]/a.primary,"--model-dir",p["fit"],
        "--out",p["discover"],"--min-targets",a.min_targets,"--max-candidates",a.max_candidates)
    if module=="validate":return py("validate.py","--matrix-root",p["matrix"],"--model-dir",p["fit"],"--discovery-dir",p["discover"],
        "--out",p["validate"],"--primary",a.primary,"--datasets",",".join(d for d in DATASETS if d!=a.primary),
        "--min-targets",a.validation_min_targets or a.min_targets,"--min-null",a.min_null,"--seed",a.seed,"--replace",a.replace)
    if module=="plot":return py("plot.py","--matrix-dir",p["matrix"]/a.primary,"--model-dir",p["fit"],"--discovery-dir",p["discover"],
        "--validation-dir",p["validate"],"--cross-check-dir",p["cross"],"--out",p["plot"],"--seed",a.seed,"--replace",a.replace)
    raise ValueError("Module has no standalone Python command: "+module)


def rebuild_commands(a):
    # Read original full GWAS, never display-thinned or COJO inputs.
    return {"full-gwas-summaries":command(a,"rebuild-summaries")}


def write_rebuild_plan(a):
    folder=a.dir_out/"logs";folder.mkdir(parents=True,exist_ok=True)
    cmds=rebuild_commands(a)
    target=folder/"rebuild-summaries.cmd"
    target.write_text("#!/usr/bin/env bash\nset -euo pipefail\n# Requires cross-check/cis-bed; exact full-GWAS roots, no path guessing.\n"+
                      "\n".join(shlex.join(cmd) for cmd in cmds.values())+"\n")
    target.chmod(0o755)
    atomic_json(folder/"rebuild-summaries.plan.json",{"schema":"phole-summary-rebuild-v2","commands":cmds,"executed":False,
                                  "summary_output_root":str(paths(a)["summaries"]),"build":"38","cis_flank":a.cis_flank})
    print(f"Summary rebuild commands: {target}",flush=True)
    return cmds


def serialize(a):
    return {k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items() if k not in {"run_config","worker"}}


def write_job_plan(a):
    folder=new_run_dir(a.dir_out/"logs")
    a.run_dir=folder
    receipt=folder/"worker.plan.json";job=folder/"worker.cmd"
    planned={m:command(a,m,"full_summary_candidates")
             for m in a.selected if m!="rebuild-plan" and not (m=="manifest" and a.assay_manifest)}
    if "manifest" in a.selected and a.assay_manifest:
        planned["manifest"]={"action":"check provided complete-summary manifest; do not overwrite", "path":str(a.assay_manifest)}
    execution_order=[]
    for module in a.selected:
        execution_order.append(module)
        if module=="cross-check" and a.refresh_summaries:
            for d,cmd in rebuild_commands(a).items():
                planned["rebuild-"+d]=cmd;execution_order.append("rebuild-"+d)
            planned["cross-check-refreshed"]=command(a,"cross-check",cross_root=paths(a)["summaries"])
            execution_order.append("cross-check-refreshed")
        if module=="rebuild-summaries":
            planned["cross-check-refreshed"]=command(a,"cross-check",cross_root=paths(a)["summaries"])
            execution_order.append("cross-check-refreshed")
        if module=="rebuild-summaries" or (module=="cross-check" and "rebuild-summaries" not in a.selected):
            for companion in ("specificity","fingerprints","readiness","review-report"):
                planned[companion]=command(a,companion);execution_order.append(companion)
        if module in {"fingerprint-plan","fingerprint-check"}:
            report_name=module+"-report"
            planned[report_name]=[a.python,str(F_DIR/"fingerprint_report.py"),"--check-dir",str(paths(a)[module]),
                                  "--cross-dir",str(paths(a)["cross"])]
            execution_order.append(report_name)
    if "rebuild-plan" in a.selected:planned["rebuild-plan"]={"action":"write commands only","commands":rebuild_commands(a)}
    plan={"schema":SCHEMA,"generated":stamp(),"arguments":serialize(a),"commands":planned,"execution_order":execution_order,
                         "worker":str(job),"planning_scope":"parameters only; no GWAS or numerical imports",
                         "primary_panel_gate":"complete current full-summary exports"}
    atomic_json(receipt,plan)
    queue,slots=a.queue,a.cores
    if a.submit=="hpc":
        from qiming_resources import cpu_slots
        queue,slots=cpu_slots(queue,slots,a.memory_gb,a.hours)
    minutes=math.ceil(a.hours*60)
    directives=["#!/usr/bin/env bash","#BSUB -J phole",f"#BSUB -q {queue}",f"#BSUB -n {slots}",
                f"#BSUB -W {minutes//60}:{minutes%60:02d}",
                f'#BSUB -R "span[hosts=1] rusage[mem={a.memory_gb}GB/host]"',
                f"#BSUB -o {shlex.quote(str(folder/'worker.out'))}",f"#BSUB -e {shlex.quote(str(folder/'worker.err'))}",
                "set -euo pipefail","export LC_ALL=C PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1",
                "export CUDA_VISIBLE_DEVICES='' WORKFLOW_DEVICE=cpu",
                f"export OPENBLAS_NUM_THREADS={a.cores} OMP_NUM_THREADS={a.cores} MKL_NUM_THREADS={a.cores}"]
    for key in ("DIR_WORK","DIR_DAT","DIR_BIG","DIR_DATA","DIR_TMP","DIR_SCRATCH","DIR_GWAS","DIR_FILES","DIR_SOFTWARE","ANALYSIS_ROOT","WORKFLOW_PYTHON","PHOLE_PYTHON","GWAS_RUNTIME_ENV"):
        if key in os.environ:directives.append("export "+key+"="+shlex.quote(os.environ[key]))
    directives += ["export submit="+a.submit+" SUBMIT="+a.submit,
                   "cd -- "+shlex.quote(str(SCRIPTS)),
                   "exec "+shlex.join([a.python,str(F_DIR/"pipeline.py"),"--run-config","-","--worker"]) + " <<'PHOLE_CONFIG'",
                   json.dumps(plan, ensure_ascii=True, allow_nan=False), "PHOLE_CONFIG"]
    fd, temporary = tempfile.mkstemp(prefix=".worker.", dir=folder)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write("\n".join(directives)+"\n")
        Path(temporary).chmod(0o750)
        os.replace(temporary, job)
    finally:
        Path(temporary).unlink(missing_ok=True)
    print(f"Worker plan: {receipt}\nWorker command: {job}",flush=True)
    return job


def require_file(path,meaning):
    if not Path(path).is_file():raise ValueError(f"Missing {meaning}: {path}")


def panel_scope(a):
    summary=paths(a)["cross"]/"summary.json"
    info=json.loads(summary.read_text()) if summary.is_file() else {}
    primary=info.get("datasets",{}).get(a.primary,{})
    if primary.get("complete_full_summaries") is True:return "full_summary_candidates"
    raise ValueError(f"{a.primary}: complete primary full-summary exports are not certified by {summary}. "
                     "Run cross-check,rebuild-plan and rebuild summaries, or use all --refresh-summaries.")


def write_build_receipt(a):
    p=paths(a);summary=p["cross"]/"summary.json"
    if not summary.is_file():return False
    data=json.loads(summary.read_text());rows=[]
    for d,info in data.get("datasets",{}).items():
        b=str(info.get("build","")).replace("GRCh","")
        if d in DATASETS and b in {"37","38"}:
            rows.append({"DATASET":d,"BUILD":b,"SOURCE":info.get("build_source","cross_check"),
                         "EVIDENCE":str(summary),"EVIDENCE_SHA256":file_hash(summary),
                         "SCOPE":"dataset inference; does not certify every assay"})
    if not rows:return False
    with p["build"].open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        writer.writeheader();writer.writerows(rows)
    return True


def check_manifest(a):
    manifest=paths(a)["manifest"];require_file(manifest,"complete-summary assay manifest")
    rows=read_tsv(manifest)
    if not rows:raise ValueError(f"No complete GWAS files found in {manifest}. Inspect {manifest}.issues.tsv and provide real full summaries; no matrix was created.")
    required={"DATASET","ASSAY","PROT","PATH","BUILD","INPUT_SCOPE"}
    if not required.issubset(rows[0]):raise ValueError("Assay manifest lacks complete-summary contract columns: "+", ".join(sorted(required)))
    missing=[];primary=0
    for row in rows:
        d="ukb" if row["DATASET"].lower()=="prot" else row["DATASET"].lower()
        primary += d==a.primary
        if row["INPUT_SCOPE"].lower()!="full":raise ValueError("Selected sig/MAGMA inputs are not complete signed summaries")
        path=Path(row["PATH"]).expanduser();path=path if path.is_absolute() else manifest.parent/path
        if not path.is_file():missing.append(str(path))
    if missing:raise ValueError("Complete GWAS files are missing; no matrix was created. First paths: "+"; ".join(missing[:10]))
    if not primary:raise ValueError("Assay manifest contains no complete primary GWAS files: "+a.primary)
    print(f"Complete-summary paths available: {len(rows)} assays ({primary} primary). Header/allele validation runs in prepare.",flush=True)
    return rows


class Runner:
    def __init__(self,a):
        self.a=a;self.logs=getattr(a,"run_dir",None) or a.dir_out/"logs";self.logs.mkdir(parents=True,exist_ok=True)
        self.receipt=self.logs/"pipeline-state.json"
        self.state={"schema":SCHEMA,"started":stamp(),"arguments":serialize(a),"stages":[],"status":"running"}
        atomic_json(self.receipt,self.state)

    def run(self,name,argv):
        logfile=self.logs/(name+".log")
        event={"module":name,"command":argv,"log":str(logfile),"started":stamp(),"status":"running"}
        self.state["stages"].append(event);atomic_json(self.receipt,self.state)
        print("["+name+"] "+shlex.join(argv),flush=True)
        env={**os.environ,"PYTHONUNBUFFERED":"1"}
        for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):env.setdefault(key,str(self.a.cores))
        try:
            with logfile.open("a",encoding="utf-8") as stream:
                stream.write("\n["+stamp()+"] "+shlex.join(argv)+"\n");stream.flush()
                process=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding="utf-8",errors="replace",env=env)
                try:
                    for line in process.stdout:sys.stdout.write(line);sys.stdout.flush();stream.write(line);stream.flush()
                    code=process.wait()
                except BaseException:
                    process.terminate();process.wait();raise
            event.update(finished=stamp(),returncode=code,status="complete" if code==0 else "failed")
            atomic_json(self.receipt,self.state)
            if code:raise RuntimeError(f"{name} failed with exit {code}; see {logfile}. Resume selected stages after correcting the input.")
        except BaseException:
            event.update(finished=stamp(),status="failed");self.state["status"]="failed";atomic_json(self.receipt,self.state);raise

    def finish(self):
        self.state.update(status="complete",finished=stamp());atomic_json(self.receipt,self.state)


def execute(a):
    if a.submit=="hpc" and not os.environ.get("LSB_JOBID"):
        raise ValueError("HPC worker requires a real LSF allocation. Submit the generated .cmd; use --submit local for an explicit local run.")
    if not shutil.which(a.python):raise ValueError("Scientific Python is unavailable: "+a.python+"; set --python or PHOLE_PYTHON")
    folder = a.dir_out / "logs"
    folder.mkdir(parents=True, exist_ok=True)
    locks=[folder/"pipeline.lock"]
    with ExitStack() as stack:
        for path in locks:
            lock=stack.enter_context(path.open("a"))
            try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError(f"PHOLE output in use: {a.dir_out}") from exc
        execute_stages(a)


def execute_stages(a):
    p=paths(a);runner=Runner(a)
    try:
        for module in a.selected:
            if module=="cross-check":
                require_file(a.protein_map,"protein mapping");require_file(a.annotation,"protein annotation")
                runner.run(module,command(a,module))
                if a.refresh_summaries:
                    summary=json.loads((p["cross"]/"summary.json").read_text())
                    incompatible=[d for d in DATASETS if str(summary.get("datasets",{}).get(d,{}).get("build"))!="38"]
                    if incompatible:
                        raise ValueError("Summary rebuild commands require a GRCh38 receipt for every cohort; unresolved/incompatible: "+
                                         ", ".join(incompatible)+". Declare verified --dataset-builds or standardize upstream; this workflow does not perform liftover.")
                    for d,cmd in write_rebuild_plan(a).items():
                        runner.run("rebuild-"+d,cmd)
                    runner.run("cross-check-refreshed",command(a,module,cross_root=paths(a)["summaries"]))
                if "rebuild-summaries" not in a.selected:
                    for companion in ("specificity","fingerprints","readiness","review-report"):
                        runner.run(companion,command(a,companion))
            elif module=="rebuild-summaries":
                require_file(p["cross"]/"coordinates-used.tsv","prior cross-check coordinates")
                runner.run(module,command(a,module))
                runner.run("cross-check-refreshed",command(a,"cross-check",cross_root=paths(a)["summaries"]))
                for companion in ("specificity","fingerprints","readiness","review-report"):
                    runner.run(companion,command(a,companion))
            elif module=="rebuild-plan":write_rebuild_plan(a)
            elif module=="manifest":
                if a.assay_manifest:print("Using explicit complete-summary manifest: "+str(a.assay_manifest),flush=True)
                else:
                    cmd=command(a,module)
                    map_file=p["cross"]/"prot.names.resolved.tsv"
                    if not map_file.is_file():cmd[cmd.index("--protein-map")+1]=str(a.protein_map)
                    require_file(Path(cmd[cmd.index("--protein-map")+1]),"protein mapping")
                    if not p["candidate"].is_file():
                        index=cmd.index("--candidate-loci");del cmd[index:index+2]
                    # Summary coordinates are useful diagnostics, but cannot
                    # certify the build of a separate complete GWAS file.
                    # Only fresh per-file evidence or an explicit declaration
                    # is accepted by prepare.export_manifest.
                    write_build_receipt(a)
                    runner.run(module,cmd)
                check_manifest(a)
            elif module in {"fingerprint-plan","fingerprint-check"}:
                require_file(p["cross"]/"coordinates-used.tsv","cross-check coordinates")
                require_file(p["fingerprints"]/"assay-evidence.tsv","summary fingerprint hypotheses")
                if module=="fingerprint-check":require_file(p["manifest"],"complete-summary assay manifest")
                runner.run(module,command(a,module))
                runner.run(module+"-report",[a.python,str(F_DIR/"fingerprint_report.py"),
                    "--check-dir",str(p[module]),"--cross-dir",str(p["cross"])])
            elif module=="prepare":
                scope=panel_scope(a);check_manifest(a)
                require_file(p["candidate"],"candidate locus table")
                runner.state["panel_completeness"]=scope
                cmd=command(a,module,scope)
                if not (p["cross"]/"coordinates-used.tsv").is_file():
                    cmd[cmd.index("--protein-annotation")+1]=str(a.annotation)
                require_file(Path(cmd[cmd.index("--protein-annotation")+1]),"effective protein annotation")
                summary=p["cross"]/"summary.json"
                if summary.is_file():cmd += ["--panel-receipt",str(summary)]
                runner.run(module,cmd)
            else:runner.run(module,command(a,module))
        runner.finish()
    except BaseException:
        runner.state.update(status="failed",finished=stamp());atomic_json(runner.receipt,runner.state);raise
    print("PHOLE requested stages completed: "+", ".join(a.selected)+"; "+str(a.dir_out),flush=True)


def main(argv=None):
    a=arguments(argv)
    if a.selected==["rebuild-plan"] and not a.plan:
        write_rebuild_plan(a)
        return
    if a.worker:
        execute(a)
        return
    folder=a.dir_out/"logs"
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/"submission.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        job=write_job_plan(a)
        submit_job(a,job,job.parent)


def submit_job(a, job, folder):
    if a.plan:
        if "rebuild-plan" in a.selected:write_rebuild_plan(a)
        print("Plan only: no input files scanned, analysis run, or cluster job submitted.")
        return
    if a.submit=="hpc" and not shutil.which("bsub"):raise ValueError("bsub is unavailable; worker plan was saved. Use --submit local to run locally, or --plan-only to review.")
    submit=[sys.executable,str(SCRIPTS/"0f/submit_jobs.py"),"--submit",a.submit,"--jobs-N","1","--journal",str(folder/"submission.jsonl")]
    if a.wait=="TRUE":submit += ["--wait"]
    subprocess.run([*submit,str(job)],check=True)
    return

if __name__=="__main__":
    try:main()
    except (ValueError,RuntimeError,OSError,subprocess.SubprocessError) as error:
        raise SystemExit("ERROR: "+str(error)) from error
