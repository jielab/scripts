#!/usr/bin/env python3
"""Cohort-oriented orchestration for outcome-free diet scores and paper analyses."""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from score_registry import SCORE_ENGINES, score_names

ROOT=Path(__file__).resolve().parents[1]
COHORTS=("ukb","nhanes","chns","charls")
CAPABILITIES={"ukb":["maha","dose"],"nhanes":["maha","dose"],"chns":["maha"],"charls":["charls"]}

def parser():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples (run from scripts/diet):
  ./install.sh --check
  ./1.disco.sh all --dry-run
  ./diet.sh all --core --dry-run

  # UKB: score first, then phenotype, association and Bayes, in order.
  ./1.disco.sh ukb --score dose,dose11,dose9,dose9half
  ./diet.sh ukb --study dose
  ./diet.sh assoc --cohort ukb --study dose
  ./diet.sh bayes --cohort ukb --study dose
  ./diet.sh final --study dose
  ./diet.sh web

  # HPC: inspect a job plan, then submit (LSF; do not run analysis on login nodes).
  ./1.disco.sh ukb --score dose,dose11,dose9,dose9half --submit hpc --jobs 8 --job-cores 8 --job-memory-gb 160 --queue 33 --plan-only
  ./1.disco.sh ukb --score dose,dose11,dose9,dose9half --submit hpc --jobs 8 --job-cores 8 --job-memory-gb 160 --queue 33
  # Wait for scoring to succeed before submitting downstream work:
  ./diet.sh all --cohort ukb --study dose --core --submit hpc --jobs 8 --job-cores 8 --job-memory-gb 160 --queue 33

  # Once MAHA dictionaries and all cohort inputs are available:
  ./1.disco.sh all
  ./diet.sh all --core
  ./diet.sh all
  ./diet.sh status

Paths default to /work/sph-huangj/{files,analysis}/diet on HPC and
/mnt/d/{files,analysis}/diet locally through ../0f/0os.sh.
Override with --dir-files, --dir-output, --dir-ukb, --foods or --config.
MAHA NHANES (--dir-nhanes) requires raw/<cycle>/; DOSE NHANES
(--dir-nhanes-dose) requires cdc/ and FPED/MPED inputs.
R/Python selection follows le8/f/runtime_env.sh and ../0f/0os.sh.
HPC defaults: LE8_RSCRIPT / LE8_PYTHON; local: system R and discovered le8 Python.
Override with DIET_RSCRIPT / DIET_PYTHON, R_BIN / PYTHON_BIN, or
--rscript (--r-bin) / --python (--python-bin).
Job names: 1.disco.sh -> diet.disco; diet.sh -> diet.
With no arguments, either entry shows help. Score first, then explicitly select all/modules.
Common scheduler options: --submit auto|local|hpc, --plan-only,
--job-cores N, --job-memory-gb N, --job-hours N, --queue NAME.
`--dry-run` prints scientific stages; `--plan-only` writes a scheduler plan.
`all` never computes scores. `all --core` omits DOSE FoodWAS/ML/mechanisms/
signatures; full `all` requires their extra inputs. CHARLS needs a declared mapping.
""")
    p.add_argument("--entry",choices=["disco","diet"],default="diet",help=argparse.SUPPRESS)
    p.add_argument("action",nargs="?",default="help")
    p.add_argument("cohort_arg",nargs="?",choices=COHORTS+("all",),help="Optional cohort after an analysis module")
    p.add_argument("--cohort",choices=COHORTS+("all",),default="all")
    p.add_argument("--study",choices=["all","maha","dose"],default="all")
    p.add_argument("--score",help="Scoring entry only: all or comma-separated diet.XXX.sum names (XXX)")
    p.add_argument("--config",type=Path,default=ROOT/"config/diet.json")
    p.add_argument("--dir-output",type=Path)
    p.add_argument("--dir-files",type=Path)
    p.add_argument("--foods",type=Path)
    p.add_argument("--dir-ukb",type=Path)
    p.add_argument("--dir-nhanes",type=Path)
    p.add_argument("--dir-nhanes-dose",type=Path,help="DOSE NHANES root containing cdc/ and FPED/MPED inputs")
    p.add_argument("--dir-chns",type=Path)
    p.add_argument("--dir-charls",type=Path)
    p.add_argument("--charls-input",type=Path)
    p.add_argument("--charls-phe",type=Path)
    p.add_argument("--charls-config",type=Path)
    p.add_argument("--rscript","--r-bin",default=os.getenv("DIET_RSCRIPT","Rscript"))
    p.add_argument("--python","--python-bin",default=os.getenv("DIET_PYTHON",sys.executable))
    p.add_argument("--jobs",type=int,default=4)
    p.add_argument("--replace",choices=["TRUE","FALSE"],default="TRUE")
    p.add_argument("--steps",help="Comma-separated original DOSE steps, in registered execution order")
    p.add_argument("--list-steps",action="store_true")
    p.add_argument("--dry-run",action="store_true")
    p.add_argument("--core",action="store_true",help="Omit DOSE derivation and molecular follow-up from all")
    p.add_argument("--collect-only",action="store_true",help="Collect existing paper figures without rerendering composites")
    p.add_argument("--extension",action="append",default=[],help="Run a named configured extension after Bayes")
    return p

def validate_args(p, a):
    if a.entry == "disco" and a.action != "help":
        modules = a.action.split(",")
        if a.action != "all" and any(c not in COHORTS for c in modules):
            p.error("1.disco.sh modules: all or comma-separated ukb,nhanes,chns,charls")
        if a.cohort_arg:
            p.error("Use a single comma-separated cohort list for 1.disco.sh")
        if a.score is not None and a.study != "all":
            p.error("Use --score for scoring; do not combine it with a restricted --study")
    elif a.score is not None:
        p.error("--score belongs to 1.disco.sh; diet.sh analyses use --study")
    try:
        score_names(a.score)
    except ValueError as exc:
        p.error(str(exc))

def csv_steps():
    with (ROOT/"config/dose_steps.tsv").open() as handle:
        return list(csv.DictReader(handle,delimiter="\t"))

def path_setting(cli,env_key,config,default):
    value=cli or os.getenv(env_key) or config or default
    return str(Path(os.path.expandvars(str(value))).expanduser().resolve())

class Workflow:
    def __init__(self,args):
        self.a=args
        self.scores=score_names(args.score)
        self.cfg=json.loads(args.config.read_text())
        # scripts/0f records the DIET_ prefix when it freezes a job plan. Retain
        # the native engines' scientific overrides in that recorded namespace.
        frozen=json.loads(os.getenv("DIET_ENGINE_ENV","{}"))
        if not isinstance(frozen,dict) or any(not isinstance(k,str) or not isinstance(v,str) or
              not k.startswith(("DOSE_","MAHA_","CHNS_","NHANES_")) for k,v in frozen.items()):
            raise ValueError("Invalid DIET_ENGINE_ENV snapshot")
        os.environ.update(frozen)
        self.env=os.environ.copy()
        defaults=self.cfg.get("paths",{})
        # Python direct calls use the same single OS table as shell/R entries.
        helpers=Path(os.getenv("DIET_HELPERS",str(ROOT.parent/"0f")))
        if str(helpers) not in sys.path: sys.path.append(str(helpers))
        from workflow_paths import resolved_table, expand_path_roots
        table=resolved_table()
        defaults=expand_path_roots(defaults)
        dat=table["DIR_DAT"]
        files=path_setting(args.dir_files,"DIET_FILES",defaults.get("files"),str(Path(table["DIR_FILES"])/"diet"))
        out=path_setting(args.dir_output,"DIET_OUTDIR",defaults.get("output"),str(Path(table["ANALYSIS_ROOT"])/"diet"))
        self.out=Path(out)
        self.env.update(DIET_ROOT=str(ROOT),DIET_OUTDIR=out,DIET_OUT=out,DIET_FILES=files,
            DIET_FOODS=path_setting(args.foods,"DIET_FOODS",defaults.get("foods"),str(Path(files)/"foods.xlsx")),
            DIET_CACHE=str(self.out/".cache"),DIET_SCORES=str(self.out/"disco"),
            DIET_HELPERS=os.getenv("DIET_HELPERS",str(ROOT.parent/"0f")),
            DIET_UKB=path_setting(args.dir_ukb,"DIET_UKB",defaults.get("ukb"),table["DIR_UKB"]),
            DIET_NHANES=path_setting(args.dir_nhanes,"DIET_NHANES",defaults.get("nhanes"),str(Path(dat)/"nhanes")),
            DIET_CHNS=path_setting(args.dir_chns,"DIET_CHNS",defaults.get("chns"),str(Path(dat)/"chns")),
            DIET_CHARLS=path_setting(args.dir_charls,"DIET_CHARLS",defaults.get("charls"),str(Path(dat)/"charls")),
            DIET_THREADS=str(args.jobs),DIET_REPLACE=args.replace,
            DOSE_PYTHON=args.python,RETICULATE_PYTHON=args.python,PYTHONDONTWRITEBYTECODE="1",
            OMP_NUM_THREADS=str(args.jobs),OPENBLAS_NUM_THREADS=str(args.jobs),MKL_NUM_THREADS=str(args.jobs))
        for arg,key in [(args.charls_input,"DIET_CHARLS_INPUT"),(args.charls_phe,"DIET_CHARLS_PHE"),(args.charls_config,"DIET_CHARLS_CONFIG")]:
            if arg: self.env[key]=str(arg.resolve())
        # Preserve explicit original engine settings; central paths are canonical.
        self.env.update(DOSE_OUTDIR=str(self.out/"dose"),DOSE_PROJECT=str(self.out/"dose"),
            DOSE_SCRIPT_DIR=str(ROOT/"f/dose"),DOSE_UKB=self.env["DIET_UKB"],
            DOSE_RDATA_DIR=str(self.out/".cache/ukb/dose/Rdata"),
            DOSE_FILES_DIR=files,DOSE_HELPER_DIR=self.env["DIET_HELPERS"],
            NHANES_DATA=path_setting(args.dir_nhanes_dose,"NHANES_DATA",defaults.get("nhanes_dose"),
                str(args.dir_nhanes) if args.dir_nhanes else files),DOSE_THREADS=str(args.jobs),DOSE_REPLACE=args.replace,
            MAHA_HELPER_DIR=self.env["DIET_HELPERS"],MAHA_OUTDIR=str(self.out/"maha"),MAHA_JOBS=str(args.jobs),
            CHNS_PROJECT_DIR=self.env["DIET_CHNS"],
            DIET_UKB_RDATA=str(self.out/".cache/ukb/maha/Rdata"))
        # DOSE engine accepts a workbook filename via its directory helper. The
        # authoritative name is foods.xlsx; avoid copying it to a second cache.
        foods=Path(self.env["DIET_FOODS"])
        if foods.name!="foods.xlsx": raise ValueError("--foods must point to the unified file named foods.xlsx")
        self.env["DOSE_FILES_DIR"]=str(foods.parent)
        self.env["DOSE_MRI_DICTIONARY"]=os.getenv("DOSE_MRI_DICTIONARY",str(Path(files)/"MRI_analysis_dictionary.csv"))
        self.status=[]

    def score_skip(self, cohort, scores, reason):
        lines=[f"SKIP {cohort}: {score}: {reason}" for score in scores]
        for line in lines: print(line, flush=True)
        if not self.a.dry_run:
            log=self.out/"logs"/f"disco.{cohort}.skip.log"
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("a", encoding="utf-8") as handle:
                handle.write(datetime.now(timezone.utc).isoformat(timespec="seconds")+"\n"+"\n".join(lines)+"\n")

    def run(self,label,command,extra=None):
        env={**self.env,**(extra or {})}
        print(f"[{label}] {shlex.join([str(v) for v in command])}",flush=True)
        if self.a.dry_run: return
        first=command[0]
        if not (Path(first).is_file() or shutil.which(first)): raise RuntimeError(f"Interpreter not found: {first}. See install.sh or --rscript/--python.")
        logdir=self.out/"logs";logdir.mkdir(parents=True,exist_ok=True)
        start=datetime.now(timezone.utc).isoformat(timespec="seconds")
        log=logdir/(label.replace("/",".")+".log")
        with log.open("w",encoding="utf-8") as handle:
            handle.write(f"{start}\n{shlex.join([str(v) for v in command])}\n");handle.flush()
            with subprocess.Popen(command,env=env,cwd=self.out,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1) as proc:
                for line in proc.stdout:
                    print(line,end="",flush=True);handle.write(line)
                code=proc.wait()
        self.status.append(dict(stage=label,started=start,finished=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                                status="complete" if code==0 else "failed",exit_code=code,log=str(log.relative_to(self.out))))
        self.save_status()
        if code: raise RuntimeError(f"{label} failed (exit {code}); see {log}")

    def save_status(self):
        from results import rows,write_rows
        path=self.out/"logs/stages.tsv"
        now={x["stage"] for x in self.status}
        old=[x for x in rows(path) if x.get("stage") not in now]
        write_rows(path,old+self.status)

    def r(self,label,file,*args,extra=None):
        self.run(label,[self.a.rscript,str(ROOT/"f"/file),*args],extra)

    def py(self,label,file,*args,extra=None):
        self.run(label,[self.a.python,str(ROOT/"f"/file),*map(str,args)],extra)

    def part_env(self,cohort,study):
        return dict(DIET_COHORT=cohort,DIET_DISCO_FILE=str(self.out/"disco"/f"{cohort}.disco.rds"),
            DIET_DISCO_PART=str(self.out/".cache/disco"/f"{cohort}.{study}.rds"))

    def selected(self):
        if self.a.entry == "disco":
            # Explicitly listing every cohort has exactly the same behavior as all.
            candidates=list(COHORTS) if self.a.action=="all" else list(dict.fromkeys(self.a.action.split(",")))
            selected=[]
            for c in candidates:
                if c=="charls" and not Path(self.env.get("DIET_CHARLS_INPUT",str(Path(self.env["DIET_CHARLS"])/"food_components.tsv"))).is_file():
                    self.score_skip(c,self.scores,"measured-component input is not configured")
                    continue
                if not self.profiles(c):
                    self.score_skip(c,self.scores,"no supported scoring engine for this selection")
                    continue
                selected.append(c)
            return selected
        cohort=self.a.cohort_arg or self.a.cohort
        if self.a.action in COHORTS: cohort=self.a.action
        candidates=list(COHORTS) if cohort=="all" else [cohort]
        if cohort=="all":
            charls=Path(self.env.get("DIET_CHARLS_INPUT",str(Path(self.env["DIET_CHARLS"])/"food_components.tsv")))
            if not charls.exists():
                candidates.remove("charls")
                print("CHARLS: no measured-component file configured; module is available when its input mapping is supplied.",flush=True)
            candidates=[c for c in candidates if self.profiles(c)]
        elif not self.profiles(cohort):
            raise ValueError(f"{cohort} has no {self.a.study} profile. See cohort coverage in README.md.")
        return candidates

    def profiles(self,cohort):
        if self.a.entry == "disco":
            return [engine for engine,names in SCORE_ENGINES[cohort].items()
                    if self.a.study in ("all",engine) and set(names).intersection(self.scores)]
        if cohort=="charls":
            if self.a.study!="all": return []
            return ["charls"]
        return [s for s in CAPABILITIES[cohort] if self.a.study in ("all",s)]

    def disco(self,cohort):
        if not self.profiles(cohort): raise ValueError(f"No selected scoring profile for {cohort}")
        for study in self.profiles(cohort):
            extra=self.part_env(cohort,study)
            if study=="maha":
                if cohort=="ukb": self.r("disco/ukb.maha","ukb_maha_score.R",extra=extra)
                else: self.r(f"disco/{cohort}.maha","maha/run.R",cohort,"disco",extra=extra)
            elif study=="dose": self.r(f"disco/{cohort}.dose","dose/adapter.R",extra={**extra,"DOSE_ACTION":f"disco_{cohort}"})
            else: self.r("disco/charls","charls.R","disco",extra=extra)
        self.r(f"disco/{cohort}.publish","common.R","disco",cohort,
               extra=dict(DIET_SCORE_SELECTION=",".join(self.scores),
                          DIET_SCORE_PROFILES=",".join(self.profiles(cohort))))

    def phe(self,cohort):
        if not self.profiles(cohort): raise ValueError(f"No selected phenotype profile for {cohort}")
        for study in self.profiles(cohort):
            extra=self.part_env(cohort,study)
            if study=="maha": self.r(f"phe/{cohort}.maha","maha/run.R",cohort,"phe",extra=extra)
            elif study=="dose": self.r(f"phe/{cohort}.dose","dose/adapter.R",extra={**extra,"DOSE_ACTION":f"cohort_{cohort}"})
            else: self.r("phe/charls","charls.R","phe",extra=extra)
        self.r(f"phe/{cohort}.publish","common.R","phe",cohort)

    def assoc(self,cohort):
        for study in self.profiles(cohort):
            extra=self.part_env(cohort,study)
            if study=="maha": self.r(f"assoc/{cohort}.maha","maha/run.R",cohort,"assoc",extra=extra)
            elif study=="dose":
                if self.a.steps:
                    if cohort!="ukb": raise ValueError("--steps for assoc is currently UKB DOSE only")
                    allowed={r["step"] for r in csv_steps() if r["module"]=="score" and r["step"]!="build"}
                    selected=self.a.steps.split(",")
                    if not set(selected)<=allowed: raise ValueError("Unknown DOSE association step")
                    extra["DOSE_ASSOC_STEPS"]=self.a.steps
                self.r(f"assoc/{cohort}.dose","dose/adapter.R",extra={**extra,"DOSE_ACTION":f"assoc_{cohort}"})
            else: self.r("assoc/charls","charls.R","assoc",extra=extra)
        self.py("assoc/collect","results.py","collect","--output",self.out)

    def bayes(self,cohorts):
        for cohort in cohorts:
            if "maha" in self.profiles(cohort):
                self.r(f"bayes/{cohort}.maha","maha/run.R",cohort,"bayes",extra=self.part_env(cohort,"maha"))
        studies=["maha","dose","charls"] if self.a.study=="all" else [self.a.study]
        settings=self.cfg.get("bayes",{})
        self.py("bayes/collect","results.py","bayes","--output",self.out,"--cohorts",",".join(cohorts),
            "--studies",",".join(studies),"--priors",",".join(map(str,settings.get("prior_95_ratios",[1.05,1.1,1.2]))),
            "--margins",",".join(map(str,settings.get("equivalence_ratios",[1.05,1.075,1.1]))))

    def dose_steps(self,module,selection=None):
        table=[r for r in csv_steps() if r["module"]==module]
        if selection:
            asked=selection.split(","); known={r["step"] for r in table}
            unknown=set(asked)-known
            if unknown: raise ValueError(f"Unknown {module} steps: {', '.join(sorted(unknown))}")
            table=[r for r in table if r["step"] in asked]
        else: table=[r for r in table if r["mode"]=="default"]
        for row in table:
            extra=dict(DOSE_MODULE=module,DOSE_STEP=row["step"],DOSE_ACTION="step")
            selected_mri=self.out/"dose/ukb/05mechanism/04mri/inputs/mri_selected.tsv"
            if row["step"]!="mri_prepare": extra["DOSE_MRI_RAW"]=self.env.get("DOSE_MRI_RAW",str(selected_mri))
            if row["interpreter"]=="python": self.py(f"dose/{module}.{row['step']}",f"dose/{row['file']}",extra=extra)
            else: self.r(f"dose/{module}.{row['step']}","dose/adapter.R",extra=extra)

    def final(self):
        studies=["maha","dose"] if self.a.study=="all" else [self.a.study]
        scope=self.a.cohort_arg or self.a.cohort
        if not self.a.collect_only and "maha" in studies:
            if scope in ("all","ukb","nhanes","chns"):
                self.r("final/maha","maha/run.R",scope,"final")
        if not self.a.collect_only and "dose" in studies and scope in ("all","ukb"):
            # These composite panels are downstream of their finished inputs.
            nh=self.out/"dose/nhanes/NHANES_M2_results.csv"
            uk=self.out/"dose/ukb/04score/pattern/pattern_continuous.csv"
            comparison_inputs=[nh,uk,self.out/"dose/nhanes/NHANES_analysis.rds",
                self.out/"dose/ukb/00data/analysis_config.rds",self.out/"dose/ukb/00data/dat.base.imp.rds",
                self.out/"dose/ukb/04score/dose_score.rds"]
            if self.a.dry_run or all(p.is_file() for p in comparison_inputs): self.dose_steps("score","stroke_comparison")
            # A directory alone may be left by an interrupted signature stage.
            # Auto-refresh the native summary only after its producers completed.
            from results import rows
            completed={r.get("stage") for r in rows(self.out/"logs/stages.tsv") if r.get("status")=="complete"}
            producers={"dose/mechanism.protein","dose/mechanism.metabolite_plot",
                "dose/mechanism.inflammation_plot","dose/mechanism.pathway_plot",
                "dose/mechanism.sem_plot","dose/signature.protein_cox","dose/signature.metabolite_plot"}
            if self.a.dry_run or producers<=completed: self.dose_steps("mechanism","summary")
        if not self.a.collect_only and "dose" in studies:
            self.py("final/dose","dose/final.py","--output",self.out)
        self.py("final/collect","results.py","final","--output",self.out,"--studies",",".join(studies))

    def web(self):
        self.py("web/collect","results.py","collect","--output",self.out)
        self.py("web/render","web.py","--output",self.out)

    def extensions(self):
        configured=self.cfg.get("extensions",{})
        for name in self.a.extension:
            if name not in configured: raise ValueError(f"Unknown extension {name}; register its command in config/diet.json")
            spec=configured[name]
            command=[str(v).replace("{diet_root}",str(ROOT)).replace("{output}",str(self.out)) for v in spec["command"]]
            self.run("extension/"+name,command)

def main(argv=None):
    p=parser();a=p.parse_args(argv)
    validate_args(p,a)
    if a.action=="help":p.print_help();return 0
    allowed=set(COHORTS)|{"all","assoc","bayes","foodwas","ml","mechanism","signature","final","web","extensions","status"}
    if a.entry!="disco" and a.action not in allowed:p.error(f"Unknown module: {a.action}")
    if a.jobs<1:p.error("--jobs must be positive")
    if a.action in {"foodwas","ml","mechanism","signature"}:
        if a.study not in {"all","dose"} or (a.cohort_arg or a.cohort) not in {"all","ukb"}:
            p.error(f"{a.action} is a UKB DOSE module; use --study dose --cohort ukb")
    if a.list_steps:
        for row in csv_steps():
            if a.action in (row["module"],"all"):
                print("\t".join(row.values()))
        return 0
    w=Workflow(a);cohorts=w.selected()
    if not cohorts:
        if a.entry=="disco": return 0
        raise ValueError("No cohorts selected")
    if a.action=="status":
        from results import rows
        for record in rows(w.out/"logs/stages.tsv"): print("\t".join(str(record.get(k,"")) for k in ("stage","status","finished","log")))
        return 0
    if a.entry=="disco":
        for cohort in cohorts:w.disco(cohort)
    elif a.action in COHORTS:
        for cohort in cohorts:w.phe(cohort)
    elif a.action=="assoc":
        for cohort in cohorts:w.assoc(cohort)
    elif a.action=="bayes":w.bayes(cohorts);w.extensions()
    elif a.action in {"foodwas","ml","mechanism","signature"}:w.dose_steps(a.action,a.steps)
    elif a.action=="final":w.final()
    elif a.action=="web":w.web()
    elif a.action=="extensions":w.extensions()
    else:
        for cohort in cohorts:w.phe(cohort)
        include_dose="ukb" in cohorts and a.study in {"all","dose"}
        if include_dose and not a.core:
            w.dose_steps("foodwas");w.dose_steps("ml")
        for cohort in cohorts:w.assoc(cohort)
        if include_dose and not a.core:
            w.dose_steps("mechanism");w.dose_steps("signature")
        w.bayes(cohorts);w.extensions();w.final();w.web()
    return 0
if __name__=="__main__":
    try:sys.exit(main())
    except (OSError,ValueError,RuntimeError) as exc:sys.exit(f"ERROR: {exc}")
