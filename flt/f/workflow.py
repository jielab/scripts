"""Dispatch the original scientific programs in dependency order."""
import argparse
import os
import json
import sys
import shlex
import subprocess
from pathlib import Path

analysis = {
    "prepare": ["s0_1_flt.R", "s0_2_outcome.R", "s0_3_lifestyle_36vars.R", "s0_4_prepare_data.R"],
    "cox": ["s1_1_cox.R"], "intake": ["s1_2_intake_spearman.R"], "behavior": ["s1_3_otherphe_lm.R"],
    "protein": ["s1_4_protein.R"], "mri": ["s1_5_mri_lm.R"], "uvmr": ["s2_1_uvmr.R"],
    "mvmr": ["s2_2_mvmr.R"], "ldsc": ["s3_1_ldsc_hdl.R"], "placo": ["s3_2_placo.R"],
    "smr": ["s4_1_reform_smr.R", "s4_2_flt_smr.sh"], "phewas": ["s3_3_placo_phewas.R"]
}
paper = {
    "fig1": ["s0_1_flt_goup_plot.R", "s9_fig1.R"], "fig2": ["s9_fig2.R"],
    "fig3": ["s1_2_intake_spearman_plot.R", "s1_3_otherphe_lm_plot.R", "s9_fig3.R"],
    "fig4": ["s1_5_mri_lm.plot.R", "s1_5_mri_lm.plot.py", "s9_fig4.R"],
    "fig5": ["s5_1_fuma_manhattan_plot.R", "s5_2_fuma_gene_plot.R", "s5_3_fuma_enrichment_plot.R",
             "s3_3_placo_phewas_plot.R", "s9_fig5.R"], "tables": ["paper_tables.R"]
}
p = argparse.ArgumentParser()
p.add_argument("entry", choices=["flt.sh","flt.paper.sh"])
p.add_argument("stage")
p.add_argument("--dry-run", action="store_true")
p.add_argument("--dir-output")
p.add_argument("--cores", type=int)
p.add_argument("--plan-only", "--write-cmd", action="store_true")
p.add_argument("--job-memory-gb", type=float, default=float(os.getenv("MEMORY_DEFAULT_GIB","32")))
p.add_argument("--job-hours", type=float, default=float(os.getenv("JOB_HOURS","120")))
p.add_argument("--dir-jobs", type=Path)
p.add_argument("--queue", default=os.getenv("HPC_CPU_QUEUE","33"))
a = p.parse_args()
stages = analysis if a.entry == "flt.sh" else paper
names = list(stages) if a.stage == "run" else [a.stage]
if not all(s in stages for s in names):
    p.error(f"Unknown module {a.stage}")
code = Path(os.environ["FLT_CODE"])
if not a.dry_run and (a.plan_only or os.getenv("WORKFLOW_WORKER") != "1"):
    sys.path.insert(0,os.environ["FLT_SHARED"])
    from workflow_jobs import Job, resolve_mode, run_plan, write_plan
    label = a.entry.removesuffix(".sh") + "." + a.stage
    entry = code.parent / a.entry
    env = {k:v for k,v in os.environ.items() if k.startswith("FLT_")}
    cores = int(os.environ["WORKFLOW_CORES"])
    env.update(WORKFLOW_CORES=str(cores), WORKFLOW_DEVICE="cpu", WORKFLOW_RESUME_POLICY="new_plan")
    job = Job(label, ["bash",str(entry),a.stage], job_name=label, cwd=str(code.parent),env=env,
              resources={"queue":a.queue,"cores":cores,"memory_gib":a.job_memory_gb,"hours":a.job_hours,"gpus":0})
    mode = resolve_mode()
    plan = write_plan([job], a.dir_jobs or Path(env["FLT_OUT"])/"logs", mode=mode, name=label)
    print(json.dumps(run_plan(plan,submit=mode,jobs=1,wait=os.getenv("SUBMIT_WAIT")=="TRUE",dry_run=a.plan_only),indent=2))
    sys.exit(0)
interpreters = {".R": [os.environ["FLT_RSCRIPT"], "--vanilla"], ".sh": ["bash"], ".py": [os.environ["FLT_PYTHON"]]}
for stage in names:
    if stage=="phewas" and not a.dry_run and not (Path(os.environ["FLT_FUMA_DIR"])/"leadSNPs.txt").is_file():
        sys.exit("PheWAS needs the FUMA return for G732 × mdd_2025. Set FLT_FUMA_DIR and resume with ./flt.sh phewas.")
    for script in stages[stage]:
        path = code / script
        cmd = interpreters[path.suffix] + [str(path)]
        print(shlex.join(cmd), flush=True)
        if not a.dry_run:
            Path(os.environ["FLT_PREP"]).mkdir(parents=True, exist_ok=True)
            if stage == "smr":
                Path(os.environ["FLT_SMR_INPUT"]).mkdir(parents=True, exist_ok=True)
            subprocess.run(cmd, check=True)
