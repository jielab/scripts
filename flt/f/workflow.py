"""Plan and run FLT stages with 0os.sh mode detection and LSF job dependencies."""
import argparse
import csv
import fcntl
import json
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

CODE = Path(__file__).resolve().parent
CONFIG = CODE.parent / "config"
sys.path.insert(0, os.getenv("FLT_SHARED", str(CODE.parents[1] / "0f")))
from lsf_submit import submit_graph
from submission_io import atomic_json

ANALYSIS = {
    "prepare": ["s0_1_flt.R", "s0_2_outcome.R", "s0_3_lifestyle_36vars.R", "s0_4_prepare_data.R"],
    "cox": ["s1_1_cox.R"], "cox-sensitivity": ["cox_sensitivity.R"],
    "intake": ["s1_2_intake_spearman.R"], "behavior": ["s1_3_otherphe_lm.R"],
    "protein": ["s1_4_protein.R"], "mri": ["s1_5_mri_lm.R"],
    "uvmr-instruments": ["s2_1_uvmr.R"], "uvmr-forward": ["s2_1_uvmr.R"],
    "uvmr-reverse": ["s2_1_uvmr.R"], "uvmr": ["uvmr_collect.R"],
    "mvmr": ["s2_2_mvmr.R"],
    "ldsc-input": ["s3_1_ldsc_hdl.R"], "ldsc-pairs": ["s3_1_ldsc_hdl.R"],
    "ldsc": ["s3_1_ldsc_hdl.R"], "placo": ["s3_2_placo.R"],
    "coloc": ["coloc.R"], "smr": ["s4_1_reform_smr.R", "s4_2_flt_smr.sh"],
    "phewas": ["s3_3_placo_phewas.R"],
}
PAPER = {
    "fig1": ["s0_1_flt_goup_plot.R", "s9_fig1.R"], "fig2": ["s9_fig2.R"],
    "fig3": ["s9_fig3.R"],
    "fig4": ["s1_5_mri_lm.plot.R", "s1_5_mri_lm.plot.py", "s9_fig4.R"],
    "fig5": [],
    "genetics": ["paper_genetics.R"], "tables": ["paper_tables.R"],
}
INSTALL = {"install": ["prepare_resources.py", "install.R", "check_ready.py"],
           "check": ["install.R", "check_ready.py"]}


def records(name):
    with (CONFIG / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def graph(entry, module, panel, cores, slots):
    jobs = {}

    def add(stage, *, variant="", deps=(), cpus=1, env=None, heavy=False):
        tag = f"{stage}-{variant}" if variant else stage
        name = ("flt-paper-" if entry == "flt.paper.sh" else "flt-") + tag
        if name in jobs or not re.fullmatch(r"[A-Za-z0-9-]+(?:_[A-Za-z0-9-]+)*", name):
            raise ValueError(f"Invalid or duplicate job name {name}")
        jobs[name] = dict(stage=stage, deps=list(deps), cores=min(cpus, cores),
                          slots=slots if heavy else min(cpus, cores), env=env or {})
        return name

    if entry == "install.sh":
        add("install" if module == "run" else "check", cpus=cores, heavy=True)
        return jobs
    if entry == "flt.paper.sh":
        for stage in PAPER if module == "run" else [module]:
            deps = ["flt-paper-fig5"] if stage == "tables" and module == "run" else []
            add(stage, deps=deps, cpus=1)
        return jobs

    outcomes = [row["id"] for row in records("panels.tsv")
                if panel == "extended" or row["panel"] != "extended"]
    full = module == "run"
    cohort = {"cox": 1, "cox-sensitivity": min(cores, 8), "intake": 1, "behavior": 1,
              "protein": 1, "mri": 4}
    prep = add("prepare", cpus=4) if full or module == "prepare" else None
    for stage, threads in cohort.items():
        if full or stage == module:
            add(stage, deps=[prep] if prep else [], cpus=threads)
    uv = None
    if full or module == "uvmr":
        instrument_jobs = [
            add("uvmr-instruments", variant=str(i), cpus=2,
                env={"FLT_SHARD": str(i), "FLT_SHARDS": "6"})
            for i in range(1, 7)
        ]
        forward = [add("uvmr-forward", variant=y, deps=instrument_jobs,
                       cpus=2, env={"FLT_OUTCOME": y}) for y in outcomes]
        reverse = add("uvmr-reverse", deps=instrument_jobs, cpus=min(cores, 8), heavy=True)
        uv = add("uvmr", deps=forward + [reverse])
    if full or module == "mvmr":
        add("mvmr", deps=[uv] if uv else [], cpus=2, heavy=True)
    if full or module == "ldsc":
        inputs = [add("ldsc-input", variant=str(i), cpus=2,
                      env={"FLT_SHARD": str(i), "FLT_SHARDS": "8"}) for i in range(1, 9)]
        fits = [add("ldsc-pairs", variant=y, deps=inputs, cpus=min(cores, 8),
                    env={"FLT_OUTCOME": y}, heavy=True) for y in outcomes]
        add("ldsc", deps=fits)
    placo = []
    if full or module == "placo":
        for row in records("placo-pairs.tsv"):
            a, b = row["flt"], row["brain"]
            name = add("placo", variant=f"{a}-{b}", cpus=min(cores, 8), heavy=True,
                       env={"FLT_PAIR_FLT": a, "FLT_PAIR_BRAIN": b,
                            "FLT_PAIR_BUILD": row["build"], "FLT_PLACO_METHOD": row["method"]})
            if (a, b) == ("G732", "mdd_2025") and row["method"] == "PLACO":
                placo.append(name)
    if full or module == "smr":
        add("smr", cpus=cores, heavy=True)
    for stage in ("coloc", "phewas"):
        if full or module == stage:
            add(stage, deps=placo, cpus=2, heavy=(stage == "coloc"))
    return jobs


def software_environment(args):
    root = Path(args.dir_output or os.environ["FLT_PROJECT"]).expanduser().resolve()
    paper = root / "paper"
    env = {key: value for key, value in os.environ.items()
           if key.startswith("FLT_") and not re.search("TOKEN|PASSWORD|JWT|SECRET|KEY", key)}
    env.update(FLT_PROJECT=str(root), FLT_CONFIG=str(CONFIG), FLT_PANEL=args.panel,
               FLT_OUT=str(root / "main/results"), FLT_REPORTS=str(root / "main/reports"),
               FLT_QC=str(root / "qc"), FLT_TMP=str(root / "tmp"), FLT_CACHE=str(root / "cache"),
               FLT_PLACO_DIR=str(root / "main/results/03placo"), FLT_PLACO_PREFIX="G732-mdd_2025_PLACO",
               FLT_PHEWAS_CACHE=str(root / "cache/phewas"),
               FLT_PREP=str(root / "tmp/prepared"), FLT_COHORT=str(root / "tmp/cohort"), FLT_SMR_INPUT=str(root / "tmp/smr-input"),
               FLT_PAPER=str(paper / "main/results"), FLT_PAPER_REPORTS=str(paper / "main/reports"),
               FLT_PAPER_TMP=str(paper / "tmp"))
    if args.dir_out_phewas_cache:
        env["FLT_PHEWAS_CACHE"] = str(Path(args.dir_out_phewas_cache).expanduser().resolve())
    return root, env


def invoke(script, env, *args):
    path = CODE / script
    command = {".R": [env["FLT_RSCRIPT"], "--vanilla"], ".sh": ["bash"],
               ".py": [env["FLT_PYTHON"]]}[path.suffix]
    subprocess.run(command + [str(path), *args], cwd=CODE.parent, env=env, check=True)


def worker(plan_path, name):
    path = Path(plan_path)
    plan = json.loads(path.read_text())
    job = plan["jobs"][name]
    stage = job["stage"]
    env = dict(os.environ, **plan["env"])
    env.update(job["env"])
    temp = Path(env["FLT_PAPER_TMP"] if stage in PAPER else env["FLT_TMP"]) / name
    temp.mkdir(parents=True, exist_ok=True)
    env.update(WORKFLOW_CORES=str(job["cores"]), R_DATATABLE_NUM_THREADS=str(job["cores"]),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               TMPDIR=str(temp), VROOM_TEMP_PATH=str(temp), VROOM_THREADS="1")
    for key in ("FLT_OUT", "FLT_QC", "FLT_REPORTS", "FLT_PAPER", "FLT_PAPER_REPORTS",
                "FLT_PREP", "FLT_COHORT", "FLT_SMR_INPUT"):
        Path(env[key]).mkdir(parents=True, exist_ok=True)
    for dep in job["deps"]:
        if stage == "tables":
            continue
        state = json.loads((path.parent / f"{dep}.status.json").read_text())["status"]
        if state not in {"done", "skipped"}:
            raise RuntimeError(f"Incomplete prerequisite {dep}")

    report = Path(plan["qc"]) / f"{name}.json"
    result = path.parent / f"{name}.status.json"

    def status(state, detail=""):
        data = dict(stage=stage, status=state, detail=detail, plan=str(path))
        atomic_json(result, data)
        atomic_json(report, data)
        print(f"[{state.upper()}] {name}" + (f" {detail}" if detail else ""), flush=True)

    # Guard against two plans writing the same scientific stage concurrently.
    with (Path(plan["qc"]) / f"{name}.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status("running")
        try:
            os.environ.update({k: v for k, v in env.items() if k.startswith("FLT_")})
            reason = run_stage(stage, env)
        except Exception as exc:
            status("failed", str(exc))
            raise
        status("skipped" if reason else "done", reason or "")


def run_stage(stage, env):
    from external_inputs import check_inputs
    from resources import inspect, allow_skips

    issues = inspect()
    env["FLT_SMR"] = os.environ["FLT_SMR"]
    env["FLT_RUN_HDL"] = "FALSE" if "hdl" in issues else "TRUE"
    env["FLT_SMR_TISSUES"] = " ".join(x for x in ("eQTLGen", "Whole_Blood") if x not in issues)
    if stage.startswith("uvmr") or stage in {"mvmr", "smr"}:
        reason = issues.get("mr_ld")
        if stage == "smr":
            reason = reason or issues.get("smr") or (
                "No cis-eQTL BESD panel available" if not env["FLT_SMR_TISSUES"] else None)
        if reason:
            if allow_skips():
                if stage == "uvmr":
                    # Publish explicit availability independently of archived results.
                    (Path(env["FLT_REPORTS"]) / "uvmr-unavailable.txt").write_text(reason + "\n")
                return reason
            raise RuntimeError(reason)

    if stage == "phewas":
        from paper_inputs import review_fuma
        if reason := review_fuma(env):
            return reason
        errors, _ = check_inputs(["phewas"])
        if errors:
            return " | ".join(errors[:8])
    if stage in PAPER:
        from paper_inputs import missing_inputs
        if missing := missing_inputs(stage, env):
            return "Missing " + ", ".join(missing)
        if stage == "fig5":
            return figure5(env)
        if stage == "fig4":
            return figure4(env)
    if stage == "mvmr" and (Path(env["FLT_REPORTS"]) / "uvmr-unavailable.txt").exists():
        return "UVMR prerequisites unavailable"

    scripts = ANALYSIS.get(stage) or PAPER.get(stage) or INSTALL[stage]
    if stage in {"ldsc-input", "ldsc-pairs", "ldsc", "uvmr-instruments", "uvmr-forward", "uvmr-reverse"}:
        env["FLT_PHASE"] = stage
    if stage == "install":
        install_root = Path(env["FLT_R_LIBRARY"]).parent
        install_root.mkdir(parents=True, exist_ok=True)
        with (install_root / "install.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            invoke("prepare_resources.py", env)
            invoke("install.R", env)
            invoke("install.R", env, "--check")
            invoke("check_ready.py", env)
    else:
        for script in scripts:
            invoke(script, env, *(['--check'] if stage == "check" and script == "install.R" else []))
    if stage == "uvmr":
        (Path(env["FLT_REPORTS"]) / "uvmr-unavailable.txt").unlink(missing_ok=True)


def figure4(env):
    invoke("s1_5_mri_lm.plot.R", env)
    fallback = {}
    try:
        invoke("s1_5_mri_lm.plot.py", env)
    except (OSError, subprocess.CalledProcessError):
        fallback["FLT_MRI_HEATMAP_FALLBACK"] = "TRUE"
    invoke("s9_fig4.R", dict(env, **fallback))
    report = Path(env["FLT_PAPER_REPORTS"]) / "fig4-anatomy.txt"
    if fallback:
        report.write_text("Anatomical renderer incomplete. The MRI coefficient heatmaps were used. See the job log.\n")
    else:
        report.unlink(missing_ok=True)


def figure5(env):
    from external_inputs import check_inputs
    from paper_inputs import available, review_fuma

    local = Path(env["FLT_OUT"])
    out = Path(env["FLT_PAPER"])
    report = Path(env["FLT_PAPER_REPORTS"]) / "fig5-coverage.txt"
    reason = review_fuma(env)
    fuma_ind = Path(env["FLT_FUMA_DIR"]) / "IndSigSNPs.txt"
    env["FLT_FUMA_MATCHED"] = "FALSE" if reason or not fuma_ind.is_file() else "TRUE"
    (out / "Fig5.png").unlink(missing_ok=True)
    invoke("s5_1_fuma_manhattan_plot.R", env)
    if reason:
        report.write_text(f"PLACO Manhattan plot generated. Full annotation pending: {reason}\n")
        return f"PLACO Manhattan generated; annotation pending ({reason})"
    errors, _ = check_inputs(["fig5"])
    if errors:
        report.write_text("\n".join(errors) + "\n")
        return "FUMA annotation is incomplete; see fig5-coverage.txt"
    analysis_state = Path(env["FLT_QC"]) / "flt-phewas.json"
    phewas_ready = (available(analysis_state) and
                    json.loads(analysis_state.read_text())["status"] == "done" and
                    available(local / "03phewas/lead_phewas_results.tsv"))
    cloud = Path(env["FLT_PAPER_TMP"]) / "phewas/plot_audit.rds"
    cloud.unlink(missing_ok=True)
    if phewas_ready:
        invoke("s3_3_placo_phewas_plot.R", env)
    smr_state = Path(env["FLT_QC"]) / "flt-smr.json"
    smr_ready = available(smr_state) and json.loads(smr_state.read_text())["status"] == "done"
    env["FLT_SMR_CURRENT"] = "TRUE" if smr_ready else "FALSE"
    invoke("s5_2_fuma_gene_plot.R", env)
    invoke("s5_3_fuma_enrichment_plot.R", env)
    invoke("s9_fig5.R", env)
    meta = Path(env["FLT_METASCAPE_RESULTS"]).is_file()
    report.write_text(f"PLACO and FUMA lead variants match.\n"
                      f"PheWAS table available: {phewas_ready}\n"
                      f"Metascape input available: {meta}\n")


def write_plan(base, env, jobs, module, queue):
    jobs_root = base / "jobs"
    jobs_root.mkdir(parents=True, exist_ok=True)
    with (jobs_root / "dispatch.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        number = 1
        while True:
            folder = jobs_root / (module if number == 1 else f"{module}-{number}")
            try:
                folder.mkdir()
                break
            except FileExistsError:
                number += 1
    log_dir = base / "logs" / folder.name
    log_dir.mkdir(parents=True)
    qc = base / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    path = folder / "plan.json"
    atomic_json(path, {"env": env, "jobs": jobs, "qc": str(qc)})
    tasks = []
    for name, job in jobs.items():
        cmd = folder / f"{name}.cmd"
        argv = [env["FLT_PYTHON"], str(Path(__file__).resolve()), "worker", str(path), name]
        cmd.write_text("#!/usr/bin/env bash\nset -euo pipefail\nexec " + shlex.join(argv) + "\n")
        lsf = ["bsub", "-J", name, "-q", queue, "-n", str(job["slots"]), "-R", "span[hosts=1]",
               "-oo", str(log_dir / f"{name}.out"), "-eo", str(log_dir / f"{name}.err")]
        tasks.append(dict(key=name, path=cmd, argv=lsf, depends=[(dep, "ended" if job["stage"] == "tables" else "done")
                                 for dep in job["deps"]]))
    return path, tasks, log_dir


def reject_active_jobs(path, names):
    previously_submitted = {}
    for receipt in path.parent.parent.glob("*/submission.json"):
        for name, row in json.loads(receipt.read_text()).items():
            if name not in names:
                continue
            local_status = receipt.parent / f"{name}.status.json"
            if local_status.is_file() and json.loads(local_status.read_text())["status"] in {"done", "skipped", "failed"}:
                continue
            if row["status"] in {"unknown", "submitting"}:
                raise RuntimeError(f"Earlier submission of {name} is uncertain; inspect {receipt}")
            if row["status"] == "submitted":
                previously_submitted[row["job_id"]] = name
    if not previously_submitted:
        return
    result = subprocess.run(["bjobs", "-a", "-noheader", "-o", "jobid stat", *previously_submitted],
                            capture_output=True, text=True)
    rows = [line.split() for line in result.stdout.splitlines() if len(line.split()) == 2]
    active = [previously_submitted[job_id] for job_id, state in rows
              if job_id in previously_submitted and state not in {"DONE", "EXIT", "ZOMBI"}]
    if active:
        raise RuntimeError("Earlier LSF jobs still active: " + ", ".join(active))


def dispatch(path, tasks, logs, wait):
    if os.environ["submit"] == "local":
        failed = set()
        for task in tasks:
            name = task["key"]
            if any(dep in failed and condition == "done" for dep, condition in task["depends"]):
                failed.add(name)
                print(f"[BLOCKED] {name}")
                continue
            with (logs / f"{name}.out").open("w") as out, (logs / f"{name}.err").open("w") as err:
                rc = subprocess.run(["bash", str(task["path"])], stdout=out, stderr=err).returncode
            if rc:
                failed.add(name)
                print(f"[FAILED] {name} {logs / (name + '.err')}", flush=True)
            else:
                state = json.loads((path.parent / f"{name}.status.json").read_text())
                print(f"[{state['status'].upper()}] {name} {state['detail']}".rstrip(), flush=True)
        return int(bool(failed))

    receipt = path.parent / "submission.json"
    receipts = {}

    def persist(name, row):
        receipts[name] = row
        atomic_json(receipt, receipts)

    # A nested LSF environment must not replace every distinct stage name.
    os.environ.pop("WORKFLOW_JOB_NAME", None)
    with (path.parent.parent / "submit.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        reject_active_jobs(path, {task["key"] for task in tasks})
        submit_graph(tasks, {}, persist=persist, receipt=receipt)
    if not wait:
        return 0
    remaining = set(receipts)
    while remaining:
        for name in list(remaining):
            state_path = path.parent / f"{name}.status.json"
            if state_path.is_file():
                status = json.loads(state_path.read_text())["status"]
                if status == "failed":
                    raise RuntimeError(f"{name} failed. See {logs}.")
                if status in {"done", "skipped"}:
                    remaining.remove(name)
        if remaining:
            ids = [receipts[name]["job_id"] for name in remaining]
            states = subprocess.run(["bjobs", "-a", "-noheader", "-o", "jobid stat", *ids],
                                    text=True, capture_output=True)
            found = dict(line.split() for line in states.stdout.splitlines() if len(line.split()) == 2)
            for name in remaining:
                stat = found.get(receipts[name]["job_id"])
                if stat in {"EXIT", "ZOMBI", "UNKW"} or stat is None:
                    raise RuntimeError(f"LSF {name} status {stat}; inspect {receipt} and {logs}")
            time.sleep(8)
    return 0


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(sys.argv[2], sys.argv[3])
        return 0
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("entry")
    parser.add_argument("action")
    parser.add_argument("module", nargs="?")
    parser.add_argument("--dir-output")
    parser.add_argument("--dir-out-phewas-cache")
    parser.add_argument("--cores", type=int, default=int(os.environ["WORKFLOW_CORES"]))
    parser.add_argument("--whole-node", action="store_true")
    parser.add_argument("--panel", choices=["primary", "extended"], default="primary")
    parser.add_argument("--queue", default=os.getenv("HPC_CPU_QUEUE", "33"))
    parser.add_argument("--wait", action="store_true")
    mode = parser.add_mutually_exclusive_group()
    for flag in ("--plan-only", "--dry-run", "--check-inputs"):
        mode.add_argument(flag, action="store_true")
    args = parser.parse_args()
    if args.cores < 1:
        parser.error("--cores must be positive")
    if args.entry == "install.sh" and args.action in {"run", "check"}:
        module = args.action
    elif args.action == "run":
        module = args.module or "run"
    else:
        parser.error("Use run [MODULE]; see --help")
    allowed = (set(PAPER) if args.entry == "flt.paper.sh" else
               set(ANALYSIS) - {"uvmr-instruments", "uvmr-forward", "uvmr-reverse", "ldsc-input", "ldsc-pairs"})
    if module not in allowed | {"run", "check"} or module == "check" and args.entry != "install.sh":
        parser.error(f"Unknown module {module}")
    root, env = software_environment(args)
    slots = args.cores if args.whole_node else None
    if args.whole_node:
        if os.environ["submit"] != "hpc" or args.queue != "33":
            parser.error("--whole-node applies to HPC queue 33")
        slots = max(args.cores, int(os.getenv("QIMING_33_CORES", "40")))
    jobs = graph(args.entry, module, args.panel, args.cores, slots)
    if slots is None:
        for job in jobs.values():
            job["slots"] = job["cores"]
    base = root / "paper" if args.entry == "flt.paper.sh" else root
    if args.check_inputs:
        from paper_inputs import check_graph
        return check_graph(jobs, env)
    if args.dry_run:
        for name, job in jobs.items():
            print(name, job["slots"], ",".join(job["deps"]) or "-", sep="\t")
        return 0
    if args.entry == "install.sh" and module == "run" and os.environ["submit"] == "hpc" and not args.plan_only:
        invoke("prefetch.R", dict(os.environ, **env))
    path, tasks, logs = write_plan(base, env, jobs, module, args.queue)
    if args.plan_only:
        print(path)
        return 0
    return dispatch(path, tasks, logs, args.wait)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, ValueError, OSError, subprocess.SubprocessError) as error:
        sys.exit(f"FLT: {error}")
