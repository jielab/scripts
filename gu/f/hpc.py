#!/usr/bin/env python3
"""Portable, receipt-checked LSF orchestration for GU; no third-party imports.

Scientific code remains in gu.sh/arg.sh. Plans are immutable, each analysis
task owns its entire output root, and retries re-enter the scientific cache
validator. Sample batches are separate frequency-estimation cohorts.
"""
from __future__ import annotations

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path

import argparse
import contextlib
import csv
import datetime as dt
import fcntl
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
from typing import Iterator

VERSION = "gu-2026-lsf-v1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "0f"))
from qiming_resources import cpu_queue, cpu_slots
SOURCE_SUFFIXES = {".sh", ".py", ".R", ".cc", ".cu", ".h", ".yml", ".env"}
LIVE_STATES = {"PEND", "RUN", "WAIT", "PSUSP", "USUSP", "SSUSP", "UNKWN"}


def die(message: str) -> None:
    raise ValueError(message)


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def digest(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def atomic_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, path)


@contextlib.contextmanager
def locked(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            die(f"Another process owns {path}; wait for it to finish.")
        yield


def read_json(path: Path):
    return json.loads(path.read_text())


def chroms(value: str) -> list[str]:
    result = []
    for token in re.split(r"[,\s]+", value.strip()):
        token = re.sub(r"^chr", "", token, flags=re.I).upper()
        if re.fullmatch(r"\d+-\d+", token):
            lo, hi = map(int, token.split("-"))
            if not 1 <= lo <= hi <= 22:
                die("Chromosome ranges must be between 1 and 22.")
            parts = map(str, range(lo, hi + 1))
        else:
            parts = ["X" if token == "23" else token]
        for c in parts:
            if c not in {str(n) for n in range(1, 23)} | {"X"}:
                die(f"Unsupported chromosome {c!r}.")
            if c not in result:
                result.append(c)
    if not result:
        die("No chromosomes selected.")
    return result


def positive(value: str) -> int:
    n = int(value)
    if n < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return n


def persistent_path(value: str) -> Path:
    path = Path(value).expanduser().resolve()
    if str(path) == "/" or any(path == p or p in path.parents for p in map(Path, ["/tmp", "/var/tmp", "/dev/shm", "/run"])):
        die(f"Use a persistent result directory under /data or /work, not {path}.")
    return path


def code_manifest() -> dict[str, str]:
    result = {}
    # Parent helpers are part of the downloadable package and execution contract.
    for base in [ROOT / "f", ROOT.parent / "0f"]:
        if base.is_dir():
            for path in sorted(base.rglob("*")):
                if path.is_file() and path.suffix in SOURCE_SUFFIXES and not any(x in path.parts for x in ("__pycache__", ".git")):
                    result[str(path)] = sha(path)
    for path in ROOT.iterdir():
        if path.is_file() and path.suffix in SOURCE_SUFFIXES:
            result[str(path)] = sha(path)
    return result


def input_stat(path: Path) -> dict:
    path = path.resolve()
    st = path.stat()
    if st.st_size == 0:
        die(f"Input file is empty: {path}")
    rec = {"path": str(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    # Hash small metadata; multi-terabyte genotype files use size/mtime plus
    # the scientific program's own receipts, never a claim of full-file hashing.
    if path.suffix in {".psam", ".txt", ".tsv", ".bed", ".cojo", ".env"} and st.st_size < 100_000_000:
        rec["sha256"] = sha(path)
    return rec


def verify_stat(record: dict) -> None:
    path = Path(record["path"])
    if not path.is_file():
        die(f"Input disappeared: {path}")
    current = input_stat(path)
    if current != record:
        die(f"Input changed since planning: {path}. Create a new plan.")


def read_psam(path: Path) -> list[tuple[str, str, str]]:
    result = []
    seen = set()
    with path.open() as f:
        line = next((s for s in f if s.strip()), "")
        header = line.strip().lstrip("#").split()
        if "IID" not in header:
            die(f"PSAM requires an IID header: {path}")
        iid = header.index("IID")
        fid = header.index("FID") if "FID" in header else None
        sex = header.index("SEX") if "SEX" in header else None
        for number, line in enumerate(f, 2):
            if not line.strip():
                continue
            row = line.split()
            if len(row) < len(header):
                die(f"Malformed PSAM row {number}: {path}")
            identifier = row[iid]
            if identifier in seen:
                die(f"Duplicate IID {identifier!r} in {path}; GU exports by IID.")
            seen.add(identifier)
            result.append((row[fid] if fid is not None else "0", identifier, row[sex] if sex is not None else "NA"))
    if not result:
        die(f"PSAM has no samples: {path}")
    return result


def freeze_batches(rows: list[tuple[str, str, str]], size: int, males: bool) -> list[dict]:
    batches = []
    for start in range(0, len(rows), size):
        selected = rows[start:start + size]
        kept = [r for r in selected if r[2] == "1"] if males else selected
        batches.append({"batch": start // size + 1, "first_row": start + 1,
                        "last_row": start + len(selected), "selected_n": len(selected),
                        "scored_n": len(kept), "rows": kept})
    return batches


def discover_genotypes(prefix: str, chromosomes: list[str]) -> list[dict]:
    records = []
    for c in chromosomes:
        candidates = [Path(prefix + c)]
        if c == "X":
            candidates += [Path(prefix + "X.male"), Path(prefix + "XY")]
        found = False
        for base in candidates:
            if Path(str(base) + ".pgen").is_file():
                pvar = next((Path(str(base) + s) for s in (".pvar", ".pvar.zst") if Path(str(base) + s).is_file()), None)
                if pvar is None or not Path(str(base) + ".psam").is_file():
                    die(f"missing {base}.pvar[.zst] or {base}.psam (incomplete PLINK2 dataset)")
                for p in [Path(str(base) + ".pgen"), pvar, Path(str(base) + ".psam")]:
                    records.append({"chrom": c, **input_stat(p)})
                found = True
            for suffix in (".vcf.gz", ".vcf.bgz", ".bcf"):
                p = Path(str(base) + suffix)
                if p.is_file():
                    idx = next((Path(str(p) + s) for s in (".csi", ".tbi") if Path(str(p) + s).is_file()), None)
                    if idx is None:
                        die(f"missing {p}.csi or {p}.tbi (input must be indexed before submission)")
                    records += [{"chrom": c, **input_stat(v)} for v in (p, idx)]
                    found = True
        if not found:
            die(f"missing {prefix}{c}.pgen (+ .pvar[.zst], .psam), or {prefix}{c}.vcf.gz (+ .tbi/.csi)")
    return records


def snapshot_file(source: str | None, run: Path, name: str) -> str | None:
    if source is None:
        return None
    p = Path(source).expanduser().resolve()
    if not p.is_file() or not p.stat().st_size:
        die(f"Required file is missing/empty: {p}")
    dest = run / "inputs" / name
    shutil.copyfile(p, dest)
    return str(dest)


def resolved_runtime_environment() -> dict[str, str]:
    """Keep super-roots/runtime paths when saved plans are submitted later."""
    prefixes = ("DIR_", "OS_", "_OS_DEFAULT_")
    exact = {"ANALYSIS_ROOT", "SCRATCH_ANALYSIS_ROOT", "PERSISTENT_ANALYSIS_ROOT", "SHARED_PYTHON",
             "GU_DATA_ROOT", "GU_REF_ROOT", "GU_ANNOT_ROOT", "GU_ANALYSIS_BASE", "GU_SOFT",
             "GU_CONDA_BASE", "GU_ENV_PREFIX", "GU_ENV_NAME", "GU_CONDA_EXE", "GU_FILES_ROOT",
             "GU_COJO_CHAIN", "GU_LIFTOVER", "GU_RESULTS_RSCRIPT", "CHECK_GRCH_SNP_LIST",
             "IBDMIX_GPU_BIN", "IBDMIX_MASK_ROOT", "IBDMIX_REFERENCE_FASTA", "IBDMIX_AXT_DIR",
             "IBDMIX_RUNTIME", "IBDMIX_PREPARE_BIN", "PHYML_TREE_SEED"}
    return {key:value for key,value in os.environ.items() if key in exact or key.startswith(prefixes)}


def common_gu(plan: dict, output: str, selection: str | None = None) -> list[str]:
    args = ["bash", str(ROOT / "gu.sh"), plan["method"], "--grch", plan["grch"],
            "--target", plan["target"], "--dir-target", plan["dir_target"],
            "--dir-out", output, "--foreground", "TRUE", "--run-cmd", "TRUE",
            "--jobs", "1", "--auto-final", "FALSE", "--memory-cap", f"{plan['memory_gb']}G"]
    if plan["method"] == "ibdmix":
        args += ["--GPU", "TRUE" if plan.get("gpu", False) else "FALSE", "--replace-ibdmix", "FALSE"]
    if selection:
        args += ["--keep", selection]
    if plan["keep_males"]:
        args += ["--keep-males"]
    if plan.get("sample_panel"):
        args += ["--sample-panel", plan["sample_panel"]]
    if plan["method"] == "ibdmix" and plan.get("archaic_path"):
        args += ["--archaic-path", plan["archaic_path"]]
    if plan["method"] == "ibdmix" and plan.get("archaic_gen"):
        args += ["--archaic-gen", plan["archaic_gen"]]
    return args


def plan_run(args) -> None:
    if not re.fullmatch(r"[A-Za-z0-9_-][A-Za-z0-9_.-]*", args.target):
        die("Target must be a plain dataset name.")
    analysis = Path(os.environ.get("GU_ANALYSIS_BASE", str(Path(os.environ.get("ANALYSIS_ROOT", configured_path('ANALYSIS_ROOT'))) / "gu")))
    default_run = analysis / args.target / args.grch / args.method if args.method == "ibdmix" else analysis / args.target / args.method
    run = persistent_path(args.run_dir or str(default_run))
    if run.exists() and any(run.iterdir()):
        die(f"Plan directory is not empty: {run}. Use a new run name or submit the existing plan.")
    chromosomes = chroms(args.chr)
    gpu = (args.GPU or (os.getenv("GU_IBDMIX_GPU_DEFAULT", "FALSE") if args.method == "ibdmix" else "FALSE")) == "TRUE"
    if gpu and args.method != "ibdmix":
        die("hpc plan --GPU TRUE is implemented for ibdmix only; other methods need their dedicated backend.")
    args.queue = args.queue or (os.getenv("GU_GPU_QUEUE", os.getenv("HPC_GPU_QUEUE", "4a100-40"))
                               if gpu else os.getenv("GU_QUEUE", "33"))
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.queue):
        die("Queue must be a plain LSF queue name.")
    if not re.fullmatch(r"\d+:[0-5]\d", args.walltime):
        die("Walltime must be HHH:MM.")
    if args.memory_gb < 4:
        die("Use at least 4 GiB per scientific worker.")
    if not gpu:
        hours, minutes = map(int, args.walltime.split(':'))
        args.queue = cpu_queue(args.queue, args.threads + 1, args.memory_gb, hours + minutes / 60)
    if args.method in {"trace", "arg", "arg-trace"} and (args.psam or args.keep_males):
        die("ARG/TRACE plans use a fixed complete cohort, not IBDmix sample batches; X is handled by its dedicated ARG contract.")
    if args.method == "phyml" and not args.loci:
        die("PhyML requires --loci FILE.jma.cojo; whole-genome tree scans are not enabled here.")
    if args.method in {"ibdmix", "phyml"} and not args.archaic_path:
        die("Supply --archaic-path explicitly for HPC reference portability.")
    if args.method == "trace" and not args.arg_dir:
        die("TRACE requires --dir-arg with certified trees; use arg-trace to build them first.")
    if args.method in {"arg", "arg-trace"} and not args.dir_gen:
        die("ARG construction requires --dir-gen.")
    if args.method == "phyml" and args.psam:
        die("HPC sample batching is limited to ibdmix; use a fixed PhyML cohort.")
    if args.method != "ibdmix" and args.batch_size != 10000:
        die("--batch-size is only valid for ibdmix.")
    records = discover_genotypes(args.dir_target, chromosomes)
    code = code_manifest()
    if not (ROOT / "gu.sh").is_file():
        die("Install gu-hpc.sh alongside gu.sh before planning.")
    run.mkdir(parents=True, exist_ok=True)
    for name in ("inputs", "jobs", "logs", "state", "tasks", "summary"):
        (run / name).mkdir(exist_ok=True)
    env_file = snapshot_file(args.env_file, run, "runtime.env")
    panel = snapshot_file(args.sample_panel, run, "samples.txt")
    loci = snapshot_file(args.loci, run, "loci.jma.cojo")
    plan = {"version": VERSION, "created": now(), "run_dir": str(run),
            "date_scratch": os.environ.get("DATE_SCRATCH", "2026-11-03"),
            "source_root": str(ROOT), "source_sha256": code, "environment": env_file,
            "resolved_environment": resolved_runtime_environment(),
            "environment_sha256": sha(Path(env_file)), "method": args.method,
            "grch": args.grch, "target": args.target, "dir_target": str(Path(args.dir_target).resolve()),
            "dir_gen": str(Path(args.dir_gen).resolve()) if args.dir_gen else None,
            "arg_dir": str(Path(args.arg_dir).resolve()) if args.arg_dir else None,
            "arg_method": args.arg_method, "queue": args.queue, "gpu": gpu,
            "cpu_queue": os.getenv("HPC_CPU_QUEUE", "33"),
            "phyml_timeout_seconds": args.phyml_timeout_seconds,
            "phyml_retry_timeouts": args.phyml_retry_timeouts,
            "dir_map": str(Path(args.dir_map or Path(configured_path('DIR_GEN')) / 'hm3' / args.grch / 'gmap').resolve()),
            "threads": args.threads, "slots": args.threads + 1,
            "memory_gb": args.memory_gb, "walltime": args.walltime,
            "max_running": args.max_running, "array_size": args.array_size,
            "rusage_memory_mode": args.rusage_memory_mode,
            "chromosomes": chromosomes, "keep_males": args.keep_males,
            "sample_panel": panel, "loci": loci,
            "snapshot_sha256": {str(p): sha(Path(p)) for p in (panel, loci) if p},
            "archaic_path": str(Path(args.archaic_path).resolve()) if args.archaic_path else None,
            "archaic_gen": args.archaic_gen, "input_files": records,
            "batch_size": args.batch_size, "batch_frequency_policy": "independent-cohorts",
            "tasks": [], "batches": []}
    if args.psam:
        psam = Path(args.psam).resolve()
        plan["psam"] = input_stat(psam)
        batches = freeze_batches(read_psam(psam), args.batch_size, args.keep_males)
        for batch in batches:
            rows = batch.pop("rows")
            keep = run / "inputs" / f"batch{batch['batch']:05d}.keep"
            if rows:
                keep.write_text("#FID\tIID\n" + "".join(f"{fid}\t{iid}\n" for fid, iid, sex in rows))
                batch["keep"] = str(keep)
                batch["keep_sha256"] = sha(keep)
            else:
                batch["keep"] = None
                batch["excluded_reason"] = "zero_selected_males"
            plan["batches"].append(batch)
        plan["batches"] = batches
    else:
        plan["batches"] = [{"batch": 1, "keep": None, "scored_n": None}]
    task_id = 0
    if args.method == "ibdmix":
        for batch in plan["batches"]:
            if batch.get("excluded_reason"):
                continue
            for c in chromosomes:
                task_id += 1
                output = run / "tasks" / f"batch{batch['batch']:05d}" / f"chr{c}"
                command = common_gu(plan, str(output), batch["keep"]) + ["--chr", c]
                plan["tasks"].append({"id": task_id, "kind": "ibdmix", "batch": batch["batch"],
                                      "chromosomes": [c], "output": str(output), "command": command})
    elif args.method in {"arg", "arg-trace"}:
        for c in chromosomes:
            task_id += 1
            output = run / "tasks" / f"chr{c}"
            command = ["bash", str(ROOT / "arg.sh"), "build", "--method", args.arg_method,
                       "--format", "trace", "--dir-gen", plan["dir_gen"],
                       "--dir-pfile", str(Path(plan["dir_target"]).parent), "--dir-map", plan["dir_map"],
                       "--dir-arg", str(output / "arg"), "--dir-gen4arg", str(output / "prepared"),
                       "--chr", c, "--grch", args.grch, "--threads", str(args.threads),
                       "--jobs", "1", "--resource-guard", "scheduler",
                       "--memory-limit-gb", str(args.memory_gb)]
            if panel:
                command += ["--sample-panel", panel]
                if args.arg_method in {"threads", "needle"}:
                    # These engines select the cohort from the ancestry panel.
                    # Its identity must match the frozen TRACE sample metadata.
                    command += ["--ancestry-file", panel]
            plan["tasks"].append({"id": task_id, "kind": "arg", "batch": 1,
                                  "chromosomes": [c], "output": str(output), "command": command})
        if args.method == "arg-trace":
            trace_plan = {**plan, "method": "trace", "keep_males": False, "archaic_path": None}
            command = common_gu(trace_plan, str(run / "trace-result")) + ["--chr", ",".join(chromosomes), "--dir-arg", str(run / "arg-assembled")]
            plan["trace_task"] = {"id": "trace", "kind": "trace", "batch": 1, "chromosomes": chromosomes,
                                  "output": str(run / "trace-result"), "command": command}
    else:
        output = run / "tasks" / args.method
        command = common_gu(plan, str(output))
        if args.method == "trace":
            command += ["--chr", ",".join(chromosomes), "--dir-arg", plan["arg_dir"]]
        else:
            command += ["--loci", loci, "--loci-format", "cojo", "--phyml-jobs", "1"]
        plan["tasks"].append({"id": 1, "kind": args.method, "batch": 1, "chromosomes": chromosomes,
                              "output": str(output), "command": command})
    if not plan["tasks"]:
        die("No nonempty sample batches remain.")
    nblocks = math.ceil(len(plan["tasks"]) / plan["array_size"])
    if nblocks > plan["max_running"]:
        die("Array block count exceeds --max-running; increase --array-size or --max-running.")
    plan["plan_id"] = digest(plan)[:20]
    atomic_json(run / "plan.json", plan)
    render_jobs(plan)
    print(json.dumps({"plan": str(run / "plan.json"), "tasks": len(plan["tasks"]),
                      "batch_count": len(plan["batches"]), "queue": plan["queue"],
                      "submission": "NOT SUBMITTED", "review": str(run / "jobs" / "commands.tsv")}, indent=2))


def load_plan(run: Path) -> dict:
    plan = read_json(run / "plan.json")
    identity = plan.pop("plan_id")
    if digest(plan)[:20] != identity:
        die("Plan was edited after creation; create a new plan instead.")
    plan["plan_id"] = identity
    if plan["version"] != VERSION or Path(plan["run_dir"]) != run.resolve():
        die("Plan version or location changed; create a new plan at its final HPC location.")
    return plan


def render_jobs(plan: dict) -> None:
    run = Path(plan["run_dir"])
    for stage in ("prepare", "worker", "trace", "merge"):
        text = "#!/bin/bash -l\nset -euo pipefail\n"
        text += '[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit through gu-hpc.sh submit; do not execute on the login node" >&2; exit 2; }\n'
        if plan.get("date_scratch"):
            text += "export DATE_SCRATCH=" + shlex.quote(plan["date_scratch"]) + "\n"
        text += "export submit=hpc SUBMIT=hpc WORKFLOW_WORKER=1\n"
        for key, value in sorted(plan.get("resolved_environment", {}).items()):
            text += "export " + key + "=" + shlex.quote(value) + "\n"
        text += "export GU_WORKER_CODE_ROOT=" + shlex.quote(str(ROOT)) + "\n"
        text += "source " + shlex.quote(plan["environment"]) + "\n"
        text += "export GU_RESOURCE_GUARD=scheduler ARG_RESOURCE_GUARD=scheduler\n"
        text += "cd " + shlex.quote(plan["source_root"]) + "\n"
        text += "exec python3 " + shlex.quote(str(ROOT / "f" / "hpc.py")) + " _execute " + shlex.quote(str(run)) + " " + stage + ' "$@"\n'
        p = run / "jobs" / (stage + ".cmd")
        p.write_text(text)
        p.chmod(0o755)
    with (run / "jobs" / "commands.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["task", "kind", "batch", "chromosomes", "output", "command"])
        for task in plan["tasks"] + ([plan["trace_task"]] if plan.get("trace_task") else []):
            w.writerow([task["id"], task["kind"], task["batch"], ",".join(task["chromosomes"]), task["output"], shlex.join(task["command"])])


def verify_plan_files(plan: dict, task: dict | None = None) -> None:
    if sha(Path(plan["environment"])) != plan["environment_sha256"]:
        die("Frozen runtime environment changed.")
    for p, expected in plan["source_sha256"].items():
        if not Path(p).is_file() or sha(Path(p)) != expected:
            die(f"Code changed after planning: {p}; use a new plan.")
    if plan.get("psam"):
        verify_stat(plan["psam"])
    for p, expected in plan.get("snapshot_sha256", {}).items():
        if sha(Path(p)) != expected:
            die(f"Frozen input metadata changed: {p}")
    for batch in plan["batches"]:
        if batch.get("keep") and sha(Path(batch["keep"])) != batch["keep_sha256"]:
            die(f"Frozen sample selection changed: {batch['keep']}")
    for record in plan["input_files"]:
        if task is None or record["chrom"] in task["chromosomes"]:
            verify_stat({k: v for k, v in record.items() if k != "chrom"})


def require_scheduler(plan: dict, stage: str) -> None:
    if not os.environ.get("LSB_JOBID"):
        die("Scientific execution requires an allocated LSF job; use submit.")
    required = 1 if stage in {"prepare", "merge"} else plan["slots"]
    slots = int(os.environ.get("LSB_DJOB_NUMPROC", "0"))
    if slots < required:
        die(f"Allocated slots {slots} are below required {required}.")


def runtime_environment(plan: dict) -> dict[str, str]:
    env = dict(os.environ)
    if plan.get("date_scratch"):
        env["DATE_SCRATCH"] = plan["date_scratch"]
        env.setdefault("SCRATCH_ANALYSIS_ROOT", env.get("ANALYSIS_ROOT", configured_path('ANALYSIS_ROOT')))
    # A request must not inherit a previous gu.sh invocation's resolved roots,
    # sample selectors, output overrides, or GPU routing.
    for key in ("GU_CMD_WORKER", "GU_UKB_DIRECT", "GU_UKB_PGEN_PREPARED", "GU_FILTERED_TARGET",
                "GU_KEEP", "GU_KEEP_PSAM", "GU_MALE_ONLY", "GU_TARGET_ROOT", "GU_TARGET_NAMESPACE",
                "GU_BATCH_ID",
                "GU_PUBLISHED_ROOT", "GU_ANALYSIS_ROOT", "GU_UKB_RESULTS_ROOT", "GU_RESULT_MODULE",
                "GU_TARGET_NATIVE_VCF_PREFIX", "IBDMIX_OUT", "TRACE_OUT", "PHYML_OUT", "GU_ENV_REEXEC", "GU_ARG_CHECKED"):
        env.pop(key, None)
    env.update({"GU_RESOURCE_GUARD": "scheduler", "ARG_RESOURCE_GUARD": "scheduler",
                "GU_PREP_THREADS": str(plan["threads"]), "IBDMIX_PREPARE_THREADS": str(plan["threads"]),
                "IBDMIX_JOB_OF_CHR": "1", "IBDMIX_JOB_IN_CHR": "1",
                "PHYML_TREE_CPUS": str(plan["threads"]), "TRACE_JOB_EXTRACT": str(plan["threads"]),
                "TRACE_JOB_INFER": str(plan["threads"]), "TRACE_JOB_SUMMARIZE": str(plan["threads"]),
                "IBDMIX_BACKEND": "cuda" if plan.get("gpu", False) else "cpu", "TRACE_BACKEND": "cpu", "PHYML_BACKEND": "cpu",
                "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1", "GU_AUTO_FINAL": "0", "GU_KEEP_SCRATCH": "1"})
    if plan.get("archaic_path"):
        env["GU_ARCHAIC_ROOT"] = plan["archaic_path"]
    if plan["method"] == "phyml":
        env["PHYML_TREE_TIMEOUT"] = str(plan.get("phyml_timeout_seconds", 7200))
        env["PHYML_RETRY_TIMEOUTS"] = "TRUE" if plan.get("phyml_retry_timeouts", False) else "FALSE"
    return env


def prepare(plan: dict) -> None:
    verify_plan_files(plan)
    env = runtime_environment(plan)
    for program in ("python3", "bash", "bcftools", "plink2", "samtools", "flock"):
        if not shutil.which(program, path=env.get("PATH")):
            die(f"Required program is not on PATH: {program}; fix the runtime environment then make a new plan.")
    if any(r["path"].endswith(".pgen") for r in plan["input_files"]):
        marker = env.get("CHECK_GRCH_SNP_LIST")
        if not marker or not Path(marker).is_file() or not Path(marker).stat().st_size:
            die("Set CHECK_GRCH_SNP_LIST to the real build-validation SNP coordinate panel in runtime.env; it is input data, not bundled code.")
    for key in ("IBDMIX_REFERENCE_FASTA", "UKB_REF_FASTA", "GU_REF_FASTA"):
        if env.get(key) and (not Path(env[key]).is_file() or not Path(env[key]).stat().st_size):
            die(f"Configured reference input is missing/empty: {key}={env[key]}")
    if plan["method"] == "ibdmix":
        subprocess.run(["bash", str(ROOT / "f" / "ibdmix.sh"), "build-prepare", "--if-needed"], env=env, check=True)
        native = Path(env.get("IBDMIX_RUNTIME", str(Path(env.get("GU_SOFT", "")) / "ibdmix"))) / "build" / "src" / "ibdmix"
        if not os.access(native, os.X_OK):
            die(f"Native IBDmix executable is missing: {native}; install before submission.")
        if plan.get("gpu", False):
            cuda = Path(env.get("IBDMIX_GPU_BIN", str(Path(env.get("GU_SOFT", "")) / "ibdmix-gpu" / "ibdmix-gpu")))
            if not os.access(cuda, os.X_OK):
                die(f"CUDA IBDmix executable is missing: {cuda}; build on an allocated GPU node before submission.")
        refs = re.split(r"[,\s]+", plan["archaic_gen"].strip())
        mask = env.get("IBDMIX_MASK_ROOT", str(Path(plan["archaic_path"]).parent / "mask"))
        cmd = ["python3", str(ROOT / "f" / "ibdmix.py"), "check-masks", "--root", mask,
               "--archaic-root", plan["archaic_path"], "--chroms", *plan["chromosomes"], "--refs", *refs]
        if env.get("IBDMIX_MASK_DIR"):
            cmd += ["--custom-masks", env["IBDMIX_MASK_DIR"]]
        subprocess.run(cmd, check=True, env=env)
    if plan["method"] in {"trace", "arg-trace"}:
        for program in ("trace-extract", "trace-infer", "trace-summarize"):
            if not shutil.which(program, path=env.get("PATH")):
                die(f"TRACE entry point is missing: {program}; install TRACE before submission.")
    atomic_json(Path(plan["run_dir"]) / "state" / "prepare.json",
                {"plan_id": plan["plan_id"], "state": "complete", "finished": now(),
                 "job_id": os.environ["LSB_JOBID"]})


def artifact_candidates(task: dict) -> list[dict]:
    root = Path(task["output"])
    kind = task["kind"]
    if kind == "arg":
        c = task["chromosomes"][0]
        candidates = list((root / "arg" / "trace").glob(f"*/chr{c}.trees"))
        if len(candidates) != 1:
            die(f"Expected one TRACE-formatted ARG tree for chr{c}; found {len(candidates)} at {root / 'arg'}.")
        tree = candidates[0]
        required = [tree, tree.with_name(f"chr{c}.sample_map.tsv"), tree.parent / "ARG_TRACE_BUILD.tsv"]
        if any(not p.is_file() or not p.stat().st_size for p in required):
            die(f"ARG tree/sample map/build receipt is incomplete for chr{c}.")
        return [{"path": str(p), "sha256": sha(p), "size": p.stat().st_size} for p in required]
    basename = {"ibdmix": re.compile(r"all_archaic_refs\.lod[^/]+\.len[^/]+\.segments\.tsv\.gz$"),
                "trace": re.compile(r"trace_haplotype_segments\.tsv\.gz$"),
                "phyml": re.compile(r"gwas_loci\.tsv$")}[kind]
    results = []
    method_root = root / kind
    # Each parent directory identifies one native analysis unit. Prefer durable
    # archives where present; never count both an archive and its read copy.
    archived_dirs = set()
    for archive in sorted(method_root.rglob(kind + ".raw.tar.gz")):
        archived_dirs.add(archive.parent)
        with tarfile.open(archive, "r:gz") as bundle:
            for member in bundle:
                if member.isfile() and basename.fullmatch(Path(member.name).name) and "final" in Path(member.name).parts:
                    payload = bundle.extractfile(member)
                    h = hashlib.sha256()
                    for block in iter(lambda: payload.read(1024 * 1024), b""):
                        h.update(block)
                    results.append({"path": str(archive), "member": member.name, "sha256": h.hexdigest(), "size": member.size})
    for file in sorted(method_root.rglob("*")):
        if not file.is_file() or not basename.fullmatch(file.name) or "final" not in file.parts:
            continue
        if any(p in archived_dirs for p in file.parents):
            continue
        results.append({"path": str(file), "sha256": sha(file), "size": file.stat().st_size})
    if not results:
        die(f"Program returned success but no canonical {kind} output exists at {root}.")
    return results


@contextlib.contextmanager
def open_artifact(artifact: dict):
    with contextlib.ExitStack() as stack:
        if artifact.get("member"):
            archive = stack.enter_context(tarfile.open(artifact["path"], "r:gz"))
            stream = archive.extractfile(artifact["member"])
            if stream is None:
                die(f"Missing archived result {artifact['member']}")
            name = artifact["member"]
        else:
            stream = stack.enter_context(open(artifact["path"], "rb"))
            name = artifact["path"]
        if name.endswith(".gz"):
            stream = stack.enter_context(gzip.GzipFile(fileobj=stream))
        text = io.TextIOWrapper(stream, encoding="utf-8")
        yield text


def validate_artifacts(task: dict, artifacts: list[dict]) -> None:
    if task["kind"] == "arg":
        return
    required = {"ibdmix": {"ID", "chrom", "start", "end", "anc"},
                "trace": {"start", "end"}, "phyml": set()}[task["kind"]]
    for a in artifacts:
        with open_artifact(a) as f:
            reader = csv.DictReader(f, delimiter="\t")
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                die(f"Malformed final result header: {a['path']} ({a.get('member', '')})")
            count = 0
            for row in reader:
                if None in row:
                    die("Result row has too many columns.")
                if task["kind"] in {"ibdmix", "trace"} and int(row["end"]) <= int(row["start"]):
                    die("Result has a nonpositive genomic interval.")
                count += 1
            a["rows"] = count
    # A valid header and zero rows is a successful negative call, not failure.


def task_receipt(plan: dict, task: dict) -> Path:
    return Path(plan["run_dir"]) / "state" / f"task.{task['id']}.json"


def prerequisite(plan: dict) -> None:
    p = Path(plan["run_dir"]) / "state" / "prepare.json"
    if not p.is_file() or read_json(p).get("plan_id") != plan["plan_id"] or read_json(p).get("state") != "complete":
        die("Preparation failed or has no matching completion receipt; scientific work was not started.")


def execute_task(plan: dict, task: dict, assembled_checked: bool = False) -> None:
    prerequisite(plan)
    verify_plan_files(plan, task)
    receipt = task_receipt(plan, task)
    with locked(receipt.with_suffix(".lock")):
        state = {"plan_id": plan["plan_id"], "task": task["id"], "state": "running", "started": now(),
                 "job_id": os.environ["LSB_JOBID"], "array_index": os.environ.get("LSB_JOBINDEX"),
                 "command": task["command"]}
        atomic_json(receipt, state)
        try:
            env = runtime_environment(plan)
            env["GU_BATCH_ID"] = f"batch{task['batch']:05d}"
            if task["kind"] == "trace":
                arg_dir = Path(task["command"][task["command"].index("--dir-arg") + 1])
                if assembled_checked or validate_arg_assembly(plan, arg_dir):
                    # The public driver otherwise assumes ARG_DIR/../.. is one
                    # native build root, which is false for an assembled array.
                    # Every original build and final artifact was verified above;
                    # TRACE's independent tree/map/span checks remain enabled.
                    env["GU_ARG_CHECKED"] = "1"
            subprocess.run(task["command"], env=env, check=True)
            artifacts = artifact_candidates(task)
            validate_artifacts(task, artifacts)
            state.update(state="complete", finished=now(), artifacts=artifacts)
            atomic_json(receipt, state)
        except BaseException as error:
            state.update(state="failed", finished=now(), error=str(error))
            atomic_json(receipt, state)
            raise


def complete_receipts(plan: dict, tasks: list[dict]) -> list[dict]:
    result, bad = [], []
    for task in tasks:
        p = task_receipt(plan, task)
        rec = read_json(p) if p.is_file() else {}
        if rec.get("state") != "complete" or rec.get("plan_id") != plan["plan_id"]:
            bad.append(str(task["id"]))
            continue
        # Revalidate the final content against the completion receipt without
        # expanding large intermediate archives onto the filesystem.
        current = artifact_candidates(task)
        old = {(a["path"], a.get("member")): (a["sha256"], a["size"]) for a in rec["artifacts"]}
        new = {(a["path"], a.get("member")): (a["sha256"], a["size"]) for a in current}
        if old != new:
            bad.append(str(task["id"]) + ":outputs_changed")
        else:
            result.append(rec)
    if bad:
        die("Incomplete or invalid tasks: " + ",".join(bad) + ". No aggregate is marked complete.")
    return result


def assemble_arg(plan: dict, receipts: list[dict] | None = None) -> None:
    if receipts is None:
        receipts = complete_receipts(plan, plan["tasks"])
    dest = Path(plan["run_dir"]) / "arg-assembled"
    dest.mkdir(exist_ok=True)
    for task, receipt in zip(plan["tasks"], receipts):
        c = task["chromosomes"][0]
        source = Path(next(a["path"] for a in receipt["artifacts"] if a["path"].endswith(".trees"))).parent
        fields = dict(line.rstrip("\n").split("\t", 1) for line in (source / "ARG_TRACE_BUILD.tsv").read_text().splitlines() if "\t" in line)
        if fields.get("method") != plan["arg_method"] or fields.get("build") != "GRCh" + plan["grch"] or fields.get("format") != "trace":
            die(f"ARG build receipt does not match plan for chr{c}.")
        # Keep method-specific provenance and QA files alongside the tree/map;
        # names beginning chrN cannot collide with another chromosome.
        for file in source.iterdir():
            if re.match(r"^chr" + re.escape(c) + r"[._]", file.name) and file.is_file():
                link = dest / file.name
                if link.is_symlink() and link.resolve() == file.resolve():
                    continue
                if link.exists() or link.is_symlink():
                    die(f"Conflicting assembled ARG artifact: {link}")
                link.symlink_to(file.resolve())
    meta = "method\t" + plan["arg_method"] + "\nformat\ttrace\nbuild\tGRCh" + plan["grch"] + "\nchromosomes\t" + " ".join(plan["chromosomes"]) + "\n"
    (dest / "ARG_TRACE_BUILD.tsv").write_text(meta)
    atomic_json(dest / "HPC_ARG_ASSEMBLY.json", {"format": "gu-hpc-arg-v1", "source_run": plan["run_dir"],
                "source_plan_id": plan["plan_id"], "grch": plan["grch"], "method": plan["arg_method"],
                "chromosomes": plan["chromosomes"], "metadata_sha256": sha(dest / "ARG_TRACE_BUILD.tsv"),
                "artifacts": {p.name: {"source": str(p.resolve()), "sha256": sha(p)}
                              for p in dest.iterdir() if p.is_symlink() and p.is_file()}})


def validate_arg_assembly(plan: dict, directory: Path) -> bool:
    marker = directory / "HPC_ARG_ASSEMBLY.json"
    if not marker.is_file():
        return False
    info = read_json(marker)
    source = load_plan(Path(info["source_run"]))
    if (info.get("format") != "gu-hpc-arg-v1" or source["method"] not in {"arg", "arg-trace"}
            or source["plan_id"] != info["source_plan_id"] or info["grch"] != plan["grch"]
            or not set(plan["chromosomes"]).issubset(info["chromosomes"])):
        die("HPC ARG assembly identity/build/chromosome contract does not match this TRACE plan.")
    verify_plan_files(source)
    complete_receipts(source, source["tasks"])
    if sha(directory / "ARG_TRACE_BUILD.tsv") != info["metadata_sha256"]:
        die("HPC ARG assembly build metadata changed.")
    for name, artifact in info["artifacts"].items():
        file = directory / name
        if file.resolve() != Path(artifact["source"]) or sha(file) != artifact["sha256"]:
            die("HPC ARG assembly artifact changed: " + name)
    return True


def merge(plan: dict) -> None:
    run = Path(plan["run_dir"])
    prerequisites = complete_receipts(plan, plan["tasks"])
    if plan["method"] == "arg":
        assemble_arg(plan, prerequisites)
        atomic_json(run / "summary" / "COMPLETE.json", {"plan_id": plan["plan_id"], "finished": now(), "arg_dir": str(run / "arg-assembled")})
        return
    tasks = [plan["trace_task"]] if plan.get("trace_task") else plan["tasks"]
    receipts = complete_receipts(plan, tasks) if plan.get("trace_task") else prerequisites
    method = tasks[0]["kind"]
    stage = run / "summary" / f".merge.{os.getpid()}"
    stage.mkdir()
    try:
        out = stage / ("segments.tsv.gz" if method in {"ibdmix", "trace"} else "phyml_loci.tsv.gz")
        fields = None
        counts = []
        with gzip.open(out, "wt", encoding="utf-8") as handle:
            writer = None
            for task, receipt in zip(tasks, receipts):
                count, length_sum = 0, 0
                for artifact in receipt["artifacts"]:
                    with open_artifact(artifact) as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        if fields is None:
                            fields = reader.fieldnames
                            if any(x in fields for x in ("hpc_task", "hpc_batch")):
                                die("Input already contains HPC aggregate columns.")
                            writer = csv.DictWriter(handle, fieldnames=["hpc_task", "hpc_batch", *fields], delimiter="\t", lineterminator="\n")
                            writer.writeheader()
                        elif reader.fieldnames != fields:
                            die("Canonical result schemas differ across tasks; aggregate refused.")
                        for row in reader:
                            writer.writerow({"hpc_task": task["id"], "hpc_batch": task["batch"], **row})
                            count += 1
                            if "start" in row and "end" in row:
                                length_sum += int(row["end"]) - int(row["start"])
                batch = next((b for b in plan["batches"] if b["batch"] == task["batch"]), {})
                counts.append({"task": task["id"], "batch": task["batch"], "chromosomes": ",".join(task["chromosomes"]),
                               "batch_cohort_samples": batch.get("scored_n"), "segments": count,
                               "sum_segment_bp": length_sum, "state": "complete", "output": task["output"]})
        with (stage / "task_summary.tsv").open("w") as f:
            w = csv.DictWriter(f, fieldnames=list(counts[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(counts)
        # Aggregate counts are not ancestry proportions: overlapping references,
        # diploid genotype support, callable lengths and batch-specific frequency
        # models require separate scientific summaries.
        (stage / "NOTES.txt").write_text(
            "All task receipts and final-output hashes were checked.\n"
            "A header-only successful task contributes zero calls; a failed/missing task blocks completion.\n"
            "hpc_batch identifies an independent IBDmix frequency-estimation cohort.\n"
            "batch_cohort_samples is before chromosome-specific selection (for example male-only nonPAR X); it is not a denominator.\n"
            "sum_segment_bp sums rows, including overlaps between archaic references; it is not genome coverage.\n"
            "No pooled percentage or unweighted mean of batch percentages is calculated.\n"
            "IBDmix unphased-genotype segments do not directly measure archaic allele dosage.\n")
        for file in stage.iterdir():
            os.replace(file, run / "summary" / file.name)
        atomic_json(run / "summary" / "COMPLETE.json", {"plan_id": plan["plan_id"], "finished": now(),
                    "tasks": len(tasks), "rows": sum(x["segments"] for x in counts),
                    "outputs": {p.name: sha(p) for p in (run / "summary").iterdir() if p.is_file() and p.name != "COMPLETE.json"}})
    finally:
        shutil.rmtree(stage)


def execute(args) -> None:
    run = Path(args.run_dir).resolve()
    plan = load_plan(run)
    require_scheduler(plan, args.stage)
    if args.stage == "prepare":
        with locked(run / "state" / "prepare.lock"):
            (run / "state" / "prepare.json").unlink(missing_ok=True)
            prepare(plan)
    elif args.stage == "worker":
        selection = read_json(Path(args.selection))
        index = int(os.environ.get("LSB_JOBINDEX", "0"))
        if not 1 <= index <= len(selection):
            die("Invalid LSB_JOBINDEX for this array.")
        task_id = selection[index - 1]
        task = next(t for t in plan["tasks"] if t["id"] == task_id)
        execute_task(plan, task)
    elif args.stage == "trace":
        if not plan.get("trace_task"):
            die("No dependent TRACE stage was planned.")
        assemble_arg(plan)
        execute_task(plan, plan["trace_task"], assembled_checked=True)
    else:
        prerequisite(plan)
        with locked(run / "state" / "merge.lock"):
            (run / "summary" / "COMPLETE.json").unlink(missing_ok=True)
            merge(plan)


def bsub_args(plan: dict, stage: str, name: str, dependency: str | None, count: int = 0, cap: int = 0) -> list[str]:
    run = Path(plan["run_dir"])
    slots = 1 if stage in {"prepare", "merge"} else plan["slots"]
    if plan["queue"] == "short":
        slots = max(40, slots)
    memory = (4 if stage in {"prepare", "merge"} else plan["memory_gb"]) * 1024
    use_gpu = plan.get("gpu", False) and stage == "worker"
    queue = plan.get("cpu_queue", plan["queue"]) if plan.get("gpu", False) and not use_gpu else plan["queue"]
    if not use_gpu:
        hours, minutes = map(int, plan["walltime"].split(':'))
        queue, slots = cpu_slots(queue, slots, memory / 1024, hours + minutes / 60)
    reservation = math.ceil(memory / slots) if plan["rusage_memory_mode"] == "per-slot" else memory
    reservation = str(reservation) + {"host": "MB/host", "per-slot": "MB/task", "per-job": "MB/job"}[plan["rusage_memory_mode"]]
    job_name = name + (f"[1-{count}]%{cap}" if count else "")
    args = ["bsub", "-q", queue, "-n", str(slots), "-R", f"span[hosts=1] rusage[mem={reservation}]",
            "-M", str(memory // 1024) + "GB",
            "-W", plan["walltime"], "-J", job_name,
            "-o", str(run / "logs" / (name + ".%J.%I.out")),
            "-e", str(run / "logs" / (name + ".%J.%I.err"))]
    if use_gpu:
        args += ["-gpu", "num=1/host:mode=exclusive_process:j_exclusive=yes"]
    if dependency:
        args += ["-w", dependency]
    return args


def submit_one(plan: dict, stage: str, name: str, dependency: str | None, selection: list[int] | None,
               cap: int, journal, attempt: str) -> str:
    run = Path(plan["run_dir"])
    extra = []
    if selection is not None:
        selected = run / "jobs" / (name + ".indices.json")
        atomic_json(selected, selection)
        extra = [str(selected)]
    args = bsub_args(plan, stage, name, dependency, len(selection) if selection is not None else 0, cap)
    script = (run / "jobs" / (stage + ".cmd")).read_text()
    # Pass array-selection path as a literal argument in the submitted script.
    if extra:
        script = script.replace(' "$@"\n', " " + " ".join(shlex.quote(x) for x in extra) + "\n")
    event = {"time": now(), "attempt": attempt, "stage": stage, "name": name, "selection": selection,
             "command": args, "state": "submitting"}
    journal.write(json.dumps(event) + "\n")
    journal.flush()
    os.fsync(journal.fileno())
    result = subprocess.run(args, input=script, text=True, capture_output=True)
    response = result.stdout + result.stderr
    match = re.search(r"Job\s*<(\d+)>", response)
    if result.returncode != 0 or not match:
        event.update(state="uncertain" if result.returncode == 0 else "rejected", response=response, exit_code=result.returncode)
        journal.write(json.dumps(event) + "\n")
        journal.flush()
        die(f"bsub failed or did not return a parseable job ID. Earlier jobs remain recorded; inspect bjobs before retry. {response.strip()}")
    job = match.group(1)
    event.update(state="submitted", job_id=job, response=response.strip())
    journal.write(json.dumps(event) + "\n")
    journal.flush()
    os.fsync(journal.fileno())
    print(response.strip(), flush=True)
    return job


def prior_submissions(run: Path) -> list[dict]:
    path = run / "submissions.jsonl"
    return [json.loads(line) for line in path.read_text().splitlines()] if path.is_file() else []


def assert_inactive(events: list[dict]) -> None:
    # A successful bsub followed by connection loss can be ambiguous. Never
    # guess that no job exists from a missing parseable identifier.
    if any(x["state"] in {"submitting", "uncertain"} and not any(y.get("name") == x["name"] and y["state"] in {"submitted", "rejected"} for y in events) for x in events):
        die("An earlier submission has an uncertain outcome; inspect bjobs and resolve the journal before retrying.")
    for event in events:
        if event["state"] != "submitted":
            continue
        # Preparation receipts are shared by the whole plan. Even an array
        # disjoint from the retry indices must end before resetting them.
        result = subprocess.run(["bjobs", "-a", "-noheader", "-o", "stat", event["job_id"]], text=True, capture_output=True)
        if result.returncode and not re.search(r"not found|No unfinished job|No matching job", result.stderr + result.stdout, re.I):
            die("Cannot establish earlier job state: " + (result.stderr or result.stdout).strip())
        states = set(result.stdout.split())
        if states & LIVE_STATES:
            die(f"Earlier job {event['job_id']} is still active ({','.join(sorted(states))}); wait or cancel it with bkill before retry.")


def parse_indices(value: str, total: int) -> list[int]:
    if value == "none":
        return []
    ids = set()
    for token in value.split(","):
        if re.fullmatch(r"\d+-\d+", token):
            first, last = map(int, token.split("-"))
            if first > last:
                die("Retry index ranges must be increasing.")
            ids.update(range(first, last + 1))
        elif token.isdigit():
            ids.add(int(token))
        else:
            die("Retry indices use 1,3,7-12 syntax.")
    if not ids or min(ids) < 1 or max(ids) > total:
        die("Retry task index is outside the plan.")
    return sorted(ids)


def submit(args, retry=False) -> None:
    run = Path(args.run_dir).resolve()
    plan = load_plan(run)
    with locked(run / "submission.lock"):
        previous = prior_submissions(run)
        if previous and not retry:
            die("This plan was already submitted (possibly partially). Inspect status; use retry for explicit task indices.")
        if not shutil.which("bsub"):
            die("bsub is not on PATH; run submit from the HPC login shell.")
        selected = parse_indices(args.indices, len(plan["tasks"])) if retry else [t["id"] for t in plan["tasks"]]
        if retry:
            assert_inactive(previous)
            for i in selected:
                rec = task_receipt(plan, plan["tasks"][i - 1])
                if rec.is_file() and read_json(rec).get("state") == "complete":
                    die(f"Task {i} is complete. Do not include it in a failed-task retry.")
        verify_plan_files(plan)
        (run / "state" / "prepare.json").unlink(missing_ok=True)
        (run / "summary" / "COMPLETE.json").unlink(missing_ok=True)
        attempt = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d%H%M%S") + str(os.getpid())
        prefix = "gu." + plan["plan_id"][:10] + "." + attempt
        with (run / "submissions.jsonl").open("a") as journal:
            prep = submit_one(plan, "prepare", prefix + ".prep", None, None, 0, journal, attempt)
            chunks = [selected[i:i + plan["array_size"]] for i in range(0, len(selected), plan["array_size"])]
            ids = []
            for index, chunk in enumerate(chunks):
                cap = plan["max_running"] // len(chunks) + (index < plan["max_running"] % len(chunks))
                ids.append(submit_one(plan, "worker", prefix + f".array{index + 1}", f"ended({prep})", chunk, min(cap, len(chunk)), journal, attempt))
            dep = " && ".join(f"ended({job})" for job in ids) if ids else f"ended({prep})"
            if plan.get("trace_task"):
                trace = submit_one(plan, "trace", prefix + ".trace", dep, None, 0, journal, attempt)
                dep = f"ended({trace})"
            submit_one(plan, "merge", prefix + ".merge", dep, None, 0, journal, attempt)
        print("Submitted. Exact job IDs and commands: " + str(run / "submissions.jsonl"))


def status(args) -> None:
    run = Path(args.run_dir).resolve()
    plan = load_plan(run)
    print("task\tkind\tchromosomes\tbatch\tstate\trows\toutput")
    for task in plan["tasks"] + ([plan["trace_task"]] if plan.get("trace_task") else []):
        p = task_receipt(plan, task)
        r = read_json(p) if p.is_file() else {}
        print("\t".join(map(str, [task["id"], task["kind"], ",".join(task["chromosomes"]), task["batch"],
                                     r.get("state", "not_started"), sum(a.get("rows", 0) for a in r.get("artifacts", [])), task["output"]])))
    print("Aggregate: " + ("COMPLETE" if (run / "summary" / "COMPLETE.json").is_file() else "not complete"))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, epilog="Shell entry point: --date-scratch YYYY-MM-DD (default: 2026-11-03).")
    sub = p.add_subparsers(dest="action", required=True)
    q = sub.add_parser("plan", help="Validate input paths and write an immutable plan; never submit jobs")
    q.add_argument("method", choices=["ibdmix", "trace", "phyml", "arg", "arg-trace"])
    q.add_argument("--run-dir", help="Override plan output; IBDmix defaults beside input genotypes in BUILD/ibdmix")
    q.add_argument("--env-file", required=True)
    q.add_argument("--queue", help="Default: GU_QUEUE, or GU_GPU_QUEUE for --GPU TRUE")
    q.add_argument("--GPU", choices=["TRUE", "FALSE"],
                   help="IBDmix CUDA caller only; validate against CPU on the same cohort before production")
    q.add_argument("--threads", type=positive, default=4)
    q.add_argument("--memory-gb", type=positive, default=32)
    q.add_argument("--walltime", default=os.getenv("HPC_WALLTIME", "120:00"))
    q.add_argument("--phyml-timeout-seconds", type=positive, default=7200,
                   help="Per-tree wall clock budget; independent of the LSF walltime")
    q.add_argument("--phyml-retry-timeouts", action="store_true",
                   help="Retry recorded PhyML timeouts; complete compatible trees are still reused")
    q.add_argument("--max-running", type=positive, default=64)
    q.add_argument("--array-size", type=positive, default=500)
    q.add_argument("--rusage-memory-mode", choices=["per-job", "per-slot", "host"], default="host",
                   help="Explicit host/job/task reservation; default host for single-node jobs")
    q.add_argument("--grch", choices=["37", "38"], default="37")
    q.add_argument("--chr", default="1-22,X")
    q.add_argument("--target", default="1kg")
    q.add_argument("--dir-target", required=True, help="Genotype file prefix ending chr, not the containing directory")
    q.add_argument("--sample-panel")
    q.add_argument("--psam", help="PSAM to split into deterministic IBDmix cohorts")
    q.add_argument("--batch-size", type=positive, default=10000)
    q.add_argument("--keep-males", action="store_true")
    q.add_argument("--archaic-path")
    q.add_argument("--archaic-gen", default="Altai,Chagyr,Denisova,Denisova25,Vindija")
    q.add_argument("--loci")
    q.add_argument("--dir-arg", dest="arg_dir")
    q.add_argument("--dir-gen")
    q.add_argument("--dir-map", help="ARG recombination map directory; default /work/sph-huangj/data/gen/hm3/GRCH/gmap")
    q.add_argument("--arg-method", choices=["threads", "needle", "tsinfer"], default="threads")
    for name in ("submit", "status"):
        q = sub.add_parser(name)
        q.add_argument("run_dir")
    q = sub.add_parser("retry", help="Submit explicit failed task indices after previous jobs have ended")
    q.add_argument("run_dir")
    q.add_argument("--indices", required=True)
    q = sub.add_parser("_execute", help="Internal allocated-job entry; refuses login-node execution")
    q.add_argument("run_dir")
    q.add_argument("stage", choices=["prepare", "worker", "trace", "merge"])
    q.add_argument("selection", nargs="?")
    return p


def main():
    args = parser().parse_args()
    try:
        {"plan": plan_run, "submit": submit, "retry": lambda a: submit(a, True),
         "status": status, "_execute": execute}[args.action](args)
    except (ValueError, OSError, subprocess.CalledProcessError, KeyError, StopIteration) as e:
        print("ERROR: " + str(e), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
