"""Common explicit LSF submission for downstream modules; no shell interpolation."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

from common import CODE_ROOT, atomic_json, load_config, resolve_mode, sha256

MODULES = {
    "plant_validate": (32, 256, 48),
    "group_summary": (2, 16, 4),
    "human_validate": (32, 128, 24),
    "sample_diagnose": (2, 16, 24),
    "species_refine": (32, 256, 48),
    "sylph_profile": (32, 128, 24),
    "split_summary": (2, 16, 24),
}


def dependency(ids, mode="done"):
    parts = [x.strip() for x in ids.split(",") if x.strip()]
    if any(not re.fullmatch(r"[1-9][0-9]*(?:\[[1-9][0-9]*\])?", x) for x in parts):
        raise ValueError("--wait-for requires comma-separated job IDs, optionally JOBID[INDEX]")
    return " && ".join(f"{mode}({x})" for x in parts)


PATH_OPTIONS = {'--out','--output-dir','--representatives','--manifest','--input-root',
                '--compare-root','--reference','--bwa-mem2','--assignments','--summary',
                '--counts','--r1','--r2','--minimap2','--sylph-candidates','--db','--sylph'}


def frozen_args(argv, snapshot, cwd):
    result=[];i=0
    while i<len(argv):
        arg=argv[i];key,sep,value=arg.partition('=')
        if key in ('--config','--root'):
            i+=1 if sep else 2;continue
        if key in PATH_OPTIONS:
            if not sep:i+=1;value=argv[i]
            if key in ('--bwa-mem2','--minimap2','--sylph') and '/' not in value:
                value=shutil.which(value) or value
            result.extend([key,str((cwd/Path(value).expanduser()).resolve())])
        elif key=='--target-fastq':
            result.append(key)
            if sep:result.append(str((cwd/Path(value).expanduser()).resolve()))
            while i+1<len(argv) and not argv[i+1].startswith('--'):
                i+=1;result.append(str((cwd/Path(argv[i]).expanduser()).resolve()))
        else:result.append(arg)
        i+=1
    # Place global options before sylph's subcommand as well.
    return ['--config',str(snapshot),*result]


def main():
    module, *argv = sys.argv[1:]
    if module not in MODULES:
        raise ValueError(f"Unknown QC module {module}")
    p = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    p.add_argument("--submit", choices=("hpc", "local", "auto"), default=None)
    p.add_argument("--plan-only", action="store_true")
    p.add_argument("--wait-for", default="")
    p.add_argument("--wait-mode", choices=("done", "ended"), default="done")
    p.add_argument("--queue", default="38")
    p.add_argument("--cores", type=int, default=MODULES[module][0])
    p.add_argument("--mem-gib", type=int, default=MODULES[module][1])
    p.add_argument("--hours", type=int, default=MODULES[module][2])
    opts, worker_args = p.parse_known_args(argv)
    worker = [sys.executable, str(CODE_ROOT / "f" / (module + ".py")), *worker_args]
    mode = resolve_mode(opts.submit)
    os.environ["submit"] = os.environ["SUBMIT"] = mode
    if min(opts.cores, opts.mem_gib, opts.hours) < 1:
        p.error("Resource limits must be positive")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", opts.queue):
        p.error("Invalid queue name")
    dep = dependency(opts.wait_for, opts.wait_mode)
    # Module parsing happens inside the compute job; inputs may not exist until
    # a dependency finishes. Validate the CLI now without touching those inputs.
    config_parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    config_parser.add_argument("--config"); config_parser.add_argument("--root")
    config_parser.add_argument("--threads", type=int)
    known, _ = config_parser.parse_known_args(worker_args)
    if known.threads is not None and not 1 <= known.threads <= opts.cores:
        p.error("--threads must be positive and no larger than the submitted --cores")
    cfg = load_config(known.config, known.root)
    from source_release import freeze_source
    release = freeze_source(cfg)
    worker[1] = str(release / "f" / (module + ".py"))
    default_threads = {"human_validate": int(cfg.get("human_host_threads", 24)),
                       "species_refine": 32, "plant_validate": 32, "sylph_profile": 32}
    if module in default_threads and known.threads is None:
        worker += ["--threads", str(min(opts.cores, default_threads[module]))]
    local = mode == "local"
    if local and dep:
        p.error("Local mode cannot wait on LSF job IDs; run dependencies first")
    if local and module != "split_summary" and "--allow-local" not in worker:
        worker.append("--allow-local")
    from storage_layout import new_run
    folder = new_run(Path(cfg["work_dir"]) / "jobs" / "qc" / module)
    folder=folder.resolve()
    snapshot=folder/'effective-config.json'
    atomic_json(snapshot,cfg)
    config_digest=sha256(snapshot)
    snapshot.chmod(0o444)
    worker=[sys.executable,str(release/'f'/(module+'.py')),
            *frozen_args(worker[2:],snapshot,Path.cwd())]
    logs = Path(cfg["work_dir"]) / "logs" / "qc" / module / folder.name
    logs.mkdir(parents=True, exist_ok=True)
    path = folder / "run.cmd"
    body = ["#!/usr/bin/env bash", "set -euo pipefail",
            "source " + shlex.quote(str(release / "f" / "env.sh")),
            "cd " + shlex.quote(str(Path.cwd())),
            shlex.join([sys.executable, str(release/'f'/'qc_submit.py'),
                        '--verify-config',str(snapshot),config_digest,str(release)]),
            "exec " + shlex.join(worker)]
    path.write_text("\n".join(body) + "\n"); path.chmod(0o755)
    command = ["bsub", "-J", "metag_" + module, "-q", opts.queue, "-n", str(opts.cores),
               "-W", f"{opts.hours}:00", "-R",
               f"span[hosts=1] rusage[mem={opts.mem_gib}GB/host]",
               "-oo", str(logs / "run.out"),
               "-eo", str(logs / "run.err")]
    if dep:
        command += ["-w", dep]
    receipt = {"command": command, "script": str(path), "worker": worker,
               "dependency": dep, "status": "planned", "source_release": str(release),
               "effective_config": str(snapshot), "config_sha256": config_digest, "submission_cwd": str(Path.cwd())}
    if not opts.plan_only:
        if local:
            result = subprocess.run(["bash", str(path)])
            receipt.update(status="complete" if result.returncode == 0 else "failed", exit_code=result.returncode, mode="local")
            atomic_json(folder / "submission.json", receipt)
            raise SystemExit(result.returncode)
        with path.open() as handle:
            run = subprocess.run(command, stdin=handle, text=True, capture_output=True)
        print(run.stdout, end=""); print(run.stderr, end="", file=sys.stderr)
        if run.returncode:
            raise RuntimeError("LSF submission failed")
        match = re.search(r"Job <(\d+)>", run.stdout)
        if not match:
            raise RuntimeError("LSF accepted submission but job ID was not recognized; inspect output before retrying")
        receipt.update(status="submitted", job_id=match[1])
    atomic_json(folder / "submission.json", receipt)
    print("Manual submit: " + ("bash " + shlex.quote(str(path)) if local else shlex.join(command) + " < " + shlex.quote(str(path))))
    print("Receipt: " + str(folder / "submission.json"))


if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1]=='--verify-config':
        _,_,snapshot,digest,release=sys.argv
        if sha256(snapshot)!=digest:raise RuntimeError('Frozen QC configuration checksum mismatch')
        root=Path(release).parent
        receipt=json.loads((root/'release.json').read_text())
        for relative,expected in receipt['files'].items():
            if sha256(root/relative)!=expected:raise RuntimeError('Frozen QC source checksum mismatch: '+relative)
    else:
        main()
