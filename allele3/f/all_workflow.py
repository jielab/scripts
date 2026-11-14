#!/usr/bin/env python3
"""One allocated CPU job for downstream allele3 analyses using prepared inputs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config/workflow.json"


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                   help="Ordered scientific steps; default config/workflow.json")
    p.add_argument("--out", type=Path,
                   help="Analysis root containing all module outputs")
    p.add_argument("--python", default=os.environ.get("ALLELE_PYTHON", sys.executable))
    p.add_argument("--plan-only", action="store_true")
    return p


def load_steps(config, out, python, cores, memory_gib):
    """Expand named paths without a shell, imports, or genomic input reads."""
    from workflow_paths import configured_path
    variables = dict(os.environ)
    variables.update(
        ALLELE_ROOT=str(ROOT), ALLELE_ANALYSIS_ROOT=str(Path(out).absolute()),
        DIR_DATA=os.environ.get("DIR_DATA", configured_path("DIR_BIG")),
        DIR_FILES=os.environ.get("DIR_FILES", configured_path("DIR_FILES")),
        ALLELE_PYTHON=str(python), CORES=str(cores),
        SCAN_JOBS=str(min(4, cores)), PCA_THREADS=str(min(4, cores)),
        PCA_MEMORY_MB=str(max(512, min(65536, int(memory_gib * 768)))),
        CENSUS_JOBS=str(min(8, cores)), CENSUS_CORES=str(cores))
    raw = json.loads(Path(config).read_text())
    if set(raw) != {"steps"} or not isinstance(raw["steps"], list) or not raw["steps"]:
        raise ValueError("Workflow config must contain a nonempty steps array")
    token = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")
    def expand(value):
        if not isinstance(value, str) or any(c in value for c in ("\n", "\r", "\0")):
            raise ValueError("Every workflow argument must be a single text value")
        def replace(match):
            key = match.group(1)
            if key not in variables:
                raise ValueError("Undefined workflow variable: " + key)
            return variables[key]
        return token.sub(replace, value)
    result = []
    for item in raw["steps"]:
        if set(item) != {"name", "argv"} or not re.fullmatch(r"[a-z][a-z0-9_]*", item["name"]):
            raise ValueError("Each step needs a unique lower-case name and argv")
        if not isinstance(item["argv"], list) or not item["argv"]:
            raise ValueError("Empty step argv: " + item["name"])
        args = [expand(value) for value in item["argv"]]
        if any(a in {"--submit", "--write-cmd", "--wait", "--job-dir", "--replace"} for a in args):
            raise ValueError("Set execution and replacement options on allele3.sh all, not individual config steps")
        result.append(dict(name=item["name"], argv=args))
    if len({s["name"] for s in result}) != len(result):
        raise ValueError("Duplicate workflow step name")
    return result


def check_prerequisites(out):
    """Check only required input paths, without size/hash or receipt audits."""
    from analysis_integrity import resolve_result_input
    out = Path(out).expanduser().absolute()
    inputs = {
        "1.pca3D.sh": [out / "pca3d" / name for name in (
            "tables/sample_metadata.tsv", "coordinates/all.sscore",
            "basis/ref.acount", "basis/ref.eigenvec.allele", "basis/ld.prune.in")],
        "3.gnomad_scan.sh": [Path(resolve_result_input(out / "gnomad_scan")) / "run.json"],
    }
    # Census is independent: downstream analyses do not consume its reports.
    failures = []
    for script, paths in inputs.items():
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            failures.append(script)
            print(f"[allele3 input missing] {script}: " + "; ".join(missing),
                  file=sys.stderr, flush=True)
            print("Run the independent " + script + " first, then rerun allele3.sh all.",
                  file=sys.stderr, flush=True)
    if failures:
        raise ValueError("Required inputs missing; no downstream analysis started.")
    print("[allele3 inputs OK] Required PCA files and scan manifest exist.",
          file=sys.stderr, flush=True)


def build_plan(argv, mode, python):
    from batch_workflows import base_plan, add_job, resources, env_paths
    a = parser().parse_args(argv)
    out = (a.out or Path(env_paths()[2])).expanduser().absolute()
    check_prerequisites(out)
    resource = resources(40 if mode == "hpc" else min(16, os.cpu_count() or 1),
                         160 if mode == "hpc" else 32, 120)
    steps = load_steps(a.config, out, python, resource["cores"], resource["memory_gib"])
    plan = base_plan("all", mode)
    plan["max_running"] = 1
    plan["local_memory_gib"] = resource["memory_gib"]
    done = out / "qc/all.json"
    add_job(plan, "all", [python, __file__, "worker", "--out", str(out),
                         "--steps-json", json.dumps(steps, ensure_ascii=False)],
            python, resource=resource, outputs=[done])
    plan["jobs"][0]["env"].update(
        ALLELE_ANALYSIS_ROOT=str(out), PCA3D_ANALYSIS_ROOT=str(out),
        ALLELE_SCAN_STORE=str(out / "gnomad_scan"),
        ALLELE_METADATA=str(out / "pca3d/tables/sample_metadata.tsv"),
        ALLELE_COORDINATES=str(out / "pca3d/coordinates/all.sscore"),
        PCA3D_METADATA=str(out / "pca3d/tables/sample_metadata.tsv"),
        PCA3D_EIGENVEC=str(out / "pca3d/coordinates/all.sscore"),
        PCA3D_REUSE_BASIS=str(out / "pca3d/basis"))
    plan["provenance"].update(execution_order="sequential_success",
        steps=[step["name"] for step in steps],
        scheduling="one CPU allocation and one all.cmd; no nested submissions")
    return plan, [], out


def run_steps(steps, out, *, runner=None):
    """Stop immediately on failure and record unstarted steps as blocked."""
    from analysis_integrity import atomic_json
    from batch_workflows import pipeline_file_lock, pipeline_worker_guard
    pipeline_worker_guard()
    out = Path(out).absolute()
    check_prerequisites(out)
    runner = runner or subprocess.run
    record = dict(status="running", steps=[], started=time.time())
    path = out / "qc/all.json"
    environment = dict(os.environ, WORKFLOW_WORKER="1", ALLELE_BATCH_WORKER="1",
                       ALLELE_IN_COMPUTE="1", ALLELE_ANALYSIS_ROOT=str(out),
                       CUDA_VISIBLE_DEVICES="", WORKFLOW_DEVICE="cpu")
    with pipeline_file_lock(out / "qc/all.lock", nonblocking=True):
        atomic_json(path, record)
        for index, step in enumerate(steps):
            print(f"[allele3 {index + 1}/{len(steps)}] {step['name']}: {shlex.join(step['argv'])}", flush=True)
            row = dict(name=step["name"], status="running", started=time.time())
            record["steps"].append(row)
            atomic_json(path, record)
            try:
                code = runner(step["argv"], cwd=ROOT, env=environment).returncode
            except OSError as exc:
                code = 1
                row["error"] = str(exc)
            row.update(status="complete" if code == 0 else "failed", returncode=code, finished=time.time())
            if code:
                record["steps"].extend(dict(name=s["name"], status="blocked", blocked_by=step["name"])
                                       for s in steps[index + 1:])
                record.update(status="failed", finished=time.time())
                atomic_json(path, record)
                return code
            atomic_json(path, record)
        record.update(status="complete", finished=time.time())
        atomic_json(path, record)
        return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["worker"]:
        p = argparse.ArgumentParser(description=__doc__)
        p.add_argument("--out", required=True)
        p.add_argument("--steps-json", required=True)
        a = p.parse_args(argv[1:])
        return run_steps(json.loads(a.steps_json), a.out)
    parser().parse_args(argv)
    raise ValueError("Use bash allele3.sh all to create or submit the one-job workflow")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        raise SystemExit("allele3-all: " + str(exc))
