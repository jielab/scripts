#!/usr/bin/env python3
"""Reviewable LSF DAG for parallel reference categories and durable native builds."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid

from hpc_common import (CODE_ROOT, atomic_json, effective_workers, engine_snapshot,
                        fingerprint, load_config, locked, require_compute, utcnow)
import hpc_cli as jobs
from hpc_lsf import write_cmd
import partitioned_build as build


def _read(path):
    return json.loads(Path(path).read_text())


def _command(folder, action, *extra):
    return jobs._command(folder / "config.json", action, folder / "jobs.json", *extra)


def _require_partitioned(cfg):
    if cfg.get("database_layout") != "partitioned":
        raise ValueError("This library workflow requires database_layout: partitioned")
    if cfg.get("integrity_seal"):
        raise ValueError("Old single-index integrity seals cannot be used for the new partitioned layout; it maintains separate SHA receipts")
    for key in ("reference_raw", "reads_raw", "metadata_dir", "library_dir", "work_dir", "output_dir", "tools_dir",
                "legacy_library_dir", "legacy_metadata_dir", "parallel_reference_dir"):
        if cfg.get(key):
            path = Path(cfg[key]).expanduser().resolve()
            if (path == Path("/scratch") or Path("/scratch") in path.parents) and cfg.get("allow_scratch") is not True:
                raise ValueError(f"{key} points into auto-purged /scratch: {path}. Copy it to persistent project storage and update hpc.yml before generating resumable production jobs.")


def generate_setup(cfg, metadata_source=None, download=False):
    _require_partitioned(cfg)
    source = str(Path(metadata_source).expanduser().resolve()) if metadata_source else None
    if download and not cfg.get("legacy_metadata_dir"):
        # Network is confined to the submit host, where the operator explicitly
        # requested it; compute nodes only consume the downloaded metadata.
        from hpc_reference import _download
        target = Path(cfg["metadata_dir"]) / "import"
        target.mkdir(parents=True, exist_ok=True)
        names = [("taxdump.tar.gz", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")]
        for namespace in ("refseq", "genbank"):
            for suffix in ("", "_historical"):
                name = f"assembly_summary_{namespace}{suffix}.txt"
                names.append((name, "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/" + name))
        for name, url in names:
            if not (target / name).is_file():
                _download(url, target / name)
        source = str(target)
    folder, manifest = jobs._new_plan(cfg, "parallel_setup", {"metadata_source": source})
    path = folder / "00.setup.cmd"
    write_cmd(cfg, path, "setup", _command(folder, "parallel-setup-worker"), "metag_setup")
    manifest["commands"] = [{"role": "setup", "path": str(path)}]
    return jobs._publish_plan(cfg, folder, manifest)


def generate_selftest(cfg):
    _require_partitioned(cfg)
    cfg = dict(cfg, setup_cores=4, setup_memory_gib=16, setup_hours=2)
    destination = Path(cfg["work_dir"]) / "selftests" / ("run_" + uuid.uuid4().hex[:12])
    folder, manifest = jobs._new_plan(cfg, "parallel_selftest", {"destination": str(destination)})
    path = folder / "00.selftest.cmd"
    write_cmd(cfg, path, "setup", _command(folder, "parallel-selftest-worker"), "metag_selftest")
    manifest["commands"] = [{"role": "selftest", "path": str(path)}]
    return jobs._publish_plan(cfg, folder, manifest)


def _sample_context(cfg, datasets, only, ready_only):
    from hpc_query import discover_samples
    samples = discover_samples(cfg, datasets=datasets, only=only, ready_only=ready_only)
    states = [jobs._sample_states(s) for s in samples]
    batches = jobs.balanced_batches(samples, cfg["sample_batch_size"])
    largest = {}
    for index, sample in enumerate(samples):
        size = sum(item["size"] for item in states[index].values())
        if size > largest.get(sample["dataset"], (-1, -1))[0]:
            largest[sample["dataset"]] = (size, index)
    return {"samples": samples, "sample_states": states, "batches": batches,
            "pilot_indices": [largest[k][1] for k in sorted(largest)]}


def _add_split(cfg, folder, manifest, parent):
    path = folder / "2.pilot.cmd"
    write_cmd(cfg, path, "pilot", _command(folder, "parallel-pilot-worker"), "metag_pilot")
    manifest["commands"].append({"role": "pilot", "path": str(path), "depends_on": parent})
    batches = manifest["context"]["split"]["batches"]
    for part, (first, count) in enumerate(jobs._arrays(len(batches), cfg["max_array_size"]), 1):
        path = folder / f"2.split.array.{part:03d}.cmd"
        write_cmd(cfg, path, "classify", _command(folder, "parallel-split-worker",
                  "--batch-id", "__LSB_JOBINDEX__", "--offset", first - 1), f"metag_split_{part}",
                  array=f"1-{count}%{min(count, cfg['classify_concurrency'])}")
        manifest["commands"].append({"role": f"split_{part}", "path": str(path),
                                     "depends_on": "pilot" if part == 1 else f"split_{part-1}"})


def generate_library(cfg, with_split=False, datasets=None, only=None, ready_only=False):
    _require_partitioned(cfg)
    from partitioned_reference import load_plan
    reference = load_plan(cfg)
    tasks = reference["tasks"]
    if not tasks:
        raise RuntimeError("Reference plan has no preparation tasks")
    categories = [name for name in build.CATEGORIES if name in reference["categories"]]
    for name in categories:
        build.resources(cfg, name)
    context = {"reference_identity": reference["identity"], "reference_tasks": len(tasks),
               "categories": categories}
    if with_split:
        cfg = dict(cfg, output_dir=str(Path(cfg["output_dir"]) / "kraken"))
        context["split"] = _sample_context(cfg, datasets, only, ready_only)
    folder, manifest = jobs._new_plan(cfg, "parallel_library", context)
    workers = min(len(tasks), cfg["prepare_concurrency"], cfg["max_array_size"])
    path = folder / "1.prepare.cmd"
    write_cmd(cfg, path, "prepare", _command(folder, "parallel-prepare-worker",
              "--task-id", "__LSB_JOBINDEX__", "--stride", workers),
              "metag_prepare", array=f"1-{workers}%{workers}")
    manifest["commands"].append({"role": "prepare", "path": str(path)})
    path = folder / "1.finalize-preparation.cmd"
    write_cmd(cfg, path, "setup", _command(folder, "parallel-finalize-worker"), "metag_prepare_final")
    manifest["commands"].append({"role": "finalize_preparation", "path": str(path), "depends_on": "prepare"})
    for category in categories:
        path = folder / f"1.build.{category}.cmd"
        allocation = build.resources(cfg, category)
        write_cmd(allocation, path, "build", _command(folder, "parallel-build-worker", "--category", category),
                  "metag_build_" + category)
        manifest["commands"].append({"role": "build_" + category, "path": str(path), "depends_on": "finalize_preparation"})
    path = folder / "1.publish.cmd"
    write_cmd(cfg, path, "setup", _command(folder, "parallel-publish-worker"), "metag_index_publish")
    manifest["commands"].append({"role": "publish", "path": str(path),
                                 "depends_on": ["build_" + category for category in categories]})
    if with_split:
        _add_split(cfg, folder, manifest, "publish")
    result = jobs._publish_plan(cfg, folder, manifest)
    print(f"{len(tasks)} reference tasks, {workers} preparation jobs, {len(categories)} independent category builds")
    print("Builds run concurrently after preparation; publish requires done() for every category.")
    print("A native checkpoint pause exits 75 and blocks downstream jobs. Use submit --retry --cancel-pending to resume after all running jobs have stopped.")
    return result


def generate_split(cfg, datasets=None, only=None, ready_only=False, queue=None):
    _require_partitioned(cfg)
    cfg = dict(cfg, output_dir=str(Path(cfg["output_dir"]) / "kraken"))
    if queue is not None:
        cfg.update(classify_queue=queue, pilot_queue=queue)
    snapshot = build.database_snapshot(cfg)
    effective_workers(cfg, snapshot)
    context = {"database_identity": snapshot["identity"], "database_runtime_states": jobs._database_states(snapshot),
               "split": _sample_context(cfg, datasets, only, ready_only)}
    folder, manifest = jobs._new_plan(cfg, "parallel_split", context)
    _add_split(cfg, folder, manifest, [])
    result = jobs._publish_plan(cfg, folder, manifest)
    print(f"Samples={len(context['split']['samples'])}; batches={len(context['split']['batches'])}; index union is queried in one FASTQ pass")
    return result


def _job_states(ids):
    if not ids:
        return {}
    process = subprocess.run(["bjobs", "-a", "-noheader", "-o", "jobid stat", *map(str, sorted(ids))],
                             capture_output=True, text=True, timeout=30, env={**os.environ, "LC_ALL": "C"})
    if process.returncode:
        raise RuntimeError("Cannot establish prior job states; inspect bjobs/bhist before recovery: " + process.stderr.strip())
    states = {job: set() for job in ids}
    for line in process.stdout.splitlines():
        fields = line.split()
        if len(fields) < 2 or not fields[0].split("[", 1)[0].isdigit():
            raise RuntimeError("Unexpected bjobs format; recovery stopped")
        job = int(fields[0].split("[", 1)[0])
        if job in states:
            states[job].add(fields[1])
    if any(not value for value in states.values()):
        raise RuntimeError("Some previous JOBIDs have expired from bjobs history; verify them with bhist/bacct before a controlled manual resubmission")
    return states


def submit_plan(path, retry=False, cancel_pending=False):
    """Submit numeric success dependencies. Retrying never treats a pause as DONE."""
    from hpc_lsf import submit_cmd, SubmissionError
    path = Path(path).resolve()
    initial = _read(path)
    cfg = load_config(initial["config"])
    manifest = jobs._validate_job_plan(cfg, path)
    target = path.parent / "submission.json"
    with locked(path.parent / "submission.lock"):
        history = _read(target) if target.exists() else {"schema": 1, "attempts": []}
        if cancel_pending and not retry:
            raise ValueError("--cancel-pending is only valid with --retry")
        if history["attempts"]:
            if not retry:
                raise RuntimeError("This plan already has submission receipts. Inspect status, then use --retry for recovery.")
            for attempt in history["attempts"]:
                for job in attempt.get("jobs", []):
                    if job.get("status") in ("intent", "unknown") and not job.get("job_id"):
                        raise RuntimeError("Prior bsub acceptance is unknown; inspect scheduler history. Automatic duplicate submission refused.")
            previous = history["attempts"][-1]
            ids = {job["job_id"] for job in previous.get("jobs", []) if isinstance(job.get("job_id"), int)}
            states = _job_states(ids)
            running = {job: sorted(value) for job, value in states.items() if value - {"DONE", "EXIT", "PEND", "PSUSP"}}
            if running:
                raise RuntimeError(f"Jobs are still running/suspended or have an uncertain state; wait for them before retry: {running}")
            pending = sorted(job for job, value in states.items() if value & {"PEND", "PSUSP"})
            if pending and not cancel_pending:
                raise RuntimeError(f"Recorded downstream jobs still wait on old IDs: {pending}. Use --retry --cancel-pending to cancel this plan's pending jobs and rebuild dependencies.")
            if pending:
                result = subprocess.run(["bkill", *map(str, pending)], capture_output=True, text=True, timeout=30)
                previous["pending_cancellation"] = {"job_ids": pending, "at": utcnow(),
                    "returncode": result.returncode, "stdout": result.stdout, "stderr": result.stderr}
                atomic_json(target, history)
                for _ in range(20):
                    states = _job_states(ids)
                    if all(not (value - {"DONE", "EXIT"}) for value in states.values()):
                        break
                    time.sleep(1)
                else:
                    raise RuntimeError("Some canceled jobs have not reached terminal state; wait and retry without submitting duplicates")
            previous["recovery_verified_terminal_at"] = utcnow()
        attempt = {"started_at": utcnow(), "status": "submitting", "jobs": []}
        history["attempts"].append(attempt)
        atomic_json(target, history)
        by_role = {}
        try:
            for command in manifest["commands"]:
                parents = command.get("depends_on", [])
                parents = [parents] if isinstance(parents, str) else parents
                if any(parent not in by_role for parent in parents):
                    raise RuntimeError("Invalid/missing preceding dependency role")
                condition = command.get("dependency_condition", "done")
                if condition not in ("done", "ended"):
                    raise ValueError("Unsupported dependency condition")
                dependency = " && ".join(f"{condition}({job})" for parent in parents for job in by_role[parent]) or None
                receipt = {"status": "intent", "role": command["role"], "command_path": command["path"], "dependency": dependency}
                attempt["jobs"].append(receipt)
                atomic_json(target, history)
                try:
                    result = submit_cmd(cfg, command["path"], dependency=dependency, expected_sha256=command["sha256"])
                except SubmissionError as exc:
                    receipt.update(exc.result)
                    raise
                receipt.update(result)
                by_role.setdefault(command["role"], []).append(receipt["job_id"])
                atomic_json(target, history)
                print(f"Submitted {command['role']}: {receipt['job_id']}; dependency={dependency or '(none)'}", flush=True)
            attempt["status"] = "submitted"
        except BaseException as exc:
            attempt.update(status="submission_failed", error=str(exc))
            raise
        finally:
            atomic_json(target, history)
    return history


def worker(args, cfg):
    require_compute(args.allow_local)
    plan = jobs._validate_job_plan(cfg, args.job_plan)
    context = plan["context"]
    if args.command == "parallel-selftest-worker":
        if plan["kind"] != "parallel_selftest":
            raise RuntimeError("Selftest requires its selftest plan")
        command = [sys.executable, str(CODE_ROOT / "f/2test/parallel_smoke.py"),
                   "--destination", context["destination"], "--tools-dir", cfg["tools_dir"],
                   "--code-root", str(CODE_ROOT)]
        if args.allow_local:
            command.append("--allow-local")
        return subprocess.run(command, check=False).returncode
    if args.command == "parallel-setup-worker":
        if plan["kind"] != "parallel_setup":
            raise RuntimeError("Setup requires its setup plan")
        from hpc_install import install
        from native_checkpoint.install import install as install_checkpoint
        from native_union.install import install as install_union
        from partitioned_reference import plan as plan_references
        install(cfg)
        install_checkpoint(cfg)
        install_union(cfg)
        if not cfg.get("legacy_metadata_dir"):
            from hpc_reference import metadata
            metadata(cfg, source=context.get("metadata_source"), download=False)
        reference = plan_references(cfg)
        atomic_json(Path(args.job_plan).parent / "setup.complete.json",
                    {"status": "complete", "reference_identity": reference["identity"], "finished_at": utcnow()})
        print("Setup complete. Generate the parallel DAG with ./1.library.sh --with-split --dat-list ...")
        return 0
    if plan["kind"] == "parallel_library":
        from partitioned_reference import load_plan, prepare_task, finalize_preparation, catalog_snapshot
        reference = load_plan(cfg)
        if reference["identity"] != context["reference_identity"]:
            raise RuntimeError("Reference plan changed after job generation")
        if args.command == "parallel-prepare-worker":
            for task_id in jobs._prepare_task_ids(args.task_id, len(reference["tasks"]), args.stride):
                print(f"Preparing/adopting task {task_id}/{len(reference['tasks'])}", flush=True)
                prepare_task(cfg, task_id)
            return 0
        if args.command == "parallel-finalize-worker":
            finalize_preparation(cfg)
            return 0
        if args.command == "parallel-build-worker":
            if args.category not in context["categories"]:
                raise RuntimeError("Category not included in the frozen job plan")
            return build.build_category(cfg, args.category)
        if args.command == "parallel-publish-worker":
            value = build.publish(cfg)
            # Do not publish a query allocation that is too small silently.
            if "split" in context:
                effective_workers(cfg, value)
            return 0
    elif plan["kind"] != "parallel_split":
        raise RuntimeError("Wrong job plan for parallel worker")
    if args.command not in ("parallel-pilot-worker", "parallel-split-worker"):
        raise RuntimeError("Worker action does not match plan")
    from hpc_query import run_sample, run_batch, BatchFailure
    snapshot = build.database_snapshot(cfg)
    if plan["kind"] == "parallel_library":
        prepared = catalog_snapshot(cfg)
        if snapshot["context"]["preparation_identity"] != prepared["identity"]:
            raise RuntimeError("Published index belongs to another preparation")
    else:
        if snapshot["identity"] != context["database_identity"] or jobs._database_states(snapshot) != context["database_runtime_states"]:
            raise RuntimeError("Database changed after split job generation")
    cfg = effective_workers(cfg, snapshot)
    splitting = context["split"]
    pilot = args.command == "parallel-pilot-worker"
    batch_id = None if pilot else args.batch_id + args.offset
    if not pilot and not 1 <= batch_id <= len(splitting["batches"]):
        raise RuntimeError("Batch ID is outside the frozen plan")
    indices = splitting["pilot_indices"] if pilot else splitting["batches"][batch_id - 1]
    samples = [splitting["samples"][i] for i in indices]
    for index in indices:
        if jobs._sample_states(splitting["samples"][index]) != splitting["sample_states"][index]:
            raise RuntimeError("FASTQ input changed after job generation")
    report = Path(args.job_plan).parent / ("pilot.complete.json" if pilot else f"batch.{batch_id:05d}.json")
    try:
        outputs = ([run_sample(cfg, sample, snapshot=snapshot, limit=cfg["pilot_fragments"]) for sample in samples]
                   if pilot else run_batch(cfg, samples, snapshot=snapshot))
    except BatchFailure as exc:
        atomic_json(report, {"status": "failed", "finished_at": utcnow(), "result": exc.summary})
        raise
    atomic_json(report, {"status": "complete", "finished_at": utcnow(), "outputs": outputs,
                        "pilot_fragments": cfg["pilot_fragments"] if pilot else None,
                        "note": "Pilot prefixes test operation, not total throughput or biological sensitivity" if pilot else None})
    return 0


def status(cfg):
    values = {}
    for category in build.CATEGORIES:
        folder = build.index_root(cfg) / category
        item = {}
        for name in ("database.json", "process.json", "resource_estimate.json", "checkpoints/progress.json"):
            path = folder / name
            if path.is_file():
                value = _read(path)
                item[name] = {k: v for k, v in value.items() if k not in ("context", "command", "files")}
        current = folder / "checkpoints/CURRENT"
        if current.is_file():
            _, item["committed_checkpoint"] = build.committed_receipt(folder / "checkpoints")
        values[category] = item or {"status": "not_started"}
    print(json.dumps(values, indent=2))
    print("Only committed_checkpoint.completed_files are durable; progress may be ahead of that. LSF status is inspected with bjobs.")


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--config")
    result.add_argument("--root")
    sub = result.add_subparsers(dest="command", required=True)
    def command(name, help_text):
        obj = sub.add_parser(name, help=help_text)
        obj.add_argument("--config", default=argparse.SUPPRESS)
        obj.add_argument("--root", default=argparse.SUPPRESS)
        return obj
    def submission(obj):
        obj.add_argument("--submit", type=jobs.submit_bool, default=False, metavar="TRUE|FALSE")
    def samples(obj):
        obj.add_argument("--dat-list")
        obj.add_argument("--sample", action="append")
        obj.add_argument("--ready-only", action="store_true")
    obj = command("setup", "Generate setup/compiler/metadata .cmd; default only generates")
    obj.add_argument("--metadata-source")
    obj.add_argument("--download-metadata", action="store_true", help="Download official metadata on this submit host before generating setup")
    obj.add_argument("--offline", action="store_true", help="Explicitly prohibit metadata download")
    submission(obj)
    obj = command("selftest", "Generate an ordinary-node real pause/resume/paired-read pipeline selftest")
    submission(obj)
    obj = command("library", "Prepare/reuse references; parallel resumable category builds; publish catalog")
    obj.add_argument("--with-split", action="store_true", help="Also submit pilot and split jobs with successful-build dependencies")
    samples(obj)
    submission(obj)
    obj = command("split", "Classify and split using a completed partitioned index")
    samples(obj)
    obj.add_argument("--queue", type=jobs.queue_name)
    submission(obj)
    obj = command("submit", "Submit a generated DAG; retry can rebuild stale dependencies")
    obj.add_argument("--manifest", required=True)
    obj.add_argument("--retry", action="store_true")
    obj.add_argument("--cancel-pending", action="store_true", help="On retry, cancel only this recorded plan's pending jobs; running jobs cause a refusal")
    command("doctor", "Inventory LSF and configured persistent filesystems")
    command("status", "Read category progress and durable checkpoint receipts")
    obj = command("audit", "Validate completed partitioned index")
    obj.add_argument("--full", action="store_true", help="Compute job: additionally read/hash every index byte")
    obj.add_argument("--allow-local", action="store_true")
    for name in ("selftest", "setup", "prepare", "finalize", "build", "publish", "pilot", "split"):
        obj = command("parallel-" + name + "-worker", "Generated compute entry")
        obj.add_argument("--job-plan", required=True)
        obj.add_argument("--allow-local", action="store_true")
        if name == "prepare":
            obj.add_argument("--task-id", type=int, required=True)
            obj.add_argument("--stride", type=int, default=0)
        if name == "build":
            obj.add_argument("--category", choices=build.CATEGORIES, required=True)
        if name == "split":
            obj.add_argument("--batch-id", type=int, required=True)
            obj.add_argument("--offset", type=int, default=0)
    return result


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # The existing benchmark implementation remains available for a separately
    # completed single index. It cannot convert/merge partitioned final taxids.
    legacy = {"compare", "compare-worker", "compare-summary-worker", "convert-worker", "pilot-worker", "split-worker"}
    if any(value in legacy for value in argv[:1]):
        return jobs.main(argv)
    args = parser().parse_args(argv or ["setup"])
    cfg = load_config(args.config, args.root)
    if args.command == "setup":
        if args.offline and args.download_metadata:
            raise ValueError("--offline and --download-metadata conflict")
        path = generate_setup(cfg, args.metadata_source, args.download_metadata)
    elif args.command == "selftest":
        path = generate_selftest(cfg)
    elif args.command == "library":
        path = generate_library(cfg, args.with_split, args.dat_list, args.sample, args.ready_only)
    elif args.command == "split":
        path = generate_split(cfg, args.dat_list, args.sample, args.ready_only, args.queue)
    elif args.command == "submit":
        submit_plan(args.manifest, args.retry, args.cancel_pending)
        return 0
    elif args.command == "doctor":
        from hpc_lsf import doctor
        report = doctor(cfg)
        atomic_json(Path(cfg["work_dir"]) / "doctor.json", report)
        print(f"Inventory saved: {Path(cfg['work_dir']) / 'doctor.json'}")
        return 0
    elif args.command == "status":
        status(cfg)
        return 0
    elif args.command == "audit":
        if args.full:
            require_compute(args.allow_local)
        value = build.database_snapshot(cfg, full=args.full)
        print(json.dumps({"identity": value["identity"], "categories": value["context"]["categories"],
                          "hash_bytes": sum(value["files"][name]["state"]["size"] for name in value["hashes"]),
                          "full_sha256": args.full}, indent=2))
        return 0
    else:
        return worker(args, cfg)
    if args.submit:
        submit_plan(path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
