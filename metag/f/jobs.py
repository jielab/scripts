#!/usr/bin/env python3
"""Generate reviewable LSF commands; run heavy steps only inside compute jobs."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from common import (split_job_name, CODE_ROOT, atomic_json, code_identity, effective_workers,
                        file_state, fingerprint, load_config, read_config_file, locked, require_compute, sha256, utcnow,
                        ENGINES, engine_resources, engine_snapshot)


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _cfg_identity(cfg):
    return fingerprint({k: v for k, v in cfg.items() if k not in ("config_path", "code_root")})


def _integrity_seal(cfg, expected=None, required_indexes=()):
    """Optional persisted reference baseline; only cheap checks at job boundaries."""
    value = cfg.get("integrity_seal")
    if not value:
        if expected is not None:
            raise RuntimeError("Frozen integrity seal is no longer configured")
        return None
    from audit_integrity import verify_seal
    path = Path(value).expanduser().resolve(strict=True)
    observed = {"path": str(path), "sha256": sha256(path)}
    if expected is not None and observed != expected:
        raise RuntimeError("Integrity seal changed after job generation")
    verified = verify_seal(path)
    for key in ("reference_raw", "metadata_dir", "library_dir", "db_dir", "kunpeng_db_dir"):
        if (Path(verified["paths"][key]).resolve() != Path(cfg[key]).resolve()):
            raise RuntimeError(f"Integrity seal belongs to a different {key}")
    for name in required_indexes:
        index = verified["indexes"].get(name, {})
        if index.get("status") != "complete":
            raise RuntimeError(f"Integrity seal does not cover the completed {name} index; create a new seal")
        if _read(Path(index["directory"]) / "database.json").get("identity") != index.get("identity"):
            raise RuntimeError(f"Integrity seal does not match the current {name} database identity")
    return observed


def _required_sealed_indexes(kind, cfg, context):
    if kind == "split":
        return ("kraken" if cfg.get("classifier_engine", "kraken") == "kraken" else "kunpeng",)
    if kind == "compare":
        return sorted({"kraken" if case["classifier_engine"] == "kraken" else "kunpeng"
                       for case in context["cases"]})
    if kind == "convert":
        return ("kraken",)
    return ()


def _new_plan(cfg, kind, context):
    from source_release import freeze_source
    release = freeze_source(cfg)
    cfg["code_root"] = str(release)
    context = dict(context, source_release=str(release))
    seal = _integrity_seal(cfg, required_indexes=_required_sealed_indexes(kind, cfg, context))
    if seal is not None:
        context["integrity_seal"] = seal
    code = code_identity()
    identity = fingerprint({"schema": 1, "kind": kind, "config": _cfg_identity(cfg),
                            "code": code, "context": context})
    from storage_layout import plan_directory, new_run
    label = kind.removeprefix("parallel_")
    if label == "split":
        label = cfg.get("run_label", "split")
    parent = Path(cfg["work_dir"]) / "jobs" / label
    if kind in ("split", "parallel_split", "compare") or "split" in context:
        folder = new_run(parent)
        atomic_json(folder / "identity.json", {"identity": identity})
    else:
        folder = plan_directory(parent, identity)
    if folder.is_symlink():
        raise RuntimeError(f"Refusing symlink at immutable job plan: {folder}")
    folder.mkdir(parents=True, exist_ok=True)
    frozen = folder / "config.json"
    if frozen.exists() and _cfg_identity(_read(frozen)) != _cfg_identity(cfg):
        raise RuntimeError("Frozen job configuration changed")
    if not frozen.exists():
        atomic_json(frozen, cfg)
    manifest = {"schema": 1, "identity": identity, "kind": kind,
                "code_identity": code, "config_identity": _cfg_identity(cfg),
                "config": str(frozen), "context": context, "commands": []}
    return folder, manifest


def _command(frozen, action, manifest, *extra):
    return ["bash", str(Path(_read(frozen)["code_root"]) / "f" / "entry.sh"), action,
            "--config", str(frozen), "--job-plan", str(manifest), *map(str, extra)]


def _publish_plan(cfg, folder, manifest):
    path = folder / "jobs.json"
    for command in manifest["commands"]:
        command["sha256"] = sha256(command["path"])
    manifest["manifest_digest"] = fingerprint(manifest)
    if path.exists() and _read(path) != manifest:
        raise RuntimeError(f"Immutable job manifest differs: {path}")
    if not path.exists():
        atomic_json(path, manifest)
    pointer = folder.parent / "current.txt"
    with locked(folder.parent / "publish.lock"):
        temporary = pointer.with_name("current.tmp")
        temporary.write_text(str(folder) + "\n")
        temporary.replace(pointer)
    print(f"Generated {manifest['kind']} command files: {folder}")
    print(f"Review: {path}")
    print("Submit from the scripts directory by omitting --plan-only and selecting --submit hpc or local.")
    return path


def generate_setup(cfg, metadata_source=None, download=True):
    from scheduler import write_cmd
    # Compute nodes have no external DNS; prefetch taxonomy on the submit host.
    meta = Path(cfg["metadata_dir"])
    if download and not (meta / "metadata.json").exists() and not Path(cfg["reference_plan"]).exists():
        from reference import _download
        meta.mkdir(parents=True, exist_ok=True)
        archive = meta / "taxdump.tar.gz"
        if not archive.is_file():
            _download("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", archive)
    source = str(Path(metadata_source).absolute()) if metadata_source else None
    folder, manifest = _new_plan(cfg, "setup", {"metadata_source": source, "download": download})
    path = folder / "00.setup.cmd"
    write_cmd(cfg, path, "setup", _command(folder / "config.json", "setup-worker", folder / "jobs.json"),
              f"metag_setup")
    manifest["commands"] = [{"role": "setup", "path": str(path)}]
    return _publish_plan(cfg, folder, manifest)


def _arrays(count, limit):
    return [(first, min(limit, count - first + 1)) for first in range(1, count + 1, limit)]


def _prepare_task_ids(first, total, stride=0):
    if not 1 <= first <= total or stride < 0:
        raise ValueError("Invalid reference preparation task range")
    return range(first, total + 1, stride) if stride else (first,)


def generate_library(cfg, prepared=False):
    from scheduler import write_cmd
    plan_path = Path(cfg["reference_plan"])
    if not plan_path.is_file():
        raise RuntimeError("Reference plan is not ready. Generate and submit 00.setup.cmd first.")
    reference = _read(plan_path)
    if (reference.get("status") != "complete" or
            fingerprint(reference.get("context")) != reference.get("identity") or
            reference.get("tasks") != reference.get("context", {}).get("tasks")):
        raise RuntimeError("Reference plan fingerprint or task list is invalid")
    tasks = reference.get("tasks", [])
    if not tasks or [t["task_id"] for t in tasks] != list(range(1, len(tasks) + 1)):
        raise RuntimeError("Reference tasks must be nonempty and consecutively numbered")
    folder, manifest = _new_plan(cfg, "library", {"reference_identity": reference["identity"],
                                                 "reference_tasks": len(tasks), "prepared": prepared})
    frozen, mf = folder / "config.json", folder / "jobs.json"
    workers = min(len(tasks), int(cfg["prepare_concurrency"]), int(cfg["max_array_size"]), 32)
    if workers < 1:
        raise ValueError("prepare_concurrency must be positive")
    if prepared:
        from reference import _load_plan, _prepared_receipt
        verified = _load_plan(cfg)
        for task in tasks:
            _prepared_receipt(cfg, verified, task)
        print(f"Verified {len(tasks)} completed preparation receipts; scheduling build directly")
    else:
        path = folder / "1.prepare.cmd"
        write_cmd(cfg, path, "prepare", _command(frozen, "prepare-worker", mf,
                  "--task-id", "__LSB_JOBINDEX__", "--stride", workers),
                  f"metag_prepare",
                  array=f"1-{workers}%{workers}")
        manifest["commands"].append({"role": "prepare", "path": str(path)})
    path = folder / "1.build.cmd"
    write_cmd(cfg, path, "build", _command(frozen, "build-worker", mf),
              f"metag_build")
    manifest["commands"].append({"role": "build", "path": str(path), "depends_on": [] if prepared else "prepare"})
    if cfg.get("build_kunpeng", True):
        path = folder / "1.kunpeng.cmd"
        write_cmd(cfg, path, "convert", _command(frozen, "convert-worker", mf),
                  f"metag_kunpeng_index")
        manifest["commands"].append({"role": "convert", "path": str(path), "depends_on": "build"})
    result = _publish_plan(cfg, folder, manifest)
    print(f"Reference preparation tasks: {len(tasks)}; worker jobs: {workers}")
    print("The preparation groups are not query databases: all groups feed one complete index.")
    return result


def _sample_states(sample):
    from query import _runtime_state
    return {key: _runtime_state(sample[key]) for key in ("r1", "r2") if sample.get(key)}


def _database_states(snapshot):
    from query import _database_runtime_states
    states = _database_runtime_states(snapshot)
    if snapshot.get("schema") == "metag-partitioned-index-v1":
        # The catalog and quick check establish library readiness. A queued
        # split must not fail just because library timestamps have changed.
        return {path: {"size": state["size"]} for path, state in states.items()}
    return states


def balanced_batches(samples, size):
    """Balance compressed input bytes, preserving every sample exactly once."""
    count = math.ceil(len(samples) / size)
    batches = [[] for _ in range(count)]
    weights = [0] * count
    weighted = [(sum(v["size"] for v in _sample_states(s).values()), i) for i, s in enumerate(samples)]
    for weight, i in sorted(weighted, reverse=True):
        dest = min((j for j in range(count) if len(batches[j]) < size), key=lambda j: (weights[j], j))
        batches[dest].append(i)
        weights[dest] += weight
    return batches


def generate_conversion(cfg):
    """Reuse an already completed Kraken database without scheduling preparation."""
    from scheduler import write_cmd
    from reference import database_snapshot
    snapshot = database_snapshot(cfg)
    folder, manifest = _new_plan(cfg, "convert", {"source_database_identity": snapshot["identity"]})
    path = folder / "1.kunpeng.cmd"
    write_cmd(cfg, path, "convert", _command(folder / "config.json", "convert-worker", folder / "jobs.json"),
              f"metag_kunpeng_index")
    manifest["commands"] = [{"role": "convert", "path": str(path)}]
    return _publish_plan(cfg, folder, manifest)


def generate_split(cfg, datasets=None, only=None, engine=None, ready_only=False, queue=None):
    from scheduler import write_cmd
    from query import discover_samples
    cfg = engine_resources(cfg, engine)
    if queue is not None:
        # Apply after engine selection so Kun-peng resource defaults cannot undo it.
        cfg["classify_queue"] = queue
        cfg["pilot_queue"] = queue
    # Every engine has a fresh namespace, including the previous Kraken baseline.
    cfg["output_dir"] = str(Path(cfg["output_dir"]) / cfg["classifier_engine"])
    snapshot = engine_snapshot(cfg)
    effective_workers(cfg, snapshot)  # Reject an impossible allocation before submission.
    samples = discover_samples(cfg, datasets=datasets, only=only, ready_only=ready_only)
    states = [_sample_states(s) for s in samples]
    batches = balanced_batches(samples, cfg["sample_batch_size"])
    # One prefix per dataset: an operational pilot, not a biological accuracy test.
    largest = {}
    for i, sample in enumerate(samples):
        amount = sum(x["size"] for x in states[i].values())
        old = largest.get(sample["dataset"])
        if old is None or amount > old[0]:
            largest[sample["dataset"]] = (amount, i)
    pilots = [largest[key][1] for key in sorted(largest)]
    context = {"database_identity": snapshot["identity"],
               "database_runtime_states": _database_states(snapshot), "samples": samples,
               "sample_states": states, "batches": batches, "pilot_indices": pilots}
    folder, manifest = _new_plan(cfg, "split", context)
    frozen, mf = folder / "config.json", folder / "jobs.json"
    with (folder / "samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "dataset", "sample", "r1", "r2"], delimiter="\t")
        writer.writeheader()
        writer.writerows(samples)
    path = folder / "2.pilot.cmd"
    write_cmd(cfg, path, "pilot", _command(frozen, "pilot-worker", mf),
              f"metag_pilot")
    manifest["commands"].append({"role": "pilot", "path": str(path)})
    for number in range(1, len(batches) + 1):
        path = folder / f"2.batch.{number:05d}.cmd"
        write_cmd(cfg, path, "classify", _command(frozen, "split-worker", mf, "--batch-id", number),
                  split_job_name(cfg, samples, number))
    for part, (first, count) in enumerate(_arrays(len(batches), cfg["max_array_size"]), 1):
        path = folder / f"2.split.array.{part:03d}.cmd"
        write_cmd(cfg, path, "classify", _command(frozen, "split-worker", mf,
                  "--batch-id", "__LSB_JOBINDEX__", "--offset", first - 1),
                  split_job_name(cfg, samples, part),
                  array=f"1-{count}%{min(cfg['classify_concurrency'], count)}")
        manifest["commands"].append({"role": "split", "path": str(path),
                                     "depends_on": ["pilot", "split"] if part > 1 else "pilot"})
    result = _publish_plan(cfg, folder, manifest)
    print(f"Samples={len(samples)}, batch jobs={len(batches)}, requested workers/job={cfg['classify_workers']}")
    print("Submission helper waits for the operational pilot to succeed, then releases batch jobs.")
    return result


def _compare_allocation(resources):
    """Canonical allocation rendered by write_cmd(), including rounded limits."""
    from scheduler import _settings, _memory
    queue, cores, memory, walltime = _settings(resources, "classify")
    return {"queue": queue, "cores": cores, "memory_gib": float(memory),
            "walltime": walltime, "memory_request": _memory(resources, memory, cores)}


def _validate_compare_allocation(context, case_index, resources, group=None):
    """Require a case's runtime budget to match its frozen submission group."""
    membership = [name for name, indices in context["groups"].items() if case_index in indices]
    if len(membership) != 1 or (group is not None and group != membership[0]):
        raise RuntimeError("Comparison case does not belong to exactly one submitted resource group")
    descriptor = context.get("group_resources", {}).get(membership[0])
    if not descriptor or descriptor.get("allocation") != _compare_allocation(resources):
        raise RuntimeError("Comparison worker allocation differs from its generated array request")
    return descriptor


def generate_compare(cfg, datasets=None, only=None, methods=None, repeats=None,
                     limit=None, profiles=None, truth=None, all_samples=False, placement=None):
    from compare import make_cases
    from query import discover_samples
    from scheduler import write_cmd
    cfg = dict(cfg)
    if placement is not None:
        cfg["compare_placement"] = placement
    samples = discover_samples(cfg, datasets, only)
    if not only and not all_samples:
        largest = {}
        for sample in samples:
            size = sum(x["size"] for x in _sample_states(sample).values())
            if sample["dataset"] not in largest or size > largest[sample["dataset"]][0]:
                largest[sample["dataset"]] = (size, sample)
        samples = [largest[k][1] for k in sorted(largest)]
    chosen = methods.split(",") if isinstance(methods, str) else methods
    if cfg.get("database_layout") == "partitioned" and chosen is None and not profiles:
        chosen = ["kraken"]
    profile_values = read_config_file(profiles) if profiles else None
    if isinstance(profile_values, dict):
        profile_values = profile_values["compare_profiles"]
    cases = make_cases(samples, methods=chosen, repeats=cfg["compare_repeats"] if repeats is None else repeats,
                       limit=cfg["compare_fragments"] if limit is None else limit, profiles=profile_values)
    methods_used = sorted({c["classifier_engine"] for c in cases})
    if cfg.get("database_layout") == "partitioned" and any(name != "kraken" for name in methods_used):
        raise ValueError("The partitioned catalog supports the CPU union engine only. Kun-peng needs its own completed single-index configuration; no final-taxid merge is substituted.")
    snapshots, database_states = {}, {}
    sources = set()
    for method in methods_used:
        resources = engine_resources(cfg, method, benchmark=True)
        snapshot = engine_snapshot(resources)
        effective_workers(resources, snapshot)
        snapshots[method] = snapshot["identity"]
        database_states[method] = _database_states(snapshot)
        sources.add(snapshot.get("source_database_identity", snapshot["identity"]))
    if len(sources) != 1:
        raise RuntimeError("Comparison engines must use the same verified source reference index")
    groups, group_resources = {}, {}
    for i, case in enumerate(cases):
        resources = engine_resources(cfg, case["classifier_engine"], benchmark=True)
        allocation = _compare_allocation(resources)
        matched = cfg.get("compare_placement", "matched") == "matched"
        cohort = "memory" if not matched and case["classifier_engine"] == "kunpeng_memory" else "full"
        descriptor = {"cohort": cohort, "allocation": allocation}
        key = "matched" if matched else next((name for name, value in group_resources.items()
                                               if value == descriptor), f"{cohort}{len(group_resources)+1}")
        if key in group_resources and group_resources[key] != descriptor:
            raise RuntimeError("Comparison resource group contains different allocations")
        groups.setdefault(key, []).append(i)
        group_resources[key] = descriptor
    truth_file = str(Path(truth).resolve(strict=True)) if truth else None
    context = {"cases": cases, "groups": groups, "group_resources": group_resources,
               "database_identities": snapshots,
               "database_runtime_states": database_states,
               "source_database_identity": next(iter(sources)),
               "sample_states": {s["sample_id"]: _sample_states(s) for s in samples},
               "truth": ({"path": truth_file, "state": file_state(truth_file), "sha256": sha256(truth_file)}
                         if truth_file else None)}
    folder, manifest = _new_plan(cfg, "compare", context)
    frozen, mf = folder / "config.json", folder / "jobs.json"
    for index, case in enumerate(cases, 1):
        resources = engine_resources(cfg, case["classifier_engine"], benchmark=True)
        path = folder / f"3.case.{index:05d}.cmd"
        write_cmd(resources, path, "classify", _command(frozen, "compare-worker", mf, "--case-id", index),
                  f"metag_compare_{index}")
    previous_role = {}
    for group_number, (group, indices) in enumerate(groups.items(), 1):
        resources = engine_resources(cfg, cases[indices[0]]["classifier_engine"], benchmark=True)
        cohort = group_resources[group]["cohort"]
        concurrency = cfg[f"compare_concurrency_{cohort}"]
        role = f"compare_{group}"
        for part, (first, count) in enumerate(_arrays(len(indices), cfg["max_array_size"]), 1):
            path = folder / f"3.compare.{group}.{part:03d}.cmd"
            write_cmd(resources, path, "classify", _command(frozen, "compare-worker", mf,
                      "--case-id", "__LSB_JOBINDEX__", "--group", group, "--offset", first - 1),
                      f"metag_compare_{cohort}_{group_number}_{part}",
                      array=f"1-{count}%{min(concurrency, count)}")
            # Distinct allocations need distinct arrays. Serializing those arrays
            # within each cohort prevents splitting groups from multiplying the
            # configured full/memory concurrency limit.
            parent = role if part > 1 else previous_role.get(cohort)
            manifest["commands"].append({"role": role, "path": str(path),
                **({"depends_on": parent, "dependency_condition": "ended"} if parent else {})})
        previous_role[cohort] = role
    path = folder / "3.summary.cmd"
    write_cmd(cfg, path, "compare_summary", _command(frozen, "compare-summary-worker", mf),
              f"metag_compare_summary")
    manifest["commands"].append({"role": "compare_summary", "path": str(path),
        "depends_on": [f"compare_{key}" for key in groups], "dependency_condition": "ended"})
    result = _publish_plan(cfg, folder, manifest)
    print(f"Comparison: {len(samples)} samples, {len(cases)} cases; placement={cfg.get('compare_placement', 'matched')}")
    print(f"Reports after jobs finish: {folder / 'report'}")
    return result


def _validate_job_plan(cfg, path, kind=None):
    if not path:
        raise RuntimeError("Compute command requires its generated --job-plan")
    plan = _read(path)
    if plan.get("schema") != 1 or (kind and plan.get("kind") != kind):
        raise RuntimeError("Invalid job plan kind/schema")
    expected = fingerprint({"schema": 1, "kind": plan["kind"], "config": plan["config_identity"],
                            "code": plan["code_identity"], "context": plan["context"]})
    if expected != plan.get("identity") or plan.get("manifest_digest") != fingerprint(
            {key: value for key, value in plan.items() if key != "manifest_digest"}):
        raise RuntimeError("Generated job manifest changed or has an invalid fingerprint")
    if plan.get("code_identity") != code_identity():
        raise RuntimeError("Code changed after .cmd generation; regenerate commands before running this job")
    if plan.get("config_identity") != _cfg_identity(cfg):
        raise RuntimeError("Frozen job configuration changed")
    seal = _integrity_seal(cfg, plan["context"].get("integrity_seal"),
                           _required_sealed_indexes(plan["kind"], cfg, plan["context"]))
    if seal != plan["context"].get("integrity_seal"):
        raise RuntimeError("Frozen job plan has no matching integrity seal")
    for command in plan["commands"]:
        if sha256(command["path"]) != command.get("sha256"):
            raise RuntimeError(f"Generated .cmd changed: {command['path']}; update config.json and regenerate")
    return plan


def _active_jobs(ids):
    if not ids:
        return []
    if not shutil.which("bjobs"):
        raise RuntimeError("Cannot verify previous jobs without bjobs; refusing automatic resubmission")
    result = subprocess.run(["bjobs", "-a", "-noheader", "-o", "jobid stat", *map(str, ids)],
                            text=True, capture_output=True, timeout=30)
    if result.returncode:
        raise RuntimeError("Cannot establish that previous jobs finished; inspect bjobs/bhist before retrying. "
                           + result.stderr.strip())
    active, observed = [], set()
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) < 2 or not fields[0].split("[", 1)[0].isdigit():
            raise RuntimeError("Unexpected bjobs status format; refusing automatic resubmission")
        observed.add(int(fields[0].split("[", 1)[0]))
        if fields[1] not in ("DONE", "EXIT"):
            active.append(line.strip())
    if set(ids) - observed:
        raise RuntimeError("Some prior JOBIDs are absent from bjobs history; inspect bhist/bacct before recovery")
    return active


def submit_plan(path, retry=False):
    from scheduler import submit_cmd, SubmissionError
    path = Path(path).absolute()
    initial = _read(path)
    cfg = load_config(initial["config"])
    manifest = _validate_job_plan(cfg, path)
    target = path.parent / "submission.json"
    with locked(path.parent / "submission.lock"):
        history = _read(target) if target.exists() else {"schema": 1, "attempts": []}
        prior_jobs = [job for attempt in history["attempts"] for job in attempt.get("jobs", [])]
        previous_ids = sorted({job["job_id"] for job in prior_jobs if isinstance(job.get("job_id"), int)})
        if history["attempts"] and not retry:
            raise RuntimeError(f"This plan was already submitted: {previous_ids}. "
                               "Inspect bjobs/status; use --retry only after failed/finished jobs have exited.")
        if any(job.get("status") in ("intent", "unknown") and not job.get("job_id") for job in prior_jobs):
            raise RuntimeError("A prior bsub acceptance is unknown. Inspect the stored submission receipt and "
                               "bjobs/bhist; no automatic retry is safe. The individual .cmd files remain available for controlled recovery.")
        if previous_ids:
            active = _active_jobs(previous_ids)
            if active:
                raise RuntimeError("Previous jobs are still active; refusing duplicate submission: " + "; ".join(active))
        attempt = {"started_at": utcnow(), "jobs": [], "status": "submitting"}
        history["attempts"].append(attempt)
        atomic_json(target, history)
        by_role = {}
        try:
            for command in manifest["commands"]:
                role = command["role"]
                dependencies = command.get("depends_on", [])
                if isinstance(dependencies, str):
                    dependencies = [dependencies]
                if any(not by_role.get(parent) for parent in dependencies):
                    raise RuntimeError(f"Missing successful dependency submission for {role}")
                parents = [job_id for parent in dependencies for job_id in by_role[parent]]
                condition = command.get("dependency_condition", "done")
                if condition not in ("done", "ended"):
                    raise RuntimeError("Invalid dependency condition")
                dependency = " && ".join(f"{condition}({job_id})" for job_id in parents) or None
                position = len(attempt["jobs"])
                attempt["jobs"].append({"status": "intent", "role": role,
                                         "command_path": command["path"], "dependency": dependency})
                atomic_json(target, history)
                try:
                    result = submit_cmd(cfg, command["path"], dependency=dependency, expected_sha256=command["sha256"])
                except SubmissionError as exc:
                    attempt["jobs"][position] = {**exc.result, "role": role}
                    atomic_json(target, history)
                    raise
                result["role"] = role
                attempt["jobs"][position] = result
                by_role.setdefault(role, []).append(result["job_id"])
                atomic_json(target, history)
                print(f"Submitted {role}: {result['job_id']}" + (f"; depends on {dependency}" if dependency else ""))
            attempt["status"] = "submitted"
        except BaseException as exc:
            attempt.update(status="submission_failed", error=str(exc))
            raise
        finally:
            atomic_json(target, history)
    return history


def run_worker(args, cfg):
    require_compute(args.allow_local)
    if args.command == "convert-worker":
        kind = _read(args.job_plan).get("kind")
        if kind not in ("library", "convert"):
            raise RuntimeError("Conversion requires a library or conversion plan")
    elif args.command in ("compare-worker", "compare-summary-worker"):
        kind = "compare"
    else:
        kind = "setup" if args.command == "setup-worker" else "library" if args.command in ("prepare-worker", "build-worker") else "split"
    plan = _validate_job_plan(cfg, args.job_plan, kind)
    if args.command == "convert-worker":
        from kunpeng_index import install, convert_database
        from reference import database_snapshot
        expected = plan["context"].get("source_database_identity")
        snapshot = database_snapshot(cfg)
        if expected and snapshot["identity"] != expected:
            raise RuntimeError("Source database changed after conversion job generation")
        if kind == "library" and _read(cfg["reference_plan"])["identity"] != plan["context"]["reference_identity"]:
            raise RuntimeError("Reference plan changed after library generation")
        install(cfg)
        convert_database(cfg)
        return
    if kind == "compare":
        from compare import run_case, summarize_cases
        context = plan["context"]
        cases = context["cases"]
        folder = Path(args.job_plan).parent
        if args.command == "compare-summary-worker":
            truth = context.get("truth")
            if truth and (file_state(truth["path"]) != truth["state"] or sha256(truth["path"]) != truth["sha256"]):
                raise RuntimeError("Truth labels changed after comparison generation")
            # Use a taxonomy that was actually frozen for these cases. A
            # Kun-Peng-only comparison need not reopen an unrecorded Kraken
            # snapshot just to calculate truth metrics.
            selected = "kraken" if "kraken" in context["database_identities"] else sorted(context["database_identities"])[0]
            source_cfg = dict(cfg, classifier_engine=selected)
            snap = engine_snapshot(source_cfg)
            if (snap["identity"] != context["database_identities"][selected] or
                    snap.get("source_database_identity", snap["identity"]) != context["source_database_identity"]):
                raise RuntimeError("Reference changed before comparison summary")
            if _database_states(snap) != context["database_runtime_states"][selected]:
                raise RuntimeError("Index/taxonomy runtime state changed before comparison summary")
            report = summarize_cases([folder / "cases" / c["case_id"] / "case.json" for c in cases],
                                     folder / "report", truth_path=truth["path"] if truth else None,
                                     taxonomy_dir=snap["taxonomy_dir"])
            print(json.dumps(report, indent=2))
            if report.get("status") != "complete":
                raise RuntimeError("Comparison is incomplete; see report for failed/missing/mismatched cases")
            from storage_layout import publish_comparison
            print(f"Published comparison: {publish_comparison(cfg, folder / 'report')}")
            return
        position = args.case_id + args.offset - 1
        indices = context["groups"].get(args.group) if args.group else list(range(len(cases)))
        if indices is None or not 0 <= position < len(indices):
            raise RuntimeError("Comparison case index outside generated plan")
        case = cases[indices[position]]
        sample = case["sample"]
        if _sample_states(sample) != context["sample_states"][sample["sample_id"]]:
            raise RuntimeError("Comparison FASTQ changed after job generation")
        resources = engine_resources(cfg, case["classifier_engine"], benchmark=True)
        _validate_compare_allocation(context, indices[position], resources, args.group)
        snapshot = engine_snapshot(resources)
        if snapshot["identity"] != context["database_identities"][case["classifier_engine"]]:
            raise RuntimeError("Comparison database changed after job generation")
        if _database_states(snapshot) != context["database_runtime_states"][case["classifier_engine"]]:
            raise RuntimeError("Comparison index/taxonomy runtime state changed after job generation")
        resources = effective_workers(resources, snapshot)
        run_case(resources, case, folder / "cases" / case["case_id"], snapshot=snapshot)
        return
    if args.command == "setup-worker":
        from software_install import install
        from reference import metadata, plan_references
        install(cfg)
        from kunpeng_index import install as install_kunpeng
        install_kunpeng(cfg)
        metadata(cfg, source=plan["context"]["metadata_source"], download=plan["context"]["download"])
        plan_references(cfg)
        generated = generate_library(cfg)
        atomic_json(Path(args.job_plan).parent / "setup.complete.json",
                    {"status": "complete", "finished_at": utcnow(), "library_plan": str(generated)})
        return
    if kind == "library":
        from reference import prepare_task, build_database
        if _read(cfg["reference_plan"])["identity"] != plan["context"]["reference_identity"]:
            raise RuntimeError("Reference plan changed after submission")
        if args.command == "prepare-worker":
            total = len(_read(cfg["reference_plan"])["tasks"])
            for task_id in _prepare_task_ids(args.task_id + args.offset, total, args.stride):
                print(f"Preparing reference task {task_id}/{total}", flush=True)
                prepare_task(cfg, task_id)
        else:
            build_database(cfg)
        return
    from query import run_sample, run_batch, BatchFailure
    snapshot = engine_snapshot(cfg)
    if snapshot["identity"] != plan["context"]["database_identity"]:
        raise RuntimeError("Database changed after job generation")
    if _database_states(snapshot) != plan["context"]["database_runtime_states"]:
        raise RuntimeError("Index/taxonomy runtime state changed after job generation")
    cfg = effective_workers(cfg, snapshot)
    context = plan["context"]
    if args.command == "pilot-worker":
        indices = context["pilot_indices"]
    else:
        batch_id = args.batch_id + args.offset
        if not 1 <= batch_id <= len(context["batches"]):
            raise RuntimeError("Batch ID is outside the generated plan")
        indices = context["batches"][batch_id - 1]
    samples = [context["samples"][i] for i in indices]
    for index in indices:
        if _sample_states(context["samples"][index]) != context["sample_states"][index]:
            raise RuntimeError(f"FASTQ changed after job generation: {context['samples'][index]['sample_id']}")
    if args.command == "pilot-worker":
        outputs = []
        for sample in samples:
            outputs.append(run_sample(cfg, sample, snapshot=snapshot, limit=cfg["pilot_fragments"]))
        atomic_json(Path(args.job_plan).parent / "pilot.complete.json",
                    {"status": "complete", "finished_at": utcnow(), "prefix_fragments": cfg["pilot_fragments"],
                     "note": "Operational prefix test, not biological sensitivity validation or whole-sample throughput",
                     "outputs": outputs})
    else:
        report = Path(args.job_plan).parent / f"batch.{batch_id:05d}.json"
        runtime = {"lsf_job_id": os.environ.get("LSB_JOBID"), "host": os.uname().nodename,
                   "requested_workers": cfg.get("requested_classify_workers"),
                   "effective_workers": cfg["classify_workers"], "accounting": cfg["memory_accounting"]}
        try:
            result = run_batch(cfg, samples, snapshot=snapshot)
        except BatchFailure as exc:
            atomic_json(report, {"status": "failed", "finished_at": utcnow(),
                                 "runtime": runtime, "result": exc.summary})
            raise
        atomic_json(report, {"status": "complete", "finished_at": utcnow(),
                             "runtime": runtime, "result": result})


def status(cfg):
    print(f"Root: {cfg['root']}")
    for label, path in (("reference plan", cfg["reference_plan"]),
                        ("Kraken database", str(Path(cfg["db_dir"]) / "database.json")),
                        ("Kun-peng database", str(Path(cfg["kunpeng_db_dir"]) / "database.json"))):
        print(f"{label}: {'present' if Path(path).is_file() else 'not ready'} — {path}")
    outputs = Path(cfg["output_dir"])
    completed = []
    for engine in ENGINES:
        markers = [p for p in (outputs / engine).glob("*/*/complete.json")
                   if ".partial." not in p.parent.name]
        completed.extend(markers)
        print(f"Published {engine} sample markers: {len(markers)}")
    print(f"Published sample markers: {len(completed)} (marker count only, not a content audit)")
    jobs = Path(cfg["work_dir"]) / "jobs"
    for path in sorted(jobs.glob("*/submission.json")) if jobs.exists() else []:
        value = _read(path)
        for attempt in value.get("attempts", []):
            print(path.parent.name, attempt["status"], "JOBIDs=", [j.get("job_id") for j in attempt["jobs"]])
    if shutil.which("bjobs"):
        subprocess.run(["bjobs"], check=False)


def boolean_value(value):
    if value.upper() == "TRUE":
        return True
    if value.upper() == "FALSE":
        return False
    raise argparse.ArgumentTypeError("Expected TRUE or FALSE")


def queue_name(value):
    import re
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.+-]*", value):
        raise argparse.ArgumentTypeError("invalid LSF queue name")
    return value


def parser():
    examples = """Optional comparison of already completed classifier databases:
  ./1.build-lib.sh check --benchmark --submit hpc --plan-only
--submit hpc|local|auto selects the execution machine; --plan-only writes commands.
Queue overrides do not change CPU or memory requests.
"""
    result = argparse.ArgumentParser(description=__doc__, epilog=examples,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    result.add_argument("--config", help="Runtime JSON configuration (default: config.json)")
    result.add_argument("--root", help="Override root for relative configured paths")
    sub = result.add_subparsers(dest="command")
    def command(name, help_text):
        obj = sub.add_parser(name, help=help_text, formatter_class=argparse.RawDescriptionHelpFormatter)
        obj.add_argument("--config", default=argparse.SUPPRESS)
        obj.add_argument("--root", default=argparse.SUPPRESS)
        return obj
    command("doctor", "Read-only LSF, node, queue and filesystem inventory")
    obj = command("setup", "Generate the ordinary-node setup .cmd; no compute on login node")
    obj.add_argument("--metadata-source")
    obj.add_argument("--offline", action="store_true", help="Import local metadata only; no NCBI downloads")
    obj.add_argument("--submit", choices=("hpc", "local", "auto"), default=None, help="Execution machine [auto]")
    obj.add_argument("--plan-only", action="store_true", help="Write commands without running or submitting")
    obj = command("library", "Generate preparation .cmd files and one uncapped build .cmd")
    obj.add_argument("--prepared", action="store_true", help="Validate completed preparation receipts and skip prepare jobs")
    obj.add_argument("--submit", choices=("hpc", "local", "auto"), default=None, help="Execution machine [auto]")
    obj.add_argument("--plan-only", action="store_true", help="Write commands without running or submitting")
    obj.add_argument("--kunpeng", action="store_true", help="Only convert an already completed Kraken index")
    obj.add_argument("--queue", help="Reference preparation queue only (not build or index conversion)")
    obj.add_argument("--build-queue", type=queue_name, help="Build queue (configured default smp, 4096 GiB)")
    obj.add_argument("--convert-queue", type=queue_name, help="Kunpeng index conversion queue (configured default short, 40 cores; memory from 0f/0os.sh: local 32 / HPC 160 GiB)")
    obj.epilog = examples
    obj = command("split", "Classify reads and split by species using pilot and batch jobs")
    obj.add_argument("--queue", type=queue_name, help="Queue for both pilot and all split batches; CPU/memory requests remain unchanged")
    obj.add_argument("--dat-list", default=None)
    obj.add_argument("--sample", action="append", help="dataset/sample or an unambiguous sample name")
    obj.add_argument("--submit", choices=("hpc", "local", "auto"), default=None, help="Execution machine [auto]")
    obj.add_argument("--plan-only", action="store_true", help="Write commands without running or submitting")
    obj.add_argument("--engine", choices=ENGINES, help="Select an actual native engine/mode")
    obj.add_argument("--replace", type=boolean_value, default=True)
    obj = command("compare", "Generate matched-input method comparison jobs and a terminal-status summary")
    obj.add_argument("--dat-list")
    obj.add_argument("--sample", action="append")
    obj.add_argument("--all-samples", action="store_true", help="Include every discovered sample instead of one per dataset")
    obj.add_argument("--methods", help="Comma-separated engines (default all three)")
    obj.add_argument("--profiles", help="JSON list of named engines and predeclared threshold overrides")
    obj.add_argument("--repeats", type=int)
    obj.add_argument("--fragments", type=int, help="Matched prefix fragments; 0 means complete samples")
    obj.add_argument("--truth", help="Optional complete sample_id/ordinal/species_taxid TSV labels")
    obj.add_argument("--placement", choices=("matched", "production"))
    obj.add_argument("--replace", type=boolean_value, default=True)
    obj.add_argument("--submit", choices=("hpc", "local", "auto"), default=None, help="Execution machine [auto]")
    obj.add_argument("--plan-only", action="store_true", help="Write commands without running or submitting")
    obj = command("submit", "Submit a generated plan with explicit success dependencies")
    obj.add_argument("--manifest", required=True)
    obj.add_argument("--retry", action="store_true")
    command("status", "Show checkpoint and submitted LSF job status")
    obj = command("audit", "Read-only database file-state and uncapped-option verification")
    obj.add_argument("--engine", choices=ENGINES, default="kraken")
    for name in ("setup-worker", "prepare-worker", "build-worker", "convert-worker", "pilot-worker", "split-worker", "compare-worker", "compare-summary-worker"):
        obj = command(name, "Compute entry used by generated .cmd files")
        obj.add_argument("--job-plan", required=True)
        obj.add_argument("--allow-local", action="store_true", help="Deliberately allow a small local test")
        if name == "prepare-worker":
            obj.add_argument("--task-id", required=True, type=int)
            obj.add_argument("--stride", type=int, default=0)
            obj.add_argument("--offset", type=int, default=0)
        if name == "split-worker":
            obj.add_argument("--batch-id", required=True, type=int)
            obj.add_argument("--offset", type=int, default=0)
        if name == "compare-worker":
            obj.add_argument("--case-id", required=True, type=int)
            obj.add_argument("--group")
            obj.add_argument("--offset", type=int, default=0)
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    if args.command is None:
        args.command, args.metadata_source, args.offline, args.submit, args.plan_only = "setup", None, False, None, False
    if hasattr(args, "submit"):
        from common import resolve_mode
        os.environ["submit"] = os.environ["SUBMIT"] = resolve_mode(args.submit)
    cfg = load_config(args.config, args.root)
    if args.command == "doctor":
        from scheduler import doctor
        report = doctor(cfg)
        atomic_json(Path(cfg["work_dir"]) / "doctor.json", report)
        print(f"Saved inventory: {Path(cfg['work_dir']) / 'doctor.json'}")
    elif args.command == "setup":
        path = generate_setup(cfg, args.metadata_source, not args.offline)
        if not args.plan_only:
            submit_plan(path)
    elif args.command == "library":
        if args.build_queue:
            cfg["build_queue"] = args.build_queue
        if args.convert_queue:
            cfg["convert_queue"] = args.convert_queue
        if args.queue:
            cfg["prepare_queue"] = args.queue
        if cfg["prepare_queue"] in ("short", "spec"):
            cfg["prepare_cores"] = max(40, int(cfg["prepare_cores"]))
        path = generate_conversion(cfg) if args.kunpeng else generate_library(cfg, args.prepared)
        if not args.plan_only:
            submit_plan(path)
    elif args.command == "split":
        cfg["replace"] = args.replace
        path = generate_split(cfg, args.dat_list, args.sample, args.engine, False, args.queue)
        if not args.plan_only:
            submit_plan(path)
    elif args.command == "compare":
        cfg["replace"] = args.replace
        path = generate_compare(cfg, args.dat_list, args.sample, args.methods, args.repeats,
                                args.fragments, args.profiles, args.truth, args.all_samples, args.placement)
        if not args.plan_only:
            submit_plan(path)
    elif args.command == "submit":
        submit_plan(args.manifest, retry=args.retry)
    elif args.command == "status":
        status(cfg)
    elif args.command == "audit":
        value = engine_snapshot(dict(cfg, classifier_engine=args.engine))
        print("Quick audit: receipt structure, recorded paths/sizes and native options; no content or timestamp audit. Use f/audit_integrity.py for a separate verified baseline.")
        print(json.dumps(value, indent=2))
    else:
        run_worker(args, cfg)
        # With a sealed reference baseline, also detect reference changes made
        # during a long compute stage. Downstream done() dependencies stay blocked.
        if cfg.get("integrity_seal"):
            _integrity_seal(cfg, _read(args.job_plan)["context"].get("integrity_seal"))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, OSError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
