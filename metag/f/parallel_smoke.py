#!/usr/bin/env python3
"""Run a tiny real five-category pipeline on one compute node, without bsub/network.

Examples (run after loading the compiler/Python environment used by metag):

  python3 f/parallel_smoke.py --destination /work/PROJECT/metag_smoke_NEW \
      --tools-dir /work/PROJECT/tools/metag_parallel

  python3 f/parallel_smoke.py --destination /tmp/metag_smoke_NEW --allow-local

The destination must be new or empty. Existing tools may be reused through their
verified installation receipts; omitted tools are compiled from bundled source.
The script generates immutable setup/library command manifests, invokes their
real worker CLI actions on the current node, verifies a native pause/resume,
publishes five indexes, and classifies/splits six original paired FASTQ records.
No scheduler commands, external downloads, or production references are used.
The final stdout is one JSON result; detailed per-step output is saved in logs/.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time


CATEGORIES = ("microbes", "viruses", "plants", "invertebrates", "vertebrates")


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def ensure(value, message):
    if not value:
        raise RuntimeError(message)


def make_fixture(root):
    """Two whole FASTA files per category and six PE fragments with exact truth."""
    source = root / "source"
    source.mkdir()
    nodes = [(1, 1, "no rank"), (2, 1, "superkingdom")]
    mapping = ["path\ttaxid\n"]
    sequences = {}
    rng = random.Random(16384)
    for taxid, category in enumerate(CATEGORIES, 10):
        nodes.append((taxid, 2, "species"))
        for part in range(2):
            path = root / "reference_raw" / category / f"reference{taxid}_{part}.fna"
            path.parent.mkdir(parents=True, exist_ok=True)
            seq = "".join(rng.choices("ACGT", k=25000))
            path.write_text(f">synthetic{taxid}_{part}\n" +
                            "\n".join(seq[i:i+80] for i in range(0, len(seq), 80)) + "\n")
            mapping.append(f"{category}/{path.name}\t{taxid}\n")
            if part == 0:
                sequences[taxid] = seq.encode()
    (root / "mapping.tsv").write_text("".join(mapping))
    (source / "nodes.dmp").write_text("".join(
        f"{t}\t|\t{p}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
        for t, p, rank in nodes))
    (source / "names.dmp").write_text("".join(
        f"{t}\t|\tSynthetic taxon {t}\t|\t\t|\tscientific name\t|\n" for t, _, _ in nodes))
    (source / "merged.dmp").write_text("")
    truth = [10, 11, 12, 13, 14, 0]
    originals = [[], []]
    reads = root / "reads_raw" / "smoke"
    reads.mkdir(parents=True)
    for ordinal, taxid in enumerate(truth):
        for mate in range(2):
            offset = 1000 + ordinal * 101 + mate * 317
            sequence = sequences[taxid][offset:offset+150] if taxid else b"N" * 150
            # Duplicate identifiers are intentional: assignments must use the
            # validated ordinal stream, not a dictionary keyed by read name.
            record = (b"@duplicate/" + str(mate+1).encode() + b" original comment\r\n" +
                      sequence + b"\r\n+original plus comment\r\n" + b"I" * 150 +
                      (b"" if ordinal == len(truth)-1 else b"\r\n"))
            originals[mate].append(record)
    for mate in range(2):
        with gzip.open(reads / f"tiny_R{mate+1}.fastq.gz", "wb", compresslevel=1) as handle:
            handle.write(b"".join(originals[mate]))
    return truth, originals


def settings(root, tools):
    # The build estimate deliberately reserves 8 GiB in production. A 16 GiB
    # declared budget passes that guard; the synthetic hash itself is tiny.
    category = {"queue": "38", "cores": 2, "threads": 2, "memory_gib": 16, "hours": 120}
    return {
        "schema": 1, "root": str(root), "reference_raw": str(root / "reference_raw"),
        "reads_raw": str(root / "reads_raw"), "metadata_dir": str(root / "metadata"),
        "library_dir": str(root / "library"), "work_dir": str(root / "work"),
        "output_dir": str(root / "output"), "tools_dir": str(tools),
        "reference_taxid_map": str(root / "mapping.tsv"),
        "database_layout": "partitioned", "classifier_engine": "kraken", "build_kunpeng": False,
        "legacy_library_dir": None, "legacy_metadata_dir": None, "parallel_reference_dir": None,
        "integrity_seal": None, "category_resources": {name: dict(category) for name in CATEGORIES},
        "masking": True, "prepare_target_gib": 0.0001, "prepare_max_files": 1,
        "prepare_threads": 2, "prepare_cores": 2, "prepare_concurrency": 3,
        "prepare_memory_gib": 1, "prepare_queue": "38", "prepare_hours": 1,
        "setup_cores": 2, "setup_memory_gib": 4, "setup_queue": "38", "setup_hours": 1,
        "build_threads": 2, "build_cores": 2, "build_memory_gib": 16, "build_queue": "38",
        "build_hours": 120, "build_soft_hours": 110, "checkpoint_headroom_hours": 10,
        "checkpoint_hours": 4, "native_stop_after_files": 0,
        "classify_threads": 2, "classify_cores": 4, "classify_workers": 1,
        "classify_memory_gib": 1, "classify_queue": "38", "classify_hours": 1,
        "classify_concurrency": 1, "worker_memory_gib": 0.25, "system_reserve_gib": 0,
        "memory_accounting": "conservative", "min_free_gib": 0, "max_array_size": 3,
        "sample_batch_size": 1, "pilot_fragments": 3, "output_format": "directory",
        "confidence": 0.1, "minimum_hit_groups": 2, "minimum_base_quality": 0,
    }


class Runner:
    def __init__(self, root, code, allow_local):
        self.root, self.code, self.allow_local = root, code, allow_local
        self.steps = []
        self.logs = root / "logs"
        self.logs.mkdir()
        self.environment = dict(os.environ, OMP_NUM_THREADS="2", K2_DB_READ_THREADS="2",
                                OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def cli(self, label, args, expected=0):
        command = [sys.executable, str(self.code / "f/parallel_cli.py"), *map(str, args)]
        log = self.logs / f"{len(self.steps)+1:02d}.{label}.log"
        print(f"SMOKE {label}", file=sys.stderr, flush=True)
        started = time.monotonic()
        with log.open("wb") as handle:
            process = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=self.environment)
        self.steps.append({"step": label, "command": command, "exit_code": process.returncode,
                           "seconds": time.monotonic()-started, "log": str(log)})
        if process.returncode != expected:
            detail = log.read_text(errors="replace")[-12000:]
            raise RuntimeError(f"{label}: expected exit {expected}, got {process.returncode}; log={log}\n{detail}")

    def manifest(self, kind):
        pointer = self.root / "work/jobs" / kind.removeprefix("parallel_") / "current.txt"
        path = Path(pointer.read_text().strip()) / "jobs.json"
        ensure(path.is_file(), f"Missing generated {kind} job manifest")
        return path

    def worker(self, manifest, action, *extra, expected=0, label=None):
        plan = read_json(manifest)
        args = [action, "--config", plan["config"], "--job-plan", manifest, *extra]
        if self.allow_local:
            args.append("--allow-local")
        self.cli(label or action, args, expected)


def verify_split(folder, truth, originals):
    payload = folder
    with gzip.open(payload / "qc/read-assignments.tsv.gz", "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assignments = [int(row["species_taxid"]) for row in rows]
    ensure(assignments == truth, f"Assignments differ from synthetic truth: {assignments} != {truth}")
    for bin_name in sorted({row["bin"] for row in rows}):
        base = payload / bin_name
        for mate in range(2):
            with gzip.open(payload / f"{bin_name}.R{mate+1}.fastq.gz", "rb") as handle:
                actual = handle.read()
            expected = b"".join(originals[mate][i] for i, row in enumerate(rows)
                                if row["bin"] == bin_name)
            ensure(actual == expected, f"Original FASTQ bytes changed for bin={bin_name}, mate={mate+1}")
    ensure(rows[-1]["bin"] == "unclassified", "True no-call was not separated")
    ensure(all(row["bin"] == "bacteria/species_other" for row in rows[:-1]),
           "Species below the storage threshold did not share their group file")
    return {"paired_fragments": len(truth), "species_taxids": truth, "exact_original_fastq_bytes": True,
            "duplicate_read_ids": True, "CRLF_plus_comments_and_missing_final_newline": True}


def smoke(args):
    code = args.code_root.expanduser().resolve(strict=True)
    sys.path.insert(0, str(code / "f"))
    from common import require_compute, sha256, utcnow
    require_compute(args.allow_local)
    root = args.destination.expanduser().absolute()
    ensure(not root.is_symlink() and (not root.exists() or (root.is_dir() and not any(root.iterdir()))),
           "--destination must be a new or empty directory; no existing files will be removed")
    root.mkdir(parents=True, exist_ok=True)
    root = root.resolve(strict=True)
    tools = args.tools_dir.expanduser().resolve(strict=True) if args.tools_dir else root / "tools"
    ensure(not args.tools_dir or tools.is_dir(), "--tools-dir must be an existing directory")
    started = time.monotonic()
    report = {"schema": 1, "status": "running", "destination": str(root), "code_root": str(code),
              "tools_dir": str(tools), "started_at": utcnow(), "scheduler_submitted": False,
              "network_requested": False}
    runner = Runner(root, code, args.allow_local)
    try:
        truth, originals = make_fixture(root)
        config = settings(root, tools)
        # The explicitly chosen fixture destination is disposable test data.
        config["allow_scratch"] = True
        config["storage_layout"] = args.storage_layout
        if args.storage_layout == "separated":
            config["output_dir"] = str(root / "output/kraken")
        base_config, pause_config = root / "config.json", root / "config.pause.json"
        write_json(base_config, config)
        write_json(pause_config, {**config, "native_stop_after_files": 1})

        runner.cli("generate-setup", ["setup", "--config", base_config, "--metadata-source", root / "source",
                                      "--offline", "--submit", "local", "--plan-only"])
        setup = runner.manifest("parallel_setup")
        runner.worker(setup, "parallel-setup-worker", label="setup-tools-and-plan")

        runner.cli("generate-pause-library", ["library", "--config", pause_config, "--submit", "local", "--plan-only"])
        paused_plan = runner.manifest("parallel_library")
        paused = read_json(paused_plan)
        ensure(paused["context"]["reference_tasks"] == 10, "Expected ten fresh preparation tasks")
        for first in range(1, 4):
            runner.worker(paused_plan, "parallel-prepare-worker", "--task-id", first, "--stride", 3,
                          label=f"prepare-stride-{first}-of-3")
        runner.worker(paused_plan, "parallel-finalize-worker", label="finalize-preparation")
        runner.worker(paused_plan, "parallel-build-worker", "--category", "microbes", expected=75,
                      label="build-microbes-durable-pause")
        from storage_layout import library_root, category_jobs, category_logs
        index = library_root(config) / "index"
        current = category_jobs(config, "microbes") / "checkpoints/CURRENT"
        generation = current.read_text().split("\t")[0]
        checkpoint = read_json(current.parent / generation / "receipt.json")
        ensure(checkpoint["completed_files"] == 1 and not checkpoint["complete"], "Native pause did not commit exactly one input")
        ensure(not (index / "microbes/database.json").exists(), "Partial hash was incorrectly published as complete")

        runner.cli("generate-resume-library", ["library", "--config", base_config, "--with-split",
                                               "--dat-list", "smoke", "--replace", "FALSE", "--submit", "local", "--plan-only"])
        resumed_plan = runner.manifest("parallel_library")
        resumed = read_json(resumed_plan)
        commands = {row["role"]: row for row in resumed["commands"]}
        dependencies = commands["publish"]["depends_on"]
        ensure(set(dependencies) == {"build_"+name for name in CATEGORIES}, "Publish does not wait for all five builds")
        ensure(all(commands["build_"+name]["depends_on"] == "finalize_preparation" for name in CATEGORIES),
               "Category builds are not independent siblings in the generated DAG")
        ensure("pilot" not in commands and commands["split_1"]["depends_on"] == "publish",
               "Default library/split DAG unexpectedly includes a pilot or lost publish dependency")
        # Exercise cache reuse through the second frozen plan, rather than
        # bypassing worker validation or changing a published config.json.
        for first in range(1, 4):
            runner.worker(resumed_plan, "parallel-prepare-worker", "--task-id", first, "--stride", 3,
                          label=f"reuse-prepare-stride-{first}-of-3")
        runner.worker(resumed_plan, "parallel-finalize-worker", label="reuse-finalize-preparation")
        for category in CATEGORIES:
            runner.worker(resumed_plan, "parallel-build-worker", "--category", category,
                          label="build-"+category)
        log = (category_logs(config, "microbes") / "build.log").read_text()
        ensure("METAG_RESUME files=1/2" in log and log.count("METAG_FILE_BEGIN 1/2") == 1,
               "Resume did not skip the committed first input")
        for name in ("taxo.k2d", "opts.k2d"):
            ensure(len({sha256(index / category / name) for category in CATEGORIES}) == 1,
                   "Independent categories have different " + name)
        runner.worker(resumed_plan, "parallel-publish-worker", label="publish-five-categories")
        catalog = read_json(index / "database.json")
        ensure(len(catalog["hashes"]) == 5, "Published catalog lacks five hashes")
        from ordinal_query import resolve_native, classifier_command
        from fastq_stream import fragments
        from partitioned_build import database_snapshot
        snapshot = database_snapshot(config)
        native = resolve_native(config)
        paired_bytes = b"".join(b"".join(line.rstrip(b"\r\n") + b"\n" for line in (*a, *b))
                                for a, b in fragments(root / "reads_raw/smoke/tiny_R1.fastq.gz",
                                                      root / "reads_raw/smoke/tiny_R2.fastq.gz"))
        outputs = []
        for mapped in (False, True):
            command = classifier_command(native, snapshot, dict(config, union_memory_mapping=mapped), True)
            result = subprocess.run(command, input=paired_bytes, capture_output=True, env=runner.environment, timeout=60)
            ensure(result.returncode == 0, result.stderr.decode(errors="replace"))
            outputs.append(sorted(result.stdout.splitlines()))
        ensure(len(outputs[0]) == 6 and outputs[0] == outputs[1], "RAM and mmap classification outputs differ")
        report["ram_mmap_exact_match"] = True
        runner.worker(resumed_plan, "parallel-pilot-worker", label="three-pair-pilot")
        runner.worker(resumed_plan, "parallel-split-worker", "--batch-id", 1, "--offset", 0,
                      label="six-pair-split")
        runner.worker(resumed_plan, "parallel-summary-worker", label="cohort-summary")
        cohort = read_json(root / "work/qc/split/smoke/summary/split-summary.json")
        ensure(cohort["state"] == "complete" and cohort["totals"]["fragments"] == 6,
               "Automatic cohort summary did not conserve synthetic fragments")
        output = root / "output/kraken/smoke/tiny"
        exact = verify_split(output, truth, originals)
        runner.worker(resumed_plan, "parallel-split-worker", "--batch-id", 1, "--offset", 0,
                      label="reuse-complete-split")
        batch = read_json(resumed_plan.parent / "batch.00001.json")
        ensure(not batch["outputs"]["failed"] and len(batch["outputs"]["succeeded"]) == 1 and
               batch["outputs"]["succeeded"][0]["status"] == "skipped", "Completed split was not reused")
        runner.cli("full-index-audit", ["audit", "--config", base_config, "--full",
                                        *(["--allow-local"] if args.allow_local else [])])
        from unittest.mock import patch
        from common import load_config
        from query import run_sample, run_batch, BatchFailure
        from split_summary import load_sample
        from assignment_io import find_assignment_source
        import query as query_module
        real_cfg = load_config(base_config)
        sample = resumed["context"]["split"]["samples"][0]
        recovery_cfg = dict(real_cfg, output_dir=str(root / "recovery"), replace=True,
                            materialize_species_min_fragments=1, run_label="recovery")
        recovery_output = root / "recovery/smoke/tiny"
        actual_replace = os.replace
        def interrupted_publish(source, destination):
            if Path(destination) == recovery_output:
                raise OSError("injected interrupted publication")
            return actual_replace(source, destination)
        with patch("query.os.replace", side_effect=interrupted_publish):
            try:
                run_sample(recovery_cfg, sample, snapshot=snapshot)
                raise AssertionError("Interrupted publication unexpectedly succeeded")
            except OSError as error:
                ensure("injected" in str(error), str(error))
        ensure(not (recovery_output / "qc/complete.json").exists(), "Partial output was marked complete")
        recovered = run_sample(recovery_cfg, sample, snapshot=snapshot)
        load_sample(recovery_output, sample["sample_id"])
        ensure(find_assignment_source(recovery_output).parent.name == "qc", "Wrong assignment metadata location")
        ensure(len(list(recovery_output.rglob("summary.json"))) == 1 and
               not list(recovery_output.rglob("classifier.log")), "Duplicate metadata/runtime log in data")
        ensure((recovery_output / "bacteria/10_Synthetic_taxon_10.R1.fastq.gz").is_file(),
               "Materialized species is outside the flat group layout")
        before = read_json(recovery_output / "qc/complete.json")["completed_utc"]
        overwritten = run_sample(recovery_cfg, sample, snapshot=snapshot)
        ensure(overwritten["status"] == "complete" and
               read_json(recovery_output / "qc/complete.json")["completed_utc"] != before,
               "Default replacement did not rebuild")
        reuse_cfg = dict(recovery_cfg, replace=False)
        ensure(run_sample(reuse_cfg, sample, snapshot=snapshot)["status"] == "skipped", "Explicit FALSE did not reuse")
        fastq = recovery_output / "bacteria/10_Synthetic_taxon_10.R1.fastq.gz"
        fastq.unlink()
        ensure(run_sample(reuse_cfg, sample, snapshot=snapshot)["status"] == "complete",
               "Missing artifact was not rebuilt under replace FALSE")
        # The same read IDs deliberately recur in two different samples.
        second = dict(sample, sample="tiny2", sample_id="smoke/tiny2")
        batch_cfg = dict(real_cfg, output_dir=str(root / "batch"), run_label="batch", replace=True,
                         classify_workers=1, batch_query=True)
        batch_started = time.monotonic()
        batch_result = run_batch(batch_cfg, [sample, second], snapshot=snapshot)
        ensure(batch_result["native_invocations"] == 1 and len(batch_result["succeeded"]) == 2,
               "Resident batch did not share one native invocation")
        for item in (sample, second):
            verify_split(root / "batch" / item["sample_id"], truth, originals)
        ensure(batch_result["succeeded"][0]["summary"]["batch_samples"] == 2, "Missing batch timing scope")
        cache_info = batch_result["succeeded"][0]["summary"]["project_preparation"]
        ensure(cache_info["taxonomy_cache_hit"], "Later batch did not reuse compiled reference taxonomy")
        cache_path = Path(cache_info["taxonomy_cache_path"])
        cache_mtime = cache_path.stat().st_mtime_ns
        with patch("project_cache._compile", side_effect=AssertionError("Unchanged project cache rebuilt")):
            next_batch = run_batch(batch_cfg, [sample, second], snapshot=snapshot)
        ensure(next_batch["native_invocations"] == 1 and
               all(item["status"] == "complete" for item in next_batch["succeeded"]),
               "Default replacement must rerun classification while retaining project preparation")
        ensure(cache_path.stat().st_mtime_ns == cache_mtime, "Sample overwrite replaced taxonomy cache")
        for item in (sample, second):
            verify_split(root / "batch" / item["sample_id"], truth, originals)
        report["project_cache"] = {
            "reused_across_batches": True, "default_overwrite_reclassified_every_sample": True,
            "taxonomy_cache_rows": cache_info["taxonomy_cache_rows"],
            "initial_build_seconds": read_json(cache_path.parent / "complete.json")["build_seconds"],
            "later_preparation_seconds": next_batch["succeeded"][0]["summary"]["project_preparation"]["project_preparation_seconds"],
            "native_invocations_per_batch": [batch_result["native_invocations"], next_batch["native_invocations"]]}
        # A failure in the second sample's split preserves the first commit.
        failed_cfg = dict(batch_cfg, output_dir=str(root / "partial"), run_label="partial")
        actual_split = query_module.split_originals
        def fail_second(cfg, current, *rest, **kwargs):
            if current["sample_id"] == second["sample_id"]:
                raise RuntimeError("injected later sample split failure")
            return actual_split(cfg, current, *rest, **kwargs)
        with patch("query.split_originals", side_effect=fail_second):
            try:
                run_batch(failed_cfg, [sample, second], snapshot=snapshot)
                raise AssertionError("Second sample split failure was ignored")
            except BatchFailure as error:
                ensure(len(error.summary["succeeded"]) == 1 and len(error.summary["failed"]) == 1,
                       "First committed sample was lost")
        repaired = run_batch(dict(failed_cfg, replace=False), [sample, second], snapshot=snapshot)
        ensure([x["status"] for x in repaired["succeeded"]] == ["skipped", "complete"],
               "Retry did not reuse complete/rebuild incomplete samples")
        report["layout_checks"] = {"interrupted_publication_recovery": True, "flat_layout": True,
                                  "single_qc_metadata": True, "missing_artifact_rebuilt": True,
                                  "default_replace": True, "explicit_reuse": True,
                                  "later_split_failure_preserves_completed_sample": True}
        report["resident_batch"] = {"samples": 2, "native_index_loads": batch_result["native_invocations"],
                                    "seconds": batch_result["seconds"],
                                    "separate_query_native_index_loads": 2,
                                    "separate_query_classification_seconds": recovered["summary"]["classification_seconds"] + overwritten["summary"]["classification_seconds"],
                                    "batch_query_classification_seconds": batch_result["succeeded"][0]["summary"]["batch_classification_seconds"],
                                    "original_bytes_and_assignments_equal": True,
                                    "index_load_seconds": batch_result["succeeded"][0]["summary"].get("batch_index_load_seconds"),
                                    "note": "Tiny fixture measures invocation reduction, not production throughput"}
        report.update(status="pass", reference_categories=list(CATEGORIES), preparation_tasks=10,
                      preparation_worker_stride=3, reference_bases=250000,
                      native_pause_exit=75, native_committed_files_before_resume=1,
                      native_resume_skips_committed_input=True, shared_taxonomy_and_options=True,
                      all_five_build_dependencies=True, complete_split_reuse="skipped",
                      split=exact, catalog_identity=catalog["identity"], output=str(output),
                      manifests={"setup": str(setup), "pause_library": str(paused_plan), "resume_library": str(resumed_plan)},
                      note="Synthetic correctness/operability smoke test; not a full-reference speed or sensitivity benchmark")
    except BaseException as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        report.update(completed_at=utcnow(), elapsed_seconds=time.monotonic()-started, steps=runner.steps)
        write_json(root / "smoke-result.json", report)
        print(json.dumps(report, indent=2, sort_keys=True), flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--tools-dir", type=Path, help="Existing verified tools directory; omitted means compile under destination/tools")
    parser.add_argument("--code-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--storage-layout", choices=("legacy", "separated"), default="separated")
    parser.add_argument("--allow-local", action="store_true", help="Explicit small local test; otherwise require an LSF compute job")
    args = parser.parse_args()
    try:
        return smoke(args)
    except (OSError, RuntimeError, ValueError, KeyError, AssertionError) as error:
        print(f"SMOKE FAILED: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
