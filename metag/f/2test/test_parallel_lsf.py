"""Exercise the real frozen DAG and submission receipts without an HPC account.

Only scheduler subprocess calls (bsub/bjobs/bkill) are mocked. Plan generation,
configuration loading, metadata checks, command rendering, fingerprints and
recovery decisions use their real implementations. No jobs leave this process.
"""
import contextlib
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_cli
import hpc_common
import hpc_lsf
import hpc_reference
import parallel_cli
import partitioned_build
import partitioned_reference


class Scheduler:
    def __init__(self):
        self.submissions = []
        self.states = {}
        self.kills = []
        self.queries = []
        self.failure_at = None
        self.failure = None
        self.next_id = 17000
        self.on_submission = None

    def run(self, argv, **kwargs):
        if argv[0] == "bsub":
            self.next_id += 1
            job = self.next_id
            script = kwargs["stdin"].read().decode()
            self.submissions.append({"job_id": job, "argv": argv[:], "script": script})
            self.states[job] = ["DONE"]
            if self.on_submission is not None:
                self.on_submission(self)
            output = f"Job <{job}> is submitted to queue <test>.\n"
            if len(self.submissions) == self.failure_at:
                if self.failure == "accepted_error":
                    return subprocess.CompletedProcess(argv, 1, output, "connection closed after acceptance")
                if self.failure == "unknown_timeout":
                    raise subprocess.TimeoutExpired(argv, kwargs["timeout"], output="", stderr="network timeout")
                if self.failure == "known_timeout":
                    raise subprocess.TimeoutExpired(argv, kwargs["timeout"], output=output, stderr="network timeout")
                raise AssertionError(f"Unexpected failure mode: {self.failure}")
            return subprocess.CompletedProcess(argv, 0, output, "")
        if argv[0] == "bjobs":
            ids = [int(value) for value in argv[5:]]
            self.queries.append(ids)
            lines = []
            for job in ids:
                states = self.states.get(job, [])
                for number, state in enumerate(states, 1):
                    label = f"{job}[{number}]" if len(states) > 1 else str(job)
                    lines.append(f"{label} {state}\n")
            return subprocess.CompletedProcess(argv, 0, "".join(lines), "")
        if argv[0] == "bkill":
            ids = [int(value) for value in argv[1:]]
            self.kills.append(ids)
            for job in ids:
                self.states[job] = ["EXIT"]
            return subprocess.CompletedProcess(argv, 0, "".join(f"Job <{job}> is being terminated\n" for job in ids), "")
        raise AssertionError(f"Unexpected external command under scheduler mock: {argv}")


class ParallelSchedulingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.environ = patch.dict(os.environ, {"submit": "hpc", "MEMORY_DEFAULT_GIB": "32"})
        self.environ.start()
        self.addCleanup(self.environ.stop)
        values = {key: str(self.root / key) for key in hpc_common.PATH_KEYS}
        values.update(schema=1, root=str(self.root), database_layout="partitioned", classifier_engine="kraken",
                      legacy_library_dir=None, legacy_metadata_dir=None, parallel_reference_dir=None,
                      reference_taxid_map=None, integrity_seal=None, build_kunpeng=False,
                      convert_memory_gib=32, min_free_gib=0, masking=False,
                      prepare_concurrency=3, max_array_size=4, sample_batch_size=1,
                      prepare_threads=2, prepare_cores=4, prepare_memory_gib=16, prepare_hours=24,
                      prepare_queue="prepare_q", setup_queue="setup_q", setup_cores=4,
                      setup_memory_gib=16, setup_hours=24,
                      build_cores=96, build_threads=88, build_memory_gib=4608, build_hours=120,
                      build_queue="build_q", classify_cores=16, classify_threads=12,
                      classify_workers=1, classify_concurrency=2, classify_queue="query_q",
                      classify_memory_gib=4608, classify_hours=120)
        values["category_resources"] = {
            "microbes": dict(queue="microbes_q", cores=12, threads=10, memory_gib=96, hours=96),
            "viruses": dict(queue="viruses_q", cores=4, threads=4, memory_gib=16, hours=48),
            "plants": dict(queue="plants_q", cores=48, threads=40, memory_gib=1024, hours=110),
            "invertebrates": dict(queue="invertebrates_q", cores=32, threads=32, memory_gib=512, hours=114),
            "vertebrates": dict(queue="vertebrates_q", cores=96, threads=88, memory_gib=4608, hours=120),
        }
        config = self.root / "config.json"
        hpc_common.atomic_json(config, values)
        self.cfg = hpc_common.load_config(config)
        self._reference_plan()

    def _reference_plan(self):
        """A structurally real tiny frozen reference plan; no native tool needed."""
        raw = Path(self.cfg["reference_raw"])
        metadata = Path(self.cfg["metadata_dir"])
        taxonomy = metadata / "taxonomy"
        taxonomy.mkdir(parents=True)
        (taxonomy / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n")
        (taxonomy / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        (taxonomy / "merged.dmp").write_text("")
        for category in partitioned_reference.CATEGORIES:
            for number in range(3):
                path = raw / category / f"ref{number}.fna"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(">test\nACGTACGTACGT\n")
        inventory = hpc_reference._inventory(self.cfg)
        mapping = metadata / "reference_taxids.tsv"
        mapping.write_text("path\taccession\ttaxid\tsource\n" + "".join(
            f"{row['path']}\t\t1\ttest\n" for row in inventory))
        records = {name: {"state": hpc_common.file_state(metadata / name)} for name in (
            "taxonomy/nodes.dmp", "taxonomy/names.dmp", "taxonomy/merged.dmp", "reference_taxids.tsv")}
        meta_context = {"raw_inventory": hpc_common.fingerprint(inventory), "custom_map": None, "files": records}
        meta = {"schema": hpc_reference.SCHEMA, "status": "complete", "context": meta_context,
                "identity": hpc_common.fingerprint(meta_context)}
        hpc_common.atomic_json(metadata / "metadata.json", meta)
        base = Path(self.cfg["library_dir"]) / "parallel"
        tasks = []
        for index, item in enumerate(inventory):
            category = item["path"].split("/", 1)[0]
            task_id = index + 1
            path = base / "manifests" / f"task{task_id:06d}.json"
            hpc_common.atomic_json(path, {"task_id": task_id,
                "files": [{**item, "index": index, "taxid": 1, "category": category}], "legacy": None})
            tasks.append({"task_id": task_id, "mode": "prepare", "categories": [category],
                          "compressed_bytes": item["state"]["size"], "manifest": {
                              "path": path.relative_to(base).as_posix(), "size": path.stat().st_size,
                              "sha256": hpc_common.sha256(path)}})
        context = {"raw_inventory": hpc_common.fingerprint(inventory), "metadata_identity": meta["identity"],
                   "mapping": {"size": mapping.stat().st_size, "sha256": hpc_common.sha256(mapping)},
                   "policy": {"masking": False, "record_policy": partitioned_reference.POLICY,
                              "source_commit": hpc_reference.SOURCE_COMMIT, "source_tag": "2.17.2"},
                   "tasks": tasks, "categories": list(partitioned_reference.CATEGORIES),
                   "category_overrides": {}, "legacy_plan_identity": None}
        self.reference = {"schema": partitioned_reference.SCHEMA, "status": "complete", "context": context,
                          "identity": hpc_common.fingerprint(context), "tasks": tasks,
                          "categories": list(partitioned_reference.CATEGORIES), "references": len(inventory)}
        hpc_common.atomic_json(base / "plan.json", self.reference)

    def generate(self, *, with_split=False):
        if with_split:
            reads = Path(self.cfg["reads_raw"]) / "dataset"
            reads.mkdir(parents=True, exist_ok=True)
            for number in range(9):
                (reads / f"sample{number}.fastq").write_text("@read\nACGT\n+\nIIII\n")
        with contextlib.redirect_stdout(io.StringIO()):
            return parallel_cli.generate_library(self.cfg, with_split=with_split)

    def submit(self, path, scheduler, **kwargs):
        with patch.object(subprocess, "run", side_effect=scheduler.run), contextlib.redirect_stdout(io.StringIO()):
            return parallel_cli.submit_plan(path, **kwargs)

    @staticmethod
    def load(path):
        return json.loads(Path(path).read_text())

    def test_five_category_builds_parallel_and_actual_resource_requests_match(self):
        path = self.generate()
        manifest = self.load(path)
        by_role = {command["role"]: command for command in manifest["commands"]}
        self.assertEqual(len(by_role), 8)
        self.assertEqual(by_role["finalize_preparation"]["depends_on"], "prepare")
        self.assertEqual(set(by_role["publish"]["depends_on"]), {"build_" + c for c in partitioned_build.CATEGORIES})
        frozen = hpc_common.load_config(manifest["config"])
        for category in partitioned_build.CATEGORIES:
            command = by_role["build_" + category]
            self.assertEqual(command["depends_on"], "finalize_preparation")
            lines = Path(command["path"]).read_text().splitlines()
            profile = self.cfg["category_resources"][category]
            self.assertIn(f"#BSUB -q {profile['queue']}", lines)
            self.assertIn(f"#BSUB -n {profile['cores']}", lines)
            self.assertIn(f"#BSUB -W {profile['hours']}:00", lines)
            self.assertTrue(any(f"rusage[mem={profile['memory_gib']}GB" in line for line in lines))
            invocation = shlex.split(next(line for line in lines if line.startswith("exec ")))[1:]
            self.assertIn("parallel-build-worker", invocation)
            self.assertEqual(invocation[invocation.index("--category") + 1], category)
            self.assertEqual(invocation[invocation.index("--config") + 1], manifest["config"])
            self.assertEqual(partitioned_build.resources(frozen, category)["build_threads"], profile["threads"])
            self.assertLessEqual(profile["threads"], profile["cores"])
            result = subprocess.run(["bash", "-n", command["path"]], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        prepare_text = Path(by_role["prepare"]["path"]).read_text()
        self.assertIn('"metag_prepare[1-3]%3"', prepare_text)
        self.assertIn('--task-id "${LSB_JOBINDEX:?}" --stride 3', prepare_text)
        assigned = [task for worker in range(1, 4) for task in hpc_cli._prepare_task_ids(worker, 15, 3)]
        self.assertEqual(sorted(assigned), list(range(1, 16)))
        self.assertEqual(len(assigned), len(set(assigned)))

    def test_pilot_and_split_wait_for_successful_publish_and_prior_arrays(self):
        path = self.generate(with_split=True)
        manifest = self.load(path)
        commands = {c["role"]: c for c in manifest["commands"]}
        self.assertEqual(commands["pilot"]["depends_on"], "publish")
        self.assertEqual(commands["split_1"]["depends_on"], "pilot")
        self.assertEqual(commands["split_2"]["depends_on"], "split_1")
        self.assertEqual(commands["split_3"]["depends_on"], "split_2")
        split = manifest["context"]["split"]
        self.assertEqual(sorted(i for batch in split["batches"] for i in batch), list(range(9)))
        self.assertEqual(len(split["pilot_indices"]), 1)
        self.assertIn("--offset 4", Path(commands["split_2"]["path"]).read_text())
        self.assertIn("--offset 8", Path(commands["split_3"]["path"]).read_text())

    def test_submit_uses_numeric_done_dependencies_and_keeps_receipts(self):
        path = self.generate()
        scheduler = Scheduler()
        history = self.submit(path, scheduler)
        records = history["attempts"][0]["jobs"]
        by_role = {record["role"]: record for record in records}
        self.assertEqual(history["attempts"][0]["status"], "submitted")
        self.assertIsNone(by_role["prepare"]["dependency"])
        self.assertEqual(by_role["finalize_preparation"]["dependency"], f"done({by_role['prepare']['job_id']})")
        for category in partitioned_build.CATEGORIES:
            self.assertEqual(by_role["build_" + category]["dependency"], f"done({by_role['finalize_preparation']['job_id']})")
        expected = " && ".join(f"done({by_role['build_' + c]['job_id']})" for c in partitioned_build.CATEGORIES)
        self.assertEqual(by_role["publish"]["dependency"], expected)
        self.assertEqual(len(scheduler.submissions), 8)
        self.assertEqual(self.load(path.parent / "submission.json"), history)
        self.assertTrue(all(Path(record["receipt_path"]).is_file() for record in records))

    def test_second_submission_requires_explicit_retry(self):
        path, scheduler = self.generate(), Scheduler()
        self.submit(path, scheduler)
        with self.assertRaisesRegex(RuntimeError, "already has submission receipts"):
            self.submit(path, scheduler)
        self.assertEqual(len(scheduler.submissions), 8)

    def test_known_id_survives_partial_acceptance_error_and_can_retry(self):
        path, scheduler = self.generate(), Scheduler()
        scheduler.failure_at, scheduler.failure = 3, "accepted_error"
        with self.assertRaises(hpc_lsf.SubmissionError):
            self.submit(path, scheduler)
        previous = self.load(path.parent / "submission.json")["attempts"][0]
        self.assertEqual(previous["status"], "submission_failed")
        self.assertEqual(len(previous["jobs"]), 3)
        self.assertEqual(previous["jobs"][-1]["job_id"], 17003)
        scheduler.states[17003] = ["EXIT"]
        history = self.submit(path, scheduler, retry=True)
        self.assertEqual(len(history["attempts"]), 2)
        self.assertEqual(scheduler.queries[0], [17001, 17002, 17003])
        self.assertEqual(len(scheduler.submissions), 11)

    def test_unknown_acceptance_blocks_any_automatic_retry(self):
        path, scheduler = self.generate(), Scheduler()
        scheduler.failure_at, scheduler.failure = 2, "unknown_timeout"
        with self.assertRaises(hpc_lsf.SubmissionError):
            self.submit(path, scheduler)
        receipt = self.load(path.parent / "submission.json")["attempts"][0]["jobs"][-1]
        self.assertEqual(receipt["status"], "unknown")
        self.assertIsNone(receipt["job_id"])
        with self.assertRaisesRegex(RuntimeError, "acceptance is unknown"):
            self.submit(path, scheduler, retry=True, cancel_pending=True)
        self.assertFalse(scheduler.queries)
        self.assertFalse(scheduler.kills)
        self.assertEqual(len(scheduler.submissions), 2)

    def test_timeout_with_known_id_is_verified_before_retry(self):
        path, scheduler = self.generate(), Scheduler()
        scheduler.failure_at, scheduler.failure = 2, "known_timeout"
        with self.assertRaises(hpc_lsf.SubmissionError):
            self.submit(path, scheduler)
        scheduler.states[17002] = ["EXIT"]
        history = self.submit(path, scheduler, retry=True)
        self.assertEqual(scheduler.queries[0], [17001, 17002])
        self.assertEqual(history["attempts"][-1]["status"], "submitted")

    def test_paused_build_exit_blocks_old_publish_until_explicit_pending_recovery(self):
        path, scheduler = self.generate(), Scheduler()
        self.submit(path, scheduler)
        scheduler.states[17003] = ["EXIT"]  # Native exit 75 is an LSF EXIT, never DONE.
        scheduler.states[17008] = ["PEND"]
        with self.assertRaisesRegex(RuntimeError, "still wait on old IDs"):
            self.submit(path, scheduler, retry=True)
        self.assertEqual(len(scheduler.submissions), 8)
        self.assertFalse(scheduler.kills)
        history = self.submit(path, scheduler, retry=True, cancel_pending=True)
        self.assertEqual(scheduler.kills, [[17008]])
        self.assertEqual(len(scheduler.submissions), 16)
        jobs = {r["role"]: r for r in history["attempts"][-1]["jobs"]}
        self.assertTrue(all(f"done({jobs['build_' + category]['job_id']})" in jobs["publish"]["dependency"]
                            for category in partitioned_build.CATEGORIES))
        self.assertNotIn("ended(", jobs["publish"]["dependency"])

    def test_running_array_child_refuses_retry_and_never_cancels_any_job(self):
        path, scheduler = self.generate(), Scheduler()
        self.submit(path, scheduler)
        scheduler.states[17001] = ["DONE", "RUN", "PEND"]
        scheduler.states[17002] = ["PEND"]
        with self.assertRaisesRegex(RuntimeError, "still running"):
            self.submit(path, scheduler, retry=True, cancel_pending=True)
        self.assertFalse(scheduler.kills)
        self.assertEqual(len(scheduler.submissions), 8)

    def test_expired_scheduler_history_refuses_blind_retry(self):
        path, scheduler = self.generate(), Scheduler()
        self.submit(path, scheduler)
        del scheduler.states[17003]
        with self.assertRaisesRegex(RuntimeError, "expired"):
            self.submit(path, scheduler, retry=True, cancel_pending=True)
        self.assertFalse(scheduler.kills)
        self.assertEqual(len(scheduler.submissions), 8)

    def test_cancel_pending_needs_retry_flag(self):
        path, scheduler = self.generate(), Scheduler()
        with self.assertRaisesRegex(ValueError, "only valid with --retry"):
            self.submit(path, scheduler, cancel_pending=True)
        self.assertFalse(scheduler.submissions)

    def test_changed_command_fails_before_bsub(self):
        path, scheduler = self.generate(), Scheduler()
        command = Path(self.load(path)["commands"][2]["path"])
        command.write_text(command.read_text() + "echo changed\n")
        with self.assertRaisesRegex(RuntimeError, r"\.cmd changed"):
            self.submit(path, scheduler)
        self.assertFalse(scheduler.submissions)

    def test_later_command_changed_during_submission_is_never_sent_to_bsub(self):
        path, scheduler = self.generate(), Scheduler()
        command = Path(self.load(path)["commands"][1]["path"])
        original = command.read_text()
        def mutate_next_script(current):
            if len(current.submissions) == 1:
                command.write_text(original + "echo concurrent-change\n")
        scheduler.on_submission = mutate_next_script
        with self.assertRaisesRegex(hpc_lsf.SubmissionError, "[Cc]hanged|[Ss][Hh][Aa]|[Cc]hecksum"):
            self.submit(path, scheduler)
        self.assertEqual(len(scheduler.submissions), 1)
        attempt = self.load(path.parent / "submission.json")["attempts"][-1]
        self.assertEqual(attempt["status"], "submission_failed")
        self.assertEqual(len(attempt["jobs"]), 2)
        self.assertEqual(attempt["jobs"][-1]["status"], "failed")
        self.assertIsNone(attempt["jobs"][-1]["job_id"])
        # This rejection happened before bsub. After restoring the reviewed
        # bytes, automatic retry may safely inspect the one accepted ID.
        command.write_text(original)
        scheduler.on_submission = None
        history = self.submit(path, scheduler, retry=True)
        self.assertEqual(scheduler.queries[0], [17001])
        self.assertEqual(len(scheduler.submissions), 9)
        self.assertEqual(history["attempts"][-1]["status"], "submitted")

    def test_changed_frozen_config_fails_before_bsub(self):
        path, scheduler = self.generate(), Scheduler()
        config = Path(self.load(path)["config"])
        value = self.load(config)
        value["category_resources"]["microbes"]["threads"] = 9
        hpc_common.atomic_json(config, value)
        with self.assertRaisesRegex(RuntimeError, "Frozen job configuration changed"):
            self.submit(path, scheduler)
        self.assertFalse(scheduler.submissions)

    def test_changed_manifest_fails_before_bsub(self):
        path, scheduler = self.generate(), Scheduler()
        value = self.load(path)
        value["commands"][-1]["depends_on"] = ["build_microbes"]
        hpc_common.atomic_json(path, value)
        with self.assertRaisesRegex(RuntimeError, "manifest changed"):
            self.submit(path, scheduler)
        self.assertFalse(scheduler.submissions)

    def test_threads_above_allocated_cores_are_refused_before_generation(self):
        self.cfg["category_resources"]["microbes"]["threads"] = 100
        with self.assertRaisesRegex(ValueError, "threads exceeds"):
            self.generate()


if __name__ == "__main__":
    unittest.main()
