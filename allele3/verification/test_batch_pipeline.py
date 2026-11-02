"""Small real payloads plus stub scientific stages; no genomic data or LSF jobs."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout, redirect_stderr
import argparse
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
import batch_pipeline as bp
import batch_runner as br


class StubPipeline:
    VERSION = "fixture-1"
    STAGES = ("primary", "archaic", "summary", "alleles", "time", "selection", "candidates", "genotypes", "report")
    now = staticmethod(bp.now)

    def __init__(self):
        self.calls = []
        self.plans = []
        self.results = {}
        self.failures = {}
        self.missing = {}
        self.metadata_calls = 0
        self.barrier = None

    @staticmethod
    def parser():
        p = argparse.ArgumentParser()
        p.add_argument("command", nargs="?", default="pipeline")
        p.add_argument("--out", required=True)
        p.add_argument("--targets")
        p.add_argument("--allele-catalog")
        p.add_argument("--label", default="fixture")
        p.add_argument("--dry-run", action="store_true")
        p.add_argument("--preflight", action="store_true")
        return p

    @staticmethod
    def resolve(args):
        out = Path(args.out).resolve()
        return dict(out=out, root=out.parent, metadata=out / "metadata" / "samples.tsv",
                    metadata_mode="prepare_from_original_supplement")

    def plan(self, args, paths, names):
        self.plans.extend(names)
        return [SimpleNamespace(name=name, inputs=[], out=paths["out"] if name == "report" else paths["out"] / name)
                for name in names]

    def availability(self, stage):
        return dict(missing_inputs=self.missing.get(stage.name, []), missing_packages=[], missing_programs=[])

    @staticmethod
    def snapshot(stage):
        return dict(schema="fixture-1", command=["fixture", stage.name],
                    inputs=[br.path_identity(path) for path in stage.inputs])

    def prepare_metadata_if_needed(self, paths, receipts):
        self.metadata_calls += 1
        paths["metadata"].parent.mkdir(exist_ok=True)
        paths["metadata"].write_text("id\tage_bp\nsample\t1000\n")

    def run_stage(self, stage, receipts, always=False):
        self.calls.append((stage.name, always, list(stage.inputs)))
        if self.barrier:
            self.barrier.wait(timeout=5)
        if stage.name in self.failures:
            raise self.failures[stage.name]
        if self.missing.get(stage.name):
            return dict(stage=stage.name, status="not_available", missing_inputs=self.missing[stage.name])
        status = self.results.get(stage.name, "complete")
        if status in {"complete", "reused"}:
            stage.out.mkdir(parents=True, exist_ok=True)
            target = stage.out / ("report.html" if stage.name == "report" else "result.txt")
            subprocess.run([sys.executable, "-c", "from pathlib import Path; import sys; Path(sys.argv[1]).write_text(sys.argv[2])",
                            str(target), stage.name], check=True)
            if stage.name == "report":
                seen = json.loads((stage.out / "pipeline_status.json").read_text())
                (stage.out / "report_seen.json").write_text(json.dumps(seen))
            else:
                (receipts / (stage.name + ".json")).write_text(json.dumps(dict(
                    identity=self.snapshot(stage), outputs={str(target): br.path_identity(target)}, result={"status": status})))
        return dict(stage=stage.name, status=status, result_status="complete", completed=bp.now())

    @staticmethod
    def overall(records):
        return "complete_requested_stages" if all(r["status"] in bp.SUCCESS for r in records) else "unexpected"

    @staticmethod
    def execute(args):
        raise AssertionError("The batch wrapper must never call serial execute()")


class BatchPipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="allele3-stage-qa-")
        self.path = Path(self.temp.name)
        self.out = self.path / "results with spaces"
        self.batch = "a" * 64
        self.stub = StubPipeline()
        self.env = patch.dict(os.environ, {"ALLELE_BATCH_WORKER": "1", "WORKFLOW_WORKER": "1",
                                          "ALLELE_EXECUTION_MODE": "hpc", "LSB_JOBID": "123"})
        self.env.start()
        self.importer = patch.object(bp, "load_pipeline", return_value=self.stub)
        self.importer.start()
        self.audit = patch.object(bp, "audit_stage_sources", return_value=dict(has_conflicts=False, datasets=[]))
        self.audit_mock = self.audit.start()

    def tearDown(self):
        self.audit.stop()
        self.importer.stop()
        self.env.stop()
        self.temp.cleanup()

    def argv(self, *extra):
        return ["pipeline", "--out", str(self.out), *extra]

    def run_stage(self, name, selected=None, batch=None, extra=()):
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            return bp.run_stage_worker(name, selected or name, batch or self.batch, self.argv(*extra))

    def record(self, name, batch=None):
        return json.loads((self.out / ".batch_stage_status" / (batch or self.batch) / (name + ".json")).read_text())

    def test_one_stage_runs_payload_and_writes_its_own_receipt(self):
        self.assertEqual(self.run_stage("primary"), 0)
        row = self.record("primary")
        self.assertEqual(row["status"], "complete")
        self.assertEqual(row["lsf_job_id"], "123")
        self.assertEqual(row["worker_exit_code"], 0)
        self.assertTrue((self.out / ".receipts/primary.json").is_file())
        self.assertEqual((self.out / "primary/result.txt").read_text(), "primary")
        self.assertFalse((self.out / "pipeline_status.json").exists())
        self.assertEqual(self.stub.plans, ["primary"])

    def test_stage_ledger_is_bound_to_actual_products_and_blocks_stale_dependency(self):
        self.assertEqual(self.run_stage("primary", "primary,report"), 0)
        ledger=self.out / ".batch_stage_status" / self.batch / "primary.json"
        br.validate_output_file(ledger)
        before=ledger.stat().st_mtime_ns
        (self.out / "primary/result.txt").unlink()
        with self.assertRaisesRegex(ValueError, "Scientific file identity changed"):
            br.validate_output_file(ledger)
        self.assertEqual(ledger.stat().st_mtime_ns,before)
        self.assertEqual(self.run_stage("report", "primary,report"), 2)
        self.assertFalse((self.out / "report.html").exists())

    def test_stage_ledger_checks_input_snapshot_and_report_output(self):
        self.assertEqual(self.run_stage("time", "time,report"), 0)
        self.assertEqual(self.run_stage("report", "time,report"), 0)
        ledger=self.out / ".batch_stage_status" / self.batch / "report.json"
        br.validate_output_file(ledger)
        (self.out / "metadata/samples.tsv").write_text("changed dates\n")
        with self.assertRaisesRegex(ValueError, "Scientific file identity changed"):
            br.validate_output_file(ledger)

    def test_old_ledger_without_scientific_binding_is_not_accepted_as_verified(self):
        path=self.path / "old_stage.json"
        path.write_text(json.dumps(dict(schema=bp.SCHEMA,status="complete")))
        with self.assertRaisesRegex(ValueError, "lacks verified scientific outputs"):
            br.validate_output_file(path)

    def test_report_sees_only_this_batch_and_final_state_contains_report(self):
        self.assertEqual(self.run_stage("primary", "primary,report"), 0)
        foreign = self.out / ".batch_stage_status" / ("b" * 64)
        foreign.mkdir()
        (foreign / "ghost.json").write_text('{"stage":"ghost","status":"complete"}')
        (self.out / "pipeline_status.json").write_text('{"stages":[{"stage":"ghost","status":"complete"}]}')
        self.assertEqual(self.run_stage("report", "primary,report"), 0)
        before = json.loads((self.out / "report_seen.json").read_text())
        after = json.loads((self.out / "pipeline_status.json").read_text())
        self.assertEqual(before["batch_id"], self.batch)
        self.assertEqual([r["stage"] for r in before["stages"]], ["primary"])
        self.assertEqual(before["status"], "complete_requested_stages")
        self.assertEqual(after["status"], "complete_requested_stages")
        self.assertEqual(after["stages"][1]["status"], "complete")
        self.assertTrue(self.stub.calls[-1][1])

    def test_failed_missing_dependency_does_not_use_stale_catalog(self):
        (self.out / "alleles").mkdir(parents=True)
        (self.out / "alleles/site_catalog.tsv.gz").write_text("old output")
        self.assertEqual(self.run_stage("time", "alleles,time"), 2)
        self.assertEqual(self.record("time")["status"], "failed")
        self.assertIn("Dependency alleles is not_started", self.record("time")["reason"])
        self.assertEqual(self.stub.calls, [])

    def test_report_fails_if_selected_stage_has_no_this_batch_receipt(self):
        self.assertEqual(self.run_stage("report", "primary,report"), 2)
        final = json.loads((self.out / "pipeline_status.json").read_text())
        self.assertEqual(final["status"], "failed")
        self.assertEqual(final["stages"][0]["status"], "not_started")
        self.assertFalse((self.out / "report.html").exists())

    def test_explicit_catalog_does_not_require_an_unrelated_selected_alleles_job(self):
        self.assertEqual(self.run_stage("selection", "alleles,selection", extra=("--allele-catalog", "/explicit/catalog")), 0)

    def test_failed_and_unavailable_results_exit_nonzero(self):
        for status, code in (("failed", 2), ("not_available", 3)):
            self.stub.results["primary"] = status
            self.assertEqual(self.run_stage("primary"), code)
            row = self.record("primary")
            self.assertEqual(row["status"], status)
            self.assertEqual(row["worker_exit_code"], code)

    def test_exception_is_recorded_and_returns_nonzero(self):
        self.stub.failures["primary"] = ValueError("scientific input changed")
        self.assertEqual(self.run_stage("primary"), 2)
        row = self.record("primary")
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["error_type"], "ValueError")
        self.assertIn("scientific input changed", row["reason"])

    def test_source_conflict_is_audited_before_the_stage_plan_or_payload(self):
        self.audit_mock.return_value = dict(has_conflicts=True, datasets=[dict(id="fixture", expected_build="38", observed_build="37")])
        self.assertEqual(self.run_stage("primary"), 2)
        row = self.record("primary")
        self.assertTrue(row["source_audit"]["has_conflicts"])
        self.assertIn("Conflicting source identities", row["reason"])
        self.assertEqual(self.stub.plans, [])
        self.assertEqual(self.stub.calls, [])

    def test_actual_source_audit_is_narrowed_and_bounded_without_mutating_argv(self):
        from unittest.mock import Mock
        source = SimpleNamespace(DEFAULT_CONFIG=Path("fixture.tsv"), load_config=Mock(return_value=[{"id": "fixture"}]),
                                 audit_entries=Mock(return_value=dict(has_conflicts=False)))
        workflow = SimpleNamespace(actual_sources=Mock(return_value=([{"id": "fixture"}], "/ref37", "/ref38")))
        argv = self.argv("--stages", "primary,time,report")
        original = list(argv)
        with patch.dict(sys.modules, data_sources=source, workflow_data_check=workflow):
            # Call the original helper, bypassing the default fixture mock.
            self.audit.stop()
            try:
                result = bp.audit_stage_sources("time", argv)
            finally:
                self.audit_mock = self.audit.start()
        self.assertFalse(result["has_conflicts"])
        self.assertEqual(argv, original)
        self.assertEqual(workflow.actual_sources.call_args.args[1][-2:], ["--stages", "time"])
        source.audit_entries.assert_called_once_with([{"id": "fixture"}], out=None, fasta37="/ref37", fasta38="/ref38",
                                                    max_records=1000, max_sites=100)

    def test_same_stage_and_legacy_pipeline_locks_refuse_duplicate_payload(self):
        for path in (self.out / ".batch_stage_locks/primary.lock", self.out / ".pipeline.lock"):
            with bp.file_lock(path):
                with self.assertRaisesRegex(RuntimeError, "no duplicate stage"):
                    self.run_stage("primary")
            self.assertEqual(self.stub.calls, [])

    def test_different_stage_payloads_can_run_simultaneously(self):
        self.stub.barrier = threading.Barrier(2)
        # A real barrier inside run_stage would time out if .pipeline.lock were
        # exclusive or the owner lock were held while either payload executes.
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()), ThreadPoolExecutor(2) as pool:
            futures = [pool.submit(bp.run_stage_worker, name, "primary,candidates", self.batch, self.argv())
                       for name in ("primary", "candidates")]
            codes = [future.result(timeout=8) for future in futures]
        self.assertEqual(codes, [0, 0])
        self.assertEqual(self.record("primary")["status"], "complete")
        self.assertEqual(self.record("candidates")["status"], "complete")

    def test_unfinished_batch_keeps_output_directory(self):
        self.assertEqual(self.run_stage("primary", "primary,report"), 0)
        with self.assertRaisesRegex(RuntimeError, "unfinished batch"):
            self.run_stage("candidates", batch="b" * 64)
        self.assertFalse((self.out / ".batch_stage_status" / ("b" * 64) / "candidates.json").exists())

    def test_completed_batch_can_be_replaced_but_not_reactivated(self):
        self.assertEqual(self.run_stage("primary"), 0)
        self.assertEqual(self.run_stage("primary", batch="b" * 64), 0)
        with self.assertRaisesRegex(RuntimeError, "superseded"):
            self.run_stage("primary")
        self.assertEqual(json.loads((self.out / ".batch_stage_owner.json").read_text())["batch_id"], "b" * 64)

    def test_same_batch_cannot_change_original_arguments_or_selected_stages(self):
        self.assertEqual(self.run_stage("primary"), 0)
        for selected, extra in (("primary,report", ()), ("primary", ("--label", "changed"))):
            with self.assertRaisesRegex(ValueError, "different arguments/stages"):
                self.run_stage("primary", selected, extra=extra)

    def test_time_prepares_metadata_and_records_it_as_stage_input(self):
        self.assertEqual(self.run_stage("time"), 0)
        self.assertEqual(self.stub.metadata_calls, 1)
        self.assertIn(self.out / "metadata/samples.tsv", self.stub.calls[0][2])
        self.stub.missing["time"] = ["missing.pgen"]
        self.assertEqual(self.run_stage("time"), 3)
        self.assertEqual(self.stub.metadata_calls, 1)

    def test_unchanged_stage_reuse_remains_success_and_archaic_always_checks_cache(self):
        self.stub.results["primary"] = "reused"
        self.assertEqual(self.run_stage("primary", "primary,archaic"), 0)
        self.assertEqual(self.record("primary")["status"], "reused")
        self.assertEqual(self.run_stage("archaic", "primary,archaic"), 0)
        self.assertTrue(self.stub.calls[-1][1])

    def test_login_guard_runs_before_scientific_import(self):
        with patch.dict(os.environ, {"LSB_JOBID": ""}), patch.object(bp, "load_pipeline") as imported:
            with self.assertRaisesRegex(RuntimeError, "login node"):
                self.run_stage("primary")
            imported.assert_not_called()
        with patch.dict(os.environ, {"ALLELE_BATCH_WORKER": "0"}):
            with self.assertRaisesRegex(RuntimeError, "batch_runner"):
                self.run_stage("primary")

    def test_cli_help_and_import_do_not_import_original_pipeline(self):
        code = "import sys;sys.path.insert(0,sys.argv[1]);import batch_pipeline;assert 'herv_pipeline' not in sys.modules;batch_pipeline.main(['--help'])"
        result = subprocess.run([sys.executable, "-I", "-c", code, str(ROOT / "f")], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--args-json", result.stdout)

    def test_invalid_request_and_planning_flags_are_rejected(self):
        for names in ("primary,primary", "primary,unknown", "primary,", "../primary"):
            with self.assertRaises(ValueError):
                self.run_stage("primary", names)
        with self.assertRaisesRegex(ValueError, "SHA prefix"):
            self.run_stage("primary", batch="not-a-hash")
        with self.assertRaisesRegex(ValueError, "--dry-run/--preflight"):
            self.run_stage("primary", extra=("--preflight",))
        with redirect_stderr(io.StringIO()):
            code = bp.main(["--stage", "primary", "--selected", "primary", "--batch-id", self.batch,
                            "--args-json", '{"not":"an argument array"}'])
        self.assertEqual(code, 2)
        self.assertEqual(self.stub.calls, [])


if __name__ == "__main__":
    unittest.main()
