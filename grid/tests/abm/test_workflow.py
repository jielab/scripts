import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from abm_core import pipeline, workflow
from abm_core import prsformer_adapter as adapter
import pandas as pd


class WorkflowTests(unittest.TestCase):
    def test_all_modules_complete_before_final(self):
        events = []
        with tempfile.TemporaryDirectory() as tmp, \
             patch("abm_core.cold_start.benchmark_device", return_value="cuda"), \
             patch("abm_core.cold_start.ensure_inputs", side_effect=lambda *a: events.append("prepare")), \
             patch("prsformer_resume.dispatch", side_effect=lambda *a: events.append("prsformer")), \
             patch.object(adapter, "main", side_effect=lambda *a: events.append("prsformer2")), \
             patch.object(pipeline, "main", side_effect=lambda *a: events.append("abm")), \
             patch("final_dispatch.main", side_effect=lambda *a: events.append("final")):
            workflow.main(["--analysis-root", tmp, "--traits", "height,ldl,t2dm"])
        self.assertEqual(events, ["prepare", "prsformer", "prsformer2", "abm", "final"])

    def test_failure_in_any_module_prevents_final(self):
        for failed in ("prepare", "prsformer", "prsformer2", "abm"):
            with self.subTest(failed=failed):
                events = []
                def record(name):
                    def call(*a):
                        events.append(name)
                        if name == failed:
                            raise RuntimeError("injected failure")
                    return call
                with tempfile.TemporaryDirectory() as tmp, \
                     patch("abm_core.cold_start.benchmark_device", return_value="cuda"), \
                     patch("abm_core.cold_start.ensure_inputs", side_effect=record("prepare")), \
                     patch("prsformer_resume.dispatch", side_effect=record("prsformer")), \
                     patch.object(adapter, "main", side_effect=record("prsformer2")), \
                     patch.object(pipeline, "main", side_effect=record("abm")), \
                     patch("final_dispatch.main") as final:
                    with self.assertRaisesRegex(RuntimeError, "injected"):
                        workflow.main(["--analysis-root", tmp])
                    final.assert_not_called()
                self.assertEqual(events[-1], failed)

    def test_module_only_does_not_train_benchmarks_or_publish(self):
        with patch("prsformer_resume.dispatch") as first, patch.object(adapter, "main") as second, \
             patch.object(pipeline, "main") as abm, patch("final_dispatch.main") as final:
            workflow.main(["--module-only", "--traits", "ldl"])
        abm.assert_called_once_with(["--traits", "ldl"])
        first.assert_not_called()
        second.assert_not_called()
        final.assert_not_called()

    def test_runtime_failure_stops_before_input_reconstruction(self):
        with patch("abm_core.cold_start.benchmark_device", side_effect=RuntimeError("NATTEN missing")), \
             patch("abm_core.cold_start.ensure_inputs") as prepare:
            with self.assertRaisesRegex(RuntimeError, "NATTEN"):
                workflow.main([])
        prepare.assert_not_called()

    def test_bare_options_and_explicit_abm_are_standalone(self):
        import twin
        for argv in (["--traits", "height,ldl,t2dm"], ["abm", "--traits", "height,ldl,t2dm"]):
            with patch("automatic_check.ensure_self_check") as check, \
                 patch.object(pipeline, "main") as run:
                twin.main(argv)
            check.assert_called_once()
            run.assert_called_once_with(["--traits", "height,ldl,t2dm"])

    def test_removed_aliases_are_rejected_before_self_check(self):
        import twin
        with patch("automatic_check.ensure_self_check") as check:
            for name in ("abm_GPU", "prsformer-calibrated", "prsformer-v2", "abm-legacy"):
                with self.assertRaisesRegex(ValueError, "Unknown module"):
                    twin.main([name])
        check.assert_not_called()

    def test_module_list_dry_run_and_help_through_shell_entry(self):
        for flag in ("--dry-run", "--help"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp) / "absent analysis"
                result = subprocess.run(
                    ["bash", str(pipeline.ROOT / "4.prsformer.sh"), "all",
                     "--traits", "height,ldl", "--analysis-root", str(root), flag],
                    env=dict(os.environ, GRID_PYTHON=sys.executable, SCRIPT_LOG_DIR=str(Path(tmp) / "logs")),
                    capture_output=True, text=True, check=True)
                self.assertFalse(root.exists())
                output = result.stdout
                if flag == "--dry-run":
                    self.assertLess(output.index("PRSformer plan:"), output.index("PRSformer2 plan:"))
                    self.assertIn("PRSformer plan:", output)
                    self.assertIn("PRSformer2 plan:", output)
                    self.assertEqual(output.count("independent traits in order: height -> ldl"), 2)
                    for module in ("prsformer", "prsformer2"):
                        for trait in ("height", "ldl"):
                            self.assertIn(str(root / module / trait), output)
                else:
                    self.assertIn("usage: ./4.prsformer.sh prsformer ", output)
                    self.assertIn("usage: ./4.prsformer.sh prsformer2 ", output)

    def test_invalid_module_list_fails_before_starting_any_module(self):
        import twin
        with patch("twin.subprocess.run") as run, patch("automatic_check.ensure_self_check") as check:
            for modules in ("prsformer,unknown", "prsformer,", ",prsformer", "prsformer,,prsformer2"):
                with self.subTest(modules=modules), self.assertRaisesRegex(ValueError, "Unknown module"):
                    twin.main([modules, "--traits", "height,ldl"])
        run.assert_not_called()
        check.assert_not_called()

    def test_module_list_preserves_arguments_order_and_failure_status(self):
        import twin
        shared = ["--traits", "height,ldl", "--analysis-root", "/tmp/grid analysis"]
        calls = []
        def run(cmd, **kwargs):
            calls.append(cmd[2:])
            if cmd[2] == "prsformer":
                raise subprocess.CalledProcessError(7, cmd)
        with patch("twin.subprocess.run", side_effect=run):
            with self.assertRaises(subprocess.CalledProcessError) as caught:
                twin.main(["prsformer2,prsformer,final", *shared])
        self.assertEqual(caught.exception.returncode, 7)
        self.assertEqual(calls, [["prsformer2", *shared], ["prsformer", *shared]])

    def test_workflow_dry_run_is_read_only_even_with_empty_root(self):
        text = io.StringIO()
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(text):
            root = Path(tmp) / "absent"
            workflow.main(["--analysis-root", str(root), "--config", str(pipeline.ROOT / "config"),
                           "--dry-run"])
            self.assertFalse(root.exists())
        plan = json.loads(text.getvalue())
        self.assertEqual(plan["modules"], ["prsformer", "prsformer2", "abm", "final"])
        self.assertEqual(plan["device"], "auto")
        self.assertEqual(plan["predictor"], "blocknet")
        self.assertEqual(plan["representation"], "csx_ldblocks_fixed")

    def test_partial_work_requires_module_only(self):
        with self.assertRaises(SystemExit):
            workflow.main(["--through", "prepare"])


class PRSformer2ResumeTests(unittest.TestCase):
    def test_cleanup_runs_after_success_and_skip_but_not_check_or_failure(self):
        cases = [(keep, check, failed) for keep in (None, "0", "1")
                 for check, failed in ((False, False), (True, False), (False, True))]
        for keep, check, failed in cases:
            with self.subTest(keep=keep, check=check, failed=failed), \
                 patch.dict(os.environ), \
                 patch.object(adapter, "benchmark_device", return_value="cuda"), \
                 patch.object(adapter, "_run", side_effect=RuntimeError("interrupted") if failed else None) as run, \
                 patch("dosage_storage.release_dosages") as cleanup:
                # Exercise both policies independently of the launching shell.
                if keep is None:
                    os.environ.pop("GRID_KEEP_SHARED_INPUTS", None)
                else:
                    os.environ["GRID_KEEP_SHARED_INPUTS"] = keep
                argv = ["--check"] if check else []
                if failed:
                    with self.assertRaisesRegex(RuntimeError, "interrupted"):
                        adapter.main(argv)
                else:
                    adapter.main(argv)
                run.assert_called_once()
                if keep == "1" or check or failed:
                    cleanup.assert_not_called()
                else:
                    cleanup.assert_called_once()
                    self.assertTrue(cleanup.call_args.kwargs["apply"])
                    self.assertTrue(cleanup.call_args.kwargs["best_effort"])

    def test_roles_preserve_bytes_and_reject_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "roles.tsv.gz"
            roles = pd.DataFrame(dict(id=["001", "002"], group_id=["01", "02"],
                                      role=["build", "test"], parent_role=["build", "test"]))
            adapter.save_roles(path, roles)
            before, stamp = path.read_bytes(), path.stat().st_mtime_ns
            adapter.save_roles(path, roles)
            self.assertEqual(before, path.read_bytes())
            self.assertEqual(stamp, path.stat().st_mtime_ns)
            roles.loc[0, "role"] = "test"
            with self.assertRaises(AssertionError):
                adapter.save_roles(path, roles)

    def test_resume_and_prediction_require_complete_calibration(self):
        import torch
        from abm_core.model import file_sha256
        with tempfile.TemporaryDirectory() as tmp:
            local = Path(tmp)
            self.assertEqual(adapter.action_for(local), "train")
            model = local / "model.pt"
            torch.save({"preprocessing": {}}, model)
            with self.assertRaisesRegex(ValueError, "Incomplete"):
                adapter.action_for(local)
            (local / "training.pt").write_bytes(b"fixture")
            self.assertEqual(adapter.action_for(local), "resume")
            receipt = local / "TRAINED.json"
            receipt.write_text(json.dumps({"model_sha256": file_sha256(model)}))
            self.assertEqual(adapter.action_for(local), "resume")
            torch.save({"preprocessing": {"prediction_calibration": {"height": [1, 0]}}}, model)
            receipt.write_text(json.dumps({"model_sha256": file_sha256(model)}))
            self.assertEqual(adapter.action_for(local), "predict")
            model.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "changed"):
                adapter.action_for(local)

    def test_dry_run_does_not_create_run_or_require_cache(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(adapter, "ensure_inputs") as prepare:
            root = Path(tmp) / "absent"
            adapter.main(["--analysis-root", str(root), "--dry-run"])
            self.assertFalse(root.exists())
            prepare.assert_not_called()
