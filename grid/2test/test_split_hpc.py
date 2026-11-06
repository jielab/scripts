"""Mixed-resource DAGs, frozen handoffs, and CUDA-free publication."""
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from test_hpc_workflow import planner
from test_count_abm import exported_fixture
from count_abm import workflow as count

sys.path.insert(0, str(Path(__file__).parent / "abm"))
from test_prsformer_independent import cache_fixture, evaluated_run, value
import prsformer_resume
from prsformer_stages import freeze_model, frozen_model
from abm_core import prsformer_adapter


class SplitStagesTests(unittest.TestCase):
    def test_gpu_jobs_do_not_inherit_cpu_resources_and_depend_on_cpu_preparation(self):
        args = planner.parser().parse_args([
            "--entry", "prsformer", "--queue", "73x", "--gpu-queue", "2v100-32-sc",
            "--job-cores", "32", "--job-memory-gb", "256", "--threads", "16",
            "--gpu-cores", "8", "--gpu-memory-gb", "160", "--gpu-threads", "2"])
        with patch.dict(os.environ, submit="hpc"):
            jobs = planner.build_jobs(args, [])
        self.assertEqual(len(jobs), 13)
        prepare = jobs[0]
        self.assertEqual(prepare.resources["gpus"], 0)
        self.assertEqual(prepare.resources["memory_gib"], 256)
        for model, report in zip(jobs[1::2], jobs[2::2]):
            self.assertEqual(model.depends, [prepare.key])
            self.assertEqual(report.depends, [model.key])
            self.assertEqual(model.resources["cores"], 8)
            self.assertEqual(model.resources["memory_gib"], 160)
            self.assertEqual(model.env["GRID_THREADS"], "2")
            self.assertEqual(report.resources["gpus"], 0)
            self.assertEqual(report.env["WORKFLOW_DEVICE"], "cpu")
            self.assertEqual(model.env["GRID_PRSFORMER_STAGE_DIR"], report.env["GRID_PRSFORMER_STAGE_DIR"])
            self.assertNotIn("GRID_TRAIT_STEPS", model.env)

    def test_count_cuda_keeps_raw_inputs_and_geometry_on_cpu(self):
        args = planner.parser().parse_args(["--entry", "grid", "--device", "cuda", "--traits", "height,ldl"])
        with patch.dict(os.environ, submit="hpc"):
            jobs = planner.build_jobs(args, [])
        for job in jobs:
            if "gpu_radius_twins" in job.key:
                self.assertEqual(job.resources["gpus"], 1)
                self.assertTrue(job.depends[0].endswith("cpu_geometry"))
            else:
                self.assertEqual(job.resources["gpus"], 0, job.key)
                self.assertEqual(job.env["WORKFLOW_DEVICE"], "cpu")

    def test_geometry_does_not_estimate_radius_and_gpu_stage_cannot_rebuild_it(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            root = Path(tmp)
            manifest, trait = exported_fixture(root / "input")
            argv = ["--input-manifest", str(manifest), "--analysis-root", str(root / "analysis"),
                    "--traits", trait, "--device", "cpu", "--radius-k-ref", "8"]
            with patch("abm_core.support.estimate_radius", side_effect=AssertionError("GPU work in CPU geometry")):
                count.main(["geometry", *argv])
            folder = root / "analysis/count_abm" / trait
            self.assertFalse((folder / "PREPARE_COMPLETE.json").exists())
            before = (folder / "GEOMETRY_COMPLETE.json").read_bytes()
            with patch.object(count, "load_export", side_effect=AssertionError("GPU stage rebuilt CPU inputs")):
                count.main(["radius_twins", *argv])
            self.assertEqual((folder / "GEOMETRY_COMPLETE.json").read_bytes(), before)
            self.assertTrue((folder / "TWINS_COMPLETE.json").is_file())
            donor = folder / "_LOCAL_ONLY/donor.z.npy"
            z = np.load(donor)
            z[0, 0] += 1
            np.save(donor, z)
            with self.assertRaisesRegex(ValueError, "missing or changed"):
                count.main(["radius", *argv])

    def test_prsformer_model_does_not_publish_and_cpu_report_needs_no_cuda(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            root = Path(tmp)
            cache = cache_fixture(root)
            calls = []
            def execute(command, **kwargs):
                mode = value(command, "--prsformer-mode")
                calls.append(mode)
                if mode == "train":
                    evaluated_run(Path(value(command, "--prsformer-run-dir")), cache, "height")
            argv = ["--analysis-root", str(root), "--traits", "height"]
            with patch.dict(os.environ, GRID_PRSFORMER_STAGE_DIR=str(root / "handoff")), \
                 patch("abm_core.cold_start.benchmark_device", return_value="cuda"), \
                 patch("abm_core.cold_start.ensure_inputs"), \
                 patch.object(prsformer_resume, "subprocess", SimpleNamespace(run=execute)):
                prsformer_resume.dispatch(planner.ROOT, [*argv, "--workflow-stage", "model"])
                self.assertEqual(calls, ["train"])
                run = frozen_model(root, "prsformer", "height", cache)
                with patch("abm_core.cold_start.benchmark_device", side_effect=AssertionError("CPU report probed CUDA")):
                    prsformer_resume.dispatch(planner.ROOT, [*argv, "--workflow-stage", "report", "--device", "cpu"])
                self.assertEqual(calls, ["train", "report"])
                (run / "model.pt").write_bytes(b"different model")
                with self.assertRaisesRegex(ValueError, "changed"):
                    frozen_model(root, "prsformer", "height", cache)

    def test_prsformer2_cpu_report_uses_recorded_run_without_cuda_or_role_refitting(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            root = Path(tmp)
            cache = cache_fixture(root)
            run = root / "prsformer2/height/frozen-run"
            local = run / "_LOCAL_ONLY"
            evaluated_run(local, cache, "height")
            (local / "roles.tsv.gz").write_bytes(b"frozen roles")
            (run / "RUN_MANIFEST.json").write_text("{}")
            with patch.dict(os.environ, GRID_PRSFORMER_STAGE_DIR=str(root / "handoff")):
                freeze_model(root, "prsformer2", "height", local, cache)
                with patch.object(prsformer_adapter, "benchmark_device", side_effect=AssertionError("CPU report probed CUDA")), \
                     patch.object(prsformer_adapter, "_run", side_effect=AssertionError("CPU report retrained")):
                    prsformer_adapter.main(["--analysis-root", str(root), "--traits", "height",
                                            "--workflow-stage", "report", "--device", "cpu"])
                self.assertTrue((run / "COMPLETE.json").is_file())
                self.assertEqual(json.loads((run.parent / "CURRENT.json").read_text())["run"], run.name)


if __name__ == "__main__":
    unittest.main()
