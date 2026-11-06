"""CPU/GPU scheduling and frozen-artifact handoff; no cluster/GPU is required."""
from contextlib import nullcontext, redirect_stdout, redirect_stderr
from dataclasses import asdict
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
sys.path.append(str(ROOT.parent / "0f"))
import c1_correlate
import c1_models
import hpc_workflow
import scope_workflow
import workflow_jobs
from test_count_hpc import synthetic


def steps(job):
    return json.loads(job.argv[job.argv.index("--steps-json") + 1])


class SplitPlanTests(unittest.TestCase):
    def test_three_allocations_per_scope_and_actual_lsf_dependencies(self):
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, submit="hpc"):
            args, extra = hpc_workflow.parse_cli(["c1_correlate", "--Y", "ad,amr,asthma",
                "--biom", "prot,met", "--analysis-root", td, "--queue", "33",
                "--gpu-queue", "2v100-32-e5", "--cores", "20", "--gpu-cores", "12",
                "--memory-limit-gb", "128", "--gpu-memory-limit-gb", "96", "--device", "cuda"])
            jobs = hpc_workflow.build_jobs(args, extra)
            self.assertEqual(len(jobs), 18)
            self.assertEqual(sum(j.resources["gpus"] for j in jobs), 6)
            by_key = {j.key: j for j in jobs}
            for scope in (f"{t}.{b}" for t in ("ad", "amr", "asthma") for b in ("prot", "met")):
                pre, gpu, post = [by_key[scope + "." + phase] for phase in
                                 ("cpu_prepare_models", "gpu_abm_deep", "cpu_evaluation")]
                self.assertEqual(pre.depends, [])
                self.assertEqual(gpu.depends, [pre.key])
                self.assertEqual(post.depends, [gpu.key])
                self.assertEqual([j.resources["cores"] for j in (pre, gpu, post)], [20, 12, 20])
                self.assertEqual([j.resources["memory_gib"] for j in (pre, gpu, post)], [128, 96, 128])
                self.assertEqual([j.resources["gpus"] for j in (pre, gpu, post)], [0, 1, 0])
                self.assertEqual([s["key"].split(".c1.")[-1] for s in steps(pre)],
                                 ["prepare", "association", "models_cpu"])
                self.assertEqual(len(steps(gpu)), 1)
                gpu_argv = steps(gpu)[0]["argv"]
                self.assertEqual(gpu_argv[gpu_argv.index("--methods") + 1], "abm_deep")
                cpu_argv = steps(pre)[-1]["argv"]
                self.assertEqual(cpu_argv[cpu_argv.index("--methods") + 1], "glmnet,lightgbm,abm_count")
                self.assertEqual(steps(post)[0]["argv"][-2:], ["--device", "cpu"])
                for job in (pre, gpu, post):
                    self.assertIn("--fail-fast", job.argv)
                manifests = {s["env"]["LE8_C1_BATCH_DIR"] for job in (pre, gpu, post)
                             for s in steps(job) if s["key"].split(".c1.")[-1] not in ("reliability", "summary")}
                self.assertEqual(len(manifests), 1)
            path = workflow_jobs.write_plan(jobs, Path(td) / "plans", mode="hpc", name="le8.c1_correlate")
            submitted = []

            def bsub(argv, **kwargs):
                text = kwargs["stdin"].read()
                key = next(line.removeprefix("#BSUB -J le8.c1_correlate.")
                           for line in text.splitlines() if line.startswith("#BSUB -J"))
                submitted.append((key, argv, text))
                return subprocess.CompletedProcess(argv, 0, f"Job <{1000 + len(submitted)}> is submitted", "")

            with patch.object(workflow_jobs, "_submit_bsub", side_effect=bsub):
                result = workflow_jobs.run_plan(path, submit="hpc", jobs=6, wait=False)
                workflow_jobs.run_plan(path, submit="hpc", jobs=6, wait=False)
            self.assertEqual(len(submitted), 18)
            for key, argv, text in submitted:
                parents = by_key[key].depends
                if parents:
                    condition = argv[argv.index("-w") + 1]
                    for parent in parents:
                        self.assertIn(f"done({result['jobs'][parent]['job_id']})", condition)
                else:
                    self.assertNotIn("-w", argv)
                self.assertEqual('#BSUB -gpu' in text, bool(by_key[key].resources["gpus"]))

    def test_phase_failure_stops_children_and_worker_checks_parent_receipt(self):
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, submit="local"):
            root = Path(td)
            output = root / "should-not-exist"
            run = [dict(key="prepare", argv=[sys.executable, "-c", "raise SystemExit(7)"]),
                   dict(key="models", argv=[sys.executable, "-c",
                        "from pathlib import Path; import sys; Path(sys.argv[1]).touch()", str(output)])]
            self.assertEqual(scope_workflow.run_steps(root, "toy.prot", run,
                phase="cpu_prepare_models", fail_fast=True), 1)
            self.assertFalse(output.exists())
            record = json.loads((root / "logs/scopes/toy.prot.cpu_prepare_models.json").read_text())
            self.assertEqual(len(record["steps"]), 1)
            jobs = [workflow_jobs.Job("cpu", [sys.executable, "-c", "pass"]),
                    workflow_jobs.Job("gpu", [sys.executable, "-c", "pass"], ["cpu"])]
            plan = workflow_jobs.write_plan(jobs, root / "plan", mode="local")
            workflow_jobs.atomic_json(plan.parent / "status/cpu.json", dict(status="EXIT"))
            with self.assertRaisesRegex(RuntimeError, "Dependency did not complete successfully"):
                workflow_jobs.worker(plan, "gpu")


class FrozenStageTests(unittest.TestCase):
    def test_cpu_evaluator_uses_exact_gpu_predictions_without_cuda_or_refit(self):
        d, x, truth, cfg = synthetic()
        with tempfile.TemporaryDirectory() as td, threadpool_limits(1), \
             patch.dict(os.environ, WORKFLOW_WORKER="1"), redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            root = Path(td)
            source = root / "inputs"
            source.mkdir()
            for name, values in (("development", d), ("test_X", x), ("test_truth", truth)):
                np.savez_compressed(source / f"{name}_LOCAL_ONLY.npz", **values)
            (source / "host_contract.json").write_text(json.dumps(dict(trait="toy", layer="prot",
                clinical_feature_names=["age"], independent_family_support_verified=True)))
            config = root / "abm.json"
            config.write_text(json.dumps(asdict(cfg)))
            batch = root / "batch"
            command = ["--Y", "toy", "--biom", "prot", "--prepared-source", str(source),
                "--analysis-root", str(root / "results"), "--abm-config", str(config),
                "--association", "skip", "--bootstrap", "20", "--min-cases-with-biom", "2", "--cores", "1"]
            env = dict(LE8_C1_BATCH_DIR=str(batch), LE8_C1_MODEL_PARTS="cpu,gpu")

            def gpu_fit(c, out, network_config, device):
                self.assertEqual(device, "cuda")
                probability = np.full(len(c.ids), .4)
                diagnostics = dict(technical_qc=c.qc, supported=c.qc)
                return c1_models.finish_projection(c, probability, diagnostics, out,
                    fitted={}, metadata=dict(method="abm_deep", actual_device="cuda"))

            with patch.dict(os.environ, env):
                self.assertEqual(c1_correlate.main(command + ["--stage", "prepare", "--device", "cpu",
                    "--methods", "abm_count,abm_deep"]), 0)
                with patch.dict(os.environ, LE8_C1_MODEL_PART="cpu"):
                    self.assertEqual(c1_correlate.main(command + ["--stage", "models", "--device", "cpu",
                        "--methods", "abm_count"]), 0)
                resources = c1_correlate.load_module("le8_c1_resources", ROOT / "f/0.resources.py")
                with patch.dict(os.environ, LE8_C1_MODEL_PART="gpu"), \
                     patch("c1_protocol._torch_cuda_available", return_value=True), \
                     patch.object(resources, "gpu_lease", return_value=nullcontext()), \
                     patch.object(c1_correlate, "fit_deep", side_effect=gpu_fit) as fit:
                    self.assertEqual(c1_correlate.main(command + ["--stage", "models", "--device", "cuda",
                        "--methods", "abm_deep"]), 0)
                    fit.assert_called_once()
                receipt = json.loads((batch / "models.gpu.json").read_text())
                artifact = Path(receipt["methods"]["abm_deep"]["path"]) / "test_predictions_LOCAL_ONLY.joblib"
                original = artifact.read_bytes()
                evaluate = command + ["--stage", "evaluation", "--device", "cpu",
                                      "--methods", "abm_count,abm_deep"]
                with patch("c1_protocol._torch_cuda_available", side_effect=AssertionError("CPU evaluation touched CUDA")), \
                     patch.object(c1_correlate, "fit_deep", side_effect=AssertionError("evaluation refitted deep")), \
                     patch("c1_count.fit_count", side_effect=AssertionError("evaluation refitted count")):
                    self.assertEqual(c1_correlate.main(evaluate), 0)
                self.assertEqual(artifact.read_bytes(), original)
                scope = root / "results/toy/prot/c1_correlate"
                final = json.loads((scope / "C1_FINAL_INDEX.json").read_text())
                self.assertEqual(set(final["methods"]), {"abm_count", "abm_deep"})
                self.assertEqual(final["methods"]["abm_deep"]["device"]["actual"], "cuda")
                self.assertTrue((scope / "evaluation/FINAL/performance.tsv").is_file())
                # Missing handoff, changed artifacts, or changed input bytes must
                # block evaluation even when old model directories still exist.
                marker = batch / "models.gpu.json"
                saved = marker.read_bytes()
                marker.unlink()
                with patch.object(c1_correlate, "evaluate_all") as evaluator:
                    self.assertEqual(c1_correlate.main(evaluate), 1)
                    evaluator.assert_not_called()
                marker.write_bytes(saved)
                artifact.write_bytes(original + b"changed")
                with patch.object(c1_correlate, "evaluate_all") as evaluator:
                    self.assertEqual(c1_correlate.main(evaluate), 1)
                    evaluator.assert_not_called()
                artifact.write_bytes(original)
                with (source / "test_X_LOCAL_ONLY.npz").open("ab") as stream:
                    stream.write(b"changed")
                with patch.object(c1_correlate, "evaluate_all") as evaluator:
                    self.assertEqual(c1_correlate.main(evaluate), 1)
                    evaluator.assert_not_called()


if __name__ == "__main__":
    unittest.main()
