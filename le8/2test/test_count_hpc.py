"""Executable synthetic correctness checks; no UKB data, R, GPU or LSF needed."""
from __future__ import annotations

from dataclasses import asdict
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
sys.path.append(str(ROOT.parent / "0f"))
import abm
import c1_count
from c1_models import context, fit_lightgbm
from c1_protocol import build_abm_config
from c1_evaluation import evaluate_methods
import hpc_workflow
import workflow_jobs


def synthetic():
    rng = np.random.default_rng(412)
    n, nt = 600, 120
    x = rng.normal(size=(n + nt, 7))
    y = (rng.uniform(size=n + nt) < 1 / (1 + np.exp(-x[:, 1]))).astype(int)
    time = np.where(y == 1, rng.uniform(1, 9, n + nt), rng.uniform(11, 15, n + nt))
    # Known early censoring must never become a donor control.
    censored = np.arange(0, n, 17)
    y[censored] = 0
    time[censored] = 3
    ids = np.array([f"p{i:05}" for i in range(n + nt)])
    groups = np.array([f"g{i // 2:05}" for i in range(n + nt)])
    features = np.array(["age", *[f"protein{i}" for i in range(6)]])
    d = dict(X=x[:n], eid=ids[:n], group=groups[:n], feature_names=features,
             role=np.repeat(abm.ROLES, n // 5), time=time[:n], event=y[:n])
    tx = x[n:].copy()
    tx[0, 1:] = 1e5  # no geometric donors
    tx[1, 1:] = np.nan  # technical-QC failure
    t = dict(X=tx, eid=ids[n:], group=groups[n:], feature_names=features)
    truth = dict(eid=ids[n:], time=time[n:], event=y[n:])
    cfg = build_abm_config("abm_count", dict(min_events=2, audit_min_cases=2, audit_min_controls=2,
        bootstrap=20, k=12, min_twins=4, min_family_ess=3, dimensions=3, prior_strength=8),
        radius_k_ref=12, cores=1, parallel_workers=1, always_include_indices=(0,))
    return d, t, truth, cfg


class TransparentCountTests(unittest.TestCase):
    def test_exact_donor_evidence_abstention_and_common_evaluation(self):
        d, t, truth, cfg = synthetic()
        c = context(d, t, cfg)
        with tempfile.TemporaryDirectory() as td, \
             patch("c1_models._fit_assessment", side_effect=AssertionError("regression head forbidden")), \
             patch.object(abm.Calibrator, "fit", side_effect=AssertionError("risk calibration forbidden")):
            out = Path(td) / "fit"
            saved = c1_count.fit_count(c, out)
            p = saved["test"]
            self.assertTrue(np.isnan(p["predicted_error"]).all())
            self.assertFalse(p["supported"][0])
            self.assertFalse(p["supported"][1])
            self.assertEqual(p["diagnostics"]["abstention_code"][0], 2)
            self.assertEqual(p["diagnostics"]["abstention_code"][1], 1)
            self.assertTrue(np.isnan(p["diagnostics"]["risk_for_decision"][:2]).all())
            self.assertGreater(p["supported"].sum(), 70)
            trace = pd.read_csv(out / "donor_evidence_LOCAL_ONLY/test_donors.tsv.gz", sep="\t", dtype={"query_eid": str})
            self.assertTrue((trace.kernel_weight == 1).all())
            self.assertFalse((trace.query_eid == trace.donor_eid).any())
            self.assertFalse((trace.query_group == trace.donor_group).any())
            self.assertFalse(set(d["eid"][np.arange(0, len(d["eid"]), 17)]) & set(trace.donor_eid))
            queries = pd.read_csv(out / "donor_evidence_LOCAL_ONLY/test_queries.tsv.gz", sep="\t").set_index("eid")
            for query, rows in trace.groupby("query_eid"):
                q = queries.loc[query]
                self.assertAlmostEqual(rows.normalized_weight.sum(), 1., places=12)
                local = float(np.sum(rows.normalized_weight * rows.donor_outcome))
                family_weights = rows.groupby("donor_group").normalized_weight.sum()
                ess = 1 / np.sum(family_weights ** 2)
                self.assertAlmostEqual(q.family_ess, ess, places=10)
                self.assertAlmostEqual(q.local_event_rate, local, places=12)
                self.assertEqual(q.matched_case_count, rows.donor_outcome.sum())
                expected = ((ess * local + cfg.prior_strength * q.prior_event_rate) / (ess + cfg.prior_strength)
                            if q.supported else q.prior_event_rate)
                self.assertAlmostEqual(q.probability, expected, places=12)
            model = joblib.load(out / "model_LOCAL_ONLY.joblib")
            found = []
            model["bank"].predict(c.z[:1], c.ids[:1], c.groups[:1], np.array([model["prior"]]),
                                   evidence=lambda i, rows: found.extend(rows["donor_group"]))
            self.assertNotIn(c.groups[0], found)
            tables = evaluate_methods({"abm_count": p}, truth, c.km, bootstrap=20, primary_coverage=.5)
            all_rows = tables["performance"].query("subset == 'all_common'")
            self.assertEqual(int(all_rows.iloc[0].evaluated_N), len(t["eid"]) - 1)
            perf = tables["performance"].set_index("subset")
            self.assertEqual(int(perf.loc["selected::abm_count::0.5", "evaluated_N"] +
                                 perf.loc["not_selected::abm_count::0.5", "evaluated_N"]), len(t["eid"]) - 1)
            self.assertFalse(any("eid" in frame.columns for frame in tables.values()))
            manifest = json.loads((out / "MODEL_FROZEN.json").read_text())
            self.assertFalse(manifest["query_regression_head"])
            self.assertFalse(manifest["test_outcomes_used"])
            replay = c1_count.predict_frozen_count(out, t["X"], t["eid"], t["group"], t["feature_names"], Path(td) / "new_queries")
            np.testing.assert_array_equal(replay["probability"], p["probability"])
            np.testing.assert_array_equal(replay["supported"], p["supported"])

    @unittest.skipUnless(importlib.util.find_spec("lightgbm"), "Install LightGBM for the conventional-baseline integration check")
    def test_real_lightgbm_and_count_share_evaluation_roster(self):
        from dataclasses import replace
        from threadpoolctl import threadpool_limits
        d, t, truth, cfg = synthetic()
        cfg = replace(cfg, tree_iterations=8, gate_iterations=3)
        c = context(d, t, cfg)
        with tempfile.TemporaryDirectory() as td, threadpool_limits(1):
            count = c1_count.fit_count(c, Path(td) / "count")["test"]
            conventional = fit_lightgbm(c, Path(td) / "lightgbm")["test"]
            tables = evaluate_methods({"lightgbm": conventional, "abm_count": count}, truth, c.km,
                                      bootstrap=20, primary_coverage=.5)
            performance = tables["performance"]
            for subset, rows in performance.groupby("subset"):
                self.assertEqual(rows.evaluated_N.nunique(), 1)
                self.assertEqual(set(rows.model), {"lightgbm", "abm_count"})
            self.assertTrue(np.isnan(count["predicted_error"]).all())
            self.assertTrue(np.isfinite(conventional["predicted_error"]).all())

    def test_nonbuild_outcomes_cannot_change_query_count_prediction(self):
        d, t, truth, cfg = synthetic()
        changed = {k: v.copy() for k, v in d.items()}
        ix = changed["role"] != "build"
        changed["event"][ix] = 1 - changed["event"][ix]
        changed["time"][ix] = np.where(changed["event"][ix] == 1, 5., 13.)
        with tempfile.TemporaryDirectory() as td:
            a = c1_count.fit_count(context(d, t, cfg), Path(td) / "a")["test"]
            b = c1_count.fit_count(context(changed, t, cfg), Path(td) / "b")["test"]
            np.testing.assert_array_equal(a["probability"], b["probability"])
            np.testing.assert_array_equal(a["supported"], b["supported"])
            for q in a["masks"]:
                np.testing.assert_array_equal(a["masks"][q], b["masks"][q])

    def test_real_restartable_models_and_evaluation_cli(self):
        d, t, truth, cfg = synthetic()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            source = root / "prepared"
            source.mkdir()
            for name, data in (("development", d), ("test_X", t), ("test_truth", truth)):
                np.savez_compressed(source / f"{name}_LOCAL_ONLY.npz", **data)
            (source / "host_contract.json").write_text(json.dumps(dict(trait="toy", layer="prot",
                clinical_feature_names=["age"], independent_family_support_verified=True)))
            config = root / "abm.json"
            config.write_text(json.dumps(asdict(cfg)))
            command = [sys.executable, str(ROOT / "f/c1_correlate.py"), "--Y", "toy", "--biom", "prot",
                "--methods", "abm_count", "--prepared-source", str(source), "--analysis-root", str(root / "results"),
                "--abm-config", str(config), "--association", "skip", "--min-cases-with-biom", "2",
                "--bootstrap", "20", "--cores", "1"]
            env = dict(os.environ, WORKFLOW_WORKER="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
            for stage in ("prepare", "models", "evaluation"):
                run = subprocess.run(command + ["--stage", stage], env=env, capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            files = list((root / "results").rglob("test_predictions_LOCAL_ONLY.joblib"))
            self.assertEqual(len(files), 1)
            frozen = files[0].read_bytes()
            run = subprocess.run(command + ["--stage", "evaluation", "--replace", "TRUE"],
                                 env=env, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertEqual(frozen, files[0].read_bytes())
            evaluation = root / "results/toy/prot/c1_correlate/evaluation/FINAL"
            self.assertTrue((evaluation / "performance.tsv").is_file())
            self.assertTrue((evaluation / "C1_unified_evaluation.xlsx").is_file())


class WorkflowTests(unittest.TestCase):
    def test_lsf_phase_requires_allocation_and_respects_budget(self):
        spec = importlib.util.spec_from_file_location("le8_test_resources", ROOT / "f/0.resources.py")
        resources = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(resources)
        command = ["--cores", "2", "--memory-gib", "4", "--", sys.executable, "-c",
                   "import os; assert os.environ['LE8_PHASE_CORES']=='2'"]
        env = dict(WORKFLOW_RESOURCE_MANAGER="lsf", LSB_JOBID="synthetic-test-only",
                   WORKFLOW_CORES="4", WORKFLOW_MEMORY_GIB="8", WORKFLOW_GPUS="0")
        with patch.dict(os.environ, env), patch.object(resources, "SharedBudget", side_effect=AssertionError("No cross-job pool")):
            self.assertEqual(resources.phase_main(command), 0)
            with self.assertRaises(ValueError):
                resources.phase_main(["--cores", "5", "--memory-gib", "4", "--", "true"])
            with self.assertRaises(ValueError):
                resources.phase_main(["--cores", "2", "--memory-gib", "9", "--", "true"])
            with self.assertRaises(ValueError):
                resources.phase_main(["--cores", "2", "--memory-gib", "4", "--gpu", "--", "true"])
            with patch.dict(os.environ, {"LSB_JOBID": ""}):
                with self.assertRaises(RuntimeError):
                    resources.phase_main(command)
        with patch.dict(os.environ, dict(WORKFLOW_RESOURCE_MANAGER="local", LE8_PHASE_ADMITTED="1",
                  LE8_PHASE_CORES="4", LE8_PHASE_MEMORY_GIB="8")), \
             patch.object(resources, "SharedBudget", side_effect=AssertionError("Nested admission must not deadlock")):
            self.assertEqual(resources.phase_main(["--cores", "2", "--memory-gib", "4", "--", "true"]), 0)

    def test_public_hpc_worker_refuses_login_and_keeps_job_tmpdir(self):
        env = dict(os.environ, submit="hpc", SUBMIT="hpc", WORKFLOW_WORKER="1", LSB_JOBID="",
                   WORKFLOW_PYTHON=sys.executable)
        run = subprocess.run(["bash", str(ROOT / "le8.sh"), "--help"], env=env, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("LSF allocation", run.stderr)
        with tempfile.TemporaryDirectory() as td:
            env.update(submit="local", SUBMIT="local", TMPDIR=td, LE8_ROOT=str(ROOT),
                       ANALYSIS_ROOT=td, DIR_WORK="/test/work", DIR_DATA="/test/data", MEMORY_DEFAULT_GIB="16")
            command = ["bash", "-c", 'source "$1"; "$2" -c "import os; assert os.environ[\'TMPDIR\']==os.environ[\'EXPECTED_TMP\']; assert os.environ[\'LE8_RESOURCE_POOL\'].startswith(os.environ[\'EXPECTED_TMP\'])"',
                       "le8-env-test", str(ROOT / "f/runtime_env.sh"), sys.executable]
            env["EXPECTED_TMP"] = td
            run = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)

    def test_all_modules_dag_and_mock_lsf_dependencies(self):
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, {"submit": "hpc", "RUN_GPU_COLOC": "FALSE"}):
            a, extra = hpc_workflow.parse_cli(["c1_correlate,c2_cause,c3_coloc,c4_connect,c5_cellulation,pgs_focus,final",
                "--Y", "toy,other", "--biom", "prot,met", "--analysis-root", td, "--eigh-backend", "cpu"])
            jobs = hpc_workflow.build_jobs(a)
            by_key = {j.key: j for j in jobs}
            self.assertIn("toy.prot.c1.prepare", by_key)
            self.assertIn("other.met.c4_explain", by_key)
            self.assertIn("other.met.c5_cellulation", by_key)
            self.assertIn("toy.prot.pgs_focus", by_key)
            self.assertEqual(by_key["toy.prot.c1.models"].depends, ["toy.prot.c1.association"])
            self.assertFalse(by_key["other.met.c1.prepare"].depends)
            self.assertFalse(any(j.resources["gpus"] for j in jobs))
            spec = importlib.util.spec_from_file_location("le8_plan_test_common", ROOT / "f/0.common.py")
            common = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(common)
            # Validate the ACTUAL worker parsers as well as scheduler records.
            # Scientific method flags must never be passed to C2/C3/Final.
            for job in jobs:
                if job.argv[0] == "bash":
                    parsed = common.pipeline_parse(job.argv[2:])
                    self.assertEqual(parsed[2], [job.argv[2]])
                    if job.argv[2] != "c1_correlate":
                        self.assertNotIn("--methods", job.argv)
            path = workflow_jobs.write_plan(jobs, Path(td) / "plans", mode="hpc", name="le8test")
            submitted = []

            def bsub(argv, **kwargs):
                self.assertEqual(argv[0], "bsub")
                text = kwargs["stdin"].read()
                self.assertIn("LSB_JOBID", text)
                self.assertIn("#BSUB -M", text)
                submitted.append((argv, text))
                return subprocess.CompletedProcess(argv, 0, f"Job <{1000 + len(submitted)}> is submitted", "")

            with patch.object(workflow_jobs.subprocess, "run", side_effect=bsub):
                result = workflow_jobs.run_plan(path, submit="hpc", jobs=128, wait=False)
                workflow_jobs.run_plan(path, submit="hpc", jobs=128, wait=False)
            self.assertEqual(len(submitted), len(jobs))  # resume does not duplicate jobs
            self.assertEqual(result["status"], "SUBMITTED")
            self.assertTrue(any("-w" in call[0] for call in submitted))
            a.device = "cuda"
            gpu_jobs = [j for j in hpc_workflow.build_jobs(a) if j.resources["gpus"]]
            self.assertEqual(len(gpu_jobs), 4)
            self.assertTrue(all(j.key.endswith("c3_coloc") for j in gpu_jobs))
            self.assertTrue(all(j.resources["cores"] >= 20 for j in gpu_jobs))
            a.fit_joint = True
            fitted_jobs = hpc_workflow.build_jobs(a)
            self.assertEqual([j.key for j in fitted_jobs if "--fit-joint" in j.argv], ["aggregate.final"])

    def test_default_modules_preserve_unknown_option_values(self):
        a, extra = hpc_workflow.parse_cli(["--ukb-phe", "/some/input", "--prepared-source", "/some/prepared",
                                          "--Y", "toy", "--biom", "prot"])
        self.assertEqual(a.modules, hpc_workflow.DEFAULT)
        self.assertEqual(extra, ["--ukb-phe", "/some/input", "--prepared-source", "/some/prepared"])
        jobs = hpc_workflow.build_jobs(a, extra)
        self.assertTrue(all("/some/prepared" in job.argv for job in jobs if ".c1." in job.key))

    def test_explicit_coloc_switch_and_gpu_allocation_agree(self):
        for switch, device, inherited, expected in (
                ("TRUE", "auto", "FALSE", 1),
                ("FALSE", "cuda", "TRUE", 0),
                ("yes", "cpu", "FALSE", 1)):
            with patch.dict(os.environ, {"submit": "hpc", "RUN_GPU_COLOC": inherited}):
                args, extra = hpc_workflow.parse_cli(["c3_coloc", "--Y", "toy", "--biom", "prot",
                    "--device", device, "--run-gpu-coloc", switch])
                self.assertNotIn("--run-gpu-coloc", extra)
                job = hpc_workflow.build_jobs(args, extra)[0]
                self.assertEqual(job.resources["gpus"], expected)
                self.assertEqual(job.env["RUN_GPU_COLOC"], "TRUE" if expected else "FALSE")

    def test_cpu_device_overrides_site_coloc_default_and_deep_auto_is_scheduled(self):
        with patch.dict(os.environ, {"submit": "hpc", "RUN_GPU_COLOC": "TRUE"}):
            args, extra = hpc_workflow.parse_cli(["c3_coloc", "--Y", "toy", "--biom", "prot", "--device", "cpu"])
            job = hpc_workflow.build_jobs(args, extra)[0]
            self.assertEqual(job.resources["gpus"], 0)
            self.assertEqual(job.env["RUN_GPU_COLOC"], "FALSE")
            args, extra = hpc_workflow.parse_cli(["c1_correlate", "--Y", "toy", "--biom", "prot", "--methods", "abm_deep"])
            jobs = hpc_workflow.build_jobs(args, extra)
            self.assertEqual([j.key for j in jobs if j.resources["gpus"]], ["toy.prot.c1.models"])

    def test_native_hpc_runtime_auto_coloc_generates_matching_gpu_plan(self):
        with tempfile.TemporaryDirectory() as td:
            env = {k:v for k,v in os.environ.items() if k not in {"RUN_GPU_COLOC", "WORKFLOW_WORKER"}}
            run = subprocess.run(["bash", str(ROOT / "le8.sh"), "c3_coloc", "--submit", "hpc",
                "--plan-only", "--Y", "toy", "--biom", "prot", "--analysis-root", td],
                env=env, capture_output=True, text=True, timeout=30)
            self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
            plans = list(Path(td).glob("logs/workflows/*/plan.json"))
            self.assertEqual(len(plans), 1)
            job = json.loads(plans[0].read_text())["jobs"][0]
            self.assertEqual(job["resources"]["gpus"], 1)
            self.assertEqual(job["env"]["RUN_GPU_COLOC"], "TRUE")

    def test_eigh_auto_on_hpc_and_explicit_cpu_choices(self):
        with patch.dict(os.environ, {"submit": "hpc", "MRLINK2_EIGH_BACKEND": "auto"}):
            for options, expected in (([], 1), (["--device", "cpu"], 0),
                                      (["--eigh-backend", "cpu"], 0),
                                      (["--device", "cpu", "--eigh-backend", "cuda"], 1)):
                args, extra = hpc_workflow.parse_cli(["c2_cause", "--Y", "toy", "--biom", "prot", *options])
                job = hpc_workflow.build_jobs(args, extra)[0]
                self.assertEqual(job.resources["gpus"], expected)
                self.assertEqual(job.env["MRLINK2_EIGH_BACKEND"], "cuda" if expected else "cpu")


if __name__ == "__main__":
    unittest.main()
