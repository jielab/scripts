"""Regression coverage for one allocation per scope and independent scarce outcomes."""
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

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
from test_count_hpc import synthetic
import abm
import c1_correlate
import c1_report
from c1_evaluation import locate_method, ReliabilityHead
from c1_models import context, InsufficientDevelopmentSupport
import hpc_workflow
import scope_workflow


def bundle(root, d, x, truth, layer="prot"):
    source = root / (layer + "_inputs")
    source.mkdir()
    for name, data in (("development", d), ("test_X", x), ("test_truth", truth)):
        np.savez_compressed(source / f"{name}_LOCAL_ONLY.npz", **data)
    contract = dict(trait="toy", layer=layer, clinical_feature_names=["age"],
                    independent_family_support_verified=True, incident_test_ids_hash=abm.fingerprint(x["eid"]))
    (source / "host_contract.json").write_text(json.dumps(contract))
    return source, contract, None


class SparseC1Tests(unittest.TestCase):
    def arguments(self, root, cfg, source, methods="glmnet,abm_count"):
        config = root / "abm.json"
        config.write_text(json.dumps(asdict(cfg)))
        return ["--Y", "toy", "--biom", "prot", "--methods", methods,
                "--prepared-source", str(source), "--analysis-root", str(root / "results"),
                "--abm-config", str(config), "--min-cases-with-biom", "2", "--bootstrap", "20", "--cores", "1"]

    def test_sparse_tuning_skips_glmnet_but_fits_and_evaluates_count(self):
        d, x, truth, cfg = synthetic()
        poor = d["role"] == "tune_model"
        d["event"][poor] = 0
        d["time"][poor] = 14
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, WORKFLOW_WORKER="1"), \
             patch.object(c1_correlate, "native_association") as association, threadpool_limits(1):
            root = Path(td)
            source, _, _ = bundle(root, d, x, truth)
            # A saved whole-scope exclusion from the old code must not block this run.
            exclusion = root / "results/logs/exclusions/toy.prot.json"
            exclusion.parent.mkdir(parents=True)
            exclusion.write_text(json.dumps(dict(status=c1_correlate.INSUFFICIENT_SCOPE_SUPPORT)))
            with patch.object(c1_correlate, "fit_glmnet", side_effect=AssertionError("unsupported glmnet was fitted")):
                rc = c1_correlate.main(self.arguments(root, cfg, source))
            self.assertEqual(rc, 0)
            association.assert_called_once()
            scope = root / "results/toy/prot/c1_correlate"
            index = json.loads((scope / "C1_FINAL_INDEX.json").read_text())
            self.assertEqual(set(index["methods"]), {"abm_count"})
            self.assertEqual(set(index["unavailable_methods"]), {"glmnet"})
            self.assertTrue((scope / "evaluation/FINAL/performance.tsv").is_file())
            self.assertEqual(json.loads(exclusion.read_text())["status"], "superseded_by_per_analysis_support")

    def test_audit_cases_do_not_block_model_fitting(self):
        d, x, _, cfg = synthetic()
        d["event"][d["role"] == "calibration_audit"] = 0
        d["time"][d["role"] == "calibration_audit"] = 14
        c = context(d, x, cfg, required_roles=set())
        c1_correlate.check_method_support(c, "glmnet", 2)
        c1_correlate.check_method_support(c, "twin_count", 2)
        head = ReliabilityHead(name="test")
        head.fit(np.zeros((120, 3)), np.full(120, .1), np.zeros(120), np.ones(120),
                 np.ones(120, bool), ids=np.array([str(i) for i in range(120)]))
        self.assertIn("insufficient_calibration_audit_support", head.status)
        self.assertIsNone(head.error_model)

    def test_zero_cases_prot_does_not_block_met_or_associations(self):
        d, x, truth, cfg = synthetic()
        prot = {k: v.copy() for k, v in d.items()}
        prot["event"][:] = 0
        prot["time"][:] = 14
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, WORKFLOW_WORKER="1"), \
             patch.object(c1_correlate, "native_association") as association, threadpool_limits(1):
            root = Path(td)
            sources = {"prot": bundle(root, prot, x, truth, "prot"), "met": bundle(root, d, x, truth, "met")}
            argv = self.arguments(root, cfg, sources["prot"][0], "abm_count")
            i = argv.index("--prepared-source")
            del argv[i:i+2]
            argv[argv.index("--biom")+1] = "prot,met"
            with patch.object(c1_correlate, "prepare", side_effect=lambda args, *unused: sources[args.biom]):
                self.assertEqual(c1_correlate.main(argv), 0)
            self.assertEqual([c.args[0].biom for c in association.call_args_list], ["prot", "met"])
            indexes = {b: json.loads((root / f"results/toy/{b}/c1_correlate/C1_FINAL_INDEX.json").read_text())
                       for b in ("prot", "met")}
            self.assertEqual(indexes["prot"]["methods"], {})
            self.assertIn("abm_count", indexes["prot"]["unavailable_methods"])
            self.assertIn("abm_count", indexes["met"]["methods"])

    def test_failed_method_remains_error_but_count_is_preserved(self):
        d, x, truth, cfg = synthetic()
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, WORKFLOW_WORKER="1"), \
             patch.object(c1_correlate, "native_association"), threadpool_limits(1):
            root = Path(td)
            source, _, _ = bundle(root, d, x, truth)
            with patch.object(c1_correlate, "fit_glmnet", side_effect=RuntimeError("synthetic genuine failure")):
                rc = c1_correlate.main(self.arguments(root, cfg, source))
            self.assertEqual(rc, 1)
            index = json.loads((root / "results/toy/prot/c1_correlate/C1_FINAL_INDEX.json").read_text())
            self.assertIn("glmnet", index["failed_methods"])
            self.assertIn("abm_count", index["methods"])

    def test_unavailable_current_method_does_not_reuse_old_model(self):
        with tempfile.TemporaryDirectory() as td:
            scope = Path(td)
            old = scope / "methods/glmnet/old"
            old.mkdir(parents=True)
            for name in ("COMPLETE.json", "test_predictions_LOCAL_ONLY.joblib", "calibration_audit_predictions_LOCAL_ONLY.joblib"):
                (old / name).write_text("old")
            (scope / "C1_FINAL_INDEX.json").write_text(json.dumps(dict(methods={}, unavailable_methods={"glmnet": {}})))
            self.assertIsNone(locate_method(scope, "glmnet"))
            tables, contract = c1_report.reliability_tables(scope)
            self.assertEqual(set(tables), {"reliability_status"})
            self.assertTrue(tables["reliability_status"].status.eq("unavailable_glmnet").all())
            self.assertEqual(c1_report.reliability_summary_row("toy",tables,contract)["primary_conclusion"],
                             "primary_unavailable")


class ScopeRunnerTests(unittest.TestCase):
    def test_native_scope_stages_do_not_retry_failed_prerequisites(self):
        script = (ROOT / "f/0.engine.sh").read_text()
        function = script[script.index("plan_job() {"):script.index("\nIFS=',' read -ra raw_traits")]
        command = ("set -eu\ndeclare -A completed_results=() job_actions=()\nreplace=FALSE\n"
                   "layer_complete() { return 1; }\n" + function +
                   "\nplan_job toy c3_coloc prot\nplan_job toy c4_explain prot\n"
                   "printf '%s\\n' \"${!job_actions[@]}\"\n")
        for scope_worker in ("0", "1"):
            result = subprocess.run(["bash", "-c", command], text=True, capture_output=True,
                                    env={**os.environ, "LE8_SCOPE_WORKER": scope_worker})
            self.assertEqual(result.returncode, 0, result.stderr)
            planned = {line.split("|")[1] for line in result.stdout.splitlines()}
            self.assertEqual(planned, {"c3_coloc", "c4_explain"} if scope_worker == "1" else
                {"c1_correlate", "c2_cause", "c3_coloc", "c4_panel_validation", "c4_explain"})

    def test_aggregate_runs_once_after_every_scope_without_extra_job(self):
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, submit="local"):
            root = Path(td)
            output = root / "aggregate.txt"
            steps = [dict(key="aggregate", aggregate=True, argv=[sys.executable, "-c",
                     "from pathlib import Path; import sys; Path(sys.argv[1]).open('a').write('done\\n')", str(output)])]
            scopes = ["toy.prot", "toy.met"]
            self.assertEqual(scope_workflow.run_steps(root, scopes[0], steps, batch_id="batch", scopes=scopes), 0)
            self.assertFalse(output.exists())
            self.assertEqual(scope_workflow.run_steps(root, scopes[1], steps, batch_id="batch", scopes=scopes), 0)
            self.assertEqual(output.read_text(), "done\n")
            self.assertFalse(scope_workflow.claim_aggregates(root, scopes[0], "batch", scopes))

    def test_cpu_only_two_jobs_and_no_nested_submission(self):
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, submit="hpc"):
            args, extra = hpc_workflow.parse_cli(["--Y", "cvd_cad", "--biom", "prot,met",
                                                 "--device", "cpu", "--analysis-root", td])
            jobs = hpc_workflow.build_jobs(args, extra)
            self.assertEqual([j.key for j in jobs], ["cvd_cad.prot", "cvd_cad.met"])
            for job in jobs:
                steps = json.loads(job.argv[job.argv.index("--steps-json")+1])
                self.assertEqual(len(steps), 12)
                self.assertTrue(all("bsub" not in s["argv"] for s in steps))

    def test_stage_failure_runs_all_later_stages_and_returns_failure(self):
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, submit="local"):
            output = Path(td) / "called.txt"
            steps = []
            for stage in ("c1", "c2", "c3", "c4", "c5"):
                code = ("import os,sys; from pathlib import Path; "
                        "assert os.environ['WORKFLOW_WORKER']=='1' and os.environ['LE8_SCOPE_WORKER']=='1'; "
                        "p=Path(sys.argv[1]); p.open('a').write(sys.argv[2]+'\\n'); "
                        "sys.exit(7 if sys.argv[2]=='c1' else 0)")
                steps.append(dict(key="toy.prot."+stage, argv=[sys.executable, "-c", code, str(output), stage]))
            self.assertEqual(scope_workflow.run_steps(Path(td), "toy.prot", steps), 1)
            self.assertEqual(output.read_text().splitlines(), ["c1", "c2", "c3", "c4", "c5"])
            state = json.loads((Path(td) / "logs/scopes/toy.prot.json").read_text())
            self.assertEqual(state["failed_steps"], ["toy.prot.c1"])


class SparseC5Tests(unittest.TestCase):
    def test_small_gene_does_not_cancel_other_genes(self):
        spec = importlib.util.spec_from_file_location("le8_c5_sparse", ROOT / "f/c5.cellulation.py")
        c5 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(c5)
        rng = np.random.default_rng(5)
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            rows = []
            for gene, n in (("SMALL", 5), ("ENOUGH", 30)):
                donors = [f"d{i}" for i in range(n)]
                cells = ["a", "b"]
                matrices = {"ctp": rng.normal(size=(n, 2)), "ctnu": np.full((n, 2), .1),
                            "P": np.full((n, 2), .5), "K": np.eye(n) + np.full((n, n), .05)}
                row = dict(gene=gene, tissue="blood", build="38", kinship_scope="cis", ctnu_definition="variance_of_pseudobulk_mean")
                for key, values in matrices.items():
                    path = root / f"{gene}_{key}.csv"
                    pd.DataFrame(values, index=donors, columns=donors if key == "K" else cells).to_csv(path, index_label="donor")
                    row[key] = path.name
                rows.append(row)
            manifest = root / "manifest.csv"
            pd.DataFrame(rows).to_csv(manifest, index=False)
            self.assertEqual(c5.run(manifest, root / "out", validate_only=True), 0,
                             (root / "out/cigma.status.csv").read_text())
            status = pd.read_csv(root / "out/cigma.status.csv").set_index("gene")
            self.assertEqual(status.loc["SMALL", "status"], "skipped_insufficient_support")
            self.assertEqual(status.loc["ENOUGH", "status"], "validated")


if __name__ == "__main__":
    unittest.main()
