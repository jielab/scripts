"""Public-entry scheduling contracts and real small-data local orchestration."""
from __future__ import annotations
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

GRID = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(GRID / "f"))
spec = importlib.util.spec_from_file_location("grid_workflow_planner", GRID / "f/hpc_workflow.py")
planner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(planner)


class PlanningTests(unittest.TestCase):
    def args(self, entry, module="all", *options):
        return planner.parser().parse_args(["--entry", entry, "--workflow-module", module,
            "--analysis-root", "/tmp/grid-unit-output", *options])

    def build(self, entry, module="all", options=(), extra=(), mode="hpc"):
        with patch.dict(os.environ, {"submit": mode}):
            return planner.build_jobs(self.args(entry, module, *options), list(extra))

    def test_csx_all_has_each_core_score_and_own_evaluation(self):
        jobs = self.build("csx")
        index = {j.key: j for j in jobs}
        self.assertEqual(len(index), 10)
        self.assertEqual(set(index["shared.prepare"].depends), {t+".score" for t in ("height", "ldl", "t2dm")})
        for trait in ("height", "ldl", "t2dm"):
            self.assertEqual(index[trait+".score"].depends, [trait+".weights"])
            self.assertEqual(index[trait+".targetfit"].depends, ["shared.prepare"])
            self.assertEqual(index[trait+".weights"].argv[2], "weights")
        self.assertNotIn("grid.sh", " ".join(" ".join(j.argv) for j in jobs))

    def test_existing_external_cohort_does_not_start_genotype_preparation(self):
        jobs = self.build("csx", extra=["--cohort", "/tmp/existing.tsv", "--validation-n", "400"])
        self.assertNotIn("shared.prepare", {j.key for j in jobs})
        for job in jobs:
            if job.key.endswith("targetfit"):
                self.assertIn("/tmp/existing.tsv", job.argv)
                self.assertEqual(len(job.depends), 3)
            else:
                self.assertNotIn("--cohort", job.argv)
                self.assertNotIn("--validation-n", job.argv)

    def test_block_and_count_have_distinct_jobs_and_boundary_options(self):
        block = self.build("csx", "blocknet", options=["--traits", "height"],
                           extra=["--blocks", "1", "--ld-block-files", "ALL.37.bed"])
        self.assertEqual([j.key for j in block], ["shared.prepare", "height.blockprepare", "height.blocknet"])
        self.assertEqual(block[-1].resources["gpus"], 1)
        self.assertEqual(block[-1].resources["cores"], 20)
        self.assertIn("ALL.37.bed", block[-2].argv)
        count = self.build("grid", options=["--traits", "height,ldl"],
                           extra=["--input-manifest", "/tmp/{trait}.json"])
        self.assertEqual(len(count), 8)
        self.assertEqual(count[-2].depends, ["height.evaluate", "ldl.evaluate"])
        self.assertEqual(count[-1].depends, ["architecture"])
        self.assertFalse(any("blocknet" in j.key for j in count))

    def test_shared_raw_options_are_not_passed_to_block_or_count_models(self):
        for entry, module in [("grid", "all"), ("csx", "blockridge")]:
            jobs = self.build(entry, module, options=["--traits", "height"],
                              extra=["--pheno-file", "/tmp/phe.rds", "--group-file", "/tmp/families.tsv"])
            self.assertIn("--pheno-file", jobs[0].argv)
            for job in jobs[1:]:
                self.assertNotIn("--pheno-file", job.argv)
                self.assertNotIn("--group-file", job.argv)

    def test_count_uses_declared_producer_config_and_serializes_its_writers(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "custom.json"
            path.write_text(json.dumps(dict(schema_version="grid.count_abm.config.1", block_config="custom-blocks.json")))
            jobs = self.build("grid", options=["--traits", "height,ldl"], extra=["--config", str(path)])
            lookup = {j.key:j for j in jobs}
            for key in ("shared.prepare", "height.blockprepare", "ldl.blockprepare", "height.prepare", "ldl.prepare"):
                argv = lookup[key].argv
                self.assertEqual(argv[argv.index("--block-config")+1], str(path.parent/"custom-blocks.json"))
            self.assertEqual(lookup["ldl.blockprepare"].depends, ["height.blockprepare"])
            self.assertEqual(lookup["height.prepare"].depends, ["ldl.blockprepare"])
            self.assertEqual(lookup["ldl.prepare"].depends, ["ldl.blockprepare"])
            self.assertEqual(planner.count_block_config(["--config", str(path), "--block-config", "grid/config/abm.blocks.json"]),
                             str((Path.cwd()/"grid/config/abm.blocks.json").resolve()))

    def test_local_auto_does_not_force_cuda_or_hpc_ram(self):
        with patch.dict(os.environ, {}, clear=True):
            jobs = self.build("csx", "blocknet", options=["--traits", "height"], mode="local")
        self.assertEqual(jobs[-1].argv[jobs[-1].argv.index("--device")+1], "auto")
        self.assertEqual(jobs[-1].resources["gpus"], 0)
        self.assertEqual(jobs[-1].resources["memory_gib"], 16)
        self.assertEqual(jobs[-1].env["WORKFLOW_DEVICE"], "auto")

    def test_explicit_gpu_and_cpu_requests(self):
        gpu = self.build("grid", "twins", options=["--traits", "height", "--device", "cuda"])[0]
        self.assertEqual(gpu.resources["queue"], "gpu")
        self.assertEqual(gpu.resources["gpus"], 1)
        cpu = self.build("csx", "blocknet", options=["--traits", "height", "--device", "cpu"])[-1]
        self.assertEqual(cpu.resources["gpus"], 0)
        self.assertEqual(cpu.argv[cpu.argv.index("--device")+1], "cpu")
        with self.assertRaisesRegex(ValueError, "requires --gpu-count"):
            self.build("grid", "twins", options=["--device", "cuda", "--gpu-count", "0"])

    def test_eight_threads_and_cuda_radius_do_not_force_gpu_reports(self):
        jobs = self.build("grid", options=["--traits", "height", "--device", "cuda", "--threads", "8"],
                          extra=["--input-manifest", "/tmp/height.json"])
        self.assertEqual(jobs[0].resources["gpus"], 1)
        self.assertEqual(jobs[0].argv[jobs[0].argv.index("--device")+1], "cuda")
        for job in jobs[-2:]:
            self.assertEqual(job.resources["gpus"], 0)
            self.assertEqual(job.env["GRID_THREADS"], "4")

    def test_within_job_parallelism_cannot_exceed_lsf_cores(self):
        with self.assertRaisesRegex(ValueError, "exceeds --job-cores"):
            self.build("csx", "weights", options=["--threads", "16", "--job-cores", "32"], extra=["--jobs", "4"])

    def test_disjoint_method_entries_and_pca_installer(self):
        for entry, module in [("disco", "all"), ("prsformer", "all"), ("pca", "all"), ("install", "all")]:
            jobs = self.build(entry, module)
            self.assertTrue(jobs)
            self.assertFalse(any(j.key.endswith(".twins") for j in jobs))
        install = self.build("install")[0]
        self.assertEqual(install.argv, ["bash", str(GRID/"install.sh")])


class PublicEntryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="grid scheduler ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.env = {**os.environ, "DIR_WORK": str(self.root/"work"), "DIR_DATA": str(self.root/"data"),
                    "DIR_SCRATCH": str(self.root/"scratch"), "SHARED_PYTHON": sys.executable,
                    "WORKFLOW_PYTHON": sys.executable, "GRID_PYTHON": sys.executable,
                    "PYTHONPATH": str(GRID/"f")+os.pathsep+os.environ.get("PYTHONPATH", "")}
        self.env.pop("WORKFLOW_WORKER", None)

    def call(self, filename, arguments, output):
        result = subprocess.run(["bash", str(GRID/filename), *arguments, "--analysis-root", str(output)],
            env=self.env, capture_output=True, text=True, timeout=90)
        self.assertEqual(result.returncode, 0, msg=result.stdout+result.stderr)
        plans = list((output/"jobs").glob("*/plan.json"))
        self.assertEqual(len(plans), 1, msg=result.stdout)
        return json.loads(plans[0].read_text()), plans[0]

    def test_all_public_entry_previews_need_no_data_or_conda(self):
        for filename, module in [("1.pca.sh", []), ("2.csx.sh", []), ("3.disco.sh", []),
                                 ("4.prsformer.sh", []), ("install.sh", []), ("grid.sh", ["all"])]:
            with self.subTest(filename=filename):
                output = self.root/filename
                plan, path = self.call(filename, [*module, "--submit", "hpc", "--plan-only"], output)
                self.assertTrue(plan["jobs"])
                self.assertFalse((path.parent/"submission.json").exists())
                for job in plan["jobs"]:
                    self.assertIn("LSB_JOBID", (path.parent/"commands"/(job["key"]+".cmd")).read_text())

    def test_legacy_stage_and_unknown_argument_value_do_not_become_module(self):
        plan, _ = self.call("2.csx.sh", ["--stage", "weights", "--trait", "height", "--submit", "hpc", "--plan-only"], self.root/"csx")
        self.assertEqual([j["key"] for j in plan["jobs"]], ["height.weights"])
        plan, _ = self.call("grid.sh", ["--input-manifest", "/tmp/{trait}.json", "--traits", "height", "--submit", "hpc", "--plan-only"], self.root/"count")
        self.assertEqual(len(plan["jobs"]), 5)

    def test_scientific_worker_cannot_be_forced_onto_login_node(self):
        env = {k:v for k,v in self.env.items() if not k.startswith("LSB_")}
        env.update(WORKFLOW_WORKER="1", submit="hpc", WORKFLOW_RESOURCE_MANAGER="lsf")
        for filename in ("1.pca.sh", "2.csx.sh", "3.disco.sh", "4.prsformer.sh", "grid.sh", "install.sh"):
            result = subprocess.run(["bash", str(GRID/filename)], env=env, capture_output=True, text=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("allocated workflow .cmd job", result.stderr)

    def test_real_public_local_count_pipeline_completes_all_five_workers(self):
        from test_count_abm import exported_fixture
        manifest, trait = exported_fixture(self.root/"input", maps=1)
        output = self.root/"actual"
        plan, path = self.call("grid.sh", ["all", "--submit", "local", "--jobs-N", "2",
            "--traits", trait, "--input-manifest", os.path.relpath(manifest, Path.cwd()), "--device", "cpu", "--threads", "1",
            "--radius-k-ref", "8", "--min-count", "5", "--min-group-ess", "3", "--bootstrap", "5"], output)
        for job in plan["jobs"]:
            state = json.loads((path.parent/"status"/(job["key"]+".json")).read_text())
            self.assertEqual(state["status"], "DONE")
        self.assertTrue((output/"count_abm/final/index.html").is_file())
        self.assertTrue((output/"count_abm"/trait/"_LOCAL_ONLY/test.predictions.tsv.gz").is_file())


if __name__ == "__main__":
    unittest.main()
