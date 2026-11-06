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

    def test_check_never_schedules_cohort_preparation_or_target_fitting(self):
        for entry in ("csx", "disco"):
            jobs = self.build(entry, extra=["--check"])
            self.assertEqual(len(jobs), 3)
            self.assertTrue(all(job.key.endswith(".check") for job in jobs))
            self.assertTrue(all("--check" in job.argv and not job.depends for job in jobs))
            self.assertTrue(all(job.argv[2] == "all" for job in jobs))
        for entry, module in (("csx", "blocknet"), ("prsformer", "all"), ("grid", "all")):
            with self.assertRaisesRegex(ValueError, "use --plan-only"):
                self.build(entry, module, extra=["--check"])

    def test_chromosome_subset_is_carried_to_its_own_target_evaluation(self):
        for specification in ("21-22", "chr21;CHR22;21", "21,22"):
            jobs = self.build("csx", options=["--traits", "height"], extra=["--chrs", specification])
            fit = next(job for job in jobs if job.key.endswith("targetfit"))
            self.assertEqual(planner.option(fit.argv, {"--score-subset"}), "chr21,22")
        fit = self.build("csx", options=["--traits", "height"])[-1]
        self.assertNotIn("--score-subset", fit.argv)

    def test_block_and_count_have_distinct_jobs_and_boundary_options(self):
        block = self.build("csx", "blocknet", options=["--traits", "height"],
                           extra=["--blocks", "1", "--ld-block-files", "ALL.37.bed"])
        self.assertEqual([j.key for j in block], ["shared.prepare", "height.blockprepare", "height.blocknet"])
        self.assertEqual(block[-1].resources["gpus"], 1)
        self.assertEqual(block[-1].resources["cores"], 8)
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
        self.assertEqual(gpu.resources["queue"], os.environ.get("HPC_GPU_QUEUE", "4a100-40"))
        self.assertEqual(gpu.resources["gpus"], 1)
        cpu = self.build("csx", "blocknet", options=["--traits", "height", "--device", "cpu"])[-1]
        self.assertEqual(cpu.resources["gpus"], 0)
        self.assertEqual(cpu.argv[cpu.argv.index("--device")+1], "cpu")
        with self.assertRaisesRegex(ValueError, "requires --gpu-count"):
            self.build("grid", "twins", options=["--device", "cuda", "--gpu-count", "0"])

    def test_hpc_count_auto_requests_gpu_but_cpu_override_is_preserved(self):
        for stage in ("prepare", "twins"):
            extra = ["--input-manifest", "/tmp/{trait}.json"] if stage == "prepare" else []
            automatic = self.build("grid", stage, options=["--traits", "height"], extra=extra)[0]
            self.assertEqual(automatic.resources["gpus"], 1)
            self.assertEqual(automatic.resources["memory_gib"], 160)
            self.assertEqual(automatic.argv[automatic.argv.index("--device")+1], "cuda")
            cpu = self.build("grid", stage, options=["--traits", "height", "--device", "cpu"], extra=extra)[0]
            self.assertEqual(cpu.resources["gpus"], 0)
            self.assertEqual(cpu.argv[cpu.argv.index("--device")+1], "cpu")

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
        check = self.build("install", extra=["--check"])[0]
        self.assertEqual(check.argv, ["bash", str(GRID/"install.sh"), "--check"])


class UpstreamFileTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="grid prerequisites ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.pca = self.root / "projection/pca.tsv.gz"
        self.ancestry = self.root / "projection/ancestry.tsv.gz"
        self.scores = self.root / "scores"
        environment = patch.dict(os.environ, {
            "submit": "hpc", "GRID_PCA_FILE": str(self.pca),
            "GRID_ANCESTRY_FILE": str(self.ancestry), "GRID_SCORE_DIR": str(self.scores),
            "GRID_GWAS_DIR": str(self.root / "gwas"), "GRID_CHRS": "1-22",
        })
        environment.start()
        self.addCleanup(environment.stop)

    def touch(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("exists\n")  # Only existence is checked on the submission host.
        return path

    def jobs(self, entry, module="all", extra=()):
        args = planner.parser().parse_args(["--entry", entry, "--workflow-module", module,
            "--analysis-root", str(self.root / "analysis"), "--traits", "height"])
        return planner.build_jobs(args, list(extra))

    def projection(self):
        self.touch(self.pca)
        self.touch(self.ancestry)

    def test_csx_all_checks_pca_but_allows_its_own_future_scores(self):
        jobs = self.jobs("csx")
        with self.assertRaisesRegex(ValueError, "Missing prerequisite files; no jobs were submitted") as error:
            planner.check_upstream_files(jobs)
        self.assertIn(str(self.pca), str(error.exception))
        self.assertIn("bash 1.pca.sh", str(error.exception))
        self.projection()
        planner.check_upstream_files(jobs)
        self.assertFalse((self.scores / "height/1csx.scores.rds").exists())

    def test_disco_respects_explicit_paths_and_lists_missing_csx(self):
        pca = self.touch(self.root / "custom/pca.gz")
        ancestry = self.touch(self.root / "custom/ancestry.gz")
        scores = self.root / "custom/scores"
        jobs = self.jobs("disco", extra=["--pca-file", str(pca), "--ancestry-file", str(ancestry),
                                        "--dir-score", str(scores)])
        self.assertEqual(jobs[0].env["GRID_PCA_FILE"], str(pca))
        with self.assertRaises(ValueError) as error:
            planner.check_upstream_files(jobs)
        score = scores / "height/1csx.scores.rds"
        self.assertIn(str(score), str(error.exception))
        self.assertIn("bash 2.csx.sh", str(error.exception))
        self.touch(score)
        planner.check_upstream_files(jobs)

    def test_prsformer_all_needs_projection_without_other_methods_outputs(self):
        jobs = self.jobs("prsformer")
        with self.assertRaises(ValueError):
            planner.check_upstream_files(jobs)
        self.projection()
        planner.check_upstream_files(jobs)
        self.assertFalse(self.scores.exists())

    def test_grid_checks_score_provenance_and_referenced_weight_files(self):
        self.projection()
        self.touch(self.scores / "height/1csx.scores.rds")
        jobs = self.jobs("grid")
        provenance = self.root / "analysis/csx/height/provenance.json"
        with self.assertRaises(ValueError) as error:
            planner.check_upstream_files(jobs)
        self.assertIn(str(provenance), str(error.exception))
        weights = {p: self.root / "weights" / (p + ".gz") for p in ("AFR", "EAS", "EUR", "SAS")}
        provenance.parent.mkdir(parents=True, exist_ok=True)
        provenance.write_text(json.dumps({"models": {"csx."+p: {"weights_file": str(w)} for p,w in weights.items()}}))
        with self.assertRaises(ValueError) as error:
            planner.check_upstream_files(jobs)
        for weight in weights.values():
            self.assertIn(str(weight), str(error.exception))
            self.touch(weight)
        planner.check_upstream_files(jobs)
        self.assertFalse((self.root / "analysis/prsformer/cache/data.tsv.gz").exists())

    def test_score_only_needs_existing_weights_and_signatures(self):
        for pop in ("AFR", "EAS", "EUR", "SAS"):
            self.touch(self.root / "gwas" / f"height.{pop}.gz")
        jobs = self.jobs("csx", "score", ["--chrs", "22"])
        with self.assertRaises(ValueError) as error:
            planner.check_upstream_files(jobs)
        self.assertIn("height.AFR.chr22.csx.gz", str(error.exception))
        for pop in ("AFR", "EAS", "EUR", "SAS"):
            self.touch(self.root / "gwas" / f"height.{pop}.chr22.csx.gz")
            self.touch(self.root / "gwas" / f"height.{pop}.chr22.csx.gz.signature")
        planner.check_upstream_files(jobs)

    def test_standalone_count_stage_reports_the_producing_command(self):
        jobs = self.jobs("grid", "twins")
        with self.assertRaises(ValueError) as error:
            planner.check_upstream_files(jobs)
        receipt = self.root / "analysis/count_abm/height/PREPARE_COMPLETE.json"
        self.assertIn(str(receipt), str(error.exception))
        self.assertIn("bash grid.sh prepare", str(error.exception))
        self.touch(receipt)
        planner.check_upstream_files(jobs)

    def test_external_manifest_is_checked_and_internal_stages_are_not_blocked(self):
        manifest = self.root / "input/height.json"
        jobs = self.jobs("grid", extra=["--input-manifest", str(manifest.parent / "{trait}.json")])
        with self.assertRaises(ValueError) as error:
            planner.check_upstream_files(jobs)
        self.assertIn(str(manifest), str(error.exception))
        self.touch(manifest)
        planner.check_upstream_files(jobs)


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

    def test_hpc_installer_resolves_an_existing_recipe_without_installing(self):
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        conda = bin_dir / "conda"
        conda.write_text("#!/usr/bin/env bash\nexit 0\n")
        conda.chmod(0o755)
        env = dict(self.env, PATH=str(bin_dir)+os.pathsep+os.environ["PATH"],
                   WORKFLOW_WORKER="1", submit="hpc", SUBMIT="hpc", LSB_JOBID="test",
                   WORKFLOW_RESOURCE_MANAGER="lsf", GRID_CONDA_ENV=str(self.root / "env"))
        result = subprocess.run(["bash", str(GRID / "install.sh"), "--dry-run"],
                                env=env, capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn(str(GRID / "environment.hpc.yml"), result.stdout)
        self.assertFalse((self.root / "env").exists())

    def test_default_hpc_submission_returns_without_querying_job_status(self):
        bin_dir = self.root / "fake-lsf"
        bin_dir.mkdir()
        submitted, polled = self.root / "submitted", self.root / "polled"
        for name, body in {
            "bsub": 'cat >/dev/null\nprintf "accepted\\n" >> "$TEST_SUBMITTED"\necho "Job <12345> is submitted to queue <33>."\n',
            "bjobs": 'touch "$TEST_POLLED"\nexit 1\n',
        }.items():
            executable = bin_dir / name
            executable.write_text("#!/usr/bin/env bash\n" + body)
            executable.chmod(0o755)
        env = dict(self.env, PATH=str(bin_dir)+os.pathsep+os.environ["PATH"], SUBMIT_WAIT="FALSE",
                   TEST_SUBMITTED=str(submitted), TEST_POLLED=str(polled))
        result = subprocess.run(["bash", str(GRID / "1.pca.sh"), "--submit", "hpc",
            "--analysis-root", str(self.root / "submit")], env=env, capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('"status": "SUBMITTED"', result.stdout)
        self.assertEqual(submitted.read_text(), "accepted\n")
        self.assertFalse(polled.exists())

    def test_missing_prerequisite_is_printed_before_any_plan_or_submission(self):
        output = self.root / "missing-output"
        pca = self.root / "missing.pca.tsv.gz"
        scores = self.root / "missing-scores"
        bin_dir = self.root / "blocked-lsf"
        bin_dir.mkdir()
        submitted = self.root / "unexpected-submission"
        bsub = bin_dir / "bsub"
        bsub.write_text('#!/usr/bin/env bash\ntouch "$TEST_SUBMITTED"\nexit 1\n')
        bsub.chmod(0o755)
        env = dict(self.env, GRID_PCA_FILE=str(pca), GRID_SCORE_DIR=str(scores),
                   GRID_ANCESTRY_FILE=str(self.root / "missing.ancestry.tsv.gz"), SUBMIT_WAIT="FALSE",
                   PATH=str(bin_dir)+os.pathsep+os.environ["PATH"], TEST_SUBMITTED=str(submitted))
        result = subprocess.run(["bash", str(GRID / "3.disco.sh"), "--traits", "height",
            "--submit", "hpc", "--analysis-root", str(output)], env=env,
            capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        self.assertIn("[MISSING] " + str(pca), result.stderr)
        self.assertIn(str(scores / "height/1csx.scores.rds"), result.stderr)
        self.assertIn("bash 1.pca.sh", result.stderr)
        self.assertIn("bash 2.csx.sh", result.stderr)
        self.assertNotIn("[SUBMITTED]", result.stdout)
        self.assertFalse(submitted.exists())
        self.assertFalse((output / "jobs").exists())

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
