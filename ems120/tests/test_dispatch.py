"""Subprocess tests for EMS120 routing; no R, CUDA, LSF or real data required."""
from __future__ import annotations

import json
import fcntl
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


PROJECT = Path(__file__).resolve().parents[1]
ENTRY = PROJECT / "ems120.sh"


class DispatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="ems120-dispatch-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.stubs = self.root / "bin"
        self.stubs.mkdir()
        self.calls = self.root / "unexpected-runtime-calls"
        self.jobs = self.root / "plans with spaces"
        self.env = dict(os.environ)
        for name in list(self.env):
            if name.startswith(("EMS120_", "DIR_", "_OS_DEFAULT_", "HPC_", "WORKFLOW_", "LSB_")):
                self.env.pop(name)
        for name in ("submit", "SUBMIT", "SUBMIT_WAIT", "jobs_N", "SCRIPTS_ROOT",
                     "ANALYSIS_ROOT", "PERSISTENT_ANALYSIS_ROOT", "DATE_SCRATCH",
                     "CONDA_EXE", "CONDA_PREFIX", "RETICULATE_PYTHON", "SHARED_PYTHON"):
            self.env.pop(name, None)
        self.env["DIR_WORK"] = str(self.root / "work")
        self.env["DIR_TMP"] = str(self.root / "disposable")
        self.env["EMS120_PLANNER_PYTHON"] = sys.executable
        self.env.pop("CUDA_VISIBLE_DEVICES", None)
        self.env["PATH"] = str(self.stubs) + os.pathsep + os.environ.get("PATH", "/usr/bin:/bin")
        self.env["EMS_TEST_RUNTIME_CALLS"] = str(self.calls)
        self.fail_conda = self._executable("conda", """#!/usr/bin/env python3
import os, pathlib, sys
with pathlib.Path(os.environ['EMS_TEST_RUNTIME_CALLS']).open('a') as out:
    out.write('conda ' + repr(sys.argv[1:]) + '\\n')
sys.exit(93)
""")
        self.fail_r = self._executable("Rscript", """#!/usr/bin/env python3
import os, pathlib, sys
with pathlib.Path(os.environ['EMS_TEST_RUNTIME_CALLS']).open('a') as out:
    out.write('Rscript ' + repr(sys.argv[1:]) + '\\n')
sys.exit(94)
""")

    def _executable(self, name, text):
        path = self.stubs / name
        path.write_text(text, encoding="utf-8")
        path.chmod(0o755)
        return path

    def call(self, *args, check=True, env=None, workflow="single"):
        command = ["bash", str(ENTRY), "--dir-jobs", str(self.jobs),
                   "--conda-exe", str(self.fail_conda), "--rscript", str(self.fail_r), *map(str, args)]
        if workflow is not None:
            command.extend(["--workflow", workflow])
        proc = subprocess.run(command, env=env or self.env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=20)
        if check:
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        return proc

    def plan(self, *args, workflow="single"):
        before = set(self.jobs.glob("*/plan.json"))
        result = self.call(*args, "--plan-only", workflow=workflow)
        after = set(self.jobs.glob("*/plan.json"))
        self.assertEqual(len(after - before), 1, result.stdout + result.stderr)
        path = (after - before).pop()
        data = json.loads(path.read_text())
        return path, data, data["jobs"][0]

    def assert_no_runtime(self):
        self.assertFalse(self.calls.exists(), self.calls.read_text() if self.calls.exists() else "")

    def fake_worker_environment(self, fail_runtime=False):
        """A real shell activation and process chain, with a recording R stub."""
        base = self.root / "fake conda"
        prefix = self.root / "fake ai environment"
        profile = base / "etc/profile.d/conda.sh"
        profile.parent.mkdir(parents=True)
        (prefix / "bin").mkdir(parents=True)
        (prefix / "lib").mkdir()
        self.env["EMS_TEST_CONDA_BASE"] = str(base)
        self.env["EMS_TEST_CONDA_PREFIX"] = str(prefix)
        self.env["EMS_TEST_R_LOG"] = str(self.root / "r-calls.jsonl")
        self.env["EMS_TEST_FAIL_RUNTIME"] = "1" if fail_runtime else "0"
        self.env["EMS120_DATA_DIR"] = str(self.root / "raw inputs")
        self.env["EMS120_ANALYSIS_DIR"] = str(self.root / "worker analysis")
        (Path(self.env["EMS120_DATA_DIR"]) / "clean").mkdir(parents=True)
        profile.write_text('''conda() {
  if [[ "$1" != activate ]]; then return 2; fi
  export CONDA_PREFIX="$EMS_TEST_CONDA_PREFIX"
  export PATH="$CONDA_PREFIX/bin:$PATH"
}
''')
        self.fail_conda.write_text('''#!/usr/bin/env bash
if [[ "${1:-}" == info && "${2:-}" == --base ]]; then
  printf '%s\\n' "$EMS_TEST_CONDA_BASE"
else
  exit 3
fi
''')
        py = prefix / "bin/python"
        # Executable alias keeps the actual parser/stdlib available without ML imports.
        py.symlink_to(sys.executable)
        self.fail_r.write_text('''#!/usr/bin/env python3
import json, os, pathlib, sys
names = ('RETICULATE_PYTHON', 'EMS120_ANALYSIS_DIR', 'EMS120_DATA_DIR',
         'EMS120_BASE_MODEL_DIR', 'EMS120_MODEL_DIR', 'WORKFLOW_CORES',
         'CUDA_VISIBLE_DEVICES', 'PNG_PYTHON', 'TMPDIR', 'TMP', 'TEMP',
         'LC_ALL', 'EMS120_PHASE', 'EMS120_CORES', 'EMS120_REQUIRE_GPU', 'EMS120_PHONE_GPU')
with pathlib.Path(os.environ['EMS_TEST_R_LOG']).open('a') as out:
    out.write(json.dumps({'argv': sys.argv[1:], 'env': {k: os.environ.get(k) for k in names}}) + '\\n')
if os.environ.get('EMS_TEST_FAIL_RUNTIME') == '1' and any(x.endswith('/check_runtime.R') for x in sys.argv):
    print('Synthetic R/reticulate CUDA failure', file=sys.stderr)
    sys.exit(42)
print('Synthetic R runtime passed')
''')
        return prefix

    def test_single_hpc_plan_uses_gpu_and_persistent_analysis(self):
        self.env["CUDA_VISIBLE_DEVICES"] = "7"
        path, plan, job = self.plan()
        self.assertEqual(plan["mode"], "hpc")
        self.assertEqual(job["resources"]["queue"], "2v100-32-sc")
        self.assertEqual(job["resources"]["cores"], 8)
        self.assertEqual(job["resources"]["memory_gib"], 96)
        self.assertEqual(job["resources"]["hours"], 120)
        self.assertEqual(job["resources"]["gpus"], 1)
        self.assertEqual(Path(job["env"]["EMS120_ANALYSIS_DIR"]), self.root / "work/analysis/ems120")
        self.assertEqual(Path(job["env"]["EMS120_DATA_DIR"]), self.root / "work/phe/ems120")
        self.assertEqual(Path(job["env"]["EMS120_BASE_MODEL_DIR"]), self.root / "work/data/AI_model/hfl")
        self.assertEqual(Path(job["env"]["EMS120_SCRIPT_DIR"]), PROJECT)
        self.assertNotIn("CUDA_VISIBLE_DEVICES", job["env"])
        text = (path.parent / "commands" / (job["key"] + ".cmd")).read_text()
        self.assertIn('#BSUB -gpu "num=1/host:', text)
        self.assertIn("span[hosts=1]", text)
        self.assertIn("rusage[mem=96GB/host]", text)
        self.assert_no_runtime()

    def test_local_root_overrides_propagate_to_plan(self):
        dat = self.root / "data files"
        big = self.root / "large data"
        scratch = self.root / "scratch area"
        _, plan, job = self.plan("--submit", "local", "--dir-dat", dat,
                                 "--dir-big", big, "--dir-tmp", scratch)
        self.assertEqual(plan["mode"], "local")
        env = job["env"]
        self.assertEqual(Path(env["EMS120_DATA_DIR"]), dat / "ems120")
        self.assertEqual(Path(env["EMS120_BASE_MODEL_DIR"]), big / "AI_model/hfl")
        self.assertEqual(Path(env["EMS120_FILES_DIR"]), dat / "files")
        self.assertEqual(Path(env["EMS120_ANALYSIS_DIR"]), dat / "analysis/ems120")
        self.assertEqual(Path(env["EMS120_MODEL_DIR"]), dat / "analysis/ems120/ml")
        self.assert_no_runtime()

    def test_project_overrides_preserve_spaces_and_stage_settings(self):
        analysis = self.root / "my analysis"
        models = self.root / "trained models"
        data = self.root / "inputs with spaces"
        base = self.root / "base model"
        files = self.root / "artwork files"
        _, _, job = self.plan("--dir-analysis", analysis, "--dir-model", models,
                              "--dir-data", data, "--dir-base-model", base,
                              "--dir-files", files, "--start-step", "fig2",
                              "--end-step", "fig3", "--run-cv", "FALSE")
        expected = dict(EMS120_ANALYSIS_DIR=analysis, EMS120_MODEL_DIR=models,
                        EMS120_DATA_DIR=data, EMS120_BASE_MODEL_DIR=base, EMS120_FILES_DIR=files)
        for name, value in expected.items():
            self.assertEqual(Path(job["env"][name]), value)
        self.assertEqual(job["env"]["START_STEP"].lower(), "fig2")
        self.assertEqual(job["env"]["END_STEP"].lower(), "fig3")
        self.assertEqual(job["env"]["EMS120_RUN_DX_CV"], "0")
        self.assert_no_runtime()

    def test_gpu_check_plan_is_a_real_gpu_job(self):
        _, plan, job = self.plan("--check-gpu")
        self.assertEqual(plan["mode"], "hpc")
        self.assertEqual(job["key"], "check_gpu")
        self.assertEqual(job["resources"]["gpus"], 1)
        self.assert_no_runtime()

    def test_preflight_checks_inputs_without_running_production_analysis(self):
        self.fake_worker_environment()
        self.env["LC_ALL"] = "C"
        result = self.call("--submit", "local", "--preflight", "--dir-analysis", self.root / "analysis")
        records = [json.loads(line) for line in Path(self.env["EMS_TEST_R_LOG"]).read_text().splitlines()]
        self.assertEqual([Path(r["argv"][1]).name for r in records], ["check_runtime.R", "preflight.R"])
        self.assertRegex(records[0]["env"]["LC_ALL"], r"(?i)utf-?8")
        plan = json.loads(next(self.jobs.glob("*/plan.json")).read_text())
        self.assertEqual(plan["jobs"][0]["key"], "preflight")
        self.assertEqual(plan["jobs"][0]["resources"]["hours"], 1)

    def test_analysis_batch_and_renderer_settings_are_saved_in_plan(self):
        settings = {
            "EMS120_SENS_N_PERM": "25",
            "EMS120_DX_BATCH_SIZE": "128",
            "EMS120_GEO_BATCH_SIZE": "256",
            "EMS120_PHONE_BATCH_SIZE": "512",
            "CHROME_BIN": str(self.root / "browser executable"),
            "CHROME_PATH": str(self.root / "alternate browser"),
            "EDGE_BIN": str(self.root / "edge executable"),
            "FIG2D_PYTHON": str(self.root / "map python"),
        }
        self.env.update(settings)
        _, _, job = self.plan()
        for name, expected in settings.items():
            self.assertEqual(job["env"].get(name), expected, name)
        self.assert_no_runtime()

    def test_hpc_command_refuses_direct_login_execution(self):
        path, _, job = self.plan()
        command = path.parent / "commands" / (job["key"] + ".cmd")
        proc = subprocess.run(["bash", str(command)], env=self.env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10)
        self.assertNotEqual(proc.returncode, 0)
        self.assertRegex(proc.stdout + proc.stderr, r"(?i)(login|bsub|allocation)")
        self.assert_no_runtime()

    def test_hpc_submission_records_receipt_without_loading_science(self):
        captured = self.root / "submitted.cmd"
        self.env["EMS_TEST_BSUB_CAPTURE"] = str(captured)
        self._executable("bsub", """#!/usr/bin/env python3
import os, pathlib, sys
pathlib.Path(os.environ['EMS_TEST_BSUB_CAPTURE']).write_text(sys.stdin.read())
print('Job <837120> is submitted to queue <4a100-40>.')
""")
        result = self.call("--check-gpu")
        self.assertTrue(captured.exists(), result.stdout + result.stderr)
        self.assertIn("#BSUB -gpu", captured.read_text())
        receipts = list(self.jobs.glob("*/submission.json"))
        self.assertEqual(len(receipts), 1)
        receipt = json.loads(receipts[0].read_text())
        self.assertEqual(receipt["jobs"]["check_gpu"]["job_id"], "837120")
        self.assert_no_runtime()

    def test_failed_bsub_does_not_fall_back_to_local(self):
        self._executable("bsub", "#!/usr/bin/env bash\necho 'queue disabled' >&2\nexit 1\n")
        result = self.call("--check-gpu", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assert_no_runtime()

    def test_multi_gpu_and_hpc_cpu_fallback_are_rejected(self):
        for args in (("--gpus", "2"), ("--gpu", "FALSE"), ("--cores", "0"),
                     ("--memory-gb", "-1"), ("--hours", "0"), ("--submit", "auto")):
            with self.subTest(args=args):
                result = self.call(*args, "--plan-only", check=False)
                self.assertNotEqual(result.returncode, 0)
        self.assert_no_runtime()

    def test_explicit_local_cpu_plan_hides_gpus_and_disables_phone_gpu(self):
        _, plan, job = self.plan("--submit", "local", "--gpu", "FALSE")
        self.assertEqual(plan["mode"], "local")
        self.assertEqual(job["resources"]["gpus"], 0)
        self.assertEqual(job["env"]["EMS120_REQUIRE_GPU"], "0")
        self.assertEqual(job["env"]["EMS120_PHONE_GPU"], "off")
        self.assert_no_runtime()

    def test_gpu_runtime_failure_stops_before_analysis(self):
        self.fake_worker_environment(fail_runtime=True)
        result = self.call("--submit", "local", "--dir-analysis", self.root / "analysis", check=False)
        self.assertNotEqual(result.returncode, 0)
        log = Path(self.env["EMS_TEST_R_LOG"])
        self.assertTrue(log.is_file(), result.stdout + result.stderr)
        records = [json.loads(line) for line in log.read_text().splitlines()]
        self.assertEqual(len(records), 1)
        self.assertFalse(any(str(PROJECT / "f/ems120.R") in r["argv"] for r in records))
        statuses = [json.loads(p.read_text()) for p in self.jobs.glob("*/status/*.json")]
        self.assertEqual(len(statuses), 1)
        self.assertEqual(statuses[0]["status"], "EXIT")

    def test_local_worker_uses_configured_environment_and_preserves_device_visibility(self):
        prefix = self.fake_worker_environment()
        self.env["CUDA_VISIBLE_DEVICES"] = "4"
        analysis = self.root / "analysis with spaces"
        scratch = self.root / "job temporary files"
        result = self.call("--submit", "local", "--dir-analysis", analysis,
                           "--dir-tmp", scratch, "--cores", "3", "--start-step", "fig2",
                           "--end-step", "fig3", "--run-cv", "FALSE")
        log = Path(self.env["EMS_TEST_R_LOG"])
        self.assertTrue(log.is_file(), result.stdout + result.stderr)
        records = [json.loads(line) for line in log.read_text().splitlines()]
        main = [r for r in records if str(PROJECT / "f/ems120.R") in r["argv"]]
        self.assertEqual(len(main), 1, records)
        self.assertEqual(main[0]["argv"][-4:], ["--start_step", "fig2", "--end_step", "fig3"])
        env = main[0]["env"]
        self.assertEqual(Path(env["RETICULATE_PYTHON"]), prefix / "bin/python")
        self.assertEqual(Path(env["EMS120_ANALYSIS_DIR"]), analysis)
        self.assertEqual(env["WORKFLOW_CORES"], "3")
        self.assertEqual(env["CUDA_VISIBLE_DEVICES"], "4")
        self.assertEqual(env["PNG_PYTHON"], env["RETICULATE_PYTHON"])
        plans = list(self.jobs.glob("*/plan.json"))
        self.assertEqual(len(plans), 1)
        self.assertEqual(json.loads(plans[0].read_text())["jobs"][0]["env"]["CUDA_VISIBLE_DEVICES"], "4")
        expected_temporary = analysis / "tmp" / plans[0].parent.name / "analysis"
        for name in ("TMPDIR", "TMP", "TEMP"):
            self.assertEqual(Path(env[name]), expected_temporary)
        self.assertTrue(expected_temporary.is_dir())
        self.assertEqual(plans[0].parent.parent, self.jobs)
        statuses = [json.loads(p.read_text()) for p in self.jobs.glob("*/status/*.json")]
        self.assertEqual(len(statuses), 1)
        self.assertEqual(statuses[0]["status"], "DONE")

    def test_local_unset_device_visibility_defaults_to_one_gpu(self):
        self.fake_worker_environment()
        self.env.pop("CUDA_VISIBLE_DEVICES", None)
        self.call("--submit", "local", "--check-gpu")
        records = [json.loads(line) for line in Path(self.env["EMS_TEST_R_LOG"]).read_text().splitlines()]
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["env"]["CUDA_VISIBLE_DEVICES"], "0")

    def test_cache_only_worker_can_redraw_without_original_excel_directory(self):
        self.fake_worker_environment()
        raw = Path(self.env["EMS120_DATA_DIR"])
        (raw / "clean").rmdir()
        raw.rmdir()
        analysis = self.root / "cached analysis"
        (analysis / "dat").mkdir(parents=True)
        # Only the shell's presence/size check runs; fake R never parses this fixture.
        (analysis / "dat/dat0.list.rds").write_bytes(b"nonempty synthetic cache marker")
        self.call("--submit", "local", "--dir-analysis", analysis,
                  "--start-step", "fig2", "--end-step", "fig2", "--run-cv", "FALSE")
        records = [json.loads(line) for line in Path(self.env["EMS_TEST_R_LOG"]).read_text().splitlines()]
        self.assertEqual(sum(str(PROJECT / "f/ems120.R") in r["argv"] for r in records), 1)

    def test_existing_analysis_or_shared_model_lock_blocks_new_worker(self):
        for which in ("analysis", "model"):
            with self.subTest(lock=which):
                analysis = self.root / (which + " analysis")
                model = self.root / (which + " shared model")
                analysis.mkdir()
                model.mkdir()
                locked = analysis if which == "analysis" else model
                with (locked / ".ems120.lock").open("w") as held:
                    fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    result = self.call("--submit", "local", "--dir-analysis", analysis,
                                       "--dir-model", model, check=False)
                self.assertNotEqual(result.returncode, 0)
                errors = "\n".join(p.read_text() for p in self.jobs.glob("*/logs/*.err"))
                self.assertIn("Another EMS120 job is using " + which + " directory", errors)
        self.assert_no_runtime()


    def test_default_split_resources_names_and_dependencies(self):
        self.env["CUDA_VISIBLE_DEVICES"] = "9"
        path, plan, _ = self.plan(workflow=None)
        prep, ml, report = plan["jobs"]
        self.assertEqual([j["key"] for j in plan["jobs"]],
                         ["data.cpu_prepare", "ml_phone_geo_dx.gpu_ml", "statistics.cpu_report"])
        self.assertEqual([j["depends"] for j in plan["jobs"]], [[], [prep["key"]], [ml["key"]]])
        self.assertEqual(len({j["env"]["EMS120_HANDOFF_DIR"] for j in plan["jobs"]}), 1)
        for j in (prep, report):
            self.assertEqual(j["resources"]["gpus"], 0)
            self.assertEqual(j["resources"]["queue"], "33")
            self.assertEqual(j["resources"]["cores"], 20)
            self.assertEqual(j["env"]["EMS120_CORES"], "8")
            self.assertEqual(j["env"]["CUDA_VISIBLE_DEVICES"], "")
            self.assertEqual(j["env"]["EMS120_REQUIRE_GPU"], "0")
            self.assertNotIn("#BSUB -gpu", (path.parent / "commands" / (j["key"] + ".cmd")).read_text())
        self.assertEqual(ml["resources"]["gpus"], 1)
        self.assertEqual(ml["env"]["EMS120_PHONE_GPU"], "required")
        self.assertNotIn("CUDA_VISIBLE_DEVICES", ml["env"])
        self.assert_no_runtime()

    def test_gpu_upgrade_does_not_change_cpu_or_scientific_settings(self):
        _, plan, _ = self.plan("--gpu-queue", "2h100-80", "--gpu-model", "H100",
                               "--gpu-cores", "16", "--gpu-memory-gb", "192", "--gpu-hours", "200", workflow="split")
        prep, ml, report = plan["jobs"]
        self.assertEqual(ml["resources"]["gpu_model"], "H100")
        self.assertEqual(ml["resources"]["cores"], 16)
        self.assertEqual(ml["resources"]["memory_gib"], 192)
        self.assertEqual(ml["resources"]["hours"], 200)
        self.assertEqual(prep["resources"], report["resources"])
        self.assertEqual(prep["resources"]["cores"], 20)
        self.assertEqual(prep["env"]["EMS120_CORES"], "8")
        self.assertEqual(ml["env"]["EMS120_RUN_DX_CV"], "1")
        self.assertNotIn("EMS120_DX_BATCH_SIZE", ml["env"])
        self.assert_no_runtime()

    def test_report_plan_and_worker_skip_gpu_check_and_raw_inputs(self):
        self.fake_worker_environment()
        self.env["CUDA_VISIBLE_DEVICES"] = "4"
        self.env["EMS_TEST_FAIL_RUNTIME"] = "1"  # CPU report must never call it.
        analysis = self.root / "reports"
        (analysis / "dat").mkdir(parents=True)
        (analysis / "dat/dat.list.rds").write_bytes(b"synthetic R stub fixture")
        (Path(self.env["EMS120_DATA_DIR"]) / "clean").rmdir()
        self.call("--submit", "local", "--dir-analysis", analysis, workflow="report")
        records = [json.loads(s) for s in Path(self.env["EMS_TEST_R_LOG"]).read_text().splitlines()]
        self.assertEqual(len(records), 1)
        self.assertTrue(any(x.endswith("/ems120.R") for x in records[0]["argv"]))
        self.assertEqual(records[0]["env"]["CUDA_VISIBLE_DEVICES"], "")
        self.assertEqual(records[0]["env"]["EMS120_PHASE"], "cpu_report")
        _, plan, job = self.plan(workflow="report")
        self.assertEqual(len(plan["jobs"]), 1)
        self.assertEqual(job["resources"]["gpus"], 0)
        self.assertEqual(job["depends"], [])

    def test_missing_report_cache_never_initializes_ml(self):
        self.fake_worker_environment()
        result = self.call("--submit", "local", "--dir-analysis", self.root / "absent", workflow="report", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(Path(self.env["EMS_TEST_R_LOG"]).exists())

    def test_gpu_failure_blocks_cpu_report_in_split_chain(self):
        self.fake_worker_environment(fail_runtime=True)
        result = self.call("--submit", "local", "--dir-analysis", self.root / "split", workflow="split", check=False)
        self.assertNotEqual(result.returncode, 0)
        records = [json.loads(s) for s in Path(self.env["EMS_TEST_R_LOG"]).read_text().splitlines()]
        self.assertEqual([r["env"]["EMS120_PHASE"] for r in records], ["cpu_prepare", "gpu_ml"])
        self.assertEqual(records[0]["env"]["CUDA_VISIBLE_DEVICES"], "")
        states = {p.stem: json.loads(p.read_text())["status"] for p in self.jobs.glob("*/status/*.json")}
        self.assertEqual(states["data.cpu_prepare"], "DONE")
        self.assertEqual(states["ml_phone_geo_dx.gpu_ml"], "EXIT")
        self.assertNotIn("statistics.cpu_report", states)

    def test_split_workers_hide_cuda_only_on_cpu_and_keep_requested_threads(self):
        self.fake_worker_environment()
        self.env["CUDA_VISIBLE_DEVICES"] = "4"
        analysis = self.root / "successful split"
        (analysis / "dat").mkdir(parents=True)
        (analysis / "dat/dat.list.rds").write_bytes(b"R stub fixture")
        self.call("--submit", "local", "--dir-analysis", analysis, "--gpu-cores", "3",
                  "--cpu-cores", "4", "--cpu-threads", "2", workflow="split")
        records = [json.loads(s) for s in Path(self.env["EMS_TEST_R_LOG"]).read_text().splitlines()]
        main = [r for r in records if any(x.endswith("/ems120.R") for x in r["argv"])]
        self.assertEqual([r["env"]["EMS120_PHASE"] for r in main], ["cpu_prepare", "gpu_ml", "cpu_report"])
        self.assertEqual([r["env"]["CUDA_VISIBLE_DEVICES"] for r in main], ["", "4", ""])
        self.assertEqual([r["env"]["EMS120_CORES"] for r in main], ["2", "3", "2"])
        self.assertEqual(sum(any(x.endswith("/check_runtime.R") for x in r["argv"]) for r in records), 1)

    def test_cpu_queue_capacity_and_gpu_separation_are_validated(self):
        for args in (("--cpu-cores", "8"), ("--cpu-cores", "10", "--cpu-memory-gb", "96"),
                     ("--cpu-threads", "21"), ("--cpu-queue", "2v100-32-sc")):
            with self.subTest(args=args):
                result = self.call(*args, "--plan-only", workflow="split", check=False)
                self.assertNotEqual(result.returncode, 0)
        self.assert_no_runtime()

    def test_lsf_split_submission_uses_success_dependencies(self):
        capture = self.root / "submissions.jsonl"
        self.env["EMS_TEST_BSUB_CAPTURE"] = str(capture)
        self._executable("bsub", """#!/usr/bin/env python3
import json, os, pathlib, sys
p = pathlib.Path(os.environ['EMS_TEST_BSUB_CAPTURE'])
n = len(p.read_text().splitlines()) if p.exists() else 0
with p.open('a') as f:
    f.write(json.dumps({'args': sys.argv[1:], 'cmd': sys.stdin.read()}) + '\\n')
print('Job <%d> is submitted to queue <test>.' % (91000 + n))
""")
        self.call(workflow="split")
        records = [json.loads(s) for s in capture.read_text().splitlines()]
        self.assertEqual(len(records), 3)
        self.assertNotIn("-w", records[0]["args"])
        self.assertEqual(records[1]["args"], ["-w", "done(91000) && ended(91000)"])
        self.assertEqual(records[2]["args"], ["-w", "done(91001) && ended(91001)"])
        self.assertIn("#BSUB -gpu", records[1]["cmd"])
        self.assertNotIn("#BSUB -gpu", records[0]["cmd"] + records[2]["cmd"])
        self.assert_no_runtime()


if __name__ == "__main__":
    unittest.main(verbosity=2)
