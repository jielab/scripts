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

SHARED = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SHARED))
from workflow_jobs import Job, normalized, resources, run_plan, write_plan
from workflow_paths import relocate_defaults


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="workflow test ")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def job(self, key, code, depends=(), env=None):
        return Job(key, [sys.executable, "-c", code], list(depends),
                   {"cores": 1, "memory_gib": 1}, env or {}, str(self.root))

    def test_dependency_validation(self):
        for jobs in ([self.job("x", "pass", ["missing"])],
                     [self.job("x", "pass", ["y"]), self.job("y", "pass", ["x"])],
                     [self.job("x", "pass"), self.job("x", "pass")]):
            with self.assertRaises(ValueError):
                normalized(jobs)

    def test_topological_order(self):
        jobs = [self.job("c", "pass", ["a", "b"]), self.job("a", "pass"), self.job("b", "pass")]
        self.assertEqual([j["key"] for j in normalized(jobs)], ["a", "b", "c"])

    def test_local_dag_and_resume(self):
        jobs = [self.job("a", "from pathlib import Path; Path('a').write_text('A')"),
                self.job("b", "from pathlib import Path; Path('b').write_text('B')"),
                self.job("merge", "from pathlib import Path; Path('result').write_text(Path('a').read_text()+Path('b').read_text())", ["a", "b"])]
        plan = write_plan(jobs, self.root, mode="local")
        result = run_plan(plan, submit="local", jobs=2)
        self.assertEqual(result["status"], "DONE")
        self.assertEqual((self.root / "result").read_text(), "AB")
        before = (self.root / "result").stat().st_mtime_ns
        run_plan(plan, submit="local", jobs=2)
        self.assertEqual(before, (self.root / "result").stat().st_mtime_ns)

    def test_failure_does_not_run_dependent(self):
        jobs = [self.job("fail", "raise SystemExit(7)"),
                self.job("merge", "from pathlib import Path; Path('wrong').touch()", ["fail"])]
        plan = write_plan(jobs, self.root, mode="local")
        with self.assertRaises(RuntimeError):
            run_plan(plan, submit="local", jobs=2)
        self.assertFalse((self.root / "wrong").exists())
        self.assertEqual(json.loads((plan.parent / "status/fail.json").read_text())["returncode"], 7)

    def test_literal_arguments_and_environment(self):
        literal = "$(touch SHOULD_NOT_EXIST) `touch ALSO_NOT` ' quote; amp&"
        job = Job("literal", [sys.executable, "-c", "import os,sys; from pathlib import Path; Path('seen').write_text(sys.argv[1]+os.environ['EXAMPLE'])", literal],
                  resources={"cores": 1}, env={"EXAMPLE": literal}, cwd=str(self.root))
        plan = write_plan([job], self.root, mode="local")
        run_plan(plan, submit="local", jobs=1)
        self.assertEqual((self.root / "seen").read_text(), literal * 2)
        self.assertFalse((self.root / "SHOULD_NOT_EXIST").exists())
        self.assertFalse((self.root / "ALSO_NOT").exists())

    def test_new_plan_never_skips_new_inputs(self):
        job = self.job("same", "pass")
        p1 = write_plan([job], self.root, mode="local")
        p2 = write_plan([job], self.root, mode="local")
        self.assertNotEqual(p1, p2)

    def test_hpc_no_login_execution(self):
        plan = write_plan([self.job("heavy", "from pathlib import Path; Path('wrong').touch()")], self.root, mode="hpc")
        env = {k: v for k, v in os.environ.items() if not k.startswith("LSB_")}
        p = subprocess.run(["bash", str(plan.parent / "commands/heavy.cmd")], capture_output=True, env=env)
        self.assertNotEqual(p.returncode, 0)
        self.assertFalse((self.root / "wrong").exists())

    def test_gpu_and_memory_header(self):
        job = Job("gpu", ["true"], resources={"queue": "gpu", "cores": 8, "memory_gib": 64, "gpus": 1, "hours": 120})
        plan = write_plan([job], self.root)
        cmd = (plan.parent / "commands/gpu.cmd").read_text()
        self.assertIn("#BSUB -M 64GB", cmd)
        self.assertIn("rusage[mem=64GB/host]", cmd)
        self.assertIn("#BSUB -W 120:00", cmd)
        self.assertIn('num=1:mode=exclusive_process', cmd)

    def test_resource_and_label_injection_rejected(self):
        for r in ({"queue": "ser\n#BSUB -q hacked"}, {"cores": 1.5}, {"memory_gib": float('nan')}, {"gpus": -1}):
            with self.assertRaises(ValueError):
                resources(r)
        with self.assertRaises(ValueError):
            write_plan([self.job("../x", "pass")], self.root)

    def test_mock_lsf_submission_dependencies_and_receipts(self):
        jobs = [self.job("a", "pass"), self.job("b", "pass"), self.job("merge", "pass", ["a", "b"])]
        plan = write_plan(jobs, self.root)
        calls = []
        def fake(argv, **kwargs):
            calls.append(argv)
            return subprocess.CompletedProcess(argv, 0, f"Job <{100+len(calls)}> is submitted to queue <ser>.\n", "")
        with patch("workflow_jobs.subprocess.run", side_effect=fake):
            run_plan(plan, submit="hpc", jobs=2)
            run_plan(plan, submit="hpc", jobs=2)
        self.assertEqual(len(calls), 3)
        self.assertEqual(calls[0], ["bsub"])
        self.assertEqual(calls[1], ["bsub"])
        self.assertIn("done(101)", calls[2][2])
        self.assertIn("done(102)", calls[2][2])
        self.assertIn("ended(101)", calls[2][2])
        journal = json.loads((plan.parent / "submission.json").read_text())
        self.assertEqual(journal["jobs"]["merge"]["job_id"], "103")

    def test_unknown_submission_not_retried(self):
        plan = write_plan([self.job("x", "pass")], self.root)
        with patch("workflow_jobs.subprocess.run", return_value=subprocess.CompletedProcess([], 1, "", "transport failed")) as call:
            with self.assertRaises(RuntimeError):
                run_plan(plan)
            with self.assertRaises(RuntimeError):
                run_plan(plan)
            self.assertEqual(call.call_count, 1)

    def test_preview_never_calls_bsub(self):
        plan = write_plan([self.job("x", "pass")], self.root)
        with patch("workflow_jobs.subprocess.run", side_effect=AssertionError("should not execute")):
            self.assertEqual(run_plan(plan, dry_run=True)["status"], "PLANNED")

    def test_portable_directory_alias(self):
        code = 'source "$1"; script_directory_args /work/sph-huangj/scripts/grid/2.csx.sh --dir-analysis-root /result; printf "%s\\n" "${SCRIPT_DIRECTORY_ARGS[@]}"'
        p = subprocess.run(["bash", "-c", code, "bash", str(SHARED / "cli_dirs.sh")], capture_output=True, text=True, check=True)
        self.assertEqual(p.stdout.splitlines(), ["--analysis-root", "/result"])

    def test_path_table_mapping(self):
        with patch.dict(os.environ, {"submit": "hpc", "DIR_WORK": "/work/me", "DIR_DATA": "/data/me", "DIR_SCRATCH": "/scratch/date/me", "ANALYSIS_ROOT": "/scratch/date/me/analysis"}):
            self.assertEqual(relocate_defaults(["/mnt/f/gen/a", "/mnt/d/analysis/grid"]), ["/data/me/gen/a", "/scratch/date/me/analysis/grid"])
            self.assertEqual(relocate_defaults("/custom/path"), "/custom/path")

    def test_scratch_date_relocates_inherited_defaults(self):
        env = {**os.environ, "submit": "hpc", "DATE_SCRATCH": "2026-10-31",
               "DIR_SCRATCH": "/scratch/2026-10-31/sph-huangj",
               "GRID_OUTPUT_ROOT": "/scratch/2026-10-31/sph-huangj/analysis/grid",
               "SCRIPT_LOG_DIR": "/my/custom/logs"}
        code = 'source "$1"; workflow_init grid --date-scratch 2026-11-02; printf "%s\\n" "$DIR_SCRATCH" "$GRID_OUTPUT_ROOT" "$SCRIPT_LOG_DIR"'
        p = subprocess.run(["bash", "-c", code, "bash", str(SHARED/"workflow_env.sh")],
                           env=env, capture_output=True, text=True, check=True)
        self.assertEqual(p.stdout.splitlines(), ["/scratch/2026-11-02/sph-huangj",
            "/scratch/2026-11-02/sph-huangj/analysis/grid", "/my/custom/logs"])

    def test_scratch_date_preserves_custom_directory(self):
        env = {**os.environ, "submit": "hpc", "DATE_SCRATCH": "2026-10-31", "DIR_SCRATCH": "/my/project/scratch"}
        code = 'source "$1"; workflow_init grid --date-scratch 2026-11-02; printf "%s\\n" "$DATE_SCRATCH" "$DIR_SCRATCH"'
        p = subprocess.run(["bash", "-c", code, "bash", str(SHARED/"workflow_env.sh")],
                           env=env, capture_output=True, text=True, check=True)
        self.assertEqual(p.stdout.splitlines(), ["2026-11-02", "/my/project/scratch"])


if __name__ == "__main__":
    unittest.main()
