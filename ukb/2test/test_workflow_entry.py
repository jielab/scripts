"""Plan/dispatch regressions with fake LSF; no real cohort or scheduler access."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class WorkflowEntry(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="ukb workflow ")
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.env = {**os.environ, "UKB_PYTHON": sys.executable,
                    "DIR_DAT": str(self.work), "ANALYSIS_ROOT": str(self.work / "analysis")}

    def run_script(self, name, *args, ok=True, env=None):
        result = subprocess.run(["bash", str(ROOT / name), *map(str, args)],
                                env=env or self.env, cwd=self.work, capture_output=True, text=True)
        if ok:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_hpc_plan_freezes_paths_resources_and_scientific_options(self):
        result = self.run_script("phe.sh", "--submit", "hpc", "--steps", "phe4gwas",
            "--phenotype-dir", self.work / "cohort's files", "--write-cmd",
            env={**self.env, "ADULT_ENABLE": "FALSE", "SRD_YEAR_CAP": "2026"})
        cmd = next(self.work.rglob("*.cmd"))
        text = cmd.read_text()
        self.assertIn("#BSUB -n 40\n", text)
        self.assertIn("rusage[mem=160GB/host]", text)
        self.assertIn("export submit=hpc", text)
        self.assertIn("export SRD_YEAR_CAP=2026", text)
        self.assertIn("--execute --steps phe4gwas", text)
        self.assertIn("LSB_JOBID", text)
        self.assertIn(str(cmd), result.stdout)
        subprocess.run(["bash", "-n", cmd], check=True)

    def test_all_genotype_modules_plan_without_running_tools(self):
        for module, args in [
            ("pgs", ["--input", "not-present.tsv", "--label", "score", "--pgs-threads", "12"]),
            ("snp", ["--data", "typ", "--snp-file", "not-present.bed", "--snp-file-type", "bed", "--label", "test"]),
            ("apoe", []),
            ("liftGen", ["--dir-in", "not-present", "--chain-file", "not-present.chain"]),
        ]:
            with self.subTest(module=module):
                out = self.work / module
                self.run_script("gen.sh", module, "--submit", "hpc", "--dir-out", out, "--write-cmd", *args)
                commands = list(out.rglob("*.cmd"))
                self.assertEqual(len(commands), 1)
                self.assertFalse(list(out.rglob("*.worker.sh")))
                self.assertIn("--execute", commands[0].read_text())

    def test_fake_bsub_receives_one_cmd_and_records_job_id(self):
        fakebin = self.work / "bin"
        fakebin.mkdir()
        fake = fakebin / "bsub"
        fake.write_text("#!/bin/sh\ncat > \"$BSUB_CAPTURE\"\nprintf 'Job <123456> is submitted to queue <33>.\\n'\n")
        fake.chmod(0o755)
        capture = self.work / "received.cmd"
        result = self.run_script("phe.sh", "--submit", "hpc", "--steps", "phe4gwas", "--memory-gb", "32",
            env={**self.env, "PATH": str(fakebin) + os.pathsep + os.environ["PATH"], "BSUB_CAPTURE": str(capture)})
        self.assertIn("123456", result.stdout)
        self.assertIn("#BSUB -n 10\n", capture.read_text())
        journals = list((self.work / "analysis").rglob("submission-*.jsonl"))
        self.assertEqual(len(journals), 1)
        records = [json.loads(line) for line in journals[0].read_text().splitlines()]
        self.assertEqual(records[-1]["job_id"], "123456")

    def test_local_worker_keeps_roots_arguments_and_failure_code(self):
        raw = self.work / "data with spaces" / "rap"
        raw.mkdir(parents=True)
        fake = self.work / "custom Rscript"
        record = self.work / "runtime.json"
        fake.write_text("#!/usr/bin/env python3\nimport os,json\n"
            "from pathlib import Path\nPath(os.environ['RECORD']).write_text(json.dumps("
            "{k:os.environ[k] for k in ['submit','DIR_DAT','PHEDIR','RUN_STEPS']}))\n")
        fake.chmod(0o755)
        self.run_script("phe.sh", "--submit", "local", "--phenotype-dir", raw.parent,
            "--steps", "phe4gwas", "--rscript", fake,
            env={**self.env, "submit": "hpc", "RECORD": str(record)})
        result = json.loads(record.read_text())
        self.assertEqual(result["submit"], "local")
        self.assertEqual(result["PHEDIR"], str(raw.parent))
        self.assertEqual(result["DIR_DAT"], str(self.work))
        self.assertEqual(result["RUN_STEPS"], "phe4gwas")

    def test_invalid_options_and_login_worker_fail_before_dispatch(self):
        self.run_script("phe.sh", "--submit", "hpc", "--steps", "wrong", "--write-cmd", ok=False)
        self.run_script("gen.sh", "apoe", "--submit", "hpc", "--memory-gb", "999", "--write-cmd", ok=False)
        self.run_script("phe.sh", "--submit", "hpc", "--execute", "--steps", "phe4gwas",
                        env={k: v for k, v in self.env.items() if k != "LSB_JOBID"}, ok=False)
        self.assertFalse(list(self.work.rglob("*.cmd")))


if __name__ == "__main__":
    unittest.main()
