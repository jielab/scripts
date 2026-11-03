"""Entrypoint routing only; never read or rebuild participant data."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class RscriptEntry(unittest.TestCase):
    def test_pipeline_honors_selected_interpreter_and_steps(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "rap").mkdir()
            fake = root / "custom Rscript"
            fake.write_text("#!/usr/bin/env python3\nimport os,sys,json\n"
                            "print(json.dumps({'args':sys.argv[1:],"
                            "'steps':os.environ['RUN_STEPS'],"
                            "'root':os.environ['PHEDIR']}))\n")
            fake.chmod(0o755)
            result = subprocess.run([
                "bash", str(ROOT / "phe.sh"), "--phenotype-dir", tmp,
                "--steps", "data_qc,vip_phe,merge,phe4gwas",
                "--rscript", str(fake)], text=True, capture_output=True,
                env={**os.environ, "SCRIPT_CONSOLE_ACTIVE": "1"}, check=True)
            record = json.loads(result.stdout.strip().splitlines()[-1])
            self.assertEqual(record["args"], [str(ROOT / "f/phe.R")])
            self.assertEqual(record["steps"], "data_qc,vip_phe,merge,phe4gwas")
            self.assertEqual(record["root"], tmp)

    def test_interpreter_failure_is_not_hidden(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "rap").mkdir()
            fake = root / "Rscript"
            fake.write_text("#!/bin/sh\nexit 17\n")
            fake.chmod(0o755)
            result = subprocess.run([
                "bash", str(ROOT / "phe.sh"), "--phenotype-dir", tmp,
                "--steps", "data_qc"], capture_output=True,
                env={**os.environ, "SCRIPT_CONSOLE_ACTIVE": "1",
                     "PATH": tmp + os.pathsep + os.environ["PATH"],
                     "RSCRIPT": "/missing/old/Windows/Rscript.exe"})
            self.assertEqual(result.returncode, 17)


if __name__ == "__main__":
    unittest.main()
