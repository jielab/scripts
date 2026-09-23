"""Launcher regression checks without loading data or starting model training."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.calls = Path(self.tmp.name) / "calls.jsonl"
        stub = Path(self.tmp.name) / "python"
        stub.write_text('''#!/usr/bin/env python3
import json, os, runpy, sys
entry = runpy.run_path(sys.argv[1])
a = entry["configure"](entry["parser"]().parse_args(sys.argv[2:]))
with open(os.environ["PANOME_TEST_CALLS"], "a") as out:
    out.write(json.dumps(vars(a)) + "\\n")
if os.environ.get("PANOME_TEST_FAIL") == a.trait + "/" + a.biom:
    sys.exit(7)
''')
        stub.chmod(0o755)
        self.env = dict(os.environ, PANOME_PYTHON=str(stub),
                        PANOME_TEST_CALLS=str(self.calls), PANOME_TEST_FAIL="")

    def launch(self, *args):
        return subprocess.run(["bash", str(ROOT / "panome.sh"), *args],
                              env=self.env, capture_output=True, text=True)

    def records(self):
        return [json.loads(line) for line in self.calls.read_text().splitlines()]

    def test_outcome_layer_product_and_defaults(self):
        result = self.launch("--Y", "cvd_cad,ra", "--biom", "prot,met")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = self.records()
        self.assertEqual([(a["trait"], a["biom"]) for a in rows],
                         [("cvd_cad", "prot"), ("cvd_cad", "met"),
                          ("ra", "prot"), ("ra", "met")])
        for a in rows:
            self.assertEqual(a["diagnosis_col"], "fod_icd10_" + a["trait"])
            self.assertEqual(a["residualize"], a["biom"] + ".plate")
            self.assertEqual(a["transform"], "log1p" if a["biom"] == "met" else "none")
            self.assertEqual(a["device"], "cuda")

    def test_equals_syntax_and_explicit_overrides(self):
        result = self.launch("--Y=ra,cvd_cad", "--biom=met", "--device", "cpu",
                             "--omics-file=/tmp/shared matrix.csv", "--run-name", "check")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([a["trait"] for a in self.records()], ["ra", "cvd_cad"])
        for a in self.records():
            self.assertEqual(a["omics_file"], "/tmp/shared matrix.csv")
            self.assertEqual(a["device"], "cpu")
            self.assertEqual(a["run_name"], "check")

    def test_single_run_and_last_option_wins(self):
        result = self.launch("--Y", "cvd_cad,ra", "--Y=ra", "--biom=met,prot",
                             "--biom", "prot", "--diagnosis-col", "custom_date",
                             "--run-dir", "/tmp/custom run")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = self.records()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["trait"], "ra")
        self.assertEqual(rows[0]["diagnosis_col"], "custom_date")
        self.assertEqual(rows[0]["run_dir"], "/tmp/custom run")

    def test_invalid_batches_fail_before_launch(self):
        cases = [("--Y", value) for value in ["", ",ra", "ra,", "ra,,cad", "ra,ra"]]
        cases += [("--biom", value) for value in ["prot,prot", "prot,", "other"]]
        cases += [("--Y=ra,cad", option) for option in
                  ["--run-dir=/tmp/run", "--output=/tmp/out.csv", "--diagnosis-col=date"]]
        cases += [("--biom=prot,met", option) for option in
                  ["--run-dir=/tmp/run", "--output=/tmp/out.csv", "--omics-file=/tmp/x.csv"]]
        for args in cases:
            with self.subTest(args=args):
                self.assertNotEqual(self.launch(*args).returncode, 0)
                self.assertFalse(self.calls.exists())

    def test_first_failure_stops_remaining_runs(self):
        self.env["PANOME_TEST_FAIL"] = "cvd_cad/met"
        result = self.launch("--Y=cvd_cad,ra", "--biom=prot,met")
        self.assertEqual(result.returncode, 7)
        self.assertEqual([(a["trait"], a["biom"]) for a in self.records()],
                         [("cvd_cad", "prot"), ("cvd_cad", "met")])


if __name__ == "__main__":
    unittest.main()
