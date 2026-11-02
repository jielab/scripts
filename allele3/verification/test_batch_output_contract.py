"""Regression checks for output selection before any scheduler submission."""
import os
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
import batch_workflows as workflows


class BatchOutputContract(unittest.TestCase):
    def setUp(self):
        self.env = patch.dict(os.environ, {
            "DIR_DATA": "/data/test", "DIR_WORK": "/work/test",
            "DIR_SCRATCH": "/scratch/test", "ALLELE_ANALYSIS_ROOT": "/scratch/test/analysis/allele3",
        }, clear=True)
        self.env.start()
        self.addCleanup(self.env.stop)

    def plan(self, entry, argv, mode="hpc"):
        return workflows.single_plan(ROOT / entry, list(argv), mode, sys.executable)

    def test_fresh_scan_never_invents_permission_to_replace_default_store(self):
        for argv in ([], ["run"], ["run", "--replace", "FALSE"]):
            with self.subTest(argv=argv):
                with self.assertRaisesRegex(ValueError, "requires --out"):
                    self.plan("3.scan.sh", argv)
        plan, _, out = self.plan("3.scan.sh", ["run", "--out", "/scratch/test/new-scan"])
        self.assertEqual(str(out), "/scratch/test/new-scan")
        self.assertIn("--out", plan["jobs"][0]["argv"])

    def test_fresh_pca_preserves_explicit_output_guard_and_environment_override(self):
        for argv in ([], ["run"]):
            with self.subTest(argv=argv):
                with self.assertRaisesRegex(ValueError, "requires --out"):
                    self.plan("1.pca3D.sh", argv)
        with patch.dict(os.environ, {"PCA3D_OUT": "/scratch/test/selected-pca"}):
            plan, _, out = self.plan("1.pca3D.sh", ["run"])
            self.assertEqual(str(out), "/scratch/test/selected-pca")
            argv = plan["jobs"][0]["argv"]
            self.assertEqual(argv[argv.index("--out") + 1], str(out))

    def test_environment_outputs_are_honored_and_cli_wins(self):
        for entry, env_key, argv in (
                ("2.herv.sh", "HERV_OUT", []),
                ("2.herv.sh", "ALLELE_HERV_OUT", []),
                ("allele3.sh", "ALLELE_HERV_FOCUS_OUT", ["primary"]),
                ("1.pca3D.sh", "PCA3D_OUT", ["render"])):
            with self.subTest(entry=entry, env_key=env_key), patch.dict(os.environ, {env_key: "/scratch/test/configured"}):
                _, _, configured = self.plan(entry, argv)
                _, _, explicit = self.plan(entry, [*argv, "--out", "/scratch/test/explicit"])
                self.assertEqual(str(configured), "/scratch/test/configured")
                self.assertEqual(str(explicit), "/scratch/test/explicit")

    def test_default_herv_outputs_keep_dataset_build_and_sample_identity(self):
        for mode, prefix in (("hpc", "/scratch/test/allele3/herv/1kg"),
                             ("local", "/data/test/gen/1kg")):
            outputs = []
            for build in ("37", "38"):
                _, _, out = self.plan("2.herv.sh", ["--build", build, "--sample", "HG00096"], mode)
                expected = f"{prefix}/{build}/HG00096" if mode == "hpc" else f"{prefix}/{build}/herv/HG00096"
                self.assertEqual(str(out), expected)
                outputs.append(out)
            self.assertNotEqual(*outputs)

    def test_existing_scan_actions_still_accept_configured_store(self):
        with patch.dict(os.environ, {"ALLELE_SCAN_STORE": "/data/test/existing-scan"}):
            for action in ("plan", "verify", "adopt"):
                with self.subTest(action=action):
                    _, _, out = self.plan("3.scan.sh", [action])
                    self.assertEqual(str(out), "/data/test/existing-scan")

    def test_cohort_plan_and_summary_keep_separate_default_directories(self):
        outputs = []
        for action in ("plan", "summarize"):
            _, _, out = self.plan("allele3.sh", ["research", "cohort", action])
            self.assertEqual(str(out), "/scratch/test/analysis/allele3/research/cohort/" + action)
            outputs.append(out)
        self.assertNotEqual(*outputs)


if __name__ == "__main__":
    unittest.main()
