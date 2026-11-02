"""Exercise the public shell/dispatcher boundary without launching large jobs."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
import herv_pipeline


class ResearchCLI(unittest.TestCase):
    def invoke(self, *args):
        return subprocess.run(["bash", str(ROOT / "allele3.sh"), "research", *args,
                               "--python", sys.executable], cwd=ROOT, text=True,
                              capture_output=True, timeout=30)

    def test_each_public_help_reaches_its_parser(self):
        for module in ("data", "allofus", "mfi", "mapped-targets", "panel", "balance", "space", "domains", "cohort", "archaic-alleles"):
            with self.subTest(module=module):
                result = self.invoke(module, "--help")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("usage:", result.stdout.lower())

    def test_default_output_is_forwarded_and_dry_run_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {"ALLELE_ANALYSIS_ROOT": "/tmp/allele3-test-default"}):
            result = self.invoke("balance", "--genotypes", str(Path(tmp) / "absent.tsv"), "--dry-run")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--out", result.stdout)
            self.assertIn("/tmp/allele3-test-default/research/balance", result.stdout)
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_cohort_defaults_keep_plan_and_summary_separate(self):
        for mode in ("plan", "summarize"):
            with self.subTest(mode=mode), patch.dict(os.environ, {"ALLELE_ANALYSIS_ROOT": "/tmp/allele3-test-default"}):
                result = self.invoke("cohort", mode, "--dry-run")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("/research/cohort/" + mode, result.stdout)

    def test_selected_references_reach_both_summary_and_archaic(self):
        args = herv_pipeline.parser().parse_args(["--archaic-references", "Altai,Vindija"])
        with tempfile.TemporaryDirectory() as tmp:
            args.dir_analysis_root = tmp
            args.dir_out = str(Path(tmp)/"result")
            args.dir_modern = str(Path(tmp)/"scan")
            paths = herv_pipeline.resolve(args)
            # Input inventories are independent of command construction; do not
            # require or forge an actual large scan/GU corpus for this check.
            with patch.object(herv_pipeline, "modern_files", return_value=[]), \
                    patch("herv_archaic.discover", return_value=[]) as discover, \
                    patch("herv_archaic.modern_inputs", return_value=[]):
                stages = herv_pipeline.plan(args, paths, ["summary", "archaic"])
            by_name = {s.name: s for s in stages}
            for name, flag in (("summary", "--archaic-references"), ("archaic", "--references")):
                command = by_name[name].command
                self.assertEqual(command[command.index(flag)+1], "Altai,Vindija")
            self.assertEqual(discover.call_args.args[-1], "Altai,Vindija")


if __name__ == "__main__":
    unittest.main()
