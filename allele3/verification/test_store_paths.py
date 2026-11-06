"""Persistent stores and primary-stage isolation from genotype uploads."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
import herv_pipeline as pipeline


class StorePaths(unittest.TestCase):
    def test_primary_needs_only_completed_scan_and_merged_counters(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            modern = root / "stores/gnomad_scan"
            (modern / "summary").mkdir(parents=True)
            (modern / "run.json").write_text(json.dumps({
                "status": "complete", "signature": "fixture", "shards": [str(modern / "shards/not-needed")],
                "options": {"fasta": str(root / "gen-still-uploading/reference.fa")},
            }))
            counters = [modern / "summary" / f"architecture.{kind}_strata.tsv.gz"
                        for kind in ("global", "exposure")]
            for path in counters:
                path.write_bytes(b"metadata-only fixture")
            args = pipeline.parser().parse_args([
                "primary", "--dir-analysis-root", str(root), "--scan-dir", str(modern),
                "--dir-ancient", str(root / "gen-still-uploading/ancient"),
                "--fasta38", str(root / "gen-still-uploading/reference.fa"),
            ])
            with patch("herv_archaic.discover", side_effect=AssertionError("primary discovered gen inputs")):
                steps = pipeline.plan(args, pipeline.resolve(args), ["primary"])
            step = next(s for s in steps if s.name == "primary")
            self.assertEqual(set(step.inputs), {modern / "run.json", *counters})
            self.assertTrue(all(p.is_file() for p in step.inputs))

    def test_archaic_store_is_separate_from_input_and_honors_override(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {"ALLELE_ARCHAIC_CACHE": ""}):
            root = Path(tmp)
            for explicit in (None, root / "custom-cache"):
                argv = ["archaic", "--ibdmix", "TRUE", "--dir-analysis-root", str(root / "analysis"),
                        "--scan-dir", str(root / "annot/source-scan")]
                if explicit:
                    argv += ["--archaic-cache", str(explicit)]
                args = pipeline.parser().parse_args(argv)
                with patch("herv_archaic.discover", return_value=[]), patch("herv_archaic.modern_inputs", return_value=[]):
                    steps = pipeline.plan(args, pipeline.resolve(args), ["archaic"])
                step = next(s for s in steps if s.name == "archaic")
                actual = step.command[step.command.index("--cache-dir") + 1]
                self.assertEqual(actual, str(explicit or root / "analysis/stores/archaic_overlap"))


if __name__ == "__main__":
    unittest.main()
