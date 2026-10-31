"""Run with: python3 -m unittest discover -s verification -p 'test_*.py'."""
import contextlib
import gzip
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
import herv_candidates
import scan

spec = importlib.util.spec_from_file_location("pipeline_console", ROOT / "runtime/pipeline_console.py")
pipeline_console = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipeline_console)


class PipelineRegressions(unittest.TestCase):
    def test_empty_final_annotation_field_and_cache_reuse(self):
        row = "22\tchr13\t114292321\t114298487\tL1P4#LINE/L1\t2610\t-\t114292338\t114292382\t0\t3\t17,44,6105\t-1,17,-1\t"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "joined.gz"
            with gzip.open(source, "wt") as f:
                f.write("# annotation\n\n" + row + "\n" + row + "record2\n")
            cache, _ = herv_candidates.split_annotation(source, {}, root / "cache")
            with gzip.open(cache / "chr13.tsv.gz", "rt") as f:
                self.assertEqual(f.read(), row + "\n" + row + "record2\n")
            self.assertEqual(len(scan.parse_joined_rmsk(cache / "chr13.tsv.gz")), 2)
            before = (cache / "chr13.tsv.gz").stat().st_mtime_ns
            herv_candidates.split_annotation(source, {}, cache)
            self.assertEqual(before, (cache / "chr13.tsv.gz").stat().st_mtime_ns)

    def test_malformed_annotation_reports_location_without_publishing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "bad.tsv"
            source.write_text("# comment\n22\tchr13\t114292321\n")
            with self.assertRaisesRegex(ValueError, r"bad.tsv:2: expected 14.*got 3"):
                herv_candidates.split_annotation(source, {}, root / "cache")
            self.assertFalse((root / "cache/manifest.json").exists())

    def test_compact_screen_preserves_full_log_and_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "full.log"
            raw = ("[stage] alleles: /long/python command --many-options\n"
                   + "[herv-candidates] read shard\n" * 100
                   + '{"status": "complete", "failed": 0}\n'
                   + "herv-candidates: expected 14 columns\n"
                   + "[stage] ld: /long/python ld\n"
                   + "[herv-focus] failed: /out/pipeline_status.json\n")
            old = pipeline_console.console.ScreenLog
            pipeline_console.console.ScreenLog = pipeline_console.PipelineScreenLog
            output = io.StringIO()
            try:
                with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
                    code = pipeline_console.console.run(
                        [sys.executable, "-c", "import sys; sys.stdout.write(" + repr(raw) + "); sys.exit(2)"],
                        log, interval=300)
            finally:
                pipeline_console.console.ScreenLog = old
            self.assertEqual(code, 2)
            self.assertEqual(log.read_text(), raw)
            screen = output.getvalue()
            self.assertEqual(screen.count("read shard"), 1)
            self.assertNotIn("/long/python", screen)
            self.assertNotIn('"status"', screen)
            self.assertIn("[stage] alleles: running", screen)
            self.assertIn("herv-candidates: expected 14 columns", screen)
            self.assertLess(len(screen.splitlines()), 12)


if __name__ == "__main__":
    unittest.main()
