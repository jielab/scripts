"""Native R checks for the reviewed PCA QC and legacy evaluation guard."""
from pathlib import Path
import importlib.util
import os
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
import openpyxl
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


class NativeRReview(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rscript = os.environ.get("GRID_TEST_RSCRIPT") or shutil.which("Rscript") or str(
            Path.home() / "miniforge3/envs/grid/bin/Rscript")
        if not Path(cls.rscript).is_file():
            raise unittest.SkipTest("Rscript required for native R checks")

    def test_grid_rds_table_with_raw_attributes(self):
        spec = importlib.util.spec_from_file_location("grid_data_rds_review", ROOT / "f/grid.data.py")
        data = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(data)
        with tempfile.TemporaryDirectory(prefix="grid-raw-rds-") as temporary:
            path = Path(temporary) / "scores.rds"
            script = Path(temporary) / "fixture.R"
            script.write_text('''x <- data.frame(eid=c("001","002","003"), score=c(1/3,NA,-2), extra=1:3, family_id=c("01","1","02"))
attr(x, "provenance_payload") <- as.raw(c(0,127,255))
saveRDS(x, commandArgs(TRUE)[1])
''')
            subprocess.run([self.rscript, "--vanilla", str(script), str(path)], check=True,
                           capture_output=True, text=True)
            frame = data.read_table(path)
            self.assertEqual(frame.eid.tolist(), ["001", "002", "003"])
            self.assertEqual(frame.columns.tolist(), ["eid", "score", "extra", "family_id"])
            self.assertEqual(frame.family_id.tolist(), ["01", "1", "02"])
            np.testing.assert_allclose(frame.score, [1/3, np.nan, -2], atol=1e-15)
            selected = data.read_table(path, ["eid", "score", "absent"])
            pd.testing.assert_frame_equal(selected, frame[["eid", "score"]])

    def test_pca_qc_reads_separate_score_directory(self):
        with tempfile.TemporaryDirectory(prefix="grid-r-pca-") as temporary:
            root = Path(temporary)
            raw = root / "scores"
            raw.mkdir()
            rng = np.random.default_rng(901)
            pcs = [f"PC{i}" for i in range(1, 6)]
            pca = pd.DataFrame(rng.normal(size=(40, 5)), columns=pcs)
            pca.insert(0, "IID", [f"synthetic{i}" for i in range(len(pca))])
            pca.to_csv(root / "pca.tsv", sep="\t", index=False)
            med = pd.DataFrame(rng.normal(size=(4, 5)), columns=pcs)
            med.insert(0, "POP", ["AFR", "EAS", "EUR", "SAS"])
            med.to_csv(root / "med.tsv", sep="\t", index=False)
            for chrom in range(1, 23):
                (raw / f"chr{chrom}.sscore.vars").write_text(
                    "".join(f"rs{chrom}_{j}\n" for j in range(chrom)))
            run = subprocess.run([self.rscript, str(ROOT / "f/0.pca.R"), "projection",
                "--pca", str(root / "pca.tsv"), "--med", str(root / "med.tsv"),
                "--outdir", str(root / "qc"), "--n-pc", "5", "--distance-pcs", "5",
                "--score-dir", str(raw)], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            book = openpyxl.load_workbook(root / "qc/pca_qc.xlsx", read_only=True, data_only=True)
            try:
                sheet = book["Variant_match"]
                # openxlsx may emit dimension="A1" despite complete cell XML.
                # Stream the actual cells instead of trusting that hint.
                sheet.reset_dimensions()
                rows = list(sheet.values)
                self.assertEqual(rows[0], ("chr", "matched_pca_variants"))
                self.assertEqual(rows[1:], [(i, i) for i in range(1, 23)])
            finally:
                book.close()

    def test_r_yeval_rejects_learned_scores_before_input_access(self):
        with tempfile.TemporaryDirectory(prefix="grid-r-yeval-") as temporary:
            output = Path(temporary) / "must-not-be-created"
            run = subprocess.run([self.rscript, str(ROOT / "f/Yeval.R"),
                "--trait", "height", "--type", "ct", "--grid-file", "/does/not/exist",
                "--run-dir", str(output)], capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("Legacy --grid-file is disabled", run.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
