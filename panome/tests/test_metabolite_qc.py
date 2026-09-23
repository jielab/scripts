"""Regression coverage for raw negative abundances, including held-out rows."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
from io_data import map_metabolites, numeric
from preprocess import MolecularPreprocessor


class MetaboliteQCTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.mapping = Path(tmp.name) / "met.lst"
        self.mapping.write_text(
            "data_field\tmet_name\n"
            "p1\tNon_HDL_C\np2\tB\np3\tC\np2/p3\tRatio\np4\tAbsent\n")
        self.frame = pd.DataFrame({
            "eid": ["a", "b", "c", "d", "held_out"],
            "p1_i0": [1., 2., 3., 4., -.073459],
            "p2_i0": [2., 3., 4., 5., 6.],
            "p3_i0": [1., 3., 2., 4., 0.],
            "p1_i1": [-10.] * 5})

    def test_held_out_negative_is_missing_before_qc_and_imputation(self):
        source = self.frame.copy(deep=True)
        mapped, audit = map_metabolites(self.frame, self.mapping, "eid", nonnegative=True)
        audit = audit.set_index("feature")
        self.assertEqual(audit.loc["Non_HDL_C", "negative_to_missing"], 1)
        self.assertEqual(audit.loc["Ratio", "nonfinite_to_missing"], 1)
        self.assertEqual(audit.loc["Absent", "status"], "missing_source")
        raw = numeric(mapped, list(mapped.columns[1:]), "test")
        people = pd.DataFrame(index=range(len(raw)))
        prep = MolecularPreprocessor(transform="log1p").fit(raw[:4], people.iloc[:4])
        missing = (~np.isfinite(prep.scale_transform(raw[:, prep.keep]))).mean(1)
        self.assertEqual(missing[-1], .5)
        x, observed = prep.transform(raw, people)
        self.assertTrue(np.isfinite(x).all())
        self.assertFalse(observed[-1, 0])
        self.assertFalse(observed[-1, 3])
        expected = (prep.median[0] - prep.mean[0]) / prep.scale[0]
        self.assertAlmostEqual(x[-1, 0], expected, places=6)
        pd.testing.assert_frame_equal(self.frame, source)

    def test_signed_input_is_preserved_without_nonnegative_policy(self):
        mapped, audit = map_metabolites(self.frame, self.mapping, "eid")
        self.assertAlmostEqual(mapped.Non_HDL_C.iloc[-1], -.073459, places=6)
        self.assertEqual(audit.negative_to_missing.sum(), 0)
        raw = numeric(mapped, list(mapped.columns[1:]), "test")
        np.testing.assert_allclose(MolecularPreprocessor(transform="none").scale_transform(raw), raw)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            MolecularPreprocessor(transform="log1p").scale_transform(raw)


if __name__ == "__main__":
    unittest.main()
