from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
import data
import prsformer


class PhenotypeTests(unittest.TestCase):
    def read(self, frame, traits, **options):
        args = SimpleNamespace(**dict({"covariates": "none"}, **options))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "phenotypes.tsv"
            frame.to_csv(path, sep="\t", index=False)
            before = path.read_bytes()
            result = data.read_phenotypes(path, traits, args)
            self.assertEqual(path.read_bytes(), before)
        return result

    def test_ldl_mapping_preserves_ids_values_and_missing_labels(self):
        frame = pd.DataFrame({"eid": ["003", "001", "002"], "bb_LDL": [2.847, np.nan, 4.1]})
        result = self.read(frame, ["ldl"])
        self.assertEqual(result.eid.tolist(), frame.eid.tolist())
        np.testing.assert_allclose(result.ldl, frame.bb_LDL, equal_nan=True)
        self.assertEqual(result.attrs["phenotype_definitions"]["ldl"]["source_column"], "bb_LDL")

    def test_existing_or_explicit_ldl_columns_take_precedence(self):
        frame = pd.DataFrame({"eid": ["001", "002"], "ldl": [8., np.nan],
                              "bb_LDL": [3., 4.], "custom": [5., 6.]})
        np.testing.assert_allclose(self.read(frame, ["ldl"]).ldl, [8., np.nan], equal_nan=True)
        np.testing.assert_array_equal(self.read(frame, ["ldl"], ldl_col="custom").ldl, [5., 6.])
        with self.assertRaisesRegex(ValueError, "Missing phenotype column: missing"):
            self.read(frame, ["ldl"], ldl_col="missing")

    def test_raw_t2dm_uses_baseline_not_incident_status(self):
        # Prior/same-day cases, later/after-censor cases, controls, undated
        # disease evidence, missing baseline and invalid follow-up.
        frame = pd.DataFrame({
            "eid": [f"{i:03}" for i in range(10)],
            "date_attend": ["2010-01-01"] * 7 + [None, "2010-01-01", "2024-01-01"],
            "date_lost": [None] * 8 + ["2009-01-01", None],
            "date_death": [None] * 10,
            "fod_icd10_t2dm": ["2009-01-01", "2010-01-01", "2011-01-01", "2025-01-01",
                               None, None, None, "2009-01-01", "2008-01-01", None],
            "cnt_icd10_t2dm": [1, 1, 1, 1, 0, 1, np.nan, 1, 1, 0],
        })
        result = self.read(frame, ["t2dm"])
        np.testing.assert_allclose(result.t2dm, [1, 1, 0, 0, 0, np.nan, 0, np.nan, np.nan, np.nan], equal_nan=True)
        definition = result.attrs["phenotype_definitions"]["t2dm"]["derived_from"]
        self.assertEqual(definition["administrative_end_date"], "2023-04-01")
        self.assertEqual(Path(definition["helper"]["path"]), data.PHENOTYPE_HELPER)
        # Reading already-derived endpoints produces identical labels without R.
        with patch.object(data, "derive_t2dm_endpoints", side_effect=AssertionError("must reuse")):
            existing = self.read(result, ["t2dm"])
        np.testing.assert_allclose(existing.t2dm, result.t2dm, equal_nan=True)

    def test_missing_partial_and_explicit_t2dm_endpoints(self):
        frame = pd.DataFrame({"eid": ["001", "002"], "baseline": [0., 1.]})
        with self.assertRaisesRegex(ValueError, "Run UKB phe4gwas"):
            self.read(frame, ["t2dm"])
        with self.assertRaisesRegex(ValueError, "Incomplete T2DM endpoint pair"):
            self.read(frame.assign(**{"t2dm.Yr2e": 0.}), ["t2dm"])
        with patch.object(data, "derive_t2dm_endpoints", side_effect=AssertionError("explicit endpoint")):
            result = self.read(frame, ["t2dm"], t2dm_col="baseline")
        np.testing.assert_array_equal(result.t2dm, [0., 1.])
        with self.assertRaisesRegex(ValueError, "baseline 0/1 outcome"):
            self.read(frame.assign(baseline=[0., 2.]), ["t2dm"], t2dm_col="baseline")

    def test_outcome_alias_and_derivation_inputs_cannot_be_covariates(self):
        for trait, column in (("ldl", "bb_LDL"), ("t2dm", "cnt_icd10_t2dm"),
                              ("t2dm", "fod_icd10_t2dm"), ("t2dm", "t2dm.Yr2e")):
            with self.subTest(column=column), self.assertRaisesRegex(ValueError, "outcome cannot"):
                self.read(pd.DataFrame({"eid": ["001"]}), [trait], covariates=column)

    def test_abm_and_prsformer_align_the_same_outcomes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ids = [f"{i:03}" for i in range(18)]
            frame = pd.DataFrame({"eid": ids, "height": np.arange(18) + 160., "bb_LDL": np.arange(18) / 10 + 2.,
                                  "date_attend": "2010-01-01", "date_lost": None, "date_death": None,
                                  "fod_icd10_t2dm": ["2009-01-01", "2011-01-01"] * 9,
                                  "cnt_icd10_t2dm": 1, "genetic_ancestry": "EUR"})
            path = root / "phenotypes.tsv"
            frame.to_csv(path, sep="\t", index=False)
            splits = root / "splits.tsv"
            pd.DataFrame({"eid": ids, "split": ["train"] * 6 + ["validation"] * 6 + ["test"] * 6}).to_csv(splits, sep="\t", index=False)
            args = prsformer.parser().parse_args(["prepare", "--pheno-file", str(path), "--cache-dir", str(root / "cache"),
                                                 "--covariates", "none", "--remove", "none", "--split-file", str(splits)])
            expected = data.read_phenotypes(path, list(data.TRAITS), args).set_index("eid")
            with patch.object(prsformer, "sample_ids", return_value=ids[::-1]):
                actual = prsformer.make_cohort(args, [(1, "genotypes", "variants", "samples", "pgen")])
            self.assertEqual(actual.eid.tolist(), ids[::-1])
            np.testing.assert_allclose(actual[list(data.TRAITS)], expected.loc[actual.eid, list(data.TRAITS)])
            self.assertEqual(actual.attrs["phenotype_definitions"], expected.attrs["phenotype_definitions"])


if __name__ == "__main__":
    unittest.main()
