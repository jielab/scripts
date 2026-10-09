"""Regression checks for usable age evidence and honestly named negative controls.

These synthetic checks need NumPy/pandas only; no network access or UKB data.
Run: python grid/f/test_grid_evolution_review.py
"""
from pathlib import Path
import importlib.util
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("grid_evolution_review", HERE / "grid.evolution.py")
evolution = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = evolution
SPEC.loader.exec_module(evolution)


class AgeEvidence(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="grid-age-review-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        n = 12
        self.weights = pd.DataFrame({
            "CHR": 1, "BP": np.arange(101, 101 + n),
            "SNP": [f"rs{j}" for j in range(n)],
            "A1": ["C", "A"] * (n // 2), "A2": ["A", "C"] * (n // 2),
            "beta_EUR": np.linspace(-.2, .3, n),
        })
        self.annotation = self.weights[["CHR", "BP"]].assign(REF="A", ALT="C", ANC="A")
        self.annotation["AGE_GEN"] = np.tile([500., 2000., 6000.], n // 3)
        self.annotation["AGE_LO"] = self.annotation.AGE_GEN * .9
        self.annotation["AGE_HI"] = self.annotation.AGE_GEN * 1.1
        self.annotation["AGE_QUAL"] = .95
        self.annotation["AGE_SOURCE"] = "synthetic_test"
        self.annotation["AGE_METHOD"] = "synthetic_known_age"
        self.annotation["AGE_UNCERTAINTY"] = "synthetic_bounds"

    def run_builder(self, *options):
        self.weights.to_csv(self.root / "weights.tsv", sep="\t", index=False)
        self.annotation.to_csv(self.root / "annotation.tsv", sep="\t", index=False)
        args = evolution.parser().parse_args([
            "build", "--weights", str(self.root / "weights.tsv"),
            "--annotations", str(self.root / "annotation.tsv"),
            "--out-dir", str(self.root / "output"), "--populations", "EUR",
            "--chromosomes", "1", "--seed", "93", "--force", *options,
        ])
        result = evolution.build(args)
        variants = pd.read_csv(self.root / "output/variants.tsv.gz", sep="\t")
        modules = pd.read_csv(self.root / "output/modules.tsv.gz", sep="\t")
        return result, variants, modules

    def add_frequency(self):
        # Each stratum contains all three age bins.
        self.annotation["AF_EUR"] = np.repeat([.02, .2], 6)
        self.annotation["FREQ_SOURCE"] = "synthetic_external_reference"

    def test_missing_bounds_cannot_establish_usable_age(self):
        self.annotation[["AGE_LO", "AGE_HI"]] = np.nan
        with self.assertRaisesRegex(ValueError, "Only 0 usable"):
            self.run_builder()
        result, variants, _ = self.run_builder("--allow-proxy-only")
        self.assertEqual(result["mode"], "proxy_only")
        self.assertEqual(result["age_quality_pass_variants"], 12)
        self.assertEqual(result["age_usable_variants"], 0)
        self.assertTrue(variants.age_bin.eq("uncertain").all())
        self.assertTrue(variants.age_quality_pass.all())
        self.assertFalse(variants.age_usable.any())
        self.assertEqual(result["permutation_scheme"], "not_applicable")

    def test_ci_crossing_and_incompatible_ci_stay_uncertain(self):
        self.annotation.loc[0, ["AGE_LO", "AGE_HI"]] = [450, 1050]
        self.annotation.loc[1, ["AGE_LO", "AGE_HI"]] = [2100, 2500]
        self.annotation.loc[2, "AGE_QUAL"] = .5
        result, variants, _ = self.run_builder()
        self.assertEqual(result["age_usable_variants"], 9)
        self.assertEqual(variants.loc[0, "age_status"], "boundary_uncertain")
        self.assertEqual(variants.loc[1, "age_status"], "point_outside_bounds")
        self.assertEqual(variants.loc[2, "age_status"], "low_quality")
        self.assertTrue(variants.loc[:2, "age_bin"].eq("uncertain").all())
        self.assertFalse(variants.loc[:2, "age_usable"].any())
        self.assertTrue((variants.loc[:2, "perm_age_bin"] == variants.loc[:2, "age_bin"]).all())

    def test_minimum_is_checked_on_usable_nonzero_weights(self):
        self.weights.loc[:10, "beta_EUR"] = 0.
        with self.assertRaisesRegex(ValueError, "requires 2"):
            self.run_builder("--min-age-variants", "2")
        result, _, _ = self.run_builder("--min-age-variants", "2", "--allow-proxy-only")
        self.assertEqual(result["age_usable_variants"], 12)
        self.assertEqual(result["age_usable_weighted_variants"], 1)
        self.assertEqual(result["mode"], "proxy_only")

    def test_zero_weight_ages_cannot_make_human_evolution_claim(self):
        self.weights["beta_EUR"] = 0.
        with self.assertRaisesRegex(ValueError, "Only 0 usable"):
            self.run_builder()

    def test_missing_reference_maf_is_explicitly_exploratory(self):
        result, variants, modules = self.run_builder()
        self.assertEqual(result["permutation_scheme"], "chr_only_exploratory")
        self.assertEqual(result["permutation_maf_covered_fraction"], 0.)
        self.assertFalse(result["permutation_ld_controlled"])
        self.assertTrue(variants.permutation_scheme.eq("chr_only_exploratory").all())
        negative = modules.loc[modules.group.eq("permuted_age"), "description"]
        self.assertTrue(negative.str.contains("not MAF-matched", regex=False).all())
        self.assertTrue(negative.str.contains("Not LD-matched", regex=False).all())

    def test_complete_maf_strata_preserve_labels_and_signs(self):
        self.add_frequency()
        result, variants, _ = self.run_builder("--age-permutation-mode", "chr-maf")
        self.assertEqual(result["permutation_scheme"], "chr_reference_maf")
        self.assertEqual(result["permutation_maf_covered_fraction"], 1.)
        self.assertFalse(result["permutation_ld_controlled"])
        for _, group in variants.groupby("perm_stratum"):
            self.assertCountEqual(group.age_bin, group.perm_age_bin)
        matrix = pd.read_csv(self.root / "output/chr1.weights.tsv.gz", sep="\t")
        np.testing.assert_allclose(matrix.total_EUR, self.weights.beta_EUR)
        for prefix in ("age", "perm_age"):
            partitions = matrix[[f"{prefix}_{label}_EUR" for label in evolution.AGE_BINS]]
            np.testing.assert_allclose(partitions.sum(axis=1), self.weights.beta_EUR)
            present = partitions.to_numpy() != 0
            expected_sign = np.broadcast_to(np.sign(self.weights.beta_EUR.to_numpy())[:, None], present.shape)
            np.testing.assert_array_equal(np.sign(partitions.to_numpy())[present], expected_sign[present])

    def test_partial_maf_is_reported_and_strict_mode_rejects(self):
        self.add_frequency()
        self.annotation.loc[:2, "AF_EUR"] = np.nan
        result, _, _ = self.run_builder()
        self.assertEqual(result["permutation_scheme"], "chr_reference_maf_partial")
        self.assertEqual(result["permutation_maf_covered_variants"], 9)
        self.assertEqual(result["permutation_maf_covered_fraction"], .75)
        with self.assertRaisesRegex(ValueError, "3 are missing"):
            self.run_builder("--age-permutation-mode", "chr-maf")

    def test_strict_maf_needs_only_permutable_sites(self):
        self.add_frequency()
        self.annotation.loc[0, ["AGE_LO", "AGE_HI", "AF_EUR"]] = np.nan
        self.annotation.loc[1, ["ANC", "AF_EUR"]] = ["", np.nan]
        result, variants, _ = self.run_builder("--age-permutation-mode", "chr-maf")
        self.assertEqual(result["permutation_eligible_variants"], 10)
        self.assertEqual(result["permutation_maf_covered_fraction"], 1.)
        self.assertEqual(variants.loc[0, "perm_age_bin"], "uncertain")
        self.assertEqual(variants.loc[1, "perm_age_bin"], "unknown")

    def test_explicit_chromosome_only_does_not_claim_maf_control(self):
        self.add_frequency()
        result, variants, _ = self.run_builder("--age-permutation-mode", "chr-only")
        self.assertEqual(result["permutation_scheme"], "chr_only_exploratory")
        self.assertEqual(result["permutation_maf_covered_fraction"], 1.)
        self.assertTrue(variants.perm_stratum.eq("1:all").all())

    def test_homogeneous_age_bins_report_no_negative_control_contrast(self):
        self.annotation["AGE_GEN"] = 500.
        self.annotation["AGE_LO"] = 450.
        self.annotation["AGE_HI"] = 550.
        result, _, _ = self.run_builder()
        self.assertEqual(result["age_bin_count"], 1)
        self.assertEqual(result["age_partition_status"], "single_dated_bin")
        self.assertEqual(result["permutation_changed_weighted_variants"], 0)
        self.assertFalse(result["permutation_has_weighted_contrast"])


if __name__ == "__main__":
    unittest.main()
