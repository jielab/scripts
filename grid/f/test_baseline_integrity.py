#!/usr/bin/env python3
"""Targeted guards for wrong-allele centering and outcome leakage.

Run: python3 grid/f/test_baseline_integrity.py
Uses synthetic public-free inputs; no R, PLINK, Torch, or UKB files required.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def module(name, filename):
	spec = importlib.util.spec_from_file_location(name, HERE / filename)
	value = importlib.util.module_from_spec(spec)
	sys.modules[name] = value
	spec.loader.exec_module(value)
	return value


common = module("baseline_common_integrity", "0.common.py")
prsformer = module("baseline_prsformer_integrity", "3.prsformer.py")
prsformer_data = module("baseline_prsformer_data_integrity", "3.prsformer_data.py")


class AlleleFrequencyTests(unittest.TestCase):
	def frequency(self, table):
		columns = {key: common.choose(table.columns, key) for key in common.ALIASES}
		return common.effect_allele_frequency(table, columns, table.A1, table.A2)

	def test_alt_frequency_is_flipped_when_effect_is_reference(self):
		table = pd.DataFrame({"A1": ["G", "A"], "A2": ["A", "G"], "ALT": ["G", "G"], "POOLED_ALT_AF": [.2, .2]})
		frequency, _ = self.frequency(table)
		np.testing.assert_allclose(frequency, [.2, .8])

	def test_reference_column_can_identify_the_alternative(self):
		table = pd.DataFrame({"A1": ["A", "G"], "A2": ["G", "A"], "REF": ["A", "A"], "ALT_FREQ": [.2, .2]})
		frequency, _ = self.frequency(table)
		np.testing.assert_allclose(frequency, [.8, .2])

	def test_unoriented_af_is_not_silently_assumed_effect_frequency(self):
		table = pd.DataFrame({"A1": ["A"], "A2": ["G"], "AF": [.2]})
		frequency, source = self.frequency(table)
		self.assertTrue(frequency.isna().all())
		self.assertIn("unspecified", source)
		table["AF_ALLELE"] = "G"
		frequency, _ = self.frequency(table)
		np.testing.assert_allclose(frequency, [.8])

	def test_explicit_effect_frequency_takes_precedence_and_rejects_invalid(self):
		table = pd.DataFrame({"A1": ["A", "A"], "A2": ["G", "G"], "EAF": [.3, 1.1], "ALT_FREQ": [.7, .2]})
		frequency, _ = self.frequency(table)
		self.assertAlmostEqual(frequency.iloc[0], .3)
		self.assertTrue(np.isnan(frequency.iloc[1]))

	def test_beta_p_inputs_pass_preflight_and_get_same_normalized_z(self):
		with tempfile.TemporaryDirectory() as td:
			root = Path(td)
			ref, source = root / "snpinfo.tsv", root / "source.tsv"
			pd.DataFrame({"CHR": [1, 1], "SNP": ["rs1", "rs2"], "BP": [100, 200]}).to_csv(ref, sep="\t", index=False)
			pd.DataFrame({"CHR": [1, 1], "SNP": ["rs1", "rs2"], "BP": [100, 200],
				"A1": ["A", "G"], "A2": ["G", "A"], "ALT": ["G", "G"],
				"BETA": [.2, -.1], "P": [.05, .01], "ALT_FREQ": [.2, .2], "N": [10000, 10000]}).to_csv(source, sep="\t", index=False)
			common.inspect(str(ref), [str(source)])
			output, metadata = root / "normalized.gz", root / "normalized.json"
			proc = subprocess.run([sys.executable, str(HERE / "0.common.py"), "prepare-sumstats", "--input", str(source),
				"--output", str(output), "--metadata", str(metadata), "--snpinfo", str(ref), "--trait", "height", "--pop", "EUR"],
				capture_output=True, text=True)
			self.assertEqual(proc.returncode, 0, proc.stderr)
			prepared = pd.read_csv(output, sep="\t")
			np.testing.assert_allclose(prepared.EAF, [.8, .2])
			np.testing.assert_allclose(prepared.BETA / prepared.SE, [1.95996398454, -2.57582930355], rtol=1e-8)
			self.assertEqual(json.loads(metadata.read_text())["eaf_unavailable"], 0)


class HoldoutTests(unittest.TestCase):
	def setUp(self):
		self.roster = pd.DataFrame({"eid": ["1", "2", "3", "4", "5", "6"],
			"split": ["train", "train", "validation", "validation", "test", "test"],
			"group": ["f1", "f2", "f3", "f4", "f5", "f6"]})
		self.checkpoint = {"split_registry": prsformer.split_registry(self.roster), "training_identity": {"data_sha256": "original"}}

	def test_original_held_out_people_pass(self):
		prsformer.validate_prediction_holdout(self.roster, {}, self.checkpoint)

	def test_relabeling_training_or_validation_people_is_rejected(self):
		for person in ["1", "3", "1.0"]:
			query = pd.DataFrame({"eid": [person], "split": ["test"], "group": ["new-family"]})
			with self.assertRaisesRegex(ValueError, "not held out"):
				prsformer.validate_prediction_holdout(query, {}, self.checkpoint)

	def test_different_person_from_a_development_family_is_rejected(self):
		query = pd.DataFrame({"eid": ["new"], "split": ["test"], "group": ["f1"]})
		with self.assertRaisesRegex(ValueError, "families overlap"):
			prsformer.validate_prediction_holdout(query, {}, self.checkpoint)
		with self.assertRaisesRegex(ValueError, "related-family column"):
			prsformer.validate_prediction_holdout(query.drop(columns="group"), {}, self.checkpoint)

	def test_legacy_checkpoint_only_accepts_exact_original_cohort(self):
		legacy = {"training_identity": {"data_sha256": "original"}}
		prsformer.validate_prediction_holdout(self.roster, {"data_sha256": "original"}, legacy)
		with self.assertRaisesRegex(ValueError, "Legacy checkpoint"):
			prsformer.validate_prediction_holdout(self.roster, {"data_sha256": "modified"}, legacy)

	def test_direct_training_reader_checks_family_split(self):
		with tempfile.TemporaryDirectory() as td:
			root = Path(td)
			d = self.roster.assign(target="EUR", height=np.arange(6, dtype=float))
			d.loc[d.eid == "5", "group"] = "f1"
			d.to_csv(root / "data.tsv", sep="\t", index=False)
			pd.DataFrame({"CHR": [1], "BP": [100], "SNP": ["rs1"], "REF": ["A"], "ALT": ["G"]}).to_csv(root / "variants.tsv", sep="\t", index=False)
			np.save(root / "genotypes.npy", np.ones((6, 1), dtype=np.float16))
			args = SimpleNamespace(data=root / "data.tsv", variants=root / "variants.tsv", genotypes=root / "genotypes.npy")
			with self.assertRaisesRegex(ValueError, "spans train"):
				prsformer.read_inputs(args, ["height"])

	def test_legacy_yeval_learned_score_path_fails_before_any_data_access(self):
		proc = subprocess.run(["bash", str(HERE.parent / "Yeval.sh"), "--grid-file", "/does/not/exist"], capture_output=True, text=True)
		self.assertEqual(proc.returncode, 2)
		self.assertIn("cannot be assigned new CV folds", proc.stderr)
		self.assertIn("grid.sh --stage report", proc.stderr)


class ScoreProvenanceTests(unittest.TestCase):
	def test_appending_meta_preserves_valid_population_lineage(self):
		with tempfile.TemporaryDirectory() as td:
			root = Path(td)
			score = root / "1csx.scores.rds"
			populations = ("AFR", "EAS", "EUR", "SAS")
			weights = {f"csx.{pop}": root / f"{pop}.weights" for pop in populations}
			for path in weights.values():
				path.write_text("SNP\tA1\tBETA\nrs1\tA\t0.1\n")
			records = common.score_source_records(weights, "joint-settings", [1, 2])
			table = pd.DataFrame({"eid": ["1", "2"], **{key: [.1, .2] for key in weights}})
			# Exercise publication/locking/hash semantics independently of R encoding.
			with patch.object(common, "write_rds", side_effect=lambda d, p: d.to_csv(p, sep="\t", index=False)), \
				patch.object(common, "read_result_table", side_effect=lambda p, **kw: pd.read_csv(p, sep="\t", dtype={"eid": str})):
				common.update_csx_table(table, score, provenance=records)
				meta_weight = root / "meta.weights"
				meta_weight.write_text("SNP\tA1\tBETA\nrs1\tA\t0.2\n")
				common.update_csx_table(pd.DataFrame({"eid": ["1", "2"], "csx.meta": [.3, .4]}), score,
					provenance=common.score_source_records({"csx.meta": meta_weight}, "meta-settings", [1, 2]))
				provenance = json.loads(score.with_suffix(".provenance.json").read_text())
				self.assertTrue(provenance["complete_population_provenance"])
				self.assertEqual(set(provenance["models"]), set(weights) | {"csx.meta"})
				self.assertEqual(provenance["score_file_sha256"], common.sha256(score))
				self.assertEqual(provenance["models"]["csx.EUR"]["weights_sha256"], common.sha256(weights["csx.EUR"]))
				# Replacing a column without attestation must clear that column only.
				common.update_csx_table(pd.DataFrame({"eid": ["1", "2"], "csx.EUR": [9., 8.]}), score)
				provenance = json.loads(score.with_suffix(".provenance.json").read_text())
				self.assertFalse(provenance["complete_population_provenance"])
				self.assertNotIn("csx.EUR", provenance["models"])
				self.assertIn("csx.AFR", provenance["models"])
				# If the source table itself was edited, no old lineage survives.
				with score.open("a") as stream:
					stream.write("\n")
				common.update_csx_table(pd.DataFrame({"eid": ["1", "2"], "csx.meta": [.5, .6]}), score)
				provenance = json.loads(score.with_suffix(".provenance.json").read_text())
				self.assertEqual(provenance["models"], {})

	def test_changed_weights_are_rejected_before_publication(self):
		with tempfile.TemporaryDirectory() as td:
			root = Path(td)
			weight = root / "weights.tsv"
			weight.write_text("rs1 A .1\n")
			records = common.score_source_records({"csx.EUR": weight}, "settings", [1])
			weight.write_text("rs1 A .2\n")
			with self.assertRaisesRegex(ValueError, "weight file changed"):
				common.update_csx_table(pd.DataFrame({"eid": ["1"], "csx.EUR": [.1]}), root / "scores.rds", provenance=records)
			self.assertFalse((root / "scores.rds").exists())


class BenchmarkFieldTests(unittest.TestCase):
	def test_yeval_metric_names_keep_correlation_and_predictive_formulas_distinct(self):
		# R is unavailable in this test environment; protect the published field contract.
		source = (HERE / "Yeval.R").read_text()
		wrapper = (HERE.parent / "Yeval.sh").read_text()
		self.assertNotIn("OOF_prediction_R2", source)
		self.assertIn('ct = "OOF_residual_correlation_R2"', source)
		self.assertRegex(source, r"OOF_predictive_R2\s*=\s*1\s*-\s*sse1/den")
		self.assertRegex(source, r"full_R2\s*=\s*1\s*-\s*sse1/den")
		self.assertIn('ct = "Residual correlation R²"', source)
		self.assertNotIn('ct = "Prediction R²"', source)
		self.assertNotIn('ct = "prediction R²"', source)
		self.assertIn("OOF_predictive_R2", wrapper)

	def test_predictive_r2_penalizes_opposite_sign_perfect_correlation(self):
		y = np.array([-1., 0., 1.])
		metrics = prsformer.metrics_for(y, -y, np.zeros(3), -y, "height")
		self.assertAlmostEqual(metrics["residual_correlation_R2"], 1.)
		self.assertAlmostEqual(metrics["predictive_R2"], -3.)
		self.assertEqual(metrics["predictive_R2"], metrics["full_R2"])
		self.assertNotIn("prediction_R2", metrics)

	def test_legacy_metric_migration_preserves_original_formula_values(self):
		old = pd.DataFrame({"trait": ["height", "height"], "target": ["ALL", "ALL"], "split": ["test", "test"],
			"metric": ["prediction_R2", "full_R2"], "value": [.8, -.2]})
		updated = prsformer_data.canonical_metrics(old).set_index("metric")
		self.assertAlmostEqual(updated.loc["residual_correlation_R2", "value"], .8)
		self.assertAlmostEqual(updated.loc["predictive_R2", "value"], -.2)
		self.assertEqual(len(updated), 3)
		self.assertEqual(len(prsformer_data.canonical_metrics(updated.reset_index())), 3)

	def test_covariate_values_follow_subject_ids_not_row_order(self):
		scores = pd.DataFrame({"eid": ["2", "1"], "prediction": [.4, .3]})
		cohort = pd.DataFrame({"eid": ["1", "2"], "age": [40., 60.], "sex": [0., 1.]})
		result = prsformer_data.attach_covariate_values(scores, cohort, "age,sex")
		self.assertEqual(result.eid.tolist(), ["2", "1"])
		np.testing.assert_allclose(result["covariate.age"], [60., 40.])
		np.testing.assert_allclose(result["covariate.sex"], [1., 0.])
		with self.assertRaisesRegex(ValueError, "no matching original covariates"):
			prsformer_data.attach_covariate_values(scores, cohort.iloc[:1], "age,sex")
		with self.assertRaisesRegex(ValueError, "missing or nonfinite"):
			prsformer_data.attach_covariate_values(scores, cohort.assign(age=[40., np.nan]), "age,sex")


if __name__ == "__main__":
	unittest.main(verbosity=2)
