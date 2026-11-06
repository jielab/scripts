"""Frozen-model, identity, label-separation and deployment-policy contracts."""
from pathlib import Path
import json
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from abm_core.model import ABMAssessor


def make_role(prefix, n_each=8):
    counts = np.repeat([50, 100], n_each)
    low_count = counts == 50
    return pd.DataFrame({
        "id": [f"{prefix}-{i:03d}" for i in range(len(counts))],
        "group_id": [f"group-{prefix}-{i:03d}" for i in range(len(counts))],
        "n_twins": counts,
        "supported": True,
        "p_base": np.zeros(len(counts)),
        "p_candidate": np.ones(len(counts)),
        # Low-count subjects have base loss 1, gain +1. High-count subjects
        # have base loss .04, gain -.60. Predictability and benefit differ.
        "y": np.where(low_count, 1.0, 0.2),
        "ancestry": np.where(np.arange(len(counts)) % 2, "EUR", "EAS"),
    })


def manifest():
    return {
        "baseline_id": "csx-linear",
        "candidate_id": "csx-ldblock-net",
        "bank_sha256": "1" * 64,
        "geometry_sha256": "2" * 64,
        "predictor_bundle_sha256": "3" * 64,
        "role_manifest_sha256": "4" * 64,
        "grouping_status": "synthetic_verified",
    }


def fit_model(**kwargs):
    options = dict(manifest=manifest(), task="continuous",
                   head_options={"n_bins": 2, "risk_bins": 2,
                                 "prior_strength": 2,
                                 "min_calibration_groups": 2,
                                 "n_boot": 20, "seed": 7})
    options.update(kwargs)
    return ABMAssessor.fit(make_role("gate"), make_role("headcal"),
                           make_role("policycal"), **options)


class ModelTests(unittest.TestCase):
    def test_ancestry_is_only_a_label_for_fit_assessment_and_control_matching(self):
        from abm_core.evaluation import make_matched_controls
        model = fit_model()
        query = make_role("query", n_each=30).drop(columns="y")
        expected = model.assess(query, manifest=manifest())
        pairs, _ = make_matched_controls(expected, control_eligible_col="control_eligible")
        for changed in (query.assign(ancestry=None), query.assign(ancestry="RARE"),
                        query.drop(columns="ancestry")):
            result = model.assess(changed, manifest=manifest())
            cols = [c for c in expected if c != "ancestry"]
            pd.testing.assert_frame_equal(expected[cols], result[cols])
            observed, _ = make_matched_controls(result, control_eligible_col="control_eligible")
            pd.testing.assert_frame_equal(pairs, observed)
        other = ABMAssessor.fit(make_role("gate").assign(ancestry=None),
            make_role("headcal").drop(columns="ancestry"), make_role("policycal").assign(ancestry="RARE"),
            manifest=manifest(), head_options={"n_bins":2,"risk_bins":2,"prior_strength":2,
                "min_calibration_groups":2,"n_boot":20,"seed":7})
        pd.testing.assert_frame_equal(expected, other.assess(query, manifest=manifest()))

    def test_roundtrip_is_strict_json_and_preserves_all_assessments(self):
        model = fit_model()
        query = make_role("query").drop(columns="y")
        expected = model.assess(query, manifest=manifest())
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "model.json"
            model.save(path)
            def reject_constant(value):
                raise AssertionError(f"Nonstandard JSON constant: {value}")
            parsed = json.loads(path.read_text(), parse_constant=reject_constant)
            self.assertEqual(parsed["rules"]["benefit_c100"]["threshold"], "-Infinity")
            restored = ABMAssessor.load(path)
            observed = restored.assess(query, manifest=manifest())
        pd.testing.assert_frame_equal(expected, observed)

    def test_unavailable_bin_nan_interval_roundtrip(self):
        gate, headcal, policycal = make_role("gate"), make_role("headcal"), make_role("policycal")
        headcal = headcal.loc[headcal.n_twins == 50].copy()
        model = ABMAssessor.fit(gate, headcal, policycal, manifest=manifest(),
                                head_options={"n_bins": 2, "risk_bins": 2,
                                              "min_calibration_groups": 2, "n_boot": 20})
        query = make_role("query").drop(columns="y")
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "model.json"
            model.save(path)
            restored = ABMAssessor.load(path)
            out = restored.assess(query, manifest=manifest())
        unknown = out.n_twins == 100
        self.assertTrue(out.loc[unknown, "expected_gain_ci_low"].isna().all())
        self.assertFalse(out.loc[unknown, "assessment_available"].any())
        np.testing.assert_array_equal(out.loc[unknown, "p_policy"], out.loc[unknown, "p_base"])

    def test_bank_and_model_change_require_recalibration(self):
        model = fit_model()
        query = make_role("query").drop(columns="y")
        for key in ("bank_sha256", "geometry_sha256", "predictor_bundle_sha256", "role_manifest_sha256"):
            changed = manifest()
            changed[key] = "a" * 64
            with self.subTest(field=key), self.assertRaisesRegex(ValueError, key):
                model.assess(query, manifest=changed)

    def test_query_outcome_rejected_even_when_all_missing(self):
        model = fit_model()
        for y in (np.nan, 1.0):
            query = make_role("query")
            query["y"] = y
            with self.subTest(y=y), self.assertRaisesRegex(ValueError, "outcome"):
                model.assess(query, manifest=manifest())

    def test_threshold_is_stable_under_query_permutation_and_batches(self):
        model = fit_model()
        query = make_role("query", n_each=20).drop(columns="y")
        whole = model.assess(query, manifest=manifest()).set_index("id").sort_index()
        shuffled = query.sample(frac=1, random_state=71)
        # Deliberately preserve non-default dataframe indices to exercise joins.
        pieces = [model.assess(shuffled.iloc[start:start + 7], manifest=manifest())
                  for start in range(0, len(shuffled), 7)]
        batched = pd.concat(pieces).set_index("id").sort_index()
        columns = [c for c in whole if c.startswith("selected_")] + [
            "selected", "switched", "p_policy", "match_stratum", "expected_gain"]
        pd.testing.assert_frame_equal(whole[columns], batched[columns])

    def test_same_family_cannot_appear_in_different_training_roles(self):
        gate, headcal, policycal = make_role("gate"), make_role("headcal"), make_role("policycal")
        headcal.loc[0, "group_id"] = gate.loc[0, "group_id"]
        with self.assertRaisesRegex(ValueError, "group_id"):
            ABMAssessor.fit(gate, headcal, policycal, manifest=manifest())

    def test_negative_expected_gain_never_routes_to_candidate(self):
        model = fit_model()
        out = model.assess(make_role("query", n_each=30).drop(columns="y"), manifest=manifest())
        negative = out.expected_gain < 0
        self.assertTrue(negative.any())
        self.assertFalse(out.loc[negative, "switched"].any())
        for name in [c for c in out if c.startswith("selected_benefit_")]:
            self.assertFalse(out.loc[negative, name].any())
        np.testing.assert_array_equal(out.loc[negative, "p_policy"], out.loc[negative, "p_base"])

    def test_negative_or_nonfinite_gain_threshold_is_rejected(self):
        for threshold in (-0.01, np.nan, np.inf):
            with self.subTest(threshold=threshold), self.assertRaises(ValueError):
                fit_model(min_gain=threshold)

    def test_candidate_missingness_falls_back_without_dropping_queries(self):
        model = fit_model()
        query = make_role("query", n_each=30).drop(columns="y")
        query.loc[query.n_twins == 50, "p_candidate"] = np.nan
        out = model.assess(query, manifest=manifest())
        self.assertEqual(len(out), len(query))
        self.assertFalse(out.switched.any())
        np.testing.assert_array_equal(out.p_policy, out.p_base)

    def test_dataset_expectation_does_not_extrapolate_unavailable_people(self):
        model = fit_model()
        query = make_role("query").drop(columns="y")
        query.loc[0, "n_twins"] = 1000
        out = model.assess(query, manifest=manifest())
        summary = model.dataset_assessment(out)
        self.assertEqual(summary["n_query"], len(query))
        self.assertEqual(summary["n_assessment_available"], len(query) - 1)
        self.assertIsNone(summary["full_dataset_expected_loss"])
        self.assertIsNotNone(summary["expected_mean_loss_assessable_subset"])


if __name__ == "__main__":
    unittest.main()
