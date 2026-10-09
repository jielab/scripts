"""Synthetic checks for exact eligible support versus the capped donor bank query.

Run: python grid/f/test_grid_matching_review.py
No individual UKB data or GPU/index dependency is needed.
"""
from pathlib import Path
import importlib.util
import pickle
import sys
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("grid_matching_review", HERE / "grid.abm.py")
abm = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = abm
SPEC.loader.exec_module(abm)


class ExactSupport(unittest.TestCase):
    def setUp(self):
        n = 24
        x = np.linspace(-2., 2., n)
        self.frame = pd.DataFrame({
            "eid": [f"D{j:02d}" for j in range(n)],
            "family_id": [f"F{j // 2:02d}" for j in range(n)],
            "ancestry": "EUR", "x": x, "y": np.sin(x),
        })
        self.config = abm._configuration({
            "abm_backend": "reference", "device": "cpu", "retrieval": "kd_tree",
            "verbose": False, "query_batch": 3, "k_grid": [3, 5],
            "radius_grid": [1., 2., 4.], "alpha_grid": [0., .5, 1.],
        })
        self.bank = abm.ReferenceBank().fit(
            self.frame, np.sin(x), np.zeros(n), np.arange(n) % 3,
            {"csx": ["x"]}, self.config,
        )

    def query(self):
        return pd.DataFrame({
            "eid": ["new", "D00", "D02", "same-family-new", "missing"],
            "family_id": ["unseen", "F00", "F00", "F03", "unseen"],
            "ancestry": "EUR", "x": [0., -2., 1., -.5, np.nan],
        })

    def brute_count(self, frame, caliper, qc=None):
        z = self.bank.geometry.transform(frame)
        # Dense distances are used only as this tiny independent test oracle.
        distance = np.linalg.norm(z[:, None, :] - self.bank.z[None, :, :], axis=2)
        eligible = distance <= caliper
        eligible &= frame.eid.to_numpy()[:, None] != self.bank.ids[None, :]
        eligible &= frame.family_id.to_numpy()[:, None] != self.bank.groups[None, :]
        if qc is not None:
            eligible &= np.asarray(qc)[:, None]
        return eligible.sum(axis=1)

    def test_exact_count_includes_all_donors_beyond_k(self):
        frame = self.query().iloc[:1]
        _, diagnostics, packed = self.bank.borrow(
            frame, np.zeros(len(frame)), k=3, alpha=1., radius_multiplier=100., retain=True,
        )
        self.assertEqual(diagnostics.matched_count.iloc[0], 3)
        self.assertEqual(diagnostics.eligible_match_count.iloc[0], len(self.frame))
        self.assertEqual(diagnostics.retrieval_k.iloc[0], 3)
        self.assertTrue(diagnostics.matches_truncated.iloc[0])
        self.assertEqual(packed["eligible_match_count"][0], 24)
        self.assertEqual(packed["neighbor_indices"].shape, (1, 3))
        # Full support is explanatory only: actual three-donor ESS still
        # determines the original borrowing fraction and decomposition.
        used_weights = packed["weights"][0]
        self.assertAlmostEqual(diagnostics.donor_ess.iloc[0], 1. / np.sum(used_weights ** 2))
        np.testing.assert_allclose(
            packed["prediction"], packed["baseline"] + packed["weighted_residual_contribution"].sum(axis=1),
        )

    def test_family_and_identity_exclusions_are_unioned_once(self):
        frame = self.query().iloc[:4]
        z = self.bank.geometry.transform(frame)
        count = self.bank.count_eligible(z, frame.eid, frame.family_id, caliper=100.)
        # New; self within the excluded family; self outside the named
        # excluded family; unseen ID but known two-person family.
        np.testing.assert_array_equal(count, [24, 22, 21, 22])
        np.testing.assert_array_equal(count, self.brute_count(frame, 100.))

    def test_radius_and_qc_match_independent_count_oracle(self):
        frame = self.query()
        qc = frame.x.notna().to_numpy()
        z = self.bank.geometry.transform(frame)
        for caliper in (.1, .5, 1., 2.):
            count = self.bank.count_eligible(z, frame.eid, frame.family_id, caliper, qc)
            np.testing.assert_array_equal(count, self.brute_count(frame, caliper, qc))
        _, diagnostics, _ = self.bank.borrow(
            frame, np.zeros(len(frame)), k=5, alpha=1., radius_multiplier=100.,
        )
        self.assertEqual(diagnostics.eligible_match_count.iloc[-1], 0)
        self.assertEqual(diagnostics.matched_count.iloc[-1], 0)
        self.assertFalse(diagnostics.matches_truncated.iloc[-1])
        self.assertEqual(diagnostics.fallback_reason.iloc[-1], "excess_missing_predictors")

    def test_untruncated_small_support_has_correct_flag(self):
        frame = self.frame.iloc[[10]].drop(columns="y")
        caliper = .22
        _, diagnostics, _ = self.bank.borrow(
            frame, np.zeros(1), k=5, alpha=1., radius_multiplier=caliper / self.bank.radius,
        )
        expected = self.brute_count(frame, caliper)
        np.testing.assert_array_equal(diagnostics.eligible_match_count, expected)
        np.testing.assert_array_equal(diagnostics.matched_count, expected)
        self.assertFalse(diagnostics.matches_truncated.any())

    def test_cached_counts_preserve_predictions_and_enforce_qc(self):
        frame = self.query()
        z = self.bank.geometry.transform(frame)
        counts_without_qc = self.bank.count_eligible(z, frame.eid, frame.family_id, self.bank.radius * 4)
        first = self.bank.borrow(frame, np.zeros(len(frame)), k=5, alpha=.5, radius_multiplier=4.)
        cached = self.bank.borrow(
            frame, np.zeros(len(frame)), k=5, alpha=.5, radius_multiplier=4.,
            eligible_counts=counts_without_qc,
        )
        np.testing.assert_array_equal(first[0], cached[0])
        pd.testing.assert_frame_equal(first[1], cached[1])
        self.assertEqual(cached[1].eligible_match_count.iloc[-1], 0)

    def test_tuning_computes_counts_once_per_radius(self):
        frame = self.query().iloc[:4]
        with patch.object(self.bank, "count_eligible", wraps=self.bank.count_eligible) as counter:
            best, rows = abm._tune_bank(
                "test", self.bank, frame, np.zeros(len(frame)), np.array([.1, -.2, .3, 0.]), self.config,
            )
        self.assertEqual(counter.call_count, len(self.config["radius_grid"]))
        self.assertEqual(len(rows), len(self.config["k_grid"]) * len(self.config["radius_grid"]) * len(self.config["alpha_grid"]))
        self.assertIn(best["k"], self.config["k_grid"])

    def test_count_index_is_rebuilt_after_serialization(self):
        frame = self.query().iloc[:4]
        z = self.bank.geometry.transform(frame)
        before = self.bank.count_eligible(z, frame.eid, frame.family_id, .8)
        self.assertIsNotNone(self.bank._count_tree)
        state = self.bank.__getstate__()
        self.assertFalse(any(name.startswith("_count_") for name in state))
        restored = pickle.loads(pickle.dumps(self.bank))
        self.assertIsNone(restored._count_tree)
        after = restored.count_eligible(z, frame.eid, frame.family_id, .8)
        np.testing.assert_array_equal(before, after)

    def test_exact_count_can_use_auxiliary_tree_without_hnsw_import(self):
        # count_eligible must be exact and independent of the ANN library.
        self.bank.config["retrieval"] = "hnsw"
        self.bank._count_tree = None
        frame = self.query().iloc[:4]
        z = self.bank.geometry.transform(frame)
        result = self.bank.count_eligible(z, frame.eid, frame.family_id, .8)
        np.testing.assert_array_equal(result, self.brute_count(frame, .8))

    def test_nonpositive_caliper_disables_eligible_support(self):
        frame = self.query().iloc[:4]
        z = self.bank.geometry.transform(frame)
        for caliper in (0., -1.):
            np.testing.assert_array_equal(
                self.bank.count_eligible(z, frame.eid, frame.family_id, caliper), np.zeros(len(frame)),
            )


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        unittest.main()
