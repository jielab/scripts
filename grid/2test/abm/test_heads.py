import unittest
import numpy as np
import pandas as pd

from abm_core.heads import CountHeads, FrozenRule


def frame(prefix, repeats=30, harmful=False):
    count = np.tile([30, 60, 120], repeats)
    # Deliberately nonmonotone count->loss: high support need not mean low noise.
    y = np.tile([1.0, 0.2, 2.0], repeats)
    candidate = -np.ones(len(y)) if harmful else np.tile([0.9, -0.4, 1.2], repeats)
    ids = [f"{prefix}{i}" for i in range(len(y))]
    return pd.DataFrame({"id": ids, "group_id": ids, "n_twins": count,
                         "supported": True, "p_base": 0.0, "p_candidate": candidate, "y": y})


class CountHeadTests(unittest.TestCase):
    def fit_head(self, harmful=False):
        return CountHeads.fit(frame("g", harmful=harmful), frame("c", harmful=harmful),
                              n_bins=3, min_calibration_groups=2, n_boot=20, prior_strength=0)

    def test_count_loss_not_forced_monotone_and_gain_is_distinct(self):
        head = self.fit_head()
        query = frame("q", repeats=1).drop(columns="y")
        out = head.predict(query)
        self.assertGreater(out.expected_loss_base.iloc[2], out.expected_loss_base.iloc[1])
        self.assertGreater(out.expected_loss_base.iloc[0], out.expected_loss_base.iloc[1])
        self.assertLess(out.expected_gain.iloc[1], 0)
        self.assertGreater(out.expected_gain.iloc[0], 0)

    def test_negative_absolute_gain_never_switches(self):
        head = self.fit_head(harmful=True)
        query = head.predict(frame("q", harmful=True).drop(columns="y"))
        rule = FrozenRule.fit(query, objective="benefit", coverage=1.0)
        self.assertFalse(rule.apply(query).any())

    def test_missing_candidate_cannot_redefine_prediction_difficulty(self):
        g = frame("g")
        g.loc[0, "p_candidate"] = np.nan
        with self.assertRaisesRegex(ValueError, "candidate missingness"):
            CountHeads.fit(g, frame("c"), n_boot=20, min_calibration_groups=2)

    def test_query_outcomes_forbidden(self):
        with self.assertRaisesRegex(ValueError, "Inference features"):
            self.fit_head().predict(frame("q"))

    def test_unsupported_means_unknown_not_zero_error(self):
        q = frame("q", repeats=1).drop(columns="y")
        q.loc[0, "supported"] = False
        result = self.fit_head().predict(q)
        self.assertFalse(result.assessment_available.iloc[0])
        self.assertTrue(np.isnan(result.expected_loss_base.iloc[0]))

    def test_role_family_overlap_rejected(self):
        c = frame("c")
        c.loc[0, "group_id"] = "g0"
        with self.assertRaisesRegex(ValueError, "group_id"):
            CountHeads.fit(frame("g"), c)

    def test_unavailable_genetic_count_is_not_fabricated_zero_or_dropped(self):
        q = frame("q", repeats=1).drop(columns="y")
        q.loc[0, "supported"] = False
        q["n_twins"] = q.n_twins.astype(float)
        q.loc[0, "n_twins"] = np.nan
        result = self.fit_head().predict(q)
        self.assertEqual(len(result), len(q))
        self.assertTrue(np.isnan(result.n_twins.iloc[0]))
        self.assertFalse(result.assessment_available.iloc[0])
        self.assertEqual(result.assessment_reason.iloc[0], "genetic_count_unavailable")

    def test_fixed_threshold_is_invariant_to_query_composition(self):
        head = self.fit_head()
        cal = head.predict(frame("r").drop(columns="y"))
        rule = FrozenRule.fit(cal, objective="predictability", coverage=0.5)
        q = head.predict(frame("q").drop(columns="y"))
        all_selected = pd.Series(rule.apply(q), index=q.id)
        subset = q.sample(15, random_state=5)
        np.testing.assert_array_equal(rule.apply(subset), all_selected.loc[subset.id])


if __name__ == "__main__":
    unittest.main()
