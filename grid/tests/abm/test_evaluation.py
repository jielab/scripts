import json
import unittest

import numpy as np
import pandas as pd

from abm_core.evaluation import evaluate_predictions, make_matched_controls, individual_loss
from abm_core.losses import row_loss


def toy():
    return pd.DataFrame({
        "id": ["s1", "s2", "c1", "c2"],
        "group_id": ["f1", "f2", "f3", "f4"],
        "y": [0., 0., 0., 0.],
        "p_base": [2., 2., 2., 2.],
        "p_candidate": [1., 1., 2., 2.],
        "supported": [True] * 4,
        "selected": [True, True, False, False],
        "match_stratum": ["EUR|risk0"] * 4,
    })


def evaluate(frame, **kwargs):
    return evaluate_predictions(
        frame, task="continuous", loss_name="mse", selector="benefit",
        candidate_id="csx-block-direct", n_boot=kwargs.pop("n_boot", 100),
        min_matched_pairs=kwargs.pop("min_matched_pairs", 2),
        min_matched_groups=kwargs.pop("min_matched_groups", 2), **kwargs,
    )


class EvaluationTests(unittest.TestCase):
    def test_gain_enrichment_is_difference_of_gains_and_policy_is_whole_population(self):
        result = evaluate(toy(), coverage=.5)
        self.assertAlmostEqual(result["matched"]["gain_enrichment"]["estimate"], 3.)
        self.assertAlmostEqual(result["matched"]["predictability_loss_difference"]["estimate"], 0.)
        self.assertAlmostEqual(result["whole_cohort"]["policy"]["mse"], 2.5)
        self.assertAlmostEqual(result["whole_cohort"]["policy_gain_vs_base"]["estimate"], 1.5)
        assert result["whole_cohort"]["base"]["prediction_r2"] is None
        if result["matched"]["gain_enrichment"]["ci_low"] is not None:
            self.assertAlmostEqual(result["matched"]["gain_enrichment"]["ci_low"], 3.)
        json.dumps(result, allow_nan=False)


    def test_matching_is_independent_of_labels_predictions_and_row_order(self):
        frame = toy()
        pairs1, _ = make_matched_controls(frame)
        frame.loc[0, "y"] = np.nan
        frame["p_candidate"] = [np.nan, 99., -9., -19.]
        pairs2, _ = make_matched_controls(frame.sample(frac=1, random_state=1))
        assignments = lambda pairs: sorted(zip(pairs.selected_id, pairs.control_id))
        assert assignments(pairs1) == assignments(pairs2)
        result = evaluate(frame, n_boot=0)
        assert result["coverage"]["n_selected_query"] == 2
        assert result["coverage"]["n_selected_evaluable"] == 1
        assert result["matching"]["n_pairs"] == 2
        assert result["matching"]["n_pairs_with_both_outcomes"] == 1
        assert result["policy"]["n_selected_without_candidate"] == 1


    def test_control_excludes_entire_selected_family(self):
        frame = toy()
        frame.loc[2, "group_id"] = frame.loc[0, "group_id"]
        pairs, diagnostics = make_matched_controls(frame)
        assert len(pairs) == 1
        assert "c1" not in set(pairs.control_id)
        assert diagnostics["n_unmatched_selected"] == 1
        assert diagnostics["n_excluded_unselected_same_group"] == 1
        result = evaluate(frame, n_boot=0)
        assert result["coverage"]["n_switched_query"] == 2  # A shortage never changes policy.


    def test_no_support_means_base_policy_and_selection_must_respect_support(self):
        frame = toy()
        frame["supported"] = False
        with self.assertRaisesRegex(ValueError, "imply supported"):
            evaluate(frame)
        frame["selected"] = False
        result = evaluate(frame, n_boot=0)
        assert result["whole_cohort"]["policy"] == result["whole_cohort"]["base"]
        assert result["coverage"]["n_switched_query"] == 0
        assert result["matched"]["gain_enrichment"]["estimate"] is None


    def test_all_selected_has_no_disjoint_control_but_valid_policy(self):
        frame = toy()
        frame["selected"] = True
        result = evaluate(frame, n_boot=0, coverage=1.)
        assert result["matching"]["status"] == "no_disjoint_matched_controls"
        assert result["matched"]["gain_enrichment"]["status"] == "not_estimable"
        assert result["whole_cohort"]["policy"] == result["whole_cohort"]["candidate"]


    def test_missing_stratum_no_expansion_or_backfill(self):
        frame = toy()
        frame.loc[0, "match_stratum"] = "EAS|risk0"
        frame.loc[1, "match_stratum"] = None
        pairs, diagnostics = make_matched_controls(frame)
        assert pairs.empty
        assert diagnostics["n_selected_missing_stratum"] == 1
        assert diagnostics["n_unmatched_selected"] == 2


    def test_binary_metrics_are_probabilities_and_two_objectives_can_disagree(self):
        frame = toy()
        frame["y"] = [0., 1., 0., 1.]
        frame["p_base"] = [.1, .9, .4, .6]
        frame["p_candidate"] = [.1, .9, .1, .9]
        result = evaluate_predictions(
            frame, task="binary", loss_name="brier", selector="predictability",
            candidate_id="csx-block-direct", n_boot=0,
        )
        # Selected persons are easier for baseline, while candidate gains are in controls.
        assert result["matched"]["predictability_loss_difference"]["estimate"] < 0
        assert result["matched"]["gain_enrichment"]["estimate"] < 0
        assert result["whole_cohort"]["base"]["auc"] == 1.
        assert result["whole_cohort"]["base"]["cases"] == 2
        assert result["whole_cohort"]["policy_gain_vs_base"]["estimate"] == 0.
        frame.loc[0, "p_candidate"] = 3.
        with self.assertRaisesRegex(ValueError, "probabilities"):
            evaluate_predictions(frame, task="binary", loss_name="brier", selector="x", candidate_id="x")


    def test_group_bootstrap_uses_original_groups_and_is_reproducible(self):
        frame = pd.concat([toy(), toy().assign(id=lambda x: x.id + "_relative")], ignore_index=True)
        a = evaluate(frame)
        b = evaluate(frame)
        assert a["bootstrap"]["n_original_groups"] == 4
        assert a["matched"]["gain_enrichment"] == b["matched"]["gain_enrichment"]
        assert a["whole_cohort"]["policy_gain_vs_base"] == b["whole_cohort"]["policy_gain_vs_base"]


    def test_invalid_duplicate_ids_and_nonboolean_selection_rejected(self):
        frame = toy()
        frame.loc[1, "id"] = frame.loc[0, "id"]
        with self.assertRaisesRegex(ValueError, "unique"):
            evaluate(frame)
        frame = toy()
        frame["selected"] = ["True", "True", "False", "False"]
        with self.assertRaisesRegex(ValueError, "booleans"):
            evaluate(frame)

    def test_tiny_matched_sets_remain_descriptive_without_inferential_ci(self):
        result = evaluate(toy(), min_matched_pairs=20)
        effect = result["matched"]["gain_enrichment"]
        self.assertEqual(effect["estimate"], 3.)
        self.assertEqual(effect["status"], "insufficient_matched_pairs")
        self.assertIsNone(effect["ci_low"])
        self.assertIsNotNone(result["whole_cohort"]["policy_gain_vs_base"]["estimate"])

    def test_policy_vs_always_candidate_has_its_own_paired_bootstrap_contrast(self):
        frame = toy()
        frame.loc[~frame.selected, "p_candidate"] = 3.
        result = evaluate(frame)
        effect = result["whole_cohort"]["policy_gain_vs_candidate"]
        self.assertAlmostEqual(effect["estimate"], 2.5)
        self.assertEqual(effect["n_left"], 4)
        self.assertEqual(effect["bootstrap_valid_replicates"], 100)
        self.assertIsNotNone(effect["ci_low"])
        frame.loc[2, "p_candidate"] = np.nan
        result = evaluate(frame, n_boot=0)
        self.assertEqual(result["whole_cohort"]["policy_gain_vs_candidate"]["n_left"], 3)
        self.assertEqual(result["whole_cohort"]["policy_on_candidate_evaluable"]["n"], 3)
        self.assertAlmostEqual(result["whole_cohort"]["policy_gain_vs_candidate"]["estimate"], 5. / 3)

    def test_logloss_uses_same_clipping_as_head_training(self):
        y, p = np.array([1., 0., 1., np.nan]), np.array([0., 1., 1e-20, .5])
        actual = individual_loss(y, p, "logloss")
        expected = row_loss(y, p, task="binary", loss_name="logloss")
        np.testing.assert_allclose(actual, expected, equal_nan=True)
        self.assertAlmostEqual(actual[0], -np.log(1e-7))

    def test_outcome_blind_control_eligibility_restricts_pool_and_records_shortage(self):
        frame = toy()
        frame["control_eligible"] = [True, True, False, True]
        pairs, diagnostics = make_matched_controls(frame, control_eligible_col="control_eligible")
        self.assertEqual(list(pairs.control_id), ["c2"])
        self.assertEqual(diagnostics["n_excluded_unselected_control_ineligible"], 1)
        self.assertEqual(diagnostics["n_unmatched_selected"], 1)
        self.assertEqual(diagnostics["control_eligible_column"], "control_eligible")
        self.assertEqual(diagnostics["n_excluded_unselected_same_group"], 0)
        result = evaluate(frame, n_boot=0, control_eligible_col="control_eligible")
        self.assertEqual(result["matching"]["n_pairs"], 1)
        self.assertEqual(result["coverage"]["n_switched_query"], 2)
        frame["y"] = np.nan
        pairs2, _ = make_matched_controls(frame, control_eligible_col="control_eligible")
        self.assertEqual(pairs.to_dict("records"), pairs2.to_dict("records"))

    def test_many_pairs_in_two_families_do_not_create_matched_inference(self):
        frame = pd.concat([
            toy().assign(id=lambda x, k=k: x.id + str(k)) for k in range(10)
        ], ignore_index=True)
        frame["group_id"] = np.where(frame.selected, "selected_family", "control_family")
        result = evaluate(frame, min_matched_pairs=20, min_matched_groups=10)
        effect = result["matched"]["gain_enrichment"]
        self.assertEqual(result["matching"]["n_pairs"], 20)
        self.assertEqual(effect["n_unique_groups_left"], 1)
        self.assertEqual(effect["n_unique_groups_right"], 1)
        self.assertEqual(effect["status"], "insufficient_matched_groups")
        self.assertEqual(result["matching"]["gain_status"], "insufficient_matched_groups")
        self.assertIsNone(effect["ci_low"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
