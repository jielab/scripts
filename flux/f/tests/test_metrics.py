"""Numerical and leakage checks for the evaluation and baseline modules.

Run: python -m unittest discover -s f/tests -p test_metrics.py -v
These tests use synthetic outcomes only and never require UKB files.
"""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from baseline import RidgeCoxPH, TabularPreprocessor, fit_predict_baselines
from metrics import (CensoringKM, binary_metrics, calibration_bins, evaluate_endpoint,
                     evaluate_predictions, harrell_c_index, ipcw_weights, macro_average)


def brute_concordance(time, event, score):
    comparable, correct, tied = 0, 0, 0
    for i in range(len(time)):
        if not event[i]:
            continue
        for j in range(len(time)):
            if time[i] < time[j] or (time[i] == time[j] and not event[j]):
                comparable += 1
                correct += score[i] > score[j]
                tied += score[i] == score[j]
    return ((correct + 0.5 * tied) / comparable if comparable else np.nan,
            comparable, correct, tied)


class MetricTests(unittest.TestCase):
    def test_uncensored_matches_sklearn_and_ordinary_brier(self):
        p = np.array([0.9, 0.4, 0.6, 0.1])
        y = np.array([1, 1, 0, 0], dtype=bool)
        result = evaluate_endpoint([1, 2, np.inf, 8], [12] * 4, p,
                                   ref_time=[1, 2, 4, 6, 10, 12], ref_event=[True] * 6,
                                   probabilities=p, probability_horizon_years=5)
        self.assertAlmostEqual(result["recorded_auroc"], roc_auc_score(y, p))
        self.assertAlmostEqual(result["recorded_ap"], average_precision_score(y, p))
        self.assertAlmostEqual(result["ipcw_auroc"], roc_auc_score(y, p))
        self.assertAlmostEqual(result["ipcw_brier"], np.mean((y - p) ** 2))
        self.assertEqual(result["ipcw_status"], "ok")
        self.assertEqual(result["probability_status"], "ok")

    def test_recorded_and_censoring_estimands_are_separate(self):
        result = evaluate_endpoint([1, 4, np.inf], [10, 2, 10], [0.9, 0.8, 0.1],
                                   ref_time=[1, 2, 8, 12], ref_event=[True] * 4)
        self.assertEqual(result["recorded_cases"], 2)
        self.assertEqual(result["ipcw_cases"], 1)
        self.assertEqual(result["ipcw_controls"], 1)

    def test_reverse_km_ties_and_left_limit_weights(self):
        km = CensoringKM.fit([2, 2, 4, 8], [True, False, False, True])
        np.testing.assert_allclose(km.predict([1, 2, 4, 5]), [1, 2 / 3, 1 / 3, 1 / 3])
        np.testing.assert_allclose(km.predict([2, 4], left_limit=True), [1, 2 / 3])
        case, control, weight, info = ipcw_weights([2, 4, 8, 1], [True, True, False, False], 5, km)
        np.testing.assert_array_equal(case, [True, True, False, False])
        np.testing.assert_array_equal(control, [False, False, True, False])
        np.testing.assert_allclose(weight, [1, 1.5, 3, 0])
        self.assertEqual(info["status"], "ok")
        self.assertTrue(np.isnan(km.predict([9])[0]))

    def test_missing_support_is_not_silently_clipped(self):
        km = CensoringKM.fit([1, 3, 5], [False] * 3)
        _, _, weights, info = ipcw_weights([2, 6], [True, False], 5, km)
        self.assertEqual(info["status"], "insufficient_censoring_support")
        np.testing.assert_array_equal(weights, [0, 0])
        _, _, _, info = ipcw_weights([2, 6], [True, False], 6, km)
        self.assertEqual(info["status"], "horizon_outside_development_followup")

    def test_positive_extreme_weights_are_reported_and_truncated(self):
        km = CensoringKM.fit([1] * 19 + [10], [False] * 19 + [True])
        _, _, weight, info = ipcw_weights([3, 10], [True, False], 5, km, weight_cap=10)
        np.testing.assert_allclose(weight, [10, 10])
        self.assertEqual(info["weights_truncated"], 2)
        _, _, untruncated, _ = ipcw_weights([3, 10], [True, False], 5, km, weight_cap=None)
        np.testing.assert_allclose(untruncated, [20, 20])

    def test_no_implicit_test_censoring_and_no_ten_year_probability(self):
        result = evaluate_endpoint([1, np.inf], [12, 12], [0.8, 0.2])
        self.assertEqual(result["ipcw_status"], "missing_development_censoring_reference")
        self.assertTrue(np.isnan(result["ipcw_auroc"]))
        result = evaluate_endpoint([1, np.inf], [12, 12], [0.8, 0.2], horizon_years=10,
                                   ref_time=[1, 12], ref_event=[True, True],
                                   probabilities=[0.8, 0.2], probability_horizon_years=5)
        self.assertEqual(result["ipcw_auroc"], 1)
        self.assertEqual(result["probability_status"], "different_or_unknown_probability_horizon")
        self.assertTrue(np.isnan(result["ipcw_brier"]))

    def test_no_events_or_no_controls_remain_undefined(self):
        self.assertTrue(np.isnan(binary_metrics([0, 0], [0.1, 0.2])["average_precision"]))
        self.assertTrue(np.isnan(harrell_c_index([1, 2], [False, False], [0.1, 0.2])))
        result = evaluate_endpoint([np.inf, np.inf], [12, 12], [0.1, 0.2],
                                   ref_time=[2, 12], ref_event=[True, True],
                                   probabilities=[0.1, 0.2], probability_horizon_years=5)
        self.assertEqual(result["ipcw_status"], "no_cases_or_controls_at_horizon")
        self.assertTrue(np.isnan(result["ipcw_auroc"]))
        self.assertAlmostEqual(result["ipcw_brier"], 0.025)

    def test_prevalence_and_explicit_risk_mask(self):
        result = evaluate_endpoint([-1, 0, 1, np.inf, 2], [12] * 5, [0.2] * 5,
                                   risk_mask=[True, True, True, True, False])
        self.assertEqual(result["n"], 2)
        self.assertEqual(result["excluded"], 3)
        self.assertEqual(result["recorded_cases"], 1)

    def test_fenwick_c_index_matches_independent_pairwise_definition(self):
        rng = np.random.default_rng(721)
        for _ in range(20):
            time = rng.integers(1, 9, 35)
            event = rng.uniform(size=35) < 0.6
            score = rng.integers(0, 6, 35)
            expected = brute_concordance(time, event, score)
            found = harrell_c_index(time, event, score, return_details=True)
            self.assertAlmostEqual(found["c_index"], expected[0])
            self.assertEqual(found["comparable"], expected[1])
            self.assertEqual(found["concordant"], expected[2])
            self.assertEqual(found["tied_scores"], expected[3])

    def test_calibration_is_descriptive_with_zero_weight_early_censor(self):
        km = CensoringKM.fit([2, 2, 4, 8], [True, False, False, True])
        rows = calibration_bins([0.2, 0.3, 0.4, 0.5], [2, 4, 8, 1],
                                [True, True, False, False], horizon_years=5,
                                reference=km, bins=1)
        self.assertEqual(rows[0]["n"], 4)
        self.assertAlmostEqual(rows[0]["predicted_mean"], 0.35)
        self.assertAlmostEqual(rows[0]["ipcw_observed"], 2.5 / 5.5)
        self.assertEqual(rows[0]["observed_controls"], 1)

    def test_prediction_table_death_censoring_and_horizon_contract(self):
        # Times are days after a shared epoch; nonfatal death cuts follow-up.
        table = pd.DataFrame({"eid": [1, 2, 3, 4], "split": ["train", "train", "test", "test"],
                              "endpoint": ["I21"] * 4, "score": [0.8, 0.1, 0.8, 0.1],
                              "probability": [0.8, 0.1, 0.8, 0.1], "prediction_horizon_years": [5] * 4,
                              "landmark_day": [0] * 4, "event_day": [365.25, np.nan, 365.25, np.nan],
                              "exit_day": [12 * 365.25] * 4, "death_day": [np.nan] * 4})
        metrics, calibration = evaluate_predictions(table, horizons=(5, 10))
        self.assertEqual(len(metrics), 2)
        self.assertEqual(metrics.iloc[0]["censor_reference_n"], 2)
        self.assertEqual(metrics.iloc[0]["probability_status"], "ok")
        self.assertEqual(metrics.iloc[1]["probability_status"], "different_or_unknown_probability_horizon")
        self.assertEqual(set(calibration["horizon_years"]), {5})
        table.loc[3, "death_day"] = 2 * 365.25
        metrics, _ = evaluate_predictions(table, horizons=(5,))
        self.assertEqual(metrics.iloc[0]["recorded_controls"], 1)
        self.assertEqual(metrics.iloc[0]["ipcw_controls"], 0)
        table["endpoint"] = "death"
        metrics, _ = evaluate_predictions(table, horizons=(5,))
        self.assertEqual(metrics.iloc[0]["ipcw_controls"], 1)
        with self.assertRaisesRegex(ValueError, "training rows"):
            evaluate_predictions(table, reference=table[table.split == "test"])

    def test_macro_reports_denominators(self):
        rows = pd.DataFrame({"model": ["flux"] * 2, "split": ["test"] * 2,
                             "horizon_years": [5] * 2, "recorded_auroc": [0.7, np.nan],
                             "ipcw_auroc": [0.6, 0.8]})
        macro = macro_average(rows).iloc[0]
        self.assertEqual(macro["endpoints"], 2)
        self.assertEqual(macro["recorded_auroc_endpoints"], 1)
        self.assertAlmostEqual(macro["recorded_auroc"], 0.7)
        self.assertEqual(macro["ipcw_auroc_endpoints"], 2)
        self.assertAlmostEqual(macro["ipcw_auroc"], 0.7)


class BaselineTests(unittest.TestCase):
    def test_preprocessing_uses_training_statistics_and_preserves_inputs(self):
        train = np.array([[1, np.nan, 10], [3, np.nan, np.nan], [np.nan, np.nan, 14]], dtype=float)
        original = train.copy()
        preprocessor = TabularPreprocessor().fit(train)
        np.testing.assert_allclose(train, original, equal_nan=True)
        np.testing.assert_allclose(preprocessor.median_, [2, 0, 12])
        test = np.array([[100, np.nan, np.nan], [np.nan, 10, 20]])
        transformed = preprocessor.transform(test)
        self.assertEqual(transformed.dtype, np.float32)
        self.assertEqual(transformed[0, 1], 0)
        self.assertEqual(transformed[0, 2], 0)
        self.assertEqual(transformed[1, 0], 0)
        self.assertGreater(transformed[0, 0], 50)
        np.testing.assert_allclose(preprocessor.median_, [2, 0, 12])

    def test_cox_objective_matches_partial_likelihood_and_finite_difference(self):
        rng = np.random.default_rng(81)
        x = rng.normal(size=(7, 3))
        time = np.array([7, 5, 5, 3, 2, 2, 1])
        event = np.array([1, 1, 0, 1, 0, 1, 1], dtype=float)
        ends = np.array([0, 2, 3, 5, 6])
        counts = np.array([1, 1, 1, 1, 1])
        beta = np.array([0.2, -0.3, 0.1])
        alpha = 0.07
        value, gradient = RidgeCoxPH._objective(beta, x, event, ends, counts, alpha)
        eta = x @ beta
        naive = sum(np.log(np.exp(eta[time >= time[i]]).sum()) - eta[i]
                    for i in range(len(time)) if event[i]) / event.sum() + alpha / 2 * (beta @ beta)
        self.assertAlmostEqual(value, naive, places=12)
        epsilon = 1e-6
        finite = []
        for j in range(len(beta)):
            plus, minus = beta.copy(), beta.copy()
            plus[j] += epsilon
            minus[j] -= epsilon
            finite.append((RidgeCoxPH._objective(plus, x, event, ends, counts, alpha)[0]
                           - RidgeCoxPH._objective(minus, x, event, ends, counts, alpha)[0]) / (2 * epsilon))
        np.testing.assert_allclose(gradient, finite, atol=1e-8, rtol=1e-6)

    def test_cox_recovers_risk_direction_and_monotone_horizon_probability(self):
        rng = np.random.default_rng(98)
        x = rng.normal(size=(700, 3)).astype(np.float32)
        failure = rng.exponential(5 * np.exp(-1.2 * x[:, 0]))
        censor = rng.uniform(3, 15, len(x))
        observed, event = np.minimum(failure, censor), failure <= censor
        model = RidgeCoxPH(alpha=0.01).fit(x, observed, event)
        self.assertGreater(model.coef_[0], 0.8)
        self.assertLess(abs(model.coef_[1]), 0.25)
        self.assertGreater(harrell_c_index(observed, event, model.predict(x)), 0.7)
        p1, p5, p10 = [model.predict_proba(x, horizon) for horizon in (1, 5, 10)]
        self.assertTrue(np.all(p1 <= p5))
        self.assertTrue(np.all(p5 <= p10))
        self.assertTrue(np.all((p10 >= 0) & (p10 <= 1)))

    @staticmethod
    def synthetic_inputs():
        rng = np.random.default_rng(41)
        x = rng.normal(size=(240, 4)).astype(np.float32)
        failure = rng.exponential(5 * np.exp(-1.5 * x[:, 0]))
        observed, event = np.minimum(failure, 12), failure <= 12
        labels = np.column_stack([failure <= 5, np.zeros(len(x), dtype=bool)])
        mask = np.ones_like(labels)
        time = np.column_stack([observed, np.full(len(x), 12)])
        status = np.column_stack([event, np.zeros(len(x), dtype=bool)])
        return x, labels, mask, time, status

    def test_baselines_fit_serialize_and_explicitly_skip_rare_endpoint(self):
        x, labels, mask, time, event = self.synthetic_inputs()
        with tempfile.TemporaryDirectory() as directory:
            result = fit_predict_baselines(x[:150], labels[:150], mask[:150],
                                           x[150:200], labels[150:200], mask[150:200], x[200:],
                                           methods=("logistic", "histgb", "cox"),
                                           survival_train=(time[:150], event[:150]),
                                           survival_val=(time[150:200], event[150:200]),
                                           endpoint_names=["I21", "rare"], output_dir=directory,
                                           options={"logistic": {"C": [0.1, 1]},
                                                    "histgb": {"max_leaf_nodes": [7]},
                                                    "cox": {"alpha": [0.01, 0.1]}},
                                           min_cases=5, max_iter=70)
            self.assertNotIn("models", result)
            self.assertEqual(len(result["status"]), 6)
            self.assertTrue((Path(directory) / "preprocess.joblib").exists())
            preprocessor = joblib.load(Path(directory) / "preprocess.joblib")
            test_x = preprocessor.transform(x[200:])
            for status in result["status"]:
                scores = result["predictions"][status["method"]]["test"]
                self.assertIsInstance(scores, np.memmap)
                if status["endpoint"] == "rare":
                    self.assertEqual(status["status"], "insufficient_training_cases_or_controls")
                    self.assertTrue(np.isnan(scores[:, 1]).all())
                    continue
                self.assertIn(status["status"], ("ok", "optimizer_not_converged", "optimizer_iteration_limit"))
                self.assertTrue(np.isfinite(scores[:, 0]).all())
                self.assertGreater(roc_auc_score(labels[200:, 0], scores[:, 0]), 0.65)
                persisted = joblib.load(status["model_path"])
                estimator = persisted["estimator"]
                reloaded = (estimator.predict_proba(test_x, 5) if status["method"] == "cox"
                            else estimator.predict_proba(test_x)[:, 1])
                np.testing.assert_allclose(reloaded, scores[:, 0], rtol=1e-6, atol=1e-7)

    def test_test_covariates_cannot_change_training_or_selection(self):
        x, labels, mask, _, _ = self.synthetic_inputs()
        arguments = (x[:150], labels[:150, :1], mask[:150, :1],
                     x[150:200], labels[150:200, :1], mask[150:200, :1])
        options = {"methods": ("logistic",), "return_models": True, "min_cases": 5}
        original = fit_predict_baselines(*arguments, x[200:], **options)
        changed = fit_predict_baselines(*arguments, 1000 * x[200:], **options)
        np.testing.assert_array_equal(original["preprocessor"].median_, changed["preprocessor"].median_)
        np.testing.assert_allclose(original["predictions"]["logistic"]["validation"],
                                   changed["predictions"]["logistic"]["validation"])
        np.testing.assert_allclose(original["models"]["logistic"]["0"].coef_,
                                   changed["models"]["logistic"]["0"].coef_)
        self.assertEqual(original["status"][0]["selected_parameter"], changed["status"][0]["selected_parameter"])


if __name__ == "__main__":
    unittest.main()
