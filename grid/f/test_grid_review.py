"""Regression tests for evaluation, score lineage and frozen-feature contracts.

Synthetic data only. Run from the shared GRID environment with
python -m unittest discover -s f -p 'test_*review.py'.
"""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd
from scipy.special import expit

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("grid_metrics_review", ROOT / "grid.py")
grid = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = grid
spec.loader.exec_module(grid)
abm = grid.load_module("grid.abm")
data = grid.load_module("grid.data")


class EvaluationReview(unittest.TestCase):
    def test_inverse_prediction_is_not_high_predictive_r2(self):
        y = np.arange(-10., 11.)
        values = grid.prediction_metrics(y, -y, np.zeros(len(y)), False)
        self.assertAlmostEqual(values["residual_correlation_R2"], 1.)
        self.assertAlmostEqual(values["predictive_R2"], -3.)
        self.assertEqual(values["predictive_R2"], values["total_R2"])
        self.assertNotIn("Prediction_R2", values)
        self.assertAlmostEqual(values["calibration_slope"], -1.)

    def test_calibration_is_diagnostic_and_does_not_correct_error(self):
        p = np.arange(1., 11.)
        y = 4. + 2. * p
        original = p.copy()
        values = grid.prediction_metrics(y, p, np.zeros(len(y)), False)
        self.assertAlmostEqual(values["calibration_intercept"], 4.)
        self.assertAlmostEqual(values["calibration_slope"], 2.)
        self.assertEqual(values["MSE"], np.mean((y-p)**2))
        np.testing.assert_array_equal(original, p)
        self.assertFalse(values["calibration_refits_prediction"])

    def test_binary_calibration_with_exact_group_frequencies(self):
        p = np.repeat([.2, .5, .8], 100)
        y = np.concatenate([np.r_[np.ones(k), np.zeros(100-k)] for k in (20,50,80)])
        result = grid.calibration_diagnostics(y, p, True)
        self.assertEqual(result["calibration_status"], "estimated")
        self.assertAlmostEqual(result["calibration_intercept"], 0., places=5)
        self.assertAlmostEqual(result["calibration_slope"], 1., places=5)

    def test_constant_prediction_calibration_is_not_identifiable(self):
        result = grid.calibration_diagnostics(np.arange(8.), np.ones(8), False)
        self.assertEqual(result["calibration_status"], "not_identifiable")
        self.assertTrue(np.isnan(result["calibration_slope"]))

    def test_separated_binary_calibration_has_no_finite_mle(self):
        y = np.array([0,0,0,1,1,1])
        for p in (np.array([.1,.2,.3,.7,.8,.9]), np.array([.1,.2,.5,.5,.8,.9])):
            for outcome in (y, 1-y):
                result = grid.calibration_diagnostics(outcome,p,True)
                self.assertEqual(result["calibration_status"],"separated")
                self.assertTrue(np.isnan(result["calibration_slope"]))
                self.assertTrue(np.isnan(result["calibration_intercept"]))

    def test_linear_and_logit_baseline_terms_reconstruct(self):
        rng = np.random.default_rng(902)
        frame = pd.DataFrame({"age":rng.normal(size=500), "csx.EUR":rng.normal(size=500),
                              "ancestry":np.repeat(["EUR","EAS"],250),
                              "eid":[str(i) for i in range(500)], "family_id":[str(i) for i in range(500)]})
        eta = .2 + .3*frame.age + .7*frame["csx.EUR"]
        for kind, y in (("continuous", eta.to_numpy()), ("binary", rng.binomial(1, expit(eta)))):
            model = abm.AncestryPredictor(["age","csx.EUR"],kind,1.).fit(frame,y)
            terms = model.decompose(frame)
            linear = terms.baseline_intercept + terms.filter(like="baseline_contribution.").sum(axis=1)
            predicted = expit(linear) if kind == "binary" else linear
            np.testing.assert_allclose(predicted,model.predict(frame),atol=1e-9)
            self.assertTrue((terms.baseline_link == ("logit" if kind == "binary" else "identity")).all())

    def test_low_error_audit_uses_same_model_and_prespecified_mask(self):
        n = 160
        frame = pd.DataFrame({"eid":[str(i) for i in range(n)]})
        table = pd.DataFrame({"selected_absolute_error":np.arange(n)<80,
                              "GRID_evolution":np.zeros(n), "family_id":[str(i//2) for i in range(n)]})
        bundle = {"primary_arm":"GRID_evolution", "config":{**abm.DEFAULTS,"audit_bootstrap":100}}
        y = np.r_[np.ones(80)*.1,np.ones(80)*2]
        with mock.patch.object(abm,"_raw_predictions",return_value=(table,{})), \
             mock.patch.object(abm,"_apply_policy",side_effect=lambda b,r:r):
            audit = abm._internal_low_error_audit(bundle,frame,y)
        self.assertEqual(audit["status"],"supported_lower_group_error")
        self.assertAlmostEqual(audit["selected_loss"],.01)
        self.assertAlmostEqual(audit["all_loss"],2.005)
        self.assertLess(audit["upper"],0.)
        self.assertFalse(audit["test_outcomes_used"])
        self.assertFalse(audit["changes_prediction_policy"])


class SourceContracts(unittest.TestCase):
    def test_missing_score_lineage_fails_before_loading_phenotypes(self):
        with tempfile.TemporaryDirectory() as temporary:
            argv = ["grid.py", "--stage", "check", "--score-dir", str(Path(temporary) / "absent-scores")]
            with mock.patch.object(sys, "argv", argv), \
                 mock.patch.object(grid, "load_module", side_effect=AssertionError("Data must not be read")):
                with self.assertRaisesRegex(ValueError, "1csx.scores.provenance.json"):
                    grid.main()

    def test_prsformer_comparison_checks_covariate_values_by_id(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/"prsformer").mkdir()
            (root/"scores/height").mkdir(parents=True)
            roster_path = root/"prsformer/3.prsformer.split.rds"
            score_path = root/"scores/height/3.prsformer.scores.rds"
            roster_path.write_bytes(b"synthetic roster identity")
            score_path.write_bytes(b"synthetic score identity")
            cohort = pd.DataFrame({"eid":[str(i) for i in range(6)],
                "family_id":[str(i) for i in range(6)], "split":["train"]*3+["test"]*3})
            roster = cohort[["eid","split"]].copy()
            roster.loc[2,"split"] = "validation"
            features = pd.DataFrame({"eid":cohort.eid, "age":np.arange(40.,46.), "sex":[0,1]*3})
            table = pd.DataFrame({"eid":["3","4","5"],"y":[170.,180.,175.]})
            scores = pd.DataFrame({"eid":["5","3","4"], "split":"test",
                "prediction":[174.,171.,178.], "outcome":[175.,170.,180.],
                "covariate_names":"age,sex", "covariate.age":[45.,43.,44.],
                "covariate.sex":[1,1,0], "endpoint_definition":"height"})
            args = SimpleNamespace(prsformer_root=str(root), covariates="age,sex")
            result = {"metadata":{"endpoint":"height"}}
            def read(path):
                return roster.copy() if Path(path) == roster_path else scores.copy()
            with mock.patch.object(grid,"named_table",side_effect=read):
                combined, status = grid.add_prsformer(args,"height",result,cohort,table,features)
                np.testing.assert_array_equal(combined.PRSformer,[171.,178.,174.])
                self.assertEqual(status["covariate_values"],"matched_by_eid")
                scores.loc[0,"covariate.age"] = 20.
                with self.assertRaisesRegex(ValueError,"covariate values differ for age"):
                    grid.add_prsformer(args,"height",result,cohort,table,features)
                scores.drop(columns="covariate.age",inplace=True)
                with self.assertRaisesRegex(ValueError,"original covariate values"):
                    grid.add_prsformer(args,"height",result,cohort,table,features)
                # Both methods support an explicitly empty covariate model.
                # The sentinel must not be treated as a column named "none".
                args.covariates = "none"
                scores["covariate_names"] = "none"
                combined, status = grid.add_prsformer(args,"height",result,cohort,table)
                np.testing.assert_array_equal(combined.PRSformer,[171.,178.,174.])

    def test_feature_manifest_rejects_wrong_contract_or_changed_matrix(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root/"features.tsv"
            source.write_text("eid\tcsx.EUR\n1\t.2\n")
            manifest = root/"features.manifest.json"
            body = {"format":"GRID-features-1", "trait":"height",
                    "feature_contract_sha256":"frozen-contract", "data_sha256":grid.sha256(source)}
            manifest.write_text(json.dumps(body))
            metadata = {"trait":"height", "feature_contract":{"contract_sha256":"frozen-contract"}}
            result = grid.validate_feature_manifest(manifest,source,metadata)
            self.assertEqual(result["status"],"provider_declared_contract_and_file_hash_matched")
            body["feature_contract_sha256"] = "different-weights"
            manifest.write_text(json.dumps(body))
            with self.assertRaisesRegex(ValueError,"different weights"):
                grid.validate_feature_manifest(manifest,source,metadata)
            body["feature_contract_sha256"] = "frozen-contract"
            manifest.write_text(json.dumps(body))
            source.write_text("eid\tcsx.EUR\n1\t.7\n")
            with self.assertRaisesRegex(ValueError,"exact data file"):
                grid.validate_feature_manifest(manifest,source,metadata)

    def test_model_runtime_is_checked_before_deserialization(self):
        metadata = {"sources":grid.source_signature(),"runtime":{"sklearn":"other"}}
        saved = {"format":"GRID-model-1","metadata_json":json.dumps(metadata),
                 "python_payload_codec":"base64-zlib-pickle5","python_payload":"invalid"}
        with mock.patch.object(grid,"read_rds",return_value=saved):
            with self.assertRaisesRegex(ValueError,"runtime differs"):
                grid.read_model("unused")

    def test_native_csx_weight_lineage_and_chromosomes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/"height").mkdir()
            score = root/"height/1csx.scores.rds"
            score.write_bytes(b"synthetic score-file identity")
            models,registry = {},{}
            for pop in data.POPS:
                weight = root/f"{pop}.csx.gz"
                weight.write_bytes(f"weight-{pop}".encode())
                registry[f"csx.{pop}"] = {"weight":weight,"signature":"joint-run"}
                models[f"csx.{pop}"] = {"weights_sha256":grid.sha256(weight),
                    "joint_inference_signature":"joint-run","chromosomes":[1,2],
                    "scoring_convention":"PLINK2_SCORESUM_no_mean_imputation"}
            manifest = score.with_suffix(".provenance.json")
            manifest.write_text(json.dumps({"schema":"grid_csx_scores_v1", "models":models,
                                            "score_file_sha256":grid.sha256(score)}))
            args = SimpleNamespace(score_dir=str(root))
            self.assertEqual(data.verify_csx_score_provenance(args,"height",registry,[1,2]),manifest)
            with self.assertRaisesRegex(ValueError,"chromosome sets differ"):
                data.verify_csx_score_provenance(args,"height",registry,[1])
            registry["csx.EUR"]["weight"].write_bytes(b"changed weights")
            with self.assertRaisesRegex(ValueError,"different CSx runs"):
                data.verify_csx_score_provenance(args,"height",registry,[1,2])


if __name__ == "__main__":
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        unittest.main(verbosity=2)
