import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from abm_core.grid_inputs import partition_calibration, validate_lineage, IndexedBlocks
from abm_core.predictors import fit_predictors, predict
from abm_core.model import ABMAssessor
from abm_core.representation import BlockGeometry
from abm_core.pipeline import Stages
from abm_core.reporting import (freeze_matches, evaluate, fit_diagnostics, diagnostic_predictions,
                             evaluate_diagnostics, evaluation_calibration, plot_calibration)


def fixture(binary=False):
    rng = np.random.default_rng(9)
    roles = np.repeat(["build", "tune_model", "tune_gate", "prediction_calibration",
                       "head_calibration", "policy_calibration", "audit", "test"], 220)
    x = rng.normal(size=(len(roles), 18))
    y = 2 * x[:, 0] + x[:, 14] + rng.normal(size=len(x))
    if binary:y = (y > np.median(y)).astype(float)
    features = [f"side{i}" for i in range(18)]
    d = pd.DataFrame(x, columns=features)
    d["id"] = [f"i{i:05d}" for i in range(len(x))]
    d["group_id"] = [f"g{i:05d}" for i in range(len(x))]
    d["role"] = roles
    d["y"] = y
    d["ancestry"] = "EUR"
    config = {"_binary": binary, "seed": 7, "moe_clusters": [4], "moe_min_n": 30}
    return d, features, config


class AdapterTests(unittest.TestCase):
    def test_role_partition_is_group_and_order_invariant(self):
        d = pd.DataFrame({"id": [str(i) for i in range(1200)],
                          "group_id": [str(i//2) for i in range(1200)], "role": "calibration_fit"})
        a = partition_calibration(d, 42).set_index("id")
        b = partition_calibration(d.sample(frac=1, random_state=4), 42).set_index("id").loc[a.index]
        pd.testing.assert_series_equal(a.role, b.role)
        self.assertEqual(a.groupby("group_id").role.nunique().max(), 1)
        self.assertEqual(set(a.role), {"prediction_calibration", "head_calibration", "policy_calibration"})

    def test_cross_role_family_is_rejected(self):
        d, _, _ = fixture()
        validate_lineage(d)
        d.loc[230, "group_id"] = d.loc[0, "group_id"]
        with self.assertRaisesRegex(ValueError, "group"):
            validate_lineage(d)

    def test_predictor_outcome_lineage_and_recomputation(self):
        d, features, cfg = fixture()
        a, pa = fit_predictors(d, features, list(range(4)), "moe", cfg)
        perturbed = d.copy()
        unseen = d.role.isin(["tune_gate", "head_calibration", "policy_calibration", "audit", "test"])
        perturbed.loc[unseen, "y"] = np.nan
        b, pb = fit_predictors(perturbed, features, list(range(4)), "moe", cfg)
        np.testing.assert_allclose(pa.p_base, pb.p_base, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(pa.p_candidate, pb.p_candidate, rtol=1e-12, atol=1e-12)
        p0, p1 = predict(a, d.drop(columns="y"))
        np.testing.assert_array_equal(p0, pa.p_base)
        np.testing.assert_array_equal(p1, pa.p_candidate)
        self.assertEqual(a["scale"].n_features_in_, 18)
        self.assertEqual(a["lineage"]["calibrator_role"], "prediction_calibration")
        self.assertEqual(a["lineage"]["selection_role"], "tune_model")

    def test_binary_predictions_use_probabilities_and_all_18_columns(self):
        d, features, cfg = fixture(True)
        b, p = fit_predictors(d, features, list(range(4)), "moe", cfg)
        self.assertTrue(p.p_base.between(0, 1).all())
        self.assertTrue(p.p_candidate.between(0, 1).all())
        self.assertEqual(b["base"].n_features_in_, 18)
        self.assertEqual(b["training_objectives"]["candidate"], "logloss")
        self.assertEqual(b["training_objectives"]["selection"], "brier")
        d.loc[0, features[-1]] = np.nan
        with self.assertRaisesRegex(ValueError, "no silent"):
            fit_predictors(d, features, list(range(4)), "moe", cfg)

    def test_memmapped_geometry_matches_in_memory_oracle(self):
        rng = np.random.default_rng(3)
        x = rng.normal(size=(70, 12, 1))
        rows = np.arange(0, 70, 2)
        maps = np.repeat(["AFR", "EAS", "EUR", "SAS"], 3)
        indexed = IndexedBlocks(x, rows)
        a = BlockGeometry.fit(indexed, maps, row_batch=9)
        b = BlockGeometry.fit(x[rows].reshape(len(rows), -1), maps, row_batch=9)
        self.assertEqual(a.fingerprint(), b.fingerprint())
        with tempfile.TemporaryDirectory() as tmp:
            a.transform_to_npy(indexed, Path(tmp)/"x.npy", dtype="float64")
            np.testing.assert_array_equal(np.load(Path(tmp)/"x.npy"), b.transform(x[rows, :, 0]))

    def test_phase_receipts_reject_artifact_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stages = Stages(root, {"source": "frozen"})
            def prepare():
                (root/"a.txt").write_text("one")
                return [root/"a.txt"]
            stages.perform("prepare", prepare, "cpu", 1)
            stages.perform("prepare", lambda: self.fail("must resume"), "cpu", 1)
            (root/"a.txt").write_text("two")
            with self.assertRaisesRegex(ValueError, "Changed phase"):
                stages.perform("prepare", prepare, "cpu", 1)

    def test_unsupported_development_yields_na_and_baseline_policy(self):
        d, features, cfg = fixture()
        _, p = fit_predictors(d, features, list(range(4)), "moe", cfg)
        d = d[["id", "group_id", "role", "y", "ancestry"]].merge(
            p[["id", "p_base", "p_candidate"]], on="id", validate="one_to_one")
        d["n_twins"] = 0
        d["supported"] = False
        manifest = {"baseline_id": "csx-linear", "candidate_id": "csx-moe",
                    **{k: "a"*64 for k in ("bank_sha256", "geometry_sha256",
                                           "predictor_bundle_sha256", "role_manifest_sha256")}}
        m = ABMAssessor.fit(d[d.role.eq("tune_gate")], d[d.role.eq("head_calibration")],
            d[d.role.eq("policy_calibration")].drop(columns="y"), manifest=manifest)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"m.json";m.save(path);m=ABMAssessor.load(path)
            q=m.assess(d[d.role.eq("test")].drop(columns="y"),manifest=manifest)
        self.assertFalse(q.assessment_available.any())
        self.assertTrue(q.expected_loss_base.isna().all())
        np.testing.assert_array_equal(q.p_policy, q.p_base)

    def test_frozen_diagnostics_and_observed_calibration_report(self):
        d, features, cfg = fixture()
        _, p = fit_predictors(d, features, list(range(4)), "moe", cfg)
        d = d.merge(p[["id", "p_base", "p_candidate"]], on="id", validate="one_to_one")
        d["n_twins"] = np.random.default_rng(21).integers(30, 300, len(d))
        d["supported"] = True
        gate = d[d.role.eq("tune_gate")]
        head = d[d.role.eq("head_calibration")]
        policy = d[d.role.eq("policy_calibration")].drop(columns="y")
        manifest = {"baseline_id": "csx-linear", "candidate_id": "csx-moe",
                    **{k: "b"*64 for k in ("bank_sha256", "geometry_sha256",
                                           "predictor_bundle_sha256", "role_manifest_sha256")}}
        model = ABMAssessor.fit(gate, head, policy, manifest=manifest,
                               head_options={"n_bins": 4, "n_boot": 30})
        diagnostics = fit_diagnostics(gate, head, features[4:], "continuous", policy)
        test = d[d.role.eq("test")]
        assessed = model.assess(test.drop(columns="y"), manifest=manifest)
        frozen = diagnostic_predictions(diagnostics, assessed)
        perturbed = diagnostic_predictions(diagnostics, assessed.assign(y=1e10))
        cols = [c for c in frozen if c.startswith("diagnostic_")]
        pd.testing.assert_frame_equal(frozen[cols], perturbed[cols])
        evaluated = frozen.merge(test[["id", "y"]], on="id")
        evaluated.loc[evaluated.index[:5], "y"] = np.nan
        options = {"evaluation": {"bootstrap_replicates": 30},
                   "head_options": {"seed": 2, "min_calibration_groups": 20}}
        table = evaluation_calibration(evaluated, model, options, "test")
        row = table.loc[table.stratum.eq("all") & table.objective.eq("expected_loss") &
                        table.n_observed.gt(0)].iloc[0]
        subset = evaluated.loc[evaluated.count_bin.eq(row.count_bin) &
                               evaluated.assessment_available & evaluated.y.notna()]
        self.assertAlmostEqual(row.observed_estimate, float(((subset.y - subset.p_base)**2).mean()))
        self.assertEqual(row.n_observed, len(subset))
        report = evaluate_diagnostics(diagnostics, evaluated, options)
        self.assertEqual(len(report), 8)
        self.assertTrue((report.selection_threshold_role == "policy_calibration_X_only").all())
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration.png"
            plot_calibration(table, path)
            self.assertGreater(path.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
