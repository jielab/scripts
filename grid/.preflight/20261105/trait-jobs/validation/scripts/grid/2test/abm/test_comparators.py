"""Pooled fits must not depend on optional ancestry labels or an R process."""
import unittest
import subprocess
from unittest.mock import patch
import numpy as np
import pandas as pd
from abm_core.comparators import fit_comparators
from abm_core.pipeline import predictor_identity


class PooledComparatorTests(unittest.TestCase):
    def test_rare_or_missing_ancestry_does_not_partition_comparator_training(self):
        rng = np.random.default_rng(29)
        features = ["csx.AFR", "csx.EAS", "csx.EUR", "csx.SAS", "age"]
        frame = pd.DataFrame(rng.normal(size=(600, 5)), columns=features)
        frame["role"] = np.repeat(["build", "tune_model", "prediction_calibration"], [240,180,180])
        frame["y"] = (frame["csx.EUR"] + rng.normal(size=600) > 0).astype(float)
        frame["ancestry"] = "EUR"
        frame.loc[0, "ancestry"] = "EAS"
        config = dict(_binary=True, seed=29, logistic_c_grid=[.1], moe_clusters=[2],
                      moe_min_n=30, moe_min_class=10)
        original_run = subprocess.run
        def guard(command, *args, **kwargs):
            if any("Rscript" in str(part) or str(part).endswith(".R") for part in command):
                raise AssertionError("pooled comparators need no R program")
            return original_run(command, *args, **kwargs)
        with patch("subprocess.run", side_effect=guard):
            expected, _, status = fit_comparators(frame, features, [0,1,2,3], config)
            self.assertEqual(set(status.method), {"Clinical", "csx-moe"})
            self.assertTrue(status.status.eq("completed").all())
            for changed in (frame.assign(ancestry=None), frame.drop(columns="ancestry")):
                observed, _, new_status = fit_comparators(changed, features, [0,1,2,3], config)
                pd.testing.assert_frame_equal(status, new_status)
                for name in expected:
                    np.testing.assert_array_equal(expected[name], observed[name])

    def test_predictor_cache_identity_does_not_include_ancestry_labels(self):
        frame = pd.DataFrame(dict(id=["a","b","c"], group_id=["a","b","c"],
            role=["build","tune_model","prediction_calibration"], y=[0.,1.,0.],
            ancestry=["EUR","EAS","AFR"], score=[1.,2.,3.]))
        first = predictor_identity(frame, ["score"], [0], "moe", {}, None, {})
        other = predictor_identity(frame.drop(columns="ancestry"), ["score"], [0], "moe", {}, None, {})
        self.assertEqual(first, other)
