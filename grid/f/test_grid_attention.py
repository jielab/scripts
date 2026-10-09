"""Real CUDA regression tests for GRID's selective-attention ABM; synthetic data only."""
from pathlib import Path
import importlib.util
import pickle
import sys
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("grid_grid_abm", HERE / "grid.abm.py")
abm = sys.modules.get(spec.name)
if abm is None:
    abm = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = abm
    spec.loader.exec_module(abm)
attention = abm._attention()


class DeviceContract(unittest.TestCase):
    def test_no_silent_cpu_fallback(self):
        with patch.object(torch.cuda, "is_available", return_value=False):
            with self.assertRaisesRegex(RuntimeError, "refuses CPU fallback"):
                attention.execution_device("cuda")

    def test_explicit_cpu_smoke(self):
        self.assertEqual(attention.preflight("cpu")["forward_backward"], "passed")


@unittest.skipUnless(torch.cuda.is_available(), "CUDA required for production backend tests")
class CudaAttention(unittest.TestCase):
    def test_production_model_forward_backward(self):
        info = attention.preflight("cuda")
        self.assertEqual(info["forward_backward"], "passed")
        self.assertTrue(info["device"].startswith("cuda"))

    def test_blocked_retrieval_matches_numpy_with_ties_families_and_folds(self):
        rng = np.random.default_rng(16)
        bank = rng.integers(-2, 3, (43, 4)).astype(float)
        bank[4:14] = bank[3]  # ties must obey hash order across block boundaries
        query = np.r_[bank[3:4], bank[:5], rng.normal(size=(5, 4))]
        ids = np.array([f"I{i}" for i in range(len(bank))])
        families = np.array([f"F{i//2}" for i in range(len(bank))])
        qids = np.r_[ids[:6], [f"Q{i}" for i in range(5)]]
        qfam = np.r_[families[:6], ["unknown"]*5]
        folds, qfolds = np.arange(43)%3, np.arange(len(query))%3
        tie = rng.permutation(len(bank))
        actual_j, actual_d, actual_n = attention.neighbors(bank, query, ids, qids, families,
            qfam, tie, 8, "cuda", 3, 7, folds, qfolds, radius=3.)
        distance = np.linalg.norm(query[:, None] - bank[None], axis=-1)
        allowed = (ids[None] != qids[:, None]) & (families[None] != qfam[:, None]) & (folds[None] != qfolds[:, None])
        distance[~allowed] = np.inf
        expected = np.array([np.lexsort((tie, row))[:8] for row in distance])
        np.testing.assert_array_equal(actual_j, expected)
        np.testing.assert_allclose(actual_d, np.take_along_axis(distance, expected, axis=1), atol=1e-12)
        np.testing.assert_array_equal(actual_n, (allowed & (distance <= 3)).sum(1))

    def test_masked_attention_empty_and_partial_neighborhoods(self):
        torch.manual_seed(81)
        model = attention.RowEncoder([0, 0, 1, 1], width=16, heads=2, layers=1, dropout=0).cuda()
        query = torch.randn(2, 16, device="cuda", requires_grad=True)
        donors = torch.randn(2, 4, 16, device="cuda")
        allowed = torch.tensor([[False]*4, [True, False, True, False]], device="cuda")
        w = model.donor_weights(query, donors, torch.zeros(2, 4, device="cuda"), allowed)
        np.testing.assert_array_equal(w[0].detach().cpu(), np.zeros(4))
        self.assertAlmostEqual(float(w[1].sum().detach()), 1., places=6)
        w.square().sum().backward()
        self.assertTrue(torch.isfinite(query.grad).all())

    def test_train_freeze_reload_and_test_outcome_isolation_both_endpoints(self):
        rng = np.random.default_rng(145)
        n = 640
        x, evolution = rng.normal(size=(2, n))
        frame = pd.DataFrame(dict(eid=[f"I{i:04}" for i in range(n)],
            family_id=[f"F{i//2:04}" for i in range(n)], ancestry="EUR",
            split=np.where(np.arange(n)<n//2, "train", "test"),
            x=x, evo=evolution, pc=rng.normal(size=n), age=rng.normal(size=n)))
        config = dict(abm_backend="selective_attention", device="cuda", retrieval="cuda",
            attention_epochs=2, attention_width=16, attention_heads=2, attention_layers=1,
            attention_dropout=0., attention_batch=64, donor_block=73,
            k_grid=[8], alpha_grid=[0.,1.], radius_grid=[4.], folds=3, min_role_n=8,
            hgb_iterations=3, gate_iterations=3, hgb_leaves_grid=[7], csx_alpha_grid=[1.],
            ridge_alpha_grid=[1.], audit_bootstrap=20, verbose=False)
        groups = dict(csx=["x"], ancestry=["pc"], covariates=["age"], evolution=["evo"])
        for kind in ("continuous", "binary"):
            frame["y"] = x + np.sin(2*evolution) + rng.normal(scale=.3, size=n)
            if kind == "binary":
                frame["y"] = rng.binomial(1, 1/(1+np.exp(-frame.y.to_numpy())))
            config["trait_type"] = kind
            with threadpool_limits(limits=1):
                result = abm.fit_predict(frame, groups, config)
                poisoned = frame.copy()
                poisoned.loc[poisoned.split=="test", "y"] = np.nan
                second = abm.fit_predict(poisoned, groups, config)
            for name, bank in result["bundle"]["banks"].items():
                status = bank.attention.fit_status
                self.assertGreater(status["gradient_steps"], 0)
                self.assertGreater(status["selected_tensor_changed_count"], 0)
                self.assertGreater(status["gpu_peak_allocated_bytes"], 0)
                for key in ("query.weight", "key.weight", "context.0.qkv.weight"):
                    self.assertIn(key, status["gradient_names"])
                for key, value in bank.attention.state.items():
                    np.testing.assert_allclose(value, second["bundle"]["banks"][name].attention.state[key], atol=1e-6)
            np.testing.assert_allclose(result["predictions"]["GRID_evolution"], second["predictions"]["GRID_evolution"], atol=1e-6)
            frozen = pickle.loads(pickle.dumps(result["bundle"]))
            with threadpool_limits(limits=1):
                projected = abm.predict_new(frozen, frame.loc[frame.split=="test"].drop(columns="y"))
            np.testing.assert_allclose(result["predictions"]["GRID_evolution"], projected["predictions"]["GRID_evolution"], atol=1e-7)
            matches = result["matches"]["GRID_evolution"]
            np.testing.assert_allclose(matches["prediction"], matches["baseline"]+
                matches["weighted_residual_contribution"].sum(1)+matches["clipping_correction"], atol=1e-10)
            if kind == "binary":
                self.assertTrue(result["predictions"]["GRID_evolution"].between(0, 1).all())


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        torch.set_num_threads(1)
        unittest.main()
