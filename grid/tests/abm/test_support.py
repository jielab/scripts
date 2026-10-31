"""Meaningful support-definition and exclusion checks; no UKB data required."""
import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from abm_core.support import count_twins, estimate_radius


class SupportTests(unittest.TestCase):
    def call(self, donors, queries, **kwargs):
        n, m = len(donors), len(queries)
        defaults = dict(donor_ids=[f"d{i:03d}" for i in range(n)],
                        query_ids=[f"q{i:03d}" for i in range(m)],
                        donor_groups=[f"g{i:03d}" for i in range(n)],
                        query_groups=[f"h{i:03d}" for i in range(m)],
                        radius=2.0, min_count=1, min_group_ess=1, device="cpu",
                        query_batch=2, donor_batch=3)
        defaults.update(kwargs)
        return count_twins(np.asarray(donors, dtype=float), np.asarray(queries, dtype=float), **defaults)

    def test_count_is_not_k_and_density_uses_eligible_bank(self):
        d = np.linspace(0.0, 0.8, 10).reshape(-1, 1)
        small = self.call(d, [[0.0]], radius=1, k_used=2)
        large = self.call(d, [[0.0]], radius=1, k_used=8)
        self.assertEqual(small.table.n_twins.item(), 10)
        self.assertEqual(large.table.n_twins.item(), 10)
        self.assertEqual(small.table.n_used.item(), 2)
        self.assertEqual(large.table.n_used.item(), 8)
        self.assertEqual(small.table.twins_per_10000.item(), 10000)

    def test_self_and_family_exclusion_even_if_id_group_disagree(self):
        r = self.call([[0], [.1], [.2], [.3]], [[0]], donor_ids=["a", "b", "c", "d"],
                      query_ids=["a"], donor_groups=["other", "family", "free", "free"],
                      query_groups=["family"], k_used=4)
        self.assertEqual(r.table.n_eligible.item(), 2)
        self.assertEqual(r.table.n_twins.item(), 2)
        self.assertEqual(r.table.n_groups_used.item(), 1)
        self.assertAlmostEqual(r.table.ess_group.item(), 1)
        self.assertGreater(r.table.ess_individual.item(), 1.9)
        np.testing.assert_array_equal(r.neighbor_indices[0], [2, 3, -1, -1])

    def test_support_requires_group_ess_not_only_headcount(self):
        r = self.call([[0], [.1], [.2]], [[0]], donor_groups=["g", "g", "g"],
                      min_count=3, min_group_ess=2)
        self.assertFalse(r.table.supported.item())
        self.assertEqual(r.table.support_reason.item(), "low_group_ess")

    def test_empty_radius_still_reports_nearest_eligible_distance(self):
        r = self.call([[10], [12]], [[0]], radius=.1, k_used=4)
        self.assertEqual(r.table.n_twins.item(), 0)
        self.assertFalse(r.table.supported.item())
        self.assertEqual(r.table.d_nearest.item(), 10)
        self.assertTrue(np.isinf(r.table.d_k_used.item()))
        self.assertTrue((r.neighbor_indices == -1).all())
        self.assertTrue(np.isinf(r.neighbor_distances).all())

    def test_no_eligible_and_empty_bank(self):
        r = self.call([[0]], [[0]], donor_groups=["g"], query_groups=["g"])
        self.assertEqual(r.table.support_reason.item(), "no_eligible_donors")
        self.assertTrue(np.isnan(r.table.twins_per_10000.item()))
        empty = self.call(np.empty((0, 1)), [[0]])
        self.assertEqual(empty.table.n_eligible.item(), 0)

    def test_ties_are_canonical_and_permutation_independent(self):
        d = np.asarray([[1, 0], [-1, 0], [0, 1], [0, -1]])
        ids = np.asarray(["z", "b", "a", "m"])
        a = self.call(d, [[0, 0]], donor_ids=ids, k_used=2)
        p = np.asarray([3, 1, 0, 2])
        b = self.call(d[p], [[0, 0]], donor_ids=ids[p], donor_groups=np.asarray(["g0", "g1", "g2", "g3"])[p],
                      k_used=2, donor_batch=1)
        self.assertEqual(ids[a.neighbor_indices[0]].tolist(), ["a", "b"])
        self.assertEqual(ids[p][b.neighbor_indices[0]].tolist(), ["a", "b"])

    def test_exact_radius_boundary(self):
        r = self.call([[1], [1 + 1e-8], [1 - 1e-8]], [[0]], radius=1)
        self.assertEqual(r.table.n_twins.item(), 2)

    def test_blocked_search_matches_direct_brute_force(self):
        rng = np.random.default_rng(12)
        d, q = rng.normal(size=(71, 5)), rng.normal(size=(7, 5))
        ids = np.asarray([f"id{i:03d}" for i in rng.permutation(len(d))])
        dg = np.asarray([f"family{i // 3}" for i in range(len(d))])
        qg = np.asarray([f"family{i}" for i in range(len(q))])
        out = self.call(d, q, donor_ids=ids, donor_groups=dg, query_groups=qg,
                        radius=2.7, k_used=9, donor_batch=13, query_batch=3)
        exact = np.sqrt(np.square(q[:, None, :] - d[None, :, :]).sum(axis=2))
        for row in range(len(q)):
            valid = np.flatnonzero((exact[row] <= 2.7) & (dg != qg[row]))
            expected = sorted(valid, key=lambda j: (exact[row, j], ids[j]))[:9]
            used = out.neighbor_indices[row]
            self.assertEqual(used[used >= 0].tolist(), expected)
            self.assertEqual(out.table.n_twins.iloc[row], len(valid))
            self.assertEqual(out.table.n_eligible.iloc[row], int((dg != qg[row]).sum()))
            np.testing.assert_allclose(out.neighbor_distances[row, :len(expected)],
                                       exact[row, expected], rtol=1e-12, atol=1e-12)

    def test_radius_fits_one_threshold_not_query_specific_k(self):
        d = np.asarray([[0.0], [1.0], [2.0], [10.0]])
        q = np.asarray([[0.2], [1.2]])
        radius = estimate_radius(d, q, donor_ids=["d0", "d1", "d2", "d3"],
                                 query_ids=["q0", "q1"], donor_groups=["d0", "d1", "d2", "d3"],
                                 query_groups=["q0", "q1"], k_ref=2, quantile=.9, device="cpu")
        self.assertAlmostEqual(radius, .8)
        out = self.call(d, [[0.2], [100]], radius=radius, k_used=3)
        self.assertEqual(out.table.n_twins.tolist(), [2, 0])

    def test_invalid_geometry_and_identifiers_fail_explicitly(self):
        with self.assertRaises(ValueError):
            self.call([[0], [1]], [[0]], donor_ids=["a", "a"])
        with self.assertRaises(ValueError):
            self.call([[0]], [[np.nan]])
        with self.assertRaises(ValueError):
            self.call([[0]], [[0]], donor_groups=[None])
        with self.assertRaises(ValueError):
            self.call([[0]], [[0]], radius=0)

    def test_provisional_group_status_is_explicit(self):
        r = self.call([[0], [.1]], [[0]], groups_verified=False)
        self.assertTrue(r.table.supported.item())
        self.assertEqual(r.table.support_reason.item(), "supported_provisional_groups")
        self.assertFalse(r.table.groups_verified.item())

    def test_explicit_cuda_never_silently_falls_back(self):
        if importlib.util.find_spec("torch") is None:
            with self.assertRaises(RuntimeError):
                self.call([[0]], [[0]], device="cuda")
        else:
            import torch
            if not torch.cuda.is_available():
                with self.assertRaises(RuntimeError):
                    self.call([[0]], [[0]], device="cuda")

    def test_cuda_parity_if_available(self):
        if importlib.util.find_spec("torch") is None:
            self.skipTest("PyTorch not installed; CUDA execution was not verified")
        import torch
        if not torch.cuda.is_available():
            self.skipTest("No CUDA device; CUDA execution was not verified")
        rng = np.random.default_rng(19)
        d, q = rng.normal(size=(45, 7)), rng.normal(size=(9, 7))
        kwargs = dict(radius=3, k_used=8, donor_batch=11, query_batch=3)
        cpu = self.call(d, q, device="cpu", **kwargs)
        gpu = self.call(d, q, device="cuda", **kwargs)
        np.testing.assert_array_equal(cpu.neighbor_indices, gpu.neighbor_indices)
        np.testing.assert_allclose(cpu.neighbor_distances, gpu.neighbor_distances, atol=1e-10, rtol=1e-10)
        for col in ["n_twins", "n_eligible", "n_used", "n_groups_used", "supported"]:
            np.testing.assert_array_equal(cpu.table[col], gpu.table[col])
        np.testing.assert_allclose(cpu.table.ess_group, gpu.table.ess_group, rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
