import unittest
import numpy as np
from abm_core.representation import BlockGeometry


class RepresentationTests(unittest.TestCase):
    def test_map_balanced_rms_distance(self):
        # One map has one variable block, the other has two; maps weigh equally.
        x = np.array([[0., 0., 0., 9.], [1., 1., 1., 9.], [2., 2., 2., 9.]])
        g = BlockGeometry.fit(x, ["AFR", "EUR", "EUR", "EUR"], floor_quantile=0)
        z = g.transform(x)
        self.assertEqual(g.keep, [0, 1, 2])
        self.assertAlmostEqual(np.sum((z[2] - z[0]) ** 2), 6.0)

    def test_query_does_not_refit_scale(self):
        x = np.array([[0., 1.], [1., 2.], [2., 4.]])
        g = BlockGeometry.fit(x, ["A", "B"])
        before = g.fingerprint()
        one = g.transform(x[:1])
        many = g.transform(np.vstack([x[:1], np.array([[1000., -9000.]])]))
        np.testing.assert_array_equal(one, many[:1])
        self.assertEqual(g.fingerprint(), before)

    def test_missing_block_values_are_not_hidden_by_imputation(self):
        x = np.array([[0., 1.], [1., 2.]])
        g = BlockGeometry.fit(x, ["A", "B"])
        with self.assertRaisesRegex(ValueError, "Invalid genotype"):
            g.transform([[np.nan, 3.]])


if __name__ == "__main__":
    unittest.main()
