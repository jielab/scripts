"""Boundaries that prevent label leakage and accidental checkpoint misuse."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"f"))
from panome_tf import parse_args, batch_configs
from tf_download import verify_model
from tf_model import context_indices, fold_labels, select_features


class TransformerWorkflowTests(unittest.TestCase):
    def test_batch_has_distinct_input_defaults_and_tf_output_root(self):
        runs = list(batch_configs(parse_args(["--Y", "cvd_cad,ra", "--biom", "prot,met"])))
        self.assertEqual([(a.trait, a.biom) for a in runs],
                         [("cvd_cad", "prot"), ("cvd_cad", "met"), ("ra", "prot"), ("ra", "met")])
        for a in runs:
            self.assertEqual(a.analysis_root, "/mnt/d/analysis/panome_TF")
            self.assertEqual(a.transform, "log1p" if a.biom == "met" else "none")
            self.assertEqual(a.residualize, a.biom+".plate")
            self.assertEqual(a.learning_rate, 1e-5)

    def test_ambiguous_batch_paths_are_rejected(self):
        for args in [["--biom", "prot,met", "--omics-file", "x.csv"],
                     ["--Y", "ra,ra"], ["--Y", "ra,cvd_cad", "--run-dir", "/tmp/x"]]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                list(batch_configs(parse_args(args)))

    def test_git_lfs_pointer_is_not_accepted_as_weights(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"model.ckpt"
            path.write_text("version https://git-lfs.github.com/spec/v1\n")
            with self.assertRaisesRegex(FileNotFoundError, "pointer"):
                verify_model(path)

    def test_context_sampling_never_reintroduces_excluded_queries_or_families(self):
        y = np.tile([0, 0, 1, 1], 50)
        groups = np.repeat(np.arange(100).astype(str), 2)
        folds = fold_labels(y, groups, 5, 42)
        for fold in range(5):
            pool = np.flatnonzero(folds != fold)
            query = np.flatnonzero(folds == fold)
            ctx = context_indices(pool, y, np.ones(len(y)), 30, np.random.default_rng(42))
            self.assertFalse(set(ctx) & set(query))
            self.assertFalse(set(groups[ctx]) & set(groups[query]))
            self.assertEqual(len(ctx), len(set(ctx)))
            self.assertEqual(set(y[ctx]), {0, 1})

    def test_zero_weight_censored_labels_cannot_change_feature_selection(self):
        rng = np.random.default_rng(42)
        x = rng.normal(size=(200, 8))
        y = rng.integers(0, 2, len(x))
        w = np.ones(len(x))
        w[::3] = 0
        ix, score = select_features(x, y, w, 4)
        y[w == 0] = 1-y[w == 0]
        changed, changed_score = select_features(x, y, w, 4)
        np.testing.assert_array_equal(ix, changed)
        np.testing.assert_array_equal(score, changed_score)


if __name__ == "__main__":
    unittest.main()
