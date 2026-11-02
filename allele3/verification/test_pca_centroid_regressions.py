"""Small regressions for explicit cohort identity and separate modern markers."""
import importlib.util
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"f"))
spec = importlib.util.spec_from_file_location("pca3d_entry", ROOT/"f/0.pca3D.py")
pca = importlib.util.module_from_spec(spec); spec.loader.exec_module(pca)


class CentroidRegressions(unittest.TestCase):
    def fixture(self, root):
        coordinates = root/"scores.sscore"
        coordinates.write_text("#IID PC1_AVG PC2_AVG\na0 0 0\na1 1 1\nm1 2 2\n")
        metadata = root/"meta.tsv"
        metadata.write_text("sample_id\tcohort\tage_bp\tinclude\na0\tancient\t0\t1\na1\tancient\t3000\t1\nm1\tmodern\t0\t1\n")
        dosages = root/"dosages.tsv"
        dosages.write_text("sample_id\tvariant_id\tdosage\tstatus\na0\ts\t1\tok\na1\ts\t1\tok\nm1\ts\t1\tok\n")
        return SimpleNamespace(coordinates=coordinates, metadata=metadata, dosages=dosages, out=root/"out", bin_years=2000)

    def test_age_zero_does_not_infer_modern(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = self.fixture(Path(tmp))
            a.metadata.write_text("sample_id\tage_bp\tinclude\na0\t0\t1\n")
            with self.assertRaisesRegex(ValueError, "explicit ancient/modern cohort"):
                pca.run_centroids(a)

    def test_modern_not_connected_to_ancient_line(self):
        import matplotlib
        matplotlib.use("Agg")
        from mpl_toolkits.mplot3d import Axes3D
        with tempfile.TemporaryDirectory() as tmp:
            a = self.fixture(Path(tmp)); seen = []; original = Axes3D.plot
            def record_plot(axis, xs, ys, *rest, **kwargs):
                seen.append((list(xs), list(ys), list(rest[0]), kwargs.get("label")))
                return original(axis, xs, ys, *rest, **kwargs)
            with patch.object(Axes3D, "plot", record_plot):
                pca.run_centroids(a)
            self.assertEqual(len(seen), 1)
            self.assertEqual(seen[0][2], [1000., 3000.])
            self.assertEqual(seen[0][3], "s ancient")
            rows = (a.out/"allele_centroids.tsv").read_text()
            self.assertIn("s\tancient\t0\t2000", rows)
            self.assertIn("s\tmodern\t0\t0", rows)


if __name__ == "__main__":
    unittest.main()
