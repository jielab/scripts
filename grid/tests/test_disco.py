"""Regression checks for the independent published-method R wrapper."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import numpy as np
import pandas as pd


def r_runtime():
    env = os.environ.copy()
    for key in ("R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_ENVIRON_USER"):
        env.pop(key, None)
    for r in (os.getenv("GRID_RSCRIPT"), str(Path.home()/"miniforge3/envs/grid/bin/Rscript"), shutil.which("Rscript")):
        if r and Path(r).is_file():
            cmd = [r, "--vanilla"]
            probe = subprocess.run([*cmd, "-e", 'stopifnot(all(vapply(c("data.table","dplyr","stringr","rio","optparse"),requireNamespace,logical(1),quietly=TRUE)))'],
                                   env=env, capture_output=True, timeout=60)
            if probe.returncode == 0:
                return cmd, env
    raise RuntimeError("No R runtime with the upstream wrapper dependencies")


class PublishedInterpolationTests(unittest.TestCase):
    def test_real_r_interpolation_preserves_ids_with_and_without_pca_regression(self):
        try:
            command, env = r_runtime()
        except RuntimeError as exc:
            self.skipTest(str(exc))
        root = Path(__file__).resolve().parents[1]
        rng = np.random.default_rng(17)
        n = 24
        ids = [f"{i:04d}" for i in range(n)]
        pops = ["AFR", "EAS", "EUR", "SAS"]
        centers = np.eye(4, 10)
        pcs = rng.normal(size=(n, 10))
        scores = rng.normal(size=(n, 4))
        # Equidistant centers give equal shrinkage: normalized inverse distances.
        weights = 1 / np.linalg.norm(pcs[:, None, :] - centers[None, :, :], axis=2)
        expected = (weights * scores).sum(axis=1) / weights.sum(axis=1)
        with tempfile.TemporaryDirectory() as td:
            folder = Path(td)
            c = pd.DataFrame(centers, columns=[f"PC{i}" for i in range(1, 11)])
            c.insert(0, "cohort", pops)
            c.to_csv(folder / "centers.tsv", sep="\t", index=False)
            p = pd.DataFrame(pcs, columns=[f"PC{i}" for i in range(1, 11)])
            p.insert(0, "id", ids)
            p.iloc[rng.permutation(n)].to_csv(folder / "pca.tsv", sep="\t", index=False)
            files = []
            for i, pop in enumerate(pops):
                path = folder / (pop + ".tsv")
                pd.DataFrame(dict(IID=ids, PRS=scores[:, i])).iloc[rng.permutation(n)].to_csv(
                    path, sep="\t", index=False)
                files.append(str(path))
            for regress in ("FALSE", "TRUE"):
                with self.subTest(regress=regress):
                    out = folder / regress
                    result = subprocess.run([*command, str(root / "f/3.disco.R"),
                        "--med.file", str(folder / "centers.tsv"), "--pca.file", str(folder / "pca.tsv"),
                        "--prs.list", ",".join(files), "--select.col", "IID,PRS",
                        "--regress.PCA", regress, "--print.coef", "FALSE", "--out", str(out)],
                        env=env, capture_output=True, text=True, timeout=60)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    actual = pd.read_csv(str(out) + ".tsv.gz", sep="\t", dtype={"IID": str})
                    self.assertEqual(set(actual.IID), set(ids))
                    self.assertTrue(actual.IID.is_unique)
                    self.assertTrue(np.isfinite(actual.PRS).all())
                    if regress == "FALSE":
                        np.testing.assert_allclose(actual.set_index("IID").loc[ids, "PRS"],
                                                   expected, rtol=1e-10, atol=1e-10)

