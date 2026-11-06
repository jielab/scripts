"""Shared preparation must not change when a benchmark selects fewer traits."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

import prsformer
from test_block_dosage import write_pgen

ROOT = Path(__file__).resolve().parents[1]


def commands(script, args):
    result = subprocess.run(
        ["bash", str(ROOT / "f" / script), *args, "--dry-run"],
        env=dict(os.environ, GRID_PYTHON=sys.executable, SCRIPT_CONSOLE_ACTIVE="1"),
        text=True, capture_output=True, check=True)
    return [shlex.split(line[4:]) for line in result.stdout.splitlines() if line.startswith("RUN ")]


class SharedCacheTests(unittest.TestCase):
    def test_phenotype_only_refresh_preserves_frozen_cache(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp)
            source = cache / "phenotypes.tsv"
            source.write_text("updated source")
            cohort = pd.DataFrame(dict(eid=["001", "002"], target=["EUR", "EUR"],
                split=["train", "test"], height=[160.12345678901234, np.nan], group=["1", "2"]))
            cohort.attrs["phenotype_definitions"] = {"height": {"source_column": "height"}}
            cohort.to_csv(cache / "data.tsv.gz", sep="\t", index=False, na_rep="NA")
            (cache / "variants.tsv.gz").write_bytes(b"variants")
            np.save(cache / "genotypes.npy", np.zeros((2, 1), dtype=np.float16))
            old = dict(format_version=1, options={"maf": .01},
                inputs=[dict(path=str(source.resolve()), size=1), dict(path="genotype", size=10)],
                phenotype_code_sha256={"helper": "old"}, signature="frozen", shape=[2, 1],
                phenotype_definitions=cohort.attrs["phenotype_definitions"],
                sha256={n: prsformer.sha256(cache / n) for n in ("data.tsv.gz", "variants.tsv.gz")})
            (cache / "prepare.json").write_text(json.dumps(old))
            current = {k: old[k] for k in ("format_version", "options", "inputs", "phenotype_code_sha256")}
            current = dict(current, inputs=[dict(path=str(source.resolve()), size=20), old["inputs"][1]],
                           phenotype_code_sha256={"helper": "new"})
            frozen = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in cache.iterdir()}
            args = SimpleNamespace(pheno_file=str(source))
            with patch.object(prsformer, "make_cohort", return_value=cohort.copy()), \
                 patch.object(prsformer, "identity", return_value=("new", current)):
                self.assertTrue(prsformer.reuse_equivalent_phenotypes(args, [], cache, current))
                changed = dict(current, inputs=[current["inputs"][0], dict(path="genotype", size=11)])
                self.assertFalse(prsformer.reuse_equivalent_phenotypes(args, [], cache, changed))
                self.assertFalse(prsformer.reuse_equivalent_phenotypes(args, [], cache,
                    dict(current, options={"maf": .02})))
            for column, value in (("height", 161.), ("split", "validation"), ("group", "9"), ("eid", "003")):
                changed = cohort.copy()
                changed.loc[0, column] = value
                with patch.object(prsformer, "make_cohort", return_value=changed):
                    with self.assertRaisesRegex(ValueError, "current phenotype/cohort values differ"):
                        prsformer.reuse_equivalent_phenotypes(args, [], cache, current)
            for path, state in frozen.items():
                self.assertEqual((path.read_bytes(), path.stat().st_mtime_ns), state)

    def test_trait_output_paths_reach_training_and_publication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "absent analysis"
            for trait in ("height", "ldl", "t2dm"):
                base = root / "prsformer" / trait
                planned = commands("prsformer.entry.sh", [
                    "--cache-dir", str(root / "cache/abm"), "--out-root", str(root / "abm"),
                    "--prsformer-root", str(root), "--traits", trait, "--no-replace",
                    "--prsformer-run-dir", str(base / "run"),
                    "--prsformer-publication-dir", str(base), "--prsformer-score-dir", str(root / "prsformer")])
                model = next(cmd for cmd in planned if str(ROOT / "f/prsformer.sh") in cmd)
                stages = commands("prsformer.sh", model[2:])
                for stage in stages:
                    expected = "height,ldl,t2dm" if stage[2] == "prepare" else trait
                    self.assertEqual(stage[stage.index("--traits") + 1], expected)
                    if stage[2] == "train":
                        self.assertEqual(stage[stage.index("--out-dir") + 1], str(base / "run"))
                        self.assertEqual(stage[stage.index("--genotypes") + 1], str(root / "prsformer/cache/genotypes.npy"))
                    elif stage[2] == "publish":
                        self.assertEqual(stage[stage.index("--run-dir") + 1], str(base / "run"))
                        self.assertEqual(stage[stage.index("--publication-dir") + 1], str(base))
                        self.assertEqual(stage[stage.index("--score-dir") + 1], str(root / "prsformer"))
            self.assertFalse(root.exists())

    def test_shared_preparation_and_model_tasks_are_separate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "absent analysis"
            planned = commands("prsformer.entry.sh", [
                "--cache-dir", str(root / "cache/abm"), "--out-root", str(root / "abm"),
                "--prsformer-root", str(root), "--traits", "ldl,height", "--no-replace"])
            cohort = next(cmd for cmd in planned if str(ROOT / "f/main.py") in cmd)
            self.assertEqual(cohort[cohort.index("--traits") + 1], "height,ldl,t2dm")
            model = next(cmd for cmd in planned if str(ROOT / "f/prsformer.sh") in cmd)
            stages = commands("prsformer.sh", model[2:])
            for stage in stages:
                expected = "height,ldl,t2dm" if stage[2] == "prepare" else "ldl,height"
                self.assertEqual(stage[stage.index("--traits") + 1], expected)
            standalone = commands("prsformer.sh", ["--traits", "ldl,height"])
            for stage in standalone:
                self.assertEqual(stage[stage.index("--traits") + 1], "ldl,height")
            self.assertFalse(root.exists())

    def test_subset_reuses_and_restores_frozen_full_trait_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ids = [f"{i:03}" for i in range(18)]
            write_pgen(root / "hap", ids, [("rs1", 100, "A", "G")], [np.arange(18) % 3])
            phenotype = pd.DataFrame(dict(eid=ids, height=np.arange(18) + 160.,
                                          ldl=np.arange(18) / 10 + 2., t2dm=np.arange(18) % 2,
                                          genetic_ancestry="EUR"))
            # This participant belongs to the full roster but has neither selected label.
            phenotype.loc[0, ["height", "ldl"]] = np.nan
            phenotype.to_csv(root / "phenotypes.tsv", sep="\t", index=False)
            pd.DataFrame(dict(eid=ids, split=["train"] * 6 + ["validation"] * 6 + ["test"] * 6)).to_csv(
                root / "split.tsv", sep="\t", index=False)
            pd.DataFrame(dict(eid=ids, group=ids)).to_csv(root / "families.tsv", sep="\t", index=False)
            (root / "keep").write_text("\n".join(ids) + "\n")
            cache = root / "cache"
            common = ["--mode", "prepare", "--require-shared-cohort", "--dir-gen", str(root / "hap"),
                      "--chrs", "1", "--pheno-file", str(root / "phenotypes.tsv"),
                      "--ancestry-file", "none", "--covariates", "none", "--t2dm-col", "t2dm",
                      "--snp-list", "none", "--remove", "none", "--keep", str(root / "keep"),
                      "--split-file", str(root / "split.tsv"), "--split-group-file", str(root / "families.tsv"),
                      "--cache-dir", str(cache)]
            full = commands("prsformer.sh", [*common, "--traits", "height,ldl,t2dm"])[0][2:]
            subset = commands("prsformer.sh", [*common, "--traits", "height,ldl"])[0][2:]
            prsformer.main(full)
            files = [cache / name for name in ("prepare.json", "data.tsv.gz", "variants.tsv.gz", "genotypes.npy")]
            frozen = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in files}
            with patch.object(prsformer, "make_cohort", side_effect=AssertionError("unexpected rebuild")):
                prsformer.main(subset)
            # Simulate the verified release marker, then rebuild from the tiny real PGEN.
            marker = cache / "genotypes.released.json"
            marker.write_text(json.dumps(dict(schema="grid.dosage_release.1",
                prepare_sha256=prsformer.sha256(cache / "prepare.json"),
                genotype_sha256=prsformer.sha256(cache / "genotypes.npy"),
                genotype_mtime_ns=(cache / "genotypes.npy").stat().st_mtime_ns)))
            (cache / "genotypes.npy").unlink()
            prsformer.main(subset)
            self.assertFalse(marker.exists())
            for path, (contents, stamp) in frozen.items():
                self.assertEqual(path.read_bytes(), contents)
                self.assertEqual(path.stat().st_mtime_ns, stamp)
            model = prsformer.parser().parse_args([
                "train", "--genotypes", str(cache / "genotypes.npy"), "--data", str(cache / "data.tsv.gz"),
                "--variants", str(cache / "variants.tsv.gz"), "--out-dir", str(root / "run"),
                "--upstream-dir", str(root / "upstream"), "--traits", "height,ldl"])
            cohort, labels, _, _ = prsformer.read_inputs(model, ["height", "ldl"])
            self.assertEqual(cohort.eid.tolist(), ids)
            self.assertEqual(labels.shape, (18, 2))
            self.assertTrue(np.isnan(labels[0]).all())
            with self.assertRaisesRegex(ValueError, "Cache input/settings mismatch"):
                prsformer.main([*subset, "--maf", "0.02"])


if __name__ == "__main__":
    unittest.main()
