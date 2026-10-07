"""Synthetic regression checks; never read or write UKB data.

Run with the GRID Python environment: python f/test_grid_integration.py
Set GRID_TEST_PLINK2 if PLINK 2 is not on PATH or in the default grid environment.
"""
from pathlib import Path
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
import rdata
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("grid_review_runner", ROOT / "f/grid.py")
grid = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = grid
spec.loader.exec_module(grid)
data = grid.load_module("grid.data")
abm = grid.load_module("grid.abm")


class Integration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="grid-regression-")
        cls.addClassCleanup(cls.temp.cleanup)
        cls.work = w = Path(cls.temp.name)
        cls.plink = os.environ.get("GRID_TEST_PLINK2") or shutil.which("plink2") or str(
            Path.home() / "miniforge3/envs/grid/bin/plink2")
        if not Path(cls.plink).is_file():
            raise unittest.SkipTest("PLINK 2 required for native scoring tests")
        rng = np.random.default_rng(20261007)
        cls.ids = ids = [str(100000 + i) for i in range(600)]
        all_ids = ids + ["-10", "999999"]
        cls.g = g = rng.binomial(2, .35, (602, 24)).astype(np.int8)
        g[4, 3] = -1
        (w / "gen").mkdir()
        cls.beta = beta = rng.normal(0, .1, (24, 4))
        cls.a1_alt = a1_alt = np.arange(24) % 3 != 0
        weights = pd.DataFrame({"CHR": np.repeat([1, 2], 12), "BP": np.arange(101, 125),
            "SNP": [f"rs{i}" for i in range(24)], "A1": np.where(a1_alt, "C", "A"),
            "A2": np.where(a1_alt, "A", "C")})
        for j, pop in enumerate(data.POPS):
            weights[f"beta_{pop}"] = beta[:, j]
        weights.to_csv(w / "weights.tsv", sep="\t", index=False)
        weights[["CHR", "BP", "SNP"]].to_csv(w / "snps.tsv", sep="\t", index=False)
        annotation = weights[["CHR", "BP"]].assign(REF="A", ALT="C", ANC="A")
        annotation["AGE_GEN"] = np.tile([500., 2000., 6000.], 8)
        annotation["AGE_LO"] = annotation.AGE_GEN * .9
        annotation["AGE_HI"] = annotation.AGE_GEN * 1.1
        annotation["AGE_QUAL"] = .95
        annotation["AGE_SOURCE"] = "synthetic_fixture"
        annotation["AGE_METHOD"] = "synthetic_known_age"
        annotation["AGE_UNCERTAINTY"] = "synthetic_bounds"
        annotation.loc[0, "AGE_QUAL"] = .1
        annotation.loc[1, "ANC"] = ""
        for pop, af in zip(data.POPS, [.15, .18, .3, .4]):
            annotation[f"AF_{pop}"] = af
        annotation["FREQ_SOURCE"] = "synthetic_reference"
        annotation.to_csv(w / "annotations.tsv", sep="\t", index=False)
        for chrom in (1, 2):
            prefix = w / "gen" / f"bed{chrom}"
            order = np.arange(602) if chrom == 1 else rng.permutation(602)
            variants = list(range((chrom - 1) * 12, chrom * 12))
            prefix.with_suffix(".fam").write_text("".join(f"0 {all_ids[i]} 0 0 0 -9\n" for i in order))
            prefix.with_suffix(".bim").write_text("".join(f"{chrom} rs{j} 0 {101+j} C A\n" for j in variants))
            packed = bytearray(b"\x6c\x1b\x01")
            bits = {-1: 1, 0: 3, 1: 2, 2: 0}
            for j in variants:
                row = g[order, j]
                for start in range(0, len(row), 4):
                    packed.append(sum(bits[int(v)] << (2*k) for k, v in enumerate(row[start:start+4])))
            prefix.with_suffix(".bed").write_bytes(packed)
            subprocess.run([cls.plink, "--bfile", str(prefix), "--make-pgen", "--out", str(w/"gen"/f"chr{chrom}")],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        families = pd.DataFrame({"eid": ids, "family_id": [f"F{i//2}" for i in range(600)]})
        families.to_csv(w / "families.tsv", sep="\t", index=False)
        (w / "remove.txt").write_text("999999\n")
        (w / "keep.txt").write_text("\n".join(all_ids) + "\n")
        cls.train = np.array([data.family_assignment(families.family_id, .5, 20260904)[v] for v in families.family_id])
        pd.DataFrame({"eid": ids, "split": np.where(cls.train, "train", "test")}).to_csv(w / "split.tsv", sep="\t", index=False)
        raw_g = np.where(g < 0, .7, g)
        raw_score = np.where(a1_alt, raw_g, 2-raw_g) @ beta
        age = rng.normal(50, 5, 602)
        pca = pd.DataFrame({"eid": all_ids, "PC1": rng.normal(size=602), "PC2": rng.normal(size=602)})
        pca.to_csv(w / "pca.tsv", sep="\t", index=False)
        ancestry = pd.DataFrame({"eid": all_ids, "genetic_ancestry": ["EUR" if i % 4 < 2 else "AFR" for i in range(602)]})
        ancestry.to_csv(w / "ancestry.tsv", sep="\t", index=False)
        pheno = pd.DataFrame({"eid": all_ids, "age": age, "sex": rng.integers(0,2,602),
            "PC1": pca.PC1, "PC2": pca.PC2,
            "height": 165 + .1*age + raw_score[:,0] + rng.normal(size=602),
            "ldl": 3 + raw_score[:,1] + rng.normal(size=602),
            "t2dm.Yr2e": np.arange(602)%2, "t2dm.Yt2e": 0})
        # The outer roster stays 300/300 despite asymmetric missing LDL outcomes.
        pheno.loc[np.flatnonzero(cls.train)[:80], "ldl"] = np.nan
        rdata.write_rds(w / "pheno.rds", pheno)
        for trait in data.TRAITS:
            folder = w / "scores" / trait
            folder.mkdir(parents=True)
            scores = pd.DataFrame(raw_score, columns=data.CSX).assign(eid=all_ids)
            rdata.write_rds(folder / "1csx.scores.rds", scores)
        cls.args = args = grid.parser().parse_args([
            "--cache-dir", str(w/"cache"), "--out-root", str(w/"output"),
            "--dir-gen", str(w/"gen"), "--pheno-file", str(w/"pheno.rds"),
            "--score-dir", str(w/"scores"), "--pca-file", str(w/"pca.tsv"),
            "--ancestry-file", str(w/"ancestry.tsv"), "--split-group-file", str(w/"families.tsv"),
            "--split-file", str(w/"split.tsv"), "--remove", str(w/"remove.txt"), "--keep", str(w/"keep.txt"),
            "--weights-file", str(w/"weights.tsv"), "--snpinfo", str(w/"snps.tsv"),
            "--annotation-file", str(w/"annotations.tsv"), "--geva-dir", str(w/"absent-geva"),
            "--chrs", "1-2", "--distance-pcs", "2", "--plink2", cls.plink,
            "--hgb-iterations", "4", "--gate-iterations", "4", "--folds", "3",
            "--k-grid", "8", "--alpha-grid", "0,1", "--radius-grid", "2",
            "--audit-bootstrap", "20", "--bootstrap", "20", "--threads", "1", "--quiet"])
        args.traits = list(data.TRAITS)
        args.prsformer_root = str(w / "absent-benchmark")
        cls.prepared = prepared = grid.evolution(args, grid.prepare(args))
        grid.fit(args, prepared)
        grid.report(args, prepared)

    def test_native_scores_and_shared_split(self):
        self.assertEqual(self.prepared["cohort"].split.value_counts().to_dict(), {"train":300,"test":300})
        self.assertEqual(self.prepared["cohort"].groupby("family_id").split.nunique().max(), 1)
        self.assertEqual(len(self.prepared["tables"]["ldl"].dropna(subset=["y"])), 520)
        freq = np.where(self.g[:600][self.train] < 0, np.nan, self.g[:600][self.train]).astype(float)
        centered = self.g[:600] - np.nanmean(freq, axis=0)
        centered[self.g[:600] < 0] = 0
        expected = (centered * np.where(self.a1_alt, 1, -1)) @ self.beta
        table = self.prepared["tables"]["height"].set_index("eid").loc[self.ids]
        np.testing.assert_allclose(table[[f"evo.total_{p}" for p in data.POPS]], expected, atol=2e-5, rtol=2e-5)
        modules = pd.read_csv(self.work/"cache/evolution/height/modules.tsv.gz",sep="\t")
        scores = table.filter(regex=r"^evo\.").rename(columns=lambda c:c[4:])
        data.check_partition(scores, modules)

    def test_age_permutation_preserves_uncertain_and_unknown(self):
        frame = pd.read_csv(self.work/"cache/evolution/height/variants.tsv.gz",sep="\t")
        self.assertEqual(frame.loc[0,"age_bin"],"uncertain")
        self.assertEqual(frame.loc[1,"age_bin"],"unknown")
        fixed = frame.age_bin.isin(["unknown","uncertain"])
        self.assertTrue((frame.loc[fixed,"age_bin"] == frame.loc[fixed,"perm_age_bin"]).all())
        for _, group in frame.groupby("perm_stratum"):
            self.assertEqual(sorted(group.age_bin),sorted(group.perm_age_bin))

    def test_invalid_model_options_checked_before_annotation_build(self):
        args = grid.parser().parse_args(["--k-grid","0"])
        with self.assertRaisesRegex(ValueError,"k_grid"):
            abm._configuration(grid.model_configuration(args,{},"height",require_features=False))

    def test_empty_withdrawal_file_can_be_fingerprinted(self):
        with tempfile.TemporaryDirectory() as temporary:
            empty = Path(temporary)/"remove.txt"
            empty.touch()
            self.assertEqual(data.id_list(empty),set())
            self.assertEqual(data.source_stamp(empty,allow_empty=True)["size"],0)

    def test_all_cohort_control_files_are_fingerprinted(self):
        stamps = {x["path"] for x in self.prepared["metadata"]["sources"]}
        for name in ("families.tsv", "split.tsv", "remove.txt", "keep.txt"):
            self.assertIn(str(self.work/name), stamps)
            path = self.work/name
            original = path.read_bytes()
            stat = path.stat()
            try:
                path.write_bytes(original + b"\n")
                with self.assertRaisesRegex(ValueError, "Prepared input changed"):
                    grid.validate_input_sources(self.prepared)
            finally:
                path.write_bytes(original)
                os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))

    def test_test_labels_do_not_change_prediction_or_selection(self):
        frame = self.prepared["tables"]["height"].copy()
        first = grid.joblib.load(self.work/"cache/height/fit.joblib")
        frame.loc[frame.split == "test", "y"] = np.arange(300)*100.
        second = abm.fit_predict(frame, self.prepared["feature_groups"]["height"], first["bundle"]["config"])
        pd.testing.assert_frame_equal(first["predictions"], second["predictions"])
        self.assertTrue((first["crossfit_audit"].family_overlap == False).all())

    def test_models_reload_without_outcomes_and_explain_contributions(self):
        for trait in data.TRAITS:
            bundle, metadata = grid.read_model(self.work/f"output/{trait}/grid.model.rds")
            source = self.prepared["tables"][trait]
            query = source[(source.split == "test") & source.y.notna()].drop(columns="y")
            result = abm.predict_new(bundle, query)
            saved = grid.read_rds(self.work/f"output/{trait}/grid.test_individuals.rds")
            np.testing.assert_allclose(result["predictions"].GRID_policy, saved.GRID_policy)
            bank = bundle["banks"][bundle["primary_arm"]]
            baseline = bundle["models"]["CSx"].predict(query)
            pred, diag, matches = bank.borrow(query, baseline, k=8, alpha=1., radius_multiplier=100., retain=True)
            self.assertTrue(np.any(np.abs(matches["weighted_residual_contribution"]) > 0))
            np.testing.assert_allclose(pred, baseline + matches["weighted_residual_contribution"].sum(1) + matches["clipping_correction"])
            valid = matches["neighbor_indices"] >= 0
            safe = np.maximum(matches["neighbor_indices"],0)
            self.assertFalse(np.any((query.family_id.to_numpy()[:,None] == matches["donor_family_ids"][safe]) & valid))

    def test_metrics_and_artifact_formats(self):
        import openpyxl
        for trait in data.TRAITS:
            folder = self.work/"output"/trait
            table = grid.read_rds(folder/"grid.test_individuals.rds")
            book = openpyxl.load_workbook(folder/"grid.performance.xlsx",read_only=True,data_only=True)
            book.close()
            for png in folder.glob("*.png"):
                self.assertTrue(png.with_suffix(".xlsx").is_file())
            for path in folder.glob("*.xlsx"):
                book = openpyxl.load_workbook(path,read_only=True,data_only=True)
                for sheet in book:
                    headers = next(sheet.values,())
                    self.assertFalse(set(headers) & {"eid","IID","donor_eid","family_id"})
                book.close()
            delta = grid.paired_loss_comparison(table, ["GRID_policy"], 20, 3)[0]
            expected = np.mean((table.y-table.GRID_policy)**2-(table.y-table.CSx_full_training)**2)
            self.assertAlmostEqual(delta["delta_loss"],expected)

    def test_scoring_tampering_rejected(self):
        item = self.prepared["scoring_files"]["height"][0]
        path = Path(item["source"])
        original = path.read_bytes()
        try:
            path.write_bytes(bytes([original[0]^1])+original[1:])
            with tempfile.TemporaryDirectory() as destination:
                with self.assertRaisesRegex(ValueError,"Frozen scoring material changed"):
                    grid.publish_scoring_files(self.prepared,"height",Path(destination))
        finally:
            path.write_bytes(original)

    def test_shell_delegates_dotted_entries(self):
        for module,name in [("pca","0.pca.sh"),("csx","1.csx.sh"),("disco","2.disco.sh")]:
            run = subprocess.run(["bash",str(ROOT/"grid.sh"),module,"--help"],capture_output=True,text=True,check=True)
            self.assertIn(name,run.stdout)


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        unittest.main(verbosity=2)
