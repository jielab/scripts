"""Regression tests plus a full CLI pipeline using synthetic data and fake Atlas."""

from __future__ import annotations

import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import pysam
import zarr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
from atlas_query import query, response_matrix, variant_id  # noqa: E402
from fake_atlas import FakeAtlas, vector  # noqa: E402
from alphagenome.data import genome  # noqa: E402
import validate  # noqa: E402

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=6000>
##contig=<ID=chr2,length=6000>
##contig=<ID=chr3,length=6000>
##FILTER=<ID=LowQual,Description="low quality">
##INFO=<ID=AC,Number=A,Type=Integer,Description="count">
##INFO=<ID=AF,Number=A,Type=Float,Description="frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="number">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def command(*args, env=None):
    return subprocess.run(
        [str(x) for x in args],
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )


def extract(directory, body, sequence="A" * 6000):
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / "ref.fa"
    fasta.write_text(">chr1\n" + sequence + "\n")
    pysam.faidx(str(fasta))
    (directory / "input.vcf").write_text(HEADER + body)
    output = directory / "candidates.parquet"
    result = command(
        sys.executable,
        ROOT / "f/extract_sites.py",
        "--vcf-dir",
        directory,
        "--fasta",
        fasta,
        "--out",
        output,
        "--sites-out",
        directory / "sites.parquet",
    )
    return result, output


class InputTests(unittest.TestCase):
    def test_indexed_cohort_region_retains_counts_and_split_alternatives(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            fasta = base / "ref.fa"
            fasta.write_text(">chr1\n" + "A" * 6000 + "\n")
            pysam.faidx(str(fasta))
            source = base / "cohort.vcf"
            header = (
                HEADER.replace(
                    "#CHROM",
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM',
                ).rstrip()
                + "\tFORMAT\tS1\tS2\n"
            )
            source.write_text(
                header
                + "chr1\t9\t.\tA\tC\t.\tPASS\tAC=1;AN=4\tGT\t0/1\t0/0\n"
                + "chr1\t10\t.\tA\tC\t.\tPASS\tAC=1;AN=4\tGT\t0/1\t0/0\n"
                + "chr1\t10\t.\tA\tG\t.\tPASS\tAC=1;AN=4\tGT\t0/0\t0/1\n"
                + "chr1\t20\t.\tA\tT\t.\tPASS\tAC=1;AN=4\tGT\t0/1\t0/0\n"
                + "chr1\t21\t.\tA\tC\t.\tPASS\tAC=1;AN=4\tGT\t0/1\t0/0\n"
            )
            compressed = str(source) + ".gz"
            pysam.tabix_compress(str(source), compressed)
            pysam.tabix_index(compressed, preset="vcf")
            out = base / "output"
            result = command(
                ROOT / "allele.sh",
                "prepare",
                "--vcf",
                compressed,
                "--region",
                "1:10-20",
                "--fasta",
                fasta,
                "--out",
                out,
                env={
                    **os.environ,
                    "ALLELE_CONFIG": "/dev/null",
                    "ALLELE_ENV_NAME": "allele",
                },
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            path = out / "work/candidates.parquet"
            data = pd.read_parquet(path)
            self.assertEqual(set(data.pos), {10, 20})
            self.assertTrue((data.loc[data.pos == 10, "occupancy"] == 2).all())
            observed = data[data.observed == 1]
            self.assertTrue((observed.ac == 1).all())
            self.assertTrue((observed.an == 4).all())
            self.assertTrue((observed.af == 0.25).all())
            manifest = json.loads(Path(str(path) + ".inputs.json").read_text())
            self.assertEqual(manifest["files"][0]["sample_count"], 2)
            Path(compressed + ".tbi").unlink()
            result = command(
                sys.executable,
                ROOT / "f/extract_sites.py",
                "--vcf",
                compressed,
                "--region",
                "chr1:10-20",
                "--fasta",
                fasta,
                "--out",
                base / "unused.parquet",
                "--sites-out",
                base / "unused-sites.parquet",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("requires a VCF .tbi or .csi index", result.stdout)

    def test_controls_remain_ref_matched_when_context_is_sparse(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            data = pd.DataFrame(
                [
                    dict(
                        site_id="multi",
                        chrom="chr1",
                        pos=1,
                        ref="A",
                        context="AAA",
                        occupancy=2,
                        primary_alt="C",
                    ),
                    dict(
                        site_id="a_control",
                        chrom="chr1",
                        pos=2,
                        ref="A",
                        context="TAT",
                        occupancy=1,
                        primary_alt="G",
                    ),
                    dict(
                        site_id="c_control",
                        chrom="chr1",
                        pos=3,
                        ref="C",
                        context="CCC",
                        occupancy=1,
                        primary_alt="G",
                    ),
                ]
            )
            data.to_parquet(base / "input.parquet", index=False)
            result = command(
                sys.executable,
                ROOT / "f/select_atlas_sites.py",
                "--candidates",
                base / "input.parquet",
                "--out",
                base / "selected.parquet",
                "--ratio",
                "2",
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            selected = pd.read_parquet(base / "selected.parquet")
            self.assertEqual(set(selected.site_id), {"multi", "a_control"})

    def test_fixed_alt_is_not_mislabeled_as_biallelic(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, _ = extract(
                Path(tmp), "chr1\t10\t.\tA\tC\t.\tPASS\tAC=100;AN=100\n"
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("No eligible", result.stdout)

    def test_mixed_indel_snv_counts_and_zero_count_alt(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, out = extract(
                Path(tmp),
                "chr1\t10\t.\tA\tAT,C,G\t.\tPASS\tAF=0.4,0.1,0;AC=40,10,0;AN=100\n"
                "chr1\t10\t.\tA\tT\t.\tPASS\tAC=5;AN=100\n"
                "chr1\t20\t.\tA\tC\t.\tLowQual\tAC=10;AN=100\n",
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            data = pd.read_parquet(out).set_index("alt")
            self.assertEqual(len(data), 3)
            self.assertAlmostEqual(data.loc["C", "af"], 0.1, places=6)
            self.assertAlmostEqual(data.loc["T", "af"], 0.05, places=6)
            self.assertEqual(data.loc["G", "observed"], 0)
            self.assertTrue((data.occupancy == 2).all())
            self.assertEqual(data.loc["C", "is_primary"], 1)

    def test_wrong_reference_fails_clearly(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, _ = extract(Path(tmp), "chr1\t10\t.\tC\tG\t.\tPASS\tAC=1;AN=100\n")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("disagree with FASTA", result.stdout)

    def test_empty_filtered_input_fails_clearly(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, _ = extract(Path(tmp), "chr1\t10\t.\tA\tC\t.\tPASS\tAC=0;AN=100\n")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("No eligible", result.stdout)

    def test_shell_validation_and_missing_environment(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = {
                **os.environ,
                "ALLELE_CONFIG": "/dev/null",
                "ALLELE_ENV_NAME": "allele_missing_test_environment",
                "ALLELE_OUT": tmp + "/unused",
            }
            result = command(ROOT / "allele.sh", "--help", env=env)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertFalse(Path(env["ALLELE_OUT"]).exists())
            result = command(ROOT / "allele.sh", "prepare", "--fasta", env=env)
            self.assertIn("missing value for --fasta", result.stdout)
            result = command(ROOT / "install.sh", "--check", env=env)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Run:", result.stdout)
            self.assertNotIn("EnvironmentLocationNotFound", result.stdout)

    def test_check_returns_failure_for_broken_import(self):
        real_import = validate.importlib.import_module

        def broken(name):
            if name == "numpy":
                raise ImportError("simulated broken binary")
            return real_import(name)

        with patch.object(
            sys, "argv", ["validate", "--root", str(ROOT), "--out", "/tmp"]
        ), patch.object(
            validate.importlib, "import_module", side_effect=broken
        ), contextlib.redirect_stdout(
            io.StringIO()
        ), contextlib.redirect_stderr(
            io.StringIO()
        ):
            self.assertEqual(validate.main(), 1)


class AtlasTests(unittest.TestCase):
    def sites(self):
        return pd.DataFrame(
            [
                dict(
                    site_id=f"1:{pos}:A",
                    chrom="1",
                    pos=pos,
                    ref="A",
                    occupancy=2,
                    atlas_group="multiallelic",
                )
                for pos in [10, 20]
            ]
        )

    def test_alignment_checks_chromosome_and_duplicates(self):
        variants = [genome.Variant("chr1", 10, "A", alt) for alt in "CGT"]
        response = FakeAtlas().query_variants(variants, requested_scorers=["ATAC"])
        expected = [variant_id(v) for v in variants]
        matrix, _ = response_matrix(response, ["ATAC"], expected)
        np.testing.assert_allclose(
            matrix, np.stack([vector(v, "ATAC") for v in variants])
        )
        with self.assertRaisesRegex(ValueError, "response variants differ"):
            response_matrix(
                response, ["ATAC"], [v.replace("chr1:", "chr2:") for v in expected]
            )
        response["ATAC"].obs["variant"] = [variants[0]] * 3
        with self.assertRaisesRegex(ValueError, "duplicate variant"):
            response_matrix(response, ["ATAC"], expected)

    def test_checkpoint_resume_and_changed_request_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            client = FakeAtlas()
            client.fail_at = 2
            with patch("alphagenome.atlas.atlas.create", return_value=client):
                with self.assertRaisesRegex(ValueError, "simulated interrupted"):
                    query(self.sites(), out, ["ATAC"], 3, 1, "test-only")
            self.assertEqual(len(pd.read_parquet(out / "rows.parquet")), 3)
            first = zarr.open(str(out / "vectors.zarr"), mode="r")[:].copy()
            client = FakeAtlas()
            with patch("alphagenome.atlas.atlas.create", return_value=client):
                query(self.sites(), out, ["ATAC"], 3, 1, "test-only")
                self.assertEqual(client.calls, 1)
                query(self.sites(), out, ["ATAC"], 3, 1, "test-only")
                self.assertEqual(client.calls, 1)
            np.testing.assert_array_equal(
                zarr.open(str(out / "vectors.zarr"), mode="r")[:3], first
            )
            self.assertTrue(
                json.loads((out / "atlas_manifest.json").read_text())["complete"]
            )
            with self.assertRaisesRegex(ValueError, "different sites/scorers/SDK"):
                query(self.sites(), out, ["DNASE"], 3, 1, "test-only")

    def test_geometry_known_vectors_and_zero_vectors(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            rows = pd.DataFrame(
                [
                    dict(
                        site_id="chr1:10:A",
                        chrom="chr1",
                        pos=10,
                        ref="A",
                        alt=alt,
                        variant_id=f"chr1:10:A>{alt}",
                        occupancy=2,
                        observed=int(alt in "CG"),
                        is_primary=int(alt == "C"),
                        primary_alt="C",
                        af={"C": 0.1, "G": 0.05, "T": 0}[alt],
                    )
                    for alt in "CGT"
                ]
            )
            rows.to_parquet(base / "candidates.parquet", index=False)
            rows[["variant_id"]].to_parquet(base / "rows.parquet", index=False)
            pd.DataFrame(
                [{"site_id": "chr1:10:A", "atlas_group": "multiallelic"}]
            ).to_parquet(base / "sites.parquet", index=False)
            pd.DataFrame({"feature": ["ATAC|x", "ATAC|y"]}).to_csv(
                base / "feature_names.tsv", sep="\t", index=False
            )
            (base / "atlas_manifest.json").write_text('{"complete":true}')
            matrix = zarr.open(
                str(base / "vectors.zarr"), mode="w", shape=(3, 2), dtype="f4"
            )
            matrix[:] = [[1, 0], [0, 1], [1, 1]]

            def calculate():
                result = command(
                    sys.executable,
                    ROOT / "f/geometry.py",
                    "--candidates",
                    base / "candidates.parquet",
                    "--atlas-sites",
                    base / "sites.parquet",
                    "--atlas-dir",
                    base,
                    "--out-dir",
                    base / "results",
                )
                self.assertEqual(result.returncode, 0, result.stdout)
                return pd.read_parquet(
                    base / "results/geometry_candidates.parquet"
                ).set_index("alt")

            g = calculate()
            self.assertAlmostEqual(
                g.loc["G", "affinity_to_primary"], (1 - np.sqrt(2)) / (1 + np.sqrt(2))
            )
            self.assertAlmostEqual(g.loc["G", "cosine_to_primary"], 0)
            self.assertAlmostEqual(g.loc["G", "novelty_perpendicular"], 1)
            self.assertAlmostEqual(
                g.loc["T", "affinity_to_primary"], (np.sqrt(2) - 1) / (np.sqrt(2) + 1)
            )
            self.assertAlmostEqual(g.loc["T", "cosine_to_primary"], 1 / np.sqrt(2))
            matrix[:] = 0
            g = calculate()
            for name in [
                "effect_magnitude",
                "cosine_to_primary",
                "affinity_to_primary",
                "novelty_perpendicular",
            ]:
                np.testing.assert_allclose(g[name], 0)


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = None
        keep = os.getenv("ALLELE_TEST_OUT")
        if keep:
            cls.base = Path(keep)
            cls.base.mkdir(parents=True, exist_ok=False)
        else:
            cls.temporary = tempfile.TemporaryDirectory(prefix="allele-test-")
            cls.base = Path(cls.temporary.name)
        inputs = cls.base / "input"
        inputs.mkdir()
        cls.out = cls.base / "output"
        rng = np.random.default_rng(20261015)
        fasta = inputs / "ref.fa"
        records, sequences = [], []
        for chrom in ["chr1", "chr2", "chr3"]:
            sequence = "".join(rng.choice(list("ACGT"), 6000))
            sequences.append(f">{chrom}\n{sequence}\n")
            for i in range(180):
                pos = 100 + 10 * i
                ref = sequence[pos - 1]
                alts = rng.permutation([v for v in "ACGT" if v != ref])
                occupancy = [1, 1, 1, 2, 2, 3][i % 6]
                chosen = alts[:occupancy]
                counts = rng.integers(1, 25, occupancy)
                records.append(
                    f'{chrom}\t{pos}\t.\t{ref}\t{",".join(chosen)}\t.\tPASS\tAC={",".join(map(str, counts))};AN=100\n'
                )
        fasta.write_text("".join(sequences))
        pysam.faidx(str(fasta))
        (inputs / "input.vcf").write_text(HEADER + "".join(records))
        inject = cls.base / "test_imports"
        inject.mkdir()
        (inject / "sitecustomize.py").write_text(
            f"import sys\nif sys.prefix == {sys.prefix!r}:\n    from alphagenome.atlas import atlas\n    from fake_atlas import FakeAtlas\n    atlas.create = lambda key: FakeAtlas()\n"
        )
        cls.env = {
            **os.environ,
            "ALLELE_CONFIG": "/dev/null",
            "ALPHA_GENOME_API_KEY": "offline-test-only",
            "PYTHONPATH": str(inject) + os.pathsep + str(ROOT / "tests"),
            "ALLELE_GNOMAD": str(inputs),
            "ALLELE_VCF_DIR": str(inputs),
            "ALLELE_FASTA": str(fasta),
            "ALLELE_OUT": str(cls.out),
        }
        result = command(ROOT / "allele.sh", "all", "--chunk", "200", env=cls.env)
        (cls.base / "pipeline.log").write_text(result.stdout)
        if result.returncode:
            raise AssertionError(result.stdout)

    @classmethod
    def tearDownClass(cls):
        if cls.temporary:
            cls.temporary.cleanup()

    def test_all_stages_complete_and_exact_occupancy(self):
        counts = pd.read_csv(self.out / "final/Table1_occupancy_counts.tsv", sep="\t")
        self.assertEqual(
            dict(zip(counts.occupancy, counts.n_sites)), {1: 270, 2: 180, 3: 90}
        )
        self.assertEqual(len(list((self.out / "final").glob("Fig*.png"))), 4)
        self.assertEqual(len(list((self.out / "final").glob("Table*.tsv"))), 6)
        for figure in (self.out / "final").glob("Fig*.png"):
            self.assertGreater(figure.stat().st_size, 10000)

    def test_vectors_are_aligned_after_shuffled_responses(self):
        rows = pd.read_parquet(self.out / "atlas/rows.parquet")
        matrix = zarr.open(str(self.out / "atlas/vectors.zarr"), mode="r")
        for i in [0, 199, len(rows) - 1]:
            row = rows.iloc[i]
            v = genome.Variant(row.chrom, int(row.pos), row.ref, row.alt)
            np.testing.assert_allclose(matrix[i, :6], vector(v, "ATAC"))
            np.testing.assert_allclose(matrix[i, 6:12], vector(v, "DNASE"))

    def test_no_primary_label_leakage_and_cross_fitted_scores(self):
        g = pd.read_parquet(self.out / "results/geometry_candidates.parquet")
        metrics = json.loads((self.out / "results/stage2_metrics.json").read_text())
        self.assertEqual(metrics["conditional"]["n_candidates"], 360)
        self.assertEqual(metrics["conditional"]["n_loci"], 180)
        self.assertEqual(metrics["conditional"]["status"], "ok", metrics)
        self.assertEqual(metrics["all_candidates"]["status"], "ok", metrics)
        self.assertEqual(metrics["site_cv_status"], "ok", metrics)
        primary = g[g.is_primary == 1]
        np.testing.assert_allclose(primary.affinity_to_primary, 1, atol=1e-10)
        np.testing.assert_allclose(primary.novelty_perpendicular, 0, atol=1e-10)
        m = pd.read_parquet(self.out / "results/stage1_candidates.parquet")
        self.assertTrue(
            m.loc[m.occupancy < 3, "mutation_prediction_source"].eq("out_of_fold").all()
        )
        self.assertTrue(np.isfinite(m.mutation_opportunity).all())

    def test_single_chromosome_stage1_uses_locus_folds(self):
        subset = self.base / "single.parquet"
        candidates = pd.read_parquet(self.out / "work/candidates.parquet")
        candidates[candidates.chrom == "chr1"].to_parquet(subset, index=False)
        dest = self.base / "single_stage1"
        result = command(
            sys.executable,
            ROOT / "f/stage1_mutation.py",
            "--candidates",
            subset,
            "--out-dir",
            dest,
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        metrics = json.loads((dest / "stage1_metrics.json").read_text())
        self.assertEqual(metrics["cv_group"], "site_id")
        self.assertEqual(metrics["cv_status"], "ok")

    def test_geometry_rejects_nonfinite_vectors(self):
        # Alter a copy, preserving the successful full pipeline artifact.
        import shutil

        target = self.base / "bad_atlas"
        shutil.copytree(self.out / "atlas", target)
        matrix = zarr.open(str(target / "vectors.zarr"), mode="a")
        matrix[0, 0] = np.nan
        result = command(
            sys.executable,
            ROOT / "f/geometry.py",
            "--candidates",
            self.out / "work/candidates.parquet",
            "--atlas-sites",
            self.out / "work/atlas_sites.parquet",
            "--atlas-dir",
            target,
            "--out-dir",
            self.base / "bad_geometry",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("nonfinite", result.stdout)

    def test_single_chromosome_site_cv_is_explicitly_unavailable(self):
        dest = self.base / "single_stage2"
        geometry = pd.read_parquet(self.out / "results/geometry_candidates.parquet")
        sites = pd.read_parquet(self.out / "results/geometry_sites.parquet")
        geometry_path = self.base / "single_geometry.parquet"
        sites_path = self.base / "single_geometry_sites.parquet"
        geometry[geometry.chrom == "chr1"].to_parquet(geometry_path, index=False)
        sites[sites.chrom == "chr1"].to_parquet(sites_path, index=False)
        result = command(
            sys.executable,
            ROOT / "f/stage2_retention.py",
            "--geometry",
            geometry_path,
            "--sites",
            sites_path,
            "--mutation",
            self.out / "results/stage1_candidates.parquet",
            "--out-dir",
            dest,
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        metrics = json.loads((dest / "stage2_metrics.json").read_text())
        self.assertIsNone(metrics["site_cv_auc"])
        self.assertEqual(metrics["site_cv_status"], "fewer than two validation groups")
        self.assertTrue(
            pd.read_parquet(dest / "stage2_sites.parquet")
            .pred_multiallelic.isna()
            .all()
        )


if __name__ == "__main__":
    unittest.main()
