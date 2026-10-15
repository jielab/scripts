"""Protect cohort completeness, streamed boundaries and exact weighted training."""

from pathlib import Path
import os
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
import pysam
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
from stage1_mutation import CAT, NUM, CountModel, site_bucket
from vcf_inputs import resolve_inputs, validate_inputs, GNOMAD_PREFIX
from test_pipeline import HEADER, command


class StreamingTests(unittest.TestCase):
    def test_weighted_feature_counts_match_expanded_logistic_model(self):
        rng = np.random.default_rng(22)
        rows = []
        for context in ("AAA", "CAC", "GAT"):
            for alt in ("C", "G", "T"):
                for observed in (0, 1):
                    rows.append(
                        dict(
                            context=context,
                            ref_alt=f"A>{alt}",
                            transition=int(alt == "G"),
                            cpg_transition=0,
                            observed=observed,
                            count=int(rng.integers(1, 40)),
                        )
                    )
        compact = pd.DataFrame(rows)
        expanded = compact.loc[compact.index.repeat(compact["count"])].reset_index(
            drop=True
        )
        reference = make_pipeline(
            ColumnTransformer(
                [
                    ("cat", OneHotEncoder(handle_unknown="ignore"), CAT),
                    ("num", StandardScaler(), NUM),
                ]
            ),
            LogisticRegression(max_iter=1000, class_weight="balanced", random_state=22),
        )
        reference.fit(expanded[CAT + NUM], expanded.observed)
        actual = CountModel(22).fit(compact).predict(compact)
        np.testing.assert_allclose(
            actual, reference.predict_proba(compact[CAT + NUM])[:, 1], atol=1e-8
        )
        ids = pd.Series(["chr1:10:A"] * 3 + ["chr1:20:C"] * 3)
        folds = site_bucket(ids, 22)
        self.assertEqual(len(set(folds[:3])), 1)
        np.testing.assert_array_equal(folds, site_bucket(ids, 22))

    def test_gnomad_requires_requested_set_and_excludes_unrelated_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            file = base / f"{GNOMAD_PREFIX}22.vcf.bgz"
            file.write_bytes(b"partial")
            (base / "another_cohort.vcf").write_text(HEADER)
            files = resolve_inputs(directory=base)
            self.assertEqual(len(files), 22)
            self.assertTrue(all(p.name.startswith(GNOMAD_PREFIX) for p in files))
            with self.assertRaisesRegex(ValueError, "Missing VCF"):
                validate_inputs(files)
            self.assertEqual(resolve_inputs(directory=base, chrs="22"), [file])
            with self.assertRaisesRegex(ValueError, "wrong size"):
                validate_inputs([file])
            self.assertEqual(len(resolve_inputs(directory=base, chrs="all")), 24)

    def test_aria2_and_truncated_bgzf_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            source = base / "chr1.vcf"
            source.write_text(HEADER + "chr1\t10\t.\tA\tC\t.\tPASS\tAC=1;AN=10\n")
            compressed = base / "chr1.vcf.gz"
            pysam.tabix_compress(str(source), str(compressed))
            pysam.tabix_index(str(compressed), preset="vcf")
            validate_inputs([compressed], require_index=True)
            marker = Path(str(compressed) + ".aria2")
            marker.touch()
            with self.assertRaisesRegex(ValueError, "Download incomplete"):
                validate_inputs([compressed])
            marker.unlink()
            compressed.write_bytes(compressed.read_bytes()[:-10])
            with self.assertRaisesRegex(ValueError, "end marker"):
                validate_inputs([compressed])

    def test_streamed_batch_boundaries_merge_files_and_refuse_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            fasta = base / "ref.fa"
            fasta.write_text(">chr1\n" + "A" * 6000 + "\n")
            pysam.faidx(str(fasta))
            (base / "a.vcf").write_text(
                HEADER
                + "chr1\t10\t.\tA\tC\t.\tPASS\tAC=1;AN=10\n"
                + "chr1\t20\t.\tA\tT\t.\tPASS\tAC=2;AN=10\n"
            )
            (base / "b.vcf").write_text(
                HEADER + "chr1\t10\t.\tA\tG\t.\tPASS\tAC=1;AN=10\n"
            )
            output = base / "candidates.parquet"
            args = [
                sys.executable,
                ROOT / "f/extract_sites.py",
                "--vcf-dir",
                base,
                "--fasta",
                fasta,
                "--batch-sites",
                "1",
                "--out",
                output,
                "--sites-out",
                base / "sites.parquet",
            ]
            result = command(*args)
            self.assertEqual(result.returncode, 0, result.stdout)
            data = pd.read_parquet(output)
            self.assertEqual(
                data.groupby("pos").occupancy.first().to_dict(), {10: 2, 20: 1}
            )
            self.assertEqual(len(data), 6)
            original = output.read_bytes()
            result = command(*args)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("already exists", result.stdout)
            self.assertEqual(output.read_bytes(), original)

    def test_unsorted_vcf_does_not_publish_partial_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            fasta = base / "ref.fa"
            fasta.write_text(">chr1\n" + "A" * 6000 + "\n")
            pysam.faidx(str(fasta))
            source = base / "input.vcf"
            source.write_text(
                HEADER
                + "".join(
                    f"chr1\t{pos}\t.\tA\tC\t.\tPASS\tAC=1;AN=10\n"
                    for pos in [10, 20, 5]
                )
            )
            output = base / "candidates.parquet"
            result = command(
                sys.executable,
                ROOT / "f/extract_sites.py",
                "--vcf",
                source,
                "--fasta",
                fasta,
                "--batch-sites",
                "1",
                "--out",
                output,
                "--sites-out",
                base / "sites.parquet",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("coordinate sorted", result.stdout)
            self.assertFalse(output.exists())
            self.assertFalse((base / "sites.parquet").exists())
            self.assertFalse(list(base.glob("*.tmp")))

    def test_check_input_does_not_create_analysis_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            fasta = base / "ref.fa"
            fasta.write_text(">chr1\nAAAA\n")
            pysam.faidx(str(fasta))
            source = base / "input.vcf"
            source.write_text(HEADER)
            result = command(
                ROOT / "allele.sh",
                "check-input",
                "--vcf",
                source,
                "--fasta",
                fasta,
                "--out",
                base / "unused",
                env={
                    **os.environ,
                    "ALLELE_CONFIG": "/dev/null",
                    "ALLELE_ENV_NAME": "allele",
                },
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("Inputs ready", result.stdout)
            self.assertFalse((base / "unused").exists())


if __name__ == "__main__":
    unittest.main()
