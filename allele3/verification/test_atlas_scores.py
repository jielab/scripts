"""Real BGZF + Tabix checks for bounded, exact-allele Atlas lookups."""
from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
import atlas_scores as atlas

HAVE_PYSAM = importlib.util.find_spec("pysam") is not None
if HAVE_PYSAM:
    import pysam


@unittest.skipUnless(HAVE_PYSAM, "pysam required: these tests use actual BGZF and TBI/CSI files")
class AtlasScoresTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.targets = self.root / "targets.tsv"
        self.data = self.root / "avi.tsv.gz"
        self.out = self.root / "analysis" / "atlas"
        self.write_targets([("s1", "1", 10, "A", "C", "GRCh38")])
        self.write_atlas([
            ("chr1", 10, "A", "C", "-0.125", "0"),
            ("chr1", 10, "A", "G", "0", "1.2"),
            ("chr1", 10, "A", "T", "0.6", "20"),
            ("chr1", 11, "C", "A", "0.1", "2"),
        ])

    def tearDown(self):
        self.temp.cleanup()

    def write_targets(self, rows, columns=("site_id", "chrom", "pos", "ref", "alts", "build")):
        with self.targets.open("w", newline="") as handle:
            out = csv.writer(handle, delimiter="\t", lineterminator="\n")
            out.writerow(columns)
            out.writerows(rows)

    def indexed(self, path, rows, header="", zero_based=False, csi=False, seq_col=0, start_col=1, end_col=1):
        plain = path.with_suffix("")
        plain.write_text(header + "".join("\t".join(map(str, row)) + "\n" for row in rows))
        for suffix in (".tbi", ".csi"):
            Path(str(path) + suffix).unlink(missing_ok=True)
        pysam.tabix_compress(str(plain), str(path), force=True)
        pysam.tabix_index(str(path), seq_col=seq_col, start_col=start_col, end_col=end_col,
                          zerobased=zero_based, force=True, csi=csi)
        return path

    def write_atlas(self, rows, **kwargs):
        return self.indexed(self.data, rows, **kwargs)

    def args(self, *more):
        return atlas.parser().parse_args(["--targets", str(self.targets), "--build", "GRCh38", "--avi-file", str(self.data),
                                         "--out", str(self.out), "--source-version", "fixture-2026", *map(str, more)])

    def run_atlas(self, *more):
        return atlas.run(self.args(*more))

    def test_default_data_root_follows_execution_environment(self):
        from unittest.mock import patch
        for root in ('/mnt/f', '/data/sph-huangj'):
            with self.subTest(root=root), patch.dict(os.environ, DIR_DATA=root):
                self.assertEqual(self.args().dir_data, root)

    def test_missing_explicit_index_reports_its_exact_path(self):
        missing=self.root/'matching-index.csi'
        with self.assertRaises(ValueError) as caught:
            self.run_atlas('--avi-index',missing)
        self.assertIn(str(missing),str(caught.exception))
        self.assertFalse(self.out.exists())

    def rows(self, filename="atlas_variants.tsv.gz"):
        with gzip.open(self.out / filename, "rt") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def variant(self, allele):
        return next(row for row in self.rows() if row["allele"] == allele)

    def custom_config(self, changes):
        with atlas.DEFAULT_RESOURCES.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = reader.fieldnames, list(reader)
        for rid, values in changes.items():
            next(row for row in rows if row["resource_id"] == rid).update(values)
        path = self.root / "resources.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_negative_zero_all_four_bases_and_default_theoretical_alts(self):
        result = self.run_atlas()
        self.assertEqual(result["status"], "complete")
        self.assertEqual([row["allele"] for row in self.rows()], list("ACGT"))
        self.assertEqual(self.variant("A")["avi_raw"], "NA")
        self.assertEqual(self.variant("A")["observed_from_targets"], "NA")
        self.assertEqual(self.variant("C")["avi_raw"], "-0.125")
        self.assertEqual(self.variant("C")["avi_phred"], "0")
        self.assertEqual(self.variant("G")["avi_raw"], "0")
        self.assertEqual(self.variant("C")["observed_from_targets"], "1")
        self.assertEqual(self.variant("G")["observed_from_targets"], "0")
        self.assertEqual(self.variant("T")["observed_from_targets"], "0")
        self.assertEqual(result["query_audits"]["avi"]["returned_rows"], 3)
        self.assertEqual(result["query_audits"]["avi"]["tabix_fetches"], 1)
        self.assertEqual(result["counts"]["metric_rows"], 8)
        self.assertNotIn("pairwise", "\t".join(self.rows()[0]))

    def test_missing_row_and_missing_metrics_are_not_zero(self):
        self.write_atlas([("1", 10, "A", "C", "0", ""), ("1", 10, "A", "G", ".", "NA")])
        result = self.run_atlas()
        self.assertEqual(result["status"], "complete_with_exclusions")
        self.assertEqual(self.variant("C")["avi_raw"], "0")
        self.assertEqual(self.variant("C")["avi_phred"], "NA")
        self.assertEqual(self.variant("C")["status"], "available_partial")
        self.assertEqual(self.variant("G")["status"], "score_missing")
        self.assertEqual(self.variant("T")["status"], "not_scored")
        self.assertEqual(self.variant("T")["avi_raw"], "NA")

    def test_target_ref_mismatch_fails_and_cannot_publish_complete_manifest(self):
        self.write_targets([("s1", "1", 10, "G", "C", "38")])
        with self.assertRaisesRegex(ValueError, "REF mismatch"):
            self.run_atlas()
        manifest = json.loads((self.out / "manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")
        self.assertFalse((self.out / "atlas_variants.tsv.gz").exists())

    def test_target_build_is_checked_before_any_publication(self):
        self.write_targets([("s1", "1", 10, "A", "C", "GRCh37")])
        with self.assertRaisesRegex(ValueError, "conflicts with GRCh38"):
            self.run_atlas()
        self.assertFalse(self.out.exists())
        with self.assertRaisesRegex(ValueError, "require GRCh38"):
            self.run_atlas("--build", "GRCh37")

    def test_all_build_columns_must_agree(self):
        self.write_targets([("s1", "1", 10, "A", "C", "hg38", "37")],
                           columns=("site_id", "chrom", "pos", "ref", "alts", "build", "assembly"))
        with self.assertRaisesRegex(ValueError, "assembly=37"):
            self.run_atlas()

    def test_comment_header_and_source_build_checked(self):
        rows = [("chr1", 10, "A", "C", "-1", "2")]
        self.write_atlas(rows, header="##reference=GRCh38\n#CHROM\tPOS\tREF\tALT\traw_score\tPHRED\n")
        result = self.run_atlas()
        self.assertEqual(result["resources"][0]["header_check"]["verification"], "GRCh38_header_and_declared")
        self.write_atlas(rows, header="##reference=GRCh37\n#CHROM\tPOS\tREF\tALT\traw_score\tPHRED\n")
        with self.assertRaisesRegex(ValueError, "source header build"):
            self.run_atlas()

    def test_header_columns_cannot_be_silently_reordered(self):
        self.write_atlas([("1", 10, "A", "C", "1", "2")],
                         header="#CHROM\tPOS\tREF\tALT\tPHRED\traw_score\n")
        with self.assertRaisesRegex(ValueError, "column header conflicts"):
            self.run_atlas()

    def test_zero_based_index_rejected_before_point_query(self):
        self.write_atlas([("1", 10, "A", "C", "1", "2")], zero_based=True)
        with self.assertRaisesRegex(ValueError, "index uses 0-based"):
            self.run_atlas()

    def test_zero_based_header_rejected_even_if_index_is_one_based(self):
        self.write_atlas([("1", 10, "A", "C", "1", "2")], header="##coordinate_system=0-based\n")
        with self.assertRaisesRegex(ValueError, "header declares zero-based"):
            self.run_atlas()

    def test_csi_index_is_supported_and_checked(self):
        self.write_atlas([("1", 10, "A", "C", "-1", "0")], csi=True)
        result = self.run_atlas()
        self.assertEqual(result["resources"][0]["index_check"]["format"], "CSI")
        self.assertEqual(self.variant("C")["avi_raw"], "-1")

    def test_ambiguous_index_requires_explicit_selection(self):
        pysam.tabix_index(str(self.data), seq_col=0, start_col=1, end_col=1, force=True, csi=True)
        self.assertTrue(Path(str(self.data) + ".tbi").exists())
        with self.assertRaisesRegex(ValueError, "Both .tbi and .csi"):
            self.run_atlas()
        self.assertEqual(self.run_atlas("--avi-index", str(self.data) + ".csi")["status"], "complete")

    def test_negative_phred_and_nonfinite_raw_rejected(self):
        for raw, phred, message in [("0", "-1", "invalid PHRED"), ("inf", "1", "invalid raw_score")]:
            with self.subTest(raw=raw, phred=phred):
                self.write_atlas([("1", 10, "A", "C", raw, phred)])
                with self.assertRaisesRegex(ValueError, message):
                    self.run_atlas()

    def test_identical_duplicates_deduplicated_but_conflicts_fail(self):
        rows = [("1", 10, "A", "C", "1", "2"), ("1", 10, "A", "C", "1.0", "2.0")]
        self.write_atlas(rows)
        result = self.run_atlas()
        self.assertEqual(result["query_audits"]["avi"]["identical_duplicate_rows"], 1)
        self.assertEqual(len(self.rows()), 4)
        self.write_atlas(rows + [("1", 10, "A", "C", "3", "2")])
        with self.assertRaisesRegex(ValueError, "conflicting duplicate"):
            self.run_atlas()

    def test_optional_generic_metrics_and_missing_optional_file(self):
        metric_path = self.root / "metric.tsv.gz"
        self.indexed(metric_path, [("1", 10, "A", "C", "-0.4", "")],
                     header="#CHROM\tPOS\tREF\tALT\tsigned_effect\toptional_effect\n")
        schema = dict(kind="metrics", path=str(metric_path), columns="CHROM,POS,REF,ALT,signed_effect,optional_effect",
                      chrom_column="CHROM", pos_column="POS", ref_column="REF", alt_column="ALT",
                      metric_columns="signed_effect,optional_effect", enabled="1", required="0")
        missing_schema = dict(schema, path=str(self.root / "not_downloaded.tsv.gz"))
        config = self.custom_config({"splicing": schema, "avi_attributions": missing_schema})
        result = self.run_atlas("--resources", config)
        self.assertEqual(result["status"], "complete_with_exclusions")
        metrics = [row for row in self.rows("atlas_metrics.tsv.gz") if row["resource_id"] == "splicing" and row["allele"] == "C"]
        self.assertEqual([row["value"] for row in metrics], ["-0.4", "NA"])
        self.assertEqual([row["status"] for row in metrics], ["available", "score_missing"])
        self.assertEqual(result["resources"][2]["status"], "deferred_missing_file")

    def test_optional_existing_resource_wrong_build_not_ignored(self):
        schema = dict(path=str(self.data), build="GRCh37", columns="CHROM,POS,REF,ALT,raw_score,PHRED",
                      chrom_column="CHROM", pos_column="POS", ref_column="REF", alt_column="ALT",
                      metric_columns="raw_score,PHRED", enabled="1")
        config = self.custom_config({"splicing": schema})
        with self.assertRaisesRegex(ValueError, "GRCh38 1-based"):
            self.run_atlas("--resources", config)

    def test_targets_can_omit_site_id_and_build_but_need_cli_build(self):
        self.write_targets([("chr1", 10, "A", "C,T")], columns=("chrom", "pos", "ref", "alts"))
        self.run_atlas()
        self.assertEqual(self.rows()[0]["site_id"], "1:10:A")
        self.assertEqual(self.variant("T")["observed_from_targets"], "1")

    def test_chromosome_region_filters_and_xy_supported(self):
        self.write_atlas([("chr1", 10, "A", "C", "0", "1"), ("chrX", 30, "T", "A", "1", "2"),
                          ("chrY", 40, "C", "A", "2", "3")])
        self.write_targets([("s1", "1", 10, "A", "C", "38"), ("sx", "X", 30, "T", "A", "38"), ("sy", "Y", 40, "C", "A", "38")])
        result = self.run_atlas("--chromosomes", "X,Y", "--region", "chrY:40-40")
        self.assertEqual(result["counts"]["target_sites"], 1)
        self.assertEqual({row["chrom"] for row in self.rows()}, {"Y"})
        self.assertEqual(self.rows()[0]["allele"], "C")
        with self.assertRaisesRegex(ValueError, "No target sites"):
            self.run_atlas("--region", "1:100-200")

    def test_absent_contig_is_explicit_and_main_returns_unavailable(self):
        self.write_targets([("s2", "2", 10, "A", "C", "38")])
        result = self.run_atlas()
        self.assertEqual(result["status"], "not_available")
        self.assertEqual(self.variant("C")["status"], "not_indexed_contig")
        rc = atlas.main(["--targets", str(self.targets), "--build", "38", "--avi-file", str(self.data), "--out", str(self.out)])
        self.assertEqual(rc, 2)

    def test_reuse_tracks_input_and_output_identity(self):
        first = self.run_atlas()
        second = self.run_atlas()
        self.assertFalse(first["reused"])
        self.assertTrue(second["reused"])
        self.assertEqual(first["run_id"], second["run_id"])
        self.write_targets([("s1", "1", 10, "A", "C,T", "38")])
        changed = self.run_atlas()
        self.assertFalse(changed["reused"])
        self.assertEqual(self.variant("T")["observed_from_targets"], "1")
        with (self.out / "target_availability.tsv").open("a") as handle:
            handle.write("tampered\n")
        self.assertFalse(self.run_atlas()["reused"])

    def test_source_change_invalidates_reuse(self):
        self.run_atlas()
        self.write_atlas([("chr1", 10, "A", "C", "8", "20")])
        result = self.run_atlas()
        self.assertFalse(result["reused"])
        self.assertEqual(self.variant("C")["avi_raw"], "8")

    def test_configured_data_root_is_not_a_protected_input_directory(self):
        path = self.root / "annot" / "allele3" / "alphagenome" / "hg38" / "alphagenome_variant_impact_score_snvs.tsv.gz"
        path.parent.mkdir(parents=True)
        self.indexed(path, [("1", 10, "A", "C", "-1", "0")])
        args = atlas.parser().parse_args(["--targets", str(self.targets), "--build", "38", "--dir-data", str(self.root), "--out", str(self.out)])
        self.assertEqual(atlas.run(args)["status"], "complete_with_exclusions")
        self.assertTrue((self.out / "manifest.json").is_file())

    def test_output_cannot_overlap_input_file_directory(self):
        with self.assertRaisesRegex(ValueError, "overlaps protected"):
            self.run_atlas("--out", self.root)
        self.assertFalse((self.root / "manifest.json").exists())
        self.assertTrue(self.targets.is_file())
        self.assertTrue(self.data.is_file())

    def test_source_ref_equals_alt_rejected(self):
        self.write_atlas([("1", 10, "A", "A", "0", "0")])
        with self.assertRaisesRegex(ValueError, "not a canonical nonreference SNV"):
            self.run_atlas()

    def test_ambiguous_contig_aliases_rejected(self):
        self.write_atlas([("1", 10, "A", "C", "0", "0"), ("chr1", 10, "A", "G", "1", "2")])
        with self.assertRaisesRegex(ValueError, "ambiguous indexed contig aliases"):
            self.run_atlas()

    def test_selected_target_cap_requires_sharding(self):
        self.write_targets([("s1", "1", 10, "A", "C", "38"), ("s2", "1", 11, "C", "A", "38")])
        with self.assertRaisesRegex(ValueError, "Target cap 1 exceeded"):
            self.run_atlas("--max-targets", "1")
        self.assertFalse(self.out.exists())

    def test_public_shell_dispatches_atlas_help(self):
        root = Path(__file__).resolve().parents[1]
        result = subprocess.run(["bash", str(root / "allele3.sh"), "research", "atlas", "--help", "--python", sys.executable],
                                cwd=root, text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--avi-file", result.stdout)
        self.assertIn("--resources", result.stdout)


if __name__ == "__main__":
    unittest.main()
