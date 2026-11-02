"""End-to-end AF-only window query + genuine indexed Atlas + safe exact join."""
from __future__ import annotations

from contextlib import redirect_stdout
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
import allofus_shards as aou
import atlas_scores as atlas
import atlas_allofus as join
from analysis_integrity import inventory

HAVE_PYSAM = importlib.util.find_spec("pysam") is not None
if HAVE_PYSAM:
    import pysam


@unittest.skipUnless(HAVE_PYSAM, "pysam required for real upstream BGZF/Tabix fixtures")
class AtlasAllOfUsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.targets = self.root / "targets.tsv"
        self.targets.write_text("site_id\tchrom\tpos\tref\talts\tbuild\ns1\t1\t10\tA\tC\tGRCh38\n")
        source = self.root / "allofus.gz"
        with gzip.open(source, "wt") as handle:
            handle.write("#" + "\t".join(aou.ANNOVAR_COLUMNS) + "\n")
            handle.write("chr1\t10\t10\tA\tC\t0.1\t0.4\t0.2\t0.1\t0.3\t0.1\t\t0.4\t0.0\n")
            handle.write("chr1\t10\t10\tA\tG\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n")
        plan = aou.make_plan(source, "1", 300000000)
        plan_path = self.root / "plan.json"
        aou.write_plan(plan, plan_path)
        prepared, normalized = self.root / "prepared", self.root / "normalized"
        self.aou_dir, self.atlas_dir, self.out = self.root / "query", self.root / "atlas", self.root / "joined"
        with redirect_stdout(io.StringIO()):
            aou.prepare(plan_path, prepared, 0)
            for shard in plan["shards"]:
                aou.normalize_shard(prepared, shard["id"], normalized / shard["id"])
            aou.query(prepared, normalized, self.targets, self.aou_dir)
        self.avi = self.root / "avi.tsv.gz"
        plain = self.root / "avi.tsv"
        plain.write_text("chr1\t10\tA\tC\t-0.125\t0\nchr1\t10\tA\tG\t0\t1\nchr1\t10\tA\tT\t0.4\t10\n")
        pysam.tabix_compress(str(plain), str(self.avi), force=True)
        pysam.tabix_index(str(self.avi), seq_col=0, start_col=1, end_col=1, force=True, zerobased=False)
        self.make_atlas(self.targets, self.atlas_dir)

    def tearDown(self):
        self.temp.cleanup()

    def make_atlas(self, targets, out):
        args = atlas.parser().parse_args(["--targets", str(targets), "--build", "38", "--avi-file", str(self.avi), "--out", str(out)])
        with redirect_stdout(io.StringIO()):
            return atlas.run(args)

    def execute(self):
        with redirect_stdout(io.StringIO()):
            return join.run(self.aou_dir, self.atlas_dir, self.out)

    def rows(self, path=None):
        with gzip.open(path or self.out / "allofus_atlas.tsv.gz", "rt") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def mutate_manifest(self, directory, action):
        path = directory / "manifest.json"
        manifest = json.loads(path.read_text())
        action(manifest)
        path.write_text(json.dumps(manifest))

    def rewrite_table(self, directory, filename, transform):
        path = directory / filename
        rows = transform(self.rows(path))
        with gzip.open(path, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        self.mutate_manifest(directory, lambda manifest: manifest.update(output_inventory=inventory(directory, "manifest.json")))

    def test_exact_join_retains_extra_alt_all_populations_and_missingness(self):
        result = self.execute()
        rows = {row["alt"]: row for row in self.rows()}
        self.assertEqual(set(rows), set("ACGT"))
        self.assertEqual(result["counts"]["rows"], 4)
        self.assertEqual(result["counts"]["substitutions_with_AF_and_AVI"], 2)
        self.assertEqual(rows["A"]["allele_role"], "REF")
        self.assertEqual(rows["A"]["avi_raw"], "NA")
        self.assertTrue(all(rows["A"][col] == "NA" for col in [*aou.AF_COLUMNS, "gvs_max_af"]))
        self.assertEqual(rows["C"]["observed_from_targets"], "1")
        self.assertEqual(rows["C"]["gvs_mid_af"], "NA")
        self.assertEqual(rows["C"]["gvs_max_af"], "0.4")
        self.assertEqual(rows["C"]["avi_raw"], "-0.125")
        self.assertEqual(rows["C"]["avi_phred"], "0")
        self.assertEqual(rows["G"]["gvs_all_af"], "0.0")
        self.assertEqual(rows["G"]["allofus_additional_reported_alt"], "1")
        self.assertEqual(rows["G"]["observed_from_targets"], "0")
        self.assertEqual(rows["T"]["allofus_availability"], "no_exact_alt_row_in_query")
        self.assertEqual(rows["T"]["gvs_all_af"], "NA")
        self.assertEqual(rows["T"]["avi_raw"], "0.4")
        self.assertEqual(rows["T"]["allofus_observed_in_target"], "NA")
        self.assertFalse({"AC", "AN", "ref_AF", "pairwise_distance"} & set(rows["C"]))
        self.assertEqual(len(aou.AF_COLUMNS), 8)

    def test_identical_target_content_at_new_path_is_allowed(self):
        copied = self.root / "targets_copy.tsv"
        copied.write_bytes(self.targets.read_bytes())
        replacement = self.root / "atlas_copy"
        self.make_atlas(copied, replacement)
        self.atlas_dir = replacement
        self.assertEqual(self.execute()["counts"]["rows"], 4)

    def test_different_target_content_rejected(self):
        other = self.root / "other_targets.tsv"
        other.write_text(self.targets.read_text().replace("\tC\tGRCh38", "\tG\tGRCh38"))
        replacement = self.root / "atlas_other"
        self.make_atlas(other, replacement)
        self.atlas_dir = replacement
        with self.assertRaisesRegex(ValueError, "Target-file content identities differ"):
            self.execute()

    def test_missing_target_sha_rejected(self):
        self.mutate_manifest(self.aou_dir, lambda m: m["targets"].pop("sha256"))
        with self.assertRaisesRegex(ValueError, "SHA256 identities"):
            self.execute()

    def test_incomplete_source_manifest_and_wrong_schema_rejected(self):
        self.mutate_manifest(self.atlas_dir, lambda m: m.update(status="running"))
        with self.assertRaisesRegex(ValueError, "Nonterminal"):
            self.execute()
        self.mutate_manifest(self.atlas_dir, lambda m: m.update(status="complete", schema="different"))
        with self.assertRaisesRegex(ValueError, "completed atlas_scores"):
            self.execute()

    def test_output_inventory_tampering_rejected(self):
        with (self.aou_dir / "target_alleles.tsv.gz").open("ab") as handle:
            handle.write(b"tamper")
        with self.assertRaisesRegex(ValueError, "Output identity changed"):
            self.execute()

    def test_build_mismatch_rejected(self):
        self.mutate_manifest(self.aou_dir, lambda m: m.update(build="GRCh37"))
        with self.assertRaisesRegex(ValueError, "Both completed sources"):
            self.execute()

    def test_source_ref_mismatch_is_not_silently_merged(self):
        def change(rows):
            rows[0]["source_ref"] = "G"
            return rows
        self.rewrite_table(self.aou_dir, "target_alleles.tsv.gz", change)
        with self.assertRaisesRegex(ValueError, "source REF conflicts"):
            self.execute()

    def test_wrong_site_id_cannot_be_joined_by_position_only(self):
        def change(rows):
            rows[0]["site_id"] = "other"
            return rows
        self.rewrite_table(self.aou_dir, "target_alleles.tsv.gz", change)
        with self.assertRaisesRegex(ValueError, "site_id or REF mismatch"):
            self.execute()

    def test_duplicate_allele_rows_rejected(self):
        self.rewrite_table(self.aou_dir, "target_alleles.tsv.gz", lambda rows: rows + [rows[0]])
        with self.assertRaisesRegex(ValueError, "unique canonical"):
            self.execute()

    def test_missing_theoretical_atlas_row_rejected(self):
        self.rewrite_table(self.atlas_dir, "atlas_variants.tsv.gz", lambda rows: [r for r in rows if r["allele"] != "T"])
        with self.assertRaisesRegex(ValueError, "required theoretical-base"):
            self.execute()

    def test_misleading_target_membership_flag_rejected(self):
        def change(rows):
            rows[0]["requested_alt"] = "0"
            rows[0]["observed_in_target"] = "0"
            rows[0]["additional_reported_alt"] = "1"
            return rows
        self.rewrite_table(self.aou_dir, "target_alleles.tsv.gz", change)
        with self.assertRaisesRegex(ValueError, "membership flags disagree"):
            self.execute()

    def test_invalid_af_and_fake_reference_score_rejected(self):
        def bad_af(rows):
            rows[0]["gvs_all_af"] = "1.01"
            return rows
        self.rewrite_table(self.aou_dir, "target_alleles.tsv.gz", bad_af)
        with self.assertRaisesRegex(ValueError, "Invalid gvs_all_af"):
            self.execute()
        def bad_reference(rows):
            next(r for r in rows if r["allele_role"] == "REF")["avi_raw"] = "0"
            return rows
        self.rewrite_table(self.atlas_dir, "atlas_variants.tsv.gz", bad_reference)
        with self.assertRaisesRegex(ValueError, "REF must retain NA"):
            self.execute()

    def test_protected_source_directory_cannot_be_output(self):
        with self.assertRaisesRegex(ValueError, "overlaps protected"):
            join.run(self.aou_dir, self.atlas_dir, self.aou_dir)
        self.assertEqual(json.loads((self.aou_dir / "manifest.json").read_text())["status"], "complete_descriptive_unverified")

    def test_completed_no_score_audit_preserves_af_and_na_predictions(self):
        plain = self.root / "avi.tsv"
        plain.write_text("chr1\t11\tC\tA\t0\t0\n")
        pysam.tabix_compress(str(plain), str(self.avi), force=True)
        pysam.tabix_index(str(self.avi), seq_col=0, start_col=1, end_col=1, force=True, zerobased=False)
        self.assertEqual(self.make_atlas(self.targets, self.atlas_dir)["status"], "not_available")
        result = self.execute()
        self.assertEqual(result["status"], "complete_with_exclusions")
        self.assertEqual(result["counts"]["substitutions_with_AF_and_AVI"], 0)
        row = next(row for row in self.rows() if row["alt"] == "C")
        self.assertEqual(row["gvs_all_af"], "0.1")
        self.assertEqual(row["avi_raw"], "NA")
        self.assertEqual(row["atlas_status"], "not_scored")


if __name__ == "__main__":
    unittest.main()
