"""Independent failure cases and an indexed-VCF end-to-end regression."""
import csv
import gzip
import io
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"f"))
import archaic_alleles as aa
import scan
from herv_archaic import inspect_dataset

try:
    import pysam
except ImportError:
    pysam = None


class Call(dict):
    def __init__(self, gt, phased=False, **kwargs):
        super().__init__(GT=gt, **kwargs); self.phased = phased


class Record:
    def __init__(self, alts, gt, phased=False, **kwargs):
        self.ref = "A"; self.alleles = ("A", *alts)
        self.samples = {"s": Call(gt, phased, **kwargs)}; self.filter = {"PASS": True}


class EvidenceTests(unittest.TestCase):
    def test_research_dispatcher_help_forwards_arguments(self):
        cp = subprocess.run([sys.executable, str(ROOT/"f/research.py"), "archaic-alleles", "--help"],
                            capture_output=True, text=True)
        self.assertEqual(cp.returncode, 0, cp.stdout+cp.stderr)
        self.assertIn("--target-to-gu-chain", cp.stdout)
        self.assertIn("--emit-all-individuals", cp.stdout)

    def test_unencoded_third_base_is_missing_not_zero(self):
        row = aa.decode_base([Record(("C",), (0, 1))], "A", "G", "s")
        self.assertEqual(row["status"], "not_encoded")
        self.assertIsNone(row["dosage"])

    def test_target_phasing_does_not_phase_ibdmix(self):
        row = aa.decode_base([Record(("C", "G"), (0, 2), True)], "A", "G", "s")
        self.assertTrue(row["target_gt_phased"])
        self.assertEqual(aa.evidence(row, True, {"dosage": 1}), "cooccurrence_and_archaic_base_match")

    def test_ploidy_and_quality_are_observed(self):
        row = aa.decode_base([Record(("C", "G"), (2,))], "A", "G", "s")
        self.assertEqual((row["dosage"], row["ploidy"]), (1, 1))
        self.assertEqual(aa.decode_base([Record(("G",), (1,))], "A", "G", "s", min_gq=20)["status"], "missing_gq")

    def test_conflicting_split_phase_and_legal_trans(self):
        bad = [Record(("C",), (1, 0), True), Record(("G",), (1, 0), True)]
        self.assertEqual(aa.decode_base(bad, "A", "G", "s")["status"], "split_phase_conflict")
        good = [Record(("C",), (1, 0), True), Record(("G",), (0, 1), True)]
        row = aa.decode_base(good, "A", "G", "s")
        self.assertEqual((row["dosage"], row["representation"]), (1, "split_alt_dosage_only"))

    def test_overlapping_tracts_count_the_same_person_once_and_half_open_end(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"normalized.tsv"
            path.write_text("1\tAltai\ts1\t9\t20\tsource\n1\tAltai\ts1\t9\t21\tsource\n1\tAltai\ts2\t0\t9\tsource\n")
            db = sqlite3.connect(":memory:")
            aa.membership_database(path, [{"gu_pos": 10}], db)
            self.assertEqual(db.execute("SELECT * FROM member").fetchall(), [(10, "Altai", "s1")])
            db.close()

    def test_missing_nan_lod_and_fractional_coordinate_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"segments.tsv"
            for value in ("NA", "nan", "inf", ""):
                path.write_text("ID\tchrom\tstart\tend\tanc\tslod\ns1\t1\t0\t60000\tAltai\t"+value+"\n")
                with self.assertRaisesRegex(ValueError, "LOD"):
                    scan.normalize_file(path, "ibdmix", io.StringIO(), 4, 50000)
            path.write_text("ID\tchrom\tstart\tend\tanc\tslod\ns1\t1\t0.5\t60000\tAltai\t5\n")
            with self.assertRaisesRegex(ValueError, "Malformed tract"):
                scan.normalize_file(path, "ibdmix", io.StringIO(), 4, 50000)

    def test_default_reference_selection_does_not_add_historical_pair(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"segments.tsv"
            path.write_text("ID\tchrom\tstart\tend\tanc\tslod\ns1\t1\t0\t60000\tAltai\t5\ns1\t1\t0\t60000\tAltai.2013\t5\n")
            buf = io.StringIO()
            scan.normalize_file(path, "ibdmix", buf, 4, 50000, allowed_references=scan.track_references())
            self.assertEqual(len(buf.getvalue().splitlines()), 1)
            self.assertNotIn(".2013", buf.getvalue())

    def test_caller_versions_cannot_be_mixed_and_untested_reference_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/"1kg"
            refs = " ".join(scan.DEFAULT_ARCHAIC_REFERENCES)
            for chrom in ("1", "2"):
                directory = root/"ibdmix"/f"chr{chrom}"
                (directory/"final").mkdir(parents=True); (directory/"samples").mkdir()
                (directory/"samples/ALL.txt").write_text("s1\ns2\n")
                (directory/".complete").touch()
                with gzip.open(directory/"final/all_archaic_refs.segments.tsv.gz", "wt") as handle:
                    handle.write("ID\tchrom\tstart\tend\tanc\tslod\n")
                (directory/"run.meta.tsv").write_text(f"pipeline_version\tversion{chrom}\ngenome_build\tb37\ncoordinate_system\t0-based-half-open\nrefs\t{refs}\nlod_cut\t4\nlen_cut\t50000\n")
            with self.assertRaisesRegex(ValueError, "contract differs"):
                inspect_dataset(root, chromosomes=("1", "2"))
            with self.assertRaisesRegex(ValueError, "not tested"):
                inspect_dataset(root, chromosomes=("1",), references=("Altai.2013",))

    def test_reverse_liftover_complements_all_bases(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); fake = root/"liftOver"
            fake.write_text("#!/usr/bin/env python3\nimport pathlib,sys\npathlib.Path(sys.argv[-2]).write_text('chr1\\t19\\t20\\tt0\\t0\\t-\\n')\npathlib.Path(sys.argv[-1]).write_text('')\n")
            fake.chmod(0o755)
            targets = [dict(chrom="1", pos=10, ref="A", alts=["C", "G"], site_id="1:10:A")]
            mapped, audit = aa.map_targets(targets, "38", "37", "chain", str(fake), root)
            self.assertEqual((mapped[0]["gu_pos"], mapped[0]["gu_ref"], mapped[0]["gu_alts"]), (20, "T", ["G", "C"]))


@unittest.skipIf(pysam is None, "pysam is required for indexed-VCF integration")
class IndexedIntegration(unittest.TestCase):
    def create_vcf(self, path, people, calls):
        path.parent.mkdir(parents=True, exist_ok=True)
        raw = path.with_suffix("")
        raw.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                       + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(people)+"\n"
                       + "1\t100001\t.\tA\tC,G\t.\tPASS\t.\tGT\t"+"\t".join(calls)+"\n")
        pysam.tabix_compress(str(raw), str(path), force=True)
        pysam.tabix_index(str(path), preset="vcf", force=True)

    def test_actual_cli_preserves_zero_tract_and_missing_genotype_people(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            modern, ancient = root/"vcf/chr1.vcf.gz", root/"archaic/Altai/chr1.vcf.gz"
            self.create_vcf(modern, ["s1", "s2", "s3"], ["0|2", "1|0", "./."])
            self.create_vcf(ancient, ["Altai"], ["2/2"])
            gu = root/"gu/1kg/ibdmix/chr1"; (gu/"samples").mkdir(parents=True); (gu/"final").mkdir()
            (gu/"samples/ALL.txt").write_text("s1\ns2\ns3\n"); (gu/".complete").touch()
            refs = " ".join(scan.DEFAULT_ARCHAIC_REFERENCES)+" Altai.2013"
            (gu/"run.meta.tsv").write_text("pipeline_version\ttest-1\ngenome_build\tb37\ncoordinate_system\t0-based-half-open\nrefs\t"+refs+"\nlod_cut\t4\nlen_cut\t50000\nmodern_vcf\t"+str(modern)+":1:1\narchaic_vcf\t"+str(ancient)+":1:1\n")
            with gzip.open(gu/"final/all_archaic_refs.lod4.len50000.segments.tsv.gz", "wt") as handle:
                handle.write("ID\tchrom\tstart\tend\tanc\tslod\ns1\t1\t100000\t160000\tAltai\t5\ns1\t1\t99000\t161000\tAltai\t6\ns2\t1\t100000\t160000\tAltai.2013\t5\n")
            targets = root/"targets.tsv"; targets.write_text("chrom\tpos\tref\talts\tfamily\tgroup\n1\t100001\tA\tC,G\tERVK\tcase\n")
            out = root/"result"
            cp = subprocess.run([sys.executable, str(ROOT/"f/research.py"), "archaic-alleles", "--gu-dataset", str(root/"gu/1kg"),
                                 "--targets", str(targets), "--target-build", "37", "--modern-build", "37", "--emit-all-individuals", "--out", str(out)], capture_output=True, text=True)
            self.assertEqual(cp.returncode, 0, cp.stdout+cp.stderr)
            manifest = json.loads((out/"manifest.json").read_text())
            self.assertEqual(manifest["n_rows"], 30)
            with gzip.open(out/"individual_evidence.tsv.gz", "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            row = next(r for r in rows if r["reference"] == "Altai" and r["sample_id"] == "s1" and r["allele"] == "G")
            self.assertEqual(row["evidence"], "cooccurrence_and_archaic_base_match")
            self.assertEqual(row["archaic_haplotype_assignment"], "unavailable_IBDmix_is_unphased")
            self.assertFalse(any(".2013" in r["reference"] for r in rows))
            with (out/"allele_tract_summary.tsv").open() as handle:
                summary = list(csv.DictReader(handle, delimiter="\t"))
            outside = next(r for r in summary if r["reference"] == "Altai" and r["allele"] == "G" and r["tract_call"] == "no_called_tract")
            self.assertEqual((outside["n_tested"], outside["n_called"], outside["AN"], outside["missing_genotypes"]), ("2", "1", "2", "1"))


if __name__ == "__main__":
    unittest.main()
