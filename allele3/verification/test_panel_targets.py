"""Candidate counts from real native/split PGEN fixtures, including QC absences."""
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"f"))
import panel_targets as panel
import ukb_mfi
from analysis_integrity import validate_manifest


def rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def fixture(root):
    import pgenlib
    pvar = [
        "1 100 native A C,G,T", "1 200 bial A G", "1 300 splitC A C", "1 300 splitG A G",
        "1 400 mixed_native A C,G", "1 400 mixed_split A T", "1 500 indel A AT", "1 500 companion A G",
        "1 600 swapped G A", "1 700 missing A G", "1 800 zeros A G",
    ]
    (root/"chr1.pvar").write_text("##reference=GRCh37\n#CHROM POS ID REF ALT\n"+"\n".join(pvar)+"\n")
    (root/"chr1.psam").write_text("#IID\nA\nB\nC\nD\n")
    writer = pgenlib.PgenWriter(str(root/"chr1.pgen").encode(), 4, variant_ct=len(pvar), allele_ct_limit=4)
    writer.append_alleles(np.array([1, 2, -9, -9, 2, 2, 0, 3], dtype=np.int32), allele_ct=4)
    for counts in ([0, 1, 2, -9], [1, 0, 0, 0], [0, 1, 0, 0]):
        writer.append_biallelic(np.array(counts, dtype=np.int8))
    writer.append_alleles(np.array([1, 2, 0, 0, 1, 1, 2, 2], dtype=np.int32), allele_ct=3)
    for counts in ([0, 1, 0, 0], [0, 0, 0, 0], [1, 1, 0, 0], [0, 1, 2, -9], [-9]*4, [0]*4):
        writer.append_biallelic(np.array(counts, dtype=np.int8))
    writer.close()
    targets = [(100, "C,G,T"), (200, "G,C"), (300, "C,G"), (400, "C,G,T"),
               (500, "G,C"), (600, "G,C"), (700, "G,C"), (800, "G,C"), (900, "G,C")]
    (root/"targets.tsv").write_text("site_id\tchrom\tpos\tref\talts\n"+"".join(f"s{pos}\t1\t{pos}\tA\t{alts}\n" for pos, alts in targets))
    return ["--targets", str(root/"targets.tsv"), "--target-build", "37", "--panel-build", "GRCh37",
            "--pgen-template", str(root/"chr{chrom}"), "--evidence-kind", "direct_gt", "--out", str(root/"out")]


def result_lookup(root):
    return {(r["site_id"], r["allele"], r["population"]): r for r in rows(root/"out/allele_counts.tsv.gz")}


def quality_file(root):
    p = root/"quality.tsv"
    p.write_text("chrom\tpos\tref\talt\tinfo\n"+"".join(f"1\t{pos}\t{ref}\t{alt}\t{value}\n" for pos, ref, alt, value in (
        (100, "A", "C", .95), (100, "A", "G", .2), (100, "A", "T", .95),
        (200, "A", "G", .95), (300, "A", "C", .95), (300, "A", "G", .2),
        (400, "A", "C", .95), (400, "A", "G", .95), (500, "A", "G", .95),
        (600, "G", "A", .95), (700, "A", "G", .95), (800, "A", "G", .95))))
    return ["--quality-info", str(p), "--quality-build", "37"]


class PanelTargets(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try: import pgenlib
        except ImportError: raise unittest.SkipTest("pgenlib required for actual PGEN integration tests")

    def test_native_four_bases_use_real_psam_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root)
            result = panel.run(panel.parser().parse_args([*args, "--chunk-size", "2"]))
            validate_manifest(root/"out", require_evidence=True)
            data = result_lookup(root)
            for base, ac in (("A", 1), ("C", 1), ("G", 3), ("T", 1)):
                row = data["s100", base, "ALL"]
                self.assertEqual(int(row["AC"]), ac); self.assertEqual(int(row["AN"]), 6)
                self.assertEqual(row["counting_scope"], "native_joint_nucleotide_count")
                self.assertEqual(float(row["call_rate"]), .75)
                self.assertEqual(row["n_panel_samples"], "4")
            self.assertEqual(result["chromosomes"][0]["n_pvar_rows"], 11)
            self.assertEqual(result["chromosomes"][0]["allele_offset_bytes"], 12*8)

    def test_split_and_missing_encoded_are_not_ref_or_unencoded_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); panel.run(panel.parser().parse_args(fixture(root)))
            data = result_lookup(root)
            self.assertEqual(data["s200", "A", "ALL"]["measurement_status"], "reference_requires_native_joint_genotype")
            self.assertEqual(data["s200", "A", "ALL"]["AN"], "NA")
            self.assertEqual(data["s200", "C", "ALL"]["measurement_status"], "not_represented_in_panel")
            self.assertEqual(data["s200", "G", "ALL"]["AC"], "3")
            self.assertEqual(data["s300", "A", "ALL"]["AN"], "NA")
            self.assertEqual(data["s300", "C", "ALL"]["AC"], "1")
            self.assertEqual(data["s300", "G", "ALL"]["AC"], "1")
            self.assertEqual(data["s700", "G", "ALL"]["AN"], "0")
            self.assertEqual(data["s700", "G", "ALL"]["AF"], "NA")
            self.assertEqual(data["s700", "C", "ALL"]["AN"], "NA")
            self.assertEqual(data["s800", "G", "ALL"]["AN"], "8")
            self.assertEqual(data["s800", "G", "ALL"]["AF"], "0.0")
            self.assertEqual(data["s900", "G", "ALL"]["measurement_status"], "not_represented_in_panel")

    def test_mixed_native_and_non_snv_companions_cannot_disappear(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root)
            panel.run(panel.parser().parse_args([*args, *quality_file(root)]))
            data = result_lookup(root)
            self.assertEqual(data["s400", "C", "ALL"]["measurement_status"], "mixed_native_and_split_representation")
            self.assertEqual(data["s500", "G", "ALL"]["measurement_status"], "non_snv_or_invalid_companion_at_locus")
            self.assertEqual(data["s300", "A", "ALL"]["measurement_status"], "reference_requires_native_joint_genotype")
            self.assertEqual(data["s300", "G", "ALL"]["measurement_status"], "not_in_high_quality")
            self.assertEqual(data["s100", "A", "ALL"]["measurement_status"], "not_in_high_quality")
            self.assertEqual(data["s100", "C", "ALL"]["AC"], "1")

    def test_target_ref_explicitly_encoded_as_alt_can_be_counted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); panel.run(panel.parser().parse_args(fixture(root)))
            row = result_lookup(root)["s600", "A", "ALL"]
            self.assertEqual(row["AC"], "3")
            self.assertEqual(row["allele_idx"], "1")
            self.assertEqual(row["counting_scope"], "panel_encoded_ALT_count")
            self.assertEqual(row["reference_label_relation"], "target_REF_explicitly_in_panel_alleles")

    def test_population_subset_joins_by_iid_and_uses_selected_callability(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root)
            p = root/"samples.tsv"; p.write_text("sample_id\tpopulation\nD\tEUR\nB\tAFR\nA\tEUR\n")
            panel.run(panel.parser().parse_args([*args, "--samples", str(p)]))
            data = result_lookup(root)
            self.assertEqual(data["s100", "G", "ALL"]["AC"], "1")
            self.assertEqual(data["s100", "G", "ALL"]["AN"], "4")
            self.assertEqual(data["s100", "G", "ALL"]["n_selected_samples"], "3")
            self.assertEqual(data["s100", "G", "EUR"]["AN"], "4")
            self.assertEqual(data["s100", "G", "AFR"]["AN"], "0")
            p.write_text("sample_id\nmissing\n")
            with self.assertRaisesRegex(ValueError, "absent from PSAM"):
                panel.run(panel.parser().parse_args([*args, "--samples", str(p)]))

    def test_imputed_requires_quality_and_labels_hardcall_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root); args[args.index("direct_gt")] = "imputed_hardcall"
            with self.assertRaisesRegex(ValueError, "requires --quality-info"):
                panel.run(panel.parser().parse_args(args))
            quality = quality_file(root)
            result = panel.run(panel.parser().parse_args([*args, *quality]))
            self.assertEqual(result["evidence_kind"], "imputed_hardcall")
            self.assertIn("no posterior", result["genotype_semantics"])
            self.assertEqual(result_lookup(root)["s200", "G", "ALL"]["evidence_kind"], "imputed_hardcall")
            with self.assertRaisesRegex(ValueError, "requires explicit --quality-build"):
                panel.run(panel.parser().parse_args([*args, *quality[:2]]))

    def test_mfi_output_contract_binds_original_pvar_and_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root); args[args.index("direct_gt")] = "imputed_hardcall"
            mf = root/"mfi.txt"
            mf.write_text("1:200_A_G rs200 200 A G 0.1 G 0.95\n1:600_A_G rs600 600 A G 0.2 G 0.95\n1:800_A_G rs800 800 A G 0.1 G 0.95\n")
            ukb_mfi.run(ukb_mfi.parser().parse_args(["--mfi", str(mf), "--pvar", str(root/"chr1.pvar"),
                "--chromosomes", "1", "--genome-build", "GRCh37", "--out", str(root/"mfi_result")]))
            args += ["--quality-info", str(root/"mfi_result/target_info.tsv.gz")]
            result = panel.run(panel.parser().parse_args(args))
            self.assertEqual(result["quality_contract"]["status"], "verified_mfi_output_and_original_PVAR")
            self.assertEqual(result_lookup(root)["s200", "G", "ALL"]["AC"], "3")
            self.assertEqual(result_lookup(root)["s100", "C", "ALL"]["measurement_status"], "not_in_high_quality")
            with (root/"chr1.pvar").open("a") as handle: handle.write("\n")
            with self.assertRaisesRegex(ValueError, "PVAR changed"):
                panel.run(panel.parser().parse_args(args))

    def test_explicit_pvar_resolves_coexisting_companions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root)
            (root/"chr1.pvar.zst").write_bytes(b"unused deliberate conflict")
            with self.assertRaisesRegex(ValueError, "Exactly one PVAR"):
                panel.run(panel.parser().parse_args(args))
            result = panel.run(panel.parser().parse_args([*args, "--pvar-template", str(root/"chr{chrom}.pvar")]))
            self.assertGreater(result["counts"]["measured_nucleotides"], 0)

    def test_explicit_build_xy_target_cap_and_dimension_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root)
            mismatched = args.copy(); mismatched[mismatched.index("GRCh37")] = "GRCh38"
            with self.assertRaisesRegex(ValueError, "builds must match"):
                panel.run(panel.parser().parse_args(mismatched))
            with self.assertRaisesRegex(ValueError, "cap"):
                panel.run(panel.parser().parse_args([*args, "--max-targets", "1"]))
            with (root/"chr1.psam").open("a") as handle: handle.write("E\n")
            with self.assertRaises((RuntimeError, ValueError)):
                panel.run(panel.parser().parse_args(args))
            (root/"targets.tsv").write_text("site_id chrom pos ref alts\nx X 100 A C,G\n")
            with self.assertRaisesRegex(ValueError, "X/Y"):
                panel.run(panel.parser().parse_args(args))

    def test_explicit_target_row_build_cannot_be_relabelled_by_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root)
            for column in ("build", "genome_build", "assembly", "target_build"):
                (root/"targets.tsv").write_text(f"site_id chrom pos ref alts {column}\ns100 1 100 A C,G,T GRCh38\n")
                with self.assertRaisesRegex(ValueError, "conflicts with --target-build"):
                    panel.run(panel.parser().parse_args(args))
            (root/"targets.tsv").write_text("site_id chrom pos ref alts genome_build\ns100 1 100 A C,G,T 37\n")
            result = panel.run(panel.parser().parse_args(args))
            self.assertEqual(result["counts"]["measured_nucleotides"], 4)

    def test_requested_bases_exclude_structural_ref_without_losing_native_encoding(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root)
            (root/"targets.tsv").write_text("site_id chrom pos ref alts build requested_bases\ns100 1 100 A C,G,T GRCh37 C,G\n")
            result = panel.run(panel.parser().parse_args(args))
            self.assertEqual(result["counts"]["requested_nucleotides"], 2)
            self.assertEqual(result["counts"]["measured_nucleotides"], 2)
            data = result_lookup(root)
            self.assertEqual({key[1] for key in data}, {"C", "G"})
            self.assertEqual(data["s100", "G", "ALL"]["AC"], "3")
            self.assertEqual(data["s100", "G", "ALL"]["counting_scope"], "native_joint_nucleotide_count")
            audit = rows(root/"out/target_availability.tsv")[0]
            self.assertEqual(audit["structural_ref_requested"], "0")
            self.assertEqual(audit["n_requested_bases"], "2")
            self.assertEqual(audit["requested_bases"], "C,G")
            self.assertEqual(audit["panel_representation"], "native_multiallelic_row")
            for requested in ("A,A", "A,X", "NA", "C,,G"):
                (root/"targets.tsv").write_text(f"site_id chrom pos ref alts requested_bases\ns100 1 100 A C,G,T {requested}\n")
                with self.assertRaisesRegex(ValueError, "requested_bases"):
                    panel.run(panel.parser().parse_args(args))

    def test_direct_cli_real_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = fixture(root)
            run = subprocess.run([sys.executable, str(ROOT/"f/panel_targets.py"), *args], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
            self.assertEqual(json.loads(run.stdout.splitlines()[-1])["counts"]["target_sites"], 9)
            self.assertIn("actual_PSAM_N=4 selected_N=4", run.stdout)
            self.assertIn("[build-check]", run.stdout)
            validate_manifest(root/"out", require_evidence=True)
            dispatched = subprocess.run([sys.executable, str(ROOT/"f/research.py"), "panel", *args], capture_output=True, text=True)
            self.assertEqual(dispatched.returncode, 0, dispatched.stdout+dispatched.stderr)
            self.assertEqual(json.loads(dispatched.stdout.splitlines()[-1])["counts"]["target_sites"], 9)
            validate_manifest(root/"out", require_evidence=True)


if __name__ == "__main__":
    unittest.main()
