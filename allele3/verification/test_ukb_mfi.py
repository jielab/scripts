"""Real TSV/BGZF integration tests; no genotype matrix or decoder mocks."""
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"f"))
import ukb_mfi as mfi
from analysis_integrity import validate_manifest
from pca_project import quality_values


def make_mfi(pos=100, a1="A", a2="G", maf="0.1", minor="G", info="0.9"):
    return dict(chrom="1", pos=pos, mfi_line=1, mfi_id=f"1:{pos}_{a1}_{a2}", mfi_rsid=f"rs{pos}",
                allele1=a1, allele2=a2, maf_raw=maf, minor_allele=minor, info_raw=info, input_layout="v3")


def make_pvar(pos=100, ref="A", alt="G", index=0):
    return dict(chrom="1", pos=pos, pvar_index=index, pvar_id=f"rs{pos}", pvar_ref=ref, pvar_alt=alt)


def read_gzip_table(path):
    import csv
    with gzip.open(path, "rt") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def files(root, unsorted=False, header="##reference=GRCh37\n"):
    # Exact, swapped, low INFO, missing INFO/MAF, duplicated MFI/PVAR,
    # native multi-ALT, true PVAR-only, true MFI-only, and strand mismatch.
    mrows = [
        "1:100_A_G rs100 100 A G 0.1 G 0.95",
        "1:110_A_C rs110 110 A C 0.2 A 0.9",
        "1:120_C_T rs120 120 C T 0.2 T 0.79",
        "1:130_A_G rs130 130 A G 0.00001 G NA",
        "1:140_A_G rs140 140 A G NA G 0.99",
        "1:150_A_G rs150 150 A G 0.1 G 0.99",
        "1:150_A_C rs150b 150 A C 0.1 C 0.99",
        "1:160_A_G rs160 160 A G 0.1 G 0.99",
        "1:170_A_G rs170 170 A G 0.1 G 0.99",
        "1:190_A_G rs190 190 A G 0.1 G 0.99",
        "1:200_A_G rs200 200 A G 0.1 G 0.99",
        "1:210_A_AC rs210 210 A AC 0.1 AC 0.99",
    ]
    prows = [
        "1 100 rs100 A G", "1 110 rs110 C A", "1 120 rs120 C T", "1 130 rs130 A G",
        "1 140 rs140 A G", "1 150 rs150 A G", "1 160 rs160 A G", "1 160 rs160b A C",
        "1 170 rs170 A G,C", "1 180 rs180 A G", "1 200 rs200 T C", "1 210 rs210 A AC",
    ]
    if unsorted:
        mrows = list(reversed(mrows)); prows = list(reversed(prows))
    (root/"input.mfi").write_text("\n".join(mrows)+"\n")
    (root/"input.pvar").write_text(header+"#CHROM POS ID REF ALT\n"+"\n".join(prows)+"\n")
    return ["--mfi", str(root/"input.mfi"), "--pvar", str(root/"input.pvar"), "--chromosomes", "1",
            "--genome-build", "GRCh37", "--out", str(root/"out")]


class MfiAdapter(unittest.TestCase):
    def test_actual_official_v3_layout_and_explicit_legacy(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"official_prefix.txt"
            p.write_text("1:10505_A_T\trs548419688\t10505\tA\tT\t1.7934e-05\tT\t0.230206\n")
            row = next(mfi.mfi_rows(p, 1, "v3"))
            self.assertEqual(row["minor_allele"], "T")
            self.assertEqual(float(row["maf_raw"]), 1.7934e-5)
            with self.assertRaisesRegex(ValueError, "invalid MAF"):
                list(mfi.mfi_rows(p, 1, "legacy_doc"))
            p.write_text("1:10505_A_T rs548419688 10505 A T T 1.7934e-05 0.230206\n")
            self.assertEqual(next(mfi.mfi_rows(p, 1, "legacy_doc"))["minor_allele"], "T")
            with self.assertRaisesRegex(ValueError, "invalid MAF"):
                list(mfi.mfi_rows(p, 1, "v3"))

    def test_exact_swap_and_no_complement_guess(self):
        exact = mfi.compare([make_mfi()], [make_pvar()], .8, 0)
        self.assertEqual(exact["match_status"], "matched_exact")
        swapped = mfi.compare([make_mfi()], [make_pvar(ref="G", alt="A")], .8, 0)
        self.assertEqual(swapped["match_status"], "matched_ref_alt_swapped")
        self.assertEqual(swapped["quality_pass"], 1)
        rejected = mfi.compare([make_mfi()], [make_pvar(ref="T", alt="C")], .8, 0)
        self.assertEqual(rejected["match_status"], "allele_pair_mismatch_no_strand_guess")
        self.assertEqual(rejected["quality_pass"], 0)

    def test_invalid_quality_never_becomes_one(self):
        for changes, reason in ((dict(info="NA"), "info_missing"), (dict(info="-1"), "info_out_of_range"),
                                (dict(info="1.01"), "info_out_of_range"), (dict(info="bad"), "info_invalid"),
                                (dict(maf="NA"), "maf_missing"), (dict(maf="0.6"), "maf_out_of_range"),
                                (dict(minor="C"), "minor_allele_mismatch"), (dict(info="0.79"), "info_below_threshold")):
            result = mfi.compare([make_mfi(**changes)], [make_pvar()], .8, 0)
            self.assertEqual(result["quality_pass"], 0)
            self.assertIn(reason, result["quality_status"])
        low_maf = mfi.compare([make_mfi(maf="0.00001")], [make_pvar()], .8, .001)
        self.assertEqual(low_maf["quality_status"], "maf_below_threshold")

    def test_duplicate_positions_and_native_multi_alt_rejected(self):
        for ms, ps, reason in (([make_mfi(), make_mfi(a2="C", minor="C")], [make_pvar()], "duplicate_mfi_position"),
                               ([make_mfi()], [make_pvar(), make_pvar(alt="C")], "duplicate_pvar_position"),
                               ([make_mfi()], [make_pvar(alt="G,C")], "pvar_not_unique_biallelic_snv")):
            result = mfi.compare(ms, ps, .8, 0)
            self.assertEqual(result["match_status"], reason)
            self.assertEqual(result["quality_pass"], 0)

    def test_coordinate_and_column_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"x"
            for content in ("1:100_A_G rs1 101 A G 0.1 G 0.9\n", "2:100_A_G rs1 100 A G 0.1 G 0.9\n",
                            "X:100_A_G rs1 100 A G 0.1 G 0.9\n",
                            "1:100_A_G rs1 100 A G 0.1 0.9\n"):
                p.write_text(content)
                with self.assertRaises(ValueError): list(mfi.mfi_rows(p, 1, "v3"))

    def test_build_conflict_and_declaration_only(self):
        for line in ("##reference=GRCh38", "##assembly=hg38", "##contig=<ID=1,assembly=GRCh38.p14>",
                     "##genome_build=38", "##contig=<ID=1,assembly=38>"):
            with self.assertRaisesRegex(ValueError, "conflicts"):
                mfi.check_build_header(line, "GRCh37")
        self.assertEqual(mfi.check_build_header("##reference=file:///GRCH37.fasta.gz", "GRCh37"), ["GRCh37"])
        with tempfile.TemporaryDirectory() as tmp:
            args = files(Path(tmp), header="")
            result = mfi.run(mfi.parser().parse_args(args))
            self.assertEqual(result["pvar_build_audits"][0]["verification"], "declaration_only_header_uninformative")

    def test_real_projection_info_contract_and_full_audits(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = files(root)
            # Genotype-like files deliberately contain poison; converter cannot
            # accidentally decode or derive AC from them.
            (root/"input.pgen").write_bytes(b"not a genotype file")
            (root/"input.acount").write_text("not observed WGS allele counts\n")
            result = mfi.run(mfi.parser().parse_args(args))
            validate_manifest(root/"out", require_evidence=True)
            self.assertEqual(result["counts"]["quality_pass"], 2)
            self.assertEqual(result["counts"]["mfi_rows"], 12)
            self.assertEqual(result["counts"]["pvar_rows"], 12)
            info = quality_values(root/"out/target_info.tsv.gz", {("1", x) for x in range(100, 211, 10)})
            self.assertEqual(info, {("1", 100, "A", "G"): .95, ("1", 110, "C", "A"): .9})
            audit = read_gzip_table(root/"out/mfi_audit.tsv.gz")
            mask = read_gzip_table(root/"out/pvar_quality_mask.tsv.gz")
            self.assertEqual(len(audit), 12); self.assertEqual(len(mask), 12)
            self.assertEqual(next(r for r in audit if r["pos"] == "190")["match_status"], "outside_pvar_subset")
            self.assertEqual(next(r for r in mask if r["pos"] == "180")["match_status"], "mfi_missing")
            self.assertEqual([int(r["pvar_index"]) for r in mask], list(range(12)))
            self.assertIn("unknown_from_MFI", result["count_evidence"])

    def test_sorted_merge_rejects_disorder_and_external_sort_preserves_indexes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = files(root, unsorted=True)
            with self.assertRaisesRegex(ValueError, "not position sorted"):
                mfi.run(mfi.parser().parse_args(args))
            result = mfi.run(mfi.parser().parse_args([*args, "--sort-inputs", "--sort-memory-mb", "16"]))
            self.assertEqual(result["counts"]["quality_pass"], 2)
            mask = read_gzip_table(root/"out/pvar_quality_mask.tsv.gz")
            self.assertEqual(next(r for r in mask if r["pos"] == "100")["pvar_index"], "11")
            self.assertEqual([int(r["pos"]) for r in mask], sorted(int(r["pos"]) for r in mask))

    def test_templates_multiple_chromosomes_and_input_protection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for chrom in (1, 2):
                (root/f"chr{chrom}.mfi").write_text(f"{chrom}:100_A_G rs100 100 A G 0.1 G 0.9\n")
                (root/f"chr{chrom}.pvar").write_text(f"#CHROM POS ID REF ALT\n{chrom} 100 rs100 A G\n")
            args = ["--mfi-template", str(root/"chr{chrom}.mfi"), "--pvar-template", str(root/"chr{chrom}.pvar"),
                    "--chromosomes", "2,1", "--genome-build", "GRCh37", "--out", str(root/"out")]
            result = mfi.run(mfi.parser().parse_args(args))
            self.assertEqual(result["counts"]["quality_pass"], 2)
            self.assertEqual([r["chrom"] for r in read_gzip_table(root/"out/target_info.tsv.gz")], ["1", "2"])
            args[-1] = str(root)
            with self.assertRaisesRegex(ValueError, "overlaps"):
                mfi.run(mfi.parser().parse_args(args))

    def test_real_bgzf_and_tabix_lookup(self):
        try: import pysam
        except ImportError: self.skipTest("pysam is required for BGZF/tabix integration")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = files(root)
            result = mfi.run(mfi.parser().parse_args([*args, "--tabix"]))
            self.assertEqual(result["tabix_index"], "created")
            with pysam.TabixFile(str(root/"out/target_info.tsv.gz")) as table:
                self.assertEqual(list(table.fetch("1", 99, 100)), ["1\t100\tA\tG\t0.95"])
            self.assertEqual(quality_values(root/"out/target_info.tsv.gz", {("1", 110)}), {("1", 110, "C", "A"): .9})

    def test_cli_main_subprocess_and_failure_has_no_usable_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = files(root)
            run = subprocess.run([sys.executable, str(ROOT/"f/ukb_mfi.py"), *args], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
            self.assertEqual(json.loads(run.stdout.splitlines()[-1])["counts"]["quality_pass"], 2)
            self.assertIn("[build-check]", run.stdout)
            dispatched = subprocess.run([sys.executable, str(ROOT/"f/research.py"), "mfi", *args], capture_output=True, text=True)
            self.assertEqual(dispatched.returncode, 0, dispatched.stdout+dispatched.stderr)
            self.assertEqual(json.loads(dispatched.stdout.splitlines()[-1])["counts"]["quality_pass"], 2)
            files(root, header="##reference=GRCh38\n")
            run = subprocess.run([sys.executable, str(ROOT/"f/ukb_mfi.py"), *args], capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("conflicts", run.stderr)
            with self.assertRaisesRegex(ValueError, "Nonterminal/unusable"):
                validate_manifest(root/"out", require_evidence=True)

    def test_no_passing_markers_is_explicit_and_never_filled(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); args = files(root)
            result = mfi.run(mfi.parser().parse_args([*args, "--min-info", "1"]))
            self.assertEqual(result["status"], "no_usable_exact_high_quality_contrast")
            self.assertEqual(result["counts"]["quality_pass"], 0)
            self.assertEqual(quality_values(root/"out/target_info.tsv.gz", {("1", 100)}), {})
            self.assertEqual(len(read_gzip_table(root/"out/pvar_quality_mask.tsv.gz")), 12)


if __name__ == "__main__":
    unittest.main()
