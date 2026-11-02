import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
from allofus_sites import open_text, parser, run


class TestAllofusSites(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="allele3-allofus-")
        self.root = Path(self.temp.name)
        self.attempt = 0

    def tearDown(self):
        self.temp.cleanup()

    def prepare(self, text, extra=(), mode="targeted", counts=False,
                target_text="chrom\tpos\tref\talts\tfamily\tgroup\n1\t10\tA\tC,G\tERV1\tcase\n",
                metrics=True, gzip_source=True):
        self.attempt += 1
        source = self.root / f"source{self.attempt}.tsv.gz"
        if gzip_source:
            with gzip.open(source, "wt") as handle:
                handle.write(text)
        else:
            source = source.with_suffix("")
            source.write_text(text)
        out = self.root / f"result{self.attempt}"
        argv = ["--source", str(source), "--build", "38", "--out", str(out), "--mode", mode]
        if mode != "preview":
            argv += ["--coordinate-system", "1-based", "--progress-every", "0"]
            if metrics:
                argv += ["--population", "ALL", "--af-column", "AF"]
                if counts:
                    argv += ["--ac-column", "AC", "--an-column", "AN"]
        if mode == "targeted":
            target = self.root / f"targets{self.attempt}.tsv"
            target.write_text(target_text)
            argv += ["--targets", str(target), "--targets-build", "GRCh38"]
        return parser().parse_args(argv + list(extra)), out

    def table(self, path):
        with open_text(path) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_preview_default_is_bounded_and_does_not_guess_schema(self):
        args, out = self.prepare("unknown\tcolumns\n" + "1\t10\tA\tC\t0.05\n" * 100,
                                 mode="preview", extra=("--preview-chars", "32", "--preview-lines", "2"))
        self.assertEqual(parser().get_default("mode"), "preview")
        result = run(args)
        preview = json.loads((out / "preview.json").read_text())
        self.assertTrue(preview["prefix_truncated"])
        self.assertLessEqual(len(preview["lines"]), 2)
        self.assertFalse((out / "allofus_alleles.tsv.gz").exists())
        self.assertEqual(result["scan_policy"], "bounded prefix")

    def test_headerless_needs_explicit_column_order(self):
        text = "1\t10\tA\tC\t0.05\n1\t10\tA\tG\t0.02\n"
        args, out = self.prepare(text)
        with self.assertRaisesRegex(ValueError, "No usable source header"):
            run(args)
        self.assertEqual(json.loads((out / "allofus_manifest.json").read_text())["status"], "failed")
        args, out = self.prepare(text, extra=("--columns", "chrom,pos,ref,alt,AF"))
        run(args)
        row = self.table(out / "target_sites.tsv")[0]
        self.assertEqual(row["frequency_supported_multi_alt"], "1")
        self.assertEqual(row["count_supported_multi_alt"], "NA")

    def test_annovar_header_exact_mapping_and_zero_based_half_open_coordinates(self):
        text = "#Chr\tStart\tEnd\tRef\tAlt\tFreq\nchr1\t9\t10\tA\tC\t0.05\n"
        args, out = self.prepare(text, metrics=False, gzip_source=False, extra=(
            "--chrom-column", "Chr", "--pos-column", "Start", "--end-column", "End",
            "--ref-column", "Ref", "--alt-column", "Alt", "--population", "ALL",
            "--af-column", "Freq", "--coordinate-system", "0-based"))
        run(args)
        row = self.table(out / "target_alleles.tsv")[0]
        self.assertEqual(row["pos"], "10")
        self.assertEqual(row["match_status"], "matched_frequency_only")

    def test_gzip_af_only_duplicates_missing_metrics_and_unreported_alt(self):
        text = ("chrom\tpos\tref\talt\tAF\ttranscript\n"
                "1\t10\tA\tC\t0.02\tTX1\n"
                "1\t10\tA\tG\t.\tTX3\n"
                "1\t10\tA\tC\t0.02\tTX2\n")
        targets = "chrom\tpos\tref\talts\n1\t10\tA\tC,G,T\n1\t20\tG\tA,T\n"
        args, out = self.prepare(text, target_text=targets)
        result = run(args)
        self.assertEqual(result["analysis"]["counts"]["identical_duplicate_rows"], 1)
        rows = self.table(out / "target_alleles.tsv")
        at10 = {r["alt"]: r for r in rows if r["pos"] == "10"}
        self.assertEqual(at10["C"]["match_status"], "matched_frequency_only")
        self.assertEqual(at10["C"]["AC"], "NA")
        self.assertEqual(at10["C"]["AN"], "NA")
        self.assertEqual(at10["G"]["match_status"], "matched_metrics_missing")
        self.assertEqual(at10["G"]["AF"], "NA")
        self.assertEqual(at10["T"]["match_status"], "not_in_summary")
        self.assertEqual(at10["T"]["AC"], "NA")
        self.assertTrue(all(r["match_status"] == "not_in_summary" for r in rows if r["pos"] == "20"))

    def test_true_counts_multi_alt_do_not_assume_reference_presence(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\n"
                "1\t10\tA\tC\t5\t100\t0.05\n1\t10\tA\tG\t7\t100\t0.07\n")
        args, out = self.prepare(text, counts=True)
        run(args)
        row = self.table(out / "target_sites.tsv")[0]
        self.assertEqual(row["count_supported_multi_alt"], "1")
        self.assertEqual(row["target_all_alts_count_supported"], "1")
        self.assertEqual(row["ref_AC"], "NA")
        self.assertEqual(row["n_supported_bases_lower_bound"], "2")
        self.assertEqual(row["at_least_three_bases_demonstrated"], "0")

    def test_explicit_reference_zero_is_not_three_bases(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\tREF_AC\n"
                "1\t10\tA\tC\t60\t100\t0.6\t0\n1\t10\tA\tG\t40\t100\t0.4\t0\n")
        args, out = self.prepare(text, counts=True, extra=("--ref-ac-column", "REF_AC"))
        run(args)
        row = self.table(out / "target_sites.tsv")[0]
        self.assertEqual(row["count_supported_multi_alt"], "1")
        self.assertEqual(row["ref_AC"], "0")
        self.assertEqual(row["at_least_three_bases_demonstrated"], "0")

    def test_three_measured_bases_can_be_established(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\tREF_AC\n"
                "1\t10\tA\tC\t5\t100\t0.05\t90\n1\t10\tA\tG\t5\t100\t0.05\t90\n")
        args, out = self.prepare(text, counts=True, extra=("--ref-ac-column", "REF_AC"))
        run(args)
        self.assertEqual(self.table(out / "target_sites.tsv")[0]["at_least_three_bases_demonstrated"], "1")

    def test_additional_alt_is_retained_even_if_not_requested(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\n"
                "1\t10\tA\tC\t5\t100\t0.05\n1\t10\tA\tG\t7\t100\t0.07\n1\t10\tA\tT\t6\t100\t0.06\n")
        args, out = self.prepare(text, counts=True)
        run(args)
        row = self.table(out / "target_sites.tsv")[0]
        self.assertEqual(row["additional_reported_alts"], "T")
        self.assertEqual(row["at_least_three_bases_demonstrated"], "1")
        self.assertEqual(len(self.table(out / "allofus_alleles.tsv.gz")), 3)

    def test_missing_an_does_not_establish_a_count_threshold(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\n"
                "1\t10\tA\tC\t5\t.\t0.05\n1\t10\tA\tG\t7\t.\t0.07\n")
        args, out = self.prepare(text, counts=True)
        run(args)
        row = self.table(out / "target_sites.tsv")[0]
        self.assertEqual(row["count_supported_multi_alt"], "NA")
        self.assertEqual(row["frequency_supported_multi_alt"], "1")

    def test_missing_an_counts_are_not_pooled_with_another_alts_denominator(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\tREF_AC\n"
                "1\t10\tA\tC\t700\t.\t0.7\t200\n1\t10\tA\tG\t50\t100\t0.5\t.\n")
        args, out = self.prepare(text, counts=True, extra=("--ref-ac-column", "REF_AC"))
        run(args)
        row = self.table(out / "target_sites.tsv")[0]
        self.assertEqual(row["ref_AC"], "NA")
        self.assertEqual(row["count_supported_multi_alt"], "NA")
        self.assertEqual(row["n_supported_bases_lower_bound"], "1")
        self.assertEqual(row["reported_AF_sum_exceeds_one"], "1")

    def test_af_only_nonjoint_frequencies_are_flagged_without_fabricating_counts(self):
        text = "chrom\tpos\tref\talt\tAF\n1\t10\tA\tC\t0.7\n1\t10\tA\tG\t0.6\n"
        args, out = self.prepare(text)
        run(args)
        row = self.table(out / "target_sites.tsv")[0]
        self.assertEqual(row["reported_AF_sum_exceeds_one"], "1")
        self.assertEqual(row["count_supported_multi_alt"], "NA")
        self.assertEqual(row["n_supported_bases_lower_bound"], "NA")

    def test_zero_callable_alleles_is_not_frequency_or_count_absence_evidence(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\n"
                "1\t10\tA\tC\t0\t0\t.\n1\t10\tA\tG\t0\t0\t.\n")
        args, out = self.prepare(text, counts=True)
        run(args)
        row = self.table(out / "target_sites.tsv")[0]
        self.assertEqual(row["count_supported_multi_alt"], "NA")
        self.assertEqual(row["frequency_supported_multi_alt"], "NA")
        self.assertEqual(self.table(out / "target_alleles.tsv")[0]["match_status"], "matched_no_callable_alleles")

    def test_af_can_be_derived_only_from_real_counts(self):
        args, out = self.prepare("chrom\tpos\tref\talt\tAC\tAN\tAF\n1\t10\tA\tC\t5\t100\t.\n", counts=True)
        run(args)
        row = self.table(out / "allofus_alleles.tsv.gz")[0]
        self.assertEqual(row["AF"], "0.05")
        self.assertEqual(row["af_source"], "derived_from_reported_AC_AN")
        self.assertEqual(row["ref_AC"], "NA")

    def test_invalid_counts_duplicate_conflicts_and_reference_conflicts_fail_closed(self):
        cases = [
            ("1\t10\tA\tC\t5\t100\t.\n1\t10\tA\tG\t5\t200\t.\n", "Inconsistent AN"),
            ("1\t10\tA\tC\t70\t100\t.\n1\t10\tA\tG\t50\t100\t.\n", "exceeds AN"),
            ("1\t10\tA\tC\t5\t100\t.\n1\t10\tA\tC\t6\t100\t.\n", "Conflicting duplicate"),
            ("1\t10\tA\tC\t5\t100\t.\n1\t10\tG\tA\t6\t100\t.\n", "Conflicting REF"),
            ("1\t10\tA\tA\t5\t100\t.\n", "REF/ALT identity"),
            ("1\t10\tA\tC\t5.5\t100\t.\n", "nonnegative integer"),
            ("1\t10\tA\tC\t5\t100\t0.8\n", "AF disagrees"),
            ("1\t10\tA\tC\t0\t0\t0.01\n", "Positive AF with AN=0"),
        ]
        for body, message in cases:
            with self.subTest(message=message):
                args, out = self.prepare("chrom\tpos\tref\talt\tAC\tAN\tAF\n" + body, counts=True)
                with self.assertRaisesRegex(ValueError, message):
                    run(args)
                self.assertEqual(json.loads((out / "allofus_manifest.json").read_text())["status"], "failed")

    def test_population_columns_never_pool_alt_counts(self):
        text = ("chrom\tpos\tref\talt\tAC_EUR\tAN_EUR\tAC_AFR\tAN_AFR\n"
                "1\t10\tA\tC\t10\t20\t0\t20\n1\t10\tA\tG\t0\t20\t10\t20\n")
        mapping = json.dumps({"EUR": {"ac": "AC_EUR", "an": "AN_EUR"}, "AFR": {"ac": "AC_AFR", "an": "AN_AFR"}})
        args, out = self.prepare(text, metrics=False, extra=("--population-map", mapping))
        run(args)
        rows = self.table(out / "target_sites.tsv")
        self.assertEqual({r["population"] for r in rows}, {"AFR", "EUR"})
        self.assertTrue(all(r["count_supported_multi_alt"] == "0" for r in rows))
        self.assertTrue(all(r["n_alts_ac_ge_min"] == "1" for r in rows))

    def test_reference_mismatch_is_not_a_target_allele_match(self):
        args, out = self.prepare("chrom\tpos\tref\talt\tAF\n1\t10\tG\tC\t0.1\n")
        run(args)
        rows = self.table(out / "target_alleles.tsv")
        self.assertTrue(all(r["match_status"] == "reference_mismatch" and r["AC"] == "NA" for r in rows))
        self.assertTrue(all(r["source_ref"] == "G" and r["ref"] == "A" for r in rows))

    def test_target_and_source_builds_must_agree(self):
        args, _ = self.prepare("chrom\tpos\tref\talt\tAF\n1\t10\tA\tC\t0.1\n",
                               extra=("--targets-build", "37"))
        args.source = str(self.root / "nonexistent-large-source.gz")
        with self.assertRaisesRegex(ValueError, "Target build differs"):
            run(args)

    def test_row_build_is_checked_even_outside_targets(self):
        args, _ = self.prepare("chrom\tpos\tref\talt\tAF\tassembly\n1\t20\tA\tC\t0.1\t37\n",
                               extra=("--build-column", "assembly"))
        with self.assertRaisesRegex(ValueError, "row build conflicts"):
            run(args)

    def test_normalize_streams_multialt_vectors_and_deduplicates(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\n"
                "1\t10\tA\tC,G\t5,7\t100\t0.05,0.07\n"
                "1\t10\tA\tC\t5\t100\t0.05\n1\t20\tC\tA\t1\t100\t0.01\n")
        args, out = self.prepare(text, mode="normalize", counts=True)
        result = run(args)
        self.assertEqual(len(self.table(out / "allofus_sites.tsv.gz")), 2)
        self.assertEqual(len(self.table(out / "allofus_alleles.tsv.gz")), 3)
        self.assertEqual(result["analysis"]["counts"]["identical_duplicate_rows"], 1)
        self.assertFalse((out / "target_sites.tsv").exists())

    def test_mixed_snv_indel_record_retains_the_actual_snv_vector_entries(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\n"
                "1\t10\tA\tC,AT,G\t5,8,7\t100\t0.05,0.08,0.07\n")
        args, out = self.prepare(text, counts=True)
        result = run(args)
        rows = self.table(out / "target_alleles.tsv")
        self.assertEqual({r["alt"]: r["AC"] for r in rows}, {"C": "5", "G": "7"})
        self.assertEqual(result["analysis"]["counts"]["excluded_non_snv_alts"], 1)
        self.assertEqual(self.table(out / "target_sites.tsv")[0]["ref_AC"], "NA")

    def test_mixed_record_checks_all_original_alt_counts_before_excluding_indels(self):
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\n"
                "1\t10\tA\tC,AT,G\t5,95,7\t100\t.\n")
        args, out = self.prepare(text, counts=True)
        with self.assertRaisesRegex(ValueError, "original ALT counts.*exceeds AN"):
            run(args)
        self.assertEqual(json.loads((out / "allofus_manifest.json").read_text())["status"], "failed")

    def test_mixed_record_ref_count_added_once_and_missing_ac_not_invented(self):
        for counts, ref_count, expected in (("5,8,7", "80", "80"), ("5,.,7", ".", "NA")):
            with self.subTest(counts=counts):
                text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\tREF_AC\n"
                        f"1\t10\tA\tC,AT,G\t{counts}\t100\t.\t{ref_count}\n")
                args, out = self.prepare(text, counts=True, extra=("--ref-ac-column", "REF_AC"))
                run(args)
                row = self.table(out / "target_sites.tsv")[0]
                self.assertEqual(row["ref_AC"], expected)
                self.assertEqual(row["count_supported_multi_alt"], "1")
        text = ("chrom\tpos\tref\talt\tAC\tAN\tAF\tREF_AC\n"
                "1\t10\tA\tC,AT,G\t5,8,7\t100\t.\t90\n")
        args, _ = self.prepare(text, counts=True, extra=("--ref-ac-column", "REF_AC"))
        with self.assertRaisesRegex(ValueError, "original ALT counts.*exceeds AN"):
            run(args)

    def test_normalize_rejects_unsorted_or_repeated_chromosome_blocks(self):
        for body in ("1\t20\tA\tC\t0.1\n1\t10\tA\tG\t0.1\n",
                     "1\t10\tA\tC\t0.1\n2\t10\tA\tC\t0.1\n1\t20\tA\tG\t0.1\n"):
            with self.subTest(body=body):
                args, _ = self.prepare("chrom\tpos\tref\talt\tAF\n" + body, mode="normalize")
                with self.assertRaisesRegex(ValueError, "normalize requires"):
                    run(args)

    def test_cli_default_preview_runs_without_population_mapping(self):
        source = self.root / "arbitrary-summary.txt"
        source.write_text("not\ta\tverified\tschema\n")
        out = self.root / "cli-preview"
        result = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / "f" / "allofus_sites.py"),
                                 "--source", str(source), "--build", "38", "--out", str(out)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((out / "preview.json").is_file())
        self.assertEqual(json.loads((out / "allofus_manifest.json").read_text())["settings"]["mode"], "preview")

    def test_actual_annovar_preset_keeps_eight_populations_max_indicator_and_blanks(self):
        from allofus_sites import ANNOVAR_COLUMNS
        text = ("#"+"\t".join(ANNOVAR_COLUMNS)+"\n"
                "chr1\t10001\t10001\tT\tA\t0.013301\t0.040816\t0.001887\t0.014706\t0.033333\t0.015748\t0.0\t0.040816\t0.0\n"
                "chr1\t10001\t10001\tT\tC\t0.003448\t0.005906\t0.001887\t0.0\t0.0\t0.005906\t0.0\t0.0\t0.0\n"
                "chr1\t10001\t10001\tT\tG\t0.001478\t0.002953\t0.0\t0.0\t0.0\t0.002953\t\t0.0\t0.0\n"
                "chr1\t10110\t10110\tA\tAAC\t0.000165\t0.000199\t0.000145\t0.000183\t0.0\t0.000199\t0.0\t0.0\t0.0\n"
                "chr1\t10110\t10110\tACCCT\tA\t0.000237\t0.000546\t0.000248\t0.000304\t0.00051\t0.000189\t0.0\t0.000254\t0.000546\n")
        targets="chrom\tpos\tref\talts\n1\t10001\tT\tA,C\n1\t10110\tA\tG\n"
        args,out=self.prepare(text,metrics=False,target_text=targets,extra=("--preset","annovar-allofus-hg38"))
        result=run(args)
        rows=self.table(out/"allofus_alleles.tsv.gz")
        self.assertEqual(len(rows),24)
        self.assertEqual({r["population"] for r in rows},{"ALL","AFR","AMR","EAS","EUR","MID","OTH","SAS"})
        self.assertTrue(all(r["AC"]=="NA" and r["AN"]=="NA" for r in rows))
        mid_g=next(r for r in rows if r["population"]=="MID" and r["alt"]=="G")
        self.assertEqual(mid_g["AF"],"NA")
        self.assertEqual(mid_g["reported_max_population_AF"],"0.002953")
        site=next(r for r in self.table(out/"allofus_sites.tsv.gz") if r["population"]=="ALL")
        self.assertEqual(site["n_alts_positive_AF"],"3")
        self.assertEqual(site["count_supported_multi_alt"],"NA")
        self.assertEqual(result["analysis"]["counts"]["selected_non_snv_rows"],2)

    def test_preset_infers_only_declared_schema_build_not_counts(self):
        from allofus_sites import apply_preset, ANNOVAR_POPULATIONS
        args=parser().parse_args(["--source","no-file.gz","--preset","annovar-allofus-hg38","--out","unused"])
        apply_preset(args)
        self.assertEqual(args.build,"GRCh38")
        self.assertEqual(args.coordinate_system,"1-based")
        self.assertEqual(json.loads(args.population_map),ANNOVAR_POPULATIONS)
        self.assertIsNone(args.ac_column)

    def test_preset_refuses_conflicting_build_schema_or_whitespace(self):
        from allofus_sites import apply_preset
        for extra in (("--build","37"),("--delimiter","whitespace"),("--coordinate-system","0-based"),
                      ("--population","MAX","--af-column","gvs_max_af"),("--chrom-column","wrong")):
            with self.subTest(extra=extra):
                args=parser().parse_args(["--source","none","--preset","annovar-allofus-hg38","--out","unused",*extra])
                with self.assertRaises(ValueError):apply_preset(args)


if __name__ == "__main__":
    unittest.main()
