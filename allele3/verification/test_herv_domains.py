"""Domain partition, original-denominator, and supported LOCO counterexamples."""
import csv
from collections import defaultdict
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
import numpy as np
import herv_domains as hd
import scan


class DomainAssignmentTests(unittest.TestCase):
    def test_bed_boundary_and_duplicate_annotations(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "domains.bed"
            p.write_text("chr1\t10\t20\tGag\t0\t+\nchr1\t10\t20\tGag\t0\t+\n")
            domains, labels = hd.read_domains(p)
            self.assertEqual((len(domains), labels), (1, ["Gag"]))
            sweep = hd.DomainSweep(domains)
            self.assertEqual(hd.assign_domain(sweep.locate("1", 10), {"+"}), "unannotated")
            self.assertEqual(hd.assign_domain(sweep.locate("1", 11), {"+"}), "Gag")
            self.assertEqual(hd.assign_domain(sweep.locate("1", 20), {"+"}), "Gag")
            self.assertEqual(hd.assign_domain(sweep.locate("1", 21), {"+"}), "unannotated")

    def test_conflicting_domains_and_strands_retained_as_residuals(self):
        gag = hd.Domain("1", 0, 100, "Gag", "+")
        pol = hd.Domain("1", 0, 100, "Pol", "+")
        self.assertEqual(hd.assign_domain([gag, pol], {"+"}), "ambiguous")
        self.assertEqual(hd.assign_domain([gag], {"-"}), "strand_mismatch")
        self.assertEqual(hd.assign_domain([gag], {"+", "-"}), "ambiguous")
        self.assertEqual(hd.assign_domain([gag], {"+"}), "Gag")

    def test_no_geometric_name_inference_and_no_reversed_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "bad.bed"
            p.write_text("chr1\t20\t10\tGag\n")
            with self.assertRaisesRegex(ValueError, "invalid domain"):
                hd.read_domains(p)
            p.write_text("chr1\t10\t20\tunannotated\n")
            with self.assertRaisesRegex(ValueError, "reserved label"):
                hd.read_domains(p)


class DomainStatisticsTests(unittest.TestCase):
    def test_expected_weighting_not_unweighted_mean(self):
        data = {("ERVK", "strict", 5, "Gag"): {"B1": np.array([20., 10., 100., 0., 0.])},
                ("ERVK", "strict", 5, "Pol"): {"B2": np.array([40., 90., 900., 0., 0.])}}
        rows, audit, _ = hd.summarize_domains(data, {}, ["Gag", "Pol"], ["ERVK"])
        total = next(r for r in rows if r["domain"] == "total")
        self.assertAlmostEqual(total["oe"], .6)
        self.assertAlmostEqual(audit[0]["expected_weighted_oe"], .6)
        self.assertAlmostEqual(audit[0]["expected_weight_sum"], 1)
        self.assertEqual(audit[0]["observed_delta"], 0)
        self.assertIsNone(total["ci_low"])

    def test_unsupported_part_invalidates_full_scope_weighted_decomposition(self):
        data = {("ERVK", "strict", 5, "Gag"): {"B1": np.array([20., 10., 100., 0., 0.])},
                ("ERVK", "strict", 5, "unannotated"): {"B2": np.array([1., 0., 2., 1., 2.])}}
        rows, audit, _ = hd.summarize_domains(data, {}, ["Gag"], ["ERVK"])
        self.assertTrue(all(r["oe"] is None and r["expected_weight"] is None for r in rows))
        total = next(r for r in rows if r["domain"] == "total")
        self.assertEqual((total["observed_multi"], total["eligible_sites"]), (21, 102))
        self.assertEqual(total["n_without_loco"], 2)
        self.assertEqual(audit[0]["status"], "raw_accounting_only_loco_unavailable")

    def test_outside_ld_is_counted_but_not_used_for_inference(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "counts.tsv.gz"
            with gzip.open(p, "wt") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(hd.FIELDS)
                w.writerow(["1", "B1", "ERVK", "Gag", "all", 5, "AAA", "bin", "an", 10, 2])
                w.writerow(["1", "outside_partition", "ERVK", "Gag", "all", 5, "AAA", "bin", "an", 100, 90])
            background = SimpleNamespace(key=lambda r: (), rate=lambda c, key: .1)
            data, raw, outside = hd.aggregate([p], background)
            rows, _, _ = hd.summarize_domains(data, outside, ["Gag"], ["ERVK"])
            total = next(r for r in rows if r["domain"] == "total")
            self.assertEqual((total["eligible_sites"], total["observed_multi"], total["outside_partition_sites"]), (10, 2, 100))
            self.assertEqual(raw[("1", "ERVK", "all", 5, "all_family_sites", "AAA", "bin", "an")].tolist(), [110, 92])


class DomainQueryTests(unittest.TestCase):
    def test_local_query_counts_single_alt_controls_and_residuals(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fasta = root / "ref.fa"
            fasta.write_text(">chr1\n" + "A" * 100 + "\n")
            Path(str(fasta) + ".fai").write_text("chr1\t100\t6\t100\t101\n")
            ld = root / "ld"
            ld.mkdir()
            (ld / "EUR.38.bed").write_text("chr1\t0\t60\n")
            annotation, vcf, domain = root / "annotation.tsv", root / "x.vcf.gz", root / "domains.bed"
            for p in (annotation, vcf, Path(str(vcf) + ".tbi")):
                p.write_text("")
            domain.write_text("chr1\t10\t30\tGag\nchr1\t50\t70\tPol\nchr1\t50\t60\tEnv\n")
            masks = {}
            for name in ("segdup", "lcr", "low_mappability"):
                p = root / (name + ".bed")
                p.write_text("")
                masks["architecture_" + name + "_bed"] = str(p)
            options = dict(fasta=str(fasta), ld_block_dir=str(ld), architecture_build="38", include_unfiltered=False,
                           exclude_flags=[], min_an=100, min_af=0, **masks)
            args = SimpleNamespace(out=str(root / "result"), fasta=None, families="ERVK", bcftools="mock", domain_bed=str(domain))
            header = hd.hm._parse_header('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            records = [dict(chrom="chr1", pos=pos, ref="A", alts=alts, filter="PASS", info=dict(AC=ac, AN="100"))
                       for pos, alts, ac in ((11, ["C"], "10"), (21, ["C", "G"], "6,6"),
                                             (51, ["C", "G"], "6,6"), (61, ["C"], "1"))]
            # Fragmented aligned blocks exclude 31..50; duplicate same-family
            # annotations must never count the same eligible position twice.
            e = scan.ERV("1", 10, 70, "ERV-test", "LTR", "ERVK", "+", "test", ((10, 30), (50, 70)))
            domains, _ = hd.read_domains(domain)
            with mock.patch.object(hd.hm, "read_indexed_header", return_value=header), mock.patch.object(scan, "parse_joined_rmsk", return_value=[e, e]), mock.patch.object(hd.hm, "query_regions", return_value=iter(records)):
                result = hd.scan_chrom(str(vcf), "1", annotation, options, args, domains)
            rows = list(hd.hm.rows(result))
            selected = [r for r in rows if r["sensitivity"] == "strict" and r["threshold"] == "5"]
            self.assertEqual(sum(int(r["n_sites"]) for r in selected), 4)
            self.assertEqual(sum(int(r["n_multi"]) for r in selected), 2)
            self.assertEqual(sum(int(r["n_sites"]) for r in selected if r["domain"] == "ambiguous"), 1)
            self.assertEqual(sum(int(r["n_sites"]) for r in selected if r["block_id"] == "outside_partition"), 1)
            self.assertEqual((root / "result/cache/chr1.regions.bed").read_text(), "chr1\t10\t30\nchr1\t50\t70\n")


class DomainCounterIntegrationTests(unittest.TestCase):
    def fixture(self, root, mismatch=False):
        modern = root / "modern"
        (modern / "summary").mkdir(parents=True)
        bed = root / "domain.bed"
        bed.write_text("chr1\t0\t100\tGag\nchr2\t0\t100\tPol\n")
        files = {}
        global_rows, erv_rows, exposure_rows = [], [], []
        for chrom in ("1", "2"):
            path = root / ("chr" + chrom + ".tsv.gz")
            files[chrom] = path
            with gzip.open(path, "wt") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(hd.FIELDS)
                for sensitivity in ("all", "strict"):
                    for threshold in (5, 10):
                        observed = 20 if threshold == 5 else 0
                        old_sens = "exclude_combined" if sensitivity == "strict" else sensitivity
                        baseline = dict(chrom=chrom, context="AAA", af_bin="bin", an_bin="an", sensitivity=old_sens,
                                        population_mode="all", threshold=threshold, n_sites=100, n_multi=observed)
                        global_rows.append(baseline)
                        erv_rows.append(dict(baseline, label="family:ERVK"))
                        exposure_rows.append(dict(baseline, label="family:ERVK", block_id="EUR.38." + chrom + ".0001"))
                        for label, n, m in (("Gag", 25, observed // 2), ("Pol", 50, observed // 2), ("unannotated", 25, 0)):
                            if mismatch and chrom == "1" and sensitivity == "all" and threshold == 5 and label == "Gag":
                                n += 1
                            w.writerow([chrom, "EUR.38." + chrom + ".0001", "ERVK", label, sensitivity, threshold, "AAA", "bin", "an", n, m])
        for kind, rows in (("global", global_rows), ("erv", erv_rows), ("exposure", exposure_rows)):
            with gzip.open(modern / "summary" / ("architecture." + kind + "_strata.tsv.gz"), "wt") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
        args = SimpleNamespace(domain_build="38", domain_source="synthetic-reviewed-test", families="ERVK",
            domain_bed=str(bed), modern_dir=str(modern), out=str(root / "result"), bootstrap=2000, seed=8,
            min_blocks=2, chromosomes="1,2", fasta=None, erv_joined=None, segdup_bed=None, lcr_bed=None, low_mappability_bed=None)
        return args, files

    def test_two_chromosome_loco_and_original_counter_reconciliation(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args, files = self.fixture(root)
            original = dict(options=dict(ld_block_dir="test", architecture_build="38", thresholds=[5, 10]))
            with mock.patch.object(hd.hm, "targeted_inputs", return_value=(original, root, [("1", "vcf1"), ("2", "vcf2")], {"1", "2"})), mock.patch.object(hd, "scan_chrom", side_effect=lambda vcf, chrom, *rest: files[chrom]):
                result = hd.run(args)
            self.assertEqual(result["status"], "complete")
            rows = list(hd.hm.rows(root / "result/domain_enrichment.tsv"))
            total = next(r for r in rows if r["domain"] == "total" and r["threshold"] == "5" and r["sensitivity"] == "strict")
            self.assertEqual(int(total["observed_multi"]), 40)
            self.assertAlmostEqual(float(total["expected_multi"]), 200 * 20.5 / 101)
            self.assertAlmostEqual(float(total["oe"]), 40 / (200 * 20.5 / 101))
            self.assertEqual(result["comparison_scope"], "selected_chromosomes_subset")
            self.assertTrue(all(r["status"] == "matched" for r in hd.hm.rows(root / "result/domain_total_reconciliation.tsv")))
            self.assertTrue(result["output_inventory"])

    def test_count_drift_prevents_publishing_success(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args, files = self.fixture(root, mismatch=True)
            original = dict(options=dict(ld_block_dir="test", architecture_build="38", thresholds=[5, 10]))
            with mock.patch.object(hd.hm, "targeted_inputs", return_value=(original, root, [("1", "vcf1"), ("2", "vcf2")], {"1", "2"})), mock.patch.object(hd, "scan_chrom", side_effect=lambda vcf, chrom, *rest: files[chrom]):
                with self.assertRaisesRegex(ValueError, "disagree with original"):
                    hd.run(args)
            manifest = json.loads((root / "result/domain_manifest.json").read_text())
            self.assertEqual(manifest["status"], "failed")


if __name__ == "__main__":
    unittest.main()
