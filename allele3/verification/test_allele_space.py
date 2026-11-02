"""Scientific regressions for allele-space; run with unittest (no pytest needed)."""
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
import allele_space as space
from analysis_integrity import identity, validate_manifest


def sample(sid, age, x=0., y=0., cohort="ancient", region="E"):
    return dict(sample_id=sid, age_bp=age, PC1=x, PC2=y, cohort=cohort, region=region, eligible=True)


def info(kind="direct_gt", chrom="1"):
    return dict(site_id="site1", chrom=chrom, pos=100, genome_build="GRCh38", allele="G",
                role="secondary", evidence_kind=kind, counting_scope="explicit_count_of_one_nucleotide")


class AlleleSpace(unittest.TestCase):
    def test_research_dispatch_space_help(self):
        r = subprocess.run([sys.executable, str(ROOT/"f/research.py"), "space", "--help"], text=True, capture_output=True)
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assertIn("--time-dir", r.stdout)
        self.assertIn("--alleles-dir", r.stdout)
        self.assertIn("--observations", r.stdout)

    def test_time_stage_cli_alias_preserves_adapter_contract(self):
        common = ["--metadata", "metadata.tsv", "--coordinates", "all.sscore", "--out", "result", "--pgen-prefix", "panel"]
        for flag in ("--time-dir", "--alleles-dir"):
            args = space.parser().parse_args([flag, "pipeline/time", *common])
            self.assertEqual(args.alleles_dir, "pipeline/time")

    def test_unassayed_is_not_zero_and_sex_ploidy_is_actual(self):
        samples = [sample("m", 100), sample("f", 200), sample("unknown", 300)]
        cells, _ = space.build_cells(samples, 1, 1000, 10000)
        rows, summary = space.summarize(info(chrom="X"), np.array([1., 0., np.nan]),
            np.array([1., 2., np.nan]), np.array([True, True, False]), samples, cells, 1, 1)
        self.assertEqual(rows[0]["AN"], 3)
        self.assertAlmostEqual(rows[0]["AF"], 1/3)
        self.assertEqual(rows[0]["n_called"], 2)
        self.assertEqual(rows[0]["n_declared_assayed"], 2)
        self.assertEqual(rows[0]["call_rate_among_declared"], 1)
        self.assertAlmostEqual(rows[0]["frame_measurement_fraction"], 2/3)
        self.assertEqual(summary[0]["origin_inference"], "not_estimated")

    def test_frequency_peak_is_not_largest_carrier_count(self):
        # A densely sampled cell has 10 copies/200 chromosomes; a sparse cell
        # has 4 copies/8 chromosomes. Frequency correctly peaks in the latter.
        samples = [sample(f"a{i}", 100, 0) for i in range(100)] + [sample(f"b{i}", 100, 1) for i in range(4)]
        counts = np.r_[np.ones(10), np.zeros(90), np.ones(4)]
        cells, _ = space.build_cells(samples, 2, 1000, 10000)
        rows, summary = space.summarize(info(), counts, np.full(104, 2.), np.ones(104, dtype=bool), samples, cells, 4, 2)
        self.assertEqual(summary[0]["peak_AF_cell"], "0:1:0")
        self.assertEqual(summary[0]["peak_AF"], .5)
        self.assertEqual(summary[0]["peak_AF_AN"], 8)
        self.assertEqual({r["AC"] for r in rows}, {10., 4.})

    def test_oldest_singleton_and_supported_interval_are_separate(self):
        samples = [sample("old", 9000)] + [sample(f"young{i}", 2300+i) for i in range(4)]
        cells, _ = space.build_cells(samples, 1, 1000, 10000)
        _, summary = space.summarize(info(), np.array([1., 1., 1., 0., 0.]), np.full(5, 2.),
            np.ones(5, dtype=bool), samples, cells, 4, 2)
        row = summary[0]
        self.assertEqual(row["oldest_observed_carrier_age_bp"], 9000)
        self.assertEqual(row["earliest_supported_time_bin"], 2)
        self.assertEqual(row["earliest_supported_AN"], 8)
        self.assertEqual(row["earliest_supported_n_carriers"], 2)

    def test_posterior_dosage_is_not_called_carriage(self):
        samples = [sample("a", 5000), sample("b", 1000)]
        cells, _ = space.build_cells(samples, 1, 1000, 10000)
        rows, summary = space.summarize(info("posterior_dosage"), np.array([.00001, .9]),
            np.full(2, 2.), np.ones(2, dtype=bool), samples, cells, 1, 1)
        self.assertIsNone(summary[0]["n_hardcall_carriers"])
        self.assertNotIn("oldest_observed_carrier_age_bp", summary[0])
        self.assertTrue(all(r["n_hardcall_carriers"] is None for r in rows))

    def test_bin_last_endpoint_and_modern_remain_distinct(self):
        samples = [sample("old", 10000), sample("ancient0", 0), sample("modern", None, cohort="modern")]
        cells, _ = space.build_cells(samples, 1, 1000, 10000)
        self.assertEqual({c["time_bin"] for c in cells}, {0, 9, "modern"})
        modern = next(c for c in cells if c["time_bin"] == "modern")
        self.assertIsNone(modern["frame_age_bp"])

    def test_observations_reject_duplicates_and_inconsistent_dosage(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"input.tsv"
            header = "site_id\tchrom\tpos\tallele\tsample_id\tcount\tploidy\n"
            row = "s\tY\t1\tA\ta\t1\t1\n"
            p.write_text(header+row+row)
            with self.assertRaisesRegex(ValueError, "Duplicate sample"):
                list(space.observation_groups(p, [sample("a", 100)], "GRCh38", "direct_gt"))
            p.write_text(header+"s\tY\t1\tA\ta\t2\t1\n")
            with self.assertRaisesRegex(ValueError, "outside ploidy"):
                list(space.observation_groups(p, [sample("a", 100)], "GRCh38", "direct_gt"))

    def test_each_allele_retains_its_own_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"input.tsv"
            p.write_text("site_id\tchrom\tpos\tallele\tsample_id\tcount\tploidy\n"
                "s\t1\t1\tA\ta\t1\t2\n"
                "s\t1\t1\tA\tb\t0\t2\n"
                "s\t1\t1\tG\ta\t1\t2\n")
            ss = [sample("a", 100), sample("b", 100)]
            cells, _ = space.build_cells(ss, 1, 1000, 10000)
            result = [space.summarize(i, c, p, d, ss, cells, 1, 1)[1][0]
                for i, c, p, d in space.observation_groups(p, ss, "GRCh38", "direct_gt")]
            self.assertEqual([r["AN"] for r in result], [4, 2])
            self.assertEqual([r["AF"] for r in result], [.25, .5])

    def test_date_reference_checks_ids_and_reports_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"dates.tsv"
            p.write_text("sample_id\tage_bp\tregion\nB\t500\tE\nA\t3000\tE\n")
            rows, audit = space.date_audit([sample("A", 3000), sample("B", 1500), sample("C", 100)], p, 1)
            self.assertEqual([r["status"] for r in rows], ["concordant", "date_mismatch", "reference_missing"])
            self.assertEqual(audit["status"], "metadata_concordance_only")

    def test_pgen_adapter_validates_panel_and_reorders_by_psam(self):
        try:
            import pgenlib
        except ImportError:
            self.skipTest("pgenlib is needed for the native multiallelic integration test")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); stage = root/"time"; stage.mkdir(); prefix = root/"panel"
            for suffix, content in ((".pvar", "#CHROM\tPOS\tID\tREF\tALT\n1\t1\tx\tA\tC,G,T\n"),
                                    (".psam", "#IID\nB\nA\nC\nD\n")):
                Path(str(prefix)+suffix).write_text(content)
            writer = pgenlib.PgenWriter(str(prefix).encode()+b".pgen", 4, variant_ct=1, allele_ct_limit=4)
            writer.append_alleles(np.array([1, 2, -9, -9, 2, 2, 0, 3], dtype=np.int32), allele_ct=4)
            writer.close()
            entries = []
            for base, index in (("C", 1), ("G", 2), ("T", 3)):
                entries.append(dict(allele_id="s:"+base, site_id="s", chrom="1", pos=1, allele=base, role="secondary",
                    counting_scope="native_joint_nucleotide_count", measurement_plan=dict(status="assayed_pass",
                        mode="pgen_allele", pvar_index=0, allele_idx=index, invert_dosage=0)))
            plans = stage/"measurement_plans.json"
            plans.write_text(json.dumps(dict(schema="nucleotide-count-1", ALT_measurements=entries)))
            inputs = {str(Path(str(prefix)+suffix).resolve()): identity(str(prefix)+suffix) for suffix in (".pgen", ".pvar", ".psam")}
            output_id = identity(plans)
            (stage/"manifest.json").write_text(json.dumps(dict(status="complete", inputs=inputs,
                output_inventory={"measurement_plans.json": {k: v for k, v in output_id.items() if k != "path"}})))
            ss = [sample("A", 100), sample("B", 200), sample("C", 300), sample("D", 400)]
            result = list(space.pgen_groups(stage, prefix, ss, 2))
            self.assertTrue(np.isnan(result[0][1][0]))
            np.testing.assert_equal(result[0][1], [np.nan, 1, 0, 0])
            np.testing.assert_equal(result[1][1], [np.nan, 1, 2, 0])
            np.testing.assert_equal(result[2][1], [np.nan, 0, 0, 1])
            Path(str(prefix)+".psam").write_text("#IID\nA\nB\nC\nD\n")
            with self.assertRaisesRegex(ValueError, "Original panel changed"):
                list(space.pgen_groups(stage, prefix, ss, 1))

    def test_cli_fixture_outputs_verified_manifest_and_no_fake_origin(self):
        fixture = ROOT/"verification/fixtures/allele_space"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/"result"
            command = [sys.executable, str(ROOT/"f/allele_space.py"), "--observations", str(fixture/"observations.tsv"),
                "--metadata", str(fixture/"metadata.tsv"), "--coordinates", str(fixture/"coordinates.sscore"),
                "--reference-metadata", str(fixture/"reference_dates.tsv"), "--genome-build", "GRCh38",
                "--evidence-kind", "direct_gt", "--allow-unverified-coordinates", "--pc-bins", "2",
                "--min-samples", "2", "--out", str(out), "--no-plots"]
            r = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
            manifest = validate_manifest(out, require_evidence=True)
            self.assertEqual(manifest["status"], "complete_descriptive_unverified")
            self.assertEqual(manifest["n_alleles"], 4)
            rows = list(space.read_tsv(out/"allele_space_summary.tsv.gz"))
            self.assertTrue(all(r["origin_inference"] == "not_estimated" for r in rows))
            y = next(r for r in rows if r["chrom"] == "Y" and r["cohort"] == "ancient")
            self.assertEqual(float(y["AN"]), 3)
            dispatched = [sys.executable, str(ROOT/"f/research.py"), "space", *command[2:]]
            r = subprocess.run(dispatched, text=True, capture_output=True)
            self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
            self.assertEqual(validate_manifest(out, require_evidence=True)["n_alleles"], 4)
            # The explicit descriptive flag is mandatory without old receipts.
            command.remove("--allow-unverified-coordinates")
            r = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(r.returncode, 0)


if __name__ == "__main__":
    unittest.main()
