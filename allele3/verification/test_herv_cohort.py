"""Counterexamples for cohort sampling, denominators and census provenance."""
import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

MODULE = Path(__file__).resolve().parents[1] / "f" / "herv_cohort.py"
spec = importlib.util.spec_from_file_location("herv_cohort", MODULE)
hc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hc)


def write_census(root, sample="S1", library="library-one", missing=0):
    root = Path(root)
    root.mkdir(parents=True)
    settings = dict(sample=sample, method="test", fasta_sha256="reference", library_sha256=library,
        repeatmasker_sha256="engine", rmblast_version="test", chunk_bp=10000, halo=1000,
        uncertain_flank=50, parameters=["-e", "rmblast"], code_sha256="code")
    rows = []
    for label in hc.LABELS:
        for category in hc.CATEGORIES:
            count = 800 if label == "diploid" else 790
            row = dict(label=label, category=category, union_bp=count, genome_bp=10000,
                       union_non_N_bp=count, genome_non_N_bp=10000 - missing)
            row["percent_genome"] = 100 * count / 10000
            row["percent_non_N"] = 100 * count / (10000 - missing)
            rows.append(row)
    table = root / "genome_fraction.tsv"
    with table.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (root / "run.json").write_text(json.dumps(dict(status="complete", scope=[f"chr{i}" for i in range(1, 23)], settings=settings)))
    (root / "analysis.done.json").write_text(json.dumps(dict(results_sha256=hc.sha256(table), settings=settings)))
    (root / "genome_fraction.json").write_text(json.dumps(dict(categories=rows, settings=settings)))
    return root


class CensusAuditTests(unittest.TestCase):
    def test_fractional_json_counts_are_rejected(self):
        row = dict(union_bp=8.5, genome_bp=100, union_non_N_bp=8, genome_non_N_bp=100,
                   percent_genome=8, percent_non_N=8)
        with self.assertRaisesRegex(ValueError, "Noninteger"):
            hc.validate_fraction(row)

    def test_valid_paired_fraction_and_no_missing_as_zero(self):
        with tempfile.TemporaryDirectory() as temp:
            root = write_census(Path(temp) / "S1", missing=2000)
            rows, comparison = hc.inspect_census(root, "S1")
            self.assertEqual(rows[0]["model_percent"], 8)
            self.assertEqual(rows[0]["model_percent_non_N"], 10)
            self.assertAlmostEqual(rows[0]["delta_percentage_points"], 0.125)
            self.assertEqual(comparison["chromosomes"], [str(i) for i in range(1, 23)])

    def test_table_corruption_fails_receipt(self):
        with tempfile.TemporaryDirectory() as temp:
            root = write_census(Path(temp) / "S1")
            table = root / "genome_fraction.tsv"
            table.write_text(table.read_text().replace("800", "900"))
            with self.assertRaisesRegex(ValueError, "SHA256"):
                hc.inspect_census(root, "S1")

    def test_json_changes_cannot_hide_behind_tsv_hash(self):
        with tempfile.TemporaryDirectory() as temp:
            root = write_census(Path(temp) / "S1")
            p = root / "genome_fraction.json"
            data = json.loads(p.read_text())
            r = data["categories"][0]
            r.update(union_bp=900, union_non_N_bp=900, percent_genome=9.0, percent_non_N=9.0)
            p.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "JSON and TSV"):
                hc.inspect_census(root, "S1")

    def test_different_library_is_not_pooled(self):
        with tempfile.TemporaryDirectory() as temp:
            a, _ = hc.inspect_census(write_census(Path(temp) / "a", "A", "one"), "A")
            b, _ = hc.inspect_census(write_census(Path(temp) / "b", "B", "two"), "B")
            self.assertNotEqual(a[0]["comparison_id"], b[0]["comparison_id"])

    def test_same_scientific_inputs_different_sample_share_comparison(self):
        with tempfile.TemporaryDirectory() as temp:
            a, _ = hc.inspect_census(write_census(Path(temp) / "a", "A"), "A")
            b, _ = hc.inspect_census(write_census(Path(temp) / "b", "B"), "B")
            self.assertEqual(a[0]["comparison_id"], b[0]["comparison_id"])

    def test_wrong_sample_and_failed_run_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = write_census(Path(temp) / "S1")
            with self.assertRaisesRegex(ValueError, "sample differs"):
                hc.inspect_census(root, "S2")
            p = root / "run.json"
            run = json.loads(p.read_text())
            run["status"] = "failed"
            p.write_text(json.dumps(run))
            with self.assertRaisesRegex(ValueError, "not complete"):
                hc.inspect_census(root, "S1")

    def test_percentage_with_wrong_denominator_fails(self):
        row = dict(union_bp=8, genome_bp=100, union_non_N_bp=8, genome_non_N_bp=80,
                   percent_genome=8, percent_non_N=8)
        with self.assertRaisesRegex(ValueError, "percent_non_N"):
            hc.validate_fraction(row)

    def test_unavailable_roster_writes_no_success(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            meta = root / "metadata.tsv"
            meta.write_text("sample_id\tpopulation\tsex\nS1\tEUR\tF\n")
            out = root / "summary"
            args = SimpleNamespace(metadata=meta, group_column="population", sex_column="sex",
                                   census_root=root / "censuses", out=out)
            with self.assertRaisesRegex(ValueError, "No verified"):
                hc.summarize(args)
            self.assertEqual(json.loads((out / "manifest.json").read_text())["status"], "unavailable")
            self.assertIn("unavailable", (out / "sample_audit.tsv").read_text())


class SamplingTests(unittest.TestCase):
    def setUp(self):
        self.rows = [dict(sample_id=f"{g}{s}{i}", population=g, sex=s, include="1", unrelated="1")
                     for g in ("AFR", "EUR") for s in ("F", "M") for i in range(3)]

    def test_sampling_is_balanced_reproducible_and_row_order_independent(self):
        a, support = hc.balanced_subset(self.rows, total=8, per_stratum=2)
        b, _ = hc.balanced_subset(list(reversed(self.rows)), total=8, per_stratum=2)
        self.assertEqual(a, b)
        self.assertEqual([r["selected"] for r in support], [2] * 4)

    def test_ineligible_forced_sample_cannot_be_selected(self):
        self.rows[0]["include"] = "0"
        with self.assertRaisesRegex(ValueError, "eligible"):
            hc.balanced_subset(self.rows, always=(self.rows[0]["sample_id"],))

    def test_forced_sample_is_counted_not_duplicated(self):
        sample = self.rows[0]["sample_id"]
        chosen, _ = hc.balanced_subset(self.rows, total=4, per_stratum=1, always=(sample,))
        self.assertEqual(len(chosen), 4)
        self.assertEqual(sum(r["sample_id"] == sample for r in chosen), 1)
        self.assertEqual(len({(r["population"], r["sex"]) for r in chosen}), 4)

    def test_known_relatives_excluded(self):
        self.rows[0]["unrelated"] = "0"
        chosen, _ = hc.balanced_subset(self.rows, total=20, per_stratum=3)
        self.assertNotIn(self.rows[0]["sample_id"], {r["sample_id"] for r in chosen})

    def test_duplicate_metadata_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "m.tsv"
            p.write_text("sample_id\tpopulation\tsex\nS1\tEUR\tF\nS1\tEUR\tF\n")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                hc.metadata(p)

    def test_existing_1kg_panel_no_reformatting(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "1kg.tsv"
            p.write_text("sample\tpop\tsuper_pop\tgender\nHG00096\tGBR\tEUR\tmale\n")
            row = hc.metadata(p)[0]
            self.assertEqual((row["sample_id"], row["population"], row["sex"]), ("HG00096", "EUR", "male"))
            self.assertEqual(hc.metadata(p, "pop")[0]["population"], "GBR")


class CensusReferenceLabelTests(unittest.TestCase):
    def test_37_and_38_report_the_actual_reference(self):
        sys.path.insert(0, str(MODULE.parent))
        spec = importlib.util.spec_from_file_location("herv_individual_test", MODULE.parent / "1.herv.py")
        engine = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = engine
        spec.loader.exec_module(engine)
        for build in (37, 38):
            with self.subTest(build=build), tempfile.TemporaryDirectory() as temp:
                root = write_census(Path(temp) / "S1")
                p = root / "genome_fraction.json"
                data = json.loads(p.read_text())
                data.update(elapsed_wall_seconds=1, measured_cpu_seconds={}, boundary_disagreement_bp={}, cohort_projection={})
                p.write_text(json.dumps(data))
                (root / "variant_summary.tsv").write_text("uncertain_reference_bp\tunhandled_records\n0\t0\n")
                fasta = Path(f"/references/GRCH{build}.fasta.gz")
                args = SimpleNamespace(out=root, fasta=fasta, sample="S1", effective_input="vcf", jobs=1)
                with mock.patch.object(engine.os, "sched_setaffinity"), mock.patch.object(engine.os, "nice"), mock.patch.object(engine.subprocess, "run"):
                    engine.execute(args, [["mock-engine"]], [0])
                report = (root / "report.txt").read_text()
                self.assertIn(str(fasta), report)
                self.assertNotIn("GRCh38 autosomal", report)


if __name__ == "__main__":
    unittest.main()
