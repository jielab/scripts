"""The default main analysis excludes both GU jobs and cached IBDmix evidence."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
from analysis_scope import DEFAULT_STAGES, STAGES, select_stages
import batch_workflows as batch
import herv_pipeline as pipeline
import herv_modern as modern
import herv_primary as primary
import herv_report as report
import scan
from verification.test_archaic_summary_filter import CommonSupportReferences


class IBDmixOptional(unittest.TestCase):
    def test_default_native_and_batch_plans_do_not_discover_gu(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "out"
            argv = ["--out", str(out), "--scan-dir", str(Path(temp) / "scan"),
                    "--archaic-root", str(Path(temp) / "absent-GU")]
            args = pipeline.parser().parse_args(argv)
            self.assertFalse(args.ibdmix)
            self.assertEqual(pipeline.stage_selection(args), list(DEFAULT_STAGES))
            with patch("herv_archaic.discover", side_effect=AssertionError("Default accessed GU")):
                steps = pipeline.plan(args, pipeline.resolve(args), pipeline.stage_selection(args))
                plan, _ = batch.pipeline_plan(argv, "hpc", out, sys.executable)
            self.assertEqual([job["id"] for job in plan["jobs"]], ['all'])
            self.assertEqual(plan['provenance']['selected'], list(DEFAULT_STAGES))
            for step in steps:
                if step.name in {"primary", "summary", "report"}:
                    self.assertEqual(step.command[step.command.index("--ibdmix") + 1], "FALSE")
            self.assertFalse(plan["provenance"]["ibdmix_enabled"])

    def test_explicit_ibdmix_opt_in_is_required_and_selects_the_requested_scope(self):
        for command, subset in (("archaic", None), ("pipeline", "primary,archaic")):
            with self.assertRaisesRegex(ValueError, "requires --ibdmix TRUE"):
                select_stages(command, subset, False)
        self.assertEqual(select_stages(ibdmix=True), list(STAGES))
        self.assertEqual(select_stages("archaic", ibdmix=True), ["archaic"])
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "exploration"
            plan, _ = batch.pipeline_plan(["archaic", "--ibdmix", "TRUE", "--archaic-datasets", "ukb009"],
                                           "hpc", out, sys.executable)
            self.assertEqual([job["id"] for job in plan["jobs"]], ["archaic"])
            self.assertTrue(plan["provenance"]["ibdmix_enabled"])

    def test_summary_default_excludes_all_ibdmix_labels_before_fdr(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = CommonSupportReferences().fixture(root)
            out = root / "summary"
            with patch.object(modern, "plot_families"), patch.object(scan, "bh_fdr", wraps=scan.bh_fdr) as bh:
                modern.summarize_common_support(source, out)
            with (out / "common_support_enrichment.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["label"] for row in rows}, {"family:ERV1", "family:Altai.2013_like"})
            self.assertEqual(len(bh.call_args.args[0]), 2)
            manifest = json.loads((out / "summary_manifest.json").read_text())
            self.assertFalse(manifest["ibdmix_enabled"])
            self.assertEqual(manifest["selected_archaic_references"], [])
            self.assertGreater(manifest["excluded_archaic_counter_rows"], 0)

    def test_primary_default_filters_cached_ibdmix_before_statistics(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = CommonSupportReferences().fixture(root)
            results = {}
            for enabled in (False, True):
                out = root / str(enabled)
                with patch.object(primary, "plot_primary"):
                    manifest = primary.compute_primary_from_counters(source, out, include_ibdmix=enabled)
                with (out / "enrichment_summary.tsv").open() as handle:
                    results[enabled] = {row["label"]: row for row in csv.DictReader(handle, delimiter="\t")}
                self.assertEqual(manifest["ibdmix_enabled"], enabled)
                self.assertEqual(manifest["excluded_ibdmix_counter_rows"] > 0, not enabled)
            self.assertEqual(set(results[False]), {"family:ERV1", "family:Altai.2013_like"})
            self.assertEqual(set(results[True]), set(CommonSupportReferences.labels))
            for label, row in results[False].items():
                for field in ("observed_multi", "expected_multi", "oe", "bootstrap_p"):
                    self.assertEqual(row[field], results[True][label][field])

    def test_default_report_does_not_load_saved_archaic_results(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            archaic = root / "archaic"
            archaic.mkdir()
            table = "archaic_enrichment.tsv"
            (archaic / table).write_text("dataset\tlabel\tthreshold\tsensitivity\nARCHAIC_SENTINEL\t1kg::archaic:Altai:any\t10\texclude_combined\n")
            seen = []
            def validate(path, *args, **kwargs):
                seen.append(Path(path).name)
                if Path(path) == archaic:
                    return dict(status="complete", output_inventory={table: {}})
                raise ValueError("No current fixture output")
            with patch.object(report, "validate_manifest", side_effect=validate), \
                    patch.dict("os.environ", {"ALLELE_COHORT_TABLE": ""}):
                report.make_report(root, root / "scan", root)
                self.assertNotIn("archaic", seen)
                html = (root / "report.html").read_text()
                self.assertNotIn("ARCHAIC_SENTINEL", html)
                self.assertNotIn(table, html)
                report.make_report(root, root / "scan", root, include_ibdmix=True)
                self.assertIn("ARCHAIC_SENTINEL", (root / "report.html").read_text())


if __name__ == "__main__":
    unittest.main()
