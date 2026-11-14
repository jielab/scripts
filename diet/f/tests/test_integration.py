#!/usr/bin/env python3
"""Small aggregate-data and coordinator regression suite; no participant data."""
from pathlib import Path
import contextlib
import io
import json
import math
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE))
import results
import pipeline
import web
from score_registry import SCORES


class Contracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def workflow(self, *args):
        with patch.dict(os.environ, {"DIET_ENGINE_ENV": "{}"}):
            return pipeline.Workflow(pipeline.parser().parse_args([
                *args, "--dir-output", str(self.out), "--dir-charls", str(self.out / "absent"), "--dry-run"]))

    def test_bayes_point_null_and_equivalence(self):
        b = results.bf_row(0, .005, 1.10, 1.05)
        tau = math.log(1.10) / 1.96
        self.assertAlmostEqual(b["BF01"], math.sqrt(1 + tau ** 2 / .005 ** 2))
        self.assertAlmostEqual(b["BF01"] * b["BF10"], 1)
        self.assertEqual(b["evidence"], "point_null")
        self.assertTrue(b["equivalent"])
        a = results.bf_row(math.log(.80), .02, 1.10, 1.05)
        self.assertEqual(a["evidence"], "alternative")
        self.assertFalse(a["equivalent"])
        with self.assertRaises(ValueError):
            results.bf_row(0, 0)

    def test_nhanes_standard_error_is_not_reconstructed_from_t_ci(self):
        r = {"estimate": ".8", "lower": ".6", "upper": "1.1"}
        self.assertIsNone(results.normalise(r, "dose", "nhanes", "source")["se"])
        self.assertIsNotNone(results.normalise(r, "dose", "ukb", "source")["se"])
        r.update(beta="-.22314355", se=".12", df="12", contrast="per_SD", family="main")
        n = results.normalise(r, "dose", "nhanes", "source")
        self.assertEqual(n["se"], .12)
        self.assertEqual(n["beta"], -.22314355)
        self.assertEqual(n["effect_type"], "OR")
        self.assertEqual(n["family"], "main")

    def test_selected_profiles_and_distinct_nhanes_roots(self):
        w = self.workflow("all", "--study", "dose", "--dir-files", str(self.out / "files"),
                          "--dir-nhanes", str(self.out / "maha"), "--dir-nhanes-dose", str(self.out / "dose"))
        self.assertEqual(w.selected(), ["ukb", "nhanes"])
        self.assertEqual(w.env["DIET_NHANES"], str(self.out / "maha"))
        self.assertEqual(w.env["NHANES_DATA"], str(self.out / "dose"))
        with self.assertRaises(ValueError):
            self.workflow("chns", "--study", "dose").selected()

    def test_native_modules_reject_incompatible_selection(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as err:
            pipeline.main(["ml", "--study", "maha", "--cohort", "nhanes", "--dry-run"])
        self.assertEqual(err.exception.code, 2)

    def test_disco_cohort_list_matches_all_and_selects_score_engines(self):
        with contextlib.redirect_stdout(io.StringIO()):
            all_cohorts = self.workflow("--entry", "disco", "all").selected()
            listed = self.workflow("--entry", "disco", "ukb,nhanes,chns,charls").selected()
        self.assertEqual(listed, all_cohorts)
        w = self.workflow("--entry", "disco", "ukb,nhanes", "--score", "dose,dose9half")
        self.assertEqual(w.selected(), ["ukb", "nhanes"])
        self.assertEqual(w.profiles("ukb"), ["dose"])
        w = self.workflow("--entry", "disco", "ukb", "--score", "upf")
        self.assertEqual(w.profiles("ukb"), ["maha"])
        with patch.object(w, "r") as run:
            w.disco("ukb")
        self.assertEqual(run.call_args.kwargs["extra"],
                         {"DIET_SCORE_SELECTION": "upf", "DIET_SCORE_PROFILES": "maha"})

    def test_unsupported_score_is_logged_without_running_an_engine(self):
        with contextlib.redirect_stdout(io.StringIO()), patch.object(pipeline.Workflow, "r") as run:
            code = pipeline.main(["--entry", "disco", "nhanes", "--score", "upf",
                                  "--dir-output", str(self.out)])
        self.assertEqual(code, 0)
        run.assert_not_called()
        self.assertIn("SKIP nhanes: upf", (self.out / "logs/disco.nhanes.skip.log").read_text())

    def test_scoring_options_reject_typos_before_workflow_creation(self):
        for argv in (["--entry", "disco", "ukb,bogus"],
                     ["--entry", "disco", "ukb", "--score", "dash,bogus"],
                     ["--entry", "disco", "ukb", "--score", "dash", "--study", "dose"],
                     ["ukb", "--score", "dash"]):
            with self.subTest(argv=argv), patch.object(pipeline, "Workflow") as workflow:
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    pipeline.main(argv)
                workflow.assert_not_called()

    def test_public_score_example_lists_every_registered_score(self):
        help_text = (CODE.parent / "1.disco.sh").read_text().split("Examples:\n")[1].split("PUBLIC_ENTRY_HELP")[0]
        line = next(line for line in help_text.splitlines() if "--score " in line)
        self.assertEqual(line.split("--score ")[1].split()[0].split(","), list(SCORES))
        self.assertNotIn("--plan-only", help_text)

    def test_worker_relative_outputs_stay_in_analysis_directory(self):
        w = self.workflow("--entry", "disco", "ukb")
        w.a.dry_run = False
        with contextlib.redirect_stdout(io.StringIO()):
            w.run("test/output_directory", [sys.executable, "-c",
                  "from pathlib import Path; Path('relative-output.txt').write_text(str(Path.cwd()))"])
        self.assertEqual((self.out / "relative-output.txt").read_text(), str(self.out))

    def test_final_scope_and_collect_only(self):
        w = self.workflow("final", "--study", "maha", "--cohort", "ukb")
        with contextlib.redirect_stdout(io.StringIO()) as output:
            w.final()
        self.assertIn("run.R ukb final", output.getvalue())
        self.assertNotIn("run.R all final", output.getvalue())
        w = self.workflow("final", "--collect-only")
        with contextlib.redirect_stdout(io.StringIO()) as output:
            w.final()
        self.assertNotIn("dose/final.py", output.getvalue())
        self.assertNotIn("maha/run.R", output.getvalue())

    def test_aggregate_guard_reads_only_header(self):
        p = self.out / "big.csv"
        p.write_text("eid,estimate\n0001,0.8\n")
        with patch.object(results, "rows", side_effect=AssertionError("must not read full participant table")):
            self.assertFalse(results.is_aggregate_table(p))
            p.write_text("term,estimate\nx,0.8\n")
            self.assertTrue(results.is_aggregate_table(p))

    def test_exact_paper_mapping_partial_and_stale_cleanup(self):
        native = self.out / "dose/ukb/05mechanism"
        native.mkdir(parents=True)
        unrelated = native / "Figure4_MRI.png"
        unrelated.write_bytes(b"synthetic test payload")
        records = results.collect_final(self.out, ["dose"])
        paper = {r["manuscript_id"]: r for r in records if r["kind"] == "manuscript_item"}
        self.assertEqual(len(paper), 5)
        self.assertTrue(all(r["status"] == "not_generated" for r in paper.values()))
        final = self.out / "final/dose"
        (final / "Fig3.DOSE.score.png").write_bytes(b"synthetic test payload")
        results.write_rows(final / "composite_status.tsv", [dict(figure="Fig3", status="partial", missing_panels="f")])
        records = results.collect_final(self.out, ["dose"])
        item = next(r for r in records if r.get("manuscript_id") == "Fig3")
        self.assertEqual(item["status"], "partial")
        self.assertIn("Missing panels: f", item["detail"])
        owned = final / "ukb__05mechanism__Figure4_MRI.png"
        self.assertTrue(owned.is_file())
        unrelated.unlink()
        results.collect_final(self.out, ["dose"])
        self.assertFalse(owned.exists())
        self.assertTrue((final / "Fig3.DOSE.score.png").exists())

    def test_bayes_update_preserves_unselected_results(self):
        effects = [dict(study="dose", cohort="ukb", beta=0.0, se=.01, status="ok"),
                   dict(study="maha", cohort="nhanes", beta=-.05, se=.02, status="ok")]
        with patch.object(results, "collect_effects", return_value=effects):
            results.bayes(self.out, ["nhanes"], ["maha"], [1.1], [1.05])
            results.bayes(self.out, ["ukb"], ["dose"], [1.1], [1.05])
        self.assertEqual(len(results.rows(self.out / "bayes/effects.tsv")), 2)

    def test_report_safe_embedding_and_empty_state(self):
        results.write_rows(self.out / "assoc/ukb.dose.tsv", [dict(score="</script><script>bad()</script>",
                          outcome="stroke", beta="-.1", se=".02", estimate=".9048")])
        target = web.render(self.out)
        doc = target.read_text()
        self.assertNotIn("</script><script>bad()", doc)
        self.assertIn("\\u003c/script\\u003e", doc)
        self.assertNotIn("__DIET_DATA__", doc)
        self.assertNotIn("<script src=", doc)
        self.assertEqual(len(web.build_data(self.out)["effects"]), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
