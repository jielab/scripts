"""Optional benchmark comparison must use identical frozen people and labels."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from count_abm.benchmarks import compare_current, load_benchmark
from count_abm.inputs import file_hash, json_write, read_table
from count_abm.workflow import main
from test_count_abm import exported_fixture


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.root = self.base/"analysis"
        self.manifest, self.trait = exported_fixture(self.base/"input")
        self.args = ["--traits", self.trait, "--analysis-root", str(self.root),
                     "--input-manifest", str(self.manifest), "--device", "cpu",
                     "--radius-k-ref", "8", "--min-count", "5", "--min-group-ess", "3",
                     "--bootstrap", "4", "--threads", "1"]
        with contextlib.redirect_stdout(io.StringIO()):
            main(["all", *self.args])
        self.count_dir = self.root/"count_abm"/self.trait
        self.count = read_table(self.count_dir/"_LOCAL_ONLY/test.predictions.tsv.gz")
        self.truth = read_table(self.count_dir/"_LOCAL_ONLY/test.truth.tsv.gz")

    def tearDown(self):
        self.temporary.cleanup()

    def benchmark(self, *, drop=False, shift_label=False, missing_label=False):
        out = self.root/"csx"/self.trait/"targetfit"
        local = out/"_LOCAL_ONLY"
        local.mkdir(parents=True, exist_ok=True)
        pred = self.count[["id", "group_id"]].copy()
        pred["role"] = "test"
        pred["populations.prediction"] = self.count.p_population+.1
        truth = self.truth.copy()
        if drop:
            pred = pred.iloc[:-1]
        if shift_label:
            truth.loc[0, "y"] += 1
        if missing_label:
            truth.loc[0, "y"] = np.nan
        p, t = local/"final_prs.tsv.gz", local/"test_truth.tsv.gz"
        pred.sample(frac=1, random_state=5).to_csv(p, sep="\t", index=False)
        truth.sample(frac=1, random_state=6).to_csv(t, sep="\t", index=False)
        receipt = dict(schema="grid.target_benchmark.1", status="completed", method="csx",
            trait=self.trait, task="continuous", endpoint=self.trait, models=["populations"],
            model_labels={"populations": "PRS-CSx validation combination"},
            predictions_file=str(p), test_truth_file=str(t),
            test_labels_used_only_after_prediction_freeze=True,
            roster={"actual_validation_N": 50}, files={str(p): file_hash(p), str(t): file_hash(t)})
        json_write(out/"COMPLETE.json", receipt)
        return p, out/"COMPLETE.json"

    def test_identical_ids_labels_align_by_identity_and_compare_both_populations(self):
        self.benchmark()
        comparisons, status = compare_current(self.root, [self.trait], repetitions=4)
        self.assertEqual(len(comparisons), 2)
        self.assertEqual(status.loc[status.benchmark_family.eq("csx"), "status"].iloc[0],
                         "comparable_exact_test_roster_and_labels")
        all_people = comparisons.loc[comparisons.population.eq("all_with_explicit_fallback")].iloc[0]
        supported = comparisons.loc[comparisons.population.eq("count_supported_same_people")].iloc[0]
        self.assertEqual(all_people.n_same_people, len(self.count))
        self.assertEqual(all_people.n_evaluable, len(self.count)-1)
        self.assertEqual(supported.n_same_people, int(self.count.supported.sum()))
        frame = self.count.merge(self.truth, on="id")
        expected = float(np.nanmean((frame.y-frame.p_population)**2-(frame.y-frame.p_population-.1)**2))
        self.assertAlmostEqual(all_people.count_minus_benchmark_loss, expected)
        with contextlib.redirect_stdout(io.StringIO()):
            main(["final", *self.args])
        saved = read_table(self.root/"count_abm/final/count_vs_benchmarks.tsv")
        self.assertEqual(len(saved), 2)
        self.assertEqual(saved.benchmark_label.iloc[0], "PRS-CSx validation combination")

    def test_never_uses_smaller_intersection_or_changed_labels(self):
        for options, reason in ((dict(drop=True), "roster differs"),
                                (dict(shift_label=True), "labels or missing-label"),
                                (dict(missing_label=True), "labels or missing-label")):
            with self.subTest(options=options):
                self.benchmark(**options)
                comparisons, statuses = compare_current(self.root, [self.trait], repetitions=4)
                self.assertTrue(comparisons.empty)
                row = statuses.loc[statuses.benchmark_family.eq("csx")].iloc[0]
                self.assertEqual(row.status, "not_comparable")
                self.assertIn(reason, row.reason)

    def test_changed_artifact_and_legacy_only_outputs_are_not_compared(self):
        prediction, receipt = self.benchmark()
        with prediction.open("ab") as stream:
            stream.write(b"changed after freeze")
        comparisons, statuses = compare_current(self.root, [self.trait], repetitions=4)
        self.assertTrue(comparisons.empty)
        self.assertEqual(statuses.loc[statuses.benchmark_family.eq("csx"), "status"].iloc[0], "not_comparable")
        receipt.rename(receipt.with_name("performance.json"))
        comparisons, statuses = compare_current(self.root, [self.trait], repetitions=4)
        self.assertTrue(comparisons.empty)
        self.assertEqual(statuses.loc[statuses.benchmark_family.eq("csx"), "status"].iloc[0], "missing_current_compatible_receipt")

    def test_production_csx_source_cannot_be_silently_changed(self):
        _, receipt = self.benchmark()
        count = self.count.merge(self.truth, on="id", validate="one_to_one")
        with self.assertRaisesRegex(ValueError, "CSx score provenance differs"):
            load_benchmark(receipt, method="csx", trait=self.trait, task="continuous",
                           endpoint=self.trait, count_frame=count, count_csx_score_sha256="a different frozen source")


if __name__ == "__main__":
    unittest.main()
