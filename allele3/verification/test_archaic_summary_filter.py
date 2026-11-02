"""Reference selection changes the tested family, never the cached site counts."""
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"f"))
import herv_modern as hm
import scan


class CommonSupportReferences(unittest.TestCase):
    labels = ("family:ERV1", "family:Altai.2013_like", "1kg::archaic:Altai:any",
              "1kg::archaic:Altai.2013:any", "ukb009::archaic:Denisova:any")

    def fixture(self, root):
        modern = root/"modern"
        summary = modern/"summary"
        summary.mkdir(parents=True)
        base = dict(population_mode="all", threshold="5", context="ACA", af_bin="common",
                    an_bin="high", sensitivity="exclude_combined", n_sites="1000", n_multi="100")
        global_rows = [dict(base, chrom=str(chrom)) for chrom in (1, 2)]
        exposure_rows = [dict(base, chrom=str(chrom), label=label, block_id=str(block),
                              n_sites="10", n_multi=str((block+index) % 3))
                         for index, label in enumerate(self.labels) for chrom in (1, 2) for block in range(10)]
        exposure_rows += [dict(base, chrom="1", label="1kg::archaic:Altai.2013:any", block_id="outside_partition")]
        # These are outside the summary's defined scope, and must not inflate the audit.
        exposure_rows += [dict(base, chrom="X", label="1kg::archaic:Altai.2013:any", block_id="outside_partition"),
                          dict(base, chrom="1", population_mode="EUR", label="1kg::archaic:Altai.2013:any", block_id="outside_partition")]
        for kind, data in (("global", global_rows), ("exposure", exposure_rows)):
            with gzip.open(summary/f"architecture.{kind}_strata.tsv.gz", "wt") as handle:
                writer = csv.DictWriter(handle, list(data[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(data)
        return modern

    def read_results(self, out):
        with (out/"common_support_enrichment.tsv").open() as handle:
            return {row["label"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def summarize(self, modern, out, references=None):
        with patch.object(hm, "plot_families"), patch.object(scan, "bh_fdr", wraps=scan.bh_fdr) as bh, \
                patch.object(scan, "load_run", side_effect=AssertionError("No VCF scan or run metadata is needed")):
            hm.summarize_common_support(modern, out, archaic_references=references)
        return self.read_results(out), bh.call_args.args[0]

    def test_default_five_filter_precedes_bh_and_has_explicit_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); modern = self.fixture(root); out = root/"result"
            results, tested = self.summarize(modern, out)
            self.assertEqual(set(results), set(self.labels)-{"1kg::archaic:Altai.2013:any"})
            self.assertEqual(len(tested), 4)
            self.assertTrue(all(row["status"] == "complete" for row in results.values()))
            manifest = json.loads((out/"summary_manifest.json").read_text())
            self.assertEqual(manifest["selected_archaic_references"], list(scan.DEFAULT_ARCHAIC_REFERENCES))
            self.assertEqual(manifest["excluded_archaic_counter_rows"], 21)
            self.assertEqual(manifest["excluded_archaic_labels"], ["1kg::archaic:Altai.2013:any"])
            with (out/"archaic_reference_selection.tsv").open() as handle:
                audit = list(csv.DictReader(handle, delimiter="\t"))
            old = next(row for row in audit if row["reference"] == "Altai.2013")
            self.assertEqual((old["decision"], old["counter_rows"]), ("excluded", "21"))

    def test_explicit_historical_selection_preserves_other_point_estimates_and_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); modern = self.fixture(root)
            before = {path: path.read_bytes() for path in (modern/"summary").glob("*.gz")}
            current, current_tests = self.summarize(modern, root/"current")
            history, history_tests = self.summarize(modern, root/"history", [*scan.DEFAULT_ARCHAIC_REFERENCES, "Altai.2013"])
            self.assertEqual(set(history), set(self.labels))
            self.assertEqual((len(current_tests), len(history_tests)), (4, 5))
            for label in current:
                for field in ("observed_multi", "expected_multi", "oe", "supported_sites", "bootstrap_p", "ci_low", "ci_high"):
                    self.assertEqual(current[label][field], history[label][field])
            self.assertEqual(before, {path: path.read_bytes() for path in before})

    def test_exact_reference_names_do_not_match_prefixes_or_nonarchaic_labels(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); modern = self.fixture(root)
            results, tested = self.summarize(modern, root/"result", "Altai.2013")
            self.assertEqual(set(results), {"family:ERV1", "family:Altai.2013_like", "1kg::archaic:Altai.2013:any"})
            self.assertEqual(len(tested), 3)

    def test_namespace_validation(self):
        self.assertEqual(hm.archaic_label_reference("ukb009::archaic:Altai.2013:carrier_ge_0.1"), "Altai.2013")
        for label in ("family:Altai.2013", "other:archaic:Altai:any", "ERVK"):
            self.assertIsNone(hm.archaic_label_reference(label))
        for label in ("::archaic:Altai:any", "1kg::archaic::any", "1kg::archaic:Altai:",
                      "1kg::archaic:Altai", "1kg::archaic:Altai:any::archaic:Denisova:any"):
            with self.assertRaises(ValueError):
                hm.archaic_label_reference(label)

    def test_cli_default_and_explicit_selection_forward_to_summary(self):
        base = ["summary", "--modern-dir", "source", "--out", "derived"]
        self.assertEqual(hm.parser().parse_args(base).archaic_references, ",".join(scan.DEFAULT_ARCHAIC_REFERENCES))
        with patch.object(hm, "summarize_common_support", return_value={}) as summary:
            hm.main(base+["--archaic-references", "Altai,Altai.2013"])
        self.assertEqual(summary.call_args.args[-1], "Altai,Altai.2013")


if __name__ == "__main__":
    unittest.main()
