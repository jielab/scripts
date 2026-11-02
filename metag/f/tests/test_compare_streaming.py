"""Compare real gzip assignment tables using one bounded streaming pass per case."""
from collections import Counter
import csv
import gzip
import itertools
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_compare as compare
import hpc_query


class StreamingComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def records(self, count=6):
        sample = dict(sample_id="dataset/sample", dataset="dataset", sample="sample", r1="unused", r2="unused2")
        profiles = [{"name": f"profile{i}", "classifier_engine": "kraken", "parameters": {}}
                    for i in range(count)]
        cases = compare.make_cases([sample], profiles=profiles, repeats=1, limit=5)
        patterns = [(11, "higher", 0, 22, 11), (11, 33, 0, "higher", 22),
                    (0, 33, 44, 22, 11)]
        records = []
        for i, case in enumerate(cases):
            folder = self.root / case["case_id"]
            output = folder / "output"
            rows = []
            for ordinal, value in enumerate(patterns[i % len(patterns)], 1):
                higher = value == "higher"
                species = 0 if higher else value
                rows.append({"ordinal": ordinal, "r1_read_id": f"read{ordinal}/1",
                             "r2_read_id": f"read{ordinal}/2", "assigned_taxid": 10 if higher else species,
                             "assigned_rank": "genus" if higher else "species" if species else "",
                             "species_taxid": species,
                             "status": "higher_rank" if higher else "assigned_species" if species else "unclassified"})
            record = {"schema": compare.SCHEMA, "identity": f"identity-{i}", "case": case,
                      "status": "complete", "source_database_identity": "source",
                      "output": str(output), "output_identity": f"output-{i}",
                      "pipeline_wall_seconds": 1.0, "timing_eligible": True,
                      "summary": {"fragments": len(rows), "reads": 2 * len(rows), "paired": True,
                                  "input_digest": "identical-input", "count_conservation": True,
                                  "statuses": dict(Counter(row["status"] for row in rows)),
                                  "classification_seconds": 0.4, "split_seconds": 0.6}}
            self.write_rows(record, rows)
            records.append(record)
        return records

    @staticmethod
    def write_rows(record, rows):
        path = Path(record["output"]) / "species" / "read-assignments.tsv.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        fields = ["ordinal", "r1_read_id", "r2_read_id", "assigned_taxid", "assigned_rank", "species_taxid", "status"]
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def read_rows(record):
        with compare._assignment_reader(record["output"]) as reader:
            return list(reader)

    def summarize(self, records):
        receipts = []
        for record in records:
            receipt = Path(record["output"]).parent / "case.json"
            receipt.write_text(json.dumps(record))
            receipts.append(receipt)
        # Publication marker verification is independently tested by the query
        # tests. Here the real assignment files and all comparison checks run.
        with patch.object(hpc_query, "_complete", return_value=True):
            return compare.summarize_cases(receipts, self.root / "report")

    def test_all_pair_counts_match_two_case_reference_with_one_scan_per_case(self):
        records = self.records()
        expected = [compare._pair(left, right) for left, right in itertools.combinations(records, 2)]
        calls = Counter()
        original = compare._assignments

        def counted(record):
            calls[record["case"]["case_id"]] += 1
            yield from original(record)

        with patch.object(compare, "_assignments", side_effect=counted):
            actual, errors = compare._all_pairs(records)
        self.assertEqual(actual, expected)
        self.assertEqual(errors, [])
        self.assertEqual(list(calls.values()), [1] * len(records))
        # Independent hand-counted fixture: each outcome occurs once, with two
        # identical normalized assignments (one species and one unclassified).
        for key in ("same_species", "different_species", "left_species_right_unresolved",
                    "right_species_left_unresolved", "both_unresolved"):
            self.assertEqual(actual[0][key], 1)
        self.assertEqual(actual[0]["same_normalized_assignment"], 2)
        self.assertEqual(actual[0]["species_agreement_among_both_species"], 0.5)
        self.assertEqual(actual[0]["species_agreement_over_all_fragments"], 0.2)

    def test_summary_performs_two_scans_per_case_instead_of_case_count_squared(self):
        records = self.records()
        calls = Counter()
        original = compare._assignments

        def counted(record):
            calls[record["case"]["case_id"]] += 1
            yield from original(record)

        with patch.object(compare, "_assignments", side_effect=counted):
            report = self.summarize(records)
        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["pairwise_comparisons"], 15)
        self.assertEqual(list(calls.values()), [2] * len(records))
        self.assertEqual(sum(calls.values()), 12)

    def test_different_inputs_and_sources_report_errors_without_losing_valid_cohorts(self):
        records = self.records(5)
        for index in (1, 3):
            records[index]["summary"]["input_digest"] = "different-input"
        records[4]["source_database_identity"] = "different-reference"
        expected = [compare._pair(records[0], records[2]), compare._pair(records[1], records[3])]
        actual, errors = compare._all_pairs(records)
        self.assertEqual(actual, expected)
        self.assertEqual(len(errors), 8)
        self.assertTrue(all("different input" in error["error"] for error in errors))

    def test_read_id_mismatch_makes_summary_incomplete_but_keeps_unaffected_pair(self):
        records = self.records(3)
        rows = self.read_rows(records[1])
        rows[2]["r1_read_id"] = "wrong-read/1"
        self.write_rows(records[1], rows)
        report = self.summarize(records)
        self.assertEqual(report["status"], "incomplete")
        self.assertEqual(report["valid_completed_cases"], 3)
        self.assertEqual(report["pairwise_comparisons"], 1)
        self.assertEqual(len(report["errors"]), 2)
        self.assertTrue(all("matching ordinal read pairs" in error["error"] for error in report["errors"]))
        with (self.root / "report" / "agreement.tsv").open() as handle:
            agreement = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(agreement[0]["left_case"], records[0]["case"]["case_id"])
        self.assertEqual(agreement[0]["right_case"], records[2]["case"]["case_id"])

    def test_truncation_during_comparison_discards_only_affected_pairs(self):
        records = self.records(3)
        self.write_rows(records[1], self.read_rows(records[1])[:-1])
        actual, errors = compare._all_pairs(records)
        self.assertEqual(actual, [compare._pair(records[0], records[2])])
        self.assertEqual(len(errors), 2)
        self.assertTrue(all("row count" in error["error"] for error in errors))

    def test_early_pair_failure_closes_assignment_iterators(self):
        records = self.records(2)
        rows = self.read_rows(records[1])
        rows[0]["r1_read_id"] = "wrong-read/1"
        self.write_rows(records[1], rows)
        closed = set()
        original = compare._assignments

        def watched(record):
            try:
                yield from original(record)
            finally:
                closed.add(record["case"]["case_id"])

        with patch.object(compare, "_assignments", side_effect=watched):
            actual, errors = compare._all_pairs(records)
        self.assertEqual(actual, [])
        self.assertEqual(len(errors), 1)
        self.assertEqual(closed, {record["case"]["case_id"] for record in records})


if __name__ == "__main__":
    unittest.main()
