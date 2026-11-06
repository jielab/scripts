#!/usr/bin/env python3
"""Small, answerable regression cases for the aggregate-only date report.

Run directly with Python 3; the report tests need neither R nor UKB data.
The fixtures deliberately distinguish missing dates, definition-only matches,
duplicate evidence, unavailable observation history, and operation context.
"""
from __future__ import annotations

import base64
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "f" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


report = load_module("aggregate_date_report", "date_check_report.py")
input_io = load_module("aggregate_date_io", "date_check_io.py")

PF = ("eid", "min1", "min2", "ref", "date_assess", "date_birth", "date_death", "date_lost")
EF = ("eid", "source", "code", "date", "date_raw", "date_status", "origin",
      "granularity", "match_local", "match_official")
META = {"Y": "cvd_cad", "as_of": "2026-01-01"}


def patient(eid, min1="", min2="", ref="", baseline="2010-01-01",
            birth="1970-01-01", death="", lost=""):
    return [eid, min1, min2, ref, baseline, birth, death, lost]


def event(eid, day, source="icd10", granularity="episode", local=True,
          official=False, code="I210"):
    return [eid, source, code, day, day, "valid" if day else "unknown",
            "synthetic test row", granularity, str(local).upper(), str(official).upper()]


def summary(people, events=(), metadata=None):
    return report.make_summary(PF, people, EF, iter(events), {**META, **(metadata or {})})


def series_of(result, source="clinical"):
    return next(row for row in result["trajectory"]["series"] if row["key"] == source)


def total_retained(series):
    return sum(series["retained_days"]) + series["underflow_days"] + series["overflow_days"]


def write_tsv(path, fields, rows):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


class DenominatorTests(unittest.TestCase):
    def setUp(self):
        day = "2010-01-01"
        self.people = [
            patient("PRIVATE_ID_A", day, day, day),
            patient("PRIVATE_ID_B", "2010-01-02", day, day),
            patient("PRIVATE_ID_C", day, "2010-01-03", day),
            patient("PRIVATE_ID_D", "2010-01-04", "2010-01-05", day),
            patient("PRIVATE_ID_E", day, "", ""),
            patient("PRIVATE_ID_F", "", day, day),
            patient("PRIVATE_ID_G", "", day, ""),
            patient("PRIVATE_ID_H", "1900-01-01", "", "1901-01-01"),
        ]
        self.result = summary(self.people)

    def test_four_way_counts_partition_every_person(self):
        first = self.result["coverage"]["min1"]
        self.assertEqual((first["both"], first["a_only"], first["b_only"], first["neither_date"]),
                         (4, 1, 1, 2))
        second = self.result["coverage"]["min2"]
        self.assertEqual((second["both"], second["a_only"], second["b_only"], second["neither_date"]),
                         (5, 1, 0, 2))
        for metrics in (first, second):
            self.assertEqual(metrics["total"], 8)
            self.assertEqual(sum(row["n"] for row in metrics["four_way"]), 8)
            self.assertEqual(sum(row["n"] for row in metrics["excluded_reasons"]), metrics["excluded"])
            self.assertEqual(metrics["both"] + metrics["excluded"], 8)

    def test_common_people_prevent_a_false_improvement(self):
        # Crude comparisons look like 50% -> 60%. In the same four people,
        # one disagreement is corrected and one agreement is lost: no gain.
        self.assertEqual(self.result["pairs"]["min1"]["ref"]["same_day_pct"], 50)
        self.assertEqual(self.result["pairs"]["min2"]["ref"]["same_day_pct"], 60)
        common = self.result["common_comparison"]
        self.assertEqual(common["n"], 4)
        self.assertEqual((common["min1_same_day"], common["min2_same_day"]), (2, 2))
        self.assertEqual((common["gained"], common["lost"], common["both_same"], common["neither_same"]),
                         (1, 1, 1, 1))

    def test_empty_date_comparison_is_unknown_not_zero_percent(self):
        result = summary([patient("PRIVATE_ID_NO_DATES", baseline="")])
        self.assertEqual(result["coverage"]["min1"]["both"], 0)
        self.assertIsNone(result["pairs"]["min1"]["ref"]["same_day_pct"])
        self.assertIsNone(result["pairs"]["min1"]["ref"]["median_days"])
        self.assertEqual(result["trajectory"]["n_baseline"], 0)


class TrajectoryTests(unittest.TestCase):
    def test_duplicate_evidence_and_different_code_systems_do_not_multiply_days(self):
        eid = "PRIVATE_ID_DUPLICATE"
        people = [patient(eid)]
        events = [
            event(eid, "2008-01-01", granularity="fod_only"),
            event(eid, "2008-01-01", granularity="cached_aod"),
            event(eid, "2008-01-01", granularity="summary_first_per_code"),
            event(eid, "2008-01-01"),
            event(eid, "2008-01-01", source="gp2", granularity="record"),
            event(eid, "2008-01-01", source="gp3", granularity="record"),
            event(eid, "2011-01-01", granularity="cached_aod"),
            event(eid, "2011-01-01"),
            event(eid, "2012-01-01", source="gp2", granularity="record"),
            event(eid, "2012-01-01", source="gp3", granularity="record"),
            event(eid, "1980-01-01", source="ref", granularity="official_ref", local=False),
            event(eid, "1995-01-01", local=False, official=True),
            event(eid, "1990-01-01", source="srd", granularity="self_report_interpolated"),
            event(eid, "2013-01-01", source="opcs4_context", granularity="record", code="K751"),
        ]
        result = summary(people, events, {"procedure_context_defined": True})
        clinical = series_of(result)
        self.assertEqual(total_retained(clinical), 3)
        self.assertEqual(result["trajectory"]["n_clinical_people"], 1)
        self.assertEqual(result["trajectory"]["persistence"]["n_repeated"], 1)
        self.assertEqual(result["procedures"]["n_people"], 1)
        self.assertEqual(result["n_patients"], 1)
        self.assertEqual(total_retained(series_of(result, "hospital")), 2)
        self.assertEqual(total_retained(series_of(result, "gp")), 2)

    def test_first_date_only_is_identified_as_limited_history(self):
        eid = "PRIVATE_ID_FOD_ONLY"
        result = summary([patient(eid)], [event(eid, "2011-01-01", granularity="fod_only")])
        provenance = result["trajectory"]["provenance"]
        self.assertGreaterEqual(sum(row["fod_fallback_only"] for row in provenance), 1)
        self.assertEqual(result["trajectory"]["n_raw_clinical_people"], 0)
        self.assertEqual(result["trajectory"]["persistence"]["n_repeated"], 0)
        self.assertEqual(sum(series_of(result)["raw_record_days"]), 0)

    def test_aod_in_one_hospital_code_system_does_not_erase_other_system_fod(self):
        eid = "PRIVATE_ID_MIXED_HISTORY"
        events = [event(eid, "1995-01-01", source="icd9", granularity="fod_only", code="410"),
                  event(eid, "2011-01-01", source="icd10", granularity="cached_aod")]
        result = summary([patient(eid)], events)
        self.assertEqual(total_retained(series_of(result, "hospital")), 2)
        self.assertEqual(total_retained(series_of(result)), 2)

    def test_unavailable_gp_does_not_become_a_zero_percent_patient_curve(self):
        result = summary([patient("PRIVATE_ID_NO_GP")], metadata={"raw_gp2_available": False,
                                                               "raw_gp3_available": False})
        gp = next((row for row in result["trajectory"]["series"] if row["key"] == "gp"), None)
        if gp is not None:
            self.assertEqual(gp["n_people"], 0)
            weights = gp["patient_weighted_pct"]
            self.assertTrue(weights is None or all(value is None for value in weights))

    def test_missing_baseline_excludes_alignment_but_keeps_disease_evidence(self):
        eid = "PRIVATE_ID_BASELINE_MISSING"
        result = summary([patient(eid, baseline="")], [event(eid, "2011-01-01")])
        self.assertEqual(result["n_patients"], 1)
        self.assertEqual(result["trajectory"]["n_baseline"], 0)
        self.assertEqual(result["trajectory"]["n_clinical_people"], 1)
        self.assertEqual(result["trajectory"]["n_baseline_clinical"], 0)
        self.assertEqual(total_retained(series_of(result)), 0)

    def test_implausible_baselines_do_not_create_aligned_events(self):
        people = [
            patient("PRIVATE_ID_BASE_BEFORE_BIRTH", baseline="1960-01-01"),
            patient("PRIVATE_ID_BASE_AFTER_DEATH", baseline="2020-01-01", death="2019-01-01"),
            patient("PRIVATE_ID_BASE_AFTER_LOSS", baseline="2020-01-01", lost="2019-01-01"),
            patient("PRIVATE_ID_BASE_FUTURE", baseline="2027-01-01"),
        ]
        result = summary(people, [event(row[0], "2011-01-01") for row in people])
        self.assertEqual(result["trajectory"]["n_baseline"], 0)
        self.assertEqual(result["quality"]["invalid_baseline"], 4)
        self.assertEqual(total_retained(series_of(result)), 0)

    def test_post_death_or_loss_and_future_dates_are_not_trajectory_days(self):
        people = [patient("PRIVATE_ID_END", death="2012-01-01"),
                  patient("PRIVATE_ID_LOST", lost="2012-01-01"),
                  patient("PRIVATE_ID_FUTURE")]
        events = [event("PRIVATE_ID_END", "2011-01-01"), event("PRIVATE_ID_END", "2013-01-01"),
                  event("PRIVATE_ID_LOST", "2013-01-01"), event("PRIVATE_ID_FUTURE", "2027-01-01")]
        result = summary(people, events)
        self.assertEqual(total_retained(series_of(result)), 1)
        quality = result["quality"]
        self.assertEqual(quality["after_known_end_rows"] + quality["invalid_or_undated_event_rows"], 3)

    def test_special_invalid_and_prebirth_dates_are_not_real_events(self):
        eid = "PRIVATE_ID_BAD_DATE"
        dates = ["1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03",
                 "1909-09-09", "2037-07-07", "2011-02-30", "1969-01-01", ""]
        result = summary([patient(eid)], [event(eid, day) for day in dates])
        self.assertEqual(result["trajectory"]["n_clinical_people"], 0)
        self.assertEqual(total_retained(series_of(result)), 0)

    def test_reference_and_procedures_alone_do_not_become_clinical_activity(self):
        eid = "PRIVATE_ID_CONTEXT"
        events = [event(eid, "2011-01-01", source="ref", granularity="official_ref", local=False),
                  event(eid, "2012-01-01", local=False, official=True),
                  event(eid, "2013-01-01", source="opcs4_context", granularity="record", code="K401")]
        result = summary([patient(eid)], events, {"procedure_context_defined": True})
        self.assertEqual(result["n_patients"], 1)
        self.assertEqual(result["trajectory"]["n_clinical_people"], 0)
        self.assertEqual(result["procedures"]["n_people"], 1)
        self.assertEqual(result["coverage"]["min1"]["both"], 0)

    def test_unknown_context_participant_cannot_silently_expand_population(self):
        with self.assertRaises(ValueError):
            summary([patient("PRIVATE_ID_ALREADY_INCLUDED")],
                    [event("PRIVATE_ID_CONTEXT_ONLY", "2011-01-01", source="opcs4_context", code="K401")])


class PrivacyBuildTests(unittest.TestCase):
    def test_real_build_contains_only_aggregates_including_embedded_payload(self):
        eid = "PRIVATE_PARTICIPANT_EVIDENCE_1000102_DO_NOT_EXPORT"
        private_meta = "PRIVATE_HOME_AND_PERSON_872341"
        people = [patient(eid, "2011-01-01", "2011-01-01", "2011-01-01")]
        events = [event(eid, "2011-01-01")]
        metadata = {**META, "debug_eid": eid, "r_inputs": [{"path": private_meta}],
                    "unexpected_patient_rows": [people], "is_demo": True}
        with tempfile.TemporaryDirectory() as td:
            folder = Path(td)
            write_tsv(folder / "patients.tsv.gz", PF, people)
            write_tsv(folder / "events.tsv.gz", EF, events)
            (folder / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
            output = folder / "report.html"
            result = report.build(folder, output, "cvd_cad", folder / "metadata.json")
            document = output.read_text(encoding="utf-8")
            written = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(written["n_patients"], 1)
            self.assertEqual(result["n_patients"], 1)
            texts = [document, json.dumps(result), json.dumps(written)]
            # Older templates used a gzip/base64 payload. Inspect any such
            # payload too, so hiding rows under compression cannot pass.
            for token in re.findall(r"[A-Za-z0-9+/]{80,}={0,2}", document):
                try:
                    raw = base64.b64decode(token, validate=True)
                    if raw.startswith(b"\x1f\x8b"):
                        texts.append(gzip.decompress(raw).decode("utf-8"))
                except (ValueError, OSError, UnicodeError):
                    pass
            for text in texts:
                self.assertNotIn(eid, text)
                self.assertNotIn(private_meta, text)
                self.assertNotIn('"patient_fields"', text)
                self.assertNotIn('"event_fields"', text)
                self.assertNotIn('"patients":', text)
                self.assertNotIn('"events":', text)


class ProcedureInputTests(unittest.TestCase):
    def make_definitions(self, folder):
        common = folder / "common"
        common.mkdir()
        (common / "icd10.lst").write_text("I2[0-5]\tcvd_cad\n", encoding="utf-8")
        (common / "opcs4.lst").write_text("K4[0-6]|k49|K50|K7[5-6]\tmi\n", encoding="utf-8")
        (common / "fod.lst").write_text("p131296|p131298\tcvd_cad\n", encoding="utf-8")
        (common / "fod.each.lst").write_text("p131296\tI20\np131298\tI21\n", encoding="utf-8")
        mapping = folder / "maps.tsv"
        mapping.write_text("source\tcode\ticd10\nsrd\t1074\tI20\nsrd\t1075\tI21\n", encoding="utf-8")
        return input_io.Definitions("cvd_cad", common, mapping)

    def test_cad_context_excludes_conduit_and_diagnostic_procedures(self):
        with tempfile.TemporaryDirectory() as td:
            folder = Path(td)
            defs = self.make_definitions(folder)
            codes = ["K40.1", "K499", "K50.4", "K75.1", "k45.1", "K761", "K519", "K75"]
            operations = folder / "operations.tsv"
            rows = [[f"PRIVATE_OP_{i}", code, "" if i in (4, 7) else "2011-01-01", "2011-01-02"]
                    for i, code in enumerate(codes)]
            write_tsv(operations, ("eid", "oper4", "opdate", "admidate"), rows)
            buffer, info = io.StringIO(), {}
            input_io.stream_operations(operations, defs, csv.writer(buffer, delimiter="\t"), info)
            output = list(csv.reader(io.StringIO(buffer.getvalue()), delimiter="\t"))
            self.assertEqual(len(output), 6)
            self.assertEqual({row[1] for row in output}, {"opcs4_context"})
            self.assertEqual({row[2] for row in output}, {"K40.1", "K499", "K50.4", "K75.1", "k45.1", "K75"})
            self.assertEqual(sum(row[3] == "" for row in output), 2)
            self.assertEqual(info["matched_records"], 6)
            self.assertFalse(defs.local_defined("opcs4"))

    def test_absent_operation_date_is_not_replaced_by_admission_date(self):
        with tempfile.TemporaryDirectory() as td:
            folder = Path(td)
            defs = self.make_definitions(folder)
            operations = folder / "operations.tsv"
            write_tsv(operations, ("eid", "oper4", "admidate"),
                      [["PRIVATE_OP_DATE_UNKNOWN", "K401", "2010-01-01"]])
            buffer, info = io.StringIO(), {}
            input_io.stream_operations(operations, defs, csv.writer(buffer, delimiter="\t"), info)
            row = next(csv.reader(io.StringIO(buffer.getvalue()), delimiter="\t"))
            self.assertEqual(row[3], "")
            self.assertIsNone(info["date_column"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
