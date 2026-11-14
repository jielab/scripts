"""Scientific data boundary tests: exact dates, provenance and no future input."""
from __future__ import annotations

import copy
import csv
import gzip
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import resolve_followup
from data import SCHEMA, Importer, CohortDataset, collate_batch, day, field_parts, prepare
from demo import make_demo


def write_tsv(path, names, records):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(names); writer.writerows(records)


class DataBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.raw = self.root / "input"; self.raw.mkdir()
        self.cfg = {"seed": 19, "paths": {"ukb": str(self.raw), "out": str(self.root / "output")},
                    "data": {"admin_end": "2025-01-01", "use_input_split": True,
                             "min_token_count": 1, "min_feature_count": 1,
                             "landmark": "baseline", "min_clinical_days": 0,
                             "endpoints": {"cad": ["I25"], "t2d": ["E11"], "death": ["DEATH"]},
                             "modalities": []}}
        fu = resolve_followup()
        fu['date'] = self.cfg['data']['admin_end']
        self.cfg['followup_provenance'] = fu
        mocked = patch('data.resolve_followup', return_value=fu)
        mocked.start(); self.addCleanup(mocked.stop)
        write_tsv(self.raw / "participants.tsv",
                  ["eid", "birth_date", "baseline_date", "exit_date", "death_date", "sex", "split"],
                  [[str(i), "1950-01-15", "2010-01-01", "2025-01-01", "", i % 2, split]
                   for i, split in [(1, "train"), (2, "validation"), (3, "test")]])

    def tearDown(self):
        self.temporary.cleanup()

    def canonical(self):
        write_tsv(self.raw / "events.tsv", ["eid", "date", "code"],
                  [[i, date, code] for i in (1, 2, 3) for date, code in [
                      ("2000-01-01", "I25.1"), ("2000-01-01", "K21"),
                      ("2005-01-01", "J06"), ("2012-01-01", "E11.9"),
                      ("2018-01-01", "N18")]])
        write_tsv(self.raw / "observations.tsv", ["eid", "date", "modality", "feature", "value"], [
            [1, "2011-01-01", "clinical", "x", 12],
            [2, "2009-01-01", "clinical", "x", 1000],
            [3, "2009-01-01", "clinical", "x", -5000],
            [1, "2009-01-01", "clinical", "x", 10],
            [2, "2009-01-01", "clinical", "validation_only", 30]])
        return prepare(self.cfg)

    def database(self):
        return sqlite3.connect(self.root / "output/data/cohort.sqlite")

    def test_shared_cap_keeps_earlier_loss_and_coverage(self):
        write_tsv(self.raw / "participants.tsv",
                  ["eid", "birth_date", "baseline_date", "exit_date", "death_date", "date_lost", "split"], [
                      [1, "1950-01-15", "2010-01-01", "2030-01-01", "", "", "train"],
                      [2, "1950-01-15", "2010-01-01", "2030-01-01", "", "2020-01-01", "validation"],
                      [3, "1950-01-15", "2010-01-01", "2030-01-01", "2021-01-01", "", "test"],
                      [4, "1950-01-15", "2010-01-01", "2019-01-01", "", "", "test"]])
        db = sqlite3.connect(":memory:"); db.executescript(SCHEMA)
        self.addCleanup(db.close)
        importer = Importer(db, self.cfg, self.raw)
        importer.participants()
        self.assertEqual([importer.people[str(i)][2] for i in range(1, 5)],
                         [day("2025-01-01"), day("2020-01-01"), day("2025-01-01"), day("2019-01-01")])
        # Death remains a separate event, rather than pre-censoring the death endpoint.
        self.assertEqual(importer.people['3'][3], day("2021-01-01"))

    def test_dates_and_actual_ukb_field_headers(self):
        self.assertIsNone(day("1900-01-01"))
        self.assertIsNone(day("2037-07-07"))
        self.assertIsNone(day("2024-02-30"))
        self.assertIsNone(day("01/02/2024"))
        self.assertEqual(day("2024-02-29T12:30:00"), day("2024-02-29"))
        self.assertEqual(field_parts("p41280_a11"), ("p41280", 0, 11))
        self.assertEqual(field_parts("p53_i2"), ("p53", 2, None))
        self.assertEqual(field_parts("p4080_i0_a1"), ("p4080", 0, 1))
        self.assertEqual(field_parts("53-2.0"), ("p53", 2, 0))
        self.assertEqual(field_parts("f.53.2.0"), ("p53", 2, 0))

    def test_training_only_schema_and_scaler(self):
        meta = self.canonical()
        self.assertEqual(meta["features"], ["clinical:x"])
        self.assertAlmostEqual(meta["numeric_scaler"]["clinical:x"]["mean"], 11)
        self.assertAlmostEqual(meta["numeric_scaler"]["clinical:x"]["std"], np.sqrt(2))
        self.assertEqual(meta["numeric_scaler"]["clinical:x"]["n"], 2)
        self.assertEqual(meta["preprocessing_fit_split"], "train")

    def test_risk_prefix_hides_future_and_retains_prevalence(self):
        meta = self.canonical()
        ds = CohortDataset(self.cfg, "train", meta=meta, mode="risk", max_nodes=1)
        item = ds[0]
        self.assertGreater(len(item["times"]), 1)  # Risk history is never windowed.
        self.assertTrue(np.all(item["times"] <= item["query_age"]))
        self.assertGreater(item["query_age"], item["times"][-1])
        self.assertEqual(item["numeric_mask_flat"].sum(), 1)
        self.assertAlmostEqual(float(item["numeric_flat"][item["numeric_mask_flat"]][0]), -1 / np.sqrt(2), places=5)
        self.assertFalse(item["label_mask"].any())
        self.assertFalse(item["next_event"].any())
        self.assertFalse(item["risk_mask"][0])  # CAD before all sequence windows.
        self.assertEqual(item["labels"].tolist(), [0, 1, 0])
        self.assertTrue(item["status"][1])
        self.assertGreater(item["event_times"][1], 0)

    def test_strictly_later_targets_and_reproducible_windows(self):
        self.canonical()
        full = CohortDataset(self.cfg, "train", max_nodes=100)[0]
        self.assertTrue(np.all(full["next_dt"][full["label_mask"]] > 0))
        # The two same-day codes predict the next dated encounter, not each other.
        repeated = np.flatnonzero(np.diff(full["times"]) == 0)
        self.assertEqual(len(repeated), 1)
        i = repeated[0]
        self.assertEqual(full["next_event"][i], full["next_event"][i + 1])
        windows = CohortDataset(self.cfg, "train", max_nodes=2)
        windows.set_epoch(7)
        np.testing.assert_array_equal(windows[0]["times"], windows[0]["times"])
        self.assertEqual(len(windows[0]["times"]), 2)

    def test_midpoint_uses_distinct_days_and_contains_whole_day(self):
        self.canonical()
        cfg = copy.deepcopy(self.cfg)
        cfg["data"].update(landmark="midpoint", min_clinical_days=3)
        item = CohortDataset(cfg, "train", mode="risk")[0]
        self.assertEqual(item["landmark_day"], day("2012-01-01"))
        self.assertFalse(item["risk_mask"][1])  # On-landmark disease is prevalent.
        with self.database() as db:
            db.execute("DELETE FROM events WHERE eid='2' AND day>?", (day("2005-01-01"),))
        self.assertEqual(len(CohortDataset(cfg, "validation", mode="risk")), 0)

    def test_exact_hospital_episode_join_and_repeated_encounters(self):
        write_tsv(self.raw / "rap/raw/hesin.tsv.gz", ["eid", "ins_index", "epistart", "admidate"], [
            [1, 1, "2000-01-01", "1999-12-30"], [2, 1, "2001-01-01", ""],
            [3, 8, "", "2002-01-01"], [1, 2, "2003-01-01", ""],
            [1, 3, "1900-01-01", ""]])
        write_tsv(self.raw / "rap/raw/hesin_diag.tsv.gz", ["eid", "ins_index", "diag_icd10"], [
            [2, 1, "E11.9"], [1, 1, "I25.1"], [1, 2, "I25.2"],
            [3, 8, "I50"], [1, 999, "N18"], [1, 3, "G30"]])
        meta = prepare(self.cfg)
        self.assertEqual(meta["event_source_mode"], "hospital_episode_records")
        with self.database() as db:
            events = db.execute("SELECT * FROM events ORDER BY eid,day").fetchall()
        self.assertEqual(events, [("1", day("2000-01-01"), "I25"), ("1", day("2003-01-01"), "I25"),
                                  ("2", day("2001-01-01"), "E11"), ("3", day("2002-01-01"), "I50")])
        self.assertEqual(meta["summary"]["diagnosis_unmatched_episode"], 2)

    def test_first_occurrence_pairing_is_by_eid_and_array_index(self):
        write_tsv(self.raw / "rap/icd10.code.tab.gz", ["eid", "p41270"], [
            [1, "I25.1|E11.9"], [2, "N18|I50"], [3, "J06|K21"]])
        write_tsv(self.raw / "rap/icd10.date.tab.gz", ["eid", "p41280_a1", "p41280_a0"], [
            [3, "2013-01-01", "2003-01-01"], [1, "2011-01-01", "2001-01-01"],
            [2, "2012-01-01", "2002-01-01"]])
        self.cfg["data"]["allow_first_occurrence"] = True
        meta = prepare(self.cfg)
        self.assertEqual(meta["event_source_mode"], "first_occurrence_only")
        self.assertTrue(any("first occurrences" in s for s in meta["warnings"]))
        with self.database() as db:
            events = db.execute("SELECT day,code FROM events WHERE eid='1' ORDER BY day").fetchall()
        self.assertEqual(events, [(day("2001-01-01"), "I25"), (day("2011-01-01"), "E11")])

    def test_raw_npx_negatives_and_imaging_instance_dates(self):
        write_tsv(self.raw / "events.tsv", ["eid", "date", "code"], [[i, "2001-01-01", "I25"] for i in (1, 2, 3)])
        write_tsv(self.raw / "rap/raw/prot.tab.gz", ["eid", "GDF15"], [[1, -2.5], [2, -1.2], [3, .4]])
        write_tsv(self.raw / "rap/img.tab.gz", ["eid", "p53_i2", "p25000_i2"], [
            [1, "2016-01-01", 1440], [2, "2017-01-01", 1500], [3, "1900-01-01", 1400]])
        self.cfg["data"].update(modalities=["prot", "img"], prot_date="baseline")
        meta = prepare(self.cfg)
        with self.database() as db:
            records = db.execute("SELECT eid,day,modality,values_json FROM observations ORDER BY eid,day").fetchall()
        self.assertEqual(json.loads(records[0][3])["GDF15"], -2.5)
        self.assertEqual(records[0][1], day("2010-01-01"))
        image = [r for r in records if r[2] == "img"]
        self.assertEqual(len(image), 2)
        self.assertEqual(image[0][1], day("2016-01-01"))
        self.assertTrue(any("explicitly assigned" in s for s in meta["warnings"]))
        item = CohortDataset(self.cfg, "train", mode="risk")[0]
        img_index = meta["features"].index("img:p25000")
        self.assertFalse(item["numeric_mask_flat"][:, img_index].any())

    def test_protein_collection_date_requires_explicit_assumption(self):
        write_tsv(self.raw / "events.tsv", ["eid", "date", "code"], [[i, "2001-01-01", "I25"] for i in (1, 2, 3)])
        write_tsv(self.raw / "rap/raw/prot.tab.gz", ["eid", "GDF15"], [[1, -2.5], [2, -1.2], [3, .4]])
        self.cfg["data"]["modalities"] = ["prot"]
        with self.assertRaisesRegex(ValueError, "Undated wide protein"):
            prepare(self.cfg)
        self.assertFalse((self.root / "output/data/cohort.sqlite").exists())

    def test_raw_vip_birth_proxy_and_repeated_assessment_dates(self):
        # Exact names used in the user's rap/vip.tab.gz, not processed Rdata.
        write_tsv(self.raw / "rap/vip.tab.gz",
                  ["eid", "p34", "p52", "p53_i0", "p53_i2", "p31", "p21001_i0", "p21001_i2"],
                  [[i, 1950, 7, "2010-01-01", "2016-01-01", i % 2, 24 + i / 10, 25 + i / 10]
                   for i in range(1, 61)])
        write_tsv(self.raw / "events.tsv", ["eid", "date", "code"],
                  [[i, "2004-01-01", "J06"] for i in range(1, 61)])
        self.cfg["data"].update(modalities=["clinical"], sources={"participants": None})
        meta = prepare(self.cfg)
        self.assertEqual(meta["summary"]["birth_midmonth_proxy"], 60)
        self.assertEqual(meta["features"], ["clinical:p21001"])
        self.assertTrue(all(meta["split_counts"].get(key, 0) > 0 for key in ("train", "validation", "test")))
        with self.database() as db:
            birth = db.execute("SELECT birth_day FROM participants WHERE eid='1'").fetchone()[0]
            dates = [r[0] for r in db.execute("SELECT day FROM observations WHERE eid='1' ORDER BY day")]
        self.assertEqual(birth, day("1950-07-15"))
        self.assertEqual(dates, [day("2010-01-01"), day("2016-01-01")])

    def test_missingness_threshold_uses_training_people_with_modality(self):
        self.canonical()
        with (self.raw / "participants.tsv").open("a") as handle:
            handle.write("4\t1950-01-15\t2010-01-01\t2025-01-01\t\t1\ttrain\n")
        with (self.raw / "events.tsv").open("a") as handle:
            handle.write("4\t2004-01-01\tJ06\n")
        with (self.raw / "observations.tsv").open("a") as handle:
            handle.write("4\t2009-01-01\tclinical\tx\t14\n")
            handle.write("1\t2009-01-01\tclinical\tsparse\t30\n")
            handle.write("2\t2009-01-01\tclinical\tsparse\t3000\n")
        self.cfg["data"]["max_missing_fraction"] = .25
        meta = prepare(self.cfg)
        self.assertEqual(meta["features"], ["clinical:x"])
        self.assertEqual(meta["modality_people_train"], {"clinical": 2})
        self.assertEqual(meta["feature_people_train"]["clinical:sparse"], 1)

    def test_pca_table_does_not_promote_numeric_metadata_to_context(self):
        self.canonical()
        write_tsv(self.raw / "rap/pca.tab.gz",
                  ["eid", "PC1", "p22009_a1", "batch_id", "assessment_day", "p21003_i2"],
                  [[i, .1 * i, -.2 * i, 55, 19000, 70] for i in (1, 2, 3)])
        meta = prepare(self.cfg)
        self.assertEqual(meta["context_features"], ["sex", "PC1", "PC2"])
        with self.database() as db:
            context = json.loads(db.execute("SELECT context FROM participants WHERE eid='1'").fetchone()[0])
        self.assertEqual(set(context), {"sex", "PC1", "PC2"})

    def test_verified_field_dictionary_and_explicit_raw_feature_filters(self):
        write_tsv(self.raw / "events.tsv", ["eid", "date", "code"], [[i, "2001-01-01", "I25"] for i in (1, 2, 3)])
        write_tsv(self.raw / "rap/bbc.tab.gz",
                  ["eid", "p30600_i0", "p30601_i0", "p30602_i0", "p30700_i0", "p30800_i0"],
                  [[i, 40, 999, 19000, 100, 25] for i in (1, 2, 3)])
        (self.raw / "common").mkdir()
        (self.raw / "common/bbc.lst").write_text(
            "p30600\tbb_albumin\np30601\tassay_qc\np30602\tdate_bb_albumin\n"
            "p30700\tbb_creatinine\np30600/p30700\tderived_ratio\n")
        self.cfg["data"].update(modalities=["bbc"], feature_exclude={"bbc": ["p30700"]})
        meta = prepare(self.cfg)
        self.assertEqual(meta["features"], ["bbc:p30600"])
        self.assertTrue(any(s["name"] == "bbc_feature_dictionary" for s in meta["sources"]))
        # Explicit, exact raw-column include overrides the local dictionary.
        self.cfg["data"]["feature_include"] = {"bbc": ["p30800_i0"]}
        meta = prepare(self.cfg)
        self.assertEqual(meta["features"], ["bbc:p30800"])
        self.assertFalse(any(s["name"] == "bbc_feature_dictionary" for s in meta["sources"]))

    def test_death_absorption_and_competing_followup(self):
        write_tsv(self.raw / "participants.tsv",
                  ["eid", "birth_date", "baseline_date", "exit_date", "death_date", "sex", "split"], [
                      [1, "1950-01-15", "2010-01-01", "2025-01-01", "2013-01-01", 1, "train"],
                      [2, "1950-01-15", "2010-01-01", "2025-01-01", "", 0, "validation"],
                      [3, "1950-01-15", "2010-01-01", "2025-01-01", "", 1, "test"]])
        self.canonical()
        with self.database() as db:
            records = db.execute("SELECT day,code FROM events WHERE eid='1' ORDER BY day").fetchall()
        self.assertEqual(records[-1], (day("2013-01-01"), "DEATH"))
        self.assertNotIn("N18", [code for _, code in records])
        item = CohortDataset(self.cfg, "train", mode="risk")[0]
        self.assertEqual(item["labels"][-1], 1)
        self.assertLess(item["death_time"], item["exit_time"])

    def test_demo_and_model_batch_contract(self):
        cfg = make_demo(self.cfg, n_people=60)
        inputs = Path(cfg["paths"]["ukb"])
        mtimes = {name: (inputs / name).stat().st_mtime_ns for name in ("participants.tsv", "events.tsv", "observations.tsv")}
        again = make_demo(self.cfg, n_people=60)
        self.assertEqual(again, cfg)
        self.assertEqual(mtimes, {name: (inputs / name).stat().st_mtime_ns for name in mtimes})
        meta = prepare(cfg)
        dataset = CohortDataset(cfg, "validation", mode="risk")
        items = [dataset[0], dataset[1]]
        batch = collate_batch(items, meta)
        self.assertEqual(len(dataset), 9)
        self.assertEqual(batch["times"].shape, batch["valid"].shape)
        self.assertTrue(all(v["values"].shape == v["mask"].shape for v in batch["numeric"].values()))
        self.assertEqual(batch["labels"].shape, (2, 4))
        for row in batch["valid"]:
            self.assertFalse(np.any(np.diff(row.numpy().astype(int)) == 1))
        self.assertTrue(np.all(batch["times"].numpy() <= batch["query_age"].numpy()[:, None]))


if __name__ == "__main__":
    unittest.main()
