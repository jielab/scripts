"""Prepared-data identity and safe reuse across expensive UKB import stages."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare import prepare_cohort
from common import resolve_followup


class PrepareReuseTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.raw = self.root / "raw"; self.raw.mkdir()
        self.source = self.raw / "events.tsv"
        self.source.write_text("eid\tdate\tcode\n1\t2000-01-01\tI25\n")
        self.folder = self.root / "out/data"; self.folder.mkdir(parents=True)
        (self.folder / "cohort.sqlite").touch()
        self.cfg = {"seed": 2026, "paths": {"ukb": str(self.raw), "out": str(self.root / "out")},
                    "data": {"admin_end": "2023-04-01", "landmark": "midpoint", "min_clinical_days": 3}}
        self.fu = resolve_followup()
        self.fu['date'] = self.cfg['data']['admin_end']
        self.cfg['followup_provenance'] = copy.deepcopy(self.fu)
        mocked = patch('prepare.resolve_followup', side_effect=lambda: dict(self.fu))
        mocked.start(); self.addCleanup(mocked.stop)
        stat = self.source.stat()
        self.meta = {"input_root": str(self.raw.resolve()), "seed": 2026,
                     "data_config": copy.deepcopy(self.cfg["data"]),
                     "sources": [{"name": "events", "path": str(self.source),
                                  "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}]}
        self.write_meta()

    def write_meta(self):
        (self.folder / "meta.json").write_text(json.dumps(self.meta))

    def tearDown(self):
        self.temporary.cleanup()

    def test_unchanged_inputs_and_runtime_landmark_reuse_prepared_data(self):
        with patch("prepare.prepare") as rebuild:
            self.assertEqual(prepare_cohort(self.cfg), self.meta)
            changed = copy.deepcopy(self.cfg)
            changed["data"].update(landmark="age", landmark_age=60, max_nodes=32, min_clinical_days=0)
            self.assertEqual(prepare_cohort(changed), self.meta)
            rebuild.assert_not_called()

    def test_input_root_change_cannot_reuse_old_participants(self):
        changed = copy.deepcopy(self.cfg)
        changed["paths"]["ukb"] = str(self.root / "other_ukb_release")
        with patch("prepare.prepare", return_value={"rebuilt": True}) as rebuild:
            self.assertEqual(prepare_cohort(changed), {"rebuilt": True})
            rebuild.assert_called_once_with(changed)

    def test_modified_source_requires_rebuild(self):
        self.source.write_text(self.source.read_text() + "2\t2001-01-01\tE11\n")
        with patch("prepare.prepare", return_value={"rebuilt": True}) as rebuild:
            self.assertEqual(prepare_cohort(self.cfg), {"rebuilt": True})
            rebuild.assert_called_once_with(self.cfg)

    def test_changed_source_with_fitted_model_requires_explicit_replacement(self):
        model = self.root / "out/model"; model.mkdir()
        (model / "flux.pt").touch()
        self.source.write_text(self.source.read_text() + "2\t2001-01-01\tE11\n")
        with patch("prepare.prepare") as rebuild:
            with self.assertRaisesRegex(ValueError, "inputs/settings changed"):
                prepare_cohort(self.cfg)
            rebuild.assert_not_called()
        with patch("prepare.prepare", return_value={"rebuilt": True}) as rebuild:
            self.assertEqual(prepare_cohort(self.cfg, replace=True), {"rebuilt": True})
            rebuild.assert_called_once_with(self.cfg)

    def test_offline_original_source_remains_usable(self):
        self.source.unlink()
        with patch("prepare.prepare") as rebuild:
            self.assertEqual(prepare_cohort(self.cfg), self.meta)
            rebuild.assert_not_called()

    def test_source_hash_only_does_not_rebuild_date_identical_cache(self):
        self.fu['source_sha256'] = 'f' * 64
        with patch('prepare.prepare') as rebuild:
            self.assertEqual(prepare_cohort(self.cfg), self.meta)
            rebuild.assert_not_called()

    def test_seed_and_administrative_end_are_part_of_preparation_identity(self):
        for key, value in (("seed", 101), ("admin_end", "2024-01-01")):
            changed = copy.deepcopy(self.cfg)
            if key == "seed":
                changed[key] = value
            else:
                self.fu['date'] = value
                changed['followup_provenance'] = dict(self.fu)
                changed["data"][key] = self.fu['date']
            with self.subTest(key=key), patch("prepare.prepare", return_value={"rebuilt": True}) as rebuild:
                self.assertEqual(prepare_cohort(changed), {"rebuilt": True})
                rebuild.assert_called_once_with(changed)


if __name__ == "__main__":
    unittest.main()
