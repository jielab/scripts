"""Current-only result admission: no historical trees or automatic repairs."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"f"))
from current_release import require_current, SCHEMA
from dx_layout import scope_for_read
from outputs import default_cohort_path


class CurrentRelease(unittest.TestCase):
    def test_current_resume_and_changed_code_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/"analysis"
            identity=dict(schema=SCHEMA,source_sha256="one")
            require_current(root,create=True,identity=identity)
            (root/"results").mkdir()
            require_current(root,identity=identity)
            with self.assertRaisesRegex(ValueError,"mismatch"):
                require_current(root,identity=dict(identity,source_sha256="two"))
            self.assertEqual(json.loads((root/"RELEASE.json").read_text()),identity)

    def test_existing_unversioned_tree_is_not_adopted(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/"old.xlsx").write_text("old")
            with self.assertRaisesRegex(ValueError,"Missing current"):
                require_current(root,create=True)
            self.assertFalse((root/"RELEASE.json").exists())
            self.assertEqual((root/"old.xlsx").read_text(),"old")

    def test_readers_never_fall_back_to_old_layout(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/"icd/qc").mkdir(parents=True)
            (root/"icd/qc/dx.json").write_text("{}")
            old=root/"logs/freeze/c1_outer_roster_seed2026.tsv"
            old.parent.mkdir(parents=True);old.touch()
            self.assertEqual(scope_for_read(root,"icd"),root/"main/results/icd")
            self.assertEqual(default_cohort_path(root),root/"tmp/freeze"/old.name)

    def test_endpoint_hash_uses_only_current_endpoint_columns(self):
        import pandas as pd
        from dx_dates import derive, frame_hash, BASE_ENDPOINTS
        cohort=pd.DataFrame(dict(eid=["1","2","3"],date_attend=["2010-01-01"]*3,
            birth_date=["1960-01-01"]*3,fod_icd10_cad=["2012-01-01",None,"2009-01-01"],
            fod_ref_cad=["2012-01-01",None,None]))
        data,counts,_=derive(cohort,pd.DataFrame(dict(eid=["1","2","3"])),"cad","2023-04-01")
        self.assertFalse(any(c.startswith("ham") for c in data))
        self.assertEqual(set(counts.dx),{"icd","ref","icd1","ref1","any1"})
        columns=["eid",*[key+suffix for key in BASE_ENDPOINTS for suffix in ("","_state","_Yt2e","_t2e")]]
        self.assertEqual(set(counts.loc[counts.dx.isin(BASE_ENDPOINTS),"endpoint_sha256"]),
                         {frame_hash(data[columns])})


if __name__=="__main__":unittest.main()
