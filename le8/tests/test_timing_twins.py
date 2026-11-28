from pathlib import Path
import sys
import unittest
import json
import tempfile
from zipfile import ZipFile
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"f"))
from paper_focus import aggregate_exports
from dx_dates import derive
from c1_timing_twins import timing, twin_counts, neighbor_risk


class TimingTwins(unittest.TestCase):
    def test_discordant_prior_early_and_late_source(self):
        base=pd.Timestamp("2010-01-01")
        d=pd.DataFrame(dict(eid=list("abcdef"),birth_date="1950-01-01",date_attend=base,
             fod_icd10_Y=[base+pd.Timedelta(days=k) for k in (700,365,366,1600,1400,1200)],
             fod_ref_Y=[base-pd.Timedelta(days=700),pd.NaT,pd.NaT,
                        base+pd.Timedelta(days=1000),base+pd.Timedelta(days=1500),base+pd.Timedelta(days=1300)]))
        p=pd.DataFrame(dict(eid=list("abcdef"),srd=[pd.NaT]*4+[base+pd.Timedelta(days=4000),pd.NaT]))
        snapshot,_,_=derive(d,p,"Y","2023-04-01")
        broad,clean,flags,_=timing(snapshot)
        self.assertTrue(broad.all())
        np.testing.assert_array_equal(clean,[False,False,True,False,True,True])
        self.assertTrue(flags["prior_vs_ICD_at_risk"][0])
        self.assertTrue(flags["first_year_any_source"][1])
        self.assertTrue(flags["source_first_span_gt365"][3])

    def test_counts_exclude_own_family_and_empty_donors_use_background(self):
        d=pd.DataFrame(dict(eid=list("abcde"),group=["f","f","g","h","i"]))
        xy=np.array([[0.],[0.],[.5],[2.],[.25]])
        sex=np.array([0,0,0,0,1]); ok=np.ones(5,bool)
        count=twin_counts(d,xy,sex,ok,np.array([0,2]),1)
        np.testing.assert_array_equal(count,[1,1,1,0,0])
        background=np.array([.1,.3]); train=np.array([True,False,False,False,False])
        for residual in (False,True):
            out=neighbor_risk(d,xy,sex,ok,train,np.array([3,4]),np.array([1,0,0,0,0]),np.ones(5),background,1,residual)
            np.testing.assert_equal(out,background)

    def test_equal_budget_source_survives_workbook_packing(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with ZipFile(root/"c4.Fig16.equal_budget.xlsx","w") as book:
                part="le8/exports/value.bin"
                book.writestr(part,"model,delta_AUC\nYS10,0.015\n")
                book.writestr("le8/manifest.json",json.dumps({"files":{"c4.focus.contrasts.csv":{"part":part}}}))
            exports=list(aggregate_exports(root,{"c4.focus.contrasts.csv":"C4_panel_contrasts"}))
            self.assertEqual(len(exports),1)
            self.assertAlmostEqual(exports[0][1].delta_AUC.iloc[0],.015)


if __name__=="__main__":unittest.main()
