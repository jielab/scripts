from pathlib import Path
from tempfile import TemporaryDirectory
import sys
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"f"))
from cad_icd_first import analyze_icd_first
from dx_dates import derive
from paper_focus import highlight_figures


class CADFirst(unittest.TestCase):
    def test_source_mixing_and_blinding(self):
        n=720; k=np.arange(n); baseline=pd.Timestamp("2010-01-01"); rng=np.random.default_rng(71)
        def dates(days):return baseline+pd.to_timedelta(days,unit="D")
        # The same canonical first-year group contains cases with a raw GP prior, a composite-only
        # prior, an earlier ICD-10 version, baseline-day GP, and no other recorded history.
        kind=k%12
        day=np.select([kind<=2,kind<=5,kind<=10],[-800,900,120],default=-800)
        raw_icd=day.copy();raw_icd[kind==8]=-20
        gp=np.select([kind==6,kind==9],[-500,0],default=np.nan)
        p=pd.DataFrame(dict(eid=k.astype(str),computed_icd10=dates(raw_icd),computed_gp2=dates(gp),
                computed_srd=dates(np.where(kind==11,60,np.nan))))
        c=pd.DataFrame(dict(eid=k.astype(str),date_attend=baseline,birth_date="1950-01-01",
                            fod_icd10_cvd_cad=dates(day),fod_ref_cvd_cad=dates(np.where(kind==7,-500,day))))
        d=derive(c,p,"cvd_cad","2023-04-01")[0]
        d["group"]="f"+d.eid;d["baseline_lipid_i0"]=k%2
        x=rng.normal(size=(n,6)); x[:,0]+=(day<0)
        test=k>=480
        with TemporaryDirectory() as temp:
            tmp=Path(temp);out=analyze_icd_first(d,x,~test,test,tmp,30,2026)
            scores=pd.read_csv(tmp/"icd_first.tsv.gz",sep="\t")
            e=scores[scores.period=="first_year_ICD"]
            self.assertEqual(set(e.eid%12),{6,7,9,10})
            self.assertEqual(set(e.history),{"raw_GP_prior","REF_prior","day_zero_only","no_prior_record_found"})
            self.assertEqual(out["icd_score_status"].status.iloc[0],"estimated")
            overlap=out["icd_any1_overlap"]
            self.assertTrue(overlap.saved_ICD_first.eq("baseline_or_earlier").any())
            self.assertTrue(overlap.ICD_any_version_prior_or_day0.any())
            tables={"focus_icd_first_scores":out["icd_first_scores"].assign(trait="cvd_cad")}
            self.assertIn("Fig11_first_ICD.png",highlight_figures(tables,tmp))
            # Change secondary-source evidence in every split. Scores and frozen threshold must be identical.
            for col in ("src_gp_days","src_raw_gp_days","src_ref_days","src_srd_days"):d[col]=np.nan
            out2=analyze_icd_first(d,x,~test,test,tmp,30,2026)
            new=pd.read_csv(tmp/"icd_first.tsv.gz",sep="\t")
            np.testing.assert_allclose(scores.score,new.score)
            self.assertEqual(out["icd_score_status"].cut.iloc[0],out2["icd_score_status"].cut.iloc[0])
            self.assertTrue(out2["icd_first_tests"].query('estimator=="record_evidence_AUC"').status.eq("insufficient_test_groups").all())


if __name__=="__main__":unittest.main()
