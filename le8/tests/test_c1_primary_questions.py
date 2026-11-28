"""Source overlap, anchor blinding and prospective individual scoring."""
from pathlib import Path
import sys
import unittest
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"f"))
from dx_dates import derive
from c1_four_groups import first_year_groups, four_groups, GROUPS
from c1_personal_predictivity import personal_predictivity
from paper_focus import primary_figures


class PrimaryQuestions(unittest.TestCase):
    def test_first_icd_365_days_overlaps_and_assays(self):
        rng=np.random.default_rng(52)
        baseline=pd.Timestamp("2010-01-01")
        def dates(days):return baseline+pd.to_timedelta(days,unit="D")
        kind=np.tile(np.arange(15),48);i=np.arange(len(kind))
        day=np.select([kind==0,kind==1,kind==2,kind==3,kind==4,kind==5],
                      [-700,366,0,365,1,90],default=240)
        earlier=np.where(kind==10,-8,day)
        ref=np.where(np.isin(kind,[5,6]),-30,day)
        gp=np.select([np.isin(kind,[6,7]), kind==12],[-90,0],default=np.nan)
        srd=np.where(np.isin(kind,[6,8]),-150,np.nan)
        p=pd.DataFrame(dict(eid=i.astype(str),computed_icd10=dates(earlier),computed_gp2=dates(gp),
                            computed_srd=dates(srd),cnt_undated_srd=(kind==11).astype(int),
                            cnt_undated_icd10=(kind==9).astype(int),
                            opcs4=dates(np.where(kind==13,-85,np.nan))))
        c=pd.DataFrame(dict(eid=i.astype(str),date_attend=baseline,birth_date="1950-01-01",
                            fod_icd10_stroke=dates(day),fod_ref_stroke=dates(ref)))
        d,_,_=derive(c,p,"stroke","2025-01-01")
        d["group"]=(i//2).astype(str);test=i>=450
        eligible,flags=first_year_groups(d)
        self.assertEqual(set(kind[eligible]),{3,4,5,6,7,8,11,12,13,14})
        for name,k in zip(GROUPS,({3,4,14},{5,6},{6,7},{6,8})):
            self.assertEqual(set(kind[flags[name]]),k)
        self.assertEqual(flags.loc[kind==6].sum(axis=1).min(),3)
        x=rng.normal(size=(len(i),3));x[:,0]+=(day<0)
        x[kind==4,1]=np.nan
        with TemporaryDirectory() as temp:
            folder=Path(temp);out=four_groups(d,x,["a","b","c"],~test,test,folder,30,72)
            self.assertEqual(out["four_anchors"].status.iloc[0],"estimated")
            self.assertEqual(len(out["four_features"]),12)
            mask=flags[GROUPS[0]]
            b=out["four_features"].query('feature=="b" and group=="no_prior_any_source"').iloc[0]
            self.assertEqual(b.N,np.isfinite(x[mask,1]).sum())
            self.assertAlmostEqual(b["median"],np.nanmedian(x[mask,1]))
            self.assertTrue(out["four_features"].FDR.dropna().between(0,1).all())
            np.testing.assert_array_equal(out["four_groups"].N_all_assayed,flags.sum().values)
            public={f"focus_{k}":v.assign(trait="stroke",layer="prot") for k,v in out.items()}
            self.assertIn("Fig12_ICD_first_year_profiles.png",primary_figures(public,folder))
            scores=pd.read_csv(folder/"four_scores.tsv.gz",sep="\t").set_index("eid")
            # Removing other-source labels can change group membership but cannot change the ICD-only classifier.
            d["ref_days"]=d.icd_days;d["src_gp_days"]=np.nan;d["src_srd_days"]=np.nan
            out=four_groups(d,x,["a","b","c"],~test,test,folder,30,72)
            new=pd.read_csv(folder/"four_scores.tsv.gz",sep="\t").set_index("eid")
            np.testing.assert_array_equal(scores.similarity_score,new.similarity_score)

    def test_scores_are_blind_to_test_outcomes(self):
        rng=np.random.default_rng(411);n=1400;i=np.arange(n)
        x=rng.normal(size=(n,6));y=rng.random(n)<1/(1+np.exp(1-x[:,0]))
        test=i>=960;pre=(i%35==0);censor=(i%19==0)&~pre
        valid=i%23!=0;sex=i%2;geo=x[:,:2]/2
        d=pd.DataFrame(dict(eid=i.astype(str),group=(i//2).astype(str),
                role=np.where(test,"test","training"),
                icd_state=np.select([pre,censor,y],["prevalent","control","incident"],default="control"),
                icd_t2e=np.select([pre,censor,y],[np.nan,2.,4.],default=14.)))
        eligible=~pre;counts=valid.astype(int)*8
        with TemporaryDirectory() as temp:
            folder=Path(temp)
            def run():
                report=personal_predictivity(d,x,geo,sex,valid,counts,eligible,folder,1,673,30)
                return report,pd.read_csv(folder/"predictivity.tsv.gz",sep="\t")
            out,sc=run()
            self.assertEqual(sc.risk_10y.notna().sum(),eligible.sum())
            self.assertTrue(sc.risk_10y[pre].isna().all())
            self.assertTrue(sc.predictivity[censor].notna().all())
            self.assertTrue((sc.donor_families[valid&eligible]>0).any())
            self.assertTrue(sc.donor_families.dropna().between(0,32).all())
            self.assertTrue(sc.predictivity.dropna().between(0,1).all())
            np.testing.assert_allclose(sc.expected_Brier[~valid&eligible],sc.risk_only_Brier[~valid&eligible])
            np.testing.assert_allclose(1-sc.expected_Brier,sc.predictivity,equal_nan=True)
            self.assertEqual(out["predictivity_test"].status.tolist(),["estimated","estimated"])
            self.assertEqual(out["predictivity"].known_N.iloc[0],(test&eligible&~censor).sum())
            self.assertIn("Fig13_individual_predictivity.png",primary_figures({"focus_predictivity_test":
                out["predictivity_test"].assign(trait="stroke",layer="prot_met")},folder))
            change=test&eligible
            d.loc[change,"icd_state"]=np.where(y[change],"control","incident")
            d.loc[change,"icd_t2e"]=np.where(y[change],15.,3.)
            out2,sc2=run()
            pd.testing.assert_frame_equal(sc,sc2)
            self.assertNotEqual(out["predictivity"].observed_Brier.iloc[0],out2["predictivity"].observed_Brier.iloc[0])
            self.assertFalse(any("eid" in frame for frame in out.values()))
            d["icd_state"]="control";d["icd_t2e"]=15.
            sparse,sc=run()
            self.assertTrue(sc.predictivity.isna().all())
            self.assertTrue(sparse["predictivity_test"].status.eq("no_known_test_outcomes").all())


if __name__=="__main__":unittest.main()
