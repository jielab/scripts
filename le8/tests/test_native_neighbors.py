import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from dataclasses import replace
import numpy as np
from sklearn.decomposition import PCA
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"f"))
import abm
from c1_glmnet_neighbors import freeze,evaluate


class NativeIncrement(unittest.TestCase):
    def test_blinded_residuals_preserve_unsupported_and_pair_by_id(self):
        rng=np.random.default_rng(21); n=260
        x=rng.normal(size=(n,4)); x[-12:]+=50
        ids=np.arange(n).astype(str)
        cfg=replace(abm.Config(),k=10,min_twins=2,min_family_ess=1,dimensions=2)
        c=SimpleNamespace(parts={"tune_gate":np.arange(80,155),"tune_model":np.arange(45,80),"test":np.arange(155,n)},
             ids=ids,groups=ids,z=x,qc=np.ones(n,bool),cfg=cfg,w=np.ones(n),y=rng.integers(0,2,n),
             common_geometry=PCA(2).fit(x[:45]))
        p=.2+.6/(1+np.exp(-x[:,0]))
        arm,audit,_=freeze(c,p)
        c.y[155:]=1-c.y[155:]
        same=freeze(c,p)[0]; np.testing.assert_equal(same["probability"],arm["probability"])
        use=arm["supported"];self.assertTrue(use.any() and (~use).any())
        base=p[155:];np.testing.assert_equal(arm["probability"][~use],base[~use])
        qids=ids[155:]
        pred=dict(eid=qids,group=qids,probability=base,technical_qc=np.ones(len(qids),bool),
                  frozen_residual_arm=arm,metadata={"residual_bank":audit})
        event=np.arange(len(qids))%2; times=np.where(event,2.,12.)
        truth=dict(eid=qids[::-1],group=qids[::-1],time=times[::-1],event=event[::-1])
        km=abm.CensoringKM().fit(times,event,10.,.01)
        out=evaluate({"glmnet":pred},truth,km,22,2026)
        self.assertTrue(out[out.subset.eq("unsupported")].difference.eq(0).all())
        self.assertEqual(out[out.subset.eq("all_common")].test_known_N.iloc[0],len(qids))


if __name__=="__main__":unittest.main()
