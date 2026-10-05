"""Generate labelled SYNTHETIC C1 smoke inputs. Never reads UKB or trains a model."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n',type=int,default=6000)
    ap.add_argument('--features',type=int,default=24)
    ap.add_argument('--seed',type=int,default=2026)
    a=ap.parse_args()
    if a.n < 1000 or a.n%2 or a.features < 8:
        ap.error('Use an even N >= 1000 and >= 8 features')
    a.out.mkdir(parents=True,exist_ok=True)
    names=['phe.csv','prot.csv','groups.csv','outer.csv','SIMULATION_NOTICE.json']
    if any((a.out/n).exists() for n in names):
        raise FileExistsError('Refusing to overwrite prior inputs; choose a new directory')
    rng=np.random.default_rng(a.seed)
    ids=[f'p{i:07d}' for i in range(a.n)]
    group=np.array([f'fam{i//2:06d}' for i in range(a.n)])
    age=rng.uniform(40,70,a.n); sex=rng.integers(0,2,a.n)
    x=rng.normal(size=(a.n,a.features))
    rare=rng.random(a.n)<.003
    x[:,5]=20.*rare
    x[:,6]=rng.integers(0,2,a.n);x[:,7]=rng.integers(0,2,a.n)
    family_effect=np.repeat(rng.normal(0,.15,a.n//2),2)
    # Enough events for role/interface smoke; not a realistic CAD prevalence model.
    rate=.10*np.exp(.03*(age-55)+.25*sex+.35*x[:,0]-.25*x[:,1]+.2*(x[:,6]==x[:,7])+family_effect+.5*rare)
    failure=rng.exponential(1/rate);censor=rng.uniform(6,9,a.n)
    event=failure<=censor
    baseline=pd.Timestamp('2010-01-01')
    end=baseline+pd.to_timedelta(np.rint(censor*365.25).astype(int),unit='D')
    # Only convert observed diagnoses: rare enormous hypothetical times are not material.
    diagnosis=np.full(a.n,'',dtype=object)
    days=np.maximum(1,np.rint(failure[event]*365.25).astype(int))
    diagnosis[event]=(baseline+pd.to_timedelta(days,unit='D')).strftime('%Y-%m-%d')
    phe=pd.DataFrame(dict(eid=ids,age=age,sex=sex,date_attend='2010-01-01',
        date_lost=end.strftime('%Y-%m-%d'),date_death='',fod_icd10_cvd_cad=diagnosis))
    omics=pd.DataFrame(x,columns=[f'F{i:03d}' for i in range(a.features)])
    omics.insert(0,'eid',ids)
    # Fold definition depends only on family IDs and seed, not failure outcomes.
    levels=sorted(set(group),key=lambda v:hashlib.blake2b(f'{a.seed}:{v}'.encode(),digest_size=8).digest())
    testing=set(levels[:round(.2*len(levels))])
    roles=np.array(['test' if v in testing else 'training' for v in group])
    phe.to_csv(a.out/'phe.csv',index=False)
    omics.to_csv(a.out/'prot.csv',index=False)
    pd.DataFrame(dict(eid=ids,group=group)).to_csv(a.out/'groups.csv',index=False)
    pd.DataFrame(dict(eid=ids,role=roles)).to_csv(a.out/'outer.csv',index=False)
    notice=dict(synthetic=True,clinical_use=False,N=a.n,features=a.features,seed=a.seed,
        events_observed=int(event.sum()),events_by_5years=int(((failure<=5)&event).sum()),
        test_N=int((roles=='test').sum()),rare_count=int(rare.sum()),
        scope='Input/interface smoke only; no UKB data and no model performance claim')
    (a.out/'SIMULATION_NOTICE.json').write_text(json.dumps(notice,indent=2)+'\n')
    print(json.dumps(notice,indent=2))

if __name__=='__main__':main()
