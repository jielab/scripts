#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np,pandas as pd
POPS=['AFR','EAS','EUR','SAS']
def read(path): return pd.read_csv(path,sep='\t',compression='infer',dtype={'eid':str})
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--scores',required=True);ap.add_argument('--ancestry',required=True);ap.add_argument('--prefix',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
 s=read(a.scores);anc=read(a.ancestry);idc=next((c for c in ['eid','IID','#IID','ID_2'] if c in anc),None)
 if idc is None: raise SystemExit('No ancestry ID')
 anc=anc.rename(columns={idc:'eid'})
 cols=[f'{a.prefix}_{p}' for p in POPS];missing=[c for c in cols if c not in s]
 if missing:raise SystemExit('Missing population scores '+','.join(missing))
 qcols=[f'posterior_{p}' for p in POPS]
 if all(c in anc for c in qcols):q=anc[['eid']+qcols].copy();q.columns=['eid']+[f'q_{p}' for p in POPS]
 else:
  ac=next((c for c in ['ancestry','genetic_ancestry','predicted_ancestry'] if c in anc),None)
  if ac is None:raise SystemExit('No ancestry/posterior columns')
  q=anc[['eid',ac]].copy();
  for p in POPS:q[f'q_{p}']=(q[ac].astype(str).str.upper()==p).astype(float)
  q=q.drop(columns=ac)
 z=s.merge(q,on='eid',how='left');Q=z[[f'q_{p}' for p in POPS]].apply(pd.to_numeric,errors='coerce').fillna(0).to_numpy();den=Q.sum(1);Q=np.divide(Q,den[:,None],out=np.full_like(Q,.25),where=den[:,None]>0);B=z[cols].apply(pd.to_numeric,errors='coerce').to_numpy()
 z[f'{a.prefix}_posterior']=np.nansum(B*Q,axis=1);ix=np.argmax(Q,axis=1);z[f'{a.prefix}_matched']=B[np.arange(len(B)),ix]
 out=z[['eid']+cols+[f'{a.prefix}_posterior',f'{a.prefix}_matched']];Path(a.out).parent.mkdir(parents=True,exist_ok=True);out.to_csv(a.out,sep='\t',index=False,compression='gzip');print(len(out))
if __name__=='__main__':main()
