#!/usr/bin/env python3
from __future__ import annotations
import argparse,gzip,re
from pathlib import Path
import pandas as pd

def read_score(path):
 d=pd.read_csv(path,sep=r'\s+',engine='python',dtype=str); idc=next((c for c in ['IID','#IID','ID_2','eid'] if c in d.columns),None)
 if idc is None: raise SystemExit(f'No IID in {path}')
 cols=[c for c in d.columns if c.endswith('_SUM') and c not in {'NAMED_ALLELE_DOSAGE_SUM'}]
 if not cols: cols=[c for c in d.columns if c.startswith('SCORE') and c not in {'SCORE_AVG'}]
 if not cols: raise SystemExit(f'No score sum in {path}: {list(d.columns)}')
 return pd.DataFrame({'eid':d[idc].astype(str),'score':pd.to_numeric(d[cols[-1]],errors='coerce').fillna(0)})
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--inputs',nargs='+',required=True); ap.add_argument('--name',required=True); ap.add_argument('--output',required=True); a=ap.parse_args()
 z=None
 for p in a.inputs:
  d=read_score(p); z=d if z is None else z.merge(d,on='eid',how='outer',suffixes=('','_x')).fillna(0).assign(score=lambda x:x['score']+x['score_x']).drop(columns='score_x')
 z=z.rename(columns={'score':a.name}); z.to_csv(a.output,sep='\t',index=False,compression='gzip'); print(len(z))
if __name__=='__main__':main()
