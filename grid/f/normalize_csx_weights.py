#!/usr/bin/env python3
from __future__ import annotations
import argparse,re
from pathlib import Path
import pandas as pd

def read(path):
  # Standard PRS-CSx posterior files are whitespace-delimited, often without a header.
  with open(path,errors='replace') as h: first=h.readline().split()
  known={x.upper() for x in first}&{'CHR','SNP','BP','A1','A2','BETA'}
  if len(known)>=3: d=pd.read_csv(path,sep=r'\s+',engine='python')
  else:
    d=pd.read_csv(path,sep=r'\s+',engine='python',header=None)
    if d.shape[1]<6: raise SystemExit(f'Unrecognized PRS-CSx output {path}: {d.shape[1]} columns')
    d=d.iloc[:,:6]; d.columns=['CHR','SNP','BP','A1','A2','BETA']
  up={str(c).upper():c for c in d.columns}
  def c(*x): return next((up[y] for y in x if y in up),None)
  sc=c('SNP','RSID','ID'); ac=c('A1','EA','EFFECT_ALLELE'); bc=c('BETA','POSTERIOR_BETA','EFFECT')
  if not sc or not ac or not bc: raise SystemExit(f'Missing SNP/A1/BETA in {path}: {list(d.columns)}')
  out=pd.DataFrame({'SNP':d[sc].astype(str),'A1':d[ac].astype(str).str.upper(),'BETA':pd.to_numeric(d[bc],errors='coerce')})
  for new,aliases in [('CHR',('CHR','CHROM')),('BP',('BP','POS')),('A2',('A2','NEA'))]:
    z=c(*aliases); out[new]=d[z] if z else pd.NA
  return out.dropna(subset=['BETA'])[['SNP','A1','BETA','CHR','BP','A2']]
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--output',required=True); a=ap.parse_args(); d=read(a.input); d.to_csv(a.output,sep='\t',index=False); print(len(d))
if __name__=='__main__':main()
