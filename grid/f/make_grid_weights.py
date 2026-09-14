#!/usr/bin/env python3
"""Shrink ancestry-specific PRS-CSx effects toward a shared effect using the GRID prior."""
from __future__ import annotations
import argparse,math,re
from pathlib import Path
import numpy as np,pandas as pd
COMP=str.maketrans('ACGT','TGCA')
def comp(s):return str(s).translate(COMP)
def align(beta,a1,a2,r1,r2):
 a1=np.asarray(a1,str);a2=np.asarray(a2,str);r1=np.asarray(r1,str);r2=np.asarray(r2,str); b=np.asarray(beta,float)
 same=(a1==r1)&(a2==r2);swap=(a1==r2)&(a2==r1);cs=np.array([comp(x) for x in a1])==r1;cs&=np.array([comp(x) for x in a2])==r2;cw=np.array([comp(x) for x in a1])==r2;cw&=np.array([comp(x) for x in a2])==r1
 out=np.full(len(b),np.nan);out[same|cs]=b[same|cs];out[swap|cw]=-b[swap|cw];return out
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--weights-dir',required=True); ap.add_argument('--conservation',required=True); ap.add_argument('--manifest',required=True); ap.add_argument('--out-dir',required=True); ap.add_argument('--pops',default='AFR,EAS,EUR,SAS'); ap.add_argument('--chrs',required=True); a=ap.parse_args()
 pops=[x.upper() for x in re.split('[,; ]+',a.pops) if x];chrs=[x for x in re.split('[,; ]+',a.chrs) if x];out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
 man=pd.read_csv(a.manifest,sep='\t'); ns={str(r['pop']).upper():float(r['n_gwas']) for _,r in man.iterrows()}; nw={p:math.sqrt(max(ns.get(p,1),1)) for p in pops}
 con=pd.read_csv(a.conservation,sep='\t',compression='infer',dtype={'SNP':str,'chr':str})[['SNP','chr','conservation']].drop_duplicates(['SNP','chr'])
 audits=[]
 for c in chrs:
  ds={}
  for p in pops:
   f=Path(a.weights_dir)/f'{p}.chr{c}.tsv'
   if not f.exists():continue
   d=pd.read_csv(f,sep='\t',dtype={'SNP':str});d['A1']=d.A1.astype(str).str.upper();d['A2']=d.A2.astype(str).str.upper();d['BETA']=pd.to_numeric(d.BETA,errors='coerce');ds[p]=d.dropna(subset=['SNP','A1','BETA']).drop_duplicates('SNP')
  if len(ds)<2: raise SystemExit(f'Need >=2 PRS-CSx populations for chr{c}')
  priority=[p for p in ['EUR','AFR','EAS','SAS'] if p in ds]
  ref=pd.concat([ds[p][['SNP','A1','A2']].assign(prio=i) for i,p in enumerate(priority)],ignore_index=True).sort_values('prio').drop_duplicates('SNP').drop(columns='prio')
  z=ref.copy(); bet=[]
  for p in pops:
   if p not in ds: z[f'beta_{p}']=np.nan;continue
   x=z[['SNP','A1','A2']].merge(ds[p][['SNP','A1','A2','BETA']],on='SNP',how='left',suffixes=('_REF',''))
   z[f'beta_{p}']=align(x.BETA,x.A1,x.A2,x.A1_REF,x.A2_REF)
  B=np.column_stack([z[f'beta_{p}'].to_numpy(float) for p in pops]);W=np.array([nw[p] for p in pops],float)[None,:]*np.isfinite(B);shared=np.nansum(B*W,axis=1)/np.where(W.sum(1)>0,W.sum(1),np.nan)
  z['beta_shared']=shared;z['chr']=str(c);z=z.merge(con,on=['SNP','chr'],how='left');z['prior_missing']=z.conservation.isna();z['conservation']=pd.to_numeric(z.conservation,errors='coerce').fillna(0).clip(0,1)
  pd.DataFrame({'SNP':z.SNP,'A1':z.A1,'BETA':z.beta_shared}).dropna().to_csv(out/f'GRID_shared.chr{c}.tsv',sep='\t',index=False)
  for p in pops:
   bp=z[f'beta_{p}'].to_numpy(float); bp=np.where(np.isfinite(bp),bp,shared); bg=z.conservation.to_numpy()*shared+(1-z.conservation.to_numpy())*bp
   z[f'beta_GRID_{p}']=bg;pd.DataFrame({'SNP':z.SNP,'A1':z.A1,'BETA':bg}).dropna().to_csv(out/f'GRID_{p}.chr{c}.tsv',sep='\t',index=False)
  audits.append(z)
 pd.concat(audits,ignore_index=True).to_csv(out/'GRID_weight_audit.tsv.gz',sep='\t',index=False,compression='gzip',float_format='%.9g')
 print(f'chromosomes={len(chrs)} variants={sum(len(x) for x in audits)} shared_weight=sqrt_n conservation_missing_is_zero')
if __name__=='__main__':main()
