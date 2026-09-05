#!/usr/bin/env python3
import argparse,gzip
from pathlib import Path
import pandas as pd
ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--out-dir',required=True); ap.add_argument('--prefix',required=True); ap.add_argument('--chrs',required=True); a=ap.parse_args()
chrs={str(x) for x in a.chrs.replace(',',' ').split()}; Path(a.out_dir).mkdir(parents=True,exist_ok=True); first={c:True for c in chrs}
for d in pd.read_csv(a.input,sep='\t',compression='infer',dtype={'CHR':str},chunksize=200000):
 d['CHR']=d.CHR.str.replace('chr','',case=False,regex=False)
 for c,x in d[d.CHR.isin(chrs)].groupby('CHR'):
  p=f'{a.out_dir}/{a.prefix}.chr{c}.tsv.gz'; x.to_csv(p,sep='\t',index=False,compression='gzip',mode='wt' if first[c] else 'at',header=first[c]); first[c]=False
for c in chrs:
 if first[c]: raise SystemExit(f'No variants for chromosome {c}')
