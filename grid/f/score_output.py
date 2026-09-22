#!/usr/bin/env python3
"""Publish compact score tables, excluding withdrawn individuals."""
import argparse
import fcntl,os
from pathlib import Path
import pandas as pd
import numpy as np


def excluded_ids(path):
    if not path:
        return set()
    if not Path(path).is_file():raise FileNotFoundError(f'Missing withdrawal list: {path}')
    return {parts[1] if len(parts)>1 else parts[0]
            for line in Path(path).read_text().splitlines()
            if (parts := line.split()) and not parts[0].startswith('#')}


def filter_samples(d, idcol, remove):
    ids = d[idcol].astype(str)
    return d.loc[~ids.str.startswith('-') & ~ids.isin(excluded_ids(remove))].copy()

def update_csx_table(d,dest,remove=''):
    """Atomically replace only supplied score columns, preserving other models."""
    d=d.rename(columns={'IID':'eid'}).copy();dest=Path(dest)
    dest.parent.mkdir(parents=True,exist_ok=True)
    with Path(str(dest)+'.lock').open('a') as lock:
      fcntl.flock(lock,fcntl.LOCK_EX)
      if dest.exists():
        old=pd.read_csv(dest,sep='\t',dtype={'eid':str,'IID':str}).rename(columns={'IID':'eid'})
        old=filter_samples(old,'eid',remove)
        if old.eid.isna().any() or old.eid.duplicated().any():raise ValueError('Invalid existing CSx IDs')
        if set(old.eid)!=set(d.eid):raise ValueError('New and existing CSx sample sets differ; use a separate --score-dir for a different cohort')
        extra=[c for c in old if c!='eid' and c not in d]
        d=d.merge(old[['eid']+extra],on='eid',validate='one_to_one')
      if d.empty or d.eid.isna().any() or d.eid.duplicated().any():raise ValueError('Empty/missing/duplicate CSx IDs')
      if not np.isfinite(d.drop(columns='eid').to_numpy(float)).all():raise ValueError('Nonfinite CSx scores')
      tmp=dest.with_name(dest.name+f'.tmp.{os.getpid()}')
      d.to_csv(tmp,sep='\t',index=False,compression='gzip');tmp.replace(dest)
    return d


def publish_table(src, dest, method, remove):
    d = pd.read_csv(src, sep='\t', dtype={'eid':str, 'IID':str})
    idcol = 'eid' if 'eid' in d else 'IID'
    before = len(d)
    d = filter_samples(d, idcol, remove)
    mapping = {f'CSX_{p}':f'csx.{p}' for p in ['AFR','EAS','EUR','SAS']}
    if method == 'disco': mapping = {'PRS':'disco'}
    d = d.rename(columns=mapping)
    if d.empty or d[idcol].isna().any() or d[idcol].duplicated().any(): raise ValueError('Empty/missing/duplicate score IDs')
    if method=='csx':
        update_csx_table(d,dest,remove)
        print(f'{method}: retained={len(d)}; excluded={before-len(d)}; output={dest}')
        return
    if not np.isfinite(d.drop(columns=idcol).to_numpy(float)).all():raise ValueError('Nonfinite scores')
    dest=Path(dest); dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_name(dest.name+'.tmp')
    d.to_csv(tmp,sep='\t',index=False,compression='gzip');tmp.replace(dest)
    print(f'{method}: retained={len(d)}; excluded={before-len(d)}; output={dest}')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('method',choices=['csx','disco'])
    p.add_argument('input');p.add_argument('output');p.add_argument('--remove',default='/mnt/d/files/ukb.exclude.id')
    a=p.parse_args();publish_table(a.input,a.output,a.method,a.remove)
