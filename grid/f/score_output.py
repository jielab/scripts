#!/usr/bin/env python3
"""Publish compact score tables, excluding withdrawn individuals."""
import argparse
from pathlib import Path
import pandas as pd


def excluded_ids(path):
    if not path or not Path(path).is_file():
        return set()
    return {parts[1] if len(parts)>1 else parts[0]
            for line in Path(path).read_text().splitlines()
            if (parts := line.split()) and not parts[0].startswith('#')}


def filter_samples(d, idcol, remove):
    ids = d[idcol].astype(str)
    return d.loc[~ids.str.startswith('-') & ~ids.isin(excluded_ids(remove))].copy()


def publish_table(src, dest, method, remove):
    d = pd.read_csv(src, sep='\t', dtype={'eid':str, 'IID':str})
    idcol = 'eid' if 'eid' in d else 'IID'
    before = len(d)
    d = filter_samples(d, idcol, remove)
    mapping = {f'CSX_{p}':f'csx.{p}' for p in ['AFR','EAS','EUR','SAS']}
    if method == 'disco': mapping = {'PRS':'disco'}
    d = d.rename(columns=mapping)
    if d.empty or d[idcol].duplicated().any(): raise ValueError('Empty/duplicate score IDs')
    dest=Path(dest); dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_name(dest.name+'.tmp')
    d.to_csv(tmp,sep='\t',index=False,compression='gzip');tmp.replace(dest)
    print(f'{method}: retained={len(d)}; excluded={before-len(d)}; output={dest}')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('method',choices=['csx','disco'])
    p.add_argument('input');p.add_argument('output');p.add_argument('--remove',default='/mnt/d/files/ukb.exclude.id')
    a=p.parse_args();publish_table(a.input,a.output,a.method,a.remove)
