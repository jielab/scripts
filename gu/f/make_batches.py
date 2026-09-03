#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib
from pathlib import Path
import pandas as pd

def stable(xs): return sorted(map(str,xs), key=lambda x: hashlib.sha256(x.encode()).hexdigest())

def main():
    ap=argparse.ArgumentParser(description='Deterministic ancestry-stratified target batches with fixed anchors')
    ap.add_argument('--panel',required=True,type=Path); ap.add_argument('--outdir',required=True,type=Path)
    ap.add_argument('--batch-size',type=int,default=1000); ap.add_argument('--sample-col',default='sample'); ap.add_argument('--group-col',default='super_pop')
    ap.add_argument('--anchor-list',type=Path); ap.add_argument('--anchors-per-group',type=int,default=0)
    a=ap.parse_args(); d=pd.read_csv(a.panel,sep=None,engine='python',dtype=str)
    if a.sample_col not in d: raise SystemExit(f'missing sample column {a.sample_col}; columns={list(d)}')
    if a.group_col not in d: d[a.group_col]='ALL'
    d=d[[a.sample_col,a.group_col]].dropna(subset=[a.sample_col]).drop_duplicates(a.sample_col)
    a.outdir.mkdir(parents=True,exist_ok=True)
    anchors=[]
    if a.anchor_list:
        anchors=[x.split()[0] for x in a.anchor_list.read_text().splitlines() if x.strip() and not x.startswith('#')]
    elif a.anchors_per_group>0:
        for _,g in d.groupby(a.group_col,dropna=False): anchors += stable(g[a.sample_col])[:a.anchors_per_group]
        anchors=sorted(set(anchors)); (a.outdir/'anchors.samples.txt').write_text('\n'.join(anchors)+'\n')
    aset=set(anchors); rows=[]
    for group,g in d.groupby(a.group_col,dropna=False):
        group='NA' if pd.isna(group) else str(group); ids=[x for x in stable(g[a.sample_col]) if x not in aset]
        for k in range(0,len(ids),a.batch_size):
            target=ids[k:k+a.batch_size]; bid=f'{group}.b{k//a.batch_size+1:04d}'
            tf=a.outdir/f'{bid}.targets.txt'; tf.write_text('\n'.join(target)+'\n')
            jf=a.outdir/f'{bid}.joint.txt'; jf.write_text('\n'.join(sorted(set(anchors+target)))+'\n')
            rows.append([bid,group,len(target),len(anchors),str(tf),str(jf)])
    pd.DataFrame(rows,columns=['batch_id','group','n_target','n_anchor','target_list','joint_list']).to_csv(a.outdir/'batch_manifest.tsv',sep='\t',index=False)
    print(f'batches={len(rows)} target_samples={sum(r[2] for r in rows)} anchors={len(anchors)} manifest={a.outdir/"batch_manifest.tsv"}')
if __name__=='__main__': main()
