#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, re
from pathlib import Path
import pandas as pd


def first(cols, names):
    low={c.lower():c for c in cols}
    for n in names:
        if n in cols:return n
        if n.lower() in low:return low[n.lower()]
    return None


def norm_one(path:Path,node:str):
    try: d=pd.read_csv(path,sep='\t')
    except pd.errors.EmptyDataError: return pd.DataFrame()
    if d.empty:return d
    c=first(d.columns,['chrom','chromosome','chr']); s=first(d.columns,['start','start_bp','left','begin']); e=first(d.columns,['end','end_bp','right','stop'])
    p=first(d.columns,['mean_posterior','posterior','posterior_mean','prob','probability']); cm=first(d.columns,['length(cM)','length_cM','length_cm','genetic_length_cm'])
    if not all([c,s,e]): raise SystemExit(f'Unrecognized TRACE summary columns in {path}: {list(d.columns)}')
    out=pd.DataFrame({'tree_node_id':str(node),'chr':d[c].astype(str).str.replace('chr','',regex=False),'start':pd.to_numeric(d[s],errors='coerce'),'end':pd.to_numeric(d[e],errors='coerce')})
    out['length_bp']=out.end-out.start
    out['posterior']=pd.to_numeric(d[p],errors='coerce') if p else pd.NA
    out['length_cM']=pd.to_numeric(d[cm],errors='coerce') if cm else pd.NA
    return out[(out.start.notna())&(out.end>out.start)].copy()


def reduce_person(g):
    g=g.sort_values(['start','end']).copy(); rows=[]
    cs=ce=None; haps=set(); posts=[]
    for r in g.itertuples(index=False):
        if cs is None or r.start>ce:
            if cs is not None: rows.append((cs,ce,len(haps),max(posts) if posts else None))
            cs,ce=float(r.start),float(r.end); haps={str(r.haplotype)}; posts=[] if pd.isna(r.posterior) else [float(r.posterior)]
        else:
            ce=max(ce,float(r.end)); haps.add(str(r.haplotype));
            if not pd.isna(r.posterior): posts.append(float(r.posterior))
    if cs is not None: rows.append((cs,ce,len(haps),max(posts) if posts else None))
    return pd.DataFrame(rows,columns=['start','end','dosage','posterior'])


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True,type=Path); ap.add_argument('--sample-map',required=True,type=Path); ap.add_argument('--build',default='GRCh37'); args=ap.parse_args()
    m=pd.read_csv(args.sample_map,sep='\t',dtype={'tree_node_id':str}); pieces=[]
    for f in sorted((args.root/'summary').glob('hap*.summary.txt')):
        mt=re.match(r'hap(.+)\.summary\.txt$',f.name)
        if not mt: continue
        x=norm_one(f,mt.group(1))
        if not x.empty: pieces.append(x)
    hap=pd.concat(pieces,ignore_index=True) if pieces else pd.DataFrame(columns=['tree_node_id','chr','start','end','length_bp','posterior','length_cM'])
    hap=hap.merge(m[['tree_node_id','sample','haplotype']],on='tree_node_id',how='left')
    hap['method']='trace'; hap['source']='ghost_or_unknown_archaic'; hap['genome_build']=args.build
    hap=hap[['sample','haplotype','tree_node_id','method','source','chr','start','end','length_bp','posterior','length_cM','genome_build']]
    final=args.root/'final'; report=args.root/'report'; final.mkdir(parents=True,exist_ok=True); report.mkdir(parents=True,exist_ok=True)
    hap.to_csv(final/'trace_haplotype_segments.tsv.gz',sep='\t',index=False,compression='gzip')
    persons=[]
    if not hap.empty:
        for (sample,chrom),g in hap.groupby(['sample','chr'],dropna=False):
            z=reduce_person(g[['start','end','haplotype','posterior']])
            if z.empty: continue
            z.insert(0,'chr',chrom); z.insert(0,'sample',sample); persons.append(z)
    person=pd.concat(persons,ignore_index=True) if persons else pd.DataFrame(columns=['sample','chr','start','end','dosage','posterior'])
    if not person.empty: person['length_bp']=person.end-person.start
    else: person['length_bp']=pd.Series(dtype=float)
    person['method']='trace'; person['source']='ghost_or_unknown_archaic'; person['genome_build']=args.build
    person.to_csv(final/'trace_person_segments.tsv.gz',sep='\t',index=False,compression='gzip')
    if hap.empty:
        burden=pd.DataFrame(columns=['sample','n_haplotype_segments','total_haplotype_bp','mean_posterior'])
    else:
        burden=hap.groupby('sample').agg(n_haplotype_segments=('start','size'),total_haplotype_bp=('length_bp','sum'),mean_posterior=('posterior','mean')).reset_index()
    burden.to_csv(report/'trace_sample_burden.tsv',sep='\t',index=False)
    print(f'TRACE combine: hap_segments={len(hap)} person_segments={len(person)} samples={burden.shape[0]}')
if __name__=='__main__': main()
