#!/usr/bin/env python3
"""Bridge permanent 1csx outputs to the GRID transport/weight input contract."""
import argparse,json
from pathlib import Path
import pandas as pd
from csx_config import find_gwas,POPS
from prepare_sumstats import preparation_signature
from pipeline_io import signature
from sumstats_cache import ensure_prepared

def main():
 p=argparse.ArgumentParser()
 for k in ('trait','gwas-dir','snpinfo','chrs','out'):p.add_argument('--'+k,required=True)
 a=p.parse_args();out=Path(a.out);chrs=list(map(int,a.chrs.split()))
 suffix='' if chrs==list(range(1,23)) else '.chr'+','.join(map(str,chrs))
 srcs=[find_gwas(a.gwas_dir,a.trait,pop) for pop in POPS]
 weights=[Path(str(src)[:-3]+suffix+'.csx.gz') for src in srcs]
 signatures=[];metas=[]
 for src,w,pop in zip(srcs,weights,POPS):
  sf=Path(str(w)+'.signature');mf=Path(str(w)+'.metadata.json')
  if not all(f.is_file() for f in (w,sf,mf)):raise ValueError(f'Run 1csx.sh --stage weights with the same chromosomes first: missing {w} or its sidecars')
  meta=json.loads(mf.read_text())
  if meta.get('preparation_signature')!=preparation_signature(src,a.snpinfo,a.trait,pop):raise ValueError(f'CSx source/preparation mismatch: {w}; rerun 1csx.sh')
  signatures.append(sf.read_text().strip());metas.append(meta)
 if len(set(signatures))!=1:raise ValueError('Population weights come from different joint CSx runs')
 files=srcs+weights+[Path(str(w)+ext) for w in weights for ext in ('.signature','.metadata.json')]+[Path(a.snpinfo),Path(__file__)]
 key=signature([a.trait,chrs],files);cache=out/'grid'/'inputs'/key;cache.mkdir(parents=True,exist_ok=True)
 dest=out/'sumstats'/'bychr';dest.mkdir(parents=True,exist_ok=True)
 wd=out/'csx'/'weights';wd.mkdir(parents=True,exist_ok=True);manifest=[]
 marker=out/'grid'/'inputs.signature'
 expected=[dest/f'{a.trait}.{pop}.chr{c}.tsv.gz' for pop in POPS for c in chrs]+[wd/f'{pop}.chr{c}.tsv' for pop in POPS for c in chrs]+[out/'csx'/'manifest.tsv']
 if marker.is_file() and marker.read_text().strip()==key and all(f.is_file() and f.stat().st_size for f in expected):
  print('SKIP GRID inputs: matching permanent CSx/GWAS inputs');return
 for pop,src,w,meta in zip(POPS,srcs,weights,metas):
  table=cache/f'{pop}.tsv.gz';info=cache/f'{pop}.json'
  ensure_prepared(src,a.snpinfo,a.trait,pop,output=table,metadata=info)
  wanted={c:[] for c in chrs}
  for chunk in pd.read_csv(table,sep='\t',chunksize=500000):
   for c in chrs:
    z=chunk.loc[chunk.CHR==c]
    if not z.empty:wanted[c].append(z)
  beta=pd.read_csv(w,sep='\t')
  for c in chrs:
   if not wanted[c] or not (beta.CHR==c).any():raise ValueError(f'No GWAS/CSx variants for {pop} chr{c}')
   pd.concat(wanted[c]).to_csv(dest/f'{a.trait}.{pop}.chr{c}.tsv.gz',sep='\t',index=False,compression='gzip')
   beta.loc[beta.CHR==c].to_csv(wd/f'{pop}.chr{c}.tsv',sep='\t',index=False)
  n=meta.get('n_gwas_used',meta['n_gwas_median'])
  if n is None or float(n)<=0:raise ValueError(f'Missing actual GWAS N for {pop}')
  manifest.append({'pop':pop,'n_gwas':n,'source':str(src),'weights':str(w),'joint_signature':signatures[0]})
 pd.DataFrame(manifest).to_csv(out/'csx'/'manifest.tsv',sep='\t',index=False)
 (out/'grid'/'inputs.signature').write_text(key+'\n')
 print(f'GRID inputs ready: {len(chrs)} chromosomes; joint CSx={signatures[0][:12]}',flush=True)

if __name__=='__main__':main()
