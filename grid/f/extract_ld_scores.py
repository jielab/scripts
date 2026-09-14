#!/usr/bin/env python3
"""Extract per-SNP LD scores (sum r^2, diagonal included) from PRS-CSx HDF5 blocks."""
from __future__ import annotations
import argparse,gzip
from pathlib import Path
import h5py,numpy as np,pandas as pd

def dec(x):
 if isinstance(x,(bytes,np.bytes_)): return x.decode()
 if isinstance(x,np.void) and x.dtype.names:
  for n in x.dtype.names:
   if 'snp' in n.lower() or 'rs' in n.lower() or 'id'==n.lower(): return dec(x[n])
 return str(x)
def groups(h):
 out=[]
 def visit(name,obj):
  if isinstance(obj,h5py.Group):
   ds={k:v for k,v in obj.items() if isinstance(v,h5py.Dataset)}
   mats=[(k,v) for k,v in ds.items() if v.ndim==2 and v.shape[0]==v.shape[1]]
   vec=[(k,v) for k,v in ds.items() if v.ndim==1]
   for mk,m in mats:
    sv=next(((k,v) for k,v in vec if len(v)==m.shape[0] and any(z in k.lower() for z in ['snp','rs','id'])),None)
    if sv is None: sv=next(((k,v) for k,v in vec if len(v)==m.shape[0]),None)
    if sv: out.append((name,mk,m,sv[0],sv[1])); break
 h.visititems(visit); return out
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--hdf5',required=True); ap.add_argument('--out',required=True); ap.add_argument('--pop',required=True); ap.add_argument('--chr',required=True); a=ap.parse_args()
 rows=[]
 with h5py.File(a.hdf5,'r') as h:
  gs=groups(h)
  if not gs: raise SystemExit(f'No square LD matrix + SNP list found in {a.hdf5}')
  for name,mk,m,sk,s in gs:
   x=np.asarray(m[...],dtype=np.float64); ld=np.einsum('ij,ij->i',x,x); ids=[dec(z) for z in s[...]]
   rows.extend(zip(ids,ld))
 d=pd.DataFrame(rows,columns=['SNP','ldscore']).drop_duplicates('SNP'); d['pop']=a.pop.upper(); d['chr']=str(a.chr)
 Path(a.out).parent.mkdir(parents=True,exist_ok=True); d.to_csv(a.out,sep='\t',index=False,compression='gzip'); print(f'blocks={len(gs)} snps={len(d)}')
if __name__=='__main__':main()
