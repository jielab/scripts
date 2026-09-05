#!/usr/bin/env python3
"""Interpolate a GRCh37 cumulative genetic map onto every Oxford HAPS variant."""
from __future__ import annotations
import argparse,gzip,re
from pathlib import Path
import numpy as np

def op(path): return gzip.open(path,"rt") if str(path).endswith(".gz") else open(path)
def numeric(x):
    try: return float(x)
    except: return None

def read_source(path):
    rows=[]; header=None
    with op(path) as h:
      for line in h:
        if not line.strip() or line.lstrip().startswith("#"): continue
        z=re.split(r"\s+",line.strip())
        if header is None and any(numeric(x) is None for x in z): header=[x.lower() for x in z]; continue
        rows.append(z)
    if not rows: raise SystemExit(f"No map rows in {path}")
    if header:
        def idx(pred): return next((i for i,x in enumerate(header) if pred(x)),None)
        ip=idx(lambda x:x in {"position","pos","bp","base_pair_location"} or "position" in x)
        ic=idx(lambda x:("cm" in x and "rate" not in x) or x in {"map","genetic_map"} or ("map" in x and "rate" not in x))
        ir=idx(lambda x:"rate" in x)
    else: ip=ic=ir=None
    vals=[]
    for z in rows:
        nums=[numeric(x) for x in z]
        if ip is not None and ip<len(nums) and nums[ip] is not None: bp=nums[ip]
        elif len(nums)>=4 and nums[-1] is not None: bp=nums[-1]
        else: bp=next((x for x in nums if x is not None and x>=1),None)
        if ic is not None and ic<len(nums): cm=nums[ic]
        elif len(nums)>=4: cm=nums[-2]
        elif len(nums)>=3: cm=nums[-1]
        else: cm=None
        rate=nums[ir] if ir is not None and ir<len(nums) else None
        if bp is not None: vals.append((float(bp),None if cm is None else float(cm),rate))
    vals=sorted({int(bp):(bp,cm,rate) for bp,cm,rate in vals}.values())
    bp=np.array([x[0] for x in vals],float); cm=np.array([np.nan if x[1] is None else x[1] for x in vals],float)
    if np.isfinite(cm).sum()<2:
        rate=np.array([np.nan if x[2] is None else x[2] for x in vals],float)
        if np.isfinite(rate).sum()<2: raise SystemExit("Map needs cumulative cM or recombination rate")
        rate=np.interp(bp,bp[np.isfinite(rate)],rate[np.isfinite(rate)])
        cm=np.zeros(len(bp)); cm[1:]=np.cumsum((bp[1:]-bp[:-1])*(rate[1:]+rate[:-1])/2/1e6)
    else: cm=np.interp(bp,bp[np.isfinite(cm)],cm[np.isfinite(cm)])
    if np.any(np.diff(bp)<0) or np.any(np.diff(cm)<-1e-6): raise SystemExit("Source map is not monotonic")
    return bp,cm

def haps_variants(path):
    with op(path) as h:
      for line in h:
        if not line.strip(): continue
        z=line.split()
        # Oxford/PLINK HAPS: chr ID BP A0 A1 ...; accept chr ID rsID BP A0 A1 ...
        if len(z)<6: raise SystemExit("Malformed HAPS")
        if numeric(z[2]) is not None: chrom,sid,pos=z[0],z[1],int(float(z[2]))
        elif numeric(z[3]) is not None: chrom,sid,pos=z[0],z[2],int(float(z[3]))
        else: raise SystemExit(f"Cannot identify HAPS position: {' '.join(z[:6])}")
        yield chrom,sid,pos

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--source",required=True); ap.add_argument("--haps",required=True); ap.add_argument("--out",required=True); ap.add_argument("--chr",required=True); a=ap.parse_args()
    bp,cm=read_source(a.source); Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    n=0; prev=-1
    with open(a.out,"w") as o:
      for chrom,sid,pos in haps_variants(a.haps):
        if pos<prev: raise SystemExit("HAPS positions are not sorted")
        prev=pos; g=float(np.interp(pos,bp,cm,left=cm[0],right=cm[-1])); o.write(f"{a.chr}\t{sid}\t{g:.9f}\t{pos}\n"); n+=1
    if n<2: raise SystemExit("Too few mapped variants")
    print(f"variants={n} source_points={len(bp)} map={a.out}")
if __name__=="__main__": main()
