#!/usr/bin/env python3
"""Extract mutation age and population-level local genealogy from a tskit ARG."""
from __future__ import annotations
import argparse,gzip,math
from pathlib import Path
import numpy as np, pandas as pd, tskit
POPS=["AFR","EAS","EUR","SAS"]
PAIRS=[("AFR","EAS"),("AFR","EUR"),("AFR","SAS"),("EAS","EUR"),("EAS","SAS"),("EUR","SAS")]

def op(path): return gzip.open(path,"rt") if str(path).endswith(".gz") else open(path)
def read_haps(path):
    d={}
    with op(path) as h:
      for line in h:
        if not line.strip(): continue
        z=line.split()
        if len(z)<6: continue
        if z[2].replace('.','',1).isdigit(): chrom,sid,bp,a0,a1=z[0],z[1],int(float(z[2])),z[3],z[4]
        elif z[3].replace('.','',1).isdigit(): chrom,sid,bp,a0,a1=z[0],z[2],int(float(z[3])),z[4],z[5]
        else: continue
        d.setdefault(bp,(sid,a0,a1))
    return d

def stat_arrays(ts,sets,windows):
    pairs=[(POPS.index(a),POPS.index(b)) for a,b in PAIRS]
    try: div=ts.divergence(sample_sets=sets,indexes=pairs,windows=windows,mode="branch",span_normalise=True)
    except TypeError: div=ts.divergence(sets,indexes=pairs,windows=windows,mode="branch")
    try: within=ts.diversity(sample_sets=sets,windows=windows,mode="branch",span_normalise=True)
    except TypeError: within=ts.diversity(sets,windows=windows,mode="branch")
    div=np.asarray(div); within=np.asarray(within)
    if div.ndim==1: div=div[:,None]
    if within.ndim==1: within=within[:,None]
    return div,within

def root_by_window(ts,windows):
    n=len(windows)-1; acc=np.zeros(n); span=np.zeros(n)
    for tree in ts.trees():
      l,r=tree.interval; roots=list(tree.roots); rt=max((ts.node(x).time for x in roots),default=np.nan)
      if not np.isfinite(rt): continue
      a=max(0,int(np.searchsorted(windows,l,side="right")-1)); b=min(n-1,int(np.searchsorted(windows,r,side="left")))
      for w in range(a,b+1):
        ov=max(0,min(r,windows[w+1])-max(l,windows[w]))
        if ov: acc[w]+=ov*rt; span[w]+=ov
    return np.divide(acc,span,out=np.full(n,np.nan),where=span>0)

def mutation_info(ts,site):
    tree=ts.at(site.position); ages=[]; lengths=[]
    for m in site.mutations:
      node=int(m.node); nt=ts.node(node).time; par=tree.parent(node); pt=ts.node(par).time if par!=tskit.NULL else np.nan
      mt=getattr(m,"time",np.nan)
      if mt is not None and np.isfinite(mt) and mt!=tskit.UNKNOWN_TIME: age=float(mt)
      elif np.isfinite(pt): age=float((nt+pt)/2)
      else: age=float(nt)
      ages.append(age); lengths.append(float(pt-nt) if np.isfinite(pt) else np.nan)
    return (max(ages) if ages else np.nan, max([x for x in lengths if np.isfinite(x)],default=np.nan), len(site.mutations))

def entropy(v):
    x=np.asarray(v,float); x=x[x>0]
    if len(x)<=1:return 0.0
    p=x/x.sum(); return float(-(p*np.log(p)).sum()/math.log(len(POPS)))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--trees",required=True); ap.add_argument("--haps",required=True); ap.add_argument("--anchors",required=True); ap.add_argument("--out",required=True); ap.add_argument("--chr",required=True); ap.add_argument("--window-bp",type=int,default=1000000); a=ap.parse_args()
    ts=tskit.load(a.trees); hv=read_haps(a.haps); anc=pd.read_csv(a.anchors,sep="\t")
    sets=[anc.loc[anc["pop"]==p,"sample_node"].astype(int).tolist() for p in POPS]
    if any(len(x)<2 for x in sets): raise SystemExit("Each population needs >=2 anchor haplotypes")
    L=float(ts.sequence_length)
    windows=np.arange(0,L,a.window_bp,dtype=float)
    if len(windows)==0 or windows[0] != 0: windows=np.insert(windows,0,0.0)
    if windows[-1] != L: windows=np.append(windows,L)
    offset=int(ts.metadata.get("offset",0)) if isinstance(ts.metadata,dict) else 0
    div,within=stat_arrays(ts,sets,windows); root=root_by_window(ts,windows)
    all_nodes=[x for s in sets for x in s]; slices=[]; k=0
    for s in sets: slices.append(slice(k,k+len(s))); k+=len(s)
    sample_nodes=np.asarray(ts.samples(),dtype=int); node_to_index={n:i for i,n in enumerate(sample_nodes)}
    anchor_sample_idx=np.array([node_to_index[n] for n in all_nodes],dtype=int)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    cols=["chr","bp","SNP","allele0","allele1","arg_age_gen","mutation_branch_gen","n_mutations","root_time_gen","carrier_frequency","lineage_breadth","lineage_entropy"]
    cols += [f"af_{p}" for p in POPS]+[f"div_{x}_{y}" for x,y in PAIRS]+[f"within_{p}" for p in POPS]
    n=0; missing_id=0
    with gzip.open(a.out,"wt") as o:
      o.write("\t".join(cols)+"\n")
      for var in ts.variants(samples=sample_nodes,isolated_as_missing=False):
        site=var.site; bp=int(round(site.position))+offset; info=hv.get(bp)
        if info is None:
          # ARG conversions can use 0-based site coordinates while HAPS is 1-based.
          info=hv.get(bp+1) or hv.get(bp-1)
        if info is None: sid,a0,a1=f"{a.chr}:{bp}","NA","NA"; missing_id+=1
        else: sid,a0,a1=info
        wi=min(len(windows)-2,max(0,int(np.searchsorted(windows,site.position,side="right")-1)))
        age,bl,nmut=mutation_info(ts,site)
        g=np.asarray(var.genotypes)[anchor_sample_idx]
        counts=[]; af=[]
        for sl in slices:
          z=g[sl]; z=z[z>=0]; c=float(np.sum(z>0)); counts.append(c); af.append(c/max(1,len(z)))
        breadth=sum(x>0 for x in counts); freq=float(np.mean(g>0))
        row=[a.chr,bp,sid,a0,a1,age,bl,nmut,root[wi],freq,breadth,entropy(counts)]
        row += af+list(np.asarray(div[wi]).ravel())+list(np.asarray(within[wi]).ravel())
        o.write("\t".join("NA" if (isinstance(x,float) and not np.isfinite(x)) else str(x) for x in row)+"\n"); n+=1
    print(f"sites={n} missing_haps_id={missing_id} windows={len(windows)-1}")
if __name__=="__main__": main()
