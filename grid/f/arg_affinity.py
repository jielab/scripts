#!/usr/bin/env python3
"""Optional chromosome-averaged genealogical nearest-neighbour affinities."""
from __future__ import annotations
import argparse,gzip
from pathlib import Path
import numpy as np,pandas as pd,tskit
POPS=["AFR","EAS","EUR","SAS"]
def main():
  ap=argparse.ArgumentParser(); ap.add_argument("--trees",required=True); ap.add_argument("--sample-map",required=True); ap.add_argument("--anchors",required=True); ap.add_argument("--out",required=True); ap.add_argument("--batch",type=int,default=10000); a=ap.parse_args()
  ts=tskit.load(a.trees); m=pd.read_csv(a.sample_map,sep="\t",dtype={"eid":str}); anc=pd.read_csv(a.anchors,sep="\t")
  refs=[anc.loc[anc["pop"]==p,"sample_node"].astype(int).tolist() for p in POPS]; focal=m.sample_node.astype(int).to_numpy()
  Path(a.out).parent.mkdir(parents=True,exist_ok=True)
  with gzip.open(a.out,"wt") as o:
    o.write("eid\thap\tq_AFR\tq_EAS\tq_EUR\tq_SAS\n")
    for st in range(0,len(focal),a.batch):
      en=min(len(focal),st+a.batch); q=np.asarray(ts.genealogical_nearest_neighbours(focal[st:en],refs))
      q=np.nan_to_num(q,nan=0.0); den=q.sum(1); q=np.divide(q,den[:,None],out=np.full_like(q,.25),where=den[:,None]>0)
      for (_,r),v in zip(m.iloc[st:en].iterrows(),q): o.write(f"{r.eid}\t{r.hap}\t"+"\t".join(map(str,v))+"\n")
  print(f"haplotypes={len(focal)}")
if __name__=="__main__": main()
