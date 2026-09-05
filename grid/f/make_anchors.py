#!/usr/bin/env python3
from __future__ import annotations
import argparse
import pandas as pd

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--panel",required=True); ap.add_argument("--sample-map",required=True); ap.add_argument("--out",required=True); ap.add_argument("--per-pop",type=int,default=1000); ap.add_argument("--prob-min",type=float,default=.999); a=ap.parse_args()
    p=pd.read_csv(a.panel,sep="\t",dtype={"eid":str}); m=pd.read_csv(a.sample_map,sep="\t",dtype={"eid":str})
    p["pop"]=p["pop"].astype(str).str.upper(); p["prob"]=pd.to_numeric(p.get("prob"),errors="coerce")
    keep=[]
    for pop in ["AFR","EAS","EUR","SAS"]:
      x=p[(p["pop"]==pop)&((p["prob"]>=a.prob_min)|p["prob"].isna())].sort_values(["prob","source_order"],ascending=[False,True]).head(a.per_pop)
      if len(x)<10: raise SystemExit(f"Only {len(x)} anchors for {pop}")
      y=m.merge(x[["eid","prob"]],on="eid",how="inner"); y["pop"]=pop; keep.append(y)
    z=pd.concat(keep,ignore_index=True); z[["eid","hap","sample_node","pop","prob"]].to_csv(a.out,sep="\t",index=False)
    print(z.groupby("pop").size().to_dict())
if __name__=="__main__": main()
