#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd, tskit

def sample_ids(path):
    lines=Path(path).read_text().splitlines(); h=lines[0].split(); rows=[x.split() for x in lines[2:] if x.strip()]
    ix={x.lower():i for i,x in enumerate(h)}; j=ix.get("id_2",ix.get("iid",ix.get("id",ix.get("id_1",min(1,len(h)-1)))))
    return [r[j] for r in rows]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--trees",required=True); ap.add_argument("--sample",required=True); ap.add_argument("--out",required=True); a=ap.parse_args()
    ids=sample_ids(a.sample); ts=tskit.load(a.trees); nodes=list(map(int,ts.samples()))
    if len(nodes)!=2*len(ids): raise SystemExit(f"tree samples={len(nodes)} != 2*individuals={2*len(ids)}")
    out=[]
    for i,eid in enumerate(ids):
        out.append((eid,0,nodes[2*i])); out.append((eid,1,nodes[2*i+1]))
    pd.DataFrame(out,columns=["eid","hap","sample_node"]).to_csv(a.out,sep="\t",index=False)
    print(f"individuals={len(ids)} sample_nodes={len(nodes)}")
if __name__=="__main__": main()
