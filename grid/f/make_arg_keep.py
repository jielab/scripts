#!/usr/bin/env python3
"""Create a deterministic, ancestry-balanced UKB ARG sample order."""
from __future__ import annotations
import argparse, gzip
from pathlib import Path
import numpy as np
import pandas as pd

def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", compression="infer", dtype=str)

def read_sample(path: str) -> pd.DataFrame:
    lines=[x for x in Path(path).read_text().splitlines() if x.strip()]
    if len(lines)<2: raise SystemExit(f"Bad Oxford sample/PSAM file: {path}")
    h=lines[0].lstrip("#").split()
    # Oxford SAMPLE has a type row after the header; PLINK PSAM does not.
    start=2 if lines[1].split() and all(x in {"0","D","B","C","P"} for x in lines[1].split()) else 1
    rows=[x.split() for x in lines[start:] if x.strip()]
    if not rows: raise SystemExit(f"No samples in {path}")
    d=pd.DataFrame(rows, columns=h[:len(rows[0])])
    lower={x.lower():x for x in d.columns}
    col=lower.get("id_2") or lower.get("iid") or lower.get("id") or lower.get("id_1") or d.columns[min(1,len(d.columns)-1)]
    out=pd.DataFrame({"eid":d[col].astype(str)})
    out=out[~out.eid.isin(["0","NA","nan",""])].drop_duplicates("eid")
    out["source_order"]=np.arange(len(out),dtype=np.int64)
    return out

def norm_pop(x: pd.Series) -> pd.Series:
    z=x.fillna("").astype(str).str.upper().str.strip().str.replace("_"," ",regex=False).str.replace("-"," ",regex=False)
    out=pd.Series("OTH",index=x.index,dtype="object")
    rules={
      "AFR":r"^(?:AFR|AFRICAN|BLACK)", "EAS":r"^(?:EAS|EAST ASIAN|CHINESE|JAPANESE|KOREAN)",
      "EUR":r"^(?:EUR|EUROPEAN|WHITE)", "SAS":r"^(?:SAS|SOUTH ASIAN|INDIAN|PAKISTANI|BANGLADESHI)",
      "AMR":r"^(?:AMR|HISPANIC|LATINO|ADMIXED AMERICAN)"}
    for p,pat in rules.items(): out[z.str.contains(pat,regex=True,na=False)]=p
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sample",required=True); ap.add_argument("--ancestry",required=True)
    ap.add_argument("--out-keep",required=True); ap.add_argument("--out-panel",required=True)
    ap.add_argument("--max-individuals",type=int,default=20000)
    ap.add_argument("--full",action="store_true"); ap.add_argument("--anchors-per-pop",type=int,default=1000)
    ap.add_argument("--prob-min",type=float,default=.999); ap.add_argument("--seed",type=int,default=20260904)
    ap.add_argument("--custom-keep",default="")
    a=ap.parse_args(); rng=np.random.default_rng(a.seed)
    s=read_sample(a.sample)
    anc=read_table(a.ancestry)
    idc=next((c for c in ["eid","IID","#IID","sample","ID_2"] if c in anc.columns),None)
    if idc is None: raise SystemExit(f"No ID column in {a.ancestry}: {list(anc.columns)}")
    ac=next((c for c in ["ancestry","genetic_ancestry","predicted_ancestry","super_pop","population"] if c in anc.columns),None)
    if ac is None: raise SystemExit(f"No ancestry column in {a.ancestry}: {list(anc.columns)}")
    pc=next((c for c in ["ancestry_prob","ancestry_probability","probability","max_posterior","posterior_prob"] if c in anc.columns),None)
    anc=anc.rename(columns={idc:"eid"}); anc["eid"]=anc.eid.astype(str); anc["pop"]=norm_pop(anc[ac])
    anc["prob"]=pd.to_numeric(anc[pc],errors="coerce") if pc else np.nan
    # If no single probability column exists, use posterior of the assigned group.
    if not np.isfinite(anc["prob"]).any():
        vals=[]
        for _,r in anc.iterrows():
            c=f"posterior_{r['pop']}"; vals.append(pd.to_numeric(r.get(c,np.nan),errors="coerce"))
        anc["prob"]=vals
    d=s.merge(anc[["eid","pop","prob"]],on="eid",how="left"); d["pop"]=d["pop"].fillna("OTH")
    if a.custom_keep:
        k=pd.read_csv(a.custom_keep,sep=None,engine="python",comment="#",header=None,dtype=str)
        wanted=[x for x in k.iloc[:,0].astype(str) if x in set(d.eid)]
        rank={x:i for i,x in enumerate(wanted)}; d=d[d.eid.isin(rank)].copy(); d["rank"]=d.eid.map(rank); d=d.sort_values("rank")
    else:
        from itertools import zip_longest
        anchor_lists=[]
        for p in ["AFR","EAS","EUR","SAS"]:
            x=d[(d["pop"]==p) & ((d["prob"]>=a.prob_min) | d["prob"].isna())].copy()
            x=x.sort_values(["prob","source_order"],ascending=[False,True])
            anchor_lists.append(x.head(a.anchors_per_pop).eid.tolist())
        # Round-robin ordering ensures the initial ARG scaffold contains all groups.
        selected=[eid for row in zip_longest(*anchor_lists) for eid in row if eid is not None]
        selected=list(dict.fromkeys(selected))
        remaining=d[~d.eid.isin(selected)].copy()
        if a.full or a.max_individuals<=0:
            fill=remaining.sort_values("source_order").eid.tolist()
        else:
            n=max(0,a.max_individuals-len(selected))
            if n>len(remaining): n=len(remaining)
            # Stratified random fill so the pilot is not almost entirely EUR.
            fill=[]
            groups=[g for g in ["AFR","EAS","EUR","SAS","AMR","OTH"] if (remaining["pop"]==g).any()]
            if groups and n:
                alloc={g:min(len(remaining[remaining["pop"]==g]), n//len(groups)) for g in groups}
                for g,m in alloc.items():
                    ids=remaining.loc[remaining["pop"]==g,"eid"].to_numpy();
                    if m: fill.extend(rng.choice(ids,size=m,replace=False).tolist())
                remn=n-len(fill)
                pool=remaining[~remaining.eid.isin(fill)].eid.to_numpy()
                if remn: fill.extend(rng.choice(pool,size=min(remn,len(pool)),replace=False).tolist())
        selected=selected+fill
        order={x:i for i,x in enumerate(selected)}; d=d[d.eid.isin(order)].copy(); d["rank"]=d.eid.map(order); d=d.sort_values("rank")
    if not len(d): raise SystemExit("No ARG samples selected")
    Path(a.out_keep).parent.mkdir(parents=True,exist_ok=True)
    with open(a.out_keep,"w") as h:
        h.write("#FID\tIID\n"); h.writelines(f"{x}\t{x}\n" for x in d.eid)
    d[["eid","pop","prob","source_order"]].to_csv(a.out_panel,sep="\t",index=False)
    print(f"selected_individuals={len(d)} pop_counts={d['pop'].value_counts().to_dict()}")
if __name__=="__main__": main()
