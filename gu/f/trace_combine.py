#!/usr/bin/env python3
"""Normalize TRACE calls and optionally apply requested loci post hoc.

TRACE inference should normally use chromosome-scale/multi-chromosome context.  The
loci BED is therefore a final segment-selection map, not an HMM training region.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path
import pandas as pd

def first(cols, names):
    low = {str(c).lower(): c for c in cols}
    return next((low[n.lower()] for n in names if n.lower() in low), None)

def norm_one(path: Path, node: str) -> pd.DataFrame:
    try: d = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError: return pd.DataFrame()
    if d.empty: return d
    c=first(d.columns,["chrom","chromosome","chr"]); s=first(d.columns,["start","start_bp","left","begin"]); e=first(d.columns,["end","end_bp","right","stop"])
    p=first(d.columns,["mean_posterior","posterior","posterior_mean","prob","probability"]); cm=first(d.columns,["length(cM)","length_cM","length_cm","genetic_length_cm"])
    if c is None or s is None or e is None: raise SystemExit(f"Unrecognized TRACE columns in {path}: {list(d.columns)}")
    out=pd.DataFrame({"tree_node_id":str(node),"chr":d[c].astype(str).str.replace("chr","",regex=False),"start":pd.to_numeric(d[s],errors="coerce"),"end":pd.to_numeric(d[e],errors="coerce")})
    out["posterior"]=pd.to_numeric(d[p],errors="coerce") if p else pd.NA; out["length_cM"]=pd.to_numeric(d[cm],errors="coerce") if cm else pd.NA
    return out[(out.start.notna())&(out.end>out.start)].copy()

def clip_loci(data: pd.DataFrame, path: Path) -> pd.DataFrame:
    loci=pd.read_csv(path,sep="\t",header=None,names=["chr","locus_start","locus_end","locus_id"],dtype={"chr":str}); loci["chr"]=loci.chr.str.replace("chr","",regex=False)
    joined=data.merge(loci,on="chr",how="inner"); joined["start"]=joined[["start","locus_start"]].max(axis=1); joined["end"]=joined[["end","locus_end"]].min(axis=1)
    return joined[joined.end>joined.start].drop(columns=["locus_start","locus_end"])

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",required=True,type=Path); ap.add_argument("--sample-map",required=True,type=Path); ap.add_argument("--build",default="GRCh37"); ap.add_argument("--loci",type=Path); args=ap.parse_args()
    sample_map=pd.read_csv(args.sample_map,sep="\t",dtype={"tree_node_id":str}); pieces=[]
    for path in sorted((args.root/"calls").glob("hap*.summary.txt")):
        match=re.match(r"hap(.+)\.summary\.txt$",path.name)
        if match:
            part=norm_one(path,match.group(1)); pieces.extend([part] if not part.empty else [])
    cols=["tree_node_id","chr","start","end","posterior","length_cM"]; hap=pd.concat(pieces,ignore_index=True) if pieces else pd.DataFrame(columns=cols)
    if args.loci and not hap.empty: hap=clip_loci(hap,args.loci)
    hap["length_bp"]=hap.end-hap.start; hap=hap.merge(sample_map[["tree_node_id","sample","haplotype"]],on="tree_node_id",how="left")
    if not hap.empty and hap["sample"].isna().any():
        missing=sorted(hap.loc[hap["sample"].isna(),"tree_node_id"].astype(str).unique())
        raise SystemExit(f"TRACE sample map is missing tree nodes: {','.join(missing[:10])}")
    hap["method"]="trace"; hap["source"]="ghost_or_unknown_archaic"; hap["genome_build"]=args.build
    order=["sample","haplotype","tree_node_id","method","source","chr","start","end","length_bp","posterior","length_cM","genome_build"]
    if "locus_id" in hap: order.insert(5,"locus_id")
    final=args.root/"final"; final.mkdir(parents=True,exist_ok=True); hap[order].to_csv(final/"trace_haplotype_segments.tsv.gz",sep="\t",index=False,compression="gzip")
    print(f"TRACE segments: {len(hap)} -> {final/'trace_haplotype_segments.tsv.gz'}")
if __name__=="__main__": main()
