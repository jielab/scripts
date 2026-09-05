#!/usr/bin/env python3
"""Normalize heterogeneous GWAS summary statistics to PRS-CSx input."""
from __future__ import annotations
import argparse,csv,gzip,json,math,re
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import norm
ALIASES={
 "SNP":["SNP","RSID","RS_ID","MARKERNAME","VARIANT_ID","ID","RS_NUMBER"],
 "CHR":["CHR","CHROM","CHROMOSOME","CHROMSOME","CHR_ID"], "BP":["BP","POS","POSITION","BASE_PAIR_LOCATION"],
 "A1":["A1","EA","EFFECT_ALLELE","EFFECTALLELE","ALT","ALLELE1","TESTED_ALLELE","CODED_ALLELE"],
 "A2":["A2","NEA","OTHER_ALLELE","NON_EFFECT_ALLELE","NONEFFECTALLELE","REF","ALLELE0","REFERENCE_ALLELE"],
 "BETA":["BETA","EFFECT","EFFECT_SIZE","ES","LOG_ODDS","B"], "OR":["OR","ODDS_RATIO","ODDSRATIO"],
 "SE":["SE","STDERR","STANDARD_ERROR","BETA_SE"], "P":["P","PVAL","PVALUE","P_VALUE","P-VALUE"],
 "N":["N","N_TOTAL","TOTAL_N","OBS_CT","N_SAMPLES","NEFF","N_EFF"],
 "NCASE":["N_CASE","NCASE","NCASES","CASES","N_CASES"], "NCTRL":["N_CONTROL","N_CONTROLS","NCTRL","NCONTROLS","CONTROLS"],
 "EAF":["EAF","AF","POOLED_ALT_AF","EFFECT_ALLELE_FREQUENCY","A1FREQ","FREQ1","ALT_FREQ","MAF"]}
def canon(s): return re.sub(r"[^A-Z0-9]+","_",str(s).upper()).strip("_")
def choose(cols,key):
    d={canon(c):c for c in cols}
    return next((d[canon(x)] for x in ALIASES[key] if canon(x) in d),None)
def read_snpinfo(path):
    d=pd.read_csv(path,sep=r"\s+",engine="python",dtype=str)
    sc=choose(d.columns,"SNP") or d.columns[1 if len(d.columns)>1 else 0]
    cc=choose(d.columns,"CHR"); bc=choose(d.columns,"BP")
    snps=set(d[sc].dropna().astype(str)); locus={}
    if cc and bc:
      for s,c,b in zip(d[sc],d[cc],d[bc]):
        try: locus[(str(c).replace("chr","").replace("CHR",""),int(float(b)))]=str(s)
        except: pass
    return snps,locus
def sep_for(path):
    op=gzip.open if str(path).endswith(".gz") else open
    with op(path,"rt",errors="replace") as h:
      line=next((x for x in h if x.strip() and not x.startswith("##")),"")
    if "\t" in line:return "\t"
    if "," in line:return ","
    return r"\s+"
def main():
  ap=argparse.ArgumentParser(); ap.add_argument("--input",required=True); ap.add_argument("--output",required=True); ap.add_argument("--metadata",required=True); ap.add_argument("--snpinfo",required=True); ap.add_argument("--trait",required=True); ap.add_argument("--pop",required=True); ap.add_argument("--chunk",type=int,default=500000); a=ap.parse_args()
  snps,locus=read_snpinfo(a.snpinfo); sep=sep_for(a.input); engine="python" if sep==r"\s+" else "c"
  Path(a.output).parent.mkdir(parents=True,exist_ok=True); seen=set(); nvals=[]; stats={"input_rows":0,"kept_rows":0,"bad_allele":0,"ambiguous":0,"missing_effect":0,"not_hm3":0,"duplicates":0}
  first=True
  try: it=pd.read_csv(a.input,sep=sep,engine=engine,compression="infer",chunksize=a.chunk,comment="#",low_memory=False)
  except TypeError: it=pd.read_csv(a.input,sep=sep,engine=engine,compression="infer",chunksize=a.chunk,comment="#")
  for d in it:
    stats["input_rows"]+=len(d); cols={k:choose(d.columns,k) for k in ALIASES}
    if not cols["A1"] or not cols["A2"]: raise SystemExit(f"Allele columns not found in {a.input}: {list(d.columns)}")
    chrom=d[cols["CHR"]].astype(str).str.replace(r"^chr","",case=False,regex=True) if cols["CHR"] else pd.Series("",index=d.index)
    bp=pd.to_numeric(d[cols["BP"]],errors="coerce") if cols["BP"] else pd.Series(np.nan,index=d.index)
    sid=d[cols["SNP"]].astype(str) if cols["SNP"] else pd.Series("",index=d.index)
    # Recover rsID from SNPINFO by build-37 position when necessary.
    sid=[s if s in snps else locus.get((str(c),int(b)),s) if np.isfinite(b) else s for s,c,b in zip(sid,chrom,bp)]
    a1=d[cols["A1"]].astype(str).str.upper(); a2=d[cols["A2"]].astype(str).str.upper()
    beta=pd.to_numeric(d[cols["BETA"]],errors="coerce") if cols["BETA"] else (np.log(pd.to_numeric(d[cols["OR"]],errors="coerce")) if cols["OR"] else pd.Series(np.nan,index=d.index))
    se=pd.to_numeric(d[cols["SE"]],errors="coerce") if cols["SE"] else pd.Series(np.nan,index=d.index)
    pv=pd.to_numeric(d[cols["P"]],errors="coerce") if cols["P"] else pd.Series(np.nan,index=d.index)
    need=se.isna()&beta.notna()&pv.notna()&(pv>0)&(pv<=1)&(beta!=0); se.loc[need]=np.abs(beta.loc[need])/norm.isf(np.clip(pv.loc[need]/2,1e-323,.5))
    if cols["N"]: n=pd.to_numeric(d[cols["N"]],errors="coerce")
    elif cols["NCASE"] and cols["NCTRL"]:
      ca=pd.to_numeric(d[cols["NCASE"]],errors="coerce"); co=pd.to_numeric(d[cols["NCTRL"]],errors="coerce"); n=4/(1/ca+1/co)
    else:n=pd.Series(np.nan,index=d.index)
    eaf=pd.to_numeric(d[cols["EAF"]],errors="coerce") if cols.get("EAF") else pd.Series(np.nan,index=d.index)
    out=pd.DataFrame({"SNP":sid,"A1":a1,"A2":a2,"BETA":beta,"SE":se,"P":pv,"N":n,"EAF":eaf,"CHR":chrom,"BP":bp})
    valid=out.A1.isin(list("ACGT"))&out.A2.isin(list("ACGT"))&(out.A1!=out.A2); stats["bad_allele"]+=int((~valid).sum()); out=out[valid]
    amb=(out.A1+out.A2).isin(["AT","TA","CG","GC"]); stats["ambiguous"]+=int(amb.sum()); out=out[~amb]
    eff=out.BETA.notna()&out.SE.notna()&(out.SE>0); stats["missing_effect"]+=int((~eff).sum()); out=out[eff]
    hm=out.SNP.isin(snps); stats["not_hm3"]+=int((~hm).sum()); out=out[hm]
    dup=out.SNP.isin(seen); stats["duplicates"]+=int(dup.sum()); out=out[~dup]; seen.update(out.SNP.tolist())
    if len(out):
      nvals.extend(out.N.dropna().to_numpy().tolist()); out.to_csv(a.output,sep="\t",index=False,mode="wt" if first else "at",header=first,compression="gzip",float_format="%.10g"); first=False; stats["kept_rows"]+=len(out)
  if first: raise SystemExit(f"No usable HapMap3 variants in {a.input}")
  nmed=float(np.nanmedian(nvals)) if nvals else None
  meta={**stats,"trait":a.trait,"pop":a.pop,"input":str(Path(a.input).resolve()),"output":str(Path(a.output).resolve()),"n_gwas_median":nmed}
  Path(a.metadata).write_text(json.dumps(meta,indent=2)+"\n")
  print(json.dumps(meta))
if __name__=="__main__":main()
