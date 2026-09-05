#!/usr/bin/env python3
from __future__ import annotations
import argparse,gzip,json
from pathlib import Path

def op(path): return gzip.open(path,"rt") if str(path).endswith(".gz") else open(path)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--haps",required=True); ap.add_argument("--sample",required=True); ap.add_argument("--out",required=True); ap.add_argument("--max-check",type=int,default=2000); a=ap.parse_args()
    sl=Path(a.sample).read_text().splitlines(); n=max(0,len([x for x in sl[2:] if x.strip()])); exp=2*n
    nv=0; prev=-1; errors=[]; first=None; last=None
    with op(a.haps) as h:
      for line in h:
        if not line.strip(): continue
        z=line.split(); nv+=1
        # PLINK/Oxford HAPS usually has five metadata columns. Six-column variants are accepted.
        meta=5 if len(z)-5==exp else (6 if len(z)-6==exp else None)
        if meta is None:
            errors.append(f"line {nv}: columns={len(z)} not metadata+2N ({exp})");
            if len(errors)>10: break
            continue
        try: pos=int(float(z[2] if meta==5 else z[3]))
        except Exception: errors.append(f"line {nv}: bad position"); continue
        if pos<prev: errors.append(f"line {nv}: positions not sorted {prev}>{pos}")
        prev=pos; first=pos if first is None else first; last=pos
        if nv<=a.max_check:
            bad={x for x in z[meta:] if x not in {"0","1","?","NA"}}
            if bad: errors.append(f"line {nv}: invalid hap values {sorted(bad)[:5]}")
    q={"individuals":n,"haplotypes":exp,"variants":nv,"first_position":first,"last_position":last,"errors":errors}
    Path(a.out).parent.mkdir(parents=True,exist_ok=True); Path(a.out).write_text(json.dumps(q,indent=2)+"\n")
    if errors or nv==0 or n==0: raise SystemExit("HAPS validation failed: "+"; ".join(errors[:3]))
    print(json.dumps(q))
if __name__=="__main__": main()
