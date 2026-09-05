#!/usr/bin/env python3
from __future__ import annotations
import argparse,gzip,json
from pathlib import Path

def op(path): return gzip.open(path,"rt") if str(path).endswith(".gz") else open(path)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--haps",required=True); ap.add_argument("--sample",required=True); ap.add_argument("--out",required=True); ap.add_argument("--max-check",type=int,default=0); ap.add_argument('--keep'); ap.add_argument('--map'); a=ap.parse_args()
    # Compatibility with older callers; new preparation validates every genotype.
    sl=Path(a.sample).read_text().splitlines(); n=max(0,len([x for x in sl[2:] if x.strip()])); exp=2*n
    nv=0; prev=-1; errors=[]; first=None; last=None
    ids=[x.split()[1] for x in sl[2:] if x.strip()]
    if len(set(ids))!=n: errors.append('duplicate sample IID')
    if a.keep:
        wanted=[x.split()[-1] for x in Path(a.keep).read_text().splitlines() if x.strip() and not x.startswith('#')]
        if ids!=wanted: errors.append('sample identities/order differ from selected panel')
    maps=open(a.map) if a.map else None
    prev_cm=-float('inf')
    with op(a.haps) as h:
      for line in h:
        if not line.strip(): continue
        z=line.split(); nv+=1
        # PLINK/Oxford HAPS usually has five metadata columns. Six-column variants are accepted.
        meta=5 if len(z)-5==exp else None
        if meta is None:
            errors.append(f"line {nv}: columns={len(z)} not metadata+2N ({exp})");
            if len(errors)>10: break
            continue
        try: pos=int(float(z[2] if meta==5 else z[3]))
        except Exception: errors.append(f"line {nv}: bad position"); continue
        if pos<=prev: errors.append(f"line {nv}: positions not strictly increasing {prev}>={pos}")
        prev=pos; first=pos if first is None else first; last=pos
        if a.max_check==0 or nv<=a.max_check:
            bad={x for x in z[meta:] if x not in {"0","1"}}
            if bad: errors.append(f"line {nv}: invalid hap values {sorted(bad)[:5]}")
        if len(z[3])!=1 or len(z[4])!=1 or z[3] not in 'ACGT' or z[4] not in 'ACGT' or z[3]==z[4]: errors.append(f'line {nv}: invalid alleles')
        if maps:
            m=maps.readline().split()
            if len(m)!=4 or m[0].removeprefix('chr')!=z[0].removeprefix('chr') or m[1]!=z[1] or int(m[3])!=pos: errors.append(f'line {nv}: map/HAPS mismatch')
            elif float(m[2])<=prev_cm: errors.append(f'line {nv}: genetic positions not strictly increasing')
            else: prev_cm=float(m[2])
        if len(errors)>10: break
    if maps:
        if maps.readline(): errors.append('extra map records')
        maps.close()
    q={"individuals":n,"haplotypes":exp,"variants":nv,"first_position":first,"last_position":last,"errors":errors}
    Path(a.out).parent.mkdir(parents=True,exist_ok=True); Path(a.out).write_text(json.dumps(q,indent=2)+"\n")
    if errors or nv==0 or n==0: raise SystemExit("HAPS validation failed: "+"; ".join(errors[:3]))
    if a.keep:
        # Inference identifies samples by IID; Oxford requires identical ID_1/ID_2.
        # Only normalize after proving the IID order matches the selected panel.
        rows=[x.split() for x in sl[2:] if x.strip()]
        if any(r[0]!=r[1] for r in rows):
            for r in rows: r[0]=r[1]
            Path(a.sample).write_text('\n'.join(sl[:2])+ '\n' + ''.join(' '.join(r)+'\n' for r in rows))
    print(json.dumps(q))
if __name__=="__main__": main()
