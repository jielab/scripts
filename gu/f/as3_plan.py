#!/usr/bin/env python3
"""Create non-overlapping AS3 target/reference panels.

If target samples overlap the African reference samples, cross-fitting keeps every
target sample while ensuring that a sample is never used as its own reference.
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path


def ids(path: Path) -> list[str]:
    out=[]
    for line in path.read_text().splitlines():
        x=line.strip().split()[0] if line.strip() else ""
        if x and not x.startswith("#"): out.append(x)
    if len(out)!=len(set(out)): raise SystemExit(f"Duplicate IDs in {path}")
    return out


def write(path: Path, values: list[str]) -> None:
    path.write_text("\n".join(values)+("\n" if values else ""))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--target', required=True, type=Path)
    p.add_argument('--reference', required=True, type=Path)
    p.add_argument('--outdir', required=True, type=Path)
    p.add_argument('--policy', choices=['error','crossfit'], default='crossfit')
    p.add_argument('--folds', type=int, default=5)
    a=p.parse_args()
    if a.folds < 2: raise SystemExit('--folds must be >=2')
    t=ids(a.target); r=ids(a.reference)
    ts=set(t); rs=set(r); overlap=sorted(ts&rs)
    a.outdir.mkdir(parents=True, exist_ok=True)
    rows=[]
    if overlap and a.policy=='error':
        raise SystemExit(f"Target/reference overlap contains {len(overlap)} samples; examples: {','.join(overlap[:10])}")
    nonoverlap=[x for x in t if x not in rs]
    if nonoverlap:
        tp=a.outdir/'main.target.txt'; rp=a.outdir/'main.reference.txt'
        write(tp,nonoverlap); write(rp,r)
        rows.append(('main',tp,rp,len(nonoverlap),len(r)))
    if overlap:
        fold_members={k:[] for k in range(a.folds)}
        # Stable hash ordering plus round-robin assignment keeps folds balanced.
        ordered=sorted(overlap, key=lambda x: hashlib.sha256(x.encode()).hexdigest())
        for i, x in enumerate(ordered):
            fold_members[i % a.folds].append(x)
        used=0
        for k in range(a.folds):
            target=sorted(fold_members[k])
            if not target: continue
            ref=[x for x in r if x not in set(target)]
            if len(ref)<2: raise SystemExit(f"Cross-fit fold {k+1} leaves fewer than 2 African reference samples")
            pid=f'overlap_fold{k+1}'
            tp=a.outdir/f'{pid}.target.txt'; rp=a.outdir/f'{pid}.reference.txt'
            write(tp,target); write(rp,ref)
            rows.append((pid,tp,rp,len(target),len(ref))); used+=len(target)
        if used!=len(overlap): raise SystemExit('Internal cross-fit accounting error')
    if not rows:
        tp=a.outdir/'all.target.txt'; rp=a.outdir/'all.reference.txt'
        write(tp,t); write(rp,r); rows.append(('all',tp,rp,len(t),len(r)))
    with (a.outdir/'panels.tsv').open('w') as h:
        h.write('panel_id\ttarget_list\tafr_reference_list\tn_target\tn_afr_reference\n')
        for row in rows: h.write('\t'.join(map(str,row))+'\n')
    write(a.outdir/'overlap.txt', overlap)
    print(f"panels={len(rows)} target={len(t)} reference={len(r)} overlap={len(overlap)}")

if __name__=='__main__': main()
