#!/usr/bin/env python3
"""Create TRACE tree-node -> sample/haplotype map.

Priority: metadata in tree sequence; fallback: VCF sample order when the tree has
exactly two sample nodes per VCF individual. The fallback is appropriate only
when the ARG converter preserved input haplotype order; the script records the
mapping mode so it is auditable.
"""
from __future__ import annotations
import argparse, json, subprocess
from pathlib import Path


def load_ts(path: Path):
    import tskit
    if str(path).endswith('.tsz'):
        import tszip
        return tszip.decompress(str(path))
    return tskit.load(str(path))


def meta(obj):
    x=getattr(obj,'metadata',None)
    if isinstance(x,dict): return x
    if isinstance(x,(bytes,bytearray)):
        try: return json.loads(x.decode())
        except Exception: return {}
    return {}


def vcf_samples(path: Path):
    out=subprocess.check_output(['bcftools','query','-l',str(path)],text=True)
    return [x for x in out.splitlines() if x]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--tree',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--vcf',type=Path)
    ap.add_argument('--panel',type=Path)
    args=ap.parse_args()
    ts=load_ts(args.tree)
    nodes=list(map(int,ts.samples()))
    vcf_ids=vcf_samples(args.vcf) if args.vcf else []
    known=set(vcf_ids)
    rows=[]; unresolved=[]; seen={}
    for j,node_id in enumerate(nodes):
        node=ts.node(node_id); cand=[]
        for obj in [node] + ([ts.individual(node.individual)] if node.individual >= 0 else []):
            m=meta(obj)
            for k in ('variant_data_sample_id','sample','sample_id','individual','individual_id','name','id'):
                if k in m: cand.append(str(m[k]))
        sample=next((x for x in cand if not known or x in known),None)
        if sample is None:
            unresolved.append((j,node_id))
        else:
            seen[sample]=seen.get(sample,0)+1
            rows.append([node_id,sample,seen[sample],'metadata'])
    if unresolved:
        if not vcf_ids or len(nodes) != 2*len(vcf_ids):
            raise SystemExit(f'Cannot map {len(unresolved)} sample nodes: tree nodes={len(nodes)}, VCF individuals={len(vcf_ids)}')
        rows=[]
        for j,node_id in enumerate(nodes):
            rows.append([node_id,vcf_ids[j//2],j%2+1,'vcf_order_fallback'])
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w') as h:
        h.write('tree_node_id\tsample\thaplotype\tmapping_mode\n')
        for r in rows: h.write('\t'.join(map(str,r))+'\n')
    print(f'Wrote {len(rows)} nodes to {args.out}; mode={rows[0][3] if rows else "none"}')

if __name__=='__main__': main()
