#!/usr/bin/env python3
"""Create a TRACE tree-node -> sample/haplotype map from tree metadata."""
from __future__ import annotations
import argparse, json
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


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--tree',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path)
    args=ap.parse_args()
    ts=load_ts(args.tree)
    nodes=list(map(int,ts.samples()))
    rows=[]; unresolved=[]; seen={}
    for j,node_id in enumerate(nodes):
        node=ts.node(node_id); cand=[]
        for obj in [node] + ([ts.individual(node.individual)] if node.individual >= 0 else []):
            m=meta(obj)
            for k in ('variant_data_sample_id','sample','sample_id','individual','individual_id','name','id'):
                if k in m: cand.append(str(m[k]))
        sample=next((x for x in cand if x),None)
        if sample is None:
            unresolved.append((j,node_id))
        else:
            seen[sample]=seen.get(sample,0)+1
            rows.append([node_id,sample,seen[sample],'metadata'])
    if unresolved:
        examples=','.join(str(node_id) for _,node_id in unresolved[:10])
        raise SystemExit(
            f'Cannot map {len(unresolved)} sample nodes from ARG metadata; '
            f'example node IDs: {examples}'
        )
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w') as h:
        h.write('tree_node_id\tsample\thaplotype\tmapping_mode\n')
        for r in rows: h.write('\t'.join(map(str,r))+'\n')
    print(f'Wrote {len(rows)} nodes to {args.out}; mode={rows[0][3] if rows else "none"}')

if __name__=='__main__': main()
