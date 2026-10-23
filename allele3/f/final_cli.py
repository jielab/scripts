#!/usr/bin/env python3
"""FINAL utilities use the same modules as the integrated scanner/pipeline."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
if __package__: sys.path.insert(0,str(Path(__file__).resolve().parent))
from ldblocks import validate_directory,load_block_sets,hotspot_consensus
from erv import parse_joined_rmsk,write_erv_bed,prepare_joined_cache
from architecture import loco_expected,write_summary
from ukb_mfi import scan_mfi_dir
from universes import build_universes,select_u3
from core import dump_json

def main(argv=None):
    p=argparse.ArgumentParser(prog='allele3-final'); sp=p.add_subparsers(dest='cmd',required=True)
    q=sp.add_parser('validate-ld'); q.add_argument('--dir',required=True); q.add_argument('--build',choices=['37','38'],default='38'); q.add_argument('--out',required=True)
    q=sp.add_parser('prepare-erv'); q.add_argument('--joined-rmsk',required=True); q.add_argument('--out',required=True); q.add_argument('--cache-dir'); q.add_argument('--chrom')
    q=sp.add_parser('summarize-architecture'); q.add_argument('--global-strata',required=True); q.add_argument('--ld-strata'); q.add_argument('--erv-strata'); q.add_argument('--out-dir',required=True); q.add_argument('--min-sites',type=int,default=10000); q.add_argument('--min-training-sites',type=int,default=1)
    q=sp.add_parser('ukb-mfi'); q.add_argument('--mfi-dir',required=True); q.add_argument('--targets',required=True); q.add_argument('--out',required=True)
    q=sp.add_parser('build-universes'); q.add_argument('--intersection',required=True); q.add_argument('--out',required=True); q.add_argument('--ld-block-dir',required=True); q.add_argument('--u1-sites',type=int,default=20000); q.add_argument('--controls-per-case',type=int,default=4); q.add_argument('--seed',type=int,default=20261022)
    q=sp.add_parser('select-u3'); q.add_argument('--universes',required=True); q.add_argument('--out',required=True); q.add_argument('--u3-sites',type=int,default=3000); q.add_argument('--seed',type=int,default=20261022)
    q=sp.add_parser('final-report'); q.add_argument('--modern-dir',required=True); q.add_argument('--out',required=True); q.add_argument('--metag-census'); q.add_argument('--no-plots',action='store_true')
    a=p.parse_args(argv)
    if a.cmd=='validate-ld':
        load_block_sets(a.dir,a.build); result=validate_directory(a.dir,a.build); dump_json(a.out,result)
    elif a.cmd=='prepare-erv':
        if a.cache_dir: result=prepare_joined_cache(a.joined_rmsk,a.cache_dir)
        else:
            e=parse_joined_rmsk(a.joined_rmsk,chromosome=a.chrom); write_erv_bed(e,a.out)
            from collections import Counter
            result=dict(n_elements=len(e),classes=dict(Counter(x.rep_class for x in e)),tiers=dict(Counter(x.tier for x in e)),source=a.joined_rmsk)
        dump_json(str(a.out)+'.manifest.json',result)
    elif a.cmd=='summarize-architecture':
        od=Path(a.out_dir);od.mkdir(parents=True,exist_ok=True);result={}
        for kind,path in [('ld',a.ld_strata),('erv',a.erv_strata)]:
            if path:
                r=loco_expected(a.global_strata,path,kind,a.min_training_sites)
                write_summary(r,od/('ld_block_enrichment.tsv' if kind=='ld' else 'erv_enrichment.tsv'));result[kind]=len(r)
                if kind=='ld': write_summary(hotspot_consensus(r,a.min_sites),od/'ld_hotspot_consensus.tsv')
    elif a.cmd=='ukb-mfi': result=scan_mfi_dir(a.mfi_dir,a.targets,a.out)
    elif a.cmd=='build-universes': result=build_universes(a.intersection,a.out,ld_block_dir=a.ld_block_dir,u1_sites=a.u1_sites,controls_per_case=a.controls_per_case,seed=a.seed)
    elif a.cmd=='select-u3': result=select_u3(a.universes,a.out,a.u3_sites,a.seed)
    elif a.cmd=='final-report':
        from final_figures import write_local_review
        result=write_local_review(a.modern_dir,a.out,a.metag_census,not a.no_plots)
    print(json.dumps(result,indent=2,ensure_ascii=False));return 0
if __name__=='__main__':
    try: raise SystemExit(main())
    except (ValueError,OSError) as e: print(f'allele3 FINAL: {e}',file=sys.stderr);raise SystemExit(2)
