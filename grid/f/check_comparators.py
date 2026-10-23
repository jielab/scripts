#!/usr/bin/env python3
"""Read-only sanity check for GRID published-method comparators.

No model is trained or modified.  The purpose is to verify that saved PRSformer
results can be reused and that comparator artefacts are present before a new ABM run.
"""
from __future__ import annotations
import argparse, hashlib, json, subprocess
from pathlib import Path

PRSFORMER_OFFICIAL_COMMIT='7dea4e1bb27975885c937f2be82f9243bc82bda7'

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--analysis-root',default='/mnt/d/analysis/grid')
    a.add_argument('--prsformer-upstream',default='/mnt/f/software/PRSformer')
    a.add_argument('--strict',action='store_true')
    x=a.parse_args();root=Path(x.analysis_root);up=Path(x.prsformer_upstream)
    out={'mode':'read_only_sanity_check','comparators_are_not_GRID_innovation':True,'PRSformer':{},'PRS-CSx':{},'DiscoDivas':{}}
    p=root/'prsformer/run/evaluation.json';e={}
    if p.is_file():
        e=json.loads(p.read_text());out['PRSformer'].update(saved_evaluation=True,evaluation_sha256=sha256(p),preparation_signature=e.get('preparation_signature'),model=e.get('model',{}))
        metrics=root/'prsformer/run/metrics.tsv'
        out['PRSformer']['metrics_present']=metrics.is_file()
    else:out['PRSformer']['saved_evaluation']=False
    if up.is_dir():
        try:
            head=subprocess.run(['git','-C',str(up),'rev-parse','HEAD'],capture_output=True,text=True,check=True,timeout=5).stdout.strip()
        except Exception:head=None
        out['PRSformer'].update(upstream_checkout=str(up),upstream_head=head,expected_official_commit=PRSFORMER_OFFICIAL_COMMIT,official_commit_match=head==PRSFORMER_OFFICIAL_COMMIT)
    else:out['PRSformer'].update(upstream_checkout=str(up),upstream_head=None,expected_official_commit=PRSFORMER_OFFICIAL_COMMIT,official_commit_match=None)
    for trait in ('height','ldl','t2dm'):
        out['PRS-CSx'][trait]=(root/'csx'/trait).is_dir()
        out['DiscoDivas'][trait]=(root/'disco'/trait).is_dir()
    out['PRSformer']['reuse_recommendation']='reuse saved PRSformer; FINAL ABM will additionally require exact outer-test-roster match before fitting/comparison' if out['PRSformer'].get('saved_evaluation') else 'PRSformer result unavailable'
    blockers=[]
    if not out['PRSformer'].get('saved_evaluation'):blockers.append('missing_saved_PRSformer_evaluation')
    if out['PRSformer'].get('official_commit_match') is not True:blockers.append('PRSformer_upstream_commit_unverified_or_changed')
    if e:
        checks={}
        source=e.get('model',{}).get('upstream_source_sha256',{})
        for name in ('utils','modules','model'):
            path=up/'src'/f'{name}.py'
            checks['official_source_'+name]=path.is_file() and sha256(path)==source.get(name)
        path=root/'prsformer/run/model.pt';identity=e.get('model',{})
        checks['saved_model']=path.is_file() and path.stat().st_size==identity.get('size') and path.stat().st_mtime_ns==identity.get('mtime_ns')
        for name in ('metrics.tsv','history.tsv','test_predictions.tsv.gz'):
            path=root/'prsformer/run'/name
            checks[name]=path.is_file() and sha256(path)==e.get('outputs_sha256',{}).get(name)
        path=root/'prsformer/cache/prepare.json'
        preparation=json.loads(path.read_text()) if path.is_file() else {}
        checks['preparation_signature']=bool(preparation.get('signature')) and preparation.get('signature')==e.get('preparation_signature')
        for name,key in (('data.tsv.gz','data_sha256'),('variants.tsv.gz','variants_sha256')):
            path=root/'prsformer/cache'/name
            checks[name]=path.is_file() and sha256(path)==e.get('inputs',{}).get(key)
        path=root/'prsformer/cache/genotypes.npy';identity=e.get('inputs',{}).get('genotypes',{})
        checks['genotype_input']=path.is_file() and path.stat().st_size==identity.get('size') and path.stat().st_mtime_ns==identity.get('mtime_ns') and preparation.get('shape')==identity.get('shape')
        if not path.exists():
            from prsformer_cleanup import validate_release
            try:
                validate_release(path.parent, root/'prsformer/run/evaluation.json')
                checks['genotype_input'] = preparation.get('shape') == identity.get('shape')
            except (OSError, ValueError, KeyError):
                pass
        out['PRSformer']['receipt_checks']=checks
        blockers.extend('PRSformer_receipt_mismatch:'+name for name,ok in checks.items() if not ok)
    out['blockers']=blockers
    print(json.dumps(out,indent=2))
    if x.strict and blockers:raise SystemExit(2)
if __name__=='__main__':main()
