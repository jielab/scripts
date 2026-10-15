#!/usr/bin/env python3
"""Run twin confidence and molecular deviation from current cohort inputs."""
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import fcntl
import hashlib
import json
import os
import shutil
import sys
import tempfile
import uuid
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from abm_twin_core import json_write
from host_inputs import ROOT,host_config,prepare_host
from run_abm_twin import load_config
from research_io import write_tables,load_pgs,merge_native_scores
import run_deviation as deviation
from c1_abm_compare import run_comparison


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--Y',required=True);p.add_argument('--biom',choices=['prot','met'],required=True)
    p.add_argument('--analysis-root',type=Path,required=True);p.add_argument('--config',required=True)
    p.add_argument('--seed',type=int,default=2026);p.add_argument('--cores',type=int,default=4)
    p.add_argument('--replace',choices=['TRUE','FALSE'],default='FALSE');p.add_argument('--r-bin',default='Rscript')
    for flag in ['ukb-phe','group-file','group-col','outer-roster','end-date','diagnosis-col',
                 'shared-covariates','deviation-config','objective-file','pgs-root','pgs-mapping']:
        p.add_argument('--'+flag)
    p.add_argument('--host-args',default='');p.add_argument('--remerge-pgs',action='store_true')
    p.add_argument('--preflight',action='store_true')
    return p


def publish(stage,target):
    """Retain the previous complete run and publish only a successful replacement."""
    previous=None
    if target.exists():
        previous=target.parent/'_history'/(target.name+'_'+uuid.uuid4().hex[:12])
        previous.parent.mkdir(exist_ok=True);target.rename(previous)
    try:
        stage.rename(target)
    except BaseException:
        if previous:
            previous.rename(target)
        raise


def completed(target):
    marker=target/'COMPLETE.json'
    if not marker.is_file():
        return False
    saved=json.loads(marker.read_text())
    for name,digest in saved.get('files',{}).items():
        path=target/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            return False
    return bool(saved.get('files'))


def check_input_reuse(target, current):
    """Do not silently reuse a raw-assay fit under the cleaned-RDS contract."""
    contract=target/'host_contract.json'
    if not contract.is_file():
        return
    previous=json.loads(contract.read_text()).get('input_config',{})
    if not previous:
        return
    keys=('phe_file','omics_file','met_input','transform','covariates','residualize','categorical')
    changed=[key for key in keys if previous.get(key)!=getattr(current,key)]
    if changed:
        raise ValueError('ABM input configuration differs from the saved fit: '+','.join(changed)+
                         '. Run c1_abm --replace TRUE to fit the new inputs and archive the previous results: '+str(target))


def finish(stage, audit):
    contract=stage/'host_contract.json'
    existing=json.loads(contract.read_text()) if contract.is_file() else {}
    json_write(contract,{**existing,**audit})
    # Only small integrity anchors; do not rehash every large local array on resume.
    files={str(p.relative_to(stage)):hashlib.sha256(p.read_bytes()).hexdigest()
           for p in stage.rglob('*') if p.is_file() and (p.name in {'MODEL_FROZEN.json','table_manifest.json','host_contract.json'})}
    json_write(stage/'COMPLETE.json',dict(status='completed',files=files))


def pgs_path(args):
    if args.pgs_root:
        return Path(args.pgs_root)
    cached=args.analysis_root/'_pgs_fullprecision'/args.biom/(args.biom+'.pgs.fullprecision.gz')
    return cached if cached.is_file() else Path('/mnt/f/gwas')/args.biom/'pgs'/(args.biom+'.pgs.gz')


def run_deviation(args,data,stage,config,audit):
    source=pgs_path(args)
    if args.remerge_pgs:
        manifest=Path('/mnt/f/gwas')/args.biom/'pgs'/(args.biom+'.pgs.files.tsv')
        source=args.analysis_root/'_pgs_fullprecision'/args.biom/(args.biom+'.pgs.fullprecision.gz')
        if not source.exists():
            temp=source.with_name('_pending_'+uuid.uuid4().hex+'.gz')
            merge_native_scores(manifest,temp);temp.rename(source)
    if not source.exists():
        stage.mkdir()
        write_tables(stage/'aggregate',{'status':pd.DataFrame([dict(status='PGS_input_unavailable',
                     source=str(source),validation_status='not_validated_for_diagnosis')])},
                     'deviation',args.biom,args.Y,'unavailable')
        print('[LE8 twin] Deviation unavailable: missing PGS '+str(source),flush=True)
        finish(stage,audit);return
    names,G,pgs_audit=load_pgs(source,args.biom,data['eid'],data['features'],args.pgs_mapping)
    if len(names)<2:
        raise ValueError('Fewer than two exactly matched molecular PGSs')
    indices=[data['features'].index(n) for n in names]
    M=data['M'][:,indices];role=data['role'];dev=role!='test'
    eligible,y,reason=deviation.objective_endpoint(data['endpoint'],config)
    known=pd.to_numeric(data['endpoint'].known_cad_at_index,errors='coerce').to_numpy()
    stage.mkdir();local=stage/'_LOCAL_ONLY';local.mkdir(mode=0o700)
    fields=dict(M=M,G=G,C=data['C'],eid=data['eid'],group=data['group'],role=role,
                reference=known==0,y=y,endpoint_eligible=eligible)
    schema=dict(feature_names=np.asarray(names),covariate_names=np.asarray(data['covariate_names']))
    np.savez_compressed(local/'development.npz',**{k:v[dev] for k,v in fields.items()},**schema)
    np.savez_compressed(local/'test_X.npz',**{k:fields[k][~dev] for k in ['M','G','C','eid','group']},**schema)
    pd.DataFrame(dict(eid=data['eid'][~dev],y=y[~dev],endpoint_eligible=eligible[~dev],
                      exclusion_reason=reason[~dev])).to_csv(local/'test_truth.tsv.gz',sep='\t',index=False)
    pd.DataFrame(dict(eid=data['eid'],group=data['group'],role=role)).to_csv(local/'roster.tsv.gz',sep='\t',index=False)
    from dataclasses import asdict
    json_write(stage/'config.json',asdict(config))
    flow=pd.DataFrame(dict(role=role,reason=reason)).groupby(['role','reason']).size().reset_index(name='N')
    write_tables(stage/'aggregate',{'pgs_quality_audit':pgs_audit,'endpoint_flow':flow,
                 'support_status':pd.DataFrame([{'family_status':audit['family_status'],
                    'validation_status':'exploratory','objective_endpoint_N':int(eligible.sum())}])},
                 'deviation',args.biom,args.Y,'host')
    A=SimpleNamespace
    deviation.fit(A(data=local/'development.npz',config=stage/'config.json',output=stage/'fit'))
    model=stage/'fit/_LOCAL_ONLY/model.joblib'
    deviation.project(A(model=model,data=local/'test_X.npz',output=local/'prediction.joblib'))
    deviation.evaluate(A(model=model,predictions=local/'prediction.joblib',truth=local/'test_truth.tsv.gz',
                         output=str(stage/'results'),anchors='L_VLDL_TG.pct,Total_TG,ApoB,PCSK9,LPA,GDF15'))
    finish(stage,audit)


def main(args):
    scope=args.analysis_root/args.Y/args.biom/'c1_abm'
    cfg=replace(load_config(args.config),trait=args.Y,layer=args.biom,seed=args.seed,cores=args.cores)
    cfg.validate()
    dc=deviation.read_config(args.deviation_config or ROOT/f'config/deviation_{args.biom}.json')
    dc=replace(dc,trait=args.Y,layer=args.biom,seed=args.seed);dc.validate()
    # Explicit alternative endpoints have distinct directories and completion
    # identities; they cannot overwrite or resume an objective-plaque analysis.
    deviation_target=(scope/'molecular_deviation' if dc.task=='objective_plaque'
                      else scope/dc.task/'molecular_deviation')
    if args.Y!='cvd_cad':
        # Do not reinterpret another disease's prevalence as known CAD.
        deviation_supported=False
    else:
        deviation_supported=True
    if args.preflight:
        host,a=host_config(args)
        if not a.demo:
            for path in [a.phe_file,a.omics_file,*([a.met_map] if a.biom=='met' and a.met_input=='raw' else [])]:
                if not Path(path).is_file():
                    raise FileNotFoundError(path)
        for path in [args.objective_file,args.pgs_mapping,args.group_file,args.outer_roster]:
            if path and not Path(path).is_file():
                raise FileNotFoundError(path)
        if args.replace!='TRUE':
            for target in (scope/'abm_twin',deviation_target):
                if completed(target):
                    check_input_reuse(target,a)
        print('[LE8 twin] Preflight passed; family support='+('supplied' if a.group_col else 'unverified; individual IDs only')+
              '; objective endpoint='+('supplied' if args.objective_file else 'unavailable; reference calibration only'),flush=True)
        print(f'[LE8] DONE ABM input preflight: biom={a.biom}; omics={a.omics_file}; '
              f'met_input={a.met_input}; transform={a.transform}; covariates={a.covariates}; '
              f'extra_adjustment={a.residualize or "none"}',flush=True)
        return
    scope.mkdir(parents=True,exist_ok=True)
    lock_key=hashlib.sha256(str(scope.resolve()).encode()).hexdigest()[:20]
    with open('/tmp/le8-abm-'+lock_key+'.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        targets={'abm_twin':scope/'abm_twin','molecular_deviation':deviation_target}
        todo={name for name,target in targets.items() if args.replace=='TRUE' or not completed(target)}
        if args.replace!='TRUE' and set(targets)-todo:
            _,current=host_config(args)
            for name in set(targets)-todo:
                check_input_reuse(targets[name],current)
        if not todo:
            print('[LE8 twin] SKIP: both twin fits completed; refreshing C1 comparison if needed',flush=True)
            with tempfile.TemporaryDirectory(prefix='_abm_staging_',dir=scope) as work:
                run_comparison(args,scope,Path(work),publish,completed,finish)
                refresh_performance(args,scope,Path(work))
            return
        work=Path(tempfile.mkdtemp(prefix='_abm_staging_',dir=scope))
        try:
            abm,data,audit=prepare_host(args,work)
            unit='family' if audit['independent_family_support_verified'] else 'individual_unverified'
            cfg=replace(cfg,grouping_unit=unit);dc=replace(dc,grouping_unit=unit)
            if 'abm_twin' in todo:
                print('[LE8 twin] START twin fit/project/evaluate',flush=True)
                from host_adapter import run_prepared_numeric_cohort
                stage=work/'abm_twin'
                result=run_prepared_numeric_cohort(**abm,output=stage,config=cfg)
                write_tables(stage/'aggregate',{'support_status':pd.DataFrame([{
                    'family_status':audit['family_status'],'independent_family_support_verified':audit['independent_family_support_verified'],
                    'validation_status':'exploratory','split_source':audit['split_source']}])},'ABM twin',args.biom,args.Y,'host')
                finish(stage,audit);publish(stage,scope/'abm_twin')
                print('[LE8 twin] DONE twin',flush=True)
            del abm
            if 'molecular_deviation' in todo:
                stage=work/'molecular_deviation'
                if deviation_supported:
                    print('[LE8 twin] START molecular deviation',flush=True)
                    run_deviation(args,data,stage,dc,audit)
                else:
                    stage.mkdir()
                    write_tables(stage/'aggregate',{'status':pd.DataFrame([dict(status='not_applicable_non_CAD_trait',trait=args.Y)])},
                                 'deviation',args.biom,args.Y,'not_applicable')
                    finish(stage,audit)
                deviation_target.parent.mkdir(parents=True,exist_ok=True)
                publish(stage,deviation_target)
                print('[LE8 twin] DONE deviation',flush=True)
            run_comparison(args,scope,work,publish,completed,finish)
            refresh_performance(args,scope,work)
        except BaseException as exc:
            # The public console recognizes this prefix before the outer
            # CalledProcessError hides the actual cause behind a long command.
            print(f'[LE8] FAIL c1_abm Y={args.Y} biom={args.biom}: {type(exc).__name__}: {exc}',
                  file=sys.stderr, flush=True)
            print('[LE8 twin] Failed; earlier published results retained. Local diagnostics: '+str(work),
                  file=sys.stderr, flush=True)
            raise
        else:
            shutil.rmtree(work)


def refresh_performance(args,scope,work):
    from performance import collect, write_snapshot
    stage=work/'performance'
    tables=collect(args.analysis_root,args.analysis_root,args.Y,args.biom,'c1_abm')
    write_snapshot(tables,stage,f'{args.Y} / {args.biom} | through c1_abm')
    publish(stage,scope/'performance')


if __name__=='__main__':
    os.umask(0o077)
    args=parser().parse_args()
    with threadpool_limits(args.cores):
        main(args)
