#!/usr/bin/env python3
"""GRID twin: run/fit/project/evaluate/demo with separated outcome access."""
from __future__ import annotations
import argparse,json,os,sys,uuid,datetime
from pathlib import Path
from dataclasses import asdict,replace
import joblib,numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from twin_core import Config,TwinModel,state_hash,id_hash,ROLE_NAMES
from twin_evaluate import evaluate
from twin_io import json_write,read_table,identified,write_aggregate,sha256
from twin_host import load_source,numeric_cohort,calibrate_disco,verified_prsformer

DEFAULT_CONFIG=Path(__file__).resolve().parents[1]/'config/twin.json'


def config_file(path,trait=None):
    d=json.loads(Path(path).read_text()) if path else {}
    for name in ('role_fractions','coverages','ridge_grid','logistic_c_grid'):
        if name in d:d[name]=tuple(d[name])
    cfg=Config(**d)
    if trait:cfg=replace(cfg,trait=trait,trait_type='binary' if trait=='t2dm' else 'continuous')
    cfg.validate();return cfg


def npz(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}


def export_projection(p,path):
    d=pd.DataFrame(dict(eid=p['eid'],group=p['group'],selected=p['selected'],reason=p['reason'],**p['diagnostics']))
    for name,v in p['predictions'].items():d['prediction::'+name]=v
    d['expected_squared_error']=p['expected_squared_error'];d['model_hash']=p['frozen_hash']
    d.to_csv(path,sep='\t',index=False,float_format='%.17g',compression='gzip')


def save_fit(model,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);private=out/'_LOCAL_ONLY';private.mkdir(mode=0o700,exist_ok=True)
    joblib.dump(model,private/'model.joblib',compress=3)
    digest=state_hash(model)
    json_write(out/'MODEL_FROZEN.json',dict(model_hash=digest,config=asdict(model.config),test_outcomes_used=False,clinical_release=False))
    pd.DataFrame(model.donors).to_csv(private/'donors.tsv.gz',sep='\t',index=False,float_format='%.17g')
    feature_rows=[]
    for pos,original in enumerate(model.prep.keep):
        is_genetic=original in set(model.prep.profile_keep)
        gp=list(model.prep.profile_keep).index(original) if is_genetic else None
        feature_rows.append(dict(feature=model.features[original],role='genetic_profile' if is_genetic else 'clinical_covariate',training_mean=model.prep.mean[pos],training_sd=model.prep.scale[pos],imputation_median=model.prep.median[pos],distance_weight=model.geometry.weights[gp] if is_genetic else 0.,full_model_coefficient=model.predictors['ABM_full'].coefficients[pos]))
    write_aggregate(out/'aggregate',dict(training_summary=pd.DataFrame([model.training_summary]),roles=pd.DataFrame(model.role_summary),tuning=pd.DataFrame(model.tuning),OOF_audit=pd.DataFrame(model.oof_audit),profile_features=pd.DataFrame(feature_rows)),dict(trait=model.config.trait,module='Twin training',run_id=digest))


def fit_cmd(args):
    d=npz(args.data);cfg=config_file(args.config)
    needed={'X','y','ids','groups','roles','features','genetic_indices','covariate_indices','csx_indices'}
    if not needed<=set(d) or 'test' in set(d['roles']):raise ValueError('Fit accepts development-only arrays')
    out=Path(args.output)
    if out.exists():raise FileExistsError('Choose a fresh fit directory; old models are never overwritten')
    model=TwinModel(cfg).fit(**{k:d[k] for k in needed})
    save_fit(model,out);print(json.dumps(dict(model=str(out/'_LOCAL_ONLY/model.joblib'),hash=state_hash(model))))


def project_cmd(args):
    model=joblib.load(args.model);d=npz(args.data)
    if set(d)!= {'X','ids','groups','features'}:raise ValueError('Project accepts only predictors/identity, never query outcomes')
    p=model.project(**d);out=Path(args.output)
    if out.exists():raise FileExistsError(out)
    out.parent.mkdir(parents=True,exist_ok=True);joblib.dump(p,out,compress=3)
    export_projection(p,out.with_suffix('.confidence.tsv.gz'));print(json.dumps(dict(query_N=len(p['eid']),selected_N=int(p['selected'].sum()))))


def benchmark_manifest(path,ids,truth,trait,endpoint):
    if not path:return {},[]
    path=Path(path);spec=json.loads(path.read_text());items=spec['methods'];external={};audit=[]
    for item in items:
        name=item['name'];source=path.parent/item['path']
        if item.get('trait')!=trait or item.get('endpoint_definition')!=endpoint or item.get('value_type')!='outcome_scale_prediction':
            raise ValueError('Comparator endpoint/units manifest mismatch: '+name)
        if item.get('test_roster_sha256')!=id_hash(ids):raise ValueError('Comparator does not declare this exact frozen test roster: '+name)
        if not item.get('predictions_frozen_without_test_labels',False):raise ValueError('Comparator must declare out-of-training prediction provenance')
        if not source.exists():audit.append(dict(method=name,status='file_unavailable'));continue
        if 'sha256' not in item or sha256(source)!=item['sha256']:raise ValueError('Comparator file checksum missing/mismatch: '+name)
        d=identified(read_table(source))
        if set(d.eid)!=set(ids):audit.append(dict(method=name,status='different_or_incomplete_test_roster'));continue
        d=d.set_index('eid').loc[ids];p=pd.to_numeric(d[item['prediction_column']],errors='raise').to_numpy(dtype=float,na_value=np.nan)
        if name in external:raise ValueError('Duplicate benchmark method')
        if item.get('outcome_column'):
            yy=pd.to_numeric(d[item['outcome_column']],errors='raise').to_numpy(dtype=float,na_value=np.nan)
            if not np.allclose(yy,truth.set_index('eid').loc[ids,'y'],equal_nan=True):raise ValueError('Comparator outcome definitions/units differ')
        external[name]=p;audit.append(dict(method=name,status='provider_declared_frozen_prediction_verified_file_and_roster',source_sha256=item['sha256'],training_budget=item.get('training_budget','not_recorded')))
    return external,audit


def evaluate_cmd(args):
    model=joblib.load(args.model);p=joblib.load(args.predictions);truth=identified(read_table(args.truth))
    external,availability=benchmark_manifest(args.comparators,p['eid'],truth,model.config.trait,args.endpoint)
    tables,status=evaluate(model,p,truth,external=external,bootstrap=args.bootstrap)
    if availability:tables['comparator_provenance']=pd.DataFrame(availability)
    out=Path(args.output)
    if out.exists():raise FileExistsError('Choose a new report directory')
    write_aggregate(out,tables,dict(trait=model.config.trait,module='Twin evaluation',run_id=p['frozen_hash']))
    json_write(out/'status.json',status)
    from twin_figures import draw_report
    draw_report(tables,out/'figures');print(json.dumps(status))


def run(args):
    traits=args.traits.split(',')
    if not traits or len(set(traits))!=len(traits) or not set(traits)<= {'height','ldl','t2dm'}:raise ValueError('Use unique height,ldl,t2dm traits')
    for trait in traits:
        cfg=config_file(args.config,trait)
        d,groups,meta=load_source(args.prepared or str(Path(args.cache_dir)/'prepared.joblib'),trait,
                                   args.data_file.replace('{trait}',trait) if args.data_file else None,
                                   args.schema.replace('{trait}',trait) if args.schema else None)
        raw,query,truth,flow,cfg,extra=numeric_cohort(d,groups,cfg,args.roles_file)
        endpoint=meta.get('endpoint_definition') or meta.get('prepared_settings',{}).get(trait+'_col',trait)
        if trait=='t2dm' and endpoint=='auto':endpoint='baseline_t2dm_from_Yr2e_Yt2e'
        endpoint=args.endpoint or endpoint
        if args.check:
            print(json.dumps(dict(trait=trait,input='validated',endpoint=endpoint,**flow)));continue
        root=Path(args.out_root)/'twin'/trait;root.mkdir(parents=True,exist_ok=True)
        runid=args.run_id or datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
        import re
        if not re.fullmatch(r'[A-Za-z0-9_-]+',runid):raise ValueError('Unsafe run name')
        dest=root/runid
        if dest.exists():raise FileExistsError('Existing run is immutable: '+str(dest))
        dest.mkdir();local=dest/'_LOCAL_ONLY';local.mkdir(mode=0o700)
        json_write(dest/'protocol.json',dict(config=asdict(cfg),input=meta,flow=flow,endpoint_definition=endpoint,test_roster_changed=False,validation_scope='exploratory_previously_inspected_test'))
        np.savez_compressed(local/'development.npz',**raw)
        np.savez_compressed(local/'test_X.npz',**query)
        truth.to_csv(local/'test_truth.tsv.gz',sep='\t',index=False,float_format='%.17g')
        pd.DataFrame(dict(eid=raw['ids'],group=raw['groups'],role=raw['roles'])).to_csv(local/'development_roles.tsv.gz',sep='\t',index=False)
        try:
            print(f'[GRID twin] FIT {trait}: {len(raw["ids"])} development; {len(query["ids"])} frozen test rows',flush=True)
            model=TwinModel(cfg).fit(**raw);save_fit(model,dest/'fit')
            audit=raw['roles']=='audit';pa=model.project(raw['X'][audit],raw['ids'][audit],raw['groups'][audit],raw['features'],audit=True)
            ta=pd.DataFrame(dict(eid=raw['ids'][audit],y=raw['y'][audit]))
            at,ast=evaluate(model,pa,ta,bootstrap=args.bootstrap)
            write_aggregate(dest/'audit',at,dict(trait=trait,module='Twin audit',run_id=pa['frozen_hash']))
            p=model.project(**query);joblib.dump(p,local/'predictions.joblib',compress=3)
            export_projection(p,local/'confidence.tsv.gz')
            # The predictor/membership was frozen and saved BEFORE test truth is used below.
            ext,availability=benchmark_manifest(args.comparators,p['eid'],truth,trait,endpoint)
            disco,disco_model,da=calibrate_disco(extra,model)
            availability.append(da)
            if disco is not None:
                ext['Disco_reference']=disco;joblib.dump(disco_model,local/'disco_calibration.joblib',compress=3)
            if args.analysis_root:
                cov=meta.get('prepared_settings',{}).get('covariates','')
                cov=cov.split(',') if isinstance(cov,str) and cov else list(groups.get('covariates',[]))
                try:
                    pp,ppstat=verified_prsformer(args.analysis_root,trait,d,truth,cfg,cov,endpoint)
                    availability.append(ppstat)
                    if pp is not None:ext['PRSformer']=pp
                except (ValueError,FileNotFoundError,RuntimeError,ImportError) as exc:
                    availability.append(dict(method='PRSformer',status='comparison_not_verified',reason=str(exc)))
            tables,status=evaluate(model,p,truth,external=ext,bootstrap=args.bootstrap)
            if ext:joblib.dump(ext,local/'comparator_predictions.joblib',compress=3)
            tables['comparator_provenance']=pd.DataFrame(availability)
            write_aggregate(dest/'test',tables,dict(trait=trait,module='Twin test',run_id=p['frozen_hash']))
            json_write(dest/'test/status.json',status)
            from twin_figures import draw_report
            draw_report(tables,dest/'test/figures')
            json_write(dest/'COMPLETE.json',dict(status='completed',model_hash=p['frozen_hash'],test_roster_hash=id_hash(query['ids']),clinical_release=False))
            json_write(root/'CURRENT.json',dict(run=runid,model_hash=p['frozen_hash']))
            print(f'[GRID twin] COMPLETE {trait}: {dest}',flush=True)
        except BaseException:
            json_write(dest/'FAILED.json',dict(status='failed_partial_run_not_promoted',previous_CURRENT_unchanged=True))
            raise


def demo(args):
    root=Path(args.output)
    if root.exists():raise FileExistsError(root)
    root.mkdir(parents=True);rng=np.random.default_rng(219)
    n=2400;ids=np.asarray([f'SYNTHETIC_{i}' for i in range(n)]);pop=np.where(np.arange(n)%4,'EUR','OTHER')
    g=rng.normal(size=(n,6));age=rng.normal(size=n);sex=rng.integers(0,2,n)
    # Heteroscedasticity exercises confidence evaluation; no tests require an ABM win.
    y=170+4*g[:,0]+.5*g[:,1]+2*sex+rng.normal(size=n)*(3+np.abs(g[:,2]))
    if args.binary:y=rng.binomial(1,1/(1+np.exp(-(-1+.6*g[:,0]+.2*age))))
    d=pd.DataFrame(g,columns=[f'csx_{i}' for i in range(6)]);d['age']=age;d['sex']=sex
    d['eid']=ids;d['family_id']=np.asarray([f'F_{i//2}' for i in range(n)]);d['ancestry']=pop;d['y']=y
    d['split']=np.where(np.arange(n)<n//2,'train','test');d.loc[n//2:n//2+4,[f'csx_{i}' for i in range(6)]]=np.nan
    data=root/'synthetic.tsv.gz';d.to_csv(data,sep='\t',index=False)
    schema=root/'schema.json';json_write(schema,dict(feature_groups=dict(csx=[f'csx_{i}' for i in range(6)],covariates=['age','sex']),metadata=dict(synthetic_only=True)))
    cfg=replace(Config(),folds=3,repeats=1,min_role_n=12,min_class_n=3,k=12,min_twins=3,min_group_ess=2,gate_iterations=8,bootstrap=25,ridge_grid=(10.,100.),logistic_c_grid=(.1,))
    config=root/'config.json';json_write(config,asdict(cfg))
    trait='t2dm' if args.binary else 'height'
    run(argparse.Namespace(traits=trait,prepared=None,cache_dir='unused',data_file=str(data),schema=str(schema),roles_file=None,config=str(config),check=False,out_root=str(root/'abm'),analysis_root=None,run_id='synthetic',endpoint='SIMULATED_ONLY',comparators=None,bootstrap=25))
    (root/'SYNTHETIC_ONLY.txt').write_text('Synthetic execution check, not UK Biobank evidence or performance validation.\n')


def main():
    a=argparse.ArgumentParser(description=__doc__);sp=a.add_subparsers(dest='command',required=True)
    r=sp.add_parser('run');r.add_argument('--traits',default='height,ldl,t2dm');r.add_argument('--cache-dir',default='/tmp/grid-cache/grid');r.add_argument('--prepared');r.add_argument('--data-file');r.add_argument('--schema');r.add_argument('--roles-file')
    r.add_argument('--config',default=str(DEFAULT_CONFIG),help='Configuration JSON (default: config/twin.json beside grid.sh)');r.add_argument('--out-root',default='/mnt/d/analysis/grid/abm');r.add_argument('--analysis-root',default='/mnt/d/analysis/grid');r.add_argument('--run-id');r.add_argument('--check',action='store_true');r.add_argument('--endpoint');r.add_argument('--comparators');r.add_argument('--bootstrap',type=int);r.add_argument('--threads',type=int,default=4,help='Numeric worker threads (default: 4)')
    f=sp.add_parser('fit');f.add_argument('--data',required=True);f.add_argument('--config',required=True);f.add_argument('--output',required=True);f.add_argument('--threads',type=int,default=4)
    p=sp.add_parser('project');p.add_argument('--model',required=True);p.add_argument('--data',required=True);p.add_argument('--output',required=True);p.add_argument('--threads',type=int,default=4)
    e=sp.add_parser('evaluate');e.add_argument('--model',required=True);e.add_argument('--predictions',required=True);e.add_argument('--truth',required=True);e.add_argument('--output',required=True);e.add_argument('--comparators');e.add_argument('--endpoint',required=True);e.add_argument('--bootstrap',type=int);e.add_argument('--threads',type=int,default=4)
    d=sp.add_parser('demo');d.add_argument('--output',required=True);d.add_argument('--binary',action='store_true');d.add_argument('--threads',type=int,default=2)
    args=a.parse_args()
    if args.threads<1:a.error('--threads must be positive')
    action=dict(run=run,fit=fit_cmd,project=project_cmd,evaluate=evaluate_cmd,demo=demo)[args.command]
    os.umask(0o077)
    with threadpool_limits(args.threads):action(args)

if __name__=='__main__':
    try:main()
    except (ValueError,FileNotFoundError,FileExistsError,RuntimeError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);sys.exit(2)
