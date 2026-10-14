#!/usr/bin/env python3
"""Prepare -> fit -> project -> evaluate. All individual artifacts stay local."""
from __future__ import annotations
import argparse,json,hashlib,os
from dataclasses import asdict,replace
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from molecular_deviation import DeviationConfig,DeviationStudy,objective_endpoint,evaluate_deviation
from research_io import id_table,read_table,load_pgs,merge_native_scores,sha256,write_tables
from abm_twin_core import model_hash,json_write


def read_config(path):
    values=json.loads(Path(path).read_text()) if path else {}
    if 'c_grid' in values:values['c_grid']=tuple(values['c_grid'])
    return DeviationConfig(**values)


def fixed_roles(groups,seed):
    # No labels, disease prevalence or future outcomes enter the allocation.
    def role(g):
        v=int(hashlib.sha256((str(seed)+'|'+g).encode()).hexdigest()[:16],16)/2**64
        return 'train' if v<.60 else 'tune' if v<.75 else 'calibrate' if v<.90 else 'test'
    return np.array([role(str(g)) for g in groups])


def prepare(args):
    cfg=read_config(args.config);cfg.validate();out=Path(args.output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Choose a new prepared directory')
    cohort=id_table(args.cohort_file);omics=id_table(args.omics_file)
    if 'known_cad_at_index' not in cohort:raise ValueError('Supply verified baseline CAD status (unknown stays missing)')
    ids=cohort.index[cohort.index.isin(omics.index)]
    if not len(ids):raise ValueError('No shared phenotype/omics IDs')
    d=cohort.loc[ids].copy();m=omics.loc[ids].drop(columns=['eid'])
    if args.feature_file:
        features=[v.strip() for v in Path(args.feature_file).read_text().splitlines() if v.strip()]
        if len(set(features))!=len(features) or not set(features)<=set(m):raise ValueError('Invalid feature allowlist')
        m=m[features]
    columns=[v.strip() for v in args.clinical.split(',') if v.strip()]
    if not columns or not set(columns)<=set(d):raise ValueError('Supply existing numeric clinical/covariate column names')
    forbidden={'objective_cad','known_cad_at_index','event','time','y','endpoint_eligible','intervening_cad'}
    if set(columns)&forbidden:raise ValueError('Outcome or disease ascertainment column cannot be a predictor')
    C=d[columns].apply(pd.to_numeric,errors='raise').to_numpy(float)
    names,G,audit=load_pgs(args.pgs_root,cfg.layer,ids,m.columns,args.pgs_mapping)
    M=m[names].apply(pd.to_numeric,errors='raise').to_numpy(float)
    if args.group_col:
        if args.group_col not in d or d[args.group_col].isna().any():raise ValueError('Missing family/group identifiers')
        group=d[args.group_col].astype(str).to_numpy()
    else:group=np.asarray(ids,str)
    support=dict(family_status='verified_supplied_groups' if args.group_col else 'unavailable_individual_ID_groups_only',
                 independent_family_support_verified=bool(args.group_col),
                 grouping_unit='family' if args.group_col else 'individual_unverified',
                 assume_unrelated=bool(args.assume_unrelated),validation_status='exploratory')
    cfg=replace(cfg,grouping_unit=support['grouping_unit'])
    if args.role_col:
        if args.role_col not in d:raise ValueError('Missing role column')
        role=d[args.role_col].astype(str).to_numpy()
    else:role=fixed_roles(group,cfg.seed)
    if set(role)!={'train','tune','calibrate','test'}:raise ValueError('All four roles required')
    rg=pd.DataFrame({'g':group,'role':role})
    if rg.groupby('g').role.nunique().max()>1:raise ValueError('Family crosses roles')
    known=pd.to_numeric(d.known_cad_at_index,errors='coerce').to_numpy()
    reference=known==0  # No use of future disease status, nor an assumption of plaque-free health.
    if args.reference_col:
        if args.reference_col not in d:raise ValueError('Missing predeclared reference eligibility')
        reference &= pd.to_numeric(d[args.reference_col],errors='coerce').to_numpy()==1
    eligible,y,reason=objective_endpoint(d,cfg)
    out.mkdir(parents=True);local=out/'_LOCAL_ONLY';local.mkdir(mode=0o700)
    group=np.asarray(group,str);role=np.asarray(role,str)
    arrays=dict(M=M,G=G,C=C,eid=np.asarray(ids,str),group=group,role=role,reference=reference,
                feature_names=np.asarray(names),covariate_names=np.asarray(columns))
    dev=role!='test'
    np.savez_compressed(local/'development.npz',**{k:v[dev] if k not in ('feature_names','covariate_names') else v for k,v in arrays.items()},
                        y=y[dev],endpoint_eligible=eligible[dev])
    np.savez_compressed(local/'test_X.npz',M=M[~dev],G=G[~dev],C=C[~dev],eid=np.asarray(ids,str)[~dev],group=group[~dev],
                        feature_names=np.asarray(names),covariate_names=np.asarray(columns))
    # Full ascertainment remains a private audit table. Labels are in a separate file.
    pd.DataFrame({'eid':np.asarray(ids,str)[~dev],'y':y[~dev],'endpoint_eligible':eligible[~dev],
                  'exclusion_reason':reason[~dev]}).to_csv(local/'test_truth.tsv.gz',sep='\t',index=False)
    pd.DataFrame({'eid':ids,'group':group,'role':role,'reference_eligible':reference,
                  'endpoint_eligible':eligible,'endpoint_reason':reason}).to_csv(local/'roster.tsv.gz',sep='\t',index=False)
    counts=pd.DataFrame({'role':role,'endpoint_reason':reason}).groupby(['role','endpoint_reason']).size().reset_index(name='N')
    status=pd.DataFrame([dict(task=cfg.task,status='objective_endpoint_available' if eligible.any() else 'objective_endpoint_unavailable',
                      molecular_N=len(M),assays=len(names),endpoint_N=int(eligible.sum()),clinical_release=False,**support)])
    write_tables(out/'aggregate',{'pgs_quality_audit':audit,'endpoint_flow':counts,'status':status},'deviation',cfg.layer,cfg.trait,'prepared')
    json_write(out/'config.json',asdict(cfg))
    json_write(out/'preparation.json',dict(config=asdict(cfg),
               group_assumption=args.group_col or ('declared_unrelated' if args.assume_unrelated else 'individual_ID_only_unverified'),
               **support,
               reference_definition='known_cad_at_index==0'+(' and '+args.reference_col+'==1' if args.reference_col else ''),
               evaluation_status='exploratory_until_independent_validation',inputs={str(p):sha256(p) for p in (args.cohort_file,args.omics_file)},
               warning='PGS GWAS overlap is not corrected by calibration cross-fitting; objective ascertainment may be selected.'))
    print(json.dumps(dict(prepared=str(out),N=len(M),features=len(names),objective_N=int(eligible.sum()))))


def read_npz(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}


def fit(args):
    cfg=read_config(args.config);d=read_npz(args.data);out=Path(args.output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Choose a new model output directory')
    roles=np.asarray(d['role'],str);tr=roles=='train'
    if set(roles)!={'train','tune','calibrate'}:raise ValueError('No test labels allowed in fit')
    model=DeviationStudy(cfg).fit_reference(d['M'][tr],d['G'][tr],d['C'][tr],d['eid'][tr],d['group'][tr],
                d['reference'][tr],d['feature_names'],d['covariate_names'])
    model.fit_disease_models(d);out.mkdir(parents=True)
    local=out/'_LOCAL_ONLY';local.mkdir(mode=0o700);joblib.dump(model,local/'model.joblib',compress=3)
    frozen=model_hash(model);json_write(out/'MODEL_FROZEN.json',{'hash':frozen,'version':'2026-10-14','config':asdict(cfg),
                        'test_outcomes_used':False,'status':model.classification_status,'clinical_release':False})
    write_tables(out/'aggregate',{'reference_calibration':pd.DataFrame(model.reference_table),
                 'reference_folds':pd.DataFrame(model.reference_audit),
                 'model_tuning':pd.DataFrame(getattr(model,'model_tuning',[])),
                 'status':pd.DataFrame([{'status':model.classification_status,'task':cfg.task,'clinical_release':False}])},
                 'deviation',cfg.layer,cfg.trait,frozen[:16])
    print(json.dumps({'model':str(local/'model.joblib'),'status':model.classification_status}))


def project(args):
    # Only trusted local model artifacts; never load model files from an untrusted source.
    model=joblib.load(args.model);d=read_npz(args.data)
    allowed={'M','G','C','eid','group','feature_names','covariate_names'}
    if set(d)!=allowed:raise ValueError('Projection accepts X and identity only, no outcome fields')
    p=model.project(d['M'],d['G'],d['C'],d['eid'],d['group'],d['feature_names'],d['covariate_names'])
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists():raise FileExistsError(out)
    joblib.dump(p,out,compress=3)
    # Arrays already provide per-person, per-feature expected and deviation values locally.
    print(json.dumps({'prediction':str(out),'N':len(p['eid']),'hash':p['frozen_hash']}))


def evaluate(args):
    model=joblib.load(args.model);p=joblib.load(args.predictions);truth=id_table(args.truth).reset_index(drop=True)
    if truth.endpoint_eligible.dtype==object:
        val=truth.endpoint_eligible.astype(str).str.lower()
        if not val.isin(['true','false','1','0']).all():raise ValueError('Invalid boolean eligibility')
        truth['endpoint_eligible']=val.isin(['true','1'])
    if set(truth.eid)!=set(p['eid']):raise ValueError('Different prediction/truth populations')
    truth=truth.set_index('eid').loc[p['eid']].reset_index()
    tables,status=evaluate_deviation(model,p,truth,args.output)
    from research_figures import plot_deviation
    plot_deviation(tables,Path(args.output)/'figures',args.anchors.split(',') if args.anchors else [])
    print(json.dumps({'output':args.output,'classification_status':model.classification_status}))


def demo(args):
    # Synthetic smoke test, not evidence of UKB validity or discovery.
    cfg=DeviationConfig(reference_min_n=20,min_class_n=5,bootstrap=25,c_grid=(.1,),folds=3)
    rng=np.random.default_rng(73);n=1600;p=12
    ids=np.asarray(['DEMO_'+str(i) for i in range(n)]);g=rng.normal(size=(n,p));C=rng.normal(size=(n,3));state=rng.binomial(1,.25,n)
    m=.5*g+.2*C[:,[0]]+rng.normal(size=(n,p));m[:,0]+=state*.7
    roles=np.repeat(['train','tune','calibrate','test'],[800,300,300,200]);ref=np.ones(n,bool)
    root=Path(args.output);root.mkdir(parents=True,exist_ok=False);local=root/'_LOCAL_ONLY';local.mkdir(mode=0o700)
    common={'M':m,'G':g,'C':C,'eid':ids,'group':ids,'reference':ref,'role':roles}
    dev=roles!='test';names=np.asarray(['MET_'+str(i) for i in range(p)]);covs=np.asarray(['age_z','sex_code','PC1'])
    np.savez_compressed(local/'dev.npz',**{k:v[dev] for k,v in common.items()},y=state[dev],endpoint_eligible=np.ones(dev.sum(),bool),feature_names=names,covariate_names=covs)
    np.savez_compressed(local/'X.npz',M=m[~dev],G=g[~dev],C=C[~dev],eid=ids[~dev],group=ids[~dev],feature_names=names,covariate_names=covs)
    pd.DataFrame({'eid':ids[~dev],'y':state[~dev],'endpoint_eligible':True}).to_csv(local/'truth.tsv',sep='\t',index=False)
    json_write(root/'config.json',asdict(cfg));A=argparse.Namespace
    fit(A(data=str(local/'dev.npz'),config=str(root/'config.json'),output=str(root/'model')))
    model=root/'model/_LOCAL_ONLY/model.joblib';project(A(model=str(model),data=str(local/'X.npz'),output=str(local/'prediction.joblib')))
    evaluate(A(model=str(model),predictions=str(local/'prediction.joblib'),truth=str(local/'truth.tsv'),output=str(root/'results'),anchors='MET_0,MET_1'))
    (root/'SYNTHETIC_ONLY.txt').write_text('Synthetic pipeline smoke test. Not a UKB result. Endpoint simulated, no actual imaging.\n')


def main():
    a=argparse.ArgumentParser(description=__doc__);s=a.add_subparsers(dest='command',required=True)
    p=s.add_parser('prepare');p.add_argument('--cohort-file',required=True);p.add_argument('--omics-file',required=True)
    p.add_argument('--pgs-root',required=True);p.add_argument('--clinical',required=True);p.add_argument('--config',required=True)
    p.add_argument('--pgs-mapping');p.add_argument('--feature-file')
    p.add_argument('--group-col',help='Optional family/group column; absent: unverified individual ID groups')
    p.add_argument('--role-col');p.add_argument('--reference-col')
    p.add_argument('--assume-unrelated',action='store_true',help='Optional user declaration, not proof of independent families; unnecessary for running')
    p.add_argument('--output',required=True)
    p=s.add_parser('fit');p.add_argument('--data',required=True);p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    p=s.add_parser('project');p.add_argument('--model',required=True);p.add_argument('--data',required=True);p.add_argument('--output',required=True)
    p=s.add_parser('evaluate');p.add_argument('--model',required=True);p.add_argument('--predictions',required=True);p.add_argument('--truth',required=True);p.add_argument('--output',required=True);p.add_argument('--anchors',default='L_VLDL_TG.pct,Total_TG,ApoB')
    p=s.add_parser('merge-pgs');p.add_argument('--manifest',required=True);p.add_argument('--output',required=True)
    p=s.add_parser('demo');p.add_argument('--output',required=True)
    args=a.parse_args()
    if args.command=='merge-pgs':print(json.dumps(merge_native_scores(args.manifest,args.output)));return
    globals()[args.command](args)

if __name__=='__main__':
    os.umask(0o077)
    main()
