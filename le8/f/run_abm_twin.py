#!/usr/bin/env python3
"""Separate fit -> project -> evaluate commands. Only fit sees development labels."""
from __future__ import annotations
import argparse
from dataclasses import asdict
from pathlib import Path
import json
import joblib
import numpy as np
from threadpoolctl import threadpool_limits
from abm_twin_core import Config,TwinABM,ROLES,VERSION,json_write,tie_hash,model_hash
from abm_twin_evaluate import evaluate


def read_npz(path):
    with np.load(path,allow_pickle=False) as f: return {k:f[k] for k in f.files}


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    s=p.add_subparsers(dest='command',required=True)
    a=s.add_parser('fit');a.add_argument('--development',required=True);a.add_argument('--output',required=True);a.add_argument('--config')
    a=s.add_parser('project');a.add_argument('--model',required=True);a.add_argument('--features',required=True);a.add_argument('--output',required=True)
    a=s.add_parser('evaluate');a.add_argument('--model',required=True);a.add_argument('--predictions',required=True);a.add_argument('--truth',required=True);a.add_argument('--output',required=True);a.add_argument('--bootstrap',type=int)
    a=s.add_parser('demo');a.add_argument('--output',required=True);a.add_argument('--seed',type=int,default=12345)
    return p


def load_config(path):
    values=json.loads(Path(path).read_text()) if path else {}
    for key in ('coverages','c_grid','always_include_indices'):
        if key in values: values[key]=tuple(values[key])
    return Config(**values)


def write_result(out,result):
    out.mkdir(parents=True,exist_ok=True)
    for name,value in result.items(): json_write(out/(name+'.json'),value)
    import pandas as pd
    from research_io import write_tables
    from research_figures import plot_abm
    tables={k:pd.DataFrame(v) for k,v in result.items() if isinstance(v,list) and (not v or isinstance(v[0],dict))}
    tables['status']=pd.DataFrame([result['status']])
    for d in tables.values():
        d['evaluation_scope']=result.get('evaluation_scope','unspecified')
        d['grouping_unit']=result.get('grouping_unit','family')
    write_tables(out,tables,'ABM twin',result.get('layer','prot'),result.get('trait','cvd_cad'),run_id=result.get('run_id',''),workbook='ABM_twin.results.xlsx')
    plot_abm(tables,out/'figures')


def fit_command(development,output,config):
    out=Path(output)
    if (out/'model_LOCAL_ONLY.joblib').exists(): raise FileExistsError('Use a new output directory; old fits are never overwritten')
    data=read_npz(development)
    required={'X','time','event','eid','group','role','feature_names'}
    if not required<=set(data): raise ValueError('Missing development fields: '+str(required-set(data)))
    model=TwinABM(config).fit(data['X'],data['time'],data['event'],data['eid'],data['group'],data['role'],data['feature_names'])
    out.mkdir(parents=True,exist_ok=True)
    joblib.dump(model,out/'model_LOCAL_ONLY.joblib',compress=3)
    json_write(out/'MODEL_FROZEN.json',dict(version=VERSION,hash=model_hash(model),config=asdict(config),test_outcomes_used=False))
    np.savez_compressed(out/'donor_quality_LOCAL_ONLY.npz',**model.donor_table)
    json_write(out/'training_audit.json',model.training_audit)
    json_write(out/'model_tuning.json',model.model_tuning)
    json_write(out/'selected_panel.json',dict(features=model.panel_names, input_indices=model.panel_indices))
    import pandas as pd
    from research_io import write_tables
    audit=model.training_audit
    summary={k:v for k,v in audit.items() if not isinstance(v,(list,dict))}
    write_tables(out/'aggregate',{
        'training_summary':pd.DataFrame([summary]),
        'training_folds':pd.DataFrame(audit['oof']),
        'pruning_audit':pd.DataFrame(audit['pruning']),
        'selected_panel':pd.DataFrame({'feature':model.panel_names,'input_index':model.panel_indices,
           'fixed_covariate':[j in config.always_include_indices for j in model.panel_indices]}),
        'model_tuning':pd.DataFrame(model.model_tuning),
        'role_counts':pd.DataFrame([dict(role=r,N=v['N'],cases=v['cases']) for r,v in audit['roles'].items()])
    },'ABM twin',config.layer,config.trait,model_hash(model))
    ix=data['role']=='calibration_audit'
    before=model_hash(model)
    pred=model.predict(data['X'][ix],data['eid'][ix],data['group'][ix],data['feature_names'],allow_development_audit=True)
    result=evaluate(model,pred,data['time'][ix],data['event'][ix],data['eid'][ix],full_curves=False)
    assert before==model_hash(model),'Audit changed the frozen model'
    write_result(out/'audit',result)
    return model


def demo_command(out,seed):
    """Synthetic mechanics demonstration only; no clinical performance claim."""
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    n=2400;p=12
    X=rng.normal(size=(n,p))
    # Ordinary stochastic survival data; deliberately not engineered to guarantee a pruning gain.
    eta=-2.3+1.1*X[:,0]-.7*X[:,1]+.4*X[:,2]
    rate=np.exp(eta)/10
    failure=rng.exponential(1/rate)
    censor=rng.uniform(11,18,n)
    time=np.minimum(failure,censor);event=(failure<=censor).astype(int)
    ids=np.array([f'synthetic_{i}' for i in range(n)])
    groups=ids.copy();names=np.array([f'X{j}' for j in range(p)])
    # Split is outcome-blind. This demo is not a UKB result and is not tuned for gains.
    order=np.argsort([tie_hash(v,seed) for v in ids]);dev=order[:n//2];test=order[n//2:]
    counts=np.array([.6,.1,.1,.1,.1]);cut=np.cumsum(np.r_[0,np.round(counts*len(dev)).astype(int)]);cut[-1]=len(dev)
    roles=np.empty(len(dev),dtype='U25')
    for j,r in enumerate(ROLES):roles[cut[j]:cut[j+1]]=r
    np.savez_compressed(out/'development.npz',X=X[dev],time=time[dev],event=event[dev],eid=ids[dev],group=groups[dev],role=roles,feature_names=names)
    np.savez_compressed(out/'test_X.npz',X=X[test],eid=ids[test],group=groups[test],feature_names=names)
    np.savez_compressed(out/'test_truth.npz',time=time[test],event=event[test],eid=ids[test])
    cfg=Config(seed=seed,folds=3,repeats=1,c_grid=(.1,),min_events=3,k=32,min_twins=8,min_family_ess=6,
               dimensions=8,gate_iterations=20,bootstrap=25,audit_min_cases=3,audit_min_controls=3,max_iter=3000)
    model=fit_command(out/'development.npz',out/'fit',cfg)
    pred=model.predict(X[test],ids[test],groups[test],names)
    joblib.dump(pred,out/'test_predictions_LOCAL_ONLY.joblib')
    result=evaluate(model,pred,time[test],event[test],ids[test],bootstrap=25)
    write_result(out/'test_report',result)
    print(json.dumps({'synthetic_only':True,'selected_n':int(pred['selected'].sum()),'test_n':len(test),'status':result['status']},indent=2))


def main():
    a=parser().parse_args()
    if a.command=='fit':
        config=load_config(a.config)
        with threadpool_limits(config.cores): fit_command(a.development,a.output,config)
    elif a.command=='demo':
        with threadpool_limits(2): demo_command(a.output,a.seed)
    else:
        model=joblib.load(a.model)  # Only load a trusted model created by this package.
        with threadpool_limits(model.config.cores):
            if a.command=='project':
                data=read_npz(a.features)
                permitted={'X','eid','group','feature_names'}
                if set(data)!=permitted: raise ValueError('Projection file must contain exactly X, eid, group, feature_names; no outcomes')
                pred=model.predict(data['X'],data['eid'],data['group'],data['feature_names'])
                output=Path(a.output)
                if output.exists(): raise FileExistsError(output)
                output.parent.mkdir(parents=True,exist_ok=True)
                joblib.dump(pred,output,compress=3)
                import pandas as pd
                table=pd.DataFrame({'eid':pred['eid'],'selected':pred['selected'],'reason':pred['reason'],
                    'candidate_risk':pred['probability'][model.config.primary],'fallback_risk':pred['probability']['elasticnet_full'],
                    'policy_risk':pred['probability']['fallback_policy'],**pred['diagnostics']})
                table.to_csv(output.with_name(output.stem+'.confidence.tsv.gz'),sep='\t',index=False)
                json_write(output.with_suffix('.manifest.json'),dict(model_hash=model_hash(model),version=VERSION,outcomes_used=False,selected_N=int(pred['selected'].sum())))
            else:
                pred=joblib.load(a.predictions)
                truth=read_npz(a.truth)
                before=model_hash(model)
                result=evaluate(model,pred,truth['time'],truth['event'],truth['eid'],bootstrap=a.bootstrap)
                assert before==model_hash(model)
                write_result(Path(a.output),result)

if __name__=='__main__':main()
