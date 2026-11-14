"""Endpoint-wise disk-backed evaluation and matched-prefix tabular baselines."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common import digest, file_digest, log, outdir
from data import CohortDataset, load_meta
from metrics import CensoringKM, calibration_bins, evaluate_endpoint, macro_average
from train import extract_states, risk_cache, state_cache_key


def censor_for(arrays, endpoint, meta):
    death = meta['endpoints'][endpoint] == ['DEATH']
    return np.asarray(arrays['exit_time']) if death else np.minimum(arrays['exit_time'],arrays['death_time'])


def score_matrix(cfg, model_name, predictions, *, probability_horizon=None, split='test',
                 horizons=None, suffix=None):
    """Never builds a people x endpoints long table or fits censoring on test."""
    training, tm = risk_cache(cfg,'train'); testing, em = risk_cache(cfg,split)
    meta = load_meta(cfg); names=tm['endpoint_names']
    if names != em['endpoint_names'] or names != meta['endpoint_names']:
        raise ValueError('Endpoint order differs between training, evaluation cache and prepared metadata')
    if tm['key'] != em['key'] or tm['key'] != state_cache_key(cfg):
        raise ValueError('Risk cache is stale for the current world model/landmark settings. Run readout first.')
    if predictions.shape != testing['event_time'].shape: raise ValueError('Prediction matrix and risk set do not align')
    rows=[]; bins=[]
    horizons = horizons or cfg['evaluation']['horizons']
    for j,name in enumerate(names):
        tr_event=np.asarray(training['event_time'][:,j]); tr_censor=censor_for(training,name,meta)
        tr_mask=np.asarray(training['risk_mask'][:,j])&(tr_censor>0)
        ref_time=np.minimum(tr_event,tr_censor)[tr_mask]
        ref_event=(tr_event<=tr_censor)[tr_mask]
        te_event=np.asarray(testing['event_time'][:,j]); te_censor=censor_for(testing,name,meta)
        mask=np.asarray(testing['risk_mask'][:,j]); scores=np.asarray(predictions[:,j])
        cached_c_index=None
        for h in horizons:
            result=evaluate_endpoint(te_event,te_censor,scores,horizon_years=float(h),
                ref_time=ref_time,ref_event=ref_event,risk_mask=mask,
                probabilities=scores if probability_horizon is not None else None,
                probability_horizon_years=probability_horizon,
                weight_cap=cfg['evaluation']['weight_cap'],compute_c_index=cached_c_index is None)
            if cached_c_index is None: cached_c_index=result['harrell_c_index']
            else: result['harrell_c_index']=cached_c_index
            result.update(model=model_name,endpoint=name,split=split)
            rows.append(result)
            if probability_horizon is not None and np.isclose(float(h),probability_horizon):
                good=mask&np.isfinite(scores)&(te_event>0)&(te_censor>0)
                try: reference=CensoringKM.fit(ref_time,ref_event)
                except ValueError: reference=None
                for item in calibration_bins(scores[good],np.minimum(te_event,te_censor)[good],
                    (te_event<=te_censor)[good],horizon_years=float(h),reference=reference,
                    weight_cap=cfg['evaluation']['weight_cap']):
                    bins.append(dict(model=model_name,endpoint=name,split=split,horizon_years=h,**item))
        if (j+1)%100==0: log(f'Evaluated {model_name}: {j+1}/{len(names)} endpoints')
    frame=pd.DataFrame(rows); folder=outdir(cfg,'results'); filename=suffix or model_name
    frame.to_csv(folder/f'{filename}.metrics.tsv',sep='\t',index=False)
    pd.DataFrame(bins).to_csv(folder/f'{filename}.calibration.tsv',sep='\t',index=False)
    macro_average(frame).to_csv(folder/f'{filename}.macro.tsv',sep='\t',index=False)
    group='readout' if model_name=='flux' else 'rollout' if model_name=='rollout' else 'baseline'
    provenance=dict(risk_key=tm['key'],endpoint_names=names,model=model_name,config=cfg[group],
                    seed=cfg['seed'],evaluation=cfg['evaluation'])
    if model_name=='flux': provenance['readout_sha256']=file_digest(Path(cfg['paths']['out'])/'model/readout.pt')
    (folder/f'{filename}.provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    return frame


def evaluate(cfg):
    folder=outdir(cfg,'results')
    description=json.loads((folder/'flux.prediction.json').read_text())
    _,manifest=risk_cache(cfg,'test')
    if (description.get('states_key') != manifest['key'] or
        description['endpoint_names'] != manifest['endpoint_names'] or
        description.get('readout_sha256') != file_digest(Path(cfg['paths']['out'])/'model/readout.pt')):
        raise ValueError('Saved prediction provenance differs from the current risk/readout cache. Run readout again.')
    predictions=np.load(folder/'flux.test.npy',mmap_mode='r')
    result=score_matrix(cfg,'flux',predictions,probability_horizon=description['horizon_years'])
    log('Five-year recorded-label metrics and censoring-aware metrics saved separately in results/flux.metrics.tsv')
    return result


def _prefix_feature_indices(meta, maximum):
    # Selection uses only training observation counts, never outcome association.
    ranked=sorted(range(len(meta['features'])), key=lambda i:(-meta['numeric_scaler'][meta['features'][i]]['n'],meta['features'][i]))
    return sorted(ranked[:int(maximum)])


def prefix_features(cfg,split,*,replace=False):
    meta=load_meta(cfg); indices=_prefix_feature_indices(meta,cfg['baseline']['max_numeric_features'])
    dataset=CohortDataset(cfg,split,mode='risk',horizon_years=cfg['readout']['horizon_years'])
    risk,manifest=risk_cache(cfg,split)
    folder=outdir(cfg,f'data/risk_{split}')
    key=digest(dict(states_key=manifest['key'],indices=indices))
    schema=folder/'baseline.json'; output=folder/'baseline.npy'
    if output.exists() and schema.exists() and not replace:
        if json.loads(schema.read_text())['key']==key: return np.load(output,mmap_mode='r')
    vocab=meta['vocab']; c=len(meta['context_features']); v=len(vocab); f=len(indices)
    width=3+2*c+v+2*f
    matrix=np.lib.format.open_memmap(output,mode='w+',dtype='float32',shape=(len(dataset),width))
    for i in range(len(dataset)):
        item=dataset[i]
        if item['eid'] != str(risk['eid'][i]): raise ValueError('Baseline and world-model risk rows do not align')
        tokens=item['tokens']; ages=item['times']
        clinical=(tokens>=3)
        counts=np.bincount(tokens[clinical],minlength=v).astype('float32')
        last_gap=float(item['query_age']-ages[clinical][-1]) if clinical.any() else np.nan
        last_values=np.full(f,np.nan,dtype='float32'); observed=np.zeros(f,dtype='float32')
        for k,j in enumerate(indices):
            positions=np.flatnonzero(item['numeric_mask_flat'][:,j])
            if len(positions):
                last_values[k]=item['numeric_flat'][positions[-1],j]; observed[k]=1
        x=np.r_[item['query_age'],last_gap,np.log1p(clinical.sum()),item['context'],item['context_mask'],
                np.log1p(counts),last_values,observed]
        matrix[i]=x
        if (i+1)%10000==0: log(f'Baseline prefix {split}: {i+1:,}/{len(dataset):,}')
    matrix.flush()
    columns=['age','time_since_last_diagnosis','log1p_diagnosis_count']
    columns+=meta['context_features']+['available:'+x for x in meta['context_features']]
    columns+=['log1p_count:'+x for x in vocab]
    columns+=[meta['features'][j] for j in indices]+['available:'+meta['features'][j] for j in indices]
    schema.write_text(json.dumps(dict(key=key,columns=columns,numeric_selection='training observation count',
        note='Independent matched-prefix representation; not the paper authors 213-dimensional comparator features.'),indent=2)+'\n')
    return matrix


def run_baselines(cfg,device='auto',*,replace=False):
    from baseline import fit_predict_baselines
    # Uses the exact rows/landmarks already frozen for HealthFlux.
    extract_states(cfg,device,replace=False)
    split_arrays={s:risk_cache(cfg,s)[0] for s in ('train','validation','test')}
    matrices={s:prefix_features(cfg,s,replace=replace) for s in split_arrays}
    meta=load_meta(cfg); names=meta['endpoint_names']; h=float(cfg['readout']['horizon_years'])
    train=split_arrays['train']; valid=split_arrays['validation']
    bc=cfg['baseline']; modeldir=outdir(cfg,'model/baseline')
    # Byte labels are much smaller than float outcome matrices at UKB scale.
    labels={s:(a['event_time']<=h).astype('uint8') for s,a in split_arrays.items() if s!='test'}
    survival={}
    if 'cox' in bc['methods']:
        for s in ('train','validation'):
            a=split_arrays[s]
            time=np.empty_like(a['event_time']); status=np.zeros(a['event_time'].shape,dtype=bool)
            for j,name in enumerate(names):
                censor=censor_for(a,name,meta)
                time[:,j]=np.minimum(a['event_time'][:,j],censor)
                status[:,j]=a['event_time'][:,j]<=censor
            survival[s]=(time,status)
    log(f'Fitting endpoint-specific baselines: {", ".join(bc["methods"])}')
    result=fit_predict_baselines(matrices['train'],labels['train'],train['risk_mask'],
        matrices['validation'],labels['validation'],valid['risk_mask'],matrices['test'],
        methods=bc['methods'],seed=int(cfg['seed']),endpoint_names=names,
        survival_train=survival.get('train'),survival_val=survival.get('validation'),
        horizon_years=h,min_cases=int(bc['min_cases']),max_iter=int(bc['max_iter']),
        n_jobs=int(bc['n_jobs']),output_dir=modeldir,options=bc.get('options'))
    pd.DataFrame(result['status']).to_csv(outdir(cfg,'results')/'baseline.status.tsv',sep='\t',index=False)
    for method,values in result['predictions'].items():
        score_matrix(cfg,method,values['test'],probability_horizon=h)
    return result
