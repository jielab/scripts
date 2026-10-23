#!/usr/bin/env python3
"""Compare saved CPU ABM/GPU results. No fit, projection, or CUDA execution."""
from pathlib import Path
import argparse
import gzip
import hashlib
import html
import json
import math
import os
import shutil
import tempfile
import uuid
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import abm
from c1_abm_compare import _metrics, benchmark_rows, conventional_scores
from c1_lock_roster import eid_hash
from le8_common import write_tables, sha256

VERSION = 'c1-final-saved-comparison-v3'


def read_json(path):
    return json.loads(Path(path).read_text())


def checked_completion(root):
    path = root / 'COMPLETE.json'
    if not path.is_file(): raise FileNotFoundError(f'Completed results unavailable: {root}')
    marker = read_json(path)
    if marker.get('status') != 'completed' or not marker.get('files'):
        raise ValueError(f'Invalid completion marker: {path}')
    for name, digest in marker['files'].items():
        p = root / name
        if not p.resolve().is_relative_to(root.resolve()): raise ValueError('Invalid receipt path')
        if not p.is_file() or sha256(p) != digest: raise ValueError(f'Completed file missing or changed: {p}')
    return [path, *(root / name for name in marker['files'])]


def ids(values, label):
    x = np.asarray(values, str)
    if len(np.unique(x)) != len(x) or np.any(x == ''): raise ValueError(f'Duplicate/empty IDs: {label}')
    return x


def boolean(values):
    x = pd.Series(values)
    if not x.astype(str).str.lower().isin(['true','false','1','0']).all(): raise ValueError('Invalid saved selection mask')
    return x.astype(str).str.lower().isin(['true','1']).to_numpy()


def paired_cohort(cpu, truth, gpu, km, conventional=None):
    ci = ids(cpu['eid'], 'CPU'); ti = ids(truth['eid'], 'truth'); gi = ids(gpu.eid, 'GPU')
    if not np.array_equal(ci, ti): raise ValueError('CPU predictions/truth ID order mismatch')
    index = pd.Index(ci).get_indexer(gi); keep = index >= 0
    ix = index[keep]; g = gpu.loc[keep].reset_index(drop=True)
    if not len(ix): raise ValueError('No common test participants')
    if not np.array_equal(np.asarray(cpu['group'], str)[ix], g.family.astype(str)):
        raise ValueError('CPU/GPU family assignments differ')
    t, e = np.asarray(truth['time'])[ix], np.asarray(truth['event'])[ix]
    if not np.allclose(t, g.time, rtol=0, atol=1e-8, equal_nan=False) or not np.array_equal(e,g.event):
        raise ValueError('CPU/GPU endpoint or follow-up differs')
    y,w = km.labels_weights(t,e)
    if 'y' in g and not np.array_equal(y,g.y): raise ValueError('CPU/GPU horizon labels differ')
    probabilities = {f'CPU_{name}':np.asarray(p,float)[ix] for name,p in cpu['probability'].items()}
    probabilities['GPU_Transformer'] = g.probability.to_numpy(float)
    if conventional is not None and not conventional.empty:
        probabilities.update(conventional_scores(conventional, g.eid.to_numpy(), t, e))
    finite = np.logical_and.reduce([np.isfinite(v) for v in probabilities.values()])
    if not finite.any(): raise ValueError('No shared finite predictions across models')
    cpu_qc = np.asarray(cpu.get('diagnostics', {}).get('technical_qc'))
    if cpu_qc.shape != ci.shape or not np.isin(cpu_qc, [True, False, 0, 1]).all():
        raise ValueError('Frozen CPU technical-QC mask is required for the shared comparison')
    qc = cpu_qc.astype(bool)[ix]
    # Legacy GPU exports lack a QC column. Their shared estimand uses the
    # documented CPU reference QC; do not infer QC from zero IPCW weights,
    # since an early-censored participant also has zero horizon weight.
    gpu_qc_status = 'not_saved_using_frozen_CPU_reference_QC'
    if 'technical_qc' in g:
        qc &= boolean(g.technical_qc)
        gpu_qc_status = 'saved_GPU_QC_intersected_with_CPU_QC'
    evaluated = finite & qc
    if not evaluated.any(): raise ValueError('No shared QC-passing finite predictions')
    for name,v in probabilities.items():
        if not name.startswith('C1_') and np.any((v[finite]<0)|(v[finite]>1)): raise ValueError('Invalid risk probability: '+name)
    requested_masks = dict(all_common=np.ones(len(ix), bool),
                           CPU_selected=np.asarray(cpu['selected'],bool)[ix],
                           GPU_selected=boolean(g.selected))
    requested_masks['both_selected'] = requested_masks['CPU_selected'] & requested_masks['GPU_selected']
    masks = {name: mask & evaluated for name, mask in requested_masks.items()}
    audit = dict(CPU_test_N=len(ci),GPU_test_N=len(gi),common_test_N=len(ix),
        CPU_only_N=len(ci)-len(ix),GPU_only_N=len(gi)-len(ix),excluded_nonfinite_N=int((~finite).sum()),
        excluded_QC_N=int((~qc).sum()),excluded_QC_or_nonfinite_N=int((~evaluated).sum()),
        comparison_N=int(evaluated.sum()),family_N=len(np.unique(g.family.to_numpy()[evaluated])),
        common_test_roster_sha256=eid_hash(g.eid),comparison_roster_sha256=eid_hash(g.eid[evaluated]),
        QC_rule='frozen CPU reference QC intersected with saved GPU QC when available',GPU_QC_status=gpu_qc_status,
        evaluation_weights='same frozen CPU training-derived censoring KM for every model',
        selection_scope='Every row compares models on exactly the same saved subset; no selection refitting',
        family_match=True,endpoint_match=True)
    if 'weight' in g:
        audit['saved_GPU_vs_shared_weight_max_difference'] = float(np.max(np.abs(g.weight.to_numpy()-w)))
    return dict(eid=g.eid.to_numpy(),time=t,event=e,y=y,w=w,groups=g.family.to_numpy(),
                probability=probabilities,masks=masks,requested_masks=requested_masks,qc_mask=qc,audit=audit)


def weighted_metrics(y, w, p, probability=True):
    use = (w > 0) & np.isfinite(p)
    result = dict(AUC_IPCW=np.nan, Brier_IPCW=np.nan, LogLoss_IPCW=np.nan)
    if not use.any(): return result
    y, w, p = y[use], w[use], p[use]
    # Conventional Cox scores are unrestricted risks, not probabilities.
    if len(np.unique(y)) == 2:
        result['AUC_IPCW'] = float(roc_auc_score(y, p, sample_weight=w))
    if probability:
        p = np.clip(p, 1e-7, 1-1e-7)
        result.update(Brier_IPCW=float(np.average((y-p)**2, weights=w)),
            LogLoss_IPCW=float(np.average(-(y*np.log(p)+(1-y)*np.log(1-p)), weights=w)))
    return result


def performance(data, km, bootstrap, seed):
    if int(bootstrap) != bootstrap or bootstrap < 0:
        raise ValueError('Bootstrap replicates must be a nonnegative integer')
    bootstrap = int(bootstrap)
    minimum_valid = max(20, math.ceil(.8 * bootstrap))
    probabilities = {name: np.asarray(p, float) for name, p in data['probability'].items()}
    if not probabilities or any(p.shape != (len(data['eid']),) for p in probabilities.values()):
        raise ValueError('Every model must provide one score per comparison participant')
    common_finite = np.logical_and.reduce([np.isfinite(p) for p in probabilities.values()])
    grouping_unit = data.get('grouping_unit', 'individual_unverified')
    if grouping_unit not in ('family', 'individual_unverified'):
        raise ValueError('Invalid comparison grouping unit')
    rows, deltas = [], []
    for subset, mask in data['masks'].items():
        # Enforce one mask even for callers supplying unfiltered subset masks.
        requested = np.asarray(data.get('requested_masks', {}).get(subset, mask), bool)
        mask = np.asarray(mask, bool)
        if requested.shape != common_finite.shape or mask.shape != common_finite.shape or np.any(mask & ~requested):
            raise ValueError('Invalid requested/evaluated comparison mask: ' + subset)
        mask = mask & common_finite
        y, w = data['y'][mask], data['w'][mask]
        t, e = data['time'][mask], data['event'][mask]
        p = {k:v[mask] for k,v in probabilities.items()}
        _, codes = np.unique(data['groups'][mask], return_inverse=True)
        families = len(np.unique(codes))
        family_available = data.get('family_available', True)
        common = dict(trait=data.get('trait', ''), layer=data.get('layer', ''), subset=subset,
            comparison_N=int(mask.sum()), requested_N=int(requested.sum()), evaluated_N=int(mask.sum()),
            excluded_N=int(requested.sum() - mask.sum()),
            excluded_nonfinite_N=int(np.sum(requested & ~common_finite)),
            excluded_QC_N=int(np.sum(requested & ~np.asarray(data.get('qc_mask', np.ones(len(requested), bool)), bool))),
            events=int(e.sum()), bootstrap_cluster_N=families if family_available else np.nan,
            families=families if family_available and grouping_unit == 'family' else np.nan,
            requested_eid_hash=eid_hash(np.asarray(data['eid'])[requested]),
            comparison_eid_hash=eid_hash(np.asarray(data['eid'])[mask]),
            outer_roster_sha256=data.get('outer_roster_sha256', ''), same_people=True,
            grouping_unit=grouping_unit, Harrell_time_scope='full_follow_up_strict_time',
            IPCW_horizon=km.horizon)
        measures = {}
        for name, risk in p.items():
            measure = _metrics(t, e, risk, km, not name.startswith('C1_'))
            measures[name] = measure
            rows.append(dict(**common, model=name, model_family=name.split('_')[0], **measure,
                status='evaluated' if np.isfinite(measure['AUC_IPCW']) and np.isfinite(measure['Harrell_C']) else 'insufficient_information'))
        pairs = [(a,b) for a in p for b in p
                 if (a.startswith('CPU_') and b.startswith('C1_')) or
                    (a == 'GPU_Transformer' and b != a)]
        # Planned CPU contrasts share the exact same subset as every external
        # benchmark; cleaning effects need both co-primary metric intervals too.
        cpu_pairs = [('elasticnet_clean', 'elasticnet_full'), ('elasticnet_clean', 'elasticnet_random'),
                     ('twin_clean', 'twin_full'), ('twin_clean', 'twin_random'),
                     ('twin_integrated_clean', 'twin_integrated_full'), ('twin_integrated_clean', 'twin_integrated_random'),
                     ('twin_clean', 'twin_integrated_clean'), ('fallback_policy', 'elasticnet_full'),
                     ('lightgbm_clean', 'lightgbm_full'), ('knn_unweighted_clean', 'knn_unweighted_full'),
                     ('mixture_experts', 'elasticnet_full')]
        pairs.extend(('CPU_' + a, 'CPU_' + b) for a, b in cpu_pairs if 'CPU_' + a in p and 'CPU_' + b in p)
        metric_names = ('Harrell_C','AUC_IPCW','Brier_IPCW','LogLoss_IPCW')
        samples = {pair:{m:[] for m in metric_names} for pair in pairs}
        # Fixed scores/time never change during bootstrap: cache rank sweeps,
        # but apply each fresh family draw identically to every model.
        harrell_plans = {name:abm._harrell_bootstrap_plan(t, e, risk) for name,risk in p.items()} if bootstrap and families and family_available else {}
        rng = np.random.default_rng(seed)
        for _ in range(bootstrap if families and family_available else 0):
            multiplicity = np.bincount(rng.integers(families, size=families), minlength=families)[codes]
            estimates = {}
            for name, risk in p.items():
                estimates[name] = weighted_metrics(y, w*multiplicity, risk, not name.startswith('C1_'))
                estimates[name]['Harrell_C'] = harrell_plans[name](multiplicity)[0]
            for (a,b), metrics in samples.items():
                for metric in metrics:
                    difference = estimates[a][metric] - estimates[b][metric]
                    if np.isfinite(difference): metrics[metric].append(difference)
        for (a,b), metrics in samples.items():
            for metric, values in metrics.items():
                enough = len(values) >= minimum_valid
                ci = np.quantile(values,[.025,.975]) if enough else [np.nan,np.nan]
                difference = measures[a][metric] - measures[b][metric]
                inference_status = ('insufficient_information' if not np.isfinite(difference) else
                    'grouping_unavailable' if not family_available else
                    'bootstrap_not_requested' if not bootstrap else
                    'evaluated' if enough else 'insufficient_bootstrap')
                deltas.append(dict(**common, model_A=a, model_B=b, CPU_model=b,
                    metric=metric, N=int(mask.sum()), delta=difference,
                    GPU_minus_CPU=difference if a == 'GPU_Transformer' and b.startswith('CPU_') else np.nan,
                    lower95=ci[0], upper95=ci[1], bootstrap_requested=bootstrap, bootstrap_valid=len(values),
                    bootstrap_minimum_valid=minimum_valid,
                    bootstrap_unit=grouping_unit if family_available else 'unavailable',
                    better_GPU_direction='positive' if metric in ('AUC_IPCW','Harrell_C') else 'negative',
                    status=inference_status, inference_status=inference_status,
                    point_status='evaluated' if np.isfinite(difference) else 'insufficient_information',
                    inference='exploratory paired fixed-prediction group bootstrap; no refitting'))
    return pd.DataFrame(rows), pd.DataFrame(deltas)


def primary_tables(perf, deltas):
    performance_table = perf.rename(columns={'AUC_IPCW':'AUC_IPCW_10y',
        'Brier_IPCW':'Brier_IPCW_10y','LogLoss_IPCW':'LogLoss_IPCW_10y'})
    rows = []
    if not deltas.empty:
        for _, d in deltas.groupby(['subset','model_A','model_B'], sort=False):
            first = d.iloc[0]
            row = {k:first[k] for k in ('trait','layer','subset','model_A','model_B',
                'comparison_N','requested_N','evaluated_N','excluded_N','excluded_nonfinite_N',
                'events','families','bootstrap_cluster_N','requested_eid_hash','comparison_eid_hash',
                'outer_roster_sha256','same_people','bootstrap_unit','Harrell_time_scope','IPCW_horizon')}
            row['status'] = 'evaluated'
            for metric, label, ci_label in [('Harrell_C','Harrell_C','Harrell_C'),('AUC_IPCW','AUC_IPCW_10y','AUC')]:
                m = d.loc[d.metric == metric].iloc[0]
                row['delta_'+label] = m['delta']
                row[ci_label+'_lower95'], row[ci_label+'_upper95'] = m.lower95, m.upper95
                row[ci_label+'_bootstrap_valid'] = m.bootstrap_valid
                row[ci_label+'_inference_status'] = m.inference_status
                row[ci_label+'_point_status'] = m.point_status
                if m.status != 'evaluated': row['status'] = m.status
            rows.append(row)
    return dict(same_subset_primary_performance=performance_table,
                same_subset_primary_deltas=pd.DataFrame(rows))


def numerical_parity(cpu, saved, network):
    ci=ids(cpu['eid'],'CPU'); gi=ids(saved['eid'],'parity')
    ix=pd.Index(gi).get_indexer(ci)
    if len(ci)!=len(gi) or np.any(ix<0): raise ValueError('Parity test roster differs from CPU')
    if not np.array_equal(np.asarray(saved['group'],str)[ix],np.asarray(cpu['group'],str)):
        raise ValueError('Parity family mismatch')
    rows=[]
    exact={'eligible_match_count','matched_count','family_count','supported','technical_qc'}
    for kind in ('probability','diagnostics','masks'):
        actual=saved.get(kind,{})
        if not actual: raise ValueError('Missing saved parity fields: '+kind)
        if kind=='probability' and not {'twin_clean','twin_full'} <= set(actual):
            raise ValueError('Saved parity lacks core twin model predictions')
        for field,a in actual.items():
            if field not in cpu[kind]: raise ValueError('Unknown parity field: '+field)
            a=np.asarray(a)[ix];b=np.asarray(cpu[kind][field]);delta=np.abs(a.astype(float)-b.astype(float))
            tolerance=0. if kind=='masks' or field in exact else float(network['parity_probability_tolerance' if kind=='probability' else 'parity_continuous_tolerance'])
            finite=np.isfinite(delta).all();maximum=float(delta.max()) if finite else float('inf')
            rows.append(dict(kind=kind,field=field,N=len(ci),max_abs_difference=maximum,
                tolerance=tolerance,mismatches=int(np.count_nonzero(a!=b)),passed=bool(finite and maximum<=tolerance)))
    a=np.asarray(saved['selected'])[ix];b=np.asarray(cpu['selected'])
    rows.append(dict(kind='selector',field='selected',N=len(ci),max_abs_difference=float(np.any(a!=b)),
                     tolerance=0.,mismatches=int(np.count_nonzero(a!=b)),passed=bool(np.array_equal(a,b))))
    return pd.DataFrame(rows)


def timing_tables(gpu_root, transformer):
    rows=[]
    parity=gpu_root/'parity/timing.tsv'
    if parity.is_file():
        for row in pd.read_csv(parity,sep='\t').to_dict('records'):
            rows.append(dict(category='same_algorithm_matching',status='recorded',**row,
                interpretation='Recorded warmed-up median on identical queries; matching speedup only, not training speedup'))
    steps=transformer/'performance/step_timings.tsv'
    if steps.is_file():
        for row in pd.read_csv(steps,sep='\t').to_dict('records'):
            rows.append(dict(category='different_model_step',status='recorded',**row,
                interpretation='Recorded distinct-model core inference/training; not algorithm-equivalent speedup'))
    history=transformer/'training_history.tsv'
    if history.is_file():
        d=pd.read_csv(history,sep='\t')
        if 'seconds' in d:
            rows.append(dict(category='GPU_recorded_epochs',status='recorded',seconds=float(d.seconds.sum()),
                interpretation='Sum of recorded epoch times; excludes preparation, reporting and unrecorded work'))
    report=gpu_root/'GPU_REPORT.json'
    if report.is_file():
        for key,value in read_json(report).get('timing',{}).items():
            if key.endswith('_seconds') and value is not None:
                rows.append(dict(category='GPU_workflow',step=key,seconds=value,status='recorded',
                    interpretation='Recorded workflow wall time; not comparable to CPU training alone'))
    rows.append(dict(category='CPU_training',status='not_recorded',seconds=np.nan,
                     interpretation='CPU artifacts have no reliable total training timer; no training speedup inferred'))
    if not parity.is_file(): rows.append(dict(category='same_algorithm_matching',status='unavailable',interpretation='No saved CPU/GPU matching benchmark'))
    if not steps.is_file(): rows.append(dict(category='different_model_step',status='unavailable',interpretation='No saved paired inference timings'))
    return pd.DataFrame(rows)


def load_inputs(scope):
    cpu_root=scope/'c1_abm/abm_twin';gpu_root=scope/'c1_abm_GPU'
    files=checked_completion(cpu_root)+checked_completion(gpu_root)
    transformer=gpu_root/'transformer' if (gpu_root/'transformer').is_dir() else gpu_root
    paths=[cpu_root/'fit/model_LOCAL_ONLY.joblib',cpu_root/'test_predictions_LOCAL_ONLY.joblib',
           cpu_root/'test_truth_LOCAL_ONLY.npz',cpu_root/'host_contract.json',gpu_root/'host_contract.json',
           transformer/'test_predictions_LOCAL_ONLY.tsv.gz',transformer/'MODEL_FROZEN.json']
    for p in paths:
        if not p.is_file(): raise FileNotFoundError('Saved comparison input unavailable: '+str(p))
    files+=paths
    ca,ga=read_json(cpu_root/'host_contract.json'),read_json(gpu_root/'host_contract.json')
    for key in ('outer_roster_sha256','group_map_sha256'):
        if not ca.get(key) or ca.get(key)!=ga.get(key): raise ValueError('CPU/GPU protocol mismatch: '+key)
    if ca.get('grouping_unit')!='family' or ga.get('grouping_unit')!='family': raise ValueError('Verified family grouping required')
    # New GPU runs retain exact CPU artifact hashes. Historical independent
    # Transformer runs are instead verified by endpoint/roster/family alignment.
    for name,digest in ga.get('cpu_source',{}).get('files',{}).items():
        path=cpu_root/name
        if not path.is_file() or sha256(path)!=digest: raise ValueError('GPU used different CPU artifacts: '+name)
        files.append(path)
    model=abm.load_model(paths[0]);cpu=joblib.load(paths[1])
    if cpu['frozen_model_hash']!=abm.model_hash(model): raise ValueError('CPU frozen model/prediction hash mismatch')
    with np.load(paths[2],allow_pickle=False) as t: truth={k:t[k] for k in t.files}
    gpu=pd.read_csv(paths[-2],sep='\t',dtype={'eid':str,'family':str})
    marker=read_json(paths[-1])
    if not np.isclose(marker['horizon'],model.config.horizon,rtol=0,atol=1e-10): raise ValueError('CPU/GPU prediction horizons differ')
    optional=[gpu_root/'parity/parity_predictions_LOCAL_ONLY.joblib',gpu_root/'parity/timing.tsv',
        gpu_root/'parity/timing_summary.json',gpu_root/'GPU_REPORT.json',
        transformer/'performance/step_timings.tsv',transformer/'performance/timing_scope.json',transformer/'training_history.tsv']
    files += [p for p in optional if p.is_file()]
    benchmark_path = scope/'c1_correlate/prediction/_LOCAL_ONLY/benchmark.json.gz'
    conventional = None
    if benchmark_path.is_file():
        with gzip.open(benchmark_path, 'rt') as stream:
            benchmark = json.load(stream)
        conventional = benchmark_rows(benchmark, scope.parent.name, scope.name,
                                     model.config.horizon, ca['outer_roster_sha256'])
        files.append(benchmark_path)
    ga = dict(ga, conventional_status='available' if conventional is not None and not conventional.empty else 'unavailable')
    fingerprint={str(p):sha256(p) for p in sorted(set(files))}
    return model,cpu,truth,gpu,ga,transformer,fingerprint,conventional


def publish(target, tables, metadata):
    target.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='_comparison_',dir=target.parent))
    try:
        write_tables(stage,tables,'CPU/GPU saved comparison',target.parent.name,metadata['trait'],workbook='CPU_GPU_comparison.xlsx')
        (stage/'comparison.json').write_text(json.dumps(metadata,indent=2,allow_nan=False))
        note=('All performance rows use the same participants within each subset and the frozen CPU censoring model. '
              'GPU Transformer is a different model; numerical equivalence applies only to the saved frozen-model parity branch. '
              'Timing entries retain their original scope. Missing timing is not estimated. No model was trained or projected.')
        body=''.join('<h2>'+html.escape(name)+'</h2>'+d.to_html(index=False,na_rep='Unavailable',escape=True) for name,d in tables.items())
        (stage/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>CPU/GPU comparison</title>'
            '<style>body{font:15px system-ui;margin:32px}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border:1px solid #ddd}h2{margin-top:32px}</style>'
            '<h1>CPU ABM / GPU comparison</h1><p>'+note+'</p>'+body)
        for path,digest in metadata.get('sources',{}).items():
            if sha256(path)!=digest: raise ValueError('Source changed during comparison: '+path)
        manifest={p.name:sha256(p) for p in stage.iterdir() if p.is_file()}
        (stage/'COMPLETE.json').write_text(json.dumps(dict(status='completed',files=manifest),indent=2))
        backup=None
        if target.exists():
            backup=target.parent/'_history'/(target.name+'_'+uuid.uuid4().hex[:10]);backup.parent.mkdir(exist_ok=True);target.rename(backup)
        try: stage.rename(target)
        except BaseException:
            if backup: backup.rename(target)
            raise
    finally:
        if stage.exists(): shutil.rmtree(stage)


def compare_scope(root, trait, layer, bootstrap=200, seed=2026, check=False):
    scope=Path(root)/trait/layer;target=scope/'c1_abm_comparison'
    model,cpu,truth,gpu,contract,transformer,sources,conventional=load_inputs(scope)
    data=paired_cohort(cpu,truth,gpu,model.km,conventional)
    data.update(trait=trait, layer=layer, outer_roster_sha256=contract['outer_roster_sha256'],
                grouping_unit=getattr(model.config, 'grouping_unit', 'individual_unverified'))
    metadata=dict(version=VERSION,trait=trait,layer=layer,bootstrap=bootstrap,seed=seed,
        sources=sources,code_sha256={name:sha256(Path(__file__).with_name(name)) for name in
            ('abm.py','c1_abm_compare.py','c1_abm_gpu_compare.py')},
        conventional_status=contract['conventional_status'],
        co_primary_metrics='Harrell C-index over follow-up; 10-year IPCW AUC',
        secondary_metrics='horizon-truncated Uno C; 10-year IPCW Brier; 10-year IPCW log loss',
        caution='selected-vs-all performance is descriptive selection composition, not algorithm superiority',
        training_performed=False,inference_performed=False)
    if check: return dict(status='inputs_verified',**data['audit'])
    if (target/'comparison.json').is_file() and read_json(target/'comparison.json')==metadata:
        try:
            checked_completion(target)
            print(f'[LE8 comparison] SKIP {trait}/{layer}: verified unchanged report',flush=True)
            return target
        except (ValueError,FileNotFoundError): pass
    metrics,deltas=performance(data,model.km,bootstrap,seed)
    parity_path=scope/'c1_abm_GPU/parity/parity_predictions_LOCAL_ONLY.joblib'
    if parity_path.is_file():
        if not contract.get('cpu_source'): raise ValueError('Parity lacks frozen CPU source receipt')
        parity=numerical_parity(cpu,joblib.load(parity_path),contract['gpu_config'])
        parity_status='passed' if parity.passed.all() else 'failed'
    else:
        parity=pd.DataFrame([dict(status='unavailable',reason='No saved frozen-model CPU/GPU parity predictions; Transformer predictions are not an equivalence test')])
        parity_status='unavailable'
    tables=dict(**primary_tables(metrics,deltas),status=pd.DataFrame([dict(status='evaluated',trait=trait,layer=layer,parity_status=parity_status,
        exploratory=True,conventional_status=contract['conventional_status'],
        co_primary_metrics=metadata['co_primary_metrics'],secondary_metrics=metadata['secondary_metrics'],
        caution=metadata['caution'],training_performed=False,inference_performed=False)]),cohort=pd.DataFrame([data['audit']]),
        performance=metrics,paired_differences=deltas,numerical_parity=parity,
        timings=timing_tables(scope/'c1_abm_GPU',transformer))
    publish(target,tables,metadata)
    print(f'[LE8 comparison] DONE {trait}/{layer}: {target}/CPU_GPU_comparison.xlsx',flush=True)
    return target


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-root',type=Path,default=Path('/mnt/d/analysis/le8'))
    p.add_argument('--Y',required=True);p.add_argument('--biom',default='prot,met')
    p.add_argument('--bootstrap',type=int,default=200);p.add_argument('--seed',type=int,default=2026)
    p.add_argument('--check',action='store_true');a=p.parse_args()
    if a.bootstrap<0: p.error('--bootstrap must be nonnegative')
    pending=[]
    from threadpoolctl import threadpool_limits
    with threadpool_limits(2):
        for trait in a.Y.split(','):
            for layer in a.biom.split(','):
                try: compare_scope(a.analysis_root,trait,layer,a.bootstrap,a.seed,a.check)
                except FileNotFoundError as exc:
                    pending.append(str(exc));print('[LE8 comparison] PENDING '+str(exc),flush=True)
    if pending: raise SystemExit(2)


if __name__=='__main__': main()
