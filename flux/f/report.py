"""Summarize the evaluated split, panel membership and provenance of each model."""
import html
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from common import file_digest, log, outdir
from data import load_meta
from metrics import macro_average


def common_comparison(frame):
    rows=[]
    for h,part in frame.groupby('horizon_years'):
        models=sorted(part.model.unique())
        for metric in ('recorded_auroc','recorded_ap','ipcw_auroc','ipcw_brier','harrell_c_index'):
            wide=part.pivot(index='endpoint',columns='model',values=metric).reindex(columns=models)
            shared=wide.replace([np.inf,-np.inf],np.nan).dropna()
            for model in models:
                rows.append(dict(horizon_years=h,metric=metric,model=model,common_endpoints=len(shared),mean=shared[model].mean()))
    return pd.DataFrame(rows)


def table(frame):
    return frame.to_html(index=False,border=0,float_format=lambda x:f'{x:.4f}',escape=True,na_rep='—')


def report(cfg, *, paper_index=None):
    from train import risk_cache, state_cache_key
    base=Path(cfg['paths']['out'])
    result=outdir(cfg,'main/results')
    dest=outdir(cfg,'main/reports')
    meta=load_meta(cfg)
    key=state_cache_key(cfg)
    readout_hash=file_digest(base/'main/results/models/readout.pt')
    methods=['flux',*cfg['baseline']['methods']]
    frames=[]
    accepted=[]
    for path in sorted(result.glob('*.metrics.tsv')):
        desc=path.with_name(path.name.replace('.metrics.tsv','.provenance.json'))
        if not desc.exists():
            continue
        p=json.loads(desc.read_text())
        method=p.get('model')
        if method not in methods+['rollout']:
            continue
        group='readout' if method=='flux' else 'rollout' if method=='rollout' else 'baseline'
        if (p['risk_key']!=key or p['config']!=cfg[group] or p.get('evaluation')!=cfg['evaluation']
                or p.get('endpoint_names')!=meta['endpoint_names'] or p.get('seed')!=cfg['seed']):
            continue
        if method=='flux' and p.get('readout_sha256')!=readout_hash:
            continue
        frames.append(pd.read_csv(path,sep='\t'))
        accepted.append(path.name)
    if 'flux.metrics.tsv' not in accepted:
        raise ValueError('No matching primary evaluation. Run evaluate before report.')
    frame=pd.concat(frames,ignore_index=True)
    common=common_comparison(frame[frame.model!='rollout'])
    common.to_csv(dest/'comparison_common.tsv',sep='\t',index=False)
    macro=macro_average(frame)
    macro.to_csv(dest/'summary.tsv',sep='\t',index=False)
    # One grouped event scan quantifies the effect of the old three-day cutoff.
    import sqlite3
    with sqlite3.connect(f"file:{base/'main/results/data/cohort.sqlite'}?mode=ro",uri=True) as db:
        rows=db.execute("""
            WITH days AS (
                SELECT eid,COUNT(DISTINCT CASE WHEN code != 'DEATH' THEN day END) AS n
                FROM events GROUP BY eid
            )
            SELECT p.split,
                   SUM(d.eid IS NOT NULL OR EXISTS (SELECT 1 FROM observations o WHERE o.eid=p.eid)),
                   SUM(COALESCE(d.n,0)=2), SUM(COALESCE(d.n,0)>=3)
            FROM participants p LEFT JOIN days d ON d.eid=p.eid GROUP BY p.split
        """)
        counts={split:(usable,two,three) for split,usable,two,three in rows}
    flow=[]
    for split,prepared in meta['split_counts'].items():
        a,m=risk_cache(cfg,split)
        usable,two,three=counts[split]
        flow.append(dict(split=split,prepared=prepared,usable_for_world=usable,exactly_two_clinical_days=two,
                         at_least_three_clinical_days=three,eligible_landmarks=m['n'],excluded=prepared-m['n']))
    flow=pd.DataFrame(flow)
    flow.to_csv(dest/'cohort_flow.tsv',sep='\t',index=False)
    evaluated=frame[(frame.model=='flux') & (frame.horizon_years==cfg['readout']['horizon_years'])].set_index('endpoint')
    description=json.loads((result/'flux.prediction.json').read_text())
    if description['states_key']!=key or description['readout_sha256']!=readout_hash:
        raise ValueError('The primary predictions do not match the fitted decoder.')
    panel=pd.DataFrame(dict(endpoint=meta['endpoint_names'],trained=description['eligible'])).set_index('endpoint')
    panel=panel.join(evaluated[['n','recorded_cases','recorded_controls','recorded_auroc','recorded_ap']])
    panel['estimable']=np.isfinite(panel.recorded_auroc)
    panel.to_csv(dest/'endpoint_panel.tsv',sep='\t')
    events=meta['summary'].get('events',0)
    synthetic=any('demo_inputs' in str(x) for x in meta['sources'])
    scope=dict(synthetic=synthetic,panel='training vocabulary' if cfg['data']['endpoints']=='all' else 'prespecified (user configuration)',
               events=events,features=len(meta['features']),context=len(meta['context_features']),endpoints=len(panel),
               event_source=meta['event_source_mode'],features_by_modality={m:len(v) for m,v in meta['modalities'].items()},
               admin_end=meta['admin_end'],landmark=cfg['data']['landmark'],minimum_clinical_days=cfg['data']['min_clinical_days'],
               warnings=meta.get('warnings',[]),metrics=accepted,risk_key=key,readout_sha256=readout_hash,
               blackout=False,rollout=False,missing_baselines=[m for m in methods[1:] if m+'.metrics.tsv' not in accepted])
    if (result/'blackout.json').exists():
        b=json.loads((result/'blackout.json').read_text())
        scope['blackout']=b['risk_key']==key and b['readout_sha256']==readout_hash and b['config']==cfg['blackout']
    if (result/'rollout.json').exists():
        r=json.loads((result/'rollout.json').read_text())
        scope['rollout']=r.get('risk_key')==key and r.get('config')==cfg['rollout']
    (dest/'scope.json').write_text(json.dumps(scope,indent=2,ensure_ascii=False))
    page='<meta charset="utf-8"><h1>HealthFlux held-out YGB evaluation</h1>'
    if synthetic: page+='<p>Synthetic test cohort</p>'
    page+=table(flow)+f'<p>{events:,} events, {len(panel)} candidate endpoints, {panel.estimable.sum()} evaluable five-year endpoints.</p>'
    index = paper_index or base/'paper/main/reports/index.html'
    link = html.escape(os.path.relpath(index, dest), quote=True)
    page+=table(common[common.horizon_years==5])+ f'<p><a href="{link}">Publication figures and tables</a></p>'
    page+='<pre>'+html.escape(json.dumps(scope,indent=2,ensure_ascii=False))+'</pre>'
    (dest/'report.html').write_text(page)
    log(f'Main report: {dest}/report.html')
    return scope
