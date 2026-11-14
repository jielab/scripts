"""Aggregate-only offline HTML report; no participant identifiers embedded."""
from __future__ import annotations
import html
import json
from pathlib import Path
import pandas as pd
import numpy as np
from common import file_digest, log, outdir
from data import load_meta
from metrics import macro_average


def table(frame):
    return frame.to_html(index=False,border=0,float_format=lambda x:f'{x:.4f}',escape=True,na_rep='—')


def common_comparison(frame):
    """Use identical endpoints across fitted models, separately for each metric."""
    rows=[]
    for horizon,part in frame.groupby('horizon_years'):
        models=sorted(part['model'].unique())
        for metric in ('recorded_auroc','recorded_ap','ipcw_auroc','ipcw_brier','harrell_c_index'):
            wide=part.pivot_table(index='endpoint',columns='model',values=metric,dropna=False)
            wide=wide.reindex(columns=models).replace([np.inf,-np.inf],np.nan).dropna()
            for model in models:
                rows.append(dict(horizon_years=horizon,metric=metric,model=model,
                                 common_endpoints=len(wide),mean=wide[model].mean()))
    return pd.DataFrame(rows)


def report(cfg):
    root=Path(cfg['paths']['out']); folder=outdir(cfg,'results'); meta=load_meta(cfg)
    from train import state_cache_key
    risk_key=state_cache_key(cfg)
    frames=[]; ignored=[]
    for path in sorted(folder.glob('*.metrics.tsv')):
        provenance=path.with_name(path.name.replace('.metrics.tsv','.provenance.json'))
        if not provenance.exists():
            ignored.append(path.name); continue
        p=json.loads(provenance.read_text()); method=p.get('model')
        group='readout' if method=='flux' else 'rollout' if method=='rollout' else 'baseline'
        current=(p.get('risk_key')==risk_key and p.get('endpoint_names')==meta['endpoint_names']
                 and p.get('config')==cfg[group] and p.get('seed')==cfg['seed']
                 and p.get('evaluation')==cfg['evaluation'])
        if method=='flux': current &= p.get('readout_sha256')==file_digest(root/'model/readout.pt')
        if group=='baseline': current &= method in cfg['baseline']['methods']
        if not current:
            ignored.append(path.name); continue
        frames.append(pd.read_csv(path,sep='\t'))
    if not frames: raise FileNotFoundError('Run evaluate before report')
    frame=pd.concat(frames,ignore_index=True)
    macro=macro_average(frame); macro.to_csv(folder/'summary.tsv',sep='\t',index=False)
    common=common_comparison(frame[frame['model']!='rollout'])
    common.to_csv(folder/'comparison_common.tsv',sep='\t',index=False)
    compact=macro[[x for x in ('model','horizon_years','endpoints','recorded_auroc','recorded_auroc_endpoints',
        'recorded_ap','ipcw_auroc','ipcw_auroc_endpoints','ipcw_brier') if x in macro]]
    warnings=meta.get('warnings',[])
    synthetic=any('demo_inputs' in str(x) for x in meta.get('sources',[]))
    label='SYNTHETIC SMOKE TEST — not a scientific benchmark' if synthetic else 'Independent UKB implementation — not the authors’ released model'
    parts=[f'<h1>HealthFlux / UKB</h1><p class="badge">{label}</p>',
        '<p>Latent state → continuous-time transition → observation update → frozen disease readout.</p>',
        '<h2>Run and input scope</h2>',
        f'<p>Event source: <strong>{html.escape(str(meta.get("event_source_mode")))}</strong>. '
        f'Train/validation/test participants: {html.escape(str(meta["split_counts"]))}. '
        f'{len(meta["vocab"]):,} tokens, {len(meta["features"]):,} numeric features, {len(meta["endpoint_names"]):,} endpoints. '
        f'Follow-up end: {html.escape(str(meta.get("admin_end")))}. Landmark: {html.escape(str(cfg["data"]["landmark"]))}.</p>',
        '<ul>'+''.join('<li>'+html.escape(str(w))+'</li>' for w in warnings)+'</ul>',
        '<h2>Held-out results</h2>',table(compact),
        '<p>Macro values average endpoints, not participants. Each metric reports its own eligible endpoint count. '
        'A five-year output used at ten years supplies ranking metrics only; its ten-year Brier score is omitted.</p>',
        '<h2>Matched endpoint comparison</h2>',table(common),
        '<p>The common panel keeps the same endpoints with finite values for every displayed model. '
        'All models use the same participant splits, landmark dates and endpoint-specific prevalent-case exclusions.</p>',
        '<h2>What the metrics mean</h2>',
        '<ul><li><b>Recorded AUROC/AP:</b> the paper-style recorded-diagnosis classification task. No recorded event by the horizon is a control, including incomplete follow-up.</li>'
        '<li><b>IPCW AUROC/Brier:</b> administrative follow-up and death censoring are explicit. The censoring distribution is fitted in the training risk set. '
        'Deaths censor nonfatal disease: this is not a competing-risk cumulative-incidence estimand.</li>'
        '<li><b>Weighted BCE readout:</b> outputs are not presumed calibrated. The calibration table and Brier score must be inspected before interpreting absolute risks.</li>'
        '<li><b>Rollout:</b> Monte Carlo conditional-incidence scores from generated events; not treatment effects or calibrated lifetime probabilities.</li></ul>',
        '<h2>Reproduction boundary</h2><p>This implementation reproduces the continuous-time latent-state core, '
        'chronological numeric/EHR inputs, frozen five-year readout, explicit baselines, no-ODE/no-update retraining, '
        'and event-feedback rollouts. It does not reproduce the authors’ participant split not supplied here, '
        'exact feature dictionary, pretrained raw-3D-MRI branch, external cohorts, or RCT-informed medication model. '
        'MRI inputs here are dated numeric image-derived phenotypes. Feature selection and architecture details '
        'are independently specified and recorded in logs/resolved.yml.</p>',
        '<h2>Files</h2><p>Full per-endpoint metrics: *.metrics.tsv. Calibration: *.calibration.tsv. '
        'Common-endpoint comparisons: comparison_common.tsv. Baseline fitting diagnostics: baseline.status.tsv. '
        'All arrays and checkpoints remain local; this page contains aggregates only.</p>']
    if ignored: parts+=['<p>Ignored stale or unversioned metric files: '+html.escape(', '.join(ignored))+'. Rerun their stage to refresh.</p>']
    rollout=folder/'rollout.json'
    if rollout.exists(): parts+=['<h2>Rollout diagnostics</h2>','<pre>'+html.escape(rollout.read_text())+'</pre>']
    ablations=[]
    for name in ('no_ode','no_update'):
        path=root/'ablation'/name/'results/flux.macro.tsv'
        if path.exists():
            a=pd.read_csv(path,sep='\t'); a.insert(0,'variant',name); ablations.append(a)
    if ablations: parts+=['<h2>Retrained ablations</h2>',table(pd.concat(ablations,ignore_index=True))]
    parts+=['<h2>Sources</h2><ul>'
        '<li><a href="https://www.medrxiv.org/content/10.64898/2026.09.19.26363460v3.full">HealthFlux v3 primary methods</a></li>'
        '<li><a href="https://github.com/jielab/scripts/tree/8edca20e15c6aaa9efa36a993f57862b3c3188a9">jielab/scripts, pinned source-layout reference</a></li></ul>']
    page='<!doctype html><html lang="en"><meta charset="utf-8"><title>HealthFlux UKB report</title><style>'
    page+='body{font:16px/1.6 system-ui,sans-serif;max-width:1260px;margin:40px auto;padding:0 24px;color:#16303b;background:#f8fafb}'
    page+='h1{font-size:40px;color:#074f58}h2{margin-top:38px}table{border-collapse:collapse;width:100%;font-size:13px;background:white}'
    page+='th,td{padding:9px 12px;border-bottom:1px solid #dce6e9;text-align:right}th:first-child,td:first-child{text-align:left}'
    page+='th{background:#e3eff0}a{color:#006675}.badge{background:#fff1d0;padding:12px 18px;border-radius:5px}pre{white-space:pre-wrap;background:#eef3f4;padding:18px}'
    page+='</style><main>'+''.join(parts)+'</main></html>'
    (folder/'report.html').write_text(page)
    log(f'Report: {folder/"report.html"}')
