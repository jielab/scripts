"""Paper-oriented aggregate tables and PNG figures. Reads completed scientific outputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from common import outdir, log
from report import report, table


def interval(values, rng):
    x=np.asarray(values,dtype=float)
    x=x[np.isfinite(x)]
    if not len(x):
        return 0,np.nan,np.nan,np.nan
    means=np.concatenate([rng.choice(x,size=(100,len(x))).mean(axis=1) for _ in range(100)])
    return len(x),x.mean(),*np.quantile(means,[.025,.975])


def paper(cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(cfg['paths']['out'])
    src=root/'main/results'
    main=root/'main/reports'
    dest=outdir(cfg,'paper/main/results')
    web=outdir(cfg,'paper/main/reports')
    scope=report(cfg)
    rng=np.random.default_rng(cfg['seed'])
    parts=['<h1>HealthFlux · UK Biobank</h1>']
    if scope['synthetic']:
        parts+=['<b>Synthetic cohort · workflow check only</b>']
    parts+=[f"<p>Endpoint panel: {scope['panel']}. Source: {scope['event_source']}. "
            f"Administrative end: {scope['admin_end']}. The study benchmark used 1,256 diagnoses and death.</p>"]
    def section(title, df, name, display=None):
        df.to_csv(dest/(name+'.tsv'),sep='\t',index=False)
        parts.extend([f'<h2>{title}</h2>',f'<p><a href="../results/{name}.tsv">Download table</a></p>',
                      table(df if display is None else display)])
    def figure(name, fig):
        fig.tight_layout()
        fig.savefig(web/(name+'.png'),dpi=155)
        plt.close(fig)
        parts.append(f'<img src="{name}.png" alt="{name}">')
    flow=pd.read_csv(main/'cohort_flow.tsv',sep='\t')
    section('Participant eligibility',flow,'cohort_flow')
    summary=pd.read_csv(main/'comparison_common.tsv',sep='\t')
    section('Five-year comparison on common endpoints',summary,'common_endpoints',summary.query('horizon_years==5'))
    frame=pd.concat([pd.read_csv(src/p,sep='\t') for p in scope['metrics']],ignore_index=True)
    paired=[]
    scores=frame[(frame.horizon_years==5)&(frame.model!='rollout')]
    for metric in ('recorded_auroc','recorded_ap'):
        wide=scores.pivot(index='endpoint',columns='model',values=metric).dropna()
        for comparator in wide.columns.difference(['flux']):
            paired.append([comparator,metric,*interval(wide['flux']-wide[comparator],rng)])
    df=pd.DataFrame(paired,columns=['comparator','metric','common_endpoints','mean_flux_minus_comparator','lower','upper'])
    section('Paired five-year differences across common endpoints',df,'paired_differences')
    single=frame.query('model=="flux" and horizon_years==5').copy()
    single['chapter']=single.endpoint.str.extract(r'^([A-Z])[0-9]{2}$',expand=False).fillna(single.endpoint)
    rarity=[]
    for min_cases in (1,20):
        supported=single.loc[single.recorded_cases.ge(min_cases),'endpoint']
        sub=scores[scores.endpoint.isin(supported)]
        from report import common_comparison
        shared=common_comparison(sub) if len(sub) else pd.DataFrame()
        if len(shared):
            shared=shared[shared.metric.isin(['recorded_auroc','recorded_ap'])].copy()
            shared.insert(0,'min_test_cases',min_cases)
            rarity.append(shared)
    if rarity:
        section('Sensitivity to rare test outcomes',pd.concat(rarity,ignore_index=True),'case_support')
    per=single.groupby('chapter').agg(candidate=('endpoint','size'),estimable=('recorded_auroc','count'),
        mean_auroc=('recorded_auroc','mean'),median_auroc=('recorded_auroc','median'),mean_ap=('recorded_ap','mean')).reset_index()
    section('Five-year discrimination by ICD-10 initial letter',per,'chapter')
    fig,ax=plt.subplots(figsize=(8,3.5))
    values=per.dropna(subset=['mean_auroc'])
    ax.scatter(values.chapter,values.mean_auroc,color='#176585')
    ax.set(ylabel='Mean recorded AUROC',xlabel='ICD-10 initial letter or endpoint');ax.tick_params(axis='x',labelrotation=55)
    figure('five_year_discrimination',fig)
    enrich=pd.read_csv(src/'flux.enrichment.tsv',sep='\t')
    focus=['DEATH','I21','I25','I50','I63','E11','N18','J44']
    section('Observed five-year incidence in the top-risk groups',enrich,'enrichment',
            enrich[enrich.endpoint.isin(focus)&enrich.top_fraction.eq(.01)])
    cal=pd.read_csv(src/'flux.validation_calibration.tsv',sep='\t')
    ok=cal[cal.status=='ok']
    section('IPCW calibration fitted on validation and scored on test',cal,'calibration',
            ok.reindex(columns=['endpoint','raw_ipcw_brier','calibrated_ipcw_brier']).query('endpoint in @focus'))
    if len(ok):
        parts.append(f'<p>Paired mean test Brier across {len(ok)} endpoints: raw {ok.raw_ipcw_brier.mean():.5f}, calibrated {ok.calibrated_ipcw_brier.mean():.5f}.</p>')
    if scope['blackout']:
        black=pd.read_csv(src/'blackout.tsv',sep='\t')
        rows=[]
        for lag,g in black.groupby('lag_years'):
            rows.append([lag,*interval(g.auroc,rng)])
        df=pd.DataFrame(rows,columns=['lag_years','endpoints','mean_auroc','lower','upper'])
        section('Five-year prediction after an observation blackout',df,'blackout')
        if len(df):
            fig,ax=plt.subplots(figsize=(6,3.5));ax.plot(df.lag_years,df.mean_auroc,'o-',color='#176585')
            ax.fill_between(df.lag_years,df.lower,df.upper,color='#176585',alpha=.2)
            ax.set(xlabel='Years since last observation',ylabel='Next-five-year recorded AUROC')
            figure('blackout',fig)
        else:
            parts.append('<p>No local endpoints met the prespecified stress-query support rule.</p>')
    else: parts.append('<p>Blackout results unavailable for this run.</p>')
    if scope['rollout']:
        status=json.loads((src/'rollout.json').read_text())
        roll=frame[frame.model=='rollout']
        rows=[]
        for h,info in status['horizon_status'].items():
            vals=roll.loc[roll.horizon_years==float(h),'ipcw_auroc']
            rows.append([float(h),info['status'],info['completed_participants'],*interval(vals,rng)])
        df=pd.DataFrame(rows,columns=['horizon','simulation_status','participants_complete','endpoints','ipcw_auroc','lower','upper'])
        df['evaluation_status']=np.where(df.endpoints>0,'estimable','no estimable endpoint')
        section('Recursive event simulation from the observed landmark',df,'recursive')
        if df.ipcw_auroc.notna().any():
            fig,ax=plt.subplots(figsize=(7,3.5));ax.plot(df.horizon,df.ipcw_auroc,'o-',color='#257152',ms=3)
            ax.fill_between(df.horizon,df.lower,df.upper,color='#257152',alpha=.2)
            ax.set(xlabel='Simulated horizon (years)',ylabel='Mean cumulative/dynamic IPCW AUROC')
            figure('recursive',fig)
    else: parts.append('<p>Recursive simulation results unavailable for this run.</p>')
    timings=[]
    def append_dynamics(folder,variant):
        if (folder/'next_event.tsv').exists():
            d=pd.read_csv(folder/'next_event.tsv',sep='\t')
            d.insert(0,'variant',variant)
            timings.append(d)
    append_dynamics(src,'full')
    ab=[]
    for variant in ('no_ode','no_update'):
        r=root/'ablation'/variant/'main/results'
        stamp=r/'ablation.json'
        if stamp.exists() and json.loads(stamp.read_text()).get('parent_risk_key')==scope['risk_key']:
            variant_cfg = {**cfg, 'paths': {**cfg['paths'], 'out': str(root/'ablation'/variant)},
                           'model': {**cfg['model'], variant: True}}
            report(variant_cfg, paper_index=web/'index.html')
            append_dynamics(r,variant)
            e=pd.read_csv(r/'flux.metrics.tsv',sep='\t').query('horizon_years==5')
            m=single[['endpoint','recorded_auroc']].merge(e[['endpoint','recorded_auroc']],on='endpoint',suffixes=('_full','_variant')).dropna()
            delta=m.recorded_auroc_variant-m.recorded_auroc_full
            ab.append([variant,*interval(delta,rng)])
    if timings:
        section('Held-out next-event cross-entropy and timing error',pd.concat(timings,ignore_index=True),'event_dynamics')
    if ab:
        section('Retrained architecture ablations: paired change in five-year AUROC',
                pd.DataFrame(ab,columns=['variant','endpoints','mean_delta','lower','upper']),'ablations')
    if scope['missing_baselines']:parts.append('<p>Unavailable baseline results: '+', '.join(scope['missing_baselines'])+'</p>')
    parts+=['<p>Intervals resample endpoints 10,000 times. Error bands describe between-endpoint variation. A five-year decoder used at ten years supplies a ranking sensitivity only. Recursive results use the particle simulator and endpoint-specific censoring. Top-group incidence uses recorded diagnoses.</p>',
        '<p>Numeric imaging phenotypes and locally fitted models. External cohorts, raw 3D MRI, strictly unseen diseases, the drug/trial analyses and official comparator checkpoints require separate data or implementations. '
        '<a href="https://www.medrxiv.org/content/10.64898/2026.09.19.26363460v3.full">HealthFlux v3</a></p>']
    page='<meta charset="utf-8"><style>body{font:15px/1.55 system-ui;max-width:1120px;margin:36px auto;color:#19343f}table{border-collapse:collapse}td,th{padding:4px 10px;border-bottom:1px solid #ccd6d9}img{max-width:100%}h2{margin-top:32px}</style>'
    (web/'index.html').write_text(page+'\n'.join(parts))
    log(f'Paper report: {web}/index.html')
