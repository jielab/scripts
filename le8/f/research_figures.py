"""A few aggregate figures, not one file per assay. Dynamic details live in Shiny."""
from pathlib import Path
import numpy as np
import pandas as pd

def _save(fig,path):
    fig.tight_layout();fig.savefig(path,dpi=160,bbox_inches='tight')
    import matplotlib.pyplot as plt
    plt.close(fig)

def plot_deviation(tables,out,anchors):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    status=tables.get('status',pd.DataFrame())
    tasks=status.task.dropna().unique() if 'task' in status else []
    labels={'objective_plaque':'Objective coronary plaque','coronary_calcium':'Coronary calcium',
            'known_cad_state':'Known CAD state (exploratory)'}
    endpoint=labels.get(tasks[0],'Recorded outcome state') if len(tasks)==1 else 'Recorded outcome state'
    t=tables.get('reference_calibration',pd.DataFrame())
    if len(t):
        v=pd.to_numeric(t.own_PGS_partial_R2,errors='coerce').dropna()
        fig,ax=plt.subplots(figsize=(9,5));ax.hist(v,bins=35);ax.axvline(0,ls='--')
        ax.set(xlabel='Held-out partial R² from own molecular PGS',ylabel='Assays',title='Genetic calibration: all eligible assays')
        _save(fig,out/'Fig1_genetic_calibration.png')
    t=tables.get('model_metrics',pd.DataFrame())
    if len(t) and 'AUC' in t:
        fig,ax=plt.subplots(figsize=(10,5));ax.plot(t.AUC,t.model,'o')
        ax.set(xlabel='Held-out AUC (same endpoint population)',title=endpoint+'; frozen research models')
        _save(fig,out/'Fig2_diagnostic_comparison.png')
    t=tables.get('roc_points',pd.DataFrame())
    if len(t):
        fig,ax=plt.subplots(figsize=(7,6))
        for name,d in t.groupby('model',sort=False):ax.plot(d.FPR,d.TPR,label=name)
        ax.plot([0,1],[0,1],ls='--');ax.set(xlabel='False-positive rate',ylabel='True-positive rate',title=endpoint+'; held-out ROC')
        ax.legend(fontsize=8);_save(fig,out/'Fig3_ROC.png')
    t=tables.get('paired_contrasts',pd.DataFrame())
    if len(t):
        d=t[t.metric=='AUC'];fig,ax=plt.subplots(figsize=(10,5))
        xx=np.arange(len(d));ax.hlines(xx,d.lower,d.upper);ax.plot(d.delta,xx,'o')
        ax.set_yticks(xx,d.contrast);ax.axvline(0,ls='--');ax.set(xlabel='Paired ΔAUC; family bootstrap, frozen fits',title='Incremental and reparameterization comparisons')
        _save(fig,out/'Fig4_incremental_value.png')
    t=tables.get('feature_state_summary',pd.DataFrame())
    if len(t):
        d=t[t.feature.isin(anchors)] if anchors else t[t.feature.isin(t.feature.unique()[:3])]
        fig,ax=plt.subplots(figsize=(10,5))
        for state,s in d.groupby('state'):
            ax.plot(s.feature,s.mean_deviation,'o',label=f'{endpoint}={state}')
        ax.axhline(0,ls='--');ax.set(ylabel='Mean genetic deviation (reference residual SD)',title='Declared illustrative assays; all assays available in browser')
        ax.legend();_save(fig,out/'Fig5_anchor_deviations.png')
    t=tables.get('decision_curve',pd.DataFrame())
    if len(t):
        fig,ax=plt.subplots(figsize=(8,5))
        for name,d in t.groupby('model',sort=False):ax.plot(d.threshold,d.net_benefit,label=name)
        d=t[t.model==t.model.iloc[0]];ax.plot(d.threshold,d.treat_all,ls='--',label='Treat all');ax.axhline(0,ls=':')
        ax.set(xlabel='Risk threshold',ylabel='Net benefit',title='Exploratory decision curve in the observed ascertainment sample')
        ax.legend(fontsize=7);_save(fig,out/'Fig6_decision_curve.png')
    if not len(tables.get('model_metrics',[])):
        fig,ax=plt.subplots(figsize=(9,4));ax.axis('off')
        ax.text(.03,.6,endpoint+': validation not completed',fontsize=16)
        ax.text(.03,.3,'Reference calibration is available. Missing endpoints are not negative results.',wrap=True)
        _save(fig,out/'Fig_status_not_validated.png')

def plot_abm(tables,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    t=tables.get('metrics',pd.DataFrame())
    if len(t):
        fig,ax=plt.subplots(figsize=(8,5))
        d=t[t.subset.str.startswith('support_quality:',na=False)]
        for name,s in d.groupby('model',sort=False):
            if name in ('elasticnet_full','elasticnet_clean','twin_clean'):
                s=s.sort_values('coverage');ax.plot(s.coverage,s.AUC_IPCW,'o-',label=name)
        ax.set(xlabel='Actual held-out population coverage',ylabel='IPCW AUC',title='Frozen twin-support thresholds: same-mask comparisons')
        if len(d):ax.legend(fontsize=8)
        _save(fig,out/'Fig1_twin_coverage_AUC.png')
    t=tables.get('confidence_validation',pd.DataFrame())
    if len(t):
        fig,ax=plt.subplots(figsize=(8,5))
        for name,d in t.groupby('model',sort=False):ax.plot(d.mean_confidence,d.Brier_IPCW,'o-',label=name)
        ax.set(xlabel='Pre-outcome confidence score (0–100)',ylabel='Observed IPCW Brier error',title='Empirical validation of twin confidence, not individual certainty');ax.legend()
        _save(fig,out/'Fig2_confidence_validation.png')
    t=tables.get('contrasts',pd.DataFrame())
    if len(t):
        d=t[(t.metric=='AUC_IPCW')&t.contrast.isin(['selection_effect','training_effect_same_subset','combined_vs_bulk','beyond_confidence_selection'])]
        fig,ax=plt.subplots(figsize=(10,5));y=np.arange(len(d))
        # Percentile intervals need not contain the point in finite samples.
        ax.hlines(y,d.lower,d.upper);ax.plot(d.delta,y,'o');ax.set_yticks(y,d.contrast);ax.axvline(0,ls='--')
        ax.set(xlabel='Paired ΔAUC',title='Selection versus training-filter effect')
        _save(fig,out/'Fig3_effect_decomposition.png')
