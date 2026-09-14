"""Descriptive molecular-state interpretation and auditable local research report."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sksurv.nonparametric import kaplan_meier_estimator


def report(x,p,z,g,features,a,out,root):
    tr=p.split.values=='train';te=p.split.values=='test';states=g['state']
    profiles=pd.DataFrame([x[tr & (states==s)].mean(axis=0) for s in range(g['soft'].shape[1])],columns=features,
                          index=[f'S{s+1}' for s in range(g['soft'].shape[1])])
    profiles.index.name='state';profiles.to_csv(out/'training_state_profiles.csv')
    top=[]
    for state,row in profiles.iterrows():
        for feature in row.abs().nlargest(20).index:top.append(dict(state=state,feature=feature,mean_standardized=float(row[feature])))
    pd.DataFrame(top).to_csv(out/'state_top_features.csv',index=False)
    zx=(x[tr]-x[tr].mean(axis=0))/np.maximum(x[tr].std(axis=0),1e-8)
    zz=(z[tr]-z[tr].mean(axis=0))/np.maximum(z[tr].std(axis=0),1e-8)
    correlations=zx.T@zz/len(zz)
    pd.DataFrame(correlations,index=features,columns=[f'AE{i+1}' for i in range(z.shape[1])]).to_csv(out/'training_axis_feature_correlations.csv')
    persons=p[[a.id_col,'split','time','event']].copy();persons['state']=states+1
    persons['novel']=g['novelty'];persons['neighbor_distance']=g['neighbor_distance']
    for j in range(g['soft'].shape[1]):persons[f'state_weight_{j+1}']=g['soft'][:,j]
    for j in range(z.shape[1]):persons[f'AE{j+1}']=z[:,j]
    for j in range(g['diffusion'].shape[1]):persons[f'diffusion{j+1}']=g['diffusion'][:,j]
    persons.to_csv(out/'person_molecular_states.csv',index=False)
    counts=persons.groupby(['split','state']).agg(n=('event','size'),events=('event','sum'),novel=('novel','mean'))
    counts.to_csv(out/'state_counts.csv')
    fig,ax=plt.subplots(figsize=(7,5));rng=np.random.default_rng(a.seed);ix=np.flatnonzero(tr)
    ix=rng.choice(ix,min(8000,len(ix)),replace=False)
    ax.scatter(z[ix,0],z[ix,1],c=states[ix],s=4,cmap='tab20',alpha=.5,rasterized=True)
    ax.set(xlabel='AE coordinate 1',ylabel='AE coordinate 2',title='Discovery molecular landscape (two coordinates)')
    fig.tight_layout();fig.savefig(out/'landscape.png',dpi=180);plt.close(fig)
    cols=list(dict.fromkeys(profiles.loc[s].abs().nlargest(5).index.tolist()[j] for s in profiles.index for j in range(min(5,len(features)))))
    fig,ax=plt.subplots(figsize=(max(8,len(cols)*.22),4))
    im=ax.imshow(profiles[cols],aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
    ax.set_xticks(range(len(cols)),cols,rotation=90,fontsize=7);ax.set_yticks(range(len(profiles)),profiles.index)
    ax.set_title('Training state profiles; descriptive standardized means');fig.colorbar(im,ax=ax)
    fig.tight_layout();fig.savefig(out/'state_heatmap.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,5))
    for state in sorted(np.unique(states)):
        ix=te&(states==state)
        if not ix.any():continue
        t,s=kaplan_meier_estimator(p.event.values[ix].astype(bool),p.time.values[ix])
        ax.step(np.r_[0,t],np.r_[1,s],where='post',label=f'S{state+1} (n={ix.sum()})')
    ax.set(xlabel='Years since baseline',ylabel='CAD-free net survival',title='Independent test set; death treated as censoring')
    ax.legend();fig.tight_layout();fig.savefig(out/'test_state_survival.png',dpi=180);plt.close(fig)
    metrics=pd.read_csv(root/'s5_predict/test_metrics.csv');m=metrics[metrics.metric=='Harrell_C'].set_index('model').value
    fig,ax=plt.subplots(figsize=(8,5));ax.barh(m.index,m.values,color='#347e91');ax.set(xlim=(.45,1),xlabel='Held-out Harrell C')
    fig.tight_layout();fig.savefig(out/'test_model_comparison.png',dpi=180);plt.close(fig)
    audit=json.loads((root/'s2_preprocess/cohort_audit.json').read_text())
    gs=json.loads((root/'s4_graph/summary.json').read_text())
    text=f'''# Panome: {a.biom} → molecular state → {a.trait}

Run type: {'SYNTHETIC SOFTWARE TEST' if a.demo else ('REAL-DATA PILOT, not final inference' if a.max_samples else 'REAL DATA')}.

Input cohort: {audit['joined']}; eligible after phenotype rules: {audit['eligible']}; after omics row QC: {len(p)}.
Retained molecular features: {len(features)}. Train / validation / test: {dict(p.split.value_counts())}.
Baseline exclusion date columns: {a.diagnosis_col}, {a.healthy_date_cols or '(no additional endpoints)'}.
This is CAD-free baseline sampling unless a broader disease panel was explicitly supplied.

Representation: denoising AE {len(features)} → {a.hidden} → 64 → {a.latent}, with observed-entry loss and validation early stopping.
Residualization: {a.residualize or 'none'}. Categorical variables: {a.categorical}.
Graph uses training AE coordinates only; held-out people use frozen neighbor projection.
States: {gs['n_states']}; resolution: {gs['resolution']}; connected components: {gs['n_components']}.
Mean ARI after 10% edge removal: {gs['perturbation_ARI_mean']:.3f}. This measures graph perturbation stability, not full-pipeline reproducibility.

![Landscape](landscape.png)
![Profiles](state_heatmap.png)
![Test survival](test_state_survival.png)
![Comparison](test_model_comparison.png)

## Held-out model comparison

{m.to_string()}

See ../s5_predict/test_metrics.csv for IPCW AUC, Brier score and Uno C at configured horizons;
metric_limitations.json records unsupported evaluations. test_calibration.csv contains net-risk calibration by quintile.
Test bootstrap confidence intervals are conditional on this split and trained models, not repeated nested-CV intervals.

## Interpretation and limits

State weights are local neighbor membership weights, not measured biological pathway percentages.
State profiles and axis-feature correlations are descriptive; cluster-defining features do not provide independent discovery P values.
Marginal PWAS is a training-only univariate Cox screen on the same corrected molecular matrix, followed by a clinically adjusted score model;
it is not a multivariable-adjusted etiological PWAS. PCA and elastic-net Cox provide additional reference models.
Test-set state Cox estimates and PH diagnostics are exploratory; sparse states or PH violations must be checked.
The {a.lag_years:g}-year landmark analysis keeps only participants still under follow-up and CAD-free at that landmark.
Death is censored: predicted 1-S is net risk, not competing-risk cumulative incidence.
Neither baseline AE coordinates nor clustering establishes disease trajectory, causal mechanism, or molecular ancestry.
A random person split does not ensure family separation unless --group-col was supplied.
Inputs may already have undergone upstream normalization/imputation; original raw-data missingness and upstream leakage cannot be recovered.
External validation, repeated split / AE-seed stability, competing-risk modeling, and pathway validation remain necessary for publication.

Sources: [original discussion](https://chatgpt.com/c/6aa01c8e-2ce0-83ec-981e-fdca59accfbd),
[Zhang et al., Genome Medicine 2026](https://doi.org/10.1186/s13073-026-01696-w).
This pipeline is an AE-based extension of the sample-graph idea, not a reproduction of that paper.
'''
    (out/'REPORT.md').write_text(text)
