"""Small aggregate figure set. No per-person/per-SNP static plot expansion."""
from pathlib import Path
import numpy as np

def draw_report(tables,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    m=tables['metrics'];c=tables['confidence_validation'];d=tables['contrasts']
    score='Brier' if 'Brier' in m else 'MSE'
    # One chart per file, with no fixed palette.
    allrows=m[(m.ancestry=='ALL')&m['subset'].isin(['all','selected','rejected'])]
    fig,ax=plt.subplots(figsize=(10,6))
    for name,z in allrows.groupby('method',sort=False):
        z=z.set_index('subset').reindex(['all','selected','rejected']);ax.plot(z.index,z[score],marker='o',label=name)
    ax.set(xlabel='Identical evaluation subsets',ylabel=score+' (lower is better)',title='Full / selected / rejected — same-person model comparison')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'same_subset.png',dpi=170);plt.close(fig)
    if len(c):
        fig,ax=plt.subplots(figsize=(9,5))
        for name,z in c[c.ancestry=='ALL'].groupby('method',sort=False):
            z=z.sort_values('mean_confidence');ax.plot(z.mean_confidence,z[score],marker='o',label=name)
        ax.set(xlabel='Pre-outcome confidence (reference support, 0–100)',ylabel=score,title='Confidence versus held-out error')
        ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'confidence_error.png',dpi=170);plt.close(fig)
    selected=d[(d.ancestry=='ALL')&(d.metric==score)&d.contrast.isin(['selection_effect','training_effect_same_subset','combined_effect','beyond_baseline_confidence','beyond_risk_ancestry_composition'])]
    fig,ax=plt.subplots(figsize=(10,5))
    for i,row in enumerate(selected.itertuples()):
        if np.isfinite(row.lower) and np.isfinite(row.upper) and row.lower<=row.delta<=row.upper:
            ax.errorbar(row.delta,i,xerr=[[row.delta-row.lower],[row.upper-row.delta]],fmt='o',capsize=3)
        else:
            ax.plot(row.delta,i,marker='o')
            if np.isfinite(row.lower) and np.isfinite(row.upper):ax.plot([row.lower,row.upper],[i,i])
    ax.axvline(0,linestyle='--');ax.set_yticks(np.arange(len(selected)),selected.contrast)
    ax.set(xlabel='Paired Δ'+score+' (negative favours the candidate)',title='Selection and training are different effects')
    fig.tight_layout();fig.savefig(out/'effect_decomposition.png',dpi=170);plt.close(fig)
    q=m[(m.ancestry=='ALL')&m['subset'].str.startswith('support_quality:')]
    fig,ax=plt.subplots(figsize=(9,5))
    for name,z in q.groupby('method',sort=False):
        z=z.sort_values('actual_coverage');ax.plot(z.actual_coverage,z[score],marker='o',label=name)
    ax.set(xlabel='Actual coverage of the original test roster',ylabel=score,title='Frozen twin-confidence coverage')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'coverage.png',dpi=170);plt.close(fig)
