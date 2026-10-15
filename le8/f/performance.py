#!/usr/bin/env python3
"""Cumulative performance snapshots from existing models, without fitting."""
import argparse
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import textwrap
import zipfile

import numpy as np
import pandas as pd

from research_io import write_xlsx
from c1_abm_compare import comparison_from_paths

STAGES = ['c1_correlate', 'c1_abm', 'c2_cause', 'c3_coloc', 'c4_connect', 'c5_cellulation', 'final']


def read_export(path):
    """Read the exact aggregate export, including the native workbook archive."""
    if path.is_file():
        return pd.read_csv(path, sep='\t' if '.tsv' in path.name else ',')
    if not path.parent.is_dir():
        return pd.DataFrame()
    for book in sorted(path.parent.glob('*.xlsx')):
        with zipfile.ZipFile(book) as archive:
            if 'le8/manifest.json' not in archive.namelist():
                continue
            entry = json.loads(archive.read('le8/manifest.json')).get('files', {}).get(path.name)
            if not entry:
                continue
            data = archive.read(entry['part'])
            if hashlib.sha256(data).hexdigest() != entry['sha256']:
                raise ValueError('Corrupt performance export: '+str(path))
            if path.name.endswith('.gz'):
                data = gzip.decompress(data)
            return pd.read_csv(io.BytesIO(data), sep='\t' if '.tsv' in path.name else ',')
    return pd.DataFrame()


def collect(root, published, trait, layer, stage):
    base = root/trait/layer
    twin = published/trait/layer/'c1_abm/abm_twin'
    max_stage = STAGES.index(stage)
    rows, availability = [], []

    def add(module, path, family, metrics, filters=None, horizon=None):
        d = read_export(path)
        for column, value in (filters or {}).items():
            if column in d:
                d = d.loc[d[column].isin(value) if isinstance(value, list) else d[column].eq(value)]
        count = 0
        for _, r in d.iterrows():
            model = str(r.get('model', r.get('rule', 'unspecified')))
            panel = str(r.get('biom_set', ''))
            context = ', '.join(f'{k}={r[k]}' for k in ['budget', 'budget_per_layer', 'ablation', 'arm', 'landmark', 'horizon', 'subset']
                                if k in r and pd.notna(r[k]))
            for column, metric in metrics.items():
                value = pd.to_numeric(r.get(column), errors='coerce')
                if not np.isfinite(value):
                    continue
                rows.append(dict(stage=module, family=family, model=model, panel=panel, context=context,
                                 metric=metric, value=value, N=r.get('N', np.nan),
                                 horizon=r.get('horizon', r.get('AUC_horizon_years', horizon)),
                                 evaluation_scope=family, source=str(published/path.relative_to(root)) if path.is_relative_to(root) else str(path)))
                count += 1
        availability.append(dict(stage=module, source=family, status='evaluated' if count else 'not_evaluated',
                                 metric_rows=count, detail='Existing model estimates only'))

    add('c1_correlate', base/'c1_correlate/prediction/prediction_summary.csv', 'C1 conventional holdout',
        dict(AUC='IPCW AUC', C_index='Harrell C'))
    add('c1_correlate', base/'c1_correlate/prediction/review_budget_metrics.csv', 'C1 fixed-budget holdout',
        dict(AUC='IPCW AUC', Brier='IPCW Brier'))
    add('c1_correlate', base/'c1_correlate/prediction/metrics_summary.csv', 'C1 optional nested cross-validation',
        dict(AUC10_mean='Mean CV IPCW AUC', Brier10_mean='Mean CV IPCW Brier', cindex_mean='Mean CV Harrell C'))
    # ABM is a C1 analysis. Include an already completed ABM even when refreshing C1.
    add('c1_abm', twin/'test_report/metrics.tsv.gz', 'ABM test (technical QC)',
        dict(AUC_IPCW='IPCW AUC', Uno_C_horizon='Uno C', Brier_IPCW='IPCW Brier', LogLoss_IPCW='IPCW LogLoss'),
        filters=dict(subset=['all', 'selected', 'rejected']))
    add('c1_abm', published/trait/layer/'c1_abm/molecular_deviation/results/model_metrics.tsv.gz',
        'Molecular deviation objective endpoint', dict(AUC='Binary AUC', AUPRC='AUPRC', Brier='Binary Brier', LogLoss='Binary LogLoss'))
    for module in ['c2_cause', 'c3_coloc', 'c5_cellulation']:
        if STAGES.index(module) <= max_stage:
            availability.append(dict(stage=module, source=module, status='not_evaluated', metric_rows=0,
                                     detail='This module currently produces evidence/annotation, not a prediction model'))
    if STAGES.index('c4_connect') <= max_stage:
        add('c4_connect', base/'c4_connect/c4.focus.metrics.csv', 'C4 panel validation',
            dict(AUC='IPCW AUC', Brier='IPCW Brier'), filters=dict(stratum='All', landmark=0), horizon=10)
        add('c4_connect', base/'c4_connect/c4.penalty_summary.csv', 'C4 penalty cross-validation',
            dict(mean_cindex='Mean CV Harrell C'))
    if stage == 'final':
        add('final', root/'final'/trait/'joint/metrics.csv', 'Joint modality out-of-fold',
            dict(AUC='IPCW AUC', Brier='IPCW Brier'), filters=dict(stratum='All', landmark=0))
        add('final', root/'final'/trait/'joint/factorial_metrics.csv', 'Joint modality within-fold concordance',
            dict(estimate='Horizon-restricted Harrell C'), filters=dict(metric='C_index', landmark=0, arm='all_assays'))
    benchmark = base/'c1_correlate/prediction/_LOCAL_ONLY/benchmark.json.gz'
    comparison = comparison_from_paths(benchmark, twin)
    availability.sort(key=lambda r: STAGES.index(r['stage']))
    return dict(performance=pd.DataFrame(rows), stages=pd.DataFrame(availability),
                c1_abm_comparison=comparison.get('metrics', pd.DataFrame()),
                c1_abm_status=comparison['status'])


def plot_tables(tables, output, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    performance = tables['performance']
    panels = []
    if len(performance):
        for (family, metric), d in performance.groupby(['family', 'metric'], sort=False):
            # Preserve every model/metric; paginate long groups instead of ranking winners.
            for start in range(0, len(d), 24):
                panels.append((family, metric, d.iloc[start:start+24]))
    if not panels:
        panels = [('No evaluated prediction models yet', '', pd.DataFrame())]
    stages = tables['stages'].groupby('stage', sort=False).metric_rows.sum()
    status = '; '.join(f'{k}: '+('available' if v else 'not evaluated') for k, v in stages.items())
    for page in range(0, len(panels), 4):
        selected = panels[page:page+4]
        suffix = '' if page == 0 else f'_{page//4+1:02d}'
        stem = 'performance_comparison'+suffix
        nr = (len(selected)+1)//2
        height = max(3.5, max(len(d) for _, _, d in selected)*.42+2.3)
        fig, axes = plt.subplots(nr, 2, figsize=(18, nr*height), squeeze=False)
        for i, (ax, (family, metric, d)) in enumerate(zip(axes.ravel(), selected)):
            export = d.copy() if len(d) else tables['stages'].copy()
            export['panel'] = chr(65+i)
            export.to_csv(output/f'{stem}.panel_{chr(65+i)}.csv', index=False)
            if len(d):
                labels = d.apply(lambda r: textwrap.fill(' / '.join(filter(None, [r.panel, r.model, r.context])), 56), axis=1)
                ax.scatter(d.value, np.arange(len(d)), color='#247b83', s=30)
                ax.set_yticks(np.arange(len(d)), labels, fontsize=8)
                ax.invert_yaxis()
                ax.grid(axis='x', alpha=.2)
                if 'Harrell C' in metric or metric.endswith('AUC') or metric in {'Uno C', 'AUPRC'}:
                    ax.set_xlim(0, 1)
                ax.set_xlabel(metric)
            else:
                ax.text(.5, .5, 'Not evaluated', ha='center', transform=ax.transAxes)
            ax.set_title(family+' | '+metric, fontsize=11)
        for ax in axes.ravel()[len(selected):]:
            ax.set_visible(False)
        fig.suptitle(title, fontsize=15)
        fig.text(.02, .012, textwrap.fill(status, 155)+'\nPanels retain their own cohorts/endpoints; cross-panel differences are descriptive.', fontsize=9)
        fig.tight_layout(rect=(0, .06, 1, .96))
        fig.savefig(output/f'{stem}.png', dpi=150)
        plt.close(fig)
    paired = tables['c1_abm_comparison']
    if not paired.empty:
        for subset, d in paired.groupby('subset', sort=False):
            d = d.loc[d.status.eq('evaluated')]
            for page in range(0, len(d), 30):
                z = d.iloc[page:page+30]
                fig, ax = plt.subplots(figsize=(14, max(4, len(z)*.45+2)))
                labels = z.apply(lambda r: textwrap.fill(f'{r.model} vs {r.baseline_panel} / {r.baseline_model} (n={r.N})', 75), axis=1)
                ax.scatter(z.delta_AUC, np.arange(len(z)), color='#386cb0')
                ax.set_yticks(np.arange(len(z)), labels, fontsize=8)
                ax.invert_yaxis(); ax.axvline(0, color='grey', linestyle='--')
                ax.set_xlabel('ABM minus C1 IPCW AUC (same held-out participants and weights)')
                ax.set_title(f'{title} | {subset} | exploratory point estimates')
                fig.tight_layout()
                stem = f'C1_ABM_comparison_{subset}_{page//30+1:02d}'
                z.to_csv(output/f'{stem}.panel_A.csv', index=False)
                fig.savefig(output/f'{stem}.png', dpi=150)
                plt.close(fig)


def generate(root, published, traits, layers, stage):
    for trait in traits:
        for layer in layers:
            output = (root/'final'/trait/layer if stage == 'final' else root/trait/layer/stage)/'performance'
            tables = collect(root, published, trait, layer, stage)
            write_snapshot(tables, output, f'{trait} / {layer} | through {stage}')


def write_snapshot(tables, output, title):
    output.mkdir(parents=True, exist_ok=True)
    # Only files owned by this report are refreshed.
    for pattern in ['performance_comparison*.png', 'C1_ABM_comparison_*.png',
                    'performance_comparison*.panel_*.csv', 'C1_ABM_comparison_*.panel_*.csv']:
        for path in output.glob(pattern):
            path.unlink()
    safe = {name: d if len(d.columns) else pd.DataFrame(columns=['status']) for name, d in tables.items()}
    for name, d in safe.items():
        d.to_csv(output/(name+'.csv'), index=False)
    write_xlsx(output/'performance_summary.xlsx', safe)
    plot_tables(safe, output, title)
    print('[LE8] DONE cumulative performance: '+str(output), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-root', type=Path, required=True)
    p.add_argument('--Y', required=True)
    p.add_argument('--biom', required=True)
    p.add_argument('--stage', choices=STAGES, required=True)
    a = p.parse_args()
    generate(a.analysis_root, Path(os.getenv('LE8_PUBLISHED_ROOT', str(a.analysis_root))),
             a.Y.split(','), a.biom.split(','), a.stage)


if __name__ == '__main__':
    main()
