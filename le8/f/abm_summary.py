"""Current ABM aggregate sources and descriptive summaries for Final."""
import os
from pathlib import Path
import pandas as pd

BOUNDARY = ('Confidence describes reference support, not the probability that a prediction is correct. '
            'Selection and training effects require their separate paired contrasts. '
            'Missing family groups remain unverified; missing coronary endpoints do not establish diagnostic accuracy.')


def source_root(root):
    # Native reports run in a temporary table workspace; ABM owns its publication.
    return Path(os.getenv('LE8_PUBLISHED_ROOT', str(root)))


def table_paths(root, trait, layer):
    base = source_root(root) / trait / layer / 'c1_abm'
    twin = base / 'abm_twin'
    deviation = base / 'molecular_deviation'
    return {
        'abm_metrics': twin / 'test_report/metrics.tsv.gz',
        'abm_coverage': twin / 'test_report/metrics.tsv.gz',
        'abm_paired': twin / 'test_report/contrasts.tsv.gz',
        'abm_confidence': twin / 'test_report/confidence_validation.tsv.gz',
        'abm_support': twin / 'aggregate/support_status.tsv.gz',
        'abm_training': twin / 'fit/aggregate/training_summary.tsv.gz',
        'deviation_metrics': deviation / 'results/model_metrics.tsv.gz',
        'deviation_status': deviation / 'results/status.tsv.gz',
        'deviation_support': deviation / 'aggregate/support_status.tsv.gz',
    }


def summarize(metrics):
    full = metrics.loc[metrics['subset'].eq('all')] if 'subset' in metrics else pd.DataFrame()
    fields = [c for c in ['AUC_IPCW', 'Uno_C_horizon', 'Brier_IPCW'] if c in full]
    rows = []
    for _, row in full.iterrows():
        values = ', '.join(f'{c}={row[c]:.4f}' for c in fields if pd.notna(row[c]))
        rows.append(str(row.get('model', 'model')) + ': ' + values)
    return dict(status='exploratory' if rows else 'unavailable',
                key_result='; '.join(rows) or 'Current twin model results unavailable.',
                full_cohort_models=len(full))
