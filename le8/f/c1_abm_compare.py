"""Compare frozen C1 and TwinABM scores on their common held-out participants."""
import gzip
import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from abm_twin_evaluate import uno_c


def compare(benchmark, prediction, truth, km, config):
    """Use the same outcomes, censoring weights and masks for each model pair.

    C1 exports risk rankings, not calibrated probabilities. Only discrimination
    can be compared; its Harrell C is not substituted for ABM's Uno C.
    """
    def unavailable(reason):
        return {'status': pd.DataFrame([dict(status=reason)]), 'metrics': pd.DataFrame()}

    if benchmark.get('format') != 'le8-c1-benchmark-v1':
        raise ValueError('Unknown C1 benchmark format')
    if benchmark.get('trait') != config.trait or benchmark.get('layer') != config.layer:
        return unavailable('trait_or_layer_mismatch')
    if not np.isclose(float(benchmark['horizon']), config.horizon, rtol=0, atol=1e-10):
        return unavailable('prediction_horizon_mismatch')
    ids = np.asarray(prediction['eid'], str)
    if len(set(ids)) != len(ids) or not np.array_equal(ids, np.asarray(truth['eid'], str)):
        raise ValueError('Invalid ABM prediction/truth IDs')
    time, event = np.asarray(truth['time'], float), np.asarray(truth['event'], int)
    c1 = pd.DataFrame(benchmark['predictions'])
    required = {'eid', 'biom_set', 'model', 'time', 'event', 'score'}
    if c1.empty:
        return unavailable('C1_predictions_unavailable')
    if not required <= set(c1):
        raise ValueError('Incomplete C1 benchmark predictions')
    c1['eid'] = c1.eid.astype(str)
    if c1.duplicated(['eid', 'biom_set', 'model']).any():
        raise ValueError('Duplicate C1 model/participant predictions')
    joined = c1.merge(pd.DataFrame(dict(eid=ids, abm_row=np.arange(len(ids)))), on='eid')
    if joined.empty:
        return unavailable('no_common_heldout_participants')
    j = joined.abm_row.to_numpy()
    if (not np.allclose(joined.time, time[j], rtol=0, atol=1e-8)
            or not np.array_equal(joined.event.to_numpy(), event[j])):
        return unavailable('outcome_definition_mismatch')
    qc = np.asarray(prediction['diagnostics']['technical_qc'], bool)
    selected = np.asarray(prediction['selected'], bool) & qc
    masks = dict(all=qc, selected=selected, rejected=qc & ~selected)
    rows = []
    # Report primary, unfiltered and deployed fallback predictions separately.
    models = list(dict.fromkeys([config.primary, 'elasticnet_full', 'fallback_policy']))
    for (panel, baseline), d in joined.groupby(['biom_set', 'model'], sort=False):
        indices = d.abm_row.to_numpy()
        scores = d.score.to_numpy(float)
        for name in models:
            if name not in prediction['probability']:
                continue
            p = np.asarray(prediction['probability'][name], float)[indices]
            for subset, mask in masks.items():
                use = mask[indices] & np.isfinite(scores) & np.isfinite(p)
                ix, baseline_score, abm_score = indices[use], scores[use], p[use]
                y, w = km.labels_weights(time[ix], event[ix])
                ok = w > 0
                both = len(np.unique(y[ok])) == 2
                ba = float(roc_auc_score(y[ok], baseline_score[ok], sample_weight=w[ok])) if both else np.nan
                aa = float(roc_auc_score(y[ok], abm_score[ok], sample_weight=w[ok])) if both else np.nan
                bc, pairs = uno_c(time[ix], event[ix], baseline_score, km)
                ac, _ = uno_c(time[ix], event[ix], abm_score, km)
                rows.append(dict(
                    baseline_panel=panel, baseline_model=baseline, model=name, subset=subset,
                    N=len(ix), cases=int(y.sum()), controls=int(np.sum(ok & (y == 0))),
                    common_holdout_N=len(d), ABM_test_N=len(ids),
                    coverage_of_ABM_test=len(ix)/len(ids) if len(ids) else 0,
                    horizon=config.horizon, C1_AUC_IPCW=ba, ABM_AUC_IPCW=aa, delta_AUC=aa-ba,
                    C1_Uno_C_horizon=bc, ABM_Uno_C_horizon=ac, delta_Uno_C=ac-bc,
                    comparable_pairs=pairs,
                    status='evaluated' if both else 'insufficient_case_control_information'))
    metrics = pd.DataFrame(rows)
    status = dict(
        status='evaluated' if len(metrics) and metrics.status.eq('evaluated').any() else 'insufficient_case_control_information',
        scope='intersection of C1 outer test and ABM test, same subset for both models',
        common_holdout_N=joined.eid.nunique(), ABM_test_N=len(ids), C1_test_N=c1.eid.nunique(),
        censoring='same frozen ABM development censoring KM for both models',
        grouping_unit=config.grouping_unit,
        inference='Exploratory frozen-fit point comparisons; no confidence interval or superiority claim',
        calibration='C1 exports rankings; Brier and LogLoss are not compared')
    return dict(metrics=metrics, status=pd.DataFrame([status]))


def run_comparison(args, scope, work, publish, completed, finish):
    """An independent publication: refreshing a comparison never refits ABM."""
    from abm_twin_core import json_write
    from research_io import write_tables, sha256

    benchmark = args.analysis_root/args.Y/args.biom/'c1_correlate/prediction/_LOCAL_ONLY/benchmark.json.gz'
    twin = scope/'abm_twin'
    inputs = dict(C1=benchmark, model=twin/'fit/MODEL_FROZEN.json',
                  predictions=twin/'test_predictions_LOCAL_ONLY.joblib', truth=twin/'test_truth_LOCAL_ONLY.npz')
    signature = {key: sha256(path) if path.is_file() else 'unavailable' for key, path in inputs.items()}
    target = scope/'baseline_comparison'
    if completed(target) and (target/'sources.json').is_file():
        if json.loads((target/'sources.json').read_text()) == signature:
            return
    stage = work/'baseline_comparison'
    tables = comparison_from_paths(benchmark, twin)
    write_tables(stage, tables, 'ABM twin', args.biom, args.Y, workbook='C1_ABM_comparison.xlsx')
    json_write(stage/'sources.json', signature)
    finish(stage, dict(source='C1 basic prediction versus frozen ABM', source_hashes=signature))
    publish(stage, target)
    print('[LE8 twin] DONE C1 baseline comparison: '+str(tables['status'].iloc[0]['status']), flush=True)


def comparison_from_paths(benchmark, twin):
    """Read trusted local frozen models; never fit, tune or change a test split."""
    from abm_twin_core import model_hash
    if not benchmark.is_file():
        return dict(status=pd.DataFrame([dict(status='C1_benchmark_unavailable')]))
    required = ['fit/model_LOCAL_ONLY.joblib', 'test_predictions_LOCAL_ONLY.joblib', 'test_truth_LOCAL_ONLY.npz']
    if not all((twin/name).is_file() for name in required):
        return dict(status=pd.DataFrame([dict(status='ABM_frozen_predictions_unavailable')]))
    import joblib
    with gzip.open(benchmark, 'rt') as stream:
        baseline = json.load(stream)
    model = joblib.load(twin/required[0])
    prediction = joblib.load(twin/required[1])
    if prediction['frozen_model_hash'] != model_hash(model):
        raise ValueError('ABM comparison predictions do not match frozen model')
    with np.load(twin/required[2], allow_pickle=False) as truth:
        return compare(baseline, prediction, truth, model.km, model.config)
