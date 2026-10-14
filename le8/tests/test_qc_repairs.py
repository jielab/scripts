"""Synthetic acceptance checks for the 2026-10-14 minimal repair list."""
from dataclasses import replace
import inspect
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
from abm_twin_core import Config, NumericPrep, ROLES, TwinABM, TwinBank, model_hash
from abm_twin_evaluate import evaluate
from molecular_deviation import (DeviationConfig, DeviationStudy, evaluate_deviation,
                                 feature_state_associations, objective_endpoint)


@pytest.fixture(scope='module')
def abm():
    rng = np.random.default_rng(817)
    n = 1300
    M = rng.normal(size=(n, 10))
    C = np.eye(50)[rng.integers(0, 50, n)]
    X = np.c_[M, C]
    event = rng.binomial(1, expit(1.5*M[:, 0] - M[:, 1]))
    time = np.where(event, 5., 15.)
    roles = np.repeat(ROLES, [700, 150, 150, 150, 150])
    bad = np.flatnonzero((roles == 'build') & (event == 0))[:30]
    X[bad, :10] = np.nan
    ids = np.array([f'dev_{i}' for i in range(n)])
    names = [f'm_{i}' for i in range(10)] + [f'c_{i}' for i in range(50)]
    config = Config(folds=3, repeats=2, min_events=3, c_grid=(.5,), pilot_c=.5,
                    k=12, min_twins=3, min_family_ess=2, gate_iterations=5,
                    audit_min_cases=3, audit_min_controls=3, bootstrap=2,
                    max_iter=8000, always_include_indices=tuple(range(10, 60)),
                    grouping_unit='individual_unverified')
    with threadpool_limits(1):
        model = TwinABM(config).fit(X, time, event, ids, ids, roles, names)
    query = np.c_[rng.normal(size=(80, 10)), np.eye(50)[rng.integers(0, 50, 80)]]
    query[:20, :10] = np.nan
    qi = np.array([f'query_{i}' for i in range(80)])
    with threadpool_limits(1):
        pred = model.predict(query, qi, qi, names)
    return model, X, time, event, ids, roles, bad, query, qi, pred


def test_twenty_missing_queries_abstain_and_stay_in_denominator(abm):
    model, _, _, _, _, _, _, _, qi, pred = abm
    d = pred['diagnostics']
    assert len(pred['eid']) == 80
    np.testing.assert_array_equal(d['molecular_missing_fraction'][:20], 1)
    assert not d['technical_qc'][:20].any()
    assert not d['supported'][:20].any()
    assert not d['matched_count'][:20].any()
    assert not d['pre_prediction_confidence'][:20].any()
    assert not pred['selected'][:20].any()
    assert np.isnan(pred['probability']['fallback_policy'][:20]).all()
    assert set(pred['reason'][:20]) == {'all_panel_molecules_missing'}
    truth_event = np.tile([0, 1], 40)
    result = evaluate(model, pred, np.where(truth_event, 5., 15.), truth_event, qi, bootstrap=2)
    assert result['status']['original_N'] == 80
    assert result['status']['technical_QC_excluded'] == 20
    assert result['status']['actual_coverage'] == pred['selected'].sum()/80
    assert result['grouping_unit'] == 'individual_unverified'


def test_thirty_missing_donors_excluded_from_all_arms_and_teacher_training(abm):
    model, _, _, _, ids, _, bad, *_ = abm
    donor = model.donor_table
    for field in ('retained', 'random_retained', 'donor_eligible', 'technical_qc'):
        assert not donor[field][bad].any()
    assert not donor['soft_weight'][bad].any()
    assert donor['technical_QC_excluded'][bad].all()
    assert not donor['profile_discordance_excluded'][bad].any()
    assert donor['horizon_known'][bad].all()
    for bank in model.banks.values():
        assert not set(bank.ids) & set(ids[bad])
    assert model.training_audit['technical_QC_excluded'] == 30
    folds = model.training_audit['oof']
    assert sum(f['train_technical_QC_excluded'] for f in folds) == 30*2*2
    assert sum(f['recipient_technical_QC_excluded'] for f in folds) == 30*2
    assert all(f['family_disjoint'] and f['selection_in_teacher_only'] for f in folds)
    assert donor['retained'].sum() == donor['random_retained'].sum()


def test_dummy_columns_cannot_dilute_molecular_missingness():
    rng = np.random.default_rng(4)
    X = rng.normal(size=(100, 10))
    query = X[:3].copy(); query[0] = np.nan; query[1, :3] = np.nan; query[2, :2] = np.nan
    base = NumericPrep().fit(X, Config())
    for count in (2, 50, 150):
        cfg = replace(Config(), always_include_indices=tuple(range(10, 10+count)))
        prep = NumericPrep().fit(np.c_[X, rng.normal(size=(100, count))], cfg)
        q = np.c_[query, np.zeros((3, count))]
        np.testing.assert_array_equal(prep.missing(q), base.missing(query))
        np.testing.assert_array_equal(prep.technical_qc(q), [False, False, True])


def test_actual_panel_qc_and_no_molecular_column_failure():
    X = np.random.default_rng(9).normal(size=(100, 12))
    cfg = replace(Config(), always_include_indices=(10, 11))
    prep = NumericPrep().fit(X, cfg)
    prep.keep = np.array([0, 1, 10, 11])
    query = X[:1].copy(); query[0, 0] = np.nan
    assert prep.missing(query, panel=False)[0] == .1
    assert prep.missing(query)[0] == .5
    assert not prep.technical_qc(query)[0]
    with pytest.raises(ValueError, match='molecular'):
        NumericPrep().fit(X, replace(cfg, always_include_indices=tuple(range(12))))


def test_frozen_projection_independent_of_batch_and_no_outcome_api(abm):
    model, *_, query, qi, pred = abm
    before = model_hash(model)
    for i in (0, 20, 31):
        with threadpool_limits(1):
            one = model.predict(query[i:i+1], qi[i:i+1], qi[i:i+1], model.feature_names)
        for name, values in pred['probability'].items():
            np.testing.assert_allclose(one['probability'][name], values[i:i+1], atol=1e-12)
        for name, values in pred['diagnostics'].items():
            np.testing.assert_allclose(one['diagnostics'][name], values[i:i+1], atol=1e-12)
        assert one['selected'][0] == pred['selected'][i]
    assert before == model_hash(model)
    assert not {'y', 'time', 'event'} & set(inspect.signature(model.predict).parameters)
    with pytest.raises(TypeError):
        model.predict(query, qi, qi, model.feature_names, y=np.zeros(len(qi)))


def test_family_exclusion_overlap_guard_and_technical_fallback(abm):
    model, *_, query, qi, pred = abm
    with pytest.raises(ValueError, match='overlaps'):
        model.predict(query[20:21], qi[20:21], model.development_groups[:1], model.feature_names)
    # A technical pass with no twins still gets the frozen full-model fallback.
    distant = query[20:21].copy(); distant[:, :10] = 1e4
    p = model.predict(distant, qi[:1], qi[:1], model.feature_names)
    assert p['diagnostics']['technical_qc'][0]
    assert not p['diagnostics']['supported'][0]
    assert p['reason'][0] == 'no_adequate_twins'
    np.testing.assert_array_equal(p['probability']['fallback_policy'], p['probability']['elasticnet_full'])
    bank = model.banks['full']
    _, diag = bank.predict(model.prep.transform(query[20:21]), qi[:1], np.array([bank.groups[0]]), np.array([.5]))
    # Direct all-family query: no donor may be borrowed from the same family.
    tiny = TwinBank().fit(np.arange(8.).reshape(4, 2), np.array([0, 1, 0, 1]), np.ones(4),
                         np.array(['a', 'b', 'c', 'd']), np.array(['F']*4), np.zeros(4),
                         replace(model.config, dimensions=2))
    tiny.radius = 100.
    _, diag = tiny.predict(np.array([[0., 1.]]), np.array(['q']), np.array(['F']), np.array([.5]))
    assert diag['matched_count'][0] == 0


@pytest.fixture
def association_data():
    rng = np.random.default_rng(123)
    n = 160
    sex = rng.integers(0, 2, n); age = rng.normal(size=n)
    y = rng.integers(0, 2, n)
    z = .7*y[:, None] + .4*sex[:, None] + age[:, None] + rng.normal(size=(n, 3))
    pred = dict(covariate_values=np.c_[1-sex, sex, age], covariate_names=['sex0', 'sex1', 'age'],
                Z=z, features=['A', 'B', 'C'], group=np.array([f'F{i//2}' for i in range(n)]), task='known_cad_state')
    return pred, y, np.ones(n, bool)


def test_full_onehot_matches_reference_encoding_with_cluster_se(association_data):
    pred, y, eligible = association_data
    full = feature_state_associations(pred, y, eligible)
    reference = feature_state_associations({**pred, 'covariate_values':pred['covariate_values'][:, 1:],
                                           'covariate_names':['sex1', 'age']}, y, eligible)
    assert full.status.eq('estimated_cluster_robust').all()
    np.testing.assert_allclose(full[['beta', 'se', 'p', 'lo', 'hi', 'FDR']],
                               reference[['beta', 'se', 'p', 'lo', 'hi', 'FDR']], rtol=1e-10, atol=1e-12)
    assert full.covariates_before.eq(3).all() and full.covariates_after.eq(2).all()
    assert full.design_rank.eq(4).all() and full.residual_df.eq(156).all()
    assert all(len(json.loads(s)) == 1 for s in full.dropped_covariates)
    assert not full.effect_measure.str.contains('objective|客观病变', regex=True).any()


def test_state_collinearity_is_not_silently_dropped(association_data):
    pred, y, eligible = association_data
    pred.update(covariate_values=np.c_[y, 1-y], covariate_names=['state', 'other'])
    result = feature_state_associations(pred, y, eligible)
    assert result.status.eq('disease_state_not_estimable').all()
    assert result[['beta', 'se', 'p', 'FDR']].isna().all().all()


def test_insufficient_groups_retains_na_and_legacy_task_is_neutral(association_data):
    pred, y, eligible = association_data
    pred['group'] = np.array(['family']*len(y)); pred.pop('task')
    result = feature_state_associations(pred, y, eligible)
    assert result.status.eq('insufficient_independent_groups').all()
    assert result.beta.isna().all()
    assert result.task_label.eq('记录状态/结局状态').all()


@pytest.mark.parametrize('task', ['objective_plaque', 'coronary_calcium', 'known_cad_state'])
def test_task_in_all_three_feature_tables_without_refit(association_data, tmp_path, task):
    from types import SimpleNamespace
    pred, y, eligible = association_data
    n, p = pred['Z'].shape
    cfg = DeviationConfig(task=task, bootstrap=1)
    model = SimpleNamespace(config=cfg, classifiers={}, reference_table=[], features=pred['features'],
                            keep=list(range(p)), classification_status='reference_only', operating_thresholds={})
    pred.update(eid=np.array([f'q{i}' for i in range(n)]), probability={}, qc=eligible,
                measured=pred['Z'], expected=np.zeros((n, p)), frozen_hash=model_hash(model))
    truth = dict(eid=pred['eid'], y=y, endpoint_eligible=eligible)
    tables, _ = evaluate_deviation(model, pred, truth, tmp_path)
    for name in ('feature_state_summary', 'feature_distribution', 'feature_state_association'):
        assert len(tables[name]) and tables[name].task.eq(task).all()
        assert 'task_label' in tables[name]
    assert pred['frozen_hash'] == model_hash(model)
    assert tables['roc_points'].empty  # No diagnostic model means no invented ROC.


def test_no_objective_file_keeps_reference_calibration_available():
    rng = np.random.default_rng(77); n = 120
    C = rng.normal(size=(n, 2)); G = rng.normal(size=(n, 3)); M = .5*G + rng.normal(size=(n, 3))
    ids = np.array([str(i) for i in range(n)])
    eligible, y, reason = objective_endpoint(pd.DataFrame({'known_cad_at_index':np.zeros(n)}), DeviationConfig())
    assert not eligible.any() and np.isnan(y).all()
    assert set(reason) == {'objective_endpoint_unavailable'}
    model = DeviationStudy(DeviationConfig(reference_min_n=20, folds=3)).fit_reference(
        M, G, C, ids, ids, np.ones(n, bool), ['A', 'B', 'C'], ['age', 'sex'])
    assert len(model.reference_table) == 3 and model.classifiers == {}


def test_explicit_outer_roster_preserves_test_people(tmp_path):
    from types import SimpleNamespace
    from cohort_inputs import assign_roles, family_groups
    p = pd.DataFrame({'eid':[str(i) for i in range(1000)]})
    outer = p.assign(role=np.where(np.arange(1000)%3 == 0, 'test', 'training'))
    path = tmp_path/'roster.tsv'; outer.to_csv(path, sep='\t', index=False)
    args = SimpleNamespace(group_col=None, id_col='eid', outer_roster=str(path), seed=2026, r_bin='Rscript')
    np.testing.assert_array_equal(family_groups(p, args), p.eid.to_numpy())
    roles = assign_roles(p, args)
    np.testing.assert_array_equal(roles == 'test', outer.role.to_numpy() == 'test')


@pytest.mark.parametrize('task', ['known_cad_state', 'coronary_calcium'])
def test_integrated_task_directories_and_unverified_family_status(tmp_path, monkeypatch, task):
    import run_host
    config = tmp_path/'abm.json'; config.write_text('{}')
    deviation_config = tmp_path/'deviation.json'; deviation_config.write_text(json.dumps({'task':task}))
    args = run_host.parser().parse_args(['--Y','cvd_cad','--biom','met','--analysis-root',str(tmp_path/'results'),
        '--config',str(config),'--deviation-config',str(deviation_config)])
    scope = args.analysis_root/'cvd_cad/met/c1_abm'
    old = scope/'molecular_deviation'; old.mkdir(parents=True)
    (old/'unchanged.txt').write_text('old objective result')
    audit = dict(independent_family_support_verified=False,family_status='unavailable_individual_ID_groups_only')
    monkeypatch.setattr(run_host,'prepare_host',lambda *args: ({},{},audit))
    monkeypatch.setattr(run_host,'completed',lambda path: path.name=='abm_twin')
    def fake_deviation(args, data, stage, config, audit):
        assert config.grouping_unit=='individual_unverified'
        assert not audit['independent_family_support_verified']
        stage.mkdir(); (stage/'task.txt').write_text(config.task)
    monkeypatch.setattr(run_host,'run_deviation',fake_deviation)
    run_host.main(args)
    assert (scope/task/'molecular_deviation/task.txt').read_text()==task
    assert (old/'unchanged.txt').read_text()=='old objective result'


def test_legacy_model_demands_explicit_refit_without_mutation(abm):
    import copy
    model, *_, query, qi, _ = abm
    legacy = copy.copy(model); del legacy.molecular_qc_version
    with pytest.raises(RuntimeError,match='explicit refit'):
        legacy.predict(query,qi,qi,model.feature_names)


def test_successful_refit_publication_preserves_previous_snapshot(tmp_path):
    from run_host import publish
    target=tmp_path/'abm_twin';target.mkdir();(target/'fit.txt').write_text('old')
    stage=tmp_path/'stage';stage.mkdir();(stage/'fit.txt').write_text('new')
    publish(stage,target)
    assert (target/'fit.txt').read_text()=='new'
    assert [p.read_text() for p in (tmp_path/'_history').glob('*/fit.txt')]==['old']
