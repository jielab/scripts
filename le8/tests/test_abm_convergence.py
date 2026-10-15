"""Convergence retries preserve the estimator and never publish an unfinished fit."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'f'))
from abm_twin_core import Config, fit_classifier


@pytest.fixture
def weighted_data():
    rng = np.random.default_rng(61)
    x = rng.normal(size=(500, 8))
    y = (rng.random(len(x)) < 1 / (1 + np.exp(-x[:, 0] + .6*x[:, 1]))).astype(int)
    w = rng.uniform(.5, 1.5, len(x))
    w[::9] = 0
    # These rows must not influence fitting or its normalization on any retry.
    x[w == 0] = 1e6
    return x, y, w


def direct_fit(data, budget, seed=2026):
    x, y, w = data
    ok = w > 0
    return LogisticRegression(C=.1, solver='saga', penalty='elasticnet', l1_ratio=.5,
                              max_iter=budget, random_state=seed).fit(
        x[ok], y[ok], sample_weight=w[ok]/w[ok].mean())


def test_converged_first_attempt_is_numerically_unchanged(weighted_data, capsys):
    cfg = Config(seed=2026, max_iter=1000)
    expected = direct_fit(weighted_data, cfg.max_iter)
    actual = fit_classifier(*weighted_data, cfg)
    np.testing.assert_array_equal(actual.coef_, expected.coef_)
    np.testing.assert_array_equal(actual.intercept_, expected.intercept_)
    assert actual.convergence_audit_ == [dict(max_iter=1000, n_iter=int(actual.n_iter_[0]), converged=True)]
    assert not capsys.readouterr().out


def test_real_nonconvergence_continues_to_convergence(weighted_data, capsys):
    cfg = Config(seed=2026, max_iter=5, max_iter_retries=4)
    # The same initial budget previously terminated the run.
    with pytest.raises(RuntimeError, match='did not converge'):
        fit_classifier(*weighted_data, replace(cfg, max_iter_retries=0))
    actual = fit_classifier(*weighted_data, cfg)
    audit = actual.convergence_audit_
    assert len(audit) > 1 and audit[-1]['converged']
    assert all(not row['converged'] for row in audit[:-1])
    assert [row['max_iter'] for row in audit] == [5*2**i for i in range(len(audit))]
    expected = direct_fit(weighted_data, 1000)
    x = weighted_data[0][weighted_data[2] > 0]
    np.testing.assert_allclose(actual.predict_proba(x), expected.predict_proba(x), atol=2e-4, rtol=0)
    assert actual.C == expected.C and actual.tol == expected.tol
    assert actual.l1_ratio == expected.l1_ratio and actual.solver == expected.solver
    log = capsys.readouterr().out
    assert 'convergence retry' in log and 'DONE ABM elastic-net converged' in log


def test_exhausted_budget_fails_explicitly(weighted_data):
    with pytest.raises(RuntimeError, match=r'iterations=7, attempts=3.*--abm-config'):
        fit_classifier(*weighted_data, Config(seed=2026, max_iter=1, max_iter_retries=2))


def test_invalid_input_is_not_retried(weighted_data, capsys):
    x, y, w = weighted_data
    x[1, 0] = np.nan
    with pytest.raises(ValueError, match='NaN'):
        fit_classifier(x, y, w, Config())
    assert not capsys.readouterr().out


@pytest.mark.parametrize('values', [dict(max_iter=0), dict(max_iter=1.5), dict(max_iter=True),
                                  dict(max_iter_retries=-1), dict(max_iter_retries=6),
                                  dict(max_iter_retries=.5), dict(max_iter_retries=True)])
def test_invalid_iteration_budget_is_rejected(values):
    with pytest.raises(ValueError, match='max_iter'):
        Config(**values).validate()


def test_host_reports_root_cause_and_preserves_published_result(tmp_path, monkeypatch, weighted_data, capsys):
    import host_adapter
    import run_host
    cfg = Config(max_iter=1, max_iter_retries=0)
    config_file = tmp_path/'abm.json'
    config_file.write_text(json.dumps(dict(max_iter=1, max_iter_retries=0)))
    args = run_host.parser().parse_args([
        '--Y', 'cvd_cad', '--biom', 'met', '--analysis-root', str(tmp_path/'results'),
        '--config', str(config_file), '--replace', 'TRUE'])
    scope = args.analysis_root/'cvd_cad/met/c1_abm'
    target = scope/'abm_twin'
    target.mkdir(parents=True)
    (target/'old_fit.txt').write_text('published result')

    def prepare(args, work):
        (work/'prepared.txt').write_text('retained diagnostics')
        return {}, {}, dict(independent_family_support_verified=False)

    def fail_fit(**kwargs):
        return fit_classifier(*weighted_data, cfg)

    monkeypatch.setattr(run_host, 'prepare_host', prepare)
    monkeypatch.setattr(host_adapter, 'run_prepared_numeric_cohort', fail_fit)
    with pytest.raises(RuntimeError, match='did not converge'):
        run_host.main(args)
    assert (target/'old_fit.txt').read_text() == 'published result'
    stages = list(scope.glob('_abm_staging_*'))
    assert len(stages) == 1 and (stages[0]/'prepared.txt').is_file()
    assert not list(stages[0].rglob('COMPLETE.json'))

    # Exercise the actual shared console parser without changing it.
    spec = importlib.util.spec_from_file_location('le8_console_test', ROOT.parent/'0f/console_run.py')
    console_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(console_module)
    console = console_module.Le8Console()
    for line in capsys.readouterr().err.splitlines():
        console.feed(line)
    console.feed("ERROR: Command ['python', 'run_host.py'] returned non-zero exit status 1.")
    assert 'ABM elastic-net did not converge' in console.failure
    assert 'biom=met' in console.failure
