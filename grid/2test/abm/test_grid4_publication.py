import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from runtime_state import group_contract, file_identity
from prepare_kinship import prepare, activation_locks
from kinship import active, compatible, publication_metadata
from final_outputs import current_abm_run

ROOT = Path(__file__).resolve().parents[2]


def cohort(root):
    cache = root / 'prsformer/cache'
    cache.mkdir(parents=True)
    data = pd.DataFrame(dict(eid=['a','b','c','d'], group=['a','b','c','d'], split=['train','train','test','test']))
    data.to_csv(cache / 'data.tsv.gz', sep='\t', index=False)
    source = root / 'source.dat'
    source.write_text('ID1 ID2 Kinship\na c .1\n')
    return source


def test_group_contract_exact_changes_and_resume(tmp_path):
    frame = pd.DataFrame(dict(id=['a','b','c','d'], group_id=['1','1','2','2'], role=['build','build','test','test']))
    meta = dict(groups_verified=True, grouping_status='verified_from_source', kinship_version='v', kinship_token='t')
    first = group_contract(frame, meta, tmp_path)
    assert first == group_contract(frame.copy(), meta, tmp_path)
    changed = frame.copy(); changed.group_id = ['2','2','1','1']
    second = group_contract(changed, meta, tmp_path)
    assert first['n_rows'] == second['n_rows'] and first['n_groups'] == second['n_groups']
    assert first['state_token'] != second['state_token']
    changed.role = ['test','test','build','build']
    third = group_contract(changed, meta, tmp_path)
    assert third['role_counts'] == second['role_counts'] and third['state_token'] != second['state_token']
    assert third == group_contract(changed.copy(), meta, tmp_path)
    assert first == group_contract(frame, meta, tmp_path)
    assert first != group_contract(frame, dict(meta, kinship_token="new"), tmp_path)


@pytest.mark.parametrize('name', ['prsformer/.resume.lock', 'prsformer2/.resume.lock', 'abm/FINAL/height/run/.resume.lock', 'cache/abm/.dosage.lock', 'cache/abm/blocks/.dosage.lock', 'cache/abm/.publication.lock', 'workflow/.resume.lock'])
def test_activation_busy_never_moves_outputs(tmp_path, name):
    source = cohort(tmp_path)
    stale = tmp_path / 'final/old.tsv'; stale.parent.mkdir(); stale.write_text('keep')
    lock = tmp_path / name; lock.parent.mkdir(parents=True, exist_ok=True)
    with lock.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError, match='Kinship activation refused'):
            prepare(source, tmp_path, activate=True)
    assert stale.read_text() == 'keep'
    assert not (tmp_path / 'local/kinship/ACTIVE.json').exists()
    assert not list((tmp_path / 'local/kinship').rglob('activation_actions.json'))


def test_activation_commit_last_and_qc_private(tmp_path):
    source = cohort(tmp_path)
    import prepare_kinship
    events = []
    real = prepare_kinship.atomic_json
    def record(path, data):
        events.append(Path(path).name)
        return real(path, data)
    with patch.object(prepare_kinship, 'atomic_json', side_effect=record):
        prepare(source, tmp_path, activate=True)
    assert events[-1] == 'ACTIVE.json'
    state = active(tmp_path)
    assert state['token'] and compatible(publication_metadata(tmp_path), state)
    from kinship_qc import public_qc, COLUMNS
    table = public_qc(tmp_path)
    assert list(table) == COLUMNS
    assert table.iloc[0].N_development_excluded == 1
    assert table.iloc[0].N_test_retained == 2
    assert not any('sha' in col or col in ['id','eid','group_id'] for col in table)
    assert not any('/' in str(v) for v in table.iloc[0])


def test_activation_failure_keeps_no_partial_active(tmp_path):
    source = cohort(tmp_path)
    final = tmp_path / 'final'; final.mkdir()
    (final / 'one.tsv').write_text('one'); (final / 'two.tsv').write_text('two')
    original = Path.rename
    calls = []
    def rename(path, target):
        calls.append(path.name)
        if path.name == 'two.tsv' and path.parent == final:
            raise OSError('injected failure')
        return original(path, target)
    with patch.object(Path, 'rename', rename), pytest.raises(OSError, match='injected'):
        prepare(source, tmp_path, activate=True)
    assert (final / 'one.tsv').read_text() == 'one'
    assert (final / 'two.tsv').read_text() == 'two'
    assert not (tmp_path / 'local/kinship/ACTIVE.json').exists()
    assert list((tmp_path / 'local/kinship').rglob('activation_failure.json'))


def test_active_kinship_requires_verified_metadata(tmp_path):
    source = cohort(tmp_path); prepare(source, tmp_path, activate=True)
    from test_final_results import abm_fixture
    abm_fixture(tmp_path)
    pointer = tmp_path / 'abm/FINAL/height/blocknet.CURRENT.json'
    assert current_abm_run(pointer) is None
    run = pointer.parent / 'blocknet_fixture'
    from abm_core.model import file_sha256
    for verified in (False, True):
        manifest = json.loads((run / 'RUN_MANIFEST.json').read_text())
        manifest['metadata'] = dict(publication_metadata(tmp_path), groups_verified=verified)
        (run / 'RUN_MANIFEST.json').write_text(json.dumps(manifest))
        receipt = json.loads((run / 'COMPLETE.json').read_text())
        receipt['manifest_sha256'] = file_sha256(run / 'RUN_MANIFEST.json')
        (run / 'COMPLETE.json').write_text(json.dumps(receipt))
        phases = json.loads((run / 'PHASES.json').read_text())
        phases['phases']['report']['COMPLETE.json'] = file_identity(run / 'COMPLETE.json')
        (run / 'PHASES.json').write_text(json.dumps(phases))
        assert (current_abm_run(pointer) is not None) == verified


def test_main_only_blocknet_and_sensitivity_separate(tmp_path):
    from test_final_results import abm_fixture
    from final_compare import run
    abm_fixture(tmp_path)
    parent = tmp_path / 'abm/FINAL/height'
    (parent / 'moe.CURRENT.json').write_text(json.dumps(dict(run='blocknet_fixture')))
    run(tmp_path, ['height'])
    assert len(pd.read_csv(tmp_path / 'final/performance.tsv', sep='\t')) == 1
    assert len(pd.read_csv(tmp_path / 'final/sensitivity_performance.tsv', sep='\t')) == 1


def shell(tmp_path, *args):
    return subprocess.run([str(ROOT / 'grid.sh'), *args], capture_output=True, text=True,
        env=dict(os.environ, GRID_PYTHON=sys.executable, SCRIPT_LOG_DIR=str(tmp_path / 'logs'), SCRIPT_VERBOSE='1'))


@pytest.mark.parametrize('module', ['blocknet','twins','assess','abm','all'])
def test_public_workflow_frozen(tmp_path, module):
    for flags in (['--predictor','moe'], ['--representation=csx_scores'], ['--through','prepare']):
        result = shell(tmp_path, module, *flags, '--dry-run')
        assert result.returncode != 0
        assert 'frozen' in result.stderr + result.stdout or '--through is internal' in result.stderr + result.stdout


def test_final_filters_training_options(tmp_path):
    result = shell(tmp_path, 'final', '--traits','height','--analysis-root',str(tmp_path / 'absent'),
        '--device','auto','--predictor','moe','--representation','csx_scores','--threads','4','--dry-run')
    assert result.returncode == 0, result.stderr
    assert 'paper_summary.py' in result.stdout
    assert not (tmp_path / 'absent').exists()


def test_invalid_disco_is_unavailable_without_fitting(tmp_path):
    from final_compare import current_tables
    out = tmp_path / 'disco/height'; out.mkdir(parents=True)
    (out / 'performance.tsv').write_text('trait\tmethod\nheight\tDiscoDivas\n')
    with patch('abm_core.predictors.fit_linear', side_effect=AssertionError('final cannot fit')):
        tables, status = current_tables(tmp_path, ['height'])
    assert not tables
    assert next(row for row in status if row['method'] == 'DiscoDivas')['status'] == 'unavailable'


@pytest.mark.parametrize('trait', ['height', 't2dm'])
def test_disco_adapter_development_only_then_readonly_final(tmp_path, trait):
    import disco_evaluate
    from disco_publication import read_performance
    from abm_core.predictors import fit_linear
    cache = tmp_path / 'prsformer/cache'; cache.mkdir(parents=True)
    n = 250
    ids = [str(i) for i in range(n)]
    roles = ['build']*80 + ['tune_model']*60 + ['prediction_calibration']*60 + ['test']*50
    roster = pd.DataFrame(dict(eid=ids, group=ids, split=['train']*200+['test']*50))
    roster.to_csv(cache / 'data.tsv.gz', sep='\t', index=False)
    source = tmp_path / 'source.dat'; source.write_text('ID1 ID2 Kinship\n0 1 .1\n')
    prepare(source, tmp_path, activate=True)
    groups = ids.copy(); groups[1] = '0'
    rng = np.random.default_rng(2)
    score = rng.normal(size=n)
    y = np.arange(n) % 2 if trait == 't2dm' else 2*score+rng.normal(size=n)
    frame = pd.DataFrame(dict(id=ids, group_id=groups, role=roles, y=y, age=rng.normal(size=n)))
    score_path = tmp_path / 'scores.tsv'
    pd.DataFrame(dict(eid=ids, disco=score)).to_csv(score_path, sep='\t', index=False)
    prepared = tmp_path / 'cache/abm/prepared.joblib'; prepared.write_bytes(b'synthetic prepared identity')
    with patch.object(disco_evaluate, 'load_grid', return_value=(frame, ['age'], [], publication_metadata(tmp_path))), \
         patch.object(disco_evaluate, 'fit_linear', wraps=fit_linear) as fit:
        disco_evaluate.evaluate(tmp_path, trait, score_path)
    assert np.isnan(fit.call_args.args[1][200:]).all()
    table, reason = read_performance(tmp_path, trait)
    assert table is not None, reason
    assert table.iloc[0]['n'] == 50
    assert 'id' not in table and 'group_id' not in table
    from final_compare import current_tables
    with patch.object(disco_evaluate, 'evaluate', side_effect=AssertionError('final must only read')):
        tables, status = current_tables(tmp_path, [trait])
    assert len(tables) == 1
    score_path.write_text('changed')
    assert read_performance(tmp_path, trait)[0] is None


def test_public_workflow_duplicate_refused_before_loading_inputs(tmp_path):
    root = tmp_path / 'analysis'
    folder = root / 'cache/abm'; folder.mkdir(parents=True)
    with (folder / '.publication.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = shell(tmp_path, 'all', '--analysis-root', str(root), '--traits', 'height')
    assert result.returncode != 0
    assert 'Another GRID publication workflow' in result.stdout + result.stderr


def test_paper_distinguishes_six_questions_and_method_classes(tmp_path):
    from test_final_results import abm_fixture
    from abm_core.model import file_sha256
    from paper_summary import main
    abm_fixture(tmp_path)
    run = tmp_path / 'abm/FINAL/height/blocknet_fixture'
    pd.DataFrame([dict(trait='height',method=m,prediction_r2=r) for m,r in
        [('csx-linear',.1),('csx-blocknet',.2),('ABM(base=csx-linear,candidate=csx-blocknet)',.15),('Clinical',.05)]]).to_csv(run/'main_performance.tsv',sep='\t',index=False)
    curve = {'coverage.nominal_development_coverage':.5}
    for prefix in ('matched.predictability_loss_difference','matched.gain_enrichment',
                   'whole_cohort.policy_gain_vs_base','whole_cohort.policy_gain_vs_candidate'):
        curve.update({prefix+'.estimate':0.,prefix+'.ci_low':-.1,prefix+'.ci_high':.1})
    for name in ('predictability_curve.tsv','benefit_curve.tsv'):
        pd.DataFrame([curve]).to_csv(run/name,sep='\t',index=False)
    pd.DataFrame([dict(contrast='D_vs_C',objective='expected_loss',mse_reduction=0.,status='ok')]).to_csv(run/'count_incremental_contrasts.tsv',sep='\t',index=False)
    phases=json.loads((run/'PHASES.json').read_text())
    phases['phases']['evaluate']['main_performance.tsv']=file_sha256(run/'main_performance.tsv')
    (run/'PHASES.json').write_text(json.dumps(phases))
    main(['--analysis-root',str(tmp_path),'--traits','height'])
    keys=pd.read_csv(tmp_path/'final/paper_key_results.tsv',sep='\t')
    assert set(keys.result)=={'blocknet_vs_csx','predictability_selected_vs_control','count_incremental_D_vs_C',
                              'gain_enrichment','policy_vs_base','policy_vs_candidate'}
    table=pd.read_csv(tmp_path/'final/paper_main_performance.tsv',sep='\t')
    assert set(table.method_class)=={'baseline','block_model','policy','benchmark'}
