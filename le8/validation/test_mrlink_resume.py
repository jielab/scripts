"""MR checkpoints must survive reuse from a different module transaction."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
ADAPTER = ROOT / 'f/c2.parallel.py'
spec = importlib.util.spec_from_file_location('mr_resume', ADAPTER)
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)


def checkpoint_task(tmp_path, status='ok'):
    out = tmp_path / 'results'
    marker = tmp_path / 'validated.json'
    request = tmp_path / 'worker.json'
    calls = tmp_path / 'calls.txt'
    # Real subprocess, tiny deterministic result and numerical audit.
    code = f'''
from pathlib import Path
out = Path({str(out)!r})
with Path({str(calls)!r}).open('a') as handle:
    handle.write('called\\n')
result = out / 'estimate.tsv'
result.write_text('alpha\\tse(alpha)\\tp(alpha)\\n0.1\\t0.03\\t0.02\\n')
(out / 'fixture_numerical_audit.json').write_text('{{"ordered_snps": ["rs1"]}}')
(out / 'mrlink2.status.tsv').write_text(
    'omics\\ttrait\\tstatus\\tout_prefix\\tmessage\\n'
    'protein\\tP\\t{status}\\t' + str(result) + '\\tcompleted\\n')
'''
    mr.atomic_json(request, dict(task_id='fixture', scientific_signature='science-v1',
                                omics='protein', trait='P', worker_outdir=str(out),
                                argv=[sys.executable, '-c', code]))
    return dict(task_id='fixture', scientific_signature='science-v1',
                argv=[sys.executable, str(ADAPTER), 'worker', '--task', str(request),
                      '--marker', str(marker)], cores=1, memory_gib=.01,
                outputs=[dict(path=str(marker), kind='mrlink_marker')])


def run_task(task, tmp_path, transaction):
    return mr.execute([task], tmp_path / transaction, cores=1, memory_gib=.1,
                      max_workers=1, pool_path=tmp_path / 'pool.json', reserve_gib=0)[0]


@pytest.mark.parametrize('status', ['ok', 'no_estimate'])
def test_checkpoint_reused_across_transaction_directories(tmp_path, status):
    task = checkpoint_task(tmp_path, status)
    assert run_task(task, tmp_path, 'c2')['status'] == 'completed'
    marker = Path(task['outputs'][0]['path'])
    before = marker.read_bytes(), marker.stat().st_mtime_ns
    second = run_task(task, tmp_path, 'c3')
    assert second['status'] == 'completed', second
    assert (tmp_path / 'c3/ALL_COMPLETE.json').is_file()
    assert (marker.read_bytes(), marker.stat().st_mtime_ns) == before
    assert (tmp_path / 'calls.txt').read_text() == 'called\n'
    assert run_task(task, tmp_path, 'c3')['cache_hit']


@pytest.mark.parametrize('corruption', [
    'result', 'audit', 'missing_audits', 'task_id', 'scientific_signature', 'status_row',
])
def test_invalid_checkpoint_fails_without_destroying_evidence(tmp_path, corruption):
    task = checkpoint_task(tmp_path)
    assert run_task(task, tmp_path, 'c2')['status'] == 'completed'
    marker = Path(task['outputs'][0]['path'])
    record = json.loads(marker.read_text())
    if corruption == 'result':
        Path(record['row']['out_prefix']).write_text('corrupt')
    elif corruption == 'audit':
        Path(record['numerical_audits'][0]['path']).write_text('corrupt')
    elif corruption == 'missing_audits':
        record['numerical_audits'] = []
    elif corruption in ('task_id', 'scientific_signature'):
        record[corruption] = 'different-task'
    else:
        path = tmp_path / 'results/mrlink2.status.tsv'
        path.write_text(path.read_text().replace('protein\tP\t', 'protein\tOTHER\t'))
    mr.atomic_json(marker, record)
    before = marker.read_bytes()
    result = run_task(task, tmp_path, 'c3')
    assert result['status'] == 'failed'
    assert not (tmp_path / 'c3/ALL_COMPLETE.json').exists()
    assert marker.read_bytes() == before
    assert (tmp_path / 'calls.txt').read_text() == 'called\n'


def test_unvalidated_output_directory_is_not_promoted(tmp_path):
    task = checkpoint_task(tmp_path)
    assert run_task(task, tmp_path, 'c2')['status'] == 'completed'
    marker = Path(task['outputs'][0]['path'])
    marker.unlink()
    assert run_task(task, tmp_path, 'c3')['status'] == 'failed'
    assert not marker.exists()
    assert (tmp_path / 'calls.txt').read_text() == 'called\n'
    log = tmp_path / 'c3' / (hashlib.sha256(task['task_id'].encode()).hexdigest() + '.log')
    assert 'create a new isolated attempt' in log.read_text()
