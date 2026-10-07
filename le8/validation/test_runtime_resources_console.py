"""Regression coverage for native phase RAM and bounded public console output."""
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONSOLE = ROOT.parent / '0f/console_run.py'


def clean_env(**values):
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(('LE8_', 'SCRIPT_', 'MRLINK2_'))}
    return {**env, **values}


def test_c2_reserves_full_native_budget(tmp_path):
    result = subprocess.run(
        ['bash', str(ROOT / 'le8.sh'), 'c2_cause', '--Y', 'cvd_cad',
         '--biom', 'prot', '--replace', 'TRUE', '--cores', '2', '--dry-run'],
        env=clean_env(PYTHON_BIN=sys.executable, LE8_RESOURCE_MEMORY_GIB='24'),
        capture_output=True, text=True, check=True)
    assert 'phase --label c2_cause --cores 2 --memory-gib 24.0' in result.stdout
    assert '--replace TRUE' in result.stdout
    assert '详细日志' not in result.stdout  # Dry-run remains directly inspectable.


@pytest.mark.parametrize('requested,expected', [('c2_cause', 'FALSE'), ('c1_correlate', 'TRUE')])
def test_engine_replace_only_selected_stage(tmp_path, requested, expected):
    source = (ROOT / 'f/0.engine.sh').read_text()
    run_one = 'run_one() {' + source.split('run_one() {', 1)[1].split('\nfor tr in ', 1)[0]
    fake_r = tmp_path / 'Rscript'
    recorded = tmp_path / 'replace.txt'
    fake_r.write_text('#!/bin/bash\nprintf "%s" "$LE8_REPLACE" > "$RECORDED"\n')
    fake_r.chmod(0o755)
    shell = '''set -euo pipefail
declare -A job_actions=(['cad|c1_correlate|prot']=analysis) trait_grch=() trait_gwas=()
biom_layers=(prot); jobs=(c1_correlate); files=(c1.correlate.R)
requested=("$REQUESTED")
job_needs_gwas() { return 1; }
replace=TRUE; dry_run=FALSE; analysis_root="$FIXTURE_ROOT"; fdir="$SOURCE_FDIR"
LE8_GRCH=auto; LE8_REFGEN_ROOT=/unused; MRLINK2_REF_POP=EUR
''' + run_one + '\nrun_one cad c1_correlate analysis\n'
    subprocess.run(['bash', '-c', shell], check=True, capture_output=True, text=True,
                   env=clean_env(REQUESTED=requested, FIXTURE_ROOT=str(tmp_path),
                                 SOURCE_FDIR=str(ROOT / 'f'), RECORDED=str(recorded),
                                 R_BIN=str(fake_r)))
    assert recorded.read_text() == expected


@pytest.mark.parametrize('mode', ['success', 'success_r', 'success_r_scope', 'success_gpu', 'failed_worker', 'oversized'])
def test_mr_subdivision_inside_real_phase(tmp_path, mode):
    """One CPU is sufficient: waiting parents must not hold a second lease."""
    pool = tmp_path / 'global.json'
    output = tmp_path / 'worker.json'
    child = tmp_path / 'child.py'
    child.write_text('''import importlib.util, json, os
from pathlib import Path
spec=importlib.util.spec_from_file_location('r',os.environ['RESOURCES'])
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
with r.gpu_lease():
    global_pool=json.loads(Path(os.environ['LE8_GLOBAL_RESOURCE_POOL']).read_text())
    leases=list(global_pool['leases'].values())
    assert len(leases)==(1 if os.environ['MODE']=='success_gpu' else 2)
    assert sum(x['gpu'] is not None for x in leases)==1
    if os.environ['MODE']=='success_gpu':
        local=json.loads(Path(os.environ['LE8_RESOURCE_POOL']).read_text())
        assert sum(x['gpu'] is not None for x in local['leases'].values())==1
    with r.gpu_lease():pass  # Reentrant use must not wait on itself.
    Path(os.environ['OUTPUT']).write_text(json.dumps(dict(
        affinity=sorted(os.sched_getaffinity(0)), budget=os.environ['LE8_RESOURCE_MEMORY_GIB'])))
raise SystemExit(7 if os.environ['MODE']=='failed_worker' else 0)
''')
    parent = tmp_path / 'parent.py'
    parent.write_text('''import importlib.util, os, sys
from pathlib import Path
spec=importlib.util.spec_from_file_location('r',os.environ['RESOURCES'])
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
global_path=os.environ['LE8_RESOURCE_POOL']
task=dict(task_id='mr', argv=[sys.executable,os.environ['CHILD']], cores=1,
          memory_gib=1 if os.environ['MODE']=='oversized' else .05,
          outputs=[dict(path=os.environ['OUTPUT'],kind='json')])
with r.phase_worker_budget(1,.5,global_path,.02) as (cores,memory,pool):
    assert pool!=global_path and 0<memory<.48 and cores==1
    status=r.execute([task],Path(os.environ['OUTPUT']).parent/'executor',cores,memory,1,pool,0)
raise SystemExit(0 if status[0]['status']=='completed' else 3)
''')
    command = [sys.executable, str(parent)]
    if mode.startswith('success_r'):
        # Evaluate the exact production R->Python invocation, including quoting.
        r_launch = '''
for (expr in parse(Sys.getenv('C2_SOURCE'))) {
  if (is.call(expr) && identical(expr[[1]],as.name('<-')) &&
      identical(expr[[2]],as.name('le8_execute_mrlink2'))) eval(expr)
}
find_launch <- function(expr) {
  if (!is.call(expr)) return(NULL)
  if (identical(expr[[1]],as.name('system2'))) return(expr)
  for (part in as.list(expr)[-1]) {
    found<-find_launch(part);if (!is.null(found)) return(found)
  }
  NULL
}
python<-Sys.getenv('PYTHON');adapter<-Sys.getenv('PARENT')
sh<-jobs_file<-ygfile<-linkdir<-'fixture argument with spaces'
launch<-find_launch(body(le8_execute_mrlink2));stopifnot(!is.null(launch))
quit(status=eval(launch))
'''
        command = ['Rscript', '-e', r_launch]
        if mode == 'success_r_scope':
            state = subprocess.run(['systemctl', '--user', 'is-system-running'], capture_output=True, text=True)
            if state.stdout.strip() not in ('running', 'degraded'):
                pytest.skip('User systemd is unavailable')
            command = ['systemd-run', '--user', '--scope', '--quiet', '--collect',
                       '-p', 'MemoryMax=512M', '-p', 'MemorySwapMax=0', '--', *command]
    result = subprocess.run(
        [sys.executable, str(ROOT / 'f/0.resources.py'), 'phase', '--label', 'c2_cause',
         '--cores', '1', '--memory-gib', '.5', *(['--gpu'] if mode=='success_gpu' else []), '--', *command],
        capture_output=True, text=True, timeout=20,
        env=clean_env(LE8_TOTAL_CORES='1', LE8_RESOURCE_MEMORY_GIB='.5',
                      LE8_RESOURCE_RESERVE_GIB='0', LE8_RESOURCE_POOL=str(pool),
                      RESOURCES=str(ROOT / 'f/0.resources.py'), OUTPUT=str(output),
                      CHILD=str(child), MODE=mode, PYTHON=sys.executable, PARENT=str(parent),
                      C2_SOURCE=str(ROOT / 'f/c2.cause.R')))
    if mode.startswith('success'):
        assert result.returncode == 0, result.stdout + result.stderr
        saved = json.loads(output.read_text())
        assert len(saved['affinity']) == 1 and float(saved['budget']) < .48
    else:
        assert result.returncode != 0
        if mode == 'oversized':
            assert 'Task exceeds resource budget' in result.stderr
            assert not output.exists()
        else:
            assert result.returncode == 3
    assert json.loads(pool.read_text())['leases'] == {}


def test_r_admission_reports_budget_without_disabling_guard():
    script = '''
for (expr in parse(commandArgs(TRUE)[1])) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      identical(expr[[2]], as.name('le8_parallel_workers'))) eval(expr)
}
stopifnot(le8_parallel_workers(2910,8,3.823,1,32,1)==8L)
err<-tryCatch(le8_parallel_workers(2910,8,3.823,1,3,1),error=conditionMessage)
stopifnot(is.character(err),grepl('budget=3.000 GiB',err,fixed=TRUE),
          grepl('parent RSS=3.823 GiB',err,fixed=TRUE))
'''
    subprocess.run(['Rscript', '-e', script, str(ROOT / 'f/0.common.R')],
                   capture_output=True, text=True, check=True)


def test_pgs_admission_accounts_for_fork_copies_and_runs_serially():
    script = '''
for (expr in parse(commandArgs(TRUE)[1])) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in%
      c('le8_parallel_workers', 'le8_worker_plan', 'le8_dynamic_map')) eval(expr)
}
N_CORES <- 16L
Sys.setenv(LE8_PHASE_CORES='16', LE8_PHASE_MEMORY_GIB='32', LE8_R_RESERVE_GIB='1')
# Replay the real failed admission without allocating a 20 GiB fixture.
readLines <- function(...) c('VmRSS: 21438136 kB', 'VmSwap: 0 kB')
for (choice in c('auto', '1', '4', '8')) {
  Sys.setenv(LE8_PGS_WORKERS=choice)
  workers <- le8_worker_plan(2582, 'pgs')
  stopifnot(workers == 1L)
  pids <- le8_dynamic_map(as.list(1:8), function(x) Sys.getpid(), workers=workers)
  stopifnot(all(unlist(pids) == Sys.getpid()))
}
# Larger budgets and smaller heaps still allow parallel execution.
Sys.setenv(LE8_PGS_WORKERS='auto', LE8_PHASE_MEMORY_GIB='128')
stopifnot(le8_worker_plan(2582, 'pgs') == 4L)
Sys.setenv(LE8_PHASE_MEMORY_GIB='32')
readLines <- function(...) c('VmRSS: 3145728 kB', 'VmSwap: 0 kB')
stopifnot(le8_worker_plan(2582, 'pgs') == 4L,
          le8_worker_plan(2910, 'pwas') == 8L)
readLines <- function(...) c('VmRSS: 8388608 kB', 'VmSwap: 0 kB')
stopifnot(le8_worker_plan(2582, 'pgs') == 2L)
readLines <- function(...) c('VmRSS: 8388608 kB', 'VmSwap: 8388608 kB')
stopifnot(le8_worker_plan(2582, 'pgs') == 1L)
# No pending work needs no admission, while even serial work must fit.
stopifnot(le8_worker_plan(0, 'pgs') == 0L)
Sys.setenv(LE8_PHASE_MEMORY_GIB='18')
err <- tryCatch(le8_worker_plan(2582, 'pgs'), error=conditionMessage)
stopifnot(is.character(err), grepl('Insufficient declared RAM', err, fixed=TRUE))
'''
    subprocess.run(['Rscript', '-e', script, str(ROOT / 'f/0.common.R')],
                   env=clean_env(), capture_output=True, text=True, check=True, timeout=30)


@pytest.mark.parametrize('fail', [False, True])
def test_console_keeps_full_log_and_one_failure_summary(tmp_path, fail):
    script = tmp_path / 'le8.sh'
    script.write_text('''#!/bin/bash
test "$SCRIPT_CONSOLE_ACTIVE" = 1 || exit 99
echo '[LE8] START C1/protein PWAS features=2910'
for i in {1..500}; do
  echo "Warning message: routine package warning $i"
  echo "[LE8] RESOURCE pwas: worker $i"
done
''' + ('''echo '[LE8] FAIL C1/protein: Insufficient declared RAM'
echo 'Error in le8_parallel_workers(...) :'
echo '  Insufficient declared RAM'
echo 'Calls: source -> ...'
echo 'Execution halted'
echo '[LE8] FAIL c1_correlate exit=1 log=/tmp/detail.log'
echo "ERROR: Command ['python', 'dispatch'] returned non-zero exit status 1."
echo '[LE8] Run failed; published results were retained. Diagnostic workspace: /tmp/example'
exit 7
''' if fail else "echo '[LE8] DONE C1/protein PWAS'\n"))
    result = subprocess.run([sys.executable, str(CONSOLE), '--script', str(script), '--'],
                            env=clean_env(SCRIPT_LOG_DIR=str(tmp_path / 'logs')),
                            capture_output=True, text=True, timeout=10)
    assert result.returncode == (7 if fail else 0)
    assert len(result.stdout.splitlines()) <= 4
    assert '警告详见日志' in result.stdout
    assert result.stderr == ''
    assert 'RESOURCE' not in result.stdout and 'Calls:' not in result.stdout
    assert 'ERROR: Command' not in result.stdout
    log = next((tmp_path / 'logs').glob('*.log')).read_text()
    assert log.count('routine package warning') == 500
    assert log.count('RESOURCE') == 500
    if fail:
        assert result.stdout.count('Insufficient declared RAM') == 1
        assert '诊断目录：/tmp/example' in result.stdout
        assert 'Calls: source' in log and 'ERROR: Command' in log


def test_console_handles_unwrapped_r_error():
    spec = importlib.util.spec_from_file_location('console', CONSOLE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    console = module.Le8Console()
    assert console.feed('Error in fit(x) :') is None
    assert console.feed('  root failure') is None
    assert console.failure.strip() == 'root failure'


def test_console_keeps_model_reuse_event():
    spec = importlib.util.spec_from_file_location('console', CONSOLE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    console = module.Le8Console()
    line = '[C1 selective] SKIP completed | /tmp/synthetic/model'
    assert console.feed(line) == line
    assert console.feed(line) is None
    for event in ['START S7_crossfit | fold=1', 'DONE S7_final_model | attention_retrieval_full']:
        assert console.feed('[ABM] ' + event) == '[ABM] ' + event
    for epoch in range(500):
        assert console.feed(f'[ABM] INFO S7_epoch | epoch={epoch} loss=.4') is None


def test_native_baseline_audit_is_not_an_analysis_workbook(tmp_path):
    root = tmp_path / 'logs' / 'analysis'
    stage = root / 'cvd_cad/prot/c1_correlate'
    stage.mkdir(parents=True)
    (stage / 'c1.fixture.csv').write_text('term,beta,p.value\nPCSK9,.2,.01\n')
    logs = root / 'logs/cvd_cad'
    logs.mkdir(parents=True)
    audit = logs / 'baseline_rebuild_audit.csv'
    audit.write_text('variable,N,changed\nbp.pts,2400,2400\n')
    original = audit.read_bytes()
    result = subprocess.run(['Rscript', str(ROOT / 'f/0.common.R'), '--tables-pack', str(root)],
                            env=clean_env(LE8_TABLE_WORKSPACE='1'), text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert list(stage.glob('*.xlsx'))
    assert audit.read_bytes() == original and not list(logs.glob('*.xlsx'))


def test_invalid_native_white_only_rejected_before_execution(tmp_path):
    result = subprocess.run(['bash', str(ROOT / 'f/0.engine.sh'), 'c1_correlate',
                             '--white-only', 'invalid', '--analysis-root', str(tmp_path / 'results')],
                            env=clean_env(SCRIPT_CONSOLE_ACTIVE='1'), text=True, capture_output=True)
    assert result.returncode == 2
    assert '--white-only requires TRUE/FALSE' in result.stderr
    assert not (tmp_path / 'results').exists()


def test_phase_cancellation_reaps_detached_workers_and_releases_budget(tmp_path):
    pool = tmp_path / 'pool.json'
    ready = tmp_path / 'ready.json'
    parent = tmp_path / 'parent.py'
    parent.write_text('''import importlib.util, os, sys
from pathlib import Path
spec=importlib.util.spec_from_file_location('r',os.environ['RESOURCES'])
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
code="import os,json,time;from pathlib import Path;Path(os.environ['READY']).write_text(json.dumps(dict(pid=os.getpid())));time.sleep(300)"
task=dict(task_id='pending',argv=[sys.executable,'-c',code],cores=1,memory_gib=.05,
          outputs=[dict(path=os.environ['READY'],kind='json')])
with r.phase_worker_budget(1,.5,os.environ['LE8_RESOURCE_POOL'],.02) as (cores,memory,local):
    r.execute([task],Path(os.environ['READY']).parent/'executor',cores,memory,1,local,0)
''')
    proc = subprocess.Popen([sys.executable, str(ROOT / 'f/0.resources.py'), 'phase',
                             '--label', 'cancel-fixture', '--cores', '1', '--memory-gib', '.5',
                             '--', sys.executable, str(parent)], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, start_new_session=True,
                            env=clean_env(LE8_TOTAL_CORES='1', LE8_RESOURCE_MEMORY_GIB='.5',
                                          LE8_RESOURCE_RESERVE_GIB='0', LE8_RESOURCE_POOL=str(pool),
                                          RESOURCES=str(ROOT / 'f/0.resources.py'), READY=str(ready)))
    child_pid = None
    try:
        deadline = time.monotonic() + 10
        while not ready.exists() and proc.poll() is None and time.monotonic() < deadline:
            time.sleep(.05)
        assert ready.exists(), proc.communicate(timeout=2)[0]
        child_pid = json.loads(ready.read_text())['pid']
        proc.send_signal(signal.SIGTERM)
        output = proc.communicate(timeout=12)[0]
        assert proc.returncode != 0, output
        status = Path(f'/proc/{child_pid}/stat')
        assert not status.exists() or status.read_text().rsplit(')', 1)[1].split()[0] == 'Z'
        assert json.loads(pool.read_text())['leases'] == {}
        assert not (tmp_path / 'executor/ALL_COMPLETE.json').exists()
    finally:
        if proc.poll() is None:
            proc.kill();proc.wait(timeout=5)
        if child_pid is not None:
            try:os.kill(child_pid,signal.SIGKILL)
            except ProcessLookupError:pass


def test_pgs_column_names_preserve_c1_c2_numerics():
    result = subprocess.run(['Rscript', str(ROOT / 'validation/test_pgs_join.R')],
                            env=clean_env(), capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'PASS C2 decomposition' in result.stdout
