"""Runtime: dependency order, resources, numerical backends and console behavior."""
import importlib.util
import json
import numpy as np
import os
import pytest
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from unittest import TestCase, mock

ROOT = Path(__file__).resolve().parents[1]

def load_module(name, filename):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / 'f' / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

CONSOLE = ROOT.parent / '0f/console_run.py'
shared = load_module('le8_validation_dispatch', '0.common.py')
abm = load_module('le8_validation_abm', 'c1.abm.py')
state = load_module('le8_validation_state', '0.run_state.py')
r = load_module('le8_validation_resources', '0.resources.py')
SharedBudget, execute, atomic_json, process_identity = r.SharedBudget, r.execute, r.atomic_json, r.process_identity
BlockedKNN, validated_eigh = r.BlockedKNN, r.validated_eigh

def result(root, module, value=b"validated result"):
    for name in state.RESULTS[module]:
        file = root / "cad/prot" / name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes(value)

def test_receipts_survive_transactions_but_not_changed_outputs(tmp_path, monkeypatch):
    ledger = tmp_path / "ledger"
    ledger.mkdir()
    monkeypatch.setenv("LE8_RUN_STATE_DIR", str(ledger))
    monkeypatch.setenv("LE8_PUBLISHED_ROOT", str(tmp_path / "published"))
    first, second = tmp_path / "first", tmp_path / "second"
    result(first, "c1_correlate")
    state.receipts("record", first, ["cad"], ["prot"], ["c1_correlate"])
    shutil.copytree(first, second)
    assert state.receipts("read", second, ["cad"], ["prot"], ["c1_correlate"]) == ["cad|c1_correlate|prot"]
    assert not state.receipts("read", second, ["ra"], ["prot"], ["c1_correlate"])
    result(second, "c1_correlate", b"different numerical result")
    assert not state.receipts("read", second, ["cad"], ["prot"], ["c1_correlate"])
    monkeypatch.delenv("LE8_RUN_STATE_DIR")
    assert not state.receipts("read", first, ["cad"], ["prot"], ["c1_correlate"])

def test_receipts_require_all_c4_results(tmp_path, monkeypatch):
    monkeypatch.setenv("LE8_RUN_STATE_DIR", str(tmp_path))
    result(tmp_path, "c4_connect")
    (tmp_path / "cad/prot/c4_connect/c4.penalty.res.rds").unlink()
    with pytest.raises(ValueError, match="Missing result"):
        state.receipts("record", tmp_path, ["cad"], ["prot"], ["c4_connect"])

@pytest.mark.parametrize("replace", ["TRUE", "FALSE"])
@pytest.mark.parametrize("done,expected", [(True, "c3_coloc"), (False, "c1_correlate c2_cause c3_coloc")])
def test_native_dependency_plan_respects_current_invocation(replace, done, expected):
    source = (ROOT / "f/0.engine.sh").read_text()
    functions = source[source.index("layer_complete() {"):source.index("\nIFS=',' read -ra raw_traits")]
    script = '''set -euo pipefail
declare -A completed_results=() run_completed=() job_actions=()
replace="$REPLACE"
if [[ "$DONE" == TRUE ]]; then
  run_completed['cad|c1_correlate|prot']=TRUE
  run_completed['cad|c2_cause|prot']=TRUE
fi
''' + functions + '''
plan_job cad c3_coloc prot
for job in c1_correlate c2_cause c3_coloc; do
  [[ -z "${job_actions["cad|$job|prot"]:-}" ]] || echo "$job"
done
'''
    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True,
                          env={**os.environ, "REPLACE": replace, "DONE": str(done).upper()})
    assert " ".join(proc.stdout.split()) == expected

@pytest.mark.parametrize("fail", [False, True])
def test_dispatch_order_and_failure_boundary(tmp_path, monkeypatch, fail):
    common = load_module("run_common", "0.common.py")
    monkeypatch.delenv("LE8_TABLE_WORKSPACE", raising=False)
    monkeypatch.delenv("LE8_RUN_STATE_DIR", raising=False)
    called, ledgers = [], []

    def transaction(argv):
        module = argv[1]
        called.append(module)
        ledger = Path(os.environ["LE8_RUN_STATE_DIR"])
        assert ledger.is_dir()
        ledgers.append(ledger)
        if module == "c2_cause" and fail:
            raise RuntimeError("upstream failed")
        return True

    monkeypatch.setattr(common, "_table_runtime", transaction)
    monkeypatch.setattr(common, "dispatch_main", lambda argv: None)  # ABM preflight only.
    modules = "c1_correlate,c1_abm,c2_cause,c3_coloc,c4_connect,c4_panel_validation,c5_cellulation"
    argv = ["dispatch", modules, "--Y", "cad", "--biom", "prot", "--analysis-root", str(tmp_path)]
    if fail:
        with pytest.raises(RuntimeError, match="upstream failed"):
            common.main(argv)
    else:
        common.main(argv)
    assert called == modules.split(",")[:3 if fail else 7]
    assert len(set(ledgers)) == 1 and not ledgers[0].exists()
    assert "LE8_RUN_STATE_DIR" not in os.environ

def test_result_audit_rejects_missing_numerics_and_figures(tmp_path):
    root = tmp_path / "results"
    script = r'''
root <- commandArgs(TRUE)[1]; d <- file.path(root,'cad/prot/c1_correlate');dir.create(d,recursive=TRUE)
x <- list(meta=list(trait='cad',layer='protein',module='c1_correlate',generated='fixture'),association=data.frame(beta=.1),prevalent=data.frame(beta=.2))
saveRDS(x,file.path(d,'c1.res.rds'))
for (n in c('c1.cohort.csv','c1.directionality_triage.csv')) write.csv(data.frame(value=1),file.path(d,n),row.names=FALSE)
write.csv(data.frame(file='figure.png'),file.path(d,'figure_manifest.csv'),row.names=FALSE)
writeBin(as.raw(1:4),file.path(d,'figure.png'))
'''
    subprocess.run(["Rscript", "-e", script, str(root)], check=True, capture_output=True)
    command = ["Rscript", str(ROOT / "f/0.common.R"), "--audit-results", str(root), "cad", "prot", "c1_correlate", str(tmp_path / "audit.csv")]
    assert subprocess.run(command, capture_output=True).returncode == 0
    figure = root / "cad/prot/c1_correlate/figure.png"
    figure.unlink()
    assert subprocess.run(command, capture_output=True).returncode != 0
    figure.write_bytes(b"png")
    (figure.parent / "c1.cohort.csv").unlink()
    assert subprocess.run(command, capture_output=True).returncode != 0

def test_fingerprint_ignores_directory_writes_but_tracks_real_inputs(tmp_path):
    script = r'''
for (x in parse(commandArgs(TRUE)[1])) if (is.call(x) && identical(x[[1]],as.name('<-')) &&
 is.symbol(x[[2]]) && as.character(x[[2]]) %in% c('le8_stage_fingerprint','le8_hash_object')) eval(x)
indir <- commandArgs(TRUE)[2];dir.create(file.path(indir,'Rdata'),recursive=TRUE)
input <- file.path(indir,'Rdata/all.rds');writeBin(as.raw(1),input)
.le8_loaded_code <- 'fixture';analysis_root <- indir
Sys.setenv(LE8_PQTL_IV_DIR=indir)
.le8_stage_source_values <- list(c1_association=data.frame(beta=.1))
before <- le8_stage_fingerprint()
writeLines('new tool index',file.path(indir,'index.txt'));Sys.setFileTime(indir,Sys.time()+10)
stopifnot(identical(before,le8_stage_fingerprint()))
writeBin(as.raw(1:3),input);stopifnot(!identical(before,le8_stage_fingerprint()))
before <- le8_stage_fingerprint();.le8_stage_source_values$c1_association$beta <- .2
stopifnot(!identical(before,le8_stage_fingerprint()))
'''
    subprocess.run(["Rscript", "-e", script, str(ROOT / "f/0.common.R"), str(tmp_path / "inputs")],
                   check=True, capture_output=True, text=True)

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
    fake_r.write_text('#!/bin/bash\nif [[ $1 == -e ]]; then printf "%s" "$LE8_REPLACE" > "$RECORDED"; fi\n')
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

def test_console_partial_is_visible_in_final_summary(tmp_path):
    script=tmp_path/'le8.sh'
    script.write_text("#!/bin/bash\necho '[LE8] START c5_cellulation'\necho '[LE8] SKIP C5/CIGMA status=unavailable'\necho '[LE8] PARTIAL c5_cellulation'\n")
    proc=subprocess.run([sys.executable,str(CONSOLE),'--script',str(script),'--'],
                        env=clean_env(SCRIPT_LOG_DIR=str(tmp_path/'logs')),capture_output=True,text=True,check=True)
    assert '[LE8] SKIP C5/CIGMA' in proc.stdout
    assert '[LE8] PARTIAL c5_cellulation' in proc.stdout
    assert '部分分析未完成' in proc.stdout

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
    result = subprocess.run(['Rscript', str(ROOT / 'validation/native_models.R'), '--pgs'],
                            env=clean_env(), capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'PASS C2 decomposition' in result.stdout

@pytest.mark.parametrize('qb,db', [(1, 3), (4, 7), (5, 100)])
def test_knn_matches_bruteforce(qb, db):
    rng = np.random.default_rng(31)
    x = rng.normal(size=(25, 3))
    q = rng.normal(size=(5, 3))
    ids = np.array(['d' + str(i) for i in range(25)])
    groups = np.array(['g' + str(i // 2) for i in range(25)])
    keys = np.arange(25)[::-1]
    qids = np.array(['d0', 'new1', 'new2', 'new3', 'new4'])
    qgroups = np.array(['g0', 'z', 'g2', 'z', 'z'])
    index = BlockedKNN(x, ids, groups, keys, query_block=qb, donor_block=db)
    ix, dist = index.query(q, qids, qgroups, 4)
    for i in range(5):
        d = np.linalg.norm(x - q[i], axis=1)
        d[(ids == qids[i]) | (groups == qgroups[i])] = np.inf
        expected = np.lexsort((keys, d))[:4]
        np.testing.assert_array_equal(ix[i], expected)
        np.testing.assert_allclose(dist[i], d[expected], atol=1e-12)

def test_knn_ties_and_fold_masks():
    x = np.zeros((12, 3))
    ids = np.array([str(i) for i in range(12)])
    groups = ids.copy()
    key = np.arange(12)[::-1]
    index = BlockedKNN(x, ids, groups, key, donor_block=3)
    ix, d = index.query(x[:1], ['out'], ['out'], 5, [0], np.arange(12) % 2)
    np.testing.assert_array_equal(ix[0], [11, 9, 7, 5, 3])
    ix, d = index.query(x[:1], ['out'], ['out'], 5, [0], np.zeros(12))
    assert (ix == -1).all() and np.isinf(d).all()

@pytest.mark.parametrize('seed', [1, 2, 3])
def test_eigh_cpu(seed):
    x = np.random.default_rng(seed).normal(size=(15, 10))
    a = x.T @ x
    lam, u, diag = validated_eigh(a)
    np.testing.assert_allclose(u * lam @ u.T, a, atol=1e-10)
    assert diag['relative_residual'] < 1e-10

def test_eigh_reject_asymmetry_and_nonfinite():
    for a in [np.array([[1, 2], [0, 1]]), np.full((3, 3), np.nan)]:
        with pytest.raises(ValueError):
            validated_eigh(a)

def tasks(tmp, n=4, fail_index=None):
    out = []
    for i in range(n):
        target = tmp / f'output{i}.json'
        script = "import json,time,pathlib;time.sleep(.1);pathlib.Path(%r).write_text(json.dumps({'i':%d}))" % (str(target), i)
        if i == fail_index:
            script = 'raise SystemExit(3)'
        out.append(dict(task_id=f't{i}', argv=[sys.executable, '-c', script], scientific_signature='fixture1', outputs=[dict(path=str(target), kind='json')], cores=1, memory_gib=0.01))
    return out

def test_scheduler_parallel_order_and_resume(tmp_path):
    t = tasks(tmp_path)
    out = tmp_path / 'run'
    r = execute(t, out, cores=2, memory_gib=1, max_workers=2, pool_path=tmp_path / 'pool.json', reserve_gib=0)
    assert [x['task_id'] for x in r] == ['t0', 't1', 't2', 't3'] and all((x['status'] == 'completed' for x in r))
    r2 = execute(t, out, cores=2, memory_gib=1, max_workers=2, pool_path=tmp_path / 'pool.json', reserve_gib=0)
    assert all((x['cache_hit'] for x in r2))
    Path(t[1]['outputs'][0]['path']).write_text('corrupted')
    r3 = execute(t, out, cores=2, memory_gib=1, max_workers=2, pool_path=tmp_path / 'pool.json', reserve_gib=0)
    assert not r3[1]['cache_hit'] and r3[0]['cache_hit']

def test_scheduler_failure_does_not_mark_complete(tmp_path):
    t = tasks(tmp_path, fail_index=1)
    out = tmp_path / 'run'
    r = execute(t, out, 2, 1, 2, tmp_path / 'pool.json', 0)
    assert r[1]['status'] == 'failed' and (not (out / 'ALL_COMPLETE.json').exists())
    t2 = tasks(tmp_path)
    r2 = execute(t2, out, 2, 1, 2, tmp_path / 'pool.json', 0)
    assert r2[0]['cache_hit'] and r2[1]['status'] == 'completed'

def test_shared_budget_cannot_oversubscribe_and_reaps_dead(tmp_path):
    pool = SharedBudget(tmp_path / 'pool.json', 2, 2, 0)
    one = pool.acquire(1, 1.5)
    assert one and pool.acquire(1, 1) is None
    two = pool.acquire(1, 0.4)
    assert two and pool.acquire(1, 0.01) is None
    pool.release(one[0])
    pool.release(two[0])
    assert process_identity(os.getpid())
    atomic_json(pool.path, {'leases': {'dead': {'pid': 99999999, 'birth': None, 'cpus': [0], 'memory_gib': 2}}})
    three = pool.acquire(1, 1)
    assert three
    pool.release(three[0])

def test_resource_oversized_task_errors(tmp_path):
    t = tasks(tmp_path, 1)
    t[0]['memory_gib'] = 5
    with pytest.raises(ValueError):
        execute(t, tmp_path / 'r', 2, 1, 2, tmp_path / 'p', 0)

@pytest.mark.parametrize('operation', ['knn', 'eigh'])
def test_actual_cuda_optional(operation):
    import torch
    if not torch.cuda.is_available():
        pytest.skip('Real CUDA not present; not counted as CUDA validation')
    rng = np.random.default_rng(2)
    if operation == 'eigh':
        x = rng.normal(size=(32, 25))
        a = x.T @ x
        lc, _, _ = validated_eigh(a, 'cpu')
        lg, _, _ = validated_eigh(a, 'cuda')
        np.testing.assert_allclose(lc, lg, rtol=1e-07, atol=1e-09)
    else:
        x = rng.normal(size=(150, 8))
        ids = np.array([str(i) for i in range(150)])
        keys = np.arange(150)
        cpu = BlockedKNN(x, ids, ids, keys)
        gpu = BlockedKNN(x, ids, ids, keys, device='cuda')
        a = cpu.query(x[:10], ids[:10], ids[:10], 5)
        b = gpu.query(x[:10], ids[:10], ids[:10], 5)
        np.testing.assert_array_equal(a[0], b[0])
        np.testing.assert_allclose(a[1], b[1], rtol=1e-07, atol=1e-09)

def test_worker_affinity_and_thread_environment(tmp_path):
    output = tmp_path / 'settings.json'
    script = "import os,json,pathlib;pathlib.Path(%r).write_text(json.dumps({'omp':os.environ['OMP_NUM_THREADS'],'blas':os.environ['OPENBLAS_NUM_THREADS'],'affinity':list(os.sched_getaffinity(0))}))" % str(output)
    task = dict(task_id='resource-check', argv=[sys.executable, '-c', script], outputs=[{'path': str(output), 'kind': 'json'}], cores=1, memory_gib=0.01)
    r = execute([task], tmp_path / 'executor', 2, 1, 1, tmp_path / 'pool', 0)
    assert r[0]['status'] == 'completed'
    d = json.loads(output.read_text())
    assert d['omp'] == d['blas'] == '1' and len(d['affinity']) == 1


def test_cuda_unavailable_is_error_subprocess():
    cmd=[sys.executable,str(ROOT/'f/c1.abm.py'),'--backend','reference','--device','cuda','--check-device']
    z=subprocess.run(cmd,env={**os.environ,'CUDA_VISIBLE_DEVICES':''},capture_output=True,text=True)
    assert z.returncode!=0 and ('CUDA' in z.stderr or 'CUDA' in z.stdout)

def test_knn_cuda_nearties_reorder_masks_batches():
    rng=np.random.default_rng(62);x=rng.normal(size=(130,7));x[1]=x[0]+1e-12
    ids=np.array([f'd{i}' for i in range(len(x))]);groups=np.array([f'g{i//2}' for i in range(len(x))]);tie=np.arange(len(x))
    q=x[:9].copy();qi=ids[:9];qg=groups[:9]
    cpu=BlockedKNN(x,ids,groups,tie,query_block=2,donor_block=17).query(q,qi,qg,8)
    order=rng.permutation(len(x))
    gpu=BlockedKNN(x[order],ids[order],groups[order],tie[order],device='cuda',query_block=4,donor_block=31).query(q[::-1],qi[::-1],qg[::-1],8)
    np.testing.assert_array_equal(ids[cpu[0]],ids[order[gpu[0][::-1]]]);np.testing.assert_allclose(cpu[1],gpu[1][::-1],atol=1e-12)

def test_gpu_lease_exclusive_and_nested_no_deadlock(tmp_path,monkeypatch):
    pool=SharedBudget(tmp_path/'pool',4,2,0)
    first=pool.acquire(1,.2,gpu=True);assert first
    second=pool.acquire(1,.2);assert second
    assert pool.acquire(0,0,gpu=True) is None
    pool.release(first[0]);third=pool.acquire(0,0,gpu=True);assert third
    pool.release(third[0]);pool.release(second[0])
    monkeypatch.setenv('LE8_GPU_LEASE_HELD','1')
    with r.gpu_lease():pass

def test_public_prior_module_published_when_later_fails(tmp_path):
    fixture=tmp_path/'code';shutil.copytree(ROOT,fixture,ignore=shutil.ignore_patterns('__pycache__','.backups','.pytest_cache'))
    # The public entry point shares the installed console runtime with 0.engine.sh.
    shutil.copytree(ROOT.parent/'0f',tmp_path/'0f',ignore=shutil.ignore_patterns('__pycache__'))
    (fixture/'f/0.engine.sh').write_text('''#!/bin/bash
set -eu
while [[ $# -gt 0 ]]; do
 if [[ "$1" == --analysis-root ]]; then root="$2"; shift 2; else shift; fi
done
mkdir -p "$root/cvd_cad/prot/c1_correlate"
printf 'term,beta,p.value\nTEST,.1,.2\n' > "$root/cvd_cad/prot/c1_correlate/c1.fixture.csv"
''')
    (fixture/'f/c1.abm.py').write_text("import sys\nraise SystemExit(0 if '--preflight' in sys.argv else 7)\n")
    out=tmp_path/'results'
    z=subprocess.run([str(fixture/'le8.sh'),'c1_correlate,c1_abm','--Y','cvd_cad','--biom','prot','--analysis-root',str(out),'--cores','2'],
        env={**os.environ,'PYTHON_BIN':sys.executable,'ABM_PYTHON':sys.executable,'LE8_RESOURCE_POOL':str(tmp_path/'pool.json')},capture_output=True,text=True)
    assert z.returncode!=0
    assert list((out/'cvd_cad/prot/c1_correlate').glob('*.xlsx')),z.stdout+z.stderr

def test_public_empty_native_pipeline_budget_and_publication(tmp_path):
    fixture=tmp_path/'code';shutil.copytree(ROOT,fixture,ignore=shutil.ignore_patterns('__pycache__','.backups','.pytest_cache'))
    shutil.copytree(ROOT.parent/'0f',tmp_path/'0f',ignore=shutil.ignore_patterns('__pycache__'))
    (fixture/'f/0.engine.sh').write_text('''#!/bin/bash
set -eu
module=$1; shift
while [[ $# -gt 0 ]]; do
 if [[ "$1" == --analysis-root ]]; then root="$2"; shift 2; else shift; fi
done
[[ "$LE8_PHASE_MEMORY_GIB" == 4.0 && "$LE8_PHASE_CORES" == 2 ]]
if [[ "$module" == c2_cause ]]; then
 test -s "$root/cvd_cad/prot/c1_correlate/c1.fixture.csv"
fi
mkdir -p "$root/cvd_cad/prot/$module"
printf 'term,beta,p.value\\nTEST,.1,.2\\n' > "$root/cvd_cad/prot/$module/${module%%_*}.fixture.csv"
echo "[LE8] DONE $module"
''')
    out=tmp_path/'new-results';assert not out.exists()
    env={k:v for k,v in os.environ.items() if not k.startswith(('LE8_','SCRIPT_'))}
    env.update(PYTHON_BIN=sys.executable,LE8_RESOURCE_POOL=str(tmp_path/'pool.json'),
               LE8_RESOURCE_MEMORY_GIB='4',SCRIPT_LOG_DIR=str(tmp_path/'logs'))
    z=subprocess.run([str(fixture/'le8.sh'),'c1_correlate,c2_cause','--Y','cvd_cad','--biom','prot',
        '--analysis-root',str(out),'--replace','TRUE','--cores','2'],env=env,capture_output=True,text=True,timeout=30)
    assert z.returncode==0,z.stdout+z.stderr
    for module in ('c1_correlate','c2_cause'):
        assert list((out/'cvd_cad/prot'/module).glob('*.xlsx'))
    assert json.loads((tmp_path/'pool.json').read_text())['leases']=={}
    assert len(z.stdout.splitlines())<12
    assert 'Resource admission' not in z.stdout and '0.common.py' not in z.stdout

def test_multi_module_native_options_do_not_damage_shared_values(tmp_path,monkeypatch):
    spec=importlib.util.spec_from_file_location('final_dispatch',ROOT/'f/0.common.py');d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
    calls=[];monkeypatch.setattr(d,'dispatch_main',lambda a:calls.append(a));monkeypatch.setattr(d,'_table_runtime',lambda a:False)
    d.main(['dispatch','c1_correlate,c1_abm','--native-test-option','cvd_cad','--Y','cvd_cad','--biom','prot','--analysis-root',str(tmp_path)])
    abm=[a for a in calls if a[0]=='c1_abm']
    assert len(abm)==2
    for a in abm:assert '--native-test-option' not in a and a[a.index('--Y')+1]=='cvd_cad'
    assert '--native-test-option' in next(a for a in calls if a[0]=='c1_correlate')

class DispatchRegression(TestCase):
	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory(prefix='le8-acceptance-', dir='/tmp')
		self.addCleanup(self.tmp.cleanup)
		self.path = Path(self.tmp.name)
		self.clean = mock.patch.dict(os.environ, {k:v for k,v in os.environ.items() if not k.startswith(('LE8_', 'C4_', 'C3_', 'PGS_', 'DATE_FOLLOW_END'))}, clear=True)
		self.clean.start()
		self.addCleanup(self.clean.stop)

	def test_dispatch_only_forwards_c5_options_to_c5(self):
		calls=[]
		with mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append((list(map(str,cmd)),env))):
			shared.dispatch_main(['c1_abm,c2_cause,c5_cellulation','--Y','cad','--biom','prot','--analysis-root',str(self.path),
				'--seed','47','--end-date','2024-02-01','--Y-date','custom','--group-file','/tmp/groups.tsv',
				'--contrasts','/tmp/contrast.csv','--matched-draws','100','--cigma-results','/tmp/results.csv',
				'--cigma-cells','/tmp/cells.csv','--allow-untested-cigma','--dry-run'])
		self.assertEqual(len(calls),4)
		for cmd,env in calls:
			self.assertEqual(env['DATE_FOLLOW_END'],'2024-02-01')
			if any('c5.cellulation.py' in x for x in cmd):
				self.assertEqual(cmd[cmd.index('--seed')+1],'47')
				self.assertIn('--contrasts',cmd);self.assertIn('--cigma-cells',cmd)
			else: self.assertNotIn('--contrasts',cmd);self.assertNotIn('--cigma-cells',cmd)
			if any('c1.abm.py' in x for x in cmd): self.assertIn('--group-file',cmd);self.assertIn('--diagnosis-col',cmd)

	def test_abm_dependencies_checked_before_native_scans(self):
		calls=[]
		with mock.patch.dict(os.environ,ABM_PYTHON='/tmp/selected-python'), mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append(list(map(str,cmd)))):
			shared.dispatch_main(['c1_correlate,c1_abm','--Y','cad,ra','--biom','prot,met','--abm-backend','both',
				'--analysis-root',str(self.path),'--r-bin','/tmp/selected-Rscript','--abm-args','--tree hist','--dry-run'])
		self.assertEqual(len(calls),17)
		for cmd in calls[:8]:
			self.assertEqual(cmd[0],'/tmp/selected-python')
			self.assertIn('--preflight',cmd)
			self.assertIn('/tmp/selected-Rscript',cmd)
			self.assertIn('hist',cmd)
		self.assertIn('c1_correlate',calls[8])
		for check,run in zip(calls[:8],calls[9:]): self.assertEqual(check,run[run.index('--')+1:]+['--preflight'])

	def test_public_abm_defaults_to_cuda_transformer_with_explicit_overrides(self):
		for backend,extra in [('reference',''),('both',''),('reference','--abm-design selective --device cpu')]:
			calls=[]
			with mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append(list(map(str,cmd)))):
				shared.dispatch_main(['c1_abm','--Y','cad','--biom','prot','--analysis-root',str(self.path),
					'--abm-backend',backend,'--abm-args',extra,'--dry-run'])
			for cmd in calls:
				kind=cmd[cmd.index('--backend')+1]
				a=abm.reference_parser().parse_args(cmd[cmd.index('--backend')+2:])
				self.assertEqual(a.device,'cpu' if extra else 'cuda')
				self.assertEqual(a.abm_design,'selective' if extra or kind=='tabicl' else 'selective_attention')

	def test_abm_management_commands_execute_once(self):
		for action in ['evaluate','project','--check-device','--download-model']:
			calls=[]
			with mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append(list(map(str,cmd)))):
				shared.dispatch_main(['c1_abm','--Y','cad','--biom','met','--analysis-root',str(self.path),'--abm-args='+action])
			self.assertEqual(len(calls),1)
			self.assertNotIn('--preflight',calls[0])

	def test_missing_abm_dependency_stops_before_any_native_analysis(self):
		calls=[]
		def fail(cmd,env,dry):
			calls.append(list(map(str,cmd)))
			raise subprocess.CalledProcessError(1,cmd)
		with mock.patch.object(shared,'call',fail),self.assertRaises(subprocess.CalledProcessError):
			shared.dispatch_main(['c1_correlate,c1_abm','--Y','cad','--biom','met','--analysis-root',str(self.path)])
		self.assertEqual(len(calls),1)
		self.assertIn('--preflight',calls[0])
		self.assertFalse(any('0.engine.sh' in c for c in calls[0]))

	def test_explicit_preflight_and_skip_abm(self):
		for flags,n in [(['--preflight'],2),(['--skip-abm','--dry-run'],1)]:
			calls=[]
			with mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append(list(map(str,cmd)))):
				shared.dispatch_main(['c1_correlate,c1_abm','--Y','cad','--biom','met','--analysis-root',str(self.path),*flags])
			self.assertEqual(len(calls),n)
			if '--preflight' in flags: self.assertTrue(all('--preflight' in cmd for cmd in calls))
			else: self.assertIn('c1_correlate',calls[0])

	def test_public_shell_seed_environment_and_help(self):
		env=dict(os.environ,PYTHON_BIN=sys.executable,LE8_REPORT_PYTHON=sys.executable,SEED='63')
		r=subprocess.run(['bash',str(ROOT/'le8.sh'),'c5_cellulation','--Y','cad','--biom','prot','--dry-run','--matched-draws','100'],env=env,text=True,capture_output=True,check=True)
		self.assertIn('--seed 63',r.stdout);self.assertIn('--matched-draws 100',r.stdout)



@pytest.mark.parametrize('profile,torch_expected,tabicl_expected,pytest_expected', [
    ('--report', False, False, False),
    ('--abm', True, False, False),
    ('--all', True, True, False),
    ('--validation', True, True, True),
])
def test_install_profiles_select_dependencies(tmp_path, profile, torch_expected, tabicl_expected, pytest_expected):
    # Stub the installer process, preserving the exact arguments and temporary file.
    recorder = tmp_path / 'python'
    record = tmp_path / 'requested.json'
    recorder.write_text('#!' + sys.executable + '\n' +
        'import json,os,sys\nfrom pathlib import Path\n' +
        'Path(os.environ["INSTALL_RECORD"]).write_text(json.dumps(dict(argv=sys.argv[1:], requirements=Path(sys.argv[-1]).read_text())))\n')
    recorder.chmod(0o755)
    subprocess.run(['bash', str(ROOT / 'install.sh'), profile, '--python', str(recorder), '--no-r'],
                   env={**os.environ, 'INSTALL_RECORD': str(record)}, check=True, capture_output=True, text=True)
    result = json.loads(record.read_text())
    assert result['argv'][:4] == ['-m', 'pip', 'install', '-r']
    requirements = result['requirements'].splitlines()
    assert 'numpy==2.3.5' in requirements and 'pandas==2.2.3' in requirements
    assert ('torch==2.12.0' in requirements) == torch_expected
    assert ('tabicl==2.2.0' in requirements) == tabicl_expected
    assert ('pytest>=8,<10' in requirements) == pytest_expected
    assert not any('cigma' in line for line in requirements)
    assert not Path(result['argv'][-1]).exists()


@pytest.mark.parametrize('environment,name', [('le8', 'le8'), ('cigma', 'le8-cigma')])
def test_install_environment_uses_named_recipe(tmp_path, environment, name):
    recorder = tmp_path / 'conda'
    record = tmp_path / 'requested.json'
    recorder.write_text('#!' + sys.executable + '\n' +
        'import json,os,sys\nfrom pathlib import Path\n' +
        'Path(os.environ["INSTALL_RECORD"]).write_text(json.dumps(dict(argv=sys.argv[1:], recipe=Path(sys.argv[-1]).read_text())))\n')
    recorder.chmod(0o755)
    subprocess.run(['bash', str(ROOT / 'install.sh'), '--env', environment],
                   env={**os.environ, 'CONDA_EXE': str(recorder), 'INSTALL_RECORD': str(record)},
                   check=True, capture_output=True, text=True)
    result = json.loads(record.read_text())
    assert result['argv'][:5] == ['env', 'update', '--name', name, '--file']
    assert '\nname: ' + name + '\n' in '\n' + result['recipe']
    assert ('cigma==1.1.0' in result['recipe']) == (environment == 'cigma')
    assert ('pandas=2.2.3' in result['recipe']) == (environment == 'le8')
    assert not Path(result['argv'][-1]).exists()


def test_install_preview_does_not_execute_installers(tmp_path):
    forbidden = tmp_path / 'forbidden'
    forbidden.write_text('#!/bin/bash\nexit 99\n')
    forbidden.chmod(0o755)
    for options in [['--all', '--python', str(forbidden), '--r-bin', str(forbidden)], ['--env', 'cigma']]:
        result = subprocess.run(['bash', str(ROOT / 'install.sh'), *options, '--dry-run'],
                                env={**os.environ, 'CONDA_EXE': str(forbidden)}, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
