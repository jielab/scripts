"""FINAL production kernels and failure contracts; synthetic fixtures, real CUDA."""
import dataclasses, hashlib, importlib.util, json, os, subprocess, sys, time, csv, shutil
from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np
import pandas as pd
import pytest
from test_c1_s7_production import m, config, fitted
ROOT=Path(__file__).resolve().parents[1]
r=m.s7_resources()
SharedBudget,execute,atomic_json,process_identity=r.SharedBudget,r.execute,r.atomic_json,r.process_identity
BlockedKNN,validated_eigh=r.BlockedKNN,r.validated_eigh
spec=importlib.util.spec_from_file_location('test_final_mr',ROOT/'f/c2.parallel.py')
mr=importlib.util.module_from_spec(spec);spec.loader.exec_module(mr)
validate_status=mr.validate_status
@dataclasses.dataclass
class Batch:
    bulk: np.ndarray
    tail_up: np.ndarray
    tail_down: np.ndarray
    observed: np.ndarray
    clinical: np.ndarray
    ids: np.ndarray
    groups: np.ndarray
    feature_names: list

    def design(self):
        return np.c_[self.bulk, self.tail_up, self.tail_down, (~self.observed).astype(float), self.clinical]

    def take(self, ix):
        return Batch(**{k: v if k == 'feature_names' else v[ix] for k, v in vars(self).items()})

def bank(n=1000, kind='binary'):
    x = np.zeros((n, 2))
    up = np.zeros_like(x)
    if kind == 'binary':
        x[:30, 0] = 1
    else:
        up[:30, 0] = 2
    return Batch(x, up, np.zeros_like(x), np.ones_like(x, dtype=bool), np.empty((n, 0)), np.array(['id' + str(i) for i in range(n)]), np.array(['g' + str(i) for i in range(n)]), ['a', 'b'])

def cfg():
    return dict(seed=19, s7_views=['learned', 'tail'], s7_radius_ref_n=50, s7_neighbors_per_view=8)

@pytest.mark.parametrize('kind', ['binary', 'continuous'])
def test_exact_rare_survives_zero_radius(kind):
    b = bank(kind=kind)
    z = np.zeros((len(b.ids), 2))
    g = m.S7RetrievalGeometry().fit(b, z, cfg())
    q = b.take([0])
    q.ids[:] = 'external'
    q.groups[:] = 'external'
    ix, r, ok, d = g.candidates(q, np.zeros((1, 2)))
    assert g.status['tail'] == 'exact_match_only'
    assert ok.sum() == 8 and np.all(ix[ok] < 30) and np.all(r[ok] == 0)
    assert d.tail_available.iloc[0] and (not d.learned_available.iloc[0])

@pytest.mark.parametrize('kind', ['binary', 'continuous'])
def test_exact_does_not_match_missing_or_other_assay(kind):
    b = bank(kind=kind)
    g = m.S7RetrievalGeometry().fit(b, np.zeros((len(b.ids), 2)), cfg())
    q = b.take([0])
    q.ids[:] = 'external'
    q.groups[:] = 'external'
    q.observed[:, 0] = False
    assert not g.candidates(q, np.zeros((1, 2)))[2].any()
    q.observed[:] = True
    q.bulk[:] = 0
    q.tail_up[:] = 0
    q.tail_up[:, 1] = 2
    assert not g.candidates(q, np.zeros((1, 2)))[2].any()

@pytest.mark.parametrize('exclude', ['family', 'fold'])
def test_exact_excludes_family_and_fold(exclude):
    b = bank()
    g = m.S7RetrievalGeometry().fit(b, np.zeros((len(b.ids), 2)), cfg())
    q = b.take([0])
    q.ids[:] = 'outside'
    if exclude == 'family':
        g.groups[:30] = 'F'
        q.groups[:] = 'F'
        args = {}
    else:
        q.groups[:] = 'new'
        bf = np.zeros(len(b.ids), int)
        bf[:30] = 1
        args = dict(query_folds=np.array([1]), bank_folds=bf)
    assert not g.candidates(q, np.zeros((1, 2)), **args)[2].any()

def test_constant_global_not_mislabeled_exact():
    b = bank()
    b.bulk[:] = 0
    b.tail_up[:] = 0
    g = m.S7RetrievalGeometry().fit(b, np.zeros((len(b.ids), 2)), cfg())
    assert g.status['learned'] == 'constant_geometry' and g.status['tail'] == 'unavailable'

def test_static_refresh_preserves_indexes_and_radii():
    b = bank()
    rng = np.random.default_rng(1)
    z = rng.normal(size=(len(b.ids), 2))
    g = m.S7RetrievalGeometry().fit(b, z, cfg())
    tail = g.tail_index
    radius = g.radii['tail']
    g.refresh_learned(b, z * 2)
    assert g.tail_index is tail and g.radii['tail'] == radius and (g.static_build_count == 1)
    with pytest.raises(ValueError):
        g.refresh_learned(b.take(np.arange(len(b.ids))[::-1]), z)

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

@pytest.mark.parametrize('status', ['failed', 'not_run'])
def test_mrlink_exit_zero_status_still_fails(tmp_path, status):
    f = tmp_path / 'status.tsv'
    f.write_text('omics\ttrait\tstatus\nout\ta\t' + status + '\n')
    with pytest.raises(ValueError):
        validate_status(f, 'out', 'a')

def test_mrlink_no_estimate_and_good_result(tmp_path):
    f = tmp_path / 'status.tsv'
    f.write_text('omics\ttrait\tstatus\nout\ta\tno_estimate\n')
    assert validate_status(f, 'out', 'a')['status'] == 'no_estimate'
    result = tmp_path / 'result'
    result.write_text('alpha\tse(alpha)\tp(alpha)\n0.1\t0.03\t0.02\n')
    f.write_text('omics\ttrait\tstatus\tout_prefix\nout\ta\tok\t' + str(result) + '\n')
    assert validate_status(f, 'out', 'a')['result_sha256']

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

def test_mrlink_marker_checks_underlying_result(tmp_path):
    validate_outputs, digest_file = r.validate_outputs, r.digest_file
    target = tmp_path / 'result'
    target.write_text('alpha\tp(alpha)\n.1\t.2\n')
    marker = tmp_path / 'marker.json'
    atomic_json(marker, {'status': 'ok', 'row': {'out_prefix': str(target), 'result_sha256': digest_file(target)}})
    task = {'outputs': [{'path': str(marker), 'kind': 'mrlink_marker'}]}
    assert validate_outputs(task)
    target.write_text('corrupt')
    with pytest.raises(ValueError):
        validate_outputs(task)

def test_worker_affinity_and_thread_environment(tmp_path):
    output = tmp_path / 'settings.json'
    script = "import os,json,pathlib;pathlib.Path(%r).write_text(json.dumps({'omp':os.environ['OMP_NUM_THREADS'],'blas':os.environ['OPENBLAS_NUM_THREADS'],'affinity':list(os.sched_getaffinity(0))}))" % str(output)
    task = dict(task_id='resource-check', argv=[sys.executable, '-c', script], outputs=[{'path': str(output), 'kind': 'json'}], cores=1, memory_gib=0.01)
    r = execute([task], tmp_path / 'executor', 2, 1, 1, tmp_path / 'pool', 0)
    assert r[0]['status'] == 'completed'
    d = json.loads(output.read_text())
    assert d['omp'] == d['blas'] == '1' and len(d['affinity']) == 1

def test_geometry_query_batch_and_order_invariance():
    b = bank()
    z = np.random.default_rng(3).normal(size=(len(b.ids), 3))
    geo = m.S7RetrievalGeometry().fit(b, z, cfg())
    q = b.take([0, 4])
    q.ids[:] = ['x', 'y']
    q.groups[:] = ['xg', 'yg']
    one = geo.candidates(q, z[[0, 4]])
    two = geo.candidates(q.take([1, 0]), z[[4, 0]])
    np.testing.assert_array_equal(one[0], two[0][::-1])
    np.testing.assert_allclose(one[1], two[1][::-1])


def test_real_cuda_3000_assays_and_peak(tmp_path):
    torch=m.torch
    if not torch.cuda.is_available():pytest.fail('CUDA is required on this acceptance host')
    torch.cuda.reset_peak_memory_stats()
    start=time.monotonic();m.s7_model_preflight(config('--device','cuda'),features=3000)
    record=dict(features=3000,actual_device='cuda',forward_backward=True,wall_seconds=time.monotonic()-start,
        peak_allocated_bytes=torch.cuda.max_memory_allocated(),peak_reserved_bytes=torch.cuda.max_memory_reserved(),gpu=torch.cuda.get_device_name())
    (ROOT.parent/'gpu-3000-evidence.json').write_text(json.dumps(record,indent=2))


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


def test_policy_technical_abstention_and_regular_fallback(fitted):
    b=fitted['b'];p=fitted['p'].iloc[:5].copy();x=fitted['x'][:5].copy();x[0,:]=np.nan
    table=m.s7_predict(b,x,p)
    assert not table.policy_prediction_available.iloc[0] and np.isnan(table.policy_risk.iloc[0])
    assert table.policy_source.iloc[0]=='abstain_technical_QC'
    ordinary=(~table.released)&table.policy_prediction_available
    assert ordinary.any()
    np.testing.assert_allclose(table.loc[ordinary,'policy_risk'],table.loc[ordinary,b['fallback']])
    metric=m.s7_metric(p.time.to_numpy(),p.event.to_numpy(),table.policy_risk.to_numpy(),p.family.to_numpy(),b['censoring'])
    assert metric['eligible_N']==5 and metric['abstained_N']>=1 and metric['policy_available_N']<5


def test_calibration_uses_frozen_cuts_and_preserves_artifact(fitted):
    b=fitted['b'];out=fitted['out'];before=hashlib.sha256((out/'model_bundle.joblib').read_bytes()).hexdigest()
    m.s7_evaluate(out)
    assert before==hashlib.sha256((out/'model_bundle.joblib').read_bytes()).hexdigest()
    d=pd.read_csv(out/'test_calibration.csv');g=pd.read_csv(out/'test_calibration_groups.csv')
    assert set(d.cutpoint_source)=={'calibration_fit'} and (d.policy_available_N+d.abstained_N==d.eligible_N).all()
    assert set(['all','research','released'])<=set(d.subset)
    for name,rows in g.groupby('model_id'):
        expected=set(np.asarray(b['calibration_cutpoints'][name]).tolist())
        for edge in rows.upper[np.isfinite(rows.upper)]:assert np.min(np.abs(np.asarray(list(expected))-edge))<1e-14


def test_eigen_cache_full_allele_identity(tmp_path,monkeypatch):
    monkeypatch.setenv('MRLINK2_EIGH_CACHE',str(tmp_path))
    x=np.random.default_rng(4).normal(size=(15,6));a=x.T@x;alleles=[[f's{i}','A','C'] for i in range(6)]
    l,u,d=r.cached_eigh(a,alleles);l2,u2,d2=r.cached_eigh(a,alleles)
    assert not d['cache_hit'] and d2['cache_hit'];np.testing.assert_array_equal(l,l2)
    changed=[v.copy() for v in alleles];changed[0][1:]=['C','A']
    assert r.cached_eigh(a,changed)[2]['cache_key']!=d['cache_key']
    assert r.cached_eigh(a[::-1,::-1],alleles[::-1])[2]['cache_key']!=d['cache_key']
    assert r.cached_eigh(a[:3,:3],alleles[:3])[2]['cache_key']!=d['cache_key']


def test_mr_real_N_and_outcome_only_initialization(tmp_path):
    n=mr.numeric_module();beta=np.array([.1,.2,.3]);ii=np.array([0,2]);value=n.le8_instrument_variance(beta,3000,ii)
    assert value==pytest.approx(np.sum(beta[ii]**2-1/3000))
    source=(ROOT/'f/c2.mr_link2.py').read_text()
    assert 'le8_instrument_variance(out_betas, n_out, outcome_instruments)' in source
    frame=pd.DataFrame(dict(SNP=['rs1'],CHR=[1],POS=[100],EA=['A'],NEA=['C'],BETA=[.1],SE=[.01],P=[.001],EAF=[.2]))
    file=tmp_path/'in.tsv';frame.to_csv(file,sep='\t',index=False)
    with pytest.raises((ValueError,SystemExit),match='missing_N'):
        n.le8_prepare_sumstats(str(file),str(tmp_path/'bad.tsv.gz'),float('nan'),'1:1-200')
    n.le8_prepare_sumstats(str(file),str(tmp_path/'ok.tsv.gz'),12345,'1:1-200')
    assert (pd.read_csv(tmp_path/'ok.tsv.gz',sep='\t').N==12345).all()


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


def test_frozen_cpu_cuda_probabilities_gate_release(fitted):
    import copy
    b=copy.deepcopy(fitted['b']);x=fitted['x'][:32];p=fitted['p'].iloc[:32]
    cpu=m.s7_predict(b,x,p)
    for learner in b['models'].values():
        if isinstance(learner,m.S7AttentionLearner):
            learner.device='cuda';learner.config['s7_retrieval_device']='cuda'
    gpu=m.s7_predict(b,x,p)
    for name in ('policy_risk',b['primary'],'gain_score'):
        np.testing.assert_allclose(cpu[name],gpu[name],atol=2e-6,rtol=0,equal_nan=True)
    for name in ('research_mask','released','technical_QC_pass','policy_prediction_available'):
        np.testing.assert_array_equal(cpu[name],gpu[name])


def test_prediction_preserves_original_cuda_error(fitted,monkeypatch):
    import copy
    learner=copy.deepcopy(fitted['b']['models'][fitted['b']['primary']]);learner.device='cuda'
    def fail(*args,**kwargs):raise RuntimeError('intentional CUDA kernel error')
    monkeypatch.setattr(learner,'predict_batch',fail)
    with pytest.raises(RuntimeError,match='intentional CUDA kernel error'):
        learner.predict(fitted['x'][:2],fitted['p'].iloc[:2])
    assert learner.device=='cuda' and next(learner.network.parameters()).device.type=='cpu'


def test_mr_raw_batch_extraction_and_task_identity(tmp_path,monkeypatch):
    import gzip
    source=tmp_path/'raw.gz'
    with gzip.open(source,'wt') as h:h.write('SNP\tCHR\tPOS\tN\nrs1\t1.0\t100\t300\nrs2\t1\t150\t300\nrs3\t1\t250\t300\n')
    ref=tmp_path/'ref';ref.write_text('identity fixture')
    worker=tmp_path/'c2.cause.sh';worker.write_text('#!/bin/bash\nexit 0\n')
    jobs=tmp_path/'jobs.tsv'
    jobs.write_text('omics\ttrait\texposure\toutcome\tregion\n'+''.join(f'protein\tP\t{source}\t{source}\t{v}\n' for v in ['1:90-200','1:200-300']))
    a=NS(jobs=jobs,worker_script=worker,cad_gwas=str(source),cache_root=tmp_path/'cache',hash_mode='sha256',identity_file=[ref],job_memory_gib=.1,job_base_gib=.01,ld_workspace_multiplier=8,inner_threads=1)
    tasks=mr.build_tasks(a,[]);assert len(set(t['task_id'] for t in tasks))==2
    one=pd.read_csv(next((tmp_path/'cache/_regions').rglob('*.csv.gz')))
    assert len(one) in (1,2)
    a.inner_threads=2;monkeypatch.setenv('MRLINK2_WORKERS','4');other=mr.build_tasks(a,[])
    assert [t['task_id'] for t in tasks]==[t['task_id'] for t in other]
    ref.write_text('changed reference');other=mr.build_tasks(a,[])
    assert [t['task_id'] for t in tasks]!=[t['task_id'] for t in other]
    monkeypatch.setenv('MRLINK2_EXPOSURE_N','500');changed=mr.build_tasks(a,[])
    assert [t['task_id'] for t in other]!=[t['task_id'] for t in changed]


@pytest.mark.parametrize('status',['failed','not_run','no_estimate'])
def test_actual_adapter_validates_exit_zero_status(tmp_path,status):
    out=tmp_path/'worker-results';marker=tmp_path/'marker.json';taskfile=tmp_path/'task.json'
    code="from pathlib import Path;p=Path(%r);(p/'mrlink2.status.tsv').write_text(%r)"%(str(out),'omics\ttrait\tstatus\tmessage\nprotein\tP\t'+status+'\tno_shared_snps\n')
    task=dict(task_id='test',omics='protein',trait='P',scientific_signature='fixture',worker_outdir=str(out),argv=[sys.executable,'-c',code])
    atomic_json(taskfile,task)
    z=subprocess.run([sys.executable,str(ROOT/'f/c2.parallel.py'),'worker','--task',str(taskfile),'--marker',str(marker)],capture_output=True,text=True)
    assert (z.returncode==0)==(status=='no_estimate')
    assert marker.exists()==(status=='no_estimate')
    if marker.exists():assert json.loads(marker.read_text())['row']['message']=='no_shared_snps'


def test_multi_module_native_options_do_not_damage_shared_values(tmp_path,monkeypatch):
    spec=importlib.util.spec_from_file_location('final_dispatch',ROOT/'f/0.common.py');d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
    calls=[];monkeypatch.setattr(d,'dispatch_main',lambda a:calls.append(a));monkeypatch.setattr(d,'_table_runtime',lambda a:False)
    d.main(['dispatch','c1_correlate,c1_abm','--native-test-option','cvd_cad','--Y','cvd_cad','--biom','prot','--analysis-root',str(tmp_path)])
    abm=[a for a in calls if a[0]=='c1_abm']
    assert len(abm)==2
    for a in abm:assert '--native-test-option' not in a and a[a.index('--Y')+1]=='cvd_cad'
    assert '--native-test-option' in next(a for a in calls if a[0]=='c1_correlate')


def test_force_replace_rebuilds_mr_attempt(tmp_path,monkeypatch):
    source=tmp_path/'input';source.write_text('SNP\tCHR\tPOS\tN\nrs1\t1\t100\t300\n')
    ref=tmp_path/'ref';ref.write_text('identity fixture');worker=tmp_path/'worker.sh';worker.write_text('exit 0\n')
    jobs=tmp_path/'jobs.tsv';jobs.write_text(f'omics\ttrait\texposure\toutcome\tregion\nprotein\tP\t{source}\t{source}\t1:90-200\n')
    a=NS(jobs=jobs,worker_script=worker,cad_gwas=str(source),cache_root=tmp_path/'cache',hash_mode='sha256',identity_file=[ref],job_memory_gib=.1,job_base_gib=.01,ld_workspace_multiplier=8,inner_threads=1)
    first=mr.build_tasks(a,[])[0]
    atomic_json(first['outputs'][0]['path'],dict(status='no_estimate',row=dict(message='no_shared_snps')))
    again=mr.build_tasks(a,[])[0];assert again['argv']==first['argv']
    monkeypatch.setenv('LE8_REPLACE','TRUE');new=mr.build_tasks(a,[])[0]
    assert new['task_id']==first['task_id'] and new['argv']!=first['argv']


def test_shiny_aggregate_reference_and_private_id_guard():
    script = r'''
args <- commandArgs(trailingOnly=TRUE)
e <- new.env(parent=baseenv())
for (expr in parse(args[1])) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c('private_columns', 'has_private_columns')) eval(expr, e)
}
aggregate <- c('run_id','model_id','primary_id','reference_id','n','uno_c_horizon')
stopifnot(!e$has_private_columns(aggregate))
stopifnot(e$has_private_columns(c(aggregate, 'eid')))
stopifnot(e$has_private_columns(c(aggregate, 'sample_id')))
stopifnot(e$has_private_columns(c('reference_id', 'risk')))
stopifnot(!e$has_private_columns(c('feature', 'beta', 'p.value')))
'''
    z=subprocess.run(['Rscript','-e',script,str(ROOT/'shiny/app.R')],capture_output=True,text=True)
    assert z.returncode==0,z.stdout+z.stderr
