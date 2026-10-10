"""C2 MR: numerical contracts, worker results and resumable checkpoints.

Regional numerical validation: python validation/c2_mr.py --regional PATH
"""
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
import pytest
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace as NS

ROOT = Path(__file__).resolve().parents[1]

def load_module(name, filename):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / 'f' / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

ADAPTER = ROOT / 'f/c2.parallel.py'
mr = load_module('le8_validation_mr', 'c2.parallel.py')
validate_status = mr.validate_status
r = load_module('le8_validation_resources', '0.resources.py')
atomic_json = r.atomic_json

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

def test_mr_boundary_diagnostics_do_not_disable_primary_validation(tmp_path):
    n=mr.numeric_module()
    result={'alpha':np.float64(.1),'se(alpha)':.01,'p(alpha)':.2,
            'se(sigma_y)':np.float64('nan'),'optim_ha_success':np.bool_(True)}
    audit=n.le8_optimization_audit(result)
    atomic_json(tmp_path/'audit.json',audit)
    saved=json.loads((tmp_path/'audit.json').read_text())
    assert saved['result']['se(sigma_y)'] is None
    assert saved['nonfinite_result_fields']==['se(sigma_y)']
    assert saved['result']['alpha']==.1 and saved['result']['optim_ha_success'] is True
    assert np.isnan(result['se(sigma_y)'])  # Fitted values were not changed.
    fit=tmp_path/'fit.tsv';status=tmp_path/'status.tsv'
    fit.write_text('alpha\tse(alpha)\tp(alpha)\tse(sigma_y)\n0.1\t0.01\t0.2\tnan\n')
    status.write_text(f'omics\ttrait\tstatus\tout_prefix\nprotein\tP\tok\t{fit}\n')
    assert validate_status(status,'protein','P')['status']=='ok'
    fit.write_text('alpha\tse(alpha)\tp(alpha)\n0.1\tnan\t0.2\n')
    with pytest.raises(ValueError,match='alpha inference'):
        validate_status(status,'protein','P')

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

def regional_numerics(directory):
    m = load_module("le8_validation_mr_numeric", "c2.mr_link2.py")
    w = Path(directory)
    rows = []

    for file in sorted((w/'serial-cache').rglob('*_numerical_audit.json')):
     audit=json.loads(file.read_text());trait=file.name.removeprefix('protein_').removesuffix('_numerical_audit.json')
     if not audit.get('optimizations'):continue
     matrices=np.load(w/'eigen'/audit['eigh']['cache_key']/'eigen.npz');lam=matrices['values'];u=matrices['vectors'];ld=(u*lam)@u.T
     def betas(kind):
      f=file.parent.parent/'prepared'/f'protein_{trait}.{kind}.ref.tsv';data=pd.read_csv(f,sep='\t').set_index('SNP');answer=[]
      for snp,a1,a2 in audit['ordered_alleles']:
       v=data.loc[snp];z=v.BETA/v.SE;b=z/np.sqrt(v.N+z*z);pair=(v.EA,v.NEA);comp=tuple(v.translate(str.maketrans('ACGT','TGCA')) for v in pair)
       if pair==(a1,a2) or comp==(a1,a2):answer.append(b)
       elif pair==(a2,a1) or comp==(a2,a1):answer.append(-b)
       else:raise ValueError('Unexpected allele mismatch')
      return np.array(answer)
     x,y=betas('exposure'),betas('outcome');nX,nY=audit['N_exposure'],audit['N_outcome'];guess=audit['sigma_initial']
     threshold=sorted(lam[lam>0],reverse=True)[np.argmin(np.cumsum(sorted(lam[lam>0],reverse=True))/sum(lam[lam>0])<=.99)]
     take=lam>=threshold
     for factor in [1.,.5,2.]:
      result=m.mr_link2(lam[take],u[:,take],x,y,nX,nY,guess[0]*factor,guess[1]*factor)
      rows.append(dict(trait=trait,backend='cpu',initial_factor=factor,**{k:(v.item() if hasattr(v,'item') else v) for k,v in result.items()}))
      if factor==1:np.testing.assert_allclose(result['alpha'],audit['optimizations'][0]['result']['alpha'],atol=1e-6,rtol=0)
     start=time.monotonic();gl,gu,diag=r.validated_eigh(ld,'cuda');seconds=time.monotonic()-start
     gthreshold=sorted(gl[gl>0],reverse=True)[np.argmin(np.cumsum(sorted(gl[gl>0],reverse=True))/sum(gl[gl>0])<=.99)];gt=gl>=gthreshold
     result=m.mr_link2(gl[gt],gu[:,gt],x,y,nX,nY,*guess)
     rows.append(dict(trait=trait,backend='cuda_eigh',eigh_seconds=seconds,initial_factor=1.,eigen_diagnostic=diag,
       **{k:(v.item() if hasattr(v,'item') else v) for k,v in result.items()}))
     np.testing.assert_allclose(result['alpha'],audit['optimizations'][0]['result']['alpha'],atol=1e-6,rtol=0)
    (w/'numerical-starts-cuda.json').write_text(json.dumps(rows,indent=2))
    print('PASS',len(rows),'actual regional fits; 3 CPU starts + CUDA float64 eigen backend per trait')
    for trait in sorted({v['trait'] for v in rows}):
     z=[v for v in rows if v['trait']==trait];print(trait,'alpha_range',max(v['alpha'] for v in z)-min(v['alpha'] for v in z),'ha_loglik_range',max(v['ha_loglik'] for v in z)-min(v['ha_loglik'] for v in z),'all_converged',all(v[k] for v in z for k in ['optim_ha_success','optim_alpha_h0_success','optim_sigma_y_h0_success']))

if __name__ == '__main__':
    if sys.argv[1:2] == ['--regional'] and len(sys.argv) == 3:
        regional_numerics(sys.argv[2])
    else:
        raise SystemExit(pytest.main([__file__, *sys.argv[1:]]))
