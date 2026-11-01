"""Regression checks for ARG identity, TRACE provenance, and safe output handling.

Run with: python -m pytest -q gu/tests/test_arg_trace.py
Scientific integration checks skip if their optional dependencies are absent.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import gzip
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import numpy as np
import pytest

GU = Path(__file__).resolve().parents[2]


def load(name, path):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


trace = load('gu_trace_test', GU / 'f/trace.py')
arg = load('gu_arg_test', GU / 'f/arg.py')


def test_map_conversion_three_columns_gzip(tmp_path, monkeypatch):
    src, out = tmp_path/'map.gz', tmp_path/'map.txt'
    with gzip.open(src, 'wt') as handle:
        handle.write('position rate genetic_map\n1 1.2 0\n1000 1.1 0.0012\n')
    monkeypatch.setattr(sys, 'argv', ['normalize-map', '--input', str(src), '--output', str(out), '--chr', 'X'])
    trace.normalize_map_cli()
    assert out.read_text().splitlines()[1] == 'chrX\t1\t1.2\t0'
    before = out.stat().st_mtime_ns
    trace.normalize_map_cli()
    assert out.stat().st_mtime_ns == before


def test_map_wrong_chromosome_rejected(tmp_path, monkeypatch):
    src = tmp_path/'chr10.txt'
    src.write_text('chr10 1 1 0\nchr10 1000 1 0.001\n')
    monkeypatch.setattr(sys, 'argv', ['normalize-map', '--input', str(src), '--output', str(tmp_path/'out'), '--chr', '1'])
    with pytest.raises(SystemExit, match='chromosome does not match'):
        trace.normalize_map_cli()


def test_map_nonmonotonic_rejected(tmp_path, monkeypatch):
    src = tmp_path/'map.txt'
    src.write_text('1 1 0\n100 1 0.5\n200 1 0.3\n')
    monkeypatch.setattr(sys, 'argv', ['normalize-map', '--input', str(src), '--output', str(tmp_path/'out'), '--chr', '1'])
    with pytest.raises(SystemExit, match='nondecreasing'):
        trace.normalize_map_cli()


def test_nonfinite_posterior_is_not_reused(tmp_path):
    path = tmp_path/'bad.chr1.xss.npz'
    np.savez(path, gammas=np.array([[np.nan], [1.]]), params=np.ones((2,9)))
    assert not trace.valid_npz(path)
    path2 = tmp_path/'good.chr1.xss.npz'
    np.savez(path2, gammas=np.array([[.1, .2], [.9, .8]]), params=np.ones((2,9)))
    assert trace.valid_npz(path2)


def test_summary_validation_accepts_zero_calls(tmp_path):
    path=tmp_path/'hap0.summary.txt'
    path.write_text('chromosome\tstart\tend\tmean_posterior\tlength(bp)\tlength(cM)\n')
    assert trace.valid_output(path)
    path.write_text('unrelated\theader\n')
    assert not trace.valid_output(path)


def test_arg_map_rejects_duplicate_biological_haplotypes():
    msprime = pytest.importorskip('msprime')
    ts=msprime.sim_ancestry(2, sequence_length=10000, recombination_rate=1e-7, population_size=10000, random_seed=1)
    ts=msprime.sim_mutations(ts, rate=1e-7, random_seed=2)
    rows=[(int(node), 'person', 1) for node in ts.samples()]
    with pytest.raises(SystemExit, match='duplicate/empty biological'):
        arg.trace_validate(ts, rows, '1')


def test_pinned_trace_initializer_is_used_directly(monkeypatch):
    upstream=pytest.importorskip('tracehmm.trace')
    cli=pytest.importorskip('tracehmm.infer_cli')
    original=upstream.TRACE.init_ncoal_gamma_params
    called=[]
    monkeypatch.setenv('TRACE_BACKEND','cpu')
    monkeypatch.setattr(cli,'main',lambda: called.append(True))
    monkeypatch.setattr(sys,'argv',['trace.py'])
    trace.trace_infer_main()
    assert called == [True]
    assert upstream.TRACE.init_ncoal_gamma_params is original


def test_summary_does_not_modify_inference_checkpoint(tmp_path):
    pytest.importorskip('tracehmm')
    binary = shutil.which('trace-summarize')
    if binary is None:
        pytest.skip('trace-summarize CLI is not available')
    xss=tmp_path/'hap0.chr1.xss.npz'
    spans=np.array([[0.,100000.],[100000.,200000.],[200000.,300000.]])
    np.savez_compressed(xss, gammas=np.array([[.01,.01,.01],[.99,.99,.99]]),
                        params=np.ones((2,9)), treespan=spans/1e6, treespan_phy=spans)
    original=xss.read_bytes()
    rc=trace.run(tmp_path/'hap0', ['.summary.txt'], [binary,'-f',str(xss),'-c','chr1'])
    assert rc==0
    assert xss.read_bytes()==original
    assert trace.valid_output(tmp_path/'hap0.summary.txt')


def make_arg(root, chrom, *, offset=0, haploid=False, seed=10):
    msprime=pytest.importorskip('msprime')
    tskit=pytest.importorskip('tskit')
    ts=msprime.sim_ancestry(8, ploidy=1 if haploid else 2, sequence_length=200000,
                          recombination_rate=2e-8, population_size=10000, random_seed=seed)
    ts=msprime.sim_mutations(ts, rate=2e-8, random_seed=seed+1)
    tables=ts.dump_tables()
    if offset:
        tables.sequence_length+=offset
        tables.edges.left+=offset
        tables.edges.right+=offset
        tables.sites.position+=offset
    tables.individuals.metadata_schema=tskit.MetadataSchema.permissive_json()
    inds=tables.individuals.copy()
    tables.individuals.clear()
    for individual_index, ind in enumerate(inds):
        tables.individuals.add_row(flags=ind.flags, location=ind.location, parents=ind.parents,
                                  metadata={'sample':f'sample{individual_index}'})
    ts=tables.tree_sequence()
    path=root/f'chr{chrom}.trees'
    ts.dump(path)
    rows=['tree_node_id\tsample\thaplotype']
    counts={}
    for node in ts.samples():
        sid=ts.node(node).individual
        counts[sid]=counts.get(sid,0)+1
        rows.append(f'{node}\tsample{sid}\t{counts[sid]}')
    (root/f'chr{chrom}.sample_map.tsv').write_text('\n'.join(rows)+'\n')
    return path


def test_trace_check_provenance_rejects_changed_arg(tmp_path):
    pytest.importorskip('msprime')
    trees=tmp_path/'trees'; trees.mkdir()
    path=make_arg(trees,'1')
    (trees/'ARG_TRACE_BUILD.tsv').write_text('method\ttsinfer\nformat\ttrace\nbuild\tGRCh37\n')
    env={**os.environ,'GU_ARG_DIR':str(trees),'TRACE_OUT':str(tmp_path/'out'),'TRACE_CHRS':'1',
         'TRACE_MIN_FULL_ARG_SITES':'2','TRACE_MIN_FULL_ARG_SPAN_FRACTION':'0.1'}
    command=['bash',str(GU/'f/trace.sh'),'check']
    first=subprocess.run(command,env=env,text=True,capture_output=True)
    assert first.returncode==0,first.stdout+first.stderr
    second=subprocess.run(command,env=env,text=True,capture_output=True)
    assert second.returncode==0,second.stdout+second.stderr
    stamp=path.stat().st_mtime_ns
    os.utime(path,ns=(stamp,stamp+1000000000))
    changed=subprocess.run(command,env=env,text=True,capture_output=True)
    assert changed.returncode!=0
    assert 'ARG inputs' in changed.stderr


def test_trace_x_and_autosome_maps_checked_separately(tmp_path):
    pytest.importorskip('msprime')
    trees=tmp_path/'trees';trees.mkdir()
    make_arg(trees,'1')
    make_arg(trees,'X',offset=3000000,haploid=True)
    env={**os.environ,'GU_ARG_DIR':str(trees),'TRACE_OUT':str(tmp_path/'out'),'TRACE_CHRS':'1 X',
         'TRACE_MIN_FULL_ARG_SITES':'2','TRACE_MIN_FULL_ARG_SPAN_FRACTION':'0.01'}
    result=subprocess.run(['bash',str(GU/'f/trace.sh'),'check'],env=env,text=True,capture_output=True)
    assert result.returncode==0,result.stdout+result.stderr
    auto=tmp_path/'out/chromosome_classes/autosomes/samples/trace_sample_map.tsv'
    x=tmp_path/'out/chromosome_classes/chrX/samples/trace_sample_map.tsv'
    assert len(auto.read_text().splitlines())==17
    assert len(x.read_text().splitlines())==9
    assert 'chromosomes\t1 X' in (tmp_path/'out/run.meta.tsv').read_text()


def test_real_trace_joint_autosomes_and_haploid_x(tmp_path):
    pytest.importorskip('msprime')
    pytest.importorskip('tracehmm')
    if not all(shutil.which(tool) for tool in ('trace-extract','trace-infer','trace-summarize','bedtools')):
        pytest.skip('TRACE CLIs and bedtools are required for the synthetic integration check')
    trees=tmp_path/'trees';trees.mkdir()
    make_arg(trees,'1',seed=13)
    make_arg(trees,'2',seed=17)
    make_arg(trees,'X',offset=3000000,haploid=True,seed=21)
    selected=tmp_path/'selected.txt'; selected.write_text('sample0\n')
    out=tmp_path/'out'
    env={**os.environ,'GU_ARG_DIR':str(trees),'TRACE_OUT':str(out),'TRACE_CHRS':'1 2 X',
         'TRACE_MIN_FULL_ARG_SITES':'2','TRACE_MIN_FULL_ARG_SPAN_FRACTION':'0.01',
         'TRACE_TARGET_SAMPLES':str(selected),'TRACE_JOB_EXTRACT':'1','TRACE_JOB_INFER':'1','TRACE_JOB_SUMMARIZE':'1'}
    command=['bash',str(GU/'f/trace.sh'),'run']
    result=subprocess.run(command,env=env,text=True,capture_output=True,timeout=180)
    (tmp_path/'integration.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stdout[-3000:]+result.stderr[-3000:]
    assert (out/'final/trace_haplotype_segments.tsv.gz').is_file()
    inference=sorted(out.glob('chromosome_classes/*/infer/*.xss.npz'))
    assert len(inference)==5
    assert all(trace.valid_npz(f) for f in inference)
    auto=out/'chromosome_classes/autosomes/infer'
    with np.load(auto/'hap0.chr1.xss.npz') as a, np.load(auto/'hap0.chr2.xss.npz') as b:
        np.testing.assert_array_equal(a['params'],b['params'])
        assert 'states' not in a.files
    stamps={str(f):f.stat().st_mtime_ns for f in inference}
    result=subprocess.run(command,env=env,text=True,capture_output=True,timeout=90)
    assert result.returncode==0,result.stdout[-2000:]+result.stderr[-2000:]
    assert stamps=={str(f):f.stat().st_mtime_ns for f in inference}


def test_arg_tsinfer_rejects_virtual_limit_that_breaks_lmdb(tmp_path):
    result=subprocess.run(['bash',str(GU/'arg.sh'),'build','--method','tsinfer','--resource-guard','ulimit','--dir-gen',str(tmp_path)],text=True,capture_output=True)
    assert result.returncode!=0
    assert '1 TiB LMDB' in result.stderr


def test_arg_scheduler_guard_requires_lsf_allocation(tmp_path):
    env={k:v for k,v in os.environ.items() if k!='LSB_JOBID'}
    result=subprocess.run(['bash',str(GU/'arg.sh'),'build','--method','threads','--resource-guard','scheduler','--dir-gen',str(tmp_path)],env=env,text=True,capture_output=True)
    assert result.returncode!=0
    assert 'LSB_JOBID' in result.stderr
