"""Real end-to-end test on small simulated PGENs; never touches reference data.

Run with ~/.venvs/gu-threads/bin/python f/tests/test_arg_threads.py
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import msprime
import numpy as np
import pgenlib
import tskit

GU = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(GU/'f'))
from arg_threads import repack


def fixture(root, chrom):
    prefix = root/'pfile'/f'chr{chrom}'
    n = 20
    ts = msprime.sim_ancestry(n, population_size=10000, sequence_length=20000,
                              recombination_rate=1e-8, random_seed=123)
    ts = msprime.sim_mutations(ts, rate=2e-8, random_seed=321, model=msprime.BinaryMutationModel())
    variants = [(int(v.site.position), np.array(v.genotypes, dtype=np.int32, copy=True))
                for v in ts.variants() if 0 < v.genotypes.sum() < 2*n]
    offset = 3000000 if chrom == 'X' else 1000000
    with pgenlib.PgenWriter(os.fsencode(str(prefix)+'.pgen'), n, variant_ct=len(variants),
                           hardcall_phase_present=True) as writer:
        for pos, g in variants:
            if chrom == 'X':
                for male in range(3): g[2*male+1] = g[2*male]
            writer.append_alleles(g, all_phased=True)
    Path(str(prefix)+'.psam').write_text('#IID\tSEX\n'+''.join(
        f's{i}\t{1 if i<3 else 2}\n' for i in range(n)))
    Path(str(prefix)+'.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(
        f'{chrom}\t{offset+pos}\tv{j}\tA\tG\n' for j,(pos,g) in enumerate(variants)))
    (root/'maps'/f'chr{chrom}.map').write_text(
        f'position\tCOMBINED_rate(cM/Mb)\tGenetic_Map(cM)\n{offset}\t1\t0\n{offset+20000}\t1\t0.02\n')


with tempfile.TemporaryDirectory(prefix='gu-threads-e2e-') as tmp:
    root = Path(tmp)
    for name in ('pfile','maps'): (root/name).mkdir()
    (root/'samples.txt').write_text('sample\tsuper_pop\n'+''.join(f's{i}\tEUR\n' for i in range(20)))
    fixture(root,'22'); fixture(root,'X')
    common = ['--dir-gen',str(root),'--method','threads','--format','trace','--chr','22,X',
              '--map-dir',str(root/'maps'),'--threads','2','--jobs','4','--scratch',str(root/'scratch')]
    env = dict(os.environ, GU_THREADS_ENV=str(Path(sys.executable).parent.parent))
    def run(action, *extra, ok=True):
        result = subprocess.run(['bash',str(GU/'arg.sh'),action,*common,*extra], env=env,
                                text=True,capture_output=True,timeout=300)
        if ok and result.returncode:
            print(result.stdout); print(result.stderr); raise AssertionError(result.returncode)
        return result
    run('prep_gen')
    result = run('build')
    assert 'BUILD PASS' in result.stdout, result.stdout
    assert 'SKIP verified Threads PGEN' in result.stdout
    output = root/'arg.threads'
    for chrom in ('22','X'):
        prep = root/'pfile.4arg.threads'/f'chr{chrom}'
        ts = tskit.load(output/'trace/threads'/f'chr{chrom}.trees')
        assert ts.time_units == 'generations'
        assert min(ts.tables.sites.position) >= (3000000 if chrom=='X' else 1000000)
        n = 36 if chrom=='X' else 40  # three males -> drop one, 17 females + 2 males
        assert ts.num_samples == n, (chrom,ts.num_samples)
        with pgenlib.PgenReader(os.fsencode(str(prep)+'.pgen')) as reader:
            g = np.empty(n,dtype=np.int32)
            assert reader.get_variant_ct() == ts.num_sites
            for i,v in enumerate(ts.variants()):
                reader.read_alleles(i,g)
                observed = np.array([v.alleles[j] for j in v.genotypes])
                expected = np.array(['A','G'])[g]
                assert np.array_equal(observed,expected), (chrom,i,'genotype/identity mismatch')
        log = (output/'log'/f'arg.infer.chr{chrom}.log').read_text()
        assert '--num_threads 2' in log
    run('check')
    precheck = subprocess.run(['bash',str(GU/'f/arg.sh'),'trace_check','check'],
        env=dict(env,GU_TARGET_ROOT=str(root),GU_ARG_DIR=str(output),GU_ARG_METHOD='threads',
                 GU_ARG_FORMAT='trace',GU_CHRS='22 X',GU_BUILD='37'),
        text=True,capture_output=True,timeout=60)
    assert precheck.returncode == 0, (precheck.stdout,precheck.stderr)
    cached = run('build')
    assert 'SKIP verified Threads inference chr22' in cached.stdout
    assert 'SKIP verified Threads inference chrX' in cached.stdout
    # Detect output corruption rather than silently reusing an existing file.
    with (output/'threads/chr22.threads').open('ab') as h: h.write(b'corrupt')
    assert run('check',ok=False).returncode != 0
    assert run('build','--action','convert',ok=False).returncode != 0
    # Unknown sex and odd X panel policy must fail before inference.
    assert run('prep_gen','--x-odd-male','error',ok=False).returncode != 0
    print('PASS: real Threads 2-worker PGEN -> inference -> conversion -> TRACE; '
          'autosomes/X, genotypes, coordinates, identity counts, cache and corruption checks')
