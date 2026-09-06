"""Real ASMC numerical regression, all three inference stages and resume.

Usage (grid Python): test_argneedle_memory.py --haps chr1.haps.gz --map chr1.map
Writes only a new temporary directory; the full dataset is read-only.
"""
import argparse
import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import arg_needle
import arg_needle_lib

p = argparse.ArgumentParser()
p.add_argument('--haps', required=True)
p.add_argument('--map', required=True)
a = p.parse_args()
root = Path(tempfile.mkdtemp(prefix='argneedle-memory-test-'))
print(f'Test artifacts: {root}', flush=True)
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           PYTHONHASHSEED='20260904')
package = Path(arg_needle.__file__).parent
script = package / 'scripts/infer_args_advanced.py'
runner = Path(__file__).resolve().parents[1] / 'run_argneedle_advanced.py'

# Retain real variants, map positions and haplotypes; ASMC expects >=300 haps.
with gzip.open(a.haps, 'rt') as haps, open(a.map) as maps, \
        gzip.open(root / 'input.haps.gz', 'wt') as hout, \
        open(root / 'input.map', 'w') as mout:
    count = 0
    for line, mapline in zip(haps, maps):
        fields = line.split()
        alleles = fields[5:305]
        if len(set(alleles)) != 2:
            continue
        hout.write(' '.join(fields[:5] + alleles) + '\n')
        mout.write(mapline)
        count += 1
        if count == 384:
            break
    assert count == 384
(root / 'input.sample').write_text('ID_1 ID_2 missing\n0 0 0\n' +
    ''.join(f's{i} s{i} 0\n' for i in range(150)))
bootstrap = ('import runpy,sys,random,numpy as np; random.seed(20260904); '
             'np.random.seed(20260904); script=sys.argv.pop(1); sys.argv[0]=script; '
             'runpy.run_path(script,run_name="__main__")')
for mode in ('sequence', 'array'):
    baseline = root / (mode + '-baseline')
    bounded = root / (mode + '-bounded')
    with open(root / (mode + '.log'), 'w') as log:
        for step in (1, 2, 3):
            cmd = [sys.executable, '-u', '-c', bootstrap, str(script),
                   '--hap_gz', str(root / 'input.haps.gz'), '--map', str(root / 'input.map'),
                   '--out', str(baseline), '--mode', mode, '--step', str(step),
                   '--chromosome', '1', '--normalize', '1', '--trim_num_snps', '0']
            if step == 1:
                cmd += ['--num_sequence_samples' if mode == 'sequence' else '--num_snp_samples', '150']
            subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=300)
        cmd = [sys.executable, str(runner), '--haps', str(root / 'input.haps.gz'),
               '--map', str(root / 'input.map'), '--out', str(bounded),
               '--chr', '1', '--seed-haplotypes', '150', '--mode', mode,
               '--log', str(root / (mode + '-bounded.log'))]
        subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=600)
        for suffix in ('.step1.argn', '.step2.argn', '.argn'):
            left = arg_needle_lib.arg_to_tskit(arg_needle_lib.deserialize_arg(str(baseline) + suffix))
            right = arg_needle_lib.arg_to_tskit(arg_needle_lib.deserialize_arg(str(bounded) + suffix))
            left.tables.assert_equals(right.tables, ignore_provenance=True)
            assert right.num_samples == (150 if suffix == '.step1.argn' else 300)
        result = subprocess.run(cmd, env=env, text=True, capture_output=True,
                                check=True, timeout=60)
        assert result.stdout.count('SKIP verified') == 3, result.stdout
    print(f'PASS {mode}: identical tree tables for steps 1/2/3 (including normalization); '
          '300 haplotypes; verified resume', flush=True)

# Unknown upstream implementations must fail before modifying either function.
sys.path.insert(0, str(runner.parent))
from argneedle_memory import _rewrite, install
try:
    _rewrite(install, [('nonexistent upstream operation', 'replacement')])
except RuntimeError:
    pass
else:
    raise AssertionError('Unknown source accepted')
subprocess.run([sys.executable, '-c',
    f'import sys; sys.path.insert(0, {str(runner.parent)!r}); '
    'from argneedle_memory import install; install(); install()'], env=env, check=True)
print('PASS compatibility rejection and idempotent install', flush=True)
