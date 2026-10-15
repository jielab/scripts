"""Verify CLI dispatch and failure propagation using a harmless Rscript stub."""
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
passed = 0


def check(label, condition):
    global passed
    assert condition, label
    passed += 1
    print('PASS:', label)


with tempfile.TemporaryDirectory(prefix='bmonth-cli-') as td:
    work = Path(td)
    stub = work / 'Rscript'
    stub.write_text('''#!/usr/bin/env python3
import json, os, pathlib, subprocess, sys
cg = pathlib.Path('/proc/self/cgroup').read_text().strip().split('::', 1)[1]
base = pathlib.Path('/sys/fs/cgroup' + cg)
limits = {name: (base / name).read_text().strip()
          for name in ('memory.max', 'memory.swap.max', 'memory.oom.group')}
child_cg = subprocess.check_output([sys.executable, '-c',
    "from pathlib import Path; print(Path('/proc/self/cgroup').read_text().strip())"], text=True).strip()
with open(os.environ['BMONTH_TEST_CALLS'], 'a') as f:
    f.write(json.dumps({'args': sys.argv[1:], 'out': os.environ.get('BMONTH_OUTDIR'),
                       'cores': os.environ.get('BMONTH_PHEWAS_CORES'),
                       'replace': os.environ.get('BMONTH_REPLACE'),
                       'cgroup': cg, 'child_cgroup': child_cg,
                       'limits': limits}) + '\\n')
if pathlib.Path(sys.argv[-1]).name == os.environ.get('BMONTH_TEST_FAIL'):
    sys.exit(17)
''')
    stub.chmod(0o755)
    # The stub never writes analysis results; avoid even creating the default
    # production output directory during the no-argument dispatch test.
    mkdir = work / 'mkdir'
    mkdir.write_text('#!/bin/sh\nif [ "$*" = "-p /mnt/d/analysis/bmonth" ]; then exit 0; fi\nexec /bin/mkdir "$@"\n')
    mkdir.chmod(0o755)
    calls = work / 'calls.jsonl'
    env = dict(os.environ, PATH=str(work) + ':' + os.environ['PATH'],
               BMONTH_TEST_CALLS=str(calls), SCRIPT_LOG_DIR=str(work / 'logs'))

    def run(args, **extra):
        calls.write_text('')
        result = subprocess.run([str(ROOT / 'bmonth.sh'), *args], cwd=work,
                                env=dict(env, **extra), capture_output=True, text=True)
        rows = [json.loads(line) for line in calls.read_text().splitlines()]
        return result, rows

    sequence = ['00.check.R', '1.prep.R', '3.phewas.R', '2.assoc.R', '4.med.R', '5.final.R']
    r, rows = run([])
    check('no-argument default sequence', r.returncode == 0 and [Path(x['args'][-1]).name for x in rows] == sequence)
    check('no-argument run has persistent log', any((work / 'logs').glob('*.log')))
    caller_cg = Path('/proc/self/cgroup').read_text().strip().split('::', 1)[1]
    check('analysis isolated from caller cgroup', all(x['cgroup'] != caller_cg and
          Path(x['cgroup']).name.startswith('bmonth-memory-') for x in rows))
    check('all stages share one scope', len({x['cgroup'] for x in rows}) == 1)
    check('kernel enforces 24 GiB RAM, zero swap, whole-group OOM', all(x['limits'] == {
          'memory.max': str(24 * 1024**3), 'memory.swap.max': '0', 'memory.oom.group': '1'} for x in rows))
    check('child workers inherit scope', all(x['child_cgroup'] == '0::' + x['cgroup'] for x in rows))
    check('default PheWAS concurrency is automatic', all(x['cores'] == 'auto' for x in rows))
    for stage, name in [('check', '00.check.R'), ('prep', '1.prep.R'), ('phewas', '3.phewas.R'),
                        ('assoc', '2.assoc.R'), ('med', '4.med.R'), ('final', '5.final.R'), ('plot', '6.replot.R')]:
        r, rows = run([stage, '--out', str(work / 'out with spaces')])
        check(stage + ' dispatch', r.returncode == 0 and len(rows) == 1 and Path(rows[0]['args'][-1]).name == name)
        check(stage + ' quoted output path', rows[0]['out'] == str(work / 'out with spaces'))
    r, rows = run(['--out', str(work / 'all'), '--phewas-cores', '2', '--replace'], BMONTH_TEST_FAIL='3.phewas.R')
    check('failed R stage exits with original status', r.returncode == 17)
    check('failed R stage stops later stages', [Path(x['args'][-1]).name for x in rows] == sequence[:3])
    check('options exported', all(x['cores']=='2' and x['replace']=='TRUE' for x in rows))
    for args in [['bad-stage'], ['--phewas-cores', '0'], ['--phewas-cores', '-1'],
                 ['--phewas-cores', 'bad'], ['--med-boot', '-1'], ['--med-boot', 'x'], ['--out']]:
        r, rows = run(args)
        check('reject invalid arguments: ' + ' '.join(args), r.returncode != 0 and not rows)
    r, rows = run(['--help'])
    check('help does not run R', r.returncode == 0 and not rows)
    r, rows = run(['med', '--med-boot', '0', '--out', str(work / 'zero')])
    check('zero bootstrap accepted', r.returncode == 0 and len(rows) == 1)
    systemctl = work / 'systemctl'
    systemctl.write_text('#!/bin/sh\nexit 1\n')
    systemctl.chmod(0o755)
    r, rows = run(['phewas', '--out', str(work / 'unavailable')])
    check('unavailable memory guard prevents analysis', r.returncode != 0 and not rows and
          'refusing to run uncapped' in r.stderr + r.stdout)
    r, rows = run(['-h'])
    check('help works without systemd', r.returncode == 0 and not rows)
print('ALL', passed, 'CLI CHECKS PASSED')
