#!/usr/bin/env bash
# Generate a standalone .cmd, then submit to LSF or run locally.
scripts_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
source "$scripts_root/0f/0os.sh"
set -euo pipefail
exec "$SHARED_PYTHON" - "${BASH_SOURCE[0]}" "$@" <<'PY'
import argparse
import os
from pathlib import Path
import shlex
import sys

script = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(script.parents[2] / '0f'))
from workflow_mode import resolve_mode
from qiming_resources import cpu_slots
import submit_jobs

p = argparse.ArgumentParser(prog='concat_ygb_ibdmix.sh', description='Concatenate ygb001–ygb049 IBDmix native calls; generate .cmd and run locally or submit with bsub.')
p.add_argument('--dir-in', required=True, type=Path)
p.add_argument('--dir-out', required=True, type=Path)
p.add_argument('--submit', choices=['auto', 'hpc', 'local'], default=None)
p.add_argument('--queue', default='33')
p.add_argument('--walltime', default='144:00')
p.add_argument('--plan-only', action='store_true', help='Write .cmd without running or submitting it')
a = p.parse_args(sys.argv[2:])
try:
    mode = resolve_mode(a.submit)
    src, dst = a.dir_in.expanduser().resolve(), a.dir_out.expanduser().resolve()
    if any(c in str(dst) for c in '\n\r"'):
        raise ValueError('Output path cannot contain newlines or double quotes')
    names = {f.name for f in (src / 'ygb001/raw').glob('chr*.txt.gz')}
    if not names:
        raise ValueError(f'No chr*.txt.gz files in {src}/ygb001/raw')
    for i in range(1, 50):
        folder = src / f'ygb{i:03d}' / 'raw'
        found = {f.name for f in folder.glob('chr*.txt.gz')}
        if found != names:
            raise ValueError(f'{folder}: missing={sorted(names-found)}, extra={sorted(found-names)}')
        for name in names:
            if not (folder / name).is_file() or not (folder / name).stat().st_size:
                raise ValueError(f'Missing or empty file: {folder / name}')
        if dst == folder.resolve():
            raise ValueError('Output directory must differ from input raw directories')
    hours, minutes = map(int, a.walltime.split(':'))
    if hours < 0 or not 0 <= minutes < 60:
        raise ValueError('Invalid --walltime; use HHH:MM')
    queue, slots = cpu_slots(a.queue, 1, 4, hours + minutes / 60) if mode == 'hpc' else (a.queue, 1)
except (ValueError, OSError) as e:
    p.error(str(e))

def chrom_key(name):
    chrom = name[3:].split('.')[0]
    return (0, int(chrom)) if chrom.isdigit() else (1, chrom)

jobdir = dst / 'jobs'
jobdir.mkdir(parents=True, exist_ok=True)
run = jobdir / 'concat_ygb_ibdmix'
run.mkdir(exist_ok=True)
cmd = run / 'concat_ygb_ibdmix.cmd'
lines = ['#!/bin/bash -l', f'#BSUB -q {queue}', f'#BSUB -n {slots}',
         '#BSUB -J gu.concat.ygb.ibdmix',
         '#BSUB -R "span[hosts=1] rusage[mem=4GB/host]"',
         f'#BSUB -o "{run}/concat.%J.out"', f'#BSUB -e "{run}/concat.%J.err"',
         'set -euo pipefail']
if mode == 'hpc':
    lines += [': "${LSB_JOBID:?Submit this .cmd with bsub}"']
lines += [f'src={shlex.quote(str(src))}', f'dst={shlex.quote(str(dst))}',
          'names=(' + ' '.join(shlex.quote(n) for n in sorted(names, key=chrom_key)) + ')']
lines += [r'''
mkdir -p "$dst"
# Prevent two concurrent merges from replacing the same output files.
exec 9>"$dst/.concat_ygb_ibdmix.lock"
flock -n 9 || { echo 'ERROR: another merge holds the output lock' >&2; exit 1; }
header=''
# Check all expected inputs again on the execution node before writing output.
for name in "${names[@]}"; do
    for i in $(seq -w 001 049); do
        f="$src/ygb$i/raw/$name"
        [[ -s "$f" ]] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
    done
done
tmp=''
trap '[[ -z $tmp ]] || rm -f -- "$tmp"' EXIT
for name in "${names[@]}"; do
    tmp=$(mktemp "$dst/.${name}.XXXXXX")
    header=$(gzip -cd -- "$src/ygb001/raw/$name" | awk 'NR==1 {print} {next}')
    echo "[START] $name"
    {
        for i in $(seq -w 001 049); do
            f="$src/ygb$i/raw/$name"
            gzip -cd -- "$f" | awk -v expected="$header" -v first="$i" -v input="$f" '
                NR == 1 {
                    if ($0 != expected) {
                        print "ERROR: header mismatch: " input > "/dev/stderr"
                        exit 1
                    }
                    if (first == "001") print
                    next
                }
                { print }
                END {
                    if (NR == 0) {
                        print "ERROR: empty TSV: " input > "/dev/stderr"
                        exit 1
                    }
                }'
        done
    } | gzip -c > "$tmp"
    gzip -t "$tmp"
    mv -- "$tmp" "$dst/$name"
    tmp=''
    echo "[DONE] $dst/$name"
done
''']
cmd.write_text('\n'.join(lines) + '\n')
cmd.chmod(0o755)
print(f'Generated: {cmd}', flush=True)
print(f'Mode: {mode}; chromosomes: {len(names)}; batches: 49', flush=True)
if a.plan_only:
    print('Plan only: no merge or submission performed.')
else:
    if mode == 'local':
        submit_jobs.local_one = submit_jobs.local_foreground_one
    submit_jobs.run([cmd], mode=mode, wait=False)
PY
