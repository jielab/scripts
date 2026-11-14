#!/usr/bin/env python3
"""Resume only unsubmitted workers in the inspected November 13 LD plan."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path('/work/sph-huangj/analysis/ld/logs')
PLAN = ROOT / 'run-20261113T150343501842Z-2358583/plan.json'
STATE = ROOT / 'submission-state.json'
ORIGINAL_JOURNAL = PLAN.parent / 'commands/submission-cad74c97ea39.jsonl'
ACTIVE = {'PEND', 'RUN', 'PSUSP', 'USUSP', 'SSUSP', 'WAIT', 'UNKWN', 'PROV'}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true', help='Submit; default only previews')
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from submit import bsub_process
    with STATE.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = json.loads(STATE.read_text())
        plan = json.loads(PLAN.read_text())
        if state.get('plan') != str(PLAN) or state.get('schema') != 'ld.submission.1':
            raise RuntimeError('The current plan changed; inspect before resuming')
        missing = [j for j in plan['jobs'] if j['key'] not in state['workers']]
        probe = subprocess.run(['bjobs', '-a', '-noheader', '-o', 'jobid stat job_name'],
                               capture_output=True, text=True, timeout=60)
        if probe.returncode and 'No unfinished job found' not in probe.stderr:
            raise RuntimeError(probe.stdout + probe.stderr)
        observed = {}
        names = set()
        for line in probe.stdout.splitlines():
            fields = line.split()
            if len(fields) >= 3 and fields[0].isdigit():
                observed[fields[0]] = fields[1]
                names.add(fields[2])
        for job in missing:
            if 'ld38.' + job['key'] in names:
                raise RuntimeError('Unrecorded LSF job already exists: ' + job['key'])
        uncertain = state.get('uncertain')
        if uncertain:
            expected = str(PLAN.parent / 'commands/AFR.chr20.cmd')
            if uncertain != {'cmd': expected, 'journal': str(ORIGINAL_JOURNAL)}:
                raise RuntimeError('New uncertain acknowledgement; inspect before retrying')
            history = subprocess.run(['bhist', '-n', '1', '-u', 'sph-huangj', '-J',
                                      'ld38.AFR.chr20', '-w'], capture_output=True, text=True, timeout=120)
            if 'No matching job found' not in history.stdout + history.stderr:
                raise RuntimeError('Cannot exclude an earlier AFR.chr20 submission: ' + history.stdout + history.stderr)
            if any((PLAN.parent / 'logs').glob('AFR.chr20.*')):
                raise RuntimeError('AFR.chr20 output exists; inspect before retrying')
        print('Missing workers:', len(missing), ', '.join(j['key'] for j in missing), flush=True)
        print('Merge:', 'already submitted' if state.get('merge') else 'will submit after workers', flush=True)
        if not args.execute:
            print('Preview only; use --execute to submit')
            return
        def save():
            temporary = STATE.with_suffix('.resume.tmp')
            with temporary.open('w') as stream:
                json.dump(state, stream, indent=2)
                stream.flush()
                os.fsync(stream.fileno())
            temporary.replace(STATE)
        journal = PLAN.parent / 'commands/resume-missing-20261114.jsonl'
        if uncertain:
            with journal.open('a') as stream:
                stream.write(json.dumps({'resolved_original_uncertain': uncertain,
                                         'evidence': 'lock acquired; no bjobs, bhist or output match'}) + '\n')
            state.pop('uncertain')
            save()
        dependencies = [i for i in state['workers'].values() if observed.get(i) in ACTIVE]
        tasks = [(j['key'], j['command']) for j in missing]
        if not state.get('merge'):
            tasks.append(('merge', plan['merge']))
        for key, command in tasks:
            argv = ['bsub']
            if key == 'merge' and dependencies:
                argv += ['-w', ' && '.join('ended(' + i + ')' for i in dependencies)]
            state['uncertain'] = {'cmd': command, 'journal': str(journal)}
            save()
            with open(command) as handle:
                result = bsub_process(argv, handle)
            response = result.stdout + result.stderr
            match = re.search(r'Job <(\d+)>', response)
            with journal.open('a') as stream:
                stream.write(json.dumps({'cmd': command, 'returncode': result.returncode,
                                         'response': response}) + '\n')
                stream.flush()
                os.fsync(stream.fileno())
            if result.returncode or not match:
                raise RuntimeError('Submission not confirmed; inspect ' + str(journal))
            job_id = match[1]
            if key == 'merge':
                state['merge'] = job_id
            else:
                state['workers'][key] = job_id
                dependencies.append(job_id)
            state.pop('uncertain')
            save()
            print(key, job_id, flush=True)

if __name__ == '__main__':
    main()
