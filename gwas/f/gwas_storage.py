"""Disposable SAIGE dosage exports shared within one invocation.

Association receipts bind durable source inputs, never disposable VCF mtimes.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def fingerprint(paths):
    return [(str(p), Path(p).stat().st_size, Path(p).stat().st_mtime_ns)
            for p in sorted(set(paths))]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def valid(path, signature, outputs):
    try:
        old = json.loads(path.read_text())
        return (old['signature'] == signature and old['outputs'] ==
                [list(x) for x in fingerprint(outputs)] and
                all(Path(x).is_file() and Path(x).stat().st_size for x in outputs))
    except (OSError, ValueError, KeyError):
        return False


def record(path, signature, outputs):
    data = dict(signature=signature, outputs=fingerprint(outputs))
    temp = path.with_suffix('.next')
    temp.write_text(json.dumps(data))
    temp.replace(path)


def execute_saige_chain(chain, spec, plan, stop=None):
    exports = [j for j in chain if j['name'].endswith(('.export', '.index'))]
    associations = [j for j in chain if j['name'].endswith('.association')]
    if len(exports) != 2 or len(associations) != 1 or len(chain) != 3:
        raise ValueError('Unexpected SAIGE chromosome chain; regenerate gwas.sh plan')
    association = associations[0]
    temporary_outputs = {x for j in exports for x in j['outputs']}
    durable = [x for j in chain for x in j['inputs'] if x not in temporary_outputs]
    # Script/executable changes must invalidate the receipt as well.
    dependencies = [__file__]
    for j in chain:
        executable = shutil.which(j['cmd'][0])
        if executable: dependencies.append(executable)
        dependencies += [x for x in j['cmd'][1:] if x.endswith('.R') and Path(x).is_file()]
    signature = digest([chain, fingerprint(durable + dependencies), 'transient-dosage-v1'])
    receipt = Path(plan).parent/(association['name']+'.transient.done.json')
    outputs = association['outputs']
    if not spec.get('replace') and valid(receipt, signature, outputs):
        print('RESUME '+association['name']+' without dosage export', flush=True)
        return
    receipt.unlink(missing_ok=True)
    # Remove the legacy association receipt before attempting replacement.
    (Path(plan).parent/(association['name']+'.done.json')).unlink(missing_ok=True)
    shared = os.environ.get('GU_SAIGE_CACHE_DIR')
    if shared:
        run_with_cache(Path(shared), exports, association, plan, stop)
    else:
        # A standalone generated run.cmd owns and cleans its dosage cache.
        with tempfile.TemporaryDirectory(prefix='.saige-dosage.', dir=Path(plan).parent) as temporary:
            run_with_cache(Path(temporary), exports, association, plan, stop)
    if not all(Path(x).is_file() and Path(x).stat().st_size for x in outputs):
        raise RuntimeError('Missing SAIGE association result')
    record(receipt, signature, outputs)


def run_with_cache(root, exports, association, plan, stop=None):
    export, index = exports
    command = export['cmd']
    old_prefix = command[command.index('--out')+1]

    def rewrite(value, prefix):
        if value == old_prefix or value.startswith(old_prefix+'.'):
            return prefix+value[len(old_prefix):]
        if '=' in value:
            left, right = value.split('=', 1)
            if right == old_prefix or right.startswith(old_prefix+'.'):
                return left+'='+prefix+right[len(old_prefix):]
        return value

    normalized = [[rewrite(x, '@dosage') for x in j['cmd']] for j in exports]
    key = digest([normalized, fingerprint(export['inputs']),
                  fingerprint([shutil.which(j['cmd'][0]) for j in exports]), 'DS-force-v1'])
    root.mkdir(parents=True, exist_ok=True)
    cache = root/key
    with (root/(key+'.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        prefix = str(cache/'dosage')
        generated = [rewrite(x, prefix) for j in exports for x in j['outputs']]
        marker = cache/'complete.json'
        if not valid(marker, key, generated):
            shutil.rmtree(cache, ignore_errors=True)
            cache.mkdir()
            try:
                for job in exports:
                    run(job, [rewrite(x, prefix) for x in job['cmd']], plan, stop)
                if not all(Path(x).is_file() and Path(x).stat().st_size for x in generated):
                    raise RuntimeError('Incomplete dosage export/index')
                record(marker, key, generated)
            except BaseException:
                shutil.rmtree(cache, ignore_errors=True)
                raise
        else:
            print('REUSE shared SAIGE dosage '+key, flush=True)
        run(association, [rewrite(x, prefix) for x in association['cmd']], plan, stop)


def run(job, command, plan, stop=None):
    log = Path(plan).parent/(job['name']+'.run.log')
    with log.open('w') as handle:
        proc = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT)
        try:
            while proc.poll() is None:
                if stop is not None and stop.is_set(): raise RuntimeError('SAIGE run cancelled')
                try: proc.wait(timeout=0.2)
                except subprocess.TimeoutExpired: pass
            if proc.returncode: raise RuntimeError(f'{job["name"]} failed; see {log}')
        finally:
            if proc.poll() is None:
                proc.terminate()
                try: proc.wait(timeout=10)
                except subprocess.TimeoutExpired: proc.kill(); proc.wait()
