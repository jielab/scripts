"""Read the shared 0os.sh CPU budget; never duplicate machine defaults here."""
from functools import lru_cache
import os
from pathlib import Path
import subprocess


@lru_cache(maxsize=8)
def _default(mode):
    env = dict(os.environ)
    if mode: env.update(submit=mode, SUBMIT=mode)
    return int(subprocess.check_output(['bash', '-c',
        'source "$1" || exit; printf "%s" "$WORKFLOW_CORES"',
        'workflow-cpu', str(Path(__file__).with_name('0os.sh'))], env=env, text=True))


def core_budget(value=None):
    raw = value if value is not None else os.environ.get('WORKFLOW_CORES')
    if raw is None: raw = _default(os.environ.get('submit', os.environ.get('SUBMIT', 'auto')))
    try: n = int(raw)
    except (TypeError, ValueError): raise ValueError('--cores / WORKFLOW_CORES must be a positive integer') from None
    if str(n) != str(raw) or not 1 <= n <= 10000:
        raise ValueError('--cores / WORKFLOW_CORES must be an integer in 1..10000')
    return n


def worker_cores(value=None):
    """Never exceed the requested budget, scheduler allocation or CPU affinity."""
    limits = [core_budget(value)]
    if os.environ.get('LSB_DJOB_NUMPROC'): limits.append(int(os.environ['LSB_DJOB_NUMPROC']))
    if hasattr(os, 'sched_getaffinity'): limits.append(len(os.sched_getaffinity(0)))
    else: limits.append(os.cpu_count() or 1)
    return max(1, min(limits))
