"""Apply a separate cgroup memory ceiling to reference download/masking only."""
import argparse
import os
from pathlib import Path
import subprocess
import time


def scope_command(memory_gb, command, parent_unit=None):
    if memory_gb < 1 or not command:
        raise ValueError('A positive preparation memory limit and a command are required')
    limit = memory_gb * 2**30
    options = ['systemd-run', '--user', '--scope', '--quiet', '--collect',
            f'--unit=metag-reference-{os.getpid()}-{time.time_ns()}',
            f'--property=MemoryHigh={limit * 9 // 10}', f'--property=MemoryMax={limit}',
            '--property=MemorySwapMax=1G', '--property=OOMPolicy=stop',
            '--property=KillMode=control-group']
    if parent_unit:
        options += [f'--property=BindsTo={parent_unit}', f'--property=After={parent_unit}']
    return options + ['--', *command]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--memory-gb', type=int, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    print(f'[metag] Reference preparation RAM cap: {args.memory_gb} GiB; swap cap: 1 GiB', flush=True)
    # Fail visibly if the cap cannot be enforced; never silently run uncapped.
    parent_unit = None
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        unit = line.rsplit('/', 1)[-1]
        if unit.startswith('metag-workflow-') and unit.endswith('.scope'):
            parent_unit = unit
    # A preparation scope is a sibling, so bind its lifetime to the workflow.
    # Cancelling the outer workflow must not leave its masker processes behind.
    result = subprocess.run(scope_command(args.memory_gb, command, parent_unit))
    if result.returncode:
        print('[metag] Reference preparation failed. Completed genome receipts are retained; '
              'inspect the log before retrying. For memory failures, reduce --mask-jobs '
              'or explicitly increase --prepare-memory-gb.', flush=True)
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
