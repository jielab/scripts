"""Write a small-library build command without downloading or executing it."""
import argparse
import os
from pathlib import Path
import shlex
from common import load_config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module', choices=('human','sylph'), required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--root')
    parser.add_argument('--submit', choices=('hpc','local'), required=True)
    args = parser.parse_args()
    cfg = load_config(args.config, args.root or None)
    out = Path(cfg['work_dir']) / 'jobs/small-db' / (args.module + '.plan.cmd')
    out.parent.mkdir(parents=True, exist_ok=True)
    names = (('HUMAN_DIR','BWA_DIR','GRCH38_URL','BWA_URL','HUMAN_QUEUE','HUMAN_CORES','HUMAN_MEM_GIB','HUMAN_HOURS')
             if args.module=='human' else ('SYLPH_DIR','SYLPH_TOOL_DIR','SYLPH_BUCKETS','SYLPH_QUEUE','SYLPH_THREADS','SYLPH_MEM_GIB','SYLPH_HOURS','AUDIT_DIR'))
    command = ['bash', str(Path(__file__).resolve().parents[1] / '1.build-lib.sh'), args.module,
               '--config', str(Path(args.config.split('#',1)[0]).expanduser().resolve()) + ('#'+args.config.split('#',1)[1] if '#' in args.config else ''),
               '--submit', args.submit]
    if args.root:
        command += ['--root', str(Path(args.root).expanduser().absolute())]
    lines = ['#!/usr/bin/env bash', 'set -euo pipefail']
    lines += ['export '+key+'='+shlex.quote(os.environ[key]) for key in names if key in os.environ]
    lines += ['exec '+shlex.join(command)]
    out.write_text('\n'.join(lines)+'\n');out.chmod(0o755)
    print('Plan only; no downloads or builds performed. Run: bash '+shlex.quote(str(out)))


if __name__ == '__main__':
    main()
