#!/usr/bin/env python3
"""Check the external annotator before starting expensive genome preparation.

Keep runtime changes outside f/: HERV receipts and active scan checkpoints
hash those scientific source files.
"""
from datetime import datetime
import importlib.util
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / 'f/1.herv.py'


def load_engine():
    sys.path.insert(0, str(ENGINE.parent))
    spec = importlib.util.spec_from_file_location('herv_preflight_engine', ENGINE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def short_error(log):
    lines = log.read_text(errors='replace').splitlines()
    errors = [line.strip() for line in lines if re.search(
        r'error|fatal|could not|cannot|not found|no such|permission denied', line, re.I)]
    return (errors[0] if errors else next((line.strip() for line in reversed(lines) if line.strip()), 'unknown failure'))[:300]


def preflight(args, engine):
    engine.require_library(args.library)
    # A verified completed run does not launch RepeatMasker again.
    if args.replace == 'FALSE' and (args.out / 'analysis.done.json').is_file():
        return
    log_dir = args.out / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f'repeatmasker_preflight_{datetime.now():%Y%m%d_%H%M%S_%f}.log'
    with tempfile.TemporaryDirectory(prefix='herv-preflight-') as directory:
        temp = Path(directory)
        # Fresh mode will archive the prior run.json before resolving the tool.
        selection = SimpleNamespace(repeatmasker=args.repeatmasker,
                                    out=temp if args.replace == 'TRUE' else args.out)
        executable = engine.repeatmasker_path(selection).absolute()
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError(f'RepeatMasker executable not found: {executable}')
        sequence = []
        with args.library.open() as source:
            for line in source:
                if line.startswith('>'):
                    if sequence:
                        break
                    continue
                sequence.extend(line.strip().upper())
                if len(sequence) >= 1000:
                    break
        sequence = ''.join(sequence[:1000])
        if not sequence or set(sequence) - set('ACGTN'):
            raise ValueError(f'Repeat library has no valid first FASTA sequence: {args.library}')
        fasta = temp / 'check.fa'
        fasta.write_text('>preflight\n' + sequence + '\n')
        command = [str(executable), '-e', 'rmblast', '-pa', '1',
                   '-lib', str(args.library.resolve()), '-dir', str(temp), str(fasta)]
        env = dict(os.environ, PATH=str(executable.parent) + os.pathsep + os.environ.get('PATH', ''))
        with log.open('w') as handle:
            result = subprocess.run(command, cwd=temp, env=env, stdin=subprocess.DEVNULL,
                                    stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode or not (temp / 'check.fa.out').is_file():
            raise ValueError(f'RepeatMasker preflight failed: {short_error(log)}; log={log}')
    print('[herv] PREFLIGHT RepeatMasker OK', flush=True)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    engine = load_engine()
    parser = engine.parser()
    parser.add_argument('--replace', choices=('TRUE', 'FALSE'), default='TRUE')
    parser.add_argument('--reuse-from', type=Path)
    args = parser.parse_args(argv)
    if not args.dry_run:
        if args.out is None:
            parser.error('The HERV launcher requires --out')
        args.out = args.out.expanduser().resolve()
        preflight(args, engine)
    os.execv(sys.executable, [sys.executable, str(ENGINE), *argv])


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print(f'[herv] ERROR {exc}', file=sys.stderr, flush=True)
        raise SystemExit(2)
