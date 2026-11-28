#!/usr/bin/env python3
"""Expand ordinary IBDmix CLI cohorts into independent LSF jobs.

Only scheduling happens here; gu.sh and ygb_prepare.py retain sample selection
and all scientific processing. No genotype preparation runs on the login node.
"""
from __future__ import annotations

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path

import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

from processor import from_argv, validate_queue

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / '0f'))
import submit_jobs
from workflow_cpu import core_budget
from qiming_resources import cpu_queue, cpu_slots


def expand(value):
    """Numeric selectors preserve padding; arbitrary single postfixes are valid."""
    raw = value
    if value.startswith('[') and value.endswith(']'):
        value = value[1:-1]
    elif '[' in value or ']' in value:
        raise ValueError(f'Malformed selector: {raw}')
    result = []
    for token in value.split(','):
        match = re.fullmatch(r'(\d+)-(\d+)', token)
        if match:
            a, b = match.groups()
            if not 1 <= int(a) <= int(b) or int(b) - int(a) > 10000:
                raise ValueError(f'Invalid range: {token}')
            width = max(len(a), len(b)) if a.startswith('0') or b.startswith('0') else 0
            result.extend(str(n).zfill(width) for n in range(int(a), int(b) + 1))
        elif re.fullmatch(r'[A-Za-z0-9_-][A-Za-z0-9._-]*', token):
            result.append(token)
        else:
            raise ValueError(f'Invalid selector: {raw}')
    if len(set(result)) != len(result):
        raise ValueError(f'Duplicate selector: {raw}')
    return result


def get(args, flag, default=None):
    positions = [i for i, value in enumerate(args) if value == flag]
    if len(positions) > 1:
        raise ValueError(f'Duplicate option: {flag}')
    if not positions:
        return default
    i = positions[0]
    if i + 1 == len(args) or args[i + 1].startswith('--'):
        raise ValueError(f'{flag} requires a value')
    return args[i + 1]


def put(args, flag, value):
    args = list(args)
    if flag in args:
        args[args.index(flag) + 1] = value
    else:
        args += [flag, value]
    return args


def cohorts(args):
    postfix = get(args, '--postfix')
    posts = expand(postfix) if postfix is not None else ['']
    keep = get(args, '--keep-psam')
    if not keep:
        if len(posts) > 1:
            raise ValueError('Multiple postfixes require matching --keep-psam indices')
        return [(posts[0], put(args, '--postfix', posts[0]) if postfix is not None else args)]
    if '--keep' in args:
        raise ValueError('--keep and --keep-psam are mutually exclusive')
    try:
        path, size, selector = keep.split(',', 2)
        indices = expand(selector)
        if int(size) < 1 or any(not i.isdigit() or int(i) < 1 for i in indices):
            raise ValueError()
    except ValueError:
        raise ValueError('--keep-psam requires FILE,POSITIVE_SIZE,POSITIVE_INDEX_OR_RANGE') from None
    if len(posts) != len(indices):
        raise ValueError('--postfix and --keep-psam must expand to the same number of entries')
    if len(set(map(int, indices))) != len(indices):
        raise ValueError('Duplicate sample chunk indices')
    if len(posts) > 1 and any(not p.isdigit() or int(p) != int(i) for p, i in zip(posts, indices)):
        raise ValueError('Batch postfixes must match their corresponding sample chunk indices')
    # Bare PSAM filenames share the OS-derived target PFILE directory.
    if path and '/' not in path:
        target = get(args, '--target', '1kg')
        build = get(args, '--grch', '37')
        path = str(Path(os.environ['GU_REF_ROOT']) / target / build / ('imp' if target == 'ygb' else 'pfile') / path)
    # Count metadata only, never load PGENs or prepare samples on the login node.
    with open(path) as stream:
        header = next((line for line in stream if line.strip()), '')
        if 'IID' not in header.lstrip('#').split():
            raise ValueError(f'PSAM requires an IID header: {path}')
        count = sum(bool(line.strip()) for line in stream)
    if max(map(int, indices)) > (count + int(size) - 1) // int(size):
        raise ValueError(f'Chunk index exceeds PSAM sample count ({count})')
    result = []
    for post, index in zip(posts, indices):
        cmd = put(args, '--keep-psam', f'{path},{int(size)},{int(index)}')
        if postfix is not None:
            cmd = put(cmd, '--postfix', post)
        result.append((post, cmd))
    return result


def preflight_axt(args):
    """Cheap reference inventory before spending hours converting each cohort."""
    if get(args, '--grch', '37') != '37':
        return
    raw = get(args, '--chr', '1-22,X')
    chromosomes = []
    for token in re.split(r'[,\s]+', raw.strip()):
        token = re.sub('^chr', '', token, flags=re.I).upper()
        if re.fullmatch(r'\d+-\d+', token):
            lo, hi = map(int, token.split('-')); chromosomes.extend(map(str, range(lo, hi + 1)))
        else:
            chromosomes.append('X' if token == '23' else token)
    if any(c not in {str(i) for i in range(1,23)} | {'X'} for c in chromosomes):
        raise ValueError('Invalid --chr selection')
    root = Path(os.environ.get('IBDMIX_AXT_DIR', str(Path(os.environ['GU_ANNOT_ROOT']) / 'axt/37')))
    missing = []
    for chrom in chromosomes:
        # Existing scientific policy permits missing chrX AXT (reported in output).
        if chrom == 'X':
            continue
        for species in ('panTro2', 'ponAbe2', 'rheMac2'):
            name = f'chr{chrom}.hg19.{species}.synNet.axt.gz'
            nested = root / ('vs' + species[0].upper() + species[1:])
            folder = nested if (nested / name).exists() or (nested.is_dir() and not (root / name).exists()) else root
            checksums = folder / 'md5sum.txt'
            if not checksums.is_file():
                checksums = folder / 'checksums' / (species + '.md5sum.txt')
            for file in (folder / name, checksums):
                if not file.is_file() or file.stat().st_size == 0:
                    missing.append(str(file))
    if missing:
        raise ValueError('Missing AXT resources (no jobs submitted): ' + ', '.join(dict.fromkeys(missing)))


def result_base(args):
    explicit = get(args, '--dir-out')
    if explicit:
        return Path(explicit).expanduser().resolve()
    analysis = Path(os.environ.get('GU_ANALYSIS_BASE', str(Path(os.environ.get('ANALYSIS_ROOT', configured_path('ANALYSIS_ROOT'))) / 'gu')))
    return (analysis / 'ibdmix').resolve()


def preflight_inputs(args):
    import hpc
    target, build = get(args, '--target', '1kg'), get(args, '--grch', '37')
    prefix = get(args, '--dir-target', str(Path(os.environ['GU_REF_ROOT']) / target / build / ('imp' if target == 'ygb' else 'pfile') / 'chr'))
    hpc.discover_genotypes(prefix, hpc.chroms(get(args, '--chr', '1-22,X')))
    for flag in ('--keep', '--sample-panel', '--loci'):
        value = get(args, flag)
        if value and (not Path(value).is_file() or not Path(value).stat().st_size):
            raise ValueError(f'missing {Path(value).expanduser().absolute()} ({flag})')


def preflight_publication(mode):
    if mode == 'hpc':
        return  # R package checks run in the allocated worker after activation.
    interpreter = os.environ.get('WORKFLOW_RSCRIPT') or os.environ.get('GU_RESULTS_RSCRIPT', str(Path(os.environ.get('GU_ENV_PREFIX', configured_path('DIR_WORK', '.conda/envs/gu'))) / 'bin/Rscript'))
    if not Path(interpreter).is_file() or not os.access(interpreter, os.X_OK):
        raise ValueError(f'missing {interpreter}; run gu/install.sh before submitting')
    for name in ('ygb_prepare.py', 'ibdmix.py', 'ibdmix.sh', '0.common.py', '0.common.R', 'results.R', 'results.py'):
        path = ROOT / 'f' / name
        if not path.is_file(): raise ValueError(f'missing {path}')
    env = {k:v for k,v in os.environ.items() if k not in ('R_HOME', 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE')}
    subprocess.run([interpreter, '--vanilla', '-e', 'p<-c("data.table","jsonlite","digest","DBI","RSQLite","R.utils"); m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]; if(length(m)) stop("missing R packages: ",paste(m,collapse=", "))'], env=env, check=True)


def worker_environment():
    """Retain resolved project overrides when a saved .cmd is submitted later.

    A copied runtime.env contains default expressions, not the caller's resolved
    DIR_*/GU_* values. LSF inheritance alone does not preserve those values when
    --plan-only and the later submission happen in different shells.
    """
    exact = {'DIR_WORK', 'DIR_DATA', 'DIR_SCRATCH', 'DIR_SOFTWARE', 'CHECK_GRCH_SNP_LIST', 'WORKFLOW_CORES'}
    prefixes = ('DIR_', 'HPC_', 'RESULTS_', '_OS_DEFAULT_', 'GU_', 'IBDMIX_', 'TRACE_', 'PHYML_', 'ARG_', 'YGB_')
    excluded = {'GU_SUBMIT_PYTHON', 'GU_HPC_PYTHON',
                'GU_CMD_WORKER', 'GU_ENV_REEXEC'}
    return {key: value for key, value in os.environ.items()
            if key not in excluded and (key in exact or key.startswith(prefixes))}


def main(argv=None):
    submission_cwd = Path.cwd().resolve()
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--queue')
    parser.add_argument('--slots', type=int)
    parser.add_argument('--walltime', default=os.environ.get('HPC_WALLTIME', '120:00'))
    raw = sys.argv[1:] if argv is None else argv
    raw = [item for arg in raw for item in (arg.split('=', 1) if arg.startswith('--') and '=' in arg else [arg])]
    parser.add_argument('--plan-only', action='store_true')
    opts, args = parser.parse_known_args(raw)
    if not args or args[0] != 'ibdmix':
        raise ValueError('Only ordinary ibdmix requests are accepted')
    # Resolve caller-relative inputs before entering the job working directory.
    for i, flag in enumerate(args[:-1]):
        if flag.startswith('--dir-') or flag in ('--archaic-path', '--sample-panel', '--loci', '--keep', '--keep-psam'):
            value, sep, suffix = args[i + 1].partition(',')
            if value and not Path(value).is_absolute():
                args[i + 1] = str(submission_cwd / value) + sep + suffix
    mode = os.environ.get('submit', 'hpc')
    processors = from_argv('ibdmix', args)
    processor = processors['ibdmix']
    use_gpu = processor == 'gpu'
    if any(x in args for x in ('--run-cmd', '--replace-ibdmix', '--replace-trace', '--replace-phyml')):
        raise ValueError('Use --plan-only and --replace TRUE|FALSE')
    opts.queue = opts.queue or (os.environ.get('GU_GPU_QUEUE', os.environ.get('HPC_GPU_QUEUE', '4a100-40'))
                               if use_gpu else os.environ.get('GU_QUEUE', os.environ.get('HPC_CPU_QUEUE', '33')))
    opts.slots = opts.slots if opts.slots is not None else (16 if use_gpu else core_budget())
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', opts.queue) or opts.slots < 1:
        raise ValueError('Invalid queue or slots')
    if mode == 'hpc':
        validate_queue(opts.queue, use_gpu)
    if not re.fullmatch(r'\d+:[0-5]\d', opts.walltime):
        raise ValueError('--walltime must be HOURS:MINUTES')
    memory_default = (os.environ.get('GU_GPU_MEMORY_GIB', os.environ.get('HPC_GPU_MEMORY_GIB', '160'))
                      if use_gpu and mode == 'hpc' else os.environ.get('MEMORY_DEFAULT_GIB', '160'))
    memory = get(args, '--memory-cap', memory_default + 'G').upper()
    m = re.fullmatch(r'([1-9]\d*)G', memory)
    if not m:
        raise ValueError('HPC --memory-cap must be a positive integer in G, e.g. 32G')
    if mode == 'hpc' and not use_gpu:
        hours, minutes = map(int, opts.walltime.split(':'))
        opts.queue, opts.slots = cpu_slots(opts.queue, opts.slots, int(m[1]), hours + minutes / 60)
    jobs = int(get(args, '--jobs', str(min(opts.slots, core_budget(), max(1, int(m[1])//8)))))
    if jobs < 1 or jobs > opts.slots:
        raise ValueError('--jobs must be between 1 and --slots')
    target = get(args, '--target', '1kg')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', target) or target in ('.', '..'):
        raise ValueError('Invalid target')
    base = result_base(args)
    pairs = cohorts(args)
    preflight_inputs(args)
    preflight_publication(mode)
    preflight_axt(args)
    env_file = ROOT / 'f' / 'hpc.runtime.env'
    if not env_file.is_file():
        raise ValueError(f'Missing runtime: {env_file}')
    # Keep every cohort's submission artifacts under its result directory.
    scripts = []
    for post, cmd in pairs:
        parent = base / (target + post) / 'jobs'
        parent.mkdir(parents=True, exist_ok=True)
        directory = parent / 'ibdmix'
        directory.mkdir(exist_ok=True)
        environment = worker_environment()
        environment["WORKFLOW_CORES"] = str(min(core_budget(), opts.slots))
        environment['GU_CODE_ROOT'] = str(ROOT.parent)
        for flag, value in [('--dir-out', str(base)), ('--ibdmix-run', processor),
                            ('--foreground', 'TRUE'),
                            ('--resource-guard', 'scheduler' if mode == 'hpc' else get(args, '--resource-guard', 'auto')), ('--memory-cap', memory),
                            ('--jobs', str(jobs))]:
            cmd = put(cmd, flag, value)
        name = f"gu.ibdmix.{target}" + (f".{post}" if post else "")
        path = directory / (name + '.cmd')
        quote = shlex.quote
        script = '\n'.join([
            '#!/bin/bash -l', f'#BSUB -q {opts.queue}', f'#BSUB -n {opts.slots}',
            f'#BSUB -J {name}',
 f'#BSUB -R "span[hosts=1] rusage[mem={m[1]}GB/host]"',
            *(['#BSUB -gpu "num=1/host:mode=exclusive_process:j_exclusive=yes"'] if use_gpu and mode == 'hpc' else []),
            f'#BSUB -o "{directory}/{name}.%J.out"',
            f'#BSUB -e "{directory}/{name}.%J.err"',
            'set -euo pipefail',
            f'export submit={quote(mode)} SUBMIT={quote(mode)} WORKFLOW_WORKER=1',
            ': "${LSB_JOBID:?Submit this .cmd with bsub}"' if mode == 'hpc' else ':',
            f'export WORKFLOW_DEVICE={processor}',
            *(["export CUDA_VISIBLE_DEVICES=''"] if not use_gpu else []),
            *[f'export {key}={quote(value)}' for key, value in sorted(environment.items())],
            f'source {quote(str(env_file))}',
            *(['source ' + quote(str(ROOT.parent / '0f' / '0os.sh')), 'export GU_RESULTS_RSCRIPT="$WORKFLOW_RSCRIPT"'] if mode == 'hpc' else []),
            f'cd {quote(str(directory))}',
            'exec ' + shlex.join(['bash', str(ROOT / 'gu.sh')] + cmd), '',
        ])
        path.write_text(script)
        path.chmod(0o755)
        scripts.append(path)
    print(f'Generated {len(scripts)} independent .cmd files:', flush=True)
    for path in scripts:
        print(path, flush=True)
    if opts.plan_only:
        print('Plan only: command files written; no scientific work executed.')
        return 0
    if mode == 'local':
        submit_jobs.local_one = submit_jobs.local_foreground_one
    submit_jobs.run(scripts, mode=mode, wait=os.environ.get('SUBMIT_WAIT') == 'TRUE')

    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(2)
