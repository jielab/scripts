#!/usr/bin/env python3
"""Expand ordinary IBDmix CLI cohorts into independent LSF jobs.

Only scheduling happens here; gu.sh and ukb_prepare.py retain sample selection
and all scientific processing. No genotype preparation runs on the login node.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import hashlib
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / '0f'))
import submit_jobs


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
        return [(posts[0], args)]
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
    target, build = get(args, '--target', '1kg'), get(args, '--grch', '37')
    prefix = get(args, '--dir-target')
    base = Path(os.environ.get('GU_REF_ROOT', '/work/sph-huangj/data/gen')) / target / build
    if prefix:
        parent = Path(prefix).expanduser().resolve().parent
        if parent.name in ('pfile', 'imp', 'hap', 'vcf', 'bfile', 'bed'):
            base = parent.parent
    # Keep every cohort inside the input's build directory (never 37013/).
    if get(args, '--postfix'):
        base = base / 'cohort'
    return base.resolve()


def preflight_inputs(args):
    import hpc
    target, build = get(args, '--target', '1kg'), get(args, '--grch', '37')
    prefix = get(args, '--dir-target', str(Path(os.environ['GU_REF_ROOT']) / target / build / 'pfile/chr'))
    hpc.discover_genotypes(prefix, hpc.chroms(get(args, '--chr', '1-22,X')))
    for flag in ('--keep', '--sample-panel', '--loci'):
        value = get(args, flag)
        if value and (not Path(value).is_file() or not Path(value).stat().st_size):
            raise ValueError(f'missing {Path(value).expanduser().absolute()} ({flag})')


def preflight_publication():
    interpreter = os.environ.get('GU_RESULTS_RSCRIPT', str(Path(os.environ.get('GU_ENV_PREFIX', '/work/sph-huangj/.conda/envs/gu')) / 'bin/Rscript'))
    if not Path(interpreter).is_file() or not os.access(interpreter, os.X_OK):
        raise ValueError(f'missing {interpreter}; run gu/install.sh before submitting')
    for name in ('ukb_prepare.py', 'ibdmix.py', 'ibdmix.sh', '0.common.py', '0.common.R', 'results.R', 'results.py'):
        path = ROOT / 'f' / name
        if not path.is_file(): raise ValueError(f'missing {path}')
    env = {k:v for k,v in os.environ.items() if k not in ('R_HOME', 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE')}
    subprocess.run([interpreter, '--vanilla', '-e', 'p<-c("data.table","jsonlite","digest","DBI","RSQLite","R.utils"); m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]; if(length(m)) stop("missing R packages: ",paste(m,collapse=", "))'], env=env, check=True)


def snapshot_code(directory):
    """Running cohorts must not lose helpers when the live scripts are updated."""
    destination = directory / 'code'
    destination.mkdir()
    ignore = shutil.ignore_patterns('__pycache__', '*.pyc', '.migration-backup', '.git', '*.pre-*')
    for name in ('0f', 'gu', 'gwas', 'grid', 'le8', 'allele3', 'ld', 'metag'):
        source = ROOT.parent / name
        if source.is_dir():
            shutil.copytree(source, destination / name, ignore=ignore)
    manifest = {str(p.relative_to(destination)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in destination.rglob('*') if p.is_file()}
    (directory / 'code.sha256.json').write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
    return destination / 'gu'


def worker_environment():
    """Retain resolved project overrides when a saved .cmd is submitted later.

    A copied runtime.env contains default expressions, not the caller's resolved
    DIR_*/GU_* values. LSF inheritance alone does not preserve those values when
    --run-cmd FALSE and the later submission happen in different shells.
    """
    exact = {'DIR_WORK', 'DIR_DATA', 'DIR_SCRATCH', 'DIR_SOFTWARE', 'CHECK_GRCH_SNP_LIST'}
    prefixes = ('DIR_', 'HPC_', 'RESULTS_', '_OS_DEFAULT_', 'GU_', 'IBDMIX_', 'TRACE_', 'PHYML_', 'ARG_', 'UKB_')
    excluded = {'GU_WORKER_CODE_ROOT', 'GU_SUBMIT_PYTHON', 'GU_HPC_PYTHON',
                'GU_CMD_WORKER', 'GU_ENV_REEXEC'}
    return {key: value for key, value in os.environ.items()
            if key not in excluded and (key in exact or key.startswith(prefixes))}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--queue', default=os.environ.get('GU_QUEUE', os.environ.get('HPC_CPU_QUEUE', '33')))
    parser.add_argument('--slots', type=int, default=40)
    parser.add_argument('--walltime', default='144:00')
    opts, args = parser.parse_known_args(sys.argv[1:] if argv is None else argv)
    if not args or args[0] != 'ibdmix':
        raise ValueError('Only ordinary ibdmix requests are accepted')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', opts.queue) or opts.slots < 1:
        raise ValueError('Invalid queue or slots')
    if not re.fullmatch(r'\d+:\d{2}', opts.walltime):
        raise ValueError('--walltime must be HOURS:MINUTES')
    mode = os.environ.get('submit', 'hpc')
    if mode == 'hpc' and get(args, '--GPU', 'FALSE').upper() != 'FALSE':
        raise ValueError('HPC batch submission is CPU-only; use --GPU FALSE')
    run = get(args, '--run-cmd', 'TRUE').upper()
    if run not in {'TRUE', 'FALSE'}:
        raise ValueError('--run-cmd must be TRUE or FALSE')
    memory = get(args, '--memory-cap', os.environ.get('MEMORY_DEFAULT_GIB', '160') + 'G').upper()
    m = re.fullmatch(r'([1-9]\d*)G', memory)
    if not m:
        raise ValueError('HPC --memory-cap must be a positive integer in G, e.g. 32G')
    jobs = int(get(args, '--jobs', '1'))
    if jobs < 1 or jobs > opts.slots:
        raise ValueError('--jobs must be between 1 and --slots')
    target = get(args, '--target', '1kg')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', target):
        raise ValueError('Invalid target')
    base = result_base(args)
    pairs = cohorts(args)
    preflight_inputs(args)
    preflight_publication()
    preflight_axt(args)
    env_file = ROOT / 'f' / 'hpc.runtime.env'
    if not env_file.is_file():
        raise ValueError(f'Missing runtime: {env_file}')
    # An invocation has its own scripts/logs, so another submission cannot
    # overwrite a queued job's command. Submission receipts survive failures.
    parent = base.parent / (base.name + '.jobs')
    parent.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='ibdmix-', dir=parent))
    frozen = snapshot_code(directory)
    (directory / 'runtime.env').write_bytes(env_file.read_bytes())
    scripts = []
    environment = worker_environment()
    code_prefixes = [str(ROOT.parent / name) for name in ('0f','gu','gwas','grid','le8','ld','metag','allele3')]
    for key, value in environment.items():
        if any(value == prefix or value.startswith(prefix + '/') for prefix in code_prefixes):
            environment[key] = str(frozen.parent) + value[len(str(ROOT.parent)):]
    environment['GU_CODE_ROOT'] = str(frozen.parent)
    for post, cmd in pairs:
        for flag, value in [('--dir-out', str(base)), ('--GPU', get(args, '--GPU', 'FALSE').upper()),
                            ('--run-cmd', 'TRUE'), ('--foreground', 'TRUE'),
                            ('--resource-guard', 'scheduler' if mode == 'hpc' else get(args, '--resource-guard', 'auto')), ('--memory-cap', memory),
                            ('--jobs', str(jobs))]:
            cmd = put(cmd, flag, value)
        name = f'ibdmix-{target}-{post or "all"}'
        path = directory / (name + '.cmd')
        quote = shlex.quote
        script = '\n'.join([
            '#!/bin/bash -l', f'#BSUB -q {opts.queue}', f'#BSUB -n {opts.slots}',
            f'#BSUB -W {opts.walltime}', f'#BSUB -J {name}',
            f'#BSUB -M {m[1]}GB',
            f'#BSUB -R "span[hosts=1] rusage[mem={m[1]}GB/host]"',
            f'#BSUB -o "{directory}/{name}.%J.out"',
            f'#BSUB -e "{directory}/{name}.%J.err"',
            'set -euo pipefail',
            f'export submit={quote(mode)} SUBMIT={quote(mode)} WORKFLOW_WORKER=1',
            f'export GU_WORKER_CODE_ROOT={quote(str(frozen))}',
            ': "${LSB_JOBID:?Submit this .cmd with bsub}"' if mode == 'hpc' else ':',
            f'export DATE_SCRATCH={quote(os.environ.get("DATE_SCRATCH", "2026-11-03"))}',
            *[f'export {key}={quote(value)}' for key, value in sorted(environment.items())],
            f'source {quote(str(directory / "runtime.env"))}',
            f'cd {quote(str(frozen))}',
            'exec ' + shlex.join(['bash', str(frozen / 'gu.sh')] + cmd), '',
        ])
        path.write_text(script)
        path.chmod(0o755)
        scripts.append(path)
    print(f'Generated {len(scripts)} independent .cmd files: {directory}', flush=True)
    if run == 'FALSE':
        print('Dry run: command files only. Use 0f/submit_jobs.py with the selected submit mode to run them.')
        return 0
    submit_jobs.run(scripts, mode=mode, wait=os.environ.get('SUBMIT_WAIT') == 'TRUE')

    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(2)
