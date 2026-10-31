#!/usr/bin/env python3
"""One-command independent rebuild of the public McManus et al. LD objective.

The numerical objective is exact for the supplied phased r2 values and tested D.
The finite-D stability rule and degenerate-coordinate conventions are explicitly
implementation choices, not a claim of identity to restricted author outputs.
"""
from __future__ import annotations
import argparse
import contextlib
import csv
import fcntl
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import stat
import time
import tempfile
from universal import compute_consensus, read_bed

POPS = ('AFR', 'EAS', 'EUR', 'SAS')
EXPECTED = dict(AFR=661, EAS=504, EUR=503, SAS=489)
HERE = Path(__file__).resolve().parent
VERSION = 'ld-rebuild-v2-20261020'
LEGACY_PGEN_IMPLEMENTATIONS = {'4f94282bfd912895c454dbb98b85b8340f9a39c9e01ec8d72c41e17ab10f4889'}
# This exact reader-only revision preserves the original arithmetic and input
# validation. Keep the producer hash in existing receipts; accept it only for
# the explicitly audited successor, never for arbitrary future source edits.
PARTITION_EQUIVALENT_SOURCES = {
    # 2026-10-24: buffer views, memrchr and a reused Zstandard decoder only.
    # Differential BED/QC checks and exhaustive small partitions preserve the
    # parser's rejection rules, floating conversion and every DP operation.
    '439e395d3758abc4c8bae7f70b1f53ae849f348f4565f904fa0d0efb48bed865': {
        '5c6254ee9d4629e46f2491c6cc8b7aa5ee36fe44f6b5321c8385e0d4ee9621fd',
        '4383c6c2ebd489224904b5c61b156ee458232cd397ac0f531818c3db57e0cd88',
        '81ab7a235173b4dceae91d3612cc390e8a3e0ea9ae10c45fa0327177eeb37929',
    },
    '5c6254ee9d4629e46f2491c6cc8b7aa5ee36fe44f6b5321c8385e0d4ee9621fd': {
        '4383c6c2ebd489224904b5c61b156ee458232cd397ac0f531818c3db57e0cd88',
        '81ab7a235173b4dceae91d3612cc390e8a3e0ea9ae10c45fa0327177eeb37929',
    },
    '4383c6c2ebd489224904b5c61b156ee458232cd397ac0f531818c3db57e0cd88': {
        '81ab7a235173b4dceae91d3612cc390e8a3e0ea9ae10c45fa0327177eeb37929',
    },
}

class Stop(RuntimeError):
    pass

class Tee:
    """Mirror Python status/errors to disk; PLINK retains its own native logs."""
    def __init__(self, terminal, logfile):
        self.terminal = terminal
        self.logfile = logfile

    def write(self, text):
        self.terminal.write(text)
        if self.logfile is not None:
            try:
                self.logfile.write(text)
            except OSError as exc:
                self.logfile = None
                self.terminal.write(f'\n[LD] Could not continue writing run log: {exc}\n')
        return len(text)

    def flush(self):
        self.terminal.flush()
        if self.logfile is not None:
            try:
                self.logfile.flush()
            except OSError as exc:
                self.logfile = None
                self.terminal.write(f'\n[LD] Could not flush run log: {exc}\n')

def note(message):
    print(f'[LD] {message}', flush=True)

def atomic(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    try:
        with temporary.open('w', encoding='utf-8', newline='') as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)

def uniquify_pvar_ids(path):
    """Disambiguate repeated IDs without changing variant order or genotypes."""
    path = Path(path)
    seen = set()
    counts = {}
    changed = 0
    temporary = path.with_name(path.name + '.unique.tmp')
    try:
        with path.open() as source, temporary.open('w') as target:
            for line in source:
                if line.startswith('#') or not line.strip():
                    target.write(line)
                    continue
                fields = line.rstrip('\n').split('\t')
                original = fields[2]
                candidate = original
                count = counts.get(original, 1)
                while candidate in seen:
                    count += 1
                    candidate = f'{original}:record{count}'
                counts[original] = count
                seen.add(candidate)
                if candidate != original:
                    changed += 1
                    fields[2] = candidate
                    line = '\t'.join(fields) + '\n'
                target.write(line)
        if changed:
            temporary.replace(path)
        return changed
    finally:
        temporary.unlink(missing_ok=True)


def write_json(path, value):
    atomic(path, json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + '\n')

def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(2**20), b''):
            digest.update(block)
    return digest.hexdigest()

def stamp(path):
    path = Path(path).resolve()
    st = path.stat()
    if not st.st_size:
        raise Stop(f'Empty input: {path}')
    return {'path': str(path), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}

def parse_chroms(value):
    result = []
    for token in value.split(','):
        token = token.strip().lower().removeprefix('chr')
        if '-' in token:
            a, b = map(int, token.split('-', 1))
            if a > b:
                raise Stop('CHRS range must increase')
            result.extend(range(a, b + 1))
        else:
            result.append(int(token))
    if not result or len(result) != len(set(result)) or any(c < 1 or c > 22 for c in result):
        raise Stop('CHRS must contain unique autosomes 1-22')
    return sorted(result)

def positive_int(name, default):
    value = int(os.environ.get(name, default))
    if value < 1:
        raise Stop(f'{name} must be positive')
    return value

def load_targets(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    result = {(r['population'], int(r['chrom'])): int(r['blocks']) for r in rows}
    if len(result) != len(rows) or any(n < 1 for n in result.values()):
        raise Stop('Invalid or duplicate published-count rows')
    return result

def table(path):
    with Path(path).open() as handle:
        lines = (line for line in handle if line.strip())
        header = next(lines).lstrip('#').split()
        for line in lines:
            fields = line.split()
            if len(fields) != len(header):
                raise Stop(f'Malformed table: {path}')
            yield dict(zip(header, fields))

def changed_boundaries(old, new):
    # Count differing block intervals, not merely differences in block counts.
    return len(set(old).symmetric_difference(new))

def ld_table_requirements(positions, distance_bp):
    """Count pairs and a strict byte lower bound for plain cols=id VCOR.

    IDs are v{zero_based_pvar_row}. Each pair needs both IDs, two tabs,
    at least one r2 character and a newline. Real r2 values usually need
    more characters, so passing this check does not guarantee enough space.
    """
    left = pairs = window_id_bytes = 0
    minimum_bytes = len('#ID_A\tID_B\tPHASED_R2\n')
    for right, pos in enumerate(positions):
        while pos - positions[left] > distance_bp:
            window_id_bytes -= len(str(left)) + 1
            left += 1
        id_bytes = len(str(right)) + 1
        count = right - left
        pairs += count
        minimum_bytes += window_id_bytes + count * (id_bytes + 4)
        window_id_bytes += id_bytes
    return pairs, minimum_bytes

def choose_distance(snapshots, grid, confirmations=2):
    """Smallest D whose complete suffix is identical, with enough larger tests.

A full-span test removes D as a restriction mathematically; otherwise this is
explicitly an empirical, finite-grid stability rule. Counts never rank D.
"""
    for i, distance in enumerate(grid):
        tail = grid[i:]
        same = all(snapshots[d]['signature'] == snapshots[distance]['signature'] for d in tail)
        full_span = any(snapshots[d]['full_span_all_populations'] for d in tail)
        if same and (len(tail) >= confirmations + 1 or full_span):
            return distance
    return None

class Pipeline:
    def __init__(self):
        self.genome_build = os.environ.get('GENOME_BUILD', '38')
        if self.genome_build not in ('37', '38'):
            raise Stop('GENOME_BUILD must be 37 or 38')
        self.input_format = os.environ.get('INPUT_FORMAT', 'pgen')
        if self.input_format not in ('vcf', 'pgen'):
            raise Stop('INPUT_FORMAT must be vcf or pgen')
        self.pfile_dir = Path(os.environ.get('PFILE_DIR', f'/data/sph-huangj/gen/1kg/{self.genome_build}/pfile')).expanduser().resolve()
        self.out = Path(os.environ.get('OUT', f'/data/sph-huangj/analysis/ld/{self.genome_build}')).expanduser().resolve()
        self.work = self.out / 'rebuild_v2'
        self.tmp_parent = Path(os.environ.get('LD_TMP_DIR', str(self.out / 'tmp'))).expanduser().resolve()
        self.resume_dir = (Path(os.environ['LD_RESUME_DIR']).expanduser().resolve()
                           if os.environ.get('LD_RESUME_DIR') else None)
        self.vcf_dir = Path(os.environ.get('VCF_DIR', f'/data/sph-huangj/gen/1kg/{self.genome_build}/vcf')).expanduser().resolve()
        self.panel = Path(os.environ.get('PANEL', f'/data/sph-huangj/gen/1kg/{self.genome_build}/samples.txt')).expanduser().resolve()
        self.chroms = parse_chroms(os.environ.get('CHRS', '1-22'))
        self.threads = positive_int('THREADS', 8)
        self.parallel_chroms = positive_int('PARALLEL_CHROMS', 1)
        self.memory_mb = positive_int('MEMORY_MB', 4096)
        self.confirmations = positive_int('STABLE_CONFIRMATIONS', 2)
        # At least two independent larger-D comparisons by default.
        raw_grid = os.environ.get('D_GRID_KB', '250,500,1000,2000,5000,10000')
        self.grid = [int(x) for x in raw_grid.split(',')]
        if not self.grid or self.grid != sorted(set(self.grid)) or min(self.grid) < 1:
            raise Stop('D_GRID_KB must be a strictly increasing list of positive integers')
        cap = positive_int('MAX_D_KB', os.environ.get('PILOT_MAX_D_KB', '40000'))
        if self.grid[-1] > cap:
            raise Stop('D_GRID_KB exceeds MAX_D_KB')
        if cap > (2**63 - 1) // 1000:
            raise Stop('MAX_D_KB exceeds the supported integer range')
        while self.grid[-1] < cap:
            self.grid.append(min(2 * self.grid[-1], cap))
        if len(self.grid) > 64:
            raise Stop('D_GRID_KB has too many states (maximum 64)')
        self.min_free = float(os.environ.get('MIN_FREE_GB', '5')) * 10**9
        self.reference_tolerance = float(os.environ.get('REFERENCE_MAX_REL_DEVIATION', '0.10'))
        self.reference_policy = os.environ.get('REFERENCE_COUNT_POLICY', 'off')
        if self.reference_policy not in ('warn', 'error', 'off'):
            raise Stop('REFERENCE_COUNT_POLICY must be warn, error or off')
        if self.genome_build == '38' and self.reference_policy != 'off':
            raise Stop('GRCh38 requires REFERENCE_COUNT_POLICY=off: GRCh37 cohort counts do not apply')
        if not math.isfinite(self.min_free) or self.min_free < 0:
            raise Stop('MIN_FREE_GB must be finite and nonnegative')
        if not math.isfinite(self.reference_tolerance) or self.reference_tolerance < 0:
            raise Stop('REFERENCE_MAX_REL_DEVIATION must be finite and nonnegative')
        self.coincident_policy = os.environ.get('COINCIDENT_POLICY', 'group')
        if self.coincident_policy not in ('group', 'error'):
            raise Stop('COINCIDENT_POLICY must be group or error')
        if os.environ.get('LD_MODE', 'strict') != 'strict':
            raise Stop('The automatic pipeline requires strict LD; sparse cannot pass production gates')
        if os.environ.get('KEEP_VCOR', '0') != '0':
            raise Stop('This automatic pipeline always deletes temporary LD tables; remove KEEP_VCOR=1')
        if os.environ.get('D_KB'):
            raise Stop('D_KB is selected automatically; remove manual D_KB and run ./ld.sh')
        target_name = os.environ.get('LD_TARGETS', '')
        if self.reference_policy != 'off' and not target_name:
            raise Stop('Set LD_TARGETS to an external reference-count TSV when enabling REFERENCE_COUNT_POLICY')
        target_file = Path(target_name) if target_name else None
        self.targets = load_targets(target_file) if self.reference_policy != 'off' else {}
        for c in (self.chroms if self.targets else []):
            for pop in (*POPS, 'UNIVERSAL'):
                if (pop, c) not in self.targets:
                    raise Stop(f'Missing reference count for {pop} chr{c}: {target_file}')
        plink = os.environ.get('PLINK2', 'plink2')
        self.plink = shutil.which(plink)
        self.cxx = shutil.which(os.environ.get('CXX', 'g++'))
        self.tool = None
        self.bin = self.out / 'work' / 'bin' / 'ld_partition'
        self.source_hash = sha(HERE / 'ld_partition.cpp')
        self.ld_storage = os.environ.get('LD_STORAGE', 'blocked-zstd')
        if self.ld_storage not in ('plain', 'blocked-zstd'):
            raise Stop('LD_STORAGE must be plain or blocked-zstd')
        self.pgen_cache = {}
        self.child = None
        self.samples = {}
        self.sample_counts = {}

    def input_pfiles(self, chrom):
        prefix = self.pfile_dir / f'chr{chrom}'
        pvar = Path(str(prefix) + '.pvar.zst')
        if not pvar.is_file():
            pvar = Path(str(prefix) + '.pvar')
        files = {'pgen': Path(str(prefix) + '.pgen'), 'pvar': pvar,
                 'psam': Path(str(prefix) + '.psam')}
        pgi = Path(str(prefix) + '.pgen.pgi')
        if pgi.is_file():
            files['pgi'] = pgi
        for path in files.values():
            if not path.is_file() or not path.stat().st_size:
                raise Stop(f'Missing or empty PGEN input: {path}')
        return files

    @contextlib.contextmanager
    def run_logging(self):
        logdir = self.work / 'logs'
        logdir.mkdir(parents=True, exist_ok=True)
        path = logdir / f'run.{time.strftime("%Y%m%d_%H%M%S")}.{os.getpid()}.{time.time_ns()}.log'
        with path.open('x', buffering=1) as handle:
            with contextlib.redirect_stdout(Tee(sys.stdout, handle)), contextlib.redirect_stderr(Tee(sys.stderr, handle)):
                note(f'Run log: {path}')
                note(f'Started {time.strftime("%Y-%m-%d %H:%M:%S %z")}; chromosomes={self.chroms}')
                started = time.monotonic()
                try:
                    yield
                except BaseException as exc:
                    # main() prints the final error after this context has closed.
                    with contextlib.suppress(OSError):
                        handle.write(f'[LD] STOP: {type(exc).__name__}: {exc}\n')
                    raise
                finally:
                    note(f'Run elapsed: {time.monotonic() - started:.1f} s; log: {path}')

    def show_error_log(self, path):
        if path is None:
            return
        with contextlib.suppress(OSError):
            with Path(path).open('rb') as handle:
                handle.seek(max(0, handle.seek(0, 2) - 16384))
                tail = handle.read().decode('utf-8', errors='replace').strip()
            if tail:
                note(f'{Path(path).name} (last 16 KiB):\n{tail}')

    def check_free_space(self, disk_path):
        free = shutil.disk_usage(disk_path).free
        if free < self.min_free:
            detail = ''
            with contextlib.suppress(OSError):
                size = sum(p.stat().st_size for p in Path(disk_path).glob('*.vcor') if not p.is_symlink())
                detail = f'; temporary .vcor size={size / 1e9:.3f} GB'
            raise Stop(f'Free disk below MIN_FREE_GB={self.min_free / 1e9:g}: '
                       f'{disk_path} had {free / 1e9:.3f} GB free before cleanup{detail}. '
                       'Deleting temporary files releases space, so free space after stopping can be much higher. '
                       'Set LD_TMP_DIR to a filesystem with more free space and rerun.')
        return free

    def run_command(self, command, disk_path=None, *, stderr_path=None):
        note(shlex.join(map(str, command)))
        # Keep the terminal session; a separate process group lets the parent
        # stop and reap the command and its descendants together on interruption.
        started = time.monotonic()
        error_handle = Path(stderr_path).open('wb') if stderr_path is not None else None
        try:
            self.child = subprocess.Popen(list(map(str, command)), process_group=0, stderr=error_handle)
            while True:
                try:
                    code = self.child.wait(timeout=3)
                    break
                except subprocess.TimeoutExpired:
                    if disk_path:
                        self.check_free_space(disk_path)
            if code:
                raise Stop(f'Command failed with exit {code}: {command[0]}')
        finally:
            if self.child is not None and self.child.poll() is None:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(self.child.pid, signal.SIGTERM)
                try:
                    self.child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(self.child.pid, signal.SIGKILL)
                    self.child.wait()
            self.child = None
            if error_handle is not None:
                error_handle.close()
                self.show_error_log(stderr_path)
            note(f'Command elapsed: {time.monotonic() - started:.1f} s')

    @contextlib.contextmanager
    def temporary(self, label):
        self.tmp_parent.mkdir(parents=True, exist_ok=True)
        self.check_free_space(self.tmp_parent)
        folder = Path(tempfile.mkdtemp(prefix='ld-auto-v2-', dir=self.tmp_parent))
        try:
            write_json(folder / 'OWNER.json', {'program': VERSION, 'output_root': str(self.out), 'label': label})
            yield folder
        finally:
            try:
                logdir = self.work / 'logs'
                logdir.mkdir(parents=True, exist_ok=True)
                for log in folder.glob('*.log'):
                    shutil.copy2(log, logdir / f'{label}.{log.name}')
            except OSError as exc:
                note(f'Could not preserve a temporary log: {exc}')
            finally:
                # Cleanup must run even when a full disk prevents log copying.
                # This directory contains only files created by this process.
                shutil.rmtree(folder)
                note(f'Temporary files removed: {label}')

    @contextlib.contextmanager
    def output_locks(self, *, wait=False):
        self.out.mkdir(parents=True, exist_ok=True)
        cache = self.out / 'work'
        cache.mkdir(parents=True, exist_ok=True)
        # Always acquire the legacy lock first, even if it does not yet exist.
        # Keep both lock files in place so all versions lock the same inodes.
        paths = (self.out / '.ld.lock', cache / '.ld.lock')
        with contextlib.ExitStack() as stack:
            locked_inodes = set()
            for path in paths:
                handle = stack.enter_context(path.open('a'))
                stat = os.fstat(handle.fileno())
                identity = (stat.st_dev, stat.st_ino)
                if identity in locked_inodes:
                    continue
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | (0 if wait else fcntl.LOCK_NB))
                except BlockingIOError as exc:
                    raise Stop(f'Another ld.sh process is using {self.out}; lock is busy: {path}') from exc
                locked_inodes.add(identity)
            yield

    def clean_abandoned(self):
        # Only after all output locks are held; never sweep unrelated files.
        self.tmp_parent.mkdir(parents=True, exist_ok=True)
        for folder in self.tmp_parent.glob('ld-auto-v2-*'):
            if folder.is_symlink() or not folder.is_dir():
                continue
            owner = read_json(folder / 'OWNER.json')
            if owner and owner.get('program') == VERSION and owner.get('output_root') == str(self.out):
                note(f'Removing an interrupted private temporary directory: {folder.name}')
                shutil.rmtree(folder)

    def prepare(self):
        if not self.plink or not self.cxx:
            raise Stop('Need PLINK2 and a C++17 compiler (g++); set PLINK2/CXX to their executable paths')
        self.work.mkdir(parents=True, exist_ok=True)
        self.tool = {'path': str(Path(self.plink).resolve()), 'sha256': sha(self.plink),
                     'version': subprocess.check_output([self.plink, '--version'], text=True).strip()}
        note(self.tool['version'])
        if not self.panel.is_file():
            raise Stop(f'Missing sample panel: {self.panel}')
        self.samples = {p: [] for p in POPS}
        seen = set()
        with self.panel.open() as handle:
            header = handle.readline().split()
            si, pi = header.index('sample'), header.index('super_pop')
            for line in handle:
                if not line.strip():
                    continue
                row = line.split()
                sample, pop = row[si], row[pi]
                if sample in seen:
                    raise Stop(f'Duplicate panel sample: {sample}')
                seen.add(sample)
                if pop in self.samples:
                    self.samples[pop].append(sample)
        for pop in POPS:
            if self.reference_policy != 'off' and len(self.samples[pop]) != EXPECTED[pop]:
                raise Stop(f'{pop}: expected {EXPECTED[pop]} people, found {len(self.samples[pop])}')
            if len(self.samples[pop]) < 2:
                raise Stop(f'{pop}: need at least two samples')
            self.sample_counts[pop] = len(self.samples[pop])
            self.samples[pop].sort()
            atomic(self.out / 'meta' / f'{pop}.keep', '#FID\tIID\n' + ''.join(f'{s}\t{s}\n' for s in self.samples[pop]))
        for chrom in self.chroms:
            if self.input_format == 'pgen':
                files = self.input_pfiles(chrom)
                rows = list(table(files['psam']))
                ids = [r['IID'] for r in rows]
                if len(ids) != len(set(ids)):
                    raise Stop(f'chr{chrom}: duplicate IID in source PSAM')
                required = {s for group in self.samples.values() for s in group}
                if not required.issubset(set(ids)):
                    raise Stop(f'chr{chrom}: PSAM is missing {len(required - set(ids))} requested samples')
                # The chosen cohort includes relatives. Keep all requested
                # samples in LD, rather than silently dropping nonfounders.
                # Current 3202-person PSAM has no pedigree columns.
                if any(r.get(k, '0') not in ('0', 'NA', '.') for r in rows for k in ('PAT', 'MAT')):
                    raise Stop('Source PSAM has parent IDs; explicitly normalize pedigree before using this all-sample workflow')
                fid = {r['IID']: r.get('FID', '0') for r in rows}
                for pop in POPS:
                    atomic(self.out / 'meta' / f'{pop}.chr{chrom}.keep', '#FID\tIID\n' +
                           ''.join(f'{fid[s]}\t{s}\n' for s in self.samples[pop]))
            else:
                vcf = self.vcf_dir / f'chr{chrom}.vcf.gz'
                if not vcf.is_file() or not Path(str(vcf) + '.tbi').is_file():
                    raise Stop(f'Missing VCF or .tbi: {vcf}')
        note(f'Input={self.input_format}; build=GRCh{self.genome_build}; sample counts={self.sample_counts}')
        self.bin.parent.mkdir(parents=True, exist_ok=True)
        build_receipt = self.bin.with_suffix('.json')
        spec = {'source_sha256': self.source_hash, 'compiler': self.cxx, 'flags': ['-O3', '-std=c++17', '-lzstd']}
        cached = read_json(build_receipt)
        if not (self.bin.is_file() and cached and cached.get('spec') == spec and cached.get('binary_sha256') == sha(self.bin)):
            temporary = self.bin.with_suffix('.tmp')
            self.run_command([self.cxx, '-O3', '-std=c++17', HERE / 'ld_partition.cpp', '-lzstd', '-o', temporary])
            temporary.replace(self.bin)
            write_json(build_receipt, {'spec': spec, 'binary_sha256': sha(self.bin)})
        note(f'D grid (kb): {self.grid}; exact-boundary confirmations: {self.confirmations}')
        note(f'Reference-count comparison: policy={self.reference_policy}; '
             f'threshold={self.reference_tolerance:.2%} (engineering alert, not a statistical test)')
        note(f'Temporary LD directory: {self.tmp_parent}; storage={self.ld_storage}; explicitly preserved checkpoints are retained')

    def pgen_spec(self, pop, chrom):
        record = {'schema': 3, 'phase_qc_version': 1, 'pop': pop, 'chrom': chrom,
                'genome_build': f'GRCh{self.genome_build}', 'input_format': self.input_format,
                'sample_count': self.sample_counts[pop],
                'panel_sha256': sha(self.panel), 'samples_sha256': hashlib.sha256('\n'.join(self.samples[pop]).encode()).hexdigest(),
                'filter': {'maf': 0.05, 'min_alleles': 2, 'max_alleles': 2, 'snps_only': False,
                           'double_id': True, 'variant_ids': '@:#:$r:$a;duplicate_suffix=:recordN', 'preserve_phase': True},
                'plink2': self.tool}
        if self.input_format == 'pgen':
            record['source_pfiles'] = {k: stamp(p) for k, p in self.input_pfiles(chrom).items()}
            record['filter']['double_id'] = False
        else:
            record['vcf'] = stamp(self.vcf_dir / f'chr{chrom}.vcf.gz')
            if self.genome_build == '37':
                record['schema'] = 2
                for key in ('genome_build', 'input_format', 'sample_count'):
                    del record[key]
        return record

    def legacy_pgen_ok(self, prefix, pop, chrom):
        if self.input_format != 'vcf' or self.genome_build != '37':
            return False
        old = read_json(Path(str(prefix) + '.READY.json'))
        if not old or old.get('kind') != 'pgen' or old.get('implementation') not in LEGACY_PGEN_IMPLEMENTATIONS:
            return False
        if old.get('pop') != pop or str(old.get('chr')) != str(chrom):
            return False
        if old.get('prefix') != str(prefix) or old.get('vcf') != str(self.vcf_dir / f'chr{chrom}.vcf.gz') or old.get('panel') != str(self.panel):
            return False
        # A changed tool is allowed only after a fresh actual genotype/phase QC;
        # preprocessing itself is the exact known legacy specification above.
        try:
            required = [str(prefix) + s for s in ('.pgen', '.pvar', '.psam')]
            required += [old['vcf'], old['panel']]
            for name in required:
                st = Path(name).stat()
                if old.get('files', {}).get(name) != [st.st_size, st.st_mtime_ns]:
                    return False
        except (OSError, KeyError):
            return False
        return True

    def pgen(self, pop, chrom, *, allow_rebuild=True):
        key = (pop, chrom)
        if key in self.pgen_cache:
            return self.pgen_cache[key]
        prefix = self.out / 'pgen' / pop / f'chr{chrom}'
        files = [Path(str(prefix) + suffix) for suffix in ('.pgen', '.pvar', '.psam')]
        ready = Path(str(prefix) + '.READY.v2.json')
        spec = self.pgen_spec(pop, chrom)
        record = read_json(ready)
        valid = False
        try:
            valid = bool(record and record.get('spec') == spec and record.get('files') == [stamp(p) for p in files]
                         and record.get('qc', {}).get('missing_genotypes') == 0 and record.get('qc', {}).get('phase_check_passed'))
        except OSError:
            pass
        if valid:
            if record['qc'].get('extra_colocated_records', 0) and self.coincident_policy == 'error':
                raise Stop(f'{pop} chr{chrom} has co-located records; COINCIDENT_POLICY=error')
            note(f'REUSE PGEN: {pop} chr{chrom}')
            self.pgen_cache[key] = (prefix, record)
            return prefix, record
        if not allow_rebuild:
            return None
        prefix.parent.mkdir(parents=True, exist_ok=True)
        legacy = self.legacy_pgen_ok(prefix, pop, chrom)
        if legacy:
            note(f'Reuse verified legacy PGEN; add fresh phase/missingness QC: {pop} chr{chrom}')
        else:
            with self.temporary(f'{pop}.chr{chrom}.pgen') as temp:
                dest = temp / 'import'
                if self.input_format == 'pgen':
                    source_args = [arg for k, v in spec['source_pfiles'].items() for arg in (f'--{k}', v['path'])]
                else:
                    source_args = ['--vcf', spec['vcf']['path'], '--double-id']
                keep = f'{pop}.chr{chrom}.keep' if self.input_format == 'pgen' else f'{pop}.keep'
                self.run_command([self.plink, *source_args, '--keep', self.out / 'meta' / keep,
                                  '--maf', '0.05', '--min-alleles', '2', '--max-alleles', '2',
                                  '--new-id-max-allele-len', '10000', '--set-all-var-ids', '@:#:$r:$a',
                                  '--make-pgen', '--threads', self.threads, '--memory', self.memory_mb, '--out', dest], temp)
                renamed = uniquify_pvar_ids(Path(str(dest) + '.pvar'))
                note(f'{pop} chr{chrom}: renamed {renamed} repeated variant IDs; all records retained')
                for suffix, final in zip(('.pgen', '.pvar', '.psam'), files):
                    shutil.move(str(dest) + suffix, final)
        samples = list(table(str(prefix) + '.psam'))
        actual = sorted(r['IID'] for r in samples)
        if actual != self.samples[pop]:
            raise Stop(f'PGEN sample identity/count mismatch: {pop} chr{chrom}')
        positions = []
        ids = set()
        duplicates = 0
        with Path(str(prefix) + '.pvar').open() as handle:
            for line in handle:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.split()
                if len(fields) < 5 or fields[0].removeprefix('chr') != str(chrom) or ',' in fields[4]:
                    raise Stop('PVAR must contain sorted biallelic variants on exactly one chromosome')
                pos = int(fields[1])
                if pos < 1 or (positions and pos < positions[-1]) or fields[2] in ids:
                    raise Stop(f'PVAR invalid coordinate/order or duplicate ID: {pop} chr{chrom}, position={pos}, ID={fields[2]}, previous_position={positions[-1] if positions else None}')
                if positions and pos == positions[-1]:
                    duplicates += 1
                ids.add(fields[2]); positions.append(pos)
        if not positions:
            raise Stop(f'No retained variants: {pop} chr{chrom}')
        if duplicates and self.coincident_policy == 'error':
            raise Stop(f'{pop} chr{chrom} has {duplicates} extra records at repeated positions; COINCIDENT_POLICY=error')
        with self.temporary(f'{pop}.chr{chrom}.genotype_qc') as temp:
            dest = temp / 'check'
            # A 2-variant window validates phase with linear-scale work. These
            # prune lists are discarded: NO variants are actually pruned.
            self.run_command([self.plink, '--pfile', prefix, '--missing', 'variant-only',
                              '--indep-pairphase', '2', '1', '0.999999', '--threads', self.threads,
                              '--memory', self.memory_mb, '--out', dest], temp)
            missing_rows = list(table(str(dest) + '.vmiss'))
            if len(missing_rows) != len(positions):
                raise Stop('Missingness report does not cover every retained variant')
            if any(int(r['MISSING_CT']) != 0 or int(r['OBS_CT']) != self.sample_counts[pop] for r in missing_rows):
                raise Stop('Fixed haplotype N requires complete genotypes; missing genotypes detected')
        record = {'spec': spec, 'files': [stamp(p) for p in files],
                  'qc': {'missing_genotypes': 0, 'phase_check_passed': True, 'variants': len(positions),
                         'first_position': positions[0], 'last_position': positions[-1],
                         'extra_colocated_records': duplicates, 'legacy_pgen_adopted': legacy}}
        write_json(ready, record)
        note(f'PGEN QC: {pop} chr{chrom}, N={self.sample_counts[pop]} people, {len(positions)} variants, {duplicates} extra co-located records')
        self.pgen_cache[key] = (prefix, record)
        return prefix, record

    def partition_prefix(self, pop, chrom, d):
        return self.work / 'pop_blocks' / pop / 'strict' / f'D{d}kb' / f'chr{chrom}'

    def partition_spec(self, pop, chrom, d, pgen_record):
        return {'schema': 2, 'validation_version': 1, 'pop': pop, 'chrom': chrom, 'D_kb': d, 'source_sha256': self.source_hash,
                'pgen_receipt_sha256': hashlib.sha256(json.dumps(pgen_record, sort_keys=True).encode()).hexdigest(),
                'ld_mode': 'strict', 'penalty': 1 / (2 * self.sample_counts[pop] - 1) + 0.0025,
                'coincident_policy': self.coincident_policy, 'plink2': self.tool,
                'ld_statistic': 'phased_hardcall_r2', 'ld_pair_threshold': -1, 'numeric_id_mapping': 'v{zero_based_pvar_row}'}

    def partition_valid(self, pop, chrom, d, pgen_record):
        prefix = self.partition_prefix(pop, chrom, d)
        ready = read_json(str(prefix) + '.PARTITION_READY.json')
        spec = self.partition_spec(pop, chrom, d, pgen_record)
        try:
            if not ready or not isinstance(ready.get('spec'), dict):
                return False
            producer = ready['spec'].get('source_sha256')
            if producer in PARTITION_EQUIVALENT_SOURCES.get(self.source_hash, set()):
                spec = dict(spec, source_sha256=producer)
            return bool(ready.get('spec') == spec and ready.get('outputs') ==
                        {s: sha(str(prefix) + s) for s in ('.bed', '.qc.tsv')})
        except OSError:
            return False

    def completed_ld(self, pop, chrom, d, compact, expected_pairs):
        """Read an explicitly preserved PLINK table; never trust a partial table.

        The checkpoint is outside auto-cleaned temporary directories. On failure
        it remains available; only successfully validated partitions consume it.
        The C++ reader still verifies every ordered within-D pair before output.
        """
        if self.resume_dir is None:
            return None
        ready = read_json(self.resume_dir / 'LD_READY.json')
        if not ready or ready.get('schema') != 1:
            raise Stop(f'Invalid LD_RESUME_DIR receipt: {self.resume_dir}')
        if ready.get('output_root') != str(self.out):
            raise Stop('LD_RESUME_DIR belongs to a different output root')
        if ready.get('state') == 'consumed':
            return None
        if ready.get('state') != 'plink_completed_unpartitioned':
            raise Stop('LD_RESUME_DIR does not contain completed PLINK output')
        if (ready.get('pop'), ready.get('chrom')) != (pop, chrom):
            return None
        receipt = self.out / 'pgen' / pop / f'chr{chrom}.READY.v2.json'
        expected = {'D_kb': d, 'expected_pairs': expected_pairs,
                    'plink_binary_sha256': self.tool['sha256'],
                    'pgen_receipt_sha256': sha(receipt),
                    'numeric_pvar_sha256': sha(compact),
                    'ld_statistic': 'phased_hardcall_r2', 'ld_pair_threshold': -1,
                    'ld_window': 1000000000, 'numeric_id_mapping': 'v{zero_based_pvar_row}'}
        if any(ready.get(k) != value for k, value in expected.items()):
            raise Stop('LD_RESUME_DIR inputs or LD settings do not match this run')
        vcor = self.resume_dir / 'ld.vcor'
        if vcor.is_symlink() or not vcor.is_file():
            raise Stop('LD_RESUME_DIR has no regular ld.vcor file')
        recorded = ready.get('table', {})
        stat = vcor.stat()
        if stat.st_size != recorded.get('size') or stat.st_mtime_ns != recorded.get('mtime_ns'):
            raise Stop('Preserved LD table size or modification time changed')
        with vcor.open('rb') as handle:
            head = handle.read(65536)
            handle.seek(max(0, stat.st_size - 65536))
            tail = handle.read(65536)
        if (hashlib.sha256(head).hexdigest() != recorded.get('head_sha256') or
                hashlib.sha256(tail).hexdigest() != recorded.get('tail_sha256')):
            raise Stop('Preserved LD table head/tail integrity check failed')
        self.check_free_space(self.resume_dir)
        note(f'REUSE completed LD table: {pop} chr{chrom} D={d} kb ({stat.st_size / 1e9:.3f} GB); '
             'skip PLINK; strict pair validation or verified successful child adoption required')
        return vcor

    def compressed_ld(self, command, temp):
        """Supervise both ends of a local FIFO; publish only after both exit 0.

        The FIFO lives on Linux /tmp (DrvFS need not support FIFOs). Its contents
        never become a large local file. Compressed frames go directly to temp.
        """
        packed = temp / 'ld.vcor.ldz'
        partial = temp / 'ld.vcor.ldz.partial'
        processes = []
        started = time.monotonic()
        with tempfile.TemporaryDirectory(prefix='ld-stream-') as folder, (temp / 'compressor.log').open('wb') as error_handle:
            prefix = Path(folder) / 'ld'
            fifo = Path(str(prefix) + '.vcor')
            os.mkfifo(fifo, 0o600)
            producer_finished = None
            last_space_check = 0.0
            try:
                packer = subprocess.Popen([str(self.bin), '--pack', str(fifo), str(partial)],
                                          process_group=0, stderr=error_handle)
                processes.append(packer)
                args = list(map(str, command)) + ['--out', str(prefix)]
                note(shlex.join(args) + ' [lossless indexed Zstandard stream]')
                producer = subprocess.Popen(args, process_group=0)
                processes.append(producer)
                self.child = producer
                while True:
                    pc, cc = producer.poll(), packer.poll()
                    if pc not in (None, 0) or cc not in (None, 0):
                        raise Stop(f'LD stream failed: PLINK exit={pc}, compressor exit={cc}')
                    if pc == 0:
                        if not stat.S_ISFIFO(fifo.lstat().st_mode):
                            raise Stop('PLINK replaced the LD FIFO; refusing an unverified stream')
                        if cc == 0:
                            break
                        if producer_finished is None:
                            producer_finished = time.monotonic()
                        if time.monotonic() - producer_finished > 30:
                            raise Stop('LD compressor did not finish after PLINK success')
                    now = time.monotonic()
                    if now - last_space_check >= 3:
                        self.check_free_space(temp)
                        last_space_check = now
                    time.sleep(0.25)
                # fsync before making the completed container visible to the DP.
                with partial.open('rb') as handle:
                    os.fsync(handle.fileno())
                partial.replace(packed)
                return packed
            finally:
                for process in reversed(processes):
                    if process.poll() is None:
                        with contextlib.suppress(ProcessLookupError):
                            os.killpg(process.pid, signal.SIGTERM)
                for process in processes:
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        with contextlib.suppress(ProcessLookupError):
                            os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                self.child = None
                error_handle.flush()
                self.show_error_log(temp / 'compressor.log')
                note(f'LD stream elapsed: {time.monotonic() - started:.1f} s')
                for log in Path(folder).glob('*.log'):
                    shutil.copy2(log, temp / log.name)
                partial.unlink(missing_ok=True)

    def completed_partitions(self, pop, chrom, distances, record):
        """Adopt a successful gated legacy child, preserving its producer hash."""
        if self.resume_dir is None:
            return None
        receipt = read_json(self.resume_dir / 'PARTITIONS_READY.json')
        if not receipt:
            return None
        producer = receipt.get('source_sha256')
        allowed = {self.source_hash} | PARTITION_EQUIVALENT_SOURCES.get(self.source_hash, set())
        if (receipt.get('schema') != 1 or receipt.get('exit_code') != 0 or
                producer not in allowed or receipt.get('pop') != pop or receipt.get('chrom') != chrom or
                receipt.get('ld_receipt_sha256') != sha(self.resume_dir / 'LD_READY.json')):
            raise Stop('Preserved partition receipt does not match the completed LD checkpoint')
        dest = self.resume_dir / 'partitions'
        for d in distances:
            spec = dict(self.partition_spec(pop, chrom, d, record), source_sha256=producer)
            item = receipt.get('partitions', {}).get(str(d), {})
            if item.get('spec') != spec:
                raise Stop('Preserved partition specification mismatch')
            for suffix in ('.bed', '.qc.tsv'):
                path = dest / f'D{d}kb' / f'chr{chrom}{suffix}'
                if path.is_symlink() or item.get('outputs', {}).get(suffix) != sha(path):
                    raise Stop('Preserved partition output checksum mismatch')
        return dest, producer

    def partitions(self, pop, chrom, distances):
        prefix, record = self.pgen(pop, chrom)
        missing = [d for d in distances if not self.partition_valid(pop, chrom, d, record)]
        if not missing:
            note(f'REUSE validated partitions: {pop} chr{chrom} D={distances} kb')
            return
        max_d = max(missing)
        with self.temporary(f'{pop}.chr{chrom}.D{max_d}kb') as temp:
            # Numeric local IDs only reduce text volume. PGEN record order,
            # coordinates, REF/ALT and genotype bits remain identical.
            compact = temp / 'numeric.pvar'
            positions = []
            with Path(str(prefix) + '.pvar').open() as src, compact.open('w') as dst:
                dst.write('#CHROM\tPOS\tID\tREF\tALT\n')
                for line in src:
                    if line.startswith('#') or not line.strip():
                        continue
                    row = line.split(); row[2] = f'v{len(positions)}'
                    dst.write('\t'.join(row[:5]) + '\n'); positions.append(int(row[1]))
            expected_pairs, minimum_bytes = ld_table_requirements(positions, max_d * 1000)
            note(f'{pop} chr{chrom}: {expected_pairs:,} required LD pairs at D={max_d} kb; one temporary table, one reverse read')
            vcor = self.completed_ld(pop, chrom, max_d, compact, expected_pairs)
            resumed = vcor is not None
            if not resumed:
                free = self.check_free_space(temp)
                if self.ld_storage == 'plain':
                    note(f'LD disk preflight: plain .vcor needs at least {minimum_bytes / 1e9:.3f} GB '
                         f'(actual size is larger); {free / 1e9:.3f} GB free in {temp}; '
                         f'reserve={self.min_free / 1e9:g} GB')
                    if free < minimum_bytes + self.min_free:
                        raise Stop(f'Insufficient space for {pop} chr{chrom} D={max_d} kb in {temp}: '
                                   f'plain .vcor alone needs at least {minimum_bytes / 1e9:.3f} GB plus '
                                   f'MIN_FREE_GB={self.min_free / 1e9:g}, but only {free / 1e9:.3f} GB is free. '
                                   'PLINK2 LD calculation was not started. Set LD_TMP_DIR to a filesystem '
                                   'with more free space and rerun; validated results will be reused.')
                ldp = temp / 'ld'
                command = [self.plink, '--pgen', str(prefix) + '.pgen', '--pvar', compact, '--psam', str(prefix) + '.psam',
                           '--r2-phased', 'yes-really', 'cols=id', '--ld-window', '1000000000',
                           '--ld-window-kb', max_d, '--ld-window-r2', '-1', '--threads', self.threads,
                           '--memory', self.memory_mb]
                if self.ld_storage == 'blocked-zstd':
                    vcor = self.compressed_ld(command, temp)
                else:
                    self.run_command(command + ['--out', ldp], temp)
                    vcor = Path(str(ldp) + '.vcor')
            if not vcor.is_file():
                raise Stop('PLINK2 did not produce the requested completed LD table')
            note(f'Temporary LD table size: {vcor.stat().st_size / 1e9:.3f} GB')
            adopted = self.completed_partitions(pop, chrom, missing, record) if resumed else None
            producer_hash = self.source_hash
            if adopted:
                dest, producer_hash = adopted
                note(f'REUSE successful preserved partition child: {pop} chr{chrom}')
            else:
                dest = temp / 'partitions'
                self.run_command([self.bin, '--pvar', compact, '--vcor', vcor, '--chrom', chrom,
                                  '--penalty', format(1 / (2 * self.sample_counts[pop] - 1) + 0.0025, '.17g'),
                                  '--ld-mode', 'strict', '--coincident-policy', self.coincident_policy,
                                  '--distance-grid-kb', ','.join(map(str, missing)), '--out-dir', dest], temp,
                                 stderr_path=temp / 'partition.log')
            for d in missing:
                produced = dest / f'D{d}kb' / f'chr{chrom}'
                blocks = read_bed(str(produced) + '.bed', chrom)
                rows = list(table(str(produced) + '.qc.tsv'))
                if len(rows) != 1:
                    raise Stop('Expected exactly one partition QC row')
                qc = rows[0]
                if (int(qc['variants']) != record['qc']['variants'] or int(qc['blocks']) != len(blocks) or
                    int(qc['missing_pairs_within_D']) != 0 or qc['ld_mode'] != 'strict' or
                    int(qc['max_distance_bp']) != d * 1000 or not math.isfinite(float(qc['objective']))):
                    raise Stop(f'Partition QC failed for {pop} chr{chrom} D={d}')
                # Verify a complete, ordered cover of PVAR records, respecting
                # repeated positions. Do not silently drop singletons or short blocks.
                i = 0
                for start, end in blocks:
                    if i == len(positions) or positions[i] != start:
                        raise Stop('Partition omits or splits the first variant of a block')
                    while i < len(positions) and positions[i] <= end:
                        i += 1
                    if positions[i-1] != end:
                        raise Stop('Partition ends outside its assigned variants')
                if i != len(positions):
                    raise Stop('Partition does not cover every retained variant')
                final = self.partition_prefix(pop, chrom, d)
                final.parent.mkdir(parents=True, exist_ok=True)
                for suffix in ('.bed', '.qc.tsv'):
                    staged = Path(str(final) + suffix + f'.{os.getpid()}.tmp')
                    shutil.copy2(str(produced) + suffix, staged)
                    with staged.open('rb') as handle:
                        os.fsync(handle.fileno())
                    staged.replace(str(final) + suffix)
                write_json(str(final) + '.PARTITION_READY.json',
                           {'spec': dict(self.partition_spec(pop, chrom, d, record), source_sha256=producer_hash),
                            'outputs': {s: sha(str(final) + s) for s in ('.bed', '.qc.tsv')}})
            # Explicitly archived tables remain intact even after adoption.
            preserve = resumed and read_json(self.resume_dir / 'LD_READY.json').get('preserve_table', False)
            if not preserve:
                vcor.unlink()
            if resumed:
                ready = read_json(self.resume_dir / 'LD_READY.json')
                write_json(self.resume_dir / 'LD_READY.json', dict(ready, state='consumed'))
            note(f'{"Retained preserved" if preserve else "Deleted temporary"} LD table for {pop} chr{chrom}')

    def summarize(self, chrom, distances, *, verify_only=False):
        snapshots = {}
        previous = None
        lines = ['D_kb\tAFR\tEAS\tEUR\tSAS\tuniversal\tcounts_match\tsame_boundaries_as_previous\tchanged_population_intervals\tchanged_universal_intervals\tmax_reference_relative_deviation\tfull_span_all_populations\n']
        for d in distances:
            pop_blocks = {p: read_bed(str(self.partition_prefix(p, chrom, d)) + '.bed', chrom) for p in POPS}
            result = compute_consensus(pop_blocks, beta=2000)
            if not result['bed']:
                raise Stop('Empty universal consensus')
            dest = self.work / 'consensus' / 'strict' / f'D{d}kb'
            outputs = {
                dest / f'chr{chrom}.bed': ''.join(f'{chrom}\t{s}\t{e}\n' for s, e in result['bed']),
                dest / f'chr{chrom}.cuts.tsv': 'start_1based\tend_1based\n' + ''.join(f'{s}\t{e}\n' for s, e in result['cuts']),
            }
            for path, content in outputs.items():
                if verify_only:
                    if path.read_text() != content:
                        raise Stop(f'Saved consensus differs from validated partitions: {path}')
                else:
                    atomic(path, content)
            diagnostic_file = dest / f'chr{chrom}.diagnostics.json'
            if verify_only:
                if read_json(diagnostic_file) != result['diagnostics']:
                    raise Stop(f'Saved consensus diagnostics differ: {diagnostic_file}')
            else:
                write_json(diagnostic_file, result['diagnostics'])
            counts = [len(pop_blocks[p]) for p in POPS] + [len(result['bed'])]
            expected = [self.targets[(p, chrom)] for p in (*POPS, 'UNIVERSAL')] if self.targets else None
            deviation = max(abs(n - ref) / ref for n, ref in zip(counts, expected)) if expected else None
            signature = [pop_blocks[p] for p in POPS] + [result['bed'], result['cuts']]
            full = all(d * 1000 >= self.pgen_cache[(p, chrom)][1]['qc']['last_position'] -
                       self.pgen_cache[(p, chrom)][1]['qc']['first_position'] for p in POPS)
            snapshot = {'D_kb': d, 'signature': signature, 'full_span_all_populations': full,
                        'counts': dict(zip((*POPS, 'UNIVERSAL'), counts)), 'counts_match': counts == expected if expected else None,
                        'max_reference_relative_deviation': deviation,
                        'diagnostics': {'path': str(diagnostic_file.relative_to(self.work)),
                                        'sha256': sha(diagnostic_file)}}
            snapshots[d] = snapshot
            same = signature == previous['signature'] if previous else None
            change_pop = sum(changed_boundaries(previous['signature'][i], signature[i]) for i in range(4)) if previous else 'NA'
            change_uni = changed_boundaries(previous['signature'][4], signature[4]) if previous else 'NA'
            lines.append('\t'.join(map(str, [d, *counts, int(counts == expected) if expected else 'NA', int(same) if same is not None else 'NA',
                                             change_pop, change_uni, f'{deviation:.8g}' if deviation is not None else 'NA', int(full)])) + '\n')
            previous = snapshot
        if not verify_only:
            target = self.work / 'stability' / f'chr{chrom}.strict.summary.tsv'
            atomic(target, ''.join(lines))
            note(f'Boundary-stability report: {target}')
            print(''.join(lines), end='', flush=True)
        return snapshots

    def reference_count_check(self, chrom, snapshot, *, policy=None, tolerance=None):
        policy = self.reference_policy if policy is None else policy
        tolerance = self.reference_tolerance if tolerance is None else tolerance
        if policy == 'off':
            return {'policy': 'off', 'status': 'not_applicable',
                    'reason': 'Bundled GRCh37/Phase3 counts do not validate this build/cohort',
                    'comparisons': {}, 'is_statistical_test': False,
                    'author_coordinate_identity_verified': False}
        comparisons = {}
        for pop in (*POPS, 'UNIVERSAL'):
            observed, reference = snapshot['counts'][pop], self.targets[(pop, chrom)]
            relative = (observed - reference) / reference
            comparisons[pop] = {'observed': observed, 'reference': reference,
                                'signed_relative_deviation': relative,
                                'exceeds_threshold': abs(relative) > tolerance}
        exceeded = any(row['exceeds_threshold'] for row in comparisons.values())
        return {'policy': policy, 'threshold': tolerance,
                'status': ('warning' if policy == 'warn' else 'error') if exceeded else 'within_threshold',
                'comparisons': comparisons, 'is_statistical_test': False,
                'author_coordinate_identity_verified': False}

    def check_reference_counts(self, chrom, snapshot, *, context='selected', enforce=True):
        check = self.reference_count_check(chrom, snapshot)
        if check['status'] in ('within_threshold', 'not_applicable'):
            return check
        detail = ', '.join(f'{pop}={row["observed"]} vs {row["reference"]} '
                           f'({row["signed_relative_deviation"]:+.2%})'
                           for pop, row in check['comparisons'].items() if row['exceeds_threshold'])
        message = (f'chr{chrom} {context}: reference-count difference exceeds '
                   f'{self.reference_tolerance:.2%}: {detail}. '
                   'Count agreement does not establish coordinate identity; '
                   'investigate input/version/endpoint conventions.')
        if enforce and self.reference_policy == 'error':
            raise Stop(message + ' REFERENCE_COUNT_POLICY=error requires stopping.')
        note(f'WARNING: {message} '
             + ('REFERENCE_COUNT_POLICY=warn; retaining this difference in QC.' if enforce else
                f'Preflight only; current policy={self.reference_policy}; D stability remains unconfirmed.'))
        return check

    def reference_preflight(self, chrom, distances):
        """Report already validated cached counts before generating another LD table."""
        if self.reference_policy == 'off':
            return
        records = {}
        for pop in POPS:
            cached = self.pgen(pop, chrom, allow_rebuild=False)
            if cached is None:
                return
            records[pop] = cached[1]
        for d in reversed(distances):
            if all(self.partition_valid(pop, chrom, d, records[pop]) for pop in POPS):
                blocks = {p: read_bed(str(self.partition_prefix(p, chrom, d)) + '.bed', chrom) for p in POPS}
                consensus = compute_consensus(blocks, beta=2000)
                counts = {p: len(rows) for p, rows in blocks.items()}
                counts['UNIVERSAL'] = len(consensus['bed'])
                self.check_reference_counts(chrom, {'counts': counts},
                                            context=f'cached D={d} kb before new LD', enforce=False)
                return

    def selection_record(self, selected, distances, snapshots, *, reference_tolerance=None):
        tolerance = self.reference_tolerance if reference_tolerance is None else reference_tolerance
        record = {k: v for k, v in snapshots[selected].items() if k != 'signature'}
        full_span_verified = any(snapshots[d]['full_span_all_populations'] for d in distances if d >= selected)
        record.update({'tested_grid_kb': distances, 'empirical_stability': not full_span_verified,
                       'full_span_verified_in_suffix': full_span_verified,
                       'stable_confirmations_required': self.confirmations, 'reference_max_rel_deviation': tolerance,
                       'criterion': 'smallest D with identical full suffix; counts are only an independent sanity guard'})
        return record

    def audited_selection_record(self, chrom, selected, distances, snapshots, *, policy=None, tolerance=None):
        record = self.selection_record(selected, distances, snapshots, reference_tolerance=tolerance)
        record['reference_count_check'] = self.reference_count_check(chrom, snapshots[selected],
                                                                     policy=policy, tolerance=tolerance)
        return record

    def reuse_selection(self, chrom, distances):
        """Revalidate a passed chromosome on its recorded grid, without new LD.

        Extending the budget does not extend the evidence for a passed chromosome.
        Keep its original selected D and tested grid on reuse. Reference-count
        policy metadata may be refreshed after validating the old evidence.
        """
        path = self.work / 'selected' / f'chr{chrom}.json'
        if not path.exists():
            return None
        record = read_json(path)
        old_grid = record.get('tested_grid_kb') if isinstance(record, dict) else None
        old_tolerance = record.get('reference_max_rel_deviation') if isinstance(record, dict) else None
        old_check = record.get('reference_count_check') if isinstance(record, dict) else None
        old_policy = old_check.get('policy') if isinstance(old_check, dict) else 'error'
        if (not isinstance(old_grid, list) or not old_grid
                or any(type(d) is not int or d < 1 for d in old_grid)
                or old_grid != sorted(set(old_grid)) or distances[:len(old_grid)] != old_grid
                or record.get('D_kb') not in old_grid
                or record.get('stable_confirmations_required') != self.confirmations
                or type(old_tolerance) not in (int, float) or not math.isfinite(old_tolerance) or old_tolerance < 0
                or old_policy not in ('warn', 'error', 'off')
                or ('reference_count_check' in record and not isinstance(old_check, dict))):
            note(f'RECHECK chr{chrom}: saved selection is invalid or grid/validation settings changed')
            return None
        for pop in POPS:
            # This validates the current VCF, panel, samples, tool and PGEN stamps.
            # Real input/QC failures must propagate, rather than become cache misses.
            cached = self.pgen(pop, chrom, allow_rebuild=False)
            if cached is None:
                note(f'RECHECK chr{chrom}: {pop} PGEN no longer validates against current inputs')
                return None
            _, pgen_record = cached
            if not all(self.partition_valid(pop, chrom, d, pgen_record) for d in old_grid):
                note(f'RECHECK chr{chrom}: {pop} partition receipts or outputs no longer validate')
                return None
        try:
            snapshots = self.summarize(chrom, old_grid, verify_only=True)
            selected = choose_distance(snapshots, old_grid, self.confirmations)
            expected = (self.audited_selection_record(chrom, selected, old_grid, snapshots,
                                                      policy=old_policy, tolerance=old_tolerance)
                        if selected is not None else None)
            if expected is not None and old_check is None:
                expected.pop('reference_count_check')  # Legacy records predate explicit warn/error policy.
            if selected is None or record != expected:
                note(f'RECHECK chr{chrom}: saved selection no longer matches current boundary/count checks')
                return None
        except (Stop, OSError, ValueError, KeyError, TypeError) as exc:
            note(f'RECHECK chr{chrom}: saved boundary evidence is invalid: {exc}')
            return None
        # Policy is applied after scientific validation, outside the cache-miss
        # handler: an explicit error must never trigger more LD computation.
        self.check_reference_counts(chrom, snapshots[selected], context=f'reused D={selected} kb')
        current = self.audited_selection_record(chrom, selected, old_grid, snapshots)
        if current != record:
            write_json(path, current)
            record = current
            note(f'Refreshed reference-count audit for chr{chrom}; retained original D/grid and BED files')
        note(f'REUSE validated chromosome chr{chrom}: D={selected} kb; '
             f'previously tested through {old_grid[-1]} kb; no additional LD calculation')
        return selected, record

    def evaluate_chromosome(self, chrom, distances):
        reused = self.reuse_selection(chrom, distances)
        if reused is not None:
            return reused
        self.reference_preflight(chrom, distances)
        for pop in POPS:
            self.partitions(pop, chrom, distances)
        snapshots = self.summarize(chrom, distances)
        selected = choose_distance(snapshots, distances, self.confirmations)
        if selected is None:
            latest = snapshots[distances[-1]]
            suffix_start = len(distances) - 1
            while suffix_start > 0 and snapshots[distances[suffix_start - 1]]['signature'] == latest['signature']:
                suffix_start -= 1
            stable_suffix = distances[suffix_start:]
            reference_hint = ''
            if (latest['max_reference_relative_deviation'] is not None and
                    latest['max_reference_relative_deviation'] > self.reference_tolerance):
                failed_counts = ', '.join(
                    f'{pop}={count} (reference={self.targets[(pop, chrom)]})'
                    for pop, count in latest['counts'].items()
                    if abs(count - self.targets[(pop, chrom)]) / self.targets[(pop, chrom)] > self.reference_tolerance)
                reference_hint = (
                    f' Latest tested results also exceed REFERENCE_MAX_REL_DEVIATION={self.reference_tolerance:.2%}: '
                    f'{failed_counts}; maximum deviation={latest["max_reference_relative_deviation"]:.2%}. '
                    f'Current REFERENCE_COUNT_POLICY={self.reference_policy}: '
                    + ('this remains a warning after D stabilizes; ' if self.reference_policy == 'warn' else
                       'if these boundaries persist, the reference-count check will still stop the run; ')
                    + 'investigate input/version/endpoint conventions.')
            raise Stop(f'chr{chrom}: D stability is not yet confirmed across the tested grid {distances}. '
                       f'The identical trailing suffix is {stable_suffix}: '
                       f'{len(stable_suffix) - 1} of {self.confirmations} required larger-D confirmations '
                       f'for D={stable_suffix[0]} kb, with no full-span confirmation. '
                       f'No result published for this chromosome. Increase MAX_D_KB (currently {self.grid[-1]}) and rerun the SAME command. '
                       'Completed validated partitions will be reused; counts alone cannot pass.' + reference_hint)
        chosen = snapshots[selected]
        self.check_reference_counts(chrom, chosen, context=f'stable D={selected} kb')
        record = self.audited_selection_record(chrom, selected, distances, snapshots)
        write_json(self.work / 'selected' / f'chr{chrom}.json', record)
        note(f'SELECTED chr{chrom}: D={selected} kb; tested through {distances[-1]} kb; counts={chosen["counts"]}')
        return selected, record

    def publish(self, selections):
        all_text = []
        population_text = {p: [] for p in POPS}
        for chrom in self.chroms:
            distance = selections[chrom][0]
            path = self.work / 'consensus' / 'strict' / f'D{distance}kb' / f'chr{chrom}.bed'
            all_text.append(path.read_text())
            for pop in POPS:
                population_text[pop].append(Path(str(self.partition_prefix(pop, chrom, distance)) + '.bed').read_text())
        suffix = '' if self.chroms == list(range(1, 23)) else '.SUBSET'
        output = self.work / f'ALL.{self.genome_build}{suffix}.bed'
        atomic(output, ''.join(all_text))
        for pop in POPS:
            atomic(self.work / f'{pop}.{self.genome_build}{suffix}.bed', ''.join(population_text[pop]))
        warnings = [c for c in self.chroms
                    if selections[c][1]['reference_count_check']['status'] == 'warning']
        metadata = {'version': VERSION, 'genome_build': f'GRCh{self.genome_build}', 'chromosomes': self.chroms,
                    'input_format': self.input_format, 'sample_counts': self.sample_counts,
                    'panel_sha256': sha(self.panel),
                    'sample_selection': 'all panel samples in AFR/EAS/EUR/SAS; no relatedness pruning',
                    'status': 'complete_with_reference_warnings' if warnings else 'complete',
                    'reference_count_policy': self.reference_policy,
                    'reference_max_rel_deviation': self.reference_tolerance,
                    'reference_count_warning_chromosomes': warnings,
                    'independent_reimplementation': True, 'coordinate_identity_to_author_files_verified': False,
                    'distance_selection': {
                        str(c): {'D_kb': selections[c][0], 'path': f'selected/chr{c}.json',
                                 'sha256': sha(self.work / 'selected' / f'chr{c}.json')}
                        for c in self.chroms},
                    'universal_convention': 'right endpoints; nonempty C_p; zero-width cut base assigned left',
                    'coincident_policy': self.coincident_policy, 'beta_bp': 2000, 'plink2': self.tool,
                    'files': {p.name: sha(p) for p in [output] + [self.work / f'{pop}.{self.genome_build}{suffix}.bed' for pop in POPS]}}
        write_json(self.work / 'COMPLETE.json', metadata)
        if warnings:
            note(f'Completed with reference-count warnings on chromosomes {warnings}; '
                 'see selected/chrN.json and COMPLETE.json; author coordinate identity remains unverified')
        note(f'COMPLETE: {output}')
        return output

    def execute(self):
        with self.output_locks(), self.run_logging():
            self.work.mkdir(parents=True, exist_ok=True)
            identity_path = self.out / 'DATASET.json'
            identity = {'genome_build': f'GRCh{self.genome_build}', 'input_format': self.input_format,
                        'panel_sha256': sha(self.panel)}
            existing = read_json(identity_path)
            previous = read_json(self.work / 'COMPLETE.json')
            if existing is not None and existing != identity:
                raise Stop('OUT belongs to a different build/input format/sample panel; use a new output directory')
            if previous and previous.get('genome_build') != identity['genome_build']:
                raise Stop('OUT contains results from another genome build; use a new output directory')
            if (existing is None and self.genome_build == '38' and
                    (self.out / 'pgen').exists() and not previous):
                raise Stop('OUT has unidentified PGEN caches; use a fresh GRCh38 output directory')
            write_json(identity_path, identity)
            # A previous successful run's files are left intact; its completion
            # certificate is invalidated until THIS run has passed every gate.
            (self.work / 'COMPLETE.json').unlink(missing_ok=True)
            self.clean_abandoned()
            self.prepare()
            selections = self.evaluate_autosomes(self.grid)
            self.publish(selections)

    def chromosome_worker(self, chrom, distances):
        # Fork inherits the prepared tools/sample lists and the parent's output
        # locks. Never run prepare/clean_abandoned/publish in a worker: those
        # operations affect files shared by all chromosomes.
        self.chroms = [chrom]
        try:
            with self.run_logging():
                self.evaluate_chromosome(chrom, distances)
        except BaseException as exc:
            note(f'Worker chr{chrom} STOP: {type(exc).__name__}: {exc}')
            raise SystemExit(2)

    def evaluate_autosomes(self, distances):
        chromosomes = self.chroms
        if self.parallel_chroms == 1:
            return {c: self.evaluate_chromosome(c, distances) for c in chromosomes}
        # Linux/WSL only, like flock and the existing POSIX process supervision.
        context = multiprocessing.get_context('fork')
        pending = iter(chromosomes)
        active = {}
        results = {}
        note(f'Parallel chromosomes={self.parallel_chroms}; PLINK threads per worker={self.threads}')
        try:
            exhausted = False
            while active or not exhausted:
                while not exhausted and len(active) < self.parallel_chroms:
                    chrom = next(pending, None)
                    if chrom is None:
                        exhausted = True
                        break
                    sys.stdout.flush(); sys.stderr.flush()
                    process = context.Process(target=self.chromosome_worker, args=(chrom, distances))
                    process.start()
                    active[chrom] = process
                for chrom, process in list(active.items()):
                    if process.exitcode is None:
                        continue
                    process.join()
                    del active[chrom]
                    if process.exitcode != 0:
                        raise Stop(f'chr{chrom} worker failed (exit {process.exitcode}); see run logs. '
                                   'Other workers will stop; rerun to reuse validated results.')
                    record = read_json(self.work / 'selected' / f'chr{chrom}.json')
                    if not isinstance(record, dict) or record.get('D_kb') not in distances:
                        raise Stop(f'chr{chrom} worker did not produce a valid selection record')
                    results[chrom] = (record['D_kb'], record)
                if active:
                    time.sleep(0.2)
        finally:
            # The existing signal handler unwinds each worker's subprocess/FIFO
            # contexts, stops PLINK/packer process groups and cleans private temp.
            for process in active.values():
                if process.is_alive():
                    process.terminate()
            for process in active.values():
                process.join()
        return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('command', nargs='?', default='run', choices=['run', 'all'])
    parser.parse_args()
    def interrupted(signum, frame):
        # A repeated Ctrl+C must not interrupt child reaping or file deletion.
        for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            signal.signal(sig, signal.SIG_IGN)
        raise KeyboardInterrupt(f'Interrupted by signal {signum}; current temporary files will be removed')
    signal.signal(signal.SIGINT, interrupted)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGHUP, interrupted)
    try:
        Pipeline().execute()
        return 0
    except KeyboardInterrupt as exc:
        note(str(exc)); return 130
    except (Stop, OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        note(f'STOP: {exc}'); return 2

if __name__ == '__main__':
    sys.exit(main())
