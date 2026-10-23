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
        self.out = Path(os.environ.get('OUT', '/mnt/f/analysis/ld')).expanduser().resolve()
        self.work = self.out / 'rebuild_v2'
        self.tmp_parent = Path(os.environ.get('LD_TMP_DIR', str(self.out / 'tmp'))).expanduser().resolve()
        self.resume_dir = (Path(os.environ['LD_RESUME_DIR']).expanduser().resolve()
                           if os.environ.get('LD_RESUME_DIR') else None)
        self.vcf_dir = Path(os.environ.get('VCF_DIR', '/mnt/f/gen/1kg/37/vcf')).expanduser().resolve()
        self.panel = Path(os.environ.get('PANEL', '/mnt/f/gen/1kg/37/samples.txt')).expanduser().resolve()
        self.chroms = parse_chroms(os.environ.get('CHRS', '1-22'))
        self.threads = positive_int('THREADS', 8)
        self.memory_mb = positive_int('MEMORY_MB', 4096)
        self.confirmations = positive_int('STABLE_CONFIRMATIONS', 2)
        # At least two independent larger-D comparisons by default.
        raw_grid = os.environ.get('D_GRID_KB', '250,500,1000,2000,5000,10000')
        self.grid = [int(x) for x in raw_grid.split(',')]
        if not self.grid or self.grid != sorted(set(self.grid)) or min(self.grid) < 1:
            raise Stop('D_GRID_KB must be a strictly increasing list of positive integers')
        cap = positive_int('PILOT_MAX_D_KB', 20000)
        if self.grid[-1] > cap:
            raise Stop('D_GRID_KB exceeds PILOT_MAX_D_KB')
        if cap > (2**63 - 1) // 1000:
            raise Stop('PILOT_MAX_D_KB exceeds the supported integer range')
        while self.grid[-1] < cap:
            self.grid.append(min(2 * self.grid[-1], cap))
        if len(self.grid) > 64:
            raise Stop('D_GRID_KB has too many states (maximum 64)')
        self.min_free = float(os.environ.get('MIN_FREE_GB', '5')) * 10**9
        self.reference_tolerance = float(os.environ.get('REFERENCE_MAX_REL_DEVIATION', '0.10'))
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
        target_file = Path(os.environ.get('LD_TARGETS', str(HERE.parent / 'reference' / 'published_counts.tsv')))
        self.targets = load_targets(target_file)
        for c in sorted(set(self.chroms + [22])):
            for pop in (*POPS, 'UNIVERSAL'):
                if (pop, c) not in self.targets:
                    raise Stop(f'Missing reference count for {pop} chr{c}: {target_file}')
        plink = os.environ.get('PLINK2', 'plink2')
        if plink == 'plink2' and not shutil.which(plink) and Path('/mnt/d/software/bin/plink2').is_file():
            plink = '/mnt/d/software/bin/plink2'
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

    def run_command(self, command, disk_path=None):
        note(shlex.join(map(str, command)))
        # Keep the terminal session; a separate process group lets the parent
        # stop and reap the command and its descendants together on interruption.
        self.child = subprocess.Popen(list(map(str, command)), process_group=0)
        try:
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
            if self.child.poll() is None:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(self.child.pid, signal.SIGTERM)
                try:
                    self.child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(self.child.pid, signal.SIGKILL)
                    self.child.wait()
            self.child = None

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
    def output_locks(self):
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
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
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
            if len(self.samples[pop]) != EXPECTED[pop]:
                raise Stop(f'{pop}: expected {EXPECTED[pop]} people, found {len(self.samples[pop])}')
            self.samples[pop].sort()
            atomic(self.out / 'meta' / f'{pop}.keep', '#FID\tIID\n' + ''.join(f'{s}\t{s}\n' for s in self.samples[pop]))
        for chrom in sorted(set(self.chroms + [22])):
            vcf = self.vcf_dir / f'chr{chrom}.vcf.gz'
            if not vcf.is_file() or not Path(str(vcf) + '.tbi').is_file():
                raise Stop(f'Missing VCF or .tbi: {vcf}')
        self.bin.parent.mkdir(parents=True, exist_ok=True)
        build_receipt = self.bin.with_suffix('.json')
        spec = {'source_sha256': self.source_hash, 'compiler': self.cxx, 'flags': ['-O3', '-std=c++17', '-lzstd']}
        cached = read_json(build_receipt)
        if not (self.bin.is_file() and cached and cached.get('spec') == spec and cached.get('binary_sha256') == sha(self.bin)):
            temporary = self.bin.with_suffix('.tmp')
            self.run_command([self.cxx, '-O3', '-std=c++17', HERE / 'ld_partition.cpp', '-lzstd', '-o', temporary])
            temporary.replace(self.bin)
            write_json(build_receipt, {'spec': spec, 'binary_sha256': sha(self.bin)})
        note(f'Pilot grid (kb): {self.grid}; exact-boundary confirmations: {self.confirmations}')
        note(f'Temporary LD directory: {self.tmp_parent}; storage={self.ld_storage}; explicitly preserved checkpoints are retained')

    def pgen_spec(self, pop, chrom):
        return {'schema': 2, 'phase_qc_version': 1, 'pop': pop, 'chrom': chrom, 'vcf': stamp(self.vcf_dir / f'chr{chrom}.vcf.gz'),
                'panel_sha256': sha(self.panel), 'samples_sha256': hashlib.sha256('\n'.join(self.samples[pop]).encode()).hexdigest(),
                'filter': {'maf': 0.05, 'min_alleles': 2, 'max_alleles': 2, 'snps_only': False,
                           'double_id': True, 'variant_ids': '@:#:$r:$a', 'preserve_phase': True},
                'plink2': self.tool}

    def legacy_pgen_ok(self, prefix, pop, chrom):
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
                self.run_command([self.plink, '--vcf', spec['vcf']['path'], '--double-id', '--keep', self.out / 'meta' / f'{pop}.keep',
                                  '--maf', '0.05', '--min-alleles', '2', '--max-alleles', '2',
                                  '--new-id-max-allele-len', '10000', '--set-all-var-ids', '@:#:$r:$a',
                                  '--make-pgen', '--threads', self.threads, '--memory', self.memory_mb, '--out', dest], temp)
                for suffix, final in zip(('.pgen', '.pvar', '.psam'), files):
                    shutil.move(str(dest) + suffix, final)
        samples = list(table(str(prefix) + '.psam'))
        actual = sorted(r['IID'] for r in samples)
        if actual != self.samples[pop] or any(r.get('FID', r['IID']) != r['IID'] for r in samples):
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
                    raise Stop('PVAR has invalid coordinates, sorting or duplicate IDs')
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
            if any(int(r['MISSING_CT']) != 0 or int(r['OBS_CT']) != EXPECTED[pop] for r in missing_rows):
                raise Stop('Fixed haplotype N requires complete genotypes; missing genotypes detected')
        record = {'spec': spec, 'files': [stamp(p) for p in files],
                  'qc': {'missing_genotypes': 0, 'phase_check_passed': True, 'variants': len(positions),
                         'first_position': positions[0], 'last_position': positions[-1],
                         'extra_colocated_records': duplicates, 'legacy_pgen_adopted': legacy}}
        write_json(ready, record)
        note(f'PGEN QC: {pop} chr{chrom}, N={EXPECTED[pop]} people, {len(positions)} variants, {duplicates} extra co-located records')
        self.pgen_cache[key] = (prefix, record)
        return prefix, record

    def partition_prefix(self, pop, chrom, d):
        return self.work / 'pop_blocks' / pop / 'strict' / f'D{d}kb' / f'chr{chrom}'

    def partition_spec(self, pop, chrom, d, pgen_record):
        return {'schema': 2, 'validation_version': 1, 'pop': pop, 'chrom': chrom, 'D_kb': d, 'source_sha256': self.source_hash,
                'pgen_receipt_sha256': hashlib.sha256(json.dumps(pgen_record, sort_keys=True).encode()).hexdigest(),
                'ld_mode': 'strict', 'penalty': 1 / (2 * EXPECTED[pop] - 1) + 0.0025,
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
        with tempfile.TemporaryDirectory(prefix='ld-stream-') as folder:
            prefix = Path(folder) / 'ld'
            fifo = Path(str(prefix) + '.vcor')
            os.mkfifo(fifo, 0o600)
            producer_finished = None
            last_space_check = 0.0
            try:
                packer = subprocess.Popen([str(self.bin), '--pack', str(fifo), str(partial)], process_group=0)
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
                                  '--penalty', format(1 / (2 * EXPECTED[pop] - 1) + 0.0025, '.17g'),
                                  '--ld-mode', 'strict', '--coincident-policy', self.coincident_policy,
                                  '--distance-grid-kb', ','.join(map(str, missing)), '--out-dir', dest], temp)
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
            expected = [self.targets[(p, chrom)] for p in (*POPS, 'UNIVERSAL')]
            deviation = max(abs(n - ref) / ref for n, ref in zip(counts, expected))
            signature = [pop_blocks[p] for p in POPS] + [result['bed'], result['cuts']]
            full = all(d * 1000 >= self.pgen_cache[(p, chrom)][1]['qc']['last_position'] -
                       self.pgen_cache[(p, chrom)][1]['qc']['first_position'] for p in POPS)
            snapshot = {'D_kb': d, 'signature': signature, 'full_span_all_populations': full,
                        'counts': dict(zip((*POPS, 'UNIVERSAL'), counts)), 'counts_match': counts == expected,
                        'max_reference_relative_deviation': deviation,
                        'diagnostics': {'path': str(diagnostic_file.relative_to(self.work)),
                                        'sha256': sha(diagnostic_file)}}
            snapshots[d] = snapshot
            same = signature == previous['signature'] if previous else None
            change_pop = sum(changed_boundaries(previous['signature'][i], signature[i]) for i in range(4)) if previous else 'NA'
            change_uni = changed_boundaries(previous['signature'][4], signature[4]) if previous else 'NA'
            lines.append('\t'.join(map(str, [d, *counts, int(counts == expected), int(same) if same is not None else 'NA',
                                             change_pop, change_uni, f'{deviation:.8g}', int(full)])) + '\n')
            previous = snapshot
        if not verify_only:
            target = self.work / 'pilot' / f'chr{chrom}.strict.summary.tsv'
            atomic(target, ''.join(lines))
            note(f'Boundary-stability report: {target}')
            print(''.join(lines), end='', flush=True)
        return snapshots

    def selection_record(self, selected, distances, snapshots):
        record = {k: v for k, v in snapshots[selected].items() if k != 'signature'}
        full_span_verified = any(snapshots[d]['full_span_all_populations'] for d in distances if d >= selected)
        record.update({'tested_grid_kb': distances, 'empirical_stability': not full_span_verified,
                       'full_span_verified_in_suffix': full_span_verified,
                       'stable_confirmations_required': self.confirmations, 'reference_max_rel_deviation': self.reference_tolerance,
                       'criterion': 'smallest D with identical full suffix; counts are only an independent sanity guard'})
        return record

    def reuse_selection(self, chrom, distances):
        """Revalidate a passed chromosome on its recorded grid, without new LD.

        Extending the budget does not extend the evidence for a passed chromosome.
        Keep its original selection record and tested grid unchanged on reuse.
        """
        path = self.work / 'selected' / f'chr{chrom}.json'
        if not path.exists():
            return None
        record = read_json(path)
        old_grid = record.get('tested_grid_kb') if isinstance(record, dict) else None
        if (not isinstance(old_grid, list) or not old_grid
                or any(type(d) is not int or d < 1 for d in old_grid)
                or old_grid != sorted(set(old_grid)) or distances[:len(old_grid)] != old_grid
                or record.get('D_kb') not in old_grid
                or record.get('stable_confirmations_required') != self.confirmations
                or record.get('reference_max_rel_deviation') != self.reference_tolerance):
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
            if (selected is None
                    or snapshots[selected]['max_reference_relative_deviation'] > self.reference_tolerance
                    or record != self.selection_record(selected, old_grid, snapshots)):
                note(f'RECHECK chr{chrom}: saved selection no longer matches current boundary/count checks')
                return None
        except (Stop, OSError, ValueError, KeyError, TypeError) as exc:
            note(f'RECHECK chr{chrom}: saved boundary evidence is invalid: {exc}')
            return None
        note(f'REUSE validated chromosome chr{chrom}: D={selected} kb; '
             f'previously tested through {old_grid[-1]} kb; no additional LD calculation')
        return selected, record

    def evaluate_chromosome(self, chrom, distances):
        reused = self.reuse_selection(chrom, distances)
        if reused is not None:
            return reused
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
            if latest['max_reference_relative_deviation'] > self.reference_tolerance:
                failed_counts = ', '.join(
                    f'{pop}={count} (reference={self.targets[(pop, chrom)]})'
                    for pop, count in latest['counts'].items()
                    if abs(count - self.targets[(pop, chrom)]) / self.targets[(pop, chrom)] > self.reference_tolerance)
                reference_hint = (
                    f' Latest tested results also exceed REFERENCE_MAX_REL_DEVIATION={self.reference_tolerance:.2%}: '
                    f'{failed_counts}; maximum deviation={latest["max_reference_relative_deviation"]:.2%}. '
                    'If these boundaries persist, extending D alone will still fail the reference-count check; '
                    'investigate input/version/endpoint conventions.')
            raise Stop(f'chr{chrom}: D stability is not yet confirmed across the tested grid {distances}. '
                       f'The identical trailing suffix is {stable_suffix}: '
                       f'{len(stable_suffix) - 1} of {self.confirmations} required larger-D confirmations '
                       f'for D={stable_suffix[0]} kb, with no full-span confirmation. '
                       f'No whole-genome continuation. Increase PILOT_MAX_D_KB (currently {self.grid[-1]}) and rerun the SAME command. '
                       'Completed validated partitions will be reused; counts alone cannot pass.' + reference_hint)
        chosen = snapshots[selected]
        if chosen['max_reference_relative_deviation'] > self.reference_tolerance:
            raise Stop(f'chr{chrom}: stable boundaries, but maximum count deviation from the paper is '
                       f'{chosen["max_reference_relative_deviation"]:.1%} > {self.reference_tolerance:.1%}. '
                       'Investigate input/version/endpoint conventions; do not tune D to force a count match.')
        record = self.selection_record(selected, distances, snapshots)
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
        output = self.work / f'ALL.37{suffix}.bed'
        atomic(output, ''.join(all_text))
        for pop in POPS:
            atomic(self.work / f'{pop}.37{suffix}.bed', ''.join(population_text[pop]))
        metadata = {'version': VERSION, 'genome_build': 'GRCh37', 'chromosomes': self.chroms,
                    'independent_reimplementation': True, 'coordinate_identity_to_author_files_verified': False,
                    'distance_selection': {
                        str(c): {'D_kb': selections[c][0], 'path': f'selected/chr{c}.json',
                                 'sha256': sha(self.work / 'selected' / f'chr{c}.json')}
                        for c in self.chroms},
                    'universal_convention': 'right endpoints; nonempty C_p; zero-width cut base assigned left',
                    'coincident_policy': self.coincident_policy, 'beta_bp': 2000, 'plink2': self.tool,
                    'files': {p.name: sha(p) for p in [output] + [self.work / f'{pop}.37{suffix}.bed' for pop in POPS]}}
        write_json(self.work / 'COMPLETE.json', metadata)
        note(f'COMPLETE: {output}')
        return output

    def execute(self):
        with self.output_locks():
            self.work.mkdir(parents=True, exist_ok=True)
            # A previous successful run's files are left intact; its completion
            # certificate is invalidated until THIS run has passed every gate.
            (self.work / 'COMPLETE.json').unlink(missing_ok=True)
            self.clean_abandoned()
            self.prepare()
            if 22 in self.chroms:
                pilot = self.evaluate_chromosome(22, self.grid)
            else:
                # An explicit subset must never start work on an excluded chr22.
                pilot = self.reuse_selection(22, self.grid)
                if pilot is None:
                    raise Stop('The requested chromosome subset excludes chr22, but no reusable chr22 pilot is valid. '
                               'Include 22 in --chr/CHRS, or run ./ld.sh with the default chromosome list first.')
            pilot_d, pilot_record = pilot
            note(f'PILOT PASSED. Automatic initial D={pilot_d} kb. Continuing without a second command.')
            distances = [d for d in self.grid if d >= pilot_d]
            selections = {22: (pilot_d, pilot_record)}
            for chrom in self.chroms:
                if chrom == 22:
                    continue
                # chr22 convergence is not a guarantee for chr6 or another
                # chromosome. Verify every chromosome against larger D too.
                selections[chrom] = self.evaluate_chromosome(chrom, distances)
            self.publish(selections)


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
