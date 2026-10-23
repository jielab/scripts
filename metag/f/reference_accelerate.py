"""Prepare untouched genomes while preserving a live downloader in a frozen scope.

The original downloader and maskers remain in memory. Only files beyond its
bounded scheduling frontier are eligible. Original code later reuses the normal
receipts and publishes the library manifest in its original order. Run this in
a memory-limited user service with ExecStopPost=systemctl --user thaw UNIT.
"""
import argparse
from contextlib import contextmanager
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
from urllib.parse import urlsplit

from k2_compat import patched_source
from reference_batches import checksum, identity, prepare_one, SCHEMA
from reference_layout import source_path
from taxonomy import atomic_json

BIN = Path('/home/huangj/anaconda3/envs/metag/share/kraken2-2.17.1/libexec')


def genome_table(folder):
    frozen = folder / 'metag.assembly_summary.txt'
    if json.loads((folder / 'metag.assembly_summary.json').read_text()) != dict(schema=SCHEMA, sha256=checksum(frozen)):
        raise RuntimeError('Frozen assembly selection changed')
    table = {}
    for line in frozen.read_text().splitlines():
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        if not re.match('chromosome|complete genome', fields[11], re.I) or fields[19] == 'na':
            continue
        ftp = fields[19].rstrip('/')
        rel = urlsplit(ftp + '/' + ftp.rsplit('/', 1)[-1] + '_genomic.fna.gz').path.lstrip('/')
        table[rel] = dict(taxid=fields[5], url='https://ftp.ncbi.nlm.nih.gov/' + rel)
    manifest = (folder / 'manifest.txt').read_text().splitlines()
    if len(manifest) != len(set(manifest)) or set(table) != set(manifest):
        raise RuntimeError('Frozen assembly and running manifest differ')
    return table


def protected_frontier(files, log, original_jobs):
    positions = {Path(name).name: i for i, name in enumerate(files)}
    done = [positions[m[1]] for m in re.finditer(r'BATCH \d+/\d+: taxid=\d+; (?:completed|reused) (\S+);', log)
            if m[1] in positions]
    if not done:
        raise RuntimeError('Cannot establish the running downloader scheduling frontier')
    # At most original_jobs futures can be submitted beyond the highest logged
    # completion. Exclude that entire prefix, including any still-running file.
    return min(len(files), max(done) + 1 + original_jobs)


def validate_receipt(folder, rel, metadata, masker_sha):
    output = (folder / rel).with_suffix('')
    receipt = output.with_name(output.name + '.metag.json')
    if not receipt.exists():
        return False
    saved = json.loads(receipt.read_text())
    source = Path(source_path(str(folder / rel)))
    expected = dict(schema=SCHEMA, source=str(source), source_identity=identity(source),
                    taxid=metadata['taxid'], accessions={}, url=metadata['url'],
                    masker_sha256=masker_sha, masking='k2mask -outfmt fasta -r x')
    if (saved.get('parameters') != expected or identity(output) != saved['output_identity']
            or identity(output.with_name(output.name + '.prelim.tsv')) != saved['map_identity']):
        raise RuntimeError(f'Existing genome receipt failed validation: {rel}')
    return True


def unit_state(unit):
    return subprocess.check_output(['systemctl', '--user', 'show', unit, '-p', 'FreezerState', '--value'], text=True).strip()


@contextmanager
def frozen_unit(unit):
    if unit_state(unit) != 'running':
        raise RuntimeError('Original scope must be active and not already frozen')
    try:
        subprocess.run(['systemctl', '--user', 'freeze', unit], check=True, timeout=30)
        if unit_state(unit) != 'frozen':
            raise RuntimeError('Original scope did not fully freeze')
        yield
    finally:
        subprocess.run(['systemctl', '--user', 'thaw', unit], check=True, timeout=30)


def worker(folder, rel):
    table = genome_table(folder)
    metadata = table[rel]
    scope = dict(__name__='accelerator_k2', __file__=str(BIN / 'k2'))
    exec(compile(patched_source(BIN / 'k2'), str(BIN / 'k2'), 'exec'), scope)
    scope.update(LOG=logging.getLogger('masking'), SCRIPT_PATHNAME=str(BIN / 'k2'))
    os.chdir(folder)
    return prepare_one(rel, metadata['taxid'], {}, metadata['url'], BIN / 'k2mask',
                       checksum(BIN / 'k2mask'), 1, scope['scan_fasta_file'])


def accelerate(args):
    folder = args.db.resolve() / 'library' / args.library
    table = genome_table(folder)
    files = sorted(table)
    proc = Path('/proc') / str(args.owner_pid)
    command = (proc / 'cmdline').read_bytes().split(b'\0')
    if (str(args.db.resolve()).encode() not in command or b'download-library' not in command
            or args.library.encode() not in command or args.unit not in (proc/'cgroup').read_text()):
        raise RuntimeError('Original PID, command, database or scope does not match')
    if (b'--masker-threads' not in command or command[command.index(b'--masker-threads')+1] != b'1'):
        raise RuntimeError('Accelerator requires the same one-thread masker settings')
    env = (proc/'environ').read_bytes().split(b'\0')
    if f'METAG_MASK_JOBS={args.original_jobs}'.encode() not in env:
        raise RuntimeError('Original scheduling window does not match')
    args.report.mkdir(parents=True, exist_ok=True)
    masker_sha = checksum(BIN/'k2mask')
    # Validate completed progress before pausing or scheduling anything.
    before = {rel: (folder/rel).with_suffix('').with_name(Path(rel).with_suffix('').name+'.metag.json').read_bytes()
              for rel in files if validate_receipt(folder, rel, table[rel], masker_sha)}
    active = {}
    finished = []
    started = reported = time.monotonic()
    with frozen_unit(args.unit):
        try:
            frontier = protected_frontier(files, args.log.read_text(), args.original_jobs)
            todo = []
            for rel in files[frontier:]:
                out = (folder/rel).with_suffix('')
                if validate_receipt(folder, rel, table[rel], masker_sha):
                    continue
                if any(p.exists() for p in (out, out.with_name(out.name+'.metag.partial'),
                                             out.with_name(out.name+'.prelim.tsv'))):
                    raise RuntimeError(f'Unreceipted output beyond the protected frontier: {rel}')
                size = Path(source_path(str(folder/rel))).stat().st_size
                # Conservative working-set allowance; the service has its own
                # hard ceiling as an additional guard. Largest tasks start first.
                weight = max(512 << 20, size * 6 + (256 << 20))
                todo.append((rel, weight))
            todo.sort(key=lambda item: item[1], reverse=True)
            atomic_json(args.report/'plan.json', dict(owner_pid=args.owner_pid, unit=args.unit,
                preserved_completed=len(before), protected_prefix=frontier, jobs=args.jobs,
                memory_budget_gb=args.memory_budget_gb, pending=[rel for rel,_ in todo]))
            print(f'[accelerate] original scope frozen in memory; {len(before)} receipts preserved; '
                  f'{len(todo)} untouched genomes; up to {args.jobs} workers', flush=True)
            total = len(todo)
            while todo or active:
                if unit_state(args.unit) != 'frozen':
                    raise RuntimeError('Original scope was externally thawed; stopping auxiliary writers')
                for pid, (child, rel, weight) in list(active.items()):
                    code = child.poll()
                    if code is None:
                        continue
                    del active[pid]
                    if code or not validate_receipt(folder, rel, table[rel], masker_sha):
                        raise RuntimeError(f'Auxiliary masking failed: {rel}; see per-genome log')
                    finished.append(rel)
                    print(f'[accelerate] {len(finished)}/{total} completed: {Path(rel).name}', flush=True)
                used = sum(item[2] for item in active.values())
                while len(active) < args.jobs and todo:
                    available = args.memory_budget_gb * 2**30 - used
                    eligible = next((i for i, (_, weight) in enumerate(todo) if weight <= available), None)
                    if eligible is None:
                        if active:
                            break
                        raise RuntimeError('A genome exceeds the requested memory budget; original process will resume')
                    rel, weight = todo.pop(eligible)
                    with (args.report/(Path(rel).name+'.log')).open('w') as log:
                        child = subprocess.Popen([sys.executable, '-B', str(Path(__file__).resolve()),
                            'worker', '--folder', str(folder), '--genome', rel], stdout=log, stderr=subprocess.STDOUT)
                    active[child.pid] = (child, rel, weight)
                    used += weight
                if time.monotonic()-reported >= 60:
                    print(f'[accelerate] {len(finished)}/{total}; active={len(active)}; '
                          f'elapsed={(time.monotonic()-started)/60:.1f} min', flush=True)
                    reported = time.monotonic()
                time.sleep(2)
            for rel, contents in before.items():
                out = (folder/rel).with_suffix('')
                if out.with_name(out.name+'.metag.json').read_bytes() != contents:
                    raise RuntimeError(f'Previously completed receipt changed: {rel}')
                validate_receipt(folder, rel, table[rel], masker_sha)
            atomic_json(args.report/'completed.json', dict(preserved_completed=len(before),
                auxiliary_completed=len(finished), seconds=time.monotonic()-started))
        finally:
            for child, _, _ in active.values():
                if child.poll() is None:
                    child.terminate()
            for child, _, _ in active.values():
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
    print('[accelerate] original scope resumed; existing active genomes continue and new receipts will be reused', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    one = sub.add_parser('worker')
    one.add_argument('--folder', type=Path, required=True)
    one.add_argument('--genome', required=True)
    run = sub.add_parser('run')
    run.add_argument('--db', type=Path, required=True)
    run.add_argument('--library', default='invertebrate')
    run.add_argument('--owner-pid', type=int, required=True)
    run.add_argument('--unit', required=True)
    run.add_argument('--log', type=Path, required=True)
    run.add_argument('--report', type=Path, required=True)
    run.add_argument('--original-jobs', type=int, default=2)
    run.add_argument('--jobs', type=int, default=6)
    run.add_argument('--memory-budget-gb', type=int, default=12)
    args = parser.parse_args()
    if args.mode == 'worker':
        print(json.dumps(worker(args.folder.resolve(), args.genome)), flush=True)
        return
    if min(args.jobs, args.original_jobs, args.memory_budget_gb) < 1:
        parser.error('Worker counts and memory budget must be positive')
    def stop(signum, frame):
        raise KeyboardInterrupt('Auxiliary preparation cancelled; restoring original scope')
    signal.signal(signal.SIGTERM, stop)
    accelerate(args)


if __name__ == '__main__':
    main()
