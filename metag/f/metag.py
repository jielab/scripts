#!/usr/bin/env python3
"""Local ENA metagenome inventory, immutable raw organization, QC and species bins.

No sample merging, no host deletion, no per-read in-memory dictionary. See README.md in the project root.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
from datetime import datetime, timezone
import csv
import fcntl
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import uuid

from taxonomy import GROUPS, Taxonomy, audit_database, sha256, slug, snapshot, prepare_sqlite
from splitter import fastq, pair_ids, GzipBins, split_stream, ASSIGNMENT_FIELDS

VERSION = '2.0.0-20261014'
SCHEMA = 'v2'
EXPERIMENT_FIELDS = ['input_nucleic_acid', 'library_method', 'reverse_transcription',
    'dna_retained', 'virus_enriched', 'processing_protocol', 'is_control',
    'metadata_evidence', 'metadata_review_status']
REPORT_METADATA = EXPERIMENT_FIELDS + ['assay_class', 'library_strategy', 'library_source',
    'library_selection', 'human_read_masking', 'pooling', 'source_matrix']
PROVENANCE_FIELDS = ['schema', 'mapping_schema', 'mapping_sha256', 'mapping_build_signature']
ENA_FIELDS = ['run_accession', 'study_accession', 'sample_accession',
              'secondary_sample_accession', 'experiment_accession', 'sample_alias',
              'sample_title', 'experiment_title', 'library_layout', 'instrument_platform',
              'instrument_model', 'library_strategy', 'library_source', 'library_selection',
              'fastq_ftp', 'fastq_md5', 'fastq_bytes']
ENA_REQUIRED_FIELDS = ['run_accession', 'fastq_ftp', 'fastq_md5', 'fastq_bytes',
                       'sample_accession', 'study_accession', 'library_layout',
                       'library_strategy', 'library_source', 'instrument_platform']
SAMPLE_FIELDS = ['dataset', 'project', 'run_accession', 'experiment_accession',
                 'sample_accession', 'secondary_sample_accession', 'sample_key',
                 'mapping_status', 'sample_alias', 'sample_title', 'experiment_title',
                 'library_layout', 'instrument_platform', 'instrument_model',
                 'library_strategy', 'library_source', 'library_selection', 'assay_class',
                 'observed_layout', 'qc_policy', 'human_read_masking', 'include', 'valid',
                 'r1', 'r2', 'single', 'source_matrix', 'sampling_date', 'site', 'country',
                 'latitude', 'longitude', 'pooling', 'notes'] + EXPERIMENT_FIELDS
FILE_FIELDS = ['run_accession', 'role', 'source_path', 'raw_path', 'bytes',
               'ena_bytes', 'ena_md5', 'ena_ftp']
COUNTS_FIELDS = ['dataset', 'sample_key', 'run_accession', 'group', 'bucket',
                 'species_taxid', 'species_name', 'stream', 'fragments', 'reads', 'scope'] + REPORT_METADATA + PROVENANCE_FIELDS
VIRUS_TAXA_FIELDS = ['dataset', 'sample_key', 'run_accession', 'scope', 'stream',
    'original_assigned_taxid', 'canonical_taxid', 'assigned_rank', 'assigned_name',
    'species_taxid', 'species_name', 'group', 'genome_class', 'genome_type_detail',
    'resolution_status', 'evidence_taxid', 'fragments', 'reads'] + REPORT_METADATA + PROVENANCE_FIELDS
AUDIT_FIELDS = ['dataset', 'sample_key', 'run_accession', 'original_assigned_taxid',
    'canonical_taxid', 'assigned_rank', 'assigned_name', 'genome_class', 'genome_type_detail',
    'resolution_status', 'evidence_taxid', 'evidence_method', 'conflict_status',
    'source_ref', 'source_version', 'source_accession', 'evidence_json'] + PROVENANCE_FIELDS
FASTQ_NAME = re.compile(r'^([SED]RR\d+)(?:_([12]))?\.(?:fastq|fq)\.gz$')


def now():
    return datetime.now(timezone.utc).isoformat()


def json_write(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.partial')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')
    tmp.replace(path)


def write_tsv(path: Path, fields: list, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.partial')
    with tmp.open('w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def read_tsv(path: Path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def selected(rows: list, runs: str = '', limit: int = 0):
    wanted = set(filter(None, runs.split(',')))
    if wanted:
        missing = wanted - {r['run_accession'] for r in rows}
        if missing:
            raise RuntimeError(f'Unknown runs: {sorted(missing)}')
        rows = [r for r in rows if r['run_accession'] in wanted]
    return rows[:limit] if limit else rows


def lock(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    f = path.open('a+')
    try:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        f.close()
        raise RuntimeError(f'Another process holds {path}')
    return f


def parse_ena_report(text: str, source: str):
    """Reject HTTP-200 error bodies and partial TSVs before caching or inventorying."""
    reader = csv.DictReader(io.StringIO(text.lstrip('\ufeff')), delimiter='\t', strict=True)
    try:
        fields = reader.fieldnames or []
        missing = set(ENA_REQUIRED_FIELDS) - set(fields)
        if missing or len(set(fields)) != len(fields):
            raise RuntimeError(f'Invalid ENA report ({source}): missing or duplicate columns; '
                               f'request a full file report. Missing: {sorted(missing)}')
        rows = []
        for row in reader:
            runid = row.get('run_accession') or ''
            if not re.fullmatch(r'[SED]RR\d+', runid):
                raise RuntimeError(f'Invalid ENA report ({source}), line {reader.line_num}: '
                                   f'expected a run accession, got {runid[:160]!r}. '
                                   'The service may have returned an error or partial report.')
            if None in row or any(value is None for value in row.values()):
                raise RuntimeError(f'Invalid ENA report ({source}), line {reader.line_num}: '
                                   'incomplete or extra TSV columns; report may be truncated.')
            rows.append(row)
    except csv.Error as error:
        raise RuntimeError(f'Invalid ENA report ({source}): malformed TSV: {error}') from error
    if not rows:
        raise RuntimeError(f'ENA returned no usable read_run records ({source})')
    return rows


def fetch_ena(project: str, destination: Path):
    query = {'accession': project, 'result': 'read_run',
             'fields': ','.join(ENA_FIELDS), 'format': 'tsv', 'limit': '0'}
    url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urllib.parse.urlencode(query)
    last = None
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': f'metag-split/{VERSION}'})
            with urllib.request.urlopen(request, timeout=180) as response:
                data = response.read()
            rows = parse_ena_report(data.decode('utf-8-sig'), project)
            pending = destination.with_name(destination.name + '.partial')
            pending.write_bytes(data)
            pending.replace(destination)
            json_write(destination.with_suffix('.source.json'),
                       {'url': url, 'retrieved_utc': now(), 'sha256': hashlib.sha256(data).hexdigest()})
            return rows
        except Exception as error:
            last = error
            if attempt < 3:
                print(f'WARNING: ENA {project} attempt {attempt + 1}/4 failed: {error}; retrying.',
                      file=sys.stderr)
                time.sleep(2 ** attempt)
    raise RuntimeError(f'ENA retrieval failed: {last}. Use --ena-tsv with a saved full report.')


def project_policy(row: dict):
    project = row.get('study_accession', '')
    strategy = row.get('library_strategy', '').upper()
    selection = row.get('library_selection', '').upper()
    notes = []
    platform = row.get('instrument_platform', '').upper()
    qc = 'fastp' if platform == 'ILLUMINA' else 'minimal'
    assay = 'shotgun_DNA' if strategy == 'WGS' else ('RNA_seq' if strategy == 'RNA-SEQ' else 'review_required')
    masked = 'unknown'
    if row.get('library_source', '').upper() != 'METAGENOMIC':
        assay = 'review_required'
        notes.append('library_source is not METAGENOMIC; review isolate/host/other-source sequencing.')
    if strategy == 'AMPLICON' or any(s in selection for s in ('HYBRID', 'CAPTURE')):
        assay = 'targeted_or_review'
    if project == 'PRJEB1787':
        notes.append('README: DNA size-fractionated Tara subset; includes single-end 454; not total ocean life.')
    if project == 'PRJEB68319':
        qc = 'minimal'
        notes.append('README: archive FASTQ already quality-trimmed; disable new adapter auto-trimming.')
    if project in ('PRJEB22068', 'PRJNA406976', 'PRJNA406977'):
        notes.append('README: check processed/merged NEON reads; file layout is not necessarily original library layout.')
    if project == 'PRJNA1247874':
        masked = 'reported_masked'
        notes.append('README: CASPER RNA/cDNA workflow; human reads reported masked; zero human is not absence.')
    if project == 'PRJEB62830':
        assay = 'targeted_or_review'
        notes.append('README: London capture/amplicon/controls and pooled samples; no automatic ecological comparison.')
    if project == 'PRJEB34634':
        notes.append('README: Oxford sample mapping unreliable; isolate every run until an external sample map is validated.')
    return qc, assay, masked, notes


def experiment_defaults(row):
    values = {k: 'unknown' for k in EXPERIMENT_FIELDS}
    if row.get('study_accession', row.get('project')) == 'PRJEB1787' and row.get('library_strategy','').upper() == 'WGS':
        values.update(input_nucleic_acid='DNA', library_method='shotgun',
                      metadata_evidence='README.md#ocean-source; ENA library_strategy=WGS',
                      metadata_review_status='provisional')
    return values


def metadata_upgrade(args):
    meta = args.root.resolve() / 'metadata'
    path = meta / 'samples.tsv'
    with lock(meta / '.metadata.lock'):
        rows = read_tsv(path)
        with path.open(newline='') as f:
            fields = next(csv.reader(f, delimiter='\t'))
        fields += [k for k in EXPERIMENT_FIELDS if k not in fields]
        changes, conflicts = [], []
        updates = read_tsv(args.updates) if args.updates else []
        known = {r['run_accession'] for r in rows}
        if len(known) != len(rows):
            raise RuntimeError('Duplicate run in samples.tsv')
        patch = {}
        for r in updates:
            runid = r.get('run_accession','')
            if runid not in known or runid in patch:
                raise RuntimeError(f'Unknown/duplicate metadata update run: {runid}')
            if not r.get('metadata_evidence') or r.get('metadata_review_status') not in ('provisional','reviewed'):
                raise RuntimeError('Metadata updates require evidence and provisional/reviewed status')
            if set(r)-set(REPORT_METADATA)-{'run_accession'}:
                raise RuntimeError('Metadata updates cannot modify include/valid/path/sample identity')
            patch[runid] = r
            fields += [k for k in r if k not in fields]
        for row in rows:
            rid = row['run_accession']
            # Upgrade is deliberately conservative. Project-level defaults only apply
            # to newly generated manifests, never replace existing/manual values here.
            for k in EXPERIMENT_FIELDS:
                if not row.get(k):
                    changes.append(dict(run_accession=rid, field=k, old=row.get(k,''), new='unknown', status='added'))
                    row[k] = 'unknown'
            pending = patch.get(rid, {})
            for k,v in pending.items():
                if k=='run_accession' or not v or v=='unknown' or v==row.get(k):
                    continue
                old = row.get(k,'')
                if old and old != 'unknown':
                    conflicts.append(dict(run_accession=rid, field=k, old=old, proposed=v,
                                          evidence=pending['metadata_evidence'], status='needs_review'))
                else:
                    row[k]=v
                    changes.append(dict(run_accession=rid, field=k, old=old, new=v, status='filled_with_evidence'))
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '.' + uuid.uuid4().hex[:8]
        backup = path.with_name(f'samples.before_v2.{stamp}.tsv')
        shutil.copy2(path, backup)
        write_tsv(path, fields, rows)
        write_tsv(meta/f'metadata_upgrade.{stamp}.changes.tsv', ['run_accession','field','old','new','status'], changes)
        write_tsv(meta/f'metadata_upgrade.{stamp}.conflicts.tsv', ['run_accession','field','old','proposed','evidence','status'], conflicts)
        print(f'Metadata upgraded; backup={backup}; changes={len(changes)}; conflicts={len(conflicts)}')


def report_context(row, db_audit):
    mapping = db_audit['virus_mapping']
    return {**{k: row.get(k,'unknown') or 'unknown' for k in REPORT_METADATA},
            'schema': SCHEMA, 'mapping_schema': mapping['mapping_schema'],
            'mapping_sha256': mapping['mapping_sha256'], 'mapping_build_signature': mapping['build_signature']}


def manifest(args):
    root = args.root.resolve()
    meta = root / 'metadata'
    meta.mkdir(parents=True, exist_ok=True)
    if (meta / 'samples.tsv').exists() and not args.refresh:
        raise RuntimeError('samples.tsv exists; do not overwrite manual curation. Use --refresh only deliberately.')
    records = []
    if args.ena_tsv:
        for p in args.ena_tsv:
            records.extend(parse_ena_report(p.read_text(encoding='utf-8-sig'), str(p)))
        json_write(meta / 'local_metadata.sources.json',
                   [{'path': str(p.resolve()), 'sha256': sha256(p)} for p in args.ena_tsv])
    else:
        if not args.project:
            raise RuntimeError('Provide --project (repeatable) or --ena-tsv.')
        for project in args.project:
            records.extend(fetch_ena(project, meta / f'{project}.ena.tsv'))
    index = {}
    for r in records:
        runid = r.get('run_accession', '')
        if not runid:
            raise RuntimeError('ENA report must contain run_accession')
        if runid in index and index[runid] != r:
            raise RuntimeError(f'Conflicting metadata for {runid}')
        index[runid] = r
    # Deliberately never recurse into classified/QC output folders.
    paths = list(root.glob('*.gz'))
    if (root / 'raw').exists():
        paths.extend((root / 'raw').rglob('*.gz'))
    files_by_run = defaultdict(dict)
    unknown = []
    for p in sorted(paths):
        m = FASTQ_NAME.fullmatch(p.name)
        if m is None:
            if p.name.endswith(('.fastq.gz', '.fq.gz')):
                unknown.append(str(p))
            continue
        runid, mate = m.groups()
        role = {'1': 'r1', '2': 'r2', None: 'single'}[mate]
        if role in files_by_run[runid]:
            if os.path.samefile(p, files_by_run[runid][role]):
                continue
            raise RuntimeError(f'Duplicate role {runid}/{role}: different files; do not infer lanes.')
        files_by_run[runid][role] = p.absolute()
    if unknown:
        raise RuntimeError(f'Unrecognized FASTQ names (manual manifest needed): {unknown[:10]}')
    if not files_by_run:
        raise RuntimeError(f'No ENA-style FASTQ files found in {root} or raw/')
    rows, file_rows, errors, warnings = [], [], [], []
    for runid, files in sorted(files_by_run.items()):
        if runid not in index:
            raise RuntimeError(f'Local run {runid} not in metadata; fetch correct project or complete table.')
        r = index[runid]
        # Required file provenance fields; a summary without FASTQ file columns is insufficient.
        for field in ENA_REQUIRED_FIELDS:
            if field not in r:
                raise RuntimeError(f'ENA table missing {field}; request a full file report.')
        qc, assay, masked, notes = project_policy(r)
        sample = r.get('sample_accession', '')
        unresolved = not sample or r.get('study_accession') == 'PRJEB34634'
        sample_key = f'UNRESOLVED_{runid}' if unresolved else slug(sample)
        mapping = 'unresolved_run_isolated' if unresolved else 'archive_mapping_not_independent_sample_count'
        sources = r['fastq_ftp'].split(';') if r['fastq_ftp'] else []
        md5s = r['fastq_md5'].split(';') if r['fastq_md5'] else [''] * len(sources)
        sizes = r['fastq_bytes'].split(';') if r['fastq_bytes'] else [''] * len(sources)
        if len(md5s) != len(sources) or len(sizes) != len(sources):
            raise RuntimeError(f'Inconsistent ENA file arrays for {runid}')
        expected = {Path(url).name: (url, md5, size) for url, md5, size in zip(sources, md5s, sizes)}
        local_names = {p.name for p in files.values()}
        run_errors = []
        if not expected:
            run_errors.append('No ENA generated FASTQ file list; cannot verify file roles')
        if local_names - set(expected):
            run_errors.append('Local filenames absent from ENA file report: ' + ','.join(local_names - set(expected)))
        if set(expected) - local_names:
            run_errors.append('Missing expected FASTQ file(s): ' + ','.join(set(expected) - local_names))
        if ('r1' in files) != ('r2' in files):
            run_errors.append('Missing mate: _1 and _2 must both exist; never treat missing mate as intentional SE')
        pe = 'r1' in files and 'r2' in files
        layout = 'PE+ARCHIVE_SINGLE' if pe and 'single' in files else ('PE' if pe else 'SINGLE_FILE')
        if not pe and r.get('library_layout', '').upper() == 'PAIRED':
            # Could be merged/interleaved: requires manual review before processing.
            assay = 'review_required'
            notes.append('PAIRED library but only one unpaired-named file: resolve merged vs interleaved before run.')
        out = {field: '' for field in SAMPLE_FIELDS}
        out.update({field: r.get(field, '') for field in SAMPLE_FIELDS if field in r})
        out.update(experiment_defaults(r))
        out.update(dataset=args.dataset or root.name, project=r.get('study_accession', ''),
                   run_accession=runid, sample_key=sample_key, mapping_status=mapping,
                   observed_layout=layout, qc_policy=qc, assay_class=assay,
                   human_read_masking=masked, include='yes', valid='yes')
        for role, p in files.items():
            url, md5, size = expected.get(p.name, ('', '', ''))
            raw = root / 'raw' / sample_key / runid / p.name
            actual_size = p.stat().st_size
            if size and actual_size != int(size):
                run_errors.append(f'{p.name}: size mismatch {actual_size} vs ENA {size}')
            if Path(str(p) + '.aria2').exists():
                run_errors.append(f'{p.name}: .aria2 control file still present')
            if actual_size == 0:
                run_errors.append(f'{p.name}: zero-byte file')
            out[role] = str(raw)
            file_rows.append(dict(run_accession=runid, role=role, source_path=str(p),
                                  raw_path=str(raw), bytes=actual_size, ena_bytes=size,
                                  ena_md5=md5, ena_ftp=url))
        if run_errors:
            out['valid'] = 'no'
            errors.extend([f'{runid}: {x}' for x in run_errors])
            notes.extend(run_errors)
        if assay not in ('shotgun_DNA', 'RNA_seq'):
            warnings.append(f'{runid}: {assay}; organized, but not automatically included in shotgun analysis.')
        out['notes'] = ' | '.join(notes)
        rows.append(out)
    write_tsv(meta / 'samples.tsv', SAMPLE_FIELDS, rows)
    write_tsv(meta / 'files.tsv', FILE_FIELDS, file_rows)
    json_write(meta / 'inventory.audit.json', {'created_utc': now(), 'version': VERSION,
               'local_runs': len(rows), 'local_files': len(file_rows), 'errors': errors, 'warnings': warnings,
               'md5_verified': False, 'biological_samples_confirmed': False})
    print(f'Manifest: {meta / "samples.tsv"}; {len(rows)} runs, {len(file_rows)} files; '
          f'{len(errors)} errors, {len(warnings)} review warnings.')
    if errors:
        raise RuntimeError('Inventory has errors. See inventory.audit.json; no files were moved.')


def source_file(row):
    # Organized file first. Original source retained for provenance and symlink operation.
    for key in ('raw_path', 'source_path'):
        p = Path(row[key])
        if p.is_file():
            if p.stat().st_size != int(row['bytes']):
                raise RuntimeError(f'File size changed: {p}')
            return p
    raise RuntimeError(f'Missing raw and original file: {row}')


def organize(args):
    root, meta = args.root.resolve(), args.root.resolve() / 'metadata'
    rows = selected(read_tsv(meta / 'samples.tsv'), args.runs)
    wanted = {r['run_accession'] for r in rows}
    if any(r['valid'] != 'yes' for r in rows):
        raise RuntimeError('Selected run has invalid inventory; fix manifest first.')
    logpath = meta / 'organize.jsonl'
    receipts = set()
    if logpath.exists():
        for line in logpath.read_text().splitlines():
            r = json.loads(line)
            receipts.add((r['source'], r['destination']))
    with lock(meta / '.organize.lock'), logpath.open('a') as log:
        for row in read_tsv(meta / 'files.tsv'):
            if row['run_accession'] not in wanted:
                continue
            src, dst = Path(row['source_path']), Path(row['raw_path'])
            # Prevent a manually corrupted manifest from writing out of raw/.
            if '..' in dst.parts or not dst.parent.resolve().is_relative_to((root / 'raw').resolve()):
                raise RuntimeError(f'Unsafe raw destination: {dst}')
            if dst.is_file():
                if src.exists() and os.path.samefile(src, dst):
                    continue
                if ((str(src), str(dst)) in receipts and not src.exists()
                        and dst.stat().st_size == int(row['bytes'])):
                    continue
                raise RuntimeError(f'Refuse existing unrelated target: {dst}')
            if dst.is_symlink():
                raise RuntimeError(f'Broken symlink: {dst}')
            if not src.is_file() or src.stat().st_size != int(row['bytes']):
                raise RuntimeError(f'Original file missing/changed: {src}')
            dst.parent.mkdir(parents=True, exist_ok=True)
            if args.mode == 'symlink':
                dst.symlink_to(src.resolve())
            elif args.mode == 'hardlink':
                os.link(src, dst)
            else:
                # Same-filesystem atomic rename; intentionally NO copy-and-delete fallback.
                os.rename(src, dst)
            entry = {'time_utc': now(), 'mode': args.mode, 'source': str(src), 'destination': str(dst)}
            log.write(json.dumps(entry, ensure_ascii=False) + '\n')
            log.flush()
            os.fsync(log.fileno())
    print(f'raw organization completed: {args.mode}. No FASTQ contents were modified.')


def verify(args):
    if not (args.md5 or args.gzip):
        raise RuntimeError('Specify --md5 and/or --gzip. Inventory already checks byte size.')
    root = args.root.resolve()
    wanted = {r['run_accession'] for r in selected(read_tsv(root / 'metadata/samples.tsv'), args.runs)}
    results, bad = [], False
    for row in read_tsv(root / 'metadata/files.tsv'):
        if row['run_accession'] not in wanted:
            continue
        status, error, reads = 'PASS', '', ''
        try:
            p = source_file(row)
            if args.md5:
                if not row['ena_md5']:
                    raise RuntimeError('No expected ENA MD5; cannot certify checksum')
                h = hashlib.md5()
                with p.open('rb') as f:
                    for chunk in iter(lambda: f.read(8 << 20), b''):
                        h.update(chunk)
                if h.hexdigest().lower() != row['ena_md5'].lower():
                    raise RuntimeError('MD5 mismatch')
            if args.gzip:
                reads = sum(1 for _ in fastq(p))
        except Exception as exc:
            status, error, bad = 'FAIL', str(exc), True
        results.append({'run': row['run_accession'], 'role': row['role'],
                        'path': str(source_file(row)) if status == 'PASS' else row['raw_path'],
                        'md5_checked': args.md5, 'gzip_fastq_checked': args.gzip,
                        'reads': reads, 'status': status, 'error': error})
        print(f'{row["run_accession"]} {row["role"]}: {status}', flush=True)
    name = f'verify_{datetime.now().strftime("%Y%m%dT%H%M%S")}.tsv'
    write_tsv(root / 'metadata' / name, list(results[0]) if results else ['status'], results)
    if bad:
        raise RuntimeError('Verification failed; see verification table.')


def tool_version(command):
    p = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode:
        raise RuntimeError(f'Cannot get version: {command}')
    return p.stdout.strip()


def run_command(cmd, logfile, cwd=None):
    print('$ ' + shlex.join(map(str, cmd)), flush=True)
    with logfile.open('w') as log:
        subprocess.run(list(map(str, cmd)), check=True, stdout=log, stderr=log, cwd=cwd)


def prepare_inputs(row, raw_files, work: Path, args):
    """QC preserves surviving mates as singleton streams. Each stream remains distinct."""
    clean, qc = work / 'clean', work / 'qc'
    clean.mkdir()
    qc.mkdir()
    paths = {r['role']: source_file(r) for r in raw_files}
    for role, p in paths.items():
        if not Path(row[role]).is_file():
            raise RuntimeError(f'Run organize first; expected raw link/file: {row[role]}')
    if args.max_records:
        head = work / 'pilot_raw'
        head.mkdir()
        for role, p in list(paths.items()):
            dest = head / f'{role}.fastq.gz'
            with gzip.open(dest, 'wb', compresslevel=1) as out:
                for i, rec in enumerate(fastq(p)):
                    if i == args.max_records:
                        break
                    out.write(b''.join(rec))
            paths[role] = dest
    jobs = []
    if 'r1' in paths:
        jobs.append(('PE', paths['r1'], paths['r2']))
    if 'single' in paths:
        jobs.append(('ARCHIVE_SINGLE', paths['single'], None))
    streams, raw_reads, failed_reads, all_passthrough = [], 0, 0, True
    qc_jobs = {}
    policy = args.qc if args.qc != 'auto' else row['qc_policy']
    if policy not in ('fastp', 'minimal', 'passthrough'):
        raise RuntimeError(f'Unknown QC policy {policy}; edit samples.tsv deliberately.')
    for label, r1, r2 in jobs:
        if policy == 'passthrough':
            # Includes input reads unchanged; count through the checked classifier splitter.
            streams.append((label, r1, r2))
            continue
        all_passthrough = False
        out1 = clean / f'{label}.R1.fastq.gz'
        out2 = clean / f'{label}.R2.fastq.gz'
        failed = qc / f'{label}.failed.fastq.gz'
        report = qc / f'{label}.fastp.json'
        cmd = ['fastp', '--in1', str(r1), '--out1', str(out1), '--thread', str(min(args.threads, 16)),
               '--json', str(report), '--html', str(qc / f'{label}.fastp.html'),
               '--failed_out', str(failed), '--qualified_quality_phred', str(args.q),
               '--unqualified_percent_limit', '40', '--n_base_limit', '5',
               '--length_required', str(args.min_length), '--compression', '1',
               '--disable_trim_poly_g', '--dont_eval_duplication', '--dont_overwrite']
        if policy == 'minimal':
            cmd.append('--disable_adapter_trimming')
        if r2:
            u1, u2 = clean / 'PE.U1.fastq.gz', clean / 'PE.U2.fastq.gz'
            cmd += ['--in2', str(r2), '--out2', str(out2), '--unpaired1', str(u1), '--unpaired2', str(u2)]
            if policy == 'fastp':
                cmd.append('--detect_adapter_for_pe')
            streams.extend([(label, out1, out2), ('QC_U1', u1, None), ('QC_U2', u2, None)])
        else:
            streams.append((label, out1, None))
        run_command(cmd, qc / f'{label}.fastp.log', cwd=work)
        data = json.loads(report.read_text())
        before = int(data['summary']['before_filtering']['total_reads'])
        after = int(data['summary']['after_filtering']['total_reads'])
        filtering = data['filtering_result']
        filtered = sum(int(v) for k,v in filtering.items() if k.endswith('_reads') and k != 'passed_filter_reads')
        if min(before, after, filtered) < 0 or int(filtering['passed_filter_reads']) != after or before != after + filtered:
            raise RuntimeError('Inconsistent fastp summary/filtering counts')
        emitted_failed = sum(1 for _ in fastq(failed)) if failed.exists() else 0
        raw_reads += before
        failed_reads += emitted_failed
        qc_jobs[label] = dict(raw_reads=before, paired=bool(r2), passed_primary_reads=after,
                              filtered_report_reads=filtered, failed_fastq_reads=emitted_failed)
    json_write(qc/'jobs.json', qc_jobs)
    return streams, raw_reads, failed_reads, all_passthrough, policy, qc_jobs


def qc_failure_count(jobs, stream_reads):
    """PE fastp filtering counts both mates of a failing pair, including survivors.

    failed_out is diagnostic only: fastp 1.4.0 omits double-failed pairs there.
    Validate the independently reported primary outputs, then subtract rescued
    singleton reads from the reported failed-pair read count. Raw stays intact.
    """
    failures = 0
    for label, job in jobs.items():
        if stream_reads[label] != job['passed_primary_reads']:
            raise RuntimeError(f'fastp report/actual primary output mismatch for {label}')
        rescued = stream_reads['QC_U1'] + stream_reads['QC_U2'] if job['paired'] else 0
        count = job['filtered_report_reads'] - rescued
        if count < 0 or job['failed_fastq_reads'] > count:
            raise RuntimeError('fastp failed/singleton counts inconsistent')
        failures += count
    return failures


def signature(row, files, db_audit, args, versions):
    scripts = {p.name: sha256(p) for p in Path(__file__).parent.glob('*.py')}
    return hashlib.sha256(json.dumps({'row': row, 'sources': [
        {'path': str(source_file(r).resolve()), 'size': source_file(r).stat().st_size,
         'mtime_ns': source_file(r).stat().st_mtime_ns} for r in files],
        'db': db_audit, 'args': {k: getattr(args, k) for k in (
            'confidence', 'hit_groups', 'q', 'min_length', 'qc', 'max_records',
            'memory_mapping', 'keep_clean', 'keep_kraken_output', 'allow_incomplete_db',
            'allow_incomplete_virus', 'allow_reviewed_targeted')},
        'version': VERSION, 'schema': SCHEMA, 'scripts': scripts, 'tools': versions}, sort_keys=True).encode()).hexdigest()


def run_one(row, files, output: Path, tax, db_audit, versions, args):
    runid, sample = row['run_accession'], row['sample_key']
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', sample) or sample in ('.', '..'):
        raise RuntimeError(f'Unsafe sample_key {sample}')
    if not files or len({r['role'] for r in files}) != len(files):
        raise RuntimeError(f'Missing/duplicate file manifest records for {runid}')
    if not re.fullmatch(r'[SED]RR[0-9]+', runid):
        raise RuntimeError(f'Unsafe run_accession {runid}')
    tax.virus.assert_unchanged()
    context = report_context(row, db_audit)
    report_dir = output / 'reports' / runid
    done = report_dir / 'done.json'
    sig = signature(row, files, db_audit, args, versions)
    if done.exists():
        status = json.loads(done.read_text())
        if status.get('schema') != SCHEMA or status.get('virus_mapping') != tax.virus.identity or status.get('signature') != sig:
            raise RuntimeError(f'{runid}: completed output has different settings/input/database. Use a NEW --out directory.')
        for entry in status['outputs']:
            p = output / entry['path']
            if not p.resolve().is_relative_to(output) or not p.is_file() or p.stat().st_size != entry['bytes']:
                raise RuntimeError(f'Completed output missing/size-changed/unsafe: {p}')
            if entry.get('sha256') and sha256(p) != entry['sha256']:
                raise RuntimeError(f'Completed report checksum changed: {p}')
        print(f'{runid}: already completed and size-checked; skip.', flush=True)
        return
    if report_dir.exists():
        raise RuntimeError(f'Incomplete final report directory {report_dir}; inspect it before retrying.')
    for group in (*GROUPS, 'virus'):
        if (output / group / sample / runid).exists():
            raise RuntimeError(f'Pre-existing output without compatible done.json: {group}/{sample}/{runid}')
    work = output / 'work' / f'{runid}.{uuid.uuid4().hex[:10]}.partial'
    work.mkdir(parents=True)
    json_write(work / 'input.json', {'schema': SCHEMA, 'sample_row': row, 'signature': sig, 'versions': versions,
               'virus_mapping': tax.virus.identity, 'db_audit': db_audit, 'started_utc': now()})
    bins = GzipBins(work / 'bins', memory_mb=args.split_buffer_mb)
    counter, totals, virus_counter = Counter(), Counter(), Counter()
    stream_reads = Counter()
    try:
        streams, raw_reads, failed_emitted, passthrough, policy, qc_jobs = prepare_inputs(row, files, work, args)
        for stream, r1, r2 in streams:
            # Empty gzip files are valid when no survivors; check both mates.
            first = next(fastq(r1), None)
            second = next(fastq(r2), None) if r2 else None
            if r2 and ((first is None) != (second is None)):
                raise RuntimeError('One clean mate file is empty while the other is not')
            if first is None:
                continue
            if r2:
                pair_ids(first[0], second[0])
            kr = work / 'kraken'
            kr.mkdir(exist_ok=True)
            cmd = ['kraken2', '--db', str(args.db.resolve()), '--threads', str(args.threads),
                   '--confidence', str(args.confidence), '--minimum-hit-groups', str(args.hit_groups),
                   '--report', str(kr / f'{stream}.kreport'), '--report-minimizer-data', '--gzip-compressed']
            if args.memory_mapping:
                cmd.append('--memory-mapping')
            if r2:
                cmd += ['--paired', str(r1), str(r2)]
            else:
                cmd.append(str(r1))
            print('$ ' + shlex.join(cmd), flush=True)
            (kr / f'{stream}.command.txt').write_text(shlex.join(cmd) + '\n')
            with ExitStack() as stack:
                err = stack.enter_context((kr / f'{stream}.log').open('w'))
                assignments = stack.enter_context(gzip.open(kr / f'{stream}.assignments.tsv.gz', 'wb', compresslevel=1))
                assignments.write(('\t'.join(ASSIGNMENT_FIELDS) + '\n').encode())
                full = stack.enter_context(gzip.open(kr / f'{stream}.kraken.gz', 'wb', compresslevel=1)) if args.keep_kraken_output else None
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, cwd=work)
                try:
                    c, t = split_stream(proc.stdout, r1, r2, stream, tax, bins, assignments, full, virus_counter)
                    rc = proc.wait()
                    if rc:
                        raise RuntimeError(f'Kraken2 failed (exit {rc}); see {kr / (stream + ".log")}')
                except BaseException:
                    proc.terminate()
                    try:
                        proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
                    raise
                finally:
                    proc.stdout.close()
                counter.update(c)
                totals.update(t)
                stream_reads[stream] += t['reads']
        bins.close()
        if passthrough:
            raw_reads = totals['reads']
            failed_reads = 0
        else:
            failed_reads = qc_failure_count(qc_jobs, stream_reads)
        if raw_reads != totals['reads'] + failed_reads:
            raise RuntimeError(f'Read conservation failed: raw={raw_reads}, clean={totals["reads"]}, '
                               f'failed={failed_reads}. Inspect fastp version/orphan behavior; do not publish.')
        counts = []
        for (group, bucket, species, name, stream, mates), fragments in sorted(counter.items()):
            counts.append(dict(dataset=row['dataset'], sample_key=sample, run_accession=runid,
                               group=group, bucket=bucket, species_taxid=species, species_name=name,
                               stream=stream, fragments=fragments, reads=fragments * mates,
                               scope='pilot_subset' if args.max_records else 'full_run', **context))
        write_tsv(work / 'counts.tsv', COUNTS_FIELDS, counts)
        group_rows = []
        for group in GROUPS:
            group_rows.append({'group': group, 'fragments': sum(x['fragments'] for x in counts if x['group'] == group),
                               'reads': sum(x['reads'] for x in counts if x['group'] == group)})
        write_tsv(work / 'groups.tsv', ['group', 'fragments', 'reads'], group_rows)
        virus_reads = sum(r['reads'] for r in group_rows if r['group'].startswith('virus_'))
        virus_fragments = sum(r['fragments'] for r in group_rows if r['group'].startswith('virus_'))
        taxa_rows, annotation_rows = [], []
        seen_taxids = set()
        for (tid, stream, mates), fragments in sorted(virus_counter.items()):
            ann = tax.annotation(tid)
            group, _, species, species_name = tax.route(tid)
            base = dict(dataset=row['dataset'], sample_key=sample, run_accession=runid,
                        original_assigned_taxid=tid, canonical_taxid=ann['canonical_taxid'],
                        assigned_rank=tax.node(tid)[1], assigned_name=tax.node(tid)[2])
            taxa_rows.append(dict(base, **context, scope='pilot_subset' if args.max_records else 'full_run',
                stream=stream, species_taxid=species, species_name=species_name, group=group,
                genome_class=ann['genome_class'], genome_type_detail=ann['genome_type_detail'],
                resolution_status=ann['resolution_status'], evidence_taxid=ann['evidence_taxid'],
                fragments=fragments, reads=fragments*mates))
            if tid not in seen_taxids:
                annotation_rows.append({**base, **context, **ann})
                seen_taxids.add(tid)
        if sum(r['reads'] for r in taxa_rows) != virus_reads or sum(r['fragments'] for r in taxa_rows) != virus_fragments:
            raise RuntimeError('Virus taxon/group conservation failed')
        write_tsv(work/'virus_taxa.tsv', VIRUS_TAXA_FIELDS, taxa_rows)
        write_tsv(work/'virus_annotation.tsv', AUDIT_FIELDS, annotation_rows)
        stats = {'raw_reads_processed': raw_reads, 'clean_reads': totals['reads'],
                 'clean_fragments': totals['fragments'], 'clean_bases': totals['bases'],
                 'qc_failed_reads': failed_reads, 'qc_failed_reads_in_failed_fastq': failed_emitted,
                 'qc_failed_reads_not_emitted': failed_reads - failed_emitted,
                 'qc_count_method': 'fastp filtering report minus rescued singleton reads; primary output independently checked',
                 'read_conservation_passed': True,
                 'virus_all_reads': virus_reads, 'virus_all_fragments': virus_fragments,
                 'virus_conservation_passed': True,
                 'missing_database_groups': ','.join(db_audit['missing_groups']),
                 'missing_virus_classes': ','.join(db_audit['missing_virus_classes']),
                 **context,
                 'scope': 'pilot_subset' if args.max_records else 'full_run', 'qc_policy': policy,
                 'warning': 'Taxonomic candidate assignments, NOT validated species presence or biomass abundance.'}
        json_write(work / 'run.stats.json', stats)
        tax.virus.assert_unchanged()
        if snapshot(args.db) != db_audit['snapshot']:
            raise RuntimeError('Database changed during run')
        # No final outputs become visible until ALL streams have passed consistency checks.
        outputs = []
        for group in GROUPS:
            src = work / 'bins' / group
            if not src.exists():
                continue
            dst = output / group / sample / runid
            if dst.exists():
                raise RuntimeError(f'Refuse pre-existing species output: {dst}')
        if not args.keep_clean:
            shutil.rmtree(work / 'clean')
            shutil.rmtree(work / 'pilot_raw', ignore_errors=True)
        for group in GROUPS:
            src = work / 'bins' / group
            if not src.exists():
                continue
            dst = output / group / sample / runid
            dst.parent.mkdir(parents=True, exist_ok=True)
            src.rename(dst)
            for p in sorted(dst.rglob('*.gz')):
                outputs.append({'path': str(p.relative_to(output)), 'bytes': p.stat().st_size})
        if (work / 'bins').exists():
            (work / 'bins').rmdir()
        report_dir.parent.mkdir(parents=True, exist_ok=True)
        work.rename(report_dir)
        for p in report_dir.rglob('*'):
            if p.is_file():
                outputs.append({'path': str(p.relative_to(output)), 'bytes': p.stat().st_size,
                                **({'sha256':sha256(p)} if p.suffix in ('.tsv','.json') else {})})
        json_write(done, {'schema': SCHEMA, 'virus_mapping': tax.virus.identity, 'signature': sig, 'completed_utc': now(), 'stats': stats, 'outputs': outputs})
        print(f'{runid}: {raw_reads:,} raw reads -> {totals["reads"]:,} clean + {failed_reads:,} QC-failed; PASS.', flush=True)
    except BaseException as error:
        bins.close()
        if work.exists():
            json_write(work / 'FAILED.json', {'error': str(error), 'time_utc': now()})
        raise


def run(args):
    root, output = args.root.resolve(), (args.out or args.root / 'results_v2').resolve()
    if args.max_records and (args.out is None or output == root):
        raise RuntimeError('Pilot --max-records requires a separate --out, never mix subsets and full outputs.')
    if (not (0 <= args.confidence <= 1) or args.threads < 1 or args.min_length < 35
            or args.split_buffer_mb < 1 or args.max_records < 0 or args.hit_groups < 1 or args.limit < 0):
        raise RuntimeError('Require confidence in [0,1], positive threads, minimum length >=35 for default Kraken2 k.')
    rows = selected(read_tsv(root / 'metadata/samples.tsv'), args.runs, args.limit)
    if len({r['run_accession'] for r in rows}) != len(rows):
        raise RuntimeError('Duplicate run_accession in samples.tsv')
    eligible, skipped = [], []
    for row in rows:
        if row['include'].lower() not in ('yes', 'true', '1'):
            skipped.append({'run': row['run_accession'], 'reason': 'include=no'})
        elif row['valid'] != 'yes':
            raise RuntimeError(f'Invalid inventory for {row["run_accession"]}')
        elif row['assay_class'] not in ('shotgun_DNA', 'RNA_seq') and not args.allow_reviewed_targeted:
            skipped.append({'run': row['run_accession'], 'reason': row['assay_class']})
        else:
            eligible.append(row)
    print(f'Eligible: {len(eligible)}; skipped: {len(skipped)}')
    for r in skipped[:20]:
        print(f'SKIP {r["run"]}: {r["reason"]}')
    if args.dry_run:
        for r in eligible:
            print(r['run_accession'], r['sample_key'], r['observed_layout'], r['qc_policy'])
        return
    if not eligible:
        raise RuntimeError('No eligible runs. Review metadata; do not relabel targeted libraries as shotgun.')
    for program in ('kraken2', 'fastp'):
        if not shutil.which(program):
            raise RuntimeError(f'Missing executable {program}; activate metag environment.')
    versions = {'kraken2': tool_version(['kraken2', '--version']), 'fastp': tool_version(['fastp', '--version'])}
    if args.threads > 1 and not re.search(r'(?<![0-9])2\.17\.1(?![0-9.])', versions['kraken2']):
        raise RuntimeError('Multi-thread output ordering verified for Kraken2 2.17.1 only; use that version or --threads 1')
    versions['classifier_order'] = '2.17.1 ordered block queue' if args.threads > 1 else 'single_thread'
    args.db = args.db.resolve()
    audit_path = args.db / 'metag.audit.json'
    if not audit_path.exists():
        raise RuntimeError('Run db-audit first. Group names alone do not prove reference coverage.')
    db_audit = json.loads(audit_path.read_text())
    if db_audit.get('schema') != SCHEMA or not db_audit.get('virus_mapping'):
        raise RuntimeError('Require v2 db-audit --virus-map before processing')
    if db_audit['snapshot'] != snapshot(args.db):
        raise RuntimeError('Database changed since audit; run db-audit again, use new output if already analyzed.')
    if db_audit['missing_groups'] and not args.allow_incomplete_db:
        raise RuntimeError(f'Database missing reference support for {db_audit["missing_groups"]}. STOP.')
    if db_audit['missing_virus_classes'] and not args.allow_incomplete_virus:
        raise RuntimeError(f'Insufficient DNA/RNA virus index support: {db_audit["missing_virus_classes"]}')
    # Database loading is the dominant memory demand, not the read splitter.
    if not args.memory_mapping and Path('/proc/meminfo').exists():
        mem = dict((line.split(':')[0], int(line.split()[1]) * 1024)
                   for line in Path('/proc/meminfo').read_text().splitlines() if ':' in line)
        available = mem.get('MemAvailable', 0)
        if available and (args.db / 'hash.k2d').stat().st_size > 0.85 * available:
            raise RuntimeError('Database hash exceeds 85% of currently available RAM. Free RAM, use a smaller/fuller-memory setup, or explicitly use --memory-mapping (may be I/O-bound).')
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'virus').exists():
        raise RuntimeError('Legacy virus/ output found; use a NEW --out directory')
    raw = read_tsv(root / 'metadata/files.tsv')
    by_run = defaultdict(list)
    for r in raw:
        by_run[r['run_accession']].append(r)
    with lock(output / '.metag.run.lock'):
        tax = Taxonomy(prepare_sqlite(args.db), args.virus_map)
        if db_audit['virus_mapping'] != tax.virus.identity:
            tax.close()
            raise RuntimeError('Virus mapping differs from db-audit; re-audit first')
        for completed in (output/'reports').glob('*/done.json'):
            state = json.loads(completed.read_text())
            if (state.get('schema') != SCHEMA or state.get('virus_mapping') != tax.virus.identity
                    or state.get('stats',{}).get('scope') != ('pilot_subset' if args.max_records else 'full_run')):
                tax.close()
                raise RuntimeError('Output root contains incompatible schema/mapping/scope; use a NEW --out')
        for group in GROUPS:
            (output / group).mkdir(exist_ok=True)
        write_tsv(output / 'skipped_runs.tsv', ['run', 'reason'], skipped)
        try:
            for row in eligible:
                run_one(row, by_run[row['run_accession']], output, tax, db_audit, versions, args)
        finally:
            tax.close()
    summarize(argparse.Namespace(out=output))


def summarize(args):
    output = args.out.resolve()
    long_rows, group_rows, stats_rows, taxa_rows, species_rows, genome_rows, audit_rows = [], [], [], [], [], [], []
    mapping_identity = None
    summary_scope = None
    for done in sorted((output / 'reports').glob('*/done.json')):
        directory = done.parent
        status = json.loads(done.read_text())
        inp = json.loads((directory / 'input.json').read_text())
        if status.get('schema') != SCHEMA or inp.get('schema') != SCHEMA:
            raise RuntimeError(f'Incompatible old output schema at {directory}; use separate v2 outputs')
        if not status.get('signature') or status['signature'] != inp.get('signature'):
            raise RuntimeError(f'Completion/input signature mismatch: {directory}')
        mapping = status.get('virus_mapping')
        if not mapping or mapping != inp.get('virus_mapping') or mapping != inp.get('db_audit',{}).get('virus_mapping'):
            raise RuntimeError(f'Completion/input mapping mismatch: {directory}')
        if mapping_identity is not None and mapping != mapping_identity:
            raise RuntimeError('Cannot combine different frozen virus/taxonomy mappings in a summary')
        mapping_identity = mapping
        for entry in status['outputs']:
            file = output / entry['path']
            if not file.resolve().is_relative_to(output) or not file.is_file() or file.stat().st_size != entry['bytes']:
                raise RuntimeError(f'Missing/changed/unsafe completed output: {file}')
            if entry.get('sha256') and sha256(file) != entry['sha256']:
                raise RuntimeError(f'Completed report checksum mismatch: {file}')
        sample = inp['sample_row']
        ctx = report_context(sample, inp['db_audit'])
        common = dict(dataset=sample['dataset'], sample_key=sample['sample_key'], run_accession=sample['run_accession'])
        counts = read_tsv(directory/'counts.tsv')
        groups = read_tsv(directory/'groups.tsv')
        taxa = read_tsv(directory/'virus_taxa.tsv')
        annotations = read_tsv(directory/'virus_annotation.tsv')
        stat = json.loads((directory/'run.stats.json').read_text())
        for r in counts + taxa + annotations:
            if any(r.get(k) != str(ctx[k]) for k in PROVENANCE_FIELDS):
                raise RuntimeError(f'Row provenance mismatch in {directory}')
        if summary_scope is not None and stat.get('scope') != summary_scope:
            raise RuntimeError('Cannot combine pilot subsets and full runs in one summary')
        summary_scope = stat.get('scope')
        if stat != status.get('stats'):
            raise RuntimeError(f'Completion/statistics mismatch: {directory}')
        if {r['group'] for r in groups} != set(GROUPS) or len(groups) != len(GROUPS):
            raise RuntimeError('Missing/duplicate/old output group')
        for r in groups:
            if any(int(r[k]) != sum(int(c[k]) for c in counts if c['group']==r['group']) for k in ('reads','fragments')):
                raise RuntimeError('Species/group conservation failed in summary')
        if (stat['raw_reads_processed'] != stat['clean_reads'] + stat['qc_failed_reads'] or
            stat['clean_reads'] != sum(int(r['reads']) for r in groups) or
            stat['clean_fragments'] != sum(int(r['fragments']) for r in groups)):
            raise RuntimeError('Read/fragment conservation failed in summary')
        viruses = [r for r in groups if r['group'].startswith('virus_')]
        for k in ('reads','fragments'):
            total = sum(int(r[k]) for r in viruses)
            if total != sum(int(r[k]) for r in taxa) or total != stat['virus_all_'+k]:
                raise RuntimeError('Virus taxon/group conservation failed in summary')
        long_rows.extend(counts)
        taxa_rows.extend(taxa)
        audit_rows.extend(annotations)
        stats_rows.append({**common, **stat})
        for r in groups:
            group_rows.append({**common, **ctx, 'scope':stat['scope'], **r})
        # Original counts already aggregate all assigned nodes within a species bin.
        # Do not sum this view with the taxon-level view.
        for r in counts:
            if r['group'].startswith('virus_') and int(r['species_taxid'])>0:
                species_rows.append(dict(r, genome_class=r['group'][6:]))
        all_nodes = {r['canonical_taxid'] for r in taxa}
        unknown_nodes = {r['canonical_taxid'] for r in taxa if r['genome_class']=='unknown'}
        virus_reads = stat['virus_all_reads']
        unknown_reads = sum(int(r['reads']) for r in viruses if r['group']=='virus_unknown')
        for r in viruses:
            genome_rows.append({**common, **ctx, 'scope':stat['scope'], **r,
                'genome_class':r['group'][6:], 'virus_all_reads':virus_reads,
                'virus_all_fragments':stat['virus_all_fragments'], 'virus_nodes':len(all_nodes),
                'unknown_nodes':len(unknown_nodes),
                'unknown_read_fraction':unknown_reads/virus_reads if virus_reads else 'NA',
                'unknown_node_fraction':len(unknown_nodes)/len(all_nodes) if all_nodes else 'NA',
                'reads_per_million_clean':int(r['reads'])*1000000/stat['clean_reads'] if stat['clean_reads'] else 'NA'})
    common_fields = ['dataset','sample_key','run_accession','scope'] + REPORT_METADATA + PROVENANCE_FIELDS
    write_tsv(output/'summary/species_by_run.tsv', COUNTS_FIELDS, long_rows)
    write_tsv(output/'summary/groups_by_run.tsv', common_fields+['group','fragments','reads'], group_rows)
    write_tsv(output/'summary/qc_by_run.tsv', list(stats_rows[0]) if stats_rows else common_fields+['raw_reads_processed','clean_reads','qc_failed_reads'], stats_rows)
    write_tsv(output/'summary/virus_taxa_by_run.tsv', VIRUS_TAXA_FIELDS, taxa_rows)
    write_tsv(output/'summary/virus_species_by_run.tsv', COUNTS_FIELDS+['genome_class'], species_rows)
    write_tsv(output/'summary/virus_genome_by_run.tsv', common_fields+['group','genome_class','fragments','reads',
              'virus_all_reads','virus_all_fragments','virus_nodes','unknown_nodes','unknown_read_fraction',
              'unknown_node_fraction','reads_per_million_clean'], genome_rows)
    write_tsv(output/'summary/virus_annotation_audit.tsv', AUDIT_FIELDS, audit_rows)
    json_write(output/'summary/summary.source.json', dict(schema=SCHEMA, virus_mapping=mapping_identity,
                completed_runs=len(stats_rows), metadata_source='per-run input.json snapshot'))
    print(f'Summaries written to {output / "summary"}; NO biological-sample pooling performed.')


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version=VERSION)
    sub = p.add_subparsers(dest='command', required=True)
    a = sub.add_parser('manifest', help='Read ENA metadata and inventory LOCAL files only')
    a.add_argument('--root', type=Path, required=True)
    a.add_argument('--project', action='append')
    a.add_argument('--ena-tsv', type=Path, action='append')
    a.add_argument('--dataset', default='')
    a.add_argument('--refresh', action='store_true')
    a.set_defaults(func=manifest)
    a = sub.add_parser('metadata-upgrade', help='Back up and append metadata; conflicts remain pending review')
    a.add_argument('--root', type=Path, required=True)
    a.add_argument('--updates', type=Path, help='Optional evidence-bearing TSV keyed by run_accession')
    a.set_defaults(func=metadata_upgrade)
    a = sub.add_parser('organize', help='Create raw/sample/run files; default is nondestructive symlinks')
    a.add_argument('--root', type=Path, required=True)
    a.add_argument('--mode', choices=['symlink', 'hardlink', 'move'], default='symlink')
    a.add_argument('--runs', default='')
    a.set_defaults(func=organize)
    a = sub.add_parser('verify', help='Full ENA MD5 and/or gzip/FASTQ validation')
    a.add_argument('--root', type=Path, required=True)
    a.add_argument('--runs', default='')
    a.add_argument('--md5', action='store_true')
    a.add_argument('--gzip', action='store_true')
    a.set_defaults(func=verify)
    a = sub.add_parser('db-audit', help='Audit same-build taxonomy and actual index coverage')
    a.add_argument('--db', type=Path, required=True)
    a.add_argument('--threads', type=int, default=8)
    a.add_argument('--inspect-file', type=Path, help='Existing SIX-column k2 inspect report from this exact DB')
    choice = a.add_mutually_exclusive_group(required=True)
    choice.add_argument('--virus-map', type=Path)
    choice.add_argument('--base-only', action='store_true', help='Seven-group audit only; NOT ready for v2 run')
    a.add_argument('--allow-incomplete-virus', action='store_true')
    a.set_defaults(func=lambda a: print(json.dumps(audit_database(a.db.resolve(), a.threads, a.inspect_file,
                        a.virus_map, a.allow_incomplete_virus), indent=2)))
    a = sub.add_parser('run', help='Sequential runs; threaded classifier; all surviving reads retained')
    a.add_argument('--root', type=Path, required=True)
    a.add_argument('--db', type=Path, required=True)
    a.add_argument('--out', type=Path, help='Default: ROOT/results_v2')
    a.add_argument('--virus-map', type=Path, required=True)
    a.add_argument('--allow-incomplete-virus', action='store_true', help='Exploratory only; missing types not adequately tested')
    a.add_argument('--runs', default='', help='Comma-separated run IDs')
    a.add_argument('--limit', type=int, default=0)
    a.add_argument('--threads', type=int, default=16)
    a.add_argument('--confidence', type=float, default=0.1)
    a.add_argument('--hit-groups', type=int, default=3)
    a.add_argument('--min-length', type=int, default=50)
    a.add_argument('--q', type=int, default=20)
    a.add_argument('--qc', choices=['auto', 'fastp', 'minimal', 'passthrough'], default='auto')
    a.add_argument('--max-records', type=int, default=0, help='Pilot: limit each raw stream; requires separate --out')
    a.add_argument('--memory-mapping', action='store_true')
    a.add_argument('--split-buffer-mb', type=int, default=64)
    a.add_argument('--keep-clean', action='store_true')
    a.add_argument('--keep-kraken-output', action='store_true')
    a.add_argument('--allow-incomplete-db', action='store_true', help='EXPLICIT exploratory override; missing groups not tested')
    a.add_argument('--allow-reviewed-targeted', action='store_true',
                   help='Classify manually reviewed targeted/control runs only; not ecological abundance inference')
    a.add_argument('--dry-run', action='store_true')
    a.set_defaults(func=run)
    a = sub.add_parser('summary', help='Summarize only completed runs, without pooling')
    a.add_argument('--out', type=Path, required=True)
    a.set_defaults(func=summarize)
    return p


def main():
    args = parser().parse_args()
    try:
        args.func(args)
    except (RuntimeError, OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
