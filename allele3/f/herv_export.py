#!/usr/bin/env python3
"""Publish complete HERV censuses as portable chromosome tables (no reanalysis)."""
from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import hashlib
import json
from pathlib import Path
import shutil

SCHEMA = 'herv-chromosome-export-1'
FIELDS = ['chrom', 'start0', 'end0', 'sample', 'haplotype', 'family', 'pop', 'super_pop']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_raw_sample(folder, chrom):
    from herv_cohort import annotation_intervals
    return {hap: annotation_intervals(folder, chrom, hap) for hap in ('hap1', 'hap2')}


def export_cohort(out, io_workers=4):
    from herv_cohort import (CHROMS, FAMILIES, annotation_intervals, complete_sample,
                             file_lock, read_json, save_json, sample_directory,
                             projection_directory, comparable_sequence_settings)
    out = Path(out).resolve()
    with file_lock(out/'tmp/export.lock'):
        plan = read_json(out/'qc/plan.json')
        summary = read_json(out/'qc/summary.json')
        samples = plan['samples']
        if plan['identity']['build'] != 'GRCh38':
            raise ValueError('Only verified native GRCh38 censuses can be exported')
        if summary['status'] != 'complete' or summary['missing_samples'] or summary['completed_samples'] != len(samples):
            raise ValueError('Export requires a complete, summarized frozen cohort')
        if len({r['sample'] for r in samples}) != len(samples):
            raise ValueError('Duplicate samples in frozen cohort')
        # Metadata waits dominate on the shared filesystem. Bounded I/O
        # concurrency preserves all checks without using extra compute workers.
        with ThreadPoolExecutor(max_workers=min(8, len(samples))) as pool:
            pending = {pool.submit(complete_sample, sample_directory(out, r['sample']), r['sample']): r['sample'] for r in samples}
            for index, future in enumerate(as_completed(pending), 1):
                sample = pending[future]
                if not future.result():
                    raise ValueError('Incomplete census: ' + sample)
                print(f'[HERV export] validated {index}/{len(samples)} {sample}', flush=True)
        settings = comparable_sequence_settings(out, samples)
        if settings != summary['sequence_annotation_settings']:
            raise ValueError('Sample settings differ from completed summary')
        for kind in ('tmp/export_tables/raw', 'tmp/export_tables/segments', 'tmp/export_tables/uncallable', 'qc/export', 'tmp/export'):
            (out/kind).mkdir(parents=True, exist_ok=True)
        # This manifest is the commit marker; stale products must not appear complete.
        manifest_path = out/'qc/export_manifest.json'
        save_json(manifest_path, dict(schema=SCHEMA, status='running', build='GRCh38'))
        receipts = []
        for chrom in CHROMS:
            sources = []
            for row in samples:
                folder = sample_directory(out, row['sample'])
                proj = projection_directory(out, folder)
                for name in ('run.json', 'analysis.done.json', 'plan.json'):
                    sources.append(sha(folder/name))
                for hap in ('hap1', 'hap2'):
                    for suffix in ('.bed.gz', '.uncallable.bed.gz', '.json'):
                        sources.append(sha(proj/(chrom+'.'+hap+suffix)))
            identity = dict(schema=SCHEMA, samples=samples, sources=sources, settings=settings,
                            exporter_sha256=sha(__file__))
            paths = {k: out/'tmp/export_tables'/k/(chrom+'.tsv.gz') for k in ('raw', 'segments', 'uncallable')}
            receipt_path = out/'qc/export'/(chrom+'.json')
            if receipt_path.is_file():
                old = read_json(receipt_path)
                if old.get('identity') == identity and all(p.is_file() and sha(p) == old['files'][k]['sha256'] for k, p in paths.items()):
                    receipts.append(old)
                    print('[HERV export] REUSE '+chrom, flush=True)
                    continue
            handles = {}
            writers = {}
            counts = {k: 0 for k in paths}
            temps = {k: out/'tmp/export'/(chrom+'.'+k+'.partial.gz') for k in paths}
            # Prefetch chunk files while the main thread serializes deterministic rows.
            pool = ThreadPoolExecutor(max_workers=max(1,min(4,io_workers)))
            raw_futures = {r['sample']: pool.submit(read_raw_sample, sample_directory(out, r['sample']), chrom) for r in samples}
            try:
                for k in paths:
                    handles[k] = gzip.open(temps[k], 'wt', compresslevel=3)
                    writers[k] = csv.writer(handles[k], delimiter='\t', lineterminator='\n')
                    writers[k].writerow(FIELDS)
                for row in samples:
                    folder = sample_directory(out, row['sample'])
                    proj = projection_directory(out, folder)
                    for hap in ('hap1', 'hap2'):
                        meta = read_json(proj/(chrom+'.'+hap+'.json'))
                        expected = {r['family']: r for r in meta['counts']}
                        raw = raw_futures[row['sample']].result()[hap]
                        for family in FAMILIES:
                            if sum(b-a for a,b in raw[family]) != expected[family]['raw_haplotype_annotation_bp']:
                                raise ValueError(f'Raw annotation changed since projection: {row["sample"]}/{chrom}/{hap}')
                            for a,b in raw[family]:
                                writers['raw'].writerow([chrom,a,b,row['sample'],hap,family,row['pop'],row['super_pop']])
                                counts['raw'] += 1
                        observed = dict.fromkeys(FAMILIES, 0)
                        for kind, suffix in (('segments','.bed.gz'), ('uncallable','.uncallable.bed.gz')):
                            mask_bp = 0
                            with gzip.open(proj/(chrom+'.'+hap+suffix), 'rt') as stream:
                                for line in stream:
                                    fields = line.rstrip().split('\t')
                                    c,a,b = fields[:3]; a,b = int(a),int(b)
                                    family = fields[3] if kind == 'segments' else 'uncallable'
                                    if c != chrom or not 0 <= a < b <= meta['reference_length']:
                                        raise ValueError('Invalid GRCh38 projection interval')
                                    if kind == 'segments': observed[family] += b-a
                                    else: mask_bp += b-a
                                    writers[kind].writerow([chrom,a,b,row['sample'],hap,family,row['pop'],row['super_pop']])
                                    counts[kind] += 1
                            if kind == 'uncallable' and mask_bp != meta['uncallable_bp']:
                                raise ValueError('Uncallable mask differs from projection receipt')
                        if any(observed[f] != expected[f]['callable_annotation_bp'] for f in FAMILIES):
                            raise ValueError('Projected annotation differs from receipt')
                    del raw_futures[row['sample']]
                for handle in handles.values(): handle.close()
                for k,p in paths.items(): temps[k].replace(p)
                receipt = dict(identity=identity, chrom=chrom, files={k: dict(path=str(p.relative_to(out)),
                    rows=counts[k], sha256=sha(p), bytes=p.stat().st_size) for k,p in paths.items()})
                save_json(receipt_path, receipt)
                receipts.append(receipt)
                print('[HERV export] '+chrom+' '+json.dumps(counts), flush=True)
            finally:
                pool.shutdown(wait=True, cancel_futures=True)
                for handle in handles.values(): handle.close()
                for path in temps.values(): path.unlink(missing_ok=True)
        manifest = dict(schema=SCHEMA, status='complete', build='GRCh38', coordinate_convention='0-based half-open',
            samples=len(samples), haplotypes=2*len(samples), chromosomes=list(CHROMS), fields=FIELDS,
            raw_coordinates='personalized haplotype coordinates; NOT reference BED coordinates',
            segments_coordinates='GRCh38 homologous reference coordinates, masked bases excluded',
            uncallable_coordinates='GRCh38; unknown spans are NOT absence',
            families=list(FAMILIES), order='sample, haplotype, family, start0; not globally coordinate-sorted',
            interpretation='Variant-conditioned reference reconstruction, not de novo insertion discovery; inserted bases without reference homology are omitted from segments.',
            sequence_annotation_settings=settings,
            files=[dict(chrom=r['chrom'],kind=k,**v) for r in receipts for k,v in r['files'].items()])
        save_json(manifest_path, manifest)
        return dict(status='complete', manifest=str(manifest_path), samples=len(samples), chromosomes=len(CHROMS))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    print(json.dumps(export_cohort(parser.parse_args().out), indent=2))
