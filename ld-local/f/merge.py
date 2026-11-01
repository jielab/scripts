#!/usr/bin/env python3
"""Merge only completed, compatible per-chromosome jobs; never compute LD."""
import os
from ld_pipeline import Pipeline, POPS, Stop, atomic, read_json, sha, write_json


def merge():
    # Final publication always covers all autosomes, even after a subset run.
    os.environ['CHRS'] = '1-22'
    pipeline = Pipeline()
    missing = [str(c) for c in pipeline.chroms
               if not (pipeline.out / f'chr{c}' / 'rebuild_v2' / 'COMPLETE.json').is_file()]
    if missing:
        raise Stop('chr1-22 results are required; missing completion records: ' + ', '.join(missing))
    with pipeline.output_locks():
        pipeline.work.mkdir(parents=True, exist_ok=True)
        (pipeline.work / 'COMPLETE.json').unlink(missing_ok=True)
        texts = {p: [] for p in ('ALL', *POPS)}
        common = None
        evidence = {}
        warnings = []
        for chrom in pipeline.chroms:
            source = pipeline.out / f'chr{chrom}' / 'rebuild_v2'
            receipt = source / 'COMPLETE.json'
            record = read_json(receipt)
            if not record or record.get('chromosomes') != [chrom] or record.get('status') not in ('complete', 'complete_with_reference_warnings'):
                raise Stop(f'chr{chrom}: missing successful single-chromosome completion record')
            if record.get('genome_build') != f'GRCh{pipeline.genome_build}' or record.get('panel_sha256') != sha(pipeline.panel):
                raise Stop(f'chr{chrom}: build or sample panel mismatch')
            signature = {k: record[k] for k in ('version', 'genome_build', 'input_format', 'panel_sha256', 'sample_counts', 'sample_selection',
                                               'coincident_policy', 'beta_bp', 'plink2', 'universal_convention',
                                               'reference_count_policy', 'reference_max_rel_deviation')}
            if common is not None and signature != common:
                raise Stop(f'chr{chrom}: incompatible analysis settings')
            common = signature
            selection = record['distance_selection'][str(chrom)]
            if sha(source / selection['path']) != selection['sha256']:
                raise Stop(f'chr{chrom}: selection receipt changed')
            for pop in texts:
                path = source / f'{pop}.{pipeline.genome_build}.SUBSET.bed'
                if sha(path) != record['files'].get(path.name):
                    raise Stop(f'chr{chrom}: BED checksum mismatch: {path}')
                texts[pop].append(path.read_text())
            evidence[str(chrom)] = {'path': str(receipt.relative_to(pipeline.out)), 'sha256': sha(receipt)}
            if record['status'] == 'complete_with_reference_warnings':
                warnings.append(chrom)
        suffix = '' if pipeline.chroms == list(range(1, 23)) else '.SUBSET'
        files = {}
        for pop, parts in texts.items():
            path = pipeline.work / f'{pop}.{pipeline.genome_build}{suffix}.bed'
            atomic(path, ''.join(parts))
            files[path.name] = sha(path)
        write_json(pipeline.work / 'COMPLETE.json', dict(common,
                   chromosomes=pipeline.chroms, chromosome_jobs=evidence, files=files,
                   status='complete_with_reference_warnings' if warnings else 'complete',
                   reference_count_warning_chromosomes=warnings,
                   independent_reimplementation=True, coordinate_identity_to_author_files_verified=False))
        print(f'[LD] MERGED: {pipeline.work}')


if __name__ == '__main__':
    try:
        merge()
    except (Stop, OSError, KeyError, ValueError) as exc:
        raise SystemExit(f'[LD] MERGE STOP: {exc}')
