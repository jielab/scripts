"""Publish one reference-coordinate table per chromosome; keep runtime in analysis."""
from pathlib import Path
import csv
import gzip
import json
import os
import shutil
import uuid
from herv_layout import product_path, ensure_layout

SCHEMA = 'herv-chromosome-export-2'
PRODUCTS = ('plan.json', 'summary.json', 'samples.tsv', 'projection_summary.tsv',
            'sample_completion.tsv', 'sample_summary.tsv', 'population_summary.tsv',
            'locus_calls.npz', 'locus_stability.tsv.gz', 'stable_loci.bed',
            'reference_herv.bed', 'family_stability.tsv', 'sample_family_stability.tsv',
            'pairwise_jaccard.tsv', 'stability_bootstrap.tsv')


def prune_duplicate_pdfs(out):
    for pdf in (Path(out)/'figures').glob('*.pdf'):
        if pdf.with_suffix('.png').is_file():
            pdf.unlink()


def archive_root():
    from workflow_paths import configured_path
    return Path(os.environ.get('ALLELE_HERV_ARCHIVE', configured_path('DIR_GEN', '1kg/38/herv_scan')))


def hydrate_runtime(out, archive=None, *, refresh=False):
    """Initialize analysis directories; restore missing inputs from analysis receipts."""
    out, archive = Path(out).resolve(), Path(archive or archive_root()).resolve()
    if out == archive or archive in out.parents or out in archive.parents:
        raise ValueError('HERV runtime and permanent archive must be separate directories')
    ensure_layout(out)
    metadata = out/'qc/archive'
    if not (metadata/'export_manifest.json').is_file():
        return
    manifest = json.loads((metadata/'export_manifest.json').read_text())
    if manifest.get('status') != 'complete':
        raise ValueError('HERV archive is incomplete')
    for name in PRODUCTS:
        if (refresh or not product_path(out,name).exists()) and (metadata/name).is_file():
            shutil.copy2(metadata/name, product_path(out,name))
    summary_path = product_path(out,'summary.json')
    if summary_path.is_file():
        from herv_cohort import save_json
        summary = json.loads(summary_path.read_text())
        summary.update(stable_bed=str(product_path(out,'stable_loci.bed')),
                       reference_bed=str(product_path(out,'reference_herv.bed')))
        save_json(summary_path,summary)


def publish(out, archive):
    from herv_export import sha, FIELDS
    from herv_cohort import file_lock, save_json
    out, archive = Path(out).resolve(), Path(archive).resolve()
    if out == archive or archive in out.parents or out in archive.parents:
        raise ValueError('HERV runtime and permanent archive must be separate directories')
    ensure_layout(out)
    archive.parent.mkdir(parents=True, exist_ok=True)
    with file_lock(out/'tmp/publish.lock'):
        source = json.loads((out/'qc/export_manifest.json').read_text())
        summary = json.loads((out/'qc/summary.json').read_text())
        if source.get('status') != 'complete' or summary.get('status') != 'complete':
            raise ValueError('Only a complete cohort can be published')
        stage = archive.with_name('.'+archive.name+'.staging-'+uuid.uuid4().hex)
        stage.mkdir()
        metadata = out/'tmp'/('archive-metadata-'+uuid.uuid4().hex)
        metadata.mkdir()
        try:
            files = []
            for chrom in source['chromosomes']:
                counts = {}
                target = stage/(chrom+'.tsv.gz')
                with gzip.open(target, 'wt', compresslevel=3, newline='') as dest:
                    writer = csv.writer(dest, delimiter='\t', lineterminator='\n')
                    writer.writerow(FIELDS+['record_type'])
                    for kind in ('segments', 'uncallable'):
                        item = next(r for r in source['files'] if r['chrom']==chrom and r['kind']==kind)
                        path = out/item['path']
                        if sha(path) != item['sha256']:
                            raise ValueError('Source checksum mismatch: '+str(path))
                        count = 0
                        with gzip.open(path, 'rt', newline='') as handle:
                            reader = csv.reader(handle, delimiter='\t')
                            if next(reader) != FIELDS:
                                raise ValueError('Unexpected table header: '+str(path))
                            for row in reader:
                                if len(row)!=len(FIELDS) or row[0]!=chrom:
                                    raise ValueError('Invalid chromosome row')
                                writer.writerow(row+[kind]); count += 1
                        if count != item['rows']:
                            raise ValueError('Source row count mismatch: '+str(path))
                        counts[kind] = count
                files.append(dict(chrom=chrom, kind='reference', path=target.name,
                                  rows=sum(counts.values()), record_counts=counts,
                                  bytes=target.stat().st_size, sha256=sha(target)))
                print('[HERV publish] '+chrom+' '+str(counts), flush=True)
            for name in PRODUCTS:
                shutil.copy2(product_path(out,name), metadata/name)
            # Freeze scientific identity in analysis, separate from mutable enrichment status.
            summary = {k:v for k,v in summary.items() if k not in ('export','enrichment','report')}
            summary.update(stable_bed='stable_loci.bed', reference_bed='reference_herv.bed')
            save_json(metadata/'summary.json', summary)
            manifest = {k:v for k,v in source.items() if k not in ('files','raw_coordinates','order')}
            manifest.update(schema=SCHEMA, fields=FIELDS+['record_type'], files=files,
                record_types=['segments','uncallable'],
                order='record_type, sample, haplotype, family, start0; not globally coordinate-sorted',
                runtime=str(out))
            save_json(metadata/'export_manifest.json', manifest)
            checks = [(sha(p),p.relative_to(stage).as_posix()) for p in sorted(stage.rglob('*')) if p.is_file()]
            checks += [(sha(p),p.relative_to(metadata).as_posix()) for p in sorted(metadata.rglob('*')) if p.is_file()]
            (metadata/'SHA256SUMS').write_text(''.join(h+'  '+p+'\n' for h,p in checks))
            # Stage and validate all files before replacing the visible archive.
            previous = out/'tmp'/('archive-before-'+uuid.uuid4().hex)
            if archive.exists():
                shutil.move(str(archive), str(previous))
            receipt = out/'qc/archive'
            receipt.parent.mkdir(parents=True,exist_ok=True)
            old_receipt = out/'tmp'/('metadata-before-'+uuid.uuid4().hex)
            try:
                stage.rename(archive)
                if receipt.exists():receipt.rename(old_receipt)
                metadata.rename(receipt)
            except BaseException:
                if archive.exists():shutil.move(str(archive),str(out/'tmp'/('failed-publication-'+uuid.uuid4().hex)))
                if previous.exists():shutil.move(str(previous),str(archive))
                if old_receipt.exists() and not receipt.exists():old_receipt.rename(receipt)
                raise
            result = dict(status='complete', archive=str(archive), chromosome_tables=len(files), schema=SCHEMA)
            save_json(out/'qc/publication.json',result)
            prune_duplicate_pdfs(out)
            return result
        finally:
            if stage.exists():
                shutil.rmtree(stage)
            if metadata.exists():
                shutil.rmtree(metadata)
