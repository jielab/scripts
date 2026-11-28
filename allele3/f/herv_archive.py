"""Validate permanent HERV products and record reuse in the analysis runtime."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys


def read(path):
    return json.loads(Path(path).read_text())


def validate(args):
    from herv_cohort import requested_populations
    from workflow_paths import configured_path
    root=args.archive.expanduser().resolve()
    metadata=args.out.expanduser().resolve()/'qc/archive'
    manifest=read(metadata/'export_manifest.json');summary=read(metadata/'summary.json');plan=read(metadata/'plan.json')
    identity=plan['identity'];samples=plan['samples']
    chroms=['chr'+str(i) for i in range(1,23)]
    if manifest.get('schema')!='herv-chromosome-export-2' or manifest.get('status')!='complete' or manifest.get('build')!='GRCh38' or manifest.get('chromosomes')!=chroms:
        raise ValueError('Incomplete GRCh38 chromosome export')
    from herv_export import FIELDS
    if manifest.get('fields') != FIELDS+['record_type'] or manifest.get('record_types') != ['segments','uncallable']:
        raise ValueError('Missing reference record types or incompatible table fields')
    if summary.get('status')!='complete' or summary.get('missing_samples') or not samples:
        raise ValueError('Incomplete cohort summary')
    n=len(samples)
    if any(x!=n for x in [manifest.get('samples'),summary.get('planned_samples'),summary.get('completed_samples')]) or manifest.get('haplotypes')!=2*n:
        raise ValueError('Sample counts disagree')
    if identity.get('build')!='GRCh38' or identity.get('seed')!=args.seed or identity.get('per_pop')!=args.per_pop or identity.get('populations')!=requested_populations(args.populations):
        raise ValueError('Requested cohort/seed differs from archive')
    if str(args.pfile_dir.expanduser().resolve())!=identity.get('pfile_dir'):
        raise ValueError('Requested genotype directory differs from archive')
    # Old archives lack a reusable input-location receipt for custom resources.
    # Never silently accept a different reference, library or supplied catalogue.
    for key,expected in [('fasta',configured_path('DIR_GEN','1hgp/38/GRCH38.fasta.gz')),
                         ('library',configured_path('DIR_ANNOT','dfam_4.0_human/human.curated.fa')),
                         ('sample_panel',configured_path('DIR_GEN','1kg/38/samples.txt'))]:
        if getattr(args,key).expanduser().resolve()!=Path(expected).expanduser().resolve():
            raise ValueError('Custom '+key+' requires native validation')
    if args.reference_bed or args.repeatmasker or args.pgen_python or args.allow_partial:
        raise ValueError('Custom reference/tool or partial-summary request requires native validation')
    settings=summary['settings'];sequence=summary['sequence_annotation_settings']
    for key in ('min_callable','min_region_coverage','stable_retention','stable_min_sample_fraction','stable_min_population_fraction','bootstrap','seed'):
        if settings.get(key)!=getattr(args,key):raise ValueError('Requested setting differs: '+key)
    for key in ('chunk_bp','halo'):
        if sequence.get(key)!=getattr(args,key):raise ValueError('Requested sequence setting differs: '+key)
    if sequence!=manifest.get('sequence_annotation_settings'):raise ValueError('Export and summary settings differ')
    with (metadata/'samples.tsv').open() as f:ids=[r['sample'] for r in csv.DictReader(f,delimiter='\t')]
    if ids!=identity.get('sample_ids') or ids!=[r['sample'] for r in samples] or len(set(ids))!=n:
        raise ValueError('Sample identities disagree')
    checks={}
    for line in (metadata/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split(None,1);name=name.lstrip(' *')
        base=root if re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])\.tsv\.gz',name) else metadata
        path=base/name
        if not re.fullmatch('[0-9a-f]{64}',digest) or Path(name).is_absolute() or '..' in Path(name).parts or not path.resolve().is_relative_to(base) or name in checks:
            raise ValueError('Invalid archive checksum entry')
        checks[name]=digest
    required={'export_manifest.json','plan.json','summary.json','samples.tsv','projection_summary.tsv',
              'sample_completion.tsv','sample_summary.tsv','population_summary.tsv',
              'locus_calls.npz','locus_stability.tsv.gz','stable_loci.bed','reference_herv.bed',
              'family_stability.tsv','sample_family_stability.tsv','pairwise_jaccard.tsv','stability_bootstrap.tsv'}
    expected={(c,'reference') for c in chroms}
    found=set()
    for item in manifest['files']:
        pair=(item['chrom'],item['kind'])
        if pair not in expected or pair in found:raise ValueError('Duplicate/unexpected chromosome table')
        found.add(pair);name=item['path'];required.add(name)
        if name != item['chrom']+'.tsv.gz' or set(item.get('record_counts',{}))!={'segments','uncallable'} or sum(item['record_counts'].values())!=item['rows']:
            raise ValueError('Invalid combined reference table receipt')
        if checks.get(name)!=item['sha256'] or (root/name).stat().st_size!=item['bytes']:
            raise ValueError('Export file size/checksum receipt mismatch: '+name)
    if found!=expected or not required.issubset(checks):raise ValueError('Archive lacks required products/checksums')
    print('[HERV archive] checking SHA256: '+str(root),flush=True)
    for name,expected_hash in checks.items():
        base=root if re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])\.tsv\.gz',name) else metadata
        p=base/name;before=p.stat();h=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
        after=p.stat()
        if h.hexdigest()!=expected_hash or (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
            raise ValueError('Archive content changed or checksum failed: '+name)
    return dict(status='reused',archive=str(root),samples=n,chromosomes=22,chromosome_tables=22,
                verified_files=len(checks),validation='full SHA256, completion, cohort and scientific settings',
                sample_jobs_submitted=0,repeatmasker_runs=0)


def try_reuse(args):
    if args.command not in {'run','summarize','export'} or args.replace=='TRUE' or os.environ.get('ALLELE_REPLACE')=='TRUE':return None
    from workflow_paths import configured_path
    default_out=Path(os.environ.get('HERV_OUT',os.environ.get('ALLELE_HERV_SCAN_STORE',configured_path('ANALYSIS_ROOT','allele3/2.herv_scan'))))
    if args.out.expanduser().resolve()!=default_out.expanduser().resolve():return None
    if args.sample or args.sample_index:return None
    try:result=validate(args)
    except FileNotFoundError:
        # First runs and deleted analysis receipts simply require normal analysis.
        return None
    except (OSError,ValueError,KeyError,TypeError) as exc:
        print('[HERV archive] not reusable: '+str(exc)+'; continuing analysis workflow',flush=True)
        return None
    from herv_publish import hydrate_runtime
    hydrate_runtime(args.out, args.archive, refresh=True)
    from herv_cohort import save_json
    save_json(args.out/'qc/archive_reuse.json',result)
    save_json(args.out/'qc/publication.json',result)
    print('[HERV archive] REUSE '+json.dumps(result),flush=True)
    return result


def main(argv=None):
    from herv_cohort import parser
    argv=list(sys.argv[1:] if argv is None else argv)
    # Shell-only execution controls do not alter scientific archive identity.
    filtered=[];i=0
    while i<len(argv):
        value=argv[i];key=value.split('=',1)[0]
        if key in {'--python','--memory-limit-gb','--screen-log-interval'}:
            i+=1 if '=' in value else 2;continue
        if key=='--verbose-screen':i+=1;continue
        filtered.append(value);i+=1
    args=parser().parse_args(filtered)
    return 0 if try_reuse(args) is not None else 1


if __name__=='__main__':sys.exit(main())
