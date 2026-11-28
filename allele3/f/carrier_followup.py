"""Measure target allele dosage within individuals with finalized GU calls."""
from collections import defaultdict
import csv
import json
import os
from pathlib import Path
import shutil

from research_common import atomic_json, write_tsv

def followup(products, gu, out):
    from archaic_alleles import metadata, metadata_vcfs, parser, run
    from herv import archaic_inspect_dataset
    from workflow_paths import configured_path
    from analysis_integrity import resolve_result_input
    out.mkdir(parents=True,exist_ok=True)
    status=out/'availability.json'
    def unavailable(reason):
        atomic_json(status,{'status':'unavailable','reason':reason})
    atomic_json(status,{'status':'running'})
    herv=products/'herv_enrichment'
    candidates=[herv/'artifacts/runtime/candidates/targets.tsv',herv/'results/candidates.targets.tsv',
                Path(resolve_result_input(str(herv/'candidates/targets.tsv')))]
    target=next((p for p in candidates if p.is_file()),None)
    if not target:return unavailable('Detailed HERV candidate panel unavailable')
    files=sorted((gu/'ibdmix').glob('chr*/run.meta.tsv'))
    if not files:return unavailable('GU chromosome metadata unavailable')
    first=metadata(files[0]);build=first.get('genome_build','').removeprefix('GRCh').removeprefix('b')
    # Restrict to sources with native evidence on all requested autosomes.
    try:
        references=os.environ.get('ALLELE_GU_REFERENCES') or first.get('references') or 'Altai,Chagyr,Denisova,Denisova25,Vindija'
        references=references.replace(' ',',')
        refs=[v for v in references.split(',') if v]
        info=archaic_inspect_dataset(gu,'ibdmix',tuple(str(i) for i in range(1,23)),refs)
        for m in files:
            modern,ancient=metadata_vcfs(m,None,None,m.parent.name.removeprefix('chr'),refs)
            for vcf in [modern,*ancient.values()]:
                if not Path(vcf).is_file() or not any(Path(str(vcf)+ext).is_file() for ext in ('.tbi','.csi')):
                    return unavailable('Source VCF or index unavailable: '+str(vcf))
    except (ValueError,FileNotFoundError) as exc:return unavailable(str(exc))
    args=['--gu-dataset',str(gu),'--modern-build',build,'--targets',str(target),'--out',str(out/'artifacts'),
          '--references',','.join(refs),'--min-lod',str(info['native_min_lod']),'--min-length',str(info['native_min_length'])]
    if build=='37':
        chain=os.environ.get('ALLELE_CHAIN38TO37') or configured_path('DIR_FILES','liftOver/hg38ToHg19.over.chain.gz')
        exe=shutil.which(os.environ.get('ALLELE_ARCHAIC_LIFTOVER','liftOver')) or shutil.which('liftOver')
        if not Path(chain).is_file() or not exe:return unavailable('GRCh38 target-to-GU chain or liftOver missing')
        args+=['--target-to-gu-chain',chain,'--liftover',exe]
    run(parser().parse_args(args))
    groups=defaultdict(lambda:[0,0,0,0])
    with (out/'artifacts/allele_tract_summary.tsv').open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=tuple(row[k] for k in ('reference','allele_role','tract_call'))
            totals=groups[key]
            totals[0]+=1
            for i,k in enumerate(('n_called','AC','AN'),1):totals[i]+=int(row[k])
    fields=['reference','allele_role','tract_call','site_allele_rows','called_person_site_records','AC','AN','pooled_AF']
    table=[(*key,*v,v[2]/v[3] if v[3] else None) for key,v in sorted(groups.items())]
    write_tsv(out/'results/role_by_tract.tsv',fields,table)
    atomic_json(status,{'status':'complete','target_file':str(target),'rows':len(table),'interpretation':'Repeated person-site records; unphased IBDmix carrier overlap.'})
