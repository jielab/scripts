#!/usr/bin/env python3
"""Merge completed chromosome MFI audits into one exact PVAR-bound INFO table."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
from analysis_integrity import InputPath,identity,transactional,validate_manifest


def checked_sources(directories):
    result=[];chroms=set();contract=None
    for directory in directories:
        directory=Path(directory)
        m=validate_manifest(directory,require_evidence=True)
        if m.get('schema')!='ukb-mfi-1.0' or m.get('genome_build')!='GRCh37' or m.get('status') not in {'complete','no_usable_exact_high_quality_contrast'}:
            raise ValueError('Not a completed GRCh37 MFI audit: '+str(directory))
        params=m.get('parameters',{})
        current={k:params.get(k) for k in ('mfi_layout','min_info','min_maf')}
        if contract is None:contract=current
        elif current!=contract:raise ValueError('Chromosome MFI quality thresholds/layout differ')
        source_files=m.get('source_files',[])
        for source in source_files:
            c=source['chrom']
            if c in chroms: raise ValueError('Overlapping chromosome audits')
            chroms.add(c)
            for key in ('mfi','pvar'):
                saved=m.get('inputs',{}).get(str(Path(source[key]).resolve()))
                if not saved or identity(source[key],full='sha256' in saved)!=saved:
                    raise ValueError('Original source differs from its MFI audit: '+source[key])
        result.append((directory,m))
    if not result:raise ValueError('No MFI audit directories')
    return sorted(result,key=lambda pair:min(int(s['chrom']) for s in pair[1]['source_files']))


def inputs(raw):
    result=[]
    for directory,m in checked_sources(raw['directories']):
        result.append(InputPath(directory,'directory'))
        result += [InputPath(s[key]) for s in m['source_files'] for key in ('mfi','pvar')]
    return result


@transactional('manifest.json',resolve_inputs=inputs,input_prefixes=())
def merge(directories,out):
    import pysam
    sources=checked_sources(directories);out=Path(out);counts=Counter();seen=set();last=None
    target=out/'target_info.tsv.gz';n=0
    with pysam.BGZFile(str(target),'wb') as handle:
        handle.write(b'chrom\tpos\tref\talt\tinfo\n')
        for directory,m in sources:
            counts.update(m.get('counts',{}))
            with gzip.open(directory/'target_info.tsv.gz','rt',newline='') as inp:
                reader=csv.DictReader(inp,delimiter='\t')
                if reader.fieldnames!=['chrom','pos','ref','alt','info']:raise ValueError('Unexpected quality table schema')
                for row in reader:
                    key=(int(row['chrom']),int(row['pos']))
                    if last is not None and key<=last:raise ValueError('Quality records overlap or are not globally sorted')
                    last=key;n+=1
                    handle.write(('\t'.join(row[c] for c in reader.fieldnames)+'\n').encode())
    if n:pysam.tabix_index(str(target),seq_col=0,start_col=1,end_col=1,line_skip=1,zerobased=False,force=True)
    if n!=counts['quality_pass']:raise ValueError('Merged quality row count does not reconcile with input audits')
    return dict(status='complete' if n else 'no_usable_exact_high_quality_contrast',schema='ukb-mfi-1.0',
        genome_build='GRCh37',counts=dict(counts),mfi_layout=sources[0][1]['mfi_layout'],
        parameters={k:sources[0][1]['parameters'][k] for k in ('mfi_layout','min_info','min_maf')},
        source_files=[s for _,m in sources for s in m['source_files']],
        chromosome_manifests=[identity(d/'manifest.json',full=True) for d,_ in sources],
        tabix_index='created' if n else 'not_created_no_pass_rows',
        evidence='Exact PVAR-bound imputation quality only; no genotype counts or WGS event claims',
        limits=['Missing contrasts remain excluded/unknown. Original MFI/PVAR inputs must match each source audit.'])


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directories-json',required=True);p.add_argument('--out',required=True)
    a=p.parse_args(argv);m=merge(json.loads(a.directories_json),a.out)
    print(json.dumps(dict(status=m['status'],counts=m['counts']),sort_keys=True));return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,RuntimeError,KeyError) as exc:raise SystemExit('mfi-merge: '+str(exc))
