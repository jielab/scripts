#!/usr/bin/env python3
"""Summarize split read counts and retain sample-level reference evidence."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import sys

from common import DIR_OUT, DIR_REF, DIR_SPLIT, atomic_json
from library import digest
from sample_state import is_complete


def write_table(path, fields, rows):
    partial=path.with_name(path.name+'.partial')
    with partial.open('w', newline='', encoding='utf-8') as handle:
        writer=csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def snapshot_references(gather, output, library):
    """Join actual gather matches to the reference inventory, never assign reads."""
    with Path(gather).open(newline='') as handle:
        reader=csv.DictReader(handle)
        fields=reader.fieldnames or []
        matches=list(reader)
    name_field='match_name' if 'match_name' in fields else 'name'
    if name_field not in fields:
        raise RuntimeError(f'Sourmash gather result lacks reference names: {gather}')
    inventory={}
    if matches:
        with (library/'reference/inventory.tsv').open(newline='') as handle:
            inventory={row['filename']:row for row in csv.DictReader(handle,delimiter='\t')}
    columns=['evidence_type','assembly_accession','reference_file','reference_path',
             'reference_source_group','taxid','species_taxid','organism_name']
    rows=[]
    for match in matches:
        name=match[name_field]
        if name not in inventory:
            raise RuntimeError(f'Gather reference absent from library inventory: {name}')
        reference=inventory[name]
        relative=Path(reference['relative_path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise RuntimeError(f'Unsafe reference inventory path: {relative}')
        row=dict(evidence_type='sample_sourmash_gather_candidate',
                 assembly_accession=reference['assembly_accession'],reference_file=name,
                 reference_path=str(DIR_REF/relative),reference_source_group=reference['source_group'],
                 taxid=reference['taxid'],species_taxid=reference['species_taxid'],
                 organism_name=reference['organism_name'])
        row.update({'sourmash_'+field:match[field] for field in fields})
        rows.append(row)
    write_table(output,columns+['sourmash_'+field for field in fields],rows)


def summarize(sample, output=None):
    sample=Path(sample).resolve()
    receipt_path=sample/'split.done'
    try:
        receipt=json.loads(receipt_path.read_text())
    except (OSError,ValueError) as exc:
        raise RuntimeError(f'Completed split receipt missing/invalid: {receipt_path}') from exc
    if not isinstance(receipt,dict) or not receipt.get('input_identity') or not is_complete(sample,receipt['input_identity']):
        raise RuntimeError(f'Split output is incomplete or has changed: {sample}')
    with (sample/'counts.tsv').open(newline='') as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    species=Counter()
    groups=Counter()
    total_fragments=total_reads=0
    multiplier=2 if receipt['paired'] else 1
    for row in rows:
        fragments,reads=int(row['fragments']),int(row['reads'])
        if fragments<0 or reads!=fragments*multiplier:
            raise RuntimeError('Inconsistent fragment/read counts in counts.tsv')
        key=(row['destination'],row['group'],row['species_taxid'],row.get('species_name',''),row['reason'])
        species[key]+=fragments
        groups[row['group']]+=fragments
        total_fragments+=fragments
        total_reads+=reads
    if total_fragments!=receipt['fragments']:
        raise RuntimeError('counts.tsv total differs from split completion receipt')
    evidence=sample/'reference-matches.tsv'
    has_sample_evidence='reference-matches.tsv' in receipt['outputs']
    if output is None:
        try:
            relative=sample.relative_to(DIR_SPLIT.resolve())
        except ValueError:
            relative=Path(sample.name+'-'+digest(str(sample))[:8])
        output=DIR_OUT/'composition'/relative
    output=Path(output).resolve()
    if output==sample or sample in output.parents:
        raise RuntimeError('Report output must be outside the split sample directory')
    output.mkdir(parents=True,exist_ok=True)

    def amounts(fragments):
        reads=fragments*multiplier
        fraction=reads/total_reads if total_reads else 0
        return dict(fragments=fragments,reads=reads,fraction_all_reads=fraction,percent_all_reads=100*fraction)

    species_rows=[dict(destination=key[0],group=key[1],species_taxid=key[2],species_name=key[3],reason=key[4],**amounts(n))
                  for key,n in sorted(species.items(),key=lambda pair:(-pair[1],pair[0]))]
    group_rows=[dict(group=group,**amounts(n)) for group,n in sorted(groups.items(),key=lambda pair:(-pair[1],pair[0]))]
    metrics=['fragments','reads','fraction_all_reads','percent_all_reads']
    write_table(output/'species.tsv',['destination','group','species_taxid','species_name','reason']+metrics,species_rows)
    write_table(output/'groups.tsv',['group']+metrics,group_rows)
    statuses=receipt.get('reference_status_fragments',{})
    write_table(output/'reference-status.tsv',['reference_status']+metrics,
                [dict(reference_status=status,**amounts(n)) for status,n in sorted(statuses.items())])
    references=None
    if has_sample_evidence:
        shutil.copyfile(evidence,output/'reference-matches.tsv.partial')
        (output/'reference-matches.tsv.partial').replace(output/'reference-matches.tsv')
        with evidence.open(newline='') as handle:
            references=sum(1 for _ in csv.DictReader(handle,delimiter='\t'))
    else:
        write_table(output/'reference-matches.tsv',
                    ['evidence_type','assembly_accession','reference_file','reference_path',
                     'reference_source_group','taxid','species_taxid','organism_name'], [])
    atomic_json(output/'summary.json',dict(sample=str(sample),input_identity=receipt['input_identity'],
                fragments=total_fragments,reads=total_reads,reference_candidates=references,
                read_assignments=str(sample/'read-assignments.tsv.gz'),
                read_reference_hits=str(sample/'read-reference-hits.tsv.gz'),
                reference_status_fragments=statuses,
                denominator='all input reads, including unclassified and unresolved reads',
                sample_reference_evidence=('Sample-level sourmash gather candidates (default threshold 50 kb)'
                                           if has_sample_evidence else 'not_recorded'),
                read_reference_evidence=receipt.get('reference_evidence','not recorded')))
    print(f'Composition: {total_reads:,} reads ({total_fragments:,} fragments) -> {output}')
    print(f'Sample sketch candidates: {references if has_sample_evidence else "not recorded"}; '
          f'per-read reference evidence: {receipt.get("reference_evidence", "not recorded")}')
    return output


def main():
    parser=argparse.ArgumentParser(prog='./metag.sh composition',description='Summarize species read fractions and reference candidates from completed split data')
    parser.add_argument('--sample',type=Path,required=True,help='completed sample directory containing split.done')
    parser.add_argument('--out',type=Path,help='report directory (default: analysis/metag/composition/<sample>)')
    args=parser.parse_args()
    summarize(args.sample,args.out)
    return 0


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (RuntimeError,OSError,ValueError,KeyError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        raise SystemExit(1)

