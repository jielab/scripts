"""Schemes shared by the per-sample evidence job and the manuscript report."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFINE_FORMAT = 'metag-refine-v2'

SCHEMA = {
    'refine_sources': 'sample group rank taxid name bin_pairs food_source split_evidence',
    'refine_blast': 'sample split_taxid split_name split_rank bin_pairs sampled_pairs reference_scope outcome comparison_taxid comparison_rank comparison_name pairs',
    'refine_amr': 'sample family gene_annotation annotation_group true_markers uniquely_supported_markers assigned_pairs both_mates_pairs input_pairs pairs_per_million',
    'refine_virus': 'sample category species_bin_pairs input_pairs library_type interpretation',
    'refine_status': 'sample state sample_type library_type input_pairs blast amr blast_db_id amr_marker_sha256 blast_pairs_per_bin amr_engine amr_tested_families amr_ambiguous_pairs split_sha256 refine_version',
}
SCHEMA = {k: v.split() for k,v in SCHEMA.items()}

def code_id():
    h = hashlib.sha256()
    for p in sorted(HERE.glob('*.py')) + [HERE/'food_sources.tsv']:
        h.update(p.name.encode() + b'\0' + p.read_bytes() + b'\0')
    return h.hexdigest()

def rows(path):
    with path.open(encoding='utf-8', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')

def write(path, fields, data):
    with path.open('w', encoding='utf-8', newline='') as handle:
        w = csv.writer(handle, delimiter='\t', lineterminator='\n')
        w.writerow(fields)
        w.writerows(data)

def collect(root, sample, split_sha, writers):
    from common import locked, sha256
    loc = root/'main'/sample
    def missing(state):
        row = dict(sample=sample, state=state, split_sha256=split_sha)
        writers['refine_status'].writerow([row.get(k,'NA') for k in SCHEMA['refine_status']])
    lock = loc.parent/('.'+loc.name+'.lock')
    if not lock.exists():
        missing('not_run')
        return
    with locked(lock, shared=True):
        if not (loc/'receipt.json').is_file():
            missing('not_run')
            return
        rec = json.loads((loc/'receipt.json').read_text())
        if rec['sample'] != sample or rec['split_sha256'] != split_sha:
            missing('split_changed')
            return
        if rec.get('format') != REFINE_FORMAT:
            raise ValueError('Expected current refine receipt in '+str(loc))
        if rec.get('state')!='complete':
            missing('incomplete')
            return
        for table, fields in SCHEMA.items():
            if rec.get('table_hashes', {}).get(table)!=sha256(loc/(table+'.tsv')):
                raise ValueError('Refine output differs from its receipt: '+str(loc/(table+'.tsv')))
            with (loc/(table+'.tsv')).open(encoding='utf-8',newline='') as h:
                r = csv.DictReader(h, delimiter='\t')
                if r.fieldnames != fields:
                    raise ValueError('Incompatible refine table: '+str(loc/table))
                for row in r:
                    if row['sample'] != sample:
                        raise ValueError('Refine sample mismatch: '+str(loc/table))
                    writers[table].writerow([row[k] for k in fields])

def summarize(folder):
    from collections import Counter, defaultdict
    status = {r['sample']:r for r in rows(folder/'refine_status.tsv')}
    denominators, input_pairs = Counter(),Counter()
    for sample,r in status.items():
        if r['state']=='complete' and int(r['input_pairs'])>0:
            cohort=(sample.split('/')[0],r['sample_type'],r['library_type'])
            denominators[cohort]+=1
            input_pairs[cohort]+=int(r['input_pairs'])
    totals, prevalence = Counter(),Counter()
    for r in rows(folder/'refine_sources.tsv'):
        sample=r['sample']
        cohort=(sample.split('/')[0],status[sample]['sample_type'],status[sample]['library_type'])
        key=(*cohort,r['group'],r['rank'],r['taxid'],r['name'],r['food_source'])
        if int(r['bin_pairs'])>0:
            prevalence[key]+=1
            totals[key]+=int(r['bin_pairs'])
    source_rows = []
    for key,n in totals.most_common():
        cohort=key[:3]
        source_rows.append([*key,n,prevalence[key],denominators[cohort],input_pairs[cohort],
                            1e6*n/input_pairs[cohort]])
    write(folder/'refine_sources_by_dataset.tsv',
          'dataset sample_type library_type group rank taxid name food_source bin_pairs samples_with_bin completed_nonempty_samples input_pairs pairs_per_million'.split(),source_rows)
    audit = defaultdict(Counter)
    outcomes = ['supports_bin_rank','higher_rank_only','different_lineage','insufficient','hit_limit','missing_taxonomy','panel_not_evaluable']
    for r in rows(folder/'refine_blast.tsv'):
        key = tuple(r[k] for k in ('sample','split_taxid','split_name','split_rank','bin_pairs','reference_scope','sampled_pairs'))
        if r['outcome'] not in outcomes:
            raise ValueError('Unknown BLAST audit outcome')
        audit[key][r['outcome']] += int(r['pairs'])
    for key,c in audit.items():
        if sum(c.values())!=int(key[-1]):
            raise ValueError('BLAST audit denominator mismatch: '+str(key))
    write(folder/'refine_blast_by_bin.tsv',
          ['sample','split_taxid','name','rank','bin_pairs','reference_scope','sampled_pairs',*outcomes],
          ([*key,*[c[x] for x in outcomes]] for key,c in sorted(audit.items())))
    return source_rows
