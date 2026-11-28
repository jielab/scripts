"""Render a portable explorer from published TSVs. Counts are always read pairs."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

def rows(folder,name):
    p = folder/(name+'.tsv')
    if p.is_file():
        with p.open(encoding='utf-8',newline='') as f:
            yield from csv.DictReader(f,delimiter='\t')

def projection(folder,name,keys):
    return [{k:r[k] for k in keys.split()} for r in rows(folder,name)]

def write(folder, preview=False):
    samples=list(rows(folder,'samples'))
    if not samples:
        raise ValueError('Explorer requires at least one completed sample.')
    by_sample={r['sample']:i for i,r in enumerate(samples)}
    taxa={}
    catalog=[]
    bins=[[] for _ in samples]
    counts=defaultdict(Counter)
    by_cohort=defaultdict(list)
    for r in samples:
        by_cohort[r['sample'].split('/')[0]].append(r)
    for r in rows(folder,'taxa'):
        sid=r['sample']
        key=(r['taxid'],r['name'],r['group'],r['rank'])
        if key not in taxa:
            taxa[key]=len(catalog)
            catalog.append(key)
        pairs=int(r['pairs'])
        bins[by_sample[sid]].extend([taxa[key],pairs])
        if r['rank']=='species' and pairs:
            counts[sid.split('/')[0],key][sid]+=pairs
    denom={d:(sum(int(r['pairs'])>0 for r in ss),sum(int(r['pairs']) for r in ss))
           for d,ss in by_cohort.items()}
    with (folder/'species_recurrence.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow('dataset taxid name group positive_samples nonempty_samples pairs input_pairs ppm top_sample top_sample_pairs top_sample_fraction'.split())
        for (dataset,(tid,name,group,_)),counter in sorted(counts.items()):
            n,total=denom[dataset]
            pairs=sum(counter.values())
            top,peak=sorted(counter.items(),key=lambda kv:(-kv[1],kv[0]))[0]
            w.writerow([dataset,tid,name,group,len(counter),n,pairs,total,1e6*pairs/total,top,peak,peak/pairs])
    qc=defaultdict(Counter)
    for r in rows(folder,'fastq_qc'):
        c=qc[r['sample'],r['mate']]
        bases=int(r['bases'])
        if bases:
            c['bases']+=bases
            c['q30']+=bases*float(r['q30_fraction'])
            c['gc']+=bases*float(r['gc_fraction'])
    compact={}
    for name,keys in {
        'seal_status':'sample state panel_id',
        'seal_bins':'sample panel_id bin_taxid bin_pairs tested_pairs panel_scope target_nuclear target_other_contig target_nonconcordant target_multilocus low_complexity shared_phocid shared_or_conflicting alternative_taxon weak_or_one_mate hit_limit panel_limited unique_nuclear_fragments nuclear_windows_10kb separated_nuclear_loci_10kb most_populated_window_pairs interpretation',
        'seal_loci':'sample bin_taxid accession contig_class window_10kb_start pairs distinct_fragment_spans',
        'sources':'sample group rank taxid name bin_pairs food_source split_evidence',
        'panel':'sample taxid name status exact_pairs anchor_windows_capped split_min_windows',
        'support':'sample group taxid name evidence pairs',
        'status':'sample state sample_type library_type blast amr',
        'blast':'sample split_taxid name rank reference_scope sampled_pairs supports_bin_rank higher_rank_only different_lineage insufficient hit_limit missing_taxonomy panel_not_evaluable',
        'amr':'sample family assigned_pairs both_mates_pairs pairs_per_million',
        'viral':'sample category species_bin_pairs library_type interpretation',
        'evidence':'sample evidence pairs',
        'truth':'sample rank label method positives called tp fp fn precision recall f1',
    }.items():
        source={'sources':'refine_sources','panel':'target_panel_qc','support':'target_species_evidence','status':'refine_status',
                'blast':'refine_blast_by_bin','amr':'refine_amr','viral':'refine_virus',
                'truth':'truth_metrics'}.get(name,name)
        compact[name]=projection(folder,source,keys)
    compact['fastq']=[dict(sample=s,mate=m,bases=int(c['bases']),q30=round(c['q30']/c['bases'],6),
                           gc=round(c['gc']/c['bases'],6)) for (s,m),c in qc.items() if c['bases']]
    data=dict(samples=samples,catalog=catalog,bins=bins,detail=compact,preview=preview)
    payload=json.dumps(data,ensure_ascii=False,separators=(',',':')).replace('&','\\u0026').replace('<','\\u003c').replace('>','\\u003e')
    template=(HERE/'explorer.html').read_text(encoding='utf-8')
    (folder/'explore.html').write_text(template.replace('/*__APP__*/',(HERE/'explorer.js').read_text(encoding='utf-8')).replace('__DATA__',payload),encoding='utf-8')
