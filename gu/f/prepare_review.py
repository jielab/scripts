#!/usr/bin/env python3
"""Read-only GU evidence review. Standard-library Python; never edits GU inputs.

The legacy support table is labelled ANY-LOCUS OVERLAP throughout. Optional
per-copy comparisons use the actual candidate interval in evidence_haplotypes.tsv.
IBDmix is unphased: a match in the same individual does not establish phase.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from typing import Iterable

csv.field_size_limit(100_000_000)
VERSION = '2026-09-05.2'
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dual_lead
MISSING = {'', 'NA', 'NaN', 'nan', 'None', 'null'}


def number(v, default=None):
    if v is None or str(v).strip() in MISSING:
        return default
    try:
        x = float(v)
        return x if math.isfinite(x) else default
    except (ValueError, TypeError):
        return default


def integer(v, default=None):
    x = number(v)
    if x is None:
        return default
    if x != int(x):
        raise ValueError(f'Expected integer, got {v!r}')
    return int(x)


def truth(v):
    return str(v).strip().lower() in {'1', 'true', 't', 'yes', 'pass'}


def chrom(v):
    c = str(v).removeprefix('chr').removeprefix('CHR')
    return 'X' if c == '23' else c


def key(r):
    return (r['dataset_id'], r['genome_build'], r['locus_id'])


def read_tsv(path: Path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'Missing header: {path}')
        for line, r in enumerate(reader, 2):
            if None in r:
                raise ValueError(f'Unexpected extra TSV fields at {path}:{line}')
            yield r


def atomic_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.'+path.name+'.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='') as f:
            f.write(text)
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def write_tsv(path: Path, rows: list[dict], fields=None):
    import io
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    if not fields:
        fields = ['status']
    s = io.StringIO(newline='')
    w = csv.DictWriter(s, fieldnames=fields, delimiter='\t', extrasaction='raise', lineterminator='\n')
    w.writeheader()
    for row in rows:
        w.writerow({k: ('' if v is None else v) for k, v in row.items()})
    atomic_text(path, s.getvalue())


def union(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    out = []
    for st, en in sorted(set(intervals)):
        if st < 0 or en <= st:
            raise ValueError(f'Invalid half-open interval: {st}, {en}')
        if out and st <= out[-1][1]:
            out[-1] = (out[-1][0], max(en, out[-1][1]))
        else:
            out.append((st, en))
    return out


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def contains_pos1(interval, pos):
    """VCF position is 1-based; GU candidate interval is 0-based half-open."""
    return pos is not None and interval[0] <= pos - 1 < interval[1]


def coverage(candidate, intervals):
    if candidate[0] < 0 or candidate[1] <= candidate[0]:
        raise ValueError(f'Invalid candidate interval {candidate}')
    ints = list(intervals)
    clipped = [(max(candidate[0], s), min(candidate[1], e)) for s, e in ints if overlap(candidate, (s, e))]
    reduced = union(clipped)
    bp = sum(e-s for s, e in reduced)
    size = candidate[1]-candidate[0]
    return dict(overlap_bp=bp, union_fraction=bp/size,
                best_single_fraction=max((overlap(candidate, x)/size for x in ints), default=0),
                n_covered_blocks=len(reduced),
                longest_block_fraction=max(((e-s)/size for s, e in reduced), default=0))


def wilson(n, total, z=1.959963984540054):
    if total == 0:
        return None, None
    if not 0 <= n <= total:
        raise ValueError(f'Invalid binomial counts {n}/{total}')
    p=n/total; den=1+z*z/total
    mid=(p+z*z/(2*total))/den
    rad=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return max(0.0, mid-rad), min(1.0, mid+rad)


def contingency(n, p, b, both):
    vals=(both, p-both, b-both, n-p-b+both)
    if min(vals) < 0:
        raise ValueError(f'Inconsistent support counts N={n}, P={p}, B={b}, both={both}')
    lo, hi = wilson(both, p)
    return dict(both=vals[0], phyml_only=vals[1], ibdmix_only=vals[2], neither=vals[3],
                ibdmix_given_phyml=both/p if p else None,
                conditional_ci_low=lo, conditional_ci_high=hi,
                phyml_given_ibdmix=both/b if b else None,
                jaccard=both/(p+b-both) if p+b-both else None)


def read_authoritative(root):
    rows=[]; files=[]
    for path in sorted((root/'phyml').glob('**/final/evidence_loci.tsv')):
        parts=path.relative_to(root).parts
        # normalize/phyml/<dataset>/<scope>/final/evidence_loci.tsv
        if len(parts)<5:
            raise ValueError(f'Cannot resolve dataset from {path}')
        for r in read_tsv(path):
            r['dataset_id']=r.get('dataset_id') or parts[1]
            r['_evidence_file']=str(path)
            if not r.get('genome_build'):
                raise ValueError(f'Genome build missing in {path}; refusing to infer it')
            rows.append(r)
        files.append(path)
    seen=set()
    for r in rows:
        region_key=dual_lead.coord_key(r)
        if region_key in seen:
            raise ValueError(f'Duplicate authoritative region key {region_key}')
        seen.add(region_key)
    if not rows:
        raise ValueError('No phyml/<dataset>/<scope>/final/evidence_loci.tsv found')
    return rows, files


def anchor_assessment(r):
    if not truth(r.get('named_anchor_found')):
        return 'not_evaluable_anchor_missing'
    if r.get('introgression_call') == 'not_evaluable':
        return 'not_evaluable_phyml'
    cs,ce,pos=(integer(r.get(x)) for x in ['candidate_start','candidate_end','named_anchor_pos'])
    if cs is None or ce is None:
        return 'no_supported_candidate_for_anchor'
    if not truth(r.get('candidate_tree_pass')):
        return 'tree_not_supported_for_anchor'
    if not contains_pos1((cs,ce),pos):
        return 'outside_candidate_envelope_allele_origin_not_established'
    status=r.get('named_anchor_support_status','')
    if status in {'named_anchor_not_phase_assigned','lineage_anchor_uninformative_or_representation_mismatch'}:
        return 'not_evaluable_anchor_phase_or_representation'
    if number(r.get('named_anchor_candidate_match')) == 0:
        return 'candidate_does_not_match_lineage_anchor'
    if truth(r.get('named_anchor_candidate_match')):
        return 'allele_linked_candidate_inside_envelope_not_final_origin_proof'
    return 'not_evaluable_anchor_allele'


def prepare_summary(authoritative, legacy_rows, support_rows):
    legacy={key(r):r for r in legacy_rows}
    idx={(key(r),r['method'],r['population']):r for r in support_rows}
    overview=[]
    label_counts=__import__('collections').Counter(key(r) for r in authoritative)
    for r in authoritative:
        k=key(r); ambiguous=label_counts[k]>1
        old={} if ambiguous else legacy.get(k,{})
        ibd={} if ambiguous else idx.get((k,'ibdmix','ALL'),{})
        n=integer(ibd.get('n_tested'),integer(old.get('n_tested'),0))
        p=integer(ibd.get('n_phyml_carriers'),integer(old.get('n_phyml_carriers'),0))
        available=truth(ibd.get('method_available')); eligible=truth(ibd.get('evidence_eligible'))
        b=integer(ibd.get('n_carriers'),0) if available else None
        a=integer(ibd.get('n_phyml_supported'),0) if available else None
        c=contingency(n,p,b,a) if available else {x:None for x in contingency(0,0,0,0)}
        call=r.get('introgression_call','not_evaluable')
        if call == 'not_evaluable':
            state='phyml_not_evaluable'
        elif truth(r.get('candidate_tree_pass')):
            state='tree_plus_broad_ibdmix_overlap' if eligible and a else 'tree_candidate_needs_tract_review'
        elif available and eligible and b:
            state='discordant_phyml_ibdmix'
        elif available and eligible and b == 0:
            state='no_call_under_current_filters_not_proven_absent'
        else:
            state='phyml_no_candidate_ibdmix_not_evaluable'
        cs,ce=integer(r.get('candidate_start')),integer(r.get('candidate_end'))
        pos=integer(r.get('named_anchor_pos'))
        row=dict(legacy_join_status='ambiguous_legacy_label_not_joined' if ambiguous else 'unique_legacy_label',dataset_id=k[0],genome_build=k[1],locus_id=k[2],chr=chrom(r['chr']),
                 core_start=integer(r.get('core_start')),core_end=integer(r.get('core_end')),
                 candidate_start=cs,candidate_end=ce,
                 candidate_envelope_bp=(ce-cs if cs is not None and ce is not None else None),
                 candidate_span_definition='min_start_to_max_end_across_candidate_types_not_shared_contiguous_tract',
                 evidence_lineage=r.get('evidence_lineage') or None,
                 ibdmix_comparison_lineage=ibd.get('source_class'),
                 phyml_call=call,review_state=state,
                 diagnostic_sites=integer(r.get('n_lineage_diagnostic_sites')),
                 sequence_candidate_types=integer(r.get('n_candidate_haplotypes')),
                 sequence_candidate_copies=integer(r.get('n_candidate_copies')),
                 tree_bootstrap=number(r.get('candidate_tree_bootstrap')),
                 candidate_purity=number(r.get('candidate_tree_candidate_purity')),
                 candidate_sensitivity=number(r.get('candidate_tree_candidate_sensitivity')),
                 tree_pass=int(truth(r.get('candidate_tree_pass'))),
                 n_tested=None if ambiguous else n,phyml_carriers=None if ambiguous else p,ibdmix_carriers=b,
                 ibdmix_available=int(available),ibdmix_eligible=int(eligible),
                 ibdmix_availability_note=ibd.get('availability_note'),
                 anchor_pos_1based=pos,
                 anchor_consensus_allele_index=integer(r.get('named_anchor_lineage_consensus_allele_index')),
                 anchor_inside_candidate_envelope=(int(contains_pos1((cs,ce),pos)) if cs is not None and ce is not None and pos is not None else None),
                 anchor_assessment=anchor_assessment(r),
                 upstream_anchor_call=r.get('anchor_linked_introgression_call'),
                 locus_wide_ibdmix_max_lod=number(old.get('ibdmix_max_lod')),
                 raw_best_lineage_descriptive_only=old.get('best_lineage'),
                 raw_sequence_identity_descriptive_only=number(old.get('prop_match')),
                 upstream_reason=r.get('evidence_reason'),
                 support_definition='legacy_same_individual_same_lineage_any_locus_overlap',**c)
        overview.append(row)
    def order(r):
        return (r['dataset_id'],r['genome_build'],99 if r['chr']=='X' else integer(r['chr'],98),r['core_start'] or 0)
    overview.sort(key=order)
    population=[]
    for r in support_rows:
        if r['method']!='ibdmix' or label_counts[key(r)]!=1:
            continue
        n,p,b,a=(integer(r.get(x),0) for x in ['n_tested','n_phyml_carriers','n_carriers','n_phyml_supported'])
        available=truth(r.get('method_available'))
        c=contingency(n,p,b,a) if available else {x:None for x in contingency(0,0,0,0)}
        population.append(dict(dataset_id=r['dataset_id'],genome_build=r['genome_build'],locus_id=r['locus_id'],
                               population=r['population'],source_class=r['source_class'],
                               method_available=int(available),evidence_eligible=int(truth(r.get('evidence_eligible'))),
                               n_tested=n,phyml_carriers=p,ibdmix_carriers=b if available else None,**c))
    return overview,population


def candidate_copies(authoritative):
    rows=[]; files=[]; seen=set()
    for r in authoritative:
        path=Path(r['_evidence_file']).with_name('evidence_haplotypes.tsv')
        if not path.exists():
            continue
        files.append(path)
        tips=set(filter(None,r.get('candidate_tree_candidates_in_clade','').split(',')))
        for h in read_tsv(path):
            if h.get('locus_id')!=r['locus_id'] or h.get('diagnostic_lineage')!=r.get('evidence_lineage'):
                continue
            if not truth(h.get('diagnostic_candidate_pass')):
                continue
            st,en=integer(h.get('candidate_start')),integer(h.get('candidate_end'))
            if st is None or en is None or en<=st or st<0:
                raise ValueError(f'Invalid per-haplotype candidate interval at {path}: {h.get("hap_id")}')
            tokens=[x.strip() for x in h.get('copies','').split(';') if x.strip()]
            if not tokens:
                raise ValueError(f'Candidate with no copy IDs: {path} {h.get("hap_id")}')
            for token in tokens:
                if ':' not in token:
                    raise ValueError(f'Invalid sample:haplotype token {token!r} in {path}')
                sample,hap=token.rsplit(':',1)
                ident=(dual_lead.coord_key(r),sample,hap,h['hap_id'],st,en)
                if ident in seen:
                    continue
                seen.add(ident)
                rows.append(dict(locus_key=dual_lead.coord_key(r),dataset_id=r['dataset_id'],genome_build=r['genome_build'],locus_id=r['locus_id'],
                                 chr=chrom(r['chr']),sample_id=sample,haplotype=hap,hap_id=h['hap_id'],
                                 candidate_start=st,candidate_end=en,source_class=r['evidence_lineage'],
                                 candidate_group=('tree_supported_copy' if truth(r.get('candidate_tree_pass')) and h['hap_id'] in tips else 'sequence_candidate_copy'),
                                 input_lead_snp=r.get('_input_lead_override',r.get('named_anchor_id',r['locus_id'])),
                                 upstream_anchor_pos_1based=integer(r.get('named_anchor_pos')),
                                 anchor_pos_1based=integer(r.get('_display_input_pos1',r.get('named_anchor_pos')))))
    return rows,files


def selected_segments(root, database, candidates, offset):
    # One pass through a compressed file. With SQLite, use only indexed locus windows.
    windows=defaultdict(list)
    for r in candidates:
        windows[(r['dataset_id'],r['genome_build'],r['chr'])].append((r['candidate_start'],r['candidate_end']))
    windows={k:union(v) for k,v in windows.items()}
    targets={(r['dataset_id'],r['genome_build'],r['chr'],r['sample_id']) for r in candidates}
    kept={}; nread=0
    def accept(r):
        nonlocal nread
        nread+=1
        if r.get('method')!='ibdmix':
            return
        unit=(r.get('dataset_id'),r.get('genome_build'),chrom(r.get('chr')))
        if (*unit,r.get('sample_id')) not in targets or unit not in windows:
            return
        st,en=integer(r.get('start')),integer(r.get('end'))
        if st is None or en is None or st+offset<0 or en<=st:
            raise ValueError('Invalid IBDmix interval; refusing silent row removal')
        st+=offset; en+=offset
        if not any(overlap((st,en),q) for q in windows[unit]):
            return
        out={x:r.get(x) for x in ['dataset_id','genome_build','sample_id','method','source','source_class','locus_id','raw_file']}
        out.update(chr=unit[2],start=st,end=en,score=number(r.get('score')))
        ident=(*unit,out['sample_id'],out['source'],st,en,out['score'])
        kept[ident]=out
    if database:
        uri=Path(database).resolve().as_uri()+'?mode=ro'
        with sqlite3.connect(uri,uri=True) as con:
            con.row_factory=sqlite3.Row
            cols={r[1] for r in con.execute('PRAGMA table_info(segments)')}
            need={'dataset_id','genome_build','chr','start','end','sample_id','method','source_class'}
            if not need<=cols:
                raise ValueError('Database segments schema is incompatible')
            for (d,b,c),ints in windows.items():
                for st,en in ints:
                    for row in con.execute("SELECT * FROM segments WHERE dataset_id=? AND genome_build=? AND chr=? AND method='ibdmix' AND end>? AND start<?",(d,b,c,st-offset,en-offset)):
                        accept(dict(row))
        source=str(Path(database).resolve())
    else:
        path=root/'summary'/'segments.tsv.gz'
        if not path.exists():
            return None,{'status':'not_computed_segments_missing'}
        for r in read_tsv(path):
            accept(r)
        source=str(path)
    return list(kept.values()),dict(status='computed',source=source,rows_read=nread,selected_calls=len(kept),ibdmix_coordinate_offset=offset)


def exact_support(candidates, segments, runs, populations, threshold):
    bysample=defaultdict(list)
    for s in segments:
        bysample[(s['dataset_id'],s['genome_build'],s['chr'],s['sample_id'])].append(s)
    run={}
    for r in runs:
        if r['method']=='ibdmix':
            k=(r['dataset_id'],r['genome_build'],chrom(r['chr']))
            prev=run.get(k)
            if prev and truth(prev.get('evidence_eligible'))!=truth(r.get('evidence_eligible')):
                raise ValueError(f'Conflicting IBDmix eligibility: {k}')
            run[k]=r
    pop={(r['dataset_id'],r['sample_id']):(r.get('population'),r.get('super_population')) for r in populations}
    out=[]
    for c in candidates:
        unit=(c['dataset_id'],c['genome_build'],c['chr'])
        pool=bysample.get((*unit,c['sample_id']),[])
        interval=(c['candidate_start'],c['candidate_end'])
        same=[s for s in pool if s['source_class']==c['source_class'] and overlap(interval,(s['start'],s['end']))]
        other=[s for s in pool if s['source_class']!=c['source_class'] and overlap(interval,(s['start'],s['end']))]
        metric=coverage(interval,[(s['start'],s['end']) for s in same])
        info=run.get(unit)
        completed=bool(info and info.get('status') == 'complete')
        eligible=bool(completed and truth(info.get('evidence_eligible')))
        if same:
            status='observed_same_individual_overlap' if eligible else 'observed_overlap_eligibility_unconfirmed_or_exploratory'
        else:
            status='not_detected_individual_callability_unknown' if completed else 'not_evaluable_no_completion_metadata'
        known=int(bool(info))
        if not same and not completed:
            metric={x:None for x in metric}
        pp,sp=pop.get((c['dataset_id'],c['sample_id']),(None,None))
        r=dict(c, population=pp,super_population=sp,**metric,
               n_matching_reference_calls=len(same),
               matching_references=','.join(sorted({s['source'] for s in same if s['source']})),
               matched_interval_max_lod=max((s['score'] for s in same if s['score'] is not None),default=None),
               other_lineage_overlap=int(bool(other)),
               other_lineage_max_lod=max((s['score'] for s in other if s['score'] is not None),default=None),
               anchor_inside_this_candidate=int(contains_pos1(interval,c['anchor_pos_1based'])) if c['anchor_pos_1based'] else None,
               ibdmix_covers_anchor_in_same_individual=(int(any(contains_pos1((s['start'],s['end']),c['anchor_pos_1based']) for s in same)) if c['anchor_pos_1based'] else None),
               completion_metadata_present=known,evidence_eligible=int(eligible),
               strict_coverage_threshold=threshold,
               strict_coverage_pass=(int(metric['best_single_fraction']>=threshold) if metric['best_single_fraction'] is not None and eligible else None),
               callability='unknown_individual_callable_bases_not_in_normalized_summary',
               phase_corroboration='individual_only_IBDmix_phase_unknown',
               status=status)
        out.append(r)
    return out


def collapse_exact(rows, thresholds=(.25,.5,.8)):
    grouped=defaultdict(list)
    for r in rows:
        grouped[(*key(r),r['candidate_group'],r.get('locus_key',''))].append(r)
    out=[]
    for k, rr in grouped.items():
        samples={r['sample_id'] for r in rr}
        eligible=[r for r in rr if r['evidence_eligible']==1 and r['best_single_fraction'] is not None]
        for t in thresholds:
            passing={r['sample_id'] for r in eligible if r['best_single_fraction']>=t}
            seen={r['sample_id'] for r in eligible}
            out.append(dict(locus_key=k[4] or None,dataset_id=k[0],genome_build=k[1],locus_id=k[2],candidate_group=k[3],
                            threshold=t,n_candidate_individuals=len(samples),n_candidate_copies=len(rr),
                            n_individuals_with_eligible_comparison=len(seen),n_individuals_passing=len(passing),
                            support_fraction=len(passing)/len(seen) if seen else None,
                            comparison='single_matching_lineage_IBDmix_call_covers_candidate_fraction',
                            phase='not_confirmed'))
    return out


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--normalize',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--database',type=Path,help='Optional existing GU SQLite, read-only')
    ap.add_argument('--summary-only',action='store_true',help='Skip strict IBDmix comparison; dual-lead analysis remains enabled')
    ap.add_argument('--loci',type=Path,help='Optional CHR START END lead file; exact coordinate match, fourth column is display metadata only')
    ap.add_argument('--loci-format',choices=['bed','1based'],default='bed')
    ap.add_argument('--skip-tags',action='store_true',help='Explicitly skip dual-lead sequence analysis')
    ap.add_argument('--min-match-sites',type=int,default=10)
    ap.add_argument('--min-match-callable',type=float,default=.8)
    ap.add_argument('--min-tag-call-rate',type=float,default=.95)
    ap.add_argument('--min-tag-target-copies',type=int,default=2)
    ap.add_argument('--strong-tag-r2',type=float,default=.8,help='In-sample proxy quality label; does not change the regional introgression call')
    ap.add_argument('--top-tags',type=int,default=20)
    ap.add_argument('--tag-population-level',choices=['none','super','all'],default='super',help='Reselect tag within superpopulations or all groups; fixed-tag metrics are always exported for all known populations')
    ap.add_argument('--coverage-threshold',type=float,default=.8)
    ap.add_argument('--ibdmix-coordinate-offset',type=int,choices=[-1,0,1],default=0,
                    help='Explicit audit-only shift of BOTH IBDmix bounds; default preserves normalized coordinates')
    args=ap.parse_args(argv)
    if not 0<args.coverage_threshold<=1:
        ap.error('--coverage-threshold must be in (0,1]')
    if args.min_match_sites<1 or args.min_tag_target_copies<1 or args.top_tags<1:
        ap.error('Site, copy and top-tag counts must be positive')
    if not all(0 < v <= 1 for v in [args.min_match_callable,args.min_tag_call_rate,args.strong_tag_r2]):
        ap.error('Match/tag fractions and r2 thresholds must be in (0,1]')
    root=args.normalize.resolve(); out=args.output.resolve()
    if not root.is_dir():
        ap.error('normalize directory does not exist')
    if out==root or root in out.parents:
        ap.error('output must be outside the normalized input tree')
    if out in root.parents:
        ap.error('output must not be an ancestor of normalized input')
    auth,files=read_authoritative(root)
    if args.loci:
        auth=dual_lead.apply_loci_file(auth,args.loci,args.loci_format)
        files.append(args.loci.resolve())
    def load(name,required=True):
        path=root/'summary'/(name+'.tsv.gz')
        if not path.exists():
            if required: raise ValueError(f'Missing required summary: {path}')
            return []
        files.append(path)
        return list(read_tsv(path))
    legacy=load('locus_evidence');support=load('locus_method_support');runs=load('method_runs',False)
    overview,population=prepare_summary(auth,legacy,support)
    pops=load('sample_populations',False)
    if args.skip_tags:
        tags={k:[] for k in ['summary','comparisons','ranked','population','by_population','haplotypes','gwas','files']}
        for r in auth:
            tags['summary'].append(dict(locus_key=dual_lead.coord_key(r),dataset_id=r['dataset_id'],genome_build=r['genome_build'],locus_id=r['locus_id'],input_lead_snp=r.get('_input_lead_override',r.get('named_anchor_id',r['locus_id'])),tag_status='skipped_by_user'))
    else:
        tags=dual_lead.analyse_all(auth,pops,args)
        files+=tags.pop('files')
    ti={t['locus_key']:t for t in tags['summary']}
    for r,a in zip(overview,sorted(auth,key=lambda z:(z['dataset_id'],z['genome_build'],99 if chrom(z['chr'])=='X' else integer(z['chr'],98),integer(z.get('core_start'),0)))):
        r['locus_key']=dual_lead.coord_key(a)
        r['input_lead_snp']=a.get('_input_lead_override',a.get('named_anchor_id',a.get('name',a['locus_id'])))
        t=ti.get(r['locus_key'],{})
        a['_display_input_pos1']=(t.get('input_pos_1based') if '_input_lead_override' in a and a['_input_lead_override']!=a.get('named_anchor_id',a['locus_id']) else integer(a.get('named_anchor_pos')))
        for field in ['best_haplotype_id','best_haplotype_tier','best_haplotype_lineage','best_tag_snp','best_tag_pos_1based','best_tag_allele','best_tag_r2','best_tag_quality','input_tag_r2','input_tag_allele','input_vs_best_snp_r2','same_input_and_best','n_equivalent_best_tags','tag_status','tag_detail']:
            r[field]=t.get(field)
        if '_input_lead_override' in a and a['_input_lead_override']!=a.get('named_anchor_id',a['locus_id']):
            r['upstream_anchor_snp']=a.get('named_anchor_id',a['locus_id'])
            r['anchor_assessment']='not_evaluable_input_changed'
            r['anchor_pos_1based']=t.get('input_pos_1based')
            r['anchor_consensus_allele_index']=None
            r['upstream_anchor_call']=None
    coordinates={key(r):r['locus_key'] for r in overview}
    for r in population:r['locus_key']=coordinates.get(key(r))
    copies=[];tracks=[];strict=[];status={'status':'not_computed_summary_only'}
    if not args.summary_only:
        candidates,cfiles=candidate_copies(auth);files+=cfiles
        if not candidates:
            status={'status':'not_computed_per_haplotype_candidates_missing_or_empty'}
        else:
            tracks,status=selected_segments(root,args.database,candidates,args.ibdmix_coordinate_offset)
            if tracks is not None:
                copies=exact_support(candidates,tracks,runs,pops,args.coverage_threshold)
                strict=collapse_exact(copies,tuple(sorted({.25,.5,.8,args.coverage_threshold})))
            else: tracks=[]
    manifest={'version':VERSION,'created_utc':datetime.now(timezone.utc).isoformat(),
              'normalize_root':str(root),'output_root':str(out),'input_files':[],
              'n_loci':len(overview),'strict_comparison':status,
              'dual_lead':{'status_counts':dict(__import__('collections').Counter(t.get('tag_status') for t in tags['summary'])),
                          'analysis_unit':'dataset + genome build + CHR:START-END (0-based half-open)',
                          'input_lead_role':'annotation/comparator only; never a target or rank constraint',
                          'parameters':{k:getattr(args,k) for k in ['min_match_sites','min_match_callable','min_tag_call_rate','min_tag_target_copies','strong_tag_r2','top_tags','tag_population_level']},
                          'selection_rule':'evidence tier, diagnostic match proportion, contiguous diagnostic matches, raw identity, informative sites, copies; tag ranked by copy-level r2 then F1/call rate; no GWAS',
                          'scope_limit':'Best among stored biallelic SNPs with phased A/C/G/T calls. Exact input indels use anchor_copies.tsv. Not a genome-wide or all-variant VCF rerun.'},
              'coordinate_note':'Candidate start/end are GU 0-based half-open; anchor is 1-based. IBDmix input coordinates are preserved by default. Verify native coordinate provenance before endpoint-specific SNP claims.',
              'inference_limits':['No introgression probability is calculated.','Legacy overlap is not tract coverage.',
                                  'IBDmix phase and individual callable-base denominator are unknown.',
                                  'Multiple archaic references are not independent methods.',
                                  'No VCF, phylogenetic tree, IBDmix, TRACE or AS3 caller is rerun.', 'Best tag and all proxy quality metrics are in-sample, not external validation or GWAS association.', 'Population tag reselection keeps the same global haplotype target.']}
    for path in sorted(set(files)):
        data=path.read_bytes()
        manifest['input_files'].append({'path':str(path.relative_to(root)) if path.is_relative_to(root) else str(path),'size':len(data),
                                       'sha256':hashlib.sha256(data).hexdigest(),
                                       'git_blob_sha1':hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()})
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/'loci_review.tsv',overview)
    write_tsv(out/'population_concordance.tsv',population)
    write_tsv(out/'method_runs.tsv',runs)
    write_tsv(out/'candidate_copy_support.tsv',copies)
    write_tsv(out/'strict_support_summary.tsv',strict)
    write_tsv(out/'selected_ibdmix_calls.tsv',tracks)
    for label,name in [('summary','locus_tag_summary'),('comparisons','lead_comparison'),('ranked','haplotype_tag_candidates'),('population','tag_population_metrics'),('by_population','best_tag_by_population'),('haplotypes','best_haplotypes'),('gwas','gwas_lookup')]:
        write_tsv(out/(name+'.tsv'),tags[label])
    payload={'dual_lead':tags,'overview':overview,'population':population,'methods':runs,'copies':copies,
             'strict':strict,'tracks':tracks,'manifest':manifest}
    atomic_text(out/'review_data.json',json.dumps(payload,ensure_ascii=False,allow_nan=False,indent=2))
    atomic_text(out/'manifest.json',json.dumps(manifest,ensure_ascii=False,indent=2))
    template=Path(__file__).with_name('review_template.html').read_text(encoding='utf-8')
    safe_json=json.dumps(payload,ensure_ascii=False,allow_nan=False).replace('<','\\u003c')
    atomic_text(out/'review.html',template.replace('/*__GU_DATA__*/',safe_json))
    print(json.dumps({'loci':len(overview),'candidate_copies_compared':len(copies),
                      'strict_status':status['status'],'output':str(out)},ensure_ascii=False,indent=2))
    return 0

if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (OSError,ValueError,KeyError,csv.Error,sqlite3.Error) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        raise SystemExit(2)
