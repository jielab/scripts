#!/usr/bin/env python3
"""Describe within-population versus pooled ALT diversity at one SNV site.

ALL is a pooled label, never a second ancestry group. A reported zero differs
from an absent field/population. Patterns concern the ALT identities actually
reported by the source, not all possible alleles, REF occupancy or selection.
The pure API consumes one site's primary rows while the scanner streams them.
"""
from __future__ import annotations

import itertools
import math

PATTERNS=('within_population_multi_ALT','between_population_only','unresolved','not_multi_ALT')
LAYER_FIELDS=('pattern','union_alt_bases','n_union_alts','max_alts_per_population',
    'multi_alt_populations','n_multi_alt_populations','n_complete_populations','n_incomplete_populations')
FIELDS=['chrom','pos','ref','reported_alt_bases','candidate','structure_status',
    'ancestry_populations','required_populations','missing_populations',
    'n_populations_reported','n_populations_required','n_populations_missing',
    'pooled_available','pooled_observed_alt_bases','pooled_qualified_alt_bases',
    'pooled_n_observed_alts','pooled_n_qualified_alts','pooled_observed_complete','pooled_qualified_complete',
    'threshold_kind','threshold','qualified_threshold_comparable',
    *[layer+'_'+field for layer in ('observed','qualified') for field in LAYER_FIELDS],
    'af_rank_status','af_rank_comparable_pairs','af_rank_reversed_pairs',
    'af_rank_comparison_populations','af_rank_tied_populations']
COUNT_FIELDS=['sites','candidate_sites','pooled_available_sites','pooled_observed_multi_alt_sites','pooled_qualified_multi_alt_sites',
    'sites_with_at_least_two_ancestry_groups','pooled_only_sites','sites_with_missing_populations','invalid_ref_sites',
    *[layer+'_'+pattern+'_sites' for layer in ('observed','qualified') for pattern in PATTERNS],
    'qualified_threshold_mismatch_sites','af_rank_comparable_sites','af_rank_reversal_sites']


def _number(value,count=False):
    if value is None:return None
    try:value=float(value)
    except (ValueError,TypeError):return None
    if not math.isfinite(value) or value<0:return None
    if count:return int(value) if value.is_integer() else None
    return value if value<=1 else None


def _evidence(row,universe):
    alts=list(row['alts']);ac=list(row['alt_acs']);af=list(row['alt_afs'])
    if len(alts)!=len(ac) or len(alts)!=len(af):raise ValueError('ALT evidence arrays must align')
    if len(set(alts))!=len(alts) or any(b not in 'ACGT' or len(b)!=1 or b==row['ref'] for b in alts):
        raise ValueError('Population patterns require distinct SNV ALTs different from REF')
    count_mode=row['count_mode']!='AF_only'
    kind=str(row['threshold_kind']).upper()
    if kind not in {'AC','AF'}:raise ValueError('Primary threshold_kind must be AC or AF')
    threshold=float(row['threshold'])
    if not math.isfinite(threshold) or threshold<0:raise ValueError('Invalid primary threshold')
    counts={b:_number(v,True) for b,v in zip(alts,ac)}
    freqs={b:_number(v) for b,v in zip(alts,af)}
    observed={b:counts.get(b) if count_mode else freqs.get(b) for b in universe}
    qualified={}
    for b in universe:
        if kind=='AC':value=counts.get(b)
        elif not count_mode:value=freqs.get(b)
        else:
            count=counts.get(b)
            value=None if count is None else 0. if count==0 else freqs.get(b)
        qualified[b]=None if value is None else bool(value>0 and value>=threshold)
    return dict(observed={b:None if v is None else v>0 for b,v in observed.items()},
        qualified=qualified,af=freqs,rule=(kind,threshold))


def _positive(evidence,layer):return {b for b,present in evidence[layer].items() if present is True}


def _layer(evidence,populations,expected,layer,comparable=True):
    positives={p:_positive(evidence[p],layer) for p in populations}
    union=set().union(*positives.values()) if positives else set()
    multi=[p for p in populations if len(positives[p])>=2]
    complete=[p for p in populations if evidence[p][layer] and all(v is not None for v in evidence[p][layer].values())]
    all_complete=len(complete)==len(expected)
    if multi:pattern='within_population_multi_ALT'
    elif len(expected)<2 or not all_complete or not comparable:pattern='unresolved'
    elif len(union)>=2:pattern='between_population_only'
    else:pattern='not_multi_ALT'
    return dict(pattern=pattern,union_alt_bases=''.join(sorted(union)),n_union_alts=len(union),
        max_alts_per_population=max(map(len,positives.values())) if positives else None,
        multi_alt_populations=','.join(multi),n_multi_alt_populations=len(multi),
        n_complete_populations=len(complete),n_incomplete_populations=len(expected)-len(complete))


def _ranks(evidence,populations,universe):
    comparable=0;reversed_pairs=0;used=set();tied=set()
    for a,b in itertools.combinations(sorted(universe),2):
        ranks={}
        for p in populations:
            e=evidence[p]
            if e['observed'].get(a) is not True or e['observed'].get(b) is not True:continue
            av,bv=e['af'].get(a),e['af'].get(b)
            if av is None or bv is None or av<=0 or bv<=0:continue
            if av==bv:tied.add(p);continue
            ranks[p]=1 if av>bv else -1
        if len(ranks)>=2:
            comparable+=1;used.update(ranks)
            reversed_pairs+=int(len(set(ranks.values()))>1)
    return dict(af_rank_status='reversed' if reversed_pairs else 'consistent' if comparable else 'unresolved',
        af_rank_comparable_pairs=comparable,af_rank_reversed_pairs=reversed_pairs,
        af_rank_comparison_populations=','.join(sorted(used)),af_rank_tied_populations=','.join(sorted(tied)))


def classify_site(rows,required_populations=None):
    """Return one compact result for typed primary rows at the same chr/pos/REF.

    Required populations may identify absent requested rows. Every supplied
    ancestry group is also considered. Completeness is evaluated over the ALT
    set reported in any supplied row, including ALL. An ALT omitted in one
    group's row is unknown there, not zero. Nonreported possible bases outside
    that source ALT set are outside this descriptive comparison.
    """
    rows=list(rows)
    if not rows:raise ValueError('classify_site requires at least one population row')
    first=rows[0];key=(str(first['chrom']),int(first['pos']),str(first['ref']))
    if key[1]<1 or key[2] not in 'ACGT' or len(key[2])!=1:raise ValueError('Invalid SNV coordinate/REF')
    invalid_ref=any('REF_conflict' in str(r.get('flags','')).split(';') for r in rows)
    by_pop={}
    for r in rows:
        if (str(r['chrom']),int(r['pos']),str(r['ref']))!=key:raise ValueError('Population rows must describe one chr/pos/REF')
        p=str(r['population']).upper()
        if not p or p in by_pop:raise ValueError('Population rows must have distinct population labels')
        if invalid_ref:
            # The scanner retains conflicting-REF anchors with unknown evidence.
            # One conflicting record's ALT may equal the retained REF identity;
            # preserve the raw union below, but never interpret those counts.
            safe_alts=[b for b in r['alts'] if b!=key[2]]
            r=dict(r,alts=safe_alts,alt_acs=[None]*len(safe_alts),alt_afs=[None]*len(safe_alts))
        by_pop[p]=r
    populations=sorted(p for p in by_pop if p!='ALL')
    supplied_required=required_populations or []
    if isinstance(supplied_required,str):supplied_required=supplied_required.split(',')
    expected=sorted(set(populations)|{str(p).strip().upper() for p in supplied_required if str(p).strip().upper() not in {'','ALL'}})
    missing=sorted(set(expected)-set(populations))
    universe=set().union(*(set(r['alts']) for r in rows))
    evidence={p:_evidence(r,universe) for p,r in by_pop.items()}
    rules={evidence[p]['rule'] for p in populations}
    comparable=len(rules)==1
    observed=_layer(evidence,populations,expected,'observed')
    qualified=_layer(evidence,populations,expected,'qualified',comparable)
    pooled=evidence.get('ALL')
    pooled_obs=_positive(pooled,'observed') if pooled else set()
    pooled_qual=_positive(pooled,'qualified') if pooled else set()
    all_rules={e['rule'] for e in evidence.values()}
    same_rule=next(iter(all_rules)) if len(all_rules)==1 else None
    structure=('invalid_REF' if invalid_ref else 'pooled_only' if not expected else 'missing_population_rows' if missing else
               'single_ancestry_group' if len(populations)==1 else 'ancestry_groups_available')
    result=dict(chrom=key[0],pos=key[1],ref=key[2],reported_alt_bases=''.join(sorted(universe)),
        candidate=bool(len(pooled_qual)>=2 or qualified['n_union_alts']>=2),structure_status=structure,
        ancestry_populations=','.join(populations),required_populations=','.join(expected),missing_populations=','.join(missing),
        n_populations_reported=len(populations),n_populations_required=len(expected),n_populations_missing=len(missing),
        pooled_available=bool(pooled),pooled_observed_alt_bases=''.join(sorted(pooled_obs)) if pooled else None,
        pooled_qualified_alt_bases=''.join(sorted(pooled_qual)) if pooled else None,
        pooled_n_observed_alts=len(pooled_obs) if pooled else None,pooled_n_qualified_alts=len(pooled_qual) if pooled else None,
        pooled_observed_complete=bool(pooled['observed']) and all(v is not None for v in pooled['observed'].values()) if pooled else None,
        pooled_qualified_complete=bool(pooled['qualified']) and all(v is not None for v in pooled['qualified'].values()) if pooled else None,
        threshold_kind=','.join(sorted({rule[0] for rule in all_rules})),threshold=same_rule[1] if same_rule else None,
        qualified_threshold_comparable=comparable)
    for layer,values in [('observed',observed),('qualified',qualified)]:result.update({layer+'_'+k:v for k,v in values.items()})
    result.update(_ranks(evidence,populations,universe))
    return result


def add_counts(counter,row):
    """Accumulate all input sites, including noncandidate/unresolved sites."""
    counter['sites']+=1;counter['candidate_sites']+=int(row['candidate'])
    counter['pooled_available_sites']+=int(row['pooled_available'])
    for layer in ('observed','qualified'):
        counter['pooled_'+layer+'_multi_alt_sites']+=int((row['pooled_n_'+layer+'_alts'] or 0)>=2)
        counter[layer+'_'+row[layer+'_pattern']+'_sites']+=1
    counter['sites_with_at_least_two_ancestry_groups']+=int(row['n_populations_reported']>=2)
    counter['pooled_only_sites']+=int(row['structure_status']=='pooled_only')
    counter['sites_with_missing_populations']+=int(row['n_populations_missing']>0)
    counter['invalid_ref_sites']+=int(row['structure_status']=='invalid_REF')
    counter['qualified_threshold_mismatch_sites']+=int(row['n_populations_reported']>=2 and not row['qualified_threshold_comparable'])
    counter['af_rank_comparable_sites']+=int(row['af_rank_comparable_pairs']>0)
    counter['af_rank_reversal_sites']+=int(row['af_rank_reversed_pairs']>0)
