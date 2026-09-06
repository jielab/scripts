#!/usr/bin/env python3
"""Region-first post-processing: original lead vs haplotype-derived tag.

Uses only saved, phased GU sequences and copy IDs. No GWAS input, risk allele,
lead-SNP LD window, environment changes, or rerunning of scientific callers.
The binary SNP--haplotype r^2 is calculated on chromosome copies, not people.
"""
from __future__ import annotations
import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from collections.abc import Mapping
from pathlib import Path
from comm import resolve_tsv_path

BASES = set('ACGT')
VERSION = '2026-09-06.1'
MISSING = {'', '.', 'NA', 'nan', 'None', 'null'}


def num(v, default=None):
    try:
        x = float(v)
        return x if math.isfinite(x) else default
    except (ValueError, TypeError):
        return default


def integer(v, default=None):
    x = num(v)
    if x is None:
        return default
    if x != int(x):
        raise ValueError(f'Nonintegral value: {v!r}')
    return int(x)


def truth(v):
    return str(v).strip().upper() in {'1', 'TRUE', 'T', 'YES', 'PASS'}


def chrom(v):
    x = re.sub(r'^chr', '', str(v), flags=re.I)
    return 'X' if x == '23' else x


def lineage(v):
    x = str(v).lower()
    if any(t in x for t in ('neander', 'altai', 'vindija', 'chagyr')):
        return 'Neanderthal'
    if 'denis' in x:
        return 'Denisovan'
    return str(v)


def rows(path):
    if path is None:
        return []
    path = resolve_tsv_path(Path(path))
    if not path.is_file():
        return []
    op = gzip.open if path.suffix == '.gz' else open
    with op(path, 'rt', encoding='utf-8-sig', newline='') as f:
        out = list(csv.DictReader(f, delimiter='\t'))
    if any(None in r for r in out):
        raise ValueError(f'Malformed TSV fields: {path}')
    return out


class SiteAlleles(Mapping):
    """Lazy shared sequence-backed map: avoids copies x SNPs dictionaries."""
    def __init__(self, copy_sequences, index, allowed):
        self.sequences=copy_sequences
        self.index=index
        self.allowed=set(allowed)
    def get(self, key, default=None):
        seq=self.sequences.get(key)
        if seq is None:return default
        base=seq[self.index]
        return base if base in self.allowed else default
    def __getitem__(self,key):
        value=self.get(key)
        if value is None:raise KeyError(key)
        return value
    def __iter__(self):
        return (cp for cp in self.sequences if self.get(cp) is not None)
    def __len__(self):
        return sum(self.get(cp) is not None for cp in self.sequences)
    def __bool__(self):
        return any(self.get(cp) is not None for cp in self.sequences)


def coord_key(r):
    st, en = integer(r.get('core_start')), integer(r.get('core_end'))
    if st is None or en is None or st < 0 or en <= st:
        raise ValueError('A valid core_start/core_end is required for a region-first locus key')
    return f"{r['dataset_id']}|{r['genome_build']}|{chrom(r['chr'])}:{st}-{en}"


def aliases(v):
    return {x.strip() for x in re.split('[;,]', str(v or '')) if x.strip() not in MISSING}


def copy_ids(h):
    tokens = [x.strip() for x in str(h.get('copies', '')).split(';') if x.strip()]
    if len(tokens) != len(set(tokens)) or any(':' not in t for t in tokens):
        raise ValueError(f'Duplicate or invalid sample:haplotype copy ID: {h.get("hap_id")}')
    if not tokens:
        raise ValueError(f'Copy IDs missing for haplotype {h.get("hap_id")}')
    return set(tokens)


def inside(st, en, pos1):
    return pos1 is not None and st <= pos1 - 1 < en


def binary_metrics(tp, fp, fn, tn, total_target=None, total=None):
    """Target=haplotype membership; predicted=specified allele present.

    PPV/sensitivity here describe a tag for this computational target, NOT
    accuracy of introgression detection and NOT tree purity.
    """
    a, b, c, d = tp, fp, fn, tn
    if min(a, b, c, d) < 0:
        raise ValueError('Negative contingency count')
    n = a+b+c+d
    denom = (a+b)*(c+d)*(a+c)*(b+d)
    signed = (a*d-b*c)/math.sqrt(denom) if denom else None
    return dict(n_callable_copies=n, n_target_callable=a+c,
                n_background_callable=b+d, tp=a, fp=b, fn=c, tn=d,
                r=signed, r2=signed*signed if signed is not None else None,
                sensitivity=a/(a+c) if a+c else None,
                specificity=d/(b+d) if b+d else None,
                ppv=a/(a+b) if a+b else None,
                f1=2*a/(2*a+b+c) if 2*a+b+c else None,
                allele_frequency=(a+b)/n if n else None,
                haplotype_frequency=(a+c)/n if n else None,
                call_rate=n/total if total else None,
                target_call_rate=(a+c)/total_target if total_target else None)


def count_variant(alleles, target, cohort, tag_allele):
    a=b=c=d=0
    for cp in cohort:
        value = alleles.get(cp)
        if value is None:
            continue
        ish = cp in target
        ist = value == tag_allele
        if ish and ist: a += 1
        elif ist: b += 1
        elif ish: c += 1
        else: d += 1
    return binary_metrics(a,b,c,d,len(target & cohort),len(cohort))


def qualify(m, min_call_rate=.95, min_target=2):
    if m['n_target_callable'] < min_target:
        return 'too_few_callable_target_copies'
    if m['n_background_callable'] < 2:
        return 'too_few_callable_background_copies'
    if m['call_rate'] is None or m['call_rate'] < min_call_rate:
        return 'low_call_rate'
    if m['target_call_rate'] is None or m['target_call_rate'] < min_call_rate:
        return 'low_target_call_rate'
    if m['r'] is None:
        return 'monomorphic_allele_or_target'
    if m['r'] <= 0:
        return 'allele_not_enriched_in_target'
    return 'eligible'


def rank_tag(r):
    # Neither input lead position, input lead identity nor GWAS appears here.
    return (-round(r['r2'],12), -round(r['f1'],12),
            -round(r['call_rate'],12), -int(r['inside_candidate']),
            -int(r['diagnostic_for_target_lineage']),
            r['pos_1based'], r['ref'], r['alt'], r['tag_allele'], r['variant_key'])


def rank_evidence(r):
    return (-r['evidence_tier'],
            -(r.get('diagnostic_match_prop') or 0),
            -r.get('n_contiguous_diagnostic_match',0),
            -r['prop_match'], -r['n_compared'], -r['n_copies'],
            r['hap_id'], r['lineage'], r['archaic'])


def make_groups(copies, population_rows):
    meta = {}
    for p in population_rows:
        sid = p.get('sample_id',p.get('sample'))
        v = (p.get('population',p.get('pop')),p.get('super_population',p.get('super_pop')))
        if sid in meta and meta[sid] != v:
            raise ValueError(f'Conflicting population metadata for {sid}')
        meta[sid] = v
    groups = {'ALL': set(copies)}
    for cp in copies:
        pop, sup = meta.get(cp.rsplit(':',1)[0],(None,None))
        for prefix, v in [('POP',pop),('SUPER',sup)]:
            if v and v not in MISSING:
                groups.setdefault(prefix+':'+v,set()).add(cp)
    return groups


def resolve_sequence_files(auth):
    f = Path(auth['_evidence_file']).parent
    locus = f.parent/'loci'
    if not (locus/'sites.tsv').is_file():
        locus = locus/auth['locus_id']
    # The region-wide catalogue is preferred; legacy tables may be LD-selected.
    region = resolve_tsv_path(f/'region_haplotypes.unfiltered.tsv.gz')
    hp = region if region.is_file() else f/'haplotypes.tsv'
    return f, locus, hp, hp == region


def load_data(auth):
    final, locus, hp, region_inventory = resolve_sequence_files(auth)
    paths = [hp,locus/'sites.tsv',locus/'archaic.tsv']
    missing = [p.name for p in paths if not p.is_file()]
    if missing:
        return None, {'status':'sequence_files_missing','detail':','.join(missing)}, []
    sites = rows(paths[1]); haps = rows(hp); arch = rows(paths[2])
    haps = [h for h in haps if not h.get('locus_id') or h['locus_id']==auth['locus_id']]
    if not sites or not haps or not arch:
        return None, {'status':'empty_sequence_catalogue'}, paths
    needed = {'chr','pos','id','ref','alt'}
    if not needed <= sites[0].keys():
        raise ValueError('sites.tsv lacks required chr/pos/id/ref/alt columns')
    nsite = len(sites)
    if any(len(str(h.get('seq',''))) != nsite for h in haps+arch):
        raise ValueError('Sequence/site length mismatch; never truncate or infer SNP positions')
    if len({h['hap_id'] for h in haps}) != len(haps):
        raise ValueError('Duplicate hap_id in the same locus catalogue')
    st,en=integer(auth['core_start']),integer(auth['core_end'])
    if region_inventory:
        if any(integer(h.get('core_start')) != st or integer(h.get('core_end')) != en for h in haps):
            raise ValueError('Unfiltered catalogue does not match the full input locus coordinates')
        scope='region_wide_unfiltered_catalogue'
    else:
        lr=[r for r in rows(final/'loci.tsv') if r.get('locus_id')==auth['locus_id']]
        if not lr:
            return None, {'status':'cannot_verify_region_wide_scope','detail':'Legacy haplotypes need loci.tsv core/selected bounds'}, paths
        sel = lr[0]
        s0=integer(sel.get('selected_start'),integer(sel.get('region_start')))
        e0=integer(sel.get('selected_end'),integer(sel.get('region_end')))
        if (s0,e0)!=(st,en):
            return None, {'status':'legacy_lead_selected_region_not_used','detail':'Run a region-wide PhyML inventory; an LD-pruned subset cannot produce a full-locus best tag'}, paths+[final/'loci.tsv']
        paths.append(final/'loci.tsv'); scope='legacy_catalogue_full_core_bounds_verified'
    all_copies=set()
    for h in haps:
        h['_copies']=copy_ids(h)
        n=integer(h.get('hap_n',h.get('n')),len(h['_copies']))
        if n != len(h['_copies']):
            raise ValueError(f'Copy count mismatch: {h["hap_id"]}')
        if all_copies & h['_copies']:
            raise ValueError('A chromosome copy belongs to multiple catalogue haplotypes')
        all_copies |= h['_copies']
        h['seq']=h['seq'].upper()
    for a in arch:
        a['seq']=a['seq'].upper()
        a['lineage']=lineage(a.get('lineage') or a.get('archaic'))
    in_core=[]; seen_sites=set()
    for i,s in enumerate(sites):
        s['pos']=integer(s['pos'])
        if s['pos'] is None or s['pos']<1:
            raise ValueError('Invalid site position')
        sid=(chrom(s['chr']),s['pos'],s['ref'],s['alt'])
        if sid in seen_sites:
            raise ValueError('Duplicate chromosome/position/REF/ALT site')
        seen_sites.add(sid)
        if chrom(s['chr'])==chrom(auth['chr']) and inside(st,en,s['pos']):
            in_core.append(i)
    if not in_core:
        return None, {'status':'no_stored_sites_in_core'}, paths
    evpath=final/'evidence_haplotypes.tsv';ev=rows(evpath)
    ep=final/'evidence_sites.tsv';diag=rows(ep)
    tp=final/'evidence_trees.tsv';trees=rows(tp)
    for p in [evpath,ep,tp,final/'anchors.tsv',final/'anchor_copies.tsv']:
        if p.is_file():paths.append(p)
    # Membership labels are validated against actual copy sets, not raw H labels.
    byid={h['hap_id']:h for h in haps}
    for e in ev:
        if e.get('locus_id')!=auth['locus_id'] or e.get('hap_id') not in byid:
            continue
        e['_copies']=copy_ids(e)
        if e['_copies']!=byid[e['hap_id']]['_copies']:
            raise ValueError('Evidence and region catalogue H labels map to different copies')
    return dict(final=final,locus=locus,sites=sites,haps=haps,arch=arch,
                in_core=in_core,copies=all_copies,copy_sequences={cp:h['seq'] for h in haps for cp in h['_copies']},evidence=ev,diagnostic=diag,
                trees=trees,scope=scope), {'status':'computed'}, paths


def choose_haplotype(auth,data,min_sites=10,min_callable=.8):
    ev={(r['hap_id'],r.get('diagnostic_lineage')):r for r in data['evidence']
        if r.get('locus_id')==auth['locus_id'] and r.get('hap_id')}
    tree_members={}
    for t in data['trees']:
        if t.get('locus_id')==auth['locus_id'] and truth(t.get('candidate_clade_pass')):
            tree_members[t.get('expected_lineage')]=set(t.get('candidate_tips_in_clade','').split(','))
    if truth(auth.get('candidate_tree_pass')):
        tree_members.setdefault(auth.get('evidence_lineage'),set(auth.get('candidate_tree_candidates_in_clade','').split(',')))
    comparisons=[]
    for h in data['haps']:
        for a in data['arch']:
            pos=[i for i in data['in_core'] if a['seq'][i] in BASES]
            called=[i for i in pos if h['seq'][i] in BASES]
            nc=len(called); nm=sum(h['seq'][i]==a['seq'][i] for i in called)
            if nc<min_sites or not pos or nc/len(pos)<min_callable:
                continue
            e=ev.get((h['hap_id'],a['lineage']),{})
            diagnostic=truth(e.get('diagnostic_candidate_pass'))
            tree=diagnostic and h['hap_id'] in tree_members.get(a['lineage'],set())
            comparisons.append(dict(hap_id=h['hap_id'],archaic=a['archaic'],lineage=a['lineage'],
                                    n_compared=nc,n_match=nm,prop_match=nm/nc,n_copies=len(h['_copies']),
                                    archaic_callable_fraction=nc/len(pos),evidence_tier=3 if tree else 2 if diagnostic else 1,
                                    diagnostic_match_prop=num(e.get('diagnostic_match_prop')),
                                    n_contiguous_diagnostic_match=integer(e.get('n_contiguous_diagnostic_match'),0),
                                    candidate_start=integer(e.get('candidate_start')),
                                    candidate_end=integer(e.get('candidate_end'))))
    if not comparisons:
        return None,None,ev,tree_members
    best=sorted(comparisons,key=rank_evidence)[0]
    raw=sorted(comparisons,key=lambda r:(-r['prop_match'],-r['n_compared'],-r['n_copies'],r['hap_id'],r['archaic']))[0]
    return best,raw,ev,tree_members


def build_variants(auth,data,best):
    h=next(x for x in data['haps'] if x['hap_id']==best['hap_id'])
    a=next(x for x in data['arch'] if x['archaic']==best['archaic'])
    passed_diag={(integer(x.get('pos')),x.get('lineage'),x.get('lineage_consensus'))
                 for x in data['diagnostic'] if truth(x.get('diagnostic_site_pass')) and x.get('locus_id')==auth['locus_id']}
    pool=[]
    for i in data['in_core']:
        s=data['sites'][i]
        ref,alt=s['ref'].upper(),s['alt'].upper()
        # An A/C/G/T alignment cannot reconstruct indels or multiallelic alleles.
        if ref not in BASES or alt not in BASES or ref==alt:
            continue
        allele=h['seq'][i]
        if allele not in {ref,alt} or a['seq'][i]!=allele:
            continue
        amap=SiteAlleles(data['copy_sequences'],i,{ref,alt})
        cs,ce=best['candidate_start'],best['candidate_end']
        pool.append(dict(chr=chrom(s['chr']),pos_1based=s['pos'],snp_id=s['id'] if s['id'] not in MISSING else None,
                         ref=ref,alt=alt,tag_allele=allele,other_allele=alt if allele==ref else ref,
                         variant_key=f"{auth['genome_build']}:{chrom(s['chr'])}:{s['pos']}:{ref}:{alt}",
                         archaic_ref=best['archaic'],archaic_matching_allele=allele,
                         inside_candidate=int(cs is not None and ce is not None and inside(cs,ce,s['pos'])),
                         diagnostic_for_target_lineage=int((s['pos'],best['lineage'],allele) in passed_diag),
                         _alleles=amap))
    return pool


def input_variant(auth,data):
    requested=auth.get('_input_lead_override',auth.get('named_anchor_id',auth.get('name',auth['locus_id'])))
    if not requested or requested in MISSING:
        return requested,None,'input_not_provided'
    hits=[s for s in data['sites'] if requested in aliases(s.get('id'))]
    # Exact CHR:POS:REF:ALT is also accepted; no nearest-SNP substitution.
    if not hits:
        m=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([^:]+):([^:]+)',requested,re.I)
        if m:
            hits=[s for s in data['sites'] if (chrom(s['chr']),s['pos'],s['ref'],s['alt'])==(chrom(m[1]),int(m[2]),m[3],m[4])]
    if len(hits)>1:
        return requested,None,'ambiguous_input_variant_id'
    if hits:
        s=hits[0];i=data['sites'].index(s);valid={s['ref'],s['alt']}
        if len(valid)==2 and valid<=BASES:
            amap=SiteAlleles(data['copy_sequences'],i,valid)
            return requested,dict(chr=chrom(s['chr']),pos_1based=s['pos'],snp_id=requested,ref=s['ref'],alt=s['alt'],
                                   variant_key=f"{auth['genome_build']}:{chrom(s['chr'])}:{s['pos']}:{s['ref']}:{s['alt']}",
                                   _alleles=amap),'observed_in_phased_alignment'
    anchors=[r for r in rows(data['final']/'anchors.tsv') if r.get('locus_id')==auth['locus_id']
             and requested in aliases(r.get('anchor_id')) | aliases(r.get('requested_anchor_id'))]
    if len(anchors)>1:
        return requested,None,'ambiguous_input_anchor_record'
    if anchors and truth(anchors[0].get('exact_anchor_found')):
        s=anchors[0];pos=integer(s.get('pos'));ref=s.get('ref','');alt=s.get('alt','')
        rec=dict(chr=chrom(s.get('chr',auth['chr'])),pos_1based=pos,snp_id=requested,ref=ref,alt=alt,
                 variant_key=f"{auth['genome_build']}:{chrom(s.get('chr',auth['chr']))}:{pos}:{ref}:{alt}",_alleles={})
        if ',' in alt or not ref or not alt:
            return requested,rec,'multiallelic_input_not_scored'
        for c in rows(data['final']/'anchor_copies.tsv'):
            if c.get('locus_id')!=auth['locus_id'] or not (aliases(c.get('anchor_id')) & (aliases(s.get('anchor_id')) | {requested})):
                continue
            if integer(c.get('anchor_pos'))!=pos or not truth(c.get('assigned_to_hap_id')):
                continue
            cp=str(c.get('sample',c.get('sample_id')))+':'+str(c.get('haplotype'))
            val=c.get('allele')
            if cp in data['copies'] and val in {ref,alt}:
                old=rec['_alleles'].get(cp)
                if old is not None and old!=val:
                    raise ValueError('Conflicting phased input-marker allele')
                rec['_alleles'][cp]=val
        return requested,rec,'observed_in_anchor_copies' if rec['_alleles'] else 'input_found_phase_unavailable'
    return requested,None,'input_not_found_no_substitution'


def clean(r):
    return {k:v for k,v in r.items() if not k.startswith('_')}


def empty_metrics():
    return {k:None for k in binary_metrics(0,0,0,0)}


def marker_metrics(record,target,cohort,allele=None):
    if record is None or not record.get('_alleles'):
        return dict(tag_allele=allele,**empty_metrics())
    options=[allele] if allele else [record['ref'],record['alt']]
    scored=[dict(tag_allele=v,**count_variant(record['_alleles'],target,cohort,v)) for v in options]
    # Freeze input orientation from ALL and re-use it in population comparisons.
    return max(scored,key=lambda x:(x['r'] if x['r'] is not None else -2,x['f1'] or 0,x['tag_allele']))


def pair_r2(a,b,cohort,aa,ba):
    if a is None or b is None or not aa or not ba:
        return None
    common=cohort & set(a['_alleles']) & set(b['_alleles'])
    first={cp for cp in common if a['_alleles'][cp]==aa}
    return count_variant(b['_alleles'],first,common,ba)['r2']


def analyse_locus(auth,population_rows,options):
    lk=coord_key(auth)
    ident=dict(locus_key=lk,dataset_id=auth['dataset_id'],genome_build=auth['genome_build'],
               locus_id=auth['locus_id'],chr=chrom(auth['chr']),core_start=integer(auth['core_start']),core_end=integer(auth['core_end']))
    requested=auth.get('_input_lead_override',auth.get('named_anchor_id',auth.get('name',auth['locus_id'])))
    summary=dict(ident,input_lead_snp=requested,best_haplotype_id=None,best_haplotype_tier=None,
                 best_haplotype_lineage=None,best_haplotype_archaic=None,best_tag_snp=None,
                 best_tag_pos_1based=None,best_tag_allele=None,best_tag_r2=None,
                 best_tag_ppv=None,best_tag_sensitivity=None,best_tag_quality=None,
                 input_tag_r2=None,input_tag_allele=None,input_pos_1based=None,
                 input_vs_best_snp_r2=None,same_input_and_best=None,n_equivalent_best_tags=0,
                 family_tag_snp=None,family_tag_r2=None,tag_status='not_computed',tag_detail='')
    result=dict(summary=summary,comparisons=[],ranked=[],population=[],by_population=[],haplotypes=[],gwas=[],files=[])
    data,state,files=load_data(auth);result['files']=files
    if data is None:
        summary.update(tag_status=state['status'],tag_detail=state.get('detail',''))
        return result
    best,raw,ev,members=choose_haplotype(auth,data,options.min_match_sites,options.min_match_callable)
    if best is None:
        summary.update(tag_status='no_haplotype_meets_sequence_information_minimum',tag_detail='No SNP or lead-SNP substitution is made')
        return result
    targets={}
    selected=next(h for h in data['haps'] if h['hap_id']==best['hap_id'])
    targets['best_haplotype']=selected['_copies']
    fam=set()
    for h in data['haps']:
        e=ev.get((h['hap_id'],best['lineage']),{})
        if best['evidence_tier']==3:
            ok=truth(e.get('diagnostic_candidate_pass')) and h['hap_id'] in members.get(best['lineage'],set())
        else:
            ok=truth(e.get('diagnostic_candidate_pass'))
        if ok:fam|=h['_copies']
    if fam:
        targets['candidate_family']=fam
    groups=make_groups(data['copies'],population_rows)
    pool=build_variants(auth,data,best)
    req,iv,input_status=input_variant(auth,data)
    summary.update(best_haplotype_id=best['hap_id'],best_haplotype_tier={3:'tree_supported_candidate',2:'sequence_candidate',1:'raw_similarity_only'}[best['evidence_tier']],
                   best_haplotype_lineage=best['lineage'],best_haplotype_archaic=best['archaic'],
                   best_haplotype_copies=best['n_copies'],best_haplotype_prop_match=best['prop_match'],
                   n_stored_region_sites=len(data['in_core']),n_archaic_matching_snp_candidates=len(pool),
                   n_tested_haplotype_copies=len(data['copies']),sequence_scope=data['scope'],
                   raw_best_haplotype_id=raw['hap_id'],raw_best_archaic=raw['archaic'],raw_best_prop_match=raw['prop_match'],
                   input_status=input_status,input_pos_1based=iv['pos_1based'] if iv else None,
                   raw_match_is_not_introgression=int(best['evidence_tier']==1))
    result['haplotypes']=[dict(ident,role='evidence_prioritized_best',**best),dict(ident,role='raw_similarity_best',**raw)]
    for target_name,target in targets.items():
        ranked=[]
        for v in pool:
            m=count_variant(v['_alleles'],target,data['copies'],v['tag_allele'])
            if qualify(m,options.min_tag_call_rate,options.min_tag_target_copies)=='eligible':
                ranked.append(dict(v,**m))
        ranked.sort(key=rank_tag)
        chosen=ranked[0] if ranked else None
        ties=[]
        if chosen:
            ties=[r for r in ranked if (round(r['r2'],12),round(r['f1'],12),round(r['call_rate'],12)) ==
                  (round(chosen['r2'],12),round(chosen['f1'],12),round(chosen['call_rate'],12))]
        im=marker_metrics(iv,target,data['copies'])
        sm=marker_metrics(chosen,target,data['copies'],chosen['tag_allele'] if chosen else None)
        inp=dict(ident,target_definition=target_name,role='input_lead',requested_snp=req,status=input_status,
                 **(clean(iv) if iv else {}),**im)
        if iv:
            inp['inside_locus']=int(inside(ident['core_start'],ident['core_end'],iv['pos_1based']))
        bst={**ident,'target_definition':target_name,'role':'best_tag','requested_snp':None,
                 'status':('strong_in_sample_proxy' if chosen and chosen['r2']>=options.strong_tag_r2 else 'weak_in_sample_proxy' if chosen else 'no_eligible_single_snp_tag'),
                 **(clean(chosen) if chosen else {}),**sm,'n_equivalent_best_tags':len(ties)}
        result['comparisons'] += [inp,bst]
        for rank,r in enumerate(ranked[:options.top_tags],1):
            result['ranked'].append(dict(ident,target_definition=target_name,rank=rank,**clean(r),
                                          tied_on_statistical_metrics=int(r in ties)))
        # Fixed ALL-selected alleles/tags: these rows quantify portability.
        for g,cohort in sorted(groups.items()):
            for role,v,allele in [('input_lead',iv,im['tag_allele']),('best_tag',chosen,chosen['tag_allele'] if chosen else None)]:
                met=marker_metrics(v,target,cohort,allele)
                result['population'].append(dict(ident,target_definition=target_name,population=g,
                    role=role,snp_id=v.get('snp_id') if v else req if role=='input_lead' else None,
                    variant_key=v.get('variant_key') if v else None,**met,
                    metric_status=qualify(met,options.min_tag_call_rate,options.min_tag_target_copies) if v and v.get('_alleles') else 'not_evaluable',
                    selection_population='ALL',tag_reselected=0))
        # Optional per-superpopulation reselection: same global haplotype target,
        # not a new ancestry-specific haplotype definition.
        for g,cohort in sorted(groups.items()):
            if g=='ALL' or (g.startswith('POP:') and options.tag_population_level!='all') or options.tag_population_level=='none':
                continue
            candidates=[]
            for v in pool:
                met=count_variant(v['_alleles'],target,cohort,v['tag_allele'])
                if qualify(met,options.min_tag_call_rate,options.min_tag_target_copies)=='eligible':
                    candidates.append(dict(v,**met))
            candidates.sort(key=rank_tag)
            bb=candidates[0] if candidates else None
            result['by_population'].append(dict(ident,target_definition=target_name,population=g,
                   status='computed' if bb else 'not_evaluable_or_no_single_snp_tag',
                   **(clean(bb) if bb else {}),tag_reselected=1))
        for role,v,met in [('input_lead',iv,im),('best_tag',chosen,sm)]:
            if v and v.get('pos_1based'):
                result['gwas'].append(dict(ident,target_definition=target_name,role=role,
                    snp_id=v.get('snp_id'),variant_key=v.get('variant_key'),pos_1based=v['pos_1based'],
                    ref=v['ref'],alt=v['alt'],tag_allele=met['tag_allele'],tag_haplotype_r2=met['r2'],
                    gwas_effect_allele=None,gwas_other_allele=None,gwas_beta=None,gwas_se=None,gwas_p=None,
                    gwas_status='not_queried_user_followup',allele_note='tag_allele_is_not_assumed_to_be_risk_or_effect_allele'))
        if target_name=='best_haplotype':
            summary.update(best_tag_snp=(chosen['snp_id'] or chosen['variant_key']) if chosen else None,
                           best_tag_pos_1based=chosen['pos_1based'] if chosen else None,
                           best_tag_allele=sm['tag_allele'],best_tag_r2=sm['r2'],best_tag_ppv=sm['ppv'],best_tag_sensitivity=sm['sensitivity'],
                           best_tag_quality=bst['status'],input_tag_r2=im['r2'],input_tag_allele=im['tag_allele'],
                           input_vs_best_snp_r2=pair_r2(iv,chosen,data['copies'],im['tag_allele'],sm['tag_allele']),
                           same_input_and_best=int(iv['variant_key']==chosen['variant_key']) if iv and chosen else None,
                           n_equivalent_best_tags=len(ties),tag_status='computed',
                           tag_detail='Full stored locus SNP scan; no GWAS or input-lead constraint. r2 is SNP--haplotype, not SNP--SNP. In-sample, not externally validated.')
        else:
            summary.update(family_tag_snp=(chosen['snp_id'] or chosen['variant_key']) if chosen else None,family_tag_r2=sm['r2'])
    return result


def analyse_all(authoritative,population_rows,options):
    out={k:[] for k in ['summary','comparisons','ranked','population','by_population','haplotypes','gwas','files']}
    for auth in authoritative:
        try:
            r=analyse_locus(auth,[p for p in population_rows if p.get('dataset_id',auth['dataset_id'])==auth['dataset_id']],options)
        except (OSError,ValueError,KeyError,IndexError,csv.Error) as exc:
            # Do not take down the other loci. Invalid inputs must not produce tags.
            r={k:[] for k in out}
            r['summary']=dict(locus_key=coord_key(auth),dataset_id=auth['dataset_id'],genome_build=auth['genome_build'],
                              locus_id=auth['locus_id'],chr=chrom(auth['chr']),core_start=integer(auth['core_start']),core_end=integer(auth['core_end']),
                              input_lead_snp=auth.get('_input_lead_override',auth.get('named_anchor_id',auth['locus_id'])),
                              tag_status='invalid_sequence_inputs',tag_detail=str(exc))
        out['summary'].append(r['summary'])
        for k in out:
            if k!='summary':out[k]+=r[k]
    out['files']=sorted(set(out['files']))
    return out


def apply_loci_file(authoritative,path,coordinate_format='bed'):
    """Exact coordinate matching only. Column 4 is output metadata, never a filter.

    Existing sequence results must already cover the requested region. An input
    file is not allowed to silently redefine or widen a previously computed locus.
    """
    requests={}
    with Path(path).open(encoding='utf-8-sig') as f:
        for ln,text in enumerate(f,1):
            if not text.strip() or text.lstrip().startswith(('#','track','browser')):continue
            fields=text.split()
            if len(fields)<3:raise ValueError(f'Invalid loci line {ln}')
            try:st,en=int(fields[1]),int(fields[2])
            except ValueError:
                if ln==1 and fields[0].lower() in {'chr','chrom','chromosome'}:continue
                raise ValueError(f'Invalid coordinates on loci line {ln}')
            if coordinate_format=='1based':st-=1
            if st<0 or en<=st:raise ValueError(f'Invalid interval on loci line {ln}')
            ck=(chrom(fields[0]),st,en)
            name=fields[3] if len(fields)>3 else ''
            if ck in requests and requests[ck]!=name:raise ValueError('Two different input lead labels for the same coordinates')
            requests[ck]=name
    result=[];matched=set()
    for r in authoritative:
        ck=(chrom(r['chr']),integer(r['core_start']),integer(r['core_end']))
        if ck in requests:
            r=dict(r);r['_input_lead_override']=requests[ck];result.append(r);matched.add(ck)
    missing=set(requests)-matched
    if missing:raise ValueError('Requested regions lack an exact saved PhyML coordinate match: '+str(sorted(missing)))
    if not result:raise ValueError('No loci in input file')
    return result
