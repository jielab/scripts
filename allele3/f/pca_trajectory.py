#!/usr/bin/env python3
"""Temporal structure in frozen shared projected PCA coordinates; never fits PCA."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import gzip
import html
import json
import math
from pathlib import Path
import sys
import numpy as np
from analysis_integrity import atomic_json, identity, transactional, code_identity
from pca_contract import read_scores, historical_identity, coordinate_join
from temporal_stats import (number, read_table, write_table, technology, age_bin, seed_for,
                            bootstrap_mean, fixed_composition, linear_model, adjust_fdr, NA)

LIMITATIONS = [
    'PC1/PC2/age BP are the original three axes; time is not PC3 and never enters the genotype PCA fit.',
    'PC changes can reflect migration, mixture and sampling. No selection, ancestral path, or evolutionary endpoint is inferred.',
    'Intervals resample unrelated individuals conditional on frozen coordinates, the supplied cohort and declared composition.',
    'Source is a measurement category, not an independent study/site cluster; no invented clusters or dating errors.',
    'PCA dispersion is not nucleotide diversity. Modern samples are a separate endpoint, never ancient age-zero observations.',
    'Shotgun labels do not establish independent, unimputed whole-genome validation. Genotype encoding needs source documentation.',
    'No missing-marker reprojection, new basis fitting, or independent read validation is performed by this coordinate module.',
]


def read_unique(path):
    result = {}
    for row in read_table(path):
        sid = row.get('sample_id')
        if not sid or sid in result: raise ValueError('Missing/duplicate metadata sample_id: ' + str(sid))
        result[sid] = row
    if not result: raise ValueError('Empty metadata table')
    return result


def training_members(args, contract):
    basis = Path(contract['basis_dir']) if contract.get('basis_dir') else None
    path = Path(args.training_ids) if args.training_ids else basis / 'basis.keep' if basis else None
    if not path or not path.is_file(): return set(), {'status': 'unavailable'}
    recorded = None
    if basis:
        plan = json.loads((basis / 'pca_plan.json').read_text())
        matched = [digest for name, digest in plan.get('plan_file_sha256', {}).items() if Path(name).name == path.name]
        if len(matched) == 1: recorded = {'sha256': matched[0]}
    receipt = historical_identity(path, args.provenance_manifest, recorded)
    ids = []
    with path.open() as handle:
        for line in handle:
            parts = line.split()
            if not parts or parts[0].startswith('#') or parts[-1] == 'IID': continue
            ids.append(parts[-1])
    if len(set(ids)) != len(ids): raise ValueError('Duplicate training ID')
    return set(ids), dict(receipt, status=receipt['historical_metadata_verification'])


def load_inputs(args):
    requested = args.pcs or 2
    coords, columns, contract = read_scores(args.coordinates, requested, basis_dir=args.basis_dir,
        receipt=args.coordinate_receipt, require_verified=False)
    # Optional prespecified 5/10-PC supplements only when every column exists.
    if args.pcs is None:
        for k in (10, 5):
            try:
                extra, names, detail = read_scores(args.coordinates, k, basis_dir=args.basis_dir,
                    receipt=args.coordinate_receipt, require_verified=False)
            except ValueError as exc:
                if 'Missing consecutive requested PCs' in str(exc): continue
                raise
            coords, columns, contract = extra, names, detail; break
    contract.update(primary_pc_count=requested, supplemental_k=[k for k in (5, 10) if k <= len(columns)])
    meta = read_unique(args.metadata)
    meta_receipt = historical_identity(args.metadata, args.provenance_manifest, contract.get('metadata_recorded'))
    training, train_receipt = training_members(args, contract)
    supplement_info = None
    if args.supplement_metadata:
        supplemental = read_unique(args.supplement_metadata)
        for sid, row in supplemental.items():
            if sid not in meta: raise ValueError('Supplement has unmatched sample_id: ' + sid)
            for key, value in row.items():
                if value in NA or key == 'sample_id': continue
                if key in meta[sid] and meta[sid][key] not in NA and meta[sid][key] != value:
                    raise ValueError('Supplement would rewrite existing metadata: ' + sid + ':' + key)
                meta[sid][key] = value
        supplement_info = identity(args.supplement_metadata, full=True)
    biological = Counter(r.get('individual_id') for r in meta.values() if r.get('individual_id') not in NA | {None})
    rows = []
    # Coordinate order is unchanged, with unmatched metadata explicitly appended.
    ids = list(coords) + [sid for sid in meta if sid not in coords]
    for sid in ids:
        raw = dict(meta.get(sid, {})); age = number(raw.get('age_bp')); cohort = raw.get('cohort', 'unknown')
        region = raw.get('region', raw.get('stratum', 'unknown'))
        if region in NA: region = 'unknown'
        source = raw.get('data_source', raw.get('source', 'unknown'))
        if source in NA: source = 'unknown'
        tech = technology(raw)
        reason = []
        if sid not in meta: reason.append('metadata_missing')
        if sid not in coords: reason.append('coordinates_missing')
        if raw.get('include') != '1': reason.append('include_not_explicit_1')
        if raw.get('unrelated') != '1': reason.append('related' if raw.get('unrelated') == '0' else 'relationship_unknown')
        if cohort not in {'ancient', 'modern'}: reason.append('cohort_unknown')
        if age is None or age < 0: reason.append('invalid_date')
        if region == 'unknown': reason.append('region_unknown')
        if biological.get(raw.get('individual_id'), 0) > 1: reason.append('duplicate_biological_individual')
        formal = not reason
        membership = 'training' if sid in training else 'projected' if train_receipt.get('status') == 'verified' else 'unknown'
        if train_receipt.get('status') != 'verified': membership = 'unknown'
        row = dict(raw, sample_id=sid, IID=sid, cohort=cohort, age_bp=age, original_age_bp=raw.get('age_bp'),
                   region=region, data_source=source, technology=tech, measurement_source=tech + '|' + source,
                   formal=formal, formal_exclusion=';'.join(reason), coordinate_present=int(sid in coords),
                   training_status=membership, individual_audit='known' if raw.get('individual_id') not in NA | {None} else 'not_provided',
                   basis_id=contract.get('basis_id'), coordinate_kind=contract['coordinate_kind'],
                   time_bin=age_bin(age, args.bin_years, args.max_age_bp) if cohort == 'ancient' else None,
                   date_uncertainty='not_provided', iqs_definition='as supplied; not inferred from coverage',
                   genotype_encoding='not_established_from_coordinates')
        for k in range(len(columns)): row[f'PC{k+1}'] = coords[sid][k] if sid in coords else None
        rows.append(row)
    # Retain FID, denominators and measurement fields from the exact original score.
    opener = gzip.open if str(args.coordinates).endswith('.gz') else open
    score_info = {}
    with opener(args.coordinates, 'rt') as handle:
        lines = (line for line in handle if line.strip() and not line.startswith('##'))
        names = next(lines).lstrip('#').split()
        for line in lines:
            record = dict(zip(names, line.split())); score_info[record['IID']] = record
    for row in rows:
        record = score_info.get(row['sample_id'], {})
        for key in ('FID', 'ALLELE_CT', 'DENOM', 'NAMED_ALLELE_DOSAGE_SUM', 'marker_count', 'call_fraction'):
            if key in record: row[key] = record[key]
    return rows, coords, columns, dict(contract, join=coordinate_join(list(meta.values()), coords),
        metadata_provenance=meta_receipt, training_provenance=train_receipt, supplement=supplement_info)


def composition(samples, width, maximum):
    groups = defaultdict(list)
    for s in samples:
        epoch = age_bin(s['age_bp'], width, maximum) if s['cohort'] == 'ancient' else 'modern' if s['cohort'] == 'modern' else 'unknown'
        if epoch is None: epoch = 'older_or_invalid'
        groups[s['region'], epoch, s['cohort'], s['measurement_source']].append(s)
    result = []
    for (region, epoch, cohort, source), ss in sorted(groups.items(), key=lambda x: str(x[0])):
        row = dict(region=region, time_bin=epoch, cohort=cohort, measurement_source=source, n=len(ss),
                   n_formal=sum(s['formal'] for s in ss), n_related=sum(s.get('unrelated') == '0' for s in ss),
                   n_relationship_unknown=sum(s.get('unrelated') not in {'0','1'} for s in ss),
                   n_training=sum(s['training_status'] == 'training' for s in ss),
                   n_projected=sum(s['training_status'] == 'projected' for s in ss),
                   n_modern=sum(s['cohort'] == 'modern' for s in ss),
                   exclusion_counts=json.dumps(dict(Counter(r for s in ss for r in s['formal_exclusion'].split(';') if r)), sort_keys=True))
        for name in ('target_coverage', 'iqs', 'marker_count', 'call_fraction', 'ALLELE_CT', 'DENOM'):
            values = [number(s.get(name)) for s in ss]; values = [v for v in values if v is not None]
            row.update({name+'_n': len(values), name+'_median': float(np.median(values)) if values else None,
                        name+'_q25': float(np.quantile(values,.25)) if values else None, name+'_q75': float(np.quantile(values,.75)) if values else None})
        result.append(row)
    return result


def scaling(samples, k):
    reference = [s for s in samples if s['coordinate_present'] and s['training_status'] == 'training']
    definition = 'verified_frozen_training_members'
    if len(reference) < 2:
        reference = [s for s in samples if s['formal'] and s['coordinate_present']]
        definition = 'fixed_entire_formal_reference_cohort; training_members_unavailable'
    if len(reference) < 2: raise ValueError('At least two fixed reference samples required for PC scaling')
    matrix = np.array([[s[f'PC{i+1}'] for i in range(k)] for s in reference])
    scale = matrix.std(axis=0, ddof=1)
    result = [dict(pc=f'PC{i+1}', mean=float(matrix[:,i].mean()), scale=float(scale[i]), n=len(reference),
                   reference_population=definition, sample_ids=';'.join(s['sample_id'] for s in reference),
                   status='usable' if scale[i] > 1e-15 else 'zero_reference_variance') for i in range(k)]
    scale[scale <= 1e-15] = np.nan
    return scale, result


def trajectory_tables(samples, k, args, scales):
    centroids, dispersions, contrasts, support, cache = [], [], [], [], {}
    regions = sorted({s['region'] for s in samples if s['region'] != 'unknown'})
    epochs = range(math.ceil(args.max_age_bp / args.bin_years))
    for region in ['all', *regions]:
        eligible = [s for s in samples if s['formal'] and (region == 'all' or s['region'] == region)]
        modern = [s for s in eligible if s['cohort'] == 'modern']
        for epoch in [*epochs, 'modern']:
            ss = modern if epoch == 'modern' else [s for s in eligible if s['cohort'] == 'ancient' and s['time_bin'] == epoch]
            n = len(ss); ok = n >= args.min_samples
            base = dict(analysis_set='formal_include_unrelated', region=region, cohort='modern' if epoch == 'modern' else 'ancient',
                        time_bin=epoch, age_start_bp=None if epoch == 'modern' else epoch*args.bin_years,
                        age_end_bp=None if epoch == 'modern' else min((epoch+1)*args.bin_years,args.max_age_bp),
                        n=n, n_eff=n, status='estimated' if ok else 'not_estimable', reason='' if ok else 'below_min_samples',
                        weights_id='equal_individual', support_id='fixed_bins', seed=args.seed)
            support.append(dict(base, minimum=args.min_samples, modern_reference_n=len(modern),
                                modern_reference_status='available' if len(modern) >= args.min_samples else 'no_matched_modern_support'))
            values = np.array([[s[f'PC{i+1}'] for i in range(k)] for s in ss]).reshape(n,k)
            se, lo, hi, draws = bootstrap_mean(values,args.bootstrap,seed_for(args.seed,(region,epoch))) if ok else (None,None,None,None)
            if ok: cache[region,epoch] = (values.mean(axis=0),draws,n)
            for i in range(k):
                centroids.append(dict(base, measure='raw_mean', pc_or_k=f'PC{i+1}', estimate=float(values[:,i].mean()) if ok else None,
                    descriptive_mean=float(values[:,i].mean()) if n else None, robust_center=float(np.median(values[:,i])) if n else None,
                    individual_q25=float(np.quantile(values[:,i],.25)) if n else None, individual_q75=float(np.quantile(values[:,i],.75)) if n else None,
                    age_median=float(np.median([s['age_bp'] for s in ss])) if n else None,
                    actual_age_min=min((s['age_bp'] for s in ss),default=None),actual_age_max=max((s['age_bp'] for s in ss),default=None),
                    se=float(se[i]) if se is not None else None,ci_low=float(lo[i]) if lo is not None else None,ci_high=float(hi[i]) if hi is not None else None,
                    sample_ids=';'.join(s['sample_id'] for s in ss), uncertainty='mean_bootstrap_CI; individual_quartiles_are_distribution'))
            for dimension in sorted({2, k, *[v for v in (5,10) if v <= k]}):
                if dimension > k: continue
                valid = ok and np.isfinite(scales[:dimension]).all()
                z = values[:,:dimension]/scales[:dimension] if valid else None
                dispersion = float(np.mean(np.sum((z-z.mean(axis=0))**2,axis=1))) if valid else None
                dbase = dict(base,pc_or_k=dimension)
                if not valid: dbase.update(status='not_estimable',reason=base['reason'] or 'zero_reference_variance')
                dispersions.append(dict(dbase,measure='mean_squared_distance_to_centroid_in_fixed_reference_scale',estimate=dispersion))
                matched = valid and epoch != 'modern' and region != 'all' and len(modern) >= args.min_samples
                distance = None
                if matched:
                    endpoint = np.array([[s[f'PC{i+1}'] for i in range(dimension)] for s in modern]).mean(axis=0)
                    distance = float(np.linalg.norm((values[:,:dimension].mean(axis=0)-endpoint)/scales[:dimension]))
                dispersions.append(dict(dbase,measure='distance_to_region_matched_modern',estimate=distance,
                    status='estimated' if matched else 'not_estimable',reason='' if matched else 'modern_separate_or_no_region_matched_support', modern_reference_n=len(modern)))
        # Adjacent bins only; gaps are never bridged.
        for epoch in list(epochs)[1:]:
            left,right = cache.get((region,epoch-1)),cache.get((region,epoch))
            for i in range(k):
                ok = left is not None and right is not None
                draws = right[1][:,i]-left[1][:,i] if ok and right[1] is not None and left[1] is not None else None
                contrasts.append(dict(region=region,cohort='ancient',analysis_set='formal_include_unrelated',time_bin=epoch,
                    age_start_bp=(epoch-1)*args.bin_years,age_end_bp=min((epoch+1)*args.bin_years,args.max_age_bp),
                    measure='older_bin_minus_adjacent_younger_bin',pc_or_k=f'PC{i+1}',estimate=float(right[0][i]-left[0][i]) if ok else None,
                    ci_low=float(np.quantile(draws,.025)) if draws is not None else None,ci_high=float(np.quantile(draws,.975)) if draws is not None else None,
                    n=left[2]+right[2] if ok else 0,n_eff=left[2]+right[2] if ok else None,
                    status='estimated' if ok else 'not_estimable',reason='' if ok else 'adjacent_bin_support_gap',seed=args.seed))
    return centroids,dispersions,contrasts,support


def standardized_tables(samples, k, args):
    primary = [s for s in samples if s['formal'] and s['cohort']=='ancient' and s['time_bin'] is not None]
    result, weights_out, support = [], [], []
    definitions = [('fixed_region','all','region',primary),('fixed_tech_source','all','measurement_source',primary)]
    definitions += [('fixed_tech_source',r,'measurement_source',[s for s in primary if s['region']==r]) for r in sorted({s['region'] for s in primary})]
    for method,region,key,ss in definitions:
        if not ss: continue
        counts = Counter(s[key] for s in ss)
        weights = {g:n/len(ss) for g,n in sorted(counts.items())}
        support_id=f'{method}:{region}:0-{args.max_age_bp}:{args.bin_years}'
        for g,w in weights.items():
            weights_out.append(dict(weights_id=support_id,stratum=g,weight=w,n_reference=counts[g],reference_population='entire_formal_ancient_window',region=region, fixed_across_all_bins=True))
        estimated = fixed_composition([[s[f'PC{i+1}'] for i in range(k)] for s in ss],[s['time_bin'] for s in ss],[s[key] for s in ss],
            weights,range(math.ceil(args.max_age_bp/args.bin_years)),args.min_samples,0 if args.command == 'trajectory-plan' else args.bootstrap,args.seed)
        for e in estimated:
            base={key:value for key,value in e.items() if key not in {'estimate','ci_low','ci_high','se','stratum_n'}}
            base.update(region=region,cohort='ancient',analysis_set='formal_include_unrelated',measure=method,weights_id=support_id,support_id=support_id,
                age_start_bp=e['time_bin']*args.bin_years,age_end_bp=min((e['time_bin']+1)*args.bin_years,args.max_age_bp),seed=args.seed)
            support.append(dict(base,stratum_n=json.dumps(e['stratum_n'],sort_keys=True)))
            for i in range(k): result.append(dict(base,pc_or_k=f'PC{i+1}',**{key:float(e[key][i]) if e[key] is not None else None for key in ('estimate','se','ci_low','ci_high')}))
    return result,weights_out,support


def models(samples, k, verified, args):
    primary = [s for s in samples if s['formal'] and s['cohort']=='ancient' and s['time_bin'] is not None]
    results = []
    for region in ['all', *sorted({s['region'] for s in primary})]:
        subset = [s for s in primary if region=='all' or s['region']==region]
        for i in range(k):
            fit=linear_model(subset,[s[f'PC{i+1}'] for s in subset],minimum=args.min_samples) if verified else dict(status='not_estimable',reason='unverified_basis_or_metadata',n=len(subset),estimate=None,p=None,q=None)
            results.append(dict(fit,region=region,cohort='ancient',measure='PC_response_time_association',pc_or_k=f'PC{i+1}',
                analysis_set='formal_include_unrelated',age_start_bp=0,age_end_bp=args.max_age_bp,seed=args.seed,
                hypothesis_family='all_prespecified_PC_by_region_and_global_linear_models', genetic_PC_covariates='none_PC_is_response'))
    return adjust_fdr(results)


def sensitivities(samples, k, args, verified):
    rows=[]
    for name,predicate in [('all_formal',lambda s:True),('coverage_ge_1',lambda s:(number(s.get('target_coverage')) or 0)>=1),
                           ('shotgun_label',lambda s:s['technology']=='shotgun_label'),('training_only',lambda s:s['training_status']=='training'),
                           ('projected_only',lambda s:s['training_status']=='projected')]:
        subset=[s for s in samples if s['formal'] and s['cohort']=='ancient' and s['time_bin'] is not None and predicate(s)]
        for i in range(min(k,2)):
            fit=linear_model(subset,[s[f'PC{i+1}'] for s in subset],minimum=args.min_samples) if verified else dict(status='not_estimable',reason='unverified_basis_or_metadata',n=len(subset),estimate=None)
            rows.append(dict(fit,analysis_set=name,measure='technical_subset_time_slope',pc_or_k=f'PC{i+1}',region='all',cohort='ancient',seed=args.seed))
    for field in ('target_coverage','marker_count','call_fraction','ALLELE_CT','DENOM'):
        for i in range(min(k,2)):
            pairs=[(number(s.get(field)),s[f'PC{i+1}']) for s in samples if s['formal'] and number(s.get(field)) is not None]
            valid=len(pairs)>=3 and np.std(np.array(pairs)[:,0])>0 and np.std(np.array(pairs)[:,1])>0
            rows.append(dict(analysis_set='formal',measure='descriptive_PC_correlation_'+field,pc_or_k=f'PC{i+1}',n=len(pairs),
                            estimate=float(np.corrcoef(np.array(pairs).T)[0,1]) if valid else None,
                            status='descriptive' if valid else 'not_estimable',reason='' if valid else 'missing_or_constant_covariate'))
    for width,offset in [(500,0),(2000,0),(args.bin_years,args.bin_years/2)]:
        groups=defaultdict(list)
        for s in samples:
            if s['formal'] and s['cohort']=='ancient' and s['age_bp'] is not None and 0<=s['age_bp']<=args.max_age_bp:
                groups[s['region'], math.floor((s['age_bp']-offset)/width)].append(s)
        for (region,epoch),ss in sorted(groups.items()):
            for i in range(min(k,2)):
                rows.append(dict(analysis_set=f'bin_{width}_offset_{offset:g}',region=region,cohort='ancient',measure='raw_mean_boundary_sensitivity',
                    age_start_bp=max(0,offset+epoch*width),age_end_bp=min(args.max_age_bp,offset+(epoch+1)*width),pc_or_k=f'PC{i+1}',n=len(ss),
                    estimate=float(np.mean([s[f'PC{i+1}'] for s in ss])) if len(ss)>=args.min_samples else None,
                    status='descriptive' if len(ss)>=args.min_samples else 'not_estimable',reason='' if len(ss)>=args.min_samples else 'below_min_samples'))
    for name,need in [('missing_marker_reprojection','actual missing-marker masks and explicit separate projection'),('dosage_vs_hardcall','actual dosage-bearing genotype files'),('independent_read_validation','independent reads and study IDs'),('alternative_basis','explicit new basis fit and its receipts')]:
        rows.append(dict(analysis_set=name,status='not_run',reason=need,measure='conditional_extension'))
    return rows


def date_sensitivity(samples, args):
    if not args.date_distributions and not args.date_interval_uniform:
        return [], {'status':'not_provided','point_ages_preserved':True}
    distributions=json.loads(Path(args.date_distributions).read_text()) if args.date_distributions else {}
    known={s['sample_id'] for s in samples}
    if set(distributions)-known: raise ValueError('Date distribution contains unmatched sample IDs')
    prepared={}; modern_ignored=[]
    for s in samples:
        sid=s['sample_id']
        if s['cohort']!='ancient':
            if sid in distributions: modern_ignored.append(sid)
            continue
        if sid in distributions:
            d=distributions[sid]; ages=np.asarray(d['ages_bp'],dtype=float); probability=np.asarray(d['probabilities'],dtype=float)
            if ages.ndim!=1 or not len(ages) or ages.shape!=probability.shape or not np.isfinite(ages).all() or np.any(ages<0) or not np.isfinite(probability).all() or np.any(probability<0) or not math.isclose(float(probability.sum()),1.,abs_tol=1e-8):
                raise ValueError('Invalid calibrated age distribution: '+sid)
            prepared[sid]=('supplied_distribution',ages,probability/probability.sum())
        elif args.date_interval_uniform:
            lo,hi=number(s.get('age_lower_bp')),number(s.get('age_upper_bp'))
            if lo is not None and hi is not None:
                if lo<0 or hi<lo: raise ValueError('Invalid supplied age interval: '+sid)
                prepared[sid]=('assumed_uniform_interval',lo,hi)
    for s in samples:
        if s['sample_id'] in prepared: s['date_uncertainty'] = prepared[s['sample_id']][0]
    rng=np.random.default_rng(args.seed); rows=[]
    for draw in range(args.date_draws):
        grouped=defaultdict(list)
        for s in samples:
            if not s['formal'] or s['cohort']!='ancient': continue
            d=prepared.get(s['sample_id']); age=s['age_bp']
            if d: age=float(rng.choice(d[1],p=d[2])) if d[0]=='supplied_distribution' else float(rng.uniform(d[1],d[2]))
            epoch=age_bin(age,args.bin_years,args.max_age_bp)
            if epoch is not None: grouped[s['region'],epoch].append(s)
        for (region,epoch),ss in sorted(grouped.items()):
            for i in (1,2): rows.append(dict(draw=draw,region=region,time_bin=epoch,pc_or_k=f'PC{i}',n=len(ss),
                estimate=float(np.mean([s[f'PC{i}'] for s in ss])) if len(ss)>=args.min_samples else None,
                status='date_sensitivity' if len(ss)>=args.min_samples else 'not_estimable',seed=args.seed))
    return rows,dict(status='sensitivity_only',n_supplied=len(prepared),methods=dict(Counter(d[0] for d in prepared.values())),
        modern_not_sampled=modern_ignored,draws=args.date_draws,seed=args.seed,point_ages_preserved=True,
        distribution_input=identity(args.date_distributions,full=True) if args.date_distributions else None)


def display_sample(samples, limit, seed):
    """Stable stratified display subset; statistical functions never receive it."""
    if len(samples)<=limit: return samples
    groups=defaultdict(list)
    for s in samples: groups[s['cohort'],s['region'],s.get('time_bin')].append(s)
    rng=np.random.default_rng(seed); selected=[]
    for key,ss in sorted(groups.items(),key=lambda x:str(x[0])):
        n=max(1,int(limit*len(ss)/len(samples)))
        selected.extend(ss[i] for i in rng.choice(len(ss),min(n,len(ss)),replace=False))
    if len(selected)>limit: selected=[selected[i] for i in rng.choice(len(selected),limit,replace=False)]
    return selected


def figures(samples, centroids, dispersion, standardized, sensitivity, composition_rows, out, args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figdir=out/'figures';figdir.mkdir(exist_ok=True)
    points=[s for s in samples if s['coordinate_present'] and s['age_bp'] is not None]
    shown=display_sample(points,args.display_limit,args.seed)
    write_table(figdir/'F1_plot_points.tsv.gz',shown)
    write_table(figdir/'F1_full_points.tsv.gz',points)
    def save(fig,name,caption):
        fig.suptitle(caption,fontsize=14)
        fig.text(.01,.012,f'N coordinates={len(points):,}; shown={len(shown):,}. Statistical tables use full eligible data. See ../manifest.json and sample_audit.tsv.gz.',fontsize=8)
        fig.tight_layout(rect=(0,.04,1,.94))
        for ext in ('png','pdf'):fig.savefig(figdir/f'{name}.{ext}',dpi=180)
        plt.close(fig)
    fig=plt.figure(figsize=(16,9));ax=fig.add_subplot(121,projection='3d');b=fig.add_subplot(122)
    for cohort,color in [('ancient','#1676a3'),('modern','#d16b38')]:
        ss=[s for s in shown if s['cohort']==cohort]
        ax.scatter([s['PC1'] for s in ss],[s['PC2'] for s in ss],[s['age_bp'] for s in ss],s=4,alpha=.4,label=f'{cohort} n={len(ss):,}',color=color)
        b.scatter([s['PC1'] for s in ss],[s['PC2'] for s in ss],s=4,alpha=.4,label=cohort,color=color)
    ax.set(xlabel='PC1 projected AVG',ylabel='PC2 projected AVG',zlabel='Original age BP (time)');b.set(xlabel='PC1 projected AVG',ylabel='PC2 projected AVG');ax.legend();b.legend()
    save(fig,'F1_shared_coordinates','Fixed genetic coordinates with original age; modern endpoint shown separately')
    write_table(figdir/'F2_composition.tsv',composition_rows)
    fig,axes=plt.subplots(2,2,figsize=(16,9));regions=sorted({r['region'] for r in composition_rows});bins=range(math.ceil(args.max_age_bp/args.bin_years))
    matrix=np.array([[sum(r['n_formal'] for r in composition_rows if r['region']==reg and r['cohort']=='ancient' and r['time_bin']==ep) for ep in bins] for reg in regions])
    axes[0,0].imshow(matrix,aspect='auto');axes[0,0].set(yticks=range(len(regions)),yticklabels=regions,xlabel='Ancient time bin',title='Formal sample N by region and time')
    for source in sorted({s['technology'] for s in samples}):
        counts=[sum(s['formal'] and s['cohort']=='ancient' and s['time_bin']==ep and s['technology']==source for s in samples) for ep in bins]
        axes[0,1].plot(list(bins),counts,label=source)
    axes[0,1].set(xlabel='Ancient time bin',ylabel='N',title='Technology labels');axes[0,1].legend(fontsize=8)
    ss=[s for s in shown if number(s.get('target_coverage')) is not None and s['cohort']=='ancient']
    axes[1,0].scatter([s['age_bp'] for s in ss],[float(s['target_coverage']) for s in ss],s=5,alpha=.3);axes[1,0].set(xlabel='Age BP',ylabel='Reported target coverage')
    reasons=Counter(r for s in samples for r in s['formal_exclusion'].split(';') if r)
    axes[1,1].barh(list(reasons),list(reasons.values()));axes[1,1].set(title='Formal exclusions (reasons may overlap)',xlabel='N')
    save(fig,'F2_sample_composition','Sampling and measurement audit; unknown relationships excluded from formal estimates')
    write_table(figdir/'F3_centroids.tsv',centroids)
    fig,axes=plt.subplots(1,2,figsize=(16,9))
    for i,ax in enumerate(axes,1):
        for region in regions:
            rows=[r for r in centroids if r['region']==region and r['cohort']=='ancient' and r['pc_or_k']==f'PC{i}']
            x=[(r['age_start_bp']+r['age_end_bp'])/2 for r in rows];y=[r['estimate'] if r['estimate'] is not None else np.nan for r in rows]
            ax.plot(x,y,marker='o',ms=3,label=region)
            ax.fill_between(x,[r['ci_low'] if r['ci_low'] is not None else np.nan for r in rows],[r['ci_high'] if r['ci_high'] is not None else np.nan for r in rows],alpha=.12)
        ss=[s for s in shown if s['formal'] and s['cohort']=='ancient' and s['time_bin'] is not None]
        ax.scatter([s['age_bp'] for s in ss],[s[f'PC{i}'] for s in ss],color='gray',s=2,alpha=.12,label='Individual distribution')
        ax.set(xlabel='Age BP',ylabel=f'PC{i} projected AVG',title=f'PC{i}: mean and individual distribution');ax.legend(fontsize=8)
    save(fig,'F3_regional_trajectories','Within-region means; shaded intervals are bootstrap CIs for means (not individual ellipses)')
    write_table(figdir/'F4_standardized.tsv',standardized);write_table(figdir/'F4_dispersion.tsv',dispersion)
    fig,axes=plt.subplots(2,2,figsize=(16,9))
    for i in (1,2):
        ax=axes[0,i-1]
        for method in ('raw','fixed_region','fixed_tech_source'):
            rows=[r for r in (centroids if method=='raw' else standardized) if r['region']=='all' and r['cohort']=='ancient' and r['pc_or_k']==f'PC{i}' and (method=='raw' or r['measure']==method)]
            ax.plot([(r['age_start_bp']+r['age_end_bp'])/2 for r in rows],[r['estimate'] if r['estimate'] is not None else np.nan for r in rows],marker='o',label=method)
        ax.set(xlabel='Age BP',ylabel=f'PC{i} projected AVG');ax.legend()
    for ax,measure in zip(axes[1],('mean_squared_distance_to_centroid_in_fixed_reference_scale','distance_to_region_matched_modern')):
        for region in regions:
            rows=[r for r in dispersion if r['region']==region and r['cohort']=='ancient' and r['pc_or_k']==2 and r['measure']==measure]
            ax.plot([(r['age_start_bp']+r['age_end_bp'])/2 for r in rows],[r['estimate'] if r['estimate'] is not None else np.nan for r in rows],label=region)
        ax.set(xlabel='Age BP',ylabel='Fixed reference scale',title='Within-bin dispersion' if measure.startswith('mean') else 'Distance to region-matched modern (NA without support)');ax.legend(fontsize=8)
    save(fig,'F4_fixed_composition','Raw and fixed-composition estimands; missing support remains a gap')
    write_table(figdir/'F5_sensitivity.tsv',sensitivity)
    fig,axes=plt.subplots(1,2,figsize=(16,9))
    for i,ax in enumerate(axes,1):
        rows=[r for r in sensitivity if r.get('measure')=='technical_subset_time_slope' and r.get('pc_or_k')==f'PC{i}']
        for j,r in enumerate(rows):
            if r.get('estimate') is not None:
                ax.plot(r['estimate'],j,'o')
                if r.get('ci_low') is not None:ax.plot([r['ci_low'],r['ci_high']],[j,j],'-')
            else:ax.text(0,j,'NA: '+str(r.get('reason')),fontsize=8)
        ax.axvline(0,color='.7');ax.set(yticks=range(len(rows)),yticklabels=[r['analysis_set'] for r in rows],xlabel=f'PC{i} per 1000 years toward present',title='Prespecified technical subsets')
    save(fig,'F5_technical_sensitivity','Conditional sensitivity; shotgun labels are not independent WGS validation')
    return [str(p.relative_to(out)) for p in sorted(figdir.glob('*.png'))]


@transactional('manifest.json')
def run(args):
    if args.min_samples<2 or args.bin_years<1 or args.max_age_bp<args.bin_years or args.bootstrap<2 or args.display_limit<1 or args.date_draws<1 or args.pcs is not None and args.pcs<2:
        raise ValueError('Require >=2 PCs/samples/bootstrap, positive bins, display limit and date draws')
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    samples,coords,columns,contract=load_inputs(args);k=len(columns)
    verified=contract['status']=='verified' and contract['metadata_provenance']['historical_metadata_verification']=='verified'
    for s in samples:
        s['basis_id']=contract.get('basis_id');s['coordinate_kind']=contract['coordinate_kind']
    counts=dict(total_coordinates=len(coords),total_metadata=contract['join']['n_metadata'],joined=contract['join']['n_joined'],
                displayed_cohorts=dict(Counter(s['cohort'] for s in samples if s['coordinate_present'])),formal=sum(s['formal'] for s in samples))
    plan=dict(status='planned',command=args.command,coordinate_contract=contract,sample_counts=counts,parameters=vars(args),
              bins='left-closed right-open; last ancient bin includes max_age_bp; modern separate',
              weights='fixed over all formal ancient samples in prespecified window; no missing-stratum extrapolation',
              inference='PC response linear HC3, region/source/technology and available coverage; all PC x region/global BH family',
              bootstrap='unrelated individual, conditional on fixed coordinates/composition',date_uncertainty='not_provided_unless_explicit_inputs',
              planned_outputs=['sample_audit.tsv.gz','sample_composition.tsv','support_audit.tsv','reference_weights.tsv','pc_scaling.tsv',
                               'temporal_centroids.tsv','temporal_dispersion.tsv','temporal_contrasts.tsv','standardized_trajectories.tsv',
                               'temporal_models.tsv','sensitivity_summary.tsv','figures','report.html'],code=code_identity())
    atomic_json(out/'trajectory_plan.json',plan)
    write_table(out/'sample_audit.tsv.gz',samples)
    comp=composition(samples,args.bin_years,args.max_age_bp);write_table(out/'sample_composition.tsv',comp)
    standardized,weights,stdsupport=standardized_tables(samples,k,args)
    write_table(out/'reference_weights.tsv',weights)
    if args.command=='trajectory-plan':
        write_table(out/'support_audit.tsv',stdsupport)
        return dict(plan,analysis_execution='not_run_plan_only',limitations=LIMITATIONS)
    scales,scalerows=scaling(samples,k);write_table(out/'pc_scaling.tsv',scalerows)
    centroid,dispersion,contrasts,support=trajectory_tables(samples,k,args,scales)
    fits=models(samples,k,verified,args);sens=sensitivities(samples,k,args,verified)
    date_rows,date_status=date_sensitivity(samples,args);write_table(out/'date_sensitivity.tsv.gz',date_rows)
    write_table(out/'sample_audit.tsv.gz',samples)
    for name,records in [('temporal_centroids',centroid),('temporal_dispersion',dispersion),('temporal_contrasts',contrasts),
                          ('standardized_trajectories',standardized),('support_audit',support+stdsupport),('temporal_models',fits),('sensitivity_summary',sens)]:
        for r in records:r.update(basis_id=contract.get('basis_id'),coordinate_kind=contract['coordinate_kind'])
        write_table(out/(name+'.tsv'),records)
    pictures=[] if args.no_plots else figures(samples,centroid,dispersion,standardized,sens,comp,out,args)
    text='<html lang="en"><meta charset="utf-8"><title>Fixed PCA temporal structure</title><body style="font:16px system-ui;max-width:1400px;margin:32px auto"><h1>Fixed PCA temporal structure</h1>'
    text+='<p>Raw, fixed-composition, and adjusted PC-response estimands are separate. HERV conditional associations require the independent time stage. Modern is a separate endpoint.</p>'
    text+='<pre>'+html.escape(json.dumps(counts,indent=2))+'</pre>'
    for picture in pictures:text+='<img style="width:100%" src="'+picture+'">'
    text+='<ul>'+''.join('<li>'+html.escape(s)+'</li>' for s in LIMITATIONS)+'</ul>'
    text+='<p>Provenance: '+html.escape('verified' if verified else 'unverified; formal adjusted models unavailable')+'</p><a href="manifest.json">Manifest and output inventory</a></body></html>'
    (out/'report.html').write_text(text)
    return dict(status='complete' if verified else 'complete_descriptive_unverified',schema='pca-trajectory-1',sample_counts=counts,
                coordinate_contract=contract,parameters=vars(args),date_uncertainty=date_status,limitations=LIMITATIONS,
                hypothesis_family='all_prespecified_PC_by_region_and_global_linear_models',n_models=len(fits),
                n_estimable_models=sum(r['status']=='estimated' for r in fits),genotype_data_read=False,
                original_points_and_files_modified=False,unexecuted=['new_reference_fits','missing_marker_reprojection','dosage_validation','independent_reads'],
                figures=pictures,interpretation='coordinate temporal structure; no selection or causal HERV evidence')


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['trajectory-plan','trajectory'])
    p.add_argument('--coordinates','--eigenvec',dest='coordinates',required=True)
    p.add_argument('--metadata',required=True)
    p.add_argument('--out',default='/mnt/d/analysis/allele3/pca3d/temporal')
    for name in ('basis-dir','coordinate-receipt','provenance-manifest','training-ids','supplement-metadata','date-distributions'):p.add_argument('--'+name)
    p.add_argument('--pcs',type=int,default=None,help='Default PC1/2 primary plus available 5/10 supplements; explicit K requires every PC')
    for name,value in [('bin-years',1000),('max-age-bp',10000),('min-samples',20),('bootstrap',500),('seed',20261028),('display-limit',20000),('date-draws',100)]:p.add_argument('--'+name,type=int,default=value)
    p.add_argument('--date-interval-uniform',action='store_true',help='Explicit assumption using real age_lower_bp/age_upper_bp fields')
    p.add_argument('--no-plots',action='store_true')
    return p


if __name__=='__main__':
    try:
        result=run(parser().parse_args());print(json.dumps({k:result.get(k) for k in ('status','sample_counts','n_models','n_estimable_models')}))
    except (ValueError,OSError,RuntimeError) as exc:raise SystemExit('pca-trajectory: '+str(exc))
