"""Bounded same-person matched HERV statistics on a frozen measurable panel.

Pass 1 freezes pair membership using callability alone. Pass 2 reads one pair
at a time and accumulates sums/counts. No sample-by-catalog matrix is retained.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import base64
import json
import math
import os
from pathlib import Path
import numpy as np
from analysis_integrity import atomic_json, identity
from allele_measurement import open_pgen_with_pvar, panel_pvar, read_counted_dosage
from temporal_stats import number, technology, linear_model, adjust_fdr, fixed_composition, write_table, read_table, seed_for


def pair_candidates(targets, selected, fixed_flags):
    import herv_allele_time as hat
    expected=defaultdict(Counter); by_target={t['site_id']:t for t in targets}; found=defaultdict(lambda:{'case':[],'control':[]})
    for target in targets:
        for base in target['alts'].split(','):expected[target['match_id'],hat.alt_role(target,base)][target['group']]+=1
    for entry,keep in zip(selected,fixed_flags):
        if keep and entry['strict_filter_complete']:found[entry['match_id'],entry['role']][entry['group']].append(entry)
    pairs=[];audit=[]
    for key,count in sorted(expected.items()):
        cases,controls=found[key]['case'],found[key]['control']
        good=len(cases)==count['case']==1 and len(controls)==count['control'] and bool(controls)
        good=good and all(e['chrom']==cases[0]['chrom'] for e in controls)
        audit.append(dict(match_id=key[0],role=key[1],expected_cases=count['case'],expected_controls=count['control'],
                          measured_cases=len(cases),measured_controls=len(controls),status='candidate' if good else 'not_estimable',
                          reason='' if good else 'missing_nonunique_or_cross_chromosome_exact_pair'))
        if not good:continue
        case=cases[0];target=by_target[case['site_id']]
        memberships=[v for v in target.get('candidate_families',case['family']).split(',') if v not in {'','.'}]
        scope='boundary_only_sensitivity' if case['partition']=='boundary_only' else 'nonoverlap_main' if len(memberships)<=1 and str(target.get('spatial_overlap','0')) in {'0','.','False','false','none','','single_geometry'} and case['partition']!='ambiguous_overlap' else 'overlap_sensitivity'
        pairs.append(dict(pair_id=key[0]+'|'+key[1],case=case,controls=controls,scope=scope,
                          group=(case['family'],case['partition'],case['role'],scope),membership=memberships))
    return pairs,audit


def iter_pair_values(prefix, n, pairs, *, audit=None):
    """Use the same nucleotide count as individual and family trajectories."""
    with open_pgen_with_pvar(str(prefix)+'.pgen', panel_pvar(prefix), n) as reader:
        def read(entry):
            counted = read_counted_dosage(reader, entry['measurement_plan'], n, audit)
            value = counted.astype(float)/2
            value[counted < 0] = np.nan
            return value
        for pair in pairs:
            if not pair['controls']:
                raise ValueError('A matched pair requires at least one control')
            case = read(pair['case'])
            total = np.zeros(n)
            valid = np.isfinite(case)
            for control in pair['controls']:
                value = read(control)
                valid &= np.isfinite(value)
                total += np.nan_to_num(value)
            controls = total/len(pair['controls'])
            case[~valid] = np.nan
            controls[~valid] = np.nan
            yield pair, case, controls, valid


def sample_cells(samples,args):
    cells=defaultdict(list);formal=[]
    for i,s in enumerate(samples):
        row=dict(s,formal=s.get('include')=='1' and s.get('unrelated')=='1' and s['region']!='unknown',
                 measurement_source=technology(s)+'|'+s.get('data_source','unknown'))
        formal.append(row)
        if not row['formal']:continue
        if s['age_scope']=='primary_ancient':
            cells['all',s['time_bin']].append(i);cells[s['region'],s['time_bin']].append(i)
        elif s['cohort']=='modern':cells[s['region'],'modern'].append(i)
    return formal,{key:np.array(v,dtype=np.int32) for key,v in cells.items()}


def finalize_panel_cells(blocks, panel_counts, cells, args):
    """Each cell must contain every frozen pair; no per-cell panel renormalization."""
    rows=[]
    for group,n_pairs in sorted(panel_counts.items()):
        for cell,indices in sorted(cells.items(),key=lambda x:str(x[0])):
            data=blocks.get((group,cell),{})
            accepted=sum(v[3] for v in data.values())
            ok=accepted==n_pairs and len(indices)>=args.min_samples
            total=np.sum(list(data.values()),axis=0) if data else np.zeros(4)
            row=dict(family=group[0],partition=group[1],role=group[2],membership_scope=group[3],region=cell[0],time_bin=cell[1],
                cohort='modern' if cell[1]=='modern' else 'ancient',n=len(indices),n_fixed_pairs=n_pairs,n_estimable_pairs=int(accepted),
                estimate=float(total[0]/n_pairs) if ok else None,case_mean=float(total[1]/n_pairs) if ok else None,
                control_mean=float(total[2]/n_pairs) if ok else None,status='estimated' if ok else 'not_estimable',
                reason='' if ok else 'some_frozen_pairs_lack_cell_support',measure='fixed_measurable_panel_same_person_pair_difference',
                weights_id='equal_frozen_pair',support_id='|'.join(group),ci_low=None,ci_high=None,
                uncertainty='whole_chromosome_bootstrap_conditional_on_people',n_genomic_blocks=len(data),seed=args.seed)
            if ok and len(data)>=3 and args.bootstrap>=2:
                a=np.array(list(data.values()));rng=np.random.default_rng(seed_for(args.seed,(group,cell)));draws=[]
                for _ in range(args.bootstrap):
                    x=a[rng.integers(len(a),size=len(a))].sum(axis=0);draws.append(x[0]/x[3])
                row.update(ci_low=float(np.quantile(draws,.025)),ci_high=float(np.quantile(draws,.975)))
            rows.append(row)
    return rows


def marker_overlap(args,targets,out):
    basis=getattr(args,'basis_dir',None) or (str(Path(args.coordinates).parent) if args.coordinates else None)
    marker=Path(basis)/'all_pruned.pvar' if basis else None
    if marker is None or not marker.is_file():return {'status':'not_available','needed':'original basis all_pruned.pvar'}
    import herv_time as ht
    positions={(str(t.get('chrom37')),int(t['pos37'])) for t in targets if t.get('mapping_status')=='mapped' and t.get('pos37') is not None}
    counts=Counter();overlap=[]
    for row in ht.rows(marker,('CHROM','POS','ID','REF','ALT'),whitespace=True):
        chrom=row['CHROM'].removeprefix('chr');counts[chrom]+=1
        if (chrom,int(row['POS'])) in positions:overlap.append(row)
    write_table(out/'pca_marker_target_overlap.tsv',overlap)
    return dict(status='audited',marker_source=identity(marker,full=True),n_original_markers=sum(counts.values()),
                markers_by_chromosome=dict(counts),n_exact_coordinate_overlap=len(overlap),
                exclusions_applied=False,flanks_added=False,LOCO='not_run_requires_separate_basis_design')


def run_joint(args,targets,selected,samples,coord_contract,coordinates,fixed_flags):
    import herv_time as ht
    import herv_allele_time as hat
    out=Path(args.out);n=len(samples)
    people,cells=sample_cells(samples,args)
    pairs,pair_audit=pair_candidates(targets,selected,fixed_flags)
    epochs=sorted({s['time_bin'] for s in samples if s['age_scope']=='primary_ancient' and s.get('include')=='1' and s.get('unrelated')=='1'})
    qualified={ep:cells.get(('all',ep),np.array([],dtype=np.int32)) for ep in epochs}
    qualified={ep:idx for ep,idx in qualified.items() if len(idx)>=args.min_samples}
    frozen=[];freeze_audit=[]
    measurement_qc=Counter()
    # Measurability-only selection. No age slope, effect sign, peak or p value is examined.
    for pair,case,control,valid in iter_pair_values(args.pgen_prefix,n,pairs,audit=measurement_qc) if pairs else []:
        rates={ep:float(valid[idx].mean()) for ep,idx in qualified.items()}
        keep=len(rates)>=2 and all(rate>=args.min_call_rate for rate in rates.values())
        freeze_audit.append(dict(pair_id=pair['pair_id'],family=pair['group'][0],partition=pair['group'][1],role=pair['group'][2],
            membership_scope=pair['scope'],candidate_families=','.join(pair['membership']),
            epoch_call_rates=json.dumps(rates,sort_keys=True),frozen=int(keep),
            case_contrast=pair['case']['contrast_id'],control_contrasts=';'.join(c['contrast_id'] for c in pair['controls']),
            status='fixed' if keep else 'not_estimable',reason='' if keep else 'insufficient_same_person_cross_epoch_callability'))
        if keep:frozen.append(pair)
    panel_counts=Counter(p['group'] for p in frozen)
    for row in freeze_audit:
        key=(row['family'],row['partition'],row['role'],row['membership_scope'])
        row['fixed_pair_weight']=1/panel_counts[key] if row['frozen'] else None
    write_table(out/'joint_pair_candidate_audit.tsv',pair_audit)
    write_table(out/'joint_fixed_pair_panel.tsv',freeze_audit)
    atomic_json(out/'joint_panel_plan.json',dict(status='planned',n_fixed_pairs=len(frozen),
        selection='HQ strict counted-nucleotide measurements; same-person pair callability >= threshold in every qualified pooled ancient epoch; frozen before effect estimates',
        qualified_epochs={str(k):len(v) for k,v in qualified.items()},min_call_rate=args.min_call_rate,
        weights='equal fixed pairs; controls equal within each fixed match',
        formal_estimand='pair-by-region-by-time means, then fixed pair weights; any missing pair makes that cell unavailable',
        personal_estimand='observed_pair_conditional unless every frozen pair is measured; formal personal model uses complete panel only'))
    totals={g:np.zeros((n,3)) for g in panel_counts} # case sum, control sum, paired count
    chromosome_totals=defaultdict(dict)
    blocks=defaultdict(lambda:defaultdict(lambda:np.zeros(4)))
    cell_fields=['pair_id','family','partition','role','membership_scope','region','time_bin','cohort','n','n_eligible','call_rate','case_mean','control_mean','difference','status']
    with hat.table(out/'joint_pair_cell_estimates.tsv.gz',cell_fields) as emit, hat.table(out/'joint_pair_call_masks.tsv.gz',['pair_id','PSAM_order_bitmask_base64','n_called','n_PSAM']) as mask_emit:
        for pair,case,control,valid in iter_pair_values(args.pgen_prefix,n,frozen,audit=measurement_qc) if frozen else []:
            g=pair['group'];chrom=pair['case']['chrom'];data=totals[g]
            data[:,0]+=np.nan_to_num(case);data[:,1]+=np.nan_to_num(control);data[:,2]+=valid
            if chrom not in chromosome_totals[g]:chromosome_totals[g][chrom]=np.zeros((n,3))
            chromosome_totals[g][chrom]+=np.column_stack((np.nan_to_num(case),np.nan_to_num(control),valid))
            mask_emit(dict(pair_id=pair['pair_id'],PSAM_order_bitmask_base64=base64.b64encode(np.packbits(valid).tobytes()).decode(),n_called=int(valid.sum()),n_PSAM=n))
            for cell,idx in sorted(cells.items(),key=lambda x:str(x[0])):
                called=idx[valid[idx]];good=len(called)>=args.min_samples and len(called)/len(idx)>=args.min_call_rate
                c=float(case[called].mean()) if good else None;z=float(control[called].mean()) if good else None
                emit(dict(pair_id=pair['pair_id'],family=g[0],partition=g[1],role=g[2],membership_scope=g[3],region=cell[0],time_bin=cell[1],
                    cohort='modern' if cell[1]=='modern' else 'ancient',n=len(called),n_eligible=len(idx),call_rate=len(called)/len(idx),case_mean=c,control_mean=z,
                    difference=c-z if good else None,status='estimated' if good else 'not_estimable'))
                if good:blocks[g,cell][chrom]+=[c-z,c,z,1]
    cell_rows=finalize_panel_cells(blocks,panel_counts,cells,args)
    write_table(out/'joint_fixed_panel_trajectories.tsv',cell_rows)
    # Fixed regional composition, with identical pair sets in every region/epoch.
    region_counts=Counter(s['region'] for s in people if s['formal'] and s['age_scope']=='primary_ancient')
    region_weights={r:n/sum(region_counts.values()) for r,n in sorted(region_counts.items())} if region_counts else {}
    write_table(out/'joint_region_weights.tsv',[dict(region=r,weight=w,weights_id='formal_ancient_region_composition') for r,w in region_weights.items()])
    standardized=[]
    for g in panel_counts:
        for ep in epochs:
            found={r['region']:r for r in cell_rows if (r['family'],r['partition'],r['role'],r['membership_scope'])==g and r['time_bin']==ep}
            missing=[r for r in region_weights if r not in found or found[r]['estimate'] is None]
            good=bool(region_weights) and not missing
            standardized.append(dict(family=g[0],partition=g[1],role=g[2],membership_scope=g[3],time_bin=ep,
                measure='fixed_region_and_pair_composition',n_fixed_pairs=panel_counts[g],
                estimate=sum(w*found[r]['estimate'] for r,w in region_weights.items()) if good else None,
                status='estimated' if good else 'not_estimable',reason='' if good else 'missing_region_pair_support:'+';'.join(missing),
                weights_id='formal_ancient_region_composition',support_id='|'.join(g)))
    write_table(out/'joint_standardized_trajectories.tsv',standardized)
    fits=[];missingness=[]
    verified=coord_contract.get('status')=='verified' and coord_contract.get('metadata_provenance',{}).get('historical_metadata_verification')=='verified'
    explicit_unadjusted=args.pcs==0
    personal_fields=['sample_id','cohort','age_bp','region','family','partition','role','membership_scope','n_fixed_pairs','n_valid_pairs','case_mean','control_mean','paired_difference','callability','estimand','formal_model_eligible','support_id']
    with hat.table(out/'joint_person_statistics.tsv.gz',personal_fields) as emit:
        for g,data in sorted(totals.items()):
            P=panel_counts[g];counts=data[:,2];complete=counts==P
            difference=np.divide(data[:,0]-data[:,1],counts,out=np.full(n,np.nan),where=counts>0)
            strict=np.where(complete,difference,np.nan)
            for i,s in enumerate(samples):
                emit(dict(sample_id=s['sample_id'],cohort=s['cohort'],age_bp=s['age_bp'],region=s['region'],family=g[0],partition=g[1],role=g[2],membership_scope=g[3],
                    n_fixed_pairs=P,n_valid_pairs=int(counts[i]),case_mean=float(data[i,0]/counts[i]) if counts[i] else None,
                    control_mean=float(data[i,1]/counts[i]) if counts[i] else None,paired_difference=float(difference[i]) if counts[i] else None,callability=float(counts[i]/P),
                    estimand='strict_complete_fixed_panel' if complete[i] else 'observed_pair_conditional',
                    formal_model_eligible=int(complete[i] and people[i]['formal'] and s['age_scope']=='primary_ancient'),support_id='|'.join(g)))
            for region in ['all',*sorted(region_counts)]:
                members=[i for i,s in enumerate(people) if s['age_scope']=='primary_ancient' and (region=='all' or s['region']==region)]
                ss=[people[i] for i in members]
                fit=linear_model(ss,strict[members],pcs=coordinates if args.pcs else None,minimum=args.min_samples) if verified or explicit_unadjusted else dict(status='not_estimable',reason='genetic_adjustment_unavailable',n=0,estimate=None,p=None)
                fit.update(family=g[0],partition=g[1],role=g[2],membership_scope=g[3],region=region,n_fixed_pairs=P,
                    measure='strict_complete_fixed_panel_pair_time_association',estimand='complete_panel_people_only; transport_to_incomplete_people_not_assumed',
                    genetic_adjustment='verified_shared_PCs' if verified and args.pcs else 'none_explicit_pcs_0' if explicit_unadjusted else 'unavailable',
                    inference_scope='person_HC3_conditional_on_fixed_loci; chromosome_sensitivity_reported_separately')
                # Whole chromosomes keep all related ALTs and both matched arms together.
                slopes=[]
                if fit.get('estimate') is not None and len(chromosome_totals[g])>=3:
                    for chromosome,removed in chromosome_totals[g].items():
                        n_removed=sum(p['group']==g and p['case']['chrom']==chromosome for p in frozen)
                        denom=P-n_removed
                        if denom<=0:continue
                        response=np.where(complete,(data[:,0]-data[:,1]-removed[:,0]+removed[:,1])/denom,np.nan)
                        loo=linear_model(ss,response[members],pcs=coordinates if args.pcs else None,minimum=args.min_samples)
                        if loo.get('estimate') is not None:slopes.append(loo['estimate'])
                fit.update(n_chromosome_sensitivity=len(slopes),chromosome_jackknife_se=float(np.std(slopes,ddof=0)*math.sqrt(len(slopes)-1)) if len(slopes)>=3 else None,
                           chromosome_sensitivity_status='estimated' if len(slopes)>=3 else 'insufficient_genomic_blocks')
                fits.append(fit)
            members=[i for i,s in enumerate(people) if s['age_scope']=='primary_ancient'];ss=[people[i] for i in members]
            fit=linear_model(ss,counts[members]/P,pcs=coordinates if args.pcs else None,minimum=args.min_samples) if verified or explicit_unadjusted else dict(status='not_estimable',reason='genetic_adjustment_unavailable')
            missingness.append(dict(fit,family=g[0],partition=g[1],role=g[2],membership_scope=g[3],measure='paired_callability_time_association'))
    adjust_fdr(fits);write_table(out/'joint_person_models.tsv',fits);write_table(out/'joint_missingness_models.tsv',missingness)
    funnel=measurability_funnel(out,selected,fixed_flags,freeze_audit,fits,args)
    overlap=marker_overlap(args,targets,out)
    figure_files=[] if args.no_plots else plot_joint(funnel,cell_rows,standardized,out,args)
    result=dict(status='complete' if frozen else 'no_usable_contrast',n_candidate_pairs=len(pairs),n_fixed_pairs=len(frozen),n_PSAM=n,
        per_person_storage='O(samples * family_partition_role * chromosomes); never samples * loci',
        formal_person_model='strict complete fixed panel only; incomplete personal means descriptive observed_pair_conditional',
        formal_cell_model='all frozen pairs with original fixed weights; any unsupported pair produces NA',
        dependencies='person HC3 conditional on fixed loci; genomic chromosome bootstrap and leave-chromosome-out sensitivity separately',
        extra_ALT='unencoded is NA; actual native allele indices are retained; unknown split REF is NA, not an ALT complement',
        measurement_schema='nucleotide-count-1',measurement_qc=dict(measurement_qc),
        measurement_qc_scope='read operations across freeze and estimation passes, not unique biological observations',
        spatial_membership='nonoverlapping main; ambiguous/multiple memberships separately labelled sensitivity; seeded assignment is not unique biology',
        basis_id=coord_contract.get('basis_id'),coordinate_kind=coord_contract.get('coordinate_kind'),
        pca_background_audit=overlap,figures=figure_files,interpretation='fixed measurable allele panel; not ancient original gnomAD O/E or selection coefficients')
    atomic_json(out/'joint_manifest.json',result)
    return result


def measurability_funnel(out,selected,fixed_flags,panel,fits,args):
    counts=defaultdict(lambda:defaultdict(set))
    for row in read_table(out/'all_allele_audit.tsv.gz'):
        if row['role']=='REF':continue
        g=(row['family'],row['partition'],row['role']);allele=row['site_id']+'|'+row['allele']
        counts[g]['catalog_ALT'].add(allele)
        if row['mapping_status']=='mapped':counts[g]['mapped_GRCh37_ALT'].add(allele)
        if row['panel_status'] in {'assayed_pass','not_in_high_quality','pvar_filter_failed'} or row.get('pvar_index') not in {'NA','.','',None}:
            counts[g]['exact_PVAR_REF_ALT'].add(allele)
        if row['panel_status']=='assayed_pass':counts[g]['HQ_exact_ALT'].add(allele)
    for entry,keep in zip(selected,fixed_flags):
        if keep:counts[entry['family'],entry['partition'],entry['role']]['fixed_epoch_callable_ALT'].add(entry['contrast_id'])
    stages_order=['catalog_ALT','mapped_GRCh37_ALT','exact_PVAR_REF_ALT','HQ_exact_ALT','fixed_epoch_callable_ALT','modeled_ALT_contrasts']
    for row in read_table(out/'adjusted_allele_trends.tsv.gz'):
        if row.get('status') in {'descriptive_adjusted_imputed_association','descriptive_no_genetic_adjustment'}:
            counts[row['family'],row['partition'],row['role']]['modeled_ALT_contrasts'].add(row['site_id']+'|'+row['counted_base38'])
    rows=[dict(family=g[0],partition=g[1],role=g[2],stage=stage,n=len(stages[stage]),unit='ALT_contrasts') for g,stages in sorted(counts.items()) for stage in stages_order]
    epochs=defaultdict(set)
    for row in read_table(out/'allele_trajectories.tsv.gz'):
        if row['role']=='REF' or row['region']!='all' or row.get('measurement_source','.') not in {'.','all','NA',''}: continue
        if (number(row.get('n_called')) or 0) >= args.min_samples and (number(row.get('call_rate')) or 0) >= args.min_call_rate:
            epochs[row['family'],row['partition'],row['role'],row['cohort'],row['time_bin']].add(row['contrast_id'])
    rows += [dict(family=g[0],partition=g[1],role=g[2],cohort=g[3],time_bin=g[4],stage='epoch_callable_ALT',n=len(values),unit='ALT_contrasts') for g,values in sorted(epochs.items())]
    for row in panel:
        if row['frozen']:rows.append(dict(family=row['family'],partition=row['partition'],role=row['role'],stage='fixed_same_person_pair',n=1,unit='matched_pair',membership_scope=row['membership_scope']))
    for row in fits:
        rows.append(dict(family=row['family'],partition=row['partition'],role=row['role'],stage='strict_person_model_complete_cases',n=row.get('n',0),unit='people',region=row['region'],membership_scope=row['membership_scope'],status=row['status']))
    write_table(out/'joint_measurability_funnel.tsv',rows)
    return rows


def plot_joint(funnel,cells,standardized,out,args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(16,9))
    stages=['catalog_ALT','mapped_GRCh37_ALT','exact_PVAR_REF_ALT','HQ_exact_ALT','fixed_epoch_callable_ALT']
    for family in sorted({r['family'] for r in funnel}):
        axes[0].plot(range(len(stages)),[sum(r['n'] for r in funnel if r['family']==family and r['stage']==stage) for stage in stages],marker='o',label=family)
    axes[0].set(xticks=range(len(stages)),xticklabels=stages,ylabel='ALT contrasts',title='Measurability funnel');axes[0].tick_params(axis='x',rotation=25);axes[0].legend()
    for family in sorted({r['family'] for r in standardized}):
        for partition in sorted({r['partition'] for r in standardized if r['family']==family}):
            rr=sorted([r for r in standardized if r['family']==family and r['partition']==partition and r['role']=='primary' and r['membership_scope']=='nonoverlap_main'],key=lambda r:r['time_bin'])
            axes[1].plot([(r['time_bin']+.5)*args.bin_years for r in rr],[r['estimate'] if r['estimate'] is not None else np.nan for r in rr],marker='o',label=family+':'+partition)
    axes[1].set(xlabel='Age BP',ylabel='Same-person case minus matched control allele fraction',title='Fixed regional and pair composition; primary ALT');axes[1].legend(fontsize=8)
    fig.suptitle('HERV conditional measurable-panel statistics; not ancient gnomAD enrichment or selection')
    fig.text(.01,.02,'Machine tables: joint_measurability_funnel.tsv / joint_standardized_trajectories.tsv. N, exclusions and provenance: manifest.json.',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,.94))
    for ext in ('png','pdf'):fig.savefig(out/('F6_joint_measurability.'+ext),dpi=180)
    plt.close(fig)
    return ['F6_joint_measurability.png','F6_joint_measurability.pdf']


def family_composition(family_rows,samples,args):
    """Additional fixed-composition versions of the same full ALT family panel."""
    primary=[s for s in samples if s['age_scope']=='primary_ancient']
    definitions={
        'fixed_region':Counter((s['region'],'all') for s in primary),
        'fixed_region_technology_source':Counter((s['region'],technology(s)+'|'+s.get('data_source','unknown')) for s in primary),
    }
    output=[];weights=[]
    keys=sorted({(r['family'],r['group'],r['partition'],r['role']) for r in family_rows})
    lookup={(r['family'],r['group'],r['partition'],r['role'],r['region'],r.get('measurement_source','all'),r['time_bin']):r for r in family_rows if r['age_scope']=='primary_ancient'}
    epochs=range(math.ceil(args.max_age_bp/args.bin_years))
    for method,counts in definitions.items():
        if not counts:continue
        total=sum(counts.values());ws={key:n/total for key,n in counts.items()}
        weights.extend(dict(method=method,region=key[0],measurement_source=key[1],weight=w,n_reference=counts[key]) for key,w in ws.items())
        for key in keys:
            for ep in epochs:
                strata={cell:lookup.get((*key,*cell,ep)) for cell in ws}
                missing=[cell for cell,r in strata.items() if not r or r['n_metadata']<args.min_samples or r['n_fixed_contrasts']!=r['expected_fixed_contrasts']]
                valid=not missing
                output.append(dict(family=key[0],group=key[1],partition=key[2],role=key[3],time_bin=ep,cohort='ancient',measure=method,
                    estimate=sum(ws[cell]*r['mean_AF'] for cell,r in strata.items()) if valid else None,
                    n=sum(r['n_metadata'] for r in strata.values() if r),n_eff=1/sum(ws[cell]**2/r['n_metadata'] for cell,r in strata.items()) if valid else None,
                    n_fixed_contrasts=next(iter(strata.values()))['expected_fixed_contrasts'] if valid else None,
                    status='descriptive_fixed_composition' if valid else 'not_estimable',
                    reason='' if valid else 'missing_stratum_or_contrast_support:'+json.dumps(missing),
                    uncertainty='not_PC_adjusted; no_person_or_locus_independence_assumed',weights_id=method))
    write_table(Path(args.out)/'family_composition_weights.tsv',weights)
    write_table(Path(args.out)/'standardized_family_trajectories.tsv',output)
    return output
