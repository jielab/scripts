#!/usr/bin/env python3
"""GRCh38 enrichment and within-site AVI tests using observed target summaries.

enrich_test uses a two-pass, bounded-memory LOCO background conditional on
observed eligible variants. allele_test compares reference-relative impact
scores of nonreference substitutions; no AVI value is invented for REF.
"""
from __future__ import annotations
import argparse, csv, gzip, heapq, itertools, json, math, os, bisect
from collections import Counter, defaultdict
from pathlib import Path

import sys as _path_sys
_path_sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "0f"))
from workflow_paths import configured_path
import numpy as np
from research_common import (Intervals,atomic_json,bed_rows,bh,chrom_name,chrom_key,
    chromosomes,fmt,html_report,html_table,identity,number,open_text,stable_hash,write_tsv)
from avi_scan import resolve_vcf,query_file
from allele_geometry import GEOMETRY_FIELDS,analyze_site

OUTCOMES=('multi_ALT','two_ALT','three_ALT','singleton_ALT','observed_biallelic','observed_triallelic','observed_quadallelic')
METRICS=('mean_avi','sd_site_mean','mean_alt_range','mean_abs_gradient','high_fraction','low_fraction')


def source_rows(a):
    from summary_scan import iter_sites
    return iter_sites(a.scan_dir,chromosomes=chromosomes(a.chr),population=a.population)


def frequency_bin(value):
    return 'NA' if value is None else str(bisect.bisect_right([1e-4,1e-3,.01,.05,.5],value))


def denominator_bin(value):
    return 'NA' if value is None else str(bisect.bisect_right([100,1000,10000,100000],value))


def observed_outcomes(row,a):
    """Each tuple is (eligible denominator, binary outcome), never imputed REF."""
    afs=row['alt_afs'];acs=row['alt_acs'];ans=row['alt_ans']
    count_mode=row['count_mode']!='AF_only'
    if a.threshold_kind=='AC' and not count_mode:
        raise ValueError('--threshold-kind AC requires reported AC counts; this population is AF-only. Use --threshold-kind auto or AF.')
    observed=[i for i in range(len(afs)) if acs[i] is not None and acs[i]>0] if count_mode else [i for i in range(len(afs)) if afs[i] is not None and afs[i]>0]
    if a.threshold_kind=='AF' or (a.threshold_kind=='auto' and not count_mode):
        qualified=[i for i in observed if afs[i] is not None and afs[i]>=a.min_af and afs[i]>0]
    else:
        qualified=[i for i in observed if acs[i] is not None and acs[i]>=a.min_ac and ans[i] is not None and ans[i]>=a.min_an]
    evidence_complete=all(x is not None for x in (acs if count_mode else afs))
    if a.threshold_kind=='AF':evidence_complete=evidence_complete and all(x is not None for x in afs)
    if count_mode and a.threshold_kind!='AF':evidence_complete=evidence_complete and all(ans[i] is not None for i in observed)
    bases=row.get('observed_bases');known=len(bases) if isinstance(bases,list) else number(bases)
    strict=known is not None and known>=2
    return {
      'multi_ALT':(bool(qualified) and evidence_complete,len(qualified)>=2 and evidence_complete),
      'two_ALT':(bool(qualified) and evidence_complete,len(qualified)==2 and evidence_complete),
      'three_ALT':(bool(qualified) and evidence_complete,len(qualified)==3 and evidence_complete),
      'singleton_ALT':(count_mode and bool(observed) and all(x is not None for x in acs) and row.get('singleton_status') in ('present','absent'),row.get('singleton_status')=='present'),
      'observed_biallelic':(strict,known==2),
      'observed_triallelic':(strict,known==3),
      'observed_quadallelic':(strict,known==4)},qualified


class Context:
    def __init__(self,path):
        self.fasta=None
        if path:
            from scan import FastaContext
            self.fasta=FastaContext(path)
    def key(self,row):
        seq='NA'
        if self.fasta:
            seq=self.fasta.context(row['chrom'],int(row['pos']))
            if seq[1]!=row['ref']:raise ValueError('REF/FASTA mismatch at '+str(row['chrom'])+':'+str(row['pos']))
            if row['ref'] in 'AG':seq=seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
        eligible=[(af,i) for i,af in enumerate(row['alt_afs']) if af is not None and af>0 and (row['count_mode']=='AF_only' or (row['alt_acs'][i] is not None and row['alt_acs'][i]>0))]
        if eligible:
            af,i=max(eligible);an=row['alt_ans'][i]
        else:af=an=None
        return (seq,frequency_bin(af),denominator_bin(an))
    def close(self):
        if self.fasta:self.fasta.close()


class WindowLandscape:
    """Load one chromosome's compact AVI windows at a time, not basewise scores."""
    def __init__(self,root):
        self.root=Path(root);manifest=self.root/'qc/manifest.json'
        if not manifest.is_file():raise ValueError('Run 4.avi_scan.sh first: '+str(manifest))
        self.manifest=json.loads(manifest.read_text())
        if self.manifest.get('status')!='complete':raise ValueError('AVI scan is not complete')
        self.chrom=None;self.data={};self.cutoffs={}
        with open_text(self.root/'results/cutoffs.tsv') as f:
            for r in csv.DictReader(f,delimiter='\t'):self.cutoffs[(int(r['window_bp']),r['metric'])]=r
    def load(self,chrom):
        chrom=chrom_name(chrom)
        if chrom==self.chrom:return
        self.data=defaultdict(dict);self._profile_index=None;p=self.root/'shards'/f'{chrom}.windows.tsv.gz'
        if p.is_file():
            with open_text(p) as f:
                for r in csv.DictReader(f,delimiter='\t'):
                    w=int(r['window_bp']);start=int(r['start']);n=int(r['n_positions'])
                    self.data[w][start//w]=(int(r['end']),n,int(r['n_complete_positions']),float(r['coverage_fraction']),*[number(r[m]) for m in METRICS])
        self.chrom=chrom
    def eligible(self,row):
        if row is None:return False
        e=self.manifest['eligibility']
        return row[1]>=e['min_positions'] and row[3]>=e['min_coverage'] and row[2]/row[1]>=e['min_complete_fraction']
    def labels(self,chrom,pos,window):
        self.load(chrom);row=self.data.get(window,{}).get((pos-1)//window)
        if not self.eligible(row):return None
        result=[]
        for i,m in enumerate(METRICS):
            cut=self.cutoffs.get((window,m));v=row[4+i]
            if not cut or cut['status']!='ready' or v is None:continue
            lo,hi=float(cut['low_cut']),float(cut['high_cut'])
            if lo==hi:continue
            if v<=lo:result.append('AVI_'+m+'_low')
            if v>=hi:result.append('AVI_'+m+'_high')
        return result
    def profile(self,chrom,pos,radii):
        """Coverage-weighted means of local tile metrics around a position.

        Prefix sums are built lazily for the smallest tile scale of the current
        chromosome. Each query excludes the entire focal tile. Effective bp is
        metric-specific: an unavailable gradient contributes no gradient weight.
        Partial-tile coverage assumes scored positions are uniform within that
        tile; it is not a population callability or exact sequence denominator.
        """
        self.load(chrom)
        if not self.data:return []
        width=min(self.data);center=pos-1;focal_start=(center//width)*width;focal_end=focal_start+width;rows=[]
        if getattr(self,'_profile_index',None) is None:
            tiles=[(k,r) for k,r in sorted(self.data[width].items()) if self.eligible(r) and r[0]>k*width]
            starts=np.asarray([k*width for k,r in tiles],dtype=np.int64)
            ends=np.asarray([r[0] for k,r in tiles],dtype=np.int64)
            n=np.asarray([r[1] for k,r in tiles],dtype=float)
            values=np.asarray([[np.nan if v is None else v for v in r[4:]] for k,r in tiles],dtype=float).reshape((-1,len(METRICS)))
            weights=np.isfinite(values)*n[:,None]
            sums=np.where(np.isfinite(values),values,0.)*weights
            prefix_weights=np.vstack((np.zeros((1,len(METRICS))),np.cumsum(weights,axis=0)))
            prefix_sums=np.vstack((np.zeros((1,len(METRICS))),np.cumsum(sums,axis=0)))
            self._profile_index=(starts,ends,prefix_sums,prefix_weights,np.r_[0.,np.cumsum(n)])
        starts,ends,prefix_sums,prefix_weights,prefix_covered=self._profile_index
        zero=np.zeros(len(METRICS))
        def integral(x):
            i=int(np.searchsorted(starts,x,side='right'))-1
            if i<0:return zero,zero,0.
            fraction=min(1.,max(0.,(x-int(starts[i]))/(int(ends[i])-int(starts[i]))))
            return (prefix_sums[i]+fraction*(prefix_sums[i+1]-prefix_sums[i]),
                    prefix_weights[i]+fraction*(prefix_weights[i+1]-prefix_weights[i]),
                    float(prefix_covered[i]+fraction*(prefix_covered[i+1]-prefix_covered[i])))
        for lo,hi in zip([0]+radii[:-1],radii):
            sums=np.zeros(len(METRICS));weights=np.zeros(len(METRICS));covered=0.;band_bp=0
            for start,end in ((max(0,center-hi),max(0,center-lo)),(center+lo+1,center+hi+1)):
                band_bp+=max(0,end-start)
                # These two pieces exclude the focal tile, including any part
                # of that tile on either side of the focal nucleotide.
                for left,right in ((start,min(end,focal_start)),(max(start,focal_end),end)):
                    if right<=left:continue
                    sn,wn,cn=integral(right);sp,wp,cp=integral(left)
                    sums+=sn-sp;weights+=wn-wp;covered+=cn-cp
            for j,m in enumerate(METRICS):
                weight=max(0.,float(weights[j]))
                rows.append(dict(inner_bp=lo,outer_bp=hi,metric=m,value=float(sums[j]/weight) if weight>1e-9 else None,
                    effective_bp=weight,covered_bp=max(0.,covered),band_bp=band_bp,tile_bp=width))
        return rows


def prepare_herv(a,out):
    if a.herv_bed:path=Path(a.herv_bed)
    else:
        scan=Path(a.herv_scan_dir);choices=[scan/'reference_herv.bed',scan/'stable_loci.bed']
        path=next((p for p in choices if p.is_file()),None)
        if path is None:
            if not a.erv_joined or not Path(a.erv_joined).is_file():raise ValueError('Provide --herv-bed, finish a HERV sample, or provide --erv-joined')
            from scan import parse_joined_rmsk
            rows=[]
            for e in parse_joined_rmsk(a.erv_joined):
                if e.rep_class=='LTR' and e.family in {'ERV1','ERVK','ERVL','ERVL-MaLR','MaLR'}:
                    family='ERVL-MaLR' if e.family=='MaLR' else e.family
                    rows.extend((chrom_name(e.chrom),lo,hi,family+'|'+e.source_id) for lo,hi in e.segments)
            path=out/'annotations/reference_herv.bed';write_tsv(path,['#chrom','start0','end0','family|id'],rows)
    return Intervals(bed_rows(path,family=True)),dict(path=str(path),definition='Family-specific GRCh38 aligned reference intervals; sample stability where the selected source supplies it')


def prepare_archaic(a,out):
    paths=sorted(Path(a.archaic_dir).glob('*.segments.tsv.gz'))
    if not paths:raise ValueError('No GU IBDmix *.segments.tsv.gz: '+a.archaic_dir)
    groups=defaultdict(list);n=0;references=set();pops=set();chosen=set(chromosomes(a.chr))
    for path in paths:
        with open_text(path) as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r.get('genome_build') not in ('b38','38','GRCh38','hg38'):raise ValueError('IBDmix segment is not GRCh38: '+str(path))
                chrom=chrom_name(r['chrom'])
                if chrom not in chosen:continue
                if a.archaic_population!='ALL' and r.get('super_pop')!=a.archaic_population:continue
                if a.archaic_reference!='all' and r['anc']!=a.archaic_reference:continue
                lo,hi=int(r['start']),int(r['end'])
                if lo<0 or hi<=lo:raise ValueError('Invalid GU BED0 segment')
                if number(r.get('slod')) is not None and float(r['slod'])<a.archaic_min_slod:continue
                if hi-lo<a.archaic_min_bp:continue
                groups[(chrom,r['anc'],r['ID'])].append((lo,hi));n+=1;references.add(r['anc']);pops.add(r.get('super_pop','NA'))
    events=defaultdict(lambda:defaultdict(int))
    for (chrom,anc,sample),segments in groups.items():
        union=Intervals((chrom,x,y,'x') for x,y in segments)
        for c,x,y,_ in union.rows():events[(c,anc)][x]+=1;events[(c,anc)][y]-=1
    rows=[];coverage=[]
    for (c,anc),ev in events.items():
        active=0;previous=None;bp=0
        for pos,change in sorted(ev.items()):
            if previous is not None and pos>previous and active>=a.archaic_min_carriers:
                rows.append((c,previous,pos,'archaic:'+anc));bp+=pos-previous
            active+=change;previous=pos
        coverage.append(dict(chrom=c,reference=anc,covered_bp=bp,min_distinct_carriers=a.archaic_min_carriers))
    index=Intervals(rows);path=out/'annotations/archaic.bed';write_tsv(path,['#chrom','start0','end0','label'],index.rows())
    write_tsv(out/'annotations/archaic_coverage.tsv',['chrom','reference','covered_bp','min_distinct_carriers'],coverage)
    return index,dict(files=[identity(p) for p in paths],input_rows=n,references=sorted(references),source_populations=sorted(pops),
        coordinate='GU IBDmix already exports BED0 half-open; no extra minus-one conversion',min_carriers=a.archaic_min_carriers,
        interpretation='Genomic overlap with 1KG IBDmix tracts; not proof that target alleles reside on the introgressed haplotype',output=str(path))


def region_provider(a,out):
    if a.ref=='alphagenome':
        landscape=WindowLandscape(a.avi_dir)
        if not any(w==a.window_bp for w,m in landscape.cutoffs):raise ValueError('Requested --window-bp missing from AVI scan')
        return lambda c,p:landscape.labels(c,p,a.window_bp),dict(source=str(a.avi_dir),window_bp=a.window_bp),landscape
    index,meta=prepare_herv(a,out) if a.ref=='herv' else prepare_archaic(a,out)
    return lambda c,p:index.labels(c,p-1),meta,None


def bootstrap_oe(blocks,replicates,seed):
    pairs=np.asarray(list(blocks.values()),float)
    if len(pairs)<2 or pairs[:,1].sum()<=0 or pairs[:,0].sum()<=0:return (None,None,None,len(pairs))
    point=pairs[:,0].sum()/pairs[:,1].sum();rng=np.random.default_rng(seed);boot=[]
    for _ in range(replicates):
        v=pairs[rng.integers(0,len(pairs),len(pairs))].sum(axis=0)
        if v[1]>0 and v[0]>0:boot.append(v[0]/v[1])
    if len(boot)<max(20,replicates//2):return (None,None,None,len(pairs))
    lo,hi=np.quantile(boot,[.025,.975]);delta=np.log(boot)-math.log(point)
    p=(1+np.count_nonzero(np.abs(delta)>=abs(math.log(point))))/(len(boot)+1)
    return float(lo),float(hi),float(p),len(pairs)


class Reservoir:
    def __init__(self,n):self.n=n;self.data=defaultdict(list);self.serial=0
    def add(self,key,rank,value):
        if self.n<=0:return
        h=self.data[key];self.serial+=1;item=(-rank,self.serial,value)
        if len(h)<self.n:heapq.heappush(h,item)
        elif rank<-h[0][0]:heapq.heapreplace(h,item)


def enrichment(a):
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    regions,ref_meta,landscape=region_provider(a,out);context=Context(None if a.no_context else a.fasta)
    from enrichment_regions import RegionHighlighter
    local_regions=RegionHighlighter(a,out,ref_meta,landscape)
    masks=Intervals(itertools.chain.from_iterable(bed_rows(p,'excluded') for p in a.exclude_bed))
    selected=a.outcomes.split(',')
    if set(selected)-set(OUTCOMES):raise ValueError('--outcomes must be selected from '+','.join(OUTCOMES))
    background_outcomes=list(dict.fromkeys([*selected,'multi_ALT']))
    sensitivities=['all','exclude_mask'] if a.exclude_bed else ['all']
    total=defaultdict(lambda:np.zeros(2));per_chrom=defaultdict(lambda:np.zeros(2));status=Counter({outcome+suffix:0 for outcome in selected for suffix in ('_eligible','_case')});examined=0
    def usable(row):
        labels=regions(row['chrom'],int(row['pos']))
        if labels is None:return None
        sens=['all']
        if a.exclude_bed and not masks.labels(row['chrom'],int(row['pos'])-1):sens.append('exclude_mask')
        return labels,sens
    from enrichment_parallel import workers_requested, run as parallel_enrichment
    workers=workers_requested(a)
    if workers>1:
        context.close()
        examined,status,observed,expected,denom,unsupported,block,cases,controls = parallel_enrichment(
            a,regions,masks,local_regions,background_outcomes,workers)
    else:
        print('[enrich] pass 1: observed-site background by chromosome/context/AF/AN',flush=True)
        for row in source_rows(a):
            examined+=1;state,qualified=observed_outcomes(row,a);r=usable(row)
            if r is None:status['outside_reference_scoring_support']+=1;continue
            labels,sens=r;key=context.key(row);c=chrom_name(row['chrom'])
            for outcome in background_outcomes:
                eligible,case=state[outcome]
                if not eligible:continue
                status[outcome+'_eligible']+=1;status[outcome+'_case']+=int(case)
                for s in sens:total[(outcome,s,key)]+=(1,int(case));per_chrom[(c,outcome,s,key)]+=(1,int(case))
        print('[enrich] pass 2: region observed/expected and genomic blocks',flush=True)
        observed=defaultdict(float);expected=defaultdict(float);denom=Counter();unsupported=Counter();block=defaultdict(lambda:defaultdict(lambda:np.zeros(2)))
        cases=Reservoir(a.profile_sites);controls=Reservoir(a.controls_per_stratum)
        for row in source_rows(a):
            state,qualified=observed_outcomes(row,a);r=usable(row)
            if r is None:continue
            labels,sens=r;key=context.key(row);c=chrom_name(row['chrom']);pos=int(row['pos'])
            if state['multi_ALT'][0]:
                for sensitivity in sens:
                    other=total[('multi_ALT',sensitivity,key)]-per_chrom[(c,'multi_ALT',sensitivity,key)]
                    rate=(other[1]+.5)/(other[0]+1.) if other[0]>=a.min_background_sites else None
                    local_regions.record(c,pos,state['multi_ALT'][1],rate,sensitivity,prior_only=bool(other[1]==0),count_mode=row['count_mode'])
            for outcome in selected:
                eligible,case=state[outcome]
                if not eligible:continue
                for s in sens:
                    other=total[(outcome,s,key)]-per_chrom[(c,outcome,s,key)]
                    for label in labels:
                        k=(outcome,s,label)
                        if other[0]<a.min_background_sites:unsupported[k]+=1;continue
                        rate=(other[1]+.5)/(other[0]+1.)
                        denom[k]+=1;observed[k]+=int(case);expected[k]+=rate
                        block[k][(c,(pos-1)//a.block_bp)]+=(int(case),rate)
            if landscape and a.profile_sites>0 and state['multi_ALT'][0] and ('exclude_mask' in sens or not a.exclude_bed):
                rank=stable_hash(a.seed,c,pos,row['ref']);record=dict(chrom=c,pos=pos,ref=row['ref'],key=key)
                if state['multi_ALT'][1]:cases.add('all',rank,record)
                else:controls.add((c,key),rank,record)
    if workers==1:context.close()
    rows=[]
    for k in sorted(set(observed)|set(unsupported)):
        o,e=observed[k],expected[k];lo,hi,p,nblocks=bootstrap_oe(block[k],a.bootstrap,a.seed)
        rows.append(dict(outcome=k[0],sensitivity=k[1],region=k[2],observed=o,expected=e,oe=o/e if e else None,
            ci_low=lo,ci_high=hi,p=p,n_eligible_supported=denom[k],n_unsupported=unsupported[k],n_blocks=nblocks,
            inference='conditional_fixed_background_block_bootstrap' if p is not None else 'insufficient_information'))
    for outcome in selected:
        for s in sensitivities:
            indexes=[i for i,r in enumerate(rows) if r['outcome']==outcome and r['sensitivity']==s]
            for i,v in zip(indexes,bh([rows[i]['p'] for i in indexes])):rows[i]['fdr']=v
    fields=['outcome','sensitivity','region','observed','expected','oe','ci_low','ci_high','p','fdr','n_eligible_supported','n_unsupported','n_blocks','inference']
    write_tsv(out/'results/enrichment.tsv',fields,rows)
    write_tsv(out/'results/status_counts.tsv',['measure','n'],[dict(measure=k,n=v) for k,v in sorted(status.items())])
    region_manifest=local_regions.finalize(rows)
    profiles=matched_profiles(a,landscape,cases,controls,out) if landscape and a.profile_sites>0 else []
    enrichment_plot(out,rows,profiles)
    manifest=dict(schema='allele3-enrichment-1',status='complete',genome_build='GRCh38',ref=a.ref,target=a.target,
        source_scan=str(a.scan_dir),population=a.population,reference=ref_meta,examined_rows=examined,
        parameters={k:getattr(a,k) for k in ('min_ac','min_an','min_af','threshold_kind','window_bp','block_bp','bootstrap','seed','outcomes')},
        background='Leave-one-chromosome-out empirical rate with Jeffreys smoothing in trinucleotide x maximum-observed-ALT-AF x corresponding-AN strata; minimum other-chromosome support enforced',
        context='unadjusted; no FASTA requested' if a.no_context else str(a.fasta),
        target_frame='Eligible observed variant positions. No invariant/all-genome callable denominator; AF-only files cannot identify singleton or joint REF occupancy.',
        uncertainty='Genomic block bootstrap conditional on fitted LOCO background. P uses centered log-ratio bootstrap; BH within outcome x sensitivity x this reference/target run.',
        profile_inference='Matched AVI profiles use whole-chromosome bootstrap to retain within-chromosome dependence from reused controls and overlapping neighborhoods. Tile-weighted approximations exclude the focal tile; one chromosome has no CI/P/FDR.',
        avi_circularity='AVI contains conservation and gnomAD-frequency-trained features. AVI-versus-gnomAD analyses are descriptive, not independent evidence of selection.',
        masks=[identity(p) for p in a.exclude_bed],outputs=['results/enrichment.tsv','results/status_counts.tsv','reports/enrich_test.html'])
    manifest['region_highlights']=dict(threshold=a.highlight_fold,**region_manifest['counts'])
    manifest['outputs']+=['qc/region_highlights.json',*region_manifest['outputs']]
    atomic_json(out/'qc/manifest.json',manifest)
    html_report(out/'reports/enrich_test.html',f'{a.target} enrichment in {a.ref} — GRCh38',[
        ('Estimand','<p>Conditional enrichment among eligible observed variant positions. Multi-ALT candidates and proven joint nucleotide states are separate. AF-only data cannot identify singleton counts.</p>'),
        ('Region estimates','<img src="../figures/enrichment.png">'+html_table(rows,fields,150)),
        ('All high-enrichment regions',f'<p><b>{region_manifest["counts"]["highlight_regions"]:,}</b> coordinate-resolved candidate rows have O/E ≥ {a.highlight_fold:g}. Every candidate is exported, including low-count candidates. <a href="highlights.html">Open the complete searchable catalogue and genome-wide plots</a>.</p>'),
        ('Near versus far AVI',('<img src="../figures/avi_distance.png">' if profiles else '')+'<p>Case/control profiles use fixed hash sampling and matching on chromosome, trinucleotide and primary-ALT AF/AN. The focal AVI tile is excluded. Profiles are coverage-weighted tile approximations; local SD is the average of tile SDs, not a ring-wide SD. Whole-chromosome bootstrap retains dependence from reused controls and overlapping neighborhoods; one chromosome has no CI/P/FDR. Per-metric effective coverage and actual control positions/weights are exported.</p>'),
        ('Interpretation','<p>AVI-versus-gnomAD is descriptive because AVI uses gnomAD-frequency proxies. HERV/archaic overlaps alone do not establish mechanism or balancing selection. Archaic intervals do not identify the ancestry of the target allele.</p><p><a href="../results/enrichment.tsv">Complete estimates</a> · <a href="../results/status_counts.tsv">State availability</a></p>')])
    print(json.dumps(dict(status='complete',out=str(out),estimates=len(rows),source_rows=examined)));return manifest


def matched_profiles(a,landscape,cases,controls,out):
    chosen=[r[2] for r in cases.data.get('all',[])];chosen.sort(key=lambda r:(chrom_key(r['chrom']),r['pos']))
    perchrom=defaultdict(lambda:np.zeros(2));audit={};radii=[int(x) for x in a.radii.split(',')]
    if radii!=sorted(set(radii)) or not radii or radii[0]<1:raise ValueError('Positive increasing --radii required')
    fields=['chrom','pos','metric','inner_bp','outer_bp','case_value','control_mean','difference','n_controls','effective_case_bp','tile_bp',
            'control_chrom','control_positions','control_weights','control_effective_bps','n_selected_controls','effective_control_mean_bp','effective_control_min_bp']
    pair_path=Path(out)/'results/avi_distance_pairs.tsv.gz';pair_path.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(pair_path,'wt',compresslevel=3) as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for case in chosen:
            pool=controls.data.get((case['chrom'],case['key']),[])
            if not pool:continue
            pool=sorted(pool,key=lambda x:stable_hash(a.seed,case['chrom'],case['pos'],x[2]['pos']))[:a.controls_per_case]
            if any(chrom_name(x[2]['chrom'])!=chrom_name(case['chrom']) for x in pool):raise ValueError('Profile controls must be matched on chromosome')
            pcase=landscape.profile(case['chrom'],case['pos'],radii)
            pcs=[landscape.profile(x[2]['chrom'],x[2]['pos'],radii) for x in pool]
            for j,r in enumerate(pcase):
                valid=[(pool[i][2],profile[j]) for i,profile in enumerate(pcs)
                       if j<len(profile) and profile[j]['value'] is not None and profile[j]['effective_bp']>0]
                if r['value'] is None or r['effective_bp']<=0 or not valid:continue
                ctrl=[v['value'] for c,v in valid];bp=[float(v['effective_bp']) for c,v in valid]
                ctrl_mean=float(np.mean(ctrl));diff=r['value']-ctrl_mean;key=(r['inner_bp'],r['outer_bp'],r['metric'])
                # Whole chromosomes preserve all dependence caused by shared
                # controls and overlapping neighborhoods within a chromosome.
                # The effect estimate still weights each matched case equally.
                perchrom[(key,case['chrom'])]+=(diff,1)
                current=audit.setdefault(key,dict(controls=Counter(),case_bp=0.,control_bp=0.,min_control_bp=math.inf,n=0))
                for c,v in valid:current['controls'][(chrom_name(c['chrom']),int(c['pos']))]+=1
                current['case_bp']+=r['effective_bp'];current['control_bp']+=float(np.mean(bp));current['min_control_bp']=min(current['min_control_bp'],min(bp));current['n']+=1
                record=dict(chrom=case['chrom'],pos=case['pos'],metric=r['metric'],inner_bp=r['inner_bp'],outer_bp=r['outer_bp'],case_value=r['value'],
                    control_mean=ctrl_mean,difference=diff,n_controls=len(valid),effective_case_bp=r['effective_bp'],tile_bp=r['tile_bp'],
                    control_chrom=case['chrom'],control_positions=','.join(str(c['pos']) for c,v in valid),
                    control_weights=','.join(format(1./len(valid),'.12g') for _ in valid),control_effective_bps=','.join(format(v,'.12g') for v in bp),
                    n_selected_controls=len(pool),effective_control_mean_bp=float(np.mean(bp)),effective_control_min_bp=min(bp))
                writer.writerow({k:fmt(record.get(k)) for k in fields})
    summary=[]
    for key in sorted(audit):
        blocks=[v for (k,c),v in perchrom.items() if k==key];estimate,lo,hi,p=bootstrap_mean(blocks,a.bootstrap,a.seed);current=audit[key]
        summary.append(dict(inner_bp=key[0],outer_bp=key[1],metric=key[2],case_minus_control=estimate,ci_low=lo,ci_high=hi,p=p,
            n_pairs=current['n'],n_blocks=len(blocks),n_chromosomes=len(blocks),n_unique_controls=len(current['controls']),max_control_reuse=max(current['controls'].values()),
            mean_effective_case_bp=current['case_bp']/current['n'],mean_effective_control_bp=current['control_bp']/current['n'],min_effective_control_bp=current['min_control_bp'],
            inference='chromosome_bootstrap' if len(blocks)>=2 and p is not None else 'insufficient_chromosomes'))
    for row,q in zip(summary,bh([r['p'] for r in summary])):row['fdr']=q
    write_tsv(out/'results/avi_distance_summary.tsv',['inner_bp','outer_bp','metric','case_minus_control','ci_low','ci_high','p','fdr','n_pairs','n_blocks',
        'n_chromosomes','n_unique_controls','max_control_reuse','mean_effective_case_bp','mean_effective_control_bp','min_effective_control_bp','inference'],summary)
    return summary


def bootstrap_mean(blocks,replicates,seed):
    v=np.asarray(blocks,float)
    if not len(v) or v[:,1].sum()<=0:return None,None,None,None
    point=v[:,0].sum()/v[:,1].sum()
    if len(v)<2:return float(point),None,None,None
    rng=np.random.default_rng(seed);boot=[]
    for _ in range(replicates):
        s=v[rng.integers(0,len(v),len(v))].sum(axis=0)
        if s[1]>0:boot.append(s[0]/s[1])
    if len(boot)<20:return float(point),None,None,None
    lo,hi=np.quantile(boot,[.025,.975]);p=(1+np.count_nonzero(np.abs(np.asarray(boot)-point)>=abs(point)))/(len(boot)+1)
    return float(point),float(lo),float(hi),float(p)


def enrichment_plot(out,rows,profiles):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    (out/'figures').mkdir(parents=True,exist_ok=True)
    chosen=[r for r in rows if r['outcome']=='multi_ALT' and r['oe'] is not None]
    fig,ax=plt.subplots(figsize=(11,max(4,min(16,len(chosen)*.4))))
    for i,r in enumerate(chosen):
        ax.plot(r['oe'],i,'o',color='#246f91')
        if r['ci_low'] is not None:ax.plot([r['ci_low'],r['ci_high']],[i,i],color='#246f91')
    ax.set_yticks(range(len(chosen)),[r['region']+' / '+r['sensitivity'] for r in chosen]);ax.axvline(1,color='#899ba3',ls='--');ax.set_xlabel('Observed / LOCO expected multi-ALT sites');fig.tight_layout();fig.savefig(out/'figures/enrichment.png',dpi=170);plt.close(fig)
    if profiles:
        fig,axes=plt.subplots(2,3,figsize=(14,8))
        for ax,m in zip(axes.flat,METRICS):
            rs=[r for r in profiles if r['metric']==m];xs=np.arange(len(rs));ys=[r['case_minus_control'] for r in rs]
            ax.plot(xs,ys,'o-',color='#278a7a')
            for x,r in zip(xs,rs):
                if r['ci_low'] is not None:ax.vlines(x,r['ci_low'],r['ci_high'],color='#278a7a')
            ax.axhline(0,color='gray',ls='--');ax.set_xticks(xs,[f"{r['inner_bp']/1000:g}–{r['outer_bp']/1000:g}" for r in rs],rotation=30);ax.set_title(m);ax.set_xlabel('Distance band (kb)');ax.set_ylabel('Case − matched control')
        fig.tight_layout();fig.savefig(out/'figures/avi_distance.png',dpi=170);plt.close(fig)


def site_contrasts(row,predicted,a):
    """Reference-relative impact similarity, not molecular direction or fitness."""
    state,indices=observed_outcomes(row,a)
    if not state['multi_ALT'][0]:return []
    available=[i for i in indices if row['alt_afs'][i] is not None and row['alt_afs'][i]>0]
    ranked=sorted(available,key=lambda i:(-row['alt_afs'][i],row['alts'][i]))
    result=[]
    for scale in ('RAW','PHRED'):
        subset=[i for i in ranked if predicted.get(row['alts'][i],{}).get(scale) is not None]
        if len(subset)<2:continue
        x=np.log10([row['alt_afs'][i] for i in subset]);y=np.asarray([predicted[row['alts'][i]][scale] for i in subset])
        dx=x-x.mean();dy=y-y.mean();sxx=float(dx@dx);sxy=float(dx@dy)
        result.append(dict(metric='within_site_log10_AF_slope',scale=scale,numerator=sxy,denominator=sxx))
        if len(ranked)==3 and len(subset)==3 and len(set(row['alt_afs'][i] for i in ranked))==3:
            v1,v2,v3=[predicted[row['alts'][i]][scale] for i in ranked]
            result.append(dict(metric='least_ALT_distance_to_middle_minus_top',scale=scale,numerator=abs(v3-v2)-abs(v3-v1),denominator=1.))
        # Two observed positive ALTs vs the third theoretically possible ALT.
        # The latter is "not reported positive", never declared extinct.
        observed_values=row['alt_afs'] if row['count_mode']=='AF_only' else row['alt_acs']
        reported_positive={b for b,value in zip(row['alts'],observed_values) if value is not None and value>0}
        remaining=[b for b in predicted if b not in reported_positive and predicted[b].get(scale) is not None]
        if len(reported_positive)==2 and len(ranked)==2 and len(subset)==2 and len(remaining)==1 and row['alt_afs'][ranked[0]]!=row['alt_afs'][ranked[1]]:
            top,extra=[predicted[row['alts'][i]][scale] for i in ranked];unreported=predicted[remaining[0]][scale]
            result.append(dict(metric='extra_ALT_distance_to_top_minus_unreported',scale=scale,numerator=abs(extra-top)-abs(unreported-top),denominator=1.))
    return result


def allele_test(a):
    out=Path(a.out)
    for sub in ('tmp/positions','tmp/avi','results','qc','figures','reports'): (out/sub).mkdir(parents=True,exist_ok=True)
    paths={};handles={};n_positions=Counter();seen={};status=Counter()
    try:
        for row in source_rows(a):
            state,idx=observed_outcomes(row,a)
            if not state['multi_ALT'][1]:continue
            c=chrom_name(row['chrom']);pos=int(row['pos'])
            if c not in handles:
                paths[c]=out/'tmp/positions'/f'{c}.txt';handles[c]=paths[c].open('w')
            if seen.get(c)!=pos:handles[c].write(str(pos)+'\n');n_positions[c]+=1;seen[c]=pos
    finally:
        for f in handles.values():f.close()
    scan_manifest=Path(a.scan_dir)/'qc/manifest.json'
    request=dict(scan=identity(scan_manifest),population=a.population,min_ac=a.min_ac,min_af=a.min_af,min_an=a.min_an,threshold_kind=a.threshold_kind)
    for c in sorted(paths,key=chrom_key):
        vcf=resolve_vcf(a.atlas_dir,c);dest=out/'tmp/avi'/f'{c}.tsv';receipt=out/'qc'/f'avi_{c}.json'
        desired=dict(request=request,atlas=identity(vcf),n_targets=n_positions[c],worker=identity(Path(__file__).with_name('avi_stream.cpp')),
            native=identity(a.native_binary) if a.native_binary else None)
        old=json.loads(receipt.read_text()) if receipt.is_file() else {}
        if a.replace=='FALSE' and dest.is_file() and old.get('status')=='complete' and old.get('output_identity')==identity(dest) and all(old.get(k)==v for k,v in desired.items()):
            print('[allele AVI reuse] '+c,flush=True)
        else:
            query_file(vcf,c,paths[c],dest,out/'tmp/bin',native=a.native_binary)
            atomic_json(receipt,dict(desired,status='complete',output_identity=identity(dest)))
    blocks=defaultdict(lambda:defaultdict(lambda:np.zeros(2)));site_n=Counter();plots=Reservoir(15000)
    fields=['chrom','pos','ref','population','alt','alt_AF_rank','rank_tie','af','ac','an','avi_raw','avi_phred','count_mode','ref_evidence']
    metric_fields=['chrom','pos','scale','metric','sensitivity','value','numerator','denominator']
    current=None;score_handle=None;score_iter=None;score=None;geometry_status=Counter()
    with gzip.open(out/'results/allele_avi.tsv.gz','wt',compresslevel=3) as afile, \
            gzip.open(out/'results/site_contrasts.tsv.gz','wt',compresslevel=3) as mfile, \
            gzip.open(out/'results/allele_geometry.tsv.gz','wt',compresslevel=3) as gfile:
        writer=csv.DictWriter(afile,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        mw=csv.DictWriter(mfile,fieldnames=metric_fields,delimiter='\t',lineterminator='\n');mw.writeheader()
        gw=csv.DictWriter(gfile,fieldnames=GEOMETRY_FIELDS,delimiter='\t',lineterminator='\n');gw.writeheader()
        try:
            for row in source_rows(a):
                state,idx=observed_outcomes(row,a)
                if not state['multi_ALT'][1]:continue
                c=chrom_name(row['chrom']);pos=int(row['pos']);status['eligible_sites']+=1
                if c!=current:
                    if score_handle:score_handle.close()
                    score_handle=(out/'tmp/avi'/f'{c}.tsv').open();score_iter=iter(csv.DictReader(score_handle,delimiter='\t'));score=next(score_iter,None);current=c
                while score is not None and int(score['pos'])<pos:score=next(score_iter,None)
                if score is None or int(score['pos'])!=pos:status['missing_atlas_position']+=1;continue
                if score['ref']!=row['ref']:status['reference_mismatch']+=1;continue
                raw=[number(x) for x in score['avi_raw'].split(',')];phred=[number(x) for x in score['avi_phred'].split(',')];alts=score['alts'].split(',')
                if len(raw)==1 and raw[0] is None:raw=[None]*len(alts)
                if len(phred)==1 and phred[0] is None:phred=[None]*len(alts)
                if len(raw)!=len(alts) or len(phred)!=len(alts):raise ValueError('Malformed Number=A AVI extraction')
                predicted={b:{'RAW':r,'PHRED':p} for b,r,p in zip(alts,raw,phred)}
                ranked=sorted([i for i in idx if row['alt_afs'][i] is not None and row['alt_afs'][i]>0],key=lambda i:(-row['alt_afs'][i],row['alts'][i]))
                frequencies=Counter(row['alt_afs'][i] for i in ranked)
                for rank,i in enumerate(ranked,1):
                    b=row['alts'][i];v=predicted.get(b,{})
                    r=dict(chrom=c,pos=pos,ref=row['ref'],population=a.population,alt=b,alt_AF_rank=rank,rank_tie=frequencies[row['alt_afs'][i]]>1,
                        af=row['alt_afs'][i],ac=row['alt_acs'][i],an=row['alt_ans'][i],avi_raw=v.get('RAW'),avi_phred=v.get('PHRED'),count_mode=row['count_mode'],ref_evidence=row.get('ref_evidence'))
                    writer.writerow({k:fmt(r.get(k)) for k in fields})
                    if v.get('PHRED') is not None:plots.add('all',stable_hash(a.seed,c,pos,b),r)
                contrasts=site_contrasts(row,predicted,a);status['matched_sites']+=1
                geometry,normalized=analyze_site(row,predicted,idx)
                for g in geometry:
                    gw.writerow({k:fmt(g.get(k)) for k in GEOMETRY_FIELDS})
                    geometry_status[(g['scale'],g['geometry_status'])]+=1
                contrasts.extend(normalized)
                reported_ans=[row['alt_ans'][i] for i in ranked]
                sensitivities=['all_available']
                if reported_ans and all(x is not None for x in reported_ans) and len(set(reported_ans))==1:sensitivities.append('equal_reported_AN')
                for v in contrasts:
                    if v['denominator']<=0:continue
                    for s in sensitivities:
                        key=(v['scale'],v['metric'],s);block=(c,(pos-1)//a.block_bp)
                        blocks[key][block]+=(v['numerator'],v['denominator']);site_n[key]+=1
                        mr=dict(chrom=c,pos=pos,scale=v['scale'],metric=v['metric'],sensitivity=s,value=v['numerator']/v['denominator'],numerator=v['numerator'],denominator=v['denominator'])
                        mw.writerow({k:fmt(mr.get(k)) for k in metric_fields})
        finally:
            if score_handle:score_handle.close()
    results=[]
    for key,b in sorted(blocks.items()):
        estimate,lo,hi,p=bootstrap_mean(list(b.values()),a.bootstrap,a.seed)
        results.append(dict(scale=key[0],metric=key[1],sensitivity=key[2],estimate=estimate,ci_low=lo,ci_high=hi,p=p,n_sites=site_n[key],n_blocks=len(b)))
    for r,q in zip(results,bh([r['p'] for r in results])):r['fdr']=q
    write_tsv(out/'results/within_site_tests.tsv',['scale','metric','sensitivity','estimate','ci_low','ci_high','p','fdr','n_sites','n_blocks'],results)
    write_tsv(out/'results/matching_status.tsv',['status','n'],[dict(status=k,n=v) for k,v in sorted(status.items())])
    geometry_counts=[dict(scale=k[0],status=k[1],n=v) for k,v in sorted(geometry_status.items())]
    write_tsv(out/'results/geometry_status.tsv',['scale','status','n'],geometry_counts)
    within_plot(out,[x[2] for x in plots.data.get('all',[])],results)
    geometry_plot(out,results)
    manifest=dict(schema='allele3-within-site-2',status='complete',genome_build='GRCh38',target=a.target,population=a.population,
        request=request,counts=dict(status),inference='Site fixed-effect slope and paired scalar impact-distance contrasts; genomic block bootstrap; BH across the reported tests.',
        rank='Positive observed ALT alleles sorted by AF; ties flagged and excluded from rank-dependent contrasts. ALT rank is neither global base rank nor birth order.',
        reference='AVI(REF) is unavailable, not filled with zero. No REF-versus-ALT impact-distance test is fabricated.',
        equal_AN='Equal reported AN is a sensitivity label, not proof of the same called chromosome set.',
        meaning='Similarity of reference-relative AVI magnitudes, not directional agreement of molecular functions or fitness. Unreported positive ALT is not proven absent or extinct.',
        geometry='Primary normalized weak-pair contrast: [d23-(d12+d13)/2]/range, where ranks refer to three observed positive, threshold-qualified ALTs; genomic-block inference. RAW and PHRED remain separate.',
        geometry_status=geometry_counts,
        permutation='Six score-to-fixed-AF-rank assignments are descriptive exchangeability diagnostics, not per-site discovery P-values or a mutation/selection null.',
        avi_circularity='gnomAD frequency contributes to AVI training; use All of Us for cohort sensitivity, while acknowledging shared population biology and possible cohort overlap.',
        outputs=['results/within_site_tests.tsv','results/allele_avi.tsv.gz','results/site_contrasts.tsv.gz',
            'results/allele_geometry.tsv.gz','results/geometry_status.tsv','figures/allele_geometry.png','figures/allele_geometry.pdf','reports/allele_test.html'])
    atomic_json(out/'qc/manifest.json',manifest)
    html_report(out/'reports/allele_test.html',f'Within-site AF and AVI — {a.target}',[
      ('Tests','<p>Positive observed ALT alleles are ranked by AF. REF is not assigned an AVI score. RAW and PHRED scales are both reported.</p>'+html_table(results,['scale','metric','sensitivity','estimate','ci_low','ci_high','p','fdr','n_sites'])),
      ('Frequency and impact','<img src="../figures/af_avi.png">'),
      ('Do the two lower-frequency ALTs have closer AVI magnitudes?',
       '<p>With three observed ALTs ranked by AF, the primary statistic is [d23 − (d12+d13)/2] / score range. Negative values mean the two lower-AF ALTs are closer to one another than the average of their distances to the highest-AF ALT. Each eligible site contributes equally to this normalized contrast; the aggregate uses genomic-block bootstrap.</p>'
       '<img src="../figures/allele_geometry.png"><p>RAW and PHRED are different scales. These distances measure similarity of reference-relative AVI magnitudes; they are not direct ALT-to-ALT molecular effects or pleiotropy. AF ties, missing predictions and zero score range have explicit status labels.</p>'
       '<p><a href="../results/allele_geometry.tsv.gz">Every matched site: identities, ranks, distances and six-assignment diagnostics</a> · <a href="../results/geometry_status.tsv">Availability and exclusion reasons</a></p>'
       '<p>The six-assignment tail fractions are descriptive rank-exchangeability diagnostics. They do not test mutation or selection and are not used to declare individual sites significant.</p>'),
      ('Interpretation','<p>A negative least_ALT_distance_to_middle_minus_top means the least frequent ALT has a more similar AVI magnitude to the middle-frequency ALT than to the highest-frequency ALT. It does not establish matching molecular directions, cooperation or balancing selection.</p><p>The extra_ALT_distance_to_top_minus_unreported contrast compares the second observed ALT to the theoretically possible ALT not reported positive. Missing detectability prevents calling the latter extinct.</p><p>Equal reported AN does not establish a common called set. gnomAD-trained AVI introduces frequency circularity.</p><p><a href="../results/within_site_tests.tsv">Tests</a> · <a href="../results/allele_avi.tsv.gz">Exact allele scores</a> · <a href="../results/matching_status.tsv">Availability</a></p>')])
    print(json.dumps(dict(status='complete',out=str(out),counts=dict(status))));return manifest


def within_plot(out,rows,results):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    if rows:
        x=np.log10([r['af'] for r in rows]);y=[r['avi_phred'] for r in rows]
        axes[0].hexbin(x,y,gridsize=45,mincnt=1,cmap='viridis');axes[0].set(xlabel='log10 observed ALT AF',ylabel='AVI PHRED',title='Hash sample of matched observed ALTs')
    else:axes[0].text(.5,.5,'No matched AVI scores',ha='center',transform=axes[0].transAxes)
    rs=[r for r in results if r['scale']=='PHRED' and r['sensitivity']=='all_available' and r['metric']!='within_site_log10_AF_slope' and not r['metric'].startswith('normalized_')]
    for i,r in enumerate(rs):
        axes[1].plot(r['estimate'],i,'o',color='#278a7a')
        if r['ci_low'] is not None:axes[1].plot([r['ci_low'],r['ci_high']],[i,i],color='#278a7a')
    axes[1].axvline(0,color='gray',ls='--');axes[1].set_yticks(range(len(rs)),[r['metric'].replace('_',' ') for r in rs]);axes[1].set_xlabel('Paired impact-distance difference (PHRED)')
    fig.tight_layout();fig.savefig(out/'figures/af_avi.png',dpi=180);plt.close(fig)


def geometry_plot(out,results):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels={
        'normalized_weak_pair_distance_minus_other_pairs':'Two lower-AF ALTs vs pairs involving highest AF',
        'normalized_least_ALT_distance_to_middle_minus_top':'Lowest AF: distance to middle minus highest',
        'normalized_extra_ALT_distance_to_top_minus_unreported':'Second observed vs non-positive/unreported ALT'}
    fig,axes=plt.subplots(1,2,figsize=(14,4.8),layout='constrained',sharey=True)
    for ax,scale in zip(axes,['RAW','PHRED']):
        selected={r['metric']:r for r in results if r['scale']==scale and r['sensitivity']=='all_available'}
        for i,key in enumerate(labels):
            r=selected.get(key)
            if r is None:continue
            ax.scatter(r['estimate'],i,color='#007d79',s=34,zorder=3)
            if r['ci_low'] is not None:ax.hlines(i,r['ci_low'],r['ci_high'],color='#7d969d',lw=2)
        ax.axvline(0,color='#805d48',ls='--',lw=1);ax.set_xlim(-1.05,1.05);ax.set_ylim(2.55,-.55)
        ax.set_title(scale+' score geometry');ax.set_xlabel('Mean normalized distance contrast (95% block interval)');ax.grid(axis='x',alpha=.18)
        if not any(k in selected for k in labels):ax.text(.5,.5,'No eligible normalized contrast',ha='center',transform=ax.transAxes)
    axes[0].set_yticks(range(len(labels)),list(labels.values()),fontsize=9)
    fig.suptitle('Within-site AF ranks and AVI magnitude similarity\nNegative values favor the specified closer-pair hypothesis',fontsize=13)
    for ext in ('png','pdf'):fig.savefig(Path(out)/'figures'/('allele_geometry.'+ext),dpi=180)
    plt.close(fig)


def parser():
    root=os.environ.get('ALLELE_ANALYSIS_ROOT',configured_path('ANALYSIS_ROOT', 'allele3'))
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['enrich_test','allele_test'])
    p.add_argument('--jobs',type=int,default=0,help='Enrichment chromosome workers; 0 uses allocated cores')
    p.add_argument('--ref',choices=['alphagenome','herv','archaic'],default='alphagenome');p.add_argument('--target',choices=['gnomad','allofus','custom'],default='gnomad')
    p.add_argument('--scan-dir');p.add_argument('--out');p.add_argument('--chr',default='autosomes');p.add_argument('--population',default='ALL')
    p.add_argument('--avi-dir',default=os.environ.get('ALLELE_AVI_STORE',root+'/avi_scan'))
    p.add_argument('--atlas-dir',default=os.environ.get('ALLELE_ATLAS_DIR',configured_path('DIR_ANNOT', 'alphagenome/38')))
    p.add_argument('--herv-scan-dir',default=os.environ.get('ALLELE_HERV_SCAN_STORE',root+'/herv_scan'));p.add_argument('--herv-bed');p.add_argument('--stable-herv-bed')
    p.add_argument('--erv-joined',default=os.environ.get('ALLELE_ERV_JOINED'))
    p.add_argument('--archaic-dir',default=os.environ.get('ALLELE_IBDMIX_SEGMENTS',configured_path('ANALYSIS_ROOT', 'gu/ibdmix/1kg/38/segments')))
    p.add_argument('--archaic-population',default='ALL');p.add_argument('--archaic-reference',default='all')
    p.add_argument('--archaic-min-carriers',type=int,default=2);p.add_argument('--archaic-min-slod',type=float,default=0);p.add_argument('--archaic-min-bp',type=int,default=0)
    p.add_argument('--fasta',default=os.environ.get('ALLELE_FASTA',configured_path('DIR_GEN', '1hgp/38/GRCH38.fasta.gz')));p.add_argument('--no-context',action='store_true')
    p.add_argument('--exclude-bed',action='append',default=[]);p.add_argument('--outcomes',default=','.join(OUTCOMES))
    p.add_argument('--threshold-kind',choices=['auto','AC','AF'],default='auto');p.add_argument('--min-ac',type=int,default=5);p.add_argument('--min-an',type=int,default=100);p.add_argument('--min-af',type=float,default=.0001)
    p.add_argument('--window-bp',type=int,default=10000);p.add_argument('--block-bp',type=int,default=1000000);p.add_argument('--min-background-sites',type=int,default=100)
    p.add_argument('--region-windows',default='10000,100000,1000000',help='Coordinate-resolved multi-ALT enrichment windows, in bp')
    p.add_argument('--highlight-fold',type=float,default=3.,help='Export every finite region O/E >= this value, including low-count candidates')
    p.add_argument('--bootstrap',type=int,default=1000);p.add_argument('--seed',type=int,default=20261115)
    p.add_argument('--profile-sites',type=int,default=10000);p.add_argument('--controls-per-stratum',type=int,default=100);p.add_argument('--controls-per-case',type=int,default=3)
    p.add_argument('--radii',default='1000,10000,100000,1000000');p.add_argument('--native-binary');p.add_argument('--replace',choices=['TRUE','FALSE'],default='FALSE')
    return p


def main(argv=None):
    a=parser().parse_args(argv);root=os.environ.get('ALLELE_ANALYSIS_ROOT',configured_path('ANALYSIS_ROOT', 'allele3'))
    if a.scan_dir is None:a.scan_dir=str(Path(os.environ.get('ALLELE_SUMMARY_STORE',root+'/allele3_scan'))/a.target)
    if a.out is None:a.out=root+'/'+a.command+'/'+(a.ref+'_'+a.target if a.command=='enrich_test' else a.target)
    if a.block_bp<1 or a.bootstrap<20 or a.min_ac<1 or a.min_an<1 or not 0<=a.min_af<=1:raise ValueError('Invalid threshold/block/bootstrap parameter')
    if a.archaic_min_carriers<1 or a.profile_sites<0 or a.controls_per_case<1 or a.controls_per_stratum<1 or a.min_background_sites<1:raise ValueError('Invalid carrier/profile/background parameter')
    source=Path(a.scan_dir)/'qc/manifest.json'
    if not source.is_file():raise ValueError('Run 3.allele3_scan.sh for this target first: '+str(source))
    meta=json.loads(source.read_text())
    if meta.get('status')!='complete':raise ValueError('Target scan is not complete')
    if meta.get('build',meta.get('genome_build'))!='GRCh38':raise ValueError('Target scan must use GRCh38')
    if meta.get('target')!=a.target:raise ValueError('--target does not match the target scan manifest')
    config=meta.get('config',{})
    columns=(config.get('schema') or {}).get('populations',{}).get(a.population)
    if a.threshold_kind=='AC' and config.get('format')=='tsv' and columns is not None and not columns.get('ac'):
        raise ValueError('--threshold-kind AC requires reported AC counts; population '+a.population+' is AF-only. Use --threshold-kind auto or AF.')
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    request=dict(parameters={k:v for k,v in vars(a).items() if k not in ('replace','out')},scan=identity(source),program=identity(__file__),
        helpers=[identity(Path(__file__).with_name(name)) for name in ('research_common.py','summary_scan.py','avi_scan.py','avi_stream.cpp','enrichment_regions.py','enrichment_parallel.py','allele_geometry.py')])
    request_path=out/'qc/request.json';manifest=out/'qc/manifest.json'
    old=json.loads(request_path.read_text()) if request_path.is_file() else {}
    prior=json.loads(manifest.read_text()) if manifest.is_file() else {}
    # Input annotations are included without large-file hashing.
    refs=[]
    for path in [a.fasta if a.command=='enrich_test' and not a.no_context else None,a.herv_bed,a.stable_herv_bed,*a.exclude_bed]:
        if path:refs.append(identity(path))
    for path in [Path(a.avi_dir)/'qc/manifest.json' if a.ref=='alphagenome' and a.command=='enrich_test' else None]:
        if path:refs.append(identity(path))
    if a.command=='allele_test':refs.extend(identity(resolve_vcf(a.atlas_dir,c)) for c in chromosomes(a.chr))
    elif a.ref=='archaic':refs.extend(identity(p) for p in sorted(Path(a.archaic_dir).glob('*.segments.tsv.gz')))
    elif a.ref=='herv' and not a.herv_bed:
        available=next((p for p in [Path(a.herv_scan_dir)/'reference_herv.bed',Path(a.herv_scan_dir)/'stable_loci.bed',Path(a.erv_joined) if a.erv_joined else None] if p and p.is_file()),None)
        if available:refs.append(identity(available))
    if a.command=='enrich_test' and a.ref=='herv' and not a.stable_herv_bed:
        stable=Path(a.herv_scan_dir)/'stable_loci.bed'
        if stable.is_file():refs.append(identity(stable))
    request['reference_inputs']=refs
    if a.replace=='FALSE' and old==request and prior.get('status')=='complete' and all((out/p).is_file() for p in prior.get('outputs',[])):
        print(json.dumps(dict(status='reused',out=str(out))));return prior
    atomic_json(request_path,request);atomic_json(manifest,dict(status='running',command=a.command))
    try:return enrichment(a) if a.command=='enrich_test' else allele_test(a)
    except Exception as e:
        atomic_json(manifest,dict(status='failed',command=a.command,reason=str(e)));raise


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,RuntimeError) as e:raise SystemExit('allele3-research: '+str(e))
