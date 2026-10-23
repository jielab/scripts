"""Independent U0 full contrasts, U1 probability sample, U2 all strong, U3 budget."""
from __future__ import annotations
from collections import Counter, defaultdict
from pathlib import Path
import csv, gzip, heapq, json, math, random, shutil, sqlite3
from ancient_join import (table_rows,number,_match_stratum,_stage_identity,completed_stage,finish_stage)
from core import stable_hash,tsv_writer,dump_json,fingerprint,open_text
from common import ac_tier
from ldblocks import read_blocks,BlockIndex

EXTRA_FIELDS=['group','match_set','match_stratum','locus_block','selection_evidence','comparison_universe',
              'universe','sampling_stratum','sampling_probability','sampling_weight']

def eligible(row,min_extra_ac=5,exclude_mhc=True,allow_filter_warnings=False):
    return (number(row.get('extra_ac'),0)>=min_extra_ac and
        (row.get('selection_filter')=='PASS' or allow_filter_warnings) and
        row.get('has_non_snv')!='1' and row.get('region_flags','.') in ('.','','none') and
        not (exclude_mhc and row['chrom37']=='6' and 25000000<=int(row['pos37'])<=34000000) and
        number(row.get('s'),0)!=0)

def sampling_stratum(r):
    p=number(r['posterior'],0); s=abs(number(r['s'],0))
    pb='strong' if p>=.99 else 'moderate' if p>=.9 else 'intermediate' if p>.5 else 'weak'
    sb='zero' if not s else str(math.floor(math.log10(s)))
    return '|'.join([pb,sb,ac_tier(r['extra_ac']),r.get('canonical_context',r['context']),
                     r.get('population_class','NA'),r['n_modern_alt'],r['chrom']])

def sample_allocations(counts,budget,seed=20261022,with_expectations=False):
    """Fixed-size stratified sample with randomized fractional allocation.

    Expected stratum counts mu_h are proportional to sqrt(N_h) times declared
    oversampling factors, capped at N_h. Systematic randomized rounding gives
    E[K_h]=mu_h and sum K_h=budget, even when strata outnumber the budget.
    Independent SRS inside each stratum gives exact marginal pi=mu_h/N_h.
    """
    budget=min(budget,sum(counts.values()))
    if budget<0: raise ValueError('Sample budget must be nonnegative')
    weights={s:math.sqrt(n)*(4 if s.split('|')[0]=='strong' else 1)*(3 if s.split('|')[-2]=='3' else 1)*(2 if s.split('|')[2] in {'1','2-4'} else 1) for s,n in counts.items()}
    remaining=set(counts); expected={}; left=float(budget)
    while remaining:
        total=sum(weights[s] for s in sorted(remaining))
        full=[s for s in sorted(remaining) if left*weights[s]/total>=counts[s]]
        if not full:
            expected.update({s:left*weights[s]/total for s in remaining});break
        for st in full: expected[st]=float(counts[st]);left-=counts[st];remaining.remove(st)
    alloc={}; offset=random.Random('allocation:'+str(seed)).random(); cumulative=0.; previous=0
    ordered=sorted(counts)
    for i,st in enumerate(ordered):
        cumulative+=expected[st]
        rounded=budget if i==len(ordered)-1 else math.floor(cumulative+offset)
        alloc[st]=rounded-previous; previous=rounded
        if not 0<=alloc[st]<=counts[st]: raise RuntimeError('Invalid sampling allocation')
    return (alloc,expected) if with_expectations else alloc


def _decorate(r,universe,group,block,stratum='NA',probability=1,match='.'):
    return dict(r,group=group,match_set=match,match_stratum=_match_stratum(r),locus_block=block,
        selection_evidence='published_strong' if group=='strong' else 'weak_signal_not_neutral' if group=='weak_signal' else 'not_selected_by_significance',
        comparison_universe='modern_multi_alt_exact_intersection',universe=universe,
        sampling_stratum=stratum,sampling_probability=probability,sampling_weight=1/probability)

def clustered_bootstrap(rows,out,seed=20261022,reps=1000):
    """Descriptive all-strong means, resampling whole GRCh37 blocks with replacement."""
    groups=defaultdict(lambda:[0,0.,0.])
    for r in rows:
        z=groups[r['locus_block']]; z[0]+=1; z[1]+=float(r['extra_af']); z[2]+=float(r['extra_ac'])
    values=list(groups.values()); rng=random.Random(seed); result=[]
    for index,name in ((1,'extra_af'),(2,'extra_ac')):
        means=[]
        if len(values)>=2:
            for _ in range(reps):
                chosen=[values[rng.randrange(len(values))] for _ in values]
                means.append(sum(z[index] for z in chosen)/sum(z[0] for z in chosen))
            means.sort()
        result.append(dict(metric=name,n_rows=len(rows),n_clusters=len(values),
            mean=sum(z[index] for z in values)/len(rows) if rows else 'NA',
            ci_low=means[int(.025*len(means))] if means else 'NA',
            ci_high=means[min(len(means)-1,int(.975*len(means)))] if means else 'NA',
            method='whole_LD_block_bootstrap; outside_partition uses labelled 1Mb spatial clusters'))
    h,w=tsv_writer(out,list(result[0]))
    with h: w.writerows(result)

def build_universes(intersection,out_dir,*,ld_block_dir,u1_sites=20000,strong_posterior=.99,
                    weak_posterior=.5,min_extra_ac=5,controls_per_case=4,seed=20261022,
                    exclude_mhc=True,allow_filter_warnings=False):
    if u1_sites<1 or not 3<=controls_per_case<=5 or not 0<=weak_posterior<strong_posterior<=1:
        raise ValueError('Invalid universe budget, posterior thresholds or controls (3-5 required)')
    block_path=Path(ld_block_dir)/'EUR.37.bed'; blocks,_=read_blocks(block_path,'EUR','37'); index=BlockIndex(blocks)
    out=Path(out_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    params=dict(schema=2,u1_sites=u1_sites,strong_posterior=strong_posterior,weak_posterior=weak_posterior,
        min_extra_ac=min_extra_ac,controls_per_case=controls_per_case,seed=seed,
        exclude_mhc=exclude_mhc,allow_filter_warnings=allow_filter_warnings)
    identity=_stage_identity('independent_universes',dict(intersection=intersection,ld=block_path),params)
    if old:=completed_stage(out,identity): return old
    dbpath=out/'sampling.sqlite'
    if dbpath.exists(): dbpath.unlink()
    db=sqlite3.connect(dbpath); db.execute('PRAGMA cache_size=-48000')
    db.execute('CREATE TABLE sites (site TEXT PRIMARY KEY,stratum TEXT,rand TEXT,rank REAL,payload TEXT,block TEXT)')
    db.execute('CREATE TABLE weak (site TEXT,stratum TEXT,rand TEXT,payload TEXT,block TEXT)')
    fields=None; strong=[]; qc=Counter()
    def block(r):
        b=index.locate(r['chrom37'],r['pos37'])
        return b.block_id if b else f'outside_partition:{r["chrom37"]}:{(int(r["pos37"])-1)//1000000}'
    # U0 is byte-preserved if compressed; otherwise losslessly recompressed.
    if str(intersection).endswith('.gz'): shutil.copyfile(intersection,out/'U0.tsv.gz')
    else:
        with open(intersection,'rb') as src,gzip.open(out/'U0.tsv.gz','wb') as dst: shutil.copyfileobj(src,dst)
    for r in table_rows(intersection):
        if fields is None: fields=list(r)
        qc['U0_rows']+=1; st=sampling_stratum(r); b=block(r)
        rank=float(r['posterior']); payload=json.dumps(r); h=stable_hash([seed,r['site_id']])
        # One deterministic representative per unique site for probability sampling.
        db.execute('INSERT INTO sites VALUES (?,?,?,?,?,?) ON CONFLICT(site) DO UPDATE SET stratum=excluded.stratum,rank=excluded.rank,payload=excluded.payload,block=excluded.block WHERE excluded.rank>sites.rank',
            (r['site_id'],st,h,rank,payload,b))
        if eligible(r,min_extra_ac,exclude_mhc,allow_filter_warnings):
            p=float(r['posterior'])
            if p>=strong_posterior: strong.append(_decorate(r,'U2','strong',b)); qc['U2_all_strong_rows']+=1
            elif p<=weak_posterior:
                db.execute('INSERT INTO weak VALUES (?,?,?,?,?)',(r['site_id'],_match_stratum(r),h,payload,b))
        if qc['U0_rows']%50000==0: db.commit()
    db.commit(); db.execute('CREATE INDEX site_sample ON sites(stratum,rand)'); db.execute('CREATE INDEX weak_sample ON weak(stratum,rand)'); db.commit()
    fields=(fields or [])+EXTRA_FIELDS
    counts=dict(db.execute('SELECT stratum,COUNT(*) FROM sites GROUP BY stratum')); alloc,expected=sample_allocations(counts,u1_sites,seed,True)
    u1=[]
    for st,k in sorted(alloc.items()):
        for payload,b in db.execute('SELECT payload,block FROM sites WHERE stratum=? ORDER BY rand,site LIMIT ?',(st,k)):
            r=json.loads(payload); u1.append(_decorate(r,'U1','architecture_sample',b,st,expected[st]/counts[st]))
    qc['U0_unique_sites']=sum(counts.values()); qc['U1_unique_sites']=len(u1)
    strong.sort(key=lambda r:(-float(r['posterior']),-abs(float(r['x'])),r['site_id'],r['selection_id'],r['extra_alt']))
    strong_sites={r['site_id'] for r in strong}; used=set(strong_sites); controls=[]; match_report=[]
    for i,case in enumerate(strong,1):
        match=f'case{i:06d}'; case['match_set']=match; chosen=[]
        for site,_,_,payload,b in db.execute('SELECT * FROM weak WHERE stratum=? ORDER BY rand,site',(case['match_stratum'],)):
            if site in used or b==case['locus_block']: continue
            chosen.append(_decorate(json.loads(payload),'U2','weak_signal',b,match=match)); used.add(site)
            if len(chosen)==controls_per_case: break
        controls.extend(chosen)
        match_report.append(dict(match_set=match,site_id=case['site_id'],selection_id=case['selection_id'],extra_alt=case['extra_alt'],
            requested=controls_per_case,matched=len(chosen),status='complete' if len(chosen)==controls_per_case else 'insufficient_exact_matched_weak_controls'))
    qc['U2_controls']=len(controls); qc['U2_cases_under_3_controls']=sum(r['matched']<3 for r in match_report)
    reps={}
    for r in strong:
        if not r['locus_block'].startswith('outside_partition'): reps.setdefault(r['locus_block'],r)
    qc['U2_representatives']=len(reps); qc['U2_outside_partition_rows']=sum(r['locus_block'].startswith('outside_partition') for r in strong)
    outputs=[]
    def write(name,rows,cols=fields):
        h,w=tsv_writer(out/name,cols)
        with h: w.writerows(rows)
        outputs.append(name)
    write('U1.tsv',u1); write('U2_all_strong.tsv',strong); write('U2_controls.tsv',controls)
    write('U2_representatives.tsv',list(reps.values()))
    write('U2_matching.tsv',match_report,['match_set','site_id','selection_id','extra_alt','requested','matched','status'])
    write('sampling_design.tsv',[dict(stratum=s,population=counts[s],sample=alloc[s],expected_sample=expected[s],probability=expected[s]/counts[s],weight=counts[s]/expected[s]) for s in sorted(counts)],['stratum','population','sample','expected_sample','probability','weight'])
    combined={}
    for r in [*strong,*controls,*u1]: combined.setdefault((r['site_id'],r['selection_id'],r['extra_alt']),r)
    write('candidates.tsv',list(combined.values()))
    clustered_bootstrap(strong,out/'U2_clustered_bootstrap.tsv',seed)
    clustered_bootstrap(list(reps.values()),out/'U2_representative_bootstrap.tsv',seed)
    db.close()
    outputs+=['U0.tsv.gz','U2_clustered_bootstrap.tsv','U2_representative_bootstrap.tsv','sampling.sqlite']
    return finish_stage(out,identity,dict(qc=dict(qc),candidate_path=str(out/'candidates.tsv'),
        U0_path=str(out/'U0.tsv.gz'),U3_status='not_selected; run select-u3 after local LD/ERV/UKB review',
        caveats=['U0 is unsampled.','U1 marginal probability is mu_h/N_h; randomized systematic stratum allocation, independent within-stratum SRS; all units have positive probability.',
        'U2 includes every eligible strong contrast; no max-sites truncation.','Matching never relaxes scientific thresholds; deficits are explicit.',
        'Outside LD partitions cannot be labelled LD representatives.','Matching extra AF precludes paired AF-outcome tests.']),outputs)

def select_u3(universe_dir,out_path,sites=3000,seed=20261022):
    if sites<1: raise ValueError('U3 size must be positive')
    directory=Path(universe_dir);src=directory/'sampling.sqlite'
    # A separate functional budget: keep U2 cases and their controls first,
    # then fill from U0. U1 sampling is never consulted.
    selected={};matches=set()
    for r in table_rows(directory/'U2_all_strong.tsv'):
        if r['site_id'] in selected:continue
        if len(selected)>=sites:break
        selected[r['site_id']]=dict(r,universe='U3',sampling_stratum='U3_priority_U2')
        matches.add(r['match_set'])
    pools=defaultdict(list)
    for r in table_rows(directory/'U2_controls.tsv'):
        if r['match_set'] in matches:pools[r['match_set']].append(r)
    for round_ in range(5):
        for match in sorted(pools):
            if len(selected)>=sites:break
            if round_<len(pools[match]):
                r=pools[match][round_]
                selected.setdefault(r['site_id'],dict(r,universe='U3',sampling_stratum='U3_priority_U2'))
    priority=len(selected);budget=sites-priority
    db=sqlite3.connect(f'file:{src}?mode=ro',uri=True);heap=[];total=0;pool_n=0
    for site,payload,block in db.execute('SELECT site,payload,block FROM sites'):
        total+=1
        if site in selected:continue
        pool_n+=1
        if not budget:continue
        h=int(stable_hash(['U3',seed,site]),16);item=(-h,site,payload,block)
        if len(heap)<budget:heapq.heappush(heap,item)
        elif item[:2]>heap[0][:2]:heapq.heapreplace(heap,item)
    db.close();result_rows=list(selected.values())
    for _,_,payload,block in sorted(heap,reverse=True):
        r=json.loads(payload);result_rows.append(_decorate(r,'U3','functional_sample',block,'U0_independent_fill',min(budget,pool_n)/pool_n))
    if result_rows:
        h,w=tsv_writer(out_path,list(result_rows[0]))
        with h:w.writerows(result_rows)
    else:
        with open_text(directory/'U0.tsv.gz') as f:fields=next(csv.reader(f,delimiter='\t'))+EXTRA_FIELDS
        h,w=tsv_writer(out_path,fields);h.close()
    result=dict(n_sites=len(result_rows),U0_unique_sites=total,requested_sites=sites,U2_priority_sites=priority,
        independent_of_U1=True,source=fingerprint(src),design='One row per site; U2 strong then round-robin matched controls; independent U0 seeded fill; targeted functional panel')
    dump_json(str(out_path)+'.json',result);return result
