"""Competitive alignment validation of plant or human candidate fragments.

Outputs are compatibility evidence relative to an explicitly bounded reference
panel. Genome representatives are reused from the classification reference
collection: the alignment method is independent, but reference contamination
can be shared. Original assignments and split FASTQs are never modified.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import re
import subprocess

from result_output import managed_outputs, prepare as prepare_result, finish as finish_result, replace_bool
from common import atomic_json, code_identity, load_config, require_compute, sha256
from bin_layout import fastq_paths
from fastq_stream import fragments, read_id, new_digest, update_digest
from qc_paths import configured_split_root
from storage_layout import metadata_dir
from split_summary import load_sample, table
from species_refine import _load_representatives, _executable, _open_fasta
from species_split import _cigar_metrics
from assignment_io import open_assignment_text
from taxonomy import Taxonomy


def file_state(path):
    s=Path(path).stat()
    return dict(path=str(Path(path).resolve()),size=s.st_size,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)


def select_panel(claims, counts, reps, tax, extra, maximum, budget, background=()):
    wanted={t:'claimed_species' for t in claims}
    for tid in claims:
        lineage=list(tax.lineage(tid))
        # Closest available genus/family relatives, with explicit omissions.
        clades=[x for x in lineage if tax.rank.get(x) in ('genus','family')]
        close=set(clades[:2])
        for other in reps:
            if other!=tid and close.intersection(tax.lineage(other)):
                wanted.setdefault(other,'close_relative')
    for tid in extra: wanted.setdefault(tid,'explicit_competitor')
    wanted.setdefault(9606,'human_competitor')
    for tid in background: wanted.setdefault(tid,'other_claimed_species')
    groups=('bacteria','archaea','fungi','plant','virus','other_eukaryote','invertebrate','other_vertebrate','mammal')
    missing=[]
    for group in groups:
        ranked=sorted((t for t,r in reps.items() if r['group']==group),key=lambda t:(-counts[t],t))
        if ranked:
            for tid in ranked[:2]: wanted.setdefault(tid,group+'_competitor')
        else: missing.append(dict(species_taxid=0,reason='missing_group:'+group,bases=0))
    selected=[]; omitted=list(missing); total=0
    priority={'claimed_species':0,'human_competitor':1,'explicit_competitor':1,
              'bacteria_competitor':1,'archaea_competitor':1,'fungi_competitor':1,'plant_competitor':1,
              'close_relative':2,'other_claimed_species':3}
    order=sorted(wanted,key=lambda t:(priority.get(wanted[t],4),-counts[t],t))
    for tid in order:
        if tid not in reps:
            omitted.append(dict(species_taxid=tid,reason='representative_missing:'+wanted[tid],bases=0));continue
        row=reps[tid]; n=int(row['bases'])
        if not Path(row['resolved_path']).is_file():
            omitted.append(dict(species_taxid=tid,reason='reference_file_missing:'+wanted[tid],bases=n));continue
        reason='max_species' if len(selected)>=maximum else ('reference_budget' if total+n>budget else '')
        if reason:
            omitted.append(dict(species_taxid=tid,reason=reason+':'+wanted[tid],bases=n));continue
        selected.append({**row,'role':wanted[tid]});total+=n
    return selected,omitted,total


def plan_panels(claims, counts, reps, tax, extra, maximum, budget, group='plant'):
    """Assign every candidate to one declared, bounded competitive panel.

    The former all-species panel could not validate even 61 claimed species
    under its default max-species=60. Reserve half the nominal budget for
    competitors, serialize panels, and report rather than discard untestable
    claims. These are within-panel checks; their MAPQs are never compared.
    """
    target_taxids=sorted({x['species_taxid'] for x in claims if x['species_taxid']},
                        key=lambda tid:(-counts[tid],tid))
    panels=[]; untested={}; pending=[]; used=0
    target_budget=max(1,budget//2); target_count=max(1,maximum//2)

    def add(targets):
        selected,omitted,bases=select_panel(targets,counts,reps,tax,extra,maximum,budget,
                                           background=target_taxids)
        kept={r['species_taxid'] for r in selected}
        available=set(targets)&kept
        if available:
            panels.append(dict(panel_id=f'{len(panels)+1:03d}',target_taxids=sorted(available),
                               selected=selected,omitted=omitted,declared_bases=bases))
        for tid in set(targets)-kept:
            untested[tid]='target_reference_not_selected'

    for tid in target_taxids:
        row=reps.get(tid)
        reason=('representative_missing' if row is None else
                'reference_file_missing' if not Path(row['resolved_path']).is_file() else
                'single_reference_exceeds_panel_budget' if int(row['bases'])>budget else '')
        if reason:
            untested[tid]=reason;continue
        bases=int(row['bases'])
        if pending and (used+bases>target_budget or len(pending)>=target_count):
            add(pending);pending=[];used=0
        pending.append(tid);used+=bases
    if pending:add(pending)
    owners={tid:panel['panel_id'] for panel in panels for tid in panel['target_taxids']}

    # A higher-rank claim is tested only against a panel containing an actual
    # descendant. Its evidence remains at the original higher rank.
    broad=sorted({x['assigned_taxid'] for x in claims if not x['species_taxid']},
                 key=lambda tid:(-sum(x['assigned_taxid']==tid for x in claims),tid))
    broad_owners={}
    for tid in broad:
        candidates=[p for p in panels if any(tid in set(tax.lineage(r['species_taxid']))
                                             for r in p['selected'] if r['group']==group)]
        if not candidates:
            choices=sorted((t for t,r in reps.items() if r['group']==group and
                            int(r['bases'])<=budget and Path(r['resolved_path']).is_file() and
                            tid in set(tax.lineage(t))),key=lambda t:(-counts[t],int(reps[t]['bases']),t))
            if choices:
                before=len(panels);add([choices[0]])
                if len(panels)>before:candidates=[panels[-1]]
        if candidates:broad_owners[tid]=candidates[0]['panel_id']
        else:untested[tid]='no_available_descendant_reference_within_budget'
    assignment={x['ordinal']:(owners.get(x['species_taxid']) if x['species_taxid'] else
                             broad_owners.get(x['assigned_taxid'])) for x in claims}
    return panels,assignment,untested


def build_panel(selected, path, max_bases):
    contigs={};total=0
    with Path(path).open('wb') as out:
        for row in selected:
            tid=int(row['species_taxid']); idx=0; name=None
            with _open_fasta(row['resolved_path']) as stream:
                for raw in stream:
                    if raw.startswith(b'>'):
                        idx+=1; name=f'T{tid}|R{idx}'
                        desc=raw[1:].decode('utf-8','replace').strip()
                        lower=desc.lower()
                        cls=('chloroplast' if 'chloroplast' in lower or 'plastid' in lower else
                             'mitochondrial' if 'mitochond' in lower else
                             'nuclear' if re.search(r'\b(chromosome|nuclear)\b',lower) else 'unknown')
                        contigs[name]=dict(contig=name,species_taxid=tid,source_path=row['resolved_path'],
                                           original_header=desc,contig_class=cls,repeat_annotation='repeat' in lower,bases=0)
                        out.write(('>'+name+'\n').encode())
                    else:
                        seq=raw.strip()
                        if not seq:continue
                        if name is None or not re.fullmatch(b'[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+',seq):
                            raise RuntimeError('Malformed reference FASTA: '+row['resolved_path'])
                        total+=len(seq);contigs[name]['bases']+=len(seq)
                        if total>max_bases: raise RuntimeError('Actual FASTA exceeds declared panel budget')
                        out.write(seq+b'\n')
            if idx==0:raise RuntimeError('Empty reference '+row['resolved_path'])
    if not total: raise RuntimeError('Empty competitive reference')
    return contigs,total


def extract_candidates(folder, rows, assignment, out, taxids, broad, maximum, expected_total, group='plant'):
    chosen=[r for r in rows if r['group']==group and
            ((r['status']=='assigned_species' and (not taxids or int(r['species_taxid']) in taxids)) or
             (broad and r['status']=='higher_rank'))]
    bins={r['bin'] for r in chosen};expected=sum(int(r['fragments']) for r in chosen)
    if not expected:raise RuntimeError('No selected '+group+' candidates')
    if expected>maximum:raise RuntimeError('Candidate count exceeds --max-fragments; no reads silently dropped')
    receipt_path=metadata_dir(folder)/'complete.json'
    states=[file_state(receipt_path)]
    outputs=json.loads(receipt_path.read_text())['outputs']
    def verify(path):
        states.append(file_state(path))
        name=str(Path(path).relative_to(folder))
        saved=outputs.get(name,{}).get('state',{}).get('size')
        if saved is None or path.stat().st_size!=saved:
            raise RuntimeError('Committed artifact missing/size mismatch: '+name)
    verify(assignment)
    keys=('assigned_taxid','species_taxid','group','status','bin')
    expected_counts=Counter()
    for row in rows:expected_counts[tuple(str(row[k]) for k in keys)]+=int(row['fragments'])
    actual_counts=Counter()
    bybin=defaultdict(list);n=0
    with open_assignment_text(assignment) as h:
        for n,row in enumerate(csv.DictReader(h,delimiter='\t'),1):
            if int(row['ordinal'])!=n:raise RuntimeError('Noncontiguous saved assignment ordinals')
            actual_counts[tuple(row[k] for k in keys)]+=1
            if row['bin'] in bins:bybin[row['bin']].append(row)
    if n!=expected_total:raise RuntimeError('Assignments count differs from committed summary')
    if actual_counts!=expected_counts:raise RuntimeError('Assignment taxon/bin/status counts differ from counts.tsv')
    claims=[]; digest=new_digest()
    paths=[out/'candidates.R1.fastq.gz',out/'candidates.R2.fastq.gz']
    with gzip.open(paths[0],'wb') as aout,gzip.open(paths[1],'wb') as bout:
        for base in sorted(bins):
            if Path(base).is_absolute() or '..' in Path(base).parts:raise RuntimeError('Unsafe saved bin path')
            pair=fastq_paths(folder,base)
            for path in pair:verify(path)
            source=iter(bybin[base]);count=0
            for count,(a,b) in enumerate(fragments(*pair),1):
                row=next(source,None)
                if row is None or row['r1_read_id']!=read_id(a[0]).decode() or row['r2_read_id']!=read_id(b[0]).decode():
                    raise RuntimeError('Bin FASTQ order/ID differs from saved assignments')
                if not ((row['status']=='assigned_species' and (not taxids or int(row['species_taxid']) in taxids)) or (broad and row['status']=='higher_rank')):continue
                update_digest(digest,a,b);idx=len(claims)+1
                seqs=[x[1].strip().upper() for x in (a,b)]
                low=any(max(Counter(x).values())/len(x)>.8 or len({x[i:i+3] for i in range(max(0,len(x)-2))})<5 for x in seqs)
                claims.append(dict(ordinal=idx,original_ordinal=row['ordinal'],r1_read_id=row['r1_read_id'],r2_read_id=row['r2_read_id'],
                    assigned_taxid=int(row['assigned_taxid']),species_taxid=int(row['species_taxid']),status=row['status'],bin=base,low_complexity=low))
                for m,(record,handle) in enumerate(((a,aout),(b,bout)),1):
                    handle.write(f'@{idx}/{m}\n'.encode()+record[1].rstrip(b'\r\n')+b'\n+\n'+record[3].rstrip(b'\r\n')+b'\n')
            if next(source,None) is not None:raise RuntimeError('Bin FASTQ truncated relative to assignments')
    if len(claims)!=expected:raise RuntimeError('Extracted candidate count differs from counts.tsv')
    return claims,paths,states,digest.hexdigest()


def cigar_reference_span(cigar):
    """SAM reference span: M/D/N/=/X consume reference; I/S/H/P do not."""
    return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar) if op in 'MDN=X')


def parse_sam(path, expected, contigs, min_mapq, min_identity, min_fraction):
    """Return zero-based half-open reference intervals [position, end)."""
    align={}; alternatives=defaultdict(list)
    with Path(path).open() as h:
        for raw in h:
            if raw.startswith('@'):continue
            f=raw.rstrip('\r\n').split('\t')
            if len(f)<11:raise RuntimeError('Malformed/incomplete mapper SAM')
            flag=int(f[1]); q=f[0].split('/')[0]
            if flag&64 and flag&128:raise RuntimeError('SAM row has both read1 and read2 flags')
            mate=1 if flag&64 else 2 if flag&128 else int(f[0][-1]) if f[0].endswith(('/1','/2')) else 0
            if not q.isdigit() or not 1<=int(q)<=expected or not mate:raise RuntimeError('Unexpected mapper fragment/mate')
            key=(int(q),mate)
            if flag&2048:continue
            tags={x[:2]:x[5:] for x in f[11:] if len(x)>=5}
            mapped=not flag&4
            if mapped and f[2] not in contigs:raise RuntimeError('Unknown mapped contig')
            aligned,consumed=_cigar_metrics(f[5]) if mapped else (0,0)
            length=len(f[9]) if f[9]!='*' else consumed
            identity=max(0,1-int(tags.get('NM',aligned))/max(aligned,1)) if mapped else 0
            fraction=aligned/max(length,1);mapq=int(f[4]); score=int(tags.get('AS',0))
            row=dict(mapped=mapped,species_taxid=contigs[f[2]]['species_taxid'] if mapped else 0,
                     contig=f[2],position=max(0,int(f[3])-1),
                     end=max(0,int(f[3])-1)+cigar_reference_span(f[5]),
                     flag=flag,cigar=f[5],proper_pair=bool(flag&2),reverse=bool(flag&16),
                     mate_reverse=bool(flag&32),rnext=f[6],pnext=int(f[7]),tlen=int(f[8]),
                     mapq=mapq,identity=identity,aligned_fraction=fraction,
                     score=score,passed=mapped and min_mapq<=mapq<255 and identity>=min_identity and fraction>=min_fraction,
                     contig_class=contigs[f[2]]['contig_class'] if mapped else 'unknown',
                     repeat_annotation=contigs[f[2]]['repeat_annotation'] if mapped else False)
            if flag&256:alternatives[key].append(row);continue
            if key in align:raise RuntimeError('Duplicate primary SAM mate (possible multipart index)')
            align[key]=row
    if len(align)!=expected*2 or any((i,m) not in align for i in range(1,expected+1) for m in (1,2)):
        raise RuntimeError('Incomplete mapper output: expected exactly two primary records per fragment')
    for key,row in align.items():
        row['near_tied']=any(x['species_taxid']!=row['species_taxid'] and x['score']>=row['score']-5 for x in alternatives[key])
    return align


def category(claim, mates, tax):
    if not any(x['mapped'] for x in mates):return 'unmapped'
    if claim['low_complexity'] or any(x['near_tied'] for x in mates):return 'ambiguous'
    good=[x for x in mates if x['passed']]
    if len(good)!=2 or len({x['species_taxid'] for x in good})!=1:return 'ambiguous'
    target=good[0]['species_taxid']; original=tax.resolve(claim['assigned_taxid'])
    if claim['species_taxid'] and target==tax.resolve(claim['species_taxid']):
        # Organellar-only sequence is compatible with the clade, insufficient
        # by itself for a confident exact-species conclusion.
        if all(x['contig_class'] in ('chloroplast','mitochondrial') for x in good):return 'higher_rank_compatible'
        return 'supported_same_species'
    if original in set(tax.lineage(target)):return 'higher_rank_compatible'
    return 'conflicting_competitor'


def pair_evidence(mates):
    a,b=mates;both=all(x['mapped'] for x in mates)
    same=both and a['contig']==b['contig']
    proper=all(x['proper_pair'] for x in mates)
    reciprocal=same and all(x['rnext'] in ('=',y['contig']) and x['pnext']==y['position']+1
                            and x['mate_reverse']==y['reverse'] for x,y in ((a,b),(b,a)))
    concordant=bool(reciprocal and proper and a['tlen']==-b['tlen'] and a['tlen']!=0)
    return dict(cross_species_mates=both and a['species_taxid']!=b['species_taxid'],
                same_species_mates=both and a['species_taxid']==b['species_taxid'],
                same_contig=same,proper_pair=proper,concordant_pair=concordant,
                orientation=''.join('R' if x['reverse'] else 'F' for x in mates),
                insert_size=abs(a['tlen']) if same else 0,
                pairing_limitation='' if concordant else 'pair concordance not established; cross-scaffold hits may be compatible')


def support_class(x,m):
    cls=x[f'r{m}_contig_class']
    if cls in ('chloroplast','mitochondrial'):return 'organellar'
    if x[f'r{m}_repeat_annotation']:return 'repeat'
    return cls


def locus_stats(entries):
    regions=defaultdict(list);loci=Counter();classes=Counter()
    for x,m in entries:
        contig=x[f'r{m}_contig'];left=x[f'r{m}_position'];right=x[f'r{m}_end']
        regions[contig].append((left,right));loci[(contig,left//10000)]+=1
        classes[x[f'r{m}_contig_class']]+=1
    independent=0;covered=0
    for intervals in regions.values():
        end=-1
        for left,right in sorted(intervals):
            # Half-open adjacent intervals are separate nonoverlapping regions.
            if left>=end:independent+=1
            covered+=max(0,right-max(left,end))
            end=max(end,right)
    return dict(mates=len(entries),nonoverlapping_regions=independent,loci_10kb=len(loci),
                distinct_contigs=len(regions),covered_reference_bases=covered,
                maximum_locus_fraction=max(loci.values(),default=0)/max(1,len(entries)),
                organellar_only=bool(entries) and all(k in ('chloroplast','mitochondrial') for k in classes),
                repeat_annotated_mates=sum(x[f'r{m}_repeat_annotation'] for x,m in entries),
                unknown_class_mates=classes['unknown'])


def write_panel_reads(reads, claims, assigned, panel_id, out):
    """Keep each original pair together; numeric mapper IDs are panel-local."""
    selected=[]
    paths=[out/f'.panel.{panel_id}.R1.fastq.gz',out/f'.panel.{panel_id}.R2.fastq.gz']
    with gzip.open(paths[0],'wb') as left,gzip.open(paths[1],'wb') as right:
        seen=0
        for seen,(a,b) in enumerate(fragments(*reads),1):
            if assigned[seen]!=panel_id:continue
            selected.append(claims[seen-1]);ordinal=len(selected)
            for mate,(record,handle) in enumerate(((a,left),(b,right)),1):
                handle.write(f'@{ordinal}/{mate}\n'.encode()+record[1]+b'+\n'+record[3])
    if seen!=len(claims):raise RuntimeError('Extracted candidate FASTQs no longer match claims')
    return selected,paths


def evidence_row(claim, mates, tax, panel_id):
    result={**claim,'panel_id':panel_id,'category':category(claim,mates,tax),
            'tested':True,'untested_reason':'','mapq_scope':'within_panel',
            'both_mates_pass':all(x['passed'] for x in mates),
            'discordant_mates':all(x['mapped'] for x in mates) and mates[0]['species_taxid']!=mates[1]['species_taxid']}
    result.update(pair_evidence(mates))
    for m,row in enumerate(mates,1):result.update({f'r{m}_{k}':v for k,v in row.items()})
    return result


def aggregate_evidence(evidence, tax):
    groups=defaultdict(list)
    for row in evidence:groups[row['assigned_taxid']].append(row)
    aggregates=[]
    for tid,items in sorted(groups.items()):
        tested=[x for x in items if x['tested']]
        mapped=[(x,m) for x in tested for m in (1,2) if x[f'r{m}_mapped']]
        def matching(x,m):
            actual=x[f'r{m}_species_taxid']
            return (actual==tax.resolve(x['species_taxid']) if x['species_taxid'] else
                    tax.resolve(x['assigned_taxid']) in set(tax.lineage(actual)))
        support=[(x,m) for x,m in mapped if x['category'] in ('supported_same_species','higher_rank_compatible')
                 and x[f'r{m}_passed'] and matching(x,m)]
        result=dict(assigned_taxid=tid,name=tax.name.get(tax.resolve(tid),str(tid)),
                    fragments=len(items),mates=len(items)*2,tested_fragments=len(tested),
                    untested_fragments=len(items)-len(tested),
                    panel_ids=','.join(sorted({x['panel_id'] for x in items if x.get('panel_id')})),
                    untested_reasons=';'.join(sorted({x['untested_reason'] for x in items if x.get('untested_reason')})),
                    all_mapped_mates=len(mapped),all_mapped_fraction=len(mapped)/(2*len(items)),
                    claimed_support_mates=len(support),claimed_support_fraction=len(support)/(2*len(items)),
                    exact_species_support_mates=sum(bool(x['species_taxid']) for x,m in support),
                    higher_rank_compatible_mates=sum(not x['species_taxid'] for x,m in support),
                    both_mates_pass=sum(x['both_mates_pass'] for x in tested),
                    cross_species_mates=sum(x['cross_species_mates'] for x in tested),
                    concordant_support_fragments=sum(x['concordant_pair'] and x['category']=='supported_same_species' for x in tested),
                    same_species_only_fragments=sum(x['same_species_mates'] and not x['concordant_pair'] for x in tested),
                    one_mate_pass_fragments=sum(sum(x[f'r{m}_passed'] for m in (1,2))==1 for x in tested),
                    low_complexity_fragments=sum(x['low_complexity'] for x in items),
                    interpretation='bounded_panel_compatibility')
        strata={'all_mapped':mapped,'claimed_support':support}
        for cls in ('nuclear','organellar','repeat','unknown'):
            strata['claimed_'+cls]=[(x,m) for x,m in support if support_class(x,m)==cls]
        for cat in ('supported_same_species','higher_rank_compatible','ambiguous','conflicting_competitor','unmapped'):
            strata[cat]=[(x,m) for x,m in mapped if x['category']==cat]
        strata['one_mate_pass']=[(x,m) for x,m in mapped if sum(x[f'r{j}_passed'] for j in (1,2))==1 and x[f'r{m}_passed']]
        for name,entries in strata.items():
            result.update({name+'_'+key:value for key,value in locus_stats(entries).items()})
        result.update(Counter(x['category'] for x in items));aggregates.append(result)
    return aggregates


@managed_outputs
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--replace',type=replace_bool,default=True);p.add_argument('--config');p.add_argument('--root');p.add_argument('--sample',required=True);p.add_argument('--out')
    p.add_argument('--group',choices=('plant','human'),default='plant')
    p.add_argument('--representatives');p.add_argument('--taxid',type=int,action='append',default=[])
    p.add_argument('--competitor-taxid',type=int,action='append',default=[]);p.add_argument('--include-broad',action='store_true')
    p.add_argument('--max-species',type=int,default=60);p.add_argument('--max-reference-gib',type=float,default=24)
    p.add_argument('--max-fragments',type=int,default=1000000);p.add_argument('--threads',type=int,default=32)
    p.add_argument('--min-mapq',type=int,default=20);p.add_argument('--min-identity',type=float,default=.95)
    p.add_argument('--min-aligned-fraction',type=float,default=.8);p.add_argument('--minimap2');p.add_argument('--allow-local',action='store_true')
    p.add_argument('--keep-alignment-files',action='store_true',
                   help='Keep per-panel temporary FASTA, SAM and paired inputs (can be large)')
    a=p.parse_args(argv);require_compute(a.allow_local);cfg=load_config(a.config,a.root)
    initial_source=code_identity()
    taxonomy_states=[file_state(x) for x in sorted(Path(cfg['taxonomy_dir']).glob('*.dmp'))]
    if Path(a.sample).is_absolute() or '..' in Path(a.sample).parts:p.error('Unsafe sample ID')
    if min(a.max_species,a.max_fragments,a.threads)<1 or not math.isfinite(a.max_reference_gib) or a.max_reference_gib<=0:p.error('Limits must be positive')
    if not 0<=a.min_mapq<255 or not 0<=a.min_identity<=1 or not 0<=a.min_aligned_fraction<=1:p.error('Invalid alignment thresholds')
    output_group='plants' if a.group=='plant' else 'human-competitive'
    report_name='plant-validation.json' if a.group=='plant' else 'human-competitive-validation.json'
    schema='metag-plant-validation-v2' if a.group=='plant' else 'metag-human-competitive-validation-v1'
    if a.group=='human':a.competitor_taxid=sorted(set(a.competitor_taxid)|{9598,9595})
    out=Path(a.out or str(Path(cfg['work_dir'])/'qc'/output_group/a.sample)).absolute()
    if not prepare_result(out,replace=a.replace,cfg=cfg,protected=[a.representatives]):return
    folder=configured_split_root(cfg)/a.sample; summary,rows,_=load_sample(folder,a.sample)
    if not summary['paired']:p.error('Competitive validation requires paired input')
    selected_rows=[r for r in rows if r['group']==a.group and
                   ((r['status']=='assigned_species' and (not a.taxid or int(r['species_taxid']) in a.taxid)) or
                    (a.include_broad and r['status']=='higher_rank'))]
    total_selected=sum(int(r['fragments']) for r in selected_rows)
    if not total_selected:
        table(out/'summary.tsv',['assigned_taxid','name','fragments','tested_fragments','untested_fragments'],[])
        atomic_json(out/report_name,dict(schema=schema,status='no_candidates',target_group=a.group,
            sample=a.sample,sample_fragments=int(summary['fragments']),candidate_fragments=0,
            tested_fragments=0,untested_fragments=0,parameters=vars(a),
            interpretation='No selected saved candidates. This does not exclude reads missed by the classifier.'))
        finish_result(out);return
    tax=Taxonomy(cfg['taxonomy_dir']);count=Counter()
    for row in rows:count[int(row['species_taxid'])]+=int(row['fragments'])
    claims,reads,states,digest=extract_candidates(folder,rows,metadata_dir(folder)/'read-assignments.tsv.gz',
        out,set(a.taxid),a.include_broad,a.max_fragments,int(summary['fragments']),group=a.group)
    repfile=Path(a.representatives or str(Path(cfg['work_dir'])/'qc/reference-audit/species-representatives.tsv'))
    reps=_load_representatives(repfile,cfg['reference_raw'],require_files=False)
    budget=int(a.max_reference_gib*2**30)
    panels,assigned,untested=plan_panels(claims,count,reps,tax,a.competitor_taxid,a.max_species,budget,group=a.group)
    states.append(file_state(repfile))
    omitted=[dict(panel_id=panel['panel_id'],**row) for panel in panels for row in panel['omitted']]
    omitted.extend(dict(panel_id='',species_taxid=tid,reason=reason,bases=int(reps.get(tid,{}).get('bases',0)))
                   for tid,reason in sorted(untested.items()))
    table(out/'omitted-references.tsv',['panel_id','species_taxid','reason','bases'],omitted)
    atomic_json(out/'panel.json',dict(panels=panels,untested=untested,max_bases_per_panel=budget,
        max_species_per_panel=a.max_species,execution='serial',
        restriction='Each claim is tested in one declared panel. MAPQ is within-panel only; omitted organisms and references in other panels are not excluded.'))
    evidence={}
    for claim in claims:
        if not assigned[claim['ordinal']]:
            evidence[claim['ordinal']]={**claim,'panel_id':'','category':'untested','tested':False,
                'untested_reason':untested.get(claim['species_taxid'] or claim['assigned_taxid'],'no_competitive_panel'),
                'mapq_scope':'not_tested'}
    executable=None;failures=[];commands=[]
    if panels:
        executable=_executable(a.minimap2 or str(Path(cfg['tools_dir'])/'minimap2/2.28/minimap2'))
        states.append(file_state(executable))
    for panel in panels:
        panel_id=panel['panel_id']
        subset,pair=write_panel_reads(reads,claims,assigned,panel_id,out)
        if not subset:continue
        reference=out/f'.panel.{panel_id}.fa';sam=out/f'.panel.{panel_id}.sam'
        print(f'{a.group} validation panel {panel_id}: {len(subset):,} candidate pairs, '
              f'{len(panel["selected"])} references, {panel["declared_bases"]:,} declared bases',flush=True)
        try:
            states.extend(file_state(x['resolved_path']) for x in panel['selected'])
            contigs,actual=build_panel(panel['selected'],reference,budget)
            table(out/f'contigs.{panel_id}.tsv',list(next(iter(contigs.values()))),contigs.values())
            # One complete index per panel: multipart minimap2 MAPQ is invalid.
            command=[executable,'-ax','sr','--secondary=yes','-N','20','-I',str(actual+1000000),
                     '-t',str(a.threads),str(reference),*[str(x) for x in pair]]
            commands.append(dict(panel_id=panel_id,command=command,bases=actual,fragments=len(subset)))
            with sam.open('w') as stdout,(out/f'minimap2.{panel_id}.log').open('w') as stderr:
                subprocess.run(command,stdout=stdout,stderr=stderr,check=True)
            align=parse_sam(sam,len(subset),contigs,a.min_mapq,a.min_identity,a.min_aligned_fraction)
            for local,claim in enumerate(subset,1):
                evidence[claim['ordinal']]=evidence_row(claim,[align[(local,m)] for m in (1,2)],tax,panel_id)
        except (OSError,ValueError,RuntimeError,subprocess.CalledProcessError) as exc:
            failures.append(dict(panel_id=panel_id,error=str(exc),fragments=len(subset)))
            for claim in subset:
                evidence[claim['ordinal']]={**claim,'panel_id':panel_id,'category':'untested',
                    'tested':False,'untested_reason':'panel_failed: '+str(exc),'mapq_scope':'not_tested'}
        finally:
            if not a.keep_alignment_files:
                for path in [reference,sam,*pair]:path.unlink(missing_ok=True)
    if len(evidence)!=len(claims):raise RuntimeError('Validation did not account for every selected fragment')
    evidence=[evidence[key] for key in sorted(evidence)]
    table(out/'fragment-evidence.tsv',list(dict.fromkeys(k for x in evidence for k in x)),evidence)
    aggregates=aggregate_evidence(evidence,tax)
    fields=sorted({k for row in aggregates for k in row});table(out/'summary.tsv',fields,aggregates)
    states.extend(taxonomy_states)
    # Stat checks protect active inputs without re-hashing large references.
    if any(file_state(x['path'])!=x for x in states):raise RuntimeError('Inputs/references changed during validation')
    if code_identity()!=initial_source:raise RuntimeError('Source changed during validation')
    tested=sum(x['tested'] for x in evidence)
    validation_status='complete' if tested==len(claims) else 'partially_tested' if tested else 'not_tested'
    interpretation=('Bounded-panel compatibility, not confirmed species. Each candidate is evaluated in one '
        'explicit panel; MAPQ is within that panel and is not compared across panels. Unknown/omitted organisms '
        'remain possible. All selected candidates are reported, including untested and failed-panel records. '
        'Contig classes/repeats are header heuristics, not interval annotations. Exact-species support requires '
        'both quality-passing mates matching the claim; higher-rank compatibility stays at its original rank. '
        'Nuclear support excludes repeat headers. Locus metrics are stratified; covered bases merge overlapping '
        'reference intervals. Paired concordance requires mapper proper-pair flags and reciprocal coordinates. '
        'Candidate-only validation cannot estimate reads missed by the original classifier. Alignment uses '
        'an independent method but reuses the classification collection of reference genomes; reference '
        'contamination is therefore not excluded by agreement alone.')
    report=dict(schema=schema,status=validation_status,target_group=a.group,sample=a.sample,
        sample_fragments=int(summary['fragments']),candidate_fragments=len(claims),
        tested_fragments=tested,untested_fragments=len(claims)-tested,
        candidate_species=len({x['species_taxid'] for x in claims if x['species_taxid']}),
        higher_rank_fragments=sum(not x['species_taxid'] for x in claims),
        categories=dict(Counter(x['category'] for x in evidence)),panels=len(commands),panel_failures=failures,
        parameters=vars(a),input_digest=digest,commands=commands,
        coordinate_convention='0-based half-open reference coordinates',interpretation=interpretation)
    atomic_json(out/report_name,report)
    atomic_json(out/'complete.json',dict(validation_status=validation_status,
        sample=a.sample,fragments=len(claims),tested_fragments=tested,source_identity=initial_source,
        input_digest=digest,input_states=states,parameters=vars(a),
        interpretation=interpretation))
    finish_result(out)
    print(json.dumps(report,indent=2,sort_keys=True))
    if failures:
        marker=json.loads((out/'complete.json').read_text())
        marker['status']='incomplete';atomic_json(out/'complete.json',marker)
        raise SystemExit('Some competitive panels failed; other panels completed. See '+report_name+'.')

if __name__=='__main__':main()

