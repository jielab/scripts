"""Independent competitive validation of plant species and higher-rank candidates.

Outputs are compatibility evidence relative to an explicitly bounded reference
panel. Original signed assignments and split FASTQs are never modified.
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


def select_panel(claims, counts, reps, tax, extra, maximum, budget):
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
    groups=('bacteria','archaea','fungi','other_eukaryote','invertebrate','other_vertebrate','mammal')
    missing=[]
    for group in groups:
        ranked=sorted((t for t,r in reps.items() if r['group']==group),key=lambda t:(-counts[t],t))
        if ranked:
            for tid in ranked[:2]: wanted.setdefault(tid,group+'_competitor')
        else: missing.append(dict(species_taxid=0,reason='missing_group:'+group,bases=0))
    selected=[]; omitted=list(missing); total=0
    order=sorted(wanted,key=lambda t:({'claimed_species':0,'human_competitor':1,'explicit_competitor':1,'close_relative':2}.get(wanted[t],3),-counts[t],t))
    for tid in order:
        if tid not in reps:
            omitted.append(dict(species_taxid=tid,reason='representative_missing:'+wanted[tid],bases=0));continue
        row=reps[tid]; n=int(row['bases'])
        reason='max_species' if len(selected)>=maximum else ('reference_budget' if total+n>budget else '')
        if reason:
            omitted.append(dict(species_taxid=tid,reason=reason+':'+wanted[tid],bases=n));continue
        selected.append({**row,'role':wanted[tid]});total+=n
    return selected,omitted,total


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


def extract_candidates(folder, rows, assignment, out, taxids, broad, maximum, expected_total):
    chosen=[r for r in rows if r['group']=='plant' and
            ((r['status']=='assigned_species' and (not taxids or int(r['species_taxid']) in taxids)) or
             (broad and r['status']=='higher_rank'))]
    bins={r['bin'] for r in chosen};expected=sum(int(r['fragments']) for r in chosen)
    if not expected:raise RuntimeError('No selected plant candidates')
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
                     score=score,passed=mapped and mapq>=min_mapq and identity>=min_identity and fraction>=min_fraction,
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
    independent=0
    for intervals in regions.values():
        end=-1
        for left,right in sorted(intervals):
            # Half-open adjacent intervals are separate nonoverlapping regions.
            if left>=end:independent+=1
            end=max(end,right)
    return dict(mates=len(entries),nonoverlapping_regions=independent,loci_10kb=len(loci),
                maximum_locus_fraction=max(loci.values(),default=0)/max(1,len(entries)),
                organellar_only=bool(entries) and all(k in ('chloroplast','mitochondrial') for k in classes),
                repeat_annotated_mates=sum(x[f'r{m}_repeat_annotation'] for x,m in entries),
                unknown_class_mates=classes['unknown'])


@managed_outputs
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--replace',type=replace_bool,default=True);p.add_argument('--config');p.add_argument('--root');p.add_argument('--sample',required=True);p.add_argument('--out')
    p.add_argument('--representatives');p.add_argument('--taxid',type=int,action='append',default=[])
    p.add_argument('--competitor-taxid',type=int,action='append',default=[]);p.add_argument('--include-broad',action='store_true')
    p.add_argument('--max-species',type=int,default=60);p.add_argument('--max-reference-gib',type=float,default=24)
    p.add_argument('--max-fragments',type=int,default=1000000);p.add_argument('--threads',type=int,default=32)
    p.add_argument('--min-mapq',type=int,default=20);p.add_argument('--min-identity',type=float,default=.95)
    p.add_argument('--min-aligned-fraction',type=float,default=.8);p.add_argument('--minimap2');p.add_argument('--allow-local',action='store_true')
    a=p.parse_args(argv);require_compute(a.allow_local);cfg=load_config(a.config,a.root)
    initial_source=code_identity()
    taxonomy_states=[file_state(x) for x in sorted(Path(cfg['taxonomy_dir']).glob('*.dmp'))]
    if Path(a.sample).is_absolute() or '..' in Path(a.sample).parts:p.error('Unsafe sample ID')
    if min(a.max_species,a.max_fragments,a.threads)<1 or not math.isfinite(a.max_reference_gib) or a.max_reference_gib<=0:p.error('Limits must be positive')
    if not 0<=a.min_mapq<=255 or not 0<=a.min_identity<=1 or not 0<=a.min_aligned_fraction<=1:p.error('Invalid alignment thresholds')
    out=Path(a.out or str(Path(cfg['work_dir'])/'qc/plants'/a.sample)).absolute()
    if not prepare_result(out,replace=a.replace,cfg=cfg,protected=[a.representatives]):return
    folder=configured_split_root(cfg)/a.sample; summary,rows,_=load_sample(folder,a.sample)
    if not summary['paired']:p.error('Plant competitive validation requires paired input')
    tax=Taxonomy(cfg['taxonomy_dir']);count=Counter()
    for row in rows:count[int(row['species_taxid'])]+=int(row['fragments'])
    claims,reads,states,digest=extract_candidates(folder,rows,metadata_dir(folder)/'read-assignments.tsv.gz',out,set(a.taxid),a.include_broad,a.max_fragments,int(summary['fragments']))
    repfile=Path(a.representatives or str(Path(cfg['work_dir'])/'qc/reference-audit/species-representatives.tsv'))
    reps=_load_representatives(repfile,cfg['reference_raw']);budget=int(a.max_reference_gib*2**30)
    claimed={x['species_taxid'] for x in claims if x['species_taxid']}
    selected,omitted,bases=select_panel(claimed,count,reps,tax,a.competitor_taxid,a.max_species,budget)
    table(out/'omitted-references.tsv',['species_taxid','reason','bases'],omitted)
    atomic_json(out/'panel.json',dict(selected=selected,omitted=omitted,max_bases=budget,max_species=a.max_species,
        restriction='MAPQ is relative to this panel; omitted organisms are not excluded'))
    if not claimed.issubset({r['species_taxid'] for r in selected}):raise RuntimeError('Claimed species missing/over budget; see panel.json; no misleading validation published')
    states.extend(file_state(x['resolved_path']) for x in selected);states.append(file_state(repfile))
    panel=out/'competitive.fa';contigs,actual=build_panel(selected,panel,budget)
    table(out/'contigs.tsv',list(next(iter(contigs.values()))),contigs.values())
    executable=_executable(a.minimap2 or str(Path(cfg['tools_dir'])/'minimap2/2.28/minimap2'));states.append(file_state(executable))
    # Use one index part exceeding actual bases. Duplicate-primary checks also
    # reject unexpected multipart SAM rather than interpreting its MAPQ.
    command=[executable,'-ax','sr','--secondary=yes','-N','20','-I',str(actual+1000000),'-t',str(a.threads),str(panel),*[str(x) for x in reads]]
    sam=out/'alignment.sam'
    with sam.open('w') as stdout,(out/'minimap2.log').open('w') as stderr:
        subprocess.run(command,stdout=stdout,stderr=stderr,check=True)
    align=parse_sam(sam,len(claims),contigs,a.min_mapq,a.min_identity,a.min_aligned_fraction)
    evidence=[];groups=defaultdict(list)
    for claim in claims:
        mates=[align[(claim['ordinal'],m)] for m in (1,2)]
        result={**claim,'category':category(claim,mates,tax),'both_mates_pass':all(x['passed'] for x in mates),
                'discordant_mates':all(x['mapped'] for x in mates) and mates[0]['species_taxid']!=mates[1]['species_taxid']}
        result.update(pair_evidence(mates))
        for m,row in enumerate(mates,1):result.update({f'r{m}_{k}':v for k,v in row.items()})
        evidence.append(result);groups[claim['assigned_taxid']].append(result)
    table(out/'fragment-evidence.tsv',list(evidence[0]),evidence)
    aggregates=[]
    for tid,items in sorted(groups.items()):
        mapped=[(x,m) for x in items for m in (1,2) if x[f'r{m}_mapped']]
        support=[(x,m) for x,m in mapped if x['category'] in ('supported_same_species','higher_rank_compatible')
                 and x[f'r{m}_passed'] and x[f'r{m}_species_taxid']==tax.resolve(x['species_taxid'])]
        result=dict(assigned_taxid=tid,fragments=len(items),mates=len(items)*2,
                    all_mapped_mates=len(mapped),all_mapped_fraction=len(mapped)/(2*len(items)),
                    claimed_support_mates=len(support),claimed_support_fraction=len(support)/(2*len(items)),
                    both_mates_pass=sum(x['both_mates_pass'] for x in items),
                    cross_species_mates=sum(x['cross_species_mates'] for x in items),
                    concordant_support_fragments=sum(x['concordant_pair'] and x['category']=='supported_same_species' for x in items),
                    same_species_only_fragments=sum(x['same_species_mates'] and not x['concordant_pair'] for x in items),
                    one_mate_pass_fragments=sum(sum(x[f'r{m}_passed'] for m in (1,2))==1 for x in items),
                    low_complexity_fragments=sum(x['low_complexity'] for x in items))
        strata={'all_mapped':mapped,'claimed_support':support}
        for cls in ('nuclear','organellar','repeat','unknown'):
            strata['claimed_'+cls]=[(x,m) for x,m in support if support_class(x,m)==cls]
        for cat in ('supported_same_species','higher_rank_compatible','ambiguous','conflicting_competitor','unmapped'):
            strata[cat]=[(x,m) for x,m in mapped if x['category']==cat]
        strata['one_mate_pass']=[(x,m) for x,m in mapped if sum(x[f'r{j}_passed'] for j in (1,2))==1 and x[f'r{m}_passed']]
        for name,entries in strata.items():
            result.update({name+'_'+key:value for key,value in locus_stats(entries).items()})
        result.update(Counter(x['category'] for x in items));aggregates.append(result)
    fields=sorted({k for row in aggregates for k in row});table(out/'summary.tsv',fields,aggregates)
    states.extend(taxonomy_states)
    if any(file_state(x['path'])!=x for x in states):raise RuntimeError('Inputs/references changed during validation')
    if code_identity()!=initial_source:raise RuntimeError('Source changed during validation')
    atomic_json(out/'complete.json',dict(status='complete',sample=a.sample,fragments=len(claims),mates=len(claims)*2,
        source_identity=initial_source,source_release=str(Path(__file__).resolve().parent.parent),input_digest=digest,input_states=states,command=command,parameters=vars(a),
        panel_bases=actual,omitted=omitted,outputs={p.name:sha256(p) for p in (out/'summary.tsv',out/'fragment-evidence.tsv',out/'contigs.tsv',out/'panel.json')},
        coordinate_convention='0-based half-open reference coordinates; CIGAR M/D/N/=/X consume reference',
        interpretation='Bounded-panel compatibility, not confirmed species. Unknown/omitted organisms remain possible. Contig classes/repeats are header heuristics, not interval annotations. Claimed support requires both quality-passing compatible mates and matching claimed species; nuclear excludes repeat headers. Strata are organellar, repeat, nuclear, unknown in that precedence. Locus fractions use mapped mates within each named stratum; overall fractions use all candidate mates. Concordance uses mapper proper-pair flags plus same contig and reciprocal mate coordinates. Cross-contig same-species hits retain compatibility with pairing limitations. Region counts are not organism abundances.'))

    finish_result(out)

if __name__=='__main__':main()
