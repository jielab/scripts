"""Read-pair competition, nuclear placement and cross-locus support in the two Weddell bins."""
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory

from common import atomic_json, fingerprint, locked, sha256
from fastq_stream import fragments, new_digest, update_digest
from ranked_storage import publish_directory, recover_directory, staging
from metag_refine_tables import rows, write
import seal_panel as panel

HERE=Path(__file__).resolve().parent
BINS=(panel.TARGET,panel.GENUS)
OUTCOMES=('target_nuclear','target_other_contig','target_nonconcordant','target_multilocus',
          'low_complexity','shared_phocid','shared_or_conflicting','alternative_taxon','weak_or_one_mate','hit_limit','panel_limited')
SUMMARY=('sample bin_taxid bin_pairs tested_pairs panel_id panel_scope '+ ' '.join(OUTCOMES)+
         ' unique_nuclear_fragments nuclear_windows_10kb separated_nuclear_loci_10kb most_populated_window_pairs interpretation').split()
LOCI='sample bin_taxid accession contig_class window_10kb_start zero_based_end_exclusive pairs distinct_fragment_spans'.split()
PAIR_COLS=('bin_taxid bin_ordinal outcome mate1_near_taxids mate2_near_taxids target_accession '
           'target_mate1_start1 target_mate1_end1 target_mate2_start1 target_mate2_end1 '
           'mate1_top_bits mate2_top_bits mate1_best_other_bits mate2_best_other_bits').split()
FILES=('seal_bins.tsv','seal_loci.tsv','seal_pairs.tsv.gz','seal_hits.tsv.gz')
LIMIT=64
HSPS=8

def location(root,sample):
    return root/'seal'/sample

def tool_version():
    return subprocess.check_output(['blastn','-version'],text=True).splitlines()[0]

def audit_key(folder):
    ident,_=panel.signature(folder)
    return fingerprint([ident,sha256(HERE/'seal_audit.py'),sha256(HERE/'seal_panel.py'),tool_version()])

def current(root,sample,split_sha,opts):
    opt=opts.get('seal')
    if not opt:
        return True
    dest=location(root,sample)
    if opts.get('replace') or not (dest/'receipt.json').is_file() or not Path(opt['path']).exists():
        return False
    with locked(dest.parent/('.'+dest.name+'.lock'),shared=True):
        return matches(dest,sample,split_sha,audit_key(Path(opt['path'])))

def matches(dest,sample,split_sha,key):
    p=dest/'receipt.json'
    if not p.is_file():
        return False
    r=json.loads(p.read_text())
    return (r.get('sample')==sample and r.get('split_sha256')==split_sha and r.get('audit_key')==key
            and all((dest/f).is_file() and (dest/f).stat().st_size==r.get('file_sizes',{}).get(f) for f in FILES))

def contig_class(title):
    s=title.casefold()
    if any(k in s for k in ('mitochond','plastid','chloroplast','ribosomal rna','unlocalized','unplaced','repeat','vector','plasmid','alternate locus')):
        return 'organelle_repeat_or_unplaced'
    if re.search(r'\bchromosome\s+[\w.-]+',s) and 'scaffold' not in s and 'contig' not in s:
        return 'nuclear_chromosome'
    return 'other_or_untyped'

@dataclass(frozen=True)
class Hit:
    sid:str
    score:float
    start:int
    end:int
    @property
    def tid(self):
        return int(self.sid.split('_')[1])

def blast_block(pairs,db,tmp,threads,raw):
    query=tmp/'queries.fa'
    out=tmp/'hits.tsv'
    with query.open('wb') as f:
        for i,(tid,pos,a,b) in enumerate(pairs):
            for m,r in enumerate((a,b),1):
                f.write(f'>{i}_{m}\n'.encode()+r[1].strip()+b'\n')
    fields='qseqid sseqid bitscore pident qcovhsp qstart qend sstart send'
    subprocess.run(['blastn','-task','blastn','-query',str(query),'-db',str(db/'competition'),
                    '-dust','yes','-num_threads',str(threads),'-max_target_seqs',str(LIMIT),
                    '-max_hsps',str(HSPS),'-outfmt','6 '+fields,'-out',str(out)],check=True)
    hits=defaultdict(list)
    counts=defaultdict(Counter)
    with out.open() as f:
        for fields in csv.reader(f,delimiter='\t'):
            q,subject,bs,pid,cov,qs,qe,ss,se=fields
            i,m=map(int,q.split('_'));sid=subject.removeprefix('lcl|')
            if not 0<=i<len(pairs) or m not in (1,2):
                raise ValueError('Bad BLAST query ordinal')
            tid,pos,*_=pairs[i]
            raw.writerow([tid,pos,m,sid,*fields[2:]])
            counts[i,m][sid]+=1
            if float(pid)>=95 and float(cov)>=80:
                start,end=int(ss),int(se)
                if int(qs)>int(qe):
                    start,end=end,start
                hits[i,m].append(Hit(sid,float(bs),start,end))
    return hits,counts

def near(hits):
    if not hits:return []
    best=max(h.score for h in hits)
    return list({h for h in hits if h.score>=best*.98})

def simple_sequence(sequence):
    seq=sequence.strip().upper()
    bases=Counter(seq)
    canonical=sum(bases[b] for b in b'ACGT')
    return (not seq or canonical<.8*len(seq) or max(bases.values())>.8*len(seq)
            or len({seq[j:j+3] for j in range(len(seq)-2)})<5)

def competition(left,right,truncated,headers,scope,low_complexity=False):
    candidates=[near(x) for x in (left,right)]
    ids=[set(h.tid for h in end) for end in candidates]
    ids_text=['|'.join(str(t) for t in sorted(v)) or 'NA' for v in ids]
    top=[max((h.score for h in end),default='NA') for end in (left,right)]
    rivals=[max((h.score for h in end if h.tid!=panel.TARGET),default='NA') for end in (left,right)]
    coords=['NA']*5
    span=None
    if low_complexity:
        outcome='low_complexity'
    elif truncated:
        outcome='hit_limit'
    elif not all(ids):
        outcome='weak_or_one_mate'
    elif ids[0]==ids[1]=={panel.TARGET}:
        if any(len(end)!=1 for end in candidates):
            outcome='target_multilocus'
        else:
            a,b=(end[0] for end in candidates)
            same=a.sid==b.sid
            acc,title=headers[a.sid]
            coords=[acc if same else acc+'|'+headers[b.sid][0],a.start,a.end,b.start,b.end]
            lo=min(a.start,a.end,b.start,b.end);hi=max(a.start,a.end,b.start,b.end)
            plus,minus=(a,b) if a.start<a.end else (b,a)
            inward=plus.start<plus.end and minus.start>minus.end and plus.start<=minus.start
            if not same or not inward or hi-lo+1>1000:
                outcome='target_nonconcordant'
            else:
                kind=contig_class(title)
                outcome='target_nuclear' if kind=='nuclear_chromosome' else 'target_other_contig'
                span=(acc,kind,lo,hi)
    elif ids[0]==ids[1] and len(ids[0])==1:
        outcome='alternative_taxon'
    elif (ids[0]|ids[1]) <= {panel.TARGET,*scope['relative_species']}:
        outcome='shared_phocid'
    else:
        outcome='shared_or_conflicting'
    if not scope['adequate'] and outcome=='target_nuclear':
        outcome='panel_limited'
    return outcome,[*ids_text,*coords,*top,*rivals],span

def analyze(folder,out,published,db,scope,ident,threads):
    s=published['sample']
    taxa={int(r['taxid']):r for r in rows(folder/'taxa.tsv') if int(r['taxid']) in BINS}
    header_path=db/'contigs.tsv.gz'
    # Only target headers are needed to label genomic coordinates.
    headers={}
    with gzip.open(header_path,'rt',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if int(r['taxid'])==panel.TARGET:
                headers[r['subject']]=(r['accession'],r['description'])
    counts={tid:Counter() for tid in BINS}
    windows=defaultdict(Counter)
    expected={t:int(taxa[t]['pairs']) if t in taxa else 0 for t in BINS}
    with TemporaryDirectory(dir=out) as wd, gzip.open(out/'seal_pairs.tsv.gz','wt',newline='') as detail, gzip.open(out/'seal_hits.tsv.gz','wt',newline='') as raw:
        tmp=Path(wd)
        pairs=csv.writer(detail,delimiter='\t',lineterminator='\n');pairs.writerow(PAIR_COLS)
        hits_writer=csv.writer(raw,delimiter='\t',lineterminator='\n')
        hits_writer.writerow('bin_taxid bin_ordinal mate subject bitscore pident qcovhsp qstart qend sstart send'.split())
        def process(block):
            hits,n=blast_block(block,db,tmp,threads,hits_writer)
            for i,(tid,pos,left,right) in enumerate(block):
                cut=any(len(n[i,m])>=LIMIT or any(v>=HSPS for v in n[i,m].values()) for m in (1,2))
                outcome,fields,span=competition(hits[i,1],hits[i,2],cut,headers,scope,
                                                simple_sequence(left[1]) or simple_sequence(right[1]))
                counts[tid][outcome]+=1
                pairs.writerow([tid,pos,outcome,*fields])
                if span:
                    acc,kind,lo,hi=span
                    win=((lo+hi)//2-1)//10000*10000
                    windows[tid,acc,kind,win][lo,hi]+=1
        block=[]
        for tid,row in taxa.items():
            rel=Path(row['bin'])
            if rel.is_absolute() or '..' in rel.parts:
                raise ValueError('Unexpected bin path: '+str(rel))
            base=folder/rel;digest=new_digest();pos=0
            for pos,(a,b) in enumerate(fragments(base/'R1.fastq.gz',base/'R2.fastq.gz'),1):
                update_digest(digest,a,b)
                block.append((tid,pos,a,b))
                if len(block)>=5000:
                    process(block);block=[]
            if pos!=expected[tid] or digest.hexdigest()!=row['pairs_sha256']:
                raise ValueError('Input bin count or digest mismatch: '+str(base))
        if block:
            process(block)
    label='adequate_local_competition' if scope['adequate'] else 'missing_'+','.join(k for k,v in scope['components'].items() if not v)
    summaries=[]
    for tid in BINS:
        c=counts[tid];n=sum(c.values())
        if n!=expected[tid]:raise ValueError('Audit denominator mismatch')
        nuclear=[spans for (t,_,kind,_),spans in windows.items() if t==tid and kind=='nuclear_chromosome']
        distinct=sum(len(x) for x in nuclear)
        by_contig=defaultdict(set)
        for (t,acc,kind,_),spans in windows.items():
            if t==tid and kind=='nuclear_chromosome':
                by_contig[acc].update((lo+hi)//2 for lo,hi in spans)
        separated=0
        for positions in by_contig.values():
            last=None
            for pos in sorted(positions):
                if last is None or pos-last>=10000:
                    separated+=1
                    last=pos
        message=('no_bin' if not n else 'panel_limited' if not scope['adequate'] else
                 'distributed_panel_support' if separated>=3 else 'localized_panel_support' if nuclear else 'no_unique_nuclear_support')
        summaries.append([s,tid,n,n,ident,label,*[c[k] for k in OUTCOMES],distinct,len(nuclear),separated,
                          max((sum(x.values()) for x in nuclear),default=0),message])
    write(out/'seal_bins.tsv',SUMMARY,summaries)
    write(out/'seal_loci.tsv',LOCI,([s,t,a,k,w,w+10000,sum(spans.values()),len(spans)]
                                  for (t,a,k,w),spans in sorted(windows.items())))

def review(folder,root,opts,expected=None):
    from rank_paper import read_sample
    opt=opts['seal'];db=panel.prepare(opt);key=audit_key(db)
    ident,meta=panel.signature(db)
    sample=folder.parent.name+'/'+folder.name;dest=location(root,sample)
    with read_sample(folder),locked(dest.parent/('.'+dest.name+'.lock')):
        src=folder/'summary.json';split_sha=sha256(src);s=json.loads(src.read_text())
        if s['sample']!=sample or (expected and expected!=split_sha):
            raise ValueError('Split changed before Weddell seal audit: '+sample)
        recover_directory(dest)
        if dest.is_symlink() or (dest.exists() and not (dest/'receipt.json').is_file()):
            raise ValueError('Unexpected seal audit destination: '+str(dest))
        if not opts.get('replace') and matches(dest,sample,split_sha,key):
            return
        if sha256(db/'contigs.tsv.gz')!=meta['headers_sha256']:
            raise ValueError('Changed seal panel contig manifest')
        with staging(dest.parent,'.building-') as out:
            analyze(folder,out,s,db,meta['scope'],ident,opts['cores'])
            if key!=audit_key(db):
                raise ValueError('Seal competition panel changed during audit')
            atomic_json(out/'receipt.json',dict(sample=sample,split_sha256=split_sha,audit_key=key,
                        panel_id=ident,panel_path=str(db),scope=meta['scope'],blast_version=tool_version(),
                        file_sizes={f:(out/f).stat().st_size for f in FILES},
                        table_hashes={f:sha256(out/f) for f in ('seal_bins.tsv','seal_loci.tsv')}))
            publish_directory(out,dest)
    print('Published Weddell seal audit '+str(dest),flush=True)

def collect(root,sample,split_sha,writers):
    dest=location(root,sample)
    def status(value,panel_id='NA'):
        writers['seal_status'].writerow([sample,value,panel_id])
    lock=dest.parent/('.'+dest.name+'.lock')
    if not lock.is_file():
        status('not_run');return
    with locked(lock,shared=True):
        p=dest/'receipt.json'
        if not p.is_file():
            status('not_run');return
        rec=json.loads(p.read_text())
        if rec.get('split_sha256')!=split_sha or rec.get('sample')!=sample:
            status('split_changed');return
        for name in ('seal_bins','seal_loci'):
            file=dest/(name+'.tsv')
            if rec['table_hashes'][file.name]!=sha256(file):
                raise ValueError('Seal audit table differs from receipt: '+str(file))
            for r in rows(file):
                if r['sample']!=sample:raise ValueError('Unexpected seal audit sample')
                writers[name].writerow([r[x] for x in (SUMMARY if name=='seal_bins' else LOCI)])
        status('complete',rec['panel_id'])

def summarize_windows(folder):
    """Count samples sharing each locus within one frozen panel and original bin."""
    status={r['sample']:r for r in rows(folder/'seal_status.tsv')}
    audited=Counter()
    occupied=defaultdict(Counter)
    for sample,r in status.items():
        if r['state']=='complete':
            audited[sample.split('/')[0],r['panel_id']]+=1
    for r in rows(folder/'seal_loci.tsv'):
        s=r['sample'];pid=status[s]['panel_id']
        k=(s.split('/')[0],pid,r['bin_taxid'],r['accession'],r['contig_class'],int(r['window_10kb_start']))
        if s in occupied[k]:
            raise ValueError('Duplicate per-sample window record')
        occupied[k][s]=int(r['pairs'])
    fields='dataset panel_id bin_taxid accession contig_class window_10kb_start samples_with_window audited_samples pairs top_sample top_sample_pairs'.split()
    result=[]
    for k,counts in occupied.items():
        top=max(sorted(counts),key=counts.get)
        result.append([*k,len(counts),audited[k[:2]],sum(counts.values()),top,counts[top]])
    write(folder/'seal_recurrent_loci.tsv',fields,sorted(result,key=lambda r:(r[0],-r[6],-r[8],r[3],r[5])))
