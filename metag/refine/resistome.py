"""Community-level AMR family support from curated ShortBRED true markers and DIAMOND."""
from collections import Counter, defaultdict
import csv
import re
import subprocess
from common import sha256
from fastq_stream import fragments, new_digest, update_digest

MARKER = re.compile(r'(.+)_(TM|JM|QM)\d*_#\d+')

def true_markers(source, path):
    seen = set()
    family_by_index = ['']
    families = Counter()
    def emit(ident, bits, out):
        if ident is None:
            return
        match = MARKER.fullmatch(ident.split('__[')[0])
        if not match:
            raise ValueError('Expected a ShortBRED marker ID, found '+ident)
        if match[2] != 'TM':
            return
        seq = ''.join(bits)
        if not seq:
            raise ValueError('Empty true marker '+ident)
        if ident in seen:
            raise ValueError('Duplicate true marker '+ident)
        seen.add(ident)
        if len(seq) < 25:
            return
        family_by_index.append(match[1])
        families[match[1]] += 1
        out.write(f'>m{len(family_by_index)-1}\n{seq}\n')
    ident, bits = None, []
    with source.open() as src, path.open('w') as out:
        for line in src:
            if line.startswith('>'):
                emit(ident,bits,out)
                ident, bits = line[1:].split()[0], []
            else:
                bits.append(line.strip())
        emit(ident,bits,out)
    if not families:
        raise ValueError('No true markers of at least 25 amino acids in the selected AMR reference')
    return family_by_index,families

def annotations(path):
    if not path:
        return {}
    labels = defaultdict(set)
    with path.open(newline='',encoding='utf-8') as f:
        rows = csv.DictReader(f,delimiter='\t')
        if not {'Family','ID_Orig','Merged.ID'} <= set(rows.fieldnames or []):
            raise ValueError('Annotation TSV requires Family, ID_Orig, Merged.ID')
        for r in rows:
            labels[r['Family']].add((r['ID_Orig'],r['Merged.ID']))
    return {k:tuple('|'.join(sorted({v[i] for v in a})) for i in (0,1)) for k,a in labels.items()}

def count_hits(path, family_by_index, n):
    # DIAMOND emits query-ordered output. Queries were written 1_1, 1_2, 2_1 ...
    # Keep only one physical fragment in memory, irrespective of dataset size.
    counts=Counter()
    both_mates=Counter()
    unique_markers=defaultdict(set)
    ambiguous=0
    previous=0
    scores=defaultdict(dict)
    def finish():
        nonlocal ambiguous
        possible=set()
        distinct=set()
        for marker_scores in scores.values():
            best=max(marker_scores.values())
            near=[m for m,score in marker_scores.items() if score>=.98*best]
            possible.update(family_by_index[m] for m in near)
            if len(near)==1:
                distinct.add(near[0])
        if len(possible)==1:
            fam=next(iter(possible))
            counts[fam]+=1
            both_mates[fam]+=len(scores)==2
            unique_markers[fam].update(distinct)
        elif possible:
            ambiguous+=1
    with path.open() as handle:
        for line in handle:
            q,subject,ident,aa,query_bp,marker_aa,score=line.split()
            ordinal,end = map(int,q.split('_'))
            if ordinal < previous or not 1 <= ordinal <= n or end not in (1,2):
                raise ValueError('Unexpected translated-search query order or ID')
            if ordinal != previous:
                finish()
                scores.clear()
                previous=ordinal
            if not subject.startswith('m') or not 0<int(subject[1:])<len(family_by_index):
                raise ValueError('Unexpected translated-search subject ID')
            aligned=int(aa)
            if float(ident)<95 or aligned<25 or aligned<.95*min(int(query_bp)//3,int(marker_aa)):
                continue
            marker_id=int(subject[1:])
            scores[end][marker_id]=max(float(score),scores[end].get(marker_id,0))
        finish()
    return counts,ambiguous,both_mates,unique_markers

def run(sample, receipt, markers, cores, tmp, annotation=None):
    query=tmp/'community.fa'
    digest=new_digest()
    n=0
    with query.open('wb') as f:
        for n,(left,right) in enumerate(fragments(*(r[0] for r in receipt['request']['raw'])),1):
            update_digest(digest,left,right)
            for end,record in enumerate((left,right),1):
                f.write(f'>{n}_{end}\n'.encode()+record[1].strip()+b'\n')
    if n!=receipt['pairs'] or digest.hexdigest()!=receipt['input_digest']:
        raise ValueError('Original FASTQ differs from split: '+sample)
    true=tmp/'true.faa'
    by_index,families=true_markers(markers,true)
    labels=annotations(annotation)
    db=tmp/'markers.dmnd'
    out=tmp/'hits.tsv'
    version=subprocess.check_output(['diamond','version'],text=True).strip()
    subprocess.run(['diamond','makedb','--in',str(true),'--db',str(db),'--threads',str(cores)],check=True)
    counts,ambiguous,both_mates,unique_markers=Counter(),0,Counter(),defaultdict(set)
    if n:
        subprocess.run(['diamond','blastx','--db',str(db),'--query',str(query),'--threads',str(cores),
                        '--sensitive','--max-target-seqs','0','--out',str(out),'--outfmt','6',
                        'qseqid','sseqid','pident','length','qlen','slen','bitscore'],check=True)
        counts,ambiguous,both_mates,unique_markers=count_hits(out,by_index,n)
    calls=[[sample,fam,*labels.get(fam,('NA','NA')),families[fam],len(unique_markers[fam]),counts[fam],both_mates[fam],n,
            1e6*counts[fam]/n if n else 'NA'] for fam in sorted(families)]
    return calls,dict(engine=version,marker_set=str(markers),tested_families=len(families),
                      ambiguous_pairs=ambiguous,assigned_pairs=sum(counts.values()),
                      annotation_sha256=sha256(annotation) if annotation else None)
