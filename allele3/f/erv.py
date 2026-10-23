"""UCSC joined repeats: aligned bodies, strand-aware observed element edges."""
from __future__ import annotations
from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path
import gzip, json, os, fcntl
try:
    from .common import canonical_chrom, open_text, chrom_sort_key
    from .core import fingerprint, dump_json
except ImportError:
    from common import canonical_chrom, open_text, chrom_sort_key
    from core import fingerprint, dump_json
CORE_FAMILIES={'ERV1','ERVK','ERVL'}
RELATED_FAMILIES=CORE_FAMILIES | {'ERVL-MaLR','Gypsy'}
PRIMARY={str(i) for i in range(1,23)}|{'X','Y'}

@dataclass(frozen=True, slots=True)
class ERV:
    chrom:str
    start:int
    end:int
    name:str
    rep_class:str
    family:str
    strand:str
    source_id:str=''
    segments:tuple=()
    @property
    def tier(self):
        if self.rep_class!='LTR': return 'non_ltr'
        if self.family in CORE_FAMILIES: return 'erv_core'
        if self.family in RELATED_FAMILIES: return 'erv_related'
        return 'ltr_other'
    @property
    def five_prime(self): return self.start if self.strand=='+' else self.end
    @property
    def three_prime(self): return self.end if self.strand=='+' else self.start

def parse_joined_rmsk(path, include_alt_contigs=False, chromosome=None):
    """Aligned blocks define observed joined edges, excluding consensus tails.

    Only blocks whose relative start is not the -1 sentinel are aligned DNA; nested-insertion gaps do
    not count toward repeat bodies. Retain all classes for negative controls.
    """
    out=[]
    with open_text(path) as f:
        for ln,line in enumerate(f,1):
            if not line.strip() or line.startswith('#'): continue
            a=line.rstrip('\r\n').split('\t') if '\t' in line else line.split()
            if len(a)!=14: raise ValueError(f'{path}:{ln}: expected 14 joined RepeatMasker columns')
            chrom=canonical_chrom(a[1])
            if chromosome and chrom!=canonical_chrom(chromosome): continue
            if not include_alt_contigs and chrom not in PRIMARY: continue
            outer_start,outer_end=int(a[2]),int(a[3])
            start,end=int(a[7]),int(a[8]); strand=a[6]
            if not (0<=outer_start<outer_end and 0<=start<end) or strand not in {'+','-'}:
                raise ValueError(f'{path}:{ln}: invalid joined coordinates/strand')
            sizes=[int(x) for x in a[11].rstrip(',').split(',')]
            offsets=[int(x) for x in a[12].rstrip(',').split(',')]
            if len(sizes)!=int(a[10]) or len(offsets)!=len(sizes): raise ValueError(f'{path}:{ln}: invalid joined block cardinality')
            segments=[]
            for size,offset in zip(sizes,offsets):
                if offset != -1 and size<0: raise ValueError('Invalid aligned joined block')
                # Unaligned consensus gaps can have negative sizes (overlap); ignore them.
                if offset!=-1 and size>0:
                    lo=outer_start+offset; hi=lo+size
                    if lo<0: raise ValueError(f'{path}:{ln}: negative aligned genomic coordinate')
                    segments.append((lo,hi))
            if not segments: raise ValueError(f'{path}:{ln}: no aligned repeat sequence')
            # A few UCSC records have overlapping, out-of-order blocks: alignEnd
            # is the final block's end, not the maximum observed aligned end.
            # Use the envelope of ALL aligned blocks, preserving every body base.
            start=min(s for s,t in segments); end=max(t for s,t in segments)
            name,sep,cf=a[4].rpartition('#')
            if not sep: raise ValueError(f'{path}:{ln}: repeat class missing from joined name')
            rep_class,_,family=cf.partition('/')
            out.append(ERV(chrom,start,end,name,rep_class,family,strand,a[13],tuple(segments)))
    out.sort(key=lambda e:(chrom_sort_key(e.chrom),e.start,e.end,e.name))
    return out

def write_erv_bed(elements,out_path):
    Path(out_path).parent.mkdir(parents=True,exist_ok=True)
    with open(out_path,'w') as f:
        f.write('chrom\tstart\tend\tname\trep_class\tfamily\tstrand\ttier\tsource_id\tsegments\n')
        for e in elements:
            f.write('\t'.join(map(str,('chr'+e.chrom,e.start,e.end,e.name,e.rep_class,e.family,e.strand,e.tier,e.source_id,','.join(f'{s}-{t}' for s,t in e.segments))))+'\n')

class ERVOverlapSweep:
    def __init__(self,elements):
        self.by_chrom=defaultdict(list)
        for e in elements: self.by_chrom[e.chrom].append(e)
        for a in self.by_chrom.values(): a.sort(key=lambda e:e.start)
        self.chrom=None; self.i=0; self.active=[]; self.last={}
    def locate(self,chrom,pos1):
        chrom=canonical_chrom(chrom); pos0=int(pos1)-1
        if pos0<self.last.get(chrom,-1): raise ValueError('Repeat sweep requires sorted positions')
        self.last[chrom]=pos0
        arr=self.by_chrom.get(chrom,[])
        if chrom!=self.chrom: self.chrom=chrom; self.i=0; self.active=[]
        self.active=[e for e in self.active if e.end>pos0]
        while self.i<len(arr) and arr[self.i].start<=pos0:
            e=arr[self.i]
            if e.end>pos0: self.active.append(e)
            self.i+=1
        return tuple(e for e in self.active if any(lo<=pos0<hi for lo,hi in (e.segments or ((e.start,e.end),))))

class BoundarySweep:
    def __init__(self,elements,tiers=('erv_core','erv_related')):
        self.points=defaultdict(lambda:{'5p':[],'3p':[]}); allowed=set(tiers)
        for e in elements:
            if e.tier not in allowed: continue
            self.points[e.chrom]['5p'].append((e.five_prime,e))
            self.points[e.chrom]['3p'].append((e.three_prime,e))
        for v in self.points.values():
            for a in v.values(): a.sort(key=lambda x:(x[0],x[1].name,x[1].source_id))
        self.state={}; self.last={}
    def nearest(self,chrom,pos1,boundary_type,max_distance=10000):
        chrom=canonical_chrom(chrom); p=int(pos1)-1; key=(chrom,boundary_type)
        if boundary_type not in {'5p','3p'}: raise ValueError('Boundary must be 5p or 3p')
        if p<self.last.get(key,-1): raise ValueError('Boundary sweep requires sorted sites')
        self.last[key]=p
        a=self.points.get(chrom,{}).get(boundary_type,[])
        if not a: return None
        i=self.state.get(key,0)
        while i+1<len(a) and a[i+1][0]<=p: i+=1
        self.state[key]=i
        point,e=min(a[i:i+2],key=lambda x:(abs(p-x[0]),x[0],x[1].name))
        distance=abs(p-point)
        if distance>=max_distance: return None
        signed=(p-point)*(1 if e.strand=='+' else -1)
        return dict(distance=distance,signed_genomic_distance=p-point,signed_transcription_distance=signed,boundary_pos0=point,element=e)

def boundary_bin(distance):
    if distance is None: return 'gt10kb'
    lo=0
    for hi in (100,200,500,1000,2000,5000,10000):
        if abs(int(distance))<hi: return f'{lo}-{hi}'
        lo=hi
    return 'ge10000'

def prepare_joined_cache(path, directory):
    """Split once, under a lock, so chromosome workers do not reread 5M repeats."""
    d=Path(directory); d.mkdir(parents=True,exist_ok=True)
    with (d/'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        source=fingerprint(path); marker=d/'manifest.json'
        if marker.exists():
            old=json.loads(marker.read_text())
            if old['source']!=source: raise ValueError('Joined annotation changed; choose a new annotation cache')
            for name,fp in old['files'].items():
                if fingerprint(d/name)!=fp: raise ValueError('Joined annotation cache changed')
            return old
        handles={}; counts=defaultdict(int)
        try:
            with open_text(path) as f:
                for line in f:
                    if not line.strip() or line.startswith('#'): continue
                    a=line.split('\t',2)
                    if len(a)<3: raise ValueError('Malformed joined annotation')
                    c=canonical_chrom(a[1])
                    if c not in PRIMARY: continue
                    if c not in handles: handles[c]=gzip.open(d/f'chr{c}.tsv.gz.partial','wt',compresslevel=2)
                    handles[c].write(line); counts[c]+=1
        finally:
            for f in handles.values(): f.close()
        for c in PRIMARY:
            if c not in handles:
                with gzip.open(d/f'chr{c}.tsv.gz.partial','wt'): pass
            os.replace(d/f'chr{c}.tsv.gz.partial',d/f'chr{c}.tsv.gz')
        result=dict(schema=1,source=source,counts=dict(counts),files={p.name:fingerprint(p) for p in d.glob('chr*.tsv.gz')})
        dump_json(marker,result); return result
