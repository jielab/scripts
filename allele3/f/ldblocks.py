from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
try:
    from .common import canonical_chrom, chrom_sort_key
except ImportError:
    from common import canonical_chrom, chrom_sort_key

@dataclass(frozen=True)
class Block:
    chrom: str
    start: int
    end: int
    ancestry: str
    build: str
    index: int

    @property
    def block_id(self):
        return f"{self.ancestry}.{self.build}.{self.chrom}.{self.index:04d}"

def read_blocks(path, ancestry, build):
    if str(build) not in {"37", "38"}:
        raise ValueError("LD build must be 37 or 38")
    rows=[]
    with open(path, "rt", encoding="utf-8") as f:
        for line_no,line in enumerate(f,1):
            if not line.strip() or line.startswith("#"): continue
            p=line.rstrip("\n").split()
            if len(p)<3:
                raise ValueError(f"{path}:{line_no}: expected >=3 columns")
            if p[:3] in (["chr","start","end"],["chrom","start","end"],["chromosome","start","end"]): continue
            chrom=canonical_chrom(p[0]); start=int(p[1]); end=int(p[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: invalid BED interval")
            rows.append((chrom,start,end))
    rows.sort(key=lambda x:(chrom_sort_key(x[0]),x[1],x[2]))
    out=[]
    per=defaultdict(int)
    prev={}
    stats={"n":0,"overlaps":0,"touching":0,"gaps":0}
    for chrom,start,end in rows:
        per[chrom]+=1
        b=Block(chrom,start,end,ancestry,str(build),per[chrom])
        if chrom in prev:
            pe=prev[chrom]
            if start < pe: raise ValueError(f"Overlapping BED intervals in {path}: {chrom}:{start}-{end}")
            elif start == pe: stats["touching"]+=1
            else: stats["gaps"]+=1
        prev[chrom]=max(prev.get(chrom,0),end)
        out.append(b)
        stats["n"]+=1
    return out,stats

class BlockSweep:
    """O(1)-amortized assignment for sites supplied in sorted chromosome/position order."""
    def __init__(self, blocks):
        self.by_chrom=defaultdict(list)
        for b in blocks: self.by_chrom[b.chrom].append(b)
        self._chrom=None
        self._i=0
        self._last={}

    def locate(self, chrom, pos1):
        chrom=canonical_chrom(chrom)
        pos0=int(pos1)-1
        if pos0 < 0 or pos0 < self._last.get(chrom, -1):
            raise ValueError("LD sweep requires nondecreasing positive positions per chromosome")
        self._last[chrom]=pos0
        arr=self.by_chrom.get(chrom,[])
        if chrom != self._chrom:
            self._chrom=chrom; self._i=0
        while self._i < len(arr) and arr[self._i].end <= pos0:
            self._i += 1
        if self._i < len(arr):
            b=arr[self._i]
            if b.start <= pos0 < b.end:
                return b
        return None

def validate_directory(directory, build="38", ancestries=("AFR","EAS","EUR","SAS")):
    d=Path(directory)
    result={}
    for anc in ancestries:
        p=d/f"{anc}.{build}.bed"
        if not p.exists():
            result[anc]={"path":str(p),"missing":True}
            continue
        blocks,stats=read_blocks(p,anc,build)
        result[anc]={"path":str(p),"missing":False,**stats}
    return result


def load_block_sets(directory, build="38"):
    result={}
    for anc in ("AFR","EAS","EUR","SAS"):
        path=Path(directory)/f"{anc}.{build}.bed"
        blocks,_=read_blocks(path,anc,build)
        if not blocks: raise ValueError(f"Empty LD partition: {path}")
        result[anc]=blocks
    return result

class BlockIndex:
    """Random-access assignment for candidate tables which are not coordinate sorted."""
    def __init__(self, blocks):
        from bisect import bisect_right
        self.bisect=bisect_right
        self.rows=defaultdict(list)
        for b in blocks: self.rows[b.chrom].append(b)
        for c in self.rows: self.rows[c].sort(key=lambda b:b.start)
        self.starts={c:[b.start for b in a] for c,a in self.rows.items()}
    def locate(self, chrom, pos1):
        c=canonical_chrom(chrom); p=int(pos1)-1
        i=self.bisect(self.starts.get(c,[]),p)-1
        if i>=0 and p<self.rows[c][i].end: return self.rows[c][i]
        return None

def hotspot_consensus(rows, min_sites=10000, top_fraction=.01):
    """Intersect independently ranked partitions; never equate their block IDs."""
    import math
    groups=defaultdict(list)
    for r in rows:
        if r.get("sensitivity","all") != "all" or str(r.get("threshold",1)) != "1" or r.get("population_mode","all") != "all": continue
        if r["block_id"]=="outside_partition" or int(r["n_sites"])<min_sites or r.get("oe") in (None,"NA"): continue
        if not math.isfinite(float(r["oe"])): continue
        groups[r["ancestry"]].append(r)
    events=defaultdict(list)
    for anc, group in groups.items():
        group.sort(key=lambda r:(-float(r['oe']),r['chrom'],int(r['start'])))
        for r in group[:max(1,math.ceil(len(group)*top_fraction))]:
            events[r['chrom']].extend([(int(r['start']),1,anc,r['block_id']),(int(r['end']),-1,anc,r['block_id'])])
    result=[]
    for chrom, ev in sorted(events.items(),key=lambda x:chrom_sort_key(x[0])):
        active={}; previous=None
        import itertools
        for pos, at in itertools.groupby(sorted(ev),key=lambda e:e[0]):
            if previous is not None and pos>previous and len(active)>=2:
                result.append(dict(chrom=chrom,start=previous,end=pos,overlap_bp=pos-previous,
                    supporting_ancestries=','.join(sorted(active)),n_ancestries=len(active),
                    block_ids=','.join(active[a] for a in sorted(active))))
            for _,direction,anc,bid in at:
                if direction<0: active.pop(anc,None)
                else: active[anc]=bid
            previous=pos
    return result
