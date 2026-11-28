"""Small shared utilities for the GRCh38 scan and enrichment workflows."""
from __future__ import annotations
import bisect, csv, gzip, hashlib, html, json, math, os, tempfile
from collections import defaultdict
from pathlib import Path


def open_text(path, mode='rt'):
    return gzip.open(path, mode, encoding='utf-8', newline='') if str(path).endswith(('.gz','.bgz')) else open(path, mode.replace('t',''), encoding='utf-8', newline='')


def chrom_name(value):
    value=str(value).removeprefix('chr')
    return 'chr'+('M' if value in ('M','MT') else value)


def chrom_key(value):
    c=chrom_name(value)[3:]
    return (int(c), '') if c.isdigit() else ({'X':23,'Y':24,'M':25}.get(c,99),c)


def chromosomes(value):
    if value=='autosomes': return ['chr'+str(i) for i in range(1,23)]
    if value=='all': return ['chr'+str(i) for i in range(1,23)]+['chrX','chrY','chrM']
    answer=[]
    for part in value.split(','):
        part=part.strip().removeprefix('chr')
        if '-' in part:
            lo,hi=map(int,part.split('-')); answer.extend('chr'+str(i) for i in range(lo,hi+1))
        else: answer.append(chrom_name(part))
    if not answer or len(set(answer))!=len(answer): raise ValueError('Select unique chromosomes')
    return sorted(answer,key=chrom_key)


def number(value):
    if value in (None,'','.','NA','nan','None'): return None
    try: n=float(value)
    except (TypeError,ValueError): return None
    return n if math.isfinite(n) else None


def fmt(value):
    if value is None: return 'NA'
    if isinstance(value,float): return format(value,'.9g') if math.isfinite(value) else 'NA'
    return str(value)


def atomic_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,temp=tempfile.mkstemp(prefix=path.name+'.',suffix='.partial',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f: json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
        os.replace(temp,path)
    finally:
        if os.path.exists(temp):os.unlink(temp)


def identity(path):
    p=Path(path).resolve();s=p.stat()
    return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)


def write_tsv(path,fields,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with open_text(path,'wt') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields)
        for row in rows:
            w.writerow([fmt(row.get(k)) for k in fields] if isinstance(row,dict) else [fmt(x) for x in row])


def bh(values):
    finite=sorted((p,i) for i,p in enumerate(values) if p is not None and math.isfinite(p))
    result=[None]*len(values);m=len(finite);best=1.
    for rank in range(m,0,-1):
        p,i=finite[rank-1];best=min(best,p*m/rank);result[i]=best
    return result


def stable_hash(*parts):
    return int.from_bytes(hashlib.blake2b('|'.join(map(str,parts)).encode(),digest_size=8).digest(),'big')


class Intervals:
    """Merged BED0 intervals indexed independently by chromosome and label."""
    def __init__(self,rows=()):
        raw=defaultdict(list)
        for chrom,start,end,label in rows:
            start,end=int(start),int(end)
            if start<0 or end<=start:raise ValueError('Invalid BED0 interval')
            raw[(chrom_name(chrom),str(label))].append((start,end))
        self.data={};self.by_chrom=defaultdict(list)
        for (c,label),items in raw.items():
            merged=[]
            for a,b in sorted(items):
                if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
                else:merged.append((a,b))
            self.data[(c,label)]=([a for a,b in merged],[b for a,b in merged])
            self.by_chrom[c].append(label)
    def labels(self,chrom,pos0):
        c=chrom_name(chrom);answer=[]
        for label in self.by_chrom.get(c,[]):
            starts,ends=self.data[(c,label)];i=bisect.bisect_right(starts,pos0)-1
            if i>=0 and pos0<ends[i]:answer.append(label)
        return answer
    def rows(self):
        for (c,label),(starts,ends) in sorted(self.data.items(),key=lambda x:(chrom_key(x[0][0]),x[0][1])):
            for a,b in zip(starts,ends):yield c,a,b,label


def bed_rows(path,label=None,family=False):
    header=None
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('track','browser')):continue
            if line.startswith('#'):
                parts=line.lstrip('#').rstrip('\n').split('\t')
                if parts[0].lower() in ('chrom','chromosome'):header=parts
                continue
            a=line.rstrip('\n').split('\t')
            if a[0].lower() in ('chrom','chromosome'):header=a;continue
            if len(a)<3:raise ValueError('BED needs chrom,start0,end0: '+str(path))
            name=label or (a[3] if len(a)>3 else Path(path).stem)
            if family:
                if header and 'family' in header:name=a[header.index('family')]
                elif len(a)>=10 and a[3] in ('erv_core','erv_related','ltr_other','non_ltr'):name=a[4]
                else:name=name.split('|')[0]
            yield a[0],int(a[1]),int(a[2]),name


def html_report(path,title,sections):
    content=['<!doctype html><meta charset="utf-8"><title>'+html.escape(title)+'</title>',
      '<style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:0 24px;color:#173244}h1,h2{color:#123f53}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:8px;border-bottom:1px solid #ddd;text-align:left}th{background:#edf5f5}img{max-width:100%}.note{background:#eef5f6;padding:16px}code{background:#eee;padding:2px}</style>', '<h1>'+html.escape(title)+'</h1>']
    for heading,body in sections:content+=['<h2>'+html.escape(heading)+'</h2>',body]
    Path(path).parent.mkdir(parents=True,exist_ok=True);Path(path).write_text('\n'.join(content))


def html_table(rows,fields,limit=100):
    return '<table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in fields)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(fmt(r.get(k)))+'</td>' for k in fields)+'</tr>' for r in rows[:limit])+'</table>'
