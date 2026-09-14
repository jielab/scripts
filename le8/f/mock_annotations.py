"""Public annotation enrichment. Input/output contain gene symbols only."""
import csv,hashlib,io,json,math,sys,time,os,tempfile
from pathlib import Path
import requests
src,out,cache=map(Path,sys.argv[1:]);cache.mkdir(parents=True,exist_ok=True)
r=list(csv.DictReader(src.open()));background={x['gene'].upper() for x in r};sig={x['gene'].upper() for x in r if x['selected'].upper()=='TRUE'}
terms=[];tf=[];ppi=[];audit=[]
def atomic_text(path,text):
 fd,tmp=tempfile.mkstemp(dir=path.parent)
 try:
  with os.fdopen(fd,'w') as f:f.write(text)
  os.replace(tmp,path)
 finally:
  if os.path.exists(tmp):os.unlink(tmp)
def write(name,rows,cols):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
def tail(k,m,N,n):
 def lc(a,b):return math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)
 logs=[lc(m,j)+lc(N-m,n-j)-lc(N,n) for j in range(max(k,n-(N-m)),min(m,n)+1)]
 if not logs:return 1.
 z=max(logs);return min(1.,math.exp(z)*sum(math.exp(x-z) for x in logs))
def bh(rows):
 rows.sort(key=lambda x:x['p_raw']);best=1.
 for i in range(len(rows)-1,-1,-1):best=min(best,rows[i]['p_raw']*len(rows)/(i+1));rows[i]['adjusted_p']=best
 return rows
def library(label,names):
 for name in names:
  p=cache/(name+'.gmt')
  try:
   if not p.exists():
    z=requests.get('https://maayanlab.cloud/Enrichr/geneSetLibrary',params={'mode':'text','libraryName':name},timeout=45);z.raise_for_status()
    if '\t' not in z.text:continue
    atomic_text(p,z.text)
   lines=p.read_text().splitlines();rows=[]
   for line in lines:
    c=line.split('\t')
    if label=='TF' and not c[0].lower().endswith(' human'):continue
    g={x.split(',')[0].upper() for x in c[2:]}&background
    if len(g)<3:continue
    hits=g&sig;k=len(hits)
    # Retain every eligible term for multiplicity, including zero overlaps.
    rows.append(dict(source=label,term_name=c[0],p_raw=tail(k,len(g),len(background),len(sig)) if k else 1.,intersection_size=k,genes=';'.join(sorted(hits))))
   if not rows:continue
   audit.append(dict(source=label,status='ok',detail=name));return bh(rows)
  except Exception as e:audit.append(dict(source=label,status='unavailable',detail=str(e)))
 return []
def mgi_high_level():
 try:
  paths=[]
  for name in ['HMD_HumanPhenotype.rpt','VOC_MammalianPhenotype.rpt']:
   p=cache/name
   if not p.exists():
    z=requests.get('https://www.informatics.jax.org/downloads/reports/'+name,timeout=45);z.raise_for_status();atomic_text(p,z.text)
   paths.append(p)
  names={c[0]:c[1] for c in csv.reader(paths[1].open(),delimiter='\t') if len(c)>1}
  groups={}
  for c in csv.reader(paths[0].open(),delimiter='\t'):
   if len(c)<5 or c[0].upper() not in background:continue
   for term in c[4].replace(',',' ').split():groups.setdefault(term,set()).add(c[0].upper())
  universe=set().union(*groups.values()) if groups else set();query=sig&universe;rows=[]
  if not universe or not query:return []
  for term,g in groups.items():
   hits=g&query
   rows.append(dict(source='MGI',term_name=names.get(term,term),p_raw=tail(len(hits),len(g),len(universe),len(query)) if hits else 1.,intersection_size=len(hits),genes=';'.join(sorted(hits))))
  audit.append(dict(source='MGI',status='ok',detail=f'MGI high-level phenotypes; annotated assayed background N={len(universe)}, selected N={len(query)}'))
  return bh(rows)
 except Exception as e:
  audit.append(dict(source='MGI',status='unavailable',detail=str(e)));return []
if sig:
 terms+=mgi_high_level()
 for label,names in [('KEGG',['KEGG_2021_Human']),('TF',['TRRUST_Transcription_Factors_2019'])]:
  rows=library(label,names);terms+=rows
  if label=='TF':
   for r in [x for x in rows if x['adjusted_p']<.05][:15]:
    regulator=r['term_name'].split('_')[0].split(' ')[0]
    for g in r['genes'].split(';'):
     if g and regulator!=g:tf.append({'from':regulator,'to':g})
 key=hashlib.sha256('\n'.join(sorted(sig)).encode()).hexdigest()[:20];p=cache/('string_physical_'+key+'.tsv')
 try:
  if not p.exists():
   z=requests.post('https://version-12-0.string-db.org/api/tsv/network',data={'identifiers':'\r'.join(sorted(sig)),'species':9606,'network_type':'physical','required_score':400},timeout=60);z.raise_for_status();atomic_text(p,z.text)
  for r in csv.DictReader(io.StringIO(p.read_text()),delimiter='\t'):
   a,b=r.get('preferredName_A',''),r.get('preferredName_B','')
   if a.upper() in sig and b.upper() in sig:ppi.append({'from':a,'to':b,'score':r['score']})
  audit.append(dict(source='STRING',status='ok',detail=f'{len(ppi)} physical edges'))
 except Exception as e:audit.append(dict(source='STRING',status='unavailable',detail=str(e)))
write('c1.mock_function_terms.csv',terms,['source','term_name','p_raw','intersection_size','genes','adjusted_p'])
write('c1.mock_tf_edges.csv',tf,['from','to']);write('c1.mock_ppi_edges.csv',ppi,['from','to','score'])
write('c1.mock_annotation_status.csv',audit,['source','status','detail'])
print('MOCK enrichment:',len(terms),'terms;',len(tf),'TF edges;',len(ppi),'physical edges',flush=True)
