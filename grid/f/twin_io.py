"""Lossless local numeric IO and aggregate-only publication utilities."""
from __future__ import annotations
import hashlib,json,os,re,uuid,gzip,io,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pandas as pd

PRIVATE_NAMES={'eid','iid','fid','participantid','individualid','subjectid','patientid','queryid','targetid','donorid','familyid','group','id','id1','id2','sampleid','ukbid','ukbeid'}


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def json_write(path,value):
    def safe(v):
        if isinstance(v,pd.DataFrame):return {'type':'dataframe_metadata','rows':len(v),'columns':list(map(str,v.columns)),'contents_omitted':True}
        if isinstance(v,pd.Series):return {'type':'series_metadata','rows':len(v),'name':str(v.name),'contents_omitted':True}
        if isinstance(v,pd.Timestamp):return v.isoformat()
        if v is pd.NA or v is pd.NaT:return None
        if isinstance(v,dict):return {str(k):safe(x) for k,x in v.items()}
        if isinstance(v,(list,tuple)):return [safe(x) for x in v]
        if isinstance(v,np.ndarray):return safe(v.tolist())
        if isinstance(v,np.generic):return safe(v.item())
        if isinstance(v,Path):return str(v)
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp.'+uuid.uuid4().hex)
    tmp.write_text(json.dumps(safe(value),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');os.replace(tmp,p)


def read_table(path):
    p=Path(path)
    if p.suffix=='.rds':
        try:import rdata
        except ImportError as e:raise RuntimeError('RDS input needs rdata>=1.1; use the existing grid environment or supply a TSV') from e
        d=rdata.read_rds(p)
        if isinstance(d,dict):
            candidates=[x for x in d.values() if isinstance(x,pd.DataFrame)]
            if len(candidates)!=1:raise ValueError('Expected exactly one data.frame in RDS')
            d=candidates[0]
        if not isinstance(d,pd.DataFrame):raise ValueError('Expected a data.frame in RDS')
        return d
    if p.suffix=='.parquet':return pd.read_parquet(p)
    sep=',' if '.csv' in p.name else '\t'
    # Explicit string IDs preserve values such as 00012; pandas otherwise infers integers.
    head=pd.read_csv(p,sep=sep,header=None,nrows=1).iloc[0].astype(str)
    if head.duplicated().any():raise ValueError('Duplicate source columns')
    dtype={c:str for c in head if re.sub(r'[^a-z0-9]','',c.lower()) in PRIVATE_NAMES or c=='ancestry'}
    return pd.read_csv(p,sep=sep,dtype=dtype,low_memory=False)


def identified(d):
    d=pd.DataFrame(d).copy()
    if 'eid' not in d:
        c=next((c for c in ('IID','#IID','ID') if c in d),None)
        if c is None:raise ValueError('Participant table needs eid or IID')
        d=d.rename(columns={c:'eid'})
    if d.eid.isna().any():raise ValueError('Missing IDs')
    if pd.api.types.is_numeric_dtype(d.eid):
        v=d.eid.to_numpy(float)
        if not np.isfinite(v).all() or np.any(v!=np.floor(v)):raise ValueError('Noninteger numeric IDs')
        d['eid']=d.eid.astype('int64').astype(str)
    else:d['eid']=d.eid.astype(str).str.strip()
    if d.eid.duplicated().any() or d.eid.eq('').any():raise ValueError('Duplicate/blank IDs')
    return d


def public_table(d):
    for c in d:
        v=re.sub(r'[^a-z0-9]','',str(c).lower())
        if v in PRIVATE_NAMES or re.search(r'(^|[_. #])(eid|iid|fid)($|[_. ])',str(c),re.I):
            raise ValueError('Participant identifiers must remain in _LOCAL_ONLY: '+str(c))
    return d


def write_aggregate(directory,tables,metadata):
    out=Path(directory);out.mkdir(parents=True,exist_ok=True);manifest=[]
    for name,value in tables.items():
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',name):raise ValueError('Unsafe table name')
        d=public_table(pd.DataFrame(value).copy())
        if not len(d.columns):d=pd.DataFrame(columns=['status'])
        for c in d.select_dtypes(include='object'):
            d[c]=d[c].map(lambda v:json.dumps(v,ensure_ascii=False,default=str) if isinstance(v,(dict,list,tuple)) else v)
        path=out/(name+'.tsv.gz');tmp=out/(name+'.tmp.gz')
        d.to_csv(tmp,sep='\t',index=False,float_format='%.17g',compression='gzip');os.replace(tmp,path)
        manifest.append(dict(name=name,path=path.name,sha256=sha256(path),rows=len(d),columns=list(map(str,d.columns)),**metadata))
    json_write(out/'table_manifest.json',manifest)
    return manifest


def xlsx_values(path,max_bytes=256*1024**2):
    """Read cached Excel *values* only. Never evaluate formula text or execute macros.
    Supports the current GRID ordinary XLSX aggregates. Empty/status sheets retained.
    This parser is for importing legacy evidence, not editing spreadsheets.
    """
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    def text(e):return ''.join(t.text or '' for t in e.findall('.//s:t',ns))
    with zipfile.ZipFile(path) as z:
        if sum(i.file_size for i in z.infolist())>max_bytes:raise ValueError('Expanded workbook exceeds configured limit')
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():strings=[text(si) for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
        relationships={e.attrib['Id']:e.attrib['Target'] for e in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        w=ET.fromstring(z.read('xl/workbook.xml'));out={}
        for sheet in w.findall('s:sheets/s:sheet',ns):
            rid=sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'];target=relationships[rid]
            target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            if '..' in Path(target).parts:raise ValueError('Unsafe XLSX relationship')
            rows=[]
            for row in ET.fromstring(z.read(target)).findall('s:sheetData/s:row',ns):
                vals=[]
                for cell in row.findall('s:c',ns):
                    letters=re.sub('[0-9]','',cell.get('r','A1'));col=0
                    for char in letters:col=26*col+ord(char)-64
                    while len(vals)<col:vals.append(None)
                    v=cell.find('s:v',ns);v=v.text if v is not None else None;typ=cell.get('t')
                    if typ=='s' and v is not None:v=strings[int(v)]
                    elif typ=='inlineStr':v=text(cell)
                    elif typ=='b' and v is not None:v=v=='1'
                    elif v is not None and typ not in ('str','e'):
                        try:v=float(v)
                        except ValueError:pass
                    vals[col-1]=v
                rows.append(vals)
            if not rows:out[sheet.attrib['name']]=pd.DataFrame();continue
            width=max(map(len,rows));header=rows[0]+[None]*(width-len(rows[0]));names=[str(v) if v is not None else f'column_{j+1}' for j,v in enumerate(header)]
            if len(set(names))!=len(names):raise ValueError('Duplicate worksheet columns')
            out[sheet.attrib['name']]=pd.DataFrame([r+[None]*(width-len(r)) for r in rows[1:]],columns=names)
        return out
