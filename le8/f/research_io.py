"""Portable tabular IO and strict score alignment. No implicit phenotype recoding."""
from __future__ import annotations
import csv, gzip, hashlib, io, json, os, re, tempfile, zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np
import pandas as pd

NS = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
PRIVATE = re.compile(r'(^|[._ #])(eid|iid|fid|participant_id|patient_id|donor_id|query_id|subject_id)($|[._ ])', re.I)


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def read_table(path, **kwargs):
    p=Path(path)
    if p.suffix.lower()=='.rds':
        import pyreadr
        out=pyreadr.read_r(str(p))
        if len(out)!=1: raise ValueError('RDS must contain exactly one data.frame')
        return next(iter(out.values()))
    if p.suffix.lower()=='.parquet': return pd.read_parquet(p,**kwargs)
    if p.suffix.lower()=='.xlsx':
        sheet=kwargs.pop('sheet',None)
        frames=read_xlsx(p)
        if sheet is None and len(frames)!=1: raise ValueError('Select an explicit XLSX sheet')
        return frames[sheet] if sheet else next(iter(frames.values()))
    return pd.read_csv(p,sep=',' if p.name.endswith(('.csv','.csv.gz')) else '\t',**kwargs)


def id_table(path, id_col='eid'):
    # Read IDs as text so leading zeros cannot be silently lost in delimited inputs.
    d=read_table(path,dtype={id_col:str}) if not str(path).endswith(('.rds','.parquet','.xlsx')) else read_table(path)
    if id_col not in d: raise ValueError(f'{path}: missing ID {id_col}')
    if d[id_col].isna().any(): raise ValueError('Missing IDs')
    
    if pd.api.types.is_float_dtype(d[id_col]) and np.equal(d[id_col],np.floor(d[id_col])).all():
        d[id_col]=d[id_col].astype('int64')
    d[id_col]=d[id_col].astype(str)
    if d[id_col].duplicated().any() or d[id_col].eq('').any(): raise ValueError('Duplicate/empty IDs')
    return d.set_index(id_col,drop=False)


def win_to_local(path):
    s=str(path)
    if os.name!='nt' and re.match(r'^[A-Za-z]:[\\/]',s):
        s='/mnt/'+s[0].lower()+'/'+s[3:].replace('\\','/')
    return Path(s)


def score_source(root, layer):
    p=win_to_local(root)
    if p.is_file(): return p
    path=p/f'{layer}.pgs.gz'
    if path.exists(): return path
    candidates=list(p.glob('*.pgs.gz'))
    if len(candidates)==1: return candidates[0]
    raise ValueError(f'Cannot unambiguously find {layer}.pgs.gz in {p}; run format.sh pgs step2 or supply an exact file')


def load_pgs(root, layer, ids, measured_features, mapping=None):
    """Reads format.sh's merged eid + <assay>.SCORE_SUM; row-order-independent.

    A per-trait #IID/ALLELE_CT/SCORE_SUM file is also accepted. SCORE_AVG is NOT
    silently substituted. Constant/header-only scores are unavailable, never zeros.
    A wide merged score cannot certify individual variant-pattern QC.
    """
    path=score_source(root,layer)
    opener=gzip.open if path.name.endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8') as h:head=next(csv.reader(h,delimiter='\t'))
    if len(set(head))!=len(head):raise ValueError('Duplicate raw PGS column names')
    id_candidates=[c for c in ('eid','#IID','IID') if c in head]
    if len(id_candidates)!=1: raise ValueError('PGS needs exactly one eid/#IID/IID column')
    id_col=id_candidates[0]
    amap={}
    for c in head:
        if c.endswith('.SCORE_SUM'): amap[c]=c[:-len('.SCORE_SUM')]
    if mapping is not None:
        mp=read_table(mapping)
        if not {'pgs_column','feature'}<=set(mp): raise ValueError('Mapping needs pgs_column,feature')
        if mp.pgs_column.duplicated().any() or mp.feature.duplicated().any(): raise ValueError('Ambiguous mapping')
        amap.update(dict(zip(mp.pgs_column,mp.feature)))
    amap={c:f for c,f in amap.items() if c in head and f in measured_features}
    if len(set(amap.values()))!=len(amap): raise ValueError('Multiple PGSs map to one assay; use an explicit source-specific run')
    if not amap: raise ValueError('No exact PGS/assay matches; do not use fuzzy gene/assay matching')
    features=[f for f in measured_features if f in set(amap.values())]
    cols=[next(c for c,f2 in amap.items() if f2==f) for f in features]
    # Keep only matched columns, and align by ID, not position.
    requested=set(map(str,ids));chunks=[];seen=set()
    for chunk in pd.read_csv(path,sep='\t',usecols=[id_col]+cols,dtype={id_col:str},chunksize=20000):
        current=chunk[id_col]
        if current.isna().any() or current.duplicated().any() or set(current)&seen:raise ValueError('PGS duplicate/missing IDs')
        seen.update(current)
        chunks.append(chunk[current.isin(requested)])
    d=pd.concat(chunks,ignore_index=True) if chunks else pd.DataFrame(columns=[id_col]+cols)
    d=d.set_index(id_col).reindex(list(map(str,ids)))
    g=d[cols].apply(pd.to_numeric,errors='raise').to_numpy(float)
    if np.isinf(g).any(): raise ValueError('Infinite PGS value')
    audit=[]
    for j,f in enumerate(features):
        x=g[:,j];v=x[np.isfinite(x)];u=len(np.unique(v))
        audit.append(dict(feature=f,pgs_column=cols[j],source=str(path),source_sha256=sha256(path) if j==0 else '',
            matched_N=len(v),unique_scores=u,rounded_3dp_fraction=float(np.mean(np.abs(v*1000-np.round(v*1000))<1e-7)) if len(v) else np.nan,
            status='available' if u>1 else 'constant_or_unavailable',variant_pattern_qc='not_inferred_from_scores'))
    return features,g,pd.DataFrame(audit)


def merge_native_scores(manifest, output):
    """Re-merge existing per-trait scores at full precision; never regenerate GWAS.
    Fails on missing people rather than inventing zero scores. Output is create-only.
    """
    manifest=win_to_local(manifest);output=win_to_local(output)
    if output.exists(): raise FileExistsError(output)
    m=read_table(manifest)
    if not {'GWAS','STATUS','PGS'}<=set(m): raise ValueError('Need format.sh <layer>.pgs.files.tsv manifest')
    incomplete=m[~m.STATUS.isin(['included','empty_no_matched_variants'])]
    if len(incomplete): raise ValueError('Manifest contains incomplete traits')
    merged=None;seen=set()
    for r in m.itertuples(index=False):
        if r.STATUS!='included': continue
        if r.GWAS in seen: raise ValueError('Duplicate trait in PGS manifest')
        seen.add(r.GWAS);p=win_to_local(r.PGS)
        if not p.is_absolute():p=manifest.parent/p
        d=pd.read_csv(p,sep='\t',dtype={'#IID':str})
        col=f'{r.GWAS}.SCORE_SUM'
        if not {'#IID',col}<=set(d) or d.empty or d['#IID'].isna().any() or d['#IID'].duplicated().any(): raise ValueError(f'Invalid native score: {p}')
        s=pd.to_numeric(d.set_index('#IID')[col],errors='raise')
        if not np.isfinite(s).all(): raise ValueError('Nonfinite score')
        if merged is None: merged=s.to_frame()
        else:
            if set(merged.index)!=set(s.index): raise ValueError('Different PGS participant sets')
            merged[col]=s.reindex(merged.index)
    if merged is None: raise ValueError('No nonempty PGSs')
    output.parent.mkdir(parents=True,exist_ok=True)
    merged.index.name='eid';merged.to_csv(output,sep='\t',float_format='%.17g',compression='infer')
    return dict(output=str(output),sha256=sha256(output),assays=merged.shape[1],N=len(merged),precision='17 significant digits; cannot recover previously rounded native inputs')


def _letter(n):
    s=''
    while n:
        n,r=divmod(n-1,26);s=chr(65+r)+s
    return s


def write_xlsx(path, tables):
    """Minimal standards-compliant OOXML writer for the user's portable code.
    All strings are inline strings, not formulas. Clinical outputs remain numeric.
    """
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    from xml.sax.saxutils import escape,quoteattr
    ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    rel='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    sheets=[];payload={};used=set()
    for i,(name,df) in enumerate(tables.items(),1):
        name=re.sub(r'[\[\]:*?/\\]','_',str(name))[:31] or f'Sheet{i}'
        if name in used: name=name[:25]+'_'+str(i)
        used.add(name)
        if len(df)>1048575: raise ValueError('Sheet exceeds Excel row limit; use table store instead')
        rows=[]
        for n,record in enumerate([list(map(str,df.columns))]+list(df.itertuples(index=False,name=None)),1):
            cells=[]
            for j,v in enumerate(record,1):
                ref=_letter(j)+str(n);style=' s="1"' if n==1 else ''
                if v is None or (isinstance(v,(float,np.floating)) and not np.isfinite(v)): continue
                if isinstance(v,(bool,np.bool_)): c=f'<c r="{ref}" t="b"{style}><v>{int(v)}</v></c>'
                elif isinstance(v,(int,float,np.number)): c=f'<c r="{ref}"{style}><v>{v}</v></c>'
                else: c=f'<c r="{ref}" t="inlineStr"{style}><is><t xml:space="preserve">{escape(str(v))}</t></is></c>'
                cells.append(c)
            rows.append(f'<row r="{n}">'+''.join(cells)+'</row>')
        cols=f'<cols><col min="1" max="{max(1,len(df.columns))}" width="20" customWidth="1"/></cols>'
        payload[f'xl/worksheets/sheet{i}.xml']=f'<worksheet xmlns="{ns}"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" state="frozen"/></sheetView></sheetViews>{cols}<sheetData>'+''.join(rows)+'</sheetData></worksheet>'
        sheets.append(f'<sheet name={quoteattr(name)} sheetId="{i}" r:id="rId{i}"/>')
    payload['xl/workbook.xml']=f'<workbook xmlns="{ns}" xmlns:r="{rel}"><sheets>'+''.join(sheets)+'</sheets></workbook>'
    rns='http://schemas.openxmlformats.org/package/2006/relationships'
    payload['xl/_rels/workbook.xml.rels']=f'<Relationships xmlns="{rns}">'+''.join(f'<Relationship Id="rId{i}" Type="{rel}/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,len(sheets)+1))+f'<Relationship Id="style" Type="{rel}/styles" Target="styles.xml"/></Relationships>'
    payload['_rels/.rels']=f'<Relationships xmlns="{rns}"><Relationship Id="rId1" Type="{rel}/officeDocument" Target="xl/workbook.xml"/></Relationships>'
    payload['xl/styles.xml']=f'<styleSheet xmlns="{ns}"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts><fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills><borders count="1"><border/></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs></styleSheet>'
    types=['<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>','<Default Extension="xml" ContentType="application/xml"/>','<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>','<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>']
    types += [f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,len(sheets)+1)]
    payload['[Content_Types].xml']='<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'+''.join(types)+'</Types>'
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
        for n,v in payload.items(): z.writestr(n,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'+v)


def read_xlsx(path, max_bytes=256*1024**2):
    """Read ordinary XLSX worksheet VALUES without executing formulas/macros.
    Raw embedded exports are indexed separately by evidence for full precision.
    """
    with zipfile.ZipFile(path) as z:
        if sum(x.file_size for x in z.infolist())>max_bytes: raise ValueError('Workbook exceeds expanded size cap')
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            root=ET.fromstring(z.read('xl/sharedStrings.xml'))
            shared=[''.join(e.itertext()) for e in root.findall('s:si',NS)]
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        rm={e.attrib['Id']:e.attrib['Target'] for e in rels}
        wb=ET.fromstring(z.read('xl/workbook.xml'));out={}
        for sh in wb.findall('s:sheets/s:sheet',NS):
            rid=sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
            target=rm[rid];target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            if '..' in Path(target).parts: raise ValueError('Unsafe worksheet path')
            root=ET.fromstring(z.read(target));rows=[]
            for row in root.findall('s:sheetData/s:row',NS):
                vals=[]
                for c in row.findall('s:c',NS):
                    ref=re.match(r'[A-Z]+',c.attrib['r']).group();j=0
                    for letter in ref:j=j*26+ord(letter)-64
                    while len(vals)<j:vals.append(None)
                    typ=c.attrib.get('t');v=c.find('s:v',NS)
                    if typ=='inlineStr':value=''.join(c.find('s:is',NS).itertext())
                    elif v is None:value=None
                    elif typ=='s':value=shared[int(v.text)]
                    elif typ=='b':value=v.text=='1'
                    elif typ in ('str','e'):value=v.text
                    else:
                        try:value=float(v.text)
                        except (TypeError,ValueError):value=v.text
                    vals[j-1]=value
                rows.append(vals)
            if not rows:out[sh.attrib['name']]=pd.DataFrame();continue
            width=max(map(len,rows));header=rows[0]+[None]*(width-len(rows[0]));seen={};names=[]
            for j,h in enumerate(header):
                h=str(h) if h is not None else f'column_{j+1}';seen[h]=seen.get(h,0)+1
                names.append(h if seen[h]==1 else f'{h}__{seen[h]}')
            out[sh.attrib['name']]=pd.DataFrame([r+[None]*(width-len(r)) for r in rows[1:]],columns=names)
        return out


def write_tables(out, tables, module, layer, trait='cvd_cad',run_id='', workbook='results.xlsx'):
    """One workbook + machine-readable aggregate TSVs; no participant identifiers."""
    out=Path(out);out.mkdir(parents=True,exist_ok=True);registry=[];safe={}
    for name,d in tables.items():
        d=pd.DataFrame(d)
        if not len(d.columns):d=pd.DataFrame(columns=['status'])
        if any(PRIVATE.search(c) for c in d.columns):raise ValueError(f'Private column in public table {name}')
        stem=re.sub('[^a-zA-Z0-9_.-]','_',name);p=out/(stem+'.tsv.gz')
        d.to_csv(p,sep='\t',index=False)
        registry.append(dict(name=name,path=p.name,module=module,layer=layer,trait=trait,run_id=run_id,
                             rows=len(d),sha256=sha256(p),status='available' if len(d) else 'empty'))
        if len(d)<=200000 and len(d)*len(d.columns)<=2000000:safe[name]=d
    if safe:write_xlsx(out/workbook,safe)
    (out/'table_manifest.json').write_text(json.dumps(registry,ensure_ascii=False,indent=2))
    return registry
