#!/usr/bin/env python3
"""Read-only result indexing for modular table/figure-first Shiny.

Reads exact exports in existing le8-workbook-v1 XLSX archives and new manifests.
Never loads pickle/RDS/model/participant objects. Missing evidence is not zero.
"""
from __future__ import annotations
import argparse,hashlib,html,io,json,os,re,sqlite3,tempfile,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from research_io import sha256,read_xlsx,read_table,write_xlsx,PRIVATE

FEATURES=('feature','assay','protein','Protein','metabolite','Metabolite','exposure','term','gene')
PRIVATE_PATH=re.compile(r'(^|[/_.-])(private|LOCAL_ONLY|input|inputs|neural|checkpoint|outcomes|predictions|participant|individual|matches|roster|split)([/_.-]|$)',re.I)
PRIVATE_EXACT={'eid','iid','fid','id','id1','id2','participantid','individualid','subjectid','patientid','queryid','targetid','referenceid','donorid','ukbid','ukbeid','sampleid','birthdate','datebirth'}


def has_private(d):
    names={re.sub('[^a-z0-9]','',str(c).lower()) for c in d.columns}
    allowed={'referenceid'} if {'runid','modelid','primaryid','referenceid','n','unochorizon'}<=names else set()
    return bool((names-allowed)&PRIVATE_EXACT) or any(PRIVATE.search(str(c)) for c in d.columns)


def private_path(path):
    return any(p.startswith('.') or p.startswith('_') for p in Path(path).parts) or bool(PRIVATE_PATH.search(str(path)))


def context(path):
    parts=Path(path).parts;layer=next((p for p in parts if p in ('prot','met')),'all')
    ix=parts.index(layer) if layer in parts else -1
    trait=parts[ix-1] if ix>0 else 'all'
    module=next((p for p in parts if re.match(r'c[1-5]_',p)),'other')
    if 'abm_twin' in parts or 'twin' in str(path).lower():module='ABM twin'
    elif 'deviation' in str(path).lower():module='Molecular deviation'
    elif module=='other' and 'final' in parts:module='Final'
    return dict(module=module,layer=layer,trait=trait)


def evidence_status(d):
    if not len(d):return 'empty','No numerical result rows; not a null hypothesis result'
    if set(map(str.lower,d.columns))<= {'status','note','message','detail','reason','caption'}:
        values='; '.join(d.fillna('').astype(str).to_numpy().ravel())[:600]
        return 'status_only',values
    return 'available','Source values; no causal or clinical claim inferred by the viewer'


def plot_kind(d):
    c=set(d)
    if {'confidence_lo','mean_confidence','Brier_IPCW'}<=c:return 'confidence'
    if 'beta' in c and ('time' in c or 'landmark_years' in c):return 'temporal'
    if {'delta_R2_vs_NS','component','model'}<=c:return 'proxy_reconstruction'
    if c & {'PP_H4','PP.H4.abf','PP.H4'}:return 'coloc'
    if {'FPR','TPR'}<=c:return 'roc'
    if {'net_benefit','threshold'}<=c:return 'decision'
    if {'mean_deviation','state'}<=c:return 'deviation_state'
    if {'bin','lower','upper','state','N'}<=c:return 'deviation_histogram'
    if {'own_PGS_partial_R2','feature'}<=c:return 'pgs_capture'
    if {'coverage','AUC_IPCW','model'}<=c:return 'coverage'
    if {'delta','lower','upper','metric'}<=c:return 'contrasts'
    if {'predicted','observed'}<=c:return 'calibration'
    if {'beta','p.value'}<=c or {'beta','p'}<=c or {'beta','se'}<=c or {'beta','std.error'}<=c:return 'effect'
    if {'pgs_full_beta','observed_beta'}<=c:return 'concordance'
    if {'AUC','model'}<=c:return 'models'
    return 'generic'


class Indexer:
    def __init__(self,root,output,max_mb=128):
        self.root=Path(root).resolve();self.out=Path(output).resolve();self.max_bytes=int(max_mb*1024**2)
        if not self.root.is_dir():raise FileNotFoundError(self.root)
        self.out.mkdir(parents=True,exist_ok=True)
        fd,p=tempfile.mkstemp(prefix='evidence-',suffix='.sqlite.tmp',dir=self.out);os.close(fd)
        self.tmp=Path(p);self.db=sqlite3.connect(p);self.rows=[];self.audit=[];self.seen={};self.owned=set()
    def warn(self,source,status,detail):self.audit.append(dict(source=source,status=status,detail=str(detail)[:800]))
    def add(self,d,name,source,source_hash,meta=None,kind='native'):
        d=pd.DataFrame(d)
        if not len(d.columns):d=pd.DataFrame(columns=['status'])
        if has_private(d):self.warn(source,'blocked_private','Participant-like columns; no rows imported');return
        if len(d)>1000000:self.warn(source,'skipped_size','Table exceeds one million rows; prepare a filtered aggregate view');return
        d.columns=list(map(str,d.columns))
        if len(set(d.columns))!=len(d.columns):self.warn(source,'skipped_schema','Duplicate columns');return
        for c in d.select_dtypes(include='object').columns:
            if any(isinstance(v,(dict,list,tuple,set)) for v in d[c].head(500)):
                d[c]=d[c].map(lambda v:json.dumps(v,default=str) if isinstance(v,(dict,list,tuple,set)) else v)
        md={**context(source),**(meta or {})}
        md['module']={'deviation':'Molecular deviation','c1_correlate':'C1 association and prediction','c2_cause':'C2 genetic evidence','c3_coloc':'C3 colocalization','c4_connect':'C4 connection','c5_cellulation':'C5 cellular context'}.get(md['module'],md['module'])
        # A content hash is exact to this exported table. Never min-P/max-H4 merge.
        content=hashlib.sha256(d.to_csv(index=False).encode()).hexdigest()
        # Different runs retain independent provenance even for identical tables.
        key=(md['trait'],md['layer'],md['module'],md.get('run_id') or str(Path(source.split('::')[0]).parent),name,content)
        if key in self.seen:
            self.warn(source,'duplicate_identical',self.seen[key]);return
        self.seen[key]=source
        tid=hashlib.sha256((source+'|'+name+'|'+content).encode()).hexdigest()[:20]
        table='data_'+tid;d.to_sql(table,self.db,index=False,if_exists='fail')
        feature=next((c for c in FEATURES if c in d),'')
        if feature:
            self.db.execute(f'CREATE INDEX idx_{tid}_feature ON "{table}" ("{feature}")')
            features=pd.DataFrame({'table_id':tid,'feature':d[feature].dropna().astype(str).unique()})
            features.to_sql('feature_index',self.db,index=False,if_exists='append')
        state,note=evidence_status(d)
        self.rows.append(dict(table_id=tid,sql_table=table,name=name,**{k:v for k,v in md.items() if k!='run_id'},run_id=md.get('run_id') or source_hash[:12],
                     source=source,source_sha256=source_hash,content_sha256=content,rows=len(d),columns_json=json.dumps(list(d.columns)),
                     feature_column=feature,status=state,note=note,source_kind=kind,plot_kind=plot_kind(d)))
    def workbook(self,path):
        rel=str(path.relative_to(self.root));h=sha256(path)
        try:
            with zipfile.ZipFile(path) as z:
                if 'le8/manifest.json' in z.namelist():
                    mf=json.loads(z.read('le8/manifest.json'))
                    if mf.get('format')!='le8-workbook-v1':raise ValueError('Unknown archive format')
                    for name,e in mf.get('files',{}).items():
                        source=rel+'::'+name
                        if Path(name).name!=name or private_path(name):self.warn(source,'blocked_private_path','not imported');continue
                        part=e.get('part','')
                        if not re.fullmatch(r'le8/exports/[a-f0-9]{64}\.bin',part):raise ValueError('Unsafe export part')
                        if e.get('kind') not in ('table','workbook'):
                            self.warn(source,'skipped_non_table',e.get('kind','unknown'));continue
                        info=z.getinfo(part)
                        if info.file_size>self.max_bytes:self.warn(source,'skipped_size','Expanded export exceeds cap');continue
                        raw=z.read(part)
                        if len(raw)!=int(e['size']) or hashlib.sha256(raw).hexdigest()!=e['sha256']:raise ValueError('Stored export checksum mismatch')
                        if e['kind']=='workbook':
                            for sh,d in read_xlsx(io.BytesIO(raw),max_bytes=self.max_bytes).items():self.add(d,name+' / '+sh,source,e['sha256'],kind='embedded_workbook')
                        else:
                            import gzip
                            if name.endswith('.gz'):raw=gzip.decompress(raw)
                            if len(raw)>self.max_bytes:self.warn(source,'skipped_size','Decompressed table exceeds cap');continue
                            d=pd.read_csv(io.BytesIO(raw),sep='\t' if '.tsv' in name else ',') if raw.strip() else pd.DataFrame()
                            self.add(d,name,source,e['sha256'],kind='exact_embedded_export')
                    return
            for sh,d in read_xlsx(path,max_bytes=self.max_bytes).items():self.add(d,sh,rel+'::'+sh,h,kind='worksheet_values')
        except Exception as exc:self.warn(rel,'workbook_read_failed',exc)
    def native(self,path):
        meta=json.loads(path.read_text());directory=path.parent
        for e in meta:
            rel=e['path']
            if Path(rel).name!=rel:raise ValueError('Native table manifest path must be a basename')
            source=directory/rel
            if not source.resolve().is_relative_to(self.root):raise ValueError('Native source escapes result root')
            if private_path(str(source.relative_to(self.root))):self.warn(rel,'blocked_private_path','not imported');continue
            if not source.exists():self.warn(rel,'missing','Manifest target missing');continue
            if source.stat().st_size>self.max_bytes:self.warn(rel,'skipped_size','File exceeds cap');continue
            if sha256(source)!=e['sha256']:self.warn(rel,'checksum_mismatch','Index not updated for changed source');continue
            d=read_table(source);md={k:e[k] for k in ('module','layer','trait','run_id') if k in e}
            self.add(d,e['name'],str(source.relative_to(self.root)),e['sha256'],md)
            self.owned.add(source.resolve())
        self.owned.update(p.resolve() for p in directory.glob('*.xlsx'))
    def run(self,traits=None,layers=None):
        def allowed(p):
            rel=p.relative_to(self.root)
            if any(p.resolve().is_relative_to(out) for out in [self.out,self.root/'evidence']) or private_path(str(rel)):return False
            cx=context(rel)
            return (not traits or cx['trait']=='all' or cx['trait'] in traits) and (not layers or cx['layer']=='all' or cx['layer'] in layers)
        try:
            for p in sorted(self.root.rglob('table_manifest.json')):
                if allowed(p):
                    try:self.native(p)
                    except Exception as exc:self.warn(str(p.relative_to(self.root)),'manifest_failed',exc)
            for p in sorted(self.root.rglob('*')):
                if not p.is_file() or not allowed(p) or p.resolve() in self.owned:continue
                if p.suffix.lower()=='.xlsx':self.workbook(p)
                elif p.name.endswith(('.csv','.tsv','.csv.gz','.tsv.gz')):
                    rel=str(p.relative_to(self.root))
                    if p.stat().st_size>self.max_bytes:self.warn(rel,'skipped_size','File exceeds cap');continue
                    try:self.add(read_table(p),p.name,rel,sha256(p),kind='delimited')
                    except Exception as exc:self.warn(rel,'read_failed',exc)
            catalogue=pd.DataFrame(self.rows,columns=['table_id','sql_table','name','module','layer','trait','run_id','source','source_sha256','content_sha256','rows','columns_json','feature_column','status','note','source_kind','plot_kind'])
            catalogue.to_sql('catalogue',self.db,index=False,if_exists='replace')
            audit=pd.DataFrame(self.audit,columns=['source','status','detail']);audit.to_sql('audit',self.db,index=False,if_exists='replace')
            if not self.db.execute("SELECT count(*) FROM sqlite_master WHERE name='feature_index'").fetchone()[0]:
                self.db.execute('CREATE TABLE feature_index(table_id TEXT,feature TEXT)')
            self.db.execute('CREATE INDEX IF NOT EXISTS feature_lookup ON feature_index(table_id,feature)')
            # Scientific states are reported separately; no automatic claim from counts.
            meta=pd.DataFrame([dict(format='le8-evidence',source_root=str(self.root),tables=len(catalogue),audit_rows=len(audit),
                 completeness='partial' if audit.status.isin(['skipped_size','workbook_read_failed','read_failed','manifest_failed','checksum_mismatch']).any() else 'indexed_available_aggregate_sources',
                 models_refit=False,private_records_imported=False)])
            meta.to_sql('metadata',self.db,index=False,if_exists='replace');self.db.commit();self.db.close()
            os.replace(self.tmp,self.out/'evidence.sqlite')
            write_xlsx(self.out/'catalogue.xlsx',{'tables':catalogue,'audit':audit,'metadata':meta})
            document='<!doctype html><meta charset="utf-8"><title>LE8 results</title><style>body{font:15px system-ui;margin:2rem}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}</style><h1>LE8 · Table and figure index</h1>'
            document+='<p>Read-only aggregate index. Choose one table and one molecule in Shiny. No new model fitting.</p>'
            document+=meta.to_html(index=False,escape=True)+catalogue[['trait','layer','module','name','rows','status','plot_kind']].to_html(index=False,escape=True)
            document+='<details><summary>Import and provenance audit</summary>'+audit.to_html(index=False,escape=True)+'</details>'
            (self.out/'index.html').write_text(document)
            return {'tables':len(catalogue),'audit_rows':len(audit),'database':str(self.out/'evidence.sqlite')}
        except BaseException:
            self.db.close();self.tmp.unlink(missing_ok=True);raise


def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--analysis-root',required=True);a.add_argument('--out');a.add_argument('--Y');a.add_argument('--biom');a.add_argument('--max-table-mb',type=int,default=128)
    x=a.parse_args();out=x.out or str(Path(x.analysis_root)/'evidence')
    print(json.dumps(Indexer(x.analysis_root,out,x.max_table_mb).run(x.Y.split(',') if x.Y else None,x.biom.split(',') if x.biom else None)))
if __name__=='__main__':main()
