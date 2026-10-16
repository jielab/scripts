#!/usr/bin/env python3
"""GRID aggregate browser: current twin runs plus existing benchmark XLSX tables.

Reads numeric source values, never fits models or loads participant/model RDS.
Historical residual-correlation R2 is NOT renamed as predictive R2.
"""
from __future__ import annotations
import argparse,json,os,re,sqlite3,tempfile,html,hashlib
from pathlib import Path
import pandas as pd
from twin import sha256,read_table,xlsx_values,public_table,json_write

SKIP={'dashboard','cache','checkpoints','scoring','work','run','runs','neural','__pycache__','evidence','shiny.upstream','shiny_legacy'}
PRIVATE_FILE=re.compile(r'(^|[._-])(individuals?|participants?|predictions?|outcomes?|matches|donors?|roster|split|families|scores\.rds|model)([._-]|$)',re.I)
TRAITS={'height','ldl','t2dm'}


def context(path):
    parts=Path(path).parts;trait=next((p for p in parts if p in TRAITS),'all')
    module=next((p for p in parts if p in ('csx','disco','prsformer','abm','pca','final')),'other')
    if 'twin' in parts:module='Twin'
    return trait,module


def kind_of(d):
    c=set(d)
    if {'contrast','delta','lower','upper','metric'}<=c:return 'contrasts'
    if {'mean_confidence','method'}<=c:return 'confidence'
    if {'actual_coverage','subset','method'}<=c:return 'metrics'
    if {'estimate','method'}<=c:return 'legacy_performance'
    if {'difference','method'}<=c:return 'legacy_contrasts'
    if {'observed_loss','confidence_bin'}<=c:return 'legacy_confidence'
    if {'actual_coverage_among_eligible','observed_loss'}<=c:return 'legacy_coverage'
    if {'reference_eligible_N','retained_N'}<=c:return 'training'
    return 'table'


class Index:
    def __init__(self,root,out,max_mb=128):
        self.root=Path(root).resolve();self.out=Path(out).resolve();self.limit=max_mb*1024**2
        if not self.root.is_dir():raise FileNotFoundError(self.root)
        self.out.mkdir(parents=True,exist_ok=True)
        fd,p=tempfile.mkstemp(prefix='grid-',suffix='.sqlite.tmp',dir=self.out);os.close(fd)
        self.tmp=Path(p);self.db=sqlite3.connect(p);self.catalogue=[];self.audit=[];self.owned=set()
        # Only this disposable database is built without per-table disk syncs.
        # Flush the complete file once before atomically publishing it below.
        self.db.execute('PRAGMA journal_mode=MEMORY')
        self.db.execute('PRAGMA synchronous=OFF')
    def note(self,source,status,detail):self.audit.append(dict(source=str(source),status=status,detail=str(detail)))
    def sources(self):
        # Prune excluded trees before traversal: UKB caches can contain hundreds
        # of thousands of files that the aggregate browser must never inspect.
        for directory,dirs,files in os.walk(self.root,followlinks=False):
            parent=Path(directory)
            dirs[:]=sorted(d for d in dirs if not (parent/d).is_symlink() and self.allowed(parent/d,within_root=True))
            for name in sorted(files):
                if name!='table_manifest.json' and (not name.endswith(('.xlsx','.tsv','.tsv.gz','.csv','.csv.gz')) or PRIVATE_FILE.search(name)):continue
                path=parent/name
                if self.allowed(path,within_root=not path.is_symlink()):yield path
    def allowed(self,p,within_root=False):
        rel=p.relative_to(self.root)
        if p.is_relative_to(self.out) or any(x.startswith(('.', '_')) or x in SKIP for x in rel.parts):return False
        # os.walk reaches only verified nonsymlink parents. Resolve individual
        # links/manifest targets, avoiding repeated ancestor stats on mounted disks.
        if not within_root:
            resolved=p.resolve()
            if not resolved.is_relative_to(self.root) or resolved.is_relative_to(self.out):return False
        if 'twin' in rel.parts:
            pos=rel.parts.index('twin')
            if len(rel.parts)<=pos+2:return p.is_dir()
            parent=self.root.joinpath(*rel.parts[:pos+2]);pointer=parent/'CURRENT.json'
            if not pointer.is_file():return False
            run=json.loads(pointer.read_text()).get('run')
            if rel.parts[pos+2]!=run or not (parent/run/'COMPLETE.json').is_file():return False
        return True
    def add(self,d,source,name,meta=None):
        d=pd.DataFrame(d)
        if not len(d.columns):d=pd.DataFrame(columns=['status'])
        try:public_table(d)
        except ValueError as e:self.note(source,'blocked_private',e);return
        if len(d)>1000000:self.note(source,'size_limit','More than one million aggregate rows');return
        if d.columns.duplicated().any():self.note(source,'invalid_schema','Duplicate columns');return
        for c in d.select_dtypes(include='object'):
            d[c]=d[c].map(lambda v:json.dumps(v,default=str,ensure_ascii=False) if isinstance(v,(list,dict,tuple)) else v)
        trait,module=context(source);md=meta or {};trait=md.get('trait',trait);module=md.get('module',module)
        tid=hashlib.sha256((str(source)+'|'+name).encode()).hexdigest()[:20];sql='t_'+tid
        d.to_sql(sql,self.db,index=False,if_exists='fail')
        state='empty' if not len(d) else 'available'
        if set(d)&{'classification_status','validation_status'}:state='status_or_validation_table'
        original=str(source).split('::')[0];path=self.root/original
        legacy='twin' not in Path(original).parts
        self.catalogue.append(dict(table_id=tid,sql_table=sql,trait=trait,module=module,name=name,rows=len(d),source=str(source),source_sha256=sha256(path),
                                   run_id=md.get('run_id','legacy_saved_result'),plot_kind=kind_of(d),status=state,
                                   legacy=legacy,columns_json=json.dumps(list(map(str,d.columns))),
                                   note='Previously generated source values; no refit and no metric relabeling' if legacy else 'Frozen genetic-profile twin experiment; selected-subset gain is separate from training gain'))
    def run(self,traits=None):
        try:
            paths=list(self.sources())
            for path in (p for p in paths if p.name=='table_manifest.json'):
                try:
                    manifest=json.loads(path.read_text())
                    for item in manifest:
                        name=item['path']
                        if Path(name).name!=name:raise ValueError('Unsafe manifest target')
                        target=path.parent/name
                        self.owned.add(target)
                        if not self.allowed(target):continue
                        if traits and item.get('trait') not in set(traits)|{'all',None}:continue
                        if not target.exists():self.note(target,'missing','Manifest target missing');continue
                        if sha256(target)!=item['sha256']:self.note(target,'checksum_mismatch','Aggregate changed since publication');continue
                        if target.stat().st_size>self.limit:self.note(target,'size_limit','Table file over cap');continue
                        self.add(read_table(target),str(target.relative_to(self.root)),item['name'],item)
                        self.owned.add(target)
                except Exception as e:self.note(path,'manifest_failed',e)
            for path in paths:
                if not path.is_file() or path in self.owned or PRIVATE_FILE.search(path.name):continue
                # New twin publications require their manifest even if it is
                # missing or malformed. Never fall back to unverified raw TSVs.
                if 'twin' in path.relative_to(self.root).parts:continue
                trait,module=context(path.relative_to(self.root))
                if traits and trait not in set(traits)|{'all'}:continue
                rel=str(path.relative_to(self.root))
                try:
                    if path.stat().st_size>self.limit:continue
                    if path.suffix.lower()=='.xlsx':
                        for sheet,frame in xlsx_values(path,self.limit).items():self.add(frame,rel+'::'+sheet,sheet)
                    elif path.name.endswith(('.tsv','.tsv.gz','.csv','.csv.gz')):
                        self.add(read_table(path),rel,path.name)
                except Exception as e:self.note(rel,'source_read_failed',e)
            columns=['table_id','sql_table','trait','module','name','rows','source','source_sha256','run_id','plot_kind','status','legacy','columns_json','note']
            cat=pd.DataFrame(self.catalogue,columns=columns);audit=pd.DataFrame(self.audit,columns=['source','status','detail'])
            cat.to_sql('catalogue',self.db,index=False,if_exists='replace');audit.to_sql('audit',self.db,index=False,if_exists='replace')
            meta=pd.DataFrame([dict(tables=len(cat),audit_items=len(audit),root=str(self.root),read_only=True,individuals_imported=False,models_refit=False,
                completeness='partial' if len(audit) else 'available_aggregate_sources_indexed')])
            meta.to_sql('metadata',self.db,index=False,if_exists='replace');self.db.commit();self.db.close()
            with self.tmp.open('rb') as stream:os.fsync(stream.fileno())
            os.replace(self.tmp,self.out/'evidence.sqlite')
            cat.to_csv(self.out/'catalogue.tsv',sep='\t',index=False);audit.to_csv(self.out/'audit.tsv',sep='\t',index=False)
            page='<!doctype html><meta charset="utf-8"><title>GRID results</title><style>body{font:15px system-ui;max-width:1280px;margin:30px auto;padding:0 20px}table{border-collapse:collapse;width:100%}th,td{padding:7px;border-bottom:1px solid #ddd;text-align:left}details{margin:20px 0}</style><h1>GRID · Results catalogue</h1><p>Select a trait, source table and metric in Shiny. Original benchmark results and current twin runs remain distinguishable. No automatic superiority claim.</p>'
            page+=meta.to_html(index=False,escape=True)+cat[['trait','module','name','rows','status','legacy']].to_html(index=False,escape=True)
            page+='<details><summary>Import audit</summary>'+audit.to_html(index=False,escape=True)+'</details>'
            (self.out/'index.html').write_text(page,encoding='utf-8')
            return dict(tables=len(cat),audit_items=len(audit),db=str(self.out/'evidence.sqlite'))
        except BaseException:
            try:self.db.close()
            finally:self.tmp.unlink(missing_ok=True)
            raise


def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--analysis-root',default='/mnt/d/analysis/grid');a.add_argument('--out');a.add_argument('--traits',default='height,ldl,t2dm');a.add_argument('--max-table-mb',type=int,default=128);a.add_argument('--summary-only',action='store_true',help='Accepted; indexing never refits models')
    x=a.parse_args()
    if x.max_table_mb<=0:a.error('Positive size cap required')
    print(json.dumps(Index(x.analysis_root,x.out or str(Path(x.analysis_root)/'evidence'),x.max_table_mb).run(x.traits.split(','))))
if __name__=='__main__':main()
