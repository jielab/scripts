#!/usr/bin/env python3
"""GRID FINAL aggregate browser; no model fitting or participant-level imports."""
from __future__ import annotations
import argparse,json,os,re,sqlite3,tempfile,html,hashlib
from pathlib import Path
import pandas as pd
from twin import sha256,read_table,xlsx_values,json_write
from final_outputs import PAPER_FILES, aggregate

SKIP={'abm_GPU','abm_v2','prsformer_v2','calibrated','dashboard','cache','checkpoints','scoring','work','run','runs','neural','__pycache__','evidence','shiny.upstream','shiny_legacy','local'}
PRIVATE_FILE=re.compile(r'(^|[._-])(individuals?|participants?|predictions?|outcomes?|matches|donors?|roster|split|families|scores\.rds|model)([._-]|$)',re.I)
TRAITS={'height','ldl','t2dm'}


def current_abm_position(parts):
    """Locate FINAL ABM publications."""
    for i, part in enumerate(parts):
        if part == 'abm' and len(parts) > i + 1 and parts[i + 1] == 'FINAL':
            return i
    return None


def context(path):
    parts=Path(path).parts;trait=next((p for p in parts if p in TRAITS),'all')
    module=next((p for p in parts if p in ('csx','disco','prsformer','prsformer2','abm','pca','architecture','final')),'other')
    return trait,module


def kind_of(d):
    c=set(d)
    if {'objective','contrast','mse_reduction','comparator_target_mse'}<=c:return 'abm_count_contrasts'
    if {'policy_gain_vs_base.estimate','method','assessment_loss'}<=c:return 'abm_policy'
    if 'matched.gain_enrichment.estimate' in c:
        return 'abm_benefit' if 'selector' in d and d.selector.astype(str).str.contains('benefit').any() else 'abm_predictability'
    if {'purpose','selector','actual_coverage','candidate_gain_same_subset','selector_delta_loss_vs_control'}<=c:return 'abm_selector_curve'
    if {'candidate_minus_fallback_loss','selector_minus_control_loss'}<=c:return 'abm_primary_contrasts'
    if {'contrast','delta_loss','lower','upper','metric'}<=c:return 'abm_reference_cleaning'
    if {'contrast','delta','lower','upper','metric'}<=c:return 'contrasts'
    if {'mean_confidence','method'}<=c:return 'confidence'
    if {'actual_coverage','subset','method'}<=c:return 'metrics'
    if {'estimate','method'}<=c:return 'performance'
    if {'difference','method'}<=c:return 'contrasts'
    if {'observed_loss','confidence_bin'}<=c:return 'confidence'
    if {'actual_coverage_among_eligible','observed_loss'}<=c:return 'coverage'
    if {'reference_eligible_N','retained_N'}<=c:return 'training'
    return 'table'


class Index:
    def __init__(self,root,out,max_mb=128):
        self.root=Path(root).resolve();self.out=Path(out).resolve();self.limit=max_mb*1024**2
        if not self.root.is_dir():raise FileNotFoundError(self.root)
        self.out.mkdir(parents=True,exist_ok=True)
        fd,p=tempfile.mkstemp(prefix='grid-',suffix='.sqlite.tmp',dir=self.out);os.close(fd)
        self.tmp=Path(p);self.db=sqlite3.connect(p);self.catalogue=[];self.audit=[];self.owned=set();self.verified_abm={}
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
        if any('SMOKE' in x for x in rel.parts) or (rel.parts and rel.parts[0] in {'logs', 'validation'}):return False
        final_files = {'performance.tsv', 'performance.xlsx', 'performance_status.tsv', 'count_incremental_contrasts.tsv', 'historical_performance.tsv', 'sensitivity_performance.tsv', 'kinship_qc.tsv', *PAPER_FILES}
        if p.is_relative_to(self.out) and p != self.out and not (p.parent == self.out and p.name in final_files):return False
        if any(x.startswith(('.', '_')) or x in SKIP for x in rel.parts):return False
        if rel.parts and rel.parts[0] == 'architecture':
            from regional_dimensionality import completed_tables, TABLES
            if not hasattr(self, 'architecture_valid'):
                try:
                    self.architecture_valid = completed_tables(self.root) is not None
                except (OSError, ValueError, KeyError):
                    self.architecture_valid = False
            if not self.architecture_valid:return False
            if len(rel.parts) > 1 and p.name not in TABLES:return False
        pos=current_abm_position(rel.parts)
        if pos is not None:
            if any("SMOKE" in x for x in rel.parts):return False
            if len(rel.parts)>pos+3:
                parent=self.root.joinpath(*rel.parts[:pos+3])
                from final_outputs import current_abm_run
                if parent not in self.verified_abm:
                    self.verified_abm[parent]={run.name for q in parent.glob("*.CURRENT.json") if (run := current_abm_run(q)) is not None}
                current=self.verified_abm[parent]
                if rel.parts[pos+3] not in current or not (parent/rel.parts[pos+3]/"COMPLETE.json").is_file():return False
        # os.walk reaches only verified nonsymlink parents. Resolve individual
        # links/manifest targets, avoiding repeated ancestor stats on mounted disks.
        if not within_root:
            resolved=p.resolve()
            if not resolved.is_relative_to(self.root):return False
            if resolved.is_relative_to(self.out) and resolved.name not in final_files:return False
        if rel.parts and rel.parts[0] == 'abm' and len(rel.parts) > 1 and rel.parts[1] not in {'FINAL', 'ld_block_report'}:
            return False
        for module in ('prsformer', 'prsformer2'):
            if module in rel.parts and rel.parts[0] != module:return False
        if rel.parts and rel.parts[0] == 'prsformer2' and len(rel.parts) > 1:
            from benchmark_paths import prsformer2_current
            if rel.parts[1] not in TRAITS:return False
            if len(rel.parts) == 2:return True
            run = prsformer2_current(self.root, rel.parts[1])
            if run is None or rel.parts[2] != run.name or not (run/'COMPLETE.json').is_file():return False
        if rel.parts and rel.parts[0] == 'prsformer' and len(rel.parts) > 1:
            if rel.parts[1] not in TRAITS:return False
            if not (self.root/'prsformer'/rel.parts[1]/'PUBLISHED.json').is_file():return False
        if rel.parts and rel.parts[0] in ('abm_reports', 'historical') and len(rel.parts)>1:
            marker=self.root/rel.parts[0]/rel.parts[1]/'COMPLETE.json'
            if not marker.is_file() or json.loads(marker.read_text()).get('smoke'):return False
        return True
    def add(self,d,source,name,meta=None):
        d=pd.DataFrame(d)
        from final_compare import LABELS
        if "method" in d:
            d["method"] = d.method.replace(LABELS)
        try:aggregate(d)
        except ValueError as exc:
            self.note(source,'private_columns',str(exc));return
        if not len(d.columns):d=pd.DataFrame(columns=['status'])
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
        self.catalogue.append(dict(table_id=tid,sql_table=sql,trait=trait,module=module,name=name,rows=len(d),source=str(source),source_sha256=sha256(path),
                                   run_id=md.get('run_id','FINAL'),plot_kind=kind_of(d),status=state,
                                   analysis_version='FINAL',columns_json=json.dumps(list(map(str,d.columns))),
                                   note='FINAL aggregate results; no model refitting'))
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
                parts=path.relative_to(self.root).parts
                if 'twin' in parts or current_abm_position(parts) is not None:continue
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
            columns=['table_id','sql_table','trait','module','name','rows','source','source_sha256','run_id','plot_kind','status','analysis_version','columns_json','note']
            cat=pd.DataFrame(self.catalogue,columns=columns);audit=pd.DataFrame(self.audit,columns=['source','status','detail'])
            cat.to_sql('catalogue',self.db,index=False,if_exists='replace');audit.to_sql('audit',self.db,index=False,if_exists='replace')
            meta=pd.DataFrame([dict(tables=len(cat),audit_items=len(audit),root=str(self.root),read_only=True,individuals_imported=False,models_refit=False,
                completeness='partial' if len(audit) else 'available_aggregate_sources_indexed')])
            meta.to_sql('metadata',self.db,index=False,if_exists='replace');self.db.commit();self.db.close()
            with self.tmp.open('rb') as stream:os.fsync(stream.fileno())
            os.replace(self.tmp,self.out/'evidence.sqlite')
            cat.to_csv(self.out/'catalogue.tsv',sep='\t',index=False);audit.to_csv(self.out/'audit.tsv',sep='\t',index=False)
            page='<!doctype html><meta charset="utf-8"><title>GRID results</title><style>body{font:15px system-ui;max-width:1280px;margin:30px auto;padding:0 20px}table{border-collapse:collapse;width:100%}th,td{padding:7px;border-bottom:1px solid #ddd;text-align:left}details{margin:20px 0}</style><h1>GRID · Results catalogue</h1><p>Select a trait, source table and metric in Shiny. FINAL results for the completed workflow.</p>'
            page+=meta.to_html(index=False,escape=True)+cat[['trait','module','name','rows','status','analysis_version']].to_html(index=False,escape=True)
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
    print(json.dumps(Index(x.analysis_root,x.out or str(Path(x.analysis_root)/'final'),x.max_table_mb).run(x.traits.split(','))))
if __name__=='__main__':main()
