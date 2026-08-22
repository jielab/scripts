#!/usr/bin/env python3
"""Normalize GU methods into a common, UKB-friendly segment database.

The database is the primary artifact for large cohorts. Exact raw tract boundaries are
kept, while a cohort-level canonical segment_code is assigned by reciprocal-overlap
clustering. This gives each participant 0-N stable-ish tract codes analogous to an
ICD event table, while preserving the original calls for sensitivity analyses.
"""
from __future__ import annotations
import argparse, gzip, hashlib, math, re, sqlite3, urllib.parse
from pathlib import Path
from typing import Iterable
import pandas as pd

SEG_COLS=['sample_id','method','source','source_class','chr','start','end','length_bp','haplotype','score','posterior','trait','locus_id','genome_build','batch_id','raw_file']

def chrom_norm(x):
    s=str(x); s=re.sub(r'^chr','',s,flags=re.I); s=re.sub(r'_part\d+$','',s,flags=re.I)
    return 'X' if s=='23' else s

def source_class(x):
    s=str(x).lower()
    if any(k in s for k in ('altai','vindija','chagyr','neander')): return 'Neanderthal'
    if 'denis' in s: return 'Denisovan'
    if 'mosaic' in s: return 'Mosaic'
    if 'trace' in s or 'ghost' in s or 'unknown' in s: return 'Ghost/Unknown'
    return 'Archaic/Other'

def empty_segments(): return pd.DataFrame(columns=SEG_COLS)

def finish(df,method,build,raw,batch=''):
    if df.empty:return empty_segments()
    for c in SEG_COLS:
        if c not in df: df[c]=pd.NA
    df['method']=method; df['genome_build']=build; df['raw_file']=str(raw); df['batch_id']=batch
    df['chr']=df['chr'].map(chrom_norm)
    df['start']=pd.to_numeric(df['start'],errors='coerce'); df['end']=pd.to_numeric(df['end'],errors='coerce')
    df=df[df.start.notna() & df.end.notna() & (df.end>df.start)].copy()
    df['start']=df.start.astype('int64'); df['end']=df.end.astype('int64'); df['length_bp']=(df.end-df.start).astype('int64')
    df['source']=df['source'].fillna('unknown').astype(str); df['source_class']=df['source'].map(source_class)
    return df[SEG_COLS]

def find_first(root,patterns):
    out=[]
    for p in patterns: out.extend(root.glob(p))
    out=[x for x in out if x.is_file() and x.stat().st_size>0]
    return sorted(set(out))

def batch_from_path(f, method):
    """Infer batch ID from UKB batch layout; return empty string for ordinary single-run outputs."""
    parts=list(f.parts)
    try:
        i=parts.index('ukb')
        if i+2 < len(parts) and parts[i+1] == method:
            return parts[i+2]
    except ValueError:
        pass
    return ''


def iter_ibdmix(root,build):
    # Prefer final caller output. Fall back to report output only when no final calls exist.
    files=find_first(root,[
        "ibdmix/final/all_archaic_refs*.segments.tsv.gz",
        "ukb/ibdmix/**/final/all_archaic_refs*.segments.tsv.gz",
    ])
    if not files:
        files=find_first(root,[
            "ibdmix/report/genomewide_segments.reduced.tsv.gz",
            "ibdmix/report/genomewide_segments.tsv.gz",
            "ukb/ibdmix/**/report/genomewide_segments.reduced.tsv.gz",
            "ukb/ibdmix/**/report/genomewide_segments.tsv.gz",
        ])
    for f in files:
        batch=batch_from_path(f,'ibdmix')
        for d in pd.read_csv(f,sep='\t',compression='infer',chunksize=200000):
            ren={}
            for a,b in [('ID','sample_id'),('chrom','chr'),('length','length_bp'),('slod','score'),('anc','source')]:
                if a in d: ren[a]=b
            d=d.rename(columns=ren)
            if 'sample_id' not in d or 'chr' not in d: continue
            d['haplotype']=pd.NA; d['posterior']=pd.NA; d['trait']=pd.NA; d['locus_id']=pd.NA
            if 'source' not in d:d['source']='IBDmix_archaic'
            yield finish(d,'ibdmix',build,f,batch=batch)


def iter_trace(root,build):
    files=find_first(root,[
        "trace/final/trace_haplotype_segments.tsv.gz",
        "ukb/trace/**/final/trace_haplotype_segments.tsv.gz",
    ])
    for f in files:
        batch=batch_from_path(f,'trace')
        for d in pd.read_csv(f,sep='\t',compression='infer',chunksize=200000):
            d=d.rename(columns={'sample':'sample_id','chrom':'chr','mean_posterior':'posterior'})
            if 'chr' not in d and 'chromosome' in d:d=d.rename(columns={'chromosome':'chr'})
            d['source']='TRACE_ghost_or_unknown'; d['score']=pd.NA; d['trait']=pd.NA; d['locus_id']=pd.NA
            yield finish(d,'trace',build,f,batch=batch)

def parse_as3_file(f,build):
    # Official AS3.1-Mamba output is BED-like. Accept headered or headerless variants.
    first=''
    with open(f,'rt',errors='replace') as h:
        for line in h:
            if line.strip() and not line.startswith('#'):
                first=line.rstrip('\n'); break
    toks=[x.strip().lower() for x in first.split('\t')]
    has_header=any(x in {'chr','chrom','chromosome'} for x in toks) and any(x=='start' for x in toks)
    d=pd.read_csv(f,sep='\t',comment='#',header=0 if has_header else None)
    low={str(c).lower():c for c in d.columns}
    def col(*names):
        for n in names:
            if n.lower() in low:return low[n.lower()]
        return None
    cc=col('Chr','chrom','chr'); ss=col('Start','start'); ee=col('End','end')
    if cc is None and d.shape[1]>=3:
        d.columns=list(range(d.shape[1])); cc,ss,ee=0,1,2
    if cc is None:return empty_segments()
    out=pd.DataFrame({'chr':d[cc],'start':d[ss],'end':d[ee]})
    hc=col('Haplotype','haplotype'); ac=col('Archaic','archaic'); sc=col('Score','score'); idc=col('SampleID_HapID','sampleid_hapid','sample','SampleID')
    # Headerless official order: Chr Start End Haplotype Archaic #SNP Score #SNP_A1 #SNP_A2 SampleID_HapID
    if isinstance(cc,int) and d.shape[1]>=10:
        hc,ac,sc,idc=3,4,6,9
    out['haplotype']=d[hc].astype(str) if hc is not None else pd.NA
    archaic=d[ac].astype(str) if ac is not None else pd.Series(['unknown']*len(d))
    amap={'1':'Denisovan','2':'Neanderthal','3':'Mosaic','1.0':'Denisovan','2.0':'Neanderthal','3.0':'Mosaic'}
    out['source']=archaic.map(lambda x: amap.get(str(x),str(x)))
    out['score']=pd.to_numeric(d[sc],errors='coerce') if sc is not None else pd.NA
    sid=d[idc].astype(str) if idc is not None else pd.Series(['']*len(d))
    def split_id(x):
        m=re.match(r'^(.*?)(?:[_:.-](?:hap)?)([12])$',x,re.I)
        return (m.group(1),m.group(2)) if m else (x,None)
    sp=[split_id(x) for x in sid]
    out['sample_id']=[x[0] for x in sp]
    out.loc[out.sample_id.eq(''),'sample_id']=pd.NA
    if hc is None: out['haplotype']=[x[1] for x in sp]
    out['posterior']=pd.NA; out['trait']=pd.NA; out['locus_id']=pd.NA
    return finish(out,'as3',build,f,batch=f.parent.name)

def iter_as3(root,build):
    files=find_first(root,[
        "as3/**/introgression_prediction.bed",
        "ukb/as3/**/introgression_prediction.bed",
    ])
    for f in files:
        d=parse_as3_file(f,build)
        if not d.empty:
            batch=batch_from_path(f,'as3')
            if batch:
                d['batch_id']=batch
            yield d

def loci_tables(root,method):
    base=root/method
    patterns=['**/report/inherited_region.tsv'] if method=='loci_avcf' else ['**/report/asnp_region_summary.tsv','**/report/asnp_inherited_segments.tsv']
    seen=set()
    for pat in patterns:
        for f in base.glob(pat):
            if not f.is_file() or f.stat().st_size==0:continue
            try:d=pd.read_csv(f,sep='\t')
            except Exception:continue
            if d.empty:continue
            carrier_counts={}
            if method=='loci_avcf':
                sm=f.parent/'haplotype_sample_map.tsv'
                if sm.is_file() and sm.stat().st_size:
                    try:
                        x=pd.read_csv(sm,sep='\t')
                        idc='id' if 'id' in x.columns else ('locus_id' if 'locus_id' in x.columns else None)
                        sc='sample' if 'sample' in x.columns else ('sample_id' if 'sample_id' in x.columns else None)
                        if idc and sc:
                            carrier_counts=x.dropna(subset=[idc,sc]).groupby(idc)[sc].nunique().to_dict()
                    except Exception: pass
            for _,r in d.iterrows():
                chrom=r.get('lead_chr',r.get('CHR',r.get('chr','')))
                st=r.get('core_start',r.get('region_start',r.get('start',r.get('POS',None))))
                en=r.get('core_end',r.get('region_end',r.get('end',r.get('POS',None))))
                try: st=int(float(st)); en=int(float(en))
                except Exception: continue
                if en<=st: en=st+1
                key=(method,str(r.get('trait','')),str(r.get('id',r.get('lead_snp',''))),chrom,st,en)
                if key in seen:continue
                seen.add(key)
                source=r.get('matched_archaics',r.get('best_lineage',r.get('asnp_lineage','archaic_map')))
                lid=r.get('id',r.get('lead_snp',None))
                yield {'method':method,'trait':r.get('trait',None),'locus_id':lid,'chr':chrom_norm(chrom),'start':st,'end':en,'source':str(source),'status':r.get('status',r.get('strict_filter_reason',None)),'n_carriers':carrier_counts.get(lid,None),'raw_file':str(f)}

def init_db(con):
    con.executescript('''
    PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL; PRAGMA temp_store=FILE;
    DROP TABLE IF EXISTS segments_raw; DROP TABLE IF EXISTS loci; DROP TABLE IF EXISTS coord_catalog;
    DROP TABLE IF EXISTS segment_catalog; DROP TABLE IF EXISTS segments; DROP TABLE IF EXISTS carriers; DROP TABLE IF EXISTS sample_burden;
    DROP TABLE IF EXISTS overview_metrics; DROP TABLE IF EXISTS overview_summary; DROP TABLE IF EXISTS segment_chr_summary; DROP TABLE IF EXISTS segment_density_1mb;
    CREATE TABLE segments_raw(sample_id TEXT,method TEXT,source TEXT,source_class TEXT,chr TEXT,start INTEGER,end INTEGER,length_bp INTEGER,haplotype TEXT,score REAL,posterior REAL,trait TEXT,locus_id TEXT,genome_build TEXT,batch_id TEXT,raw_file TEXT);
    CREATE TABLE loci(method TEXT,trait TEXT,locus_id TEXT,chr TEXT,start INTEGER,end INTEGER,source TEXT,status TEXT,n_carriers INTEGER,raw_file TEXT);
    CREATE TABLE coord_catalog(source_class TEXT,chr TEXT,start INTEGER,end INTEGER,segment_code TEXT);
    CREATE TABLE segment_catalog(segment_code TEXT PRIMARY KEY,source_class TEXT,chr TEXT,start INTEGER,end INTEGER,length_bp INTEGER,n_distinct_boundaries INTEGER,genome_build TEXT);
    ''')

def insert_df(con,table,df):
    if df is None or df.empty:return 0
    df.to_sql(table,con,if_exists='append',index=False,method='multi',chunksize=1000); return len(df)

def reciprocal(a,b,c,d):
    ov=max(0,min(b,d)-max(a,c));
    if ov<=0:return 0.0
    return min(ov/max(1,b-a),ov/max(1,d-c))

def code_for(build,src,chrom,start,end):
    tag={'Neanderthal':'NEA','Denisovan':'DEN','Mosaic':'MOS','Ghost/Unknown':'GST','Archaic/Other':'ARC'}.get(src,'ARC')
    dig=hashlib.sha1(f'{build}|{src}|{chrom}|{round(start,-3)}|{round(end,-3)}'.encode()).hexdigest()[:10].upper()
    return f'GU{build.replace("GRCh","")}-{tag}-{chrom}-{dig}'

def build_catalog(con,build,threshold):
    groups=con.execute('SELECT DISTINCT source_class,chr FROM segments_raw ORDER BY source_class, CASE WHEN chr GLOB "[0-9]*" THEN CAST(chr AS INT) ELSE 99 END, chr').fetchall()
    ins_coord=[]; ins_cat=[]
    for src,chrom in groups:
        cur=con.execute('SELECT DISTINCT start,end FROM segments_raw WHERE source_class=? AND chr=? ORDER BY start,end',(src,chrom))
        coords=cur.fetchall()
        cluster=[]; mean_s=mean_e=None
        def flush():
            nonlocal cluster,mean_s,mean_e,ins_coord,ins_cat
            if not cluster:return
            rs=int(round(mean_s)); re=int(round(mean_e)); code=code_for(build,src,chrom,rs,re)
            ins_cat.append((code,src,chrom,rs,re,re-rs,len(cluster),build))
            ins_coord.extend((src,chrom,s,e,code) for s,e in cluster)
            if len(ins_coord)>=100000:
                con.executemany('INSERT INTO coord_catalog VALUES(?,?,?,?,?)',ins_coord); ins_coord=[]
            if len(ins_cat)>=10000:
                con.executemany('INSERT OR REPLACE INTO segment_catalog VALUES(?,?,?,?,?,?,?,?)',ins_cat); ins_cat=[]
            cluster=[]; mean_s=mean_e=None
        for s,e in coords:
            if not cluster:
                cluster=[(s,e)]; mean_s=float(s); mean_e=float(e); continue
            if reciprocal(s,e,mean_s,mean_e) >= threshold:
                n=len(cluster); cluster.append((s,e)); mean_s=(mean_s*n+s)/(n+1); mean_e=(mean_e*n+e)/(n+1)
            else:
                flush(); cluster=[(s,e)]; mean_s=float(s); mean_e=float(e)
        flush()
    if ins_coord: con.executemany('INSERT INTO coord_catalog VALUES(?,?,?,?,?)',ins_coord)
    if ins_cat: con.executemany('INSERT OR REPLACE INTO segment_catalog VALUES(?,?,?,?,?,?,?,?)',ins_cat)
    con.commit()
    con.executescript('''
      CREATE INDEX idx_coord ON coord_catalog(source_class,chr,start,end);
      CREATE TABLE segments AS SELECT s.*, c.segment_code FROM segments_raw s JOIN coord_catalog c USING(source_class,chr,start,end);
      CREATE INDEX idx_segments_sample ON segments(sample_id); CREATE INDEX idx_segments_region ON segments(chr,start,end); CREATE INDEX idx_segments_code ON segments(segment_code);
      CREATE INDEX idx_segments_method_source ON segments(method,source_class,chr,start);
      CREATE TABLE carriers AS
        SELECT sample_id,segment_code,source_class,
               MIN(2, CASE WHEN COUNT(DISTINCT CASE WHEN haplotype NOT IN ('','nan','None') THEN haplotype END)>0 THEN COUNT(DISTINCT CASE WHEN haplotype NOT IN ('','nan','None') THEN haplotype END) ELSE 1 END) AS dosage,
               COUNT(DISTINCT method) AS n_methods,
               GROUP_CONCAT(DISTINCT method) AS methods_support,
               MAX(score) AS max_score, MAX(posterior) AS max_posterior
        FROM segments WHERE sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None') GROUP BY sample_id,segment_code,source_class;
      CREATE INDEX idx_carriers_sample ON carriers(sample_id); CREATE INDEX idx_carriers_code ON carriers(segment_code);
      CREATE TABLE sample_burden AS
        SELECT sample_id,method,source_class,COUNT(*) AS n_segments,SUM(length_bp) AS total_bp,AVG(length_bp) AS mean_length_bp,MAX(score) AS max_score,MAX(posterior) AS max_posterior
        FROM segments WHERE sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None') GROUP BY sample_id,method,source_class;
      CREATE INDEX idx_burden_sample ON sample_burden(sample_id);
      CREATE TABLE overview_metrics(metric TEXT PRIMARY KEY,value INTEGER);
      INSERT INTO overview_metrics VALUES
        ('n_samples',(SELECT COUNT(DISTINCT sample_id) FROM segments WHERE sample_id IS NOT NULL)),
        ('n_segments',(SELECT COUNT(*) FROM segments)),
        ('n_codes',(SELECT COUNT(*) FROM segment_catalog)),
        ('n_methods',(SELECT COUNT(*) FROM (SELECT method FROM segments UNION SELECT method FROM loci)));
      CREATE TABLE overview_summary AS
        SELECT 'individual_segments' layer,method,source_class source,COUNT(*) n_records,
               COUNT(DISTINCT sample_id) n_samples,SUM(length_bp) total_bp
        FROM segments GROUP BY method,source_class;
      INSERT INTO overview_summary
        SELECT 'locus_evidence',method,source,COUNT(*),NULL,SUM(end-start)
        FROM loci GROUP BY method,source;
      CREATE TABLE segment_chr_summary AS
        SELECT method,source_class source,chr,COUNT(*) n_records,
               COUNT(DISTINCT sample_id) n_samples,SUM(length_bp) total_bp,AVG(length_bp) mean_length_bp
        FROM segments GROUP BY method,source_class,chr;
      CREATE TABLE segment_density_1mb AS
        SELECT method,source_class source,chr,CAST(start/1000000 AS INTEGER)*1000000 bin_start,
               COUNT(*) n_calls,COUNT(DISTINCT sample_id) n_carriers,SUM(length_bp) total_bp
        FROM segments GROUP BY method,source_class,chr,CAST(start/1000000 AS INTEGER);
      CREATE INDEX idx_density_method_source ON segment_density_1mb(method,source,chr,bin_start);
    ''')
    con.commit()

def export_query(con,query,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    cur=con.execute(query); cols=[x[0] for x in cur.description]
    with gzip.open(path,'wt') as h:
        h.write('\t'.join(cols)+'\n')
        while True:
            rows=cur.fetchmany(50000)
            if not rows:break
            for r in rows:h.write('\t'.join('' if x is None else str(x) for x in r)+'\n')

def export_browser_loci(con, outdir, build):
    """Precompute browser-ready BED and stable IGV/UCSC links for every locus."""
    db='hg38' if '38' in str(build) else 'hg19'
    rows=con.execute('SELECT method,trait,locus_id,chr,start,end,source,status FROM loci ORDER BY chr,start,end').fetchall()
    bed=outdir/'loci.browser.bed'; links=outdir/'loci.browser_links.tsv'
    with bed.open('w',encoding='utf-8') as bh, links.open('w',encoding='utf-8') as lh:
        bh.write('track name="GU_loci" description="GU archaic-introgression loci" itemRgb="On"\n')
        lh.write('method\ttrait\tlocus_id\tgenome\tlocus\tigv_url\tucsc_url\n')
        for method,trait,lid,chrom,start,end,source,status in rows:
            c='chr'+re.sub(r'^chr','',str(chrom),flags=re.I); start0=max(0,int(start)-1); end=int(end)
            name='|'.join('' if x is None else str(x) for x in (trait,lid,source))
            bh.write(f'{c}\t{start0}\t{end}\t{name}\t0\t.\t{start0}\t{end}\t36,113,163\n')
            flank=250000; locus=f'{c}:{max(1,int(start)-flank)}-{end+flank}'; enc=urllib.parse.quote(locus,safe='')
            igv=f'https://igv.org/app/?genome={db}&locus={enc}'
            ucsc=f'https://genome.ucsc.edu/cgi-bin/hgTracks?db={db}&position={enc}'
            lh.write('\t'.join('' if x is None else str(x) for x in (method,trait,lid,db,locus,igv,ucsc))+'\n')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--analysis-root',required=True,type=Path); ap.add_argument('--output-dir',required=True,type=Path); ap.add_argument('--build',default='GRCh37'); ap.add_argument('--reciprocal-overlap',type=float,default=.5); args=ap.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True); db=args.output_dir/'gu.sqlite'
    if db.exists(): db.unlink()
    con=sqlite3.connect(db); init_db(con); counts={}
    for name,it in [('ibdmix',iter_ibdmix(args.analysis_root,args.build)),('trace',iter_trace(args.analysis_root,args.build)),('as3',iter_as3(args.analysis_root,args.build))]:
        n=0
        if it is not None:
            for d in it:n+=insert_df(con,'segments_raw',d)
        counts[name]=n; con.commit()
    loci=[]
    for method in ('loci_avcf','loci_asnp'): loci.extend(list(loci_tables(args.analysis_root,method)))
    if loci: insert_df(con,'loci',pd.DataFrame(loci))
    con.execute('CREATE INDEX idx_loci_region ON loci(chr,start,end)'); con.commit()
    build_catalog(con,args.build,args.reciprocal_overlap)
    export_query(con,'SELECT * FROM segments ORDER BY chr,start,end,sample_id',args.output_dir/'segments.tsv.gz')
    export_query(con,'SELECT * FROM segment_catalog ORDER BY chr,start,end',args.output_dir/'segment_catalog.tsv.gz')
    export_query(con,'SELECT * FROM carriers ORDER BY sample_id,segment_code',args.output_dir/'carriers.tsv.gz')
    export_query(con,'SELECT * FROM sample_burden ORDER BY sample_id,method,source_class',args.output_dir/'sample_burden.tsv.gz')
    export_query(con,'SELECT * FROM loci ORDER BY chr,start,end',args.output_dir/'loci.tsv.gz')
    export_browser_loci(con,args.output_dir,args.build)
    con.close()
    print('Normalized:',counts,'loci=',len(loci),'database=',db)

if __name__=='__main__':main()
