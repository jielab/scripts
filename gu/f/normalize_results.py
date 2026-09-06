#!/usr/bin/env python3
"""Normalize GU methods into a common, UKB-friendly segment database.

The database is the primary artifact for large cohorts. Exact raw tract boundaries are
kept, while a cohort-level canonical segment_code is assigned by reciprocal-overlap
clustering. This gives each participant 0-N stable-ish tract codes analogous to an
ICD event table, while preserving the original calls for sensitivity analyses.
"""
from __future__ import annotations
import argparse, bisect, gzip, hashlib, math, os, re, shutil, sqlite3, tempfile
from pathlib import Path
import pandas as pd
from comm import resolve_tsv_path

SEG_COLS=['dataset_id','sample_id','method','source','source_class','chr','start','end','length_bp','haplotype','method_haplotype_index','score','posterior','trait','locus_id','genome_build','batch_id','raw_file']
REFERENCE_COLS=['dataset_id','population','genome_build','chr','start','end','source_class','reference_role','raw_file']
REFERENCE_RAW_COLS=REFERENCE_COLS+['reference_sample_id','reference_super_population']
POPS_1KG={
    'ACB','ASW','ESN','GWD','LWK','MSL','YRI','CLM','MXL','PEL','PUR',
    'CDX','CHB','CHS','JPT','KHV','CEU','FIN','GBR','IBS','TSI',
    'BEB','GIH','ITU','PJL','STU'
}

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

def build_norm(x, fallback):
    if pd.isna(x) or str(x).strip() in ('', 'nan', 'None'): x=fallback
    s=str(x).strip()
    if re.fullmatch(r'b?3[78]',s,re.I): return f'GRCh{s[-2:]}'
    if re.fullmatch(r'GRCh3[78]',s,re.I): return f'GRCh{s[-2:]}'
    return s

def dataset_from_path(path,method,default='1kg'):
    """Resolve the target namespace from GU's <method>/<target>/<scope> layout."""
    parts=list(Path(path).parts)
    if 'ukb' in [x.lower() for x in parts]:return 'ukb'
    indexes=[i for i,x in enumerate(parts) if x.lower()==method.lower()]
    if not indexes:return default
    i=indexes[-1]
    if i+2>=len(parts):
        raise ValueError(f'target-aware GU path required for {method}: {path}')
    if parts[i+2].lower() in {'final','loci','segments','results','raw','calls','samples','log'}:
        raise ValueError(f'target namespace missing from GU path: {path}')
    return parts[i+1]

def finish(df,method,build,raw,batch='',dataset=None):
    if df.empty:return empty_segments()
    for c in SEG_COLS:
        if c not in df: df[c]=pd.NA
    df['dataset_id']=dataset or dataset_from_path(raw,method)
    df['method']=method; df['genome_build']=df['genome_build'].map(lambda x:build_norm(x,build)); df['raw_file']=str(raw); df['batch_id']=batch
    df['chr']=df['chr'].map(chrom_norm)
    df['start']=pd.to_numeric(df['start'],errors='coerce'); df['end']=pd.to_numeric(df['end'],errors='coerce')
    df=df[df.start.notna() & df.end.notna() & (df.end>df.start)].copy()
    df['start']=df.start.astype('int64'); df['end']=df.end.astype('int64'); df['length_bp']=(df.end-df.start).astype('int64')
    df['source']=df['source'].fillna('unknown').astype(str); df['source_class']=df['source'].map(source_class)
    return df[SEG_COLS]

def is_transient_output(path):
    """Return True for AS3 atomic-work directories such as .chr22.run.<pid>."""
    return any(
        part.startswith('.') and ('.run.' in part or '.old.' in part)
        for part in path.parts
    )

def find_first(root,patterns):
    out=[]
    for p in patterns: out.extend(root.glob(p))
    out=[x for x in out if not is_transient_output(x) and x.is_file() and x.stat().st_size>0]
    return sorted(set(out))

CORE_FILE_PATTERNS={
    'ibdmix':(),
    'trace':(),
    'as3':(),
    'phyml':(
        '**/final/*.tsv','**/final/evidence_parameters.json','**/final/*unfiltered.tsv.gz','**/loci/**/sites.tsv','**/loci/**/search_sites.tsv',
        '**/loci/**/ld.tsv','**/loci/**/selected_region.tsv','**/loci/**/archaic.tsv',
        '**/loci/**/ancestral.tsv',
        '**/loci/**/haplotypes.tsv','**/loci/**/haplotypes.phy','**/loci/**/haplotypes.phy.meta.tsv',
        '**/loci/**/*_phyml_tree.txt','**/loci/**/*_phyml_tree.png','**/loci/**/*_phyml_stats.txt',
        '**/run.meta.tsv','**/request.loci.analysis.bed','**/request.loci.core.bed','**/request.loci.map.tsv',
    ),
}

def package_core_results(analysis_root,output_dir):
    """Copy the small, reusable result subset used by Shiny into normalize/<method>."""
    analysis_root=analysis_root.resolve(); output_dir=output_dir.resolve()
    stage=output_dir/f'.core.part.{os.getpid()}'
    if stage.exists():shutil.rmtree(stage)
    copied={}
    try:
        for method,patterns in CORE_FILE_PATTERNS.items():
            dst_root=stage/method; dst_root.mkdir(parents=True,exist_ok=True); copied[method]=0
            for src_root,prefix in ((analysis_root/method,Path()),(analysis_root/'ukb'/method,Path('ukb'))):
                files=[]
                if src_root.is_dir():
                    for pattern in patterns:files.extend(src_root.glob(pattern))
                files=sorted(set(resolve_tsv_path(x) for x in files if x.is_file() and not is_transient_output(x) and x.stat().st_size>0))
                for src in files:
                    dst=dst_root/prefix/src.relative_to(src_root)
                    dst.parent.mkdir(parents=True,exist_ok=True)
                    if src.name.endswith('unfiltered.tsv'):
                        # Legacy runs are packaged in the same compressed format
                        # as new runs, without changing the original source.
                        dst=dst.with_name(dst.name+'.gz')
                        with src.open('rb') as reader, gzip.open(dst,'wb') as writer:
                            shutil.copyfileobj(reader,writer,length=1024*1024)
                    else:
                        shutil.copy2(src,dst)
                copied[method]+=len(files)
        for method in CORE_FILE_PATTERNS:
            dst=output_dir/method; backup=output_dir/f'.{method}.old.{os.getpid()}'
            if backup.exists():shutil.rmtree(backup)
            if dst.exists():os.replace(dst,backup)
            try:os.replace(stage/method,dst)
            except Exception:
                if backup.exists() and not dst.exists():os.replace(backup,dst)
                raise
            if backup.exists():shutil.rmtree(backup)
    finally:
        if stage.exists():shutil.rmtree(stage)
    return copied

def portable_paths(df,package_root,columns):
    """Store paths inside the share package relative to normalize/, leaving external provenance intact."""
    root=package_root.resolve()
    for column in columns:
        if column not in df:continue
        def convert(value):
            if pd.isna(value) or not str(value):return value
            try:return Path(str(value)).resolve().relative_to(root).as_posix()
            except ValueError:return str(value)
        df[column]=df[column].map(convert)
    return df

def phyml_artifact_in_run(value,run_dir):
    """Map an absolute path recorded by the raw run to its packaged counterpart."""
    if value is None or pd.isna(value) or not str(value):return ''
    path=Path(str(value)); parts=path.parts
    for marker in ('loci','final'):
        indexes=[i for i,x in enumerate(parts) if x==marker]
        if indexes:
            candidate=run_dir.joinpath(*parts[indexes[-1]:])
            if candidate.is_file():return str(candidate)
    hits=list((run_dir/'loci').rglob(path.name)) if (run_dir/'loci').is_dir() else []
    return str(hits[0]) if len(hits)==1 else str(value)

def batch_from_path(f, method):
    """Infer batch ID from UKB batch layout; return empty string for ordinary single-run outputs."""
    parts=list(f.parts)
    try:
        i=parts.index('ukb')
        if i+2 < len(parts) and parts[i+1] == method:
            return parts[i+2]
    except ValueError:
        pass
    try:
        i=parts.index(method)
        if i+2 < len(parts) and parts[i+1] == 'ukb':
            return parts[i+2]
    except ValueError:
        pass
    return ''


def iter_ibdmix(root,build):
    files=find_first(root,[
        "ibdmix/**/final/all_archaic_refs*.segments.tsv.gz",
        "ukb/ibdmix/**/final/all_archaic_refs*.segments.tsv.gz",
    ])
    for f in files:
        batch=batch_from_path(f,'ibdmix')
        for d in pd.read_csv(f,sep='\t',compression='infer',chunksize=200000):
            ren={}
            for a,b in [('ID','sample_id'),('chrom','chr'),('length','length_bp'),('slod','score'),('anc','source')]:
                if a in d: ren[a]=b
            d=d.rename(columns=ren)
            if 'sample_id' not in d or 'chr' not in d: continue
            d['haplotype']=pd.NA; d['posterior']=pd.NA; d['trait']=pd.NA
            if 'locus_id' not in d:d['locus_id']=pd.NA
            if 'source' not in d:d['source']='IBDmix_archaic'
            yield finish(d,'ibdmix',build,f,batch=batch,dataset=dataset_from_path(f,'ibdmix'))


def iter_trace(root,build):
    files=find_first(root,[
        "trace/**/final/trace_haplotype_segments.tsv.gz",
        "ukb/trace/**/final/trace_haplotype_segments.tsv.gz",
    ])
    for f in files:
        batch=batch_from_path(f,'trace')
        for d in pd.read_csv(f,sep='\t',compression='infer',chunksize=200000):
            d=d.rename(columns={'sample':'sample_id','chrom':'chr','mean_posterior':'posterior'})
            if 'chr' not in d and 'chromosome' in d:d=d.rename(columns={'chromosome':'chr'})
            d['source']='TRACE_ghost_or_unknown'; d['score']=pd.NA; d['trait']=pd.NA
            if 'locus_id' not in d:d['locus_id']=pd.NA
            yield finish(d,'trace',build,f,batch=batch,dataset=dataset_from_path(f,'trace'))

def find_parent_file(f,name):
    for parent in [f.parent,*f.parents]:
        candidate=parent/name
        if candidate.is_file() and candidate.stat().st_size:return candidate
    return None

def as3_build(f,fallback):
    meta=find_parent_file(f,'run.meta.tsv')
    if not meta:return fallback
    try:
        d=pd.read_csv(meta,sep='\t',dtype=str)
        row=d.loc[d.iloc[:,0].eq('genome_build')]
        return build_norm(row.iloc[0,1],fallback) if not row.empty else fallback
    except Exception:return fallback

def clip_to_loci(df,bed):
    if not bed or df.empty:return df
    try:loci=pd.read_csv(bed,sep='\t',comment='#',header=None,usecols=[0,1,2,3],names=['chr','locus_start','locus_end','clip_locus_id'])
    except Exception:return df
    loci['chr']=loci['chr'].map(chrom_norm); loci['locus_start']=pd.to_numeric(loci.locus_start,errors='coerce'); loci['locus_end']=pd.to_numeric(loci.locus_end,errors='coerce')
    z=df.copy(); z['chr']=z['chr'].map(chrom_norm); z['start']=pd.to_numeric(z.start,errors='coerce'); z['end']=pd.to_numeric(z.end,errors='coerce')
    z=z.merge(loci,on='chr',how='inner'); z=z[(z.end>z.locus_start)&(z.start<z.locus_end)].copy()
    if z.empty:return z
    z['start']=z[['start','locus_start']].max(axis=1); z['end']=z[['end','locus_end']].min(axis=1); z['locus_id']=z['clip_locus_id']
    return z.drop(columns=['locus_start','locus_end','clip_locus_id'])

def parse_as3_file(f,build):
    # Official AS3.1-Mamba output is BED-like but Start/End are VCF positions:
    # both are 1-based inclusive. Convert once to the internal 0-based half-open
    # convention so interval lengths and clipping agree with every other method.
    first=''; skiprows=0
    with open(f,'rt',errors='replace') as h:
        for line in h:
            if line.strip() and not line.startswith('#'):
                first=line.rstrip('\n'); break
            skiprows+=1
    toks=[x.strip().lower() for x in first.split('\t')]
    has_header=any(x in {'chr','chrom','chromosome'} for x in toks) and any(x=='start' for x in toks)
    # Do not use pandas comment='#': the official header names include #SNP,
    # #SNP_A1, and #SNP_A2, which would otherwise truncate the header itself.
    d=pd.read_csv(f,sep='\t',skiprows=skiprows,header=0 if has_header else None)
    low={str(c).lower():c for c in d.columns}
    def col(*names):
        for n in names:
            if n.lower() in low:return low[n.lower()]
        return None
    cc=col('Chr','chrom','chr'); ss=col('Start','start'); ee=col('End','end')
    if cc is None and d.shape[1]>=3:
        d.columns=list(range(d.shape[1])); cc,ss,ee=0,1,2
    if cc is None:return empty_segments()
    out=pd.DataFrame({'chr':d[cc],'start':pd.to_numeric(d[ss],errors='coerce')-1,'end':pd.to_numeric(d[ee],errors='coerce')})
    hc=col('Haplotype','haplotype'); ac=col('Archaic','archaic'); sc=col('Score','score'); idc=col('SampleID_HapID','sampleid_hapid','sample','SampleID')
    # Headerless official order: Chr Start End Haplotype Archaic #SNP Score #SNP_A1 #SNP_A2 SampleID_HapID
    if isinstance(cc,int) and d.shape[1]>=10:
        hc,ac,sc,idc=3,4,6,9
    # Derive stable haplotype 1/2 from SampleID_HapID.
    out['method_haplotype_index']=d[hc].astype(str) if hc is not None else pd.NA
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
    out['haplotype']=[x[1] for x in sp]
    out['posterior']=pd.NA; out['trait']=pd.NA; out['locus_id']=pd.NA
    out=clip_to_loci(out,find_parent_file(f,'input.loci.bed'))
    return finish(out,'as3',as3_build(f,build),f,batch='',dataset=dataset_from_path(f,'as3'))

def iter_as3(root,build):
    files=find_first(root,[
        "as3/**/introgression.bed",
        "ukb/as3/**/introgression.bed",
    ])
    for f in files:
        d=parse_as3_file(f,build)
        if not d.empty:
            batch=batch_from_path(f,'as3')
            if batch:
                d['batch_id']=batch
            yield d


def read_run_meta(path):
    values={}
    if not path or not Path(path).is_file():return values
    try:
        with Path(path).open(errors='replace') as handle:
            for line in handle:
                fields=line.rstrip('\n').split('\t',1)
                if len(fields)==2 and fields[0] not in values:values[fields[0]]=fields[1]
    except OSError:pass
    return values


def discover_method_runs(root,fallback_build):
    """Record completed chromosome units, including valid zero-call runs.

    Segment presence is not a completion signal: a successful caller can emit
    only a header.  These rows let Shiny distinguish a true negative result
    from a method that was never run (or is unsupported on chrX).
    """
    rows={}
    def add(method,path,build,chrom,status='complete'):
        chrom=chrom_norm(chrom); build=build_norm(build,fallback_build)
        if not chrom:return
        eligible=1; note='completed native analysis'
        if chrom=='X' and method=='ibdmix':
            eligible=0; note='experimental male-X pseudo-diploid adaptation; excluded from GU evidence score'
        elif chrom=='X' and method=='as3':
            eligible=0; note='AS3 Ref1028/model is autosomes-only; chrX is unsupported'
        elif method=='trace':
            note='haplotype-node call; validity depends on the full-chromosome ARG backend and calibration'
        key=(dataset_from_path(path,method),build,chrom,method)
        rows[key]=(*key,status,eligible,note,str(path))

    for marker in sorted(root.glob('ibdmix/**/.complete'))+sorted(root.glob('ukb/ibdmix/**/.complete')):
        meta=marker.parent/'run.meta.tsv'; values=read_run_meta(meta); text=meta.read_text(errors='replace') if meta.is_file() else ''
        build=values.get('genome_build',fallback_build)
        chroms=set(re.findall(r'(?:^|[/\t_.-])chr([0-9]+|X)(?=$|[/\t_.-])',text,flags=re.I|re.M))
        if not chroms:
            chroms=set(re.findall(r'chr([0-9]+|X)',str(marker.parent),flags=re.I))
        for chrom in chroms:add('ibdmix',meta if meta.is_file() else marker,build,chrom)

    trace_files=find_first(root,["trace/**/final/trace_haplotype_segments.tsv.gz","ukb/trace/**/final/trace_haplotype_segments.tsv.gz"])
    for output in trace_files:
        meta=find_parent_file(output,'run.meta.tsv'); values=read_run_meta(meta)
        for chrom in re.split(r'[, ]+',values.get('chromosomes','').strip()):
            if chrom:add('trace',meta or output,values.get('build',fallback_build),chrom)

    for marker in sorted(root.glob('as3/**/results/**/chr*/.complete'))+sorted(root.glob('ukb/as3/**/results/**/chr*/.complete')):
        meta=find_parent_file(marker,'run.meta.tsv'); values=read_run_meta(meta)
        match=re.search(r'^chr([0-9]+|X)',marker.parent.name,flags=re.I)
        if match:add('as3',meta or marker,values.get('genome_build',fallback_build),match.group(1))
    return list(rows.values())

def population_from_text(text):
    toks={x for x in re.split(r'[^A-Za-z0-9]+',str(text).upper()) if x}
    hit=sorted(toks & POPS_1KG)
    return hit[0] if len(hit)==1 else None

def reference_build_from_text(text, fallback='GRCh38'):
    s=str(text).lower()
    if any(x in s for x in ('chm13','t2t')): return 'CHM13'
    if any(x in s for x in ('grch38','hg38','b38')): return 'GRCh38'
    if any(x in s for x in ('grch37','hg19','b37')): return 'GRCh37'
    return fallback

def parse_reference_handle(handle, raw_name, dataset_id='AS3_1KG', fallback_build='GRCh38', chunksize=100000):
    """Parse BED-like population callsets without treating them as sample calls."""
    path_pop=population_from_text(raw_name)
    build=reference_build_from_text(raw_name,fallback_build)
    header=None; rows=[]
    for raw in handle:
        if isinstance(raw,bytes): raw=raw.decode('utf-8','replace')
        line=raw.strip()
        if not line or line.startswith(('#','track ','browser ')): continue
        fields=line.split('\t') if '\t' in line else line.split()
        if len(fields)<3: continue
        if header is None and not (re.fullmatch(r'-?\d+',fields[1]) and re.fullmatch(r'-?\d+',fields[2])):
            header={re.sub(r'[^a-z0-9]+','_',x.strip().lower()).strip('_'):i for i,x in enumerate(fields)}
            continue
        if header:
            ci=next((header[x] for x in ('chr','chrom','chromosome') if x in header),0)
            si=next((header[x] for x in ('start','start0','chromstart') if x in header),1)
            ei=next((header[x] for x in ('end','stop','chromend') if x in header),2)
            pi=next((header[x] for x in ('population','pop','population_id') if x in header),None)
            spi=next((header[x] for x in ('super_population','super_pop','superpopulation') if x in header),None)
            idi=next((header[x] for x in ('individual_id','sample_id','sample','individual','iid') if x in header),None)
        else: ci,si,ei,pi,spi,idi=0,1,2,None,None,None
        try: chrom=chrom_norm(fields[ci]); start=int(float(fields[si])); end=int(float(fields[ei]))
        except (ValueError,IndexError): continue
        if end<=start: continue
        pop=path_pop
        if pi is not None and pi<len(fields) and fields[pi].upper() in POPS_1KG: pop=fields[pi].upper()
        if pop is None:
            for value in fields[3:]:
                if value.upper() in POPS_1KG: pop=value.upper(); break
        super_pop=fields[spi].upper() if spi is not None and spi<len(fields) else None
        reference_sample=fields[idi] if idi is not None and idi<len(fields) else None
        # Official database_introgression_call headerless order is:
        # chr,start,end,LOD,length,population,super_population,individual_id.
        if header is None and len(fields)>=8 and fields[5].upper() in POPS_1KG:
            pop=fields[5].upper(); super_pop=fields[6].upper(); reference_sample=fields[7]
        rows.append((dataset_id,pop or 'ALL',build,chrom,start,end,'Neanderthal','external_reference',raw_name,reference_sample,super_pop))
        if len(rows)>=chunksize:
            yield pd.DataFrame(rows,columns=REFERENCE_RAW_COLS); rows=[]
    if rows: yield pd.DataFrame(rows,columns=REFERENCE_RAW_COLS)

def iter_reference_callsets(path,dataset_id='AS3_1KG'):
    if not path:return
    path=Path(path)
    if not path.is_dir():
        raise ValueError(f'reference callset directory is missing: {path}')
    files=sorted(x for x in path.rglob('*') if x.is_file() and re.search(r'\.(bed|tsv)(\.gz)?$',x.name,re.I) and x.stat().st_size>0)
    for f in files:
        opener=gzip.open if f.name.lower().endswith('.gz') else open
        with opener(f,'rt',encoding='utf-8',errors='replace') as h:
            yield from parse_reference_handle(h,str(f),dataset_id)

def sample_population_table(path):
    cols=['dataset_id','sample_id','population','super_population','raw_file']
    if not path or not Path(path).is_file():return pd.DataFrame(columns=cols)
    try:d=pd.read_csv(path,sep=r'\s+',dtype=str,comment='#')
    except Exception:return pd.DataFrame(columns=cols)
    low={str(c).lower().lstrip('#'):c for c in d.columns}
    sample=next((low[x] for x in ('sample','sample_id','iid') if x in low),None)
    pop=next((low[x] for x in ('population','pop') if x in low),None)
    super_pop=next((low[x] for x in ('super_population','super_pop','superpopulation') if x in low),None)
    if sample is None or pop is None:return pd.DataFrame(columns=cols)
    out=pd.DataFrame({'dataset_id':'1kg','sample_id':d[sample].astype(str),'population':d[pop].astype(str).str.upper()})
    out['super_population']=d[super_pop].astype(str).str.upper() if super_pop else pd.NA
    out['raw_file']=str(path)
    return out[out.sample_id.ne('') & out.population.isin(POPS_1KG)].drop_duplicates('sample_id')[cols]

def truth_int(value):
    if pd.isna(value):return 0
    return int(str(value).strip().upper() in {'1','TRUE','T','YES','Y','PASS','DIRECT_MATCH_PASS'})

def iter_phyml_carriers(root,fallback_build):
    """Normalize every tested PhyML copy and keep tree and pairwise flags separate."""
    for hp in sorted((root/'phyml').glob('**/final/haplotypes.tsv')):
        if not hp.is_file() or not hp.stat().st_size:continue
        try:
            haps=pd.read_csv(hp,sep='\t',dtype=str)
        except Exception:continue
        if not {'locus_id','hap_id'}.issubset(haps.columns):continue
        dataset=dataset_from_path(hp,'phyml')
        # The copies column records every tested sample/haplotype and is the
        # prevalence denominator.  A header-only table is a valid result for a
        # locus with no selected haplotypes (for example, low LD information).
        if 'copies' not in haps:
            raise ValueError(f'PhyML haplotypes table lacks required copies column: {hp}')
        if haps.empty:continue
        if not haps['copies'].notna().any():
            raise ValueError(f'PhyML haplotypes table has no populated copies values: {hp}')
        tree_tips={}; tree_lineages={}; tree_pass={}
        trees=hp.parent/'trees.tsv'
        if trees.is_file() and trees.stat().st_size:
            try:
                for _,row in pd.read_csv(trees,sep='\t',dtype=str,keep_default_na=False).iterrows():
                    lid=str(row.get('locus_id',''))
                    tips={x for x in str(row.get('candidate_clade_tips','')).split(',') if x.startswith('H')}
                    tree_tips[lid]=tips
                    tree_lineages[lid]=str(row.get('candidate_lineage',''))
                    tree_pass[lid]=truth_int(row.get('candidate_clade_pass',pd.NA))
            except Exception:
                tree_tips={};tree_lineages={};tree_pass={}
        buffer=[]
        for _,hap in haps.iterrows():
            lid=str(hap['locus_id']); hap_id=str(hap['hap_id'])
            in_supported_clade=int(bool(tree_pass.get(lid,0) and hap_id in tree_tips.get(lid,set())))
            copies='' if pd.isna(hap.get('copies',pd.NA)) else str(hap.get('copies',''))
            for copy in copies.split(';'):
                if not copy or ':' not in copy:continue
                sample,haplotype=copy.rsplit(':',1)
                if not sample:continue
                buffer.append({
                    'dataset_id':dataset,
                    'genome_build':build_norm(hap.get('genome_build',None),fallback_build),
                    'locus_id':lid,'sample_id':sample,'haplotype':haplotype,
                    'hap_id':hap_id,'best_archaic':hap.get('best_archaic',pd.NA),
                    'best_lineage':hap.get('best_lineage',pd.NA),
                    'candidate_lineage':tree_lineages.get(lid,''),
                    'n_compared':pd.to_numeric(hap.get('n_compared',pd.NA),errors='coerce'),
                    'n_match':pd.to_numeric(hap.get('n_match',pd.NA),errors='coerce'),
                    'prop_match':pd.to_numeric(hap.get('prop_match',pd.NA),errors='coerce'),
                    'direct_match_pass':truth_int(hap.get('direct_match_pass',pd.NA)),
                    'tree_candidate_pass':in_supported_clade,
                    'raw_file':str(hp),
                })
                if len(buffer)>=100000:
                    yield pd.DataFrame(buffer);buffer=[]
        if buffer:yield pd.DataFrame(buffer)

def loci_tables(root,fallback_build):
    def scalar_int(value):
        try:
            if value is None or pd.isna(value):return None
            return int(float(value))
        except (TypeError,ValueError):return None
    def scalar_float(value):
        try:
            if value is None or pd.isna(value):return None
            return float(value)
        except (TypeError,ValueError):return None
    def scalar_text(value):
        return None if value is None or pd.isna(value) or not str(value).strip() else str(value)
    seen=set()
    for base in (root/'phyml',):
      for pat in ['**/final/loci.tsv']:
        for f in base.glob(pat):
            if not f.is_file() or f.stat().st_size==0:continue
            try:d=pd.read_csv(f,sep='\t')
            except Exception:continue
            if d.empty:continue
            trees={}
            tree_file=f.parent/'trees.tsv'
            if tree_file.is_file() and tree_file.stat().st_size:
                try:
                    td=pd.read_csv(tree_file,sep='\t',dtype=str,keep_default_na=False)
                    trees={str(row.get('locus_id','')):row for _,row in td.iterrows()}
                except Exception:pass
            dataset=dataset_from_path(f,'phyml')
            for _,r in d.iterrows():
                chrom=r.get('lead_chr',r.get('CHR',r.get('chr','')))
                input_st=scalar_int(r.get('core_start',r.get('region_start',r.get('start',r.get('POS',None)))))
                input_en=scalar_int(r.get('core_end',r.get('region_end',r.get('end',r.get('POS',None)))))
                if input_st is None or input_en is None:continue
                if input_en<=input_st:input_en=input_st+1
                ast=scalar_int(r.get('analysis_start',r.get('region_start',input_st)))
                aen=scalar_int(r.get('analysis_end',r.get('region_end',input_en)))
                if ast is None or aen is None or aen<=ast:ast,aen=input_st,input_en
                st=scalar_int(r.get('selected_start',input_st))
                en=scalar_int(r.get('selected_end',input_en))
                if st is None or en is None or en<=st:st,en=input_st,input_en
                method='phyml'; lid=r.get('locus_id',r.get('name',None))
                genome_build=build_norm(r.get('genome_build',None),fallback_build)
                key=(dataset,method,genome_build,str(lid),chrom,st,en)
                if key in seen:continue
                seen.add(key)
                tree=trees.get(str(lid),{}); run_dir=f.parent.parent
                candidate_lineage=scalar_text(tree.get('candidate_lineage',None))
                source=candidate_lineage or next((v for v in (r.get('best_lineage'),r.get('matched_archaics'),r.get('asnp_lineage')) if not pd.isna(v) and str(v).strip()),'archaic_map')
                yield {
                    'dataset_id':dataset,'method':method,'trait':None,'locus_id':lid,'chr':chrom_norm(chrom),
                    'start':st,'end':en,'input_start':input_st,'input_end':input_en,
                    'analysis_start':ast,'analysis_end':aen,
                    'flank_bp':max(0,max(input_st-ast,aen-input_en)),
                    'anchor_pos':scalar_int(r.get('anchor_pos',None)),
                    'selection_method':scalar_text(r.get('selection_method',None)),
                    'ld_r2_threshold':scalar_float(r.get('ld_r2_threshold',None)),
                    'n_search_sites':scalar_int(r.get('n_search_sites',None)),
                    'n_ld_sites':scalar_int(r.get('n_ld_sites',None)),
                    'n_ancestral_sites':scalar_int(r.get('n_ancestral_sites',None)),
                    'n_sites':scalar_int(r.get('n_sites',None)),
                    'best_archaic':scalar_text(r.get('best_archaic',None)),
                    'best_lineage':scalar_text(r.get('best_lineage',None)),
                    'n_compared':scalar_int(r.get('n_compared',None)),
                    'n_match':scalar_int(r.get('n_match',None)),
                    'prop_match':scalar_float(r.get('prop_match',None)),
                    'source':str(source),'source_class':source_class(source),'status':r.get('status',None),
                    'direct_match_pass':r.get('direct_match_pass',None),
                    'n_carriers':None,'genome_build':genome_build,
                    'tree_status':tree.get('tree_status','not_run'),
                    'tree_file':phyml_artifact_in_run(tree.get('tree_file',''),run_dir),
                    'stats_file':phyml_artifact_in_run(tree.get('stats_file',''),run_dir),
                    'plot_file':phyml_artifact_in_run(tree.get('plot_file',''),run_dir),
                    'n_bootstrap_nodes':tree.get('n_bootstrap_nodes',None),
                    'bootstrap_min':tree.get('bootstrap_min',None),
                    'bootstrap_median':tree.get('bootstrap_median',None),
                    'bootstrap_max':tree.get('bootstrap_max',None),
                    'candidate_lineage':candidate_lineage,
                    'tree_has_ancestral_outgroup':truth_int(tree.get('tree_has_ancestral_outgroup',None)),
                    'candidate_clade_pass':truth_int(tree.get('candidate_clade_pass',None)),
                    'candidate_clade_rule':scalar_text(tree.get('candidate_clade_rule',None)),
                    'candidate_clade_bootstrap':scalar_float(tree.get('candidate_clade_bootstrap',None)),
                    'candidate_clade_node':scalar_text(tree.get('candidate_clade_node',None)),
                    'candidate_clade_side':scalar_text(tree.get('candidate_clade_side',None)),
                    'candidate_clade_n_tips':scalar_int(tree.get('candidate_clade_n_tips',None)),
                    'candidate_clade_modern_tips':scalar_int(tree.get('candidate_clade_modern_tips',None)),
                    'candidate_clade_archaic_tips':scalar_int(tree.get('candidate_clade_archaic_tips',None)),
                    'candidate_clade_specificity':scalar_float(tree.get('candidate_clade_specificity',None)),
                    'candidate_clade_tips':scalar_text(tree.get('candidate_clade_tips',None)),
                    'tree_newick':tree.get('tree_newick',''),'raw_file':str(f)
                }

def init_db(con):
    con.executescript('''
    PRAGMA journal_mode=DELETE; PRAGMA synchronous=FULL; PRAGMA temp_store=FILE;
    DROP TABLE IF EXISTS segments_raw; DROP TABLE IF EXISTS loci; DROP TABLE IF EXISTS coord_catalog;
    DROP TABLE IF EXISTS segment_catalog; DROP TABLE IF EXISTS segments; DROP TABLE IF EXISTS carriers;
    DROP TABLE IF EXISTS reference_callsets; DROP TABLE IF EXISTS reference_callset_overlaps; DROP TABLE IF EXISTS database_summary;
    DROP TABLE IF EXISTS phyml_haplotype_carriers; DROP TABLE IF EXISTS locus_sample_support;
    DROP TABLE IF EXISTS locus_method_support; DROP TABLE IF EXISTS locus_trajectory; DROP TABLE IF EXISTS locus_evidence;
    DROP TABLE IF EXISTS method_runs;
    DROP TABLE IF EXISTS sample_populations; DROP TABLE IF EXISTS sample_burden; DROP TABLE IF EXISTS burden_union_chr;
    DROP TABLE IF EXISTS temp.reference_callsets_raw;
    CREATE TABLE segments_raw(dataset_id TEXT,sample_id TEXT,method TEXT,source TEXT,source_class TEXT,chr TEXT,start INTEGER,end INTEGER,length_bp INTEGER,haplotype TEXT,method_haplotype_index TEXT,score REAL,posterior REAL,trait TEXT,locus_id TEXT,genome_build TEXT,batch_id TEXT,raw_file TEXT);
    CREATE TABLE loci(dataset_id TEXT,method TEXT,trait TEXT,locus_id TEXT,chr TEXT,start INTEGER,end INTEGER,input_start INTEGER,input_end INTEGER,analysis_start INTEGER,analysis_end INTEGER,flank_bp INTEGER,anchor_pos INTEGER,selection_method TEXT,ld_r2_threshold REAL,n_search_sites INTEGER,n_ld_sites INTEGER,n_ancestral_sites INTEGER,n_sites INTEGER,best_archaic TEXT,best_lineage TEXT,n_compared INTEGER,n_match INTEGER,prop_match REAL,source TEXT,source_class TEXT,status TEXT,direct_match_pass INTEGER,n_carriers INTEGER,genome_build TEXT,tree_status TEXT,tree_file TEXT,stats_file TEXT,plot_file TEXT,n_bootstrap_nodes INTEGER,bootstrap_min REAL,bootstrap_median REAL,bootstrap_max REAL,candidate_lineage TEXT,tree_has_ancestral_outgroup INTEGER,candidate_clade_pass INTEGER,candidate_clade_rule TEXT,candidate_clade_bootstrap REAL,candidate_clade_node TEXT,candidate_clade_side TEXT,candidate_clade_n_tips INTEGER,candidate_clade_modern_tips INTEGER,candidate_clade_archaic_tips INTEGER,candidate_clade_specificity REAL,candidate_clade_tips TEXT,tree_newick TEXT,raw_file TEXT);
    CREATE TABLE coord_catalog(dataset_id TEXT,genome_build TEXT,source_class TEXT,chr TEXT,start INTEGER,end INTEGER,segment_code TEXT);
    CREATE TABLE segment_catalog(segment_code TEXT PRIMARY KEY,dataset_id TEXT,source_class TEXT,chr TEXT,start INTEGER,end INTEGER,length_bp INTEGER,n_distinct_boundaries INTEGER,genome_build TEXT);
    CREATE TABLE reference_callsets(dataset_id TEXT,population TEXT,genome_build TEXT,chr TEXT,start INTEGER,end INTEGER,source_class TEXT,reference_role TEXT,raw_file TEXT);
    CREATE TEMP TABLE reference_callsets_raw(dataset_id TEXT,population TEXT,genome_build TEXT,chr TEXT,start INTEGER,end INTEGER,source_class TEXT,reference_role TEXT,raw_file TEXT,reference_sample_id TEXT,reference_super_population TEXT);
    CREATE TABLE sample_populations(dataset_id TEXT,sample_id TEXT,population TEXT,super_population TEXT,raw_file TEXT,PRIMARY KEY(dataset_id,sample_id));
    CREATE TABLE phyml_haplotype_carriers(dataset_id TEXT,genome_build TEXT,locus_id TEXT,sample_id TEXT,haplotype TEXT,hap_id TEXT,best_archaic TEXT,best_lineage TEXT,candidate_lineage TEXT,n_compared INTEGER,n_match INTEGER,prop_match REAL,direct_match_pass INTEGER,tree_candidate_pass INTEGER,raw_file TEXT);
    CREATE TABLE method_runs(dataset_id TEXT,genome_build TEXT,chr TEXT,method TEXT,status TEXT,evidence_eligible INTEGER,availability_note TEXT,raw_file TEXT);
    CREATE TABLE database_summary(dataset_id TEXT,genome_build TEXT,metric TEXT,value INTEGER,PRIMARY KEY(dataset_id,genome_build,metric));
    ''')

def insert_df(con,table,df):
    if df is None or df.empty:return 0
    df.to_sql(table,con,if_exists='append',index=False,method='multi',chunksize=1000); return len(df)

def insert_reference_df(con,df):
    if df is None or df.empty:return 0
    con.executemany(
        'INSERT INTO temp.reference_callsets_raw VALUES(?,?,?,?,?,?,?,?,?,?,?)',
        df[REFERENCE_RAW_COLS].itertuples(index=False,name=None),
    )
    return len(df)

def reciprocal(a,b,c,d):
    ov=max(0,min(b,d)-max(a,c));
    if ov<=0:return 0.0
    return min(ov/max(1,b-a),ov/max(1,d-c))

def code_for(dataset,build,src,chrom,start,end):
    tag={'Neanderthal':'NEA','Denisovan':'DEN','Mosaic':'MOS','Ghost/Unknown':'GST','Archaic/Other':'ARC'}.get(src,'ARC')
    dig=hashlib.sha1(f'{dataset}|{build}|{src}|{chrom}|{round(start,-3)}|{round(end,-3)}'.encode()).hexdigest()[:10].upper()
    return f'GU{build.replace("GRCh","")}-{tag}-{chrom}-{dig}'

def build_catalog(con,build,threshold):
    groups=con.execute('SELECT DISTINCT dataset_id,genome_build,source_class,chr FROM segments_raw ORDER BY dataset_id,genome_build,source_class, CASE WHEN chr GLOB "[0-9]*" THEN CAST(chr AS INT) ELSE 99 END, chr').fetchall()
    ins_coord=[]; ins_cat=[]
    for dataset,group_build,src,chrom in groups:
        cur=con.execute('SELECT DISTINCT start,end FROM segments_raw WHERE dataset_id=? AND genome_build=? AND source_class=? AND chr=? ORDER BY start,end',(dataset,group_build,src,chrom))
        coords=cur.fetchall()
        cluster=[]; mean_s=mean_e=None
        def flush():
            nonlocal cluster,mean_s,mean_e,ins_coord,ins_cat
            if not cluster:return
            rs=int(round(mean_s)); re=int(round(mean_e)); code=code_for(dataset,group_build,src,chrom,rs,re)
            ins_cat.append((code,dataset,src,chrom,rs,re,re-rs,len(cluster),group_build))
            ins_coord.extend((dataset,group_build,src,chrom,s,e,code) for s,e in cluster)
            if len(ins_coord)>=100000:
                con.executemany('INSERT INTO coord_catalog VALUES(?,?,?,?,?,?,?)',ins_coord); ins_coord=[]
            if len(ins_cat)>=10000:
                con.executemany('INSERT OR REPLACE INTO segment_catalog VALUES(?,?,?,?,?,?,?,?,?)',ins_cat); ins_cat=[]
            cluster=[]; mean_s=mean_e=None
        for s,e in coords:
            if not cluster:
                cluster=[(s,e)]; mean_s=float(s); mean_e=float(e); continue
            if reciprocal(s,e,mean_s,mean_e) >= threshold:
                n=len(cluster); cluster.append((s,e)); mean_s=(mean_s*n+s)/(n+1); mean_e=(mean_e*n+e)/(n+1)
            else:
                flush(); cluster=[(s,e)]; mean_s=float(s); mean_e=float(e)
        flush()
    if ins_coord: con.executemany('INSERT INTO coord_catalog VALUES(?,?,?,?,?,?,?)',ins_coord)
    if ins_cat: con.executemany('INSERT OR REPLACE INTO segment_catalog VALUES(?,?,?,?,?,?,?,?,?)',ins_cat)
    con.commit()
    con.executescript('''
      CREATE INDEX idx_coord ON coord_catalog(dataset_id,genome_build,source_class,chr,start,end);
      CREATE TABLE segments AS SELECT s.*, c.segment_code FROM segments_raw s JOIN coord_catalog c USING(dataset_id,genome_build,source_class,chr,start,end);
      CREATE INDEX idx_segments_sample ON segments(dataset_id,sample_id); CREATE INDEX idx_segments_region ON segments(dataset_id,genome_build,chr,start,end); CREATE INDEX idx_segments_code ON segments(segment_code);
      CREATE INDEX idx_segments_method_source ON segments(dataset_id,genome_build,method,source_class,chr,start);
      CREATE INDEX idx_segments_build_region ON segments(genome_build,chr,start,end);
      CREATE INDEX idx_segments_build_method ON segments(genome_build,dataset_id,method,source_class,chr,start);
      CREATE TABLE carriers AS
        SELECT dataset_id,sample_id,segment_code,source_class,
               MIN(2, CASE WHEN COUNT(DISTINCT CASE WHEN haplotype NOT IN ('','nan','None') THEN haplotype END)>0 THEN COUNT(DISTINCT CASE WHEN haplotype NOT IN ('','nan','None') THEN haplotype END) ELSE 1 END) AS dosage,
               COUNT(DISTINCT method) AS n_methods,
               GROUP_CONCAT(DISTINCT method) AS methods_support,
               MAX(score) AS max_score, MAX(posterior) AS max_posterior
        FROM segments WHERE sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None') GROUP BY dataset_id,sample_id,segment_code,source_class;
      CREATE INDEX idx_carriers_sample ON carriers(dataset_id,sample_id); CREATE INDEX idx_carriers_code ON carriers(segment_code);
    ''')
    con.commit()

def build_reference_overlaps(con,raw_table='temp.reference_callsets_raw'):
    if raw_table not in ('temp.reference_callsets_raw','refcache.reference_callsets_raw'):
        raise ValueError(f'unsupported reference raw table: {raw_table}')
    con.executescript('''
      CREATE INDEX IF NOT EXISTS idx_reference_callsets_region ON reference_callsets(genome_build,source_class,chr,start,end);
      CREATE INDEX IF NOT EXISTS idx_reference_callsets_population ON reference_callsets(dataset_id,population,genome_build);
      CREATE INDEX IF NOT EXISTS idx_reference_callsets_build_population ON reference_callsets(genome_build,population);
      CREATE INDEX IF NOT EXISTS idx_reference_callsets_shiny_region ON reference_callsets(genome_build,chr,start,end,population);
      DROP TABLE IF EXISTS reference_callset_overlaps;
      CREATE TABLE reference_callset_overlaps(
        result_rowid INTEGER,result_dataset_id TEXT,sample_id TEXT,sample_population TEXT,method TEXT,source_class TEXT,segment_code TEXT,
        dataset_id TEXT,reference_population TEXT,reference_sample_id TEXT,genome_build TEXT,chr TEXT,
        segment_start INTEGER,segment_end INTEGER,reference_start INTEGER,reference_end INTEGER,
        overlap_bp INTEGER,result_overlap REAL,reference_overlap REAL,reciprocal_overlap REAL,
        reference_role TEXT,result_raw_file TEXT,reference_raw_file TEXT
      );
    ''')
    n_results=con.execute("SELECT COUNT(*) FROM segments WHERE method='as3' AND source_class='Neanderthal'").fetchone()[0]
    has_reference_samples=con.execute(f"SELECT COUNT(*) FROM {raw_table} WHERE reference_sample_id IS NOT NULL AND reference_sample_id!=''").fetchone()[0]
    if n_results and has_reference_samples:
        if raw_table.startswith('temp.'):
            con.execute('CREATE INDEX temp.idx_reference_raw_sample_region ON reference_callsets_raw(reference_sample_id,genome_build,source_class,chr,start,end)')
        con.execute(f'''
          INSERT INTO reference_callset_overlaps
          SELECT s.rowid,s.dataset_id,s.sample_id,COALESCE(p.population,r.population),s.method,s.source_class,s.segment_code,
                 r.dataset_id,r.population,r.reference_sample_id,s.genome_build,s.chr,
                 s.start,s.end,r.start,r.end,
                 MIN(s.end,r.end)-MAX(s.start,r.start) AS overlap_bp,
                 CAST(MIN(s.end,r.end)-MAX(s.start,r.start) AS REAL)/MAX(1,s.end-s.start) AS result_overlap,
                 CAST(MIN(s.end,r.end)-MAX(s.start,r.start) AS REAL)/MAX(1,r.end-r.start) AS reference_overlap,
                 MIN(
                   CAST(MIN(s.end,r.end)-MAX(s.start,r.start) AS REAL)/MAX(1,s.end-s.start),
                   CAST(MIN(s.end,r.end)-MAX(s.start,r.start) AS REAL)/MAX(1,r.end-r.start)
                 ) AS reciprocal_overlap,
                 r.reference_role,s.raw_file,r.raw_file
          FROM segments s
          JOIN {raw_table} r
            ON r.reference_sample_id=s.sample_id
           AND r.genome_build=s.genome_build AND r.source_class=s.source_class AND r.chr=s.chr
           AND r.end>s.start AND r.start<s.end
          LEFT JOIN sample_populations p ON p.dataset_id=s.dataset_id AND p.sample_id=s.sample_id
          WHERE s.method='as3' AND s.source_class='Neanderthal'
            AND (p.population IS NULL OR r.population='ALL' OR r.population=p.population)
        ''')
    elif n_results and con.execute('SELECT COUNT(*) FROM reference_callsets').fetchone()[0]:
        # Generic population callsets without individual IDs still support a
        # population-matched comparison against the unioned overlay table.
        con.execute('''
          INSERT INTO reference_callset_overlaps
          SELECT s.rowid,s.dataset_id,s.sample_id,p.population,s.method,s.source_class,s.segment_code,
                 r.dataset_id,r.population,NULL,s.genome_build,s.chr,
                 s.start,s.end,r.start,r.end,
                 MIN(s.end,r.end)-MAX(s.start,r.start),
                 CAST(MIN(s.end,r.end)-MAX(s.start,r.start) AS REAL)/MAX(1,s.end-s.start),
                 CAST(MIN(s.end,r.end)-MAX(s.start,r.start) AS REAL)/MAX(1,r.end-r.start),
                 MIN(
                   CAST(MIN(s.end,r.end)-MAX(s.start,r.start) AS REAL)/MAX(1,s.end-s.start),
                   CAST(MIN(s.end,r.end)-MAX(s.start,r.start) AS REAL)/MAX(1,r.end-r.start)
                 ),
                 r.reference_role,s.raw_file,r.raw_file
          FROM segments s
          JOIN reference_callsets r
            ON r.genome_build=s.genome_build AND r.source_class=s.source_class AND r.chr=s.chr
           AND r.end>s.start AND r.start<s.end
          LEFT JOIN sample_populations p ON p.dataset_id=s.dataset_id AND p.sample_id=s.sample_id
          WHERE s.method='as3' AND s.source_class='Neanderthal'
            AND (r.population='ALL' OR r.population=p.population)
        ''')
    con.executescript('''
      CREATE INDEX idx_reference_overlaps_region ON reference_callset_overlaps(genome_build,chr,segment_start,segment_end);
      CREATE INDEX idx_reference_overlaps_population ON reference_callset_overlaps(dataset_id,reference_population,method);
      CREATE INDEX idx_reference_overlaps_sample ON reference_callset_overlaps(result_dataset_id,sample_id);
    ''')
    if raw_table.startswith('temp.'):con.execute('DELETE FROM temp.reference_callsets_raw')
    con.commit()

def collapse_reference_callsets(con):
    """Union carrier-level official calls into population callset intervals.

    The published BED contains individual IDs, while the requested external
    reference schema is population-level and intentionally has no individual
    carrier column.  Collapsing overlapping/touching rows here avoids storing
    millions of anonymous duplicate carrier intervals or multiplying every
    GU-AS3 overlap by the number of reference individuals.
    """
    con.executescript('''
      DELETE FROM reference_callsets;
      INSERT INTO reference_callsets
      WITH preceding AS (
        SELECT *,
          MAX(end) OVER (
            PARTITION BY dataset_id,population,genome_build,chr,source_class,reference_role
            ORDER BY start,end ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
          ) AS preceding_max_end
        FROM reference_callsets_raw
      ), marked AS (
        SELECT *,CASE WHEN preceding_max_end IS NULL OR start>preceding_max_end THEN 1 ELSE 0 END AS new_island
        FROM preceding
      ), islands AS (
        SELECT *,
          SUM(new_island) OVER (
            PARTITION BY dataset_id,population,genome_build,chr,source_class,reference_role
            ORDER BY start,end ROWS UNBOUNDED PRECEDING
          ) AS island_id
        FROM marked
      )
      SELECT dataset_id,population,genome_build,chr,MIN(start),MAX(end),source_class,reference_role,MIN(raw_file)
      FROM islands
      GROUP BY dataset_id,population,genome_build,chr,source_class,reference_role,island_id;
    ''')
    con.commit()

def interval_summary(intervals,weighted=False):
    if not intervals:return 0,0,0,0
    x=sorted(intervals,key=lambda z:(z[0],z[1]))
    merged=0; total=0; s=e=None
    for start,end,*_ in x:
        if s is None:s,e=start,end
        elif start<=e:e=max(e,end)
        else:merged+=1; total+=e-s; s,e=start,end
    merged+=1; total+=e-s
    if not weighted:return len(x),merged,total,total
    events={}
    for start,end,dose in x:
        dose=max(1,min(2,int(dose or 1)))
        events.setdefault(start,[]).append((1,dose)); events.setdefault(end,[]).append((-1,dose))
    active={1:0,2:0}; prev=None; dosage_bp=0
    for pos in sorted(events):
        if prev is not None:
            dose=2 if active[2] else (1 if active[1] else 0)
            dosage_bp+=(pos-prev)*dose
        for delta,dose in events[pos]:active[dose]+=delta
        prev=pos
    return len(x),merged,total,dosage_bp

def insert_interval_burden(con,query,burden_type,weighted=False):
    cur=con.execute(query); current=None; intervals=[]; rows=[]
    def flush():
        nonlocal intervals,rows
        if current is None:return
        n,merged,total,dosage=interval_summary(intervals,weighted)
        rows.append((*current,burden_type,n,merged,total,dosage))
        intervals=[]
        if len(rows)>=10000:
            con.executemany('INSERT INTO burden_union_chr VALUES(?,?,?,?,?,?,?,?,?,?,?)',rows); rows=[]
    for row in cur:
        key=tuple(row[:6]); interval=tuple(row[6:])
        if current is not None and key!=current:flush()
        current=key; intervals.append(interval)
    flush()
    if rows:con.executemany('INSERT INTO burden_union_chr VALUES(?,?,?,?,?,?,?,?,?,?,?)',rows)

def build_burden(con):
    con.executescript('''
      DROP TABLE IF EXISTS sample_burden;
      CREATE TABLE sample_burden(
        dataset_id TEXT,sample_id TEXT,method TEXT,source_class TEXT,genome_build TEXT,burden_type TEXT,
        n_input_segments INTEGER,n_merged_intervals INTEGER,n_chromosomes INTEGER,total_bp INTEGER,dosage_bp INTEGER
      );
      INSERT INTO sample_burden
        SELECT dataset_id,sample_id,method,source_class,genome_build,'raw_call',COUNT(*),NULL,COUNT(DISTINCT chr),SUM(length_bp),SUM(length_bp)
        FROM segments WHERE sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None')
        GROUP BY dataset_id,sample_id,method,source_class,genome_build;
      CREATE TEMP TABLE burden_union_chr(
        dataset_id TEXT,sample_id TEXT,method TEXT,source_class TEXT,genome_build TEXT,chr TEXT,burden_type TEXT,
        n_input_segments INTEGER,n_merged_intervals INTEGER,total_bp INTEGER,dosage_bp INTEGER
      );
    ''')
    # The biological source_class collapses Altai/Chagyrskaya/Vindija support,
    # preventing the same Neanderthal tract from being counted repeatedly.
    insert_interval_burden(con,'''
      SELECT dataset_id,sample_id,method,source_class,genome_build,chr,start,end
      FROM segments WHERE sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None')
      ORDER BY dataset_id,sample_id,method,source_class,genome_build,chr,start,end
    ''','nonredundant_union',False)
    insert_interval_burden(con,'''
      SELECT c.dataset_id,c.sample_id,'consensus',c.source_class,s.genome_build,s.chr,s.start,s.end,c.dosage
      FROM carriers c JOIN segment_catalog s
        ON s.dataset_id=c.dataset_id AND s.segment_code=c.segment_code AND s.source_class=c.source_class
      WHERE c.n_methods>=2
      ORDER BY c.dataset_id,c.sample_id,c.source_class,s.genome_build,s.chr,s.start,s.end
    ''','consensus_catalog',True)
    con.execute('''
      INSERT INTO sample_burden
      SELECT dataset_id,sample_id,method,source_class,genome_build,burden_type,
             SUM(n_input_segments),SUM(n_merged_intervals),COUNT(DISTINCT chr),SUM(total_bp),SUM(dosage_bp)
      FROM burden_union_chr GROUP BY dataset_id,sample_id,method,source_class,genome_build,burden_type
    ''')
    con.executescript('''
      CREATE INDEX idx_sample_burden_sample ON sample_burden(dataset_id,sample_id);
      CREATE INDEX idx_sample_burden_filter ON sample_burden(genome_build,dataset_id,burden_type,method,source_class);
    ''')
    con.commit()

def wilson_interval(successes,total,z=1.959963984540054):
    if not total:return (None,None,None)
    p=successes/total; z2=z*z; denominator=1+z2/total
    centre=(p+z2/(2*total))/denominator
    half=z*math.sqrt((p*(1-p)+z2/(4*total))/total)/denominator
    return p,max(0.0,centre-half),min(1.0,centre+half)

def build_locus_support(con,trajectory_bins=120):
    """Align PhyML carriers with individual tract calls and build scalable curves."""
    con.executescript('''
      CREATE INDEX idx_phyml_carrier_locus ON phyml_haplotype_carriers(dataset_id,genome_build,locus_id,sample_id);
      CREATE INDEX idx_phyml_carrier_match ON phyml_haplotype_carriers(dataset_id,genome_build,locus_id,tree_candidate_pass);
      DROP TABLE IF EXISTS locus_sample_support;
      CREATE TABLE locus_sample_support(
        dataset_id TEXT,genome_build TEXT,locus_id TEXT,chr TEXT,source_class TEXT,
        sample_id TEXT,population TEXT,super_population TEXT,
        phyml INTEGER,ibdmix INTEGER,trace INTEGER,as3 INTEGER,n_methods INTEGER,
        methods_support TEXT,matched_haplotypes INTEGER,best_lineages TEXT,matched_source_classes TEXT,max_prop_match REAL,
        max_segment_score REAL,max_posterior REAL
      );
      WITH tested AS (
        SELECT pc.dataset_id,pc.genome_build,pc.locus_id,l.chr,l.source_class,pc.sample_id,
               p.population,p.super_population,
               MAX(pc.tree_candidate_pass) AS phyml,
               COUNT(DISTINCT CASE WHEN pc.tree_candidate_pass=1 THEN pc.hap_id END) AS matched_haplotypes,
               GROUP_CONCAT(DISTINCT CASE WHEN pc.tree_candidate_pass=1 THEN pc.candidate_lineage END) AS best_lineages,
               GROUP_CONCAT(DISTINCT CASE WHEN pc.tree_candidate_pass=1 THEN
                 CASE
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%neander%' OR LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%altai%' OR LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%vindija%' OR LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%chagyr%' THEN 'Neanderthal'
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%denis%' THEN 'Denisovan'
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%mosaic%' THEN 'Mosaic'
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%ghost%' OR LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%unknown%' THEN 'Ghost/Unknown'
                   ELSE 'Archaic/Other'
                 END END) AS matched_source_classes,
               MAX(CASE WHEN pc.tree_candidate_pass=1 THEN pc.prop_match END) AS max_prop_match
        FROM phyml_haplotype_carriers pc
        JOIN loci l ON l.dataset_id=pc.dataset_id AND l.genome_build=pc.genome_build
                   AND l.locus_id=pc.locus_id AND l.method='phyml'
        LEFT JOIN sample_populations p ON p.dataset_id=pc.dataset_id AND p.sample_id=pc.sample_id
        GROUP BY pc.dataset_id,pc.genome_build,pc.locus_id,l.chr,l.source_class,pc.sample_id,p.population,p.super_population
      ), flags AS (
        SELECT t.dataset_id,t.genome_build,t.locus_id,t.sample_id,
               MAX(s.method='ibdmix') AS ibdmix,MAX(s.method='trace') AS trace,MAX(s.method='as3') AS as3,
               MAX(s.score) AS max_segment_score,MAX(s.posterior) AS max_posterior
        FROM tested t
        JOIN loci l ON l.dataset_id=t.dataset_id AND l.genome_build=t.genome_build AND l.locus_id=t.locus_id
        JOIN segments s ON s.dataset_id=t.dataset_id AND s.genome_build=t.genome_build AND s.sample_id=t.sample_id
                       AND s.chr=l.chr AND (
                         s.method='trace'
                         OR (t.phyml=1 AND INSTR(','||COALESCE(t.matched_source_classes,'')||',',','||s.source_class||',')>0)
                         OR (t.phyml=0 AND s.source_class=l.source_class)
                       )
                       AND s.end>l.start AND s.start<l.end
                       AND s.method IN ('ibdmix','trace','as3')
        GROUP BY t.dataset_id,t.genome_build,t.locus_id,t.sample_id
      )
      INSERT INTO locus_sample_support
      SELECT t.dataset_id,t.genome_build,t.locus_id,t.chr,t.source_class,t.sample_id,t.population,t.super_population,
             t.phyml,COALESCE(f.ibdmix,0),COALESCE(f.trace,0),COALESCE(f.as3,0),
             t.phyml+COALESCE(f.ibdmix,0)+COALESCE(f.trace,0)+COALESCE(f.as3,0),
             TRIM((CASE WHEN t.phyml=1 THEN 'phyml,' ELSE '' END)||
                  (CASE WHEN COALESCE(f.ibdmix,0)=1 THEN 'ibdmix,' ELSE '' END)||
                  (CASE WHEN COALESCE(f.trace,0)=1 THEN 'trace,' ELSE '' END)||
                  (CASE WHEN COALESCE(f.as3,0)=1 THEN 'as3,' ELSE '' END),','),
             t.matched_haplotypes,t.best_lineages,t.matched_source_classes,t.max_prop_match,f.max_segment_score,f.max_posterior
      FROM tested t LEFT JOIN flags f USING(dataset_id,genome_build,locus_id,sample_id);
      CREATE UNIQUE INDEX idx_locus_sample_support_pk ON locus_sample_support(dataset_id,genome_build,locus_id,sample_id);
      CREATE INDEX idx_locus_sample_support_filter ON locus_sample_support(dataset_id,genome_build,locus_id,population,n_methods);
      UPDATE loci SET n_carriers=(
        SELECT COUNT(*) FROM locus_sample_support s
        WHERE s.dataset_id=loci.dataset_id AND s.genome_build=loci.genome_build
          AND s.locus_id=loci.locus_id AND s.phyml=1
      ) WHERE method='phyml';
      DROP TABLE IF EXISTS locus_method_support;
      CREATE TABLE locus_method_support(
        dataset_id TEXT,genome_build TEXT,locus_id TEXT,chr TEXT,population TEXT,method TEXT,source_class TEXT,
        method_available INTEGER,evidence_eligible INTEGER,availability_note TEXT,
        n_tested INTEGER,n_carriers INTEGER,carrier_fraction REAL,ci_low REAL,ci_high REAL,
        n_phyml_carriers INTEGER,n_phyml_supported INTEGER,phyml_support_fraction REAL
      );
      DROP TABLE IF EXISTS locus_trajectory;
      CREATE TABLE locus_trajectory(
        dataset_id TEXT,genome_build TEXT,locus_id TEXT,chr TEXT,population TEXT,denominator_mode TEXT,method TEXT,source_class TEXT,
        bin_start INTEGER,bin_end INTEGER,bin_mid INTEGER,n_denominator INTEGER,n_carriers INTEGER,
        prevalence REAL,ci_low REAL,ci_high REAL
      );
    ''')

    available={(d,b,c,m) for d,b,c,m in con.execute(
        "SELECT DISTINCT dataset_id,genome_build,chr,method FROM segments WHERE method IN ('ibdmix','trace','as3') UNION SELECT dataset_id,genome_build,chr,method FROM method_runs"
    )}
    run_info={(d,b,c,m):(int(eligible or 0),note or '') for d,b,c,m,eligible,note in con.execute(
        "SELECT dataset_id,genome_build,chr,method,MAX(evidence_eligible),MAX(availability_note) FROM method_runs GROUP BY dataset_id,genome_build,chr,method"
    )}
    aggregate_sql='''
      SELECT dataset_id,genome_build,locus_id,chr,source_class,'ALL' AS population,COUNT(*) AS n_tested,
             SUM(phyml),SUM(ibdmix),SUM(trace),SUM(as3),
             SUM(phyml*ibdmix),SUM(phyml*trace),SUM(phyml*as3)
      FROM locus_sample_support GROUP BY dataset_id,genome_build,locus_id,chr,source_class
      UNION ALL
      SELECT dataset_id,genome_build,locus_id,chr,source_class,population,COUNT(*),
             SUM(phyml),SUM(ibdmix),SUM(trace),SUM(as3),
             SUM(phyml*ibdmix),SUM(phyml*trace),SUM(phyml*as3)
      FROM locus_sample_support WHERE population IS NOT NULL AND population!=''
      GROUP BY dataset_id,genome_build,locus_id,chr,source_class,population
    '''
    support_rows=[]
    for d,b,lid,chrom,src,pop,n,phy,ibd,trc,a3,phy_ibd,phy_trc,phy_a3 in con.execute(aggregate_sql):
        values={'phyml':(phy,phy),'ibdmix':(ibd,phy_ibd),'trace':(trc,phy_trc),'as3':(a3,phy_a3)}
        for method,(carriers,intersection) in values.items():
            carriers=int(carriers or 0); intersection=int(intersection or 0); n=int(n or 0); phy=int(phy or 0)
            fraction,lo,hi=wilson_interval(carriers,n)
            method_source='Ghost/Unknown' if method=='trace' else src
            key=(d,b,chrom,method); method_available=int(method=='phyml' or key in available)
            if method=='phyml':eligible,note=1,'bootstrap-supported tree clade; pairwise match is descriptive only'
            elif key in run_info:eligible,note=run_info[key]
            elif chrom=='X' and method=='ibdmix':eligible,note=0,'experimental male-X pseudo-diploid adaptation; excluded from GU evidence score'
            elif chrom=='X' and method=='as3':eligible,note=0,'AS3 Ref1028/model is autosomes-only; chrX is unsupported'
            elif method=='trace':eligible,note=1,'haplotype-node call; validity depends on the full-chromosome ARG backend and calibration'
            else:eligible,note=1,'segment result observed; completion metadata unavailable'
            support_rows.append((d,b,lid,chrom,pop,method,method_source,
                method_available,int(eligible),note,n,carriers,fraction,lo,hi,
                phy,intersection,(intersection/phy if phy else None)))
    con.executemany('INSERT INTO locus_method_support VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',support_rows)
    con.executescript('''
      CREATE UNIQUE INDEX idx_locus_method_support_pk ON locus_method_support(dataset_id,genome_build,locus_id,population,method);
      CREATE INDEX idx_locus_method_support_filter ON locus_method_support(dataset_id,genome_build,locus_id,population);
    ''')

    trajectory_rows=[]
    loci=con.execute('''
      SELECT dataset_id,genome_build,locus_id,chr,source_class,start,end,analysis_start,analysis_end
      FROM loci WHERE method='phyml' ORDER BY dataset_id,genome_build,locus_id
    ''').fetchall()
    for d,b,lid,chrom,src,core_start,core_end,analysis_start,analysis_end in loci:
        ast=int(analysis_start if analysis_start is not None else core_start)
        aen=int(analysis_end if analysis_end is not None else core_end)
        if aen<=ast:continue
        bins=max(1,min(int(trajectory_bins),aen-ast))
        edges=[ast+(aen-ast)*i//bins for i in range(bins+1)]
        denominators={'ALL':con.execute('''
          SELECT COUNT(*) FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=?
        ''',(d,b,lid)).fetchone()[0]}
        phy_counts={'ALL':con.execute('''
          SELECT COALESCE(SUM(phyml),0) FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=?
        ''',(d,b,lid)).fetchone()[0]}
        for pop,n,carriers in con.execute('''
          SELECT population,COUNT(*),SUM(phyml) FROM locus_sample_support
          WHERE dataset_id=? AND genome_build=? AND locus_id=? AND population IS NOT NULL AND population!=''
          GROUP BY population
        ''',(d,b,lid)):
            denominators[pop]=n; phy_counts[pop]=carriers or 0
        diffs={}
        def diff_for(method,pop,mode):return diffs.setdefault((method,pop,mode),[0]*(bins+1))
        current=None; ranges=[]; current_pop=''; current_phyml=0
        def flush_ranges():
            nonlocal ranges
            if current is None or not ranges:return
            method,_sample=current
            merged=[]
            for lo,hi in sorted(ranges):
                if merged and lo<=merged[-1][1]+1:merged[-1][1]=max(merged[-1][1],hi)
                else:merged.append([lo,hi])
            for pop in ('ALL',current_pop):
                if pop!='ALL' and (not pop or pop not in denominators):continue
                modes=['tested']+(['phyml_positive'] if current_phyml else [])
                for mode in modes:
                    diff=diff_for(method,pop,mode)
                    for lo,hi in merged:diff[lo]+=1;diff[hi+1]-=1
            ranges=[]
        query='''
          SELECT s.method,s.sample_id,COALESCE(ls.population,''),ls.phyml,s.start,s.end
          FROM segments s JOIN locus_sample_support ls
            ON ls.dataset_id=s.dataset_id AND ls.genome_build=s.genome_build AND ls.sample_id=s.sample_id
           AND ls.locus_id=?
          WHERE s.dataset_id=? AND s.genome_build=? AND s.chr=? AND (
              s.method='trace'
              OR (ls.phyml=1 AND INSTR(','||COALESCE(ls.matched_source_classes,'')||',',','||s.source_class||',')>0)
              OR (ls.phyml=0 AND s.source_class=?)
            )
            AND s.method IN ('ibdmix','trace','as3') AND s.end>? AND s.start<?
          ORDER BY s.method,s.sample_id,s.start,s.end
        '''
        for method,sample,pop,is_phyml,start,end in con.execute(query,(lid,d,b,chrom,src,ast,aen)):
            key=(method,sample)
            if current is not None and key!=current:flush_ranges()
            current=key;current_pop=pop or '';current_phyml=int(is_phyml or 0)
            start=max(ast,int(start));end=min(aen,int(end))
            if end<=start:continue
            lo=max(0,min(bins-1,bisect.bisect_right(edges,start)-1))
            hi=max(0,min(bins-1,bisect.bisect_left(edges,end)-1))
            if hi>=lo:ranges.append((lo,hi))
        flush_ranges()
        locus_methods={m for dd,bb,cc,m in available if (dd,bb,cc)==(d,b,chrom)}
        for pop,tested_total in denominators.items():
            phy_n=int(phy_counts.get(pop,0))
            for mode,total in (('tested',tested_total),('phyml_positive',phy_n)):
                if not total:continue
                for method in sorted(locus_methods):diff_for(method,pop,mode)
                for method in sorted(locus_methods):
                    running=0;diff=diff_for(method,pop,mode)
                    for i in range(bins):
                        running+=diff[i];fraction,lo,hi=wilson_interval(running,total)
                        method_source='Ghost/Unknown' if method=='trace' else src
                        trajectory_rows.append((d,b,lid,chrom,pop,mode,method,method_source,edges[i],edges[i+1],(edges[i]+edges[i+1])//2,total,running,fraction,lo,hi))
                phy_carriers=phy_n if mode=='tested' else total
                fraction,lo,hi=wilson_interval(phy_carriers,total)
                for i in range(bins):
                    if edges[i+1]<=core_start or edges[i]>=core_end:continue
                    trajectory_rows.append((d,b,lid,chrom,pop,mode,'phyml',src,edges[i],edges[i+1],(edges[i]+edges[i+1])//2,total,phy_carriers,fraction,lo,hi))
        if len(trajectory_rows)>=100000:
            con.executemany('INSERT INTO locus_trajectory VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',trajectory_rows);trajectory_rows=[]
    if trajectory_rows:con.executemany('INSERT INTO locus_trajectory VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',trajectory_rows)
    con.executescript('''
      CREATE INDEX idx_locus_trajectory_filter ON locus_trajectory(dataset_id,genome_build,locus_id,population,denominator_mode,method,bin_start);
    ''')
    con.commit()


def build_locus_evidence(con):
    """Build a QC-first evidence row without combining incomparable scores."""
    con.executescript('''
      DROP TABLE IF EXISTS locus_evidence;
      CREATE TABLE locus_evidence(
        dataset_id TEXT,genome_build TEXT,locus_id TEXT,chr TEXT,source_class TEXT,
        n_tested INTEGER,n_phyml_carriers INTEGER,carrier_fraction REAL,
        best_archaic TEXT,best_lineage TEXT,n_compared INTEGER,n_match INTEGER,prop_match REAL,
        sequence_information TEXT,candidate_clade_bootstrap REAL,candidate_clade_specificity REAL,
        independent_methods_available INTEGER,independent_methods_supporting INTEGER,independent_method_names TEXT,
        tract_support_mean REAL,ibdmix_max_lod REAL,trace_max_posterior REAL,as3_max_score REAL,
        evidence_score REAL,evidence_completeness REAL,evidence_tier TEXT,
        probability_calibrated INTEGER,score_interpretation TEXT,
        candidate_clade_pass INTEGER,region_coverage_fraction REAL,phyml_qc TEXT,evidence_state TEXT,evidence_summary TEXT
      );
    ''')
    rows=[]
    loci=con.execute('''
      SELECT dataset_id,genome_build,locus_id,chr,source_class,start,end,input_start,input_end,
             best_archaic,best_lineage,n_compared,n_match,prop_match,tree_status,candidate_lineage,tree_has_ancestral_outgroup,
             candidate_clade_pass,candidate_clade_bootstrap,candidate_clade_specificity
      FROM loci WHERE method='phyml' ORDER BY dataset_id,genome_build,locus_id
    ''').fetchall()
    for d,b,lid,chrom,src,start,end,input_start,input_end,best_archaic,best_lineage,n_compared,n_match,prop_match,tree_status,candidate_lineage,has_outgroup,clade_pass,clade_boot,clade_specificity in loci:
        if candidate_lineage:src=source_class(candidate_lineage)
        n_tested,n_carriers=con.execute('''
          SELECT COUNT(*),COALESCE(SUM(phyml),0) FROM locus_sample_support
          WHERE dataset_id=? AND genome_build=? AND locus_id=?
        ''',(d,b,lid)).fetchone()
        carrier_fraction=(n_carriers/n_tested) if n_tested else None
        method_rows=con.execute('''
          SELECT method,phyml_support_fraction FROM locus_method_support
          WHERE dataset_id=? AND genome_build=? AND locus_id=? AND population='ALL'
            AND method!='phyml' AND method_available=1 AND evidence_eligible=1
          ORDER BY method
        ''',(d,b,lid)).fetchall()
        fractions=[float(x[1]) for x in method_rows if x[1] is not None]
        tract_support=sum(fractions)/len(fractions) if fractions else None
        supporting=sum(float(x[1] or 0)>0 for x in method_rows)
        names=','.join(x[0] for x in method_rows)
        native=con.execute('''
          SELECT MAX(CASE WHEN method='ibdmix' THEN score END),
                 MAX(CASE WHEN method='trace' THEN posterior END),
                 MAX(CASE WHEN method='as3' THEN score END)
          FROM segments WHERE dataset_id=? AND genome_build=? AND chr=? AND end>? AND start<?
            AND (method='trace' OR source_class=?)
        ''',(d,b,chrom,start,end,src)).fetchone()
        n_compared=int(n_compared or 0); n_match=int(n_match or 0)
        sequence_information='adequate' if n_compared>=10 else ('low' if n_compared>0 else 'missing')
        input_bp=max(0,int(input_end or 0)-int(input_start or 0))
        selected_bp=max(0,int(end or 0)-int(start or 0))
        region_coverage=(selected_bp/input_bp) if input_bp else None
        if tree_status!='complete':phyml_qc='not_evaluable'
        elif sequence_information!='adequate':phyml_qc='low_sequence_information'
        elif region_coverage is not None and region_coverage<.25:phyml_qc='narrow_selected_region'
        elif int(clade_pass or 0)!=1:phyml_qc='no_supported_archaic_clade'
        elif int(has_outgroup or 0)!=1:phyml_qc='unrooted_tree'
        else:phyml_qc='tree_supported'
        if phyml_qc in {'not_evaluable','low_sequence_information'}:state='not_evaluable'
        elif phyml_qc=='narrow_selected_region':state='requires_core_region_rerun'
        elif phyml_qc=='no_supported_archaic_clade':state='phyml_not_supported'
        elif phyml_qc=='unrooted_tree':state='provisional_tree_support_with_independent_overlap' if supporting else 'provisional_tree_support'
        elif not method_rows:state='tree_supported_no_independent_method'
        elif supporting:state='tree_supported_with_independent_overlap'
        else:state='tree_supported_without_independent_overlap'
        summary=(f'{phyml_qc}; tree carriers={n_carriers}/{n_tested}; '
                 f'eligible callers with overlap={supporting}/{len(method_rows)}')
        rows.append((d,b,lid,chrom,src,n_tested,n_carriers,carrier_fraction,best_archaic,best_lineage,
            n_compared,n_match,prop_match,sequence_information,clade_boot,clade_specificity,
            len(method_rows),supporting,names,tract_support,*native,None,None,state,0,
            'No composite score: inspect tree QC and each caller in its native scale',
            int(clade_pass or 0),region_coverage,phyml_qc,state,summary))
    con.executemany('INSERT INTO locus_evidence VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',rows)
    con.executescript('''
      CREATE UNIQUE INDEX idx_locus_evidence_pk ON locus_evidence(dataset_id,genome_build,locus_id);
      CREATE INDEX idx_locus_evidence_state ON locus_evidence(dataset_id,genome_build,evidence_state,phyml_qc);
    ''')
    con.commit()

def build_database_summary(con):
    con.executescript('''
      DELETE FROM database_summary;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'segment_calls',COUNT(*) FROM segments GROUP BY dataset_id,genome_build;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'samples',COUNT(DISTINCT sample_id) FROM segments
        WHERE sample_id IS NOT NULL AND sample_id!='' GROUP BY dataset_id,genome_build;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'methods',COUNT(DISTINCT method) FROM (
          SELECT dataset_id,genome_build,method FROM segments UNION ALL SELECT dataset_id,genome_build,method FROM loci
        ) GROUP BY dataset_id,genome_build;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'phyml_loci',COUNT(*) FROM loci WHERE method='phyml' GROUP BY dataset_id,genome_build;
      INSERT INTO database_summary
        SELECT 'reference',genome_build,'external_reference_intervals',COUNT(*) FROM reference_callsets
        WHERE reference_role='external_reference' GROUP BY genome_build;
    ''')
    con.commit()

def export_query(con,query,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(f'.{path.name}.part.{os.getpid()}')
    cur=con.execute(query); cols=[x[0] for x in cur.description]
    with gzip.open(tmp,'wt') as h:
        h.write('\t'.join(cols)+'\n')
        while True:
            rows=cur.fetchmany(50000)
            if not rows:break
            for r in rows:h.write('\t'.join('' if x is None else str(x) for x in r)+'\n')
    os.replace(tmp,path)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--analysis-root',required=True,type=Path); ap.add_argument('--output-dir',required=True,type=Path); ap.add_argument('--build',default='GRCh37'); ap.add_argument('--reciprocal-overlap',type=float,default=.5); ap.add_argument('--trajectory-bins',type=int,default=120)
    ap.add_argument('--database',type=Path,help='SQLite destination; defaults to <analysis-root>/gu.sqlite')
    ap.add_argument('--reference-callsets',type=Path,help='Extracted database_introgression_call directory; external validation only')
    ap.add_argument('--reference-cache',type=Path,help='Prepared persistent published-callset SQLite cache')
    ap.add_argument('--reference-dataset-id',default='AS3_1KG')
    ap.add_argument('--sample-panel',type=Path,help='1KG sample metadata with sample/pop/super_pop columns')
    args=ap.parse_args()
    if args.trajectory_bins<20 or args.trajectory_bins>1000:ap.error('--trajectory-bins must be between 20 and 1000')
    if args.reference_callsets and args.reference_cache:ap.error('use --reference-cache or --reference-callsets, not both')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for obsolete in ('README.txt','manifest.tsv'):
        path=args.output_dir/obsolete
        if path.exists():path.unlink()
    packaged=package_core_results(args.analysis_root,args.output_dir)
    parse_root=args.output_dir
    db=(args.database or (args.analysis_root/'gu.sqlite')).resolve(); db.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp_name=tempfile.mkstemp(prefix='gu-normalize-',suffix='.sqlite'); os.close(fd)
    tmp_db=Path(tmp_name); db_part=db.with_name(f'.{db.name}.part.{os.getpid()}')
    con=None; loci=[]; reference_raw_table='temp.reference_callsets_raw'; reference_cache_attached=False
    local_reference_cache=None
    try:
        con=sqlite3.connect(tmp_db); init_db(con); counts={}
        method_runs=discover_method_runs(args.analysis_root,args.build)
        if method_runs:
            con.executemany('INSERT INTO method_runs VALUES(?,?,?,?,?,?,?,?)',method_runs)
        con.execute('CREATE INDEX idx_method_runs_unit ON method_runs(dataset_id,genome_build,chr,method)')
        counts['completed_method_units']=len(method_runs)
        populations=sample_population_table(args.sample_panel)
        if not populations.empty:insert_df(con,'sample_populations',populations)
        nref=0
        if args.reference_cache:
            if not args.reference_cache.is_file():raise ValueError(f'reference cache is missing: {args.reference_cache}')
            cache_fd,cache_name=tempfile.mkstemp(prefix='gu-reference-cache-',suffix='.sqlite'); os.close(cache_fd)
            local_reference_cache=Path(cache_name); shutil.copy2(args.reference_cache,local_reference_cache)
            con.execute('ATTACH DATABASE ? AS refcache',(str(local_reference_cache),)); reference_cache_attached=True
            required={'metadata','reference_callsets','reference_callsets_raw'}
            tables={row[0] for row in con.execute("SELECT name FROM refcache.sqlite_master WHERE type='table'")}
            if not required.issubset(tables):raise ValueError(f'invalid reference cache schema: {args.reference_cache}')
            nref=con.execute('SELECT COUNT(*) FROM refcache.reference_callsets_raw').fetchone()[0]
            con.execute('INSERT INTO reference_callsets SELECT * FROM refcache.reference_callsets')
            reference_raw_table='refcache.reference_callsets_raw'; con.commit()
        elif args.reference_callsets:
            for d in iter_reference_callsets(args.reference_callsets,args.reference_dataset_id):nref+=insert_reference_df(con,d)
        counts['external_reference_raw']=nref; con.commit()
        if not args.reference_cache:collapse_reference_callsets(con)
        counts['external_reference']=con.execute('SELECT COUNT(*) FROM reference_callsets').fetchone()[0]
        for name,it in [('ibdmix',iter_ibdmix(args.analysis_root,args.build)),('trace',iter_trace(args.analysis_root,args.build)),('as3',iter_as3(args.analysis_root,args.build))]:
            n=0
            if it is not None:
                for d in it:n+=insert_df(con,'segments_raw',portable_paths(d,db.parent,['raw_file']))
            counts[name]=n; con.commit()
        loci=list(loci_tables(parse_root,args.build))
        if loci: insert_df(con,'loci',portable_paths(pd.DataFrame(loci),db.parent,['raw_file','tree_file','stats_file','plot_file']))
        n_phyml_carriers=0
        for d in iter_phyml_carriers(parse_root,args.build):
            n_phyml_carriers+=insert_df(con,'phyml_haplotype_carriers',portable_paths(d,db.parent,['raw_file']))
        counts['phyml_haplotype_copies']=n_phyml_carriers
        con.execute('CREATE INDEX idx_loci_region ON loci(dataset_id,genome_build,chr,start,end)'); con.commit()
        build_catalog(con,args.build,args.reciprocal_overlap)
        build_reference_overlaps(con,reference_raw_table)
        if reference_cache_attached:
            con.execute('DETACH DATABASE refcache'); reference_cache_attached=False
        build_burden(con)
        build_locus_support(con,args.trajectory_bins)
        build_locus_evidence(con)
        build_database_summary(con)
        summary=args.output_dir/'summary'
        export_query(con,'SELECT * FROM segments ORDER BY dataset_id,genome_build,chr,start,end,sample_id',summary/'segments.tsv.gz')
        export_query(con,'SELECT * FROM segment_catalog ORDER BY dataset_id,genome_build,chr,start,end',summary/'segment_catalog.tsv.gz')
        export_query(con,'SELECT * FROM carriers ORDER BY dataset_id,sample_id,segment_code',summary/'carriers.tsv.gz')
        export_query(con,'SELECT * FROM loci ORDER BY dataset_id,genome_build,chr,start,end',summary/'loci.tsv.gz')
        export_query(con,'SELECT * FROM phyml_haplotype_carriers ORDER BY dataset_id,genome_build,locus_id,sample_id,haplotype',summary/'phyml_haplotype_carriers.tsv.gz')
        export_query(con,'SELECT * FROM locus_sample_support ORDER BY dataset_id,genome_build,locus_id,n_methods DESC,sample_id',summary/'locus_sample_support.tsv.gz')
        export_query(con,'SELECT * FROM locus_method_support ORDER BY dataset_id,genome_build,locus_id,population,method',summary/'locus_method_support.tsv.gz')
        export_query(con,'SELECT * FROM locus_evidence ORDER BY dataset_id,genome_build,locus_id',summary/'locus_evidence.tsv.gz')
        export_query(con,'SELECT * FROM locus_trajectory ORDER BY dataset_id,genome_build,locus_id,population,denominator_mode,method,bin_start',summary/'locus_trajectory.tsv.gz')
        export_query(con,'SELECT * FROM method_runs ORDER BY dataset_id,genome_build,chr,method',summary/'method_runs.tsv.gz')
        export_query(con,'SELECT * FROM reference_callsets ORDER BY genome_build,population,chr,start,end',summary/'reference_callsets.tsv.gz')
        export_query(con,'SELECT * FROM reference_callset_overlaps ORDER BY genome_build,reference_population,chr,segment_start,sample_id',summary/'reference_callset_overlaps.tsv.gz')
        export_query(con,'SELECT * FROM sample_populations ORDER BY dataset_id,population,sample_id',summary/'sample_populations.tsv.gz')
        export_query(con,'SELECT * FROM sample_burden ORDER BY dataset_id,genome_build,burden_type,sample_id,method,source_class',summary/'sample_burden.tsv.gz')
        for method in ('as3','trace','ibdmix'):
            export_query(con,f"SELECT * FROM segments WHERE method='{method}' ORDER BY dataset_id,genome_build,chr,start,end,sample_id",args.output_dir/method/'segments.tsv.gz')
        export_query(con,"SELECT * FROM loci WHERE method='phyml' ORDER BY dataset_id,genome_build,chr,start,end",args.output_dir/'phyml/loci.tsv.gz')
        check=con.execute('PRAGMA quick_check').fetchone()
        if not check or check[0] != 'ok':raise RuntimeError(f'SQLite quick_check failed: {check}')
        con.close(); con=None
        if db_part.exists():db_part.unlink()
        shutil.copy2(tmp_db,db_part); os.replace(db_part,db)
    finally:
        if con is not None:con.close()
        if tmp_db.exists():tmp_db.unlink()
        if db_part.exists():db_part.unlink()
        if local_reference_cache is not None and local_reference_cache.exists():local_reference_cache.unlink()
    print('Normalized:',counts,'loci=',len(loci),'packaged_files=',packaged,'database=',db)

if __name__=='__main__':main()
