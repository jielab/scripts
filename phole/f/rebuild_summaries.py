#!/usr/bin/env python3
"""Regenerate certified PHOLE summaries by reading every supplied full-GWAS row.
No LD-independence claim: trans loci are groups of significant positions with
successive gaps <=1 Mb. Original numeric strings are retained; no COJO values.
"""
from __future__ import annotations
import argparse, concurrent.futures, csv, gzip, hashlib, json, math, os, subprocess, tempfile
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
import cross_check as cc
from summary_metadata import read_metadata, metadata_path

FIELDS=['SNP','CHR','POS','BETA','SE','P','EA','NEA','EAF','N','LOG10P','GENE','CLASS']
AWK=r'''
BEGIN { FS="[ \t]+"; OFS="\t" }
FILENAME!="-" { if ($0 !~ /^#/ && NF>=3) { nc++; c[nc]=$1; lo[nc]=$2; hi[nc]=$3 }; next }
!header { sub(/^#/,"",$1); for(i=1;i<=NF;i++) h[toupper($i)]=i;
 for(i=1;i<=6;i++) { split("SNP CHR POS BETA SE P",need," "); if(!h[need[i]]) {print "Missing column " need[i] > "/dev/stderr"; exit 2} }
 header=1; print $0; next }
/^#/ || NF==0 { next }
{ rows++; chr=$(h["CHR"]); sub(/^[cC][hH][rR]/,"",chr); if(chr==23)chr="X"; if(chr==24)chr="Y"; if(chr==25 || chr=="M")chr="MT";
 pos=$(h["POS"]); p=$(h["P"]);
 if(NF!=length(h) || pos !~ /^[0-9]+$/ || pos+0<1 || p !~ /^[+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ || p+0<0 || p+0>1) {bad++;next}
 cis=0; for(i=1;i<=nc;i++)if(chr==c[i] && pos+0>=lo[i] && pos+0<=hi[i])cis=1;
 if(cis)cisrows++;
 if(p+0<=0.00000005) { if(cis)cissig++; print $0 }
}
END { print rows+0,cisrows+0,cissig+0,bad+0 > stats; if(!header || bad || !rows)exit 3 }
'''
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def identity(p):
 s=p.stat();return [str(p.resolve()),s.st_size,s.st_mtime_ns]
def table(path,rows,fields):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def work(task):
 dataset,assay,source,sig,intervals,flank,code,old,*original=task
 source=Path(source);sig=Path(sig)
 before=identity(source)
 key=hashlib.sha256(json.dumps([before,intervals,flank,code],sort_keys=True).encode()).hexdigest()
 if sig.is_file() and old:
  if old.get('summary_rebuild_key')==key and old.get('sig_md5')==hashlib.md5(sig.read_bytes()).hexdigest():return dataset,assay,'reused',old
 # Everything is prepared before the old summary is replaced. Metadata is installed last.
 sig.parent.mkdir(parents=True,exist_ok=True)
 with tempfile.TemporaryDirectory(prefix='.phole-rebuild-',dir=sig.parent) as tmp:
  tmp=Path(tmp);bounds=tmp/'bounds';stats=tmp/'stats';hits=tmp/'hits';errors=tmp/'stderr'
  bounds.write_text('# 1-based inclusive\n'+''.join(f'{c}\t{lo}\t{hi}\n' for c,lo,hi in intervals))
  with hits.open('w') as out,errors.open('w') as err:
   z=subprocess.Popen(['gzip','-dc','--',str(source)],stdout=subprocess.PIPE,stderr=err)
   try:
    result=subprocess.run(['awk','-v','stats='+str(stats),AWK,str(bounds),'-'],stdin=z.stdout,stdout=out,stderr=err)
   finally:z.stdout.close()
   zrc=z.wait()
  if result.returncode or zrc:raise ValueError(f'{dataset}/{assay}: full scan failed: gzip={zrc}, awk={result.returncode}; '+errors.read_text()[-1500:])
  counts=list(map(int,stats.read_text().split()))
  if before!=identity(source):raise ValueError(f'{dataset}/{assay}: input changed during full scan')
  rows=[]
  with hits.open() as f:
   header=f.readline().strip().lstrip('#').upper().split()
   for line in f:
    fields=line.split();raw=dict(zip(header,fields))
    p=Decimal(raw['P']);pos=int(raw['POS']);chrom=cc.chrom(raw['CHR'])
    if not p.is_finite() or not 0<=p<=Decimal('5e-8'):continue # recheck awk floating boundary
    if chrom not in cc.PRIMARY_CHROMS:continue
    beta,se=Decimal(raw['BETA']),Decimal(raw['SE'])
    if not beta.is_finite() or not se.is_finite() or se<=0:raise ValueError(f'{dataset}/{assay}: invalid significant BETA/SE')
    row={k:raw.get(k,'.') for k in FIELDS};row.update(CHR=chrom,POS=str(pos),GENE='.')
    row['_key']=(p,pos,raw['SNP']);row['_pos']=pos;rows.append(row)
  # Select an own-cis lead per distinct component window, and all trans position clusters.
  selected={};cis_ids=set()
  for c,lo,hi in intervals:
   indices=[i for i,r in enumerate(rows) if r['CHR']==c and lo<=r['_pos']<=hi]
   cis_ids.update(indices)
   if indices:
    i=min(indices,key=lambda i:rows[i]['_key']);selected[i]='cis'
  groups=defaultdict(list)
  for i,r in enumerate(rows):
   if i not in cis_ids:groups[r['CHR']].append(i)
  for c,indices in groups.items():
   indices.sort(key=lambda i:(rows[i]['_pos'],rows[i]['_key']))
   group=[];last=None
   for i in indices:
    if last is not None and rows[i]['_pos']-last>1000000:
     selected[min(group,key=lambda j:rows[j]['_key'])]='trans';group=[]
    group.append(i);last=rows[i]['_pos']
   if group:selected[min(group,key=lambda j:rows[j]['_key'])]='trans'
  result_rows=[]
  for i,kind in sorted(selected.items(),key=lambda kv:(cc.chrom_key(rows[kv[0]]['CHR']),rows[kv[0]]['_pos'],rows[kv[0]]['SNP'])):
   r={k:rows[i][k] for k in FIELDS};r.update(CLASS=kind,GENE=assay if kind=='cis' else '.');result_rows.append(r)
  magma='magma:\nGENE\tCHR\tREGION\tZSTAT\tP\tTRAIT\tCLASS\n'
  previous=sig if sig.exists() or not original else Path(original[0])
  if previous.exists():
   text=previous.read_text()
   if '\nmagma:' in text:magma=text.split('\nmagma:',1)[1];magma='magma:'+magma
  new=tmp/'summary'
  with new.open('w',newline='') as f:
   f.write('common:\n');w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(result_rows);f.write('\n'+magma)
  md={'ASSAY':assay,'schema_version':2,'trait':assay,'dataset':dataset,'build':'38','build_source':'user_confirmed_GRCh38_2026-11-07',
      'common_rows':len(result_rows),'sig_md5':hashlib.md5(new.read_bytes()).hexdigest(),'p_threshold':'5e-8',
      'cis_flank_bp':flank,'cis_variant_rows':counts[1],'cis_significant_rows':len(cis_ids),
      'cis_intervals':';'.join(f'{c}:{lo}-{hi}' for c,lo,hi in intervals),'cis_annotation_available':str(bool(intervals)).lower(),
      'cis_annotation_build':'38','cis_window':'gene_body_flank','trans_selection':'all_position_clumped_loci','trans_gap_bp':1000000,
      'association_scope':'full_gwas','association_input':str(source.resolve()),'p_source':'raw_gwas','beta_se_source':'raw_gwas',
      'full_input_rows':counts[0],'malformed_rows':counts[3],'source_size':before[1],'source_mtime_ns':before[2],
      'source_identity_policy':'absolute path,size,mtime_ns; full stream read to gzip EOF; not a full content checksum',
      'summary_rebuild_key':key,'generator_sha256':code,'ld_clumping':'false','preserved_magma_section':'not regenerated; not used by PHOLE'}
  os.replace(new,sig)
 return dataset,assay,'rebuilt',md
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--gwas-root',type=Path,required=True);p.add_argument('--cross-dir',type=Path,required=True)
 p.add_argument('--summary-root',type=Path,help='Destination for regenerated cohort summaries')
 p.add_argument('--out',type=Path,required=True);p.add_argument('--cores',type=int,default=40);p.add_argument('--cis-flank',type=int,default=500000)
 p.add_argument('--datasets',default='ukb,ckb,decode,scot',help='Explicit cohort subset, for partial recovery or bounded validation')
 p.add_argument('--dataset-builds',required=True);p.add_argument('--limit',type=int,default=0,help='Explicit smoke-test limit; never certifies unprocessed summaries')
 a=p.parse_args();builds=dict(x.split('=',1) for x in a.dataset_builds.split(',') if x)
 if any(builds.get(d)!='38' for d in cc.DATASETS):raise ValueError('Rebuilding requires explicit GRCh38 for all four cohorts; no liftover')
 selected=set(a.datasets.split(','))
 if not selected or not selected.issubset(cc.DATASETS):raise ValueError('Invalid dataset subset')
 a.out.mkdir(parents=True,exist_ok=True);(a.out/'summary.json').unlink(missing_ok=True);tasks=[];missing=[];code=sha(__file__);cohort_meta={};folders={}
 pending=[]
 for d,folder in cc.DATASETS.items():
  if d not in selected:continue
  beds=defaultdict(list)
  for line in (a.cross_dir/'cis-bed'/f'{d}.GRCh38.bed').read_text().splitlines():
   if line.startswith('#') or not line.strip():continue
   c,lo,hi,assay=line.split()[:4];beds[assay].append((cc.chrom(c),max(1,int(lo)-a.cis_flank),int(hi)+a.cis_flank))
  root=a.gwas_root/folder;folders[d]=(a.summary_root or a.gwas_root)/folder;cohort_meta[d]=dict(read_metadata(folders[d]))
  existing={x.name.removesuffix('.sig.txt') for x in (root/'mplot').glob('*.sig.txt')}
  aliases=set(existing)
  if (root/'common').is_dir():aliases.update(x.name for x in (root/'common').iterdir())
  if (root/'clean').is_dir():aliases.update(x.name.removesuffix('.gz') for x in (root/'clean').glob('*.gz'))
  print(f'{d}: checking {len(aliases)} assays, {len(existing)} existing summaries',flush=True)
  for assay in sorted(aliases):pending.append((d,assay,root,assay in existing,sorted(set(beds[assay]))))
 def resolve(item):
  d,assay,root,existing,intervals=item
  choices=[x for x in [root/'common'/assay/'gwas'/(assay+'.gz'),root/'clean'/(assay+'.gz')] if x.is_file()]
  if len(choices)!=1:
   return None,({'dataset':d,'assay':assay,'reason':'missing_or_ambiguous_complete_GWAS'} if existing or choices else None)
  return (d,assay,str(choices[0]),str(folders[d]/'mplot'/(assay+'.sig.txt')),intervals,a.cis_flank,code,cohort_meta[d].get(assay,{}),str(root/'mplot'/(assay+'.sig.txt'))),None
 with concurrent.futures.ThreadPoolExecutor(max_workers=min(a.cores,40)) as pool:
  for task,error in pool.map(resolve,pending):
   if task:tasks.append(task)
   if error:missing.append(error)
 (a.out/'missing-inputs.json').write_text(json.dumps(missing,indent=2))
 if missing:raise ValueError(f'{len(missing)} summaries have missing/ambiguous full inputs; see {a.out}/missing-inputs.json')
 def flush():
  for d,rows in cohort_meta.items():
   if not rows:continue
   dest=metadata_path(folders[d]);temporary=dest.with_name(dest.name+'.tmp')
   fields=list(dict.fromkeys(k for row in rows.values() for k in row))
   table(temporary,[rows[k] for k in sorted(rows)],fields);os.replace(temporary,dest)
 if a.limit:tasks=tasks[:a.limit]
 (a.out/'input-plan.json').write_text(json.dumps({'tasks':len(tasks),'compressed_bytes':sum(Path(t[2]).stat().st_size for t in tasks),'limit':a.limit,'workers':max(1,a.cores//2)},indent=2))
 print(f'Rebuilding {len(tasks)} full files with {max(1,a.cores//2)} gzip/awk worker pairs',flush=True)
 counts=defaultdict(int);failures=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=max(1,a.cores//2)) as pool:
  futures={pool.submit(work,t):t[:2] for t in tasks}
  for i,future in enumerate(concurrent.futures.as_completed(futures),1):
   try:
    d,assay,status,md=future.result();counts[status]+=1;cohort_meta[d][assay]=md
   except Exception as e:failures.append({'assay':futures[future],'error':str(e)});print('ERROR',failures[-1],flush=True)
   if i%25==0 or i==len(tasks):
    flush();print(f'{i}/{len(tasks)} {dict(counts)} failed={len(failures)}',flush=True)
 report={'summary_root':str((a.summary_root or a.gwas_root).resolve()),'datasets':sorted(selected),'tasks':len(tasks),'counts':dict(counts),'failures':failures,'limit':a.limit}
 (a.out/'summary.json').write_text(json.dumps(report,indent=2))
 if failures:raise ValueError(f'{len(failures)} rebuild failures; no certification was emitted for these files')
if __name__=='__main__':main()
