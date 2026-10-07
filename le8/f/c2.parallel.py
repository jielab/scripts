#!/usr/bin/env python3
"""Parallel adapter around the existing LE8 serial MR-link-2 runner.

Deploy beside c2.cause.sh + bounded_jobs.py + mrlink_task_adapter.py. Invoke this
adapter instead of the serial shell loop, or add the documented delegation hook.
Every child invokes the existing scientific implementation in an isolated folder.
No changes to the MR likelihood, harmonization, SNPs, or p-value corrections.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,os,re,shutil,sys,tempfile,time,uuid
from pathlib import Path
import importlib.util
spec=importlib.util.spec_from_file_location('le8_resources',Path(__file__).with_name('0.resources.py'))
resources=importlib.util.module_from_spec(spec);sys.modules.setdefault('le8_resources',resources);spec.loader.exec_module(resources)
execute,atomic_json,digest_file=resources.execute,resources.atomic_json,resources.digest_file


def fingerprint(path,mode,memo):
    p=Path(path).resolve(strict=True);st=p.stat();key=(str(p),st.st_size,st.st_mtime_ns)
    if key not in memo:
        memo[key]={'path':str(p),'size':st.st_size,'mtime_ns':st.st_mtime_ns,'mode':mode}
        if mode=='sha256':memo[key]['sha256']=digest_file(p)
    return memo[key]


def key_of(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()


def numeric_module():
    spec=importlib.util.spec_from_file_location('le8_mr_numeric',Path(__file__).with_name('c2.mr_link2.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def reference_identity_files(rows,forward):
    def option(flag,variable,default=''):
        return forward[forward.index(flag)+1] if flag in forward else os.getenv(variable,default)
    bed=option('--reference-bed','MRLINK2_REF_BED')
    if bed:
        files=[Path(bed+'.'+ext) for ext in ('bed','bim','fam')]
        if not all(p.is_file() for p in files):raise ValueError('Incomplete explicit reference BED')
        return list(map(str,files))
    directory=Path(option('--reference-pfile-dir','MRLINK2_REF_PFILE_DIR',f"/mnt/f/gen/1kg/{os.getenv('LE8_GRCH','37')}/pfile"))
    pop=option('--reference-pop','MRLINK2_REF_POP','EUR').upper();files=[]
    for chrom in sorted({r['region'].split(':')[0].removeprefix('chr') for r in rows}):
        candidates=[(directory.parent/'bfile'/pop/f'chr{chrom}','bed'),
            (directory/f'{pop}.chr{chrom}','bed'),(directory/f'chr{chrom}','bed')] if pop=='ALL' else [
            (directory.parent/'bfile'/pop/f'chr{chrom}','bed'),(directory/f'{pop}.chr{chrom}','bed')]
        candidates.insert(0,(directory/pop/f'chr{chrom}','bed'))
        candidates += [(directory/f'{pop}.chr{chrom}','pgen'),(directory/f'chr{chrom}','pgen'),(directory/f'chr{chrom}','bed')]
        found=False
        for prefix,kind in candidates:
            if kind=='bed':parts=[Path(str(prefix)+'.'+ext) for ext in ('bed','bim','fam')]
            else:
                pvar=Path(str(prefix)+'.pvar');pvar=pvar if pvar.exists() else Path(str(prefix)+'.pvar.zst')
                parts=[Path(str(prefix)+'.pgen'),Path(str(prefix)+'.psam'),pvar]
            if all(p.is_file() for p in parts):files+=parts;found=True;break
        if not found:raise ValueError('No complete reference for chromosome '+chrom)
    ids=Path(option('--reference-id-dir','MRLINK2_REF_ID_DIR') or directory.parent/'id')
    samples=Path(option('--reference-samples','MRLINK2_REF_SAMPLES') or directory.parent/'samples.txt')
    keep=ids/(pop+'.id.2col')
    if pop!='ALL':
        if keep.is_file():files.append(keep)
        elif samples.is_file():files.append(samples)
    return sorted(set(map(str,files)))


def batch_region_inputs(rows,a):
    """One streaming pass per unindexed input for all missing requested regions.

    Raw columns and every in-region row are retained. No SNP or allele selection.
    The scientific task identity continues to reference the original input.
    """
    import gzip,bisect,collections,fcntl
    selected={};requests={};memo={}
    for row in rows:
        for field in ('exposure','outcome'):
            file=str(Path(row[field] or a.cad_gwas).resolve())
            if not Path(file).is_file():continue
            # Fresh indexes are left to the existing region-aware standardizer.
            fresh=any(Path(file+ext).is_file() and Path(file+ext).stat().st_mtime_ns>=Path(file).stat().st_mtime_ns for ext in ('.tbi','.csi'))
            if not fresh:requests.setdefault(file,set()).add(row['region'])
    for file,regions in requests.items():
        identity=key_of(dict(source=fingerprint(file,a.hash_mode,memo),extractor='raw-region-v1'))
        cache=a.cache_root/'_regions'/identity;cache.mkdir(parents=True,exist_ok=True)
        with open(cache/'producer.lock','a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            missing=[]
            for region in sorted(regions):
                dest=cache/(key_of(region)+'.csv.gz');proof=dest.with_suffix('.json')
                valid=False
                try:valid=json.loads(proof.read_text())['sha256']==digest_file(dest)
                except (OSError,ValueError,KeyError):pass
                if not valid:missing.append(region)
                selected[file,region]=str(dest)
            if not missing:continue
            intervals={};targets={};writers=collections.OrderedDict();counts={region:0 for region in missing}
            def writer(region,header):
                if region in writers:writers.move_to_end(region);return writers[region][1]
                if len(writers)>=32:_,(handle,_)=writers.popitem(last=False);handle.close()
                path=targets[region];empty=not path.exists()
                handle=gzip.open(path,'at',newline='');obj=csv.DictWriter(handle,fieldnames=header)
                if empty:obj.writeheader()
                writers[region]=(handle,obj);return obj
            for region in missing:
                chrom,span=region.split(':');lo,hi=map(int,span.split('-'))
                if hi<lo:raise ValueError('Invalid region bounds')
                intervals.setdefault(chrom.removeprefix('chr'),[]).append((lo,hi,region))
                targets[region]=cache/(key_of(region)+'.tmp-'+uuid.uuid4().hex+'.gz')
            spans={chrom:max(hi-lo for lo,hi,_ in items) for chrom,items in intervals.items()}
            for items in intervals.values():items.sort()
            starts={chrom:[v[0] for v in items] for chrom,items in intervals.items()}
            opener=gzip.open if file.endswith('.gz') else open
            print('[LE8 MR] Single input pass: '+file+'; regions='+str(len(missing)),flush=True)
            try:
                with opener(file,'rt',newline='') as handle:
                    headerline=handle.readline();delimiter='\t' if headerline.count('\t')>=headerline.count(',') else ','
                    header=next(csv.reader([headerline],delimiter=delimiter));upper={v.upper():v for v in header}
                    cc=next((upper[k] for k in ('CHROMOSOME','CHR','CHROM') if k in upper),None)
                    pc=next((upper[k] for k in ('POSITION','POS','BP','BASE_PAIR_LOCATION','BASEPAIR') if k in upper),None)
                    if cc is None or pc is None:raise ValueError('Missing chromosome/position for region extraction: '+file)
                    for region in missing:writer(region,header)
                    for row in csv.DictReader(handle,fieldnames=header,delimiter=delimiter):
                        chrom=re.sub(r'\.0$','',str(row[cc]).lower().removeprefix('chr'))
                        if chrom not in intervals:continue
                        try:pos=int(float(row[pc]))
                        except (TypeError,ValueError):continue
                        lo=bisect.bisect_left(starts[chrom],pos-spans[chrom]);hi=bisect.bisect_right(starts[chrom],pos)
                        for _,end,region in intervals[chrom][lo:hi]:
                            if pos<=end:writer(region,header).writerow(row);counts[region]+=1
                for handle,_ in writers.values():handle.close()
                writers.clear()
                for region,path in targets.items():
                    dest=Path(selected[file,region]);os.replace(path,dest)
                    atomic_json(dest.with_suffix('.json'),dict(source=fingerprint(file,a.hash_mode,memo),region=region,rows=counts[region],sha256=digest_file(dest)))
            finally:
                for handle,_ in writers.values():handle.close()
                for path in targets.values():path.unlink(missing_ok=True)
    return selected


def match_main(argv):
    import fcntl,gzip
    if len(argv)!=4:raise ValueError('match reference input output audit')
    reference,source,output,audit=argv
    phe=Path(os.getenv('PHE_F','/mnt/d/scripts/0f/phenotype.sh')).resolve(strict=True)
    h=hashlib.sha256()
    with gzip.open(source,'rb') if source.endswith('.gz') else open(source,'rb') as handle:
        for block in iter(lambda:handle.read(2**20),b''):h.update(block)
    identity=dict(summary_sha256=h.hexdigest(),reference=fingerprint(reference+'.bim','stat',{}),
        matcher=dict(path=str(phe),sha256=digest_file(phe)),build=os.getenv('LE8_GRCH','37'),policy='match_GRCH-v1')
    cache=Path(os.getenv('MRLINK2_PREPARE_CACHE','/tmp/le8-mrlink2-prepared'))/'aligned'/key_of(identity)
    cache.mkdir(parents=True,exist_ok=True);dest=cache/'matched.tsv';audit_dest=cache/'audit.tsv'
    with open(cache/'lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);valid=False
        try:
            proof=json.loads((cache/'complete.json').read_text())
            valid=proof['sha256']==digest_file(dest) and proof['audit_sha256']==digest_file(audit_dest)
        except (OSError,ValueError,KeyError):pass
        if not valid:
            with tempfile.TemporaryDirectory(prefix='match-',dir=cache) as tmp:
                target=Path(tmp)/'matched.tsv';record=Path(tmp)/'audit.tsv'
                command=['bash','-c','set -euo pipefail; source "$1"; match_GRCH --reference "$2" --output "$4" --audit "$5" "$3"',
                         'le8-match',str(phe),reference,source,str(target),str(record)]
                subprocess.run(command,check=True)
                if not target.is_file() or not record.is_file():raise ValueError('Missing matched outputs')
                os.replace(target,dest);os.replace(record,audit_dest)
                atomic_json(cache/'complete.json',dict(identity=identity,sha256=digest_file(dest),audit_sha256=digest_file(audit_dest)))
        shutil.copyfile(dest,output);shutil.copyfile(audit_dest,audit)
    return 0

def prepare_main(argv):
    import fcntl,gzip
    if len(argv) not in (3,4):raise ValueError('prepare input output verified_N_or_nan [region]')
    source,output,n=argv[:3];region=argv[3] if len(argv)>3 else None
    code=Path(__file__).with_name('c2.mr_link2.py')
    identity=dict(input=fingerprint(source,'stat',{}),region=region,verified_default_N=n,numeric_code=digest_file(code))
    cache=Path(os.getenv('MRLINK2_PREPARE_CACHE','/tmp/le8-mrlink2-prepared'))/key_of(identity)
    cache.mkdir(parents=True,exist_ok=True);dest=cache/'sumstats.tsv.gz'
    with open(cache/'lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);valid=False
        try:valid=json.loads((cache/'complete.json').read_text())['sha256']==digest_file(dest)
        except (OSError,ValueError,KeyError):pass
        if not valid:
            tmp=cache/('prepared-'+uuid.uuid4().hex+'.gz')
            try:
                numeric_module().le8_prepare_sumstats(source,str(tmp),float(n),region)
                with gzip.open(tmp,'rb') as f:
                    while f.read(1024*1024):pass
                os.replace(tmp,dest);atomic_json(cache/'complete.json',dict(identity=identity,sha256=digest_file(dest)))
            finally:tmp.unlink(missing_ok=True)
        shutil.copyfile(dest,output)
    return 0


def build_tasks(a,forward):
    with a.jobs.open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t');header=reader.fieldnames;rows=list(reader)
    if not header or not {'omics','trait','exposure','outcome','region'}<=set(header):raise ValueError('Invalid jobs header')
    script=a.worker_script.resolve(strict=True);implementation=script.with_name('c2.mr_link2.py')
    code={str(p):digest_file(p) for p in (script,implementation,Path(__file__),Path(__file__).with_name('0.resources.py')) if p.is_file()}
    memo={};tasks=[];seen=set();root=a.cache_root.resolve();root.mkdir(parents=True,exist_ok=True)
    if not a.identity_file:a.identity_file=reference_identity_files(rows,forward)
    extracted=batch_region_inputs(rows,a)
    for i,row in enumerate(rows):
        if not row['trait']:continue
        exposure=Path(row['exposure']);outcome=Path(row['outcome'] or a.cad_gwas)
        if not re.fullmatch(r'(?:chr)?[^:]+:[0-9]+-[0-9]+',row['region']):raise ValueError('Prespecified region required for bounded MR jobs')
        raw=dict(row);raw['outcome']=str(outcome.resolve());raw['exposure']=str(exposure.resolve())
        source=[fingerprint(exposure,a.hash_mode,memo),fingerprint(outcome,a.hash_mode,memo)]
        for file in a.identity_file:source.append(fingerprint(file,a.hash_mode,memo))
        identity=dict(job=raw,inputs=source,worker_code=code,forwarded_options=forward,
                      environment={k:v for k,v in os.environ.items() if k.startswith(('MRLINK2_','LE8_GRCH','PHE_F','LE8_REFGEN')) and k not in {'MRLINK2_WORKERS','MRLINK2_INNER_THREADS','MRLINK2_JOB_MEMORY_GIB'}},
                      signature_policy='explicit_manifest_v1')
        key=key_of(identity)
        if key in seen:raise ValueError('Duplicate complete MR job identity')
        seen.add(key);job=root/key;job.mkdir(exist_ok=True)
        taskpath=job/'executor_task.json';existing=None
        if taskpath.is_file() and os.getenv('LE8_REPLACE','FALSE').upper() not in ('TRUE','1','YES'):
            candidate=json.loads(taskpath.read_text())
            try:
                validate_outputs=resources.validate_outputs
                validate_outputs(candidate);existing=candidate
            except (OSError,ValueError,KeyError):pass
        if existing is not None:
            existing['cores']=a.inner_threads
            tasks.append(existing);continue
        attempt=job/('attempt-'+uuid.uuid4().hex);attempt.mkdir()
        single=attempt/'job.tsv'
        with single.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=header,delimiter='\t');writer.writeheader();writer.writerow({**raw,'exposure':extracted.get((str(exposure.resolve()),row['region']),str(exposure.resolve())),
                'outcome':extracted.get((str(outcome.resolve()),row['region']),str(outcome.resolve()))})
        worker_out=attempt/'results'
        argv=['bash',str(script),'mr-link2','--jobs',str(single),'--cad-gwas',str(outcome.resolve()),
              '--outdir',str(worker_out),*forward]
        marker=attempt/'validated.json';request=attempt/'worker.json'
        atomic_json(request,dict(task_id=key,omics=row['omics'],trait=row['trait'],scientific_signature=key,
                                worker_outdir=str(worker_out),argv=argv))
        variants=0
        try:variants=int(float(row.get('region_variants','0')))
        except (ValueError,OverflowError):pass
        # Conservative accounting hint, NOT an LD shape cap. Never discard SNPs.
        estimated=max(a.job_memory_gib, a.ld_workspace_multiplier*8*variants*variants/2**30+a.job_base_gib)
        task=dict(task_id=key,argv=[sys.executable,str(Path(__file__)),'worker',
                      '--task',str(request),'--marker',str(marker)],scientific_signature=key,
                  outputs=[dict(path=str(marker),kind='mrlink_marker')],cores=a.inner_threads,memory_gib=estimated,
                  env={'LE8_MRLINK2_WORKER':'1','LE8_REPLACE':'FALSE',**({'CUDA_VISIBLE_DEVICES':''} if os.getenv('MRLINK2_EIGH_BACKEND','cpu')=='cpu' else {}),
                       'MRLINK2_EXPOSURE_N':row.get('exposure_n') or os.getenv('MRLINK2_EXPOSURE_N','nan'),
                       'MRLINK2_OUTCOME_N':row.get('outcome_n') or os.getenv('MRLINK2_OUTCOME_N','nan')},
                  job=raw,source_position=i)
        atomic_json(taskpath,task);tasks.append(task)
    return tasks


def collect(tasks,status,outdir):
    outdir.mkdir(parents=True,exist_ok=True);(outdir/'results').mkdir(exist_ok=True)
    records=[];all_rows=[];columns=None;complete=True
    for task,run in zip(tasks,status):
        row=task['job'];entry=dict(omics=row['omics'],trait=row['trait'],exposure=row['exposure'],outcome=row['outcome'],
                                  status='failed',out_prefix='',message=run.get('reason',''),job_id=task['task_id'],region=row['region'])
        if run['status']=='completed':
            marker=json.loads(Path(task['outputs'][0]['path']).read_text());entry['status']=marker['status'];entry['message']=marker['row'].get('message','validated')
            if marker['status']=='ok':
                raw=Path(marker['row']['out_prefix'])
                if digest_file(raw)!=marker['row']['result_sha256']:raise ValueError('MR source changed during collection')
                dest=outdir/'results'/(task['task_id']+'.mrlink2');shutil.copyfile(raw,dest);entry['out_prefix']=str(dest.resolve())
                with raw.open(newline='') as f:
                    reader=csv.DictReader(f,delimiter='\t');cols=reader.fieldnames;data=list(reader)
                if columns is None:columns=cols
                if cols!=columns:raise ValueError('Inconsistent MR output schema')
                for d in data:all_rows.append({'omics':row['omics'],'trait':row['trait'],**d,'job_id':task['task_id']})
        else:complete=False
        records.append(entry)
    def tsv(name,fields,rows):
        target=outdir/name;fd,tmp=tempfile.mkstemp(prefix='.'+name,dir=outdir)
        with os.fdopen(fd,'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
        os.replace(tmp,target)
    tsv('mrlink2.status.tsv',['omics','trait','exposure','outcome','status','out_prefix','message','job_id','region'],records)
    tsv('mrlink2.all.tsv',['omics','trait',*(columns or []),'job_id'],all_rows)
    marker=outdir/'mrlink2.complete'
    if complete:atomic_json(marker,{'task_signatures':[t['scientific_signature'] for t in tasks],'count':len(tasks),'all_terminal_validated':True})
    else:marker.unlink(missing_ok=True)
    return complete


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--jobs',type=Path,required=True);p.add_argument('--cad-gwas',required=True);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--worker-script',type=Path,default=Path(__file__).with_name('c2.cause.sh'))
    p.add_argument('--cores',type=int,default=int(os.getenv('LE8_TOTAL_CORES',os.getenv('N_CORES','16'))))
    p.add_argument('--workers',type=int,default=int(os.getenv('MRLINK2_WORKERS','4')))
    p.add_argument('--inner-threads',type=int,default=int(os.getenv('MRLINK2_INNER_THREADS','1')))
    p.add_argument('--memory-gib',type=float,default=float(os.getenv('LE8_RESOURCE_MEMORY_GIB','32')))
    p.add_argument('--job-memory-gib',type=float,default=float(os.getenv('LE8_CPU_TASK_MEMORY_GIB',os.getenv('MRLINK2_JOB_MEMORY_GIB','4'))))
    p.add_argument('--job-base-gib',type=float,default=1);p.add_argument('--ld-workspace-multiplier',type=float,default=8)
    p.add_argument('--cache-root',type=Path,default=Path('/tmp/le8-mrlink2-validated'))
    p.add_argument('--hash-mode',choices=['stat','sha256'],default='stat')
    p.add_argument('--identity-file',action='append',default=[],help='LD/keep/build manifest or actual files to fingerprint; REQUIRED')
    p.add_argument('--pool');p.add_argument('--reserve-gib',type=float,default=4)
    a,forward=p.parse_known_args()
    if a.inner_threads<1 or min(a.job_memory_gib,a.job_base_gib,a.ld_workspace_multiplier)<=0:raise ValueError('Invalid resource values')
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/'mrlink2.complete').unlink(missing_ok=True)
    a.pool=a.pool or os.getenv('LE8_RESOURCE_POOL')
    tasks=build_tasks(a,forward)
    if not tasks:raise ValueError('No MR jobs')
    executor=a.cache_root/('_executor_'+key_of(str(a.outdir.resolve()))[:20]);atomic_json(executor/'tasks.json',tasks)
    with resources.phase_worker_budget(a.cores,a.memory_gib,a.pool,a.reserve_gib) as (cores,memory,pool):
        status=execute(tasks,executor,cores,memory,a.workers,pool,a.reserve_gib)
    return 0 if collect(tasks,status,a.outdir) else 2

import subprocess,math
def validate_status(status_file,omics,trait):
    with open(status_file,newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    rows=[r for r in rows if r.get('omics')==omics and r.get('trait')==trait]
    if len(rows)!=1:raise ValueError('Expected exactly one fresh per-job status row')
    row=rows[0]
    if row['status'] not in ('ok','no_estimate'):raise ValueError('Scientific worker status: '+row['status'])
    if row['status']=='ok':
        path=Path(row['out_prefix'])
        if not path.is_file():raise ValueError('ok without a result file')
        with path.open(newline='') as f:
            reader=csv.DictReader(f,delimiter='\t');header=reader.fieldnames;records=list(reader)
        if not records or not {'alpha','se(alpha)','p(alpha)'}<=set(header or []):
            raise ValueError('Invalid MR-link-2 estimate output')
        for record in records:
            beta,se,p=map(float,(record['alpha'],record['se(alpha)'],record['p(alpha)']))
            if not all(map(math.isfinite,(beta,se,p))) or se<=0 or not 0<=p<=1:
                raise ValueError('Nonfinite/invalid MR-link-2 alpha inference')
        row['result_sha256']=digest_file(path)
    return row


def worker_main():
    p=argparse.ArgumentParser();p.add_argument('--task',type=Path,required=True);p.add_argument('--marker',type=Path,required=True)
    a=p.parse_args(sys.argv[2:]);task=json.loads(a.task.read_text());a.marker.unlink(missing_ok=True)
    # Each attempt MUST use a fresh private worker outdir. Do not reuse its old
    # .complete/status files: cache reuse is decided by the parent fingerprints.
    worker_out=Path(task['worker_outdir'])
    if worker_out.exists():raise ValueError('Worker outdir exists; create a new isolated attempt before retry')
    worker_out.mkdir(parents=True)
    proc=subprocess.run(task['argv'],cwd=task.get('cwd'),check=False)
    if proc.returncode:return proc.returncode
    row=validate_status(worker_out/'mrlink2.status.tsv',task['omics'],task['trait'])
    audits=[dict(path=str(p.resolve()),sha256=digest_file(p)) for p in sorted(worker_out.rglob('*_numerical_audit.json'))]
    if row['status']=='ok' and not audits:raise ValueError('Estimate missing ordered SNP/allele numerical audit')
    atomic_json(a.marker,{'numerical_audits':audits,'status':row['status'],'row':row,'task_id':task['task_id'],
                         'scientific_signature':task['scientific_signature']})
    return 0

if __name__=='__main__':
    try: sys.exit(worker_main() if len(sys.argv)>1 and sys.argv[1]=='worker' else match_main(sys.argv[2:]) if len(sys.argv)>1 and sys.argv[1]=='match' else prepare_main(sys.argv[2:]) if len(sys.argv)>1 and sys.argv[1]=='prepare' else main())
    except (OSError,ValueError,KeyError) as exc: print('MR scheduler: '+str(exc),file=sys.stderr);sys.exit(2)
