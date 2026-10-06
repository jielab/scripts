#!/usr/bin/env python3
"""Bounded external-process executor, Python stdlib only.

This is an execution layer, not a statistical estimator. A task is successful
only when exit status AND all declared output files pass validation. Global
admission leases coordinate cooperating invocations. RAM claims are estimates;
retain the existing whole-process-tree cgroup hard limit in LE8.
"""
from __future__ import annotations
import argparse, contextlib, csv, hashlib, json, math, os, signal, subprocess, sys, tempfile, time, uuid
from pathlib import Path
import fcntl


def digest_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def atomic_json(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name+'-',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f:
            json.dump(value,f,indent=2,allow_nan=False);f.flush();os.fsync(f.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)


def process_identity(pid):
    try:
        raw=Path(f'/proc/{pid}/stat').read_text();tail=raw[raw.rfind(')')+2:].split()
        if tail[0]=='Z':return None
        return tail[19]  # field 22: starttime; field 3 is tail[0]
    except (OSError,IndexError):return None


def memory_available_gib():
    try:
        fields=dict(line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())
        return int(fields['MemAvailable'].split()[0])/1024**2
    except (OSError,KeyError,ValueError):return float('inf')


class SharedBudget:
    def __init__(self,path,cores,memory_gib,reserve_gib=4):
        if cores<1 or memory_gib<=0 or reserve_gib<0:raise ValueError('Invalid resource capacity')
        self.path=Path(path);self.path.parent.mkdir(parents=True,exist_ok=True)
        self.lock=self.path.with_suffix(self.path.suffix+'.lock')
        self.cores=int(cores);self.memory=float(memory_gib);self.reserve=float(reserve_gib)
        inherited=os.getenv('LE8_POOL_CPUS','')
        self.cpus=([int(v) for v in inherited.split(',')] if inherited else sorted(os.sched_getaffinity(0)))[:self.cores] if hasattr(os,'sched_getaffinity') else list(range(self.cores))
        self.cores=len(self.cpus)
        self.gpu_jobs=int(os.getenv('LE8_GPU_JOBS','1'))
        if self.gpu_jobs!=1:raise ValueError('This release requires --gpu-jobs 1 for exclusive GPU phases')

    @contextlib.contextmanager
    def transaction(self):
        with open(self.lock,'a') as handle:
            fcntl.flock(handle,fcntl.LOCK_EX)
            try:
                value=json.loads(self.path.read_text()) if self.path.exists() else {'leases':{}}
                value['leases']={k:v for k,v in value['leases'].items() if process_identity(v['pid']) is not None and process_identity(v['pid'])==v['birth']}
                active=bool(value['leases'])
                definition={'cpus':self.cpus,'memory_gib':self.memory,'gpu_jobs':self.gpu_jobs}
                if active and value.get('definition')!=definition:
                    raise ValueError('Other active LE8 jobs use different global CPU/RAM capacity; use the same budget settings')
                value['definition']=definition
                yield value
                atomic_json(self.path,value)
            finally:fcntl.flock(handle,fcntl.LOCK_UN)

    def acquire(self,cores,memory_gib,gpu=False):
        if cores>self.cores or memory_gib>self.memory:raise ValueError('One task exceeds total resource capacity')
        with self.transaction() as d:
            used={c for v in d['leases'].values() for c in v['cpus']}
            free=[c for c in self.cpus if c not in used]
            available=self.memory-sum(v['memory_gib'] for v in d['leases'].values())
            gpu_key=os.getenv('LE8_GPU_KEY',os.getenv('CUDA_VISIBLE_DEVICES','0') or '0')
            if gpu and any(v.get('gpu')==gpu_key for v in d['leases'].values()):return None
            if len(free)<cores or available+1e-9<memory_gib or memory_available_gib()<memory_gib+self.reserve:
                return None
            token=uuid.uuid4().hex
            d['leases'][token]={'pid':os.getpid(),'birth':process_identity(os.getpid()),
                                'cpus':free[:cores],'memory_gib':memory_gib,'gpu':gpu_key if gpu else None}
            return token,free[:cores]

    def transfer(self,token,pid):
        with self.transaction() as d:
            if token in d['leases']:
                d['leases'][token].update(pid=pid,birth=process_identity(pid))

    def release(self,token):
        with self.transaction() as d:d['leases'].pop(token,None)


def validate_outputs(task):
    result=[]
    for item in task.get('outputs',[]):
        cfg={'path':item} if isinstance(item,str) else item
        path=Path(cfg['path'])
        if not path.is_file() or path.stat().st_size<cfg.get('min_bytes',1):raise ValueError('Missing/empty output: '+str(path))
        if cfg.get('kind') in ('json','mrlink_marker'):
            record=json.loads(path.read_text())
            if cfg.get('kind')=='mrlink_marker':
                if record.get('status') not in ('ok','no_estimate'):raise ValueError('Invalid MR marker')
                for audit in record.get('numerical_audits',[]):
                    if digest_file(audit['path'])!=audit['sha256']:raise ValueError('MR alignment audit changed')
                if record['status']=='ok':
                    row=record['row'];target=Path(row['out_prefix'])
                    if not target.is_file() or digest_file(target)!=row['result_sha256']:
                        raise ValueError('MR result changed after completion')
        result.append({'path':str(path.resolve()),'sha256':digest_file(path)})
    if not result:raise ValueError('Declare outputs; exit 0 alone is not completion')
    return result


def task_signature(task):
    # Caller supplies verified scientific/data identities. Execution order and
    # requested worker count deliberately do not enter the scientific identity.
    value={k:task[k] for k in ('task_id','argv','cwd','scientific_signature','outputs') if k in task}
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()


def execute(tasks,outdir,cores=4,memory_gib=16,max_workers=4,pool_path=None,reserve_gib=4):
    outdir=Path(outdir).resolve();outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'ALL_COMPLETE.json').unlink(missing_ok=True)
    ids=[str(t['task_id']) for t in tasks]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate task IDs')
    if max_workers<1:raise ValueError('max_workers must be positive')
    pool=SharedBudget(pool_path or f'/tmp/le8-resources-{os.getuid()}.json',cores,memory_gib,reserve_gib)
    pending=[];result={};running={}
    for position,task in enumerate(tasks):
        task=dict(task);task['_position']=position
        if not isinstance(task['argv'],list) or not all(isinstance(s,str) for s in task['argv']) or not task['argv']:
            raise ValueError('argv must be a nonempty string list; shell strings are not accepted')
        task['_signature']=task_signature(task)
        task['_cache']=outdir/(hashlib.sha256(str(task['task_id']).encode()).hexdigest()+'.complete.json')
        try:
            old=json.loads(task['_cache'].read_text())
            if old['signature']==task['_signature'] and old['outputs']==validate_outputs(task):
                result[position]={**old,'cache_hit':True};continue
        except (OSError,ValueError,KeyError):pass
        task['_cache'].unlink(missing_ok=True)
        if int(task.get('cores',1))<1 or float(task.get('memory_gib',1))<=0:raise ValueError('Invalid task resources')
        if int(task.get('cores',1))>pool.cores or float(task.get('memory_gib',1))>pool.memory:raise ValueError('Task exceeds resource budget')
        pending.append(task)
    # Largest tasks first avoids leaving one big LD job as the final straggler.
    pending.sort(key=lambda t:(-float(t.get('memory_gib',1)),t['_position']))
    old_handlers={}
    def interrupted(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    for sig in (signal.SIGTERM,signal.SIGINT):old_handlers[sig]=signal.signal(sig,interrupted)
    try:
        while pending or running:
            launched=False
            for task in list(pending):
                if len(running)>=max_workers:break
                threads=int(task.get('cores',1));mem=float(task.get('memory_gib',1))
                acquired=pool.acquire(threads,mem)
                if acquired is None:continue
                token,cpus=acquired
                logpath=outdir/(hashlib.sha256(str(task['task_id']).encode()).hexdigest()+'.log')
                handle=open(logpath,'w')
                env=os.environ.copy();env.update({k:str(v) for k,v in task.get('env',{}).items()})
                for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
                    env[key]=str(threads)
                env['LE8_IN_WORKER']='1';env['LE8_INNER_THREADS']=str(threads)
                env['LE8_POOL_CPUS']=','.join(map(str,pool.cpus))
                env['LE8_TOTAL_CORES']=str(pool.cores)
                env['LE8_RESOURCE_MEMORY_GIB']=str(pool.memory)
                env['LE8_RESOURCE_POOL']=str(pool.path)
                def affinity():
                    if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,cpus)
                try:
                    proc=subprocess.Popen(task['argv'],cwd=task.get('cwd'),env=env,stdout=handle,stderr=subprocess.STDOUT,
                                          start_new_session=True,preexec_fn=affinity)
                    pool.transfer(token,proc.pid)
                except BaseException:
                    handle.close();pool.release(token);raise
                running[proc.pid]=(proc,task,token,handle,time.monotonic(),cpus)
                pending.remove(task);launched=True
            for pid,(proc,task,token,handle,start,cpus) in list(running.items()):
                code=proc.poll()
                if code is None:continue
                handle.close();pool.release(token);running.pop(pid)
                rec={'task_id':task['task_id'],'signature':task['_signature'],'exit_code':code,
                     'wall_seconds':time.monotonic()-start,'cores':cpus,'cache_hit':False,'status':'failed'}
                rec['memory_claim_gib']=float(task.get('memory_gib',1))
                try:
                    if code:raise ValueError('Worker exited with '+str(code))
                    rec['outputs']=validate_outputs(task);rec['status']='completed';atomic_json(task['_cache'],rec)
                except (OSError,ValueError) as exc:rec['reason']=str(exc)
                result[task['_position']]=rec
                atomic_json(outdir/'progress.json',{'finished':len(result),'total':len(tasks),'failed':sum(v['status']!='completed' for v in result.values())})
            if pending or running:time.sleep(.03 if launched else .1)
    finally:
        for proc,task,token,handle,start,cpus in running.values():
            try:os.killpg(proc.pid,signal.SIGTERM)
            except ProcessLookupError:pass
        for proc,task,token,handle,start,cpus in running.values():
            try:proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:os.killpg(proc.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                proc.wait()
            handle.close();pool.release(token)
        for sig,old in old_handlers.items():signal.signal(sig,old)
    ordered=[result[i] for i in range(len(tasks))]
    atomic_json(outdir/'run_status.json',ordered)
    complete=outdir/'ALL_COMPLETE.json'
    if all(v['status']=='completed' for v in ordered):atomic_json(complete,{'tasks':[v['signature'] for v in ordered]})
    else:complete.unlink(missing_ok=True)
    return ordered



class BlockedKNN:
    def __init__(self,coordinates,ids,groups,tie_keys,device='cpu',dtype='float64',query_block=32,donor_block=2048):
        global np, torch
        import numpy as np
        import torch
        x=np.asarray(coordinates)
        if x.ndim!=2 or not len(x) or not np.isfinite(x).all():raise ValueError('Finite donor matrix required')
        self.ids=np.asarray(ids,str);self.groups=np.asarray(groups,str)
        if len(set(self.ids))!=len(x) or len(self.groups)!=len(x) or len(tie_keys)!=len(x):raise ValueError('Invalid donor identities')
        if query_block<1 or donor_block<1:raise ValueError('Positive block sizes required')
        if device.startswith('cuda') and not torch.cuda.is_available():raise RuntimeError('CUDA explicitly requested, no CPU fallback')
        self.device=device;self.dtype={'float64':torch.float64,'float32':torch.float32}[dtype]
        self.order=np.argsort(np.asarray(tie_keys),kind='stable')
        self.x=torch.as_tensor(x[self.order].copy(),dtype=self.dtype,device=device)
        self.sid=self.ids[self.order];self.sgroup=self.groups[self.order]
        self.query_block=query_block;self.donor_block=donor_block

    def query(self,coordinates,ids,groups,k,query_folds=None,bank_folds=None):
        q=np.ascontiguousarray(coordinates);ids=np.asarray(ids,str);groups=np.asarray(groups,str)
        if q.ndim!=2 or q.shape[1]!=self.x.shape[1] or len(ids)!=len(q) or len(groups)!=len(q) or not np.isfinite(q).all():raise ValueError('Invalid query')
        if k<1:raise ValueError('Positive k required')
        if (query_folds is None)!=(bank_folds is None):raise ValueError('Supply both query and bank folds')
        bf=None if bank_folds is None else np.asarray(bank_folds)[self.order]
        if bf is not None and (len(bf)!=len(self.x) or len(query_folds)!=len(q)):raise ValueError('Fold shape mismatch')
        keep=min(k,len(self.x));all_index=[];all_distance=[]
        for begin in range(0,len(q),self.query_block):
            end=min(begin+self.query_block,len(q));xx=torch.as_tensor(q[begin:end],dtype=self.dtype,device=self.device)
            best_d=torch.empty((len(xx),0),dtype=self.dtype,device=self.device)
            best_i=torch.empty((len(xx),0),dtype=torch.int64,device=self.device)
            for start in range(0,len(self.x),self.donor_block):
                stop=min(start+self.donor_block,len(self.x))
                # Direct Euclidean accumulation avoids norm-subtraction cancellation.
                d=torch.cdist(xx,self.x[start:stop],p=2,compute_mode='donot_use_mm_for_euclid_dist')
                valid=(ids[begin:end,None]!=self.sid[None,start:stop])&(groups[begin:end,None]!=self.sgroup[None,start:stop])
                if bf is not None:valid &= np.asarray(query_folds)[begin:end,None]!=bf[None,start:stop]
                d.masked_fill_(~torch.as_tensor(valid,device=self.device),torch.inf)
                ix=torch.arange(start,stop,device=self.device).expand(len(xx),-1)
                di=torch.cat([best_d,d],1);ii=torch.cat([best_i,ix],1)
                # Stable lexicographic (distance, ID rank), including block boundaries.
                rank=torch.argsort(ii,dim=1,stable=True);di=di.gather(1,rank);ii=ii.gather(1,rank)
                rank=torch.argsort(torch.round(di/1e-10),dim=1,stable=True)[:,:keep]
                best_d=di.gather(1,rank);best_i=ii.gather(1,rank)
            distance=best_d.cpu().numpy();index=self.order[best_i.cpu().numpy()]
            index[~np.isfinite(distance)]=-1
            all_index.append(index);all_distance.append(distance)
        if not len(q):return np.empty((0,keep),int),np.empty((0,keep),float)
        return np.concatenate(all_index),np.concatenate(all_distance)


@contextlib.contextmanager
def gpu_lease():
    """GPU-only admission inside an already admitted MR CPU job."""
    if os.getenv('LE8_GPU_LEASE_HELD')=='1':
        yield
        return
    pool=SharedBudget(os.getenv('LE8_RESOURCE_POOL',f'/tmp/le8-resources-{os.getuid()}.json'),
        int(os.getenv('LE8_TOTAL_CORES','16')),float(os.getenv('LE8_RESOURCE_MEMORY_GIB','24')),0)
    slot=None
    while slot is None:
        slot=pool.acquire(0,0,gpu=True)
        if slot is None:time.sleep(.2)
    try:yield
    finally:pool.release(slot[0])


def phase_main(argv):
    p=argparse.ArgumentParser(description='Shared phase-level CPU/RAM/GPU admission; cgroup remains the hard cap')
    p.add_argument('--cores',type=int,required=True);p.add_argument('--memory-gib',type=float,required=True)
    p.add_argument('--gpu',action='store_true');p.add_argument('--label',default='phase')
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args(argv)
    command=a.command[1:] if a.command[:1]==['--'] else a.command
    if not command:raise ValueError('Missing phase command')
    if a.cores<1 or a.memory_gib<=0:raise ValueError('Invalid phase request')
    pool=SharedBudget(os.getenv('LE8_RESOURCE_POOL',f'/tmp/le8-resources-{os.getuid()}.json'),
        int(os.getenv('LE8_TOTAL_CORES','16')),float(os.getenv('LE8_RESOURCE_MEMORY_GIB','24')),
        float(os.getenv('LE8_RESOURCE_RESERVE_GIB','2')))
    slot=None;started=time.monotonic();proc=None
    while slot is None:
        slot=pool.acquire(a.cores,a.memory_gib,gpu=a.gpu)
        if slot is None:
            if time.monotonic()-started<.5:print('[LE8] Waiting for shared resources: '+a.label+' (retry every 10 minutes)',flush=True)
            time.sleep(600)
    token,cpus=slot
    env=dict(os.environ,LE8_POOL_CPUS=','.join(map(str,pool.cpus)),LE8_PHASE_MEMORY_GIB=str(a.memory_gib),
        LE8_PHASE_CORES=str(a.cores),LE8_PHASE_ADMITTED='1',LE8_TOTAL_CORES=str(pool.cores),
        LE8_RESOURCE_MEMORY_GIB=str(pool.memory),LE8_RESOURCE_POOL=str(pool.path))
    if a.gpu:env['LE8_GPU_LEASE_HELD']='1'
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    print('[LE8] Resource admission '+json.dumps(dict(phase=a.label,cpus=cpus,memory_claim_gib=a.memory_gib,
        exclusive_gpu=a.gpu,pool=str(pool.path),wait_seconds=round(time.monotonic()-started,3))),flush=True)
    old={}
    def stop(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    try:
        for sig in (signal.SIGTERM,signal.SIGINT):old[sig]=signal.signal(sig,stop)
        def affinity():os.sched_setaffinity(0,cpus)
        proc=subprocess.Popen(command,env=env,start_new_session=True,preexec_fn=affinity)
        # Keep ownership in this supervisor, which cleans the complete child group.
        return proc.wait()
    finally:
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid,signal.SIGTERM)
            try:proc.wait(timeout=5)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
        pool.release(token)
        for sig,handler in old.items():signal.signal(sig,handler)


def validated_eigh(matrix,backend='cpu',probes=8,tolerance=1e-8):
    import numpy as np
    a=np.asarray(matrix,dtype=np.float64)
    if a.ndim!=2 or a.shape[0]!=a.shape[1] or not len(a) or not np.isfinite(a).all():raise ValueError('Finite square matrix required')
    if not np.allclose(a,a.T,rtol=0,atol=1e-10):raise ValueError('LD is not symmetric; do not silently alter it')
    if backend=='cpu':values,vectors=np.linalg.eigh(a)
    elif backend=='cuda':
        import torch
        if not torch.cuda.is_available():raise RuntimeError('CUDA requested; CPU fallback is not silent')
        x=torch.as_tensor(a,device='cuda',dtype=torch.float64)
        lam,u=torch.linalg.eigh(x);torch.cuda.synchronize()
        values,vectors=lam.cpu().numpy(),u.cpu().numpy()
    else:raise ValueError('backend must be cpu or cuda')
    # Random deterministic probes avoid another full O(m^3) residual evaluation.
    r=np.random.default_rng(581).normal(size=(len(a),min(probes,len(a))))
    residual=np.linalg.norm(a@(vectors@r)-vectors@(values[:,None]*r))/max(np.linalg.norm(a)*np.linalg.norm(r),1e-30)
    orthogonality=np.linalg.norm(vectors.T@(vectors@r)-r)/max(np.linalg.norm(r),1e-30)
    trace_error=abs(values.sum()-np.trace(a))/max(abs(np.trace(a)),1e-30)
    if not np.isfinite([residual,orthogonality,trace_error]).all() or max(residual,orthogonality,trace_error)>tolerance:
        raise ValueError('Eigendecomposition validation failed; keep failure explicit')
    return values,vectors,dict(backend=backend,dtype='float64',relative_residual=float(residual),
                               orthogonality_probe=float(orthogonality),relative_trace_error=float(trace_error))

def cached_eigh(matrix, ordered_alleles, backend='cpu'):
    """Reuse only byte-identical float64 LD with the same ordered SNP/alleles."""
    import numpy as np
    import contextlib
    a=np.ascontiguousarray(matrix,dtype=np.float64)
    identity=dict(ordered_alleles=ordered_alleles,ld_sha256=hashlib.sha256(a.tobytes()).hexdigest(),
        backend=backend,dtype='float64',shape=list(a.shape),algorithm='validated-eigh-v1',numpy=np.__version__)
    if backend=='cuda':
        import torch
        identity['torch']=torch.__version__;identity['cuda']=torch.version.cuda
    key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    root=Path(os.getenv('MRLINK2_EIGH_CACHE','/tmp/le8-mrlink2-eigh'))/key
    root.mkdir(parents=True,exist_ok=True)
    with open(root/'lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        target=root/'eigen.npz';proof=root/'complete.json'
        try:
            saved=json.loads(proof.read_text())
            if saved['identity']!=identity or saved['sha256']!=digest_file(target):raise ValueError('Invalid eigen cache')
            with np.load(target,allow_pickle=False) as z:values=z['values'];vectors=z['vectors']
            return values,vectors,{**saved['audit'],'cache_hit':True,'cache_key':key}
        except (OSError,ValueError,KeyError):pass
        with gpu_lease() if backend=='cuda' else contextlib.nullcontext():values,vectors,audit=validated_eigh(a,backend)
        temporary=root/('build-'+uuid.uuid4().hex+'.npz')
        try:
            np.savez(temporary,values=values,vectors=vectors)
            os.replace(temporary,target)
            atomic_json(proof,dict(identity=identity,audit=audit,sha256=digest_file(target)))
        finally:temporary.unlink(missing_ok=True)
        return values,vectors,{**audit,'cache_hit':False,'cache_key':key}

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='phase':sys.exit(phase_main(sys.argv[2:]))
    raise SystemExit('Use: 0.resources.py phase --cores N --memory-gib N [--gpu] -- command ...')
