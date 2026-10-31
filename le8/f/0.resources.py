#!/usr/bin/env python3
"""Bounded external-process executor, Python stdlib only.

This is an execution layer, not a statistical estimator. A task is successful
only when exit status AND all declared output files pass validation. Global
admission leases coordinate cooperating invocations. RAM claims are estimates;
retain the existing whole-process-tree cgroup hard limit in LE8.
"""
from __future__ import annotations
import argparse, contextlib, csv, re, hashlib, json, math, os, signal, subprocess, sys, tempfile, time, uuid
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


def process_tree(pid):
    """Snapshot descendants by PID and birth time before stopping their parent."""
    entries={}
    for path in Path('/proc').glob('[0-9]*/stat'):
        try:
            raw=path.read_text();fields=raw[raw.rfind(')')+2:].split()
            entries[int(path.parent.name)]=(int(fields[1]),fields[19])
        except (OSError,ValueError,IndexError):pass
    selected={pid}
    while True:
        children={child for child,(parent,_) in entries.items() if parent in selected}
        if children<=selected:break
        selected|=children
    return {child:entries[child][1] for child in selected if child in entries}


def signal_tree(snapshot,signum):
    for pid,birth in snapshot.items():
        if process_identity(pid)==birth:
            try:os.kill(pid,signum)
            except ProcessLookupError:pass


def stop_worker(proc):
    # Workers share the terminal's process group: never signal that whole group
    # here, since it also contains the waiting supervisor and console logger.
    snapshot=process_tree(proc.pid)
    signal_tree(snapshot,signal.SIGTERM)
    try:proc.wait(timeout=5)
    except subprocess.TimeoutExpired:pass
    signal_tree(snapshot,signal.SIGKILL)
    proc.wait()


def memory_available_gib():
    try:
        fields=dict(line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())
        return int(fields['MemAvailable'].split()[0])/1024**2
    except (OSError,KeyError,ValueError):return float('inf')


def phase_reservation(label, hard_limit):
    # Measured C1 peaks: CPU ~10 GiB, GPU input preparation ~5 GiB.
    # Leave room for larger layers and training/evaluation allocations.
    estimate = {"c1_abm": 16.0, "c1_abm_GPU": 16.0, "c1_correlate": 20.0}.get(label, hard_limit)
    return min(float(hard_limit), estimate)


def process_tree_memory_gib(pid):
    """Use bounded cgroup memory (includes page cache), RSS fallback otherwise."""
    pending, seen, scopes, rss = [int(pid)], set(), {}, 0
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        root = Path('/proc') / str(current)
        try:
            pending.extend(map(int, (root / 'task' / str(current) / 'children').read_text().split()))
            for line in (root / 'status').read_text().splitlines():
                if line.startswith('VmRSS:'):
                    rss += int(line.split()[1]) * 1024
            for line in (root / 'cgroup').read_text().splitlines():
                if not line.startswith('0::'):
                    continue
                path = Path('/sys/fs/cgroup') / line[3:].lstrip('/')
                if (path / 'memory.max').read_text().strip() != 'max':
                    scopes[path] = int((path / 'memory.current').read_text())
        except (OSError, ValueError):
            continue
    # Avoid counting nested bounded cgroups twice.
    if scopes:
        return sum(v for p,v in scopes.items() if not any(q in p.parents for q in scopes)) / 2**30
    return rss / 2**30


class SharedBudget:
    def __init__(self,path,cores,memory_gib,reserve_gib=4,cpus=None):
        if cores<1 or memory_gib<=0 or reserve_gib<0:raise ValueError('Invalid resource capacity')
        self.path=Path(path);self.path.parent.mkdir(parents=True,exist_ok=True)
        self.lock=self.path.with_suffix(self.path.suffix+'.lock')
        self.cores=int(cores);self.memory=float(memory_gib);self.reserve=float(reserve_gib)
        inherited=','.join(map(str,cpus)) if cpus is not None else os.getenv('LE8_POOL_CPUS','')
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

    def acquire(self,cores,memory_gib,gpu=False,hard_limit_gib=None,label=None):
        if cores>self.cores or memory_gib>self.memory:raise ValueError('One task exceeds total resource capacity')
        with self.transaction() as d:
            used={c for v in d['leases'].values() for c in v['cpus']}
            free=[c for c in self.cpus if c not in used]
            usage = {key: (process_tree_memory_gib(v['pid'])
                           if v.get('label') in ('c1_abm','c1_abm_GPU','c1_correlate') else v['memory_gib'])
                     for key,v in d['leases'].items()}
            charged = sum(max(v['memory_gib'], usage[key]) for key,v in d['leases'].items())
            available=self.memory-charged
            # MemAvailable already excludes current allocations. Reserve only
            # the unconsumed part of existing estimates, avoiding double counting.
            headroom = sum(max(0, v['memory_gib']-usage[key]) for key,v in d['leases'].items())
            gpu_key=os.getenv('LE8_GPU_KEY',os.getenv('CUDA_VISIBLE_DEVICES','0') or '0')
            if gpu and any(v.get('gpu')==gpu_key for v in d['leases'].values()):return None
            if len(free)<cores or available+1e-9<memory_gib or memory_available_gib()<memory_gib+headroom+self.reserve:
                return None
            token=uuid.uuid4().hex
            d['leases'][token]={'pid':os.getpid(),'birth':process_identity(os.getpid()),
                                'cpus':free[:cores],'memory_gib':memory_gib,'gpu':gpu_key if gpu else None,
                                'hard_limit_gib':hard_limit_gib if hard_limit_gib is not None else memory_gib,'label':label}
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


def execute(tasks,outdir,cores=4,memory_gib=32,max_workers=4,pool_path=None,reserve_gib=4):
    outdir=Path(outdir).resolve();outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'ALL_COMPLETE.json').unlink(missing_ok=True)
    ids=[str(t['task_id']) for t in tasks]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate task IDs')
    if max_workers<1:raise ValueError('max_workers must be positive')
    pool=SharedBudget(pool_path or f'/tmp/le8-resources-{os.getuid()}.json',cores,memory_gib,reserve_gib)
    pending=[];result={};running={}
    try:
        previous={str(record['task_id']):record for record in json.loads((outdir/'run_status.json').read_text())}
    except (OSError, ValueError, KeyError, TypeError):previous={}
    for position,task in enumerate(tasks):
        task=dict(task);task['_position']=position
        if not isinstance(task['argv'],list) or not all(isinstance(s,str) for s in task['argv']) or not task['argv']:
            raise ValueError('argv must be a nonempty string list; shell strings are not accepted')
        task['_signature']=task_signature(task)
        task['_cache']=outdir/(hashlib.sha256(str(task['task_id']).encode()).hexdigest()+'.complete.json')
        try:
            old=previous.get(str(task['task_id']))
            if old is None:old=json.loads(task['_cache'].read_text())
            if old.get('status')=='completed' and old['signature']==task['_signature'] and old['outputs']==validate_outputs(task):
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
    for sig in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):old_handlers[sig]=signal.signal(sig,interrupted)
    try:
        while pending or running:
            launched=False
            for task in list(pending):
                if len(running)>=max_workers:break
                threads=int(task.get('cores',1));mem=float(task.get('memory_gib',1))
                acquired=pool.acquire(threads,mem)
                if acquired is None:continue
                token,cpus=acquired
                label=re.sub(r'[^A-Za-z0-9_.-]+','_',str(task['task_id']))[:64] or 'task'
                logpath=outdir/f"{task['_position']+1:04d}-{label}.log"
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
                                          preexec_fn=affinity)
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
                    rec['outputs']=validate_outputs(task);rec['status']='completed'
                except (OSError,ValueError) as exc:rec['reason']=str(exc)
                result[task['_position']]=rec
                atomic_json(outdir/'run_status.json',[result[i] for i in sorted(result)])
            if pending or running:time.sleep(.03 if launched else .1)
    finally:
        # Ctrl+C reaches both supervisors and workers. Ignore a second signal
        # from the console's cancellation handler until leases are released.
        for sig in old_handlers:signal.signal(sig,signal.SIG_IGN)
        for proc,task,token,handle,start,cpus in running.values():
            stop_worker(proc)
            handle.close();pool.release(token)
        for sig,old in old_handlers.items():signal.signal(sig,old)
    ordered=[result[i] for i in range(len(tasks))]
    atomic_json(outdir/'run_status.json',ordered)
    complete=outdir/'ALL_COMPLETE.json'
    if all(v['status']=='completed' for v in ordered):atomic_json(complete,{'tasks':[v['signature'] for v in ordered]})
    else:complete.unlink(missing_ok=True)
    return ordered



class BlockedKNN:
    """Compatibility constructor; all numerical code lives in ../0f/abm.py."""
    def __new__(cls, *args, **kwargs):
        # Supervisors must stay stdlib-only until retrieval is actually requested.
        import importlib.util
        if "shared_abm" not in sys.modules:
            path = Path(__file__).resolve().parents[2] / "0f" / "abm.py"
            spec = importlib.util.spec_from_file_location("shared_abm", path)
            module = importlib.util.module_from_spec(spec)
            sys.modules["shared_abm"] = module
            spec.loader.exec_module(module)
        return sys.modules["shared_abm"].BlockedKNN(*args, **kwargs)


def phase_parent_rss_gib(supervisor_pid):
    """Charge the waiting R parent and launchers, stopping at the phase owner."""
    pid=os.getpid();seen=set();rss_kib=0
    while pid>0 and pid not in seen:
        seen.add(pid)
        fields=dict(line.split(':',1) for line in Path(f'/proc/{pid}/status').read_text().splitlines())
        rss_kib+=int(fields.get('VmRSS','0 kB').split()[0])
        if pid==supervisor_pid:return rss_kib/1024**2
        pid=int(fields['PPid'].strip())
    raise ValueError('Phase supervisor is not an ancestor; cannot account for waiting parent RAM')


@contextlib.contextmanager
def phase_worker_budget(cores,memory_gib,pool_path,reserve_gib):
    """Subdivide an admitted phase for synchronous MR jobs without double admission.

    The outer lease remains in the global pool, protecting all of these workers
    from other runs. The waiting parent RSS and headroom are charged locally.
    """
    if os.getenv('LE8_PHASE_ADMITTED')!='1':
        yield cores,memory_gib,pool_path
        return
    parent=phase_parent_rss_gib(int(os.environ['LE8_PHASE_SUPERVISOR_PID']))
    budget=min(memory_gib,float(os.environ['LE8_PHASE_MEMORY_GIB']))
    remaining=budget-parent-reserve_gib
    if remaining<=0:
        raise ValueError(f'No RAM for MR workers: phase={budget:.3f} GiB, parent RSS={parent:.3f} GiB, reserve={reserve_gib:.3f} GiB')
    cpus=sorted(os.sched_getaffinity(0))[:min(cores,int(os.environ['LE8_PHASE_CORES']))]
    if not cpus:raise ValueError('No CPUs within the admitted phase')
    # The parent is blocked until execute() finishes, including child cleanup.
    names=('LE8_POOL_CPUS','LE8_GPU_LEASE_HELD','LE8_GPU_LOCAL_LEASE')
    previous={name:os.environ.get(name) for name in names}
    with tempfile.TemporaryDirectory(prefix='le8-phase-workers-',dir='/tmp') as directory:
        os.environ['LE8_POOL_CPUS']=','.join(map(str,cpus))
        if os.getenv('LE8_GPU_LEASE_HELD')=='1':
            # A GPU phase owns the device globally, but parallel MR children
            # still serialize their GPU operations in the local worker pool.
            os.environ['LE8_GPU_LOCAL_LEASE']='1'
            os.environ['LE8_GPU_LEASE_HELD']='0'
        try:
            print(f'[LE8] RESOURCE MR: cores={len(cpus)} parent_RSS_GiB={parent:.3f} worker_budget_GiB={remaining:.3f} phase_GiB={budget:.3f}',flush=True)
            yield len(cpus),remaining,str(Path(directory)/'pool.json')
        finally:
            for name,value in previous.items():
                if value is None:os.environ.pop(name,None)
                else:os.environ[name]=value


@contextlib.contextmanager
def gpu_lease():
    """GPU-only admission inside an already admitted MR CPU job."""
    if os.getenv('LE8_GPU_LEASE_HELD')=='1':
        yield
        return
    # MR workers may use a local CPU/RAM pool; GPU exclusion stays global.
    prefix='LE8_' if os.getenv('LE8_GPU_LOCAL_LEASE')=='1' else 'LE8_GLOBAL_'
    gpu_cpus=os.getenv(prefix+'POOL_CPUS')
    pool=SharedBudget(os.getenv(prefix+'RESOURCE_POOL',os.getenv('LE8_RESOURCE_POOL',f'/tmp/le8-resources-{os.getuid()}.json')),
        int(os.getenv(prefix+'TOTAL_CORES',os.getenv('LE8_TOTAL_CORES','16'))),
        float(os.getenv(prefix+'RESOURCE_MEMORY_GIB',os.getenv('LE8_RESOURCE_MEMORY_GIB','32'))),0,
        cpus=None if gpu_cpus is None else list(map(int,gpu_cpus.split(','))))
    slot=None
    while slot is None:
        slot=pool.acquire(0,0,gpu=True)
        if slot is None:time.sleep(.2)
    previous=os.environ.get('LE8_GPU_LEASE_HELD')
    os.environ['LE8_GPU_LEASE_HELD']='1'
    try:yield
    finally:
        pool.release(slot[0])
        if previous is None:os.environ.pop('LE8_GPU_LEASE_HELD',None)
        else:os.environ['LE8_GPU_LEASE_HELD']=previous


def phase_main(argv):
    p=argparse.ArgumentParser(description='Shared phase-level CPU/RAM/GPU admission; cgroup remains the hard cap')
    p.add_argument('--cores',type=int,required=True);p.add_argument('--memory-gib',type=float,required=True)
    p.add_argument('--gpu',action='store_true');p.add_argument('--label',default='phase')
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args(argv)
    command=a.command[1:] if a.command[:1]==['--'] else a.command
    if not command:raise ValueError('Missing phase command')
    if a.cores<1 or a.memory_gib<=0:raise ValueError('Invalid phase request')
    pool=SharedBudget(os.getenv('LE8_RESOURCE_POOL',f'/tmp/le8-resources-{os.getuid()}.json'),
        int(os.getenv('LE8_TOTAL_CORES','16')),float(os.getenv('LE8_RESOURCE_MEMORY_GIB','32')),
        float(os.getenv('LE8_RESOURCE_RESERVE_GIB','2')))
    reservation=phase_reservation(a.label,a.memory_gib)
    # Existing leases already reserve their unconsumed RAM in acquire(). Keep
    # 4 GiB for the rest of the system in addition to those full phase claims.
    # Two 8-core/16-GiB C1 runs can then share the 16-core/32-GiB pool without
    # relaxing either process tree's cgroup limit. Explicit larger reserves win.
    if a.label in ('c1_abm','c1_abm_GPU','c1_correlate'):pool.reserve=max(pool.reserve,4.0)
    slot=None;started=time.monotonic();proc=None
    while slot is None:
        slot=pool.acquire(a.cores,reservation,gpu=a.gpu,hard_limit_gib=a.memory_gib,label=a.label)
        if slot is None:
            if time.monotonic()-started<.5:print('[LE8] Waiting for shared resources: '+a.label+' (retry every 10 seconds)',flush=True)
            time.sleep(10)
    token,cpus=slot
    env=dict(os.environ,LE8_POOL_CPUS=','.join(map(str,pool.cpus)),LE8_PHASE_MEMORY_GIB=str(a.memory_gib),
        LE8_PHASE_CORES=str(a.cores),LE8_PHASE_ADMITTED='1',LE8_TOTAL_CORES=str(pool.cores),
        LE8_RESOURCE_MEMORY_GIB=str(pool.memory),LE8_RESOURCE_POOL=str(pool.path),
        LE8_PHASE_SUPERVISOR_PID=str(os.getpid()),LE8_GLOBAL_RESOURCE_POOL=str(pool.path),
        LE8_GLOBAL_TOTAL_CORES=str(pool.cores),LE8_GLOBAL_RESOURCE_MEMORY_GIB=str(pool.memory),
        LE8_GLOBAL_POOL_CPUS=','.join(map(str,pool.cpus)))
    if a.gpu:env['LE8_GPU_LEASE_HELD']='1'
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    print('[LE8] Resource admission '+json.dumps(dict(phase=a.label,cpus=cpus,memory_claim_gib=reservation,memory_limit_gib=a.memory_gib,system_reserve_gib=pool.reserve,
        exclusive_gpu=a.gpu,pool=str(pool.path),wait_seconds=round(time.monotonic()-started,3))),flush=True)
    old={}
    def stop(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    try:
        for sig in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):old[sig]=signal.signal(sig,stop)
        def affinity():os.sched_setaffinity(0,cpus)
        proc=subprocess.Popen(command,env=env,preexec_fn=affinity)
        # Share terminal job control; this supervisor cleans its own descendants.
        return proc.wait()
    finally:
        for sig in old:signal.signal(sig,signal.SIG_IGN)
        if proc is not None and proc.poll() is None:
            stop_worker(proc)
        pool.release(token)
        for sig,handler in old.items():signal.signal(sig,handler)


def resolve_eigh_backend(request='auto'):
    if request not in ('auto','cpu','cuda'):
        raise ValueError('eigh backend must be auto|cpu|cuda')
    if request=='cpu':return 'cpu'
    try:
        import torch
        available=torch.cuda.is_available()
    except ImportError:
        available=False
    if request=='cuda' and not available:
        raise RuntimeError('MR-link-2 CUDA eigendecomposition requested but unavailable')
    return 'cuda' if available else 'cpu'


def validated_eigh(matrix,backend='cpu',probes=8,tolerance=1e-8):
    import numpy as np
    backend=resolve_eigh_backend(backend)
    started=time.perf_counter();gpu_name=None;peak=0
    a=np.asarray(matrix,dtype=np.float64)
    if a.ndim!=2 or a.shape[0]!=a.shape[1] or not len(a) or not np.isfinite(a).all():raise ValueError('Finite square matrix required')
    if not np.allclose(a,a.T,rtol=0,atol=1e-10):raise ValueError('LD is not symmetric; do not silently alter it')
    if backend=='cpu':values,vectors=np.linalg.eigh(a)
    elif backend=='cuda':
        import torch
        if not torch.cuda.is_available():raise RuntimeError('CUDA requested; CPU fallback is not silent')
        torch.cuda.reset_peak_memory_stats()
        gpu_name=torch.cuda.get_device_name()
        x=torch.as_tensor(a,device='cuda',dtype=torch.float64)
        lam,u=torch.linalg.eigh(x);torch.cuda.synchronize()
        values,vectors=lam.cpu().numpy(),u.cpu().numpy()
        peak=torch.cuda.max_memory_allocated()/2**20
    else:raise ValueError('backend must be cpu or cuda')
    # Random deterministic probes avoid another full O(m^3) residual evaluation.
    r=np.random.default_rng(581).normal(size=(len(a),min(probes,len(a))))
    residual=np.linalg.norm(a@(vectors@r)-vectors@(values[:,None]*r))/max(np.linalg.norm(a)*np.linalg.norm(r),1e-30)
    orthogonality=np.linalg.norm(vectors.T@(vectors@r)-r)/max(np.linalg.norm(r),1e-30)
    trace_error=abs(values.sum()-np.trace(a))/max(abs(np.trace(a)),1e-30)
    if not np.isfinite([residual,orthogonality,trace_error]).all() or max(residual,orthogonality,trace_error)>tolerance:
        raise ValueError('Eigendecomposition validation failed; keep failure explicit')
    return values,vectors,dict(backend=backend,dtype='float64',wall_seconds=time.perf_counter()-started,
                               gpu_name=gpu_name,gpu_peak_memory_mb=peak,n=len(a),relative_residual=float(residual),
                               orthogonality_probe=float(orthogonality),relative_trace_error=float(trace_error))

def cached_eigh(matrix, ordered_alleles, backend='cpu'):
    """Reuse only byte-identical float64 LD with the same ordered SNP/alleles."""
    import numpy as np
    import contextlib
    backend=resolve_eigh_backend(backend)
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
