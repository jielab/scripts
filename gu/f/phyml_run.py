#!/usr/bin/env python3
"""Run one PhyML tree; a pre-bootstrap main tree is never a completion marker."""
import argparse, fcntl, hashlib, json, os, re, shutil, signal, subprocess, sys, time
from pathlib import Path

SUFFIXES=['_phyml_tree.txt','_phyml_stats.txt','_phyml_boot_trees.txt','_phyml_boot_stats.txt',
          '_phyml_tree.png','_phyml_tree.pdf','.phyml.log','.phyml.complete.json']

def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def outputs(phy,boot):
    names=['_phyml_tree.txt','_phyml_stats.txt','.phyml.log']
    if boot>0:names+=['_phyml_boot_trees.txt','_phyml_boot_stats.txt']
    return [Path(str(phy)+s) for s in names]

def completion_error(phy,boot):
    files=outputs(phy,boot)
    if any(not p.is_file() or not p.stat().st_size for p in files):return 'missing/empty tree, stats, log or bootstrap outputs'
    log=Path(str(phy)+'.phyml.log').read_text(errors='replace')
    finish=list(re.finditer(r'\. Time used\s+\d+h\d+m\d+s',log))
    if not finish:return 'no final runtime footer (interrupted run)'
    if not re.search(r'Printing the most likely tree',log):return 'missing final tree report'
    if boot>0:
        progress=list(re.finditer(r'\b'+str(boot)+r'\s*/\s*'+str(boot)+r'\b',log))
        if not progress or progress[-1].end()>finish[-1].start():return f'bootstrap did not finish {boot}/{boot}'
        with Path(str(phy)+'_phyml_boot_trees.txt').open() as f:
            count=0
            for line in f:
                line=line.strip()
                if not line:continue
                if not line.endswith(';') or line.count('(')!=line.count(')'):return 'invalid bootstrap Newick'
                count+=1
        if count!=boot:return f'bootstrap trees={count}, expected={boot}'
    tree=Path(str(phy)+'_phyml_tree.txt').read_text().strip()
    if not tree.endswith(';') or tree.count('(')!=tree.count(')'):return 'invalid final Newick'
    return None

def request(phy,boot,binary,seed=None):
    data=dict(schema=1,input_sha256=digest(phy),binary=str(binary),binary_sha256=digest(binary),
                model='HKY85',categories=4,alpha='e',invariant='e',bootstrap=boot)
    if seed is not None:data['seed']=seed
    return data

def seal(phy,req,execution=None):
    p=Path(str(phy)+'.phyml.complete.json');q=p.with_suffix('.next')
    data=dict(request=req,outputs={str(x):digest(x) for x in outputs(phy,req['bootstrap'])})
    if execution is not None:data['execution']=execution
    q.write_text(json.dumps(data,indent=2)+'\n');q.replace(p)

def reusable(phy,req,adopt=True):
    why=completion_error(phy,req['bootstrap'])
    if why:return False,why
    receipt=Path(str(phy)+'.phyml.complete.json')
    if receipt.exists():
        try:
            data=json.loads(receipt.read_text())
            if data['request']!=req:return False,'input/software/options changed'
            if set(data['outputs'])!={str(x) for x in outputs(phy,req['bootstrap'])}:return False,'incomplete output receipt'
            if any(digest(p)!=d for p,d in data['outputs'].items()):return False,'completed output changed'
        except (OSError,ValueError,KeyError):return False,'invalid completion receipt'
        return True,'verified completion receipt'
    # Adopt historical results only when their own full log proves completion
    # and the input predates the main tree. Subsequent reuse uses content hashes.
    log=Path(str(phy)+'.phyml.log').read_text(errors='replace')
    import shlex
    command=re.search(r'\. Command line:\s*(.*)',log)
    try:args=shlex.split(command.group(1)) if command else []
    except ValueError:args=[]
    if not args or Path(args[0]).name!=Path(req['binary']).name:return False,'historical executable does not match'
    options=[('-i',str(phy)),('-m','HKY85'),('-c','4'),('-a','e'),('-v','e'),('-b',str(req['bootstrap']))]
    if 'seed' in req:options.append(('--r_seed',str(req['seed'])))
    for key,val in options:
        if key not in args or args.index(key)+1>=len(args) or args[args.index(key)+1]!=val:return False,'historical command does not match requested input/options'
    if phy.stat().st_mtime_ns>Path(str(phy)+'_phyml_tree.txt').stat().st_mtime_ns:return False,'input newer than historical tree'
    if adopt:seal(phy,req)
    return True,'verified complete historical bootstrap'

def clean(phy):
    for suffix in SUFFIXES:
        Path(str(phy)+suffix).unlink(missing_ok=True)
    # Only products for this exact input, never other loci/trees.
    for p in list(phy.parent.glob(phy.name+'_phyml_tree.panelB*'))+list(phy.parent.glob(phy.name+'_phyml_tree.*.panelB*')):
        if p.is_file():p.unlink()

def reuse_identical_sibling(phy,req):
    """Only byte-identical alignments and identical settings may share a tree."""
    for source in sorted(phy.parent.glob('haplotypes.evidence.*.phy')):
        if source==phy or digest(source)!=req['input_sha256']:continue
        ok,_=reusable(source,req)
        if not ok:continue
        clean(phy)
        for src,dst in zip(outputs(source,req['bootstrap']),outputs(phy,req['bootstrap'])):
            shutil.copy2(src,dst)
        with Path(str(phy)+'.phyml.log').open('a') as f:
            f.write(f'\n[GU PHYML] Reused identical alignment and model from {source}\n')
        seal(phy,req)
        return source
    return None

def state(phy,scope,status,rc,elapsed=0):
    p=Path(str(phy)+'.phyml.run.status.tsv');q=p.with_suffix('.next')
    q.write_text('scope\tstatus\trc\telapsed_seconds\tphy\ttree\tlog\n'+f'{scope}\t{status}\t{rc}\t{elapsed}\t{phy}\t{phy}_phyml_tree.txt\t{phy}.phyml.log\n');q.replace(p)

def stop_process_group(proc):
    # Distribution launchers may spawn MPI/worker children.
    try: os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError: return
    try: proc.wait(timeout=10)
    except subprocess.TimeoutExpired: pass
    try: os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError: pass
    proc.wait()

def interrupted(signum, frame):
    raise RuntimeError(f'PhyML interrupted by signal {signum}')

def discard_partial(phy):
    for suffix in SUFFIXES:
        if suffix!='.phyml.log':Path(str(phy)+suffix).unlink(missing_ok=True)

def run_attempt(phy,cmd,deadline):
    remaining=None if deadline is None else deadline-time.monotonic()
    proc=None
    try:
        with Path(str(phy)+'.phyml.log').open('w') as f:
            if remaining is not None and remaining<=0:
                f.write('[GU PHYML] Timeout budget exhausted before starting this attempt\n')
                return 124
            # MPI otherwise consumes the caller's remaining task-list lines.
            proc=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=f,
                                  stderr=subprocess.STDOUT,start_new_session=True)
            while True:
                remaining=None if deadline is None else deadline-time.monotonic()
                if remaining is not None and remaining<=0:
                    stop_process_group(proc)
                    return 124
                try:return proc.wait(timeout=60 if remaining is None else min(60,remaining))
                except subprocess.TimeoutExpired:
                    if deadline is not None and time.monotonic()>=deadline:
                        stop_process_group(proc)
                        return 124
                    # Continue enforcing the deadline silently.
    finally:
        if proc is not None and proc.poll() is None:stop_process_group(proc)

def archive_failure(phy,req,attempt,why):
    prefix=str(phy)+f'.phyml.failed.{time.time_ns()}'
    log=Path(prefix+'.log')
    shutil.copy2(str(phy)+'.phyml.log',log)
    Path(prefix+'.json').write_text(json.dumps(dict(request=req,attempt=attempt,
        completion_error=why,log=str(log)),indent=2)+'\n')
    return str(log)

def main():
    p=argparse.ArgumentParser();p.add_argument('--phy',type=Path,required=True);p.add_argument('--scope',required=True)
    p.add_argument('--bootstrap',type=int,default=100);p.add_argument('--timeout',type=float,default=0)
    p.add_argument('--cpus',type=int,default=1)
    p.add_argument('--seed',type=int)
    p.add_argument('--mpi-fallback',choices=['serial','error'],default='serial')
    p.add_argument('--verify-only',action='store_true')
    a=p.parse_args();phy=a.phy.resolve()
    if a.bootstrap<0 or a.timeout<0 or a.cpus<1:p.error('bootstrap/timeout must be nonnegative and cpus positive')
    serial=shutil.which('phyml');mpi=shutil.which('phyml-mpi');mpirun=shutil.which('mpirun')
    if not serial:raise RuntimeError('phyml executable unavailable')
    lock=Path(str(phy)+'.phyml.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    # Changing the execution CPU count alone must not replace a complete tree.
    for binary in dict.fromkeys(x for x in (serial,mpi) if x):
        req=request(phy,a.bootstrap,Path(binary).resolve(),a.seed)
        ok,why=reusable(phy,req,adopt=not a.verify_only)
        if ok:
            if not a.verify_only:state(phy,a.scope,'COMPLETE_VERIFIED',0)
            print(f'[GU PHYML] tree={"VERIFIED" if a.verify_only else "SKIP"} reason={why} input={phy}',flush=True)
            return 0
    if a.verify_only:
        print(f'ERROR: planned tree is not complete: {why}; input={phy}',file=sys.stderr)
        return 1
    binary=serial
    launcher=[]
    if a.cpus > 1 and a.bootstrap > 0:
        if mpi and mpirun:
            # PhyML rounds replicates UP to a multiple of MPI ranks. Choose a
            # divisor so requesting 100 always produces exactly 100, not 112.
            a.cpus=max(n for n in range(1,min(a.cpus,a.bootstrap)+1) if a.bootstrap%n==0)
            binary=mpi;launcher=[mpirun,'--bind-to','none','-np',str(a.cpus)]
        else: print('[GU PHYML] MPI unavailable; using one CPU',flush=True)
    req=request(phy,a.bootstrap,Path(binary).resolve(),a.seed)
    reused=reuse_identical_sibling(phy,req)
    if reused:
        state(phy,a.scope,'COMPLETE_REUSED',0)
        print(f'[GU PHYML] tree=REUSE identical_input={reused} input={phy}',flush=True);return 0
    print(f'[GU PHYML] tree=REPLACE reason={why} input={phy}',flush=True)
    clean(phy);state(phy,a.scope,'RUNNING','');start=time.monotonic()
    deadline=start+a.timeout if a.timeout else None
    cmd=launcher+[binary,'-i',str(phy),'-m','HKY85','-c','4','-a','e','-v','e','-b',str(a.bootstrap)]
    if a.seed is not None:cmd+=['--r_seed',str(a.seed)]
    print(f'[GU PHYML] bootstrap_cpus={a.cpus if launcher else 1} input={phy}',flush=True)
    attempts=[];rc=1
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    try:
        rc=run_attempt(phy,cmd,deadline)
        why=completion_error(phy,a.bootstrap)
        log=Path(str(phy)+'.phyml.log').read_text(errors='replace')
        match=re.search(r'\. Random seed:\s*(\d+)',log)
        seed=int(match.group(1)) if match else a.seed
        attempts.append(dict(backend='mpi' if launcher else 'serial',command=cmd,rc=rc,seed=seed))
        # One computational-error recovery, never a retry for weak support,
        # incomplete zero-exit output, timeout, or interruption. MPI and serial
        # have different random streams even with the same initial seed.
        mpi_error = rc > 0 and rc not in (124, 130, 137, 143) and bool(
            re.search(r'MPI_ERR_[A-Z_]+|MPI_ERRORS_ARE_FATAL', log))
        if launcher and (rc==1 or mpi_error) and a.mpi_fallback=='serial' and seed is not None:
            failed_log=archive_failure(phy,req,attempts[-1],why)
            attempts[-1]['failed_log']=failed_log
            if deadline is not None and time.monotonic()>=deadline:
                rc=124
            else:
                print(f'[GU PHYML] mpi=FAILED rc={rc} fallback=serial seed={seed} failed_log={failed_log}',flush=True)
                clean(phy)
                req=request(phy,a.bootstrap,Path(serial).resolve(),a.seed)
                cmd=[serial,'-i',str(phy),'-m','HKY85','-c','4','-a','e','-v','e',
                     '-b',str(a.bootstrap),'--r_seed',str(seed)]
                rc=run_attempt(phy,cmd,deadline)
                why=completion_error(phy,a.bootstrap)
                attempts.append(dict(backend='serial',command=cmd,rc=rc,seed=seed))
        if rc==0 and why is None:
            seal(phy,req,dict(attempts=attempts));state(phy,a.scope,'COMPLETE',0,int(time.monotonic()-start));return 0
        print(f'ERROR: incomplete PhyML rc={rc}: {why}; input={phy}',file=sys.stderr)
        # Prevent this invocation's summary stages from reading a partial tree.
        discard_partial(phy)
        state(phy,a.scope,'FAILED',rc or 1,int(time.monotonic()-start));return rc or 1
    except Exception:
        discard_partial(phy)
        state(phy,a.scope,'FAILED',1,int(time.monotonic()-start))
        raise

if __name__=='__main__':
    try:sys.exit(main())
    except (OSError,RuntimeError,ValueError) as e:print(f'ERROR: {e}',file=sys.stderr);sys.exit(1)
