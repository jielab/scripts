#!/usr/bin/env python3
"""Request-bound, validated ARG-Needle checkpoints and visible progress."""
from __future__ import annotations
import argparse, fcntl, hashlib, json, os, shlex, signal, subprocess, sys, time
from pathlib import Path

def locate(home):
    paths=list(Path(home).glob('**/infer_args_advanced.py')) if home else []
    try:
        import importlib.util
        s=importlib.util.find_spec('arg_needle.scripts.infer_args_advanced')
        if s and s.origin: paths.append(Path(s.origin))
    except (ImportError, AttributeError, ValueError): pass
    return next((p.resolve() for p in paths if p.is_file()), None)

def stamp(p):
    p=Path(p).resolve(strict=True); s=p.stat()
    return [str(p),s.st_size,s.st_mtime_ns]

def main():
    def stopped(signum, frame): raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM, stopped)
    p=argparse.ArgumentParser()
    p.add_argument('--home',default='')
    for k in ('haps','map','out','log'): p.add_argument('--'+k,required=True)
    p.add_argument('--chr',type=int,required=True)
    p.add_argument('--seed-haplotypes',type=int,required=True)
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--normalize',type=int,choices=(0,1),default=1)
    p.add_argument('--random-seed',type=int,default=20260904)
    p.add_argument('--mode',choices=('array','sequence'),default='array')
    p.add_argument('--normalize-demography',default='')
    p.add_argument('--replace',action='store_true')
    a=p.parse_args(); script=locate(a.home)
    if script is None: raise SystemExit('Cannot find infer_args_advanced.py')
    env=os.environ.copy()
    if a.home: env['PYTHONPATH']=a.home+os.pathsep+env.get('PYTHONPATH','')
    env['PYTHONUNBUFFERED']='1'; env['PYTHONHASHSEED']=str(a.random_seed)
    helptext=subprocess.check_output([sys.executable,str(script),'--help'],env=env,text=True)
    prefix=Path(a.out); prefix.parent.mkdir(parents=True,exist_ok=True)
    Path(a.log).parent.mkdir(parents=True,exist_ok=True)
    lock=open(str(prefix)+'.lock','a')
    try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: raise SystemExit(f'Another inference uses {prefix}')
    sample=a.haps.removesuffix('.gz').removesuffix('.haps').removesuffix('.hap')+'.sample'
    n=2*len(Path(sample).read_text().splitlines()[2:])
    if not 2<=a.seed_haplotypes<=n: raise SystemExit('Invalid scaffold size')
    import importlib.metadata
    import arg_needle
    package=Path(arg_needle.__file__).parent
    memory_adapter=Path(__file__).resolve().with_name('argneedle_memory.py')
    from argneedle_memory import backend_path
    backend=backend_path()
    if not backend.is_file(): raise SystemExit(f'Missing bounded ASMC backend: {backend}; run f/build_asmc_memory.py --source /path/to/ASMC-1.4.0')
    dependencies=[memory_adapter,backend,package/'inference.py',package/'decoders.py',package/'resources/30-100-2000_CEU.decodingQuantities.gz',package/'resources/CEU.demo']
    request={'schema':3,'inputs':[stamp(x) for x in (a.haps,a.map,sample,script,__file__,*dependencies) if Path(x).is_file()],
      'versions':{name:importlib.metadata.version(name) for name in ('arg-needle','arg-needle-lib')},
      'mode':a.mode,'scaffold':a.seed_haplotypes,'seed':a.random_seed,'normalize':a.normalize,
      'chr':a.chr,'trim':0,'demography':stamp(a.normalize_demography) if a.normalize_demography else 'default CEU'}
    rid=hashlib.sha256(json.dumps(request,sort_keys=True).encode()).hexdigest()
    paths={1:Path(str(prefix)+'.step1.argn'),2:Path(str(prefix)+'.step2.argn'),3:Path(str(prefix)+'.argn')}
    bootstrap="import runpy,sys,random,numpy as np; adapter=sys.argv.pop(1); runpy.run_path(adapter)['install'](); seed=int(sys.argv.pop(1)); random.seed(seed); np.random.seed(seed); script=sys.argv.pop(1); sys.argv[0]=script; runpy.run_path(script,run_name='__main__')"
    dirty=a.replace
    with open(a.log,'a',buffering=1) as log:
      th=next((x for x in ('--num_threads','--threads','--n_threads') if x in helptext),None)
      if th is None:
        msg='ARG-Needle inference is single-threaded in this version; --threads controls preparation; --jobs controls chromosome concurrency.'
        print(msg,flush=True); log.write(msg+'\n')
      for step in (1,2,3):
        output=paths[step]; marker=Path(str(prefix)+f'.step{step}.done'); valid=False
        if not dirty and marker.is_file() and output.is_file():
          try: valid=json.loads(marker.read_text())=={'request':rid,'output':stamp(output)} and output.stat().st_size>0
          except (ValueError,OSError): pass
        if valid:
          print(f'SKIP verified chr{a.chr} step={step}',flush=True); continue
        dirty=True
        for ds in range(step,4):
          Path(str(prefix)+f'.step{ds}.done').unlink(missing_ok=True); paths[ds].unlink(missing_ok=True)
        cmd=[sys.executable,'-u','-c',bootstrap,str(memory_adapter),str(a.random_seed),str(script),
          '--hap_gz',a.haps,'--map',a.map,'--out',str(prefix),'--mode',a.mode,
          '--step',str(step),'--chromosome',str(a.chr),'--verbose','1','--trim_num_snps','0']
        if step==1: cmd+=['--num_snp_samples' if a.mode=='array' else '--num_sequence_samples',str(a.seed_haplotypes)]
        if step==3: cmd+=['--normalize',str(a.normalize)]
        if a.normalize_demography: cmd+=['--normalize_demography',a.normalize_demography]
        if th: cmd+=[th,str(a.threads)]
        log.write('COMMAND '+shlex.join(cmd)+'\n')
        print(f'INFER chr{a.chr} step={step}/3 haplotypes={n} scaffold={a.seed_haplotypes} log={a.log}',flush=True)
        start=time.monotonic(); proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
        try:
          while True:
            try: code=proc.wait(timeout=60); break
            except subprocess.TimeoutExpired:
              with open(a.log,'rb') as tail:
                tail.seek(max(0,tail.seek(0,2)-1500)); last=tail.read().decode(errors='replace').splitlines()
              print(f'PROGRESS chr{a.chr} step={step}/3 elapsed={int(time.monotonic()-start)}s '+(last[-1][:180] if last else ''),flush=True)
          if code: raise subprocess.CalledProcessError(code,cmd)
        finally:
          if proc.poll() is None:
            proc.terminate()
            try: proc.wait(timeout=10)
            except subprocess.TimeoutExpired: proc.kill(); proc.wait()
        if not output.is_file() or output.stat().st_size==0: raise SystemExit(f'Step {step} did not create {output}')
        import arg_needle_lib
        arg=arg_needle_lib.deserialize_arg(str(output))
        if arg.num_samples()!=(a.seed_haplotypes if step==1 else n): raise SystemExit(f'Step {step}: wrong sample count')
        del arg
        tmp=Path(str(marker)+'.next'); tmp.write_text(json.dumps({'request':rid,'output':stamp(output)})+'\n'); tmp.replace(marker)
    print(paths[3],flush=True)

if __name__=='__main__': main()
