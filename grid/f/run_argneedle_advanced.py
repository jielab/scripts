#!/usr/bin/env python3
"""Version-tolerant wrapper around ARG-Needle's official three-step large-data script."""
from __future__ import annotations
import argparse,os,re,shlex,subprocess,sys
from pathlib import Path

def locate(home:str)->Path:
    c=[]
    if home:
      h=Path(home)
      c += [h/"arg_needle/infer_args_advanced.py",h/"infer_args_advanced.py"]
      c += list(h.glob("**/infer_args_advanced.py"))
    try:
      import importlib.util
      spec=importlib.util.find_spec("arg_needle.scripts.infer_args_advanced")
      if spec is not None and spec.origin: c.append(Path(spec.origin))
    except (ImportError, AttributeError, ValueError):
      pass
    for p in os.environ.get("PATH","").split(os.pathsep): c.append(Path(p)/"infer_args_advanced.py")
    return next((x.resolve() for x in c if x.is_file()),None)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--home",default=""); ap.add_argument("--haps",required=True); ap.add_argument("--map",required=True); ap.add_argument("--out",required=True); ap.add_argument("--chr",type=int,required=True); ap.add_argument("--seed-haplotypes",type=int,required=True); ap.add_argument("--threads",type=int,default=8); ap.add_argument("--normalize",type=int,choices=[0,1],default=1); ap.add_argument("--random-seed",type=int,default=20260904); ap.add_argument("--replace",action="store_true"); ap.add_argument("--log",required=True); a=ap.parse_args()
    script=locate(a.home)
    if script is None: raise SystemExit("Cannot find arg_needle/infer_args_advanced.py; set --arg-needle-home")
    env=os.environ.copy()
    if a.home: env["PYTHONPATH"]=a.home+(os.pathsep+env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    helptext=subprocess.check_output([sys.executable,str(script),"--help"],text=True,stderr=subprocess.STDOUT,env=env)
    def opt(*names): return next((x for x in names if x in helptext),None)
    th=opt("--num_threads","--threads","--n_threads")
    rs=opt("--random_seed","--random-seed","--seed")
    prefix=Path(a.out); prefix.parent.mkdir(parents=True,exist_ok=True); Path(a.log).parent.mkdir(parents=True,exist_ok=True)
    import shutil,time
    with open(a.log,"a") as log:
      for step in (1,2,3):
        marker=Path(str(prefix)+f".step{step}.done")
        if marker.is_file() and not a.replace:
            log.write(f"SKIP step={step} marker={marker}\n"); log.flush(); continue
        before={str(x):x.stat().st_mtime_ns for x in prefix.parent.glob(prefix.name+"*.argn")}
        cmd=[sys.executable,str(script),"--hap_gz",a.haps,"--map",a.map,"--out",str(prefix),"--mode","array","--step",str(step),"--chromosome",str(a.chr),"--verbose","1"]
        if step==1: cmd += ["--num_snp_samples",str(a.seed_haplotypes)]
        if step==3 and "--normalize" in helptext: cmd += ["--normalize",str(a.normalize)]
        if th: cmd += [th,str(a.threads)]
        if rs: cmd += [rs,str(a.random_seed)]
        log.write("COMMAND "+" ".join(shlex.quote(x) for x in cmd)+"\n"); log.flush()
        subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT,env=env)
        files=[x for x in prefix.parent.glob(prefix.name+"*.argn") if x.is_file() and x.stat().st_size>0]
        if not files: raise SystemExit(f"ARG-Needle step {step} created no nonempty .argn file for prefix {prefix}")
        marker.write_text("completed\n")
    final=Path(str(prefix)+".argn")
    if not final.is_file() or final.stat().st_size==0:
        files=sorted([x for x in prefix.parent.glob(prefix.name+"*.argn") if x.is_file() and x.stat().st_size>0],key=lambda x:(x.stat().st_mtime_ns,x.stat().st_size))
        if not files: raise SystemExit("No final ARG-Needle .argn file")
        shutil.copy2(files[-1],final)
    print(final)
if __name__=="__main__": main()
