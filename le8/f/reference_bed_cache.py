#!/usr/bin/env python3
"""Publish one verified BED conversion per source/population for all MR-link jobs."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile


def stamp(path):
    p=Path(path).resolve(strict=True);s=p.stat()
    return [str(p),s.st_size,s.st_mtime_ns]


def outputs(root):
    return [[ext,(root/('reference.'+ext)).stat().st_size,(root/('reference.'+ext)).stat().st_mtime_ns]
            for ext in ('bed','bim','fam')]


def prepare(source, kind, keep, root, plink='plink2'):
    binary=shutil.which(plink)
    if binary is None:raise ValueError('PLINK2 is required')
    args=[binary]
    if kind=='pfile':
        pvar=Path(str(source)+'.pvar')
        args+=['--pfile',str(source)]
        if not pvar.exists():pvar=Path(str(source)+'.pvar.zst');args+=['vzs']
        inputs=[str(source)+'.pgen',str(source)+'.psam',pvar]
    else:
        args+=['--bfile',str(source)];inputs=[str(source)+'.'+ext for ext in ('bed','bim','fam')]
    request=dict(inputs=[stamp(x) for x in inputs],binary=stamp(binary),kind=kind,
                 keep=hashlib.sha256(Path(keep).read_bytes()).hexdigest() if keep else None,
                 policy='max-alleles=2',code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    key=hashlib.sha256(json.dumps(request,sort_keys=True).encode()).hexdigest()
    root=Path(root);root.mkdir(parents=True,exist_ok=True);cache=root/key
    with (root/(key+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        try:
            record=json.loads((cache/'complete.json').read_text())
            valid=record['request']==request and record['outputs']==outputs(cache) and all(row[1]>0 for row in record['outputs'])
        except (OSError,ValueError,KeyError):valid=False
        if not valid:
            with tempfile.TemporaryDirectory(prefix='.build.',dir=root) as temporary:
                staging=Path(temporary)/'ready';staging.mkdir()
                if keep:args+=['--keep',str(keep)]
                args+=['--max-alleles','2','--make-bed','--out',str(staging/'reference')]
                with (staging/'build.log').open('w') as log:
                    proc=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT)
                    try:
                        if proc.wait():raise RuntimeError('PLINK reference conversion failed')
                    finally:
                        if proc.poll() is None:
                            proc.terminate()
                            try:proc.wait(timeout=10)
                            except subprocess.TimeoutExpired:proc.kill();proc.wait()
                result=outputs(staging)
                if not all(row[1]>0 for row in result):raise ValueError('Incomplete reference BED')
                (staging/'complete.json').write_text(json.dumps(dict(request=request,outputs=result)))
                if cache.exists():shutil.rmtree(cache)
                staging.replace(cache)
    return cache/'reference'


def main():
    def stopped(signum,frame):raise SystemExit(128+signum)
    for sig in (signal.SIGTERM,signal.SIGHUP):signal.signal(sig,stopped)
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--kind',choices=['pfile','bfile'],required=True)
    p.add_argument('--keep');p.add_argument('--cache-root',required=True)
    a=p.parse_args();print(prepare(a.source,a.kind,a.keep,a.cache_root))

if __name__=='__main__':main()
