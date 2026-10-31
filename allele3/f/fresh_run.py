"""Fresh output preparation without deleting previous calculations or inputs."""
from contextlib import contextmanager
import datetime
import fcntl
import os
from pathlib import Path
import uuid


def boolean(value):
    value=str(value).upper()
    if value not in {'TRUE','FALSE'}:raise ValueError('--replace must be TRUE or FALSE')
    return value=='TRUE'


@contextmanager
def output_run(output, replace, protected=()):
    output=Path(output).expanduser().resolve()
    analysis=Path(os.environ.get('ALLELE_ANALYSIS_ROOT','/mnt/d/analysis/allele3')).resolve()
    if output==analysis or output in analysis.parents:raise ValueError('Choose a separate task output, not the analysis root or its parent')
    for original in protected:
        if not original:continue
        original=Path(original).expanduser().resolve()
        if output==original or original in output.parents or output in original.parents:
            raise ValueError(f'Output overlaps a preserved input/cache: {original}; choose a new --out')
    output.parent.mkdir(parents=True,exist_ok=True)
    with (output.parent/('.'+output.name+'.recalculation.lock')).open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc:raise RuntimeError('Another calculation owns this output directory') from exc
        held=[]
        try:
            for name in ('job.lock','run.lock','.allele3.lock','.pca3d.lock'):
                p=output/name
                if p.exists():
                    h=p.open('a');held.append(h)
                    try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
                    except BlockingIOError as exc:raise RuntimeError('Another job is active: '+str(p)) from exc
            old=[p for p in output.iterdir() if p.name not in {'logs','job.lock','run.lock','.allele3.lock','.pca3d.lock'}] if output.is_dir() else []
            if replace and old:
                stamp=datetime.datetime.now().strftime('%Y%m%d_%H-%M-%S')+'_'+uuid.uuid4().hex[:6]
                archive=output.parent/(output.name+'.previous_'+stamp);archive.mkdir()
                for p in old:p.rename(archive/p.name)
                print('[preserved previous output] '+str(archive),flush=True)
            output.mkdir(parents=True,exist_ok=True)
        finally:
            for handle in held:handle.close()
        print('[mode] '+('FRESH RECALCULATION (--replace TRUE)' if replace else 'CHECKPOINT RESUME (--replace FALSE)'),flush=True)
        print('[new output] '+str(output),flush=True)
        for original in protected:
            if original:print('[permanent input/cache; retained] '+str(original),flush=True)
        yield output
