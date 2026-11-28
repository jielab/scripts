"""Plan split and downstream evidence work for exactly the selected raw samples."""
import argparse
from contextlib import contextmanager
import threading
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import hashlib
import stat
import time
from functools import cached_property
from itertools import islice
import os
from pathlib import Path
import json
import sys

@contextmanager
def progress_phase(name):
    """Report slow startup/metadata work even while its main thread waits on I/O."""
    if __name__ != '__main__':
        yield
        return
    started = time.monotonic()
    stopped = threading.Event()
    print(f'PHASE {name} start', flush=True)

    def heartbeat():
        while not stopped.wait(5):
            print(f'PHASE {name} running elapsed={time.monotonic()-started:.1f}s', flush=True)

    reporter = threading.Thread(target=heartbeat, daemon=True)
    reporter.start()
    status = 'done'
    try:
        yield
    except BaseException:
        status = 'failed'
        raise
    finally:
        stopped.set()
        reporter.join()
        print(f'PHASE {name} {status} elapsed={time.monotonic()-started:.1f}s', flush=True)


with progress_phase('imports'):
    ROOT = Path(__file__).resolve().parent
    sys.path[:0] = [str(ROOT/'f'),str(ROOT/'refine')]
    from common import atomic_json, locked, resolve_mode
    from rank_workflow import FORMAT, effective_config, release, stamp, executable_path, _mapper, positive
    from rank_plan import automatic_batches
    from rank_paper import read_sample
    from ranked_storage import recover_directory
    from query import discover_samples
    from project_cache import batch_plan
    from scheduler import write_cmd, submit_cmd
    from storage_layout import new_run
    from taxonomy import safe_name
    import run as refine
    import cache

def check_file_sizes(folder, files, sample, workers=8, progress_seconds=5):
    """Check every recorded file with one stat, bounding shared-storage requests."""
    total = len(files)
    started = last = time.monotonic()
    checked = 0

    def report():
        print(f'CHECK sample={sample} files={checked}/{total} '
              f'elapsed={time.monotonic()-started:.1f}s', flush=True)

    def check(item):
        name, expected = item
        try:
            info = (folder/name).stat()
        except (FileNotFoundError, NotADirectoryError):
            return name
        return None if stat.S_ISREG(info.st_mode) and info.st_size == expected else name

    report()
    items = iter(files.items())
    # Keep only workers requests in flight, rather than queueing thousands of futures.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(check, item) for item in islice(items, workers)}
        while pending:
            done, pending = wait(pending, timeout=progress_seconds,
                                 return_when=FIRST_COMPLETED)
            for future in done:
                mismatch = future.result()
                checked += 1
                if mismatch is not None:
                    report()
                    return mismatch
                item = next(items, None)
                if item is not None:
                    pending.add(pool.submit(check, item))
            now = time.monotonic()
            if now-last >= progress_seconds:
                report()
                last = now
    report()
    return None


class Inventory:
    def __init__(self,cfg,opts):
        self.cfg = cfg
        self.opts = opts
        self.observed = {}

    @cached_property
    def common_context(self):
        c = self.cfg
        opts = c.get('stage_alignment') or {}
        keys = ('rank_panel_gib','rank_max_refs','rank_min_windows','rank_panel_taxids',
                'confidence','minimum_hit_groups','minimum_base_quality')
        return dict(code=c['rank_method'],database=c['source_full_index_identity'],
                    programs=[stamp(executable_path(c)),stamp(_mapper(c,opts))],
                    catalog=[stamp(opts['representatives'])],
                    settings={k:c[k] for k in keys},reference_raw=c['reference_raw'],
                    metadata_dir=c['metadata_dir'], metadata=[stamp(p) for p in
                    sorted((Path(c['metadata_dir'])/'summaries').glob('*')) if p.is_file()])

    @cached_property
    def assay_keys(self):
        return {a:cache.assay_key(a,self.opts) for a in cache.requested(self.opts)}

    def inspect(self,sample,folder,refined,condition):
        if folder.with_name('.'+folder.name+'.metag-rotation').exists():
            with locked(folder.parent/('.'+folder.name+'.runs')/'publish.lock'):
                recover_directory(folder)
        with read_sample(folder):
            if not folder.exists():
                return True,True
            if folder.is_symlink():
                raise ValueError('Expected a regular ranked sample directory: '+str(folder))
            receipt = folder/'summary.json'
            receipt_bytes = receipt.read_bytes()
            old = json.loads(receipt_bytes)
            key = sample['sample_id']
            if old.get('sample')!=key or old.get('format')!=FORMAT:
                raise ValueError('Directory belongs to another sample or format: '+str(folder))
            if self.opts['replace']:
                return True,True
            previous=dict(old.get('request',{}))
            context=dict(self.common_context,raw=[stamp(sample[k]) for k in ('r1','r2')])
            files=old['file_sizes']
            mismatch=[k for k in sorted(set(previous)|set(context)) if previous.get(k)!=context.get(k)]
            if not mismatch:
                for ref in old['reference_states']:
                    if not Path(ref[0]).is_file() or stamp(ref[0])!=ref:
                        mismatch.append('reference '+ref[0])
                        break
                if not files or not {'assignments.tsv.gz','taxa.tsv','taxon_evidence.tsv'} <= files.keys():
                    mismatch.append('file inventory')
                else:
                    changed = check_file_sizes(folder, files, key)
                    if changed is not None:
                        mismatch.append(changed)
            if mismatch:
                raise ValueError(f'{key}: existing split differs ({", ".join(mismatch)}). '
                                 f'Inspect it before using run --sample {key} --replace TRUE')
            self.observed[str(folder)]=hashlib.sha256(receipt_bytes).hexdigest()
            if not (refined/'receipt.json').is_file():
                return False,True
            with locked(refined.parent/('.'+refined.name+'.lock'),shared=True):
                base_done=cache.complete(refined,cache.receipt(refined),key,self.observed[str(folder)],
                                         self.opts,condition,self.assay_keys)
                audit_done=base_done and refine.seal_audit.current(
                    Path(self.cfg['index_base_work_dir'])/'refine',key,self.observed[str(folder)],self.opts)
                return False,not audit_done

def job(cfg, folder, name, args, worker, job_name):
    spec=folder/(name+'.json')
    atomic_json(spec,args)
    cmd=write_cmd(cfg,folder/(name+'.cmd'),'classify',
                  [sys.executable,'-B',str(worker),'--worker',str(spec)],job_name)
    return cmd

def split_jobs(cfg, folder, label, members, staged, replace):
    commands=[]
    for stage in (('query','finish') if staged else ('combined',)):
        resources=dict(cfg)
        if stage=='finish':
            resources.update(classify_queue=cfg.get('rank_finish_queue','33'),
                             classify_memory_gib=cfg.get('rank_finish_memory_gib',80))
        name=label+('.'+stage if staged else '')
        payload=dict(config=resources,samples=members,replace=replace,stage=stage,
                     prepared_dir=str(folder/(label+'.prepared')))
        commands.append(job(resources,folder,name,payload,ROOT/'f/rank_workflow.py','metag.'+name))
    return commands

def main():
    p=argparse.ArgumentParser(prog='./metag.sh run')
    p.add_argument('--sample',action='append')
    for name in ('dat-list','batch','config','root','representatives','panel-taxids','queue','metadata','refine-queue'):
        p.add_argument('--'+name)
    for name in ('cores','threads','mem-gib','hours','batch-size','min-windows','max-refs','refine-cores','refine-mem-gib'):
        p.add_argument('--'+name,type=positive)
    p.add_argument('--submit',choices=('auto','local','hpc'))
    p.add_argument('--replace',choices=('TRUE','FALSE'),default='FALSE')
    p.add_argument('--confidence',type=float)
    p.add_argument('--panel-gib',type=float)
    p.add_argument('--seal-audit',action='store_true')
    p.add_argument('--seal-db')
    p.add_argument('--blast-pairs',type=positive,default=100)
    for name in ('blast-db','amr-markers','amr-annotation'):
        p.add_argument('--'+name,default=os.getenv('METAG_'+name.upper().replace('-','_')))
    p.add_argument('--plan-only',action='store_true')
    a=p.parse_args()
    if a.sample and (a.batch or a.batch_size):
        p.error('--sample excludes --batch and --batch-size')
    if a.amr_annotation and not a.amr_markers:
        p.error('--amr-annotation requires --amr-markers')
    with progress_phase('configuration'):
        cfg=effective_config(a)
    with progress_phase('code-fingerprint'):
        key,_=release(cfg)
    cfg.update(rank_release=key,code_root=str(ROOT))
    with progress_phase('sample-discovery'):
        samples=discover_samples(cfg,datasets=a.dat_list,only=a.sample,probe=False)
    if not samples:
        p.error('No complete paired-end samples selected')
    automatic=not (a.sample or a.batch or a.batch_size)
    with progress_phase('batch-planning'):
        plans=(automatic_batches(samples,cfg) if automatic else
               batch_plan(cfg,samples,1 if a.sample else a.batch_size or cfg['sample_batch_size'],a.batch or 'all'))
    staged=resolve_mode(a.submit)=='hpc' and automatic
    rebuild=a.replace=='TRUE'
    root=Path(cfg['index_base_work_dir'])/'refine'
    cores=a.refine_cores or (4 if a.blast_db or a.amr_markers or a.seal_audit or a.seal_db else 1)
    opts=dict(replace=rebuild,cores=cores,blast_pairs=a.blast_pairs,
              blast_db=str(Path(a.blast_db).expanduser().absolute()) if a.blast_db else None,
              amr_markers=str(Path(a.amr_markers).expanduser().resolve(strict=True)) if a.amr_markers else None,
              amr_annotation=str(Path(a.amr_annotation).expanduser().resolve(strict=True)) if a.amr_annotation else None)
    opts['seal']=refine.seal_panel.configure(a,cfg,root)
    small=dict(cfg,work_dir=str(root),classify_memory_gib=a.refine_mem_gib or 32,
               classify_cores=cores,classify_hours=24,classify_queue=a.refine_queue or cfg['convert_queue'])
    with progress_phase('refine-metadata'):
        labels=refine.meta(a.metadata)
        inv=Inventory(cfg,opts)
        split_folder=ref_folder=None
        ref_version=refine.ref_version()
    chains=[]
    decisions=[]
    total_samples = sum(len(members) for _, _, members in plans)
    print(f'CHECK selected_samples={total_samples} batches={len(plans)}', flush=True)
    for dataset,index,members in plans:
        label=safe_name(dataset)+f'.batch{index:03d}'
        pending=[]
        independent=[]
        after=[]
        for s in members:
            sid=s['sample_id']
            if a.metadata and sid not in labels:
                p.error('No metadata for '+sid)
            path=Path(cfg['output_dir'])/sid
            cond=labels.get(sid)
            started = time.monotonic()
            print(f'CHECK sample={sid} sample_index={len(decisions)+1}/{total_samples} '
                  f'batch={label} phase=inspect', flush=True)
            needs_split,needs_ref=inv.inspect(s,path,root/'main'/sid,cond)
            decisions.append((sid,needs_split,needs_ref))
            print(f'CHECK sample={sid} split={"run" if needs_split else "reuse"} '
                  f'refine={"run" if needs_ref else "reuse"} '
                  f'elapsed={time.monotonic()-started:.1f}s', flush=True)
            if needs_split:
                pending.append(s)
            if needs_ref:
                (after if needs_split else independent).append([str(path),cond])
        split=[]
        if pending:
            if split_folder is None:
                split_folder=new_run(Path(cfg['work_dir'])/'jobs/split')
            split=split_jobs(cfg,split_folder,label,pending,staged,rebuild)
        for suffix,inputs,parents in [('ready',independent,[]),('after',after,split)]:
            if not inputs:
                continue
            if ref_folder is None:
                ref_folder=new_run(root/'jobs')
            chunk=4 if cache.requested(opts) or opts.get('seal') else len(inputs)
            tasks=[]
            for first in range(0,len(inputs),chunk):
                name=label+'.'+suffix+(f'.{first//chunk+1}' if len(inputs)>chunk else '')
                tasks.append(job(small,ref_folder,name,dict(inputs=inputs[first:first+chunk],root=str(root),
                                  version=ref_version,opts=opts,expected_splits={
                                      path:inv.observed[path] for path,_ in inputs[first:first+chunk] if path in inv.observed}),ROOT/'refine/run.py','metag.refine.'+name))
            chains.append((parents,tasks))
    for sid,sp,rf in decisions:
        print(sid,'split='+('run' if sp else 'reuse'),'refine='+('run' if rf else 'reuse'),sep='\t')
    print(f'Selected {len(decisions)} samples; split {sum(s for _,s,_ in decisions)}; refine {sum(r for _,_,r in decisions)}.',flush=True)
    for split,tasks in chains:
        print(' -> '.join(map(str,split))+(' -> ' if split else '')+', '.join(map(str,tasks)),flush=True)
        if a.plan_only:
            continue
        parent=None
        if split:
            parent=submit_cmd(cfg,split[0])['job_id']
            if staged:
                parent=submit_cmd(cfg,split[1],dependency=parent)['job_id']
        for task in tasks:
            submit_cmd(small,task,dependency=parent if resolve_mode(a.submit)=='hpc' else None)

if __name__=='__main__':
    main()
