"""Submit one ordered workflow per sample or per selected dataset batch."""
import os
from pathlib import Path
import subprocess
import sys

from common import CODE_ROOT, load_config, resolve_mode, default_cores
from hierarchical_split import parser as worker_parser
from query import discover_samples
from workflow_selection import positive
from project_cache import batch_plan, probe_samples
from index_selection import INDEX_NAMES, select_index


def batch_range_label(numbers):
    ranges = []
    for number in sorted(set(numbers)):
        if ranges and number == ranges[-1][1] + 1:
            ranges[-1][1] = number
        else:
            ranges.append([number, number])
    return '_'.join(str(first) if first == last else f'{first}-{last}' for first,last in ranges)


def main(argv=None):
    p=worker_parser()
    p.set_defaults(threads=None)
    p.add_argument('--batch', help='Batch indices per dataset: 1-2, 1,3-5, 3-, all')
    p.add_argument('--batch-size', type=positive, default=None)
    p.add_argument('--index-use', choices=INDEX_NAMES, default=None, help='Index profile [optimized]; legacy preserves the old group index')
    p.add_argument('--submit', choices=('hpc','local','auto'), default=None)
    p.add_argument('--plan-only', action='store_true')
    p.add_argument('--queue', default='38')
    p.add_argument('--cores', type=positive, default=None, help='CPU budget [HPC: 40; local: 16]')
    p.add_argument('--mem-gib', type=positive, default=448, help='LSF memory reservation in GiB [448; one large group job per 512-GiB node]')
    p.add_argument('--hours', type=positive, default=120)
    a=p.parse_args(argv)
    mode=resolve_mode(a.submit)
    if a.cores is None:a.cores=default_cores(mode)
    if a.sample and (a.batch is not None or a.batch_size is not None):
        p.error('--sample cannot be combined with --batch or --batch-size')
    if a.include_uncertain and 'species' not in a.action:
        p.error('--include-uncertain requires split-species')
    if a.groups and 'species' not in a.action:
        p.error('--groups requires split-species')
    if a.out and 'report' not in a.action:
        p.error('--out requires report')
    if a.threads is not None and not 1 <= a.threads <= a.cores:
        p.error('--threads must be positive and no larger than --cores')
    if a.evaluation_pairs is not None and a.evaluation_pairs < 1:
        p.error('--evaluation-pairs must be positive')
    os.environ['submit']=os.environ['SUBMIT']=mode
    cfg=load_config(a.config,a.root)
    cfg.update(classify_memory_gib=a.mem_gib, classify_cores=a.cores,
               classify_threads=a.threads or max(1,a.cores-2), classify_workers=1)
    if 'library_dir' in cfg:
        if a.routing=='full' and a.index_use not in (None,'full'):p.error('Use --index-use full for full routing')
        if a.routing=='compact' and a.index_use=='full':p.error('--routing compact conflicts with --index-use full')
        cfg=select_index(cfg,'full' if a.routing=='full' else a.index_use)
        print(f"Library: {cfg['library_dir']}; index={cfg['index_use']}; group={cfg['group_library_dir']}",flush=True)
        print(f"Output: {cfg['output_dir']}; state: {cfg['work_dir']}",flush=True)
        if a.action != ['report']:
            from partitioned_build import database_snapshot
            from group_library import compact_snapshot
            if not (Path(cfg['library_dir'])/'catalog.json').is_file():
                p.error('Selected library is not ready; run 1.lib.sh build first')
            database_snapshot(cfg)
            if 'group' in a.action and cfg['group_routing']=='compact':
                try:
                    from ordinal_query import group_query_config
                    group_query_config(cfg, compact_snapshot(cfg))
                except (OSError, ValueError, RuntimeError) as exc:
                    p.error('INPUT NOT READY: ' + str(exc) + '; complete the selected group index first; no job submitted')
    print(f"Discovering cached inputs: datasets={a.dat_list or 'selected/all'}",flush=True)
    samples=discover_samples(cfg,datasets=a.dat_list,only=a.sample,probe=False)
    if not samples:
        p.error('No complete selected input sample pairs')
    if a.truth and len(samples)!=1:
        p.error('--truth requires exactly one selected sample')
    if a.sample:
        jobs=[('metag.sample', s['sample_id'], [s]) for s in samples]
    elif a.action==['report'] and a.batch is None and a.batch_size is None:
        jobs=[('metag.sample', 'selected cohort report', samples)]
    else:
        size=a.batch_size or int(cfg.get('sample_batch_size',10))
        try:
            plan=batch_plan(cfg,samples,size,a.batch)
            jobs=[(f'metag.{dataset}.batch{number}', f'{dataset}/batch{number:04d}', members)
                  for dataset,number,members in plan]
        except ValueError as exc:
            p.error(str(exc))
    # Prefix checks belong only to selected batches; unrelated raw FASTQs stay closed.
    if 'reads_raw' in cfg:
        selected={s['sample_id']:s for _,_,members in jobs for s in members}
        probe_samples(cfg,list(selected.values()))
    if 'library_dir' in cfg:
        from submission_inputs import check_hierarchy_stages
        try:
            check_hierarchy_stages(cfg, a.action, [s for _,_,members in jobs for s in members], a.stage)
        except (OSError, ValueError, RuntimeError) as exc:
            p.error('INPUT NOT READY: ' + str(exc) + '; no job submitted')
    if a.out and len(jobs)>1:
        p.error('--out requires a single job; omit it for separate per-job reports')
    common=['--submit',mode,'--queue',a.queue,'--cores',str(a.cores),
            '--hours',str(a.hours)]
    if 'library_dir' in cfg:
        common+=['--index-use',cfg['index_use']]
    if a.mem_gib is not None:common+=['--mem-gib',str(a.mem_gib)]
    if a.plan_only:common+=['--plan-only']
    for flag,value in [('--config',a.config),('--root',a.root),('--groups',a.groups),
                       ('--truth',a.truth),('--out',a.out),('--threads',a.threads),
                       ('--evaluation-pairs',a.evaluation_pairs)]:
        if value is not None:common += [flag,str(value)]
    common+=['--replace','TRUE' if a.replace else 'FALSE','--stage',a.stage]
    common+=['--evaluation',a.evaluation,'--routing',cfg.get('group_routing') or a.routing or 'compact']
    if a.include_uncertain:common+=['--include-uncertain']
    if a.allow_local:common+=['--allow-local']
    # All selection validation finishes before the first job is submitted.
    for name,label,members in jobs:
        print(f'{name}: {label}; samples={len(members)}; stages={",".join(a.action)}',flush=True)
        selected=[x for s in members for x in ('--sample',s['sample_id'])]
        command=[sys.executable,str(CODE_ROOT/'f/qc_submit.py'),'hierarchical_split',
                 ','.join(a.action),'--job-name',name,*common,*selected]
        subprocess.run(command,check=True)


if __name__=='__main__':
    main()
