#!/usr/bin/env python3
"""Pure joint selection/publication after the exact planned ancestry×chromosome set."""
from __future__ import annotations
from contextlib import ExitStack
import os
from pathlib import Path
from ld_pipeline import Pipeline, POPS, Stop, choose_distance, read_json, sha, stamp, write_json


def invalidate_publication(root, build):
    """An in-progress replacement must not leave apparently current final BEDs."""
    root=Path(root)
    for pop in ('ALL',*POPS):
        for suffix in ('','.SUBSET'):
            (root/f'{pop}.{build}{suffix}.bed').unlink(missing_ok=True)
    for name in ('COMPLETE.json','NEEDS_LARGER_GRID.json'):
        (root/name).unlink(missing_ok=True)


def activate_plan(plan_path):
    plan=read_json(plan_path)
    if not plan or plan.get('schema')!='ld.plan.1':raise Stop('Invalid LD plan')
    pipeline=Pipeline();root=Path(plan['output_root'])
    with pipeline.output_locks(out=root):
        invalidate_publication(root,plan['genome_build'])
        write_json(root/'CURRENT_PLAN.json',{'run_id':plan['run_id'],'plan':str(Path(plan_path).resolve())})


def merge(plan_path=None):
    plan_path=plan_path or os.environ.get('LD_PLAN_FILE')
    plan=read_json(plan_path or '')
    if not plan or plan.get('schema')!='ld.plan.1':raise Stop('Merge requires its exact LD plan')
    pipeline=Pipeline();root=Path(plan['output_root']).resolve()
    if pipeline.out!=root or pipeline.genome_build!=plan['genome_build'] or pipeline.chroms!=plan['chromosomes']:
        raise Stop('Merge environment does not match its planned output/build/chromosomes')
    expected=[f'{p}.chr{c}' for c in pipeline.chroms for p in POPS]
    if plan.get('populations')!=list(POPS) or [j['key'] for j in plan['jobs']]!=expected:
        raise Stop('Plan must contain exactly the four populations for every requested chromosome in canonical order')
    if plan['grid_kb']!=pipeline.grid or plan['confirmations']!=pipeline.confirmations:
        raise Stop('Merge D grid differs from its population workers')
    with pipeline.output_locks(),ExitStack() as locks:
        active=read_json(root/'CURRENT_PLAN.json')
        if not active or active.get('run_id')!=plan['run_id']:
            raise Stop('This merge plan has been superseded; current results were not changed')
        invalidate_publication(root,pipeline.genome_build)
        records={};common=None;references={};evidence={}
        panel_identity=sha(pipeline.panel)
        pipeline.population_roots={}
        for job in plan['jobs']:
            pop,chrom=job['population'],job['chromosome'];source=Path(job['output']).resolve()
            if source!=root/'parts'/pop/f'chr{chrom}' or job['key']!=f'{pop}.chr{chrom}':
                raise Stop('Wrong worker output or ancestry/chromosome identity')
            locks.enter_context(pipeline.output_locks(out=source))
            done=read_json(Path(plan['run_dir'])/'receipts'/(job['key']+'.json'))
            record=read_json(source/'COMPLETE.json')
            if not done or done!={'run_id':plan['run_id'],'key':job['key'],'status':'complete',
                                 'output':str(source),'population':pop,'chromosome':chrom}:
                raise Stop('Missing current-plan successful worker receipt: '+job['key'])
            if (not record or record.get('schema')!='ld.population.1' or record.get('status')!='complete'
                    or record.get('run_id')!=plan['run_id'] or record.get('population')!=pop
                    or record.get('chromosome')!=chrom or record.get('grid_kb')!=pipeline.grid):
                raise Stop('Incomplete/stale/wrong worker outputs: '+job['key'])
            signature={k:record[k] for k in ('version','genome_build','input_format','sample_counts','panel_sha256',
                       'grid_kb','confirmations','coincident_policy','reference_count_policy','reference_max_rel_deviation','beta_bp')}
            if common is not None and signature!=common:raise Stop('Incompatible population/chromosome settings: '+job['key'])
            common=signature
            if record['panel_sha256']!=panel_identity:raise Stop('Sample panel identity differs from population results')
            pipeline.sample_counts=record['sample_counts'];pipeline.tool=record['plink2']
            pgen=record['pgen_record'];spec=pgen['spec'];qc=pgen['qc']
            if (spec['pop']!=pop or spec['chrom']!=chrom or spec['sample_count']!=pipeline.sample_counts[pop]
                    or qc.get('missing_genotypes')!=0 or not qc.get('phase_check_passed')):
                raise Stop('Invalid population sample/genotype-QC identity: '+job['key'])
            reference=spec.get('source_pfiles',spec.get('vcf'))
            if chrom in references and references[chrom]!=reference:raise Stop('Reference input changed between ancestry workers')
            if chrom not in references:
                current=({name:stamp(path) for name,path in pipeline.input_pfiles(chrom).items()}
                         if pipeline.input_format=='pgen' else stamp(pipeline.vcf_dir/f'chr{chrom}.vcf.gz'))
                if current!=reference:raise Stop('Reference input changed during or after population computation')
            references[chrom]=reference
            pipeline.population_roots[(pop,chrom)]=source
            pipeline.pgen_cache[(pop,chrom)]=(None,pgen)
            if not all(pipeline.partition_valid(pop,chrom,d,pgen) for d in pipeline.grid):
                raise Stop('Missing/incomplete strict partition: '+job['key'])
            records[(pop,chrom)]=record
            evidence[job['key']]={'path':str((source/'COMPLETE.json').relative_to(root)),
                                   'run_id':record['run_id']}
        selections={};unstable=[]
        for chrom in pipeline.chroms:
            snapshots=pipeline.summarize(chrom,pipeline.grid)
            selected=choose_distance(snapshots,pipeline.grid,pipeline.confirmations)
            if selected is None:
                last=pipeline.grid[-1];signature=snapshots[last]['signature'];tail=[]
                for d in reversed(pipeline.grid):
                    if snapshots[d]['signature']!=signature:break
                    tail.append(d)
                span=max((record['pgen_record']['qc']['last_position']-record['pgen_record']['qc']['first_position']+999)//1000
                         for (pop,c),record in records.items() if c==chrom)
                unstable.append(dict(chromosome=chrom,tested_grid_kb=pipeline.grid,identical_suffix_kb=sorted(tail),
                                     required_larger_confirmations=pipeline.confirmations,
                                     suggested_max_d_kb=min(2*last,span) if span>last else last))
                continue
            pipeline.check_reference_counts(chrom,snapshots[selected],context=f'stable common D={selected} kb')
            selected_record=pipeline.audited_selection_record(chrom,selected,pipeline.grid,snapshots)
            write_json(root/'selected'/f'chr{chrom}.json',selected_record)
            selections[chrom]=(selected,selected_record)
        if unstable:
            write_json(root/'NEEDS_LARGER_GRID.json',dict(status='joint_grid_not_stable',run_id=plan['run_id'],
                       chromosomes=unstable,action='Increase --max-d-kb and rerun the same chromosome set with --replace FALSE; no final BED was published'))
            raise Stop('Joint four-population boundary stability failed; see NEEDS_LARGER_GRID.json. No final BED was published')
        pipeline.publish_context=dict(run_id=plan['run_id'],population_jobs=evidence,
                                      computation='independent_population_chromosome_CPU_jobs',
                                      distance_policy='same explicit common grid; no hidden extension')
        try:
            pipeline.publish(selections)
        except BaseException:
            invalidate_publication(root,pipeline.genome_build)
            raise
        print(f'[LD] MERGED {len(expected)} successful population/chromosome jobs in chromosome order: {root}',flush=True)


if __name__=='__main__':
    try:merge()
    except (Stop,OSError,KeyError,ValueError,TypeError) as exc:raise SystemExit(f'[LD] MERGE STOP: {exc}')
