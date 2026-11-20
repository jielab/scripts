"""Chromosome workers for the two-pass enrichment estimator.

Fork shares read-only reference indexes. Each worker opens its own input/FASTA.
All chromosomes finish pass 1 before any leave-one-chromosome-out rate is used.
Results are reduced in chromosome order, independently of scheduling order.
"""
import copy
import multiprocessing as mp
import os
import time
from collections import Counter, defaultdict

import numpy as np

_SHARED = None


def workers_requested(a):
    allocated = int(os.environ.get('ALLELE_ALLOCATED_CORES') or
                    os.environ.get('LSB_DJOB_NUMPROC') or os.cpu_count() or 1)
    requested = getattr(a, 'jobs', 0)
    if requested < 0: raise ValueError('--jobs must be nonnegative')
    requested = requested or int(os.environ.get('WORKFLOW_CORES') or os.environ.get('ALLELE_ALLOCATED_CORES') or os.environ.get('LSB_DJOB_NUMPROC') or 1)
    if requested < 1:
        raise ValueError('--jobs must be positive (or 0 for allocated cores)')
    # The inherited reference catalogue is shared; reserve room for private
    # Python objects, chromosome statistics, and parent-side report generation.
    memory = float(os.environ.get('ALLELE_ALLOCATED_MEMORY_GIB') or
                   os.environ.get('ALLELE_MEMORY_LIMIT_GB') or 64)
    return max(1, min(requested, allocated, int(max(1, memory - 4) // 2)))


def _worker(task):
    from research_tests import (Context, Reservoir, source_rows, observed_outcomes,
                                chrom_name, stable_hash)
    from enrichment_regions import LocusSweep
    chrom, phase = task
    a, regions, masks, template, background, total, per_chrom = _SHARED
    a = copy.copy(a)
    a.chr = chrom
    context = Context(None if a.no_context else a.fasta)
    selected = a.outcomes.split(',')
    examined = 0
    started = time.monotonic()
    print(f'[enrich] pass {phase} {chrom} START pid={os.getpid()}', flush=True)
    zero = np.zeros(2)
    local_total = defaultdict(lambda: np.zeros(2))
    status = Counter()
    observed = defaultdict(float); expected = defaultdict(float)
    denom = Counter(); unsupported = Counter()
    block = defaultdict(lambda: defaultdict(lambda: np.zeros(2)))
    cases = Reservoir(a.profile_sites); controls = Reservoir(a.controls_per_stratum)
    collector = None
    if phase == 2:
        collector = copy.copy(template)
        collector.regions = {rid:r for rid,r in template.regions.items() if r.chrom == chrom}
        collector.counts = {k:list(v) for k,v in template.counts.items() if k[0] in collector.regions}
        collector.count_modes = set()
        collector.sweep = LocusSweep(collector.regions)
    try:
        for row in source_rows(a):
            examined += 1
            if examined % 1000000 == 0:
                print(f'[enrich] pass {phase} {chrom} rows={examined:,} elapsed={time.monotonic()-started:.0f}s', flush=True)
            state, qualified = observed_outcomes(row, a)
            c = chrom_name(row['chrom']); pos = int(row['pos'])
            labels = regions(c, pos)
            if labels is None:
                if phase == 1: status['outside_reference_scoring_support'] += 1
                continue
            sens = ['all']
            if a.exclude_bed and not masks.labels(c, pos-1): sens.append('exclude_mask')
            key = context.key(row)
            if phase == 1:
                for outcome in background:
                    eligible, case = state[outcome]
                    if not eligible: continue
                    status[outcome+'_eligible'] += 1; status[outcome+'_case'] += int(case)
                    for s in sens: local_total[(outcome,s,key)] += (1,int(case))
                continue
            if state['multi_ALT'][0]:
                for s in sens:
                    other = total.get(('multi_ALT',s,key), zero) - per_chrom.get((c,'multi_ALT',s,key), zero)
                    rate = (other[1]+.5)/(other[0]+1.) if other[0]>=a.min_background_sites else None
                    collector.record(c,pos,state['multi_ALT'][1],rate,s,prior_only=bool(other[1]==0),count_mode=row['count_mode'])
            for outcome in selected:
                eligible, case = state[outcome]
                if not eligible: continue
                for s in sens:
                    other = total.get((outcome,s,key), zero) - per_chrom.get((c,outcome,s,key), zero)
                    for label in labels:
                        k = (outcome,s,label)
                        if other[0]<a.min_background_sites: unsupported[k] += 1; continue
                        rate = (other[1]+.5)/(other[0]+1.)
                        denom[k] += 1; observed[k] += int(case); expected[k] += rate
                        block[k][(c,(pos-1)//a.block_bp)] += (int(case),rate)
            if a.ref == 'alphagenome' and a.profile_sites>0 and state['multi_ALT'][0] and ('exclude_mask' in sens or not a.exclude_bed):
                rank=stable_hash(a.seed,c,pos,row['ref']); record=dict(chrom=c,pos=pos,ref=row['ref'],key=key)
                if state['multi_ALT'][1]: cases.add('all',rank,record)
                else: controls.add((c,key),rank,record)
    finally:
        context.close()
    print(f'[enrich] pass {phase} {chrom} DONE rows={examined:,} elapsed={time.monotonic()-started:.1f}s', flush=True)
    if phase == 1: return chrom, examined, dict(local_total), status
    return (dict(observed),dict(expected),denom,unsupported,{k:dict(v) for k,v in block.items()},
            dict(cases.data),dict(controls.data),collector.regions,collector.counts,collector.count_modes)


def run(a, regions, masks, collector, background, workers):
    from research_tests import chromosomes, Reservoir
    global _SHARED
    chroms = chromosomes(a.chr)
    workers = min(workers, len(chroms))
    if 'fork' not in mp.get_all_start_methods():
        raise RuntimeError('Parallel enrichment requires fork; use --jobs 1 on this platform')
    total = defaultdict(lambda: np.zeros(2)); per_chrom = {}
    status = Counter({k+suffix:0 for k in a.outcomes.split(',') for suffix in ('_eligible','_case')})
    examined = 0
    context = mp.get_context('fork')
    _SHARED = (a,regions,masks,collector,background,None,None)
    try:
        print(f'[enrich] pass 1: {workers} workers / {len(chroms)} chromosomes; global background', flush=True)
        with context.Pool(workers) as pool:
            for chrom,n,counts,counts_status in pool.imap(_worker,[(c,1) for c in chroms]):
                examined += n; status.update(counts_status)
                for k,v in counts.items(): total[k] += v; per_chrom[(chrom,*k)] = v
        _SHARED = (a,regions,masks,collector,background,dict(total),per_chrom)
        observed=defaultdict(float); expected=defaultdict(float); denom=Counter(); unsupported=Counter()
        block=defaultdict(dict); cases=Reservoir(a.profile_sites); controls=Reservoir(a.controls_per_stratum)
        print(f'[enrich] pass 2: {workers} workers; shared leave-one-chromosome-out background', flush=True)
        with context.Pool(workers) as pool:
            for obs,exp,den,uns,blocks,case_data,control_data,rs,counts,modes in pool.imap(_worker,[(c,2) for c in chroms]):
                for k,v in obs.items(): observed[k] += v
                for k,v in exp.items(): expected[k] += v
                denom.update(den); unsupported.update(uns)
                for k,v in blocks.items(): block[k].update(v)
                for target,data in ((cases,case_data),(controls,control_data)):
                    for k,items in data.items():
                        for neg_rank,serial,value in sorted(items,key=lambda x:x[1]): target.add(k,-neg_rank,value)
                collector.regions.update(rs); collector.counts.update(counts); collector.count_modes.update(modes)
        return examined,status,observed,expected,denom,unsupported,block,cases,controls
    finally:
        _SHARED = None
