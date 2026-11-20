"""Parallel reduction of existing chromosome counters for the primary estimator."""
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import csv
import gzip
import json
import os
from scan import bootstrap_statistics, bh_fdr

DIMS = ('context','af_bin','an_bin','sensitivity','population_mode','threshold')

def workers_for(n):
    limits = [n, len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)]
    for key in ('ALLELE_ALLOCATED_CORES', 'LSB_DJOB_NUMPROC', 'WORKFLOW_CORES'):
        if os.environ.get(key):
            limits.append(max(1, int(os.environ[key])))
    return max(1, min(limits))

def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as f:
        for line, row in enumerate(csv.DictReader(f, delimiter='\t'), 2):
            if None in row or any(value is None or value == '' for value in row.values()):
                raise ValueError(f'Incomplete counter row in {path}:{line}')
            try:
                n, m = int(row['n_sites']), int(row['n_multi'])
                if not 0 <= m <= n:
                    raise ValueError('require 0 <= n_multi <= n_sites')
            except (KeyError, ValueError) as exc:
                raise ValueError(f'Invalid counter counts in {path}:{line}: {exc}') from exc
            yield row

def global_counts(path):
    total = defaultdict(lambda: [0, 0])
    chroms = defaultdict(lambda: [0, 0])
    for r in read(path):
        if r['threshold'] not in {'5','10'} or r['population_mode'] != 'all':
            continue
        st = tuple(r[x] for x in DIMS)
        for table, key in ((total, st), (chroms, (r['chrom'], *st))):
            table[key][0] += int(r['n_sites'])
            table[key][1] += int(r['n_multi'])
    print('[primary global complete] ' + str(path), flush=True)
    return dict(total), dict(chroms)

_BACKGROUND = None

def init_background(total, chroms, minimum, exclude_ibdmix):
    global _BACKGROUND
    _BACKGROUND = total, chroms, minimum, exclude_ibdmix

def exposure_counts(path):
    total, chroms, minimum, exclude_ibdmix = _BACKGROUND
    blocks = defaultdict(lambda: defaultdict(lambda: [0, 0., 0, 0]))
    outside = Counter()
    excluded = 0
    for r in read(path):
        if exclude_ibdmix and '::archaic:' in r.get('label', ''):
            from herv import modern_archaic_label_reference
            if modern_archaic_label_reference(r['label']) is not None:
                excluded += 1
                continue
        if r['threshold'] not in {'5','10'} or r['population_mode'] != 'all':
            continue
        key = (r['label'], r['threshold'], r['sensitivity'])
        if r['block_id'] == 'outside_partition':
            outside[key] += int(r['n_sites'])
            continue
        st = tuple(r[x] for x in DIMS)
        n, observed = int(r['n_sites']), int(r['n_multi'])
        t, c = total.get(st, (0, 0)), chroms.get((r['chrom'], *st), (0, 0))
        tn, tm = t[0]-c[0], t[1]-c[1]
        z = blocks[key][r['block_id']]
        z[0] += observed
        z[2] += n
        if tn < minimum:
            z[3] += n
        else:
            z[1] += n*(tm+.5)/(tn+1)
    print('[primary exposure complete] ' + str(path), flush=True)
    return {k: dict(v) for k, v in blocks.items()}, outside, excluded

def summarize(global_paths, exposure_paths, out, reps=2000, seed=20261024,
              min_training_sites=1, exclude_ibdmix=False):
    if reps < 2000:
        raise ValueError('At least 2000 bootstrap replicates required')
    for kind, paths in (('global', global_paths), ('exposure', exposure_paths)):
        if not paths or len({Path(p).resolve() for p in paths}) != len(paths):
            raise ValueError('Counter paths must be nonempty and unique')
        header = None
        for path in paths:
            opener = gzip.open if str(path).endswith('.gz') else open
            with opener(path, 'rt', newline='') as handle:
                fields = next(csv.reader(handle, delimiter='\t'), None)
            if not fields or (header is not None and fields != header):
                raise ValueError('Counter tables require the same nonempty schema')
            required = {'chrom', 'n_sites', 'n_multi', *DIMS}
            if kind == 'exposure':
                required.update(('label', 'block_id'))
            missing = required - set(fields)
            if missing:
                raise ValueError(f'{kind} counter {path} is missing required fields: {sorted(missing)}')
            header = fields
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    total = defaultdict(lambda: [0, 0])
    chroms = defaultdict(lambda: [0, 0])
    workers = workers_for(max(len(global_paths), len(exposure_paths)))
    print(f'[primary parallel] workers={workers}; global_files={len(global_paths)}; exposure_files={len(exposure_paths)}', flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for part_total, part_chroms in pool.map(global_counts, global_paths):
            for target, part in ((total, part_total), (chroms, part_chroms)):
                for key, values in part.items():
                    target[key][0] += values[0]
                    target[key][1] += values[1]
    blocks = defaultdict(lambda: defaultdict(lambda: [0, 0., 0, 0]))
    outside = Counter()
    excluded = 0
    with ProcessPoolExecutor(max_workers=workers, initializer=init_background,
            initargs=(dict(total), dict(chroms), min_training_sites, exclude_ibdmix)) as pool:
        # Ordered reduction preserves source and block order for seeded bootstrap.
        for part, part_outside, n_excluded in pool.map(exposure_counts, exposure_paths):
            excluded += n_excluded
            outside.update(part_outside)
            for key, values in part.items():
                for bid, value in values.items():
                    if bid not in blocks[key]:
                        blocks[key][bid] = value
                    else:
                        for i in range(4):
                            blocks[key][bid][i] += value[i]
    print('[primary reduction complete] bootstrap and output', flush=True)
    fields=['label','threshold','sensitivity','observed_multi','expected_multi','oe','ci_low','ci_high','bootstrap_p','fdr','block_count','n_eligible','n_without_loco','outside_partition_sites','status']
    rows=[]
    with gzip.open(out/'enrichment_block_counts.tsv.gz','wt',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['label','threshold','sensitivity','block_id','observed_multi','expected_multi','n_eligible','n_without_loco'])
        for key,b in sorted(blocks.items()):
            # Stable bootstrap input order for either merged tables or shards.
            b = dict(sorted(b.items()))
            for bid,vals in sorted(b.items()):w.writerow([*key,bid,*vals])
            o=sum(x[0] for x in b.values());e=sum(x[1] for x in b.values());n=sum(x[2] for x in b.values());missing=sum(x[3] for x in b.values())
            ci=(None,None,None) if missing else bootstrap_statistics({bid:vals[:2] for bid,vals in b.items()},reps,seed)
            status='insufficient_loco_training' if missing else 'insufficient_blocks' if len(b)<2 else 'zero_expected' if not e else 'complete'
            rows.append(dict(zip(('label','threshold','sensitivity'),key))|dict(observed_multi=o,expected_multi=e if not missing else 'NA',oe=o/e if e and not missing else 'NA',ci_low=ci[0],ci_high=ci[1],bootstrap_p=ci[2],fdr=None,block_count=len(b),n_eligible=n,n_without_loco=missing,outside_partition_sites=outside[key],status=status))
    groups=defaultdict(list)
    for r in rows:
        if r['bootstrap_p'] is not None:groups[r['threshold'],r['sensitivity']].append(r)
    for group in groups.values():
        for r,q in zip(group,bh_fdr([r['bootstrap_p'] for r in group])):r['fdr']=q
    with open(out/'enrichment_summary.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
    manifest=dict(status='complete',bootstrap_replicates=reps,seed=seed,primary_thresholds=[5,10],ld_ancestry='EUR',inference='LOCO context/primary-AF/AN expected; percentile LD-block bootstrap CI; centered two-sided bootstrap P; BH FDR within threshold/sensitivity',missing_loco='NA if any exposed stratum is unsupported',outside_partition='excluded from inference; count retained',n_tests=len(rows),difficult_region_complete=any(r['sensitivity']=='exclude_combined' for r in rows))
    (out/'enrichment_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest, excluded

