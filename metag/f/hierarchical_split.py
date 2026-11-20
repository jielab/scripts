"""Disk-backed group routing, selected species refinement and automatic evaluation.

Stage 1 competes all reference categories in one compact Kraken index. Stage 2
queries only the detailed category needed for a selected group. The latter is
conditional evidence, never a replacement for all-category competition.
Original paired FASTQ bytes and original ordinals are retained at both stages.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import gzip
import hashlib
import html
import itertools
import json
import os
from pathlib import Path
import random
import re
import shutil
import struct
import sys
import tempfile
import time

from common import (atomic_json, code_identity, file_state, fingerprint,
                    load_config, locked, require_compute, utcnow)
from compressed_output import BufferedWriters, gzip_text
from fastq_stream import fragments, new_digest, read_id, update_digest
from ordinal_query import input_states, query_once, group_query_config
from query import discover_samples
from result_output import replace_bool
from species_split import taxids
from storage_layout import library_root
from taxonomy import GROUPS, Taxonomy
from project_cache import routing_taxonomy

SCHEMA = 'metag-hierarchy-v1'
FIELDS = ('ordinal', 'read_id', 'r1_read_id', 'r2_read_id', 'group', 'taxid',
          'species_taxid', 'status', 'bin', 'source_group', 'candidate_taxid')
BIO_GROUPS = tuple(g for g in GROUPS if g != 'unresolved')
ROUTES = {'bacteria':'microbes', 'archaea':'microbes', 'fungi':'microbes',
          'virus':'viruses', 'plant':'plants', 'human':'vertebrates',
          'mammal':'vertebrates', 'other_vertebrate':'vertebrates',
          'invertebrate':'invertebrates'}
UNCERTAIN = ('ambiguous', 'unclassified')


def state(path):
    p = Path(path).resolve(); s = p.stat()
    return dict(path=str(p), size=s.st_size, mtime_ns=s.st_mtime_ns,
                ctime_ns=s.st_ctime_ns)


def read_json(path):
    return json.loads(Path(path).read_text())


def write_table(path, fields, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def canonical_id(header):
    name = read_id(header).decode('utf-8', 'strict')
    if name.endswith(('/1', '/2')):
        return name[:-2]
    m = re.fullmatch(r'([SED]RR\d+\.\d+)\.[12]', name)
    return m[1] if m else name


def biological_group(info):
    if not info['taxid']:
        return 'unclassified'
    return 'ambiguous' if info['group'] == 'unresolved' else info['group']


def qc_folder(sample_root, stage):
    return Path(sample_root)/'group'/('qc' if stage == 'group' else 'species.qc')


class Reservoir:
    """Uniform sampling within declared strata, independent of file prefix."""
    def __init__(self, cfg):
        self.seed = int(cfg.get('evaluation_seed', 20261115))
        self.general = int(cfg.get('evaluation_pairs_per_stratum', 1000))
        self.target = int(cfg.get('evaluation_target_pairs_per_stratum', 10000))
        if min(self.general, self.target) < 1:
            raise ValueError('Evaluation sample sizes must be positive')
        self.counts = Counter(); self.rows = defaultdict(list); self.rng = {}

    def add(self, row, a, b):
        key = str(row['group']) + ':' + str(row['status'])
        limit = self.target if row['group'] in ('human', 'plant') else self.general
        if key not in self.rng:
            value = int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], 'little')
            self.rng[key] = random.Random(self.seed ^ value)
        self.counts[key] += 1; n = self.counts[key]
        item = (dict(row), a, b)
        if len(self.rows[key]) < limit:
            self.rows[key].append(item)
        else:
            j = self.rng[key].randrange(n)
            if j < limit:
                self.rows[key][j] = item

    def publish(self, qc):
        qc = Path(qc); qc.mkdir(parents=True, exist_ok=True)
        selected = []
        for key, rows in self.rows.items():
            for row, a, b in rows:
                row.update(stratum=key, stratum_pairs=self.counts[key],
                           sampled_pairs=len(rows))
                selected.append((row, a, b))
        selected.sort(key=lambda x:int(x[0]['ordinal']))
        with gzip.open(qc/'sample.R1.fastq.gz', 'wb', compresslevel=3) as aout, \
             gzip.open(qc/'sample.R2.fastq.gz', 'wb', compresslevel=3) as bout:
            for _, a, b in selected:
                aout.write(b''.join(a)); bout.write(b''.join(b))
        write_table(qc/'sample.tsv', (*FIELDS, 'stratum', 'stratum_pairs', 'sampled_pairs'),
                    [x[0] for x in selected])
        return dict(method='uniform_reservoir_within_group_and_status', seed=self.seed,
                    pairs=len(selected), strata=dict(self.counts),
                    target_limit=self.target, other_limit=self.general)


def stage_request(cfg, sample, stage, extra):
    request=dict(schema=SCHEMA, stage=stage, code=code_identity(),
                sample=sample['sample_id'], inputs=[state(sample[k]) for k in ('r1','r2')],
                gzip_level=cfg.get('split_gzip_level', 3),
                minimum_base_quality=cfg.get('minimum_base_quality',0),
                evaluation_seed=cfg.get('evaluation_seed', 20261115),
                evaluation_pairs_per_stratum=cfg.get('evaluation_pairs_per_stratum', 1000),
                evaluation_target_pairs_per_stratum=cfg.get('evaluation_target_pairs_per_stratum', 10000))
    request.update(extra)
    return request


def reuse_stage(sample_root, stage, request, replace):
    qc = qc_folder(sample_root, stage)
    marker = qc/'summary.json'
    if not marker.is_file():
        return False
    previous = read_json(marker)
    if replace:
        return False
    from cache_compatibility import compatible_request
    if not compatible_request(previous.get('request'), request) or previous.get('status') != 'complete':
        raise RuntimeError(f'{stage} results differ from this request; use --replace TRUE to recompute this stage')
    for path, size in previous.get('artifacts', {}).items():
        p = Path(sample_root)/path
        if not p.is_file() or p.stat().st_size != size:
            raise RuntimeError(f'Incomplete stage artifact: {p}; use --replace TRUE')
        expected=previous.get('artifact_observations',{}).get(path)
        if expected and {k:v for k,v in state(p).items() if k!='path'}!=expected:
            raise RuntimeError(f'Stage artifact changed: {p}; use --replace TRUE')
    print(f'REUSE {stage}: {sample_root}', flush=True)
    return True


def artifact_sizes(root, paths):
    root=Path(root)
    return {str(p.relative_to(root)):p.stat().st_size for p in paths if p.is_file()}


def finish_counts(qc, physical_root, bins, taxa, prefix, population, request, sampling,
                  queries, input_digest, extra=None):
    rows=[]
    for name, n in sorted(bins.items()):
        r1=physical_root/name/'R1.fastq.gz'; r2=physical_root/name/'R2.fastq.gz'
        rows.append(dict(bin=f'{prefix}/{name}', pairs=n,
                         r1=f'{prefix}/{name}/R1.fastq.gz', r2=f'{prefix}/{name}/R2.fastq.gz',
                         bytes_r1=r1.stat().st_size, bytes_r2=r2.stat().st_size))
    write_table(qc/'counts.tsv', ('bin','pairs','r1','r2','bytes_r1','bytes_r2'), rows)
    write_table(qc/'taxa.tsv', ('group','taxid','species_taxid','status','pairs'),
                [dict(zip(('group','taxid','species_taxid','status','pairs'), (*k,n)))
                 for k,n in sorted(taxa.items())])
    count=sum(bins.values())
    report=dict(schema=SCHEMA, status='complete', stage=request['stage'],
                sample_id=request['sample'], request=request, paired=True, input_pairs=count, written_pairs=count,
                assignment_rows=count, expected_pairs=population, input_digest=input_digest,
                sampling=sampling, queries=queries, completed_at=utcnow(),
                evaluation_status='pending', **(extra or {}))
    atomic_json(qc/'summary.json', report)
    return report


def _publish_directory(stage, destination, replace):
    """Publish within one filesystem; replace only an explicitly owned stage."""
    stage, destination = Path(stage), Path(destination)
    backup = destination.with_name(destination.name+'.previous')
    if destination.exists():
        if not replace:
            raise RuntimeError(f'Output already exists: {destination}')
        if not (destination/'qc/summary.json').exists() and not (destination/'summary.json').exists():
            raise RuntimeError(f'Refusing to replace unowned output: {destination}')
        if backup.exists():
            raise RuntimeError(f'Previous interrupted publication needs inspection: {backup}')
        destination.rename(backup)
    try:
        stage.rename(destination)
    except BaseException:
        if backup.exists() and not destination.exists(): backup.rename(destination)
        raise
    if backup.exists(): shutil.rmtree(backup)


def split_group(cfg, sample, *, replace=False):
    root=Path(cfg['output_dir'])/sample['sample_id']; root.mkdir(parents=True, exist_ok=True)
    previous=root/'group.previous'
    if previous.exists():
        if not (previous/'qc/summary.json').is_file():raise RuntimeError('Unowned group.previous requires inspection')
        if not (root/'group').exists():previous.rename(root/'group')
        elif (root/'group/qc/summary.json').is_file():shutil.rmtree(previous)
        else:raise RuntimeError('Interrupted group publication requires inspection')
    from group_library import compact_snapshot
    if cfg.get('group_routing', 'compact') == 'full':
        from partitioned_build import database_snapshot
        snapshot=database_snapshot(cfg)
        method='full_competitive_kraken'
    else:
        snapshot=compact_snapshot(cfg)
        method='compact_competitive_kraken'
    extra=dict(method=method, database=snapshot['identity'],
               confidence=float(cfg.get('group_confidence',0.0)),
               minimum_hit_groups=int(cfg.get('group_minimum_hit_groups',2)))
    taxonomy_dir=snapshot['taxonomy_dir']
    tax=routing_taxonomy(cfg,taxonomy_dir)
    request=stage_request(cfg,sample,'group',extra)
    if reuse_stage(root,'group',request,replace): return root
    destination=root/'group'
    if destination.exists() and not replace:
        raise RuntimeError('Incomplete/unowned group directory exists; inspect it or use --replace TRUE for an owned completed stage')
    stage=Path(tempfile.mkdtemp(prefix='.group-', dir=root)); qc=stage/'qc';qc.mkdir()
    work=Path(cfg['work_dir'])/'tmp'/'hierarchy'/sample['sample_id']
    work.mkdir(parents=True, exist_ok=True)
    query_work=Path(tempfile.mkdtemp(prefix='group-',dir=work))
    try:
        queries=[]
        qcfg=dict(group_query_config(cfg, snapshot),
                  confidence=extra['confidence'],minimum_hit_groups=extra['minimum_hit_groups'])
        q=query_once(qcfg,sample,snapshot,query_work,taxonomy=tax);queries=[q]
        calls=iter(taxids(q['taxids'],q['fragments']))
        expected=q['fragments']; expected_digest=q['input_digest']
        bins=Counter();taxa=Counter(); reservoir=Reservoir(cfg); digest=new_digest(); n=0
        group_digests=defaultdict(new_digest); ordinal_digests=defaultdict(hashlib.sha256)
        input_before=input_states(sample)
        writers=BufferedWriters(stage,threads=max(1,int(cfg.get('classify_threads',4))-3),
                                compresslevel=int(cfg.get('split_gzip_level',3)))
        with ExitStack() as stack:
            out=stack.enter_context(gzip_text(qc/'assignments.tsv.gz',threads=min(2,max(0,int(cfg.get('classify_threads',4))-2))))
            aw=csv.DictWriter(out,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');aw.writeheader()
            sidecars={}
            try:
                source=fragments(sample['r1'],sample['r2'])
                for n,(a,b) in enumerate(source,1):
                    call=next(calls,None)
                    if call is None: raise RuntimeError('More FASTQ pairs than assignments')
                    tid=call
                    info=tax.describe(tid); group=biological_group(info)
                    row=dict(ordinal=n,read_id=canonical_id(a[0]),r1_read_id=read_id(a[0]).decode(),
                             r2_read_id=read_id(b[0]).decode(),group=group,taxid=info['taxid'],
                             species_taxid=info['species_taxid'],status=info['status'],
                             bin=f'group/{group}',source_group='',candidate_taxid=info['taxid'])
                    aw.writerow(row); bins[group]+=1
                    taxa[(group,row['taxid'],row['species_taxid'],row['status'])]+=1
                    for mate,rec in ((1,a),(2,b)):writers.write(f'{group}/R{mate}.fastq.gz',rec)
                    if group not in sidecars:
                        (stage/group).mkdir(parents=True,exist_ok=True)
                        sidecars[group]=stack.enter_context(gzip.open(stage/group/'ordinals.tsv.gz','wt',compresslevel=3))
                        sidecars[group].write('ordinal\ttaxid\n')
                    ordinal_line=f'{n}\t{row["taxid"]}\n'
                    sidecars[group].write(ordinal_line);ordinal_digests[group].update(ordinal_line.encode())
                    update_digest(group_digests[group],a,b)
                    reservoir.add(row,a,b);update_digest(digest,a,b)
                    if n%1000000==0:print(f'group {sample["sample_id"]}: wrote {n:,} pairs',flush=True)
                if next(calls,None) is not None or n!=expected:raise RuntimeError('FASTQ/assignment count mismatch')
                if digest.hexdigest()!=expected_digest:raise RuntimeError('Original FASTQ stream identity differs')
                if input_before!=input_states(sample):raise RuntimeError('Original input changed during routing')
                writers.close()
            except BaseException:
                writers.close(flush=False);raise
        sampling=reservoir.publish(qc)
        report=finish_counts(qc,stage,bins,taxa,'group',n,request,sampling,queries,digest.hexdigest(),
                             dict(group_counts=dict(bins),scope='all_original_pairs',taxonomy_dir=taxonomy_dir,
                                  group_input_digests={k:v.hexdigest() for k,v in group_digests.items()},
                                  group_ordinal_digests={k:v.hexdigest() for k,v in ordinal_digests.items()},
                                  note='Stored species labels are evidence only; no species FASTQs are materialized in this stage.'))
        # Paths recorded relative to the final sample root; stage is atomic.
        report['artifacts']={str(Path('group')/p.relative_to(stage)):p.stat().st_size
                             for p in stage.rglob('*') if p.is_file() and p!=qc/'summary.json'}
        report['artifact_observations']={str(Path('group')/p.relative_to(stage)):
            {k:v for k,v in state(p).items() if k!='path'} for p in stage.rglob('*')
            if p.is_file() and p!=qc/'summary.json'}
        atomic_json(qc/'summary.json',report)
        _publish_directory(stage,destination,replace)
        return root
    finally:
        if stage.exists():shutil.rmtree(stage)
        # On failure keep native log/taxids for diagnosis; on success remove only taxids.
        if (destination/'qc/summary.json').is_file():
            for p in query_work.glob('*.bin'):p.unlink(missing_ok=True)


def category_snapshot(cfg, category):
    """Use an existing category as an ordinary single index, with the global taxonomy."""
    from partitioned_build import database_snapshot
    full=database_snapshot(cfg)
    db=Path(full['db_dir'])/category
    if category not in full['context']['categories']:
        raise RuntimeError(f'No existing detailed category index: {category}')
    files={name:dict(state=file_state(db/name)) for name in ('hash.k2d','taxo.k2d','opts.k2d')}
    # A category has its own published database.json; reuse its real identity,
    # because the ordinary native-query verifier also checks that marker.
    marker=read_json(db/'database.json')
    if marker.get('status')!='complete' or marker.get('identity')!=fingerprint(marker.get('context')):
        raise RuntimeError('Incomplete detailed category: '+category)
    identity=marker['identity']
    return dict(status='complete',identity=identity,layout='single',db_dir=str(db),
                taxonomy_dir=full['taxonomy_dir'],files=files,taxonomy_files=full['taxonomy_files'],
                parent_identity=full['identity'],category=category,
                scope='conditional_on_coarse_group; competitors_inside_this_category_only')


def species_route(info, source_group, original_taxid, counts, minimum):
    actual=biological_group(info)
    if source_group in UNCERTAIN:
        if info['species_taxid']:
            tid=info['species_taxid'];leaf=str(tid) if tid==9606 or counts.get(tid,0)>=minimum else 'other'
            return dict(group=actual,taxid=info['taxid'],species_taxid=tid,status='assigned_species',leaf=leaf)
        return dict(group=actual,taxid=info['taxid'],species_taxid=0,
                    status='taxon_unresolved' if info['taxid'] else 'unclassified',leaf='unresolved')
    if info['taxid'] and actual not in (source_group,'ambiguous','unclassified'):
        return dict(group='ambiguous',taxid=0,species_taxid=0,status='group_conflict',leaf='ambiguous')
    if info['species_taxid'] and actual==source_group:
        tid=info['species_taxid']
        leaf=str(tid) if tid==9606 or counts.get(tid,0)>=minimum else 'other'
        return dict(group=source_group,taxid=info['taxid'],species_taxid=tid,status='assigned_species',leaf=leaf)
    if info['taxid'] and actual==source_group:
        return dict(group=source_group,taxid=info['taxid'],species_taxid=0,status='taxon_unresolved',leaf='unresolved')
    # A broad ancestor is not a contrary group. Retain the original coarse claim,
    # while marking the detailed stage as unresolved; do not manufacture a species.
    return dict(group=source_group,taxid=int(original_taxid),species_taxid=0,
                status='unclassified_within_group' if not info['taxid'] else 'taxon_unresolved',leaf='unresolved')


def recover_species_publication(parent):
    """Recover only this workflow's journalled, same-filesystem publication."""
    parent=Path(parent)
    for backup in parent.glob('.species-backup-*'):
        journal=read_json(backup/'transaction.json')
        if journal.get('schema')!=SCHEMA:raise RuntimeError('Unowned species transaction: '+str(backup))
        final=parent/'species.qc/summary.json'
        if final.is_file() and fingerprint(read_json(final).get('request'))==journal['request_identity']:
            shutil.rmtree(backup);continue
        prior=set(journal['prior']);chosen=set(journal['chosen'])
        for group in sorted(prior|chosen):
            saved=backup/group;dest=parent/group/'species'
            if saved.exists():
                if dest.exists():shutil.rmtree(dest)
                saved.rename(dest)
            elif group not in prior and dest.exists():
                shutil.rmtree(dest)
        saved=backup/'species.qc'
        if saved.exists():
            if (parent/'species.qc').exists():shutil.rmtree(parent/'species.qc')
            saved.rename(parent/'species.qc')
        shutil.rmtree(backup)


def publish_species(stage,parent,chosen,prior,request,replace):
    """Preflight every destination and roll back a partial multi-child update."""
    stage,parent=Path(stage),Path(parent)
    for group in set(prior)|set(chosen):
        dest=parent/group/'species'
        if dest.is_symlink():raise RuntimeError('Refusing symlink species destination: '+str(dest))
        if dest.exists() and (not replace or group not in prior):
            raise RuntimeError('Refusing to overwrite unowned species output: '+str(dest))
    oldqc=parent/'species.qc'
    if oldqc.exists() and (not replace or not (oldqc/'summary.json').is_file()):
        raise RuntimeError('Refusing to overwrite unowned species metadata')
    backup=Path(tempfile.mkdtemp(prefix='.species-backup-',dir=parent))
    atomic_json(backup/'transaction.json',dict(schema=SCHEMA,prior=sorted(prior),chosen=sorted(chosen),
                request_identity=fingerprint(request)))
    try:
        if oldqc.exists():oldqc.rename(backup/'species.qc')
        for group in sorted(prior):
            dest=parent/group/'species'
            if dest.exists():dest.rename(backup/group)
        for group in chosen:(stage/group/'species').rename(parent/group/'species')
        (stage/'species.qc').rename(oldqc)
    except BaseException:
        recover_species_publication(parent);raise
    shutil.rmtree(backup)


def split_species(cfg,sample,*,groups=None,include_uncertain=False,replace=False):
    root=Path(cfg['output_dir'])/sample['sample_id'];parent=root/'group'
    recover_species_publication(parent)
    source=read_json(parent/'qc/summary.json')
    if source.get('status')!='complete':raise RuntimeError('Run split-group successfully first')
    groups=groups or list(cfg.get('species_groups',['human','plant']))
    if groups==['all']:groups=[g for g in BIO_GROUPS if source.get('group_counts',{}).get(g,0)]
    if len(set(groups))!=len(groups) or any(g not in BIO_GROUPS for g in groups):
        raise ValueError('Choose biological groups from: '+','.join(BIO_GROUPS))
    chosen=[g for g in groups if source.get('group_counts',{}).get(g,0)]
    if include_uncertain:chosen += [g for g in UNCERTAIN if source.get('group_counts',{}).get(g,0)]
    group_inputs={}
    for group in chosen:
        for name in ('R1.fastq.gz','R2.fastq.gz','ordinals.tsv.gz'):
            relative=f'group/{group}/{name}';p=root/relative
            if not p.is_file() or p.stat().st_size!=source['artifacts'].get(relative):
                raise RuntimeError('Coarse input no longer matches its completion record: '+str(p))
            group_inputs[relative]=state(p)
    tax=routing_taxonomy(cfg,source.get('taxonomy_dir',cfg['taxonomy_dir']))
    snapshots={}
    for group in chosen:
        if group=='other_eukaryote' or group in UNCERTAIN:
            from partitioned_build import database_snapshot
            snapshots[group]=database_snapshot(cfg)
        else:snapshots[group]=category_snapshot(cfg,ROUTES[group])
    request=stage_request(cfg,sample,'species',dict(
        method='detailed_category_kraken',groups=groups,parent=fingerprint(source['request']),
        include_uncertain=include_uncertain,
        confidence=cfg['confidence'],minimum_hit_groups=cfg['minimum_hit_groups'],
        minimum_species_file_pairs=cfg.get('materialize_species_min_fragments',20),
        databases={g:s['identity'] for g,s in snapshots.items()},group_inputs=group_inputs))
    if reuse_stage(root,'species',request,replace):return root
    if (parent/'species.qc').exists() and not replace:raise RuntimeError('Species result exists; use --replace TRUE')
    stage=Path(tempfile.mkdtemp(prefix='.species-',dir=parent));qc=stage/'species.qc';qc.mkdir()
    work=Path(cfg['work_dir'])/'tmp'/'hierarchy'/sample['sample_id'];work.mkdir(parents=True,exist_ok=True)
    reservoir=Reservoir(cfg);bins=Counter();taxa=Counter();queries=[];total=0
    all_digest=new_digest();minimum=int(cfg.get('materialize_species_min_fragments',20))
    try:
        with gzip_text(qc/'assignments.tsv.gz',threads=min(2,max(0,int(cfg['classify_threads'])-2))) as assignments:
            aw=csv.DictWriter(assignments,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');aw.writeheader()
            for group in chosen:
                folder=parent/group
                current=dict(sample_id=sample['sample_id']+'/'+group,r1=str(folder/'R1.fastq.gz'),r2=str(folder/'R2.fastq.gz'))
                qwork=Path(tempfile.mkdtemp(prefix=group+'-',dir=work))
                snapshot=snapshots[group]
                qcfg=dict(cfg,database_layout=snapshot['layout'],union_memory_mapping=True)
                q=query_once(qcfg,current,snapshot,qwork,taxonomy=tax);queries.append(q)
                expected=int(source['group_counts'][group])
                if q['fragments']!=expected:raise RuntimeError('Coarse bin count changed: '+group)
                if source.get('group_input_digests',{}).get(group)!=q['input_digest']:
                    raise RuntimeError('Coarse input content differs from the committed group stream: '+group)
                observed=Counter(tax.species(t) for t in taxids(q['taxids'],expected))
                calls=taxids(q['taxids'],expected);digest=new_digest();ordinal_digest=hashlib.sha256();n=0;previous=0
                writers=BufferedWriters(stage,threads=max(1,int(cfg['classify_threads'])-3),compresslevel=int(cfg['split_gzip_level']))
                try:
                    with gzip.open(folder/'ordinals.tsv.gz','rt',newline='') as idx:
                        indexes=csv.DictReader(idx,delimiter='\t')
                        for n,(a,b) in enumerate(fragments(current['r1'],current['r2']),1):
                            saved=next(indexes,None);tid=next(calls,None)
                            if saved is None or tid is None:raise RuntimeError('Missing group ordinal/classification')
                            original=int(saved['ordinal'])
                            if not previous<original<=int(source['expected_pairs']):raise RuntimeError('Invalid original ordinal in coarse bin')
                            previous=original
                            ordinal_digest.update(f'{original}\t{int(saved["taxid"])}\n'.encode())
                            route=species_route(tax.describe(tid),group,int(saved['taxid']),observed,minimum)
                            name=f'{group}/species/{route.pop("leaf")}'
                            row=dict(ordinal=original,read_id=canonical_id(a[0]),r1_read_id=read_id(a[0]).decode(),
                                     r2_read_id=read_id(b[0]).decode(),bin='group/'+name,source_group=group,
                                     candidate_taxid=tid,**route)
                            aw.writerow(row);bins[name]+=1
                            taxa[(row['group'],row['taxid'],row['species_taxid'],row['status'])]+=1
                            for mate,rec in ((1,a),(2,b)):writers.write(f'{name}/R{mate}.fastq.gz',rec)
                            reservoir.add(row,a,b);update_digest(digest,a,b);update_digest(all_digest,a,b);total+=1
                        if next(indexes,None) is not None or next(calls,None) is not None or n!=expected:
                            raise RuntimeError('Detailed pair/ordinal count mismatch')
                    if digest.hexdigest()!=q['input_digest']:raise RuntimeError('Coarse FASTQ changed during species routing')
                    if ordinal_digest.hexdigest()!=source.get('group_ordinal_digests',{}).get(group):
                        raise RuntimeError('Coarse original-ordinal mapping changed: '+group)
                    writers.close()
                except BaseException:
                    writers.close(flush=False);raise
                Path(q['taxids']).unlink(missing_ok=True)
        sampling=reservoir.publish(qc)
        report=finish_counts(qc,stage,bins,taxa,'group',int(source['expected_pairs']),request,sampling,queries,
                             all_digest.hexdigest(),dict(selected_groups=groups,processed_groups=chosen,
                             unprocessed_pairs=int(source['expected_pairs'])-total,
                             scope='selected_groups_and_uncertain' if include_uncertain else 'selected_coarse_groups_only',
                             taxonomy_dir=source.get('taxonomy_dir',cfg['taxonomy_dir']),
                             note='Known groups use conditional category evidence; optional uncertain inputs use full all-category competition. Unselected pairs remain in coarse bins.'))
        old_qc=parent/'species.qc'
        old=read_json(old_qc/'summary.json') if (old_qc/'summary.json').is_file() else {}
        prior=set(old.get('processed_groups',[]))
        report['artifacts']={}
        report['artifact_observations']={}
        for name in bins:
            for mate in (1,2):
                p=stage/name/f'R{mate}.fastq.gz';key=str(Path('group')/p.relative_to(stage))
                report['artifacts'][key]=p.stat().st_size
                report['artifact_observations'][key]={k:v for k,v in state(p).items() if k!='path'}
        for p in qc.rglob('*'):
            if p.is_file() and p.name!='summary.json':
                report['artifacts'][str(Path('group/species.qc')/p.relative_to(qc))]=p.stat().st_size
                report['artifact_observations'][str(Path('group/species.qc')/p.relative_to(qc))]={k:v for k,v in state(p).items() if k!='path'}
        atomic_json(qc/'summary.json',report)
        publish_species(stage,parent,chosen,prior,request,replace)
        return root
    finally:
        if stage.exists():shutil.rmtree(stage)


def classifier_agreement(sample_root,stage,cfg):
    from alignment_cache import observation, read_json as cache_json
    from partitioned_build import database_snapshot
    qc=qc_folder(sample_root,stage);out=qc/'agreement'
    receipt=out/'cache.json'
    try:
        summary=read_json(qc/'summary.json')
        full=database_snapshot(cfg)
        from common import sha256
        context=dict(schema=1,stage=stage,implementation={name:sha256(Path(__file__).with_name(name))
            for name in ('hierarchical_split.py','query.py','ordinal_query.py','taxonomy.py','common.py')},
            inputs={name:observation(qc/name) for name in ('sample.tsv','sample.R1.fastq.gz','sample.R2.fastq.gz')},
            database=full['identity'],seed=cfg.get('evaluation_seed',20261115),
            pairs=cfg.get('evaluation_agreement_pairs_per_stratum',200),
            parameters={k:summary['request'][k] for k in ('confidence','minimum_hit_groups','minimum_base_quality')})
        key=fingerprint(context)
        saved=cache_json(receipt)
        valid=False
        try:
            valid=(saved.get('key')==key and bool(saved.get('files')) and all(
                observation(out/name)==stamp for name,stamp in saved['files'].items()))
        except OSError:
            pass
        if valid:
            result=read_json(out/'summary.json')
            if result.get('status')=='complete':
                print(f'CACHE agreement HIT: {sample_root} {stage}',flush=True)
                return result
    except (OSError,ValueError,KeyError,RuntimeError):
        context=None
    if context is not None and not receipt.exists():
        from cache_compatibility import legacy_agreement_valid
        previous=cache_json(out/'summary.json')
        if previous.get('status')=='complete' and legacy_agreement_valid(qc,summary,context,previous,code_identity()):
            atomic_json(receipt,dict(key=key,files={name:observation(out/name)
                for name in ('summary.json','pairs.tsv','transitions.tsv')}))
            print(f'CACHE agreement VERIFIED LEGACY HIT: {sample_root} {stage}',flush=True)
            return previous
    print(f'CACHE agreement BUILD: {sample_root} {stage}',flush=True)
    result=_classifier_agreement_uncached(sample_root,stage,cfg)
    if context is not None and result.get('status')=='complete':
        current={name:observation(qc/name) for name in context['inputs']}
        if current==context['inputs']:
            atomic_json(receipt,dict(key=key,files={name:observation(out/name)
                for name in ('summary.json','pairs.tsv','transitions.tsv')}))
    return result


def _classifier_agreement_uncached(sample_root,stage,cfg):
    """A bounded larger-library comparison; agreement is not biological accuracy."""
    qc=qc_folder(sample_root,stage);out=qc/'agreement';out.mkdir(exist_ok=True)
    summary=read_json(qc/'summary.json')
    try:
        with (qc/'sample.tsv').open(newline='') as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        per=int(cfg.get('evaluation_agreement_pairs_per_stratum',200));selected=[];counts=Counter()
        # Reservoir ordering is random within strata; reshuffle indices with a
        # fixed seed before selecting an additional bounded subset.
        indices=defaultdict(list)
        for i,row in enumerate(rows):indices[row['stratum']].append(i)
        rng=random.Random(int(cfg.get('evaluation_seed',20261115))+1)
        keep=set()
        for key,values in sorted(indices.items()):
            chosen=rng.sample(values,min(per,len(values)));keep.update(chosen);counts[key]=len(chosen)
        with gzip.open(out/'sample.R1.fastq.gz','wb') as left,gzip.open(out/'sample.R2.fastq.gz','wb') as right:
            for i,(a,b) in enumerate(fragments(qc/'sample.R1.fastq.gz',qc/'sample.R2.fastq.gz')):
                if i in keep:left.write(b''.join(a));right.write(b''.join(b));selected.append(rows[i])
        if not selected:
            result=dict(status='no_candidates',pairs=0);atomic_json(out/'summary.json',result);return result
        from partitioned_build import database_snapshot
        full=database_snapshot(cfg);qcfg=dict(cfg,database_layout='partitioned',union_memory_mapping=True,
                                             min_free_gib=min(float(cfg['min_free_gib']),1.0))
        for key in ('confidence','minimum_hit_groups','minimum_base_quality'):
            qcfg[key]=summary['request'][key]
        work=out/'query';work.mkdir(exist_ok=True)
        sample=dict(sample_id='evaluation/'+stage,r1=str(out/'sample.R1.fastq.gz'),r2=str(out/'sample.R2.fastq.gz'))
        tax=routing_taxonomy(cfg,full['taxonomy_dir']);query=query_once(qcfg,sample,full,work,taxonomy=tax)
        transition=defaultdict(lambda:[0,0.0]);detail=[]
        for row,tid in zip(selected,taxids(query['taxids'],len(selected)),strict=True):
            info=tax.describe(tid);group=biological_group(info)
            weight=int(row['stratum_pairs'])/counts[row['stratum']]
            key=(row['group'],group);transition[key][0]+=1;transition[key][1]+=weight
            detail.append(dict(ordinal=row['ordinal'],stratum=row['stratum'],weight=weight,
                               stage_group=row['group'],stage_taxid=row['taxid'],full_group=group,
                               full_taxid=tid,full_species_taxid=info['species_taxid']))
        write_table(out/'transitions.tsv',('stage_group','full_group','sample_pairs','weighted_pairs'),
                    [dict(stage_group=k[0],full_group=k[1],sample_pairs=v[0],weighted_pairs=v[1]) for k,v in sorted(transition.items())])
        write_table(out/'pairs.tsv',tuple(detail[0]),detail)
        result=dict(status='complete',pairs=len(selected),sampling='subsample_of_uniform_stratified_reservoir',
                    full_database=full['identity'],query=query,
                    matched_parameters={key:qcfg[key] for key in ('confidence','minimum_hit_groups','minimum_base_quality')},
                    interpretation='Larger-library agreement and routing sensitivity; neither call is ground truth.')
        atomic_json(out/'summary.json',result)
        for path in (out/'sample.R1.fastq.gz',out/'sample.R2.fastq.gz',Path(query['taxids'])):path.unlink(missing_ok=True)
        return result
    except Exception as exc:
        result=dict(status='failed',error=str(exc));atomic_json(out/'summary.json',result);return result


def classify_controls(sample_root,stage,cfg):
    """Measure optimistic in-reference recovery, separately from real accuracy."""
    qc=qc_folder(sample_root,stage);saved=read_json(qc/'summary.json')
    truth=qc/'controls.truth.tsv'
    if not truth.is_file():return dict(status='not_available',reason='No reference control fragments')
    with truth.open(newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows:return dict(status='no_candidates',reason='No usable reference control fragments')
    work=qc/'controls-query';work.mkdir(exist_ok=True)
    destinations=[];tasks=[];streams=defaultdict(list)
    for row,(a,b) in zip(rows,fragments(qc/'controls.R1.fastq.gz',qc/'controls.R2.fastq.gz'),strict=True):
        if stage=='group' or row['true_group'] in saved.get('selected_groups',[]):
            streams[row['true_group'] if stage=='species' else 'all'].append((row,a,b))
    qcfg=dict(cfg,union_memory_mapping=True,min_free_gib=min(float(cfg.get('min_free_gib',0)),1.0))
    for group,entries in streams.items():
        folder=work/group;folder.mkdir(exist_ok=True)
        for mate in (1,2):
            with gzip.open(folder/f'R{mate}.fastq.gz','wb') as handle:
                for item in entries:handle.write(b''.join(item[mate]))
        if stage=='group' and saved['request']['method']=='compact_competitive_kraken':
            from group_library import compact_snapshot
            snapshot=compact_snapshot(cfg)
        elif stage=='group' or group=='other_eukaryote':
            from partitioned_build import database_snapshot
            snapshot=database_snapshot(cfg)
        else:snapshot=category_snapshot(cfg,ROUTES[group])
        local=dict(qcfg,database_layout=snapshot['layout'])
        for key in ('confidence','minimum_hit_groups','minimum_base_quality'):
            local[key]=saved['request'][key]
        expected=(saved['request'].get('database') or saved['request'].get('original_database')) if stage=='group' \
                 else saved['request'].get('databases',{}).get(group)
        if expected and snapshot['identity']!=expected:
            raise RuntimeError('Reference controls would use a different classification database than the saved stage: '+group)
        tax=routing_taxonomy(cfg,snapshot['taxonomy_dir']);query_dir=folder/'query';query_dir.mkdir(exist_ok=True)
        query=query_once(local,dict(sample_id='reference-controls/'+group,r1=str(folder/'R1.fastq.gz'),
                         r2=str(folder/'R2.fastq.gz')),snapshot,query_dir,taxonomy=tax)
        tasks.append(query)
        for item,tid in zip(entries,taxids(query['taxids'],len(entries)),strict=True):
            row,a,b=item;info=tax.describe(tid)
            destinations.append(dict(ordinal=row['ordinal'],read_id=canonical_id(a[0]),
                r1_read_id=read_id(a[0]).decode(),r2_read_id=read_id(b[0]).decode(),
                group=biological_group(info),taxid=info['taxid'],species_taxid=info['species_taxid'],
                status=info['status'],bin='',source_group=group if stage=='species' else '',candidate_taxid=tid))
        for p in (folder/'R1.fastq.gz',folder/'R2.fastq.gz',Path(query['taxids'])):p.unlink(missing_ok=True)
    destinations.sort(key=lambda r:int(r['ordinal']))
    with gzip.open(qc/'controls.assignments.tsv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(destinations)
    result=dict(status='complete',queries=tasks,pairs=len(destinations),
                interpretation='Optimistic in-reference prefix-window recovery; species controls are conditional on known source group, not end-to-end accuracy.')
    atomic_json(qc/'controls-query.json',result)
    return result


def publish_evaluation_link(sample_root,stage,cfg):
    saved=read_json(qc_folder(sample_root,stage)/'summary.json')
    target=Path(cfg['work_dir'])/'qc/evaluate'/saved['sample_id']/stage
    target.parent.mkdir(parents=True,exist_ok=True)
    source=qc_folder(sample_root,stage).resolve()
    if target.is_symlink():
        if target.resolve()!=source:raise RuntimeError('Existing evaluation link targets another sample: '+str(target))
    elif target.exists():raise RuntimeError('Unowned evaluation destination: '+str(target))
    else:target.symlink_to(source,target_is_directory=True)
    print('EVALUATION_REPORT '+str(target/'evaluation.html'),flush=True)
    return target


def evaluate(sample_root,stage,cfg,*,truth_path=None):
    from stage_evaluate import (evaluate_stage, evaluate_controls, render_evaluation,
                                export_diagnostic_package)
    from stage_alignment import evaluate_alignment
    qc=qc_folder(sample_root,stage)
    saved=read_json(qc/'summary.json')
    cfg=dict(cfg,taxonomy_dir=saved.get('taxonomy_dir',cfg['taxonomy_dir']))
    # All independent methods are attempted even if one fails. The consolidated
    # report records the limitation and the command exits nonzero on actual failure.
    profile=cfg.get('evaluation_profile','full')
    if profile == 'full':
        agreement=classifier_agreement(sample_root,stage,cfg)
        try:alignment=evaluate_alignment(sample_root,stage,cfg)
        except Exception as exc:
            alignment=dict(status='failed',error=str(exc));atomic_json(qc/'alignment/summary.json',alignment)
    else:
        agreement=dict(status='skipped',reason='basic_profile_no_full_library_agreement')
        alignment=dict(status='skipped',reason='basic_profile_no_independent_alignment')
    atomic_json(qc/'evaluation-profile.json',dict(profile=profile,agreement=agreement.get('status'),alignment=alignment.get('status')))
    try:
        report=evaluate_stage(sample_root,stage,cfg,truth_path=truth_path)
        try:
            control_execution=classify_controls(sample_root,stage,cfg)
            report['controls_execution']=control_execution
            report['controls']=evaluate_controls(qc,stage,cfg)
        except Exception as exc:
            report['controls']=dict(status='failed',error=str(exc))
            atomic_json(qc/'controls-evaluation.json',report['controls'])
        report['evaluation_profile']=profile
        report=render_evaluation(qc,report)
    except Exception as exc:
        report=dict(status='failed',error=str(exc));atomic_json(qc/'evaluation.json',report)
    result=dict(stage=stage,profile=profile,agreement=agreement,alignment=alignment,evaluation=report)
    status_values=[x.get('status','') for x in (agreement,alignment,report)]
    result['status']='failed' if any(x in ('failed','error','invalid','incomplete','dependency_missing') for x in status_values) else \
                     ('partial' if any(x in ('partial','not_available','unavailable','skipped','not_tested','complete_with_limitations') for x in status_values) else 'complete')
    atomic_json(qc/'evaluate.json',result)
    summary=read_json(qc/'summary.json');summary['evaluation_status']=result['status']
    atomic_json(qc/'summary.json',summary)
    # The report is rendered before the orchestration status is known. Publish
    # the final bundle only after summary.json no longer says "pending".
    bundle=export_diagnostic_package(qc,report)
    atomic_json(qc/'diagnostic-bundle.json',bundle)
    publish_evaluation_link(sample_root,stage,cfg)
    return result


def report_samples(cfg,samples,out=None):
    out=Path(out) if out else Path(cfg.get('hierarchy_report_dir') or Path(cfg['work_dir'])/'qc'/'hierarchy')
    out.mkdir(parents=True,exist_ok=True);rows=[]
    for sample in samples:
        root=Path(cfg['output_dir'])/sample['sample_id']
        for stage in ('group','species'):
            qc=qc_folder(root,stage);p=qc/'summary.json'
            if p.is_file():
                s=read_json(p);rows.append(dict(sample=sample['sample_id'],stage=stage,
                    input_pairs=s['input_pairs'],original_pairs=s['expected_pairs'],
                    evaluation=s.get('evaluation_status','pending'),report=str(Path(cfg['work_dir'])/'qc/evaluate'/sample['sample_id']/stage/'evaluation.html')))
    write_table(out/'samples.tsv',('sample','stage','input_pairs','original_pairs','evaluation','report'),rows)
    body=''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in r)+'</tr>' for r in rows)
    (out/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>Metag hierarchy</title>'
        '<h1>Metag group and species stages</h1><p>All counts are paired fragments. Evidence agreement is not truth.</p>'
        '<table border="1"><tr><th>Sample</th><th>Stage</th><th>Stage input pairs</th><th>Original pairs</th><th>Evaluation</th><th>Report path</th></tr>'+body+'</table>')
    return rows


def parser():
    from workflow_selection import actions
    p=argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument('action',type=actions,help='all or ordered comma-separated stages')
    p.add_argument('--config');p.add_argument('--root');p.add_argument('--allow-local',action='store_true')
    p.add_argument('--sample',action='append');p.add_argument('--dat-list')
    p.add_argument('--groups',help='Selected biological groups; default human,plant; or all')
    p.add_argument('--include-uncertain',action='store_true',help='Also refine ambiguous/unclassified with the full union; potentially expensive')
    p.add_argument('--stage',choices=('group','species','both'),default='both')
    p.add_argument('--truth',type=Path,help='Optional real truth TSV; ordinal,true_taxid/true_group')
    p.add_argument('--threads',type=int,default=22)
    p.add_argument('--evaluation-pairs',type=int,help='Uniform reservoir limit per group/status stratum')
    p.add_argument('--evaluation',choices=('basic','full'),default='basic',help='basic QC/controls; full adds union agreement and competitive alignment')
    p.add_argument('--routing',choices=('compact','full'),default=None,help=argparse.SUPPRESS)
    p.add_argument('--replace',type=replace_bool,default=False)
    p.add_argument('--out',type=Path,help='Cohort report destination only')
    return p


def main(argv=None):
    p=parser();a=p.parse_args(argv);cfg=load_config(a.config,a.root);require_compute(a.allow_local)
    if a.threads<1: p.error('--threads must be positive')
    cfg.update(classify_threads=a.threads,classify_cores=a.threads+2,
               evaluation_profile=a.evaluation,group_routing=a.routing or cfg.get('group_routing','compact'))
    if a.evaluation_pairs is not None:
        if a.evaluation_pairs<1:p.error('--evaluation-pairs must be positive')
        cfg.update(evaluation_pairs_per_stratum=a.evaluation_pairs,
                   evaluation_target_pairs_per_stratum=a.evaluation_pairs)
    if a.include_uncertain and 'species' not in a.action:p.error('--include-uncertain is only for split-species')
    if a.out and 'report' not in a.action:p.error('--out is only for report; data paths come from config')
    samples=discover_samples(cfg,datasets=a.dat_list,only=a.sample)
    if not samples:raise RuntimeError('No complete selected input sample pairs')
    if a.truth and len(samples)!=1:p.error('--truth requires exactly one sample')
    groups=a.groups.split(',') if a.groups else None
    if a.action==['report']:report_samples(cfg,samples,a.out);return
    failures=[];blocked=set()
    # Stage-major execution keeps one scheduler job while making every split
    # available before any sample starts the optional expensive evaluation.
    for action in a.action:
        if action=='report':continue
        if len(blocked)==len(samples):break
        print(f'BATCH STAGE {action}: {len(samples)} samples; evaluation={a.evaluation}',flush=True)
        for number,sample in enumerate(samples,1):
            if sample['sample_id'] in blocked:continue
            root=Path(cfg['output_dir'])/sample['sample_id'];root.mkdir(parents=True,exist_ok=True)
            started=time.monotonic()
            try:
                with locked(root/'.hierarchy.lock'):
                    print(f'{utcnow()} START {sample["sample_id"]}: {action} ({number}/{len(samples)})',flush=True)
                    if action=='group':split_group(cfg,sample,replace=a.replace)
                    elif action=='species':split_species(cfg,sample,groups=groups,include_uncertain=a.include_uncertain,replace=a.replace)
                    elif action=='evaluate':
                        stages=[s for s in ('group','species') if a.stage in ('both',s) and (qc_folder(root,s)/'summary.json').is_file()]
                        if not stages:raise RuntimeError('No completed hierarchy stage; run split first')
                        failed=[]
                        for stage in stages:
                            result=evaluate(root,stage,cfg,truth_path=a.truth)
                            print(f'EVALUATE {sample["sample_id"]} {stage}: {result["status"]}',flush=True)
                            if result['status']=='failed':failed.append(stage)
                        if failed:raise RuntimeError('Evaluation failed: '+','.join(failed))
                    duration=time.monotonic()-started
                    print(f'{utcnow()} DONE {sample["sample_id"]}: {action}; seconds={duration:.1f}',flush=True)
                    progress=root/'workflow-progress.json';previous=read_json(progress) if progress.exists() else {}
                    previous[action]=dict(status='complete',seconds=duration,completed_at=utcnow(),evaluation_profile=a.evaluation)
                    atomic_json(progress,previous)
            except Exception as exc:
                blocked.add(sample['sample_id']);failures.append(sample['sample_id']+': '+str(exc))
                print('FAILED '+failures[-1],file=sys.stderr,flush=True)
        # Publish useful completion state even if the later evaluation is stopped.
        report_samples(cfg,samples,a.out)
    if not a.action or a.action==['report']:report_samples(cfg,samples,a.out)
    if failures:raise SystemExit('Some requested stages failed; independent results retained:\n'+'\n'.join(failures))


if __name__=='__main__':main()
