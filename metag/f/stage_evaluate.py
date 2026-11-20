#!/usr/bin/env python3
"""Automatic, bounded evaluation after group or species splitting.

Composition and agreement are diagnostics, not accuracy estimates. Accuracy
metrics require a supplied true-label table. Every sample row retains its
original ordinal; optional truth matching uses bounded SQLite storage, so
species assignments may be emitted in group order.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import gzip
import html
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import sqlite3
import statistics
import time
import zipfile

from common import atomic_json, load_config, locked
from fastq_stream import fragments, read_id
from sample_diagnose import Metrics

SCHEMA = 'metag-stage-evaluation-v1'
UNCLASSIFIED = '__unclassified__'
MISSING = '__not_processed__'
UNKNOWN_TRUTH = '__truth_not_available__'
_BASE = {65: 0, 67: 1, 71: 2, 84: 3}
_SRA_MATE = re.compile(r'^([SED]RR\d+\.\d+)\.[12]$')
_PALETTE = ('#147d92','#bd5b26','#7774ad','#54a36c','#c94567','#c0931a',
            '#2e6eb5','#7f7f3c','#9a60a0','#409e9d','#a76151','#67767b')


def stage_qc(sample_root, stage):
    if stage not in ('group','species'):
        raise ValueError('stage must be group or species')
    return Path(sample_root)/'group'/('qc' if stage=='group' else 'species.qc')


@contextmanager
def text_table(path):
    path=Path(path)
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8',newline='') as handle:
        yield csv.DictReader(handle,delimiter='\t')


def write_table(path, rows, fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    rows=list(rows)
    fields=fields or list(dict.fromkeys(k for row in rows for k in row)) or ['status']
    with path.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def json_file(path, default=None):
    try:return json.loads(Path(path).read_text(encoding='utf-8'))
    except FileNotFoundError:return default


def pair_id(value):
    fields=str(value).lstrip('@').split()
    if not fields:raise ValueError('Empty read identifier')
    value=fields[0]
    if value.endswith(('/1','/2')):return value[:-2]
    match=_SRA_MATE.fullmatch(value)
    return match.group(1) if match else value


def state(path):
    path=Path(path)
    try:
        stat=path.stat()
        return dict(path=str(path.resolve()),size=stat.st_size,mtime_ns=stat.st_mtime_ns)
    except OSError:return dict(path=str(path.absolute()),exists=False)


def accounting(sample_root, qc):
    summary=json_file(qc/'summary.json',{})
    checks=[];errors=[];notes=[];rows=[]
    path=qc/'counts.tsv'
    if not path.is_file():
        return dict(status='failed',errors=['counts.tsv missing'],checks=[]),summary,rows
    with text_table(path) as reader:
        for row in reader:rows.append(row)
    def check(name,observed,expected):
        ok=observed==expected
        checks.append(dict(check=name,observed=observed,expected=expected,passed=ok))
        if not ok:errors.append(f'{name}: {observed} != {expected}')
    try:
        counts=sum(int(row.get('pairs',row.get('fragments',''))) for row in rows)
        if any(int(row.get('pairs',row.get('fragments','')))<0 for row in rows):raise ValueError('negative pairs')
    except (TypeError,ValueError):
        counts=None;errors.append('counts.tsv requires nonnegative integer pairs')
    input_pairs=summary.get('input_pairs')
    written=summary.get('written_pairs')
    assignments=summary.get('assignment_rows')
    if input_pairs is None or written is None:
        errors.append('summary.json requires input_pairs and written_pairs')
    else:
        check('writer_pair_conservation',int(written),int(input_pairs))
        check('counts_vs_writer',counts,int(written))
    if assignments is not None and input_pairs is not None:
        check('assignment_rows_vs_input',int(assignments),int(input_pairs))
    else:notes.append('assignment_rows not recorded; assignment total not independently recounted')
    if qc.name=='qc' and summary.get('expected_pairs') is not None and input_pairs is not None:
        check('group_input_vs_original_population',int(input_pairs),int(summary['expected_pairs']))
    if qc.name=='species.qc' and input_pairs is not None:
        parent_counts=Path(sample_root)/'group/qc/counts.tsv'
        if parent_counts.is_file() and summary.get('processed_groups') is not None:
            # selected_groups is the biological target scope. processed_groups
            # also includes ambiguous/unclassified inputs when their recovery
            # is requested; conservation must use the latter population.
            processed=set(summary['processed_groups'])
            with text_table(parent_counts) as parent_rows:
                processed_total=sum(int(r.get('pairs',0)) for r in parent_rows if Path(r.get('bin','')).name in processed)
            check('species_input_vs_processed_parent_bins',int(input_pairs),processed_total)
        if summary.get('unprocessed_pairs') is not None and summary.get('expected_pairs') is not None:
            check('species_processed_plus_preserved_parent_pairs',int(input_pairs)+int(summary['unprocessed_pairs']),
                  int(summary['expected_pairs']))
    if not (qc/'assignments.tsv.gz').is_file():errors.append('assignments.tsv.gz missing')
    paths={};missing_size=0
    for row in rows:
        n=int(row.get('pairs',row.get('fragments','0')) or 0)
        if n==0:continue
        # Species refinement may resolve an uncertain source into a biological
        # group or demote a conflicting claim. The source describes the input
        # branch, not a required final classification.
        if qc.name!='species.qc' and row.get('source_group') and row.get('group') and row['source_group']!=row['group']:
            errors.append('source_group differs from final group in counts: '+row.get('bin',''))
        for mate in ('r1','r2'):
            if not row.get(mate):
                errors.append('counts.tsv missing '+mate+' for nonempty bin '+row.get('bin',''));continue
            file=Path(row[mate]);file=file if file.is_absolute() else Path(sample_root)/file
            expected=row.get('bytes_'+mate)
            if file in paths and paths[file]!=expected:
                errors.append('Conflicting expected file sizes: '+str(file))
            paths[file]=expected
    for file,expected in paths.items():
        if not file.is_file():errors.append('Expected bin file missing: '+str(file));continue
        if expected not in (None,''):
            check('artifact_bytes:'+str(file.relative_to(sample_root) if file.is_relative_to(sample_root) else file),
                  file.stat().st_size,int(expected))
        else:missing_size+=1
    if missing_size:notes.append(f'{missing_size} artifacts had existence checks only; expected sizes were not recorded')
    return dict(status='failed' if errors else 'complete',input_pairs=input_pairs,written_pairs=written,
                counts_pairs=counts,assignment_rows=assignments,checked_artifacts=len(paths),
                checks=checks,errors=errors,notes=notes,
                scope='full writer/count accounting plus artifact existence/size; output FASTQs are not reread in full'),summary,rows


def canonical_feature_map():
    canonical=[]
    for code in range(256):
        current=code;reverse=0
        for _ in range(4):reverse=(reverse<<2)|(3-(current&3));current>>=2
        canonical.append(min(code,reverse))
    columns={code:i for i,code in enumerate(sorted(set(canonical)))}
    return [columns[x] for x in canonical]


_CANONICAL=canonical_feature_map()


def tetranucleotide_features(sequences):
    """136 reverse-complement canonical 4-mers; never cross mate boundaries."""
    counts=[0]*136;total=0
    for sequence in sequences:
        code=valid=0
        for base in sequence.upper():
            if base not in _BASE:code=valid=0;continue
            code=((code<<2)|_BASE[base])&255;valid+=1
            if valid>=4:counts[_CANONICAL[code]]+=1;total+=1
    return ([math.sqrt(value/total) for value in counts] if total else None),total


def combine_metrics(strata, sizes):
    result=Metrics()
    attributes=[x for x in Metrics.__slots__ if x!='lengths']
    for name,metrics in strata.items():
        weight=sizes[name]/metrics.fragments if metrics.fragments else 0
        for key in attributes:setattr(result,key,getattr(result,key)+getattr(metrics,key)*weight)
        for length,count in metrics.lengths.items():result.lengths[length]+=count*weight
    return result.report()


def sampled_data(qc):
    metadata=qc/'sample.tsv'
    left,right=qc/'sample.R1.fastq.gz',qc/'sample.R2.fastq.gz'
    if not metadata.is_file() or not left.is_file() or not right.is_file():
        return dict(status='not_available',reason='sample.tsv and both sampled FASTQs are required'),[],[]
    with text_table(metadata) as reader:rows=list(reader)
    strata=defaultdict(Metrics);sizes={};declared={};observations=[];errors=[];seen=set();lengths=[[],[]]
    actual=0
    for actual,pair in enumerate(fragments(left,right),1):
        if actual>len(rows):raise RuntimeError('Sample FASTQ has more pairs than sample.tsv')
        row=rows[actual-1]
        ordinal=int(row['ordinal'])
        if ordinal<1 or ordinal in seen:raise RuntimeError('Duplicate/invalid sample original ordinal')
        seen.add(ordinal)
        ids=[pair_id(read_id(rec[0]).decode('utf-8','backslashreplace')) for rec in pair]
        if ids[0]!=ids[1] or ids[0]!=pair_id(row['read_id']):raise RuntimeError('sample.tsv / FASTQ read ID mismatch')
        name=row.get('stratum') or row.get('group') or 'all'
        population=int(row.get('stratum_pairs') or 0)
        selected=int(row.get('sampled_pairs') or 0)
        if population<1 or selected<1 or selected>population:
            raise RuntimeError('Each sample row needs positive stratum_pairs >= sampled_pairs')
        if name in sizes and (sizes[name]!=population or declared[name]!=selected):
            raise RuntimeError('Sampling stratum N/n values are inconsistent')
        sizes[name]=population;declared[name]=selected;strata[name].add_fragment(pair)
        seqs=[record[1].strip() for record in pair]
        for mate,seq in enumerate(seqs):lengths[mate].append(len(seq))
        features,windows=tetranucleotide_features(seqs)
        sequence=b''.join(seqs).upper()
        observations.append(dict(kind='read',ordinal=ordinal,read_id=row['read_id'],group=row.get('group',''),
            species_taxid=row.get('species_taxid','0'),taxid=row.get('taxid','0'),status=row.get('status',''),
            stratum=name,weight=population/selected,gc_fraction=(sequence.count(b'G')+sequence.count(b'C'))/len(sequence),
            valid_4mer_windows=windows,features=features,reference=''))
    if actual!=len(rows):raise RuntimeError('Sample FASTQ has fewer pairs than sample.tsv')
    records=[]
    for name,metrics in sorted(strata.items()):
        if metrics.fragments!=declared[name]:raise RuntimeError('Sampled stratum count differs from sampled_pairs')
        records.append(dict(stratum=name,stratum_pairs=sizes[name],sampled_pairs=metrics.fragments,
                            inverse_sampling_weight=sizes[name]/metrics.fragments,**metrics.report()))
    write_table(qc/'sample-metrics.tsv',records)
    report=dict(status='complete' if actual else 'empty',sampled_pairs=actual,
        represented_pairs=sum(sizes.values()),strata=len(strata),
        weighting='inverse sampling fraction N_h/n_h within each represented stratum',
        weighted_population_estimates=combine_metrics(strata,sizes),
        median_r1_length=statistics.median(lengths[0]) if lengths[0] else None,
        median_r2_length=statistics.median(lengths[1]) if lengths[1] else None,
        scope='QC metrics estimate represented strata only; absent strata and species are not inferred',errors=errors)
    return report,observations,records


def reference_catalog(cfg,qc=None):
    from group_library import group_directory
    library=group_directory(cfg) if cfg.get('library_dir') or cfg.get('group_library_dir') else Path('.')/'group'
    candidates=[library/'references.tsv',Path(cfg.get('work_dir','.'))/'qc/reference-audit/species-representatives.tsv']
    if qc is not None:
        # Existing flat results need no compact group index. The preceding
        # alignment step can have created the representative table by its
        # metadata-only audit; reuse it for anchors and actual control reads.
        alignment=Path(qc)/'alignment'
        try:metadata=json_file(alignment/'summary.json',{})
        except (OSError,ValueError):metadata={}
        if metadata.get('reference_table'):
            candidates.append(Path(metadata['reference_table']))
        candidates.append(alignment/'representatives.tsv')
    return next((path for path in candidates if path.is_file()),None),library


def bounded_reference_windows(path, *, bases_limit, read_lengths, insert_size, maximum, seed):
    """Deterministic sampling within a bounded reference prefix, not genome-wide."""
    left,right=read_lengths;span=max(insert_size,left+right)
    opener=gzip.open if str(path).endswith('.gz') else open
    contigs=[];name=None;parts=[];total=0
    with opener(path,'rb') as handle:
        while total<bases_limit:
            raw=handle.readline(1<<20)
            if not raw:break
            if raw.startswith(b'>'):
                if name is not None and parts:contigs.append((name,b''.join(parts)))
                name=raw[1:].split()[0].decode('utf-8','replace');parts=[]
            else:
                seq=b''.join(raw.split()).upper()[:bases_limit-total]
                if not seq:continue
                if name is None:raise RuntimeError('FASTA sequence precedes first header: '+str(path))
                parts.append(seq);total+=len(seq)
    if name is not None and parts:contigs.append((name,b''.join(parts)))
    regions=[];possible=0
    for contig,sequence in contigs:
        for match in re.finditer(rb'[ACGT]+',sequence):
            length=match.end()-match.start()
            if length>=span:
                count=length-span+1
                regions.append((possible,possible+count,contig,sequence,match.start()))
                possible+=count
    chosen=sorted(random.Random(seed).sample(range(possible),min(maximum,possible)))
    result=[];idx=0;complement=bytes.maketrans(b'ACGT',b'TGCA')
    for value in chosen:
        while value>=regions[idx][1]:idx+=1
        start,end,contig,sequence,offset=regions[idx]
        position=offset+value-start
        result.append(dict(contig=contig,position=position,
            sequence1=sequence[position:position+left].decode(),
            sequence2=sequence[position+span-right:position+span].translate(complement)[::-1].decode(),
            insert_size=span,sampling_scope='bounded_reference_prefix_windows'))
    return result,total


def _export_controls(qc,anchors):
    truth=[]
    with gzip.open(qc/'controls.R1.fastq.gz','wt') as left,gzip.open(qc/'controls.R2.fastq.gz','wt') as right:
        for ordinal,row in enumerate(anchors,1):
            name=f'control{ordinal}'
            for mate,handle in ((1,left),(2,right)):
                seq=row['sequence'+str(mate)]
                handle.write(f'@{name}/{mate}\n{seq}\n+\n'+('I'*len(seq))+'\n')
            truth.append(dict(ordinal=ordinal,read_id=name,true_group=row['group'],
                true_taxid=row.get('species_taxid') or row.get('taxid') or '0',
                true_species_taxid=row.get('species_taxid') or '0',
                reference=row.get('reference',''),control_kind='reference_prefix_positive_control',
                source_in_classification_reference='true'))
    write_table(qc/'controls.truth.tsv',truth,
                ['ordinal','read_id','true_group','true_taxid','true_species_taxid','reference','control_kind','source_in_classification_reference'])


def prepare_reference_anchors(qc,cfg):
    """Return qc/anchors.tsv and export coherent paired positive controls, or None."""
    qc=Path(qc);qc.mkdir(parents=True,exist_ok=True)
    catalog,library=reference_catalog(cfg,qc)
    prepared=next((p for p in (library/'qc/anchors.tsv',library/'qc/anchors.tsv.gz') if p.is_file()),None)
    limits=dict(references_per_group=int(cfg.get('evaluation_anchor_references',5)),
        bases_per_reference=int(cfg.get('evaluation_anchor_bases',1<<20)),
        pairs_per_group=int(cfg.get('evaluation_anchor_pairs',100)),
        read_lengths=[int(cfg.get('evaluation_read_length1',150)),int(cfg.get('evaluation_read_length2',150))],
        insert_size=int(cfg.get('evaluation_anchor_insert',350)),seed=int(cfg.get('evaluation_seed',20261115)))
    if min(limits['references_per_group'],limits['bases_per_reference'],limits['pairs_per_group'],*limits['read_lengths'])<1:
        raise ValueError('Reference anchor limits must be positive')
    anchors=[];issues=[];selected=[]
    if prepared:
        with text_table(prepared) as reader:
            group_counts=Counter()
            for row in reader:
                if not all(row.get(k) for k in ('group','sequence1','sequence2')):
                    issues.append('prepared anchor lacks group or paired sequences');continue
                if group_counts[row['group']]>=limits['pairs_per_group']:continue
                if not re.fullmatch('[ACGTacgt]+',row['sequence1']+row['sequence2']):continue
                group_counts[row['group']]+=1;anchors.append(row)
        report=dict(status='complete' if anchors else 'not_available',source=str(prepared),
                    sampling_scope='supplied_reference_anchors',issues=issues,limits=limits)
    elif not catalog:
        report=dict(status='not_available',reason='No group references.tsv or species-representatives.tsv',limits=limits)
    else:
        # Keep at most five distinct species per group; deterministic selection
        # is deliberately not claimed to represent the entire reference library.
        best=defaultdict(dict)
        with text_table(catalog) as reader:
            for row in reader:
                group=row.get('group') or 'unresolved'
                tid=int(row.get('species_taxid') or row.get('taxid') or row.get('source_taxid') or 0)
                raw=Path(row.get('path') or row.get('resolved_path') or '')
                if not raw.is_absolute():raw=Path(cfg.get('reference_raw','.'))/raw
                value=dict(row,path=str(raw),group=group,species_taxid=tid,
                           taxid=int(row.get('taxid') or row.get('source_taxid') or tid))
                previous=best[group].get(tid)
                if previous is None or str(raw)<previous['path']:best[group][tid]=value
                if len(best[group])>limits['references_per_group']:
                    del best[group][max(best[group])]
        selected=[row for group in sorted(best) for _,row in sorted(best[group].items())]
        request=dict(schema='reference-prefix-anchors-v1',catalog=state(catalog),limits=limits,
                     selected=[{**{k:r[k] for k in ('group','taxid','species_taxid')},**state(r['path'])} for r in selected])
        cache=Path(cfg.get('work_dir',qc))/'cache/evaluation';cache.mkdir(parents=True,exist_ok=True)
        with locked(cache/'anchors.lock'):
            try:saved=json_file(cache/'anchors.json',{})
            except (OSError,ValueError):saved={}
            if saved.get('request')==request and (cache/'anchors.tsv').is_file():
                with text_table(cache/'anchors.tsv') as reader:anchors=list(reader)
                report=saved['report'];report={**report,'cache_reused':True}
            else:
                totals=Counter(r['group'] for r in selected);scanned=0
                for index,row in enumerate(selected):
                    maximum=max(1,limits['pairs_per_group']//totals[row['group']])
                    try:
                        windows,bases=bounded_reference_windows(row['path'],bases_limit=limits['bases_per_reference'],
                            read_lengths=limits['read_lengths'],insert_size=limits['insert_size'],maximum=maximum,
                            seed=limits['seed']+index)
                        scanned+=bases
                    except (OSError,ValueError,RuntimeError,EOFError) as exc:
                        issues.append(dict(reference=row['path'],reason=str(exc)));continue
                    for window in windows:
                        anchors.append(dict(anchor_id=f'anchor{len(anchors)+1}',group=row['group'],taxid=row['taxid'],
                            species_taxid=row['species_taxid'],reference=row['path'],**window))
                report=dict(status='complete' if anchors and not issues else 'partial' if anchors else 'not_available',
                    source=str(catalog),sampling_scope='bounded_reference_prefix_windows',
                    selection='smallest species taxids, one deterministic assembly per species, up to the per-group cap',
                    reference_files_considered=len(selected),reference_prefix_bases_scanned=scanned,
                    limits=limits,issues=issues,cache_reused=False,
                    interpretation='Exploratory reference anchors and optimistic positive controls; not genome-wide random sampling or held-out accuracy')
                write_table(cache/'anchors.tsv',anchors)
                atomic_json(cache/'anchors.json',dict(request=request,report=report))
    report.update(anchor_pairs=len(anchors),groups=dict(Counter(x['group'] for x in anchors)))
    atomic_json(qc/'anchors.json',report)
    if not anchors:
        for name in ('anchors.tsv','controls.R1.fastq.gz','controls.R2.fastq.gz','controls.truth.tsv'):
            (qc/name).unlink(missing_ok=True)
        return None
    write_table(qc/'anchors.tsv',anchors);_export_controls(qc,anchors)
    return qc/'anchors.tsv'


def _anchor_observations(path):
    result=[]
    if path is None:return result
    with text_table(path) as reader:
        for row in reader:
            seqs=[row['sequence1'].encode(),row['sequence2'].encode()]
            features,windows=tetranucleotide_features(seqs)
            seq=b''.join(seqs).upper()
            result.append(dict(kind='reference',ordinal='',read_id=row.get('anchor_id',''),group=row['group'],
                species_taxid=row.get('species_taxid','0'),taxid=row.get('taxid','0'),status='reference_anchor',
                stratum='reference:'+row['group'],weight=1,gc_fraction=(seq.count(b'G')+seq.count(b'C'))/len(seq),
                valid_4mer_windows=windows,features=features,reference=row.get('reference','')))
    return result


def composition(qc, observations, anchors):
    rows=[x for x in observations+anchors if x['features'] is not None]
    excluded=len(observations)+len(anchors)-len(rows)
    try:import numpy as np
    except ImportError:
        return dict(status='dependency_missing',reason='NumPy is required for PCA; install the declared metag environment')
    if len(rows)<2:return dict(status='insufficient_data',usable_points=len(rows),excluded_no_valid_4mers=excluded)
    reference_indices=[i for i,x in enumerate(rows) if x['kind']=='reference']
    reference_groups={rows[i]['group'] for i in reference_indices}
    reference_fit=len(reference_indices)>=3 and len(reference_groups)>=2
    fit_indices=reference_indices if reference_fit else [i for i,x in enumerate(rows) if x['kind']=='read']
    if len(fit_indices)<2:return dict(status='insufficient_data',reason='Neither sufficient reference anchors nor at least two reads')
    matrix=np.asarray([x['features'] for x in rows],dtype=np.float64)
    fitted=matrix[fit_indices];center=fitted.mean(axis=0);fitted=fitted-center
    _,singular,components=np.linalg.svd(fitted,full_matrices=False)
    total=float((singular**2).sum())
    if total<=1e-20:return dict(status='zero_variance',usable_points=len(rows),fit_points=len(fit_indices))
    for i in range(min(2,len(components))):
        if components[i,np.abs(components[i]).argmax()]<0:components[i]*=-1
    scores=(matrix-center)@components[:2].T
    variances=[float(singular[i]**2/total) if i<len(singular) else 0 for i in (0,1)]
    coordinates=[]
    for index,row in enumerate(rows):
        coordinates.append({k:v for k,v in row.items() if k!='features'}|dict(pc1=float(scores[index,0]),
                           pc2=float(scores[index,1]) if scores.shape[1]>1 else 0))
    write_table(qc/'composition.tsv',coordinates)
    plot_composition(qc/'composition.svg',coordinates,variances,reference_fit)
    return dict(status='complete',features=136,transformation='sqrt(canonical 4-mer frequency); centered PCA',
        fit_scope='reference_anchors_only' if reference_fit else 'sample_only',
        reference_anchor_groups=len(reference_groups),reference_anchor_points=len(reference_indices),
        usable_read_points=sum(x['kind']=='read' for x in rows),excluded_no_valid_4mers=excluded,
        explained_variance_fraction=variances,
        interpretation='Short-read nucleotide composition is exploratory. Separation, overlap and distance are not species accuracy, ancestry inference, or proof of contamination.')


def plot_composition(path,rows,variance,reference_fit):
    width,height=1100,700;left,top,plotw,ploth=80,70,780,530
    xs=[x['pc1'] for x in rows];ys=[x['pc2'] for x in rows]
    xmin,xmax=min(xs),max(xs);ymin,ymax=min(ys),max(ys)
    xpad=max((xmax-xmin)*.06,1e-6);ypad=max((ymax-ymin)*.06,1e-6)
    xmin-=xpad;xmax+=xpad;ymin-=ypad;ymax+=ypad
    labels=sorted({x['group'] for x in rows});colors={g:_PALETTE[i%len(_PALETTE)] for i,g in enumerate(labels)}
    def px(v):return left+(v-xmin)/(xmax-xmin)*plotw
    def py(v):return top+ploth-(v-ymin)/(ymax-ymin)*ploth
    elements=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#233746} .small{font-size:12px}</style>',
        '<text x="80" y="30" font-size="21">Paired-read composition: exploratory PCA</text>',
        '<text x="80" y="51" class="small">'+('Fit on bounded reference anchors; reads projected' if reference_fit else 'Sample-only fit: reference anchors unavailable or insufficient')+'</text>',
        f'<path d="M {left} {top} V {top+ploth} H {left+plotw}" fill="none" stroke="#405467"/>']
    for i in range(6):
        x=xmin+(xmax-xmin)*i/5;y=ymin+(ymax-ymin)*i/5
        elements.append(f'<text x="{px(x):.2f}" y="{top+ploth+22}" text-anchor="middle" class="small">{x:.3g}</text>')
        elements.append(f'<text x="{left-10}" y="{py(y)+4:.2f}" text-anchor="end" class="small">{y:.3g}</text>')
    for row in sorted(rows,key=lambda x:x['kind']=='reference'):
        x,y=px(row['pc1']),py(row['pc2']);color=colors[row['group']]
        title=html.escape(f"{row['kind']} | {row['group']} | species {row['species_taxid']} | {row['read_id']}")
        if row['kind']=='reference':
            elements.append(f'<rect x="{x-3:.2f}" y="{y-3:.2f}" width="6" height="6" fill="white" stroke="{color}" stroke-width="1.2"><title>{title}</title></rect>')
        else:elements.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.4" fill="{color}" opacity=".48"><title>{title}</title></circle>')
    for i,label in enumerate(labels):
        y=90+i*24
        elements.append(f'<circle cx="885" cy="{y}" r="5" fill="{colors[label]}"/><text x="900" y="{y+4}" class="small">{html.escape(label)}</text>')
    elements.extend([f'<text x="{left+plotw/2}" y="650" text-anchor="middle">PC1 ({variance[0]:.1%} of fit variance)</text>',
        f'<text x="20" y="{top+ploth/2}" transform="rotate(-90 20 {top+ploth/2})" text-anchor="middle">PC2 ({variance[1]:.1%})</text>',
        '<text x="880" y="520" class="small">● sampled read pair</text><text x="880" y="542" class="small">□ reference pair</text>',
        '<text x="80" y="680" class="small">136 canonical 4-mer features. Cluster geometry is not a measure of classification accuracy.</text></svg>'])
    Path(path).write_text(''.join(elements),encoding='utf-8')


def _predicted_label(row,stage):
    if stage=='group':
        group=row.get('group','')
        return group if group and group not in ('unresolved','unclassified','unknown','ambiguous') else UNCLASSIFIED
    tid=int(row.get('species_taxid') or 0)
    return str(tid) if tid>0 else UNCLASSIFIED


def evaluate_truth(assignments,truth_path,out,stage,cfg=None,*,prefix='truth',selected_groups=()):
    """Evaluate supplied labels with an on-disk ordinal join and honest denominators."""
    cfg=cfg or {};out=Path(out);out.mkdir(parents=True,exist_ok=True)
    database=out/f'.{prefix}.sqlite';database.unlink(missing_ok=True)
    connection=sqlite3.connect(database)
    connection.execute('PRAGMA journal_mode=OFF');connection.execute('PRAGMA synchronous=OFF')
    connection.execute('PRAGMA cache_size=-16384');connection.execute('PRAGMA temp_store=FILE')
    connection.execute('CREATE TABLE truth (ordinal INTEGER PRIMARY KEY, read_id TEXT, label TEXT, true_group TEXT, seen INTEGER DEFAULT 0)')
    truth_count=excluded=assignment_count=matched=missing=unlabeled=0
    confusion=Counter();actual=Counter();predicted=Counter();correct=Counter();target=Counter()
    taxonomy=None;target_groups=set(selected_groups);normalization=Counter()
    def tax_info(tid):
        nonlocal taxonomy
        if not cfg.get('taxonomy_dir'):return None
        if taxonomy is None:
            from taxonomy import Taxonomy
            taxonomy=Taxonomy(cfg['taxonomy_dir'])
        return taxonomy.describe(tid)
    def observe(label,pred,truth_group,pred_group):
        nonlocal excluded
        confusion[label,pred]+=1
        if label==UNKNOWN_TRUTH:excluded+=1;return
        actual[label]+=1
        if pred not in (UNCLASSIFIED,MISSING):predicted[pred]+=1
        if pred==label:correct[label]+=1
        if target_groups:
            true_target=truth_group in target_groups;pred_target=pred_group in target_groups and pred not in (UNCLASSIFIED,MISSING)
            target['true_target']+=int(true_target);target['predicted_target']+=int(pred_target)
            target['correct_target']+=int(true_target and pred_target and pred==label)
    try:
        with text_table(truth_path) as reader:
            if not reader.fieldnames or 'ordinal' not in reader.fieldnames:raise ValueError('Truth TSV needs ordinal')
            if stage=='group' and not ({'true_group','true_taxid'}&set(reader.fieldnames)):
                raise ValueError('Group truth needs true_group or true_taxid with taxonomy')
            if stage=='species' and 'true_taxid' not in reader.fieldnames:raise ValueError('Species truth needs true_taxid')
            for row in reader:
                ordinal=int(row['ordinal'])
                if ordinal<1:raise ValueError('Truth ordinals must be positive')
                group=row.get('true_group','');tid=int(row.get('true_taxid') or 0)
                info=None
                if not group and tid:
                    info=tax_info(tid);group=info['group'] if info else ''
                if stage=='group':
                    label=group if group not in ('unresolved','unclassified','unknown','ambiguous') else ''
                elif row.get('true_species_taxid') not in (None,''):
                    species=int(row['true_species_taxid']);label=str(species) if species>0 else ''
                    normalization['explicit_true_species_taxid']+=1
                elif tid and cfg.get('taxonomy_dir'):
                    info=info or tax_info(tid);species=info['species_taxid']
                    label=str(species) if species else '';normalization['frozen_taxonomy_species_ancestor']+=1
                else:
                    label=str(tid) if tid>0 else '';normalization['supplied_taxid_asserted_as_species']+=1
                if stage=='group' and not label and tid and not cfg.get('taxonomy_dir'):
                    raise ValueError('true_taxid-only group truth requires taxonomy_dir')
                label=label or UNKNOWN_TRUTH
                try:connection.execute('INSERT INTO truth(ordinal,read_id,label,true_group) VALUES (?,?,?,?)',
                                       (ordinal,row.get('read_id',''),label,group))
                except sqlite3.IntegrityError as exc:raise ValueError('Duplicate truth ordinal '+str(ordinal)) from exc
                truth_count+=1
        connection.commit()
        with text_table(assignments) as reader:
            for row in reader:
                assignment_count+=1;ordinal=int(row['ordinal'])
                if ordinal<1:raise ValueError('Assignment ordinals must be positive')
                found=connection.execute('SELECT read_id,label,true_group,seen FROM truth WHERE ordinal=?',(ordinal,)).fetchone()
                if found is None:unlabeled+=1;continue
                identity,label,group,seen=found
                if seen:raise ValueError('Duplicate assignment for a truth-labeled original ordinal '+str(ordinal))
                if identity and pair_id(identity)!=pair_id(row.get('read_id','')):
                    raise ValueError('Truth/assignment read ID mismatch at ordinal '+str(ordinal))
                connection.execute('UPDATE truth SET seen=1 WHERE ordinal=?',(ordinal,));matched+=1
                observe(label,_predicted_label(row,stage),group,row.get('group',''))
        for label,group in connection.execute('SELECT label,true_group FROM truth WHERE seen=0'):
            missing+=1;observe(label,MISSING,group,'')
        metrics=[]
        for label in sorted(set(actual)|set(predicted)):
            tp=correct[label];fp=predicted[label]-tp;fn=actual[label]-tp
            precision=tp/(tp+fp) if tp+fp else None;recall=tp/(tp+fn) if tp+fn else None
            f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None
            metrics.append(dict(label=label,true_pairs=actual[label],predicted_pairs=predicted[label],tp=tp,fp=fp,fn=fn,
                                precision=precision,recall=recall,f1=f1))
        confusion_rows=[dict(true_label=a,predicted_label=b,pairs=n) for (a,b),n in sorted(confusion.items())]
        write_table(out/f'{prefix}-confusion.tsv',confusion_rows,['true_label','predicted_label','pairs'])
        write_table(out/f'{prefix}-metrics.tsv',metrics,['label','true_pairs','predicted_pairs','tp','fp','fn','precision','recall','f1'])
        tp=sum(correct.values());known=sum(actual.values());positive=sum(predicted.values())
        report=dict(status='complete',stage=stage,truth_rows=truth_count,known_label_rows=known,
            unknown_truth_rows=excluded,assignment_rows=assignment_count,matched_truth_rows=matched,
            missing_predictions=missing,assignment_rows_without_truth=unlabeled,correct_pairs=tp,
            assigned_pairs_on_labeled_truth=positive,
            precision_assigned=tp/positive if positive else None,recall=tp/known if known else None,
            f1=2*tp/(positive+known) if positive+known else None,
            classification_rate_on_labeled_truth=positive/known if known else None,
            truth_source=state(truth_path),assignment_source=state(assignments),
            taxid_normalization=dict(normalization),
            label_policy='group labels or species-level truth (explicit true_species_taxid, otherwise frozen-taxonomy ancestor when available); unclassified/missing calls contribute false negatives. Without taxonomy the supplied true_taxid must already be a species ID.',
            interpretation='Metrics apply only to the supplied true-label population. Unlabeled real reads are not declared correct or incorrect. Precision has labeled assigned predictions as denominator; recall includes labeled unclassified and missing calls.')
        if target_groups:
            hit=target['correct_target'];calls=target['predicted_target'];total=target['true_target']
            report['selected_group_scope']=dict(groups=sorted(target_groups),true_pairs=total,predicted_pairs=calls,
                correct_pairs=hit,precision=hit/calls if calls else None,recall=hit/total if total else None,
                interpretation='Target-group recall uses all supplied target-group truths, including those missed during group routing; non-target reads called into target groups contribute false positives.')
        atomic_json(out/f'{prefix}-evaluation.json',report)
        return report
    finally:
        connection.close();database.unlink(missing_ok=True)


def evaluate_controls(qc,stage,cfg):
    qc=Path(qc);truth=qc/'controls.truth.tsv';assignments=qc/'controls.assignments.tsv.gz'
    def pending(status,reason):
        report=dict(status=status,reason=reason)
        atomic_json(qc/'controls-evaluation.json',report)
        return report
    if not truth.is_file():return pending('not_available','No reference anchors generated')
    if not assignments.is_file():return pending('generated_not_run','Classifier has not yet been run on the exported positive controls')
    if assignments.stat().st_mtime_ns<truth.stat().st_mtime_ns:
        return pending('generated_not_run','Previous control assignments predate the currently generated control sequences')
    selected_groups=(json_file(qc/'summary.json',{}) or {}).get('selected_groups',[])
    if isinstance(selected_groups,str):selected_groups=[x for x in selected_groups.split(',') if x]
    control_scope='all_exported_reference_groups'
    if stage=='species' and selected_groups:
        with text_table(truth) as reader:rows=[r for r in reader if r.get('true_group') in selected_groups]
        target=qc/'controls.target.truth.tsv'
        write_table(target,rows,['ordinal','read_id','true_group','true_taxid','true_species_taxid','reference','control_kind','source_in_classification_reference'])
        truth=target;control_scope='conditional_on_known_source_group'
    report=evaluate_truth(assignments,truth,qc,stage,cfg,prefix='controls')
    report.update(control_kind='reference_prefix_positive_control',source_in_classification_reference=True,
                  control_scope=control_scope,selected_groups=selected_groups,
                  interpretation='Optimistic recovery control using windows from classification reference genomes. This checks routing/implementation and in-reference recovery; it is not held-out performance on real metagenomic samples.')
    atomic_json(qc/'controls-evaluation.json',report)
    return report


def export_diagnostic_package(qc,report):
    qc=Path(qc)
    description=('Metagenomics stage diagnostic bundle\n\n'
        'Contains the bounded stratified paired-read sample, original ordinals, sampling fractions, full counts, '
        'and evaluation evidence. No full raw sample or full split FASTQ bin is included.\n'
        'Composition and sampled alignment/classifier agreement are not accuracy. True-label metrics, if '
        'present, apply only to their declared labels. Reference controls are optimistic in-reference checks.\n')
    (qc/'diagnostic-README.txt').write_text(description)
    names=['diagnostic-README.txt','sample.R1.fastq.gz','sample.R2.fastq.gz','sample.tsv','counts.tsv','summary.json',
           'stage-evaluation.json','evaluation-profile.json','evaluation.html','sample-metrics.tsv','composition.tsv','composition.svg',
           'anchors.json','anchors.tsv','controls.truth.tsv','controls.target.truth.tsv','controls-evaluation.json','controls-confusion.tsv',
           'controls-metrics.tsv','truth-evaluation.json','truth-confusion.tsv','truth-metrics.tsv','agreement.json']
    if report.get('evaluation_profile','full')=='basic':names.remove('agreement.json')
    paths=[qc/name for name in names if (qc/name).is_file()]
    if report.get('evaluation_profile','full')=='full':
        for name in ('summary.json','sampled-evidence.tsv','reference-panel.tsv',
                     'omitted-references.tsv','panels.json'):
            path=qc/'alignment'/name
            if path.is_file():paths.append(path)
        for name in ('summary.json','transitions.tsv','pairs.tsv'):
            path=qc/'agreement'/name
            if path.is_file():paths.append(path)
    destination=qc/'diagnostic.zip'
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4) as archive:
        for path in paths:
            archive.write(path,str(path.relative_to(qc)),compress_type=zipfile.ZIP_STORED if path.suffix=='.gz' else zipfile.ZIP_DEFLATED)
    return dict(path='diagnostic.zip',bytes=destination.stat().st_size,files=len(paths),sampled_sequences_only=True)


def render_evaluation(qc,report):
    qc=Path(qc)
    def evidence(path):
        try:return json_file(path)
        except (OSError,ValueError) as exc:
            return dict(status='failed',reason='Cannot read '+str(path.relative_to(qc))+': '+str(exc))
    profile=json_file(qc/'evaluation-profile.json',{}).get('profile', 'full')
    report['evaluation_profile']=profile
    if profile=='basic':
        alignment=dict(status='skipped',reason='Not requested in basic evaluation profile')
        agreement=dict(status='skipped',reason='Not requested in basic evaluation profile')
    else:
        alignment=evidence(qc/'alignment/summary.json')
        agreement=evidence(qc/'agreement/summary.json') or evidence(qc/'agreement.json')
    controls=evidence(qc/'controls-evaluation.json')
    if alignment is not None:report['alignment']=alignment
    if agreement is not None:report['agreement']=agreement
    if controls is not None:report['controls']=controls
    components=('accounting','sample_qc','anchors','composition','alignment','agreement','truth','controls')
    failures=[k for k in components if (report.get(k) or {}).get('status') in ('failed','dependency_missing')]
    limited={'not_available','unavailable','skipped','insufficient_data','zero_variance','partial',
             'partially_tested','not_tested','generated_not_run','empty'}
    if failures:
        report['status']='incomplete';report['failed_components']=failures
    elif any((report.get(k) or {}).get('status') in limited for k in components):
        report['status']='complete_with_limitations';report.pop('failed_components',None)
    else:
        report['status']='complete';report.pop('failed_components',None)
    # Keep exact archive bytes outside the archive itself; otherwise this JSON
    # would necessarily contain a stale/self-referential ZIP size.
    report['diagnostic_bundle']=dict(path='diagnostic.zip',sampled_sequences_only=True,
                                     metadata_path='diagnostic-bundle.json')
    atomic_json(qc/'stage-evaluation.json',report)
    def text(value):return html.escape(str(value))
    statuses=[]
    for key in components:
        value=report.get(key) or {}
        statuses.append('<tr><td>'+text(key)+'</td><td>'+text(value.get('status','not_available'))+'</td><td>'+
                        text(value.get('reason') or value.get('interpretation') or value.get('scope',''))+'</td></tr>')
    image='<img src="composition.svg" alt="Exploratory paired-read composition PCA" style="width:100%;height:auto">' if (qc/'composition.svg').is_file() else ''
    details=json.dumps({k:report[k] for k in ('accounting','sample_qc','composition','truth','controls') if k in report},indent=2)
    content='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Metag stage evaluation</title>
<style>body{font:16px/1.6 system-ui,sans-serif;max-width:1200px;margin:32px auto;padding:0 22px;color:#203644;background:#f7fafb}h1,h2{line-height:1.3}table{border-collapse:collapse;width:100%;background:white}td,th{padding:10px;border-bottom:1px solid #d9e2e8;text-align:left;vertical-align:top}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#fff;padding:18px;font-size:12px}.note{background:#fff1d1;padding:18px;border-left:4px solid #d39e36}a{color:#09637c}</style>'''
    content+='<h1>'+text(report['stage'])+' split: automatic evaluation</h1>'
    content+='<p>Sample: '+text(report.get('sample_id',''))+' · Execution status: <strong>'+text(report['status'])+'</strong></p>'
    content+='<div class="note">A clean split and a clear PCA plot do not establish biological accuracy. Composition describes sequence similarity; alignment and classifier agreement describe concordance under their reference panels. Precision/recall are shown only when true labels are supplied. Reference-prefix controls are optimistic in-reference checks.</div>'
    content+='<h2>Checks and evidence</h2><table><tr><th>Component</th><th>Status</th><th>Meaning / limitation</th></tr>'+''.join(statuses)+'</table>'
    content+='<h2>Composition with reference anchors</h2>'+image
    content+='<p>Reads are projected into the reference-anchor PCA when enough reference groups are available. Canonical tetranucleotide features count within each mate, never across their artificial join. Reference windows are bounded prefix samples, not a representative genome-wide sample. Short reads from different taxa can overlap.</p>'
    content+='<h2>Files</h2><p><a href="diagnostic.zip">Bounded diagnostic ZIP</a> · <a href="stage-evaluation.json">Exact evaluation JSON</a> · <a href="sample-metrics.tsv">Sampling-weighted QC strata</a>'
    if (qc/'composition.tsv').is_file():content+=' · <a href="composition.tsv">PCA coordinates</a>'
    if profile=='full' and alignment is not None:content+=' · <a href="alignment/summary.json">Competitive alignment evidence</a>'
    if profile=='full' and agreement is not None:
        target='agreement/summary.json' if (qc/'agreement/summary.json').is_file() else 'agreement.json'
        content+=' · <a href="'+target+'">Sampled classifier agreement</a>'
    if controls is not None:content+=' · <a href="controls-evaluation.json">Reference control recovery</a>'
    content+='</p><details><summary>Numerical details and denominators</summary><pre>'+text(details)+'</pre></details></html>'
    (qc/'evaluation.html').write_text(content,encoding='utf-8')
    bundle=export_diagnostic_package(qc,report)
    atomic_json(qc/'diagnostic-bundle.json',bundle)
    return report


def evaluate_stage(sample_root,stage,cfg,*,truth_path=None):
    """Evaluate committed stage output; always publish explicit component states."""
    started=time.monotonic();qc=stage_qc(sample_root,stage);qc.mkdir(parents=True,exist_ok=True);cfg=dict(cfg)
    # Only replace evaluator-owned derived files. Never remove stage reads,
    # assignments, counts, external alignment evidence, or classifier controls.
    for name in ('composition.tsv','composition.svg','sample-metrics.tsv','truth-evaluation.json',
                 'truth-confusion.tsv','truth-metrics.tsv','controls-evaluation.json',
                 'controls-confusion.tsv','controls-metrics.tsv','controls.target.truth.tsv'):
        (qc/name).unlink(missing_ok=True)
    report=dict(schema=SCHEMA,stage=stage,sample_id=Path(sample_root).name,
                created_at=datetime.now(timezone.utc).isoformat(),status='complete',errors=[])
    try:
        audit,summary,counts=accounting(Path(sample_root),qc)
        report['accounting']=audit;report['sample_id']=summary.get('sample_id',report['sample_id'])
        if summary.get('taxonomy_dir'):cfg['taxonomy_dir']=summary['taxonomy_dir']
        report['stage_scope']={k:summary[k] for k in ('input_pairs','expected_pairs','selected_groups','processed_groups',
            'unprocessed_pairs','input_scope','scope') if k in summary}
    except (OSError,ValueError,TypeError,KeyError) as exc:
        summary={};report['accounting']=dict(status='failed',reason=str(exc));report['errors'].append(str(exc))
    observations=[]
    try:
        sample,observations,_=sampled_data(qc);report['sample_qc']=sample
        n=summary.get('input_pairs')
        if n:
            sample['represented_population_fraction']=sample.get('represented_pairs',0)/int(n)
            if sample.get('represented_pairs',0)>int(n):raise ValueError('Sampling strata claim more pairs than stage input')
    except (OSError,ValueError,RuntimeError,KeyError,EOFError) as exc:
        report['sample_qc']=dict(status='failed',reason=str(exc));report['errors'].append(str(exc))
    anchor_path=None
    try:
        settings=dict(cfg)
        for mate in (1,2):
            length=(report.get('sample_qc') or {}).get(f'median_r{mate}_length')
            if length:settings[f'evaluation_read_length{mate}']=int(length)
        anchor_path=prepare_reference_anchors(qc,settings)
        report['anchors']=json_file(qc/'anchors.json',dict(status='not_available'))
    except (OSError,ValueError,RuntimeError,KeyError,EOFError) as exc:
        report['anchors']=dict(status='failed',reason=str(exc));report['errors'].append(str(exc))
    try:report['composition']=composition(qc,observations,_anchor_observations(anchor_path))
    except (OSError,ValueError,RuntimeError,KeyError) as exc:
        report['composition']=dict(status='failed',reason=str(exc));report['errors'].append(str(exc))
    if truth_path:
        try:report['truth']=evaluate_truth(qc/'assignments.tsv.gz',truth_path,qc,stage,cfg,
                                          selected_groups=summary.get('selected_groups',()))
        except (OSError,ValueError,RuntimeError,KeyError,sqlite3.Error) as exc:
            report['truth']=dict(status='failed',reason=str(exc));report['errors'].append(str(exc))
    else:report['truth']=dict(status='not_supplied',interpretation='No real-sample accuracy estimate without true labels')
    try:report['controls']=evaluate_controls(qc,stage,cfg)
    except (OSError,ValueError,RuntimeError,KeyError,sqlite3.Error) as exc:
        report['controls']=dict(status='failed',reason=str(exc));report['errors'].append(str(exc))
    failures=[k for k,v in report.items() if isinstance(v,dict) and v.get('status') in ('failed','dependency_missing')]
    if failures:report['status']='incomplete';report['failed_components']=failures
    elif any((report.get(k) or {}).get('status') in ('not_available','insufficient_data','zero_variance','partial')
             for k in ('sample_qc','anchors','composition')):report['status']='complete_with_limitations'
    report['evaluation_seconds']=round(time.monotonic()-started,3)
    return render_evaluation(qc,report)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-root',required=True);parser.add_argument('--stage',choices=('group','species'),required=True)
    parser.add_argument('--config');parser.add_argument('--root');parser.add_argument('--truth')
    args=parser.parse_args(argv)
    result=evaluate_stage(args.sample_root,args.stage,load_config(args.config,args.root),truth_path=args.truth)
    print(json.dumps(result,indent=2))
    if result['status']=='incomplete':raise SystemExit(1)


if __name__=='__main__':main()
