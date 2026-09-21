#!/usr/bin/env python3
"""Non-destructive Cell 2020 postprocessing for GU final.

The published text does not specify the DAF window construction. Accept a
provenanced exclusion BED; never substitute segment LOD or carrier frequency
for archaic derived-allele proportion. All intervals here are BED0 half-open.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, math, os, tempfile
from collections import defaultdict, Counter
from pathlib import Path
from bisect import bisect_right

VERSION = '2026-09-15.1'
AFR = {'ESN', 'GWD', 'LWK', 'MSL', 'YRI'}
NEA = {'Altai', 'Vindija', 'Chagyr', 'Chagyrskaya'}
DEN = {'Denisova', 'Denisova25'}


def rows(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        yield from csv.DictReader(f, delimiter='\t')


def union(intervals):
    result = []
    for lo, hi in sorted(intervals):
        if lo < 0 or hi <= lo: raise ValueError('Invalid interval')
        if result and lo <= result[-1][1]: result[-1] = (result[-1][0], max(hi, result[-1][1]))
        else: result.append((lo, hi))
    return result


def subtract(lo, hi, mask):
    i = max(0, bisect_right(mask, (lo, math.inf)) - 1)
    for j in range(i,len(mask)):
        left,right=mask[j]
        if left >= hi: break
        if right <= lo: continue
        if left > lo: yield lo, left
        lo = max(lo, right)
        if lo >= hi: return
    if lo < hi: yield lo, hi


def carrier_mask(calls, members, threshold=.30):
    """Pointwise pooled AFR5 carrier frequency, one vote per person.

    Include zero-call participants in the denominator. This is a carrier
    fraction, not an inferred phased allele frequency (IBDmix is unphased).
    """
    if not members: raise ValueError('No African denominator')
    changes = Counter()
    for sample, intervals in calls.items():
        if sample not in members: raise ValueError('Control sample outside denominator')
        for lo, hi in union(intervals):
            changes[lo] += 1; changes[hi] -= 1
    count = 0; previous = None; result = []
    minimum = math.ceil(threshold * len(members))
    for pos, change in sorted(changes.items()):
        if previous is not None and count >= minimum and pos > previous:
            result.append((previous, pos))
        count += change; previous = pos
    return union(result)


def checksum(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while data := f.read(1024*1024): h.update(data)
    return h.hexdigest()


def load_daf(path):
    """BED plus .json declaring build, reference, percentile and provenance."""
    if path is None: return {}, {'status': 'not_applied_missing_author_mask'}
    path = Path(path)
    meta = json.loads(Path(str(path)+'.json').read_text())
    if meta.get('genome_build') != 'GRCh37' or meta.get('reference') != 'Altai':
        raise ValueError('DAF mask requires genome_build=GRCh37 and reference=Altai')
    if meta.get('percentile') != 99.9 or not meta.get('source'):
        raise ValueError('DAF mask requires percentile=99.9 and a source citation')
    masks = defaultdict(list)
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else path.open()) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')): continue
            ch, lo, hi, *_ = line.split(); ch = ch.removeprefix('chr')
            if ch not in set(map(str,range(1,23))): raise ValueError('DAF mask must be autosomal')
            masks[ch].append((int(lo), int(hi)))
    if not masks: raise ValueError('DAF mask is empty')
    return {ch:union(iv) for ch,iv in masks.items()}, dict(meta, status='applied', path=str(path.resolve()), sha256=checksum(path))


def raw_masks(run, chrom, refs):
    """Use original Neanderthal calls, before the AFR-Denisova subtraction."""
    nea = []; controls = {ref:defaultdict(list) for ref in refs if ref in DEN}
    members = set(); seen_pops = set(); input_paths = []
    units = sorted((run/'samples').glob('*/populations.tsv'))
    if not units: raise ValueError(f'IBDmix sample denominators missing: {run}')
    for table in units:
        populations = list(rows(table)); unit = table.parent.name
        for pop in populations:
            name = pop['population']; sample_file = table.parent/(name+'.txt')
            ids = set(sample_file.read_text().split()); input_paths.append(sample_file)
            if len(ids) != int(pop['n']): raise ValueError(f'Population denominator mismatch: {sample_file}')
            if name in AFR: members.update(ids); seen_pops.add(name)
            # Altai defines the paper comparison. Denisova25 uses the same
            # Altai overlap mask as an explicitly labelled extension.
            for ref in ['Altai'] + sorted(controls):
                if ref != 'Altai' and name not in AFR: continue
                path = run/'raw'/unit/f'{ref}.{name}.raw.txt.gz'
                if not path.is_file(): raise ValueError(f'Original IBDmix calls required by final filter: {path}')
                input_paths.append(path)
                for r in rows(path):
                    ch = r['chrom'].removeprefix('chr'); ch = 'X' if ch=='23' else ch
                    if ch != chrom: continue
                    if r['ID'] not in ids: raise ValueError(f'Wrong population sample: {path}')
                    lo, hi = int(r['start'])-1, int(r['end'])-1
                    score = float(r.get('LOD',r.get('slod','nan')))
                    if not math.isfinite(score): raise ValueError(f'Invalid LOD: {path}')
                    if score < 4 or hi-lo < 50000: continue
                    if ref == 'Altai': nea.append((lo,hi))
                    else: controls[ref][r['ID']].append((lo,hi))
        input_paths.append(table)
    if seen_pops != AFR: raise ValueError(f'All five African controls required: {run}')
    return union(nea), {ref:carrier_mask(data,members) for ref,data in controls.items()}, len(members)


def filter_rows(source, masks, daf, chrom, counts):
    for row in rows(source):
        ref = row['anc']; lo, hi = int(row['start']), int(row['end'])
        if row['chrom'].removeprefix('chr') != chrom: raise ValueError('Mixed chromosome final output')
        counts[ref]['input_segments'] += 1; counts[ref]['input_bp'] += hi-lo
        mask = daf.get(chrom, []) if ref=='Altai' else masks.get(ref, [])
        for left, right in subtract(lo,hi,mask):
            if right-left < 50000:
                counts[ref]['short_fragment_bp'] += right-left
                continue
            result = dict(row, start=left, end=right, length=right-left)
            if (left,right)!=(lo,hi):
                result['score_scope'] = 'parent_native_call_after_final_mask'
                # Parent LOD is retained; site counts do not describe fragments.
                for key in ('sites','positive_lods','negative_lods'):
                    if key in result: result[key] = ''
            counts[ref]['output_segments'] += 1; counts[ref]['output_bp'] += right-left
            yield result


def prepare(source, output_root, meta, daf_path=None):
    """Return a filtered snapshot and audit; native chromosome files stay intact."""
    source = Path(source); run = source.parent.parent
    daf, daf_info = load_daf(daf_path)
    with gzip.open(source, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t'); fields = reader.fieldnames
        first = next(reader, None)
    if first is None: return source, {'status':'empty_native_result', 'daf':daf_info}
    chrom = first['chrom'].removeprefix('chr')
    if chrom not in set(map(str,range(1,23))):
        return source, {'status':'extension_chrX_not_filtered', 'daf':{'status':'not_applicable'}}
    if meta.get('genome_build') not in ('b37','37','GRCh37') or meta.get('coordinate_system')!='0-based-half-open':
        raise ValueError(f'Cell final filters require GRCh37 BED0 calls: {source}')
    refs = set(meta.get('refs','').split())
    if refs & DEN and 'Altai' not in refs:
        raise ValueError('Denisovan filtering requires original Altai calls')
    dependencies = [source, run/'run.meta.tsv', *sorted((run/'raw').glob('*/*.raw.txt.gz')), *sorted((run/'samples').glob('*/*.txt')), *sorted((run/'samples').glob('*/populations.tsv'))]
    signature = dict(version=VERSION, code_sha256=checksum(__file__), daf=daf_info,run_meta=meta,
                     inputs=[[str(p.resolve()),p.stat().st_size,p.stat().st_mtime_ns] for p in dependencies])
    key = hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
    dest = Path(output_root)/key
    result = dest/'segments.tsv.gz'; manifest = dest/'audit.json'
    if result.is_file() and manifest.is_file():
        saved = json.loads(manifest.read_text())
        if saved['signature']==signature and saved['output_sha256']==checksum(result): return result, saved
    Path(output_root).mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.filter-',dir=output_root) as tmp:
        stage = Path(tmp)
        overlap, high, denominator = raw_masks(run,chrom,refs) if refs & DEN else ([],{},0)
        masks = {ref:union(overlap+iv) for ref,iv in high.items()}
        for name,iv in [('Altai_overlap',overlap), *[(r+'_AFR30',iv) for r,iv in high.items()], ('Altai_DAF',daf.get(chrom,[]))]:
            (stage/(name+'.bed')).write_text(''.join(f'{chrom}\t{lo}\t{hi}\n' for lo,hi in iv))
        counts = defaultdict(Counter)
        with gzip.open(stage/'segments.tsv.gz','wt') as f:
            w = csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
            w.writerows(filter_rows(source,masks,daf,chrom,counts))
        audit = dict(signature=signature,status='filtered',daf=daf_info,chrom=chrom,
                     denisovan_overlap='union_of_original_Altai_calls_across_all_populations',
                     african_frequency='pointwise_distinct_carriers_in_pooled_ESN_GWD_LWK_MSL_YRI / all_tested_members >= 0.30',
                     n_african_controls=denominator,reference_scope='Altai/Denisova paper pair; Denisova25 extension; other Neanderthals unchanged',
                     counts=dict(counts),output_sha256=checksum(stage/'segments.tsv.gz'))
        (stage/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
        # A process-specific cache generation avoids replacing active readers.
        if dest.exists():
            dest = Path(str(dest)+f'.{os.getpid()}')
        os.replace(stage,dest)
    return dest/'segments.tsv.gz',audit
