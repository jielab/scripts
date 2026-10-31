#!/usr/bin/env python3
"""Independent HERV assembly, PGEN preparation, cohort analysis and plots."""
from __future__ import annotations

# The shared utility filename is sorted first; use a stable import name in Python.
import importlib.util as _common_import
import sys as _common_sys
from pathlib import Path as _CommonPath
if "allele_common" not in _common_sys.modules:
    _common_spec = _common_import.spec_from_file_location(
        "allele_common", _CommonPath(__file__).with_name("0.common.py"))
    _common_module = _common_import.module_from_spec(_common_spec)
    _common_sys.modules[_common_spec.name] = _common_module
    try:
        _common_spec.loader.exec_module(_common_module)
    except BaseException:
        _common_sys.modules.pop(_common_spec.name, None)
        raise


# ========================================================================
# herv_assembly
# ========================================================================

# Chunked, sequence-based ERV census of a VCF-personalized diploid reference.
# 
# No individual chromosome or whole-genome FASTA is saved. Unresolvable alleles
# are audited and masked with N in reference-coordinate spans (plus a declared
# flank); their unknown sequence lengths cannot be inferred from this model.
import argparse
import bisect
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

METHOD = 'vcf-chunk-repeatmasker-v1'
ERV = {'ERV1', 'ERVK', 'ERVL'}
GROUPS = ['ERV_without_MaLR', 'ERV_with_MaLR', 'all_LTR']


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4 << 20), b''):
            h.update(b)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def union(intervals):
    result = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if result and a <= result[-1][1]:
            result[-1][1] = max(b, result[-1][1])
        else:
            result.append([a, b])
    return result


def intersection(left, right):
    result, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        a, b = max(left[i][0],right[j][0]), min(left[i][1],right[j][1])
        if a < b:
            result.append([a,b])
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return result


def fasta_region(fasta, chrom, start, end):
    if start == end:
        return ''
    p = subprocess.run(['samtools', 'faidx', str(fasta), f'{chrom}:{start+1}-{end}'],
                       check=True, capture_output=True, text=True)
    seq = ''.join(p.stdout.splitlines()[1:]).upper()
    if len(seq) != end-start:
        raise ValueError(f'Incorrect indexed FASTA length: {chrom}:{start}-{end}')
    return seq


def minimal_edit(start, ref, alt):
    while ref and alt and ref[0] == alt[0]:
        start += 1
        ref, alt = ref[1:], alt[1:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return start, start+len(ref), alt


def intersects(a, b, c, d):
    # A zero-width insertion belongs to the base on its right.
    return a < d and max(b, a+1) > c


def resolve_edits(edits, uncertain, length, flank, report):
    """Reject overlapping selected edits together, then close masks over edits."""
    edits.sort(key=lambda e: (e[0], e[1], e[2]))
    masks = list(uncertain)
    accepted = []
    cluster, end = [], -1

    def flush():
        if len(cluster) == 1:
            accepted.append(cluster[0])
        elif cluster:
            masks.append([max(0, cluster[0][0]-flank), min(length, end+flank)])
            for e in cluster:
                report(e[3], 'overlapping_selected_alleles')

    for e in edits:
        if cluster and e[0] >= end:
            flush()
            cluster = []
        cluster.append(e)
        end = max(max(x[1], x[0]+1) for x in cluster)
    flush()
    masks = union(masks)
    # Expanding only to the edit's extent avoids cascades of flanking masks.
    while True:
        kept, expanded = [], list(masks)
        j = 0
        for e in accepted:
            while j < len(masks) and masks[j][1] <= e[0]:
                j += 1
            if j < len(masks) and intersects(e[0], e[1], *masks[j]):
                expanded.append([e[0], min(length, max(e[1], e[0]+1))])
                report(e[3], 'overlaps_uncertain_region')
            else:
                kept.append(e)
        new_masks = union(expanded)
        if new_masks == masks:
            return [e[:3] for e in kept], masks
        accepted, masks = kept, new_masks


def prepare_chrom(args, chrom, length):
    root = args.out / 'variants' / chrom
    root.mkdir(parents=True, exist_ok=True)
    vcf = args.vcf_dir / (chrom + '.vcf.gz')
    stat = vcf.stat()
    inputs = dict(method=METHOD, sample=args.sample, chromosome=chrom, length=length,
                  vcf=str(vcf.resolve()), vcf_bytes=stat.st_size, vcf_mtime_ns=stat.st_mtime_ns,
                  fasta_sha256=args.fasta_sha256, uncertain_flank=args.uncertain_flank)
    receipt = root / 'done.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old['inputs'] != inputs or any(digest(root/n) != h for n, h in old['files'].items()):
            raise ValueError(f'Changed prepared inputs/outputs: {root}')
        return old
    t0 = time.monotonic()
    reference = fasta_region(args.fasta, chrom, 0, length)
    edits, masks = [[], []], [[], []]
    counts = Counter()
    audit_path = root / 'unhandled.tsv.gz'
    # BGZF does many small reads. A bounded temporary compressed VCF copy on
    # Linux avoids the Windows-mount syscall penalty; no sample FASTA is saved.
    vcf_temp = tempfile.TemporaryDirectory(prefix='herv-vcf-', dir=getattr(args,'vcf_scratch',None))
    local_vcf = Path(vcf_temp.name)/(chrom+'.vcf.gz')
    shutil.copyfile(vcf, local_vcf)
    query = ['bcftools', 'query', '-s', args.sample, '-i', 'GT="alt" || GT="mis"',
             '-f', '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\t%FILTER\t[%GT]\n', str(local_vcf)]
    selected_hash = hashlib.sha256()
    with gzip.open(audit_path, 'wt') as audit, open(root/'query.log', 'w') as log:
        writer = csv.writer(audit, delimiter='\t', lineterminator='\n')
        writer.writerow(['chromosome', 'haplotype', 'position_1based', 'id', 'REF', 'ALT', 'GT', 'reason'])
        def report(hap, row, reason):
            counts[f'hap{hap+1}:unhandled:{reason}'] += 1
            writer.writerow([chrom, hap+1, row[1], row[2], row[3], row[4], row[7], reason])
        process = subprocess.Popen(query, stdout=subprocess.PIPE, stderr=log, text=True)
        try:
            for line in process.stdout:
                selected_hash.update(line.encode())
                row = line.rstrip('\n').split('\t')
                if len(row) != 8 or row[0] != chrom:
                    raise ValueError(f'Invalid query record: {row[:3]}')
                counts['selected_records'] += 1
                pos, ref, alts, gt = int(row[1])-1, row[3].upper(), row[4].split(','), row[7]
                info = dict(x.split('=', 1) if '=' in x else (x, True) for x in row[5].split(';'))
                end = max(pos+len(ref), int(info.get('END', pos+len(ref))))
                if pos < 0 or end > length:
                    raise ValueError(f'Out of range variant {chrom}:{pos+1}-{end}')
                alleles = re.split(r'[|/]', gt)
                unphased_het = '/' in gt and len(set(alleles)) > 1
                for hap in range(2):
                    allele = alleles[hap] if len(alleles) == 2 else '.'
                    reason = None
                    if len(alleles) != 2:
                        reason = 'non_diploid_autosomal_GT'
                    elif unphased_het:
                        reason = 'unphased_heterozygote'
                    elif allele == '.':
                        reason = 'missing_GT'
                    elif allele == '0':
                        continue
                    elif row[6] not in ('.', 'PASS'):
                        reason = 'record_FILTER'
                    elif not allele.isdigit() or int(allele) > len(alts):
                        reason = 'invalid_allele_index'
                    elif reference[pos:pos+len(ref)] != ref:
                        reason = 'REF_mismatch'
                    else:
                        alt = alts[int(allele)-1].upper()
                        if re.fullmatch('[ACGT]+', alt) and re.fullmatch('[ACGT]+', ref):
                            a, b, replacement = minimal_edit(pos, ref, alt)
                            if a == b and not replacement:
                                continue
                            edits[hap].append([a, b, replacement, row])
                            counts[f'hap{hap+1}:candidate_literal'] += 1
                            continue
                        if alt == '<DEL>' and 'END' in info and 'IMPRECISE' not in info and \
                                info.get('CIPOS', '0,0') == '0,0' and info.get('CIEND', '0,0') == '0,0' and end > pos+1:
                            # VCF symbolic deletion retains its anchor base.
                            edits[hap].append([pos+1, end, '', row])
                            counts[f'hap{hap+1}:candidate_symbolic_DEL'] += 1
                            continue
                        reason = 'unsupported_or_imprecise_allele'
                    report(hap, row, reason)
                    masks[hap].append([max(0, pos-args.uncertain_flank),
                                       min(length, end+args.uncertain_flank)])
            if process.wait() != 0:
                raise RuntimeError(f'bcftools query failed; see {root}/query.log')
        except BaseException:
            process.kill()
            process.wait()
            raise
        finally:
            process.stdout.close()
            vcf_temp.cleanup()
        for hap in range(2):
            edits[hap], masks[hap] = resolve_edits(edits[hap], masks[hap], length,
                args.uncertain_flank, lambda row, reason: report(hap, row, reason))
            counts[f'hap{hap+1}:applied_edits'] = len(edits[hap])
            counts[f'hap{hap+1}:uncertain_reference_bp'] = sum(b-a for a, b in masks[hap])
            counts[f'hap{hap+1}:uncertain_reference_non_N_bp'] = sum(b-a-reference.count('N', a, b) for a, b in masks[hap])
            with gzip.open(root/f'hap{hap+1}.json.gz', 'wt') as out:
                json.dump(dict(edits=edits[hap], masks=masks[hap]), out, separators=(',', ':'))
    result = dict(inputs=inputs, counts=dict(counts), query_command=query,
                  selected_records_sha256=selected_hash.hexdigest(), seconds=time.monotonic()-t0,
                  completed=now(), files={n: digest(root/n) for n in ['unhandled.tsv.gz', 'hap1.json.gz', 'hap2.json.gz']})
    save(receipt, result)
    print(json.dumps(dict(event='prepared', chromosome=chrom, seconds=round(result['seconds'], 1),
                          records=counts['selected_records'])), flush=True)
    return result


def load_model(root, chrom, hap):
    if hap == 'reference':
        return dict(edits=[], masks=[])
    with gzip.open(root/'variants'/chrom/(hap+'.json.gz'), 'rt') as f:
        return json.load(f)


class Model:
    def __init__(self, edits, masks):
        self.edits = edits
        self.masks = masks
        self.starts = [e[0] for e in edits]
        self.ends = [e[1] for e in edits]
        self.delta = [0]
        for a, b, alt in edits:
            self.delta.append(self.delta[-1]+len(alt)-(b-a))

    def boundary(self, x):
        """Map a reference cut before any insertion at that exact cut."""
        i = bisect.bisect_left(self.starts, x)
        if i and self.ends[i-1] > x:
            raise ValueError('Boundary inside a replacement')
        return x+self.delta[i]

    def expand(self, a, b):
        i = bisect.bisect_left(self.starts, a)
        if i and self.ends[i-1] > a:
            a = self.starts[i-1]
        i = bisect.bisect_left(self.starts, b)
        if i and self.ends[i-1] > b:
            b = self.ends[i-1]
        return a, b

    def render(self, sequence, a, b, chromosome_length):
        # Mask unresolved reference spans without inventing sequence or length.
        seq = bytearray(sequence.encode())
        for x, y in self.masks:
            if x >= b:
                break
            if y > a:
                lo, hi = max(x, a)-a, min(y, b)-a
                seq[lo:hi] = b'N'*(hi-lo)
        sequence = seq.decode()
        out, cursor = [], a
        i = bisect.bisect_left(self.starts, a)
        j = bisect.bisect_left(self.starts, b)
        if b == chromosome_length:
            j = bisect.bisect_right(self.starts, b)
        for x, y, alt in self.edits[i:j]:
            if y > b:
                raise ValueError('Rendering split an edit')
            out.extend((sequence[cursor-a:x-a], alt))
            cursor = y
        out.append(sequence[cursor-a:])
        return ''.join(out)


def chunks(length, models, chunk_size):
    a = 0
    while a < length:
        b = min(a+chunk_size, length)
        while True:
            new_b = max(m.expand(a, b)[1] for m in models)
            if b == new_b:
                break
            b = new_b
        yield a, b
        a = b


def groups(category):
    if not category.startswith('LTR/') and category != 'LTR':
        return []
    family = category.partition('/')[2]
    result = ['all_LTR', 'family:'+family]
    if family in ERV:
        result.extend(['ERV_without_MaLR', 'ERV_with_MaLR'])
    elif family == 'ERVL-MaLR':
        result.append('ERV_with_MaLR')
    return result


def parse_hits(path, seq_length, core_start, core_end):
    hits = defaultdict(list)
    records = 0
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            fields = line.split()
            if not fields or not fields[0].isdigit():
                continue
            if len(fields) < 14:
                raise ValueError('Malformed RepeatMasker output')
            a, b, category = int(fields[5])-1, int(fields[6]), fields[10]
            if not (0 <= a < b <= seq_length):
                raise ValueError('RepeatMasker annotation outside input sequence')
            a, b = max(a, core_start), min(b, core_end)
            if a >= b:
                continue
            records += 1
            for group in groups(category):
                hits[group].append((a-core_start, b-core_start))
    for group in GROUPS:
        hits[group] = hits[group]
    return {g: union(v) for g, v in hits.items()}, records


def run_chunk(args, descriptor, model, signature, verify_only=False):
    chrom, label, a, b, length = descriptor
    folder = args.out/'chunks'/chrom/label/f'{a:09d}-{b:09d}'
    receipt = folder/'done.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old['signature'] != signature or any(digest(folder/n) != h for n, h in old['files'].items()):
            raise ValueError(f'Changed chunk receipt: {folder}')
        return old
    if verify_only:
        raise ValueError(f'Completed analysis is missing a chunk receipt: {receipt}')
    folder.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    start, end = model.expand(max(0, a-args.halo), min(length, b+args.halo))
    # A long deletion near a boundary must not consume the entire halo.
    while start and model.boundary(a)-model.boundary(start) < args.halo:
        start, _ = model.expand(max(0, start-args.halo), end)
    while end < length and model.boundary(end)-model.boundary(b) < args.halo:
        _, end = model.expand(start, min(length, end+args.halo))
    reference = fasta_region(args.fasta, chrom, start, end)
    sequence = model.render(reference, start, end, length)
    offset = model.boundary(start)
    core_start, core_end = model.boundary(a)-offset, model.boundary(b)-offset
    if b == length:
        core_end = len(sequence)
    core = sequence[core_start:core_end]
    if core_end < core_start or len(core) != core_end-core_start:
        raise ValueError('Invalid core mapping')
    records, hits, command, seconds = 0, {g: [] for g in GROUPS}, [], 0
    resources = dict(elapsed_seconds=0, user_seconds=0, system_seconds=0, max_rss_kib=0)
    with tempfile.TemporaryDirectory(prefix='herv-rm-', dir=args.scratch) as temp:
        temp = Path(temp)
        input_fa = temp/'chunk.fa'
        with input_fa.open('w') as f:
            f.write('>chunk\n')
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60]+'\n')
        if sequence and sequence.count('N') != len(sequence):
            command = [str(args.repeatmasker), '-e', 'rmblast', '-pa', '1',
                       '-lib', str(args.library), '-dir', str(temp), str(input_fa)]
            env = dict(os.environ, PATH=str(args.repeatmasker.parent)+os.pathsep+os.environ.get('PATH', ''))
            t_rm = time.monotonic()
            with (folder/'repeatmasker.log').open('w') as log:
                timed = ['/usr/bin/time', '-f',
                         '{"elapsed_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M}',
                         '-o', str(temp/'resources.json')] + command
                # -dir only controls final outputs. RepeatMasker places RM_*
                # working files in cwd, so cwd must also be the Linux scratch.
                subprocess.run(timed, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=temp, check=True)
            seconds = time.monotonic()-t_rm
            resources = json.loads((temp/'resources.json').read_text())
            output = temp/'chunk.fa.out'
            if not output.is_file():
                raise ValueError(f'RepeatMasker did not create .out: {folder}')
            hits, records = parse_hits(output, len(sequence), core_start, core_end)
            with output.open('rb') as src, gzip.open(folder/'repeatmasker.out.gz', 'wb') as dst:
                shutil.copyfileobj(src, dst)
        else:
            (folder/'repeatmasker.log').write_text('All-N chunk: no searchable sequence.\n')
            with gzip.open(folder/'repeatmasker.out.gz', 'wt') as out:
                out.write('No searchable sequence\n')
        rows = []
        with gzip.open(folder/'unions.bed.gz', 'wt') as out:
            for group, intervals in sorted(hits.items()):
                covered = sum(y-x for x, y in intervals)
                non_n = sum(y-x-core.count('N', x, y) for x, y in intervals)
                rows.append(dict(category=group, union_bp=covered, union_non_N_bp=non_n))
                for x, y in intervals:
                    out.write(f'{chrom}\t{x+offset+core_start}\t{y+offset+core_start}\t{group}\n')
        result = dict(signature=signature, chromosome=chrom, label=label, reference_core=[a,b],
                      reference_context=[start,end], haplotype_core=[offset+core_start,offset+core_end],
                      context_haplotype_start=offset, core_offsets=[core_start,core_end],
                      sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                      core_bp=len(core), core_non_N_bp=len(core)-core.count('N'),
                      input_bp=len(sequence), input_non_N_bp=len(sequence)-sequence.count('N'),
                      categories=rows, annotation_records=records, command=command,
                      repeatmasker_seconds=seconds, resources=resources,
                      wall_seconds=time.monotonic()-t0, completed=now())
    result['files'] = {n: digest(folder/n) for n in ['repeatmasker.log','repeatmasker.out.gz','unions.bed.gz']}
    save(receipt, result)
    print(json.dumps(dict(event='annotated', chromosome=chrom, label=label, reference_core=[a,b],
                          core_bp=len(core), seconds=round(result['wall_seconds'],1))), flush=True)
    return result


def prepare(args):
    lengths = {line.split('\t')[0]: int(line.split('\t')[1])
               for line in Path(str(args.fasta)+'.fai').read_text().splitlines()}
    args.fasta_sha256 = digest(args.fasta)
    args.vcf_scratch.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(prepare_chrom, args, chrom, lengths[chrom]) for chrom in args.chromosomes]
        for f in as_completed(futures):
            f.result()
    save(args.out/'prepared.json', dict(method=METHOD, sample=args.sample, chromosomes=args.chromosomes,
         fasta=str(args.fasta), fasta_sha256=args.fasta_sha256, uncertain_flank=args.uncertain_flank,
         policy='Apply phased literal alleles and precise symbolic DEL. Reject selected overlaps; '
                'mask unresolved spans and flanks with N; reference-coordinate lengths of unresolved SVs '
                'are placeholders, not inferred individual sequence lengths.', completed=now()))


def require_library(path):
    if not path.is_file():
        raise ValueError(f'Repeat library not found: {path}; set --library FILE or HERV_LIBRARY '
                         'to the existing repeat-library FASTA')


def verify_completed_results(out, settings, descriptors):
    done = json.loads((out/'analysis.done.json').read_text())
    plan = json.loads((out/'plan.json').read_text())
    result = json.loads((out/'genome_fraction.json').read_text())
    if done['settings'] != settings or result['settings'] != settings or \
       done['chunks'] != len(descriptors) or plan['chunks'] != len(descriptors) or \
       plan['descriptors'] != [list(d) for d in descriptors] or \
       result['completed_chunks'] != dict(Counter(d[1] for d in descriptors)):
        raise ValueError('Completed analysis metadata differs from the verified run')
    if digest(out/'genome_fraction.tsv') != done['results_sha256']:
        raise ValueError('Completed genome_fraction.tsv checksum changed')
    with (out/'genome_fraction.tsv').open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        rows = list(reader)
    if rows != [{key: '' if row[key] is None else str(row[key]) for key in fields}
                for row in result['categories']]:
        raise ValueError('Completed genome_fraction.json differs from the verified TSV')
    for name in ['variant_summary.tsv', 'chromosome_fraction.tsv', 'boundary_audit.tsv']:
        if not (out/name).is_file():
            raise ValueError(f'Completed analysis is missing {name}')


def run(args):
    require_library(args.library)
    prepared = json.loads((args.out/'prepared.json').read_text())
    if digest(args.fasta) != prepared['fasta_sha256']:
        raise ValueError('Reference changed after variant preparation')
    # Preserve the environment's bin/ directory even when the executable is a
    # symlink into share/RepeatMasker; its Perl dependencies live in this env.
    args.repeatmasker = args.repeatmasker.absolute()
    args.library = args.library.resolve()
    lengths = {line.split('\t')[0]: int(line.split('\t')[1])
               for line in Path(str(args.fasta)+'.fai').read_text().splitlines()}
    settings = dict(method=METHOD, sample=prepared['sample'], fasta_sha256=prepared['fasta_sha256'],
                    library=str(args.library), library_sha256=digest(args.library),
                    repeatmasker=str(args.repeatmasker), repeatmasker_sha256=digest(args.repeatmasker),
                    chunk_bp=args.chunk_bp, halo=args.halo, uncertain_flank=prepared['uncertain_flank'],
                    parameters=['-e','rmblast','-pa','1'],
                    code_sha256=digest(__file__))
    env = dict(os.environ, PATH=str(args.repeatmasker.parent)+os.pathsep+os.environ.get('PATH',''))
    settings['rmblast_version'] = subprocess.check_output(
        [str(args.repeatmasker.parent/'rmblastn'),'-version'],text=True,env=env).strip()
    manifest = args.out/'run.json'
    reuse_complete = False
    if manifest.exists():
        previous = json.loads(manifest.read_text())
        reuse_complete = previous.get('status') == 'complete' and (args.out/'analysis.done.json').is_file()
        # A relocated library is identical only when its SHA-256 still matches.
        # Completed results can survive code packaging changes, but must pass
        # summary, preparation and every chunk's signature/output checks below.
        ignored = {'library'} | ({'code_sha256'} if reuse_complete else set())
        changed = sorted(k for k in settings.keys() | previous['settings'].keys()
                         if k not in ignored and settings.get(k) != previous['settings'].get(k))
        if changed or previous['scope'] != prepared['chromosomes']:
            raise ValueError('Run settings changed ('+', '.join(changed or ['scope'])+
                             '); use a new output directory')
        # Retain the historical settings for receipt signatures and provenance;
        # any new chunks still read args.library at its current location.
        settings = previous['settings']
    if not manifest.exists():
        save(manifest, dict(settings=settings, started=now(), started_unix=time.time(), jobs=args.jobs,
                           cpu_affinity=sorted(os.sched_getaffinity(0)),
                           max_search_threads=4*args.jobs, nice=os.getpriority(os.PRIO_PROCESS,0),
                           scope=prepared['chromosomes'], status='running'))
    if not reuse_complete:
        args.scratch.mkdir(parents=True, exist_ok=True)
    jobs = []
    for chrom in prepared['chromosomes']:
        receipt = json.loads((args.out/'variants'/chrom/'done.json').read_text())
        if any(receipt['inputs'][key] != prepared[key]
               for key in ['sample','fasta_sha256','uncertain_flank']) or receipt['inputs']['length'] != lengths[chrom]:
            raise ValueError('Inconsistent chromosome preparation metadata')
        if any(digest(args.out/'variants'/chrom/n) != h for n,h in receipt['files'].items()):
            raise ValueError('Prepared variants changed')
        models = {label: Model(**load_model(args.out, chrom, label)) for label in ['reference','hap1','hap2']}
        for a,b in chunks(lengths[chrom], list(models.values()), args.chunk_bp):
            for label, model in models.items():
                desc = (chrom,label,a,b,lengths[chrom])
                sig = hashlib.sha256(json.dumps([settings,receipt['files'],desc],sort_keys=True).encode()).hexdigest()
                jobs.append((desc, model, sig))
    if reuse_complete:
        verify_completed_results(args.out, settings, [j[0] for j in jobs])
        print(f'[herv] VERIFY completed analysis: {len(jobs)} chunk receipts and outputs', flush=True)
    else:
        save(args.out/'plan.json', dict(chunks=len(jobs), reference_autosomal_bp=sum(lengths[c] for c in prepared['chromosomes']),
             descriptors=[j[0] for j in jobs]))
    with ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_chunk, args, *job, verify_only=reuse_complete) for job in jobs]
        completed = 0
        try:
            for f in as_completed(futures):
                f.result()
                completed += 1
                if not reuse_complete:
                    save(args.out/'progress.json', dict(completed_chunks=completed, total_chunks=len(jobs), updated=now()))
        except BaseException:
            for f in futures:
                f.cancel()
            raise
    if reuse_complete:
        print(f'[herv] REUSE completed analysis: verified {len(jobs)} chunks; output={args.out}', flush=True)
        return
    summarize(args)


def summarize(args):
    plan = json.loads((args.out/'plan.json').read_text())
    run_info = json.loads((args.out/'run.json').read_text())
    totals = defaultdict(Counter)
    chromosome_totals = defaultdict(Counter)
    annotation_seconds = Counter()
    cpu_seconds = Counter()
    counts = Counter()
    receipts = []
    by_contig = defaultdict(list)
    for chrom,label,a,b,length in plan['descriptors']:
        receipt = args.out/'chunks'/chrom/label/f'{a:09d}-{b:09d}'/'done.json'
        if not receipt.is_file():
            raise ValueError(f'Incomplete analysis: {receipt}')
        r = json.loads(receipt.read_text())
        if any(digest(receipt.parent/n) != h for n,h in r['files'].items()):
            raise ValueError(f'Changed annotation outputs: {receipt.parent}')
        receipts.append(r)
        by_contig[(chrom,label)].append((r,receipt.parent))
        counts[label] += 1
        annotation_seconds[label] += r['repeatmasker_seconds']
        cpu_seconds[label] += r['resources']['user_seconds']+r['resources']['system_seconds']
        for row in r['categories']:
            t = totals[(label,row['category'])]
            t['union_bp'] += row['union_bp']
            t['union_non_N_bp'] += row['union_non_N_bp']
            chromosome_totals[(chrom,label,row['category'])].update(
                union_bp=row['union_bp'],union_non_N_bp=row['union_non_N_bp'])
        totals[(label,'denominator')]['genome_bp'] += r['core_bp']
        totals[(label,'denominator')]['genome_non_N_bp'] += r['core_non_N_bp']
        chromosome_totals[(chrom,label,'denominator')].update(genome_bp=r['core_bp'],genome_non_N_bp=r['core_non_N_bp'])
    # Matched reference controls remove the same uncertainty spans as each
    # haplotype, separating masking effects from the effects of resolved edits.
    variant_rows = []
    for chrom in run_info['scope']:
        prep = json.loads((args.out/'variants'/chrom/'done.json').read_text())
        length = prep['inputs']['length']
        sequence = fasta_region(prep['inputs'].get('fasta',
            json.loads((args.out/'prepared.json').read_text())['fasta']),chrom,0,length)
        ref_hits = defaultdict(list)
        for folder in sorted((args.out/'chunks'/chrom/'reference').iterdir()):
            if not folder.is_dir():
                continue
            with gzip.open(folder/'unions.bed.gz','rt') as f:
                for line in f:
                    _,a,b,group = line.rstrip().split('\t')
                    ref_hits[group].append([int(a),int(b)])
        for group in ref_hits:
            ref_hits[group] = union(ref_hits[group])
        for group in GROUPS:
            ref_hits.setdefault(group,[])
        for hap in ['hap1','hap2']:
            masks = load_model(args.out,chrom,hap)['masks']
            label = 'reference_callable_'+hap
            denom = totals[(label,'denominator')]
            denom['genome_bp'] += length
            denom['genome_non_N_bp'] += length-sequence.count('N')-sum(b-a-sequence.count('N',a,b) for a,b in masks)
            with gzip.open(args.out/'variants'/chrom/(hap+'.uncertain.bed.gz'),'wt') as f:
                for a,b in masks:
                    f.write(f'{chrom}\t{a}\t{b}\n')
            for group,intervals in ref_hits.items():
                removed = intersection(intervals,masks)
                t = totals[(label,group)]
                t['union_bp'] += sum(b-a for a,b in intervals)-sum(b-a for a,b in removed)
                t['union_non_N_bp'] += sum(b-a-sequence.count('N',a,b) for a,b in intervals)-sum(b-a-sequence.count('N',a,b) for a,b in removed)
            variant_rows.append(dict(chromosome=chrom,haplotype=hap,
                applied_edits=prep['counts'].get(hap+':applied_edits',0),
                unhandled_records=sum(v for k,v in prep['counts'].items() if k.startswith(hap+':unhandled:')),
                uncertain_reference_bp=prep['counts'].get(hap+':uncertain_reference_bp',0),
                uncertain_reference_non_N_bp=prep['counts'].get(hap+':uncertain_reference_non_N_bp',0)))
    with (args.out/'variant_summary.tsv').open('w') as f:
        writer = csv.DictWriter(f,list(variant_rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(variant_rows)
    # Quantify annotation disagreement in the shared inner half of each halo.
    # This reports sensitivity to batching rather than claiming equivalence to
    # an unfragmented RepeatMasker run (RepeatMasker also fragments internally).
    boundary_rows = []
    for (chrom,label), pieces in by_contig.items():
        pieces.sort(key=lambda x:x[0]['reference_core'][0])
        previous = None
        previous_hits = None
        for current,folder in pieces:
            h0,_ = parse_hits(folder/'repeatmasker.out.gz',current['input_bp'],0,current['input_bp'])
            shifted = {g:[[a+current['context_haplotype_start'],b+current['context_haplotype_start']]
                          for a,b in h0[g]] for g in GROUPS}
            if previous:
                boundary = current['haplotype_core'][0]
                if previous['haplotype_core'][1] != boundary:
                    raise ValueError('Gap or duplicate sequence between reconstructed cores')
                radius = run_info['settings']['halo']//2
                lo = max(boundary-radius,current['context_haplotype_start'],previous['context_haplotype_start'])
                hi = min(boundary+radius,current['context_haplotype_start']+current['input_bp'],
                         previous['context_haplotype_start']+previous['input_bp'])
                for group in GROUPS:
                    left = intersection(previous_hits[group],[[lo,hi]]) if lo < hi else []
                    right = intersection(shifted[group],[[lo,hi]]) if lo < hi else []
                    common = intersection(left,right)
                    boundary_rows.append(dict(chromosome=chrom,label=label,haplotype_boundary=boundary,
                        category=group,compared_bp=max(0,hi-lo),
                        disagreement_bp=sum(b-a for a,b in left+right)-2*sum(b-a for a,b in common)))
            previous,previous_hits = current,shifted
    with (args.out/'boundary_audit.tsv').open('w') as f:
        keys=['chromosome','label','haplotype_boundary','category','compared_bp','disagreement_bp']
        writer=csv.DictWriter(f,keys,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(boundary_rows)
    for label in ['hap1','hap2','reference_callable_hap1','reference_callable_hap2']:
        for (source,category), values in list(totals.items()):
            if source == label:
                target = 'reference_callable_diploid' if source.startswith('reference_') else 'diploid'
                totals[(target,category)].update(values)
    rows = []
    for (label,category), values in sorted(totals.items()):
        if category == 'denominator':
            continue
        row = dict(label=label, category=category, **dict(values), **dict(totals[(label,'denominator')]))
        row['percent_genome'] = 100*row['union_bp']/row['genome_bp'] if row['genome_bp'] else None
        row['percent_non_N'] = 100*row['union_non_N_bp']/row['genome_non_N_bp'] if row['genome_non_N_bp'] else None
        rows.append(row)
    keys = ['label','category','union_bp','genome_bp','percent_genome','union_non_N_bp','genome_non_N_bp','percent_non_N']
    with (args.out/'genome_fraction.tsv').open('w') as f:
        writer = csv.DictWriter(f,keys,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    chromosome_rows=[]
    for (chrom,label,category),values in sorted(chromosome_totals.items()):
        if category=='denominator':
            continue
        row=dict(chromosome=chrom,label=label,category=category,**dict(values),
                 **dict(chromosome_totals[(chrom,label,'denominator')]))
        row['percent_genome']=100*row['union_bp']/row['genome_bp'] if row['genome_bp'] else None
        row['percent_non_N']=100*row['union_non_N_bp']/row['genome_non_N_bp'] if row['genome_non_N_bp'] else None
        chromosome_rows.append(row)
    with (args.out/'chromosome_fraction.tsv').open('w') as f:
        writer=csv.DictWriter(f,['chromosome']+keys,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(chromosome_rows)
    ended_unix = run_info.get('completed_unix',time.time())
    elapsed = ended_unix-run_info['started_unix']
    cohort_job_seconds = annotation_seconds['reference'] + 3202*(annotation_seconds['hap1']+annotation_seconds['hap2'])
    result = dict(scope='autosomes only; VCF-conditioned reference models, not independent assemblies',
                  categories=rows, completed_chunks=dict(counts), repeatmasker_job_seconds=dict(annotation_seconds),
                  measured_cpu_seconds=dict(cpu_seconds),
                  boundary_disagreement_bp={label:sum(r['disagreement_bp'] for r in boundary_rows
                      if r['label']==label and r['category']=='ERV_with_MaLR') for label in ['reference','hap1','hap2']},
                  cohort_projection=dict(samples=3202, reference_runs=1, haplotypes_per_sample=2,
                      annotation_job_hours=cohort_job_seconds/3600,
                      annotation_wall_days_at_same_jobs=cohort_job_seconds/(run_info['jobs']*86400),
                      annotation_cpu_core_hours=(cpu_seconds['reference']+3202*(cpu_seconds['hap1']+cpu_seconds['hap2']))/3600,
                      caveat='Linear projection from HG00096, assuming similar workloads and resource contention. '
                             'Excludes installation, downloads, VCF preparation and scheduling overhead.'),
                  elapsed_wall_seconds=elapsed, completed=now(), settings=run_info['settings'],
                  interpretation='Unresolved variation was masked; unknown inserted sequence/length is not recovered. '
                    'Percent including N uses model lengths; percent_non_N excludes reference gaps and masked uncertainty.')
    save(args.out/'genome_fraction.json', result)
    save(args.out/'analysis.done.json', dict(completed=now(), chunks=len(receipts),
         results_sha256=digest(args.out/'genome_fraction.tsv'), settings=run_info['settings']))
    run_info['status']='complete'
    run_info['completed']=now()
    run_info['completed_unix']=ended_unix
    save(args.out/'run.json',run_info)
    print(json.dumps(dict(event='complete', seconds=elapsed, output=str(args.out))), flush=True)


def main_assembly():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command',required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--fasta',type=Path,required=True)
    prep.add_argument('--vcf-dir',type=Path,required=True)
    prep.add_argument('--sample',required=True)
    prep.add_argument('--chromosomes',default=','.join('chr'+str(i) for i in range(1,23)))
    prep.add_argument('--uncertain-flank',type=int,default=1000)
    prep.add_argument('--jobs',type=int,default=2)
    prep.add_argument('--vcf-scratch',type=Path,default=Path('/tmp/metag-herv-vcf'))
    prep.add_argument('--out',type=Path,required=True)
    ann = sub.add_parser('run')
    ann.add_argument('--fasta',type=Path,required=True)
    ann.add_argument('--library',type=Path,required=True)
    ann.add_argument('--repeatmasker',type=Path,required=True)
    ann.add_argument('--chunk-bp',type=int,default=5000000)
    ann.add_argument('--halo',type=int,default=20000)
    ann.add_argument('--jobs',type=int,default=3)
    ann.add_argument('--scratch',type=Path,default=Path('/tmp/metag-herv-assembly'))
    ann.add_argument('--out',type=Path,required=True)
    summary = sub.add_parser('summarize')
    summary.add_argument('--out',type=Path,required=True)
    args = p.parse_args()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/(args.command+'.lock')).open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.command == 'prepare':
            args.chromosomes = args.chromosomes.split(',')
            if len(set(args.chromosomes)) != len(args.chromosomes) or any(not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])',c) for c in args.chromosomes):
                p.error('Specify unique autosomes chr1..chr22')
            prepare(args)
        elif args.command == 'run':
            if args.chunk_bp <= 0 or args.halo < 0 or args.jobs < 1:
                p.error('Invalid chunk size, halo or jobs')
            run(args)
        else:
            summarize(args)


# ========================================================================
# herv_pgen
# ========================================================================

# Prepare the existing HERV sequence models directly from phased PLINK 2 files.
# 
# Only a compact, temporary single-sample variant stream bridges to the unchanged
# VCF preparation engine. No individual chromosome/genome FASTA is materialized.
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import fcntl
import gzip
import importlib.metadata
import io
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time
from types import SimpleNamespace

engine = sys.modules[__name__]


def source_files(root, chrom):
    base = root/chrom
    pvar = Path(str(base)+'.pvar.zst')
    if not pvar.is_file():
        pvar = Path(str(base)+'.pvar')
    files = {'pgen': Path(str(base)+'.pgen'), 'pvar': pvar,
             'psam': Path(str(base)+'.psam')}
    for path in files.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    pgi = Path(str(base)+'.pgen.pgi')
    if pgi.exists():
        files['pgi'] = pgi
    return files


def file_state(path):
    stat = path.stat()
    return dict(path=str(path.resolve()), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def sample_index(psam, sample):
    found, total = [], 0
    with psam.open() as stream:
        header = stream.readline().lstrip('#').split()
        if 'IID' not in header:
            raise ValueError(f'PSAM requires an IID header: {psam}')
        column = header.index('IID')
        for line in stream:
            if not line.strip():
                continue
            row = line.split()
            if len(row) != len(header):
                raise ValueError(f'Malformed PSAM: {psam}')
            if row[column] == sample:
                found.append(total)
            total += 1
    if len(found) != 1:
        raise ValueError(f'Expected exactly one IID={sample} in {psam}; found {len(found)}')
    return found[0], total


@contextmanager
def pvar_text(path):
    if path.suffix == '.zst':
        import zstandard
        with path.open('rb') as raw, zstandard.ZstdDecompressor().stream_reader(raw) as reader:
            with io.TextIOWrapper(reader) as stream:
                yield stream
    else:
        with path.open() as stream:
            yield stream


def compact_info(info):
    # These are all INFO fields consulted by the unchanged reconstruction policy.
    return ';'.join(x for x in info.split(';')
                    if x.partition('=')[0] in {'END', 'IMPRECISE', 'CIPOS', 'CIEND'}) or '.'


def export_selected(files, sample, chrom, length, output):
    import numpy as np
    import pgenlib
    index, total = sample_index(files['psam'], sample)
    counts = Counter()
    # The minimal native PVAR reader supplies multiallelic offsets, while the
    # streaming text pass preserves FILTER and structural-variant metadata.
    with pgenlib.PvarReader(os.fsencode(files['pvar']), omit_chrom=True, omit_pos=True) as variants, \
         pgenlib.PgenReader(os.fsencode(files['pgen']), raw_sample_ct=total,
             sample_subset=np.asarray([index], dtype=np.uint32), pvar=variants) as genotypes, \
         pvar_text(files['pvar']) as source, gzip.open(output, 'wt', compresslevel=1) as dest:
        if not genotypes.hardcall_phase_present():
            raise ValueError(f'PGEN contains no hardcall phase information: {files["pgen"]}; use the phased VCF')
        n = variants.get_variant_ct()
        if n != genotypes.get_variant_ct():
            raise ValueError('PVAR/PGEN variant count mismatch')
        columns = None
        for line in source:
            if line.startswith('#CHROM'):
                columns = line.lstrip('#').rstrip('\n').split('\t')
                break
            if not line.startswith('##'):
                raise ValueError('PVAR requires a #CHROM header')
        required = {'CHROM', 'POS', 'ID', 'REF', 'ALT'}
        if columns is None or not required <= set(columns):
            raise ValueError('Missing required PVAR columns')
        lookup = {name: columns.index(name) for name in columns}
        dest.write('##fileformat=VCFv4.3\n'+f'##contig=<ID={chrom},length={length}>\n')
        dest.write('##FILTER=<ID=PGEN_FILTERED,Description="Non-PASS FILTER retained from PVAR">\n')
        for name, number, kind in [('END', '1', 'Integer'), ('IMPRECISE', '0', 'Flag'),
                                   ('CIPOS', '2', 'Integer'), ('CIEND', '2', 'Integer')]:
            dest.write(f'##INFO=<ID={name},Number={number},Type={kind},Description="Preserved from PVAR">\n')
        dest.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased PGEN hardcalls">\n')
        dest.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')
        batch = 8192
        alleles = np.empty((batch, 2), dtype=np.int32)
        phase = np.empty((batch, 1), dtype=np.uint8)
        previous_pos = 0
        for ordinal, line in enumerate(source):
            if ordinal >= n or not line.strip() or line.startswith('#'):
                raise ValueError('PVAR row count/order differs from the native reader')
            if ordinal % batch == 0:
                genotypes.read_alleles_and_phasepresent_range(ordinal, min(ordinal+batch, n), alleles, phase)
            row = line.rstrip('\n').split('\t')
            if len(row) != len(columns):
                raise ValueError(f'Malformed PVAR row {ordinal+1}')
            position = int(row[lookup['POS']])
            row_chrom = row[lookup['CHROM']]
            if ('chr'+row_chrom.removeprefix('chr')) != chrom or not previous_pos <= position <= length:
                raise ValueError(f'PVAR chromosome/position/order mismatch at row {ordinal+1}')
            previous_pos = position
            counts['total_variants'] += 1
            a, b = map(int, alleles[ordinal % batch])
            if a == b == 0:
                continue
            phased = bool(phase[ordinal % batch, 0])
            if a < 0 or b < 0:
                counts['missing_GT'] += 1
            elif a != b:
                counts['phased_heterozygotes' if phased else 'unphased_heterozygotes'] += 1
            else:
                counts['homozygous_alternate'] += 1
            # Homozygous calls have no haplotype ambiguity even when PGEN does
            # not retain their original VCF separator. Missing calls stay missing.
            separator = '|' if phased or (a == b and a >= 0) else '/'
            gt = separator.join('.' if x < 0 else str(x) for x in (a, b))
            info = compact_info(row[lookup['INFO']]) if 'INFO' in lookup else '.'
            filt = row[lookup['FILTER']] if 'FILTER' in lookup else '.'
            # FILTER names can be arbitrary; the existing engine only needs
            # PASS/dot vs failure, so declare one deterministic failure code.
            if filt not in ('.', 'PASS'):
                filt = 'PGEN_FILTERED'
            dest.write('\t'.join([chrom, str(position), row[lookup['ID']], row[lookup['REF']],
                row[lookup['ALT']], '.', filt, info, 'GT', gt])+'\n')
            counts['selected_records'] += 1
        if counts['total_variants'] != n:
            raise ValueError('PVAR ended before all PGEN variants were read')
    return dict(counts)


def prepare_pgen_chrom(args, chrom, length):
    files = source_files(args.pfile_dir, chrom)
    sources = {name: file_state(path) for name, path in files.items()}
    inputs = dict(method=engine.METHOD, input_source='pgen', sample=args.sample,
        chromosome=chrom, length=length, fasta_sha256=args.fasta_sha256,
        uncertain_flank=args.uncertain_flank, sources=sources,
        preparer_sha256=engine.digest(__file__), pgenlib_version=importlib.metadata.version('pgenlib'))
    root = args.out/'variants'/chrom
    receipt = root/'done.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old['inputs'] != inputs or any(engine.digest(root/n) != h for n,h in old['files'].items()):
            raise ValueError(f'Changed PGEN inputs/outputs: {root}')
        return old
    t0 = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='herv-pgen-', dir=args.scratch) as folder:
        temp = Path(folder)
        # PGEN performs many small reads; staging also avoids WSL/NTFS latency.
        # At most --jobs chromosome filesets are staged and each is removed
        # after its compact models are published.
        staged = {}
        for name, source in files.items():
            staged[name] = temp/source.name
            shutil.copyfile(source, staged[name])
        staging_seconds = time.monotonic()-t0
        selected = temp/(chrom+'.vcf.gz')
        export_counts = export_selected(staged, args.sample, chrom, length, selected)
        export_seconds = time.monotonic()-t0
        if sources != {name: file_state(path) for name, path in files.items()}:
            raise ValueError('PGEN input files changed during reading')
        shared = SimpleNamespace(out=temp/'models', vcf_dir=temp, fasta=args.fasta,
            fasta_sha256=args.fasta_sha256, sample=args.sample, uncertain_flank=args.uncertain_flank,
            vcf_scratch=args.scratch)
        result = engine.prepare_chrom(shared, chrom, length)
        result.update(inputs=inputs, pgen_export_counts=export_counts, export_seconds=export_seconds,
                      staging_seconds=staging_seconds,
                      seconds=time.monotonic()-t0, completed=engine.now())
        result['query_command_note'] = 'Temporary compact PGEN hardcall export; removed after preparation.'
        root.mkdir(parents=True, exist_ok=True)
        source_root = shared.out/'variants'/chrom
        for name in [*result['files'], 'query.log']:
            shutil.copy2(source_root/name, root/name)
        # Completion is published only after every model/audit file is copied.
        engine.save(receipt, result)
    print(json.dumps(dict(event='prepared_pgen', chromosome=chrom,
        seconds=round(result['seconds'], 2), export_counts=export_counts)), flush=True)
    return result


def main_pgen():
    os.environ['PATH'] = str(Path(sys.executable).parent)+os.pathsep+os.environ.get('PATH', '')
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['prepare'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--pfile-dir', type=Path, required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--chromosomes', default=','.join('chr'+str(i) for i in range(1,23)))
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--uncertain-flank', type=int, default=1000)
    p.add_argument('--scratch', type=Path, default=Path('/tmp/metag-herv-pgen'))
    args = p.parse_args()
    args.chromosomes = args.chromosomes.split(',')
    if args.jobs < 1 or args.uncertain_flank < 0 or len(set(args.chromosomes)) != len(args.chromosomes) or \
       any(not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])', c) for c in args.chromosomes):
        p.error('Specify positive jobs, nonnegative flank and unique chromosomes chr1..chr22')
    args.out = args.out.resolve()
    args.fasta = args.fasta.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    args.scratch.mkdir(parents=True, exist_ok=True)
    with (args.out/'prepare.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        lengths = {r.split('\t')[0]: int(r.split('\t')[1])
                   for r in Path(str(args.fasta)+'.fai').read_text().splitlines()}
        args.fasta_sha256 = engine.digest(args.fasta)
        # PVAR parsing is CPU-bound Python work; separate processes let the
        # default two chromosome preparations use both allocated CPUs.
        with ProcessPoolExecutor(args.jobs) as pool:
            futures = [pool.submit(prepare_pgen_chrom, args, c, lengths[c]) for c in args.chromosomes]
            for future in futures:
                future.result()
        engine.save(args.out/'prepared.json', dict(method=engine.METHOD, input_source='pgen',
            sample=args.sample, chromosomes=args.chromosomes, fasta=str(args.fasta),
            fasta_sha256=args.fasta_sha256, uncertain_flank=args.uncertain_flank,
            pfile_dir=str(args.pfile_dir.resolve()), preparer_sha256=engine.digest(__file__),
            policy='Same reconstruction/masking policy as herv_assembly.prepare_chrom; '
                   'PGEN hardcalls and phase, PVAR REF/ALT/FILTER/SV INFO. No dosage rounding.',
            completed=engine.now()))


# ========================================================================
# herv
# ========================================================================

# Individual-genome HERV workflow, cohort burden and plots.





# ========================================================================
# herv assembly job
# ========================================================================

# Run/resume individual autosomal ERV annotation with reference and both haplotypes.
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys



def positive(value):
    number=int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('Must be a positive integer')
    return number


def component(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',value):
        raise argparse.ArgumentTypeError('Use a single name containing letters, digits, dots, underscores or hyphens')
    return value


def parser():
    p=argparse.ArgumentParser(description=__doc__,allow_abbrev=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='Run: ./2.herv.sh. Advanced annotation commands: '
               'python f/1.herv.py assembly {prepare,run,summarize} --help.')
    p.add_argument('--out','--dir-out',type=Path,
        help='Override results directory; otherwise /mnt/d/analysis/allele3/herv/<gen>/<sample>')
    p.add_argument('--fasta',type=Path,default=Path('/mnt/f/gen/1hgp/38/GRCH38.fasta.gz'),help='Indexed reference FASTA')
    p.add_argument('--library',type=Path,default=Path('/mnt/f/annot/allele3/dfam_4.0_human/human.curated.fa'),help='Repeat library')
    p.add_argument('--repeatmasker',type=Path,
        help='RepeatMasker executable; reuse checkpoint path, then discover allele3-repeatmasker or legacy metag-repeatmasker')
    p.add_argument('--gen',type=component,default='1kg',help='Dataset name used in the default output directory')
    p.add_argument('--sample',type=component,default='HG00096',help='VCF sample or PGEN IID; scope is chr1-22')
    p.add_argument('--input-format',choices=['auto','vcf','pgen'],default='auto',
        help='Reuse the prepared source; for a new sample auto prefers a complete PGEN fileset and available reader')
    p.add_argument('--vcf-dir','--dir-vcf',type=Path,default=Path('/mnt/f/gen/1kg/38/vcf'),help='VCFs used if preparation is missing')
    p.add_argument('--pfile-dir','--dir-pfile',type=Path,default=Path('/mnt/f/gen/1kg/38/pfile'),help='Chromosome PGEN/PVAR/PSAM filesets')
    p.add_argument('--pgen-python',type=Path,help='Python with numpy, pgenlib and zstandard; autodetected when omitted')
    p.add_argument('--jobs',type=positive,default=8,help='Concurrent RepeatMasker chunks')
    p.add_argument('--cores',type=positive,default=16,help='Limit CPU affinity to the first available CPUs')
    p.add_argument('--chunk-bp',type=positive,default=5000000,help='Reference core chunk length')
    p.add_argument('--halo',type=int,default=20000,help='Minimum flanking context in bp')
    p.add_argument('--scratch','--dir-scratch',type=Path,default=Path('/tmp/allele3-herv-assembly'),help='Temporary chunk directory')
    p.add_argument('--dry-run',action='store_true',help='Print commands without writing files or running analysis')
    return p


def pgen_runtime(explicit=None):
    if explicit:
        candidates=[str(explicit)]
    elif os.environ.get('HERV_PGEN_PYTHON') or os.environ.get('ALLELE_HERV_PGEN_PYTHON') or os.environ.get('METAG_PGEN_PYTHON'):
        candidates=[os.environ.get('HERV_PGEN_PYTHON') or os.environ.get('ALLELE_HERV_PGEN_PYTHON') or os.environ['METAG_PGEN_PYTHON']]
    else:
        root=conda_root()
        candidates=[sys.executable,str(root/'envs/herv-pgen/bin/python'),str(root/'envs/allele3-herv-pgen/bin/python'),
                    str(root/'envs/allele3-evolution/bin/python'),str(root/'envs/metag-herv-pgen/bin/python')]
    for candidate in dict.fromkeys(candidates):
        if not Path(candidate).is_file():
            continue
        try:
            subprocess.run([candidate,'-c',
                'import numpy,pgenlib,zstandard; '
                'assert hasattr(pgenlib,"PvarReader"); '
                'assert hasattr(pgenlib.PgenReader,"read_alleles_and_phasepresent_range")'],
                check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=15)
            return candidate
        except (OSError,subprocess.SubprocessError):
            continue
    return None


def pgen_available(root):
    return all((root/f'chr{i}.pgen').is_file() and (root/f'chr{i}.psam').is_file()
        and ((root/f'chr{i}.pvar.zst').is_file() or (root/f'chr{i}.pvar').is_file()) for i in range(1,23))


def conda_root():
    return Path(os.environ.get('HERV_CONDA_ROOT') or os.environ.get('ALLELE_CONDA_ROOT') or
                os.environ.get('METAG_CONDA_ROOT') or str(Path.home()/'anaconda3'))


def repeatmasker_path(args):
    if args.repeatmasker:
        return args.repeatmasker
    if os.environ.get('HERV_REPEATMASKER') or os.environ.get('ALLELE_HERV_REPEATMASKER'):
        return Path(os.environ.get('HERV_REPEATMASKER') or os.environ['ALLELE_HERV_REPEATMASKER'])
    # Executable paths are part of the existing run signature. Moving results
    # must not silently select a different environment and invalidate receipts.
    checkpoint=args.out/'run.json'
    if checkpoint.is_file():
        return Path(json.loads(checkpoint.read_text())['settings']['repeatmasker'])
    installed = shutil.which('RepeatMasker')
    if installed:
        return Path(installed)
    candidates=[conda_root()/'envs'/env/'bin/RepeatMasker'
                for env in ['herv-repeatmasker','allele3-repeatmasker','metag-repeatmasker']]
    return next((p for p in candidates if p.is_file()),candidates[0])


def commands(args):
    engine=[sys.executable,str(Path(__file__)), 'assembly']
    common=['--out',str(args.out),'--fasta',str(args.fasta.resolve())]
    prepared=args.out/'prepared.json'
    steps=[]
    if prepared.exists():
        saved=json.loads(prepared.read_text())
        if saved['sample'] != args.sample or saved['chromosomes'] != ['chr'+str(i) for i in range(1,23)]:
            raise ValueError('Prepared sample/scope differs from requested sample and chr1-22; use a different --out')
        args.effective_input=saved.get('input_source','vcf')
        if args.input_format not in ('auto',args.effective_input):
            raise ValueError('Prepared input format differs; use a separate --out for a PGEN/VCF comparison')
    else:
        complete=pgen_available(args.pfile_dir)
        runtime=pgen_runtime(args.pgen_python) if args.input_format!='vcf' and complete else None
        if args.input_format=='pgen' and (not complete or not runtime):
            raise ValueError('PGEN requires chr1-22 .pgen/.pvar[.zst]/.psam and a --pgen-python '
                             'with numpy, pgenlib and zstandard')
        args.effective_input='pgen' if args.input_format!='vcf' and complete and runtime else 'vcf'
        if args.effective_input=='pgen':
            steps.append([runtime,str(Path(__file__)), 'pgen', 'prepare']+common+
                ['--pfile-dir',str(args.pfile_dir.resolve()),'--sample',args.sample,'--jobs','2',
                 '--scratch',str(args.scratch.resolve()/'pgen-prepare')])
        else:
            steps.append(engine+['prepare']+common+['--vcf-dir',str(args.vcf_dir.resolve()),
                '--sample',args.sample,'--jobs','2','--vcf-scratch',str(args.scratch.resolve()/'vcf-prepare')])
    steps.append(engine+['run']+common+[
        '--library',str(args.library.resolve()),'--repeatmasker',str(repeatmasker_path(args).absolute()),
        '--jobs',str(args.jobs),'--chunk-bp',str(args.chunk_bp),'--halo',str(args.halo),
        '--scratch',str(args.scratch.resolve())])
    return steps


def main_job():
    p=parser()
    p.add_argument('--replace', choices=['TRUE','FALSE'], default='TRUE', help='TRUE recalculates in explicit --out; FALSE resumes/reuses')
    p.add_argument('--reuse-from', type=Path, help='Completed F-drive cache; new reports stay in --out')
    args=p.parse_args()
    if args.replace == 'TRUE' and args.out is None:
        p.error('Fresh calculation requires explicit --out NEW_DIR')
    if args.replace == 'TRUE' and args.reuse_from:
        p.error('--reuse-from requires --replace FALSE')
    if args.replace == 'TRUE' and not args.dry_run:
        from fresh_run import output_run
        protected = [os.environ.get('HERV_REUSE_DIR','/mnt/f/gen/1kg/38/allele3/herv/HG00096'),args.fasta,args.vcf_dir,args.pfile_dir,args.library]
        with output_run(args.out, True, protected):
            return run_job(args, p)
    return run_job(args, p)


def run_job(args, p):
    if args.halo < 0:
        p.error('--halo must be nonnegative')
    args.out=(args.out or Path('/mnt/d/analysis/allele3/herv')/args.gen/args.sample).resolve()
    if args.reuse_from and args.dry_run:
        print(f'[HERV] planned completed-cache reuse: {args.reuse_from}; new reports: {args.out}; no files written')
        return
    if args.reuse_from and (not (args.out/'prepared.json').exists() or (args.out/'reuse_source.json').exists()):
        from reuse_results import seed_herv
        seed_herv(args.reuse_from, args.out)
    steps=commands(args)
    cpus=sorted(os.sched_getaffinity(0))[:args.cores]
    if args.dry_run:
        print(f'[HERV] sample={args.sample}; input={args.effective_input}; chr1-22; '
              f'reference/hap1/hap2; jobs={args.jobs}; CPUs={cpus}')
        print(f'[HERV] output={args.out}; completed chunks are verified and reused')
        for command in steps:
            print(shlex.join(command))
        return
    require_library(args.library)
    legacy=Path('/mnt/d/analysis/metag/herv')/args.gen/args.sample
    if args.out == Path('/mnt/d/analysis/allele3/herv')/args.gen/args.sample and not (args.out/'prepared.json').exists() and (legacy/'prepared.json').exists():
        p.error(f'Existing results are still at {legacy}; move the entire sample directory to {args.out}, '
                f'or use --out {legacy} to resume in place. See 0.README.md.')
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/'job.lock').open('a') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            print(f'[herv] SKIP herv: another job is already running in {args.out}',flush=True)
            return
        # Refuse an already active direct engine run before overwriting job status.
        with (args.out/'run.lock').open('a') as core_lock:
            try:
                fcntl.flock(core_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:
                print(f'[herv] SKIP herv: annotation is already running in {args.out}',flush=True)
                return
        execute(args,steps,cpus)


def execute(args,steps,cpus):
    os.sched_setaffinity(0,cpus)
    os.nice(5)
    state=dict(status='running',started=now(),pid=os.getpid(),command=steps[-1],
               commands=steps,cpu_affinity=cpus,jobs=args.jobs,input_source=args.effective_input)
    save(args.out/'job.state.json',state)
    print(f'[herv] START herv sample={args.sample}; input={args.effective_input}; chr1-22; reference/hap1/hap2; '
          f'jobs={args.jobs}; output={args.out}',flush=True)
    try:
        for command in steps:
            subprocess.run(command,check=True)
        result=json.loads((args.out/'genome_fraction.json').read_text())
        with (args.out/'variant_summary.tsv').open() as f:variants=list(csv.DictReader(f,delimiter='\t'))
        lines=[result['settings']['sample']+' / GRCh38 autosomal RepeatMasker census',
               'Sequence-based annotation; variant-conditioned haplotypes, not complete independent assemblies.',
               'Unresolved spans were masked with N. Unknown inserted sequence lengths remain unknown.',
               '', 'label\tERV_with_MaLR_percent_model_length\tERV_with_MaLR_percent_non_N']
        for row in result['categories']:
            if row['category']=='ERV_with_MaLR':
                lines.append(f"{row['label']}\t{row['percent_genome']:.6f}\t{row['percent_non_N']:.6f}")
        lines.extend(['',f"Calendar hours from first annotation start (includes pauses): {result['elapsed_wall_seconds']/3600:.3f}",
            f"Measured CPU core-hours: {sum(result['measured_cpu_seconds'].values())/3600:.3f}",
            'Uncertain reference bases across both haplotypes: '+str(sum(int(r['uncertain_reference_bp']) for r in variants)),
            'Unhandled/audited allele records across both haplotypes: '+str(sum(int(r['unhandled_records']) for r in variants)),
            'Boundary disagreement (bp; ERV_with_MaLR): '+json.dumps(result['boundary_disagreement_bp']),
            '', 'Linear projection for 3202 individuals; one shared reference control:',
            json.dumps(result['cohort_projection'],indent=2),
            '', 'All results: genome_fraction.tsv/json; chromosome_fraction.tsv; variant_summary.tsv;',
            'boundary_audit.tsv; variants/<chromosome>/unhandled.tsv.gz; chunks/*/*/*/repeatmasker.out.gz.',
            'No complete individual FASTA was saved. Temporary sequence chunks are removed after each job.'])
        (args.out/'report.txt').write_text('\n'.join(lines)+'\n')
        subprocess.run([sys.executable,str(Path(__file__)), 'plot','--run',str(args.out),'--out',str(args.out/'plots')],check=True)
        state.update(status='complete',completed=now(),report=str(args.out/'report.txt'))
        save(args.out/'job.state.json',state)
        reuse = args.out/'reuse_source.json'
        if reuse.is_file():
            provenance = json.loads(reuse.read_text())
            provenance.update(status='verified_reused', verified=now())
            save(reuse, provenance)
        print(f'[herv] DONE herv sample={args.sample}; results={args.out}',flush=True)
    except BaseException as error:
        state.update(status='failed',failed=now(),error=repr(error))
        save(args.out/'job.state.json',state)
        print(f'[herv] FAIL herv: {error}',file=sys.stderr,flush=True)
        raise


# ========================================================================
# herv plot
# ========================================================================

# Publication plots for a completed 2.herv.sh individual run.

import argparse, csv, gzip, json, math
from collections import defaultdict
from pathlib import Path

def merge(rs):
    out=[]
    for a,b in sorted(rs):
        if not out or a>out[-1][1]:out.append([a,b])
        else:out[-1][1]=max(out[-1][1],b)
    return out

def collect(root,category="ERV_with_MaLR",labels=("hap1","hap2")):
    rows=[]
    for hap in labels:
        by=defaultdict(list)
        for p in Path(root).glob(f"chunks/chr*/{hap}/*/unions.bed.gz"):
            with gzip.open(p,"rt") as f:
                for line in f:
                    a=line.rstrip("\n").split("\t")
                    if len(a)>=4 and a[3]==category:by[a[0]].append((int(a[1]),int(a[2])))
        for c,rs in by.items():
            for a,b in merge(rs):rows.append((hap,c,a,b,b-a))
    return rows

def q(x,p):
    x=sorted(x)
    if not x:return None
    z=(len(x)-1)*p;i=int(z);f=z-i
    return x[i]*(1-f)+x[min(i+1,len(x)-1)]*f

def main_plot(argv=None):
    p=argparse.ArgumentParser();p.add_argument("--run",required=True);p.add_argument("--out");a=p.parse_args(argv)
    out=Path(a.out or Path(a.run)/"plots");out.mkdir(parents=True,exist_ok=True)
    rows=collect(a.run); refrows=collect(a.run,labels=("reference",)); vals=[r[-1] for r in rows if r[-1]>0]
    with gzip.open(out/"herv_union_segments.tsv.gz","wt") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(["haplotype","chrom","start","end","length_bp","coordinate_system"]);w.writerows([list(r)+["personalized_haplotype"] for r in rows])
    with gzip.open(out/"reference_herv_segments.GRCh38.tsv.gz","wt") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(["label","chrom","start","end","length_bp","coordinate_system"]);w.writerows([list(r)+["GRCh38_reference"] for r in refrows])
    with open(out/"herv_length_summary.tsv","wt") as f:
        w=csv.writer(f,delimiter="\t");w.writerow(["n","median","q25","q75","q95","max"])
        w.writerow([len(vals),q(vals,.5),q(vals,.25),q(vals,.75),q(vals,.95),max(vals) if vals else "NA"])
    if vals:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        bins=np.logspace(math.floor(math.log10(min(vals))),max(math.floor(math.log10(min(vals)))+1,math.ceil(math.log10(max(vals)))),60)
        fig,ax=plt.subplots(figsize=(8,5));ax.hist(vals,bins=bins,density=True)
        ax.set_xscale("log");ax.set(xlabel="Merged HERV/ERV segment length (bp, log scale)",ylabel="Density",title="HERV/ERV segment-length density")
        fig.tight_layout();fig.savefig(out/"herv_length_density.png",dpi=220,bbox_inches="tight");fig.savefig(out/"herv_length_density.pdf",bbox_inches="tight");plt.close(fig)
    (out/"manifest.json").write_text(json.dumps({"personalized_segments":len(rows),"reference_segments":len(refrows),"definition":"merged ERV_with_MaLR union segments; hap1/hap2 are personalized coordinates, reference is GRCh38"},indent=2)+"\n")


# ========================================================================
# herv cohort
# ========================================================================

# Reference HERV atlas and cohort genotype burden without per-person RepeatMasker.

import argparse, csv, gzip, json, random, os, sys, io
from pathlib import Path
from collections import defaultdict, Counter
from allele_common import fingerprint


def atlas(joined,out):
    from scan import prepare_joined_cache
    from scan import parse_joined_rmsk
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    cache=out/'joined_cache';prepare_joined_cache(joined,cache)
    dest=out/'herv_atlas.GRCh38.bed.gz';n=bp=0;families=Counter()
    with gzip.open(dest,'wt') as f:
        for c in [str(i) for i in range(1,23)]:
            source=cache/f'chr{c}.tsv.gz'
            if not source.exists():continue
            for e in parse_joined_rmsk(source,chromosome=c):
                if e.rep_class!='LTR':continue
                for start,end in e.segments:
                    f.write(f'chr{c}\t{start}\t{end}\t{e.tier}\t{e.family or "NA"}\t{e.strand}\t{e.start}\t{e.end}\t{e.name}\t{e.source_id}\n')
                    n+=1;bp+=end-start
                families[e.family]+=1
    result=dict(status='complete',output=str(dest),aligned_segments=n,raw_aligned_bp=bp,families=dict(families),coordinate='GRCh38 reference BED0 half-open; aligned blocks only; nested gaps and unaligned tails excluded',source=fingerprint(joined),columns=['chrom','start','end','tier','family','strand','element_start','element_end','name','source_id'])
    (out/'atlas_manifest.json').write_text(json.dumps(result,indent=2)+'\n');return result


def subset(panel,n,seed,out):
    from scan import read_panel
    groups=defaultdict(list)
    for s,g in read_panel(panel).items():groups[g].append(s)
    rng=random.Random(seed)
    for x in groups.values():rng.shuffle(x)
    chosen=[]
    while len(chosen)<n and any(groups.values()):
        for g in sorted(groups):
            if groups[g] and len(chosen)<n:chosen.append(groups[g].pop())
    Path(out).parent.mkdir(parents=True,exist_ok=True);Path(out).write_text('\n'.join(chosen)+'\n')
    return dict(n=len(chosen),output=str(out),method='round-robin genetic ancestry groups; seeded within-group shuffle')


def read_atlas(path,chrom):
    from scan import ERV
    elements={}
    with gzip.open(path,'rt') if str(path).endswith('.gz') else open(path) as f:
        for line in f:
            r=line.rstrip().split('\t')
            if r[0].removeprefix('chr')!=chrom:continue
            if len(r)<10:raise ValueError('Use herv-cohort atlas output, including aligned-block and element-edge columns')
            key=tuple(r[4:]);elements.setdefault(key,[]).append((int(r[1]),int(r[2])))
    return [ERV(chrom,int(k[2]),int(k[3]),k[4],'LTR',k[0],k[1],k[5],tuple(segments)) for k,segments in elements.items()]


def text_pvar(path):
    if str(path).endswith('.zst'):
        import zstandard
        return io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(open(path,'rb')))
    return open(path)


def burden(args):
    from scan import ERVOverlapSweep
    from scan import BoundarySweep
    import numpy as np
    import pgenlib
    from evolve import read_psam_ids
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    prefixes=([Path(args.pgen_prefix)] if args.pgen_prefix else [Path(args.pgen_dir)/f'chr{i}' for i in range(1,23)])
    merged=defaultdict(lambda:np.zeros(5,dtype=np.int64));inputs=[];expected_samples=None
    for prefix in prefixes:
        pvar=Path(str(prefix)+'.pvar.zst')
        if not pvar.exists():pvar=Path(str(prefix)+'.pvar')
        pgen=Path(str(prefix)+'.pgen');psam=Path(str(prefix)+'.psam')
        samples=read_psam_ids(psam)
        if expected_samples is None:expected_samples=set(samples)
        elif set(samples)!=expected_samples:raise ValueError('Cohort sample set differs between chromosome PGEN files')
        identity={k:fingerprint(p) for k,p in [('pgen',pgen),('pvar',pvar),('psam',psam),('atlas',args.atlas)]};inputs.append(identity)
        done=out/(prefix.name+'.burden.json');table=out/(prefix.name+'.burden.tsv')
        if done.exists() and table.exists() and (receipt:=json.loads(done.read_text())).get('inputs')==identity and receipt.get('table')==fingerprint(table):
            print(f'[herv-cohort] reuse verified {prefix.name}',flush=True)
            with table.open() as f:
                for r in csv.DictReader(f,delimiter='\t'): merged[r['sample_id'],r['label']]+=np.array([int(r[x]) for x in ('n_alt','n_indel_alt','n_nonmissing','n_assayed','n_missing')])
            continue
        print(f'[herv-cohort] genotype burden {prefix}',flush=True)
        stats={};current=None;ordinal=0
        alleles=np.empty(2*len(samples),dtype=np.int32)
        with pgenlib.PvarReader(os.fsencode(pvar),omit_chrom=True,omit_pos=True) as variants, pgenlib.PgenReader(os.fsencode(pgen),raw_sample_ct=len(samples),pvar=variants) as reader, text_pvar(pvar) as f:
            for line in f:
                if line.startswith('##'):continue
                if line.startswith('#CHROM'):header=line.lstrip('#').split();continue
                r=dict(zip(header,line.split()));c=r['CHROM'].removeprefix('chr');pos=int(r['POS']);index=ordinal;ordinal+=1
                if ordinal % 1000000 == 0:print(f'[herv-cohort] {prefix.name}: {ordinal:,} variants inspected',flush=True)
                if c not in {str(i) for i in range(1,23)} or r.get('FILTER','PASS') not in {'PASS','.'}:continue
                if c!=current:
                    elements=read_atlas(args.atlas,c);body=ERVOverlapSweep(elements);boundary=BoundarySweep(elements,('erv_core','erv_related','ltr_other'));current=c
                labels=set()
                overlap=body.locate(c,pos)
                if overlap:
                    labels.add('body:all_ltr')
                    if any(e.tier=='erv_core' for e in overlap):labels.add('body:erv_core')
                    if any(e.tier in {'erv_core','erv_related'} for e in overlap):labels.add('body:erv_related')
                for end in ('5p','3p'):
                    near=boundary.nearest(c,pos,end)
                    if near:
                        d=near['distance'];b='0-100' if d<100 else '100-500' if d<500 else '500-2000' if d<2000 else '2000-10000'
                        labels.add(end+':'+b)
                if any(label.startswith(('5p:','3p:')) for label in labels):labels.add('boundary:any')
                if not labels:continue
                bases=[r['REF']]+r['ALT'].split(',')
                if any(any(b not in 'ACGT' for b in a.upper()) for a in bases):continue
                reader.read_alleles(index,alleles)
                gt=alleles.reshape(-1,2);called=(gt>=0).all(axis=1)
                alt=((gt>0)&called[:,None]).sum(axis=1)
                indel_index=np.array([len(a)!=len(bases[0]) for a in bases])
                indels=(indel_index[np.maximum(gt,0)]&called[:,None]).sum(axis=1)
                for label in labels:
                    z=stats.setdefault(label,np.zeros((len(samples),5),dtype=np.int64))
                    z[:,0]+=alt;z[:,1]+=indels;z[:,2]+=called;z[:,3]+=1;z[:,4]+=~called
            if ordinal!=reader.get_variant_ct():raise ValueError('PVAR/PGEN variant count mismatch')
        with table.open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','label','n_alt','n_indel_alt','n_nonmissing','n_assayed','n_missing'])
            for label,z in sorted(stats.items()):
                for sample,counts in zip(samples,z):w.writerow([sample,label,*counts]);merged[sample,label]+=counts
        done.write_text(json.dumps(dict(status='complete',inputs=identity,table=fingerprint(table)),indent=2)+'\n')
    dest=out/'herv_cohort_burden_by_label.tsv'
    with dest.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','label','n_alt','n_indel_alt','n_nonmissing','n_assayed','n_missing','call_rate'])
        for (sample,label),z in sorted(merged.items()):w.writerow([sample,label,*z,z[2]/z[3] if z[3] else 'NA'])
    wide=out/'herv_cohort_burden.tsv'
    with wide.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','n_alt_erv','n_alt_erv_boundary','n_indel_erv','call_rate','n_nonmissing_erv','n_assayed_erv'])
        for sample in sorted({key[0] for key in merged}):
            body=merged.get((sample,'body:erv_related'),np.zeros(5,dtype=np.int64));boundary=merged.get((sample,'boundary:any'),np.zeros(5,dtype=np.int64))
            w.writerow([sample,body[0],boundary[0],body[1],body[2]/body[3] if body[3] else 'NA',body[2],body[3]])
    result=dict(status='complete',output=str(wide),by_label=str(dest),chromosome_inputs=len(prefixes),inputs=inputs,measure='diploid nonreference hardcall burden in shared reference aligned bodies and nearest observed boundaries; not insertion/deletion presence of entire HERV elements')
    (out/'burden_manifest.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main_cohort(argv=None):
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='cmd',required=True)
    q=sp.add_parser('atlas');q.add_argument('--joined-rmsk',required=True);q.add_argument('--out',required=True)
    q=sp.add_parser('subset');q.add_argument('--sample-panel',required=True);q.add_argument('--n',type=int,default=20);q.add_argument('--seed',type=int,default=20261024);q.add_argument('--out',required=True)
    q=sp.add_parser('burden');q.add_argument('--atlas',required=True);g=q.add_mutually_exclusive_group(required=True);g.add_argument('--pgen-prefix');g.add_argument('--pgen-dir');q.add_argument('--out',required=True)
    a=p.parse_args(argv);result=atlas(a.joined_rmsk,a.out) if a.cmd=='atlas' else subset(a.sample_panel,a.n,a.seed,a.out) if a.cmd=='subset' else burden(a)
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    try:
        mode = sys.argv[1] if len(sys.argv) > 1 else ''
        if mode in {'assembly', 'pgen'}:
            sys.argv.pop(1)
            result = main_assembly() if mode == 'assembly' else main_pgen()
        elif mode == 'plot':
            result = main_plot(sys.argv[2:])
        elif mode == 'cohort':
            result = main_cohort(sys.argv[2:])
        else:
            result = main_job()
        raise SystemExit(result if isinstance(result, int) else 0)
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f'herv: {exc}', file=sys.stderr)
        raise SystemExit(2)
