#!/usr/bin/env python3
"""Chunked, sequence-based ERV census of a VCF-personalized diploid reference.

No individual chromosome or whole-genome FASTA is saved. Unresolvable alleles
are audited and masked with N in reference-coordinate spans (plus a declared
flank); their unknown sequence lengths cannot be inferred from this model.
"""
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


def run_chunk(args, descriptor, model, signature):
    chrom, label, a, b, length = descriptor
    folder = args.out/'chunks'/chrom/label/f'{a:09d}-{b:09d}'
    receipt = folder/'done.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old['signature'] != signature or any(digest(folder/n) != h for n, h in old['files'].items()):
            raise ValueError(f'Changed chunk receipt: {folder}')
        return old
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


def run(args):
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
    if manifest.exists() and json.loads(manifest.read_text())['settings'] != settings:
        raise ValueError('Run settings changed; use a new output directory')
    if not manifest.exists():
        save(manifest, dict(settings=settings, started=now(), started_unix=time.time(), jobs=args.jobs,
                           cpu_affinity=sorted(os.sched_getaffinity(0)),
                           max_search_threads=4*args.jobs, nice=os.getpriority(os.PRIO_PROCESS,0),
                           scope=prepared['chromosomes'], status='running'))
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
    save(args.out/'plan.json', dict(chunks=len(jobs), reference_autosomal_bp=sum(lengths[c] for c in prepared['chromosomes']),
         descriptors=[j[0] for j in jobs]))
    with ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_chunk, args, *job) for job in jobs]
        completed = 0
        try:
            for f in as_completed(futures):
                f.result()
                completed += 1
                save(args.out/'progress.json', dict(completed_chunks=completed, total_chunks=len(jobs), updated=now()))
        except BaseException:
            for f in futures:
                f.cancel()
            raise
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


def main():
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


if __name__ == '__main__':
    main()
