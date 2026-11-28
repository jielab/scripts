"""Competitive human/green-plant evidence. One panel and one ordinal per original pair."""
import itertools
import re
import shlex
import subprocess
import threading
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

from fastq_stream import fragments, new_digest, update_digest
from alignment_cache import panel_index
from stage_alignment import _mapper, _references

TARGETS = {'human', 'plant'}
RANKS = {'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom', 'domain'}
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


class Lineages:
    def __init__(self, tax):
        self.tax = tax
        self.chain = lru_cache(maxsize=131072)(lambda tid: tuple(tax.lineage(tid)))

    def common(self, ids):
        ids = set(ids)
        if not ids:
            return 0
        if len(ids) == 1:
            return next(iter(ids))
        chains = [self.chain(t) for t in ids]
        shared = set.intersection(*(set(c) for c in chains))
        return next((t for t in chains[0] if t in shared), 1)

    def ancestor(self, tid, ranks, skip=False):
        return next((t for t in self.chain(tid)[int(skip):] if self.tax.rank[t] in ranks), 0)

    def parent_rank(self, tid):
        return self.ancestor(tid, RANKS, skip=True) or 1

    def group(self, tid):
        if not tid:
            return 'unclassified'
        g = self.tax.describe(tid)['group']
        return 'shared' if g == 'unresolved' else g


def choose_panel(cfg, line, votes, tmp):
    opts = cfg.get('stage_alignment') or {}
    reps, table = _references(cfg, opts, tmp, line.tax)
    budget = int(cfg['rank_panel_gib'] * 2**30)
    present = {t: r for t, r in reps.items() if Path(r['resolved_path']).is_file()}
    def priority(t): return (-votes[t], t)
    candidates = sorted((t for t in votes if t in present and line.group(t) == 'plant'), key=priority)
    selected, omitted, used = {}, {}, 0

    def add(t, role):
        nonlocal used
        if t in selected:
            return
        if t not in present:
            omitted[t] = 'reference_missing'
            return
        n = int(present[t]['bases'])
        if len(selected) >= cfg['rank_max_refs'] or used+n > budget:
            omitted[t] = 'panel_budget'
            return
        selected[t] = dict(present[t], role=role)
        used += n

    def add_first(options, role, key=priority):
        # A large, highly ranked genome may not fit the remaining base budget.
        for t in sorted(options, key=key):
            add(t, role)
            if t in selected:
                return True
            if len(selected) >= cfg['rank_max_refs']:
                break
        return False

    for t in (9606, 9598):
        add(t, 'human_or_close_relative')
    if 9598 not in selected and 9595 in present:
        add(9595, 'human_close_relative')
    add_first((t for t in present if line.group(t) == 'bacteria'), 'background')
    # Give a voted non-target eukaryote a chance to compete for human/plant reads.
    # Admit the first preferred plant and its relative before spending this slot.
    euk_groups = {'fungi', 'other_eukaryote', 'mammal', 'other_vertebrate', 'invertebrate'}
    euk_background = [t for t in votes if t in present and votes[t] > 0
                      and t not in selected and line.group(t) in euk_groups]
    # Select a near relative together with each candidate, before admitting the next target.
    relatives = defaultdict(list)
    for t in present:
        if line.group(t) == 'plant':
            for node in line.chain(t):
                if line.tax.rank[node] in ('genus', 'family'):
                    relatives[node].append(t)
    forced = cfg.get('rank_panel_taxids') or []
    for t in forced:
        if t in present and line.group(t) != 'plant':
            raise ValueError('Configured candidate must be a green-plant species: '+str(t))
    # A genus/family-only Kraken vote can seed the panel even with zero species calls.
    broad = sorted((t for t in votes if t in relatives and votes[t]), key=priority)
    broad_reps = [min(relatives[t], key=priority) for t in broad]
    for t in dict.fromkeys([*forced, *candidates, *broad_reps]):
        add(t, 'plant_candidate')
        if t not in selected:
            continue
        for node in line.chain(t):
            if node in relatives:
                options = [other for other in relatives[node] if other != t]
                if add_first(options, 'plant_relative', key=lambda x: (
                        x not in selected, -votes[x], int(present[x]['bases']), x)):
                    break
        if euk_background and (not forced or t == forced[-1] or t not in forced):
            add_first(euk_background, 'eukaryote_background')
            euk_background = []
    if euk_background:
        add_first(euk_background, 'eukaryote_background')
    for group in ('fungi', 'other_eukaryote', 'mammal', 'plant'):
        options = [t for t in present if line.group(t) == group and t not in selected]
        if options:
            add(min(options, key=priority), 'background')
    for t in votes:
        if line.tax.rank.get(t) == 'species' and line.group(t) in TARGETS and t not in selected:
            omitted.setdefault(t, 'reference_missing' if t not in present else 'panel_budget')
    eligible = set()
    has_background = any(line.group(t) == 'bacteria' for t in selected)
    for t in selected:
        group = line.group(t)
        clade = line.ancestor(t, {'genus'} if group == 'plant' else {'family'})
        if group in TARGETS and clade and has_background and any(
            u != t and clade in line.chain(u) and (group != 'human' or u in (9598, 9595, 9597)) for u in selected
        ):
            eligible.add(t)
    return list(selected.values()), eligible, omitted, str(table)


def mate_hits(entries, contigs):
    """High-identity alternatives near each primary imply a shared lineage."""
    primary = None
    alternatives = []
    for f in entries:
        flag = int(f[1])
        if flag & 2048:
            continue
        if not flag & 256:
            if primary is not None:
                raise ValueError('Multiple primary SAM records for a mate')
            primary = f
        if not flag & 4:
            tags = {x[:2]: x[5:] for x in f[11:]}
            score = int(tags['AS'])
            ops = [(int(n), op) for n, op in CIGAR.findall(f[5])]
            aligned = sum(n for n, op in ops if op in 'MI=X')
            columns = sum(n for n, op in ops if op in 'MID=X')
            length = sum(n for n, op in ops if op in 'MISH=X')
            if aligned >= .80*length and columns and 1-int(tags['NM'])/columns >= .95:
                alternatives.append((score, f[2], not flag & 256))
    if primary is None:
        raise ValueError('Missing primary SAM mate')
    best = next((s for s, t, p in alternatives if p), None)
    possible = {ref for score, ref, _ in alternatives if best is not None and score >= best-5}
    return primary, possible


def pair_signal(records, contigs, line, kraken):
    ends = ([], [])
    for f in records:
        flag = int(f[1])
        if bool(flag & 64) == bool(flag & 128):
            raise ValueError('Invalid SAM mate flags')
        ends[int(bool(flag & 128))].append(f)
    (a, ra), (b, rb) = [mate_hits(end, contigs) for end in ends]
    if not ra or not rb:
        return None
    nearby = [contigs[r] for r in ra | rb]
    ids = {c['species_taxid'] for c in nearby}
    if line.group(kraken) not in TARGETS and not any(line.group(t) in TARGETS for t in ids):
        return None
    if a[2] != b[2] or not (int(a[1]) & int(b[1]) & 2) or int(a[8]) == 0 or int(a[8]) != -int(b[8]):
        return None
    if any(f[6] not in ('=', other[2]) or int(f[7]) != int(other[3]) for f, other in ((a, b), (b, a))):
        return None
    for seq in (a[9], b[9]):
        if seq == '*' or max(Counter(seq).values()) > .8*len(seq) or len({seq[i:i+3] for i in range(len(seq)-2)}) < 5:
            return None
    if any(c['contig_class'] == 'excluded' for c in nearby):
        return None
    common = line.common(ids)
    target_ids = [t for t in ids if line.group(t) in TARGETS]
    if target_ids and any(c['contig_class'] in ('chloroplast', 'mitochondrial') or c['repeat_annotation'] for c in nearby):
        ceilings = [33090 if line.group(t) == 'plant' else line.ancestor(
            t, {'family'}) or 1 for t in target_ids]
        return line.common([common, *ceilings]), 'organelle_or_repeat', ''
    if len(ids) > 1:
        return common, 'shared_alignment', ''
    if any(not 20 <= int(f[4]) < 255 for f in (a, b)):
        return None
    tid = common
    if line.group(tid) not in TARGETS:
        return tid, 'competing_non_target', ''
    if any(c['contig_class'] != 'nuclear' for c in nearby):
        return line.parent_rank(tid), 'untyped_contig', ''
    window = f"{a[2]}:{(min(int(a[3]), int(b[3]))-1)//10000}"
    return tid, 'nuclear', window


def map_all(cfg, sample, refs, calls, line, work, count, expected_digest, progress):
    """Stream interleaved ordinal reads; keep only informative pairs on disk."""
    path = work/'panel_hits.tsv'
    seeds = defaultdict(set)
    options = cfg.get('stage_alignment') or {}
    mapper = _mapper(cfg, options)
    version = subprocess.check_output([mapper, '--version'], text=True).strip()
    threads = cfg['classify_threads']
    if not refs:
        raise RuntimeError('No usable human/plant or background references; check the reference catalog')
    progress.stage('panel_lookup', references=len(refs))
    with panel_index(cfg, refs, mapper, version, threads, int(cfg['rank_panel_gib']*2**30),
                     work/'panel.fa', progress=progress) as (index, contigs, bases, cache):
        # Pin the immutable index inode. Other samples may start using the same
        # panel as soon as it is published, including while this mapper runs.
        panel_handle = Path(index).open('rb')
    with panel_handle:
        progress.stage('mapper_setup', index=index, cached=cache['index_cached'],
                       hit=cache['index_hit'], reason=cache.get('index_fallback_reason', 'none'),
                       bases=bases)
        fd = panel_handle.fileno()
        for contig in contigs.values():
            desc = contig['original_header'].lower()
            if re.search(r'\b(?:chrebv|decoy|vector|synthetic construct|contaminant)\b', desc):
                contig['contig_class'] = 'excluded'
            elif re.search(r'\b(?:chr(?:omosome)?\s*(?:m|mt|pt|pltd)|mt|pltd)\b', desc):
                contig['contig_class'] = 'mitochondrial' if 'pt' not in desc and 'pltd' not in desc else 'chloroplast'
            elif re.search(r'\b(?:unplaced|unlocalized|unlocalised|unanchored)\b', desc):
                contig['contig_class'] = 'unknown'
            elif contig['contig_class'] == 'unknown' and re.match(r'^chr(?:[1-9]\d*|x|y)\b', desc):
                contig['contig_class'] = 'nuclear'
        command = [mapper, '-ax', 'sr', '-v', '2', '--frag=yes', '--secondary=yes', '-N', '20', '-I', str(bases+1000000),
                   '-t', str(threads), f'/dev/fd/{fd}', '-']
        errors = []
        with (work/'mapping.log').open('w') as log, path.open('w') as sparse, subprocess.Popen(
                command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log, pass_fds=(fd,)) as proc:
            print(f'minimap2 {version}\n{shlex.join(command)}', file=log, flush=True)
            progress.stage('mapper_to_first_sam', pid=proc.pid, threads=threads)
            def feed():
                try:
                    digest = new_digest()
                    with proc.stdin as sink:
                        for i, (a, b) in enumerate(fragments(sample['r1'], sample['r2']), 1):
                            update_digest(digest, a, b)
                            for m, r in enumerate((a, b), 1):
                                sink.write(f'@{i}/{m}\n'.encode()+r[1].rstrip()+b'\n+\n'+r[3].rstrip()+b'\n')
                    if digest.hexdigest() != expected_digest:
                        raise RuntimeError('Raw reads changed between classification and mapping')
                except BaseException as e:
                    errors.append(e)
            feed_thread = threading.Thread(target=feed, daemon=True)
            feed_thread.start()
            try:
                def sam_rows():
                    # The first byte observed on stdout may arrive after stdio buffering.
                    first = True
                    for raw in proc.stdout:
                        if first:
                            progress.stage('sam_to_first_pair', pid=proc.pid)
                            first = False
                        if not raw.startswith(b'@'):
                            yield raw.decode().rstrip('\r\n').split('\t')
                groups = itertools.groupby(sam_rows(), key=lambda f: int(f[0].split('/')[0]))
                for i, ((kraken, mate_review), pair) in enumerate(zip(calls, groups, strict=True), 1):
                    ordinal, records = pair
                    if ordinal != i:
                        raise ValueError('Mapper ordinal/order mismatch')
                    if i == 1:
                        progress.stage('mapping_pairs', pid=proc.pid, expected=count)
                    if i % 100000 == 0:
                        progress.counted(i)
                    hit = pair_signal(records, contigs, line, kraken)
                    if hit:
                        t, kind, window = hit
                        sparse.write(f'{i}\t{t}\t{kind}\n')
                        if kind == 'nuclear' and kraken == mate_review == t and len(seeds[t]) < cfg['rank_min_windows']:
                            seeds[t].add(window)
                if proc.wait():
                    raise RuntimeError('minimap2 failed: '+log.name)
                progress.counted(count)
                print(f'Completed {count} pairs', file=log, flush=True)
            finally:
                if proc.poll() is None:
                    proc.kill()
                proc.stdout.close()
                feed_thread.join()
            if errors:
                raise errors[0]
    return path, {t: len(s) for t, s in seeds.items()}, dict(version=version, command=command, index_source=str(index), pairs=count, cache=cache)


def reconcile(kraken, signal, line, eligible, seeds, minimum):
    if signal:
        mapped, method = signal
        if method == 'nuclear' and (mapped not in eligible or seeds.get(mapped, 0) < minimum):
            mapped = line.parent_rank(mapped)
            method = 'limited_comparison_or_breadth'
        if method == 'competing_non_target':
            return line.common([kraken, mapped]), 'competing_non_target_lca'
        if not kraken:
            if method in ('organelle_or_repeat', 'untyped_contig', 'limited_comparison_or_breadth'):
                return 0, 'weak_alignment_only'
            return mapped, 'alignment_'+method
        if mapped in line.chain(kraken):
            return mapped, method
        if kraken in line.chain(mapped):
            return mapped, 'refined_'+method
        return line.common([kraken, mapped]), 'conflict_lca'
    if line.group(kraken) in TARGETS and line.tax.rank.get(kraken) == 'species':
        return line.parent_rank(kraken), 'unconfirmed_species'
    return kraken, 'classifier_only' if kraken else 'unclassified'
