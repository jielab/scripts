"""Local, paired evidence audit for eukaryotic source bins. No FASTQ is republished."""
from collections import Counter, defaultdict
import hashlib
from pathlib import Path
import random
import subprocess
from fastq_stream import fragments, new_digest, update_digest

LIMIT = 100

def resolve_mate(hits, taxonomy):
    if len(hits) >= LIMIT:
        return None,'hit_limit'
    keep = [(ids,score) for s,ids,score,identity,coverage in hits if identity >= 95 and coverage >= 80]
    if not keep:
        return None,'insufficient'
    best = max(score for _,score in keep)
    taxa = set()
    for ids,score in keep:
        if score < 0.98*best:
            continue
        for item in ids.split(';'):
            if not item.isdecimal() or int(item) == 0:
                return None,'missing_taxonomy'
            tid = int(item)
            if tid not in taxonomy.parent and tid not in taxonomy.merged:
                return None,'missing_taxonomy'
            taxa.add(taxonomy.resolve(tid))
    chains = [list(taxonomy.lineage(t)) for t in taxa]
    shared = set.intersection(*(set(c) for c in chains))
    return next(t for t in chains[0] if t in shared),'ok'

def paired(left, right, original, taxonomy):
    calls = [resolve_mate(end,taxonomy) for end in (left,right)]
    for label in ('hit_limit','missing_taxonomy','insufficient'):
        if any(state == label for _,state in calls):
            return label,0
    a,b = (t for t,_ in calls)
    left_chain = list(taxonomy.lineage(a))
    right_chain = set(taxonomy.lineage(b))
    lca = next(t for t in left_chain if t in right_chain)
    lineage = set(taxonomy.lineage(lca))
    if original in lineage:
        return 'supports_bin_rank',lca
    if lca in taxonomy.lineage(original):
        return 'higher_rank_only',lca
    return 'different_lineage',lca

def scope(tid, group, rank, panel, taxonomy):
    if panel is None:
        return 'external_inventory_unverified'
    if not any(tid in lineage for lineage in panel.values()):
        return 'target_absent'
    if rank == 'genus':
        close = tid
        relative = sum(close in c for c in panel.values()) > 1
    else:
        level = 'genus' if group == 'plant' else 'family'
        close = next((t for t in taxonomy.lineage(tid) if taxonomy.rank[t] == level), None)
        relative = close and any(close in c and tid not in c for c in panel.values())
    if not relative:
        return 'close_relative_absent'
    if all(close in lineage for lineage in panel.values()):
        return 'background_absent'
    return 'target_relative_background_present'

def audit(sample, folder, src, taxa, db, taxonomy, cap, threads, tmp, reference_ids=None):
    query, result = tmp/'audit.fa', tmp/'audit.tsv'
    # Use a deterministic uniform sample in each bin, independent of FASTQ headers.
    manifest = {r["taxid"]:r for r in taxa}
    panel = ({t:set(taxonomy.lineage(t)) for t in reference_ids} if reference_ids is not None else None)
    chosen = []
    excluded = []
    with query.open('wb') as out:
        for r in src:
            _,group,rank,tid,name,n,food,_ = r
            if rank != 'species' and not (rank == 'genus' and food != 'unlisted'):
                continue
            coverage = scope(int(tid),group,rank,panel,taxonomy)
            if coverage not in ('external_inventory_unverified','target_relative_background_present'):
                excluded.append([sample,tid,name,rank,n,0,coverage,'panel_not_evaluable','NA','NA','NA',0])
                continue
            declared = manifest[tid]
            rel = Path(declared['bin'])
            if rel.is_absolute() or '..' in rel.parts:
                raise ValueError('Invalid bin manifest path: '+str(rel))
            base = folder/rel
            digest = new_digest()
            seed = hashlib.sha256((sample+'\t'+tid).encode()).digest()
            targets = set(random.Random(seed).sample(range(n),min(cap,n)))
            seen = 0
            for pos,(a,b) in enumerate(fragments(base/'R1.fastq.gz',base/'R2.fastq.gz')):
                seen += 1
                update_digest(digest,a,b)
                if pos in targets:
                    number = len(chosen)
                    chosen.append((int(tid), name, rank, n, coverage))
                    for end,record in enumerate((a,b),1):
                        out.write(f'>{number}_{end}\n'.encode()+record[1].strip()+b'\n')
            if seen != n or digest.hexdigest() != declared['pairs_sha256']:
                raise ValueError(f'Bin count or FASTQ checksum differs for {sample}, taxid {tid}')
    hits = defaultdict(list)
    if chosen:
        subprocess.run(['blastn','-task','blastn','-query',str(query),'-db',db,
                        '-num_threads',str(threads),'-dust','yes','-max_target_seqs',str(LIMIT),
                        '-max_hsps','1','-outfmt','6 qseqid sseqid staxids bitscore pident qcovhsp',
                        '-out',str(result)],check=True)
        with result.open() as f:
            for line in f:
                qry,subject,ids,score,pid,cov = line.rstrip('\n').split('\t')
                num,end = map(int,qry.split('_'))
                if not (0 <= num < len(chosen) and end in (1,2)):
                    raise ValueError('Unexpected BLAST query label')
                hits[num,end].append((subject,ids,float(score),float(pid),float(cov)))
    counts = Counter()
    for i,bin_info in enumerate(chosen):
        outcome,tid = paired(hits[i,1],hits[i,2],bin_info[0],taxonomy)
        counts[bin_info,outcome,tid] += 1
    return excluded + [[sample,orig,name,rank,total,min(cap,total),coverage,outcome,tid or 'NA',
             taxonomy.rank.get(tid,'NA'),taxonomy.name.get(tid,'NA'),n]
            for ((orig,name,rank,total,coverage),outcome,tid),n in sorted(counts.items())]
