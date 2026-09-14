#!/usr/bin/env python3
"""Read PRS-CSx 1KG LD by SNP ID; relocate coordinates for GRCh38 display."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import h5py

RACES = ('EUR', 'AFR', 'EAS', 'SAS', 'AMR')


def identity(path):
    path = Path(path)
    if not path.is_file():
        return [str(path.resolve()), None]
    stat = path.stat()
    return [str(path.resolve()), stat.st_size, stat.st_mtime_ns]


def key(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:24]


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(value, allow_nan=False, separators=(',', ':')))
    tmp.replace(path)


def positions(root, cache, build, chain_dir, liftover):
    """snpinfo_mult_1kg_hm3 BP is GRCh37; lifted positions retain source chr."""
    source = root / 'snpinfo_mult_1kg_hm3'
    chain = chain_dir / 'hg19ToHg38.over.chain.gz'
    signature = [identity(source), build, identity(chain) if build == '38' else None,
                 identity(liftover) if build == '38' else None, identity(__file__)]
    path = cache / ('positions-' + key(signature) + '.npz')
    if path.is_file():
        with np.load(path, allow_pickle=False) as data:
            return {k: data[k] for k in data.files}
    import pandas as pd
    data = pd.read_csv(source, sep=r'\s+', usecols=['CHR', 'SNP', 'BP'],
                       dtype={'CHR': np.int16, 'SNP': str, 'BP': np.int64})
    source_chr = data.CHR.to_numpy()
    chr_ = source_chr.copy()
    pos = data.BP.to_numpy().copy()
    snps = data.SNP.to_numpy(dtype=str)
    if build == '38':
        if not chain.is_file() or not liftover.is_file():
            raise FileNotFoundError('GRCh38 LD display requires hg19ToHg38 chain and liftOver')
        with tempfile.TemporaryDirectory(prefix='ld-map-') as work:
            bed, mapped, unmapped = [Path(work) / name for name in ('input.bed', 'mapped.bed', 'unmapped.bed')]
            with bed.open('w') as handle:
                for i, (ch, bp) in enumerate(zip(source_chr, pos)):
                    handle.write(f'chr{ch}\t{bp-1}\t{bp}\t{i}\n')
            subprocess.run([str(liftover), str(bed), str(chain), str(mapped), str(unmapped)], check=True)
            chr_.fill(0)
            pos.fill(0)
            seen = set()
            with mapped.open() as handle:
                for line in handle:
                    ch, start, end, i = line.split()[:4]
                    i = int(i)
                    if i in seen:
                        chr_[i] = 0
                        continue
                    seen.add(i)
                    ch = ch.removeprefix('chr')
                    if ch.isdigit() and 1 <= int(ch) <= 22 and int(end) - int(start) == 1:
                        chr_[i], pos[i] = int(ch), int(end)
    result = dict(snp=snps, chr=chr_, pos=pos, source_chr=source_chr)
    stage = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with stage.open('wb') as handle:
        np.savez_compressed(handle, **result)
    stage.replace(path)
    return result


def index(path, cache):
    cache_file = cache / ('index-' + key([identity(path), 1]) + '.json')
    if cache_file.is_file():
        return json.loads(cache_file.read_text())
    result = {}
    with h5py.File(path, 'r') as handle:
        for group in handle:
            if 'snplist' not in handle[group] or 'ldblk' not in handle[group]:
                continue
            snps = handle[group]['snplist'].asstr()[:]
            if handle[group]['ldblk'].shape != (len(snps), len(snps)):
                raise ValueError(f'{path}:{group}: LD matrix and snplist differ in size')
            for row, snp in enumerate(snps):
                if snp in result:
                    raise ValueError(f'{path}: duplicate reference SNP {snp}')
                result[snp] = [group, row]
    write_json(cache_file, result)
    return result


def matrix(axis, mappings):
    """Unknown SNPs and cross-HDF5-block pairs are NA, never zero."""
    result = np.full((len(axis), len(axis)), np.nan)
    groups = {}
    for col, snp in enumerate(axis):
        if snp not in mappings:
            continue
        path, group, row = mappings[snp]
        groups.setdefault((path, group), []).append((int(row), col))
    for (path, group), rows in groups.items():
        rows.sort()
        source, dest = zip(*rows)
        with h5py.File(path, 'r') as handle:
            ld = handle[group]['ldblk'][list(source), :][:, list(source)]
        diag = np.diag(ld)
        scale = np.sqrt(np.outer(diag, diag))
        good = np.isfinite(scale) & (scale > 0) & (diag[:, None] > 0) & (diag[None, :] > 0)
        corr = np.divide(ld, scale, out=np.full_like(ld, np.nan), where=good)
        corr = (corr + corr.T) / 2
        r2 = np.clip(corr, -1, 1) ** 2
        result[np.ix_(dest, dest)] = r2
    return result, len(groups)


def query(a):
    root, cache = Path(a.root), Path(a.cache)
    cache.mkdir(parents=True, exist_ok=True)
    if str(a.chr) not in map(str, range(1, 23)):
        return dict(build=a.grch, chr=a.chr, start=a.start, end=a.end, snps=[], pos=[], total=0,
                    populations=[dict(race=r, note='1KG LD reference has autosomes 1–22 only', available=0, shown=0) for r in RACES])
    data = positions(root, cache, a.grch, Path(a.chain_dir), Path(a.liftover_bin))
    indices = np.flatnonzero((data['chr'] == int(a.chr)) & (data['pos'] >= a.start) & (data['pos'] <= a.end))
    indices = indices[np.argsort(data['pos'][indices], kind='stable')]
    candidates = set(data['snp'][indices])
    chroms = sorted(set(data['source_chr'][indices].tolist()))
    files = {race: [root / f'ldblk_1kg_{race}' / f'ldblk_1kg_chr{ch}.hdf5' for ch in chroms] for race in RACES}
    signature = [a.grch, a.chr, a.start, a.end, a.max_snps, identity(__file__),
                 identity(root / 'snpinfo_mult_1kg_hm3'), identity(Path(a.chain_dir) / 'hg19ToHg38.over.chain.gz'),
                 [[identity(p) for p in files[r]] for r in RACES]]
    saved = cache / ('region-' + key(signature) + '.json')
    if saved.is_file():
        return json.loads(saved.read_text())
    maps, notes = {}, {}
    for race in RACES:
        maps[race] = {}
        missing = []
        try:
            for path in files[race]:
                if not path.is_file():
                    missing.append(path.name)
                    continue
                rows = index(path, cache)
                for snp in candidates.intersection(rows):
                    group, row = rows[snp]
                    maps[race][snp] = (str(path), group, row)
            if missing:
                notes[race] = 'Missing: ' + ', '.join(missing)
        except (OSError, ValueError, KeyError) as error:
            maps[race] = {}
            notes[race] = str(error)
    present = set().union(*(set(v) for v in maps.values()))
    indices = np.array([i for i in indices if data['snp'][i] in present], dtype=int)
    total = len(indices)
    if total > a.max_snps:
        indices = indices[np.linspace(0, total - 1, a.max_snps).round().astype(int)]
    axis = data['snp'][indices].tolist()
    populations = []
    for race in RACES:
        item = dict(race=race, available=len(maps[race]), shown=sum(s in maps[race] for s in axis),
                    note=notes.get(race, ''))
        if item['shown']:
            try:
                r2, groups = matrix(axis, maps[race])
                item.update(groups=groups, r2=[[None if not np.isfinite(v) else round(float(v), 4) for v in row] for row in r2])
            except (OSError, ValueError, KeyError) as error:
                item['note'] = str(error)
        elif not item['note']:
            item['note'] = 'No reference SNPs in this region'
        populations.append(item)
    result = dict(build=a.grch, chr=a.chr, start=a.start, end=a.end, total=total,
                  snps=axis, pos=data['pos'][indices].tolist(), populations=populations)
    write_json(saved, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'cache', 'chr', 'output', 'chain-dir', 'liftover-bin'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--grch', choices=('37', '38'), required=True)
    parser.add_argument('--start', type=int, required=True)
    parser.add_argument('--end', type=int, required=True)
    parser.add_argument('--max-snps', type=int, default=300)
    a = parser.parse_args()
    if not (1 <= a.start < a.end) or not 2 <= a.max_snps <= 1000:
        parser.error('Require 1 <= start < end and 2 <= max-snps <= 1000')
    write_json(a.output, query(a))


if __name__ == '__main__':
    main()
