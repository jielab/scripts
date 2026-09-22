#!/usr/bin/env python3
"""Validation, provenance and atomic-output staging for the three launchers."""
import argparse, gzip, hashlib, json, os
from pathlib import Path
import numpy as np
import pandas as pd
from prepare_sumstats import MIN_COORDINATE_MATCH,choose,reader_options
from score_output import filter_samples


def stamp(path):
    p = Path(path)
    if not p.is_file():
        raise ValueError(f'Missing file: {p}')
    z = p.stat()
    return [str(p.resolve()), z.st_size, z.st_mtime_ns]


def inspect(snpinfo, files):
    ref = pd.read_csv(snpinfo, sep=r'\s+', usecols=['SNP', 'CHR', 'BP'], dtype={'SNP': str})
    ref = ref.drop_duplicates('SNP').set_index('SNP')
    for path in files:
        marker = Path(path + '.grch')
        if marker.exists() and marker.read_text().strip() != '37':
            raise ValueError(f'GRCh37 required: {marker}')
        d = pd.read_csv(path, nrows=100000, **reader_options(path))
        columns = {key: choose(d.columns, key) for key in ('SNP','CHR','BP','A1','A2','BETA','SE','N')}
        missing = [key for key in ('A1','A2','SE') if columns[key] is None]
        if columns['BETA'] is None and choose(d.columns,'OR') is None:missing.append('BETA or OR')
        if columns['SNP'] is None and (columns['CHR'] is None or columns['BP'] is None):missing.append('rsID or CHR/BP')
        if missing:
            raise ValueError(f'{path}: missing {missing}')
        if columns['SNP'] is None or columns['CHR'] is None or columns['BP'] is None:
            print(f'input GWAS: {path}; coordinate completion uses the explicitly supplied GRCh37 SNPINFO')
            continue
        check = pd.DataFrame({'SNP': d[columns['SNP']].astype(str),
                              'CHR': pd.to_numeric(d[columns['CHR']].astype(str).str.strip().str.replace(r'^chr','',regex=True),errors='coerce'),
                              'POS': pd.to_numeric(d[columns['BP']],errors='coerce')})
        z = check.merge(ref, on='SNP', suffixes=('', '_ref'))
        ok = (z.CHR == z.CHR_ref) & (z.POS == z.BP)
        n = len(z); rate = float(ok.mean()) if n else 0
        if n < min(100,len(ref)) or rate < MIN_COORDINATE_MATCH:
            raise ValueError(f'{path}: GRCh37 reference-coordinate check failed ({int(ok.sum())}/{n})')
        print(f'input GWAS: {path}\n  GRCh37 coordinate sample: {int(ok.sum())}/{n}; N is median of usable HM3 variants (override: --n-gwas)')


def coverage(chrs, files):
    """Read the BGZF index when available; plain raw gzip is checked after preparation."""
    import struct
    wanted = set(chrs.split())
    for path in files:
        index = Path(path + '.tbi')
        if not index.is_file():
            continue
        with gzip.open(index, 'rb') as stream:
            if stream.read(4) != b'TBI\x01':
                raise ValueError(f'Invalid tabix index: {index}')
            header = stream.read(32)
            if len(header) != 32:
                raise ValueError(f'Truncated tabix index: {index}')
            length = struct.unpack('<8i', header)[7]
            names = stream.read(length).decode().strip('\x00').split('\x00')
        present = {s.removeprefix('chr') for s in names}
        missing = sorted(wanted - present, key=int)
        if missing:
            raise ValueError(f'{path}: source GWAS is missing chromosomes {",".join(missing)} (tabix index). '
                             'Repair the formatted common GWAS from the original source data first; '
                             'do not silently omit missing chromosomes.')


def signature(values, files):
    obj = {'settings': values, 'files': [stamp(p) for p in files]}
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def weights(inputs, output):
    d = pd.concat([pd.read_csv(p, sep='\t') for p in inputs], ignore_index=True)
    if d.empty or d.SNP.duplicated().any() or not np.isfinite(d.BETA).all():
        raise ValueError('Invalid/duplicate posterior weights')
    d.to_csv(output, sep='\t', index=False, compression='gzip')
    print(f'posterior variants: {len(d)}')


def disco_inputs(pca, centers, score_dir, outdir, npc, remove='/mnt/d/files/ukb.exclude.id'):
    if not 5 <= npc <= 20: raise ValueError('Disco distance PCs must be 5..20')
    pc = [f'PC{i}' for i in range(1,npc+1)]
    d = pd.read_csv(pca, sep='\t', dtype={'IID': str, 'eid': str, '#IID': str})
    idcol = next(x for x in ('IID','eid','#IID') if x in d)
    d = d.rename(columns={idcol:'IID'})[['IID']+pc]
    if d.IID.isna().any() or d.IID.duplicated().any() or not np.isfinite(d[pc].to_numpy()).all():
        raise ValueError('Invalid PCA IDs or PCs')
    med = pd.read_csv(centers, sep='\t'); med = med[[med.columns[0]]+pc]
    pops = ['AFR','EAS','EUR','SAS']
    if list(med.iloc[:,0]) != pops or not np.isfinite(med[pc].to_numpy()).all():
        raise ValueError('Centers must have finite PCs in AFR,EAS,EUR,SAS order')
    scores = []
    merged = Path(score_dir)/'csx.pgs.gz'
    combined = pd.read_csv(merged, sep='\t', dtype={'eid':str}) if merged.is_file() else None
    for pop in pops:
        p = Path(score_dir)/f'csx.{pop}.tsv.gz'
        source = combined if combined is not None else pd.read_csv(p, sep='\t', dtype={'eid':str, 'IID':str})
        z = source.rename(columns={'eid':'IID', f'CSX_{pop}':'PRS', f'csx.{pop}':'PRS'})[['IID','PRS']]
        if z.IID.isna().any() or z.IID.duplicated().any() or not np.isfinite(z.PRS).all(): raise ValueError(f'Invalid scores: {p}')
        if z.PRS.std() == 0: raise ValueError(f'Constant scores: {p}')
        if scores and set(z.IID) != set(scores[0].IID): raise ValueError('Population PRS sample sets differ')
        scores.append(z)
    before = len(scores[0])
    scores = [filter_samples(z, 'IID', remove) for z in scores]
    print(f'Disco withdrawn filter: excluded={before-len(scores[0])}; retained={len(scores[0])}')
    ids = set(scores[0].IID)
    retained = ids.intersection(d.IID)
    print(f'Disco sample filter: scored={len(ids)}; missing PCA={len(ids-retained)}; retained={len(retained)}')
    if len(retained) < 2: raise ValueError('Fewer than two scored samples remain after PCA filtering')
    d = d[d.IID.isin(retained)].sort_values('IID')
    scores = [z[z.IID.isin(retained)].sort_values('IID') for z in scores]
    for pop,z in zip(pops,scores):
        if z.PRS.std() == 0: raise ValueError(f'Constant {pop} scores after PCA filtering')
    for _, row in med.iterrows():
        dist = np.linalg.norm(d[pc].to_numpy()-row[pc].to_numpy(dtype=float), axis=1)
        if (dist <= 0).any(): raise ValueError('Sample exactly at reference center; official interpolation is undefined')
    target = Path(outdir); target.mkdir(parents=True, exist_ok=True)
    d.to_csv(target/'pca.tsv', sep='\t', index=False)
    med.to_csv(target/'centers.tsv', sep='\t', index=False)
    for pop,z in zip(pops,scores): z.to_csv(target/f'{pop}.tsv',sep='\t',index=False)
    print(f'Disco aligned samples: {len(d)}; distance PCs: {npc}')


def validate_disco(path, inputs):
    d = pd.read_csv(path, sep='\t', dtype={'IID':str})
    ref = pd.read_csv(Path(inputs)/'AFR.tsv', sep='\t', dtype={'IID':str})
    if d.empty or d.IID.duplicated().any() or set(d.IID) != set(ref.IID) or not np.isfinite(pd.to_numeric(d.PRS)).all():
        raise ValueError('Official DiscoDivas returned invalid/incomplete scores')
    co = pd.read_csv(path.replace('.tsv.gz','.coef.tsv.gz'), sep='\t', dtype={'IID':str})
    vals = co.drop(columns='IID').to_numpy()
    if co.IID.duplicated().any() or set(co.IID) != set(ref.IID) or not np.isfinite(vals).all() or not np.allclose(vals.sum(axis=1),1):
        raise ValueError('Invalid interpolation coefficients')
    print(f'validated Disco scores: {len(d)}')


if __name__ == '__main__':
    import sys
    action, *args = sys.argv[1:]
    try:
        if action == 'inspect': inspect(args[0],args[1:])
        elif action == 'coverage': coverage(args[0],args[1:])
        elif action == 'signature':
            split = args.index('--files'); print(signature(args[:split], args[split+1:]))
        elif action == 'weights': weights(args[1:],args[0])
        elif action == 'disco-inputs': disco_inputs(*args[:4],int(args[4]), *args[5:])
        elif action == 'disco-output': validate_disco(*args)
        else: raise ValueError(f'Unknown action: {action}')
    except (ValueError, KeyError, OSError) as e:
        raise SystemExit(f'ERROR: {e}')
