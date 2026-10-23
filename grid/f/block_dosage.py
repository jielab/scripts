"""PGEN supplements for weighted SNPs absent from the PRSformer dosage cache."""
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd


class SupplementedDosage:
    """Read batches from two matrices without copying the full PRSformer cache."""
    def __init__(self, original, supplement):
        self.original = original
        self.supplement = supplement
        self.shape = (original.shape[0], original.shape[1] + supplement.shape[1])

    def __getitem__(self, key):
        if isinstance(key, slice):
            return np.concatenate((self.original[key], self.supplement[key]), axis=1)
        rows, cols = key  # np.ix_ used by independent reconstruction QC
        rows = np.asarray(rows).reshape(-1)
        cols = np.asarray(cols).reshape(-1)
        out = np.empty((len(rows), len(cols)), dtype=np.float32)
        old = cols < self.original.shape[1]
        out[:, old] = self.original[np.ix_(rows, cols[old])]
        out[:, ~old] = self.supplement[np.ix_(rows, cols[~old] - self.original.shape[1])]
        return out


def supplement(config, prep_meta, cohort, variants, geno, missing, sha256):
    from prsformer import genotype_files, sample_ids, scan_variants
    import pgenlib

    # CSx may use a broader target panel (e.g. typ) than PRSformer (e.g. hap).
    # Only the missing columns come from this source; retain the frozen cache.
    directory = config.get('dosage_supplement_dir') or prep_meta.get('options', {}).get('dir_gen')
    if not directory:
        raise ValueError('Set dosage_supplement_dir or provide PRSformer prepare.json options.dir_gen')
    directory = str(Path(directory).expanduser().resolve())
    print(f'[GRID ABM_GPU] Dosage supplement source: {directory}', flush=True)
    files = genotype_files(directory, sorted(missing.CHR.unique()))
    if any(mode != 'pgen' for _, _, _, _, mode in files):
        raise ValueError(f'Exact ABM_GPU dosage supplementation requires PGEN files: {directory}')
    allowed = {}
    for w in missing.itertuples():
        coord = (int(w.CHR), int(w.BP))
        if w.SNP in allowed and allowed[w.SNP] != coord:
            raise ValueError(f'Conflicting PRS-CSx coordinates for {w.SNP}')
        allowed[w.SNP] = coord
    extra = scan_variants(files, allowed, 0)
    absent = sorted(set(allowed) - set(extra.SNP))
    if absent:
        raise ValueError(f'{len(absent)} weighted SNPs unavailable in biallelic PGEN at {directory}: '
                         f'{absent[:20]}; set dosage_supplement_dir to the matching CSx target PGEN directory')
    # Fail before writing a large cache on allele or coordinate disagreements.
    aligned = missing.merge(extra, on=['SNP', 'CHR', 'BP'], validate='many_to_one')
    alt = aligned.A1 == aligned.ALT
    valid = ((alt | (aligned.A1 == aligned.REF)) &
             ((aligned.A2 == '') | (aligned.A2 == np.where(alt, aligned.REF, aligned.ALT))))
    if not valid.all():
        raise ValueError(f'PGEN allele mismatch: {aligned.loc[~valid, "SNP"].head(20).tolist()}')
    sources = []
    for _, gp, vp, sp, _ in files:
        for path in (gp, vp, sp):
            stat = path.stat()
            sources.append(dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                                ctime_ns=stat.st_ctime_ns,
                                metadata_sha256=sha256(path) if path != gp else None))
    identity = dict(version=2, source_directory=directory, sources=sources, variants=extra.to_dict('records'),
                    cohort=hashlib.sha256(cohort.eid.to_csv(index=False).encode()).hexdigest(),
                    source_sha256=sha256(Path(__file__)))
    signature = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    root = Path(config['cache_root']) / 'dosage_supplements'
    root.mkdir(parents=True, exist_ok=True)
    target = root / signature
    shape = (len(cohort), len(extra))
    if not target.exists() or ((target / 'genotypes.released.json').is_file()
                               and not (target / 'genotypes.npy').exists()):
        print(f'[GRID ABM_GPU] Supplement {len(extra):,} exact weighted SNPs from PGEN', flush=True)
        with tempfile.TemporaryDirectory(prefix='prepare-', dir=root) as td:
            work = Path(td)
            # Fortran layout allows sequential variant-major PGEN writes.
            out = np.lib.format.open_memmap(work / 'genotypes.npy', mode='w+',
                                           dtype=np.float32, shape=shape, fortran_order=True)
            for chrom, gp, vp, sp, mode in files:
                vc = extra[extra.CHR == chrom]
                ids = sample_ids(sp, mode)
                lookup = {eid: i for i, eid in enumerate(ids)}
                absent_ids = [eid for eid in cohort.eid if eid not in lookup]
                if absent_ids:
                    raise ValueError(f'chr{chrom}: {len(absent_ids)} cohort IDs missing from PGEN')
                original = np.asarray([lookup[eid] for eid in cohort.eid], dtype=np.uint32)
                subset = np.sort(original)
                reorder = np.searchsorted(subset, original)
                with pgenlib.PvarReader(os.fsencode(vp)) as pv:
                    with pgenlib.PgenReader(os.fsencode(gp), raw_sample_ct=len(ids),
                                           sample_subset=subset, pvar=pv) as pg:
                        for start in range(0, len(vc), 32):
                            chunk = vc.iloc[start:start + 32]
                            buf = np.empty((len(chunk), len(cohort)), dtype=np.float32)
                            pg.read_dosages_list(np.asarray(chunk.source_index, dtype=np.uint32), buf, allele_idx=1)
                            unknown = (buf == -9) | ~np.isfinite(buf)
                            if np.any(~unknown & ((buf < 0) | (buf > 2))):
                                raise ValueError(f'chr{chrom}: invalid PGEN dosage')
                            buf[unknown] = -1
                            out[:, chunk.index.to_numpy()] = buf[:, reorder].T
                print(f'[GRID ABM_GPU] Supplement chr{chrom}: {len(vc):,} SNPs', flush=True)
            out.flush()
            del out
            extra.to_csv(work / 'variants.tsv.gz', sep='\t', index=False)
            receipt = dict(signature=signature, identity=identity, shape=list(shape),
                           genotype_sha256=sha256(work / 'genotypes.npy'),
                           dosage='ALT float32; missing=-1; exact weighted SNPs without PRSformer MAF/call-rate filtering')
            (work / 'COMPLETE.json').write_text(json.dumps(receipt, indent=2) + '\n')
            target.mkdir(exist_ok=True)
            for name in ('genotypes.npy', 'variants.tsv.gz', 'COMPLETE.json'):
                os.replace(work / name, target / name)
            (target / 'genotypes.released.json').unlink(missing_ok=True)
    receipt = json.loads((target / 'COMPLETE.json').read_text())
    path = target / 'genotypes.npy'
    if receipt['signature'] != signature or sha256(path) != receipt['genotype_sha256']:
        raise ValueError(f'Invalid dosage supplement: {target}')
    added = np.load(path, mmap_mode='r', allow_pickle=False)
    if added.shape != shape:
        raise ValueError(f'Dosage supplement shape mismatch: {target}')
    print(f'[GRID ABM_GPU] Exact dosage supplement ready: {len(extra):,} SNPs', flush=True)
    return pd.concat([variants, extra], ignore_index=True), SupplementedDosage(geno, added), receipt
