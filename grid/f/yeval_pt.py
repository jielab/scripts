#!/usr/bin/env python3
"""Ancestry-matched COJO scores; one sparse genotype pass per chromosome."""
import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

POPS = ('EUR', 'AFR', 'EAS', 'SAS')


def prepare_scoring_index(path, wanted):
    """Give each extracted PGEN row a unique scoring ID, preserving row order.

    Imputed data can contain split multiallelic records with the same rsID.
    Only the temporary PVAR is rewritten; the source genotype files are untouched.
    """
    found = {}
    tmp = Path(str(path) + '.tmp')
    with path.open() as handle, tmp.open('w') as output:
        columns = None
        row_number = 0
        for line in handle:
            if line.startswith('##'):
                output.write(line)
                continue
            bits = line.split()
            if not bits:
                output.write(line)
                continue
            if columns is None:
                columns = bits
                chrcol, poscol, idcol, refcol, altcol = (
                    columns.index(z) for z in ('#CHROM', 'POS', 'ID', 'REF', 'ALT'))
                output.write(line)
                continue
            sid = bits[idcol]
            row_number += 1
            score_id = f'PTv{row_number}'
            if sid in wanted:
                found.setdefault(sid, []).append(dict(
                    chromosome=int(bits[chrcol]), position=int(bits[poscol]),
                    ref=bits[refcol], alt=bits[altcol], score_id=score_id,
                    alleles=set([bits[refcol], *bits[altcol].split(',')]),
                    variant=':'.join(bits[c] for c in (chrcol, poscol, refcol, altcol))))
            bits[idcol] = score_id
            output.write('\t'.join(bits) + '\n')
    tmp.replace(path)
    return found


def match_variant(candidates, chromosome, position, allele):
    """Select exactly one COJO-compatible row; never choose an arbitrary duplicate."""
    if not candidates:
        return None, 'absent'
    at_position = [v for v in candidates if v['chromosome'] == chromosome
                   and v['position'] == position]
    if not at_position:
        return None, 'position_mismatch'
    compatible = [v for v in at_position if allele in v['alleles']]
    if not compatible:
        return None, 'allele_mismatch'
    if len(compatible) != 1:
        return None, 'ambiguous'
    return compatible[0], 'scored'


def score(a):
    dest = Path(a.output)
    sources, weights = {}, {}
    for pop in POPS:
        tag = f'{a.trait}.{pop}'
        src = Path(a.dir_gwas) / tag / 'gwas' / f'{tag}.jma.cojo'
        if not src.is_file() and a.trait == 't2dm' and pop == 'AFR':
            src = Path(a.dir_gwas) / 't2dm.AFA/gwas/t2dm.AFA.jma.cojo'
        if not src.is_file():
            raise FileNotFoundError(f'PT {pop}: {src}; supply --pt-file for an existing table')
        sources[pop] = src
        w = pd.read_csv(src, sep=r'\s+')
        cols = ['Chr', 'SNP', 'bp', 'refA', a.effect]
        if not set(cols).issubset(w):
            raise ValueError(f'Missing COJO columns {cols}: {src}')
        w = w[cols].copy()
        if w.empty or w.SNP.duplicated().any() or w.isna().any().any():
            raise ValueError(f'Empty, duplicated or missing COJO weights: {src}')
        w[a.effect] = pd.to_numeric(w[a.effect], errors='coerce')
        bad_effect = ~np.isfinite(w[a.effect])
        if bad_effect.any():
            examples = ', '.join(w.loc[bad_effect, 'SNP'].astype(str).head(5))
            raise ValueError(f'Invalid/nonfinite {a.effect} for SNPs {examples}: {src}')
        # COJO can include sequence alleles for indels, e.g. C/CAA. Keep the
        # full sequence: match_variant still requires an exact genotype allele.
        bad_allele = ~w.refA.astype(str).str.fullmatch(r'[ACGT]+')
        if bad_allele.any():
            examples = ', '.join(f'{row.SNP}={row.refA}' for row in w.loc[bad_allele].head(5).itertuples())
            raise ValueError(f'Invalid refA (expected A/C/G/T sequence): {examples}: {src}')
        if not w.Chr.isin(range(1, 23)).all():
            raise ValueError(f'Non-autosomal COJO variants: {src}')
        if not (np.isfinite(w.bp) & (w.bp > 0) & (w.bp == np.floor(w.bp))).all():
            raise ValueError(f'Invalid COJO positions: {src}')
        weights[pop] = w
    chromosomes = sorted(set(int(c) for w in weights.values() for c in w.Chr))
    inputs = list(sources.values()) + [Path(__file__)]
    gen, variants = {}, {}
    for ch in chromosomes:
        prefix = Path(a.dir_gen) / f'chr{ch}'
        if prefix.with_suffix('.pgen').is_file():
            mode = ['--pfile', str(prefix)]
            var = prefix.with_suffix('.pvar')
            if not var.is_file():
                var = prefix.with_suffix('.pvar.zst'); mode += ['vzs']
            files = [prefix.with_suffix('.pgen'), prefix.with_suffix('.psam'), var]
        else:
            mode = ['--bfile', str(prefix)]
            files = [prefix.with_suffix(x) for x in ('.bed', '.bim', '.fam')]
            var = prefix.with_suffix('.bim')
        for f in files:
            if not f.is_file():
                raise FileNotFoundError(f)
        inputs += files
        gen[ch], variants[ch] = mode, var
    if a.remove:
        inputs.append(Path(a.remove))
    stamp = [(str(p.resolve()), p.stat().st_size, p.stat().st_mtime_ns) for p in inputs]
    sig = hashlib.sha256(json.dumps([stamp, a.effect, 'SUM-no-mean-imputation'], sort_keys=True).encode()).hexdigest()
    sf = Path(str(dest) + '.signature')
    audit_file = Path(str(dest) + '.variants.tsv')
    matches_file = Path(str(dest) + '.matches.tsv')
    dest.parent.mkdir(parents=True, exist_ok=True)
    if all(p.is_file() for p in (dest, sf, audit_file, matches_file)) and sf.read_text().strip() == sig:
        print(f'PT: reuse matching cache {dest}', flush=True)
        return
    if not shutil.which('plink2'):
        raise RuntimeError('plink2 not found')
    # Invalidate old scores before updating their audit sidecars, including on failure.
    sf.unlink(missing_ok=True)
    audit, matches, total = [], [], None
    with tempfile.TemporaryDirectory(prefix='yeval-pt-') as scratch:
        scratch = Path(scratch)
        for ch in chromosomes:
            wanted = set(s for w in weights.values() for s in w.loc[w.Chr == ch, 'SNP'])
            prefix = scratch / f'chr{ch}'
            extracted = scratch / f'gen.chr{ch}'
            extract = scratch / 'extract.txt'; extract.write_text('\n'.join(sorted(wanted))+'\n')
            extract_cmd = ['plink2', *gen[ch], '--extract', str(extract), '--make-pgen',
                           '--threads', str(a.threads), '--memory', '4096', '--out', str(extracted)]
            if a.remove:
                extract_cmd += ['--remove', a.remove]
            print('PLINK command: ' + json.dumps(extract_cmd), flush=True)
            print(f'PT chr{ch}: extracting {len(wanted)} distinct COJO variants from imputed genotypes', flush=True)
            subprocess.run(extract_cmd, stderr=subprocess.STDOUT, check=True)
            found = prepare_scoring_index(Path(str(extracted)+'.pvar'), wanted)
            selected, fnames = set(), []
            resolved, ambiguous = 0, 0
            for pop in POPS:
                w = weights[pop].loc[weights[pop].Chr == ch].copy()
                statuses, score_ids = [], []
                duplicates, duplicates_resolved = 0, 0
                for row in w.itertuples(index=False):
                    candidates = found.get(row.SNP, [])
                    variant, status = match_variant(candidates, row.Chr, row.bp, row.refA)
                    statuses.append(status)
                    score_ids.append(variant['score_id'] if variant else None)
                    duplicate = len(candidates) > 1
                    duplicates += duplicate
                    duplicates_resolved += duplicate and status == 'scored'
                    matches.append(dict(population=pop, chromosome=ch, SNP=row.SNP,
                                        bp=row.bp, refA=row.refA, status=status,
                                        genotype_records=len(candidates),
                                        candidates=';'.join(v['variant'] for v in candidates),
                                        selected_variant=variant['variant'] if variant else '',
                                        score_id=variant['score_id'] if variant else ''))
                w['score_id'] = score_ids
                wc = w.loc[w.score_id.notna()]
                audit.append(dict(population=pop, chromosome=ch, source=str(sources[pop]),effect=a.effect,
                                  requested=len(w), scored=len(wc), absent=statuses.count('absent'),
                                  position_mismatch=statuses.count('position_mismatch'),
                                  allele_mismatch=statuses.count('allele_mismatch'),
                                  ambiguous=statuses.count('ambiguous'), duplicate_ids=duplicates,
                                  duplicates_resolved=duplicates_resolved))
                resolved += duplicates_resolved
                ambiguous += statuses.count('ambiguous')
                if wc.empty:
                    continue
                selected.update(wc.SNP)
                wf = scratch / f'pt.{pop}.chr{ch}.tsv'
                wc[['score_id', 'refA', a.effect]].rename(columns={'score_id':'SNP', a.effect:f'pt.{pop}'}).to_csv(wf, sep='\t', index=False)
                fnames.append(str(wf))
            pd.DataFrame(audit).to_csv(audit_file, sep='\t', index=False)
            pd.DataFrame(matches).to_csv(matches_file, sep='\t', index=False)
            if resolved or ambiguous:
                print(f'PT chr{ch}: {resolved} ancestry-specific duplicate-ID weights resolved; '
                      f'{ambiguous} ambiguous weights excluded (see {matches_file})', flush=True)
            if not selected:
                continue
            flist = scratch / 'score-list.txt'; flist.write_text('\n'.join(fnames)+'\n')
            cmd = ['plink2', '--pfile', str(extracted), '--score-list', str(flist), '1', '2', '3',
                   'header-read', 'no-mean-imputation', 'cols=maybefid,scoresums',
                   '--threads', str(a.threads), '--memory', '4096', '--out', str(prefix)]
            print('PLINK command: ' + json.dumps(cmd), flush=True)
            print(f'PT chr{ch}: {len(selected)}/{len(wanted)} distinct COJO variants; {len(fnames)} ancestry scores', flush=True)
            subprocess.run(cmd, stderr=subprocess.STDOUT, check=True)
            sc = pd.read_csv(str(prefix)+'.sscore', sep=r'\s+', dtype={'IID':str,'#IID':str})
            idcol = 'IID' if 'IID' in sc else '#IID'
            if sc[idcol].duplicated().any():
                raise ValueError('Duplicated PLINK score IDs')
            sc = sc.set_index(idcol)
            z = pd.DataFrame(0., index=sc.index, columns=[f'pt.{p}' for p in POPS])
            for pop in POPS:
                column = f'pt.{pop}_SUM'
                if column in sc:
                    z[f'pt.{pop}'] = sc[column]
                elif any(f'pt.{pop}.' in f for f in fnames):
                    raise ValueError(f'Missing PLINK score {column}; got {list(sc.columns)}')
            if total is None:
                total = z
            else:
                if set(total.index) != set(z.index):
                    raise ValueError('Chromosome sample sets differ')
                total += z.reindex(total.index)
    audit = pd.DataFrame(audit)
    audit.to_csv(audit_file, sep='\t', index=False)
    if total is None or any(audit.groupby('population').scored.sum().reindex(POPS).fillna(0)==0):
        raise ValueError(f'An ancestry has no matched PT variants; see {audit_file}')
    total = total.loc[~total.index.str.startswith('-')]
    if total.isna().any().any() or not np.isfinite(total.to_numpy()).all():
        raise ValueError('Nonfinite final PT scores')
    total.index.name = 'eid'
    tmp = Path(str(dest)+'.tmp'); total.to_csv(tmp, sep='\t', compression='gzip'); tmp.replace(dest)
    sf.write_text(sig+'\n')
    print(f'PT published: {dest}; N={len(total)}', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('trait','dir-gwas','dir-gen','output'):
        p.add_argument('--'+name, required=True)
    p.add_argument('--out', help=argparse.SUPPRESS)  # Legacy; audits stay beside the score cache.
    p.add_argument('--effect',choices=['bJ','b'],default='bJ')
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--remove',default='/mnt/d/files/ukb.exclude.id')
    score(p.parse_args())
