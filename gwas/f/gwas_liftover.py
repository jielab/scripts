#!/usr/bin/env python3
"""Stream a complete standardized GWAS through UCSC liftOver, retaining X/Y."""
import argparse
from collections import Counter
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from gwas_alleles import IndexedFasta, id_alleles, normalize

PRIMARY_CHROMOSOMES = frozenset(str(c) for c in range(1, 26))

# GRC PAR intervals, inclusive and 1-based. Keep source-X PAR associations on
# X when a chain picks the homologous Y representation. True source-Y and
# non-PAR mappings are not changed.
PAR = {
    37: ((60001, 2699520, 10001, 2649520), (154931044, 155260560, 59034050, 59363566)),
    38: ((10001, 2781479, 10001, 2781479), (155701383, 156030895, 56887903, 57217415)),
}


def source_x_par_target(source_chr, source_pos, source_build, target_build, mapping):
    if chrom(source_chr) != '23' or chrom(mapping[0]) != '24':
        return mapping
    start, end = int(mapping[1]) + 1, int(mapping[2])
    for old, new in zip(PAR[source_build], PAR[target_build]):
        if old[0] <= int(source_pos) <= old[1] and new[2] <= start <= end <= new[3]:
            result = list(mapping)
            offset = new[0] - new[2]
            result[0], result[1], result[2] = 'chrX', str(start - 1 + offset), str(end + offset)
            return result
    return mapping


def identity(path):
    p = Path(path).resolve()
    s = p.stat()
    return [str(p), s.st_size, s.st_mtime_ns]


def chrom(value):
    value = value.upper().removeprefix('CHR')
    return {'X': '23', 'Y': '24', 'M': '25', 'MT': '25'}.get(value, value)


def ucsc(value):
    return 'chr' + {'23': 'X', '24': 'Y', '25': 'M'}.get(chrom(value), chrom(value))


def reverse_complement(value):
    if not value or any(c not in 'ACGT' for c in value.upper()):
        raise ValueError('non-DNA allele on reverse strand')
    return value.upper().translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def run(a):
    src, dst = Path(a.input), Path(a.output)
    if src.resolve() == dst.resolve():
        raise ValueError('Source and target must be distinct; retain the source build')
    qc = Path(a.qc_prefix)
    qc.parent.mkdir(parents=True, exist_ok=True)
    dst.parent.mkdir(parents=True, exist_ok=True)
    state = Path(str(qc) + '.liftover.state.json')
    # Reference files are optional for legacy SNV-only use, but mandatory for
    # preserving sequence-resolved indels. Never interpret a point lift as an
    # indel lift without checking its full reference interval.
    fasta_paths = []
    for option, build in (('source_fasta', a.source_build), ('target_fasta', a.target_build)):
        path = getattr(a, option, None)
        if not path:
            candidate = Path(f'/mnt/e/refGen/fasta/GRCH{build}.fasta')
            path = str(candidate) if candidate.is_file() and Path(str(candidate)+'.fai').is_file() else None
        fasta_paths.append(path)
    signature = dict(version=5, source=identity(src), chain=identity(a.chain),
                     source_build=a.source_build, target_build=a.target_build,
                     fastas=[identity(p) if p else None for p in fasta_paths])
    if a.replace == 'FALSE' and dst.exists() and state.exists():
        saved = json.loads(state.read_text())
        if saved.get('signature') == signature and saved.get('output') == identity(dst):
            print('Reuse verified liftover: ' + str(dst), flush=True)
            return
    # Invalidate before work; a failed run must not leave a success state.
    state.unlink(missing_ok=True)
    env = dict(os.environ, LC_ALL='C')
    counts = Counter()
    by_chr = Counter()
    source_fasta, target_fasta = [IndexedFasta(p) if p else None for p in fasta_paths]
    sequence_rows = {}
    with tempfile.TemporaryDirectory(prefix='liftover.', dir=qc.parent) as temp:
        work = Path(temp)
        bed, old = work / 'input.bed', work / 'input.tsv'
        opener = gzip.open if str(src).endswith(('.gz', '.bgz')) else open
        with opener(src, 'rt') as inp, bed.open('w') as bp, old.open('w') as op:
            header = inp.readline().rstrip('\r\n').split('\t')
            ci, pi, ai, bi = [header.index(c) for c in ('CHR', 'POS', 'EA', 'NEA')]
            for n, line in enumerate(inp, 1):
                row = line.rstrip('\r\n').split('\t')
                if len(row) != len(header):
                    raise ValueError(f'Malformed source row {n}')
                pos = int(row[pi])
                if pos < 1:
                    raise ValueError(f'Invalid position at row {n}')
                if chrom(row[ci]) == '25' and pos > 16569:
                    raise ValueError('CHR=25 exceeds mitochondrial length. If the source uses PLINK '
                                     '25=XY/PAR, standardize it to X (23) before liftover.')
                key = f'{n:012d}'
                length = 1
                if max(len(row[ai]), len(row[bi])) > 1 and source_fasta and target_fasta:
                    pair = id_alleles(row[header.index('SNP')], pos)
                    if pair and set(pair) == {row[ai], row[bi]}:
                        ref, alt = pair
                        if source_fasta.fetch(row[ci], pos, len(ref)) == ref:
                            sequence_rows[key] = (ref, alt, row[ai] == ref)
                            length = len(ref)
                        else:
                            sequence_rows[key] = 'source_reference_mismatch'
                bp.write(f'{ucsc(row[ci])}\t{pos-1}\t{pos+length-1}\t{key}\t0\t+\n')
                op.write(key + '\t' + '\t'.join(row) + '\n')
                counts['input'] += 1
                by_chr[(chrom(row[ci]), 'input')] += 1
        lifted, unmapped = work / 'lifted.bed', work / 'unmapped.bed'
        with Path(str(qc) + '.liftover.log').open('w') as log:
            subprocess.run([a.liftOver, str(bed), a.chain, str(lifted), str(unmapped)],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        ordered = work / 'mapped.sorted.bed'
        with ordered.open('w') as out:
            subprocess.run(['sort', '-T', temp, '-S', '256M', '-k4,4', str(lifted)],
                           stdout=out, env=env, check=True)
        mapped = work / 'mapped.tsv'
        with ordered.open() as lp, old.open() as op, mapped.open('w') as out, \
                gzip.open(str(qc) + '.liftover.par.tsv.gz', 'wt') as par_audit, \
                gzip.open(str(qc) + '.liftover.nonprimary.tsv.gz', 'wt') as nonprimary, \
                gzip.open(str(qc) + '.liftover.alleles.tsv.gz', 'wt') as allele_audit:
            allele_audit.write('SNP\tCHR_SOURCE\tPOS_SOURCE\tEA_SOURCE\tNEA_SOURCE\tREASON\n')
            nonprimary.write('\t'.join(header + ['TARGET_CONTIG', 'TARGET_START_0', 'TARGET_END_0',
                                                'TARGET_STRAND', 'REASON']) + '\n')
            par_audit.write('SNP\tCHR_SOURCE\tPOS_SOURCE\tCHAIN_TARGET\tCHAIN_POS\tFINAL_TARGET\tFINAL_POS\tREASON\n')
            current = next(lp, '').split()
            for line in op:
                key, rest = line.rstrip('\n').split('\t', 1)
                row = rest.split('\t')
                source_chr = chrom(row[ci])
                if not current or current[3] != key:
                    counts['unmapped'] += 1
                    by_chr[(source_chr, 'unmapped')] += 1
                    continue
                if len(current) != 6:
                    raise ValueError('Invalid mapped BED6 interval')
                canonical = source_x_par_target(source_chr, row[pi], a.source_build, a.target_build, current)
                if canonical != current:
                    par_audit.write('\t'.join([row[header.index('SNP')], source_chr, row[pi], current[0],
                        str(int(current[1])+1), canonical[0], str(int(canonical[1])+1),
                        'source_X_PAR_retained_on_homologous_target_X_PAR'])+'\n')
                    current = canonical
                    counts['source_x_par_retained_on_x'] += 1
                if chrom(current[0]) not in PRIMARY_CHROMOSOMES:
                    # The project schema/LD references support primary chromosomes.
                    # Numeric sorting would interleave chr1 with chr1_* contigs.
                    # Preserve the entire source association and mapping separately;
                    # do not mislabel these mappings as an autosomal position.
                    nonprimary.write('\t'.join(row + [current[0], current[1], current[2], current[5],
                                      'nonprimary_target_unsupported_by_project_schema_and_LD_references']) + '\n')
                    counts['nonprimary_excluded'] += 1
                    by_chr[(source_chr, 'nonprimary_excluded')] += 1
                    previous = current[3]
                    current = next(lp, '').split()
                    if current and current[3] == previous:
                        raise ValueError('Multiple mappings for the same source variant')
                    continue
                seq = sequence_rows.get(key)
                reject = None
                new_pos = int(current[1]) + 1
                expected_length = len(seq[0]) if isinstance(seq, tuple) else 1
                if int(current[2]) - int(current[1]) != expected_length:
                    reject = 'reference_interval_changed_length'
                if isinstance(seq, str):
                    reject = seq
                if isinstance(seq, tuple) and reject is None:
                    ref, alt, effect_is_ref = seq
                    if current[5] == '-':
                        ref, alt = reverse_complement(ref), reverse_complement(alt)
                    target_chr = chrom(current[0])
                    if target_fasta.fetch(target_chr, new_pos, len(ref)) != ref:
                        reject = 'target_reference_mismatch'
                    else:
                        new_pos, ref, alt = normalize(new_pos, ref, alt, target_fasta, target_chr)
                        row[ai], row[bi] = (ref, alt) if effect_is_ref else (alt, ref)
                        counts['sequence_alleles_lifted'] += 1
                if reject:
                    allele_audit.write('\t'.join([row[header.index('SNP')], row[ci], row[pi], row[ai], row[bi],
                                                  reject + '; association_excluded']) + '\n')
                    counts['alleles_excluded'] += 1
                    by_chr[(source_chr, 'alleles_excluded')] += 1
                    previous = current[3]
                    current = next(lp, '').split()
                    if current and current[3] == previous:
                        raise ValueError('Multiple mappings for the same source variant')
                    continue
                # Unknown allele codes or indels lacking explicit REF/ALT and
                # FASTA validation retain association statistics with NA alleles.
                unresolved = any(len(row[i]) != 1 or row[i].upper() not in 'ACGT' for i in (ai, bi))
                if isinstance(seq, tuple):
                    unresolved = False
                if unresolved:
                    allele_audit.write('\t'.join([row[header.index('SNP')], row[ci], row[pi], row[ai], row[bi],
                                                  'alleles_not_resolved_by_point_liftover; association_retained']) + '\n')
                    row[ai] = row[bi] = 'NA'
                    counts['alleles_unresolved'] += 1
                    by_chr[(source_chr, 'alleles_unresolved')] += 1
                try:
                    if current[5] == '-' and not unresolved and not isinstance(seq, tuple):
                        row[ai], row[bi] = reverse_complement(row[ai]), reverse_complement(row[bi])
                        counts['reverse_complemented'] += 1
                    row[ci], row[pi] = chrom(current[0]), str(new_pos)
                    out.write('\t'.join(row) + '\n')
                    counts['lifted'] += 1
                    by_chr[(source_chr, 'lifted')] += 1
                except ValueError:
                    counts['unsupported_reverse_allele'] += 1
                    by_chr[(source_chr, 'unsupported_reverse_allele')] += 1
                previous = current[3]
                current = next(lp, '').split()
                if current and current[3] == previous:
                    raise ValueError('Multiple mappings for the same source variant')
            if current:
                raise ValueError('Unconsumed mapping rows')
        if not counts['lifted']:
            raise ValueError('No variants lifted; inspect chain/build and liftover.log')
        tmp = Path(str(dst) + '.tmp')
        with tmp.open('wb') as out:
            bg = subprocess.Popen(['bgzip', '-@', '4', '-c'], stdin=subprocess.PIPE, stdout=out)
            try:
                bg.stdin.write(('\t'.join(header) + '\n').encode())
                with subprocess.Popen(['sort', '-T', temp, '-S', '512M', '-t', '\t',
                                       f'-k{ci+1},{ci+1}n', f'-k{pi+1},{pi+1}n', '-k1,1',
                                       str(mapped)], stdout=subprocess.PIPE, env=env) as sorter:
                    shutil.copyfileobj(sorter.stdout, bg.stdin, 1024 * 1024)
                    if sorter.wait():
                        raise RuntimeError('Coordinate sort failed')
                bg.stdin.close()
                if bg.wait():
                    raise RuntimeError('BGZF compression failed')
            except BaseException:
                bg.kill(); bg.wait()
                raise
        subprocess.run(['tabix', '-s', str(ci+1), '-b', str(pi+1), '-e', str(pi+1),
                        '-S', '1', str(tmp)], check=True)
        os.replace(tmp, dst)
        os.replace(str(tmp) + '.tbi', str(dst) + '.tbi')
        Path(str(dst) + '.csi').unlink(missing_ok=True)
        with gzip.open(str(qc) + '.liftover.unmapped.bed.gz', 'wb') as out, unmapped.open('rb') as inp:
            shutil.copyfileobj(inp, out)
    with Path(str(qc) + '.liftover.n.tsv').open('w') as out:
        out.write('CHR\tN_INPUT\tN_LIFTED\tN_UNMAPPED\tN_ALLELES_UNRESOLVED\tN_ALLELES_EXCLUDED\tN_NONPRIMARY_EXCLUDED\n')
        for c in sorted({c for c, _ in by_chr}, key=lambda c: int(c)):
            out.write(c + '\t' + '\t'.join(str(by_chr[c, k]) for k in
                      ('input', 'lifted', 'unmapped', 'alleles_unresolved', 'alleles_excluded', 'nonprimary_excluded')) + '\n')
    for fasta in (source_fasta, target_fasta):
        if fasta:
            fasta.close()
    Path(str(src) + '.grch').write_text(str(a.source_build) + '\n')
    Path(str(dst) + '.grch').write_text(str(a.target_build) + '\n')
    Path(str(qc) + '.grch').write_text(str(a.target_build) + '\n')
    state.write_text(json.dumps(dict(signature=signature, output=identity(dst), counts=counts), indent=2))
    print('Liftover completed: ' + json.dumps(counts), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'output', 'chain', 'qc-prefix'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--liftOver', default='liftOver')
    p.add_argument('--source-fasta', help='Indexed source FASTA for sequence-resolved indels')
    p.add_argument('--target-fasta', help='Indexed target FASTA for validation and left alignment')
    p.add_argument('--source-build', type=int, choices=(37, 38), required=True)
    p.add_argument('--target-build', type=int, choices=(37, 38), required=True)
    p.add_argument('--replace', choices=('TRUE', 'FALSE'), default='FALSE')
    run(p.parse_args())
