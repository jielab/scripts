#!/usr/bin/env python3
"""Lift BIM SNP coordinates without changing BED row or allele encoding."""
import collections
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    source, outdir, chain, fasta = map(Path, sys.argv[1:])
    source = source.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    output = outdir / source.name
    if output.resolve() == source:
        raise ValueError('Input and output must differ')
    for path in (source, chain, fasta, Path(str(fasta)+'.fai')):
        if not path.is_file():
            raise FileNotFoundError(path)
    companions = [source.with_suffix(ext) for ext in ('.bed', '.fam', '.nosex')]
    for path in companions:
        if not path.exists():
            raise FileNotFoundError(path)
        dst = outdir / path.name
        if dst.exists() and not dst.is_symlink():
            raise FileExistsError(dst)
    rows = [line.split() for line in source.read_text().splitlines()]
    if not rows or any(len(row) != 6 for row in rows):
        raise ValueError('Expected nonempty six-column BIM')
    lift = os.environ.get('LIFTOVER_BIN') or shutil.which('liftOver') or '/mnt/d/software/bin/liftOver'
    with tempfile.TemporaryDirectory(prefix='.liftGen.bim.', dir=outdir) as tmp:
        tmp = Path(tmp)
        bed, mapped, unmapped = [tmp / name for name in ('source.bed', 'mapped.bed', 'unmapped.bed')]
        expected = {}
        with bed.open('w') as handle:
            for index, row in enumerate(rows):
                chrom = row[0].removeprefix('chr')
                chrom = {'23':'X', '24':'Y', '25':'XY', '26':'MT', 'M':'MT'}.get(chrom, chrom)
                chrom = 'chr'+chrom
                # Non-SNP alleles need normalization; do not guess their encoding.
                if int(row[3]) > 0 and all(a.upper() in ('A','C','G','T') for a in row[4:6]):
                    pos = int(row[3])
                    expected[index] = chrom
                    handle.write(f'{chrom}\t{pos-1}\t{pos}\t{index}\t0\t+\n')
        subprocess.run([lift, str(bed), str(chain), str(mapped), str(unmapped)], check=True)
        groups = collections.defaultdict(list)
        for line in mapped.read_text().splitlines():
            fields = line.split()
            groups[int(fields[3])].append(fields)
        positions = {}
        candidates = tmp / 'candidate.bed'
        with candidates.open('w') as handle:
            for index, records in groups.items():
                if len(records) != 1:
                    continue
                chrom, start, end, _, _, strand = records[0]
                if chrom == expected[index] and strand == '+' and int(end)-int(start) == 1:
                    positions[index] = int(start)+1
                    handle.write(f'{chrom}\t{start}\t{end}\t{index}\n')
        verified = {}
        if positions:
            result = subprocess.run(['bedtools','getfasta','-fi',str(fasta),'-bed',str(candidates),'-nameOnly','-tab'],
                                    check=True, capture_output=True, text=True)
            for line in result.stdout.splitlines():
                index, allele = line.split('\t')
                index = int(index)
                if allele.upper() in [a.upper() for a in rows[index][4:6]]:
                    verified[index] = positions[index]
        staged = tmp / source.name
        with staged.open('w') as handle:
            for index, row in enumerate(rows):
                row[3] = str(verified.get(index, -1))
                handle.write('\t'.join(row)+'\n')
        os.replace(staged, output)
        for path in companions:
            dst = outdir / path.name
            link = tmp / path.name
            link.symlink_to(os.path.relpath(path, outdir.resolve()))
            os.replace(link, dst)
        report = (f'input\t{source}\noutput\t{output}\ntotal_variants\t{len(rows)}\n'
                  f'unique_forward_same_chr_length\t{len(positions)}\n'
                  f'verified_either_allele\t{len(verified)}\nfailed_pos_minus_1\t{len(rows)-len(verified)}\n')
        log = tmp / 'complete.log'
        log.write_text(report)
        os.replace(log, str(output)+'.liftGen.log')
        print(report, end='')


if __name__ == '__main__':
    main()
