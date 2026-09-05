#!/usr/bin/env python3
"""Full, streaming phase/AA/position validation before committing prepared VCF."""
import argparse
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('vcf')
parser.add_argument('--chr', default='')
parser.add_argument('--build', choices=('37', '38'), default='37')
parser.add_argument('--x-psam')
args = parser.parse_args()
with pysam.VariantFile(args.vcf) as vcf:
    ploidies = None
    if args.chr == 'X':
        from arg_x import nonpar, read_sexes
        if args.x_psam:
            sexes = read_sexes(args.x_psam)
            if not set(vcf.header.samples) <= sexes.keys():
                raise SystemExit('X VCF contains samples absent from PSAM')
            ploidies = {n: 1 if sexes[n] == 1 else 2 for n in vcf.header.samples}
    previous = None
    count = 0
    for r in vcf:
        count += 1
        if previous and (r.contig != previous[0] or r.pos <= previous[1]):
            raise SystemExit(f'Non-unique/unsorted positions: {r.contig}:{r.pos}')
        previous = (r.contig, r.pos)
        if args.chr == 'X':
            if r.contig.removeprefix('chr') not in ('X','23') or not nonpar(r.pos, args.build):
                raise SystemExit(f'Expected non-PAR X: {r.contig}:{r.pos}')
            if ploidies is None:
                ploidies = {n: len(c['GT']) for n, c in r.samples.items()}
        aa = r.info.get('AA')
        if len(r.alleles) != 2 or aa not in r.alleles or any(x not in ('A','C','G','T') for x in r.alleles):
            raise SystemExit(f'Invalid alleles/AA: {r.contig}:{r.pos}')
        for name, call in r.samples.items():
            gt = call.get('GT')
            if gt is None or len(gt) not in (1,2): raise SystemExit(f'Invalid ploidy at {r.pos}: {name}')
            if ploidies is not None and len(gt) != ploidies[name]:
                raise SystemExit(f'Incorrect/changing non-PAR X ploidy at {r.pos}: {name}; expected {ploidies[name]}, got {len(gt)}. Supply a VCF with haploid males and diploid females.')
            if any(x is not None and x not in (0,1) for x in gt): raise SystemExit('Invalid GT allele')
            if len(gt)==2 and None not in gt and gt[0]!=gt[1] and not call.phased:
                raise SystemExit(f'Unphased heterozygote at {r.contig}:{r.pos}: {name}')
    if count == 0: raise SystemExit('No variants')
print(f'Validated all {count} VCF records: phased GT, valid AA, unique positions')
