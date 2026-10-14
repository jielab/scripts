"""Validate a completed chromosome before reuse or creation of its done marker."""
import argparse
import json
import math
from pathlib import Path

import h5py
import numpy as np


def read_weights(path, chrom):
    rows = []
    with path.open() as handle:
        for line in handle:
            row = line.split()
            if len(row) != 6 or int(row[0]) != chrom or not math.isfinite(float(row[5])):
                raise ValueError(f'Invalid posterior effect row in {path}')
            rows.append(row)
    if not rows or len({row[1] for row in rows}) != len(rows):
        raise ValueError(f'Empty or duplicate posterior SNPs in {path}')
    return rows


def validate_output(root, chrom, seed, iterations, burnin, thin, posterior=True):
    root = Path(root)
    pops = ('AFR', 'EAS', 'EUR', 'SAS')
    weights = {pop: read_weights(root/pop/f'{pop}.chr{chrom}.pst_eff.txt', chrom) for pop in pops}
    meta = list(root.glob('*META*pst_eff*.txt'))
    if len(meta) != 1:
        raise ValueError('Expected exactly one META posterior')
    meta_rows = read_weights(meta[0], chrom)
    if {row[1] for row in meta_rows} != {row[1] for rows in weights.values() for row in rows}:
        raise ValueError('META/population SNP sets differ')
    if posterior:
        expected = np.arange(thin, iterations+1, thin)
        expected = expected[expected > burnin]
        with h5py.File(root/'joint_posterior.h5', 'r') as handle:
            if not handle.attrs.get('complete', False) or handle.attrs.get('schema') != 'grid_csx_draws_v1':
                raise ValueError('Posterior HDF5 is incomplete or has an unknown schema')
            if int(handle.attrs.get('chromosome', -1)) != chrom or int(handle.attrs.get('seed', -1)) != seed:
                raise ValueError('Posterior chromosome/seed differs from requested chain')
            if not np.array_equal(handle['iteration'][:], expected):
                raise ValueError('Posterior retained iterations differ from requested chain')
            for pop, rows in weights.items():
                group = handle[pop]
                if group['beta'].shape != (len(rows), len(expected)):
                    raise ValueError(f'{pop} posterior shape is incomplete')
                for name, column in (('SNP', 1), ('A1', 3), ('A2', 4)):
                    if list(group[name].asstr()[:]) != [row[column] for row in rows]:
                        raise ValueError(f'{pop} {name} differs between weights and draws')
    return {'chromosome': chrom, 'population_snps': {p: len(w) for p,w in weights.items()},
            'posterior': posterior, 'status': 'valid'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    for name in ('chrom', 'seed', 'iterations', 'burnin', 'thin'):
        p.add_argument('--'+name, type=int, required=True)
    p.add_argument('--posterior', choices=('TRUE', 'FALSE'), default='TRUE')
    p.add_argument('--quiet', action='store_true')
    a = p.parse_args()
    try:
        result = validate_output(a.root, a.chrom, a.seed, a.iterations, a.burnin, a.thin, a.posterior == 'TRUE')
    except (OSError, ValueError, KeyError) as exc:
        p.exit(1, f'Invalid saved chromosome output: {exc}\n')
    if not a.quiet:
        print(json.dumps(result))
