#!/usr/bin/env python3
"""Pack real non-PAR X copies into ASMC's two-column storage containers."""
import argparse
import csv
import json
from pathlib import Path

from arg_x import nonpar, read_sexes


def panel(keep, psam, policy, qc):
    sexes = read_sexes(psam)
    path = Path(keep)
    rows = path.read_text().splitlines()
    names = [r.split()[-1] for r in rows[1:] if r.strip()]
    if not set(names) <= sexes.keys():
        raise ValueError('Selected X panel has samples absent from PSAM')
    males = [name for name in names if sexes[name] == 1]
    dropped = []
    if len(males) % 2:
        if policy != 'drop-last':
            raise ValueError('ASMC requires an even number of real haplotypes. '
                             'Use --x-odd-male drop-last to exclude the last selected male, '
                             'or provide an even-male sample panel with --keep.')
        dropped = males[-1:]
        rows = [rows[0]] + [r for r in rows[1:] if r.split()[-1] not in dropped]
        path.write_text('\n'.join(rows) + '\n')
    result = dict(policy=policy, excluded_males=dropped, individuals=len(names)-len(dropped),
                  haplotypes=sum(1 if sexes[n] == 1 else 2 for n in names)-len(dropped))
    qc = Path(qc)
    text = json.dumps(result, indent=2) + '\n'
    if not qc.is_file() or qc.read_text() != text:
        temporary = qc.with_name(qc.name + '.next')
        temporary.write_text(text)
        temporary.replace(qc)
    print('Needle X panel: ' + json.dumps(result))


def pack(prefix, psam, build, keep):
    import numpy as np
    prefix = str(prefix)
    sexes = read_sexes(psam)
    sample = Path(prefix + '.sample')
    names = [line.split()[1] for line in sample.read_text().splitlines()[2:] if line.strip()]
    wanted = [line.split()[-1] for line in Path(keep).read_text().splitlines() if line.strip() and not line.startswith('#')]
    if names != wanted:
        raise ValueError('X export sample order differs from selected panel')
    identities = [(name, h) for name in names for h in range(1 if sexes[name] == 1 else 2)]
    if len(identities) % 2:
        raise ValueError('Odd real X haplotype count: select an even-male panel before export')
    male = np.array([sexes[name] == 1 for name in names])
    columns = np.array([2*i+h for i, name in enumerate(names) for h in range(1 if sexes[name] == 1 else 2)])
    haps = Path(prefix + '.haps')
    tmp = Path(prefix + '.haps.x.next')
    retained = removed = 0
    try:
        with haps.open() as src, tmp.open('w') as dest:
            for line in src:
                z = line.split()
                if len(z) != 5 + 2 * len(names) or z[0].removeprefix('chr') not in ('X', '23'):
                    raise ValueError('Malformed X HAPS row')
                if not nonpar(int(z[2]), build):
                    removed += 1
                    continue
                # Vectorize across individuals: full-density X contains billions
                # of calls. ASCII 0/1 are 48/49; PLINK's absent male copy is '-'.
                values = np.frombuffer(''.join(z[5:]).encode('ascii'), dtype=np.uint8)
                if len(values) != 2 * len(names):
                    raise ValueError(f'Missing/unphased X HAPS at {z[2]}')
                pairs = values.reshape(-1, 2)
                a, b = pairs[:, 0], pairs[:, 1]
                valid_a = (a == 48) | (a == 49)
                valid_b = np.where(male, (b == 45) | (b == a), (b == 48) | (b == 49))
                if not np.all(valid_a & valid_b):
                    bad = int(np.flatnonzero(~(valid_a & valid_b))[0])
                    raise ValueError(f'Invalid/missing X HAPS: {names[bad]} at {z[2]}')
                gt = values[columns].tobytes().decode('ascii')
                dest.write(' '.join(['X'] + z[1:5]) + ' ' + ' '.join(gt) + '\n')
                retained += 1
        if retained < 2:
            raise ValueError('Too few non-PAR X HAPS variants')
        tmp.replace(haps)
    finally:
        tmp.unlink(missing_ok=True)
    # Pairs are file-format containers only; biological identities live in the sidecar.
    storage_ids = [f'X_pair_{i}' for i in range(len(identities)//2)]
    sample.write_text('ID_1 ID_2 missing\n0 0 0\n' + ''.join(f'{n} {n} 0\n' for n in storage_ids))
    Path(prefix + '.packed.keep').write_text('#IID\n' + '\n'.join(storage_ids) + '\n')
    with open(prefix + '.haplotypes.tsv', 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('eid', 'hap'))
        writer.writerows(identities)
    print(f'Needle X packed: individuals={len(names)} real_haplotypes={len(identities)} '
          f'nonpar_variants={retained} excluded_outside_nonpar={removed}')


def sample_map(tree, identities, output):
    import tskit
    from arg_x import restore_individuals

    with open(identities) as handle:
        rows = [(r['eid'], int(r['hap'])) for r in csv.DictReader(handle, delimiter='\t')]
    ts = tskit.load(tree)
    if len(rows) != ts.num_samples:
        raise ValueError('X tree sample count differs from biological haplotype map')
    restored = restore_individuals(ts, rows)
    temp = Path(str(tree) + '.x.next')
    restored.dump(str(temp))
    temp.replace(tree)
    with open(output, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('eid', 'hap', 'sample_node'))
        writer.writerows((name, hap, int(node)) for (name, hap), node in zip(rows, restored.samples()))


def annotate_qc(qc, identities):
    from collections import Counter
    with open(identities) as handle:
        counts = Counter(r['eid'] for r in csv.DictReader(handle, delimiter='\t'))
    path = Path(qc)
    result = json.loads(path.read_text())
    if result['haplotypes'] != sum(counts.values()):
        raise ValueError('X HAPS QC count differs from real haplotypes')
    result['storage_pairs'] = result['individuals']
    result['individuals'] = len(counts)
    result['haploid_individuals'] = sum(n == 1 for n in counts.values())
    result['diploid_individuals'] = sum(n == 2 for n in counts.values())
    result['storage'] = 'ASMC-pairs; biological identities in chrX.haplotypes.tsv'
    path.write_text(json.dumps(result, indent=2) + '\n')


def check(tree, identities, output, psam, keep):
    import tskit
    sexes = read_sexes(psam)
    names = [r.split()[-1] for r in Path(keep).read_text().splitlines()[1:] if r.strip()]
    expected = [(name, h) for name in names for h in range(1 if sexes[name] == 1 else 2)]
    with open(identities) as handle:
        actual = [(r['eid'], int(r['hap'])) for r in csv.DictReader(handle, delimiter='\t')]
    if actual != expected or len(set(actual)) != len(actual):
        raise ValueError('X biological haplotype identities disagree with panel/sex')
    ts = tskit.load(tree)
    with open(output) as handle:
        mapped = [(r['eid'], int(r['hap']), int(r['sample_node'])) for r in csv.DictReader(handle, delimiter='\t')]
    if len(expected) != ts.num_samples or mapped != [(n, h, int(node)) for (n, h), node in zip(expected, ts.samples())]:
        raise ValueError('X sample map does not match tree samples and real ploidy')
    for name, hap, node in mapped:
        individual = ts.node(node).individual
        if individual < 0 or ts.individual(individual).metadata.get('sample') != name:
            raise ValueError('X tree individual identity disagrees with sample map')
    print(f'Needle X identity CHECK PASS: individuals={len(names)} real_haplotypes={len(expected)}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('panel')
    for key in ('keep', 'psam', 'qc'):
        p.add_argument('--' + key, required=True)
    p.add_argument('--policy', choices=('error', 'drop-last'), default='drop-last')
    p = sub.add_parser('pack')
    for key in ('prefix', 'psam', 'keep', 'build'):
        p.add_argument('--' + key, required=True)
    p = sub.add_parser('map')
    for key in ('tree', 'identities', 'output'):
        p.add_argument('--' + key, required=True)
    p = sub.add_parser('check')
    for key in ('tree', 'identities', 'output', 'psam', 'keep'):
        p.add_argument('--' + key, required=True)
    p = sub.add_parser('qc')
    for key in ('qc', 'identities'):
        p.add_argument('--' + key, required=True)
    args = vars(parser.parse_args())
    command = args.pop('command')
    {'panel': panel, 'pack': pack, 'map': sample_map, 'check': check, 'qc': annotate_qc}[command](**args)


if __name__ == '__main__':
    main()
