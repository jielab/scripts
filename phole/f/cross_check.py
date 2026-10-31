"""Compare existing mplot signal summaries; never read full GWAS or MAGMA results."""
import argparse
import csv
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

DATASETS = {'ukb': 'prot', 'ckb': 'ckb', 'decode': 'decode', 'scot': 'scot'}


def read_tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_tsv(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def parse_signals(path):
    signals = []
    header = None
    with path.open(encoding='utf-8-sig') as stream:
        for number, line in enumerate(stream, 1):
            line = line.strip()
            if line.lower().startswith('magma:'):
                break
            if not line or line.lower() == 'common:':
                continue
            fields = line.split('\t')
            if fields[0] == 'SNP':
                header = fields
                if not {'CHR', 'POS', 'CLASS'}.issubset(header):
                    raise ValueError(f'{path}:{number}: missing CHR/POS/CLASS header')
                continue
            if header is None or len(fields) != len(header):
                raise ValueError(f'{path}:{number}: malformed signal row')
            row = dict(zip(header, fields))
            chrom = row['CHR'].removeprefix('chr').upper()
            chrom = {'23': 'X', '24': 'Y', '25': 'MT', 'M': 'MT'}.get(chrom, chrom)
            pos = int(row['POS'])
            kind = row['CLASS'].lower()
            if pos < 1 or kind not in ('cis', 'trans'):
                raise ValueError(f'{path}:{number}: invalid POS or CLASS')
            signals.append({'chr': chrom, 'pos': pos, 'kind': kind, 'snp': row['SNP']})
    if header is None:
        raise ValueError(f'{path}: missing signal header (incomplete file)')
    return signals


def index_names(folder):
    """Exact assay/file names only: gene labels must not conflate distinct assays."""
    names = set()
    for filename in ('filename.tsv', 'pheno.tsv'):
        path = folder / filename
        if path.exists():
            for row in read_tsv(path):
                for field in ('filename', 'phenocode', 'phenotype', 'Assay'):
                    value = row.get(field, '')
                    if value:
                        names.add(value[:-3] if value.endswith('.gz') else value)
    return names


def collect(root, mapping):
    rows = read_tsv(mapping)
    if not rows or not {'prot', *DATASETS}.issubset(rows[0]):
        raise ValueError('Mapping requires prot, ukb, ckb, decode, scot columns')
    if len({r['prot'] for r in rows}) != len(rows):
        raise ValueError('Duplicate prot identifiers in mapping')
    indices = {d: index_names(root / folder) for d, folder in DATASETS.items()}
    cache, proteins, coverage, long_rows = {}, [], [], []
    for row in rows:
        protein = {'prot': row['prot'], 'data': {}}
        for dataset, folder in DATASETS.items():
            aliases = list(dict.fromkeys(x.strip() for x in row[dataset].split(';') if x.strip() not in ('', '.')))
            # NA is a valid literal identifier after the one-time mapping migration.
            found, missing, errors, merged = [], [], [], {}
            for alias in aliases:
                if '/' in alias or '\\' in alias or alias in ('.', '..'):
                    raise ValueError(f'Invalid assay name: {alias}')
                path = root / folder / 'mplot' / f'{alias}.sig.txt'
                if not path.is_file():
                    missing.append(alias)
                    continue
                if path not in cache:
                    try:
                        cache[path] = (parse_signals(path), '')
                    except (ValueError, OSError) as exc:
                        cache[path] = ([], str(exc))
                signals, error = cache[path]
                if error:
                    errors.append(error)
                    continue
                found.append(alias)
                for signal in signals:
                    key = (signal['chr'], signal['pos'], signal['kind'])
                    if key not in merged:
                        merged[key] = {**signal, 'assays': [], 'snps': []}
                    merged[key]['assays'].append(alias)
                    merged[key]['snps'].append(signal['snp'])
            status = ('not_mapped' if not aliases else 'invalid' if errors and not found else
                      'missing' if not found else 'partial' if missing or errors else 'complete')
            signals = sorted(merged.values(), key=lambda s: (chrom_key(s['chr']), s['pos'], s['kind']))
            for signal in signals:
                signal['assays'] = sorted(set(signal['assays']))
                signal['snp'] = ';'.join(sorted(set(signal.pop('snps'))))
                long_rows.append({'prot': row['prot'], 'dataset': dataset, 'CHR': signal['chr'],
                                  'POS': signal['pos'], 'CLASS': signal['kind'], 'SNP': signal['snp'],
                                  'assays': ';'.join(signal['assays'])})
            protein['data'][dataset] = {'aliases': aliases, 'status': status, 'signals': signals,
                                       'cis': int(any(s['kind'] == 'cis' for s in signals)),
                                       'trans': sum(s['kind'] == 'trans' for s in signals)}
            coverage.append({'prot': row['prot'], 'dataset': dataset, 'status': status,
                             'aliases': ';'.join(aliases) or '.', 'files_read': ';'.join(found) or '.',
                             'files_missing': ';'.join(missing) or '.',
                             'not_in_index': ';'.join(a for a in aliases if a not in indices[dataset]) or '.',
                             'errors': ';'.join(errors) or '.'})
        proteins.append(protein)
    return proteins, coverage, long_rows


def chrom_key(chrom):
    return (0, int(chrom)) if chrom.isdigit() else (1, {'X': '01', 'Y': '02', 'MT': '03'}.get(chrom, chrom))


def static_plots(out, proteins, axes):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    offsets = {a['chr']: a['offset'] for a in axes}
    total = sum(a['span'] for a in axes)
    ticks = [a['offset'] + a['span'] / 2 for a in axes]
    labels = [a['chr'] for a in axes]
    for primary in DATASETS:
        selected = [p for p in proteins if p['data'][primary]['cis']]
        fig, panels = plt.subplots(1, 4, figsize=(13, 4), sharey=True)
        for ax, dataset in zip(panels, DATASETS):
            counts = Counter('cis present' if p['data'][dataset]['cis'] else
                             'no cis' if p['data'][dataset]['status'] == 'complete' else
                             'unavailable / partial' for p in selected)
            keys = ['cis present', 'no cis', 'unavailable / partial']
            bars = ax.bar(range(3), [counts[k] for k in keys], color=['#dc4b52', '#78b9aa', '#c8d2dc'])
            ax.bar_label(bars)
            ax.set_xticks(range(3), keys, rotation=25, ha='right', fontsize=8)
            ax.set_title(dataset.upper())
        panels[0].set_ylabel('Proteins')
        fig.suptitle(f'Primary {primary.upper()}: {len(selected):,} proteins with observed cis signals')
        fig.tight_layout()
        fig.savefig(out / f'primary-{primary}.overview.png', dpi=160)
        plt.close(fig)
        if not selected:
            # Remove only our own obsolete preview if a primary became empty.
            (out / f'primary-{primary}.first-protein.png').unlink(missing_ok=True)
            continue
        protein = selected[0]
        fig, panels = plt.subplots(4, 1, figsize=(14, 6), sharex=True)
        for ax, dataset in zip(panels, DATASETS):
            data = protein['data'][dataset]
            for kind, color in [('trans', '#15976b'), ('cis', '#dc4b52')]:
                positions = [offsets[s['chr']] + s['pos'] for s in data['signals'] if s['kind'] == kind]
                ax.vlines(positions, 0, 1, color=color, label=kind)
            ax.set_ylabel(dataset.upper())
            ax.set_yticks([])
            ax.set_title(data['status'], fontsize=9, loc='right')
            ax.set_xlim(0, total)
            for tick in axes:
                ax.axvline(tick['offset'], color='#dde4eb', linewidth=.5)
        panels[-1].set_xticks(ticks, labels)
        panels[-1].set_xlabel('Chromosome; original CHR/POS, no liftover; red = cis, green = trans')
        panels[0].legend(loc='upper left', ncol=2)
        fig.suptitle(f'{protein["prot"]} | Primary {primary.upper()} | first protein preview')
        fig.tight_layout()
        fig.savefig(out / f'primary-{primary}.first-protein.png', dpi=160)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gwas-root', type=Path, default=Path('/mnt/f/gwas'))
    parser.add_argument('--protein-map', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--primary', choices=DATASETS, default='ukb')
    args = parser.parse_args()
    root = args.gwas_root
    mapping = args.protein_map or root / 'prot.names.tsv'
    proteins, coverage, signals = collect(root, mapping)
    args.out.mkdir(parents=True, exist_ok=True)
    fields = ['prot'] + [f'{kind}-{d}' for kind in ('cis', 'trans') for d in DATASETS]
    tables = [{ 'prot': p['prot'], **{f'{kind}-{d}': p['data'][d][kind]
               for kind in ('cis', 'trans') for d in DATASETS}} for p in proteins]
    write_tsv(args.out / 'cross-check.tsv', fields, tables)
    for primary in DATASETS:
        write_tsv(args.out / f'primary-{primary}.tsv', fields, [r for r in tables if r[f'cis-{primary}']])
    write_tsv(args.out / 'coverage.tsv', ['prot', 'dataset', 'status', 'aliases', 'files_read', 'files_missing', 'not_in_index', 'errors'], coverage)
    write_tsv(args.out / 'signals.tsv', ['prot', 'dataset', 'CHR', 'POS', 'CLASS', 'SNP', 'assays'], signals)
    spans = {}
    for signal in signals:
        spans[signal['CHR']] = max(spans.get(signal['CHR'], 0), signal['POS'])
    axes, offset = [], 0
    for chrom in sorted(spans, key=chrom_key):
        span = int(spans[chrom] * 1.02) + 1
        axes.append({'chr': chrom, 'span': span, 'offset': offset})
        offset += span
    payload = {'proteins': proteins, 'axes': axes, 'primary': args.primary,
               'generated': datetime.now().isoformat(timespec='seconds'), 'root': str(root)}
    template = Path(__file__).with_name('cross_check.html').read_text(encoding='utf-8')
    html = template.replace('__PAYLOAD__', json.dumps(payload, ensure_ascii=False, separators=(',', ':')).replace('<', '\\u003c'))
    (args.out / 'cross-check.html').write_text(html, encoding='utf-8')
    static_plots(args.out, proteins, axes)
    summary = {d: {'primary_cis_proteins': sum(p['data'][d]['cis'] for p in proteins),
                   'coverage': dict(Counter(p['data'][d]['status'] for p in proteins))} for d in DATASETS}
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))
    print(f'Report: {args.out / "cross-check.html"}')


if __name__ == '__main__':
    main()
