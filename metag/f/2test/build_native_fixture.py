#!/usr/bin/env python3
"""Build a fresh, tiny two-species fixture with the bundled real CPU programs.

Run on a compute node, or explicitly --allow-local for a small developer test.
The destination must be new or empty. No production data is read or modified.
The companion native_query_smoke.py consumes the resulting config.json.
"""
import argparse
import json
from pathlib import Path
import random
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--code-root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--allow-local', action='store_true')
    args = parser.parse_args()
    sys.path.insert(0, str(args.code_root.resolve() / 'f'))
    from hpc_common import load_config, require_compute, atomic_json
    import hpc_install
    import hpc_kunpeng_index
    import hpc_reference
    require_compute(args.allow_local)
    base = args.destination.expanduser().absolute()
    if base.is_symlink() or (base.exists() and (not base.is_dir() or any(base.iterdir()))):
        raise RuntimeError('Native test destination must be new or empty')
    base.mkdir(parents=True, exist_ok=True)
    for name in ('raw', 'reads', 'source', 'metadata', 'lib', 'work', 'output'):
        (base / name).mkdir()
    rng = random.Random(41)
    for taxid in (10, 20):
        sequence = ''.join(rng.choices('ACGT', k=100000))
        (base / 'raw' / f'species{taxid}.fna').write_text(
            f'>ref{taxid}\n' + '\n'.join(sequence[i:i+80] for i in range(0, len(sequence), 80)) + '\n')
    nodes = [(1, 1, 'no rank'), (2, 1, 'superkingdom'), (3, 2, 'genus'),
             (10, 3, 'species'), (20, 3, 'species')]
    (base / 'source/nodes.dmp').write_text(''.join(
        f'{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n'
        for taxid, parent, rank in nodes))
    (base / 'source/names.dmp').write_text(''.join(
        f'{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n'
        for taxid, name in [(1, 'root'), (2, 'Bacteria'), (3, 'Fixture'),
                            (10, 'Fixture one'), (20, 'Fixture two')]))
    (base / 'source/merged.dmp').write_text('')
    (base / 'map.tsv').write_text('path\ttaxid\nspecies10.fna\t10\nspecies20.fna\t20\n')
    settings = {'root': str(base), 'library_dir': str(base / 'lib'),
                'metadata_dir': str(base / 'metadata'), 'reference_raw': str(base / 'raw'),
                'reads_raw': str(base / 'reads'), 'work_dir': str(base / 'work'),
                'output_dir': str(base / 'output'), 'tools_dir': str(base / 'tools'),
                'reference_taxid_map': str(base / 'map.tsv'), 'min_free_gib': 0,
                'masking': True, 'prepare_threads': 2, 'prepare_target_gib': 0.0001,
                'prepare_max_files': 1, 'build_threads': 2, 'build_memory_gib': 16,
                'kunpeng_hash_capacity': 1024, 'convert_memory_gib': 4,
                'setup_cores': 2, 'kunpeng_install_mode': 'binary',
                'database_layout': 'single', 'legacy_library_dir': None, 'legacy_metadata_dir': None}
    config = base / 'config.json'
    atomic_json(config, settings)
    cfg = load_config(config)
    hpc_install.install(cfg)
    hpc_kunpeng_index.install(cfg)
    hpc_reference.metadata(cfg, source=base / 'source', download=False)
    plan = hpc_reference.plan_references(cfg)
    for task in plan['tasks']:
        hpc_reference.prepare_task(cfg, task['task_id'])
    original = hpc_reference.build_database(cfg)
    converted = hpc_kunpeng_index.convert_database(cfg)
    report = {'status': 'pass', 'config': str(config), 'reference_bases': 200000,
              'kraken_identity': original['identity'], 'kunpeng_identity': converted['identity'],
              'kraken_hash_bytes': original['files']['hash.k2d']['state']['size'],
              'note': 'Small synthetic correctness fixture; no HPC throughput claim.'}
    atomic_json(base / 'fixture-result.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
