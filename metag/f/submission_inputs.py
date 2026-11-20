"""Metadata-only admission checks for independent submissions (no scheduler waits)."""
import argparse
import json
from pathlib import Path


def require_files(paths):
    missing = [str(p) for p in map(Path, paths) if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise RuntimeError('Missing or empty required input(s): ' + ', '.join(missing))


def check_human_index(cfg, reference=None):
    default = reference is None
    fasta = Path(reference).expanduser() if reference else Path(cfg['human_reference_dir']) / 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa'
    files = [fasta, *[Path(str(fasta) + suffix) for suffix in ('.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac')]]
    if default:
        files.append(fasta.parent / 'BUILD.COMPLETE')
    require_files(files)


def check_stage(root, stage):
    from hierarchical_split import qc_folder
    path = qc_folder(root, stage) / 'summary.json'
    require_files([path])
    record = json.loads(path.read_text())
    if record.get('status') != 'complete':
        raise RuntimeError('Required stage is not complete: ' + str(path))
    for name, size in record.get('artifacts', {}).items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise RuntimeError('Invalid artifact path in ' + str(path))
        artifact = root / relative
        if not artifact.is_file() or artifact.stat().st_size != size:
            raise RuntimeError('Missing or changed required input: ' + str(artifact))


def check_hierarchy_stages(cfg, actions, samples, stage='both'):
    """Require earlier outputs only when this invocation will not produce them."""
    from hierarchical_split import qc_folder
    for sample in samples:
        root = Path(cfg['output_dir']) / sample['sample_id']
        if 'species' in actions and 'group' not in actions:
            check_stage(root, 'group')
        if 'evaluate' in actions:
            requested = [s for s in ('group', 'species') if stage in ('both', s)]
            available = [s for s in requested if s in actions or (qc_folder(root, s) / 'summary.json').is_file()]
            if not available:
                raise RuntimeError('No completed requested hierarchy stage for ' + str(root) + '; run split first')
            for existing in available:
                if existing not in actions:
                    check_stage(root, existing)


def check_submission_inputs(module, cfg, argv):
    if module == 'hierarchical_split':
        from hierarchical_split import parser
        from query import discover_samples
        from partitioned_build import database_snapshot
        from group_library import compact_snapshot
        args, _ = parser().parse_known_args(argv)
        if args.action == ['report']:
            return
        database_snapshot(cfg)
        if 'group' in args.action and (args.routing or cfg.get('group_routing', 'compact')) == 'compact':
            from ordinal_query import group_query_config
            group_query_config(cfg, compact_snapshot(cfg))
        samples = discover_samples(cfg, datasets=args.dat_list, only=args.sample, probe=False)
        if not samples:
            raise RuntimeError('No complete selected input sample pairs')
        require_files([s[k] for s in samples for k in ('r1', 'r2')])
        check_hierarchy_stages(cfg, args.action, samples, args.stage)
    elif module == 'group_library':
        from partitioned_build import database_snapshot
        from group_library import source_catalog
        database_snapshot(cfg)
        source_catalog(cfg)
    elif module in ('human_validate', 'validation_suite'):
        parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
        parser.add_argument('--reference')
        parser.add_argument('--sample')
        parser.add_argument('--output-dir')
        parser.add_argument('--target-fastq', nargs='+')
        parser.add_argument('--summary')
        parser.add_argument('--assignments')
        args, _ = parser.parse_known_args(argv)
        check_human_index(cfg, args.reference)
        require_files([p for p in [args.summary, args.assignments, *(args.target_fastq or [])] if p])
        if module == 'validation_suite':
            from qc_paths import configured_split_root
            from storage_layout import metadata_dir
            from validation_suite import sample_parts
            if not args.sample:
                raise ValueError('--sample is required')
            meta = metadata_dir(configured_split_root(cfg).joinpath(*sample_parts(args.sample)))
            require_files([meta / 'complete.json', meta / 'summary.json'])
        elif args.sample or args.output_dir:
            from qc_paths import configured_split_root
            from storage_layout import metadata_dir
            from human_validate import _safe_sample_id
            root = Path(args.output_dir) if args.output_dir else configured_split_root(cfg) / _safe_sample_id(args.sample)
            if not root.is_dir():
                raise RuntimeError('Missing required split directory: ' + str(root))
            if not args.target_fastq:
                require_files([metadata_dir(root) / 'summary.json'])
