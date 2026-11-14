#!/usr/bin/env python3
"""Run or inspect the current gnomAD scan format; no output migration interface."""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '0f'))
from workflow_paths import configured_path
if 'allele_common' not in sys.modules:
    spec = importlib.util.spec_from_file_location('allele_common', Path(__file__).with_name('0.common.py'))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
from allele_common import fingerprint

SCHEMA = 'allele3-scan-store-final'
DEFAULT_STORE = str(Path(os.environ.get('ALLELE_ANALYSIS_ROOT', configured_path('ANALYSIS_ROOT', 'allele3'))) / 'gnomad_scan')
MASK_OPTIONS = {'lcr': 'architecture_lcr_bed', 'segdup': 'architecture_segdup_bed',
                'low_mappability': 'architecture_low_mappability_bed'}


def absolute(path, root=None):
    from analysis_integrity import resolve_result_input
    value = Path(path).expanduser()
    return Path(resolve_result_input(value if value.is_absolute() else Path(root or Path.cwd()) / value)).resolve()


def load_run_metadata(out):
    root = absolute(out)
    run = json.loads((root / 'run.json').read_text())
    if run.get('status') not in {'complete', 'scanned'}:
        raise ValueError('Scan is incomplete: ' + str(root))
    if not run.get('signature') or not isinstance(run.get('shards'), list) or not run['shards']:
        raise ValueError('Scan manifest lacks its shard list: ' + str(root))
    if len(set(run['shards'])) != len(run['shards']):
        raise ValueError('Duplicate scan shards')
    run['_store'] = {'root': str(root)}
    run['shards'] = [str(resolve_run_path(run, name)) for name in run['shards']]
    return run


def resolve_run_path(run, path):
    value = Path(path).expanduser()
    return value if value.is_absolute() else Path(run.get('_store', {}).get('root', '.')) / value


def verify_recorded_input(run, expected):
    if not isinstance(expected, dict) or 'path' not in expected:
        return None
    return verify_cached_fingerprint(resolve_run_path(run, expected['path']), expected)


def verify_cached_fingerprint(path, expected, *, context=None):
    path = Path(path)
    if not path.is_file() or path.stat().st_size != expected.get('size'):
        raise ValueError('Required scan file is missing or incomplete: ' + str(path))
    # No content read, timestamp gate, or old-output relocation fallback.
    return dict(path=str(path.resolve()), size=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns)


def verify_shard(shard, signature):
    root = Path(shard)
    complete = json.loads((root / 'COMPLETE.json').read_text())
    if complete.get('signature') != signature or not complete.get('outputs'):
        raise ValueError('Incomplete or mismatched scan shard: ' + str(root))
    for name, expected in complete['outputs'].items():
        if Path(name).name != name:
            raise ValueError('Invalid shard output name: ' + name)
        verify_cached_fingerprint(root / name, expected)
    return complete


def verify(store):
    root = absolute(store)
    run = load_run_metadata(root)
    for shard in run['shards']:
        verify_shard(shard, run['signature'])
    return dict(schema=SCHEMA, status='verified', store=str(root), shards=len(run['shards']),
                vcf_reads=0, validation='completion markers and file sizes')


def scan_command(store, extra):
    """Current scan profile; generated files belong to the selected output."""
    extra = [part for word in extra for part in (word.split('=', 1) if word.startswith('--') and '=' in word else [word])]
    def has(*names):
        return any(word.split("=", 1)[0] in names for word in extra)
    command = [sys.executable, str(Path(__file__).with_name("scan.py")), "run", "--out", str(store)]
    defaults = [("--fasta", os.environ.get("ALLELE_FASTA", configured_path('DIR_GEN', '1hgp/38/GRCH38.fasta.gz'))),
                ("--jobs", os.environ.get("ALLELE_JOBS", "2")), ("--min-an", "100"), ("--min-af", "0"),
                ("--thresholds", "1,2,5,10"), ("--backend", "auto"), ("--candidate-storage", "deferred"),
                ("--ld-block-dir", os.environ.get("ALLELE_LD_BLOCK_DIR", configured_path('DIR_LD', 'block'))),
                ("--architecture-build", "38"),
                ("--erv-joined", os.environ.get("ALLELE_ERV_JOINED", configured_path('DIR_ANNOT', 'allele3/repeats/hg38/rmskJoinedCurrent.txt.gz'))),
                ("--architecture-annotation-cache", str(store / "annotations")),
                ("--architecture-low-mappability-bed", os.environ.get("ALLELE_LOW_MAPPABILITY_BED", configured_path('DIR_FILES', 'allele3/umap/hg38.k100.lt1.bed.gz')))]
    if not has("--vcf", "--dir-vcf", "--vcf-dir"):
        command += ["--dir-vcf", os.environ.get("ALLELE_VCF_DIR", configured_path('DIR_ANNOT', 'allele3/gnomad/sites'))]
    if has("--dir-cache", "--cache-dir"):
        raise ValueError("Scan shards are owned by --out; external shared scan caches are not supported")
    for option, value in defaults:
        if not has(option): command += [option, value]
    if not has("--with-architecture"): command += ["--with-architecture"]
    if has("--out", "--dir-out"):
        raise ValueError("Set the scan output once with --out")
    return command + list(extra)



def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('scan', 'verify', 'scan,verify'), nargs='?', default='scan')
    parser.add_argument('--out', dest='store', default=DEFAULT_STORE)
    parser.add_argument('--replace', choices=('TRUE', 'FALSE'), default=os.environ.get('ALLELE_REPLACE', 'TRUE'))
    parser.add_argument('--dry-run', action='store_true')
    args, extra = parser.parse_known_args(argv)
    if extra and args.command not in ('scan', 'scan,verify'):
        parser.error('Scan arguments are valid only for scan')
    store = absolute(args.store)
    if args.command == 'verify':
        print(json.dumps(verify(store), indent=2))
        return 0
    command = scan_command(store, extra)
    if args.dry_run:
        print(json.dumps(dict(status='planned', out=str(store), argv=command), indent=2))
        return 0
    if args.replace == 'FALSE':
        try:
            result = verify(store)
            print(json.dumps(dict(result, status='reused'), indent=2))
            return 0
        except (OSError, ValueError, KeyError, TypeError):
            pass  # Missing/partial output is rebuilt, never a dead-end checkpoint.
    from analysis_integrity import output_run, publish_scan
    input_options = {'--dir-vcf','--vcf-dir','--vcf','--fasta','--erv-joined','--ld-block-dir','--architecture-low-mappability-bed'}
    protected = [command[i+1] for i,x in enumerate(command[:-1]) if x in input_options]
    with output_run(store, True, protected):
        code = subprocess.run(command, check=False).returncode
        if code == 0:
            verify(store)
            publish_scan(store)
        return code


def resolve_mask_options(run, overrides=None, *, verify=False):
    """Explicit mask relocation retains scientific mask identity checks only."""
    options = dict(run.get('options', run.get('identity', {}).get('options', {})))
    records = (options.get('architecture_annotations') or {}).get('files', {})
    for name, destination in (overrides or {}).items():
        if not destination:
            continue
        if name not in MASK_OPTIONS:
            raise ValueError('Unknown mask: ' + name)
        key = MASK_OPTIONS[name]
        expected = records.get(key)
        if not options.get(key) or not isinstance(expected, dict):
            raise ValueError('Cannot add an unrecorded mask to an existing scan')
        path = Path(destination).expanduser().resolve()
        if verify:
            actual = fingerprint(path)
            # Bounded edge sampling for an explicitly moved scientific mask;
            # ordinary completed-output reuse never reads file contents.
            if any(actual.get(k) != expected.get(k) for k in ('size', 'edge_sha256')):
                raise ValueError('Explicit mask differs from the mask used by the scan: ' + str(path))
        options[key] = str(path)
    return options


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError) as exc:
        raise SystemExit('scan-store: ' + str(exc))
