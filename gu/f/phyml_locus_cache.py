#!/usr/bin/env python3
"""Validate whole-locus completion before launching a GU worker.

Large reference files use path/size/mtime fingerprints; small results use SHA256.
Historical runs require an explicit successful worker footer, never just TSVs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def command(path):
    env, argv = {}, []
    for line in path.read_text().splitlines():
        words = shlex.split(line)
        if words[:1] == ['export']:
            for word in words[1:]:
                key, value = word.split('=', 1)
                env[key] = value
        elif words[:1] == ['exec']:
            argv = words[1:]
    if len(argv) < 2 or argv[1] != 'phyml':
        raise ValueError('not a phyml worker')
    args = argv[2:]
    if args[:1] == ['run']:
        args = args[1:]
    if len(args) % 2 or any(not k.startswith('--') for k in args[::2]):
        raise ValueError('unsupported worker arguments')
    opts = dict(zip(args[::2], args[1::2]))
    return env, opts


def request(cmd, archaic_root):
    env, opts = command(cmd)
    bed = Path(opts['--loci']).read_text().split()
    if len(bed) != 4:
        raise ValueError('expected one locus')
    lead = [r for r in rows(Path(env['GU_PHYML_LEAD_TABLE'])) if r['locus_id'] == bed[3]]
    if len(lead) != 1:
        raise ValueError('missing original lead')
    prefix = Path(opts['--target-dir'])
    ch = bed[0].removeprefix('chr')
    sources = [p for p in prefix.parent.glob(prefix.name + ch + '.*')
               if p.name.endswith(('.pgen', '.pvar', '.pvar.zst', '.psam', '.bed', '.bim', '.fam', '.vcf.gz', '.bcf', '.tbi', '.csi'))]
    if not sources:
        raise ValueError('missing target input')
    panel = opts.get('--sample-panel') or os.environ.get('GU_SAMPLE_PANEL')
    if not panel:
        panel = next((str(p) for p in (prefix.parent/'samples.txt', prefix.parent.parent/'samples.txt') if p.is_file()), None)
    if not panel:
        raise ValueError('sample panel unavailable')
    sources.append(Path(panel))
    root = Path(archaic_root)
    if not root.is_dir():
        raise ValueError('archaic root unavailable')
    sources += [p for p in root.rglob('*') if re.search(r'(?<![0-9])(?:chr)?' + re.escape(ch) + r'[_.]', p.name) and p.is_file()]
    stamps = {}
    for p in sorted(set(sources)):
        if p.is_file():
            s = p.stat()
            stamps[str(p.resolve())] = [s.st_size, s.st_mtime_ns]
    ignored = {'--memory-cap', '--foreground', '--auto-final', '--replace-phyml', '--phyml-jobs', '--loci', '--sample-panel'}
    code = {name: digest(Path(__file__).with_name(name)) for name in
            ('phyml_gwas.py', 'phyml_core.py', 'phyml_thresholds.py', 'phyml_tree_summary.py', 'phyml_panel_b.R', 'phyml_run.py')}
    return dict(schema=1, lead=lead[0], bed=bed, options={k:v for k,v in opts.items() if k not in ignored},
                sources=stamps, sample_panel=str(Path(panel).resolve()), archaic_root=str(root.resolve()), code=code,
                environment={k:os.environ.get(k, '') for k in ('GU_CHRX_MALE_ONLY', 'GU_CHRX_PAR_DIPLOID')})


def outputs(out):
    required = ('gwas_loci.tsv', 'gwas_haplotypes.tsv', 'gwas_copies.tsv', 'gwas_lead.tsv',
                'gwas_parameters.json', 'loci.tsv', 'trees.tsv', 'evidence_trees.tsv',
                'haplotypes.tsv', 'haplotype_samples.tsv', 'skipped_loci.tsv')
    if any(not (out/'final'/name).is_file() or not (out/'final'/name).stat().st_size for name in required):
        raise ValueError('missing final outputs')
    return {str(p.relative_to(out)): digest(p) for folder in ('final', 'loci')
            for p in sorted((out/folder).glob('*')) if p.is_file()
            and not p.name.endswith('.lock') and '.failed.' not in p.name}


def successful(cmd, req):
    out = cmd.parent
    log = cmd.with_suffix('.log').read_text()
    if cmd.with_suffix('.err').exists() or f'analysis_unit={cmd.stem} status=complete' not in log:
        raise ValueError('no successful worker completion')
    if 'status=failed' in log or 'plot export failed' in log:
        raise ValueError('failed worker')
    if f'archaic reference root={req["archaic_root"]}\n' not in log:
        raise ValueError('reference root changed')
    if f'sample metadata={req["sample_panel"]}\n' not in log:
        raise ValueError('sample panel changed')
    if rows(out/'final/gwas_lead.tsv') != [req['lead']]:
        raise ValueError('lead changed')
    params = json.loads((out/'final/gwas_parameters.json').read_text())
    if params['workflow'] != 'gwas_lead_ld_core_archaic5_v2':
        raise ValueError('historical workflow changed')
    summaries = rows(out/'final/gwas_loci.tsv')
    if len(summaries) != 2 or any(r['status'] in ('tree_failed', 'not_evaluable') for r in summaries):
        raise ValueError('incomplete locus')
    trees = rows(out/'final/trees.tsv')
    if len(trees) != 1 or trees[0]['tree_status'] not in ('complete', 'not_run'):
        raise ValueError('incomplete tree summary')
    if any(r['status'] in ('tree_supported', 'tree_not_supported') for r in summaries) and trees[0]['tree_status'] != 'complete':
        raise ValueError('tree result without completed tree')
    for tree in trees:
        if tree['tree_status'] == 'complete':
            from phyml_run import completion_error
            if completion_error(out/'loci/haplotypes.phy', 100):
                raise ValueError('incomplete tree')
            for lineage in ('Neanderthal', 'Denisovan'):
                for suffix in ('.png', '.pdf', '.full.png', '.full.pdf'):
                    plot = out/'loci'/f'haplotypes.phy_phyml_tree.{lineage}.panelB{suffix}'
                    if not plot.is_file() or not plot.stat().st_size:
                        raise ValueError('missing plot')


def process(mode, cmd, archaic_root):
    receipt = cmd.parent/'.phyml.locus.complete.json'
    if mode == 'adopt' and receipt.exists():
        return False
    req = request(cmd, archaic_root)
    if mode == 'check':
        if command(cmd)[1].get('--replace-phyml', 'FALSE') == 'TRUE':
            return False
        data = json.loads(receipt.read_text())
        return data['request'] == req and data['outputs'] == outputs(cmd.parent)
    successful(cmd, req)
    if mode == 'adopt':
        # Old runs have no source fingerprints. Only adopt sources older than
        # their successful log, with the original command and exact saved lead.
        finished = cmd.with_suffix('.log').stat().st_mtime_ns
        if any(stamp[1] > finished for stamp in req['sources'].values()):
            return False
        if cmd.stat().st_mtime_ns > finished:
            return False
    data = dict(request=req, outputs=outputs(cmd.parent))
    tmp = receipt.with_name(receipt.name + f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(receipt)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('check', 'seal', 'adopt'))
    p.add_argument('cmd', type=Path)
    p.add_argument('--archaic-root', required=True)
    a = p.parse_args()
    try:
        ok = process(a.mode, a.cmd, a.archaic_root)
    except (OSError, ValueError, KeyError, TypeError):
        ok = False
    raise SystemExit(0 if ok else 1)


if __name__ == '__main__':
    main()
