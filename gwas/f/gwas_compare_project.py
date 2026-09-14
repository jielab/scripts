#!/usr/bin/env python3
"""Resolve project GWAS files and create non-destructive build-aligned comparison inputs."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def build_of(file):
    trait = file.name.removesuffix('.gz').removesuffix('.thin')
    candidates = [Path(str(file) + '.grch'), file.parent.parent / 'qc' / (trait + '.grch')]
    builds = {p.read_text().strip() for p in candidates if p.is_file()}
    if len(builds) != 1 or not builds.issubset({'37', '38'}):
        raise ValueError('Missing/conflicting genome-build metadata for ' + str(file))
    return next(iter(builds))


def resolve(project, category, anchor):
    base = Path(project) / category
    files = [p / 'gwas' / (p.name + '.gz') for p in base.iterdir() if p.is_dir()]
    files = sorted((p for p in files if p.is_file()),
                   key=lambda p: [int(x) if x.isdigit() else x for x in re.split(r'(\d+)', p.parent.parent.name)])
    if anchor:
        first = [p for p in files if p.parent.parent.name == anchor]
        if not first:
            raise ValueError(f'Anchor {anchor} not found; finish processing that GWAS first')
        files = first + [p for p in files if p not in first]
    return files


def require_build(files, required):
    builds = {str(f): build_of(f) for f in files}
    wrong = [f'{Path(f).stem}=GRCh{b}' for f, b in builds.items() if b != required]
    if wrong:
        raise ValueError(f'All inputs must be GRCh{required}; finish gwas_format.sh liftover first: '
                         + ', '.join(wrong))
    return builds


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('compare', 'ldsc', 'shiny'))
    p.add_argument('--dir-gwas', '--project-dir', dest='project_dir')
    p.add_argument('--category', default='common')
    p.add_argument('--anchor')
    p.add_argument('--gwas-files')
    p.add_argument('--grch', choices=('37', '38'))
    p.add_argument('--require-grch', choices=('37', '38'),
                   help='Require every finalized input to have this build; do not lift comparison copies')
    p.add_argument('--dir-out', '--output-dir', dest='output_dir', help='Default: /mnt/d/analysis/gwas/<mode>')
    p.add_argument('--check-only', choices=('TRUE', 'FALSE'), default='FALSE')
    a, other = p.parse_known_args()
    if bool(a.project_dir) == bool(a.gwas_files):
        p.error('Specify one of --dir-gwas or --gwas-files')
    files = resolve(a.project_dir, a.category, a.anchor) if a.project_dir else [Path(x).resolve() for x in a.gwas_files.split(',')]
    minimum = 1 if a.mode == 'shiny' else 2
    if len(files) < minimum or any(not f.is_file() for f in files):
        p.error(f'{a.mode} requires at least {minimum} existing GWAS file(s)')
    if a.mode == 'shiny':
        files = [f.with_name(f.name[:-3] + '.thin.gz')
                 if not f.name.endswith('.thin.gz') and f.with_name(f.name[:-3] + '.thin.gz').is_file() else f
                 for f in files]
    labels = [f.name.removesuffix('.gz').removesuffix('.thin') if a.mode == 'shiny'
              else f.name.removesuffix('.gz') for f in files]
    checked_builds = {}
    if a.require_grch:
        if a.grch and a.grch != a.require_grch:
            p.error('--grch and --require-grch must agree')
        checked_builds = require_build(files, a.require_grch)
        a.grch = a.require_grch
    out = Path(a.output_dir or ('/mnt/d/analysis/gwas/' + a.mode)).resolve()
    if a.mode == 'shiny':
        a.grch = a.grch or '37'
        builds = {str(f): build_of(f) for f in files}
        if len(set(builds.values())) != 1:
            p.error('shiny inputs must share one source GRCh build; harmonize the inputs first')
        if len(set(labels)) != len(labels):
            p.error('shiny requires unique GWAS filenames')
        source_build = next(iter(builds.values()))
        out.mkdir(parents=True, exist_ok=True)
        config = out / 'shiny.inputs.json'
        config.write_text(json.dumps(dict(source_build=source_build,
            tracks=[dict(file=str(f.resolve()), trait=label) for f, label in zip(files, labels)]), indent=2))
        print(f'shiny: {len(files)} GWAS; source GRCh{source_build}; initial GRCh{a.grch or source_build}', flush=True)
        if a.check_only == 'TRUE':
            print(config)
            return
        os.execvp('Rscript', ['Rscript', '--vanilla', str(Path(__file__).resolve().parent.parent / 'shiny' / 'app.R'),
                        '--manifest', str(config), '--output-dir', str(out),
                        '--grch', a.grch or source_build] + other)
        return
    manifest = []
    aligned = []
    for f in files:
        build = checked_builds.get(str(f)) or (build_of(f) if a.mode == 'compare' else 'not_required_for_rsID_LDSC')
        target = f
        if a.mode == 'compare' and a.grch and build != a.grch:
            target = out / 'inputs' / ('grch' + a.grch) / f.name
            if a.check_only != 'TRUE':
                chain = {'37': 'hg19ToHg38.over.chain.gz', '38': 'hg38ToHg19.over.chain.gz'}[build]
                subprocess.run([sys.executable, str(Path(__file__).with_name('gwas_liftover.py')),
                                '--input', str(f), '--output', str(target),
                                '--qc-prefix', str(target.parent / 'qc' / f.stem),
                                '--source-build', build, '--target-build', a.grch,
                                '--chain', '/mnt/d/files/liftOver/' + chain,
                                '--liftOver', '/mnt/d/software/bin/liftOver'], check=True)
        aligned.append(target)
        manifest.append(dict(trait=f.stem, input=str(f), source_build=build, comparison_input=str(target)))
    if a.mode == 'compare' and not a.grch and len({x['source_build'] for x in manifest}) > 1:
        p.error('Mixed builds; supply --grch 38 for cached harmonization')
    out.mkdir(parents=True, exist_ok=True)
    import csv
    with (out / 'inputs.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(manifest[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(manifest)
    print(f'Discovered {len(files)} GWAS: ' + ','.join(labels), flush=True)
    if a.check_only == 'TRUE':
        print(out / 'inputs.tsv'); return
    args = ['--gwas-files', ','.join(map(str, aligned)), '--output-dir', str(out)] + other
    if a.mode == 'compare':
        if '--labels' not in other:
            args += ['--labels', ','.join(labels)]
        if a.grch:
            args += ['--grch', a.grch]
        cmd = ['Rscript', str(Path(__file__).with_name('gwas_compare.R'))] + args
    else:
        cmd = [sys.executable, str(Path(__file__).with_name('gwas_ldsc.py'))] + args
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
