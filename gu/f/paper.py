#!/usr/bin/env python3
"""Generate PhyML summaries, figures and audits from completed GU results."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / '0f'))
from workflow_paths import configured_path


def main():
    parser = argparse.ArgumentParser(prog='gu.paper.sh', description=__doc__)
    modules = parser.add_subparsers(dest='module', required=True)
    report = modules.add_parser('phyml', help=__doc__)
    for name, help_text in (
        ('normalize', 'Normalized results containing phyml/'),
        ('database', 'Existing results database for carrier validation'),
        ('dataset-dir', 'Saved dataset containing chr*/phyml.raw.tar.gz'),
        ('lead-table', 'Original lead table for the audit'),
        ('dir-phyml', 'Completed locus run or parent containing locus runs')):
        report.add_argument('--'+name, required=True, type=Path, help=help_text)
    report.add_argument('--dir-out', '--output', dest='output', type=Path,
        default=Path(configured_path('ANALYSIS_ROOT', 'gu', 'gu.paper', 'phyml')))
    args = parser.parse_args()
    for path in (args.database, args.lead_table):
        if not path.is_file(): parser.error('Missing input file: '+str(path))
    if not list((args.normalize/'phyml').glob('**/final/gwas_loci.tsv')):
        parser.error('--normalize contains no PhyML locus summaries')
    if not list(args.dataset_dir.glob('chr*/phyml.raw.tar.gz')):
        parser.error('--dataset-dir contains no saved PhyML archives')
    runs = sorted({p.parent.parent for p in args.dir_phyml.glob('**/final/evidence_trees.tsv')})
    if not runs or any(not (run/'loci').is_dir() for run in runs):
        parser.error('--dir-phyml must contain completed runs with final/evidence_trees.tsv and loci/')
    output = args.output.resolve()
    inputs = [args.normalize, args.dataset_dir, args.dir_phyml]
    if any(path.resolve().is_relative_to(output) or output.is_relative_to(path.resolve()) for path in inputs):
        parser.error('--dir-out must be separate from the input directories')
    engine = [sys.executable, str(ROOT/'phyml.py')]
    rscript = os.environ.get('WORKFLOW_RSCRIPT') or str(ROOT.parents[1]/'0f/0os.r.sh')
    commands = [('summary', engine+['report', '--normalize', str(args.normalize),
        '--database', str(args.database), '--output', str(output/'summary')])]
    for run in runs:
        relative = run.relative_to(args.dir_phyml)
        commands.append(('figures/'+str(relative), [rscript, '--vanilla', str(ROOT/'phyml.R'),
            '--out', str(run), '--dir-out', str(output/'figures'/relative)]))
    commands.append(('audit', engine+['audit', '--dataset-dir', str(args.dataset_dir),
        '--lead-table', str(args.lead_table), '--output', str(output/'audit')]))
    output.mkdir(parents=True, exist_ok=True)
    record = dict(status='running', stages=[])
    def save():
        temporary = output/'report.json.partial'
        temporary.write_text(json.dumps(record, indent=2)+'\n')
        temporary.replace(output/'report.json')
    save()
    for stage, command in commands:
        print('[GU paper] '+stage, flush=True)
        try:
            result = subprocess.run(command, check=False)
        except OSError as exc:
            record.update(status='failed', error=str(exc), failed_stage=stage); save()
            raise
        record['stages'].append(dict(stage=stage, returncode=result.returncode))
        if result.returncode:
            record.update(status='failed', failed_stage=stage); save()
            return result.returncode
        save()
    record['status'] = 'complete'; save()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
