#!/usr/bin/env python3
"""Create a separate projection basis using unchanged axes and unique coordinates."""
import argparse
from collections import Counter
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pca_project as pp


def subset_table(source, destination, keep):
    """Keep exact original header/row bytes, including frequency orientation."""
    with source.open('rb') as src, destination.open('wb') as dest:
        header = src.readline()
        index = header.decode().lstrip('#').split().index('ID')
        dest.write(header)
        for line in src:
            if line.decode().split()[index] in keep:
                dest.write(line)


def prepare(args):
    source, out = args.basis_dir.resolve(), args.out.resolve()
    if source == out or source in out.parents or out in source.parents:
        raise ValueError('New basis must be separate from the original basis')
    if args.threads < 1 or args.memory_mb < 640:
        raise ValueError('Positive threads and at least 640 MiB memory required')
    # Only this derivation step may read duplicate coordinates. IDs, alleles,
    # frequencies, loadings, receipts and sample order remain strictly checked.
    original = pp._basis(source, allow_duplicate_coordinates=True)
    executable = shutil.which(str(args.plink2))
    if not executable:
        raise ValueError('PLINK2 executable unavailable')
    tool = dict(executable=pp.state(executable), version=subprocess.check_output([executable,'--version'],text=True).strip())
    coordinate_counts = Counter((m['chrom'],m['pos']) for m in original['markers'].values())
    excluded = [m for m in original['markers'].values() if coordinate_counts[m['chrom'],m['pos']] > 1]
    keep = [m['id'] for m in original['markers'].values() if coordinate_counts[m['chrom'],m['pos']] == 1]
    if not keep:
        raise ValueError('No unique-coordinate reference markers remain')
    identity = dict(source=original['inputs'], tool=tool, policy='exclude_all_records_at_repeated_coordinates',
                    version=1, code_sha256=pp.sha(__file__))
    signature = hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    out.mkdir(parents=True,exist_ok=True)
    with (out/'.basis.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        receipt = out/'derivation.json'
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if saved['signature'] != signature:
                raise ValueError('Basis derivation inputs/code/tool changed; use a new directory')
            if saved['status'] == 'complete':
                for name,digest in saved['outputs_sha256'].items():
                    if pp.sha(out/name) != digest:
                        raise ValueError('Derived basis file changed: '+name)
                pp.validate_basis(out)
                return saved
        elif any(p.name != '.basis.lock' for p in out.iterdir()):
            raise ValueError('New basis requires an empty output directory')
        record = dict(status='preparing',signature=signature,identity=identity,
                      original_basis=str(source),original_signature=original['plan']['signature'],
                      n_original_markers=len(original['markers']),n_retained_markers=len(keep),
                      n_excluded_markers=len(excluded),n_reference_samples=len(original['samples']),
                      axes_refitted=False,original_basis_modified=False)
        pp.save_json(receipt,record)
        pp.atomic_text(out/'retained_markers.ids',''.join(x+'\n' for x in keep))
        pp.write_table(out/'excluded_markers.tsv',excluded,['chrom','pos','id','ref','alt'])
        for name in ('ref.acount','ref.eigenvec.allele'):
            subset_table(source/name,out/name,set(keep))
        shutil.copy2(source/'ref.eigenval',out/'ref.eigenval')
        files=original['files']
        common=[executable,'--threads',str(args.threads),'--memory',str(args.memory_mb)]
        commands=[common+['--pgen',files['pgen'],'--pvar',files['pvar'],'--psam',files['psam'],
                  '--extract',str(out/'retained_markers.ids'),'--make-pgen','--out',str(out/'all_pruned')],
                  common+['--pfile',str(out/'all_pruned'),'--read-freq',str(out/'ref.acount'),
                  '--score',str(out/'ref.eigenvec.allele'),'2','5','header-read','variance-standardize','list-variants',
                  'cols=fid,nallele,denom,dosagesum,scoreavgs,scoresums','--score-col-nums',f"6-{5+len(original['pcs'])}",
                  '--out',str(out/'all')]]
        for label,command in zip(('extract','score'),commands):
            print('[projection basis] '+label,flush=True)
            with (out/(label+'.console.log')).open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        used=(out/'all.sscore.vars').read_text().splitlines()
        if set(used)!=set(keep) or len(used)!=len(keep):
            raise ValueError('PLINK did not score exactly the retained reference panel')
        for path,saved in original['inputs'].items():
            if pp.state(saved['path'])!=saved:
                raise ValueError('Original basis changed during derivation')
        plan=dict(schema_version='derived-frozen-projection-basis-1',signature=signature,build='GRCh37',
                  inputs=original['plan'].get('inputs',{}),settings=original['plan'].get('settings',{}),
                  counts={**original['plan'].get('counts',{}),'retained_projection_markers':len(keep)},
                  coordinate_kind='unchanged_reference_axes_rescored_on_unique_coordinate_subset',
                  source_basis=str(source),source_signature=original['plan']['signature'],axes_refitted=False)
        pp.save_json(out/'pca_plan.json',plan)
        expected=('all.sscore','ref.acount','ref.eigenvec.allele','ref.eigenval')
        pp.atomic_text(out/'completed.sha256',''.join(f'{pp.sha(out/name)}  {name}\n' for name in expected))
        pp.atomic_text(out/'PCA_COMPLETE',signature+'\n')
        validation=pp.validate_basis(out)
        geometry=pp.geometry_change(source/'all.sscore',out/'all.sscore',original['pcs'],out/'reference_geometry_change.tsv')
        record.update(status='complete',commands=commands,validation=validation,geometry_change=geometry,
                      outputs_sha256={p.name:pp.sha(p) for p in out.iterdir() if p.is_file() and p.name not in {'.basis.lock','derivation.json'}})
        pp.save_json(receipt,record)
        return record


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--basis-dir',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--plink2',default='plink2');p.add_argument('--threads',type=int,default=4)
    p.add_argument('--memory-mb',type=int,default=8192)
    args=p.parse_args(argv)
    try:
        result=prepare(args)
        print(json.dumps({k:result[k] for k in ('status','n_original_markers','n_retained_markers','n_excluded_markers','axes_refitted')},indent=2))
        return 0
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
        print('pca-basis: '+str(exc),file=sys.stderr)
        return 2


if __name__=='__main__':raise SystemExit(main())
