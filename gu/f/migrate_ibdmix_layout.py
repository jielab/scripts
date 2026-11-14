#!/usr/bin/env python3
"""Verify legacy archives, publish plain files, then move each verified cohort."""
import concurrent.futures
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import tarfile
from ibdmix_layout import publish_run_layout, destination_for

OLD=Path('/work/sph-huangj/analysis/gu-OLD/ibdmix')
NEW=Path('/work/sph-huangj/analysis/gu/ibdmix')

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def migrate(i):
    name=f'ukb{i:03d}'; source=OLD/name; target=NEW/name
    if target.exists():
        if (target/'jobs/migration.complete.json').exists() and not source.exists(): return name+' already complete'
        raise FileExistsError(target)
    stage=NEW/('.migrate-'+name)
    stage.mkdir(parents=True,exist_ok=True)
    records=[]
    archives=sorted(source.glob('chr*.raw.tar.gz'))
    if len(archives)!=23: raise ValueError(f'{name}: expected 23 archives, got {len(archives)}')
    for archive in archives:
        unit=re.sub(r'^chrX\.nonpar-v2\.mac[0-9]+$', 'chrX', archive.name.removesuffix('.raw.tar.gz'))
        run=stage/'tmp/ibdmix'/unit
        # A failed staging cohort can safely restart from its unchanged archives.
        if run.exists(): shutil.rmtree(run)
        run.mkdir(parents=True)
        with tarfile.open(archive,'r:gz') as tf:
            manifest=json.load(tf.extractfile('GU-MANIFEST.json'))
            expected={r['path']:r for r in manifest['files']}
            observed=set()
            for member in tf:
                p=Path(member.name)
                if p.is_absolute() or '..' in p.parts: raise ValueError(member.name)
                if member.isdir(): continue
                if not member.isfile(): raise ValueError(f'Unsupported archive link: {archive}:{member.name}')
                out=run/p;out.parent.mkdir(parents=True,exist_ok=True)
                with tf.extractfile(member) as inp, out.open('wb') as dest: shutil.copyfileobj(inp,dest,8*1024*1024)
                if member.name!='GU-MANIFEST.json':
                    row=expected[member.name]
                    if out.stat().st_size!=row['bytes'] or digest(out)!=row['sha256']: raise ValueError(f'Checksum: {out}')
                    observed.add(member.name)
            if observed!=set(expected): raise ValueError(f'Missing members: {archive}')
        if not publish_run_layout(run,stage,replace=True): raise ValueError(run)
        for rel,row in expected.items():
            out=stage/destination_for(rel,unit)
            if digest(out)!=row['sha256']: raise ValueError(f'Published checksum: {out}')
        records.append(dict(archive=archive.name,files=len(expected)))
    for p in source.rglob('*'):
        if p.is_dir() or p.name.endswith('.raw.tar.gz'): continue
        rel=p.relative_to(source)
        if len(rel.parts)==1 and p.name.endswith('.run.meta.tsv'):
            unit=re.sub(r'^chrX\.nonpar-v2\.mac[0-9]+$', 'chrX', p.name.removesuffix('.run.meta.tsv')); dest=stage/'jobs'/unit/'run.meta.tsv'
        elif len(rel.parts)==1 and p.name.endswith('.tracts.xlsx'):
            dest=stage/destination_for('ibdmix.tracts.xlsx',p.name.removesuffix('.tracts.xlsx'))
        elif len(rel.parts)==1 and p.suffix=='.sample':
            dest=stage/destination_for('samples/ALL.txt',p.name.removesuffix('.sample'))
        elif rel.as_posix()=='gu-target.env': dest=stage/rel
        elif rel.parts[0] in ('logs','jobs','qc','tmp'): dest=stage/rel
        else: dest=stage/'jobs'/rel
        dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists():
            if digest(p)!=digest(dest): raise ValueError(f'Conflicting export: {p}')
        else: os.link(p,dest)
    receipt=stage/'jobs/migration.complete.json'
    receipt.write_text(json.dumps(dict(source=str(source),archives=records),indent=2)+'\n')
    os.rename(stage,target)
    shutil.rmtree(source)
    return name+' verified and moved'

if __name__=='__main__':
    NEW.mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for result in pool.map(migrate,range(3,50)): print(result,flush=True)
