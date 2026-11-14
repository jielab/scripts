#!/usr/bin/env python3
"""Move flat IBDmix exports into category/unit directories, preserving native links.

Dry-run by default. Preflight every cohort before --apply; never overwrite files.
A per-cohort JSON receipt records all renames and link changes for recovery.
"""
import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import re
import stat
from ibdmix_layout import destination_for


def normalize(path):
    return Path(re.sub(r'chrX\.nonpar-v2\.mac[0-9]+(?=/|\.|$)', 'chrX', str(path)))


def identity(path):
    s=path.stat()
    return [s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns]


def inventory(root):
    files,links,dirs={}, {}, []
    for directory,subdirs,names in os.walk(root,followlinks=False):
        dirs.append(Path(directory))
        for name in subdirs[:]:
            p=Path(directory)/name
            if p.is_symlink():
                raise ValueError(f'Unexpected directory symlink: {p}')
        for name in names:
            p=Path(directory)/name
            if p.is_symlink():links[p]=os.readlink(p)
            elif p.is_file():files[p]=None
            else:raise ValueError(f'Unsupported file: {p}')
    return files,links,dirs


def plan(root):
    files,links,dirs=inventory(root)
    mapping={p:root/normalize(p.relative_to(root)) for p in files}
    native=root/'tmp/ibdmix'
    # Native links unambiguously recover dots and directory components flattened
    # by the old publisher, including reference names such as Altai.2013.
    for p,link in links.items():
        if not p.is_relative_to(native):continue
        rel=p.relative_to(native)
        if len(rel.parts)<2:continue
        old=Path(os.path.abspath(p.parent/link))
        if old in files:
            mapping[old]=root/destination_for(Path(*rel.parts[1:]).as_posix(),rel.parts[0])
    for p in files:
        if p.parent!=root:continue
        if p.name.endswith('.tracts.xlsx'):
            mapping[p]=root/destination_for('ibdmix.tracts.xlsx',p.name[:-len('.tracts.xlsx')])
        elif p.name.endswith('.sample'):
            mapping[p]=root/destination_for('samples/ALL.txt',p.name[:-len('.sample')])
        elif p.name.startswith('chr') and mapping[p]==p:
            raise ValueError(f'Unmapped flat export: {p}')
    moves=[];destinations={}
    for old,new in mapping.items():
        if new in destinations:raise ValueError(f'Collision: {old}, {destinations[new]} -> {new}')
        destinations[new]=old
        if old!=new:
            if os.path.lexists(new):raise FileExistsError(new)
            moves.append(dict(old=str(old.relative_to(root)),new=str(new.relative_to(root)),identity=identity(old)))
    rewrites=[]
    for p,link in links.items():
        old_target=Path(os.path.abspath(p.parent/link))
        if not p.exists():raise ValueError(f'Existing broken link: {p} -> {link}')
        new_path=root/normalize(p.relative_to(root))
        target=mapping.get(old_target,root/normalize(old_target.relative_to(root)) if old_target.is_relative_to(root) else old_target)
        new_link=os.path.relpath(target,new_path.parent)
        if new_path!=p and os.path.lexists(new_path):raise FileExistsError(new_path)
        if new_path!=p or new_link!=link:
            rewrites.append(dict(old=str(p.relative_to(root)),new=str(new_path.relative_to(root)),before=link,after=new_link))
    additions=[]
    normalized_links={normalize(p.relative_to(root)) for p in links}
    # Legacy archives omitted workbooks from their native read view. Restore the
    # same lightweight links that new runs produce, without duplicating data.
    for row in moves:
        dest=Path(row['new'])
        if dest.parts[0]=='workbooks':
            link=Path('tmp/ibdmix')/dest.name.removesuffix('.tracts.xlsx')/'ibdmix.tracts.xlsx'
            if not os.path.lexists(root/link) and link not in normalized_links:
                additions.append(dict(path=str(link),target=os.path.relpath(root/dest,(root/link).parent)))
    return dict(root=str(root),moves=moves,links=rewrites,additions=additions,
                regular_files=len(files),symlinks=len(links),directories=[str(p.relative_to(root)) for p in dirs])



def promote_context(root):
    """Keep the cohort context as one real top-level file, without legacy storage."""
    root=Path(root); context=root/'gu-target.env'; legacy=root/'jobs/legacy'
    if not legacy.exists():
        if context.exists():
            assert context.is_file() and not context.is_symlink(),context
        return False
    source=legacy/'gu-target.env'
    assert sorted(p.name for p in legacy.iterdir())==['gu-target.env'],legacy
    assert context.is_symlink() and context.resolve()==source.resolve(),context
    before=identity(source)
    os.replace(source,context)
    assert not context.is_symlink() and identity(context)==before,context
    legacy.rmdir()
    return True


def apply(p):
    root=Path(p['root']); receipt=root/'jobs/layout-v2.migration.json'
    if receipt.exists():raise FileExistsError(receipt)
    receipt.parent.mkdir(parents=True,exist_ok=True)
    receipt.write_text(json.dumps(dict(p,status='planned'),indent=2)+'\n')
    for row in p['moves']:
        old,new=root/row['old'],root/row['new']
        assert identity(old)==row['identity'],old
        assert not os.path.lexists(new),new
        new.parent.mkdir(parents=True,exist_ok=True)
        os.rename(old,new)
    for row in p['links']:
        old,new=root/row['old'],root/row['new']
        assert os.readlink(old)==row['before'],old
        new.parent.mkdir(parents=True,exist_ok=True)
        tmp=new.with_name(new.name+'.layout-link')
        assert not os.path.lexists(tmp),tmp
        tmp.symlink_to(row['after']);os.replace(tmp,new)
        if old!=new:old.unlink()
    for row in p['additions']:
        dest=root/row['path'];dest.parent.mkdir(parents=True,exist_ok=True)
        dest.symlink_to(row['target'])
    for rel in sorted(p['directories'],key=lambda x:len(Path(x).parts),reverse=True):
        path=root/rel
        if path!=root and path.is_dir() and not any(path.iterdir()):path.rmdir()
    for row in p['moves']:
        assert identity(root/row['new'])==row['identity'],row
        assert not os.path.lexists(root/row['old']),row
    files,links,_=inventory(root)
    assert len(files)==p['regular_files']+1,(root,len(files),p['regular_files'])
    assert len(links)==p['symlinks']+len(p['additions'])
    for link in links:assert link.exists(),link
    assert not any(q.name.startswith('chr') for q in root.iterdir()),root
    assert not any('nonpar-v2.mac' in str(q.relative_to(root)) for q in list(files)+list(links))
    context_promoted=promote_context(root)
    receipt.write_text(json.dumps(dict(p,status='verified',context_promoted=context_promoted),indent=2)+'\n')
    return f'{root.name}: verified {len(p["moves"])} renames, {len(p["links"])} link updates; all data files retain inode, size and mtime'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cohorts',nargs='+',type=Path)
    parser.add_argument('--apply',action='store_true')
    a=parser.parse_args()
    roots=[p.resolve() for p in a.cohorts]
    assert len(set(roots))==len(roots)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        plans=[]
        for p in pool.map(plan,roots):
            plans.append(p)
            print(f'{Path(p["root"]).name}: plan {len(p["moves"])} renames, {len(p["links"])} link updates',flush=True)
        if a.apply:
            for result in pool.map(apply,plans):print(result,flush=True)

if __name__=='__main__':main()
