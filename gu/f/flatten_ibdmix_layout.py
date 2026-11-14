#!/usr/bin/env python3
"""Flatten physical IBDmix exports while preserving native read links.

Preflight all cohorts before renaming anything. No file contents are changed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
from ibdmix_layout import destination_for

CATEGORIES={'raw','segments','lod','samples','workbooks','qc','logs','jobs'}

def identity(p):
    s=p.stat()
    return [s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns]

def inventory(root):
    files=set();links={};dirs=[]
    for directory,subdirs,names in os.walk(root):
        dirs.append(Path(directory))
        for name in subdirs:
            assert not (Path(directory)/name).is_symlink(),(directory,name)
        for name in names:
            p=Path(directory)/name
            if p.is_symlink():links[p]=os.readlink(p)
            else:files.add(p)
    return files,links,dirs

def plan(root):
    files,links,dirs=inventory(root);mapping={};additions=[]
    native=root/'tmp/ibdmix'
    for p,link in links.items():
        if not p.is_relative_to(native):continue
        rel=p.relative_to(native)
        assert len(rel.parts)>1,rel
        source=Path(os.path.abspath(p.parent/link))
        assert source in files,(p,source)
        dest=root/destination_for(Path(*rel.parts[1:]).as_posix(),rel.parts[0])
        assert source not in mapping or mapping[source]==dest,(source,dest,mapping.get(source))
        mapping[source]=dest
    for p in files:
        if p.is_relative_to(native):
            rel=p.relative_to(native)
            dest=root/destination_for(Path(*rel.parts[1:]).as_posix(),rel.parts[0])
            mapping[p]=dest
            additions.append(dict(path=str(p.relative_to(root)),target=os.path.relpath(dest,p.parent)))
        else:mapping.setdefault(p,p)
    destinations={};moves=[]
    for old,new in mapping.items():
        assert new not in destinations,(old,new,destinations.get(new))
        destinations[new]=old
        relative=new.relative_to(root)
        if relative.parts[0] in CATEGORIES:
            assert len(relative.parts)==2,('Unmapped nested output',old,new)
        if relative.parts[0]=='tmp':
            assert len(relative.parts)==2,('Unmapped intermediate',old,new)
        if old!=new:
            assert not os.path.lexists(new),new
            moves.append(dict(old=str(old.relative_to(root)),new=str(new.relative_to(root)),identity=identity(old)))
    rewrites=[]
    for p,link in links.items():
        target=Path(os.path.abspath(p.parent/link))
        assert p.exists(),('Broken original link',p)
        newtarget=mapping.get(target,target)
        value=os.path.relpath(newtarget,p.parent)
        if value!=link:rewrites.append(dict(path=str(p.relative_to(root)),before=link,after=value))
    return dict(root=str(root),moves=moves,links=rewrites,additions=additions,
                regular_files=len(files),symlinks=len(links),directories=[str(d.relative_to(root)) for d in dirs])

def apply(p):
    root=Path(p['root']);receipt=root/'jobs/layout-flat.migration.json'
    assert not receipt.exists(),receipt
    receipt.write_text(json.dumps(dict(p,status='planned'),indent=2)+'\n')
    for row in p['moves']:
        old,new=root/row['old'],root/row['new']
        assert identity(old)==row['identity'],old
        assert not os.path.lexists(new),new
        new.parent.mkdir(parents=True,exist_ok=True);os.rename(old,new)
    for row in p['links']:
        link=root/row['path'];assert os.readlink(link)==row['before']
        tmp=link.with_name(link.name+'.flat-link');assert not os.path.lexists(tmp)
        tmp.symlink_to(row['after']);os.replace(tmp,link)
    for row in p['additions']:
        (root/row['path']).symlink_to(row['target'])
    for d in sorted(p['directories'],key=lambda s:len(Path(s).parts),reverse=True):
        directory=root/d
        if directory!=root and directory.is_dir() and not any(directory.iterdir()):directory.rmdir()
    for row in p['moves']:
        assert identity(root/row['new'])==row['identity'],row
        assert not (root/row['old']).exists() or (root/row['old']).is_symlink(),row
    files,links,dirs=inventory(root)
    assert len(files)==p['regular_files']+1,(root,len(files),p['regular_files'])
    assert len(links)==p['symlinks']+len(p['additions'])
    for link in links:assert link.exists(),link
    for category in CATEGORIES:
        if (root/category).is_dir():assert not any(q.is_dir() for q in (root/category).iterdir()),category
    assert all(q.name=='ibdmix' for q in (root/'tmp').iterdir() if q.is_dir())
    assert not (root/'final').exists()
    assert (root/'gu-target.env').is_file() and not (root/'gu-target.env').is_symlink()
    receipt.write_text(json.dumps(dict(p,status='verified'),indent=2)+'\n')
    return f'{root.name}: verified {len(p["moves"])} renames and {len(p["links"])} native links; all data retain inode, size and mtime'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cohorts',nargs='+',type=Path)
    parser.add_argument('--apply',action='store_true')
    args=parser.parse_args();roots=[p.resolve() for p in args.cohorts]
    assert len(roots)==len(set(roots))
    with ThreadPoolExecutor(max_workers=4) as pool:
        plans=[]
        for p in pool.map(plan,roots):
            plans.append(p)
            print(f'{Path(p["root"]).name}: preflight {len(p["moves"])} file renames',flush=True)
        if args.apply:
            for result in pool.map(apply,plans):print(result,flush=True)

if __name__=='__main__':main()
