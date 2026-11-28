"""Prepare a frozen local Weddell seal competition panel from available assemblies."""
import argparse
import json
from pathlib import Path
import sys
import shutil
import tempfile

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'f'))
from common import load_config, locked, sha256
from stage_alignment import _references
from taxonomy import Taxonomy
from build_blast_db import build, read_manifest

TARGET=9713
GENUS=9712
# Extra backgrounds available in many local full libraries.
BACKGROUNDS=(9606,9612,9615,9685,9823,9913,10090,9031,1423,303,562,70863,10847)

def configure(args,cfg,root):
    if not args.seal_audit and not args.seal_db:
        return None
    dest=Path(args.seal_db).expanduser().absolute() if args.seal_db else root/'references'/'weddell'
    if dest.name=='competition' and not dest.is_dir():
        dest=dest.parent
    return dict(path=str(dest),auto=not bool(args.seal_db),library_dir=cfg['library_dir'],
                reference_raw=cfg['reference_raw'],work_dir=cfg.get('index_base_work_dir',cfg['work_dir']),
                metadata_dir=cfg['metadata_dir'],root=cfg['root'])

def scope(refs,tax):
    chain=list(tax.lineage(TARGET))
    family=next(t for t in chain if tax.rank[t]=='family')
    ids={int(r['taxid']) for r in refs}
    lineages={i:set(tax.lineage(i)) for i in ids}
    relatives=sorted({tax.describe(i)['species_taxid'] for i,c in lineages.items()
                      if family in c and TARGET not in c and tax.describe(i)['species_taxid']})
    present=dict(target=TARGET in ids,phocid_relative=bool(relatives),human=9606 in ids,
                 other_mammal=any(40674 in c and 9606 not in c and family not in c for c in lineages.values()),
                 microbial_background=any(2 in c for c in lineages.values()))
    return dict(family=family,relative_species=relatives,components=present,
                adequate=all(present.values()),taxids=sorted(ids))

def from_library(cfg,tax):
    lib=Path(cfg['library_dir'])
    catalogs=[lib/'group/optimized/references.tsv',lib/'group/references.tsv']
    cat=next((p for p in catalogs if p.is_file()),None)
    if cat is None:
        raise FileNotFoundError('No local representative catalog under '+str(lib))
    with tempfile.TemporaryDirectory() as tmp:
        reps,_=_references(cfg,{'representatives':str(cat)},Path(tmp),tax)
    def available(i):
        return i in reps and Path(reps[i]['resolved_path']).is_file()
    if not available(TARGET):
        raise ValueError('No Weddell seal assembly in local catalog. Supply a separate seal panel manifest.')
    target_line=list(tax.lineage(TARGET))
    family=next(t for t in target_line if tax.rank[t]=='family')
    relatives=[]
    for i,r in reps.items():
        if i==TARGET or r['group']!='mammal':
            continue
        lineage=set(tax.lineage(i))
        if family in lineage and available(i):
            distance=next(j for j,t in enumerate(target_line) if t in lineage)
            relatives.append((distance,i))
    chosen=[TARGET,*[i for _,i in sorted(relatives)[:4]],*[i for i in BACKGROUNDS if available(i)]]
    return [dict(taxid=i,fasta=reps[i]['resolved_path'],role=reps[i]['group']) for i in dict.fromkeys(chosen)]

def signature(folder):
    marker=folder/'seal_panel.json'
    data=json.loads(marker.read_text())
    if data['format']!='metag-seal-panel-v1' or not data['index_files']:
        raise ValueError('Expected a named seal competition panel in '+str(folder))
    for name,(size,mtime) in data['index_files'].items():
        st=(folder/name).stat()
        if st.st_size!=size or st.st_mtime_ns!=mtime:
            raise ValueError('Changed BLAST volume: '+str(folder/name))
    if data['refs_sha256']!=sha256(folder/'references.json') or not (folder/'contigs.tsv.gz').is_file():
        raise ValueError('Incomplete seal panel: '+str(folder))
    return sha256(marker), data

def prepare(opt):
    dest=Path(opt['path'])
    needed=['blastn'] if dest.exists() else ['blastn','makeblastdb','blastdb_aliastool','blastdbcmd']
    missing=[x for x in needed if shutil.which(x) is None]
    if missing:
        raise RuntimeError('Load NCBI BLAST+ in the worker environment; missing: '+', '.join(missing))
    if dest.exists():
        signature(dest)
        return dest
    if not opt['auto']:
        raise FileNotFoundError('Run seal_panel.py to build '+str(dest))
    with locked(dest.parent/'.weddell-build.lock'):
        if not dest.exists():
            tax=Taxonomy(Path(opt['library_dir'])/'global/taxonomy')
            refs=from_library(opt,tax)
            build(refs,dest,seal_scope=scope(refs,tax))
    signature(dest)
    return dest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,help='TSV: taxid, fasta, role. One labeled assembly per species.')
    p.add_argument('--out',required=True,type=Path,help='New version directory.')
    p.add_argument('--config');p.add_argument('--root')
    a=p.parse_args()
    cfg=load_config(a.config,a.root)
    tax=Taxonomy(Path(cfg['library_dir'])/'global/taxonomy')
    refs=read_manifest(a.manifest.resolve()) if a.manifest else from_library(cfg,tax)
    ids=[int(r['taxid']) for r in refs]
    if TARGET not in ids or len(ids)!=len(set(ids)):
        p.error('Include 9713 and supply one assembly per taxid in this panel. Use separate panels for alternative target assemblies.')
    build(refs,a.out.absolute(),seal_scope=scope(refs,tax))

if __name__=='__main__':
    main()
