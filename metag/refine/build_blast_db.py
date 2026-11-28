"""Build a separately named competitive nucleotide DB from local reference assemblies."""
from contextlib import ExitStack
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import sys

def read_manifest(manifest):
    with manifest.open(encoding='utf-8',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not {'taxid','fasta','role'} <= set(reader.fieldnames or []):
            raise ValueError('Provide tab-separated taxid, fasta, role')
        refs=list(reader)
    for r in refs:
        r['fasta']=str((manifest.parent/Path(r['fasta']).expanduser()).resolve(strict=True))
    return refs

def library_refs(taxids,config,root):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'f'))
    from common import load_config
    from taxonomy import Taxonomy
    from stage_alignment import _references
    cfg=load_config(config,root)
    lib=Path(cfg['library_dir'])
    cat=next((p for p in (lib/'group/optimized/references.tsv',lib/'group/references.tsv') if p.is_file()),None)
    if cat is None:
        raise FileNotFoundError('No representative catalog in '+str(lib))
    taxonomy=Taxonomy(lib/'global/taxonomy')
    with tempfile.TemporaryDirectory() as working:
        reps,_ = _references(cfg,{'representatives':str(cat)},Path(working),taxonomy)
    tids=list(dict.fromkeys(int(x.strip()) for x in taxids.split(',')))
    missing=[t for t in tids if t not in reps or not Path(reps[t]['resolved_path']).is_file()]
    if missing:
        raise ValueError('Requested taxa absent from available assemblies: '+str(missing))
    return [dict(taxid=t,fasta=reps[t]['resolved_path'],role=reps[t]['group']) for t in tids]

def build(refs, dest, *, seal_scope=None):
    if len({int(r['taxid']) for r in refs})<2:
        raise ValueError('A competitive database requires at least two taxa')
    dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists() or dest.is_symlink():
        raise FileExistsError('Choose a new reference version directory: '+str(dest))
    with tempfile.TemporaryDirectory(dir=dest.parent) as tmp, ExitStack() as streams:
        folder=Path(tmp)/'db'
        folder.mkdir()
        headers = None
        if seal_scope is not None:
            header_stream = streams.enter_context(gzip.open(folder/'contigs.tsv.gz','wt',encoding='utf-8',newline=''))
            headers = csv.writer(header_stream,delimiter='\t',lineterminator='\n')
            headers.writerow(['subject','taxid','accession','description'])
        names=[]
        records=[]
        for i,r in enumerate(refs,1):
            tid=int(r['taxid'])
            if tid<1 or not r['role'].strip():
                raise ValueError('Invalid taxid or role')
            fasta=Path(r['fasta']).resolve(strict=True)
            name=f'ref{i:03d}_{tid}'
            names.append(name)
            cmd=['makeblastdb','-in','-','-dbtype','nucl','-parse_seqids','-taxid',str(tid),
                 '-out',name,'-title',name]
            digest=hashlib.sha256()
            n=0
            op=gzip.open if fasta.suffix.lower()=='.gz' else open
            with op(fasta,'rb') as h, subprocess.Popen(cmd,cwd=folder,stdin=subprocess.PIPE) as proc:
                with proc.stdin as sink:
                    for line in h:
                        digest.update(line)
                        if line.startswith(b'>'):
                            n+=1
                            if headers is not None:
                                title=line[1:].decode('utf-8','replace').strip().replace('\t',' ')
                                token=title.split(maxsplit=1)
                                if not token:
                                    raise ValueError('Empty FASTA header in '+str(fasta))
                                headers.writerow([f'{name}_{n}',tid,token[0],title])
                            sink.write(f'>lcl|{name}_{n}\n'.encode())
                        else:
                            sink.write(line)
            if not n or proc.returncode:
                raise ValueError('makeblastdb failed or empty reference: '+str(fasta))
            records.append(dict(taxid=tid,role=r['role'],file=str(fasta),records=n,
                                decompressed_fasta_sha256=digest.hexdigest()))
        subprocess.run(['blastdb_aliastool','-dbtype','nucl','-dblist',' '.join(names),
                        '-out','competition','-title','Source audit competition'],cwd=folder,check=True)
        subprocess.run(['blastdbcmd','-db','competition','-info'],cwd=folder,check=True)
        (folder/'references.json').write_text(json.dumps(records,indent=2)+'\n')
        streams.close()
        if seal_scope is not None:
            from common import sha256
            inventory={p.name:[p.stat().st_size,p.stat().st_mtime_ns] for p in folder.iterdir()
                       if p.suffix.startswith('.n')}
            meta=dict(format='metag-seal-panel-v1',scope=seal_scope,index_files=inventory,
                      headers_sha256=sha256(folder/'contigs.tsv.gz'),refs_sha256=sha256(folder/'references.json'))
            (folder/'seal_panel.json').write_text(json.dumps(meta,indent=2)+'\n')
        folder.rename(dest)
    print(dest/'competition')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--library-taxids',help='Comma-separated species to select from the installed METAG representative catalog.')
    parser.add_argument('--config');parser.add_argument('--root')
    source.add_argument('--manifest',type=Path,help='TSV: taxid, fasta, role. Relative FASTA paths resolve next to the TSV.')
    parser.add_argument('--out',type=Path,required=True,help='New output directory; final DB prefix is DIR/competition.')
    args=parser.parse_args()
    refs=read_manifest(args.manifest.absolute()) if args.manifest else library_refs(args.library_taxids,args.config,args.root)
    build(refs,args.out.absolute())
