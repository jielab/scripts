"""Refine published splits using community-level evidence. See metag.sh refine -h."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'f'))
from common import atomic_json, fingerprint, load_config, locked, resolve_mode, sha256
from rank_workflow import FORMAT, release
from rank_paper import discover, read_sample
from ranked_storage import staging, recover_directory, publish_directory
from scheduler import submit_cmd, write_cmd
from storage_layout import new_run
from taxonomy import Taxonomy
from metag_refine_tables import SCHEMA, code_id, rows, write
from signals import sources, virus
import blast_audit
import resistome
import cache
import seal_panel
import seal_audit

def ref_version():
    return fingerprint([release({})[0],code_id()])

def meta(path):
    if path is None:
        return {}
    out = {}
    with open(path,encoding='utf-8',newline='') as f:
        r = csv.DictReader(f,delimiter='\t')
        if not {'sample','sample_type','library_type'} <= set(r.fieldnames or []):
            raise ValueError('Metadata requires sample, sample_type, library_type (tab-separated)')
        for row in r:
            key,mat,lib = (row[k].strip() for k in ('sample','sample_type','library_type'))
            if (key in out or len(key.split('/'))!=2 or not all(key.split('/')) or
                mat not in ('stool','wastewater','environment','blank','mock','other','unknown') or
                lib not in ('DNA','RNA','DNA_RNA','unknown')):
                raise ValueError('Invalid or duplicate metadata row for '+key)
            out[key] = [mat,lib]
    return out

def run_sample(folder,root,opts,condition,version,tax_cache,assay_keys=None,expected_split=None):
    sample = folder.parent.name+'/'+folder.name
    dest = root/'main'/sample
    with read_sample(folder), locked(dest.parent/('.'+dest.name+'.lock')):
        s = json.loads((folder/'summary.json').read_text())
        if s['sample'] != sample or s['format'] != FORMAT:
            raise ValueError('Split receipt differs from selected sample '+sample)
        split_sha = sha256(folder/'summary.json')
        if expected_split is not None and split_sha!=expected_split:
            raise ValueError('Split changed after planning. Resubmit '+sample)
        recover_directory(dest)
        if dest.is_symlink() or (dest.exists() and not (dest/'receipt.json').is_file()):
            raise ValueError('Unexpected refine destination '+str(dest))
        old = json.loads((dest/'receipt.json').read_text()) if dest.exists() else {}
        if dest.exists() and old.get('format') != cache.REFINE_FORMAT:
            raise ValueError('Expected current refine receipt in '+str(dest))
        keys = assay_keys if assay_keys is not None else {a:cache.assay_key(a,opts) for a in cache.requested(opts)}
        if cache.complete(dest,old,sample,split_sha,opts,condition,keys):
            print('Reuse refine '+sample,flush=True)
            return
        reuse = old.get('split_sha256') == split_sha and old.get('sample') == sample
        if condition is None:
            condition = ([old['sample_type'],old['library_type']] if old.get('split_sha256')==split_sha
                         else ['unknown','unknown'])
        rec = dict(format=cache.REFINE_FORMAT,sample=sample,state='complete',sample_type=condition[0],library_type=condition[1],
                   input_pairs=s['pairs'],split_sha256=split_sha,refine_version=version,baseline_key=cache.baseline_key(),
                   blast='not_run',amr='not_run',blast_db_id='NA',amr_marker_sha256='NA',
                   split_method=s['request']['code'],split_database=s['request']['database'])
        taxa = list(rows(folder/'taxa.tsv'))
        if sum(int(r['pairs']) for r in taxa) != s['pairs']:
            raise ValueError('Taxon total differs from receipt: '+sample)
        with staging(dest.parent,'.building-') as out:
            with TemporaryDirectory(dir=out) as scratch:
                tmp = Path(scratch)
                source = sources(sample,taxa,rows(folder/'taxon_evidence.tsv'))
                data = {'refine_sources':source, 'refine_virus':virus(sample,taxa,s['pairs'],condition[1])}
                for assay,option in [('blast','blast_db'),('amr','amr_markers')]:
                    table = 'refine_'+assay
                    if reuse and old.get(assay)=='complete':
                        previous = dest/(table+'.tsv')
                        valid = cache.valid_table(dest,old,table)
                        keep = (opts[option] is None or
                                (not opts.get('replace') and old.get(assay+'_key')==keys.get(assay)))
                        if keep and valid:
                            with previous.open(encoding='utf-8') as f:
                                if f.readline().rstrip('\r\n').split('\t') != SCHEMA[table]:
                                    raise ValueError('Unexpected table header: '+str(previous))
                            shutil.copyfile(previous,out/(table+'.tsv'))
                            rec.update({k:v for k,v in old.items() if k==assay or k.startswith(assay+'_')})
                            continue
                        if keep and opts[option] is None:
                            raise ValueError('Damaged optional assay table. Supply its references to repeat: '+str(previous))
                    if opts[option] is None:
                        data[table] = []
                        continue
                    rec[assay+'_key'] = keys[assay]
                    if assay=='blast':
                        db = opts[option]
                        info = subprocess.check_output(['blastdbcmd','-db',db,'-info'],text=True)
                        tool = subprocess.check_output(['blastn','-version'],text=True).splitlines()[0]
                        inventory = Path(db).parent/'references.json'
                        local_inventory = Path(db).name == 'competition' and inventory.is_file()
                        reference_sha = sha256(inventory) if local_inventory else None
                        reference_ids = {int(r['taxid']) for r in json.loads(inventory.read_text())} if local_inventory else None
                        rec.update(blast_db_id=fingerprint([db,info,reference_sha]),blast_reference_manifest_sha256=reference_sha,blast_db=db,blast_db_info=info,
                                   blast_tool=tool,blast_pairs_per_bin=opts['blast_pairs'])
                        tax_key = s['taxonomy'],s['request']['database']
                        if tax_key not in tax_cache:
                            tax_cache.clear()
                            tax_cache[tax_key] = Taxonomy(s['taxonomy'])
                        data[table] = blast_audit.audit(sample,folder,source,taxa,db,tax_cache[tax_key],
                                                       opts['blast_pairs'],opts['cores'],tmp,reference_ids)
                        if info != subprocess.check_output(['blastdbcmd','-db',db,'-info'],text=True):
                            raise ValueError('BLAST database changed during audit')
                    else:
                        markers = Path(opts[option])
                        rec['amr_marker_sha256'] = sha256(markers)
                        data[table],rec['amr_settings'] = resistome.run(sample,s,markers,opts['cores'],tmp,
                                                                       Path(opts['amr_annotation']) if opts['amr_annotation'] else None)
                        if sha256(markers) != rec['amr_marker_sha256']:
                            raise RuntimeError('AMR markers changed during analysis')
                    if cache.assay_key(assay,opts) != keys[assay]:
                        raise RuntimeError('Assay inputs changed during '+assay)
                    rec[assay] = 'complete'
                amr = rec.get('amr_settings',{})
                rec.update(amr_tested_families=amr.get('tested_families','NA'),
                           amr_ambiguous_pairs=amr.get('ambiguous_pairs','NA'),
                           amr_engine=amr.get('engine','NA'))
                rec.setdefault('blast_pairs_per_bin','NA')
                data['refine_status'] = [[rec[k] for k in SCHEMA['refine_status']]]
                for name, items in data.items():
                    write(out/(name+'.tsv'),SCHEMA[name],items)
                rec['table_hashes'] = {name:sha256(out/(name+'.tsv')) for name in SCHEMA}
                atomic_json(out/'receipt.json',rec)
            publish_directory(out,dest)
    print('Published '+str(dest),flush=True)

def main():
    if len(sys.argv)==3 and sys.argv[1]=='--worker':
        spec=json.loads(Path(sys.argv[2]).read_text())
        if spec['version']!=ref_version():
            raise RuntimeError('Installed scripts changed after submission. Resubmit refine.')
        tax_cache = {}
        keys={a:cache.assay_key(a,spec['opts']) for a in cache.requested(spec['opts'])}
        for path,condition in spec['inputs']:
            folder=Path(path)
            expected=spec.get('expected_splits',{}).get(path)
            run_sample(folder,Path(spec['root']),spec['opts'],condition,spec['version'],tax_cache,keys,expected)
            if spec['opts'].get('seal'):
                seal_audit.review(folder,Path(spec['root']),spec['opts'],expected)
        return
    p=argparse.ArgumentParser(prog='./metag.sh refine',description=__doc__)
    p.add_argument('--dat-list',metavar='CSV',help='All completed ranked splits in these datasets.')
    p.add_argument('--sample',action='append',metavar='DATASET/SAMPLE',help='Select completed sample (repeatable).')
    p.add_argument('--seal-audit',action='store_true')
    p.add_argument('--seal-db',help='Version directory made by seal_panel.py; enables Weddell audit.')
    p.add_argument('--metadata',metavar='TSV',help='Columns sample, sample_type, library_type.')
    p.add_argument('--blast-db',default=os.getenv('METAG_BLAST_DB'),help='Local competitive BLAST nucleotide DB prefix; optional.')
    p.add_argument('--blast-pairs',type=int,default=100,help='Uniform pairs per eligible bin [100].')
    p.add_argument('--amr-markers',default=os.getenv('METAG_AMR_MARKERS'),help='Curated ShortBRED AMR protein markers FASTA; optional.')
    p.add_argument('--amr-annotation',default=os.getenv('METAG_AMR_ANNOTATION'),help='Optional ShortBRED family annotation TSV.')
    p.add_argument('--config');p.add_argument('--root')
    p.add_argument('--submit',choices=['auto','hpc','local'],default='auto')
    p.add_argument('--queue',default=None);p.add_argument('--cores',type=int,default=None)
    p.add_argument('--mem-gib',type=int,default=32);p.add_argument('--hours',type=int,default=24)
    p.add_argument('--batch-size',type=int,help='Samples per cohort job [all for tables, 4 with sequence searches].')
    p.add_argument('--replace',choices=('FALSE','TRUE'),default='FALSE')
    p.add_argument('--plan-only',action='store_true')
    a=p.parse_args()
    if a.amr_annotation and not a.amr_markers:
        p.error('--amr-annotation requires --amr-markers')
    if not a.dat_list and not a.sample:
        p.error('Select --dat-list or --sample')
    cores=a.cores or (4 if a.amr_markers or a.blast_db or a.seal_audit or a.seal_db else 1)
    if a.cores==0 or min(cores,a.blast_pairs,a.mem_gib,a.hours,a.batch_size if a.batch_size is not None else 1)<1:
        p.error('Resources and --blast-pairs must be positive')
    mode=resolve_mode(a.submit)
    os.environ['submit']=os.environ['SUBMIT']=mode
    cfg=load_config(a.config,a.root)
    root=Path(cfg['work_dir'])/'refine'
    folders=discover(cfg,vars(a))
    labels=meta(a.metadata)
    opts=dict(replace=a.replace=='TRUE',cores=cores,blast_pairs=a.blast_pairs,blast_db=str(Path(a.blast_db).expanduser().absolute()) if a.blast_db else None,
              amr_markers=str(Path(a.amr_markers).expanduser().resolve(strict=True)) if a.amr_markers else None,amr_annotation=str(Path(a.amr_annotation).expanduser().resolve(strict=True)) if a.amr_annotation else None)
    opts['seal']=seal_panel.configure(a,cfg,root)
    cfg.update(work_dir=str(root),classify_queue=a.queue or cfg['classify_queue'],classify_cores=cores,
               classify_hours=a.hours,classify_memory_gib=a.mem_gib)
    version=ref_version()
    expected_splits={}
    keys=None
    # Bound sequence-search jobs while keeping lightweight summaries in one job per cohort.
    cohorts={}
    for path in folders:
        sample=path.parent.name+'/'+path.name
        if a.metadata and sample not in labels:
            p.error('Supplied metadata missing selected sample '+sample)
        dest=root/'main'/sample
        with read_sample(path):
            expected_splits[str(path)]=sha256(path/'summary.json')
            if not opts['replace'] and (dest/'receipt.json').is_file():
                with locked(dest.parent/('.'+dest.name+'.lock'),shared=True):
                    if keys is None:
                        keys={assay:cache.assay_key(assay,opts) for assay in cache.requested(opts)}
                    if (cache.complete(dest,cache.receipt(dest),sample,expected_splits[str(path)],opts,labels.get(sample),keys)
                        and seal_audit.current(root,sample,expected_splits[str(path)],opts)):
                        print('Reuse refine '+sample,flush=True)
                        continue
        cohorts.setdefault(path.parent.name,[]).append([str(path),labels.get(sample)])
    if not cohorts:
        print('All selected evidence steps are complete.',flush=True)
        return
    jobs=new_run(root/'jobs')
    for cohort,inputs in cohorts.items():
        size = a.batch_size or (4 if opts['blast_db'] or opts['amr_markers'] or opts.get('seal') else len(inputs))
        for first in range(0,len(inputs),size):
            unit = cohort if len(inputs)<=size else f'{cohort}.batch{first//size+1:03d}'
            job=jobs/(unit+'.json')
            atomic_json(job,dict(inputs=inputs[first:first+size],opts=opts,version=version,root=str(root),
                                 expected_splits={path:expected_splits[path] for path,_ in inputs[first:first+size]}))
            cmd_args=[sys.executable,'-B',str(Path(__file__).resolve()),'--worker',str(job)]
            cmd=write_cmd(cfg,jobs/(unit+'.cmd'),'classify',cmd_args,'metag.refine.'+unit)
            print(cmd,flush=True)
            if not a.plan_only:
                if mode=='local':
                    subprocess.run(cmd_args,check=True)
                else:
                    submit_cmd(cfg,cmd)


if __name__=='__main__':
    main()
