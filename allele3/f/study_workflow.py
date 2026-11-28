"""Input inventory and prespecified downstream contrasts (metadata-only planning)."""
import csv
import json
import os
from pathlib import Path
import sys

from product_paths import defaults, configured_path
from output_layout import data_path
from gu_segments import COHORTS, GU_CHROMS, inventory, lift_program

MODULE = {'alphagenome':'avi_window_enrichment','herv':'herv_region_enrichment',
          'archaic':'archaic_tract_enrichment','alleles':'alt_impact_similarity'}
RETAIN = {'1.ancient_pca','2.herv_scan','3.allele3_scan','4.avi_scan'}
AUTOSOMES = tuple('chr'+str(n) for n in range(1,23))

def pathenv(key, default):
    return Path(os.environ.get(key) or default).expanduser().absolute()

def restorable(path, analysis):
    path = Path(path).resolve()
    return not path.is_relative_to(analysis) or (path!=analysis and path.relative_to(analysis).parts[0] in RETAIN)

def shard_path(store, recorded):
    recorded = Path(recorded)
    local = Path(store)/'shards'/recorded.name
    return local if local.is_file() else Path(store)/recorded

def manifest(store):
    return data_path(Path(store),'qc')/'manifest.json'

def inputs(out, chroms=AUTOSOMES):
    out=Path(out).resolve()
    d=defaults()
    gen=Path(configured_path('DIR_GEN'))
    summary=Path(d['summary'])
    stores={target:summary/target for target in ('gnomad','allofus')}
    pca=pathenv('PCA3D_OUT',gen/'ancient/38/pca3d')
    detail=Path(d['detailed'])
    avi=Path(d['avi'])
    atlas=pathenv('ALLELE_ATLAS_DIR',configured_path('DIR_ANNOT','alphagenome/38'))
    pca_files=[pca/name for name in ('tables/sample_metadata.tsv','coordinates/all.sscore','basis/ref.acount',
                                      'basis/ref.eigenvec.allele','basis/ld.prune.in')]
    records={}
    def add(name,path,files,required=True,detail=''):
        files=list(files)
        missing=[str(p) for p in files if not p.is_file()]
        unsafe=[str(p) for p in files if p.is_file() and not restorable(p,out)]
        reason=detail or ('missing '+', '.join(missing[:3]) if missing else '')
        if unsafe:reason='depends on a directory marked for deletion: '+', '.join(unsafe[:3])
        records[name]=dict(path=str(path),status='unavailable' if reason else 'ready',reason=reason,required=required)
    add('pca',pca,pca_files)
    run=detail/'run.json'
    files=[run];problem=''
    if run.is_file():
        from scan_store import load_run_metadata
        try:
            for shard in load_run_metadata(detail)['shards']:
                files.extend(Path(shard)/name for name in ('COMPLETE.json','scan.json','multiallelic_sites.tsv.gz','multiallelic_alleles.tsv.gz','strata.tsv.gz'))
        except (ValueError,KeyError) as exc:problem=str(exc)
    add('detailed',detail,files,detail=problem)
    for target,store in stores.items():
        receipt=manifest(store);files=[receipt];populations=[];reason=''
        if receipt.is_file():
            meta=json.loads(receipt.read_text())
            if meta.get('status')!='complete' or meta.get('build')!='GRCh38' or meta.get('target')!=target:
                reason='expected a complete GRCh38 '+target+' summary'
            shards={r['chrom']:r for r in meta.get('shards',[])}
            for c in chroms:
                files.append(shard_path(store,shards[c]['file']) if c in shards else store/'shards'/(c+'.sites.tsv.gz'))
                if c not in shards:reason='missing '+c+' in the scan manifest'
            populations=meta.get('populations',[])
            if 'ALL' not in populations:reason='pooled ALL population is absent from the scan manifest'
        add(target,store,files,target=='gnomad',reason)
        records[target]['populations']=populations
    afile=manifest(avi)
    add('avi',avi,[afile,avi/'results/cutoffs.tsv',avi/'results/hotspots.bed',
                   *(avi/'shards'/(c+'.windows.tsv.gz') for c in chroms)],
        detail='incomplete AVI manifest' if afile.is_file() and json.loads(afile.read_text()).get('status')!='complete' else '')
    def vcf(c):
        return next((atlas/(s+ext) for s in (c,c[3:]) for ext in ('.vcf.bgz','.vcf.gz','.vcf') if (atlas/(s+ext)).is_file()),atlas/(c+'.vcf.bgz'))
    add('atlas',atlas,[vcf(c) for c in chroms])
    fasta=pathenv('ALLELE_FASTA',gen/'1hgp/38/GRCH38.fasta.gz')
    indexes=[fasta,Path(str(fasta)+'.fai')]
    if fasta.name.endswith(('.gz','.bgz')):indexes.append(Path(str(fasta)+'.gzi'))
    add('fasta',fasta,indexes)
    ancient=gen/'ancient/38/ancient'
    variants=next((Path(str(ancient)+ext) for ext in ('.pvar','.pvar.zst') if Path(str(ancient)+ext).is_file()),Path(str(ancient)+'.pvar'))
    add('ancient_genotypes',ancient,[Path(str(ancient)+'.pgen'),variants,Path(str(ancient)+'.psam'),gen/'ancient/files/HighQualityVariants.tsv'])
    aux=[Path(configured_path('DIR_FILES','liftOver',f)) for f in ('hg19ToHg38.over.chain.gz','hg38ToHg19.over.chain.gz')]
    fasta37=gen/'1hgp/37/GRCH37.fasta.gz'
    add('ancient_liftover',fasta37,[*aux,fasta37,Path(str(fasta37)+'.fai'),Path(str(fasta37)+'.gzi')])
    reference=pathenv('ALLELE_HERV_REFERENCE_BED',pathenv('ALLELE_HERV_SCAN_STORE',out/'2.herv_scan')/'results/reference_herv.bed')
    joined=pathenv('ALLELE_ERV_JOINED',configured_path('DIR_ANNOT','repeats/hg38/rmskJoinedCurrent.txt.gz'))
    herv=reference if reference.is_file() or not joined.is_file() else joined
    add('herv',herv,[herv]);records['herv']['option']='--herv-bed' if herv==reference else '--erv-joined'
    for k,env,fallback in [('low_mappability','ALLELE_LOW_MAPPABILITY_BED',configured_path('DIR_FILES','allele3/umap/hg38.k100.lt1.bed.gz')),
                           ('repeat_mask','ALLELE_REPEAT_MASK_BED',configured_path('DIR_FILES','allele3/hg38.segmental_duplications.bed.gz'))]:
        file=pathenv(env,fallback);add(k,file,[file],False)
    records['repeat_mask']['required']=True
    gu_root=pathenv('ALLELE_IBDMIX_ROOT',configured_path('ANALYSIS_ROOT','gu/ibdmix'))
    cache=pathenv('ALLELE_ARCHAIC_CACHE',out/'archaic_overlap')
    chain=pathenv('ALLELE_CHAIN37TO38',configured_path('DIR_FILES','liftOver/hg19ToHg38.over.chain.gz'))
    lift=lift_program()
    for cohort in COHORTS:
        raw=gu_root/cohort/'raw'
        units,skipped=inventory(raw,GU_CHROMS,chain,lift) if raw.is_dir() else ({},{})
        reason='' if units else ('Missing '+str(raw) if not raw.is_dir() else next(iter(skipped.values()),'No GU inputs'))
        records['archaic/'+cohort]=dict(path=str(raw),status='ready' if units else 'skipped',required=False,
            reason=reason,chain=str(chain),lift=lift,track=str(cache/cohort/'track.json'),
            skipped_chromosomes=skipped,chromosomes=list(units))
    return records

def designs(rec, results):
    rows=[]
    masks=[rec[k]['path'] for k in ('low_mappability','repeat_mask') if rec[k]['status']=='ready']
    def add(reference,target,population='ALL',design='primary',cohort=''):
        mod=MODULE[reference]; suffix={'primary':'','shared_af':'_shared_af_'+population,'mean_adjusted':'_mean_adjusted'}[design]
        dest=Path(results)/mod/cohort/(target+suffix)
        keys=[target]+({'alphagenome':['avi','fasta'],'alleles':['atlas'],'herv':['herv','fasta'],'archaic':['archaic/'+cohort,'fasta']}[reference])
        peer='allofus' if target=='gnomad' else 'gnomad'
        if design=='shared_af':keys.append(peer)
        issues=[k+': '+rec[k]['reason'] for k in keys if rec[k]['status']!='ready']
        for k in keys:
            if k in ('gnomad','allofus') and rec[k]['status']=='ready' and population not in rec[k]['populations']:
                issues.append(k+' has no '+population+' stratum')
        args=[mod,'--target',target,'--scan-dir',rec[target]['path'],'--population',population,'--out',str(dest)]
        if design=='shared_af':
            args+=['--threshold-kind','AF','--min-af','0.0001','--peer-scan-dir',rec[peer]['path']]
        if design=='mean_adjusted':args+=['--match-window-mean']
        if reference=='alphagenome':args+=['--avi-dir',rec['avi']['path']]
        if reference=='alleles':args+=['--atlas-dir',rec['atlas']['path']]
        else:args+=['--fasta',rec['fasta']['path']]
        if reference=='herv':args += [rec['herv']['option'],rec['herv']['path']]
        if reference=='archaic':
            a=rec['archaic/'+cohort];args+=['--archaic-dir',a['path'],'--archaic-track',a['track']]
        for mask in masks:args+=['--exclude-bed',mask]
        rows.append(dict(module=mod,reference=reference,target=target,gu_cohort=cohort,population=population,design=design,
                         out=str(dest),argv=args,status=('skipped' if reference=='archaic' and rec['archaic/'+cohort]['status']=='skipped' else
                                 'unavailable' if issues else 'ready'),reason=' | '.join(issues),
                         masks=' | '.join(masks) or 'none supplied'))
    for target in ('gnomad','allofus'):
        for ref in MODULE:
            if ref=='archaic':
                for cohort in COHORTS:add(ref,target,cohort=cohort)
            else:add(ref,target)
    add('alphagenome','gnomad',design='mean_adjusted')
    pops=list(dict.fromkeys(s.strip().upper() for s in os.environ.get('ALLELE_REPLICATION_POPULATIONS','ALL,AFR').split(',') if s.strip()))
    for pop in pops:
        if not pop.isalnum():raise ValueError('Invalid population label '+pop)
        for ref in ('alphagenome','alleles','herv'):
            for target in ('gnomad','allofus'):add(ref,target,pop,'shared_af')
    return rows

def show(rec):
    for key,r in rec.items():
        info=r['reason']
        if key.startswith('archaic/') and r.get('skipped_chromosomes'):
            info+=' skipped chromosomes='+','.join(r['skipped_chromosomes'])
        if key in ('gnomad','allofus') and r.get('populations'):
            info=(info+'  ' if info else '')+'populations='+','.join(r['populations'])
        print(f"{key:18} {r['status']:11} {r['path']}  {info}")
    fixes={'repeat_mask':'./0.mask.sh run',
           'pca':'./1.ancient_pca.sh run --wait TRUE',
           'detailed':'./3.allele3_scan.sh scan --target gnomad --mode detailed --wait TRUE',
           'gnomad':'./3.allele3_scan.sh scan --target gnomad --mode summary --populations all --wait TRUE',
           'allofus':'./3.allele3_scan.sh scan --target allofus --mode summary --populations all --wait TRUE',
           'avi':'./4.avi_scan.sh run --wait TRUE',
           'herv':'./2.herv_scan.sh summarize --wait TRUE (when per-sample products exist)'}
    for key,command in fixes.items():
        if rec[key]['status']!='ready':
            print('If the original product cannot be located: '+command)
        elif key in ('gnomad','allofus'):
            pops=set(os.environ.get('ALLELE_REPLICATION_POPULATIONS','ALL,AFR').upper().replace(' ','').split(','))
            if pops-set(rec[key]['populations']):
                print('For ancestry replication only, if the source supplies these strata: '+command)
    if rec['detailed']['status']!='ready':
        print('Detailed native store: set ALLELE_SCAN_STORE to its existing directory with run.json and shards; if absent, run ./3.allele3_scan.sh scan --target gnomad --mode detailed.')
    return any(r['required'] and r['status']!='ready' for r in rec.values())

def publish(out, rows, sources):
    from research_common import write_tsv
    import tempfile
    root=Path(out); dest=root/'results';dest.mkdir(parents=True,exist_ok=True)
    for row in rows:
        if row['design']=='core' and row['status']=='ready':
            if json.loads(Path(row['out']).read_text()).get('status')!='complete':raise ValueError('Missing core completion: '+row['module'])
            row['status']='complete'
        elif row['status']=='ready':
            receipt=manifest(row['out'])
            info=json.loads(receipt.read_text())
            if row['reference']=='archaic' and info.get('status')=='skipped':
                row.update(status='skipped',reason=info.get('reason','No usable GU intervals'))
            elif info.get('status')=='complete':
                row['status']='complete'
                if row['reference']=='archaic':
                    skipped=info.get('reference',{}).get('skipped_chromosomes',{})
                    row['reason']='; '.join(c+' '+r for c,r in skipped.items())
            else:
                raise ValueError('Incomplete '+str(receipt))
    with tempfile.TemporaryDirectory(dir=dest) as tmp:
        for name,fields,data in [('study_status',['module','reference','target','population','gu_cohort','design','status','reason','masks','out'],rows),
                                  ('study_inputs',['source','path','status','reason','required'],[dict(v,source=k) for k,v in sources.items()])]:
            path=Path(tmp)/(name+'.tsv')
            write_tsv(path,fields,[{k:r.get(k,'') for k in fields} for r in data]);os.replace(path,dest/path.name)
    reports=root/'reports';reports.mkdir(exist_ok=True)
    from research_common import html_report,html_table
    html_report(reports/'study_index.html','Multi-ALT study',[('Planned analyses', html_table(rows,['module','target','gu_cohort','population','design','status','reason'],limit=len(rows)))])

if __name__=='__main__':
    if sys.argv[1]=='check':
        sys.exit(int(show(inputs(defaults()['analysis']))))
    if sys.argv[1]=='publish':
        plan=json.loads(Path(sys.argv[3]).read_text())
        publish(sys.argv[2],plan['rows'],plan['sources'])
