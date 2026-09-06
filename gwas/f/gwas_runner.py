#!/usr/bin/env python3
"""Independent WSL GWAS modules, resumable command plans and canonical outputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
SELF=str(Path(__file__).resolve())
SAIGE_ENV=next((Path.home()/p for p in ('anaconda3/envs/saige','miniforge3/envs/saige')
                if (Path.home()/p/'bin/createSparseGRM.R').exists()),Path.home()/'anaconda3/envs/saige')
TRAITS={'height':'qt','bald':'ordinal','bald12':'bt','cvd_cad.Yt2e':'bt','cvd_cad.t2e':'t2e'}

def boolean(x):
    if x.upper() not in ('TRUE','FALSE'): raise argparse.ArgumentTypeError('Use TRUE or FALSE')
    return x.upper()=='TRUE'

def tool(name):
    return shutil.which(name) or (f'/mnt/d/software/bin/{name}' if Path(f'/mnt/d/software/bin/{name}').exists() else name)

def parse():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,epilog='''
Examples (WSL; these are independent modules):
  # Generate the phenotype file first:
  bash /mnt/d/scripts/ukb/phe.sh --steps phe4gwas

  # Once only: prune typed genotypes and create the shared sparse GRM.
  bash /mnt/d/scripts/gwas/gwas.sh prep_gwas --sparse-grm TRUE --threads 16

  # If only the REGENIE pruning input is needed:
  bash /mnt/d/scripts/gwas/gwas.sh prep_gwas --sparse-grm FALSE --threads 16
  # Later, re-run prep_gwas with --sparse-grm TRUE to add the GRM; pruning resumes.

  # Formal REGENIE run: all five traits, auto-selecting qt/ordinal/bt/t2e.
  # Activate your REGENIE >=4 environment first, or pass --regenie /path/to/regenie.
  bash /mnt/d/scripts/gwas/gwas.sh run_regenie --grch 37 --threads 16

  # Equivalent explicit trait list:
  bash /mnt/d/scripts/gwas/gwas.sh run_regenie --grch 37 \\
    --pheno /mnt/d/data/ukb/phe/common/ukb.phe --pheno-name height,bald,bald12,cvd_cad.Yt2e,cvd_cad.t2e --type auto

  # Review commands without starting GWAS; each method writes run.cmd.sh/run.plan.json.
  bash /mnt/d/scripts/gwas/gwas.sh run_regenie --grch 38 --run FALSE \\
    --output-dir /mnt/d/data.BIG/gwas/self/grch38/common

  # Alternative methods (T2E is explicitly skipped):
  bash /mnt/d/scripts/gwas/gwas.sh run_plink2 --grch 37
  bash /mnt/d/scripts/gwas/gwas.sh run_saige --grch 37

  # Compare the original HPC bald GWAS to the new REGENIE result:
  bash /mnt/d/scripts/0data/gwas_post.sh compare \\
    --gwas-files /mnt/d/data.BIG/gwas/main/common/bald/gwas/bald.gz,/mnt/d/data.BIG/gwas/self/common/bald/gwas/bald.regenie.gz \\
    --mplot TRUE --compare-beta TRUE --compare-EAF TRUE \\
    --output-dir /mnt/d/data.BIG/gwas/self/common/bald/qc/compare

Output: <output-dir>/<trait>/gwas/<trait>.<method>.gz
Raw chromosome files: <output-dir>/<trait>/gwas/<method>/
Default chromosomes: 1-22; add --chromosomes 1,2,...,22,X if chrX is needed.
Ordinal bald is modeled as a quantitative score (no rank normalization).
T2E duration cvd_cad.t2e is paired with event cvd_cad.Yt2e.
Shared prep files always default to /mnt/h/ukbGen/37/typ, including GRCh38 runs.
Completed jobs resume only when their input/configuration fingerprint matches.
''')
    p.add_argument('module',choices=['prep_gwas','run_plink2','run_regenie','run_saige'])
    p.add_argument('--grch',choices=['37','38'],default='37')
    p.add_argument('--threads',type=int,default=16)
    p.add_argument('--memory',type=int,default=24000,help='PLINK memory in MiB')
    p.add_argument('--typed-dir',default='/mnt/h/ukbGen/37/typ',help='Shared build-37 pruning/GRM cache')
    p.add_argument('--imputed-dir',help='Default /mnt/h/ukbGen/<grch>/imp')
    p.add_argument('--pheno','--phenotype-file',dest='phenotype_file',metavar='FILE',default='/mnt/d/data/ukb/phe/common/ukb.phe',help='Phenotype file (FID IID and trait columns)')
    p.add_argument('--covariate-file')
    p.add_argument('--pheno-name','--phenotypes','--traits',dest='phenotypes',metavar='CSV',default=','.join(TRAITS),help='Comma-separated trait names; default: the five example traits')
    p.add_argument('--type',choices=['auto','qt','bt','ordinal','t2e'],default='auto')
    p.add_argument('--event-col',help='Default <trait stem>.Yt2e for *.t2e')
    p.add_argument('--covariates',default='age,sex,tdi,PC1,PC2,PC3,PC4')
    p.add_argument('--categorical-covariates',default='center')
    p.add_argument('--output-dir',default='/mnt/d/data.BIG/gwas/self/common',help='Trait-root directory')
    p.add_argument('--chromosomes',default=','.join(map(str,range(1,23))))
    p.add_argument('--extract',help='Optional variant IDs for test/subset association; does not alter prep')
    p.add_argument('--keep',help='Optional FID IID subset for test runs (all model stages)')
    p.add_argument('--min-mac',type=int,default=100)
    p.add_argument('--plink2',default=tool('plink2'))
    p.add_argument('--plink',default=tool('plink'))
    p.add_argument('--regenie',default=tool('regenie'))
    p.add_argument('--rscript',default=tool('Rscript'))
    p.add_argument('--saige-rscript',default=str(HERE/'saige_Rscript.sh') if (SAIGE_ENV/'bin/Rscript').exists() else tool('Rscript'))
    p.add_argument('--saige-dir',default=str(SAIGE_ENV/'bin') if (SAIGE_ENV/'bin/createSparseGRM.R').exists() else '/mnt/d/software/SAIGE/extdata')
    p.add_argument('--sparse-grm',type=boolean,default=True,help='prep_gwas: FALSE prunes only; TRUE also creates sparse GRM')
    p.add_argument('--run',type=boolean,default=True,help='FALSE: generate reviewable command plans only')
    p.add_argument('--replace',type=boolean,default=False,help='Re-run completed jobs with identical configuration')
    a=p.parse_args()
    if a.threads<1 or a.memory<640 or a.min_mac<1: p.error('Invalid threads/memory/min-mac')
    a.chromosomes=a.chromosomes.replace(',',' ').split()
    if not a.chromosomes or len(set(a.chromosomes))!=len(a.chromosomes) or any(x not in list(map(str,range(1,23)))+['X'] for x in a.chromosomes): p.error('Invalid chromosomes')
    a.chromosomes.sort(key=lambda c:23 if c=='X' else int(c))
    if a.module=='prep_gwas' and 'X' in a.chromosomes: p.error('Prune/GRM uses autosomes only')
    a.imputed_dir=a.imputed_dir or f'/mnt/h/ukbGen/{a.grch}/imp'
    a.covariate_file=a.covariate_file or a.phenotype_file
    return a

def pfile(prefix):
    return ['--pfile',str(prefix)]+(['vzs'] if Path(str(prefix)+'.pvar.zst').exists() else [])

def node(name,cmd,outputs,inputs=(),covars=None):
    return dict(name=name,cmd=list(map(str,cmd)),outputs=list(map(str,outputs)),inputs=list(map(str,inputs)),covars=covars)

def execute(plan):
    # Serialize jobs in one plan; never let two invocations overwrite a shared cache.
    import fcntl
    lock=open(str(plan)+'.lock','w')
    try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: raise RuntimeError(f'Already running: {plan}')
    spec=json.loads(Path(plan).read_text())
    for job in spec['jobs']:
        executable=job['cmd'][0]
        if not shutil.which(executable): raise FileNotFoundError('Required executable: '+executable)
        if len(job['cmd'])>1 and job['cmd'][1].endswith('.R') and not Path(job['cmd'][1]).is_file():
            raise FileNotFoundError('Required R script: '+job['cmd'][1]+' (set --saige-dir/--saige-rscript)')
    for job in spec['jobs']:
        for f in job['inputs']:
            if not Path(f).is_file(): raise FileNotFoundError(f)
        cmd=job['cmd'][:]
        if job['covars']:
            path,method=job['covars']; names=Path(path).read_text().strip()
            if names:
                cmd+= {'plink2':['--covar',path.removesuffix('.covars'),'--covar-name',names,'--covar-variance-standardize'],
                       'regenie':['--covarFile',path.removesuffix('.covars'),'--covarColList',names],
                       'saige':['--covarColList='+names]}[method]
        signature=hashlib.sha256(json.dumps([cmd,[(f,Path(f).stat().st_size,Path(f).stat().st_mtime_ns) for f in job['inputs']]],sort_keys=True).encode()).hexdigest()
        stamp=Path(plan).parent/(job['name']+'.done.json')
        if not spec['replace'] and stamp.exists() and json.loads(stamp.read_text()).get('signature')==signature and all(Path(f).is_file() and (Path(f).stat().st_size or str(f).endswith('.covars')) for f in job['outputs']):
            print('RESUME '+job['name'],flush=True); continue
        # Invalidate success before execution so failures cannot be mistaken for completion.
        stamp.unlink(missing_ok=True)
        print(shlex.join(cmd),flush=True)
        with open(Path(plan).parent/(job['name']+'.run.log'),'w') as log:
            result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:
            logpath=Path(plan).parent/(job['name']+'.run.log')
            print('\n'.join(logpath.read_text(errors='replace').splitlines()[-12:]),file=sys.stderr)
            raise RuntimeError(f'{job["name"]} failed (exit {result.returncode}); see {logpath}')
        if not all(Path(f).is_file() and (Path(f).stat().st_size or str(f).endswith('.covars')) for f in job['outputs']): raise RuntimeError('Missing output: '+job['name'])
        stamp.write_text(json.dumps({'signature':signature,'command':cmd}))

def save_plan(directory,jobs,a):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    plan=directory/'run.plan.json'
    spec={'config':vars(a),'replace':a.replace,'jobs':jobs}
    # Do not silently reuse results across different builds, subsets or model settings.
    if plan.exists():
        old=json.loads(plan.read_text()); oc=old.get('config',{}).copy(); nc=vars(a).copy()
        for key in ('run','replace','sparse_grm','saige_dir','saige_rscript'): oc.pop(key,None); nc.pop(key,None)
        if oc!=nc and not a.replace and any(directory.glob('*.done.json')): raise ValueError(f'Configuration changed at {directory}; choose a new output directory or --replace TRUE')
    plan.write_text(json.dumps(spec,indent=2)+'\n')
    lines=['#!/usr/bin/env bash','set -euo pipefail',f'# Commands and outputs: {plan}']
    lines.extend('# '+shlex.join(j['cmd']) for j in jobs)
    lines.append(shlex.join(['python3',SELF,'execute-plan',str(plan)]))
    script=directory/'run.cmd.sh'; script.write_text('\n'.join(lines)+'\n')
    print(script,flush=True)
    if a.run: execute(plan)

def genotype_inputs(prefix):
    return [str(prefix)+'.pgen',str(prefix)+('.pvar.zst' if Path(str(prefix)+'.pvar.zst').exists() else '.pvar'),str(prefix)+'.psam']

def prep(a):
    root=Path(a.typed_dir).resolve(); cache=root/'.prep'; cache.mkdir(parents=True,exist_ok=True)
    jobs=[]; beds=[]
    base=[a.plink2,'--threads',a.threads,'--memory',a.memory]
    for c in a.chromosomes:
        src=root/f'chr{c}'; dst=cache/f'chr{c}.prune'
        # Explicit k=0 preserves the original fixed 1e-15 HWE cutoff and
        # acknowledges PLINK2's large-sample strictness check for this GRM panel.
        prune=base+pfile(src)+['--maf','.01','--mac','100','--max-alleles','2','--geno','.1','--hwe','1e-15','0','--indep-pairwise','1000kb','0.1','--out',dst]
        if c=='6': prune+=['--exclude','range',cache/'mhc.range']
        (cache/'mhc.range').write_text('6\t25000000\t35000000\tMHC\n')
        jobs.append(node(f'chr{c}.prune',prune,[str(dst)+'.prune.in'],genotype_inputs(src)))
        outputs=[str(dst)+ext for ext in ('.bed','.bim','.fam')]; beds.append(str(dst))
        jobs.append(node(f'chr{c}.bed',base+pfile(src)+['--extract',str(dst)+'.prune.in','--make-bed','--out',dst],outputs,genotype_inputs(src)+[str(dst)+'.prune.in']))
    merged=root/'typ.prune'
    mergefile=cache/'merge.list'; mergefile.write_text('\n'.join(beds[1:])+'\n')
    cmd=[a.plink,'--bfile',beds[0],'--make-bed','--out',merged,'--threads',a.threads,'--memory',a.memory]
    if len(beds)>1: cmd+=['--merge-list',mergefile]
    jobs.append(node('merge',cmd,[str(merged)+x for x in ('.bed','.bim','.fam')],[b+e for b in beds for e in ('.bed','.bim','.fam')]))
    if not a.sparse_grm:
        save_plan(cache,jobs,a)
        return
    prefix=root/'typ'; actual=str(prefix)+'_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx'
    jobs.append(node('sparseGRM',[a.saige_rscript,Path(a.saige_dir)/'createSparseGRM.R','--plinkFile='+str(merged),f'--nThreads={a.threads}','--numRandomMarkerforSparseKin=2000','--relatednessCutoff=0.125','--outputPrefix='+str(prefix)],
                     [actual,actual+'.sampleIDs.txt'],[str(merged)+e for e in ('.bed','.bim','.fam')]))
    jobs.append(node('grm-links',['python3',SELF,'grm-links',actual,str(prefix)], [str(prefix)+'.sparseGRM.mtx',str(prefix)+'.sparseGRM.mtx.sampleIDs.txt'],[actual,actual+'.sampleIDs.txt']))
    save_plan(cache,jobs,a)

def association(a):
    method=a.module.removeprefix('run_')
    for trait in a.phenotypes.split(','):
        if not trait or any(c in trait for c in '/\\\n\r') or trait in ('.','..'): raise ValueError('Invalid trait name')
        typ=TRAITS.get(trait) if a.type=='auto' else a.type
        if typ is None: raise ValueError('Specify --type for '+trait)
        root=Path(a.output_dir).resolve()/trait/'gwas'; work=root/method
        work.mkdir(parents=True,exist_ok=True)
        if typ=='t2e' and method!='regenie':
            message=f'{method}: time-to-event is unsupported; use run_regenie --type t2e'
            (work/'UNSUPPORTED.txt').write_text(message+'\n'); print(message,flush=True); continue
        event=a.event_col or (trait.removesuffix('.t2e')+'.Yt2e')
        phe=work/'analysis.phe'; covars=str(phe)+'.covars'
        jobs=[node('phenotype',[a.rscript,HERE/'gwas_inputs.R',a.phenotype_file,a.covariate_file,trait,typ,event,a.covariates,a.categorical_covariates,phe],
                   [phe,covars],[a.phenotype_file,a.covariate_file,str(HERE/'gwas_inputs.R')])]
        pruned=Path(a.typed_dir)/'typ.prune'; grm=Path(a.typed_dir)/'typ.sparseGRM.mtx'
        keep=['--keep',a.keep] if a.keep else []
        extra=['--extract',a.extract] if a.extract else []
        common_inputs=[phe,covars]+([a.keep] if a.keep else [])
        if method=='regenie':
            flags=['--bt'] if typ=='bt' else ['--t2e'] if typ=='t2e' else ['--qt']
            phflags=['--phenoColList',trait,'--eventColList',event] if typ=='t2e' else ['--phenoCol',trait]
            common=[a.regenie,'--threads',a.threads,'--bsize','1000','--phenoFile',phe]+phflags+flags+keep
            step1=work/'step1'
            jobs.append(node('null-model',common+['--step','1','--bed',pruned,'--lowmem','--lowmem-prefix',work/'tmp','--out',step1],
                             [str(step1)+'_pred.list',str(step1)+'_1.loco'],common_inputs+[str(pruned)+e for e in ('.bed','.bim','.fam')],(covars,method)))
        if method=='saige':
            if a.keep: raise ValueError('SAIGE --keep: provide a subset phenotype-file instead')
            step1=work/'step1'; sparse=['--sparseGRMFile='+str(grm),'--sparseGRMSampleIDFile='+str(grm)+'.sampleIDs.txt']
            cmd=[a.saige_rscript,Path(a.saige_dir)/'step1_fitNULLGLMM.R','--plinkFile='+str(pruned),
                 '--phenoFile='+str(phe),'--phenoCol='+trait,'--sampleIDColinphenoFile=IID',
                 '--traitType='+('binary' if typ=='bt' else 'quantitative'),'--invNormalize=FALSE',
                 '--useSparseGRMtoFitNULL=TRUE','--skipVarianceRatioEstimation=FALSE','--IsOverwriteVarianceRatioFile=TRUE',
                 '--LOCO=FALSE',f'--nThreads={a.threads}','--outputPrefix='+str(step1)]+sparse
            jobs.append(node('null-model',cmd,[str(step1)+'.rda',str(step1)+'.varianceRatio.txt'],common_inputs+[str(pruned)+e for e in ('.bed','.bim','.fam')]+[grm,str(grm)+'.sampleIDs.txt'],(covars,method)))
        raw=[]
        for c in a.chromosomes:
            prefix=Path(a.imputed_dir)/f'chr{c}'; dest=work/f'chr{c}'; inputs=genotype_inputs(prefix)+common_inputs+([a.extract] if a.extract else [])
            if method=='plink2':
                result=str(dest)+'.'+trait+'.glm.'+('logistic.hybrid' if typ=='bt' else 'linear')
                glm=['hide-covar','allow-no-covars','omit-ref','no-x-sex','cols=+omitted,+a1freq,+nobs']
                if typ=='bt': glm+=['firth-fallback']
                cmd=[a.plink2,'--threads',a.threads,'--memory',a.memory]+pfile(prefix)+['--glm']+glm+['--pheno',phe,'--pheno-name',trait,'--no-psam-pheno','--mac',a.min_mac,'--max-alleles','2','--out',dest]+keep+extra
                if typ=='bt': cmd+=['--1']
                if c=='X': cmd+=['--xchr-model','2']
            elif method=='regenie':
                result=str(dest)+'_'+trait+'.regenie'
                # REGENIE needs a plain .pvar. Keep genotype data in place; create
                # a per-analysis view with symlinks and a decompressed variant table.
                view=work/f'chr{c}.input'
                jobs.append(node(f'chr{c}.pgen-view',['python3',SELF,'pgen-view',prefix,view],
                                 [str(view)+e for e in ('.pgen','.pvar','.psam')],genotype_inputs(prefix)))
                cmd=common+['--step','2','--pgen',view,'--pred',str(step1)+'_pred.list','--minMAC',a.min_mac,'--out',dest]+extra
                if typ=='bt': cmd+=['--firth','--approx','--pThresh','0.01']
                if c=='X': cmd+=['--par-region','b'+a.grch]
                inputs+=[str(step1)+'_pred.list',str(step1)+'_1.loco']+[str(view)+e for e in ('.pgen','.pvar','.psam')]
            else:
                # Preserve imputed dosage: PGEN -> bgzipped VCF with DS, never hard-call BED.
                vcf=str(dest)+'.dosage.vcf.gz'
                jobs.append(node(f'chr{c}.export',[a.plink2,'--threads',a.threads,'--memory',a.memory]+pfile(prefix)+['--mac',a.min_mac,'--max-alleles','2','--export','vcf','bgz','vcf-dosage=DS-force','id-paste=iid','--out',str(dest)+'.dosage']+extra,[vcf],genotype_inputs(prefix)+([a.extract] if a.extract else [])))
                jobs.append(node(f'chr{c}.index',[tool('tabix'),'-f','-C','-p','vcf',vcf],[vcf+'.csi'],[vcf]))
                result=str(dest)+'.saige.txt'
                cmd=[a.saige_rscript,Path(a.saige_dir)/'step2_SPAtests.R','--vcfFile='+vcf,'--vcfFileIndex='+vcf+'.csi','--vcfField=DS','--chrom='+c,
                     '--GMMATmodelFile='+str(step1)+'.rda','--varianceRatioFile='+str(step1)+'.varianceRatio.txt',
                     '--SAIGEOutputFile='+result,'--LOCO=FALSE',f'--minMAC={a.min_mac}',
                     '--is_output_moreDetails=TRUE','--is_overwrite_output=TRUE']
                if typ=='bt': cmd+=['--is_Firth_beta=TRUE','--pCutoffforFirth=0.01']
                inputs+=[vcf,vcf+'.csi',str(step1)+'.rda',str(step1)+'.varianceRatio.txt']
            jobs.append(node(f'chr{c}.association',cmd,[result],inputs,(covars,method) if method!='saige' else None)); raw.append(result)
        final=root/(trait+'.'+method+'.gz')
        jobs.append(node('merge-results',['python3',HERE/'gwas_io.py',method,final]+raw,[final],raw+[str(HERE/'gwas_io.py')]))
        jobs.append(node('build-metadata',['python3',SELF,'build-metadata',final,a.grch],[str(final)+'.grch'],[final]))
        save_plan(work,jobs,a)

if __name__=='__main__':
    try:
        if len(sys.argv)>1 and sys.argv[1]=='execute-plan': execute(sys.argv[2])
        elif len(sys.argv)>1 and sys.argv[1]=='build-metadata': Path(sys.argv[2]+'.grch').write_text(sys.argv[3]+'\n')
        elif len(sys.argv)>1 and sys.argv[1]=='pgen-view':
            src,dst=sys.argv[2:]
            for ext in ('.pgen',):
                target=Path(dst+ext)
                if target.is_symlink(): target.unlink()
                elif target.exists(): raise FileExistsError(target)
                target.symlink_to(src+ext)
            with open(src+'.psam') as inp, open(dst+'.psam.tmp','w') as out:
                header=inp.readline()
                iid_only=header.startswith('#IID')
                if not iid_only and not header.startswith('#FID'): raise ValueError('Unsupported PSAM header')
                out.write(('#FID\t'+header.lstrip('#')) if iid_only else header)
                for line in inp: out.write(('0\t' if iid_only else '')+line)
            os.replace(dst+'.psam.tmp',dst+'.psam')
            with open(dst+'.pvar.tmp','wb') as out:
                if Path(src+'.pvar.zst').exists(): subprocess.run(['zstd','-dc',src+'.pvar.zst'],stdout=out,check=True)
                else:
                    with open(src+'.pvar','rb') as inp: shutil.copyfileobj(inp,out)
            os.replace(dst+'.pvar.tmp',dst+'.pvar')
        elif len(sys.argv)>1 and sys.argv[1]=='grm-links':
            actual,prefix=sys.argv[2:]
            for src,dst in [(actual,prefix+'.sparseGRM.mtx'),(actual+'.sampleIDs.txt',prefix+'.sparseGRM.mtx.sampleIDs.txt')]:
                p=Path(dst)
                if p.is_symlink(): p.unlink()
                elif p.exists(): raise FileExistsError('Refusing to replace existing GRM: '+dst)
                p.symlink_to(src)
        else:
            args=parse(); prep(args) if args.module=='prep_gwas' else association(args)
    except (ValueError,RuntimeError,FileNotFoundError,subprocess.CalledProcessError) as e:
        print('ERROR: '+str(e),file=sys.stderr); sys.exit(1)
