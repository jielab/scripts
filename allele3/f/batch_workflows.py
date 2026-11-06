#!/usr/bin/env python3
"""Plan allele3 work without opening large inputs; execute only via batch_runner.

The shell loads the repository's shared 0f/0os.sh before this module. This module
uses only the standard library during planning. Genomic reads belong to jobs.
"""
from __future__ import annotations

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from analysis_scope import STAGES, boolean_option, select_stages

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / 'f'


def absolute(path):
    return os.path.abspath(os.path.expanduser(str(path)))


def option(argv, *names, default=None):
    value = default
    for i, arg in enumerate(argv):
        if arg in names and i+1 < len(argv): value = argv[i+1]
        for name in names:
            if arg.startswith(name+'='): value = arg[len(name)+1:]
    return value


def remove_options(argv, names, flags=()):
    result=[]; i=0
    while i < len(argv):
        arg=argv[i]
        if arg in flags: i+=1; continue
        if arg in names:
            if i+1 >= len(argv): raise ValueError(arg+' needs a value')
            i+=2; continue
        if any(arg.startswith(n+'=') for n in names): i+=1; continue
        result.append(arg); i+=1
    return result


def env_paths():
    data=os.environ.get('DIR_DATA',configured_path('DIR_BIG'))
    work=os.environ.get('DIR_WORK',configured_path('DIR_DAT'))
    analysis=os.environ.get('ALLELE_ANALYSIS_ROOT',os.environ.get('ANALYSIS_ROOT',work+'/analysis')+'/allele3')
    return data,work,analysis


def worker_python(argv, mode):
    given=option(argv,'--python',default=os.environ.get('ALLELE_PYTHON'))
    if given: return given
    if mode=='hpc': return env_paths()[1]+'/.conda/envs/allele3-evolution/bin/python'
    candidates=[]
    if os.environ.get('CONDA_PREFIX') and os.environ.get('CONDA_DEFAULT_ENV')!='base':
        candidates.append(Path(os.environ['CONDA_PREFIX'])/'bin/python')
    candidates += [Path.home()/'anaconda3/envs/allele3-evolution/bin/python',
                   Path.home()/'miniconda3/envs/allele3-evolution/bin/python']
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate,os.X_OK): return str(candidate)
    return sys.executable


def worker_env(python, resource):
    # A saved .cmd may be submitted from another shell. Freeze the full shared
    # path table, including the super-roots and the markers used by 0os.sh,
    # instead of relying on bsub's ambient environment inheritance.
    prefixes=('ALLELE_','PCA3D_','HERV_','DIR_','OS_','_OS_DEFAULT_')
    exact={'DIR_WORK','DIR_DATA','DIR_SCRATCH','DIR_SOFTWARE','DATE_SCRATCH',
           'SHARED_PYTHON','ANALYSIS_ROOT','SCRATCH_ANALYSIS_ROOT','OS_FILES_ROOT',
           'OS_GU_CONDA','SCRIPTS_ROOT','ALLELE_OS_CONFIG','jobs_N'}
    excluded={'ALLELE_JOB_DIR','ALLELE_WRITE_CMD','ALLELE_RETRY_FAILED',
              'ALLELE_MAX_RUNNING','ALLELE_PLANNER_PYTHON'}
    result={k:v for k,v in os.environ.items() if (k in exact or k.startswith(prefixes)) and k not in excluded and not k.startswith('ALLELE_JOB_')}
    result.update(ALLELE_PYTHON=python, PCA3D_PYTHON=python, HERV_PYTHON=python,
                  ALLELE_MEMORY_LIMIT_GB=str(resource['memory_gib']),
                  PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
    if '/' in python: result['PATH']=str(Path(python).parent)+os.pathsep+os.environ.get('PATH','')
    return result


def resources(cores=1, memory=8, hours=24, queue='33'):
    r=dict(queue=os.environ.get('ALLELE_JOB_QUEUE',queue),
           cores=int(os.environ.get('ALLELE_JOB_CORES',cores)),
           memory_gib=int(os.environ.get('ALLELE_JOB_MEMORY_GIB',memory)),
           hours=int(os.environ.get('ALLELE_JOB_HOURS',hours)))
    if r['queue']=='short': r['cores']=max(r['cores'],40)
    if any(r[k]<1 for k in ('cores','memory_gib','hours')): raise ValueError('Job cores, GiB and hours must be positive integers')
    return r


def code_inputs():
    # Small source-code fingerprints are cheap; no checksum of the billion-row
    # source is computed on the login node.
    paths=[*ROOT.glob('*.sh'),ROOT/'allele3.env',*F.glob('*.py'),*F.glob('*.sh'),
           * (ROOT/'runtime').glob('*.py'),* (ROOT/'config').glob('*.tsv')]
    if os.environ.get('ALLELE_CONFIG'): paths.append(Path(os.environ['ALLELE_CONFIG']))
    if os.environ.get('ALLELE_OS_CONFIG'): paths.append(Path(os.environ['ALLELE_OS_CONFIG']))
    elif (ROOT.parent/'0f/0os.sh').is_file(): paths.append(ROOT.parent/'0f/0os.sh')
    return [dict(path=absolute(p),sha256=True) for p in sorted(set(paths))]


def base_plan(name, mode):
    return dict(name='allele3_'+name,jobs=[],inputs=code_inputs(),
        max_running=int(os.environ.get('ALLELE_MAX_RUNNING',32)),
        lsf_memory_scope=os.environ.get('ALLELE_LSF_MEMORY_SCOPE','host'),
        local_memory_gib=32, driver_python=sys.executable,
        provenance=dict(entry='batch_workflows.py',roots=dict(zip(('data','work','analysis'),env_paths())),
            scheduling='numeric done(jobid) for data dependencies; ended(jobid) only for concurrency lanes',
            source_reads='genomic parsing and annotation run only inside allocated workers in hpc mode'))


def add_job(plan,name,argv,python,deps=(),resource=None,outputs=()):
    resource=resource or resources()
    plan['jobs'].append(dict(id=name,argv=[str(a) for a in argv],deps=list(deps),cwd=str(ROOT),
        resources=resource,env=worker_env(python,resource),outputs=[absolute(p) for p in outputs]))


def selected_chromosomes(text):
    result=[]
    for token in text.split(','):
        token=token.strip().removeprefix('chr')
        if re.fullmatch(r'\d+-\d+',token):
            a,b=map(int,token.split('-'))
            if not 1<=a<=b<=22: raise ValueError('Invalid chromosome range')
            result += [str(i) for i in range(a,b+1)]
        else: result.append(token)
    if not result or len(set(result))!=len(result) or set(result)-set([str(i) for i in range(1,23)]+['X','Y','MT']):
        raise ValueError('Chromosomes must be unique 1..22, X, Y or MT')
    return result


def atomic_plan_file(path,value):
    from batch_runner import atomic_json
    path=Path(path)
    if path.exists():
        if json.loads(path.read_text())!=value: raise ValueError('Existing plan attachment differs: '+str(path))
    else: atomic_json(path,value)


def allofus_plan(args, mode, out, jobdir, python):
    from allofus_shards import make_plan
    plan=base_plan('allofus',mode)
    windows=make_plan(args.source,args.chromosomes,args.window_bp,args.source_version)
    wp=Path(jobdir)/'allofus_windows.json'
    prepared=out/'prepared'; normalized=out/'normalized'; summary=out/'summary'
    plan['inputs'].append(absolute(args.source))
    add_job(plan,'allofus_prepare',[python,F/'allofus_shards.py','prepare','--plan',wp,'--out',prepared],
            python,resource=resources(1,16,72),outputs=[prepared/'prepare_manifest.json'])
    if args.with_atlas:
        if not args.targets: raise ValueError('--with-atlas requires --targets')
        bind_atlas_inputs(plan,args)
        add_job(plan,'atlas_targets',atlas_argv(args,python,out/'atlas'),python,
                resource=resources(1,8,6),outputs=[out/'atlas/manifest.json'])
    ids=[]
    for shard in windows['shards']:
        name='aou_'+shard['id'];ids.append(name)
        add_job(plan,name,[python,F/'allofus_shards.py','normalize','--prepared',prepared,
                '--shard',shard['id'],'--out',normalized/shard['id']],python,
                deps=['allofus_prepare'],resource=resources(1,8,24),outputs=[normalized/shard['id']/'manifest.json'])
    # Limit each final dependency expression to 32 IDs. Large plans can have
    # hundreds/thousands of shards without exceeding an LSF expression limit.
    deps=barriers(plan,ids,python,out/'barriers')
    add_job(plan,'allofus_summary',[python,F/'allofus_shards.py','summarize','--plan',wp,
            '--normalized-root',normalized,'--out',summary],python,deps=deps,
            resource=resources(1,8,6),outputs=[summary/'manifest.json'])
    if args.targets:
        plan['inputs'].append(absolute(args.targets))
        add_job(plan,'allofus_query',[python,F/'allofus_shards.py','query','--prepared',prepared,
            '--normalized-root',normalized,'--targets',absolute(args.targets),'--out',out/'targets'],
            python,deps=['allofus_summary'],resource=resources(1,8,6),outputs=[out/'targets/manifest.json'])
        if args.with_atlas:
            add_job(plan,'allofus_atlas_join',[python,F/'atlas_allofus.py','--allofus-dir',out/'targets',
                '--atlas-dir',out/'atlas','--out',out/'functional'],python,
                deps=['allofus_query','atlas_targets'],resource=resources(1,4,2),outputs=[out/'functional/manifest.json'])
    elif args.with_atlas: raise ValueError('--with-atlas requires --targets')
    plan['provenance'].update(window_bp=args.window_bp,windows=len(windows['shards']),
        evidence='AF-only. No AC/AN/REF count is inferred; raw gzip is read once in the prepare job.')
    return plan,[(wp,windows)]


def barriers(plan,ids,python,directory):
    level=0
    while len(ids)>32:
        next_ids=[]
        for offset in range(0,len(ids),32):
            name=f'barrier_{level}_{offset//32:04d}'; receipt=Path(directory)/(name+'.json')
            add_job(plan,name,[python,F/'batch_workflows.py','worker','barrier','--out',receipt],python,
                    deps=ids[offset:offset+32],resource=resources(1,1,1),outputs=[receipt])
            next_ids.append(name)
        ids=next_ids;level+=1
    return ids


def atlas_argv(args,python,out):
    argv=[python,F/'atlas_scores.py','--targets',absolute(args.targets),'--build','GRCh38','--out',out,
          '--dir-data',env_paths()[0],'--source-version',args.atlas_source_version,'--max-targets',str(args.max_targets)]
    for attr,flag in (('avi_file','--avi-file'),('avi_index','--avi-index'),('atlas_resources','--resources')):
        value=getattr(args,attr,None)
        if value: argv += [flag,absolute(value)]
    return argv


def bind_atlas_inputs(plan,args):
    """Bind default and configured source identities without opening the atlas.

    Only the small resource descriptor is read on the submission host. Index
    discovery records both potential companions if none was explicitly chosen,
    so a later replacement or additional index changes the frozen plan.
    """
    from atlas_scores import load_resources
    descriptor=Path(args.atlas_resources or ROOT/'config/alphagenome_resources.tsv').expanduser().resolve()
    specs=load_resources(argparse.Namespace(resources=str(descriptor),dir_data=env_paths()[0],
        avi_file=args.avi_file,avi_index=args.avi_index,source_version=args.atlas_source_version))
    if not Path(args.targets).expanduser().is_file():
        raise ValueError('Required Atlas targets file missing: '+absolute(args.targets))
    plan['inputs'].append(dict(path=absolute(descriptor),sha256=True))
    for row in specs:
        if not row['enabled'] or not row['path']: continue
        path=row['path']
        plan['inputs'].append(path)
        explicit=args.avi_index if row['kind']=='avi' else None
        if row['index_path']: plan['inputs'].append(row['index_path'])
        if not explicit: plan['inputs'] += [path+'.tbi',path+'.csi']


def atlas_plan(args,mode,out,python):
    plan=base_plan('atlas',mode);plan['inputs'].append(absolute(args.targets))
    bind_atlas_inputs(plan,args)
    add_job(plan,'atlas_targets',atlas_argv(args,python,out),python,
            resource=resources(1,8,6),outputs=[out/'manifest.json'])
    return plan,[]


def mfi_plan(args,mode,out,python):
    plan=base_plan('ukb_mfi',mode);ids=[]
    chromosomes=selected_chromosomes(args.chromosomes)
    if set(chromosomes)-set(str(i) for i in range(1,23)): raise ValueError('MFI audit currently supports autosomes; X/Y require a separate ploidy-aware contract')
    for c in chromosomes:
        name='mfi_chr'+c;ids.append(name);target=out/('chr'+c)
        argv=[python,F/'ukb_mfi.py','--mfi-template',args.mfi_template,'--pvar-template',args.pvar_template,
              '--chromosomes',c,'--genome-build','GRCh37','--min-info',str(args.min_info),'--tabix','--out',target]
        if args.sort_inputs: argv+=['--sort-inputs']
        add_job(plan,name,argv,python,resource=resources(2,16,24),outputs=[target/'manifest.json'])
    add_job(plan,'mfi_summary',[python,F/'mfi_merge.py',
        '--directories-json',json.dumps([str(out/('chr'+c)) for c in chromosomes]),'--out',out/'summary'],
        python,deps=ids,resource=resources(1,4,2),outputs=[out/'summary/manifest.json'])
    return plan,[]


def herv_plan(args,mode,out,python):
    samples=[s.strip() for s in args.samples.split(',') if s.strip()]
    if not samples or len(samples)!=len(set(samples)) or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',s) for s in samples):
        raise ValueError('--samples must be distinct comma-separated sample IDs')
    plan=base_plan('herv_samples',mode)
    for sample in samples:
        resource=resources(16,64 if mode=='hpc' else 16,72)
        argv=['bash',ROOT/'2.herv.sh','--sample',sample,'--dat',args.dat,'--build',args.build,
              '--out',out/sample,'--replace',args.replace,'--jobs',str(min(8,resource['cores'])),
              '--cores',str(resource['cores']),'--python',python,'--memory-limit-gb',str(resource['memory_gib'])]
        for key,flag in (('dir_vcf','--dir-vcf'),('dir_pfile','--dir-pfile'),('fasta','--fasta'),('library','--library')):
            if getattr(args,key,None): argv += [flag,absolute(getattr(args,key))]
        add_job(plan,'herv_'+sample,argv,python,resource=resource,outputs=[out/sample/'job.state.json',out/sample/'genome_fraction.json'])
    # Personal HERV outputs retain their original census/checkpoint format.
    # Cohort inference must use research cohort's balanced plan and audit;
    # a submitted sample list is not labelled population representative here.
    plan['provenance']['sampling']='Explicit sample list; no population-representativeness claim'
    return plan,[]


def pipeline_plan(argv,mode,out,python):
    argv=remove_options(argv,('--python','--memory-limit-gb','--screen-log-interval'),('--verbose-screen',))
    if not argv or argv[0].startswith('-'): argv=['pipeline',*argv]
    requested=option(argv,'--stages')
    ibdmix=boolean_option(option(argv,'--ibdmix',default=False))
    names=select_stages(argv[0],requested,ibdmix)
    argv=remove_options(argv,('--ibdmix',))+['--ibdmix','TRUE' if ibdmix else 'FALSE']
    argv=remove_options(argv,('--out','--dir-out'))+['--out',str(out)]
    if not option(argv,'--metadata') and os.environ.get('ALLELE_METADATA'): argv+=['--metadata',os.environ['ALLELE_METADATA']]
    if not option(argv,'--coordinates') and os.environ.get('ALLELE_COORDINATES'): argv+=['--coordinates',os.environ['ALLELE_COORDINATES']]
    key=hashlib.sha256(json.dumps(dict(argv=argv,names=names,source=(F/'batch_pipeline.py').read_text()),sort_keys=True).encode()).hexdigest()
    plan=base_plan('pipeline',mode)
    # One command file and one allocation for every selected main module.
    # Native stage planning, audits, reuse and failure checks stay in compute.
    mode_args=['--all-stages'] if len(names)>1 else ['--stage',names[0]]
    cmd=[python,F/'batch_pipeline.py',*mode_args,'--selected',','.join(names),
         '--batch-id',key,'--args-json',json.dumps(argv,ensure_ascii=False)]
    add_job(plan,'all' if len(names)>1 else names[0],cmd,python,
            resource=resources(4,64 if mode=='hpc' else 16,72),
            outputs=[out/'.batch_stage_status'/key/(name+'.json') for name in names])
    plan['provenance'].update(selected=names,batch_id=key,ibdmix_enabled=ibdmix,
        disabled_stages=[] if ibdmix else ['archaic'],
        execution_order='sequential_success',scheduling='one job; selected modules run sequentially inside one allocation')
    return plan,[]


def single_plan(entry,argv,mode,python):
    name=Path(entry).name
    data,work,analysis=env_paths()
    analysis=option(argv,'--dir-analysis-root',default=analysis)
    if name=='install.sh':
        out=Path(absolute(analysis+'/environment_check'))
        plan=base_plan('environment_check',mode)
        cmd=remove_options(argv,('--python',))+['--python',python]
        add_job(plan,'environment_check',['bash',absolute(entry),*cmd],python,resource=resources(1,4,1))
        return plan,[],out
    explicit_out=option(argv,'--out','--dir-out','--store')
    action=argv[0] if argv and not argv[0].startswith('-') else 'run'
    # Keep the public shell's fresh-run protections before adding a resolved
    # --out to the worker. An invented default would otherwise satisfy that
    # shell guard and authorize --replace TRUE against a completed store.
    if name=='3.scan.sh' and action=='run' and not explicit_out:
        raise ValueError('Fresh scan requires --out DIR (or --store DIR); select the intended output explicitly')
    out=explicit_out
    if not out:
        if name=='1.pca3D.sh':
            out=os.environ.get('PCA3D_OUT')
            if action=='run' and not out:
                raise ValueError('PCA calculation requires --out DIR (or PCA3D_OUT); select the intended output explicitly')
            out=out or analysis+'/pca3d'
        elif name=='2.herv.sh':
            out=os.environ.get('HERV_OUT') or os.environ.get('ALLELE_HERV_OUT')
            sample=option(argv,'--sample',default=os.environ.get('HERV_SAMPLE','HG00096'))
            dataset=option(argv,'--dat','--gen',default=os.environ.get('HERV_GEN','1kg'))
            build=option(argv,'--build',default=os.environ.get('HERV_BUILD','38'))
            selected_root=option(argv,'--dir-analysis-root',default=os.environ.get('HERV_ANALYSIS_ROOT'))
            if not out:
                if selected_root: out=selected_root+'/herv/'+dataset+'/'+sample
                elif mode=='hpc': out=os.environ.get('DIR_SCRATCH',work)+'/allele3/herv/'+dataset+'/'+build+'/'+sample
                else: out=data+'/gen/'+dataset+'/'+build+'/herv/'+sample
        elif name=='3.scan.sh': out=os.environ.get('ALLELE_SCAN_STORE',analysis+'/stores/gnomad_scan')
        elif argv[:1]==['research']:
            out=analysis+'/research/'+(argv[1] if len(argv)>1 else 'unknown')
            if len(argv)>2 and argv[1]=='cohort' and argv[2] in {'plan','summarize'}: out+='/'+argv[2]
        else: out=os.environ.get('ALLELE_HERV_FOCUS_OUT') or analysis+'/final/herv_focus_20261028'
        argv += ['--out',out]
    out=Path(absolute(out))
    if name=='allele3.sh' and (not argv or argv[0]!='research'):
        analysis=Path(absolute(option(argv,'--dir-analysis-root',default=analysis)))
        if out==analysis: out=analysis/'final/herv_focus_20261028'
        return (*pipeline_plan(argv,mode,out,python),out)
    label=re.sub(r'[^A-Za-z0-9_]','_',name)
    plan=base_plan(label,mode)
    explicit_memory=option(argv,'--memory-limit-gb')
    default_cores=16 if name=='2.herv.sh' else 4
    requested_cores=int(option(argv,'--pca-threads','--threads','--cores',default=default_cores))
    resource=resources(requested_cores, int(explicit_memory) if explicit_memory else (64 if mode=='hpc' else 16),72)
    if explicit_memory and int(explicit_memory)>resource['memory_gib']: raise ValueError('--memory-limit-gb exceeds --job-memory-gb')
    argv=remove_options(argv,('--python','--memory-limit-gb'))
    argv+=['--python',python,'--memory-limit-gb',str(resource['memory_gib'])]
    add_job(plan,label,['bash',absolute(entry),*argv],python,resource=resource)
    return plan,[],out


def batch_parser():
    data,work,analysis=env_paths()
    p=argparse.ArgumentParser(description=__doc__,epilog='Execution flags: --submit local|hpc --jobs-N 1..8 --max-running 200 --write-cmd --wait TRUE|FALSE --date-scratch YYYY-MM-DD --job-dir DIR --queue NAME --job-cores N --job-memory-gb N --job-hours N --retry-failed')
    sub=p.add_subparsers(dest='command',required=True)
    def shared(q):
        q.add_argument('--out');q.add_argument('--python');q.add_argument('--dry-run',action='store_true')
    def atlas_opts(q):
        q.add_argument('--avi-file');q.add_argument('--avi-index');q.add_argument('--atlas-resources')
        q.add_argument('--atlas-source-version',default='unverified')
        q.add_argument('--max-targets',type=int,default=5000)
    q=sub.add_parser('allofus',help='One gzip preparation job -> window jobs -> audited summary / optional targets + Atlas')
    shared(q);atlas_opts(q)
    q.add_argument('--source',default=os.environ.get('ALLELE_ALL_OF_US',data+'/annot/allele3/allofus/hg38_allofus.txt.gz'))
    q.add_argument('--chromosomes',default='1-22,X,Y,MT');q.add_argument('--window-bp',type=int,default=10000000)
    q.add_argument('--source-version',default='unverified');q.add_argument('--targets');q.add_argument('--with-atlas',action='store_true')
    q=sub.add_parser('atlas',help='Targeted local AVI Tabix annotation; CPU only')
    shared(q);atlas_opts(q);q.add_argument('--targets',required=True)
    q=sub.add_parser('mfi',help='One exact PVAR/MFI quality audit per autosome')
    shared(q);q.add_argument('--chromosomes',default='1-22')
    q.add_argument('--mfi-template',default=data+'/annot/allele3/ukb_imp_mfi/ukb_mfi_chr{chrom}_v3.txt')
    q.add_argument('--pvar-template',default=data+'/gen/ukb/37/imp/chr{chrom}.pvar')
    q.add_argument('--min-info',type=float,default=.8);q.add_argument('--sort-inputs',action='store_true')
    q=sub.add_parser('herv',help='One individual HERV job per explicit sample ID')
    shared(q);q.add_argument('--samples',required=True);q.add_argument('--dat',default='1kg')
    q.add_argument('--build',choices=('37','38'),default='38');q.add_argument('--replace',choices=('TRUE','FALSE'),default='FALSE')
    for opt in ('dir-vcf','dir-pfile','fasta','library'):q.add_argument('--'+opt)
    q=sub.add_parser('pipeline',help='Stage DAG with independent stage receipts and final report')
    shared(q);q.add_argument('arguments',nargs=argparse.REMAINDER,help='Use -- before ordinary pipeline options')
    return p


def worker(argv):
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('barrier');q.add_argument('--out',required=True)
    q=sub.add_parser('manifest-summary');q.add_argument('--directories-json',required=True);q.add_argument('--out',required=True)
    a=p.parse_args(argv)
    if os.environ.get('ALLELE_EXECUTION_MODE')=='hpc' and not os.environ.get('LSB_JOBID'):
        raise RuntimeError('HPC worker requires an LSF allocation')
    from batch_runner import atomic_json
    if a.action=='barrier': atomic_json(a.out,dict(status='complete',lsf_job_id=os.environ.get('LSB_JOBID')));return 0
    from analysis_integrity import validate_manifest,transactional,identity,InputPath
    directories=json.loads(a.directories_json)
    @transactional('manifest.json',resolve_inputs=lambda raw:[InputPath(p,'directory') for p in raw['directories']],input_prefixes=())
    def summarize(directories,out):
        rows=[]
        for directory in directories:
            m=validate_manifest(directory,require_evidence=True)
            if m.get('status') not in {'complete','no_usable_exact_high_quality_contrast'}: raise ValueError('Incomplete quality audit: '+directory)
            rows.append(dict(directory=directory,status=m['status'],counts=m.get('counts',{}),input_manifest=identity(Path(directory)/'manifest.json',full=True)))
        Path(out,'cohorts.json').write_text(json.dumps(rows,indent=2)+'\n')
        return dict(status='complete',sources=rows,evidence='Quality audit only; no WGS or genotype counts inferred')
    summarize(directories,a.out);return 0


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if argv[:1]==['worker']: return worker(argv[1:])
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--entry',required=True);parser.add_argument('--mode',choices=('local','hpc'),required=True)
    parser.add_argument('payload',nargs=argparse.REMAINDER)
    outer=parser.parse_args(argv)
    args=outer.payload[1:] if outer.payload[:1]==['--'] else outer.payload
    is_batch=args[:1]==['batch']
    py=worker_python(args,outer.mode)
    if is_batch:
        a=batch_parser().parse_args(args[1:])
        default=env_paths()[2]+'/research/'+a.command
        if a.command=='allofus': default=env_paths()[0]+'/annot/allele3/allofus/cache_hg38'
        rest=[]
        if a.command=='pipeline':
            rest=a.arguments[1:] if a.arguments[:1]==['--'] else a.arguments
            default=option(rest,'--out','--dir-out',default=default)
        out=Path(absolute(a.out or default))
        if a.command=='pipeline':
            analysis=Path(absolute(option(rest,'--dir-analysis-root',default=env_paths()[2])))
            if out==analysis:out=analysis/'final/herv_focus_20261028'
        jobdir=Path(absolute(os.environ.get('ALLELE_JOB_DIR',str(out)+'.jobs')))
        if a.command=='allofus': plan,attachments=allofus_plan(a,outer.mode,out,jobdir,py)
        elif a.command=='atlas': plan,attachments=atlas_plan(a,outer.mode,out,py)
        elif a.command=='mfi': plan,attachments=mfi_plan(a,outer.mode,out,py)
        elif a.command=='herv': plan,attachments=herv_plan(a,outer.mode,out,py)
        else:
            plan,attachments=pipeline_plan(['pipeline',*rest],outer.mode,out,py)
        dry=a.dry_run
    else:
        dry='--dry-run' in args or '--preflight' in args
        args=[x for x in args if x not in {'--dry-run','--preflight'}]
        plan,attachments,out=single_plan(outer.entry,args,outer.mode,py)
        jobdir=Path(absolute(os.environ.get('ALLELE_JOB_DIR',str(out)+'.jobs')))
    if dry:
        print(json.dumps(dict(status='plan_only',mode=outer.mode,job_directory=str(jobdir),
            output=str(out),plan=plan,attachments=[str(p) for p,v in attachments]),indent=2))
        return 0
    for path,value in attachments: atomic_plan_file(path,value)
    # Bind attachments too: workers must not consume a changed shard plan.
    plan['inputs'] += [dict(path=str(path),sha256=True) for path,value in attachments]
    from batch_runner import run_plan
    print(f"[batch] mode={outer.mode} jobs={len(plan['jobs'])} max_running={plan['max_running']} job_dir={jobdir} output={out}",flush=True)
    result=run_plan(plan,jobdir,mode=outer.mode,write_only=os.environ.get('ALLELE_WRITE_CMD')=='1',
        jobs_n=int(os.environ.get('jobs_N',1)),wait=os.environ.get('SUBMIT_WAIT','FALSE')=='TRUE',
        retry_failed=os.environ.get('ALLELE_RETRY_FAILED')=='1')
    print(json.dumps(result,indent=2,default=str))
    return 0


if __name__=='__main__':
    try: raise SystemExit(main())
    except (ValueError,OSError,RuntimeError,KeyError) as exc: raise SystemExit('batch-workflows: '+str(exc))
