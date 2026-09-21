#!/usr/bin/env python3
"""Run the user's bulik/ldsc checkout in its Python 2 conda environment."""
import argparse
import gzip
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import json
import csv
import hashlib
import tempfile
import fcntl
from contextlib import contextmanager
import sys
from collections import Counter

DEFAULT_REF='/mnt/e/refLD/ldsc/1000G/1000G_Phase3_ldscores/LDscore.'
DEFAULT_WEIGHTS='/mnt/e/refLD/ldsc/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'
DEFAULT_ALLELES='/mnt/e/refLD/ldsc/hm3/w_hm3.snplist'
AUTOSOMES={str(c) for c in range(1,23)}

def validate_references(ref, weights):
    paths=[]
    for c in range(1,23):
        for prefix,suffix in ((ref,'.l2.ldscore.gz'),(ref,'.l2.M_5_50'),(weights,'.l2.ldscore.gz')):
            stem=prefix.replace('@',str(c)) if '@' in prefix else prefix+str(c)
            f=Path(stem+suffix)
            if not f.is_file() or not f.stat().st_size: raise ValueError('Missing LDSC reference: '+str(f))
            paths.append(f)
    return paths

def boolean(x):
    if x.upper() not in ('TRUE','FALSE'): raise argparse.ArgumentTypeError('Use TRUE/FALSE')
    return x.upper() == 'TRUE'

@contextmanager
def directory_lock(path):
    fd=os.open(path,os.O_RDONLY)
    try:
        fcntl.flock(fd,fcntl.LOCK_EX)
        yield
    finally:
        os.close(fd)


def cache_paths(source):
    source=Path(source)
    trait=source.name.removesuffix('.gz')
    qc=source.parent.parent/'qc' if source.parent.name=='gwas' else source.parent/'qc'
    return (source.parent/(trait+'.sumstats.gz'), qc/(trait+'.sumstats.cache.tsv'),
            qc/(trait+'.sumstats.chromosomes.tsv'))


def file_identity(path):
    path=Path(path).resolve(); st=path.stat()
    return [str(path),st.st_size,st.st_mtime_ns]


def validate_sumstats(cache):
    if not cache.is_file():
        raise ValueError('Missing pre-generated sumstats: '+str(cache))
    with gzip.open(cache,'rt') as handle:
        if not {'SNP','A1','A2','Z','N'}.issubset(handle.readline().split()) or not handle.readline().strip():
            raise ValueError('Empty/invalid pre-generated sumstats: '+str(cache))


def prepare(source, merge_alleles, software, interpreter, fallback_n=None, allow_build='True'):
    source=Path(source).resolve()
    cache,meta,auditfile=cache_paths(source)
    if allow_build != 'True':
        validate_sumstats(cache)
        print('Reuse pre-generated sumstats: '+str(cache),flush=True)
        return
    meta.parent.mkdir(parents=True,exist_ok=True)
    # Lock the directory inode, avoiding persistent per-trait lock files.
    with directory_lock(source.parent):
        signature=hashlib.sha256(json.dumps(dict(version=2,source=file_identity(source),
            alleles=file_identity(merge_alleles),munge=file_identity(Path(software)/'munge_sumstats.py'),
            interpreter=interpreter,N=fallback_n),sort_keys=True).encode()).hexdigest()
        prior={}
        if meta.exists():
            with meta.open() as handle: prior=dict(csv.reader(handle,delimiter='\t'))
        if (cache.is_file() and auditfile.is_file() and prior.get('signature')==signature
                and prior.get('output')==str(file_identity(cache))):
            print('Reuse munged GWAS: '+str(cache),flush=True); return
        print('Prepare munged GWAS: '+str(cache),flush=True)
        with tempfile.TemporaryDirectory(prefix='.ldsc-',dir=source.parent) as work:
            prefix=Path(work)/source.name.removesuffix('.gz'); fixed=str(prefix)+'.input.gz'
            fix_p(str(source),fixed,merge_alleles)
            audit=json.loads(Path(fixed+'.chromosomes.json').read_text())
            with gzip.open(fixed,'rt') as handle: header=handle.readline().split()
            if 'N' not in header and fallback_n is None: raise ValueError('Missing N: '+str(source))
            n=['--N-col','N'] if 'N' in header else ['--N',str(fallback_n)]
            cmd=interpreter+[str(Path(software)/'munge_sumstats.py'),'--sumstats',fixed,
                '--out',str(prefix),'--merge-alleles',merge_alleles,'--snp','SNP','--a1','EA',
                '--a2','NEA','--signed-sumstats','BETA,0','--p','P','--chunksize','10000']+n
            log=meta.parent/(source.name.removesuffix('.gz')+'.sumstats.log')
            with log.open('w') as handle:
                handle.write(audit['reason']+'\n'); handle.flush()
                subprocess.run(cmd,stdout=handle,stderr=subprocess.STDOUT,check=True)
            munged=Path(str(prefix)+'.sumstats.gz')
            with gzip.open(munged,'rt') as handle:
                if not {'SNP','A1','A2','Z','N'}.issubset(handle.readline().split()) or not handle.readline():
                    raise ValueError('Empty/invalid munged GWAS: '+str(munged))
            with auditfile.open('w') as handle:
                writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
                writer.writerow(['CHR','INPUT','KEPT','REASON'])
                for chrom,count in audit['input'].items():
                    writer.writerow([chrom,count,audit['kept'].get(chrom,0),
                        audit['reason'] if chrom not in AUTOSOMES else 'HM3 and valid P filter before munge QC'])
            os.replace(munged,cache)
            tmp=meta.with_suffix('.tmp')
            with tmp.open('w') as handle:
                csv.writer(handle,delimiter='\t',lineterminator='\n').writerows(
                    [('signature',signature),('output',str(file_identity(cache)))])
            os.replace(tmp,meta)


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''Examples (replace /path/... with your files):
  gwas_compare.sh ldsc --gwas-files /path/A.gz,/path/B.gz --output-dir /path/ldsc-results \\
    --merge-alleles /path/w_hm3.snplist \\
    --ref-ld-chr /path/eur_w_ld_chr/ --w-ld-chr /path/eur_w_ld_chr/
  Add --run FALSE to write commands only, --run-rg FALSE for h2 only,
  or --N 10000 if the actual sample size is 10000 and no N column is present.
''')
    p.add_argument('--gwas-files',required=True,help='Comma-separated standardized GWAS files')
    p.add_argument('--output-dir',default='/mnt/d/analysis/gwas/ldsc')
    p.add_argument('--ldsc-software-dir',default='/mnt/d/software/ldsc')
    p.add_argument('--conda',default=shutil.which('conda') or str(Path.home()/'anaconda3/bin/conda'))
    p.add_argument('--conda-env',default='ldsc')
    p.add_argument('--python',help='Explicit LDSC-compatible interpreter; bypass conda')
    p.add_argument('--merge-alleles',default=DEFAULT_ALLELES,help='w_hm3.snplist')
    p.add_argument('--ref-ld-chr',default=DEFAULT_REF,help='Reference LD-score prefix, including trailing / or dot')
    p.add_argument('--w-ld-chr',default=DEFAULT_WEIGHTS,help='Regression-weight prefix')
    p.add_argument('--N',type=float,help='Explicit fallback only when N is absent')
    p.add_argument('--missing-n',choices=('error','skip'),default='error',help='skip writes an explicit unavailable status; never invents N')
    p.add_argument('--run-munge',type=boolean,default=True)
    p.add_argument('--run-h2',type=boolean,default=True)
    p.add_argument('--run-rg',type=boolean,default=True)
    p.add_argument('--run',type=boolean,default=True,help='FALSE writes commands only')
    a=p.parse_args()
    files=[Path(x.strip()).resolve() for x in a.gwas_files.split(',')]
    if any(not f.is_file() for f in files): p.error('Input file missing')
    if a.N is not None and a.N<=0: p.error('--N must be positive')
    if a.run_munge and not Path(a.merge_alleles).is_file(): p.error('Missing --merge-alleles file')
    validate_references(a.ref_ld_chr,a.w_ld_chr)
    out=Path(a.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    python=[a.python] if a.python else [a.conda,'run','--no-capture-output','-n',a.conda_env,'python']
    script=['#!/usr/bin/env bash','set -euo pipefail',
            shlex.join(['rm','-f','--']+[str(out/name) for name in ('h2.tsv','rg.tsv','rg.png')])]
    names=[];statuses=[]
    labels=[f.name.removesuffix('.gz') for f in files]
    if len(set(labels)) != len(labels): p.error('Trait names must be unique')
    helper=str(Path(__file__).resolve())
    for folder in ('h2.log','rg.log'):
        (out/folder).mkdir(exist_ok=True)
    for f,trait in zip(files,labels):
        cache,_,_=cache_paths(f)
        if not a.run_munge:
            validate_sumstats(cache)
            names.append((trait,str(cache)))
            statuses.append((str(f),'SCHEDULED','Pre-generated adjacent sumstats; standard LDSC/reference covers 1-22 only'))
            continue
        with (gzip.open(f,'rt') if f.suffix=='.gz' else f.open()) as handle:
            header=handle.readline().strip().split()
        if 'N' not in header and a.N is None:
            if a.missing_n!='skip': p.error(str(f)+': missing N; provide --N (no assumed sample size)')
            statuses.append((str(f),'UNAVAILABLE','Missing N: supply a verified sample size; excluded from h2 and rg'))
            print('UNAVAILABLE LDSC (missing N): '+str(f),flush=True)
            continue
        if not set(('SNP','CHR','EA','NEA','BETA','P')).issubset(header):
            p.error(str(f)+': need SNP CHR EA NEA BETA P')
        statuses.append((str(f),'SCHEDULED','Standard LDSC/reference covers 1-22 only; X/Y/MT exclusions recorded in trait qc chromosome audit'))
        cache,_,_=cache_paths(f)
        names.append((trait,str(cache)))
        cmd=[sys.executable,helper,'prepare','--source',str(f),'--merge-alleles',a.merge_alleles,
             '--software',a.ldsc_software_dir,'--allow-build',str(a.run_munge),
             '--interpreter']+python
        if a.N is not None: cmd[3:3]=['--fallback-n',str(a.N)]
        script.append(shlex.join(cmd))
    ref=['--ref-ld-chr',a.ref_ld_chr,'--w-ld-chr',a.w_ld_chr]
    ldsc=python+[str(Path(a.ldsc_software_dir)/'ldsc.py')]
    if a.run_h2:
        for trait,f in names:
            script.append(shlex.join(ldsc+['--h2',f,'--out',str(out/'h2.log'/(trait+'.h2'))]+ref))
    if a.run_rg:
        for i in range(len(names)-1):
            script.append(shlex.join(ldsc+['--rg',','.join(f for _,f in names[i:]),
                                         '--out',str(out/'rg.log'/(names[i][0]+'.rg'))]+ref))
    script.append(shlex.join([sys.executable,str(Path(__file__).with_name('gwas_ldsc_summary.py')),
                             '--output-dir',str(out),'--h2',str(a.run_h2),'--rg',str(a.run_rg)]))
    cmdfile=out/'ldsc.cmd.sh'; cmdfile.write_text('\n'.join(script)+'\n')
    (out/'inputs.status.tsv').write_text('FILE\tSTATUS\tREASON\n'+'\n'.join('\t'.join(row) for row in statuses)+'\n')
    if not names: p.error('No GWAS with usable N; see inputs.status.tsv')
    print(cmdfile,flush=True)
    if a.run: subprocess.run(['bash',str(cmdfile)],check=True)

def fix_p(src,dst,merge_alleles=None):
    import math
    tmp=dst+'.tmp'
    chromosomes=Counter(); kept=Counter(); invalid=0
    allowed=None
    if merge_alleles:
        with open(merge_alleles) as inp: allowed={line.split()[0] for line in inp if line.strip()}
    with (gzip.open(src,'rt') if src.endswith('.gz') else open(src)) as inp, gzip.open(tmp,'wt') as out:
        header=inp.readline().split(); pi=header.index('P'); ci=header.index('CHR'); si=header.index('SNP'); out.write('\t'.join(header)+'\n')
        for line in inp:
            row=line.split()
            if len(row)!=len(header): raise ValueError('Malformed GWAS row')
            chr=row[ci].upper().removeprefix('CHR')
            chr={'X':'23','Y':'24','MT':'25','M':'25'}.get(chr,chr)
            chromosomes[chr]+=1
            # Standard LDSC is autosomal. Audit omitted X/Y explicitly.
            if chr not in AUTOSOMES: continue
            if allowed is not None and row[si] not in allowed: continue
            try: pv=float(row[pi])
            except ValueError: invalid+=1;continue
            if not math.isfinite(pv) or pv<0 or pv>1: invalid+=1;continue
            if pv<1e-300: row[pi]='1e-300'
            out.write('\t'.join(row)+'\n')
            kept[chr]+=1
    os.replace(tmp,dst)
    Path(dst+'.chromosomes.json').write_text(json.dumps(dict(input=dict(chromosomes),kept=dict(kept),invalid_p=invalid,scope='autosomes',reason='Standard LDSC and installed LD scores cover chromosomes 1-22; X/Y/MT are not estimated'),indent=2))

if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='fix-p': fix_p(*sys.argv[2:])
    elif len(sys.argv)>1 and sys.argv[1]=='prepare':
        parser=argparse.ArgumentParser()
        parser.add_argument('--source',required=True)
        parser.add_argument('--merge-alleles',required=True)
        parser.add_argument('--software',required=True)
        parser.add_argument('--fallback-n',type=float)
        parser.add_argument('--allow-build',choices=('True','False'),default='True')
        parser.add_argument('--interpreter',nargs=argparse.REMAINDER,required=True)
        args=parser.parse_args(sys.argv[2:])
        prepare(**vars(args))
    else: main()
