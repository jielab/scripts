#!/usr/bin/env python3
"""Run the user's bulik/ldsc checkout in its Python 2 conda environment."""
import argparse
import gzip
import os
from pathlib import Path
import shlex
import shutil
import subprocess

def boolean(x):
    if x.upper() not in ('TRUE','FALSE'): raise argparse.ArgumentTypeError('Use TRUE/FALSE')
    return x.upper() == 'TRUE'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gwas-files',required=True,help='Comma-separated standardized GWAS files')
    p.add_argument('--output-dir',default='ldsc')
    p.add_argument('--ldsc-software-dir',default='/mnt/d/software/ldsc')
    p.add_argument('--conda',default=shutil.which('conda') or str(Path.home()/'anaconda3/bin/conda'))
    p.add_argument('--conda-env',default='ldsc')
    p.add_argument('--python',help='Explicit LDSC-compatible interpreter; bypass conda')
    p.add_argument('--merge-alleles',required=True,help='w_hm3.snplist')
    p.add_argument('--ref-ld-chr',required=True,help='Reference LD-score prefix, including trailing / or dot')
    p.add_argument('--w-ld-chr',required=True,help='Regression-weight prefix')
    p.add_argument('--N',type=float,help='Explicit fallback only when N is absent')
    p.add_argument('--run-munge',type=boolean,default=True)
    p.add_argument('--run-h2',type=boolean,default=True)
    p.add_argument('--run-rg',type=boolean,default=True)
    p.add_argument('--run',type=boolean,default=True,help='FALSE writes commands only')
    a=p.parse_args()
    files=[Path(x.strip()).resolve() for x in a.gwas_files.split(',')]
    if any(not f.is_file() for f in files): p.error('Input file missing')
    if a.N is not None and a.N<=0: p.error('--N must be positive')
    if not Path(a.merge_alleles).is_file(): p.error('Missing --merge-alleles file')
    out=Path(a.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    python=[a.python] if a.python else [a.conda,'run','--no-capture-output','-n',a.conda_env,'python']
    script=['#!/usr/bin/env bash','set -euo pipefail']
    names=[]
    for i,f in enumerate(files,1):
        name=out/(str(i).zfill(2)+'.'+f.name.removesuffix('.gz'))
        names.append(str(name)+'.sumstats.gz')
        if not a.run_munge:
            if a.run and not Path(names[-1]).is_file(): p.error('Missing munged file: '+names[-1])
            continue
        with (gzip.open(f,'rt') if f.suffix=='.gz' else f.open()) as handle:
            header=handle.readline().strip().split()
        if not set(('SNP','EA','NEA','BETA','P')).issubset(header): p.error(str(f)+': need SNP EA NEA BETA P')
        n=['--N-col','N'] if 'N' in header else ['--N',str(a.N)]
        if 'N' not in header and a.N is None: p.error(str(f)+': missing N; provide --N (no assumed sample size)')
        # Original GWAS is never modified. munge rejects P=0, so stream a capped copy.
        fixed=str(name)+'.input.gz'
        helper=str(Path(__file__).resolve())
        script.append(shlex.join(['python3',helper,'fix-p',str(f),fixed]))
        cmd=python+[str(Path(a.ldsc_software_dir)/'munge_sumstats.py'),'--sumstats',fixed,
                    '--out',str(name),'--merge-alleles',a.merge_alleles,'--snp','SNP','--a1','EA',
                    '--a2','NEA','--signed-sumstats','BETA,0','--p','P','--chunksize','10000']+n
        script.append(shlex.join(cmd))
    ref=['--ref-ld-chr',a.ref_ld_chr,'--w-ld-chr',a.w_ld_chr]
    ldsc=python+[str(Path(a.ldsc_software_dir)/'ldsc.py')]
    if a.run_h2:
        for f in names: script.append(shlex.join(ldsc+['--h2',f,'--out',f.removesuffix('.sumstats.gz')+'.h2']+ref))
    if a.run_rg:
        for i in range(len(names)-1):
            script.append(shlex.join(ldsc+['--rg',','.join(names[i:]),'--out',str(out/f'{i+1:02d}.rg')]+ref))
    cmdfile=out/'ldsc.cmd.sh'; cmdfile.write_text('\n'.join(script)+'\n')
    print(cmdfile,flush=True)
    if a.run: subprocess.run(['bash',str(cmdfile)],check=True)

def fix_p(src,dst):
    import math
    tmp=dst+'.tmp'
    with (gzip.open(src,'rt') if src.endswith('.gz') else open(src)) as inp, gzip.open(tmp,'wt') as out:
        header=inp.readline().split(); pi=header.index('P'); out.write('\t'.join(header)+'\n')
        for line in inp:
            row=line.split()
            if len(row)!=len(header): raise ValueError('Malformed GWAS row')
            try: pv=float(row[pi])
            except ValueError: continue
            if not math.isfinite(pv) or pv<0 or pv>1: continue
            if pv<1e-300: row[pi]='1e-300'
            out.write('\t'.join(row)+'\n')
    os.replace(tmp,dst)

if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='fix-p': fix_p(*sys.argv[2:])
    else: main()
