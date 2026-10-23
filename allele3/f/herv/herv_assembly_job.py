#!/usr/bin/env python3
"""Run/resume individual autosomal ERV annotation with reference and both haplotypes."""
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

from herv_assembly import save, now


def positive(value):
    number=int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('Must be a positive integer')
    return number


def component(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',value):
        raise argparse.ArgumentTypeError('Use a single name containing letters, digits, dots, underscores or hyphens')
    return value


def parser():
    p=argparse.ArgumentParser(description=__doc__,allow_abbrev=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='Run: ./1.herv.sh. Advanced annotation commands: '
               'python f/herv/herv_assembly.py {prepare,run,summarize} --help.')
    p.add_argument('--out','--dir-out',type=Path,
        help='Override results directory; otherwise /mnt/d/analysis/allele3/herv/<gen>/<sample>')
    p.add_argument('--fasta',type=Path,default=Path('/mnt/f/gen/1hgp/38/GRCH38.fasta.gz'),help='Indexed reference FASTA')
    p.add_argument('--library',type=Path,default=Path('/mnt/f/annot/dfam_4.0_human/human.curated.fa'),help='Repeat library')
    p.add_argument('--repeatmasker',type=Path,
        help='RepeatMasker executable; reuse checkpoint path, then discover allele3-repeatmasker or legacy metag-repeatmasker')
    p.add_argument('--gen',type=component,default='1kg',help='Dataset name used in the default output directory')
    p.add_argument('--sample',type=component,default='HG00096',help='VCF sample or PGEN IID; scope is chr1-22')
    p.add_argument('--input-format',choices=['auto','vcf','pgen'],default='auto',
        help='Reuse the prepared source; for a new sample auto prefers a complete PGEN fileset and available reader')
    p.add_argument('--vcf-dir','--dir-vcf',type=Path,default=Path('/mnt/f/gen/1kg/38/vcf'),help='VCFs used if preparation is missing')
    p.add_argument('--pfile-dir','--dir-pfile',type=Path,default=Path('/mnt/f/gen/1kg/38/pfile'),help='Chromosome PGEN/PVAR/PSAM filesets')
    p.add_argument('--pgen-python',type=Path,help='Python with numpy, pgenlib and zstandard; autodetected when omitted')
    p.add_argument('--jobs',type=positive,default=8,help='Concurrent RepeatMasker chunks')
    p.add_argument('--cores',type=positive,default=16,help='Limit CPU affinity to the first available CPUs')
    p.add_argument('--chunk-bp',type=positive,default=5000000,help='Reference core chunk length')
    p.add_argument('--halo',type=int,default=20000,help='Minimum flanking context in bp')
    p.add_argument('--scratch','--dir-scratch',type=Path,default=Path('/tmp/allele3-herv-assembly'),help='Temporary chunk directory')
    p.add_argument('--dry-run',action='store_true',help='Print commands without writing files or running analysis')
    return p


def pgen_runtime(explicit=None):
    if explicit:
        candidates=[str(explicit)]
    elif os.environ.get('HERV_PGEN_PYTHON') or os.environ.get('ALLELE_HERV_PGEN_PYTHON') or os.environ.get('METAG_PGEN_PYTHON'):
        candidates=[os.environ.get('HERV_PGEN_PYTHON') or os.environ.get('ALLELE_HERV_PGEN_PYTHON') or os.environ['METAG_PGEN_PYTHON']]
    else:
        root=conda_root()
        candidates=[sys.executable,str(root/'envs/herv-pgen/bin/python'),str(root/'envs/allele3-herv-pgen/bin/python'),
                    str(root/'envs/allele3-evolution/bin/python'),str(root/'envs/metag-herv-pgen/bin/python')]
    for candidate in dict.fromkeys(candidates):
        if not Path(candidate).is_file():
            continue
        try:
            subprocess.run([candidate,'-c',
                'import numpy,pgenlib,zstandard; '
                'assert hasattr(pgenlib,"PvarReader"); '
                'assert hasattr(pgenlib.PgenReader,"read_alleles_and_phasepresent_range")'],
                check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=15)
            return candidate
        except (OSError,subprocess.SubprocessError):
            continue
    return None


def pgen_available(root):
    return all((root/f'chr{i}.pgen').is_file() and (root/f'chr{i}.psam').is_file()
        and ((root/f'chr{i}.pvar.zst').is_file() or (root/f'chr{i}.pvar').is_file()) for i in range(1,23))


def conda_root():
    return Path(os.environ.get('HERV_CONDA_ROOT') or os.environ.get('ALLELE_CONDA_ROOT') or
                os.environ.get('METAG_CONDA_ROOT') or str(Path.home()/'anaconda3'))


def repeatmasker_path(args):
    if args.repeatmasker:
        return args.repeatmasker
    if os.environ.get('HERV_REPEATMASKER') or os.environ.get('ALLELE_HERV_REPEATMASKER'):
        return Path(os.environ.get('HERV_REPEATMASKER') or os.environ['ALLELE_HERV_REPEATMASKER'])
    # Executable paths are part of the existing run signature. Moving results
    # must not silently select a different environment and invalidate receipts.
    checkpoint=args.out/'run.json'
    if checkpoint.is_file():
        return Path(json.loads(checkpoint.read_text())['settings']['repeatmasker'])
    installed = shutil.which('RepeatMasker')
    if installed:
        return Path(installed)
    candidates=[conda_root()/'envs'/env/'bin/RepeatMasker'
                for env in ['herv-repeatmasker','allele3-repeatmasker','metag-repeatmasker']]
    return next((p for p in candidates if p.is_file()),candidates[0])


def commands(args):
    engine=[sys.executable,str(Path(__file__).with_name('herv_assembly.py'))]
    common=['--out',str(args.out),'--fasta',str(args.fasta.resolve())]
    prepared=args.out/'prepared.json'
    steps=[]
    if prepared.exists():
        saved=json.loads(prepared.read_text())
        if saved['sample'] != args.sample or saved['chromosomes'] != ['chr'+str(i) for i in range(1,23)]:
            raise ValueError('Prepared sample/scope differs from requested sample and chr1-22; use a different --out')
        args.effective_input=saved.get('input_source','vcf')
        if args.input_format not in ('auto',args.effective_input):
            raise ValueError('Prepared input format differs; use a separate --out for a PGEN/VCF comparison')
    else:
        complete=pgen_available(args.pfile_dir)
        runtime=pgen_runtime(args.pgen_python) if args.input_format!='vcf' and complete else None
        if args.input_format=='pgen' and (not complete or not runtime):
            raise ValueError('PGEN requires chr1-22 .pgen/.pvar[.zst]/.psam and a --pgen-python '
                             'with numpy, pgenlib and zstandard')
        args.effective_input='pgen' if args.input_format!='vcf' and complete and runtime else 'vcf'
        if args.effective_input=='pgen':
            steps.append([runtime,str(Path(__file__).with_name('herv_pgen.py')),'prepare']+common+
                ['--pfile-dir',str(args.pfile_dir.resolve()),'--sample',args.sample,'--jobs','2',
                 '--scratch',str(args.scratch.resolve()/'pgen-prepare')])
        else:
            steps.append(engine+['prepare']+common+['--vcf-dir',str(args.vcf_dir.resolve()),
                '--sample',args.sample,'--jobs','2','--vcf-scratch',str(args.scratch.resolve()/'vcf-prepare')])
    steps.append(engine+['run']+common+[
        '--library',str(args.library.resolve()),'--repeatmasker',str(repeatmasker_path(args).absolute()),
        '--jobs',str(args.jobs),'--chunk-bp',str(args.chunk_bp),'--halo',str(args.halo),
        '--scratch',str(args.scratch.resolve())])
    return steps


def main():
    p=parser()
    args=p.parse_args()
    if args.halo < 0:
        p.error('--halo must be nonnegative')
    args.out=(args.out or Path('/mnt/d/analysis/allele3/herv')/args.gen/args.sample).resolve()
    steps=commands(args)
    cpus=sorted(os.sched_getaffinity(0))[:args.cores]
    if args.dry_run:
        print(f'[HERV] sample={args.sample}; input={args.effective_input}; chr1-22; '
              f'reference/hap1/hap2; jobs={args.jobs}; CPUs={cpus}')
        print(f'[HERV] output={args.out}; completed chunks are verified and reused')
        for command in steps:
            print(shlex.join(command))
        return
    legacy=Path('/mnt/d/analysis/metag/herv')/args.gen/args.sample
    if args.out == Path('/mnt/d/analysis/allele3/herv')/args.gen/args.sample and not (args.out/'prepared.json').exists() and (legacy/'prepared.json').exists():
        p.error(f'Existing results are still at {legacy}; move the entire sample directory to {args.out}, '
                f'or use --out {legacy} to resume in place. See 0.README.md.')
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/'job.lock').open('a') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            print(f'[herv] SKIP herv: another job is already running in {args.out}',flush=True)
            return
        # Refuse an already active direct engine run before overwriting job status.
        with (args.out/'run.lock').open('a') as core_lock:
            try:
                fcntl.flock(core_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:
                print(f'[herv] SKIP herv: annotation is already running in {args.out}',flush=True)
                return
        execute(args,steps,cpus)


def execute(args,steps,cpus):
    os.sched_setaffinity(0,cpus)
    os.nice(5)
    state=dict(status='running',started=now(),pid=os.getpid(),command=steps[-1],
               commands=steps,cpu_affinity=cpus,jobs=args.jobs,input_source=args.effective_input)
    save(args.out/'job.state.json',state)
    print(f'[herv] START herv sample={args.sample}; input={args.effective_input}; chr1-22; reference/hap1/hap2; '
          f'jobs={args.jobs}; output={args.out}',flush=True)
    try:
        for command in steps:
            subprocess.run(command,check=True)
        result=json.loads((args.out/'genome_fraction.json').read_text())
        with (args.out/'variant_summary.tsv').open() as f:variants=list(csv.DictReader(f,delimiter='\t'))
        lines=[result['settings']['sample']+' / GRCh38 autosomal RepeatMasker census',
               'Sequence-based annotation; variant-conditioned haplotypes, not complete independent assemblies.',
               'Unresolved spans were masked with N. Unknown inserted sequence lengths remain unknown.',
               '', 'label\tERV_with_MaLR_percent_model_length\tERV_with_MaLR_percent_non_N']
        for row in result['categories']:
            if row['category']=='ERV_with_MaLR':
                lines.append(f"{row['label']}\t{row['percent_genome']:.6f}\t{row['percent_non_N']:.6f}")
        lines.extend(['',f"Calendar hours from first annotation start (includes pauses): {result['elapsed_wall_seconds']/3600:.3f}",
            f"Measured CPU core-hours: {sum(result['measured_cpu_seconds'].values())/3600:.3f}",
            'Uncertain reference bases across both haplotypes: '+str(sum(int(r['uncertain_reference_bp']) for r in variants)),
            'Unhandled/audited allele records across both haplotypes: '+str(sum(int(r['unhandled_records']) for r in variants)),
            'Boundary disagreement (bp; ERV_with_MaLR): '+json.dumps(result['boundary_disagreement_bp']),
            '', 'Linear projection for 3202 individuals; one shared reference control:',
            json.dumps(result['cohort_projection'],indent=2),
            '', 'All results: genome_fraction.tsv/json; chromosome_fraction.tsv; variant_summary.tsv;',
            'boundary_audit.tsv; variants/<chromosome>/unhandled.tsv.gz; chunks/*/*/*/repeatmasker.out.gz.',
            'No complete individual FASTA was saved. Temporary sequence chunks are removed after each job.'])
        (args.out/'report.txt').write_text('\n'.join(lines)+'\n')
        state.update(status='complete',completed=now(),report=str(args.out/'report.txt'))
        save(args.out/'job.state.json',state)
        print(f'[herv] DONE herv sample={args.sample}; results={args.out}',flush=True)
    except BaseException as error:
        state.update(status='failed',failed=now(),error=repr(error))
        save(args.out/'job.state.json',state)
        print(f'[herv] FAIL herv: {error}',file=sys.stderr,flush=True)
        raise


if __name__=='__main__':
    main()
