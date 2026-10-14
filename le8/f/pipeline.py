#!/usr/bin/env python3
"""Public C1/final bridge. Native statistics and ABM publish independent transactions."""
from pathlib import Path
import argparse
import os
import re
import runpy
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / 'f'


def parse(argv):
    options = argparse.ArgumentParser(add_help=False)
    options.add_argument('--abm-config', default=str(ROOT/'config/abm.json'))
    options.add_argument('--deviation-config')
    options.add_argument('--objective-file')
    options.add_argument('--pgs-root')
    options.add_argument('--pgs-mapping')
    options.add_argument('--remerge-pgs', action='store_true')
    options.add_argument('--abm-host-args', default='')
    research, rest = options.parse_known_args(argv)
    common = runpy.run_path(str(ROOT/'f/0.common.py'), run_name='le8_dispatch_config')
    args, unknown = common['parser']().parse_known_args(rest)
    modules = common['canonical_modules'](args.modules)
    if set(modules)-set(common['ORDER']):
        raise ValueError('Unknown module: '+','.join(sorted(set(modules)-set(common['ORDER']))))
    if args.run_abm:
        modules.append('c1_abm')
    if args.run_abm is False:
        modules = [m for m in modules if m != 'c1_abm']
    modules = [m for m in common['ORDER'] if m in modules]
    if not args.Y or any(not re.fullmatch(r'[A-Za-z0-9_]+',y) for y in args.Y.split(',')):
        raise ValueError('Invalid outcome names')
    if not set(args.biom.split(',')) <= {'prot','met'}:
        raise ValueError('Invalid omic layer')
    if args.cores is not None and args.cores < 1:
        raise ValueError('--cores must be positive')
    if args.host not in {'127.0.0.1','localhost','::1'} and os.getenv('LE8_SHINY_ALLOW_NETWORK') != 'TRUE':
        raise ValueError('Network binding requires LE8_SHINY_ALLOW_NETWORK=TRUE')
    if not 1024 <= args.port <= 65535 or args.max_table_mb < 1:
        raise ValueError('Invalid port or table limit')
    # The shell supplies a positional module list. Direct Python use also works.
    tail = rest[1:] if rest and not rest[0].startswith('-') else rest
    tail = [v for v in tail if v not in ('--run-abm','--skip-abm')]
    return research, args, modules, tail, common, unknown


def worker_command(research, a, trait, layer):
    cmd = [os.getenv('ABM_PYTHON') or sys.executable, str(HERE/'run_host.py'),
           '--Y',trait,'--biom',layer,'--analysis-root',str(a.analysis_root),
           '--config',research.abm_config,'--seed',str(a.seed),'--cores',str(a.cores or 4),
           '--r-bin',a.r_bin,'--replace',a.replace,
           '--host-args='+research.abm_host_args]
    for key in ['ukb_phe','group_file','group_col','outer_roster','end_date','diagnosis_col','shared_covariates']:
        value = getattr(a,key)
        if value is not None:
            cmd += ['--'+key.replace('_','-'),str(value)]
    for key in ['deviation_config','objective_file','pgs_root','pgs_mapping']:
        value = getattr(research,key)
        if value:
            cmd += ['--'+key.replace('_','-'),value.replace('{biom}',layer)]
    if research.remerge_pgs:
        cmd += ['--remerge-pgs']
    return cmd


def main(argv=None):
    research,a,modules,tail,common,unknown = parse(list(sys.argv[1:] if argv is None else argv))
    def call(cmd, dry=None):
        cmd=list(map(str,cmd));print('[LE8 research] '+shlex.join(cmd),flush=True)
        if not (a.dry_run if dry is None else dry):
            subprocess.run(cmd,check=True)
    if any(t.startswith(('--legacy','--abm-backend','--abm-engine','--abm-args','--s7-','--abm-tf-root','--abm-root')) for t in unknown):
        raise ValueError('Obsolete ABM options are unsupported; use --abm-config and --abm-host-args')
    if unknown and not (set(modules)&common['NATIVE'] or a.fit_joint or a.fit_reference or a.fit_genetic or a.details):
        raise ValueError('Unrecognized options: '+' '.join(unknown))
    # Fail before a native C1 scan if research configuration or files are unavailable.
    if 'c1_abm' in modules:
        for trait in a.Y.split(','):
            for layer in a.biom.split(','):
                call([*worker_command(research,a,trait,layer),'--preflight'])
    for module in modules:
        if module == 'share':
            if modules != ['share']:
                raise ValueError('Run share separately after preparing final/shiny')
            if a.preflight:
                call([a.r_bin,'-e',"p<-c('shiny','bslib','DT','DBI','RSQLite','ggplot2','plotly','openxlsx','jsonlite');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste(m,collapse=', '))"])
            else:
                call([sys.executable,HERE/'share.py','--analysis-root',a.analysis_root,
                      '--out',a.share_out or a.analysis_root.parent/'le8-share.zip',
                      '--r-bin',a.r_bin,'--max-table-mb',a.max_table_mb])
        elif module == 'c1_abm':
            if a.preflight:
                continue
            for trait in a.Y.split(','):
                for layer in a.biom.split(','):
                    worker=worker_command(research,a,trait,layer)
                    budget=a.memory_limit_gb if a.memory_limit_gb is not None else 32
                    if budget < 0:
                        raise ValueError('Invalid memory budget')
                    # Same cgroup guard and CPU lease as the existing dispatcher.
                    if budget:
                        if os.getenv('LE8_MEMORY_SCOPE_ACTIVE')!='1':
                            worker=['systemd-run','--user','--scope','--quiet','--collect',
                                    '-p',f'MemoryMax={budget}G','-p',f'MemorySwapMax={a.memory_swap_gb if a.memory_swap_gb is not None else 4}G',
                                    '-p','OOMPolicy=kill','--','env','LE8_MEMORY_SCOPE_ACTIVE=1',*worker]
                        worker=[sys.executable,ROOT/'f/0.resources.py','phase','--label','c1_abm',
                                '--cores',str(a.cores or 4),'--memory-gib',str(budget),'--',*worker]
                    call(worker)
        elif module in {'final','shiny'}:
            if module == 'final':
                call([sys.executable,ROOT/'f/0.common.py','dispatch','final',*tail])
            if not a.preflight and (module=='final' or not a.no_reindex):
                call([sys.executable,HERE/'evidence_index.py','--analysis-root',a.analysis_root,
                      '--Y',a.Y,'--biom',a.biom,'--max-table-mb',a.max_table_mb])
            if module == 'shiny':
                if a.preflight:
                    call([a.r_bin,'-e',"p<-c('shiny','bslib','DT','DBI','RSQLite','ggplot2','plotly');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste(m,collapse=', '))"])
                elif not a.prepare_only:
                    call([a.r_bin,ROOT/'shiny/app.R','--db',a.analysis_root/'evidence/evidence.sqlite',
                          '--port',a.port,'--host',a.host,*(['--review'] if a.shiny_review else [])])
        else:
            call([sys.executable,ROOT/'f/0.common.py','dispatch',module,*tail])


if __name__ == '__main__':
    try:
        main()
    except (ValueError,OSError,subprocess.CalledProcessError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr)
        sys.exit(2)
