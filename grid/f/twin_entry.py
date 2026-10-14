#!/usr/bin/env python3
"""Thin public dispatcher; preserved upstream benchmark scripts remain unchanged."""
import os,sys,subprocess,shlex,shutil
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
HELP='''GRID — genetic-profile twin confidence

  ./grid.sh abm --check
  ./grid.sh abm --traits height,ldl,t2dm
  ./grid.sh final --analysis-root /mnt/d/analysis/grid
  ./grid.sh shiny --analysis-root /mnt/d/analysis/grid
  ./grid.sh demo --output /tmp/grid-smoke
  ./grid.sh project --model MODEL --data test_X.npz --output local_prediction.joblib
  ./grid.sh evaluate --help

ABM consumes the EXISTING prepared.joblib and preserves its outer test roster.
ABM defaults: 4 threads; config/twin.json beside grid.sh; /tmp/grid-cache/grid/prepared.joblib.
--data-file plus --schema can provide an explicit equivalent numeric dataset.
New results: <out-root>/twin/<trait>/<run>. Old results are not overwritten.
The direct fit/project/evaluate entry points keep test labels out of projection.

Upstream PCA/CSx/Disco scripts remain byte-identical. PRSformer is delegated.
  ./grid.sh prsformer [original options]
  ./grid.sh legacy-abm [original options]
  ./grid.sh legacy-final [original options]

final indexes saved results only; no model is re-fitted by opening Shiny.
Sequences abm,final,shiny accept only --traits and --dry-run. Use separate commands
for custom input/output paths or benchmark-specific flags. Shiny is foreground.
'''


def r_runtime(explicit=None):
    """Find an R with the viewer dependencies, keeping system/Conda libs separate."""
    override=explicit or os.getenv('GRID_RSCRIPT') or os.getenv('R_BIN')
    candidates=[override] if override else [shutil.which('Rscript'),'/usr/bin/Rscript']
    probe="p<-c('shiny','DT','DBI','RSQLite','ggplot2','jsonlite');quit(status=if(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))0L else 1L)"
    for candidate in dict.fromkeys(x for x in candidates if x):
        command=[candidate]
        resolved=shutil.which(candidate)
        if resolved and Path(resolved).resolve()==Path('/usr/bin/Rscript').resolve():
            command=['env','-u','R_ENVIRON_USER','-u','R_LIBS_USER',candidate]
        try:
            if subprocess.run([*command,'--vanilla','-e',probe],capture_output=True).returncode==0:return command
        except OSError:
            pass
    raise RuntimeError('No usable Shiny R runtime. Install R_PACKAGES.txt dependencies or set GRID_RSCRIPT.')


def main(argv=None):
    args=list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ('-h','--help','help'):print(HELP);return
    command=args.pop(0)
    if ',' in command:
        commands=command.split(',')
        if len(set(commands))!=len(commands) or not set(commands)<={'abm','final','shiny','prsformer'}:raise ValueError('Invalid/repeated sequence module')
        if 'shiny' in commands and commands[-1]!='shiny':raise ValueError('Shiny must be the last module')
        shared=[];dry=False;i=0
        while i<len(args):
            if args[i]=='--traits' and i+1<len(args):shared+=args[i:i+2];i+=2
            elif args[i]=='--dry-run':dry=True;i+=1
            else:raise ValueError('Only --traits and --dry-run are shared sequence options')
        for cmd in commands:main([cmd,*shared,*(['--dry-run'] if dry else [])])
        return
    dry='--dry-run' in args;args=[x for x in args if x!='--dry-run'];env=None
    if command in ('prsformer','legacy-abm','legacy-final'):
        upstream=ROOT/'grid.upstream.sh'
        if not upstream.is_file():raise FileNotFoundError('Run install.py into the original scripts/grid first')
        cmd=['bash',str(upstream),command.replace('legacy-',''),*args]
        # The preserved entry's logger hardcodes grid.sh. Invoke the logger here
        # with the preserved entry so legacy-abm cannot jump into the new ABM.
        logger=ROOT.parent/'0f/console_run.py'
        quiet=any(x in ('--help','-h','help','--python-help') for x in args)
        if not dry and not quiet and os.getenv('SCRIPT_CONSOLE_ACTIVE')!='1' and os.getenv('SCRIPT_VERBOSE')!='1' and logger.is_file():
            cmd=[sys.executable,str(logger),'--script',str(upstream),'--',*cmd[2:]]
            env=dict(os.environ,SCRIPT_LOG_DIR=os.getenv('SCRIPT_LOG_DIR','/tmp/grid-logs/runs'))
        else:
            env=dict(os.environ,SCRIPT_CONSOLE_ACTIVE='1')
    elif command in ('abm','fit','project','evaluate','demo'):
        cmd=[sys.executable,str(ROOT/'f/twin_cli.py'),'run' if command=='abm' else command,*args]
    elif command=='final':cmd=[sys.executable,str(ROOT/'f/browser_index.py'),*args]
    elif command=='shiny':
        import argparse
        p=argparse.ArgumentParser();p.add_argument('--analysis-root',default='/mnt/d/analysis/grid');p.add_argument('--db');p.add_argument('--traits',default='height,ldl,t2dm');p.add_argument('--port',type=int,default=3839);p.add_argument('--host',default='127.0.0.1');p.add_argument('--r-bin');p.add_argument('--review',action='store_true');p.add_argument('--reindex',action='store_true')
        a=p.parse_args(args);db=Path(a.db or str(Path(a.analysis_root)/'evidence/evidence.sqlite'))
        if not 1024<=a.port<=65535:raise ValueError('Invalid port')
        if a.host not in ('127.0.0.1','::1','localhost') and os.getenv('GRID_SHINY_ALLOW_NETWORK')!='TRUE':raise ValueError('Network binding requires GRID_SHINY_ALLOW_NETWORK=TRUE and your own authentication')
        if a.reindex or not db.is_file():
            ix=[sys.executable,str(ROOT/'f/browser_index.py'),'--analysis-root',a.analysis_root,'--out',str(db.parent),'--traits',a.traits]
            print(shlex.join(ix),flush=True)
            if not dry:subprocess.run(ix,check=True)
        runtime=[a.r_bin or os.getenv('GRID_RSCRIPT') or os.getenv('R_BIN','Rscript')] if dry else r_runtime(a.r_bin)
        cmd=[*runtime,str(ROOT/'shiny/app.R'),'--db',str(db),'--port',str(a.port),'--host',a.host,*(['--review'] if a.review else [])]
    else:raise ValueError('Unknown command: '+command)
    print('[GRID] '+shlex.join(cmd),flush=True)
    if not dry:subprocess.run(cmd,check=True,env=env)

if __name__=='__main__':
    try:main()
    except subprocess.CalledProcessError as exc:
        print('ERROR: '+str(exc),file=sys.stderr);sys.exit(exc.returncode if exc.returncode>0 else 128-exc.returncode)
    except (ValueError,OSError,RuntimeError) as exc:print('ERROR: '+str(exc),file=sys.stderr);sys.exit(2)
