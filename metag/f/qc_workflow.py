"""Common explicit LSF submission for downstream modules; no shell interpolation."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from workflow_entry import ScriptWorkflow
import sys

from common import CODE_ROOT, atomic_json, load_config, resolve_mode, sha256, default_cores

MODULES = {
    "library_check": (2, 16, 120),
    "hierarchical_split": (24, 192, 120),
    "group_library": (32, 192, 120),
    "validation_suite": (32, 256, 120),
    "threshold_diagnose": (24, 128, 48),
    "plant_validate": (32, 256, 48),
    "group_summary": (2, 16, 4),
    "human_validate": (32, 128, 24),
    "sample_diagnose": (2, 16, 24),
    "species_refine": (32, 256, 48),
    "sylph_profile": (32, 128, 24),
    "split_summary": (2, 16, 24),
}


def active_job(job):
    if not re.fullmatch(r'[1-9][0-9]*',str(job)):return False
    run=subprocess.run(['bjobs','-noheader','-o','stat',str(job)],text=True,capture_output=True,timeout=20)
    return run.returncode==0 and bool(set(run.stdout.split()) & {'RUN','PEND','PSUSP','USUSP','SSUSP','WAIT'})


def hierarchy_scope(module,worker_args,cfg,known):
    from group_library import group_directory
    group_dir=str(Path(known.out).expanduser().resolve()) if module=='group_library' and known.out else str(group_directory(cfg))
    return dict(module=module,action=worker_args[0] if module=='hierarchical_split' and worker_args else '',
                samples=known.sample or [],datasets=(known.dat_list or '').split(',') if known.dat_list else [],
                output_dir=cfg['output_dir'],library_dir=cfg['library_dir'],
                group_library_dir=group_dir,index_use=cfg.get('index_use'))


PATH_OPTIONS = {'--out','--output-dir','--representatives','--manifest','--input-root',
                '--compare-root','--truth','--reference','--bwa-mem2','--assignments','--summary',
                '--counts','--r1','--r2','--minimap2','--sylph-candidates','--db','--sylph'}


def frozen_args(argv, snapshot, cwd):
    result=[];i=0
    while i<len(argv):
        arg=argv[i];key,sep,value=arg.partition('=')
        if key in ('--config','--root','--library','--index-use'):
            i+=1 if sep else 2;continue
        if key in PATH_OPTIONS:
            if not sep:i+=1;value=argv[i]
            if key in ('--bwa-mem2','--minimap2','--sylph') and '/' not in value:
                value=shutil.which(value) or value
            result.extend([key,str((cwd/Path(value).expanduser()).resolve())])
        elif key=='--target-fastq':
            result.append(key)
            if sep:result.append(str((cwd/Path(value).expanduser()).resolve()))
            while i+1<len(argv) and not argv[i+1].startswith('--'):
                i+=1;result.append(str((cwd/Path(argv[i]).expanduser()).resolve()))
        else:result.append(arg)
        i+=1
    # Place global options before sylph's subcommand as well.
    return ['--config',str(snapshot),*result]


def build_workflow(argv=None):
    module, *argv = list(sys.argv[1:] if argv is None else argv)
    if module not in MODULES:
        raise ValueError(f"Unknown QC module {module}")
    p = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    p.add_argument("--submit", choices=("hpc", "local", "auto"), default=None)
    p.add_argument("--plan-only", action="store_true")
    if any(x.split("=", 1)[0] in ("--wait-for", "--wait-mode") for x in argv):
        p.error("External job dependencies are not supported; complete the required input files before submitting")
    p.add_argument("--queue", default="38")
    p.add_argument("--cores", type=int, default=None)
    p.add_argument("--mem-gib", type=int, default=None if module == "hierarchical_split" else MODULES[module][1])
    p.add_argument("--job-name", default=None)
    p.add_argument("--hours", type=int, default=MODULES[module][2])
    opts, worker_args = p.parse_known_args(argv)
    if opts.job_name and not re.fullmatch(r'[A-Za-z0-9_.-]+', opts.job_name):
        p.error('Invalid job name')
    worker = [sys.executable, str(CODE_ROOT / "f" / (module + ".py")), *worker_args]
    mode = resolve_mode(opts.submit)
    if opts.cores is None:
        opts.cores = default_cores(mode) if MODULES[module][0] > 2 else MODULES[module][0]
    os.environ["submit"] = os.environ["SUBMIT"] = mode
    if min(opts.cores, opts.hours) < 1 or (opts.mem_gib is not None and opts.mem_gib < 1):
        p.error("Resource limits must be positive")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", opts.queue):
        p.error("Invalid queue name")
    # Independent submissions require their external inputs to be ready now.
    config_parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    config_parser.add_argument("--config"); config_parser.add_argument("--root")
    config_parser.add_argument('--library')
    config_parser.add_argument('--index-use')
    config_parser.add_argument("--threads", type=int)
    config_parser.add_argument('--sample',action='append')
    config_parser.add_argument('--dat-list')
    config_parser.add_argument('--out')
    known, _ = config_parser.parse_known_args(worker_args)
    if known.threads is not None and not 1 <= known.threads <= opts.cores:
        p.error("--threads must be positive and no larger than the submitted --cores")
    cfg = load_config(known.config, known.root)
    if known.library:
        from library_paths import select_library
        cfg = select_library(cfg, known.library)
    if known.index_use:
        if known.library:p.error('--library and --index-use cannot be combined')
        from index_selection import select_index
        cfg=select_index(cfg,known.index_use)
    if module == 'group_library':
        # The builder's preflight budget must match this job's allocation.
        # Otherwise --mem-gib alone could request a larger job but retain a
        # stale smaller guard (or request less RAM than the guard assumes).
        cfg['group_build_memory_gib'] = opts.mem_gib
    if module == 'hierarchical_split':
        # Admission and frozen worker budgets must match the LSF reservation.
        if opts.mem_gib is None:
            opts.mem_gib = cfg['classify_memory_gib']
        cfg.update(classify_memory_gib=opts.mem_gib, classify_cores=opts.cores,
                   classify_threads=known.threads or max(1, opts.cores - 2), classify_workers=1)
    scope=hierarchy_scope(module,worker_args,cfg,known)
    from submission_inputs import check_submission_inputs
    try:
        check_submission_inputs(module, cfg, worker_args)
    except (OSError, ValueError, RuntimeError) as exc:
        p.error("INPUT NOT READY: " + str(exc) + "; no job submitted")
    if mode == 'hpc' and module == 'group_library':
        # Refuse to replace an index with active readers; never queue behind them.
        for path in (Path(cfg['work_dir'])/'jobs/qc/hierarchical_split').glob('run*/submission.json'):
            rec = json.loads(path.read_text())
            if (rec.get('status') == 'submitted'
                    and rec.get('scope', {}).get('group_library_dir') == scope['group_library_dir']
                    and active_job(rec.get('job_id', ''))):
                p.error('Index is in use by an independent reader; no job submitted: ' + str(path))
    from source_release import freeze_source
    release = freeze_source(cfg)
    worker[1] = str(release / "f" / (module + ".py"))
    # Streaming classifiers reserve two cores for input/output helpers.
    default_threads = {name: opts.cores for name in (
        'group_library', 'human_validate', 'validation_suite',
        'species_refine', 'plant_validate', 'sylph_profile')}
    default_threads.update(hierarchical_split=max(1, opts.cores - 2),
                           threshold_diagnose=max(1, opts.cores - 2))
    if module in default_threads and known.threads is None:
        worker += ["--threads", str(min(opts.cores, default_threads[module]))]
    local = mode == "local"
    if local and module != "split_summary" and "--allow-local" not in worker:
        worker.append("--allow-local")
    from storage_layout import new_run
    folder = new_run(Path(cfg["work_dir"]) / "jobs" / "qc" / module)
    folder=folder.resolve()
    if module == 'hierarchical_split' and scope['action'] != 'report':
        cfg['hierarchy_report_dir'] = str(folder / 'report')
    snapshot=folder/'effective-config.json' 
    atomic_json(snapshot,cfg)
    config_digest=sha256(snapshot)
    snapshot.chmod(0o444)
    worker=[sys.executable,str(release/'f'/(module+'.py')),
            *frozen_args(worker[2:],snapshot,Path.cwd())]
    logs = Path(cfg["work_dir"]) / "logs" / "qc" / module / folder.name
    logs.mkdir(parents=True, exist_ok=True)
    path = folder / "run.cmd"
    body = ["#!/usr/bin/env bash", "set -euo pipefail",
            "source " + shlex.quote(str(release / "f" / "env.sh")),
            "cd " + shlex.quote(str(Path.cwd())),
            shlex.join([sys.executable, str(release.parent/'0f'/'submit.py'), '--project', 'metag-qc',
                        '--verify-config',str(snapshot),config_digest,str(release)]),
            "exec " + shlex.join(worker)]
    receipt = {"worker": worker,
               "scope": scope,
               "dependency": "", "status": "planned", "source_release": str(release),
               "effective_config": str(snapshot), "config_sha256": config_digest, "submission_cwd": str(Path.cwd())}
    return ScriptWorkflow(path, body,
                          dict(queue=opts.queue, cores=opts.cores, hours=opts.hours, memory_gib=opts.mem_gib),
                          opts.job_name or ('metag.sample' if module == 'hierarchical_split' else 'metag_' + module),
                          logs, mode, receipt, opts.plan_only)


def verify_config(argv):
    snapshot, digest, release = argv
    if sha256(snapshot) != digest: raise RuntimeError('Frozen QC configuration checksum mismatch')
    root = Path(release).parent
    receipt = json.loads((root/'release.json').read_text())
    for relative, expected in receipt['files'].items():
        if sha256(root/relative) != expected: raise RuntimeError('Frozen QC source checksum mismatch: ' + relative)
