#!/usr/bin/env python3
"""Schedule downstream multi-ALT analyses and inventory the available tests."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config/workflow.json"


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                   help="Ordered scientific steps; default config/workflow.json")
    p.add_argument("--out", type=Path,
                   help="Analysis root containing all module outputs")
    p.add_argument("--python", default=os.environ.get("ALLELE_PYTHON", sys.executable))
    p.add_argument("--plan-only", action="store_true")
    return p


def load_steps(config, out, python, cores, memory_gib, *, pca=None, scan=None):
    """Expand named paths without a shell, imports, or genomic input reads."""
    from workflow_paths import configured_path
    variables = dict(os.environ)
    variables.update(
        ALLELE_RESULT_ROOT=str(__import__('product_paths').defaults()['results']),
        ALLELE_ROOT=str(ROOT), ALLELE_ANALYSIS_ROOT=str(Path(out).absolute()),
        DIR_DATA=os.environ.get("DIR_DATA", configured_path("DIR_BIG")),
        DIR_GEN=str(Path(configured_path("DIR_GEN")).expanduser().absolute()),
        DIR_FILES=os.environ.get("DIR_FILES", configured_path("DIR_FILES")),
        PCA3D_OUT=str(pca if pca is not None else pca_root()),
        ALLELE_SCAN_STORE=str(scan if scan is not None else scan_root(out)),
        ALLELE_PYTHON=str(python), CORES=str(cores),
        SCAN_JOBS=str(min(4, cores)), PCA_THREADS=str(min(4, cores)),
        PCA_MEMORY_MB=str(max(512, min(65536, int(memory_gib * 768)))),
        CENSUS_JOBS=str(min(8, cores)), CENSUS_CORES=str(cores))
    raw = json.loads(Path(config).read_text())
    if set(raw) != {"steps"} or not isinstance(raw["steps"], list) or not raw["steps"]:
        raise ValueError("Workflow config must contain a nonempty steps array")
    token = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")
    def expand(value):
        if not isinstance(value, str) or any(c in value for c in ("\n", "\r", "\0")):
            raise ValueError("Every workflow argument must be a single text value")
        def replace(match):
            key = match.group(1)
            if key not in variables:
                raise ValueError("Undefined workflow variable: " + key)
            return variables[key]
        return token.sub(replace, value)
    result = []
    for item in raw["steps"]:
        if set(item) != {"name", "argv"} or not re.fullmatch(r"[a-z][a-z0-9_]*", item["name"]):
            raise ValueError("Each step needs a unique lower-case name and argv")
        if not isinstance(item["argv"], list) or not item["argv"]:
            raise ValueError("Empty step argv: " + item["name"])
        args = [expand(value) for value in item["argv"]]
        if any(a in {"--submit", "--write-cmd", "--wait", "--job-dir", "--replace"} for a in args):
            raise ValueError("Set execution and replacement options on allele3.sh run, not individual config steps")
        result.append(dict(name=item["name"], argv=args))
    if len({s["name"] for s in result}) != len(result):
        raise ValueError("Duplicate workflow step name")
    return result


def pca_root():
    from workflow_paths import configured_path
    return Path(os.environ.get('PCA3D_OUT', configured_path('DIR_GEN','ancient/38/pca3d'))).expanduser().absolute()


def scan_root(out):
    from workflow_paths import configured_path
    return Path(os.environ.get('ALLELE_SCAN_STORE', str(__import__('product_paths').defaults()['detailed']))).expanduser().absolute()


def check_prerequisites(out, *, pca=None, scan=None):
    """Check only required input paths, without size/hash or receipt audits."""
    from analysis_integrity import resolve_result_input
    out = Path(out).expanduser().absolute()
    pca = Path(pca) if pca is not None else pca_root()
    inputs = {
        "1.ancient_pca.sh": [pca / name for name in (
            "tables/sample_metadata.tsv", "coordinates/all.sscore",
            "basis/ref.acount", "basis/ref.eigenvec.allele", "basis/ld.prune.in")],
        "3.allele3_scan.sh": [Path(resolve_result_input(scan if scan is not None else scan_root(out))) / "run.json"],
    }
    scan_store = Path(scan if scan is not None else scan_root(out))
    if not (scan_store / 'run.json').is_file() and (scan_store / 'metadata/run.json').is_file():
        print('[allele3 input incomplete] Detailed scan contains published reports and metadata only; native runtime inputs must be regenerated with 3.allele3_scan.sh scan --target gnomad --mode detailed.', file=sys.stderr)
    # Census is independent: downstream analyses do not consume its reports.
    failures = []
    for script, paths in inputs.items():
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            failures.append(script)
            print(f"[allele3 input missing] {script}: " + "; ".join(missing),
                  file=sys.stderr, flush=True)
            command = (script + " scan --target gnomad --mode detailed"
                       if script == "3.allele3_scan.sh" else script)
            print("Run the independent " + command + " first, then rerun allele3.sh run.",
                  file=sys.stderr, flush=True)
    if failures:
        raise ValueError("Required inputs missing; no downstream analysis started.")
    print("[allele3 inputs OK] Required PCA files and scan manifest exist.",
          file=sys.stderr, flush=True)


def build_plan(argv, mode, python):
    from batch_workflows import base_plan, add_job, resources, env_paths
    a = parser().parse_args(argv)
    out = (a.out or Path(env_paths()[2])).expanduser().absolute()
    from study_workflow import inputs, designs, show, manifest
    attachments=[]
    automatic = a.config.resolve() == DEFAULT_CONFIG.resolve()
    sources = inputs(out) if automatic else {}
    if automatic:
        show(sources)
        if not a.plan_only and sources['repeat_mask']['status']!='ready':
            raise ValueError('Prepare the segmental-duplication mask first: ./0.mask.sh run')
    pca = Path(sources['pca']['path']) if automatic else pca_root()
    scan = Path(sources['detailed']['path']) if automatic else scan_root(out)
    if not automatic:check_prerequisites(out, pca=pca, scan=scan)
    resource = resources(40 if mode == "hpc" else min(16, os.cpu_count() or 1),
                         160 if mode == "hpc" else 32, 120)
    steps = load_steps(a.config, out, python, resource["cores"], resource["memory_gib"], pca=pca, scan=scan)
    plan = base_plan("all", mode)
    plan["max_running"] = 1
    plan["local_memory_gib"] = resource["memory_gib"]
    common_env=dict(ALLELE_ANALYSIS_ROOT=str(out),PCA3D_ANALYSIS_ROOT=str(out),PCA3D_OUT=str(pca),
                    ALLELE_SCAN_STORE=str(scan),ALLELE_METADATA=str(pca/'tables/sample_metadata.tsv'),
                    ALLELE_COORDINATES=str(pca/'coordinates/all.sscore'),PCA3D_REUSE_BASIS=str(pca/'basis'),
                    PCA3D_METADATA=str(pca/'tables/sample_metadata.tsv'),PCA3D_EIGENVEC=str(pca/'coordinates/all.sscore'))
    if not automatic:
        done=out/'qc/all.json'
        add_job(plan,'configured_analyses',[python,__file__,'worker','--out',out,'--steps-json',json.dumps(steps)],
                python,resource=resource,outputs=[done])
        plan['jobs'][-1]['env'].update(common_env)
    else:
        products = __import__('product_paths').defaults()['results']
        studies = designs(sources, products)
        requirements={'environment':(), 'hotspots':('detailed','atlas'), 'density':('detailed',),
                      'herv_enrichment':('detailed','pca','herv','ancient_genotypes','ancient_liftover'),
                      'allele_time':('pca','ancient_genotypes','ancient_liftover'), 'space':('pca','ancient_genotypes'),
                      'base_simulation':()}
        parents={'allele_time':['hotspots'],'space':['allele_time']}
        core=[]
        ready=set()
        plan['max_running'] = int(os.environ.get('ALLELE_MAX_RUNNING', 0))
        for step in steps:
            name=step['name'];deps=parents.get(name,[])
            missing=[k+': '+sources[k]['reason'] for k in requirements[name] if sources[k]['status']!='ready']
            missing += [dep+' unavailable' for dep in deps if dep not in ready]
            receipt=out/'qc'/'steps'/(name+'.json')
            row=dict(module=name,reference='',target='',population='',design='core',out=str(receipt),
                     status='unavailable' if missing else 'ready',reason=' | '.join(missing),masks='')
            core.append(row)
            if missing:continue
            ready.add(name)
            add_job(plan,name,[python,__file__,'core-step',receipt,name,json.dumps(step['argv'])],python,
                    deps=(deps or (['environment'] if name!='environment' else [])),resource=resource,outputs=[receipt])
            plan['jobs'][-1]['env'].update(common_env)
        cohorts=sorted({row['gu_cohort'] for row in studies if row['gu_cohort'] and row['status']=='ready'})
        for cohort in cohorts:
            source=sources['archaic/'+cohort]
            dest=Path(source['track']).parent
            add_job(plan,'gu_track.'+cohort,[python,ROOT/'f/gu_segments.py','--raw',source['path'],
                    '--out',dest,'--chain',source['chain'],'--liftover',source['lift']],python,
                    resource=resources(1,8,24),outputs=[source['track'],dest/'annotations/archaic.bed',
                                                       dest/'annotations/archaic_coverage.tsv'])
        for row in studies:
            if row['status']!='ready':continue
            cohort=row['gu_cohort']
            job_id='.'.join(x for x in (row['module'],cohort,Path(row['out']).name) if x)
            add_job(plan,job_id,['bash',ROOT/'allele3.sh',*row['argv']],python,
                    deps=['gu_track.'+cohort] if cohort else [],resource=resource,outputs=[manifest(row['out'])])
        # Individual haplotype assignments require evidence beyond native unphased raw calls.
        studies.append(dict(module='archaic_carrier_genotypes',reference='archaic',target='1kg',population='ALL',
                       design='individual_followup',out=str(products/'archaic_carrier_genotypes'/'1kg'),status='skipped',
                       reason='GU native unphased tract overlap only',masks=''))
        parents_used={d for job in plan['jobs'] for d in job['deps']}
        dependencies=[job['id'] for job in plan['jobs'] if job['id'] not in parents_used]
        inventory = [{k:v for k,v in row.items() if k!='argv'} for row in [*core,*studies]]
        inventory_data=dict(rows=inventory,sources=sources)
        digest=hashlib.sha256(json.dumps(inventory_data,sort_keys=True).encode()).hexdigest()[:24]
        inventory_file=out/'jobs'/('study_inputs_'+digest+'.json')
        attachments.append((inventory_file,inventory_data))
        add_job(plan,'study_inventory',[python,ROOT/'f/study_workflow.py','publish',products,
                    inventory_file],python,deps=dependencies,
                resource=resources(1,2,1),outputs=[products/'results/study_status.tsv',products/'reports/study_index.html'])
        plan['provenance'].update(studies=inventory,scheduling='explicit dependencies and independent cohorts; inventory after successful scheduled jobs')
    if automatic:
        inactive=[s['module']+'/'+(s['target'] or 'core')+' ('+s['design']+')' for s in inventory if s['status']=='unavailable']
        if inactive:print('[allele3] Unavailable studies: '+', '.join(inactive),file=sys.stderr,flush=True)
    return plan, attachments, out


def run_steps(steps, out, *, runner=None):
    """Stop immediately on failure and record unstarted steps as blocked."""
    from analysis_integrity import atomic_json
    from batch_workflows import pipeline_file_lock, pipeline_worker_guard
    pipeline_worker_guard()
    out = Path(out).absolute()
    check_prerequisites(out)
    # Check coordinate provenance before hours of downstream computation.
    # Full receipt/hash checks belong to the allocated worker, not the planner.
    if runner is None:
        from pca import frame_read_scores
        frame_read_scores(pca_root() / 'coordinates/all.sscore', 2,
                          basis_dir=pca_root() / 'basis')
    runner = runner or subprocess.run
    record = dict(status="running", steps=[], started=time.time())
    path = out / "qc/all.json"
    environment = dict(os.environ, WORKFLOW_WORKER="1", ALLELE_BATCH_WORKER="1",
                       ALLELE_IN_COMPUTE="1", ALLELE_ANALYSIS_ROOT=str(out),
                       CUDA_VISIBLE_DEVICES="", WORKFLOW_DEVICE="cpu")
    with pipeline_file_lock(out / "qc/all.lock", nonblocking=True):
        atomic_json(path, record)
        for index, step in enumerate(steps):
            print(f"[allele3 {index + 1}/{len(steps)}] {step['name']}: {shlex.join(step['argv'])}", flush=True)
            row = dict(name=step["name"], status="running", started=time.time())
            record["steps"].append(row)
            atomic_json(path, record)
            try:
                code = runner(step["argv"], cwd=ROOT, env=environment).returncode
            except OSError as exc:
                code = 1
                row["error"] = str(exc)
            row.update(status="complete" if code == 0 else "failed", returncode=code, finished=time.time())
            if code:
                record["steps"].extend(dict(name=s["name"], status="blocked", blocked_by=step["name"])
                                       for s in steps[index + 1:])
                record.update(status="failed", finished=time.time())
                atomic_json(path, record)
                return code
            atomic_json(path, record)
        record.update(status="complete", finished=time.time())
        atomic_json(path, record)
        return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1]==['core-step']:
        from analysis_integrity import atomic_json
        from batch_workflows import pipeline_worker_guard
        pipeline_worker_guard()
        dest,name,args=Path(argv[1]),argv[2],json.loads(argv[3])
        atomic_json(dest,dict(module=name,status='running'))
        try:
            if name in ('herv_enrichment','allele_time','space'):
                from pca import frame_read_scores
                frame_read_scores(pca_root()/'coordinates/all.sscore',2,basis_dir=pca_root()/'basis')
            code=subprocess.run(args,cwd=ROOT).returncode
        except Exception as exc:
            atomic_json(dest,dict(module=name,status='failed',returncode=1,error=str(exc)))
            raise
        atomic_json(dest,dict(module=name,status='complete' if code==0 else 'failed',returncode=code))
        return code
    if argv[:1] == ["worker"]:
        p = argparse.ArgumentParser(description=__doc__)
        p.add_argument("--out", required=True)
        p.add_argument("--steps-json", required=True)
        a = p.parse_args(argv[1:])
        return run_steps(json.loads(a.steps_json), a.out)
    parser().parse_args(argv)
    raise ValueError("Use bash allele3.sh run to create or submit the workflow")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        raise SystemExit("allele3-all: " + str(exc))
