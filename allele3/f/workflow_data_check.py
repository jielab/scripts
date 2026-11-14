#!/usr/bin/env python3
"""Audit actual workflow inputs, then execute within the existing cgroup/log.

Bounded header/reference checks do not certify a complete dataset or silently
convert genome builds. Module-specific allele harmonization remains mandatory.
"""
from __future__ import annotations

# Workflow Data Check
# ======================================================================
# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path
import argparse
import glob
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def options(argv):
    values = {}
    for i, arg in enumerate(argv):
        if arg.startswith('--'):
            if '=' in arg:
                key, value = arg.split('=', 1); values[key] = value
            elif i+1 < len(argv) and not argv[i+1].startswith('--'):
                values[arg] = argv[i+1]
    return values


def actual_sources(workflow, argv, registry):
    """Bind paths to the current command, not just to registry defaults."""
    opt = options(argv)
    by_id = {e['id']: e for e in registry}
    chosen = []
    def add(name, path=None, build=None, kind=None, suffix=''):
        e = dict(by_id[name])
        if path is not None: e['path'] = str(path)
        if build is not None: e['build'] = str(build).replace('GRCh', '')
        if kind is not None: e['type'] = kind
        e['id'] += suffix
        chosen.append(e)
    fasta37 = opt.get('--fasta37', os.environ.get('ALLELE_FASTA37', configured_path('DIR_GEN', '1hgp/37/GRCH37.fasta.gz')))
    fasta38 = opt.get('--fasta38', os.environ.get('ALLELE_FASTA', configured_path('DIR_GEN', '1hgp/38/GRCH38.fasta.gz')))
    ordinary = {'pipeline','primary','archaic','summary','spatial','windows','alleles','time','selection','candidates','genotypes','ld','report'}
    if workflow in ordinary:
        import herv as herv_pipeline
        a = herv_pipeline.parser().parse_args(argv)
        paths = herv_pipeline.resolve(a); stages = set(herv_pipeline.stage_selection(a))
        if stages & {'primary','archaic','summary','spatial','windows','alleles','candidates'}:
            add('gnomad_scan', paths['modern'])
        if stages & {'spatial','windows','alleles','candidates','ld'}:
            add('herv_annotation', a.erv_joined)
        if 'time' in stages:
            add('ancient_pgen', paths['prefix'])
            add('ancient_dates', paths['metadata'])
            add('ancient_hq', paths['hq'])
            print('[build-check] genomic joins: modern labels GRCh38 <-> ancient GRCh37; use exact nucleotide mapping, not rsID-only joins', flush=True)
        if 'genotypes' in stages:
            add('1kg_vcf', a.kg_vcf_dir, kind='vcf_directory')
        if 'archaic' in stages:
            # Discovery audits source caller contracts. Do not invent a build
            # for newly named UKB cohorts or guess which module produced them.
            from herv import archaic_discover as discover
            for item in discover(a.archaic_root, a.archaic_datasets, a.archaic_module,
                                 a.archaic_chromosomes, a.archaic_references):
                if item.get('status') == 'ready':
                    name = item.get('name', item.get('dataset', 'cohort'))
                    print(f"[dataset] GU {name}: build=GRCh{item.get('build', 'unknown')}; evidence=IBDmix unphased called tracts", flush=True)
        fasta37, fasta38 = a.fasta37, a.fasta38
    elif workflow == 'pca3d':
        command = next((x for x in argv if not x.startswith('-')), '')
        if command in {'plan','project-plan','atlas-plan'}:
            print('[dataset] plan only: genomic input sampling is deferred to execution', flush=True)
        elif command in {'render','atlas','centroid','centroids'}:
            print('[dataset] PCA render/atlas uses frozen coordinates; original genotype identities are checked by the PCA contract', flush=True)
        elif command == 'project':
            pattern = opt.get('--target-glob', os.environ.get('PCA3D_PROJECT_TARGET_GLOB', configured_path('DIR_GEN', 'ukb/37/imp/chr[0-9]*.pgen')))
            files = sorted(glob.glob(pattern))
            if files:
                for i, path in enumerate(files): add('ukb_imputed', path[:-5] if path.endswith('.pgen') else path, suffix=f':{i+1}')
            else: add('ukb_imputed', pattern)
            print('[build-check] UKB projection requires exact shared markers, allele orientation and matched INFO; MFI is imputation quality, not WGS allele counts', flush=True)
        else:
            prefix = opt.get('--pgen-prefix', os.environ.get('PCA3D_PGEN_PREFIX', configured_path('DIR_GEN', 'ancient/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes')))
            add('ancient_pgen', prefix)
    elif workflow == 'pca-trajectory':
        print('[dataset] PCA trajectory uses frozen PC coordinates and dated metadata; age is supplied metadata, not PC3', flush=True)
    elif workflow == 'herv':
        # The census chooses a prepared source or an available PGEN runtime.
        # Audit inside that selection, never block on an unused alternative.
        print('[dataset] HERV census will audit its selected PGEN/VCF source after resolving auto/resume; this is a variant-conditioned sequence model', flush=True)
    elif workflow == 'scan-run':
        source = opt.get('--vcf', opt.get('--vcf-dir', opt.get('--dir-vcf', os.environ.get('ALLELE_VCF_DIR', configured_path('DIR_ANNOT', 'allele3/gnomad/sites')))))
        build = opt.get('--build', 'GRCh38')
        add('gnomad_sites', source, build, 'vcf' if '--vcf' in opt else 'vcf_directory')
        if '--fasta' in opt:
            if build in {'37','GRCh37'}: fasta37 = opt['--fasta']
            else: fasta38 = opt['--fasta']
    elif workflow.startswith('base_occupancy-'):
        module = workflow.removeprefix('base_occupancy-')
        if module == 'space' and '--pgen-prefix' in opt and ('--time-dir' in opt or '--alleles-dir' in opt):
            add('ancient_pgen', opt['--pgen-prefix'], '37')
            print('[build-check] space: measurement-plan nucleotide labels are GRCh38; original encoded ancient PGEN is GRCh37 and is identity-bound', flush=True)
        elif module == 'domains':
            add('gnomad_scan', opt.get('--modern-dir', by_id['gnomad_scan']['path']))
            add('herv_annotation', opt.get('--domain-bed', ''), opt.get('--domain-build', 'unknown'))
            fasta38 = opt.get('--fasta', fasta38)
        elif module == 'archaic-alleles':
            add('gu_1kg', opt.get('--gu-dataset', by_id['gu_1kg']['path']), opt.get('--modern-build', 'unknown'))
            chosen[-1]['module'] = opt.get('--module', 'ibdmix')
            if '--modern-vcf-template' in opt:
                add('1kg_vcf', opt['--modern-vcf-template'], opt.get('--modern-build','unknown'))
                chosen[-1].update(id='archaic_link_modern_vcf', evidence='source_genotypes_to_confirm',
                                 role='Actual modern genotypes for within-person allele/tract linkage')
            print(f"[build-check] target={opt.get('--target-build','38')} -> GU/modern={opt.get('--modern-build','unknown')}; unique mapping and exact VCF REF are required", flush=True)
        elif module in {'data','mfi','allofus'}:
            print(f'[dataset] {module}: the module reports its explicit input schema/build checks', flush=True)
        else:
            print(f'[dataset] {module}: source-specific count/measurement provenance is checked by the module; no implicit build conversion', flush=True)
    return chosen, fasta37, fasta38


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workflow', required=True)
    p.add_argument('--config', default=str(ROOT/'config/data_sources.tsv'))
    p.add_argument('--module-out', help='Public result directory; native artifacts go to sibling tmp/')
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args(argv)
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    if len(command) < 2: raise ValueError('Expected Python executable and workflow entry')
    from analysis_integrity import resolve_cli_inputs
    command = resolve_cli_inputs(command)
    from data_sources import load_config, audit_entries
    entries, fasta37, fasta38 = actual_sources(a.workflow, command[2:], load_config(a.config))
    if entries:
        report = audit_entries(entries, out=None, fasta37=fasta37, fasta38=fasta38,
                               max_records=1000, max_sites=100)
        if report.get('has_conflicts'):
            print('[build-check] conflicting source identities; workload was not started', file=sys.stderr, flush=True)
            return 2
    sys.stdout.flush(); sys.stderr.flush()
    if a.workflow == 'scan-verify':
        # Verification is read-only; never enter the output replacement path.
        os.execvp(command[0], command)
    if a.module_out and os.environ.get('ALLELE_LAYOUT_ACTIVE') != '1':
        return execute_layout(command, a.module_out)
    os.execvp(command[0], command)


# Census Check
# ======================================================================
_description_census_check = 'Check the external annotator before starting expensive genome preparation.\n\nRuns inside the same allocated census worker.\n'

import importlib.util
import re
import subprocess
import tempfile
from types import SimpleNamespace

CENSUS_ROOT = Path(__file__).resolve().parents[1]
CENSUS_ENGINE = CENSUS_ROOT / 'f/census.py'


def census_load_engine():
    sys.path.insert(0, str(CENSUS_ENGINE.parent))
    spec = importlib.util.spec_from_file_location('herv_preflight_engine', CENSUS_ENGINE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def census_short_error(log):
    lines = log.read_text(errors='replace').splitlines()
    errors = [line.strip() for line in lines if re.search(
        r'error|fatal|could not|cannot|not found|no such|permission denied', line, re.I)]
    return (errors[0] if errors else next((line.strip() for line in reversed(lines) if line.strip()), 'unknown failure'))[:300]


def census_preflight(args, engine):
    engine.require_library(args.library)
    # A verified completed run does not launch RepeatMasker again.
    if args.replace == 'FALSE' and (args.out / 'analysis.done.json').is_file():
        return
    log_dir = args.out / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / 'repeatmasker_preflight.log'
    with tempfile.TemporaryDirectory(prefix='herv-preflight-') as directory:
        temp = Path(directory)
        # Fresh mode will archive the prior run.json before resolving the tool.
        selection = SimpleNamespace(repeatmasker=args.repeatmasker,
                                    out=temp if args.replace == 'TRUE' else args.out)
        executable = engine.repeatmasker_path(selection).absolute()
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError(f'RepeatMasker executable not found: {executable}')
        sequence = []
        with args.library.open() as source:
            for line in source:
                if line.startswith('>'):
                    if sequence:
                        break
                    continue
                sequence.extend(line.strip().upper())
                if len(sequence) >= 1000:
                    break
        sequence = ''.join(sequence[:1000])
        if not sequence or set(sequence) - set('ACGTN'):
            raise ValueError(f'Repeat library has no valid first FASTA sequence: {args.library}')
        fasta = temp / 'check.fa'
        fasta.write_text('>preflight\n' + sequence + '\n')
        command = [str(executable), '-e', 'rmblast', '-pa', '1',
                   '-lib', str(args.library.resolve()), '-dir', str(temp), str(fasta)]
        env = dict(os.environ, PATH=str(executable.parent) + os.pathsep + os.environ.get('PATH', ''))
        with log.open('w') as handle:
            result = subprocess.run(command, cwd=temp, env=env, stdin=subprocess.DEVNULL,
                                    stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode or not (temp / 'check.fa.out').is_file():
            raise ValueError(f'RepeatMasker preflight failed: {census_short_error(log)}; log={log}')
    print('[herv] PREFLIGHT RepeatMasker OK', flush=True)


def census_main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    engine = census_load_engine()
    parser = engine.parser()
    parser.add_argument('--replace', choices=('TRUE', 'FALSE'), default='FALSE')
    parser.add_argument('--reuse-from', type=Path)
    args = parser.parse_args(argv)
    if not args.dry_run:
        if args.out is None:
            parser.error('The HERV launcher requires --out')
        args.out = args.out.expanduser().resolve()
        census_preflight(args, engine)
    os.execv(sys.executable, [sys.executable, str(CENSUS_ENGINE), *argv])


# Check Environment
# ======================================================================
_description_check_environment = "Check only the selected stage's actual interpreter and executables."

import importlib
import importlib.util
import json
import shutil


def environment_main(argv=None):
    p=argparse.ArgumentParser(description=_description_check_environment)
    p.add_argument('--stage',choices=['render','trajectory','project','baseline','time','scan','herv','assoc','all'],default='trajectory')
    p.add_argument('--plink2',default='plink2')
    a=p.parse_args(argv)
    packages={'render':['numpy','matplotlib','plotly'],'trajectory':['numpy','scipy','matplotlib'],
              'project':['numpy','pgenlib'],'baseline':['numpy','matplotlib'],
              'time':['numpy','scipy','matplotlib','pgenlib','pysam','pyliftover','zstandard'],
              'scan':['numpy','scipy','matplotlib','pysam'],'herv':['numpy','matplotlib','pysam'],
              'assoc':['numpy','scipy','pandas','pysam','pgenlib','statsmodels','pyreadr','matplotlib']}
    binaries={'project':[a.plink2],'scan':['bcftools'],'herv':['bcftools','samtools']}
    stages=list(packages) if a.stage=='all' else [a.stage]
    result={'python':str(Path(sys.executable).resolve()),'python_version':sys.version,'stage':a.stage,'packages':{},'executables':{},'errors':[]}
    if 'herv' in stages:
        spec=importlib.util.spec_from_file_location('herv_dependency_check',Path(__file__).with_name('census.py'))
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
        executable=module.repeatmasker_path(SimpleNamespace(repeatmasker=None,out=Path('/nonexistent-allele3-check-only')))
        binaries['herv'] += [str(executable),str(executable.parent/'rmblastn')]
        runtime=module.pgen_runtime()
        result['herv_pgen_python']=runtime
        if not runtime:result['errors'].append('No HERV PGEN interpreter with phased allele support')
    for name in sorted({x for s in stages for x in packages[s]}):
        try:
            module=importlib.import_module(name)
            result['packages'][name]={'version':str(getattr(module,'__version__','not_exposed')),'path':getattr(module,'__file__',None)}
        except ImportError as exc:result['errors'].append(str(exc))
    for name in sorted({x for s in stages for x in binaries.get(s,[])}):
        found=shutil.which(name)
        if not found:result['errors'].append('Executable not found: '+name);continue
        path=str(Path(found).absolute())
        try:
            option={'RepeatMasker':'-v','rmblastn':'-version'}.get(Path(name).name,'--version')
            env=dict(os.environ,PATH=str(Path(path).parent)+os.pathsep+os.environ.get('PATH',''))
            version=subprocess.run([path,option],capture_output=True,text=True,errors='replace',timeout=20,env=env)
            result['executables'][name]={'path':path,'version':(version.stdout+version.stderr).splitlines()[:3],'exit_code':version.returncode}
            if version.returncode:result['errors'].append('Version check failed: '+path)
        except (OSError,subprocess.SubprocessError) as exc:result['errors'].append(str(exc))
    result['status']='failed' if result['errors'] else 'verified_selected_dependencies'
    print(json.dumps(result,indent=2));return bool(result['errors'])


def _layout_files(native):
    """Inspect published products and small receipts, never traverse genotype caches."""
    excluded = {'shards', 'variants', 'chunks', 'annotations', 'candidate_cache',
                'runtime', '__pycache__', 'history', 'previous'}
    for parent, directories, names in os.walk(native):
        directories[:] = [n for n in directories if n not in excluded and not n.endswith('.partial')]
        for name in names:
            yield Path(parent) / name


def _layout_destination(layout, relative):
    """Explicit publication policy: tables/figures/reports versus support files."""
    name = relative.name
    parts = relative.parts
    text = str(relative)
    if name.startswith('.') or name.endswith(('.lock', '.partial', '.tmp')):
        return 'qc', relative
    if 'logs' in parts or name.endswith(('.log', '.out', '.err')) and not name.startswith('ld.prune'):
        return 'logs', Path(*[p for p in parts if p != 'logs'])
    if name.endswith(('.sh', '.cmd')) or name == 'commands.sh':
        return 'jobs', relative
    if ('qc' in parts or name.endswith(('.json', '.sha256')) or name in {'PCA_COMPLETE', 'COMPLETE'}
            or 'audit' in name or name.startswith('metadata_qc')):
        return 'qc', Path(*[p for p in parts if p != 'qc'])
    module = layout.key.parts[0]
    if module == 'pca3d':
        if text == 'metadata/sample_metadata.tsv': return 'public', Path('tables/sample_metadata.tsv')
        if text in {'pca/all.sscore', 'pca/ref.eigenvec', 'pca/ref.eigenval'}:
            return 'public', Path('coordinates') / name
        if text in {'pca/ref.eigenvec.allele', 'pca/ref.acount', 'pca/ld.prune.in'}:
            return 'public', Path('basis') / name
        tail = Path(*parts[1:]) if parts[0] == 'pca' else relative
        if name.endswith(('.png', '.pdf', '.svg')): return 'public', Path('figures') / tail
        if name.endswith('.html'): return 'public', Path('reports') / tail
        if name.endswith(('.tsv', '.tsv.gz', '.csv')) and name not in {'overlay_targets.tsv'}:
            return 'public', Path('tables') / tail
        return None
    if module in {'gnomad_scan', 'hotspots', 'density', 'herv_enrichment', 'base_simulation'}:
        return ('public', relative) if parts[0] in {'results', 'reports'} else None
    if module == 'herv_census':
        tables = {'genome_fraction.tsv', 'variant_summary.tsv', 'chromosome_fraction.tsv', 'boundary_audit.tsv'}
        if name in tables: return 'public', Path('results') / name
        if name.endswith(('.html', '.png', '.pdf', '.svg')): return 'public', Path('reports') / relative
        return None
    if parts[0] in {'results', 'reports', 'figures', 'tables'}: return 'public', relative
    if name.endswith(('.html', '.png', '.pdf', '.svg')): return 'public', Path('reports') / relative
    if name.endswith(('.tsv', '.tsv.gz', '.csv', '.bed', '.bed.gz', '.xlsx')):
        return 'public', Path('results') / relative
    return None


def publish_layout(output, *, status='complete'):
    """Publish a clean result view; native scientific receipts remain unchanged."""
    import html
    import json
    import re
    import shutil
    import tempfile
    import time
    from analysis_integrity import module_layout, atomic_json
    layout = module_layout(output)
    for path in (layout.output.parent, layout.qc, layout.logs, layout.jobs, layout.runtime.parent):
        path.mkdir(parents=True, exist_ok=True)
    receipt_path = layout.qc / 'layout.json'
    previous = json.loads(receipt_path.read_text()) if receipt_path.is_file() else {}
    public_files, support_files, destinations = {}, {}, {}
    staging = Path(tempfile.mkdtemp(prefix='publish-', dir=layout.runtime.parent))
    try:
        for source in _layout_files(layout.runtime):
            if source.is_symlink() or not source.is_file():
                continue
            relative = source.relative_to(layout.runtime)
            destination = _layout_destination(layout, relative)
            if destination is None: continue
            kind, target_relative = destination
            if kind == 'public':
                if status != 'complete': continue
                target = staging / target_relative
            else:
                target = getattr(layout, kind) / target_relative
            if target in destinations:
                raise ValueError('Colliding published artifact: ' + str(target))
            destinations[target] = source
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            spec = dict(native=str(relative), size=source.stat().st_size)
            if kind == 'public':
                public_files[str(target_relative)] = spec
            else:
                support_files[str(relative)] = str(target)
        # Generated reports may link to native tables/receipts. Rewrite only the
        # displayed copy; original HTML and its historical checksum stay intact.
        by_native = {str(layout.runtime / spec['native']): str(layout.output / name)
                     for name, spec in public_files.items()}
        by_native.update({str(layout.runtime / name): value for name, value in support_files.items()})
        for name, spec in public_files.items():
            target = staging / name
            if target.suffix == '.html':
                native_parent = (layout.runtime / spec['native']).parent
                def link(match):
                    attribute, quote, url = match.groups()
                    if re.match(r'^(?:[a-zA-Z]+:|#|//)', url): return match.group(0)
                    path, separator, fragment = html.unescape(url).partition('#')
                    absolute = os.path.normpath(str(native_parent / path))
                    if absolute not in by_native: return match.group(0)
                    mapped = os.path.relpath(by_native[absolute], (layout.output / name).parent)
                    if separator: mapped += '#' + fragment
                    return attribute + '=' + quote + html.escape(mapped, quote=True) + quote
                document = re.sub(r'''(href|src)=(["'])(.*?)\2''', link, target.read_text())
                target.write_text(document)
            spec['size'] = target.stat().st_size
        if status == 'complete' and public_files:
            if layout.output.exists() and any(layout.output.iterdir()):
                archive = layout.runtime.parent / 'previous_published'
                if archive.exists(): shutil.rmtree(archive)
                layout.output.rename(archive)
            elif layout.output.exists():
                layout.output.rmdir()
            staging.rename(layout.output)
            archive = layout.runtime.parent / 'previous_published'
            if archive.exists(): shutil.rmtree(archive)
        else:
            layout.output.mkdir(parents=True, exist_ok=True)
            public_files = previous.get('public_files', {})
        aliases = {'.': '.'}
        if layout.key.parts[0] == 'pca3d':
            aliases.update(basis='pca', coordinates='pca', metadata='metadata', pca='pca')
        receipt = dict(schema='allele3-module-layout-1', status=status,
                       public_root=str(layout.output), runtime_root=str(layout.runtime),
                       directory_map=aliases,
                       public_files=public_files, support_files=support_files,
                       provenance='Native scientific receipts retained unchanged; this receipt describes file layout only')
        atomic_json(receipt_path, receipt)
        return receipt
    finally:
        if staging.exists(): shutil.rmtree(staging)


def layout_command(command, output):
    """Keep the public output argument stable while executing in tmp/MODULE/work."""
    import json
    from analysis_integrity import module_layout, resolve_cli_inputs
    layout = module_layout(output)
    command = list(command)
    def rewrite(arguments):
        arguments = list(arguments)
        for i, value in enumerate(arguments):
            if value in {'--out', '--dir-out', '--store'} and i+1 < len(arguments):
                if Path(arguments[i+1]).expanduser().absolute() == layout.output:
                    arguments[i+1] = str(layout.runtime)
            elif value == '--args-json' and i+1 < len(arguments):
                arguments[i+1] = json.dumps(rewrite(resolve_cli_inputs(json.loads(arguments[i+1]))))
            elif value.startswith(('--out=', '--dir-out=', '--store=')):
                flag, path = value.split('=', 1)
                if Path(path).expanduser().absolute() == layout.output:
                    arguments[i] = flag + '=' + str(layout.runtime)
        return arguments
    env = dict(os.environ, ALLELE_LAYOUT_ACTIVE='1', ALLELE_PUBLIC_OUT=str(layout.output),
               ALLELE_RUNTIME_OUT=str(layout.runtime), WORKFLOW_LOG_DIR=str(layout.logs))
    return resolve_cli_inputs(rewrite(command)), env, layout


def result_command(command):
    """Execution control does not change the scientific request."""
    result, index = [], 0
    while index < len(command):
        value = command[index]
        if value == '--replace':
            index += 2
            continue
        if value.startswith('--replace='):
            index += 1
            continue
        result.append(value)
        index += 1
    return result


def complete_layout(layout, command):
    """Reuse only a completed invocation with all recorded products present."""
    import json
    try:
        saved = json.loads((layout.qc / 'execution.json').read_text())
        if saved.get('status') != 'complete' or result_command(saved.get('command', [])) != result_command(command):
            return False
        inventory = saved.get('output_files')
        if not inventory:
            return False
        for name, size in inventory.items():
            path = Path(name)
            if not path.is_file() or path.stat().st_size != size:
                return False
        return True
    except (OSError, ValueError, TypeError):
        return False


def _execute_layout(command, output):
    import subprocess
    from analysis_integrity import atomic_json
    command, env, layout = layout_command(command, output)
    layout.qc.mkdir(parents=True, exist_ok=True)
    replace = os.environ.get('ALLELE_REPLACE', 'TRUE')
    if replace not in {'TRUE', 'FALSE'}:
        raise ValueError('--replace must be TRUE or FALSE')
    if replace == 'FALSE' and complete_layout(layout, command):
        print('[reuse] complete module: ' + str(layout.output), flush=True)
        return 0
    # Only this module's generated runtime is reset. Source/reference data and
    # the public result survive until successful publication. Partial results
    # always rebuild; stale receipts must not prevent a new attempt.
    if layout.runtime.is_symlink():
        raise ValueError('Module runtime must not be a symlink: ' + str(layout.runtime))
    input_flags = {'--eigenvec', '--coordinates', '--basis-dir', '--reuse-basis-dir', '--source-html'}
    reads_own_artifacts = any(flag in input_flags and i+1 < len(command)
        and Path(command[i+1]).expanduser().absolute().is_relative_to(layout.runtime)
        for i, flag in enumerate(command))
    from analysis_integrity import check_output_paths, InputPath
    protected = []
    output_flags = {'--out', '--dir-out', '--store', '--analysis-root', '--dir-analysis-root'}
    for flag, value in options(command).items():
        if flag in output_flags:
            continue
        path = Path(value).expanduser()
        if path.exists():
            path = path.resolve()
            if reads_own_artifacts and path.is_relative_to(layout.runtime):
                continue
            protected.append(InputPath(path, 'input'))
        if flag == '--pgen-prefix':
            protected.extend(InputPath(str(path) + suffix) for suffix in ('.pgen', '.pvar', '.psam'))
    check_output_paths(layout.output, protected)
    check_output_paths(layout.runtime, protected)
    if layout.runtime.exists() and not reads_own_artifacts:
        shutil.rmtree(layout.runtime)
    layout.runtime.mkdir(parents=True, exist_ok=True)
    env['ALLELE_REPLACE'] = 'TRUE'
    atomic_json(layout.qc / 'execution.json', dict(status='running', command=command,
                public_output=str(layout.output), runtime=str(layout.runtime)))
    code = 1
    try:
        code = subprocess.run(command, env=env).returncode
        return code
    finally:
        publish_layout(output, status='complete' if code == 0 else 'failed')
        products = {}
        if code == 0:
            for folder in (layout.runtime, layout.output):
                for path in folder.rglob('*'):
                    if path.is_file() and not path.is_symlink() and path.suffix not in {'.lock', '.log'}:
                        products[str(path)] = path.stat().st_size
        atomic_json(layout.qc / 'execution.json', dict(status='complete' if code == 0 else 'failed',
                    exit_code=code, command=command, public_output=str(layout.output),
                    runtime=str(layout.runtime), output_files=products,
                    validation='completion plus file presence and size; no content hashes'))


def execute_layout(command, output):
    import fcntl
    from analysis_integrity import module_layout
    layout = module_layout(output)
    layout.qc.mkdir(parents=True, exist_ok=True)
    with (layout.qc / 'execution.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Another worker is publishing this module: ' + str(output)) from exc
        return _execute_layout(command, output)


def layout_main(argv=None):
    parser = argparse.ArgumentParser(description='Run an analysis with separate public results and support directories')
    parser.add_argument('--module-out', required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command: parser.error('A worker command is required')
    return execute_layout(command, args.module_out)




def cli(argv=None):
    """Dispatch auxiliary commands owned by this module."""
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    routes = {'census': census_main, 'environment': environment_main, 'layout': layout_main}
    entry = routes[argv.pop(0)] if argv and argv[0] in routes else main
    return entry(argv)



if __name__ == '__main__':
    try: raise SystemExit(cli())
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        raise SystemExit(f'data-check: {exc}')

