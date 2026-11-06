#!/usr/bin/env python3
"""Audit actual workflow inputs, then execute within the existing cgroup/log.

Bounded header/reference checks do not certify a complete dataset or silently
convert genome builds. Module-specific allele harmonization remains mandatory.
"""
from __future__ import annotations

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
        import herv_pipeline
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
            from herv_archaic import discover
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
    elif workflow.startswith('research-'):
        module = workflow.removeprefix('research-')
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
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args(argv)
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    if len(command) < 2: raise ValueError('Expected Python executable and workflow entry')
    from data_sources import load_config, audit_entries
    entries, fasta37, fasta38 = actual_sources(a.workflow, command[2:], load_config(a.config))
    if entries:
        report = audit_entries(entries, out=None, fasta37=fasta37, fasta38=fasta38,
                               max_records=1000, max_sites=100)
        if report.get('has_conflicts'):
            print('[build-check] conflicting source identities; workload was not started', file=sys.stderr, flush=True)
            return 2
    sys.stdout.flush(); sys.stderr.flush()
    os.execvp(command[0], command)


if __name__ == '__main__':
    try: raise SystemExit(main())
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        raise SystemExit(f'data-check: {exc}')
