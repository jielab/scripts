"""Discover published tables and reports without entering genotype or scheduler stores."""
import csv
import json
import os
from study_workflow import MODULE
from gu_segments import COHORTS
from pathlib import Path

MODULES = ('hotspots', 'density', 'herv_enrichment', 'base_occupancy', 'base_simulation',
           'enrich_test', 'allele_test', *MODULE.values())
EXTENSIONS = ('.tsv', '.tsv.gz', '.csv', '.bed', '.bed.gz', '.html', '.pdf', '.png', '.svg',
              '.json', '.sscore', '.eigenval')
SKIP = {'tmp', 'jobs', 'logs', 'qc', 'metadata', 'artifacts', 'shards', 'samples', 'scenarios',
        'variants', 'raw', 'segments', 'uncallable', 'basis', 'runtime', 'history', 'previous',
        'enrichment', 'study'}


def has_modules(path):
    return any((path / name).is_dir() for name in MODULES)


def envdir(name, default):
    return Path(os.environ.get(name) or default).expanduser().resolve()


def public_files(root):
    if not root.is_dir():
        return
    for parent, dirs, names in os.walk(root):
        dirs[:] = sorted(d for d in dirs if d not in SKIP and not d.startswith('.')
                         and not (Path(parent) / d).is_symlink())
        for name in sorted(names):
            if name.startswith('.') or not name.endswith(EXTENSIONS):
                continue
            path = Path(parent) / name
            if path.is_file():
                yield path


class Source:
    def __init__(self, key, root, *, target=None, reference=None, gu_cohort='', receipts=()):
        self.key, self.root = key, root.resolve()
        self.target, self.reference = target, reference
        self.gu_cohort = gu_cohort
        self.files = {}
        self.receipt = None
        self.meta = {}
        self.status = 'missing'
        self.problem = ''
        self.excluded_stages = set()
        for path in public_files(self.root):
            self.files[str(path.relative_to(self.root))] = path
        for name in ('metadata/sample_audit.tsv.gz', 'metadata/sample_metadata.tsv', 'AVAILABILITY.json',
                     'metadata/audit.tsv', 'metadata/calibration.tsv', 'metadata/manifest.json',
                     'metadata/pipeline_status.json', 'qc/manifest.json', 'qc/request.json'):
            path = self.root / name
            if path.is_file():
                self.files[name] = path
        # Native snapshots mirror public files. Only index scientific rows which
        # were published solely under artifacts, using the producer's own map.
        layout = self.root / 'metadata/layout.json'
        if layout.is_file():
            try:
                info = json.loads(layout.read_text())
            except (ValueError,OSError) as ex:
                info = {}
                self.problem = f'Unreadable layout: {ex}'
            for native, mapped in info.get('native_to_public', {}).items():
                path = self.root / mapped
                if (mapped.startswith('artifacts/') and native.startswith(('results/', 'reports/', 'figures/'))
                        and path.name.endswith(EXTENSIONS) and path.is_file()):
                    self.files.setdefault(native, path)
        names = (*receipts, 'metadata/manifest.json', 'qc/manifest.json', 'PRODUCT.json')
        for name in names:
            path = self.root / name
            if path.is_file():
                self.receipt = path
                try:
                    self.meta = json.loads(path.read_text())
                    if not isinstance(self.meta, dict):
                        raise ValueError('Expected an object')
                except (OSError, ValueError) as ex:
                    self.meta = {}
                    self.problem = f'Unreadable receipt {path}: {ex}'
                break
        if self.receipt or self.files:
            self.status = str(self.meta.get('status', 'files_only'))
        elif self.root.is_dir():
            self.status = 'no_public_products'
        build = self.meta.get('genome_build', self.meta.get('build'))
        if isinstance(build, (str, int)) and str(build).lower() not in ('38', 'b38', 'hg38', 'grch38'):
            self.problem = f'Build in receipt: {build}'
        for field, expected in [('target', target), ('ref', reference)]:
            if expected and self.meta.get(field) and self.meta[field] != expected:
                self.problem = f'{field} in receipt is {self.meta[field]}, expected {expected}'
        if self.status != 'files_only' and not self.status.startswith('complete'):
            self.problem = self.problem or self.status
        if key == 'herv_enrichment':
            stage_file = self.root / 'metadata/pipeline_status.json'
            if stage_file.is_file():
                try:
                    pipeline = json.loads(stage_file.read_text())
                    stages = pipeline.get('stages', [])
                    self.excluded_stages.update(pipeline.get('configuration', {}).get('disabled_stages', []))
                    if isinstance(stages, dict):
                        stages = [dict(v,stage=k) for k,v in stages.items()]
                except (OSError,ValueError,TypeError) as ex:
                    self.problem = f'Unreadable pipeline status: {ex}'
                    stages = []
                self.excluded_stages.update(s['stage'] for s in stages if s.get('status') not in ('complete','reused'))
        self.population = self.meta.get('population') or 'unrecorded'
        self.profile_disabled = False
        if reference == 'alphagenome' and self.receipt:
            request = self.receipt.with_name('request.json')
            if request.is_file():
                try:
                    self.profile_disabled = json.loads(request.read_text()).get('parameters',{}).get('profile_sites') == 0
                except (OSError,ValueError) as ex:
                    self.problem = f'Unreadable request: {ex}'

    def get(self, *names):
        if self.problem:
            return None
        for name in names:
            if '_exclude_mask.' in name and 'outputs' in self.meta and name not in self.meta['outputs']:
                continue
            if (name=='results/gradient_heterogeneity.tsv' and 'outputs' in self.meta
                    and name not in self.meta['outputs']):
                continue
            if self.profile_disabled and Path(name).name.startswith('avi_distance'):
                continue
            stage = Path(name).name.split('.')[0]
            if name in self.files and stage not in self.excluded_stages:
                return self.files[name]
        return None


def discover(args):
    analysis = args.analysis_root
    nested = analysis / 'main'
    default = envdir('ALLELE_RESULT_ROOT', nested)
    notes = []
    roots = [args.result_root or default]
    main = roots[0]

    def product(relative):
        return main / relative

    notes.append('Source roots, in priority order: ' + ', '.join(map(str, roots)) +
                 '. Each cohort/reference uses a single source product; see the module inventory.')
    scans = args.scan_root or envdir('ALLELE_SUMMARY_STORE', analysis / '3.allele3_scan')
    avi = args.avi_root or envdir('ALLELE_AVI_STORE', analysis / '4.avi_scan')
    gen = envdir('DIR_GEN', analysis.parent / 'gen')
    pca = args.pca_root or envdir('PCA3D_OUT', gen / 'ancient/38/pca3d')
    herv = args.herv_root or envdir('ALLELE_HERV_SCAN_STORE', analysis / '2.herv_scan')
    herv_scan = herv
    used={}
    inventory=main/'results/study_inputs.tsv'
    if inventory.is_file():
        with inventory.open() as handle:
            used={r['source']:Path(r['path']) for r in csv.DictReader(handle,delimiter='\t')
                  if r['status']=='ready' and Path(r['path']).exists()}
        if not args.pca_root:pca=used.get('pca',pca)
        if not args.avi_root:avi=used.get('avi',avi)
    sources = {}
    def add(key, path, **kwargs):
        sources[key] = Source(key, path, **kwargs)
    add('pca', pca)
    add('pca_study', envdir('PCA3D_STUDY_OUT', pca / 'study'), receipts=('study.json',))
    add('herv_scan', herv, receipts=('qc/summary.json',))
    add('avi_scan', avi)
    add('detailed_scan', used.get('detailed',envdir('ALLELE_SCAN_STORE', scans / 'gnomad_detailed')))
    for key in ('hotspots', 'density', 'herv_enrichment', 'base_simulation'):
        add(key, product(key))
    for key in ('space', 'allele_time'):
        add(key, product(Path('base_occupancy') / key))
    occupancy = {path.name for root in roots if (root / 'base_occupancy').is_dir()
                 for path in (root / 'base_occupancy').iterdir() if path.is_dir()}
    for name in sorted(occupancy - SKIP - {'space','allele_time'}):
        if not name.startswith('.'):
            add('occupancy/' + name, product(Path('base_occupancy') / name))
    add('study_inventory',main/'results')
    add('archaic_carriers',product(Path('archaic_carrier_genotypes')/'1kg'),target='1kg',reference='archaic',receipts=('availability.json',))
    for target in args.targets:
        add('scan/' + target, scans/target if args.scan_root else used.get(target,scans/target), target=target)
        add('allele_test/' + target, product(Path(MODULE['alleles']) / target), target=target)
        for cohort in COHORTS:
            add('archaic/'+cohort+'/'+target,product(Path(MODULE['archaic'])/cohort/target),
                target=target,reference='archaic',gu_cohort=cohort)
        for reference in ('alphagenome', 'herv'):
            path = product(Path(MODULE[reference]) / target)
            add(reference + '/' + target, path, target=target, reference=reference)
        for ref, module in MODULE.items():
            folders={p.name for root in roots if (root/module).is_dir()
                     for p in (root/module).glob(target+'_*') if p.is_dir()}
            for folder in sorted(folders):
                if folder==target+'_mean_adjusted' or folder.startswith(target+'_shared_af_'):
                    add('study/'+ref+'/'+folder,product(Path(module)/folder),target=target,reference=ref)
    status_file=main/'results/study_status.tsv'
    if status_file.is_file():
        with status_file.open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                if row['design']=='core':key=row['module']
                elif row['design']=='primary' and row['reference']=='archaic':key='archaic/'+row['gu_cohort']+'/'+row['target']
                elif row['design']=='primary':key=('allele_test' if row['reference']=='alleles' else row['reference'])+'/'+row['target']
                elif row['design']=='individual_followup':key='archaic_carriers'
                else:key='study/'+row['reference']+'/'+Path(row['out']).name
                if key in sources and row['status']!='complete':
                    sources[key].status=row['status']
                    sources[key].problem=row['reason'] or 'Unavailable in this run'
    return main, sources, notes
