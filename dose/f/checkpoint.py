"""DOSE step checkpoints. Only successful, unchanged results may be skipped."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

SCHEMA = 1
DIETS = ('medito medi24 dash mind hpdi ahei phdi modern digm maha hlcd hlfd redii fds fudan upf').split()
PREPROCESS_FILES = ['diet0.rds', 'diet.rds', 'all.rds', 'diet.upf_2.rds'] + [f'diet.{x}.rds' for x in DIETS]
# Direct analysis dependencies; plotting dependencies also protect upstream tables.
DEPENDENCIES = {
    'foodwas/preprocess': [], 'foodwas/cohort': ['foodwas/preprocess'],
    'foodwas/patterns': ['foodwas/cohort'], 'foodwas/items': ['foodwas/cohort'],
    'foodwas/groups': ['foodwas/items'], 'foodwas/trend': ['foodwas/groups'],
    'foodwas/rcs': ['foodwas/groups'],
    'foodwas/plot': ['foodwas/items', 'foodwas/groups', 'foodwas/trend', 'foodwas/rcs'],
    'foodwas/plot_myth': ['foodwas/plot'],
    'ml/compare': ['foodwas/groups'], 'ml/importance': ['ml/compare'],
    'ml/select': ['ml/importance'], 'ml/cox': ['ml/select', 'foodwas/plot'],
    'ml/validate': ['ml/cox'], 'ml/bootstrap': ['ml/validate'], 'ml/plot': ['ml/bootstrap'],
    'score/build': ['foodwas/groups', 'ml/select'], 'score/distribution': ['score/build'],
    'score/stroke': ['foodwas/cohort', 'score/build'], 'score/sensitivity': ['score/stroke'],
    'score/plot': ['score/stroke', 'score/sensitivity', 'score/distribution'],
    'score/outcomes': ['foodwas/cohort', 'score/build'],
    'score/stroke_comparison': ['score/stroke', 'nhanes/stroke'],
    'mechanism/protein': ['foodwas/cohort', 'score/build'],
    'mechanism/pathways': ['mechanism/protein'], 'mechanism/network': ['mechanism/pathways'],
    'mechanism/pathway_plot': ['mechanism/network'],
    'mechanism/metabolite': ['foodwas/cohort', 'score/build'],
    'mechanism/metabolite_plot': ['mechanism/metabolite'],
    'mechanism/inflammation': ['foodwas/cohort', 'score/build'],
    'mechanism/inflammation_plot': ['mechanism/inflammation'],
    'mechanism/mri_prepare': ['foodwas/cohort'],
    'mechanism/mri': ['mechanism/mri_prepare', 'score/build'],
    'mechanism/mri_refit': ['mechanism/mri'], 'mechanism/mri_plot': ['mechanism/mri_refit'],
    'mechanism/mri_check': ['mechanism/mri_plot'], 'mechanism/mri_wmh': ['mechanism/mri'],
    'mechanism/sem': ['mechanism/metabolite', 'mechanism/inflammation'],
    'mechanism/sem_plot': ['mechanism/sem'], 'mechanism/sem_sensitivity': ['mechanism/sem'],
    'mechanism/summary': ['mechanism/sem_plot', 'mechanism/mri_plot', 'mechanism/pathway_plot',
                          'mechanism/network', 'mechanism/metabolite_plot', 'mechanism/inflammation_plot',
                          'signature/protein_cox', 'signature/metabolite_plot'],
    'signature/prepare': ['mechanism/protein', 'mechanism/metabolite', 'score/build'],
    'signature/protein_lasso': ['signature/prepare'],
    'signature/protein_cox': ['signature/protein_lasso', 'foodwas/cohort'],
    'signature/metabolite_lasso': ['signature/prepare'],
    'signature/metabolite_cox': ['signature/metabolite_lasso', 'foodwas/cohort'],
    'signature/metabolite_plot': ['signature/metabolite_cox'],
    'nhanes/score': ['foodwas/cohort', 'foodwas/groups', 'score/build'],
    'nhanes/stroke': ['nhanes/score'], 'nhanes/plot': ['nhanes/stroke'],
    'nhanes/flow': ['nhanes/stroke'],
}
# Legacy artifacts without a manifest cannot be assumed current. Globs are also
# used when a native cache returned successfully without rewriting its results.
MARKERS = {
    'foodwas/cohort': ['ukb/00data/dat.base*.rds', 'ukb/00data/analysis_config.rds'],
    'foodwas/patterns': ['ukb/01pattern/pattern_*.csv'],
    'foodwas/items': ['ukb/02foodwas/item/*.csv'],
    'foodwas/groups': ['ukb/02foodwas/dat.group.rds', 'ukb/03ml/ML_input.csv'],
    'foodwas/trend': ['ukb/02foodwas/group/Group_trend_*.csv'],
    'foodwas/rcs': ['ukb/02foodwas/group/RCS/RCS_summary.csv'],
    'foodwas/plot': ['ukb/02foodwas/group/Group_proper_reference.csv'],
    'foodwas/plot_myth': ['ukb/02foodwas/group/figures/food_group_MYTH_plot_data.csv'],
    'ml/compare': ['ukb/03ml/Model_Comparison_*.csv'],
    'ml/importance': ['ukb/03ml/Food_Importance.csv'],
    'ml/select': ['ukb/03ml/Primary_Selected_Features.csv'],
    'ml/cox': ['ukb/03ml/Cox_Imp_food_pattern.csv'],
    'ml/validate': ['ukb/03ml/SelectionValidation_*.csv'],
    'ml/bootstrap': ['ukb/03ml/*Bootstrap*.csv'],
    'ml/plot': ['ukb/03ml/figures/Fig_ML.pdf'],
    'score/build': ['ukb/04score/dose_score.rds'],
    'score/distribution': ['ukb/04score/figures/DOSE_food_intake_summary.csv'],
    'score/stroke': ['ukb/04score/pattern/pattern_continuous.csv'],
    'score/sensitivity': ['ukb/04score/IS_sensitivity/*.csv'],
    'score/outcomes': ['ukb/04score/otherOutcomes/*.csv'],
    'score/stroke_comparison': ['ukb/04score/overall_stroke_comparison/*.csv'],
    'score/plot': ['ukb/04score/pattern/Table_Pattern_Stroke_M1_M3.csv'],
    'mechanism/pathways': ['ukb/05mechanism/protein/M2/PPI_pathway/00_PPI_pathway_summary.csv'],
    'mechanism/network': ['ukb/05mechanism/protein/M2/PPI_pathway/STRING_network/**/00_network_summary.csv'],
    'mechanism/pathway_plot': ['ukb/05mechanism/protein/M2/PPI_pathway/Top_significant_pathways/Bubble_KEGG_Reactome_GO_BP/Displayed_terms.csv'],
    'mechanism/mri_prepare': ['ukb/05mechanism/04mri/inputs/mri_selected.tsv', 'ukb/05mechanism/04mri/inputs/mri_visit_covariates.csv'],
    'mechanism/mri': ['ukb/05mechanism/04mri/01_all_linear_associations.csv'],
    'mechanism/mri_refit': ['ukb/05mechanism/04mri/verification_refit.txt'],
    'mechanism/mri_plot': ['ukb/05mechanism/04mri/figures/main/Figure4_MRI.pdf'],
    'mechanism/mri_check': ['ukb/05mechanism/04mri/verification.json'],
    'mechanism/mri_wmh': ['ukb/05mechanism/04mri/figures/audit/WMH_adjusted_predictions.csv'],
    'mechanism/sem': ['ukb/05mechanism/sem/M2_Cox/three_model_comparison.csv'],
    'mechanism/sem_plot': ['ukb/05mechanism/sem/M2_Cox/**/plots/SEM_main_indicators.csv'],
    'mechanism/sem_sensitivity': ['ukb/05mechanism/sem/sensitivity/*/analysis_plan.csv'],
    'mechanism/summary': ['ukb/05mechanism/summary_figures/input_manifest.csv'],
    'signature/prepare': ['ukb/06signature/*/*_LASSO_input.csv.gz'],
    'signature/metabolite_plot': ['ukb/06signature/metabolite/figures_R/*.pdf'],
    'nhanes/score': ['.cache/nhanes/nhanes_*_score.rds'],
    'nhanes/stroke': ['nhanes/NHANES_analysis.rds'],
    'nhanes/plot': ['nhanes/Fig_NHANES_M2.png', 'nhanes/Fig_NHANES_M2_tertiles.png',
                    'nhanes/Sfig1_nhanes.png'],
    'nhanes/flow': ['nhanes/Sfig1_nhanes_primary.png'],
}
for domain in ('metabolite', 'inflammation'):
    MARKERS[f'mechanism/{domain}_plot'] = [f'ukb/05mechanism/{domain}/M2/Fig_DOSE_{domain}_mediation_beta_bubble_plot_data.csv']
for domain in ('protein', 'metabolite', 'inflammation'):
    MARKERS[f'mechanism/{domain}'] = [f'ukb/05mechanism/{domain}/M2/01a_DOSE_{domain}_XM_all.csv']
for domain in ('protein', 'metabolite'):
    MARKERS[f'signature/{domain}_lasso'] = [f'ukb/06signature/{domain}/{domain}_signature_scores.csv']
    MARKERS[f'signature/{domain}_cox'] = [f'ukb/06signature/{domain}/association_mediation/{domain}_signature_*.csv']


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + '.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as stream:
            json.dump(data, stream, sort_keys=True, ensure_ascii=False, indent=2)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding='utf-8'))
    except (FileNotFoundError, ValueError):
        return None


def stat_file(path):
    try:
        s = Path(path).stat()
        return {'size': s.st_size, 'mtime_ns': s.st_mtime_ns, 'ctime_ns': s.st_ctime_ns}
    except FileNotFoundError:
        return None


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def output_identity(path, expected):
    current = stat_file(path)
    if current is None:
        return None
    sha = expected['sha256'] if current == expected['stat'] else digest(path)
    return {'size': current['size'], 'sha256': sha}


def identity_changes(before, after, prefix=''):
    """Describe changed fields without printing environment values or input data."""
    if before == after:
        return []
    if isinstance(before, dict) and isinstance(after, dict):
        changes = []
        for key in sorted(before.keys() | after.keys()):
            path = f'{prefix} -> {key}' if prefix else str(key)
            if key not in before or key not in after:
                changes.append(path)
            else:
                changes.extend(identity_changes(before[key], after[key], path))
        return changes
    if isinstance(before, list) and isinstance(after, list) and len(before) == len(after):
        return [path for i, (old, new) in enumerate(zip(before, after))
                for path in identity_changes(old, new, f'{prefix} -> {i}')]
    return [prefix]


class Checkpoint:
    def __init__(self, key, command, output, ukb, scripts):
        if key not in DEPENDENCIES:
            raise ValueError(f'Unregistered checkpoint step: {key}')
        self.key, self.command = key, command
        self.output, self.ukb, self.scripts = map(lambda p: Path(p).resolve(), (output, ukb, scripts))
        # Script arguments define a separate invocation (e.g. plot-only).
        tag = hashlib.sha256(json.dumps(command).encode()).hexdigest()[:16]
        self.path = self.output / '.checkpoints' / key / (tag + '.json')

    def is_analysis_output(self, path):
        # CatBoost overwrites these training diagnostics across fits and steps.
        # Implicit R graphics devices can also leave a disposable root Rplots.pdf.
        # Also filter old manifests that accidentally registered these as results.
        try:
            relative = Path(path).relative_to(self.output)
        except ValueError:
            return True
        return 'catboost_info' not in relative.parts and relative != Path('Rplots.pdf')

    def markers(self):
        paths = set()
        for pattern in MARKERS.get(self.key, []):
            paths.update(p for p in self.output.glob(pattern) if p.is_file())
        if self.key == 'foodwas/preprocess':
            paths.update(self.ukb / 'Rdata' / name for name in PREPROCESS_FILES
                         if (self.ukb / 'Rdata' / name).is_file())
        return sorted(paths)

    def inventory(self):
        paths = []
        if self.output.exists():
            for root, dirs, files in os.walk(self.output):
                dirs[:] = [d for d in dirs if d not in {'.checkpoints', '.tools', '.r-lib', 'logs', 'tmp', '__pycache__', 'catboost_info'}
                           and not d.endswith('_cache_v2') and d != 'cache']
                paths.extend(Path(root) / name for name in files if not name.endswith(('.tmp', '.lock')))
        if self.key == 'foodwas/preprocess':
            paths.extend(self.markers())
        return {str(p): stat_file(p) for p in sorted(set(paths))
                if p.is_file() and self.is_analysis_output(p)}

    def identity(self):
        module = self.key.split('/')[0]
        code = [self.scripts / 'dose.sh', self.scripts / 'f/checkpoint.py',
                self.scripts / 'f/0.common.R', self.scripts / 'f/0.common.py']
        # R modules source definitions from other modules (including plots).
        code += list((self.scripts / 'f').glob('*.R')) + list((self.scripts / 'f').glob('*.py'))
        helper = Path(os.getenv('DOSE_HELPER_DIR', '/mnt/d/scripts/0f'))
        code += [helper / f for f in ('phenotype.R', 'association.R', 'plotting.R')]
        # RDS inputs are compressed and raw UKB files can be hundreds of GB.
        # Use size/mtime for inputs; generated output integrity additionally uses SHA256.
        inputs = list((self.ukb / 'Rdata').glob('*.rds'))
        if self.key == 'foodwas/preprocess':
            inputs = [p for p in inputs if p.name not in PREPROCESS_FILES]
        for directory in (self.ukb / 'common', Path(os.getenv('DOSE_FILES_DIR', '/mnt/d/files/diet'))):
            if directory.is_dir():
                inputs.extend(p for p in directory.rglob('*') if p.is_file())
        inputs += [self.ukb / 'rap/diet.tab.gz', self.ukb / 'rap/raw/pheno.tab.gz', self.ukb / 'rap/img.tab.gz']
        for variable in ('DOSE_MRI_RAW', 'DOSE_MRI_VISIT', 'DOSE_MRI_DICTIONARY',
                         'DOSE_NUTRIENT_FILE', 'DOSE_OUTCOME_DEFINITIONS'):
            if os.getenv(variable):
                inputs.append(Path(os.environ[variable]))
        nhanes = Path(os.getenv('NHANES_DATA', os.getenv('DOSE_FILES_DIR', '/mnt/d/files/diet')))
        if module == 'nhanes' and nhanes.is_dir():
            inputs.extend(p for p in nhanes.rglob('*') if p.is_file())
        upstream = {}
        ancestors = set()
        def visit(key):
            for dep in DEPENDENCIES[key]:
                if dep not in ancestors:
                    ancestors.add(dep)
                    visit(dep)
        visit(self.key)
        for dep in sorted(ancestors):
            states = sorted((self.output / '.checkpoints' / dep).glob('*.json'))
            upstream[dep] = []
            for path in states:
                state = read_json(path)
                if state:
                    upstream[dep].append({'status': state.get('status'),
                        'outputs': {p: output_identity(p, expected) for p, expected in state.get('outputs', {}).items()
                                    if self.is_analysis_output(p)}})
            if not states:  # Legacy upstream results still have to invalidate downstream caches.
                for pattern in MARKERS.get(dep, []):
                    upstream[dep].append({str(p): stat_file(p) for p in sorted(self.output.glob(pattern)) if p.is_file()})
        settings = {k: v for k, v in os.environ.items()
                    if (k.startswith('DOSE_') or k.startswith('NHANES_') or k.startswith('S4_'))
                    and k not in {'DOSE_REPLACE', 'DOSE_STEP'}}
        return {'schema': SCHEMA, 'key': self.key, 'command': self.command,
                'code': {str(p): digest(p) for p in sorted(set(code)) if p.is_file()},
                'inputs': {str(p): stat_file(p) for p in sorted(set(inputs))},
                'dependencies': upstream, 'settings': settings}

    def begin(self, replace=False):
        identity, previous = self.identity(), read_json(self.path)
        if previous and not replace:
            if previous.get('identity') != identity:
                raise RuntimeError('Inputs, code or settings changed; use --replace TRUE to rebuild this step.')
            if previous.get('status') == 'complete':
                missing = False
                for path, expected in previous['outputs'].items():
                    if not self.is_analysis_output(path):
                        continue
                    current = stat_file(path)
                    if current is None:
                        missing = True
                    elif current['size'] != expected['stat']['size'] or digest(path) != expected['sha256']:
                        raise RuntimeError(f'Output was modified: {path}; use --replace TRUE to rebuild.')
                if not missing:
                    return 'SKIP'
        if not previous and not replace:
            legacy = self.markers()
            status = self.output / 'logs' / (self.key + '.status')
            lines = status.read_text().splitlines() if status.is_file() else []
            was_complete = bool(lines) and lines[0] in {'COMPLETE', 'SKIPPED'}
            if legacy or was_complete:
                raise RuntimeError('Existing results have no verified checkpoint; use --replace TRUE once to rebuild and register them.')
        baseline = previous.get('baseline') if previous and previous.get('status') == 'running' and not replace else self.inventory()
        atomic_json(self.path, {'identity': identity, 'status': 'running', 'baseline': baseline,
                               'outputs': previous.get('outputs', {}) if previous else {}})
        return 'RUN'

    def finish(self):
        state = read_json(self.path)
        if not state or state.get('status') != 'running':
            raise RuntimeError('No matching running checkpoint.')
        changes = identity_changes(state['identity'], self.identity())
        if changes:
            details = '; '.join(changes[:8])
            if len(changes) > 8:
                details += f'; and {len(changes) - 8} more changed fields'
            raise RuntimeError('Inputs/code/settings or upstream results changed during analysis; '
                               f'checkpoint not committed. Changed: {details}. '
                               'Use --replace TRUE after resolving the change.')
        current = self.inventory()
        changed = {p for p, value in current.items() if state['baseline'].get(p) != value}
        changed.update(state.get('outputs', {}))
        changed.update(map(str, self.markers()))
        outputs = {p: {'stat': stat_file(p), 'sha256': digest(p)} for p in sorted(changed)
                   if self.is_analysis_output(p) and Path(p).is_file()}
        if not outputs:
            raise RuntimeError('Step produced no verifiable output; no success checkpoint was written.')
        atomic_json(self.path, {'identity': state['identity'], 'status': 'complete', 'outputs': outputs})
        print(f'Checkpoint saved: {self.key} ({len(outputs)} files)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['begin', 'finish'])
    parser.add_argument('module'); parser.add_argument('step')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    obj = Checkpoint(f'{args.module}/{args.step}', args.command, os.environ['DOSE_OUTDIR'],
                     os.environ['DOSE_UKB'], Path(__file__).resolve().parent.parent)
    try:
        if args.action == 'begin':
            print(obj.begin(os.getenv('DOSE_REPLACE', 'FALSE') == 'TRUE'))
        else:
            obj.finish()
    except (RuntimeError, OSError) as exc:
        print(f'Checkpoint {obj.key}: {exc}', file=sys.stderr)
        raise SystemExit(2)


if __name__ == '__main__':
    main()
