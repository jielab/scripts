"""Output directories and table names for the ranked report."""
from common import locked

TABLES = ('samples', 'ranks', 'taxa', 'groups', 'evidence', 'transitions',
          'panels', 'panel_omitted', 'taxon_evidence')
OPTIONAL = ('target_species_evidence', 'fastq_qc', 'truth_metrics', 'truth_confusion')
NAMES = tuple(k + '.tsv' for k in TABLES + OPTIONAL)


def flat_ranked(path):
    return all((path / (k + '.tsv')).is_file() for k in TABLES)


def paper_root(root):
    if root.is_symlink():
        raise ValueError(f'Choose a real report directory with --out: {root}')
    root.mkdir(parents=True, exist_ok=True)
    for name in ('main', 'jobs', 'logs'):
        if (root / name).is_symlink():
            raise ValueError(f'Expected a real directory: {root/name}')
    with locked(root / 'logs/publish.lock'):
        (root / 'main').mkdir(exist_ok=True)
    return root
