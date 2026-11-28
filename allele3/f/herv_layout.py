"""Shared HERV analysis layout; the module root contains directories only."""
from pathlib import Path


def product_path(out, name):
    if Path(name).name != name:
        raise ValueError('Expected a HERV product basename: '+name)
    folder = 'reports' if name.endswith('.html') else 'qc' if name.endswith('.json') or name == 'sample_completion.tsv' else 'results'
    return Path(out)/folder/name


def ensure_layout(out):
    for name in ('results','reports','figures','qc','jobs','logs','tmp'):
        (Path(out)/name).mkdir(parents=True,exist_ok=True)
