"""Keep source downloads separate from Kraken2 metadata and derived files."""
from pathlib import Path


DATA_DIRS = {
    'archaea': 'microbes/archaea',
    'bacteria': 'microbes/bacteria',
    'fungi': 'microbes/fungi',
    'protozoa': 'microbes/protozoa',
    'UniVec_Core': 'microbes',
    'viral': 'viruses',
    'plant': 'plants',
    'human': 'vertebrates/human',
    'vertebrate_mammalian': 'vertebrates/mammals',
    'vertebrate_other': 'vertebrates/other',
    'invertebrate': 'invertebrates',
}


def download_target(db, library, remote_path):
    db, remote_path = Path(db), Path(remote_path)
    if db.parent.name == '0.metadata' and library in DATA_DIRS:
        return db.parent.parent / DATA_DIRS[library] / remote_path.name
    return db / 'library' / library / remote_path


def source_path(filepath):
    """Map only source archives; processed FASTA and manifests stay in the DB."""
    path = Path(filepath).absolute()
    if not (path.name.endswith('_genomic.fna.gz') or path.name == 'UniVec_Core'):
        return str(path)
    for parent in path.parents:
        if parent.name != 'library':
            continue
        relative = path.relative_to(parent)
        if len(relative.parts) < 2:
            break
        return str(download_target(parent.parent, relative.parts[0], Path(*relative.parts[1:])))
    return str(path)
