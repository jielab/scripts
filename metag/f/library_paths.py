"""Select and validate standalone current reference libraries."""
import json
from pathlib import Path

from common import fingerprint



def read_json(path):
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError('Expected a library metadata object: ' + str(path))
    return value


def checked_catalog(path):
    value = read_json(path)
    if value.get('status') != 'complete' or value.get('identity') != fingerprint(value.get('context')):
        raise RuntimeError('Invalid published catalog: ' + str(path))
    return value


def select_library(cfg, choice=None):
    cfg = dict(cfg)
    configured = Path(cfg.get('selected_library_root') or cfg['library_dir']).absolute()
    selected = Path(choice).expanduser() if choice else configured
    if choice and not selected.is_absolute():
        if len(selected.parts) != 1 or selected.name in ('.', '..'):
            raise ValueError('--library must be a sibling name (lib/lib0) or an absolute path')
        selected = configured.parent / selected
    selected = selected.absolute()
    alias = selected.name
    # Separate results for a nondefault library; frozen configs retain this scope.
    if cfg.get('selected_library_root') != str(selected):
        if choice and alias != 'lib':
            cfg['output_dir'] = str(Path(cfg['output_dir']).with_name('species.' + alias))
        cfg['work_dir'] = str(Path(cfg['work_dir']) / 'libraries' / alias)
    cfg.update(selected_library_root=str(selected), library_alias=alias,
               library_dir=str(selected), parallel_reference_dir=None,
               group_library_dir=str(selected / 'group'))
    manifest = selected / 'library.json'
    if manifest.exists():
        raise RuntimeError('Only standalone current libraries are supported; build a new library')
    cfg.update(db_dir=str(Path(cfg['library_dir']) / 'index'),
               prepared_dir=str(Path(cfg['library_dir']) / 'prepared'),
               taxonomy_dir=str(Path(cfg['library_dir']) / 'global/taxonomy'),
               database_layout='partitioned', storage_layout='separated')
    return cfg


def prepared_catalog(base):
    base = Path(base).absolute()
    value = checked_catalog(base / 'catalog.json')
    if Path(value['taxonomy_dir']) != base / 'global/taxonomy':
        raise RuntimeError('Prepared catalog paths changed; rebuild the library')
    return value


def union_snapshot(base, work_dir, *, full=False):
    from partitioned_build import validate_snapshot
    base = Path(base).absolute()
    value = checked_catalog(base / 'index/database.json')
    if Path(value['db_dir']) != base / 'index':
        raise RuntimeError('Index paths changed; rebuild the library')
    return validate_snapshot(value, full=full)
