"""Select library bundles and read moved legacy libraries without rewriting them.

A bundle owns its group index and refers to an immutable prepared/full library.
Catalog identities and file baselines remain those of the original publication.
Only runtime paths and a tiny union hash list are relocated.
"""
import copy
import json
from pathlib import Path

from common import atomic_json, fingerprint, locked, sha256

SCHEMA = 'metag-library-bundle-v1'


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


def rebase(path, original, actual):
    p, original, actual = Path(path), Path(original), Path(actual)
    try:
        relative = p.relative_to(original)
    except ValueError:
        raise RuntimeError(f'Catalog path is outside its library: {p}') from None
    if '..' in relative.parts:
        raise RuntimeError('Unsafe relative library path')
    return str(actual / relative)


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
        bundle = read_json(manifest)
        if bundle.get('schema') != SCHEMA or bundle.get('identity') != fingerprint(bundle.get('context')):
            raise RuntimeError('Invalid library bundle: ' + str(manifest))
        source = Path(bundle['context']['source_library'])
        catalog = checked_catalog(source / 'catalog.json')
        if catalog['identity'] != bundle['context']['preparation_identity']:
            raise RuntimeError('Bundle source preparation changed')
        full = checked_catalog(source / 'index/database.json')
        if full['identity'] != bundle['context']['full_index_identity']:
            raise RuntimeError('Bundle source full index changed')
        cfg['library_dir'] = str(source)
        cfg['library_bundle_identity'] = bundle['identity']
        cfg['source_preparation_identity'] = catalog['identity']
        cfg['source_full_index_identity'] = full['identity']
    cfg.update(db_dir=str(Path(cfg['library_dir']) / 'index'),
               prepared_dir=str(Path(cfg['library_dir']) / 'prepared'),
               taxonomy_dir=str(Path(cfg['library_dir']) / 'global/taxonomy'),
               database_layout='partitioned', storage_layout='separated')
    return cfg


def create_bundle(destination, source):
    destination, source = Path(destination).absolute(), Path(source).absolute()
    if destination == source or destination in source.parents or source in destination.parents:
        raise ValueError('New and source libraries must be separate, non-nested directories')
    if destination.is_symlink() or source.is_symlink():
        raise ValueError('Library roots must not be symlinks')
    prepared = checked_catalog(source / 'catalog.json')
    full = checked_catalog(source / 'index/database.json')
    context = dict(source_library=str(source), preparation_identity=prepared['identity'],
                   full_index_identity=full['identity'])
    value = dict(schema=SCHEMA, context=context, identity=fingerprint(context))
    destination.mkdir(parents=True, exist_ok=True)
    with locked(destination / '.bundle.lock'):
        marker = destination / 'library.json'
        if marker.exists() and read_json(marker) != value:
            raise RuntimeError('Existing bundle refers to different inputs; choose another --library')
        if not marker.exists():
            if any(p.name not in ('.bundle.lock',) for p in destination.iterdir()):
                raise RuntimeError('Refusing to adopt a nonempty library directory')
            atomic_json(marker, value)
    return value


def prepared_catalog(base):
    base = Path(base).absolute()
    value = checked_catalog(base / 'catalog.json')
    original = Path(value['taxonomy_dir']).parents[1]
    result = copy.deepcopy(value)
    for group in result['groups'].values():
        for row in group['inputs']:
            row['path'] = rebase(row['path'], original, base)
        if group.get('inputs_tsv'):
            group['inputs_tsv'] = rebase(group['inputs_tsv'], original, base)
    result['taxonomy_dir'] = rebase(value['taxonomy_dir'], original, base)
    result['seqid_map'] = rebase(value['seqid_map'], original, base)
    return result


def union_snapshot(base, work_dir, *, full=False):
    """Make a checked runtime view; never rewrite lib0 catalogs or index bytes."""
    base = Path(base).absolute()
    marker = base / 'index/database.json'
    original = checked_catalog(marker)
    if Path(original['db_dir']) == base / 'index':
        from partitioned_build import validate_snapshot
        return validate_snapshot(original, full=full)
    result = copy.deepcopy(original)
    old_root = Path(original['db_dir']).parent
    result['db_dir'] = str(base / 'index')
    result['taxonomy_dir'] = rebase(original['taxonomy_dir'], old_root, base)
    result['_catalog_path'] = str(marker)
    result['_catalog_sha256'] = sha256(marker)
    result['_runtime_work_dir'] = str(Path(work_dir).absolute())
    hashes = ''.join(str(base / 'index' / name) + '\n' for name in result['hashes'])
    folder = Path(work_dir) / 'cache/library-paths' / fingerprint(dict(catalog=result['_catalog_sha256'],root=str(base)))
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / 'hashes.txt'
    with locked(folder / 'paths.lock'):
        if not path.exists():
            tmp = folder / 'hashes.tmp'
            tmp.write_text(hashes)
            tmp.replace(path)
        elif path.read_text() != hashes:
            raise RuntimeError('Runtime union hash list changed: ' + str(path))
    result['_hashes_path'] = str(path)
    return validate_relocated_union(result, full=full)


def validate_relocated_union(value, *, full=False):
    from partitioned_build import CATEGORIES
    from reference import _native_index_structure
    marker = Path(value['_catalog_path'])
    if sha256(marker) != value['_catalog_sha256']:
        raise RuntimeError('Relocated source catalog changed')
    saved = checked_catalog(marker)
    if saved.get('schema') != 'metag-partitioned-index-v1' or saved.get('layout') != 'partitioned':
        raise RuntimeError('Invalid partitioned library schema/layout')
    context = saved['context']
    for name in ('layout','db_dir','taxonomy_dir','hashes','files','taxonomy_files','parameters','verification_path'):
        if saved.get(name) != context.get(name):
            raise RuntimeError('Inconsistent source catalog: ' + name)
    categories = context.get('categories', {})
    expected = [name + '/hash.k2d' for name in CATEGORIES if name in categories]
    if not categories or set(categories) - set(CATEGORIES) or saved['hashes'] != expected:
        raise RuntimeError('Invalid category hash inventory')
    if sorted(context.get('absent_categories', [])) != sorted(set(CATEGORIES)-set(categories)):
        raise RuntimeError('Invalid absent-category inventory')
    required = {'taxo.k2d','opts.k2d','hashes.txt'} | {c+'/'+f for c in categories for f in ('hash.k2d','taxo.k2d','opts.k2d')}
    if set(saved['files']) != required or not {'nodes.dmp','names.dmp','merged.dmp'} <= set(saved['taxonomy_files']):
        raise RuntimeError('Incomplete source catalog inventory')
    db = marker.parent
    tax = rebase(saved['taxonomy_dir'], Path(saved['db_dir']).parent, db.parent)
    expected_view = copy.deepcopy(saved)
    expected_view.update(db_dir=str(db), taxonomy_dir=tax)
    for key in ('_catalog_path','_catalog_sha256','_runtime_work_dir','_hashes_path'):
        expected_view[key] = value[key]
    if expected_view != value:
        raise RuntimeError('Runtime library view differs from its source catalog')
    for directory, records in ((db,saved['files']), (Path(tax),saved['taxonomy_files'])):
        for name, record in records.items():
            if Path(name).is_absolute() or '..' in Path(name).parts:
                raise RuntimeError('Unsafe library file path')
            path = directory / name
            if path.is_symlink() or not path.is_file() or path.stat().st_size != record['state']['size']:
                raise RuntimeError('Missing or truncated library file: ' + str(path))
            if full:
                before = path.stat()
                digest = sha256(path)
                after = path.stat()
                if (before.st_size, before.st_mtime_ns, before.st_ctime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino):
                    raise RuntimeError('Library changed during verification: ' + str(path))
                if digest != record['sha256']:
                    raise RuntimeError('Library SHA256 mismatch: ' + str(path))
    original_list = ''.join(str(Path(saved['db_dir'])/n)+'\n' for n in saved['hashes'])
    runtime_list = ''.join(str(db/n)+'\n' for n in saved['hashes'])
    if (db/'hashes.txt').read_text() != original_list or Path(value['_hashes_path']).read_text() != runtime_list:
        raise RuntimeError('Union hash list differs from catalog')
    for name in saved['hashes']:
        _native_index_structure((db/name).parent)
    return value
