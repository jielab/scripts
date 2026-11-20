"""Named group indexes inside one reference library; separate run outputs."""
from pathlib import Path
from library_paths import checked_catalog

INDEX_NAMES = ('legacy', 'optimized', 'full')


def select_index(cfg, name=None):
    cfg = dict(cfg)
    name = name or cfg.get('index_use', 'optimized')
    if name not in INDEX_NAMES:
        raise ValueError('--index-use must be legacy, optimized, or full')
    base = Path(cfg['library_dir']).expanduser().absolute()
    if (base / 'library.json').exists():
        raise RuntimeError('A standalone lib is required; the former lib/lib0 bundle must be consolidated first')
    work = Path(cfg.get('index_base_work_dir') or cfg['work_dir'])
    output = Path(cfg.get('index_base_output_dir') or cfg['output_dir'])
    cfg.update(index_use=name, index_base_work_dir=str(work), index_base_output_dir=str(output),
               library_dir=str(base), parallel_reference_dir=None, storage_layout='separated',
               database_layout='partitioned', db_dir=str(base/'index'),
               prepared_dir=str(base/'prepared'), taxonomy_dir=str(base/'global/taxonomy'),
               group_library_dir=str(base/'group'/('legacy' if name=='full' else name)),
               group_estimate_source=str(base/'group/legacy'),
               group_routing='full' if name=='full' else 'compact',
               work_dir=str(work/'indices'/name), output_dir=str(output/name))
    for key, path in [('source_preparation_identity',base/'catalog.json'),
                      ('source_full_index_identity',base/'index/database.json')]:
        if path.exists():
            identity=checked_catalog(path)['identity']
            if cfg.get(key,identity)!=identity:
                raise RuntimeError('Selected library changed since its configuration was frozen: '+key)
            cfg[key]=identity
    return cfg
