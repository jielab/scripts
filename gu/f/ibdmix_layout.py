"""Unarchived IBDmix publication; native paths remain as lightweight links."""
import os
import re
from pathlib import Path


def destination_for(name, unit):
    """One physical file per artifact; category directories have flat names."""
    unit = re.sub(r'^chrX\.nonpar-v2\.mac[0-9]+$', 'chrX', unit)
    p = Path(name)
    parts = p.parts

    def prefixed(items):
        value = '.'.join(items).lstrip('.')
        return value if value.startswith(unit + '.') else unit + '.' + value

    if name == 'ibdmix.tracts.xlsx':
        return Path('workbooks') / (unit + '.tracts.xlsx')
    if name == 'samples/ALL.txt':
        return Path('samples') / (unit + '.txt')
    if parts[0] == 'raw':
        rest = parts[1:]
        match = re.fullmatch(r'chr([0-9]+)', unit)
        if match and len(rest) > 1 and rest[0] == 'C' + match[1]:
            rest = rest[1:]
        return Path('raw') / prefixed(rest)
    if parts[0] == 'segments':
        value = '.'.join(parts[1:])
        # Native C1 and X_MALE IDs already disambiguate whole chromosomes.
        call = 'C' + unit[3:] if re.fullmatch(r'chr[0-9]+', unit) else 'L00001_X_MALE' if unit == 'chrX' else ''
        if not call or not value.startswith(call + '.'):
            value = prefixed(parts[1:])
        return Path('segments') / value
    if parts[0] == 'final':
        value = '.'.join(parts[1:])
        match = re.fullmatch(r'all_archaic_refs\.lod[0-9.]+\.(len[0-9]+)\.segments\.tsv\.gz', value)
        if match:
            value = match[1] + '.tsv.gz'
        return Path('lod') / prefixed([value])
    if parts[0] == 'qc' or '.qc.' in p.name:
        return Path('qc') / prefixed(parts[1:] if parts[0] == 'qc' else parts)
    if p.suffix in ('.log', '.err') or parts[0] in ('log', 'logs'):
        return Path('logs') / prefixed(parts[1:] if parts[0] in ('log', 'logs') else parts)
    if p.suffix == '.cmd' or p.name in ('run.meta.tsv', 'cache.meta.tsv', 'GU-MANIFEST.json', '.complete'):
        return Path('jobs') / prefixed(parts)
    # Keep only the native read-link hierarchy under tmp/ibdmix. Physical
    # intermediate files live directly in tmp, with their origin in the name.
    return Path('tmp') / prefixed([part for i, part in enumerate(parts) if i == len(parts)-1 or part != unit])


def publish_run_layout(run, cohort, replace=False):
    run, cohort = Path(run), Path(cohort)
    if not (run / 'final').is_dir() or not (run / 'run.meta.tsv').is_file():
        return False
    files = [p for p in run.rglob('*') if p.is_file() and not p.is_symlink()]
    destinations = {}
    for source in files:
        if source.name.endswith('.raw.tar.gz'):
            raise ValueError(f'Unexpected legacy archive: {source}')
        target = cohort / destination_for(source.relative_to(run).as_posix(), run.name)
        if target in destinations:
            raise ValueError(f'Output collision: {target}')
        if target.exists() and not replace:
            raise FileExistsError(target)
        destinations[target] = source
    for target, source in destinations.items():
        target.parent.mkdir(parents=True, exist_ok=True)
        os.replace(source, target)
        source.symlink_to(os.path.relpath(target, source.parent))
    return True
