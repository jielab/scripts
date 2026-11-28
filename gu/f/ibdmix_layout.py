"""Flat IBDmix outputs; only auxiliary metadata retains native read links."""
import os
import re
from pathlib import Path


def chrx_output_name(value):
    """Use one chrX prefix for current caller partitions."""
    tokens = [token for token in value.split('.')
              if token != 'chrX' and not re.fullmatch(r'(?:L[0-9]+_)?X_MALE', token)]
    value = '.'.join(tokens).lstrip('.')
    return 'chrX.' + value


def destination_for(name, unit):
    """One physical file per artifact; category directories have flat names."""
    p = Path(name)
    parts = p.parts

    def prefixed(items):
        value = '.'.join(items).lstrip('.')
        if unit == 'chrX':
            return chrx_output_name(value)
        return value if value.startswith(unit + '.') else unit + '.' + value

    if name == 'samples/ALL.txt':
        return Path('samples') / (unit + '.txt')
    if parts[0] == 'raw':
        rest = parts[1:]
        match = re.fullmatch(r'chr([0-9]+)', unit)
        if match and len(rest) > 1 and rest[0] == 'C' + match[1]:
            rest = rest[1:]
        value = prefixed(rest).replace('.YGB.raw.txt.gz', '.txt.gz').replace('.raw.txt.gz', '.txt.gz')
        return Path('raw') / value
    if parts[0] == 'qc' or '.qc.' in p.name:
        return Path('qc') / prefixed(parts[1:] if parts[0] == 'qc' else parts)
    if p.suffix in ('.log', '.err') or parts[0] in ('log', 'logs'):
        return Path('logs') / prefixed(parts[1:] if parts[0] in ('log', 'logs') else parts)
    if p.suffix == '.cmd' or p.name in ('run.meta.tsv', 'cache.meta.tsv', 'GU-MANIFEST.json', '.complete'):
        return Path('jobs') / prefixed(parts)
    # Keep only the native read-link hierarchy under tmp/ibdmix. Physical
    # intermediate files live directly in tmp, with their origin in the name.
    return Path('tmp') / prefixed([part for i, part in enumerate(parts) if i == len(parts)-1 or part != unit])


def unit_metadata_name(native_unit, run_unit):
    raw_name = destination_for(f"raw/{native_unit}/__metadata__.raw.txt.gz", run_unit).name
    return raw_name.removesuffix(".__metadata__.txt.gz") + ".json"


def publish_run_layout(run, cohort, replace=False):
    run, cohort = Path(run), Path(cohort)
    if not (run / '.complete').is_file() or not (run / 'run.meta.tsv').is_file():
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
        if source.relative_to(run).parts[0] not in ('raw',):
            source.symlink_to(os.path.relpath(target, source.parent))
    for category in ('raw',):
        directory = run / category
        if directory.is_dir():
            for parent, _, _ in os.walk(directory, topdown=False):
                if not any(Path(parent).iterdir()):
                    Path(parent).rmdir()
    return True


def cohort_for(run):
    run = Path(run)
    return run.parent.parent.parent if run.parent.name == 'ibdmix' and run.parent.parent.name == 'tmp' else None


def run_file(run, relative):
    """Raw calls have one physical published location, no read links."""
    run = Path(run)
    cohort = cohort_for(run)
    return cohort / destination_for(relative, run.name) if cohort else run / relative


def native_file(path):
    path = Path(path)
    for run in path.parents:
        if cohort_for(run):
            return run_file(run, path.relative_to(run).as_posix())
    return path


def raw_files(run):
    run = Path(run)
    cohort = cohort_for(run)
    return sorted((cohort / 'raw').glob(run.name + '.*.txt.gz')) if cohort else sorted((run / 'raw').glob('*/*.raw.txt.gz'))


def clear_calls(run):
    run = Path(run)
    cohort = cohort_for(run)
    if cohort:
        for category in ('raw',):
            for path in (cohort / category).glob(run.name + '.*'):
                if path.is_file() or path.is_symlink():
                    path.unlink()


if __name__ == '__main__':
    import sys
    if sys.argv[1] == 'path':
        print(run_file(sys.argv[2], sys.argv[3]))
    elif sys.argv[1] == 'existing-calls':
        print('yes' if raw_files(sys.argv[2]) else '')
    elif sys.argv[1] == 'clear-calls':
        clear_calls(sys.argv[2])
    else:
        raise SystemExit('Expected path RUN RELATIVE or clear-calls RUN')
