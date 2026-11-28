"""Direct central support paths without filesystem compatibility aliases."""
from pathlib import Path
import os
import fcntl
from product_paths import support as product_support, permanent

MULTI_TARGET = { 'enrich_test', 'allele_test'}


def support_path(output, kind):
    output = Path(output).absolute()
    product = product_support(output, kind)
    if product is not None: return product
    module = next((p for p in output.parents if p.name in MULTI_TARGET), None)
    return module/kind/output.relative_to(module) if module else output/kind


def data_path(output, kind):
    """Direct destination for scientific QC or centralized runtime; no aliases."""
    output = Path(output).expanduser().absolute()
    if not permanent(output):
        return output / kind
    if kind == 'qc' and not support_path(output, 'tmp').is_relative_to(output):
        return output / 'metadata'
    return support_path(output, kind)


def prepare(output):
    output = Path(output).expanduser().absolute()
    if not permanent(output):
        return
    output.mkdir(parents=True, exist_ok=True)
    lockdir = support_path(output, 'qc')
    lockdir.mkdir(parents=True, exist_ok=True)
    with (lockdir / 'layout.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for kind in ('jobs', 'logs', 'qc', 'tmp'):
            source, destination = output / kind, data_path(output, kind)
            destination.mkdir(parents=True, exist_ok=True)
            if source == destination:
                continue
            if source.is_symlink() or source.exists():
                raise ValueError('Unexpected noncanonical support path: ' + str(source))


# Kept as an internal callable name for existing importers.
prepare_product = prepare


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    p.add_argument('--permanent-only', action='store_true')
    p.add_argument('--kind', choices=('jobs','logs','qc','tmp'))
    a = p.parse_args()
    if a.kind: print(support_path(a.output, a.kind) if not a.permanent_only or permanent(a.output) else '')
    else: prepare(a.output)
