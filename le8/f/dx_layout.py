"""Diagnosis-scoped kernel paths and the physical analysis directory layout."""
from pathlib import Path
import os

PRIMARY_DX = "icd"
DATES = ("icd", "ref", "icd1", "ref1", "any1")
READ_DATES = DATES
RUNTIME = ("logs", "jobs", "qc", "tmp")
MODULES = ("c1_correlate", "c2_cause", "c3_coloc", "c4_connect", "c5_cellulation", "pgs_focus")


def bind(link, target):
    target.mkdir(parents=True, exist_ok=True)
    link.parent.mkdir(parents=True, exist_ok=True)
    try:
        link.symlink_to(os.path.relpath(target, link.parent), target_is_directory=True)
    except FileExistsError:
        if not link.is_symlink() or link.resolve() != target.resolve():
            raise ValueError(f"Conflicting layout path: {link}")


def scope_for_read(layer, date):
    return Path(layer)/"main/results"/date


def init_scope(layer, date, modules):
    scope = layer/"main/results"/date
    (layer/"main/reports").mkdir(parents=True, exist_ok=True)
    for bucket in RUNTIME:
        bind(scope/bucket, layer/bucket/date)
    for name in modules:
        for bucket in RUNTIME:
            bind(scope/name/bucket, layer/bucket/date/name)
    return scope


def create_view(root, date, traits, layers, modules=MODULES, *, publication=False, deep=False):
    root = Path(root).absolute()
    view = root/"tmp"/("report_views" if publication else "dx_views")/date
    for trait in traits:
        for name in layers:
            layer = root/trait/name
            scope = scope_for_read(layer, date) if publication else init_scope(layer, date, modules)
            bind(view/trait/name, scope)
            if deep:
                prefix = Path("exploratory/abm_deep")
                layer = root/prefix/trait/name
                bind(view/prefix/trait/name, init_scope(layer, date, ("c1_correlate",)))
    for name in RUNTIME:
        bind(view/name, root/name/"by_date"/date)
    bind(view/"reports", root/"main/reports"/date)
    if publication:
        bind(view/"le8.paper", root/"paper/main/results/details"/date)
    return view


def bindings(root):
    """Validate only known kernel-view links before following them."""
    root = Path(root).absolute()
    if root.name not in (*READ_DATES, "default"):
        return
    kind = root.parent.name
    if (root.parent.parent.name, kind) not in {("tmp","dx_views"), ("tmp","report_views")}:
        return
    base = root.parents[2]
    for link in (root/"le8.paper", root/"reports", *root.glob("*/prot"), *root.glob("*/met"),
                 *root.glob("exploratory/abm_deep/*/prot"), *root.glob("exploratory/abm_deep/*/met")):
        if not link.is_symlink():
            continue
        relative = link.relative_to(root)
        if relative == Path("reports"):
            targets = [base/"main/reports"/root.name]
        elif relative == Path("le8.paper"):
            targets = [base/"paper/main/results/details"/root.name]
        else:
            layer = base/relative
            targets = [layer/"main/results"/root.name]
        if all(link.resolve() != target.resolve() for target in targets):
            raise ValueError(f"Unexpected result link: {link}")
        yield link


def result_files(root, pattern="*"):
    yield from root.rglob(pattern)
    for link in bindings(root):
        yield from link.rglob(pattern)


def contained(path, root):
    if path.resolve().is_relative_to(root.resolve()):
        return True
    return any(path.is_relative_to(link) and path.resolve().is_relative_to(link.resolve())
               for link in bindings(root))
