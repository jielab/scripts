"""Numbered physical analysis directories shared by all GRID entry points."""
from pathlib import Path

MODULES = {"ld": "0.ld_uni", "pca": "1.pca", "csx": "2.csx",
           "disco": "3.disco", "prsformer": "4.prsformer"}
TRAITS = {"height", "ldl", "t2dm"}
CSX_FILES = {"scores.tsv.gz", "baseline.predictions.tsv.gz", "predictions.tsv.gz",
             "benchmark.json", "benchmark_models.joblib", "oof_folds.tsv.gz"}


def upstream_root(root):
    root = Path(root)
    if root.name == "main":
        return root.parent
    return root.parent.parent if root.parent.name == "main" else root


def module_root(root, entry):
    root = Path(root)
    return root if entry not in MODULES else upstream_root(root) / MODULES[entry]


def module_results(root, entry):
    """Return the canonical results directory for each numbered module."""
    folder = module_root(root, entry)
    if entry == "csx":
        return folder
    return folder / "main"


def project_path(root, relative):
    """Map each project-relative artifact to its owning workflow directory."""
    root, relative = Path(root), Path(relative)
    parts = relative.parts
    if not parts or relative.anchor or ".." in parts:
        raise ValueError("Expected a safe project-relative output path")
    if parts[0] == "grid.paper":
        destination = upstream_root(root) / "paper"
        run = (root.name,) if root.parent.name == "main" else ()
        if len(parts)>1 and parts[1] in {"jobs","logs","main"}:
            kind = "main" if parts[1]=="main" else "logs"
            return destination.joinpath(kind,*run,*parts[2:])
        main = destination.joinpath("main", *run)
        if len(parts) > 1 and parts[1] == "traits":
            return main.joinpath(*parts[2:])
        if len(parts)==2:
            filename=parts[1]
            trait=filename.split(".",1)[0]
            if trait in TRAITS and "." in filename:
                return main/trait/filename
            if filename.startswith(("table_","all.")) and filename.endswith(".tsv"):
                return main/"tables"/filename
            if filename.endswith(".png"):
                return main/"figures"/filename
        return main.joinpath(*parts[1:])
    entry = None
    if parts[0] == "tmp" and len(parts) > 1:
        name = parts[1]
        if name == "pca": entry = "pca"
        elif name in {"cohort", "gwas_inputs", "gwas_manifest.tsv"} or name.startswith("gwas_manifest."):
            return module_results(root, "csx").joinpath(*parts[1:])
        elif name == "prsformer":
            return module_root(root, "prsformer").joinpath("tmp", *parts[2:])
        elif name in TRAITS:
            tail = parts[2:]
            if not tail:
                return root / name
            item = tail[0]
            if item == "disco":
                return module_results(root, "disco").joinpath(name, *tail[1:])
            if item == "predictions.tsv.gz":
                return module_results(root, "disco") / name / item
            if item.startswith("prsformer"):
                return module_root(root, "prsformer") / "main" / name / item
            if item == "inference":
                rest = tail[1:]
                if rest and rest[0] == "partitions": rest = ("chromosomes", *rest[1:])
                elif rest and rest[0] == "inference": rest = ("raw", *rest[1:])
                return module_results(root, "csx").joinpath(name, *rest)
            if item in CSX_FILES:
                return module_results(root, "csx").joinpath(name, *tail)
            return root.joinpath(name, *tail)
    elif parts[0] == "qc" and len(parts) > 1:
        name = parts[1]
        if name == "pca" or name.startswith("pca."): entry = "pca"
        elif name.startswith("prsformer"): entry = "prsformer"
        elif ".disco_" in name: entry = "disco"
        elif name.startswith(("cohort.", "gwas_inputs")) or name.endswith((".benchmark_tuning.tsv", ".coefficients.tsv")): entry = "csx"
    if entry is None and root.name == "main" and parts[0] in {"jobs", "logs", "tmp", "qc"}:
        return upstream_root(root).joinpath(*parts)
    if entry is None and root.parent.name == "main" and parts[0] in {"jobs", "logs", "tmp", "qc"}:
        return upstream_root(root).joinpath(parts[0], root.name, *parts[1:])
    return module_root(root, entry) / relative


def trait_path(root, trait, relative):
    return project_path(root, Path("tmp") / str(trait) / relative)
