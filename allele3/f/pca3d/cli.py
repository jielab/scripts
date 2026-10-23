#!/usr/bin/env python3
"""Independent GRCh37 PCA workflow and offline PC1/PC2/age visualization.

Reuses calculation libraries without importing the allele3 workflow dispatcher.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pca_time import _atomic_text, _fingerprint, _rows, _sha, prepare_pca_plan, render_pca_time

PREFIX = "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes"


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("command", choices=("run", "plan", "render"), nargs="?", default="run")
    p.add_argument("--out", "--dir-out", type=Path, default=Path("/mnt/d/analysis/pca3d"))
    p.add_argument("--ancient-dir", "--dir-ancient", type=Path, default=Path("/mnt/f/gen/ancient"))
    p.add_argument("--pgen-prefix", type=Path)
    p.add_argument("--metadata", type=Path, help="Canonical TSV; required for custom panels / render")
    p.add_argument("--hq-variants", type=Path, help="GRCh37 background variant TSV")
    p.add_argument("--supplement-xlsx", type=Path, help="Original ancient dataset metadata adapter")
    p.add_argument("--qpadm", type=Path)
    p.add_argument("--ancestry-model", choices=("3way", "4way"))
    p.add_argument("--dir-pgen-cache", type=Path)
    p.add_argument("--eigenvec", type=Path, help="Existing common-basis .sscore / .eigenvec for render")
    p.add_argument("--html", type=Path)
    p.add_argument("--selected-variants", type=Path, help="Optional indexed PVAR targets for allele overlays")
    p.add_argument("--target-variants", type=Path, help="Optional separate GRCh37 background exclusions")
    p.add_argument("--pca-basis", choices=("joint_balanced", "modern_projection"), default="joint_balanced")
    p.add_argument("--pca-basis-label")
    p.add_argument("--pca-maf", type=float, default=.05)
    p.add_argument("--pca-geno", type=float, default=.1)
    p.add_argument("--pca-components", type=int, default=10)
    p.add_argument("--pca-max-per-group", type=int, default=200)
    p.add_argument("--pca-time-bin-years", type=int, default=2000)
    p.add_argument("--pca-target-flank-bp", type=int, default=1000000)
    p.add_argument("--pca-exclude-bed", type=Path)
    p.add_argument("--pca-threads", type=int, default=4)
    p.add_argument("--pca-memory-mb", type=int, default=8192)
    p.add_argument("--seed", type=int, default=20261021)
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--pca-max-variants", type=int, choices=range(1, 13), default=12)
    p.add_argument("--count-kind", choices=("hardcall", "dosage"), default="hardcall")
    p.add_argument("--title")
    p.add_argument("--simulated", action="store_true")
    return p


def metadata_path(a, prefix):
    if a.metadata:
        return a.metadata.expanduser().resolve()
    if a.command == "render":
        raise ValueError("render requires --metadata; no panel or ancient data are needed without overlays")
    if a.pgen_prefix and not a.supplement_xlsx:
        raise ValueError("Custom --pgen-prefix requires --metadata (or an explicit compatible --supplement-xlsx)")
    from ancient_metadata import prepare_metadata
    directory = a.out / "metadata"
    metadata = directory / "sample_metadata.tsv"
    receipt = directory / "pca3d_metadata.json"
    supplement = a.supplement_xlsx or a.ancient_dir / "41586_2026_10358_MOESM3_ESM.xlsx"
    psam = Path(str(prefix) + ".psam")
    inputs = {"psam": _fingerprint(psam), "supplement": _fingerprint(supplement),
              "qpadm": _fingerprint(a.qpadm) if a.qpadm else None, "ancestry_model": a.ancestry_model}
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved["inputs"] != inputs or any(not Path(path).is_file() or _sha(path) != digest
                for path, digest in saved["outputs"].items()):
            raise ValueError("Metadata inputs/output changed; use a new --out or supply verified --metadata")
        return metadata
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Existing metadata has no pca3d receipt; pass --metadata explicitly or use a new --out")
    result = prepare_metadata(psam, supplement, directory, qpadm_path=a.qpadm, ancestry_model=a.ancestry_model)
    outputs = {result[key]: _sha(result[key]) for key in ("metadata_path", "qc_path", "manifest_path")}
    _atomic_text(receipt, json.dumps(dict(inputs=inputs, outputs=outputs), indent=2) + "\n")
    return metadata


def overlay_options(a, prefix):
    if not a.selected_variants:
        return {}
    required = {"variant_id", "pvar_index", "chrom", "pos", "ref", "alt", "status"}
    selected, total, seen = [], 0, set()
    for row in _rows(a.selected_variants, required):
        total += 1
        if row["variant_id"] in seen:
            raise ValueError("Duplicate overlay variant_id: " + row["variant_id"])
        seen.add(row["variant_id"])
        if len(selected) < a.pca_max_variants:
            selected.append(row)
    if not selected:
        raise ValueError("--selected-variants is empty")
    folder = a.out / "pca"
    folder.mkdir(parents=True, exist_ok=True)
    bounded = folder / "overlay_targets.tsv"
    if bounded.resolve() == a.selected_variants.expanduser().resolve():
        raise ValueError("Overlay input must not be the generated overlay_targets.tsv; use a separate input file")
    with bounded.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(selected[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    _atomic_text(folder / "overlay_selection.json", json.dumps(dict(
        rule="first rows in explicit input order; bounded visualization only", total_available=total,
        displayed=len(selected), max_variants=a.pca_max_variants, source=_fingerprint(a.selected_variants)), indent=2) + "\n")
    return dict(pgen_path=Path(str(prefix) + ".pgen"), psam_path=Path(str(prefix) + ".psam"),
                selected_variants_path=bounded)


def workflow(a):
    if a.command == "render" and not a.eigenvec:
        raise ValueError("render requires --eigenvec with shared PCA coordinates")
    if a.command != "render" and a.eigenvec:
        raise ValueError("Use render for existing --eigenvec coordinates")
    prefix = (a.pgen_prefix or a.ancient_dir / PREFIX).expanduser().resolve()
    metadata = metadata_path(a, prefix)
    plan = None
    if a.command != "render":
        if a.pgen_prefix and not a.hq_variants:
            raise ValueError("Custom panels require --hq-variants; background QC cannot be inferred from a PGEN")
        if a.dir_pgen_cache:
            from input_cache import cache_pgen_prefix
            prefix, cache_record = cache_pgen_prefix(prefix, a.dir_pgen_cache)
            _atomic_text(a.out / "pca_input_cache.json", json.dumps(cache_record, indent=2) + "\n")
        plan = prepare_pca_plan(prefix, metadata, a.hq_variants or a.ancient_dir / "HighQualityVariants.tsv",
            a.out / "pca", basis=a.pca_basis, maf=a.pca_maf, geno=a.pca_geno,
            n_components=a.pca_components, max_per_group=a.pca_max_per_group,
            time_bin_years=a.pca_time_bin_years, seed=a.seed, exclude_bed=a.pca_exclude_bed,
            target_variants=a.target_variants or a.selected_variants, target_flank_bp=a.pca_target_flank_bp,
            plink2=a.plink2, threads=a.pca_threads, memory_mb=a.pca_memory_mb)
        if a.command == "plan":
            return dict(status="planned", manifest_path=plan["manifest_path"], script_path=plan["script_path"])
        subprocess.run(["bash", plan["script_path"]], check=True)
    source = plan["projected_path"] if plan else a.eigenvec
    label = a.pca_basis_label or (a.pca_basis if plan else "User-supplied shared PCA basis")
    result = render_pca_time(source, metadata, a.html or a.out / "pca" / "pca_time.html",
        basis_label=label, max_variants=a.pca_max_variants, simulated=a.simulated, title=a.title,
        count_kind=a.count_kind, **overlay_options(a, prefix))
    return dict(status="complete", **{k: result[k] for k in ("html_path", "manifest_path", "samples", "variants")})


def main(argv=None):
    a = parser().parse_args(argv)
    a.out = a.out.expanduser().resolve()
    a.ancient_dir = a.ancient_dir.expanduser().resolve()
    try:
        a.out.mkdir(parents=True, exist_ok=True)
        with (a.out / ".pca3d.lock").open("a") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ValueError("Another pca3d job is already running in this output directory") from None
            result = workflow(a)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"pca3d: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
