#!/usr/bin/env python3
"""CIGMA-HE adapter, verified against Minhui-Chen/CIGMA 5813e4a.

Input is one gene per manifest row; ctp, ctnu, P, K are labelled CSV matrices.
No GWAS P value is treated as an expression variance or colocalization posterior.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPOSITORY = "https://github.com/Minhui-Chen/CIGMA"
VERIFIED_COMMIT = "5813e4ae84d7b3733dfcd938fe42d12c6b30a8aa"


def matrix(path):
    # Preserve donor identifiers such as '001'; never align by row position.
    with Path(path).open(newline="") as stream:
        header = next(csv.reader(stream), [])
    # pandas otherwise silently renames duplicate column headers.
    if len(header) < 2 or len(set(header)) != len(header) or any(not h.strip() for h in header[1:]):
        raise ValueError(f"Missing/duplicate matrix column labels: {path}")
    d = pd.read_csv(path, dtype=str)
    if d.shape[1] < 2:
        raise ValueError(f"Expected labelled matrix: {path}")
    d = d.set_index(d.columns[0])
    if (d.index.has_duplicates or d.columns.has_duplicates or d.index.isna().any()
            or any(not str(v).strip() for v in d.index)):
        raise ValueError(f"Missing/duplicate identifiers: {path}")
    d = d.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(d.to_numpy()).all():
        raise ValueError(f"Nonfinite values: {path}; missing donor/cell pairs need preprocessing")
    return d


def aligned(path, donors, columns=None):
    d = matrix(path)
    if set(d.index) != set(donors):
        raise ValueError(f"Donor set mismatch: {path}")
    if columns is not None:
        if set(d.columns) != set(columns):
            raise ValueError(f"Column/cell-type set mismatch: {path}")
        d = d.loc[:, columns]
    return d.loc[donors]


def read_inputs(row, base):
    if row.get("ctnu_definition") != "variance_of_pseudobulk_mean":
        raise ValueError("ctnu_definition must be variance_of_pseudobulk_mean: "
                         "upstream preprocess.pseudobulk returns SEM squared, not cell variance or SD")
    paths = {k: (base / row[k]).resolve() for k in ("ctp", "ctnu", "P", "K")}
    y = matrix(paths["ctp"])
    if len(y) < 20 or y.shape[1] < 2:
        raise ValueError("CIGMA adapter requires >=20 donors and >=2 cell types")
    donors, cells = y.index, y.columns
    nu = aligned(paths["ctnu"], donors, cells)
    prop = aligned(paths["P"], donors, cells)
    kin = aligned(paths["K"], donors, donors)
    if (nu.to_numpy() < 0).any() or (prop.to_numpy() < 0).any():
        raise ValueError("Negative noise variance or cell proportion")
    # Selected cell types need not exhaust the source tissue, but totals <= 1.
    sums = prop.sum(axis=1).to_numpy()
    if (sums <= 0).any() or (sums > 1 + 1e-6).any():
        raise ValueError("Invalid cell proportions; values must sum to (0,1]")
    k = kin.to_numpy()
    if not np.allclose(k, k.T, atol=1e-7) or np.linalg.eigvalsh(k).min() < -1e-6:
        raise ValueError("Kinship is not symmetric positive semidefinite")
    if np.allclose(k, np.eye(len(k)) * np.trace(k) / len(k)):
        raise ValueError("Identity kinship cannot separate genetic and environmental covariance")
    if np.diag(k).min() <= 0 or np.var(y.to_numpy(), axis=0).min() <= 0:
        raise ValueError("Degenerate kinship or unexpressed/constant cell type")
    args = dict(Y=y.to_numpy(), K=k, ctnu=nu.to_numpy(), P=prop.to_numpy())
    for name in ("fixed_covars", "random_covars"):
        if row.get(name):
            path = (base / row[name]).resolve()
            cov = aligned(path, donors)
            args[name] = {"design": cov.to_numpy()}
            paths[name] = path
    return args, list(cells), paths


def bh(values):
    p = np.asarray(values, dtype=float)
    if ((p[np.isfinite(p)] < 0) | (p[np.isfinite(p)] > 1)).any():
        raise ValueError("P values must lie in [0,1]")
    result = np.full(len(p), np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    ix = ok[np.argsort(p[ok])]
    if len(ix):
        # Keep failed/unavailable genes in the declared manifest family.
        result[ix] = np.minimum(1, np.minimum.accumulate((p[ix] * len(p) / np.arange(1, len(ix)+1))[::-1])[::-1])
    return result


def run(manifest, outdir, validate_only=False, seed=2026):
    outdir.mkdir(parents=True, exist_ok=True)
    with manifest.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or not {"gene", "ctp", "ctnu", "P", "K", "tissue", "build", "kinship_scope", "ctnu_definition"} <= rows[0].keys():
        raise ValueError("Manifest requires gene, ctp, ctnu, P, K, tissue, build, kinship_scope, ctnu_definition")
    keys = [(r["gene"], r["tissue"]) for r in rows]
    if any(not r[k].strip() for r in rows for k in ("gene", "tissue", "build", "kinship_scope")):
        raise ValueError("Manifest provenance fields must not be blank")
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate gene/tissue jobs")
    version = "not loaded (validation only)"
    if not validate_only:
        from cigma import fit
        version = importlib.metadata.version("cigma")
    status, results, cells_out, provenance = [], [], [], []
    for row in rows:
        base = {k: row[k] for k in ("gene", "tissue", "build", "kinship_scope", "ctnu_definition")}
        try:
            args, cells, paths = read_inputs(row, manifest.parent)
            for role, path in paths.items():
                provenance.append(dict(**base, role=role, path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            if not validate_only:
                np.random.seed(seed)
                est, p = fit.free_HE(**args, jk=True)
                shared = float(est["hom_g2"])
                v = np.diag(est["V"]).astype(float)
                specific = float(v.mean())
                admissible = shared >= 0 and bool((v >= 0).all()) and shared + specific > 0
                results.append(dict(**base, shared_variance=shared, mean_specific_variance=specific,
                    specificity=specific/(shared+specific) if admissible else np.nan,
                    variance_admissible=admissible, specific_p=float(p.get("V", np.nan)),
                    shared_p=float(p.get("hom_g2", np.nan)), N=args["Y"].shape[0], cell_types=len(cells),
                    method="CIGMA free_HE with donor jackknife", package_version=version,
                    evidence_type="expression genetic variance; not disease colocalization"))
                pvc = np.asarray(p.get("vc", np.full(len(cells), np.nan)))
                for i, cell in enumerate(cells):
                    cells_out.append(dict(**base, cell_type=cell, specific_variance=v[i], p=float(pvc[i])))
            status.append(dict(**base, status="validated" if validate_only else "ok", message=""))
        except Exception as e:
            status.append(dict(**base, status="failed", message=f"{type(e).__name__}: {e}"))
    result = pd.DataFrame(results)
    # Include all requested gene/tissue tests, not only successful fits.
    all_p = [next((r["specific_p"] for r in results if (r["gene"],r["tissue"]) == k), np.nan) for k in keys]
    q = dict(zip(keys, bh(all_p)))
    if len(result):
        result["specific_FDR_manifest"] = [q[(r["gene"],r["tissue"])] for r in results]
        result.to_csv(outdir / "cigma.results.csv", index=False)
        pd.DataFrame(cells_out).to_csv(outdir / "cigma.cell_types.csv", index=False)
    else:
        # Never leave a previous successful result alongside a failed rerun.
        pd.DataFrame(columns=["gene","tissue","specific_p","specific_FDR_manifest","specificity"]).to_csv(outdir / "cigma.results.csv", index=False)
        pd.DataFrame(columns=["gene","tissue","cell_type","specific_variance","p"]).to_csv(outdir / "cigma.cell_types.csv", index=False)
    pd.DataFrame(status).to_csv(outdir / "cigma.status.csv", index=False)
    pd.DataFrame(provenance).to_csv(outdir / "cigma.inputs.csv", index=False)
    (outdir / "cigma.provenance.json").write_text(json.dumps(dict(repository=REPOSITORY,
        adapter_verified_commit=VERIFIED_COMMIT, installed_version=version, manifest=str(manifest),
        tests=len(rows), seed=seed, validate_only=validate_only), indent=2))
    return 1 if any(s["status"] == "failed" for s in status) else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    raise SystemExit(run(a.manifest.resolve(), a.outdir.resolve(), a.validate_only, a.seed))
