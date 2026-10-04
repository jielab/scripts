#!/usr/bin/env python3


# 🚩 c5.cellulation
"""Cell-expression enrichment and native CIGMA, without fabricated cell labels.

Required resources are explicit. Expression enrichment is not secretion tracing,
cell abundance deconvolution, a cellular-aging clock, or disease colocalization.
Use the ``cigma`` subcommand for the CIGMA-HE adapter.
"""

from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import hashlib
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import csv
import importlib.metadata

HERE = Path(__file__).resolve().parent


REPOSITORY = "https://github.com/Minhui-Chen/CIGMA"
VERIFIED_COMMIT = "5813e4ae84d7b3733dfcd938fe42d12c6b30a8aa"


def matrix(path):
	# Preserve donor identifiers such as '001'; never align by row position.
	with Path(path).open(newline="") as stream:
		header = next(csv.reader(stream), [])
	# pandas otherwise silently renames duplicate column headers.
	if (
		len(header) < 2
		or len(set(header)) != len(header)
		or any(not h.strip() for h in header[1:])
	):
		raise ValueError(f"Missing/duplicate matrix column labels: {path}")
	d = pd.read_csv(path, dtype=str)
	if d.shape[1] < 2:
		raise ValueError(f"Expected labelled matrix: {path}")
	d = d.set_index(d.columns[0])
	if (
		d.index.has_duplicates
		or d.columns.has_duplicates
		or d.index.isna().any()
		or any(not str(v).strip() for v in d.index)
	):
		raise ValueError(f"Missing/duplicate identifiers: {path}")
	d = d.apply(pd.to_numeric, errors="raise")
	if not np.isfinite(d.to_numpy()).all():
		raise ValueError(
			f"Nonfinite values: {path}; missing donor/cell pairs need preprocessing"
		)
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
		raise ValueError(
			"ctnu_definition must be variance_of_pseudobulk_mean: "
			"upstream preprocess.pseudobulk returns SEM squared, not cell variance or SD"
		)
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
		raise ValueError(
			"Identity kinship cannot separate genetic and environmental covariance"
		)
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
		result[ix] = np.minimum(
			1,
			np.minimum.accumulate((p[ix] * len(p) / np.arange(1, len(ix) + 1))[::-1])[
				::-1
			],
		)
	return result


def run(manifest, outdir, validate_only=False, seed=2026):
	outdir.mkdir(parents=True, exist_ok=True)
	with manifest.open() as stream:
		rows = list(csv.DictReader(stream, delimiter="\t"))
	if (
		not rows
		or not {
			"gene",
			"ctp",
			"ctnu",
			"P",
			"K",
			"tissue",
			"build",
			"kinship_scope",
			"ctnu_definition",
		}
		<= rows[0].keys()
	):
		raise ValueError(
			"Manifest requires gene, ctp, ctnu, P, K, tissue, build, kinship_scope, ctnu_definition"
		)
	keys = [(r["gene"], r["tissue"]) for r in rows]
	if any(
		not r[k].strip()
		for r in rows
		for k in ("gene", "tissue", "build", "kinship_scope")
	):
		raise ValueError("Manifest provenance fields must not be blank")
	if len(set(keys)) != len(keys):
		raise ValueError("Duplicate gene/tissue jobs")
	version = "not loaded (validation only)"
	if not validate_only:
		from cigma import fit

		version = importlib.metadata.version("cigma")
	status, results, cells_out, provenance = [], [], [], []
	for row in rows:
		base = {
			k: row[k]
			for k in ("gene", "tissue", "build", "kinship_scope", "ctnu_definition")
		}
		try:
			args, cells, paths = read_inputs(row, manifest.parent)
			for role, path in paths.items():
				provenance.append(
					dict(
						**base,
						role=role,
						path=str(path),
						sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
					)
				)
			if not validate_only:
				np.random.seed(seed)
				est, p = fit.free_HE(**args, jk=True)
				shared = float(est["hom_g2"])
				v = np.diag(est["V"]).astype(float)
				specific = float(v.mean())
				admissible = (
					shared >= 0 and bool((v >= 0).all()) and shared + specific > 0
				)
				results.append(
					dict(
						**base,
						shared_variance=shared,
						mean_specific_variance=specific,
						specificity=(
							specific / (shared + specific) if admissible else np.nan
						),
						variance_admissible=admissible,
						specific_p=float(p.get("V", np.nan)),
						shared_p=float(p.get("hom_g2", np.nan)),
						N=args["Y"].shape[0],
						cell_types=len(cells),
						method="CIGMA free_HE with donor jackknife",
						package_version=version,
						evidence_type="expression genetic variance; not disease colocalization",
					)
				)
				pvc = np.asarray(p.get("vc", np.full(len(cells), np.nan)))
				for i, cell in enumerate(cells):
					cells_out.append(
						dict(
							**base,
							cell_type=cell,
							specific_variance=v[i],
							p=float(pvc[i]),
						)
					)
			status.append(
				dict(**base, status="validated" if validate_only else "ok", message="")
			)
		except Exception as e:
			status.append(
				dict(**base, status="failed", message=f"{type(e).__name__}: {e}")
			)
	result = pd.DataFrame(results)
	# Include all requested gene/tissue tests, not only successful fits.
	all_p = [
		next(
			(r["specific_p"] for r in results if (r["gene"], r["tissue"]) == k), np.nan
		)
		for k in keys
	]
	q = dict(zip(keys, bh(all_p)))
	if len(result):
		result["specific_FDR_manifest"] = [q[(r["gene"], r["tissue"])] for r in results]
		result.to_csv(outdir / "cigma.results.csv", index=False)
		pd.DataFrame(cells_out).to_csv(outdir / "cigma.cell_types.csv", index=False)
	else:
		# Never leave a previous successful result alongside a failed rerun.
		pd.DataFrame(
			columns=[
				"gene",
				"tissue",
				"specific_p",
				"specific_FDR_manifest",
				"specificity",
			]
		).to_csv(outdir / "cigma.results.csv", index=False)
		pd.DataFrame(
			columns=["gene", "tissue", "cell_type", "specific_variance", "p"]
		).to_csv(outdir / "cigma.cell_types.csv", index=False)
	pd.DataFrame(status).to_csv(outdir / "cigma.status.csv", index=False)
	pd.DataFrame(provenance).to_csv(outdir / "cigma.inputs.csv", index=False)
	(outdir / "cigma.provenance.json").write_text(
		json.dumps(
			dict(
				repository=REPOSITORY,
				adapter_verified_commit=VERIFIED_COMMIT,
				installed_version=version,
				manifest=str(manifest),
				tests=len(rows),
				seed=seed,
				validate_only=validate_only,
			),
			indent=2,
		)
	)
	return 1 if any(s["status"] == "failed" for s in status) else 0


def cigma_main(argv):
	ap = argparse.ArgumentParser(
		prog="c5.cellulation.py cigma",
		description="CIGMA-HE adapter for labelled donor/cell matrices",
	)
	ap.add_argument("--manifest", type=Path, required=True)
	ap.add_argument("--outdir", type=Path, required=True)
	ap.add_argument("--validate-only", action="store_true")
	ap.add_argument("--seed", type=int, default=2026)
	a = ap.parse_args(argv)
	raise SystemExit(
		run(a.manifest.resolve(), a.outdir.resolve(), a.validate_only, a.seed)
	)


def default_cell_atlas():
	"""Locate the shared CellAge annotation outside the script directory."""
	default = (
		"F:/annot/cellage/atlas.csv"
		if os.name == "nt"
		else "/mnt/f/annot/cellage/atlas.csv"
	)
	return Path(os.getenv("C5_CELL_ATLAS", default))


def import_cellage(repo, output):
	import shutil

	source = repo / "preprocessing/cell_type_mapping_update.csv"
	version = subprocess.check_output(
		["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
	).strip()
	table = pd.read_csv(source)
	rows = []
	for _, row in table.iterrows():
		for gene in str(row["Marker Genes"]).split(","):
			gene = gene.strip().upper()
			if gene and gene != "NAN":
				rows.append(
					dict(
						gene=gene,
						cell_type=row["Original Cell types"],
						source="Ding et al. CellAge author-published marker mapping; HPA-derived",
						source_version=version,
					)
				)
	out = pd.DataFrame(rows).drop_duplicates(["gene", "cell_type"])
	if out.empty or out.isna().any().any():
		raise ValueError("Empty or malformed CellAge mapping")
	output.mkdir(parents=True, exist_ok=True)
	out.to_csv(output / "atlas.csv", index=False)
	shutil.copy2(repo / "LICENSE", output / "LICENSE")
	provenance = dict(
		repository="https://github.com/dingdaisy/cellage",
		commit=version,
		source_file=str(source.relative_to(repo)),
		sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
		paper="https://doi.org/10.1038/s41591-026-04446-y",
		rule="Authors' supplied Marker Genes, not their SomaScan-only Somamer list",
		interpretation="Cell-enriched expression annotation, not secretion source or CIGMA inference",
		rows=len(out),
		genes=int(out.gene.nunique()),
		cell_types=int(out.cell_type.nunique()),
	)
	(output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
	print(json.dumps(provenance, indent=2))


def import_cellage_cli():
	parser = argparse.ArgumentParser(
		description="Import the authors' CellAge marker mapping"
	)
	parser.add_argument("--repo", type=Path, required=True)
	parser.add_argument("--outdir", type=Path, default=default_cell_atlas().parent)
	args = parser.parse_args(sys.argv[2:])
	import_cellage(args.repo, args.outdir)


def sha(path: Path) -> str:
	return hashlib.sha256(path.read_bytes()).hexdigest()


def annotate(universe, atlas, panels, outdir, prefix="cell"):
	u, a, p = (pd.read_csv(x, dtype=str) for x in (universe, atlas, panels))
	for frame, required in (
		(u, {"assay", "gene"}),
		(a, {"gene", "cell_type", "source", "source_version"}),
		(p, {"model", "feature"}),
	):
		if not required <= set(frame.columns):
			raise ValueError(f"Missing columns: {required - set(frame.columns)}")
	if u.assay.isna().any() or u.assay.duplicated().any():
		raise ValueError(
			"Universe must have unique nonmissing assays; leave unknown gene blank"
		)
	if a[["gene", "cell_type", "source", "source_version"]].isna().any().any():
		raise ValueError("Atlas labels and provenance must be nonmissing")
	if len(a[["source", "source_version"]].drop_duplicates()) != 1:
		raise ValueError(
			"Use one prespecified atlas/version per run; do not merge incompatible labels"
		)
	a = a.drop_duplicates(["gene", "cell_type"])
	p = p[p.feature.notna()].drop_duplicates(["model", "feature"])
	if not set(p.feature) <= set(u.assay):
		raise ValueError(
			"Every panel assay must be in the full measured assay universe"
		)
	background = set(u.gene.dropna()) - {""}
	if not background:
		raise ValueError("No mapped background genes")
	annotated = p.merge(u, left_on="feature", right_on="assay", how="left").merge(
		a, on="gene", how="left"
	)
	rows, coverage = [], []
	for model, group in p.groupby("model"):
		m = group.merge(u, left_on="feature", right_on="assay", how="left")
		genes = set(m.gene.dropna()) - {""}
		coverage.append(
			dict(
				model=model,
				assays=len(m),
				mapped_assays=int(m.gene.notna().sum()),
				unique_genes=len(genes),
				cell_labelled_genes=len(genes & set(a.gene)),
				background_genes=len(background),
			)
		)
		for cell, table in a.groupby("cell_type"):
			cellgenes = set(table.gene) & background
			hit = len(genes & cellgenes)
			tab = [
				[hit, len(genes) - hit],
				[len(cellgenes - genes), len(background - genes - cellgenes)],
			]
			odds, pv = (
				fisher_exact(tab, alternative="greater") if genes else (np.nan, np.nan)
			)
			rows.append(
				dict(
					model=model,
					cell_type=cell,
					hits=hit,
					panel_genes=len(genes),
					atlas_genes_in_assay_background=len(cellgenes),
					background_genes=len(background),
					odds_ratio=odds,
					p=pv,
				)
			)
	result = pd.DataFrame(rows)
	if len(result):
		result["FDR_all_panel_cell_tests"] = bh(result.p)
	outdir.mkdir(parents=True, exist_ok=True)
	annotated.to_csv(outdir / f"{prefix}.panel_annotation.csv", index=False)
	pd.DataFrame(coverage).to_csv(outdir / f"{prefix}.coverage.csv", index=False)
	result.to_csv(outdir / f"{prefix}.enrichment.csv", index=False)
	a[["source", "source_version"]].drop_duplicates().assign(
		interpretation="External gene-expression annotation; putative cellular relevance, not plasma protein origin",
		denominator="Unique mapped genes across all measured assays; NPPB and NTPROBNP counted once for enrichment",
	).to_csv(outdir / f"{prefix}.provenance.csv", index=False)
	if len(result):
		import matplotlib

		matplotlib.use("Agg")
		import matplotlib.pyplot as plt

		top = result.groupby("cell_type").hits.sum().nlargest(16).index
		z = result[result.cell_type.isin(top)].copy()
		z["evidence"] = -np.log10(z.FDR_all_panel_cell_tests.clip(lower=1e-300))
		# Average descriptive evidence across folds; never call this a combined P value.
		z["display_model"] = z.model.str.replace(r"\|[0-9]+$", "", regex=True)
		grid = z.pivot_table(
			index="cell_type",
			columns="display_model",
			values="evidence",
			aggfunc="mean",
		)
		fig, ax = plt.subplots(figsize=(max(9, min(28, len(grid.columns) * 0.28)), 7))
		im = ax.imshow(grid.fillna(0), aspect="auto", cmap="Blues")
		ax.set_yticks(range(len(grid.index)), grid.index)
		ax.set_xticks(range(len(grid.columns)), grid.columns, rotation=90, fontsize=6)
		ax.set_title(
			"External cell-expression annotation; descriptive mean across folds"
		)
		fig.colorbar(im, ax=ax, label="Mean -log10 adjusted P (not a combined test)")
		fig.tight_layout()
		fig.savefig(outdir / f"{prefix}.enrichment.png", dpi=180)
		plt.close(fig)


def annotation_cli():
	ap = argparse.ArgumentParser(description="Annotate frozen assay panels")
	for name in ["universe", "atlas", "panels", "outdir"]:
		ap.add_argument("--" + name, type=Path, required=True)
	ap.add_argument("--prefix", default="c5.cell")
	arg = ap.parse_args(sys.argv[2:])
	annotate(arg.universe, arg.atlas, arg.panels, arg.outdir, arg.prefix)


def default_resources(root, atlas):
	# Reuse explicit Final inputs when present. A clean default run reaches C5
	# before Final, so it must also work directly from the saved native tables.
	for directory, universe, panels in [
		(
			root.parent / "final" / root.name / "systematic",
			"c5.systematic.cell_universe.csv",
			"c5.systematic.cell_panels.csv",
		),
		(
			root.parent / "final" / root.name / "joint",
			"c5.cell_assay_universe.csv",
			"c5.cell_panels.csv",
		),
	]:
		if (directory / universe).is_file() and (directory / panels).is_file():
			return [atlas, directory / universe, directory / panels]
	from final import cell_annotation_inputs

	project = os.getenv("LE8_PUBLISHED_ROOT", str(root.parent))
	scratch = (
		Path("/tmp/le8-cache")
		/ hashlib.sha256(project.encode()).hexdigest()[:12]
		/ root.name
		/ "prot"
		/ "c5_inputs"
	)
	paths = cell_annotation_inputs(root, scratch)
	if paths is not None:
		return [atlas, *paths]
	return [atlas, None, None]


def main() -> None:
	ap = argparse.ArgumentParser(description=__doc__)
	ap.add_argument("--Y", default="cvd_cad,ra")
	ap.add_argument("--biom", default="prot,met")
	ap.add_argument("--analysis-root", type=Path, required=True)
	for name in ["atlas", "universe", "panels", "cigma-manifest", "cigma-results"]:
		ap.add_argument("--" + name, type=Path)
	ap.add_argument("--replace", action="store_true")
	ap.add_argument("--python", default=os.getenv("C5_CIGMA_PYTHON", sys.executable))
	arg = ap.parse_args()
	if arg.cigma_manifest and arg.cigma_results:
		ap.error("Use a native CIGMA manifest OR precomputed native results, not both.")
	resources = [arg.atlas, arg.universe, arg.panels]
	if (arg.universe is None) != (arg.panels is None):
		ap.error(
			"Provide --universe and --panels together; --atlas defaults to the shared CellAge annotation (C5_CELL_ATLAS)."
		)
	if arg.atlas is None:
		arg.atlas = default_cell_atlas()
	resources = [arg.atlas, arg.universe, arg.panels]
	if arg.panels and ("," in arg.Y or arg.biom != "prot"):
		ap.error(
			"A declared panel file must be annotated for exactly one --Y and --biom prot."
		)
	if (arg.cigma_manifest or arg.cigma_results) and "," in arg.Y:
		ap.error(
			"CIGMA resources must be attached to one explicitly declared outcome per invocation."
		)
	for path in resources + [arg.cigma_manifest, arg.cigma_results]:
		if path is not None and not path.is_file():
			raise FileNotFoundError(path)
	for trait in arg.Y.split(","):
		for layer in arg.biom.split(","):
			if (
				layer not in {"prot", "met"}
				or not trait
				or "/" in trait
				or ".." in trait
			):
				raise ValueError("Invalid trait or layer")
			out = arg.analysis_root / trait / layer / "c5_cellulation"
			out.mkdir(parents=True, exist_ok=True)
			status = []
			inputs = []
			resources = [arg.atlas, arg.universe, arg.panels]
			if layer == "prot" and arg.universe is None:
				root = arg.analysis_root / trait
				resources = default_resources(root, arg.atlas)
			cache_file = out / "c5.completed.json"
			signature = dict(
				trait=trait,
				layer=layer,
				code_sha256=sha(Path(__file__)),
				inputs=[
					(
						dict(path=str(p.resolve()), sha256=sha(p))
						if p is not None
						else None
					)
					for p in resources + [arg.cigma_manifest, arg.cigma_results]
				],
			)
			if not arg.replace and not arg.cigma_manifest and cache_file.is_file():
				saved = json.loads(cache_file.read_text())
				outputs = saved.get("outputs", {})
				if (
					saved.get("signature") == signature
					and outputs
					and all(
						(out / n).is_file() and sha(out / n) == digest
						for n, digest in outputs.items()
					)
				):
					print(
						json.dumps(
							dict(
								trait=trait,
								layer=layer,
								mode="reuse-completed",
								status=saved["status"],
							)
						),
						flush=True,
					)
					continue
			# Replace disposable driver outputs; reusable RDS data stay in place.
			for item in out.glob("c5.*"):
				if item.is_file() and item.suffix != ".rds":
					item.unlink()
			pd.DataFrame(
				[
					dict(
						analysis="cell_expression",
						status="running",
						detail="Not yet completed",
					),
					dict(
						analysis="CIGMA", status="running", detail="Not yet completed"
					),
				]
			).to_csv(out / "c5.cellulation_status.csv", index=False)
			if layer == "met":
				status.append(
					dict(
						analysis="cell_expression",
						status="not_applicable",
						detail="A metabolite has no unique encoding gene. No direct protein-style cell mapping is performed.",
					)
				)
			elif all(resources):
				annotate(resources[1], resources[0], resources[2], out, "c5.cell")
				status.append(
					dict(
						analysis="cell_expression",
						status="completed",
						detail="Full measured-assay background; one declared atlas/version; BH across all panel-cell tests.",
					)
				)
				for label, path in zip(["atlas", "universe", "panels"], resources):
					inputs.append(
						dict(role=label, path=str(path.resolve()), sha256=sha(path))
					)
			else:
				status.append(
					dict(
						analysis="cell_expression",
						status="unavailable",
						detail="Provide an external cell-expression atlas, full assay-gene universe and frozen panel file. Absence is not a null result.",
					)
				)
			results = arg.cigma_results
			if arg.cigma_manifest and layer == "prot":
				native_out = out / "cigma_native"
				native_out.mkdir(exist_ok=True)
				native = HERE / "c5.cellulation.py"
				with (native_out / "runner.log").open("w") as log:
					process = subprocess.run(
						[
							arg.python,
							str(native),
							"cigma",
							"--manifest",
							str(arg.cigma_manifest),
							"--outdir",
							str(native_out),
						],
						stdout=log,
						stderr=subprocess.STDOUT,
					)
				inputs.append(
					dict(
						role="cigma_manifest",
						path=str(arg.cigma_manifest.resolve()),
						sha256=sha(arg.cigma_manifest),
					)
				)
				if process.returncode:
					status.append(
						dict(
							analysis="CIGMA",
							status="failed",
							detail="Native runner failed; partial results not promoted. See cigma_native/runner.log.",
						)
					)
					results = None
				else:
					results = native_out / "cigma.results.csv"
			if results and layer == "prot":
				d = pd.read_csv(results)
				required = {
					"gene",
					"tissue",
					"specific_p",
					"specific_FDR_manifest",
					"specificity",
				}
				if not required <= set(d):
					raise ValueError(
						"Native CIGMA results missing "
						+ ", ".join(sorted(required - set(d)))
					)
				if (
					d[["gene", "tissue"]].isna().any().any()
					or d.duplicated(["gene", "tissue"]).any()
				):
					raise ValueError(
						"CIGMA gene/tissue identifiers must be unique and nonmissing."
					)
				d.assign(
					interpretation="Cell-specific expression regulation; no disease causal direction, plasma origin or PP(H4) is inferred"
				).to_csv(out / "c5.CIGMA_results.csv", index=False)
				inputs.append(
					dict(
						role="CIGMA_results",
						path=str(results.resolve()),
						sha256=sha(results),
					)
				)
				status.append(
					dict(
						analysis="CIGMA",
						status="completed",
						detail=f"{len(d)} native gene/tissue rows; FDR family retained, not recalculated on selected proteins.",
					)
				)
				if resources[1]:
					u = pd.read_csv(resources[1], dtype=str)
					if not {"assay", "gene"} <= set(u):
						raise ValueError("Universe requires assay,gene")
					# All measured mapped assays, not just the coloc-positive subset.
					annotation = u.dropna(subset=["gene"]).merge(
						d, on="gene", how="inner"
					)
					annotation.to_csv(out / "c5.CIGMA_annotation.csv", index=False)
			elif not any(x["analysis"] == "CIGMA" for x in status):
				status.append(
					dict(
						analysis="CIGMA",
						status="unavailable" if layer == "prot" else "not_applicable",
						detail="Needs donor-level pseudobulk expression/noise, cell proportions and kinship or native CIGMA output; UKB plasma proteomics/PGS alone are insufficient.",
					)
				)
			pd.DataFrame(status).to_csv(out / "c5.cellulation_status.csv", index=False)
			pd.DataFrame(inputs, columns=["role", "path", "sha256"]).to_csv(
				out / "c5.cellulation_inputs.csv", index=False
			)
			if not any(x["status"] == "failed" for x in status):
				saved = dict(
					signature=signature,
					status=status,
					outputs={
						p.name: sha(p)
						for p in out.glob("c5.*")
						if p.is_file() and p != cache_file
					},
				)
				temporary = cache_file.with_suffix(".tmp")
				temporary.write_text(json.dumps(saved, indent=2) + "\n")
				temporary.replace(cache_file)
			print(
				json.dumps(
					dict(trait=trait, layer=layer, status=status), ensure_ascii=False
				),
				flush=True,
			)


if __name__ == "__main__":
	try:
		if sys.argv[1:2] == ["cigma"]:
			cigma_main(sys.argv[2:])
		elif sys.argv[1:2] == ["--import-cellage"]:
			import_cellage_cli()
		elif sys.argv[1:2] == ["--annotate"]:
			annotation_cli()
		else:
			main()
	except Exception as exc:
		print("ERROR: " + str(exc), file=sys.stderr)
		sys.exit(2)
