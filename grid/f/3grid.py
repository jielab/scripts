#!/usr/bin/env python3
"""3grid.py: inputs, ld, transport, fit, weights, mix. Use a subcommand followed by --help."""

from __future__ import annotations
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

if "grid_0_common" not in sys.modules:
	_spec = spec_from_file_location("grid_0_common", Path(__file__).with_name("0.common.py"))
	_common = module_from_spec(_spec)
	sys.modules[_spec.name] = _common
	try:
		_spec.loader.exec_module(_common)
	except BaseException:
		sys.modules.pop(_spec.name, None)
		raise
from grid_0_common import load_module


# 🚩 inputs: grid_inputs
"""Bridge permanent 1csx outputs to the GRID transport/weight input contract."""
import argparse, json
from pathlib import Path
import pandas as pd

load_module("1csx.py")
from grid_1csx import find_gwas, POPS
from grid_0_common import preparation_signature
from grid_0_common import signature
from grid_0_common import ensure_prepared


def grid_inputs_main():
	p = argparse.ArgumentParser()
	for k in ("trait", "gwas-dir", "snpinfo", "chrs", "out"):
		p.add_argument("--" + k, required=True)
	a = p.parse_args()
	out = Path(a.out)
	chrs = list(map(int, a.chrs.split()))
	suffix = "" if chrs == list(range(1, 23)) else ".chr" + ",".join(map(str, chrs))
	srcs = [find_gwas(a.gwas_dir, a.trait, pop) for pop in POPS]
	weights = [Path(str(src)[:-3] + suffix + ".csx.gz") for src in srcs]
	signatures = []
	metas = []
	for src, w, pop in zip(srcs, weights, POPS):
		sf = Path(str(w) + ".signature")
		mf = Path(str(w) + ".metadata.json")
		if not all(f.is_file() for f in (w, sf, mf)):
			raise ValueError(
				f"Run 1csx.sh --stage weights with the same chromosomes first: missing {w} or its sidecars"
			)
		meta = json.loads(mf.read_text())
		if meta.get("preparation_signature") != preparation_signature(src, a.snpinfo, a.trait, pop):
			raise ValueError(f"CSx source/preparation mismatch: {w}; rerun 1csx.sh")
		signatures.append(sf.read_text().strip())
		metas.append(meta)
	if len(set(signatures)) != 1:
		raise ValueError("Population weights come from different joint CSx runs")
	files = (
		srcs
		+ weights
		+ [Path(str(w) + ext) for w in weights for ext in (".signature", ".metadata.json")]
		+ [Path(a.snpinfo), Path(__file__)]
	)
	key = signature([a.trait, chrs], files)
	cache = out / "grid" / "inputs" / key
	cache.mkdir(parents=True, exist_ok=True)
	dest = out / "sumstats" / "bychr"
	dest.mkdir(parents=True, exist_ok=True)
	wd = out / "csx" / "weights"
	wd.mkdir(parents=True, exist_ok=True)
	manifest = []
	marker = out / "grid" / "inputs.signature"
	expected = (
		[dest / f"{a.trait}.{pop}.chr{c}.tsv.gz" for pop in POPS for c in chrs]
		+ [wd / f"{pop}.chr{c}.tsv" for pop in POPS for c in chrs]
		+ [out / "csx" / "manifest.tsv"]
	)
	if (
		marker.is_file()
		and marker.read_text().strip() == key
		and all(f.is_file() and f.stat().st_size for f in expected)
	):
		print("SKIP GRID inputs: matching permanent CSx/GWAS inputs")
		return
	for pop, src, w, meta in zip(POPS, srcs, weights, metas):
		table = cache / f"{pop}.tsv.gz"
		info = cache / f"{pop}.json"
		ensure_prepared(src, a.snpinfo, a.trait, pop, output=table, metadata=info)
		wanted = {c: [] for c in chrs}
		for chunk in pd.read_csv(table, sep="\t", chunksize=500000):
			for c in chrs:
				z = chunk.loc[chunk.CHR == c]
				if not z.empty:
					wanted[c].append(z)
		beta = pd.read_csv(w, sep="\t")
		for c in chrs:
			if not wanted[c] or not (beta.CHR == c).any():
				raise ValueError(f"No GWAS/CSx variants for {pop} chr{c}")
			pd.concat(wanted[c]).to_csv(
				dest / f"{a.trait}.{pop}.chr{c}.tsv.gz", sep="\t", index=False, compression="gzip"
			)
			beta.loc[beta.CHR == c].to_csv(wd / f"{pop}.chr{c}.tsv", sep="\t", index=False)
		n = meta.get("n_gwas_used", meta["n_gwas_median"])
		if n is None or float(n) <= 0:
			raise ValueError(f"Missing actual GWAS N for {pop}")
		manifest.append(
			{"pop": pop, "n_gwas": n, "source": str(src), "weights": str(w), "joint_signature": signatures[0]}
		)
	pd.DataFrame(manifest).to_csv(out / "csx" / "manifest.tsv", sep="\t", index=False)
	(out / "grid" / "inputs.signature").write_text(key + "\n")
	print(f"GRID inputs ready: {len(chrs)} chromosomes; joint CSx={signatures[0][:12]}", flush=True)


def grid_inputs_cli():
	grid_inputs_main()


# 🚩 ld: extract_ld_scores
"""Extract per-SNP LD scores (sum r^2, diagonal included) from PRS-CSx HDF5 blocks."""
import argparse, gzip
from pathlib import Path
import h5py, numpy as np, pandas as pd


def dec(x):
	if isinstance(x, (bytes, np.bytes_)):
		return x.decode()
	if isinstance(x, np.void) and x.dtype.names:
		for n in x.dtype.names:
			if "snp" in n.lower() or "rs" in n.lower() or "id" == n.lower():
				return dec(x[n])
	return str(x)


def groups(h):
	out = []

	def visit(name, obj):
		if isinstance(obj, h5py.Group):
			ds = {k: v for k, v in obj.items() if isinstance(v, h5py.Dataset)}
			mats = [(k, v) for k, v in ds.items() if v.ndim == 2 and v.shape[0] == v.shape[1]]
			vec = [(k, v) for k, v in ds.items() if v.ndim == 1]
			for mk, m in mats:
				sv = next(
					(
						(k, v)
						for k, v in vec
						if len(v) == m.shape[0] and any(z in k.lower() for z in ["snp", "rs", "id"])
					),
					None,
				)
				if sv is None:
					sv = next(((k, v) for k, v in vec if len(v) == m.shape[0]), None)
				if sv:
					out.append((name, mk, m, sv[0], sv[1]))
					break

	h.visititems(visit)
	return out


def extract_ld_scores_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--hdf5", required=True)
	ap.add_argument("--out", required=True)
	ap.add_argument("--pop", required=True)
	ap.add_argument("--chr", required=True)
	a = ap.parse_args()
	rows = []
	with h5py.File(a.hdf5, "r") as h:
		gs = groups(h)
		if not gs:
			raise SystemExit(f"No square LD matrix + SNP list found in {a.hdf5}")
		for name, mk, m, sk, s in gs:
			x = np.asarray(m[...], dtype=np.float64)
			ld = np.einsum("ij,ij->i", x, x)
			ids = [dec(z) for z in s[...]]
			rows.extend(zip(ids, ld))
	d = pd.DataFrame(rows, columns=["SNP", "ldscore"]).drop_duplicates("SNP")
	d["pop"] = a.pop.upper()
	d["chr"] = str(a.chr)
	Path(a.out).parent.mkdir(parents=True, exist_ok=True)
	d.to_csv(a.out, sep="\t", index=False, compression="gzip")
	print(f"blocks={len(gs)} snps={len(d)}")


def extract_ld_scores_cli():
	extract_ld_scores_main()


# 🚩 transport: build_transport_table
"""Build pairwise cross-ancestry effect-transportability observations."""
import argparse, gzip, itertools, math, re
from pathlib import Path
import numpy as np, pandas as pd

POPS0 = ["AFR", "EAS", "EUR", "SAS"]
COMP = str.maketrans("ACGT", "TGCA")


def transport_complement(x):
	return str(x).translate(COMP)


def read_transport_table(path):
	return pd.read_csv(path, sep="\t", compression="infer", low_memory=False, dtype={"SNP": str, "CHR": str})


def centers(path, pops, npc=10):
	d = pd.read_csv(path, sep=None, engine="python", compression="infer")
	pc = sorted(
		[c for c in d.columns if re.fullmatch(r"PC\d+", str(c), re.I)], key=lambda x: int(re.findall(r"\d+", x)[0])
	)[:npc]
	pcoll = next((c for c in ["pop", "POP", "super_pop", "ancestry", "population", "race"] if c in d.columns), None)
	if pcoll is None or len(pc) < 2:
		raise SystemExit(f"Cannot read population PC centers from {path}: {list(d.columns)}")
	d[pcoll] = d[pcoll].astype(str).str.upper()
	out = {}
	for p in pops:
		x = d[d[pcoll] == p]
		if len(x):
			out[p] = x[pc].apply(pd.to_numeric, errors="coerce").median().to_numpy(float)
	if len(out) < 2:
		raise SystemExit(f"Fewer than two population centers in {path}")
	return out


def align_transport(d, ref):
	x = d.merge(ref, on="SNP", how="inner", suffixes=("", "_REF"))
	a1 = x.A1.astype(str).str.upper()
	a2 = x.A2.astype(str).str.upper()
	r1 = x.A1_REF.astype(str).str.upper()
	r2 = x.A2_REF.astype(str).str.upper()
	same = (a1 == r1) & (a2 == r2)
	swap = (a1 == r2) & (a2 == r1)
	cs = (a1.map(transport_complement) == r1) & (a2.map(transport_complement) == r2)
	cw = (a1.map(transport_complement) == r2) & (a2.map(transport_complement) == r1)
	ok = (same | swap | cs | cw) & (pd.to_numeric(x.BP, errors="coerce") == pd.to_numeric(x.BP_REF, errors="coerce"))
	x = x[ok].copy()
	flip = (swap | cw)[ok].to_numpy()
	x["BETA_ALIGNED"] = pd.to_numeric(x.BETA, errors="coerce").to_numpy() * np.where(flip, -1, 1)
	e = pd.to_numeric(x.get("EAF"), errors="coerce")
	x["EAF_ALIGNED"] = np.where(flip, 1 - e, e)
	return x


def external_age(path):
	if not path:
		return None
	d = pd.read_csv(path, sep=None, engine="python", compression="infer")
	sc = next((c for c in ["SNP", "rsid", "RSID", "id", "variant_id"] if c in d.columns), None)
	ac = next(
		(
			c
			for c in d.columns
			if any(z in str(c).lower() for z in ["age_gen", "allele_age", "geva_age", "mean_age", "age"])
		),
		None,
	)
	if sc is None or ac is None:
		raise SystemExit(f"Cannot find SNP/age in {path}")
	age = pd.to_numeric(d[ac], errors="coerce")
	unit = "generations"
	if "year" in ac.lower():
		age = age / 29.0
		unit = "years_to_generations"
	return pd.DataFrame({"SNP": d[sc].astype(str), "external_age_gen": age}).drop_duplicates("SNP"), unit


def build_transport_table_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--trait", required=True)
	ap.add_argument("--pops", default="AFR,EAS,EUR,SAS")
	ap.add_argument("--chrs", required=True)
	ap.add_argument("--sumstats-dir", required=True)
	ap.add_argument("--arg-dir", required=True)
	ap.add_argument("--ldscore-dir", required=True)
	ap.add_argument("--centers", required=True)
	ap.add_argument("--external-age", default="")
	ap.add_argument("--max-snps-per-chr", type=int, default=0)
	ap.add_argument("--out", required=True)
	a = ap.parse_args()
	pops = [x.upper() for x in re.split("[,; ]+", a.pops) if x]
	chrs = [x for x in re.split("[,; ]+", a.chrs) if x]
	cen = centers(a.centers, pops)
	ext = external_age(a.external_age)
	extdf = ext[0] if ext else None
	Path(a.out).parent.mkdir(parents=True, exist_ok=True)
	first = True
	total = 0
	for chrom in chrs:
		ds = {}
		for p in pops:
			f = Path(a.sumstats_dir) / f"{a.trait}.{p}.chr{chrom}.tsv.gz"
			if not f.exists():
				continue
			d = read_transport_table(f)
			d["SNP"] = d.SNP.astype(str)
			d["A1"] = d.A1.astype(str).str.upper()
			d["A2"] = d.A2.astype(str).str.upper()
			d["CHR"] = str(chrom)
			d["BP"] = pd.to_numeric(d.BP, errors="coerce")
			d = d.dropna(subset=["SNP", "A1", "A2", "BETA", "SE"]).drop_duplicates("SNP")
			ds[p] = d
		if len(ds) < 2:
			continue
		priority = [p for p in ["EUR", "AFR", "EAS", "SAS"] if p in ds]
		long = pd.concat(
			[ds[p][["SNP", "A1", "A2", "BP"]].assign(_prio=i) for i, p in enumerate(priority)], ignore_index=True
		).sort_values("_prio")
		ref = long.drop_duplicates("SNP")[["SNP", "A1", "A2", "BP"]]
		wide = ref.rename(columns={"A1": "REF_A1", "A2": "REF_A2", "BP": "REF_BP"}).copy()
		refalign = ref.rename(columns={"A1": "A1_REF", "A2": "A2_REF", "BP": "BP_REF"})
		for p, d in ds.items():
			z = align_transport(d, refalign)
			z = z[["SNP", "BETA_ALIGNED", "SE", "EAF_ALIGNED", "P", "N", "BP"]].rename(
				columns={c: f"{c}_{p}" for c in ["BETA_ALIGNED", "SE", "EAF_ALIGNED", "P", "N", "BP"]}
			)
			wide = wide.merge(z, on="SNP", how="left")
		if a.max_snps_per_chr > 0:
			wide = wide.sort_values("REF_BP").head(a.max_snps_per_chr)
		# Join ARG annotations by exact rsID and coordinate; never borrow a nearby site.
		afile = Path(a.arg_dir) / f"chr{chrom}.variants.tsv.gz"
		arg = read_transport_table(afile) if afile.exists() else pd.DataFrame()
		if len(arg):
			arg["SNP"] = arg.SNP.astype(str)
			arg["bp"] = pd.to_numeric(arg.bp, errors="coerce")
			arg = arg[~arg.SNP.duplicated(keep=False)]
			keep = [c for c in arg.columns if c not in {"chr", "allele0", "allele1"}]
			wide = wide.merge(arg[keep], on="SNP", how="left")
			if "bp" in wide:
				bad = wide.bp.notna() & (pd.to_numeric(wide.bp, errors="coerce") != wide.REF_BP)
				for field in [c for c in arg.columns if c not in {"chr", "SNP", "bp", "allele0", "allele1"}]:
					if field in wide:
						wide.loc[bad, field] = np.nan
			miss = wide["arg_age_gen"].isna() if "arg_age_gen" in wide else pd.Series(True, index=wide.index)
			if miss.any():
				apos = arg[~arg.bp.duplicated(keep=False)].set_index("bp")
				ix = pd.to_numeric(wide.loc[miss, "REF_BP"], errors="coerce")
				for c in [x for x in arg.columns if x not in {"chr", "SNP", "bp", "allele0", "allele1"}]:
					if c not in wide:
						wide[c] = np.nan
					wide.loc[miss, c] = ix.map(apos[c])
		if extdf is not None:
			wide = wide.merge(extdf, on="SNP", how="left")
		# LD score joins.
		for p in pops:
			lf = Path(a.ldscore_dir) / p / f"chr{chrom}.ldscore.tsv.gz"
			if lf.exists():
				ld = (
					read_transport_table(lf)[["SNP", "ldscore"]]
					.drop_duplicates("SNP")
					.rename(columns={"ldscore": f"ldscore_{p}"})
				)
				wide = wide.merge(ld, on="SNP", how="left")
		rows = []
		for pa, pb in itertools.combinations(pops, 2):
			ba = f"BETA_ALIGNED_{pa}"
			bb = f"BETA_ALIGNED_{pb}"
			sa = f"SE_{pa}"
			sb = f"SE_{pb}"
			if not all(c in wide for c in [ba, bb, sa, sb]):
				continue
			x = wide[wide[ba].notna() & wide[bb].notna() & wide[sa].notna() & wide[sb].notna()].copy()
			if not len(x):
				continue
			va = pd.to_numeric(x[sa], errors="coerce") ** 2
			vb = pd.to_numeric(x[sb], errors="coerce") ** 2
			svar = va + vb
			q = (pd.to_numeric(x[ba]) - pd.to_numeric(x[bb])) ** 2 / svar.replace(0, np.nan)
			y = np.log1p(np.maximum(q - 1, 0))
			ea = pd.to_numeric(x.get(f"EAF_ALIGNED_{pa}"), errors="coerce")
			eb = pd.to_numeric(x.get(f"EAF_ALIGNED_{pb}"), errors="coerce")
			# ARG frequencies are orientation-free after conversion to MAF.
			if f"af_{pa}" in x:
				ea = ea.fillna(pd.to_numeric(x[f"af_{pa}"], errors="coerce"))
			if f"af_{pb}" in x:
				eb = eb.fillna(pd.to_numeric(x[f"af_{pb}"], errors="coerce"))
			ma = np.minimum(ea, 1 - ea)
			mb = np.minimum(eb, 1 - eb)
			la = pd.to_numeric(x.get(f"ldscore_{pa}"), errors="coerce")
			lb = pd.to_numeric(x.get(f"ldscore_{pb}"), errors="coerce")
			dc = float(np.linalg.norm(cen[pa] - cen[pb])) if pa in cen and pb in cen else np.nan
			dcol = f"div_{pa}_{pb}" if f"div_{pa}_{pb}" in x else f"div_{pb}_{pa}"
			local = pd.to_numeric(x.get(dcol), errors="coerce")
			age = pd.to_numeric(x.get("arg_age_gen"), errors="coerce")
			if "external_age_gen" in x:
				age = age.fillna(pd.to_numeric(x.external_age_gen, errors="coerce"))
			z = pd.DataFrame(
				{
					"trait": a.trait,
					"chr": str(chrom),
					"bp": x.REF_BP,
					"SNP": x.SNP,
					"A1": x.REF_A1,
					"A2": x.REF_A2,
					"pop_a": pa,
					"pop_b": pb,
					"beta_a": x[ba],
					"beta_b": x[bb],
					"se_a": x[sa],
					"se_b": x[sb],
					"transport_heterogeneity": y,
					"sampling_var": svar,
					"log_sampling_var": np.log(svar),
					"maf_a": ma,
					"maf_b": mb,
					"mean_maf": (ma + mb) / 2,
					"delta_maf": abs(ma - mb),
					"ldscore_a": la,
					"ldscore_b": lb,
					"mean_ldscore": (la + lb) / 2,
					"delta_ldscore": abs(la - lb),
					"global_pca_distance": dc,
					"local_arg_divergence": local,
					"local_global_ratio": local / (dc + 1e-8),
					"arg_age_gen": age,
					"log_arg_age": np.log1p(age),
				}
			)
			for c in [
				"mutation_branch_gen",
				"root_time_gen",
				"carrier_frequency",
				"lineage_breadth",
				"lineage_entropy",
			]:
				z[c] = pd.to_numeric(x.get(c), errors="coerce")
			rows.append(z)
		if rows:
			z = pd.concat(rows, ignore_index=True)
			z.to_csv(
				a.out,
				sep="\t",
				index=False,
				compression="gzip",
				mode="wt" if first else "at",
				header=first,
				float_format="%.9g",
			)
			first = False
			total += len(z)
	if first:
		raise SystemExit("No pairwise transport observations were generated")
	print(f"rows={total} output={a.out}")


def build_transport_table_cli():
	build_transport_table_main()


# 🚩 fit: fit_transport_model
"""Fit baseline vs evolutionary transportability models with blocked out-of-fold prediction."""
import argparse, json, math
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

BASE = ["global_pca_distance", "mean_maf", "delta_maf", "mean_ldscore", "delta_ldscore", "log_sampling_var"]
EVOL = [
	"log_arg_age",
	"mutation_branch_gen",
	"root_time_gen",
	"carrier_frequency",
	"lineage_breadth",
	"lineage_entropy",
	"local_arg_divergence",
	"local_global_ratio",
]


def prepare(train, valid, features):
	med = train[features].replace([np.inf, -np.inf], np.nan).median(numeric_only=True).fillna(0)
	tr = train[features].replace([np.inf, -np.inf], np.nan).fillna(med).to_numpy(float)
	va = valid[features].replace([np.inf, -np.inf], np.nan).fillna(med).to_numpy(float)
	mu = np.nanmean(tr, 0)
	sd = np.nanstd(tr, 0)
	sd[~np.isfinite(sd) | (sd < 1e-12)] = 1
	return (tr - mu) / sd, (va - mu) / sd, med.to_dict(), mu, sd


def ridge(X, y, w, alpha):
	X = np.column_stack([np.ones(len(X)), X])
	sw = np.sqrt(w)
	Xw = X * sw[:, None]
	yw = y * sw
	pen = np.eye(X.shape[1]) * alpha
	pen[0, 0] = 0
	try:
		return np.linalg.solve(Xw.T @ Xw + pen, Xw.T @ yw)
	except np.linalg.LinAlgError:
		return np.linalg.pinv(Xw.T @ Xw + pen) @ (Xw.T @ yw)


def pred(X, b):
	return b[0] + X @ b[1:]


def metrics(y, p):
	ok = np.isfinite(y) & np.isfinite(p)
	y = y[ok]
	p = p[ok]
	rmse = float(np.sqrt(np.mean((y - p) ** 2)))
	mae = float(np.mean(abs(y - p)))
	r = float(spearmanr(y, p).statistic) if len(y) > 2 else np.nan
	den = np.sum((y - y.mean()) ** 2)
	r2 = float(1 - np.sum((y - p) ** 2) / den) if den > 0 else np.nan
	return {"n": int(len(y)), "rmse": rmse, "mae": mae, "spearman": r, "r2": r2}


def fit_oof(d, features, groups, alpha):
	p = np.full(len(d), np.nan)
	rec = []
	for g in sorted(pd.unique(groups)):
		te = np.asarray(groups == g)
		tr = ~te
		if tr.sum() < max(50, 3 * len(features)) or te.sum() < 10:
			continue
		Xtr, Xte, _, _, _ = prepare(d.loc[tr], d.loc[te], features)
		y = d["transport_heterogeneity"].to_numpy(float)
		w = d["model_weight"].to_numpy(float)
		b = ridge(Xtr, np.minimum(y[tr], np.quantile(y[tr], 0.995)), training_weights(d.loc[tr]), alpha)
		p[te] = pred(Xte, b)
		rec.append({"validation_group": str(g), **metrics(y[te], p[te])})
	return p, rec


def full_fit(d, features, alpha):
	X, _, med, mu, sd = prepare(d, d, features)
	y = d.transport_heterogeneity.to_numpy(float)
	w = training_weights(d)
	b = ridge(X, np.minimum(y, np.quantile(y, 0.995)), w, alpha)
	return b, med, mu, sd


def training_weights(d):
	sv = pd.to_numeric(d.sampling_var, errors="coerce").to_numpy(float)
	if np.any(~np.isfinite(sv) | (sv <= 0)):
		raise ValueError("sampling_var must be finite and positive")
	w = 1 / np.sqrt(np.maximum(sv, np.quantile(sv, 0.01)))
	w = np.clip(w, *np.quantile(w, [0.01, 0.99]))
	return w / w.mean()


def fit_transport_model_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--input", required=True)
	ap.add_argument("--out-dir", required=True)
	ap.add_argument("--ridge-alpha", type=float, default=10)
	ap.add_argument("--conservation-min", type=float, default=0.05)
	ap.add_argument("--conservation-max", type=float, default=0.995)
	ap.add_argument("--seed", type=int, default=20260904)
	ap.add_argument("--model", choices=["baseline", "evolutionary_full"], default="evolutionary_full")
	a = ap.parse_args()
	out = Path(a.out_dir)
	out.mkdir(parents=True, exist_ok=True)
	use = (
		["trait", "chr", "bp", "SNP", "A1", "A2", "pop_a", "pop_b", "transport_heterogeneity", "sampling_var"]
		+ BASE
		+ EVOL
	)
	d = pd.read_csv(
		a.input,
		sep="\t",
		compression="infer",
		usecols=lambda x: x in use,
		low_memory=False,
		dtype={"chr": str, "SNP": str},
	)
	d["transport_heterogeneity"] = pd.to_numeric(d.transport_heterogeneity, errors="coerce")
	d = d[np.isfinite(d.transport_heterogeneity)].reset_index(drop=True)
	if len(d) < 200:
		raise SystemExit(f"Too few transport observations: {len(d)}")
	d["model_weight"] = 1.0
	if a.ridge_alpha <= 0 or not 0 <= a.conservation_min < a.conservation_max <= 1:
		raise ValueError("Invalid ridge/conservation bounds")
	# Remove unusable columns; require at least one ARG-specific predictor.
	bfeat = [x for x in BASE if x in d and d[x].notna().mean() > 0.01]
	efeat = [x for x in EVOL if x in d and d[x].notna().mean() > 0.01]
	if len(bfeat) < 3:
		raise SystemExit(f"Baseline MAF/LD/global-distance features incomplete: {bfeat}")
	if not any(x in efeat for x in ["log_arg_age", "local_arg_divergence", "root_time_gen"]):
		raise SystemExit(f"ARG age/local genealogy unavailable: {efeat}")
	ffeat = bfeat + efeat
	chroms = d.chr.astype(str).nunique()
	if chroms >= 3:
		groups = d.chr.astype(str).to_numpy()
		validation = "leave_one_chromosome_out"
	else:
		bp = pd.to_numeric(d.bp, errors="coerce")
		if not np.isfinite(bp).all():
			raise ValueError("Missing genomic positions; cannot form blocked folds")
		groups = (d.chr.astype(str) + "_" + ((bp // 5_000_000).astype(int) % 5).astype(str)).to_numpy()
		validation = "five_genomic_block_folds_pilot"
	pb, rb = fit_oof(d, bfeat, groups, a.ridge_alpha)
	pf, rf = fit_oof(d, ffeat, groups, a.ridge_alpha)
	y = d.transport_heterogeneity.to_numpy(float)
	mb = metrics(y, pb)
	mf = metrics(y, pf)
	selected = a.model  # Pre-specified; do not choose a model using its own reported OOF labels.
	ps = pf if selected == "evolutionary_full" else pb
	if not np.isfinite(ps).all():
		raise ValueError(
			"Incomplete blocked OOF predictions; use more SNPs/genomic blocks. Missing predictions must not become maximum conservation."
		)
	d["pred_baseline_oof"] = pb
	d["pred_full_oof"] = pf
	d["pred_selected_oof"] = ps
	# Variant prior is based exclusively on held-out predictions.
	v = d.groupby(["trait", "chr", "bp", "SNP", "A1", "A2"], dropna=False, as_index=False).agg(
		predicted_heterogeneity=("pred_selected_oof", "median"),
		observed_heterogeneity=("transport_heterogeneity", "median"),
		n_pairs=("transport_heterogeneity", "size"),
	)
	ph = np.maximum(pd.to_numeric(v.predicted_heterogeneity, errors="raise"), 0)
	v["conservation"] = np.exp(-0.5 * ph).clip(a.conservation_min, a.conservation_max)
	v["selected_model"] = selected
	v.to_csv(out / "variant_conservation.tsv.gz", sep="\t", index=False, compression="gzip", float_format="%.9g")
	d.to_csv(out / "pair_predictions.tsv.gz", sep="\t", index=False, compression="gzip", float_format="%.9g")
	cv = []
	for name, rec in [("baseline", rb), ("evolutionary_full", rf)]:
		for x in rec:
			cv.append({"model": name, "validation": validation, **x})
	cv += [
		{"model": "baseline", "validation": "all_oof", **mb},
		{"model": "evolutionary_full", "validation": "all_oof", **mf},
	]
	pd.DataFrame(cv).to_csv(out / "model_cv.tsv", sep="\t", index=False)
	models = {}
	coeff = []
	for name, features in [("baseline", bfeat), ("evolutionary_full", ffeat)]:
		b, med, mu, sd = full_fit(d, features, a.ridge_alpha)
		models[name] = {
			"features": features,
			"intercept": float(b[0]),
			"coefficients": [float(x) for x in b[1:]],
			"medians": {k: float(v) for k, v in med.items()},
			"means": [float(x) for x in mu],
			"scales": [float(x) for x in sd],
		}
		coeff.append(pd.DataFrame({"model": name, "feature": ["intercept"] + features, "standardized_coefficient": b}))
	pd.concat(coeff).to_csv(out / "coefficients.tsv", sep="\t", index=False)
	summary = {
		"n_rows": len(d),
		"n_variants": len(v),
		"validation": validation,
		"baseline_features": bfeat,
		"evolutionary_features": efeat,
		"ridge_alpha": a.ridge_alpha,
		"baseline_oof": mb,
		"evolutionary_oof": mf,
		"selected_model": selected,
		"rmse_improvement": mb["rmse"] - mf["rmse"],
	}
	(out / "transport_model.json").write_text(json.dumps({"summary": summary, "models": models}, indent=2) + "\n")
	print(json.dumps(summary))


def fit_transport_model_cli():
	fit_transport_model_main()


# 🚩 weights: make_grid_weights
"""Shrink ancestry-specific PRS-CSx effects toward a shared effect using the GRID prior."""
import argparse, math, re
from pathlib import Path
import numpy as np, pandas as pd

COMP = str.maketrans("ACGT", "TGCA")


def weights_complement(s):
	return str(s).translate(COMP)


def align_weights(beta, a1, a2, r1, r2):
	a1 = np.asarray(a1, str)
	a2 = np.asarray(a2, str)
	r1 = np.asarray(r1, str)
	r2 = np.asarray(r2, str)
	b = np.asarray(beta, float)
	same = (a1 == r1) & (a2 == r2)
	swap = (a1 == r2) & (a2 == r1)
	cs = np.array([weights_complement(x) for x in a1]) == r1
	cs &= np.array([weights_complement(x) for x in a2]) == r2
	cw = np.array([weights_complement(x) for x in a1]) == r2
	cw &= np.array([weights_complement(x) for x in a2]) == r1
	out = np.full(len(b), np.nan)
	out[same | cs] = b[same | cs]
	out[swap | cw] = -b[swap | cw]
	return out


def make_grid_weights_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--weights-dir", required=True)
	ap.add_argument("--conservation", required=True)
	ap.add_argument("--manifest", required=True)
	ap.add_argument("--out-dir", required=True)
	ap.add_argument("--pops", default="AFR,EAS,EUR,SAS")
	ap.add_argument("--chrs", required=True)
	a = ap.parse_args()
	pops = [x.upper() for x in re.split("[,; ]+", a.pops) if x]
	chrs = [x for x in re.split("[,; ]+", a.chrs) if x]
	out = Path(a.out_dir)
	out.mkdir(parents=True, exist_ok=True)
	man = pd.read_csv(a.manifest, sep="\t")
	ns = {str(r["pop"]).upper(): float(r["n_gwas"]) for _, r in man.iterrows()}
	nw = {p: math.sqrt(max(ns.get(p, 1), 1)) for p in pops}
	con = pd.read_csv(a.conservation, sep="\t", compression="infer", dtype={"SNP": str, "chr": str})[
		["SNP", "chr", "conservation"]
	].drop_duplicates(["SNP", "chr"])
	audits = []
	for c in chrs:
		ds = {}
		for p in pops:
			f = Path(a.weights_dir) / f"{p}.chr{c}.tsv"
			if not f.exists():
				continue
			d = pd.read_csv(f, sep="\t", dtype={"SNP": str})
			d["A1"] = d.A1.astype(str).str.upper()
			d["A2"] = d.A2.astype(str).str.upper()
			d["BETA"] = pd.to_numeric(d.BETA, errors="coerce")
			ds[p] = d.dropna(subset=["SNP", "A1", "BETA"]).drop_duplicates("SNP")
		if len(ds) < 2:
			raise SystemExit(f"Need >=2 PRS-CSx populations for chr{c}")
		priority = [p for p in ["EUR", "AFR", "EAS", "SAS"] if p in ds]
		ref = (
			pd.concat([ds[p][["SNP", "A1", "A2"]].assign(prio=i) for i, p in enumerate(priority)], ignore_index=True)
			.sort_values("prio")
			.drop_duplicates("SNP")
			.drop(columns="prio")
		)
		z = ref.copy()
		bet = []
		for p in pops:
			if p not in ds:
				z[f"beta_{p}"] = np.nan
				continue
			x = z[["SNP", "A1", "A2"]].merge(
				ds[p][["SNP", "A1", "A2", "BETA"]], on="SNP", how="left", suffixes=("_REF", "")
			)
			z[f"beta_{p}"] = align_weights(x.BETA, x.A1, x.A2, x.A1_REF, x.A2_REF)
		B = np.column_stack([z[f"beta_{p}"].to_numpy(float) for p in pops])
		W = np.array([nw[p] for p in pops], float)[None, :] * np.isfinite(B)
		shared = np.nansum(B * W, axis=1) / np.where(W.sum(1) > 0, W.sum(1), np.nan)
		z["beta_shared"] = shared
		z["chr"] = str(c)
		z = z.merge(con, on=["SNP", "chr"], how="left")
		z["prior_missing"] = z.conservation.isna()
		z["conservation"] = pd.to_numeric(z.conservation, errors="coerce").fillna(0).clip(0, 1)
		pd.DataFrame({"SNP": z.SNP, "A1": z.A1, "BETA": z.beta_shared}).dropna().to_csv(
			out / f"GRID_shared.chr{c}.tsv", sep="\t", index=False
		)
		for p in pops:
			bp = z[f"beta_{p}"].to_numpy(float)
			bp = np.where(np.isfinite(bp), bp, shared)
			bg = z.conservation.to_numpy() * shared + (1 - z.conservation.to_numpy()) * bp
			z[f"beta_GRID_{p}"] = bg
			pd.DataFrame({"SNP": z.SNP, "A1": z.A1, "BETA": bg}).dropna().to_csv(
				out / f"GRID_{p}.chr{c}.tsv", sep="\t", index=False
			)
		audits.append(z)
	pd.concat(audits, ignore_index=True).to_csv(
		out / "GRID_weight_audit.tsv.gz", sep="\t", index=False, compression="gzip", float_format="%.9g"
	)
	print(
		f"chromosomes={len(chrs)} variants={sum(len(x) for x in audits)} shared_weight=sqrt_n conservation_missing_is_zero"
	)


def make_grid_weights_cli():
	make_grid_weights_main()


# 🚩 mix: mix_population_scores
import argparse
from pathlib import Path
import numpy as np, pandas as pd

POPS = ["AFR", "EAS", "EUR", "SAS"]


def read_population_scores(path):
	return pd.read_csv(path, sep="\t", compression="infer", dtype={"eid": str, "IID": str, "#IID": str, "ID_2": str})


def mix_population_scores_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--scores", required=True)
	ap.add_argument("--ancestry", required=True)
	ap.add_argument("--prefix", required=True)
	ap.add_argument("--out", required=True)
	a = ap.parse_args()
	s = read_population_scores(a.scores)
	anc = read_population_scores(a.ancestry)
	idc = next((c for c in ["eid", "IID", "#IID", "ID_2"] if c in anc), None)
	if idc is None:
		raise SystemExit("No ancestry ID")
	anc = anc.rename(columns={idc: "eid"})
	cols = [f"{a.prefix}_{p}" for p in POPS]
	missing = [c for c in cols if c not in s]
	if missing:
		raise SystemExit("Missing population scores " + ",".join(missing))
	qcols = [f"posterior_{p}" for p in POPS]
	if all(c in anc for c in qcols):
		q = anc[["eid"] + qcols].copy()
		q.columns = ["eid"] + [f"q_{p}" for p in POPS]
	else:
		ac = next((c for c in ["ancestry", "genetic_ancestry", "predicted_ancestry"] if c in anc), None)
		if ac is None:
			raise SystemExit("No ancestry/posterior columns")
		q = anc[["eid", ac]].copy()
		for p in POPS:
			q[f"q_{p}"] = (q[ac].astype(str).str.upper() == p).astype(float)
		q = q.drop(columns=ac)
	if s.eid.isna().any() or q.eid.isna().any() or s.eid.duplicated().any() or q.eid.duplicated().any():
		raise ValueError("Missing/duplicate sample IDs")
	z = s.merge(q, on="eid", how="left", validate="one_to_one")
	Q = z[[f"q_{p}" for p in POPS]].apply(pd.to_numeric, errors="coerce").to_numpy()
	den = Q.sum(1)
	valid = np.isfinite(Q).all(1) & (Q >= 0).all(1) & (den > 0)
	Q = np.divide(Q, den[:, None], out=np.full_like(Q, np.nan), where=valid[:, None])
	B = z[cols].apply(pd.to_numeric, errors="coerce").to_numpy()
	z[f"{a.prefix}_posterior"] = np.sum(B * Q, axis=1)
	ix = np.argmax(np.nan_to_num(Q, nan=-1), axis=1)
	z[f"{a.prefix}_matched"] = np.where(valid, B[np.arange(len(B)), ix], np.nan)
	out = z[["eid"] + cols + [f"{a.prefix}_posterior", f"{a.prefix}_matched"]]
	Path(a.out).parent.mkdir(parents=True, exist_ok=True)
	out.to_csv(a.out, sep="\t", index=False, compression="gzip")
	print(len(out))


def mix_population_scores_cli():
	mix_population_scores_main()


COMMANDS = {
	"inputs": grid_inputs_cli,
	"ld": extract_ld_scores_cli,
	"transport": build_transport_table_cli,
	"fit": fit_transport_model_cli,
	"weights": make_grid_weights_cli,
	"mix": mix_population_scores_cli,
}


def main():
	import sys

	if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help", "help"):
		print(__doc__)
		print("Commands: " + ", ".join(COMMANDS))
		return
	command = sys.argv.pop(1)
	if command not in COMMANDS:
		raise SystemExit("Unknown command: " + command)
	sys.argv[0] += " " + command
	COMMANDS[command]()


if __name__ == "__main__":
	main()
