#!/usr/bin/env python3
"""1.csx.py: signature, score-config, workspace, combined, posterior, normalize-weights."""

from __future__ import annotations
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

sys.dont_write_bytecode = True

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


# 🚩 signature: csx_config
"""Shared inference identity and validation for all PRS-CSx output modes."""
import json, math
from decimal import Decimal
from pathlib import Path
from grid_0_common import stamp, cache_directory, read_result_table, write_rds, rds_metadata

ROOT = Path(__file__).resolve().parent
POPS = ["AFR", "EAS", "EUR", "SAS"]


def find_gwas(directory, trait, pop):
	tags = [f"{trait}.{pop}"] + (["t2dm.AFA"] if trait == "t2dm" and pop == "AFR" else [])
	for tag in tags:
		for p in (Path(directory) / tag / "gwas" / f"{tag}.gz", Path(directory) / f"{tag}.gz"):
			if p.is_file() and p.stat().st_size:
				return p.resolve()
	raise FileNotFoundError(f"Missing GWAS: {trait}.{pop} in {directory}")


def validate_mcmc(phi, n, b, t, seed):
	if not (n > b >= 0 and t > 0 and n // t - b // t >= 2 and seed >= 0):
		raise ValueError("Require >=2 retained MCMC samples, positive thin, and a nonnegative seed")
	if phi != "auto" and not (math.isfinite(float(phi)) and float(phi) > 0):
		raise ValueError("Invalid phi")


def phi_label(phi):
	if phi == "auto":
		return "auto"
	value = Decimal(str(phi))
	if not value.is_finite() or value <= 0:
		raise ValueError("phi must be auto or positive")
	mantissa, exponent = format(value.normalize(), "e").split("e")
	return f"{mantissa}e{int(exponent)}"


def inference_name(phi, chrs):
	name = "auto" if phi == "auto" else "phi-" + phi_label(phi)
	chrs = sorted(map(int, chrs))
	if chrs != list(range(1, 23)):
		name += ".chr" + "-".join(map(str, chrs))
	return name


def configuration(settings, files):
	"""Readable cache identity; code edits alone do not invalidate computed results."""
	return json.dumps({"settings": settings, "files": [stamp(p) for p in files]}, sort_keys=True)


def inference_signature(phi, n, b, t, seed, override, chrs, gwas, snpinfo, bim, ld):
	validate_mcmc(phi, n, b, t, seed)
	return configuration(
		{"phi": phi_label(phi), "iterations": n, "burnin": b, "thin": t, "seed": seed,
		 "n_gwas": str(override), "chromosomes": sorted(map(int, chrs))},
		list(gwas) + [snpinfo, str(bim) + ".bim"] + list(ld),
	)


def workspace(parent, name, config):
	"""Reuse matching settings, keeping other runs under readable numbered names.

	Callers hold the trait's run.lock while selecting and using the workspace.
	An unlabelled or interrupted directory is never accepted as a cache hit.
	"""
	parent = Path(parent)
	if parent.name in {"scores", "posterior_scores"} or parent.parent.name == "combined_scores":
		parent = cache_directory(parent)
	expected = json.loads(config)
	# Keep run numbers distinct from chromosome lists such as .chr1-2.
	prefix = "run-" if name == "run" else name + ".run-"
	first_free = None
	candidates = [(1, parent / name)]
	if parent.is_dir():
		for path in parent.glob(prefix + "*"):
			suffix = path.name[len(prefix):]
			if suffix.isdigit() and int(suffix) >= 2:
				candidates.append((int(suffix), path))
	for number, path in sorted(candidates):
		if not path.exists():
			first_free = first_free or path
			continue
		marker = path / "config.json"
		if marker.is_file():
			try:
				if json.loads(marker.read_text()) == expected:
					return path
			except (ValueError, OSError):
				pass
	path = first_free or parent / f"{prefix}{max(n for n, _ in candidates) + 1}"
	path.mkdir(parents=True, exist_ok=False)
	marker = path / "config.json"
	tmp = marker.with_suffix(".tmp")
	tmp.write_text(json.dumps(expected, indent=2) + "\n")
	tmp.replace(marker)
	return path


def csx_workspace_cli():
	parent, name, config = sys.argv[1:]
	if name == "inference":
		settings = json.loads(config)["settings"]
		name = inference_name(settings["phi"], settings["chromosomes"])
	print(workspace(parent, name, config))


def csx_score_config_cli():
	args = sys.argv[1:]
	split = args.index("--files")
	print(configuration(args[:split], args[split + 1:]))


def sample_size(meta, override, pop):
	value = meta["n_gwas_median"]
	if override:
		value = (
			dict(x.split("=", 1) for x in override.replace(";", ",").replace(" ", ",").split(",") if x).get(pop, value)
			if "=" in override
			else override
		)
	if value is None or not math.isfinite(float(value)) or float(value) < 2:
		raise ValueError(f"Missing/invalid GWAS N: {pop}")
	return round(float(value))


def csx_config_cli():
	import sys

	phi, n, b, t, seed, override, chrs, snpinfo, bim, *paths = sys.argv[1:]
	print(
		inference_signature(
			phi, int(n), int(b), int(t), int(seed), override, chrs.split(), paths[:4], snpinfo, bim, paths[4:]
		)
	)


# 🚩 combined: csx_combined
"""Official PRS-CSx posterior meta scores; reuse the matching joint MCMC run."""
import argparse, concurrent.futures, fcntl, json, os, shutil, subprocess, threading
from pathlib import Path
import numpy as np
import pandas as pd
from grid_0_common import filter_samples, update_csx_table, score_source_records, sha256
from grid_0_common import inspect, coverage, stamp


def genotype(prefix):
	prefix = Path(prefix)
	if Path(str(prefix) + ".pgen").is_file():
		var = Path(str(prefix) + ".pvar")
		extra = []
		if not var.is_file():
			var = Path(str(prefix) + ".pvar.zst")
			extra = ["vzs"]
		files = [Path(str(prefix) + ".pgen"), Path(str(prefix) + ".psam"), var]
		command = ["--pfile", str(prefix), *extra]
	else:
		files = [Path(str(prefix) + ext) for ext in (".bed", ".bim", ".fam")]
		command = ["--bfile", str(prefix)]
	for f in files:
		stamp(f)
	return command, files


def csx_combined_main():
	p = argparse.ArgumentParser()
	for key in ["trait", "mode", "gwas-dir", "target-dir", "snpinfo", "ref-dir", "bim", "work", "score-home", "chrs"]:
		p.add_argument("--" + key, required=True)
	for key, default in [
		("jobs", 4),
		("threads", 1),
		("iterations", 4000),
		("burnin", 2000),
		("thin", 5),
		("seed", 20260904),
	]:
		p.add_argument("--" + key, type=int, default=default)
	p.add_argument("--score-jobs", type=int)
	p.add_argument("--score-memory", type=int, default=2048)
	p.add_argument("--phi", default="1e-2")
	p.add_argument("--n-gwas", default="")
	p.add_argument("--remove", default="")
	p.add_argument("--keep", default="")
	p.add_argument("--stage", choices=["all", "weights", "score"], default="all")
	p.add_argument("--replace", choices=["TRUE", "FALSE"], default="FALSE")
	p.add_argument("--check", action="store_true")
	a = p.parse_args()
	if a.mode not in ("auto", "meta") or min(a.jobs, a.threads) < 1:
		raise ValueError("Invalid mode/jobs/threads")
	if a.score_jobs is None:
		a.score_jobs = a.jobs
	if a.score_jobs < 1 or a.score_memory < 640:
		raise ValueError("Invalid score-jobs/score-memory")
	phi = "auto" if a.mode == "auto" else a.phi
	if a.mode == "meta" and phi == "auto":
		raise ValueError("meta requires fixed phi")
	validate_mcmc(phi, a.iterations, a.burnin, a.thin, a.seed)
	chrs = list(map(int, a.chrs.split()))
	if not chrs or len(set(chrs)) != len(chrs) or not set(chrs) <= set(range(1, 23)):
		raise ValueError("Invalid chromosomes")
	inputs = [find_gwas(a.gwas_dir, a.trait, pop) for pop in POPS]
	name = Path(a.snpinfo).name
	if name not in ("snpinfo_mult_1kg_hm3", "snpinfo_mult_ukbb_hm3"):
		raise ValueError("Unrecognized SNPINFO reference type")
	ref_type = "1kg" if name == "snpinfo_mult_1kg_hm3" else "ukbb"
	refs = []
	ld = []
	for pop in POPS:
		d = Path(a.ref_dir) / f"ldblk_{ref_type}_{pop.lower()}"
		if not d.is_dir():
			d = Path(a.ref_dir) / f"ldblk_{ref_type}_{pop}"
		refs.append(d.resolve())
		ld += [d / f"ldblk_{ref_type}_chr{c}.hdf5" for c in chrs]
	sig = inference_signature(phi, a.iterations, a.burnin, a.thin, a.seed, a.n_gwas, chrs, inputs, a.snpinfo, a.bim, ld)
	inspect(a.snpinfo, list(map(str, inputs)))
	coverage(a.chrs, list(map(str, inputs)))
	gen = {}
	targetfiles = []
	if a.stage != "weights":
		if not shutil.which("plink2"):
			raise ValueError("plink2 not found")
		for c in chrs:
			gen[c], files = genotype(Path(a.target_dir) / f"chr{c}")
			targetfiles += files
		for f in (a.remove, a.keep):
			if f and not Path(f).is_file():
				raise FileNotFoundError(f)
	home = Path(a.score_home)
	work = Path(a.work)
	weight = home / ".weights" / f"csx.{a.mode}.gz"
	meta = weight.with_suffix(".json")
	ready = weight.is_file() and meta.is_file() and json.loads(meta.read_text()).get("signature") == sig
	if a.stage == "score" and not ready:
		raise ValueError("Matching combined weights missing; run --stage all or weights")
	print(f"CSx {a.mode}: phi={phi}; joint run={inference_name(phi, chrs)}", flush=True)
	if a.check:
		return
	(work / a.trait).mkdir(parents=True, exist_ok=True)
	lock_dir = cache_directory(work / a.trait)
	lock_dir.mkdir(parents=True, exist_ok=True)
	lock = (lock_dir / "run.lock").open("a")
	fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
	run = workspace(work / a.trait, inference_name(phi, chrs), sig)
	scratch = cache_directory(run)
	logs = scratch / "logs"
	logs.mkdir(parents=True,exist_ok=True)
	command_lock = threading.Lock()
	env = dict(
		os.environ, OMP_NUM_THREADS=str(a.threads), OPENBLAS_NUM_THREADS=str(a.threads), MKL_NUM_THREADS=str(a.threads)
	)

	def execute(cmd, name):
		cmd = list(map(str, cmd))
		with command_lock:
			with (scratch / "commands.jsonl").open("a") as f:
				f.write(json.dumps(cmd) + "\n")
		with (logs / (name + ".log")).open("w") as f:
			rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env).returncode
		if rc:
			raise RuntimeError(f"Command exit={rc}; see {logs / (name + '.log')}")

	def posterior(c):
		found = list((run / "raw" / f"chr{c}").glob("*META*pst_eff*.txt"))
		return found[0] if len(found) == 1 and found[0].stat().st_size else None

	if a.stage != "score" and (not ready or a.replace == "TRUE"):
		raw = run / "raw"
		raw.mkdir(exist_ok=True)
		done = all(posterior(c) is not None and (raw / f"chr{c}" / "done").is_file() for c in chrs)
		if done and a.replace == "FALSE":
			print("SKIP inference: reuse matching joint population/META posteriors", flush=True)
		else:
			ref = scratch / "reference"
			ref.mkdir(exist_ok=True)
			for src, nm in [(Path(a.snpinfo).resolve(), name)] + [
				(d, f"ldblk_{ref_type}_{pop.lower()}") for d, pop in zip(refs, POPS)
			]:
				dest = ref / nm
				if not dest.exists():
					dest.symlink_to(src)
			prep = cache_directory(run / "sumstats")
			prep.mkdir(parents=True,exist_ok=True)
			sizes = []
			for pop, src in zip(POPS, inputs):
				dst = prep / f"{pop}.tsv.gz"
				info = prep / f"{pop}.json"
				execute(
					[
						"python3",
						ROOT / "0.common.py",
						"sumstats-cache",
						"--input",
						src,
						"--output",
						dst,
						"--metadata",
						info,
						"--snpinfo",
						a.snpinfo,
						"--trait",
						a.trait,
						"--pop",
						pop,
						"--work",
						work,
						"--replace",
						a.replace,
					],
					f"prepare.{pop}",
				)
				print((logs / f"prepare.{pop}.log").read_text().splitlines()[-1], flush=True)
				sizes.append(str(sample_size(json.loads(info.read_text()), a.n_gwas, pop)))
				execute(
					[
						"python3",
						ROOT / "0.common.py",
						"split-sumstats",
						"--input",
						dst,
						"--out-dir",
						prep,
						"--prefix",
						pop,
						"--chrs",
						a.chrs,
					],
					f"split.{pop}",
				)

			def infer(c):
				dest = raw / f"chr{c}"
				dest.mkdir(exist_ok=True)
				marker = dest / "done"
				if marker.exists() and posterior(c) is not None and a.replace == "FALSE":
					return
				marker.unlink(missing_ok=True)
				sst = ",".join(str(prep / f"{pop}.chr{c}.tsv") for pop in POPS)
				cmd = [
					"python3",
					ROOT / "csx/PRScsx.py",
					f"--ref_dir={ref}",
					f"--bim_prefix={a.bim}",
					f"--sst_file={sst}",
					f"--n_gwas={','.join(sizes)}",
					f"--pop={','.join(POPS)}",
					f"--chrom={c}",
					f"--n_iter={a.iterations}",
					f"--n_burnin={a.burnin}",
					f"--thin={a.thin}",
					f"--seed={a.seed + c}",
					f"--out_dir={dest}",
					f"--out_name={a.trait}",
					"--meta=TRUE",
				]
				if phi != "auto":
					cmd.append(f"--phi={phi}")
				print(f"RUN {a.trait} {a.mode} chr{c}", flush=True)
				execute(cmd, f"infer.chr{c}")
				if posterior(c) is None:
					raise ValueError(f"Missing META posterior chr{c}")
				marker.touch()

			with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
				list(ex.map(infer, chrs))
		z = pd.concat([read(posterior(c)) for c in chrs], ignore_index=True)
		if z.empty or z.SNP.duplicated().any() or not np.isfinite(z.BETA).all():
			raise ValueError("Invalid combined weights")
		weight.parent.mkdir(parents=True, exist_ok=True)
		tmp = weight.with_suffix(".tmp")
		z.to_csv(tmp, sep="\t", index=False, compression="gzip")
		tmp.replace(weight)
		meta.write_text(json.dumps({"signature": sig, "phi": phi, "inputs": list(map(str, inputs))}, indent=2) + "\n")
	if a.stage == "weights":
		return
	scoring_sig = configuration(
		[sig, a.mode, {"weights_sha256": sha256(weight)}],
		targetfiles + [weight] + [Path(x) for x in (a.remove, a.keep) if x],
	)
	score = workspace(run / "combined_scores" / a.mode, "run", scoring_sig)

	def scoring(c):
		out = score / f"chr{c}"
		done = Path(str(out) + ".done")
		if done.exists() and Path(str(out) + ".sscore").is_file() and a.replace == "FALSE":
			return
		done.unlink(missing_ok=True)
		cmd = [
			"plink2",
			*gen[c],
			"--score",
			weight,
			"1",
			"2",
			"3",
			"header-read",
			"no-mean-imputation",
			"list-variants",
			"cols=+scoresums",
			"--threads",
			a.threads,
			"--memory",
			a.score_memory,
			"--out",
			out,
		]
		if a.remove and Path(a.remove).stat().st_size:
			cmd += ["--remove", a.remove]
		if a.keep:
			cmd += ["--keep", a.keep]
		execute(cmd, f"{a.mode}.score.chr{c}")
		if not Path(str(out) + ".sscore").is_file():
			raise ValueError(f"Missing scores: {out}")
		done.touch()

	with concurrent.futures.ThreadPoolExecutor(a.score_jobs) as ex:
		list(ex.map(scoring, chrs))
	combined = score / "combined.gz"
	execute(
		[
			"python3",
			ROOT / "0.common.py",
			"combine-scores",
			"--inputs",
			*[score / f"chr{c}.sscore" for c in chrs],
			"--name",
			f"csx.{a.mode}",
			"--output",
			combined,
		],
		f"{a.mode}.combine",
	)
	z = filter_samples(pd.read_csv(combined, sep="\t", dtype={"eid": str}), "eid", a.remove)
	output = home / "1csx.scores.rds"
	update_csx_table(z, output, a.remove, provenance=score_source_records({f"csx.{a.mode}": weight}, sig, chrs))
	print(f"DONE {a.trait}: csx.{a.mode}; N={len(z)}; {output}", flush=True)


def csx_combined_cli():
	csx_combined_main()


# 🚩 posterior: csx_posterior
"""Score synchronized PRS-CSx draws; retain individual joint covariance.

One chromosome at a time. Cross-chromosome covariances are exactly zero in
the fitted chromosome-factorized model; iteration numbers from independent
chromosomes must NOT be treated as coupled joint draws.
"""
import argparse, hashlib, json, os, subprocess
from pathlib import Path
import h5py
import numpy as np
import pandas as pd

PAIRS = [(j, k) for j in range(4) for k in range(j, 4)]
COMP = str.maketrans("ACGT", "TGCA")


def align_frequencies(g, table):
	"""EAF refers to table A1; orient it to the scored A1, never use MAF."""
	if not {"SNP", "A1", "A2", "EAF"} <= set(table):
		raise ValueError("Frequency table needs SNP,A1,A2,EAF")
	if table.SNP.duplicated().any():
		raise ValueError("Duplicate frequency SNPs")
	ids = g["SNP"].asstr()[:]
	d = table.set_index("SNP").reindex(ids)
	af = pd.to_numeric(d.EAF, errors="coerce").to_numpy()
	a1, a2 = g["A1"].asstr()[:], g["A2"].asstr()[:]
	f1, f2 = d.A1.fillna("").str.upper().to_numpy(), d.A2.fillna("").str.upper().to_numpy()
	c1 = np.array([s.translate(COMP) for s in f1])
	c2 = np.array([s.translate(COMP) for s in f2])
	same = ((a1 == f1) & (a2 == f2)) | ((a1 == c1) & (a2 == c2))
	flip = ((a1 == f2) & (a2 == f1)) | ((a1 == c2) & (a2 == c1))
	valid = (same ^ flip) & np.isfinite(af) & (af > 0) & (af < 1)
	if not valid.all():
		raise ValueError(
			f"{np.sum(~valid)} posterior SNPs lack unambiguous GWAS EAF; first: {ids[~valid][:5].tolist()}"
		)
	return np.where(flip, 1 - af, af)


def pack_scores(sscore, dest, ndraw):
	"""Convert PLINK score sums to row-chunked HDF5, with strict ID checks."""
	header = pd.read_csv(sscore, sep=r"\s+", nrows=0).columns
	idcol = "IID" if "IID" in header else "#IID"
	cols = [f"DRAW_{i:04d}_SUM" for i in range(ndraw)]
	if not set(cols + [idcol]) <= set(header):
		raise ValueError("PLINK posterior score columns are incomplete")
	n = 0
	seen = set()
	with h5py.File(str(dest) + ".tmp", "w") as h:
		values = h.create_dataset(
			"scores", (0, ndraw), maxshape=(None, ndraw), chunks=(256, ndraw), dtype="f8", compression="lzf"
		)
		ids = h.create_dataset("eid", (0,), maxshape=(None,), dtype=h5py.string_dtype())
		for d in pd.read_csv(sscore, sep=r"\s+", usecols=[idcol] + cols, dtype={idcol: str}, chunksize=512):
			ix = d[idcol]
			if ix.isna().any() or (ix == "").any() or ix.duplicated().any() or seen.intersection(ix):
				raise ValueError("Missing/duplicate posterior IDs")
			ix = ix.astype(str)
			seen.update(ix)
			v = d[cols].to_numpy(float)
			if not np.isfinite(v).all():
				raise ValueError("Nonfinite posterior scores")
			values.resize(n + len(d), axis=0)
			ids.resize(n + len(d), axis=0)
			values[n : n + len(d)] = v
			ids[n : n + len(d)] = ix.to_numpy()
			n += len(d)
		if not n:
			raise ValueError("No posterior score participants")
		h.attrs["complete"] = True
	os.replace(str(dest) + ".tmp", dest)


def moments(pop_files, output):
	handles = [h5py.File(p, "r") for p in pop_files]
	try:
		shape = handles[0]["scores"].shape
		if shape[1] < 20 or any(h["scores"].shape != shape or not h.attrs.get("complete") for h in handles):
			raise ValueError("Need >=20 aligned draws and identical score dimensions")
		n, b = shape
		with h5py.File(str(output) + ".tmp", "w") as out:
			out.create_dataset("eid", shape=(n,), dtype=h5py.string_dtype())
			out.create_dataset("mean", shape=(n, 4), dtype="f8", chunks=True, compression="lzf")
			out.create_dataset("cov", shape=(n, 10), dtype="f8", chunks=True, compression="lzf")
			for start in range(0, n, 512):
				sl = slice(start, min(n, start + 512))
				ids = handles[0]["eid"].asstr()[sl]
				if any(not np.array_equal(ids, h["eid"].asstr()[sl]) for h in handles[1:]):
					raise ValueError("Population score sample order differs")
				v = np.stack([h["scores"][sl] for h in handles], axis=2)
				mean = v.mean(axis=1)
				u = v - mean[:, None, :]
				cov = np.einsum("nbj,nbk->njk", u, u) / (b - 1)
				out["eid"][sl] = ids
				out["mean"][sl] = mean
				out["cov"][sl] = np.column_stack([cov[:, j, k] for j, k in PAIRS])
			out.attrs.update(complete=True, ndraw=b)
		os.replace(str(output) + ".tmp", output)
	finally:
		for h in handles:
			h.close()


def combine_chromosomes(files):
	handles = [h5py.File(p, "r") for p in files]
	try:
		n = len(handles[0]["eid"])
		if any(len(h["eid"]) != n or not h.attrs.get("complete") for h in handles):
			raise ValueError("Incomplete chromosome moments")
		frames = []
		for start in range(0, n, 4096):
			sl = slice(start, min(start + 4096, n))
			ids = handles[0]["eid"].asstr()[sl]
			if any(not np.array_equal(ids, h["eid"].asstr()[sl]) for h in handles[1:]):
				raise ValueError("Chromosome score sample order differs")
			mu = sum(h["mean"][sl] for h in handles)
			cov = sum(h["cov"][sl] for h in handles)
			d = pd.DataFrame({"eid": ids})
			for j, pop in enumerate(POPS):
				d[f"csx.{pop}"] = mu[:, j]
			for j, (pop, other) in enumerate(PAIRS):
				d[f"cov.{POPS[pop]}.{POPS[other]}"] = cov[:, j]
			frames.append(d)
		return pd.concat(frames, ignore_index=True)
	finally:
		for handle in handles:
			handle.close()


def csx_posterior_main():
	p = argparse.ArgumentParser(description=__doc__)
	for key in ("raw-dir", "sumstats-dir", "target-dir", "output", "chrs"):
		p.add_argument("--" + key, required=True)
	p.add_argument(
		"--frequency-dir", help="Optional POP.tsv.gz files with discovery SNP,A1,A2,EAF; defaults to normalized GWAS"
	)
	p.add_argument("--threads", type=int, default=4)
	p.add_argument("--memory", type=int, default=8192)
	p.add_argument("--remove", default="")
	p.add_argument("--keep", default="")
	p.add_argument("--replace", choices=["TRUE", "FALSE"], default="FALSE")
	p.add_argument("--min-draws", type=int, default=100)
	a = p.parse_args()
	chrs = list(map(int, a.chrs.split()))
	if not chrs or len(set(chrs)) != len(chrs) or not set(chrs) <= set(range(1, 23)):
		raise ValueError("Invalid chromosomes")
	if min(a.threads, a.memory, a.min_draws) < 1:
		raise ValueError("Invalid resources/draw count")
	output = Path(a.output)
	output.parent.mkdir(parents=True, exist_ok=True)
	frequency = {}
	freq_paths = []
	for pop in POPS:
		path = Path(a.frequency_dir or a.sumstats_dir) / f"{pop}.tsv.gz"
		freq_paths.append(path)
		frequency[pop] = pd.read_csv(
			path, sep="\t", usecols=["SNP", "A1", "A2", "EAF"], dtype={"SNP": str, "A1": str, "A2": str}
		)
	raw = [Path(a.raw_dir) / f"chr{c}" / "joint_posterior.h5" for c in chrs]
	gen = {}
	gen_files = []
	for c in chrs:
		gen[c], fs = genotype(Path(a.target_dir) / f"chr{c}")
		gen_files.extend(fs)
	sig = configuration(
		[chrs, a.min_draws],
		raw + freq_paths + gen_files + [Path(x) for x in (a.keep, a.remove) if x],
	)
	if a.replace == "FALSE" and output.is_file() and rds_metadata(output).get("signature") == sig:
		print("SKIP posterior scores: matching individual covariance")
		return
	work = workspace(Path(a.raw_dir).parent / "posterior_scores", "run", sig)
	chr_files = []
	draw_counts = []
	for c, source in zip(chrs, raw):
		dst = work / f"chr{c}"
		dst.mkdir(exist_ok=True)
		moment = dst / "moments.h5"
		chr_files.append(moment)
		with h5py.File(source, "r") as h:
			if not h.attrs.get("complete") or h.attrs.get("schema") != "grid_csx_draws_v1":
				raise ValueError(f"Incomplete draws: {source}")
			ndraw = len(h["iteration"])
			draw_counts.append(ndraw)
			if ndraw < a.min_draws:
				raise ValueError(f"chr{c}: {ndraw} draws < {a.min_draws}")
			if moment.is_file() and a.replace == "FALSE":
				continue
			print(f"START posterior scoring chr{c}: {ndraw} joint draws", flush=True)
			score_files = []
			for pop in POPS:
				g = h[pop]
				eaf = align_frequencies(g, frequency[pop])
				ids = g["SNP"].asstr()[:]
				w = dst / f"{pop}.weights.tsv"
				af = dst / f"{pop}.afreq"
				pref = dst / pop
				cache = dst / f"{pop}.scores.h5"
				score_files.append(cache)
				if cache.is_file() and a.replace == "FALSE":
					continue
				# Fixed discovery allele means; never re-centre within the target group.
				pd.DataFrame(
					{
						"#CHROM": c,
						"ID": ids,
						"REF": g["A2"].asstr()[:],
						"ALT": g["A1"].asstr()[:],
						"ALT_FREQS": eaf,
						"OBS_CT": 2,
					}
				).to_csv(af, sep="\t", index=False)
				columns = ["SNP", "A1"] + [f"DRAW_{i:04d}" for i in range(ndraw)]
				for st in range(0, len(ids), 2048):
					sl = slice(st, min(st + 2048, len(ids)))
					v = pd.DataFrame(g["beta"][sl], columns=columns[2:])
					v.insert(0, "A1", g["A1"].asstr()[sl])
					v.insert(0, "SNP", ids[sl])
					v.to_csv(
						w, sep="\t", index=False, mode="w" if st == 0 else "a", header=st == 0, float_format="%.12g"
					)
				cmd = [
					"plink2",
					*gen[c],
					"--extract",
					str(w),
					"--read-freq",
					str(af),
					"--error-on-freq-calc",
					"--score",
					str(w),
					"1",
					"2",
					"header-read",
					"center",
					"no-mean-imputation",
					"list-variants",
					"cols=maybefid,scoresums",
					"--score-col-nums",
					f"3-{ndraw + 2}",
					"--threads",
					str(a.threads),
					"--memory",
					str(a.memory),
					"--out",
					str(pref),
				]
				# Centred missing genotypes contribute zero, mathematically equal to
				# discovery-mean imputation. Read SUM, never AVG (its denominator
				# varies with missingness). Explicit no-mean-imputation also avoids
				# old PLINK builds adding an uncentred mean for missing calls.
				# --extract must contain only IDs, without a header or effect columns.
				snps = dst / f"{pop}.snps"
				snps.write_text("\n".join(ids) + "\n")
				cmd[cmd.index("--extract") + 1] = str(snps)
				if a.keep:
					cmd += ["--keep", a.keep]
				if a.remove and Path(a.remove).stat().st_size:
					cmd += ["--remove", a.remove]
				with (dst / f"{pop}.command.json").open("w") as f:
					json.dump(cmd, f)
				with (dst / f"{pop}.run.log").open("w") as f:
					subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)
				used = Path(str(pref) + ".sscore.vars").read_text().split()
				if len(used) != len(ids) or set(used) != set(ids):
					raise ValueError(f"{pop} chr{c}: PLINK skipped posterior SNPs; cannot claim complete uncertainty")
				pack_scores(str(pref) + ".sscore", cache, ndraw)
				for tmp in (w, af, snps, Path(str(pref) + ".sscore")):
					tmp.unlink(missing_ok=True)
			moments(score_files, moment)
			for file in score_files:
				file.unlink(missing_ok=True)
		print(f"DONE posterior scoring chr{c}", flush=True)
	combined = combine_chromosomes(chr_files)
	metadata = {
		"schema": "grid_csx_moments_v1",
		"signature": sig,
		"populations": POPS,
		"chromosomes": chrs,
		"draw_counts": draw_counts,
		"centering": "discovery_EAF",
		"covariance": "joint within chromosome; summed across independent chromosomes",
		"means": "posterior mean scores centred at discovery EAF",
		"effect_scale": "standardized_phenotype_per_allele",
	}
	write_rds(combined, output, metadata, float_format="%.12g")
	print(f"DONE individual posterior moments: {output}", flush=True)


def csx_posterior_cli():
	csx_posterior_main()


# 🚩 normalize-weights: normalize_csx_weights
import argparse, re
from pathlib import Path
import pandas as pd
import numpy as np


def read(path):
	# Standard PRS-CSx posterior files are whitespace-delimited, often without a header.
	with open(path, errors="replace") as h:
		first = h.readline().split()
	known = {x.upper() for x in first} & {"CHR", "SNP", "BP", "A1", "A2", "BETA"}
	if len(known) >= 3:
		d = pd.read_csv(path, sep=r"\s+", engine="python")
	else:
		d = pd.read_csv(path, sep=r"\s+", engine="python", header=None)
		if d.shape[1] < 6:
			raise SystemExit(f"Unrecognized PRS-CSx output {path}: {d.shape[1]} columns")
		d = d.iloc[:, :6]
		d.columns = ["CHR", "SNP", "BP", "A1", "A2", "BETA"]
	up = {str(c).upper(): c for c in d.columns}

	def c(*x):
		return next((up[y] for y in x if y in up), None)

	sc = c("SNP", "RSID", "ID")
	ac = c("A1", "EA", "EFFECT_ALLELE")
	bc = c("BETA", "POSTERIOR_BETA", "EFFECT")
	if not sc or not ac or not bc:
		raise SystemExit(f"Missing SNP/A1/BETA in {path}: {list(d.columns)}")
	out = pd.DataFrame(
		{"SNP": d[sc].astype(str), "A1": d[ac].astype(str).str.upper(), "BETA": pd.to_numeric(d[bc], errors="coerce")}
	)
	for new, aliases in [("CHR", ("CHR", "CHROM")), ("BP", ("BP", "POS")), ("A2", ("A2", "NEA"))]:
		z = c(*aliases)
		out[new] = d[z] if z else pd.NA
	if out.empty or out.SNP.duplicated().any() or not np.isfinite(out.BETA).all():
		raise ValueError(f"Invalid posterior weights: {path}")
	return out[["SNP", "A1", "BETA", "CHR", "BP", "A2"]]


def normalize_csx_weights_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--input", required=True)
	ap.add_argument("--output", required=True)
	a = ap.parse_args()
	d = read(a.input)
	d.to_csv(a.output, sep="\t", index=False)
	print(len(d))


def normalize_csx_weights_cli():
	normalize_csx_weights_main()


COMMANDS = {
	"signature": csx_config_cli,
	"score-config": csx_score_config_cli,
	"workspace": csx_workspace_cli,
	"combined": csx_combined_cli,
	"posterior": csx_posterior_cli,
	"normalize-weights": normalize_csx_weights_cli,
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
