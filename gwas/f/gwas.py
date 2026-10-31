#!/usr/bin/env python3
"""gwas workflow utilities; use --help for commands."""

from __future__ import annotations


# 🚩 gwas_storage
"""Disposable SAIGE dosage exports shared within one invocation.

Association receipts bind durable source inputs, never disposable VCF mtimes.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def fingerprint(paths):
	return [(str(p), Path(p).stat().st_size, Path(p).stat().st_mtime_ns) for p in sorted(set(paths))]


def digest(value):
	return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def valid(path, signature, outputs):
	try:
		old = json.loads(path.read_text())
		return (
			old["signature"] == signature
			and old["outputs"] == [list(x) for x in fingerprint(outputs)]
			and all(Path(x).is_file() and Path(x).stat().st_size for x in outputs)
		)
	except (OSError, ValueError, KeyError):
		return False


def record(path, signature, outputs):
	data = dict(signature=signature, outputs=fingerprint(outputs))
	temp = path.with_suffix(".next")
	temp.write_text(json.dumps(data))
	temp.replace(path)


def execute_saige_chain(chain, spec, plan, stop=None):
	exports = [j for j in chain if j["name"].endswith((".export", ".index"))]
	associations = [j for j in chain if j["name"].endswith(".association")]
	if len(exports) != 2 or len(associations) != 1 or len(chain) != 3:
		raise ValueError("Unexpected SAIGE chromosome chain; regenerate gwas.sh plan")
	association = associations[0]
	temporary_outputs = {x for j in exports for x in j["outputs"]}
	durable = [x for j in chain for x in j["inputs"] if x not in temporary_outputs]
	# Script/executable changes must invalidate the receipt as well.
	dependencies = [__file__]
	for j in chain:
		executable = shutil.which(j["cmd"][0])
		if executable:
			dependencies.append(executable)
		dependencies += [x for x in j["cmd"][1:] if x.endswith(".R") and Path(x).is_file()]
	signature = digest([chain, fingerprint(durable + dependencies), "transient-dosage-v1"])
	receipt = Path(plan).parent / (association["name"] + ".transient.done.json")
	outputs = association["outputs"]
	if not spec.get("replace") and valid(receipt, signature, outputs):
		print("RESUME " + association["name"] + " without dosage export", flush=True)
		return
	receipt.unlink(missing_ok=True)
	# Remove the legacy association receipt before attempting replacement.
	(Path(plan).parent / (association["name"] + ".done.json")).unlink(missing_ok=True)
	shared = os.environ.get("GU_SAIGE_CACHE_DIR")
	if shared:
		run_with_cache(Path(shared), exports, association, plan, stop)
	else:
		# A standalone generated run.cmd owns and cleans its dosage cache.
		with tempfile.TemporaryDirectory(prefix=".saige-dosage.", dir=Path(plan).parent) as temporary:
			run_with_cache(Path(temporary), exports, association, plan, stop)
	if not all(Path(x).is_file() and Path(x).stat().st_size for x in outputs):
		raise RuntimeError("Missing SAIGE association result")
	record(receipt, signature, outputs)


def run_with_cache(root, exports, association, plan, stop=None):
	export, index = exports
	command = export["cmd"]
	old_prefix = command[command.index("--out") + 1]

	def rewrite(value, prefix):
		if value == old_prefix or value.startswith(old_prefix + "."):
			return prefix + value[len(old_prefix) :]
		if "=" in value:
			left, right = value.split("=", 1)
			if right == old_prefix or right.startswith(old_prefix + "."):
				return left + "=" + prefix + right[len(old_prefix) :]
		return value

	normalized = [[rewrite(x, "@dosage") for x in j["cmd"]] for j in exports]
	key = digest(
		[
			normalized,
			fingerprint(export["inputs"]),
			fingerprint([shutil.which(j["cmd"][0]) for j in exports]),
			"DS-force-v1",
		]
	)
	root.mkdir(parents=True, exist_ok=True)
	cache = root / key
	with (root / (key + ".lock")).open("a") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		prefix = str(cache / "dosage")
		generated = [rewrite(x, prefix) for j in exports for x in j["outputs"]]
		marker = cache / "complete.json"
		if not valid(marker, key, generated):
			shutil.rmtree(cache, ignore_errors=True)
			cache.mkdir()
			try:
				for job in exports:
					run(job, [rewrite(x, prefix) for x in job["cmd"]], plan, stop)
				if not all(Path(x).is_file() and Path(x).stat().st_size for x in generated):
					raise RuntimeError("Incomplete dosage export/index")
				record(marker, key, generated)
			except BaseException:
				shutil.rmtree(cache, ignore_errors=True)
				raise
		else:
			print("REUSE shared SAIGE dosage " + key, flush=True)
		run(association, [rewrite(x, prefix) for x in association["cmd"]], plan, stop)


def run(job, command, plan, stop=None):
	log = Path(plan).parent / (job["name"] + ".run.log")
	with log.open("w") as handle:
		proc = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT)
		try:
			while proc.poll() is None:
				if stop is not None and stop.is_set():
					raise RuntimeError("SAIGE run cancelled")
				try:
					proc.wait(timeout=0.2)
				except subprocess.TimeoutExpired:
					pass
			if proc.returncode:
				raise RuntimeError(f"{job['name']} failed; see {log}")
		finally:
			if proc.poll() is None:
				proc.terminate()
				try:
					proc.wait(timeout=10)
				except subprocess.TimeoutExpired:
					proc.kill()
					proc.wait()


# 🚩 gwas_io
"""Streaming, explicit effect-allele normalization; preserve chromosome ordering."""
import gzip
import math
import os
import re
from pathlib import Path

HEADER = "SNP CHR POS EA NEA EAF N BETA SE P".split()


def reader(path):
	return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def normalize(files, method, out):
	count = 0
	rejected = 0
	tmp = str(out) + ".tmp"
	with gzip.open(tmp, "wt") as target:
		target.write("\t".join(HEADER) + "\n")
		for file in files:
			with reader(file) as source:
				header = source.readline().strip().lstrip("#").split()
				for line in source:
					vals = line.split()
					if len(vals) != len(header):
						raise ValueError(f"Malformed row in {file}")
					d = dict(zip(header, vals))
					if d.get("TEST", "ADD") != "ADD":
						continue
					if method == "plink2":
						beta = d.get("BETA")
						if beta is None and "OR" in d:
							try:
								beta = str(math.log(float(d["OR"])))
							except (ValueError, OverflowError):
								beta = "NA"
						row = [
							d["ID"],
							d["CHROM"],
							d["POS"],
							d["A1"],
							d.get("OMITTED", d.get("AX", "NA")),
							d["A1_FREQ"],
							d["OBS_CT"],
							beta,
							d.get("SE", d.get("LOG(OR)_SE", "NA")),
							d["P"],
						]
					elif method == "regenie":
						# REGENIE ALLELE1 is effect allele; ALLELE0 is reference.
						try:
							pv = str(10 ** (-min(float(d["LOG10P"]), 300)))
						except ValueError:
							pv = "NA"
						row = [
							d["ID"],
							d["CHROM"],
							d["GENPOS"],
							d["ALLELE1"],
							d["ALLELE0"],
							d["A1FREQ"],
							d["N"],
							d["BETA"],
							d["SE"],
							pv,
						]
						# REGENIE collapses sex chromosomes to 23 internally.
						# Per-chromosome jobs retain the true Y identity here.
						if re.match(r"^chr(?:Y|24)[._]", Path(file).name) and row[1] in ("23", "X", "Y", "24"):
							row[1] = "24"
					elif method == "saige":
						n = d.get("N")
						if n is None:
							try:
								n = str(int(d["N_case"]) + int(d["N_ctrl"]))
							except (KeyError, ValueError):
								n = "NA"
						row = [
							d.get("MarkerID", d.get("SNPID")),
							d["CHR"],
							d["POS"],
							d["Allele2"],
							d["Allele1"],
							d["AF_Allele2"],
							n,
							d["BETA"],
							d["SE"],
							d["p.value"],
						]
					else:
						raise ValueError(method)
					try:
						nums = [float(row[j]) for j in (2, 5, 6, 7, 8, 9)]
						valid = (
							all(math.isfinite(x) for x in nums)
							and 0 <= nums[1] <= 1
							and nums[2] > 0
							and nums[4] > 0
							and 0 <= nums[5] <= 1
						)
						valid = valid and row[3] != row[4] and all(x not in (None, "NA", ".") for x in row[:5])
					except (ValueError, TypeError):
						valid = False
					if not valid:
						rejected += 1
						continue
					p = max(float(row[9]), 1e-300)
					if method == "regenie":
						row[9] = f"{p:.3f}" if p >= 0.001 else f"{p:.3E}"
					else:
						row[9] = str(p)
					target.write("\t".join(row) + "\n")
					count += 1
	if not count:
		os.unlink(tmp)
		raise ValueError(f"No valid additive results for {out}")
	os.replace(tmp, out)
	Path(str(out) + ".qc.txt").write_text(f"written\t{count}\nrejected\t{rejected}\n")


def gwas_io_cli():
	import sys

	normalize(sys.argv[3:], sys.argv[1], sys.argv[2])


# 🚩 gwas_runner
"""Generic plan/checkpoint and file-view helpers for gwas.sh.

GWAS options and commands are defined in the Bash entry point. Keep checkpoint
serialization stable so previously completed jobs remain resumable.
"""
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import re
import signal
import threading
import time
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = Path(__file__).resolve().parent
SELF = str(Path(__file__).resolve())
MODULE_ALIASES = {
	"account": "acount",
	"prep_gwas": "sparseGRM",
	"prep_regenie": "acount",
	"filter_variants": "acount",
	"run_regenie": "regenie",
}


def execute(plan):
	# Serialize jobs in one plan; never let two invocations overwrite a shared cache.
	import fcntl

	lock = open(str(plan) + ".lock", "w")
	try:
		fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
	except BlockingIOError:
		raise RuntimeError(f"Already running: {plan}")
	spec = json.loads(Path(plan).read_text())
	stop = threading.Event()
	failure_lock = threading.Lock()
	first_failure = None
	config = spec.get("config", {})
	is_regenie = MODULE_ALIASES.get(config.get("module"), config.get("module")) == "regenie"
	def job_environment(job):
		child_env = os.environ.copy()
		if not is_regenie:
			return child_env
		# Bound auxiliary BLAS/OpenMP/R pools as well as REGENIE's --threads.
		# Environment-only limits preserve existing command/checkpoint signatures.
		threads = config.get("threads", 1)
		cmd = job["cmd"]
		if "--step" in cmd and cmd[cmd.index("--step") + 1] == "2":
			threads = config.get("step2_threads", threads)
		for key in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
		            "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
		            "NUMEXPR_NUM_THREADS", "R_DATATABLE_NUM_THREADS"):
			child_env[key] = str(threads)
		return child_env

	def progress(event, label, detail=""):
		if is_regenie:
			# START/DONE/FAIL are forwarded by the shared screen-log filter.
			trait = Path(plan).parent.parent.parent.name
			when = datetime.now().astimezone().isoformat(timespec="seconds")
			print(f"[REGENIE] {event} {when} {trait} {label}{detail}", flush=True)

	def job_label(job):
		if "--step" in job["cmd"]:
			step = job["cmd"][job["cmd"].index("--step") + 1]
			return "Step 1" if step == "1" else "Step 2 " + job["name"].split(".")[0]
		return job["name"]

	if is_regenie:
		threads, workers = config.get("threads", 1), config.get("jobs", 1)
		step2_threads = config.get("step2_threads", threads)
		progress("START", f"任务；Step 1={threads} threads；Step 2 最多 {workers} 个染色体并行 × {step2_threads} threads；表型依次运行")
	for job in spec["jobs"]:
		executable = job["cmd"][0]
		if not shutil.which(executable):
			raise FileNotFoundError("Required executable: " + executable)
		if len(job["cmd"]) > 1 and job["cmd"][1].endswith(".R") and not Path(job["cmd"][1]).is_file():
			raise FileNotFoundError("Required R script: " + job["cmd"][1] + " (set --saige-dir/--saige-rscript)")

	def execute_job(job):
		for f in job["inputs"]:
			if not Path(f).is_file():
				raise FileNotFoundError(f)
		cmd = job["cmd"][:]
		if job["covars"]:
			path, method = job["covars"]
			names = Path(path).read_text().strip()
			if names:
				cmd += {
					"plink2": [
						"--covar",
						path.removesuffix(".covars"),
						"--covar-name",
						names,
						"--covar-variance-standardize",
					],
					"regenie": ["--covarFile", path.removesuffix(".covars"), "--covarColList", names],
					"saige": ["--covarColList=" + names],
				}[method]
		signature = hashlib.sha256(
			json.dumps(
				[cmd, [(f, Path(f).stat().st_size, Path(f).stat().st_mtime_ns) for f in job["inputs"]]], sort_keys=True
			).encode()
		).hexdigest()
		stamp = Path(plan).parent / (job["name"] + ".done.json")
		if (
			not spec["replace"]
			and stamp.exists()
			and json.loads(stamp.read_text()).get("signature") == signature
			and all(
				Path(f).is_file() and (Path(f).stat().st_size or str(f).endswith(".covars")) for f in job["outputs"]
			)
		):
			print("RESUME " + job["name"], flush=True)
			return True
		# Invalidate success before execution so failures cannot be mistaken for completion.
		stamp.unlink(missing_ok=True)
		print(shlex.join(cmd), flush=True)
		logpath = Path(plan).parent / (job["name"] + ".run.log")
		# Pgenlib code 4 is ReadFail, not malformed input. Reopen the file in a
		# fresh process after a bounded delay; never skip variants or append to a
		# partial association result. Keep commands/signatures stable for resume.
		retry_pgen = (is_regenie and "--pgen" in cmd and "--step" in cmd
		              and cmd[cmd.index("--step") + 1] == "2")
		for attempt in range(3 if retry_pgen else 1):
			if stop.is_set():
				raise RuntimeError("GWAS run cancelled")
			with logpath.open("w") as log:
				result = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=job_environment(job))
				try:
					while result.poll() is None:
						if stop.is_set():
							raise RuntimeError("GWAS run cancelled")
						try:
							result.wait(timeout=0.2)
						except subprocess.TimeoutExpired:
							pass
				finally:
					if result.poll() is None:
						result.terminate()
						try:
							result.wait(timeout=10)
						except subprocess.TimeoutExpired:
							result.kill()
							result.wait()
			if stop.is_set():
				raise RuntimeError("GWAS run cancelled")
			if result.returncode <= 0 or not retry_pgen or attempt == 2:
				break
			with logpath.open("rb") as log:
				log.seek(max(0, logpath.stat().st_size - 65536))
				tail = log.read().decode(errors="replace")
			if not re.search(r"PgrGet1D\(\) error 4\s*$", tail, re.MULTILINE):
				break
			archive = logpath.with_name(f"{job['name']}.retry-{time.time_ns()}-{attempt + 1}.run.log")
			shutil.copy2(logpath, archive)
			progress("START", job_label(job),
			         f" PGEN 读取失败；5s 后从头重试 {attempt + 1}/2；失败日志：{archive}")
			if stop.wait(5):
				raise RuntimeError("GWAS run cancelled")
		if result.returncode:
			print("\n".join(logpath.read_text(errors="replace").splitlines()[-12:]), file=sys.stderr)
			raise RuntimeError(f"{job['name']} failed (exit {result.returncode}); see {logpath}")
		if not all(
			Path(f).is_file() and (Path(f).stat().st_size or str(f).endswith(".covars")) for f in job["outputs"]
		):
			raise RuntimeError("Missing output: " + job["name"])
		stamp.write_text(json.dumps({"signature": signature, "command": cmd}))
		return False

	def run_job(job):
		nonlocal first_failure
		if stop.is_set():
			raise RuntimeError("GWAS run cancelled")
		started = time.monotonic()
		label = job_label(job)
		logpath = Path(plan).parent / (job["name"] + ".run.log")
		progress("START", label, f" 检查/执行；日志：{logpath}")
		try:
			reused = execute_job(job)
		except BaseException as error:
			# Signal queued workers here, before this worker can take another
			# chromosome while the main thread is still collecting the failure.
			with failure_lock:
				if first_failure is None:
					first_failure = error
				stop.set()
			progress("FAIL", label, f" 失败/中断；耗时 {time.monotonic() - started:.1f}s；日志：{logpath}")
			raise
		state = "复用已有完成结果" if reused else "完成"
		progress("DONE", label, f" {state}；耗时 {time.monotonic() - started:.1f}s")

	def run_chromosome(chain):
		if stop.is_set():
			return
		if spec.get("config", {}).get("module") == "run_saige":
			execute_saige_chain(chain, spec, plan, stop)
			return
		for job in chain:
			run_job(job)

	def flush_chromosomes(groups):
		if not groups:
			return
		started = time.monotonic()
		progress("START", f"Step 2 全部染色体（{len(groups)} 条）")
		# Each chromosome keeps its view/export/index/association dependency order.
		# Wait for every chromosome before merge-results or the next serial stage.
		with ThreadPoolExecutor(max_workers=spec.get("config", {}).get("jobs", 1)) as pool:
			futures = [pool.submit(run_chromosome, chain) for chain in groups.values()]
			try:
				for future in as_completed(futures):
					future.result()
			except BaseException as error:
				stop.set()
				for future in futures:
					future.cancel()
				progress("FAIL", "Step 2 全部染色体", f" 失败/中断；耗时 {time.monotonic() - started:.1f}s")
				# A cancelled sibling may finish first; retain the causal error.
				if isinstance(error, Exception) and first_failure is not None:
					raise first_failure
				raise
		progress("DONE", "Step 2 全部染色体", f" 完成（含复用）；耗时 {time.monotonic() - started:.1f}s")
		groups.clear()

	groups = {}
	try:
		for job in spec["jobs"]:
			match = re.match(r"^(chr(?:[0-9]+|X|Y))\.", job["name"])
			if match:
				groups.setdefault(match[1], []).append(job)
			else:
				flush_chromosomes(groups)
				run_job(job)
		flush_chromosomes(groups)
		progress("DONE", "任务完成")
	finally:
		lock.close()


def save_plan(directory, jobs, a):
	directory = Path(directory)
	directory.mkdir(parents=True, exist_ok=True)
	plan = directory / "run.plan.json"
	spec = {"config": vars(a), "replace": a.replace, "jobs": jobs}
	# Do not silently reuse results across different builds, subsets or model settings.
	if plan.exists():
		old = json.loads(plan.read_text())
		oc = old.get("config", {}).copy()
		nc = vars(a).copy()
		# Renaming a module does not change its analysis or invalidate checkpoints.
		for config in (oc, nc):
			config["module"] = MODULE_ALIASES.get(config.get("module"), config.get("module"))
			for key in ("continuous_covariate_file", "categorical_covariate_file"):
				if not config.get(key):
					config.pop(key, None)
		ignored = ["run", "replace", "sparse_grm", "saige_dir", "saige_rscript", "jobs"]
		# These filters are explicit commands/file inputs: signatures invalidate
		# Step 2 while preserving phenotype and Step 1 checkpoints on migration.
		ignored += ["global_maf"]
		if nc.get("module") == "regenie":
			ignored += ["min_mac", "step2_threads"]
		if nc.get("module") != "sparseGRM":
			# Per-trait checkpoints validate commands and actual input files;
			# selecting more chromosomes/traits can safely reuse completed jobs.
			ignored += ["chromosomes", "phenotypes"]
		for key in ignored:
			oc.pop(key, None)
			nc.pop(key, None)
		if oc != nc and not a.replace and any(directory.glob("*.done.json")):
			raise ValueError(f"Configuration changed at {directory}; choose a new output directory or --replace TRUE")
	plan.write_text(json.dumps(spec, indent=2) + "\n")
	lines = ["#!/usr/bin/env bash", "set -euo pipefail", f"# Commands and outputs: {plan}"]
	lines.extend("# " + shlex.join(j["cmd"]) for j in jobs)
	lines.append(shlex.join(["python3", SELF, "execute-plan", str(plan)]))
	script = directory / "run.cmd"
	script.write_text("\n".join(lines) + "\n")
	print(script, flush=True)
	if a.run:
		execute(plan)


def plan_node(args):
	name, covar_file, method = args[:3]
	sep = args.index("--", 3)
	out = args.index("--outputs", 3, sep)
	if args[3] != "--inputs":
		raise ValueError("Invalid plan-node arguments")
	job = dict(
		name=name,
		cmd=args[sep + 1 :],
		inputs=args[4:out],
		outputs=args[out + 1 : sep],
		covars=[covar_file, method] if covar_file else None,
	)
	print(json.dumps(job))


def finalize_plan(args):
	from types import SimpleNamespace

	directory, jobs_file = args[:2]
	config = dict(item.split("=", 1) for item in args[2:])
	config["module"] = MODULE_ALIASES.get(config["module"], config["module"])
	for key in ("threads", "memory", "min_mac", "jobs"):
		config[key] = int(config[key])
	if "step2_threads" in config:
		config["step2_threads"] = int(config["step2_threads"])
	for key in ("run", "replace", "sparse_grm"):
		config[key] = config[key] == "TRUE"
	for key in ("event_col", "extract", "keep"):
		if config[key] == "":
			config[key] = None
	config["chromosomes"] = config["chromosomes"].split(",")
	if config["module"] == "acount":
		config = {
			key: config[key]
			for key in (
				"module",
				"grch",
				"imputed_dir",
				"threads",
				"memory",
				"jobs",
				"plink2",
				"global_maf",
				"chromosomes",
				"run",
				"replace",
			)
		}
	with open(jobs_file) as source:
		jobs = [json.loads(line) for line in source if line.strip()]
	save_plan(directory, jobs, SimpleNamespace(**config))


def gwas_runner_cli():
	def interrupted(signum, frame):
		raise SystemExit(128 + signum)

	for sig in (signal.SIGTERM, signal.SIGHUP):
		signal.signal(sig, interrupted)
	try:
		if len(sys.argv) > 1 and sys.argv[1] == "execute-plan":
			execute(sys.argv[2])
		elif len(sys.argv) > 1 and sys.argv[1] == "build-metadata":
			Path(sys.argv[2] + ".grch").write_text(sys.argv[3] + "\n")
		elif len(sys.argv) > 1 and sys.argv[1] == "pgen-view":
			raise ValueError("pgen-view was removed; regenerate the plan with gwas.sh to use shared .pvar/.psam files")
		elif len(sys.argv) > 1 and sys.argv[1] == "create-grm":
			(actual, prefix) = sys.argv[2:4]
			subprocess.run(sys.argv[4:], check=True)
			pairs = [
				(Path(actual + suffix), Path(prefix + ".sparseGRM.mtx" + suffix)) for suffix in ("", ".sampleIDs.txt")
			]
			for src, dst in pairs:
				if not src.is_file() or not src.stat().st_size:
					raise RuntimeError("Missing GRM output: " + str(src))
			for src, dst in pairs:
				os.replace(src, dst)
		elif len(sys.argv) > 1 and sys.argv[1] == "grm-links":
			(actual, prefix) = sys.argv[2:]
			for src, dst in [
				(actual, prefix + ".sparseGRM.mtx"),
				(actual + ".sampleIDs.txt", prefix + ".sparseGRM.mtx.sampleIDs.txt"),
			]:
				p = Path(dst)
				if p.is_symlink():
					p.unlink()
				elif p.exists():
					raise FileExistsError("Refusing to replace existing GRM: " + dst)
				p.symlink_to(src)
		elif len(sys.argv) > 1 and sys.argv[1] == "plan-node":
			plan_node(sys.argv[2:])
		elif len(sys.argv) > 1 and sys.argv[1] == "finalize-plan":
			finalize_plan(sys.argv[2:])
		else:
			sys.exit(subprocess.call(["bash", str(HERE.parent / "gwas.sh")] + sys.argv[1:]))
	except (ValueError, RuntimeError, FileNotFoundError, subprocess.CalledProcessError) as e:
		print("ERROR: " + str(e), file=sys.stderr)
		sys.exit(1)


SUBCOMMANDS = {
	"merge": gwas_io_cli,
}


def main():
	import sys

	if len(sys.argv) > 1 and sys.argv[1] in SUBCOMMANDS:
		return SUBCOMMANDS[sys.argv.pop(1)]()
	if sys.argv[1:2] in (["--help"], ["-h"]):
		print("Helper commands: " + ", ".join(SUBCOMMANDS))
		print("Runner commands: execute-plan, plan-node, finalize-plan, create-grm, grm-links, build-metadata")
		return
	return gwas_runner_cli()


if __name__ == "__main__":
	main()
