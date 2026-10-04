#!/usr/bin/env python3
"""Read exact workbook exports and exchange named RDS tables through /tmp."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import sys
import zipfile


sys.dont_write_bytecode = True


# 🚩 R table storage
def run_r(*arguments):
	environment = os.environ.copy()
	rscript = environment.get("RESULTS_RSCRIPT", "/usr/bin/Rscript")
	if Path(rscript).resolve() == Path("/usr/bin/Rscript").resolve():
		for name in ("R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_ENVIRON_USER"):
			environment.pop(name, None)
	with tempfile.TemporaryDirectory(prefix = "result-code-", dir = "/tmp") as directory:
		script = Path(directory) / "results.R"
		script.write_bytes(Path(__file__).with_suffix(".R").read_bytes())
		return subprocess.run(
			[rscript, str(script), *map(str, arguments)],
			check = True, capture_output = True, text = True, env = environment,
		).stdout.strip()


def write_rds(table, path, metadata = None, float_format = "%.17g"):
	path = Path(path)
	if path.suffix.lower() != ".rds":
		raise ValueError("Participant results require a descriptive .rds filename")
	with tempfile.TemporaryDirectory(prefix = "result-write-", dir = "/tmp") as directory:
		if isinstance(table, dict):
			source = Path(directory) / "model.json"
			source.write_text(json.dumps(table), encoding = "utf-8")
			run_r("import-object", source, path)
			return path
		source = Path(directory) / (path.stem + ".tsv.gz")
		table.to_csv(source, sep = "\t", index = False, compression = "gzip", float_format = float_format)
		arguments = ["import-table", source, path]
		if metadata is not None:
			meta = Path(directory) / "metadata.json"
			meta.write_text(json.dumps(metadata), encoding = "utf-8")
			arguments.append(meta)
		run_r(*arguments)
	return path


def rds_metadata(path):
	return json.loads(run_r("metadata", path))


def write_workbook(tables, path):
	with tempfile.TemporaryDirectory(prefix = "result-book-", dir = "/tmp") as directory:
		sources = {}
		for name, table in tables.items():
			if table.empty:
				continue
			source = Path(directory) / (name + ".tsv")
			if source.parent != Path(directory):
				raise ValueError("Worksheet names cannot contain directory separators")
			table.to_csv(source, sep = "\t", index = False, float_format = "%.17g")
			sources[name] = str(source)
		spec = Path(directory) / "workbook.json"
		spec.write_text(json.dumps({"tables": sources}), encoding = "utf-8")
		run_r("workbook", spec, path)
	return Path(path)


def materialize_rds(path):
	path = Path(path).resolve()
	stat = path.stat()
	key = hashlib.sha256(f"{path}:{stat.st_size}:{stat.st_mtime_ns}".encode()).hexdigest()
	directory = Path("/tmp/analysis-results") / key
	directory.mkdir(parents = True, exist_ok = True, mode = 0o700)
	with (directory / "lock").open("a") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		marker = directory / "table-name"
		if marker.is_file():
			table = directory / marker.read_text()
			if table.is_file():
				return table
		table = Path(run_r("export-table", path, directory))
		if table.parent != directory or not table.is_file():
			raise ValueError("RDS exchange escaped its temporary directory")
		marker.write_text(table.name)
	return table


# 🚩 Exact aggregate workbook sources
def workbook_sources(path):
	with zipfile.ZipFile(path) as archive:
		if "results/manifest.json" not in archive.namelist():
			return {}
		manifest = json.loads(archive.read("results/manifest.json"))
		if manifest.get("format") != "analysis-tables-v1":
			raise ValueError(f"Unsupported result workbook: {path}")
		return manifest.get("files", {})


def resolve_table(path):
	path = Path(path)
	if path.is_file():
		return materialize_rds(path) if path.suffix.lower() == ".rds" else path
	stem = path.name.removesuffix(".gz").removesuffix(".tsv").removesuffix(".csv")
	private = path.with_name(stem + ".rds")
	if private.is_file():
		return materialize_rds(private)
	for workbook in sorted(path.parent.glob("*.xlsx")):
		entry = workbook_sources(workbook).get(path.name)
		if entry is None:
			continue
		part = f"results/exports/{entry['sha256']}.bin"
		if entry["part"] != part or len(entry["sha256"]) != 64:
			raise ValueError("Invalid workbook export path")
		directory = Path("/tmp/analysis-results") / entry["sha256"]
		directory.mkdir(parents = True, exist_ok = True, mode = 0o700)
		target = directory / path.name
		if not target.is_file():
			with zipfile.ZipFile(workbook) as archive:
				data = archive.read(part)
			if hashlib.sha256(data).hexdigest() != entry["sha256"]:
				raise ValueError(f"Workbook export checksum mismatch: {workbook}")
			with tempfile.NamedTemporaryFile(dir = directory, delete = False) as stream:
				stream.write(data)
				name = stream.name
			os.replace(name, target)
		return target
	return path


def read_table(path, **options):
	import pandas as pd

	options.setdefault("sep", "\t")
	options.setdefault("float_precision", "round_trip")
	return pd.read_csv(resolve_table(path), **options)
