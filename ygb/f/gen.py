#!/usr/bin/env python3
"""YGB genotype liftover and streamed PLINK exports."""

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path


# 🚩 lift_bim
"""Lift BIM SNP coordinates without changing BED row or allele encoding."""
import collections
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def lift_bim():
	source, outdir, chain, fasta = map(Path, sys.argv[1:])
	source = source.resolve()
	outdir.mkdir(parents=True, exist_ok=True)
	output = outdir / source.name
	if output.resolve() == source:
		raise ValueError("Input and output must differ")
	reference_indexes = [Path(str(fasta) + ".fai")]
	if str(fasta).endswith((".gz", ".bgz", ".bgzf")):
		reference_indexes.append(Path(str(fasta) + ".gzi"))
	for path in (source, chain, fasta, *reference_indexes):
		if not path.is_file():
			raise FileNotFoundError(path)
	companions = [source.with_suffix(ext) for ext in (".bed", ".fam")]
	if source.with_suffix(".nosex").exists():
		companions.append(source.with_suffix(".nosex"))
	for path in companions:
		if not path.exists():
			raise FileNotFoundError(path)
		dst = outdir / path.name
		if dst.exists() and not dst.is_symlink():
			raise FileExistsError(dst)
	rows = [line.split() for line in source.read_text().splitlines()]
	if not rows or any(len(row) != 6 for row in rows):
		raise ValueError("Expected nonempty six-column BIM")
	lift = os.environ.get("LIFTOVER_BIN") or shutil.which("liftOver") or configured_path('DIR_SOFTWARE', 'bin/liftOver', default_mode="local")
	with tempfile.TemporaryDirectory(prefix=".liftGen.bim.", dir=outdir) as tmp:
		tmp = Path(tmp)
		bed, mapped, unmapped = [tmp / name for name in ("source.bed", "mapped.bed", "unmapped.bed")]
		expected = {}
		with bed.open("w") as handle:
			for index, row in enumerate(rows):
				chrom = row[0].removeprefix("chr")
				chrom = {"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"}.get(chrom, chrom)
				chrom = "chr" + chrom
				# Non-SNP alleles need normalization; do not guess their encoding.
				if int(row[3]) > 0 and all(a.upper() in ("A", "C", "G", "T") for a in row[4:6]):
					pos = int(row[3])
					expected[index] = chrom
					handle.write(f"{chrom}\t{pos - 1}\t{pos}\t{index}\t0\t+\n")
		subprocess.run([lift, str(bed), str(chain), str(mapped), str(unmapped)], check=True)
		groups = collections.defaultdict(list)
		for line in mapped.read_text().splitlines():
			fields = line.split()
			groups[int(fields[3])].append(fields)
		positions = {}
		candidates = tmp / "candidate.bed"
		with candidates.open("w") as handle:
			for index, records in groups.items():
				if len(records) != 1:
					continue
				chrom, start, end, _, _, strand = records[0]
				if chrom == expected[index] and strand == "+" and int(end) - int(start) == 1:
					positions[index] = int(start) + 1
					handle.write(f"{chrom}\t{start}\t{end}\t{index}\n")
		verified = {}
		if positions:
			result = subprocess.run(
				["bedtools", "getfasta", "-fi", str(fasta), "-bed", str(candidates), "-nameOnly", "-tab"],
				check=True,
				capture_output=True,
				text=True,
			)
			for line in result.stdout.splitlines():
				index, allele = line.split("\t")
				index = int(index)
				if allele.upper() in [a.upper() for a in rows[index][4:6]]:
					verified[index] = positions[index]
		staged = tmp / source.name
		with staged.open("w") as handle:
			for index, row in enumerate(rows):
				row[3] = str(verified.get(index, -1))
				handle.write("\t".join(row) + "\n")
		os.replace(staged, output)
		for path in companions:
			dst = outdir / path.name
			link = tmp / path.name
			link.symlink_to(os.path.relpath(path, outdir.resolve()))
			os.replace(link, dst)
		report = (
			f"input\t{source}\noutput\t{output}\ntotal_variants\t{len(rows)}\n"
			f"unique_forward_same_chr_length\t{len(positions)}\n"
			f"verified_either_allele\t{len(verified)}\nfailed_pos_minus_1\t{len(rows) - len(verified)}\n"
		)
		log = tmp / "complete.log"
		log.write_text(report)
		os.replace(log, str(output) + ".liftGen.log")
		print(report, end="")


# 🚩 export_cli
"""Stream PLINK's sample-major .raw output into gzip without a plain disk copy."""
import argparse
import gzip
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile


def rename_raw_header(header, alias_file):
	aliases = {}
	for line in Path(alias_file).read_text().splitlines():
		row = line.split()
		if not row or row[0].startswith("#") or row[0].upper() in ("SNP", "RSID", "VARIANT", "ID"):
			continue
		alias = row[1] if len(row) > 1 else row[0]
		if row[0] in aliases or alias in aliases.values():
			raise ValueError("Duplicate SNP or alias in " + str(alias_file))
		aliases[row[0]] = alias
	columns = header.decode().split()
	for i, name in enumerate(columns):
		if name in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"):
			continue
		snp, sep, allele = name.rpartition("_")
		if not sep or snp not in aliases:
			raise ValueError("Exported SNP absent from alias map: " + name)
		columns[i] = aliases[snp] + "_" + allele
	if len(set(columns)) != len(columns):
		raise ValueError("Duplicate exported column names")
	return ("\t".join(columns) + "\n").encode()


def export_plink_raw(output, command, alias_file=None):
	output = Path(output)
	output.parent.mkdir(parents=True, exist_ok=True)
	# Keep the FIFO on Linux: Windows-mounted output drives may not support it.
	# The compressed staging file stays beside the output for atomic publication.
	with (
		tempfile.TemporaryDirectory(prefix=".raw-pipe.", dir="/tmp") as pipe_directory,
		tempfile.TemporaryDirectory(prefix=".raw-export.", dir=output.parent) as temporary,
	):
		prefix = Path(pipe_directory) / "data"
		fifo = Path(str(prefix) + ".raw")
		os.mkfifo(fifo)
		# Hold the pipe open until PLINK exits, including failure before its first write.
		keeper = os.open(fifo, os.O_RDWR)
		reader = fifo.open("rb")
		packed = Path(temporary) / "raw.gz"
		gzip_proc = plink_proc = None
		try:
			with packed.open("wb") as dest:
				gzip_proc = subprocess.Popen(["gzip", "-1c"], stdin=reader, stdout=dest)
				plink_proc = subprocess.Popen([*command, "--out", str(prefix)])
				while plink_proc.poll() is None:
					if gzip_proc.poll() is not None:
						raise RuntimeError("gzip exited before PLINK completed")
					try:
						plink_proc.wait(timeout=0.2)
					except subprocess.TimeoutExpired:
						pass
				rc = plink_proc.returncode
				os.close(keeper)
				keeper = None
				if rc:
					raise subprocess.CalledProcessError(rc, command)
				if gzip_proc.wait():
					raise RuntimeError("gzip failed while exporting PLINK .raw")
			with gzip.open(packed, "rb") as check:
				if not check.read(1):
					raise RuntimeError("Empty compressed PLINK output")
			if alias_file:
				renamed = Path(temporary) / "named.raw.gz"
				with gzip.open(packed, "rb") as src, gzip.open(renamed, "wb", compresslevel=1) as dst:
					dst.write(rename_raw_header(src.readline(), alias_file))
					shutil.copyfileobj(src, dst)
				packed = renamed
			packed.replace(output)
		finally:
			if keeper is not None:
				os.close(keeper)
			reader.close()
			for proc in (plink_proc, gzip_proc):
				if proc is not None and proc.poll() is None:
					proc.terminate()
					try:
						proc.wait(timeout=10)
					except subprocess.TimeoutExpired:
						proc.kill()
						proc.wait()
			log = Path(str(prefix) + ".log")
			if log.exists():
				shutil.copyfile(log, str(output) + ".log")


def export_cli():
	def stop(signum, frame):
		raise SystemExit(128 + signum)

	for sig in (signal.SIGHUP, signal.SIGTERM):
		signal.signal(sig, stop)
	p = argparse.ArgumentParser()
	p.add_argument("--output", required=True)
	p.add_argument("--alias-file")
	p.add_argument("command", nargs=argparse.REMAINDER)
	a = p.parse_args()
	command = a.command[1:] if a.command[:1] == ["--"] else a.command
	if not command or "--out" in command:
		p.error("supply PLINK command without --out")
	export_plink_raw(a.output, command, a.alias_file)


def write_apoe(vcf_path, output_path):
	"""Read the two PLINK-exported variants in their actual REF/ALT order."""
	haplotypes = {"TC": "e3", "TT": "e2", "CC": "e4", "CT": "e1"}
	samples = []
	calls = {}
	with gzip.open(vcf_path, "rt") as vcf:
		for line in vcf:
			if line.startswith("##"):
				continue
			row = line.rstrip("\r\n").split("\t")
			if row[0] == "#CHROM":
				samples = row[9:]
				if not samples or len(set(samples)) != len(samples):
					raise ValueError("Missing or duplicate APOE sample IDs")
				continue
			chrom, _, variant, ref, alt = row[:5]
			if (chrom.removeprefix("chr") != "19" or variant not in {"rs429358", "rs7412"}
					or variant in calls or {ref, alt} != {"C", "T"} or len(row) != 9 + len(samples)):
				raise ValueError("Unexpected APOE variant or sample count: " + variant)
			gt_index = row[8].split(":").index("GT")
			alleles = str.maketrans({"0": ref, "1": alt})
			calls[variant] = [v.split(":")[gt_index].translate(alleles) for v in row[9:]]
	if not samples or set(calls) != {"rs429358", "rs7412"}:
		raise ValueError("APOE requires both rs429358 and rs7412")
	with open(output_path, "w") as out:
		print("IID apoe", file=out)
		for iid, a, b in zip(samples, calls["rs429358"], calls["rs7412"]):
			result = "NA"
			if (len(a) == len(b) == 3 and a[1] in "|/" and b[1] in "|/"
					and all(base in "CT" for base in (a[0], a[2], b[0], b[2]))):
				if a[0] == a[2] or b[0] == b[2] or a[1] == b[1] == "|":
					pair = [haplotypes[a[0] + b[0]], haplotypes[a[2] + b[2]]]
					result = "".join(sorted(pair))
			print(iid, result, file=out)


# 🚩 Command dispatch
if __name__ == "__main__":
	if len(sys.argv) < 2 or sys.argv[1] in ("--help", "-h"):
		print("Usage: gen.py lift-bim BIM OUTPUT_DIR CHAIN FASTA | apoe VCF.GZ OUTPUT | export --output FILE -- PLINK_COMMAND")
	else:
		operation = sys.argv.pop(1)
		if operation == "lift-bim":
			lift_bim()
		elif operation == "apoe":
			write_apoe(*sys.argv[1:])
		elif operation == "export":
			export_cli()
		else:
			raise SystemExit("Unknown genotype operation: " + operation)
