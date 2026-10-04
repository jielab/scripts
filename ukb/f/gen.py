#!/usr/bin/env python3
"""UKB genotype liftover and streamed PLINK exports."""


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
	companions = [source.with_suffix(ext) for ext in (".bed", ".fam", ".nosex")]
	for path in companions:
		if not path.exists():
			raise FileNotFoundError(path)
		dst = outdir / path.name
		if dst.exists() and not dst.is_symlink():
			raise FileExistsError(dst)
	rows = [line.split() for line in source.read_text().splitlines()]
	if not rows or any(len(row) != 6 for row in rows):
		raise ValueError("Expected nonempty six-column BIM")
	lift = os.environ.get("LIFTOVER_BIN") or shutil.which("liftOver") or "/mnt/d/software/bin/liftOver"
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


def export_plink_raw(output, command):
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
	p.add_argument("command", nargs=argparse.REMAINDER)
	a = p.parse_args()
	command = a.command[1:] if a.command[:1] == ["--"] else a.command
	if not command or "--out" in command:
		p.error("supply PLINK command without --out")
	export_plink_raw(a.output, command)


# 🚩 Command dispatch
if __name__ == "__main__":
	if len(sys.argv) < 2 or sys.argv[1] in ("--help", "-h"):
		print("Usage: gen.py lift-bim BIM OUTPUT_DIR CHAIN FASTA | export --output FILE -- PLINK_COMMAND")
	else:
		operation = sys.argv.pop(1)
		if operation == "lift-bim":
			lift_bim()
		elif operation == "export":
			export_cli()
		else:
			raise SystemExit("Unknown genotype operation: " + operation)
