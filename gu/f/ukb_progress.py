#!/usr/bin/env python3
"""Report Linux file-read progress without changing UKB preparation receipts.

Run a command after --, or inspect an existing process with --pid PID.
Read-position ETA describes only the current file pass, not pipeline completion.
"""

import argparse
import os
from pathlib import Path
import signal
import subprocess
import time


def process_tree(pid):
	"""Include children of every thread (some native tools spawn from workers)."""
	pending, found = [pid], []
	while pending:
		current = pending.pop()
		if current in found:
			continue
		found.append(current)
		for task in Path(f"/proc/{current}/task").glob("*/children"):
			try:
				pending.extend(map(int, task.read_text().split()))
			except (OSError, ValueError):
				pass
	return found


def input_position(pid):
	files = {}
	for fd in Path(f"/proc/{pid}/fd").glob("*"):
		try:
			path = str(fd.readlink())
			if not path.endswith((".pgen", ".pvar", ".pvar.zst", ".vcf.gz", ".vcf", ".bcf")):
				continue
			info = dict(line.split(":", 1) for line in Path(f"/proc/{pid}/fdinfo/{fd.name}").read_text().splitlines())
			if int(info["flags"], 8) & os.O_ACCMODE != os.O_RDONLY:
				continue
			size, position = fd.stat().st_size, int(info["pos"])
			if size > 0:
				files[path] = (size, max(position, files.get(path, (0, 0))[1]))
		except (OSError, ValueError, KeyError):
			pass
	if not files:
		return None
	path = max(files, key=lambda x: files[x][0])
	return path, *files[path]


class Progress:
	def __init__(self, pid):
		self.pid = pid
		self.started = time.monotonic()
		self.previous = None
		self.chrom = "?"
		try:
			args = Path(f"/proc/{pid}/cmdline").read_bytes().decode().split("\0")
			if "--chr" in args:
				self.chrom = args[args.index("--chr") + 1].replace(" ", ",")
		except (OSError, UnicodeError, IndexError):
			pass

	def report(self):
		now = time.monotonic()
		for pid in reversed(process_tree(self.pid)):
			try:
				args = Path(f"/proc/{pid}/cmdline").read_bytes().decode().strip("\0").split("\0")
				if not args[0]:
					continue
				stage = Path(args[0]).name
				if "--out" in args:
					stage += ":" + Path(args[args.index("--out") + 1]).name
				elif stage == "bcftools" and len(args) > 1:
					stage += ":" + args[1]
				position = input_position(pid)
				message = f"[GU UKB] PROGRESS chr={self.chrom} stage={stage} pid={pid} observed_min={(now-self.started)/60:.1f}"
				if position:
					path, size, offset = position
					message += f" input={Path(path).name} read_position={offset/1e9:.2f}/{size/1e9:.2f}GB ({100*offset/size:.1f}%)"
					key = (pid, path)
					if self.previous and self.previous[0] == key:
						_, previous_time, previous_offset = self.previous
						rate = (offset-previous_offset) / (now-previous_time)
						if rate > 0 and offset < size:
							message += f" read_MB_s={rate/1e6:.1f} current_read_eta_min={(size-offset)/rate/60:.1f}"
					self.previous = (key, now, offset)
				else:
					self.previous = None
				print(message, flush=True)
				return
			except (OSError, UnicodeError, IndexError):
				continue


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--pid", type=int)
	parser.add_argument("--interval", type=float, default=30)
	parser.add_argument("command", nargs=argparse.REMAINDER)
	args = parser.parse_args()
	command = args.command[1:] if args.command[:1] == ["--"] else args.command
	if args.interval <= 0 or bool(args.pid) == bool(command):
		parser.error("use either --pid PID or -- COMMAND, with a positive interval")
	if args.pid:
		monitor = Progress(args.pid)
		# A start-time identity prevents following an unrelated recycled PID.
		stat = Path(f"/proc/{args.pid}/stat")
		try:
			identity = stat.read_text().rsplit(")", 1)[1].split()[19]
			while stat.read_text().rsplit(")", 1)[1].split()[19] == identity:
				monitor.report()
				time.sleep(args.interval)
		except FileNotFoundError:
			pass
		return 0
	process = subprocess.Popen(command)
	monitor = Progress(process.pid)
	try:
		while True:
			try:
				code = process.wait(timeout=args.interval)
				return code if code >= 0 else 128-code
			except subprocess.TimeoutExpired:
				monitor.report()
	finally:
		if process.poll() is None:
			for pid in reversed(process_tree(process.pid)):
				try:
					os.kill(pid, signal.SIGTERM)
				except ProcessLookupError:
					pass
			process.wait()


if __name__ == "__main__":
	raise SystemExit(main())
