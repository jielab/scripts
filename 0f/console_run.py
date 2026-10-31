#!/usr/bin/env python3


# 🚩 Console runner
"""Run a CLI in the foreground with full file logs and event-driven console output (no periodic time reports)."""

import argparse
import csv
from contextlib import contextmanager
from collections import deque
import datetime
import os
from pathlib import Path
import re
import selectors
import signal
import subprocess
import sys
import time
from process_lifetime import arm_parent_death, cleanup as cleanup_descendants, prctl


def descendants(parent):
	entries = {}
	for path in Path("/proc").glob("[0-9]*/stat"):
		try:
			fields = path.read_text().rsplit(")", 1)[1].split()
			entries[int(path.parent.name)] = (int(fields[1]), fields[19])
		except (OSError, ValueError, IndexError):
			pass
	selected = {parent}
	while True:
		new = {pid for pid, (ppid, _) in entries.items() if ppid in selected}
		if new <= selected:
			break
		selected |= new
	return {pid: entries[pid][1] for pid in selected if pid in entries}


def send_remaining(pids, sig):
	for pid, start in pids.items():
		try:
			fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
			if fields[19] == start and fields[0] != "Z":
				os.kill(pid, sig)
		except (OSError, IndexError):
			pass


@contextmanager
def managed_process(command, env, capture=True):
	# Logging mode never controls workload lifetime. Even verbose/help entry points
	# get the same terminal-bound supervisor, including parent SIGKILL protection.
	guard = [sys.executable, str(Path(__file__).with_name("process_lifetime.py")),
	         "--parent-pid", str(os.getpid()), "--", *command]
	# If the inner supervisor itself dies, adopt and clean its workers here.
	prctl(36, 1)
	proc = subprocess.Popen(guard, env=env, stdout=subprocess.PIPE if capture else None,
	                        stderr=subprocess.STDOUT if capture else None, start_new_session=False)
	try:
		yield proc
	finally:
		cleanup_descendants(proc)
		if proc.stdout is not None:
			proc.stdout.close()


def passthrough(script, args):
	"""Preserve raw stdout/stdin for verbose mode, viewers and machine-readable CLIs."""
	background = explicit_background()
	env = dict(os.environ, SCRIPT_CONSOLE_ACTIVE="1", SCRIPT_FOREGROUND="1",
	           SCRIPT_BACKGROUND="1" if background else "0")
	with managed_process(["bash", str(Path(script).resolve()), *args], env, capture=False) as proc:
		def cancel(sig, _frame):
			if proc.poll() is None:
				proc.send_signal(sig)
		for sig in ((signal.SIGINT, signal.SIGTERM) if background else (signal.SIGINT, signal.SIGHUP, signal.SIGTERM)):
			signal.signal(sig, cancel)
		if not background:
			arm_parent_death()
		rc = proc.wait()
		return rc if rc >= 0 else 128 - rc


def explicit_background():
	# nohup communicates its explicit user intent through inherited SIGHUP=IGN.
	return os.getenv("SCRIPT_BACKGROUND") == "1" or signal.getsignal(signal.SIGHUP) == signal.SIG_IGN


def console_label(script, args):
	label = script.stem
	if label != "gwas_format":
		return label
	# GWAS modules may follow options. Skip option values so names such as
	# --label magma or comma-separated --gwas values cannot become the module.
	modules = {"format", "thin", "magma", "liftover", "cis", "lead", "mplot", "h2", "pgs", "all"}
	index = 0
	while index < len(args):
		arg = args[index]
		index += 1
		if arg.startswith("-"):
			# Every GWAS option takes a value except --fill-n, whose value is
			# optional when the next argument is another long option.
			if index < len(args) and not (arg == "--fill-n" and args[index].startswith("--")):
				index += 1
			continue
		requested = arg.lower().replace("small", "format").split(",")
		if all(module in modules for module in requested):
			return f"{label}.{'_'.join(requested)}"
	return f"{label}.all"


def locus_message(event, unit, line):
	"""Describe a finished locus from its saved results, without guessing counts."""
	match = re.search(r"\b(?:log|detail)=(.*)$", line)
	if not match:
		return f"{event} {unit}"
	path = Path(match.group(1))

	def rows(name):
		with (path.parent / name).open() as f:
			return list(csv.DictReader(f, delimiter="\t"))

	if event == "FAIL":
		try:
			for row in rows("final/gwas_loci.tsv"):
				if row.get("status") == "tree_failed":
					return f"FAIL {unit}: {row['reason']}"
		except (OSError, KeyError):
			pass
		try:
			errors = [
				s.strip()
				for s in path.with_suffix(".log").read_text(errors="replace").splitlines()
				if re.search(r"(?i)\b(?:error|fatal)\b|\b\w+(?:Error|Exception):|plot export failed", s)
			]
			if errors:
				return f"FAIL {unit}: {errors[-1][:240]}"
		except OSError:
			pass
		code = re.search(r"\bexit=(\d+)", line)
		return f"FAIL {unit}: exit {code.group(1) if code else 'unknown'}"
	try:
		skipped = rows("final/skipped_loci.tsv")
		if skipped:
			return f"SKIP {unit}: {skipped[0].get('reason') or skipped[0].get('status', 'not evaluable')}"
		sites = rows("loci/sites.tsv")
		target = rows("final/haplotypes.tsv")
		archaic = rows("loci/archaic.tsv")
		return (
			f"DONE {unit}, {len(sites)} SNPs in haplotype, "
			f"{len(target)} and {len(archaic)} haplotypes in target and archaic reference"
		)
	except OSError:
		return f"DONE {unit}"


class Le8Console:
	"""Keep full diagnostics in the log and emit bounded stage summaries."""
	def __init__(self):
		self.failure = ""
		self.pending_r_error = ""
		self.warnings = False
		self.partial = False
		self.workspace = ""
		self.last_stage = ""

	def feed(self, line):
		if re.match(r"^\[LE8\] PARTIAL\b", line):
			self.partial = True
		if "Diagnostic workspace: " in line:
			self.workspace = line.split("Diagnostic workspace: ", 1)[1]
		if re.search(r"(?i)\bwarning\b|^警告", line):
			self.warnings = True
		if self.pending_r_error and line and not re.match(r"^(Calls:|In addition|Execution halted)", line):
			self.failure = line[:320]
			self.pending_r_error = ""
		if re.match(r"^\[(LE8|ABM)\] (FAIL|ERROR)\b", line):
			if not self.failure:
				self.failure = line[:320]
			return None
		if re.match(r"^Error in\b", line) and not self.failure:
			self.pending_r_error = line[:320]
		elif re.match(r"(?i)^(error|fatal)(?:\b|:)|^\w+(?:Error|Exception):", line) and not self.failure:
			if not line.startswith("ERROR: Command ["):
				self.failure = line[:320]
		# Package warnings, traceback continuations, commands and resource JSON
		# are never forwarded. The enclosing logger retains their exact bytes.
		major = re.match(r"^\[LE8(?: twin| GPU| comparison)?\] (START|DONE|SKIP|PROGRESS|PARTIAL)\b", line)
		reuse = re.match(r"^\[(ABM|C1 selective)\] SKIP\b", line)
		abm = re.match(r"^\[(ABM|C1 selective)\] (START|DONE) (S7_crossfit|S7_final_model|selective_crossfit|final_fit|input|project)\b", line)
		info = re.match(r"^\[LE8\] (Module transaction:|Waiting |.*; waiting for |Shiny is already|Share package:)", line)
		if (major or reuse or abm or info) and line != self.last_stage:
			self.last_stage = line
			return line[:320]
		return None


class GridConsole(Le8Console):
	"""GRID stages and work milestones; raw commands/JSON stay in the full log."""
	def __init__(self):
		super().__init__()
		self.milestones = {}

	def feed(self, line):
		visible = super().feed(line)
		if visible:
			return visible
		if line.startswith("[GRID ABM] ") or re.match(r"^\[GRID\] (?:MODULE|PRSformer2?|SKIP)\b", line):
			return line[:320]
		if line.startswith("[GRID] Running scientific self-check"):
			return line
		if line.startswith("RUN "):
			command = re.search(r"/prsformer\.py (shared-roster|runtime|prepare|train|predict|publish)\b", line)
			if not command:
				return None
			line = f"PRSformer {command.group(1)}"
		elif line.startswith("DEPENDENCIES OK:"):
			fields = re.findall(r"\b(?:torch|device|attention)=[^;]+", line)
			if not fields:
				return None
			line = "PRSformer " + "; ".join(fields)
		elif line.startswith('"status": "RUNTIME OK:'):
			line = "PRSformer runtime: official model forward/backward passed"
		elif line.startswith('"gpu": "'):
			line = "GPU: " + line.split('"', 3)[3].removesuffix(",").removesuffix('"')
		else:
			training = re.match(r"^TRAIN epoch (\d+/\d+): ([\d,]+)/([\d,]+) samples;", line)
			transpose = re.match(r"^TRANSPOSE: ([\d,]+)/([\d,]+) rows$", line)
			if training or transpose:
				key, current, total = training.groups() if training else ("transpose", *transpose.groups())
				bucket = 4 * int(current.replace(",", "")) // max(1, int(total.replace(",", "")))
				if self.milestones.get(key) == bucket:
					return None
				self.milestones[key] = bucket
			elif not re.match(
				r"^(MODULE |(?:START|DONE|SKIP|FAIL) GRID |SHARED ROSTER OK:|"
				r"PREPARE(?: OK)?:|GENOTYPES chr\S+: cumulative |TRANSPOSE: writing |"
				r"REUSE verified prepared cache:|VALIDATION epoch |Epoch \d+:|"
				r"Completed:|Scored \d+ |PUBLISH OK:|CHECK(?: OK| complete|:)|"
				r"WAIT for exclusive |Frozen comparison |Evaluated |DONE: |FINAL SUMMARY:|"
				r"GRID Shiny:|Listening on |\[GRID twin\] (?:FIT|COMPLETE)\b|"
				r"\[GRID (?:prepare|preflight)\] |"
				r"\[GRID\] (?:Existing PRSformer outputs reused|Updating the final evidence index)|"
				r"\[GRID\] (?:START|DONE|SKIP|FAIL|SUMMARY)\b)", line
			):
				return None
		if line == self.last_stage:
			return None
		self.last_stage = line
		return (line if line.startswith("[GRID") else f"[GRID] {line}")[:320]


class PrsformerConsole(Le8Console):
	"""Compact preparation and independent-trait training progress."""
	def __init__(self):
		super().__init__()
		self.milestones = {}
		self.trait = None
		self.last_visible = ""

	def feed(self, line):
		super().feed(line)
		trait = re.match(r"^\[GRID\] (PRSformer2?) trait (\d+/\d+): (\S+)", line)
		if trait and trait.groups() != self.trait:
			self.trait = trait.groups()
			self.milestones.clear()
		if line.startswith(("PREPARE:", "TRANSPOSE: writing")):
			self.milestones.pop("transpose", None)
		training = re.match(r"^TRAIN epoch (\d+/\d+): ([\d,]+)/([\d,]+) samples;", line)
		transpose = re.match(r"^TRANSPOSE: ([\d,]+)/([\d,]+) rows$", line)
		if training or transpose:
			key, current, total = training.groups() if training else ("transpose", *transpose.groups())
			current, total = int(current.replace(",", "")), int(total.replace(",", ""))
			if total <= 0:
				return None
			# Training already emits at most once per minute. Keep that heartbeat
			# visible even when a full epoch takes hours. Transpose stays compact.
			if not training:
				bucket = 4 * current // total
				if self.milestones.get(key) == bucket:
					return None
				self.milestones[key] = bucket
			line += f" ({100 * current / total:.1f}%)"
		elif line.startswith("RUN "):
			command = re.search(r"/prsformer\.py (shared-roster|runtime|prepare|train|predict|publish)\b", line)
			if not command:
				return None
			line = f"PRSformer {command.group(1)}"
		elif re.match(r"^\[GRID\] PRSformer2 (train|resume|predict): ", line):
			# PRSformer2 prints a full command; keep only its stage on screen.
			line = line.split(":", 1)[0]
		elif line.startswith("DEPENDENCIES OK:"):
			fields = re.findall(r"\b(?:torch|device|attention)=[^;]+", line)
			if not fields:
				return None
			line = "PRSformer " + "; ".join(fields)
		elif re.match(r'^"(?:device|gpu)": "', line):
			key, value = re.match(r'^"(device|gpu)": "([^"]+)"', line).groups()
			line = ("运行设备: " if key == "device" else "GPU: ") + value
		elif not re.match(
			r"^(COHORT:|SHARED ROSTER OK:|PREPARE(?: OK)?:|"
			r"GENOTYPES chr\S+: cumulative |TRANSPOSE: writing |"
			r"REUSE verified prepared cache:|RESUME:|RUNTIME:|EARLY STOP:|Epoch \d+:|"
			r"Completed:|Scored \d+ |PUBLISH OK:|CHECK OK:|WAIT for exclusive |"
			r"\[GRID\] (?:PRSformer2?\b|SKIP PRSformer2?\b)|"
			r"\[GRID prepare\] |\[3\.prsformer\.sh prepare\] READY:)", line
		):
			return None
		line = re.sub(r"^\[(?:GRID(?: prepare)?|3\.prsformer\.sh prepare)\] ", "", line)
		visible = f"[3.prsformer] {line}"[:320]
		if visible == self.last_visible:
			return None
		self.last_visible = visible
		return visible


class AlleleConsole(Le8Console):
	"""Stage boundaries plus at most one work-progress report per minute."""
	def __init__(self):
		super().__init__()
		self.last_progress = float("-inf")

	def feed(self, line):
		super().feed(line)
		if re.search(r"(?i)\b(?:warning|notice)\b", line) and not self.last_stage == "warning":
			self.last_stage = "warning"
			return "[allele3] 有警告，详见完整日志"
		if line.startswith("[allele3] stage="):
			return line[:320]
		if re.match(r"^\[atlas\] (?:retry |checkpoint=)", line):
			now = time.monotonic()
			if now - self.last_progress >= 60:
				self.last_progress = now
				return line[:320]
		if re.match(r"^\[(?:geometry|compute)\] (?:device=|CUDA unavailable)", line):
			return line[:320]
		if re.match(r"^\[[^]]+\] (?:complete, skip|total=|rows=)", line) or line.startswith("Final outputs:"):
			return line[:320]
		if re.match(r"^\[[^]]+\] [\d,]+(?:/[\d,]+)? .*elapsed=", line):
			now = time.monotonic()
			if now - self.last_progress >= 60 or "(100.0%)" in line:
				self.last_progress = now
				return line[:320]
		return None


class MetagConsole(Le8Console):
	"""Module/dataset boundaries only; per-read/run/download details stay in the log."""
	def __init__(self):
		super().__init__()
		self.warning_messages = set()
		self.reference_audit = False
		self.reference_step = None

	def reference_transition(self, step):
		lines = []
		if self.reference_step and self.reference_step != step:
			lines.append("[metag] DONE reference/" + self.reference_step)
		if step and self.reference_step != step:
			lines.append("[metag] START reference/" + step)
		self.reference_step = step
		return "\n".join(lines) or None

	def finish(self, returncode):
		if self.reference_step:
			step = self.reference_step
			self.reference_step = None
			return f"[metag] {'DONE' if returncode == 0 else 'FAIL'} reference/{step}"
		return None

	def warning(self, message):
		self.warnings = True
		message = message[:280]
		if message in self.warning_messages:
			return None
		self.warning_messages.add(message)
		if len(self.warning_messages) <= 5:
			return "[metag] WARNING " + message
		if len(self.warning_messages) == 6:
			return "[metag] WARNING 另有警告；后续详情保留在完整日志"
		return None

	def feed(self, line):
		super().feed(line)
		if re.match(r"^\[metag\] (START|DONE|SKIP|FAIL|ERROR)\b", line):
			if re.match(r"^\[metag\] (FAIL|ERROR)\b", line):
				self.failure = line[:320]
			return line[:320]
		# Strip timestamps so repeated warnings from download workers collapse.
		diagnostic = re.match(r"^\[(WARNING|WARN|ERROR|CRITICAL)(?:\s+-[^\]]*)?\]:?\s*(.*)", line, re.I)
		if diagnostic:
			level, message = diagnostic.groups()
			if level.upper() in ("ERROR", "CRITICAL"):
				self.failure = message[:320]
				return "[metag] ERROR " + message[:280]
			return self.warning(message)
		if line.startswith("流程停止："):
			self.failure = line[:320]
			return "[metag] FAIL " + line[:280]
		if re.search(r"(?i)\b(?:warning|warn)\b|^警告|^Virus types NOT AUDITED", line):
			return self.warning(line)
		if re.match(r"(?i)^(?:error|fatal)\b|^\w+(?:Error|Exception):", line):
			return "[metag] ERROR " + line[:280]
		if line.startswith("[raw.bad] moved "):
			return self.warning("QC 发现异常文件并移至 raw.bad；文件及原因详见日志和 qc/checks.tsv")
		if line.startswith("[qc/reference]"):
			return "[metag] START reference preparation"
		if line.startswith("[qc/audit]"):
			self.reference_audit = True
			return "[metag] DONE reference preparation; START reference/audit"
		if line == "所选流程全部完成。" and self.reference_audit:
			self.reference_audit = False
			return "[metag] DONE reference/audit"
		# build_db.sh prints each command before executing it under set -e/pipefail.
		# Reaching the next step (or the frozen marker) confirms the prior success.
		# Derive boundaries here so logging edits never invalidate analysis receipts.
		command = re.search(r"(?:^|/)k2_compat\.py (download-taxonomy|download-library)\b", line)
		if command:
			library = re.search(r"--library\s+(\S+)", line)
			step = "taxonomy" if command[1] == "download-taxonomy" else library[1] if library else "library"
			return self.reference_transition(step)
		if line.startswith("kraken2-build --build "):
			return self.reference_transition("build")
		if re.search(r"(?:^|/)metag\.py db-audit\b", line):
			return self.reference_transition("audit")
		if line.startswith("Database frozen."):
			return self.reference_transition(None)
		if line.startswith("SKIP completed step: "):
			completed = self.reference_transition(None)
			return ((completed + "\n") if completed else "") + "[metag] SKIP reference/" + line.split(": ", 1)[1]
		if line.startswith("Retrying step with missing/empty products:"):
			return self.warning(line)
		if line == "[scan] checking saved input inventory":
			return "[metag] START verify saved input inventory"
		if line.startswith("[qc] FASTQ 检查已完成"):
			return "[metag] SKIP FASTQ QC: verified unchanged inputs"
		boundary = re.match(r"^\[scan\] (start|finish) (.+)", line)
		if boundary:
			return f"[metag] {'START' if boundary[1]=='start' else 'DONE'} scan {boundary[2]}"[:320]
		boundary = re.match(r"^\[qc\] prepare ([^:]+): checking file types and pairing$", line)
		if boundary:
			return f"[metag] START input pairing {boundary[1]}"
		if line.startswith("[qc] prepared "):
			return "[metag] DONE input pairing " + line.removeprefix("[qc] prepared ")[:280]
		if line.startswith("[qc] summary "):
			return "[metag] DONE QC " + line.removeprefix("[qc] summary ")[:280]
		start = re.match(r"^\[(qc|split)\] (\d+) runs(?: pending)?; (.*)", line)
		if start:
			return f"[metag] START {start[1]}: {start[2]} runs; {start[3]}"[:320]
		end = re.match(r"^\[(qc|split)\] completed (\d+)/(\d+)", line)
		if end and end[2] == end[3]:
			return f"[metag] DONE {end[1]}: {end[3]} runs"
		return None


class IbdmixConsole(Le8Console):
	"""One shared preparation and one completion per reference; retain full file logs."""
	def feed(self, line):
		# Reuse diagnostic capture, but only display IBDmix milestones below.
		super().feed(line)
		if line.startswith("[GU UKB] PROGRESS "):
			return line
		created = re.match(r"^\[GU CMD\] created=(\d+)", line)
		if created:
			return f"[gu] 染色体任务数={created.group(1)}"
		event = re.match(r"^\[GU CMD\] (START|DONE|SKIP|FAIL) unit=(\S+)", line)
		if event:
			return f"[gu] {event.group(1)} {event.group(2)}"
		stage = re.search(r"\] ((?:PREPARE|REFERENCE DONE|REUSE) unit=.+)$", line)
		if stage:
			return f"[gu] {stage.group(1)}"
		return None


def run(script, args):
	script = Path(script).resolve()
	label = console_label(script, args)
	stage_console = {"le8": Le8Console, "grid": GridConsole, "3.prsformer": PrsformerConsole,
	                 "4.prsformer": PrsformerConsole,
	                 "allele3": AlleleConsole, "metag": MetagConsole}.get(label, lambda: None)()
	if label == "gu" and args[:1] in (["ibdmix"], ["ibdmix-GPU"]):
		stage_console = IbdmixConsole()
	completion_only_gu = label == "gu" and args[:1] in (["phyml"], ["ibdmix"])
	quiet_gu = completion_only_gu
	default_log_dir = Path("/mnt/d/analysis") / script.parent.name / "logs"
	directory = Path(os.environ.get("SCRIPT_LOG_DIR", default_log_dir))
	directory.mkdir(parents=True, exist_ok=True)
	log = directory / f"{label}.{datetime.datetime.now():%Y%m%d-%H%M%S}.{os.getpid()}.log"
	background = explicit_background()
	env = dict(os.environ, SCRIPT_CONSOLE_ACTIVE="1", SCRIPT_FOREGROUND="1",
	           SCRIPT_BACKGROUND="1" if background else "0")
	with managed_process(["bash", str(script), *args], env) as proc:
		cancelled = []

		def cancel(sig, frame):
			if not cancelled:
				cancelled.append((sig, time.monotonic(), descendants(proc.pid)))
				if proc.poll() is None:
					proc.send_signal(sig)

		for sig in ((signal.SIGINT, signal.SIGTERM) if background else (signal.SIGINT, signal.SIGHUP, signal.SIGTERM)):
			signal.signal(sig, cancel)
		if not background:
			arm_parent_death()
		print(f"[{label}] 开始；详细日志：{log}", flush=True)
		done = failed = 0
		stage = ""
		diagnostic = False
		recent = deque(maxlen=8)
		pending = b""
		active_units = set()
		skipped_units = 0
		unit_stages = {}
		total_units = "?"

		def progress(event):
			if quiet_gu:
				print(f"[{label}] {event}", flush=True)
				return
			active = (
				", ".join(
					unit + (":" + unit_stages[unit] if unit in unit_stages else "") for unit in sorted(active_units)
				)
				or "无"
			)
			print(
				f"[{label}] {event} | 完成={done}/{total_units}（复用={skipped_units}）"
				f" 失败={failed} | 运行({len(active_units)})：{active}",
				flush=True,
			)

		def console_text(text):
			if completion_only_gu:
				return text.removeprefix("[GU CMD] ")
			return text

		def show_line(raw):
			nonlocal done, failed, stage, diagnostic, skipped_units, total_units
			original = raw.decode(errors="replace")
			line = original.strip()
			if not line:
				return
			if stage_console is not None:
				visible = stage_console.feed(line)
				if visible:
					print(visible, flush=True)
				return
			if completion_only_gu:
				created = re.match(r"^\[GU CMD\] created=(\d+)", line)
				if created:
					total_units = created.group(1)
					print(f"[{label}] 任务总数={total_units}", flush=True)
					return
				if re.match(r"^\[GU CMD\] (CHECK|RESUME)\b", line):
					if line.startswith("[GU CMD] RESUME"):
						counts = dict(re.findall(r"(\w+)=(\d+)", line))
						skip = int(counts.get("skipped", 0))
						reuse = int(counts.get("reused", 0)) - skip
						print(f"[{label}] Skip {skip}, reuse {reuse}", flush=True)
					return
				detail = re.search(r"\b(GENOTYPE|CALL|REUSE) unit=C(\d+|X) ref=(\S+)(.*)", line)
				if detail:
					action, chrom, ref, rest = detail.groups()
					unit_stages["chr" + chrom] = (
						ref + ":" + ("基因型" if action == "GENOTYPE" else "复用" if action == "REUSE" else "群体分析")
					)
					return
				tree = re.match(r"^\[GU PHYML\] tree=(\S+) (.*?) input=(.*)$", line)
				if tree:
					action, reason, path = tree.groups()
					unit = Path(path).parent.parent.name
					unit_stages[unit] = "建树" if action == "REPLACE" else "复用树/整理输出"
					if action == "REPLACE" and not quiet_gu:
						print(f"[{label}] 建树 {unit} {reason}", flush=True)
					return
			unit_event = re.match(r"^\[GU CMD\] (START|DONE|SKIP|FAIL) unit=(\S+)", line)
			if completion_only_gu and unit_event:
				event, unit = unit_event.groups()
				if event == "START":
					active_units.add(unit)
				else:
					active_units.discard(unit)
					unit_stages.pop(unit, None)
				if event in ("DONE", "SKIP"):
					done += 1
					if event == "SKIP":
						skipped_units += 1
				elif event == "FAIL":
					failed += 1
				recent.append(line)
				if quiet_gu:
					if event in ("FAIL", "DONE"):
						progress(locus_message(event, unit, line))
				elif event != "SKIP" or skipped_units % 50 == 0:
					progress(f"{event} {unit}" if event != "SKIP" else "复用检查")
				if event == "FAIL" and not quiet_gu:
					print(console_text(line), flush=True)
				diagnostic = False
				return
			recent.append(line)
			if re.search(r"\[GU CMD\] (DONE|SKIP) unit=", line):
				done += 1
			if "[GU CMD] FAIL unit=" in line:
				failed += 1
			# Only real stage transitions, never a wall-clock heartbeat.
			major = re.match(
				r"^(JOB=|MODE=|Completed\s|Done\b|FINAL:|"
				r"LE8 .*preflight|\[LE8\] PROGRESS\b|MAHA pipeline:|\[MAHA(?: FINAL)?\]|"
				r"\[[^\]]+\] (?:START|DONE|FAIL|ERROR)\b|"
				r"(?:Starting|Finished|Running) (?:module|step|stage)\b)",
				line,
			)
			# Border emits bounded milestones; retain per-exposure and mediation
			# progress in the console without exposing package/model chatter.
			if label == "border":
				major = major or re.match(
					r"^\[(?:mediation|PheWAS|PheWAS plot)\] |"
					r"^\[border\] PROGRESS\b|"
					r"^(?:Follow-up end:|R packages OK|Input .*|Final figures:)", line,
				)
			# LE8 console: stage boundaries and completed-work progress, plus diagnostics.
			if label == "metag":
				major = major or re.match(r"^\[(?:scan|qc(?:/[^\]]+)?|split|raw\.bad|split\.bad)\]", line)
			if label == "le8":
				major = bool(re.match(r"^\[LE8\] (START|DONE|FAIL|ERROR|PROGRESS)\b", line))
			error = re.search(
				r"(?i)(?:^|[\] :])(?:error|fatal|warning|traceback|exception)\b|"
				r"\b[A-Za-z]+(?:Error|Exception):|Execution halted|^Calls:|^停止执行",
				line,
			)
			continuation = diagnostic and (
				original[:1].isspace() or re.match(r"^(In |[0-9]+:|Calls:|During handling|The above exception)", line)
			)
			if error or continuation:
				if not quiet_gu or not active_units:
					print(console_text(original), flush=True)
				diagnostic = True
			elif (
				major
				and not (label == "gu" and args[:1] == ["phyml"])
				and not (quiet_gu and re.match(r"^\[\d{4}-\d{2}-\d{2} [^\]]+\] (?:START|DONE)\b", line))
			):
				if line != stage:
					print(console_text(line), flush=True)
					stage = line
				diagnostic = False
			else:
				diagnostic = False

		selector = selectors.DefaultSelector()
		selector.register(proc.stdout, selectors.EVENT_READ)
		with log.open("wb") as out:
			while selector.get_map() or proc.poll() is None:
				if proc.poll() is not None:
					# An adopted worker can still hold the log pipe open after the
					# inner supervisor dies. Clean it before waiting for pipe EOF.
					cleanup_descendants(proc)
				for key, _ in selector.select(timeout=1):
					chunk = os.read(key.fileobj.fileno(), 65536)
					if not chunk:
						selector.unregister(key.fileobj)
						break
					out.write(chunk)
					out.flush()
					pending += chunk
					lines = pending.split(b"\n")
					pending = lines.pop()
					if len(pending) > 65536:
						lines.append(pending)
						pending = b""
					for raw in lines:
						show_line(raw)
				now = time.monotonic()
				if cancelled and now - cancelled[0][1] >= 10:
					send_remaining(cancelled[0][2], signal.SIGKILL)
			if pending:
				show_line(pending)
			rc = proc.wait()
		proc.stdout.close()
		selector.close()
		if cancelled:
			# Cover workers that closed stdout before the leader exited.
			send_remaining(cancelled[0][2], signal.SIGKILL)
			rc = 128 + cancelled[0][0]
		if stage_console is not None and hasattr(stage_console, "finish"):
			final_stage = stage_console.finish(rc)
			if final_stage:
				print(final_stage, flush=True)
		state = "已停止" if cancelled else "完成" if rc == 0 else f"失败（退出码 {rc}）"
		if rc == 0 and stage_console is not None and stage_console.partial:
			state = "完成（部分分析未完成，见 SKIP/PARTIAL）"
		counts = f"；完成 {done}，失败 {failed}" if (done or failed) and not quiet_gu else ""
		detail = ""
		if stage_console is not None:
			if rc and not cancelled:
				reason = stage_console.failure or stage_console.pending_r_error
				if reason:
					detail += f"；原因：{reason}"
			if stage_console.warnings:
				detail += "；警告详见日志"
			if stage_console.workspace:
				detail += f"；诊断目录：{stage_console.workspace}"
		print(f"[{label}] {state}{counts}{detail}；日志：{log}", flush=True)
		if rc and not cancelled and not quiet_gu and stage_console is None:
			print("\n".join(console_text(line) for line in recent), file=sys.stderr, flush=True)
		return rc if rc >= 0 else 128 - rc


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--script", required=True)
	parser.add_argument("--passthrough", action="store_true")
	parser.add_argument("args", nargs=argparse.REMAINDER)
	a = parser.parse_args()
	runner = passthrough if a.passthrough else run
	raise SystemExit(runner(a.script, a.args[1:] if a.args[:1] == ["--"] else a.args))
