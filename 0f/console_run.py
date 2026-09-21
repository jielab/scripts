#!/usr/bin/env python3
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


def descendants(parent):
    entries = {}
    for path in Path('/proc').glob('[0-9]*/stat'):
        try:
            fields = path.read_text().rsplit(')', 1)[1].split()
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
            fields = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
            if fields[19] == start and fields[0] != 'Z':
                os.kill(pid, sig)
        except (OSError, IndexError):
            pass


@contextmanager
def managed_process(command, env):
    proc = subprocess.Popen(command, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, start_new_session=True)
    try:
        yield proc
    finally:
        if proc.poll() is None or sys.exc_info()[0] is not None:
            # A logging/console failure must not orphan the computation.
            pids = descendants(proc.pid)
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            send_remaining(pids, signal.SIGKILL)
            proc.wait()
        proc.stdout.close()


def console_label(script, args):
    label = script.stem
    if label != 'gwas_format':
        return label
    # GWAS modules may follow options. Skip option values so names such as
    # --label magma or comma-separated --gwas values cannot become the module.
    modules = {'format', 'thin', 'magma', 'liftover', 'cis', 'lead', 'mplot', 'h2', 'pgs', 'all'}
    index = 0
    while index < len(args):
        arg = args[index]
        index += 1
        if arg.startswith('-'):
            # Every GWAS option takes a value except --fill-n, whose value is
            # optional when the next argument is another long option.
            if index < len(args) and not (arg == '--fill-n' and args[index].startswith('--')):
                index += 1
            continue
        requested = arg.lower().replace('small', 'format').split(',')
        if all(module in modules for module in requested):
            return f'{label}.{"_".join(requested)}'
    return f'{label}.all'


def locus_message(event, unit, line):
    """Describe a finished locus from its saved results, without guessing counts."""
    match = re.search(r'\b(?:log|detail)=(.*)$', line)
    if not match:
        return f'{event} {unit}'
    path = Path(match.group(1))
    def rows(name):
        with (path.parent / name).open() as f:
            return list(csv.DictReader(f, delimiter='\t'))
    if event == 'FAIL':
        try:
            for row in rows('final/gwas_loci.tsv'):
                if row.get('status') == 'tree_failed':
                    return f'FAIL {unit}: {row["reason"]}'
        except (OSError, KeyError):
            pass
        try:
            errors = [s.strip() for s in path.with_suffix('.log').read_text(errors='replace').splitlines()
                      if re.search(r'(?i)\b(?:error|fatal)\b|\b\w+(?:Error|Exception):|plot export failed', s)]
            if errors:
                return f'FAIL {unit}: {errors[-1][:240]}'
        except OSError:
            pass
        code = re.search(r'\bexit=(\d+)', line)
        return f'FAIL {unit}: exit {code.group(1) if code else "unknown"}'
    try:
        skipped = rows('final/skipped_loci.tsv')
        if skipped:
            return f'SKIP {unit}: {skipped[0].get("reason") or skipped[0].get("status", "not evaluable")}'
        sites = rows('loci/sites.tsv')
        target = rows('final/haplotypes.tsv')
        archaic = rows('loci/archaic.tsv')
        return (f'DONE {unit}, {len(sites)} SNPs in haplotype, '
                f'{len(target)} and {len(archaic)} haplotypes in target and archaic reference')
    except OSError:
        return f'DONE {unit}'


def run(script, args):
    script = Path(script).resolve()
    label = console_label(script, args)
    completion_only_gu = label == 'gu' and args[:1] in (['phyml'], ['ibdmix'])
    quiet_gu = completion_only_gu
    default_log_dir = Path('/mnt/d/analysis') / script.parent.name / 'logs'
    directory = Path(os.environ.get('SCRIPT_LOG_DIR', default_log_dir))
    directory.mkdir(parents=True, exist_ok=True)
    log = directory / f'{label}.{datetime.datetime.now():%Y%m%d-%H%M%S}.{os.getpid()}.log'
    env = dict(os.environ, SCRIPT_CONSOLE_ACTIVE='1')
    with managed_process(['bash', str(script), *args], env) as proc:
        cancelled = []
        def cancel(sig, frame):
            if not cancelled:
                cancelled.append((sig, time.monotonic(), descendants(proc.pid)))
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                send_remaining(cancelled[0][2], signal.SIGTERM)
        for sig in (signal.SIGINT, signal.SIGHUP, signal.SIGTERM):
            signal.signal(sig, cancel)
        print(f'[{label}] 开始；详细日志：{log}', flush=True)
        done = failed = 0
        stage = ""
        diagnostic = False
        recent = deque(maxlen=8)
        pending = b''
        active_units = set()
        skipped_units = 0
        unit_stages = {}
        total_units = '?'
        def progress(event):
            if quiet_gu:
                print(f'[{label}] {event}', flush=True)
                return
            active = ', '.join(unit + (':' + unit_stages[unit] if unit in unit_stages else '')
                               for unit in sorted(active_units)) or '无'
            print(f'[{label}] {event} | 完成={done}/{total_units}（复用={skipped_units}）'
                  f' 失败={failed} | 运行({len(active_units)})：{active}', flush=True)
        def console_text(text):
            if completion_only_gu:
                return text.removeprefix('[GU CMD] ')
            return text

        def show_line(raw):
            nonlocal done, failed, stage, diagnostic, skipped_units, total_units
            original = raw.decode(errors='replace')
            line = original.strip()
            if not line:
                return
            if completion_only_gu:
                created = re.match(r'^\[GU CMD\] created=(\d+)', line)
                if created:
                    total_units = created.group(1)
                    print(f'[{label}] 任务总数={total_units}', flush=True)
                    return
                if re.match(r'^\[GU CMD\] (CHECK|RESUME)\b', line):
                    if line.startswith('[GU CMD] RESUME'):
                        counts = dict(re.findall(r'(\w+)=(\d+)', line))
                        skip = int(counts.get('skipped', 0))
                        reuse = int(counts.get('reused', 0)) - skip
                        print(f'[{label}] Skip {skip}, reuse {reuse}', flush=True)
                    return
                detail = re.search(r'\b(GENOTYPE|CALL|REUSE) unit=C(\d+|X) ref=(\S+)(.*)', line)
                if detail:
                    action, chrom, ref, rest = detail.groups()
                    unit_stages['chr'+chrom] = ref + ':' + ('基因型' if action == 'GENOTYPE' else '复用' if action == 'REUSE' else '群体分析')
                    return
                tree = re.match(r'^\[GU PHYML\] tree=(\S+) (.*?) input=(.*)$', line)
                if tree:
                    action, reason, path = tree.groups()
                    unit = Path(path).parent.parent.name
                    unit_stages[unit] = '建树' if action == 'REPLACE' else '复用树/整理输出'
                    if action == 'REPLACE' and not quiet_gu:
                        print(f'[{label}] 建树 {unit} {reason}', flush=True)
                    return
            unit_event = re.match(r'^\[GU CMD\] (START|DONE|SKIP|FAIL) unit=(\S+)', line)
            if completion_only_gu and unit_event:
                event, unit = unit_event.groups()
                if event == 'START':
                    active_units.add(unit)
                else:
                    active_units.discard(unit)
                    unit_stages.pop(unit, None)
                if event in ('DONE', 'SKIP'):
                    done += 1
                    if event == 'SKIP':
                        skipped_units += 1
                elif event == 'FAIL':
                    failed += 1
                recent.append(line)
                if quiet_gu:
                    if event in ('FAIL', 'DONE'):
                        progress(locus_message(event, unit, line))
                elif event != 'SKIP' or skipped_units % 50 == 0:
                    progress(f'{event} {unit}' if event != 'SKIP' else '复用检查')
                if event == 'FAIL' and not quiet_gu:
                    print(console_text(line), flush=True)
                diagnostic = False
                return
            recent.append(line)
            if re.search(r'\[GU CMD\] (DONE|SKIP) unit=', line):
                done += 1
            if '[GU CMD] FAIL unit=' in line:
                failed += 1
            # Only real stage transitions, never a wall-clock heartbeat.
            major = re.match(
                r'^(JOB=|MODE=|Completed\s|Done\b|FINAL:|'
                r'LE8 .*preflight|MAHA pipeline:|\[MAHA(?: FINAL)?\]|'
                r'\[[^\]]+\] (?:START|DONE|FAIL|ERROR)\b|'
                r'(?:Starting|Finished|Running) (?:module|step|stage)\b)', line)
            # LE8 console: stage boundaries only, plus diagnostics below.
            if label == 'le8':
                major = bool(re.match(r'^\[LE8\] (START|DONE|FAIL|ERROR)\b', line))
            error = re.search(
                r'(?i)(?:^|[\] :])(?:error|fatal|warning|traceback|exception)\b|'
                r'\b[A-Za-z]+(?:Error|Exception):|Execution halted|^Calls:|^停止执行', line)
            continuation = diagnostic and (
                original[:1].isspace() or re.match(r'^(In |[0-9]+:|Calls:|During handling|The above exception)',line))
            if error or continuation:
                if not quiet_gu or not active_units:
                    print(console_text(original), flush=True)
                diagnostic = True
            elif major and not (label == 'gu' and args[:1] == ['phyml']) and not (quiet_gu and re.match(r'^\[\d{4}-\d{2}-\d{2} [^\]]+\] (?:START|DONE)\b', line)):
                if line != stage:
                    print(console_text(line), flush=True)
                    stage = line
                diagnostic = False
            else:
                diagnostic = False
        selector = selectors.DefaultSelector()
        selector.register(proc.stdout, selectors.EVENT_READ)
        with log.open('wb') as out:
            while selector.get_map() or proc.poll() is None:
                for key, _ in selector.select(timeout=1):
                    chunk = os.read(key.fileobj.fileno(), 65536)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        break
                    out.write(chunk)
                    out.flush()
                    pending += chunk
                    lines = pending.split(b'\n')
                    pending = lines.pop()
                    if len(pending) > 65536:
                        lines.append(pending)
                        pending = b''
                    for raw in lines:
                        show_line(raw)
                now = time.monotonic()
                if cancelled and now-cancelled[0][1] >= 10:
                    send_remaining(cancelled[0][2], signal.SIGKILL)
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
            if pending:
                show_line(pending)
            rc = proc.wait()
        proc.stdout.close()
        selector.close()
        if cancelled:
            # Cover workers that closed stdout before the leader exited.
            send_remaining(cancelled[0][2], signal.SIGKILL)
            rc = 128 + cancelled[0][0]
        state = '已停止' if cancelled else '完成' if rc == 0 else f'失败（退出码 {rc}）'
        counts = f'；完成 {done}，失败 {failed}' if (done or failed) and not quiet_gu else ''
        print(f'[{label}] {state}{counts}；日志：{log}', flush=True)
        if rc and not cancelled and not quiet_gu:
            print('\n'.join(console_text(line) for line in recent), file=sys.stderr, flush=True)
        return rc if rc >= 0 else 128-rc


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--script', required=True)
    parser.add_argument('args', nargs=argparse.REMAINDER)
    a = parser.parse_args()
    raise SystemExit(run(a.script, a.args[1:] if a.args[:1] == ['--'] else a.args))
