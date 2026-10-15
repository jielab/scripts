"""Exercise terminal job control without loading analysis inputs or models."""
import json
import os
from pathlib import Path
import pty
import select
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]


def process(pid):
    try:
        raw = Path(f'/proc/{pid}/stat').read_text()
        fields = raw[raw.rfind(')') + 2:].split()
        return dict(state=fields[0], parent=int(fields[1]), group=int(fields[2]),
                    session=int(fields[3]), birth=fields[19])
    except (OSError, IndexError):
        return None


@unittest.skipUnless(sys.platform == 'linux', 'Linux terminal and /proc test')
class ForegroundRuntimeTest(unittest.TestCase):
    def exercise(self, verbose=False, terminate_logger=False):
        with tempfile.TemporaryDirectory(prefix='le8-foreground-test-') as directory:
            work = Path(directory)
            ready, pulse, pool = (work / name for name in ('ready.json', 'pulse', 'pool.json'))
            worker = work / 'worker.py'
            worker.write_text(
                'import json,os,time\nfrom pathlib import Path\n'
                f'Path({str(ready)!r}).write_text(json.dumps(dict(pid=os.getpid())))\n'
                'while True:\n'
                f' Path({str(pulse)!r}).write_text(str(time.monotonic_ns()))\n'
                ' time.sleep(.05)\n'
            )
            executor = work / 'executor.py'
            executor.write_text(
                'import importlib.util,os,sys\n'
                f'spec=importlib.util.spec_from_file_location("resources",{str(ROOT / "f/0.resources.py")!r})\n'
                'mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)\n'
                f'mod.execute([dict(task_id="terminal",argv=[sys.executable,{str(worker)!r}],'
                f'cores=1,memory_gib=.125,outputs=[{str(work / "done")!r}])],{str(work / "output")!r},'
                f'cores=2,memory_gib=1,pool_path={str(pool)!r},reserve_gib=0)\n'
            )
            scope = []
            if shutil.which('systemctl') and shutil.which('systemd-run'):
                state = subprocess.run(['systemctl', '--user', 'is-system-running'],
                                       capture_output=True, text=True, timeout=5)
                if state.stdout.strip() in ('running', 'degraded'):
                    scope = ['systemd-run', '--user', '--scope', '--quiet', '--collect',
                             '-p', 'MemoryMax=256M', '--']
            shim = work / 'python-shim'
            shim.write_text(
                f'#!{sys.executable}\nimport os,sys\n'
                f'assert sys.argv[1] == {str(ROOT / "f/pipeline.py")!r}\n'
                'assert os.environ["SCRIPT_FOREGROUND"] == "1"\n'
                f'os.execv(sys.executable,[sys.executable,{str(ROOT / "f/0.resources.py")!r},'
                '"phase","--label","terminal-test","--cores","1","--memory-gib",".125","--",'
                f'*{scope!r},sys.executable,{str(executor)!r}])\n'
            )
            shim.chmod(0o700)
            env = dict(os.environ, PS1='LE8_TEST_PROMPT> ', PROMPT_COMMAND='',
                       SCRIPT_LOG_DIR=str(work / 'logs'), SCRIPT_VERBOSE=str(int(verbose)),
                       LE8_REPORT_PYTHON=str(shim), LE8_RESOURCE_POOL=str(pool),
                       LE8_TOTAL_CORES='2', LE8_RESOURCE_MEMORY_GIB='1', LE8_RESOURCE_RESERVE_GIB='0')
            for key in ('SCRIPT_CONSOLE_ACTIVE', 'LE8_POOL_CPUS', 'LE8_MEMORY_SCOPE_ACTIVE'):
                env.pop(key, None)
            shell, master = pty.fork()
            if shell == 0:
                os.execve('/bin/bash', ['bash', '--noprofile', '--norc', '-i'], env)
            output = bytearray()
            tracked = {}
            job_group = None

            def pump():
                if select.select([master], [], [], .05)[0]:
                    try:
                        output.extend(os.read(master, 65536))
                    except OSError:
                        pass

            def wait_for(predicate, seconds=12):
                deadline = time.monotonic() + seconds
                while time.monotonic() < deadline:
                    pump()
                    if predicate():
                        return
                self.fail('Timed out waiting for terminal state:\n' + output.decode(errors='replace'))

            def alive(pid, birth):
                entry = process(pid)
                return entry and entry['birth'] == birth and entry['state'] != 'Z'

            try:
                wait_for(lambda: b'LE8_TEST_PROMPT>' in output)
                command = shlex.join(['bash', str(ROOT / 'le8.sh'), 'c1_correlate',
                                      '--Y', 'cvd_cad', '--biom', 'prot'])
                os.write(master, (command + '\n').encode())
                wait_for(lambda: ready.exists() and pulse.exists())
                pid = json.loads(ready.read_text())['pid']
                while pid != shell:
                    entry = process(pid)
                    self.assertIsNotNone(entry)
                    tracked[pid] = entry['birth']
                    if job_group is None:
                        job_group = entry['group']
                    self.assertEqual(entry['group'], job_group, (pid, entry))
                    self.assertEqual(entry['session'], shell, (pid, entry))
                    pid = entry['parent']
                self.assertGreaterEqual(len(tracked), 4)
                self.assertNotEqual(job_group, os.getpgrp())
                if terminate_logger:
                    os.kill(job_group, signal.SIGTERM)
                else:
                    os.write(master, b'\x1a')  # Ctrl+Z, delivered by the terminal.
                    wait_for(lambda: all(process(p) and process(p)['state'] == 'T' for p in tracked))
                    paused = pulse.read_text()
                    time.sleep(.2)
                    self.assertEqual(pulse.read_text(), paused)
                    os.write(master, b'fg\n')
                    wait_for(lambda: pulse.read_text() != paused)
                    os.write(master, b'\x03')  # Ctrl+C, delivered by the terminal.
                wait_for(lambda: not any(alive(p, birth) for p, birth in tracked.items()))
                self.assertEqual(json.loads(pool.read_text())['leases'], {})
                if not verbose:
                    self.assertTrue(list((work / 'logs').glob('le8.*.log')))
            finally:
                for pid, birth in tracked.items():
                    if alive(pid, birth):
                        os.kill(pid, signal.SIGKILL)
                os.kill(shell, signal.SIGKILL)
                os.close(master)
                os.waitpid(shell, 0)

    def test_ctrl_z_and_ctrl_c_with_logger(self):
        self.exercise()

    def test_ctrl_z_and_ctrl_c_without_logger(self):
        self.exercise(verbose=True)

    def test_terminating_logger_cleans_workers(self):
        self.exercise(terminate_logger=True)


if __name__ == '__main__':
    unittest.main()
