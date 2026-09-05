"""Small Linux subprocess checks; never launch inference or touch source data."""
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]


class RunLogging(unittest.TestCase):
    def test_stdout_stderr_status(self):
        for status in (0, 7):
            with self.subTest(status=status), tempfile.TemporaryDirectory() as directory:
                p = Path(directory)
                (p/'child.sh').write_text(f'echo stdout-marker\necho stderr-marker >&2\nexit {status}\n')
                result = subprocess.run(['bash', str(ROOT/'f/arg.run_logged.sh'), 'build', str(p/'child.sh'), '--arg-dir', str(p)], capture_output=True, text=True)
                self.assertEqual(result.returncode, status)
                logs = list((p/'log').glob('arg.build.*.log'))
                self.assertEqual(len(logs), 1)
                text = logs[0].read_text()
                for expected in ('stdout-marker', 'stderr-marker', '[ARG COMMAND]', f'exit_status={status}'):
                    self.assertIn(expected, text)

    def test_terminal_hangup(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p/'child.sh').write_text('echo before-hangup\nsleep 30\n')
            process = subprocess.Popen(['bash', str(ROOT/'f/arg.run_logged.sh'), 'build', str(p/'child.sh'), '--arg-dir', str(p)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
            try:
                for _ in range(200):
                    logs = list((p/'log').glob('*.log'))
                    if logs and 'before-hangup' in logs[0].read_text():
                        break
                    time.sleep(.01)
                else:
                    self.fail('Logger never received child output')
                os.killpg(process.pid, signal.SIGHUP)
                process.communicate(timeout=5)
                self.assertEqual(process.returncode, 129)
                self.assertIn('exit_status=129', logs[0].read_text())
                self.assertIn('before-hangup', logs[0].read_text())
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.communicate()

    def test_public_defaults_and_startup_error(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p/'vcf').mkdir()
            for extra in ([], ['--unknown-option']):
                result = subprocess.run(['bash', str(ROOT/'arg.sh'), 'check', '--dir-gen', str(p), '--method', 'tsinfer', *extra], capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                if not extra:
                    self.assertIn('chr='+' '.join(map(str, range(1,23)))+' X', result.stdout)
                else:
                    self.assertIn('unknown build option', result.stdout)
            logs = list((p/'arg/log').glob('arg.check.*.log'))
            self.assertEqual(len(logs), 2)
            self.assertTrue(all('[ARG RUN END]' in log.read_text() for log in logs))


if __name__ == '__main__':
    unittest.main()
