"""Sample-level wall times and low-frequency LSF progress for the ranked finish step."""
from datetime import datetime, timezone
import os
from pathlib import Path
import threading
from time import monotonic


def mapper_counters(pid):
    """Cumulative child CPU and storage reads sampled every five minutes on Linux."""
    if pid is None:
        return ''
    folder = Path('/proc') / str(pid)
    fields = []
    try:
        stat = (folder/'stat').read_text().rsplit(') ', 1)[1].split()
        cpu = (int(stat[11]) + int(stat[12])) / os.sysconf('SC_CLK_TCK')
        fields.extend((f'pid={pid}', f'cpu_s={cpu:.1f}', f'major_faults={stat[9]}'))
        resident = int(stat[21]) * os.sysconf('SC_PAGE_SIZE')
        fields.append(f'rss_gib={resident / 2**30:.2f}')
        io = dict(row.split(': ', 1) for row in (folder/'io').read_text().splitlines())
        fields.append(f'read_bytes={int(io["read_bytes"])}')
    except (OSError, ValueError, IndexError, KeyError):
        pass
    return ' '.join(fields)


class SampleProgress:
    def __init__(self, sample, interval=300):
        self.sample = sample
        self.interval = interval
        self.phase = None
        self.timing = {}
        self.pairs = 0
        self.lock = threading.Lock()
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._heartbeat, daemon=True)

    def __enter__(self):
        self.thread.start()
        return self

    def _write(self, event, now, details=''):
        stamp = datetime.now(timezone.utc).isoformat(timespec='seconds')
        print(f'RANKED {stamp} sample={self.sample} phase={self.phase} {event} '
              f'elapsed_s={now - self.start:.1f} pairs={self.pairs} {details}'.rstrip(),
              flush=True)

    def _finish(self, now, event):
        if self.phase is not None:
            self.timing[self.phase] = round(now - self.start, 3)
            self._write(event, now)

    def stage(self, phase, **details):
        with self.lock:
            now = monotonic()
            self._finish(now, 'done')
            self.phase, self.start, self.pairs = phase, now, 0
            self.mapper_pid = details.get('pid')
            self._write('start', now, ' '.join(f'{k}={v}' for k, v in details.items()))

    def counted(self, pairs):
        with self.lock:
            self.pairs = pairs

    def _heartbeat(self):
        while not self.stop.wait(self.interval):
            with self.lock:
                if self.phase is not None:
                    self._write('running', monotonic(), mapper_counters(self.mapper_pid))

    def __exit__(self, exc_type, exc, traceback):
        self.stop.set()
        self.thread.join()
        with self.lock:
            self._finish(monotonic(), 'failed' if exc_type else 'done')
