"""Bounded parallel input preparation, inherited into the workflow's memory scope."""
import multiprocessing
import queue
import traceback
from pathlib import Path

from metacache_inputs import shards, SHARD_BASES


def prepare_worker(part, cache, identity, messages, max_bases=SHARD_BASES):
    try:
        for folder, key in shards(Path(part), Path(cache), identity, max_bases=max_bases):
            messages.put(('shard', str(folder), key))
        messages.put(('done',))
    except BaseException:
        messages.put(('error', traceback.format_exc()))


class PreparationPool:
    """At most workers source groups and four ready input shards per worker.

    Workers prepare while the parent builds on the GPU. Only small path messages
    cross queues; DNA stays in bounded buffers or in persistent files on disk.
    """
    def __init__(self, tasks, cache, identity, workers=2, max_bases=SHARD_BASES):
        if not 1 <= workers <= 4:
            raise ValueError('input preparation workers must be between 1 and 4')
        self.tasks = iter(tasks)
        self.cache, self.identity, self.workers = cache, identity, workers
        self.max_bases = max_bases
        self.ctx = multiprocessing.get_context('spawn')
        self.active = []

    def start_next(self):
        task = next(self.tasks, None)
        if task is None:
            return
        index, part = task
        messages = self.ctx.Queue(maxsize=4)
        process = self.ctx.Process(target=prepare_worker,
                                   args=(str(part), str(self.cache), self.identity, messages, self.max_bases))
        self.active.append({'index': index, 'process': process, 'queue': messages, 'done': False})
        process.start()

    def __enter__(self):
        try:
            for _ in range(self.workers):
                self.start_next()
        except BaseException:
            self.__exit__(None, None, None)
            raise
        return self

    def records(self, work):
        while True:
            try:
                message = work['queue'].get(timeout=0.5)
            except queue.Empty:
                if not work['process'].is_alive():
                    raise RuntimeError(f"Input preparation worker {work['index']} exited without completion")
                continue
            if message[0] == 'error':
                raise RuntimeError(f"Input preparation failed for group {work['index']}:\n{message[1]}")
            if message[0] == 'done':
                work['done'] = True
                return
            if message[0] != 'shard':
                raise RuntimeError('Invalid input preparation message')
            yield Path(message[1]), message[2]

    def __iter__(self):
        while self.active:
            work = self.active[0]
            yield work['index'], self.records(work)
            if not work['done']:
                raise RuntimeError('Input group was not completely consumed')
            work['process'].join(timeout=5)
            if work['process'].exitcode != 0:
                raise RuntimeError(f"Input preparation worker {work['index']} failed to exit cleanly")
            work['queue'].close()
            work['queue'].join_thread()
            self.active.pop(0)
            self.start_next()

    def __exit__(self, *exc):
        for work in self.active:
            if work['process'].is_alive():
                work['process'].terminate()
        for work in self.active:
            if work['process'].pid is not None:
                work['process'].join(timeout=5)
                if work['process'].is_alive():
                    work['process'].kill()
                    work['process'].join(timeout=5)
            work['queue'].cancel_join_thread()
            work['queue'].close()
        self.active.clear()
