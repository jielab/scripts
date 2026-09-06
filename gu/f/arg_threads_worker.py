#!/usr/bin/env python3
"""Launch the installed Threads CLI with a bounded, private Ray runtime."""
import importlib.metadata
import os
import sys

if sys.argv[1:2] == ['infer']:
    import ray
    original_init = ray.init

    def bounded_init(*args, **kwargs):
        kwargs.setdefault('num_cpus', int(os.environ['REFGEN_THREADS']))
        kwargs.setdefault('object_store_memory', 512 * 1024**2)
        kwargs.setdefault('include_dashboard', False)
        return original_init(*args, **kwargs)

    ray.init = bounded_init

entry, = [e for e in importlib.metadata.distribution('threads_arg').entry_points
          if e.group == 'console_scripts' and e.name == 'threads']
entry.load()()
