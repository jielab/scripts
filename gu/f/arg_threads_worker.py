#!/usr/bin/env python3
"""Launch the installed Threads CLI with a bounded, private Ray runtime."""
import importlib.metadata
import os
import sys

if sys.argv[1:2] == ['infer']:
    # Ray otherwise rewrites even an explicit loopback IP to an external IP.
    # Despite its name, this switch controls that behavior on Linux too.
    os.environ['RAY_ENABLE_WINDOWS_OR_OSX_CLUSTER'] = '0'
    import ray
    original_init = ray.init

    def bounded_init(*args, **kwargs):
        # This pipeline is single-host. Avoid WSL/VPN interface changes breaking
        # Ray heartbeats, and never attach to a cluster from RAY_ADDRESS.
        kwargs['address'] = 'local'
        kwargs['_node_ip_address'] = '127.0.0.1'
        kwargs.setdefault('num_cpus', int(os.environ['REFGEN_THREADS']))
        kwargs.setdefault('object_store_memory', 512 * 1024**2)
        kwargs.setdefault('include_dashboard', False)
        return original_init(*args, **kwargs)

    ray.init = bounded_init

entry, = [e for e in importlib.metadata.distribution('threads_arg').entry_points
          if e.group == 'console_scripts' and e.name == 'threads']
entry.load()()
