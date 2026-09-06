"""Memory-bounded adapter for ARG-Needle's hashed threading path.

Keep the installed package untouched. Transform the dense-posterior operations
and native batch size, refusing unknown upstream code rather than falling back.
For each site, smoothing reads only the posterior of the selected cousin.
Across a smooth interval the selected cousin is constant, so storing those
selected values gives the same float64 means without a samples-by-sites matrix.
Hashing, ASMC pair requests, random draws, MAP boundaries and threading stay
upstream. Native ASMC workspaces use batches of 16 instead of the default 64.
"""
import inspect
import logging
import textwrap
import importlib.util
from pathlib import Path
import sys
import sysconfig


def backend_path():
    return (Path.home() / '.cache/gu/asmc-memory/v1.4.0-p2' / sys.implementation.cache_tag /
            ('asmc_python_bindings' + sysconfig.get_config_var('EXT_SUFFIX')))


def _load_backend():
    name = 'asmc.asmc.asmc_python_bindings'
    if name in sys.modules:
        if getattr(sys.modules[name], '_gu_bounded_workspace', 0) == 2:
            return
        raise RuntimeError('Load argneedle_memory before importing arg_needle or asmc.asmc')
    path = backend_path()
    if not path.is_file():
        raise RuntimeError('Missing bounded ASMC backend. Run grid Python on '
                           f'{Path(__file__).with_name("build_asmc_memory.py")} '
                           '--source /path/to/ASMC-1.4.0')
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if getattr(module, '_gu_bounded_workspace', 0) != 2:
        raise RuntimeError('ASMC backend does not contain the workspace fix')
    sys.modules[name] = module


def _rewrite(function, replacements):
    source = textwrap.dedent(inspect.getsource(function))
    for old, new in replacements:
        if source.count(old) != 1:
            raise RuntimeError(
                f"Unsupported ARG-Needle implementation in {function.__name__}: "
                f"expected exactly one {old!r}")
        source = source.replace(old, new, 1)
    namespace = {}
    exec(compile(source, inspect.getfile(function) + "[bounded-memory]", "exec"),
         function.__globals__, namespace)
    return namespace[function.__name__]


def install():
    _load_backend()
    import arg_needle.inference as inference
    import arg_needle.decoders as decoders
    from arg_needle.decoders import ASMCDecoder

    if getattr(inference.thread_samples, "_bounded_memory", False):
        return
    thread = _rewrite(inference.thread_samples, [
        ("np.full((start_thread_id + num_next_samples - 1, len(posterior_phys_pos)), np.nan)",
         "np.full(len(posterior_phys_pos), np.nan)"),
        ("np.mean(tmrca_mean[indices[begin], begin:end])",
         "np.mean(tmrca_mean[begin:end] if hash_topk > 0 else "
         "tmrca_mean[indices[begin], begin:end])"),
    ])
    decode = _rewrite(ASMCDecoder.compute_with_hashing, [
        ("tmrca_mean[other_ids, from_pos:to_pos] = batch_mean\n"
         "            foo = np.argmin(batch_mean, axis=0)",
         "foo = np.argmin(batch_mean, axis=0)\n"
         "            tmrca_mean[from_pos:to_pos] = "
         "batch_mean[foo, np.arange(to_pos - from_pos)]"),
    ])
    make_decoder = _rewrite(decoders.make_asmc_decoder, [
        # Keep all 64 candidates. Only the internal SIMD batch/workspace shrinks.
        ("asmc_obj = ASMC(params)",
         "params.batchSize = 16\n    asmc_obj = ASMC(params)"),
    ])
    # Commit together only after all upstream functions pass compatibility checks.
    thread._bounded_memory = True
    inference.thread_samples = thread
    ASMCDecoder.compute_with_hashing = decode
    decoders.make_asmc_decoder = make_decoder
    inference.make_asmc_decoder = make_decoder
    logging.info("Bounded-memory threading enabled: one float64 posterior per site; "
                 "ASMC batches <=16; native scratch pages discarded after each window")
