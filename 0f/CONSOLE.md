# Console output policy

Screen output is event-driven: startup/log location, major module transitions,
completion, warnings and errors. Do not add periodic wall-clock or "still running"
messages, even for analyses that run for several days.

`console_run.py` retains the complete child output in `SCRIPT_LOG_DIR` (default
`/mnt/d/analysis/<module>/logs`, where module is the entry script's parent directory,
such as `le8`, `gwas`, or `gu`). Warnings/errors are not limited to the first three. Process
exit codes and cancellation are preserved. Timing used for deadlines, cancellation,
file provenance or numerical results is separate from console progress output.

The shared wrapper and the dedicated ARG-Needle, THREADS, SINGER and PhyML runners
have no timed console progress reports. Existing in-memory processes must be
restarted to use edited Python code. `SCRIPT_VERBOSE=1` bypasses the screen filter
and exposes the underlying program's output.
