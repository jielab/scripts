# Foreground console

Analysis entry points source `console.sh`. The driver stays attached to the
terminal; existing worker counts, thread counts and parallel scheduling remain
unchanged. The console shows startup, major stage transitions, diagnostics, and
completion or failure. Full output is saved under `/mnt/d/analysis/<module>/logs/`,
using the entry script's parent directory as the module name (for example, `le8`,
`gwas`, or `gu`); existing module logs are retained. GU summaries count completed
or skipped units and failures.

`gwas_format.sh` includes the requested module in its console label and log
filename, for example `gwas_format.mplot.<timestamp>.<pid>.log`. Combined modules
use underscores (`thin,mplot` becomes `thin_mplot`), and an omitted module is
labelled `all`. Per-GWAS START/DONE/FAIL messages also include the module.

Ctrl+C, SIGHUP and SIGTERM stop the driver and its workers. Termination escalates
after 10 seconds. GU also retains its memory-scope cleanup. Explicit background
and cluster submissions retain their separate lifetimes. Prefer Ctrl+C before
closing a terminal: terminal applications do not all deliver SIGHUP identically.

Set `SCRIPT_VERBOSE=1` to use the original console output, or `SCRIPT_LOG_DIR`
to change the new log directory. Help and dry-run output are shown directly.
