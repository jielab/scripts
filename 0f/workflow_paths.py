"""Resolve historical repository defaults through the single 0os.sh table.

Use before applying explicit CLI values; custom paths are never searched for or
silently redirected. Python helpers can use the same table without shell setup.
"""
from __future__ import annotations
from functools import lru_cache
import os
import re
from pathlib import Path
import subprocess

@lru_cache(maxsize=8)
def _fallback_table(environment_key):
    table=Path(__file__).with_name('0os.sh')
    command='source "$1" || exit; for name in "${OS_PATH_KEYS[@]}"; do printf "%s\\0%s\\0" "$name" "${!name}"; done'
    raw=subprocess.check_output(['bash','-c',command,'workflow-paths',str(table)],text=True)
    parts=raw.rstrip('\0').split('\0')
    return dict(zip(parts[::2],parts[1::2]))

def resolved_table():
    key=tuple(sorted((k,v) for k,v in os.environ.items() if k.startswith(('DIR_','_OS_DEFAULT_')) or k.startswith(('OS_', 'ANALYSIS_', 'SCRIPTS_')) or k in ('submit','SUBMIT','DATE_SCRATCH','PATH','HOME','LSB_JOBID','LSF_ENVDIR','LSF_SERVERDIR')))
    return _fallback_table(key)

def relocate_defaults(value):
    if isinstance(value,dict): return {k:relocate_defaults(v) for k,v in value.items()}
    if isinstance(value,list): return [relocate_defaults(v) for v in value]
    if isinstance(value,tuple): return tuple(relocate_defaults(v) for v in value)
    if not isinstance(value,str): return value
    t=resolved_table()
    value = expand_path_roots(value)
    mappings=[('/mnt/d/scripts',t['SCRIPTS_ROOT']),('/mnt/d/analysis',t['ANALYSIS_ROOT']),
              ('/mnt/f/analysis',t['ANALYSIS_ROOT']),('/mnt/d/software',t['DIR_SOFTWARE']),
              ('/mnt/f/software',t['DIR_SOFTWARE']),
              ('/mnt/d/phe/ukb',t['DIR_UKB']),('/mnt/d/data/ukb',t['DIR_UKB']),('/mnt/d/ukb',t['DIR_UKB']),
              ('/mnt/d/files',t['DIR_FILES'])]
    mappings += [('/mnt/f/'+cat,t[key]) for cat,key in [('annot','DIR_ANNOT'),('gen','DIR_GEN'),('gwas','DIR_GWAS'),('LD','DIR_LD'),('metag','DIR_METAG'),('ukb','DIR_UKB')]]
    mappings += [('/mnt/f',t['DIR_BIG']),('/mnt/d',t['DIR_DAT'])]
    for old,new in mappings:
        if value==old or value.startswith(old+'/'): return new.rstrip('/')+value[len(old):]
    return value


def configured_path(key, *parts, default_mode=None):
    """Get a named root from 0os.sh; explicit environment overrides survive."""
    if default_mode and not os.environ.get("submit") and not os.environ.get("SUBMIT"):
        environment = dict(os.environ, submit=default_mode)
        table = Path(__file__).with_name('0os.sh')
        command = 'source "$1" || exit; printf "%s" "${!2}"'
        root = subprocess.check_output(['bash', '-c', command, 'workflow-path', str(table), key], env=environment, text=True)
    else:
        root = resolved_table()[key]
    return str(Path(root).joinpath(*parts))


def expand_path_roots(value):
    """Expand declared root tokens, preserving literal/custom paths and non-path tokens."""
    if isinstance(value, dict): return {k: expand_path_roots(v) for k,v in value.items()}
    if isinstance(value, list): return [expand_path_roots(v) for v in value]
    if not isinstance(value, str) or "${" not in value: return value
    roots = resolved_table()
    return re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}", lambda m: roots.get(m.group(1), m.group(0)), value)
