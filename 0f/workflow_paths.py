"""Resolve historical repository defaults through the single 0os.sh table.

Use before applying explicit CLI values; custom paths are never searched for or
silently redirected. Python helpers can use the same table without shell setup.
"""
from __future__ import annotations
from functools import lru_cache
import os
from pathlib import Path
import subprocess

KEYS = ('SCRIPTS_ROOT','DIR_WORK','DIR_DAT','DIR_BIG','DIR_TMP','DIR_SOFTWARE',
        'DIR_FILES','DIR_UKB','DIR_ANNOT','DIR_GEN','DIR_GWAS','DIR_LD','DIR_METAG','ANALYSIS_ROOT')

@lru_cache(maxsize=8)
def _fallback_table(environment_key):
    table=Path(__file__).with_name('0os.sh')
    command='source "$1" || exit; shift; for name in "$@"; do printf "%s\\0%s\\0" "$name" "${!name}"; done'
    raw=subprocess.check_output(['bash','-c',command,'workflow-paths',str(table),*KEYS],text=True)
    parts=raw.rstrip('\0').split('\0')
    return dict(zip(parts[::2],parts[1::2]))

def resolved_table():
    key=tuple(sorted((k,v) for k,v in os.environ.items() if k.startswith(('DIR_','_OS_DEFAULT_')) or k in (*KEYS,'submit','SUBMIT','DATE_SCRATCH','PATH','HOME')))
    return _fallback_table(key)

def relocate_defaults(value):
    if isinstance(value,dict): return {k:relocate_defaults(v) for k,v in value.items()}
    if isinstance(value,list): return [relocate_defaults(v) for v in value]
    if isinstance(value,tuple): return tuple(relocate_defaults(v) for v in value)
    if not isinstance(value,str): return value
    t=resolved_table()
    mappings=[('/mnt/d/scripts',t['SCRIPTS_ROOT']),('/mnt/d/analysis',t['ANALYSIS_ROOT']),
              ('/mnt/f/analysis',t['ANALYSIS_ROOT']),('/mnt/d/software',t['DIR_SOFTWARE']),
              ('/mnt/f/software',t['DIR_SOFTWARE']),
              ('/mnt/d/data/ukb',t['DIR_UKB']),('/mnt/d/ukb',t['DIR_UKB']),
              ('/mnt/d/files',t['DIR_FILES'])]
    mappings += [('/mnt/f/'+cat,t[key]) for cat,key in [('annot','DIR_ANNOT'),('gen','DIR_GEN'),('gwas','DIR_GWAS'),('LD','DIR_LD'),('metag','DIR_METAG'),('ukb','DIR_UKB')]]
    mappings += [('/mnt/f',t['DIR_BIG']),('/mnt/d',t['DIR_DAT'])]
    for old,new in mappings:
        if value==old or value.startswith(old+'/'): return new.rstrip('/')+value[len(old):]
    return value
