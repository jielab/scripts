"""GWAS paths use the same 0os.sh table from shell and direct Python entries."""
import importlib.util
import os
from pathlib import Path
p=Path(__file__).resolve().parents[2]/'0f/workflow_paths.py'
spec=importlib.util.spec_from_file_location('gwas_shared_paths',p)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
table=module.resolved_table()
WORK=Path(table['DIR_WORK'])
DATA=Path(os.environ.get('GWAS_DATA_ROOT',table['DIR_BIG']))
FILES=Path(os.environ.get('GWAS_FILES_ROOT',table['DIR_FILES']))
SOFTWARE=Path(os.environ.get('GWAS_SOFTWARE_ROOT',table['DIR_SOFTWARE']))
ANALYSIS=Path(os.environ.get('GWAS_ANALYSIS_ROOT',str(Path(table['ANALYSIS_ROOT'])/'gwas')))
