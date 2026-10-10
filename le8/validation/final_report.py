"""Final report and Shiny: evidence eligibility and aggregate-table contracts."""
import importlib.util
import pandas as pd
import pytest
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def test_shiny_aggregate_reference_and_private_id_guard():
    script = r'''
args <- commandArgs(trailingOnly=TRUE)
e <- new.env(parent=baseenv())
for (expr in parse(args[1])) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c('private_columns', 'has_private_columns')) eval(expr, e)
}
aggregate <- c('run_id','model_id','primary_id','reference_id','n','uno_c_horizon')
stopifnot(!e$has_private_columns(aggregate))
stopifnot(e$has_private_columns(c(aggregate, 'eid')))
stopifnot(e$has_private_columns(c(aggregate, 'sample_id')))
stopifnot(e$has_private_columns(c('reference_id', 'risk')))
stopifnot(!e$has_private_columns(c('feature', 'beta', 'p.value')))
'''
    z=subprocess.run(['Rscript','-e',script,str(ROOT/'shiny/app.R')],capture_output=True,text=True)
    assert z.returncode==0,z.stdout+z.stderr

@pytest.mark.parametrize('layout', ['metrics', 'columns'])
@pytest.mark.parametrize('eligible', ['FALSE', 'TRUE'])
def test_final_dandelion_eligibility_audit(layout, eligible):
    spec=importlib.util.spec_from_file_location('final_dandelion_audit',ROOT/'f/final.py')
    module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    report=module.Report.__new__(module.Report)
    values=dict(primary_eligible=eligible,analysis_class='exploratory fixture')
    table=(pd.DataFrame(dict(metric=list(values),value=list(values.values())))
           if layout=='metrics' else pd.DataFrame([values]))
    report.tables={'cad.prot.dandelion':table}
    report.findings=[]
    report.audit('cad','prot')
    if eligible=='TRUE':
        assert report.findings==[]
    else:
        assert len(report.findings)==1
        finding=report.findings[0]
        assert finding['code']=='DANDELION_EXPLORATORY_ONLY'
        assert finding['severity']=='warning'
        assert 'exclude from independent causal confirmation' in finding['detail']
        assert set(finding['source_roles'].split(';'))=={
            'dandelion','dandelion_lolo','dandelion_native','state_projection','age_models'}
