#!/usr/bin/env python3
"""Describe existing frozen outcome-scale benchmark predictions; no model fitting."""
import argparse
from pathlib import Path
from twin import identified,read_table,sha256,json_write
from twin_core import id_hash

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--predictions',required=True,type=Path);p.add_argument('--trait',required=True);p.add_argument('--endpoint',required=True);p.add_argument('--methods',required=True,help='Comma separated method:column entries');p.add_argument('--outcome-column');p.add_argument('--output',required=True,type=Path);p.add_argument('--declare-frozen-without-test-labels',action='store_true')
    a=p.parse_args()
    if not a.declare_frozen_without_test_labels:p.error('This exporter cannot establish training provenance. Explicitly declare previously frozen predictions only after checking their producer.')
    d=identified(read_table(a.predictions));items=[]
    import os
    for pair in a.methods.split(','):
        method,col=pair.split(':',1)
        if col not in d:raise ValueError('Missing prediction column: '+col)
        items.append(dict(name=method,path=os.path.relpath(a.predictions.resolve(),a.output.parent.resolve()),trait=a.trait,endpoint_definition=a.endpoint,value_type='outcome_scale_prediction',prediction_column=col,outcome_column=a.outcome_column,test_roster_sha256=id_hash(d.eid),sha256=sha256(a.predictions),predictions_frozen_without_test_labels=True,provenance='provider declaration; file and roster verified on evaluation'))
    json_write(a.output,dict(methods=items))
if __name__=='__main__':main()
