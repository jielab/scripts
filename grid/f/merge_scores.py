#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
ap=argparse.ArgumentParser(); ap.add_argument('--inputs',nargs='+',required=True); ap.add_argument('--output',required=True); a=ap.parse_args(); z=None
for p in a.inputs:
 d=pd.read_csv(p,sep='\t',compression='infer',dtype={'eid':str}); z=d if z is None else z.merge(d,on='eid',how='outer')
Path(a.output).parent.mkdir(parents=True,exist_ok=True); z.to_csv(a.output,sep='\t',index=False,compression='gzip'); print(len(z))
