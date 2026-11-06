import os,json,pathlib,pandas as pd
from types import SimpleNamespace
from shared_cohort import preparation_options,prepare_frozen_roster,table
root=pathlib.Path('/work/sph-huangj/scripts/grid/.preflight/20261105/roster-fixed')
a,block=preparation_options(SimpleNamespace(analysis_root=str(root)))
prepare_frozen_roster(a,root)
roster=pd.read_csv(a.split_file,sep='\t',dtype=str)
ids=set(roster.eid)
assert ids <= set(table(os.environ['GRID_PCA_FILE']).id)
for trait in ['height','ldl','t2dm']:
 for name in ['1csx.scores.rds','2disco.scores.rds']:
  scores=table(pathlib.Path(block['score_root'])/trait/name)
  assert ids <= set(scores.id),(trait,name,len(ids-set(scores.id)))
  print('COVERAGE PASS',trait,name,len(ids),flush=True)
groups=pd.read_csv(a.split_group_file,sep='\t',dtype=str)
joined=roster.merge(groups,on='eid',validate='one_to_one')
assert joined.groupby('group').split.nunique().max()==1
print('FAMILY SPLITS PASS',len(ids),flush=True)
print((root/'cache/shared_preparation/PCA_COHORT.json').read_text())
