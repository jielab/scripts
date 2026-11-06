import json,os,pathlib,pandas as pd
from types import SimpleNamespace
from shared_cohort import align_projected_cohort,table,sha256
import prsformer
base=pathlib.Path('/work/sph-huangj/scripts/grid/.preflight/20261105')
root=base/'roster-metadata-verified'
source=base/'roster-only/cache/abm/split.rds'
data=prsformer.grid_data_helpers()
cohort=data.read_table(source)[['eid','ancestry','family_id']]
assert len(cohort)==487713
cohort=align_projected_cohort(cohort,root)
assert len(cohort)==486745
settings=SimpleNamespace(cache_dir=str(root/'cache/abm'),out_root=str(root/'abm'),prsformer_root=str(root),split_file='none',seed=20260904)
cohort,_=data.freeze_split(cohort,settings)
roster=root/'prsformer/cache/roster'
prsformer.shared_roster(SimpleNamespace(cache_dir=settings.cache_dir,grid_split=str(root/'abm/grid.split.rds'),roster_dir=str(roster),seed=20260904,check=False))
frame=table(roster/'prsformer.split.tsv.gz');align_projected_cohort(frame,root,frozen=True)
ids=set(frame.id); records=[]
for trait in ['height','ldl','t2dm']:
 for name in ['1csx.scores.rds','2disco.scores.rds']:
  path=pathlib.Path(os.environ['GRID_SCORE_DIR'])/trait/name
  scores=table(path)
  missing=len(ids-set(scores.id)); assert missing==0,(trait,name,missing)
  record=dict(trait=trait,file=str(path),cohort_N=len(ids),missing_ids=missing)
  records.append(record);print('COVERAGE PASS',trait,name,len(ids),flush=True)
groups=table(roster/'families.tsv.gz')
joined=frame.merge(groups,on='id',validate='one_to_one')
assert joined.groupby('group').split.nunique().max()==1
report=dict(status='PASS',validation='Real candidate IDs/families from the completed raw-input preflight; no genotype matrix or training',candidate_source=str(source),candidate_source_sha256=sha256(source),candidate_N=487713,common_N=len(ids),missing_pca_N=968,family_splits_disjoint=True,split_counts=frame.split.value_counts().to_dict(),coverage=records)
(base/'real-cohort-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
