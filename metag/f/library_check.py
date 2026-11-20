"""Quick library checks; full payload SHA verification is always explicit."""
import argparse
import json
from common import load_config, require_compute
from partitioned_build import database_snapshot
from group_library import compact_snapshot


def check(cfg, full=False):
    union=database_snapshot(cfg,full=full)
    if cfg.get('index_use')=='full':
        result=dict(status='complete',index_use='full',full_union_sha256=full,full_identity=union['identity'])
        print(json.dumps(result,indent=2));return result
    group=compact_snapshot(cfg)
    # Compact native files have structural/observation baselines, no published
    # payload SHA. Do not pretend --full can verify a nonexistent SHA baseline.
    result=dict(status='complete',full_union_sha256=full,group_check='headers_sizes_observations',
                full_identity=union['identity'],group_identity=group['identity'],
                selected_hash_fraction=group['downsampling']['native_options']['selected_hash_fraction'])
    print(json.dumps(result,indent=2));return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config');p.add_argument('--root')
    p.add_argument('--full',action='store_true');p.add_argument('--allow-local',action='store_true');a=p.parse_args()
    if a.full:require_compute(a.allow_local)
    check(load_config(a.config,a.root),a.full)

if __name__=='__main__':main()
