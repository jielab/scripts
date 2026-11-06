"""Cold cohort preparation is independent of all PRS score/fitted-model files."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'f'))
import prsformer
from grid_paths import read_config
from shared_cohort import prepare_frozen_roster, load_cohort, preparation_options


def test_portable_builtin_paths_honor_explicit_analysis_root(tmp_path,monkeypatch):
    monkeypatch.setenv('submit','hpc')
    monkeypatch.setenv('ANALYSIS_ROOT','/scratch/jobs/analysis')
    monkeypatch.setenv('GRID_OUTPUT_ROOT','/scratch/custom_grid')
    monkeypatch.setenv('DIR_DATA','/data/lab')
    monkeypatch.setenv('DIR_WORK','/work/lab')
    c=read_config(ROOT/'config/abm.blocks.json')
    assert c['prsformer_cache']=='/scratch/custom_grid/prsformer/cache'
    assert c['ld_block_dir']=='/data/lab/LD/block'
    c=read_config(ROOT/'config/abm.json',tmp_path/'other')
    assert c['block_config']==str(ROOT/'config/abm.blocks.json')
    c=read_config(ROOT/'config/abm.blocks.json',tmp_path/'other')
    assert c['prsformer_cache']==str(tmp_path/'other/prsformer/cache')


def test_raw_grouped_roster_uses_relationship_components_without_any_scores(tmp_path):
    n=600; rng=np.random.default_rng(37)
    gen=tmp_path/'gen';gen.mkdir()
    (gen/'chr1.pgen').write_bytes(b'metadata-only-fixture')
    (gen/'chr1.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n1\t20\trs1\tA\tG\n')
    pd.DataFrame({'#IID':[str(i) for i in range(n)]}).to_csv(gen/'chr1.psam',sep='\t',index=False)
    frame=pd.DataFrame(dict(eid=[str(i) for i in range(n)],height=rng.normal(170,8,n),
        ldl=rng.normal(3,.3,n),t2dm=rng.binomial(1,.2,n),age=rng.normal(55,4,n),
        sex=rng.integers(0,2,n),PC1=rng.normal(size=n),PC2=rng.normal(size=n),genetic_ancestry='EUR'))
    pheno=tmp_path/'pheno.tsv';frame.to_csv(pheno,sep='\t',index=False)
    rel=tmp_path/'ukb_rel.dat'
    # Include a path through an individual outside the genotype roster: the
    # full relation graph, not just retained edges, defines family components.
    pd.DataFrame(dict(ID1=['0','outside','2','4'],ID2=['outside','1','3','5'],Kinship=[.2]*4)).to_csv(rel,sep='\t',index=False)
    root=tmp_path/'analysis';cache=root/'prsformer/cache'
    args=prsformer.parser().parse_args(['prepare','--cache-dir',str(cache),'--dir-gen',str(gen),
        '--pheno-file',str(pheno),'--chrs','1','--t2dm-col','t2dm','--remove','none'])
    prepare_frozen_roster(args,root,rel)
    groups=pd.read_csv(args.split_group_file,sep='\t',dtype=str).set_index('eid')
    split=pd.read_csv(args.split_file,sep='\t',dtype=str).set_index('eid')
    assert groups.loc['0','group']==groups.loc['1','group']
    joined=split.join(groups)
    assert joined.groupby('group').split.nunique().max()==1
    assert set(joined.split)=={'train','validation','test'}
    assert len(joined)==n
    assert args.require_shared_cohort is True
    assert not (root/'csx').exists()
    assert not (root/'cache/abm/prepared.joblib').exists()
    assert not list(root.rglob('model.pt'))
    evidence=json.loads((root/'cache/shared_preparation/GROUPS_VERIFIED.json').read_text())
    assert evidence['split_uses_outcome_values'] is False
    # Export a small standalone cohort from exactly this fixed split to exercise
    # verification against the actual relation source, without producing dosage.
    d=frame.merge(split.reset_index(),on='eid').merge(groups.reset_index(),on='eid')
    path=tmp_path/'cohort.tsv';d.to_csv(path,sep='\t',index=False)
    loaded,meta=load_cohort(root,'height',cohort=path)
    assert meta['groups_verified'] is True
    assert set(loaded.loc[loaded.split.eq('test'),'group_id']).isdisjoint(
        set(loaded.loc[loaded.split.eq('train'),'group_id']))
    rel.write_text(rel.read_text()+'6\t7\t0.2\n')
    with pytest.raises(ValueError,match='family source changed'):
        load_cohort(root,'height',cohort=path)


def test_trait_feature_preparation_keeps_missing_y_without_other_traits_scores(tmp_path):
    import joblib
    from data import write_rds
    from shared_cohort import complete_score_features, prepared_for_trait, sha256
    from twin import load_source
    cache=tmp_path/'prsformer/cache';cache.mkdir(parents=True)
    n=40
    d=pd.DataFrame(dict(eid=[str(i) for i in range(n)], group=[f'g{i//2}' for i in range(n)],
        split=['train']*20+['test']*20, height=np.arange(n,dtype=float), ldl=np.arange(n,dtype=float)+3,
        age=np.arange(n,dtype=float),sex=np.arange(n)%2,PC1=np.arange(n)/2,PC2=np.arange(n)/3))
    d.loc[[2,30],'height']=np.nan
    d.to_csv(cache/'data.tsv.gz',sep='\t',index=False)
    (cache/'prepare.json').write_text(json.dumps(dict(options=dict(covariates='age,sex,PC1,PC2'))))
    cfg=json.loads((ROOT/'config/abm.blocks.json').read_text())
    cfg.update(prsformer_cache=str(cache),score_root=str(tmp_path/'scores'))
    config=tmp_path/'blocks.json';config.write_text(json.dumps(cfg))
    score=pd.DataFrame({'eid':d.eid,**{'csx.'+p:np.arange(n)*(i+1.) for i,p in enumerate(['AFR','EAS','EUR','SAS'])}})
    target=tmp_path/'scores/height';target.mkdir(parents=True)
    write_rds(score,target/'1csx.scores.rds')
    args=SimpleNamespace(analysis_root=str(tmp_path),traits='height',block_config=str(config))
    result=complete_score_features(args,{})
    prepared=prepared_for_trait(tmp_path,'height')
    assert prepared==result['height']
    frame,groups,metadata=load_source(prepared,'height')
    assert len(frame)==n and frame.y.isna().sum()==2
    assert len(groups['csx'])==4 and metadata['labels_kept_missing_until_evaluation']
    assert not (tmp_path/'scores/ldl').exists()
    before=sha256(prepared)
    target=tmp_path/'scores/ldl';target.mkdir(parents=True)
    write_rds(score,target/'1csx.scores.rds')
    args.traits='ldl';complete_score_features(args,{})
    assert sha256(prepared)==before
    assert prepared_for_trait(tmp_path,'ldl')!=prepared


@pytest.mark.parametrize('cuda',['12.6','13.0'])
def test_prsformer_accepts_only_matching_official_pinned_wheels(monkeypatch,cuda):
    natten=SimpleNamespace(__version__='0.21.7+torch2120cu'+cuda.replace('.',''),HAS_LIBNATTEN=True,na1d=lambda:None)
    torch=SimpleNamespace(__version__='2.12.0',version=SimpleNamespace(cuda=cuda),
        cuda=SimpleNamespace(is_available=lambda:False))
    monkeypatch.setattr(prsformer,'import_torch',lambda:torch)
    original=prsformer.importlib.import_module
    monkeypatch.setattr(prsformer.importlib,'import_module',lambda name: natten if name=='natten' else original(name))
    assert prsformer.neighborhood_runtime()==(natten,'0.21.7')
    natten.HAS_LIBNATTEN=False
    with pytest.raises(RuntimeError,match='compiled extension'):
        prsformer.neighborhood_runtime()
    natten.HAS_LIBNATTEN=True;torch.__version__='2.13.0'
    with pytest.raises(RuntimeError,match='Unsupported'):
        prsformer.neighborhood_runtime()
