"""ABM-only feature preparation from an existing immutable comparison roster."""
from pathlib import Path
import json
import os
import fcntl
import joblib
import pandas as pd


def prepare_abm(args):
    from twin import load_source, identified, read_table, verify_existing_prsformer_roster, frozen_comparison_roster, PreparedSourceChanged
    import main as preparation
    root=Path(args.analysis_root);cache=Path(args.cache_dir)
    target=Path(args.prepared) if args.prepared else cache/'prepared.joblib'
    r=frozen_comparison_roster(root)
    def validate(path,allow_unfinished=False):
        for trait in args.traits.split(','):
            d,_,_=load_source(path,trait,allow_unfinished_features=allow_unfinished)
            verify_existing_prsformer_roster(root,d,d.loc[d.split.eq('test'),'eid'])
    def finish_existing(path):
        validate(path,allow_unfinished=True)
        value=joblib.load(path)
        if value.get('state')!='prepared':
            if path!=target:preparation.save_cache(value,target)
            return
        # The initial PRSformer shared-cohort stage deliberately stops before
        # ABM genetic feature preparation. Complete that stage without training
        # or generating another cohort, using the original preparation settings.
        options=preparation.parser().parse_args([])
        for key,item in value.get('settings',{}).items():
            if hasattr(options,key):setattr(options,key,item)
        options.cache_dir=str(cache);options.out_root=args.out_root;options.prsformer_root=str(root)
        if isinstance(options.traits,str):options.traits=options.traits.split(',')
        value=preparation.load_module('data').evolution_data(options,value)
        value['state']='features_ready'
        value['metadata']['abm_preparation']='completed_features_on_existing_frozen_roster_without_training'
        temporary=target.with_name(target.name+'.validating')
        preparation.save_cache(value,temporary)
        try:validate(temporary);os.replace(temporary,target)
        finally:temporary.unlink(missing_ok=True)
    cache.mkdir(parents=True,exist_ok=True)
    with (cache/'abm.prepare.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        rebuild_settings={}
        # An outdated cache is optional acceleration, not a reason to block a run.
        # Keep it intact until the replacement passes all frozen-roster checks.
        candidates=[target] if target.exists() else [root/'cache/abm/prepared.joblib']
        for candidate in candidates:
            if not candidate.is_file():continue
            try:
                finish_existing(candidate)
            except PreparedSourceChanged as exc:
                rebuild_settings=joblib.load(candidate).get('settings',{})
                print('[GRID prepare] Source updated; rebuilding ABM features on the frozen roster: '+str(exc),flush=True)
                break
            print('[GRID prepare] ABM features verified/reused: '+str(candidate),flush=True);return
        options=preparation.parser().parse_args([])
        metadata=root/'prsformer/cache/prepare.json'
        if metadata.is_file():
            saved=json.loads(metadata.read_text()).get('options',{})
            for key in ('pheno_file','ancestry_file','group_col','dir_gen','remove','covariates','height_col','ldl_col','t2dm_col','seed'):
                if key in saved:setattr(options,key,saved[key])
        for key,value in rebuild_settings.items():
            if hasattr(options,key):setattr(options,key,value)
        if getattr(args,'prepare_config',None):
            settings=json.loads(Path(args.prepare_config).read_text())
            allowed={'pheno_file','score_dir','pca_file','ancestry_file','group_col','dir_gen','gwas_dir','snpinfo','weights_file','remove','covariates','distance_pcs','height_col','ldl_col','t2dm_col','annotation_file','build','plink2','threads','seed'}
            if set(settings)-allowed:raise ValueError('Unsupported ABM preparation option: '+str(set(settings)-allowed))
            for key,value in settings.items():setattr(options,key,value)
        frozen=r[['eid','group','split']].rename(columns={'group':'family_id'}).copy()
        frozen['split']=frozen.split.replace({'validation':'train'})
        if 'target' in r:frozen['ancestry']=r.target.astype(str)
        split_path=cache/'abm.frozen_outer.tsv.gz';family_path=cache/'abm.frozen_families.tsv.gz';keep_path=cache/'abm.frozen.keep'
        frozen.to_csv(split_path,sep='\t',index=False)
        frozen[['eid','family_id']].rename(columns={'family_id':'group'}).to_csv(family_path,sep='\t',index=False)
        frozen[['eid']].rename(columns={'eid':'#IID'}).to_csv(keep_path,sep='\t',index=False)
        options.traits=args.traits.split(',');options.cache_dir=str(cache);options.out_root=args.out_root
        options.prsformer_root=str(root);options.split_file=str(split_path);options.split_group_file=str(family_path);options.keep=str(keep_path)
        data=preparation.load_module('data')
        prepared=data.prepare_data(options)
        for trait,d in prepared['tables'].items():
            verify_existing_prsformer_roster(root,d,d.loc[d.split.eq('test'),'eid'])
        prepared=data.evolution_data(options,prepared)
        prepared['settings']=preparation.settings(options);prepared['state']='features_ready'
        prepared['source_signature']=preparation.source_signature()
        prepared['metadata']['abm_preparation']='frozen_PRSformer_roster_no_training_no_new_outer_split'
        # No partial prepared.joblib becomes visible on failure.
        temporary=target.with_name(target.name+'.validating')
        preparation.save_cache(prepared,temporary)
        try:
            validate(temporary);os.replace(temporary,target)
        finally:temporary.unlink(missing_ok=True)
        print('[GRID prepare] ABM features ready; PRSformer was read-only',flush=True)
