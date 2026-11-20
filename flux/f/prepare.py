"""Cheap prepared-cohort reuse; never hashes multi-terabyte raw inputs."""
import json
from pathlib import Path
from common import log, resolve_followup, validate_record
from data import prepare


def _preparation_settings(settings):
    runtime={'landmark','landmark_age','landmark_date','min_clinical_days','max_nodes'}
    return {k:v for k,v in settings.items() if k not in runtime}


def prepare_cohort(cfg,*,replace=False):
    fu=resolve_followup()
    validate_record(cfg, 'data', 'admin_end', fu)
    if cfg['data'].get('admin_end') != fu['date']:
        raise ValueError('Load configuration with the shared follow-up date before preparing data')
    folder=Path(cfg['paths']['out'])/'data' 
    if not replace and (folder/'cohort.sqlite').exists() and (folder/'meta.json').exists():
        meta=json.loads((folder/'meta.json').read_text())
        same=(meta.get('input_root')==str(Path(cfg['paths']['ukb']).resolve())
              and _preparation_settings(meta.get('data_config',{}))==_preparation_settings(cfg['data'])
              and meta.get('seed')==int(cfg['data'].get('seed',cfg['seed'])))
        for source in meta.get('sources',[]):
            p=Path(source['path'])
            if not p.exists(): continue  # A prepared cohort remains usable when the source disk is offline.
            stat=p.stat()
            same &= stat.st_size==source['bytes'] and stat.st_mtime_ns==source['mtime_ns']
        if same:
            log('Reusing data/cohort.sqlite; use prepare --replace after adding new raw source files.')
            return meta
        if (Path(cfg['paths']['out'])/'model/flux.pt').exists():
            raise ValueError('Prepared inputs/settings changed and fitted models exist. Use a new --dir-out or --replace to rebuild and retrain.')
    return prepare(cfg)
