"""Event-free state queries with a five-year outcome window starting at each query."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from common import outdir, to_device, torch_load, file_digest
from data import CohortDataset
from metrics import binary_metrics
from model import RiskHead
from train import load_world, make_loader, risk_cache, state_cache_key


def labels_at(risk, j, lag, window):
    """Alive, followed, target-free at the query; subsequent recorded diagnosis."""
    t = risk['event_time'][:,j]
    keep = (risk['risk_mask'][:,j] & (risk['exit_time'] > lag)
            & (risk['death_time'] > lag) & (t > lag))
    return keep, (t <= lag + window)


@torch.inference_mode()
def blackout(cfg, requested_device):
    model, meta, device = load_world(cfg, requested_device)
    risk, manifest = risk_cache(cfg, 'test')
    path = outdir(cfg, 'main/results/models') / 'readout.pt'
    saved = torch_load(path)
    names = manifest['endpoint_names']
    if manifest['key'] != state_cache_key(cfg) or saved['states_key'] != manifest['key'] or names != saved['endpoint_names']:
        raise ValueError('Blackout needs the readout and test landmarks from this world model.')
    conf = cfg['blackout']
    window = float(saved['horizon_years'])
    if window != 5:
        raise ValueError('The five-year blackout requires a five-year decoder.')
    lags = list(conf['lags'])
    if not lags or len(set(lags)) != len(lags) or min(lags) <= 0:
        raise ValueError('Blackout lags must be distinct positive years.')
    decoder = RiskHead(len(saved['state_mean']), len(names), saved['config']['hidden_dims'], saved['config']['dropout']).to(device)
    decoder.load_state_dict(saved['model'])
    decoder.eval()
    mean = torch.as_tensor(saved['state_mean'], device=device)
    std = torch.as_tensor(saved['state_std'], device=device)
    panel = []
    for j,name in enumerate(names):
        keep, y = labels_at(risk, j, conf['support_lag'], window)
        cases, controls = int(y[keep].sum()), int((~y[keep]).sum())
        include = bool(saved['eligible'][j] and cases >= conf['min_cases'] and controls >= conf['min_controls'])
        panel.append(dict(endpoint=name, stress_cases=cases, stress_controls=controls, selected=include))
    folder = outdir(cfg, 'main/results')
    pd.DataFrame(panel).to_csv(folder/'blackout.panel.tsv', sep='\t', index=False)
    idx = np.flatnonzero([x['selected'] for x in panel])
    columns = ['lag_years','window_years','endpoint','n','cases','controls','auroc','average_precision']
    rows = []
    if len(idx):
        cohort = CohortDataset(cfg,'test',mode='risk',horizon_years=window,compact_numeric=True)
        loader = make_loader(cohort, cfg['readout']['state_batch_size'], cfg['train']['workers'])
        tmp = outdir(cfg, 'tmp/blackout')
        pred = {lag:np.lib.format.open_memmap(tmp/f'lag{lag}.npy',mode='w+',dtype='float32',
                                             shape=(len(cohort),len(idx))) for lag in lags}
        offset = 0
        for batch in loader:
            end = offset + len(batch['eid'])
            if batch['eid'] != list(risk['eid'][offset:end]):
                raise ValueError('Blackout risk rows differ from test landmarks')
            state = model.encode_history(to_device(batch,device), sample=False, chunk_size=cfg['readout']['history_chunk_size'])
            for lag in lags:
                future = model.advance_state(state, lag, sample=False)
                logits = decoder((model.state_features(future)-mean)/std)
                pred[lag][offset:end] = logits[:,idx].cpu().numpy()
            offset = end
        if offset != len(risk['eid']):
            raise ValueError('Incomplete blackout inference')
        for lag in lags:
            for col,j in enumerate(idx):
                keep, labels = labels_at(risk,j,lag,window)
                metric = binary_metrics(labels[keep], pred[lag][keep,col])
                rows.append(dict(lag_years=lag,window_years=window,endpoint=names[j],**metric))
            del pred[lag]
            (tmp/f'lag{lag}.npy').unlink()
    pd.DataFrame(rows, columns=columns).to_csv(folder/'blackout.tsv',sep='\t',index=False)
    (folder/'blackout.json').write_text(json.dumps(dict(risk_key=manifest['key'], readout_sha256=file_digest(path),
        config=conf, support_endpoints=len(idx), window_years=window),indent=2))
