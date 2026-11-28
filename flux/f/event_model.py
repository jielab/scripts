"""Local EHR initialization and held-out next-event / interval evaluation."""
import itertools
import sqlite3

import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F
from common import log, outdir, to_device


def initialize_embeddings(base,cfg,meta):
    """Factor training-only adjacent-day transition probabilities for EHR initialization."""
    from sklearn.utils.extmath import randomized_svd
    token={name:i for i,name in enumerate(meta['vocab'])}
    counts=np.zeros((len(token),len(token)),np.float64)
    db=outdir(cfg,'main/results/data')/'cohort.sqlite'
    transitions=0
    with sqlite3.connect(f'file:{db}?mode=ro',uri=True) as conn:
        rows=conn.execute("""SELECT p.eid,e.day,e.code FROM participants p
            JOIN events e ON e.eid=p.eid WHERE p.split='train' ORDER BY p.eid,e.day,e.code""")
        for _,person in itertools.groupby(rows,key=lambda r:r[0]):
            previous=[token['<START>']]
            for _,day in itertools.groupby(person,key=lambda r:r[1]):
                current=sorted({token.get(r[2],token['<UNK>']) for r in day})
                counts[np.ix_(previous,current)] += 1.0/(len(previous)*len(current))
                previous=current
                transitions+=1
    source=counts.sum(axis=1)
    target=counts.sum(axis=0)
    if transitions==0:
        raise ValueError('No training clinical events for EHR initialization')
    # Truncated SVD of the positive pointwise mutual information matrix.
    lift=np.divide(counts*counts.sum(),source[:,None]*target[None,:],
                   out=np.zeros_like(counts),where=(source[:,None]*target[None,:])>0)
    width=min(len(token),base.embedding_dim)
    left,sv,right=randomized_svd(np.log(np.maximum(lift,1)),n_components=width,random_state=int(cfg['seed']))
    # A terminal token can have only incoming events. Include its context embedding.
    emb=(left+right.T)*np.sqrt(sv)
    norms=np.linalg.norm(emb,axis=1)
    used=(norms>1e-8) & ((source+target)>0)
    used[0]=False
    emb[used] *= np.sqrt(width)/norms[used,None]
    indices=np.flatnonzero(used)
    with torch.no_grad():
        base.token_embedding.weight[indices,:width]=torch.as_tensor(emb[used],device=base.device,dtype=base.dtype)
    log(f'EHR initialization: {transitions:,} training clinical days, {len(indices)} token embeddings')


@torch.inference_mode()
def heldout_dynamics(cfg,device):
    from train import load_world,make_loader
    from data import CohortDataset
    model,meta,dev=load_world(cfg,device)
    ds=CohortDataset(cfg,'test',mode='pretrain',max_nodes=cfg['data']['max_nodes'],compact_numeric=True)
    loader=make_loader(ds,cfg['readout']['state_batch_size'],cfg['train']['workers'])
    totals=np.zeros((3,3))
    for batch in loader:
        b=to_device(batch,dev)
        out=model(b,sample=False)
        m=out['label_mask']
        if not m.any(): continue
        logits=out['logits'][m].clone()
        logits[:,[0,2]]=torch.finfo(logits.dtype).min
        ce=F.cross_entropy(logits,b['next_event'][m],reduction='none')
        actual=b['next_dt'][m]
        error=(torch.expm1(out['time_prediction'][m])-actual).abs()*365.25
        masks=[torch.ones_like(actual,dtype=torch.bool),actual<1,actual>=5]
        for i,mask in enumerate(masks):
            totals[i]+=np.array([int(mask.sum()),float(ce[mask].sum()),float(error[mask].sum())])
    rows=[dict(gap=tag,targets=int(n),next_event_ce=ce/n if n else np.nan,
               timing_mae_days=mae/n if n else np.nan) for tag,(n,ce,mae) in zip(['all','under_1y','at_least_5y'],totals)]
    pd.DataFrame(rows).to_csv(outdir(cfg,'main/results')/'next_event.tsv',sep='\t',index=False)
