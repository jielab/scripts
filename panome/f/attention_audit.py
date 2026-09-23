"""Export actual attention tensors and perturbation evidence, without using test Y."""
import copy
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from neural import encode, device_for
from evidence import subset_indices
from borrowing import risk_from_weights, stable_seed
from common import dump, log


def rollout(maps):
    """Mean-head + residual rollout; omits FFN/LN/value mixing, NOT attribution."""
    b, layers, heads, tokens, _ = maps.shape
    joint = np.broadcast_to(np.eye(tokens), (b,tokens,tokens)).copy()
    for level in range(layers):
        flow = maps[:,level].mean(1)+np.eye(tokens)[None]
        flow /= flow.sum(-1,keepdims=True)
        joint = flow@joint
    importance = joint[:,0,1:]
    return importance/np.maximum(importance.sum(1,keepdims=True),1e-30)


def audit_attention(bundle,x,observed,ids,groups,out,a):
    indices = subset_indices(ids,a.attention_samples,a.seed+3)
    if not len(indices): return
    x, observed, ids = x[indices], observed[indices], ids[indices]
    groups = None if groups is None else groups[indices]
    out = Path(out)/'attention'; out.mkdir(exist_ok=True)
    name = bundle['primary']; spec = bundle['specs'][name]
    kw = {k:v for k,v in spec.items() if k not in ['encoder','space']}
    model = bundle['encoders'][spec['encoder']]
    bank = bundle['banks'][name]; cal = bundle['calibrators'][name]
    dev = device_for(a.device); model.to(dev).eval()
    encoded, maps = [], []
    with torch.no_grad():
        for begin in range(0,len(ids),a.batch_size):
            z,_,_,att = model(torch.as_tensor(x[begin:begin+a.batch_size],device=dev),
                torch.as_tensor(observed[begin:begin+a.batch_size],device=dev),return_attention=True)
            encoded.append(z.cpu().numpy()); maps.append(att.cpu().numpy())
    model.cpu()
    z, maps = np.concatenate(encoded), np.concatenate(maps)
    if not np.allclose(maps.sum(-1),1,atol=2e-6):
        raise AssertionError('Attention rows do not sum to one')
    # Capturing the maps must preserve the deployed model output.
    normal_z = encode(model,x,observed,dev,a.batch_size)[0]; model.cpu()
    np.testing.assert_allclose(z,normal_z,rtol=2e-5,atol=2e-6)
    # Use normal execution for every risk; captured maps are evidence only.
    raw,detail = bank.match(normal_z,ids,groups,**kw); baseline = cal.predict(raw)
    tokens = ['CLS']+list(bundle['token_names'])
    importance = rollout(maps)
    np.savez_compressed(out/'self_attention.npz',eid=ids.astype(str),token=np.array(tokens),
        attention=maps,rollout=importance,
        axis_order=np.array(['person','layer','head','query_token','key_token']))
    table = []
    for i,eid in enumerate(ids):
        for level in range(maps.shape[1]):
            for head in range(maps.shape[2]):
                matrix = maps[i,level,head]
                for query in range(len(tokens)):
                    for key in np.argsort(-matrix[query],kind='stable')[:3]:
                        table.append((eid,level+1,head+1,tokens[query],tokens[key],matrix[query,key]))
    pd.DataFrame(table,columns=['eid','layer','head','query_token','key_token','attention_weight']).to_csv(
        out/'self_attention_top_edges.csv.gz',index=False)
    summary = []
    for i,eid in enumerate(ids):
        for level in range(maps.shape[1]):
            for head in range(maps.shape[2]):
                matrix = maps[i,level,head]
                summary.append(dict(eid=eid,layer=level+1,head=head+1,
                    mean_query_entropy=float(-np.sum(matrix*np.log(np.maximum(matrix,1e-30)),-1).mean()),
                    mean_self_weight=float(np.trace(matrix)/len(matrix)),
                    across_query_weight_SD=float(matrix.std(0).mean())))
    pd.DataFrame(summary).to_csv(out/'self_attention_diagnostics.csv',index=False)
    pd.DataFrame(importance,index=ids,columns=tokens[1:]).rename_axis('eid').to_csv(out/'rollout_heuristic.csv')
    experiments = []
    def record(label,changed,changed_detail):
        for i,eid in enumerate(ids):
            experiments.append(dict(eid=eid,intervention=label,original_raw=raw[i],changed_raw=changed[i],
                original_calibrated=baseline[i],changed_calibrated=cal.predict(changed)[i],
                absolute_calibrated_change=abs(baseline[i]-cal.predict(changed)[i]),
                selected_reference_Jaccard=len(set(detail['jj'][i])&set(changed_detail['jj'][i]))/
                    len(set(detail['jj'][i])|set(changed_detail['jj'][i]))))
    # Real head weights reconstruct the marginal donor weight, then raw risk.
    if 'head_weights' in detail:
        heads, gate = detail['head_weights'], detail['head_gate']
        np.testing.assert_allclose((heads*gate[:,None,:]).sum(-1),detail['weights'],atol=2e-7)
        rows = []
        for i,eid in enumerate(ids):
            for rank,j in enumerate(detail['jj'][i]):
                for h in range(gate.shape[1]):
                    coefficient = (1-detail['prior_weight'][i])*gate[i,h]*heads[i,rank,h]
                    rows.append((eid,h+1,rank+1,bank.ids[j],detail['labels'][i,rank],gate[i,h],
                        heads[i,rank,h],coefficient,coefficient*detail['labels'][i,rank]))
        pd.DataFrame(rows,columns=['eid','head','reference_rank','reference_eid','observed_Y','head_gate',
            'within_head_copy_weight','risk_coefficient','raw_risk_contribution']).to_csv(
            out/'reference_head_contributions.csv.gz',index=False)
        head_table = []
        for i,eid in enumerate(ids):
            for h in range(gate.shape[1]):
                head_table.append(dict(eid=eid,head=h+1,gate=gate[i,h],borrowed_head_risk=detail['head_risk'][i,h],
                    effective_donors=1/max(np.sum(heads[i,:,h]**2),1e-30),prior_weight=detail['prior_weight'][i]))
        pd.DataFrame(head_table).to_csv(out/'reference_heads.csv',index=False)
        if gate.shape[1]>1:
            for h in range(gate.shape[1]):
                reduced = gate.copy(); reduced[:,h] = 0; reduced /= reduced.sum(1,keepdims=True)
                weights = (heads*reduced[:,None,:]).sum(-1)
                changed = risk_from_weights(weights,detail['labels'],bank.prior,kw['strength'])[0]
                record(f'drop_cross_head_{h+1}_fixed_neighbors',changed,detail)
    # Uniform/identity attention is an intervention, not an independently retrained comparator.
    old = [layer.intervention for layer in model.context.layers]
    try:
        for mode in ['uniform','identity']:
            for layer in model.context.layers: layer.intervention = mode
            changed_z = encode(model,x,observed,dev,a.batch_size)[0]; model.cpu()
            changed,other = bank.match(changed_z,ids,groups,**kw)
            record(mode+'_query_only_frozen_bank',changed,other)
            changed_bank = copy.copy(bank)
            changed_bank.z = encode(model,bank.x,bank.observed,dev,a.batch_size)[0]; model.cpu()
            changed,other = changed_bank.match(changed_z,ids,groups,**kw)
            record(mode+'_query_and_bank',changed,other)
    finally:
        for layer, value in zip(model.context.layers,old): layer.intervention = value
        model.cpu()
    # Does high attention predict sensitivity better than masking random modules?
    n_mask = min(2,len(tokens)-1)
    sets = {'top_rollout':np.argsort(-importance,axis=1)[:,:n_mask],
            'bottom_rollout':np.argsort(importance,axis=1)[:,:n_mask],
            'random_modules':np.array([np.random.default_rng(stable_seed(eid,a.seed+931)).choice(
                len(tokens)-1,n_mask,replace=False) for eid in ids])}
    for label, selected in sets.items():
        mask = observed.copy()
        for i in range(len(ids)): mask[i,np.isin(bundle['membership'],selected[i])] = False
        zz = encode(model,x,mask,dev,a.batch_size)[0]; model.cpu()
        changed,other = bank.match(zz,ids,groups,**kw)
        record('mask_'+label,changed,other)
    frame = pd.DataFrame(experiments)
    frame.to_csv(out/'attention_interventions.csv.gz',index=False)
    frame.groupby('intervention').agg(n=('eid','size'),mean_absolute_change=('absolute_calibrated_change','mean'),
        median_absolute_change=('absolute_calibrated_change','median'),
        mean_match_Jaccard=('selected_reference_Jaccard','mean')).to_csv(out/'intervention_summary.csv')
    # All-module perturbations from the existing module audit provide a stronger
    # person-wise correspondence check; no manufactured correlation if constant.
    module_file = Path(out).parent/'module_perturbations.csv.gz'
    correspondence = []
    if module_file.exists():
        module = pd.read_csv(module_file,dtype={'eid':str})
        lookup = {str(eid):i for i,eid in enumerate(ids)}
        for eid,sub in module.groupby('eid'):
            if eid not in lookup: continue
            ordered = sub.set_index('token').reindex(tokens[1:])
            sensitivity = ordered.delta_original_minus_masked.abs().to_numpy()
            score = importance[lookup[eid]]
            ok = np.isfinite(sensitivity)
            corr = spearmanr(score[ok],sensitivity[ok]).statistic if ok.sum()>2 and np.std(sensitivity[ok])>1e-12 and np.std(score[ok])>1e-12 else np.nan
            correspondence.append(dict(eid=eid,rollout_vs_masking_spearman=corr,n_tokens=int(ok.sum())))
    pd.DataFrame(correspondence,columns=['eid','rollout_vs_masking_spearman','n_tokens']).to_csv(out/'attention_vs_sensitivity.csv',index=False)
    dump(out/'audit.json',dict(n_people=len(ids),selection='outcome-blind stable ID hash',
        self_attention_shape=maps.shape,axis_order=['person','layer','head','query_token','key_token'],
        maps='actual pre-dropout Q/K softmax in eval mode; exported for every layer/head',
        cross_values='observed build donor outcomes, IPCW-adjusted within each head',
        caveats=['Attention and rollout are not causal effects or complete feature attributions.',
        'Interventions keep fitted parameters/calibration frozen; may move inputs out of distribution.',
        'Module count is matched in masking controls, not number of assays; modules differ in size.',
        'Retrained uniform/metric controls in main comparison answer different questions from these interventions.']))
    log('DONE','attention_audit',f'n={len(ids)}, layers={maps.shape[1]}, heads={maps.shape[2]}')
