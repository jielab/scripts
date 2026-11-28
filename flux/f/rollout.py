"""Event/time sampling with state feedback, no observed post-landmark inputs."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import torch
from torch.nn import functional as F
from common import atomic_checkpoint, digest, file_digest, log, outdir, seed_all, to_device, torch_load
from data import CohortDataset
from train import load_world, make_loader, risk_cache
from evaluate import score_matrix


@torch.no_grad()
def calibrate_generation(cfg,model,device,*,replace=False):
    """Validation-only temperature and empirical log-interval residuals."""
    path=outdir(cfg,'main/results/models')/'generation.pt'; rc=cfg['rollout']
    key=digest(dict(world=file_digest(Path(cfg['paths']['out'])/'main/results/models/flux.pt'),
        batches=rc['calibration_batches'],events=rc['max_calibration_events']))
    if path.exists() and not replace:
        result=torch_load(path)
        if result['key']==key: return result
    dataset=CohortDataset(cfg,'validation',mode='pretrain',max_nodes=cfg['data']['max_nodes'],compact_numeric=True)
    loader=make_loader(dataset,cfg['train']['batch_size'],cfg['train']['workers'])
    logits=[]; labels=[]; residuals=[]; count=0; cap=int(rc['max_calibration_events'])
    for i,batch in enumerate(loader):
        if rc['calibration_batches'] and i>=int(rc['calibration_batches']): break
        b=to_device(batch,device); result=model(b,sample=False)
        mask=result['label_mask']; remaining=cap-count
        if not mask.any(): continue
        x=result['logits'][mask][:remaining].float().cpu()
        y=b['next_event'][mask][:remaining].cpu()
        residual=(torch.log1p(b['next_dt'][mask])-result['time_prediction'][mask])[:remaining].float().cpu()
        logits.append(x); labels.append(y); residuals.append(residual); count+=len(y)
        if count>=cap: break
    if not count: raise ValueError('No validation next-event labels for generation calibration')
    x=torch.cat(logits); y=torch.cat(labels); residual=torch.cat(residuals).numpy()
    x[:,0]=-1e9; x[:,2]=-1e9
    grid=[0.5,0.75,1.0,1.5,2.0,3.0]
    losses=[float(F.cross_entropy(x/t,y)) for t in grid]
    temperature=grid[int(np.argmin(losses))]
    result=dict(key=key,temperature=temperature,time_residuals=residual,validation_n=count,
        temperature_grid=grid,validation_ce=losses,
        note='Independent validation grid and empirical residual scheme; not official calibration parameters.')
    atomic_checkpoint(result,path)
    log(f'Generation calibration: validation n={count:,}; temperature={temperature}')
    return result


def _repeat_state(state,particles):
    return {k:v.repeat_interleave(particles,dim=0) for k,v in state.items()}


def _remove_horizon_evaluation(folder, horizon):
    """Remove only this rollout horizon's stale aggregate evaluation files."""
    prefix = f'rollout.{horizon:g}y'
    for suffix in ('metrics.tsv', 'macro.tsv', 'calibration.tsv', 'provenance.json'):
        (Path(folder) / f'{prefix}.{suffix}').unlink(missing_ok=True)


@torch.no_grad()
def rollout(cfg,requested_device='auto',*,replace=False):
    model,meta,device=load_world(cfg,requested_device)
    seed_all(int(cfg['seed'])); rc=cfg['rollout']
    calibration=calibrate_generation(cfg,model,device,replace=replace)
    risk,manifest=risk_cache(cfg,'test')
    dataset=CohortDataset(cfg,'test',mode='risk',horizon_years=cfg['readout']['horizon_years'],compact_numeric=True)
    from train import state_cache_key
    if manifest['key'] != state_cache_key(cfg):
        raise ValueError('Rollout needs the test risk cache from the current world model')
    names=dataset.endpoint_names; particles=int(rc['particles']); horizons=sorted(float(x) for x in rc['horizons'])
    if particles<1 or not horizons or not np.isfinite(horizons).all() or min(horizons)<=0 or (rc['max_events'] is not None and int(rc['max_events'])<1):
        raise ValueError('rollout particles, horizons and max_events must be positive and finite')
    if len({f'{h:g}' for h in horizons}) != len(horizons):
        raise ValueError('rollout horizons must be unique at output-filename precision')
    if names != manifest['endpoint_names']: raise ValueError('Rollout endpoint order differs from saved risk rows')
    maximum=max(horizons); folder=outdir(cfg,'main/results')
    matrices={h:np.lib.format.open_memmap(folder/f'rollout.{h:g}y.test.npy',mode='w+',
        dtype='float32',shape=(len(dataset),len(names))) for h in horizons}
    # One vocabulary token may belong to several explicitly grouped endpoints.
    mapping=np.zeros((len(meta['vocab']),len(names)),bool)
    for v,code in enumerate(meta['vocab']):
        for j,name in enumerate(names):
            mapping[v,j]=any(code==prefix if prefix=='DEATH' else code.startswith(prefix[:3])
                             for prefix in meta['endpoints'][name])
    # The endpoint probabilities at each event time remain continuous even for
    # diseases rarely drawn in eight trajectories. Aggregate mutually exclusive
    # token probabilities for grouped endpoints.
    token_idx, end_idx = np.nonzero(mapping)
    token_idx=torch.as_tensor(token_idx,device=device)
    end_idx=torch.as_tensor(end_idx,device=device)
    boundaries=torch.tensor(horizons,device=device)
    
    supported=torch.as_tensor(mapping.any(axis=0),device=device)
    unsupported_endpoints=[name for j,name in enumerate(names) if not mapping[:,j].any()]
    horizon_status={f'{h:g}':dict(horizon_years=h,completed_trajectories=0,
        incomplete_trajectories=0,completed_participants=0,affected_participants=0,
        status='pending',performance_reported=False) for h in horizons}
    for matrix in matrices.values(): matrix[:]=np.nan
    death_token=meta['vocab'].index('DEATH')
    residual=torch.as_tensor(calibration['time_residuals'],device=device)
    offset=0; truncated=dead=0; simulated_events=0; unknown_events=0
    loader=make_loader(dataset,rc['batch_size'],cfg['train']['workers'])
    minimum=float(rc['min_interval_days'])/365.25
    if not np.isfinite(minimum) or minimum<=0: raise ValueError('min_interval_days must be positive and finite')
    for batch in loader:
        n=len(batch['eid'])
        if batch['eid']!=list(risk['eid'][offset:offset+n]): raise ValueError('Rollout risk-set rows differ from saved readout rows')
        initial=model.encode_history(to_device(batch,device),sample=False,chunk_size=cfg['readout']['history_chunk_size'])
        state=_repeat_state(initial,particles)
        total=n*particles
        elapsed=torch.zeros(total,device=device)
        active=torch.ones(total,dtype=torch.bool,device=device)
        log_survival=torch.zeros((len(horizons),total,len(names)),device=device)
        particle_ids=torch.arange(total,device=device)
        event_steps=0
        while active.any() and (rc['max_events'] is None or event_steps < rc['max_events']):
            event_steps+=1
            logits,predicted=model.predict_next(state)
            logits=logits.float()/float(calibration['temperature']); logits[:,0]=-1e9; logits[:,2]=-1e9
            prob=torch.softmax(logits,dim=-1)
            token=torch.multinomial(prob,1).squeeze(-1)
            noise=residual[torch.randint(len(residual),(total,),device=device)]
            interval=torch.expm1((predicted+noise).clamp(min=float(np.log1p(minimum)),max=float(np.log1p(maximum+1))))
            arrives=active&(elapsed+interval<=maximum)
            jump=torch.where(arrives,interval,torch.zeros_like(interval))
            state=model.advance_state(state,jump,sample=True)
            elapsed=elapsed+jump
            emitted=torch.where(arrives,token,torch.zeros_like(token))
            state=model.observe_event(state,emitted,sample=True)
            # At each arrival multiply conditional survival by (1 - p_endpoint).
            # The sampled token/time updates the subsequent trajectory and death stops it.
            p=torch.zeros((total,len(names)),device=device)
            p.index_add_(1,end_idx,prob.index_select(1,token_idx))
            bucket=torch.searchsorted(boundaries,elapsed[arrives])
            log_survival[bucket,particle_ids[arrives]] += torch.log1p(-p[arrives].clamp(0,1))
            died=arrives&(token==death_token)
            dead+=int(died.sum()); simulated_events+=int(arrives.sum())
            unknown_events+=int((arrives&(token==1)).sum())
            active=arrives&~died
        truncated+=int((active&(elapsed<maximum)).sum())
        cumulative=-torch.expm1(log_survival.cumsum(dim=0))
        for k,(h,matrix) in enumerate(matrices.items()):
            complete=((~active)|(elapsed>=h)).reshape(n,particles)
            good=complete.all(dim=1)
            scores=cumulative[k].reshape(n,particles,-1).mean(dim=1)
            scores[~good]=float('nan')
            scores[:,~supported]=float('nan')
            matrix[offset:offset+n]=scores.cpu().numpy()
            summary=horizon_status[f'{h:g}']
            summary['completed_trajectories']+=int(complete.sum())
            summary['incomplete_trajectories']+=int((~complete).sum())
            summary['completed_participants']+=int(good.sum())
            summary['affected_participants']+=int((~good).sum())
        offset+=n
        if offset%1000<n: log(f'Rollout: {offset:,}/{len(dataset):,}; {particles} particles each')
    if offset!=len(dataset): raise ValueError('Rollout loader did not cover every saved risk row')
    for value in matrices.values(): value.flush()
    for summary in horizon_status.values():
        if summary['incomplete_trajectories']:
            summary['status']='insufficient_simulation_budget'
        elif len(unsupported_endpoints)==len(names):
            summary['status']='no_supported_endpoints'
        else:
            summary['status']='complete'
    info=dict(particles=particles,horizons=horizons,max_events=rc['max_events'],
        trajectories=len(dataset)*particles,truncated_trajectories=truncated,
        absorbing_deaths=dead,simulated_events=simulated_events,unknown_events=unknown_events,
        horizon_status=horizon_status,unsupported_endpoints=unsupported_endpoints,
        supported_endpoint_count=len(names)-len(unsupported_endpoints),
        unsupported_endpoint_reason='No generatable vocabulary token maps to this endpoint.',
        truncation_policy='Incomplete person/horizon scores are NaN; cohort performance is withheld whenever any trajectory is incomplete at that horizon. Completed shorter horizons remain eligible.',
        calibration_validation_n=calibration['validation_n'],temperature=calibration['temperature'],
        no_observed_post_landmark_input=True,risk_key=manifest["key"],config=rc,
        score_definition="Mean across particles of 1 minus product of (1 minus conditional endpoint probability) at generated arrivals through horizon; death is absorbing.",
        interpretation='Particle-average conditional-incidence ranking scores; not calibrated lifetime or treatment-effect probabilities.')
    (folder/'rollout.json').write_text(json.dumps(info,indent=2)+'\n')
    for h,matrix in matrices.items():
        summary=horizon_status[f'{h:g}']
        if summary['status']!='complete':
            _remove_horizon_evaluation(folder,h)
            log(f'Rollout {h:g}y performance withheld: {summary["status"]}; '
                f'{summary["incomplete_trajectories"]} incomplete trajectories, '
                f'{summary["affected_participants"]} affected participants.')
            continue
        # Deliberately supply no probability horizon: these are ranking scores.
        score_matrix(cfg,'rollout',matrix,horizons=[h],suffix=f'rollout.{h:g}y')
        summary['performance_reported']=True
    (folder/'rollout.json').write_text(json.dumps(info,indent=2)+'\n')
    if truncated: log(f'Rollout limit reached in {truncated}/{len(dataset)*particles} trajectories; increase max_events to resolve unsupported horizons.')
    log('Rollout complete; generated events feed back into the state, death is absorbing.')
