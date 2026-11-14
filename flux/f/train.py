"""World-model fitting, complete-prefix state extraction and frozen risk readout."""
from __future__ import annotations
import functools
from datetime import timedelta
import json
import math
import os
from pathlib import Path
import random
import time

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F
from torch.utils.data import DataLoader
from common import (atomic_checkpoint, choose_device, digest, file_digest, log,
                    outdir, seed_all, to_device, torch_load)
from data import CohortDataset, collate_batch, load_meta
from model import FluxModel, RiskHead


def model_spec(cfg, meta):
    return dict(model_cfg=cfg['model'], n_tokens=len(meta['vocab']),
                context_dim=len(meta['context_features']),
                modality_dims={k: len(v) for k, v in meta['modalities'].items()})


def data_signature(meta):
    return digest({k: meta.get(k) for k in ('schema_version','vocab','features','numeric_scaler',
        'feature_modalities','modalities','context_features','context_scaler','endpoints',
        'split_counts','sources','admin_end','seed','normalization_clip','event_source_mode',
        'same_day_order')})


def resume_contract(cfg):
    """Settings that affect the saved minibatch cursor or optimization trajectory.

    epochs/max_steps/patience may be increased to continue a run. Changing the
    actual training algorithm requires a new run rather than reusing its RNG
    state and optimizer under different semantics.
    """
    keys = ('batch_size','accumulation_steps','learning_rate','weight_decay',
            'gradient_clip','validate_every','validation_batches','workers')
    return dict(seed=cfg['seed'], max_nodes=cfg['data']['max_nodes'],
                cores=cfg['resources']['cores'], train={k: cfg['train'][k] for k in keys},
                schedule={k:cfg['train'].get(k,v) for k,v in
                          (('warmup_steps',0),('schedule_steps',0),('final_lr_ratio',1.0))})


def scheduled_learning_rate(tc, optimizer_step):
    """Absolute-step warmup/cosine schedule, independent of a resumed run's stop."""
    base=float(tc['learning_rate']); warmup=int(tc.get('warmup_steps',0))
    total=int(tc.get('schedule_steps',0)); ratio=float(tc.get('final_lr_ratio',1.0))
    if warmup<0 or total<0 or not 0<=ratio<=1 or (total and total<=warmup):
        raise ValueError('Invalid warmup/cosine schedule: schedule_steps must exceed warmup_steps, final_lr_ratio in [0,1]')
    if warmup and optimizer_step<=warmup: return base*optimizer_step/warmup
    if not total: return base
    progress=min(1.,max(0.,(optimizer_step-warmup)/(total-warmup)))
    return base*(ratio+(1-ratio)*0.5*(1+math.cos(math.pi*progress)))


def make_loader(dataset, batch_size, workers=0, shuffle=False, sampler=None, seed=2026):
    generator = torch.Generator().manual_seed(seed)
    return DataLoader(dataset, batch_size=int(batch_size), shuffle=shuffle and sampler is None,
        sampler=sampler, num_workers=int(workers), persistent_workers=False,
        pin_memory=torch.cuda.is_available(), collate_fn=functools.partial(collate_batch, meta=dataset.meta),
        generator=generator)


def rng_state():
    return dict(python=random.getstate(), numpy=np.random.get_state(), torch=torch.get_rng_state(),
                cuda=torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None)


def restore_rng(value):
    random.setstate(value['python']); np.random.set_state(value['numpy']); torch.set_rng_state(value['torch'])
    if value.get('cuda') is not None and torch.cuda.is_available(): torch.cuda.set_rng_state_all(value['cuda'])


@torch.no_grad()
def validation_loss(model, loader, device, max_batches, *, distributed=False):
    was_training = model.training
    model.eval()
    totals = dict(event_sum=0., time_sum=0., kl_balanced_sum=0., event_count=0., kl_count=0.)
    try:
        for index, batch in enumerate(loader):
            if max_batches and index >= max_batches: break
            result = model(to_device(batch, device), sample=False, return_sequence=False)
            for key in totals: totals[key] += float(result[key])
    finally:
        model.train(was_training)
    if distributed:
        import torch.distributed as dist
        packed=torch.tensor(list(totals.values()),dtype=torch.float64,device=device)
        dist.all_reduce(packed,op=dist.ReduceOp.SUM)
        totals=dict(zip(totals,packed.cpu().tolist()))
    if totals['event_count'] == 0:
        raise ValueError('No next-event targets in validation. Check event dates, vocabulary and cohort size.')
    cfg = model.model_cfg
    loss = (totals['event_sum'] + cfg['time_weight']*totals['time_sum'])/totals['event_count']
    loss += cfg['kl_weight']*totals['kl_balanced_sum']/max(totals['kl_count'],1)
    if not math.isfinite(loss): raise RuntimeError('Non-finite world-model validation loss')
    return loss


def train_world(cfg, requested_device='auto', *, resume=False, replace=False):
    if resume and replace: raise ValueError('--resume and --replace cannot be combined for training')
    device = choose_device(requested_device)
    torch.set_num_threads(int(cfg['resources']['cores']))
    world = int(os.getenv('WORLD_SIZE','1')); rank = int(os.getenv('RANK','0'))
    distributed = world > 1
    if distributed:
        import torch.distributed as dist
        local_rank = int(os.environ['LOCAL_RANK'])
        torch.cuda.set_device(local_rank); device = torch.device('cuda',local_rank)
        dist.init_process_group('nccl',timeout=timedelta(minutes=60))
    seed_all(int(cfg['seed']) + rank)
    meta = load_meta(cfg); signature = data_signature(meta)
    spec = model_spec(cfg, meta)
    folder = outdir(cfg, 'model'); best_file = folder/'flux.pt'; last_file = folder/'flux.last.pt'
    if best_file.exists() and not replace and not resume:
        saved = torch_load(best_file)
        if saved['data_signature'] != signature or saved['spec'] != spec:
            raise ValueError('Existing world model uses different data/model settings. Use a different --dir-out or --replace.')
        previous=resume_contract(saved['config'])
        # CPU thread count does not invalidate reuse of an already fitted model.
        current=resume_contract(cfg); previous.pop('cores',None); current.pop('cores',None)
        budgets=('max_steps','epochs','patience')
        if previous!=current or any(saved['config']['train'][k]!=cfg['train'][k] for k in budgets):
            raise ValueError('Training settings changed. Use --resume for a compatible longer run, or a new --dir-out / --replace for a different experiment.')
        if rank == 0: log('Reusing model/flux.pt; --resume continues training, --replace retrains.')
        if distributed: dist.destroy_process_group()
        return
    base = FluxModel(**spec).to(device)
    tc = cfg['train']
    scheduled_learning_rate(tc,1)
    optimizer = torch.optim.AdamW(base.parameters(), lr=float(tc['learning_rate']), weight_decay=float(tc['weight_decay']))
    start_epoch = start_batch = step = 0
    best = float('inf'); stale = 0
    if resume:
        saved = torch_load(last_file)
        if saved['data_signature'] != signature or saved['spec'] != spec:
            raise ValueError('Cannot resume with different input schema or model configuration')
        if saved.get('world_size',1) != world:
            raise ValueError('Exact resume requires the original GPU count; use --replace to retrain at another count')
        previous_contract = saved.get('resume_contract')
        if previous_contract is None and 'config' in saved:
            previous_contract = resume_contract(saved['config'])
        if previous_contract != resume_contract(cfg):
            raise ValueError('Cannot exactly resume with changed seed, windows, minibatches or optimizer settings; use --replace or a new output directory')
        base.load_state_dict(saved['model']); optimizer.load_state_dict(saved['optimizer'])
        for state in optimizer.state.values():
            for key, val in state.items():
                if isinstance(val, torch.Tensor): state[key] = val.to(device)
        start_epoch, start_batch, step = saved['epoch'], saved['next_batch'], saved['step']
        best, stale = saved['best'], saved['stale']
        restore_rng(saved['rng'][rank])
    wrapped = torch.nn.parallel.DistributedDataParallel(base, device_ids=[device.index], find_unused_parameters=True) if distributed else base
    train = CohortDataset(cfg, 'train', mode='pretrain', max_nodes=cfg['data']['max_nodes'])
    validation = CohortDataset(cfg, 'validation', mode='pretrain', max_nodes=cfg['data']['max_nodes'])
    if not len(train) or not len(validation): raise ValueError('Training or validation cohort is empty')
    sampler = torch.utils.data.distributed.DistributedSampler(train, num_replicas=world, rank=rank, seed=int(cfg['seed'])) if distributed else None
    validation_limit=len(validation)
    if int(tc['validation_batches']):
        validation_limit=min(validation_limit,int(tc['validation_batches'])*int(tc['batch_size']))
    # Non-padding shards: every selected validation person is counted exactly once.
    valid_loader = make_loader(validation,tc['batch_size'],tc['workers'],
        sampler=range(rank,validation_limit,world),seed=int(cfg['seed']))
    max_steps = int(tc['max_steps']); accumulation = int(tc['accumulation_steps'])
    if max_steps < 1 or accumulation < 1: raise ValueError('max_steps and accumulation_steps must be positive')
    if min(int(tc['validate_every']), int(tc['epochs']), int(tc['patience'])) < 1:
        raise ValueError('validate_every, epochs and patience must be positive')
    if rank == 0:
        log(f'World model: {len(train):,} train participants; {sum(p.numel() for p in base.parameters()):,} parameters; {device}; effective batch {int(tc["batch_size"])*accumulation*world}')
        log('Recurrent state, attained age and RK4 integration use float32 to preserve small time increments.')
    history = outdir(cfg,'logs')/'train.tsv'
    if rank == 0 and (replace or not history.exists()): history.write_text('step\tepoch\ttrain_loss\tvalidation_loss\tseconds\n')
    tick = time.time(); stop = step >= max_steps or stale >= int(tc['patience']); recent = []
    for epoch in range(start_epoch, int(tc['epochs'])):
        if stop: break
        train.set_epoch(epoch)
        if sampler is not None: sampler.set_epoch(epoch)
        loader = make_loader(train, tc['batch_size'], tc['workers'], shuffle=True, sampler=sampler,
                             seed=int(cfg['seed'])+epoch)
        optimizer.zero_grad(set_to_none=True); base.train()
        for batch_index, batch in enumerate(loader):
            if epoch == start_epoch and batch_index < start_batch: continue
            batch = to_device(batch, device)
            result = wrapped(batch, sample=True, return_sequence=False)
            loss = result['loss']
            if not bool(torch.isfinite(loss)): raise RuntimeError(f'Non-finite world-model loss at optimizer step {step}')
            # The last partial group is normalized by its actual microbatch count.
            group_start = (batch_index//accumulation)*accumulation
            group_size = min(accumulation, len(loader)-group_start)
            (loss/group_size).backward()
            recent.append(float(loss.detach()))
            update = (batch_index+1)%accumulation == 0 or batch_index+1 == len(loader)
            if not update: continue
            current_lr=scheduled_learning_rate(tc,step+1)
            for group in optimizer.param_groups: group['lr']=current_lr
            nn.utils.clip_grad_norm_(base.parameters(), float(tc['gradient_clip']))
            optimizer.step(); optimizer.zero_grad(set_to_none=True); step += 1
            validate = (step%int(tc['validate_every']) == 0 or step >= max_steps or
                        (epoch+1==int(tc['epochs']) and batch_index+1==len(loader)))
            if validate:
                if distributed: dist.barrier()
                val = validation_loss(base,valid_loader,device,0,distributed=distributed)
                improved = val < best
                best, stale = (val,0) if improved else (best,stale+1)
                states = [None]*world
                if distributed: dist.all_gather_object(states, rng_state())
                else: states[0] = rng_state()
                saved = dict(model=base.state_dict(), optimizer=optimizer.state_dict(), spec=spec,
                    data_signature=signature, step=step, epoch=epoch, next_batch=batch_index+1,
                    best=best, stale=stale, rng=states, world_size=world, config=cfg,
                    resume_contract=resume_contract(cfg))
                if rank == 0:
                    atomic_checkpoint(saved, last_file)
                    if improved: atomic_checkpoint(saved, best_file)
                    average = float(np.mean(recent))
                    with history.open('a') as handle: handle.write(f'{step}\t{epoch}\t{average:.8g}\t{val:.8g}\t{time.time()-tick:.2f}\n')
                    log(f'train step {step}/{max_steps}; loss={average:.4f}; validation={val:.4f}; best={best:.4f}')
                recent = []
                if distributed: dist.barrier()
            stop = step >= max_steps or stale >= int(tc['patience'])
            if stop: break
        start_batch = 0
    if not best_file.exists(): raise RuntimeError('No checkpoint produced; increase train.epochs or inspect next-event targets')
    if distributed: dist.destroy_process_group()


def load_world(cfg, requested_device='auto'):
    device = choose_device(requested_device)
    torch.set_num_threads(int(cfg['resources']['cores']))
    saved = torch_load(Path(cfg['paths']['out'])/'model/flux.pt')
    meta = load_meta(cfg)
    if saved['data_signature'] != data_signature(meta): raise ValueError('World-model checkpoint does not match the prepared data')
    model = FluxModel(**saved['spec']).to(device)
    model.load_state_dict(saved['model']); model.eval()
    return model, meta, device


def risk_cache(cfg, split):
    folder = Path(cfg['paths']['out'])/'data'/f'risk_{split}'
    manifest = json.loads((folder/'index.json').read_text())
    arrays = {name: np.load(folder/f'{name}.npy', mmap_mode='r') for name in manifest['arrays']}
    return arrays, manifest


def state_cache_key(cfg):
    return digest(dict(model=file_digest(Path(cfg['paths']['out'])/'model/flux.pt'),
        data=cfg['data'], horizon=cfg['readout']['horizon_years']))


@torch.no_grad()
def extract_states(cfg, requested_device='auto', *, replace=False):
    model, meta, device = load_world(cfg, requested_device)
    key = state_cache_key(cfg)
    for split in ('train','validation','test'):
        folder = outdir(cfg, f'data/risk_{split}')
        if (folder/'index.json').exists() and not replace:
            previous = json.loads((folder/'index.json').read_text())
            if previous['key'] == key: continue
        dataset = CohortDataset(cfg,split,mode='risk',horizon_years=cfg['readout']['horizon_years'])
        n, d = len(dataset), len(dataset.endpoint_names)
        if not n: raise ValueError(f'No eligible {split} risk landmarks; inspect min_clinical_days and landmark settings')
        (folder/'index.json').unlink(missing_ok=True)
        shapes = {'features':((n,model.state_dim),'float32'), 'event_time':((n,d),'float32'),
                  'risk_mask':((n,d),'bool'), 'exit_time':((n,),'float32'), 'death_time':((n,),'float32'),
                  'landmark_day':((n,),'int32'), 'eid':((n,),'U32')}
        arrays = {name: np.lib.format.open_memmap(folder/f'{name}.npy',mode='w+',dtype=dtype,shape=shape)
                  for name,(shape,dtype) in shapes.items()}
        loader = make_loader(dataset,cfg['readout']['state_batch_size'],cfg['train']['workers'])
        offset = 0
        for batch in loader:
            count = len(batch['eid'])
            state = model.encode_history(to_device(batch,device),sample=False,
                                         chunk_size=cfg['readout']['history_chunk_size'])
            arrays['features'][offset:offset+count] = model.state_features(state).cpu().numpy()
            arrays['eid'][offset:offset+count] = batch['eid']
            for source,target in (('event_times','event_time'),('risk_mask','risk_mask'),('exit_time','exit_time'),
                                 ('death_time','death_time'),('landmark_day','landmark_day')):
                arrays[target][offset:offset+count] = batch[source].numpy()
            offset += count
            if offset%10000 < count: log(f'Frozen states {split}: {offset:,}/{n:,}')
        for value in arrays.values(): value.flush()
        manifest = dict(key=key,n=n,endpoint_names=dataset.endpoint_names,arrays=list(arrays),
                        horizon_years=cfg['readout']['horizon_years'],complete_history=True)
        (folder/'index.json').write_text(json.dumps(manifest,indent=2)+'\n')
        log(f'Frozen states {split}: {n:,}; every observation at/before each landmark included')


def _state_moments(array, block=8192):
    total = np.zeros(array.shape[1],np.float64); squares = total.copy(); n = len(array)
    for i in range(0,n,block):
        x = np.asarray(array[i:i+block],dtype=np.float64)
        total += x.sum(0); squares += np.square(x).sum(0)
    mean = total/max(n,1)
    std = np.sqrt(np.maximum(squares/max(n,1)-mean**2,1e-8))
    return mean.astype('float32'), std.astype('float32')


def _outcome_counts(arrays, horizon):
    d = arrays['event_time'].shape[1]
    cases = np.zeros(d,np.int64); counts=cases.copy()
    for i in range(0,len(arrays['features']),8192):
        mask = arrays['risk_mask'][i:i+8192]
        cases += ((arrays['event_time'][i:i+8192] <= horizon)&mask).sum(0)
        counts += mask.sum(0)
    return cases, counts


def fit_readout(cfg, requested_device='auto', *, replace=False):
    extract_states(cfg,requested_device,replace=replace)
    device = choose_device(requested_device); seed_all(int(cfg['seed']))
    train, tm = risk_cache(cfg,'train'); valid, vm = risk_cache(cfg,'validation')
    rc = cfg['readout']; horizon=float(rc['horizon_years'])
    if min(int(rc[k]) for k in ('batch_size','epochs','patience','min_cases')) < 1:
        raise ValueError('Readout batch_size, epochs, patience and min_cases must be positive')
    if not math.isfinite(horizon) or horizon <= 0 or float(rc['positive_weight_cap']) < 1:
        raise ValueError('Readout horizon must be positive and positive_weight_cap must be at least 1')
    key = digest(dict(states=tm['key'],readout=rc,seed=cfg['seed'])); target = outdir(cfg,'model')/'readout.pt'
    if target.exists() and not replace:
        saved = torch_load(target)
        if saved.get('key') != key: raise ValueError('Readout configuration changed. Use --replace or a new output directory.')
        log('Reusing model/readout.pt')
    else:
        mean,std = _state_moments(train['features'])
        positive,total = _outcome_counts(train,horizon)
        eligible = (positive >= int(rc['min_cases'])) & (total-positive >= int(rc['min_cases']))
        if not eligible.any(): raise ValueError('No endpoints have enough training cases and controls for the risk readout')
        ratio=np.clip((total-positive)/np.maximum(positive,1),1,float(rc['positive_weight_cap']))
        d = len(positive)
        head=RiskHead(len(mean),d,hidden_dim=rc['hidden_dims'],dropout=rc['dropout']).to(device)
        prevalence=np.clip(positive/np.maximum(total,1),1e-5,1-1e-5)
        with torch.no_grad(): head.net[-1].bias.copy_(torch.as_tensor(np.log(prevalence/(1-prevalence)),dtype=torch.float32,device=device))
        optimizer=torch.optim.AdamW(head.parameters(),lr=float(rc['learning_rate']),weight_decay=float(rc['weight_decay']))
        pos_weight=torch.as_tensor(ratio,dtype=torch.float32,device=device)
        active=torch.as_tensor(eligible,device=device)
        rng=np.random.default_rng(int(cfg['seed'])); batch_size=int(rc['batch_size'])
        best=float('inf'); stale=0
        history=outdir(cfg,'logs')/'readout.tsv'; history.write_text('epoch\ttrain_loss\tvalidation_loss\n')

        def batches(arrays,indices,training):
            loss_sums=torch.zeros(d,dtype=torch.float64,device=device)
            counts_sum=torch.zeros(d,dtype=torch.long,device=device)
            for start in range(0,len(indices),batch_size):
                ix=indices[start:start+batch_size]
                x=torch.as_tensor((arrays['features'][ix]-mean)/std,device=device)
                y=torch.as_tensor((arrays['event_time'][ix]<=horizon).astype('float32'),device=device)
                mask=torch.as_tensor(np.array(arrays['risk_mask'][ix]),device=device)&active
                logits=head(x)
                loss_element=F.binary_cross_entropy_with_logits(logits,y,pos_weight=pos_weight,reduction='none')
                # Endpoint-balanced mean: each supported endpoint has equal weight.
                counts=mask.sum(0); support=counts>0
                if not bool(support.any()): continue
                loss=((loss_element*mask).sum(0)[support]/counts[support]).mean()
                if not bool(torch.isfinite(loss)):
                    raise RuntimeError('Non-finite readout loss; inspect cached features and normalization')
                if training:
                    optimizer.zero_grad(set_to_none=True); loss.backward()
                    nn.utils.clip_grad_norm_(head.parameters(),1.0); optimizer.step()
                loss_sums+=(loss_element.detach()*mask).sum(0).to(torch.float64)
                counts_sum+=counts
            supported=counts_sum>0
            if not bool(supported.any()):
                raise ValueError('No eligible readout labels in this split after training-endpoint filtering')
            # Aggregate once across people: validation must not depend on where
            # batches happened to divide an endpoint's eligible participants.
            return float((loss_sums[supported]/counts_sum[supported]).mean())

        for epoch in range(int(rc['epochs'])):
            head.train(); train_loss=batches(train,rng.permutation(len(train['features'])),True)
            head.eval()
            with torch.no_grad(): val_loss=batches(valid,np.arange(len(valid['features'])),False)
            with history.open('a') as handle: handle.write(f'{epoch}\t{train_loss:.8g}\t{val_loss:.8g}\n')
            log(f'readout epoch {epoch+1}; train={train_loss:.4f}; validation={val_loss:.4f}; endpoints={eligible.sum()}/{d}')
            if val_loss<best:
                best=val_loss; stale=0
                saved=dict(model=head.state_dict(),state_mean=mean,state_std=std,eligible=eligible,
                    config=rc,endpoint_names=tm['endpoint_names'],key=key,states_key=tm['key'],positive_weight=ratio,
                    epoch=epoch,validation_loss=val_loss,horizon_years=horizon,
                    estimand='recorded_incidence',
                    label_policy='No recorded endpoint event within horizon is a control; short follow-up is retained')
                atomic_checkpoint(saved,target)
            else: stale+=1
            if stale>=int(rc['patience']): break
        saved=torch_load(target)
    head=RiskHead(len(saved['state_mean']),len(saved['endpoint_names']),hidden_dim=saved['config']['hidden_dims'],dropout=saved['config']['dropout']).to(device)
    head.load_state_dict(saved['model']); head.eval()
    folder=outdir(cfg,'results')
    for split in ('validation','test'):
        arrays,manifest=risk_cache(cfg,split)
        pred=np.lib.format.open_memmap(folder/f'flux.{split}.npy',mode='w+',dtype='float32',
                                      shape=arrays['event_time'].shape)
        with torch.no_grad():
            for start in range(0,len(pred),int(rc['batch_size'])):
                x=(arrays['features'][start:start+int(rc['batch_size'])]-saved['state_mean'])/saved['state_std']
                prob=torch.sigmoid(head(torch.as_tensor(x,device=device))).cpu().numpy()
                prob[:,~saved['eligible']]=np.nan
                pred[start:start+len(prob)]=prob
        pred.flush()
    (folder/'flux.prediction.json').write_text(json.dumps(dict(endpoint_names=saved['endpoint_names'],
        states_key=saved['states_key'],readout_key=saved['key'],readout_sha256=file_digest(target),
        horizon_years=horizon,score_kind='weighted_bce_probability',eligible=saved['eligible'].tolist(),
        note='Positive-weighted BCE outputs are not assumed calibrated. Brier/calibration are reported without post-hoc correction.'),indent=2)+'\n')
