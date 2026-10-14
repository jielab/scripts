#!/usr/bin/env python

"""
Markov Chain Monte Carlo (MCMC) sampler for cross-ethnic polygenic prediction with continuous shrinkage (CS) priors - PRS-CSx.

"""


import os
import json
import time
import fcntl
from functools import wraps
from gpu_linalg import BlockLinalg, gpu_min_block, linalg_identity
import h5py
import numpy as np
from numpy import random
import gigrnd


def _gpu_serialized(function):
    @wraps(function)
    def run(*args, **kwargs):
        backend = kwargs.get('backend') or (args[20] if len(args)>20 else None) or os.getenv('GRID_CSX_BACKEND', 'cpu')
        device = kwargs.get('gpu_device')
        if device is None: device = args[21] if len(args)>21 else os.getenv('GRID_CSX_GPU', '0')
        if backend == 'cpu': return function(*args, **kwargs)
        sizes = kwargs.get('blk_size', args[9] if len(args)>9 else {})
        # Pure CPU execution needs neither CUDA memory nor its interprocess lock.
        # The launcher still checks requested CUDA availability. This is a
        # deterministic size policy, never a fallback after a numerical failure.
        if gpu_min_block() and max((n for group in sizes.values() for n in group), default=0) < gpu_min_block():
            kwargs['_cpu_only'] = True
            return function(*args, **kwargs)
        with open('/tmp/grid-csx-gpu-%d.lock' % int(device), 'a') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                print('... Waiting for PRS-CSx GPU device %d lock ...' % int(device), flush=True)
                fcntl.flock(lock, fcntl.LOCK_EX)
                print('... PRS-CSx GPU lock acquired ...', flush=True)
            return function(*args, **kwargs)
    return run


@_gpu_serialized
def mcmc(a, b, phi, snp_dict, beta_mrg, frq_dict, idx_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, pop, chrom, out_dir, out_name, meta, write_pst, seed, backend=None, gpu_device=None, _cpu_only=False):
    started = time.perf_counter()
    backend = backend or os.getenv('GRID_CSX_BACKEND', 'cpu')
    gpu_device = int(gpu_device if gpu_device is not None else os.getenv('GRID_CSX_GPU', '0'))
    engine = BlockLinalg('cpu' if _cpu_only else backend, gpu_device)
    identity = linalg_identity(backend)
    if _cpu_only:
        print('... requested=%s; all LD blocks < %d; executing on CPU ...' % (backend, gpu_min_block()), flush=True)
    sizes = [size for group in blk_size.values() for size in group if size > 0]
    print('... chr%d LD blocks=%d; max size=%d; GPU blocks per iteration=%d ...' %
          (chrom, len(sizes), max(sizes, default=0),
           sum(size >= gpu_min_block() for size in sizes) if engine.is_cuda else 0), flush=True)
    print('... PRS-CSx linear algebra:', engine.info(), '...', flush=True)
    write_seconds = 0.0
    print('... MCMC ...', flush=True)

    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    # Samples are retained at absolute iterations divisible by thin.
    n_pst = n_iter // thin - n_burnin // thin
    if n_pst < 2:
        raise ValueError('Need at least two retained posterior samples')
    n_pop = len(pop)
    p_tot = len(snp_dict['SNP'])
    # Keep the same ordered indices, without rebuilding NumPy index arrays
    # from Python lists on every population update.
    idx_dict = {pp: np.asarray(idx_dict[pp], dtype=int) for pp in range(n_pop)}

    p = {}
    n_blk = {}
    het = {}
    for pp in range(n_pop):
        p[pp] = len(beta_mrg[pp])
        n_blk[pp] = len(ld_blk[pp])
        het[pp] = np.sqrt(2.0*frq_dict[pp]*(1.0-frq_dict[pp]))

    n_grp = np.zeros((p_tot,1))
    for pp in range(n_pop):
        n_grp[np.asarray(idx_dict[pp],dtype=int)] += 1

    # initialization
    beta = {}
    sigma = {}
    for pp in range(n_pop):
        beta[pp] = np.zeros((p[pp],1))
        sigma[pp] = 1.0

    psi = np.ones((p_tot,1))

    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    # Save each population at the SAME retained iteration. Disk-backed datasets
    # avoid keeping four SNP-by-draw matrices in RAM. META stays a mean-weight file.
    draw_file = None
    if write_pst == 'TRUE':
        draw_path = os.path.join(out_dir, 'joint_posterior.h5')
        draw_file = h5py.File(draw_path + '.tmp', 'w')
        draw_file.attrs.update(schema='grid_csx_draws_v1', chromosome=chrom,
                               backend=engine.backend, requested_backend=backend, linalg=identity,
                               gpu_min_block=gpu_min_block(), seed=-1 if seed is None else seed, complete=False,
                               effect_scale='standardized_phenotype_per_allele')
        draw_file.create_dataset('iteration', data=np.arange(thin, n_iter+1, thin)[np.arange(thin,n_iter+1,thin)>n_burnin])
        for pp, population in enumerate(pop):
            group = draw_file.create_group(population)
            ii = idx_dict[pp]
            for name in ('SNP','A1','A2'):
                group.create_dataset(name, data=np.asarray([snp_dict[name][j] for j in ii], dtype=h5py.string_dtype()))
            group.create_dataset('reference_eaf', data=frq_dict[pp].ravel())
            group.create_dataset('beta', shape=(p[pp],n_pst), dtype='f8',
                                 chunks=(min(p[pp],4096),1), compression='gzip',compression_opts=1)

    # space allocation
    beta_est = {}
    beta_sq_est = {}
    sigma_est = {}
    for pp in range(n_pop):
        beta_est[pp] = np.zeros((p[pp],1))
        beta_sq_est[pp] = np.zeros((p[pp],1))
        sigma_est[pp] = 0.0

    psi_est = np.zeros((p_tot,1))
    phi_est = 0.0

    # MCMC
    mcmc_started = time.perf_counter()
    last_progress = mcmc_started
    progress_path = os.path.join(out_dir, 'mcmc.progress.json')
    def progress(iteration):
        elapsed = time.perf_counter() - mcmc_started
        state = dict(chromosome=chrom, iteration=iteration, iterations=n_iter,
                     retained_draws=qq, elapsed_seconds=elapsed, requested_backend=backend,
                     backend=engine.backend, cpu_block_calls=engine.cpu_calls,
                     cuda_block_calls=engine.cuda_calls, updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
        with open(progress_path+'.tmp', 'w') as handle: json.dump(state, handle, indent=2)
        os.replace(progress_path+'.tmp', progress_path)
        print('... chr%d iteration %d/%d; elapsed=%.1fs; CPU/GPU blocks=%d/%d ...' %
              (chrom, iteration, n_iter, elapsed, engine.cpu_calls, engine.cuda_calls), flush=True)
    qq = 0
    progress(0)
    for itr in range(1,n_iter+1):

        for pp in range(n_pop):
            mm = 0; quad = 0.0
            psi_pp = psi[idx_dict[pp]]
            for kk in range(n_blk[pp]):
                if blk_size[pp][kk] == 0:
                    continue
                else:
                    idx_blk = slice(mm,mm+blk_size[pp][kk])
                    noise = random.randn(blk_size[pp][kk], 1)
                    updated, q = engine.beta_block(ld_blk[pp][kk], psi_pp[idx_blk],
                        beta_mrg[pp][idx_blk], noise, sigma[pp]/n[pp])
                    beta[pp][idx_blk] = updated
                    quad += q
                    mm += blk_size[pp][kk]

            err = max(n[pp]/2.0*(1.0-2.0*sum(beta[pp]*beta_mrg[pp])+quad), n[pp]/2.0*sum(beta[pp]**2/psi_pp))
            sigma[pp] = 1.0/random.gamma((n[pp]+p[pp])/2.0, 1.0/err)

        delta = random.gamma(a+b, 1.0/(psi+phi))

        xx = np.zeros((p_tot,1))
        for pp in range(n_pop):
            xx[idx_dict[pp]] += n[pp]*beta[pp]**2/sigma[pp]

        for jj in range(p_tot):
            for attempt in range(100):
                try:
                    psi[jj] = gigrnd.gigrnd(a-0.5*n_grp[jj], 2.0*delta[jj], xx[jj])
                except (ArithmeticError, ValueError):
                    continue
                else:
                    break
            else:
                raise FloatingPointError('GIG sampling failed 100 times: chr%d iteration %d SNP %s' %
                                         (chrom, itr, snp_dict['SNP'][jj]))
        psi[psi>1] = 1.0

        if phi_updt == True:
            w = random.gamma(1.0, 1.0/(phi+1.0))
            phi = random.gamma(p_tot*b+0.5, 1.0/(sum(delta)+w))

        # posterior
        if (itr > n_burnin) and (itr % thin == 0):
            for pp in range(n_pop):
                beta_est[pp] = beta_est[pp] + beta[pp]/n_pst
                beta_sq_est[pp] = beta_sq_est[pp] + beta[pp]**2/n_pst
                sigma_est[pp] = sigma_est[pp] + sigma[pp]/n_pst

            psi_est = psi_est + psi/n_pst
            phi_est = phi_est + phi/n_pst

            if draw_file is not None:
                write_started = time.perf_counter()
                for pp, population in enumerate(pop):
                    draw_file[population]['beta'][:,qq] = (beta[pp]/het[pp]).ravel()
                write_seconds += time.perf_counter()-write_started
            qq += 1

        if itr % 100 == 0 or itr == n_iter or time.perf_counter()-last_progress >= 60:
            progress(itr)
            last_progress = time.perf_counter()

    mcmc_seconds = time.perf_counter()-mcmc_started-write_seconds
    output_started = time.perf_counter()
    # convert standardized beta to per-allele beta
    for pp in range(n_pop):
        beta_est[pp] /= het[pp]
        beta_sq_est[pp] /= het[pp]**2

    if draw_file is not None:
        if qq != n_pst: raise RuntimeError('Retained draw count mismatch')
        draw_file.attrs['complete'] = True
        draw_file.close()
        os.replace(draw_path + '.tmp', draw_path)

    # meta
    if meta == 'TRUE':
        vv = np.zeros((p_tot,1))
        zz = np.zeros((p_tot,1))
        for pp in range(n_pop):
            variance = np.maximum(beta_sq_est[pp]-beta_est[pp]**2, np.finfo(float).tiny)
            vv[idx_dict[pp]] += 1.0/variance
            zz[idx_dict[pp]] += beta_est[pp]/variance
        mu = zz/vv

    # write posterior effect sizes
    for pp in range(n_pop):
        if phi_updt == True:
            eff_file = os.path.join(out_dir, pop[pp], '%s.chr%d.pst_eff.txt' % (pop[pp], chrom))
        else:
            eff_file = os.path.join(out_dir, pop[pp], '%s.chr%d.pst_eff.txt' % (pop[pp], chrom))

        snp_pp = [snp_dict['SNP'][ii] for ii in idx_dict[pp]]
        bp_pp = [snp_dict['BP'][ii] for ii in idx_dict[pp]]
        a1_pp = [snp_dict['A1'][ii] for ii in idx_dict[pp]]
        a2_pp = [snp_dict['A2'][ii] for ii in idx_dict[pp]]

        with open(eff_file, 'w') as ff:
            for snp, bp, a1, a2, beta in zip(snp_pp, bp_pp, a1_pp, a2_pp, beta_est[pp]):
                ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta[0]))

    if meta == 'TRUE':
        if phi_updt == True:
            eff_file = out_dir + '/' + '%s_META_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (out_name, a, b, chrom)
        else:
            eff_file = out_dir + '/' + '%s_META_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (out_name, a, b, phi, chrom)

        with open(eff_file, 'w') as ff:
            for snp, bp, a1, a2, beta in zip(snp_dict['SNP'], snp_dict['BP'], snp_dict['A1'], snp_dict['A2'], mu):
                ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta[0]))


    # print estimated phi
    if phi_updt == True:
        print('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )

    engine.synchronize()
    report = dict(module='grid-csx', **engine.info(), requested_backend=backend, linalg=identity,
        execution_policy='small-block-cpu-v1' if backend != 'cpu' and gpu_min_block() else 'uniform-v1',
        seed=seed, retained_draws=n_pst,
        n=p_tot, rank_or_blocks=sum(sum(size>0 for size in blk_size[pp]) for pp in range(n_pop)), mcmc_seconds=mcmc_seconds,
        posterior_write_seconds=write_seconds+time.perf_counter()-output_started,
        wall_seconds=time.perf_counter()-started)
    path = os.path.join(out_dir, 'mcmc.timing.json')
    with open(path+'.tmp', 'w') as handle: json.dump(report, handle, indent=2)
    os.replace(path+'.tmp', path)
    print('... Done ...', report, flush=True)
