"""Bounded row reads and reusable device indices for dosage block aggregation."""
from __future__ import annotations

import numpy as np

READ_BUFFER_BYTES = 256 * 1024**2


class BatchReader:
    """Coalesce strided supplement reads without staging a full dosage matrix."""

    def __init__(self, genotype, batch_size, buffer_bytes=READ_BUFFER_BYTES):
        if batch_size < 1 or buffer_bytes < 0:
            raise ValueError('Invalid dosage batch size or read buffer budget')
        self.genotype = genotype
        self.original = getattr(genotype, 'original', genotype)
        self.supplement = getattr(genotype, 'supplement', None)
        self.batch_size = batch_size
        self.buffer = None
        self.buffer_start = self.buffer_end = 0
        self.read_rows = 0
        if self.supplement is not None and not self.supplement.flags.c_contiguous:
            row_bytes = self.supplement.shape[1] * np.dtype(np.float32).itemsize
            if row_bytes:
                self.read_rows = (buffer_bytes // row_bytes // batch_size) * batch_size
        self.host = np.empty((min(batch_size, genotype.shape[0]), genotype.shape[1]), dtype=np.float32)

    def read(self, start, end):
        if not 0 <= start <= end <= self.genotype.shape[0] or end-start > self.batch_size:
            raise ValueError('Dosage read outside bounded batch')
        host = self.host[:end-start]
        split = self.original.shape[1]
        host[:, :split] = self.original[start:end]
        if self.supplement is not None:
            if self.read_rows >= self.batch_size:
                if start < self.buffer_start or end > self.buffer_end:
                    # Release the old window before allocating its replacement.
                    self.buffer = None
                    self.buffer_start = start
                    self.buffer_end = min(self.genotype.shape[0], start+self.read_rows)
                    self.buffer = np.array(self.supplement[start:self.buffer_end],
                                           dtype=np.float32, order='C', copy=True)
                host[:, split:] = self.buffer[start-self.buffer_start:end-self.buffer_start]
            else:
                host[:, split:] = self.supplement[start:end]
        return host


def score_batches(genotype, maps, tokens, channels, batch_size, device,
                  buffer_bytes=READ_BUFFER_BYTES):
    """Yield the original float32 aggregation with bounded host/device storage."""
    import torch

    device = torch.device(device)
    reader = BatchReader(genotype, batch_size, buffer_bytes)
    constants = [dict(channel=m['channel'],
                      cache_col=torch.as_tensor(m['cache_col'], dtype=torch.long, device=device),
                      token_index=torch.as_tensor(m['token_index'], dtype=torch.long, device=device),
                      coef=torch.as_tensor(m['coef'], device=device),
                      offset=torch.as_tensor(m['offset'], device=device)) for m in maps]
    for start in range(0, genotype.shape[0], batch_size):
        end = min(genotype.shape[0], start+batch_size)
        g = torch.as_tensor(reader.read(start, end), device=device)
        block_out = torch.zeros((end-start, tokens, channels), dtype=torch.float32, device=device)
        for m in constants:
            gp = g.index_select(1, m['cache_col'])
            obs = gp >= 0
            contrib = torch.where(obs, gp*m['coef']+m['offset'], torch.zeros_like(gp))
            target = torch.zeros((end-start, tokens), dtype=torch.float32, device=device)
            target.index_add_(1, m['token_index'], contrib)
            block_out[:, :, m['channel']] += target
        yield start, end, block_out.cpu().numpy()
