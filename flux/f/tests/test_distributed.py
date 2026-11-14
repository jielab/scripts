"""Actual CPU/Gloo validation reductions; no CUDA or DDP training claim."""
from __future__ import annotations

from datetime import timedelta
import json
from pathlib import Path
import sys
import tempfile
import unittest

import torch
import torch.distributed as dist
import torch.multiprocessing as mp
from torch.utils.data import DataLoader

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from train import validation_loss


class KnownValidationSums(torch.nn.Module):
    """Unequal event/KL counts expose incorrect averaging of rank means."""

    model_cfg = {'time_weight': 2., 'kl_weight': 0.125}

    def forward(self, people, **kwargs):
        values = people.to(torch.float64)
        event_count = (people.to(torch.long).remainder(3) + 1).to(torch.float64)
        kl_count = values + 1.
        return dict(event_sum=((values + .5) * event_count).sum(),
                    time_sum=((values / 4. + .25) * event_count).sum(),
                    kl_balanced_sum=((values + .25) * kl_count).sum(),
                    kl_sum=torch.tensor(1000.),
                    event_count=event_count.sum(), kl_count=kl_count.sum())


def _validation_worker(rank, rendezvous, folder):
    """Module-level worker is importable under multiprocessing spawn."""
    torch.set_num_threads(1)
    dist.init_process_group('gloo', init_method=rendezvous, rank=rank,
                            world_size=2, timeout=timedelta(seconds=30))
    try:
        people = torch.arange(5, dtype=torch.float64)
        model = KnownValidationSums().train()
        results = {}
        # Five people make unequal shards and unequal batch counts. With one
        # person, rank 1 runs no forward calls but must still join all_reduce.
        for limit in (5, 1):
            sampler = range(rank, limit, 2)
            loader = DataLoader(people, batch_size=2, sampler=sampler)
            results[str(limit)] = dict(
                loss=validation_loss(model, loader, torch.device('cpu'), 0,
                                     distributed=True),
                indices=list(sampler), mode_restored=model.training)
        model.eval()
        try:
            validation_loss(model, [], torch.device('cpu'), 0, distributed=True)
        except ValueError as error:
            results['empty_error'] = str(error)
        else:
            results['empty_error'] = None
        results['eval_mode_restored'] = not model.training
        (Path(folder) / f'rank_{rank}.json').write_text(json.dumps(results))
    finally:
        dist.destroy_process_group()


@unittest.skipUnless(dist.is_available() and dist.is_gloo_available(),
                     'This PyTorch build has no Gloo distributed backend')
class DistributedValidationTests(unittest.TestCase):
    def test_nonpadding_shards_and_empty_rank_match_serial(self):
        with tempfile.TemporaryDirectory(prefix='flux_gloo_') as temporary:
            root = Path(temporary)
            rendezvous = (root / 'rendezvous').resolve().as_uri()
            try:
                mp.spawn(_validation_worker, args=(rendezvous, str(root)),
                         nprocs=2, join=True)
            except mp.ProcessRaisedException as error:
                message = str(error)
                # Some hosted runtimes provide Gloo but prohibit its local
                # interface/socket setup. Report that unavailable capability;
                # do not replace the real collective with a mock or a pass.
                if ('gloo/transport/tcp/device.cc' in message and
                        'Operation not permitted' in message):
                    self.skipTest('Execution environment forbids local Gloo networking')
                raise
            ranks = [json.loads((root / f'rank_{rank}.json').read_text())
                     for rank in range(2)]
            people = torch.arange(5, dtype=torch.float64)
            for limit in (5, 1):
                serial = validation_loss(KnownValidationSums(),
                    DataLoader(people, batch_size=2, sampler=range(limit)),
                    torch.device('cpu'), 0)
                selected = []
                for result in ranks:
                    actual = result[str(limit)]
                    self.assertAlmostEqual(actual['loss'], serial, places=13)
                    self.assertTrue(actual['mode_restored'])
                    selected.extend(actual['indices'])
                self.assertEqual(sorted(selected), list(range(limit)))
                self.assertEqual(len(selected), len(set(selected)))
            self.assertEqual(ranks[1]['1']['indices'], [])
            for result in ranks:
                self.assertIn('No next-event targets', result['empty_error'])
                self.assertTrue(result['eval_mode_restored'])


if __name__ == '__main__':
    unittest.main()
