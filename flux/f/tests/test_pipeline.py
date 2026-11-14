"""Checkpoint/RNG continuation and frozen-readout lifecycle on synthetic data.

This exercises execution and leakage boundaries, not predictive performance.
No real UKB files are needed and every output stays in a temporary directory.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import file_digest, load_config, to_device, torch_load
from data import CohortDataset, prepare
from demo import make_demo
from train import (extract_states, fit_readout, load_world, make_loader,
                   risk_cache, train_world, validation_loss)


class PipelineTests(unittest.TestCase):
    def test_checkpoint_resume_and_frozen_complete_prefix_readout(self):
        with tempfile.TemporaryDirectory(prefix="flux_pipeline_") as temporary:
            root = Path(temporary)
            cfg = load_config()
            cfg['paths']['out'] = str(root / 'continued')
            cfg['resources']['cores'] = 1
            cfg['data']['max_nodes'] = 4
            cfg['model'].update(deterministic_dim=8, stochastic_dim=4,
                embedding_dim=8, context_hidden_dim=6, age_dim=4,
                hidden_dim=12, dropout=0.1, max_ode_step_years=3.)
            cfg['train'].update(batch_size=4, accumulation_steps=2,
                epochs=2, max_steps=1, validate_every=1,
                validation_batches=1, workers=0, patience=20)
            cfg['readout'].update(hidden_dims=[12, 6], batch_size=16,
                state_batch_size=8, history_chunk_size=3, epochs=2,
                min_cases=2, patience=2)
            cfg = make_demo(cfg, n_people=120)
            prepare(cfg)

            # Both runs use precisely the same prepared schema and source files.
            uninterrupted = copy.deepcopy(cfg)
            uninterrupted['paths']['out'] = str(root / 'uninterrupted')
            (root / 'uninterrupted' / 'data').mkdir(parents=True)
            for name in ('meta.json', 'cohort.sqlite'):
                shutil.copy2(root / 'continued' / 'data' / name,
                             root / 'uninterrupted' / 'data' / name)

            train_world(cfg, 'cpu')
            first = torch_load(root / 'continued' / 'model' / 'flux.last.pt')
            self.assertEqual(first['step'], 1)
            self.assertEqual(first['next_batch'], 2)
            cfg['train']['max_steps'] = 3
            train_world(cfg, 'cpu', resume=True)
            uninterrupted['train']['max_steps'] = 3
            train_world(uninterrupted, 'cpu')
            resumed = torch_load(root / 'continued' / 'model' / 'flux.last.pt')
            full = torch_load(root / 'uninterrupted' / 'model' / 'flux.last.pt')
            self.assertEqual(resumed['step'], 3)
            self.assertEqual(resumed['next_batch'], 6)
            self.assertEqual(resumed['best'], full['best'])
            for name in full['model']:
                torch.testing.assert_close(resumed['model'][name], full['model'][name],
                                           rtol=0, atol=0, msg=name)

            changed = copy.deepcopy(cfg)
            changed['train']['batch_size'] = 8
            with self.assertRaisesRegex(ValueError, 'Cannot exactly resume'):
                train_world(changed, 'cpu', resume=True)

            world_file = root / 'continued' / 'model' / 'flux.pt'
            frozen_digest = file_digest(world_file)
            extract_states(cfg, 'cpu')
            model, meta, device = load_world(cfg, 'cpu')
            dataset = CohortDataset(cfg, 'train', mode='risk', horizon_years=5.)
            self.assertGreater(len(dataset[0]['times']), cfg['data']['max_nodes'])
            loader = make_loader(dataset, cfg['readout']['state_batch_size'])
            batch = next(iter(loader))
            with torch.no_grad():
                state = model.encode_history(to_device(batch, device), chunk_size=3)
                expected = model.state_features(state).cpu().numpy()
            arrays, manifest = risk_cache(cfg, 'train')
            self.assertTrue(manifest['complete_history'])
            np.testing.assert_allclose(arrays['features'][:len(expected)], expected,
                                       rtol=1e-6, atol=1e-6)
            np.testing.assert_array_equal(arrays['eid'][:len(expected)], batch['eid'])
            self.assertTrue((batch['query_age'] >= state['age'].cpu() - 1e-5).all())

            fit_readout(cfg, 'cpu')
            self.assertEqual(file_digest(world_file), frozen_digest)
            readout_file = root / 'continued' / 'model' / 'readout.pt'
            fitted = torch_load(readout_file)
            self.assertEqual(fitted['estimand'], 'recorded_incidence')
            self.assertTrue(fitted['eligible'].any())
            # Normalization is fitted only on the training-state cache.
            np.testing.assert_allclose(fitted['state_mean'],
                np.asarray(arrays['features'], dtype=np.float64).mean(0), rtol=1e-6, atol=1e-6)
            predictions = np.load(root / 'continued' / 'results' / 'flux.test.npy')
            test, test_manifest = risk_cache(cfg, 'test')
            self.assertEqual(predictions.shape, test['event_time'].shape)
            self.assertTrue(np.isfinite(predictions[:, fitted['eligible']]).all())
            self.assertTrue(np.isnan(predictions[:, ~fitted['eligible']]).all())
            self.assertTrue(((predictions[:, fitted['eligible']] >= 0) &
                             (predictions[:, fitted['eligible']] <= 1)).all())
            # Reusing a fitted readout reproduces its held-out predictions.
            checkpoint_digest = file_digest(readout_file)
            fit_readout(cfg, 'cpu')
            self.assertEqual(file_digest(readout_file), checkpoint_digest)
            np.testing.assert_allclose(
                np.load(root / 'continued' / 'results' / 'flux.test.npy'), predictions,
                rtol=0, atol=0, equal_nan=True)
            metadata = json.loads((root / 'continued' / 'results' /
                                   'flux.prediction.json').read_text())
            self.assertEqual(metadata['endpoint_names'], test_manifest['endpoint_names'])

    def test_validation_uses_balanced_kl_and_restores_mode(self):
        class KnownSums(torch.nn.Module):
            model_cfg = {'time_weight': 2., 'kl_weight': 0.5}

            def forward(self, batch, **kwargs):
                return batch

        model = KnownSums().eval()
        batches = [dict(event_sum=6., time_sum=2., kl_sum=1., kl_balanced_sum=8.,
                        event_count=2., kl_count=4.),
                   dict(event_sum=3., time_sum=1., kl_sum=1., kl_balanced_sum=4.,
                        event_count=1., kl_count=2.)]
        self.assertAlmostEqual(validation_loss(model, batches, torch.device('cpu'), 0), 6.)
        self.assertFalse(model.training)


if __name__ == '__main__':
    unittest.main()
