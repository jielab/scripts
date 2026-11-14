"""Regression checks for incomplete generated trajectories and outcome support."""
from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import rollout as generation


class DeterministicEvents:
    """Person zero repeatedly has I21; person one dies on the first event."""
    def encode_history(self, batch, **kwargs):
        n = len(batch['eid'])
        return {'age': torch.zeros(n), 'person': torch.arange(n)}

    def predict_next(self, state):
        tokens = torch.where(state['person'] == 0, 4, 3)
        logits = torch.full((len(tokens), 5), -1e9)
        logits.scatter_(1, tokens[:, None], 0)
        return logits, torch.log1p(torch.full((len(tokens),), 0.25))

    def advance_state(self, state, elapsed, **kwargs):
        return {**state, 'age': state['age'] + elapsed}

    def observe_event(self, state, tokens, **kwargs):
        return state


class RolloutTests(unittest.TestCase):
    def test_one_incomplete_particle_makes_person_score_unknown(self):
        first = torch.tensor([[[0.2, float('inf')], [float('inf'), float('inf')]]])
        elapsed = torch.tensor([0.2, 0.1])
        active = torch.tensor([False, True])
        supported = torch.tensor([True, True])
        scores, counts = generation._horizon_scores(first, elapsed, active, 1, supported)
        self.assertTrue(torch.isnan(scores).all())
        self.assertEqual(counts, dict(completed_trajectories=1, incomplete_trajectories=1,
                                     completed_participants=0, affected_participants=1))
        short, counts = generation._horizon_scores(first, elapsed, active, 0.05, supported)
        torch.testing.assert_close(short, torch.zeros(1, 2))
        self.assertEqual(counts['incomplete_trajectories'], 0)

    def test_absorbed_particles_and_unsupported_endpoint_are_distinguished(self):
        first = torch.tensor([[[float('inf'), 0.1, float('inf')],
                               [float('inf'), 0.2, float('inf')]]])
        scores, counts = generation._horizon_scores(first, torch.tensor([0.1, 0.2]),
            torch.tensor([False, False]), 10, torch.tensor([True, True, False]))
        torch.testing.assert_close(scores[:, :2], torch.tensor([[0., 1.]]))
        self.assertTrue(torch.isnan(scores[:, 2]).all())
        self.assertEqual(counts['completed_trajectories'], 2)
        self.assertEqual(counts['affected_participants'], 0)

    def run_toy_rollout(self, root, names):
        cfg = {'paths': {'out': str(root)}, 'seed': 2026,
               'rollout': {'particles': 2, 'horizons': [0.25, 1], 'max_events': 2,
                           'batch_size': 2, 'min_interval_days': 1},
               'readout': {'horizon_years': 5, 'history_chunk_size': 10},
               'train': {'workers': 0}}
        meta = {'vocab': ['<PAD>', '<UNK>', '<START>', 'DEATH', 'I21'],
                'endpoints': {name: [name] for name in names}}

        class Dataset:
            endpoint_names = names

            def __len__(self):
                return 2

        batch = {'eid': ['person0', 'person1']}
        cache = {'eid': np.array(batch['eid'])}
        manifest = {'endpoint_names': names, 'key': 'synthetic'}
        calibration = {'time_residuals': np.array([0], dtype=np.float32),
                       'temperature': 1., 'validation_n': 8}
        with patch.object(generation, 'load_world', return_value=(DeterministicEvents(), meta, torch.device('cpu'))), \
             patch.object(generation, 'calibrate_generation', return_value=calibration), \
             patch.object(generation, 'risk_cache', return_value=(cache, manifest)), \
             patch.object(generation, 'CohortDataset', return_value=Dataset()), \
             patch.object(generation, 'make_loader', return_value=[batch]), \
             patch.object(generation, 'score_matrix') as score, redirect_stdout(io.StringIO()):
            generation.rollout(cfg, 'cpu')
            calls = score.call_args_list
        return json.loads((root / 'results/rollout.json').read_text()), calls

    def test_incomplete_long_horizon_is_withheld_and_stale_files_are_removed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            results = root / 'results'
            results.mkdir()
            owned = ['rollout.1y.' + suffix for suffix in
                     ('metrics.tsv', 'macro.tsv', 'calibration.tsv', 'provenance.json')]
            unrelated = ['flux.metrics.tsv', 'rollout.10y.metrics.tsv', 'rollout.1y.custom.txt']
            for filename in owned + unrelated:
                (results / filename).write_text('old result')
            info, calls = self.run_toy_rollout(root, ['I21', 'DEATH', 'E11'])
            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0].kwargs['horizons'], [0.25])
            self.assertEqual(calls[0].kwargs['suffix'], 'rollout.0.25y')
            self.assertEqual(info['unsupported_endpoints'], ['E11'])
            self.assertEqual(info['supported_endpoint_count'], 2)
            short, long = info['horizon_status']['0.25'], info['horizon_status']['1']
            self.assertEqual(short['status'], 'complete')
            self.assertEqual(short['completed_trajectories'], 4)
            self.assertTrue(short['performance_reported'])
            self.assertEqual(long['status'], 'insufficient_simulation_budget')
            self.assertEqual(long['completed_trajectories'], 2)
            self.assertEqual(long['incomplete_trajectories'], 2)
            self.assertEqual(long['completed_participants'], 1)
            self.assertEqual(long['affected_participants'], 1)
            self.assertFalse(long['performance_reported'])
            short_scores = np.load(results / 'rollout.0.25y.test.npy')
            long_scores = np.load(results / 'rollout.1y.test.npy')
            np.testing.assert_array_equal(short_scores[:, :2], [[1, 0], [0, 1]])
            self.assertTrue(np.isnan(short_scores[:, 2]).all())
            self.assertTrue(np.isnan(long_scores[0]).all())
            np.testing.assert_array_equal(long_scores[1, :2], [0, 1])
            for filename in owned:
                self.assertFalse((results / filename).exists(), filename)
            for filename in unrelated:
                self.assertEqual((results / filename).read_text(), 'old result')

    def test_no_generatable_endpoint_does_not_report_constant_zero_performance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            info, calls = self.run_toy_rollout(root, ['E11'])
            self.assertEqual(calls, [])
            self.assertEqual(info['horizon_status']['0.25']['status'], 'no_supported_endpoints')
            self.assertEqual(info['supported_endpoint_count'], 0)
            self.assertTrue(np.isnan(np.load(root / 'results/rollout.0.25y.test.npy')).all())


if __name__ == '__main__':
    unittest.main()
