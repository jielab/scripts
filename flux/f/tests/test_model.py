"""Meaningful state/time/gradient checks using a small CPU model."""
import copy
import sys
import unittest
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from model import FluxModel, RiskHead


SMALL = dict(deterministic_dim=16, stochastic_dim=8, embedding_dim=12,
             context_hidden_dim=10, age_dim=6, hidden_dim=24, dropout=0.0,
             max_ode_step_years=0.5, kl_weight=0.1)


def batch():
    return {
        "times": torch.tensor([[50., 50.5, 51.2, 52.], [45., 46., 0., 0.]]),
        "valid": torch.tensor([[True, True, True, True], [True, True, False, False]]),
        "tokens": torch.tensor([[2, 0, 3, 4], [5, 6, 0, 0]]),
        "context": torch.tensor([[0., 1.], [1., 0.]]),
        "context_mask": torch.tensor([[True, True], [True, False]]),
        "numeric": {
            "blood": {"values": torch.arange(24, dtype=torch.float32).reshape(2, 4, 3) / 10,
                      "mask": torch.tensor([[[True, True, False]] * 4, [[True, False, True]] * 4])},
            "prot": {"values": torch.arange(16, dtype=torch.float32).reshape(2, 4, 2) / 10,
                     "mask": torch.tensor([[[False, False], [True, True], [False, False], [False, False]],
                                           [[False, False], [True, False], [False, False], [False, False]]])},
        },
        "next_event": torch.tensor([[3, 3, 4, 0], [6, 0, 0, 0]]),
        "next_dt": torch.tensor([[1.2, 0.7, 0.8, 0.], [1., 0., 0., 0.]]),
        "label_mask": torch.tensor([[True, True, True, False], [True, False, False, False]]),
    }


class ModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)

    def setUp(self):
        torch.manual_seed(7)
        self.model = FluxModel(SMALL, 8, 2, {"blood": 3, "prot": 2}).eval()
        self.data = batch()

    def assertStateClose(self, one, two, atol=1e-6):
        for key in ("h", "z", "age", "context"):
            torch.testing.assert_close(one[key], two[key], rtol=1e-6, atol=atol)

    def test_training_objective_has_finite_nonzero_gradients(self):
        self.model.train()
        output = self.model(self.data, sample=True)
        self.assertEqual(int(output["event_count"]), 4)
        self.assertEqual(tuple(output["features"].shape), (2, 4, 24))
        self.assertTrue(torch.isfinite(output["loss"]))
        output["loss"].backward()
        for name in ("token_embedding.weight", "numeric_encoders.0.0.weight",
                     "vector_output.1.weight", "posterior_net.0.weight",
                     "prior_net.0.weight", "jump.weight_ih", "event_head.weight"):
            param = dict(self.model.named_parameters())[name]
            self.assertIsNotNone(param.grad, name)
            self.assertTrue(torch.isfinite(param.grad).all(), name)
            self.assertGreater(float(param.grad.abs().sum()), 0., name)

    def test_padding_does_not_update_state_or_loss(self):
        a = self.model(self.data, sample=False)
        other = copy.deepcopy(self.data)
        other["times"][1, 2:] = float("nan")
        other["tokens"][1, 2:] = 7
        other["numeric"]["blood"]["values"][1, 2:] = 1e9
        other["numeric"]["blood"]["mask"][1, 2:] = True
        b = self.model(other, sample=False)
        self.assertStateClose(a["final_state"], b["final_state"])
        torch.testing.assert_close(a["loss"], b["loss"])

    def test_masked_values_are_not_imputed_as_observations(self):
        a = self.model(self.data, sample=False)
        other = copy.deepcopy(self.data)
        for value in other["numeric"].values():
            value["values"][~value["mask"]] = float("nan")
        other["context"][~other["context_mask"]] = float("nan")
        b = self.model(other, sample=False)
        self.assertStateClose(a["final_state"], b["final_state"])

    def test_zero_interval_is_identity_even_when_sampling(self):
        state = self.model.encode_history(self.data)
        advanced = self.model.advance_state(state, 0., sample=True)
        self.assertStateClose(state, advanced, atol=0.)

    def test_future_observation_cannot_change_earlier_state(self):
        a = self.model(self.data, sample=False)
        other = copy.deepcopy(self.data)
        other["tokens"][0, 3] = 7
        other["numeric"]["blood"]["values"][0, 3] = -99.
        b = self.model(other, sample=False)
        torch.testing.assert_close(a["features"][:, :3], b["features"][:, :3])
        self.assertFalse(torch.allclose(a["features"][0, 3], b["features"][0, 3]))

    def test_chunking_preserves_full_history(self):
        a = self.model(self.data, sample=False)["final_state"]
        b = self.model.encode_history(self.data, chunk_size=1)
        self.assertStateClose(a, b)

    def test_no_observation_ode_rollout_changes_state_with_gradients(self):
        state = self.model.encode_history(self.data)
        advanced = self.model.advance_state(state, torch.tensor([2., 0.5]))
        self.assertFalse(torch.allclose(state["h"], advanced["h"]))
        torch.testing.assert_close(advanced["age"], torch.tensor([54., 46.5]))
        advanced["h"].square().mean().backward()
        self.assertGreater(float(self.model.vector_output[1].weight.grad.abs().sum()), 0.)

    def test_query_age_is_a_forecast_after_prefix(self):
        a = self.model.encode_history(self.data)
        other = dict(self.data, query_age=torch.tensor([55., 50.]))
        b = self.model.encode_history(other)
        c = self.model.advance_state(a, torch.tensor([3., 4.]))
        self.assertStateClose(b, c)
        with self.assertRaises(ValueError):
            self.model.encode_history(dict(self.data, query_age=torch.tensor([51., 50.])))

    def test_same_day_targets_do_not_create_artificial_order(self):
        other = copy.deepcopy(self.data)
        other["next_dt"][0, 1] = 0.
        output = self.model(other, sample=False)
        self.assertEqual(int(output["event_count"]), 3)
        self.assertFalse(output["label_mask"][0, 1])

    def test_no_next_targets_keeps_zero_head_gradients_for_ddp(self):
        other = copy.deepcopy(self.data)
        other["label_mask"][:] = False
        output = self.model(other, sample=False)
        self.assertEqual(int(output["event_count"]), 0)
        output["loss"].backward()
        self.assertIsNotNone(self.model.event_head.weight.grad)
        self.assertIsNotNone(self.model.time_head[0].weight.grad)
        self.assertEqual(float(self.model.event_head.weight.grad.abs().sum()), 0.)

    def test_time_reversal_and_nonfinite_interval_fail(self):
        other = copy.deepcopy(self.data)
        other["times"][0, 2] = 40.
        with self.assertRaises(ValueError):
            self.model(other)
        state = self.model.encode_history(self.data)
        for bad in (-1., float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                self.model.advance_state(state, bad)

    def test_no_update_retains_instantaneous_posterior_and_kl(self):
        model = FluxModel({**SMALL, "no_update": True}, 8, 2, {"blood": 3, "prot": 2}).eval()
        original = model(self.data, sample=False)
        a = original["final_state"]
        other = copy.deepcopy(self.data)
        other["tokens"][other["valid"]] = 7
        other["numeric"]["blood"]["values"] *= 100
        b = model.encode_history(other)
        torch.testing.assert_close(a["h"], b["h"])
        self.assertFalse(torch.allclose(a["z"], b["z"]))
        self.assertEqual(int(original["kl_count"]), int(self.data["valid"].sum()))
        self.assertGreater(float(original["loss_kl"]), 0.)
        original["loss"].backward()
        self.assertIsNotNone(model.posterior_net[0].weight.grad)
        self.assertGreater(float(model.posterior_net[0].weight.grad.abs().sum()), 0.)
        self.assertIsNone(model.jump.weight_ih.grad)
        # Without recurrent revision the observations never enter h. At a later
        # unobserved time, both identical deterministic states imply the same
        # prior; the instantaneous posterior has not persisted through a GRU.
        advanced = model.advance_state(a, 5.)
        alternate = model.advance_state(b, 5.)
        self.assertFalse(torch.allclose(a["h"], advanced["h"]))
        self.assertStateClose(advanced, alternate)

    def test_no_ode_skips_drift_but_keeps_observation_revision(self):
        model = FluxModel({**SMALL, "no_ode": True}, 8, 2, {"blood": 3, "prot": 2}).eval()
        state = model.initial_state(self.data["context"], self.data["times"][:, 0])
        revised = model.observe_event(state, self.data["tokens"][:, 0])
        self.assertFalse(torch.allclose(state["h"], revised["h"]))
        advanced = model.advance_state(revised, 5.)
        torch.testing.assert_close(revised["h"], advanced["h"])

    def test_generated_event_can_update_but_pad_cannot(self):
        state = self.model.encode_history(self.data)
        pad = self.model.observe_event(state, torch.tensor([0, 0]))
        self.assertStateClose(state, pad)
        actual = self.model.observe_event(state, torch.tensor([2, 3]))
        self.assertFalse(torch.allclose(state["h"], actual["h"]))
        logits, time = self.model.predict_next(actual)
        self.assertEqual(tuple(logits.shape), (2, 8))
        self.assertTrue((time >= 0).all())

    def test_empty_history_needs_explicit_age(self):
        empty = {"times": torch.empty(2, 0), "valid": torch.empty(2, 0, dtype=torch.bool),
                 "tokens": torch.empty(2, 0, dtype=torch.long), "context": torch.zeros(2, 2),
                 "numeric": {}}
        with self.assertRaises(ValueError):
            self.model.encode_history(empty)
        empty["start_age"] = torch.tensor([50., 55.])
        state = self.model.encode_history(empty)
        self.assertTrue(torch.isfinite(self.model.state_features(state)).all())

    def test_risk_head_consumes_state_features(self):
        states = self.model.encode_history(self.data)
        for hidden_dim in (12, [12, 6]):
            head = RiskHead(self.model.state_dim, 4, hidden_dim=hidden_dim, dropout=0.)
            logits = head(self.model.state_features(states).detach())
            self.assertEqual(tuple(logits.shape), (2, 4))
            logits.square().mean().backward()
            self.assertGreater(float(head.net[0].weight.grad.abs().sum()), 0.)

    def test_bad_config_and_modality_fail_explicitly(self):
        with self.assertRaises(ValueError):
            FluxModel({"unrecognized": 1}, 8, 2, {})
        other = copy.deepcopy(self.data)
        other["numeric"]["future_unknown"] = other["numeric"]["blood"]
        with self.assertRaises(ValueError):
            self.model(other)


if __name__ == "__main__":
    unittest.main()
