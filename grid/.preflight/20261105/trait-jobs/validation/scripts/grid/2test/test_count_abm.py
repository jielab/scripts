"""Scientific invariants of the transparent count workflow (synthetic only)."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from count_abm.core import count_predict, summarize_neighbors, freeze_policy, apply_policy
from count_abm.inputs import assign_roles, write_input_manifest, load_export, load_production, read_table
from count_abm.workflow import main, ROOT, verify_stage


class CountTests(unittest.TestCase):
    def test_all_neighbors_contribute_not_a_128_cap(self):
        donor = np.arange(200, dtype=float).reshape(-1, 1)/200
        traces = []
        result = count_predict(donor, np.array([[.5], [10.]]),
            donor_ids=[f"d{i}" for i in range(200)], query_ids=["q", "far"],
            donor_groups=[f"g{i}" for i in range(200)], query_groups=["qg", "fg"],
            donor_y=np.arange(200.), task="continuous", radius=1., min_count=20,
            min_group_ess=15, device="cpu", query_batch=1, donor_batch=37,
            trace_sink=lambda *args: traces.append(args))
        near, far = result.table.iloc[0], result.table.iloc[1]
        self.assertEqual(near.n_twins, 200)
        self.assertEqual(near.p_abm, 99.5)
        self.assertEqual(len(traces[0][3]), 200)
        self.assertEqual(near.outcome_sum, 19900)
        self.assertFalse(far.supported)
        self.assertTrue(np.isnan(far.p_abm))
        self.assertEqual(far.p_population, 99.5)
        self.assertEqual(far.support_reason, "no_twins_in_radius")

    def test_self_family_exclusion_binary_case_preservation_and_exact_count(self):
        x = np.array([[0.], [0.], [.1], [.2], [.9]])
        y = np.array([1., 1., 1., 0., 1.])
        result = count_predict(x, np.array([[0.]]), donor_ids=["q", "relative", "case", "control", "outside"],
            query_ids=["q"], donor_groups=["family", "family", "g1", "g2", "g3"],
            query_groups=["family"], donor_y=y, task="binary", radius=.2,
            min_count=2, min_group_ess=2, device="cpu").table.iloc[0]
        self.assertEqual(result.n_eligible, 3)
        self.assertEqual(result.n_twins, 2)
        self.assertEqual(result.n_cases, 1)
        self.assertEqual(result.n_controls, 1)
        self.assertEqual(result.p_count, .5)
        self.assertTrue(result.supported)

    def test_declared_group_ess_and_no_independence_claim(self):
        stats = summarize_neighbors([0., 1., 1., 0.], ["a", "a", "a", "b"], "binary")
        self.assertAlmostEqual(stats["ess_group"], 16/10)
        result = count_predict(np.zeros((4, 1)), np.zeros((1, 1)),
            donor_ids=list("abcd"), query_ids=["q"], donor_groups=["g", "g", "g", "h"],
            query_groups=["qg"], donor_y=[0, 1, 1, 0], task="binary", radius=1,
            min_count=2, min_group_ess=2, device="cpu").table.iloc[0]
        self.assertEqual(result.support_reason, "low_group_ess")
        self.assertFalse(result.groups_verified)

    def test_membership_and_results_are_donor_order_batch_invariant(self):
        rng = np.random.default_rng(4)
        x, q = rng.normal(size=(60, 4)), rng.normal(size=(8, 4))
        ids = np.asarray([f"d{i:03d}" for i in range(60)])
        y = rng.normal(size=60)
        options = dict(query_ids=[f"q{i}" for i in range(8)], query_groups=[f"qg{i}" for i in range(8)],
                       task="continuous", radius=2.4, min_count=2, min_group_ess=2, device="cpu")
        a = count_predict(x, q, donor_ids=ids, donor_groups=ids, donor_y=y, donor_batch=9, query_batch=3, **options).table
        order = rng.permutation(len(x))
        b = count_predict(x[order], q, donor_ids=ids[order], donor_groups=ids[order], donor_y=y[order], donor_batch=17, query_batch=5, **options).table
        pd.testing.assert_frame_equal(a, b, check_exact=False, rtol=1e-12, atol=1e-12)

    def test_policy_refuses_query_y_and_separates_support_from_heterogeneity(self):
        q = pd.DataFrame(dict(id=list("abcd"), supported=[False, True, True, True], ess_group=[0, 20, 30, 40],
                              neighbor_variance=[np.nan, .1, .2, .9], p_population=[.1, .2, .3, .4]))
        policy = freeze_policy(q, [.5, 1.], heterogeneity_quantile=.5)
        pred = apply_policy(q, policy)
        self.assertEqual(list(pred.evidence_state), ["low_support", "supported_lower_dispersion", "supported_lower_dispersion", "supported_heterogeneous"])
        with self.assertRaisesRegex(ValueError, "query outcomes"):
            freeze_policy(q.assign(y=0))
        with self.assertRaisesRegex(ValueError, "before attaching"):
            apply_policy(q.assign(y=0), policy)

    def test_roles_do_not_depend_on_outcomes_or_row_order(self):
        d = pd.DataFrame(dict(id=[str(i) for i in range(300)], group_id=[str(i//2) for i in range(300)],
                              split=["train"]*240+["test"]*60, y=np.arange(300.)))
        a = assign_roles(d).set_index("id")
        b = assign_roles(d.assign(y=np.nan).sample(frac=1, random_state=7)).set_index("id").loc[a.index]
        pd.testing.assert_series_equal(a.role, b.role)
        d.loc[240, "group_id"] = d.loc[0, "group_id"]
        with self.assertRaisesRegex(ValueError, "crosses"):
            assign_roles(d)

    def test_auto_reports_actual_backend_and_explicit_cuda_never_silently_falls_back(self):
        kwargs = dict(donor_ids=["a", "b"], query_ids=["q"], donor_groups=["ga", "gb"],
                      query_groups=["gq"], donor_y=[0, 1], task="binary", radius=1., min_count=1, min_group_ess=1)
        result = count_predict(np.zeros((2, 1)), np.zeros((1, 1)), device="auto", **kwargs)
        self.assertIn(result.device_used, {"cpu", "cuda"})
        try:
            import torch
            cuda = torch.cuda.is_available()
        except ImportError:
            cuda = False
        if not cuda:
            self.assertEqual(result.device_used, "cpu")
            with self.assertRaisesRegex(RuntimeError, "cuda"):
                count_predict(np.zeros((2, 1)), np.zeros((1, 1)), device="cuda", **kwargs)


def exported_fixture(folder, *, test_shift=0., binary=False, maps=4):
    folder.mkdir(parents=True)
    rng = np.random.default_rng(42)
    x = rng.normal(size=(200, 8))
    y = (x[:, 0] > 0).astype(float) if binary else x[:, 0]+rng.normal(scale=.3, size=len(x))
    if not binary:
        y[160:] += test_shift
    else:
        if test_shift:
            y[160:] = 1-y[160:]
    y[-1] = np.nan
    d = pd.DataFrame(dict(id=[f"{i:05d}" for i in range(len(x))], group_id=[f"g{i//2}" for i in range(len(x))],
                          split=["train"]*120+["validation"]*40+["test"]*40,
                          ancestry=["EUR", "AFR"]*(len(x)//2), y=y))
    d.to_csv(folder/"cohort.tsv", sep="\t", index=False)
    np.save(folder/"blocks.npy", x)
    features = pd.DataFrame(dict(feature_id=[f"block{i}" for i in range(x.shape[1])],
                                 map_id=[f"map{i%maps}" for i in range(x.shape[1])]))
    features.to_csv(folder/"features.tsv", sep="\t", index=False)
    trait = "t2dm" if binary else "height"
    write_input_manifest(folder/"cohort.tsv", folder/"blocks.npy", folder/"features.tsv", folder/"input.json",
                         trait=trait, task="binary" if binary else "continuous")
    return folder/"input.json", trait


class WorkflowTests(unittest.TestCase):
    def invoke(self, root, manifest, trait, phase="all"):
        with contextlib.redirect_stdout(io.StringIO()):
            main([phase, "--traits", trait, "--analysis-root", str(root), "--input-manifest", str(manifest),
                  "--device", "cpu", "--radius-k-ref", "8", "--min-count", "5", "--min-group-ess", "3",
                  "--bootstrap", "8", "--threads", "1"])

    def test_full_pipeline_frozen_memberships_recompute_prediction_and_retains_missing_y(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            manifest, trait = exported_fixture(base/"input", binary=True, maps=1)
            # Poison old output: the current final must never consume it.
            (base/"analysis/final").mkdir(parents=True)
            (base/"analysis/final/paper_main_performance.tsv").write_text("broken old output must not be parsed")
            self.invoke(base/"analysis", manifest, trait)
            target = base/"analysis/count_abm"/trait
            pred = read_table(target/"_LOCAL_ONLY/test.predictions.tsv.gz")
            donor = read_table(target/"_LOCAL_ONLY/donor.tsv.gz")
            self.assertEqual(len(pred), 40)
            self.assertNotIn("y", pred)
            self.assertEqual(pred.id.iloc[0], "00160")
            for path in (target/"_LOCAL_ONLY/test.memberships").glob("*.npz"):
                with np.load(path) as membership:
                    for j, query_id in enumerate(membership["query_ids"]):
                        indices = membership["donor_indices"][membership["indptr"][j]:membership["indptr"][j+1]]
                        row = pred.loc[pred.id.eq(query_id)].iloc[0]
                        self.assertEqual(len(indices), row.n_twins)
                        if len(indices):
                            self.assertAlmostEqual(float(donor.iloc[indices].y.mean()), row.p_count)
            perf = read_table(target/"performance.tsv")
            full = perf.loc[perf.role.eq("test") & perf.population.eq("all_with_explicit_fallback")].iloc[0]
            self.assertEqual(full.n_evaluable, 39)
            self.assertEqual(full.n_roster, 40)
            self.assertTrue((base/"analysis/count_abm/final/index.html").is_file())
            refs = pd.read_csv(target/"_LOCAL_ONLY/test.reference_examples.tsv.gz", sep="\t")
            blocks = pd.read_csv(target/"_LOCAL_ONLY/test.block_examples.tsv.gz", sep="\t")
            self.assertTrue(len(refs))
            self.assertTrue({"case_reference", "control_reference"} <= set(refs.reference_type))
            self.assertTrue(blocks.fraction_of_pair_distance.between(0, 1).all())
            self.invoke(base/"analysis", manifest, trait)

    def test_test_outcome_perturbation_leaves_geometry_radius_policy_predictions_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            outputs = []
            for i, shift in enumerate([0., 1000.]):
                manifest, trait = exported_fixture(base/f"input{i}", test_shift=shift)
                analysis = base/f"analysis{i}"
                self.invoke(analysis, manifest, trait, "prepare")
                self.invoke(analysis, manifest, trait, "twins")
                target = analysis/"count_abm"/trait
                outputs.append((json.loads((target/"geometry.json").read_text()),
                                json.loads((target/"policy.json").read_text()),
                                json.loads((target/"model.json").read_text())["radius"],
                                read_table(target/"_LOCAL_ONLY/test.predictions.tsv.gz")))
            self.assertEqual(outputs[0][:3], outputs[1][:3])
            pd.testing.assert_frame_equal(outputs[0][3], outputs[1][3])

    def test_input_binding_and_frozen_stage_tampering_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            manifest, trait = exported_fixture(base/"input")
            self.invoke(base/"analysis", manifest, trait, "prepare")
            target = base/"analysis/count_abm"/trait
            with (target/"_LOCAL_ONLY/test.tsv.gz").open("ab") as stream:
                stream.write(b"tampered")
            with self.assertRaisesRegex(ValueError, "artifact missing or changed"):
                verify_stage(target, "prepare")
            x = np.load(base/"input/blocks.npy")
            np.save(base/"input/blocks.npy", x[::-1])
            with self.assertRaisesRegex(ValueError, "Input binding changed"):
                load_export(manifest, trait)

    def test_dry_run_does_not_prepare_or_touch_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/"absent"
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                main(["all", "--analysis-root", str(target), "--dry-run"])
            plan = json.loads(out.getvalue())
            self.assertFalse(plan["fits_regression"])
            self.assertFalse(plan["trains_prsformer"])
            self.assertFalse(target.exists())

    def test_production_adapter_one_and_four_map_binding_and_safe_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            prepared, marker = folder/"prepared.joblib", folder/"blocks.complete"
            prepared.write_bytes(b"mocked trusted source")
            marker.write_text("complete")
            d = pd.DataFrame(dict(eid=[f"{i:03d}" for i in range(8)], family_id=[f"g{i}" for i in range(8)],
                                  split=["train"]*6+["test"]*2, y=np.arange(8.)))
            args = SimpleNamespace(analysis_root=str(folder), prepared=str(prepared), device="cpu",
                                   ld_block_files=None, ld_block_dir=None)
            cfg = {"block_config": str(folder/"config.json")}
            bc = {"ld_block_dir": "mock", "ld_block_files": "mock", "prsformer_cache": str(folder)}
            receipt = dict(array_binding={"array_sha256": "bound"}, shared_weight_fusion={"alpha": {"EUR": .5}})
            for maps in (["ALL"]*4, ["AFR", "EAS", "EUR", "SAS"]):
                tokens = pd.DataFrame(dict(block_id=[f"b{i}" for i in range(4)], source_population=maps,
                                          CHR=[1]*4, start=np.arange(4), end=np.arange(4)+1))
                with patch("twin.load_source", return_value=(d, {}, {"evolution": pd.DataFrame({"private": [1]})})), \
                     patch("abm_core.grid_inputs.locate_blocks", return_value=(np.ones((8, 4, 1)), d, tokens, receipt, folder/"blocks.npy")), \
                     patch("kinship.correct_grid", side_effect=lambda frame, metadata, root: (frame, metadata)), \
                     patch("kinship.publication_metadata", return_value={}), \
                     patch("block_inputs.read_config", return_value=bc), \
                     patch("block_inputs.resolve_block_files", return_value={}), \
                     patch("block_inputs.cache_paths", return_value=(None, None, None, None, marker)):
                    frame, array, rows, features, task, metadata = load_production(args, cfg, "height")
                self.assertEqual(list(features.map_id), maps)
                self.assertEqual(array.shape, (8, 4, 1))
                self.assertEqual(task, "continuous")
                self.assertNotIn("evolution", metadata)
                self.assertFalse(metadata["groups_verified"])
                self.assertEqual(metadata["block_weight_fusion"], receipt["shared_weight_fusion"])
                json.dumps(metadata)


if __name__ == "__main__":
    unittest.main()
