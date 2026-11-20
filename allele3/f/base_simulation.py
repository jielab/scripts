#!/usr/bin/env python3
"""Four-allele mechanism simulations; these are NOT fitted human selection tests.

Alleles are nucleotide identities A/C/G/T. An externally specified introduction
is a mutation/migration pulse, not a response to two resident frequencies.
Genotype fitness, frequency feedback, mutation and drift are separate inputs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from analysis_integrity import transactional, identity

BASES = "ACGT"


def additive_matrix(values):
    a = np.asarray(values, dtype=float)
    return (a[:, None] + a[None, :]) / 2


def builtin_scenarios():
    common = dict(initial=[.51, .49, 0, 0],
                  introductions=[dict(generation=200, base="G", frequency=.001),
                                 dict(generation=400, base="T", frequency=.001)])
    def item(name, description, **settings):
        return dict(common, name=name, description=description, **settings)
    over = np.full((4, 4), .75)
    over[0, 0], over[1, 1], over[0, 1], over[1, 0] = .9, .8, 1., 1.
    # The same antagonistic fitness table is used under hard and soft selection.
    environments = [dict(weight=.55, fitness=additive_matrix([1.2, .8, .85, .8]).tolist()),
                    dict(weight=.45, fitness=additive_matrix([.8, 1.2, .85, .8]).tolist())]
    return [
        item("neutral", "Neutral drift; no restoring force", fitness=np.ones((4, 4)).tolist()),
        item("directional", "Constant net fitness advantage", fitness=additive_matrix([1.04, 1, .97, .96]).tolist()),
        item("positive_feedback", "Common alleles receive a frequency advantage", frequency_feedback=.4),
        item("negative_feedback", "Rare alleles receive a frequency advantage", frequency_feedback=-.4),
        item("overdominance", "A/C equilibrium does not protect a less-fit G/T entrant", fitness=over.tolist()),
        item("pleiotropy_directional", "Multiple concordant effects can still produce fixation",
             trait_effects=[[.025, .015, .01], [0, 0, 0], [-.01, -.01, 0], [-.01, -.01, -.01]]),
        item("antagonism_hard", "Opposing effects with globally weighted reproductive output",
             environments=environments, environment_mode="hard"),
        item("antagonism_soft", "Same opposing effects with fixed habitat contributions",
             environments=environments, environment_mode="soft"),
    ]


def probability_vector(value, label):
    a = np.asarray(value, dtype=float)
    if a.shape != (4,) or not np.all(np.isfinite(a)) or np.any(a < 0) or not np.isclose(a.sum(), 1):
        raise ValueError(label + " must be four nonnegative probabilities summing to 1")
    return a / a.sum()


def fitness_matrix(value):
    a = np.asarray(value, dtype=float)
    if a.shape != (4, 4) or not np.all(np.isfinite(a)) or np.any(a <= 0) or not np.allclose(a, a.T):
        raise ValueError("Fitness must be a positive, finite, symmetric 4x4 genotype matrix")
    return a


def validate_scenario(raw, generations):
    s = dict(raw)
    if not s.get("name") or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in s["name"]):
        raise ValueError("Scenario name must be a simple unique identifier")
    s["initial_array"] = probability_vector(s.get("initial", [.51, .49, 0, 0]), "initial")
    specified = sum(k in s for k in ("fitness", "trait_effects", "environments"))
    if specified > 1:
        raise ValueError("Choose fitness, trait_effects or environments; do not combine them")
    if "trait_effects" in s:
        effects = np.asarray(s["trait_effects"], dtype=float)
        if effects.ndim != 2 or effects.shape[0] != 4 or effects.shape[1] < 1 or not np.all(np.isfinite(effects)):
            raise ValueError("trait_effects must contain four rows of finite log-fitness effects")
        # An explicit toy assumption, not a biological claim about GWAS beta values.
        s["matrices"] = [fitness_matrix(additive_matrix(np.exp(effects.sum(axis=1))))]
        s["weights"] = np.array([1.])
    elif "environments" in s:
        env = s["environments"]
        if not isinstance(env, list) or not env:
            raise ValueError("environments must be a nonempty list")
        s["matrices"] = [fitness_matrix(e["fitness"]) for e in env]
        w = np.asarray([e["weight"] for e in env], float)
        if not np.all(np.isfinite(w)) or np.any(w <= 0) or not np.isclose(w.sum(), 1):
            raise ValueError("Environment weights must be positive and sum to 1")
        s["weights"] = w / w.sum()
    else:
        s["matrices"] = [fitness_matrix(s.get("fitness", np.ones((4, 4))))]
        s["weights"] = np.array([1.])
    s["environment_mode"] = s.get("environment_mode", "hard")
    if s["environment_mode"] not in {"hard", "soft"}:
        raise ValueError("environment_mode must be hard or soft")
    feedback = float(s.get("frequency_feedback", 0))
    if not math.isfinite(feedback) or abs(feedback) > 100:
        raise ValueError("frequency_feedback must be finite in [-100,100]")
    s["frequency_feedback"] = feedback
    mutation = np.asarray(s.get("mutation", np.eye(4)), float)
    if mutation.shape != (4, 4) or not np.all(np.isfinite(mutation)) or np.any(mutation < 0) or not np.allclose(mutation.sum(axis=1), 1):
        raise ValueError("Mutation matrix rows (source->destination) must sum to 1")
    s["mutation_array"] = mutation
    pulses = {}
    for pulse in s.get("introductions", []):
        g, b, f = pulse["generation"], pulse["base"], float(pulse["frequency"])
        if not isinstance(g, int) or g < 1 or not isinstance(b, str) or len(b) != 1 or b not in BASES or not math.isfinite(f) or not 0 < f < 1:
            raise ValueError("Invalid introduction generation/base/frequency")
        if g in pulses:
            raise ValueError("At most one introduction is allowed per generation")
        pulses[g] = (BASES.index(b), f)
    s["pulses"] = pulses
    return s


def expected_after_selection(p, scenario):
    """Random mating, diploid viability selection, then explicit frequency feedback."""
    marginals = np.array([matrix @ p for matrix in scenario["matrices"]])
    means = marginals @ p
    if scenario["environment_mode"] == "soft":
        q = np.sum(scenario["weights"][:, None] * (p * marginals / means[:, None]), axis=0)
    else:
        marginal = scenario["weights"] @ marginals
        q = p * marginal / (p @ marginal)
    # This is a separate viability component with an explicit feedback sign.
    q *= np.exp(scenario["frequency_feedback"] * p)
    q /= q.sum()
    return q


def rare_invasion_ratio(p, index, scenario):
    """Local deterministic growth ratio when an allele is infinitesimally rare.

    Mutation, introduction and drift are intentionally excluded. A ratio > 1
    permits initial invasion in this toy model, not guaranteed establishment.
    """
    p = np.array(p, float)
    p[index] = 0
    if p.sum() <= 0:
        return None
    p /= p.sum()
    marginal = np.array([matrix @ p for matrix in scenario["matrices"]])
    means = marginal @ p
    if scenario["environment_mode"] == "soft":
        factors = np.sum(scenario["weights"][:, None] * marginal / means[:, None], axis=0)
    else:
        m = scenario["weights"] @ marginal
        factors = m / (p @ m)
    factors *= np.exp(scenario["frequency_feedback"] * p)
    return float(factors[index] / (p @ factors))


def simulate(scenario, generations, population_size, rng, record_every):
    p = scenario["initial_array"].copy()
    rows, pulses = [], []
    for generation in range(generations + 1):
        event = ""
        if generation:
            if generation in scenario["pulses"]:
                index, frequency = scenario["pulses"][generation]
                ratio = rare_invasion_ratio(p, index, scenario)
                p *= 1 - frequency
                p[index] += frequency
                event = "external_introduction_" + BASES[index]
                pulses.append(dict(generation=generation, base=BASES[index], introduced_fraction=frequency,
                                   rare_invasion_ratio=ratio, actual_post_pulse_frequency=float(p[index])))
            q = expected_after_selection(p, scenario) @ scenario["mutation_array"]
            q /= q.sum()
            p = rng.multinomial(2 * population_size, q) / (2 * population_size) if population_size else q
        if generation % record_every == 0 or generation == generations or event:
            rows.append(dict(generation=generation, event=event,
                             **{b: float(p[i]) for i, b in enumerate(BASES)},
                             heterozygosity=float(1 - p @ p)))
    return rows, pulses, p


def write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as f:
        writer = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if row.get(k) is None else row.get(k, "") for k in fields})


@transactional("qc/manifest.json", input_files=("config",))
def run(args):
    if args.generations < 1 or args.population_size < 0 or args.replicates < 1 or args.record_every < 1:
        raise ValueError("Invalid generations, population size, replicates or recording interval")
    if not 0 < args.persistence_frequency < .25:
        raise ValueError("persistence-frequency must be in (0,.25)")
    before = identity(args.config, full=True) if args.config else None
    raw = json.loads(Path(args.config).read_text())["scenarios"] if args.config else builtin_scenarios()
    if args.scenarios:
        chosen = set(args.scenarios.split(","))
        known = {s["name"] for s in raw}
        if chosen - known:
            raise ValueError("Unknown scenarios: " + ",".join(sorted(chosen - known)))
        raw = [s for s in raw if s["name"] in chosen]
    if not raw or len({s["name"] for s in raw}) != len(raw):
        raise ValueError("Scenarios must be nonempty and uniquely named")
    checked = [validate_scenario(s, args.generations) for s in raw]
    out = Path(args.out)
    for folder in ("results", "reports", "qc"):
        (out/folder).mkdir(parents=True, exist_ok=True)
    seed = np.random.SeedSequence(args.seed)
    seeds = iter(seed.spawn(len(raw) * args.replicates))
    summaries, introductions, plotting = [], [], {}
    path = out / "results/trajectories.tsv"
    fields = ["scenario", "replicate", "generation", "event", *BASES, "heterozygosity"]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n"); writer.writeheader()
        for scenario in checked:
            for replicate in range(args.replicates):
                rows, pulses, p = simulate(scenario, args.generations, args.population_size,
                                          np.random.default_rng(next(seeds)), args.record_every)
                for row in rows:
                    writer.writerow(dict(scenario=scenario["name"], replicate=replicate, **row))
                introductions += [dict(scenario=scenario["name"], replicate=replicate, **r) for r in pulses]
                summaries.append(dict(scenario=scenario["name"], replicate=replicate,
                    **{b+"_final": float(p[i]) for i, b in enumerate(BASES)},
                    n_above_persistence_threshold=int((p >= args.persistence_frequency).sum()),
                    persistence_frequency=args.persistence_frequency,
                    fixation_at_least_99pct=int(p.max() >= .99),
                    interpretation="mechanism_simulation_not_human_selection_inference"))
                if replicate == 0:
                    plotting[scenario["name"]] = rows
    write_tsv(out/"results/summary.tsv", summaries, list(summaries[0]))
    write_tsv(out/"results/introductions.tsv", introductions, ["scenario", "replicate", "generation", "base", "introduced_fraction", "rare_invasion_ratio", "actual_post_pulse_frequency"])
    (out/"qc/scenarios.json").write_text(json.dumps(dict(scenarios=raw), indent=2)+"\n")
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = len(plotting); fig, axes = plt.subplots(math.ceil(n/2), 2, figsize=(13, 3.1*math.ceil(n/2)), squeeze=False)
        colors = ["#2864A2", "#D86542", "#2D956D", "#9160AC"]
        for ax, (name, rows) in zip(axes.flat, plotting.items()):
            for b, color in zip(BASES, colors):
                ax.plot([r["generation"] for r in rows], [r[b] for r in rows], label=b, color=color)
            ax.set(title=name.replace("_", " "), ylim=(-.02, 1.02), xlabel="Generation", ylabel="Allele frequency")
            ax.spines[["top", "right"]].set_visible(False)
        for ax in list(axes.flat)[n:]: ax.set_visible(False)
        axes.flat[0].legend(ncol=4)
        fig.suptitle("Four-allele mechanism simulations (replicate 0; illustrative parameters)", fontsize=14)
        fig.tight_layout(); fig.savefig(out/"reports/base_simulation.png", dpi=180); plt.close(fig)
    if before and before != identity(args.config, full=True):
        raise ValueError("Simulation config changed while running")
    manifest = dict(status="complete", settings=vars(args), input=before,
                    hypotheses="Fixed net advantage, frequency feedback, dominance and environmental trade-offs are distinct mechanisms",
                    model="diploid random mating; viability selection; explicit feedback; mutation; optional Wright-Fisher drift",
                    interpretation="Simulation only; no estimated human parameters, no proof of pleiotropy or balancing selection",
                    introduction="Scheduled external mutation/migration pulses; never triggered by 50:50 frequencies",
                    unobserved_future_pulses="Scheduled introductions after --generations remain in scenarios.json and do not occur in this run",
                    assumptions="trait_effects are additive log-fitness contributions for a toy model, never GWAS or AVI values",
                    population_size="diploid Ne; 0 means deterministic infinite-population recursion",
                    output_rows=len(summaries))
    (out/"qc/manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", help="Optional JSON with scenarios; defaults compare eight mechanisms")
    p.add_argument("--scenarios", help="Comma-separated subset of scenario names")
    p.add_argument("--out", required=True)
    p.add_argument("--generations", type=int, default=1000)
    p.add_argument("--population-size", type=int, default=0, help="Diploid Ne; 0 for deterministic recursion")
    p.add_argument("--replicates", type=int, default=1)
    p.add_argument("--record-every", type=int, default=5)
    p.add_argument("--seed", type=int, default=20261110)
    p.add_argument("--persistence-frequency", type=float, default=.01)
    p.add_argument("--plot", action="store_true")
    return p


def main(argv=None):
    result = run(parser().parse_args(argv))
    print(json.dumps(dict(status=result["status"], output_rows=result["output_rows"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
