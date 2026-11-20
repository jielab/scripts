"""Synthetic, non-UKB smoke-test cohort; never a reproduction result."""
from __future__ import annotations

import copy
import csv
import datetime as dt
import json
from pathlib import Path

import numpy as np


def make_demo(config: dict, n_people: int = 320) -> dict:
    """Write canonical TSV inputs and return their fully resolved configuration.

    Latent variables generate correlated dated observations and ICD-10 events.
    This small cohort exercises split isolation, delayed measurements, competing
    death and five-year labels. It has no patient information or scientific
    performance interpretation.
    """
    cfg = copy.deepcopy(config)
    if n_people < 60:
        raise ValueError("Use at least 60 synthetic people so all three splits are populated")
    seed = int(cfg.get("seed", 2026))
    rng = np.random.default_rng(seed)
    root = Path(cfg.setdefault("paths", {}).get("out", "demo_out")).resolve()
    inputs = root / "data" / "demo_inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    cfg["paths"]["ukb"] = str(inputs)
    d = cfg.setdefault("data", {})
    d.update(min_token_count=2, min_feature_count=2,
             use_input_split=True, allow_first_occurrence=False,
             endpoints={"cad": ["I25"], "t2d": ["E11"], "ckd": ["N18"], "death": ["DEATH"]})
    d["sources"] = {"participants": str(inputs / "participants.tsv"),
                    "events": str(inputs / "events.tsv"),
                    "observations": str(inputs / "observations.tsv")}
    d.setdefault("landmark", "midpoint")
    d.setdefault("min_clinical_days", 3)
    marker = inputs / "SYNTHETIC.json"
    if marker.exists():
        try:
            previous = json.loads(marker.read_text())
        except (ValueError, OSError):
            previous = {}
        if (previous.get("synthetic") is True and previous.get("seed") == seed
                and previous.get("n_people") == n_people
                and previous.get("generator_version") == 2
                and all((inputs / (name + ".tsv")).is_file() for name in ("participants", "events", "observations"))):
            print(f"[flux] reusing synthetic demo: {n_people} people in {inputs}", flush=True)
            return cfg
    files = {name: (inputs / (name + ".tsv")).open("w", newline="")
             for name in ("participants", "events", "observations")}
    writers = {name: csv.writer(handle, delimiter="\t", lineterminator="\n")
               for name, handle in files.items()}
    writers["participants"].writerow(["eid", "birth_date", "baseline_date", "exit_date", "death_date", "sex", "PC1", "split"])
    writers["events"].writerow(["eid", "date", "code"])
    writers["observations"].writerow(["eid", "date", "modality", "feature", "value"])
    try:
        for i in range(n_people):
            eid = str(9000001 + i)
            z, renal, cardio = rng.normal(size=3)
            sex = int(rng.integers(0, 2))
            baseline = dt.date(2010, 1, 1) + dt.timedelta(days=int(rng.integers(-180, 181)))
            birth = baseline - dt.timedelta(days=int((50 + 12 * rng.random()) * 365.25))
            death = None
            if rng.random() < .18:
                # Early deaths also exercise a positive five-year death label:
                # their shortened EHR histories move the ordinal midpoint back.
                death_days = int(rng.integers(30, 330)) if rng.random() < .6 else int(rng.integers(7 * 365, 12 * 365))
                death = baseline + dt.timedelta(days=death_days)
            # Ensure the small smoke-test splits exercise the mortality head.
            # This artificial stratification is another reason demo metrics
            # must never be interpreted as scientific validation.
            if i % 60 in {2, 14, 17}:
                death = baseline + dt.timedelta(days=180)
            split = "train" if i % 20 < 14 else "validation" if i % 20 < 17 else "test"
            writers["participants"].writerow([eid, birth, baseline, d["admin_end"], death or "", sex, round(float(z + rng.normal(scale=.2)), 6), split])

            # Ten distinct EHR days, with some same-day comorbid codes. A fixed
            # ordinal midpoint remains different from the recruitment visit.
            offsets = [-3650, -2555, -1460, -365, 730, 1460, 2190, 2920, 3650, 4380]
            outcome = {"I25": cardio + .5 * z, "E11": z, "N18": renal + .25 * z}
            future_codes = {}
            for code, score in outcome.items():
                probability = 1 / (1 + np.exp(-1.2 * score))
                if rng.random() < probability:
                    slot = int(rng.integers(6, 10))
                    future_codes.setdefault(slot, []).append(code)
            for j, offset in enumerate(offsets):
                when = baseline + dt.timedelta(days=offset)
                if death is not None and when > death:
                    continue
                background = ["J06", "M54", "K21", "I10", "E78"][int(rng.integers(0, 5))]
                writers["events"].writerow([eid, when, background])
                for code in future_codes.get(j, []):
                    writers["events"].writerow([eid, when, code])
            # Explicit prevalent disease exercises incident-case exclusion.
            if rng.random() < .08:
                writers["events"].writerow([eid, baseline - dt.timedelta(days=100), "E11"])

            visits = [baseline, baseline + dt.timedelta(days=3 * 365), baseline + dt.timedelta(days=8 * 365)]
            for j, when in enumerate(visits):
                if death is not None and when > death:
                    continue
                age_effect = j * .1
                values = [("clinical", "bmi", 27 + 3 * z + rng.normal()),
                          ("clinical", "sbp", 125 + 10 * cardio + 3 * sex + rng.normal(scale=3)),
                          ("bbc", "crp", np.exp(.3 * z + age_effect + rng.normal(scale=.15))),
                          ("bbc", "creatinine", 75 + 12 * renal + rng.normal(scale=2)),
                          ("met", "glucose", 5.5 + .7 * z + rng.normal(scale=.2)),
                          ("prot", "GDF15", .6 * renal + age_effect + rng.normal(scale=.2)),
                          ("prot", "APOE", .5 * cardio + rng.normal(scale=.3))]
                for modality, feature, value in values:
                    writers["observations"].writerow([eid, when, modality, feature, round(float(value), 6)])
            # Half the imaging is after the risk landmark and must stay hidden.
            image_day = baseline + dt.timedelta(days=(2 if i % 2 else 8) * 365)
            if death is None or image_day <= death:
                writers["observations"].writerow([eid, image_day, "img", "brain_volume", round(float(1450 - 30 * renal + rng.normal(scale=12)), 6)])
    finally:
        for handle in files.values():
            handle.close()
    marker.write_text(json.dumps({"synthetic": True, "generator_version": 2, "seed": seed, "n_people": n_people,
        "warning": "Smoke test only; no real UKB data and no scientific benchmark claim."}, indent=2) + "\n")
    print(f"[flux] synthetic demo: {n_people} people in {inputs}", flush=True)
    return cfg
