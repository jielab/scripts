"""Dated UKB adapters and leakage-aware sequences for flux.

Raw files are read in a stream and consolidated in SQLite. Model preprocessing
uses training participants only. No processed/imputed Rdata/*.rds is consumed.
The paired UKB 41270/41280 export contains FIRST occurrences, not encounters.
Its explicitly enabled fallback is recorded in meta.json and the run report.
"""
from __future__ import annotations

import bisect
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
from typing import Any, Iterable
import warnings

import numpy as np
from common import resolve_followup, validate_record

EPOCH = dt.date(1970, 1, 1)
MISSING = {"", "NA", "N/A", "NAN", "NULL", "NONE", "."}
BAD_DATES = {"1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"}
CLINICAL_FIELDS = ["p21001", "p21002", "p48", "p49", "p50", "p46", "p47", "p4080", "p4079", "p102", "p874", "p884", "p904"]
DEFAULT_ENDPOINTS = {"cvd_cad": ["I20", "I21", "I22", "I23", "I24", "I25"], "cvd_hfail": ["I50"], "stroke_isch": ["I63"], "t2d": ["E11"], "ckd": ["N18"], "death": ["DEATH"]}


def day(value: Any) -> int | None:
    """ISO date or integer epoch day. Sentinel and ambiguous dates are missing."""
    if value is None:
        return None
    if isinstance(value, (int, np.integer)):
        return int(value)
    s = str(value).strip()
    if s.upper() in MISSING or s[:10] in BAD_DATES:
        return None
    # ISO timestamps are permitted, while locale-dependent day/month parsing is not.
    if not re.match(r"^\d{4}-\d{2}-\d{2}(?:$|[ T])", s):
        return None
    try:
        return (dt.date.fromisoformat(s[:10]) - EPOCH).days
    except ValueError:
        return None


def date_string(value: int | None) -> str:
    return "" if value is None else (EPOCH + dt.timedelta(days=int(value))).isoformat()


def number(value: Any, nonnegative: bool = False) -> float | None:
    try:
        v = float(value)
    except (ValueError, TypeError):
        return None
    if not math.isfinite(v) or (nonnegative and v < 0):
        return None
    return v


def eid_string(value: Any) -> str:
    s = str(value or "").strip()
    return re.sub(r"\.0$", "", s) if re.fullmatch(r"\d+\.0", s) else s


def code3(value: Any) -> str | None:
    s = str(value or "").upper().replace(".", "").replace(" ", "")
    if s == "DEATH":
        return s
    return s[:3] if re.fullmatch(r"[A-Z]\d{2}[A-Z0-9]*", s) else None


def field_parts(name: str) -> tuple[str, int, int | None] | None:
    """Understand RAP p53_i0, p41280_a0 and legacy 53-0.0/f.53.0.0."""
    s = str(name)
    m = re.fullmatch(r"p?(\d+)(?:_i(\d+))?(?:_a(\d+))?", s)
    if m:
        return "p" + m[1], int(m[2] or 0), int(m[3]) if m[3] is not None else None
    m = re.fullmatch(r"(?:f\.)?(\d+)[-.](\d+)\.(\d+)", s)
    if m:
        return "p" + m[1], int(m[2]), int(m[3])
    return None


def open_text(path: str | Path):
    p = Path(path)
    return gzip.open(p, "rt", encoding="utf-8-sig", newline="") if p.suffix == ".gz" else p.open("r", encoding="utf-8-sig", newline="")


def rows(path: str | Path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No TSV header: {path}")
        reader.fieldnames = [x.strip() for x in reader.fieldnames]
        for i, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Malformed TSV row {i}: {path}")
            yield row


def header(path: str | Path) -> list[str]:
    with open_text(path) as f:
        return next(csv.reader(f, delimiter="\t"), [])


def choose(row: dict, names: Iterable[str], default=None):
    for k in names:
        if k in row and str(row[k]).strip().upper() not in MISSING:
            return row[k]
    return default


def field_value(row: dict, field: str, instance: int = 0):
    aliases = [f"{field}_i{instance}", f"{field}_i{instance}_a0"]
    if instance == 0:
        aliases = [field] + aliases
    result = choose(row, aliases)
    if result is not None:
        return result
    for k, value in row.items():
        parts = field_parts(k)
        if parts and parts[0] == field and parts[1] == instance and parts[2] in (None, 0):
            return value
    return None


def split_for(eid: str, seed: int, fractions: list[float]) -> str:
    if len(fractions) != 3 or not math.isclose(sum(fractions), 1.0, abs_tol=1e-8) or min(fractions) <= 0:
        raise ValueError("data.split must contain three positive fractions summing to 1")
    u = int.from_bytes(hashlib.blake2b(f"{seed}:{eid}".encode(), digest_size=8).digest(), "big") / 2**64
    return "train" if u < fractions[0] else "validation" if u < fractions[0] + fractions[1] else "test"


def source_path(root: Path, setting: Any, defaults: Iterable[str], required=False) -> Path | None:
    if setting is None or str(setting).lower() in {"none", "false", "disabled"}:
        if required:
            raise ValueError("Required input was disabled")
        return None
    if str(setting).lower() != "auto":
        p = Path(os.path.expandvars(str(setting))).expanduser()
        p = p if p.is_absolute() else root / p
        if not p.is_file():
            raise FileNotFoundError(p)
        return p
    for rel in defaults:
        p = root / rel
        if p.is_file():
            return p
    if required:
        raise FileNotFoundError(f"No input found under {root}; tried: {', '.join(defaults)}")
    return None


def linked_defaults(stem: str) -> list[str]:
    return [f"{folder}/{stem}{suffix}" for folder in ("rap/raw", "raw", "rap", ".") for suffix in (".tsv.gz", ".tsv", ".tab.gz", ".tab")]


def _paths(config):
    paths = config.get("paths", {})
    ukb = Path(paths.get("ukb", config.get("dir_ukb", os.environ.get("DIR_UKB", "/mnt/d/phe/ukb")))).expanduser()
    out = Path(paths.get("out", config.get("out_dir", "analysis/flux"))).expanduser()
    return ukb, out


SCHEMA = """
CREATE TABLE participants(eid TEXT PRIMARY KEY,birth_day INTEGER NOT NULL,baseline_day INTEGER NOT NULL,exit_day INTEGER NOT NULL,death_day INTEGER,split TEXT NOT NULL,context TEXT NOT NULL);
CREATE TABLE visits(eid TEXT,instance INTEGER,day INTEGER,PRIMARY KEY(eid,instance));
CREATE TABLE events(eid TEXT,day INTEGER,code TEXT,PRIMARY KEY(eid,day,code)) WITHOUT ROWID;
CREATE TABLE observations(eid TEXT,day INTEGER,modality TEXT,values_json TEXT,PRIMARY KEY(eid,day,modality)) WITHOUT ROWID;
"""


class Importer:
    def __init__(self, db, config, root):
        self.db, self.config, self.root = db, config, root
        self.d = config.get("data", {})
        self.sources = self.d.get("sources", {})
        self.followup_provenance = resolve_followup()
        validate_record(config, "data", "admin_end", self.followup_provenance)
        if self.d.get("admin_end") != self.followup_provenance["date"]:
            raise ValueError("Load configuration with the shared 0f/phenotype.R follow-up date first")
        self.admin_end = day(self.followup_provenance["date"])
        self.seed = int(self.d.get("seed", config.get("seed", 2026)))
        self.fractions = list(self.d.get("split", [0.7, 0.15, 0.15]))
        self.people = {}
        self.visits = {}
        self.count = collections.Counter()
        self.source_info = []
        self.messages = []
        self.event_mode = "unknown"
        self.feature_rules = {}

    def find(self, key, defaults, required=False):
        p = source_path(self.root, self.sources.get(key, "auto"), defaults, required)
        if p:
            stat = p.stat()
            self.source_info.append({"name": key, "path": str(p.resolve()), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})
        return p

    def warn(self, message):
        self.messages.append(message)
        print("[flux] " + message, flush=True)

    def add_visits(self, eid, row):
        for key, value in row.items():
            parts = field_parts(key)
            if parts and parts[0] == "p53":
                d = day(value)
                if d is not None:
                    self.visits[(eid, parts[1])] = d

    def participants(self):
        canonical = self.find("participants", ["participants.tsv.gz", "participants.tsv"])
        vip = canonical or self.find("vip", ["rap/vip.tab.gz", "rap/vip.tab"], required=True)
        if not canonical and self.admin_end is None:
            raise ValueError("Missing resolved 0f/phenotype.R::date_follow_end; load the shared configuration first.")
        for row in rows(vip):
            eid = eid_string(choose(row, ["eid", "participant.eid"]))
            if not eid or eid.upper() in MISSING:
                self.count["participant_missing_eid"] += 1
                continue
            if eid in self.people:
                raise ValueError(f"Duplicate participant eid: {eid}")
            birth = day(choose(row, ["birth_date", "date_birth"]))
            if birth is None:
                year = number(choose(row, ["birth_year"], field_value(row, "p34")))
                month = number(choose(row, ["birth_month"], field_value(row, "p52")))
                if year is not None and month is not None:
                    birth = day(f"{int(year):04d}-{int(month):02d}-15")
                    self.count["birth_midmonth_proxy"] += 1
            baseline = day(choose(row, ["baseline_date", "date_attend", "date_assess"], field_value(row, "p53")))
            death = day(choose(row, ["death_date", "date_death"], field_value(row, "p40000")))
            loss = day(choose(row, ["date_lost", "loss_date"], field_value(row, "p191")))
            exit_day = day(choose(row, ["exit_date", "admin_end", "followup_end"]))
            limits = [d for d in (exit_day, loss, self.admin_end) if d is not None]
            exit_day = min(limits) if limits else None
            if birth is None or baseline is None or exit_day is None or baseline <= birth:
                self.count["participant_invalid_dates"] += 1
                continue
            context = {}
            sex = number(choose(row, ["sex"], field_value(row, "p31")))
            if sex in (0, 1):
                context["sex"] = sex
            for key, value in row.items():
                if re.fullmatch(r"PC\d+|context_[A-Za-z0-9_]+", key):
                    v = number(value)
                    if v is not None:
                        context[key] = v
            split = split_for(eid, self.seed, self.fractions)
            if canonical and self.d.get("use_input_split", False):
                supplied = str(row.get("split", "")).lower()
                if supplied not in {"train", "validation", "test"}:
                    raise ValueError("use_input_split requires train, validation or test for every participant")
                split = supplied
            self.people[eid] = (birth, baseline, exit_day, death, split, context)
            self.visits[(eid, 0)] = baseline
            self.add_visits(eid, row)
        if not self.people:
            raise ValueError("No valid participants; inspect birth, baseline and administrative dates")
        # Death and PC exports enrich immutable participant context before event loading.
        death_file = self.find("death", ["rap/death.tab.gz", "rap/death.tab"])
        if death_file:
            for row in rows(death_file):
                eid = eid_string(choose(row, ["eid", "participant.eid"]))
                if eid not in self.people:
                    continue
                dates = [day(v) for k, v in row.items() if k in {"date_death", "death_date"} or (field_parts(k) and field_parts(k)[0] == "p40000")]
                dates = [v for v in dates if v is not None]
                if dates:
                    p = list(self.people[eid]); p[3] = min(dates + ([p[3]] if p[3] is not None else [])); self.people[eid] = tuple(p)
        pca = self.find("pca", ["rap/pca.tab.gz", "rap/pca.tab"])
        if pca:
            for row in rows(pca):
                eid = eid_string(choose(row, ["eid", "participant.eid"]))
                if eid not in self.people:
                    continue
                for key, value in row.items():
                    parts = field_parts(key)
                    is_pc = re.fullmatch(r"PC[1-9]\d*", key)
                    if not is_pc and not (parts and parts[0] == "p22009" and parts[2] is not None):
                        continue
                    v = number(value)
                    if v is not None:
                        label = key if is_pc else "PC" + str(parts[2] + 1)
                        self.people[eid][5][label] = v
        self.db.executemany("INSERT INTO participants VALUES (?,?,?,?,?,?,?)", ((eid, *p[:5], json.dumps(p[5], separators=(",", ":"))) for eid, p in self.people.items()))
        self.db.executemany("INSERT OR REPLACE INTO visits VALUES (?,?,?)", ((eid, instance, d) for (eid, instance), d in self.visits.items()))
        self.db.commit()
        self.count["participants"] = len(self.people)
        return canonical, vip

    def valid_date(self, eid, d, role):
        if eid not in self.people:
            self.count[role + "_unknown_eid"] += 1
            return False
        if d is None:
            self.count[role + "_missing_date"] += 1
            return False
        birth, _, end, death, _, _ = self.people[eid]
        if d < birth or d > end or (death is not None and d > death):
            self.count[role + "_outside_life_or_followup"] += 1
            return False
        return True

    def add_events(self, batch):
        if batch:
            self.db.executemany("INSERT OR IGNORE INTO events VALUES (?,?,?)", batch)
            batch.clear()

    def events(self):
        long = self.find("events", ["events.tsv.gz", "events.tsv"])
        if long:
            self.event_mode = "dated_long_events"
            batch = []
            death_events = {}
            for row in rows(long):
                eid = eid_string(choose(row, ["eid", "participant.eid"]))
                d = day(choose(row, ["date", "event_date"]))
                code = code3(choose(row, ["code", "icd10", "diag_icd10"]))
                if code and self.valid_date(eid, d, "event"):
                    batch.append((eid, d, code))
                    if code == "DEATH":
                        death_events[eid] = min(d, death_events.get(eid, d))
                if len(batch) >= 20000:
                    self.add_events(batch)
            self.add_events(batch)
            # Canonical DEATH rows define the same life boundary as death_date.
            for eid, d in death_events.items():
                person = list(self.people[eid])
                if person[3] is None or d < person[3]:
                    person[3] = d
                    self.people[eid] = tuple(person)
                    self.db.execute("UPDATE participants SET death_day=? WHERE eid=?", (d, eid))
                    removed = self.db.execute("DELETE FROM events WHERE eid=? AND day>?", (eid, d))
                    self.count["event_outside_life_or_followup"] += removed.rowcount
        else:
            hes = self.find("hesin", linked_defaults("hesin"))
            diag = self.find("hesin_diag", linked_defaults("hesin_diag"))
            if hes and diag:
                self.event_mode = "hospital_episode_records"
                self.hospital(hes, diag)
            elif self.d.get("allow_first_occurrence", False):
                c = self.find("icd10_codes", ["rap/icd10.code.tab.gz", "rap/icd10.code.tab"])
                t = self.find("icd10_dates", ["rap/icd10.date.tab.gz", "rap/icd10.date.tab"])
                if not c or not t:
                    raise FileNotFoundError("Need hesin + hesin_diag, canonical events.tsv[.gz], or both rap/icd10.code.tab.gz and rap/icd10.date.tab.gz")
                self.event_mode = "first_occurrence_only"
                self.warn("Using p41270/p41280 first occurrences: repeated clinical encounters are unavailable in this source. This is a limited-data reproduction.")
                self.first_occurrence(c, t)
            else:
                raise FileNotFoundError("Full hospital episodes not found. Set data.sources.hesin/hesin_diag, supply canonical events.tsv, or explicitly enable data.allow_first_occurrence.")
        death_batch = [(eid, p[3], "DEATH") for eid, p in self.people.items() if p[3] is not None and self.valid_date(eid, p[3], "death")]
        self.add_events(death_batch)
        self.db.commit()
        self.count["events"] = self.db.execute("SELECT count(*) FROM events").fetchone()[0]
        if self.count["events"] == 0:
            raise ValueError("No valid dated clinical events")

    def hospital(self, hes, diag):
        self.db.execute("CREATE TABLE episodes(eid TEXT,ins_index TEXT,day INTEGER,PRIMARY KEY(eid,ins_index)) WITHOUT ROWID")
        batch = []
        for row in rows(hes):
            eid = eid_string(choose(row, ["eid", "participant.eid"]))
            idx = choose(row, ["ins_index"])
            d = next((d for d in (day(row.get("epistart")), day(row.get("admidate"))) if d is not None), None)
            if idx is not None and self.valid_date(eid, d, "episode"):
                batch.append((eid, str(idx), d))
            if len(batch) >= 20000:
                self.db.executemany("INSERT OR IGNORE INTO episodes VALUES (?,?,?)", batch); batch.clear()
        if batch:
            self.db.executemany("INSERT OR IGNORE INTO episodes VALUES (?,?,?)", batch)
        if "diag_icd10" not in header(diag):
            raise ValueError("hesin_diag requires diag_icd10. ICD-9 must be mapped explicitly before supplying canonical events.")
        batch = []
        previous_key = None; previous_day = None
        for row in rows(diag):
            eid = eid_string(choose(row, ["eid", "participant.eid"]))
            idx = choose(row, ["ins_index"])
            code = code3(row.get("diag_icd10"))
            if eid not in self.people or idx is None or code is None:
                continue
            key = (eid, str(idx))
            if key != previous_key:
                hit = self.db.execute("SELECT day FROM episodes WHERE eid=? AND ins_index=?", key).fetchone()
                previous_day = hit[0] if hit else None
                previous_key = key
            if previous_day is not None:
                batch.append((eid, previous_day, code))
            else:
                self.count["diagnosis_unmatched_episode"] += 1
            if len(batch) >= 20000:
                self.add_events(batch)
        self.add_events(batch)
        self.db.execute("DROP TABLE episodes")

    def first_occurrence(self, codes, dates):
        self.db.execute("CREATE TABLE widecodes(eid TEXT PRIMARY KEY,payload TEXT)")
        batch = []
        for row in rows(codes):
            eid = eid_string(choose(row, ["eid", "participant.eid"]))
            if eid not in self.people:
                continue
            mapped = {}
            for name, value in row.items():
                parts = field_parts(name)
                if not parts or parts[0] != "p41270":
                    continue
                pieces = re.split(r"[|;]", str(value))
                start = parts[2] or 0
                for offset, code in enumerate(pieces):
                    code = code3(code)
                    if code:
                        mapped[f"{parts[1]}:{start + offset}"] = code
            batch.append((eid, json.dumps(mapped, separators=(",", ":"))))
            if len(batch) >= 10000:
                self.db.executemany("INSERT INTO widecodes VALUES (?,?)", batch); batch.clear()
        if batch:
            self.db.executemany("INSERT INTO widecodes VALUES (?,?)", batch)
        batch = []
        for row in rows(dates):
            eid = eid_string(choose(row, ["eid", "participant.eid"]))
            if eid not in self.people:
                continue
            found = self.db.execute("SELECT payload FROM widecodes WHERE eid=?", (eid,)).fetchone()
            if not found:
                continue
            mapping = json.loads(found[0])
            for name, value in row.items():
                parts = field_parts(name)
                if not parts or parts[0] != "p41280":
                    continue
                start = parts[2] or 0
                for offset, raw_date in enumerate(str(value).split("|")):
                    d = day(raw_date)
                    code = mapping.get(f"{parts[1]}:{start + offset}")
                    if code and self.valid_date(eid, d, "first_occurrence"):
                        batch.append((eid, d, code))
            if len(batch) >= 20000:
                self.add_events(batch)
        self.add_events(batch)
        self.db.execute("DROP TABLE widecodes")

    def add_observation(self, eid, d, modality, values):
        if not values or not self.valid_date(eid, d, "observation"):
            return
        payload = json.dumps(values, separators=(",", ":"), allow_nan=False)
        # json_patch preserves independently supplied features on the same visit.
        self.db.execute("INSERT INTO observations VALUES (?,?,?,?) ON CONFLICT(eid,day,modality) DO UPDATE SET values_json=json_patch(observations.values_json,excluded.values_json)", (eid, d, modality, payload))

    def observations(self, vip):
        long = self.find("observations", ["observations.tsv.gz", "observations.tsv"])
        if long:
            # Long imports need not be presorted; SQLite does the external grouping.
            self.db.execute("CREATE TABLE measures(eid TEXT,day INTEGER,modality TEXT,feature TEXT,value REAL,PRIMARY KEY(eid,day,modality,feature)) WITHOUT ROWID")
            batch = []
            for row in rows(long):
                eid = eid_string(choose(row, ["eid", "participant.eid"]))
                d = day(choose(row, ["date", "measurement_date"]))
                modality = str(row.get("modality", "clinical")).strip()
                feature = str(choose(row, ["feature", "variable"], "")).strip()
                value = number(row.get("value"))
                if feature and modality and value is not None and self.valid_date(eid, d, "observation"):
                    batch.append((eid, d, modality, feature, value))
                if len(batch) >= 20000:
                    self.db.executemany("INSERT OR REPLACE INTO measures VALUES (?,?,?,?,?)", batch); batch.clear()
            if batch:
                self.db.executemany("INSERT OR REPLACE INTO measures VALUES (?,?,?,?,?)", batch)
            previous = None; values = {}
            for eid, d, modality, feature, value in self.db.execute("SELECT * FROM measures ORDER BY eid,day,modality"):
                key = (eid, d, modality)
                if previous is not None and key != previous:
                    self.add_observation(*previous, values); values = {}
                values[feature] = value; previous = key
            if previous is not None:
                self.add_observation(*previous, values)
            self.db.execute("DROP TABLE measures")
        else:
            enabled = self.d.get("modalities", ["clinical", "bbc", "prot", "met", "img"])
            paths = {}
            if "clinical" in enabled:
                paths["clinical"] = vip
            for key, defaults in {"bbc": ["rap/bbc.tab.gz"], "prot": ["rap/raw/prot.tab.gz"], "met": ["rap/met.tab.gz"], "img": ["rap/img.tab.gz"]}.items():
                if key in enabled:
                    paths[key] = self.find(key, defaults)
            # Imaging extra exports can hold the dates even if primary IDPs do not.
            extra = self.find("img_extra", ["rap/img.extra.tab.gz"]) if "img" in enabled else None
            for p in [paths.get("img"), extra]:
                if p:
                    for row in rows(p):
                        eid = eid_string(choose(row, ["eid", "participant.eid"]))
                        if eid in self.people:
                            self.add_visits(eid, row)
            for modality, p in paths.items():
                if p:
                    self.wide_observations(p, modality)
                else:
                    self.messages.append(f"Optional modality {modality} not found; skipped")
            if extra:
                existing = set(header(paths["img"])) if paths.get("img") else set()
                self.wide_observations(extra, "img", exclude=existing)
        self.db.executemany("INSERT OR REPLACE INTO visits VALUES (?,?,?)", ((eid, instance, d) for (eid, instance), d in self.visits.items()))
        self.db.commit()
        self.count["observation_visits"] = self.db.execute("SELECT count(*) FROM observations").fetchone()[0]

    def feature_selection(self, modality):
        """Explicit selectors, or verified local UKB field dictionaries."""
        if modality in self.feature_rules:
            return self.feature_rules[modality]
        def configured(key):
            mapping = self.d.get(key, {})
            if not isinstance(mapping, dict):
                raise ValueError(f"data.{key} must map modality names to feature lists")
            if modality not in mapping:
                return None
            values = mapping[modality]
            if isinstance(values, str):
                values = [values]
            if not isinstance(values, list) or not all(isinstance(x, str) for x in values):
                raise ValueError(f"data.{key}.{modality} must be a string list")
            return set(values)
        include, exclude = configured("feature_include"), configured("feature_exclude") or set()
        if include is None and modality in {"bbc", "met"}:
            dictionary = self.root / "common" / (modality + ".lst")
            if dictionary.is_file():
                allowed = set()
                ignored = 0
                with open_text(dictionary) as handle:
                    for record in csv.reader(handle, delimiter="\t"):
                        if not record:
                            continue
                        parts = field_parts(record[0].strip())
                        label = record[1].strip() if len(record) > 1 else ""
                        metadata = re.search(r"(?:^|[_\W])(date|time|batch|plate|well|assay|visit|sample|aliquot|qc|quality|dilution|device)(?:$|[_\W])", label, re.I)
                        if parts and not metadata:
                            allowed.add(parts[0])
                        else:
                            ignored += 1
                if not allowed:
                    raise ValueError(f"No standalone biological UKB field IDs in {dictionary}; inspect the dictionary or supply data.feature_include.{modality}")
                include = allowed
                stat = dictionary.stat()
                self.source_info.append({"name": modality + "_feature_dictionary", "path": str(dictionary.resolve()), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})
                self.messages.append(f"{modality}: restricted raw fields to {len(allowed)} standalone UKB IDs from {dictionary}; {ignored} metadata/derived-expression rows ignored")
            else:
                self.messages.append(f"{modality}: no common/{modality}.lst field dictionary; inspect raw numeric fields and use data.feature_include/feature_exclude to remove coded assay/QC metadata")
        self.feature_rules[modality] = (include, exclude)
        return include, exclude

    def wide_observations(self, path, modality, exclude=None):
        names = header(path)
        columns = []
        clinical = set(self.d.get("clinical_fields", CLINICAL_FIELDS))
        explicitly_dated = self.d.get("date_columns", {}).get(modality)
        include_features, exclude_features = self.feature_selection(modality)
        for name in names:
            if name in {"eid", "participant.eid", "date", "sample_date", "measurement_date", "instance"} or name in (exclude or set()) or name == explicitly_dated:
                continue
            if re.search(r"date|time|batch|plate|well|assay|visit|sample", name, re.I):
                continue
            parts = field_parts(name)
            if parts and parts[0] in {"p53", "p54", "p31", "p34", "p52", "p191", "p21003", "p40000"}:
                continue
            if modality == "clinical" and (not parts or parts[0] not in clinical):
                continue
            instance = parts[1] if parts else None
            feature = parts[0] + (f"_a{parts[2]}" if parts[2] not in (None, 0) else "") if parts else name.upper() if modality == "prot" else name
            identifiers = {name, feature, modality + ":" + feature}
            if parts:
                identifiers.add(parts[0])
            if include_features is not None and not identifiers.intersection(include_features):
                continue
            if identifiers.intersection(exclude_features):
                continue
            columns.append((name, instance, feature))
        if modality == "prot" and not explicitly_dated and not any(k in names for k in ("date", "sample_date", "measurement_date")) and self.d.get("prot_date") != "baseline":
            raise ValueError("Undated wide protein table: explicitly set data.prot_date: baseline only for baseline samples, or provide dated canonical observations")
        if modality == "prot" and self.d.get("prot_date") == "baseline":
            self.messages.append("Undated protein rows explicitly assigned to baseline collection date by data.prot_date=baseline")
        for n, row in enumerate(rows(path), 1):
            eid = eid_string(choose(row, ["eid", "participant.eid"]))
            if eid not in self.people:
                continue
            self.add_visits(eid, row)
            groups = collections.defaultdict(dict)
            for name, instance, feature in columns:
                value = number(row[name], nonnegative=modality in {"clinical", "bbc", "met"})
                if value is None:
                    continue
                if instance is None:
                    inst = int(number(row.get("instance")) or 0)
                else:
                    inst = instance
                explicit = day(row.get(explicitly_dated)) if explicitly_dated else day(choose(row, ["date", "sample_date", "measurement_date"]))
                if explicit is not None:
                    obs_day = explicit
                elif modality == "prot" and instance is None:
                    obs_day = self.people[eid][1] if self.d.get("prot_date") == "baseline" else None
                else:
                    obs_day = day(field_value(row, "p53", inst)) or self.visits.get((eid, inst))
                if obs_day is not None:
                    groups[obs_day][feature] = value
                else:
                    self.count[f"{modality}_values_without_visit_date"] += 1
            for obs_day, values in groups.items():
                self.add_observation(eid, obs_day, modality, values)
            if n % 20000 == 0:
                self.db.commit()
        self.db.commit()


def _moments_update(stats, name, value):
    n, mean, m2 = stats.get(name, (0, 0.0, 0.0))
    n += 1
    delta = value - mean
    mean += delta / n
    m2 += delta * (value - mean)
    stats[name] = (n, mean, m2)


def prepare(config: dict) -> dict:
    """Build cohort.sqlite and its train-only schema, replacing atomically at end."""
    root, out = _paths(config)
    data_dir = out / "data"; data_dir.mkdir(parents=True, exist_ok=True)
    temporary = data_dir / "cohort.sqlite.building"
    if temporary.exists():
        temporary.unlink()
    db = sqlite3.connect(temporary)
    try:
        db.execute("PRAGMA journal_mode=OFF"); db.execute("PRAGMA synchronous=OFF"); db.execute("PRAGMA temp_store=FILE")
        db.execute("PRAGMA cache_size=-65536")
        db.executescript(SCHEMA)
        imp = Importer(db, config, root)
        _, vip = imp.participants()
        imp.events()
        imp.observations(vip)
        data_cfg = config.get("data", {})
        count_min = int(data_cfg.get("min_token_count", 20))
        token_counts = dict(db.execute("SELECT e.code,count(*) FROM events e JOIN participants p ON p.eid=e.eid WHERE p.split='train' GROUP BY e.code"))
        vocab = ["<PAD>", "<UNK>", "<START>", "DEATH"] + sorted(k for k, count in token_counts.items() if k != "DEATH" and count >= count_min)
        feature_stats, context_stats = {}, {}
        feature_people, modality_people = collections.Counter(), collections.Counter()
        last_eid, seen_features, seen_modalities = None, set(), set()
        for eid, modality, payload in db.execute("SELECT o.eid,o.modality,o.values_json FROM observations o JOIN participants p ON p.eid=o.eid WHERE p.split='train' ORDER BY o.eid"):
            if eid != last_eid:
                feature_people.update(seen_features); modality_people.update(seen_modalities)
                seen_features, seen_modalities, last_eid = set(), set(), eid
            seen_modalities.add(modality)
            for name, value in json.loads(payload).items():
                full = modality + ":" + name
                _moments_update(feature_stats, full, float(value))
                seen_features.add(full)
        feature_people.update(seen_features); modality_people.update(seen_modalities)
        for (payload,) in db.execute("SELECT context FROM participants WHERE split='train'"):
            for name, value in json.loads(payload).items():
                _moments_update(context_stats, name, float(value))
        min_feature = int(data_cfg.get("min_feature_count", 20))
        max_missing = float(data_cfg.get("max_missing_fraction", 1.0))
        if not 0 <= max_missing <= 1:
            raise ValueError("data.max_missing_fraction must lie between 0 and 1")
        # Missingness is measured among training people with that modality.
        # A sparse modality (e.g. imaging) is not discarded simply because most
        # cohort members were never assayed/imaged.
        features = sorted(k for k, (n, _, _) in feature_stats.items() if n >= min_feature and feature_people[k] / modality_people[k.split(":", 1)[0]] >= 1 - max_missing)
        context_features = sorted(context_stats, key=lambda x: (x != "sex", x))
        modalities = {}
        for feature in features:
            modalities.setdefault(feature.split(":", 1)[0], []).append(feature)
        def pack(keys, stats):
            return {k: {"n": stats[k][0], "mean": stats[k][1], "std": max(math.sqrt(max(0.0, stats[k][2]) / max(1, stats[k][0] - 1)), 1e-6)} for k in keys}
        endpoints = data_cfg.get("endpoints", DEFAULT_ENDPOINTS)
        if endpoints == "all":
            endpoints = {code: [code] for code in vocab if code3(code)}
        elif isinstance(endpoints, list):
            endpoints = {str(code): [str(code)] for code in endpoints}
        elif not isinstance(endpoints, dict):
            raise ValueError("data.endpoints must be a name-to-code-list mapping, a code list, or 'all'")
        endpoints = {str(name): [str(c).upper().replace(".", "") for c in (codes if isinstance(codes, list) else [codes])] for name, codes in endpoints.items()}
        if not endpoints:
            raise ValueError("No outcome endpoints selected")
        for name, codes in endpoints.items():
            if any(code3(code) is None for code in codes):
                raise ValueError(f"Endpoint {name}: use explicit ICD-10 codes/prefixes or DEATH")
        splits = dict(db.execute("SELECT split,count(*) FROM participants GROUP BY split"))
        if not all(splits.get(s, 0) for s in ("train", "validation", "test")):
            raise ValueError("A data split is empty; add participants or adjust the deterministic split seed")
        meta = {"schema_version": 1, "input_root": str(root.resolve()), "vocab": vocab, "token_counts_train": token_counts, "features": features, "feature_modalities": [f.split(":", 1)[0] for f in features], "modalities": modalities, "context_features": context_features, "numeric_scaler": pack(features, feature_stats), "context_scaler": pack(context_features, context_stats), "feature_people_train": dict(feature_people), "modality_people_train": dict(modality_people), "endpoints": endpoints, "endpoint_names": list(endpoints), "split_counts": splits, "summary": dict(imp.count), "sources": imp.source_info, "event_source_mode": imp.event_mode, "warnings": imp.messages, "admin_end": data_cfg.get("admin_end"), "followup_provenance": imp.followup_provenance, "seed": imp.seed, "data_config": data_cfg, "preprocessing_fit_split": "train", "same_day_order": "numeric observation node first, then lexical ICD-10; targets strictly after that date", "normalization_clip": float(data_cfg.get("normalization_clip", 10.0))}
        db.execute("CREATE INDEX observations_eid ON observations(eid,day)")
        db.execute("CREATE INDEX participants_split ON participants(split,eid)")
        db.commit(); db.execute("PRAGMA optimize"); db.close()
        # No partial database is ever published on an import error.
        os.replace(temporary, data_dir / "cohort.sqlite")
        meta_path = data_dir / "meta.json"
        meta_tmp = data_dir / "meta.json.tmp"
        meta_tmp.write_text(json.dumps(meta, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
        os.replace(meta_tmp, meta_path)
        print(f"[flux] prepared {sum(splits.values()):,} people, {imp.count['events']:,} events, {len(features):,} numeric features; source={imp.event_mode}", flush=True)
        return meta
    except Exception:
        db.close()
        temporary.unlink(missing_ok=True)
        raise


def load_meta(path_or_config) -> dict:
    if isinstance(path_or_config, dict):
        _, out = _paths(path_or_config)
        path = out / "data" / "meta.json"
    else:
        path = Path(path_or_config)
        if path.suffix == ".sqlite":
            path = path.parent / "meta.json"
        elif path.is_dir():
            path = path / "meta.json" if (path / "meta.json").exists() else path / "data" / "meta.json"
    return json.loads(path.read_text())


class CohortDataset:
    """One participant per item, with no validation/test fitting or future input.

    Pretraining samples one contiguous window. Risk inference retains the entire
    prefix. Targets are the first code on the next STRICTLY later clinical day;
    lexical code ordering resolves simultaneous records without inventing times.
    """
    def __init__(self, config_or_db, split="train", meta=None, mode="pretrain", max_nodes=96, landmark=None, seed=None, endpoints=None, horizon_years=5.0):
        if mode not in {"pretrain", "risk"}:
            raise ValueError("mode must be pretrain or risk")
        if isinstance(config_or_db, dict):
            self.config = config_or_db
            _, out = _paths(config_or_db)
            self.db_path = str(out / "data" / "cohort.sqlite")
        else:
            self.db_path = str(config_or_db)
            self.config = {}
        self.meta = meta or load_meta(self.db_path)
        self.d = {**self.meta.get("data_config", {}), **self.config.get("data", {})}
        self.mode, self.split, self.max_nodes = mode, split, int(max_nodes)
        self.seed = int(seed if seed is not None else self.meta.get("seed", 2026))
        self.epoch = 0
        self.landmark = landmark or self.d.get("landmark", "baseline")
        self.horizon = float(horizon_years)
        self.endpoints = endpoints or self.meta["endpoints"]
        self.endpoint_names = list(self.endpoints)
        self.features = self.meta["features"]
        self.feature_index = {f: i for i, f in enumerate(self.features)}
        self.token_index = {k: i for i, k in enumerate(self.meta["vocab"])}
        self._connection, self._pid = None, None
        conn = sqlite3.connect(self.db_path)
        query = "SELECT * FROM participants" + (" WHERE split=?" if split not in {None, "all"} else "") + " ORDER BY eid"
        values = conn.execute(query, (split,) if split not in {None, "all"} else ())
        self.people = []
        for p in values:
            eid, birth, baseline, end, death, part, context = p
            anchor = self._landmark_for(conn, eid, birth, baseline) if mode == "risk" else baseline
            if mode == "risk" and (anchor is None or anchor < birth or anchor >= end or (death is not None and death <= anchor)):
                continue
            self.people.append({"eid": eid, "birth_day": birth, "baseline_day": baseline, "exit_day": end, "death_day": death, "split": part, "context": json.loads(context), "landmark_day": anchor})
        conn.close()
        self.modalities = self.meta["modalities"]

    def _landmark_for(self, conn, eid, birth, baseline):
        required = int(self.d.get("min_clinical_days", 3 if self.landmark == "midpoint" else 0))
        clinical = None
        if self.landmark == "midpoint" or required > 0:
            clinical = [row[0] for row in conn.execute("SELECT DISTINCT day FROM events WHERE eid=? AND code!='DEATH' ORDER BY day", (eid,))]
            if len(clinical) < required:
                return None
        if self.landmark == "midpoint":
            if len(clinical) < 2:
                return None
            return clinical[min(len(clinical) // 2, len(clinical) - 2)]
        if self.landmark == "baseline":
            return baseline
        if self.landmark == "age":
            return birth + int(round(float(self.d.get("landmark_age", 60.0)) * 365.25))
        if self.landmark == "date":
            anchor = day(self.d.get("landmark_date"))
            if anchor is None:
                raise ValueError("data.landmark=date requires landmark_date")
            return anchor
        if self.landmark == "last_event":
            # Fixed reference baseline prevents target-driven choice of a later anchor.
            hit = conn.execute("SELECT max(day) FROM events WHERE eid=? AND day<=?", (eid, baseline)).fetchone()
            return hit[0] if hit and hit[0] is not None else baseline
        raise ValueError("data.landmark must be midpoint, baseline, age, date or last_event (on/before baseline)")

    def _db(self):
        pid = os.getpid()
        if self._connection is None or self._pid != pid:
            if self._connection is not None:
                self._connection.close()
            self._connection = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True)
            self._pid = pid
        return self._connection

    def __getstate__(self):
        state = dict(self.__dict__); state["_connection"] = None; state["_pid"] = None; return state

    def __len__(self):
        return len(self.people)

    def set_epoch(self, epoch):
        self.epoch = int(epoch)

    def _context(self, values):
        features = self.meta["context_features"]
        x, mask = np.zeros(len(features), np.float32), np.zeros(len(features), bool)
        clip = self.meta.get("normalization_clip", 10.0)
        for i, f in enumerate(features):
            if f in values:
                stat = self.meta["context_scaler"][f]
                x[i] = np.clip((values[f] - stat["mean"]) / stat["std"], -clip, clip); mask[i] = True
        return x, mask

    def __getitem__(self, index):
        person = self.people[index]
        eid, birth = person["eid"], person["birth_day"]
        events = self._db().execute("SELECT day,code FROM events WHERE eid=? ORDER BY day,code", (eid,)).fetchall()
        limit = person["landmark_day"] if self.mode == "risk" else person["exit_day"]
        numeric_rows = self._db().execute("SELECT day,modality,values_json FROM observations WHERE eid=? AND day<=? ORDER BY day,modality", (eid, limit)).fetchall()
        daily = collections.defaultdict(dict)
        for obs_day, mod, payload in numeric_rows:
            for feature, value in json.loads(payload).items():
                full = mod + ":" + feature
                if full in self.feature_index:
                    daily[obs_day][full] = value
        nodes = [(d, 1, code, {}) for d, code in events if d <= limit]
        nodes += [(d, 0, None, values) for d, values in daily.items()]
        nodes.sort(key=lambda x: (x[0], x[1], x[2] or ""))
        start = max(birth, min([n[0] for n in nodes] + [limit]) - 1)
        nodes.insert(0, (start, -1, "<START>", {}))
        clinical_days = [d for d, _ in events]
        if self.mode == "pretrain" and len(nodes) > self.max_nodes:
            if self.split == "train":
                hashed = int.from_bytes(hashlib.blake2b(f"{self.seed}:{self.epoch}:{eid}".encode(), digest_size=8).digest(), "little")
                offset = int(np.random.default_rng(hashed).integers(0, len(nodes) - self.max_nodes + 1))
            else:
                offset = max(0, len(nodes) - self.max_nodes)
            nodes = nodes[offset:offset + self.max_nodes]
        size, width = len(nodes), len(self.features)
        numeric = np.zeros((size, width), np.float32); mask = np.zeros((size, width), bool)
        tokens = np.zeros(size, np.int64); times = np.zeros(size, np.float32)
        target = np.zeros(size, np.int64); gap = np.zeros(size, np.float32); target_mask = np.zeros(size, bool)
        clip = self.meta.get("normalization_clip", 10.0)
        for i, (when, _, code, values) in enumerate(nodes):
            times[i] = (when - birth) / 365.25
            tokens[i] = self.token_index.get(code, 1) if code else 0
            for feature, value in values.items():
                j = self.feature_index[feature]; stat = self.meta["numeric_scaler"][feature]
                numeric[i, j] = np.clip((value - stat["mean"]) / stat["std"], -clip, clip); mask[i, j] = True
            next_index = bisect.bisect_right(clinical_days, when)
            if self.mode == "pretrain" and next_index < len(events) and code != "DEATH":
                next_day, next_code = events[next_index]
                target[i] = self.token_index.get(next_code, 1)
                gap[i] = (next_day - when) / 365.25
                target_mask[i] = True
        context, context_mask = self._context(person["context"])
        item = {"eid": eid, "times": times, "valid": np.ones(size, bool), "tokens": tokens, "numeric_flat": numeric, "numeric_mask_flat": mask, "context": context, "context_mask": context_mask, "next_event": target, "next_dt": gap, "label_mask": target_mask, "feature_modalities": self.meta["feature_modalities"], "split": person["split"], "query_age": np.float32((limit - birth) / 365.25)}
        if self.mode == "risk":
            item.update(self._risk(person, events))
        return item

    def _risk(self, p, events):
        landmark = p["landmark_day"]
        first = {}
        for when, code in events:
            first.setdefault(code, when)
        d = len(self.endpoint_names)
        labels = np.zeros(d, np.float32); risk_mask = np.ones(d, bool)
        event_times = np.full(d, np.inf, np.float32)
        time = np.zeros(d, np.float32); status = np.zeros(d, bool)
        administrative = max(0, (p["exit_day"] - landmark) / 365.25)
        death_time = np.inf if p["death_day"] is None else (p["death_day"] - landmark) / 365.25
        for i, name in enumerate(self.endpoint_names):
            codes = self.endpoints[name]
            hits = [when for code, when in first.items() if any(code.startswith(prefix[:3]) if prefix != "DEATH" else code == "DEATH" for prefix in codes)]
            first_day = min(hits) if hits else None
            if first_day is not None:
                event_times[i] = (first_day - landmark) / 365.25
                risk_mask[i] = first_day > landmark
                labels[i] = float(first_day > landmark and event_times[i] <= self.horizon)
            is_death = codes == ["DEATH"]
            censor = administrative if is_death else min(administrative, death_time)
            time[i] = max(0.0, min(censor, float(event_times[i])))
            status[i] = bool(risk_mask[i] and np.isfinite(event_times[i]) and event_times[i] <= censor)
        return {"landmark_day": landmark, "labels": labels, "risk_mask": risk_mask, "event_times": event_times, "event": np.isfinite(event_times) & risk_mask, "time": time, "status": status, "exit_time": np.float32(administrative), "death_time": np.float32(death_time)}


def collate_batch(items, meta=None):
    """Pad sequences and split numeric features into model modality encoders."""
    import torch
    if not items:
        raise ValueError("Cannot collate an empty batch")
    batch, length = len(items), max(len(item["times"]) for item in items)
    modalities = meta["feature_modalities"] if meta is not None else items[0]["feature_modalities"]
    width = len(modalities)
    arrays = {"times": np.zeros((batch, length), np.float32), "valid": np.zeros((batch, length), bool), "tokens": np.zeros((batch, length), np.int64), "next_event": np.zeros((batch, length), np.int64), "next_dt": np.zeros((batch, length), np.float32), "label_mask": np.zeros((batch, length), bool)}
    values = np.zeros((batch, length, width), np.float32); masks = np.zeros((batch, length, width), bool)
    for i, item in enumerate(items):
        n = len(item["times"])
        for k in arrays:
            arrays[k][i, :n] = item[k]
        if n < length:
            arrays["times"][i, n:] = item["times"][-1]
        values[i, :n] = item["numeric_flat"]; masks[i, :n] = item["numeric_mask_flat"]
    result = {k: torch.from_numpy(v) for k, v in arrays.items()}
    result["numeric"] = {}
    for mod in dict.fromkeys(modalities):
        indices = [i for i, m in enumerate(modalities) if m == mod]
        result["numeric"][mod] = {"values": torch.from_numpy(values[:, :, indices]), "mask": torch.from_numpy(masks[:, :, indices])}
    for k in ("context", "context_mask", "query_age", "labels", "risk_mask", "event_times", "event", "time", "status", "exit_time", "death_time", "landmark_day"):
        if k in items[0]:
            result[k] = torch.as_tensor(np.stack([item[k] for item in items]))
    result["eid"] = [x["eid"] for x in items]
    result["split"] = [x["split"] for x in items]
    return result
