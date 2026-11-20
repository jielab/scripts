#!/usr/bin/env python3
"""Summarize exported IBDmix intervals without loading the calling environment."""

import argparse
from collections import Counter, defaultdict
from contextlib import closing
import csv
import gzip
from html import escape
from importlib.util import module_from_spec, spec_from_file_location
from itertools import groupby
import json
import math
from pathlib import Path
import re
import sqlite3
import statistics
import tempfile


ROOT = Path(__file__).resolve().parents[1]
AUTOSOMES = tuple(map(str, range(1, 23)))
CHROMOSOMES = (*AUTOSOMES, "X")
common_spec = spec_from_file_location("gu_0_common", ROOT / "f/0.common.py")
common = module_from_spec(common_spec)
common_spec.loader.exec_module(common)
LENGTHS = common.CHROM_LENGTHS
NEANDERTHAL = {"Altai", "Altai.2013", "Chagyr", "Vindija"}
DENISOVAN = set(common.SUPPORTED_ARCHAIC_REFERENCES) - NEANDERTHAL
AFR5 = {"ESN", "GWD", "LWK", "MSL", "YRI"}
CELL = {"AFR5": 17, "EUR": 51, "EAS": 55, "SAS": 55}
CORE = (45859650, 45909024)  # hg19 1-based inclusive 45,859,651–45,909,024.
FILE = re.compile(r"chr([1-9]|1[0-9]|2[0-2]|X)(?:\.segments)?\.tsv(?:\.(?:gz|bgz))?$")
TOP = 100


def arguments():
    p = argparse.ArgumentParser(prog="gu.paper.sh ibdmix", description=__doc__)
    p.add_argument("--target", required=True, help="Output prefix, e.g. 1kg or ukb")
    p.add_argument("--dir-ibdmix", required=True, type=Path, help="Directory of whole-chromosome segment exports")
    p.add_argument("--dir-out", type=Path, default=Path.cwd(), help="Report directory [current directory]")
    p.add_argument("--grch", choices=LENGTHS, default="37")
    p.add_argument("--sample-panel", type=Path, help="Actual shared tested roster, as TSV or PLINK PSAM; 1kg defaults to Phase 3")
    p.add_argument("--x-males", action="store_true", help="Use only male roster members for chrX [automatic for bundled 1KG]")
    p.add_argument("--window-bp", type=int, default=100000)
    p.add_argument("--source", default="", help="Input provenance displayed in the report")
    a = p.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", a.target):
        p.error("invalid target")
    if a.window_bp < 10000:
        p.error("--window-bp must be at least 10000")
    return a


def panel(path):
    if path is None:
        return {}
    with path.open() as handle:
        fields = handle.readline().split()
        aliases = {"sample": ("sample", "ID", "IID", "#IID", "sample_id"),
                   "pop": ("pop", "population"), "super": ("super_pop", "super_population"),
                   "sex": ("gender", "sex", "SEX")}
        idx = {key: next((fields.index(name) for name in names if name in fields), None)
               for key, names in aliases.items()}
        if idx["sample"] is None:
            raise ValueError("sample panel needs sample, ID or IID")
        people = {}
        for line in handle:
            items = line.split()
            if not items:
                continue
            data = {key: items[i] if i is not None else "" for key, i in idx.items()}
            iid = data["sample"]
            if iid in people:
                raise ValueError(f"Duplicate IID in roster: {iid}")
            groups = {"ALL", data["pop"], data["super"]} - {"", ".", "NA", "0"}
            if data["pop"] in AFR5:
                groups.add("AFR5")
            people[iid] = (tuple(sorted(groups)), data["sex"].lower())
    if not people:
        raise ValueError(f"Empty sample panel: {path}")
    return people


def inputs(folder):
    paths = {}
    for path in sorted(folder.iterdir()):
        match = FILE.fullmatch(path.name)
        if not match or not path.is_file():
            continue
        chrom = match[1]
        if chrom in paths:
            raise ValueError(f"Multiple exports for chr{chrom}: {paths[chrom]} and {path}")
        paths[chrom] = path
    if not paths:
        raise ValueError(f"No chr*.tsv[.gz] exports in {folder}")
    return paths


def calls(path, chrom, length, people):
    opener = gzip.open if path.suffix in {".gz", ".bgz"} else open
    with opener(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ID", "chrom", "anc", "start", "end", "length", "slod"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"Unexpected segment header: {path}")
        for row in reader:
            iid, ref = row["ID"], row["anc"]
            row_chrom = row["chrom"].removeprefix("chr")
            left, right = int(row["start"]), int(row["end"])
            if (not iid or (people and iid not in people)
                    or ("X" if row_chrom == "23" else row_chrom) != chrom
                    or ref not in NEANDERTHAL | DENISOVAN
                    or not 0 <= left < right <= length
                    or int(row["length"]) != right - left
                    or not math.isfinite(float(row["slod"]))):
                raise ValueError(f"Invalid interval or sample roster mismatch: {path}:{reader.line_num}")
            yield iid, left, right, ref


def union(intervals):
    left, right = -1, -1
    for lo, hi in intervals:
        if lo > right:
            if left >= 0:
                yield left, right
            left, right = lo, hi
        else:
            right = max(right, hi)
    if left >= 0:
        yield left, right


def span_bp(intervals):
    return sum(hi - lo for lo, hi in intervals)


def shared_bp(first, second):
    total = 0
    j = 0
    for lo, hi in first:
        while j < len(second) and second[j][1] <= lo:
            j += 1
        while j < len(second) and second[j][0] < hi:
            left, right = second[j]
            total += min(hi, right) - max(lo, left)
            if right > hi:
                break
            j += 1
    return total


class Windows:
    """Interval updates for physical coverage and distinct people per window."""

    def __init__(self, size, width):
        self.width = width
        n = (size + width - 1) // width
        self.edge_bp = [0] * n
        self.full = [0] * (n + 1)
        self.carriers = [0] * (n + 1)

    def add(self, intervals):
        last_counted = -1
        for lo, hi in intervals:
            first, last = lo // self.width, (hi - 1) // self.width
            if first == last:
                self.edge_bp[first] += hi - lo
            else:
                self.edge_bp[first] += (first + 1) * self.width - lo
                self.edge_bp[last] += hi - last * self.width
                self.full[first + 1] += 1
                self.full[last] -= 1
            first = max(first, last_counted + 1)
            if first <= last:
                self.carriers[first] += 1
                self.carriers[last + 1] -= 1
                last_counted = last

    def rows(self, chrom, length, n):
        full = carriers = 0
        for i, bp in enumerate(self.edge_bp):
            full += self.full[i]
            carriers += self.carriers[i]
            if carriers:
                lo, hi = i * self.width, min(length, (i + 1) * self.width)
                equivalent = (bp + full * self.width) / (hi - lo)
                yield [chrom, lo, hi, carriers, equivalent, n]


def ranking(row):
    _, lo, hi, carriers, equivalent, n = row
    return (-(equivalent / n if n else equivalent), -carriers, CHROMOSOMES.index(row[0]), lo)


def population_rows(totals, groups, span, n):
    for group, members in sorted(groups.items()):
        values = totals.get(group, {})
        if n:
            sizes = sorted(values.get(iid, 0) / 1e6 for iid in members)
            mean = statistics.mean(sizes)
            yield [group, len(members), len(values), sum(values.values()) / 1e6, mean,
                   statistics.median(sizes), sizes[int(0.1 * (len(sizes) - 1))],
                   sizes[int(0.9 * (len(sizes) - 1))], mean * 1e8 / span]
        else:
            yield [group, None, len(values), sum(values.values()) / 1e6, None, None, None, None, None]


def summarize(a):
    paths = inputs(a.dir_ibdmix)
    builtin = a.target == "1kg" and a.sample_panel is None
    roster = panel(ROOT / "config/1kg.phase3.panel.tsv" if builtin else a.sample_panel)
    x_males = builtin or a.x_males
    if x_males and roster and "X" in paths and any(
        sex not in {"male", "female", "m", "f", "1", "2"} for _, sex in roster.values()
    ):
        raise ValueError("chrX male denominator requires known sex for all roster IDs")
    lengths = LENGTHS[a.grch]
    groups = {}
    for scope in ("Autosomes", "chrX"):
        groups[scope] = defaultdict(list)
        for iid, (labels, sex) in roster.items():
            if scope == "chrX" and x_males and sex not in {"male", "m", "1"}:
                continue
            for label in labels:
                groups[scope][label].append(iid)
        if not roster:
            groups[scope]["ALL"] = []
    totals = defaultdict(lambda: defaultdict(Counter))
    peaks = defaultdict(list)
    control = defaultdict(lambda: [0, 0])
    audit = []
    ref_chroms = defaultdict(list)
    references = set()
    unique_bp = afr_bp = shared = 0

    for chrom in CHROMOSOMES:
        if chrom not in paths:
            continue
        scope = "chrX" if chrom == "X" else "Autosomes"
        people = roster
        if chrom == "X" and x_males and roster:
            people = {iid: roster[iid] for iid in groups[scope]["ALL"]}
            if not people:
                raise ValueError("No male samples in chrX roster")
        windows = {}
        ref_counts = Counter()
        ref_bp = Counter()
        # Disk-backed sort keeps peak memory independent of the number of calls.
        # Each chromosome has a private working database, removed on exit.
        with tempfile.TemporaryDirectory(prefix="gu-paper-") as tmp:
            with closing(sqlite3.connect(Path(tmp) / "chr.sqlite")) as db:
                db.execute("PRAGMA journal_mode=OFF")
                db.execute("PRAGMA synchronous=OFF")
                db.execute("PRAGMA temp_store=FILE")
                db.execute("CREATE TABLE calls (iid TEXT, lo INTEGER, hi INTEGER, ref TEXT)")
                db.executemany("INSERT INTO calls VALUES (?,?,?,?)", calls(paths[chrom], chrom, lengths[chrom], people))
                ordered = db.execute("SELECT iid,lo,hi,ref FROM calls ORDER BY iid,lo,hi,ref")
                for iid, records in groupby(ordered, key=lambda row: row[0]):
                    tracks = defaultdict(list)
                    labels = people[iid][0] if people else ("ALL",)
                    for _, lo, hi, ref in records:
                        ref_counts[ref] += 1
                        ref_bp[ref] += hi - lo
                        tracks[ref].append((lo, hi))
                        lineage = "Neanderthal" if ref in NEANDERTHAL else "Denisovan"
                        tracks[lineage].append((lo, hi))
                    for track, intervals in tracks.items():
                        merged = list(union(intervals))
                        bp = span_bp(merged)
                        core_bp = 0
                        if chrom == "3" and a.grch == "37" and track not in DENISOVAN | {"Denisovan"}:
                            core_bp = sum(max(0, min(hi, CORE[1]) - max(lo, CORE[0])) for lo, hi in merged)
                        for group in labels:
                            totals[scope, track][group][iid] += bp
                            key = track, group
                            if key not in windows:
                                windows[key] = Windows(lengths[chrom], a.window_bp)
                            windows[key].add(merged)
                            if core_bp:
                                control[track, group][0] += 1
                                control[track, group][1] += 5 * core_bp >= 4 * (CORE[1] - CORE[0])
                if scope == "Autosomes" and a.target == "1kg" and roster and ref_counts["Altai.2013"]:
                    intervals = [[], [], []]
                    for iid, lo, hi in db.execute(
                        "SELECT iid,lo,hi FROM calls WHERE ref='Altai.2013' ORDER BY lo,hi"
                    ):
                        labels = roster[iid][0]
                        intervals[0].append((lo, hi))
                        if "AFR5" in labels:
                            intervals[1].append((lo, hi))
                        elif any(group in labels for group in CELL if group != "AFR5") or "AMR" in labels:
                            intervals[2].append((lo, hi))
                    all_people, african, non_african = map(lambda arr: list(union(arr)), intervals)
                    unique_bp += span_bp(all_people)
                    afr_bp += span_bp(african)
                    shared += shared_bp(african, non_african)
        for ref, count in sorted(ref_counts.items()):
            references.add(ref)
            ref_chroms[ref].append(chrom)
            audit.append([chrom, paths[chrom].name, ref, count, round(ref_bp[ref] / count / 1000, 2)])
        if not ref_counts:
            audit.append([chrom, paths[chrom].name, "—", 0, None])
        for (track, group), window in windows.items():
            n = len(groups[scope][group])
            key = scope, track, group
            candidates = peaks[key] + list(window.rows(chrom, lengths[chrom], n))
            peaks[key] = sorted(candidates, key=ranking)[:TOP]
        print(f"[GU PAPER] chr{chrom}: {sum(ref_counts.values()):,} exported rows", flush=True)

    summaries, coverage, top = {}, {}, {}
    for scope in ("Autosomes", "chrX"):
        expected = AUTOSOMES if scope == "Autosomes" else ("X",)
        present = [chrom for chrom in expected if chrom in paths]
        for track in sorted(references | {"Neanderthal", "Denisovan"}):
            refs = ([track] if track in references else
                    sorted(references & (NEANDERTHAL if track == "Neanderthal" else DENISOVAN)))
            track_present = [chrom for chrom in present if any(chrom in ref_chroms[r] for r in refs)]
            if not track_present:
                continue
            key = f"{scope}|{track}"
            span = sum(lengths[chrom] for chrom in track_present)
            coverage[key] = [track_present, span, refs,
                             [chrom for chrom in track_present if any(chrom not in ref_chroms[r] for r in refs)]]
            summaries[key] = list(population_rows(totals[scope, track], groups[scope], span, bool(roster)))
            top[key] = {group: peaks[scope, track, group] for group in groups[scope]}
    tests = []
    if "3" in paths and a.grch == "37":
        for track in ["Neanderthal", *sorted(references & NEANDERTHAL)]:
            used = {track} if track != "Neanderthal" else references & NEANDERTHAL
            if not any("3" in ref_chroms[r] for r in used):
                continue
            for group, members in sorted(groups["Autosomes"].items()):
                near = [(i, row) for i, row in enumerate(peaks["Autosomes", track, group], 1)
                        if row[0] == "3" and row[1] < CORE[1] and row[2] > CORE[0]]
                any_bp, core80 = control[track, group]
                tests.append([track, group, len(members) or None, any_bp, core80,
                              min((i for i, _ in near), default=None)])
    benchmark = []
    benchmark_total = []
    if a.target == "1kg" and a.grch == "37" and "Autosomes|Altai.2013" in summaries:
        key = "Autosomes|Altai.2013"
        if coverage[key][0] == list(AUTOSOMES) and roster:
            benchmark = [[group, n, mean, CELL[group], mean - CELL[group], mean / CELL[group]]
                         for group, n, _, _, mean, *_ in summaries[key] if group in CELL]
            total = next(row[3] for row in summaries[key] if row[0] == "ALL")
            benchmark_total = [total / 1000, unique_bp / 1e9, 100 * shared / afr_bp if afr_bp else None]
    return dict(target=a.target, build=a.grch, window=a.window_bp, panels=bool(roster),
                roster_source="1KG Phase 3 public panel (assumed tested roster)" if builtin else
                (f"{a.sample_panel.name} (provided tested roster)" if roster else "No tested roster supplied"),
                nx=(len(groups["chrX"].get("ALL", [])) or None) if "X" in paths else None,
                na=len(roster) or None, x_males=x_males, chromosomes=list(ch for ch in CHROMOSOMES if ch in paths),
                references=sorted(references), source=a.source, audit=audit,
                coverage=coverage, summary=summaries, top=top, control=tests, benchmark=benchmark, benchmark_total=benchmark_total)


def main():
    args = arguments()
    try:
        payload = summarize(args)
        template = (ROOT / "f/ibdmix_paper.html").read_text(encoding="utf-8")
        encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c")
        rendered = template.replace("__TITLE__", escape(args.target)).replace("__DATA__", encoded)
        args.dir_out.mkdir(parents=True, exist_ok=True)
        output = args.dir_out / f"{args.target}.segments.html"
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=args.dir_out, delete=False) as tmp:
            temporary = Path(tmp.name)
            try:
                tmp.write(rendered)
                tmp.close()
                temporary.chmod(0o644)
                temporary.replace(output)
            finally:
                temporary.unlink(missing_ok=True)
        print(output)
    except (OSError, ValueError, sqlite3.Error) as exc:
        raise SystemExit(f"GU PAPER: {exc}") from exc


if __name__ == "__main__":
    main()
