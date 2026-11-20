#!/usr/bin/env python3
"""Independent-method alignment evidence for a hierarchical stage's QC sample.

Only the saved stratified sample is evaluated.  Alignments reuse the available
reference collection; they are not independent reference truth.  Each pair is
tested against one explicitly bounded competitive panel.  MAPQ is never compared
between panels, and missing organisms remain possible.  This module does not
change assignments or move reads between bins.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
import math
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import time
import traceback

from common import atomic_json, fingerprint, sha256, code_identity
from alignment_cache import panel_index, checkpoint_read, checkpoint_write
from fastq_stream import fragments, new_digest, read_id, update_digest
from plant_validate import build_panel, evidence_row, file_state, plan_panels, parse_sam
from reference_audit import audit, choose_species_representatives, enrich_assembly_metadata
from species_refine import _load_representatives
from storage_layout import library_root
from taxonomy import Taxonomy


METHOD = "minimap2_competitive_paired_alignment"
REQUIRED = {"ordinal", "read_id", "group", "taxid", "species_taxid", "status",
            "bin", "source_group", "stratum", "stratum_pairs", "sampled_pairs"}
INTERPRETATION = (
    "Stratified-sample, bounded-panel compatibility evidence; accuracy is not estimated. "
    "A different alignment method reuses available reference genomes, so shared reference "
    "contamination is not excluded. Each pair is evaluated in one panel; MAPQ is within "
    "that panel only. Unknown organisms, omitted relatives and references in other panels "
    "are not excluded. Group-stage support does not validate species. Off-target human "
    "compatibility is a possible missed signal, not a proven misclassification. Weighted "
    "category totals extrapolate sampled evidence only and do not establish biological truth. "
    "Contig classes are header heuristics, not interval-level organelle/repeat annotations."
)


class Unavailable(RuntimeError):
    """An optional evaluation prerequisite is absent."""


def _table(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _group(value):
    return {"plants": "plant", "viruses": "virus", "archaeal": "archaea"}.get(value, value)


def _canonical_id(value):
    if value.endswith(("/1", "/2")):
        return value[:-2]
    matched = re.fullmatch(r"([SED]RR\d+\.\d+)\.[12]", value)
    return matched[1] if matched else value


def _options(cfg):
    options = dict(cfg.get("stage_alignment") or {})
    options.setdefault("enabled", True)
    options.setdefault("max_references", 60)
    options.setdefault("max_reference_gib", 24)
    options.setdefault("offtarget_per_stratum", 100)
    options.setdefault("seed", 20261115)
    options.setdefault("threads", min(int(cfg.get("classify_threads", 24)), 24))
    options.setdefault("min_mapq", 20)
    options.setdefault("min_identity", .95)
    options.setdefault("min_aligned_fraction", .8)
    options.setdefault("keep_alignment_files", False)
    for name in ("max_references", "threads"):
        if isinstance(options[name], bool) or int(options[name]) != options[name] or int(options[name]) < 1:
            raise ValueError(name + " must be a positive integer")
        options[name] = int(options[name])
    for name in ("offtarget_per_stratum", "seed", "min_mapq"):
        if isinstance(options[name], bool) or int(options[name]) != options[name] or int(options[name]) < 0:
            raise ValueError(name + " must be a nonnegative integer")
        options[name] = int(options[name])
    if options["min_mapq"] >= 255:
        raise ValueError("min_mapq must be below 255")
    for name in ("min_identity", "min_aligned_fraction"):
        options[name] = float(options[name])
        if not math.isfinite(options[name]) or not 0 < options[name] <= 1:
            raise ValueError(name + " must be in (0, 1]")
    options["max_reference_gib"] = float(options["max_reference_gib"])
    if not math.isfinite(options["max_reference_gib"]) or options["max_reference_gib"] <= 0:
        raise ValueError("max_reference_gib must be positive")
    return options


def _read_sample(qc):
    paths = [qc / "sample.tsv", qc / "sample.R1.fastq.gz", qc / "sample.R2.fastq.gz"]
    if any(not path.is_file() for path in paths):
        raise Unavailable("Missing stratified sample table or paired sample FASTQs")
    states = [file_state(path) for path in paths]
    rows, seen, strata = [], set(), {}
    with paths[0].open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not REQUIRED <= set(reader.fieldnames):
            raise RuntimeError("sample.tsv missing required columns: " + ", ".join(sorted(REQUIRED)))
        for sample_ordinal, row in enumerate(reader, 1):
            for key in ("ordinal", "taxid", "species_taxid", "stratum_pairs", "sampled_pairs"):
                row[key] = int(row[key])
            ordinal = row["ordinal"]
            if ordinal <= 0 or ordinal in seen:
                raise RuntimeError("Sample original ordinals must be unique positive integers")
            seen.add(ordinal)
            if row["taxid"] < 0 or row["species_taxid"] < 0:
                raise RuntimeError("Negative sample taxid")
            if not row["stratum"] or not 0 < row["sampled_pairs"] <= row["stratum_pairs"]:
                raise RuntimeError("Invalid sample stratum population/sample count")
            current = (row["stratum_pairs"], row["sampled_pairs"])
            if row["stratum"] in strata and strata[row["stratum"]] != current:
                raise RuntimeError("Inconsistent N/n within sample stratum")
            strata[row["stratum"]] = current
            row.update(sample_ordinal=sample_ordinal, group=_group(row["group"]),
                       source_taxid=row["taxid"], source_species_taxid=row["species_taxid"])
            rows.append(row)
    observed = Counter(row["stratum"] for row in rows)
    if any(observed[name] != counts[1] for name, counts in strata.items()):
        raise RuntimeError("sample.tsv rows do not match declared sampled_pairs per stratum")
    digest = new_digest()
    count = 0
    for count, (left, right) in enumerate(fragments(*paths[1:]), 1):
        if count > len(rows):
            raise RuntimeError("Sample FASTQs contain more pairs than sample.tsv")
        row = rows[count - 1]
        first_id, second_id = read_id(left[0]).decode(), read_id(right[0]).decode()
        if row["read_id"] not in (first_id, _canonical_id(first_id)):
            raise RuntimeError("Sample FASTQ order/ID differs from sample.tsv at row " + str(count))
        row.update(r1_read_id=first_id, r2_read_id=second_id)
        sequences = [record[1].strip().upper() for record in (left, right)]
        row["low_complexity"] = any(
            max(Counter(seq).values()) / len(seq) > .8 or
            len({seq[i:i + 3] for i in range(max(0, len(seq) - 2))}) < 5
            for seq in sequences)
        update_digest(digest, left, right)
    if count != len(rows):
        raise RuntimeError("Sample FASTQs contain fewer pairs than sample.tsv")
    return rows, paths[1:], states, digest.hexdigest()


def _select(rows, options):
    chosen, other = [], defaultdict(list)
    for source in rows:
        row = dict(source)
        if row["group"] in ("human", "plant"):
            row.update(role="positive_candidate", target_group=row["group"],
                       evaluation_sample_probability=1.0,
                       evidence_weight=row["stratum_pairs"] / row["sampled_pairs"])
            chosen.append(row)
        else:
            other[row["stratum"]].append(row)
    rng = random.Random(options["seed"])
    for name in sorted(other):
        available = sorted(other[name], key=lambda row: row["ordinal"])
        count = min(len(available), options["offtarget_per_stratum"])
        if not count:
            continue
        for row in rng.sample(available, count):
            probability = count / len(available)
            row.update(role="offtarget_human_screen", target_group="human",
                       evaluation_sample_probability=probability,
                       evidence_weight=row["stratum_pairs"] / row["sampled_pairs"] / probability)
            chosen.append(row)
    for row in chosen:
        row["sample_probability"] = row["sampled_pairs"] / row["stratum_pairs"]
        row["joint_inclusion_probability"] = row["sample_probability"] * row["evaluation_sample_probability"]
        row["original_ordinal"] = row["ordinal"]
        row["assigned_taxid"] = (9606 if row["role"] == "offtarget_human_screen" else
                                 row["taxid"] or (9606 if row["group"] == "human" else 33090))
        if row["role"] == "offtarget_human_screen":
            row["species_taxid"] = 9606
        elif row["group"] == "human" and not row["species_taxid"]:
            row["species_taxid"] = 9606
    return sorted(chosen, key=lambda row: row["ordinal"])


def _mapper(cfg, options):
    explicit = options.get("minimap2") or cfg.get("stage_alignment_minimap2") or cfg.get("minimap2")
    candidates = [explicit] if explicit else []
    if not explicit:
        if cfg.get("tools_dir"):
            candidates.append(str(Path(cfg["tools_dir"]) / "minimap2/2.28/minimap2"))
        candidates.append("minimap2")
    for candidate in candidates:
        candidate = str(candidate)
        path = Path(candidate).expanduser() if os.sep in candidate else None
        resolved = str(path.resolve()) if path is not None else shutil.which(candidate)
        if resolved and Path(resolved).is_file() and os.access(resolved, os.X_OK):
            return resolved
    raise Unavailable("minimap2 executable unavailable: " + ", ".join(str(x) for x in candidates))


def _references(cfg, options, out, taxonomy):
    explicit = options.get("representatives") or cfg.get("stage_alignment_representatives")
    def absolute(value):
        path = Path(value).expanduser()
        return path if path.is_absolute() else Path(cfg.get("root", ".")) / path
    candidates = [absolute(explicit)] if explicit else []
    if not explicit and cfg.get("work_dir"):
        candidates.append(Path(cfg["work_dir"]) / "qc/reference-audit/species-representatives.tsv")
    if not explicit:
        for key in ("group_library", "group_library_dir"):
            if cfg.get(key):
                candidates.append(absolute(cfg[key]) / "references.tsv")
        if cfg.get("library_dir"):
            candidates.append(Path(library_root(cfg)) / "group/references.tsv")
    table_path = next((path for path in candidates if path.is_file()), None)
    raw_root = cfg.get("reference_raw", cfg.get("root", "."))
    if table_path is None:
        if explicit:
            raise Unavailable("Configured representative table not found: " + str(explicit))
        try:
            _, records, _, _ = audit(cfg)
            frozen_species = [row["species_taxid"] for row in records]
            enrich_assembly_metadata(records, cfg["metadata_dir"])
            # Quality metadata may contain a different species-taxid field;
            # keep the species derived from the frozen classification taxonomy.
            for row, species in zip(records, frozen_species):
                row["species_taxid"] = species
            representatives = choose_species_representatives(records)
        except (OSError, RuntimeError, KeyError, ValueError) as exc:
            raise Unavailable("No representative table; metadata-only reference audit unavailable: " + str(exc)) from exc
        if not representatives:
            raise Unavailable("Reference audit found no species representatives")
        table_path = out / "representatives.tsv"
        _table(table_path, representatives)
    with table_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        compact_catalog = {"relative_path", "source_taxid", "records"} <= fields
        records = list(reader) if compact_catalog else None
    if compact_catalog:
        if not {"species_taxid", "name", "group", "path", "bases"} <= fields:
            raise RuntimeError("Compact references.tsv lacks representative-selection fields")
        frozen_species = [row["species_taxid"] for row in records]
        if cfg.get("metadata_dir"):
            enrich_assembly_metadata(records, cfg["metadata_dir"])
        for row, species in zip(records, frozen_species):
            row["species_taxid"] = species
            row.setdefault("accession", "")
            row.setdefault("compressed_bytes", 0)
        representatives = choose_species_representatives(records)
        reps = {}
        for row in representatives:
            tid, bases = int(row["species_taxid"]), int(row["bases"])
            if bases <= 0:
                raise RuntimeError("Representative base count must be positive")
            path = Path(row["path"])
            path = path if path.is_absolute() else Path(raw_root) / path
            reps[tid] = {**row, "species_taxid": tid, "bases": bases,
                         "resolved_path": str(path.resolve())}
    else:
        reps = _load_representatives(table_path, raw_root, require_files=False)
    if not reps:
        raise Unavailable("Representative table is empty: " + str(table_path))
    normalized = {}
    for tid, row in reps.items():
        resolved = taxonomy.resolve(tid)
        if taxonomy.species(resolved) != resolved:
            raise RuntimeError("Representative taxid is not a species: " + str(tid))
        if resolved in normalized:
            raise RuntimeError("Multiple representative rows resolve to species " + str(resolved))
        normalized[resolved] = {**row, "species_taxid": resolved, "group": _group(row["group"])}
    return normalized, table_path


def _plan(selected, reps, taxonomy, options):
    counts = Counter()
    for row in selected:
        if row["source_species_taxid"]:
            counts[taxonomy.resolve(row["source_species_taxid"])] += row["evidence_weight"]
    extra = [9606, 9598, 9595]
    extra.extend(int(x) for x in options.get("competitor_taxids", []))
    budget = int(options["max_reference_gib"] * 2**30)
    panels, assignments, missing = [], {}, {}
    for group in ("human", "plant"):
        claims = [row for row in selected if row["target_group"] == group]
        if not claims:
            continue
        for row in claims:
            row["assigned_taxid"] = taxonomy.resolve(row["assigned_taxid"])
            if row["species_taxid"]:
                row["species_taxid"] = taxonomy.resolve(row["species_taxid"])
        planned, mapping, untested = plan_panels(claims, counts, reps, taxonomy, extra,
                                                options["max_references"], budget, group=group)
        for panel in planned:
            panel["panel_id"] = group + "-" + panel["panel_id"]
            panel["target_group"] = group
            panels.append(panel)
        for row in claims:
            assigned = mapping.get(row["ordinal"])
            assignments[row["ordinal"]] = group + "-" + assigned if assigned else None
            if not assigned:
                tid = row["species_taxid"] or row["assigned_taxid"]
                missing[row["ordinal"]] = untested.get(tid, "no_competitive_panel")
    return panels, assignments, missing


def _write_panel_reads(reads, selected, panel_id, assignments, out):
    subset = [row for row in selected if assignments.get(row["ordinal"]) == panel_id]
    by_sample = {row["sample_ordinal"]: row for row in subset}
    ordered = []
    paths = [out / (".panel." + panel_id + ".R1.fastq.gz"),
             out / (".panel." + panel_id + ".R2.fastq.gz")]
    with gzip.open(paths[0], "wb") as first, gzip.open(paths[1], "wb") as second:
        for index, (left, right) in enumerate(fragments(*reads), 1):
            if index not in by_sample:
                continue
            ordered.append(by_sample[index])
            for mate, (record, handle) in enumerate(((left, first), (right, second)), 1):
                handle.write(f"@{len(ordered)}/{mate}\n".encode() + record[1].rstrip(b"\r\n") +
                             b"\n+\n" + record[3].rstrip(b"\r\n") + b"\n")
    if len(ordered) != len(subset):
        raise RuntimeError("Sample pairs disappeared while preparing a panel")
    return ordered, paths


def _interpret(claim, mates, taxonomy, panel_id, stage):
    row = evidence_row(claim, mates, taxonomy, panel_id)
    original = row["category"]
    row["alignment_category"] = original
    row["claim_level"] = "group" if stage == "group" else ("species" if claim["species_taxid"] else "higher_rank")
    strong = (all(mate["passed"] for mate in mates) and not claim["low_complexity"] and
              not any(mate["near_tied"] for mate in mates) and
              mates[0]["species_taxid"] == mates[1]["species_taxid"])
    compatible_group = (taxonomy.describe(mates[0]["species_taxid"])["group"]
                        if strong else None)
    if claim["role"] == "offtarget_human_screen":
        row["claim_level"] = "offtarget_human_screen"
        if strong and row["concordant_pair"]:
            row["category"] = ("potential_missed_human_compatible" if compatible_group == "human"
                               else "other_panel_compatible")
        elif strong:
            row["category"] = "ambiguous_pairing"
        else:
            row["category"] = "unmapped" if original == "unmapped" else "ambiguous"
    elif stage == "group":
        if strong and row["concordant_pair"]:
            row["category"] = ("group_compatible" if compatible_group == claim["target_group"]
                               else "conflicting_competitor")
        elif strong:
            row["category"] = "ambiguous_pairing"
        else:
            row["category"] = "unmapped" if original == "unmapped" else "ambiguous"
    elif original in ("supported_same_species", "higher_rank_compatible", "conflicting_competitor"):
        row["category"] = ("ambiguous_pairing" if not row["concordant_pair"] else
                           "species_compatible" if original == "supported_same_species" else original)
    # Keep the public sampled label unchanged, including off-target bacteria
    # screened against human. The mapper's tested hypothesis is a separate field.
    row["evaluation_assigned_taxid"] = row["assigned_taxid"]
    row["evaluation_species_taxid"] = row["species_taxid"]
    row["species_taxid"] = row["source_species_taxid"]
    return row


def _untested(claim, reason, panel_id=""):
    return {**claim, "panel_id": panel_id, "category": "untested", "alignment_category": "untested",
            "tested": False, "untested_reason": reason, "mapq_scope": "not_tested",
            "evaluation_assigned_taxid": claim["assigned_taxid"],
            "evaluation_species_taxid": claim["species_taxid"],
            "species_taxid": claim["source_species_taxid"]}


def _reuse_completed(out, options, states, panels, selected, context):
    previous = context.get('previous', {})
    scientific = lambda value: {k:v for k,v in value.items() if k not in
        ('cache_enabled','reference_cache_gib','index_cache_gib','keep_alignment_files')}
    if (not options.get('cache_enabled', True) or previous.get('status') != 'complete'
            or scientific(previous.get('parameters', {})) != scientific(options)
            or context.get('previous_panels', {}).get('panels') != panels
            or previous.get('input_states', [])[:len(states)] != states):
        return False
    try:
        if any(file_state(row['path']) != row for row in previous['input_states']):
            return False
        names = ['sampled-evidence.tsv','reference-panel.tsv','omitted-references.tsv','panels.json']
        checksums = {name:sha256(out/name) for name in names}
        implementation = {name:sha256(Path(__file__).with_name(name)) for name in
            ('stage_alignment.py','plant_validate.py','taxonomy.py')}
        if previous.get('artifact_checksums'):
            if previous.get('cache_implementation') != implementation or previous['artifact_checksums'] != checksums:
                return False
        else:
            # Adopt only outputs from an explicitly audited previous code version.
            from cache_compatibility import compatible_request
            from alignment_cache import read_json
            request = read_json(out.parent/'summary.json').get('request', {})
            if not compatible_request(request, dict(request, code=code_identity())):
                return False
        with (out/'sampled-evidence.tsv').open(newline='') as handle:
            evidence = list(csv.DictReader(handle, delimiter='\t'))
        if sorted(int(row['ordinal']) for row in evidence) != sorted(row['ordinal'] for row in selected):
            return False
        previous['artifact_checksums'] = checksums
        previous['cache_implementation'] = implementation
        previous['cache_reused'] = True
        context['reused_complete'] = previous
        print('CACHE alignment COMPLETE HIT: ' + str(out), flush=True)
        return True
    except (OSError, ValueError, KeyError):
        return False


def _run(qc, out, stage, cfg, report, context):
    options = _options(cfg)
    report["parameters"] = options
    rows, reads, states, digest = _read_sample(qc)
    context.update(rows=rows, selected=_select(rows, options), states=states)
    report.update(input_sample_pairs=len(rows), sample_input_digest=digest,
                  selected_pairs=len(context["selected"]))
    if options["enabled"] is False:
        report.update(status="skipped", reason="alignment_disabled_by_configuration")
        return
    if not context["selected"]:
        report.update(status="skipped", reason="no_target_candidates_or_selected_offtarget_pairs")
        return
    executable = _mapper(cfg, options)
    metadata = cfg.get("taxonomy_dir") or str(Path(cfg.get("metadata_dir", ".")) / "taxonomy")
    if not (Path(metadata) / "nodes.dmp").is_file() or not (Path(metadata) / "names.dmp").is_file():
        raise Unavailable("NCBI taxonomy files unavailable: " + str(metadata))
    taxonomy = Taxonomy(metadata)
    states.extend(file_state(path) for path in sorted(Path(metadata).glob("*.dmp")))
    reps, repfile = _references(cfg, options, out, taxonomy)
    states.extend([file_state(repfile), file_state(executable)])
    report.update(reference_table=str(repfile), mapper=executable,
                  mapper_version=subprocess.check_output([executable, "--version"], text=True, timeout=30).strip())
    panels, assignments, missing = _plan(context["selected"], reps, taxonomy, options)
    report["planned_panels"] = len(panels)
    if _reuse_completed(out, options, states, panels, context['selected'], context):
        return
    context['tables_written'] = True
    selected_refs, omitted = [], []
    for panel in panels:
        for row in panel["selected"]:
            selected_refs.append({"panel_id": panel["panel_id"], "target_group": panel["target_group"], **row})
        omitted.extend({"panel_id": panel["panel_id"], **row} for row in panel["omitted"])
    for claim in context["selected"]:
        if claim["ordinal"] in missing:
            context["evidence"][claim["ordinal"]] = _untested(claim, missing[claim["ordinal"]])
            omitted.append(dict(panel_id="", species_taxid=claim["species_taxid"] or claim["assigned_taxid"],
                                reason=missing[claim["ordinal"]], bases=0))
    _table(out / "reference-panel.tsv", selected_refs,
           list(dict.fromkeys(key for row in selected_refs for key in row)) or
           ["panel_id", "target_group", "species_taxid", "resolved_path", "bases", "role"])
    unique_omitted = [dict(zip(("panel_id", "species_taxid", "reason", "bases"), values))
                      for values in sorted({(str(row["panel_id"]), int(row["species_taxid"]),
                                             str(row["reason"]), int(row["bases"])) for row in omitted})]
    _table(out / "omitted-references.tsv", unique_omitted,
           ["panel_id", "species_taxid", "reason", "bases"])
    atomic_json(out / "panels.json", dict(panels=panels, execution="serial", interpretation=INTERPRETATION))
    report["omitted_reference_entries"] = len(unique_omitted)
    if unique_omitted:
        report["warnings"].append("References or competitor groups were omitted; see omitted-references.tsv")
    budget = int(options["max_reference_gib"] * 2**30)
    checkpoint_base = dict(schema=1, stage=stage, sample_digest=digest,
        sample_rows=rows, parameters={k:v for k,v in options.items()
            if k not in ('cache_enabled', 'reference_cache_gib', 'index_cache_gib', 'keep_alignment_files')},
        inputs=list(states), implementation={name:sha256(Path(__file__).with_name(name))
            for name in ('stage_alignment.py', 'plant_validate.py', 'taxonomy.py')})
    report['completed_panels'] = 0
    report['cached_panels'] = 0
    for panel_number, panel in enumerate(panels, 1):
        panel_started=time.monotonic()
        panel_id = panel["panel_id"]
        subset, pair = _write_panel_reads(reads, context["selected"], panel_id, assignments, out)
        if not subset:
            continue
        reference, sam = out / (".panel." + panel_id + ".fa"), out / (".panel." + panel_id + ".sam")
        try:
            reference_states = [file_state(row['resolved_path']) for row in panel['selected']]
            states.extend(reference_states)
            key = fingerprint(dict(checkpoint_base, panel=panel, references=reference_states, subset=subset))
            checkpoint = out / 'checkpoints' / (panel_id + '.json')
            saved = checkpoint_read(checkpoint, key) if options.get('cache_enabled', True) else None
            print(f'ALIGNMENT {stage} {panel_id}: ' + ('REUSE' if saved else 'START') +
                  f' ({panel_number}/{len(panels)})', flush=True)
            if saved:
                for item in saved['evidence']:
                    context['evidence'][item['ordinal']] = item
                report['commands'].append(saved['command'])
                report['cached_panels'] += 1
            else:
                with panel_index(cfg, panel['selected'], executable, report['mapper_version'],
                                 options['threads'], budget, reference) as (target, contigs, actual, cache):
                    preparation_seconds=time.monotonic()-panel_started
                    _table(out / ('contigs.' + panel_id + '.tsv'), contigs.values())
                    command = [executable, '-ax', 'sr', '--secondary=yes', '-N', '20', '-I', str(actual + 1000000),
                               '-t', str(options['threads']), str(target), *(str(path) for path in pair)]
                    command_record = dict(panel_id=panel_id, command=command, reference_bases=actual,
                        fragments=len(subset), preparation_seconds=round(preparation_seconds,3), log='alignment/minimap2.' + panel_id + '.log', cache=cache)
                    report['commands'].append(command_record)
                    print(f'CACHE alignment {panel_id}: {cache}', flush=True)
                    with sam.open('w') as stdout, (out / ('minimap2.' + panel_id + '.log')).open('w') as stderr:
                        subprocess.run(command, stdout=stdout, stderr=stderr, check=True)
                    alignments = parse_sam(sam, len(subset), contigs, options['min_mapq'],
                                           options['min_identity'], options['min_aligned_fraction'])
                    evidence = [_interpret(claim, [alignments[(local,mate)] for mate in (1,2)],
                                           taxonomy, panel_id, stage) for local,claim in enumerate(subset,1)]
                    if any(file_state(x['path']) != x for x in states):
                        raise RuntimeError('Alignment inputs changed during panel execution')
                    for item in evidence:
                        context['evidence'][item['ordinal']] = item
                    command_record['total_seconds']=round(time.monotonic()-panel_started,3)
                    if options.get('cache_enabled', True):
                        checkpoint_write(checkpoint, key, evidence, command_record)
            report['completed_panels'] += 1
            report['current_panel'] = panel_id
            atomic_json(out / 'summary.json', report)
            print(f'ALIGNMENT {stage} {panel_id}: DONE ({report["completed_panels"]}/{len(panels)})', flush=True)
        except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as exc:
            if not any(item['panel_id'] == panel_id for item in report['commands']):
                report['commands'].append(dict(panel_id=panel_id, command=getattr(exc, 'cmd', []),
                    fragments=len(subset), status='preparation_failed', error=str(exc)))
            report["errors"].append(dict(panel_id=panel_id, error=str(exc), fragments=len(subset)))
            atomic_json(out / "summary.json", report)
            print(f"ALIGNMENT FAILED {stage} {panel_id}: {exc}",file=sys.stderr,flush=True)
            for claim in subset:
                context["evidence"][claim["ordinal"]] = _untested(claim, "panel_failed: " + str(exc), panel_id)
        finally:
            if not options["keep_alignment_files"]:
                for path in [reference, sam, *pair]:
                    path.unlink(missing_ok=True)
    if any(file_state(state["path"]) != state for state in states):
        raise RuntimeError("Sample, taxonomy, mapper or reference changed during alignment evaluation")
    tested = sum(row["tested"] for row in context["evidence"].values())
    report["status"] = ("failed" if report["errors"] else "complete" if tested == len(context["selected"])
                        else "partial" if tested else "unavailable")


def _finish(out, report, context):
    reason = report.get("reason") or report["status"]
    for row in context["selected"]:
        context["evidence"].setdefault(row["ordinal"], _untested(row, reason))
    evidence = [context["evidence"][key] for key in sorted(context["evidence"])]
    fields = list(dict.fromkeys(key for row in evidence for key in row)) or [
        "ordinal", "read_id", "group", "stratum", "stratum_pairs", "sampled_pairs", "category", "tested"]
    _table(out / "sampled-evidence.tsv", evidence, fields)
    report.update(tested_pairs=sum(row["tested"] for row in evidence),
                  untested_pairs=sum(not row["tested"] for row in evidence),
                  categories=dict(sorted(Counter(row["category"] for row in evidence).items())),
                  roles=dict(sorted(Counter(row["role"] for row in evidence).items())))
    by_stratum = {}
    for row in context["rows"]:
        by_stratum.setdefault(row["stratum"], dict(stratum=row["stratum"], stratum_pairs=row["stratum_pairs"],
            sampled_pairs=row["sampled_pairs"], selected_for_alignment=0, tested_pairs=0,
            categories=Counter(), weighted_category_totals=Counter(),
            weighting="N/n times inverse conditional alignment-subsample probability; evidence only"))
    for row in evidence:
        summary = by_stratum[row["stratum"]]
        summary["selected_for_alignment"] += 1
        summary["tested_pairs"] += int(row["tested"])
        summary["categories"][row["category"]] += 1
        summary["weighted_category_totals"][row["category"]] += row["evidence_weight"]
    report["by_stratum"] = [by_stratum[name] for name in sorted(by_stratum)]
    report["input_states"] = context.get("states", [])
    if report.get('status') == 'complete':
        report['cache_implementation'] = {name:sha256(Path(__file__).with_name(name)) for name in
            ('stage_alignment.py','plant_validate.py','taxonomy.py')}
        report['artifact_checksums'] = {name:sha256(out/name) for name in
            ('sampled-evidence.tsv','reference-panel.tsv','omitted-references.tsv','panels.json')}
    atomic_json(out / "summary.json", report)


def evaluate_alignment(sample_root, stage, cfg):
    """Write optional evaluation artifacts; return a JSON-serializable summary.

    Missing optional tools/references return ``unavailable``. Malformed samples,
    failed mapper panels and other execution errors return ``failed`` with their
    diagnostics; successfully evaluated panels remain visible. Paths in the
    return value are relative to the stage QC directory, so staging-directory
    renames do not invalidate report links.
    """
    stage = {"split-group": "group", "split-species": "species"}.get(stage, stage)
    if stage not in ("group", "species"):
        raise ValueError("stage must be group or species")
    qc = Path(sample_root) / "group" / ("qc" if stage == "group" else "species.qc")
    out = qc / "alignment"
    out.mkdir(parents=True, exist_ok=True)
    report = dict(schema="metag-stage-alignment-v1", stage=stage, status="running", method=METHOD,
                  scope="stratified_sample_only", sample_only=True, accuracy_estimated=False,
                  input_sample_pairs=0, selected_pairs=0, tested_pairs=0, planned_panels=0,
                  summary_path="alignment/summary.json", evidence_path="alignment/sampled-evidence.tsv",
                  reference_panel_path="alignment/reference-panel.tsv",
                  omitted_references_path="alignment/omitted-references.tsv",
                  panel_scope="one bounded competitive panel per pair; within-panel MAPQ only",
                  interpretation=INTERPRETATION, errors=[], warnings=[], commands=[])
    from alignment_cache import read_json
    context = dict(rows=[], selected=[], evidence={}, states=[],
                   previous=read_json(out/'summary.json'), previous_panels=read_json(out/'panels.json'))
    atomic_json(out / 'summary.json', report)
    try:
        _run(qc, out, stage, cfg, report, context)
    except Unavailable as exc:
        report.update(status="unavailable", reason=str(exc))
        report["warnings"].append(str(exc))
    except Exception as exc:
        report.update(status="failed", reason=str(exc))
        report["errors"].append(dict(error=str(exc), type=type(exc).__name__))
        (out / "failure.log").write_text(traceback.format_exc(), encoding="utf-8")
    if context.get('reused_complete'):
        atomic_json(out/'summary.json', context['reused_complete'])
        return context['reused_complete']
    if not context.get('tables_written'):
        _table(out/'reference-panel.tsv', [], ['panel_id','species_taxid','resolved_path','bases','role'])
        _table(out/'omitted-references.tsv', [], ['panel_id','species_taxid','reason','bases'])
    _finish(out, report, context)
    return report
