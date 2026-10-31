#!/usr/bin/env python3
"""Modern scanning, annotations, population analysis and summaries."""
from __future__ import annotations

# The shared utility filename is sorted first; use a stable import name in Python.
import importlib.util as _common_import
import sys as _common_sys
from pathlib import Path as _CommonPath
if "allele_common" not in _common_sys.modules:
    _common_spec = _common_import.spec_from_file_location(
        "allele_common", _CommonPath(__file__).with_name("0.common.py"))
    _common_module = _common_import.module_from_spec(_common_spec)
    _common_sys.modules[_common_spec.name] = _common_module
    try:
        _common_spec.loader.exec_module(_common_module)
    except BaseException:
        _common_sys.modules.pop(_common_spec.name, None)
        raise


# ========================================================================
# scan
# ========================================================================

# Modern gnomAD scanning, spectra, summaries and command handling.





# ========================================================================
# Scan
# ========================================================================

# One-pass summary-first scan of split gnomAD sites VCFs.
# 
# Only genuinely observed, PASS, counted SNV alternatives enter occupancy.
# Reference counts, a reference frequency, joint entropy and allele age are NOT
# inferred. Failed/missing alternatives remain distinct from unlisted candidates.


from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import itertools
import math
from pathlib import Path
import re
import time

from allele_common import VERSION, dump_json, fmt, tsv_writer
from allele_common import FastaContext, allele_info_values, canonical_chrom, chrom_key, iterate_records, parse_region, read_header, validate_count_header, _file_chrom
from allele_common import IntervalMask

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
ANCESTRIES = ("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas", "oth")
STRATA_FIELDS = ["chrom", "start0", "end0", "threshold", "context", "primary_change",
                 "primary_af_bin", "an_bin", "n_sites", "n_multi", "n_three"]
CANDIDATE_FIELDS = ["chrom", "pos", "ref", "alt", "status", "is_primary",
                    "primary_tie", "n_alt", "ac", "an", "af", "context",
                    "primary_alt", "transition", "cpg_transition"]
SITE_FIELDS = ["chrom", "pos", "ref", "context", "canonical_context", "cpg_site", "n_alt",
               "alts_by_af", "alt_ac", "alt_an", "alt_af", "primary_alt", "secondary_alt",
               "primary_tie", "an_min", "an_max", "an_relative_spread", "has_failed_alt",
               "has_missing_alt", "has_non_snv", "region_flags", "primary_change",
               "secondary_ref_relation", "secondary_primary_relation",
               "population_class", "n_pop_evaluable", "pops_with_multi_alt",
               "pop_presence_jaccard", "pop_profile_cosine", "pop_profile_n",
               "same_consequence_terms", "ref_frequency_status"]
ALLELE_FIELDS = ["chrom", "pos", "ref", "alt", "status", "alt_frequency_rank",
                 "is_primary", "ac", "an", "af", "original_alt_indices",
                 "consequence_terms", "gene_symbols", "region_flags"]
SITE_AC_SPECTRUM_FIELDS = ["chrom", "n_alt", "n_singleton_alts", "min_alt_ac_bin"]


def canonical_context(context, ref, alt=None):
    """Pyrimidine-centered context and mutation class; no ancestral inference."""
    if ref in "AG":
        context = context.translate(COMPLEMENT)[::-1]
        ref = ref.translate(COMPLEMENT)
        alt = alt.translate(COMPLEMENT) if alt else alt
    return context, (f"{ref}>{alt}" if alt else None)


def relation(a, b):
    if not a or not b:
        return "."
    if a == b:
        return "same"
    return "transition" if a + b in {"AG", "GA", "CT", "TC"} else "transversion"


def af_bin(af):
    for edge, name in [(1e-5, "lt_1e-5"), (1e-4, "1e-5_to_1e-4"),
                       (1e-3, "1e-4_to_1e-3"), (.01, "0.001_to_0.01"),
                       (.05, "0.01_to_0.05"), (.5, "0.05_to_0.5")]:
        if af < edge:
            return name
    return "ge_0.5"


def an_bin(an):
    for edge in (100, 1000, 10000, 50000, 100000, 120000, 140000, 150000, 200000):
        if an < edge:
            return "lt_" + str(edge)
    return "ge_200000"


def ac_bin(ac):
    if ac == 1:
        return "1"
    for edge, name in [(5, "2_to_4"), (10, "5_to_9"), (100, "10_to_99"),
                       (1000, "100_to_999")]:
        if ac < edge:
            return name
    return "ge_1000"


def _valid_counts(ac, an, where):
    if ac is not None and ac < 0 or an is not None and an < 0:
        raise ValueError(f"Negative count at {where}")
    if ac is not None and an is not None and ac > an:
        raise ValueError(f"AC>AN at {where}: AC={ac}, AN={an}")
    return ac is not None and an is not None and an > 0


def discover_populations(header, requested="auto"):
    if requested in ("none", "", None):
        return []
    info = header["info"]
    wanted = ANCESTRIES if requested == "auto" else tuple(requested.split(","))
    present = []
    for pop in wanted:
        if not re.fullmatch(r"[a-z][a-z0-9_]*", pop):
            raise ValueError(f"Invalid ancestry suffix: {pop}")
        if f"AC_{pop}" in info and f"AN_{pop}" in info:
            validate_count_header(header, f"AC_{pop}", f"AN_{pop}", f"AF_{pop}")
            present.append(pop)
        elif requested != "auto":
            raise ValueError(f"Requested ancestry {pop} lacks AC_{pop}/AN_{pop} header fields")
    return present


def vep_schema(header, enabled):
    if not enabled:
        return None, []
    for tag in ("vep", "CSQ"):
        meta = header["info"].get(tag)
        if not meta:
            continue
        text = meta["description"]
        match = re.search(r"Format:\s*(.+)", text, flags=re.I)
        if not match:
            raise ValueError(f"INFO/{tag} has no parseable Format description")
        names = match[1].strip().strip("\"' ").split("|")
        if not {"Allele", "Consequence"}.issubset(names):
            raise ValueError(f"INFO/{tag} needs Allele and Consequence in Format")
        return tag, names
    return None, []


def parse_vep(record, tag, fields):
    if not tag:
        return {}
    value = record["info"].get(tag, ".")
    if value in ("", ".", None):
        return {}
    ai, ci = fields.index("Allele"), fields.index("Consequence")
    si = fields.index("SYMBOL") if "SYMBOL" in fields else None
    result = {}
    transcripts = []
    for fragment in value.split(","):
        cells = fragment.split("|")
        # gnomAD v4.1 exports comma lists inside LOFTEE fields, including
        # LoF_flags (SINGLE_EXON,PHYLOCSF_WEAK) and terminal LoF_info.
        # Rejoin only within these explicitly named columns, and never past
        # the header's field count. A new full transcript cannot be a tail.
        if transcripts:
            previous = transcripts[-1]
            index = len(previous) - 1
            if (0 <= index < len(fields)
                    and fields[index] in {"LoF_filter", "LoF_flags", "LoF_info"}
                    and previous[-1]
                    and len(previous) + len(cells) - 1 <= len(fields)):
                previous[-1] += "," + cells[0]
                previous.extend(cells[1:])
                continue
        transcripts.append(cells)
    for cells in transcripts:
        if len(cells) <= max(ai, ci):
            raise ValueError(f"Malformed {tag} entry at {record['chrom']}:{record['pos']}")
        allele = cells[ai]
        # SNVs do not need indel-minimization conventions for VEP matching.
        if allele not in BASES or len(allele) != 1:
            continue
        item = result.setdefault(allele, [set(), set()])
        item[0].update(t for t in cells[ci].split("&") if t)
        if si is not None and si < len(cells) and cells[si]:
            item[1].add(cells[si])
    return result


def allele_filter_values(record, header):
    """Respect explicit AS_FilterStatus when supplied, without guessing shapes."""
    raw = record["info"].get("AS_FilterStatus")
    n = len(record["alts"])
    if raw in (None, "", "."):
        return [None] * n
    values = str(raw).split(",")
    if len(values) != n:
        values = str(raw).split("|")
        if len(values) == n + 1 and not values[0]:
            values = values[1:]
    if len(values) != n:
        raise ValueError(f"AS_FilterStatus/ALT cardinality mismatch at {record['chrom']}:{record['pos']}")
    return values


def group_alleles(records, header, populations, options, qc, vep_tag, vep_fields):
    """Decode a complete coordinate group, retaining failed alternatives."""
    alts, refs = {}, set()
    has_non_snv = False
    region_flags = set()
    for record in records:
        qc["records_scanned"] += 1
        ref = record["ref"].upper()
        if len(ref) != 1 or ref not in BASES:
            qc["non_snv_ref_records"] += 1
            has_non_snv = True
            continue
        refs.add(ref)
        if len(refs) > 1:
            raise ValueError(f"Conflicting SNV REF at {record['chrom']}:{record['pos']}: {refs}")
        acs = allele_info_values(record, "AC", header, kind="int", required=True)
        ans = allele_info_values(record, "AN", header, kind="int", required=True, allow_number_one=True)
        afs = allele_info_values(record, "AF", header)
        as_filters = allele_filter_values(record, header)
        pcounts = {}
        for pop in populations:
            pcounts[pop] = (
                allele_info_values(record, "AC_" + pop, header, kind="int"),
                allele_info_values(record, "AN_" + pop, header, kind="int", allow_number_one=True))
        annotations = parse_vep(record, vep_tag, vep_fields)
        flags = {tag for tag in ("lcr", "segdup") if record["info"].get(tag) not in (None, "", ".", "0", "False", "false")}
        region_flags.update(flags)
        pass_record = record["filter"] == "PASS" or (
            options["include_unfiltered"] and record["filter"] == ".")
        for index, alt in enumerate(record["alts"]):
            alt = alt.upper()
            if len(alt) != 1 or alt not in BASES or alt == ref:
                has_non_snv = True
                qc["non_snv_alt_entries"] += 1
                continue
            ac, an, given_af = acs[index], ans[index], afs[index]
            where = f"{record['chrom']}:{record['pos']}:{ref}>{alt}"
            valid = _valid_counts(ac, an, where)
            af = ac / an if valid else None
            if given_af is not None and not 0 <= given_af <= 1:
                raise ValueError(f"Invalid INFO/AF at {where}")
            if given_af is not None and af is not None and abs(given_af - af) > max(1e-6, af * 1e-3):
                qc["reported_af_disagrees_with_ac_an"] += 1
            allele_pass = pass_record and as_filters[index] in (None, ".", "PASS")
            if flags.intersection(options["exclude_flags"]):
                status = "present_failed_filter"
            elif not allele_pass:
                status = "present_failed_filter"
            elif not valid:
                status = "present_missing_counts"
            elif ac < 1 or an < options["min_an"] or af < options["min_af"]:
                status = "present_below_threshold"
            else:
                status = "observed_pass"
            pops = {}
            for pop, (pac, pan) in pcounts.items():
                a, n = pac[index], pan[index]
                ok = _valid_counts(a, n, where + ":" + pop)
                pops[pop] = (a, n, a / n if ok else None)
            terms, genes = annotations.get(alt, (set(), set()))
            entry = dict(alt=alt, ac=ac, an=an, af=af, status=status, pops=pops,
                         original_indices={index + 1}, terms=set(terms), genes=set(genes), flags=flags)
            if alt in alts:
                old = alts[alt]
                compare = ("ac", "an", "status", "pops", "flags")
                if any(old[k] != entry[k] for k in compare):
                    raise ValueError(f"Conflicting duplicate ALT at {where}; do not pool releases/cohorts or select maximum AF")
                old["original_indices"].add(index + 1)
                old["terms"].update(terms)
                old["genes"].update(genes)
                qc["identical_duplicate_alt_entries"] += 1
            else:
                alts[alt] = entry
                qc["alt_" + status] += 1
    return (next(iter(refs)) if refs else None), alts, has_non_snv, region_flags


def population_site(alts, populations, min_an, min_ac=1, min_af=0.0):
    """Describe coexistence only where all currently retained ALT counts exist."""
    evaluations, multi_pops = {}, []
    for pop in populations:
        counts = [a["pops"].get(pop, (None, None, None)) for a in alts]
        if any(ac is None or an is None or an <= 0 or an < min_an or af is None for ac, an, af in counts):
            continue
        n = sum(ac >= min_ac and af is not None and af >= min_af for ac, an, af in counts)
        evaluations[pop] = n
        if n >= 2:
            multi_pops.append(pop)
    if len(alts) < 2:
        return "single_alt", evaluations, [], None, None, 0
    if multi_pops:
        classification = "within_group_observed"
    elif len(evaluations) >= 2:
        observed_sets = [{a["alt"] for a in alts
                          if (a["pops"][p][0] or 0) >= min_ac
                          and a["pops"][p][2] is not None
                          and a["pops"][p][2] >= min_af} for p in evaluations]
        union = set().union(*observed_sets)
        classification = "distributed_across_groups" if len(union) >= 2 else "not_resolved_with_available_groups"
    else:
        classification = "insufficient_group_callability"
    jaccard = cosine = None
    profile_n = 0
    if len(alts) >= 2:
        common = list(evaluations)
        one = {p for p in common if alts[0]["pops"][p][0] >= min_ac and alts[0]["pops"][p][2] >= min_af}
        two = {p for p in common if alts[1]["pops"][p][0] >= min_ac and alts[1]["pops"][p][2] >= min_af}
        union = one | two
        if union:
            jaccard = len(one & two) / len(union)
        # Equal ancestry weighting, no imputed missing group values.
        profile_n = len(common)
        if profile_n >= 3:
            x = [alts[0]["pops"][p][2] for p in common]
            y = [alts[1]["pops"][p][2] for p in common]
            denominator = math.sqrt(sum(z*z for z in x) * sum(z*z for z in y))
            if denominator:
                cosine = sum(a*b for a, b in zip(x, y)) / denominator
    return classification, evaluations, multi_pops, jaccard, cosine, profile_n


def scan_file(path, shard_dir, options, expected_chrom=None):
    """Produce one atomic-commit-ready shard; orchestrator owns directory commit."""
    candidate_storage = options.get("candidate_storage", "eager")
    if candidate_storage not in ("eager", "deferred"):
        raise ValueError("candidate_storage must be eager or deferred")
    path, shard_dir = Path(path), Path(shard_dir)
    shard_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    last_report = t0
    header = read_header(path, options["backend"])
    validate_count_header(header)
    populations = discover_populations(header, options["populations"])
    vep_tag, vep_fields = vep_schema(header, options["vep"])
    fields = ["AC", "AN"] + (["AF"] if "AF" in header["info"] else [])
    fields += [x + "_" + pop for pop in populations for x in ("AC", "AN")]
    fields += [x for x in ("AS_FilterStatus", "lcr", "segdup") if x in header["info"]]
    if vep_tag:
        fields.append(vep_tag)
    thresholds = options["thresholds"]
    site_fields = SITE_FIELDS + ["n_alt_ac" + str(t) for t in thresholds]
    allele_fields = ALLELE_FIELDS + [x + "_" + pop for pop in populations for x in ("ac", "an", "af")]
    scope_region = parse_region(options.get("region"))
    expected_chrom = canonical_chrom(expected_chrom) if expected_chrom else None
    qc = Counter()
    freq_counts, pop_counts, population_classes = Counter(), Counter(), Counter()
    site_annotation_counts, relationship_counts = Counter(), Counter()
    site_ac_spectrum = Counter()
    strata = Counter()
    current_window = None
    scope_chroms = set()
    n_sites = n_multi = 0

    with ExitStack() as stack:
        pass
        arch = from_options(options, expected_chrom, shard_dir, header)
        if arch: stack.callback(arch.close)
        fa = stack.enter_context(FastaContext(options["fasta"]))
        if scope_region:
            rc, rs, re_ = scope_region
            if rc not in fa.lengths or re_ > fa.lengths[rc]:
                raise ValueError("Region is outside the reference FASTA bounds")
        allowed_mask = [expected_chrom] if expected_chrom else ([scope_region[0]] if scope_region else None)
        mask = IntervalMask(options["callable_bed"], allowed_mask) if options.get("callable_bed") else None
        writers = {}
        output_schemas = {"strata": STRATA_FIELDS, "multiallelic_sites": site_fields,
                          "multiallelic_alleles": allele_fields}
        if candidate_storage == "eager":
            output_schemas["candidates"] = CANDIDATE_FIELDS
        if options.get("write_all_sites"):
            output_schemas["all_sites"] = ["chrom", "pos", "ref", "context", "n_alt", "primary_alt",
                                           "primary_ac", "primary_an", "primary_af"]
        for name, columns in output_schemas.items():
            handle, writer = tsv_writer(shard_dir / (name + ".tsv.gz"), columns)
            stack.callback(handle.close)
            writers[name] = writer

        def flush_strata():
            for key, counts in sorted(strata.items()):
                chrom, start0, end0, threshold, context, change, fbin, nbin = key
                writers["strata"].writerow(dict(zip(STRATA_FIELDS, (
                    chrom, start0, end0, threshold, context, change, fbin, nbin,
                    counts[0], counts[1], counts[2]))))
            strata.clear()

        stream = iterate_records(path, fields, options.get("region"), options["backend"], 1)
        for coord, group in itertools.groupby(stream, key=lambda r: (canonical_chrom(r["chrom"]), r["pos"])):
            chrom, pos = coord
            if expected_chrom and chrom != expected_chrom:
                raise ValueError(f"Filename/CHROM mismatch in {path}: expected {expected_chrom}, observed {chrom}")
            if scope_region and (chrom != scope_region[0] or not scope_region[1] <= pos <= scope_region[2]):
                raise ValueError("I/O returned a record outside the requested region")
            scope_chroms.add(chrom)
            ref, all_alts, non_snv, flags = group_alleles(group, header, populations, options, qc, vep_tag, vep_fields)
            if ref is None or not all_alts:
                continue
            if chrom not in fa.lengths or pos > fa.lengths[chrom]:
                raise ValueError(f"VCF position absent from FASTA: {chrom}:{pos}")
            if mask and not mask.contains(chrom, pos):
                qc["loci_outside_callable_mask"] += 1
                continue
            observed = sorted((a for a in all_alts.values() if a["status"] == "observed_pass"),
                              key=lambda a: (-a["af"], -a["ac"], a["alt"]))
            if not observed:
                qc["loci_without_qualifying_alt"] += 1
                continue
            context = fa.context(chrom, pos).upper()
            if len(context) != 3 or context[1] != ref:
                raise ValueError(f"REF/FASTA mismatch at {chrom}:{pos}:{ref}, context={context}; check build")
            canonical_ctx, _ = canonical_context(context, ref)
            if any(x not in BASES for x in context):
                qc["loci_with_ambiguous_context"] += 1
            cpg = int((ref == "C" and context[2] == "G") or (ref == "G" and context[0] == "C"))
            primary = observed[0]
            n = len(observed)
            if not 1 <= n <= 3:
                raise ValueError("More than three distinct SNV ALTs at a nucleotide")
            n_sites += 1
            n_multi += int(n >= 2)
            # One joint observation per base-eligible site, before the AC
            # sensitivity thresholds. AF rank need not agree with AC rank when
            # split records have different AN; inspect the retained AC values.
            # The common single-ALT case needs no extra allele loop.
            smallest_ac = primary["ac"]
            n_singletons = int(smallest_ac == 1)
            if n > 1:
                for i in range(1, n):
                    ac = observed[i]["ac"]
                    n_singletons += int(ac == 1)
                    if ac < smallest_ac:
                        smallest_ac = ac
            site_ac_spectrum[(chrom, n, n_singletons, ac_bin(smallest_ac))] += 1
            if n_sites % 100000 == 0 and time.monotonic() - last_report >= 30:
                last_report = time.monotonic()
                print(f"[scan {path.name}] {n_sites:,} qualifying loci; {n_multi:,} multi-ALT; "
                      f"{(last_report-t0)/60:.1f} min", flush=True)
            window_start = ((pos - 1) // options["window_size"]) * options["window_size"]
            window_end = min(window_start + options["window_size"], fa.lengths[chrom])
            if scope_region:
                window_start = max(window_start, scope_region[1] - 1)
                window_end = min(window_end, scope_region[2])
            window = chrom, window_start, window_end
            if current_window is not None and window != current_window:
                flush_strata()
            current_window = window
            counts_by_threshold = {}
            for threshold in thresholds:
                retained = [a for a in observed if a["ac"] >= threshold]
                m = len(retained)
                counts_by_threshold[threshold] = m
                if not retained:
                    continue
                anchor = retained[0]
                _, change = canonical_context(context, ref, anchor["alt"])
                key = (*window, threshold, canonical_ctx, change, af_bin(anchor["af"]), an_bin(anchor["an"]))
                if key not in strata:
                    strata[key] = [0, 0, 0]
                strata[key][0] += 1
                strata[key][1] += int(m >= 2)
                strata[key][2] += int(m == 3)
                for rank, allele in enumerate(retained, 1):
                    freq_counts[(chrom, threshold, rank, ac_bin(allele["ac"]), af_bin(allele["af"]))] += 1
                if populations:
                    classification, evaluations, multis, _, _, _ = population_site(
                        retained, populations, options["min_pop_an"], threshold, options["min_pop_af"])
                    if m >= 2:
                        population_classes[(chrom, threshold, classification)] += 1
                    for pop, count in evaluations.items():
                        # n_evaluable includes known zeros but no missing-population fields.
                        pop_counts[(chrom, threshold, pop, "n_evaluable")] += 1
                        if count >= 1:
                            pop_counts[(chrom, threshold, pop, "n_sites")] += 1
                        if count >= 2:
                            pop_counts[(chrom, threshold, pop, "n_multi")] += 1
                        if count == 3:
                            pop_counts[(chrom, threshold, pop, "n_three")] += 1
                if vep_tag:
                    # This is the high-frequency ALT annotation, not an independent site category.
                    terms = "&".join(sorted(anchor["terms"])) or "unannotated"
                    site_annotation_counts[(chrom, threshold, terms, "n_sites")] += 1
                    site_annotation_counts[(chrom, threshold, terms, "n_multi")] += int(m >= 2)
            if arch:
                arch_mode = population_site(observed, populations, options["min_pop_an"], 1, options["min_pop_af"])[0] if n >= 2 and populations else "NA"
                arch.observe(chrom=chrom, pos1=pos, n_alt=n, context=canonical_ctx,
                    primary_alt_af=primary["af"], primary_alt_an=primary["an"],
                    threshold_counts=counts_by_threshold, population_mode=arch_mode, region_flags=flags)
            if "all_sites" in writers:
                writers["all_sites"].writerow(dict(chrom=chrom, pos=pos, ref=ref, context=context,
                    n_alt=n, primary_alt=primary["alt"], primary_ac=primary["ac"],
                    primary_an=primary["an"], primary_af=fmt(primary["af"])))
            if n < 2:
                continue
            qc["multi_alt_loci"] += 1
            if n == 3:
                qc["three_alt_loci"] += 1
            an_values = [a["an"] for a in observed]
            if len(set(an_values)) > 1:
                qc["multi_alt_loci_discordant_an"] += 1
            primary_tie = int(abs(observed[0]["af"] - observed[1]["af"]) < 1e-12)
            pop_class, evaluations, multi_pops, jac, cos, npop = population_site(
                observed, populations, options["min_pop_an"], 1, options["min_pop_af"])
            secondary = observed[1]
            ref_relation = relation(ref, secondary["alt"])
            primary_relation = relation(primary["alt"], secondary["alt"])
            relationship_counts[(chrom, ref_relation, primary_relation, primary_tie)] += 1
            _, pchange = canonical_context(context, ref, primary["alt"])
            same_terms = None
            if primary["terms"] and secondary["terms"]:
                same_terms = int(primary["terms"] == secondary["terms"])
            row = dict(chrom=chrom, pos=pos, ref=ref, context=context, canonical_context=canonical_ctx,
                       cpg_site=cpg, n_alt=n, alts_by_af=",".join(a["alt"] for a in observed),
                       alt_ac=",".join(str(a["ac"]) for a in observed),
                       alt_an=",".join(str(a["an"]) for a in observed),
                       alt_af=",".join(fmt(a["af"]) for a in observed), primary_alt=primary["alt"],
                       secondary_alt=secondary["alt"], primary_tie=primary_tie, an_min=min(an_values),
                       an_max=max(an_values), an_relative_spread=fmt((max(an_values)-min(an_values))/max(an_values)),
                       has_failed_alt=int(any(a["status"] == "present_failed_filter" for a in all_alts.values())),
                       has_missing_alt=int(any(a["status"] == "present_missing_counts" for a in all_alts.values())),
                       has_non_snv=int(non_snv), region_flags=",".join(sorted(flags)) or ".",
                       primary_change=pchange, secondary_ref_relation=ref_relation,
                       secondary_primary_relation=primary_relation, population_class=pop_class,
                       n_pop_evaluable=len(evaluations), pops_with_multi_alt=",".join(multi_pops) or ".",
                       pop_presence_jaccard=fmt(jac), pop_profile_cosine=fmt(cos), pop_profile_n=npop,
                       same_consequence_terms=fmt(same_terms),
                       ref_frequency_status="not_inferable_from_sites_summary")
            row.update({"n_alt_ac" + str(t): counts_by_threshold[t] for t in thresholds})
            writers["multiallelic_sites"].writerow(row)
            ranks = {a["alt"]: i + 1 for i, a in enumerate(observed)}
            for alt, a in sorted(all_alts.items()):
                allele_row = dict(chrom=chrom, pos=pos, ref=ref, alt=alt, status=a["status"],
                    alt_frequency_rank=ranks.get(alt, "."), is_primary=int(alt == primary["alt"]),
                    ac=fmt(a["ac"]), an=fmt(a["an"]), af=fmt(a["af"]),
                    original_alt_indices=",".join(map(str, sorted(a["original_indices"]))),
                    consequence_terms="&".join(sorted(a["terms"])) or ".",
                    gene_symbols=",".join(sorted(a["genes"])) or ".", region_flags=",".join(sorted(a["flags"])) or ".")
                for pop in populations:
                    pac, pan, paf = a["pops"][pop]
                    allele_row.update({"ac_"+pop: fmt(pac), "an_"+pop: fmt(pan), "af_"+pop: fmt(paf)})
                writers["multiallelic_alleles"].writerow(allele_row)
            if "candidates" in writers:
                for alt in BASES:
                    if alt == ref:
                        continue
                    a = all_alts.get(alt)
                    trans = int(relation(ref, alt) == "transition")
                    writers["candidates"].writerow(dict(
                        chrom=chrom, pos=pos, ref=ref, alt=alt,
                        status=a["status"] if a else "not_listed_in_release",
                        is_primary=int(alt == primary["alt"]), primary_tie=primary_tie, n_alt=n,
                        ac=fmt(a["ac"] if a else None), an=fmt(a["an"] if a else None),
                        af=fmt(a["af"] if a else None), context=context, primary_alt=primary["alt"],
                        transition=trans, cpg_transition=int(cpg and trans)))
        flush_strata()
        if arch: arch.write(shard_dir / options.get("architecture_prefix", "architecture"))
        if expected_chrom:
            scope_chroms.add(expected_chrom)
        if scope_region:
            scope_chroms = {scope_region[0]}
        # For a named chr file even an empty chromosome has a valid scan scope.
        # Never create zeros for the other chromosomes merely listed in its header.
        if not scope_chroms:
            conventional = [canonical_chrom(c) for c in header["contigs"] if canonical_chrom(c) in fa.lengths]
            if len(conventional) == 1:
                scope_chroms.update(conventional)
        scope_handle, scope_writer = tsv_writer(shard_dir / "scan_scope.tsv", ["chrom", "start0", "end0"])
        with scope_handle:
            for chrom in sorted(scope_chroms, key=chrom_key):
                if chrom not in fa.lengths:
                    continue
                start, end = (scope_region[1]-1, scope_region[2]) if scope_region else (0, fa.lengths[chrom])
                scope_writer.writerow(dict(chrom=chrom, start0=start, end0=end))

    def write_counter(name, columns, counter):
        handle, writer = tsv_writer(shard_dir / name, columns + ["count"])
        with handle:
            for key, value in sorted(counter.items()):
                writer.writerow(dict(zip(columns + ["count"], (*key, value))))

    write_counter("frequency_counts.tsv", ["chrom", "threshold", "alt_rank", "ac_bin", "af_bin"], freq_counts)
    write_counter("site_ac_spectrum.tsv", SITE_AC_SPECTRUM_FIELDS, site_ac_spectrum)
    write_counter("population_counts.tsv", ["chrom", "threshold", "population", "metric"], pop_counts)
    write_counter("population_classes.tsv", ["chrom", "threshold", "classification"], population_classes)
    write_counter("annotation_counts.tsv", ["chrom", "threshold", "primary_alt_terms", "metric"], site_annotation_counts)
    write_counter("relationship_counts.tsv", ["chrom", "secondary_ref_relation", "secondary_primary_relation", "primary_tie"], relationship_counts)
    qc["eligible_sites"] = n_sites
    qc["multi_alt_sites"] = n_multi
    result = dict(version=VERSION, input=str(path.resolve()), populations=populations, vep_tag=vep_tag,
                  header_sample_count=header["samples"], qc=dict(qc), elapsed_seconds=time.monotonic()-t0,
                  fields_read=fields, scope=sorted(scope_chroms, key=chrom_key),
                  site_ac_spectrum=dict(schema=1, status="measured", counted_sites=n_sites,
                      scope="base_eligible_sites_before_ac_threshold_sensitivity",
                      filename="site_ac_spectrum.tsv"),
                  candidate_storage=candidate_storage, expected_candidate_rows=3*n_multi,
                  reference_frequency="not inferred", complete=True)
    dump_json(shard_dir / "scan.json", result)
    return result


# ========================================================================
# Spectrum
# ========================================================================

# Small, fixed-cohort summaries of ALT-count spectra from verified scan shards.
# 
# Older caches remain usable but cannot manufacture a joint AC distribution that
# was not saved. Counts from partial shard coverage are explicitly labelled.


from collections import Counter
import csv
import json
from pathlib import Path

from allele_common import dump_json, tsv_writer

FIELDS = ["chrom", "n_alt", "n_singleton_alts", "min_alt_ac_bin", "count"]
AC_BINS = {"1", "2_to_4", "5_to_9", "10_to_99", "100_to_999", "ge_1000"}


def merge_site_spectra(shards, out_dir):
    """Merge only measured spectra; callers first validate their scan run.

    A complete result means complete for the supplied run, not for all human
    chromosomes. Old caches are not modified, reinterpreted or rescanned.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    available, missing = [], []
    for shard in map(Path, shards):
        record = json.loads((shard / "scan.json").read_text())
        meta = record.get("site_ac_spectrum")
        if not meta:
            missing.append(str(shard))
            continue
        if (meta.get("schema") != 1 or meta.get("status") != "measured"
                or meta.get("filename") != "site_ac_spectrum.tsv"):
            raise ValueError(f"Unsupported site AC spectrum metadata: {shard}")
        subtotal, seen = 0, set()
        with (shard / "site_ac_spectrum.tsv").open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != FIELDS:
                raise ValueError(f"Site AC spectrum columns disagree: {shard}")
            for row in reader:
                try:
                    n, singletons, count = (int(row[x]) for x in ("n_alt", "n_singleton_alts", "count"))
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Malformed site AC spectrum: {shard}") from exc
                acbin = row["min_alt_ac_bin"]
                if (not 1 <= n <= 3 or not 0 <= singletons <= n or count < 1
                        or acbin not in AC_BINS or (acbin == "1") != (singletons > 0)):
                    raise ValueError(f"Inconsistent site AC spectrum category: {shard}")
                key = (row["chrom"], n, singletons, acbin)
                if key in seen:
                    raise ValueError(f"Repeated site AC spectrum category: {shard}")
                seen.add(key)
                counts[key] += count
                subtotal += count
        if subtotal != meta["counted_sites"] or subtotal != record["qc"].get("eligible_sites", 0):
            raise ValueError(f"Site AC spectrum denominator disagrees with scan: {shard}")
        available.append(str(shard))
    status = "complete" if available and not missing else "partial" if available else "not_available"
    handle, writer = tsv_writer(out / "site_ac_spectrum.tsv", FIELDS)
    with handle:
        for key, value in sorted(counts.items()):
            writer.writerow(dict(zip(FIELDS, (*key, value))))

    # These contrasts retain the ORIGINAL base-eligible multi-ALT denominator.
    # They are not the n_multi/n_sites ratios of independently filtered runs.
    groups = {}
    for (chrom, n, singletons, acbin), count in counts.items():
        targets = [chrom]
        if chrom in {str(i) for i in range(1, 23)}:
            targets.append("autosomes_observed")
        for group in targets:
            c = groups.setdefault(group, Counter())
            c["base_eligible_sites"] += count
            if n < 2:
                continue
            c["base_multi_alt_sites"] += count
            c["multi_with_singleton_alt"] += count * (singletons > 0)
            c["multi_all_alt_ac_ge2"] += count * (singletons == 0)
            c["multi_still_multi_ac_ge2"] += count * (n - singletons >= 2)
            c["multi_lost_multi_ac_ge2"] += count * (n - singletons < 2)
            c["multi_all_alt_ac_ge5"] += count * (acbin not in {"1", "2_to_4"})
            c["multi_all_alt_ac_ge10"] += count * (acbin not in {"1", "2_to_4", "5_to_9"})
    metrics = ["base_eligible_sites", "base_multi_alt_sites", "multi_with_singleton_alt",
               "multi_all_alt_ac_ge2", "multi_still_multi_ac_ge2", "multi_lost_multi_ac_ge2",
               "multi_all_alt_ac_ge5", "multi_all_alt_ac_ge10"]
    fields = ["chrom", "coverage_status"] + metrics + ["fraction_multi_lost_ac_ge2", "fraction_all_alt_ac_ge2"]
    handle, writer = tsv_writer(out / "singleton_sensitivity.tsv", fields)
    with handle:
        for group, values in sorted(groups.items()):
            denominator = values["base_multi_alt_sites"]
            row = dict(chrom=group, coverage_status=status, **{key: values[key] for key in metrics})
            row["fraction_multi_lost_ac_ge2"] = values["multi_lost_multi_ac_ge2"] / denominator if denominator else "."
            row["fraction_all_alt_ac_ge2"] = values["multi_all_alt_ac_ge2"] / denominator if denominator else "."
            writer.writerow(row)
    result = dict(schema=1, status=status, total_shards=len(shards), measured_shards=available,
                  unavailable_shards=missing, counted_sites=sum(counts.values()) if available else None,
                  counted_chromosomes=sorted({key[0] for key in counts}),
                  denominator="original base-eligible multi-ALT sites among measured shards",
                  missing_coverage="not measured; not zero; no automatic VCF rescan",
                  interpretation="Removing AC=1 ALT is a count sensitivity, not longitudinal allele loss or proof of equilibrium.")
    dump_json(out / "site_ac_spectrum_availability.json", result)
    return result


# ========================================================================
# Summarize
# ========================================================================

# Stream compact occupancy strata into descriptive, mutation-context summaries.
# 
# The expected values are chromosome-held-out *conditional observation* rates.
# They are neither mutation rates nor a selection test. No network calls are made.
# Only Python's standard library is required; matplotlib is an optional renderer.


import argparse
import csv
import gzip
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable


REQUIRED = {
    "chrom", "start0", "end0", "threshold", "context", "primary_change",
    "primary_af_bin", "an_bin", "n_sites", "n_multi", "n_three",
}
LEVELS = ("exact", "context_change", "threshold")
AUTOSOMES = {str(i) for i in range(1, 23)}


def summary_chrom(value):
    value = str(value).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value.isdigit():
        return str(int(value))
    return "MT" if value in {"M", "MT"} else value


def summary_chrom_key(chrom):
    if chrom in AUTOSOMES:
        return (0, int(chrom))
    return (1, {"X": 0, "Y": 1, "MT": 2}.get(chrom, 3), chrom)


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith((".gz", ".bgz", ".bgzip")) else open(path, encoding="utf-8")


def summary_region(region):
    """Return (canonical chromosome, zero-based start, exclusive end).

    A string uses the usual VCF/samtools one-based inclusive chr:start-end form.
    A tuple/list is explicitly zero-based half-open. A chromosome-only string
    selects that chromosome without bounding its coordinate span.
    """
    if region is None or region == "":
        return None
    if isinstance(region, (tuple, list)) and len(region) == 3:
        chrom, start, end = region
        start, end = int(start), int(end)
    elif isinstance(region, str):
        if ":" not in region:
            return summary_chrom(region), 0, None
        chrom, coordinates = region.rsplit(":", 1)
        if "-" not in coordinates:
            raise ValueError("region must be chr:start-end (one-based inclusive)")
        start, end = coordinates.replace(",", "").split("-", 1)
        start, end = int(start) - 1, int(end)
    else:
        raise ValueError("region must be a chr:start-end string or (chrom,start0,end0)")
    if start < 0 or end <= start:
        raise ValueError("region has invalid coordinates")
    return summary_chrom(chrom), start, end


def merge_intervals(intervals):
    out = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError(f"invalid zero-based interval: {start}, {end}")
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out


def intersect_intervals(left, right):
    """Intersection of two sorted disjoint interval collections."""
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        if max(a, c) < min(b, d):
            result.append((max(a, c), min(b, d)))
        if b <= d:
            i += 1
        else:
            j += 1
    return result


def load_bed(path):
    intervals = defaultdict(list)
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: BED requires at least 3 columns")
            chrom = summary_chrom(fields[0])
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid BED coordinates") from exc
            intervals[chrom].append((start, end))
    return {chrom: merge_intervals(rows) for chrom, rows in intervals.items()}


def _find_strata(shard_dirs):
    files = []
    for item in shard_dirs:
        item = Path(item)
        if item.is_file():
            path = item
        else:
            candidates = [item / "strata.tsv.gz", item / "strata.tsv"]
            path = next((p for p in candidates if p.is_file()), None)
            if path is None:
                raise FileNotFoundError(f"No strata.tsv.gz or strata.tsv in {item}")
        files.append(path)
    if not files:
        raise ValueError("at least one strata file is required")
    if len({str(p.resolve()) for p in files}) != len(files):
        raise ValueError("duplicate strata input file would double-count observations")
    return files


def read_strata(files, window_size, region):
    for path in files:
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = REQUIRED - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{path}: missing columns: {','.join(sorted(missing))}")
            for line_number, row in enumerate(reader, 2):
                try:
                    for column in ("start0", "end0", "threshold", "n_sites", "n_multi", "n_three"):
                        row[column] = int(row[column])
                except (ValueError, TypeError) as exc:
                    raise ValueError(f"{path}:{line_number}: count/coordinate must be an integer") from exc
                row["chrom"] = summary_chrom(row["chrom"])
                start, end = row["start0"], row["end0"]
                if start < 0 or end <= start or start // window_size != (end - 1) // window_size:
                    raise ValueError(f"{path}:{line_number}: interval invalid or crosses a {window_size}-bp summary window")
                if not 0 <= row["n_three"] <= row["n_multi"] <= row["n_sites"]:
                    raise ValueError(f"{path}:{line_number}: require 0 <= n_three <= n_multi <= n_sites")
                if row["threshold"] < 1:
                    raise ValueError(f"{path}:{line_number}: threshold must be positive")
                if region:
                    rchrom, rstart, rend = region
                    if row["chrom"] != rchrom or end <= rstart or (rend is not None and start >= rend):
                        continue
                    # Counts cannot be clipped after aggregation. Root scan must
                    # already enforce region; bins may extend across its edge.
                row["window_start0"] = (start // window_size) * window_size
                yield row


def _count():
    return {"n_sites": 0, "n_multi": 0, "n_three": 0}


def _add_count(target, row):
    for field in ("n_sites", "n_multi", "n_three"):
        target[field] += row[field]


def _expected_count():
    return {"expected_multi": 0.0, "scored_sites": 0, "scored_observed_multi": 0,
            "exact_sites": 0, "context_change_sites": 0, "threshold_sites": 0,
            "unscored_sites": 0, "training_sites_min": None, "training_sites_max": None,
            "unavailable_reasons": set()}


def _add_expected(target, row, prediction):
    if prediction["rate"] is None:
        target["unscored_sites"] += row["n_sites"]
        target["unavailable_reasons"].add(prediction["reason"])
        return
    target["expected_multi"] += row["n_sites"] * prediction["rate"]
    target["scored_sites"] += row["n_sites"]
    target["scored_observed_multi"] += row["n_multi"]
    target[prediction["level"] + "_sites"] += row["n_sites"]
    ntrain = prediction["n_train"]
    target["training_sites_min"] = ntrain if target["training_sites_min"] is None else min(target["training_sites_min"], ntrain)
    target["training_sites_max"] = ntrain if target["training_sites_max"] is None else max(target["training_sites_max"], ntrain)


def _keys(row):
    exact = (row["threshold"], row["context"], row["primary_change"], row["primary_af_bin"], row["an_bin"])
    return exact, exact[:3], exact[:1]


class HeldOutRates:
    """Compact sufficient counts; target chromosome is subtracted at all levels."""

    def __init__(self, min_training_sites, restricted=False):
        self.min_training_sites = min_training_sites
        self.restricted = restricted
        self.total = defaultdict(lambda: [0, 0])
        self.by_chrom = defaultdict(lambda: [0, 0])
        self.support = defaultdict(set)
        self.autosomes = set()
        self.cache = {}

    def observe(self, row):
        chrom = row["chrom"]
        if chrom not in AUTOSOMES or row["n_sites"] == 0:
            return
        self.autosomes.add(chrom)
        for level, key in zip(LEVELS, _keys(row)):
            gkey = (level, key)
            ckey = (chrom, level, key)
            self.total[gkey][0] += row["n_sites"]
            self.total[gkey][1] += row["n_multi"]
            self.by_chrom[ckey][0] += row["n_sites"]
            self.by_chrom[ckey][1] += row["n_multi"]
            self.support[gkey].add(chrom)

    def predict(self, row):
        chrom = row["chrom"]
        cache_key = (chrom, _keys(row)[0])
        if cache_key in self.cache:
            return self.cache[cache_key]
        reason = None
        if chrom not in AUTOSOMES:
            reason = "non_autosome_descriptive_only"
        elif self.restricted:
            reason = "region_restricted_no_cross_chromosome_background"
        elif len(self.autosomes) < 2:
            reason = "fewer_than_two_represented_autosomes"
        if reason is None:
            for level, key in zip(LEVELS, _keys(row)):
                total = self.total.get((level, key), (0, 0))
                own = self.by_chrom.get((chrom, level, key), (0, 0))
                ntrain, multitrain = total[0] - own[0], total[1] - own[1]
                if ntrain >= self.min_training_sites:
                    result = {"rate": (multitrain + 0.5) / (ntrain + 1), "level": level,
                              "n_train": ntrain, "n_train_multi": multitrain,
                              "n_training_autosomes": len(self.support[(level, key)] - {chrom}), "reason": ""}
                    self.cache[cache_key] = result
                    return result
            reason = "insufficient_held_out_training_sites"
        result = {"rate": None, "level": "unavailable", "n_train": None,
                  "n_train_multi": None, "n_training_autosomes": 0, "reason": reason}
        self.cache[cache_key] = result
        return result


def wilson_interval(success, total, z=1.959963984540054):
    if total == 0:
        return None, None
    p = success / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    width = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    return max(0.0, center - width), min(1.0, center + width)


def _raw_fields(counts):
    n, multi, three = counts["n_sites"], counts["n_multi"], counts["n_three"]
    low, high = wilson_interval(multi, n)
    return {**counts, "n_single_alt": n - multi, "n_two_alt": multi - three,
            "multi_fraction": multi / n if n else None, "multi_wilson_low": low, "multi_wilson_high": high,
            "three_alt_fraction": three / n if n else None}


def _expected_fields(raw, expected):
    expected = expected or _expected_count()
    scored, unscored = expected["scored_sites"], expected["unscored_sites"]
    e = expected["expected_multi"] if scored else None
    complete = scored > 0 and unscored == 0
    return {
        "expected_status": "complete" if complete else "partial" if scored else "unavailable",
        "expected_multi": e,
        "observed_over_expected": raw["n_multi"] / e if complete and e else None,
        "scored_sites": scored,
        "scored_observed_multi": expected["scored_observed_multi"],
        "scored_observed_over_expected": expected["scored_observed_multi"] / e if e else None,
        "expected_coverage": scored / raw["n_sites"] if raw["n_sites"] else None,
        "training_exact_sites": expected["exact_sites"],
        "training_context_change_sites": expected["context_change_sites"],
        "training_threshold_sites": expected["threshold_sites"],
        "unscored_sites": unscored,
        "training_sites_min": expected["training_sites_min"],
        "training_sites_max": expected["training_sites_max"],
        "expected_unavailable_reason": ";".join(sorted(expected["unavailable_reasons"])),
    }


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else row.get(key, "") for key in fieldnames})


RAW_COLUMNS = ["n_sites", "n_multi", "n_three", "n_single_alt", "n_two_alt", "multi_fraction",
               "multi_wilson_low", "multi_wilson_high", "three_alt_fraction"]
EXPECTED_COLUMNS = ["expected_status", "expected_multi", "observed_over_expected", "scored_sites",
                    "scored_observed_multi", "scored_observed_over_expected", "expected_coverage",
                    "training_exact_sites", "training_context_change_sites", "training_threshold_sites",
                    "unscored_sites", "training_sites_min", "training_sites_max", "expected_unavailable_reason"]


def _load_scan_scopes(files):
    """Optional root-written scope records prevent unscanned intervals becoming zeros."""
    intervals = defaultdict(list)
    found = 0
    for directory in dict.fromkeys(p.parent for p in files):
        path = directory / "scan_scope.tsv"
        if not path.is_file():
            continue
        found += 1
        with open(path, encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"chrom", "start0", "end0"}.issubset(reader.fieldnames or []):
                raise ValueError(f"{path}: expected chrom,start0,end0 columns")
            for row in reader:
                intervals[summary_chrom(row["chrom"])].append((int(row["start0"]), int(row["end0"])))
    total = len(set(p.parent for p in files))
    if found not in (0, total):
        raise ValueError("scan_scope.tsv must be supplied for every shard or for none")
    return ({chrom: merge_intervals(values) for chrom, values in intervals.items()}, bool(found))


def _callable_windows(path, scopes, has_scopes, seen_chroms, window_size, region, warnings):
    bed = load_bed(path)
    allowed_chroms = set(scopes) if has_scopes else set(seen_chroms)
    if region:
        allowed_chroms = {region[0]} & (allowed_chroms or {region[0]})
    if not has_scopes:
        warnings.append("No scan_scope.tsv: callable zero windows are limited to chromosomes represented in strata (or the explicit region). Root scan must have covered each such chromosome/region completely.")
    by_window = defaultdict(int)
    by_chrom = defaultdict(int)
    for chrom in sorted(allowed_chroms, key=summary_chrom_key):
        intervals = bed.get(chrom, [])
        if has_scopes:
            intervals = intersect_intervals(intervals, scopes.get(chrom, []))
        if region:
            _, start, end = region
            limit = end if end is not None else max((b for _, b in intervals), default=start)
            intervals = intersect_intervals(intervals, [(start, limit)]) if limit > start else []
        for start, end in intervals:
            by_chrom[chrom] += end - start
            while start < end:
                wstart = (start // window_size) * window_size
                stop = min(end, wstart + window_size)
                by_window[(chrom, wstart)] += stop - start
                start = stop
    return by_window, by_chrom


def _make_plots(out_dir, chromosome_rows, window_rows, context_rows, threshold_rows):
    try:
        os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "allele3-matplotlib"))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return [], ["matplotlib is unavailable; all machine-readable tables were produced"]
    paths = []
    if not threshold_rows:
        return paths, ["No strata rows; no figures were rendered"]
    thresholds = [r["threshold"] for r in threshold_rows]
    selected = thresholds[0]
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

    def save(fig, name):
        fig.tight_layout()
        for extension in ("png", "pdf"):
            path = out_dir / f"{name}.{extension}"
            fig.savefig(path, dpi=180, bbox_inches="tight")
            paths.append(str(path))
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    chroms = sorted({r["chrom"] for r in chromosome_rows}, key=summary_chrom_key)
    xmap = {chrom: i for i, chrom in enumerate(chroms)}
    for threshold in thresholds:
        rows = [r for r in chromosome_rows if r["threshold"] == threshold and r["multi_fraction"] is not None]
        ax.plot([xmap[r["chrom"]] for r in rows], [100 * r["multi_fraction"] for r in rows], marker="o", markersize=3, linewidth=1, label=f"AC ≥ {threshold}")
    ax.set_xticks(range(len(chroms)), chroms)
    ax.set(xlabel="Chromosome", ylabel="Sites with ≥2 observed ALT alleles (%)", title="Multi-ALT fraction among eligible SNV sites")
    ax.legend(frameon=False)
    save(fig, "chromosome_occupancy")

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True, gridspec_kw={"height_ratios": [1.2, 1]})
    rows = [r for r in window_rows if r["threshold"] == selected]
    offsets = {}
    cursor = 0
    ticks, ticklabels = [], []
    spans = [max(r["end0"] for r in rows if r["chrom"] == c) -
             min(r["start0"] for r in rows if r["chrom"] == c)
             for c in chroms if any(r["chrom"] == c for r in rows)]
    gap = max(1, min(1_000_000, int(max(spans, default=1) * .02)))
    for index, chrom in enumerate(chroms):
        subset = [r for r in rows if r["chrom"] == chrom]
        if not subset:
            continue
        low, high = min(r["start0"] for r in subset), max(r["end0"] for r in subset)
        offsets[chrom] = cursor - low
        ticks.append(cursor + (high - low) / 2)
        ticklabels.append(chrom)
        color = "#306FA2" if index % 2 == 0 else "#67A998"
        valid = [r for r in subset if r["multi_fraction"] is not None]
        top.scatter([r["start0"] + offsets[chrom] for r in valid], [100 * r["multi_fraction"] for r in valid], s=7, color=color)
        modeled = [r for r in subset if r["observed_over_expected"] is not None]
        bottom.scatter([r["start0"] + offsets[chrom] for r in modeled], [r["observed_over_expected"] for r in modeled], s=7, color=color)
        cursor += high - low + gap
    top.set(ylabel="Multi-ALT sites (%)", title=f"Conditional genome distribution · AC ≥ {selected}")
    bottom.axhline(1, color="#777777", linewidth=0.8, linestyle="--")
    bottom.set(ylabel="Observed / held-out expected", xlabel="Chromosome")
    bottom.set_xticks(ticks, ticklabels)
    if not any(r["observed_over_expected"] is not None for r in rows):
        bottom.text(0.5, 0.55, "Cross-chromosome expectation unavailable", transform=bottom.transAxes, ha="center")
    save(fig, "window_distribution")

    context_aggregates = defaultdict(_count)
    for row in context_rows:
        if row["threshold"] == selected:
            _add_count(context_aggregates[row["context"]], row)
    contexts = sorted(context_aggregates)
    fig, ax = plt.subplots(figsize=(max(8, len(contexts) * 0.28), 4.6))
    ax.bar(contexts, [100 * context_aggregates[c]["n_multi"] / context_aggregates[c]["n_sites"] if context_aggregates[c]["n_sites"] else 0 for c in contexts], color="#306FA2")
    ax.tick_params(axis="x", labelrotation=90)
    ax.set(ylabel="Multi-ALT sites (%)", xlabel="Pyrimidine-canonical trinucleotide", title=f"Context distribution · AC ≥ {selected}")
    save(fig, "context_occupancy")
    return paths, []


def summarize(shard_dirs, out_dir, window_size=1_000_000, min_training_sites=100,
              make_plots=True, callable_bed=None, region=None, thresholds=None):
    """Summarize compact strata in two passes without loading a full dataframe.

    Root scanning must already enforce the requested region and callable BED;
    aggregated counts cannot be filtered by individual site coordinates here.
    Optional per-shard scan_scope.tsv bounds zero-window materialization.
    """
    if window_size <= 0 or min_training_sites < 1:
        raise ValueError("window_size and min_training_sites must be positive")
    if isinstance(shard_dirs, (str, Path)):
        shard_dirs = [shard_dirs]
    files = _find_strata(shard_dirs)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    parsed_region = summary_region(region)
    warnings = []
    raw_chr, raw_window = defaultdict(_count), defaultdict(_count)
    raw_context, raw_threshold = defaultdict(_count), defaultdict(_count)
    if thresholds is not None:
        for threshold in thresholds:
            threshold = int(threshold)
            if threshold < 1:
                raise ValueError("thresholds must be positive integers")
            raw_threshold[threshold]
    rates = HeldOutRates(min_training_sites, restricted=parsed_region is not None)
    seen_chroms, exact_rows = set(), defaultdict(_count)
    n_strata_rows = 0
    for row in read_strata(files, window_size, parsed_region):
        n_strata_rows += 1
        chrom, threshold = row["chrom"], row["threshold"]
        seen_chroms.add(chrom)
        _add_count(raw_chr[(chrom, threshold)], row)
        _add_count(raw_window[(chrom, row["window_start0"], threshold)], row)
        _add_count(raw_context[(threshold, row["context"], row["primary_change"])], row)
        _add_count(raw_threshold[threshold], row)
        _add_count(exact_rows[(chrom, _keys(row)[0])], row)
        rates.observe(row)
    scopes, has_scopes = _load_scan_scopes(files)
    if has_scopes:
        for chrom in scopes:
            if parsed_region and chrom != parsed_region[0]:
                continue
            for threshold in raw_threshold:
                raw_chr[(chrom, threshold)]
    callable_window, callable_chr = {}, {}
    if callable_bed is not None:
        callable_window, callable_chr = _callable_windows(callable_bed, scopes, has_scopes, seen_chroms, window_size, parsed_region, warnings)
        # Every scanned callable window appears even if no variant passed QC.
        for (chrom, start), bases in callable_window.items():
            if bases:
                for threshold in raw_threshold:
                    raw_window[(chrom, start, threshold)]
                    raw_chr[(chrom, threshold)]
        for (chrom, start, threshold), counts in raw_window.items():
            if counts["n_sites"] > 0 and callable_window.get((chrom, start), 0) == 0:
                raise ValueError(f"observed sites in {chrom}:{start} have zero callable denominator; scan and summary masks/scopes differ")
            if counts["n_sites"] > callable_window.get((chrom, start), 0):
                raise ValueError(f"more eligible SNV sites than callable bases in {chrom}:{start}; check duplicate shards and scan/summary mask agreement")
    expected_chr, expected_window, expected_threshold = defaultdict(_expected_count), defaultdict(_expected_count), defaultdict(_expected_count)
    for row in read_strata(files, window_size, parsed_region):
        prediction = rates.predict(row)
        _add_expected(expected_chr[(row["chrom"], row["threshold"])], row, prediction)
        _add_expected(expected_window[(row["chrom"], row["window_start0"], row["threshold"])], row, prediction)
        _add_expected(expected_threshold[row["threshold"]], row, prediction)

    chromosome_rows, window_rows, context_rows, threshold_rows = [], [], [], []
    for (chrom, threshold), counts in sorted(raw_chr.items(), key=lambda item: (summary_chrom_key(item[0][0]), item[0][1])):
        bases = callable_chr.get(chrom, 0) if callable_bed is not None else None
        chromosome_rows.append({"chrom": chrom, "threshold": threshold, **_raw_fields(counts),
                                "callable_bases": bases, "multi_per_callable_mb": counts["n_multi"] * 1e6 / bases if bases else None,
                                **_expected_fields(counts, expected_chr.get((chrom, threshold)))})
    for (chrom, start, threshold), counts in sorted(raw_window.items(), key=lambda item: (summary_chrom_key(item[0][0]), item[0][1], item[0][2])):
        bases = callable_window.get((chrom, start), 0) if callable_bed is not None else None
        end = start + window_size
        # Coordinates label the complete bin; the callable denominator is clipped
        # to scan scope/region and may be smaller than the displayed bin width.
        window_rows.append({"chrom": chrom, "start0": start, "end0": end, "threshold": threshold,
                            **_raw_fields(counts), "callable_bases": bases,
                            "multi_per_callable_mb": counts["n_multi"] * 1e6 / bases if bases else None,
                            **_expected_fields(counts, expected_window.get((chrom, start, threshold)))})
    for (threshold, context, change), counts in sorted(raw_context.items()):
        context_rows.append({"threshold": threshold, "context": context, "primary_change": change, **_raw_fields(counts)})
    for threshold, counts in sorted(raw_threshold.items()):
        threshold_rows.append({"threshold": threshold, **_raw_fields(counts), **_expected_fields(counts, expected_threshold.get(threshold))})
    expected_context_rows = []
    for (chrom, key), counts in sorted(exact_rows.items(), key=lambda item: (summary_chrom_key(item[0][0]), item[0][1])):
        threshold, context, change, af_bin, an_bin = key
        row = {"chrom": chrom, "threshold": threshold, "context": context, "primary_change": change,
               "primary_af_bin": af_bin, "an_bin": an_bin}
        prediction = rates.predict(row)
        expected_context_rows.append({**row, "observed_sites": counts["n_sites"], "observed_multi": counts["n_multi"],
                                     "training_level": prediction["level"], "training_sites": prediction["n_train"],
                                     "training_multi": prediction["n_train_multi"], "training_autosomes": prediction["n_training_autosomes"],
                                     "held_out_rate": prediction["rate"],
                                     "expected_multi": counts["n_sites"] * prediction["rate"] if prediction["rate"] is not None else None,
                                     "unavailable_reason": prediction["reason"]})
    specs = [
        ("chromosome_summary.tsv", chromosome_rows, ["chrom", "threshold"] + RAW_COLUMNS + ["callable_bases", "multi_per_callable_mb"] + EXPECTED_COLUMNS),
        ("window_summary.tsv", window_rows, ["chrom", "start0", "end0", "threshold"] + RAW_COLUMNS + ["callable_bases", "multi_per_callable_mb"] + EXPECTED_COLUMNS),
        ("context_summary.tsv", context_rows, ["threshold", "context", "primary_change"] + RAW_COLUMNS),
        ("occupancy_thresholds.tsv", threshold_rows, ["threshold"] + RAW_COLUMNS + EXPECTED_COLUMNS),
        ("expected_context.tsv", expected_context_rows, ["chrom", "threshold", "context", "primary_change", "primary_af_bin", "an_bin", "observed_sites", "observed_multi", "training_level", "training_sites", "training_multi", "training_autosomes", "held_out_rate", "expected_multi", "unavailable_reason"]),
    ]
    paths = []
    for filename, rows, columns in specs:
        path = out_dir / filename
        _write_tsv(path, rows, columns)
        paths.append(str(path))
    if make_plots:
        try:
            plot_paths, plot_warnings = _make_plots(out_dir, chromosome_rows, window_rows, context_rows, threshold_rows)
            paths.extend(plot_paths)
            warnings.extend(plot_warnings)
        except Exception as exc:
            warnings.append(f"Optional figure rendering failed; tables remain valid: {type(exc).__name__}: {exc}")
    manifest = {
        "status": "complete", "inputs": [str(p) for p in files], "outputs": paths,
        "n_strata_rows": n_strata_rows, "n_windows": len({(r["chrom"], r["start0"]) for r in window_rows}),
        "thresholds": sorted(raw_threshold), "represented_autosomes": sorted(rates.autosomes, key=summary_chrom_key),
        "window_size": window_size, "min_training_sites": min_training_sites,
        "region": region, "callable_bed": str(callable_bed) if callable_bed is not None else None,
        "scan_scope_supplied": has_scopes,
        "expected_method": "leave-one-autosome-out conditional observation probability; exact context/change/AF-bin/AN-bin, then context/change, then threshold; (multi+0.5)/(sites+1)",
        "expected_interpretation": "Conditional on an eligible polymorphic SNV site. Not a mutation rate, retention estimate, significance test, or evidence of balancing selection.",
        "wilson_interpretation": "Descriptive binomial intervals only; neighboring sites are not assumed independent for hypothesis testing.",
        "zero_window_policy": "Explicit callable zero windows only within supplied scan scope, or represented chromosomes/explicit region when scan scope is absent.",
        "window_coordinates": "0-based half-open complete bins. Callable bases are clipped to scan scope and region.",
        "warnings": warnings,
    }
    manifest_path = out_dir / "summarize_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    manifest["manifest_path"] = str(manifest_path)
    return manifest


def main_summary():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shards", nargs="+")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--window-size", type=int, default=1_000_000)
    parser.add_argument("--min-training-sites", type=int, default=100)
    parser.add_argument("--callable-bed")
    parser.add_argument("--region")
    parser.add_argument("--thresholds", help="comma-separated AC thresholds; also defines rows for an empty scan")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    thresholds = [int(x) for x in args.thresholds.split(",")] if args.thresholds else None
    result = summarize(args.shards, args.out_dir, args.window_size, args.min_training_sites,
                       not args.no_plots, args.callable_bed, args.region, thresholds)
    print(json.dumps(result, indent=2))


# ========================================================================
# Cli
# ========================================================================

# Portable local entry point. No AlphaGenome API, key or network dependency.


import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import fcntl
import json
import os
from pathlib import Path
import shutil
import sys
import time

from allele_common import VERSION, code_hash, dump_json, fingerprint, stable_hash, tsv_writer
from allele_common import FastaContext, _backend, _chromosomes, _file_chrom, discover_inputs, read_header, validate_count_header, parse_region
from allele_common import make_callable

DEFAULT_VCF = "/mnt/f/annot/gnomad/sites"
DEFAULT_FASTA = "/mnt/f/gen/1hgp/38/GRCH38.fasta.gz"
DEFAULT_OUT = "/mnt/d/analysis/allele3/gnomad_local"


def parser():
    p = argparse.ArgumentParser(
        description="allele3 local: gnomAD SNV ALT occupancy, distribution and ancestry summaries.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("command", nargs="?", default="run",
                   choices=["run", "scan", "summarize", "check", "annotate", "make-callable"])
    source = p.add_mutually_exclusive_group()
    source.add_argument("--vcf", help="One sites VCF/BCF; --region is optional")
    source.add_argument("--dir-vcf", "--vcf-dir", dest="vcf_dir",
                        help="Directory of disjoint chromosome VCF files")
    p.add_argument("--chr", dest="chromosomes", help="Directory selection: autosomes, all, or 1,2,X")
    p.add_argument("--region", help="1-based inclusive CHR:START-END; requires --vcf")
    p.add_argument("--fasta", default=os.environ.get("ALLELE_FASTA", DEFAULT_FASTA))
    p.add_argument("--dir-out", "--out", dest="out",
                   default=os.environ.get("ALLELE_LOCAL_OUT", os.environ.get("ALLELE_OUT", DEFAULT_OUT)))
    p.add_argument("--dir-cache", "--cache-dir", dest="cache_dir",
                   default=os.environ.get("ALLELE_CACHE_DIR", ""), help="Default: OUT/shards")
    p.add_argument("--jobs", type=int, default=4, help="Concurrent chromosome scans; I/O often limits speed")
    p.add_argument("--backend", choices=["auto", "bcftools", "python"], default="auto")
    p.add_argument("--window-size", type=int, default=1000000)
    p.add_argument("--thresholds", default="1,2,5,10", help="ALT count thresholds; 1 is always included")
    p.add_argument("--min-an", type=int, default=100, help="Per-ALT AN floor; not a universal callability guarantee")
    p.add_argument("--min-af", type=float, default=0.0, help="Per-ALT AF floor (AC/AN); default retains singletons")
    p.add_argument("--populations", default="auto", help="auto, none, or ancestry suffixes e.g afr,eas,nfe")
    p.add_argument("--min-pop-an", type=int, default=100)
    p.add_argument("--min-pop-af", type=float, default=0.0,
                   help="Ancestry-specific AF sensitivity; unequal detection power still remains")
    p.add_argument("--include-unfiltered", action="store_true", help="Explicitly also accept FILTER=.; default PASS only")
    p.add_argument("--exclude-flags", default="", help="Optional INFO region flags: lcr,segdup")
    p.add_argument("--vep", action="store_true", help="Also read existing INFO/vep or CSQ; increases VCF transfer/CPU")
    p.add_argument("--write-all-sites", action="store_true", help="Optional all-site compact TSV; default only multi-ALT details")
    p.add_argument("--candidate-storage", choices=["deferred", "eager"], default="deferred",
                   help="Defer the redundant three-candidate-per-site table until explicit full offline annotate; sparse Atlas does not need it")
    p.add_argument("--callable-bed", help="Optional same-build BED mask; applied to both numerator and denominator")
    p.add_argument("--build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    p.add_argument("--source-label", default="gnomAD-v4.1-genomes")
    p.add_argument("--min-training-sites", type=int, default=100)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--scores", help="Local precomputed SNV scores for annotate only")
    p.add_argument("--score-column", default="phred")
    p.add_argument("--score-format", choices=["tsv", "avi"], default="tsv",
                   help="avi explicitly selects the official six-column headerless format")
    p.add_argument("--score-label", default="local")
    p.add_argument("--score-version", default="unspecified")
    p.add_argument("--score-build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    p.add_argument("--score-backend", choices=["auto", "stream", "tabix"], default="auto")
    p.add_argument("--all-sites-an", help="Local official overall AN TSV for make-callable")
    p.add_argument("--callable-min-an", type=int, default=100000)
    p.add_argument("--output-bed", help="New plain BED output for make-callable")
    p.add_argument("--with-architecture", action="store_true")
    p.add_argument("--enrichment-track", action="append", default=[], help="NAME=BED4; uses the original scanner QC")
    p.add_argument("--enrichment-bootstrap", type=int, default=2000)
    p.add_argument("--ld-block-dir", default="/mnt/f/LD/block")
    p.add_argument("--architecture-build", choices=["37", "38"], default="38")
    p.add_argument("--erv-joined", help="Joined all-repeat table; absent annotation is NA")
    p.add_argument("--architecture-prefix", default="architecture", help="Shard-local filename prefix")
    p.add_argument("--architecture-min-sites", type=int, default=10000)
    p.add_argument("--architecture-top-fraction", type=float, default=.01)
    p.add_argument("--architecture-annotation-cache")
    p.add_argument("--architecture-segdup-bed")
    p.add_argument("--architecture-lcr-bed")
    p.add_argument("--architecture-low-mappability-bed")
    p.add_argument("--version", action="version", version=VERSION)
    return p


def validated_options(a):
    if a.jobs < 1 or a.window_size < 1 or a.min_an < 1 or a.min_pop_an < 1 or a.min_training_sites < 1:
        raise ValueError("jobs, window size, AN floors, and minimum training sites must be positive")
    if not 0 <= a.min_af <= 1 or not 0 <= a.min_pop_af <= 1:
        raise ValueError("AF thresholds must be between 0 and 1")
    try:
        thresholds = sorted({1, *(int(t) for t in a.thresholds.split(","))})
    except ValueError as e:
        raise ValueError("--thresholds requires positive comma-separated integers") from e
    if thresholds[0] < 1:
        raise ValueError("AC thresholds must be positive")
    flags = [x for x in a.exclude_flags.split(",") if x]
    if set(flags) - {"lcr", "segdup"}:
        raise ValueError("--exclude-flags accepts lcr and segdup")
    if a.region:
        if not a.vcf:
            raise ValueError("--region requires one --vcf")
        parse_region(a.region)
    if a.vcf and a.chromosomes:
        raise ValueError("Use --region with one --vcf; --chr selects directory inputs")
    options = dict(fasta=str(Path(a.fasta).expanduser().resolve()), backend=a.backend,
                window_size=a.window_size, thresholds=thresholds, min_an=a.min_an,
                min_af=a.min_af, populations=a.populations, min_pop_an=a.min_pop_an,
                min_pop_af=a.min_pop_af, include_unfiltered=a.include_unfiltered,
                exclude_flags=flags, vep=a.vep, write_all_sites=a.write_all_sites,
                candidate_storage=a.candidate_storage,
                callable_bed=str(Path(a.callable_bed).resolve()) if a.callable_bed else None,
                build=a.build, region=a.region, source_label=a.source_label)
    if a.enrichment_track and not a.with_architecture: raise ValueError('--enrichment-track requires --with-architecture')
    if a.with_architecture:
        if a.enrichment_bootstrap < 2000: raise ValueError('Formal enrichment requires at least 2000 bootstrap replicates')
        tracks={}
        for spec in a.enrichment_track:
            name,sep,path=spec.partition('=')
            if not sep or not name or name in tracks: raise ValueError('Use unique NAME=BED enrichment tracks')
            tracks[name]=str(Path(path).resolve())
        options.update(enrichment_tracks=tracks,enrichment_bootstrap=a.enrichment_bootstrap)
        if a.build != "GRCh"+a.architecture_build:
            raise ValueError("Architecture annotation build must match the VCF/reference build")
        if a.erv_joined and a.architecture_build != "38":
            raise ValueError("Current joined RepeatMasker input is GRCh38; use an explicitly prepared matching-build track")
        if not a.architecture_prefix or Path(a.architecture_prefix).name!=a.architecture_prefix or a.architecture_prefix in {".",".."}:
            raise ValueError("--architecture-prefix must be a shard-local filename, not a path")
        if a.architecture_min_sites<1 or not 0<a.architecture_top_fraction<=1:
            raise ValueError("Invalid architecture ranking thresholds")
        options.update(with_architecture=True, architecture_schema=3,
            architecture_build=a.architecture_build, ld_block_dir=str(Path(a.ld_block_dir).resolve()),
            erv_joined=str(Path(a.erv_joined).resolve()) if a.erv_joined else None,
            architecture_prefix=a.architecture_prefix, architecture_min_sites=a.architecture_min_sites,
            architecture_top_fraction=a.architecture_top_fraction)
        for key in ('architecture_annotation_cache','architecture_segdup_bed','architecture_lcr_bed','architecture_low_mappability_bed'):
            value=getattr(a,key)
            options[key]=str(Path(value).resolve()) if value else None
        pass
        options['architecture_annotations']=annotation_identity(options)
    return options


def load_run_metadata(out):
    """Path-aware metadata only; never read source VCF/FASTA or scan tables."""
    from scan_store import load_run_metadata as reader
    return reader(out)


def resolve_run_path(run, path):
    from scan_store import resolve_run_path as resolve
    return resolve(run, path)


def verify_recorded_input(run, expected):
    from scan_store import verify_recorded_input as verify
    return verify(run, expected)


def load_run(out):
    data = load_run_metadata(out)
    for shard in data["shards"]:
        verify_shard(shard, data["signature"])
    return data


def verify_shard(shard, signature):
    path = Path(shard)
    marker = path / "COMPLETE.json"
    if not marker.is_file():
        raise ValueError(f"Missing shard completion marker: {marker}")
    complete = json.loads(marker.read_text())
    if complete["signature"] != signature:
        raise ValueError(f"Cache identity differs at {path}; use a separate output/cache directory")
    for name, metadata in complete["outputs"].items():
        f = path / name
        if not f.is_file() or f.stat().st_size != metadata["size"] or fingerprint(f)["edge_sha256"] != metadata["edge_sha256"]:
            raise ValueError(f"Completed cache file is missing or changed: {f}; restore it or use a new cache")
    return complete


def _scan_task(task):
    """Serialize identical shared-cache shards across distinct output folders."""
    target = Path(task[1])
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.with_name(target.name + ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _scan_task_locked(task)


def _scan_task_locked(task):
    source, target, options, expected, signature = task
    target = Path(target)
    if (target / "COMPLETE.json").exists():
        verify_shard(target, signature)
        print(f"[resume] {Path(source).name}: completed shard reused", flush=True)
        return str(target)
    if target.exists():
        raise ValueError(f"Unrecognized incomplete directory: {target}; use a new --dir-cache")
    staging = target.with_name(target.name + ".partial")
    # This precisely named staging directory is owned by this workflow.
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    try:
        result = scan_file(source, staging, options, expected)
        outputs = {p.name: fingerprint(p) for p in staging.iterdir() if p.is_file()}
        dump_json(staging / "COMPLETE.json", dict(signature=signature, outputs=outputs,
                                                 qc=result["qc"], elapsed_seconds=result["elapsed_seconds"]))
        os.replace(staging, target)
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    print(f"[done] {Path(source).name}: {result['qc']['eligible_sites']:,} eligible loci; "
          f"{result['qc']['multi_alt_sites']:,} multi-ALT", flush=True)
    return str(target)


def merge_small_reports(shards, out):
    for name in ("frequency_counts.tsv", "population_counts.tsv", "population_classes.tsv",
                 "annotation_counts.tsv", "relationship_counts.tsv"):
        counts, fields = Counter(), None
        for shard in shards:
            with (Path(shard) / name).open(newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if fields is None:
                    fields = reader.fieldnames
                elif fields != reader.fieldnames:
                    raise ValueError(f"Inconsistent {name} schemas")
                for row in reader:
                    counts[tuple(row[x] for x in fields[:-1])] += int(row["count"])
        handle, writer = tsv_writer(Path(out) / name, fields)
        with handle:
            for key, value in sorted(counts.items()):
                writer.writerow(dict(zip(fields, (*key, value))))
    pop = defaultdict(Counter)
    with (Path(out) / "population_counts.tsv").open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pop[(r["chrom"], r["threshold"], r["population"])][r["metric"]] += int(r["count"])
    cols = ["chrom", "threshold", "population", "n_evaluable", "n_sites", "n_multi", "n_three", "multi_alt_fraction"]
    handle, writer = tsv_writer(Path(out) / "population_summary.tsv", cols)
    with handle:
        for key, values in sorted(pop.items()):
            row = dict(zip(cols[:3], key))
            row.update({m: values[m] for m in cols[3:7]})
            row["multi_alt_fraction"] = values["n_multi"]/values["n_sites"] if values["n_sites"] else "."
            writer.writerow(row)
    pass
    return merge_site_spectra(shards, out)


def summarize_run(out, min_training_sites=100, make_plots=True):
    out = Path(out)
    run = load_run(out)
    if (out / "relocation.json").is_file():
        raise ValueError("Adopted scan stores preserve the original summary unchanged. Run allele3.sh for new analyses in a separate output directory.")
    options = run["options"]
    mask = options.get("callable_bed")
    if mask and fingerprint(mask) != run["identity"]["callable_bed"]:
        raise ValueError("Callable mask changed since scanning; use the original file or rescan to a new output")
    stage = out / "summary.partial"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir()
    try:
        result = summarize(run["shards"], stage, window_size=options["window_size"],
                           min_training_sites=min_training_sites, make_plots=make_plots,
                           callable_bed=mask, region=options.get("region"), thresholds=options["thresholds"])
        spectrum = merge_small_reports(run["shards"], stage)
        if options.get('with_architecture'):
            pass
            if annotation_identity(options)!=options.get('architecture_annotations'):
                raise ValueError('Architecture annotations changed since scanning')
            summarize_architecture(run['shards'],stage,options)
        qc = Counter()
        for shard in run["shards"]:
            record = json.loads((Path(shard) / "scan.json").read_text())
            qc.update(record["qc"])
        dump_json(stage / "qc_summary.json", dict(qc))
        report = f"""# allele3 local analysis

Version: {VERSION}

Source: {options['source_label']}; genome build: {options['build']}.

- Qualifying SNV positions: {qc['eligible_sites']:,}
- Positions with at least two qualifying ALT nucleotides: {qc['multi_alt_sites']:,}
- Positions with three qualifying ALT nucleotides: {qc['three_alt_loci']:,}
- Multi-ALT positions with different ALT AN values: {qc['multi_alt_loci_discordant_an']:,}
- AC thresholds: {options['thresholds']}
- Per-ALT AN floor: {options['min_an']}
- Callable mask: {mask or 'none; proportions are conditional on observed qualifying SNV positions'}
- Site AC spectrum coverage: {spectrum['status']} ({len(spectrum['measured_shards'])}/{spectrum['total_shards']} shards)
- Theoretical candidate storage: {options.get('candidate_storage', 'eager (legacy)')}

The primary quantity is observed non-reference nucleotide occupancy. REF is not
assumed ancestral, major, present, or reconstructable from split AN values.
ALT order is not a mutation timeline. Failed/missing ALT records are not used as
unlisted counterfactuals. Frequencies are each ALT's AC/AN.

The context background uses other scanned autosomes and reports support and
fallback levels. It estimates a conditional observation frequency, not mutation
rate or evolutionary retention. X/Y/MT and regional/single-chromosome scans are
descriptive. All ancestry comparisons remain sensitive to sample size, coverage,
ancestry composition and genotype filtering.

No AlphaGenome network requests were made by this pipeline. Optional offline AVI
annotation is exploratory: AVI used gnomAD frequency proxy labels in training.
This workflow cannot infer haplotype ancestry, allele age, pleiotropy, or balancing
selection from sites-only counts.

Detailed multi-ALT sites and original per-ALT counts are kept in the per-chromosome
shards. With deferred candidate storage, the redundant three-substitution table
is reconstructed losslessly only if full offline annotate is requested. Sparse
Atlas uses its separate bounded plan. Unlisted candidate AC/AN/AF remain missing.

site_ac_spectrum.tsv and singleton_sensitivity.tsv distinguish multi-ALT sites
containing AC=1 ALTs from sites still having two ALTs after AC>=2. The denominator
is the original eligible multi-ALT cohort. This is a count sensitivity, not loss
over evolutionary time. Old caches lacking this spectrum are explicitly unavailable;
mixed cache summaries are partial and cannot be treated as genome-wide zero values.
"""
        (stage / "REPORT.md").write_text(report, encoding="utf-8")
        target, backup = out / "summary", out / "summary.previous"
        # The summarizer writes absolute paths before its directory is committed.
        # Rebase metadata while it is still private, then publish the whole set.
        def rebase(value):
            if isinstance(value, str) and value.startswith(str(stage.resolve()) + os.sep):
                return str(target.resolve()) + value[len(str(stage.resolve())):]
            if isinstance(value, list):
                return [rebase(x) for x in value]
            if isinstance(value, dict):
                return {k: rebase(v) for k, v in value.items()}
            return value
        result = rebase(result)
        for metadata in stage.glob("*.json"):
            dump_json(metadata, rebase(json.loads(metadata.read_text())))
        if backup.exists():
            shutil.rmtree(backup)
        if target.exists():
            os.replace(target, backup)
        os.replace(stage, target)
        if backup.exists():
            shutil.rmtree(backup)
        run["status"] = "complete"
        run["summary"] = str(target.resolve())
        run["summary_parameters"] = dict(min_training_sites=min_training_sites, make_plots=make_plots)
        run.pop("_store", None)  # Runtime path state is not part of a scientific manifest.
        dump_json(out / "run.json", run)
    except BaseException:
        if stage.exists():
            shutil.rmtree(stage)
        raise
    print(f"[complete] {out / 'summary'}", flush=True)
    return result


def run_scan(a, options):
    out = Path(a.out).expanduser().resolve()
    if (out / "relocation.json").is_file():
        raise ValueError("Adopted scan stores are read-only. Use 3.scan.sh verify and allele3.sh; use a new output for a new scan.")
    out.mkdir(parents=True, exist_ok=True)
    directory = None if a.vcf else (a.vcf_dir or os.environ.get("ALLELE_VCF_DIR", os.environ.get("ALLELE_GNOMAD", DEFAULT_VCF)))
    chromosomes = None if a.vcf else (a.chromosomes or "autosomes")
    paths = discover_inputs(a.vcf, directory, chromosomes)
    with FastaContext(options["fasta"]):
        pass
    if options["build"] != "GRCh38" and "gnomad-v4" in options["source_label"].lower():
        raise ValueError("gnomAD v4 uses GRCh38; supply a correct build/source label")
    identity = dict(version=VERSION, code_sha256=code_hash(), options=options,
                    inputs=[fingerprint(p) for p in paths], fasta=fingerprint(options["fasta"]),
                    fasta_index=fingerprint(options["fasta"] + ".fai"),
                    callable_bed=fingerprint(options["callable_bed"]) if options["callable_bed"] else None)
    if Path(options["fasta"] + ".gzi").exists():
        identity["fasta_gzi"] = fingerprint(options["fasta"] + ".gzi")
    signature = stable_hash(identity)
    manifest = out / "run.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if old["signature"] != signature:
            raise ValueError("Inputs, scan options or code differ from this output. Use a new --dir-out; old results are preserved.")
    elif any(out.iterdir()):
        unexpected = [p.name for p in out.iterdir() if p.name != ".allele3.lock" and not (p.name == "logs" and p.is_dir())]
        if unexpected:
            raise ValueError(f"Output is not an empty local-workflow directory: {out}; choose a new --dir-out")
    cache = Path(a.cache_dir).expanduser().resolve() if a.cache_dir else out / "shards"
    cache.mkdir(parents=True, exist_ok=True)
    tasks, shards = [], []
    for i, path in enumerate(paths):
        expected = _file_chrom(path)
        if options.get("region"):
            expected = parse_region(options["region"])[0]
        label = f"chr{expected}" if expected else f"input{i+1}"
        target = cache / (label + "." + signature[:16])
        shards.append(str(target))
        tasks.append((str(path), str(target), options, expected, signature))
    run = dict(status="running", signature=signature, identity=identity, options=options,
               shards=shards, version=VERSION, cache=str(cache), inputs=[str(p.resolve()) for p in paths])
    dump_json(manifest, run)
    if options.get('with_architecture') and options.get('architecture_annotation_cache') and options.get('erv_joined'):
        pass
        prepare_joined_cache(options['erv_joined'],options['architecture_annotation_cache'])
    started = time.monotonic()
    print(f"[scan] {len(paths)} input(s); {min(a.jobs, len(paths))} concurrent scan(s); "
          f"backend={_backend(paths[0], a.backend)}; no API calls", flush=True)
    if a.jobs == 1 or len(tasks) == 1:
        for task in tasks:
            _scan_task(task)
    else:
        with ProcessPoolExecutor(max_workers=min(a.jobs, len(tasks))) as pool:
            futures = [pool.submit(_scan_task, task) for task in tasks]
            for future in as_completed(futures):
                future.result()
    run["status"] = "scanned"
    run["scan_elapsed_seconds_this_invocation"] = time.monotonic() - started
    dump_json(manifest, run)
    return run


def main_modern(argv=None):
    a = parser().parse_args(argv)
    if a.command == "make-callable":
        if not a.all_sites_an or not a.output_bed:
            raise ValueError("make-callable requires --all-sites-an and --output-bed")
        if str(a.output_bed).endswith((".gz", ".bgz")):
            raise ValueError("--output-bed must be a plain .bed path")
        data = make_callable(a.all_sites_an, a.output_bed, a.callable_min_an, _chromosomes(a.chromosomes))
        print(json.dumps(data, ensure_ascii=False, indent=2))
        return 0
    if a.command == "check":
        options = validated_options(a)
        directory = None if a.vcf else (a.vcf_dir or os.environ.get("ALLELE_VCF_DIR", DEFAULT_VCF))
        paths = discover_inputs(a.vcf, directory, None if a.vcf else (a.chromosomes or "autosomes"))
        with FastaContext(options["fasta"]) as fa:
            print(f"FASTA: {len(fa.lengths)} indexed contigs; build={a.build}")
        for path in paths:
            h = read_header(path, a.backend)
            validate_count_header(h)
            print(f"{path.name}: backend={_backend(path,a.backend)}; "
                  f"INFO AC/AN valid; samples_in_header={h['samples']}; "
                  f"ancestries={','.join(discover_populations(h,a.populations)) or 'none'}")
        print("Ready. No API key is read. Header sample count 0 is expected for sites-only files.")
        return 0
    out = Path(a.out).expanduser().resolve()
    if (out / "relocation.json").is_file():
        raise ValueError("Adopted scan stores are read-only. Use 3.scan.sh verify and allele3.sh; use a new output for a new scan.")
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".allele3.lock").open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as e:
            raise ValueError(f"Another allele3 process is using {out}") from e
        if a.command in ("run", "scan"):
            run_scan(a, validated_options(a))
            if a.command == "run":
                summarize_run(out, a.min_training_sites, not a.no_plots)
        elif a.command == "summarize":
            summarize_run(out, a.min_training_sites, not a.no_plots)
        elif a.command == "annotate":
            if not a.scores:
                raise ValueError("annotate requires --scores pointing to an existing local file")
            run = load_run(out)
            if a.score_build != run["options"]["build"]:
                raise ValueError("Score genome build must match the scanned VCF/FASTA; no automatic liftover")
            pass
            pass
            preflight_scores(a.scores, a.score_column, a.score_build, a.score_backend, a.score_format)
            candidate_cache = candidate_files_for_run(run, out / "candidate_cache")
            result = annotate_scores(
                candidate_cache["files"],
                a.scores, out / "scores", a.score_column, a.score_build,
                a.score_label, a.score_version, a.score_backend, 1000000, a.score_format)
            print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


# ========================================================================
# scan_architecture
# ========================================================================

# LD blocks, repeat annotation, architecture, UKB and analysis universes.





# ========================================================================
# Ldblocks
# ========================================================================

from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
from allele_common import canonical_chrom, chrom_sort_key

@dataclass(frozen=True)
class Block:
    chrom: str
    start: int
    end: int
    ancestry: str
    build: str
    index: int

    @property
    def block_id(self):
        return f"{self.ancestry}.{self.build}.{self.chrom}.{self.index:04d}"

def read_blocks(path, ancestry, build):
    if str(build) not in {"37", "38"}:
        raise ValueError("LD build must be 37 or 38")
    rows=[]
    with open(path, "rt", encoding="utf-8") as f:
        for line_no,line in enumerate(f,1):
            if not line.strip() or line.startswith("#"): continue
            p=line.rstrip("\n").split()
            if len(p)<3:
                raise ValueError(f"{path}:{line_no}: expected >=3 columns")
            if p[:3] in (["chr","start","end"],["chrom","start","end"],["chromosome","start","end"]): continue
            chrom=canonical_chrom(p[0]); start=int(p[1]); end=int(p[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: invalid BED interval")
            rows.append((chrom,start,end))
    rows.sort(key=lambda x:(chrom_sort_key(x[0]),x[1],x[2]))
    out=[]
    per=defaultdict(int)
    prev={}
    stats={"n":0,"overlaps":0,"touching":0,"gaps":0}
    for chrom,start,end in rows:
        per[chrom]+=1
        b=Block(chrom,start,end,ancestry,str(build),per[chrom])
        if chrom in prev:
            pe=prev[chrom]
            if start < pe: raise ValueError(f"Overlapping BED intervals in {path}: {chrom}:{start}-{end}")
            elif start == pe: stats["touching"]+=1
            else: stats["gaps"]+=1
        prev[chrom]=max(prev.get(chrom,0),end)
        out.append(b)
        stats["n"]+=1
    return out,stats

class BlockSweep:
    """O(1)-amortized assignment for sites supplied in sorted chromosome/position order."""
    def __init__(self, blocks):
        self.by_chrom=defaultdict(list)
        for b in blocks: self.by_chrom[b.chrom].append(b)
        self._chrom=None
        self._i=0
        self._last={}

    def locate(self, chrom, pos1):
        chrom=canonical_chrom(chrom)
        pos0=int(pos1)-1
        if pos0 < 0 or pos0 < self._last.get(chrom, -1):
            raise ValueError("LD sweep requires nondecreasing positive positions per chromosome")
        self._last[chrom]=pos0
        arr=self.by_chrom.get(chrom,[])
        if chrom != self._chrom:
            self._chrom=chrom; self._i=0
        while self._i < len(arr) and arr[self._i].end <= pos0:
            self._i += 1
        if self._i < len(arr):
            b=arr[self._i]
            if b.start <= pos0 < b.end:
                return b
        return None

def validate_directory(directory, build="38", ancestries=("AFR","EAS","EUR","SAS")):
    d=Path(directory)
    result={}
    for anc in ancestries:
        p=d/f"{anc}.{build}.bed"
        if not p.exists():
            result[anc]={"path":str(p),"missing":True}
            continue
        blocks,stats=read_blocks(p,anc,build)
        result[anc]={"path":str(p),"missing":False,**stats}
    return result


def load_block_sets(directory, build="38"):
    result={}
    for anc in ("AFR","EAS","EUR","SAS"):
        path=Path(directory)/f"{anc}.{build}.bed"
        blocks,_=read_blocks(path,anc,build)
        if not blocks: raise ValueError(f"Empty LD partition: {path}")
        result[anc]=blocks
    return result

class BlockIndex:
    """Random-access assignment for candidate tables which are not coordinate sorted."""
    def __init__(self, blocks):
        from bisect import bisect_right
        self.bisect=bisect_right
        self.rows=defaultdict(list)
        for b in blocks: self.rows[b.chrom].append(b)
        for c in self.rows: self.rows[c].sort(key=lambda b:b.start)
        self.starts={c:[b.start for b in a] for c,a in self.rows.items()}
    def locate(self, chrom, pos1):
        c=canonical_chrom(chrom); p=int(pos1)-1
        i=self.bisect(self.starts.get(c,[]),p)-1
        if i>=0 and p<self.rows[c][i].end: return self.rows[c][i]
        return None

def hotspot_consensus(rows, min_sites=10000, top_fraction=.01, threshold=1):
    """Intersect independently ranked partitions; never equate their block IDs."""
    import math
    groups=defaultdict(list)
    for r in rows:
        if r.get("sensitivity","all") != "all" or str(r.get("threshold",1)) != str(threshold) or r.get("population_mode","all") != "all": continue
        if r["block_id"]=="outside_partition" or int(r["n_sites"])<min_sites or r.get("oe") in (None,"NA"): continue
        if not math.isfinite(float(r["oe"])): continue
        groups[r["ancestry"]].append(r)
    events=defaultdict(list)
    for anc, group in groups.items():
        group.sort(key=lambda r:(-float(r['oe']),r['chrom'],int(r['start'])))
        for r in group[:max(1,math.ceil(len(group)*top_fraction))]:
            events[r['chrom']].extend([(int(r['start']),1,anc,r['block_id']),(int(r['end']),-1,anc,r['block_id'])])
    result=[]
    for chrom, ev in sorted(events.items(),key=lambda x:chrom_sort_key(x[0])):
        active={}; previous=None
        import itertools
        for pos, at in itertools.groupby(sorted(ev),key=lambda e:e[0]):
            if previous is not None and pos>previous and len(active)>=2:
                result.append(dict(chrom=chrom,start=previous,end=pos,overlap_bp=pos-previous,
                    supporting_ancestries=','.join(sorted(active)),n_ancestries=len(active),
                    block_ids=','.join(active[a] for a in sorted(active))))
            for _,direction,anc,bid in at:
                if direction<0: active.pop(anc,None)
                else: active[anc]=bid
            previous=pos
    return result


def consensus_regions(fragments, required_ancestries=4):
    """Merge adjacent four-way fragments without hard-coded genomic coordinates."""
    result=[]
    for row in sorted(fragments,key=lambda r:(chrom_sort_key(r['chrom']),int(r['start']))):
        if int(row['n_ancestries']) != required_ancestries: continue
        r=dict(row); r['start']=int(r['start']); r['end']=int(r['end'])
        if result and result[-1]['chrom']==r['chrom'] and r['start']<=result[-1]['end']:
            old=result[-1]; old['end']=max(old['end'],r['end']); old['overlap_bp']=old['end']-old['start']
            old['block_ids']=','.join(sorted(set(old['block_ids'].split(','))|set(r['block_ids'].split(','))))
        else: result.append(r)
    return result


# ========================================================================
# Erv
# ========================================================================

# UCSC joined repeats: aligned bodies, strand-aware observed element edges.

from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path
import gzip, json, os, fcntl
from allele_common import canonical_chrom, open_text as architecture_open_text, chrom_sort_key
from allele_common import fingerprint, dump_json
CORE_FAMILIES={'ERV1','ERVK','ERVL'}
RELATED_FAMILIES=CORE_FAMILIES | {'ERVL-MaLR','Gypsy'}
PRIMARY={str(i) for i in range(1,23)}|{'X','Y'}

@dataclass(frozen=True, slots=True)
class ERV:
    chrom:str
    start:int
    end:int
    name:str
    rep_class:str
    family:str
    strand:str
    source_id:str=''
    segments:tuple=()
    @property
    def tier(self):
        if self.rep_class!='LTR': return 'non_ltr'
        if self.family in CORE_FAMILIES: return 'erv_core'
        if self.family in RELATED_FAMILIES: return 'erv_related'
        return 'ltr_other'
    @property
    def five_prime(self): return self.start if self.strand=='+' else self.end
    @property
    def three_prime(self): return self.end if self.strand=='+' else self.start

def parse_joined_rmsk(path, include_alt_contigs=False, chromosome=None):
    """Aligned blocks define observed joined edges, excluding consensus tails.

    Only blocks whose relative start is not the -1 sentinel are aligned DNA; nested-insertion gaps do
    not count toward repeat bodies. Retain all classes for negative controls.
    """
    out=[]
    with architecture_open_text(path) as f:
        for ln,line in enumerate(f,1):
            if not line.strip() or line.startswith('#'): continue
            a=line.rstrip('\r\n').split('\t') if '\t' in line else line.split()
            if len(a)!=14: raise ValueError(f'{path}:{ln}: expected 14 joined RepeatMasker columns')
            chrom=canonical_chrom(a[1])
            if chromosome and chrom!=canonical_chrom(chromosome): continue
            if not include_alt_contigs and chrom not in PRIMARY: continue
            outer_start,outer_end=int(a[2]),int(a[3])
            start,end=int(a[7]),int(a[8]); strand=a[6]
            if not (0<=outer_start<outer_end and 0<=start<end) or strand not in {'+','-'}:
                raise ValueError(f'{path}:{ln}: invalid joined coordinates/strand')
            sizes=[int(x) for x in a[11].rstrip(',').split(',')]
            offsets=[int(x) for x in a[12].rstrip(',').split(',')]
            if len(sizes)!=int(a[10]) or len(offsets)!=len(sizes): raise ValueError(f'{path}:{ln}: invalid joined block cardinality')
            segments=[]
            for size,offset in zip(sizes,offsets):
                if offset != -1 and size<0: raise ValueError('Invalid aligned joined block')
                # Unaligned consensus gaps can have negative sizes (overlap); ignore them.
                if offset!=-1 and size>0:
                    lo=outer_start+offset; hi=lo+size
                    if lo<0: raise ValueError(f'{path}:{ln}: negative aligned genomic coordinate')
                    segments.append((lo,hi))
            if not segments: raise ValueError(f'{path}:{ln}: no aligned repeat sequence')
            # A few UCSC records have overlapping, out-of-order blocks: alignEnd
            # is the final block's end, not the maximum observed aligned end.
            # Use the envelope of ALL aligned blocks, preserving every body base.
            start=min(s for s,t in segments); end=max(t for s,t in segments)
            name,sep,cf=a[4].rpartition('#')
            if not sep: raise ValueError(f'{path}:{ln}: repeat class missing from joined name')
            rep_class,_,family=cf.partition('/')
            out.append(ERV(chrom,start,end,name,rep_class,family,strand,a[13],tuple(segments)))
    out.sort(key=lambda e:(chrom_sort_key(e.chrom),e.start,e.end,e.name))
    return out

def write_erv_bed(elements,out_path):
    Path(out_path).parent.mkdir(parents=True,exist_ok=True)
    with open(out_path,'w') as f:
        f.write('chrom\tstart\tend\tname\trep_class\tfamily\tstrand\ttier\tsource_id\tsegments\n')
        for e in elements:
            f.write('\t'.join(map(str,('chr'+e.chrom,e.start,e.end,e.name,e.rep_class,e.family,e.strand,e.tier,e.source_id,','.join(f'{s}-{t}' for s,t in e.segments))))+'\n')

class ERVOverlapSweep:
    def __init__(self,elements):
        self.by_chrom=defaultdict(list)
        for e in elements: self.by_chrom[e.chrom].append(e)
        for a in self.by_chrom.values(): a.sort(key=lambda e:e.start)
        self.chrom=None; self.i=0; self.active=[]; self.last={}
    def locate(self,chrom,pos1):
        chrom=canonical_chrom(chrom); pos0=int(pos1)-1
        if pos0<self.last.get(chrom,-1): raise ValueError('Repeat sweep requires sorted positions')
        self.last[chrom]=pos0
        arr=self.by_chrom.get(chrom,[])
        if chrom!=self.chrom: self.chrom=chrom; self.i=0; self.active=[]
        self.active=[e for e in self.active if e.end>pos0]
        while self.i<len(arr) and arr[self.i].start<=pos0:
            e=arr[self.i]
            if e.end>pos0: self.active.append(e)
            self.i+=1
        return tuple(e for e in self.active if any(lo<=pos0<hi for lo,hi in (e.segments or ((e.start,e.end),))))

class BoundarySweep:
    def __init__(self,elements,tiers=('erv_core','erv_related')):
        self.points=defaultdict(lambda:{'5p':[],'3p':[]}); allowed=set(tiers)
        for e in elements:
            if e.tier not in allowed: continue
            self.points[e.chrom]['5p'].append((e.five_prime,e))
            self.points[e.chrom]['3p'].append((e.three_prime,e))
        for v in self.points.values():
            for a in v.values(): a.sort(key=lambda x:(x[0],x[1].name,x[1].source_id))
        self.state={}; self.last={}
    def nearest(self,chrom,pos1,boundary_type,max_distance=10000):
        chrom=canonical_chrom(chrom); p=int(pos1)-1; key=(chrom,boundary_type)
        if boundary_type not in {'5p','3p'}: raise ValueError('Boundary must be 5p or 3p')
        if p<self.last.get(key,-1): raise ValueError('Boundary sweep requires sorted sites')
        self.last[key]=p
        a=self.points.get(chrom,{}).get(boundary_type,[])
        if not a: return None
        i=self.state.get(key,0)
        while i+1<len(a) and a[i+1][0]<=p: i+=1
        self.state[key]=i
        point,e=min(a[i:i+2],key=lambda x:(abs(p-x[0]),x[0],x[1].name))
        distance=abs(p-point)
        if distance>=max_distance: return None
        signed=(p-point)*(1 if e.strand=='+' else -1)
        return dict(distance=distance,signed_genomic_distance=p-point,signed_transcription_distance=signed,boundary_pos0=point,element=e)

def boundary_bin(distance):
    if distance is None: return 'gt10kb'
    lo=0
    for hi in (100,200,500,1000,2000,5000,10000):
        if abs(int(distance))<hi: return f'{lo}-{hi}'
        lo=hi
    return 'ge10000'

def validate_existing_joined_cache(path, directory):
    """Read-only validation; never creates a directory, lock or cache file."""
    from scan_store import verify_cached_fingerprint
    d = Path(directory)
    old = json.loads((d / 'manifest.json').read_text())
    verify_cached_fingerprint(path, old['source'])
    if not old.get('files'):
        raise ValueError('Completed annotation cache has no recorded outputs')
    for name, fp in old['files'].items():
        if Path(name).name != name:
            raise ValueError('Invalid annotation cache output name')
        verify_cached_fingerprint(d / name, fp)
    return old


def prepare_joined_cache(path, directory):
    """Split once, under a lock, so chromosome workers do not reread 5M repeats."""
    d=Path(directory)
    if (d/'manifest.json').exists():
        return validate_existing_joined_cache(path, d)
    d.mkdir(parents=True,exist_ok=True)
    with (d/'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        source=fingerprint(path); marker=d/'manifest.json'
        if marker.exists():
            old=json.loads(marker.read_text())
            from scan_store import verify_cached_fingerprint
            verify_cached_fingerprint(path, old['source'])
            for name,fp in old['files'].items():
                verify_cached_fingerprint(d/name, fp)
            return old
        handles={}; counts=defaultdict(int)
        try:
            with architecture_open_text(path) as f:
                for line in f:
                    if not line.strip() or line.startswith('#'): continue
                    a=line.split('\t',2)
                    if len(a)<3: raise ValueError('Malformed joined annotation')
                    c=canonical_chrom(a[1])
                    if c not in PRIMARY: continue
                    if c not in handles: handles[c]=gzip.open(d/f'chr{c}.tsv.gz.partial','wt',compresslevel=2)
                    handles[c].write(line); counts[c]+=1
        finally:
            for f in handles.values(): f.close()
        for c in PRIMARY:
            if c not in handles:
                with gzip.open(d/f'chr{c}.tsv.gz.partial','wt'): pass
            os.replace(d/f'chr{c}.tsv.gz.partial',d/f'chr{c}.tsv.gz')
        result=dict(schema=1,source=source,counts=dict(counts),files={p.name:fingerprint(p) for p in d.glob('chr*.tsv.gz')})
        dump_json(marker,result); return result


# ========================================================================
# Architecture
# ========================================================================

# Sidecar counters; original scanner alone defines eligible loci and ALT QC.

from collections import defaultdict
from pathlib import Path
import csv, gzip, math, json, shutil
from allele_common import canonical_chrom, af_bin as architecture_af_bin, an_bin as architecture_an_bin, open_text as architecture_open_text
from allele_common import dump_json, fingerprint

SCHEMA=3
STRATA=['chrom','context','af_bin','an_bin','sensitivity','population_mode']
COUNTS=['threshold','n_sites','n_multi','n_three']
HEADERS={'global':STRATA+COUNTS,'ld':['ancestry','block_id','start','end']+STRATA+COUNTS,
         'erv':['label']+STRATA+COUNTS,
         'exposure':['label','block_id']+STRATA+COUNTS}

class ArchitectureAccumulator:
    def __init__(self,ld_block_sets=None,erv_elements=None,*,thresholds=(1,),out_prefix=None,
                 masks=None,flag_annotations=(),tracks=None,track_only=False):
        self.track_only = track_only
        self.thresholds=tuple(thresholds)
        self.width=1+2*len(self.thresholds)
        self.ld={a:BlockSweep(b) for a,b in (ld_block_sets or {}).items()}
        self.erv_available=erv_elements is not None
        elements=erv_elements or []
        self.erv_sweep=ERVOverlapSweep(elements)
        self.boundaries={} if track_only else {
            'erv_core':BoundarySweep(elements,('erv_core',)),
            'erv_related':BoundarySweep(elements,('erv_core','erv_related')),
            'ltr_all':BoundarySweep(elements,('erv_core','erv_related','ltr_other'))}
        self.global_strata=defaultdict(self.zeros)
        self.ld_strata=defaultdict(self.zeros)
        self.erv_strata=defaultdict(self.zeros)
        self.ld_current={}; self.ld_pending={a:defaultdict(self.zeros) for a in self.ld}
        self.out_prefix=out_prefix; self.ld_handle=None
        self.exposure_handle=None; self.exposure_pending=defaultdict(self.zeros); self.exposure_current=None
        pass
        self.track_sweeps={name:Sweep(intervals) for name,intervals in (tracks or {}).items()}
        self.masks=masks or {}; self.flag_annotations=set(flag_annotations)
        self.sensitivities={'all'}
        self.classes=sorted({e.rep_class for e in elements})
        if out_prefix:
            Path(out_prefix).parent.mkdir(parents=True,exist_ok=True)
            if not track_only:
                self.ld_handle=gzip.open(str(out_prefix)+'.ld_strata.tsv.gz','wt',compresslevel=2,newline='')
                self.ld_writer=csv.writer(self.ld_handle,delimiter='\t'); self.ld_writer.writerow(HEADERS['ld'])
            self.exposure_handle=gzip.open(str(out_prefix)+'.exposure_strata.tsv.gz','wt',compresslevel=2,newline='')
            self.exposure_writer=csv.writer(self.exposure_handle,delimiter='\t'); self.exposure_writer.writerow(HEADERS['exposure'])
    def zeros(self): return [0]*self.width
    @staticmethod
    def stratum(context,af,an): return (str(context or 'NA'),architecture_af_bin(af),architecture_an_bin(an))
    def _add(self,table,key,values):
        z=table[key]
        for i,v in enumerate(values): z[i]+=v
    def _rows(self,data):
        for k,v in sorted(data.items()):
            for i,t in enumerate(self.thresholds): yield [*k,t,v[0],v[1+2*i],v[2+2*i]]
    def _flush_ld(self,anc):
        self.ld_writer.writerows(self._rows(self.ld_pending[anc])); self.ld_pending[anc].clear()
    def _flush_exposure(self):
        if self.exposure_handle:
            self.exposure_writer.writerows(row for row in self._rows(self.exposure_pending) if row[-4] in (5,10))
            self.exposure_pending.clear()
    def close(self):
        if self.exposure_handle: self.exposure_handle.close(); self.exposure_handle=None
        if self.ld_handle: self.ld_handle.close(); self.ld_handle=None
    def observe(self,chrom,pos1,n_alt,context,primary_alt_af,primary_alt_an,*,
                threshold_counts=None,population_mode='NA',region_flags=()):
        c=canonical_chrom(chrom); n=int(n_alt)
        counts=threshold_counts or {1:n}
        if any(t not in counts for t in self.thresholds): raise ValueError('Persistence counts must come from the original scanner')
        values=[1]+[v for t in self.thresholds for v in (int(counts[t]>=2),int(counts[t]>=3))]
        st=self.stratum(context,primary_alt_af,primary_alt_an)
        overlaps=self.erv_sweep.locate(c,pos1) if self.erv_available else ()
        sensitivities=['all']; flags=set(region_flags)
        exclusions={}
        for flag in ('lcr','segdup'):
            if flag in self.masks: exclusions[flag]=self.masks[flag].contains(c,pos1)
            elif flag in self.flag_annotations: exclusions[flag]=flag in flags
        if all(x in exclusions for x in ('lcr','segdup')) and not (exclusions['lcr'] or exclusions['segdup']):
            sensitivities.append('exclude_lcr_segdup')
        if self.erv_available:
            exclusions['simple_repeat']=any(e.rep_class in {'Simple_repeat','Low_complexity'} for e in overlaps)
            if not exclusions['simple_repeat']: sensitivities.append('exclude_simple_low_complexity')
        if 'low_mappability' in self.masks:
            exclusions['low_mappability']=self.masks['low_mappability'].contains(c,pos1)
            if not exclusions['low_mappability']: sensitivities.append('exclude_low_mappability')
        if len(exclusions)==4 and not any(exclusions.values()): sensitivities.append('exclude_combined')
        self.sensitivities.update(sensitivities)
        modes=['all']
        if population_mode!='NA' and n>=2: modes.append(population_mode)
        labels=[]
        if self.erv_available and not self.track_only:
            ltr=[e for e in overlaps if e.rep_class=='LTR']
            labels=['all_sites','inside_any_ltr' if ltr else 'outside_ltr']
            if any(e.tier=='erv_core' for e in ltr): labels.append('inside_erv_core')
            if any(e.tier in {'erv_core','erv_related'} for e in ltr): labels.append('inside_erv_related')
            labels += ['family:'+fam for fam in sorted({e.family or 'unspecified' for e in ltr})]
            labels += ['negative_control:'+cl for cl in sorted({e.rep_class for e in overlaps}&{'LINE','SINE','DNA','Simple_repeat','Low_complexity'})]
            for tier,sweep in self.boundaries.items():
                for bt in ('5p','3p'):
                    near=sweep.nearest(c,pos1,bt)
                    if near:
                        label=f'{tier}:{bt}_boundary:{boundary_bin(near["distance"])}'
                        labels.append(label)
                        labels.append(label+(':upstream' if near['signed_transcription_distance']<0 else ':downstream'))
        for name,sweep in self.track_sweeps.items():
            labels.extend(name+'::'+x.label for x in sweep.locate(c,pos1))
        labels=sorted(set(labels))
        locations=[]; eur_block='outside_partition' 
        for anc,sweep in self.ld.items():
            b=sweep.locate(c,pos1)
            key=(anc,b.block_id,b.start,b.end) if b else (anc,'outside_partition',-1,-1)
            if anc=='EUR': eur_block=key[1]
            if self.ld_handle:
                if self.ld_current.get(anc)!=(c,key): self._flush_ld(anc); self.ld_current[anc]=(c,key)
                table=self.ld_pending[anc]
            else: table=self.ld_strata
            locations.append((table,key))
        if self.exposure_current != (c,eur_block):
            self._flush_exposure(); self.exposure_current=(c,eur_block)
        for sens in sensitivities:
            for mode in modes:
                key=(c,*st,sens,mode)
                if not self.track_only:
                    self._add(self.global_strata,key,values)
                    for table,prefix in locations: self._add(table,(*prefix,*key),values)
                for label in labels:
                    if not self.track_only: self._add(self.erv_strata,(label,*key),values)
                    # Formal inference uses the fixed eligible-locus denominator.
                    if mode=='all': self._add(self.exposure_pending,(label,eur_block,*key),values)
    def write(self,out_prefix):
        Path(out_prefix).parent.mkdir(parents=True,exist_ok=True)
        if self.exposure_handle:
            self._flush_exposure()
        else:
            with gzip.open(str(out_prefix)+'.exposure_strata.tsv.gz','wt',compresslevel=2,newline='') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(HEADERS['exposure'])
                w.writerows(row for row in self._rows(self.exposure_pending) if row[-4] in (5,10))
        if self.track_only:
            self.close()
            return
        for kind,data in [('global',self.global_strata),('ld',self.ld_strata),('erv',self.erv_strata)]:
            if kind=='ld' and self.ld_handle:
                for anc in self.ld: self._flush_ld(anc)
                self.close(); continue
            with gzip.open(str(out_prefix)+f'.{kind}_strata.tsv.gz','wt',compresslevel=2,newline='') as f:
                w=csv.writer(f,delimiter='\t'); w.writerow(HEADERS[kind]); w.writerows(self._rows(data))
        dump_json(str(out_prefix)+'.availability.json',dict(schema=SCHEMA,
            erv_status='measured' if self.erv_available else 'NA',repeat_classes=self.classes,
            ld_ancestries=sorted(self.ld),sensitivity_analyses=sorted(self.sensitivities),
            missing_sensitivity_annotations=sorted({'lcr','segdup','low_mappability'}-self.flag_annotations-set(self.masks)),
            denominator='Fixed original eligible loci at every AC tier; population profiles conditional on base multi-ALT population class',
            boundaries='Observed joined alignStart/alignEnd, BED edges, absolute distance <10000; signed distance follows repeat strand',
            body='Union of aligned joined blocks, excluding unaligned tails and nested insertion gaps'))

def _read_tsv(path):
    with architecture_open_text(path) as f: yield from csv.DictReader(f,delimiter='\t')

def annotation_identity(options):
    if not options.get('with_architecture'): return None
    files={}
    directory=options.get('ld_block_dir')
    if directory:
        load_block_sets(directory,options['architecture_build'])
        files.update({anc:fingerprint(Path(directory)/f'{anc}.{options["architecture_build"]}.bed') for anc in ('AFR','EAS','EUR','SAS')})
    for name,path in options.get('enrichment_tracks',{}).items(): files['track:'+name]=fingerprint(path)
    for key in ('erv_joined','architecture_segdup_bed','architecture_lcr_bed','architecture_low_mappability_bed'):
        if options.get(key): files[key]=fingerprint(options[key])
    return dict(schema=SCHEMA,files=files)

def from_options(options,chrom,shard_dir,header):
    if not options.get('with_architecture'): return None
    blocks=load_block_sets(options['ld_block_dir'],options['architecture_build']) if options.get('ld_block_dir') else None
    elements=None
    if options.get('erv_joined'):
        source=options['erv_joined']
        if options.get('architecture_annotation_cache') and chrom:
            source=Path(options['architecture_annotation_cache'])/f'chr{chrom}.tsv.gz'
        elements=parse_joined_rmsk(source,chromosome=chrom)
    from allele_common import IntervalMask
    masks={name:IntervalMask(options[key],[chrom] if chrom else None) for name,key in
        [('lcr','architecture_lcr_bed'),('segdup','architecture_segdup_bed'),('low_mappability','architecture_low_mappability_bed')]
        if options.get(key)}
    pass
    tracks={name:read_label_bed(path,chromosome=chrom) for name,path in options.get('enrichment_tracks',{}).items()}
    return ArchitectureAccumulator(blocks,elements,thresholds=options['thresholds'],tracks=tracks,
        out_prefix=Path(shard_dir)/options.get('architecture_prefix','architecture'),masks=masks,
        flag_annotations={'lcr','segdup'}&set(header['info']))

def loco_expected(global_path,target_path,kind,min_training_sites=1):
    """Streaming targets; report total observations even when LOCO training is missing."""
    globals_=defaultdict(lambda:[0,0,0]); by_chr=defaultdict(lambda:[0,0,0])
    dims=['context','af_bin','an_bin','sensitivity','population_mode','threshold']
    def stratum(r): return tuple(r.get(x,{'sensitivity':'all','population_mode':'all','threshold':'1'}.get(x,'NA')) for x in dims)
    for r in _read_tsv(global_path):
        st=stratum(r); vals=[int(r[x]) for x in ('n_sites','n_multi','n_three')]
        for table,key in ((globals_,st),(by_chr,(r['chrom'],*st))):
            for i,v in enumerate(vals): table[key][i]+=v
    ids=(['ancestry','block_id','chrom','start','end'] if kind=='ld' else ['label','chrom'])+['sensitivity','population_mode','threshold']
    grouped=defaultdict(lambda:[0,0,0,0.,0.,0.,0])
    for r in _read_tsv(target_path):
        gid=tuple(r.get(k,{'sensitivity':'all','population_mode':'all','threshold':'1'}.get(k,'NA')) for k in ids)
        z=grouped[gid]; n,m,t=[int(r[k]) for k in ('n_sites','n_multi','n_three')]
        z[0]+=n; z[1]+=m; z[2]+=t
        st=stratum(r); train=[a-b for a,b in zip(globals_[st],by_chr[(r['chrom'],*st)])]
        if train[0]<min_training_sites: z[6]+=n; continue
        pm=(train[1]+.5)/(train[0]+1); pt=(train[2]+.5)/(train[0]+1)
        z[3]+=n*pm; z[4]+=n*pt; z[5]+=n*pm*(1-pm)
    rows=[]
    for gid,(n,m,t,em,et,var,missing) in grouped.items():
        rows.append(dict(zip(ids,gid))|dict(n_sites=n,observed_multi=m,observed_three=t,
            expected_multi=em if not missing else 'NA',expected_three=et if not missing else 'NA',
            expected_multi_supported=em,n_sites_without_loco=missing,
            oe=m/em if em and not missing else 'NA',oe_three=t/et if et and not missing else 'NA',
            z_descriptive=(m-em)/math.sqrt(var) if var and not missing else 'NA',
            status='complete' if not missing else 'insufficient_loco_training'))
    return sorted(rows,key=lambda r:(r['oe']=='NA',-r['oe'] if r['oe']!='NA' else 0,str(r)))

def write_summary(rows,path,fields=None):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    if fields is None: fields=list(rows[0]) if rows else ['status']
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)

def summarize_architecture(shards,out,options,min_training_sites=1):
    out=Path(out); prefix=options.get('architecture_prefix','architecture')
    for kind,header in HEADERS.items():
        with gzip.open(out/f'architecture.{kind}_strata.tsv.gz','wt',compresslevel=2,newline='') as f:
            w=csv.DictWriter(f,fieldnames=header,delimiter='\t'); w.writeheader()
            for shard in shards:
                source=Path(shard)/f'{prefix}.{kind}_strata.tsv.gz'
                if not source.exists(): raise ValueError(f'Architecture counter missing: {source}; old scans cannot be reused')
                w.writerows(_read_tsv(source))
    global_path=out/'architecture.global_strata.tsv.gz'
    for kind in ('ld','erv'):
        rows=loco_expected(global_path,out/f'architecture.{kind}_strata.tsv.gz',kind,min_training_sites)
        write_summary(rows,out/('ld_block_enrichment.tsv' if kind=='ld' else 'erv_enrichment.tsv'))
        if kind=='ld':
            ranked=[r for r in rows if r['block_id']!='outside_partition' and r['n_sites']>=options.get('architecture_min_sites',10000) and r['oe']!='NA']
            write_summary(ranked,out/'ld_block_ranked.tsv',list(rows[0]) if rows else None)
            consensus=hotspot_consensus(rows,options.get('architecture_min_sites',10000),options.get('architecture_top_fraction',.01))
            write_summary(consensus,out/'ld_hotspot_consensus.tsv',['chrom','start','end','overlap_bp','supporting_ancestries','n_ancestries','block_ids'])
    pass
    pass
    regions=[]
    ld_rows=list(_read_tsv(out/'ld_block_enrichment.tsv'))
    for threshold in (5,10):
        fragments=hotspot_consensus(ld_rows,options.get('architecture_min_sites',10000),options.get('architecture_top_fraction',.01),threshold=threshold)
        regions.extend(dict(r,threshold=threshold) for r in consensus_regions(fragments))
    write_summary(regions,out/'ld_hotspot_consensus_regions.tsv',['chrom','start','end','overlap_bp','supporting_ancestries','n_ancestries','block_ids','threshold'])
    summarize_exposure(global_path,out/'architecture.exposure_strata.tsv.gz',out,options.get('enrichment_bootstrap',2000),min_training_sites=min_training_sites)
    availability=[json.loads((Path(s)/f'{prefix}.availability.json').read_text()) for s in shards]
    dump_json(out/'architecture_manifest.json',dict(schema=SCHEMA,status='complete',shards=availability,
        inputs=annotation_identity(options),ld_primary='EUR',min_sites=options.get('architecture_min_sites',10000),
        missing_annotation='NA; unavailable annotation is never an observed zero',
        inference='Formal AC5/10 exposure inference: LOCO expected counts, EUR LD-block bootstrap, centered bootstrap P and BH FDR; legacy z is descriptive'))


# ========================================================================
# Ukb mfi
# ========================================================================

# Exact GRCh37 allele matching with explicit unavailable/absent/swapped states.

from collections import Counter, defaultdict
from pathlib import Path
import csv, re, math, json, sqlite3, tempfile, os
from allele_common import canonical_chrom, open_text as architecture_open_text, chrom_sort_key
from allele_common import dump_json, fingerprint, tsv_writer

VAR_RE=re.compile(r'^(?:chr)?(?P<chr>[0-9XYMT]+):(?P<pos>[0-9]+)[_:](?P<a1>[ACGT]+)[_:](?P<a2>[ACGT]+)$',re.I)
COORD_RE=re.compile(r'^(?:chr)?([0-9XYMT]+):([0-9]+)$',re.I)
ARCHITECTURE_COMPLEMENT=str.maketrans('ACGT','TGCA')

def _float(value,label,lo,hi,missing=False):
    if missing and value in {'.','NA','nan','-1'}: return None
    try: v=float(value)
    except ValueError as exc: raise ValueError(f'Invalid MFI {label}: {value}') from exc
    if not math.isfinite(v) or not lo<=v<=hi: raise ValueError(f'Invalid MFI {label}: {value}')
    return v

def parse_mfi_row(fields,chrom_hint=None):
    """7-col, 8-col (minor allele), or local 9-col split alternate ID.

    The local first real row verifies order; never silently skip malformed data.
    A 7-column MAF cannot identify ALT AF. Both minor/MAF orders are accepted
    only when the allele token is unambiguous.
    """
    n=len(fields)
    if n==9:
        coord,pair,rsid,pos,a1,a2,v1,v2,info=fields
        alt_id=coord+'_'+pair
    elif n==8: alt_id,rsid,pos,a1,a2,v1,v2,info=fields
    elif n==7:
        alt_id,rsid,pos,a1,a2,maf,info=fields; v1=maf; v2=None
    else: raise ValueError(f'Expected 7, 8 or 9 MFI columns, found {n}')
    a1=a1.upper(); a2=a2.upper()
    if not re.fullmatch('[ACGT]+',a1) or not re.fullmatch('[ACGT]+',a2) or a1==a2:
        raise ValueError('Invalid MFI alleles')
    try: pos=int(pos)
    except ValueError as exc: raise ValueError('Invalid MFI position') from exc
    if pos<1: raise ValueError('Invalid MFI position')
    minor=None
    if v2 is not None:
        if v1.upper() in {a1,a2}: minor=v1.upper(); maf=v2
        elif v2.upper() in {a1,a2}: minor=v2.upper(); maf=v1
        elif v2 in {'.','NA'}: maf=v1
        elif v1 in {'.','NA'}: maf=v2
        else: raise ValueError('MFI minor allele disagrees with REF/ALT')
    else: maf=v1
    maf=_float(maf,'MAF',0,.5,True); info=_float(info,'INFO',-1,1,True)
    match=VAR_RE.fullmatch(alt_id)
    coord=COORD_RE.fullmatch(alt_id)
    inferred=match.group('chr') if match else coord.group(1) if coord else None
    if chrom_hint is None and inferred is None: raise ValueError('MFI chromosome requires filename hint')
    chrom=canonical_chrom(chrom_hint or inferred)
    id_alleles=None
    if match:
        id_alleles=(match.group('a1').upper(),match.group('a2').upper())
        if canonical_chrom(match.group('chr'))!=chrom or int(match.group('pos'))!=pos or id_alleles!=(a1,a2):
            raise ValueError(f'MFI alternate_id disagrees with coordinates/alleles: {alt_id}')
    elif coord and (canonical_chrom(coord.group(1))!=chrom or int(coord.group(2))!=pos):
        raise ValueError('MFI alternate_id coordinate mismatch')
    alt_af=None if minor is None or maf is None else maf if minor==a2 else 1-maf
    return dict(chrom=chrom,pos=pos,ref=a1,alt=a2,maf=maf,info=info,minor_allele=minor,
        alt_af=alt_af,alternate_id=alt_id,rsid=rsid,id_alleles=id_alleles,columns=n)

def _records(path,chrom):
    with architecture_open_text(path) as f:
        for ln,line in enumerate(f,1):
            if not line.strip() or line.startswith('#'): continue
            fields=line.split()
            if fields[0].lower() in {'alternate_id','alternate_ids'}: continue
            try: yield parse_mfi_row(fields,chrom)
            except ValueError as exc: raise ValueError(f'{path}:{ln}: {exc}') from exc

def find_mfi(directory,c):
    names=[f'ukb_mfi_chr{c}_v3.txt']
    if c.isdigit(): names.append(f'ukb_mfi_chr{int(c):02d}_v3.txt')
    for name in names:
        for suffix in ('','.gz'):
            p=Path(directory)/(name+suffix)
            if p.is_file(): return p
    return None

def target_rows(path):
    with architecture_open_text(path) as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not {'chrom37','pos37','ref37'}.issubset(reader.fieldnames or []):
            raise ValueError('MFI targets require explicit GRCh37 chrom37,pos37,ref37 columns')
        for row in reader:
            c=canonical_chrom(row['chrom37']); pos=int(row['pos37']); ref=row['ref37'].upper()
            if row.get('alt37') not in (None,'','.','NA'):
                yield c,pos,ref,row['alt37'].upper(),'anchor',row
            if row.get('extra_alt') not in (None,'','.','NA'):
                strand=row.get('strand37to38')
                if strand not in {'+','-'}: raise ValueError('extra_alt needs strand37to38 for GRCh37 matching')
                alt=row['extra_alt'].upper()
                if strand=='-': alt=alt.translate(ARCHITECTURE_COMPLEMENT)
                yield c,pos,ref,alt,'extra',row

def scan_mfi_dir(mfi_dir,target_tsv,out_tsv,info_thresholds=(.3,.8)):
    """Disk-backed U0 targets; preserves every contrast/role and every absence.

    Matching is pos+REF+ALT. A swapped row is flagged but never exact-present.
    Missing chromosome files mean unavailable, not absent. Output includes both
    MAF and ALT AF where a minor-allele column identifies frequency orientation.
    """
    out=Path(out_tsv); out.parent.mkdir(parents=True,exist_ok=True)
    marker=Path(str(out)+'.manifest.json')
    if marker.exists():
        old=json.loads(marker.read_text())
        if old.get('schema')!=2 or old.get('inputs',{}).get('targets')!=fingerprint(target_tsv) or old.get('mfi_directory')!=str(Path(mfi_dir).resolve()) or old.get('info_thresholds')!=list(info_thresholds):
            raise ValueError('UKB MFI inputs/settings changed; choose a new output')
        for c,audit in old['chromosomes'].items():
            path=find_mfi(mfi_dir,c)
            if audit['status']=='unavailable_file':
                if path is not None: raise ValueError('Previously unavailable MFI chromosome now exists; choose a new output')
            elif path is None or fingerprint(path)!=audit['source']: raise ValueError('MFI source changed; choose a new output')
        if not out.exists() or fingerprint(out)!=old.get('output'): raise ValueError('Completed MFI output changed')
        return old
    with tempfile.TemporaryDirectory(prefix='mfi-',dir=out.parent) as tmp:
        db=sqlite3.connect(str(Path(tmp)/'targets.sqlite'))
        db.execute('PRAGMA temp_store=FILE'); db.execute('PRAGMA cache_size=-32000')
        db.execute('CREATE TABLE target (chrom TEXT,pos INT,ref TEXT,alt TEXT,role TEXT,payload TEXT)')
        batch=[]; fields=None; target_count=0
        for c,p,r,a,role,row in target_rows(target_tsv):
            if fields is None: fields=list(row)
            if len(r)!=1 or len(a)!=1 or r==a or r not in 'ACGT' or a not in 'ACGT': raise ValueError('U0 target is not an exact SNV contrast')
            batch.append((c,p,r,a,role,json.dumps(row))); target_count+=1
            if len(batch)>=10000: db.executemany('INSERT INTO target VALUES (?,?,?,?,?,?)',batch); batch=[]
        if batch: db.executemany('INSERT INTO target VALUES (?,?,?,?,?,?)',batch)
        db.execute('CREATE INDEX target_key ON target(chrom,pos,ref,alt)')
        db.execute('CREATE TABLE matches (chrom TEXT,pos INT,ref TEXT,alt TEXT,status TEXT,payload TEXT,PRIMARY KEY(chrom,pos,ref,alt))')
        db.commit()
        qc=Counter(); audits={}; missing=[]
        for (c,) in db.execute('SELECT DISTINCT chrom FROM target ORDER BY chrom'):
            path=find_mfi(mfi_dir,c)
            if path is None: missing.append(c); audits[c]=dict(status='unavailable_file'); continue
            want={(p,r,a) for p,r,a in db.execute('SELECT DISTINCT pos,ref,alt FROM target WHERE chrom=?',(c,))}
            positions={k[0] for k in want}; found={}; first=None; n=0
            for x in _records(path,c):
                if first is None: first=x
                n+=1
                if x['pos'] not in positions: continue
                for key,status in [((x['pos'],x['ref'],x['alt']),'present'),((x['pos'],x['alt'],x['ref']),'swapped_alleles')]:
                    if key not in want: continue
                    if key in found and found[key][0]=='present':
                        if status=='present': raise ValueError(f'Duplicate exact MFI key: {c}:{key}')
                        continue
                    found[key]=(status,x)
            if first is None: raise ValueError(f'Empty MFI file: {path}')
            audits[c]=dict(status='scanned',source=fingerprint(path),first_real_row=first,records=n)
            db.executemany('INSERT INTO matches VALUES (?,?,?,?,?,?)',[(c,*key,status,json.dumps(x)) for key,(status,x) in found.items()]); db.commit()
            print(f'[UKB MFI] chr{c}: {n:,} rows, {len(found):,}/{len(want):,} target keys located',flush=True)
        added=['target_role','ukb_chrom','ukb_pos','ukb_ref','ukb_alt','ukb_status','ukb_allele_swap','ukb_maf','ukb_info','ukb_minor_allele','ukb_alt_af','ukb_alternate_id','ukb_rsid','gnomad_nfe_af','af_difference']+[f'ukb_info_ge_{str(t).replace(".","p")}' for t in info_thresholds]
        fields=(fields or ['chrom37','pos37','ref37','alt37'])+added
        temp=out.with_name(out.name+'.partial')
        sums=[0,0.,0.,0.,0.,0.]
        with temp.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); writer.writeheader()
            query='SELECT t.chrom,t.pos,t.ref,t.alt,t.role,t.payload,m.status,m.payload FROM target t LEFT JOIN matches m USING(chrom,pos,ref,alt) ORDER BY t.rowid'
            for c,p,r,a,role,payload,status,match in db.execute(query):
                row=json.loads(payload); x=json.loads(match) if match else {}
                status=status or ('unavailable_file' if c in missing else 'absent_from_MFI')
                qc[status]+=1
                nfe=row.get('extra_nfe_af') if role=='extra' else None
                try: nfe=float(nfe)
                except (ValueError,TypeError): nfe=None
                if nfe is not None and (not math.isfinite(nfe) or not 0<=nfe<=1): nfe=None
                alt_af=x.get('alt_af') if status=='present' else None
                diff=alt_af-nfe if alt_af is not None and nfe is not None else None
                row.update(dict(target_role=role,ukb_chrom=c,ukb_pos=p,ukb_ref=x.get('ref'),ukb_alt=x.get('alt'),ukb_status=status,
                    ukb_allele_swap=int(status=='swapped_alleles') if match else None,
                    ukb_maf=x.get('maf'),ukb_info=x.get('info'),ukb_minor_allele=x.get('minor_allele'),ukb_alt_af=alt_af,
                    ukb_alternate_id=x.get('alternate_id'),ukb_rsid=x.get('rsid'),gnomad_nfe_af=nfe,af_difference=diff))
                for t in info_thresholds: row[f'ukb_info_ge_{str(t).replace(".","p")}']=int(x['info']>=t) if x.get('info') is not None else None
                if diff is not None and (x.get('info') or 0)>=.8:
                    sums[0]+=1; sums[1]+=nfe; sums[2]+=alt_af; sums[3]+=nfe*nfe; sums[4]+=alt_af*alt_af; sums[5]+=nfe*alt_af
                writer.writerow({k:'NA' if v is None else v for k,v in row.items()})
        os.replace(temp,out); db.close()
    n,sx,sy,sxx,syy,sxy=sums; denom=(n*sxx-sx*sx)*(n*syy-sy*sy)
    result=dict(schema=2,target_rows=target_count,matched_rows=qc['present'],counts=dict(qc),chromosomes=audits,
        inputs=dict(targets=fingerprint(target_tsv)),mfi_directory=str(Path(mfi_dir).resolve()),info_thresholds=list(info_thresholds),output=fingerprint(out),concordance=dict(n=n,info_floor=.8,
        pearson_r=(n*sxy-sx*sy)/math.sqrt(denom) if n>1 and denom>0 else None),
        caveats=['Absence from MFI is not extinction.','Swapped alleles are not exact matches.','Seven-column MAF has unknown ALT direction.','Concordance conditional on exact orientation and INFO>=0.8.'])
    dump_json(str(out)+'.manifest.json',result); return result


# ========================================================================
# Universes
# ========================================================================

# Independent U0 full contrasts, U1 probability sample, U2 all strong, U3 budget.

from collections import Counter, defaultdict
from pathlib import Path
import csv, gzip, heapq, json, math, random, shutil, sqlite3
from evolve import table_rows, number, _match_stratum, _stage_identity, completed_stage, finish_stage
from allele_common import stable_hash, tsv_writer, dump_json, fingerprint, open_text as architecture_open_text
from allele_common import ac_tier

EXTRA_FIELDS=['group','match_set','match_stratum','locus_block','selection_evidence','comparison_universe',
              'universe','sampling_stratum','sampling_probability','sampling_weight']

def eligible(row,min_extra_ac=5,exclude_mhc=True,allow_filter_warnings=False):
    return (number(row.get('extra_ac'),0)>=min_extra_ac and
        (row.get('selection_filter')=='PASS' or allow_filter_warnings) and
        row.get('has_non_snv')!='1' and row.get('region_flags','.') in ('.','','none') and
        not (exclude_mhc and row['chrom37']=='6' and 25000000<=int(row['pos37'])<=34000000) and
        number(row.get('s'),0)!=0)

def sampling_stratum(r):
    p=number(r['posterior'],0); s=abs(number(r['s'],0))
    pb='strong' if p>=.99 else 'moderate' if p>=.9 else 'intermediate' if p>.5 else 'weak'
    sb='zero' if not s else str(math.floor(math.log10(s)))
    return '|'.join([pb,sb,ac_tier(r['extra_ac']),r.get('canonical_context',r['context']),
                     r.get('population_class','NA'),r['n_modern_alt'],r['chrom']])

def sample_allocations(counts,budget,seed=20261022,with_expectations=False):
    """Fixed-size stratified sample with randomized fractional allocation.

    Expected stratum counts mu_h are proportional to sqrt(N_h) times declared
    oversampling factors, capped at N_h. Systematic randomized rounding gives
    E[K_h]=mu_h and sum K_h=budget, even when strata outnumber the budget.
    Independent SRS inside each stratum gives exact marginal pi=mu_h/N_h.
    """
    budget=min(budget,sum(counts.values()))
    if budget<0: raise ValueError('Sample budget must be nonnegative')
    weights={s:math.sqrt(n)*(4 if s.split('|')[0]=='strong' else 1)*(3 if s.split('|')[-2]=='3' else 1)*(2 if s.split('|')[2] in {'1','2-4'} else 1) for s,n in counts.items()}
    remaining=set(counts); expected={}; left=float(budget)
    while remaining:
        total=sum(weights[s] for s in sorted(remaining))
        full=[s for s in sorted(remaining) if left*weights[s]/total>=counts[s]]
        if not full:
            expected.update({s:left*weights[s]/total for s in remaining});break
        for st in full: expected[st]=float(counts[st]);left-=counts[st];remaining.remove(st)
    alloc={}; offset=random.Random('allocation:'+str(seed)).random(); cumulative=0.; previous=0
    ordered=sorted(counts)
    for i,st in enumerate(ordered):
        cumulative+=expected[st]
        rounded=budget if i==len(ordered)-1 else math.floor(cumulative+offset)
        alloc[st]=rounded-previous; previous=rounded
        if not 0<=alloc[st]<=counts[st]: raise RuntimeError('Invalid sampling allocation')
    return (alloc,expected) if with_expectations else alloc


def _decorate(r,universe,group,block,stratum='NA',probability=1,match='.'):
    return dict(r,group=group,match_set=match,match_stratum=_match_stratum(r),locus_block=block,
        selection_evidence='published_strong' if group=='strong' else 'weak_signal_not_neutral' if group=='weak_signal' else 'not_selected_by_significance',
        comparison_universe='modern_multi_alt_exact_intersection',universe=universe,
        sampling_stratum=stratum,sampling_probability=probability,sampling_weight=1/probability)

def clustered_bootstrap(rows,out,seed=20261022,reps=1000):
    """Descriptive all-strong means, resampling whole GRCh37 blocks with replacement."""
    groups=defaultdict(lambda:[0,0.,0.])
    for r in rows:
        z=groups[r['locus_block']]; z[0]+=1; z[1]+=float(r['extra_af']); z[2]+=float(r['extra_ac'])
    values=list(groups.values()); rng=random.Random(seed); result=[]
    for index,name in ((1,'extra_af'),(2,'extra_ac')):
        means=[]
        if len(values)>=2:
            for _ in range(reps):
                chosen=[values[rng.randrange(len(values))] for _ in values]
                means.append(sum(z[index] for z in chosen)/sum(z[0] for z in chosen))
            means.sort()
        result.append(dict(metric=name,n_rows=len(rows),n_clusters=len(values),
            mean=sum(z[index] for z in values)/len(rows) if rows else 'NA',
            ci_low=means[int(.025*len(means))] if means else 'NA',
            ci_high=means[min(len(means)-1,int(.975*len(means)))] if means else 'NA',
            method='whole_LD_block_bootstrap; outside_partition uses labelled 1Mb spatial clusters'))
    h,w=tsv_writer(out,list(result[0]))
    with h: w.writerows(result)

def build_universes(intersection,out_dir,*,ld_block_dir,u1_sites=20000,strong_posterior=.99,
                    weak_posterior=.5,min_extra_ac=5,controls_per_case=4,seed=20261022,
                    exclude_mhc=True,allow_filter_warnings=False):
    if u1_sites<1 or not 3<=controls_per_case<=5 or not 0<=weak_posterior<strong_posterior<=1:
        raise ValueError('Invalid universe budget, posterior thresholds or controls (3-5 required)')
    block_path=Path(ld_block_dir)/'EUR.37.bed'; blocks,_=read_blocks(block_path,'EUR','37'); index=BlockIndex(blocks)
    out=Path(out_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    params=dict(schema=2,u1_sites=u1_sites,strong_posterior=strong_posterior,weak_posterior=weak_posterior,
        min_extra_ac=min_extra_ac,controls_per_case=controls_per_case,seed=seed,
        exclude_mhc=exclude_mhc,allow_filter_warnings=allow_filter_warnings)
    identity=_stage_identity('independent_universes',dict(intersection=intersection,ld=block_path),params)
    if old:=completed_stage(out,identity): return old
    dbpath=out/'sampling.sqlite'
    if dbpath.exists(): dbpath.unlink()
    db=sqlite3.connect(dbpath); db.execute('PRAGMA cache_size=-48000')
    db.execute('CREATE TABLE sites (site TEXT PRIMARY KEY,stratum TEXT,rand TEXT,rank REAL,payload TEXT,block TEXT)')
    db.execute('CREATE TABLE weak (site TEXT,stratum TEXT,rand TEXT,payload TEXT,block TEXT)')
    fields=None; strong=[]; qc=Counter()
    def block(r):
        b=index.locate(r['chrom37'],r['pos37'])
        return b.block_id if b else f'outside_partition:{r["chrom37"]}:{(int(r["pos37"])-1)//1000000}'
    # U0 is byte-preserved if compressed; otherwise losslessly recompressed.
    if str(intersection).endswith('.gz'): shutil.copyfile(intersection,out/'U0.tsv.gz')
    else:
        with open(intersection,'rb') as src,gzip.open(out/'U0.tsv.gz','wb') as dst: shutil.copyfileobj(src,dst)
    for r in table_rows(intersection):
        if fields is None: fields=list(r)
        qc['U0_rows']+=1; st=sampling_stratum(r); b=block(r)
        rank=float(r['posterior']); payload=json.dumps(r); h=stable_hash([seed,r['site_id']])
        # One deterministic representative per unique site for probability sampling.
        db.execute('INSERT INTO sites VALUES (?,?,?,?,?,?) ON CONFLICT(site) DO UPDATE SET stratum=excluded.stratum,rank=excluded.rank,payload=excluded.payload,block=excluded.block WHERE excluded.rank>sites.rank',
            (r['site_id'],st,h,rank,payload,b))
        if eligible(r,min_extra_ac,exclude_mhc,allow_filter_warnings):
            p=float(r['posterior'])
            if p>=strong_posterior: strong.append(_decorate(r,'U2','strong',b)); qc['U2_all_strong_rows']+=1
            elif p<=weak_posterior:
                db.execute('INSERT INTO weak VALUES (?,?,?,?,?)',(r['site_id'],_match_stratum(r),h,payload,b))
        if qc['U0_rows']%50000==0: db.commit()
    db.commit(); db.execute('CREATE INDEX site_sample ON sites(stratum,rand)'); db.execute('CREATE INDEX weak_sample ON weak(stratum,rand)'); db.commit()
    fields=(fields or [])+EXTRA_FIELDS
    counts=dict(db.execute('SELECT stratum,COUNT(*) FROM sites GROUP BY stratum')); alloc,expected=sample_allocations(counts,u1_sites,seed,True)
    u1=[]
    for st,k in sorted(alloc.items()):
        for payload,b in db.execute('SELECT payload,block FROM sites WHERE stratum=? ORDER BY rand,site LIMIT ?',(st,k)):
            r=json.loads(payload); u1.append(_decorate(r,'U1','architecture_sample',b,st,expected[st]/counts[st]))
    qc['U0_unique_sites']=sum(counts.values()); qc['U1_unique_sites']=len(u1)
    strong.sort(key=lambda r:(-float(r['posterior']),-abs(float(r['x'])),r['site_id'],r['selection_id'],r['extra_alt']))
    strong_sites={r['site_id'] for r in strong}; used=set(strong_sites); controls=[]; match_report=[]
    for i,case in enumerate(strong,1):
        match=f'case{i:06d}'; case['match_set']=match; chosen=[]
        for site,_,_,payload,b in db.execute('SELECT * FROM weak WHERE stratum=? ORDER BY rand,site',(case['match_stratum'],)):
            if site in used or b==case['locus_block']: continue
            chosen.append(_decorate(json.loads(payload),'U2','weak_signal',b,match=match)); used.add(site)
            if len(chosen)==controls_per_case: break
        controls.extend(chosen)
        match_report.append(dict(match_set=match,site_id=case['site_id'],selection_id=case['selection_id'],extra_alt=case['extra_alt'],
            requested=controls_per_case,matched=len(chosen),status='complete' if len(chosen)==controls_per_case else 'insufficient_exact_matched_weak_controls'))
    qc['U2_controls']=len(controls); qc['U2_cases_under_3_controls']=sum(r['matched']<3 for r in match_report)
    reps={}
    for r in strong:
        if not r['locus_block'].startswith('outside_partition'): reps.setdefault(r['locus_block'],r)
    qc['U2_representatives']=len(reps); qc['U2_outside_partition_rows']=sum(r['locus_block'].startswith('outside_partition') for r in strong)
    outputs=[]
    def write(name,rows,cols=fields):
        h,w=tsv_writer(out/name,cols)
        with h: w.writerows(rows)
        outputs.append(name)
    write('U1.tsv',u1); write('U2_all_strong.tsv',strong); write('U2_controls.tsv',controls)
    write('U2_representatives.tsv',list(reps.values()))
    write('U2_matching.tsv',match_report,['match_set','site_id','selection_id','extra_alt','requested','matched','status'])
    write('sampling_design.tsv',[dict(stratum=s,population=counts[s],sample=alloc[s],expected_sample=expected[s],probability=expected[s]/counts[s],weight=counts[s]/expected[s]) for s in sorted(counts)],['stratum','population','sample','expected_sample','probability','weight'])
    combined={}
    for r in [*strong,*controls,*u1]: combined.setdefault((r['site_id'],r['selection_id'],r['extra_alt']),r)
    write('candidates.tsv',list(combined.values()))
    clustered_bootstrap(strong,out/'U2_clustered_bootstrap.tsv',seed)
    clustered_bootstrap(list(reps.values()),out/'U2_representative_bootstrap.tsv',seed)
    db.close()
    outputs+=['U0.tsv.gz','U2_clustered_bootstrap.tsv','U2_representative_bootstrap.tsv','sampling.sqlite']
    return finish_stage(out,identity,dict(qc=dict(qc),candidate_path=str(out/'candidates.tsv'),
        U0_path=str(out/'U0.tsv.gz'),U3_status='not_selected; run select-u3 after local LD/ERV/UKB review',
        caveats=['U0 is unsampled.','U1 marginal probability is mu_h/N_h; randomized systematic stratum allocation, independent within-stratum SRS; all units have positive probability.',
        'U2 includes every eligible strong contrast; no max-sites truncation.','Matching never relaxes scientific thresholds; deficits are explicit.',
        'Outside LD partitions cannot be labelled LD representatives.','Matching extra AF precludes paired AF-outcome tests.']),outputs)

def select_u3(universe_dir,out_path,sites=3000,seed=20261022):
    if sites<1: raise ValueError('U3 size must be positive')
    directory=Path(universe_dir);src=directory/'sampling.sqlite'
    # A separate functional budget: keep U2 cases and their controls first,
    # then fill from U0. U1 sampling is never consulted.
    selected={};matches=set()
    for r in table_rows(directory/'U2_all_strong.tsv'):
        if r['site_id'] in selected:continue
        if len(selected)>=sites:break
        selected[r['site_id']]=dict(r,universe='U3',sampling_stratum='U3_priority_U2')
        matches.add(r['match_set'])
    pools=defaultdict(list)
    for r in table_rows(directory/'U2_controls.tsv'):
        if r['match_set'] in matches:pools[r['match_set']].append(r)
    for round_ in range(5):
        for match in sorted(pools):
            if len(selected)>=sites:break
            if round_<len(pools[match]):
                r=pools[match][round_]
                selected.setdefault(r['site_id'],dict(r,universe='U3',sampling_stratum='U3_priority_U2'))
    priority=len(selected);budget=sites-priority
    db=sqlite3.connect(f'file:{src}?mode=ro',uri=True);heap=[];total=0;pool_n=0
    for site,payload,block in db.execute('SELECT site,payload,block FROM sites'):
        total+=1
        if site in selected:continue
        pool_n+=1
        if not budget:continue
        h=int(stable_hash(['U3',seed,site]),16);item=(-h,site,payload,block)
        if len(heap)<budget:heapq.heappush(heap,item)
        elif item[:2]>heap[0][:2]:heapq.heapreplace(heap,item)
    db.close();result_rows=list(selected.values())
    for _,_,payload,block in sorted(heap,reverse=True):
        r=json.loads(payload);result_rows.append(_decorate(r,'U3','functional_sample',block,'U0_independent_fill',min(budget,pool_n)/pool_n))
    if result_rows:
        h,w=tsv_writer(out_path,list(result_rows[0]))
        with h:w.writerows(result_rows)
    else:
        with architecture_open_text(directory/'U0.tsv.gz') as f:fields=next(csv.reader(f,delimiter='\t'))+EXTRA_FIELDS
        h,w=tsv_writer(out_path,fields);h.close()
    result=dict(n_sites=len(result_rows),U0_unique_sites=total,requested_sites=sites,U2_priority_sites=priority,
        independent_of_U1=True,source=fingerprint(src),design='One row per site; U2 strong then round-robin matched controls; independent U0 seeded fill; targeted functional panel')
    dump_json(str(out_path)+'.json',result);return result


# ========================================================================
# Final cli
# ========================================================================

# FINAL utilities use the same modules as the integrated scanner/pipeline.

import argparse, json, sys
from pathlib import Path
from allele_common import dump_json

def main_architecture(argv=None):
    p=argparse.ArgumentParser(prog='allele3-final'); sp=p.add_subparsers(dest='cmd',required=True)
    q=sp.add_parser('validate-ld'); q.add_argument('--dir',required=True); q.add_argument('--build',choices=['37','38'],default='38'); q.add_argument('--out',required=True)
    q=sp.add_parser('prepare-erv'); q.add_argument('--joined-rmsk',required=True); q.add_argument('--out',required=True); q.add_argument('--cache-dir'); q.add_argument('--chrom')
    q=sp.add_parser('summarize-architecture'); q.add_argument('--global-strata',required=True); q.add_argument('--ld-strata'); q.add_argument('--erv-strata'); q.add_argument('--out-dir',required=True); q.add_argument('--min-sites',type=int,default=10000); q.add_argument('--min-training-sites',type=int,default=1)
    q=sp.add_parser('ukb-mfi'); q.add_argument('--mfi-dir',required=True); q.add_argument('--targets',required=True); q.add_argument('--out',required=True)
    q=sp.add_parser('build-universes'); q.add_argument('--intersection',required=True); q.add_argument('--out',required=True); q.add_argument('--ld-block-dir',required=True); q.add_argument('--u1-sites',type=int,default=20000); q.add_argument('--controls-per-case',type=int,default=4); q.add_argument('--seed',type=int,default=20261022)
    q=sp.add_parser('select-u3'); q.add_argument('--universes',required=True); q.add_argument('--out',required=True); q.add_argument('--u3-sites',type=int,default=3000); q.add_argument('--seed',type=int,default=20261022)
    q=sp.add_parser('final-report'); q.add_argument('--modern-dir',required=True); q.add_argument('--out',required=True); q.add_argument('--metag-census'); q.add_argument('--no-plots',action='store_true')
    a=p.parse_args(argv)
    if a.cmd=='validate-ld':
        load_block_sets(a.dir,a.build); result=validate_directory(a.dir,a.build); dump_json(a.out,result)
    elif a.cmd=='prepare-erv':
        if a.cache_dir: result=prepare_joined_cache(a.joined_rmsk,a.cache_dir)
        else:
            e=parse_joined_rmsk(a.joined_rmsk,chromosome=a.chrom); write_erv_bed(e,a.out)
            from collections import Counter
            result=dict(n_elements=len(e),classes=dict(Counter(x.rep_class for x in e)),tiers=dict(Counter(x.tier for x in e)),source=a.joined_rmsk)
        dump_json(str(a.out)+'.manifest.json',result)
    elif a.cmd=='summarize-architecture':
        od=Path(a.out_dir);od.mkdir(parents=True,exist_ok=True);result={}
        for kind,path in [('ld',a.ld_strata),('erv',a.erv_strata)]:
            if path:
                r=loco_expected(a.global_strata,path,kind,a.min_training_sites)
                write_summary(r,od/('ld_block_enrichment.tsv' if kind=='ld' else 'erv_enrichment.tsv'));result[kind]=len(r)
                if kind=='ld': write_summary(hotspot_consensus(r,a.min_sites),od/'ld_hotspot_consensus.tsv')
    elif a.cmd=='ukb-mfi': result=scan_mfi_dir(a.mfi_dir,a.targets,a.out)
    elif a.cmd=='build-universes': result=build_universes(a.intersection,a.out,ld_block_dir=a.ld_block_dir,u1_sites=a.u1_sites,controls_per_case=a.controls_per_case,seed=a.seed)
    elif a.cmd=='select-u3': result=select_u3(a.universes,a.out,a.u3_sites,a.seed)
    elif a.cmd=='final-report':
        from report import write_local_review
        result=write_local_review(a.modern_dir,a.out,a.metag_census,not a.no_plots)
    print(json.dumps(result,indent=2,ensure_ascii=False));return 0


# Mappability annotation preparation
# Convert published GRCh38 Umap multi-read scores into a declared exclusion BED.
# 
# Source: https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bedgraph.gz
# Documentation: https://bismap.hoffmanlab.org/ (omitted intervals have score zero).
# Only primary autosomes are used; the reference FASTA index defines their lengths.
import argparse, gzip, json, hashlib, sys
from pathlib import Path
from allele_common import fingerprint, dump_json


def prepare_mappability(source,fai,out,threshold=1.0):
    if not 0<threshold<=1:raise ValueError('Threshold must be in (0,1]')
    lengths={}
    for line in Path(fai).read_text().splitlines():
        c,n,*_=line.split();c=c.removeprefix('chr')
        if c.isdigit() and 1<=int(c)<=22:lengths[c]=int(n)
    if not lengths:raise ValueError('No primary autosomes in FASTA index')
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);tmp=out.with_suffix(out.suffix+'.partial')
    seen=set();current=None;end=0;pending=None;intervals=bp=0
    def records():
        if str(source).endswith('.bw'):
            import pyBigWig
            import numpy as np
            with pyBigWig.open(str(source)) as track:
                for c in sorted(lengths,key=int):
                    name='chr'+c
                    if track.chroms().get(name)!=lengths[c]:raise ValueError('Umap/FASTA chromosome length mismatch: '+name)
                    for start in range(0,lengths[c],1000000):
                        end=min(lengths[c],start+1000000)
                        values=track.values(name,start,end,numpy=True)
                        if np.any((values<0)|(values>1)):raise ValueError('Invalid Umap score')
                        # Expand in bounded chunks and reduce to passing runs in NumPy.
                        # This avoids Python iteration over hundreds of millions of
                        # distinct fractional scores. Unwritten bases are NaN/zero.
                        passing=np.isfinite(values)&(values>=threshold)
                        edges=np.flatnonzero(np.diff(np.concatenate(([False],passing,[False]))))
                        for s,e in zip(edges[::2],edges[1::2]):
                            yield f'{name} {start+s} {start+e} 1'
        else:
            with gzip.open(source,'rt') as f:yield from f
    with gzip.open(tmp,'wt') as w:
        f=records()
        def emit(c,s,e):
            nonlocal pending,intervals,bp
            if e<=s:return
            if pending and pending[0]==c and s<=pending[2]:pending=(c,pending[1],max(e,pending[2]));return
            if pending:
                w.write(f'chr{pending[0]}\t{pending[1]}\t{pending[2]}\n');intervals+=1;bp+=pending[2]-pending[1]
            pending=(c,s,e)
        for line in f:
            if not line.strip() or line.startswith(('#','track','browser')):continue
            c,s,e,v=line.split()[:4];c=c.removeprefix('chr')
            if c not in lengths:continue
            s,e,v=int(s),int(e),float(v)
            if not 0<=s<e<=lengths[c] or not 0<=v<=1:raise ValueError('Invalid Umap coordinate or score')
            if c!=current:
                if current is not None:emit(current,end,lengths[current])
                if c in seen:raise ValueError('Umap chromosomes must be contiguous')
                seen.add(c);current=c;end=0
            if s<end:raise ValueError('Overlapping or unsorted Umap bedGraph')
            emit(c,end,s)
            if v<threshold:emit(c,s,e)
            end=e
        if current is not None:emit(current,end,lengths[current])
        for c in sorted(set(lengths)-seen,key=int):emit(c,0,lengths[c])
        if pending:w.write(f'chr{pending[0]}\t{pending[1]}\t{pending[2]}\n');intervals+=1;bp+=pending[2]-pending[1]
    tmp.replace(out)
    digest=hashlib.sha256()
    with open(source,'rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):digest.update(part)
    result=dict(status='complete',source=fingerprint(source),source_sha256=digest.hexdigest(),documentation='https://bismap.hoffmanlab.org/',download_url=('https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw' if str(source).endswith('.bw') else 'https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bedgraph.gz'),build='GRCh38',k=100,measure='multi-read mappability',exclusion=f'score < {threshold:g}; omitted intervals treated as zero per publisher documentation',threshold=threshold,intervals=intervals,excluded_bp=bp,autosomes=sorted(lengths,key=int),output=str(out))
    dump_json(str(out)+'.json',result);return result


def main_mappability():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True);p.add_argument('--fai',required=True);p.add_argument('--out',required=True);p.add_argument('--threshold',type=float,default=1.)
    a=p.parse_args();print(json.dumps(prepare_mappability(a.source,a.fai,a.out,a.threshold),indent=2))


# ========================================================================
# scan_scores
# ========================================================================

# Candidate cache and local precomputed score annotation.





# ========================================================================
# Candidate cache
# ========================================================================

# Derive theoretical ALT candidates only when an annotation command needs them.
# 
# The input run must already have passed CLI load_run/verify_shard. Original
# shards are read-only. Committed eager tables remain usable without rebuilding;
# deferred tables are reconstructed with bounded memory from matching site and
# allele groups. Raw formatted AC/AN/AF strings are copied, never recomputed.
# 
# Derived outputs use the workflow's explicitly labelled size/mtime/edge hash
# fingerprints, not a claimed full-file checksum. Each shard is committed
# atomically and verified on reuse; interrupted shards can be regenerated while
# completed ones are reused. A completed but modified cache is never overwritten.


from contextlib import ExitStack
import csv
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import tempfile

from allele_common import dump_json, fingerprint, open_text as scores_open_text, stable_hash, tsv_writer

SCORES_SCHEMA = 1
CANDIDATES = "candidates.tsv.gz"
SITES = "multiallelic_sites.tsv.gz"
ALLELES = "multiallelic_alleles.tsv.gz"
SITE_REQUIRED = {"chrom", "pos", "ref", "context", "cpg_site", "n_alt", "primary_alt", "primary_tie"}
ALLELE_REQUIRED = {"chrom", "pos", "ref", "alt", "status", "is_primary", "ac", "an", "af"}
PRESENT_STATUSES = {"observed_pass", "present_failed_filter", "present_missing_counts", "present_below_threshold"}


def _rows(path, required):
    with scores_open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not required.issubset(fields):
            raise ValueError(f"Invalid candidate source schema: {path}; required {sorted(required)}")
        for row in reader:
            if None in row or any(row.get(key) is None for key in required):
                raise ValueError(f"Malformed candidate source row at {path}:{reader.line_num}")
            yield row


def _coordinate(row):
    try:
        pos = int(row["pos"])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid candidate source position: {row.get('pos')}") from exc
    if pos < 1 or not row["chrom"] or len(row["ref"]) != 1 or row["ref"] not in BASES:
        raise ValueError("Candidate sources require positive positions and A/C/G/T REF")
    return row["chrom"], pos, row["ref"]


def _write_candidates(site_path, allele_path, output_path, expected_rows):
    """Join the two already synchronized shard tables, never a genome-size dict."""
    with ExitStack() as stack:
        sites = _rows(site_path, SITE_REQUIRED)
        alleles = _rows(allele_path, ALLELE_REQUIRED)
        stack.callback(sites.close)
        stack.callback(alleles.close)
        groups = itertools.groupby(alleles, key=_coordinate)
        handle, writer = tsv_writer(output_path, CANDIDATE_FIELDS)
        stack.enter_context(handle)
        rows_written = sites_written = 0
        previous = None
        completed_chroms = set()
        for site in sites:
            key = _coordinate(site)
            chrom, pos, ref = key
            if previous:
                if chrom == previous[0] and pos <= previous[1]:
                    raise ValueError(f"Duplicate or unsorted source site: {chrom}:{pos}")
                if chrom != previous[0]:
                    completed_chroms.add(previous[0])
                    if chrom in completed_chroms:
                        raise ValueError(f"Non-contiguous source chromosome: {chrom}")
            previous = key
            try:
                allele_key, group = next(groups)
            except StopIteration as exc:
                raise ValueError(f"Missing allele group for candidate site {key}") from exc
            if key != allele_key:
                raise ValueError(f"Site/allele candidate sources disagree: {key} versus {allele_key}")
            entries = {}
            for allele in group:
                alt = allele["alt"]
                if len(alt) != 1 or alt not in BASES or alt == ref or alt in entries:
                    raise ValueError(f"Invalid or duplicate source ALT at {key}: {alt}")
                if allele["status"] not in PRESENT_STATUSES:
                    raise ValueError(f"Unknown source ALT status at {key}: {allele['status']}")
                if allele["is_primary"] != str(int(alt == site["primary_alt"])):
                    raise ValueError(f"Source primary ALT labels disagree at {key}")
                entries[alt] = allele
            observed = {alt for alt, entry in entries.items() if entry["status"] == "observed_pass"}
            if site["n_alt"] not in ("2", "3") or len(observed) != int(site["n_alt"]):
                raise ValueError(f"Source ALT occupancy disagrees at {key}")
            if site["primary_alt"] not in observed or site["primary_tie"] not in ("0", "1"):
                raise ValueError(f"Invalid primary ALT/tie at {key}")
            if len(site["context"]) != 3 or site["context"][1] != ref or site["cpg_site"] not in ("0", "1"):
                raise ValueError(f"Invalid source context/CpG flag at {key}")
            for alt in BASES:
                if alt == ref:
                    continue
                entry = entries.get(alt)
                transition = int(relation(ref, alt) == "transition")
                writer.writerow(dict(chrom=site["chrom"], pos=site["pos"], ref=ref, alt=alt,
                    status=entry["status"] if entry else "not_listed_in_release",
                    is_primary=int(alt == site["primary_alt"]), primary_tie=site["primary_tie"],
                    n_alt=site["n_alt"], ac=entry["ac"] if entry else ".",
                    an=entry["an"] if entry else ".", af=entry["af"] if entry else ".",
                    context=site["context"], primary_alt=site["primary_alt"],
                    transition=transition, cpg_transition=int(site["cpg_site"] == "1" and transition)))
                rows_written += 1
            sites_written += 1
        if next(groups, None) is not None:
            raise ValueError("Unmatched trailing allele groups in candidate sources")
        if rows_written != expected_rows:
            raise ValueError(f"Candidate row count differs from scan: expected {expected_rows}, wrote {rows_written}")
    return dict(rows=rows_written, sites=sites_written)


def _verified_source(shard, marker, name):
    expected = marker.get("outputs", {}).get(name)
    path = shard / name
    if not isinstance(expected, dict) or not path.is_file():
        raise ValueError(f"Missing committed candidate source file: {path}")
    actual = fingerprint(path)
    if any(actual[key] != expected.get(key) for key in ("size", "edge_sha256")):
        raise ValueError(f"Committed candidate source changed: {path}")
    return actual


def _expected_rows(scan):
    expected = scan.get("expected_candidate_rows")
    n_multi = scan.get("qc", {}).get("multi_alt_sites")
    if expected is None and type(n_multi) is int:
        expected = 3*n_multi  # Old eager scans recorded this QC count.
    if expected is not None and (type(expected) is not int or expected < 0 or expected % 3):
        raise ValueError("Invalid expected_candidate_rows in scan metadata")
    if type(n_multi) is int and expected is not None and expected != 3*n_multi:
        raise ValueError("Candidate row expectation disagrees with multi-ALT scan count")
    return expected


def _verify_derived(target, signature, identity, expected_rows):
    marker = target / "COMPLETE.json"
    if not marker.is_file():
        raise ValueError(f"Unrecognized/incomplete derived candidate cache preserved at {target}; use a new cache directory")
    saved = json.loads(marker.read_text())
    if (saved.get("complete") is not True or saved.get("signature") != signature
            or saved.get("identity") != identity or saved.get("rows") != expected_rows):
        raise ValueError(f"Derived candidate cache identity/count changed: {target}")
    path = target / CANDIDATES
    if not path.is_file():
        raise ValueError(f"Derived candidate cache output is missing: {path}")
    actual, expected = fingerprint(path), saved.get("output", {})
    if any(actual[key] != expected.get(key) for key in ("size", "mtime_ns", "edge_sha256")):
        raise ValueError(f"Derived candidate cache output changed: {path}; preserve it and use a new cache directory")
    return path


def candidate_files_for_run(run, out_dir):
    """Return reusable candidate paths for a previously verified local scan run.

    ``out_dir`` is a separate derived-cache directory, e.g. OUT/candidate_cache.
    It is not created when all shards already contain committed eager tables.
    No original shard, completed cache, or unknown pre-existing output is removed.
    """
    if run.get("status") not in ("scanned", "complete") or not run.get("signature"):
        raise ValueError("candidate_files_for_run requires a verified scanned/complete run")
    shards = [Path(p).expanduser().resolve() for p in run["shards"]]
    if len(set(shards)) != len(shards):
        raise ValueError("Duplicate shards in candidate source run")
    cache_root = Path(out_dir).expanduser().resolve()
    if any(cache_root == shard or cache_root.is_relative_to(shard) for shard in shards):
        raise ValueError("Derived candidate cache must be outside original scan shards")
    generator_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(files=[], source_shards=[], expected_candidate_rows=0,
                  reused_eager_shards=0, generated_shards=0, reused_derived_shards=0,
                  cache_dir=str(cache_root), schema=SCORES_SCHEMA,
                  validation="metadata/edge fingerprints; not full-file cryptographic checksums")
    all_expected_known = True
    for shard in shards:
        marker_path = shard / "COMPLETE.json"
        if not marker_path.is_file():
            raise ValueError(f"Missing original shard completion marker: {marker_path}")
        marker = json.loads(marker_path.read_text())
        if marker.get("signature") != run["signature"]:
            raise ValueError(f"Original shard/run identity mismatch at {shard}")
        scan_fingerprint = _verified_source(shard, marker, "scan.json")
        scan = json.loads((shard / "scan.json").read_text())
        expected_rows = _expected_rows(scan)
        if expected_rows is None:
            all_expected_known = False
        else:
            result["expected_candidate_rows"] += expected_rows
        if CANDIDATES in marker.get("outputs", {}):
            if scan.get("candidate_storage", "eager") != "eager":
                raise ValueError(f"Candidate storage metadata contradicts committed eager table at {shard}")
            _verified_source(shard, marker, CANDIDATES)
            candidate_path = shard / CANDIDATES
            mode = "reused_eager"
            result["reused_eager_shards"] += 1
            cache_manifest = None
        else:
            if scan.get("candidate_storage") != "deferred":
                raise ValueError(f"No committed eager candidate table and no explicit deferred mode at {shard}")
            if expected_rows is None:
                raise ValueError(f"Deferred scan lacks expected candidate count at {shard}")
            sources = {name: _verified_source(shard, marker, name) for name in (SITES, ALLELES)}
            sources["scan.json"] = scan_fingerprint
            identity = dict(schema=SCORES_SCHEMA, columns=CANDIDATE_FIELDS, generator_sha256=generator_sha256,
                            source_run_signature=run["signature"], source_shard=str(shard),
                            sources=sources, expected_rows=expected_rows)
            signature = stable_hash(identity)
            cache_root.mkdir(parents=True, exist_ok=True)
            target = cache_root / signature
            with (cache_root / (signature + ".lock")).open("a") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                if target.exists():
                    candidate_path = _verify_derived(target, signature, identity, expected_rows)
                    mode = "reused_derived"
                    result["reused_derived_shards"] += 1
                else:
                    # A new uniquely owned staging directory means hard-kill
                    # leftovers are preserved, not adopted or silently deleted.
                    with tempfile.TemporaryDirectory(prefix=signature + ".partial-", dir=cache_root) as temporary:
                        stage = Path(temporary)
                        stats = _write_candidates(shard / SITES, shard / ALLELES, stage / CANDIDATES, expected_rows)
                        if any(fingerprint(shard / name) != before for name, before in sources.items()):
                            raise ValueError(f"Original candidate sources changed during reconstruction: {shard}")
                        output = fingerprint(stage / CANDIDATES)
                        output["path"] = str(target / CANDIDATES)
                        dump_json(stage / "COMPLETE.json", dict(signature=signature, identity=identity,
                                  output=output, rows=stats["rows"], sites=stats["sites"], complete=True))
                        os.replace(stage, target)
                    candidate_path = target / CANDIDATES
                    mode = "generated"
                    result["generated_shards"] += 1
            cache_manifest = str(target / "COMPLETE.json")
        result["files"].append(str(candidate_path))
        result["source_shards"].append(dict(shard=str(shard), mode=mode,
                    candidate_file=str(candidate_path), expected_rows=expected_rows, cache_manifest=cache_manifest))
    if not all_expected_known:
        result["expected_candidate_rows"] = None
    return result


# ========================================================================
# Local scores
# ========================================================================

# Offline, optional annotation of same-site SNV candidates with local scores.
# 
# No prediction, download, API request, P value or selection inference is performed.
# The streaming backend requires candidate shards and scores in natural
# chromosome/position order; alleles at the same position may occur in any order.
# The tabix backend also accepts other chromosome orders when each chromosome
# forms one contiguous block with sorted positions. It streams one non-overlapping
# coordinate block at a time, so it does not materialize a million-base block.
# 
# Candidate coordinates are 1-based. ``score_build`` explicitly asserts the
# shared assembly of candidate and score inputs; the CLI verifies this against
# the scan manifest. Incompatible input declarations are rejected. TSV score
# input requires a header. Explicit ``score_format=
# 'avi'`` additionally accepts the official six-column, headerless AVI format:
# chrom, pos, ref, alt, raw_score, phred. Official AVI input requires GRCh38.


import csv
import gzip
import importlib
import json
import math
import os
import re
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Iterator


REQUIRED_CANDIDATE_COLUMNS = {
    "chrom", "pos", "ref", "alt", "status", "is_primary", "n_alt",
    "ac", "an", "af", "context",
}
SCORE_COLUMNS = ["local_score", "score_status", "score_reason"]
PAIR_COLUMNS = [
    "chrom", "pos", "ref", "context", "n_alt", "primary_alt",
    "observed_alt", "unobserved_alt", "primary_score", "observed_score",
    "unobserved_score", "delta_observed_minus_unobserved",
    "observed_ac", "observed_an", "observed_af",
]
ALLOWED_STATUSES = {
    "observed_pass", "not_listed_in_release", "present_failed_filter",
    "present_missing_counts", "present_below_threshold",
}
AVI_CAVEAT = (
    "AVI was trained using gnomAD v4.1 allele-frequency proxy labels. Its "
    "association with gnomAD frequency or occupancy is not independent evidence "
    "for selection. AVI, scalar score similarity and AVI feature attributions "
    "are not counts of biological functions or measures of pleiotropy. This "
    "module provides exploratory annotations and arithmetic score differences "
    "only; it performs no hypothesis test or causal inference."
)
_MISSING = {"", ".", "na", "n/a", "null", "none"}
_ALIASES = {
    "chrom": {"chrom", "chromosome", "chr"},
    "pos": {"pos", "position"},
    "ref": {"ref", "reference", "reference_allele"},
    "alt": {"alt", "alternate", "alternative", "alternate_allele"},
}


def scores_canonical_chrom(value: str) -> str:
    """Normalize only standard chromosome aliases; never liftover coordinates."""
    value = str(value).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    if not value:
        raise ValueError("Empty chromosome name")
    if value.isdigit():
        return "chr" + str(int(value))
    if value.upper() in {"X", "Y"}:
        return "chr" + value.upper()
    if value.upper() in {"M", "MT"}:
        return "chrM"
    return "chr" + value


def _chrom_order(chrom: str) -> tuple:
    name = scores_canonical_chrom(chrom)[3:]
    if name.isdigit():
        return (0, int(name))
    if name in {"X", "Y", "M"}:
        return (0, {"X": 1000001, "Y": 1000002, "M": 1000003}[name])
    return (1, tuple((0, int(p)) if p.isdigit() else (1, p)
                     for p in re.split(r"(\d+)", name)))


def _site_key(chrom: str, pos: int) -> tuple:
    return _chrom_order(chrom), pos


def _build(value: str) -> str:
    compact = str(value).strip().lower().replace("_", "").replace("-", "")
    if compact in {"grch38", "hg38", "38", "b38"}:
        return "GRCh38"
    if compact in {"grch37", "hg19", "37", "b37"}:
        return "GRCh37"
    return str(value).strip()


def _check_build_comment(line: str, source: str, expected_build: str) -> None:
    match = re.match(r"^#+\s*(build|assembly|genome_build|reference_build|reference)\s*[:=]\s*(.+)$", line, re.I)
    if not match:
        return
    value = match.group(2).strip().strip('"\'')
    embedded = re.search(r"GRCh(?:37|38)|hg(?:19|38)", value, re.I)
    if embedded:
        value = embedded.group(0)
    elif match.group(1).lower() == "reference":
        return  # A reference FASTA pathname without an assembly is not a declaration.
    if _build(value) != expected_build:
        raise ValueError(f"Build mismatch in {source}: declared {value!r}; candidates require {expected_build}")


@contextmanager
def _text_input(path: Path):
    opener = gzip.open if path.suffix.lower() in {".gz", ".bgz", ".bgzf"} else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        yield handle


def _header_names(line: str) -> list[str]:
    return line.lstrip("#").rstrip("\r\n").split("\t")


def _normalized_names(names: list[str]) -> list[str]:
    return [x.strip().lstrip("#").lower() for x in names]


def _coordinate_indices(names: list[str]) -> dict[str, int] | None:
    normalized = _normalized_names(names)
    indices = {}
    for key, aliases in _ALIASES.items():
        matches = [i for i, name in enumerate(normalized) if name in aliases]
        if len(matches) > 1:
            raise ValueError(f"Ambiguous score header: multiple columns for {key}")
        if not matches:
            return None
        indices[key] = matches[0]
    return indices


def _score_column_index(names: list[str], chosen: str, score_format: str) -> int:
    normalized = _normalized_names(names)
    requested = chosen.strip().lower()
    if score_format == "avi":
        requested = {"raw": "raw_score", "avi_raw": "raw_score",
                     "avi_phred": "phred"}.get(requested, requested)
    if normalized.count(requested) != 1:
        raise ValueError(f"Score column {chosen!r} must appear exactly once; available columns: {names}")
    return normalized.index(requested)


def _variant_fields(chrom: str, pos: str, ref: str, alt: str, source: str) -> tuple[str, int, str, str]:
    try:
        position = int(pos)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid 1-based position {pos!r} in {source}") from exc
    ref, alt = ref.strip().upper(), alt.strip().upper()
    if position < 1 or ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"} or ref == alt:
        raise ValueError(f"Expected a single-base substitution in {source}: {chrom}:{pos}:{ref}>{alt}")
    return scores_canonical_chrom(chrom), position, ref, alt


class _ScoreTable:
    def __init__(self, path: Path, chosen: str, score_format: str, expected_build: str):
        self.path = path
        self.chosen = chosen
        self.score_format = score_format
        self.expected_build = expected_build
        self.skip_through_line = 0
        self.stats = Counter()
        self.names: list[str] = []
        self.indices: dict[str, int] = {}
        self.score_index = -1
        self.build_indices: list[int] = []
        self._inspect()

    def _inspect(self) -> None:
        if self.score_format == "avi":
            self.names = ["chrom", "pos", "ref", "alt", "raw_score", "phred"]
            self.indices = {name: i for i, name in enumerate(self.names[:4])}
            self.score_index = _score_column_index(self.names, self.chosen, "avi")
        with _text_input(self.path) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                if line.startswith("#"):
                    _check_build_comment(line.strip(), str(self.path), self.expected_build)
                names = _header_names(line)
                indices = _coordinate_indices(names)
                if indices is not None:
                    if self.score_format == "tsv":
                        self.names, self.indices = names, indices
                        self.score_index = _score_column_index(names, self.chosen, "tsv")
                        self.build_indices = [i for i, name in enumerate(_normalized_names(names))
                                              if name in {"build", "assembly", "genome_build", "reference_build"}]
                    else:
                        if indices != {"chrom": 0, "pos": 1, "ref": 2, "alt": 3}:
                            raise ValueError("Explicit AVI format requires chrom,pos,ref,alt in the first four columns")
                        # Official positional format is fixed, even when a header exists.
                        if len(names) < 6:
                            raise ValueError("AVI input requires six columns")
                        normalized = _normalized_names(names)
                        if normalized[4] not in {"raw", "raw_score", "rawscore", "avi_raw"} or normalized[5] not in {"phred", "avi_phred"}:
                            raise ValueError("Explicit AVI header must place raw_score and phred in columns five and six")
                    self.skip_through_line = line_number
                    return
                if line.startswith("#"):
                    continue
                if self.score_format == "tsv":
                    raise ValueError(f"Score TSV requires a header with chrom,pos,ref,alt and {self.chosen!r}: {self.path}. Use score_format='avi' only for official six-column AVI input.")
                self.skip_through_line = line_number - 1
                self.parse(line, f"{self.path}:{line_number}", count=False)
                return
        if self.score_format == "tsv":
            raise ValueError(f"No score TSV header found: {self.path}")

    def parse(self, line: str, source: str, count: bool = True):
        values = line.rstrip("\r\n").split("\t")
        highest = max(*self.indices.values(), self.score_index, *self.build_indices)
        if len(values) <= highest:
            raise ValueError(f"Truncated score row in {source}")
        chrom, pos, ref, alt = _variant_fields(*(values[self.indices[x]] for x in ("chrom", "pos", "ref", "alt")), source)
        for index in self.build_indices:
            if _build(values[index]) != self.expected_build:
                raise ValueError(f"Score row build mismatch in {source}: {values[index]!r}")
        raw = values[self.score_index].strip()
        if raw.lower() in _MISSING:
            score = None
        else:
            try:
                score = float(raw)
            except ValueError as exc:
                raise ValueError(f"Nonnumeric score {raw!r} in {source}") from exc
            if not math.isfinite(score):
                raise ValueError(f"Nonfinite score {raw!r} in {source}")
        if count:
            self.stats["score_rows_read"] += 1
        return chrom, pos, ref, alt, score

    def records(self):
        with _text_input(self.path) as handle:
            for line_number, line in enumerate(handle, 1):
                if line_number <= self.skip_through_line or not line.strip():
                    continue
                if line.startswith("#"):
                    _check_build_comment(line.strip(), str(self.path), self.expected_build)
                    continue
                yield self.parse(line, f"{self.path}:{line_number}")


def _score_sites(records: Iterable, stats: Counter) -> Iterator:
    previous = None
    current = None
    scores = {}
    for chrom, pos, ref, alt, score in records:
        key = _site_key(chrom, pos)
        if previous is not None and key < previous:
            raise ValueError("Score input is not in natural chromosome/position order (chr1,chr2,...,chr10,...,chrX). Supply sorted local scores or an indexed bgzip table with backend='tabix'.")
        if current is not None and key != current:
            yield current, scores
            scores = {}
        allele = (ref, alt)
        if allele in scores:
            if scores[allele] != score:
                raise ValueError(f"Conflicting duplicate score: {chrom}:{pos}:{ref}>{alt}")
            stats["identical_duplicate_score_rows"] += 1
        else:
            scores[allele] = score
        previous = current = key
    if current is not None:
        yield current, scores


class _Candidates:
    def __init__(self, paths: list[Path], expected_build: str, natural_order: bool = True):
        self.paths = paths
        self.expected_build = expected_build
        self.natural_order = natural_order
        self.header: list[str] | None = None

    def records(self):
        previous = None
        previous_chrom = None
        finished_chroms = set()
        for path in self.paths:
            with _text_input(path) as handle:
                header = None
                for line_number, line in enumerate(handle, 1):
                    if not line.strip():
                        continue
                    if line.startswith("#"):
                        _check_build_comment(line.strip(), str(path), self.expected_build)
                    if header is None:
                        names = [x.strip() for x in _header_names(line)]
                        if REQUIRED_CANDIDATE_COLUMNS.issubset(names):
                            if len(set(names)) != len(names):
                                raise ValueError(f"Duplicate candidate header fields in {path}")
                            if set(names) & set(SCORE_COLUMNS):
                                raise ValueError(f"Candidate input already contains reserved score columns: {path}")
                            header = names
                            if self.header is None:
                                self.header = names
                            elif self.header != names:
                                raise ValueError("Candidate shard headers differ; preserve the same column order in every shard")
                            continue
                        if line.startswith("#"):
                            continue
                        missing = sorted(REQUIRED_CANDIDATE_COLUMNS - set(names))
                        raise ValueError(f"Missing candidate columns in {path}: {missing}")
                    if line.startswith("#"):
                        continue
                    values = line.rstrip("\r\n").split("\t")
                    if len(values) != len(header):
                        raise ValueError(f"Candidate row has wrong column count: {path}:{line_number}")
                    row = dict(zip(header, values))
                    chrom, pos, ref, alt = _variant_fields(row["chrom"], row["pos"], row["ref"], row["alt"], f"{path}:{line_number}")
                    row.update(chrom=chrom, pos=str(pos), ref=ref, alt=alt)
                    if row["status"] not in ALLOWED_STATUSES:
                        raise ValueError(f"Unknown candidate status {row['status']!r} in {path}")
                    if row["is_primary"] not in {"0", "1"} or row.get("primary_tie", "0") not in {"0", "1"}:
                        raise ValueError(f"is_primary and primary_tie must be 0 or 1 in {path}")
                    if row["n_alt"] not in {"0", "1", "2", "3"}:
                        raise ValueError(f"n_alt must be 0, 1, 2, or 3 in {path}")
                    for field in ("build", "assembly", "genome_build", "reference_build"):
                        if field in row and _build(row[field]) != self.expected_build:
                            raise ValueError(f"Candidate build mismatch: {row[field]!r} in {path}")
                    key = _site_key(chrom, pos)
                    if self.natural_order and previous is not None and key < previous:
                        raise ValueError("Candidate shards must be supplied in natural chromosome/position order, without overlapping or reversed shards")
                    if not self.natural_order:
                        if chrom != previous_chrom:
                            if chrom in finished_chroms:
                                raise ValueError("Candidate input revisits a chromosome; each chromosome must form one contiguous block")
                            if previous_chrom is not None:
                                finished_chroms.add(previous_chrom)
                        elif previous is not None and key < previous:
                            raise ValueError("Candidate positions must be sorted within each chromosome")
                    previous = key
                    previous_chrom = chrom
                    yield key, row
                if header is None:
                    raise ValueError(f"No candidate header found in {path}")

    def sites(self):
        current = None
        rows = []
        for key, row in self.records():
            if current is not None and key != current:
                yield current, rows
                rows = []
            current = key
            if any(x["ref"] == row["ref"] and x["alt"] == row["alt"] for x in rows):
                raise ValueError(f"Duplicate candidate: {row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}")
            if rows and rows[0]["ref"] != row["ref"]:
                raise ValueError(f"Conflicting candidate reference bases: {row['chrom']}:{row['pos']}")
            rows.append(row)
        if current is not None:
            yield current, rows


def _next(iterator):
    return next(iterator, None)


def _matched_sites(candidate_sites, score_sites):
    current_score = _next(score_sites)
    for key, candidates in candidate_sites:
        while current_score is not None and current_score[0] < key:
            current_score = _next(score_sites)
        yield candidates, current_score[1] if current_score is not None and current_score[0] == key else {}


def _tabix_matched_sites(candidate_sites, table: _ScoreTable, tabix, block_size: int, stats: Counter):
    contigs = {}
    for name in tabix.contigs:
        canonical = scores_canonical_chrom(name)
        if canonical in contigs:
            raise ValueError(f"Ambiguous tabix contig aliases for {canonical}: {contigs[canonical]}, {name}")
        contigs[canonical] = name
    block = None
    score_sites = iter(())
    current_score = None
    for key, candidates in candidate_sites:
        chrom, pos = candidates[0]["chrom"], int(candidates[0]["pos"])
        start = ((pos - 1) // block_size) * block_size
        next_block = (chrom, start)
        if next_block != block:
            # Exhaust the prior queried block to validate all its duplicates/order.
            for _ in score_sites:
                pass
            block = next_block
            current_score = None
            if chrom not in contigs:
                score_sites = iter(())
                stats["blocks_without_score_contig"] += 1
            else:
                stats["tabix_regions_fetched"] += 1
                lines = tabix.fetch(contigs[chrom], start, start + block_size)
                records = (table.parse(line, f"{table.path}:{chrom}:{start + 1}-{start + block_size}")
                           for line in lines if line.strip() and not line.startswith("#"))
                score_sites = iter(_score_sites(records, table.stats))
                current_score = _next(score_sites)
        while current_score is not None and current_score[0] < key:
            current_score = _next(score_sites)
        yield candidates, current_score[1] if current_score is not None and current_score[0] == key else {}
    for _ in score_sites:
        pass


def _number(value: float | None) -> str:
    return "." if value is None else format(value, ".17g")


def _pair(rows: list[dict[str, str]], values: dict[tuple[str, str], float | None], stats: Counter):
    if len(rows) != 3 or any(row["n_alt"] != "2" for row in rows):
        return None
    if any(row.get("primary_tie", "0") == "1" for row in rows):
        stats["pair_sites_excluded_primary_tie"] += 1
        return None
    primary = [r for r in rows if r["is_primary"] == "1" and r["status"] == "observed_pass"]
    observed = [r for r in rows if r["is_primary"] == "0" and r["status"] == "observed_pass"]
    absent = [r for r in rows if r["status"] == "not_listed_in_release" and r["is_primary"] == "0"]
    if len(primary) != 1 or len(observed) != 1 or len(absent) != 1:
        stats["pair_sites_excluded_status"] += 1
        return None
    p, o, u = primary[0], observed[0], absent[0]
    primary_score = values.get((p["ref"], p["alt"]))
    observed_score = values.get((o["ref"], o["alt"]))
    unobserved_score = values.get((u["ref"], u["alt"]))
    if observed_score is None or unobserved_score is None:
        stats["pair_sites_excluded_missing_score"] += 1
        return None
    difference = observed_score - unobserved_score
    if not math.isfinite(difference):
        raise ValueError(f"Score difference overflow at {o['chrom']}:{o['pos']}")
    return dict(
        chrom=o["chrom"], pos=o["pos"], ref=o["ref"], context=o["context"], n_alt="2",
        primary_alt=p["alt"], observed_alt=o["alt"], unobserved_alt=u["alt"],
        primary_score=_number(primary_score), observed_score=_number(observed_score),
        unobserved_score=_number(unobserved_score), delta_observed_minus_unobserved=_number(difference),
        observed_ac=o["ac"], observed_an=o["an"], observed_af=o["af"],
    )


def preflight_scores(
    scores_path,
    score_column: str,
    score_build: str = "GRCh38",
    backend: str = "auto",
    score_format: str = "tsv",
) -> dict[str, Any]:
    """Check local score prerequisites before expensive candidate construction.

    Only the header (or first AVI row) and optional tabix index are opened. No
    score-table scan, candidate read, region fetch or output write is performed.
    The returned runtime objects are reused by ``annotate_scores``; the CLI may
    discard them after this early check. Row values/order remain validated by
    the actual streaming annotation.
    """
    expected_build = _build(score_build)
    if expected_build not in {"GRCh38", "GRCh37"}:
        raise ValueError("score_build must explicitly identify GRCh38 or GRCh37; no liftover is performed")
    if score_format not in {"tsv", "avi"}:
        raise ValueError("score_format must be 'tsv' or 'avi'")
    if score_format == "avi" and expected_build != "GRCh38":
        raise ValueError(f"score_build={score_build!r} does not match official AVI's GRCh38 assembly; no liftover is performed")
    if backend == "streaming":
        backend = "stream"
    if backend not in {"auto", "stream", "tabix"}:
        raise ValueError("backend must be 'auto', 'stream', or 'tabix'")
    score_path = Path(scores_path).resolve()
    if not score_path.is_file():
        raise FileNotFoundError(f"Local input file not found: {score_path}")
    table = _ScoreTable(score_path, score_column, score_format, expected_build)
    has_index = any(Path(str(score_path) + ext).is_file() for ext in (".tbi", ".csi"))
    pysam = None
    if backend != "stream" and has_index:
        try:
            pysam = importlib.import_module("pysam")
        except ImportError:
            pass
    if backend == "tabix" and (pysam is None or not has_index):
        raise ValueError("backend='tabix' requires pysam and an adjacent .tbi or .csi index; use backend='stream' for sorted local input")
    chosen_backend = "tabix" if backend != "stream" and pysam is not None and has_index else "stream"
    if chosen_backend == "tabix":
        # Presence alone does not establish a usable BGZF/index pair. Opening
        # and closing checks it without fetching genomic intervals.
        index = pysam.TabixFile(str(score_path))
        index.close()
    return dict(expected_build=expected_build, score_path=score_path, table=table,
                backend_requested=backend, backend=chosen_backend, pysam=pysam)


def annotate_scores(
    candidate_files: list,
    scores_path,
    out_dir,
    score_column: str,
    score_build: str = "GRCh38",
    score_label: str = "local",
    score_version: str = "unspecified",
    backend: str = "auto",
    block_size: int = 1_000_000,
    score_format: str = "tsv",
) -> dict[str, Any]:
    """Annotate local candidates, returning counts, provenance and output paths.

    ``backend`` is auto, stream (streaming is an alias), or tabix. Auto uses
    pysam/tabix only when pysam and an adjacent .tbi/.csi index are available.
    Otherwise it performs a single bounded-memory sorted merge scan. Existing
    final outputs are replaced only after annotation succeeds. ``pairs.tsv.gz``
    is always written with a header and contains only eligible numeric contrasts;
    an empty file is a valid result. Failure removes this call's partial files.
    """
    if not isinstance(block_size, int) or isinstance(block_size, bool) or block_size < 1:
        raise ValueError("block_size must be a positive integer")
    paths = [Path(p).resolve() for p in candidate_files]
    if not paths:
        raise ValueError("candidate_files must contain at least one local TSV shard")
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Local input file not found: {path}")
    prepared = preflight_scores(scores_path, score_column, score_build, backend, score_format)
    expected_build, score_path, table = (prepared[x] for x in ("expected_build", "score_path", "table"))
    backend, chosen_backend, pysam = (prepared[x] for x in ("backend_requested", "backend", "pysam"))
    destination = Path(out_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    outputs = {
        "candidate_scores": destination / "candidate_scores.tsv.gz",
        "pairs": destination / "pairs.tsv.gz",
        "metadata": destination / "score_annotation.json",
    }
    if set(outputs.values()) & set([*paths, score_path]):
        raise ValueError("Output paths must not overwrite input candidate or score files")
    partials = {key: Path(str(path) + ".partial") for key, path in outputs.items()}
    if any(path.exists() for path in partials.values()):
        raise ValueError("A .partial score output already exists; use another output directory or remove the incomplete files before retrying")
    stats = Counter(total_candidates=0, available=0, absent=0, candidate_sites=0, pairs=0)
    tabix = None
    try:
        candidates = _Candidates(paths, expected_build, natural_order=chosen_backend != "tabix")
        candidate_sites = iter(candidates.sites())
        first = _next(candidate_sites)  # Establish and validate the common header.
        from itertools import chain
        candidate_sites = chain(() if first is None else (first,), candidate_sites)
        if chosen_backend == "tabix":
            tabix = pysam.TabixFile(str(score_path))
            matches = _tabix_matched_sites(candidate_sites, table, tabix, block_size, stats)
        else:
            matches = _matched_sites(candidate_sites, iter(_score_sites(table.records(), table.stats)))
        with gzip.open(partials["candidate_scores"], "wt", encoding="utf-8", newline="", compresslevel=6) as candidate_out, gzip.open(partials["pairs"], "wt", encoding="utf-8", newline="", compresslevel=6) as pair_out:
            annotation_writer = csv.DictWriter(candidate_out, fieldnames=candidates.header + SCORE_COLUMNS, delimiter="\t", lineterminator="\n")
            pair_writer = csv.DictWriter(pair_out, fieldnames=PAIR_COLUMNS, delimiter="\t", lineterminator="\n")
            annotation_writer.writeheader()
            pair_writer.writeheader()
            for rows, values in matches:
                stats["candidate_sites"] += 1
                for row in rows:
                    variant = (row["ref"], row["alt"])
                    score = values.get(variant)
                    reason = "matched" if score is not None else ("missing_score_value" if variant in values else "no_matching_record")
                    status = "available" if score is not None else "absent"
                    annotation_writer.writerow({**row, "local_score": _number(score), "score_status": status, "score_reason": reason})
                    stats["total_candidates"] += 1
                    stats[status] += 1
                    stats[f"reason_{reason}"] += 1
                pair = _pair(rows, values, stats)
                if pair is not None:
                    pair_writer.writerow(pair)
                    stats["pairs"] += 1
        metadata = {
            "schema_version": 1,
            "score_build": expected_build, "score_label": str(score_label),
            "score_version": str(score_version), "score_column": score_column,
            "score_format": score_format, "scores_path": str(score_path),
            "candidate_files": [str(p) for p in paths],
            "backend_requested": backend, "backend": chosen_backend,
            "candidate_order": "contiguous_chromosomes_sorted_positions" if chosen_backend == "tabix" else "natural_chromosome_position",
            "block_size": block_size if chosen_backend == "tabix" else None,
            "score_validation_scope": "queried_coordinate_blocks" if chosen_backend == "tabix" else "streamed_prefix_through_last_candidate",
            "counts": {**dict(stats), **dict(table.stats)},
            "outputs": {key: str(path) for key, path in outputs.items()},
            "warnings": [AVI_CAVEAT, "not_listed_in_release means absent from this released callset, not proven never mutated or extinct. Score distances do not identify ancestry or the allele from which a mutation arose."],
            "interpretation": "Exploratory local annotation and within-site observed-minus-not-listed arithmetic contrast; no P values or causal claims.",
            "sources": ["https://doi.org/10.64898/2026.09.16.26363192", "https://github.com/Ensembl/VEP_plugins/blob/release/116/AVI.pm"],
        }
        with open(partials["metadata"], "wt", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
        if tabix is not None:
            tabix.close()
            tabix = None
        for key in ("candidate_scores", "pairs", "metadata"):
            os.replace(partials[key], outputs[key])
        return metadata
    finally:
        if tabix is not None:
            tabix.close()
        for path in partials.values():
            if path.exists():
                path.unlink()


# ========================================================================
# scan_population
# ========================================================================

# Balanced population architecture and 1KG replication.



# ========================================================================
# Population architecture
# ========================================================================

# 
# Ancestry-balanced population architecture of persistent multi-ALT sites.
# 
# Primary scientific question:
# After standardizing ancestry-group allele-number (AN), are the >=2 persistent
# ALT nucleotides at a global multi-ALT locus present together within the same
# genetic-ancestry group, or does global multi-allelism arise by population
# partitioning?
# 
# AC=1 is not used for the primary claim. Defaults:
#   global per-ALT AC >=5  : primary
#   global per-ALT AC >=10 : strict replication
# Ancestry presence is standardized to common AN by exact hypergeometric
# rarefaction, with P(X>=k)>=0.95 required for "present".
# 
# Terminology: genetic ancestry groups, not race.

import argparse, csv, gzip, json, math, random, hashlib, heapq
from collections import defaultdict, Counter
from functools import lru_cache
from bisect import bisect_right
from pathlib import Path

DEFAULT_GROUPS = "AFR=afr;EAS=eas;EUR=nfe,fin;SAS=sas"

def population_open_text(path, mode="rt"):
    return gzip.open(path, mode, encoding=None if "b" in mode else "utf-8") if str(path).endswith(".gz") else open(path, mode, encoding=None if "b" in mode else "utf-8")

def parse_groups(text):
    out={}
    for item in text.split(";"):
        name, members = item.split("=",1)
        name=name.strip()
        if not name or name in out:raise ValueError('Ancestry group names must be unique and nonempty')
        out[name] = [x.strip().lower() for x in members.split(",") if x.strip()]
        if not out[name]:raise ValueError('Empty ancestry group')
    if len(out)<2: raise ValueError("Need >=2 ancestry groups")
    members=[x for v in out.values() for x in v]
    if len(members)!=len(set(members)):raise ValueError('Ancestry groups must not overlap')
    return out

def logchoose(n,k):
    if k<0 or k>n: return float("-inf")
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

@lru_cache(maxsize=100000)
def hypergeom_ge_k(N,K,n,k):
    """P[X>=k] under the exact hypergeometric distribution; stable upper tail."""
    N=int(N);K=int(K);n=int(n);k=int(k)
    if not (0<=K<=N and 0<=n<=N and k>=1): return None
    lo=max(0,n-(N-K)); hi=min(n,K)
    if hi<k: return 0.0
    if lo>=k: return 1.0
    # Stable survival function avoids cancellation for probabilities far below
    # machine precision (1-CDF via lgamma can return spurious ~1e-8 tails).
    from scipy.stats import hypergeom
    return float(hypergeom.sf(k-1,N,K,n))

def combine_counts(row, members):
    ac = an = 0
    for pop in members:
        try:
            av, nv = float(row["ac_"+pop]), float(row["an_"+pop])
            if not (av.is_integer() and nv.is_integer() and 0 <= av <= nv):
                return None, None
            ac += int(av); an += int(nv)
        except (KeyError, TypeError, ValueError):
            return None, None
    return ac, an

def cosine(a,b):
    pairs=[(x,y) for x,y in zip(a,b) if x is not None and y is not None]
    if len(pairs)<2:return None
    na=math.sqrt(sum(x*x for x,y in pairs)); nb=math.sqrt(sum(y*y for x,y in pairs))
    if not na or not nb:return None
    return sum(x*y for x,y in pairs)/(na*nb)

class population_BlockIndex:
    def __init__(self, path):
        pass
        self.by = defaultdict(list)
        for block in read_blocks(path, "EUR", "38")[0]:
            self.by[block.chrom].append(block)
        self.starts = {c: [b.start for b in blocks] for c, blocks in self.by.items()}

    def locate(self, c, pos1):
        c = str(c).removeprefix("chr"); p = int(pos1)-1
        i = bisect_right(self.starts.get(c, []), p)-1
        blocks = self.by.get(c, [])
        return blocks[i].block_id if i >= 0 and p < blocks[i].end else "outside_partition"


def locus_key(chrom, pos):
    c = str(chrom).removeprefix("chr")
    return (int(c) if c.isdigit() else {"X": 23, "Y": 24, "MT": 25}.get(c, 99), int(pos))

def iter_grouped(path,key=("chrom","pos")):
    with population_open_text(path) as f:
        r=csv.DictReader(f,delimiter="\t")
        current=None;bucket=[]
        for row in r:
            k=locus_key(row[key[0]], row[key[1]])
            if current is not None and k < current: raise ValueError(f"Unsorted allele file: {path}")
            if current is None:current=k
            if k!=current:
                yield current,bucket
                current=k;bucket=[]
            bucket.append(row)
        if current is not None:yield current,bucket

def site_allele_groups(site_path, allele_path):
    ag=iter(iter_grouped(allele_path)); akey=arows=None
    try:akey,arows=next(ag)
    except StopIteration:pass
    previous = None
    with population_open_text(site_path) as f:
        for site in csv.DictReader(f,delimiter="\t"):
            sk=locus_key(site["chrom"], site["pos"])
            if previous is not None and sk <= previous: raise ValueError(f"Unsorted or duplicate sites: {site_path}")
            previous = sk
            while akey is not None and akey<sk:
                try:akey,arows=next(ag)
                except StopIteration:akey=arows=None;break
            if akey != sk: raise ValueError(f"Site missing from allele table: {sk}")
            if any(r.get('ref',site.get('ref')) != site.get('ref',r.get('ref')) for r in arows):raise ValueError('Site/allele REF mismatch')
            yield site, arows

def find_shards(modern_dir):
    root = Path(modern_dir)
    manifest = root / "run.json"
    if manifest.exists():
        run = load_run_metadata(root)
        directories = [Path(p) for p in run["shards"]]
    else:
        directories = sorted((root / "shards").glob("*"))
    out = []
    for d in directories:
        s = d / "multiallelic_sites.tsv.gz"; a = d / "multiallelic_alleles.tsv.gz"
        if not s.is_file() or not a.is_file(): raise ValueError(f"Incomplete modern shard: {d}")
        out.append((s, a))
    if not out: raise ValueError(f"No modern shards under {root}")
    return out

def classify_site(site, alleles, groups, global_ac, target_an, presence_p=.95):
    kept=[]
    for r in alleles:
        if r.get("status")!="observed_pass":continue
        try:ac=int(float(r["ac"]))
        except Exception:continue
        if ac>=global_ac:kept.append(r)
    if len(kept)<2:return None
    present={g:set() for g in groups}
    af_profiles={}
    unknown=False
    for r in kept:
        alt=r["alt"]; vec=[]
        for g,members in groups.items():
            ac,an=combine_counts(r,members)
            if ac is None or an is None or an<target_an:
                vec.append(None); unknown=True; continue
            p=hypergeom_ge_k(an,ac,target_an,global_ac)
            vec.append(ac/an if an else 0.0)
            if p is not None and p>=presence_p:present[g].add(alt)
        af_profiles[alt]=vec
    union=set().union(*present.values()) if present else set()
    within=max((len(x) for x in present.values()),default=0)
    if within>=2: cls="within_group_coexistence"
    elif len(union)>=2 and not unknown: cls="population_partitioned"
    else: cls="unresolved_after_balance"
    # two highest global AC persistent ALTs
    kept.sort(key=lambda r:(-int(float(r["ac"])),r["alt"]))
    sim=cosine(af_profiles[kept[0]["alt"]],af_profiles[kept[1]["alt"]])
    return cls, sim, present, kept

def bootstrap(block_counts, reps=1000, seed=20261024):
    rng=random.Random(seed); blocks=[b for b in block_counts if b != "outside_partition"]
    if len(blocks)<2:return {}
    classes=sorted({c for b in blocks for c in block_counts[b]})
    vals={c:[] for c in classes}
    for _ in range(reps):
        total=Counter()
        for _j in range(len(blocks)):
            b=rng.choice(blocks); total.update(block_counts[b])
        den=sum(total.values())
        if den:
            for c in classes:vals[c].append(total[c]/den)
    out={}
    for c,x in vals.items():
        x.sort(); n=len(x)
        out[c]=(x[int(.025*(n-1))],x[int(.975*(n-1))]) if n else (None,None)
    return out

def run_population(args):
    groups=parse_groups(args.groups); thresholds=[int(x) for x in args.thresholds.split(",")]
    target_ans=[int(x) for x in args.target_an.split(",")]
    bidx=population_BlockIndex(args.ld_bed)
    agg=Counter(); block=defaultdict(lambda: defaultdict(Counter)); sim_bins=Counter(); patterns=Counter()
    example_rows=[]; n_candidate=0
    rep_heaps=defaultdict(list)
    per_rep=max(1,args.replication_sites//max(1,len(thresholds)*3))
    if not thresholds or min(thresholds)<1 or min(target_ans)<1 or not 0 < args.presence_probability <= 1 or args.bootstrap<1:
        raise ValueError("Invalid population thresholds, AN, probability or bootstrap count")
    for site_path,allele_path in find_shards(args.modern_dir):
        print(f"[population] {site_path.parent.name}", flush=True)
        for site,alleles in site_allele_groups(site_path,allele_path):
            chrom,pos=site["chrom"],int(site["pos"]); block_id=bidx.locate(chrom,pos)
            for t in thresholds:
                for N0 in target_ans:
                    result=classify_site(site,alleles,groups,t,N0,args.presence_probability)
                    if not result:continue
                    cls,sim,present,kept=result; key=(t,N0,cls)
                    agg[key]+=1; block[(t,N0)][block_id][cls]+=1; n_candidate+=1
                    if sim is not None:
                        sb=min(9,max(0,int(sim*10))); sim_bins[(t,N0,sb)]+=1
                    pat="|".join(f"{g}:{','.join(sorted(v)) or '-'}" for g,v in present.items())
                    patterns[(t,N0,pat)]+=1
                    if N0==target_ans[0]:
                        rr=dict(chrom=chrom,pos=pos,ref=site.get("ref","NA"),alts=",".join(r["alt"] for r in kept),threshold=t,target_an=N0,gnomad_class=cls)
                        hv=int.from_bytes(hashlib.sha256(f"{chrom}:{pos}:{t}:{args.seed}".encode()).digest()[:8],"big")
                        hp=rep_heaps[(t,cls)]
                        item=(-hv,rr)
                        if len(hp)<per_rep: heapq.heappush(hp,item)
                        elif hv < -hp[0][0]: heapq.heapreplace(hp,item)
                    if args.write_examples and len(example_rows)<args.max_examples:
                        example_rows.append(dict(chrom=chrom,pos=pos,threshold=t,target_an=N0,
                            classification=cls,cosine="NA" if sim is None else sim,pattern=pat,
                            alts=",".join(r["alt"] for r in kept)))
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    summary=[]
    for t in thresholds:
        for N0 in target_ans:
            counts=Counter({c:n for (tt,nn,c),n in agg.items() if tt==t and nn==N0})
            den=sum(counts.values()); ci=bootstrap(block[(t,N0)],args.bootstrap,args.seed+t+N0)
            covered=Counter()
            for bid,values in block[(t,N0)].items():
                if bid!='outside_partition':covered.update(values)
            covered_den=sum(covered.values())
            for c,n in sorted(counts.items()):
                lo,hi=ci.get(c,(None,None))
                summary.append(dict(global_ac_threshold=t,target_an=N0,classification=c,n=covered[c],n_all=n,
                    status="complete" if covered_den else "no_ld_coverage",fraction=covered[c]/covered_den if covered_den else None,fraction_all=n/den if den else None,ci_low=lo,ci_high=hi,bootstrap_blocks=len([b for b in block[(t,N0)] if b != "outside_partition"]),
                    outside_partition=sum(block[(t,N0)].get("outside_partition", {}).values())))
    def write(name,rows,fields=None):
        p=out/name
        if not rows:
            p.write_text("\t".join(fields or ["status"])+"\n");return
        fields=fields or list(rows[0])
        with p.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(rows)
    write("population_balanced_summary.tsv",summary)
    write("population_similarity_bins.tsv",
          [dict(global_ac_threshold=t,target_an=N0,bin_low=b/10,bin_high=(b+1)/10,n=n)
           for (t,N0,b),n in sorted(sim_bins.items())])
    write("population_presence_patterns.tsv",
          [dict(global_ac_threshold=t,target_an=N0,pattern=p,n=n) for (t,N0,p),n in patterns.most_common()])
    if example_rows:write("population_examples.tsv",example_rows)
    replication=[]
    for hp in rep_heaps.values(): replication.extend([item[1] for item in hp])
    replication.sort(key=lambda r:(int(r["threshold"]),r["gnomad_class"],int(r["chrom"]) if str(r["chrom"]).isdigit() else 99,int(r["pos"])))
    write("population_replication_sites.tsv",replication,["chrom","pos","ref","alts","threshold","target_an","gnomad_class"])

    # figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        combos=[(t,n) for t in thresholds for n in target_ans]
        classes=["within_group_coexistence","population_partitioned","unresolved_after_balance"]
        fig,ax=plt.subplots(figsize=(10,5))
        x=list(range(len(combos))); bottom=[0.0]*len(combos)
        for c in classes:
            vals=[]
            for t,N0 in combos:
                row=next((r for r in summary if r["global_ac_threshold"]==t and r["target_an"]==N0 and r["classification"]==c),None)
                vals.append((row["fraction"] or 0) if row else 0)
            ax.bar(x,vals,bottom=bottom,label=c)
            bottom=[a+b for a,b in zip(bottom,vals)]
        ax.set_xticks(x,[f"AC≥{t}\nAN={N0}" for t,N0 in combos])
        ax.set_ylim(0,1);ax.set_ylabel("Fraction of persistent multi-ALT loci")
        ax.set_title("Ancestry-balanced decomposition of persistent multi-ALT loci")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out/"population_partition.png",dpi=220,bbox_inches="tight")
        fig.savefig(out/"population_partition.pdf",bbox_inches="tight")
        plt.close(fig)

        fig,ax=plt.subplots(figsize=(8,5))
        for t,N0 in combos:
            bins=[0]*10
            for (tt,nn,b),v in sim_bins.items():
                if tt==t and nn==N0:bins[b]+=v
            tot=sum(bins)
            if tot:ax.plot([(i+.5)/10 for i in range(10)],[v/tot for v in bins],marker="o",label=f"AC≥{t}, AN={N0}")
        ax.set(xlabel="Cosine similarity of ALT1 vs ALT2 ancestry AF profiles",
               ylabel="Fraction",title="Do persistent ALT alleles occupy the same ancestry background?")
        ax.legend(fontsize=8);fig.tight_layout()
        fig.savefig(out/"population_profile_similarity.png",dpi=220,bbox_inches="tight")
        fig.savefig(out/"population_profile_similarity.pdf",bbox_inches="tight")
        plt.close(fig)
    except ImportError:
        raise RuntimeError("Population figures require matplotlib")

    manifest=dict(
        status="complete", primary_threshold=5,strict_threshold=10,
        bootstrap_replicates=args.bootstrap, replication_sites=len(replication),
        fraction_scope="EUR LD-covered eligible multi-ALT loci; n_all/fraction_all retain all loci; outside_partition is audited",
        class_evaluations=n_candidate, source=str(Path(args.modern_dir).resolve()),
        groups=groups,target_an=target_ans,presence_probability=args.presence_probability,
        method="Exact hypergeometric rarefaction to common AN; block-bootstrap CIs over EUR GRCh38 LD blocks",
        terminology="genetic ancestry groups; not race",
        caveats=[
            "A global gnomAD sites record is not proof that REF and every ALT co-occur in the same individuals.",
            "Ancestry-group summary counts cannot reconstruct haplotypes.",
            "Primary claim requires replication in individual-level 1KG WGS with the same allele definitions.",
            "AC=1 is descriptive emergence/QC only, not primary biological evidence."
        ])
    (out/"population_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest

def main_population(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument("--modern-dir",required=True)
    p.add_argument("--ld-bed",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--thresholds",default="5,10")
    p.add_argument("--target-an",default="5000,10000")
    p.add_argument("--groups",default=DEFAULT_GROUPS)
    p.add_argument("--presence-probability",type=float,default=.95)
    p.add_argument("--bootstrap",type=int,default=1000)
    p.add_argument("--seed",type=int,default=20261024)
    p.add_argument("--write-examples",action="store_true")
    p.add_argument("--max-examples",type=int,default=10000)
    p.add_argument("--replication-sites",type=int,default=10000)
    a=p.parse_args(argv)
    print(json.dumps(run_population(a),indent=2))


# ========================================================================
# Population 1kg
# ========================================================================

# Exact-allele, individual-genotype replication with balanced AN and samples.

import argparse, csv, json, subprocess, tempfile, re
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np

CLASSES = ('within_group_coexistence', 'population_partitioned', 'unresolved_after_balance')


def read_panel(path):
    with population_open_text(path) as f:
        lines = [line.split() for line in f if line.strip() and not line.startswith('##')]
    if not lines: raise ValueError('Empty sample panel')
    header = [s.lower().lstrip('#') for s in lines[0]]
    sample_names = ('sample', 'sample_id', 'iid', 'id')
    pop_names = ('super_pop', 'superpop', 'super_population')
    if any(k in header for k in sample_names):
        si = next(header.index(k) for k in sample_names if k in header)
        pi = next((header.index(k) for k in pop_names if k in header), None)
        if pi is None: raise ValueError('Panel requires exact super_pop labels')
        lines = lines[1:]
    else:
        si, pi = 0, 2 if len(lines[0]) > 2 else 1
    out = {}
    for row in lines:
        if len(row) <= max(si, pi): raise ValueError('Malformed sample panel row')
        if row[si] in out: raise ValueError(f'Duplicate sample in panel: {row[si]}')
        out[row[si]] = row[pi].upper()
    return out


def read_targets(path):
    out = defaultdict(lambda: defaultdict(list)); seen = set()
    with population_open_text(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            c = r['chrom'].removeprefix('chr'); pos = int(r['pos']); ref = r['ref']
            alts = tuple(sorted(set(r['alts'].split(','))))
            key = (c, pos, ref, int(r['threshold']))
            if key in seen: raise ValueError(f'Duplicate target: {key}')
            if len(alts) < 2 or any(len(x) != 1 or x not in 'ACGT' for x in (ref, *alts)):
                raise ValueError(f'Invalid target allele: {key}')
            seen.add(key); out[c][(pos, ref)].append(dict(r, alts=alts))
    return out


def gt_counts(gt, alleles):
    values = re.split(r'[|/]', gt)
    if len(values) != 2 or '.' in values: return []
    try: return [alleles[int(x)] for x in values]
    except (ValueError, IndexError): return []


def classify(counts, targets, groups, threshold, target_an, presence_p=.95):
    rows = []
    for alt in targets:
        row = dict(status='observed_pass', alt=alt, ac=max(threshold, sum(counts.get((g, alt), (0, 0))[0] for g in groups)))
        for g in groups:
            ac, an = counts.get((g, alt), (0, 0)); row['ac_'+g.lower()] = ac; row['an_'+g.lower()] = an
        rows.append(row)
    result = classify_site({}, rows, {g: [g.lower()] for g in groups}, threshold, target_an, presence_p)
    return result[0] if result else CLASSES[2]


def write_table(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t'); w.writeheader(); w.writerows(rows)


def query_chrom(a, chrom, sites, panel, groups):
    """Read only requested positions, retaining missing ALT as unavailable."""
    vcf = Path(a.vcf_dir)/a.vcf_pattern.format(chrom=chrom)
    if not vcf.is_file(): raise FileNotFoundError(vcf)
    samples = subprocess.check_output([a.bcftools, 'query', '-l', str(vcf)], text=True).split()
    absent = set(samples)-set(panel)
    if absent: raise ValueError(f'{len(absent)} VCF samples missing from exact panel: {vcf}')
    keep = [s for s in samples if panel[s] in groups]
    members = {g: np.array([i for i,s in enumerate(keep) if panel[s] == g], dtype=int) for g in groups}
    if any(len(x) == 0 for x in members.values()): raise ValueError('A requested ancestry group has no VCF samples')
    header = subprocess.check_output([a.bcftools, 'view', '-h', str(vcf)], text=True)
    contig = 'chr'+chrom if '##contig=<ID=chr'+chrom+',' in header else chrom
    result = defaultdict(dict); duplicate = set()
    with tempfile.TemporaryDirectory(prefix='allele3-1kg-') as td:
        td = Path(td); sf = td/'samples.txt'; bed = td/'sites.bed'
        sf.write_text('\n'.join(keep)+'\n')
        bed.write_text(''.join(f'{contig}\t{p-1}\t{p}\n' for p in sorted({p for p,r in sites})))
        with tempfile.TemporaryFile() as err1, tempfile.TemporaryFile() as err2:
            view = subprocess.Popen([a.bcftools,'view','-R',str(bed),'-S',str(sf),'-Ou',str(vcf)], stdout=subprocess.PIPE, stderr=err1)
            query = subprocess.Popen([a.bcftools,'query','-f','%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT]\n'], stdin=view.stdout, stdout=subprocess.PIPE, stderr=err2, text=True)
            view.stdout.close()
            try:
                for line in query.stdout:
                    z = line.rstrip('\n').split('\t'); key = (int(z[1]), z[2])
                    if key not in sites or z[4] not in {'PASS', '.'}: continue
                    alleles = [z[2]]+z[3].split(',')
                    wanted = set().union(*(set(t['alts']) for t in sites[key])) & set(alleles[1:])
                    if not wanted: continue
                    calls = [gt_counts(gt, alleles) for gt in z[5:]]
                    if len(calls) != len(keep): raise ValueError('VCF sample cardinality changed')
                    for alt in wanted:
                        if alt in result[key]: duplicate.add((key,alt)); continue
                        result[key][alt] = np.array([c.count(alt) if c else -1 for c in calls], dtype=np.int8)
            finally:
                query.stdout.close()
                rc2 = query.wait(); rc1 = view.wait()
            if rc1 or rc2:
                err1.seek(0); err2.seek(0)
                raise RuntimeError((err1.read()+err2.read()).decode(errors='replace'))
    for key,alt in duplicate: result[key].pop(alt, None)
    return result, members, keep


def run_replication(a):
    panel = read_panel(a.sample_panel); targets = read_targets(a.targets)
    groups = a.groups.split(','); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.bootstrap < 1 or a.subsamples < 1 or a.target_an < 1: raise ValueError('Invalid replication budgets')
    rows=[]; sample_rows=[]; rng=np.random.default_rng(a.seed); sample_order=None
    sampled_classes=defaultdict(lambda: np.zeros((a.subsamples,3),dtype=np.int64))
    for chrom,sites in sorted(targets.items(),key=lambda kv:int(kv[0])):
        print(f'[population-1kg] chr{chrom}: {len(sites)} exact loci', flush=True)
        data, members, order = query_chrom(a,chrom,sites,panel,groups)
        if sample_order is not None and order != sample_order: raise ValueError('1KG sample order/set differs between chromosomes')
        if sample_order is None:
            sample_order=order; m=a.equal_samples or min(map(len,members.values()))
            if m<1 or any(m>len(x) for x in members.values()): raise ValueError('Invalid equal sample count')
            # Same balanced people in each replicate across every locus/chromosome.
            subsets={g:np.array([rng.choice(ids,m,replace=False) for _ in range(a.subsamples)]) for g,ids in members.items()}
        for key, requests in sorted(sites.items()):
            observed = data.get(key,{})
            for t in requests:
                threshold=int(t['threshold'])
                base=dict(chrom=chrom,pos=key[0],ref=key[1],alts=','.join(t['alts']),threshold=threshold,gnomad_class=t['gnomad_class'])
                missing=set(t['alts'])-set(observed)
                if missing:
                    rows.append(base|dict(kg_class='not_assayed_exact_alt',concordant='NA',status='missing_or_ambiguous_exact_alt'))
                    continue
                counts={}
                for alt in t['alts']:
                    for g,ids in members.items():
                        calls=observed[alt][ids]; counts[g,alt]=(int(calls[calls>=0].sum()),int((calls>=0).sum())*2)
                cls=classify(counts,t['alts'],groups,threshold,a.target_an,a.presence_probability)
                rows.append(base|dict(kg_class=cls,concordant=int(cls==t['gnomad_class']),status='measured'))
                # Equal sample sensitivity: fully called diploid samples are required
                # at each ALT; missing calls retain unresolved status in that repeat.
                presence=[]; unknown=np.zeros(a.subsamples,dtype=bool)
                for g in groups:
                    gp=[]
                    for alt in t['alts']:
                        calls=observed[alt][subsets[g]]; callable_=(calls>=0).all(axis=1)
                        unknown |= ~callable_; gp.append(callable_ & (np.maximum(calls,0).sum(axis=1)>=threshold))
                    presence.append(np.array(gp).T)
                within=np.any([x.sum(axis=1)>=2 for x in presence],axis=0)
                union=np.any(presence,axis=0).sum(axis=1)>=2
                codes=np.where(within,0,np.where(union & ~unknown,1,2))
                for i in range(3): sampled_classes[threshold][:,i] += codes==i
                freq=np.bincount(codes,minlength=3)/a.subsamples
                sample_rows.append(base|dict(equal_samples=m,within_fraction=freq[0],partitioned_fraction=freq[1],unresolved_fraction=freq[2]))
    fields=['chrom','pos','ref','alts','threshold','gnomad_class','kg_class','concordant','status']
    write_table(out/'population_1kg_replication.tsv', rows, fields)
    write_table(out/'population_1kg_balanced_samples.tsv',sample_rows,['chrom','pos','ref','alts','threshold','gnomad_class','equal_samples','within_fraction','partitioned_fraction','unresolved_fraction'])
    cross=Counter((r['threshold'],r['gnomad_class'],r['kg_class']) for r in rows)
    write_table(out/'population_1kg_summary.tsv',[dict(threshold=t,gnomad_class=g,kg_class=k,n=n) for (t,g,k),n in sorted(cross.items())],['threshold','gnomad_class','kg_class','n'])
    effects=[]
    for t in sorted({r['threshold'] for r in rows}):
        measured=[r for r in rows if r['threshold']==t and r['status']=='measured']
        resolved=[r for r in measured if r['gnomad_class'] in CLASSES[:2] and r['kg_class'] in CLASSES[:2]]
        byblock=defaultdict(Counter)
        if a.ld_bed:
            BlockIndex = population_BlockIndex
            idx=BlockIndex(a.ld_bed)
            for r in resolved: byblock[idx.locate(r['chrom'],r['pos'])]['agree' if r['concordant'] else 'disagree']+=1
        ci=bootstrap(byblock,a.bootstrap,a.seed).get('agree',(None,None))
        kg=sum(1 if r['kg_class']==CLASSES[1] else -1 for r in resolved)
        gn=sum(1 if r['gnomad_class']==CLASSES[1] else -1 for r in resolved)
        effects.append(dict(threshold=t,n_requested=sum(r['threshold']==t for r in rows),n_exact=len(measured),n_resolved=len(resolved),direction_agreement=(kg*gn>0) if resolved and kg and gn else 'unresolved',class_agreement=sum(r['concordant'] for r in resolved)/len(resolved) if resolved else 'NA',ci_low=ci[0],ci_high=ci[1]))
    write_table(out/'effect_direction_agreement.tsv',effects,['threshold','n_requested','n_exact','n_resolved','direction_agreement','class_agreement','ci_low','ci_high'])
    balanced=[]
    for t,values in sampled_classes.items():
        den=values.sum(axis=1)
        for i,c in enumerate(CLASSES):
            fractions=np.divide(values[:,i],den,out=np.zeros(a.subsamples),where=den>0)
            lo,hi=np.quantile(fractions,[.025,.975]);balanced.append(dict(threshold=t,classification=c,mean=float(fractions.mean()),ci_low=float(lo),ci_high=float(hi)))
    write_table(out/'population_1kg_subsample_summary.tsv',balanced,['threshold','classification','mean','ci_low','ci_high'])
    manifest=dict(status='complete',n_sites=len(rows),groups=groups,target_an=a.target_an,individual_genotypes=True,exact_ref_alt_required=True,thresholds='same as gnomAD target AC (5/10)',equal_samples=a.equal_samples or (m if sample_order else 0),balanced_replicates=a.subsamples,bootstrap_replicates=a.bootstrap,seed=a.seed,sample_panel=str(a.sample_panel),note='Stratified discovery panel; agreement is conditional on this panel, not a genome-wide prevalence estimate. Missing ALTs are unavailable, never zero.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n'); return manifest


def main_replication(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--targets',required=True); p.add_argument('--vcf-dir',required=True); p.add_argument('--vcf-pattern',default='chr{chrom}.vcf.gz')
    p.add_argument('--sample-panel',required=True);p.add_argument('--groups',default='AFR,EAS,EUR,SAS')
    p.add_argument('--target-an',type=int,default=500);p.add_argument('--equal-samples',type=int)
    p.add_argument('--presence-probability',type=float,default=.95);p.add_argument('--bcftools',default='bcftools');p.add_argument('--out',required=True)
    p.add_argument('--bootstrap',type=int,default=1000);p.add_argument('--subsamples',type=int,default=1000);p.add_argument('--seed',type=int,default=20261024);p.add_argument('--ld-bed')
    a=p.parse_args(argv);print(json.dumps(run_replication(a),indent=2))


# ========================================================================
# scan_enrichment
# ========================================================================

# Archaic carrier tracks and scanner-based enrichment.



# ========================================================================
# Archaic tracks
# ========================================================================

# 
# Scalable conversion of IBDmix individual tracts to population carrier-frequency tracks.
# 
# Design:
#   per-person GRCh37 IBDmix tracts
#     -> normalize
#     -> external sort by chrom/reference/sample/start
#     -> merge overlapping tracts within person
#     -> event sweep across people
#     -> carrier-prevalence intervals
#     -> liftOver aggregated intervals once to GRCh38
# 
# IBDmix is unphased: output is carrier prevalence, NOT archaic allele frequency.

import argparse, csv, gzip, json, os, subprocess, tempfile, tarfile, io, hashlib
from pathlib import Path

def enrichment_open_text(path,mode="rt"):
    return gzip.open(path,mode,encoding=None if "b" in mode else "utf-8") if str(path).endswith(".gz") else open(path,mode,encoding=None if "b" in mode else "utf-8")

ALIASES={
 "sample":["ID","id","sample","sample_id","IID","iid"],
 "chrom":["chrom","chr","chromosome","CHR"],
 "start":["start","Start","START"],
 "end":["end","End","END"],
 "anc":["anc","archaic","reference","ref","source"],
 "lod":["LOD","lod","slod","score"],
}
def pick(header,key,required=True):
    for x in ALIASES[key]:
        if x in header:return header.index(x)
    if required:raise ValueError(f"Missing {key}: {header}")
    return None

def tract_files(root,fmt):
    p=Path(root)
    if p.is_file(): return [p]
    chrom_dirs=[p/f'chr{i}' for i in (*range(1,23), 'X', 'Y') if (p/f'chr{i}').is_dir()]
    if chrom_dirs:
        files=[]
        for d in chrom_dirs:
            native=sorted((d/'final').glob('all_archaic_refs*.segments.tsv.gz'))
            if len(native)>1: raise ValueError(f'Ambiguous native tract versions: {d}')
            if fmt=='ibdmix' and native: files+=native
            elif fmt=='ibdmix' and (d/'ibdmix.raw.tar.gz').is_file(): files.append(d/'ibdmix.raw.tar.gz')
            elif fmt=='ibdmix-xlsx' and (d/'ibdmix.tracts.xlsx').is_file(): files.append(d/'ibdmix.tracts.xlsx')
            else: raise ValueError(f'Missing selected tract format: {d}')
        return files
    if fmt=='ibdmix-xlsx': return sorted(p.rglob('ibdmix.tracts.xlsx'))
    preferred=sorted(p.rglob('segments.tsv.gz'))
    if preferred: return preferred
    return sorted(p.rglob('*.segments.tsv.gz'))

def rows_text(path):
    with enrichment_open_text(path) as f:
        r=csv.reader(f,delimiter="\t");header=next(r)
        yield header,r

def normalize_file(path,fmt,out,min_lod,min_length,allowed_samples=None,build=None,coordinates="bed0",chromosomes=None):
    def emit(header, rows, source):
        isamp,ic,is_,ie=pick(header,"sample"),pick(header,"chrom"),pick(header,"start"),pick(header,"end")
        ia,il=pick(header,"anc",False),pick(header,"lod",False)
        for row in rows:
            if row is None:continue
            try:
                sample=str(row[isamp]).strip(); chrom=str(row[ic]).replace("chr","")
                start,end=int(float(row[is_])),int(float(row[ie]))
            except Exception as exc: raise ValueError(f'Malformed tract in {source}: {row}') from exc
            if chromosomes is not None and chrom not in chromosomes:
                raise ValueError(f'Tract chromosome outside selected scope: {chrom}')
            if chrom not in {*(str(i) for i in range(1,23)), 'X', 'Y'}: continue
            if allowed_samples is not None and sample not in allowed_samples: raise ValueError(f'Tract sample outside analyzed panel: {sample}')
            if coordinates=='one-based-closed': start-=1
            if build and 'genome_build' in header and str(row[header.index('genome_build')]).replace('GRCh','').replace('b','') != str(build): raise ValueError('Tract build mismatch')
            if start<0: raise ValueError('Negative tract start')
            if not sample or end<=start or end-start<min_length:continue
            if il is not None and il<len(row) and row[il] not in (None,"","NA"):
                try:
                    if float(row[il])<min_lod:continue
                except Exception:continue
            anc=str(row[ia]).strip() if ia is not None and ia<len(row) and row[ia] not in (None,"") else "archaic"
            out.write(f"{chrom}\t{anc}\t{sample}\t{start}\t{end}\t{source}\n")
    if fmt=="ibdmix-xlsx":
        import openpyxl
        wb=openpyxl.load_workbook(path,read_only=True,data_only=True)
        for ws in wb.worksheets:
            it=ws.iter_rows(values_only=True)
            try:header=list(next(it))
            except StopIteration:continue
            if not any(x in header for x in ALIASES['sample']): continue
            emit(header,it,f"{path}#{ws.title}")
        wb.close()
    elif str(path).endswith('.tar.gz'):
        # Stream the native final member; never unpack arbitrary archive paths.
        found=False
        with tarfile.open(path,'r|gz') as archive:
            for member in archive:
                if member.isfile() and member.name.startswith('final/all_archaic_refs') and member.name.endswith('.segments.tsv.gz'):
                    with archive.extractfile(member) as raw, gzip.GzipFile(fileobj=raw) as gz, io.TextIOWrapper(gz) as text:
                        reader=csv.reader(text,delimiter='\t'); emit(next(reader),reader,f'{path}#{member.name}')
                    found=True; break
        if not found: raise ValueError(f'No native final segment member in {path}')
    else:
        with enrichment_open_text(path) as f:
            r=csv.reader(f,delimiter="\t")
            try:header=next(r)
            except StopIteration:return
            emit(header,r,str(path))

def panel_ids(path):
    ids=set()
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('##'): continue
            cells=line.split(); first=cells[0].lstrip('#')
            if first.lower() in {'sample','sample_id','iid','id','fid'}: continue
            if first in ids: raise ValueError(f'Duplicate analyzed sample: {first}')
            ids.add(first)
    if not ids: raise ValueError('Empty analyzed sample panel')
    return ids


def denominator(args):
    if args.sample_panel: return len(panel_ids(args.sample_panel)), 'sample_panel'
    if args.n_samples and args.n_samples>0: return args.n_samples,'explicit_n_samples'
    raise ValueError('Carrier prevalence requires complete analyzed --sample-panel or positive --n-samples, including zero-tract people')

def external_sort(src,dst,keys):
    cmd=["sort","-T",str(Path(dst).parent),"-S","512M"]+keys+[str(src)]
    with open(dst,"wt") as out:subprocess.run(cmd,stdout=out,check=True)

def build_tracks(args):
    N,den_source=denominator(args)
    files=[Path(p) for p in args.files] if getattr(args, 'files', None) is not None else tract_files(args.dir,args.format)
    if not files:raise ValueError("No IBDmix tract files")
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="archaic-",dir=out) as td:
        td=Path(td);norm=td/"norm.tsv";sorted1=td/"sorted.tsv";events=td/"events.tsv";sorted_events=td/"events.sorted.tsv"
        with norm.open("wt") as h:
            for p in files:
                print(f'[archaic] normalize {p}',flush=True)
                normalize_file(p,args.format,h,args.min_lod,args.min_length,panel_ids(args.sample_panel) if args.sample_panel else None,args.build,getattr(args,'coordinates','bed0'),getattr(args,'chromosomes',None))
        external_sort(norm,sorted1,["-k1,1V","-k2,2","-k3,3","-k4,4n","-k5,5n"])
        with sorted1.open() as f,events.open("wt") as e:
            key=None; lo=hi=None
            def flush(k,lo,hi):
                if k is None:return
                chrom,anc,sample=k;e.write(f"{chrom}\t{anc}\t{lo}\t1\n{chrom}\t{anc}\t{hi}\t-1\n")
            for line in f:
                chrom,anc,sample,s,e0,*_=line.rstrip("\n").split("\t");s=int(s);e0=int(e0);k=(chrom,anc,sample)
                if k!=key:
                    flush(key,lo,hi);key=k;lo,hi=s,e0
                elif s<=hi:hi=max(hi,e0)
                else:flush(key,lo,hi);lo,hi=s,e0
            flush(key,lo,hi)
        external_sort(events,sorted_events,["-k1,1V","-k2,2","-k3,3n"])
        thresholds=[float(x) for x in args.carrier_thresholds.split(",")]
        bed37=out/f"archaic_carrier.GRCh{args.build}.bed"
        with sorted_events.open() as f,bed37.open("wt") as o:
            current=None;count=0;prev=None;delta=0
            def emit(key,a,b,c):
                if key is None or a is None or b<=a or c<=0:return
                chrom,anc=key;frac=c/N
                if c>N: raise ValueError("Carrier count exceeds analyzed denominator")
                o.write(f"chr{chrom}\t{a}\t{b}\tarchaic:{anc}:any\t{frac:.10g}\n")
                for t in thresholds:
                    if frac>=t:o.write(f"chr{chrom}\t{a}\t{b}\tarchaic:{anc}:carrier_ge_{str(t).replace('.','p')}\t{frac:.10g}\n")
            for line in f:
                chrom,anc,p,d=line.rstrip("\n").split("\t");p=int(p);d=int(d);k=(chrom,anc)
                if k!=current:
                    current=k;count=0;prev=None
                if prev is not None and p>prev:emit(current,prev,p,count)
                count+=d;prev=p

    result_bed=bed37
    lift=None
    if args.target_build!=args.build:
        if not args.chain:raise ValueError("--chain required across builds")
        mapped=out/f"archaic_carrier.GRCh{args.target_build}.bed";unmapped=out/"liftover.unmapped.bed"
        cmd=[args.liftover,"-minMatch=0.95","-bedPlus=3",str(bed37),args.chain,str(mapped),str(unmapped)]
        subprocess.run(cmd,check=True)
        result_bed=mapped
        def bp(path):
            n=b=0
            with open(path) as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):continue
                    a=line.split();n+=1;b+=max(0,int(a[2])-int(a[1]))
            return n,b
        ni,bi=bp(bed37);nm,bm=bp(mapped);nu,bu=bp(unmapped)
        lift={"input_intervals":ni,"input_bp":bi,"mapped_intervals":nm,"mapped_bp":bm,
              "mapped_bp_fraction":(bi-bu)/bi if bi else None,"mapped_source_bp":bi-bu,"target_span_change_bp":bm-(bi-bu),"target_to_source_bp_ratio":bm/bi if bi else None,"unmapped_intervals":nu,"unmapped_bp":bu,"counting":"label-weighted interval bp","command":cmd}
    manifest={"status":"complete","format":args.format,"files":len(files),"denominator":N,"denominator_source":den_source,
              "coordinate_input":f"GRCh{args.build} BED0 half-open","coordinate_output":f"GRCh{args.target_build}",
              "carrier_thresholds":thresholds,"output":str(result_bed),"liftover":lift,
              "interpretation":"IBDmix is unphased; score is fraction of analyzed people carrying a tract at the locus, not archaic allele frequency."}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest

def main_tracks(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument("--dir",required=True)
    p.add_argument("--format",choices=["ibdmix","ibdmix-xlsx"],default="ibdmix")
    p.add_argument("--build",choices=["37","38"],default="37")
    p.add_argument("--target-build",choices=["37","38"],default="38")
    p.add_argument("--sample-panel")
    p.add_argument("--n-samples",type=int)
    p.add_argument("--min-lod",type=float,default=4)
    p.add_argument("--min-length",type=int,default=50000)
    p.add_argument("--carrier-thresholds",default="0.01,0.05,0.10,0.20")
    p.add_argument("--coordinates",choices=["bed0","one-based-closed"],default="bed0",help="GU finalized native segments and XLSX are BED0")
    p.add_argument("--chain")
    p.add_argument("--liftover",default="liftOver")
    p.add_argument("--out",required=True)
    a=p.parse_args(argv);print(json.dumps(build_tracks(a),indent=2))


# ========================================================================
# Enrich
# ========================================================================

# 
# Generic interval-enrichment components for integration with the existing gnomAD scanner.
# 
# Primary inference thresholds are AC>=5 and AC>=10. AC=1 is QC/emergence only.
# 
# The existing ArchitectureAccumulator.observe() hook in f/scan.py supplies
# exposure counts after the original site QC and threshold_counts calculation.
# This module normalizes annotations and performs inference from those counts.

import argparse, csv, gzip, json, math, random
from collections import defaultdict, Counter
from dataclasses import dataclass
from pathlib import Path

def canon(c): return str(c).removeprefix("chr")

@dataclass(frozen=True)
class Iv:
    chrom:str;start:int;end:int;label:str;score:float|None=None

class Sweep:
    def __init__(self,intervals):
        self.by=defaultdict(list)
        for x in intervals:self.by[x.chrom].append(x)
        for c in self.by:self.by[c].sort(key=lambda x:x.start)
        self.c=None;self.i=0;self.active=[];self.last={}
    def locate(self,c,pos1):
        c=canon(c);p=int(pos1)-1
        if p<self.last.get(c,-1):raise ValueError("Sorted site order required")
        self.last[c]=p
        a=self.by.get(c,[])
        if c!=self.c:self.c=c;self.i=0;self.active=[]
        self.active=[x for x in self.active if x.end>p]
        while self.i<len(a) and a[self.i].start<=p:
            if a[self.i].end>p:self.active.append(a[self.i])
            self.i+=1
        return tuple(x for x in self.active if x.start<=p<x.end)

def read_label_bed(path, chromosome=None):
    out=[]
    opener=gzip.open if str(path).endswith(".gz") else open
    with opener(path,"rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):continue
            a=line.split()
            if a[0].lower()=="chrom":continue
            if len(a)<4:raise ValueError("BED requires chrom,start,end,label")
            if chromosome and canon(a[0]) != canon(chromosome): continue
            if int(a[1]) < 0 or int(a[2]) <= int(a[1]): raise ValueError('Invalid BED interval')
            out.append(Iv(canon(a[0]),int(a[1]),int(a[2]),a[3],float(a[4]) if len(a)>4 and a[4]!="NA" else None))
    return out


def bootstrap_statistics(blocks, reps=2000, seed=20261024):
    """Cluster bootstrap of O/E; centered two-sided bootstrap test of O/E=1."""
    import numpy as np
    if len(blocks)<2: return None,None,None
    values=np.asarray(list(blocks.values()),dtype=float)
    total=values.sum(axis=0)
    if total[1]<=0: return None,None,None
    point=total[0]/total[1]; rng=np.random.default_rng(seed); draws=[]
    for start in range(0,reps,100):
        ids=rng.integers(0,len(values),size=(min(100,reps-start),len(values)))
        sums=values[ids].sum(axis=1); valid=sums[:,1]>0
        draws.extend((sums[valid,0]/sums[valid,1]).tolist())
    if len(draws)<reps*.95: return None,None,None
    x=np.asarray(draws);lo,hi=np.quantile(x,[.025,.975])
    p=(1+int((np.abs(x-point)>=abs(point-1)).sum()))/(len(x)+1)
    return float(lo),float(hi),float(p)


def bootstrap_oe(block_oe,reps=2000,seed=20261024):
    return bootstrap_statistics(block_oe,reps,seed)[:2]


def bh_fdr(values):
    order=sorted(range(len(values)),key=lambda i:values[i]); out=[1.]*len(values);bound=1.
    for rank,i in reversed(list(enumerate(order,1))):
        bound=min(bound,values[i]*len(values)/rank);out[i]=bound
    return out


def summarize_exposure(global_path, exposure_path, out, reps=2000, seed=20261024, min_training_sites=1):
    """LOCO background and exposure-specific EUR block counts from the QC hook."""
    if reps<2000: raise ValueError('At least 2000 bootstrap replicates required')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    def read(paths):
        for path in paths if isinstance(paths, (list, tuple)) else [paths]:
            with gzip.open(path,'rt') as f: yield from csv.DictReader(f,delimiter='\t')
    dims=('context','af_bin','an_bin','sensitivity','population_mode','threshold')
    total=defaultdict(lambda:[0,0]);chroms=defaultdict(lambda:[0,0])
    for r in read(global_path):
        if r['threshold'] not in {'5','10'} or r['population_mode']!='all':continue
        st=tuple(r[x] for x in dims)
        for table,key in ((total,st),(chroms,(r['chrom'],*st))):
            table[key][0]+=int(r['n_sites']);table[key][1]+=int(r['n_multi'])
    # Number of keys is labels x sensitivities x LD blocks, independent of loci.
    blocks=defaultdict(lambda:defaultdict(lambda:[0,0.,0,0]))
    outside=Counter()
    for r in read(exposure_path):
        if r['threshold'] not in {'5','10'} or r['population_mode']!='all':continue
        key=(r['label'],r['threshold'],r['sensitivity'])
        if r['block_id']=='outside_partition': outside[key]+=int(r['n_sites']);continue
        st=tuple(r[x] for x in dims);n=int(r['n_sites']);observed=int(r['n_multi'])
        tn=total[st][0]-chroms[(r['chrom'],*st)][0]
        tm=total[st][1]-chroms[(r['chrom'],*st)][1]
        z=blocks[key][r['block_id']];z[0]+=observed;z[2]+=n
        if tn<min_training_sites: z[3]+=n
        else: z[1]+=n*(tm+.5)/(tn+1)
    fields=['label','threshold','sensitivity','observed_multi','expected_multi','oe','ci_low','ci_high','bootstrap_p','fdr','block_count','n_eligible','n_without_loco','outside_partition_sites','status']
    rows=[]
    with gzip.open(out/'enrichment_block_counts.tsv.gz','wt',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['label','threshold','sensitivity','block_id','observed_multi','expected_multi','n_eligible','n_without_loco'])
        for key,b in sorted(blocks.items()):
            for bid,vals in sorted(b.items()):w.writerow([*key,bid,*vals])
            o=sum(x[0] for x in b.values());e=sum(x[1] for x in b.values());n=sum(x[2] for x in b.values());missing=sum(x[3] for x in b.values())
            ci=(None,None,None) if missing else bootstrap_statistics({bid:vals[:2] for bid,vals in b.items()},reps,seed)
            status='insufficient_loco_training' if missing else 'insufficient_blocks' if len(b)<2 else 'zero_expected' if not e else 'complete'
            rows.append(dict(zip(('label','threshold','sensitivity'),key))|dict(observed_multi=o,expected_multi=e if not missing else 'NA',oe=o/e if e and not missing else 'NA',ci_low=ci[0],ci_high=ci[1],bootstrap_p=ci[2],fdr=None,block_count=len(b),n_eligible=n,n_without_loco=missing,outside_partition_sites=outside[key],status=status))
    groups=defaultdict(list)
    for r in rows:
        if r['bootstrap_p'] is not None:groups[r['threshold'],r['sensitivity']].append(r)
    for group in groups.values():
        for r,q in zip(group,bh_fdr([r['bootstrap_p'] for r in group])):r['fdr']=q
    with open(out/'enrichment_summary.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
    manifest=dict(status='complete',bootstrap_replicates=reps,seed=seed,primary_thresholds=[5,10],ld_ancestry='EUR',inference='LOCO context/primary-AF/AN expected; percentile LD-block bootstrap CI; centered two-sided bootstrap P; BH FDR within threshold/sensitivity',missing_loco='NA if any exposed stratum is unsupported',outside_partition='excluded from inference; count retained',n_tests=len(rows),difficult_region_complete=any(r['sensitivity']=='exclude_combined' for r in rows))
    (out/'enrichment_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def main_enrichment(argv=None):
    p=argparse.ArgumentParser(description='Normalize annotation and run the original scanner with an enrichment sidecar')
    p.add_argument('--dir',required=True);p.add_argument('--type',choices=['herv','archaic','bed'],required=True)
    p.add_argument('--format',choices=['rmsk-joined','ibdmix','ibdmix-xlsx','bed'],required=True);p.add_argument('--out',required=True)
    p.add_argument('--vcf-dir');p.add_argument('--vcf');p.add_argument('--region');p.add_argument('--fasta',default='/mnt/f/gen/1hgp/38/GRCH38.fasta.gz')
    p.add_argument('--ld-block-dir',default='/mnt/f/LD/block');p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--sample-panel');p.add_argument('--build',default='37');p.add_argument('--chain');p.add_argument('--liftover',default='liftOver')
    p.add_argument('--segdup-bed');p.add_argument('--lcr-bed');p.add_argument('--low-mappability-bed')
    a=p.parse_args(argv)
    if not (a.vcf_dir or a.vcf):p.error('Supply --vcf-dir or --vcf; enrichment uses original eligible loci from the scanner')
    import scan as cli
    out=Path(a.out); source=Path(a.dir)
    extra=[]
    if a.type=='herv' and a.format=='rmsk-joined':
        if source.is_dir():
            candidates=[source/name for name in ('rmskJoinedCurrent.txt.gz','rmskJoinedCurrent.txt') if (source/name).is_file()]
            if len(candidates)!=1:raise ValueError('Directory must contain exactly one rmskJoinedCurrent table; otherwise supply the file')
            source=candidates[0]
        extra=['--erv-joined',str(source),'--architecture-annotation-cache',str(out/'annotations')]
    elif a.type=='archaic' and a.format in {'ibdmix','ibdmix-xlsx'}:
        pass
        from types import SimpleNamespace
        prepared=build_tracks(SimpleNamespace(dir=a.dir,format=a.format,out=str(out/'annotations'),sample_panel=a.sample_panel,n_samples=None,build=a.build,target_build='38',chain=a.chain,liftover=a.liftover,min_lod=4.,min_length=50000,carrier_thresholds='0.01,0.05,0.10,0.20',coordinates='bed0'))
        extra=['--enrichment-track',a.type+'='+prepared['output']]
    elif a.format=='bed':
        read_label_bed(source);extra=['--enrichment-track',a.type+'='+str(source)]
    else:p.error('Unsupported annotation type/format combination')
    for name in ('segdup','lcr','low-mappability'):
        path=getattr(a,name.replace('-','_')+'_bed')
        if path:extra += ['--architecture-'+name+'-bed',path]
    command=['run','--out',str(out/'scan'),'--fasta',a.fasta,'--with-architecture','--ld-block-dir',a.ld_block_dir,'--jobs',str(a.jobs)]+extra
    command += ['--vcf',a.vcf] if a.vcf else ['--vcf-dir',a.vcf_dir]
    if a.region: command += ['--region',a.region]
    return cli.main_modern(command)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'mappability':
        sys.argv.pop(1)
        main_mappability()
    else:
        raise SystemExit(main_modern())
