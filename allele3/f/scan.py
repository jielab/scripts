"""One-pass summary-first scan of split gnomAD sites VCFs.

Only genuinely observed, PASS, counted SNV alternatives enter occupancy.
Reference counts, a reference frequency, joint entropy and allele age are NOT
inferred. Failed/missing alternatives remain distinct from unlisted candidates.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import itertools
import math
from pathlib import Path
import re
import time

from core import VERSION, dump_json, fmt, tsv_writer
from gnomad_io import (FastaContext, allele_info_values, canonical_chrom,
                       chrom_key, iterate_records, parse_region, read_header,
                       validate_count_header, _file_chrom)
from intervals import IntervalMask

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
        from architecture import from_options
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
