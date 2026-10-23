"""Small, fixed-cohort summaries of ALT-count spectra from verified scan shards.

Older caches remain usable but cannot manufacture a joint AC distribution that
was not saved. Counts from partial shard coverage are explicitly labelled.
"""
from __future__ import annotations

from collections import Counter
import csv
import json
from pathlib import Path

from core import dump_json, tsv_writer

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
