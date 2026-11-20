#!/usr/bin/env python3
"""Observed A/C/G/T occupancy, counts, functional strata and ancient trajectories."""
from __future__ import annotations
import importlib
import sys

MODULES = {
    "hotspots": ("hotspot_pipeline", "Genome-wide window discovery, frozen post-discovery characterization and report"),
    "genome-windows": ("genome_windows", "Genome-wide multi-ALT windows from completed count tiles or a labelled legacy summary"),
    "density": ("allele_density", "SDS-inspired continuous multi-ALT density scores from existing count tiles; no recent-selection claim"),
    "hotspot-characterize": ("hotspot_pipeline", "Annotate frozen hotspots and matched windows after discovery; prepare exact target panels"),
    "window-selection": ("hotspot_pipeline", "Summarize exact published ancient selection contrasts in a frozen window target panel"),
    "data": ("data_sources", "Source registry, bounded schema checks and GRCh37/38 audit"),
    "allofus": ("allofus_sites", "Explicit All of Us summary schema and targeted WGS evidence"),
    "allofus-cache": ("allofus_shards", "Checkpointed wide All of Us window cache and indexed candidate query"),
    "mfi": ("ukb_mfi", "Match UKB MFI to actual PVAR alleles and produce an INFO quality panel"),
    "mapped-targets": ("chronology", "Export base-resolved GRCh37 targets from a verified existing time mapping"),
    "panel": ("panel_targets", "Count exact candidate nucleotides from a quality-audited individual PGEN panel"),
    "atlas-allofus": ("atlas_scores", "Exact merge of verified All of Us AF and Atlas predictions"),
    "atlas": ("atlas_scores", "Query local GRCh38 AlphaGenome Atlas Tabix scores for all three theoretical ALT bases"),
    "opportunity": ("allele_opportunity", "Join all four nucleotide labels to reported ALT counts and Atlas predictions; no missing-as-extinct inference"),
    "opportunity-model": ("allele_opportunity", "Within-site conditional association of molecular predictors with count support; requires joint callable WGS counts"),
    "functional": ("allele_functional", "Matched signed molecular profiles, molecular breadth and measured current-rank geometry; rejects AVI/SHAP proxies"),
    "joint-counts": ("joint_counts", "Stream native joint-call WGS GTs into actual A/C/G/T counts, preserving biallelic controls"),
    "regions": ("joint_counts", "Attach same-build regulatory/HERV/archaic BED membership to prespecified sites"),
    "strata": ("allele_functional", "Observed three/four-base fractions by AVI and region with explicit denominators and block uncertainty"),
    "base_simulation": ("base_simulation", "Four-allele mechanism simulations; not fitted human evolutionary inference"),
    "review": ("study_summary", "Re-summarize existing scan/hotspots/association result status without rescanning genomes"),
    "decompose": ("archaic_alleles", "Exact head/body/tail contributions to total ERV excess, using existing partitions without rescanning"),
    "balance": ("joint_counts", "Actual base counts, frequencies and exact sample-size rarefaction"),
    "space": ("allele_space", "Allele occupancy and callable frequencies in fixed PC1/PC2/time space"),
    "allele-time": ("chronology", "Exact ancient-panel trajectories for an explicit genome-window or HERV target catalog"),
    "domains": ("herv", "Curated ERV-domain partition enrichment using indexed eligible-site queries"),
    "cohort": ("herv", "Balanced personal HERV census planning and completed-result audit"),
    "archaic-alleles": ("archaic_alleles", "Within-person target-allele and IBDmix-tract co-occurrence"),
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"help", "-h", "--help"}:
        print("Usage: ./allele3.sh base_occupancy MODULE [OPTIONS]\n")
        for name, (_, description) in MODULES.items():
            print(f"  {name:17s} {description}")
        print("\nUse MODULE --help for exact inputs. No module reconstructs mutation birth dates.")
        return 0
    name = argv.pop(0)
    if name == 'dynamics': name = 'base_simulation'
    if name not in MODULES:
        raise SystemExit("Unknown base_occupancy module: " + name)
    module = importlib.import_module(MODULES[name][0])
    entry = {
        'allele-time': 'allele_main', 'domains': 'domains_main', 'cohort': 'cohort_main',
        'hotspot-characterize': 'characterize_main',
        'window-selection': 'selection_main',
        'atlas-allofus': 'allofus_main',
        'opportunity-model': 'model_main',
        'regions': 'regions_main',
        'balance': 'balance_main',
        'strata': 'strata_main',
        'mapped-targets': 'mapped_main',
        'decompose': 'decompose_main',
    }.get(name, 'main')
    result = getattr(module, entry)(argv)
    return result if isinstance(result, int) else 0


if __name__ == "__main__":
    raise SystemExit(main())
