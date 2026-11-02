#!/usr/bin/env python3
"""Optional analyses of third-base evidence; existing primary stages stay explicit."""
from __future__ import annotations
import importlib
import sys

MODULES = {
    "data": ("data_sources", "Source registry, bounded schema checks and GRCh37/38 audit"),
    "allofus": ("allofus_sites", "Explicit All of Us summary schema and targeted WGS evidence"),
    "allofus-cache": ("allofus_shards", "Checkpointed wide All of Us window cache and indexed candidate query"),
    "mfi": ("ukb_mfi", "Match UKB MFI to actual PVAR alleles and produce an INFO quality panel"),
    "mapped-targets": ("mapped_targets", "Export base-resolved GRCh37 targets from a verified existing time mapping"),
    "panel": ("panel_targets", "Count exact candidate nucleotides from a quality-audited individual PGEN panel"),
    "atlas-allofus": ("atlas_allofus", "Exact merge of verified All of Us AF and Atlas predictions"),
    "atlas": ("atlas_scores", "Query local GRCh38 AlphaGenome Atlas Tabix scores for all three theoretical ALT bases"),
    "balance": ("allele_balance", "Actual base counts, frequencies and exact sample-size rarefaction"),
    "space": ("allele_space", "Allele occupancy and callable frequencies in fixed PC1/PC2/time space"),
    "domains": ("herv_domains", "Curated ERV-domain partition enrichment using indexed eligible-site queries"),
    "cohort": ("herv_cohort", "Balanced personal HERV census planning and completed-result audit"),
    "archaic-alleles": ("archaic_alleles", "Within-person target-allele and IBDmix-tract co-occurrence"),
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"help", "-h", "--help"}:
        print("Usage: ./allele3.sh research MODULE [OPTIONS]\n")
        for name, (_, description) in MODULES.items():
            print(f"  {name:17s} {description}")
        print("\nUse MODULE --help for exact inputs. No module reconstructs mutation birth dates.")
        return 0
    name = argv.pop(0)
    if name not in MODULES:
        raise SystemExit("Unknown research module: " + name)
    module = importlib.import_module(MODULES[name][0])
    result = module.main(argv)
    return result if isinstance(result, int) else 0


if __name__ == "__main__":
    raise SystemExit(main())
