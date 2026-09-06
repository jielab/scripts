"""Post-discovery GWAS allele concordance; never used to select candidates."""
from __future__ import annotations
import argparse
from pathlib import Path
from phyml_evidence import read_tsv, write_tsv


def summarize_risk(out: Path, marker_map: Path):
    markers = {(r['genome_build'], r['marker_id']): r for r in read_tsv(marker_map)}
    loci = read_tsv(out / 'final/evidence_loci.tsv')
    copies = read_tsv(out / 'final/anchor_copies.tsv')
    haps = read_tsv(out / 'final/evidence_haplotypes.tsv')
    trees = {(r['locus_id'], r['expected_lineage']): r
             for r in read_tsv(out / 'final/evidence_trees.tsv')}
    results = []
    for locus in loci:
        lid = locus['locus_id']
        marker = markers.get((locus['genome_build'], locus.get('name', lid)))
        if marker is None:
            continue
        risk = marker['risk_allele']
        own = [r for r in copies if r['locus_id'] == lid]
        risk_copies = {(r['sample'], r['haplotype']) for r in own if r['allele'] == risk}
        for lineage in sorted({r['diagnostic_lineage'] for r in haps if r['locus_id'] == lid}):
            tree = trees.get((lid, lineage), {})
            candidate_haps = {r['hap_id'] for r in haps if r['locus_id'] == lid
                              and r['diagnostic_lineage'] == lineage and r['diagnostic_candidate_pass'] == '1'}
            candidate_copies = {(r['sample'], r['haplotype']) for r in own if r['hap_id'] in candidate_haps}
            matched = len(candidate_copies & risk_copies)
            results.append(dict(locus_id=lid, genome_build=locus['genome_build'],
                lineage=lineage, marker_id=marker['marker_id'], risk_allele=risk,
                trait=marker['trait'], n_risk_copies=len(risk_copies) if own else '',
                n_candidate_haplotype_types=len(candidate_haps),
                n_candidate_copies=sum(int(r.get('n', 0)) for r in haps if r['locus_id'] == lid
                                       and r['diagnostic_lineage'] == lineage and r['hap_id'] in candidate_haps),
                n_candidate_anchor_called_copies=len(candidate_copies),
                n_candidate_risk_copies=matched if own else '',
                candidate_risk_fraction=matched / len(candidate_copies) if candidate_copies else '',
                status='not_evaluable' if not own else 'no_candidate' if not candidate_haps else
                       'risk_concordant' if candidate_copies and matched == len(candidate_copies) else
                       'mixed_or_risk_discordant',
                lineage_call=('no_candidate' if not candidate_haps else 'supported_candidate'
                              if tree.get('candidate_clade_pass') == '1' else 'candidate_sequence_only'),
                lineage_tree_bootstrap=tree.get('candidate_clade_bootstrap', ''),
                lineage_tree_pass=tree.get('candidate_clade_pass', '0'),
                primary_locus_lineage=locus['evidence_lineage'],
                primary_locus_call=locus['introgression_call'],
                interpretation='post_discovery_allele_concordance;not_GWAS_colocalization_or_causality',
                source=marker['source']))
    fields = list(results[0]) if results else ['locus_id', 'status']
    write_tsv(out / 'final/gwas_risk_links.tsv', fields, results)
    return results


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--marker-map', type=Path, default=Path(__file__).with_name('phyml_markers.tsv'))
    a = p.parse_args()
    summarize_risk(a.out, a.marker_map)
