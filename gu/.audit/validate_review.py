from pathlib import Path
import csv, hashlib, json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
from phyml_risk import summarize_risk

root = Path('/mnt/d/analysis/gu/threshold_review/20260906/legacy_compatible')
source = Path('/mnt/d/analysis/gu/phyml/1kg')
marker_map = Path(__file__).resolve().parents[1] / 'f/phyml_markers.tsv'
expected_risk = {'rs12405323':737, 'rs9349320':55, 'rs1041791':978, 'rs417915':1250, 'rs6525167':837}
checks = dict(stage1_unchanged_tables=0, loci=0, candidate_trees=0, phyml_version='3.3.20260528', risk_counts={})
for d in sorted(root.glob('chr*')):
    with (d/'final/evidence_loci.tsv').open() as f:
        locus = next(csv.DictReader(f, delimiter='\t'))
    assert locus['candidate_tree_status'] not in ('pending', 'not_run'), d
    checks['loci'] += 1
    for p in (source/d.name/'final').glob('region*.gz'):
        assert hashlib.sha256(p.read_bytes()).digest() == hashlib.sha256((d/'final'/p.name).read_bytes()).digest(), p
        checks['stage1_unchanged_tables'] += 1
    for phy in (d/'loci').glob('haplotypes.evidence.*.phy'):
        stats = Path(str(phy)+'_phyml_stats.txt').read_text()
        assert '3.3.20260528' in stats and 'HKY85' in stats, phy
        assert Path(str(phy)+'_phyml_tree.txt').stat().st_size > 0, phy
        with Path(str(phy)+'_phyml_boot_trees.txt').open() as f:
            assert sum(bool(line.strip()) for line in f) == 100, phy
        assert Path(str(phy)+'_phyml_tree.png').stat().st_size > 0, phy
        checks['candidate_trees'] += 1
    risks = summarize_risk(d, marker_map)
    if locus['locus_id'] in expected_risk:
        assert risks, d
        for row in risks:
            assert row['n_risk_copies'] == expected_risk[locus['locus_id']], row
            assert row['n_candidate_risk_copies'] <= row['n_candidate_anchor_called_copies'] <= row['n_candidate_copies'], row
        checks['risk_counts'][locus['locus_id']] = risks[0]['n_risk_copies']
assert checks['loci'] == 8 and checks['stage1_unchanged_tables'] == 64 and checks['candidate_trees'] == 8, checks
(root.parent/'validation.json').write_text(json.dumps(checks, indent=2)+'\n')
print(json.dumps(checks, indent=2))
