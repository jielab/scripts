#!/usr/bin/env python3
"""Outer join UKB assay targets and CKB download identifiers; retain source audit."""
import csv
import json
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse, unquote

ROOT = Path('/mnt/e/gwas/prot')
ALIASES = {
    'CERT': ('CERT1', 'UKB CERT → CERT1；UKB HGNC.symbol=CERT1，UniProt Q9Y5P4'),
    'WARS': ('WARS1', 'UKB WARS → WARS1；UKB HGNC.symbol=WARS1，UniProt P23381'),
    'NTproBNP': ('NT-proBNP', 'UKB NTproBNP → NT-proBNP；同一 N 端前体片段，保留片段名称，不与 NPPB 靶标合并'),
    'MICB_MICA': ('MICA_MICB', 'UKB MICB_MICA → MICA_MICB；同一组合靶标，仅名称顺序不同，UniProt Q29980/Q29983'),
    'SPACA5_SPACA5B': ('SPACA5', 'UKB SPACA5_SPACA5B → SPACA5（CKB名称）；UniProt Q96QH8 对应 SPACA5/SPACA5B，共享蛋白条目；按蛋白靶标合并，不代表检测特异性相同；https://www.uniprot.org/uniprotkb/Q96QH8/entry'),
}
OTHER_TRAITS = {
    f'{trait}_bmi_{adj}' for trait in ('sbp', 'dbp', 'map', 'pp') for adj in ('adj', 'unadj')
} | {'standing_height', 'sitting_height', 'leg_length', 'shr'}


def write_table(name, rows, fields, delimiter='\t'):
    with (ROOT / name).open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def main():
    with (ROOT / 'ppp_3k.38.tsv').open() as handle:
        ukb_rows = list(csv.DictReader(handle, delimiter='\t'))
    ukb_names = Counter(row['Assay'] for row in ukb_rows)
    urls = [line.strip() for line in Path('/mnt/e/gwas/ckb/ckb.url.txt').read_text().splitlines() if line.strip()]
    ckb_names = [unquote(urlparse(url).path.rstrip('/').rsplit('/', 1)[-1]) for url in urls]
    assert len(ckb_names) == len(set(ckb_names)), 'Duplicated CKB identifiers'
    metadata = {row['phenocode']: row for row in json.loads((ROOT / 'protein_list_sources/ckb_phenotypes_20260917.json').read_text())}
    assert set(ckb_names) <= metadata.keys()
    excluded = {name for name in ckb_names if 'num_cases' in metadata[name] or name in OTHER_TRAITS}
    ukb = {ALIASES.get(name, (name, ''))[0] for name in ukb_names}
    ckb = set(ckb_names) - excluded
    # All protein targets in this snapshot must be matched or explicitly reviewed.
    assert not ckb - ukb, sorted(ckb - ukb)
    assert all(metadata[name].get('num_samples') in (3968, 3974) for name in ckb)
    notes = {canonical: note for canonical, note in ALIASES.values()}
    protein_rows = []
    for name in sorted(ukb | ckb):
        note = notes.get(name, '')
        if name not in ckb:
            note = 'CKB URL 清单未列出此靶标；不据此推断研究中未测量'
        protein_rows.append({'protein': name, 'UKB': 'Y' if name in ukb else 'N', 'CKB': 'Y' if name in ckb else 'N', '备注': note})
    fields = ['protein', 'UKB', 'CKB', '备注']
    write_table('UKB_CKB_protein_list.tsv', protein_rows, fields)
    write_table('UKB_CKB_protein_list.csv', protein_rows, fields, ',')
    excluded_rows = [{'protein': name, 'UKB': 'N', 'CKB': 'Y', '备注': '非蛋白：' + ('疾病表型（官方元数据含病例/对照数）' if 'num_cases' in metadata[name] else '血压或人体测量性状')} for name in sorted(excluded)]
    write_table('UKB_CKB_all_entries.tsv', sorted(protein_rows + excluded_rows, key=lambda row: row['protein']), fields)
    write_table('CKB_nonprotein_entries.tsv', excluded_rows, fields)
    assert sum(row['UKB'] == 'Y' for row in protein_rows) == len(ukb_names)
    assert sum(row['CKB'] == 'Y' for row in protein_rows) == len(ckb)
    assert len(protein_rows) == len({row['protein'] for row in protein_rows})
    summary = {'UKB_source_rows': len(ukb_rows), 'UKB_unique_targets': len(ukb), 'CKB_URLs': len(ckb_names), 'CKB_protein_targets': len(ckb), 'CKB_nonprotein_entries': len(excluded), 'both': len(ukb & ckb), 'UKB_only': sorted(ukb - ckb), 'CKB_only': sorted(ckb - ukb), 'protein_union': len(ukb | ckb), 'all_entries_union': len(protein_rows + excluded_rows), 'renamed_target_groups': len(ALIASES)}
    (ROOT / 'UKB_CKB_protein_list.summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
