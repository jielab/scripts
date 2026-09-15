"""Scientific workflow identity shared by generation and resume checks."""
WORKFLOW = 'gwas_lead_ld_core_archaic5_v2'
LINEAGE_REFS = {'Neanderthal': ('Altai', 'Chagyr', 'Vindija'),
                'Denisovan': ('Denisova', 'Denisova25')}
REFS = tuple(ref for refs in LINEAGE_REFS.values() for ref in refs)
