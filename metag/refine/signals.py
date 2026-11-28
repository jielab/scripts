"""Describe exact published assignments, retaining competing explanations for source."""
from collections import Counter, defaultdict
from metag_refine_tables import HERE, rows

# The split normalizes labels below species up to their species. Match only published
# species bins. Including ICTV and previous names permits both taxonomy generations.
WATCH = {
    'HBV-related species': {'Orthohepadnavirus hominoidei','Hepatitis B virus'},
    'HIV-1-related species': {'Lentivirus humimdef1','Human immunodeficiency virus 1'},
    'HIV-2-related species': {'Lentivirus humimdef2','Human immunodeficiency virus 2'},
    'influenza A': {'Alphainfluenzavirus influenzae','Influenza A virus'},
    'influenza B': {'Betainfluenzavirus influenzae','Influenza B virus'},
    'influenza C': {'Gammainfluenzavirus influenzae','Influenza C virus'},
    'influenza D': {'Deltainfluenzavirus influenzae','Influenza D virus'},
}
def sources(sample, taxa, evidence):
    foods = {(r['group'],r['rank'],r['name']):r['food_source'] for r in rows(HERE/'food_sources.tsv')}
    wanted = {}
    for r in taxa:
        key = (r['group'], r['rank'], r['name'])
        if r['group'] in ('human','plant','mammal') or key in foods:
            wanted[r['taxid']] = r, foods.get(key, 'unlisted')
    by = defaultdict(Counter)
    for r in evidence:
        tid = r['final_taxid']
        if tid in wanted:
            by[tid][r['evidence']] += int(r['pairs'])
    out = []
    for tid,(r,food) in wanted.items():
        n = int(r['pairs'])
        if sum(by[tid].values()) != n:
            raise ValueError(f'Taxon/evidence counts differ for {sample} taxid {tid}')
        details = '|'.join(f'{key}={v}' for key,v in sorted(by[tid].items()))
        out.append([sample,r['group'],r['rank'],tid,r['name'],n,food,details])
    return out

def virus(sample, taxa, total, library):
    observed = Counter()
    for r in taxa:
        if r['group']=='virus' and r['rank']=='species':
            observed[r['name'].casefold()] += int(r['pairs'])
    out = []
    for target, names in WATCH.items():
        n = sum(observed[name.casefold()] for name in names)
        if n:
            status = 'classifier_signal_requires_confirmation'
        elif library == 'unknown':
            status = 'no_matching_species_bin_library_unknown'
        elif library == 'DNA' and target != 'HBV-related species':
            status = 'no_matching_species_bin_RNA_not_assessed'
        else:
            status = 'no_matching_species_bin'
        out.append([sample,target,n,total,library,status])
    return out
