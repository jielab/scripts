#!/usr/bin/env bash
# Exercise both public backends on a small, real mixed-sex 1KG X subset.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
SOURCE=${1:-/mnt/d/data.BIG/refGen/1kg/37}
WORK=$(mktemp -d /tmp/gu-arg-x-smoke.XXXXXX)
echo "SMOKE_WORK=$WORK"
export PATH="$HOME/anaconda3/envs/gu/bin:$PATH"
mkdir -p "$WORK/pfile" "$WORK/vcf"
python3 - "$SOURCE/pfile/chrX.psam" "$WORK" <<'PY'
import sys
from pathlib import Path
lines=Path(sys.argv[1]).read_text().splitlines(); header=lines[0].lstrip('#').split()
rows=[dict(zip(header, r.split())) for r in lines[1:]]
names=[r['IID'] for r in rows if r['SEX']=='1'][:20]+[r['IID'] for r in rows if r['SEX']=='2'][:142]
p=Path(sys.argv[2]); (p/'keep.txt').write_text('#FID\tIID\n'+''.join(f'{r.get("FID", "0")}\t{r["IID"]}\n' for r in rows if r['IID'] in names))
(p/'vcf.samples').write_text('\n'.join(names)+'\n')
(p/'ancestry.tsv').write_text('IID\tancestry\n'+''.join(f'{n}\tEUR\n' for n in names))
PY
"$HOME/miniforge3/envs/grid/bin/plink2" --pfile "$SOURCE/pfile/chrX" vzs \
  --keep "$WORK/keep.txt" --chr X --from-bp 3000000 --to-bp 3100000 \
  --snps-only just-acgt --max-alleles 2 --maf 0.05 --geno 0 --thin-count 128 --seed 42 \
  --make-pgen vzs --threads 2 --out "$WORK/pfile/chrX" > "$WORK/plink.log" 2>&1
bcftools view -S "$WORK/vcf.samples" -r X:3000000-3020000 "$SOURCE/vcf/chrX.vcf.gz" -Oz -o "$WORK/vcf/chrX.vcf.gz"
tabix -p vcf "$WORK/vcf/chrX.vcf.gz"
for method in needle tsinfer; do
  echo "BUILD $method mixed-sex non-PAR X"
  bash "$ROOT/arg.sh" build --dir-gen "$WORK" --method "$method" --chr X --format trace \
    --features FALSE --threads 2 --jobs 1 --arg-dir "$WORK/arg.$method" \
    --ancestry-file "$WORK/ancestry.tsv" --scratch "$WORK/scratch" \
    --needle-mode sequence --seed-haplotypes 32 --normalize FALSE \
    > "$WORK/$method.log" 2>&1 || { tail -60 "$WORK/$method.log"; exit 1; }
  tail -5 "$WORK/$method.log"
done
python3 - "$WORK" <<'PY'
import csv,gzip,sys,tskit,pysam
from pathlib import Path
p=Path(sys.argv[1])
for method in ('needle','tsinfer'):
    output=p/f'arg.{method}'/'trace'/method
    ts=tskit.load(output/'chrX.trees')
    assert ts.num_samples==304, (method,ts.num_samples)
    rows=list(csv.DictReader((output/'chrX.sample_map.tsv').open(),delimiter='\t'))
    counts={}
    for row in rows: counts[row['sample']]=counts.get(row['sample'],0)+1
    assert sum(v==1 for v in counts.values())==20, (method,counts)
    assert sum(v==2 for v in counts.values())==142
    if method == 'needle':
        expected={}
        with gzip.open(p/'pfile.4arg/chrX.haps.gz','rt') as src:
            for line in src:
                z=line.split(); expected[int(z[2])]=[z[3+int(g)] for g in z[5:]]
    else:
        expected={}
        with pysam.VariantFile(str(p/'vcf.4arg/chrX.vcf.gz')) as src:
            for record in src:
                expected[record.pos]=[record.alleles[a] for call in record.samples.values() for a in call['GT']]
    for variant in ts.variants():
        observed=[variant.alleles[g] for g in variant.genotypes]
        assert observed==expected[int(variant.site.position)], (method,variant.site.position)
    print(f'{method}: 20 males + 142 females = 304 real haplotypes; TRACE PASS')
PY
echo "SMOKE PASS: $WORK"
