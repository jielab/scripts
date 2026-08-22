#!/usr/bin/env bash
# CHANGE LOG v45 (2026-07-31):
# - Write modern VCFs under vcf/ and sample lists under id/ (the old script later
#   referenced these directories but wrote files in the working directory).
# - Create both male and female lists for correct chrX PAR/non-PAR processing.
# - Use shared b37/b38 PAR coordinates from 0phe.f.sh.
# - Re-bgzip Chagyr files in the Chagyr directory, not the unrelated current directory.
set -euo pipefail

dir0=${dir0:-/mnt/d}; dir_ref=${dir_ref:-$dir0/data.BIG/refGen}; dirfunc=${dirfunc:-$dir0/scripts/0f}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
phe_sh=${phe_sh:-$dirfunc/0phe.f.sh}; [[ -s $phe_sh ]] || phe_sh=$script_dir/f/0phe.f.sh
[[ -s $phe_sh ]] || { echo "ERROR: missing 0phe.f.sh" >&2; exit 1; }; source "$phe_sh"
for x in wget bcftools tabix bgzip plink2 awk sort; do command -v "$x" >/dev/null || { echo "ERROR: missing command $x" >&2; exit 1; }; done

GRCH=${GRCH:-37}; split=b${GRCH}; [[ $split == b37 || $split == b38 ]] || { echo "ERROR GRCH must be 37 or 38" >&2; exit 1; }
modern_root=${MODERN_ROOT:-$PWD}; archaic_root=${ARCHAIC_ROOT:-$dir_ref/archaic/avcf}
vcf_dir=$modern_root/vcf; id_dir=$modern_root/id; pfile_dir=$modern_root/pfile
mkdir -p "$vcf_dir" "$id_dir" "$pfile_dir" "$archaic_root"/{Vindija,Altai,Chagyr,Denisova,Denisova25}

: "${sample_file:?Set sample_file to samples_v3.ALL.panel before running gen.clean.sh}"
[[ -s $sample_file ]] || { echo "ERROR sample_file missing: $sample_file" >&2; exit 1; }

if [[ $GRCH == 37 ]]; then
    baseurl=https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3
else
    echo "ERROR: automatic UCSC download is configured only for the hg19/GRCh37 1KG Phase 3 source. Supply GRCh38 VCFs separately." >&2; exit 1
fi
index=$modern_root/.1000g.index.html; files=$modern_root/.1000g.files.txt
wget -q -O "$index" "$baseurl/"
grep -oP 'href="\K[^"]+' "$index" | grep -vE '^\?|^/|\.\./|/$' > "$files"
while IFS= read -r f; do wget -c --inet4-only --timeout=20 --tries=5 --waitretry=5 -P "$vcf_dir" "$baseurl/$f"; done < "$files"
for chr in {1..22} X Y MT; do
    f=$(find "$vcf_dir" -maxdepth 1 -name "ALL.chr${chr}.*.vcf.gz" -printf '%f\n' | sort | sed -n '1p')
    [[ -n $f ]] || { echo "ERROR: missing chr$chr VCF under $vcf_dir" >&2; exit 1; }
    [[ $f == chr${chr}.vcf.gz ]] || mv -f "$vcf_dir/$f" "$vcf_dir/chr${chr}.vcf.gz"
    [[ -f $vcf_dir/$f.tbi ]] && mv -f "$vcf_dir/$f.tbi" "$vcf_dir/chr${chr}.vcf.gz.tbi"
    tabix -f -p vcf "$vcf_dir/chr${chr}.vcf.gz"
done

# 1KG panel columns: sample, pop, super_pop, gender.
for race in ALL EUR AFR EAS SAS AMR; do
    awk -v race="$race" 'NR>1 && (race=="ALL" || $3==race){print $1}' "$sample_file" > "$id_dir/$race.1id"
    awk -v race="$race" 'NR>1 && (race=="ALL" || $3==race){print $1,$1}' "$sample_file" > "$id_dir/$race.2id"
    awk -v race="$race" 'NR>1 && (race=="ALL" || $3==race){s=(tolower($4)=="male"||tolower($4)=="m")?1:(tolower($4)=="female"||tolower($4)=="f")?2:0; print $1,$1,s}' "$sample_file" > "$id_dir/$race.2id_sex"
    for sex in male female; do
        code=1; [[ $sex == female ]] && code=2
        awk -v race="$race" -v sex="$sex" 'NR>1 && (race=="ALL" || $3==race) && tolower($4)==sex{print $1}' "$sample_file" > "$id_dir/$race.$sex.1id"
        awk -v race="$race" -v sex="$sex" 'NR>1 && (race=="ALL" || $3==race) && tolower($4)==sex{print $1,$1}' "$sample_file" > "$id_dir/$race.$sex.2id"
        awk -v race="$race" -v sex="$sex" -v code="$code" 'NR>1 && (race=="ALL" || $3==race) && tolower($4)==sex{print $1,$1,code}' "$sample_file" > "$id_dir/$race.$sex.2id_sex"
    done
done
bcftools query -l "$vcf_dir/chrX.vcf.gz" | sort -u > "$modern_root/.chrX.samples.txt"
sort -u "$id_dir/ALL.1id" > "$modern_root/.panel.samples.txt"
[[ $(sort "$id_dir/ALL.1id" | uniq -d | wc -l) -eq 0 ]] || { echo "ERROR: duplicate sample IDs in sample_file" >&2; exit 1; }
cmp -s "$modern_root/.panel.samples.txt" "$modern_root/.chrX.samples.txt" || { echo "ERROR: chrX sample set differs from sample_file" >&2; exit 1; }

for chr in {1..22} X Y; do
    extra=()
    if [[ $chr == X ]]; then extra=(--update-sex "$id_dir/ALL.2id_sex" --split-par "$split" --set-missing-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing)
    elif [[ $chr == Y ]]; then extra=(--update-sex "$id_dir/ALL.male.2id_sex" --keep "$id_dir/ALL.male.2id"); fi
    plink2 --vcf "$vcf_dir/chr$chr.vcf.gz" --double-id --allow-extra-chr "${extra[@]}" --make-pgen pvar-cols=xheader --out "$pfile_dir/chr$chr"
    plink2 --pfile "$pfile_dir/chr$chr" --keep "$id_dir/EUR.2id" --make-pgen --out "$pfile_dir/EUR.chr$chr"
done

# Explicit chrX diagnostic subsets used to verify ploidy handling.
Xvcf=$vcf_dir/chrX.vcf.gz; par_var=X_PAR_${split}; nonpar_var=X_NONPAR_${split}
for sex in male female; do
    bcftools view -S "$id_dir/EUR.$sex.1id" -r "${!par_var}" -m2 -M2 -v snps -Oz -o "$vcf_dir/EUR.$sex.chrX.par.vcf.gz" "$Xvcf"; tabix -f -p vcf "$vcf_dir/EUR.$sex.chrX.par.vcf.gz"
    bcftools view -S "$id_dir/EUR.$sex.1id" -r "${!nonpar_var}" -m2 -M2 -v snps -Oz -o "$vcf_dir/EUR.$sex.chrX.nonPAR.vcf.gz" "$Xvcf"; tabix -f -p vcf "$vcf_dir/EUR.$sex.chrX.nonPAR.vcf.gz"
done

if [[ ${DOWNLOAD_ARCHAIC:-1} == 1 ]]; then
    for c in {1..22} X; do
        wget -c -P "$archaic_root/Vindija" "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
        wget -c -P "$archaic_root/Altai" "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
        wget -c -P "$archaic_root/Chagyr" "https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz"{,.tbi}
        wget -c -P "$archaic_root/Denisova" "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
        wget -c -P "$archaic_root/Denisova25" "https://cdna.eva.mpg.de/denisova/Den25/VCF/chr${c}.Den25.L35MQ25.B30.map35_100.vcf.gz"{,.tbi}
    done
    find "$archaic_root/Chagyr" -maxdepth 1 -name '*.vcf.gz' -print0 | while IFS= read -r -d '' f; do
        if ! gzip -t "$f" 2>/dev/null || ! tabix -l "$f" >/dev/null 2>&1; then mv "$f" "$f.old.gz"; gunzip -c "$f.old.gz" | bgzip -@ 4 -c > "$f"; fi
        tabix -f -p vcf "$f"
    done
    find "$archaic_root" -mindepth 2 -maxdepth 2 -name '*.vcf.gz' -print0 | xargs -0 -P 8 -I {} tabix -f -p vcf '{}'
fi

echo "ALL DONE: modern=$modern_root archaic=$archaic_root genome_build=$split"
