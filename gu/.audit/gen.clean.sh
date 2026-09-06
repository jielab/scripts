dir0=/mnt/d
source $dir0/scripts/f/0phe.f.sh

GRCH=37
split=b${GRCH}

mkdir -p vcf id pfile


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 下载1000G数据，plink2格式化
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
baseurl=https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3
baseurl=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502 # ⚠️这个版本没有rsID🕳

wget -q -O index.html "$baseurl/"
grep -oP 'href="\K[^"]+' index.html | grep -v '^?' | grep -v '^/' | grep -v '../' | grep -v '/$' > files.txt
cat files.txt | while read f; do wget -c --inet4-only --timeout=20 --tries=5 --waitretry=5 "$baseurl/$f"; done

for chr in {1..22} X Y MT; do
	f=$(find . -maxdepth 1 -name "ALL.chr${chr}.*.vcf.gz" -printf "%f\n" | sort | head -1)
	[[ -n "$f" ]] || { echo "ERROR: missing chr$chr VCF"; exit 1; }
	echo "Rename: $f -> chr${chr}.vcf.gz"
	mv "$f" "chr${chr}.vcf.gz"
	[[ -f "$f.tbi" ]] && mv "$f.tbi" "chr${chr}.vcf.gz.tbi"
done

# 拆分样本文件
sample_file=$dir0/refGen/1kg/phase3.ebi/vcf/samples_v3.ALL.panel
for race in ALL EUR AFR EAS SAS AMR; do
	awk -v race=$race 'NR>1 && (race == "ALL" || $3 == race) {print $1}' $sample_file > $race.1id
	awk -v race=$race 'NR>1 && (race == "ALL" || $3 == race) {print $1, $1}' $sample_file > $race.2id
	awk -v race=$race 'NR>1 && (race == "ALL" || $3 == race) {s=(tolower($4)== "male"||$4=="m")?1:(tolower($4)=="female"||$4=="f")?2:0; print $1,$1,s}' $sample_file > $race.2id_sex
	awk -v race=$race 'NR>1 && (race == "ALL" || $3 == race) && tolower($4)== "male" {print $1}' $sample_file > $race.male.1id
	awk -v race=$race 'NR>1 && (race == "ALL" || $3 == race) && tolower($4)== "male" {print $1,$1}' $sample_file > $race.male.2id
	awk -v race=$race 'NR>1 && (race == "ALL" || $3 == race) && tolower($4)== "male" {print $1,$1,1}' $sample_file > $race.male.2id_sex
done

# 确认chrX 里面的样本是一样的
vcf=vcf/chrX.vcf.gz
bcftools query -l $vcf > id.tmp; diff id.tmp id/ALL.1id | head
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $vcf | head

# 转换为
for chr in {1..22} X Y; do
	# tabix -f -p vcf "chr${chr}.vcf.gz"
	if [[ "$chr" == "X" ]]; then
		str="--update-sex id/ALL.2id_sex --split-par $split --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 missing" # 🏮
	elif [[ "$chr" == "Y" ]]; then
		str="--update-sex id/ALL.male.2id_sex --keep id/ALL.male.2id"
	else
		str=""
	fi
	plink2 --vcf vcf/chr$chr.vcf.gz --double-id --allow-extra-chr $str --make-pgen --out pfile/chr$chr
	plink2 --pfile pfile/EUR.chr$chr --max-alleles 2 --make-bed --out pfile/EUR.chr$chr # 🏮
done

Xvcf=vcf/chrX.vcf.gz
	par_var=X_PAR_${split}; nonpar_var=X_NONPAR_${split}
	bcftools view -S id/EUR.male.1id -r ${!par_var} -m2 -M2 -v snps -Oz -o vcf/EUR.male.chrX.par.vcf.gz "$Xvcf"; tabix -f -p vcf vcf/EUR.male.chrX.par.vcf.gz
	bcftools view -S id/EUR.male.1id -r ${!nonpar_var} -m2 -M2 -v snps -Oz -o vcf/EUR.male.chrX.nonPar.vcf.gz "$Xvcf"; tabix -f -p vcf vcf/EUR.male.chrX.nonPar.vcf.gz
	plink2 --vcf vcf/EUR.male.chrX.par.vcf.gz --double-id --allow-extra-chr --update-sex id/EUR.male.2id_sex --split-par $split --make-pgen --out pfile/EUR.male.chrX.par
	plink2 --vcf vcf/EUR.male.chrX.nonPar.vcf.gz --double-id --allow-extra-chr --update-sex id/EUR.male.2id_sex --make-pgen --out pfile/EUR.male.chrX.nonPar


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 下载古基因数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir -p Vindija Altai Chagyr Denisova Denisova25 

for c in {1..22} X; do
	wget -c -P Vindija	"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Altai	"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Chagyr	"https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz"{,.tbi}
	wget -c -P Denisova		"https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr${c}_mq25_mapab100.vcf.gz"{,.tbi}
	wget -c -P Denisova25	"https://cdna.eva.mpg.de/denisova/Den25/VCF/chr${c}.Den25.L35MQ25.B30.map35_100.vcf.gz"{,.tbi}
done

# 对 Chagyr, 重新 gzip 
ls -1 *.vcf.gz | parallel -j 4 'echo {}; mv {} {}.old.gz; gunzip -c {}.old.gz | bgzip -@ 4 -c > {} && tabix -f -p vcf {}'

ls -1 */*.vcf.gz | xargs -P 8 -I {} sh -c 'echo "{}"; tabix -f -p vcf "{}"'