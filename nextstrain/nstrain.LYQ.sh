#!/usr/bin/env bash


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Shell options, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GISAID: jiehuang Abc12345
# Web version: https://www.ebi.ac.uk/jdispatcher/msa
# Download hchikv.ref.gb (https://www.ncbi.nlm.nih.gov/nuccore/NC_004162.2).
# Download https://developers.google.com/public-data/docs/canonical/countries_csv.
# You can drag color.vip.tsv into the Nextstrain web UI.
dir0=/mnt/d
dirdat=$dir0/data/viome/raw
dirout=$dir0/analysis/viome; mkdir -p $dirout

dos2unix *
awk -F '\t' -v OFS='\t' 'NR>1 {print "country", $4, $2, $3}' countries.tsv | sed "s/ /_/g; s/'//g" > countries.centroids.tsv

sed -E 's/\|[^|]*$//' gisaid.fasta > sequences.fasta
paste <(cut -f 1,10 seq_tech.tsv) <(cut -f 1-5 patient_status.tsv) | awk -F '\t' -v OFS='\t' '{gsub(" ", "_", $1); gsub(" ", "_", $3); print}' | \
	sed '1s/Virus_name/strain/; 1s/Sequencing technology/seq_tech/; 1s/Accession ID/accession/; 1s/Collection date/date/' | \
	awk -F'\t' -v OFS='\t' 'NR==1{$3="virus"; $NF=""; print $0, "region", "country", "province", "city"; next} {$1=$1"|"$5; $3="hChikV"; split($NF, loc, " / ")
	region = (loc[1] != "") ? loc[1] : "unknown"; country = (loc[2] != "") ? loc[2] : "unknown"; gsub(" ", "_", region); gsub(" ", "_", country)
	province = (loc[3] != "") ? loc[3] : "unknown"; city = (loc[4] != "") ? loc[4] : "unknown"; gsub(" ", "_", province); gsub(" ", "_", city)
	$NF=""; print $0, region, country, province, city }' | sed "s/\t\t/\t/; s/'//g" | \
	awk -F'\t' -v OFS='\t' 'NR==1{for(i=1; i<=NF; i++) if($i == "date") dcol=i; print; next} {d=$dcol; if(d=="unknown"){} else if (d ~ /^[0-9]{4}$/){$dcol=d"-07-01"} else if(d ~ /^[0-9]{4}-[0-9]{2}$/){$dcol=d"-15"}; print}' 
> metadata.tsv

# Date or region unknown.
awk -v OFS="\t" '{if($10=="SZ") print $1, "red"; else if ($8=="China") print $1, "blue"}' metadata.tsv > vip.color
cut -f 8 metadata.tsv | sort | uniq | awk '{ printf "country\t%s\t#%06X\n", $0, rand()*16777215 }' | sed "s/ /_/g; s/'//g" > colors.tsv


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Preparation.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
conda activate nextstrain
mkdir -p auspice

# Sequence deduplication, already done.
python3 << 'EOF'
import pandas as pd
from Bio import SeqIO
meta = pd.read_csv('metadata.tsv', sep='\t')
meta_strains = set(meta['strain'].astype(str).str.strip())
seen = set()
unique = []
for rec in SeqIO.parse('sequences.fasta', 'fasta'):
    name = rec.id.strip()
    if name in meta_strains and name not in seen:
        seen.add(name)
        unique.append(rec)
SeqIO.write(unique, 'sequences.fasta', 'fasta')
print(f"保留序列数: {len(unique)}")
EOF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 2. Build sequence index.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur index --sequences sequences.fasta --output sequence_index.tsv


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 3. Filter low-quality sequences.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note: avoid --exclude-where N_ratio>0.05 because that option does not support numeric comparison.
# Use only --min-length for length filtering.
augur filter \
  --sequences sequences.fasta \
  --metadata metadata.tsv \
  --min-length 11000 \
  --output-sequences filtered.fasta \
  --output-metadata filtered_metadata.tsv


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 4. Sequence alignment.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur align \
  --sequences filtered.fasta \
  --reference-sequence chikv_ref.fasta \
  --output alignment.fasta \
  --fill-gaps


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 5. Build phylogenetic tree.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur tree \
  --alignment alignment.fasta \
  --tree-builder-args '-ninit 10 -n 4' \
  --output tree_raw.nwk


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 6. Time calibration and tree optimization.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur refine \
  --tree tree_raw.nwk \
  --alignment alignment.fasta \
  --metadata filtered_metadata.tsv \
  --timetree \
  --output-tree tree.nwk \
  --output-node-data branch_lengths.json


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 7. Add geographic coordinates.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If filtered_metadata.tsv already has latitude and longitude columns, skip this step.
python3 << 'EOF'
import pandas as pd
meta = pd.read_csv('filtered_metadata.tsv', sep='\t')
coords = {
    'China': (35.86, 104.19),
    'India': (20.59, 78.96),
    'Thailand': (15.87, 100.99),
    'Brazil': (-14.24, -51.93),
    'Madagascar': (-18.77, 46.87),
    'France': (46.23, 2.21),
    'USA': (37.09, -95.71),
    'Italy': (41.87, 12.57),
    'Singapore': (1.35, 103.82),
    'Malaysia': (4.21, 101.98),
    'Indonesia': (-0.79, 113.92),
    'Philippines': (12.88, 121.77),
    'Sri Lanka': (7.87, 80.77),
    'Maldives': (3.20, 73.22),
    'Myanmar': (21.91, 95.96),
    'Vietnam': (14.06, 108.28),
    'Cambodia': (12.57, 104.99),
    'Laos': (19.86, 102.50),
    'Bangladesh': (23.68, 90.36),
    'Pakistan': (30.38, 69.35),
    'Nepal': (28.39, 84.12),
    'Bhutan': (27.51, 90.43),
    'Japan': (36.20, 138.25),
    'South Korea': (35.91, 127.77),
    'Australia': (-25.27, 133.78),
    'New Zealand': (-40.90, 174.89),
    'Kenya': (-0.02, 37.91),
    'Nigeria': (9.08, 8.68),
    'South Africa': (-30.56, 22.94),
    'Uganda': (1.37, 32.29),
    'Angola': (-11.20, 17.87),
    'Congo': (-0.23, 15.83),
    'Gabon': (-0.80, 11.61),
    'Senegal': (14.50, -14.45),
    'Mauritius': (-20.35, 57.55),
    'Seychelles': (-4.68, 55.49),
    'Comoros': (-11.88, 43.87),
    'Reunion': (-21.115141, 55.536384),
    'Mayotte': (-12.83, 45.17),
    'Turks and Caicos Islands': (21.69, -71.79),
    'Trinidad and Tobago': (10.69, -61.22),
    'Jamaica': (18.11, -77.30),
    'Haiti': (18.97, -72.29),
    'Dominican Republic': (18.74, -70.16),
    'Puerto Rico': (18.22, -66.59),
    'Venezuela': (6.42, -66.59),
    'Colombia': (4.57, -74.30),
    'Ecuador': (-1.83, -78.18),
    'Peru': (-9.19, -75.02),
    'Bolivia': (-16.29, -63.59),
    'Paraguay': (-23.44, -58.44),
    'Argentina': (-38.42, -63.62),
    'Chile': (-35.68, -71.54),
    'Uruguay': (-32.52, -55.77),
    'Guyana': (4.86, -58.93),
    'Suriname': (3.92, -56.03),
    'French Guiana': (3.93, -53.13),
    'Panama': (8.54, -80.78),
    'Costa Rica': (9.75, -83.75),
    'Nicaragua': (12.87, -85.21),
    'Honduras': (15.20, -86.24),
    'El Salvador': (13.79, -88.90),
    'Guatemala': (15.78, -90.23),
    'Belize': (17.19, -88.50),
    'Mexico': (23.63, -102.55),
    'Canada': (56.13, -106.35),
    'Russia': (61.52, 105.32),
    'Germany': (51.17, 10.45),
    'United Kingdom': (55.38, -3.44),
    'Spain': (40.46, -3.75),
    'Portugal': (39.40, -8.22),
    'Netherlands': (52.13, 5.29),
    'Belgium': (50.50, 4.47),
    'Switzerland': (46.82, 8.23),
    'Austria': (47.52, 14.55),
    'Sweden': (60.13, 18.64),
    'Norway': (60.47, 8.47),
    'Denmark': (56.26, 9.50),
    'Finland': (61.92, 25.75),
    'Poland': (51.92, 19.15),
    'Czech Republic': (49.82, 15.47),
    'Slovakia': (48.67, 19.70),
    'Hungary': (47.16, 19.50),
    'Romania': (45.94, 24.97),
    'Bulgaria': (42.73, 25.49),
    'Greece': (39.07, 21.82),
    'Turkey': (38.96, 35.24),
    'Israel': (31.05, 34.85),
    'Saudi Arabia': (23.89, 45.08),
    'UAE': (23.42, 53.85),
    'Qatar': (25.35, 51.18),
    'Kuwait': (29.31, 47.48),
    'Oman': (21.51, 55.92),
    'Yemen': (15.55, 48.52),
    'Iran': (32.43, 53.69),
    'Iraq': (33.22, 43.68),
    'Afghanistan': (33.94, 67.71),
    'Kazakhstan': (48.02, 66.92),
    'Uzbekistan': (41.38, 64.59),
    'Turkmenistan': (38.97, 59.56),
    'Kyrgyzstan': (41.20, 74.77),
    'Tajikistan': (38.86, 71.28),
    'Mongolia': (46.86, 103.82),
    'North Korea': (40.34, 127.51),
    'Taiwan': (23.70, 121.00),
    'Hong Kong': (22.32, 114.17),
    'Macau': (22.20, 113.55)
}

meta['latitude'] = meta['country'].map(lambda x: coords.get(x, (0, 0))[0])
meta['longitude'] = meta['country'].map(lambda x: coords.get(x, (0, 0))[1])
meta.to_csv('filtered_metadata_with_coords.tsv', sep='\t', index=False)
EOF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 8. Trait inference, transmission route.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run only country and division; skip location because too many states can cause errors.
augur traits \
  --tree tree.nwk \
  --metadata filtered_metadata_with_coords.tsv \
  --columns country \
  --output-node-data traits_country.json

augur traits \
  --tree tree.nwk \
  --metadata filtered_metadata_with_coords.tsv \
  --columns division \
  --output-node-data traits_division.json


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 9. Reconstruct ancestral sequences and translate amino-acid mutations.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur ancestral \
  --tree tree.nwk \
  --alignment alignment.fasta \
  --output-node-data nt_muts.json \
  --inference joint \
  --root-sequence chikv_ref.fasta

augur translate \
  --tree tree.nwk \
  --ancestral-sequences nt_muts.json \
  --reference-sequence chikv_ref.fasta \
  --output-node-data aa_muts.json


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 10. Prepare latitude/longitude file, strictly tab-delimited.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
python3 << 'EOF'
import pandas as pd, csv
meta = pd.read_csv('filtered_metadata_with_coords.tsv', sep='\t')
lat_longs = meta[['country', 'latitude', 'longitude']].drop_duplicates()
with open('lat_longs.tsv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['country', 'latitude', 'longitude'])
    for _, row in lat_longs.iterrows():
        writer.writerow([row['country'], row['latitude'], row['longitude']])
print("lat_longs.tsv generated")
EOF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 11. Prepare visualization config file.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat > auspice_config.json << 'EOF'
{
  "title": "CHIKV Analysis",
  "maintainers": [{"name": "CDC"}],
  "colorings": [
    {"key": "country", "title": "Country", "type": "categorical"},
    {"key": "division", "title": "Division", "type": "categorical"},
    {"key": "region", "title": "Region", "type": "categorical"},
    {"key": "Genotype", "title": "Genotype", "type": "categorical"},
    {"key": "date", "title": "Date", "type": "continuous"}
  ],
  "display_defaults": {
    "color_by": "country"
  }
}
EOF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 12. Export visualization files.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur export v2 \
  --tree tree.nwk \
  --metadata filtered_metadata_with_coords.tsv \
  --metadata-id-columns strain \
  --node-data branch_lengths.json traits_country.json traits_division.json nt_muts.json aa_muts.json \
  --geo-resolutions country division \
  --lat-longs lat_longs.tsv \
  --auspice-config auspice_config.json \
  --output auspice/chikv_with_map.json


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 13. Start local viewer.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
auspice view --datasetDir auspice
# Visit http://localhost:4000 in a browser.