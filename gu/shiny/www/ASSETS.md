# Genome browser assets

These files keep the Shiny page's JavaScript, hg19 gene annotations, chromosome sizes,
and cytobands available without a request to a third-party server on each page load.
Reference bases use the existing indexed `GRCH37.fasta` / `GRCH38.fasta` under
`$GU_REF_ROOT/fasta` (default `/mnt/e/refGen/fasta`) when available, with
the chromosome 1 length checked against the requested assembly. Shiny serves
byte ranges directly; no multi-GB files are copied. Otherwise the browser falls
back to UCSC. GRCh38 and CHM13 annotations remain remote; an annotation failure
does not prevent reference navigation.

Downloaded on 2026-09-12, without modification:

| Local file | Source | SHA-256 |
| --- | --- | --- |
| `vendor/igv-3.0.0.min.js` | https://cdn.jsdelivr.net/npm/igv@3.0.0/dist/igv.min.js | `9262b079af087eb8c55836640b0b4e80251ab30e26edb9244d1645ad3e487d89` |
| `genomes/hg19/ncbiRefSeq.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz | `25d7be987073d8ecadb042da87837bfc469e432ec044353a5f28c8e9fa893238` |
| `genomes/hg19/cytoBand.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz | `f9b82309b2bca1eb9d91a5cb2c6aa0528351158e6e20b51d82cca36d01735cba` |
| `genomes/hg19/hg19.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes | `b404927655a4aada254ea94ad4da0c8901ed0737e67a0dcabedf673354b1f505` |

IGV.js is distributed under the MIT license; see `vendor/igv-LICENSE`.
Keep the hg19 annotations tied to GRCh37/hg19. When refreshing a cached file,
validate gzip integrity and update its checksum here. No analysis results are
computed from these display assets.

## IBDmix population map

`maps/ne_110m_land.tsv` contains exterior land polygons, rounded to four decimal
places and excluding Antarctica, from Natural Earth v5.1.2:
[ne_110m_land.geojson](https://github.com/nvkelso/natural-earth-vector/blob/v5.1.2/geojson/ne_110m_land.geojson).
Source GeoJSON SHA-256: `9e0729ee253ca7d7a5c4ae9395fb1902264c5377c52e224d13dd85010e2835d9`.
Natural Earth data are [public domain](https://www.naturalearthdata.com/about/terms-of-use/).
The bundled coordinates let R draw the static map without a mapping package or
an external map service.

`maps/1kg_populations.tsv` lists the 26 1000 Genomes population codes, approximate
sampling/ancestral locations, and manually spaced chart positions. GIH, ITU and
STU are placed at their South Asian ancestral locations; these coordinates are
display annotations, not participant residences or inputs to the analysis.
The population percentages and counts always come from the selected dataset.

The visual reference is [Zeberg & Pääbo, Nature 2020, Fig. 3](https://www.nature.com/articles/s41586-020-2818-3/figures/3).
That figure shows a single haplotype's allele frequency; this map instead shows
mean per-individual autosomal Neanderthal physical coverage. Benchmark Mb/person
values come from [Chen et al., Cell 2020, Fig. 2A / Table S4](https://doi.org/10.1016/j.cell.2020.01.012).
Their Altai callset, masks and population-specific calling differ from GU's
multi-reference union. The figures therefore provide scale/trend context, not
a claim of exact reproduction or diploid ancestry dosage.
