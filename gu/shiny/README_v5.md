# GU Shiny haplotype viewer v5

This patch changes only the Shiny viewer (`gu/shiny/server.R` and `gu/shiny/ui.R`).
The non-Shiny analysis pipeline, `gu.sh`, `loci.R`, `loci.sh`, normalization, and environment YAML are unchanged.

## What v5 fixes

1. **Real 1000 Genomes sample IDs**
   - Internal matrix columns/copies such as `s1418_2` are used only to extract the phased allele.
   - Display labels are mapped back to the real 1KG sample using, in order:
     1. `report/haplotype_sample_map.tsv` when present;
     2. `mat/<locus>/kg.samples.tsv` as the fallback.
   - If neither mapping is available, the viewer explicitly shows `UNMAPPED(s1418)` rather than pretending that `s1418` is a biological sample ID.

2. **Horizontal haplotypes, not a coloured heatmap**
   - The first five rows are always: Altai, Chagyr, Vindija, Denisova, Denisova25.
   - Below them are N source-matching modern phased haplotypes.
   - Every row runs left-to-right across polymorphic SNPs.
   - Cells print only `A`, `C`, `G`, or `T`; missing/ambiguous calls are blank.
   - The panel has horizontal scrolling instead of compressing hundreds of SNPs into vertical colour stripes.

3. **Mouse-over match information**
   - Hover a modern row label or any base in that row to see:
     - real 1KG sample and haplotype;
     - population/super-population when available;
     - haplotype class and frequency (`n`);
     - `best_arch` and `best_match`;
     - per-archaic `*_match` values already calculated by `loci_avcf`;
     - `carry_risk`;
     - exact SNP position and REF/ALT for the hovered base.

4. **Only polymorphic SNPs are displayed**
   - The viewer keeps biallelic SNPs for which both alleles occur in the modern 1KG matrix inside the inherited interval.
   - `Max polymorphic SNPs` defaults to 150. If there are more sites, diagnostic sites are retained first and the remainder are spread across the interval.

5. **No random control rows in the default panel**
   - v5 follows the requested layout: five archaic rows followed by N matched modern haplotypes.

## Why `./gu.sh normalize` does NOT need to change for v5

`loci_avcf` already precomputes the statistics needed by the viewer in `report/hap_match.tsv`:
`*_match`, `best_arch`, `best_match`, `carry_risk`, `hap_id`, `n`, and `copies`.
It also writes/uses the matrix sample ordering needed to map `s1`, `s2`, ... back to real 1KG samples.
Therefore v5 reads the existing locus-level outputs directly instead of duplicating the same values in `gu.sqlite`.

This is preferable for the current locus browser because only the selected locus matrix is read on demand. If GU later needs to browse thousands of loci simultaneously, these metadata can be added to the normalization database as a performance optimization.

## Install

Replace only:

```text
gu/shiny/server.R
gu/shiny/ui.R
```

Then restart:

```bash
./gu.sh shiny
```

No rerun of `loci_avcf` and no `./gu.sh normalize` are required **provided the existing locus output contains `kg.samples.tsv` (or `haplotype_sample_map.tsv`) and `hap_match.tsv`**.

If the viewer shows `UNMAPPED(s####)`, inspect the selected locus output, for example:

```bash
ls -lh /mnt/d/analysis/gu/loci_avcf/<trait>/mat/<locus_id>/kg.samples.tsv
ls -lh /mnt/d/analysis/gu/loci_avcf/<trait>/report/haplotype_sample_map.tsv
ls -lh /mnt/d/analysis/gu/loci_avcf/<trait>/report/hap_match.tsv
```

The first or second file is needed for real sample labels; `hap_match.tsv` is needed for the precomputed match statistics.
