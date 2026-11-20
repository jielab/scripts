# Food-liking traits

Place flt beside diet and 0f in scripts. The original modeling and figure programs are in f. The two entry points share diet's local/LSF interpreter selection and path table. Use Linux or WSL on the local machine.

```bash
./flt.sh run --submit local
./flt.paper.sh run --submit local
```

Use `--submit hpc --wait TRUE` on both entries for sequential LSF execution. `--plan-only` writes a plan without running the analysis, and `--dry-run` prints stage commands. Each stage can be rerun independently. Use the same `--dir-output` for analysis and publication. The default is `$ANALYSIS_ROOT/flt`.

| Input | Default location |
|---|---|
| Field lists and source RDS | `$DIR_UKB/common/*.lst`, `$DIR_UKB/Rdata` |
| RAP data | `$DIR_UKB/rap`. Existing extract files take precedence. Missing extracts are generated from `rap/raw/pheno.tab.gz` in the FLT output directory. |
| Food-liking dictionary | `$DIR_FILES/diet/flt.xlsx`, sheets Sheet1, Sheet2, Sheet4 and item182. item182 supplies the phenotype names and GWAS IDs. |
| GWAS | `$DIR_GWAS/{flt,brain,covs}/common/<id>/gwas/<id>.gz`. Covariates are income, education and physical_activity. |
| Intake / proteins / imaging | `Rdata/diet.fudan.rds` from UKB DISCO scoring, `Rdata/prot.rds`, `rap/img.tab.gz` and `$DIR_FILES/diet/MRI_analysis_dictionary.csv`. |
| Genomic references | Set `FLT_MR_LD_PREFIX` to the PLINK 1KG EUR bfile prefix if different from `$FLT_LDREF/opengwas/1kg.v3/EUR`. `FLT_LDSC_REF` and `FLT_HDL_REF` can point directly to the LD-score and HDL reference directories. |
| SMR | `FLT_SMR` selects the executable. `FLT_EQTLGEN` and `FLT_GTEX` select the two cis-BESD prefixes under `$DIR_GWAS/eqtl`. |
| MRI graphics | `FLT_MRI_ATLAS` points to the shared atlas directory containing ukb_tracts27 and subcortical assets. `FLT_MRI_TOOLS` can locate existing ENIGMA/Python extensions. |
| FUMA, Metascape and PheWAS | See below. |

**Analysis sequence.** Prepare phenotypes, then fit Cox, intake, behavioral, protein and MRI models. The five prospective outcomes are all-cause dementia, stroke, sleep disorders, depression and anxiety disorders. The inherited observational model starts at UKB attendance and uses the shared administrative cutoff. UVMR runs in both directions. The reported primary FDR is BH over actual valid IVW estimates within each disease and direction, across the common 182-trait panel. Wald ratios and sensitivity estimators have no primary FDR. The workbook retains the planned-family correction separately. Forward IVW hits enter MVMR with income, education and physical activity. The manuscript's three named follow-up traits are retained for intake, protein and MRI. LDSC/HDL, PLACO and the paired tissue-specific cis-SMR analyses concern chicken liking (G732) and depression (mdd_2025).

**External annotation.** A first `run` computes PLACO and SMR before the PheWAS stage. Use its `03placo/G732__mdd_2025_PLACO_FUMA.txt.gz` with FUMA, or reuse the matching archived return. Set `FLT_FUMA_DIR` to the returned directory containing leadSNPs.txt, genes.txt, MAGMA and tissue outputs. The default is `03placo/FUMA_job787285` under the FLT output root. Put the corresponding GWAS Atlas downloads under `03phewas/GWAS_Atlas/<rsID>.csv`. Run `./flt.sh phewas` to resume. OpenGWAS and GWAS Catalog queries use the inherited implementation and caches. Supply OpenGWAS credentials through your own environment or personal config. For an already populated cache set `FLT_PHEWAS_ONLINE=false`.

The enrichment stage in `./flt.paper.sh fig5` writes the Metascape gene lists and upload instructions under `03fuma/metascape`. Supply the matching `Enrichment_GO/_FINAL_GO.csv` via `FLT_METASCAPE_RESULTS` or the default result directory, then rerun fig5. Existing returns can be supplied before the initial run. Reusing matching, archived releases is necessary for comparison with the manuscript. Stage completion depends on the available service results.

**Runtime.** The existing diet environment supplies standard R and graphics packages. FLT also uses limma, TwoSampleMR, MVMR, MendelianRandomization, MRPRESSO, GenomicSEM, HDL, ieugwasr, robustbase, deTS, wordcloud2 and htmlwidgets. The native UVMR stage can read an offline BWMR implementation via `FLT_BWMR`. Use `FLT_PLINK`, `FLT_RSCRIPT` and `FLT_PYTHON` when the binaries are outside the configured paths. FLT and diet must see the same phenotype and reference releases.

**Output.** `main_plot/` contains the five composite manuscript PDFs and `FLT.source-tables.xlsx`. Primary Cox, both MR directions, MVMR, intake, behavioral, protein, MRI, SMR and gene/pathway evidence tables contain aggregate results. Manuscript checks include the 182-by-5 grid, the 15 concordant protein pairs, the COL18A1 depression result and the 23-gene panel. Each completed MR run records its run directory. Set `FLT_UVMR_RUN` or `FLT_MVMR_RUN` to reuse a particular prior full run.
