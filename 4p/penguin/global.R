options(shiny.maxRequestSize = 3000*1024^2)
options(timeout = 1000)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
options(warn = -1)

# ---- Single source of truth for data path ----
PATHS <- list(
  GWAS_DIR = file.path("data", "gwas")  # all GWAS live here
)

# ---- Packages ----
pacman::p_load(
  shiny, shinycustomloader, devtools, data.table, CMplot, gson, BiocManager, aplot,
  LDlinkR, DT, golem, ggrepel, reshape2, shinyjs, MendelianRandomization, corrplot,
  shinyWidgets, dplyr, ggplot2, ggfun, igraph, ggraph, MendelianRandomization,
  metap, gprofiler2, tidyverse
)

# GitHub packages
if (!require("GenomicSEM", character.only=TRUE))      install_github("GenomicSEM/GenomicSEM")
if (!require("ieugwasr", character.only=TRUE))        devtools::install_github("mrcieu/ieugwasr")
if (!require("TwoSampleMR", character.only=TRUE))     devtools::install_github("MRCIEU/TwoSampleMR")
if (!require("LDheatmap", character.only=TRUE))       devtools::install_github("SFUStatgen/LDheatmap")
if (!require("eQTpLot", character.only=TRUE))         devtools::install_github("RitchieLab/eQTpLot")
if (!require("hyprcoloc", character.only=TRUE))       devtools::install_github(
  "jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE
)
if (!require("IntAssoPlot", character.only=TRUE))     devtools::install_github("whweve/IntAssoPlot")
if (!require("locuscomparer", character.only=TRUE))   devtools::install_github("boxiangliu/locuscomparer")
pacman::p_load(GenomicSEM, ieugwasr, TwoSampleMR, eQTpLot, IntAssoPlot, hyprcoloc, locuscomparer)

# Bioconductor packages
if (!require("org.Hs.eg.db", character.only=TRUE))    BiocManager::install("org.Hs.eg.db")
if (!require("SNPRelate", character.only=TRUE))       BiocManager::install("SNPRelate")
if (!require("gdsfmt", character.only=TRUE))          BiocManager::install("gdsfmt")
if (!require("HDO.db", character.only=TRUE))          BiocManager::install("HDO.db")
if (!require("DOSE", character.only=TRUE))            BiocManager::install("DOSE")
if (!require("enrichplot", character.only=TRUE))      BiocManager::install("enrichplot")
if (!require("clusterProfiler", character.only=TRUE)) BiocManager::install("clusterProfiler")
if (!require("msigdb", character.only=TRUE))          BiocManager::install("msigdb")
if (!require("multiMiR", character.only=TRUE))        BiocManager::install("multiMiR")
pacman::p_load(org.Hs.eg.db, SNPRelate, gdsfmt, HDO.db, DOSE, enrichplot, clusterProfiler, msigdb, multiMiR)

quiet <- function(x) { sink(tempfile()); on.exit(sink()); invisible(force(x)) }
