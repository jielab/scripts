#!/usr/bin/env Rscript
# Deliberate setup only. No analysis script sources this installer.
# Existing requested R packages are retained; there is no blanket update step.
args <- commandArgs(trailingOnly=TRUE)
allowed <- c("--check", "--install", "--with-optional")
if (length(setdiff(args, allowed))) stop("Unknown installer option: ", paste(setdiff(args, allowed), collapse=", "))
if (all(c("--check", "--install") %in% args)) stop("Choose --check or --install.")
install <- "--install" %in% args
optional <- "--with-optional" %in% args
if (getRversion() < "4.5.0") stop("Use R >= 4.5; env.yml selects R 4.5.")

# Direct imports, dynamic package lists, and packages needed by the shared 0f
# functions actually used by the MAHA/DOSE workflows. Base/recommended packages
# are checked but not installed from CRAN here.
recommended <- c("MASS", "Matrix", "lattice", "survival")
cran <- c(
    "tidyverse", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "forcats", "ggplot2",
    "data.table", "R.utils", "readxl", "writexl", "openxlsx", "haven", "jsonlite", "digest",
    "survey", "broom", "rms", "pROC", "irr", "matrixStats", "Rcpp",
    "patchwork", "cowplot", "ggrepel", "scales", "flextable", "gtsummary",
    "ggnewscale", "corrplot", "pheatmap", "svglite",
    "mediation", "lavaan", "enrichR", "httr", "igraph", "sandwich",
    "sf", "raster", "rosm", "ggspatial", "rnaturalearth", "rnaturalearthdata"
)
bioc <- "impute"  # Shared phenotype.R clean_biom()/drop_10pct() when imputation is selected.
if (optional) {
    # Optional SEM diagram and other lazy-loaded 0f plotting/prediction helpers.
    cran <- c(cran, "semPlot", "circlize", "glmnet", "randomForestSRC", "xgboost")
    bioc <- c(bioc, "ComplexHeatmap")
}
github <- c(
    PheWAS="PheWAS/PheWAS@eae96fa5ed3f42e6dc5d22fad09fd7734f0df48f",
    dietaryindex="jamesjiadazhan/dietaryindex@9830a68bcfdc16204b7c2a85668548cb0b5f22d8"
)
expected_version <- c(PheWAS="0.99.6-2", dietaryindex="2.1")
cran <- unique(cran)
available <- function(x) vapply(x, requireNamespace, logical(1), quietly=TRUE)
missing_packages <- function(x) setdiff(x, rownames(utils::installed.packages()))

if (install) {
    lib <- .libPaths()[1L]
    if (file.access(lib, 2L) != 0L) stop("R library is not writable: ", lib,
        ". Activate the conda environment, or set R_LIBS_USER to an existing writable directory.")
    options(repos=c(CRAN=Sys.getenv("DIET_CRAN", "https://cloud.r-project.org")))
    workers <- suppressWarnings(as.integer(Sys.getenv("DIET_INSTALL_JOBS", "2")))
    if (is.na(workers) || workers < 1L) workers <- 2L
    options(Ncpus=workers, timeout=max(300, getOption("timeout", 60)))
    need <- missing_packages(cran)
    if (length(need)) utils::install.packages(need, lib=lib, dependencies=NA)
    need <- missing_packages(bioc)
    if (length(need)) {
        if (length(missing_packages("BiocManager"))) utils::install.packages("BiocManager", lib=lib)
        BiocManager::install(need, lib=lib, update=FALSE, ask=FALSE, dependencies=NA)
    }
    need <- missing_packages(names(github))
    if (length(need)) {
        if (length(missing_packages("remotes"))) utils::install.packages("remotes", lib=lib)
        for (pkg in need) {
            cache <- Sys.getenv("DIET_SOURCE_CACHE", "")
            archive <- file.path(cache, paste0(pkg, ".tar.gz"))
            if (nzchar(cache) && file.exists(archive)) {
                remotes::install_local(archive, lib=lib, dependencies=NA,
                    upgrade="never", build_vignettes=FALSE)
            } else remotes::install_github(github[[pkg]], lib=lib,
                dependencies=NA, upgrade="never", build_vignettes=FALSE)
        }
    }
}

required <- unique(c(recommended, cran, bioc, names(github)))
absent <- required[!available(required)]
cat(sprintf("R %s: %d/%d required packages loadable (%s).\n", getRversion(),
    length(required)-length(absent), length(required), if (optional) "core + optional" else "core"))
if (length(absent)) cat("Missing/unloadable: ", paste(absent, collapse=", "), "\n", sep="")
for (pkg in setdiff(names(github), absent)) {
    desc <- utils::packageDescription(pkg)
    actual <- desc$RemoteSha
    expected <- sub("^.*@", "", github[[pkg]])
    status <- if (is.null(actual)) "revision metadata unavailable" else if (identical(actual, expected)) "pinned revision" else "different revision retained"
    cat(sprintf("  %s %s: %s (reference version %s).\n", pkg, desc$Version, status, expected_version[[pkg]]))
    if (is.null(actual) || !identical(actual, expected)) cat("    Reference: ", github[[pkg]], "\n", sep="")
}
if (length(absent)) {
    cat("Run diet/install.sh --install", if (optional) " --with-optional" else "", " in the intended environment.\n", sep="")
    quit(save="no", status=1L)
}
