#!/usr/bin/env Rscript
# Compatibility entry point: plotting/QC now lives in gwas_post.sh compare.
a <- commandArgs(FALSE)
script <- sub('^--file=', '', a[grepl('^--file=', a)][1])
source(file.path(dirname(normalizePath(script)), '..', '0data', 'f', 'gwas_compare.R'))
