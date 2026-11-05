# Missing disease GWAS must preserve independent C2 work and never become a
# permanently completed genetic fit. Uses synthetic objects, no cohort reads.
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
suppressPackageStartupMessages({library(dplyr); library(tibble)})
e <- new.env(parent=globalenv())
wanted <- c("run_c2_without_summary_genetics", "run_c3_layer", "le8_empty_mr", "le8_completed_results")
extract <- function(x) {
    if (missing(x) || !is.call(x)) return()
    if (identical(x[[1]],as.name("<-")) && length(x)==3 && is.symbol(x[[2]]) && as.character(x[[2]]) %in% wanted) {
        eval(x,e); return()
    }
    if (identical(x[[1]],as.name("function"))) return()
    for (z in as.list(x)[-1]) extract(z)
}
for (n in c("0.common.R","c2.cause.R","c3.coloc.R")) for (x in parse(file.path(root,"f",n))) extract(x)
stopifnot(all(vapply(wanted,exists,logical(1),envir=e,inherits=FALSE)))
tmp <- tempfile("missing-genetics-"); dir.create(tmp)
e$`%||%` <- function(x,y) if (is.null(x)) y else x
e$le8_job_dir <- function(outdir,module) file.path(outdir,module)
e$write_raw_csv <- function(x,file,dir) write.csv(x,file.path(dir,file),row.names=FALSE)
e$write_xlsx2 <- function(...) NULL
e$save_plot <- e$plot_individual_genetic_decomposition <- e$plot_component_leadtime <- function(...) NULL
e$finalize_outputs <- e$le8_begin_analysis <- e$le8_finish_analysis <- function(...) NULL
e$le8_table_archive_files <- function(...) character()
e$setwd2 <- function(...) NULL
e$C2_CODE_VERSION <- e$C3_CODE_VERSION <- "synthetic"
e$RUN_MRlink2 <- e$RUN_Dandelion <- "Top"
e$module_meta <- function(layer,extra) c(list(trait="toy",layer=layer,module=e$LE8_JOB,generated="synthetic"),extra)
e$LE8_REUSE_RESULTS <- FALSE
e$C3_MR_FDR <- .05; e$MAX_FEATURES <- 0; e$MAX_LOCI_PER_FEATURE <- 3
e$out.prot <- file.path(tmp,"toy","prot"); e$out.met <- file.path(tmp,"toy","met")
e$dir.X <- e$dir.met.gwas <- file.path(tmp,"absent_qtl")
e$Y <- "toy"
e$get_y_gwas_file <- function(...) NA_character_
e$decomp_calls <- 0L
e$run_individual_genetic_decomposition <- function(...) {
    e$decomp_calls <- e$decomp_calls+1L
    list(status=tibble(status="ok"),summary=tibble(feature="protein",effect=.2),trajectory=tibble())
}
e$LE8_JOB <- "c2_cause"
rd <- file.path(e$out.prot,e$LE8_JOB); dir.create(rd,recursive=TRUE)
c2 <- e$run_c2_without_summary_genetics("protein",e$out.prot,rd,"disease GWAS missing","mode","pgs")
stopifnot(e$decomp_calls==1, nrow(c2$individual_decomposition$summary)==1, nrow(c2$MR)==0,
          c2$subanalysis_status$status[c2$subanalysis_status$analysis=="individual_PGS_decomposition"]=="completed")
e$LE8_JOB <- "c3_coloc"
c3 <- e$run_c3_layer("protein")
stopifnot(nrow(c3$summary)==0,nrow(c3$subanalysis_status)==3,
          all(c3$subanalysis_status$status=="unavailable_genetic_inputs"))
audit <- e$le8_completed_results(tmp,"toy","prot",c("c2_cause","c3_coloc"))
stopifnot(all(audit$status=="completed_with_unavailable_genetics"), !any(audit$status %in% c("completed","cached")))
# A partial result without its explicit availability export fails validation.
unlink(file.path(rd,"c2.subanalysis_status.csv"))
audit <- e$le8_completed_results(tmp,"toy","prot","c2_cause")
stopifnot(audit$status=="incomplete")
unlink(tmp,recursive=TRUE)
cat("Missing-GWAS checks passed: independent C2 work retained, C3 unavailable, future inputs rechecked.\n")
