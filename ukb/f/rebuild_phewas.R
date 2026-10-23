#!/usr/bin/env Rscript
# Standalone rebuild from existing code0/date0 exports; no regeneration of all.rds.
args <- commandArgs(TRUE)
usage <- function() cat(paste(
  "Rscript --vanilla rebuild_phewas.R --code FILE --date FILE --roster FILE --dir-out DIR",
  "  [--coverage-definition-file FILE] [--uncertainty FILE] [--replace]",
  "Roster: id, phewas_observable (logical), sex (F/M/NA).",
  "Coverage definition is required via the text file or roster RDS attribute.",
  "See docs/PHENOTYPE_REBUILD.md before constructing linkage observability.",sep="\n"),"\n")
if (!length(args) || any(args %in% c("-h","--help"))) {usage();quit(save="no",status=0L)}
opt <- list(replace=FALSE)
i <- 1L
while (i<=length(args)) {
  nm <- args[i]
  if (nm=="--dir-out") nm <- "--out"
  if (nm=="--replace") {opt$replace<-TRUE;i<-i+1L;next}
  if (!nm %in% c("--code","--date","--roster","--out","--coverage-definition-file","--uncertainty") || i==length(args))
    stop("Unknown or incomplete argument: ",nm)
  opt[[substring(nm,3)]] <- args[i+1L];i<-i+2L
}
for (nm in c("code","date","roster","out")) if (is.null(opt[[nm]])) stop("Missing --",nm)
script <- sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])[1]
source(file.path(dirname(normalizePath(script)),"phewas_icd10.R"))
coverage <- if (!is.null(opt[["coverage-definition-file"]]))
  paste(readLines(opt[["coverage-definition-file"]],warn=FALSE),collapse="\n") else NULL
manifest <- ukb_rebuild_phewas(opt$code,opt$date,opt$roster,opt$out,coverage,overwrite=opt$replace,
                              uncertainty_file=opt$uncertainty)
cat("Validated ICD10 panel saved: ",file.path(opt$out,"PheWAS.rds"),
    " | N=",manifest$n_roster,"; phenotypes=",manifest$n_phenotypes,"\n",sep="")
