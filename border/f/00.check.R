required <- c("dplyr", "tidyr", "tibble", "survival", "PheWAS", "ggplot2", "ggrepel", "digest", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "))
cat("R packages OK\n")
source(file.path(Sys.getenv("BORDER_FDIR"), "0.common.R"))
cat("Follow-up end:", as.character(follow_end()), "\n")
files <- file.path(indir,"Rdata",c("all.rds","PheWAS.rds","prot.rds","met.rds"))
if(any(!file.exists(files))) stop("Missing input files: ",paste(files[!file.exists(files)],collapse=", "))
cat("Input files and follow-up settings OK\n")
if (!"packages" %in% commandArgs(TRUE)) {
  dat <- read_all() %>% prep_covars() %>% make_birth_order()
  for(y in OUTCOMES) dat <- ensure_outcome(dat,y)
  cat("Input cohort OK:",nrow(dat),"rows; eligible birth order:",sum(dat$eligible_border),"\n")
  for(x in EXPOSURES) cat(x, sum(!is.na(dat[[x]])),"nonmissing\n")
  cat("Input and outcome checks passed; no models fitted.\n")
}
