source(file.path(Sys.getenv("BORDER_FDIR"), "0.common.R"))
suppressPackageStartupMessages(library(PheWAS))
dat <- load_cohort()
pf <- file.path(indir,"Rdata","PheWAS.rds")
if (!file.exists(pf)) stop("Missing ",pf,". Run ukb/phe.sh with MAKE_PHEWAS=TRUE.")
phenos <- validate_ids(readRDS(pf), "id", "PheWAS.rds")
if (ncol(phenos) < 2 || any(!vapply(phenos[setdiff(names(phenos), "id")], function(v) all(is.na(v) | v %in% c(0,1)), logical(1))))
  stop("PheWAS outcomes must be binary (0/1/NA)")
# PheWAS selects logistic regression only for logical outcomes; integer 0/1
# would silently invoke linear regression. Preserve NA control exclusions.
phenos[setdiff(names(phenos), "id")] <- lapply(phenos[setdiff(names(phenos), "id")],
  function(v) if(is.logical(v)) v else as_num(v) == 1)
cores <- as.integer(Sys.getenv("BORDER_PHEWAS_CORES","4"))
# Hash actual file contents, code and model options, not merely the exposure name.
inputs <- c(file.path(rawdir,"border.cohort.rds"), pf,
            file.path(Sys.getenv("BORDER_FDIR"), c("0.common.R", "3.phewas.R")))
signature <- digest::digest(list(hashes=tools::md5sum(inputs), covars=covars,
  package=as.character(packageVersion("PheWAS"))), algo="sha256")
run_one <- function(x) {
  cache <- file.path(rawdir,paste0("phewas_",x,".rds"))
  if (file.exists(cache) && !REPLACE) {
    old <- tryCatch(readRDS(cache), error=function(e) NULL)
    if(is.list(old) && identical(old$signature, signature)) {
      message("[PheWAS] matched cache: ",x); return(old$result)
    }
  }
  message("[PheWAS] fitting ",x)
  d <- complete_model_data(dat, c("eid",x,covars))
  if(nrow(d) < 500 || length(unique(d[[x]])) < 2) stop("Insufficient complete data for PheWAS: ",x)
  g <- d[,c("eid",x)]; names(g)<-c("id",x)
  cv <- d[,c("eid",covars)]; names(cv)[1] <- "id"
  ph <- phenos[phenos$id %in% d$eid, , drop=FALSE]
  if(nrow(ph) < 500) stop("Insufficient shared IDs for PheWAS: ",x)
  r <- PheWAS::phewas(phenotypes=ph, genotypes=g, covariates=cv, cores=cores,
                     additive.genotypes=FALSE, significance.threshold="bonferroni")
  r <- as.data.frame(PheWAS::addPhecodeInfo(r))
  names(r)[names(r)=="group"] <- "category"
  r$exposure <- x; r$p_raw <- r$p; r$FDR <- p.adjust(r$p,"BH")
  r$p_underflow <- !is.na(r$p) & r$p==0
  r$logp <- -log10(pmax(r$p, .Machine$double.xmin))
  r$bonferroni_threshold <- .05/max(1,sum(is.finite(r$p)))
  saveRDS(list(signature=signature,result=r),cache); r
}
res <- bind_rows(lapply(EXPOSURES,run_one))
write_csv2(res,"phewas.csv"); saveRDS(res,file.path(rawdir,"phewas.rds"))
source(file.path(Sys.getenv("BORDER_FDIR"), "phewas_plot.R"))
plot_phewas_figures(res)
