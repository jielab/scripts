#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
root <- '/mnt/d/data.BIG/gwas/self/common'
read_ids <- function(f) {
  d <- fread(f,select=c('#FID','IID'),colClasses='character')
  paste(d[[1]],d[[2]],sep=':')
}
typed <- read_ids('/mnt/h/ukbGen/37/typ/chr1.psam')
imputed <- read_ids('/mnt/h/ukbGen/37/imp/chr1.psam')
traits <- c('height','bald','bald12','cvd_cad.Yt2e','cvd_cad.t2e')
res <- rbindlist(lapply(traits,function(t) {
  d <- fread(file.path(root,t,'gwas','regenie','analysis.phe'),colClasses=c(FID='character',IID='character'))
  ids <- paste(d$FID,d$IID,sep=':')
  data.table(trait=t,complete_case_n=nrow(d),typed_chr1_overlap=sum(ids %in% typed),imputed_chr1_overlap=sum(ids %in% imputed))
}))
stopifnot(all(res$typed_chr1_overlap>0),all(res$imputed_chr1_overlap>0))
fwrite(res,file.path(root,'input_validation.tsv'),sep='\t',quote=FALSE)
print(res)
