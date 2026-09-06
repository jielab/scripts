#!/usr/bin/env Rscript
# A common numeric design matrix gives all methods identical covariates.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
if(length(a)!=8) stop('Expected phe cov trait type event covs cats output')
d <- fread(a[1],colClasses=c(FID='character',IID='character'))
c <- fread(a[2],colClasses=c(FID='character',IID='character'))
trait <- a[3]; type <- a[4]; event <- a[5]
csv <- function(x) if(x=='') character() else strsplit(x,',',fixed=TRUE)[[1]]
covs <- unique(c(csv(a[6]),csv(a[7]))); cats <- csv(a[7])
need <- c('FID','IID',trait,if(type=='t2e') event)
if(length(setdiff(need,names(d)))) stop('Missing phenotype columns: ',paste(setdiff(need,names(d)),collapse=','))
if(length(setdiff(c('FID','IID',covs),names(c)))) stop('Missing covariate columns')
if(anyDuplicated(d[,.(FID,IID)]) || anyDuplicated(c[,.(FID,IID)])) stop('Duplicate participant keys')
if(anyDuplicated(d$IID)) stop('IID must be unique for SAIGE')
d <- merge(d[,..need],c[,c('FID','IID',covs),with=FALSE],by=c('FID','IID'))
for(v in c(trait,if(type=='t2e') event,setdiff(covs,cats))) {
  old <- d[[v]]; num <- suppressWarnings(as.numeric(old))
  if(any(!is.na(old) & is.na(num))) stop('Nonnumeric phenotype/covariate: ',v)
  set(d,j=v,value=num)
}
d <- d[complete.cases(d)]
for(v in c(trait,if(type=='t2e') event,setdiff(covs,cats))) d <- d[is.finite(get(v))]
if(type=='bt' && !all(d[[trait]] %in% 0:1)) stop('Binary phenotype must use 0/1/NA')
if(type=='t2e' && (!all(d[[event]] %in% 0:1) || any(d[[trait]]<=0))) stop('T2E requires positive duration and 0/1 event')
if(nrow(d)<10 || uniqueN(d[[trait]])<2) stop('Insufficient phenotype variation/sample size')
for(v in cats) set(d,j=v,value=factor(d[[v]]))
variable <- covs[vapply(covs,function(v) uniqueN(d[[v]])>1,logical(1))]
if(length(variable)) {
  z <- as.data.frame(d[,..variable]); names(z)<-paste0('V',seq_along(variable))
  m <- model.matrix(~.,z)[,-1,drop=FALSE]
  # Drop dependent columns after complete-case filtering (e.g. male-only bald).
  full <- cbind(Intercept=1,m); qr0<-qr(full); keep<-sort(qr0$pivot[seq_len(qr0$rank)])
  m <- full[,setdiff(keep,1L),drop=FALSE]; colnames(m)<-paste0('C',seq_len(ncol(m)))
} else m <- matrix(nrow=nrow(d),ncol=0)
out <- cbind(d[,..need],as.data.table(m))
fwrite(out,a[8],sep='\t',na='NA',quote=FALSE)
writeLines(paste(colnames(m),collapse=','),paste0(a[8],'.covars'))
writeLines(c(paste('trait',trait),paste('type',type),paste('complete_case_n',nrow(d)),paste('design_columns',ncol(m))),paste0(a[8],'.qc.txt'))
