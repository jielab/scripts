# Run from scripts/: Rscript le8/tests/test_factorial_prediction.R
# Synthetic tests only; no UKB inputs and no clinical-performance claims.
suppressPackageStartupMessages({library(survival);library(dplyr);library(tibble);library(purrr)})
source("le8/f/c5_revision_validation.R")
source("le8/f/c5_joint_helpers.R")
source("le8/f/c5_factorial.R")
source("le8/f/c5_factorial_outputs.R")
Sys.unsetenv(c("C5_PRS_MANIFEST","C5_DISEASE_PRS_OVERLAP","C5_RISK_CUTS","C5_FACTORIAL_BOOT_ALL"))
ph<-data.frame(eid=1:10,cvd_cad.pgs=c(1:8,NA,Inf),cad.pgs=11:20)
z<-c5_prs_input("cvd_cad",ph)
stopifnot(z$enabled,nrow(z$data)==8,z$manifest$column=="cvd_cad.pgs",z$manifest$sample_overlap=="unknown")
stopifnot(!c5_prs_input("ra",ph)$enabled)
no_exact<-ph;no_exact$cvd_cad.pgs<-NULL
stopifnot(!c5_prs_input("cvd_cad",no_exact)$enabled) # never infer cad.pgs
constant<-ph;constant$cvd_cad.pgs<-1
stopifnot(!c5_prs_input("cvd_cad",constant)$enabled)
Sys.setenv(C5_DISEASE_PRS_OVERLAP="yes")
stopifnot(!c5_prs_input("cvd_cad",ph)$enabled)
Sys.unsetenv("C5_DISEASE_PRS_OVERLAP")
td<-tempfile();dir.create(td)
saveRDS(data.frame(eid=1:10,score=101:110),file.path(td,"scores.rds"))
mf<-file.path(td,"manifest.csv")
write.csv(data.frame(Y="cvd_cad",file="scores.rds",column="score",source="synthetic",build="b38",ancestry="synthetic",sample_overlap="none"),mf,row.names=FALSE)
Sys.setenv(C5_PRS_MANIFEST=mf)
stopifnot(identical(c5_prs_input("cvd_cad",ph)$data$disease_PRS,101:110))
Sys.unsetenv("C5_PRS_MANIFEST")
design<-c5_factorial_design(TRUE);small<-c5_factorial_design(FALSE)
stopifnot(nrow(design)==15,sum(design$available)==15,sum(small$available)==7)
ed<-c5_factorial_edges(design)
stopifnot(nrow(ed)==28,sum(ed$focus)==3,
  setequal(ed$reference[ed$focus],c("F_CPM","F_CPG","F_CMG")))
# Predictive contribution is exactly additive and sums to full minus Clinical.
v<-setNames(.6+.07*design$P+.02*design$M+.01*design$G+.02*(design$P&design$M),design$model)
ss<-c5_factorial_shapley(v)
stopifnot(max(abs(ss-c(P=.08,M=.03,G=.01)))<1e-12,abs(sum(ss)-(v['F_CPMG']-v['F_C']))<1e-12)
a<-c(.1,.2,.7,.9);y<-c(0,0,1,1);w<-c(1,2,3,1)
stopifnot(all(c5_factorial_reclassification(a,a,y,w,c(.2,.5))==0))
rc<-c5_factorial_reclassification(c(.05,.1,.8,.95),a,y,w)
stopifnot(rc['continuous_NRI']==2,rc['NRI_event']==1,rc['NRI_nonevent']==1,rc['IDI']>0)
Sys.setenv(C5_RISK_CUTS="0.1,0.05")
stopifnot(identical(c5_factorial_cutpoints(),c(.05,.1)))
Sys.unsetenv("C5_RISK_CUTS")
# All models deliberately have identical predictions, but differently ordered
# rows. Pairing must succeed by ID and all gains/intervals must be zero.
n<-240;d<-data.frame(eid=as.character(1:n),fold=rep(1:3,80),time=c(rep(2,90),rep(8,150)),
  event=c(rep(1,90),rep(0,150)),group=as.character(1:n),risk=seq(.01,.5,length.out=n),
  lp=seq(-2,2,length.out=n),landmark=0,horizon=5,arm="all_assays",high_training_q75=seq_len(n)>180)
pred<-bind_rows(lapply(seq_len(nrow(small)),function(i){
  if(!small$available[i])return(NULL)
  x<-d;x$model<-small$model[i];x$budget<-if(small$P[i]||small$M[i])5 else 0
  x[rev(seq_len(n)),,drop=FALSE]
}))
z<-c5_factorial_data(pred,small,5,0,"all_assays",d$eid,3)
stopifnot(ncol(z$risk)==7,identical(z$data$eid,d$eid))
stopifnot(c5_within_fold_c(d$time,d$event,rep(0,n),d$fold,5)==.5)
scaled<-d$lp+100*d$fold
stopifnot(c5_within_fold_c(d$time,d$event,scaled,d$fold,5)==c5_within_fold_c(d$time,d$event,d$lp,d$fold,5))
yin<-data.frame(eid=d$eid,.time=d$time,.event=d$event)
out<-c5_factorial_outputs(pred,yin,small,3,20,12345,5,0)
dd<-out$factorial_contrasts
stopifnot(all(dd$estimate==0),all(dd$lo==0),all(dd$hi==0),all(dd$N==n))
# An incomplete prediction vector never silently reduces comparison samples.
broken<-pred[-1,];z<-c5_factorial_data(broken,small,5,0,"all_assays",d$eid,3)
stopifnot(sum(z$audit$status=="complete")==6,ncol(z$risk)==6)
cat("Endpoint PRS / factorial synthetic numerical tests passed.\n")
