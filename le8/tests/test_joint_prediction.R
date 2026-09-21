# Run from scripts/: Rscript le8/tests/test_joint_prediction.R
# Tests require R + survival/glmnet/dplyr/tibble/purrr; no UKB data are used.
suppressPackageStartupMessages({library(survival);library(glmnet);library(dplyr);library(tibble);library(purrr)})
source("le8/f/c0_baseline.R")
source("le8/f/c5_revision_validation.R")
source("le8/f/c5_joint_helpers.R")
Sys.unsetenv(c("LE8_BASELINE_MAP","LE8_BASELINE_MED_COLUMNS"))
d<-data.frame(eid=1:4,date_attend=as.Date("2010-01-01"),bmi=25,bb_TC=5,bb_HDL=1.5,
  bb_HBA1C=35,sbp=120,dbp=70,drug.big3_i0=c("-7","-1","1","3"),drug.big3_i1="1;2;3",
  fod_icd10_t2dm=as.Date(c("2015-01-01",NA,NA,"2008-01-01")))
a<-le8_rebuild_baseline(d)
stopifnot(a$drug.lipid[1]==0,is.na(a$drug.lipid[2]),a$drug.lipid[3]==1,
  a$dm.yes[1]==0,a$dm.yes[4]==1,is.na(a$hba1c.pts[2]))
# Alter future visits and future diagnosis dates: baseline predictors must not change.
future<-d;future$drug.big3_i1<-"-7";future$fod_icd10_t2dm[1]<-as.Date("2020-01-01")
b<-le8_rebuild_baseline(future)
for(v in c("drug.lipid","drug.dm","drug.htn","dm.yes","hba1c.pts","nonhdl.pts","bp.pts"))stopifnot(identical(a[[v]],b[[v]]))
# Remove first-visit medication data: unknown is never silently treated as no use.
d$drug.big3_i0<-NULL;stopifnot(all(is.na(le8_rebuild_baseline(d)$drug.lipid)))
# Landmark evaluation retains controls beyond horizon and excludes earlier events.
dd<-data.frame(eid=1:120,.time=c(rep(2,20),rep(7,40),rep(12,60)),.event=c(rep(1,60),rep(0,60)))
ll<-c5_landmark(dd,5,10);iw<-le8_ipcw(ll$.time,ll$.event,5)
stopifnot(nrow(ll)==100,iw$status=="ok",iw$N_case==40,iw$N_control==60)
# Tied predictions give AUC .5; correctly ordered predictions give 1.
stopifnot(abs(le8_weighted_auc(rep(.1,100),iw$y,iw$w)-.5)<1e-12,
  le8_weighted_auc(iw$y,iw$y,iw$w)==1)
# Group assignment and modality budgets remain identical across panel paradigms.
g<-rep(letters,each=3);f<-c5_group_folds(g,5,10)
stopifnot(all(vapply(split(f,g),function(x)length(unique(x))==1,logical(1))))
q<-list(paste0("prot__",1:10),paste0("met__",1:10))
for(k in c(5,10))for(ns in list(NULL,unlist(q))){
  z<-c5_ys_panel(q,k,ns);stopifnot(length(z)==k,sum(startsWith(z,"prot__"))==ceiling(k/2),!anyDuplicated(z))
}
# Preprocessing uses training data only, including missing-value imputation.
tr<-data.frame(x=c(1,2,3,NA),z=c(0,1,0,1));te<-data.frame(x=c(NA,1e9),z=c(1,0))
x<-le8_prepare_prediction_matrix(tr,te,c("x","z"))
stopifnot(x$audit$center[x$audit$variable=="x"]==2,x$test[1,"x"]==0)
# Frozen identical risk vectors must yield zero paired increments.
pp<-data.frame(eid=1:100,fold=rep(1:5,20),time=ll$.time,event=ll$.event,
  group=as.character(1:100),risk=seq(.01,.5,length.out=100))
z<-c5_paired_delta(pp,pp,5,20,2026)
stopifnot(z$delta_AUC==0,z$delta_Brier==0,z$AUC_lo==0,z$AUC_hi==0)
# Proxy held-out targets cannot affect learned coefficients or predictions.
set.seed(720);n<-650
tr<-data.frame(.group=paste0("t",1:n),x1=rnorm(n),x2=rnorm(n))
tr$bmi<-25+2*tr$x1+rnorm(n)
te<-data.frame(.group=paste0("v",1:100),x1=rnorm(100),x2=rnorm(100),bmi=0)
p1<-c5_gaussian_proxy(tr,te,"bmi",c("x1","x2"),2026)
te$bmi<-1e8
p2<-c5_gaussian_proxy(tr,te,"bmi",c("x1","x2"),2026)
stopifnot(identical(p1$pred,p2$pred),identical(p1$coefficient,p2$coefficient))
cat("Joint prediction integrity tests passed. These are synthetic tests, not UKB validation.\n")
