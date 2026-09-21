# Run on the analysis host: Rscript le8/tests/test_prediction_audit.R
# Synthetic checks, not evidence for CAD model performance.
suppressPackageStartupMessages({library(dplyr);library(purrr);library(survival);library(glmnet)})
arg<-commandArgs(trailingOnly=FALSE)
script<-sub("^--file=","",arg[grepl("^--file=",arg)])[1]
fdir<-file.path(dirname(normalizePath(script)),"..","f")
source(file.path(fdir,"c5_revision_validation.R"))
source(file.path(fdir,"c4_focus_extensions.R"))
set.seed(17)
n<-1000;x<-rnorm(n);age<-rnorm(n)
d<-data.frame(age,x,x2=x+rnorm(n,sd=1e-7))
t<-rexp(n,rate=exp(.4*x+.2*age)/10);censor<-runif(n,3,15)
d$time<-pmin(t,censor);d$event<-as.numeric(t<=censor)
tr<-d[1:800,];te<-d[801:1000,]
fit<-le8_fit_budget_model(tr,te,"age",c("x","x2"),"time","event","ridge",seed=19)
stopifnot(fit$status=="ok",all(is.finite(le8_risk_at(fit,5))))
# Changing validation outcomes cannot change coefficients or preprocessing.
changed<-te;changed$event<-1-changed$event;changed$time<-changed$time*3
fit2<-le8_fit_budget_model(tr,changed,"age",c("x","x2"),"time","event","ridge",seed=19)
stopifnot(isTRUE(all.equal(fit$coefficient,fit2$coefficient)),identical(fit$preprocess,fit2$preprocess))
# With zero linear predictor, a tied-risk-set Breslow increment is hand-verifiable.
bh<-le8_breslow_hazard(c(1,1,2,3),c(1,1,0,1),rep(0,4))
stopifnot(isTRUE(all.equal(bh$hazard,c(.5,.5,1.5))))
# Affine copies have no distinct Cox information; include intercept in rank check.
tr$affine<-2*tr$x+3;te$affine<-2*te$x+3
mx<-le8_prepare_prediction_matrix(tr,te,c("x","affine"))
stopifnot(ncol(mx$train)==1L)
# YSbalanced never invents an unvalidated sleep marker to fill a pillar.
mm<-tibble(feature=c("A","B","C","D"),component=c("bmi","bmi","smoke","sleep"),
  selected=c(TRUE,TRUE,TRUE,FALSE),r1=c(.5,.4,.3,.9),r2=c(.4,.3,.2,.8))
stopifnot(identical(focus_balanced_panel(mm,3,c("bmi","smoke","sleep")),c("A","C","B")))
cat("Synthetic prediction and LE8 selection checks passed\n")
