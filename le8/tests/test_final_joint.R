# Synthetic numerical tests. Run natively with Rscript; no UKB data required.
suppressPackageStartupMessages({library(dplyr);library(tibble);library(purrr);library(tidyr);library(survival)})
source("le8/f/c5_revision_validation.R")
source("le8/f/c5_joint_helpers.R")
source("le8/f/c5_joint_extensions.R")
set.seed(18);n<-1200
tr<-data.frame(eid=paste0("t",seq_len(n)),.group=paste0("g",seq_len(n)),age=rnorm(n),
  pgs__prot__X=rnorm(n))
tr$prot__X<-2*tr$pgs__prot__X+rnorm(n)
tr$bmi<-tr$prot__X+rnorm(n)
te<-tr[1:200,];te$eid<-paste0("v",1:200);te$.group<-te$eid
mp<-data.frame(feature="prot__X",pgs_feature="pgs__prot__X")
# Validation measurements/outcomes cannot change genetic decomposition fitting.
a<-c5_pgs_decompose(tr,te,"prot__X",mp,2026)
te$prot__X<-te$prot__X+1e6
b<-c5_pgs_decompose(tr,te,"prot__X",mp,2026)
stopifnot(identical(a$coefficient,b$coefficient),identical(a$test$genpart__prot__X,b$test$genpart__prot__X),
  identical(a$train$genpart__prot__X,b$train$genpart__prot__X))
stopifnot(max(abs(a$train$genpart__prot__X+a$train$remainder__prot__X-tr$prot__X))<1e-12)
# A single supported domain is enough; absent domains never manufacture proxies.
d<-c5_replicated_domains(tr,"prot__X",c("bmi","not_measured"),"age",2026)
stopifnot(length(d$ranks$bmi$feature)==1,!"not_measured"%in%names(d$ranks),d$audit$eligible[d$audit$domain=="not_measured"]==0)
# Empty proxy queues must not become an all-NS panel labelled YSplus.
stopifnot(length(c5_ys_panel(list(character()),5,c(paste0("prot__",1:4),paste0("met__",1:4))))==0)
# Unsupported or wholly missing inflammation markers stay unavailable, not low.
m<-data.frame(feature="prot__X",assay="X",layer="prot")
g<-c5_marker_strata(tr,te,m)
stopifnot(all(g$groups$GDF15_IL6_median=="unavailable"))
# Participants not covered at the landmark cannot act as event-free controls.
q<-data.frame(.time=c(8,8,3),.event=c(0,1,1),.entry=c(6,2,0))
l<-c5_landmark(q,5,10)
stopifnot(nrow(l)==1,l$.event==1,l$.time==3)
cat("FINAL synthetic numerical tests passed; this is not validation on UKB.\n")
