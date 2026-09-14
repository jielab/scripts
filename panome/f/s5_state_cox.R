# Independent-test association estimates for discovery-defined states; exploratory.
a<-commandArgs(trailingOnly=TRUE);d<-read.csv(a[1],check.names=FALSE);out<-a[2]
if(!requireNamespace("survival",quietly=TRUE))stop("R package survival is required")
d<-d[d$split=="test",,drop=FALSE]
d$state<-factor(d$state)
covs<-grep("^cov_",names(d),value=TRUE)
# Remove constants in this subset; fail rather than fabricate unstable HRs.
covs<-covs[vapply(d[,covs,drop=FALSE],function(x)sd(x)>1e-8,logical(1))]
f<-as.formula(paste("survival::Surv(time,event) ~ state",if(length(covs))paste("+",paste(covs,collapse="+"))else ""))
fit<-survival::coxph(f,data=d,x=TRUE)
s<-summary(fit)
r<-data.frame(term=rownames(s$coefficients),s$coefficients,s$conf.int,check.names=FALSE)
r$fdr<-p.adjust(r[["Pr(>|z|)"]],method="BH")
write.csv(r,file.path(out,"test_state_cox.csv"),row.names=FALSE)
z<-survival::cox.zph(fit)
write.csv(data.frame(term=rownames(z$table),z$table),file.path(out,"test_state_ph_diagnostics.csv"),row.names=FALSE)
saveRDS(fit,file.path(out,"test_state_cox.rds"))
