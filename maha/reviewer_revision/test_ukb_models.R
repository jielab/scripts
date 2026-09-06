suppressPackageStartupMessages({library(dplyr);library(tidyr);library(survival);library(ggplot2)})
source("D:/scripts/0f/0phe.f.R")
source("D:/scripts/0f/assoc.f.R")
for (x in parse("f/MAHA_ukb.R", encoding="UTF-8")) {
  if (is.call(x) && identical(x[[1]],as.name("<-")) && is.call(x[[3]]) && identical(x[[3]][[1]],as.name("function"))) eval(x)
}
d <- readRDS("reviewer_revision/ukb_preflight_data.rds")
covs <- c("age.landmark", "sex.f", "center", "tdi", "PC1", "PC2", "bmi.pts", "bp.pts", "nonhdl.pts", "smoke.pts", "pa.pts", "sleep.pts")
maha_3c <- "diet.maha.3c"
for (nm in c("maha","dash","mind","medi24")) {
  d[[paste0("diet.",nm,".pts")]] <- std_10_90(d[[paste0("diet.",nm,".sum")]])
  d[[paste0("diet.",nm,".3c")]] <- f3c(d[[paste0("diet.",nm,".sum")]])
}
cp <- readRDS("D:/data/ukb/phe/Rdata/diet.maha.rds")
cp$eid <- as.character(cp$eid)
compcols <- paste0(c("protein","dairy","veg","fruit","wholegrain","fat","upf","alcohol","sodium"),".maha")
d$maha_n_components <- rowSums(!is.na(cp[match(d$eid,cp$eid),compcols]))
models <- assoc_reg(d,paste0("diet.",c("maha","dash","mind","medi24"),".pts"),covs,"death",type="t2e")
stopifnot(all(is.finite(models$estimate)),all(is.na(models$error)),all(models$N_event > 100))
write.csv(models,"reviewer_revision/ukb_model_preflight.csv",row.names=FALSE)
Y <- "death"
# Numeric sex is reconstructed solely for this preflight from its stored factor.
stopifnot(all(tolower(as.character(na.omit(d$sex.f))) %in% c("male","female")))
d$sex <- ifelse(tolower(as.character(d$sex.f))=="male",1,0)
basic <- c("age.landmark","sex","tdi","PC1","PC2")
a <- extract_joint_risk(d,"death.t2e","death.Yt2e","diet.dash.3c","diet.maha.3c",basic,panel="preflight")
b <- run_joint_risk_one(d,"death",covs)
stopifnot(nrow(a)==9,nrow(b)==9,all(is.finite(a$risk)),all(is.finite(b$risk10)),
  length(unique(a$N_model))==1,length(unique(b$N_model))==1,all(a$risk>=0 & a$risk<=1))
cc <- complete.cases(d[,c("death.t2e","death.Yt2e",covs,"diet.maha.pts")])
cat("Real-data mortality and both joint-risk functions PASSED.\n")
cat("Common MAHA model N:",sum(cc),"; complete-nine N:",sum(cc & d$maha_n_components==9),"\n")
write.csv(a,"reviewer_revision/ukb_basic_risk_preflight.csv",row.names=FALSE)
write.csv(b,"reviewer_revision/ukb_full_risk_preflight.csv",row.names=FALSE)
