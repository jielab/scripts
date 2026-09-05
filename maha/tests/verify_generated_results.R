suppressPackageStartupMessages({library(dplyr);library(readxl)})
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)[1]
root <- dirname(dirname(normalizePath(sub("^--file=","",arg))))
drive <- if(Sys.info()[["sysname"]]=="Windows") "D:" else "/mnt/d"
out <- Sys.getenv("MAHA_OUTDIR",unset=file.path(drive,"analysis","maha"))
read_sheet <- function(path,sheet) read_excel(file.path(out,path),sheet=sheet)
for(sh in c("Fig4B_risk_DASH","Fig4C_risk_MAHA")) {
  d <- read_sheet("nhanes/Fig4.mortality_validation.out.xlsx",sh)
  stopifnot(!"sd" %in% names(d),isTRUE(all.equal(d$mean,d$risk*100)),
            isTRUE(all.equal(d$lower_plot,d$lower*100)),isTRUE(all.equal(d$upper_plot,d$upper*100)))
  message(sh,": original risk confidence limits verified for ",nrow(d)," rows")
}
for(sh in c("primary_results","last_contact_results","maha_algorithm_results")) {
  d <- read_sheet("chns/mortality_models_final.out.xlsx",sh)
  if("Model" %in% names(d)) d <- d %>% filter(Model %in% c("Model4_socioeconomic","Model5_bmi"))
  stopifnot(nrow(d)>0,all(c("beta","se","CI_low","CI_high","P") %in% names(d)))
  stopifnot(all(is.finite(d$beta)),all(is.finite(d$se)),all(d$se>0),
            isTRUE(all.equal(d$CI_low,exp(d$beta-1.96*d$se))),
            isTRUE(all.equal(d$CI_high,exp(d$beta+1.96*d$se))),
            isTRUE(all.equal(d$P,2*pnorm(-abs(d$beta/d$se)))))
  message(sh,": consistent SE/CI/P verified for ",nrow(d)," cumulative-average results")
}
source_hash <- unname(tools::md5sum(file.path(root,"f","MAHA_chns.R")))
sigs <- list.files(file.path(out,"chns","mortality_sweep_runs"),pattern="run_signature.txt",recursive=TRUE,full.names=TRUE)
stopifnot(length(sigs)>0)
for(f in sigs) stopifnot(startsWith(readLines(f,n=1),paste0(source_hash,"|")))
message("Verified ",length(sigs)," wave-cache signatures against current CHNS source")
dat <- readRDS(file.path(out,"chns","analysis_dataset.rds"))
tab <- read_sheet("chns/Table1.out.xlsx","Table1")
calc <- dat %>% group_by(Diet_group=.data[["diet.maha.3c"]]) %>%
  summarise(recomputed=mean(education_years,na.rm=TRUE),.groups="drop")
check <- tab %>% mutate(Diet_group=as.character(Diet_group)) %>%
  left_join(calc %>% mutate(Diet_group=as.character(Diet_group)),by="Diet_group")
stopifnot(isTRUE(all.equal(check$education_mean,check$recomputed)))
print(tab %>% select(Diet_group,N,education_mean))
for(stem in c("Fig3.nhanes.mortality","Fig4.chns.mortality",
              "FigS6.chns.mortality_sensitivity","FigS7.chns.high_event_sensitivity")) {
  for(suffix in c(".png",".out.xlsx")) stopifnot(file.info(file.path(out,paste0(stem,suffix)))$size>1000)
}
print(read_sheet("nhanes/Fig4.mortality_validation.out.xlsx","Fig4A_assoc") %>%
      select(any_of(c("Diet","N_total","N_event"))))
cat("PASS: regenerated NHANES and CHNS results, schooling means, source signatures, publication outputs\n")
