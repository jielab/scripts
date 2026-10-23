source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))

dat <- read_all() %>% merge_supplemental_covariates() %>%
  merge_move(required=TRUE) %>% make_birth_month() %>% add_retained_covariates()
birth_year_source <- attr(dat,"birth_year_source")
birth_year_approx <- attr(dat,"birth_year_approx")
dat <- prep_covars(dat) %>% make_birth_order()
outcome_sources <- setNames(vapply(OUTCOMES,function(y)
  if(all(paste0(y,c(".t2e",".Yt2e"))%in%names(dat))) "existing_time_status" else "reconstructed_dates",
  character(1)),OUTCOMES)
for(y in OUTCOMES) {
  prepared <- tryCatch(ensure_outcome(dat,y),error=identity)
  if(inherits(prepared,"error")) {
    if(bm_truth(Sys.getenv("BMONTH_CHECK_SURVIVAL","FALSE")))stop(prepared)
    outcome_sources[y] <- "deferred_until_Cox_requested"
    message("Cox outcome ",y," deferred: ",conditionMessage(prepared),
      "; exposure/PheWAS preparation continues")
  } else dat <- prepared
}

keep <- unique(c(
  "eid","birth_month","birth_year_bm","birth_north","birth_east",
  "date_attend","birth_date","birth_year","date_lost","date_death",
  "age_fractional","assessment_year","assessment_month",
  "ethnicity","ethnic.c","sex.b","sex.f","smoke_status","smoke.c","smoke_quit_age",
  "bmi","waist","hip","whr","ldl","nonhdl.c","hba1c_ngsp","drug.lipid","drug.dm","drug.htn",
  csv_env("BMONTH_EXTRA_COVARS"), csv_env("BMONTH_FOLLOWUP_EXTRA_COVARS"),
  csv_env("BMONTH_OMICS_COVARIATES"),
  ALL_EXPOSURES,
  "cnt_older","cnt_brother","cnt_sister","adopted","adopted_num",
  "sibling_n","sibship_size","sibship_size_cap","sibship_size_cat",
  "birth_rank","impossible_rank","eligible_border",
  "tdi.birth","tdi.birth.nn","move_km","move_log1p",
  "lon.birth","lat.birth","lon.home","lat.home",
  "bmonth_cos_x_north","bmonth_sin_x_north",
  "bmonth_cos_x_east","bmonth_sin_x_east",
  "bmonth_cos_x_move","bmonth_sin_x_move",
  unlist(lapply(OUTCOMES,function(y)outcome_date_columns(dat,y))),
  unlist(covars_stage,use.names=FALSE),
  unlist(lapply(OUTCOMES,function(y)paste0(y,c(".t2e",".Yt2e"))))
))
requested_extra<-unique(c(csv_env("BMONTH_EXTRA_COVARS"),csv_env("BMONTH_FOLLOWUP_EXTRA_COVARS"),csv_env("BMONTH_OMICS_COVARIATES")))
if(length(setdiff(requested_extra,names(dat))))stop("Requested supplemental columns are missing: ",paste(setdiff(requested_extra,names(dat)),collapse=", "))
keep <- intersect(keep,names(dat))
cohort <- dat[,keep,drop=FALSE]
attr(cohort,"bmonth_metadata") <- list(
  code_version=BMONTH_CODE_VERSION,
  input_md5=unname(tools::md5sum(file.path(indir,"Rdata","all.rds"))),
  move_md5=unname(tools::md5sum(file.path(indir,"Rdata","move.rds"))),
  extra_covariates=requested_extra,
  supplemental_input=Sys.getenv("BMONTH_SUPPLEMENTAL_COVARIATES",""),
  supplemental_md5=if(nzchar(Sys.getenv("BMONTH_SUPPLEMENTAL_COVARIATES","")))
    unname(tools::md5sum(Sys.getenv("BMONTH_SUPPLEMENTAL_COVARIATES"))) else NA_character_,
  fractional_age_note="Derived from existing birth_date (UKB month-15 proxy); not exact birth date",
  input=file.path(indir,"Rdata","all.rds"),
  move=file.path(indir,"Rdata","move.rds"),
  outcomes=OUTCOMES,
  outcome_sources=outcome_sources,
  follow_end=if(any(outcome_sources=="reconstructed_dates") ||
    nzchar(Sys.getenv("BMONTH_FOLLOW_END","")))as.character(follow_end())else NA_character_,
  birth_year_source=birth_year_source,
  birth_year_approx=birth_year_approx,
  covars_primary=covars,
  covars_stage=covars_stage,
  exploratory_traits=EXPLORATORY_TRAITS,
  exploratory_covariates=exploratory_covars,
  exploratory_family_size_covariates=exploratory_family_covars,
  note_birth_order=paste("Exploratory birth-order traits exclude adopted==1, invalid sibling counts,",
    "impossible ranks and singletons; unknown adoption status is retained.",
    "Eligibility does not restrict the birth-month or geography analyses.",
    "Primary D2 birth-order models do not adjust family size.",
    "Sensitivity models for border1n, border_rank and border_relative additionally adjust",
    "categorical sibship size (2,3,4,5,6+); border12 does not require this because family size is fixed at 2."),
  note_tdi_birth=paste(
    "tdi.birth is the value created by ukb/f/phe.R: current-home TDI spatial surface",
    "mapped by nearest neighbour to birthplace. It is a birthplace deprivation proxy,",
    "not a historical census Townsend score measured in the participant's birth year.")
)
saveRDS(cohort,file.path(rawdir,"bmonth.cohort.rds"))

nn_km <- if("tdi.birth.nn"%in%names(dat)) dat$tdi.birth.nn/1000 else rep(NA_real_,nrow(dat))
safe_median <- function(x){x<-x[is.finite(x)];if(length(x))median(x) else NA_real_}
safe_q95 <- function(x){x<-x[is.finite(x)];if(length(x))as.numeric(quantile(x,.95)) else NA_real_}
qc <- bind_rows(
  tibble(metric="all.rds N",value=nrow(dat)),
  tibble(metric="birth_month nonmissing",value=sum(!is.na(dat$birth_month))),
  tibble(metric="birth_north nonmissing",value=sum(!is.na(dat$birth_north))),
  tibble(metric="birth_east nonmissing",value=sum(!is.na(dat$birth_east))),
  tibble(metric="tdi.birth nonmissing",value=sum(!is.na(dat$tdi.birth))),
  tibble(metric="move_km nonmissing",value=sum(!is.na(dat$move_km))),
  tibble(metric="tdi.birth nearest-neighbour distance median km",
         value=safe_median(nn_km)),
  tibble(metric="tdi.birth nearest-neighbour distance p95 km",
         value=safe_q95(nn_km)),
  bind_rows(lapply(BASE_EXPOSURES,function(x)
    tibble(metric=paste0(x," nonmissing"),value=sum(!is.na(dat[[x]])))))
)
write_csv2(qc,"cohort_qc.csv")
write_csv2(birth_order_qc(dat),"birth_order_qc.csv")
write_csv2(count(dat,birth_month,name="N"),"birth_month_counts.csv")
writeLines(capture.output(str(attr(cohort,"bmonth_metadata"))),
           file.path(rawdir,"cohort_metadata.txt"))
cat("Saved ",file.path(rawdir,"bmonth.cohort.rds"),"\n",sep="")
print(qc,n=Inf)
print(birth_order_qc(dat),n=Inf)
