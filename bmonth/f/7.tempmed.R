source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
source(file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R"))
source(file.path(Sys.getenv("BMONTH_FDIR"),"0.decomp_helpers.R"))
dat<-load_cohort()
tf<-file.path(rawdir,"bmonth_temperature.rds")
if(!file.exists(tf))stop("Run `bmonth.sh temp --dir-temp ...` first.")
temp<-validate_ids(readRDS(tf),"eid","bmonth_temperature.rds")
cohort_md5<-unname(tools::md5sum(file.path(rawdir,"bmonth.cohort.rds")))
if(!identical(attr(temp,"temperature_metadata")$cohort_md5,cohort_md5))
  stop("Temperature file does not match this cohort; rerun bmonth.sh temp first")
tvars<-intersect(c("temp_conception_mean","temp_birth_mean",
                   "temp_conception_anom","temp_birth_anom"),names(temp))
tvars<-tvars[vapply(temp[tvars],function(x)sum(is.finite(as_num(x)))>=1000,logical(1))]
if(!length(tvars))stop("No usable temperature variables")
dat<-left_join(dat,temp[,c("eid",tvars),drop=FALSE],by="eid")
pf<-phewas_input_path(validate=FALSE);provenance<-phewas_provenance(pf)
phenos<-phewas_packed_panel(provenance)

# All selected analyses remain exploratory. Top-ten candidates are not silently
# represented as statistically established birth-month associations.
catalog<-tibble(phenotype=character(),selection_reason=character())
if(length(TEMPMED_PHENOS)){
  catalog<-tibble(phenotype=TEMPMED_PHENOS,selection_reason="explicit_requested_candidate")
} else for(filename in c("phewas.rds","phewas_joint.rds")) {
  screen<-bm_decomp_screen(file.path(rawdir,filename),provenance,cohort_md5)
  if(!nrow(screen))next
  allowed<-if(filename=="phewas.rds")c("bmonth_cos","bmonth_sin","bmonth_winter_summer","bmonth_spring_autumn") else
    c("cycle2df","month11df")
  screen<-screen[screen$exposure %in% allowed,,drop=FALSE]
  screen<-screen[order(screen$p),,drop=FALSE]
  hit<-as.character(screen$phenotype[screen$p<=screen$bonferroni_threshold])
  selected<-unique(c(hit,head(as.character(screen$phenotype),10)))
  if(length(selected))catalog<-bind_rows(catalog,tibble(phenotype=selected,
    selection_reason=ifelse(selected %in% hit,"screen_Bonferroni_hit","exploratory_top10")))
}
catalog<-catalog %>% group_by(phenotype) %>%
  summarise(selection_reason=paste(unique(selection_reason),collapse="|"),.groups="drop") %>%
  mutate(exploratory=TRUE,available=phenotype %in% phenos$phenotypes,
         provenance_status=provenance$validation_status,
         phenotype_panel_md5=provenance$panel_md5,cohort_md5=cohort_md5)
write_csv2(catalog,"temperature_candidates.csv")
candidates<-catalog$phenotype[catalog$available]
if(!length(candidates))stop("No valid/requested birth-month phenotypes available for temperature decomposition")
joined_outcome<-function(phenotype) {
  yy<-bm_panel_column(phenos,phenotype)
  if(any(!is.na(yy$.y) & !yy$.y %in% c(0,1)))stop("Non-binary phenotype: ",phenotype)
  inner_join(dat,yy,by="eid")
}
add_temperature_meta<-function(r,phenotype,mediator) {
  r$phenotype<-phenotype;r$mediator<-mediator
  r$selection_reason<-catalog$selection_reason[match(phenotype,catalog$phenotype)]
  r$exploratory<-TRUE;r$provenance_status<-provenance$validation_status
  r$phenotype_panel_md5<-provenance$panel_md5;r$cohort_md5<-cohort_md5
  r
}

# Same participant/risk-set sample in the before/after temperature models.
# Keep the actual joint coefficients/covariance and descriptive cycle amplitude.
cycle_models<-list();cycle_curves<-list()
cycle_one<-function(yname,m) {
  z<-bm_decomp_complete(joined_outcome(yname),c(".y","bmonth_cos","bmonth_sin",m,covars))
  a<-bm_decomp_joint(z,covars,yname)
  b<-bm_decomp_joint(z,unique(c(covars,m)),yname)
  a$temperature_adjustment<-"before";b$temperature_adjustment<-"after"
  models<-bind_rows(a,b);models<-add_temperature_meta(models,yname,m)
  cycle_models[[paste(yname,m,sep="|")]]<<-models
  cycle_curves[[paste(yname,m,sep="|")]]<<-bind_rows(lapply(seq_len(nrow(models)),function(i){
    out<-bm_decomp_curve(models[i,]);out$temperature_adjustment<-models$temperature_adjustment[i]
    out$sample_id_hash<-models$sample_id_hash[i];out$method<-models$method[i];out$status<-models$status[i]
    add_temperature_meta(out,yname,m)
  }))
  comparable<-isTRUE(a$test_method=="ML_LRT") && isTRUE(b$test_method=="ML_LRT") &&
    isTRUE(a$df==b$df) && identical(a$sample_id_hash,b$sample_id_hash)
  r<-tibble(N=nrow(z),cases=sum(z$.y),sample_id_hash=bm_decomp_sample_key(z),
    cycle_chisq_before=a$chisq,cycle_p_before=a$p,cycle_chisq_after=b$chisq,cycle_p_after=b$p,
    method_before=a$method,method_after=b$method,status_before=a$status,status_after=b$status,
    test_method_before=a$test_method,test_method_after=b$test_method,
    cycle_amplitude_before=a$cycle_amplitude_logOR,cycle_amplitude_after=b$cycle_amplitude_logOR,
    statistic_ratio_comparable=comparable,
    chisq_attenuation=if(comparable && is.finite(a$chisq) && a$chisq>0)1-b$chisq/a$chisq else NA_real_,
    amplitude_attenuation=bm_decomp_ratio(b$cycle_amplitude_logOR,a$cycle_amplitude_logOR,b$method,a$method),
    interpretation="Descriptive same-sample adjustment; test statistics are not effect sizes or causal mediated proportions")
  add_temperature_meta(r,yname,m)
}
cycle<-bind_rows(lapply(candidates,function(y)bind_rows(lapply(tvars,function(m)cycle_one(y,m)))))
write_csv2(cycle,"temperature_cycle_attenuation.csv")
bm_atomic_rds(cycle,file.path(rawdir,"temperature_cycle_attenuation.rds"))
write_csv2(bind_rows(cycle_models),"temperature_cycle_models.csv")
write_csv2(bind_rows(cycle_curves),"temperature_cycle_curves.csv")

# One common sample for linear X->M and logistic total/direct models. Bootstrap
# repeats the same deterministic MLE/BR fitting rule, without rescreening.
B<-as.numeric(Sys.getenv("BMONTH_MED_BOOT","100"));set.seed(12345)
scalar<-tibble()
for(yname in candidates)for(x in c("bmonth_cos","bmonth_sin"))for(m in tvars) {
  r<-bm_temp_product(joined_outcome(yname),x,m,covars,B)
  r$exposure<-x;r<-add_temperature_meta(r,yname,m)
  r$interpretation<-"Exploratory linear-by-logistic coefficient product; not a natural causal indirect effect"
  scalar<-bind_rows(scalar,r)
  # Publish completed combinations before starting the next expensive bootstrap.
  write_csv2(scalar,"temperature_mediation_partial.csv")
  bm_atomic_rds(scalar,file.path(rawdir,"temperature_mediation_partial.rds"))
}
if(nrow(scalar)) {
  scalar$eligible_bootstrap_test<-scalar$eligible & B>=50L
  scalar$FDR_indirect<-bm_planned_adjust(scalar$p_boot,scalar$eligible_bootstrap_test)
  scalar$correction_scope<-"Selected exploratory candidate combinations; not adjusted for phenotype selection"
}
write_csv2(scalar,"temperature_mediation.csv")
bm_atomic_rds(scalar,file.path(rawdir,"temperature_mediation.rds"))
writeLines(c(
  "All temperature results here are exploratory descriptive adjustment or coefficient-product analyses, not causal mediation.",
  "Candidate provenance and selection reasons are saved; top-ten candidates may have no multiple-testing-corrected birth-month association.",
  "Before/after annual-cycle models use the same sample. Joint cosine/sine coefficients and covariance are exported with pointwise monthly contrasts.",
  "A chi-square ratio is reported only when BOTH tests are valid ordinary ML LRTs with the same df. It is not an effect-size attenuation or mediated fraction.",
  "Cycle amplitude and selected peak month are descriptive; near a null cycle the phase/peak is unstable.",
  "Scalar X->M and total/direct outcome fits share one sample; MLE/BR handling is identical in bootstrap replicates. Bootstrap does not repeat phenotype selection.",
  "The bootstrap FDR is confined to selected candidate combinations; it does not establish post-selection or causal inference.",
  "Temperature derives from month/year and birthplace; birthplace is an imperfect gestational-residence proxy and conception timing is approximate.",
  "Within-season weather anomalies, prespecified birth-year/geography control and gestational-window sensitivity remain a staged scientific follow-up, not an implemented causal design.",
  "Postbirth migration remains outside primary covariates."
),file.path(rawdir,"temperature_methods.txt"))
message("[tempmed] saved robust exploratory cycle and coefficient-product results, with component-model diagnostics")
