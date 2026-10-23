source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
source(file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R"))
source(file.path(Sys.getenv("BMONTH_FDIR"),"0.decomp_helpers.R"))
dat<-load_cohort()
pf<-phewas_input_path(validate=FALSE);provenance<-phewas_provenance(pf)
phenos<-phewas_packed_panel(provenance)
cohort_md5<-unname(tools::md5sum(file.path(rawdir,"bmonth.cohort.rds")))

# Selection is explicitly post-screening/exploratory, even for Bonferroni hits.
# Invalid fits and stale cohort/panel screens cannot become candidates.
catalog<-tibble(exposure=character(),phenotype=character(),selection_reason=character())
for(filename in c("phewas.rds","phewas_joint.rds")) {
  screen<-bm_decomp_screen(file.path(rawdir,filename),provenance,cohort_md5)
  if(!nrow(screen))next
  allowed<-if(filename=="phewas.rds")PHEWAS_EXPOSURES else c("cycle2df","month11df")
  for(x in intersect(unique(screen$exposure),allowed)) {
    d<-screen[screen$exposure==x,,drop=FALSE];d<-d[order(d$p),,drop=FALSE]
    hit<-as.character(d$phenotype[d$p<=d$bonferroni_threshold])
    top<-head(as.character(d$phenotype),10)
    selected<-unique(c(hit,top))
    catalog<-bind_rows(catalog,tibble(exposure=x,phenotype=selected,
      selection_reason=ifelse(selected %in% hit,"screen_Bonferroni_hit","exploratory_top10")))
  }
}
for(x in intersect(c("birth_north_100km","birth_east_100km"),PHEWAS_EXPOSURES))
  catalog<-bind_rows(catalog,tibble(exposure=x,phenotype=c("317","317.1","317.11"),
    selection_reason="prespecified_geography_alcohol_audit"))
catalog<-catalog %>% group_by(exposure,phenotype) %>%
  summarise(selection_reason=paste(unique(selection_reason),collapse="|"),.groups="drop") %>%
  mutate(exploratory=TRUE,available=phenotype %in% phenos$phenotypes,
         phenotype_panel_md5=provenance$panel_md5,cohort_md5=cohort_md5,
         provenance_status=provenance$validation_status)
write_csv2(catalog,"decomposition_candidates.csv")
if(!nrow(catalog))stop("No valid PheWAS results found; run phewas first")

joined_outcome<-function(phenotype) {
  yy<-bm_panel_column(phenos,phenotype)
  if(any(!is.na(yy$.y) & !yy$.y %in% c(0,1)))stop("Non-binary phenotype: ",phenotype)
  inner_join(dat,yy,by="eid")
}
attach_selection<-function(r,item) {
  r$selection_reason<-item$selection_reason;r$exploratory<-TRUE
  r$provenance_status<-provenance$validation_status
  r$phenotype_panel_md5<-provenance$panel_md5;r$cohort_md5<-cohort_md5
  r
}
fit_scalar_stage<-function(item) {
  exposure<-item$exposure;phenotype<-item$phenotype
  zz<-joined_outcome(phenotype)
  common<-bm_decomp_complete(zz,unique(c(".y",exposure,
    setdiff(covars_stage$D5_migration,exposure))))
  rows<-lapply(names(covars_stage),function(stage){
    rr<-bm_decomp_scalar(common,exposure,covars_stage[[stage]],phenotype)
    rr$stage<-stage;rr$sample_set<-"common_D5";rr$comparison_role<-"common_adjustment";rr
  })
  full<-bm_decomp_scalar(zz,exposure,covars_stage$D2_adultSES,phenotype)
  full$stage<-"D2_adultSES";full$sample_set<-"full_D2";full$comparison_role<-"sample_selection"
  r<-bind_rows(rows,full);r$exposure<-exposure
  for(stage in c("D1_core","D2_adultSES","D3_geography")) {
    base<-r[r$sample_set=="common_D5" & r$stage==stage,,drop=FALSE]
    name<-paste0("attenuation_vs_",sub("_.*$","",stage))
    r[[name]]<-bm_decomp_ratio(r$beta,base$beta,r$method,base$method)
    r[[name]][r$sample_set!="common_D5"]<-NA_real_
  }
  r$interpretation<-"Descriptive conditional log-OR change on common D5 sample; not causal mediation"
  attach_selection(r,item)
}
scalar_items<-catalog %>% filter(available,exposure %in% ALL_EXPOSURES)
scalar<-bind_rows(lapply(seq_len(nrow(scalar_items)),function(i)fit_scalar_stage(scalar_items[i,])))
if(!nrow(scalar))scalar<-tibble(stage=character(),sample_set=character(),N=integer(),cases=integer(),
  beta=double(),SE=double(),OR=double(),p=double(),exposure=character(),phenotype=character(),
  status=character(),method=character(),attenuation_vs_D1=double(),attenuation_vs_D2=double(),
  attenuation_vs_D3=double())
scalar$analysis_track<-analysis_track(scalar$exposure);scalar$exposure_unit<-exposure_unit(scalar$exposure)
write_csv2(scalar,"decomposition_scalar.csv");bm_atomic_rds(scalar,file.path(rawdir,"decomposition_scalar.rds"))

selection<-scalar %>% filter(stage=="D2_adultSES") %>%
  select(any_of(c("exposure","phenotype","sample_set","N","cases","beta","OR","SE","p",
    "method","status","sample_id_hash","selection_reason","provenance_status")))
write_csv2(selection,"decomposition_sample_selection.csv")
geo_alcohol<-scalar %>% filter(exposure %in% c("birth_north_100km","birth_east_100km"),
  phenotype %in% c("317","317.1","317.11")) %>% arrange(phenotype,exposure,sample_set,stage)
write_csv2(geo_alcohol,"audit_geography_alcohol.csv")

cycle_items<-catalog %>% filter(available,exposure %in% c("cycle2df","month11df")) %>%
  group_by(phenotype) %>% summarise(selection_reason=paste(paste(exposure,selection_reason,sep=":"),collapse="|"),
                                  .groups="drop")
fit_cycle_stage<-function(item) {
  phenotype<-item$phenotype;zz<-joined_outcome(phenotype)
  common<-bm_decomp_complete(zz,unique(c(".y","bmonth_cos","bmonth_sin",covars_stage$D5_migration)))
  rows<-lapply(names(covars_stage),function(stage){
    rr<-bm_decomp_joint(common,covars_stage[[stage]],phenotype)
    rr$stage<-stage;rr$sample_set<-"common_D5";rr$comparison_role<-"common_adjustment";rr
  })
  full<-bm_decomp_joint(zz,covars_stage$D2_adultSES,phenotype)
  full$stage<-"D2_adultSES";full$sample_set<-"full_D2";full$comparison_role<-"sample_selection"
  r<-bind_rows(rows,full);r$exposure<-"cycle2df"
  base<-r[r$sample_set=="common_D5" & r$stage=="D1_core",,drop=FALSE]
  comparable<-r$sample_set=="common_D5" & !is.na(r$test_method) &
    r$test_method=="ML_LRT" & isTRUE(base$test_method=="ML_LRT")
  r$chisq_D1<-base$chisq
  r$chisq_attenuation_vs_D1<-ifelse(comparable & is.finite(base$chisq) & base$chisq>0,
                                  1-r$chisq/base$chisq,NA_real_)
  r$statistic_ratio_comparable<-comparable
  r$interpretation<-"Joint-test statistics are not effect sizes; ratios only for matched ML LRTs; descriptive only"
  attach_selection(r,item)
}
cycle<-bind_rows(lapply(seq_len(nrow(cycle_items)),function(i)fit_cycle_stage(cycle_items[i,])))
write_csv2(cycle,"decomposition_cycle.csv");bm_atomic_rds(cycle,file.path(rawdir,"decomposition_cycle.rds"))
curves<-bind_rows(lapply(seq_len(nrow(cycle)),function(i){
  z<-bm_decomp_curve(cycle[i,]);z$phenotype<-cycle$phenotype[i];z$stage<-cycle$stage[i]
  z$sample_set<-cycle$sample_set[i];z$sample_id_hash<-cycle$sample_id_hash[i]
  z$method<-cycle$method[i];z$status<-cycle$status[i];z$exploratory<-TRUE;z
}))
write_csv2(curves,"decomposition_cycle_curves.csv")
message("[decomp] common_D5 rows retain an identical risk set across D1-D5; full_D2 is an additional sample-selection comparison. ",
        "Post-screening checks and coefficient changes are exploratory, not causal mediation.")
