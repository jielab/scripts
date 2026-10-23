source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
source(file.path(Sys.getenv("BMONTH_FDIR"),"mediation.R"))
source(file.path(Sys.getenv("BMONTH_FDIR"),"med_checkpoint.R"))
base<-load_cohort()
if(!length(EXPLORATORY_TRAITS)){
  empty<-data.frame(exposure=character(),mediator=character(),profile=character(),N=integer(),
    a=double(),p_a=double(),FDR=double(),status_a=character())
  for(nm in c("protwas_family_size.csv","immune_score_results.csv","protwas_model_comparison.csv"))write_csv2(empty,nm)
  med_atomic(empty,file.path(rawdir,"protwas_family_size.rds"))
  write_csv2(data.frame(exposure=NA_character_,status="skipped_birth_order_disabled"),"protwas_followup_completion.csv")
  write_csv2(data.frame(code_version=BMONTH_CODE_VERSION,status="skipped_birth_order_disabled"),"protwas_followup_manifest.csv")
  message("[omics-followup] skipped: birth-order traits disabled")
  quit(save="no",status=0L)
}
pf<-file.path(indir,"Rdata","prot.rds")
if(!file.exists(pf))stop("Missing ",pf)
om<-validate_ids(readRDS(pf),"eid","prot.rds")
features<-setdiff(names(om),"eid")
features<-features[vapply(om[features],is.numeric,logical(1))]
if(length(intersect(features,names(base))))stop("Protein names overlap retained cohort columns")
d0<-inner_join(base,om[,c("eid",features),drop=FALSE],by="eid")
if(!nrow(d0)||!length(features))stop("No overlapping participants/numeric proteins")
set.seed(20261020)
code_files<-file.path(Sys.getenv("BMONTH_FDIR"),c("0.common.R","0.analysis_contract.R",
  "9.omics_followup.R","mediation.R","med_checkpoint.R","phewas_engine.R"))
extra<-csv_env("BMONTH_OMICS_COVARIATES")
if(length(setdiff(extra,names(d0))))stop("Requested omics covariates are missing: ",paste(setdiff(extra,names(d0)),collapse=", "))
signature<-list(version=2L,prot_md5=unname(tools::md5sum(pf)),
  cohort_md5=unname(tools::md5sum(file.path(rawdir,"bmonth.cohort.rds"))),
  code=unname(tools::md5sum(code_files)),exposures=EXPLORATORY_TRAITS,extra=extra,
  packages=vapply(c("survival","dplyr"),function(p)as.character(packageVersion(p)),character(1)))
cache<-file.path(rawdir,"protwas_followup_cache_v2",med_hash(signature))
med_atomic(signature,file.path(cache,"signature.rds"))

# Fixed from the 2026-10-20 audit, hence an exploratory score, not an independent
# replication. Require all seven proteins; do not change membership to fit data.
markers<-c("GZMH","FCRL6","KIR2DL2","KLRD1","KIR2DL3","GZMA","GNLY")
score_name<-"cytotoxic_immune_score"
score_ok<-all(markers%in%features)
score_note<-if(score_ok)"seven-protein membership fixed; standardised proteins averaged" else
  paste("not run; missing markers:",paste(setdiff(markers,features),collapse=", "))
if(score_ok){
  finite_marker<-vapply(d0[markers],function(z)sum(is.finite(z))>=1000 && sd(z[is.finite(z)])>0,logical(1))
  score_ok<-all(finite_marker)
  if(score_ok){
    scaled<-lapply(d0[markers],function(z){z[!is.finite(z)]<-NA_real_;as.numeric(scale(z))})
    score_mat<-as.matrix(as.data.frame(scaled))
    d0[[score_name]]<-rowMeans(score_mat,na.rm=FALSE)
    d0[[score_name]]<-as.numeric(scale(d0[[score_name]]))
  } else score_note<-"not run; at least one fixed marker has insufficient variation"
}
run_features<-c(features,if(score_ok)score_name)
all_rows<-list()
completion<-data.frame(exposure=EXPLORATORY_TRAITS,status="pending",stringsAsFactors=FALSE)
write_csv2(completion,"protwas_followup_completion.csv")
write_csv2(data.frame(code_version=BMONTH_CODE_VERSION,proteins=length(features),
  shared_N=nrow(d0),score_status=if(score_ok)"available"else"unavailable",score_note=score_note,
  extra_covariates=paste(extra,collapse=","),
  technical_adjustment_note=if(length(extra))"Only explicitly listed fields adjusted; verify biological and technical meaning"else
    "Batch, sampling season/time, cell counts, LOD and ancestry not added by this default run",
  selection_status="Exploratory hypotheses selected using the 2026-10-20 results"),
  "protwas_followup_manifest.csv")

pair_models<-function(x,m){
  cv0<-analysis_covars_for(x)
  cv1<-unique(c(cv0,if(x%in%c("border1n","border_rank","border_relative"))"sibship_size_cat"))
  needed<-unique(c("eid",x,m,cv1))
  ok<-complete.cases(d0[,needed,drop=FALSE])
  for(v in needed)if(is.numeric(d0[[v]]))ok<-ok&is.finite(d0[[v]])
  z<-droplevels(d0[ok,needed,drop=FALSE])
  key<-bm_hash(sort(as.character(z$eid)))
  a<-screen_a(z,x,m,covar_names=cv0);a$profile<-"D2_common_family_sample"
  b<-screen_a(z,x,m,covar_names=cv1);b$profile<-"D2_plus_family_size"
  r<-bind_rows(a,b);r$sample_key<-key
  if(length(extra)){
    n2<-unique(c(needed,extra));ok2<-complete.cases(d0[,n2,drop=FALSE])
    for(v in n2)if(is.numeric(d0[[v]]))ok2<-ok2&is.finite(d0[[v]])
    z2<-droplevels(d0[ok2,n2,drop=FALSE]);key2<-bm_hash(sort(as.character(z2$eid)))
    c0<-screen_a(z2,x,m,covar_names=cv1);c0$profile<-"family_on_extra_common_sample";c0$sample_key<-key2
    c1<-screen_a(z2,x,m,covar_names=unique(c(cv1,extra)));c1$profile<-"family_plus_explicit_extra";c1$sample_key<-key2
    r<-bind_rows(r,c0,c1)
  }
  r$exposure<-x
  r$added_family_covariate<-x!="border12" & r$profile!="D2_common_family_sample"
  r$family_design<-if(x=="border12")"fixed_sibship_size_2" else ifelse(r$added_family_covariate,"categorical_sibship_adjustment","unadjusted_for_family_size")
  r$exposure_unit<-if(x%in%c("border12","border1n"))"0 to 1 contrast"else"per SD on paired sample"
  r$family_id<-paste(x,r$profile,if(m==score_name)"locked_exploratory_score"else"proteome_screen",sep="/")
  r
}
for(x in EXPLORATORY_TRAITS){
  message("[omics-followup] START ",x,"; ",length(run_features)," proteins/score")
  completion$status[completion$exposure==x]<-"running";write_csv2(completion,"protwas_followup_completion.csv")
  z<-med_cached(cache,list("paired_protwas",x),function()
    bind_rows(lapply(run_features,function(m)pair_models(x,m))),REPLACE)$value
  z$FDR<-NA_real_
  for(ix in split(seq_len(nrow(z)),z$family_id))z$FDR[ix]<-bm_planned_adjust(z$p_a[ix],z$eligible[ix])
  z$selection_status<-"exploratory; score and targets motivated by the audited results"
  all_rows[[x]]<-z
  results<-bind_rows(all_rows)
  write_csv2(results,"protwas_family_size.csv")
  med_atomic(results,file.path(rawdir,"protwas_family_size.rds"))
  write_csv2(results[results$mediator%in%c(markers,score_name),,drop=FALSE],"immune_score_results.csv")
  completion$status[completion$exposure==x]<-"complete";write_csv2(completion,"protwas_followup_completion.csv")
}
r<-bind_rows(all_rows)
if(nrow(r)){
  a<-r%>%filter(profile=="D2_common_family_sample")%>%
    select(exposure,mediator,N_base=N,a_base=a,p_base=p_a,FDR_base=FDR,sample_base=sample_key)
  b<-r%>%filter(profile=="D2_plus_family_size")%>%
    select(exposure,mediator,N_family=N,a_family=a,p_family=p_a,FDR_family=FDR,sample_family=sample_key)
  cmp<-left_join(a,b,by=c("exposure","mediator"))%>%
    mutate(same_sample=sample_base==sample_family,delta_beta=a_family-a_base,
      attenuation_beta=if_else(is.finite(a_base)&abs(a_base)>1e-8,1-a_family/a_base,NA_real_),
      interpretation="paired descriptive coefficient change; not a causal mediated fraction")
  if(any(!cmp$same_sample))stop("Internal error: family-size comparison samples differ")
  write_csv2(cmp,"protwas_model_comparison.csv")
}
message("[omics-followup] complete; inspect immune_score_results.csv and protwas_model_comparison.csv")
