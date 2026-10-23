source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
dat<-if(length(OUTCOMES))load_cohort() else NULL

res<-bind_rows(empty_associations(),lapply(OUTCOMES,function(y)
  bind_rows(lapply(ANALYSIS_EXPOSURES,function(x)fit_cox(dat,x,y)))))
if(nrow(res)){
  res<-res%>%group_by(exposure)%>%mutate(FDR=bm_planned_adjust(p,eligible))%>%ungroup()
  res$FDR_all_prespecified<-bm_planned_adjust(res$p,res$eligible)
  res$bonferroni_all_prespecified<-.05/max(1,sum(res$eligible))
  res$multiplicity_note<-"FDR within exposure/outcomes; FDR_all_prespecified across all selected scalar exposure-outcome pairs; unresolved eligible fits count as p=1"
}
write_csv2(res,"main_associations.csv")
saveRDS(res,file.path(rawdir,"main_associations.rds"))
print(res,n=Inf)

sens<-bind_rows(res[0,],lapply(OUTCOMES,function(y)
  bind_rows(lapply(FAMILY_SIZE_SENS_TRAITS,function(x)
    fit_cox(dat,x,y,covar_names=exploratory_family_covars,
            model_label="exploratory_D2_plus_family_size")))))
if(nrow(sens)){
  sens<-sens%>%group_by(exposure)%>%mutate(FDR=bm_planned_adjust(p,eligible))%>%ungroup()
  sens$FDR_all_sensitivity<-bm_planned_adjust(sens$p,sens$eligible)
}
write_csv2(sens,"birth_order_family_size_sensitivity.csv")
saveRDS(sens,file.path(rawdir,"birth_order_family_size_sensitivity.rds"))
print(sens,n=Inf)

write_csv2(data.frame(module="assoc",status=if(length(OUTCOMES))"completed"else"skipped_no_outcomes"),"association_status.csv")
