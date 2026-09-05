# Cross-fitted PGS calibration on baseline disease-free training participants.
# No future event-free threshold is used for calibration.
le8_oof_decompose <- function(d,feature,gcol,covars,k=5,seed=2026) {
  d<-as.data.frame(d);d$.omic<-as.numeric(d[[feature]]);d$.grs<-as.numeric(d[[gcol]])
  set.seed(seed);ord<-order(as.character(d$eid));fold<-integer(nrow(d))
  fold[ord]<-sample(rep(seq_len(k),length.out=nrow(d)))
  for(nm in c(".omic_z",".genetic_z",".residual_z",".full_prediction",".reduced_prediction"))d[[nm]]<-NA_real_
  reports<-list()
  for(f in seq_len(k)) {
    train<-d[fold!=f&!is.na(d$.prevalent)&d$.prevalent==0&complete.cases(d[,c(".omic",".grs",covars),drop=FALSE]),,drop=FALSE]
    test<-which(fold==f&complete.cases(d[,c(".omic",".grs",covars),drop=FALSE]))
    if(nrow(train)<200||!length(test)||sd(train$.omic)<=0||sd(train$.grs)<=0)next
    om<-mean(train$.omic);os<-sd(train$.omic);gm<-mean(train$.grs);gs<-sd(train$.grs)
    train$.y<-(train$.omic-om)/os;train$.g<-(train$.grs-gm)/gs
    full<-tryCatch(lm(reformulate(c(".g",covars),".y"),train),error=function(e)NULL)
    reduced<-tryCatch(lm(reformulate(covars,".y"),train),error=function(e)NULL)
    if(is.null(full)||is.null(reduced)||!is.finite(coef(full)[".g"]))next
    te<-d[test,,drop=FALSE];te$.g<-(te$.grs-gm)/gs
    bg<-unname(coef(full)[".g"])
    d$.omic_z[test]<-(te$.omic-om)/os
    d$.genetic_z[test]<-bg*te$.g
    d$.residual_z[test]<-d$.omic_z[test]-d$.genetic_z[test]
    d$.full_prediction[test]<-tryCatch(as.numeric(predict(full,te)),error=function(e)rep(NA_real_,length(test)))
    d$.reduced_prediction[test]<-tryCatch(as.numeric(predict(reduced,te)),error=function(e)rep(NA_real_,length(test)))
    reports[[length(reports)+1L]]<-tibble(feature,fold=f,N_train=nrow(train),N_test=length(test),
      genetic_beta=bg,omic_mean=om,omic_sd=os,pgs_mean=gm,pgs_sd=gs,
      calibration_rule="baseline disease-free training fold; future outcomes not consulted")
  }
  d$.calibration_fold<-factor(fold)
  ii<-!is.na(d$.prevalent)&d$.prevalent==0&is.finite(d$.full_prediction)&is.finite(d$.reduced_prediction)&is.finite(d$.omic_z)
  mse0<-sum((d$.omic_z[ii]-d$.reduced_prediction[ii])^2)
  r2<-if(mse0>0)1-sum((d$.omic_z[ii]-d$.full_prediction[ii])^2)/mse0 else NA_real_
  list(data=d,folds=bind_rows(reports),partial_R2=r2,N_calibration=sum(ii))
}
risk_set_component_scan <- function(dd,feature,covars,tvar,evar,cuts=c(0,.5,1,2,5,10,16)) {
  components<-c(Observed=".omic_z",`PGS-predicted`=".genetic_z",Residual=".residual_z")
  map_dfr(seq_len(length(cuts)-1L),function(i)imap_dfr(components,function(x,label){
    z<-le8_interval_cox(dd,x,tvar,evar,covars,cuts[i],cuts[i+1],scale_x=FALSE)
    z|>transmute(feature,component=label,lead_lo=window_lo,lead_hi=window_hi,lead_mid=time,
      beta,se=std.error,conf.low,conf.high,p=p.value,N_case,N_control,N_risk=N_total,
      person_years,status,effect_measure,unit="common training-fold omic scale")
  }))
}
run_individual_genetic_decomposition <- function(layer,outdir,top_candidates) {
  sf<-find_c2_score_file(layer);empty<-list(status=tibble(status="PGS unavailable"),summary=tibble(),trajectory=tibble(),folds=tibble())
  if(length(sf)!=1L||is.na(sf))return(empty)
  le8_dependency(sf);scores<-read_c2_scores(sf);biom<-if(layer=="protein")read_prot()else read_met()
  features<-unique(c(C2_FIXED_TOP,as.character(top_candidates$feature)))
  features<-head(intersect(features,names(biom)),le8_num_env("C2_DECOMP_MAX",100))
  mp<-map_c2_score_columns(features,names(scores));if(!length(mp))return(empty)
  covars<-unique(c(vars.basic,le8_csv_env("C2_LE4_COVARS","diet.pts,pa.pts,smoke.pts,sleep.pts"),
    le8_csv_env("C2_TREATMENT_VARS")))
  ph<-read_all(unique(c("eid","ethnic.c",covars,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y))))|>
    filter_analysis_cohort()|>make_outcome(Y)
  ph$.prevalent<-make_prevalent_status(ph,Y)
  for(nm in c("ph","scores","biom")){z<-get(nm);z$eid<-as.character(z$eid);assign(nm,z)}
  d<-inner_join(ph,biom[,c("eid",names(mp)),drop=FALSE],by="eid")|>
    inner_join(scores[,unique(c("eid",unname(mp))),drop=FALSE],by="eid")
  rm(ph,scores,biom);invisible(gc());covars<-intersect(covars,names(d))
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e")
  results<-lapply(names(mp),function(feature){
    z<-le8_oof_decompose(d[,unique(c("eid",feature,mp[[feature]],covars,tvar,evar,".prevalent")),drop=FALSE],
      feature,mp[[feature]],covars,k=as.integer(le8_num_env("C2_DECOMP_FOLDS",5)),seed=SEED)
    dd<-z$data;cv<-c(covars,".calibration_fold")
    co<-component_association(dd,".omic_z",cv,tvar,evar);cg<-component_association(dd,".genetic_z",cv,tvar,evar)
    cr<-component_association(dd,".residual_z",cv,tvar,evar)
    po<-component_association(dd,".omic_z",cv,tvar,evar,TRUE);pg<-component_association(dd,".genetic_z",cv,tvar,evar,TRUE)
    pr<-component_association(dd,".residual_z",cv,tvar,evar,TRUE)
    ji<-joint_component_association(dd,cv,tvar,evar);jp<-joint_component_association(dd,cv,tvar,evar,TRUE)
    take<-function(x,term,col){i<-match(term,x$component);if(is.na(i))NA_real_ else x[[col]][i]}
    out<-tibble(feature,score_column=mp[[feature]],status=if(any(is.finite(dd$.genetic_z)))"ok"else"calibration unavailable",
      N_calibration=z$N_calibration,N_incident=co[["N"]],incident_events=co[["events"]],
      genetic_beta=mean(z$folds$genetic_beta),genetic_partial_R2=z$partial_R2,pgs_partial_R2=z$partial_R2,
      incident_beta_observed=co[["beta"]],incident_se_observed=co[["se"]],incident_p_observed=co[["p"]],
      incident_beta_genetic=cg[["beta"]],incident_se_genetic=cg[["se"]],incident_p_genetic=cg[["p"]],
      incident_beta_residual=cr[["beta"]],incident_se_residual=cr[["se"]],incident_p_residual=cr[["p"]],
      incident_joint_beta_genetic=take(ji,".genetic_z","beta"),incident_joint_p_genetic=take(ji,".genetic_z","p"),
      incident_joint_beta_residual=take(ji,".residual_z","beta"),incident_joint_p_residual=take(ji,".residual_z","p"),
      prevalent_beta_observed=po[["beta"]],prevalent_p_observed=po[["p"]],
      prevalent_beta_genetic=pg[["beta"]],prevalent_p_genetic=pg[["p"]],
      prevalent_beta_residual=pr[["beta"]],prevalent_p_residual=pr[["p"]],
      prevalent_joint_beta_genetic=take(jp,".genetic_z","beta"),prevalent_joint_p_genetic=take(jp,".genetic_z","p"),
      prevalent_joint_beta_residual=take(jp,".residual_z","beta"),prevalent_joint_p_residual=take(jp,".residual_z","p"),
      cross_fitted=TRUE,pgs_scope="existing external/overlapping COJO score; pQTL GWAS overlap must be audited separately",
      causal_fraction_identifiable=FALSE,incident_genetic_signal_fraction_abs=NA_real_,
      interpretation="Residual includes uncaptured genetics, environment, treatment, assay error and disease processes; not a disease-consequence fraction",
      unit="common training-fold omic scale; components not separately standardized")
    trajectory<-if(feature%in%head(names(mp),le8_num_env("C2_DECOMP_TRAJECTORY_MAX",20)))
      risk_set_component_scan(dd,feature,cv,tvar,evar)else tibble()
    list(summary=out,folds=z$folds,trajectory=trajectory)
  })
  summary<-bind_rows(lapply(results,`[[`,"summary"));folds<-bind_rows(lapply(results,`[[`,"folds"))
  trajectory<-bind_rows(lapply(results,`[[`,"trajectory"));h2<-read_c2_heritability(layer)
  summary<-left_join(summary,h2$data,by="feature")
  summary$pgs_h2_coverage<-NA_real_ # partial R2 and external h2 are not a captured-causal-fraction estimator.
  for(p in grep("^(incident|prevalent)(_joint)?_p_",names(summary),value=TRUE))summary[[paste0("FDR_",p)]]<-p.adjust(summary[[p]],"BH")
  rd<-le8_job_dir(outdir,"c2_cause")
  write_raw_csv(summary,"c2.individual_genetic_decomposition.csv",rd)
  write_raw_csv(trajectory,"c2.genetic_component_leadtime.csv",rd)
  write_raw_csv(folds,"c2.decomposition_calibration_folds.csv",rd)
  list(status=tibble(status="ok",cross_fitted=TRUE,calibration="baseline disease-free training folds only"),
    summary=summary,trajectory=trajectory,folds=folds,heritability_status=h2$status)
}
plot_individual_genetic_decomposition <- function(summary) {
  z<-as_tibble(summary);if(!nrow(z)||!"pgs_partial_R2"%in%names(z))return(blank_plot("PGS decomposition","No calibrated scores"))
  z<-z|>filter(status=="ok")|>slice_head(n=30)
  a<-ggplot(z,aes(100*pgs_partial_R2,reorder(feature,pgs_partial_R2)))+geom_col()+
    labs(title="a. Out-of-fold partial R² (negative values retained)",x="PGS incremental explained variance (%)",y=NULL)+theme_5c(9)
  eff<-bind_rows(z|>transmute(feature,component="Observed",beta=incident_beta_observed,se=incident_se_observed),
    z|>transmute(feature,component="PGS-predicted",beta=incident_beta_genetic,se=incident_se_genetic),
    z|>transmute(feature,component="Residual",beta=incident_beta_residual,se=incident_se_residual))
  b<-ggplot(eff,aes(beta,feature,color=component))+geom_vline(xintercept=0)+
    geom_errorbarh(aes(xmin=beta-1.96*se,xmax=beta+1.96*se),height=.1,position=position_dodge(width=.5))+
    geom_point(position=position_dodge(width=.5))+labs(title="b. Associations on a common omic scale",x="Log HR per training-fold omic SD",y=NULL,color=NULL)+theme_5c(9)
  c<-ggplot(z,aes(N_incident,feature))+geom_point()+labs(title="c. Disease-model sample sizes (not calibration N)",x="Incident-model N",y=NULL)+theme_5c(9)
  (a|b)/c+plot_annotation(caption="No causal/consequence percentage is calculated. Cross-fitting calibration does not remove overlap in the original pQTL GWAS.")
}
