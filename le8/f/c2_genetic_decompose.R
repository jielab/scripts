# C2 individual-level decomposition of an observed omic trait into a simple
# COJO-weighted PGS component and a non-genetic residual.
#
# Auto-discovered inputs:
#   <UKB_PHE>/Rdata/prot.pgs.rds  (eid, FEATURE.pgs, ...)
#   <UKB_PHE>/Rdata/met.pgs.rds   (eid, FEATURE.pgs, ...)
# C2_GENETIC_SCORE_FILE is an optional explicit override. This is an
# exploratory in-sample decomposition, not external/cross-fitted prediction.

find_c2_score_file <- function(layer){
  explicit<-Sys.getenv("C2_GENETIC_SCORE_FILE",unset="")
  automatic<-file.path(indir,"Rdata",if(layer=="protein")"prot.pgs.rds"else"met.pgs.rds")
  z<-unique(c(explicit,automatic));z<-z[nzchar(z)&file.exists(z)&file.size(z)>0]
  if(length(z))normalizePath(z[[1]],winslash="/",mustWork=FALSE)else NA_character_
}

read_c2_scores <- function(file){
  x<-if(grepl("\\.rds$",file,ignore.case=TRUE))readRDS(file)else
    data.table::fread(file,showProgress=FALSE,check.names=FALSE)
  x<-as_tibble(x)
  if(!"eid"%in%names(x))stop("PGS file must contain eid: ",file,call.=FALSE)
  x
}

map_c2_score_columns <- function(features,nms){
  one<-function(f){
    z<-c(paste0(f,".pgs"),paste0(f,"_pgs"),paste0(f,".PGS"),
      paste0(f,"_PGS"),paste0(f,"_GRS"),paste0("GRS_",f),f)
    hit<-z[z%in%nms];if(length(hit))hit[[1]]else NA_character_
  }
  x<-setNames(vapply(features,one,character(1)),features);x[!is.na(x)]
}

component_association <- function(dd,x,covars,tvar,evar,prevalent=FALSE){
  covars<-intersect(covars,names(dd))
  if(prevalent){
    need<-unique(c(".prevalent",x,covars));z<-dd[,need,drop=FALSE]
    z<-z[complete.cases(z),,drop=FALSE];events<-sum(z$.prevalent==1)
    if(nrow(z)<500||events<20||length(unique(z$.prevalent))<2)
      return(c(beta=NA,se=NA,p=NA,N=nrow(z),events=events))
    fit<-tryCatch(glm(reformulate(c(x,covars),".prevalent"),z,family=binomial()),error=function(e)NULL)
  }else{
    need<-unique(c(tvar,evar,x,covars));z<-dd[,need,drop=FALSE]
    z<-z[complete.cases(z),,drop=FALSE]
    z<-z[is.finite(z[[tvar]])&z[[tvar]]>0&z[[evar]]%in%c(0,1),,drop=FALSE]
    events<-sum(z[[evar]]==1)
    if(nrow(z)<500||events<20)return(c(beta=NA,se=NA,p=NA,N=nrow(z),events=events))
    ff<-as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",
      paste(bt(c(x,covars)),collapse=" + ")))
    fit<-tryCatch(coxph(ff,z,ties="efron"),error=function(e)NULL)
  }
  if(is.null(fit))return(c(beta=NA,se=NA,p=NA,N=nrow(z),events=events))
  sm<-coef(summary(fit));if(!x%in%rownames(sm))
    return(c(beta=NA,se=NA,p=NA,N=nrow(z),events=events))
  if(prevalent)c(beta=sm[x,"Estimate"],se=sm[x,"Std. Error"],p=sm[x,"Pr(>|z|)"],
    N=nrow(z),events=events) else
    c(beta=sm[x,"coef"],se=sm[x,"se(coef)"],p=sm[x,"Pr(>|z|)"],
      N=nrow(z),events=events)
}

joint_component_association <- function(dd,covars,tvar,evar,prevalent=FALSE){
  covars<-intersect(covars,names(dd));xs<-c(".genetic_z",".residual_z")
  if(prevalent){
    need<-unique(c(".prevalent",xs,covars));z<-dd[,need,drop=FALSE]
    z<-z[complete.cases(z),,drop=FALSE]
    fit<-if(nrow(z)>=500&&sum(z$.prevalent==1)>=20)
      tryCatch(glm(reformulate(c(xs,covars),".prevalent"),z,family=binomial()),error=function(e)NULL) else NULL
  }else{
    need<-unique(c(tvar,evar,xs,covars));z<-dd[,need,drop=FALSE]
    z<-z[complete.cases(z),,drop=FALSE]
    z<-z[is.finite(z[[tvar]])&z[[tvar]]>0&z[[evar]]%in%c(0,1),,drop=FALSE]
    ff<-as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",
      paste(bt(c(xs,covars)),collapse=" + ")))
    fit<-if(nrow(z)>=500&&sum(z[[evar]]==1)>=20)
      tryCatch(coxph(ff,z,ties="efron"),error=function(e)NULL) else NULL
  }
  if(is.null(fit))return(tibble(component=xs,beta=NA_real_,se=NA_real_,p=NA_real_))
  sm<-coef(summary(fit))
  tibble(component=xs,
    beta=vapply(xs,function(x)if(x%in%rownames(sm))sm[x,if(prevalent)"Estimate"else"coef"]else NA_real_,numeric(1)),
    se=vapply(xs,function(x)if(x%in%rownames(sm))sm[x,if(prevalent)"Std. Error"else"se(coef)"]else NA_real_,numeric(1)),
    p=vapply(xs,function(x)if(x%in%rownames(sm))sm[x,"Pr(>|z|)"]else NA_real_,numeric(1)))
}

risk_set_component_scan <- function(dd,feature,covars,tvar,evar,
  cuts=c(0,.5,1,2,5,10,16)){
  covars<-intersect(covars,names(dd))
  components<-c(Observed=".omic_z",PGS=".genetic_z",Residual=".residual_z")
  map_dfr(seq_len(length(cuts)-1L),function(i){
    lo<-cuts[[i]];hi<-cuts[[i+1L]]
    # Later cases are valid controls for an earlier window. Participants
    # censored inside the window are excluded because disease-free status at
    # the upper boundary is unknown.
    is_case<-dd[[evar]]==1&is.finite(dd[[tvar]])&dd[[tvar]]>=lo&dd[[tvar]]<hi
    is_control<-is.finite(dd[[tvar]])&dd[[tvar]]>=hi
    keep<-is_case|is_control;base<-dd[keep,,drop=FALSE]
    base$.window_case<-as.integer(is_case[keep])
    map_dfr(names(components),function(label){
      x<-components[[label]];need<-unique(c(".window_case",x,covars))
      z<-base[,need,drop=FALSE];z<-z[complete.cases(z),,drop=FALSE]
      nc<-sum(z$.window_case==1);nn<-sum(z$.window_case==0)
      empty<-tibble(feature,component=label,lead_lo=lo,lead_hi=hi,
        lead_mid=(lo+hi)/2,beta=NA_real_,se=NA_real_,conf.low=NA_real_,
        conf.high=NA_real_,p=NA_real_,N_case=nc,N_control=nn)
      if(nc<20||nn<100||!is.finite(sd(z[[x]]))||sd(z[[x]])<=0)return(empty)
      fit<-tryCatch(glm(reformulate(c(x,covars),".window_case"),z,family=binomial()),error=function(e)NULL)
      if(is.null(fit))return(empty);sm<-coef(summary(fit))
      if(!x%in%rownames(sm))return(empty)
      b<-sm[x,"Estimate"];se<-sm[x,"Std. Error"]
      tibble(feature,component=label,lead_lo=lo,lead_hi=hi,
        lead_mid=(lo+hi)/2,beta=b,se=se,conf.low=b-1.96*se,
        conf.high=b+1.96*se,p=sm[x,"Pr(>|z|)"],N_case=nc,N_control=nn)
    })
  })
}

plot_individual_genetic_decomposition <- function(summary){
  if(is.null(summary)||!is.data.frame(summary)||!nrow(summary)||
     !all(c("status","genetic_partial_R2")%in%names(summary)))
    return(blank_plot("PGS calibration and disease association",
      "PGS input unavailable; individual genetic decomposition was skipped"))
  d<-as_tibble(summary)|>filter(status=="ok",is.finite(genetic_partial_R2))|>
    arrange(desc(genetic_partial_R2))|>slice_head(n=30)|>
    mutate(feature=factor(feature,levels=rev(feature)),r2=100*genetic_partial_R2)
  if(!nrow(d))return(blank_plot("PGS calibration and disease association",
    "No matched prot.pgs.rds/met.pgs.rds score could be calibrated"))
  pa<-ggplot(d,aes(r2,feature,fill=genetic_beta>=0))+geom_col(width=.72)+
    scale_fill_manual(values=c(`TRUE`="#3F78A8",`FALSE`="#D95F02"),guide="none")+
    labs(title="a. Observed omic variance explained by its COJO PGS",
      subtitle="Partial R² calibrated among participants event-free and observed for at least 10 years",
      x="Partial R² (%)",y=NULL)+theme_5c(9)
  eff<-bind_rows(
    d|>transmute(feature,component="Observed",beta=incident_beta_observed,se=incident_se_observed),
    d|>transmute(feature,component="PGS",beta=incident_beta_genetic,se=incident_se_genetic),
    d|>transmute(feature,component="Residual",beta=incident_beta_residual,se=incident_se_residual))|>
    filter(is.finite(beta),is.finite(se))|>mutate(lo=beta-1.96*se,hi=beta+1.96*se)
  keep<-as.character(tail(levels(d$feature),12));eff<-eff|>filter(as.character(feature)%in%keep)
  pb<-if(!nrow(eff))blank_plot("b. Component-specific incident associations","No estimable Cox model")else
    ggplot(eff,aes(beta,feature,color=component))+geom_vline(xintercept=0,color="grey75")+
      geom_errorbarh(aes(xmin=lo,xmax=hi),height=.12,position=position_dodge(width=.55))+
      geom_point(position=position_dodge(width=.55),size=2)+
      labs(title="b. Separate incident Cox models",
        subtitle="Residual includes environment, treatment and assay noise; it is not pure preclinical disease",
        x="Log HR per 1-SD component",y=NULL,color=NULL)+theme_5c(9)+theme(legend.position="top")
  pa|pb
}

plot_component_leadtime <- function(x){
  if(is.null(x)||!is.data.frame(x)||!nrow(x)||!all(c("beta","se")%in%names(x)))
    return(blank_plot("PGS and residual lead-time associations",
      "PGS input unavailable; individual genetic decomposition was skipped"))
  d<-as_tibble(x)|>filter(is.finite(beta),is.finite(se))
  if(!nrow(d))return(blank_plot("PGS and residual lead-time associations",
    "No risk-set window estimate was available"))
  ord<-d|>group_by(feature)|>summarise(best_p=safe_min_finite(p),.groups="drop")|>
    arrange(best_p)|>slice_head(n=12)|>pull(feature)
  d<-d|>filter(feature%in%ord)|>mutate(feature=factor(feature,levels=rev(ord)))
  ggplot(d,aes(lead_mid,beta,color=component,group=component))+
    geom_hline(yintercept=0,color="grey72")+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high,fill=component),alpha=.10,color=NA)+
    geom_line(linewidth=.75)+geom_point(aes(size=N_case),alpha=.88)+
    facet_wrap(~feature,scales="free_y",ncol=3)+
    scale_x_continuous(breaks=c(.5,1,2,5,10,16),limits=c(0,16))+
    labs(title="Risk-set lead-time associations of observed, PGS and residual components",
      subtitle="Cases are diagnosed in each window; controls remain observed and disease-free through its upper boundary",
      x="Years from baseline to diagnosis window",y="Log odds ratio per 1-SD",
      color=NULL,fill=NULL,size="Cases")+theme_5c(8)+theme(legend.position="bottom")
}

run_individual_genetic_decomposition <- function(layer,outdir,candidates=tibble(),
  max_features=as.integer(Sys.getenv("C2_DECOMP_MAX_FEATURES",unset="100"))){
  score_file<-find_c2_score_file(layer)
  expected<-file.path(indir,"Rdata",if(layer=="protein")"prot.pgs.rds"else"met.pgs.rds")
  empty<-list(status=tibble(status="not run",detail=paste0("Score file not found: ",expected)),
    summary=tibble(),trajectory=tibble())
  if(is.na(score_file))return(empty)
  scores<-read_c2_scores(score_file)
  biom<-if(layer=="protein")read_prot()else read_met();features<-setdiff(names(biom),"eid")
  score_cols<-map_c2_score_columns(features,names(scores))
  if(!length(score_cols))return(modifyList(empty,list(status=tibble(status="not run",
    detail="No FEATURE.pgs column matched an assayed feature"))))
  preferred<-unique(c(as.character(candidates$feature%||%character()),names(score_cols)))
  preferred<-preferred[preferred%in%names(score_cols)]
  if(!is.finite(max_features)||max_features<1L)max_features<-length(preferred)
  score_cols<-score_cols[head(preferred,max_features)]
  ph<-read_all(unique(c("eid","ethnic.c",vars.adj2,"birth_date","date_attend",
    "date_lost","date_death",paste0("fod_icd10_",Y))))|>
    filter_analysis_cohort()|>make_outcome(Y)
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");ydate<-paste0("fod_icd10_",Y)
  yd<-as.Date(ph[[ydate]]);ba<-as.Date(ph$date_attend)
  ph$.prevalent<-case_when(!is.na(yd)&!is.na(ba)&yd<ba~1,
    is.na(yd)|(!is.na(yd)&!is.na(ba)&yd>ba)~0,TRUE~NA_real_)
  ph$eid<-as.character(ph$eid);biom$eid<-as.character(biom$eid);scores$eid<-as.character(scores$eid)
  d<-ph|>inner_join(biom[,c("eid",names(score_cols)),drop=FALSE],by="eid")|>
    inner_join(scores[,unique(c("eid",unname(score_cols))),drop=FALSE],by="eid")
  covars<-intersect(vars.adj2,names(d))
  # Calibration uses everyone known to be event-free at year 10, including
  # participants diagnosed later. Restricting to lifelong non-cases would
  # condition the omic~PGS calibration on the eventual outcome.
  distal_controls<-is.finite(d[[tvar]])&d[[tvar]]>=10
  trajectory_max<-suppressWarnings(as.integer(Sys.getenv("C2_DECOMP_TRAJECTORY_MAX",unset="20")))
  if(!is.finite(trajectory_max)||trajectory_max<1L)trajectory_max<-20L
  trajectory_features<-head(names(score_cols),trajectory_max)
  rows<-parallel_map(names(score_cols),function(feature){
    gcol<-score_cols[[feature]];need<-unique(c(feature,gcol,covars,tvar,evar,".prevalent"))
    dd<-d[,need,drop=FALSE]
    dd$.omic<-suppressWarnings(as.numeric(dd[[feature]]))
    dd$.grs<-suppressWarnings(as.numeric(dd[[gcol]]))
    train<-dd[distal_controls&complete.cases(dd[,c(".omic",".grs",covars),drop=FALSE]),,drop=FALSE]
    failed<-function(status)list(summary=tibble(feature,score_column=gcol,status),trajectory=tibble())
    if(nrow(train)<500||!is.finite(sd(train$.omic))||sd(train$.omic)<=0||
       !is.finite(sd(train$.grs))||sd(train$.grs)<=0)return(failed("insufficient distal controls"))
    om_m<-mean(train$.omic);om_s<-sd(train$.omic);g_m<-mean(train$.grs);g_s<-sd(train$.grs)
    dd$.omic_z<-(dd$.omic-om_m)/om_s;dd$.grs_z<-(dd$.grs-g_m)/g_s
    tr<-train|>mutate(.omic_z=(.omic-om_m)/om_s,.grs_z=(.grs-g_m)/g_s)
    fit<-tryCatch(lm(reformulate(c(".grs_z",covars),".omic_z"),tr),error=function(e)NULL)
    reduced<-tryCatch(lm(reformulate(covars,".omic_z"),tr),error=function(e)NULL)
    if(is.null(fit)||is.null(reduced)||!".grs_z"%in%names(coef(fit)))return(failed("calibration failed"))
    bg<-unname(coef(fit)[[".grs_z"]]);dd$.genetic<-bg*dd$.grs_z
    dd$.residual<-dd$.omic_z-dd$.genetic
    dd$.genetic_z<-as.numeric(scale(dd$.genetic));dd$.residual_z<-as.numeric(scale(dd$.residual))
    co<-component_association(dd,".omic_z",covars,tvar,evar)
    cg<-component_association(dd,".genetic_z",covars,tvar,evar)
    cr<-component_association(dd,".residual_z",covars,tvar,evar)
    po<-component_association(dd,".omic_z",covars,tvar,evar,TRUE)
    pg<-component_association(dd,".genetic_z",covars,tvar,evar,TRUE)
    pr<-component_association(dd,".residual_z",covars,tvar,evar,TRUE)
    ji<-joint_component_association(dd,covars,tvar,evar,FALSE)
    jp<-joint_component_association(dd,covars,tvar,evar,TRUE)
    pick<-function(z,component,column){i<-match(component,z$component);if(is.na(i))NA_real_ else z[[column]][[i]]}
    denom<-deviance(reduced);r2<-if(is.finite(denom)&&denom>0)
      max(0,min(1,(denom-deviance(fit))/denom)) else NA_real_
    smry<-tibble(feature,score_column=gcol,status="ok",N_calibration=nrow(train),
      genetic_beta=bg,genetic_partial_R2=r2,
      incident_beta_observed=co[["beta"]],incident_se_observed=co[["se"]],incident_p_observed=co[["p"]],
      incident_beta_genetic=cg[["beta"]],incident_se_genetic=cg[["se"]],incident_p_genetic=cg[["p"]],
      incident_beta_residual=cr[["beta"]],incident_se_residual=cr[["se"]],incident_p_residual=cr[["p"]],
      incident_joint_beta_genetic=pick(ji,".genetic_z","beta"),incident_joint_p_genetic=pick(ji,".genetic_z","p"),
      incident_joint_beta_residual=pick(ji,".residual_z","beta"),incident_joint_p_residual=pick(ji,".residual_z","p"),
      prevalent_beta_observed=po[["beta"]],prevalent_p_observed=po[["p"]],
      prevalent_beta_genetic=pg[["beta"]],prevalent_p_genetic=pg[["p"]],
      prevalent_beta_residual=pr[["beta"]],prevalent_p_residual=pr[["p"]],
      prevalent_joint_beta_genetic=pick(jp,".genetic_z","beta"),prevalent_joint_p_genetic=pick(jp,".genetic_z","p"),
      prevalent_joint_beta_residual=pick(jp,".residual_z","beta"),prevalent_joint_p_residual=pick(jp,".residual_z","p"),
      pgs_scope=Sys.getenv("C2_PGS_SCOPE",unset="COJO lead SNPs; scope supplied by score builder"),
      cross_fitted=FALSE,
      interpretation="Residual contains environment, treatment, assay error and latent disease; it is not a pure preclinical-disease component")
    trj<-if(feature%in%trajectory_features)
      risk_set_component_scan(dd,feature,covars,tvar,evar) else tibble()
    list(summary=smry,trajectory=trj)
  })
  summary<-bind_rows(lapply(rows,`[[`,"summary"))
  trajectory<-bind_rows(lapply(rows,`[[`,"trajectory"))
  if(nrow(summary)){
    pcols<-grep("^(incident|prevalent)(_joint)?_p_",names(summary),value=TRUE)
    for(nm in pcols)summary[[paste0("FDR_",nm)]]<-p.adjust(summary[[nm]],"BH")
  }
  rawdir<-le8_job_dir(outdir,"c2_cause")
  write_raw_csv(summary,"c2.individual_genetic_decomposition.csv",rawdir)
  write_raw_csv(trajectory,"c2.genetic_component_leadtime.csv",rawdir)
  list(status=tibble(status="ok",detail=paste(length(score_cols),"matched scores"),
    score_file=score_file,score_construction="COJO-weighted PGS; exploratory, not cross-fitted"),
    summary=summary,trajectory=trajectory)
}
