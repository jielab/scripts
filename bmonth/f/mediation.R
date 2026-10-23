# Exploratory linear-by-Cox coefficient products; these are not causal NIEs.
# Family size is included for the three variable-size birth-order traits.
if (!exists("bm_prune_design",mode="function"))
  source(file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R"))

omics_covars_for <- function(x, family=TRUE) {
  out <- analysis_covars_for(x)
  if (family && x %in% c("border1n","border_rank","border_relative"))
    out <- unique(c(out,"sibship_size_cat"))
  unique(c(out,csv_env("BMONTH_OMICS_COVARIATES")))
}
omics_model_label <- function(x) {
  paste0(if(x %in% c("border1n","border_rank","border_relative"))
    "exploratory_D2_plus_family_size" else analysis_model_label(x),
    ";cox_center=",Sys.getenv("BMONTH_MED_CENTER","strata"),
    if(length(csv_env("BMONTH_OMICS_COVARIATES"))) ";plus_explicit_omics_covariates" else "")
}
mediation_data <- function(d,x,m,y=NULL,covar_names=omics_covars_for(x)) {
  vars <- unique(c(x,m,setdiff(covar_names,x),if(!is.null(y))paste0(y,c(".t2e",".Yt2e"))))
  if(!is.null(y)) d<-ensure_outcome(d,y)
  z<-complete_model_data(d,vars)
  if(nrow(z)<1000 || length(unique(z[[x]]))<2 || length(unique(z[[m]]))<2) return(NULL)
  if(!is.null(y) && sum(z[[paste0(y,".Yt2e")]]==1)<50) return(NULL)
  binary<-x %in% c("bmonth_winter_summer","bmonth_spring_autumn","border12","border1n")
  z$.x<-if(binary) z[[x]] else as.numeric(scale(z[[x]]))
  z$.m<-as.numeric(scale(z[[m]]))
  if(!is.null(y)){z$.time<-z[[paste0(y,".t2e")]];z$.event<-z[[paste0(y,".Yt2e")]]}
  attr(z,"med_covariates")<-setdiff(covar_names,x)
  attr(z,"exposure_sd")<-if(binary)1 else sd(z[[x]])
  z
}

# The same nuisance-column policy is used in every fit and every bootstrap.
# Centre-stratified Cox is prespecified, not a significance-based fallback.
med_model <- function(z,type=c("a","b","total")) {
  type<-match.arg(type)
  bad<-function(status,note="")list(ok=FALSE,status=status,note=note,table=NULL,
    eventless_centers=NA_integer_,informative_N=NA_integer_)
  if(is.null(z))return(bad("insufficient_data"))
  cv<-attr(z,"med_covariates"); if(is.null(cv))stop("Missing mediation covariate contract")
  protect<-if(type=="b")c(".x",".m")else".x"
  center_mode<-Sys.getenv("BMONTH_MED_CENTER","strata")
  if(!center_mode %in% c("strata","coefficient"))stop("BMONTH_MED_CENTER must be strata or coefficient")
  stratified<-type!="a" && center_mode=="strata" && "center"%in%cv
  if(stratified)cv<-setdiff(cv,"center")
  matrix_try<-tryCatch(bm_formula_matrix(z,c(protect,cv),tested=protect),error=function(e)e)
  if(inherits(matrix_try,"error"))return(bad("invalid_design",conditionMessage(matrix_try)))
  pr<-bm_prune_design(matrix_try,tested=protect)
  if(!pr$ok)return(bad(pr$status,pr$status))
  xx<-pr$x
  notes<-paste0("Removed nuisance columns: ",paste(unique(c(attr(matrix_try,"dropped_terms"),pr$dropped)),collapse="|"))
  if(type=="a") {
    ff<-model_try(lm.fit(x=xx,y=z$.m))
    if(is.null(ff$fit))return(bad("failed_linear",ff$note))
    f<-ff$fit
    if(f$rank!=ncol(xx) || f$df.residual<=0)return(bad("invalid_linear_rank"))
    inv<-chol2inv(qr.R(f$qr)[seq_len(f$rank),seq_len(f$rank),drop=FALSE])
    v<-matrix(NA_real_,ncol(xx),ncol(xx))
    v[f$qr$pivot,f$qr$pivot]<-inv*sum(f$residuals^2)/f$df.residual
    b<-f$coefficients;s<-sqrt(diag(v));p<-2*pt(abs(b/s),df=f$df.residual,lower.tail=FALSE)
    eventless<-0L; informative<-nrow(z)
  } else {
    xx<-xx[,colnames(xx)!="(Intercept)",drop=FALSE]
    if(stratified) {
      .center<-droplevels(factor(z$center))
      events_by_center<-tapply(z$.event,.center,sum)
      eventless<-sum(events_by_center==0)
      informative<-sum(as.character(.center)%in%names(events_by_center)[events_by_center>0])
      ff<-model_try(survival::coxph(survival::Surv(z$.time,z$.event)~xx+strata(.center),
        ties="efron",singular.ok=FALSE,control=survival::coxph.control(iter.max=50)))
    } else {
      eventless<-0L;informative<-nrow(z)
      ff<-model_try(survival::coxph(survival::Surv(z$.time,z$.event)~xx,
        ties="efron",singular.ok=FALSE,control=survival::coxph.control(iter.max=50)))
    }
    if(is.null(ff$fit))return(bad("failed_cox",ff$note))
    b<-setNames(as.numeric(coef(ff$fit)),colnames(xx));v<-vcov(ff$fit)
    s<-sqrt(diag(v));p<-2*pnorm(abs(b/s),lower.tail=FALSE)
    names(s)<-names(p)<-names(b)
    er<-exp(c(b[protect],b[protect]-qnorm(.975)*s[protect],b[protect]+qnorm(.975)*s[protect]))
    if(any(!is.finite(er)) || any(er<=0))return(bad("invalid_cox_effect_scale"))
  }
  if(any(!is.finite(c(b,s,p))) || any(s<=0) || any(p<0|p>1))return(bad("invalid_numeric"))
  tab<-cbind(beta=b,se=s,p=p);rownames(tab)<-names(b)
  list(ok=TRUE,status="ok",table=tab,note=paste(notes,ff$note,sep="; "),
       eventless_centers=as.integer(eventless),informative_N=as.integer(informative))
}
model_coef <- function(fit,term,type=NULL) {
  if(is.null(fit)||!isTRUE(fit$ok)||!term%in%rownames(fit$table))return(NULL)
  unname(fit$table[term,c("beta","se","p")])
}
screen_a <- function(d,x,m,covar_names=omics_covars_for(x)) {
  z<-mediation_data(d,x,m,covar_names=covar_names)
  r<-tibble(mediator=m,N=if(is.null(z))NA_integer_ else nrow(z),a=NA_real_,se_a=NA_real_,
    p_a=NA_real_,eligible=!is.null(z),status_a="insufficient_data",note="insufficient data/variation")
  if(is.null(z))return(r)
  fit<-med_model(z,"a");a<-model_coef(fit,".x")
  r$note<-fit$note;r$status_a<-fit$status
  if(!is.null(a)){r$a<-a[1];r$se_a<-a[2];r$p_a<-a[3]}
  r
}
screen_b <- function(d,x,m,y) {
  z<-mediation_data(d,x,m,y)
  r<-tibble(mediator=m,N_b=if(is.null(z))NA_integer_ else nrow(z),b=NA_real_,se_b=NA_real_,
    p_b=NA_real_,eligible=!is.null(z),status_b="insufficient_data",note_b="insufficient data/events/variation",
    eventless_centers=NA_integer_,informative_N=NA_integer_)
  if(is.null(z))return(r)
  fit<-med_model(z,"b");b<-model_coef(fit,".m")
  r$note_b<-fit$note;r$status_b<-fit$status;r$eventless_centers<-fit$eventless_centers;r$informative_N<-fit$informative_N
  if(!is.null(b)){r$b<-b[1];r$se_b<-b[2];r$p_b<-b[3]}
  r
}
empty_mediation <- function() tibble(exposure=character(),mediator=character(),outcome=character(),
  analysis_track=character(),covariate_model=character(),N=integer(),events=integer(),
  a_beta=double(),b_beta=double(),total_beta=double(),direct_beta=double(),indirect_beta=double(),indirect_se=double(),
  p_sobel=double(),p_boot=double(),indirect_lo=double(),indirect_hi=double(),prop_mediated=double(),
  bootstrap_requested=integer(),bootstrap_valid=integer(),bootstrap_sufficient=logical(),exposure_unit=character(),
  status=character(),note=character())
mediate_feature <- function(d,x,m,y,boot=100) {
  z<-mediation_data(d,x,m,y)
  row<-tibble(exposure=x,mediator=m,outcome=y,analysis_track=analysis_track(x),covariate_model=omics_model_label(x),
    N=if(is.null(z))NA_integer_ else nrow(z),events=if(is.null(z))NA_integer_ else sum(z$.event==1),
    a_beta=NA_real_,b_beta=NA_real_,total_beta=NA_real_,direct_beta=NA_real_,indirect_beta=NA_real_,indirect_se=NA_real_,
    p_sobel=NA_real_,p_boot=NA_real_,indirect_lo=NA_real_,indirect_hi=NA_real_,prop_mediated=NA_real_,
    bootstrap_requested=as.integer(boot),bootstrap_valid=0L,bootstrap_sufficient=FALSE,
    exposure_unit=if(x%in%c("bmonth_winter_summer","bmonth_spring_autumn","border12","border1n"))"0 to 1 contrast" else "per SD in common mediation sample",
    status="skipped",note="insufficient data/events/variation")
  if(is.null(z))return(row)
  fa<-med_model(z,"a");ft<-med_model(z,"total");fd<-med_model(z,"b")
  a<-model_coef(fa,".x");t<-model_coef(ft,".x");b<-model_coef(fd,".m");direct<-model_coef(fd,".x")
  row$note<-paste(fa$note,ft$note,fd$note,sep="; ")
  if(any(vapply(list(a,t,b,direct),is.null,logical(1)))){row$status<-"failed";return(row)}
  row$a_beta<-a[1];row$b_beta<-b[1];row$total_beta<-t[1];row$direct_beta<-direct[1]
  row$indirect_beta<-a[1]*b[1];row$indirect_se<-sqrt(b[1]^2*a[2]^2+a[1]^2*b[2]^2)
  if(row$indirect_se>0)row$p_sobel<-2*pnorm(abs(row$indirect_beta/row$indirect_se),lower.tail=FALSE)
  # Descriptive proportions only when total and indirect effects agree and
  # the total Wald CI excludes zero. No causal interpretation is supported.
  if(t[3]<.05 && abs(t[1])>1e-6 && sign(row$indirect_beta)==sign(t[1]) &&
     abs(row$indirect_beta)<=abs(t[1]))row$prop_mediated<-row$indirect_beta/t[1]
  bb<-rep(NA_real_,boot)
  for(i in seq_len(boot)) {
    q<-z[sample.int(nrow(z),replace=TRUE),,drop=FALSE]
    attr(q,"med_covariates")<-attr(z,"med_covariates")
    aa<-model_coef(med_model(q,"a"),".x");bt<-model_coef(med_model(q,"b"),".m")
    if(!is.null(aa)&&!is.null(bt))bb[i]<-aa[1]*bt[1]
  }
  bb<-bb[is.finite(bb)];row$bootstrap_valid<-length(bb)
  row$bootstrap_sufficient<-boot>0 && length(bb)>=max(50,ceiling(.8*boot))
  if(row$bootstrap_sufficient){
    row$p_boot<-min(1,2*min((sum(bb<=0)+1)/(length(bb)+1),(sum(bb>=0)+1)/(length(bb)+1)))
    row$indirect_lo<-quantile(bb,.025,names=FALSE);row$indirect_hi<-quantile(bb,.975,names=FALSE)
  }
  row$status<-"ok";row
}
