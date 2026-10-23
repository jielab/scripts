# Shared descriptive decomposition helpers. Source phewas_engine.R first.
# These functions load no UKB files and can be tested with synthetic inputs.
bm_decomp_complete <- function(d,variables) {
  variables<-unique(variables)
  if(length(setdiff(variables,names(d))))stop('Missing decomposition variables: ',
    paste(setdiff(variables,names(d)),collapse=', '))
  keep<-complete.cases(d[,variables,drop=FALSE])
  for(v in variables)if(is.numeric(d[[v]]))keep<-keep & is.finite(d[[v]])
  droplevels(d[keep,,drop=FALSE])
}
bm_decomp_sample_key <- function(d) {
  if(!'eid' %in% names(d))return(NA_character_)
  bm_hash(sort(as.character(d$eid)))
}
bm_decomp_scalar <- function(d,exposure,covariates,phenotype='') {
  covariates<-setdiff(unique(covariates),exposure)
  z<-bm_decomp_complete(d,c('.y',exposure,covariates))
  x<-if(nrow(z))bm_formula_matrix(z,c(exposure,covariates),tested=exposure) else
    matrix(numeric(),0,2,dimnames=list(NULL,c('(Intercept)',exposure)))
  r<-bm_scalar(list(name=exposure,rows=seq_len(nrow(z)),x=x,
                    dropped_terms=attr(x,'dropped_terms')),as.integer(z$.y),phenotype)
  r$N<-r$n_total;r$cases<-r$n_cases;r$sample_id_hash<-bm_decomp_sample_key(z)
  if(nrow(z)<500L) {
    r$eligible<-FALSE;r$status<-'skipped_low_sample_size'
    r[,c('beta','SE','OR','conf.low','conf.high','p')]<-NA_real_
  }
  r
}
bm_decomp_joint <- function(d,covariates,phenotype='') {
  tested<-c('bmonth_cos','bmonth_sin')
  covariates<-setdiff(unique(covariates),tested)
  z<-bm_decomp_complete(d,c('.y',tested,covariates))
  x<-if(nrow(z))bm_formula_matrix(z,c(covariates,tested),tested=tested) else
    matrix(numeric(),0,3,dimnames=list(NULL,c('(Intercept)',tested)))
  spec<-list(name='cycle2df',rows=seq_len(nrow(z)),x=x,
    null_cols=which(!colnames(x) %in% tested),null_key='decomp',
    dropped_terms=attr(x,'dropped_terms'))
  r<-bm_joint(spec,as.integer(z$.y),phenotype,new.env(parent=emptyenv()))
  r$chisq<-r$statistic;r$sample_id_hash<-bm_decomp_sample_key(z)
  r$cycle_amplitude_logOR<-sqrt(r$beta_cos^2+r$beta_sin^2)
  r$peak_month_descriptive<-NA_integer_
  if(is.finite(r$cycle_amplitude_logOR) && r$cycle_amplitude_logOR>0) {
    vals<- -cos(2*pi*(1:12)/12)*r$beta_cos+sin(2*pi*(1:12)/12)*r$beta_sin
    r$peak_month_descriptive<-which.max(vals)
  }
  r
}
bm_decomp_curve <- function(row,reference=1L) {
  out<-data.frame(month=1:12,reference_month=reference,log_OR=NA_real_,SE=NA_real_,
    OR=NA_real_,conf.low=NA_real_,conf.high=NA_real_,
    interval='95% pointwise from joint covariance; not simultaneous')
  values<-unlist(row[,c('beta_cos','beta_sin','var_cos','var_sin','cov_cos_sin')],use.names=FALSE)
  if(!grepl('^ok_',row$status) || any(!is.finite(values)))return(out)
  b<-c(row$beta_cos,row$beta_sin)
  v<-matrix(c(row$var_cos,row$cov_cos_sin,row$cov_cos_sin,row$var_sin),2,2)
  if(min(eigen(v,symmetric=TRUE,only.values=TRUE)$values)<=0)return(out)
  basis<-cbind(-cos(2*pi*(1:12)/12),sin(2*pi*(1:12)/12))
  cmat<-sweep(basis,2,basis[reference,],'-')
  out$log_OR<-drop(cmat %*% b)
  out$SE<-sqrt(pmax(0,rowSums((cmat %*% v)*cmat)))
  out$OR<-exp(out$log_OR);out$conf.low<-exp(out$log_OR-qnorm(.975)*out$SE)
  out$conf.high<-exp(out$log_OR+qnorm(.975)*out$SE)
  out
}
bm_decomp_ratio <- function(numerator,denominator,method,reference_method) {
  ifelse(is.finite(numerator) & is.finite(denominator) & abs(denominator)>1e-8 &
         !is.na(method) & !is.na(reference_method) & method==reference_method,
         1-numerator/denominator,NA_real_)
}
bm_decomp_screen <- function(path,provenance,cohort_md5) {
  if(!file.exists(path))return(data.frame())
  d<-as.data.frame(readRDS(path))
  req<-c('phenotype','exposure','p','status','eligible','bonferroni_threshold',
         'engine_version','phenotype_panel_md5','cohort_md5')
  if(length(setdiff(req,names(d))))stop('Stale/incomplete PheWAS screen: ',path,'; rerun phewas')
  if(nrow(d) && (anyNA(d[,c('phenotype_panel_md5','cohort_md5','engine_version')]) ||
                any(d$engine_version!=BM_ENGINE_VERSION) ||
                any(d$phenotype_panel_md5!=provenance$panel_md5) ||
                any(d$cohort_md5!=cohort_md5)))
    stop('PheWAS screen input hashes do not match the current cohort/panel: ',path)
  d[!is.na(d$eligible) & d$eligible & grepl('^ok_',d$status) &
      is.finite(d$p) & d$p>=0 & d$p<=1 &
      is.finite(d$bonferroni_threshold) & d$bonferroni_threshold>0,,drop=FALSE]
}

bm_decomp_lm <- function(x,y,tested='.x') {
  prep<-bm_prune_design(x,tested)
  out<-list(ok=FALSE,status='failed_linear',note='',coef=NULL)
  if(!prep$ok){out$status<-prep$status;out$note<-prep$status;return(out)}
  notes<-character()
  fit<-tryCatch(withCallingHandlers(lm.fit(prep$x,y),warning=function(w){
    notes<<-c(notes,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e){notes<<-c(notes,conditionMessage(e));NULL})
  if(is.null(fit) || fit$df.residual<=0L) {out$note<-paste(notes,collapse='; ');return(out)}
  # lm.fit has no terms component; compute its ordinary QR-based covariance
  # directly rather than constructing an incomplete fake lm object.
  ranked<-seq_len(fit$rank)
  sigma2<-sum(fit$residuals^2)/fit$df.residual
  vv<-tryCatch(chol2inv(qr.R(fit$qr)[ranked,ranked,drop=FALSE])*sigma2,
               error=function(e)NULL)
  if(is.null(vv) || any(!is.finite(vv)) || any(diag(vv)<=0))return(out)
  names<-colnames(prep$x)[fit$qr$pivot[ranked]]
  bb<-fit$coefficients[names];ss<-sqrt(diag(vv));tt<-bb/ss
  co<-cbind(Estimate=bb,`Std. Error`=ss,`t value`=tt,
            `Pr(>|t|)`=2*pt(abs(tt),df=fit$df.residual,lower.tail=FALSE))
  rownames(co)<-names
  if(!all(tested %in% rownames(co)) || any(!is.finite(co[tested,,drop=FALSE])))return(out)
  out$ok<-TRUE;out$status<-'ok_lm';out$coef<-co[tested,,drop=FALSE]
  out$note<-paste(notes,collapse='; ');out
}

# Coefficient product on one common sample. This is NOT a natural causal
# indirect effect; logistic non-collapsibility also precludes a mediated share.
bm_temp_product <- function(d,exposure,mediator,covariates,B=100L) {
  if(length(B)!=1L || !is.finite(B) || B<0 || B!=floor(B))stop('Invalid bootstrap count')
  z<-bm_decomp_complete(d,c('.y',exposure,mediator,covariates))
  r<-data.frame(N=nrow(z),cases=sum(z$.y),eligible=FALSE,status='skipped',note='',
    sample_id_hash=bm_decomp_sample_key(z),a_beta=NA_real_,b_beta=NA_real_,
    total_beta=NA_real_,direct_beta=NA_real_,indirect_beta=NA_real_,p_sobel=NA_real_,
    bootstrap_requested=B,bootstrap_attempted=0L,bootstrap_valid=0L,bootstrap_failed=0L,bootstrap_br_fits=0L,
    p_boot=NA_real_,indirect_lo=NA_real_,indirect_hi=NA_real_,attenuation_beta=NA_real_,
    a_status=NA_character_,total_status=NA_character_,direct_status=NA_character_,
    total_method=NA_character_,direct_method=NA_character_,stringsAsFactors=FALSE)
  if(nrow(z)<1000L || r$cases<50L || nrow(z)-r$cases<50L ||
     length(unique(z[[exposure]]))<2L || length(unique(z[[mediator]]))<2L)return(r)
  r$eligible<-TRUE
  z$.x<-as.numeric(scale(z[[exposure]]));z$.m<-as.numeric(scale(z[[mediator]]))
  covariates<-setdiff(unique(covariates),c(exposure,mediator,'.x','.m'))
  xa<-bm_formula_matrix(z,c('.x',covariates),tested='.x')
  xd<-bm_formula_matrix(z,c('.x','.m',covariates),tested=c('.x','.m'))
  fa<-bm_decomp_lm(xa,z$.m)
  ft<-bm_fit(xa,as.integer(z$.y),coefficient='.x')
  fd<-bm_fit(xd,as.integer(z$.y),tested=c('.x','.m'))
  r$a_status<-fa$status;r$total_status<-ft$status;r$direct_status<-fd$status
  r$total_method<-ft$method;r$direct_method<-fd$method
  r$note<-paste(c(fa$note,ft$note,fd$note),collapse='; ')
  if(!fa$ok || !ft$ok || !fd$ok){r$status<-'failed_component_model';return(r)}
  a<-fa$coef['.x',1];ase<-fa$coef['.x',2]
  b<-fd$coefficients['.m'];bse<-sqrt(fd$vcov['.m','.m'])
  total<-ft$coefficients['.x'];direct<-fd$coefficients['.x']
  r$status<-'ok';r$a_beta<-a;r$b_beta<-b;r$total_beta<-total;r$direct_beta<-direct
  r$indirect_beta<-a*b;sobelse<-sqrt(b*b*ase*ase+a*a*bse*bse)
  if(is.finite(sobelse) && sobelse>0)r$p_sobel<-2*pnorm(abs(a*b/sobelse),lower.tail=FALSE)
  r$attenuation_beta<-bm_decomp_ratio(direct,total,fd$method,ft$method)
  products<-rep(NA_real_,B)
  r$bootstrap_attempted<-B
  for(i in seq_len(B)) {
    index<-sample.int(nrow(z),replace=TRUE)
    aa<-bm_decomp_lm(xa[index,,drop=FALSE],z$.m[index])
    bb<-bm_fit(xd[index,,drop=FALSE],as.integer(z$.y[index]),tested=c('.x','.m'))
    if(!aa$ok || !bb$ok)next
    products[i]<-aa$coef['.x',1]*bb$coefficients['.m']
    if(bb$method=='brglm2_AS_mean')r$bootstrap_br_fits<-r$bootstrap_br_fits+1L
  }
  good<-products[is.finite(products)];r$bootstrap_valid<-length(good)
  r$bootstrap_failed<-B-length(good)
  if(B>=50L && length(good)>=max(50L,ceiling(.8*B))) {
    r$p_boot<-min(1,2*min((sum(good<=0)+1)/(length(good)+1),
                         (sum(good>=0)+1)/(length(good)+1)))
    r$indirect_lo<-quantile(good,.025,names=FALSE);r$indirect_hi<-quantile(good,.975,names=FALSE)
  } else if(B>0L)r$note<-paste(r$note,'Bootstrap insufficient; raw bootstrap p remains NA',sep='; ')
  r
}
