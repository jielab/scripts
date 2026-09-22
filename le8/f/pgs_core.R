# Statistical primitives for measured / biomarker-PGS comparisons.
# No individual records are exported. G is captured genetic prediction; R is
# the remainder, not an environmental, lifestyle, or disease-consequence fraction.
PGS_FOCUS_VERSION <- "2026-09-21.matched-locus-v1"
pgs_num <- function(key,default) {x<-suppressWarnings(as.numeric(Sys.getenv(key,as.character(default))));if(length(x)!=1L||is.na(x))stop("Invalid ",key);x}
pgs_csv <- function(key,default="") unique(Filter(nzchar,trimws(strsplit(Sys.getenv(key,default),",",fixed=TRUE)[[1]])))
pgs_bind <- function(x) as.data.frame(data.table::rbindlist(x,fill=TRUE,use.names=TRUE))
pgs_hash <- function(x) {f<-tempfile();on.exit(unlink(f));saveRDS(x,f,compress=FALSE,version=2);unname(tools::md5sum(f))}
pgs_stamp <- function(files) {files<-sort(unique(files[!is.na(files)&nzchar(files)]));data.frame(file=files,size=file.info(files)$size,mtime=as.numeric(file.info(files)$mtime))}
pgs_ids <- function(d,label) {
  if(!"eid"%in%names(d))stop(label," lacks eid")
  d$eid<-as.character(d$eid)
  if(anyNA(d$eid)||any(!nzchar(d$eid))||anyDuplicated(d$eid))stop(label," has missing/duplicate eid")
  d
}
pgs_complete <- function(d,cols) {
  if(length(setdiff(cols,names(d))))return(rep(FALSE,nrow(d)))
  ok<-complete.cases(d[,cols,drop=FALSE])
  for(v in cols)if(is.numeric(d[[v]]))ok<-ok&is.finite(d[[v]])
  ok
}
pgs_folds <- function(group,k=5L,seed=2026L) {
  group<-as.character(group);if(anyNA(group)||any(!nzchar(group)))stop("Missing PGS fold group")
  u<-sort(unique(group));if(length(u)<k)stop("Too few independent fold groups")
  set.seed(seed);f<-sample(rep(seq_len(k),length.out=length(u)));f[match(group,u)]
}
pgs_classify <- function(bm,bg,qm,qg,alpha=.05) {
  out<-rep("Neither supported",length(bm));valid<-is.finite(bm)&is.finite(bg)&is.finite(qm)&is.finite(qg)
  m<-is.finite(qm)&qm<alpha;g<-is.finite(qg)&qg<alpha
  out[m&!g]<-"Measured only";out[g&!m]<-"PGS only"
  out[which(valid&m&g&bm*bg>=0)]<-"Both supported: concordant"
  out[which(valid&m&g&bm*bg<0)]<-"Both supported: opposite"
  out[!valid]<-"Unavailable comparison";out
}
pgs_adjust <- function(d,p="p",groups=character(),out="FDR") {
  if(!nrow(d)||!p%in%names(d))return(d)
  key<-if(length(groups))do.call(paste,c(d[groups],sep="\r"))else rep("all",nrow(d))
  d[[out]]<-NA_real_
  for(g in unique(key)){ii<-which(key==g);d[[out]][ii]<-p.adjust(d[[p]][ii],"BH")};d
}
pgs_cox <- function(d,xs,covars,fold=FALSE,min_events=pgs_num("PGS_MIN_EVENTS",20)) {
  covars<-unique(setdiff(covars,xs));need<-c(".time",".event",xs,covars)
  empty<-data.frame(term=xs,beta=NA_real_,se=NA_real_,p=NA_real_,lo=NA_real_,hi=NA_real_,
    N=nrow(d),events=sum(d$.event==1,na.rm=TRUE),status="insufficient data",warning="",dropped_constant_covariates="")
  result<-list(effects=empty,contrast=data.frame(),loglik=NA_real_,vcov=NULL)
  if(length(setdiff(need,names(d)))){result$effects$status<-paste("missing covariates:",paste(setdiff(need,names(d)),collapse=","));return(result)}
  # Do not silently change the common sample between nested models.
  if(!all(pgs_complete(d,need))||any(d$.time<=0)||any(!d$.event%in%c(0,1))){result$effects$status<-"invalid common model sample";return(result)}
  if(nrow(d)<pgs_num("PGS_MIN_N",200)||sum(d$.event)<min_events)return(result)
  if(any(vapply(d[xs],function(x)!is.finite(sd(x))||sd(x)<=1e-12,logical(1)))){result$effects$status<-"constant exposure";return(result)}
  drop<-covars[vapply(d[covars],function(x)length(unique(x))<2L,logical(1))];covars<-setdiff(covars,drop)
  rhs<-paste(sprintf("`%s`",c(xs,covars)),collapse=" + ")
  if(fold&&length(unique(d$.fold))>1)rhs<-paste(rhs,"strata(.fold)",sep=" + ")
  ff<-as.formula(paste("survival::Surv(.time,.event) ~",rhs));environment(ff)<-environment()
  strata<-survival::strata;warnings<-character()
  args<-list(formula=ff,data=d,ties="efron",model=FALSE,x=FALSE,y=FALSE,
    control=survival::coxph.control(iter.max=40))
  if(".group"%in%names(d)&&anyDuplicated(d$.group)){args$cluster<-d$.group;args$robust<-TRUE}
  fit<-tryCatch(withCallingHandlers(do.call(survival::coxph,args),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
  if(inherits(fit,"error")){result$effects$status<-paste("fit failed:",conditionMessage(fit));return(result)}
  cc<-coef(fit);V<-vcov(fit);ii<-match(xs,names(cc));b<-unname(cc[ii]);s<-sqrt(diag(V)[ii])
  valid<-is.finite(b)&is.finite(s)&s>0
  if(any(grepl("infinite|converg",warnings,ignore.case=TRUE)))valid[]<-FALSE
  result$effects$beta<-ifelse(valid,b,NA_real_);result$effects$se<-ifelse(valid,s,NA_real_)
  result$effects$p<-ifelse(valid,2*pnorm(-abs(b/s)),NA_real_)
  result$effects$lo<-ifelse(valid,b-1.96*s,NA_real_);result$effects$hi<-ifelse(valid,b+1.96*s,NA_real_)
  result$effects$status<-ifelse(valid,"ok","non-estimable or convergence warning")
  result$effects$warning<-paste(unique(warnings),collapse="; ")
  result$effects$dropped_constant_covariates<-paste(drop,collapse=";")
  result$loglik<-fit$loglik[2];result$vcov<-V
  if(all(c(".G",".R")%in%xs)&&all(valid)) {
    delta<-unname(cc[".G"]-cc[".R"]);v<-V[".G",".G"]+V[".R",".R"]-2*V[".G",".R"]
    result$contrast<-data.frame(beta_difference=delta,se_difference=if(v>0)sqrt(v)else NA_real_,
      covariance_GR=V[".G",".R"],p=if(v>0)2*pnorm(-abs(delta/sqrt(v)))else NA_real_,
      lo=if(v>0)delta-1.96*sqrt(v)else NA_real_,hi=if(v>0)delta+1.96*sqrt(v)else NA_real_,
      N=nrow(d),events=sum(d$.event),uncertainty="Conditional on cross-fitted calibration; see refitted bootstrap")
  }
  result
}
pgs_incident <- function(d,cols) d[pgs_complete(d,c(cols,".time",".event",".prev"))&
  d$.prev%in%0&is.finite(d$.time)&d$.time>0&d$.event%in%c(0,1),,drop=FALSE]
pgs_pair <- function(d,covars,feature,adjustment="basic",scope="existing_full") {
  z<-pgs_incident(d,c(".m",".g",covars));id<-pgs_hash(sort(z$eid))
  z$.M<-as.numeric(scale(z$.m));z$.P<-as.numeric(scale(z$.g))
  fits<-list(measured=pgs_cox(z,".M",covars),pgs=pgs_cox(z,".P",covars),joint=pgs_cox(z,c(".M",".P"),covars))
  rows<-lapply(names(fits),function(model){a<-fits[[model]]$effects;a$feature<-feature;a$model<-model;a$adjustment<-adjustment;a$scope<-scope;a$sample_hash<-id;a$covariates<-paste(covars,collapse=";");a$unit<-"log HR per own SD in the identical complete-case sample";a})
  a<-pgs_bind(rows);m<-fits$measured$effects;g<-fits$pgs$effects;j<-fits$joint$effects
  lr<-2*(fits$joint$loglik-fits$measured$loglik)
  s<-data.frame(feature,scope,adjustment,N=nrow(z),events=sum(z$.event),sample_hash=id,
    measured_beta=m$beta,measured_se=m$se,measured_p=m$p,pgs_beta=g$beta,pgs_se=g$se,pgs_p=g$p,
    measured_joint_beta=j$beta[match(".M",j$term)],pgs_joint_beta=j$beta[match(".P",j$term)],
    pgs_joint_p=j$p[match(".P",j$term)],pgs_increment_LRT_p=if(!anyDuplicated(z$.group)&&is.finite(lr))pchisq(max(0,lr),1,lower.tail=FALSE)else NA_real_,
    LRT_status=if(anyDuplicated(z$.group))"withheld for clustered observations; use robust joint Wald"else"iid partial-likelihood diagnostic",
    correlation=if(nrow(z)>2&&all(is.finite(z$.M))&&all(is.finite(z$.P)))cor(z$.M,z$.P)else NA_real_,
    covariates=paste(covars,collapse=";"),comparison="identical people, follow-up and adjustment; own-SD associations, not MR")
  list(effects=a,summary=s)
}
pgs_calibrate <- function(d,covars,k=5L,seed=2026L,diagnostics=TRUE) {
  d<-as.data.frame(d)
  if(!".group"%in%names(d))d$.group<-d$eid
  if(!".fold"%in%names(d))d$.fold<-pgs_folds(d$.group,k,seed)
  for(nm in c(".M",".G",".R",".pred",".pred0"))d[[nm]]<-NA_real_
  reports<-list();ok<-pgs_complete(d,c(".m",".g",covars))
  for(f in sort(unique(d$.fold))) {
    train<-d[d$.fold!=f&d$.prev%in%0&ok,,drop=FALSE];ii<-which(d$.fold==f&ok)
    if(nrow(train)<100||!length(ii)||sd(train$.m)<=1e-12||sd(train$.g)<=1e-12)next
    om<-mean(train$.m);os<-sd(train$.m);gm<-mean(train$.g);gs<-sd(train$.g)
    train$.y<-(train$.m-om)/os;train$.p<-(train$.g-gm)/gs
    cv<-covars[vapply(train[covars],function(x)length(unique(x))>1,logical(1))]
    fit<-tryCatch(lm(reformulate(c(".p",cv),".y"),train),error=function(e)NULL)
    fit0<-if(diagnostics)tryCatch(lm(reformulate(if(length(cv))cv else "1",".y"),train),error=function(e)NULL)else NULL
    if(is.null(fit)||(diagnostics&&is.null(fit0))||!is.finite(coef(fit)[".p"]))next
    test<-d[ii,,drop=FALSE];test$.p<-(test$.g-gm)/gs;b<-unname(coef(fit)[".p"])
    p<-if(diagnostics)tryCatch(as.numeric(predict(fit,test)),error=function(e)rep(NA_real_,nrow(test)))else rep(NA_real_,nrow(test))
    p0<-if(diagnostics)tryCatch(as.numeric(predict(fit0,test)),error=function(e)rep(NA_real_,nrow(test)))else rep(NA_real_,nrow(test))
    d$.M[ii]<-(test$.m-om)/os;d$.G[ii]<-b*test$.p;d$.R[ii]<-d$.M[ii]-d$.G[ii]
    d$.pred[ii]<-p;d$.pred0[ii]<-p0
    reports[[length(reports)+1L]]<-data.frame(fold=f,N_train=nrow(train),N_test=length(ii),slope=b,
      omic_mean=om,omic_sd=os,pgs_mean=gm,pgs_sd=gs,
      calibration="baseline disease-free training groups; future outcomes not used")
  }
  ii<-d$.prev%in%0&pgs_complete(d,c(".M",".pred",".pred0"))
  mse0<-sum((d$.M[ii]-d$.pred0[ii])^2);r2<-if(mse0>0)1-sum((d$.M[ii]-d$.pred[ii])^2)/mse0 else NA_real_
  folds<-pgs_bind(reports);b<-folds$slope
  audit<-data.frame(N_calibration=sum(ii),partial_R2=r2,
    slope_mean=if(length(b))mean(b)else NA_real_,slope_min=if(length(b))min(b)else NA_real_,slope_max=if(length(b))max(b)else NA_real_,
    orientation=if(!length(b))"unavailable"else if(all(b>0))"positive in all folds"else if(all(b<0))"negative in all folds"else "unstable across folds",
    weak_capture=!is.finite(r2)||r2<pgs_num("PGS_MIN_PARTIAL_R2",.005),
    identity_max_error=if(any(is.finite(d$.M)))max(abs(d$.M-d$.G-d$.R),na.rm=TRUE)else NA_real_,
    calibration_folds=length(b),unit="training-fold whole-biomarker SD; G/R not separately standardized")
  list(data=d,folds=folds,audit=audit)
}
pgs_component_models <- function(d,covars,feature,scope="existing_full",landmark=0,end=Inf) {
  z<-pgs_incident(d,c(".M",".G",".R",covars));z<-z[z$.time>landmark,,drop=FALSE]
  z$.event<-as.integer(z$.event==1&z$.time<=end);z$.time<-pmin(z$.time,end)-landmark
  fits<-list(measured=pgs_cox(z,".M",covars,TRUE),genetic=pgs_cox(z,".G",covars,TRUE),
    remaining=pgs_cox(z,".R",covars,TRUE),joint=pgs_cox(z,c(".G",".R"),covars,TRUE))
  effects<-pgs_bind(lapply(names(fits),function(nm){a<-fits[[nm]]$effects;a$model<-nm;a}))
  for(nm in c("effects","contrast")) {
    a<-if(nm=="effects")effects else fits$joint$contrast
    if(nrow(a)){a$feature<-feature;a$scope<-scope;a$landmark<-landmark;a$end<-end;a$sample_hash<-pgs_hash(sort(z$eid));a$covariates<-paste(covars,collapse=";")}
    if(nm=="effects")effects<-a else contrast<-a
  }
  list(effects=effects,contrast=contrast)
}
pgs_bootstrap <- function(d,covars,feature,B=100L,k=5L,seed=2026L) {
  if(B<1)return(data.frame())
  if(!".group"%in%names(d))d$.group<-d$eid
  if(!".fold"%in%names(d))d$.fold<-pgs_folds(d$.group,k,seed)
  groups<-split(seq_len(nrow(d)),d$.group);rows<-vector("list",B)
  for(b in seq_len(B)) {
    set.seed(seed+b);ii<-unlist(groups[sample(seq_along(groups),length(groups),replace=TRUE)],use.names=FALSE)
    # Repeated copies of an individual/family always remain in the same fold.
    cal<-pgs_calibrate(d[ii,,drop=FALSE],covars,k,seed,diagnostics=FALSE)
    fitdata<-pgs_incident(cal$data,c(".G",".R",covars))
    res<-pgs_cox(fitdata,c(".G",".R"),covars,TRUE)
    z<-res$effects;delta<-res$contrast$beta_difference
    rows[[b]]<-data.frame(replicate=b,genetic=if(nrow(z))z$beta[match(".G",z$term)]else NA_real_,
      remaining=if(nrow(z))z$beta[match(".R",z$term)]else NA_real_,difference=if(length(delta))delta[1]else NA_real_)
  }
  draws<-pgs_bind(rows)
  pgs_bind(lapply(c("genetic","remaining","difference"),function(term){x<-draws[[term]];x<-x[is.finite(x)];n<-length(x);valid<-n>=max(30,ceiling(.8*B))
    data.frame(feature,term,requested=B,successful=n,lo=if(valid)unname(quantile(x,.025))else NA_real_,
      hi=if(valid)unname(quantile(x,.975))else NA_real_,
      sign_tail=if(valid)min(1,2*min((sum(x<=0)+1)/(n+1),(sum(x>=0)+1)/(n+1)))else NA_real_,
      status=if(valid)"ok"else"too few successful refits",
      uncertainty="Participant/family resampling with calibration refit; conditional on selected candidate and original GWAS weights")
  }))
}
pgs_lifestyle <- function(d,covars,components,feature) {
  pgs_bind(lapply(components,function(component)pgs_bind(lapply(c(".M",".G",".R"),function(part){
    cols<-unique(c(part,component,covars));z<-d[d$.prev%in%0&pgs_complete(d,cols),,drop=FALSE]
    ans<-data.frame(feature,component,part,N=nrow(z),beta=NA_real_,se=NA_real_,p=NA_real_,status="unavailable")
    if(nrow(z)<200||!component%in%names(z)||sd(z[[component]])<=0||sd(z[[part]])<=0)return(ans)
    z$.x<-as.numeric(scale(z[[component]]));z$.y<-z[[part]]
    cv<-setdiff(covars,component);cv<-cv[vapply(z[cv],function(x)length(unique(x))>1,logical(1))]
    f<-tryCatch(lm(reformulate(c(".x",cv),".y"),z),error=function(e)NULL)
    if(is.null(f))return(ans);s<-coef(summary(f));if(!".x"%in%rownames(s))return(ans)
    ans$beta<-s[".x",1];ans$se<-s[".x",2];ans$p<-s[".x",4]
    if(anyDuplicated(z$.group)) {
      X<-model.matrix(f);B<-tryCatch(solve(crossprod(X)),error=function(e)NULL)
      if(!is.null(B)) {
        u<-rowsum(X*as.numeric(residuals(f)),z$.group,reorder=FALSE);ng<-nrow(u)
        V<-B%*%crossprod(u)%*%B
        ans$se<-sqrt(V[".x",".x"]*ng/(ng-1)*(nrow(z)-1)/(nrow(z)-ncol(X)))
        ans$p<-2*pt(-abs(ans$beta/ans$se),df=ng-1)
      }else {ans$se<-NA_real_;ans$p<-NA_real_}
    }
    ans$status<-"association only; not mediation or environmental fraction";ans
  }))))
}
pgs_prevalent <- function(d,covars,feature) {
  # Baseline case-control association is kept separate from incident hazards.
  z<-d[pgs_complete(d,c(".M",".G",".R",".prev",covars))&d$.prev%in%c(0,1),,drop=FALSE]
  if(nrow(z)<pgs_num("PGS_MIN_N",200)||min(table(factor(z$.prev,levels=0:1)))<pgs_num("PGS_MIN_EVENTS",20))return(data.frame())
  cv<-covars[vapply(z[covars],function(x)length(unique(x))>1,logical(1))]
  pgs_bind(lapply(list(measured=".M",joint=c(".G",".R")),function(xs){
    if(any(vapply(z[xs],function(x)sd(x)<=1e-12,logical(1))))return(data.frame())
    warns<-character();f<-tryCatch(withCallingHandlers(glm(reformulate(c(xs,cv),".prev"),data=z,family=binomial()),
      warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)NULL)
    if(is.null(f)||!f$converged||length(warns))return(data.frame(feature,status="prevalent logistic non-estimable"))
    V<-vcov(f)
    if(anyDuplicated(z$.group)) {
      X<-model.matrix(f);score<-X*as.numeric(z$.prev-fitted(f));u<-rowsum(score,z$.group,reorder=FALSE)
      V<-V%*%crossprod(u)%*%V
    }
    b<-coef(f)[xs];se<-sqrt(diag(V)[xs]);data.frame(feature,term=xs,model=if(length(xs)==1)"measured"else"joint",
      beta=unname(b),se=unname(se),p=2*pnorm(-abs(b/se)),lo=b-1.96*se,hi=b+1.96*se,
      N=nrow(z),cases=sum(z$.prev),status="ok",unit="log OR per training-fold biomarker SD; not comparable in magnitude with log HR")
  }))
}
