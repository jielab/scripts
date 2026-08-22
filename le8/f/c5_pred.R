# C5-specific nested cross-validation and omics prediction utilities.
# Generic prediction helpers (including fcidx and RE_pred) come from
# /mnt/d/scripts/0f/pred.f.R, sourced by comm.f.R.

impute_train_test <- function(tr,te,vars){
  vars<-intersect(vars,intersect(names(tr),names(te)));if(!length(vars))return(list(tr=matrix(nrow=nrow(tr),ncol=0),te=matrix(nrow=nrow(te),ncol=0),vars=character()))
  a<-tr[,vars,drop=FALSE];b<-te[,vars,drop=FALSE]
  a[]<-lapply(a,function(x)suppressWarnings(as.numeric(x)));b[]<-lapply(b,function(x)suppressWarnings(as.numeric(x)))
  keep<-vapply(a,function(x)sum(is.finite(x))>=100&&is.finite(sd(x,na.rm=TRUE))&&sd(x,na.rm=TRUE)>0,logical(1));a<-a[,keep,drop=FALSE];b<-b[,keep,drop=FALSE]
  if(!ncol(a))return(list(tr=matrix(nrow=nrow(tr),ncol=0),te=matrix(nrow=nrow(te),ncol=0),vars=character()))
  cen<-vapply(a,function(x)median(x[is.finite(x)],na.rm=TRUE),numeric(1));sca<-vapply(a,function(x)sd(x,na.rm=TRUE),numeric(1));sca[!is.finite(sca)|sca==0]<-1
  for(j in seq_along(cen)){a[!is.finite(a[,j]),j]<-cen[j];b[!is.finite(b[,j]),j]<-cen[j]}
  list(tr=scale(as.matrix(a),center=cen,scale=sca),te=scale(as.matrix(b),center=cen,scale=sca),vars=colnames(a))
}

auc_binary <- function(y,score){
  id<-which(is.finite(y)&is.finite(score));y<-y[id];score<-score[id];if(length(y)<50||length(unique(y))<2)return(NA_real_)
  r<-rank(score);n1<-sum(y==1);n0<-sum(y==0);(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)
}

ipcw_at_horizon <- function(time,event,horizon){
  ok<-is.finite(time)&is.finite(event)&time>0;time<-time[ok];event<-event[ok]
  sf<-tryCatch(survival::survfit(Surv(time,1-event)~1),error=function(e)NULL)
  if(is.null(sf))return(list(ok=ok,weight=rep(NA_real_,length(time)),case=rep(FALSE,length(time)),control=rep(FALSE,length(time))))
  Gfun<-stats::stepfun(sf$time,c(1,sf$surv),right=TRUE)
  Gt<-pmax(Gfun(horizon),.02);case<-event==1&time<=horizon;control<-time>horizon
  w<-rep(0,length(time));w[case]<-1/pmax(Gfun(pmax(time[case]-1e-8,0)),.02);w[control]<-1/Gt
  list(ok=ok,weight=w,case=case,control=control)
}
weighted_time_auc <- function(time,event,score,horizon){
  id<-which(is.finite(time)&is.finite(event)&is.finite(score)&time>0);if(length(id)<100)return(NA_real_)
  ip<-ipcw_at_horizon(time[id],event[id],horizon);d<-data.frame(score=score[id],cw=ifelse(ip$case,ip$weight,0),tw=ifelse(ip$control,ip$weight,0))
  d<-d[d$cw>0|d$tw>0,,drop=FALSE];if(sum(d$cw)==0||sum(d$tw)==0)return(NA_real_)
  g<-dplyr::as_tibble(d)|>dplyr::group_by(score)|>dplyr::summarise(cw=sum(cw),tw=sum(tw),.groups="drop")|>dplyr::arrange(score)|>dplyr::mutate(tw_before=dplyr::lag(cumsum(tw),default=0))
  sum(g$cw*(g$tw_before+.5*g$tw))/(sum(g$cw)*sum(g$tw))
}
ipcw_brier <- function(time,event,risk,horizon){
  id<-which(is.finite(time)&is.finite(event)&is.finite(risk)&time>0);if(length(id)<100)return(NA_real_)
  ip<-ipcw_at_horizon(time[id],event[id],horizon);target<-as.numeric(ip$control);w<-ip$weight;sum(w*(target-(1-risk[id]))^2)/sum(w)
}

calibration_metrics <- function(y,risk){
  id<-which(is.finite(y)&is.finite(risk));if(length(id)<100||length(unique(y[id]))<2)return(c(intercept=NA_real_,slope=NA_real_))
  lp<-qlogis(pmin(pmax(risk[id],1e-6),1-1e-6));slope<-tryCatch(coef(glm(y[id]~lp,family=binomial()))[[2]],error=function(e)NA_real_)
  intercept<-tryCatch(coef(glm(y[id]~offset(lp),family=binomial()))[[1]],error=function(e)NA_real_);c(intercept=intercept,slope=slope)
}

fit_gaussian_glmnet <- function(tr,te,vars,yvar,nfolds=5,alpha=.5){
  if(!requireNamespace("glmnet",quietly=TRUE))stop("Prediction requires glmnet.",call.=FALSE)
  x<-impute_train_test(tr,te,vars);if(ncol(x$tr)<1)return(NULL);y<-suppressWarnings(as.numeric(tr[[yvar]]));ok<-is.finite(y)
  if(sum(ok)<200||sd(y[ok])==0)return(NULL)
  fit<-tryCatch(glmnet::cv.glmnet(x$tr[ok,,drop=FALSE],y[ok],family="gaussian",alpha=alpha,nfolds=min(nfolds,max(3,floor(sum(ok)/50))),standardize=FALSE,keep=TRUE),error=function(e)NULL);if(is.null(fit))return(NULL)
  b<-as.matrix(coef(fit,s="lambda.1se"));sel<-rownames(b)[b[,1]!=0];sel<-setdiff(sel,"(Intercept)")
  # Use inner-fold out-of-fold predictions for the outer-training rows. This prevents
  # the second-stage Cox model from being calibrated on in-sample omics scores.
  li<-which.min(abs(fit$lambda-fit$lambda.1se));trp<-rep(NA_real_,nrow(tr))
  if(!is.null(fit$fit.preval)&&length(li)) trp[ok]<-as.numeric(fit$fit.preval[,li])
  if(sum(is.finite(trp))<sum(ok)*.8) trp[ok]<-as.numeric(predict(fit,x$tr[ok,,drop=FALSE],s="lambda.1se"))
  list(tr=trp,te=as.numeric(predict(fit,x$te,s="lambda.1se")),selected=sel)
}
fit_cox_glmnet <- function(tr,te,vars,tvar,evar,nfolds=5,alpha=.5,penalty_factor=NULL){
  if(!requireNamespace("glmnet",quietly=TRUE))stop("Prediction requires glmnet.",call.=FALSE)
  x<-impute_train_test(tr,te,vars);if(ncol(x$tr)<2)return(NULL);time<-tr[[tvar]];event<-tr[[evar]];ok<-is.finite(time)&is.finite(event)&time>0
  if(sum(ok)<300||sum(event[ok]==1)<30)return(NULL)
  pf<-rep(1,ncol(x$tr));names(pf)<-colnames(x$tr);if(!is.null(penalty_factor)){hit<-intersect(names(penalty_factor),names(pf));pf[hit]<-penalty_factor[hit]}
  fit<-tryCatch(glmnet::cv.glmnet(x$tr[ok,,drop=FALSE],Surv(time[ok],event[ok]),family="cox",alpha=alpha,nfolds=min(nfolds,max(3,floor(sum(event[ok])/10))),standardize=FALSE,penalty.factor=pf,keep=TRUE),error=function(e)NULL);if(is.null(fit))return(NULL)
  b<-as.matrix(coef(fit,s="lambda.1se"));sel<-rownames(b)[b[,1]!=0]
  li<-which.min(abs(fit$lambda-fit$lambda.1se));trp<-rep(NA_real_,nrow(tr))
  if(!is.null(fit$fit.preval)&&length(li)) trp[ok]<-as.numeric(fit$fit.preval[,li])
  if(sum(is.finite(trp))<sum(ok)*.8) trp[ok]<-as.numeric(predict(fit,x$tr[ok,,drop=FALSE],s="lambda.1se",type="link"))
  list(tr=trp,te=as.numeric(predict(fit,x$te,s="lambda.1se",type="link")),selected=sel,beta=b[sel,1])
}

discover_ys_training <- function(tr,biom_vars,le8_vars,basic_vars,block=100,fdr=.05,specificity_cut=.35,max_per_component=40,seed=2026){
  biom_vars<-intersect(biom_vars,names(tr));le8_vars<-intersect(le8_vars,names(tr));basic_vars<-intersect(basic_vars,names(tr))
  if(!length(biom_vars)||!length(le8_vars))return(list())
  d<-tr[complete.cases(tr[,unique(c(le8_vars,basic_vars)),drop=FALSE]),,drop=FALSE];if(nrow(d)<1000)return(list())
  set.seed(seed);split<-sample(rep(c("discovery","replication"),length.out=nrow(d)))
  scan_half<-function(dd){
    if(nrow(dd)<400)return(tibble::tibble())
    blocks<-split(biom_vars,ceiling(seq_along(biom_vars)/block));rows<-list();k<-0L
    for(cmp in le8_vars){
      mm<-model.matrix(reformulate(unique(c(basic_vars,setdiff(le8_vars,cmp)))),dd);q<-qr(mm);yr<-qr.resid(q,as.numeric(scale(as.numeric(dd[[cmp]]))));sy<-sqrt(sum(yr^2));df<-max(3,nrow(dd)-ncol(mm)-2)
      for(bb in blocks){
        x<-as.matrix(dd[,bb,drop=FALSE]);storage.mode(x)<-"double"
        for(j in seq_len(ncol(x))){m<-median(x[,j],na.rm=TRUE);if(!is.finite(m))m<-0;x[!is.finite(x[,j]),j]<-m;sx<-sd(x[,j],na.rm=TRUE);x[,j]<-if(is.finite(sx)&&sx>0)as.numeric(scale(x[,j]))else 0}
        xr<-qr.resid(q,x);r<-as.numeric(crossprod(xr,yr)/pmax(sqrt(colSums(xr^2))*sy,1e-12));r<-pmax(pmin(r,.999),-.999);z<-r*sqrt(df/pmax(1-r^2,1e-9));k<-k+1L
        rows[[k]]<-tibble::tibble(feature=bb,component=sub("\\.pts$","",cmp),component_var=cmp,r=r,z=z,p=2*pt(abs(z),df=df,lower.tail=FALSE))
      }
    }
    dplyr::bind_rows(rows)|>dplyr::group_by(component)|>dplyr::mutate(FDR=p.adjust(p,"BH"))|>dplyr::ungroup()
  }
  ad<-scan_half(d[split=="discovery",,drop=FALSE]);ar<-scan_half(d[split=="replication",,drop=FALSE])
  if(!nrow(ad)||!nrow(ar))return(list())
  spec<-ad|>dplyr::group_by(feature)|>dplyr::mutate(absz=abs(z),specificity=ifelse(sum(absz,na.rm=TRUE)>0,absz/sum(absz,na.rm=TRUE),NA_real_))|>
    dplyr::slice_max(absz,n=1,with_ties=FALSE)|>dplyr::ungroup()|>dplyr::select(feature,primary_component=component,specificity)
  pri<-ad|>dplyr::select(feature,component,component_var,r_disc=r,z_disc=z,FDR_disc=FDR)|>
    dplyr::left_join(ar|>dplyr::select(feature,component,r_rep=r,z_rep=z,FDR_rep=FDR),by=c("feature","component"))|>
    dplyr::left_join(spec,by="feature")|>dplyr::filter(component==primary_component)|>
    dplyr::mutate(same_direction=is.finite(r_rep)&sign(r_disc)==sign(r_rep),strict=is.finite(FDR_disc)&is.finite(FDR_rep)&FDR_disc<fdr&FDR_rep<fdr&same_direction&specificity>=specificity_cut,
                  combined_fdr=dplyr::case_when(is.finite(FDR_disc)&is.finite(FDR_rep)~pmax(FDR_disc,FDR_rep),is.finite(FDR_disc)~FDR_disc,TRUE~Inf))
  ans<-lapply(le8_vars,function(cmp){
    nm<-sub("\\.pts$","",cmp);z<-pri|>dplyr::filter(component==nm,strict)|>dplyr::arrange(combined_fdr,dplyr::desc(abs(r_disc)))|>dplyr::slice_head(n=max_per_component)|>dplyr::pull(feature)
    if(length(z)<2)z<-pri|>dplyr::filter(component==nm)|>dplyr::arrange(dplyr::desc(same_direction),combined_fdr,dplyr::desc(specificity),dplyr::desc(abs(r_disc)))|>dplyr::slice_head(n=min(10,max_per_component))|>dplyr::pull(feature)
    unique(z)
  });names(ans)<-sub("\\.pts$","",le8_vars)
  chosen <- dplyr::bind_rows(lapply(names(ans), function(nm) {
    if (!length(ans[[nm]])) return(NULL)
    tibble::tibble(component = nm, feature = ans[[nm]])
  }))
  if (!nrow(chosen)) chosen <- tibble::tibble(component=character(),feature=character(),strict=logical(),FDR_disc=numeric(),FDR_rep=numeric(),specificity=numeric(),proxy_level=character())
  else chosen <- chosen |> dplyr::left_join(pri |> dplyr::select(feature,component,strict,FDR_disc,FDR_rep,specificity),by=c("feature","component")) |>
    dplyr::mutate(proxy_level=ifelse(strict,"YS_strict_fold","YS_exploratory_fold"))
  attr(ans,"selection_info")<-chosen;ans
}

cox_predict <- function(tr,te,vars,tvar,evar,horizon=10){
  vars<-intersect(vars,intersect(names(tr),names(te)));if(!length(vars))return(list(lp=rep(NA_real_,nrow(te)),risk=rep(NA_real_,nrow(te))))
  f<-as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(bt(vars),collapse=" + ")))
  dtr<-tr[,c(tvar,evar,vars),drop=FALSE];dtr<-dtr[complete.cases(dtr),,drop=FALSE];if(nrow(dtr)<200||sum(dtr[[evar]]==1)<20)return(list(lp=rep(NA_real_,nrow(te)),risk=rep(NA_real_,nrow(te))))
  fit<-tryCatch(coxph(f,dtr,x=TRUE),error=function(e)NULL);if(is.null(fit))return(list(lp=rep(NA_real_,nrow(te)),risk=rep(NA_real_,nrow(te))))
  id<-complete.cases(te[,vars,drop=FALSE]);lp<-risk<-rep(NA_real_,nrow(te));lp[id]<-tryCatch(as.numeric(predict(fit,newdata=te[id,,drop=FALSE],type="lp")),error=function(e)NA_real_)
  bh<-basehaz(fit,centered=TRUE)
  H0<-if(!nrow(bh)||horizon<min(bh$time,na.rm=TRUE)) 0 else bh$hazard[max(which(bh$time<=horizon))]
  if(!is.finite(H0)) H0<-tail(bh$hazard[is.finite(bh$hazard)],1)
  risk[id]<-1-exp(-H0*exp(lp[id]));list(lp=lp,risk=pmin(pmax(risk,0),1),fit=fit)
}

nested_cv_omics <- function(dat,biom_vars,le8_vars,basic_vars,prs_vars,tvar,evar,k=5,inner=5,horizon=10,seed=2026){
  dat$fold_id<-make_folds(dat,evar,k,seed);pred_rows<-list();sel_rows<-list();ii<-0L
  for(fd in seq_len(k)){
    tr<-dat[dat$fold_id!=fd,,drop=FALSE];te<-dat[dat$fold_id==fd,,drop=FALSE];idx<-which(dat$fold_id==fd)
    ys_sets<-discover_ys_training(tr,biom_vars,le8_vars,basic_vars,seed=seed+fd);ys_info<-attr(ys_sets,"selection_info")%||%tibble::tibble();ys_selected<-unique(unlist(ys_sets));ys_tr<-ys_te<-data.frame(row.names=seq_len(nrow(tr)));ys_te<-data.frame(row.names=seq_len(nrow(te)))
    for(cmp in names(ys_sets)){v<-ys_sets[[cmp]];fit<-fit_gaussian_glmnet(tr,te,v,paste0(cmp,".pts"),inner);vn<-paste0("YS_",cmp);ys_tr[[vn]]<-if(is.null(fit))NA_real_ else fit$tr;ys_te[[vn]]<-if(is.null(fit))NA_real_ else fit$te
      if(!is.null(fit)&&length(fit$selected)){si<-tibble::tibble(fold=fd,set="YS",component=cmp,feature=fit$selected)|>dplyr::left_join(ys_info|>dplyr::select(feature,component,proxy_level),by=c("feature","component"));sel_rows[[length(sel_rows)+1]]<-si}}
    tr2<-dplyr::bind_cols(tr,ys_tr);te2<-dplyr::bind_cols(te,ys_te);ys_score_vars<-grep("^YS_",names(tr2),value=TRUE);ys_score_vars<-ys_score_vars[vapply(tr2[ys_score_vars],function(x)sum(is.finite(x))>=200,logical(1))]
    ns<-fit_cox_glmnet(tr,te,biom_vars,tvar,evar,inner);tr2$NS_score<-if(is.null(ns))NA_real_ else ns$tr;te2$NS_score<-if(is.null(ns))NA_real_ else ns$te;if(!is.null(ns)&&length(ns$selected))sel_rows[[length(sel_rows)+1]]<-data.frame(fold=fd,set="NS",component=NA,feature=ns$selected)
    plus_pool<-setdiff(biom_vars,ys_selected);plus<-fit_cox_glmnet(tr,te,plus_pool,tvar,evar,inner);plus_sel<-if(is.null(plus))character() else plus$selected
    union<-unique(c(ys_selected,plus_sel));pf<-setNames(rep(1,length(union)),union);pf[ys_selected]<-.6;ysp<-fit_cox_glmnet(tr,te,union,tvar,evar,inner,penalty_factor=pf);tr2$YSP_score<-if(is.null(ysp))NA_real_ else ysp$tr;te2$YSP_score<-if(is.null(ysp))NA_real_ else ysp$te
    ysp_plus_selected<-if(is.null(ysp)) character() else intersect(ysp$selected,plus_sel)
    if(length(ysp_plus_selected))sel_rows[[length(sel_rows)+1]]<-data.frame(fold=fd,set="YSP_plus",component=NA,feature=ysp_plus_selected)
    models<-list(Base=basic_vars,PRS=c(basic_vars,prs_vars),LE8=c(basic_vars,le8_vars),YS=c(basic_vars,ys_score_vars),LE8_YS=c(basic_vars,le8_vars,ys_score_vars),PRS_LE8=c(basic_vars,prs_vars,le8_vars),PRS_LE8_YS=c(basic_vars,prs_vars,le8_vars,ys_score_vars),PRS_LE8_YSP=c(basic_vars,prs_vars,le8_vars,"YSP_score"),PRS_LE8_NS=c(basic_vars,prs_vars,le8_vars,"NS_score"))
    for(nm in names(models)){z<-cox_predict(tr2,te2,models[[nm]],tvar,evar,horizon);ii<-ii+1;pred_rows[[ii]]<-data.frame(row_id=idx,fold=fd,model=nm,lp=z$lp,risk=z$risk,time=te[[tvar]],event=te[[evar]])}
  }
  pred<-dplyr::bind_rows(pred_rows);selected<-dplyr::bind_rows(sel_rows)
  if(!nrow(selected)) selected<-tibble::tibble(fold=integer(),set=character(),component=character(),feature=character(),proxy_level=character())
  metrics<-pred|>dplyr::group_by(fold,model)|>dplyr::group_modify(function(d,key){eligible<-d$time>=horizon|(d$event==1&d$time<=horizon);y10<-as.integer(d$event==1&d$time<=horizon);cal<-calibration_metrics(y10[eligible],d$risk[eligible]);tibble::tibble(cindex=fcidx(Surv(d$time,d$event),d$lp),AUC10=weighted_time_auc(d$time,d$event,d$risk,horizon),Brier10=ipcw_brier(d$time,d$event,d$risk,horizon),cal_intercept=cal[["intercept"]],cal_slope=cal[["slope"]])})|>dplyr::ungroup()
  summary<-metrics|>dplyr::group_by(model)|>dplyr::summarise(dplyr::across(c(cindex,AUC10,Brier10,cal_intercept,cal_slope),list(mean=~mean(.x,na.rm=TRUE),sd=~sd(.x,na.rm=TRUE))),.groups="drop")
  stability<-selected|>dplyr::count(set,component,feature,name="fold_count")|>dplyr::mutate(selection_frequency=fold_count/k)|>dplyr::arrange(set,component,dplyr::desc(selection_frequency))
  list(pred=pred,metrics=metrics,summary=summary,selected=selected,stability=stability,folds=dat$fold_id)
}

risk_deciles <- function(pred,model){
  d<-pred|>dplyr::filter(.data$model==.env$model,is.finite(lp),is.finite(time),is.finite(event),time>0)|>dplyr::mutate(decile=dplyr::ntile(lp,10),decile=factor(decile))
  fit<-tryCatch(coxph(Surv(time,event)~decile,d),error=function(e)NULL);hr<-if(is.null(fit))tibble::tibble() else {sm<-coef(summary(fit));tibble::tibble(term=rownames(sm),HR=exp(sm[,"coef"]),lo=exp(sm[,"coef"]-1.96*sm[,"se(coef)"]),hi=exp(sm[,"coef"]+1.96*sm[,"se(coef)"]),p=sm[,"Pr(>|z|)"])}
  cum<-d|>dplyr::group_by(decile)|>dplyr::summarise(N=dplyr::n(),events=sum(event),risk10=mean(risk,na.rm=TRUE),.groups="drop");list(rows=d,HR=hr,cum=cum)
}

# Nature Aging-style sequential panel: C1-significant biomarkers are ordered by
# Wald P value, added one at a time, and compared with a paired DeLong test on a
# held-out subset. Stop after two consecutive non-significant AUC gains.
sequential_forward_panel <- function(dat, ranked, tvar, evar, horizon = 10, seed = 2026,
                                     alpha = .05, max_extended = 30) {
  ranked <- intersect(unique(ranked), names(dat)); extended <- head(ranked, max_extended)
  empty_log <- tibble::tibble(step=integer(),feature=character(),auc=numeric(),delta_auc=numeric(),delong_p=numeric(),accepted=logical())
  if (!length(extended) || !requireNamespace("pROC", quietly=TRUE))
    return(list(parsimonious=head(extended,min(5,length(extended))),extended=extended,log=empty_log))
  set.seed(seed); test_id <- sample(seq_len(nrow(dat)), max(1L, round(.2*nrow(dat))))
  eligible <- is.finite(dat[[tvar]]) & !is.na(dat[[evar]]) & (dat[[tvar]] >= horizon | dat[[evar]] == 1)
  tr <- dat[-test_id,,drop=FALSE]; te <- dat[test_id,,drop=FALSE]
  tr <- tr[eligible[-test_id],,drop=FALSE]; te <- te[eligible[test_id],,drop=FALSE]
  ytr <- as.integer(tr[[evar]] == 1 & tr[[tvar]] <= horizon); yte <- as.integer(te[[evar]] == 1 & te[[tvar]] <= horizon)
  if (length(unique(ytr))<2 || length(unique(yte))<2) return(list(parsimonious=head(extended,min(5,length(extended))),extended=extended,log=empty_log))
  med <- vapply(tr[,extended,drop=FALSE],function(x)median(as.numeric(x),na.rm=TRUE),numeric(1))
  prep <- function(d){x<-as.data.frame(lapply(d[,extended,drop=FALSE],as.numeric));names(x)<-extended;for(v in extended)x[[v]][!is.finite(x[[v]])]<-med[[v]];x}
  xtr <- prep(tr); xte <- prep(te); old_roc <- NULL; last_sig <- 0L; nonsig <- 0L; stop_at <- NA_integer_; rows <- list()
  for(i in seq_along(extended)){
    use <- extended[seq_len(i)]; dd <- data.frame(y=ytr,xtr[,use,drop=FALSE]); fit <- tryCatch(glm(y~.,dd,family=binomial()),error=function(e)NULL)
    score <- if(is.null(fit)) rep(NA_real_,nrow(xte)) else tryCatch(as.numeric(predict(fit,newdata=xte[,use,drop=FALSE],type="response")),error=function(e)rep(NA_real_,nrow(xte)))
    roc <- tryCatch(pROC::roc(yte,score,quiet=TRUE,direction="<"),error=function(e)NULL); auc <- if(is.null(roc))NA_real_ else as.numeric(roc$auc)
    dp <- if(is.null(old_roc)||is.null(roc))NA_real_ else tryCatch(as.numeric(pROC::roc.test(old_roc,roc,paired=TRUE,method="delong")$p.value),error=function(e)NA_real_)
    gain <- if(is.null(old_roc)) TRUE else is.finite(dp)&&dp<alpha&&auc>as.numeric(old_roc$auc)
    if(gain){last_sig<-i;nonsig<-0L}else nonsig<-nonsig+1L
    rows[[i]]<-tibble::tibble(step=i,feature=extended[[i]],auc=auc,delta_auc=if(is.null(old_roc))NA_real_ else auc-as.numeric(old_roc$auc),delong_p=dp,accepted=gain)
    old_roc<-roc
    if(nonsig>=2L){stop_at<-i;break}
  }
  if(last_sig<1L)last_sig<-min(1L,length(extended))
  # The stopping rule defines the evaluated prefix as the parsimonious panel.
  # Returning only the last individually significant step collapses a sequence
  # such as one baseline feature plus two small gains to N=1, even though all
  # three features were required to establish the stopping point.
  panel_n <- if(is.finite(stop_at)) stop_at else last_sig
  list(parsimonious=extended[seq_len(panel_n)],extended=extended,log=dplyr::bind_rows(rows))
}

# Pradeep-style evidence-set prediction: outer 80/20 splits and 10-fold CV for lambda.
evidence_biom_lists <- function(outdir, biom_vars, layer, pred.prot.use = NULL, pred.met.use = NULL) {
  rr <- function(job) { f <- file.path(le8_job_dir(outdir, job), paste0(sub("_.*$","",job), ".res.rds")); if (file.exists(f)) readRDS(f) else list() }
  r1 <- rr("c1_correlate"); r2 <- rr("c2_cause")
  c1 <- r1[["association"]] %||% tibble::tibble()
  c2 <- r2[["MR"]] %||% tibble::tibble()
  pwas <- if (nrow(c1) && "term" %in% names(c1)) {
    q <- if ("FDR" %in% names(c1)) c1$FDR else p.adjust(c1$p.value, "BH")
    unique(as.character(c1$term[is.finite(q) & q < .05]))
  } else character()
  mr <- if (nrow(c2) && "exposure" %in% names(c2)) {
    q <- if ("FDR_all" %in% names(c2)) c2$FDR_all else p.adjust(c2$pval, "BH")
    unique(as.character(c2$exposure[is.finite(q) & q < .05]))
  } else character()
  user <- if (layer == "protein") pred.prot.use else pred.met.use
  z <- list(All = biom_vars, `PWAS/MWAS significant` = base::intersect(pwas, biom_vars), `MR significant` = base::intersect(mr, biom_vars))
  if (!is.null(user) && length(base::intersect(user, biom_vars))) z$`User specified` <- base::intersect(user, biom_vars)
  empty <- names(z)[lengths(z) == 0L]
  if (length(empty)) warning("No biomarkers found for: ", paste(empty, collapse=", "), call.=FALSE)
  z
}

split_re_prediction <- function(dat, biom_lists, basic, tvar, evar, method = "cox", test_frac = .20, inner = 10, seed = 2026, checkpoint_dir = NULL, force = FALSE) {
  set.seed(seed)
  test <- sample(seq_len(nrow(dat)),max(1L,round(nrow(dat)*test_frac)))
  dat$fold_id <- 0L; dat$fold_id[test] <- 1L
  configs <- list(`Model 0`=list(fx="RE",method=method,varX=basic,vars.basic=character(),opt=list(inner_cv_n=inner,lambda_rule="lambda.min",alpha=1)))
  for (sn in names(biom_lists)) if (length(biom_lists[[sn]])) {
    configs[[paste(sn,"Biomarkers",sep="||")]] <- list(fx="RE",method=method,varX=biom_lists[[sn]],vars.basic=character(),opt=list(inner_cv_n=inner,lambda_rule="lambda.1se",alpha=1))
    configs[[paste(sn,"Combined",sep="||")]] <- list(fx="RE",method=method,varX=biom_lists[[sn]],vars.basic=basic,opt=list(inner_cv_n=inner,lambda_rule="lambda.1se",alpha=1))
  }
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
  zz <- vector("list",length(configs))
  predict_single <- function(nm, var) {
    tr <- dat[dat$fold_id != 1L,,drop=FALSE]; te <- dat[dat$fold_id == 1L,,drop=FALSE]
    med <- median(suppressWarnings(as.numeric(tr[[var]])),na.rm=TRUE); if(!is.finite(med)) med <- 0
    xtr <- suppressWarnings(as.numeric(tr[[var]])); xte <- suppressWarnings(as.numeric(te[[var]]))
    xtr[!is.finite(xtr)] <- med; xte[!is.finite(xte)] <- med
    dtr <- data.frame(time=tr[[tvar]],event=tr[[evar]],x=xtr)
    ok <- is.finite(dtr$time)&!is.na(dtr$event)&is.finite(dtr$x)
    fit <- switch(method,
      binomial=tryCatch(glm(I(event==1)~x,data=dtr[ok,,drop=FALSE],family=binomial()),error=function(e)NULL),
      gaussian=tryCatch(lm(event~x,data=dtr[ok,,drop=FALSE]),error=function(e)NULL),
      cox=tryCatch(coxph(Surv(time,event)~x,data=dtr[ok,,drop=FALSE]),error=function(e)NULL))
    if(is.null(fit)) return(tibble::tibble())
    score <- tryCatch(as.numeric(predict(fit,newdata=data.frame(x=xte),type=if(method=="binomial")"response" else if(method=="cox")"lp" else "response")),error=function(e)rep(NA_real_,nrow(te)))
    ci <- if(method=="binomial") tryCatch(as.numeric(pROC::auc(pROC::roc(te[[evar]],score,quiet=TRUE))),error=function(e)NA_real_) else fcidx(Surv(te[[tvar]],te[[evar]]),score)
    data.frame(eid=te$eid,fold=1L,model=nm,fx="RE",method=method,time=te[[tvar]],event=te[[evar]],risk=score,cidx=ci,n_sel=1L)
  }
  for (i in seq_along(configs)) {
    nm <- names(configs)[i]; cf <- configs[[i]]
    tag <- gsub("[^A-Za-z0-9]+","_",nm)
    cov_tag <- paste0("c",length(basic),"_",sum(utf8ToInt(paste(sort(basic),collapse="|"))))
    cp <- if (is.null(checkpoint_dir)) NA_character_ else file.path(checkpoint_dir,sprintf("%02d_%s_n%d_%s_%s.rds",i,tag,length(cf$varX),cov_tag,method))
    if (!force && !is.na(cp) && file.exists(cp) && file.size(cp)>0) {
      message("Reuse prediction checkpoint: ",basename(cp)); zz[[i]] <- readRDS(cp)
    } else {
      zz[[i]] <- if(length(cf$varX)==1L && !length(cf$vars.basic)) predict_single(nm,cf$varX[[1]]) else
        pred_risk_by_fold(1L,dat,tvar,evar,basic,stats::setNames(list(cf),nm),id.var="eid")
      if (!is.na(cp) && nrow(zz[[i]])) saveRDS(zz[[i]],cp,compress=FALSE)
    }
    invisible(gc())
  }
  z <- dplyr::bind_rows(zz); rm(zz); invisible(gc())
  if (!nrow(z)) return(tibble::tibble())
  z |> tidyr::separate(model,c("biom_set","model"),sep="\\|\\|",fill="left",extra="merge") |>
    dplyr::mutate(biom_set=ifelse(model=="Model 0",NA_character_,biom_set)) |>
    dplyr::select(eid,biom_set,model,time,event,score=risk,cindex=cidx,n_selected=n_sel)
}

# Fit one final Pradeep/glmnet biomarker model in the eligible prediction
# cohort, then apply its fixed coefficients to every participant with an omics
# assay. These population scores are descriptive and must not be used for
# validation metrics; held-out predictions above remain the performance data.
pradeep_population_scores <- function(train, target, biom_lists, tvar, evar,
                                      method = "binomial", inner = 10, seed = 2026) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Prediction requires glmnet.", call. = FALSE)
  out <- vector("list", length(biom_lists)); names(out) <- names(biom_lists)
  for (i in seq_along(biom_lists)) {
    nm <- names(biom_lists)[i]; vv <- intersect(biom_lists[[i]], intersect(names(train), names(target)))
    vv <- vv[vapply(train[, vv, drop=FALSE], function(z) {
      z <- suppressWarnings(as.numeric(z)); is.finite(stats::var(z, na.rm=TRUE)) && stats::var(z, na.rm=TRUE) > 0
    }, logical(1))]
    if (!length(vv)) next
    xtr <- as.matrix(data.frame(lapply(train[, vv, drop=FALSE], function(z) suppressWarnings(as.numeric(z))), check.names=FALSE))
    xall <- as.matrix(data.frame(lapply(target[, vv, drop=FALSE], function(z) suppressWarnings(as.numeric(z))), check.names=FALSE))
    storage.mode(xtr) <- storage.mode(xall) <- "double"
    med <- apply(xtr, 2, function(z) { m <- median(z[is.finite(z)], na.rm=TRUE); ifelse(is.finite(m),m,0) })
    for (j in seq_along(med)) { xtr[!is.finite(xtr[,j]),j] <- med[j]; xall[!is.finite(xall[,j]),j] <- med[j] }
    y <- switch(method, gaussian=as.numeric(train[[evar]]), binomial=as.integer(train[[evar]]==1),
                cox=survival::Surv(train[[tvar]],train[[evar]]))
    ok <- is.finite(train[[tvar]]) & !is.na(train[[evar]])
    if (method == "binomial" && length(unique(y[ok])) < 2L) next
    set.seed(seed+i)
    nf <- if (method == "binomial") max(3L,min(inner,min(table(y[ok])))) else inner
    fit <- glmnet::cv.glmnet(xtr[ok,,drop=FALSE],y[ok],family=method,alpha=1,nfolds=nf,
                             type.measure=ifelse(method=="binomial","auc","deviance"))
    bb <- as.matrix(coef(fit,s="lambda.1se")); nsel <- sum(bb[-1,1]!=0)
    score <- as.numeric(predict(fit,newx=xall,s="lambda.1se",type=ifelse(method=="binomial","link","link")))
    out[[i]] <- tibble::tibble(eid=target$eid,biom_set=nm,model="Biomarkers",score=score,n_selected=nsel)
  }
  dplyr::bind_rows(out)
}

# Yu-fair comparison: publication-aligned LightGBM using the same held-out rows
# as split_re_prediction. Missing biomarker values are left for LightGBM's
# native missing-value handling, as in the reference implementation.
split_yu_prediction <- function(dat, biom_vars, basic, tvar, evar, horizon = 10,
                                test_frac = .20, seed = 2026,
                                checkpoint_file = NULL, force = FALSE) {
  if (!requireNamespace("lightgbm", quietly = TRUE))
    stop("Yu-fair prediction requires the R package 'lightgbm' (included in environment.yml).", call. = FALSE)
  if (!force && !is.null(checkpoint_file) && file.exists(checkpoint_file) && file.size(checkpoint_file) > 0)
    return(readRDS(checkpoint_file))
  biom_vars <- intersect(biom_vars, names(dat)); basic <- intersect(basic, names(dat))
  if (!length(biom_vars)) return(tibble::tibble())
  set.seed(seed); test_id <- sample(seq_len(nrow(dat)), max(1L, round(nrow(dat) * test_frac)))
  tr_id <- setdiff(seq_len(nrow(dat)), test_id)
  known <- (dat[[evar]][tr_id] == 1 & dat[[tvar]][tr_id] <= horizon) | dat[[tvar]][tr_id] > horizon
  tr_id <- tr_id[which(known)]; y <- as.integer(dat[[evar]][tr_id] == 1 & dat[[tvar]][tr_id] <= horizon)
  if (length(unique(y)) < 2L) stop("Yu-fair training set has fewer than two horizon-status classes.", call. = FALSE)

  numeric_matrix <- function(rows, vars) {
    x <- as.matrix(data.frame(lapply(dat[rows, vars, drop = FALSE], function(z) suppressWarnings(as.numeric(z))), check.names = FALSE))
    storage.mode(x) <- "double"; x
  }
  xb_tr <- numeric_matrix(tr_id, biom_vars); xb_te <- numeric_matrix(test_id, biom_vars)
  clinical_matrices <- function() {
    if (!length(basic)) return(list(tr = matrix(nrow=length(tr_id),ncol=0), te = matrix(nrow=length(test_id),ncol=0)))
    # C5 clinical predictors are numeric scores/covariates. Construct the two
    # matrices separately so missing values can never make model.matrix drop
    # participants and desynchronise rows from tr_id/test_id.
    a <- numeric_matrix(tr_id,basic); b <- numeric_matrix(test_id,basic)
    med <- apply(a,2,function(z){v<-median(z[is.finite(z)],na.rm=TRUE);ifelse(is.finite(v),v,0)})
    for(j in seq_along(med)){a[!is.finite(a[,j]),j]<-med[j];b[!is.finite(b[,j]),j]<-med[j]}
    if(nrow(a)!=length(tr_id)||nrow(b)!=length(test_id))stop("Yu-fair clinical matrix row contract failed.",call.=FALSE)
    list(tr=a,te=b)
  }
  xc <- clinical_matrices()
  params <- list(objective="binary", metric="auc", max_depth=15L, num_leaves=10L,
                 bagging_fraction=.70, bagging_freq=1L, learning_rate=.01,
                 feature_fraction=.70, max_bin=63L, verbosity=-1L,
                 num_threads=max(1L,as.integer(Sys.getenv("N_CORES",unset="4"))), seed=as.integer(seed))
  fit_one <- function(xtr, xte, model_name) {
    ds <- lightgbm::lgb.Dataset(data=xtr,label=y)
    fit <- lightgbm::lgb.train(params=params,data=ds,nrounds=500L,verbose=-1)
    score <- as.numeric(predict(fit,xte)); imp <- tryCatch(lightgbm::lgb.importance(fit),error=function(e)NULL)
    nsel <- if(is.null(imp)) ncol(xtr) else sum(imp$Gain > 0,na.rm=TRUE)
    tibble::tibble(eid=dat$eid[test_id],biom_set="Yu / LightGBM",model=model_name,
                   time=dat[[tvar]][test_id],event=dat[[evar]][test_id],score=score,
                   cindex=fcidx(Surv(dat[[tvar]][test_id],dat[[evar]][test_id]),score),n_selected=nsel)
  }
  ans <- dplyr::bind_rows(fit_one(xb_tr,xb_te,"Biomarkers"),
                           fit_one(cbind(xb_tr,xc$tr),cbind(xb_te,xc$te),"Combined"))
  if (!is.null(checkpoint_file)) { dir.create(dirname(checkpoint_file),recursive=TRUE,showWarnings=FALSE); saveRDS(ans,checkpoint_file,compress=FALSE) }
  ans
}

plot_evidence_prediction <- function(pred, horizon = 10) {
  if (!nrow(pred)) return(blank_plot("Evidence-based prediction", "No out-of-fold prediction was available"))
  sets <- unique(stats::na.omit(pred$biom_set)); clinical <- pred |> dplyr::filter(model=="Model 0")
  pred <- dplyr::bind_rows(pred |> dplyr::filter(model!="Model 0"), lapply(sets, function(s) clinical |> dplyr::mutate(biom_set=s))) |>
    dplyr::mutate(biom_set=factor(biom_set,levels=sets)) |>
    dplyr::group_by(biom_set, model) |> dplyr::mutate(z=as.numeric(scale(score)), quintile=dplyr::ntile(score,5), decile=dplyr::ntile(score,10)) |> dplyr::ungroup()
  pal <- c(`Model 0`="#333333", Biomarkers="#A6DDA0", Combined="#32B43C")
  selected_n <- pred |> dplyr::filter(model=="Biomarkers") |> dplyr::group_by(biom_set) |>
    dplyr::summarise(n_biom=as.integer(stats::median(n_selected,na.rm=TRUE)),.groups="drop")
  biom_legend <- paste0("Biom (N=",paste(selected_n$n_biom[match(sets,as.character(selected_n$biom_set))],collapse="/"),")")
  pd <- pred |> dplyr::filter(model=="Biomarkers") |> dplyr::mutate(status=factor(event,0:1,c("Controls","Cases"))) |>
    dplyr::group_by(biom_set) |> dplyr::mutate(threshold=quantile(z[event==0],.95,na.rm=TRUE),DR=mean(z[event==1]>=dplyr::first(threshold),na.rm=TRUE)) |> dplyr::ungroup()
  ann <- pd |> dplyr::distinct(biom_set,threshold,DR)
  a <- ggplot2::ggplot(pd,ggplot2::aes(z,fill=status))+ggplot2::geom_density(alpha=.48)+ggplot2::geom_vline(data=ann,ggplot2::aes(xintercept=threshold))+ggplot2::geom_text(data=ann,ggplot2::aes(x=Inf,y=Inf,label=sprintf("DR = %.1f%%\nFPR = 5.0%%",100*DR)),hjust=1.05,vjust=1.2,inherit.aes=FALSE,size=2.7)+ggplot2::scale_fill_manual(values=c(Controls="#F8E8D2",Cases="#C89D6B"))+ggplot2::labs(x="Biomarker score (s.d.)",y="Density",fill=NULL)+theme_5c(9)+ggplot2::theme(legend.position="top",legend.justification="center",legend.box.just="center")
  cum <- pred |> dplyr::group_by(biom_set,model,quintile) |> dplyr::group_modify(~{sf<-survival::survfit(survival::Surv(time,event)~1,data=.x); tibble::tibble(time=sf$time,cuminc=1-sf$surv)}) |> dplyr::ungroup()
  b <- ggplot2::ggplot(cum |> dplyr::filter(model=="Biomarkers"),ggplot2::aes(time,cuminc,color=factor(quintile)))+ggplot2::geom_step()+ggplot2::coord_cartesian(xlim=c(0,horizon))+ggplot2::scale_color_brewer(palette="YlOrBr",direction=1)+ggplot2::labs(x="Follow-up time (years)",y="Cumulative incidence",color="Score quintile")+theme_5c(9)+ggplot2::theme(legend.position="top",legend.justification="center",legend.box.just="center")
  rate <- pred |> dplyr::filter(model=="Biomarkers") |> dplyr::group_by(biom_set,decile) |> dplyr::summarise(rate=1000*sum(event)/sum(time),.groups="drop") |> dplyr::filter(rate>0)
  c <- ggplot2::ggplot(rate,ggplot2::aes(decile*10-5,rate,color=decile))+ggplot2::geom_point(size=2)+ggplot2::scale_color_gradient(low="#FFD29B",high="#A65E00",guide="none")+ggplot2::scale_y_log10()+ggplot2::labs(x="Biomarker score percentile",y="Incidence rate per 1,000 years")+theme_5c(9)
  roc <- pred |> dplyr::group_by(biom_set,model) |> dplyr::group_modify(~{r<-pROC::roc(.x$event,.x$score,quiet=TRUE); tibble::tibble(fpr=1-r$specificities,tpr=r$sensitivities,auc=as.numeric(r$auc))}) |> dplyr::ungroup()
  labs <- roc |> dplyr::group_by(biom_set,model) |> dplyr::summarise(auc=dplyr::first(auc),.groups="drop") |>
    dplyr::left_join(selected_n,by="biom_set") |>
    dplyr::mutate(model_label=dplyr::case_when(model=="Model 0"~"Clinical",model=="Biomarkers"~paste0("Biom (N=",n_biom,")"),TRUE~"Combined"),label=sprintf("%s: %.3f",model_label,auc),x=.98,y=c(.12,.06,.18)[match(model,c("Model 0","Biomarkers","Combined"))])
  d <- ggplot2::ggplot(roc,ggplot2::aes(fpr,tpr,color=model))+ggplot2::geom_abline(slope=1,intercept=0,linetype=3,color="grey80")+ggplot2::geom_step()+ggplot2::geom_text(data=labs,ggplot2::aes(x=x,y=y,label=label,color=model),hjust=1,inherit.aes=FALSE,size=2.5)+ggplot2::scale_color_manual(values=pal,breaks=c("Biomarkers","Combined","Model 0"),labels=c(biom_legend,"Combined","Clinical"))+ggplot2::labs(x="FPR",y="True positive rate",color=NULL)+theme_5c(9)+ggplot2::theme(legend.position="top",legend.justification="center",legend.box.just="center")
  (a/b/c/d)+patchwork::plot_layout(guides="keep")+patchwork::plot_annotation(title="Prediction by biomarker evidence set: 80% training / 20% testing",theme=ggplot2::theme(plot.title=ggplot2::element_text(hjust=.5))) & ggplot2::facet_wrap(~biom_set,nrow=1,scales="free")
}
