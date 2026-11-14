# Transparent prognostic refinement of CKM stages. All fitting, preprocessing,
# molecular screening and reference percentiles are learned inside training folds.

ckm_make_folds <- function(d,cfg) {
  k<-as.integer(cfg$clock$folds %||% 5L)
  if(k<2L)stop("clock.folds must be >=2")
  group_col<-cfg$clock$group_column %||% "center"
  if(nzchar(cfg$paths$groups %||% "")) {
    gr<-ckm_table(cfg$paths$groups,"CV groups")
    if(!"group" %in% names(gr))stop("groups file needs eid,group")
    g<-as.character(gr$group[match(d$eid,gr$eid)])
  } else if(group_col=="eid")g<-as.character(d$eid)
  else {
    if(!group_col %in% names(d))stop("Missing CV grouping column ",group_col,"; set group_column='eid' explicitly for random participant CV")
    g<-as.character(d[[group_col]])
  }
  if(anyNA(g)||any(!nzchar(g)))stop("Every participant needs a nonmissing CV group")
  z<-data.table::data.table(group=g)
  sizes<-z[,.(n=.N),by=group]
  if(nrow(sizes)<k)stop("Fewer CV groups than folds")
  set.seed(cfg$seed);sizes[,tie:=runif(.N)];data.table::setorder(sizes,-n,tie)
  load<-rep(0,k);assign<-integer(nrow(sizes))
  for(i in seq_len(nrow(sizes))) {f<-which.min(load);assign[i]<-f;load[f]<-load[f]+sizes$n[i]}
  data.table::data.table(eid=as.character(d$eid),group=g,fold=assign[match(g,sizes$group)])
}

ckm_apply_transform <- function(x,kind) {
  if(length(kind)!=1L||is.na(kind)||!kind %in% c("none","log1p","log2"))
    stop("Clock transform must be none, log1p, or log2")
  x<-suppressWarnings(as.numeric(x));x[!is.finite(x)]<-NA_real_
  if(kind=="log1p") {x[x<0]<-NA_real_;x<-log1p(x)}
  if(kind=="log2") {x[x<=0]<-NA_real_;x<-log2(x)}
  x
}

ckm_design_fit <- function(d,features,transforms,cfg) {
  spec<-list();parts<-list();dropped<-character()
  for(v in c("stage","sex")) {
    lev<-if(v=="stage")as.character(0:3) else c("F","M")
    observed<-unique(as.character(d[[v]]));observed<-lev[lev %in% observed]
    if(length(observed)<2L)stop("Training set has <2 ",v," categories")
    spec[[v]]<-list(kind="factor",levels=observed,columns=paste0(v,"_",seq_len(length(observed)-1L)))
    parts[[v]]<-matrix(vapply(observed[-1L],function(z)as.numeric(as.character(d[[v]])==z),numeric(nrow(d))),
                       nrow=nrow(d),ncol=length(observed)-1L)
    colnames(parts[[v]])<-spec[[v]]$columns
  }
  for(v in unique(c("age",features))) {
    if(!v %in% names(d)) {dropped<-c(dropped,paste(v,"missing_column"));next}
    kind<-transforms[[v]] %||% "none"
    a<-ckm_apply_transform(d[[v]],kind)
    min_distinct<-if(v=="age")3L else 2L
    if(mean(is.na(a))>cfg$clock$max_missing||length(unique(a[is.finite(a)]))<min_distinct) {
      if(v=="age")stop("Insufficient age information")
      dropped<-c(dropped,paste(v,"missing_or_constant"));next
    }
    med<-median(a,na.rm=TRUE);limits<-quantile(a,c(.005,.995),na.rm=TRUE,names=FALSE)
    missing<-is.na(a);a[missing]<-med;a<-pmin(limits[2],pmax(limits[1],a))
    # A prespecified spline for age avoids a purely chronological linear clock.
    if(v=="age") {
      b<-splines::ns(a,df=3);s<-list(kind="spline",transform=kind,median=med,limits=limits,
        knots=attr(b,"knots"),boundary=attr(b,"Boundary.knots"),missing=FALSE)
    } else {
      sc<-sd(a);if(!is.finite(sc)||sc<1e-10){dropped<-c(dropped,paste(v,"constant"));next}
      b<-matrix((a-mean(a))/sc,ncol=1);s<-list(kind="numeric",transform=kind,median=med,limits=limits,mean=mean(a),sd=sc,missing=any(missing))
      if(s$missing)b<-cbind(b,as.numeric(missing))
    }
    s$columns<-paste0("v",length(spec)+1L,"_",seq_len(ncol(b)))
    colnames(b)<-s$columns;spec[[v]]<-s;parts[[v]]<-b
  }
  x<-do.call(cbind,parts);good<-apply(x,2,function(z)sd(z)>1e-10)
  x<-x[,good,drop=FALSE]
  qr0<-qr(x,tol=1e-8);keep<-sort(qr0$pivot[seq_len(qr0$rank)])
  removed<-setdiff(colnames(x),colnames(x)[keep]);x<-x[,keep,drop=FALSE]
  list(x=x,spec=spec,columns=colnames(x),dropped=c(dropped,removed))
}

ckm_design_apply <- function(d,prep) {
  parts<-list()
  for(v in names(prep$spec)) {
    s<-prep$spec[[v]]
    if(s$kind=="factor") {
      a<-as.character(d[[v]])
      b<-matrix(vapply(s$levels[-1L],function(z)as.numeric(a==z),numeric(nrow(d))),
                nrow=nrow(d),ncol=length(s$levels)-1L)
      b[!a %in% s$levels,]<-NA_real_
    } else {
      a<-ckm_apply_transform(if(v %in% names(d))d[[v]] else rep(NA_real_,nrow(d)),s$transform)
      missing<-is.na(a);a[missing]<-s$median;a<-pmin(s$limits[2],pmax(s$limits[1],a))
      if(s$kind=="spline")b<-splines::ns(a,knots=s$knots,Boundary.knots=s$boundary)
      else {b<-matrix((a-s$mean)/s$sd,ncol=1);if(s$missing)b<-cbind(b,as.numeric(missing))}
    }
    colnames(b)<-s$columns;parts[[v]]<-b
  }
  do.call(cbind,parts)[,prep$columns,drop=FALSE]
}

ckm_fit_cs <- function(d,x,cause,min_events=20L) {
  ev<-as.integer(d$status==cause)
  if(sum(ev)<min_events)stop("Insufficient cause ",cause," events: ",sum(ev))
  dat<-data.frame(.time=d$time,.event=ev,x,check.names=FALSE)
  f<-stats::as.formula(paste("survival::Surv(.time,.event) ~",paste(colnames(x),collapse=" + ")))
  warns<-character()
  fit<-withCallingHandlers(survival::coxph(f,data=dat,ties="breslow",x=TRUE,y=TRUE,model=TRUE,
      singular.ok=FALSE,control=survival::coxph.control(iter.max=50)),
      warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
  if(any(!is.finite(coef(fit)))||any(!is.finite(diag(vcov(fit))))||any(grepl("infinite|converg",warns,ignore.case=TRUE)))
    stop("Non-estimable Cox fit: ",paste(warns,collapse="; "))
  fit
}

ckm_baseline <- function(fit,horizon,steps) {
  grid<-seq(0,horizon,length.out=steps+1L)
  if(is.null(fit))return(list(grid=grid,increment=rep(0,steps)))
  b<-survival::basehaz(fit,centered=FALSE)
  j<-findInterval(grid,b$time);haz<-c(0,b$hazard)[j+1L]
  list(grid=grid,increment=pmax(0,diff(haz)))
}

ckm_cif <- function(lp1,lp2,b1,b2) {
  ok<-is.finite(lp1)&is.finite(lp2)
  n<-length(lp1);risk<-surv<-free<-rep(NA_real_,n)
  r1<-exp(pmin(30,pmax(-30,lp1[ok])));r2<-exp(pmin(30,pmax(-30,lp2[ok])))
  ss<-rep(1,length(r1));ff<-rep(0,length(r1));rr<-rep(0,length(r1))
  for(i in seq_along(b1$increment)) {
    h1<-b1$increment[i]*r1;h2<-b2$increment[i]*r2;total<-h1+h2
    inc<- -expm1(-total);portion<-ifelse(total>0,h1/total,0)
    rr<-rr+ss*inc*portion
    dt<-b1$grid[i+1L]-b1$grid[i]
    ff<-ff+ss*ifelse(total>1e-12,inc/total,1)*dt
    ss<-ss*exp(-total)
  }
  risk[ok]<-rr;surv[ok]<-ss;free[ok]<-ff
  list(risk=risk,event_free_years=free,survival=surv)
}

ckm_fit_bundle <- function(d,features,transforms,cfg) {
  if(sum(d$status==1L)<cfg$clock$min_events)stop("Too few incident events for clock fitting")
  h<-cfg$clock$horizon
  if(sum(d$time>=h)<cfg$clock$min_at_horizon)stop("Insufficient observed follow-up at the configured horizon")
  prep<-ckm_design_fit(d,features,transforms,cfg)
  if(sum(d$status==1L)<cfg$clock$events_per_parameter*ncol(prep$x))
    stop("Incident events do not cover the prespecified parameter budget")
  fit1<-ckm_fit_cs(d,prep$x,1L,cfg$clock$min_events)
  # Death is modelled separately, rather than treated as independent censoring for absolute risk.
  dp<-ckm_design_fit(d,intersect(c("bmi","egfr"),features),transforms,cfg)
  death_n<-sum(d$status==2L)
  if(death_n==0L)fit2<-NULL
  else {
    if(death_n<cfg$clock$events_per_parameter*ncol(dp$x))dp<-ckm_design_fit(d,character(),list(),cfg)
    if(death_n<max(10L,ncol(dp$x)*2L))stop("Too few deaths to estimate competing mortality reliably")
    fit2<-ckm_fit_cs(d,dp$x,2L,10L)
  }
  out<-list(prep=prep,death_prep=dp,event=fit1,death=fit2,
      b1=ckm_baseline(fit1,h,cfg$clock$integration_steps),b2=ckm_baseline(fit2,h,cfg$clock$integration_steps),
      horizon=h,features=features,death_events=death_n)
  out
}

ckm_predict_bundle <- function(fit,d) {
  x<-ckm_design_apply(d,fit$prep);xd<-ckm_design_apply(d,fit$death_prep)
  lp1<-as.numeric(x %*% coef(fit$event));lp2<-if(is.null(fit$death))rep(0,nrow(d)) else as.numeric(xd %*% coef(fit$death))
  z<-ckm_cif(lp1,lp2,fit$b1,fit$b2);z$lp<-lp1;z
}

ckm_screen <- function(d,domain,clinical_fit,cfg) {
  man<-domain$manifest
  man<-man[man$type %in% c("continuous","binary") & man$column %in% names(domain$data),,drop=FALSE]
  if(!nrow(man))return(man)
  residual<-as.numeric(residuals(clinical_fit$event,type="martingale"))
  score<-rep(NA_real_,nrow(man));vectors<-vector("list",nrow(man));ii<-match(d$eid,domain$data$eid)
  for(j in seq_len(nrow(man))) {
    a<-ckm_apply_transform(domain$data[[man$column[j]]][ii],man$transform[j])
    ok<-is.finite(a)
    if(mean(!ok)>cfg$clock$max_missing||sum(ok)<30L||sd(a[ok])<1e-10)next
    a[!ok]<-median(a[ok]);a<-pmax(quantile(a,.005),pmin(quantile(a,.995),a))
    a<-as.numeric(scale(a));score[j]<-abs(cor(a,residual));vectors[[j]]<-a
  }
  order0<-order(-score,man$feature,na.last=NA)
  selected<-integer()
  max_k<-min(cfg$clock$max_omics,max(0L,floor(sum(d$status==1L)/cfg$clock$events_per_parameter)-ncol(clinical_fit$prep$x)))
  if(max_k<1L)return(man[FALSE,,drop=FALSE])
  for(j in order0) {
    if(!length(selected)||all(vapply(selected,function(k)abs(cor(vectors[[j]],vectors[[k]]))<cfg$clock$max_correlation,logical(1))))selected<-c(selected,j)
    if(length(selected)>=max_k)break
  }
  ans<-man[selected,,drop=FALSE];ans$screen_score<-score[selected];ans
}

ckm_attach_selected <- function(d,domain,sel) {
  if(!nrow(sel))return(d)
  idx<-match(d$eid,domain$data$eid)
  for(i in seq_len(nrow(sel)))d[[paste0("omics_",i)]]<-domain$data[[sel$column[i]]][idx]
  d
}

ckm_censor_weights <- function(train,test,h) {
  w<-rep(NA_real_,nrow(test));y<-as.integer(test$status==1L & test$time<=h)
  for(st in sort(unique(test$stage))) {
    tr<-train[train$stage==st,,drop=FALSE];te<-which(test$stage==st)
    if(nrow(tr)<10L)next
    sf<-survival::survfit(survival::Surv(time,status==0L)~1,data=tr)
    gs<-c(1,as.numeric(sf$surv));gt<-sf$time
    # A participant censored exactly at h has an observed horizon outcome.
    # Use G(h-) consistently, so a common administrative cutoff at h does not
    # make every weight unavailable merely because G(h) drops to zero there.
    gh<-gs[findInterval(h,gt,left.open=TRUE)+1L]
    if(!is.finite(gh)||gh<.05)next
    tt<-test$time[te];ss<-test$status[te]
    early<-tt<=h & ss!=0L;complete<-tt>=h
    wt<-rep(0,length(te));left<-findInterval(tt,gt,left.open=TRUE)
    wt[early]<-1/gs[left[early]+1L];wt[complete & !early]<-1/gh
    w[te]<-wt
  }
  list(weight=w,target=y)
}

ckm_weighted_auc <- function(p,y,w) {
  ok<-is.finite(p)&is.finite(w)&w>0
  if(!any(ok))return(NA_real_)
  z<-data.table::data.table(p=p[ok],case=w[ok]*y[ok],control=w[ok]*(1-y[ok]))
  z<-z[,.(case=sum(case),control=sum(control)),by=p][order(p)]
  if(sum(z$case)<=0||sum(z$control)<=0)return(NA_real_)
  sum(z$case*(cumsum(z$control)-.5*z$control))/(sum(z$case)*sum(z$control))
}

ckm_metrics <- function(z) {
  ok<-is.finite(z$risk)&is.finite(z$ipcw)
  z<-z[ok,,drop=FALSE]
  if(!nrow(z))return(list(n=0L,events=0L,auc=NA_real_,brier=NA_real_,cindex=NA_real_,ipcw_mass=NA_real_,status="unsupported_horizon"))
  cc<-tryCatch(survival::concordance(survival::Surv(time,status==1L)~risk,data=z,reverse=TRUE)$concordance,error=function(e)NA_real_)
  list(n=nrow(z),events=sum(z$status==1L),auc=ckm_weighted_auc(z$risk,z$target,z$ipcw),
    brier=sum(z$ipcw*(z$target-z$risk)^2)/nrow(z),cindex=cc,ipcw_mass=sum(z$ipcw)/nrow(z),status="ok")
}

ckm_compare <- function(pred,cfg) {
  # Paired comparisons use exactly the same held-out participants and censor weights.
  res<-list();counter<-0L
  for(co in unique(pred$cohort)) {
    available<-unique(pred$model[pred$cohort==co])
    pairs<-lapply(setdiff(available,"stage"),function(mod)list(model=mod,reference="stage"))
    if(all(c("clinical","molecular") %in% available))pairs<-c(pairs,list(list(model="molecular",reference="clinical")))
    for(pair in pairs)for(st in c("all",as.character(0:3))) {
    mod<-pair$model;reference<-pair$reference
    q<-pred[pred$cohort==co & (st=="all"|as.character(pred$stage)==st),,drop=FALSE]
    a<-q[q$model==reference,];b<-q[q$model==mod,];ii<-match(a$eid,b$eid)
    good<-!is.na(ii)&is.finite(a$risk)&is.finite(a$ipcw)
    a<-a[good,,drop=FALSE];b<-b[ii[good],,drop=FALSE]
    good<-is.finite(b$risk)&is.finite(b$ipcw);a<-a[good,,drop=FALSE];b<-b[good,,drop=FALSE]
    if(nrow(a)<30L)next
    ma<-ckm_metrics(a);mb<-ckm_metrics(b);da<-mb$auc-ma$auc;db<-mb$brier-ma$brier
    bs<-matrix(NA_real_,nrow=cfg$clock$bootstrap,ncol=2)
    if(nrow(bs)>0L) {
      set.seed(cfg$seed+counter);groups<-unique(a$group)
      for(j in seq_len(nrow(bs))) {
        sampled<-sample(groups,length(groups),replace=TRUE);counts<-tabulate(match(sampled,groups),nbins=length(groups))
        mult<-counts[match(a$group,groups)];w<-a$ipcw*mult;nn<-sum(mult)
        if(nn>0L) {
          bs[j,1]<-ckm_weighted_auc(b$risk,b$target,w)-ckm_weighted_auc(a$risk,a$target,w)
          bs[j,2]<-sum(w*((b$target-b$risk)^2-(a$target-a$risk)^2))/nn
        }
      }
    }
    ci<-function(v)if(sum(is.finite(v))>=10L)quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE) else c(NA_real_,NA_real_)
    ca<-ci(bs[,1]);cb<-ci(bs[,2]);counter<-counter+1L
    res[[counter]]<-data.table::data.table(cohort=co,model=mod,reference=reference,stage=st,n=nrow(a),horizon=cfg$clock$horizon,
      delta_auc=da,delta_auc_lower=ca[1],delta_auc_upper=ca[2],delta_brier=db,delta_brier_lower=cb[1],delta_brier_upper=cb[2],
      bootstrap=cfg$clock$bootstrap,uncertainty="paired CV-group bootstrap, conditional on fitted models")
    }
  }
  data.table::rbindlist(res,fill=TRUE)
}

ckm_run_clock <- function(d,cfg,domains=list()) {
  d<-as.data.frame(d);d$eid<-as.character(d$eid)
  elig<-d$stage %in% 0:3 & is.finite(d$time)&d$time>0 & d$status %in% 0:2 & is.finite(d$age)&!is.na(d$sex)
  if("riskset_eligible" %in% names(d))elig<-elig & d$riskset_eligible %in% TRUE
  z<-d[elig,,drop=FALSE]
  if(nrow(z)<100L)stop("Insufficient eligible baseline stage 0-3 participants")
  folds<-ckm_make_folds(z,cfg);z$fold<-folds$fold;z$group<-folds$group
  populations<-list(clinical=NULL)
  for(nm in intersect(c("prot","met"),names(domains)))populations[[nm]]<-domains[[nm]]
  if(all(c("prot","met") %in% names(domains)) && isTRUE(cfg$clock$joint)) {
    p<-data.table::copy(domains$prot);m<-data.table::copy(domains$met)
    # Prefix assay columns to prevent protein/metabolite name collisions.
    for(nm in c("p","m")) {
      ob<-get(nm);old<-setdiff(names(ob$data),"eid");new<-paste0(nm,"__",old)
      data.table::setnames(ob$data,old,new);ob$manifest$column<-paste0(nm,"__",ob$manifest$column)
      assign(nm,ob)
    }
    populations$joint<-list(data=merge(p$data,m$data,by="eid",all=FALSE),manifest=data.table::rbindlist(list(p$manifest,m$manifest),fill=TRUE))
  }
  predictions<-list();fits<-list();selection<-list();fit_status<-list();nr<-ns<-0L
  clinical_features<-cfg$clock$clinical_features
  transforms<-as.list(cfg$clock$transforms %||% list())
  for(co in names(populations)) {
    domain<-populations[[co]]
    zz<-if(is.null(domain))z else z[z$eid %in% domain$data$eid,,drop=FALSE]
    models<-if(is.null(domain))c("stage","clinical") else c("stage","clinical","molecular")
    if(nrow(zz)<100L) {fit_status[[length(fit_status)+1L]]<-data.table::data.table(cohort=co,fold=NA_integer_,model="all",status="insufficient_cohort",n=nrow(zz));next}
    for(f in sort(unique(z$fold))) {
      tr<-zz[zz$fold!=f,,drop=FALSE];te<-zz[zz$fold==f,,drop=FALSE]
      if(!nrow(te))next
      cw<-ckm_censor_weights(tr,te,cfg$clock$horizon)
      clin<-NULL
      for(mod in models) {
        message("Clock: ",co," / ",mod," / held-out fold ",f)
        sel<-data.frame();features<-if(mod=="stage")character() else clinical_features
        fit_result<-tryCatch({
          train<-tr;test<-te;tf<-transforms
          if(mod=="molecular") {
            if(is.null(clin))stop("Clinical training fit failed; molecular screen unavailable")
            sel<-ckm_screen(tr,domain,clin,cfg)
            if(!nrow(sel))stop("No molecular features pass training-only screen and parameter budget")
            train<-ckm_attach_selected(tr,domain,sel);test<-ckm_attach_selected(te,domain,sel)
            features<-c(clinical_features,paste0("omics_",seq_len(nrow(sel))))
            tf<-c(transforms,stats::setNames(as.list(sel$transform),paste0("omics_",seq_len(nrow(sel)))))
          }
          fit<-ckm_fit_bundle(train,features,tf,cfg)
          if(mod=="clinical")clin<-fit
          pp<-ckm_predict_bundle(fit,test);rp<-ckm_predict_bundle(fit,train)
          percentile<-position<-accel<-rep(NA_real_,nrow(test))
          for(s in 0:3) {
            rr<-rp$risk[train$stage==s];rr<-rr[is.finite(rr)]
            ii<-which(test$stage==s)
            if(length(rr)<20L||!length(ii))next
            # Fixed training reference. Test ranks and outcomes never set these cutpoints.
            percentile[ii]<-stats::ecdf(rr)(pp$risk[ii])
            percentile[ii]<-pmin(1-1/(2*length(rr)),pmax(1/(2*length(rr)),percentile[ii]))
            position[ii]<-s+percentile[ii]
          }
          # Age/sex/stage residual is descriptive, not a true longitudinal speed.
          # This auxiliary description must not discard otherwise valid risk
          # predictions when a held-out fold contains an unseen category.
          accel<-tryCatch({
            axis<-function(p)log(-log1p(-pmin(1-1e-12,pmax(1e-12,p))))
            tr_axis<-axis(rp$risk)
            ref<-stats::lm(tr_axis~splines::ns(age,3)+sex+factor(stage),data=train)
            scale0<-sd(residuals(ref),na.rm=TRUE);expected<-rep(NA_real_,nrow(test))
            supported<-is.finite(pp$risk) & test$stage %in% train$stage & as.character(test$sex) %in% as.character(train$sex)
            if(any(supported))expected[supported]<-predict(ref,newdata=test[supported,,drop=FALSE])
            if(is.finite(scale0)&&scale0>1e-8)(axis(pp$risk)-expected)/scale0 else rep(NA_real_,nrow(test))
          },error=function(e)rep(NA_real_,nrow(test)))
          nr<-nr+1L
          predictions[[nr]]<-data.table::data.table(eid=test$eid,group=test$group,fold=f,cohort=co,model=mod,stage=test$stage,
            time=test$time,status=test$status,horizon=cfg$clock$horizon,risk=pp$risk,percentile=percentile,
            position=position,grade=ifelse(is.finite(percentile),paste0("S",test$stage,".G",pmin(3L,1L+floor(percentile*3))),NA_character_),
            age_adjusted_residual=accel,event_free_years=pp$event_free_years,ipcw=cw$weight,target=cw$target)
          if(nrow(sel)){ns<-ns+1L;selection[[ns]]<-data.table::as.data.table(sel);selection[[ns]][,`:=`(cohort=co,fold=f)]}
          list(status=if(all(is.finite(pp$risk)))"ok" else "partial: unsupported held-out categories",
               n_train=nrow(train),n_test=nrow(test),n_supported=sum(is.finite(pp$risk)),
               n_unsupported=sum(!is.finite(pp$risk)),events=sum(train$status==1L),n_parameters=length(coef(fit$event)),
               n_selected=nrow(sel),omitted_predictors=paste(fit$prep$dropped,collapse=";"))
        },error=function(e)list(status=paste0("failed: ",conditionMessage(e)),n_train=nrow(tr),n_test=nrow(te),events=sum(tr$status==1L)))
        fit_status[[length(fit_status)+1L]]<-data.table::as.data.table(c(list(cohort=co,fold=f,model=mod),fit_result))
      }
    }
    # Frozen deployment model is separate from OOF predictions, and never evaluated in-sample.
    final_clin<-NULL
    for(mod in models) {
      final<-tryCatch({
        features<-if(mod=="stage")character() else clinical_features;train<-zz;tf<-transforms;sel<-data.frame()
        if(mod=="molecular") {
          if(is.null(final_clin))stop("No full clinical fit")
          sel<-ckm_screen(zz,domain,final_clin,cfg);if(!nrow(sel))stop("No molecular features")
          train<-ckm_attach_selected(zz,domain,sel);features<-c(clinical_features,paste0("omics_",seq_len(nrow(sel))))
          tf<-c(transforms,stats::setNames(as.list(sel$transform),paste0("omics_",seq_len(nrow(sel)))))
        }
        ff<-ckm_fit_bundle(train,features,tf,cfg);if(mod=="clinical")final_clin<-ff
        rr<-ckm_predict_bundle(ff,train)$risk
        ff$reference<-lapply(0:3,function(s)sort(rr[train$stage==s]));names(ff$reference)<-as.character(0:3)
        ff$selection<-sel;ff$training_eids<-zz$eid;ff
      },error=function(e)list(error=conditionMessage(e)))
      fits[[paste(co,mod,sep=".")]]<-final
    }
  }
  pred<-data.table::rbindlist(predictions,fill=TRUE)
  statuses<-data.table::rbindlist(fit_status,fill=TRUE)
  if(!nrow(pred))stop("Every clock fit failed: ",paste(unique(statuses$status),collapse="; "))
  metrics<-calibration<-list();mm<-cc<-0L
  for(co in unique(pred$cohort))for(mod in unique(pred$model[pred$cohort==co]))for(st in c("all",as.character(0:3))) {
    q<-as.data.frame(pred[pred$cohort==co & pred$model==mod & (st=="all"|as.character(pred$stage)==st),])
    mm<-mm+1L;metrics[[mm]]<-data.table::as.data.table(c(list(cohort=co,model=mod,stage=st,horizon=cfg$clock$horizon),ckm_metrics(q)))
    q<-q[is.finite(q$risk)&is.finite(q$ipcw),,drop=FALSE]
    if(nrow(q)<30L)next
    # Quantile bins are descriptive validation bins, never used for model fitting.
    q$bin<-pmin(10L,ceiling(rank(q$risk,ties.method="average")/nrow(q)*10))
    cal<-data.table::as.data.table(q)[,.(n=.N,predicted=mean(risk),observed=sum(ipcw*target)/.N),by=bin]
    cal[,`:=`(cohort=co,model=mod,stage=st,horizon=cfg$clock$horizon)]
    cc<-cc+1L;calibration[[cc]]<-cal
  }
  distribution<-pred[cohort=="clinical" & model=="clinical" & is.finite(percentile),.(n=.N),by=.(stage,bin=pmin(9L,floor(percentile*10)))]
  list(predictions=pred,metrics=data.table::rbindlist(metrics,fill=TRUE),calibration=data.table::rbindlist(calibration,fill=TRUE),
    distribution=distribution,comparison=ckm_compare(as.data.frame(pred),cfg),selection=data.table::rbindlist(selection,fill=TRUE),
    fit_status=statuses,models=fits,folds=folds,
    notes=c("Position = observed stage + percentile against that stage's training reference; it is not a validated clinical substage or elapsed time.",
      "Absolute risk is cumulative incidence from two cause-specific Cox models, including competing death; numerical integration uses a fixed fine time grid.",
      "Molecular selection uses training-only clinical martingale residual screening, correlation pruning and a fixed feature budget; held-out outcomes never choose features.",
      "AUC/Brier use stage-specific training censor weights. Harrell C treats competing death as censoring and has a different estimand.",
      "Model comparisons and uncertainty use the same held-out IDs; bootstrap intervals are conditional on fitted models.",
      "The model named stage adjusts for age and sex as well as categorical stage; molecular models are also compared with clinical models in the same assay cohort.",
      "Unsupported held-out stage/sex categories receive missing predictions and are counted in fit_status; metrics describe only evaluable predictions with censoring support.",
      "Age-adjusted residuals and restricted event-free years must not be interpreted as stage-to-stage speed or years."))
}
