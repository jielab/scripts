# Explicit assay budgets, training-only feature selection and held-out evaluation.
# No GDF15 or NT-proBNP penalty is hard coded. Ablations are separately refitted.
le8_prepare_prediction_matrix <- function(train,test,vars) {
  train<-as.data.frame(train);test<-as.data.frame(test)
  vars<-intersect(vars,intersect(names(train),names(test)));a<-b<-list();audit<-list()
  for(v in vars) {
    x<-train[[v]];y<-test[[v]]
    if(is.numeric(x)||is.integer(x)||is.logical(x)) {
      x<-as.numeric(x);y<-as.numeric(y);m<-median(x[is.finite(x)],na.rm=TRUE)
      if(!is.finite(m))next
      x[!is.finite(x)]<-m;y[!is.finite(y)]<-m;s<-sd(x)
      if(!is.finite(s)||s<=0)next
      a[[v]]<-(x-m)/s;b[[v]]<-(y-m)/s
      audit[[length(audit)+1L]]<-tibble(variable=v,type="numeric",center=m,scale=s,reference="training median",unseen_test_levels=0L)
    }else {
      x<-as.character(x);y<-as.character(y);x[x==""]<-NA;y[y==""]<-NA
      tab<-table(x);if(!length(tab))next
      mode<-names(tab)[which.max(tab)];lev<-sort(names(tab));x[is.na(x)]<-mode
      unseen<-!is.na(y)&!y%in%lev;y[is.na(y)|unseen]<-mode
      if(length(lev)<2L)next
      # Reference and levels are learned exclusively from the training data.
      for(l in setdiff(lev,mode)) {
        n<-paste0(v,"__",l);a[[n]]<-as.numeric(x==l);b[[n]]<-as.numeric(y==l)
      }
      audit[[length(audit)+1L]]<-tibble(variable=v,type="categorical",center=NA_real_,scale=NA_real_,reference=mode,unseen_test_levels=sum(unseen))
    }
  }
  if(!length(a))return(list(train=matrix(nrow=nrow(train),ncol=0),test=matrix(nrow=nrow(test),ncol=0),audit=bind_rows(audit)))
  xx<-as.matrix(as.data.frame(a,check.names=FALSE));yy<-as.matrix(as.data.frame(b,check.names=FALSE))
  # Remove exactly redundant columns based on training predictors, not outcomes.
  qr0<-qr(xx);keep<-sort(qr0$pivot[seq_len(qr0$rank)])
  list(train=xx[,keep,drop=FALSE],test=yy[,keep,drop=FALSE],audit=bind_rows(audit))
}
le8_fit_budget_model <- function(train,test,clinical,features,tvar,evar) {
  vars<-unique(c(clinical,features));xx<-le8_prepare_prediction_matrix(train,test,vars)
  if(!ncol(xx$train))return(list(status="no usable predictors"))
  orig<-colnames(xx$train);safe<-paste0("x",seq_len(ncol(xx$train)))
  tr<-as.data.frame(xx$train);te<-as.data.frame(xx$test);names(tr)<-names(te)<-safe
  tr$.time<-train[[tvar]];tr$.event<-train[[evar]]
  f<-reformulate(safe,response="survival::Surv(.time,.event)")
  fit<-tryCatch(survival::coxph(f,tr,ties="efron",x=TRUE,y=TRUE,model=TRUE,singular.ok=FALSE),error=function(e)e)
  if(inherits(fit,"condition"))return(list(status=conditionMessage(fit)))
  if(any(!is.finite(coef(fit))))return(list(status="non-finite Cox coefficient"))
  lp<-drop(as.matrix(te)%*%coef(fit));bh<-survival::basehaz(fit,centered=FALSE)
  count<-sum(vapply(features,function(v)any(orig==v|startsWith(orig,paste0(v,"__"))),logical(1)))
  list(status="ok",fit=fit,lp=lp,baseline_hazard=bh,N_selected=count,
    coefficient=tibble(variable=orig,beta=as.numeric(coef(fit))),preprocess=xx$audit)
}
le8_risk_at <- function(obj,horizon) {
  bh<-obj$baseline_hazard;idx<-which(bh$time<=horizon)
  H<-if(length(idx))bh$hazard[max(idx)]else 0
  pmin(1-1e-8,pmax(1e-8,-expm1(-H*exp(pmin(30,pmax(-30,obj$lp))))))
}
le8_ipcw <- function(time,event,horizon) {
  sf<-survival::survfit(survival::Surv(time,1-event)~1)
  at<-function(t,left=FALSE) {
    i<-findInterval(t,sf$time)
    if(left){hit<-i>0L&sf$time[pmax(i,1L)]==t;i[hit]<-i[hit]-1L}
    z<-rep(1,length(t));z[i>0]<-sf$surv[i[i>0]];z
  }
  y<-as.numeric(event==1&time<=horizon);w<-numeric(length(time))
  case<-y==1;ctrl<-time>horizon
  gcase<-at(time[case],TRUE);gt<-at(horizon)
  valid<-length(gt)==1L&&is.finite(gt)&&gt>.01&&sum(case)>=20&&sum(ctrl)>=20
  if(valid){w[case]<-1/pmax(gcase,.01);w[ctrl]<-1/gt}
  list(y=y,w=w,status=if(valid)"ok"else"insufficient supported follow-up/cases/controls",
    N_case=sum(case),N_control=sum(ctrl),G_horizon=gt,
    assumption="marginal independent censoring; death is censored, so predicted risk is not a competing-risk CIF")
}
le8_weighted_auc <- function(p,y,w) {
  ii<-is.finite(p)&is.finite(y)&is.finite(w)&w>0
  p<-p[ii];y<-y[ii];w<-w[ii];o<-order(p);p<-p[o];y<-y[o];w<-w[o]
  if(!length(p)||sum(w*y)<=0||sum(w*(1-y))<=0)return(NA_real_)
  g<-cumsum(c(TRUE,diff(p)!=0));a<-as.numeric(rowsum(w*y,g,reorder=FALSE))
  b<-as.numeric(rowsum(w*(1-y),g,reorder=FALSE));cum<-cumsum(b)-b
  sum(a*(cum+.5*b))/(sum(a)*sum(b))
}
le8_evaluate_risk <- function(time,event,p,horizon,model,budget,paradigm,ablation="none",B=0,seed=2026) {
  iw<-le8_ipcw(time,event,horizon)
  base<-tibble(model,budget,paradigm,ablation,horizon,N=length(time),events_by_horizon=iw$N_case,
    controls_at_horizon=iw$N_control,status=iw$status,AUC=NA_real_,AUC_lo=NA_real_,AUC_hi=NA_real_,
    Brier=NA_real_,Brier_lo=NA_real_,Brier_hi=NA_real_,calibration_intercept=NA_real_,calibration_slope=NA_real_,
    censoring_survival=iw$G_horizon,estimand=iw$assumption)
  if(iw$status!="ok")return(list(metrics=base,calibration=tibble(),decision=tibble(),bootstrap=tibble()))
  y<-iw$y;w<-iw$w;n<-length(y);base$AUC<-le8_weighted_auc(p,y,w);base$Brier<-sum(w*(y-p)^2)/n
  lp<-qlogis(p);ok<-w>0
  fit<-tryCatch(glm(y[ok]~lp[ok],family=quasibinomial(),weights=w[ok]),error=function(e)NULL)
  ic<-tryCatch(glm(y[ok]~1+offset(lp[ok]),family=quasibinomial(),weights=w[ok]),error=function(e)NULL)
  if(!is.null(fit))base$calibration_slope<-unname(coef(fit)[2])
  if(!is.null(ic))base$calibration_intercept<-unname(coef(ic)[1])
  group<-pmin(10L,ceiling(rank(p,ties.method="first")/n*10L))
  cal<-tibble(group,p,y,w)|>group_by(group)|>summarise(N=n(),predicted=mean(p),
    observed_ipcw=sum(w*y)/sum(w),.groups="drop")|>mutate(model,horizon,budget,paradigm,ablation)
  dec<-map_dfr(seq(.01,.20,by=.01),function(th){tr<-p>=th
    tibble(model,horizon,threshold=th,net_benefit=(sum(w*y*tr)-sum(w*(1-y)*tr)*th/(1-th))/n,
      treat_all=(sum(w*y)-sum(w*(1-y))*th/(1-th))/n,treat_none=0)})
  boot<-tibble()
  if(B>0L) {
    set.seed(seed)
    boot<-map_dfr(seq_len(B),function(i){ix<-sample.int(n,n,replace=TRUE)
      z<-le8_ipcw(time[ix],event[ix],horizon)
      tibble(replicate=i,AUC=if(z$status=="ok")le8_weighted_auc(p[ix],z$y,z$w)else NA_real_,
        Brier=if(z$status=="ok")mean(z$w*(z$y-p[ix])^2)else NA_real_)
    })|>mutate(model,horizon)
    for(metric in c("AUC","Brier")) {
      v<-boot[[metric]];v<-v[is.finite(v)]
      if(length(v)>=max(20,ceiling(B*.8))){base[[paste0(metric,"_lo")]]<-quantile(v,.025,names=FALSE);base[[paste0(metric,"_hi")]]<-quantile(v,.975,names=FALSE)}
    }
  }
  list(metrics=base,calibration=cal,decision=dec,bootstrap=boot)
}
le8_c5_additions <- function(dat,biom_vars,ranked,training_screen,distalset,geneticset,
                             clinical,tvar,evar,split,layer,outdir) {
  budgets<-sort(unique(as.integer(le8_csv_env("C5_ASSAY_BUDGETS","3,5,10"))))
  if(!length(budgets)||anyNA(budgets)||any(budgets<1L|budgets>30L))stop("C5_ASSAY_BUDGETS must contain integers 1..30")
  horizons<-as.numeric(le8_csv_env("C5_REVIEW_HORIZONS","2,5,10"))
  if(any(!is.finite(horizons)|horizons<=0))stop("Invalid C5_REVIEW_HORIZONS")
  B<-as.integer(le8_num_env("C5_REVIEW_BOOT",100));if(B<0L)stop("C5_REVIEW_BOOT must be nonnegative")
  cv<-le8_csv_env("C5_REVIEW_CLINICAL_VARS",paste(clinical,collapse=","))
  missing<-setdiff(cv,names(dat));if(length(missing))stop("Specified clinical predictors are missing: ",paste(missing,collapse=","))
  train<-dat[split=="training",,drop=FALSE];test<-dat[split=="validation",,drop=FALSE]
  if(nrow(train)<500||nrow(test)<100||sum(train[[evar]])<40)return(list(status=tibble(status="insufficient training/validation data")))
  # Near-term candidate ranking uses only training participants; early censoring
  # is retained by a time-truncated Cox model, not converted into non-cases.
  near<-train;near$.near_time<-pmin(near[[tvar]],2);near$.near_event<-as.integer(near[[evar]]==1&near[[tvar]]<=2)
  near_scan<-cox_scan(near,head(ranked,as.integer(le8_num_env("C5_NEAR_SCREEN_MAX",100))),cv,Y,
    time_var=".near_time",event_var=".near_event")
  near_rank<-near_scan|>filter(is.finite(p.value))|>arrange(p.value)|>pull(term)
  paradigms<-list(`Overall-selected`=ranked,`Distal-selected`=ranked[ranked%in%distalset],`Near-term-selected`=near_rank)
  if(truthy(Sys.getenv("C5_GENETIC_EVIDENCE_INDEPENDENT",unset="FALSE")))
    paradigms[["Genetic-region-selected"]]<-ranked[ranked%in%geneticset]
  designs<-list(list(name="Clinical",features=character(),budget=0L,paradigm="Clinical",ablation="none"))
  for(pa in names(paradigms))for(k in budgets)if(length(paradigms[[pa]])>=k)
    designs[[length(designs)+1L]]<-list(name=paste(pa,k,sep="_"),features=head(paradigms[[pa]],k),budget=k,paradigm=pa,ablation="none")
  # Fixed-panel ablation: do not fill the removed slot with a substitute assay.
  # This answers incremental contribution of the removed marker(s).
  full<-head(ranked,max(budgets))
  if(layer=="protein")for(nm in c("GDF15","NTPROBNP","GDF15+NTPROBNP")) {
    drop<-strsplit(nm,"+",fixed=TRUE)[[1]]
    designs[[length(designs)+1L]]<-list(name=paste0("Overall_minus_",nm),features=setdiff(full,drop),
      budget=length(full),paradigm="Overall-selected",ablation=nm)
  }
  rd<-le8_job_dir(outdir,"c5_consolidate");dir.create(file.path(rd,"review_models"),showWarnings=FALSE)
  metrics<-cal<-dca<-boot<-members<-prep<-coefs<-list();pred_cache<-list()
  for(i in seq_along(designs)) {
    ds<-designs[[i]];obj<-le8_fit_budget_model(train,test,cv,ds$features,tvar,evar)
    members[[i]]<-tibble(model=ds$name,feature=if(length(ds$features))ds$features else NA_character_,
      requested_budget=ds$budget,actual_panel_size=length(ds$features),ablation=ds$ablation,status=obj$status)
    if(obj$status!="ok")next
    prep[[i]]<-obj$preprocess|>mutate(model=ds$name)
    coefs[[i]]<-obj$coefficient|>mutate(model=ds$name)
    saveRDS(obj,file.path(rd,"review_models",paste0(gsub("[^A-Za-z0-9_.-]","_",ds$name),".rds")))
    ev<-lapply(horizons,function(h) {
      p<-le8_risk_at(obj,h)
      # Identical seed at a given horizon gives paired validation resamples.
      z<-le8_evaluate_risk(test[[tvar]],test[[evar]],p,h,ds$name,ds$budget,ds$paradigm,ds$ablation,
        B=B,seed=SEED+as.integer(h*100))
      z$metrics$actual_n_assays<-length(ds$features);z$metrics$fitted_n_assays<-obj$N_selected
      z
    })
    metrics[[i]]<-bind_rows(lapply(ev,`[[`,"metrics"));cal[[i]]<-bind_rows(lapply(ev,`[[`,"calibration"))
    dca[[i]]<-bind_rows(lapply(ev,`[[`,"decision"));boot[[i]]<-bind_rows(lapply(ev,`[[`,"bootstrap"))
  }
  met<-bind_rows(metrics);ca<-bind_rows(cal);dc<-bind_rows(dca);bo<-bind_rows(boot)
  delta<-tibble()
  if(nrow(bo)) {
    ref<-bo|>filter(model=="Clinical")|>select(replicate,horizon,AUC_ref=AUC,Brier_ref=Brier)
    delta<-bo|>filter(model!="Clinical")|>inner_join(ref,by=c("replicate","horizon"))|>
      mutate(delta_AUC=AUC-AUC_ref,delta_Brier=Brier-Brier_ref)|>group_by(model,horizon)|>
      summarise(valid=sum(is.finite(delta_AUC)),delta_AUC_lo=quantile(delta_AUC,.025,na.rm=TRUE),
        delta_AUC_hi=quantile(delta_AUC,.975,na.rm=TRUE),delta_Brier_lo=quantile(delta_Brier,.025,na.rm=TRUE),
        delta_Brier_hi=quantile(delta_Brier,.975,na.rm=TRUE),.groups="drop")
  }
  tables<-list(budget_metrics=met,panel_members=bind_rows(members),calibration=ca,decision_curves=dc,
    paired_delta_CI=delta,model_coefficients=bind_rows(coefs),training_preprocessing=bind_rows(prep),near_training_scan=near_scan,
    design=tibble(item=c("Clinical comparator","Genetic evidence used in prediction","Uncertainty","Prediction risk"),
      value=c(paste(cv,collapse=";"),Sys.getenv("C5_GENETIC_EVIDENCE_INDEPENDENT",unset="FALSE"),
        "Validation bootstrap of frozen fits; not uncertainty of the entire training-selection procedure",
        "Cause-specific/net risk with deaths censored; not a competing-risk cumulative incidence")))
  for(nm in names(tables))write_raw_csv(tables[[nm]],paste0("c5.review_",nm,".csv"),rd)
  # The row split is saved only in the local raw directory, not the workbook.
  data.table::fwrite(tibble(eid=dat$eid,split=split),file.path(rd,"c5.review_split.csv.gz"),compress="gzip")
  if(nrow(met)) {
    h<-max(horizons);mm<-met|>filter(horizon==h,status=="ok")
    a<-mm|>filter(ablation=="none")|>ggplot(aes(actual_n_assays,AUC,color=paradigm))+
      geom_line()+geom_point()+geom_errorbar(aes(ymin=AUC_lo,ymax=AUC_hi),width=.2)+
      labs(title=paste0("a. Assay budget and ",h,"-year discrimination"),x="Actual measured assays",y="IPCW AUC",color=NULL)+theme_5c(9)
    b<-ca|>filter(horizon==h,model%in%c("Clinical",paste0("Overall-selected_",max(budgets))))|>
      ggplot(aes(predicted,observed_ipcw,color=model))+geom_abline(slope=1,intercept=0,linetype=2)+geom_line()+geom_point()+
      labs(title="b. Frozen-model calibration",x="Predicted risk",y="Observed IPCW risk",color=NULL)+theme_5c(9)
    c<-mm|>filter(ablation!="none"|model==paste0("Overall-selected_",max(budgets)))|>
      ggplot(aes(AUC,model))+geom_errorbarh(aes(xmin=AUC_lo,xmax=AUC_hi),height=.15)+geom_point()+
      labs(title="c. Fixed-panel ablations, refitted in training",x="IPCW AUC",y=NULL)+theme_5c(9)
    d<-dc|>filter(horizon==h,model%in%c("Clinical",paste0("Overall-selected_",max(budgets))))|>
      ggplot(aes(threshold,net_benefit,color=model))+geom_line()+geom_hline(yintercept=0,linetype=2)+
      labs(title="d. Exploratory decision curves",x="Risk threshold",y="Net benefit",color=NULL)+theme_5c(9)
    save_plot((a|b)/(c|d)+plot_annotation(caption="Same held-out participants. No validation outcomes are used to select features, tune budgets or recalibrate predictions. Genetic-region selection requires independent evidence declaration."),
      "c5.Fig13.budget_ablation_calibration.png",18,12,outdir=outdir)
  }
  tables
}

# Relearn connection-based feature membership inside the outer training split.
# No validation LE8 measurements or outcomes enter this procedure.
le8_training_connection_set <- function(train,features,components,covars,rawdir) {
  components<-intersect(components,names(train));covars<-intersect(covars,names(train))
  if(length(components)<2)return(character())
  d<-as.data.frame(train);d<-d[complete.cases(d[,unique(c(components,covars)),drop=FALSE]),,drop=FALSE]
  set.seed(SEED+411)
  cap<-as.integer(le8_num_env("C5_CONNECTION_MAX_N",60000))
  if(nrow(d)>cap)d<-d[sample.int(nrow(d),cap),,drop=FALSE]
  if(nrow(d)<1000)return(character())
  half<-sample(rep(1:2,length.out=nrow(d)));features<-intersect(features,names(d))
  scans<-map_dfr(1:2,function(h) {
    dd<-d[half==h,,drop=FALSE]
    map_dfr(components,function(cmp){
      cv<-unique(c(covars,setdiff(components,cmp)))
      M<-model.matrix(reformulate(cv),dd);q<-qr(M)
      yr<-qr.resid(q,as.numeric(scale(dd[[cmp]])));sy<-sqrt(sum(yr^2));df<-nrow(dd)-q$rank-1L
      if(!is.finite(sy)||sy<=0||df<10)return(tibble())
      blocks<-split(features,ceiling(seq_along(features)/64))
      map_dfr(blocks,function(bb){
        x<-as.matrix(dd[,bb,drop=FALSE]);storage.mode(x)<-"double"
        for(j in seq_len(ncol(x))){z<-x[,j];m<-median(z[is.finite(z)],na.rm=TRUE)
          if(!is.finite(m))m<-0;z[!is.finite(z)]<-m;s<-sd(z)
          x[,j]<-if(is.finite(s)&&s>0)(z-mean(z))/s else 0}
        xr<-qr.resid(q,x);den<-sqrt(colSums(xr^2))*sy
        r<-as.numeric(crossprod(xr,yr))/pmax(den,1e-20);r<-pmax(-.999999,pmin(.999999,r))
        z<-r*sqrt(df/(1-r*r));p<-2*pt(abs(z),df,lower.tail=FALSE)
        tibble(feature=bb,component=cmp,r,z,p,N=nrow(dd),half=h)
      })
    })
  })
  if(!nrow(scans))return(character())
  scans<-scans|>group_by(half,component)|>mutate(FDR=p.adjust(p,"BH"))|>ungroup()
  a<-scans|>filter(half==1)|>group_by(feature)|>mutate(specificity=abs(z)/sum(abs(z)))|>
    slice_max(abs(z),n=1,with_ties=FALSE)|>ungroup()|>select(feature,component,r1=r,FDR1=FDR,specificity)
  b<-scans|>filter(half==2)|>select(feature,component,r2=r,FDR2=FDR)
  joined<-left_join(a,b,by=c("feature","component"))|>
    mutate(selected=is.finite(FDR1)&is.finite(FDR2)&FDR1<.05&FDR2<.05&
      sign(r1)==sign(r2)&is.finite(specificity)&specificity>=.35,
      scope="outer training participants only; replicated LE8 association, not intervention evidence")
  write_raw_csv(joined,"c5.connection_membership_training_only.csv",rawdir)
  joined$feature[joined$selected%in%TRUE]
}
