# C5: unified score phenotyping and prediction figures.

suppressPackageStartupMessages({
  fdir<-Sys.getenv("LE8_FDIR",unset=file.path(Sys.getenv("DIRSCRIPT"),"f"))
  source(file.path(fdir,"comm.f.R"));source(file.path(fdir,"c5_pred.R"))
})
LE8_JOB <- "c5_consolidate"
C5_CODE_VERSION <- "2026-09-05.5c-audit-v1"
C5_MODEL_VERSION <- "2026-09-05.5c-audit-v1"
OUTER <- as.integer(Sys.getenv("C5_OUTER_FOLDS",unset="5"))
INNER <- as.integer(Sys.getenv("C5_INNER_FOLDS",unset="10"))
HORIZON <- as.numeric(Sys.getenv("C5_HORIZON",unset="10"))
C5_NESTED_CV <- truthy(Sys.getenv("C5_NESTED_CV",unset="FALSE"))
C5_TEST_FRAC <- as.numeric(Sys.getenv("C5_TEST_FRAC",unset="0.20"))
C5_LGB_FOLDS <- as.integer(Sys.getenv("C5_LGB_FOLDS",unset="5"))
C5_LGB_NROUND <- as.integer(Sys.getenv("C5_LGB_NROUND",unset="500"))
C5_YU_PRESELECT_N <- as.integer(Sys.getenv("C5_YU_PRESELECT_N",unset="257"))
C5_YU_LIST <- Sys.getenv("C5_YU_LIST",unset="")
C5_EXTENDED_N <- as.integer(Sys.getenv("C5_EXTENDED_N",unset="30"))
C5_YY_BINS <- as.integer(Sys.getenv("C5_YY_BINS",unset="18"))
C5_TEMPORAL_STEP <- as.numeric(Sys.getenv("C5_TEMPORAL_STEP",unset="1"))
C5_TIME_CAP <- as.numeric(Sys.getenv("C5_TIME_CAP",unset="16"))
C5_BOOT <- as.integer(Sys.getenv("C5_AUC_BOOT",unset="150"))
C5_MIN_BIN_N <- as.integer(Sys.getenv("C5_MIN_BIN_N",unset="15"))
C5_TOPN_MAX <- as.integer(Sys.getenv("C5_TOPN_MAX",unset="30"))
C5_TOPN_BOOT <- as.integer(Sys.getenv("C5_TOPN_BOOT",unset="150"))
C5_LEAD_MIN_CASES <- as.integer(Sys.getenv("C5_LEAD_MIN_CASES",unset="100"))
C5_MECHANISM_N <- as.integer(Sys.getenv("C5_MECHANISM_N",unset="10"))
C5_BOOT_MAX_CASES <- as.integer(Sys.getenv("C5_BOOT_MAX_CASES",unset="1500"))
C5_BOOT_MAX_CONTROLS <- as.integer(Sys.getenv("C5_BOOT_MAX_CONTROLS",unset="6000"))
C5_PAIR_MAX_CASES <- as.integer(Sys.getenv("C5_PAIR_MAX_CASES",unset="1500"))
C5_DISTAL_LANDMARK <- as.numeric(Sys.getenv("C5_DISTAL_LANDMARK",unset="5"))
C5_DISTAL_TOP <- as.integer(Sys.getenv("C5_DISTAL_TOP",unset="500"))
C5_INC_PROT <- Sys.getenv("C5_INC_PROT",unset="")
C5_INC_MET <- Sys.getenv("C5_INC_MET",unset="")
C5_LEAD_ANCHORS <- unique(trimws(strsplit(Sys.getenv("C1_DIRECTION_ANCHORS",
  unset="PCSK9,LPA,GDF15,NTPROBNP,MMP12"),",",fixed=TRUE)[[1]]))

# Publication typography for C5 only.
theme_5c <- function(base_size=12){
  theme_classic(base_size=base_size)+theme(
    plot.title=element_text(face="bold",size=base_size*1.12,hjust=0),
    plot.subtitle=element_text(face="bold",size=base_size*.92,color="grey30"),
    axis.title=element_text(face="bold",size=base_size*1.03),axis.text=element_text(face="bold",size=base_size,color="black"),
    legend.title=element_text(face="bold"),legend.text=element_text(face="bold"),strip.background=element_blank(),strip.text=element_text(face="bold"),
    panel.grid.major.y=element_line(color="grey91",linewidth=.25),panel.grid.minor=element_blank(),plot.margin=margin(7,10,7,10))
}


# 🚩 Fixed outer split and score fitting
make_outer_split <- function(dat,evar=NULL,test_frac=C5_TEST_FRAC,seed=SEED){
  set.seed(seed)
  if(!is.null(evar)&&evar%in%names(dat)){
    strata<-split(seq_len(nrow(dat)),as.character(dat[[evar]]),drop=TRUE)
    val<-unlist(lapply(strata,function(ii)sample(ii,max(1L,round(length(ii)*test_frac)))),use.names=FALSE)
  }else val<-sample(seq_len(nrow(dat)),max(1L,round(nrow(dat)*test_frac)))
  ifelse(seq_len(nrow(dat))%in%val,"validation","training")
}

fit_glmnet_score <- function(dat,vars,tvar,evar,split,label,inner=INNER,alpha=1,
                             lambda_rule=c("lambda.1se","lambda.min"),penalty_factor=NULL){
  lambda_rule<-match.arg(lambda_rule)
  if(!requireNamespace("glmnet",quietly=TRUE))stop("C5 requires glmnet.",call.=FALSE)
  vars<-intersect(unique(vars),names(dat));trid<-which(split=="training");vaid<-which(split=="validation")
  if(!length(vars)||length(trid)<500||length(vaid)<100)return(NULL)
  tr<-dat[trid,,drop=FALSE];va<-dat[vaid,,drop=FALSE];x<-impute_train_test(tr,va,vars)
  if(ncol(x$tr)<1)return(NULL)
  ok<-is.finite(tr[[tvar]])&!is.na(tr[[evar]])&tr[[tvar]]>0
  y<-Surv(tr[[tvar]][ok],tr[[evar]][ok]);events<-as.integer(tr[[evar]][ok]==1)
  if(sum(ok)<300||sum(events)<30)return(NULL)
  k<-min(inner,max(3L,min(table(events))))
  foldid<-make_folds(tr,evar,k,SEED+17)
  pf<-setNames(rep(1,ncol(x$tr)),x$vars)
  if(!is.null(penalty_factor)){
    supplied<-suppressWarnings(as.numeric(penalty_factor));names(supplied)<-names(penalty_factor)
    if(!is.null(names(penalty_factor))){
      hit<-intersect(names(penalty_factor),names(pf));pf[hit]<-supplied[match(hit,names(penalty_factor))]
    }else if(length(supplied)==length(pf))pf[]<-supplied
  }
  pf[!is.finite(pf)|pf<=0]<-1
  fit<-tryCatch(glmnet::cv.glmnet(x$tr[ok,,drop=FALSE],y,family="cox",alpha=alpha,foldid=foldid[ok],
    type.measure="C",standardize=FALSE,keep=TRUE,penalty.factor=unname(pf)),error=function(e)NULL)
  if(is.null(fit))return(NULL)
  lam<-if(lambda_rule=="lambda.min")fit$lambda.min else fit$lambda.1se
  li<-which.min(abs(fit$lambda-lam));oof<-rep(NA_real_,nrow(tr))
  if(!is.null(fit$fit.preval)&&length(li))oof[ok]<-as.numeric(fit$fit.preval[,li])
  val<-tryCatch(as.numeric(predict(fit,x$te,s=lambda_rule,type="link")),error=function(e)rep(NA_real_,nrow(va)))
  b<-as.matrix(coef(fit,s=lambda_rule));sel<-setdiff(rownames(b)[b[,1]!=0],"(Intercept)")
  rows<-bind_rows(
    tibble(eid=tr$eid,row_id=trid,split="training",time=tr[[tvar]],event=tr[[evar]],score=oof),
    tibble(eid=va$eid,row_id=vaid,split="validation",time=va[[tvar]],event=va[[evar]],score=val))|>
    mutate(method=label,n_features=length(x$vars),n_selected=length(sel),engine="glmnet_cox",lambda_rule=lambda_rule)
  list(rows=rows,selected=sel,available=x$vars,fit=fit,engine="glmnet_cox",lambda_rule=lambda_rule,
       penalty_factor=pf,preprocess=list(vars=x$vars,center=x$center,scale=x$scale))
}

fit_lightgbm_score <- function(dat,vars,tvar,evar,split,label="Yu-style / LightGBM",k=C5_LGB_FOLDS,nrounds=C5_LGB_NROUND){
  if(!requireNamespace("lightgbm",quietly=TRUE)){
    warning("R package 'lightgbm' is not installed; the Yu-style LightGBM score will be skipped.",call.=FALSE);return(NULL)
  }
  vars<-intersect(unique(vars),names(dat));trid<-which(split=="training");vaid<-which(split=="validation")
  if(!length(vars)||length(trid)<500||length(vaid)<100)return(NULL)
  tr<-dat[trid,,drop=FALSE];va<-dat[vaid,,drop=FALSE]
  # Yu-style LightGBM keeps missing values natively. Only remove invariant columns.
  keep<-vapply(tr[,vars,drop=FALSE],function(z){z<-suppressWarnings(as.numeric(z));sum(is.finite(z))>=100&&is.finite(sd(z,na.rm=TRUE))&&sd(z,na.rm=TRUE)>0},logical(1))
  vars<-vars[keep];if(length(vars)<2)return(NULL)
  to_mat<-function(d){x<-as.matrix(data.frame(lapply(d[,vars,drop=FALSE],function(z)suppressWarnings(as.numeric(z))),check.names=FALSE));storage.mode(x)<-"double";x}
  xtr<-to_mat(tr);xva<-to_mat(va)
  eligible<-is.finite(tr[[tvar]])&!is.na(tr[[evar]])&tr[[tvar]]>0&
    ((tr[[evar]]==1&tr[[tvar]]<=HORIZON)|tr[[tvar]]>HORIZON)
  y<-as.integer(tr[[evar]][eligible]==1&tr[[tvar]][eligible]<=HORIZON);if(sum(eligible)<300||length(unique(y))<2)return(NULL)
  nfold<-min(k,max(3L,min(table(y))));set.seed(SEED+23);foldid<-integer(length(y))
  for(value in 0:1){ii<-which(y==value);foldid[ii]<-sample(rep(seq_len(nfold),length.out=length(ii)))}
  oof<-rep(NA_real_,nrow(tr))
  params<-list(objective="binary",metric="auc",max_depth=15L,num_leaves=10L,
    bagging_fraction=.70,bagging_freq=1L,learning_rate=.01,feature_fraction=.70,max_bin=63L,
    lambda_l2=1,verbosity=-1L,num_threads=max(1L,N_CORES),seed=SEED)
  for(fd in sort(unique(foldid))){
    it<-which(foldid!=fd);iv<-which(foldid==fd);if(length(unique(y[it]))<2)next
    ds<-lightgbm::lgb.Dataset(data=xtr[eligible,,drop=FALSE][it,,drop=FALSE],label=y[it])
    ff<-tryCatch(lightgbm::lgb.train(params=params,data=ds,nrounds=nrounds,verbose=-1),error=function(e)NULL)
    if(!is.null(ff))oof[which(eligible)[iv]]<-as.numeric(predict(ff,xtr[eligible,,drop=FALSE][iv,,drop=FALSE]))
  }
  ds<-lightgbm::lgb.Dataset(data=xtr[eligible,,drop=FALSE],label=y);fit<-tryCatch(lightgbm::lgb.train(params=params,data=ds,nrounds=nrounds,verbose=-1),error=function(e)NULL)
  if(is.null(fit))return(NULL);val<-as.numeric(predict(fit,xva))
  imp<-tryCatch(lightgbm::lgb.importance(fit),error=function(e)NULL);sel<-if(is.null(imp))vars else as.character(imp$Feature[imp$Gain>0])
  rows<-bind_rows(
    tibble(eid=tr$eid,row_id=trid,split="training",time=tr[[tvar]],event=tr[[evar]],score=oof),
    tibble(eid=va$eid,row_id=vaid,split="validation",time=va[[tvar]],event=va[[evar]],score=val))|>
    mutate(method=label,n_features=length(vars),n_selected=length(sel),engine="lightgbm",lambda_rule=NA_character_)
  list(rows=rows,selected=sel,available=vars,fit=fit,engine="lightgbm",horizon=HORIZON,preprocess=list(vars=vars))
}

empty_score_rows <- function(){
  tibble(eid=numeric(),row_id=integer(),split=character(),time=numeric(),event=numeric(),
    score=numeric(),method=character(),n_features=integer(),n_selected=integer(),engine=character(),
    lambda_rule=character(),score_z=numeric(),score_native=numeric(),score_scale=character())
}

empty_prediction_rows <- function(){
  tibble(eid=numeric(),row_id=integer(),time=numeric(),event=numeric(),score=numeric(),
    biom_set=character(),model=character(),n_selected=integer(),cindex=numeric())
}

standardize_scores <- function(rows){
  if(is.null(rows)||!nrow(rows)||!all(c("method","score","split")%in%names(rows)))return(empty_score_rows())
  rows|>group_by(method)|>group_modify(function(d,key){
    m<-mean(d$score[d$split=="training"],na.rm=TRUE);s<-sd(d$score[d$split=="training"],na.rm=TRUE);if(!is.finite(s)||s==0)s<-1
    # A Cox linear predictor is not a probability; plogis(lp) was therefore a
    # misleading "native 0-1" scale. Preserve raw model output and use the
    # training-standardized score for cross-model temporal displays.
    d|>mutate(score_z=(score-m)/s,score_native=score,
      score_scale=ifelse(engine=="lightgbm","binary model output","Cox linear predictor"))
  })|>ungroup()
}

predict_prevalent_rows <- function(obj,prevalent,method,bvar,training_rows){
  if(is.null(obj)||!nrow(prevalent)||is.null(obj$fit))return(tibble())
  expected<-obj$preprocess$vars%||%obj$available
  missing<-setdiff(expected,names(prevalent));if(length(missing)){
    warning(method,": prevalent scoring missing ",length(missing)," model variables; no partial-column prediction was attempted.",call.=FALSE)
    return(tibble())
  }
  vars<-expected;if(!length(vars))return(tibble())
  if(identical(obj$engine,"lightgbm")){
    x<-as.matrix(data.frame(lapply(prevalent[,vars,drop=FALSE],function(z)suppressWarnings(as.numeric(z))),check.names=FALSE));storage.mode(x)<-"double"
    score<-tryCatch(as.numeric(predict(obj$fit,x)),error=function(e)rep(NA_real_,nrow(prevalent)))
  }else{
    cen<-obj$preprocess$center;sca<-obj$preprocess$scale
    if(is.null(cen)||is.null(sca)||!all(vars%in%names(cen))||!all(vars%in%names(sca)))return(tibble())
    cen<-cen[vars];sca<-sca[vars];x<-as.data.frame(prevalent[,vars,drop=FALSE]);x[]<-lapply(x,function(z)suppressWarnings(as.numeric(z)))
    for(j in seq_along(vars))x[!is.finite(x[,j]),j]<-cen[[j]]
    xm<-scale(as.matrix(x),center=cen,scale=sca)
    score<-tryCatch(as.numeric(predict(obj$fit,newx=xm,s=obj$lambda_rule,type="link")),error=function(e)rep(NA_real_,nrow(prevalent)))
  }
  tr<-training_rows|>filter(method==.env$method,split=="training",is.finite(score));m<-mean(tr$score,na.rm=TRUE);s<-sd(tr$score,na.rm=TRUE);if(!is.finite(s)||s==0)s<-1
  tibble(eid=prevalent$eid,row_id=NA_integer_,split="prevalent",time=prevalent[[bvar]],event=1,score=score,
    method=method,n_features=length(vars),n_selected=length(obj$selected%||%character()),engine=obj$engine,
    lambda_rule=obj$lambda_rule%||%NA_character_,score_z=(score-m)/s,
    score_native=score,score_scale=ifelse(obj$engine=="lightgbm","binary model output","Cox linear predictor"))
}

score_coverage_audit <- function(score_rows,expected_methods,prevalent_n){
  expected<-tidyr::crossing(method=expected_methods,split=c("training","validation","prevalent"))|>
    mutate(expected_n=ifelse(split=="prevalent",prevalent_n,NA_integer_))
  got<-score_rows|>group_by(method,split)|>summarise(rows=n(),finite_score=sum(is.finite(score)),finite_z=sum(is.finite(score_z)),.groups="drop")
  expected|>left_join(got,by=c("method","split"))|>mutate(across(c(rows,finite_score,finite_z),~replace_na(.x,0L)),
    coverage=ifelse(rows>0,finite_z/rows,0),status=case_when(rows==0~"missing rows",finite_z==0~"all predictions non-finite",coverage<.95~"partial coverage",TRUE~"ok"))
}

fit_combined_score <- function(dat,score_obj,clinical,tvar,evar,split,label){
  if(is.null(score_obj))return(NULL)
  vars<-unique(c(score_obj$available,clinical))
  if(identical(score_obj$engine,"lightgbm"))fit_lightgbm_score(dat,vars,tvar,evar,split,paste0(label," + clinical"))
  else {
    pf<-score_obj$penalty_factor%||%setNames(rep(1,length(score_obj$available)),score_obj$available)
    pf<-c(pf,setNames(rep(1,length(setdiff(clinical,names(pf)))),setdiff(clinical,names(pf))))
    fit_glmnet_score(dat,vars,tvar,evar,split,paste0(label," + clinical"),
      lambda_rule="lambda.1se",penalty_factor=pf)
  }
}

make_validation_bundle <- function(dat,score_obj,clinical_obj,clinical,tvar,evar,split,label){
  if(is.null(score_obj))return(tibble())
  bio<-score_obj$rows|>filter(split=="validation")|>transmute(eid,row_id,time,event,score,biom_set=label,model="Biomarkers",n_selected=first(score_obj$rows$n_selected))
  comb<-fit_combined_score(dat,score_obj,clinical,tvar,evar,split,label)
  cc<-if(is.null(comb))tibble()else comb$rows|>filter(split=="validation")|>transmute(eid,row_id,time,event,score,biom_set=label,model="Combined",n_selected=first(score_obj$rows$n_selected))
  cl<-if(is.null(clinical_obj))tibble()else clinical_obj$rows|>filter(split=="validation")|>transmute(eid,row_id,time,event,score,biom_set=label,model="Model 0",n_selected=first(score_obj$rows$n_selected))
  bind_rows(bio,cc,cl)|>group_by(biom_set,model)|>mutate(cindex=fcidx(Surv(time,event),score))|>ungroup()
}


# 🚩 Common temporal helpers
yy_score_summary <- function(rows,method,bins=C5_YY_BINS,min_bin_n=C5_MIN_BIN_N){
  d<-rows|>filter(method==.env$method,event==1,is.finite(time),is.finite(score_z));if(!nrow(d))return(list(lines=tibble(),hist=tibble()))
  lim<-max(d$time,na.rm=TRUE);br<-seq(0,lim,length.out=bins+1);mids<-(head(br,-1)+tail(br,-1))/2
  lines<-d|>mutate(bin=cut(time,br,include.lowest=TRUE,labels=FALSE))|>filter(!is.na(bin))|>group_by(split,bin)|>
    summarise(mean=mean(score_z),sd=sd(score_z),se=sd/sqrt(n()),N=n(),N_total=n_distinct(eid),.groups="drop")|>filter(N>=min_bin_n)|>mutate(year=mids[bin])
  hh<-hist(d$time,breaks=br,plot=FALSE);histdf<-tibble(xmin=head(br,-1),xmax=tail(br,-1),mid=mids,n=hh$counts)
  list(lines=lines,hist=histdf)
}

risk_set_pairs <- function(rows,dat,method,score_col=c("score_z","score_native"),seed=SEED,max_cases=C5_PAIR_MAX_CASES){
  score_col<-match.arg(score_col)
  keepcov<-intersect(c("eid","age","sex","ethnic.c"),names(dat))
  d<-rows|>filter(method==.env$method,is.finite(.data[[score_col]]),is.finite(time))|>left_join(dat|>select(all_of(keepcov)),by="eid")
  out<-list();ii<-0L;set.seed(seed)
  for(sp in intersect(c("training","validation"),unique(d$split))){
    ca<-d|>filter(split==sp,event==1)|>arrange(time);co<-d|>filter(split==sp,event==0)
    if(!nrow(ca)||!nrow(co))next
    if(nrow(ca)>max_cases)ca<-ca|>slice_sample(n=max_cases)|>arrange(time)
    for(i in seq_len(nrow(ca))){
      pool<-which(is.finite(co$time)&co$time>=ca$time[i]);if(!length(pool))next
      if("sex"%in%names(d)&&!is.na(ca$sex[i])){same<-pool[co$sex[pool]==ca$sex[i]];if(length(same))pool<-same}
      if("ethnic.c"%in%names(d)&&!is.na(ca$ethnic.c[i])){same<-pool[as.character(co$ethnic.c[pool])==as.character(ca$ethnic.c[i])];if(length(same))pool<-same}
      if("age"%in%names(d)&&is.finite(ca$age[i])){ad<-abs(co$age[pool]-ca$age[i]);pool<-pool[order(ad,na.last=TRUE)][seq_len(min(20,length(pool)))]}
      j<-sample(pool,1);ii<-ii+1L
      out[[ii]]<-tibble(pair_id=paste(sp,i,sep="_"),split=sp,index_time=ca$time[i],case_eid=ca$eid[i],control_eid=co$eid[j],
        case_score=ca[[score_col]][i],control_score=co[[score_col]][j])
    }
  }
  bind_rows(out)
}

pair_long <- function(pairs){
  if(!nrow(pairs))return(tibble())
  bind_rows(
    pairs|>transmute(pair_id,split,index_time,year=-index_time,eid=case_eid,group=ifelse(split=="training","Cases: training","Cases: validation"),score=case_score),
    pairs|>transmute(pair_id,split,index_time,year=-index_time,eid=control_eid,group="Controls",score=control_score))
}

trapz_mean <- function(x,y){
  ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok];if(length(x)<2)return(NA_real_);o<-order(x);x<-x[o];y<-y[o];span<-max(x)-min(x);if(span<=0)return(NA_real_);sum(diff(x)*(head(y,-1)+tail(y,-1))/2)/span
}
temporal_implied_auc <- function(pairs,step=C5_TEMPORAL_STEP){
  if(!nrow(pairs))return(c(dbar=NA_real_,implied_auc=NA_real_))
  z<-pairs|>mutate(bin=floor(index_time/step)*step+step/2,diff=case_score-control_score)|>group_by(bin)|>summarise(d=mean(diff,na.rm=TRUE),N=n(),.groups="drop")|>filter(N>=5)
  dbar<-trapz_mean(z$bin,z$d);c(dbar=dbar,implied_auc=ifelse(is.finite(dbar),pnorm(dbar/sqrt(2)),NA_real_))
}
match_bootstrap_pairs <- function(v,seed,max_cases=C5_PAIR_MAX_CASES){
  set.seed(seed);ca<-v|>filter(event==1);co<-v|>filter(event==0);if(!nrow(ca)||!nrow(co))return(tibble())
  if(nrow(ca)>max_cases)ca<-ca|>slice_sample(n=max_cases)
  out<-list();for(i in seq_len(nrow(ca))){pool<-which(is.finite(co$time)&co$time>=ca$time[i]);if(!length(pool))next
    if("sex"%in%names(v)&&!is.na(ca$sex[i])){same<-pool[co$sex[pool]==ca$sex[i]];if(length(same))pool<-same}
    if("ethnic.c"%in%names(v)&&!is.na(ca$ethnic.c[i])){same<-pool[as.character(co$ethnic.c[pool])==as.character(ca$ethnic.c[i])];if(length(same))pool<-same}
    if("age"%in%names(v)&&is.finite(ca$age[i])){ad<-abs(co$age[pool]-ca$age[i]);pool<-pool[order(ad,na.last=TRUE)][seq_len(min(20,length(pool)))]}
    j<-sample(pool,1);out[[length(out)+1]]<-tibble(index_time=ca$time[i],case_score=ca$score_z[i],control_score=co$score_z[j])}
  bind_rows(out)
}
bootstrap_auc_link <- function(method_rows,dat,B=C5_BOOT){
  keepcov<-intersect(c("eid","age","sex","ethnic.c"),names(dat))
  v<-method_rows|>filter(split=="validation")|>left_join(dat|>select(all_of(keepcov)),by="eid")|>filter(is.finite(score_z),!is.na(event),is.finite(time))
  if(nrow(v)<200)return(tibble())
  ca<-v|>filter(event==1);co<-v|>filter(event==0);if(!nrow(ca)||!nrow(co))return(tibble())
  base_pairs<-match_bootstrap_pairs(v,SEED+1999,C5_PAIR_MAX_CASES)
  map_dfr(seq_len(B),function(b){
    set.seed(SEED+1000+b)
    vb<-bind_rows(ca[sample(seq_len(nrow(ca)),min(nrow(ca),C5_BOOT_MAX_CASES),replace=TRUE),,drop=FALSE],
      co[sample(seq_len(nrow(co)),min(nrow(co),C5_BOOT_MAX_CONTROLS),replace=TRUE),,drop=FALSE])
    roc_auc<-weighted_time_auc(vb$time,vb$event,vb$score_z,HORIZON)
    pairs<-if(nrow(base_pairs))base_pairs[sample(seq_len(nrow(base_pairs)),nrow(base_pairs),replace=TRUE),,drop=FALSE]else tibble();ta<-temporal_implied_auc(pairs)
    tibble(bootstrap=b,roc_auc=roc_auc,temporal_d=ta[["dbar"]],trajectory_implied_auc=ta[["implied_auc"]])
  })
}

# Held-out mechanism benchmark.  Feature ranking, signs, Cox weights, means and
# scales are learned only in the outer training split.  The comparators are a
# locked training-selected top-1 biomarker, an unweighted top-N mean and a
# training marginal-Cox-weighted top-N score.  We also report the distribution of the N
# individual biomarkers; no "best" biomarker is selected on validation data.
topn_score_benchmark <- function(dat,ranked,screen,tvar,evar,split,clinical=character(),
                                 max_n=C5_TOPN_MAX,B=C5_TOPN_BOOT){
  empty<-list(rows=tibble(),summary=tibble(),bootstrap=tibble(),individuals=tibble(),individual_range=tibble())
  ranked<-head(intersect(ranked,names(dat)),max(1L,max_n));if(!length(ranked))return(empty)
  it<-which(split=="training");iv<-which(split=="validation");if(length(it)<300||length(iv)<100)return(empty)
  center<-vapply(dat[it,ranked,drop=FALSE],function(x)mean(suppressWarnings(as.numeric(x)),na.rm=TRUE),numeric(1))
  scale0<-vapply(dat[it,ranked,drop=FALSE],function(x)sd(suppressWarnings(as.numeric(x)),na.rm=TRUE),numeric(1))
  keep<-is.finite(center)&is.finite(scale0)&scale0>0;ranked<-ranked[keep];center<-center[keep];scale0<-scale0[keep]
  if(!length(ranked))return(empty)
  make_x<-function(ii){
    x<-as.matrix(data.frame(lapply(dat[ii,ranked,drop=FALSE],function(z)suppressWarnings(as.numeric(z))),check.names=FALSE));storage.mode(x)<-"double"
    for(j in seq_along(ranked))x[!is.finite(x[,j]),j]<-center[j]
    sweep(sweep(x,2,center,"-"),2,scale0,"/")
  }
  xtr<-make_x(it);xva<-make_x(iv);w<-screen$beta[match(ranked,screen$term)];w[!is.finite(w)|w==0]<-1
  clinical<-intersect(clinical,names(dat));base<-tibble(eid=dat$eid[iv],time=dat[[tvar]][iv],event=dat[[evar]][iv],row_id=iv)
  standardize_pair<-function(tr,va){m<-mean(tr,na.rm=TRUE);s<-sd(tr,na.rm=TRUE);if(!is.finite(s)||s==0)s<-1;list(tr=(tr-m)/s,va=(va-m)/s)}
  metric_vec<-function(score){
    d<-base|>mutate(score=score);auc<-weighted_time_auc(d$time,d$event,d$score,HORIZON)
    hz<-d|>filter(is.finite(score),is.finite(time),!is.na(event))|>
      mutate(class=case_when(event==1&time<=HORIZON~"Case",time>HORIZON~"Control",TRUE~NA_character_))|>filter(!is.na(class))
    ca<-hz$score[hz$class=="Case"];co<-hz$score[hz$class=="Control"]
    delta<-if(length(ca)&&length(co))mean(ca)-mean(co)else NA_real_
    pooled<-if(length(ca)>1&&length(co)>1)sqrt(((length(ca)-1)*var(ca)+(length(co)-1)*var(co))/(length(ca)+length(co)-2))else NA_real_
    dd<-dat[iv,unique(c(tvar,evar,clinical)),drop=FALSE];dd$.score<-score;dd<-dd[complete.cases(dd),,drop=FALSE];dd<-dd[dd[[tvar]]>0,,drop=FALSE]
    rhs<-c(".score",clinical);fit<-if(nrow(dd)>=100&&sum(dd[[evar]]==1)>=20)tryCatch(coxph(as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(bt(rhs),collapse=" + "))),dd,ties="efron"),error=function(e)NULL)else NULL
    sm<-if(is.null(fit))NULL else coef(summary(fit));beta<-if(!is.null(sm)&&".score"%in%rownames(sm))sm[".score","coef"]else NA_real_;se<-if(!is.null(sm)&&".score"%in%rownames(sm))sm[".score","se(coef)"]else NA_real_
    inc<-d|>filter(event==1,is.finite(time),is.finite(score));slope<-if(nrow(inc)>=20)tryCatch(unname(coef(lm(score~time,data=inc))[2]),error=function(e)NA_real_)else NA_real_
    tibble(AUC=auc,beta=beta,beta_se=se,beta_lo=beta-1.96*se,beta_hi=beta+1.96*se,
      mean_separation=delta,within_group_SD=pooled,cohen_d=delta/pooled,score_SD=sd(score,na.rm=TRUE),
      prediagnostic_steepness=-slope,cases=sum(d$event==1,na.rm=TRUE),N_validation=nrow(d))
  }
  boot_vec<-function(score,n,model_index){
    if(B<=0)return(tibble());ca<-which(base$event==1&is.finite(score));co<-which(base$event==0&is.finite(score))
    if(length(ca)<20||length(co)<30)return(tibble())
    map_dfr(seq_len(B),function(b){
      set.seed(SEED+50000+1000*n+100*model_index+b)
      ii<-c(sample(ca,min(length(ca),C5_BOOT_MAX_CASES),replace=TRUE),sample(co,min(length(co),C5_BOOT_MAX_CONTROLS),replace=TRUE))
      d<-base[ii,,drop=FALSE]|>mutate(score=score[ii]);auc<-weighted_time_auc(d$time,d$event,d$score,HORIZON)
      hz<-d|>mutate(class=case_when(event==1&time<=HORIZON~"Case",time>HORIZON~"Control",TRUE~NA_character_))|>filter(!is.na(class))
      ca0<-hz$score[hz$class=="Case"];co0<-hz$score[hz$class=="Control"]
      delta<-if(length(ca0)&&length(co0))mean(ca0)-mean(co0)else NA_real_;pooled<-if(length(ca0)>1&&length(co0)>1)sqrt(((length(ca0)-1)*var(ca0)+(length(co0)-1)*var(co0))/(length(ca0)+length(co0)-2))else NA_real_
      inc<-d|>filter(event==1,is.finite(time),is.finite(score));slope<-if(nrow(inc)>=20)tryCatch(unname(coef(lm(score~time,data=inc))[2]),error=function(e)NA_real_)else NA_real_
      tibble(bootstrap=b,AUC=auc,mean_separation=delta,within_group_SD=pooled,cohen_d=delta/pooled,prediagnostic_steepness=-slope)
    })
  }
  extra_n<-if(length(ranked)>=15L)seq(15L,length(ranked),by=5L)else integer();ns<-unique(pmin(length(ranked),c(seq_len(min(10L,length(ranked))),extra_n,length(ranked))))
  individual<-map_dfr(seq_along(ranked),function(j)metric_vec(xva[,j]*sign(w[j]))|>mutate(rank=j,feature=ranked[j]))
  summaries<-list();boots<-list();trajectory<-list();kk<-0L;target_n<-ns[which.min(abs(ns-min(C5_MECHANISM_N,length(ranked))))]
  for(n in ns){
    ww<-w[seq_len(n)];sgn<-sign(ww);sgn[sgn==0]<-1
    weighted<-standardize_pair(as.numeric(xtr[,seq_len(n),drop=FALSE]%*%ww),as.numeric(xva[,seq_len(n),drop=FALSE]%*%ww))
    unweighted<-standardize_pair(rowMeans(sweep(xtr[,seq_len(n),drop=FALSE],2,sgn,"*")),rowMeans(sweep(xva[,seq_len(n),drop=FALSE],2,sgn,"*")))
    top1<-list(tr=xtr[,1]*sign(w[1]),va=xva[,1]*sign(w[1]))
    models<-list(`Top-1 biomarker (locked)`=top1,`Unweighted top-N mean`=unweighted,
      `Marginal-Cox-weighted top-N`=weighted)
    cm<-suppressWarnings(cor(xtr[,seq_len(n),drop=FALSE],use="pairwise.complete.obs"));avgcor<-if(n>1)mean(abs(cm[row(cm)!=col(cm)]),na.rm=TRUE)else 0
    for(mi in seq_along(models)){
      nm<-names(models)[mi];sc<-models[[mi]]$va;weights0<-if(mi==1)c(1,rep(0,n-1))else if(mi==2)sgn else ww
      effn<-if(sum(abs(weights0))>0)1/sum((abs(weights0)/sum(abs(weights0)))^2)else NA_real_
      kk<-kk+1L;summaries[[kk]]<-metric_vec(sc)|>mutate(N=n,model=nm,feature=if(nm=="Top-1 biomarker (locked)")ranked[1]else paste0("top ",n),average_abs_correlation=avgcor,effective_N=effn)
      boots[[kk]]<-boot_vec(sc,n,mi)|>mutate(N=n,model=nm)
      if(n==target_n)trajectory[[kk]]<-base|>transmute(eid,time,event,N=n,model=nm,score=sc)
    }
  }
  summary<-bind_rows(summaries);boot<-bind_rows(boots)
  if(nrow(boot))summary<-summary|>left_join(boot|>group_by(N,model)|>summarise(AUC_boot_SD=sd(AUC,na.rm=TRUE),separation_boot_SD=sd(mean_separation,na.rm=TRUE),
    slope_boot_SD=sd(prediagnostic_steepness,na.rm=TRUE),AUC_lo=quantile(AUC,.025,na.rm=TRUE),AUC_hi=quantile(AUC,.975,na.rm=TRUE),.groups="drop"),by=c("N","model"))
  for(nm in c("AUC_boot_SD","separation_boot_SD","slope_boot_SD","AUC_lo","AUC_hi"))if(!nm%in%names(summary))summary[[nm]]<-NA_real_
  irange<-map_dfr(ns,function(n)individual|>filter(rank<=n)|>summarise(N=n,AUC_q25=quantile(AUC,.25,na.rm=TRUE),AUC_median=median(AUC,na.rm=TRUE),AUC_q75=quantile(AUC,.75,na.rm=TRUE),AUC_max=max(AUC,na.rm=TRUE),
    beta_q25=quantile(beta,.25,na.rm=TRUE),beta_median=median(beta,na.rm=TRUE),beta_q75=quantile(beta,.75,na.rm=TRUE),cohen_d_median=median(cohen_d,na.rm=TRUE)))
  list(rows=bind_rows(trajectory),summary=summary,bootstrap=boot,individuals=individual,individual_range=irange)
}


# 🚩 The five reusable row functions.  All C5 score figures call these directly.
c5_row1_score_profile <- function(score_rows,label,n_selected=NA_integer_){
  d<-score_rows|>filter(method==.env$label,split=="validation",is.finite(score_z),is.finite(time),time>0,!is.na(event))
  ttl<-paste0(label," (N=",ifelse(is.finite(n_selected),format(as.integer(n_selected),big.mark=","),"NA"),")")
  if(nrow(d)<100)return(blank_plot(ttl,"Insufficient validation scores"))
  status<-factor(d$event,levels=0:1,labels=c("Controls","Cases"));xr<-range(d$score_z,na.rm=TRUE);if(diff(xr)<=0)xr<-xr+c(-1,1)
  den<-map_dfr(levels(status),function(g){x<-d$score_z[status==g];if(length(x)<5)return(tibble());z<-density(x,from=xr[1],to=xr[2],n=256,na.rm=TRUE);tibble(x=z$x,density=z$y,status=g)})
  rate<-d|>mutate(bin=ntile(score_z,10))|>group_by(bin)|>summarise(x=mean(score_z),rate=1000*sum(event)/sum(time),N=n(),.groups="drop")|>filter(is.finite(rate))
  dmax<-max(den$density,na.rm=TRUE);if(!is.finite(dmax)||dmax<=0)dmax<-1
  rr<-range(rate$rate,na.rm=TRUE);if(!all(is.finite(rr)))rr<-c(0,1);if(diff(rr)<1e-9)rr<-rr+c(-.5,.5)
  a<-.78*dmax/diff(rr);offset<-.08*dmax-a*rr[1];rate<-rate|>mutate(y=a*rate+offset)
  ggplot()+geom_area(data=den,aes(x,density,fill=status,group=status),alpha=.38,position="identity")+
    geom_line(data=den,aes(x,density,color=status,group=status),linewidth=.55)+
    geom_line(data=rate,aes(x,y),color="black",linewidth=.85)+geom_point(data=rate,aes(x,y),shape=21,fill="white",size=2.2,stroke=.75)+
    scale_y_continuous(name="Density",sec.axis=sec_axis(~(.x-offset)/a,name="Incidence / 1,000 person-years"))+
    scale_fill_manual(values=c(Controls="grey78",Cases="#D9A066"))+scale_color_manual(values=c(Controls="grey48",Cases="#A8612B"))+
    coord_cartesian(xlim=xr,clip="off")+labs(title=ttl,x="Biomarker score (s.d.)",fill=NULL,color=NULL)+theme_5c(9)+
    theme(legend.position="top",legend.justification="center")
}

ipcw_roc_curve <- function(time,event,score,horizon=HORIZON){
  ok<-is.finite(time)&!is.na(event)&is.finite(score)&time>0
  time<-time[ok];event<-event[ok];score<-score[ok];if(length(score)<100)return(tibble())
  ip<-ipcw_at_horizon(time,event,horizon);cw<-ifelse(ip$case,ip$weight,0);tw<-ifelse(ip$control,ip$weight,0)
  keep<-cw>0|tw>0;if(sum(cw[keep])<=0||sum(tw[keep])<=0)return(tibble())
  d<-tibble(score=score[keep],cw=cw[keep],tw=tw[keep])|>group_by(score)|>summarise(cw=sum(cw),tw=sum(tw),.groups="drop")|>arrange(desc(score))|>
    mutate(tpr=cumsum(cw)/sum(cw),fpr=cumsum(tw)/sum(tw))
  bind_rows(tibble(score=Inf,tpr=0,fpr=0),d|>select(score,tpr,fpr),tibble(score=-Inf,tpr=1,fpr=1))
}

c5_row2_roc <- function(pred,label,boot=tibble()){
  d<-pred|>filter(biom_set==.env$label);if(!nrow(d))return(blank_plot("Validation ROC","No validation predictions"))
  roc<-d|>group_by(model)|>group_modify(~{
    z<-ipcw_roc_curve(.x$time,.x$event,.x$score,HORIZON);if(!nrow(z))return(tibble())
    z|>mutate(auc=weighted_time_auc(.x$time,.x$event,.x$score,HORIZON))
  })|>ungroup()
  if(!nrow(roc))return(blank_plot("Validation ROC","ROC unavailable"))
  al<-roc|>group_by(model)|>summarise(auc=first(auc),.groups="drop")|>mutate(txt=sprintf("%s: %.3f",case_when(model=="Model 0"~"Clinical",model=="Biomarkers"~"Biomarkers",TRUE~"Combined"),auc),x=.03,y=c(.97,.90,.83)[match(model,c("Model 0","Biomarkers","Combined"))])
  p<-ggplot(roc,aes(fpr,tpr,color=model))+geom_abline(slope=1,intercept=0,linetype=3,color="grey78")+geom_step(linewidth=.9)+
    geom_text(data=al,aes(x=x,y=y,label=txt,color=model),hjust=0,vjust=1,inherit.aes=FALSE,show.legend=FALSE,size=2.8,fontface="bold")+
    scale_color_manual(values=c(`Model 0`="grey35",Biomarkers="#C77732",Combined="#287A58"),breaks=c("Biomarkers","Combined","Model 0"),labels=c("Biomarkers","Combined","Clinical"))+
    labs(x="False positive rate",y="True positive rate",color=NULL,subtitle=paste0(HORIZON,"-year IPCW cumulative/dynamic ROC"))+theme_5c(9)+theme(legend.position="top")
  if(nrow(boot)>=20){
    rr<-suppressWarnings(cor(boot$trajectory_implied_auc,boot$roc_auc,use="complete.obs"))
    ins<-ggplot(boot,aes(trajectory_implied_auc,roc_auc))+geom_abline(slope=1,intercept=0,linetype=2,color="grey65")+
      geom_point(alpha=.28,size=.7)+geom_smooth(method="lm",se=FALSE,linewidth=.45,color="grey25")+
      annotate("text",x=Inf,y=-Inf,label=sprintf("r=%.2f",rr),hjust=1.08,vjust=-.4,size=2.4,fontface="bold")+
      labs(x="Temporal AUC",y="ROC AUC")+theme_classic(base_size=6)+theme(axis.title=element_text(face="bold"),axis.text=element_text(size=5),plot.background=element_rect(fill="white",color="grey75"))
    p<-p+patchwork::inset_element(ins,left=.38,bottom=.10,right=.84,top=.55,align_to="panel",on_top=TRUE)
  }
  p+labs(title=NULL)
}

c5_row3_followup <- function(score_rows,label,horizon=HORIZON){
  allm<-score_rows|>filter(method==.env$label,is.finite(score_z));tr<-allm|>filter(split=="training")
  if(nrow(tr)<100)return(blank_plot("Diagnosis burden across baseline","No training score distribution"))
  br<-unique(as.numeric(quantile(tr$score_z,probs=seq(0,1,.2),na.rm=TRUE)));if(length(br)<6)br<-seq(min(tr$score_z),max(tr$score_z),length.out=6)
  addq<-function(x)x|>mutate(quintile=factor(cut(score_z,breaks=br,include.lowest=TRUE,labels=1:5),levels=1:5,labels=paste0("Q",1:5)))
  pos<-addq(allm|>filter(split=="validation",is.finite(time),time>0,!is.na(event)))
  cp<-if(!nrow(pos))tibble()else pos|>filter(!is.na(quintile))|>group_by(quintile)|>group_modify(~{sf<-survfit(Surv(time,event)~1,data=.x);tibble(time=sf$time,cuminc=1-sf$surv)})|>ungroup()|>mutate(phase="Incident follow-up")
  p_inc<-if(!nrow(cp))blank_plot("Incident follow-up","No validation survival curve")else ggplot(cp,aes(time,cuminc,color=quintile))+
    geom_line(linewidth=.78)+coord_cartesian(xlim=c(0,C5_TIME_CAP))+scale_color_brewer(palette="YlOrBr",direction=1)+
    labs(title="Incident follow-up",x="Years after baseline",y="Cumulative incidence",color="Score quintile")+theme_5c(8)+theme(legend.position="top")
  # Do not draw a pseudo-survival curve before baseline. Prevalent cases are a
  # cross-sectional disease-state comparison, so report their proportion in
  # each score quintile using the full scored cohort.
  base<-addq(allm|>filter(split%in%c("training","validation","prevalent")))|>filter(!is.na(quintile))
  prev<-base|>group_by(quintile)|>summarise(N=n(),prevalent=sum(split=="prevalent"),proportion=prevalent/N,.groups="drop")
  p_prev<-if(!nrow(prev))blank_plot("Baseline prevalence","No prevalent score")else ggplot(prev,aes(quintile,proportion,fill=quintile))+
    geom_col(width=.72,show.legend=FALSE)+geom_text(aes(label=scales::percent(proportion,accuracy=.1)),vjust=-.25,size=2.6,fontface="bold")+
    scale_fill_brewer(palette="YlOrBr",direction=1)+scale_y_continuous(labels=scales::percent,expand=expansion(mult=c(0,.16)))+
    labs(title="Baseline prevalence",subtitle="Cross-sectional; not a reverse-time survival curve",x="Score quintile",y="Prevalent proportion")+theme_5c(8)
  p_prev|p_inc
}

c5_row4_yy <- function(score_rows,label){
  d<-score_rows|>filter(method==.env$label,event==1,is.finite(time),time>=-C5_TIME_CAP,time<=C5_TIME_CAP,is.finite(score_z))|>
    mutate(group=case_when(split=="prevalent"~"Prevalent (Yang)",split=="training"~"Incident: training",TRUE~"Incident: validation"),
      bin=floor(time/C5_TEMPORAL_STEP)*C5_TEMPORAL_STEP+C5_TEMPORAL_STEP/2)
  ln<-d|>group_by(group,bin)|>summarise(mean=mean(score_z),N=n(),N_total=n_distinct(eid),.groups="drop")|>mutate(reliable=N>=C5_MIN_BIN_N)
  h<-d|>count(bin,name="n")|>transmute(xmin=bin-C5_TEMPORAL_STEP/2,xmax=bin+C5_TEMPORAL_STEP/2,n)
  if(!nrow(ln))return(blank_plot("Temporal score trajectory","No sufficiently populated bins"))
  ln<-ln|>group_by(group)|>mutate(group_label=sprintf("%s (N=%s)",first(group),format(first(N_total),big.mark=",")))|>ungroup()
  yr<-range(ln$mean,na.rm=TRUE);if(diff(yr)<.2)yr<-yr+c(-.5,.5);pad<-.13*diff(yr);yr<-yr+c(-pad,pad)
  count_max<-max(h$n,1,na.rm=TRUE);plot_max<-count_max*1.15;scale_factor<-plot_max/diff(yr)
  ln<-ln|>mutate(y_plot=(mean-yr[1])*scale_factor)
  lev<-unique(ln|>arrange(factor(group,levels=c("Prevalent (Yang)","Incident: training","Incident: validation")))|>pull(group_label));base<-sub(" \\(N=.*$","",lev);pal0<-c(`Prevalent (Yang)`="grey45",`Incident: training`="#2C7FB8",`Incident: validation`="#D95F02");pal<-setNames(unname(pal0[base]),lev)
  missing_side<-tibble(x=c(-C5_TIME_CAP/2,C5_TIME_CAP/2),label=c(if(any(d$split=="prevalent"))NA_character_ else "Prevalence trajectory unavailable",if(any(d$split!="prevalent"))NA_character_ else "Incidence trajectory unavailable"))|>filter(!is.na(label))
  ggplot()+geom_rect(data=h,aes(xmin=xmin,xmax=xmax,ymin=0,ymax=n,fill=bin>=0),alpha=.31,inherit.aes=FALSE)+
    scale_fill_manual(values=c(`TRUE`="#78B7E5",`FALSE`="grey60"),guide="none")+
    geom_vline(xintercept=0,linetype=2,color="grey38")+
    geom_hline(yintercept=(0-yr[1])*scale_factor,linetype=3,color="grey62")+geom_smooth(data=ln|>filter(reliable),aes(bin,y_plot,color=group_label,group=group_label),method="loess",se=FALSE,span=.65,linewidth=1.0)+geom_point(data=ln,aes(bin,y_plot,color=group_label,alpha=reliable),size=1.45)+
    geom_text(data=missing_side,aes(x=x,y=plot_max*.52,label=label),inherit.aes=FALSE,color="grey45",fontface="bold",size=3)+
    scale_alpha_manual(values=c(`TRUE`=.78,`FALSE`=.25),guide="none")+scale_y_continuous(name="Patient count",limits=c(0,plot_max),expand=expansion(mult=c(0,.02)),sec.axis=sec_axis(~.x/scale_factor+yr[1],name="Biomarker score (s.d.)"))+scale_color_manual(values=pal)+
    coord_cartesian(xlim=c(-C5_TIME_CAP,C5_TIME_CAP))+labs(x="Recorded diagnosis relative to baseline blood draw (years)",color=NULL)+theme_5c(9)+
    guides(color=guide_legend(nrow=1,byrow=TRUE))+theme(legend.position="top")
}

c5_row5_preclinical <- function(score_rows,label){
  d<-score_rows|>filter(method==.env$label,is.finite(score_z))
  cases<-d|>filter(event==1,is.finite(time),time>=-C5_TIME_CAP,time<=C5_TIME_CAP)|>
    mutate(group=case_when(split=="prevalent"~"Cases: prevalent",split=="training"~"Cases: training",TRUE~"Cases: validation"),
      bin=floor(time/C5_TEMPORAL_STEP)*C5_TEMPORAL_STEP+C5_TEMPORAL_STEP/2)
  sm<-cases|>group_by(group,bin)|>summarise(mean=mean(score_z),se=sd(score_z)/sqrt(n()),N=n(),.groups="drop")|>filter(N>=C5_MIN_BIN_N)
  if(!nrow(sm))return(blank_plot("Baseline-anchored native score","No sufficiently populated case bins"))
  ctrl<-d|>filter(event==0,split=="validation")|>summarise(mean=mean(score_z,na.rm=TRUE),se=sd(score_z,na.rm=TRUE)/sqrt(sum(is.finite(score_z))))
  nlab<-cases|>group_by(group)|>summarise(N=n_distinct(eid),.groups="drop")|>mutate(group_label=paste0(group," (N=",format(N,big.mark=","),")"))
  sm<-sm|>left_join(nlab,by="group")
  # Standardized scores are comparable across Cox and binary engines. They are
  # deliberately not labelled as calibrated probabilities.
  ylo<-min(sm$mean-sm$se,na.rm=TRUE);yhi<-max(sm$mean+sm$se,na.rm=TRUE)
  if(!all(is.finite(c(ylo,yhi))))c(ylo,yhi)<-c(0,1)
  if(yhi-ylo<.02){mid<-(ylo+yhi)/2;ylo<-mid-.01;yhi<-mid+.01}
  ypad<-.12*(yhi-ylo);ylim_native<-c(ylo-ypad,yhi+ypad)
  gl<-unique(nlab$group_label);base_group<-sub(" \\(N=.*$","",gl);pal0<-c(`Cases: prevalent`="grey45",`Cases: training`="#2C7FB8",`Cases: validation`="#D95F02");pal<-setNames(unname(pal0[base_group]),gl)
  missing_side<-tibble(x=c(-C5_TIME_CAP/2,C5_TIME_CAP/2),label=c(if(any(cases$split=="prevalent"))NA_character_ else "Prevalence score unavailable",if(any(cases$split!="prevalent"))NA_character_ else "Incidence score unavailable"))|>filter(!is.na(label))
  ggplot(sm,aes(bin,mean,color=group_label,group=group_label))+
    geom_hline(data=ctrl,aes(yintercept=mean),inherit.aes=FALSE,linetype=3,color="grey45")+
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=.14,linewidth=.35,alpha=.55)+geom_point(size=1.8)+
    geom_smooth(method="loess",se=FALSE,span=.75,linewidth=1.0)+geom_vline(xintercept=0,linetype=2,color="grey50")+
    geom_text(data=missing_side,aes(x=x,y=mean(ylim_native),label=label),inherit.aes=FALSE,color="grey45",fontface="bold",size=3)+
    scale_color_manual(values=pal)+coord_cartesian(xlim=c(-C5_TIME_CAP,C5_TIME_CAP),ylim=ylim_native)+
    labs(x="Recorded diagnosis relative to baseline blood draw (years)",y="Model score (training s.d.)",color=NULL,
      subtitle="Dotted line: validation-control mean")+theme_5c(9)+
    guides(color=guide_legend(nrow=1,byrow=TRUE))+theme(legend.position="top",legend.text=element_text(size=7.4,face="bold"))
}

c5_make_5row_grid <- function(labels,pred,score_rows,pairs_native,boot_by_method,title=NULL){
  if(!length(labels))return(blank_plot(title%||%"C5 prediction","No score was available"))
  # Preserve the requested grid even when a model is unavailable; otherwise an
  # empty C4 cache silently collapses Fig3 into one completely blank canvas.
  nsel<-score_rows|>group_by(method)|>summarise(n_selected=as.integer(median(n_selected,na.rm=TRUE)),
    n_features=as.integer(median(n_features,na.rm=TRUE)),.groups="drop")
  getn<-function(lb){col<-if(lb=="User specified")"n_features"else"n_selected";z<-nsel[[col]][match(lb,nsel$method)];if(length(z)&&is.finite(z))z else NA_integer_}
  method_rows<-lapply(labels,function(lb)wrap_plots(list(
    c5_row1_score_profile(score_rows,lb,getn(lb)),
    c5_row2_roc(pred,lb,boot_by_method[[lb]]%||%tibble()),
    c5_row3_followup(score_rows,lb),
    c5_row4_yy(score_rows,lb),
    c5_row5_preclinical(score_rows,lb)),nrow=1,widths=c(1.15,1.05,1,1,1.08)))
  spaced<-list();for(i in seq_along(method_rows)){spaced[[length(spaced)+1]]<-method_rows[[i]];if(i<length(method_rows))spaced[[length(spaced)+1]]<-plot_spacer()}
  p<-wrap_plots(spaced,ncol=1,heights=rep(c(1,.08),length.out=length(spaced)))
  if(!is.null(title))p<-p+plot_annotation(title=title,theme=theme(plot.title=element_text(face="bold",size=16,hjust=.5)))
  p
}


# 🚩 Lead-time and consolidation figures (Figs 4-10)
make_individual_lead_rows <- function(dat,features,screen,split,tvar,evar){
  features<-intersect(unique(features),names(dat));it<-which(split=="training");iv<-which(split=="validation")
  map_dfr(features,function(x){
    tr<-suppressWarnings(as.numeric(dat[[x]][it]));va<-suppressWarnings(as.numeric(dat[[x]][iv]));m<-mean(tr,na.rm=TRUE);s<-sd(tr,na.rm=TRUE)
    if(!is.finite(s)||s==0)return(tibble());sgn<-sign(screen$beta[match(x,screen$term)]);if(!is.finite(sgn)||sgn==0)sgn<-1
    tibble(eid=dat$eid[iv],split="validation",time=dat[[tvar]][iv],event=dat[[evar]][iv],score_z=(va-m)/s*sgn,
      method=x,kind="Single biomarker",n_features=1L,n_selected=1L)
  })
}

leadtime_auc_table <- function(score_rows,methods,max_h=C5_TIME_CAP){
  if(!requireNamespace("pROC",quietly=TRUE))return(tibble())
  if(!"kind"%in%names(score_rows))score_rows$kind<-"Omic score"
  hs<-sort(unique(c(.5,seq(1,max(1,floor(max_h)),by=1))))
  map_dfr(intersect(methods,unique(score_rows$method)),function(mm){
    d0<-score_rows|>filter(method==mm,split=="validation",is.finite(time),is.finite(score_z),!is.na(event))
    kind0<-first(d0$kind)%||%"Omic score"
    map_dfr(hs,function(h){
      # Landmark-style discrimination: incident cases diagnosed at least h years
      # after baseline versus controls observed event-free for at least h years.
      d<-d0|>filter(time>=h)
      nc<-sum(d$event==1);nn<-sum(d$event==0)
      if(nc<20||nn<50)return(tibble(method=mm,kind=kind0,horizon=h,cases=nc,controls=nn,AUC=NA_real_,lo=NA_real_,hi=NA_real_,
        AUC_definition="minimum-lead-time case/control AUC"))
      roc<-tryCatch(pROC::roc(d$event,d$score_z,quiet=TRUE,direction="<"),error=function(e)NULL)
      if(is.null(roc))return(tibble(method=mm,kind=kind0,horizon=h,cases=nc,controls=nn,AUC=NA_real_,lo=NA_real_,hi=NA_real_,
        AUC_definition="minimum-lead-time case/control AUC"))
      ci<-tryCatch(as.numeric(pROC::ci.auc(roc,method="delong")),error=function(e)c(NA_real_,NA_real_,NA_real_))
      tibble(method=mm,kind=kind0,horizon=h,cases=nc,controls=nn,AUC=as.numeric(roc$auc),lo=ci[1],hi=ci[3],
        AUC_definition="minimum-lead-time case/control AUC")
    })
  })
}

leadtime_headline <- function(lead,min_cases=C5_LEAD_MIN_CASES){
  if(!nrow(lead)||!all(c("method","kind","horizon","AUC","lo","hi","cases")%in%names(lead)))
    return(tibble(method=character(),kind=character(),last_supported_horizon=numeric(),
      AUC_at_horizon=numeric(),lo_at_horizon=numeric(),hi_at_horizon=numeric(),cases_at_horizon=integer()))
  lead|>group_by(method,kind)|>group_modify(function(z,key){
    z<-z|>arrange(horizon);supported<-is.finite(z$AUC)&is.finite(z$lo)&z$lo>.5&z$cases>=min_cases
    first_fail<-match(FALSE,supported,nomatch=length(supported)+1L);last<-first_fail-1L
    if(last<1)return(tibble());tibble(last_supported_horizon=z$horizon[last],AUC_at_horizon=z$AUC[last],lo_at_horizon=z$lo[last],hi_at_horizon=z$hi[last],cases_at_horizon=z$cases[last])
  })|>ungroup()
}

prediction_window_auc_table <- function(score_rows,methods){
  if(!requireNamespace("pROC",quietly=TRUE))return(tibble())
  if(!"kind"%in%names(score_rows))score_rows$kind<-"Omic score"
  windows<-tibble(window=c("Within 5 y","Within 10 y","Over 10 y","Over 12 y"),type=c("within","within","over","over"),cut=c(5,10,10,12))
  ans<-map_dfr(intersect(methods,unique(score_rows$method)),function(mm){
    d0<-score_rows|>filter(method==mm,split=="validation",is.finite(time),is.finite(score_z),!is.na(event));kind0<-first(d0$kind)%||%"Omic score"
    map_dfr(seq_len(nrow(windows)),function(i){
      w<-windows[i,];type0<-w$type[[1]];cut0<-w$cut[[1]];window0<-w$window[[1]]
      d<-if(type0=="within")d0|>mutate(class=case_when(event==1&time<=cut0~1,time>cut0~0,TRUE~NA_real_))|>filter(!is.na(class)) else d0|>filter(time>cut0)|>mutate(class=event)
      nc<-sum(d$class==1);nn<-sum(d$class==0);if(nc<20||nn<50)return(tibble(method=mm,kind=kind0,window=window0,cases=nc,controls=nn,AUC=NA_real_,lo=NA_real_,hi=NA_real_,
        AUC_definition="prespecified binary prediction-window AUC"))
      roc<-tryCatch(pROC::roc(d$class,d$score_z,quiet=TRUE,direction="<"),error=function(e)NULL);if(is.null(roc))return(tibble())
      ci<-tryCatch(as.numeric(pROC::ci.auc(roc,method="delong")),error=function(e)c(NA_real_,NA_real_,NA_real_))
      tibble(method=mm,kind=kind0,window=window0,cases=nc,controls=nn,AUC=as.numeric(roc$auc),lo=ci[1],hi=ci[3],
        AUC_definition="prespecified binary prediction-window AUC")
    })
  })
  if(!nrow(ans))return(tibble(method=character(),kind=character(),window=factor(levels=windows$window),
    cases=integer(),controls=integer(),AUC=numeric(),lo=numeric(),hi=numeric(),AUC_definition=character()))
  ans|>mutate(window=factor(window,levels=windows$window))
}

plot_leadtime_prediction <- function(lead,windows=tibble()){
  if(!nrow(lead)||!any(is.finite(lead$AUC)))return(blank_plot("Prediction before onset","No lead-time estimate met the case/control requirement"))
  ok<-lead|>filter(is.finite(AUC));hmax<-leadtime_headline(lead)
  curve_panel<-function(kind0,title){
    d<-ok|>filter(kind==kind0);hm<-hmax|>filter(kind==kind0)
    if(!nrow(d))return(blank_plot(title,"No eligible method"))
    ggplot(d,aes(horizon,AUC))+geom_hline(yintercept=.5,linetype=2,color="grey55")+
      geom_ribbon(aes(ymin=lo,ymax=hi),fill="#6BAED6",alpha=.20)+geom_line(color="#2C7FB8",linewidth=.85)+geom_point(color="#2C7FB8",size=1.25)+
      geom_vline(data=hm,aes(xintercept=last_supported_horizon),color="#D7301F",linetype=2,linewidth=.8)+
      geom_point(data=hm,aes(last_supported_horizon,AUC_at_horizon),color="#D7301F",size=3)+
      geom_label(data=hm,aes(last_supported_horizon,AUC_at_horizon,label=paste0("up to ",format(last_supported_horizon,trim=TRUE)," y\n",cases_at_horizon," cases")),
        color="#D7301F",fill="white",label.size=.18,size=2.5,fontface="bold",nudge_y=.075)+
      facet_wrap(~method,ncol=2)+scale_x_continuous(breaks=c(.5,2,5,8,10,12,14,16))+scale_y_continuous(limits=c(.45,1))+
      labs(title=title,subtitle="Red marker = farthest consecutive horizon passing the prespecified criterion",
        x="Minimum years from blood draw to recorded diagnosis",y="Held-out AUC (95% CI)")+theme_5c(8)
  }
  pA<-curve_panel("Omic score","A. How many years ahead do scores discriminate?")
  pB<-curve_panel("Single biomarker","B. How many years ahead do individual biomarkers discriminate?")
  pC<-if(!nrow(hmax))blank_plot("C. Supported lead time",paste0("No method had lower 95% AUC CI > 0.50 with at least ",C5_LEAD_MIN_CASES," cases")) else
    ggplot(hmax,aes(last_supported_horizon,fct_reorder(method,last_supported_horizon)))+
      geom_segment(aes(x=0,xend=last_supported_horizon,yend=fct_reorder(method,last_supported_horizon)),color="grey75",linewidth=2.2)+
      geom_point(color="#D7301F",size=3.2)+geom_text(aes(label=paste0(" ",format(last_supported_horizon,trim=TRUE)," y; AUC ",sprintf("%.2f",AUC_at_horizon),"; n=",cases_at_horizon)),hjust=0,size=2.8,fontface="bold")+
      facet_grid(kind~.,scales="free_y",space="free_y")+coord_cartesian(xlim=c(0,max(hmax$last_supported_horizon)*1.55),clip="off")+
      labs(title="C. Headline lead time (red = supported, not maximum follow-up)",x="Years before recorded diagnosis",y=NULL)+theme_5c(9)
  counts<-lead|>group_by(horizon)|>summarise(cases=max(cases,na.rm=TRUE),controls=max(controls,na.rm=TRUE),.groups="drop")|>
    pivot_longer(c(cases,controls),names_to="sample",values_to="N")
  pD<-ggplot(counts,aes(horizon,N,color=sample))+geom_vline(xintercept=c(.5,1,2,5,10),linetype=3,color="grey82")+geom_line(linewidth=.9)+geom_point(size=1.5)+
    scale_color_manual(values=c(cases="#D7301F",controls="#2C7FB8"),labels=c(cases="Incident cases",controls="Eligible controls"))+
    labs(title="E. Information remaining at each lead time",subtitle=paste0("Headline requires lower 95% AUC CI > 0.50, at least ",C5_LEAD_MIN_CASES," cases, and no earlier failure"),x="Minimum lead time (years)",y="Participants",color=NULL)+theme_5c(9)+theme(legend.position="top")
  ww<-windows|>filter(kind=="Omic score",is.finite(AUC))
  pW<-if(!nrow(ww))blank_plot("D. Prespecified prediction windows","Window-specific AUC was unavailable")else ggplot(ww,aes(AUC,fct_reorder(method,AUC),color=window))+
    geom_errorbarh(aes(xmin=lo,xmax=hi),height=.15,position=position_dodge(width=.55))+geom_point(position=position_dodge(width=.55),size=2)+
    geom_vline(xintercept=.5,linetype=3,color="grey65")+scale_x_continuous(limits=c(.45,1))+
    labs(title="D. MASLD-style prespecified windows",subtitle="Within 5/10 years and diagnosis more than 10/12 years after baseline",x="Held-out AUC (95% CI)",y=NULL,color=NULL)+theme_5c(8)+theme(legend.position="top")
  (pA|pB)/plot_spacer()/(pC|pW|pD)+plot_layout(heights=c(1.35,.07,1))
}

topn_trajectory_summary <- function(rows,max_cases=C5_PAIR_MAX_CASES,step=C5_TEMPORAL_STEP){
  if(!nrow(rows))return(tibble())
  map_dfr(unique(rows$model),function(mm){
    d<-rows|>filter(model==mm,is.finite(score),is.finite(time),!is.na(event));ca<-d|>filter(event==1,time<=C5_TIME_CAP)
    if(!nrow(ca)||nrow(d)<2)return(tibble());set.seed(SEED+810+match(mm,unique(rows$model)))
    if(nrow(ca)>max_cases)ca<-ca|>slice_sample(n=max_cases)
    # Incidence-density sampling: a comparator may be event-free or may become
    # a case later, provided they are still under observation at the index
    # case's diagnosis time.  This estimates score separation within the same
    # risk set instead of contrasting only with never-cases.
    pairs<-map_dfr(seq_len(nrow(ca)),function(i){pool<-which(d$eid!=ca$eid[i]&d$time>=ca$time[i]);if(!length(pool))return(tibble());j<-sample(pool,1)
      tibble(years_before=ca$time[i],difference=ca$score[i]-d$score[j])})
    pairs|>mutate(bin=floor(years_before/step)*step+step/2)|>group_by(bin)|>
      summarise(mean_difference=mean(difference),se=sd(difference)/sqrt(n()),N=n(),.groups="drop")|>
      filter(N>=10)|>mutate(model=mm,lo=mean_difference-1.96*se,hi=mean_difference+1.96*se)
  })
}

plot_topn_mechanism <- function(topn){
  z<-topn$summary;ir<-topn$individual_range
  if(!nrow(z))return(blank_plot("Why combine biomarkers?","Top-N benchmark was unavailable"))
  pal<-c(`Top-1 biomarker (locked)`="#595959",`Unweighted top-N mean`="#2C7FB8",`Marginal-Cox-weighted top-N`="#D95F02")
  pA<-ggplot(z,aes(N,AUC,color=model,fill=model))+geom_hline(yintercept=.5,linetype=3,color="grey65")+
    {if(nrow(ir))geom_ribbon(data=ir,aes(N,ymin=AUC_q25,ymax=AUC_q75),inherit.aes=FALSE,fill="grey72",alpha=.35)else geom_blank()}+
    {if(nrow(ir))geom_line(data=ir,aes(N,AUC_max),inherit.aes=FALSE,color="grey45",linetype=2,linewidth=.65)else geom_blank()}+
    geom_ribbon(aes(ymin=AUC_lo,ymax=AUC_hi),alpha=.08,color=NA)+geom_line(linewidth=.9)+geom_point(size=1.45)+
    scale_color_manual(values=pal)+scale_fill_manual(values=pal)+labs(title="A. Held-out discrimination",
      subtitle="Grey band: IQR of the N individual biomarkers; dashed grey: their optimistic validation maximum",
      x="Training-ranked biomarkers included (N)",y="10-year IPCW AUC",color=NULL,fill=NULL)+theme_5c(9)+theme(legend.position="top")
  pB<-ggplot(z,aes(N,beta,color=model,group=model))+geom_hline(yintercept=0,linetype=3,color="grey65")+
    geom_ribbon(aes(ymin=beta_lo,ymax=beta_hi,fill=model),alpha=.08,color=NA)+geom_line(linewidth=.9)+geom_point(size=1.45)+
    scale_color_manual(values=pal)+scale_fill_manual(values=pal)+labs(title="B. Effect per training SD",
      subtitle="Validation Cox beta adjusted for the same clinical covariates",x="N",y="Log HR per 1-SD score",color=NULL,fill=NULL)+theme_5c(9)+theme(legend.position="none")
  sep<-z|>select(N,model,mean_separation,within_group_SD,cohen_d)|>
    pivot_longer(c(mean_separation,within_group_SD,cohen_d),names_to="metric",values_to="value")|>
    mutate(metric=recode(metric,mean_separation="Case-control mean separation",within_group_SD="Pooled within-group SD",cohen_d="Signal/noise (Cohen d)"))
  pC<-ggplot(sep,aes(N,value,color=model))+geom_line(linewidth=.85)+geom_point(size=1.3)+facet_wrap(~metric,ncol=1,scales="free_y")+
    scale_color_manual(values=pal)+labs(title="C. What changes when signals are combined?",
      subtitle="A score helps when separation grows relative to within-group noise; raw SD alone is not the mechanism",
      x="N",y=NULL,color=NULL)+theme_5c(8)+theme(legend.position="none")
  stab<-z|>select(N,model,AUC_boot_SD,separation_boot_SD)|>pivot_longer(c(AUC_boot_SD,separation_boot_SD),names_to="metric",values_to="value")|>filter(is.finite(value))|>
    mutate(metric=recode(metric,AUC_boot_SD="Bootstrap SD of AUC",separation_boot_SD="Bootstrap SD of mean separation"))
  pD<-if(!nrow(stab))blank_plot("D. Sampling stability","Bootstrap estimates unavailable")else ggplot(stab,aes(N,value,color=model))+geom_line(linewidth=.85)+geom_point(size=1.3)+facet_wrap(~metric,ncol=1,scales="free_y")+
    scale_color_manual(values=pal)+labs(title="D. Sampling stability",x="N",y="Bootstrap SD",color=NULL)+theme_5c(8)+theme(legend.position="none")
  tr<-topn_trajectory_summary(topn$rows)
  pE<-if(!nrow(tr))blank_plot("E. Diagnosis-anchored separation","Trajectory rows unavailable")else ggplot(tr,aes(bin,mean_difference,color=model,fill=model))+
    geom_hline(yintercept=0,linetype=3,color="grey60")+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.10,color=NA)+geom_line(linewidth=.95)+geom_point(size=1.3)+
    scale_color_manual(values=pal)+scale_fill_manual(values=pal)+scale_x_reverse(limits=c(C5_TIME_CAP,0),breaks=c(16,12,10,5,2,0))+
    labs(title=paste0("E. Score separation up to diagnosis (N=",unique(topn$rows$N)[1],")"),
      subtitle="Incidence-density risk-set matching; cross-sectional baselines in different people, not within-person trajectories",
      x="Years before the recorded diagnosis",y="Risk-set case − comparator score (95% CI)",color=NULL,fill=NULL)+theme_5c(9)+theme(legend.position="top")
  div<-z|>filter(model=="Marginal-Cox-weighted top-N")|>select(N,average_abs_correlation,effective_N)|>
    pivot_longer(c(average_abs_correlation,effective_N),names_to="metric",values_to="value")|>
    mutate(metric=recode(metric,average_abs_correlation="Average |correlation|",effective_N="Weight-based effective N"))
  pF<-ggplot(div,aes(N,value))+geom_line(color="#6A3D9A",linewidth=.9)+geom_point(color="#6A3D9A",size=1.5)+facet_wrap(~metric,ncol=1,scales="free_y")+
    labs(title="F. Redundancy versus diversification",subtitle="Highly correlated biomarkers add less independent information",x="N",y=NULL)+theme_5c(8)
  (pA|pB|pC)/plot_spacer()/(pD|pE|pF)+plot_layout(heights=c(1,.07,1))
}

attained_age_score_sensitivity <- function(score_rows,dat,methods,clinical,tvar,evar){
  age_covs<-unique(c("age",grep("^age($|[._])",clinical,value=TRUE,ignore.case=TRUE)));att_covs<-setdiff(clinical,age_covs)
  # Join by participant identifier rather than cached row position: an older
  # reusable score-model cache remains valid even if phenotype columns change.
  base_dat<-dat|>select(eid,all_of(unique(c(tvar,evar,".attained_entry",".attained_exit",clinical))))
  eligible_methods<-intersect(methods,unique(score_rows$method))
  if(!length(eligible_methods))return(tibble(method=character(),time_scale=character(),beta=numeric(),
    se=numeric(),N=integer(),events=integer(),lo=numeric(),hi=numeric()))
  map_dfr(eligible_methods,function(mm){
    d<-score_rows|>filter(method==mm,split=="validation",is.finite(score_z))|>
      select(eid,score_z)|>inner_join(base_dat,by="eid")
    fit_one<-function(kind){
      covs<-if(kind=="Baseline time")clinical else att_covs;surv<-if(kind=="Baseline time")paste0("Surv(",bt(tvar),",",bt(evar),")") else paste0("Surv(",bt(".attained_entry"),",",bt(".attained_exit"),",",bt(evar),")")
      rhs<-c("score_z",covs);dd<-d[,unique(c("score_z",covs,tvar,evar,".attained_entry",".attained_exit")),drop=FALSE];dd<-dd[complete.cases(dd),,drop=FALSE]
      if(kind!="Baseline time")dd<-dd[dd$.attained_exit>dd$.attained_entry,,drop=FALSE]
      fit<-if(nrow(dd)>=100&&sum(dd[[evar]]==1)>=20)tryCatch(coxph(as.formula(paste0(surv," ~ ",paste(bt(rhs),collapse=" + "))),dd,ties="efron"),error=function(e)NULL)else NULL
      sm<-if(is.null(fit))NULL else coef(summary(fit));if(is.null(sm)||!"score_z"%in%rownames(sm))return(tibble())
      tibble(method=mm,time_scale=kind,beta=sm["score_z","coef"],se=sm["score_z","se(coef)"],N=nrow(dd),events=sum(dd[[evar]]==1))
    }
    bind_rows(fit_one("Baseline time"),fit_one("Attained age, delayed entry"))
  })|>mutate(lo=beta-1.96*se,hi=beta+1.96*se)
}

plot_attained_age_score_sensitivity <- function(z){
  if(!nrow(z))return(blank_plot("Score time-scale sensitivity","No attained-age score model was estimable"))
  ggplot(z,aes(beta,fct_reorder(method,beta),color=time_scale))+geom_vline(xintercept=0,linetype=3,color="grey65")+
    geom_errorbarh(aes(xmin=lo,xmax=hi),height=.16,position=position_dodge(width=.45))+geom_point(position=position_dodge(width=.45),size=2.2)+
    scale_color_manual(values=c(`Baseline time`="#2C7FB8",`Attained age, delayed entry`="#D95F02"))+
    labs(title="Prediction-score sensitivity to the survival time scale",
      subtitle="Attained-age follow-up starts at baseline age; adult omics are not treated as birth measurements",
      x="Validation log HR per training-SD score",y=NULL,color=NULL)+theme_5c(9)+theme(legend.position="top")
}

plot_evidence_matrix <- function(ev){
  if(!nrow(ev))return(blank_plot("5C evidence matrix","No upstream evidence table was available"))
  grade_levels<-c("A","B","C","D")
  keep<-ev|>mutate(.grade=match(evidence_grade,grade_levels))|>
    arrange(.grade,desc(evidence_count),obs_p)|>slice_head(n=40)|>pull(feature)
  z0<-ev|>filter(feature%in%keep)|>transmute(feature,
    `C1 correlation`=C1,`Distal support`=as.integer(distal_support),`C2 MR`=C2_MR,`Reverse MR`=C2_REVERSE_MR,`C2 DPG`=C2_DANDELION,
    `C3 coloc`=C3,`C4 connection`=C4,`C5 stability`=ifelse(C5_available,C5,NA_integer_),
    `Reactive signal`=-as.integer(coalesce(reactive_compatible,FALSE)))|>
    pivot_longer(-feature,names_to="stage",values_to="support")
  grade<-ev|>filter(feature%in%keep)|>transmute(feature,stage="Grade",support=NA_integer_,grade_label=evidence_grade)
  z<-bind_rows(z0|>mutate(grade_label=NA_character_),grade)|>
    mutate(feature=factor(feature,levels=rev(keep)),
      stage=factor(stage,levels=c("C1 correlation","Distal support","C2 MR","Reverse MR","C2 DPG","C3 coloc","C4 connection","C5 stability","Reactive signal","Grade")))
  ggplot(z,aes(stage,feature,fill=factor(support)))+geom_tile(color="white",linewidth=.35)+
    geom_text(data=z|>filter(stage=="Grade"),aes(label=grade_label),fontface="bold",size=3.1)+
    scale_fill_manual(values=c(`-1`="#D9A441",`0`="grey92",`1`="#3B8064"),na.value="white",
      name="Evidence",labels=c(`-1`="Reactive-compatible",`0`="No",`1`="Yes"))+
    labs(title="5C evidence consolidation with causal safeguards",
      subtitle="Reactive evidence changes the biomarker role; it no longer vetoes concordant MR + robust colocalization",
      x=NULL,y=NULL)+theme_5c(8)+
    theme(axis.text.x=element_text(angle=35,hjust=1),panel.grid=element_blank(),legend.position="top")
}

plot_causal_reactive_map <- function(ev){
  if(!nrow(ev))return(blank_plot("Causal–reactive map","No evidence table"))
  z<-ev|>mutate(genetic_strength=coalesce(pmin(12,-log10(pmax(MR_p,1e-300))),0)+4*pmax(0,pmin(1,coalesce(PP.H4_robust_min,PP.H4,0))),
    reactive_strength=pmax(0,abs(coalesce(prevalent_score,0)))+pmax(0,abs(coalesce(duration_score,0))),
    label=ifelse(feature%in%c("GDF15","NTPROBNP","NPPB","PCSK9")|evidence_grade%in%c("A","B")&min_rank(obs_p)<=18,feature,NA_character_))
  ggplot(z,aes(genetic_strength,reactive_strength,color=consolidation_role))+
    geom_point(aes(size=pmin(12,-log10(pmax(obs_p,1e-300)))),alpha=.72)+
    geom_text_repel(aes(label=label),size=2.8,max.overlaps=25,seed=SEED,fontface="bold")+
    labs(title="Causal–reactive evidence map",subtitle="Upper-right biomarkers may be both causal and disease-responsive; neither axis alone proves mechanism",
      x="Genetic causal/locus evidence",y="Established-disease/reactive evidence",color=NULL,size="Incident\nevidence")+theme_5c(9)+theme(legend.position="bottom")
}

plot_performance_benchmark <- function(psum){
  if(!nrow(psum))return(blank_plot("Held-out prediction benchmark","No validation summary was available"))
  z<-psum|>select(biom_set,model,AUC,C_index)|>pivot_longer(c(AUC,C_index),names_to="metric",values_to="value")|>filter(is.finite(value))
  ggplot(z,aes(value,fct_reorder(biom_set,value,max),color=model,shape=model))+geom_vline(xintercept=.5,linetype=3,color="grey65")+
    geom_point(size=2.4,position=position_dodge(width=.38))+facet_wrap(~metric,scales="free_x")+
    scale_color_manual(values=c(`Model 0`="grey35",Biomarkers="#C77732",Combined="#287A58"),labels=c(`Model 0`="Clinical"))+
    labs(title="Held-out prediction benchmark",subtitle="All comparisons use the same outer validation split",x="Performance",y=NULL,color=NULL,shape=NULL)+theme_5c(9)+theme(legend.position="top")
}

plot_incremental_performance <- function(psum){
  if(!nrow(psum))return(blank_plot("Incremental discrimination","No validation summary was available"))
  z<-psum|>group_by(biom_set)|>summarise(clinical_AUC=AUC[match("Model 0",model)],biom_AUC=AUC[match("Biomarkers",model)],combined_AUC=AUC[match("Combined",model)],
      clinical_C=C_index[match("Model 0",model)],biom_C=C_index[match("Biomarkers",model)],combined_C=C_index[match("Combined",model)],.groups="drop")|>
    mutate(delta_AUC_combined=combined_AUC-clinical_AUC,delta_C_combined=combined_C-clinical_C,
      delta_AUC_biom=biom_AUC-clinical_AUC,delta_C_biom=biom_C-clinical_C)|>
    select(biom_set,starts_with("delta_"))|>pivot_longer(-biom_set,names_to="contrast",values_to="delta")|>filter(is.finite(delta))|>
    mutate(metric=ifelse(str_detect(contrast,"AUC"),"ΔAUC","ΔC-index"),model=ifelse(str_detect(contrast,"combined"),"Clinical + omics","Omics only"))
  if(!nrow(z))return(blank_plot("Incremental discrimination","Could not align clinical, omics, and combined models"))
  ggplot(z,aes(delta,fct_reorder(biom_set,delta,max),color=model,shape=model))+geom_vline(xintercept=0,color="grey60")+
    geom_point(size=2.5,position=position_dodge(width=.35))+facet_wrap(~metric,scales="free_x")+
    labs(title="Incremental discrimination over Model 0",x="Performance difference",y=NULL,color=NULL,shape=NULL)+theme_5c(9)+theme(legend.position="top")
}

plot_complexity_performance <- function(psum){
  z<-psum|>filter(model%in%c("Biomarkers","Combined"),is.finite(AUC),is.finite(n_selected),n_selected>0)
  if(!nrow(z))return(blank_plot("Parsimony versus discrimination","No finite feature-count/performance pairs"))
  z<-z|>mutate(label=ifelse(min_rank(desc(AUC))<=6|n_selected<=5,biom_set,NA_character_))
  ggplot(z,aes(n_selected,AUC,color=model))+geom_hline(yintercept=.5,linetype=3,color="grey65")+geom_point(size=2.6)+
    ggrepel::geom_text_repel(aes(label=label),size=2.8,max.overlaps=15,seed=SEED,show.legend=FALSE)+scale_x_log10()+
    labs(title="Model parsimony versus validation AUC",x="Selected biomarkers (log scale)",y="AUC",color=NULL)+theme_5c(10)+theme(legend.position="top")
}

score_correlation_table <- function(score_rows){
  w<-score_rows|>filter(split=="validation",is.finite(score_z))|>select(eid,method,score_z)|>distinct()|>pivot_wider(names_from=method,values_from=score_z)
  if(nrow(w)<20||ncol(w)<3)return(tibble())
  m<-suppressWarnings(cor(as.data.frame(w[,-1,drop=FALSE]),use="pairwise.complete.obs",method="spearman"))
  as.data.frame(as.table(m),stringsAsFactors=FALSE)|>as_tibble()|>rename(method1=Var1,method2=Var2,rho=Freq)
}

plot_score_correlation <- function(corr){
  if(!nrow(corr))return(blank_plot("Prediction-score concordance","Fewer than two validation scores were available"))
  ggplot(corr,aes(method1,method2,fill=rho))+geom_tile(color="white")+geom_text(aes(label=sprintf("%.2f",rho)),size=2.4)+
    scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,limits=c(-1,1),name="Spearman ρ")+
    labs(title="Concordance among prediction paradigms",x=NULL,y=NULL)+theme_5c(8)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
}

subgroup_auc_table <- function(score_rows,dat,methods){
  keep<-intersect(c("eid","age","sex"),names(dat));if(!all(c("eid","age","sex")%in%keep))return(tibble())
  if(!length(intersect(methods,unique(score_rows$method))))return(tibble())
  dd<-score_rows|>filter(method%in%methods,split=="validation",is.finite(score_z))|>left_join(dat|>select(all_of(keep)),by="eid")
  agecuts<-quantile(dd$age,c(0,1/3,2/3,1),na.rm=TRUE);if(anyDuplicated(agecuts))return(tibble())
  dd<-dd|>mutate(age_group=cut(age,agecuts,include.lowest=TRUE,dig.lab=3),sex_group=paste0("Sex category ",sex))
  bind_rows(dd|>mutate(domain="Age tertile",subgroup=as.character(age_group)),dd|>mutate(domain="Reported sex",subgroup=sex_group))|>
    group_by(method,domain,subgroup)|>group_modify(~{
      if(nrow(.x)<100||sum(.x$event==1)<15||sum(.x$event==0)<30)return(tibble())
      auc<-weighted_time_auc(.x$time,.x$event,.x$score_z,HORIZON);if(!is.finite(auc))return(tibble())
      set.seed(SEED+nrow(.x));bs<-replicate(100,{ii<-sample(seq_len(nrow(.x)),nrow(.x),replace=TRUE);weighted_time_auc(.x$time[ii],.x$event[ii],.x$score_z[ii],HORIZON)})
      ci<-quantile(bs,c(.025,.975),na.rm=TRUE);tibble(N=nrow(.x),events=sum(.x$event==1),AUC=auc,lo=ci[[1]],hi=ci[[2]])
    })|>ungroup()
}

plot_subgroup_auc <- function(sg){
  if(!nrow(sg))return(blank_plot("Subgroup discrimination","Age and sex subgroup estimates were unavailable"))
  ggplot(sg,aes(AUC,interaction(method,subgroup,lex.order=TRUE),color=method))+geom_vline(xintercept=.5,linetype=3,color="grey65")+
    geom_errorbarh(aes(xmin=lo,xmax=hi),height=.12)+geom_point(size=2)+facet_wrap(~domain,scales="free_y",ncol=1)+
    labs(title="Validation discrimination across prespecified subgroups",subtitle=paste0(HORIZON,"-year IPCW AUC; reported sex is retained in its source coding"),x="AUC (bootstrap 95% CI)",y=NULL,color=NULL)+theme_5c(8)+theme(legend.position="top")
}


# 🚩 Evidence table: C2 includes MR and DANDELION, C3 coloc, C4 connection.
make_evidence_table <- function(layer,outdir,stability=tibble()){
  rr<-function(job){
    f<-file.path(le8_job_dir(outdir,job),paste0(sub("_.*$","",job),".res.rds"))
    if(!file.exists(f))return(list(data=list(),available=FALSE,status="missing"))
    z<-tryCatch(readRDS(f),error=function(e)e)
    if(inherits(z,"condition"))return(list(data=list(),available=FALSE,status=paste0("unreadable: ",conditionMessage(z))))
    declared<-str_to_lower(as.character(z$meta$status%||%"ok")[[1]])
    usable<-!declared%in%c("unavailable","no qtl locus constructed","failed","not run")
    list(data=if(usable)z else list(),available=usable,
      status=if(usable)"available"else paste0("available file; ",declared))
  }
  u1<-rr("c1_correlate");u2<-rr("c2_cause");u3<-rr("c3_coloc");u4<-rr("c4_connect")
  c1<-u1$data;c2<-u2$data;c3<-u3$data;c4<-u4$data
  a0<-c1$association%||%tibble();a<-if(nrow(a0))a0|>transmute(feature=term,obs_beta=beta,obs_p=p.value,obs_FDR=FDR)else tibble(feature=character())
  g0<-c1$pgs_incident%||%tibble();g<-if(nrow(g0))g0|>transmute(feature=term,
    pgs_beta=beta,pgs_p=p.value,pgs_FDR=FDR,pgs_N=N_total)else tibble(feature=character())
  di0<-c1$directionality%||%tibble();di<-if(nrow(di0)&&all(c("term","direction_class")%in%names(di0)))di0|>
    transmute(feature=term,temporal_class=direction_class,reactive_compatible=coalesce(reactive_compatible,FALSE),
      distal_support=coalesce(distal_support,FALSE),FDR_landmark5=FDR_landmark5,
      prevalent_score=prevalent_score,duration_score=duration_score,landmark5_score=landmark5_score) else
    tibble(feature=character(),temporal_class=character(),reactive_compatible=logical(),distal_support=logical(),FDR_landmark5=numeric(),prevalent_score=numeric(),duration_score=numeric(),landmark5_score=numeric())
  m0<-c2$MR%||%tibble();m<-if(nrow(m0))m0|>filter(is.finite(pval))|>group_by(exposure)|>slice_min(pval,n=1,with_ties=FALSE)|>ungroup()|>transmute(feature=exposure,MR_analysis=analysis,MR_beta=b,MR_p=pval,MR_FDR=FDR_all)else tibble(feature=character())
  r0<-c2$MR_reverse%||%tibble();rmr<-if(nrow(r0))r0|>filter(is.finite(pval))|>transmute(feature,reverse_MR_beta=b,reverse_MR_p=pval,reverse_MR_FDR=FDR_reverse)else tibble(feature=character())
  dan<-c2$DANDELION%||%list();d0<-as_tibble(dan$targets%||%tibble())
  if(nrow(d0)&&!"consolidation_eligible"%in%names(d0)){
    primary_input<-if("gene_evidence_type"%in%names(d0))
      !str_detect(str_to_lower(coalesce(d0$gene_evidence_type,"")),"magma|adapted") else FALSE
    tested_n<-suppressWarnings(as.numeric((dan$ptrans_dimensions%||%c(NA_real_))[[1]]))
    frac<-if(is.finite(tested_n)&&tested_n>0)nrow(d0)/tested_n else NA_real_
    d0$consolidation_eligible<-primary_input&!(is.finite(frac)&&frac>.25)
  }
  dd<-if(nrow(d0))d0|>transmute(feature=gene2,DANDELION_p,DANDELION_loci=n_distal_loci,
      DANDELION_sensitivity=TRUE,DANDELION_primary=as.logical(consolidation_eligible))else
    tibble(feature=character(),DANDELION_p=numeric(),DANDELION_loci=integer(),DANDELION_sensitivity=logical(),DANDELION_primary=logical())
  co0<-as_tibble(c3$summary%||%tibble())
  if(nrow(co0)&&!"PP.H4_robust_min"%in%names(co0))co0$PP.H4_robust_min<-co0$PP.H4
  co<-if(nrow(co0)&&all(c("feature","PP.H4")%in%names(co0)))co0|>filter(status=="ok")|>
    mutate(.robust=coalesce(PP.H4_robust_min,PP.H4))|>group_by(feature)|>slice_max(.robust,n=1,with_ties=FALSE)|>ungroup()|>
    select(feature,PP.H4,PP.H4_robust_min,everything())else tibble(feature=character())
  cx0<-c4$membership%||%tibble();cx<-if(nrow(cx0)&&"feature"%in%names(cx0))cx0 else tibble(feature=character())
  st<-if(nrow(stability)&&all(c("feature","selection_frequency")%in%names(stability)))stability|>group_by(feature)|>summarise(max_selection_frequency=max(selection_frequency,na.rm=TRUE),.groups="drop")else tibble(feature=character(),max_selection_frequency=numeric())
  z<-Reduce(function(x,y)full_join(x,y,by="feature"),list(a,g,di,m,rmr,dd,co,cx,st));if(!nrow(z))return(tibble())
  defaults<-list(obs_FDR=NA_real_,obs_p=NA_real_,MR_p=NA_real_,MR_FDR=NA_real_,DANDELION_sensitivity=FALSE,
    DANDELION_primary=FALSE,PP.H4=NA_real_,PP.H4_robust_min=NA_real_,strict_YS=FALSE,
    max_selection_frequency=NA_real_,reactive_compatible=FALSE,distal_support=FALSE,reverse_MR_FDR=NA_real_,
    prevalent_score=NA_real_,duration_score=NA_real_,landmark5_score=NA_real_,
    pgs_beta=NA_real_,pgs_p=NA_real_,pgs_FDR=NA_real_,pgs_N=NA_real_)
  for(nm in names(defaults))if(!nm%in%names(z))z[[nm]]<-defaults[[nm]]
  z|>mutate(C1=as.integer(is.finite(obs_FDR)&obs_FDR<.05),
    PGS_tested=is.finite(pgs_p),C1_PGS=ifelse(PGS_tested,as.integer(pgs_FDR<.05),NA_integer_),
    PGS_observed_sign_match=PGS_tested&is.finite(obs_beta)&sign(pgs_beta)==sign(obs_beta),
    C2_MR=if(u2$available)as.integer(is.finite(MR_FDR)&MR_FDR<.05)else NA_integer_,
    C2_REVERSE_MR=if(u2$available)as.integer(is.finite(reverse_MR_FDR)&reverse_MR_FDR<.05)else NA_integer_,
    C2_DANDELION=if(u2$available)as.integer(coalesce(DANDELION_primary,FALSE))else NA_integer_,
    C2_DANDELION_sensitivity=if(u2$available)as.integer(coalesce(DANDELION_sensitivity,FALSE))else NA_integer_,
    C2=if(u2$available)as.integer(coalesce(C2_MR,0L)==1L|coalesce(C2_DANDELION,0L)==1L)else NA_integer_,
    C3=if(u3$available)as.integer(is.finite(coalesce(PP.H4_robust_min,PP.H4))&coalesce(PP.H4_robust_min,PP.H4)>=.7)else NA_integer_,
    C4=if(u4$available)as.integer(coalesce(strict_YS,FALSE))else NA_integer_,
    C5=as.integer(is.finite(max_selection_frequency)&max_selection_frequency>=.6),
    C5_available=is.finite(max_selection_frequency),reactive_compatible=coalesce(reactive_compatible,FALSE),
    distal_support=coalesce(distal_support,FALSE),temporal_warning=reactive_compatible,
    C2_status=u2$status,C3_status=u3$status,C4_status=u4$status,
    evidence_count=rowSums(cbind(C1,C2,C3,C4,ifelse(C5_available,C5,NA_integer_)),na.rm=TRUE),
    causal_locus_evidence=coalesce(C2_MR,0L)==1L|coalesce(C2_DANDELION,0L)==1L|coalesce(C3,0L)==1L,
    evidence_grade=case_when(
      coalesce(C2_MR,0L)==1L&coalesce(C3,0L)==1L&C1==1~"A",
      causal_locus_evidence&(C1==1|distal_support)~"B",
      evidence_count>=2~"C",TRUE~"D"),
    consolidation_role=case_when(evidence_grade%in%c("A","B")&reactive_compatible~"Causal + reactive mixed biomarker",
      evidence_grade=="A"~"High-priority causal candidate",
      coalesce(C2_REVERSE_MR,0L)==1L&!causal_locus_evidence~"Disease-liability-responsive biomarker",
      distal_support&!causal_locus_evidence~"Distal predictive biomarker",
      reactive_compatible&!causal_locus_evidence~"Reactive/diagnostic-compatible biomarker",
      evidence_grade=="B"~"Multi-domain causal/locus candidate",
      evidence_grade=="C"~"Multi-domain supported candidate",TRUE~"Single-domain or insufficient"),
    # This is an auditable prior for transparent penalized prediction, not a
    # fitted causal effect. Reactive-compatible features without locus evidence
    # receive the smallest weight; a name such as GDF15 is never hard-coded.
    evidence_weight_prior=case_when(
      causal_locus_evidence&distal_support&!reactive_compatible~1.00,
      causal_locus_evidence&reactive_compatible~.75,
      causal_locus_evidence~.90,
      coalesce(C1_PGS,0L)==1L&distal_support~.80,
      distal_support~.65,
      reactive_compatible~.25,
      C1==1~.50,TRUE~.35),
    prediction_role=case_when(
      causal_locus_evidence~"Mechanism-supported candidate",
      distal_support~"Distal prediction candidate",
      reactive_compatible~"Downweight; diagnostic/reactive-compatible",
      coalesce(C1_PGS,0L)==1L~"Inherited-score candidate; audit pleiotropy",
      TRUE~"Exploratory"))|>
    arrange(factor(evidence_grade,levels=c("A","B","C","D")),desc(evidence_count),obs_p)
}


# 🚩 Main C5
run_c5_layer <- function(layer=c("protein","metabolite")){
  # LE8_REVISION_UPDATES: guarded caches, definitions and additions.
  layer <- match.arg(layer)
  le8_load_revision(layer, "c5_consolidate")
  .le8_review_env <- environment()
  on.exit(le8_finish_revision(layer, "c5_consolidate", .le8_review_env), add=TRUE)

  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir);rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
  module_file<-function(job)file.path(le8_job_dir(outdir,job),paste0(sub("_.*$","",job),".res.rds"))
  safe_module<-function(job,required=FALSE){
    f<-module_file(job)
    if(!file.exists(f)){
      if(required)stop("C5 requires C1 for the same layer; missing: ",f,call.=FALSE)
      return(list(data=list(),available=FALSE,status="missing",file=f))
    }
    z<-tryCatch(readRDS(f),error=function(e)e)
    if(inherits(z,"condition")){
      if(required)stop("C5 could not read required C1 result: ",conditionMessage(z),call.=FALSE)
      warning("C5: optional ",job," result is unreadable; corresponding panels will be unavailable: ",conditionMessage(z),call.=FALSE)
      return(list(data=list(),available=FALSE,status="unreadable",file=f))
    }
    declared<-str_to_lower(as.character(z$meta$status%||%"ok")[[1]])
    usable<-!declared%in%c("unavailable","no qtl locus constructed","failed","not run")
    if(!usable)return(list(data=list(),available=FALSE,status=paste0("available file; ",declared),file=f))
    list(data=z,available=TRUE,status="available",file=f)
  }
  u1<-safe_module("c1_correlate",TRUE);u2<-safe_module("c2_cause");u3<-safe_module("c3_coloc");u4<-safe_module("c4_connect")
  upstream_signature<-paste0("c1=",u1$status,";c2=",u2$status,";c3=",u3$status,";c4=",u4$status)
  upstream_audit<-tibble(module=c("C1","C2","C3","C4"),
    required=c(TRUE,FALSE,FALSE,FALSE),available=c(u1$available,u2$available,u3$available,u4$available),
    status=c(u1$status,u2$status,u3$status,u4$status),file=c(u1$file,u2$file,u3$file,u4$file))
  write_raw_csv(upstream_audit,"c5.upstream_availability.csv",rawdir)
  message("C5/",layer,": ",upstream_signature)
  final_cache<-file.path(rawdir,"c5.res.rds")
  if(cache_valid(final_cache)){
    old<-tryCatch(readRDS(final_cache),error=function(e)NULL)
    if(!is.null(old)&&identical(old$meta$code_version%||%NA_character_,C5_CODE_VERSION)&&
       identical(old$meta$upstream_signature%||%NA_character_,upstream_signature)){
      cache_message(paste0("C5/",layer),final_cache);return(old)
    }
    message("C5/",layer,": cache predates ",C5_CODE_VERSION,"; recomputing score coverage and IPCW figures")
  }
  all_glmnet_label<-if(layer=="protein")"Pradeep-style / glmnet"else"All-metabolite / glmnet"
  lightgbm_label<-if(layer=="protein")"Yu-style / LightGBM"else"MWAS-ranked / LightGBM"
  biom<-if(layer=="protein")read_prot()else read_met();biom_vars<-setdiff(names(biom),"eid")
  incfile<-if(layer=="protein")C5_INC_PROT else C5_INC_MET
  if(nzchar(incfile)&&file.exists(incfile)){req<-unique(scan(incfile,what="character",quiet=TRUE));biom_vars<-intersect(biom_vars,req);biom<-biom[,c("eid",biom_vars),drop=FALSE]}
  all0<-read_all();prs_vars<-find_prs_vars(all0,Y,4)
  need<-unique(c("eid","ethnic.c",vars.basic,vars.le8,prs_vars,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
  dat_all<-all0[,intersect(need,names(all0)),drop=FALSE]|>filter_analysis_cohort()|>make_outcome(Y)|>add_attained_age_time(Y)|>inner_join(biom,by="eid");rm(all0,biom);gc()
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");bvar<-paste0(Y,".b2e");bivar<-paste0(Y,".bi2e")
  clinical<-intersect(covs_use,names(dat_all));biom_vars<-intersect(biom_vars,names(dat_all))
  prevalent_dat<-dat_all|>filter(is.finite(.data[[bvar]]),.data[[bvar]]<=0)
  dat<-dat_all[complete.cases(dat_all[,c(tvar,evar),drop=FALSE])&dat_all[[tvar]]>0,,drop=FALSE];rm(dat_all);invisible(gc())
  split<-make_outer_split(dat,evar=evar);dat$c5_split<-split
  message("C5/",layer,": N=",nrow(dat),", events=",sum(dat[[evar]]==1),", training=",sum(split=="training"),", validation=",sum(split=="validation"),", biomarkers=",length(biom_vars))

  c1f<-u1$file;c1<-u1$data;a1<-c1$association%||%tibble()
  if(!nrow(a1))stop("C5 requires a non-empty C1 association table for the same layer: ",c1f,call.=FALSE)
  # Re-screen C1 associations inside the outer training set.  The full-cohort
  # C1 table remains part of 5C evidence consolidation but is not allowed to
  # choose features for held-out prediction.
  screen_cache<-file.path(rawdir,"c5.training_association_screen.rds");training_screen<-read_stage_cache(screen_cache)
  if(is.null(training_screen)){
    message("C5/",layer,": training-only PWAS/MWAS screen")
    training_screen<-cox_scan(dat[split=="training",,drop=FALSE],biom_vars,clinical,Y,time_var=tvar,event_var=evar)
    write_stage_cache(training_screen,screen_cache)
  }else message("C5/",layer,": reuse training-only association screen")
  write_raw_csv(training_screen,"c5.training_association_screen.csv",rawdir)
  ranked<-if(nrow(training_screen))training_screen|>filter(is.finite(p.value))|>arrange(p.value)|>pull(term)|>intersect(biom_vars)else biom_vars
  pwas<-if(nrow(training_screen))training_screen|>filter(is.finite(FDR),FDR<.05)|>arrange(p.value)|>pull(term)|>intersect(biom_vars)else character()
  c2f<-u2$file;c2<-u2$data;m2<-c2$MR%||%tibble()
  mrset<-if(nrow(m2)&&truthy(Sys.getenv("C5_GENETIC_EVIDENCE_INDEPENDENT",unset="FALSE")))m2|>filter(is.finite(FDR_all),FDR_all<.05)|>arrange(pval)|>pull(exposure)|>unique()|>intersect(biom_vars)else character()
  local_class<-if(layer=="protein")"cis"else"local"
  local_mrset<-if(nrow(m2))m2|>filter(analysis==local_class,is.finite(FDR_all),FDR_all<.05)|>arrange(pval)|>pull(exposure)|>unique()|>intersect(biom_vars)else character()
  # Training-only five-year landmark screen: this is the predictive component
  # least compatible with a purely imminent-diagnosis signal.
  distal_cache<-file.path(rawdir,paste0("c5.distal_landmark_screen.",C5_CODE_VERSION,".rds"));distal_screen<-read_stage_cache(distal_cache)
  if(is.null(distal_screen)){
    distal_candidates<-head(ranked,C5_DISTAL_TOP);td<-dat[split=="training"&is.finite(dat[[tvar]])&dat[[tvar]]>C5_DISTAL_LANDMARK,,drop=FALSE]
    if(nrow(td)>=300&&sum(td[[evar]]==1,na.rm=TRUE)>=30&&length(distal_candidates)){
      td$.distal_time<-td[[tvar]]-C5_DISTAL_LANDMARK;td$.distal_event<-td[[evar]]
      distal_screen<-cox_scan(td,distal_candidates,clinical,Y,time_var=".distal_time",event_var=".distal_event")|>mutate(landmark_years=C5_DISTAL_LANDMARK)
    }else distal_screen<-tibble(term=character(),beta=numeric(),p.value=numeric(),FDR=numeric(),landmark_years=numeric())
    write_stage_cache(distal_screen,distal_cache)
  }
  write_raw_csv(distal_screen,"c5.training_distal_landmark_screen.csv",rawdir)
  distalset<-distal_screen|>filter(is.finite(FDR),FDR<.05)|>arrange(p.value)|>pull(term)|>intersect(biom_vars)
  c3f<-u3$file;c3<-u3$data;co3<-as_tibble(c3$summary%||%tibble())
  if(nrow(co3)&&!"PP.H4_robust_min"%in%names(co3))co3$PP.H4_robust_min<-co3$PP.H4
  if(nrow(co3)&&!"status"%in%names(co3))co3$status<-"ok"
  colocset<-if(nrow(co3))co3|>filter(status=="ok",coalesce(PP.H4_robust_min,PP.H4)>=.7)|>pull(feature)|>unique()|>intersect(biom_vars)else character()
  # The causal score is deliberately stricter than a union of all MR and all
  # coloc hits: require local/cis MR and robust colocalization for the same
  # feature. Distal/trans-only MR remains in the separate MR score.
  geneticset<-if(u2$available&&u3$available&&truthy(Sys.getenv("C5_GENETIC_EVIDENCE_INDEPENDENT",unset="FALSE"))) le8_same_locus_candidates(m2,co3,layer) else character()
  hybridset<-if(u2$available&&u3$available&&length(geneticset)&&length(distalset))unique(c(distalset,geneticset))else character()
  # A compact, transparent alternative to a black-box all-omic score. Candidate
  # selection uses only the outer training screen plus genetic/locus evidence.
  # Lower glmnet penalty factors favor distal and causal/locus-supported features;
  # training-ranked features lacking either support are retained but downweighted.
  mechanismset<-head(intersect(unique(c(geneticset,distalset,head(ranked,C5_MECHANISM_N))),biom_vars),C5_MECHANISM_N)
  mechanism_prior<-tibble(feature=mechanismset,
    distal_training=feature%in%distalset,genetic_locus=feature%in%geneticset)|>
    mutate(prior_class=case_when(genetic_locus&distal_training~"genetic + distal",
      genetic_locus~"genetic/locus",distal_training~"distal training signal",TRUE~"training-ranked only"),
      penalty_factor=1.0,
      interpretation="Uniform glmnet penalty; evidence defines candidates, not hard-coded biomarker weights")
  mechanism_pf<-setNames(mechanism_prior$penalty_factor,mechanism_prior$feature)
  write_raw_csv(mechanism_prior,"c5.mechanism_weight_prior.csv",rawdir)
  # For this analysis, "User specified" means the ten strongest training-screen biomarkers.
  user<-head(ranked,10L)
  yuvars<-if(nzchar(C5_YU_LIST)&&file.exists(C5_YU_LIST)) intersect(unique(scan(C5_YU_LIST,what="character",quiet=TRUE)),biom_vars) else head(ranked,C5_YU_PRESELECT_N)
  if(!nzchar(C5_YU_LIST)&&layer=="protein") message("C5: C5_YU_LIST not supplied; Yu column is a Yu-style 257-protein proxy using the C1 ranking, not an exact reproduction of the published 257-protein panel.")
  if(!nzchar(C5_YU_LIST)&&layer=="metabolite") message("C5: metabolite LightGBM uses the C1 MWAS ranking (up to C5_YU_PRESELECT_N available metabolites).")

  # Nature Aging-style sequential screen, learned only inside the outer training sample.
  seq0<-sequential_forward_panel(dat[split=="training",,drop=FALSE],pwas,tvar,evar,HORIZON,SEED+71,max_extended=C5_EXTENDED_N)
  pars<-intersect(seq0$parsimonious,biom_vars);extended<-intersect(seq0$extended,biom_vars)
  write_raw_csv(seq0$log,"c5.sequential_forward_training.csv",rawdir)

  # C4 connection-defined candidate sets for Fig3.
  c4f<-u4$file;c4<-u4$data;l4<-c4$lists%||%list()
  required_c4_lists<-c("NS","YS_all","YSP_all")
  c4_contract_ok<-u4$available&&all(required_c4_lists%in%names(l4))
  if(u4$available&&!c4_contract_ok)warning("C5: C4 result lacks the final list contract; C4 panels will be unavailable.",call.=FALSE)
  if(!c4_contract_ok)l4<-list(NS=character(),YS_all=character(),YSP_all=character())
  ys4<-if(c4_contract_ok)le8_training_connection_set(dat[split=="training",,drop=FALSE],biom_vars,intersect(vars.le8,names(dat)),intersect(vars.basic,names(dat)),rawdir)else character()
  # YS is selected without using disease outcome. Rebuild the outcome-ranked
  # NS and YSP-plus portions from the outer training screen to avoid leakage.
  ns4<-if(c4_contract_ok)head(ranked,as.integer(le8_num_env("C5_CONNECTION_NS_N",160)))else character()
  plus_n<-if(c4_contract_ok)as.integer(le8_num_env("C5_CONNECTION_PLUS_N",80))else 0L
  plus4<-if(plus_n>0L)head(setdiff(ranked,ys4),plus_n)else character()
  c4sets<-list(`C4 NS`=ns4,`C4 YS`=ys4,`C4 YSplus`=unique(c(ys4,plus4)))
  write_raw_csv(bind_rows(imap(c4sets,~tibble(set=.y,feature=.x))),"c5.connection_sets_training.csv",rawdir)
  empty_c4<-names(c4sets)[lengths(c4sets)==0L]
  if(length(empty_c4))message("C5: optional C4 score panels unavailable for: ",paste(empty_c4,collapse=", "))

  module_tag<-paste0("c2",as.integer(u2$available),"c3",as.integer(u3$available),"c4",as.integer(c4_contract_ok))
  cache<-file.path(rawdir,paste0("c5.score_models.",C5_MODEL_VERSION,".",module_tag,".rds"));sc<-if(cache_valid(cache))tryCatch(readRDS(cache),error=function(e)NULL)else NULL
  if(is.null(sc)){
    message("C5: fit ",all_glmnet_label," all-biomarker score")
    pr<-fit_glmnet_score(dat,biom_vars,tvar,evar,split,all_glmnet_label,lambda_rule="lambda.1se")
    message("C5: fit ",lightgbm_label," on ",length(yuvars)," C1-ranked biomarkers")
    yu<-fit_lightgbm_score(dat,yuvars,tvar,evar,split,lightgbm_label)
    us<-fit_glmnet_score(dat,user,tvar,evar,split,"User specified",lambda_rule="lambda.1se")
    pw<-fit_glmnet_score(dat,pwas,tvar,evar,split,"PWAS/MWAS significant",lambda_rule="lambda.1se")
    mr<-fit_glmnet_score(dat,mrset,tvar,evar,split,"MR significant",lambda_rule="lambda.1se")
    ds<-fit_glmnet_score(dat,distalset,tvar,evar,split,"Distal antecedent",lambda_rule="lambda.1se")
    gc<-fit_glmnet_score(dat,geneticset,tvar,evar,split,"Genetic-region evidence",lambda_rule="lambda.1se")
    hy<-fit_glmnet_score(dat,hybridset,tvar,evar,split,"Hybrid triangulated",lambda_rule="lambda.1se")
    mw<-fit_glmnet_score(dat,mechanismset,tvar,evar,split,"Evidence-selected compact",
      lambda_rule="lambda.1se",penalty_factor=mechanism_pf)
    pa<-fit_glmnet_score(dat,pars,tvar,evar,split,"Parsimonious",lambda_rule="lambda.1se")
    ex<-fit_glmnet_score(dat,extended,tvar,evar,split,"Extended",lambda_rule="lambda.1se")
    ns<-fit_glmnet_score(dat,c4sets[["C4 NS"]],tvar,evar,split,"C4 NS",lambda_rule="lambda.1se")
    ys<-fit_glmnet_score(dat,c4sets[["C4 YS"]],tvar,evar,split,"C4 YS",lambda_rule="lambda.1se")
    yp<-fit_glmnet_score(dat,c4sets[["C4 YSplus"]],tvar,evar,split,"C4 YSplus",lambda_rule="lambda.1se")
    gp<-fit_glmnet_score(dat,prs_vars,tvar,evar,split,"Genetic PRS",lambda_rule="lambda.min")
    clinical_obj<-fit_glmnet_score(dat,clinical,tvar,evar,split,"Clinical",lambda_rule="lambda.min")
    objs<-list(Pradeep=pr,Yu=yu,User=us,PWAS=pw,MR=mr,Distal=ds,Genetic=gc,Hybrid=hy,
      Mechanism= mw,Parsimonious=pa,Extended=ex,C4_NS=ns,C4_YS=ys,C4_YSplus=yp,PRS=gp)
    score_rows<-standardize_scores(bind_rows(lapply(Filter(Negate(is.null),objs),`[[`,"rows")))
    sc<-c(objs,list(Clinical=clinical_obj,rows=score_rows,split=split,
      sets=list(PWAS=pwas,MR=mrset,Distal=distalset,Genetic=geneticset,Hybrid=hybridset,
        Mechanism=mechanismset,Parsimonious=pars,Extended=extended,Yu=yuvars,User=user,C4=c4sets)))
    saveRDS(sc,cache,compress="xz")
  } else {message("C5: reuse score-model cache");score_rows<-sc$rows;clinical_obj<-sc$Clinical}

  object_map<-setNames(list(sc$Pradeep,sc$Yu,sc$User,sc$PWAS,sc$MR,sc$Distal,sc$Genetic,sc$Hybrid,sc$Mechanism,sc$Parsimonious,sc$Extended,
    sc$C4_NS,sc$C4_YS,sc$C4_YSplus,sc$PRS),
    c(all_glmnet_label,lightgbm_label,"User specified","PWAS/MWAS significant","MR significant","Distal antecedent","Genetic-region evidence","Hybrid triangulated","Evidence-selected compact","Parsimonious","Extended",
      "C4 NS","C4 YS","C4 YSplus","Genetic PRS"))
  base_score_rows<-score_rows|>filter(split!="prevalent")
  prev_rows<-imap_dfr(object_map,function(obj,nm)predict_prevalent_rows(obj,prevalent_dat,nm,bvar,base_score_rows))
  score_rows<-bind_rows(base_score_rows,prev_rows)
  score_rows$kind<-"Omic score"
  score_audit<-score_coverage_audit(score_rows,names(object_map),nrow(prevalent_dat))
  bad_prev<-score_audit|>filter(split=="prevalent",status!="ok")
  if(nrow(bad_prev))warning("C5 prevalent scoring coverage issue for: ",paste(paste0(bad_prev$method," [",bad_prev$status,"]"),collapse="; "),call.=FALSE)
  write_raw_csv(score_audit,"c5.score_coverage_audit.csv",rawdir)
  write_raw_csv(score_rows,"c5.person_scores.csv",rawdir)
  prior<-NULL
  pred<-prior$prediction%||%tibble()
  if(!nrow(pred)){
    bundles<-imap(object_map,function(obj,nm)if(is.null(obj))tibble()else make_validation_bundle(dat,obj,clinical_obj,clinical,tvar,evar,split,nm))
    pred<-bind_rows(bundles)
  }
  if(!nrow(pred))pred<-empty_prediction_rows()
  write_raw_csv(pred,"c5.prediction_validation_rows.csv",rawdir)
  psum<-if(nrow(pred))pred|>group_by(biom_set,model)|>summarise(N=n(),events=sum(event),C_index=first(cindex),
    AUC=weighted_time_auc(time,event,score,HORIZON),AUC_horizon_years=HORIZON,
    AUC_definition="IPCW cumulative/dynamic",n_selected=median(n_selected,na.rm=TRUE),.groups="drop") else
    tibble(biom_set=character(),model=character(),N=integer(),events=integer(),C_index=numeric(),
      AUC=numeric(),AUC_horizon_years=numeric(),AUC_definition=character(),n_selected=numeric())
  write_raw_csv(psum,"c5.prediction_summary.csv",rawdir)

  method_names<-unique(score_rows$method)
  fig1_order<-c(all_glmnet_label,lightgbm_label,"User specified")
  fig2_order<-c("Distal antecedent","Genetic-region evidence","Hybrid triangulated","Evidence-selected compact")
  fig3_order<-c("C4 NS","C4 YS","C4 YSplus")
  # Write the primary panels before optional matching/bootstrap diagnostics.
  # If a very large metabolite run is interrupted later, the principal C5
  # figures still exist and the output log identifies the unfinished step.
  save_plot(c5_make_5row_grid(fig1_order,pred,score_rows,list(),list(),"Prediction paradigms"),"c5.Fig1.simple_pred.png",24,13.5,outdir=outdir)
  save_plot(c5_make_5row_grid(fig2_order,pred,score_rows,list(),list(),"Evidence-screened prediction"),"c5.Fig2.screened_pred.png",24,13.5,outdir=outdir)
  save_plot(c5_make_5row_grid(fig3_order,pred,score_rows,list(),list(),"C4 connection-guided prediction"),"c5.Fig3.connection_pred.png",24,13.5,outdir=outdir)
  lead_methods<-intersect(c(all_glmnet_label,lightgbm_label,"Distal antecedent","Genetic-region evidence","Hybrid triangulated","Evidence-selected compact","C4 YSplus"),method_names)
  individual_features<-unique(c(head(ranked,3L),intersect(C5_LEAD_ANCHORS,biom_vars)))|>head(6L)
  individual_lead_rows<-make_individual_lead_rows(dat,individual_features,training_screen,split,tvar,evar)
  lead_input<-bind_rows(score_rows|>filter(method%in%lead_methods),individual_lead_rows)
  lead<-leadtime_auc_table(lead_input,unique(lead_input$method),C5_TIME_CAP)
  lead_headline<-leadtime_headline(lead,C5_LEAD_MIN_CASES)
  lead_windows<-prediction_window_auc_table(lead_input,unique(lead_input$method))
  write_raw_csv(lead,"c5.leadtime_discrimination.csv",rawdir);write_raw_csv(lead_headline,"c5.leadtime_headline.csv",rawdir);write_raw_csv(lead_windows,"c5.leadtime_prediction_windows.csv",rawdir)
  save_plot(plot_leadtime_prediction(lead,lead_windows),"c5.Fig4.leadtime_prediction.png",20,11.25,outdir=outdir)
  topn<-topn_score_benchmark(dat,ranked,training_screen,tvar,evar,split,clinical,C5_TOPN_MAX,C5_TOPN_BOOT)
  write_raw_csv(topn$rows,"c5.score_vs_topN_trajectory_rows.csv",rawdir)
  write_raw_csv(topn$summary,"c5.score_vs_topN.csv",rawdir);write_raw_csv(topn$bootstrap,"c5.score_vs_topN_bootstrap.csv",rawdir)
  write_raw_csv(topn$individuals,"c5.score_vs_topN_individuals.csv",rawdir);write_raw_csv(topn$individual_range,"c5.score_vs_topN_individual_range.csv",rawdir)
  save_plot(plot_topn_mechanism(topn),"c5.Fig5.score_vs_topN_mechanism.png",19,10.7,outdir=outdir)
  attained_score<-attained_age_score_sensitivity(score_rows,dat,lead_methods,clinical,tvar,evar)
  write_raw_csv(attained_score,"c5.attained_age_score_sensitivity.csv",rawdir)
  save_plot(plot_attained_age_score_sensitivity(attained_score),"c5.Fig12.attained_age_sensitivity.png",14,8.5,outdir=outdir)
  # Matching/trajectory bootstrap is diagnostic and is not consumed by the
  # five-row grid.  Limit it to the main paradigm scores so the 88k-person
  # metabolite branch cannot stall before any figure is written.
  pair_methods<-if(nrow(dat)>100000L)character()else intersect(fig1_order,method_names)
  if(!length(pair_methods)&&nrow(dat)>100000L)message("C5/",layer,": skip optional matched-pair bootstrap for N > 100,000; lead-time and top-N uncertainty are still estimated")
  pairs_native<-prior$preclinical_pairs%||%setNames(vector("list",length(pair_methods)),pair_methods)
  boot_by_method<-prior$yy_auc_bootstrap%||%setNames(vector("list",length(pair_methods)),pair_methods)
  for(i in seq_along(pair_methods)){
    nm<-pair_methods[[i]]
    if(is.null(pairs_native[[nm]])||!nrow(pairs_native[[nm]]))pairs_native[[nm]]<-tryCatch(risk_set_pairs(score_rows,dat,nm,"score_z",SEED+100+i),error=function(e){warning(nm,": pair construction failed: ",conditionMessage(e),call.=FALSE);tibble()})
    if(C5_BOOT>0&&(is.null(boot_by_method[[nm]])||!nrow(boot_by_method[[nm]])))boot_by_method[[nm]]<-tryCatch(bootstrap_auc_link(score_rows|>filter(method==nm),dat,C5_BOOT),error=function(e){warning(nm,": bootstrap diagnostic failed: ",conditionMessage(e),call.=FALSE);tibble()})
    if(is.null(boot_by_method[[nm]]))boot_by_method[[nm]]<-tibble()
  }
  write_raw_csv(bind_rows(imap(pairs_native,~.x|>mutate(method=.y))),"c5.preclinical_pairs_all.csv",rawdir)
  write_raw_csv(bind_rows(imap(boot_by_method,~.x|>mutate(method=.y))),"c5.yy_auc_bootstrap_all.csv",rawdir)

  f1<-c5_make_5row_grid(fig1_order,pred,score_rows,pairs_native,boot_by_method,"Prediction paradigms")
  f2<-c5_make_5row_grid(fig2_order,pred,score_rows,pairs_native,boot_by_method,"Evidence-screened prediction")
  f3<-c5_make_5row_grid(fig3_order,pred,score_rows,pairs_native,boot_by_method,"C4 connection-guided prediction")
  save_plot(f1,"c5.Fig1.simple_pred.png",24,13.5,outdir=outdir)
  save_plot(f2,"c5.Fig2.screened_pred.png",24,13.5,outdir=outdir)
  save_plot(f3,"c5.Fig3.connection_pred.png",24,13.5,outdir=outdir)

  # Optional nested-CV consolidation remains available; it is not required for Figs 1–3.
  cv<-list(pred=tibble(),metrics=tibble(),summary=tibble(),selected=tibble(),stability=tibble())
  if(C5_NESTED_CV){
    core<-unique(c(tvar,evar,vars.basic,vars.le8));dn<-dat[complete.cases(dat[,intersect(core,names(dat)),drop=FALSE]),,drop=FALSE]
    cv<-nested_cv_omics(dn,biom_vars,intersect(vars.le8,names(dn)),intersect(vars.basic,names(dn)),character(),tvar,evar,OUTER,INNER,HORIZON,SEED)
    write_raw_csv(cv$pred,"c5.prediction_rows_nested.csv",rawdir);write_raw_csv(cv$metrics,"c5.metrics_by_fold.csv",rawdir);write_raw_csv(cv$summary,"c5.metrics_summary.csv",rawdir);write_raw_csv(cv$selected,"c5.selected_by_fold.csv",rawdir);write_raw_csv(cv$stability,"c5.selection_stability.csv",rawdir)
  }
  ev<-make_evidence_table(layer,outdir,cv$stability);write_raw_csv(ev,"c5.evidence_consolidation.csv",rawdir)
  set_sizes<-tibble(set=c(paste0(all_glmnet_label," (all)"),paste0(lightgbm_label," (preselected)"),"User","PWAS/MWAS significant","MR significant",paste0(local_class," MR significant"),"Robust colocalization","Distal antecedent","Genetic-region evidence","Hybrid triangulated","Evidence-selected compact","Parsimonious","Extended","C4 NS","C4 YS","C4 YSplus","Genetic PRS"),
    N=c(length(biom_vars),length(yuvars),length(user),length(pwas),length(mrset),length(local_mrset),length(colocset),length(distalset),length(geneticset),length(hybridset),length(mechanismset),length(pars),length(extended),length(c4sets[[1]]),length(c4sets[[2]]),length(c4sets[[3]]),length(prs_vars)))
  candidate_sets<-bind_rows(tibble(set="PWAS/MWAS significant",feature=pwas),tibble(set="MR significant",feature=mrset),
    tibble(set=paste0(local_class," MR significant"),feature=local_mrset),tibble(set="Robust colocalization",feature=colocset),
    tibble(set="Distal antecedent",feature=distalset),tibble(set="Genetic-region evidence",feature=geneticset),
    tibble(set="Hybrid triangulated",feature=hybridset),tibble(set="Evidence-selected compact",feature=mechanismset),
    tibble(set="Parsimonious",feature=pars),
    tibble(set="Extended",feature=extended),bind_rows(imap(c4sets,~tibble(set=.y,feature=.x))))|>distinct()
  write_raw_csv(set_sizes,"c5.input_set_sizes.csv",rawdir)
  write_raw_csv(candidate_sets,"c5.candidate_sets.csv",rawdir)
  corr<-score_correlation_table(score_rows);write_raw_csv(corr,"c5.score_correlations.csv",rawdir)
  subgroup_methods<-intersect(c(all_glmnet_label,"Parsimonious","C4 YSplus"),method_names)
  sg<-subgroup_auc_table(score_rows,dat,subgroup_methods);write_raw_csv(sg,"c5.subgroup_auc.csv",rawdir)
  save_plot(plot_evidence_matrix(ev)|plot_causal_reactive_map(ev),"c5.Fig6.evidence_matrix.png",18,10.5,outdir=outdir)
  save_plot(plot_performance_benchmark(psum),"c5.Fig7.performance_benchmark.png",13,9,outdir=outdir)
  save_plot(plot_incremental_performance(psum),"c5.Fig8.incremental_performance.png",13,9,outdir=outdir)
  save_plot(plot_complexity_performance(psum),"c5.Fig9.parsimony_performance.png",10.5,8,outdir=outdir)
  save_plot(plot_score_correlation(corr),"c5.Fig10.score_concordance.png",11,9.5,outdir=outdir)
  save_plot(plot_subgroup_auc(sg),"c5.Fig11.subgroup_discrimination.png",12,10,outdir=outdir)
  write_xlsx2(list(upstream_availability=upstream_audit,prediction_summary=psum,input_set_sizes=set_sizes,candidate_sets=candidate_sets,mechanism_weight_prior=mechanism_prior,training_screen=training_screen,
    distal_landmark_screen=distal_screen,score_coverage_audit=score_audit,sequential_forward=seq0$log,
    score_method_summary=score_rows|>group_by(method,kind,split)|>summarise(rows=n(),finite_scores=sum(is.finite(score_z)),mean_score=mean(score_z,na.rm=TRUE),sd_score=sd(score_z,na.rm=TRUE),.groups="drop"),
    preclinical_pairs=bind_rows(imap(pairs_native,~.x|>mutate(method=.y))),yy_auc_bootstrap=bind_rows(imap(boot_by_method,~.x|>mutate(method=.y))),
    leadtime_discrimination=lead,leadtime_headline=lead_headline,leadtime_windows=lead_windows,
    score_vs_topN=topn$summary,topN_individuals=topn$individuals,topN_individual_range=topn$individual_range,
    attained_age_sensitivity=attained_score,evidence_consolidation=ev,score_correlations=corr,subgroup_AUC=sg,nested_cv_summary=cv$summary),"c5.prediction_panels.xlsx")
  out<-list(meta=module_meta(layer,extra=list(N=nrow(dat),events=sum(dat[[evar]]==1),training=sum(split=="training"),validation=sum(split=="validation"),code_version=C5_CODE_VERSION,upstream_signature=upstream_signature)),
    upstream_availability=upstream_audit,
    scores=sc,prediction=pred,summary=psum,set_sizes=set_sizes,sequential=seq0,preclinical_pairs=pairs_native,yy_auc_bootstrap=boot_by_method,
    score_coverage_audit=score_audit,candidate_sets=candidate_sets,mechanism_weight_prior=mechanism_prior,distal_screen=distal_screen,leadtime=lead,leadtime_headline=lead_headline,leadtime_windows=lead_windows,
    individual_lead_rows=individual_lead_rows,score_vs_topN=topn,attained_age_sensitivity=attained_score,
    evidence=ev,score_correlations=corr,subgroup_AUC=sg,nested_cv=cv)
  saveRDS(out,final_cache,compress="xz");finalize_outputs(LE8_JOB,outdir);out
}

if(prot_DO){invisible(run_c5_layer("protein"));gc(full=TRUE)}
if(met_DO){invisible(run_c5_layer("metabolite"));gc(full=TRUE)}
