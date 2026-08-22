# C5: Consolidation — unified five-row score phenotyping and prediction.
#
# Figure contract (all three figures use the SAME five row functions):
#   Row 1  score distribution + incidence rate on Biomarker score (s.d.); dual y axes
#   Row 2  validation ROC (Clinical / Biomarkers / Combined) + small temporal-AUC inset
#   Row 3  cumulative incidence by biomarker-score quintile
#   Row 4  Yin–Yang incident-case score trajectory (training vs validation; z score centered at training mean)
#   Row 5  Nature Aging-style preclinical trajectory (controls / case training / case validation;
#          native 0–1 model score, NOT centered at the population mean)
#
# Fig1 columns: Pradeep / glmnet, Yu / LightGBM, User specified
# Fig2 columns: PWAS/MWAS significant, MR significant, Parsimonious, Extended
# Fig3 columns: C4 NS, C4 YS, C4 YSplus
#
# Pradeep-style models use logistic LASSO, random 80/20 outer split, 10-fold tuning,
# lambda.1se for proteomic/combined models, and lambda.min for the clinical model.
# Yu-style prediction uses LightGBM with C1-ranked preselection (default top 257),
# matching the reported 257-protein input size when enough proteins are available.

suppressPackageStartupMessages({
  fdir<-Sys.getenv("LE8_FDIR",unset=file.path(Sys.getenv("DIRSCRIPT"),"f"))
  source(file.path(fdir,"comm.f.R"));source(file.path(fdir,"c5_pred.R"))
})
LE8_JOB <- "c5_consolidate"
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
C5_TEMPORAL_STEP <- as.numeric(Sys.getenv("C5_TEMPORAL_STEP",unset="0.5"))
C5_TIME_CAP <- as.numeric(Sys.getenv("C5_TIME_CAP",unset="16"))
C5_BOOT <- as.integer(Sys.getenv("C5_AUC_BOOT",unset="150"))
C5_MIN_BIN_N <- as.integer(Sys.getenv("C5_MIN_BIN_N",unset="15"))
C5_INC_PROT <- Sys.getenv("C5_INC_PROT",unset="")
C5_INC_MET <- Sys.getenv("C5_INC_MET",unset="")

# Publication typography for C5 only.
theme_5c <- function(base_size=12){
  theme_classic(base_size=base_size)+theme(
    plot.title=element_text(face="bold",size=base_size*1.12,hjust=0),
    plot.subtitle=element_text(face="bold",size=base_size*.92,color="grey30"),
    axis.title=element_text(face="bold",size=base_size*1.03),axis.text=element_text(face="bold",size=base_size,color="black"),
    legend.title=element_text(face="bold"),legend.text=element_text(face="bold"),strip.background=element_blank(),strip.text=element_text(face="bold"),
    panel.grid.major.y=element_line(color="grey91",linewidth=.25),panel.grid.minor=element_blank(),plot.margin=margin(7,10,7,10))
}

# ----------------------------------------------------------------------------
# Fixed outer split and score fitting
# ----------------------------------------------------------------------------
make_outer_split <- function(dat,test_frac=C5_TEST_FRAC,seed=SEED){
  set.seed(seed);val<-sample(seq_len(nrow(dat)),max(1L,round(nrow(dat)*test_frac)))
  ifelse(seq_len(nrow(dat))%in%val,"validation","training")
}

fit_glmnet_score <- function(dat,vars,tvar,evar,split,label,inner=INNER,alpha=1,
                             lambda_rule=c("lambda.1se","lambda.min")){
  lambda_rule<-match.arg(lambda_rule)
  if(!requireNamespace("glmnet",quietly=TRUE))stop("C5 requires glmnet.",call.=FALSE)
  vars<-intersect(unique(vars),names(dat));trid<-which(split=="training");vaid<-which(split=="validation")
  if(!length(vars)||length(trid)<500||length(vaid)<100)return(NULL)
  tr<-dat[trid,,drop=FALSE];va<-dat[vaid,,drop=FALSE];x<-impute_train_test(tr,va,vars)
  if(ncol(x$tr)<1)return(NULL);y<-as.integer(tr[[evar]]==1);if(length(unique(y))<2)return(NULL)
  k<-min(inner,max(3L,min(table(y))))
  foldid<-make_folds(tr,evar,k,SEED+17)
  fit<-tryCatch(glmnet::cv.glmnet(x$tr,y,family="binomial",alpha=alpha,foldid=foldid,
    type.measure="auc",standardize=FALSE,keep=TRUE),error=function(e)NULL)
  if(is.null(fit))return(NULL)
  lam<-if(lambda_rule=="lambda.min")fit$lambda.min else fit$lambda.1se
  li<-which.min(abs(fit$lambda-lam));oof<-rep(NA_real_,nrow(tr))
  if(!is.null(fit$fit.preval)&&length(li))oof<-as.numeric(fit$fit.preval[,li])
  val<-tryCatch(as.numeric(predict(fit,x$te,s=lambda_rule,type="link")),error=function(e)rep(NA_real_,nrow(va)))
  b<-as.matrix(coef(fit,s=lambda_rule));sel<-setdiff(rownames(b)[b[,1]!=0],"(Intercept)")
  rows<-bind_rows(
    tibble(eid=tr$eid,row_id=trid,split="training",time=tr[[tvar]],event=tr[[evar]],score=oof),
    tibble(eid=va$eid,row_id=vaid,split="validation",time=va[[tvar]],event=va[[evar]],score=val))|>
    mutate(method=label,n_features=length(x$vars),n_selected=length(sel),engine="glmnet",lambda_rule=lambda_rule)
  list(rows=rows,selected=sel,available=x$vars,fit=fit,engine="glmnet",lambda_rule=lambda_rule)
}

fit_lightgbm_score <- function(dat,vars,tvar,evar,split,label="Yu / LightGBM",k=C5_LGB_FOLDS,nrounds=C5_LGB_NROUND){
  if(!requireNamespace("lightgbm",quietly=TRUE)){
    warning("R package 'lightgbm' is not installed; Yu/LightGBM score will be skipped.",call.=FALSE);return(NULL)
  }
  vars<-intersect(unique(vars),names(dat));trid<-which(split=="training");vaid<-which(split=="validation")
  if(!length(vars)||length(trid)<500||length(vaid)<100)return(NULL)
  tr<-dat[trid,,drop=FALSE];va<-dat[vaid,,drop=FALSE]
  # Yu-style LightGBM keeps missing values natively. Only remove invariant columns.
  keep<-vapply(tr[,vars,drop=FALSE],function(z){z<-suppressWarnings(as.numeric(z));sum(is.finite(z))>=100&&is.finite(sd(z,na.rm=TRUE))&&sd(z,na.rm=TRUE)>0},logical(1))
  vars<-vars[keep];if(length(vars)<2)return(NULL)
  to_mat<-function(d){x<-as.matrix(data.frame(lapply(d[,vars,drop=FALSE],function(z)suppressWarnings(as.numeric(z))),check.names=FALSE));storage.mode(x)<-"double";x}
  xtr<-to_mat(tr);xva<-to_mat(va);y<-as.integer(tr[[evar]]==1);if(length(unique(y))<2)return(NULL)
  foldid<-make_folds(tr,evar,min(k,max(3L,min(table(y)))),SEED+23);oof<-rep(NA_real_,nrow(tr))
  params<-list(objective="binary",metric="auc",max_depth=15L,num_leaves=10L,
    bagging_fraction=.70,bagging_freq=1L,learning_rate=.01,feature_fraction=.70,max_bin=63L,
    lambda_l2=1,verbosity=-1L,num_threads=max(1L,N_CORES),seed=SEED)
  for(fd in sort(unique(foldid))){
    it<-which(foldid!=fd);iv<-which(foldid==fd);if(length(unique(y[it]))<2)next
    ds<-lightgbm::lgb.Dataset(data=xtr[it,,drop=FALSE],label=y[it])
    ff<-tryCatch(lightgbm::lgb.train(params=params,data=ds,nrounds=nrounds,verbose=-1),error=function(e)NULL)
    if(!is.null(ff))oof[iv]<-as.numeric(predict(ff,xtr[iv,,drop=FALSE]))
  }
  ds<-lightgbm::lgb.Dataset(data=xtr,label=y);fit<-tryCatch(lightgbm::lgb.train(params=params,data=ds,nrounds=nrounds,verbose=-1),error=function(e)NULL)
  if(is.null(fit))return(NULL);val<-as.numeric(predict(fit,xva))
  imp<-tryCatch(lightgbm::lgb.importance(fit),error=function(e)NULL);sel<-if(is.null(imp))vars else as.character(imp$Feature[imp$Gain>0])
  rows<-bind_rows(
    tibble(eid=tr$eid,row_id=trid,split="training",time=tr[[tvar]],event=tr[[evar]],score=oof),
    tibble(eid=va$eid,row_id=vaid,split="validation",time=va[[tvar]],event=va[[evar]],score=val))|>
    mutate(method=label,n_features=length(vars),n_selected=length(sel),engine="lightgbm",lambda_rule=NA_character_)
  list(rows=rows,selected=sel,available=vars,fit=fit,engine="lightgbm")
}

standardize_scores <- function(rows){
  rows|>group_by(method)|>group_modify(function(d,key){
    m<-mean(d$score[d$split=="training"],na.rm=TRUE);s<-sd(d$score[d$split=="training"],na.rm=TRUE);if(!is.finite(s)||s==0)s<-1
    d|>mutate(score_z=(score-m)/s,score_native=ifelse(engine=="lightgbm",pmin(pmax(score,0),1),plogis(score)))
  })|>ungroup()
}

fit_combined_score <- function(dat,score_obj,clinical,tvar,evar,split,label){
  if(is.null(score_obj))return(NULL)
  vars<-unique(c(score_obj$available,clinical))
  if(identical(score_obj$engine,"lightgbm"))fit_lightgbm_score(dat,vars,tvar,evar,split,paste0(label," + clinical"))
  else fit_glmnet_score(dat,vars,tvar,evar,split,paste0(label," + clinical"),lambda_rule="lambda.1se")
}

make_validation_bundle <- function(dat,score_obj,clinical_obj,clinical,tvar,evar,split,label){
  if(is.null(score_obj))return(tibble())
  bio<-score_obj$rows|>filter(split=="validation")|>transmute(eid,row_id,time,event,score,biom_set=label,model="Biomarkers",n_selected=first(score_obj$rows$n_selected))
  comb<-fit_combined_score(dat,score_obj,clinical,tvar,evar,split,label)
  cc<-if(is.null(comb))tibble()else comb$rows|>filter(split=="validation")|>transmute(eid,row_id,time,event,score,biom_set=label,model="Combined",n_selected=first(score_obj$rows$n_selected))
  cl<-if(is.null(clinical_obj))tibble()else clinical_obj$rows|>filter(split=="validation")|>transmute(eid,row_id,time,event,score,biom_set=label,model="Model 0",n_selected=first(score_obj$rows$n_selected))
  bind_rows(bio,cc,cl)|>group_by(biom_set,model)|>mutate(cindex=fcidx(Surv(time,event),score))|>ungroup()
}

# ----------------------------------------------------------------------------
# Common temporal helpers
# ----------------------------------------------------------------------------
yy_score_summary <- function(rows,method,bins=C5_YY_BINS,min_bin_n=C5_MIN_BIN_N){
  d<-rows|>filter(method==.env$method,event==1,is.finite(time),is.finite(score_z));if(!nrow(d))return(list(lines=tibble(),hist=tibble()))
  lim<-max(d$time,na.rm=TRUE);br<-seq(0,lim,length.out=bins+1);mids<-(head(br,-1)+tail(br,-1))/2
  lines<-d|>mutate(bin=cut(time,br,include.lowest=TRUE,labels=FALSE))|>filter(!is.na(bin))|>group_by(split,bin)|>
    summarise(mean=mean(score_z),sd=sd(score_z),se=sd/sqrt(n()),N=n(),N_total=n_distinct(eid),.groups="drop")|>filter(N>=min_bin_n)|>mutate(year=mids[bin])
  hh<-hist(d$time,breaks=br,plot=FALSE);histdf<-tibble(xmin=head(br,-1),xmax=tail(br,-1),mid=mids,n=hh$counts)
  list(lines=lines,hist=histdf)
}

risk_set_pairs <- function(rows,dat,method,score_col=c("score_z","score_native"),seed=SEED){
  score_col<-match.arg(score_col)
  keepcov<-intersect(c("eid","age","sex","ethnic.c"),names(dat))
  d<-rows|>filter(method==.env$method,is.finite(.data[[score_col]]),is.finite(time))|>left_join(dat|>select(all_of(keepcov)),by="eid")
  out<-list();ii<-0L;set.seed(seed)
  for(sp in intersect(c("training","validation"),unique(d$split))){
    ca<-d|>filter(split==sp,event==1)|>arrange(time);co<-d|>filter(split==sp,event==0)
    if(!nrow(ca)||!nrow(co))next
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
match_bootstrap_pairs <- function(v,seed){
  set.seed(seed);ca<-v|>filter(event==1);co<-v|>filter(event==0);if(!nrow(ca)||!nrow(co))return(tibble())
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
  if(!requireNamespace("pROC",quietly=TRUE)||nrow(v)<200)return(tibble())
  map_dfr(seq_len(B),function(b){
    set.seed(SEED+1000+b);ca<-v|>filter(event==1);co<-v|>filter(event==0);if(!nrow(ca)||!nrow(co))return(tibble())
    vb<-bind_rows(ca[sample(seq_len(nrow(ca)),nrow(ca),replace=TRUE),,drop=FALSE],co[sample(seq_len(nrow(co)),nrow(co),replace=TRUE),,drop=FALSE])
    roc_auc<-tryCatch(as.numeric(pROC::auc(pROC::roc(vb$event,vb$score_z,quiet=TRUE,direction="<"))),error=function(e)NA_real_)
    pairs<-match_bootstrap_pairs(vb,SEED+2000+b);ta<-temporal_implied_auc(pairs)
    tibble(bootstrap=b,roc_auc=roc_auc,temporal_d=ta[["dbar"]],trajectory_implied_auc=ta[["implied_auc"]])
  })
}

# ----------------------------------------------------------------------------
# The five reusable row functions.  All C5 score figures call these directly.
# ----------------------------------------------------------------------------
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

c5_row2_roc <- function(pred,label,boot=tibble()){
  d<-pred|>filter(biom_set==.env$label);if(!nrow(d))return(blank_plot("Validation ROC","No validation predictions"))
  roc<-d|>group_by(model)|>group_modify(~{r<-tryCatch(pROC::roc(.x$event,.x$score,quiet=TRUE,direction="<"),error=function(e)NULL);if(is.null(r))return(tibble());tibble(fpr=1-r$specificities,tpr=r$sensitivities,auc=as.numeric(r$auc))})|>ungroup()
  if(!nrow(roc))return(blank_plot("Validation ROC","ROC unavailable"))
  al<-roc|>group_by(model)|>summarise(auc=first(auc),.groups="drop")|>mutate(txt=sprintf("%s: %.3f",case_when(model=="Model 0"~"Clinical",model=="Biomarkers"~"Biomarkers",TRUE~"Combined"),auc),x=.03,y=c(.97,.90,.83)[match(model,c("Model 0","Biomarkers","Combined"))])
  p<-ggplot(roc,aes(fpr,tpr,color=model))+geom_abline(slope=1,intercept=0,linetype=3,color="grey78")+geom_step(linewidth=.9)+
    geom_text(data=al,aes(x=x,y=y,label=txt,color=model),hjust=0,vjust=1,inherit.aes=FALSE,show.legend=FALSE,size=2.8,fontface="bold")+
    scale_color_manual(values=c(`Model 0`="grey35",Biomarkers="#C77732",Combined="#287A58"),breaks=c("Biomarkers","Combined","Model 0"),labels=c("Biomarkers","Combined","Clinical"))+
    labs(x="False positive rate",y="True positive rate",color=NULL)+theme_5c(9)+theme(legend.position="top")
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
  d<-score_rows|>filter(method==.env$label,split=="validation",is.finite(score_z),is.finite(time),time>0,!is.na(event))|>mutate(quintile=factor(ntile(score_z,5),levels=1:5,labels=paste0("Q",1:5)))
  if(!nrow(d))return(blank_plot("Cumulative incidence","No validation scores"))
  cum<-d|>group_by(quintile)|>group_modify(~{sf<-survfit(Surv(time,event)~1,data=.x);tibble(time=sf$time,cuminc=1-sf$surv)})|>ungroup()
  ymax<-max(cum$cuminc[cum$time<=C5_TIME_CAP],na.rm=TRUE);if(!is.finite(ymax)||ymax<=0)ymax<-max(cum$cuminc,na.rm=TRUE);ymax<-min(1,max(.01,ymax*1.10))
  ggplot(cum,aes(time,cuminc,color=quintile))+geom_step(linewidth=.75)+coord_cartesian(xlim=c(-C5_TIME_CAP,C5_TIME_CAP),ylim=c(0,ymax))+
    scale_color_brewer(palette="YlOrBr",direction=1)+labs(x="Time from baseline (years)",y="Cumulative incidence",color="Score quintile")+
    theme_5c(9)+theme(legend.position="top")
}

c5_row4_yy <- function(score_rows,label){
  d<-score_rows|>filter(method==.env$label,is.finite(time),time>=0,time<=C5_TIME_CAP,is.finite(score_z),!is.na(event))|>
    mutate(group=case_when(event==0~"Controls",split=="training"~"Cases: training",TRUE~"Cases: validation"),bin=floor(time/C5_TEMPORAL_STEP)*C5_TEMPORAL_STEP+C5_TEMPORAL_STEP/2)
  ln<-d|>group_by(group,bin)|>summarise(mean=mean(score_z),N=n(),N_total=n_distinct(eid),.groups="drop")|>filter(N>=C5_MIN_BIN_N)
  h<-d|>count(bin,name="n")|>transmute(xmin=bin-C5_TEMPORAL_STEP/2,xmax=bin+C5_TEMPORAL_STEP/2,n)
  if(!nrow(ln))return(blank_plot("Temporal score trajectory","No sufficiently populated bins"))
  ln<-ln|>group_by(group)|>mutate(group_label=sprintf("%s (N=%s)",first(group),format(first(N_total),big.mark=",")))|>ungroup()
  yr<-range(ln$mean,na.rm=TRUE);if(diff(yr)<.2)yr<-yr+c(-.5,.5);pad<-.13*diff(yr);yr<-yr+c(-pad,pad)
  count_max<-max(h$n,1,na.rm=TRUE);plot_max<-count_max*1.15;scale_factor<-plot_max/diff(yr)
  ln<-ln|>mutate(y_plot=(mean-yr[1])*scale_factor)
  lev<-unique(ln|>arrange(factor(group,levels=c("Controls","Cases: training","Cases: validation")))|>pull(group_label));base<-sub(" \\(N=.*$","",lev);pal0<-c(Controls="grey45",`Cases: training`="#2C7FB8",`Cases: validation`="#D95F02");pal<-setNames(unname(pal0[base]),lev)
  ggplot()+geom_rect(data=h,aes(xmin=xmin,xmax=xmax,ymin=0,ymax=n),fill="#7FB8DF",alpha=.31,inherit.aes=FALSE)+
    geom_hline(yintercept=(0-yr[1])*scale_factor,linetype=3,color="grey62")+geom_smooth(data=ln,aes(bin,y_plot,color=group_label,group=group_label),method="loess",se=FALSE,span=.65,linewidth=1.0)+geom_point(data=ln,aes(bin,y_plot,color=group_label),size=1.45,alpha=.7)+
    scale_y_continuous(name="Patient count",limits=c(0,plot_max),expand=expansion(mult=c(0,.02)),sec.axis=sec_axis(~.x/scale_factor+yr[1],name="Biomarker score (s.d.)"))+scale_color_manual(values=pal)+
    coord_cartesian(xlim=c(-C5_TIME_CAP,C5_TIME_CAP))+labs(x="Time from baseline (years)",color=NULL)+theme_5c(9)+
    guides(color=guide_legend(nrow=1,byrow=TRUE))+theme(legend.position="top")
}

c5_row5_preclinical <- function(pairs_native){
  d<-pair_long(pairs_native);if(!nrow(d))return(blank_plot("Preclinical score trajectory","No matched case-control trajectory"))
  sm<-d|>mutate(bin=floor(year/C5_TEMPORAL_STEP)*C5_TEMPORAL_STEP+C5_TEMPORAL_STEP/2)|>group_by(group,bin)|>
    summarise(mean=mean(score),se=sd(score)/sqrt(n()),N=n(),.groups="drop")|>filter(N>=C5_MIN_BIN_N)
  if(!nrow(sm))return(blank_plot("Preclinical score trajectory","No sufficiently populated matched bins"))
  nlab<-d|>group_by(group)|>summarise(N=n_distinct(eid),.groups="drop")|>mutate(group_label=paste0(group," (N=",format(N,big.mark=","),")"))
  sm<-sm|>left_join(nlab,by="group")
  # Native predictions are bounded by 0 and 1, but forcing every panel to that
  # full range compresses low-risk trajectories against the x axis. Use the
  # observed mean +/- SE range with modest padding, still clipped to [0, 1].
  ylo<-min(sm$mean-sm$se,na.rm=TRUE);yhi<-max(sm$mean+sm$se,na.rm=TRUE)
  if(!all(is.finite(c(ylo,yhi))))c(ylo,yhi)<-c(0,1)
  if(yhi-ylo<.02){mid<-(ylo+yhi)/2;ylo<-mid-.01;yhi<-mid+.01}
  ypad<-.12*(yhi-ylo);ylim_native<-c(max(0,ylo-ypad),min(1,yhi+ypad))
  gl<-unique(nlab$group_label);base_group<-sub(" \\(N=.*$","",gl);pal0<-c(Controls="grey45",`Cases: training`="#2C7FB8",`Cases: validation`="#D95F02");pal<-setNames(unname(pal0[base_group]),gl)
  ggplot(sm,aes(bin,mean,color=group_label,group=group_label))+
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=.14,linewidth=.35,alpha=.55)+geom_point(size=1.8)+
    geom_smooth(method="loess",se=FALSE,span=.75,linewidth=1.0)+geom_vline(xintercept=0,linetype=2,color="grey50")+
    scale_color_manual(values=pal)+coord_cartesian(xlim=c(-C5_TIME_CAP,C5_TIME_CAP),ylim=ylim_native)+
    labs(x="Time from baseline (years)",y="Model score (native 0–1 scale)",color=NULL)+theme_5c(9)+
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
    c5_row5_preclinical(pairs_native[[lb]]%||%tibble())),nrow=1,widths=c(1.15,1.05,1,1,1.08)))
  spaced<-list();for(i in seq_along(method_rows)){spaced[[length(spaced)+1]]<-method_rows[[i]];if(i<length(method_rows))spaced[[length(spaced)+1]]<-plot_spacer()}
  p<-wrap_plots(spaced,ncol=1,heights=rep(c(1,.08),length.out=length(spaced)))
  if(!is.null(title))p<-p+plot_annotation(title=title,theme=theme(plot.title=element_text(face="bold",size=16,hjust=.5)))
  p
}

# ----------------------------------------------------------------------------
# Evidence table: C2 includes MR and DANDELION, C3 coloc, C4 connection.
# ----------------------------------------------------------------------------
make_evidence_table <- function(layer,outdir,stability=tibble()){
  rr<-function(job){f<-file.path(le8_job_dir(outdir,job),paste0(sub("_.*$","",job),".res.rds"));if(file.exists(f))readRDS(f)else list()}
  c1<-rr("c1_correlate");c2<-rr("c2_cause");c3<-rr("c3_coloc");c4<-rr("c4_connect")
  a0<-c1$association%||%tibble();a<-if(nrow(a0))a0|>transmute(feature=term,obs_beta=beta,obs_p=p.value,obs_FDR=FDR)else tibble(feature=character())
  m0<-c2$MR%||%tibble();m<-if(nrow(m0))m0|>filter(is.finite(pval))|>group_by(exposure)|>slice_min(pval,n=1,with_ties=FALSE)|>ungroup()|>transmute(feature=exposure,MR_analysis=analysis,MR_beta=b,MR_p=pval,MR_FDR=FDR_all)else tibble(feature=character())
  d0<-c2$DANDELION$targets%||%tibble();dd<-if(nrow(d0))d0|>transmute(feature=gene2,DANDELION_p,DANDELION_loci=n_distal_loci,DANDELION=TRUE)else tibble(feature=character(),DANDELION_p=numeric(),DANDELION_loci=integer(),DANDELION=logical())
  co0<-c3$summary%||%tibble();co<-if(nrow(co0)&&all(c("feature","PP.H4")%in%names(co0)))co0|>filter(status=="ok")|>group_by(feature)|>slice_max(PP.H4,n=1,with_ties=FALSE)|>ungroup()|>select(feature,PP.H4,everything())else tibble(feature=character())
  cx0<-c4$membership%||%tibble();cx<-if(nrow(cx0)&&"feature"%in%names(cx0))cx0 else tibble(feature=character())
  st<-if(nrow(stability)&&all(c("feature","selection_frequency")%in%names(stability)))stability|>group_by(feature)|>summarise(max_selection_frequency=max(selection_frequency,na.rm=TRUE),.groups="drop")else tibble(feature=character(),max_selection_frequency=numeric())
  z<-Reduce(function(x,y)full_join(x,y,by="feature"),list(a,m,dd,co,cx,st));if(!nrow(z))return(tibble())
  defaults<-list(obs_FDR=NA_real_,obs_p=NA_real_,MR_FDR=NA_real_,DANDELION=FALSE,PP.H4=NA_real_,strict_YS=FALSE,max_selection_frequency=NA_real_)
  for(nm in names(defaults))if(!nm%in%names(z))z[[nm]]<-defaults[[nm]]
  z|>mutate(C1=as.integer(is.finite(obs_FDR)&obs_FDR<.05),C2_MR=as.integer(is.finite(MR_FDR)&MR_FDR<.05),C2_DANDELION=as.integer(coalesce(DANDELION,FALSE)),
    C2=as.integer(C2_MR==1|C2_DANDELION==1),C3=as.integer(is.finite(PP.H4)&PP.H4>=.7),C4=as.integer(coalesce(strict_YS,FALSE)),
    C5=as.integer(is.finite(max_selection_frequency)&max_selection_frequency>=.6),evidence_count=C1+C2+C3+C4+C5)|>arrange(desc(evidence_count),obs_p)
}

# ----------------------------------------------------------------------------
# Main C5
# ----------------------------------------------------------------------------
run_c5_layer <- function(layer=c("protein","metabolite")){
  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir);rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
  # Delete obsolete standalone figures; their information is integrated into Figs 1–3.
  unlink(file.path(outdir,c("c5.Fig1.prediction.png","c5.Fig2.incremental.png","c5.Fig3.yy_score.png","c5.Fig4.case_control_yy.png","c5.Fig5.yy_vs_roc_auc.png","c5.Fig3.connection_pred.png")),force=TRUE)
  biom<-if(layer=="protein")read_prot()else read_met();biom_vars<-setdiff(names(biom),"eid")
  incfile<-if(layer=="protein")C5_INC_PROT else C5_INC_MET
  if(nzchar(incfile)&&file.exists(incfile)){req<-unique(scan(incfile,what="character",quiet=TRUE));biom_vars<-intersect(biom_vars,req);biom<-biom[,c("eid",biom_vars),drop=FALSE]}
  all0<-read_all();need<-unique(c("eid","ethnic.c",vars.basic,vars.le8,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
  dat<-all0[,intersect(need,names(all0)),drop=FALSE]|>filter_analysis_cohort()|>make_outcome(Y)|>inner_join(biom,by="eid");rm(all0,biom);gc()
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");clinical<-intersect(covs_use,names(dat));biom_vars<-intersect(biom_vars,names(dat))
  dat<-dat[complete.cases(dat[,c(tvar,evar),drop=FALSE])&dat[[tvar]]>0,,drop=FALSE]
  split<-make_outer_split(dat);dat$c5_split<-split
  message("C5/",layer,": N=",nrow(dat),", events=",sum(dat[[evar]]==1),", training=",sum(split=="training"),", validation=",sum(split=="validation"),", biomarkers=",length(biom_vars))

  c1f<-file.path(le8_job_dir(outdir,"c1_correlate"),"c1.res.rds");c1<-if(file.exists(c1f))readRDS(c1f)else list();a1<-c1$association%||%tibble()
  ranked<-if(nrow(a1))a1|>filter(is.finite(p.value))|>arrange(p.value)|>pull(term)|>intersect(biom_vars)else biom_vars
  pwas<-if(nrow(a1))a1|>filter(is.finite(FDR),FDR<.05)|>arrange(p.value)|>pull(term)|>intersect(biom_vars)else character()
  c2f<-file.path(le8_job_dir(outdir,"c2_cause"),"c2.res.rds");c2<-if(file.exists(c2f))readRDS(c2f)else list();m2<-c2$MR%||%tibble()
  mrset<-if(nrow(m2))m2|>filter(is.finite(FDR_all),FDR_all<.05)|>arrange(pval)|>pull(exposure)|>unique()|>intersect(biom_vars)else character()
  # For this analysis, "User specified" means the ten strongest C1 biomarkers.
  user<-head(ranked,10L)
  yuvars<-if(nzchar(C5_YU_LIST)&&file.exists(C5_YU_LIST)) intersect(unique(scan(C5_YU_LIST,what="character",quiet=TRUE)),biom_vars) else head(ranked,C5_YU_PRESELECT_N)
  if(!nzchar(C5_YU_LIST)) message("C5: C5_YU_LIST not supplied; Yu column is a Yu-style 257-protein proxy using the C1 ranking, not an exact reproduction of the published 257-protein panel.")

  # Nature Aging-style sequential screen, learned only inside the outer training sample.
  seq0<-sequential_forward_panel(dat[split=="training",,drop=FALSE],pwas,tvar,evar,HORIZON,SEED+71,max_extended=C5_EXTENDED_N)
  pars<-intersect(seq0$parsimonious,biom_vars);extended<-intersect(seq0$extended,biom_vars)
  write_raw_csv(seq0$log,"c5.sequential_forward_training.csv",rawdir)

  # C4 connection-defined candidate sets for Fig3.
  c4f<-file.path(le8_job_dir(outdir,"c4_connect"),"c4.res.rds");c4<-if(file.exists(c4f))readRDS(c4f)else list();l4<-c4$lists%||%list()
  required_c4_lists<-c("NS","YS_all","YSP_all")
  if(!all(required_c4_lists%in%names(l4)))stop("C5 requires the final C4 list contract: ",paste(required_c4_lists,collapse=", "),call.=FALSE)
  c4sets<-list(`C4 NS`=intersect(l4$NS,biom_vars),
               `C4 YS`=intersect(l4$YS_all,biom_vars),
               `C4 YSplus`=intersect(l4$YSP_all,biom_vars))
  empty_c4<-names(c4sets)[lengths(c4sets)==0L]
  if(length(empty_c4))stop("C5 connection prediction requires non-empty final C4 sets: ",paste(empty_c4,collapse=", "),call.=FALSE)

  cache<-file.path(rawdir,"c5.score_models.rds");sc<-if(cache_valid(cache))tryCatch(readRDS(cache),error=function(e)NULL)else NULL
  input_signature<-list(user=user,C4=c4sets)
  if(!is.null(sc)&&!identical(sc$input_signature,input_signature)){
    message("C5: model inputs changed; refitting model collection.");sc<-NULL
  }
  if(is.null(sc)){
    message("C5: fit Pradeep / glmnet all-biomarker score")
    pr<-fit_glmnet_score(dat,biom_vars,tvar,evar,split,"Pradeep / glmnet",lambda_rule="lambda.1se")
    message("C5: fit Yu / LightGBM on ",length(yuvars)," C1-ranked biomarkers")
    yu<-fit_lightgbm_score(dat,yuvars,tvar,evar,split,"Yu / LightGBM")
    us<-fit_glmnet_score(dat,user,tvar,evar,split,"User specified",lambda_rule="lambda.1se")
    pw<-fit_glmnet_score(dat,pwas,tvar,evar,split,"PWAS/MWAS significant",lambda_rule="lambda.1se")
    mr<-fit_glmnet_score(dat,mrset,tvar,evar,split,"MR significant",lambda_rule="lambda.1se")
    pa<-fit_glmnet_score(dat,pars,tvar,evar,split,"Parsimonious",lambda_rule="lambda.1se")
    ex<-fit_glmnet_score(dat,extended,tvar,evar,split,"Extended",lambda_rule="lambda.1se")
    ns<-fit_glmnet_score(dat,c4sets[["C4 NS"]],tvar,evar,split,"C4 NS",lambda_rule="lambda.1se")
    ys<-fit_glmnet_score(dat,c4sets[["C4 YS"]],tvar,evar,split,"C4 YS",lambda_rule="lambda.1se")
    yp<-fit_glmnet_score(dat,c4sets[["C4 YSplus"]],tvar,evar,split,"C4 YSplus",lambda_rule="lambda.1se")
    clinical_obj<-fit_glmnet_score(dat,clinical,tvar,evar,split,"Clinical",lambda_rule="lambda.min")
    objs<-list(Pradeep=pr,Yu=yu,User=us,PWAS=pw,MR=mr,Parsimonious=pa,Extended=ex,C4_NS=ns,C4_YS=ys,C4_YSplus=yp)
    score_rows<-standardize_scores(bind_rows(lapply(Filter(Negate(is.null),objs),`[[`,"rows")))
    sc<-c(objs,list(Clinical=clinical_obj,rows=score_rows,split=split,input_signature=input_signature,
      sets=list(PWAS=pwas,MR=mrset,Parsimonious=pars,Extended=extended,Yu=yuvars,User=user,C4=c4sets)))
    saveRDS(sc,cache,compress="xz")
  } else {message("C5: reuse score-model cache");score_rows<-sc$rows;clinical_obj<-sc$Clinical}
  write_raw_csv(score_rows,"c5.person_scores.csv",rawdir)

  object_map<-list(`Pradeep / glmnet`=sc$Pradeep,`Yu / LightGBM`=sc$Yu,`User specified`=sc$User,
    `PWAS/MWAS significant`=sc$PWAS,`MR significant`=sc$MR,Parsimonious=sc$Parsimonious,Extended=sc$Extended,
    `C4 NS`=sc$C4_NS,`C4 YS`=sc$C4_YS,`C4 YSplus`=sc$C4_YSplus)
  prior_file<-file.path(rawdir,"c5.res.rds");prior<-if(file.exists(prior_file))tryCatch(readRDS(prior_file),error=function(e)NULL)else NULL
  pred<-prior$prediction%||%tibble()
  if(!nrow(pred)){
    bundles<-imap(object_map,function(obj,nm)if(is.null(obj))tibble()else make_validation_bundle(dat,obj,clinical_obj,clinical,tvar,evar,split,nm))
    pred<-bind_rows(bundles)
  }
  write_raw_csv(pred,"c5.prediction_validation_rows.csv",rawdir)
  psum<-pred|>group_by(biom_set,model)|>summarise(N=n(),events=sum(event),C_index=first(cindex),
    AUC=tryCatch(as.numeric(pROC::auc(pROC::roc(event,score,quiet=TRUE,direction="<"))),error=function(e)NA_real_),n_selected=median(n_selected,na.rm=TRUE),.groups="drop")
  write_raw_csv(psum,"c5.prediction_summary.csv",rawdir)

  # Build reusable temporal objects for every score once; all three figures consume them.
  method_names<-unique(score_rows$method)
  pairs_native<-prior$preclinical_pairs%||%setNames(vector("list",length(method_names)),method_names)
  boot_by_method<-prior$yy_auc_bootstrap%||%setNames(vector("list",length(method_names)),method_names)
  for(i in seq_along(method_names)){
    nm<-method_names[[i]]
    if(is.null(pairs_native[[nm]])||!nrow(pairs_native[[nm]]))pairs_native[[nm]]<-risk_set_pairs(score_rows,dat,nm,"score_native",SEED+100+i)
    if(C5_BOOT>0&&(is.null(boot_by_method[[nm]])||!nrow(boot_by_method[[nm]])))boot_by_method[[nm]]<-bootstrap_auc_link(score_rows|>filter(method==nm),dat,C5_BOOT)
    if(is.null(boot_by_method[[nm]]))boot_by_method[[nm]]<-tibble()
  }
  write_raw_csv(bind_rows(imap(pairs_native,~.x|>mutate(method=.y))),"c5.preclinical_pairs_all.csv",rawdir)
  write_raw_csv(bind_rows(imap(boot_by_method,~.x|>mutate(method=.y))),"c5.yy_auc_bootstrap_all.csv",rawdir)

  fig1_order<-c("Pradeep / glmnet","Yu / LightGBM","User specified")
  fig2_order<-c("PWAS/MWAS significant","MR significant","Parsimonious")
  fig3_order<-c("C4 NS","C4 YS","C4 YSplus")
  f1<-c5_make_5row_grid(fig1_order,pred,score_rows,pairs_native,boot_by_method,"Prediction paradigms")
  f2<-c5_make_5row_grid(fig2_order,pred,score_rows,pairs_native,boot_by_method,"Evidence-screened prediction")
  f3<-c5_make_5row_grid(fig3_order,pred,score_rows,pairs_native,boot_by_method,"C4 connection-guided prediction")
  save_plot(f1,"c5.Fig1.simple_pred.png",26,18,outdir=outdir)
  save_plot(f2,"c5.Fig2.screened_pred.png",26,18,outdir=outdir)
  save_plot(f3,"c5.Fig3.connection_pred.png",26,18,outdir=outdir)

  # Optional nested-CV consolidation remains available; it is not required for Figs 1–3.
  cv<-list(pred=tibble(),metrics=tibble(),summary=tibble(),selected=tibble(),stability=tibble())
  if(C5_NESTED_CV){
    if(!file.exists(c4f))stop("Run C4 first for --nested-cv: ",c4f,call.=FALSE)
    core<-unique(c(tvar,evar,vars.basic,vars.le8));dn<-dat[complete.cases(dat[,intersect(core,names(dat)),drop=FALSE]),,drop=FALSE]
    cv<-nested_cv_omics(dn,biom_vars,intersect(vars.le8,names(dn)),intersect(vars.basic,names(dn)),character(),tvar,evar,OUTER,INNER,HORIZON,SEED)
    write_raw_csv(cv$pred,"c5.prediction_rows_nested.csv",rawdir);write_raw_csv(cv$metrics,"c5.metrics_by_fold.csv",rawdir);write_raw_csv(cv$summary,"c5.metrics_summary.csv",rawdir);write_raw_csv(cv$selected,"c5.selected_by_fold.csv",rawdir);write_raw_csv(cv$stability,"c5.selection_stability.csv",rawdir)
  }
  ev<-make_evidence_table(layer,outdir,cv$stability);write_raw_csv(ev,"c5.evidence_consolidation.csv",rawdir)
  set_sizes<-tibble(set=c("Pradeep_all","Yu_preselected","User","PWAS/MWAS significant","MR significant","Parsimonious","Extended","C4 NS","C4 YS","C4 YSplus"),
    N=c(length(biom_vars),length(yuvars),length(user),length(pwas),length(mrset),length(pars),length(extended),length(c4sets[[1]]),length(c4sets[[2]]),length(c4sets[[3]])))
  write_raw_csv(set_sizes,"c5.input_set_sizes.csv",rawdir)
  write_xlsx2(list(prediction_summary=psum,input_set_sizes=set_sizes,sequential_forward=seq0$log,person_scores=score_rows,
    preclinical_pairs=bind_rows(imap(pairs_native,~.x|>mutate(method=.y))),yy_auc_bootstrap=bind_rows(imap(boot_by_method,~.x|>mutate(method=.y))),
    evidence_consolidation=ev,nested_cv_summary=cv$summary),"c5.prediction_panels.xlsx")
  out<-list(meta=module_meta(layer,extra=list(N=nrow(dat),events=sum(dat[[evar]]==1),training=sum(split=="training"),validation=sum(split=="validation"))),
    scores=sc,prediction=pred,summary=psum,set_sizes=set_sizes,sequential=seq0,preclinical_pairs=pairs_native,yy_auc_bootstrap=boot_by_method,evidence=ev,nested_cv=cv)
  saveRDS(out,file.path(rawdir,"c5.res.rds"),compress="xz");finalize_outputs(LE8_JOB,outdir);out
}

if(prot_DO){invisible(run_c5_layer("protein"));gc(full=TRUE)}
if(met_DO){invisible(run_c5_layer("metabolite"));gc(full=TRUE)}
