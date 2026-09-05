# Valid interval-specific Cox risk sets; one baseline measurement per person.
le8_interval_cox <- function(dat,x,tvar,evar,covars,lo,hi,scale_x=TRUE,min_event=20L) {
  covars<-intersect(covars,names(dat));cols<-unique(c(x,tvar,evar,covars))
  d<-as.data.frame(dat[,cols,drop=FALSE]);d<-d[complete.cases(d),,drop=FALSE]
  d<-d[is.finite(d[[tvar]])&d[[tvar]]>lo&d[[evar]]%in%c(0,1),,drop=FALSE]
  d$.stop<-pmin(d[[tvar]],hi)-lo;d$.event<-as.integer(d[[evar]]==1&d[[tvar]]<=hi)
  ans<-tibble(term=x,beta=NA_real_,std.error=NA_real_,conf.low=NA_real_,conf.high=NA_real_,
    p.value=NA_real_,N_total=nrow(d),N_event=sum(d$.event),N_case=sum(d$.event),
    N_control=sum(d$.event==0),person_years=sum(d$.stop),window_lo=lo,window_hi=hi,
    time=(lo+hi)/2,effect_measure="log HR; interval-specific Cox",
    status="insufficient events or variation")
  if(!nrow(d)||sum(d$.event)<min_event||!is.finite(sd(d[[x]]))||sd(d[[x]])<=0)return(ans)
  # Fixed baseline reference SD, not a different SD in every late risk set.
  if(scale_x){ref<-as.numeric(dat[[x]]);s<-sd(ref,na.rm=TRUE);m<-mean(ref,na.rm=TRUE);d[[x]]<-(d[[x]]-m)/s}
  ff<-as.formula(paste0("survival::Surv(.stop,.event) ~ ",paste(bt(c(x,covars)),collapse=" + ")))
  fit<-tryCatch(survival::coxph(ff,d,ties="efron"),error=function(e)NULL)
  if(is.null(fit))return(ans)
  sm<-coef(summary(fit));if(!x%in%rownames(sm))return(ans)
  b<-sm[x,"coef"];se<-sm[x,"se(coef)"]
  ans$beta<-b;ans$std.error<-se;ans$conf.low<-b-1.96*se;ans$conf.high<-b+1.96*se
  ans$p.value<-sm[x,"Pr(>|z|)"];ans$status<-"ok";ans
}
risk_window_scan <- function(dat,xs,tvar,evar,bvar,covars,cuts=C1_RISK_CUTS) {
  xs<-intersect(xs,names(dat));cuts<-sort(unique(cuts))
  inc<-map_dfr(seq_len(length(cuts)-1L),function(i)map_dfr(xs,function(x)
    le8_interval_cox(dat,x,tvar,evar,covars,cuts[i],cuts[i+1])) )|>
    mutate(side="Post-baseline incident")
  # Prevalence remains logistic, explicitly separate from prospective HR.
  prev<-map_dfr(seq_len(length(cuts)-1L),function(i) {
    lo<-cuts[i];hi<-cuts[i+1]
    pc<-dat$prevalent_status==1&is.finite(dat[[bvar]])&(-dat[[bvar]])>=lo&(-dat[[bvar]])<hi
    cc<-dat$prevalent_status==0;keep<-replace_na(pc,FALSE)|replace_na(cc,FALSE)
    d<-dat[keep,,drop=FALSE];d$.window_case<-as.integer(replace_na(pc[keep],FALSE))
    z<-logistic_scan(d,xs,covars,".window_case")
    z|>transmute(term,beta,std.error,conf.low=beta-1.96*std.error,conf.high=beta+1.96*std.error,
      p.value,N_total,N_event,N_case=N_event,N_control=N_total-N_event,person_years=NA_real_,
      window_lo=lo,window_hi=hi,time=-(lo+hi)/2,effect_measure="log OR; baseline prevalence",
      status=ifelse(is.finite(p.value),"ok","insufficient data"),side="Pre-baseline prevalent")
  })
  bind_rows(prev,inc)|>group_by(side)|>mutate(FDR=p.adjust(p.value,"BH"))|>ungroup()
}
le8_time_heterogeneity <- function(dat,xs,tvar,evar,covars,cuts=c(0,2,5,10,16)) {
  if(!length(xs))return(tibble(feature=character(),N=integer(),events=integer(),p_time_interaction=numeric(),status=character(),FDR=numeric()))
  map_dfr(xs,function(x){
    cols<-unique(c(x,tvar,evar,covars));d<-as.data.frame(dat[,cols,drop=FALSE]);d<-d[complete.cases(d),,drop=FALSE]
    d<-d[d[[tvar]]>0&d[[evar]]%in%c(0,1),,drop=FALSE]
    ans<-tibble(feature=x,N=nrow(d),events=sum(d[[evar]]==1),p_time_interaction=NA_real_,status="insufficient data")
    if(nrow(d)<500||sum(d[[evar]]==1)<40||sd(d[[x]])<=0)return(ans)
    d$.z<-as.numeric(scale(d[[x]]))
    pieces<-lapply(seq_len(length(cuts)-1L),function(j){
      z<-d[d[[tvar]]>cuts[j],,drop=FALSE];z$.start<-cuts[j];z$.end<-pmin(z[[tvar]],cuts[j+1])
      z$.ev<-as.integer(z[[evar]]==1&z[[tvar]]<=cuts[j+1]);z$.window<-j;z
    })
    long<-bind_rows(pieces);long$.window<-factor(long$.window)
    cv<-paste(bt(covars),collapse=" + ");cv<-if(nzchar(cv))paste0(" + ",cv)else""
    f0<-as.formula(paste0("survival::Surv(.start,.end,.ev) ~ strata(.window) + .z",cv))
    f1<-as.formula(paste0("survival::Surv(.start,.end,.ev) ~ strata(.window) + .z:factor(.window)",cv))
    a<-tryCatch(survival::coxph(f0,long),error=function(e)NULL);b<-tryCatch(survival::coxph(f1,long),error=function(e)NULL)
    if(is.null(a)||is.null(b))return(ans)
    df<-sum(is.finite(coef(b)))-sum(is.finite(coef(a)))
    if(df>0){ans$p_time_interaction<-pchisq(max(0,2*(b$loglik[2]-a$loglik[2])),df,lower.tail=FALSE);ans$status<-"interval-slope likelihood-ratio test"}
    ans
  })|>mutate(FDR=p.adjust(p_time_interaction,"BH"))
}
le8_c1_additions <- function(dat,layer,covars,tvar,evar) {
  anchors<-intersect(le8_csv_env("C1_DIRECTION_ANCHORS","PCSK9,LPA,GDF15,NTPROBNP,MMP12,L_VLDL_TG.pct,L_VLDL_TG,Total_TG,ApoB"),names(dat))
  heterogeneity<-le8_time_heterogeneity(dat,anchors,tvar,evar,covars)
  rawdir<-le8_job_dir(if(layer=="protein")out.prot else out.met,"c1_correlate")
  paired<-tibble();status<-tibble(status="PGS input unavailable")
  sf<-find_c1_pgs_file(layer)
  if(length(sf)==1L&&!is.na(sf)&&length(anchors)) {
    le8_dependency(sf);scores<-read_c1_pgs(sf)
    mp<-map_c1_pgs_columns(anchors,names(scores))
    if(length(mp)){
    scores<-scores[,unique(c("eid",unname(mp))),drop=FALSE]
    d<-dat;d$eid<-as.character(d$eid)
    d<-inner_join(d,scores,by="eid");rm(scores);invisible(gc())
    paired<-map_dfr(names(mp),function(x){
      z<-d[,unique(c(x,mp[[x]],tvar,evar,covars)),drop=FALSE];z<-z[complete.cases(z),,drop=FALSE]
      z<-z[z[[tvar]]>0&z[[evar]]%in%c(0,1),,drop=FALSE]
      z$.measured<-as.numeric(z[[x]]);z$.pgs<-as.numeric(z[[mp[[x]]]])
      singles<-cox_scan(z,c(".measured",".pgs"),covars,Y,time_var=tvar,event_var=evar)|>
        mutate(feature=x,model="separate; identical people and covariates",component=ifelse(term==".measured","Measured","PGS"),
          unit="per own SD; not interchangeable concentration units")
      if(nrow(z)>=500&&sum(z[[evar]]==1)>=20&&sd(z$.measured)>0&&sd(z$.pgs)>0){
        z$.measured<-as.numeric(scale(z$.measured));z$.pgs<-as.numeric(scale(z$.pgs))
        f<-as.formula(paste0("survival::Surv(",bt(tvar),",",bt(evar),") ~ ",paste(bt(c(".measured",".pgs",covars)),collapse=" + ")))
        fit<-tryCatch(survival::coxph(f,z),error=function(e)NULL)
        if(!is.null(fit)) {
          sm<-coef(summary(fit));j<-intersect(c(".measured",".pgs"),rownames(sm))
          joint<-tibble(term=j,beta=sm[j,"coef"],std.error=sm[j,"se(coef)"],p.value=sm[j,"Pr(>|z|)"],
            N_total=nrow(z),N_event=sum(z[[evar]]==1),feature=x,model="joint; identical people and covariates",
            component=ifelse(j==".measured","Measured","PGS"),unit="per own SD; not interchangeable concentration units")|>
            mutate(conf.low=exp(beta-1.96*std.error),conf.high=exp(beta+1.96*std.error),estimate=exp(beta))
          singles<-bind_rows(singles,joint)
        }
      }
      singles
    })|>mutate(FDR=p.adjust(p.value,"BH"))
    status<-tibble(status="ok",N_overlap=nrow(d),matched_anchors=length(mp),
      note="This paired comparison is not MR; significance is not compared across unequal sample sizes")
    }else status<-tibble(status="No configured anchor matched a PGS column")
  }
  write_raw_csv(paired,"c1.paired_pgs_measured.csv",rawdir)
  write_raw_csv(heterogeneity,"c1.temporal_heterogeneity.csv",rawdir)
  if(nrow(paired)) {
    p1<-paired|>mutate(lo=beta-1.96*std.error,hi=beta+1.96*std.error)|>
      ggplot(aes(beta,feature,color=component))+geom_vline(xintercept=0)+
      geom_errorbarh(aes(xmin=lo,xmax=hi),height=.1,position=position_dodge(width=.4))+
      geom_point(position=position_dodge(width=.4))+facet_wrap(~model)+
      labs(title="a. Paired measured and PGS associations",x="Log HR per own SD",y=NULL,color=NULL)+theme_5c(9)
  }else p1<-blank_plot("a. Paired measured and PGS associations","Matched PGS unavailable")
  p2<-heterogeneity|>ggplot(aes(-log10(pmax(p_time_interaction,1e-300)),feature))+
    geom_point()+labs(title="b. Does the baseline association vary across follow-up intervals?",x="-log10(time-interaction P)",y=NULL)+theme_5c(9)
  save_plot(p1/p2,"c1.Fig16.paired_temporal_validation.png",16,10,
    outdir=if(layer=="protein")out.prot else out.met)
  list(paired_associations=paired,time_heterogeneity=heterogeneity,paired_status=status)
}

# Loaded after the legacy C1 definitions, so same-day status is consistent.
make_prevalent_status <- function(dat,outcome=Y) {
  y<-as.Date(dat[[paste0("fod_icd10_",outcome)]]);b<-as.Date(dat$date_attend)
  ifelse(is.na(b),NA_real_,as.numeric(!is.na(y)&y<=b))
}

landmark_incident_scan <- function(dat,xs,tvar,evar,covars,landmarks=C1_LANDMARK_YEARS) {
  xs<-intersect(xs,names(dat));d0<-dat
  for(x in xs){z<-as.numeric(d0[[x]]);sd0<-sd(z,na.rm=TRUE)
    d0[[x]]<-if(is.finite(sd0)&&sd0>0)(z-mean(z,na.rm=TRUE))/sd0 else NA_real_}
  map_dfr(landmarks,function(h){
    d<-d0|>filter(is.finite(.data[[tvar]]),!is.na(.data[[evar]]),.data[[tvar]]>h)|>
      mutate(.landmark_time=.data[[tvar]]-h,.landmark_event=.data[[evar]])
    z<-cox_scan(d,xs,covars,Y,scale_x=FALSE,time_var=".landmark_time",event_var=".landmark_event")
    z|>mutate(landmark_years=h,N_risk=N_total,events_after_landmark=N_event,
      exposure_unit="SD fixed in the baseline omics cohort, not restandardized by landmark")
  })
}
