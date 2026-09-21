# Endpoint ascertainment is separate from biomarker lead time. In particular,
# zero early AMR events must not be interpreted as biological latency.
le8_endpoint_audit<-function(dat,features,covars,Y,root){
  tv<-paste0(Y,".t2e");ev<-paste0(Y,".Yt2e");bv<-paste0(Y,".b2e")
  cuts<-c(0,.5,1,2,5,10,Inf)
  windows<-bind_rows(lapply(seq_len(length(cuts)-1L),function(i){
    at_risk<-is.finite(dat[[tv]])&dat[[tv]]>cuts[i]
    cases<-at_risk&dat[[ev]]==1&dat[[tv]]<=cuts[i+1]
    data.frame(window_lo=cuts[i],window_hi=cuts[i+1],N_at_risk=sum(at_risk),
      events=sum(cases,na.rm=TRUE),interpretation="Observed event distribution; verify recording coverage")
  }))
  write.csv(windows,file.path(root,"c1.endpoint_event_windows.csv"),row.names=FALSE)
  calendar<-if("date_attend"%in%names(dat))as.Date(dat$date_attend)else rep(as.Date(NA),nrow(dat))
  dy<-paste0("fod_icd10_",Y)
  if(dy%in%names(dat)){
    dates<-as.Date(dat[[dy]]);z<-table(format(dates[!is.na(dates)],"%Y"))
    write.csv(data.frame(year=names(z),recorded_diagnoses=as.integer(z)),file.path(root,"c1.endpoint_calendar_years.csv"),row.names=FALSE)
  }
  mf<-Sys.getenv("LE8_ENDPOINT_MANIFEST","")
  out<-data.frame(Y,prevalent=sum(is.finite(dat[[bv]])&dat[[bv]]<=0),
    early_2y_events=sum(dat[[ev]]==1&dat[[tv]]<=2,na.rm=TRUE),
    registry_status="unknown; absence of a diagnosis record is not confirmed absence of infection/resistance")
  if(nzchar(mf)){
    m<-read.csv(mf,stringsAsFactors=FALSE)
    required<-c("Y","source","definition","coverage_start","coverage_end","prebaseline_capture")
    if(!all(required%in%names(m)))stop("Endpoint manifest missing fields: ",paste(setdiff(required,names(m)),collapse=","))
    m<-m[m$Y==Y,,drop=FALSE]
    if(nrow(m)==1){
      write.csv(m,file.path(root,"c1.endpoint_manifest_used.csv"),row.names=FALSE)
      start<-as.Date(m$coverage_start);end<-as.Date(m$coverage_end)
      if(is.na(start)||is.na(end)||end<=start)stop("Registry dates must be valid, verified dates")
      d<-dat;d$.registry_entry<-pmax(0,as.numeric(start-calendar)/365.25)
      stop_time<-as.numeric(end-calendar)/365.25
      d$.registry_event<-as.integer(d[[ev]]==1&d[[tv]]<=stop_time)
      d$.registry_exit<-pmin(d[[tv]],stop_time)
      # Explicitly exclude any already observed disease at/before entry.
      eligible<-is.finite(d$.registry_entry)&is.finite(d$.registry_exit)&d$.registry_exit>d$.registry_entry
      d<-d[eligible,,drop=FALSE]
      z<-cox_scan_delayed_entry(d,features,covars,Y,entry_var=".registry_entry",exit_var=".registry_exit",event_var=".registry_event")
      z$time_scale<-"Years after blood draw, delayed entry at verified registry start"
      z$estimand<-"First recorded Y during covered follow-up; no claim of first lifetime disease if prebaseline capture incomplete"
      write.csv(z,file.path(root,"c1.registry_delayed_entry.csv"),row.names=FALSE)
      out$registry_status<-paste("Sensitivity completed;",m$source,"; prebaseline capture:",m$prebaseline_capture)
    }else out$registry_status<-"No unique manifest row for this Y; no registry correction assumed"
  }
  write.csv(out,file.path(root,"c1.endpoint_audit.csv"),row.names=FALSE)
  out
}
