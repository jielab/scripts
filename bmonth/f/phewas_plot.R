suppressPackageStartupMessages(source(file.path(dir0,"scripts","0f","plotting.R")))
plot_bmonth_phewas_figures<-function(res,joint=NULL,cohort=NULL,phenotype_ids=NULL,xmap=NULL,only=NULL){
  res<-as_tibble(res)
  if(is.null(cohort))cohort<-readRDS(file.path(rawdir,"bmonth.cohort.rds"))
  if(is.null(phenotype_ids))phenotype_ids<-readRDS(file.path(indir,"Rdata","PheWAS.rds"))$id
  types<-c(
    bmonth_cos="continuous",bmonth_sin="continuous",
    bmonth_winter_summer="binary",bmonth_spring_autumn="binary",
    birth_north_100km="continuous",birth_east_100km="continuous",
    move_log1p="continuous",bmonth_rank="ordinal",
    border12="binary",border1n="binary",border_rank="ordinal",border_relative="continuous")
  if(is.null(xmap))xmap<-make_phewas_xmap(res)
  summaries<-samples<-list()
  for(x in unique(res$exposure)){
    d<-res[res$exposure==x,,drop=FALSE]
    cv<-phewas_covars_for(x)
    si<-phewas_sample_info(cohort,x,cv,phenotype_ids,type=types[[x]])
    is_context <- x %in% c("birth_north_100km","birth_east_100km")
    is_exploratory <- x %in% EXPLORATORY_ALL
    title <- if(is_exploratory) paste0(x," exploratory birth-order PheWAS") else
      if(is_context) paste0(x," birthplace-context PheWAS") else paste0(x," birth-month PheWAS")
    prefix <- if(is_exploratory) "FigS_exploratory_phewas_" else
      if(is_context) "FigS_context_phewas_" else "Fig1_phewas_"
    if(is.null(only)||x%in%only){
      p<-plot_phewas_results(d,title=title,sample_info=si,xmap=xmap)
      save_result_figure(p,paste0(prefix,x),d,width=13,height=6.5,dpi=320)
    }
    thr<-unique(d$bonferroni_threshold)
    summaries[[x]]<-d%>%filter(is.finite(p))%>%
      mutate(category=tolower(trimws(category)))%>%group_by(category)%>%
      summarise(n_tested=n(),n_bonf=sum(p<thr),.groups="drop")%>%
      mutate(exposure=x,analysis_track=analysis_track(x),covariate_model=phewas_model_label(x),.before=1)
    samples[[x]]<-tibble(exposure=x,analysis_track=analysis_track(x),covariate_model=phewas_model_label(x),
                         type=si$type,N=si$N,cases=si$cases,
                         controls=si$controls,title_counts=si$label)
  }
  if(length(summaries))write_csv2(bind_rows(summaries),"phewas_category_summary.csv")
  if(length(samples))write_csv2(bind_rows(samples),"phewas_plot_samples.csv")

  if(!is.null(joint)&&nrow(joint)){
    ref<-res%>%select(any_of(c("phenotype","category")))%>%distinct()
    jj<-as_tibble(joint)
    if(!"category"%in%names(jj))jj<-left_join(jj,ref,by="phenotype")
    jj$category[is.na(jj$category)]<-"other"
    for(x in unique(jj$exposure)){
      if(!is.null(only)&&!x%in%only)next
      d<-jj[jj$exposure==x,,drop=FALSE]
      si<-list(type="joint",N=max(d$N,na.rm=TRUE),cases=NA,controls=NA,
               label=paste0("N up to ",format(max(d$N,na.rm=TRUE),big.mark=",")))
      p<-plot_phewas_results(d,title=paste0(x," PheWAS"),sample_info=si,xmap=xmap)
      save_result_figure(p,paste0("Fig1_phewas_",x),d,width=13,height=6.5,dpi=320)
    }
  }
  invisible(NULL)
}
plot_phewas_figures<-plot_bmonth_phewas_figures
