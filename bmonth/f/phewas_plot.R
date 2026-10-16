suppressPackageStartupMessages(source(file.path(dir0,"scripts","0f","plotting.R")))
plot_bmonth_phewas_figures<-function(res,joint=NULL,cohort=NULL,phenotype_ids=NULL,xmap=NULL,only=NULL){
  res<-as_tibble(res)
  if(is.null(cohort))cohort<-readRDS(file.path(rawdir,"bmonth.cohort.rds"))
  if(is.null(phenotype_ids))phenotype_ids<-readRDS(file.path(indir,"Rdata","PheWAS.rds"))$id
  types<-c(
    bmonth_cos="continuous",bmonth_sin="continuous",
    bmonth_winter_summer="binary",bmonth_spring_autumn="binary",
    birth_north_100km="continuous",birth_east_100km="continuous",
    move_log1p="continuous",bmonth_rank="ordinal")
  if(is.null(xmap))xmap<-make_phewas_xmap(res)
  summaries<-samples<-list()
  for(x in unique(res$exposure)){
    d<-res[res$exposure==x,,drop=FALSE]
    si<-phewas_sample_info(cohort,x,covars,phenotype_ids,type=types[[x]])
    if(is.null(only)||x%in%only){
      p<-plot_phewas_results(d,title=paste0(x," PheWAS"),sample_info=si,xmap=xmap)
      if(!is.null(p))ggplot2::ggsave(file.path(figdir,paste0("Fig1_phewas_",x,".png")),
                                    p,width=13,height=6.5,dpi=320,bg="white")
    }
    thr<-unique(d$bonferroni_threshold)
    summaries[[x]]<-d%>%filter(is.finite(p))%>%
      mutate(category=tolower(trimws(category)))%>%group_by(category)%>%
      summarise(n_tested=n(),n_bonf=sum(p<thr),.groups="drop")%>%
      mutate(exposure=x,.before=1)
    samples[[x]]<-tibble(exposure=x,type=si$type,N=si$N,cases=si$cases,
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
      if(!is.null(p))ggplot2::ggsave(file.path(figdir,paste0("Fig1_phewas_",x,".png")),
                                    p,width=13,height=6.5,dpi=320,bg="white")
    }
  }
  invisible(NULL)
}
plot_phewas_figures<-plot_bmonth_phewas_figures
