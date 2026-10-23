plot_helper<-Sys.getenv("BMONTH_PLOT_HELPER",
  file.path(dirname(dirname(Sys.getenv("BMONTH_FDIR"))),"0f","plotting.R"))
if(!file.exists(plot_helper))stop("Missing repository plotting helper: ",plot_helper)
suppressPackageStartupMessages(source(plot_helper))
source(file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R"))

# A multi-df test has significance but no single signed beta. Keep beta=NA in
# exported results and use neutral hollow circles for joint-test discoveries.
plot_bmonth_joint<-function(d,title,si,xmap){
  p<-plot_phewas_results(d,title=title,sample_info=si,xmap=xmap)
  if(is.null(p))return(NULL)
  hit<-p$data[is.finite(p$data$p)&p$data$p<=unique(d$bonferroni_threshold),,drop=FALSE]
  if(nrow(hit))p<-p+ggplot2::geom_point(data=hit,
    ggplot2::aes(phe_x,logp_display,color=category),inherit.aes=FALSE,
    shape=21,fill="white",size=2.5,stroke=.9,show.legend=FALSE)
  p+ggplot2::labs(subtitle="Joint test: hollow circles mark Bonferroni discoveries; no single effect direction")
}
plot_bmonth_phewas_figures<-function(res,joint=NULL,cohort=NULL,phenotype_ids=NULL,xmap=NULL,only=NULL,provenance=NULL,cohort_md5=NULL){
  res<-as_tibble(res)
  if(is.null(cohort))cohort<-load_cohort()
  if(is.null(provenance))provenance<-phewas_provenance(phewas_input_path(validate=FALSE))
  if(is.null(cohort_md5))cohort_md5<-unname(tools::md5sum(file.path(rawdir,"bmonth.cohort.rds")))
  for(results in list(res,joint)){
    if(is.null(results)||!nrow(results))next
    required<-c("engine_version","phenotype_panel_md5","cohort_md5")
    if(length(setdiff(required,names(results))) ||
       anyNA(results[,required,drop=FALSE]) ||
       any(results$engine_version!=BM_ENGINE_VERSION) ||
       any(results$phenotype_panel_md5!=provenance$panel_md5) ||
       is.na(cohort_md5) || any(results$cohort_md5!=cohort_md5))
      stop("PheWAS plot results do not match the current engine/cohort/panel; rerun phewas")
  }
  if(is.null(phenotype_ids))phenotype_ids<-phewas_packed_panel(provenance)$ids
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
    if(!all(c("status","eligible","n_tests_planned") %in% names(d)))
      stop("Old PheWAS results cannot be redrawn as validated v2 results; rerun phewas")
    d$p[!grepl("^ok_",d$status)]<-NA_real_
    cv<-phewas_covars_for(x)
    si<-phewas_sample_info(cohort,x,cv,phenotype_ids,type=types[[x]])
    # These counts describe exposure groups, not disease cases and controls.
    labels<-list(border12=c("first-born","second-born"),
      border1n=c("oldest","youngest"),bmonth_winter_summer=c("winter","summer"),
      bmonth_spring_autumn=c("spring","autumn"))
    if(x %in% names(labels)) {
      lab<-labels[[x]]
      si$label<-paste0(format(si$cases,big.mark=",")," ",lab[2]," vs ",
                        format(si$controls,big.mark=",")," ",lab[1],
                        "; phenotype-specific disease N varies")
    } else si$label<-paste0(si$label,"; phenotype-specific disease N varies")
    si$label<-paste0(si$label,"; eligible ",unique(d$n_tests_planned),
                     ", unresolved ",unique(d$n_tests_unresolved))
    is_context <- x %in% c("birth_north_100km","birth_east_100km")
    is_exploratory <- x %in% EXPLORATORY_ALL
    title <- if(is_exploratory) paste0(x," exploratory birth-order PheWAS") else
      if(is_context) paste0(x," birthplace-context PheWAS") else paste0(x," birth-month PheWAS")
    prefix <- if(is_exploratory) "FigS_exploratory_phewas_" else
      if(is_context) "FigS_context_phewas_" else "Fig1_phewas_"
    if(is.null(only)||x%in%only){
      p<-plot_phewas_results(d,title=title,sample_info=si,xmap=xmap)
      if(!is.null(p))p<-p+ggplot2::labs(title=paste(strwrap(p$labels$title,width=95),collapse="\n"))
      save_result_figure(p,paste0(prefix,x),d,width=13,height=6.5,dpi=320)
    }
    thr<-unique(d$bonferroni_threshold)
    summaries[[x]]<-d%>%
      mutate(category=tolower(trimws(category)))%>%group_by(category)%>%
      summarise(n_eligible=sum(eligible),n_tested=sum(is.finite(p)),
                n_unresolved=sum(eligible & !is.finite(p)),
                n_bonf=sum(p<thr,na.rm=TRUE),.groups="drop")%>%
      mutate(exposure=x,analysis_track=analysis_track(x),covariate_model=phewas_model_label(x),.before=1)
    samples[[x]]<-tibble(exposure=x,analysis_track=analysis_track(x),covariate_model=phewas_model_label(x),
                         type=si$type,N=si$N,exposure_1_n=si$cases,
                         exposure_0_n=si$controls,title_counts=si$label)
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
      if(!all(c("status","eligible","n_tests_planned") %in% names(d)))
        stop("Old joint PheWAS results need rerunning before v2 plots")
      d$p[!grepl("^ok_",d$status)]<-NA_real_
      si<-list(type="joint",N=max(d$N,na.rm=TRUE),cases=NA,controls=NA,
               label=paste0("N up to ",format(max(d$N,na.rm=TRUE),big.mark=","),
                 "; eligible ",unique(d$n_tests_planned),", unresolved ",unique(d$n_tests_unresolved)))
      p<-plot_bmonth_joint(d,title=paste0(x," PheWAS"),si=si,xmap=xmap)
      if(!is.null(p))p<-p+ggplot2::labs(title=paste(strwrap(p$labels$title,width=95),collapse="\n"))
      save_result_figure(p,paste0("Fig1_phewas_",x),d,width=13,height=6.5,dpi=320)
    }
  }
  invisible(NULL)
}
plot_phewas_figures<-plot_bmonth_phewas_figures
