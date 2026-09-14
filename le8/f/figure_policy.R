# Presentation policy: retain all available views and source layouts. Grouping
# and pagination never authorize deletion of a scientific view or old version.
if(!exists(".le8_figure_queue",inherits=FALSE)).le8_figure_queue<-new.env(parent=emptyenv())

le8_figure_rule<-function(file) {
  key<-sub("[.]png$","",sub("^[^.]+[.]Fig[0-9]+[.]","",basename(file)))
  module<-sub("[.].*$","",basename(file))
  omit<-character()
  group<-key
  groups<-switch(module,
    c1=list(temporal_profiles=c("yy_top","gradient"),temporal_evidence=c("directionality_triage","directionality_detail")),
    c2=list(instrument_diagnostics=c("pQTL_R2"),dandelion_sensitivity=c("dandelion","dandelion_evidence")),
    c3=list(posterior_evidence=c("evidence_triage","posterior_diagnostics","prior_sensitivity")),
    c4=list(supervised_connections=c("group_pillar_flow","supervised_atlas"),mediation=c("mediation_forest","mediation_diagnostics"),state_remodeling=c("state_network_remodeling","state_network_edges"),imaging_context=c("imaging_associations","imaging_overview")),
    c5=list(prediction_benchmark=c("performance_benchmark","incremental_performance","parsimony_performance"),prediction_sensitivity=c("subgroup_discrimination","attained_age_sensitivity")),list())
  for(nm in names(groups))if(key%in%groups[[nm]])group<-nm
  list(module=module,key=key,group=group,omit=key%in%omit)
}

le8_queue_figure<-function(p,target,w,h,dpi) {
  rd<-dirname(target);q<-.le8_figure_queue[[rd]]
  if(is.null(q))q<-list()
  rule<-le8_figure_rule(target)
  q[[basename(target)]]<-list(plot=if(rule$omit)NULL else p,rule=rule,width=w,height=h,dpi=dpi)
  .le8_figure_queue[[rd]]<-q
  invisible(target)
}

le8_plot_leaves<-function(p) {
  if(is.null(p))return(list())
  if(inherits(p,"patchwork"))return(unlist(lapply(seq_len(length(p)),function(i)le8_plot_leaves(p[[i]])),recursive=FALSE))
  if(inherits(p,"spacer"))return(list())
  list(p)
}

le8_clone_layer<-function(x){force(x);ggplot2::ggproto(NULL,x)}

# Shared themes precede layer-specific additions. Presentation numbers are
# independent of legacy source-file numbers and remain consecutive.
le8_figure_order<-function(module,groups) {
  preferred<-switch(module,
    c1=c('mh','circular','vc','temporal_profiles','diagnosis_timed_profiles',
      'gradient_cluster','quantile_top','enrich_sig','temporal_evidence',
      'landmark_birthline_sensitivity','diagnosis_window_riskset','paired_temporal_validation'),
    c2=c('instrument_diagnostics','effect_concordance','mr_incident_prevalent','mrlink2','bidirectional_mr',
      'genetic_decomposition','genetic_component_leadtime','evidence_grades',
      'dandelion_sensitivity','dandelion_mr_integration'),
    c4=c('supervised_connections','mediation','imaging_context','lifestyle_omics_risk','state_remodeling'),
    c5=c('prediction_benchmark','score_concordance','prediction_sensitivity',
      'budget_ablation_calibration','heldout_ROC','leadtime_prediction'),character())
  order(match(groups,preferred,nomatch=length(preferred)+1L),seq_along(groups))
}

# Version old presentation files before any overwrite or renumber. Versions are
# content-addressed, so an unchanged PNG is not copied on every run.
le8_archive_figure_files<-function(files,root) {
  files<-files[file.exists(files)];if(!length(files))return(invisible(character()))
  dest<-file.path(root,'_history');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  hashes<-unname(tools::md5sum(files))
  targets<-file.path(dest,paste0(hashes,'__',basename(files)))
  take<-!file.exists(targets)
  if(any(take)&&!all(file.copy(files[take],targets[take])))stop('Cannot preserve previous figure versions')
  invisible(targets)
}

# Merge freshly rendered themes with existing PNGs, then renumber via the same
# policy as a full module run. This permits aggregate-only presentation updates.
le8_refresh_figure_files<-function(rawdir,incoming=NULL) {
  mf<-file.path(rawdir,'figure_manifest.csv')
  old<-if(file.exists(mf))as.data.frame(data.table::fread(mf))else data.frame()
  fresh<-if(!is.null(incoming))as.data.frame(data.table::fread(file.path(incoming,'figure_manifest.csv')))else old[FALSE,]
  retained<-if(nrow(old))old[!old$group%in%fresh$group,,drop=FALSE]else old
  revised<-rbind(retained,fresh)
  if(!nrow(revised)) {
    data.table::fwrite(revised,mf)
    if(!is.null(incoming))stopifnot(file.copy(file.path(incoming,'figure_omission_audit.csv'),
      file.path(rawdir,'figure_omission_audit.csv'),overwrite=TRUE))
    return(invisible(revised))
  }
  sources<-c(file.path(rawdir,retained$file),if(nrow(fresh))file.path(incoming,fresh$file)else character())
  stopifnot(all(file.exists(sources)))
  prefix<-sub('[.].*$','',revised$file[1])
  ix<-le8_figure_order(prefix,revised$group);revised<-revised[ix,,drop=FALSE];sources<-sources[ix]
  revised$file<-paste0(prefix,'.Fig',seq_len(nrow(revised)),'.',revised$group,'.png')
  stage<-tempfile('.renumber-',tmpdir=rawdir);dir.create(stage)
  on.exit(unlink(stage,recursive=TRUE),add=TRUE)
  stopifnot(all(file.copy(sources,file.path(stage,revised$file))))
  le8_archive_figure_files(list.files(rawdir,pattern='[.]Fig[0-9]+[.].*[.]png$',full.names=TRUE),rawdir)
  le8_archive_figure_files(c(mf,file.path(rawdir,'figure_omission_audit.csv')),rawdir)
  # Stage every source before overwriting names that may be reused by another theme.
  stopifnot(all(file.copy(file.path(stage,revised$file),file.path(rawdir,revised$file),overwrite=TRUE)))
  stale<-if(nrow(old))setdiff(old$file,revised$file)else character()
  if(length(stale))unlink(file.path(rawdir,stale))
  data.table::fwrite(revised,mf)
  if(!is.null(incoming)) {
    af<-file.path(rawdir,'figure_omission_audit.csv')
    a<-if(file.exists(af))as.data.frame(data.table::fread(af))else data.frame()
    b<-as.data.frame(data.table::fread(file.path(incoming,'figure_omission_audit.csv')))
    if(nrow(a))a<-a[!a$group%in%b$group,,drop=FALSE]
    data.table::fwrite(rbind(a,b),af)
  }
  invisible(revised)
}

le8_layer_figure_correspondence<-function(traitdir) {
  read_layer<-function(layer) {
    files<-list.files(file.path(traitdir,layer),pattern='^figure_manifest[.]csv$',recursive=TRUE,full.names=TRUE)
    files<-files[!grepl('/_history/',files,fixed=TRUE)]
    ans<-lapply(files,function(f) {
      d<-as.data.frame(data.table::fread(f));if(!nrow(d))return(NULL)
      d$module<-basename(dirname(f));d$page<-ave(seq_len(nrow(d)),d$group,FUN=seq_along)
      # Genomic protein order and metabolite biochemical order share one theme.
      d$group[d$group%in%c('mh','circular')]<-'association_overview'
      d$file<-file.path(d$module,d$file);d[,c('module','group','page','file','panels')]
    })
    ans<-do.call(rbind,ans)
    if(is.null(ans))ans<-data.frame(module=character(),group=character(),page=integer(),file=character(),panels=integer())
    names(ans)[4:5]<-paste0(layer,c('_file','_panels'));ans
  }
  d<-merge(read_layer('prot'),read_layer('met'),by=c('module','group','page'),all=TRUE)
  d$status<-ifelse(!is.na(d$prot_file)&!is.na(d$met_file),'paired theme',
    ifelse(is.na(d$prot_file),'met only / no usable prot panel','prot only / no usable met panel'))
  number<-function(x)as.integer(sub('.*[.]Fig([0-9]+).*','\\1',x))
  d$same_number<-ifelse(d$status=='paired theme',number(d$prot_file)==number(d$met_file),NA)
  data.table::fwrite(d,file.path(traitdir,'figure_layer_correspondence.csv'))
  invisible(d)
}

le8_readable_plot<-function(p) {
  # Limit labels only; every observation remains in points/lines and tables.
  for(i in seq_along(p$layers))if(inherits(p$layers[[i]]$geom,"GeomTextRepel")||inherits(p$layers[[i]]$geom,"GeomLabelRepel")) {
    layer<-le8_clone_layer(p$layers[[i]])
    d<-if(is.data.frame(layer$data))layer$data else p$data
    label<-layer$mapping$label %||% p$mapping$label
    if(is.data.frame(d)&&!is.null(label)) {
      labs<-tryCatch(rlang::eval_tidy(label,data=d),error=function(e)rep(NA_character_,nrow(d)))
      keep<-head(which(!is.na(labs)&nzchar(as.character(labs))),12L)
      layer$data<-d[keep,,drop=FALSE]
      layer$aes_params$colour<-"grey20";layer$geom_params$max.overlaps<-12L
      p$layers[[i]]<-layer
    }
  }
  for(nm in c("title","subtitle","caption"))if(is.character(p$labels[[nm]])&&length(p$labels[[nm]])==1) {
    value<-p$labels[[nm]]
    if(nm=="title")value<-sub("^[A-Za-z][.] +","",value)
    p<-p+do.call(ggplot2::labs,setNames(list(stringr::str_wrap(value,if(nm=="title")60 else 85)),nm))
  }
  p+ggplot2::theme(plot.title=ggplot2::element_text(size=11,face="bold"),
    plot.subtitle=ggplot2::element_text(size=9),plot.caption=ggplot2::element_text(size=8,hjust=0),
    axis.text=ggplot2::element_text(size=9),axis.title=ggplot2::element_text(size=10),
    legend.text=ggplot2::element_text(size=8),legend.title=ggplot2::element_text(size=9),
    legend.position="bottom",legend.box="vertical",legend.box.just="center",
    plot.margin=ggplot2::margin(10,15,10,12))
}

le8_expand_facets<-function(p) {
  if(!is.null(attr(p,"le8_unavailable")))return(list())
  if(inherits(p$facet,"FacetNull"))return(list(le8_readable_plot(p)))
  layout<-ggplot2::ggplot_build(p)$layout$layout
  vars<-setdiff(names(layout),c("PANEL","ROW","COL","SCALE_X","SCALE_Y","COORD"))
  if(nrow(layout)<=1||!length(vars))return(list(le8_readable_plot(p)))
  # Each facet becomes an independent axis, so the six-panel limit is real.
  lapply(seq_len(nrow(layout)),function(j) {
    values<-layout[j,vars,drop=FALSE]
    subset_data<-function(d) {
      if(!is.data.frame(d))return(d)
      common<-intersect(names(d),vars);ok<-rep(TRUE,nrow(d))
      for(v in common)ok<-ok & (if(is.na(values[[v]]))is.na(d[[v]]) else !is.na(d[[v]])&as.character(d[[v]])==as.character(values[[v]]))
      droplevels(d[ok,,drop=FALSE])
    }
    z<-p;z$data<-subset_data(z$data)
    for(i in seq_along(z$layers)){layer<-le8_clone_layer(z$layers[[i]]);layer$data<-subset_data(layer$data);z$layers[[i]]<-layer}
    z<-z+ggplot2::facet_null()+ggplot2::labs(title=paste(p$labels$title,paste(unlist(values),collapse=" | "),sep=" — "))
    le8_readable_plot(z)
  })
}

le8_figure_composition<-function(plots,group) {
  n<-length(plots);nc<-if(n==1)1 else 2;nr<-ceiling(n/nc)
  if(group=="diagnosis_timed_profiles"&&n==4L) {
    # Restore the original c1.MOCK.Fig2 geometry: one full-height heatmap
    # alongside three aligned cluster axes. Do not flatten this into a 2x2.
    return(list(plot=patchwork::wrap_plots(plots,design="AB\nAC\nAD",widths=c(1.3,1),heights=c(1,1,1)),
      width=22,height=12,policy="full-height trajectory heatmap + three aligned cluster axes; all numeric rows retained"))
  }
  if(group=='circular')return(list(
    plot=patchwork::wrap_plots(plots,ncol=nc,guides='collect') & ggplot2::theme(legend.position='bottom'),
    width=if(nc==1)12 else 24,height=nr*12+.8,
    policy='biochemical-category radial associations; shared effect scale; <=6 axes; all significant features labelled'))
  if(group=='effect_concordance'&&n==2L)return(list(
    plot=patchwork::wrap_plots(plots,ncol=2,widths=c(1.05,1)),width=22,height=10,
    policy='paired genetic effects with 95% CIs; fixed-anchor forest; no fitted regression line'))
  if(group=='mr_incident_prevalent'&&n==5L)return(list(
    plot=patchwork::wrap_plots(plots,design='AAABBB\nCCDDEE',heights=c(.55,1.5)),width=22,height=14,
    policy='two support summaries above three row-aligned MR/incident/prevalent forests; five axes'))
  if(group=='lifestyle_omics_risk'&&n==1L)return(list(plot=plots[[1]],width=20,height=14,
    policy='LE8-selected proxy association ribbons; original nodes and signed links retained; one axis'))
  if(group=='prots.top'&&n==6L)return(list(
    plot=patchwork::wrap_plots(plots,ncol=3,heights=c(.65,1.6)),width=22,height=17,
    policy='restored three overlap panels above three aligned effect forests; six axes'))
  if(group=='directionality_causal'&&n==3L)return(list(
    plot=patchwork::wrap_plots(plots,design='AA\nBC',heights=c(1.25,1)),width=20,height=15,
    policy='restored full-width temporal/genetic evidence matrix above two diagnostic views'))
  if(group%in%c('mrlink2','genetic_decomposition')&&n==3L)return(list(
    plot=patchwork::wrap_plots(plots,design='AB\nAC',widths=c(1,1.1)),width=20,height=12,
    policy='full-height estimate panel + two aligned diagnostic panels; three axes'))
  list(plot=patchwork::wrap_plots(plots,ncol=nc,widths=rep(1,nc),heights=rep(1,nr)),
    width=if(nc==1)10 else 19,height=nr*6+.6,
    policy="equal panel cells; <=6 axes; max 12 labels in source order; all numeric rows retained")
}

le8_flush_figures<-function(rawdir) {
  q<-.le8_figure_queue[[rawdir]];if(is.null(q))return(invisible(NULL))
  nms<-names(q);ord<-order(as.integer(sub(".*[.]Fig([0-9]+).*","\\1",nms)),nms,na.last=TRUE)
  q<-q[ord];groups<-list();audit<-list();notes<-list()
  for(src in names(q)) {
    item<-q[[src]];rule<-item$rule
    audit[[src]]<-data.frame(source=src,group=rule$group,status=if(rule$omit)"omitted: redundant, obsolete or misleading view"else"queued")
    if(rule$omit)next
    leaves<-le8_plot_leaves(item$plot)
    # Keep the approved first six C1 figures stable while retaining the old
    # cluster diagnostics on additional pages, rather than discarding them.
    extra<-if(rule$key=='gradient_cluster'&&length(leaves)>2L)leaves[-c(1,2)]else list()
    if(length(extra))leaves<-leaves[1:2]
    missing<-vapply(leaves,function(p)!is.null(attr(p,"le8_unavailable")),logical(1))
    if(any(missing))notes[[rule$group]]<-c(notes[[rule$group]],vapply(leaves[missing],attr,character(1),"le8_unavailable"))
    pp<-unlist(lapply(leaves[!missing],le8_expand_facets),recursive=FALSE)
    if(!length(pp)){audit[[src]]$status<-"omitted: no usable panel";next}
    caption<-if(inherits(item$plot,"patchwork"))item$plot$patches$annotation$caption else NULL
    if(is.character(caption))notes[[rule$group]]<-c(notes[[rule$group]],caption)
    if(rule$group%in%c("temporal_profiles","gradient_cluster","diagnosis_window_riskset"))
      notes[[rule$group]]<-c(notes[[rule$group]],"Each participant contributes one baseline omic measurement. Diagnosis-time bins compare different people; these are not within-person longitudinal trajectories.")
    for(p in pp)groups[[rule$group]]<-c(groups[[rule$group]],list(list(plot=p,source=src)))
    if(length(extra)) {
      diagnostics<-unlist(lapply(extra,le8_expand_facets),recursive=FALSE)
      for(p in diagnostics)groups[['gradient_cluster_diagnostics']]<-c(groups[['gradient_cluster_diagnostics']],list(list(plot=p,source=src)))
      notes[['gradient_cluster_diagnostics']]<-'Restored cluster diagnostics are descriptive heuristics; they do not establish causal direction.'
    }
    audit[[src]]$status<-"rendered in grouped figures"
  }
  staging<-tempfile(".figures-",tmpdir=rawdir);dir.create(staging)
  on.exit(unlink(staging,recursive=TRUE),add=TRUE)
  manifest<-list();number<-0L;prefix<-q[[1]]$rule$module
  dpi<-as.integer(Sys.getenv("LE8_FIGURE_DPI",unset="220"))
  for(group in names(groups)[le8_figure_order(prefix,names(groups))]) {
    pp<-groups[[group]]
    for(start in seq(1L,length(pp),by=6L)) {
      page<-pp[start:min(length(pp),start+5L)];number<-number+1L
      file<-sprintf("%s.Fig%d.%s.png",prefix,number,group)
      n<-length(page);composition<-le8_figure_composition(lapply(page,`[[`,"plot"),group)
      caption<-paste(unique(notes[[group]]),collapse="; ")
      p<-composition$plot+
        patchwork::plot_annotation(caption=if(nzchar(caption))stringr::str_wrap(caption,150)else NULL,
          tag_levels="A",theme=ggplot2::theme(plot.caption=ggplot2::element_text(size=9,hjust=0),plot.tag=ggplot2::element_text(face="bold")))
      ggplot2::ggsave(file.path(staging,file),p,width=composition$width,height=composition$height,dpi=dpi,bg="white",limitsize=FALSE)
      manifest[[length(manifest)+1L]]<-data.frame(file=file,group=group,panels=n,sources=paste(unique(vapply(page,`[[`,character(1),"source")),collapse=";"),
        policy=composition$policy)
    }
  }
  # Preserve source geometry separately from the aligned/paginated figures.
  originals<-file.path(rawdir,'_source_figures');dir.create(originals,showWarnings=FALSE)
  for(src in names(q)) {
    item<-q[[src]];target<-file.path(originals,src)
    if(is.null(item$plot))next
    le8_archive_figure_files(target,originals)
    tryCatch(ggplot2::ggsave(target,item$plot,width=item$width,height=item$height,
      dpi=dpi,bg='white',limitsize=FALSE),error=function(e) {
        audit[[src]]$status<-paste(audit[[src]]$status,'; source-layout rendering failed:',conditionMessage(e))
      })
  }
  mf<-if(length(manifest))do.call(rbind,manifest)else data.frame(file=character(),group=character(),panels=integer(),sources=character(),policy=character())
  data.table::fwrite(mf,file.path(staging,"figure_manifest.csv"))
  data.table::fwrite(do.call(rbind,audit),file.path(staging,"figure_omission_audit.csv"))
  mf<-le8_refresh_figure_files(rawdir,staging)
  .le8_figure_queue[[rawdir]]<-NULL
  invisible(mf)
}
