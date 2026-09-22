# Aggregate-only figures. Every rendered page has a companion workbook.
pgs_theme <- function() ggplot2::theme_classic(base_size=11)+ggplot2::theme(
  plot.title=ggplot2::element_text(face="bold",size=12),plot.subtitle=ggplot2::element_text(size=9,color="#526172"),
  legend.position="bottom",legend.title=ggplot2::element_blank(),plot.margin=ggplot2::margin(10,12,10,10))
pgs_palette<-c("Both supported: opposite"="#B34359","Both supported: concordant"="#287D8E",
  "Measured only"="#B39557","PGS only"="#7964A5","Neither supported"="#C1C8CE","Unavailable comparison"="#E0E3E5")
pgs_ok<-function(d,cols) is.data.frame(d)&&nrow(d)>0&&all(cols%in%names(d))
pgs_pick<-function(x,n=8) {
  if(!pgs_ok(x$paired,c("feature","evidence_pattern")))return(character())
  d<-x$paired;anchor<-if(pgs_ok(x$selection,c("feature","selection")))x$selection$feature[x$selection$selection=="declared anchor"]else character()
  ranked<-d$feature[order(d$evidence_pattern!="Both supported: opposite",pmax(d$measured_FDR,d$pgs_FDR),d$measured_p,na.last=TRUE)]
  # Ranked discordants plus declared biological controls; selection is exploratory.
  unique(head(c(head(ranked,n%/%2),anchor,ranked),n))
}
pgs_forest_plot<-function(d,title,subtitle="",unit="Log HR",label="feature",color="series") {
  if(!pgs_ok(d,c(label,"beta","lo","hi",color)))return(NULL)
  d<-d[is.finite(d$beta)&is.finite(d$lo)&is.finite(d$hi),,drop=FALSE];if(!nrow(d))return(NULL)
  d$.label<-factor(d[[label]],levels=rev(unique(d[[label]])));d$.series<-d[[color]]
  ggplot2::ggplot(d,ggplot2::aes(beta,.label,color=.series))+ggplot2::geom_vline(xintercept=0,linetype=2,color="grey65")+
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo,xmax=hi),height=.15,position=ggplot2::position_dodge(.55))+
    ggplot2::geom_point(position=ggplot2::position_dodge(.55),size=2)+
    ggplot2::scale_color_manual(values=rep(c("#287D8E","#B34359","#7964A5","#B39557","#526172","#59A596","#A57761","#555555"),length.out=max(1,length(unique(d$.series)))))+
    ggplot2::labs(title=title,subtitle=subtitle,x=unit,y=NULL,color=NULL)+pgs_theme()
}
pgs_locus_panels<-function(features,loci) {
  pp<-list()
  if(pgs_ok(features,c("observed_beta","pgs_beta","triangulation_class"))) {
    a<-features[is.finite(features$observed_beta)&is.finite(features$pgs_beta),,drop=FALSE]
    if(nrow(a))pp$comparison<-ggplot2::ggplot(a,ggplot2::aes(pgs_beta,observed_beta,color=triangulation_class))+
      ggplot2::geom_hline(yintercept=0,color="grey80")+ggplot2::geom_vline(xintercept=0,color="grey80")+
      ggplot2::geom_point(alpha=.65)+ggplot2::scale_color_manual(values=pgs_palette)+
      ggplot2::labs(title="Measured and PGS associations",subtitle=if(all(a$matched_comparison))"Identical participants and adjustment"else"Contains legacy unmatched screening; rerun pgs_focus",
        x="PGS log HR / own SD",y="Measured log HR / own SD",color=NULL)+pgs_theme()
  }
  if(pgs_ok(loci,c("feature","locus","locus_class","PP.H4","PP.H4_robust_min","PP.H3"))) {
    rank<-match(loci$feature,features$feature);a<-head(loci[order(rank,-loci$PP.H4,na.last=TRUE),],14)
    a$.label<-paste(a$feature,a$locus_class,a$locus,sep=" | ")
    z<-pgs_bind(lapply(c("PP.H3","PP.H4","PP.H4_robust_min"),function(n)data.frame(label=a$.label,metric=n,posterior=a[[n]])))
    z$label<-factor(z$label,levels=rev(unique(a$.label)))
    pp$loci<-ggplot2::ggplot(z,ggplot2::aes(metric,label,fill=posterior))+ggplot2::geom_tile(color="white")+
      ggplot2::geom_text(ggplot2::aes(label=ifelse(is.finite(posterior),sprintf("%.2f",posterior),"NA")),size=2.7)+
      ggplot2::scale_fill_gradient(low="#F0F4F5",high="#247F8C",limits=c(0,1),na.value="#DEDEDE")+
      ggplot2::scale_x_discrete(labels=c("PP.H3"="H3","PP.H4"="H4 default","PP.H4_robust_min"="H4 min prior"))+
      ggplot2::labs(title="Locus support and prior sensitivity",subtitle="All regions retained in workbook; NA means not tested",x=NULL,y=NULL,fill=NULL)+pgs_theme()+
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=7),legend.position="none")
    z<-unique(loci[,c("region_cluster","cluster_features","cluster_robust_features")]);z<-head(z[order(-z$cluster_robust_features,-z$cluster_features),],10)
    pp$clusters<-ggplot2::ggplot(z,ggplot2::aes(cluster_robust_features,reorder(region_cluster,cluster_robust_features)))+
      ggplot2::geom_col(fill="#7964A5",width=.7)+ggplot2::labs(title="Shared-region concentration",subtitle="Overlapping windows, not independent causal signals",x="Biomarkers with robust H4",y=NULL)+pgs_theme()
  }
  pp
}
pgs_main_panels<-function(x,tri=data.frame(),loci=data.frame()) {
  pp<-list();chosen<-pgs_pick(x);d<-x$paired
  if(pgs_ok(d,c("feature","pgs_beta","measured_beta","evidence_pattern"))) {
    a<-d[is.finite(d$pgs_beta)&is.finite(d$measured_beta),,drop=FALSE];a$label<-ifelse(a$feature%in%chosen,a$feature,NA_character_)
    if(nrow(a))pp$paired<-ggplot2::ggplot(a,ggplot2::aes(pgs_beta,measured_beta,color=evidence_pattern))+
      ggplot2::geom_hline(yintercept=0,color="grey80")+ggplot2::geom_vline(xintercept=0,color="grey80")+
      ggplot2::geom_point(alpha=.65,size=1.7)+ggrepel::geom_text_repel(ggplot2::aes(label=label),size=2.8,seed=2026,max.overlaps=20,na.rm=TRUE)+
      ggplot2::scale_color_manual(values=pgs_palette)+ggplot2::labs(title="Measured level versus biomarker PGS",
        subtitle="Identical people and covariates; BH FDR within each full scan",x="PGS log HR / own SD",y="Measured log HR / own SD",color=NULL)+pgs_theme()
  }
  a<-x$matched_models
  if(pgs_ok(a,c("feature","model","beta","lo","hi"))) {
    a<-a[a$feature%in%chosen&a$model%in%c("measured","pgs"),];a$series<-ifelse(a$model=="measured","Measured","Biomarker PGS")
    pp$effects<-pgs_forest_plot(a,"Direction and uncertainty","Separate models, identical sample","Log HR per own SD")
  }
  a<-x$components
  if(pgs_ok(a,c("feature","term","model","beta","lo","hi"))) {
    a<-a[a$feature%in%chosen&a$model=="joint",];a$series<-ifelse(a$term==".G","PGS-captured G","Remaining R")
    pp$components<-pgs_forest_plot(a,"Captured and remaining components","Joint model; conditional 95% CI; refitted bootstrap in workbook","Log HR per whole-biomarker SD")
  }
  a<-x$calibration
  if(pgs_ok(a,c("feature","partial_R2","weak_capture"))) {
    a<-a[a$feature%in%chosen&is.finite(a$partial_R2),]
    if(nrow(a))pp$calibration<-ggplot2::ggplot(a,ggplot2::aes(partial_R2,reorder(feature,partial_R2),color=weak_capture))+
      ggplot2::geom_vline(xintercept=0,color="grey75")+ggplot2::geom_point(size=2.7)+
      ggplot2::scale_color_manual(values=c("FALSE"="#287D8E","TRUE"="#B39557"),labels=c("Adequate capture","Weak capture"))+
      ggplot2::labs(title="Does PGS capture the measured biomarker?",subtitle="Held-out incremental prediction; negative values retained",x="Cross-fitted partial R²",y=NULL,color=NULL)+pgs_theme()
  }
  a<-x$window_contrasts
  if(pgs_ok(a,c("feature","beta_difference","lo","hi","landmark","end"))) {
    a<-a[a$feature%in%head(chosen,4)&is.finite(a$beta_difference),]
    if(nrow(a)) {
      a$window<-paste0(a$landmark,"–",ifelse(is.finite(a$end),a$end,"end"));a$window<-factor(a$window,levels=unique(a$window[order(a$landmark)]))
      pp$time<-ggplot2::ggplot(a,ggplot2::aes(window,beta_difference,color=feature,group=feature))+
        ggplot2::geom_hline(yintercept=0,linetype=2,color="grey65")+ggplot2::geom_line()+
        ggplot2::geom_errorbar(ggplot2::aes(ymin=lo,ymax=hi),width=.12,position=ggplot2::position_dodge(.25))+
        ggplot2::geom_point(position=ggplot2::position_dodge(.25))+
        ggplot2::labs(title="Does discordance persist with follow-up?",subtitle="Risk-set intervals; joint contrast includes G/R covariance",x="Years since baseline",y="βG − βR (conditional 95% CI)",color=NULL)+pgs_theme()
    }
  }
  if(nrow(loci)) {
    # Order by the same declared figure candidates, not a second best-H4 ranking.
    tr<-if(nrow(tri))tri[order(match(tri$feature,chosen),na.last=TRUE),,drop=FALSE]else data.frame(feature=chosen)
    lp<-pgs_locus_panels(tr,loci);pp$loci<-lp$loci
  }
  Filter(Negate(is.null),pp)
}
pgs_save_panels<-function(pp,rd,stem,tables,title,caption="") {
  pp<-Filter(Negate(is.null),pp);if(!length(pp))return(invisible(NULL))
  dir.create(rd,recursive=TRUE,showWarnings=FALSE);manifest<-list()
  for(start in seq(1,length(pp),by=6L)) {
    page<-pp[start:min(length(pp),start+5L)];name<-if(length(pp)<=6)stem else paste0(stem,".page",ceiling(start/6))
    p<-patchwork::wrap_plots(page,ncol=2)+patchwork::plot_annotation(title=title,
      caption=stringr::str_wrap(caption,140),tag_levels="A",theme=ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=17),plot.caption=ggplot2::element_text(hjust=0,size=9)))
    height<-ceiling(length(page)/2)*4.7+1
    ggplot2::ggsave(file.path(rd,paste0(name,".png")),p,width=17,height=height,dpi=as.integer(Sys.getenv("LE8_FINAL_DPI","300")),bg="white",limitsize=FALSE)
    ggplot2::ggsave(file.path(rd,paste0(name,".pdf")),p,width=17,height=height,bg="white",limitsize=FALSE,device=if(capabilities("cairo"))grDevices::cairo_pdf else grDevices::pdf)
    provenance<-data.frame(panel=LETTERS[seq_along(page)],key=names(page),title=vapply(page,function(p)as.character(p$labels$title),character(1)))
    pgs_focus_export(c(list(panels=provenance,caption=data.frame(caption=caption)),tables),rd,name)
    manifest[[length(manifest)+1]]<-data.frame(file=paste0(name,".png"),group=stem,panels=length(page),source="PGS focus aggregate results")
  }
  mf<-file.path(rd,"figure_manifest.csv");old<-if(file.exists(mf))as.data.frame(data.table::fread(mf))else data.frame()
  new<-pgs_bind(manifest);if(nrow(old)&&"file"%in%names(old))old<-old[!old$file%in%new$file,]
  data.table::fwrite(pgs_bind(list(old,new)),mf);invisible(pp)
}
pgs_plot_focus<-function(x,outdir,tri=data.frame(),loci=data.frame()) {
  caption<-"Exploratory association decomposition. G denotes the part captured by the supplied PGS; R includes uncaptured genetics, lifestyle, disease and measurement. Opposite associations do not establish antagonistic pleiotropy. Conditional CIs do not include calibration uncertainty; anchor bootstrap refits calibration."
  rd<-file.path(outdir,"c1_correlate")
  pgs_save_panels(pgs_main_panels(x,tri,loci),rd,"c1.FigPGS.overview",x,"Measured levels and genetic prediction",caption)
  chosen<-pgs_pick(x,12);extra<-list()
  a<-x$landmarks
  if(pgs_ok(a,c("feature","model","term","landmark","beta","lo","hi"))) {
    a<-a[a$feature%in%head(chosen,6)&a$model=="joint",];a$series<-paste(a$term,"after",a$landmark,"y")
    extra$landmarks<-pgs_forest_plot(a,"Minimum lead-time sensitivity",unit="Joint log HR / whole-biomarker SD")
  }
  a<-x$bootstrap
  if(pgs_ok(a,c("feature","term","lo","hi","status"))) {
    a<-a[a$term=="difference"&a$status=="ok",];b<-x$contrasts
    if(nrow(a)&&pgs_ok(b,c("feature","beta_difference"))){a$beta<-b$beta_difference[match(a$feature,b$feature)];a$series<-"Refitted bootstrap";
      extra$bootstrap<-pgs_forest_plot(a,"Refitted bootstrap of βG − βR","Participant / family resampling; selected anchors only",unit="Joint coefficient difference")}
  }
  a<-x$lifestyle
  if(pgs_ok(a,c("feature","component","part","beta"))) {
    a<-a[a$feature%in%head(chosen,8)&a$part==".R"&is.finite(a$beta),]
    if(nrow(a))extra$lifestyle<-ggplot2::ggplot(a,ggplot2::aes(component,feature,fill=beta))+ggplot2::geom_tile(color="white")+
      ggplot2::scale_fill_gradient2(low="#B34359",mid="white",high="#287D8E")+
      ggplot2::labs(title="LE8 connections with remaining component R",subtitle="Adjusted contemporaneous associations; not intervention effects",x=NULL,y=NULL,fill=NULL)+pgs_theme()+
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
  }
  a<-x$source_models
  if(pgs_ok(a,c("feature","scope","model","beta","lo","hi"))) {
    a<-a[a$model=="pgs",];a$series<-a$scope
    extra$sources<-pgs_forest_plot(a,"Actual rescoring by genetic source","Common genotype-complete participants; see allele audit",unit="PGS log HR per own SD")
  }
  a<-x$prevalent
  if(pgs_ok(a,c("feature","model","term","beta","lo","hi"))) {
    a<-a[a$feature%in%chosen&a$model=="joint",];a$series<-ifelse(a$term==".G","PGS-captured G","Remaining R")
    extra$prevalent<-pgs_forest_plot(a,"Baseline prevalent cases: G and R","Descriptive case-control analysis; separate from incident hazards",unit="Log OR / whole-biomarker SD")
  }
  pgs_save_panels(extra,rd,"c1.FigPGS.sensitivity",x,"Discordance sensitivity analyses",caption)
}
pgs_publication_joint<-function(trait,analysis_root,layers=c("prot","met")) {
  if(!all(c("prot","met")%in%layers))return(invisible(NULL))
  root<-file.path(analysis_root,trait);objects<-list();tri<-list();loci<-list();tables<-list();pp<-list()
  for(layer in c("prot","met")) {
    f<-file.path(root,layer,"c1_correlate","c1.pgs_focus.rds")
    if(!file.exists(f))return(invisible(NULL));x<-readRDS(f)
    if(!pgs_ok(x$paired,c("feature","evidence_pattern")))return(invisible(NULL))
    cf<-file.path(root,layer,"c3_coloc","c3.coloc_summary.csv")
    co<-if(file.exists(cf))as.data.frame(data.table::fread(cf))else data.frame()
    tr<-make_c3_pgs_integration(coloc_summary=co,focus=x);lo<-attr(tr,"loci")
    parts<-pgs_main_panels(x,tr,lo)
    if(!is.null(parts$paired))pp[[paste0(layer,"_paired")]]<-parts$paired+ggplot2::labs(title=paste(toupper(layer),"| Measured versus PGS"))
    if(!is.null(parts$components))pp[[paste0(layer,"_components")]]<-parts$components+ggplot2::labs(title=paste(toupper(layer),"| Captured and remaining"))
    for(n in names(x))if(is.data.frame(x[[n]]))tables[[paste0(layer,"_",n)]]<-x[[n]]
    tables[[paste0(layer,"_loci")]]<-lo
    pick<-pgs_pick(x,6);cc<-x$calibration
    if(pgs_ok(cc,c("feature","partial_R2"))) {
      cc<-cc[cc$feature%in%pick,];cc$feature<-paste(toupper(layer),cc$feature,sep=": ");objects[[layer]]<-cc
    }
    tr$feature<-paste(toupper(layer),tr$feature,sep=": ");tri[[layer]]<-tr
    if(nrow(lo)){lo$feature<-paste(toupper(layer),lo$feature,sep=": ");loci[[layer]]<-lo}
  }
  cc<-pgs_bind(objects)
  if(nrow(cc)) {
    # Use all selected cross-layer calibration rows, not a new discovery ranking.
    cc<-cc[is.finite(cc$partial_R2),]
    if(nrow(cc))pp$calibration<-ggplot2::ggplot(cc,ggplot2::aes(partial_R2,reorder(feature,partial_R2),color=weak_capture))+
      ggplot2::geom_vline(xintercept=0,color="grey75")+ggplot2::geom_point(size=2.5)+
      ggplot2::scale_color_manual(values=c("FALSE"="#287D8E","TRUE"="#B39557"))+
      ggplot2::labs(title="PGS capture across both omic layers",x="Cross-fitted partial R²",y=NULL,color="Weak capture")+pgs_theme()
  }
  loc<-pgs_locus_panels(pgs_bind(tri),pgs_bind(loci));if(!is.null(loc$loci))pp$loci<-loc$loci
  if(length(pp))pgs_save_panels(pp,root,"Fig1.PGS_integrated",tables,paste(trait,"| Measured and genetically predicted biomarkers"),
    "Exploratory matched-cohort comparisons within each omic layer. G and R share the whole-biomarker scale; R is not pure lifestyle. The protein and metabolite cohorts need not contain the same people. Colocalization remains locus-specific. See each layer's Fig1 for follow-up contrasts and its workbook for refitted bootstrap and source-score sensitivity.")
  invisible(pp)
}
