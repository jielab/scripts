for(.pgs_file in c("pgs_core.R","pgs_focus_analysis.R","pgs_figures.R"))
  source(file.path(Sys.getenv("LE8_FDIR",unset="."),.pgs_file))
# Locus-resolved triangulation. A genome-wide PGS is never an input to coloc.
# Keep every tested region; overlapping regions are not independent discoveries.
pgs_coloc_loci <- function(x,h4=.70) {
  d<-as.data.frame(x);if(!nrow(d)||!all(c("feature","PP.H4")%in%names(d)))return(data.frame())
  numeric_cols<-c("chr","start","end","PP.H3","PP.H4_robust_min","credible_set_n")
  for(n in numeric_cols)if(!n%in%names(d))d[[n]]<-NA_real_
  for(n in c("locus","locus_class","status","lead_shared"))if(!n%in%names(d))d[[n]]<-NA_character_
  d$chr<-sub("^chr","",as.character(d$chr),ignore.case=TRUE)
  d$prior_tested<-is.finite(d$PP.H4_robust_min)
  d$robust_coloc<-d$status%in%"ok"&d$prior_tested&d$PP.H4_robust_min>=h4
  d$locus_evidence<-ifelse(!d$status%in%"ok","Unavailable / failed",
    ifelse(d$robust_coloc,"H4 robust across tested priors",ifelse(d$PP.H4>=h4&!d$prior_tested,
      "H4 default; prior sensitivity unavailable",ifelse(d$PP.H4>=h4,"H4 sensitive to priors",
      ifelse(is.finite(d$PP.H3)&d$PP.H3>=h4,"H3: distinct signals favored","No strong H3/H4 support")))))
  d$region_cluster<-NA_character_
  for(ch in unique(d$chr[!is.na(d$chr)])) {
    ix<-which(d$chr==ch&is.finite(d$start)&is.finite(d$end)&d$start<=d$end)
    ix<-ix[order(d$start[ix],d$end[ix])];right<--Inf;k<-0L
    for(i in ix){if(d$start[i]>right){k<-k+1L;right<-d$end[i]}else right<-max(right,d$end[i]);d$region_cluster[i]<-paste0("chr",ch,"_region",k)}
  }
  # Unknown coordinates are not merged based only on a reused lead label.
  missing<-which(is.na(d$region_cluster));d$region_cluster[missing]<-paste0("unmapped_",missing)
  d$cluster_features<-vapply(d$region_cluster,function(k)length(unique(d$feature[d$region_cluster==k])),integer(1))
  d$cluster_robust_features<-vapply(d$region_cluster,function(k)length(unique(d$feature[d$region_cluster==k&d$robust_coloc])),integer(1))
  d$interpretation<-"Overlapping windows define a region cluster, not one fine-mapped causal signal; H4 is locus-specific"
  d
}
make_c3_pgs_integration <- function(c1=list(),c2=list(),coloc_summary=data.frame(),h4=.70,focus=NULL) {
  loci<-pgs_coloc_loci(coloc_summary,h4)
  paired<-if(is.list(focus))focus$paired else NULL
  if(is.data.frame(paired)&&nrow(paired)) {
    z<-as.data.frame(paired);names(z)[names(z)=="measured_beta"]<-"observed_beta"
    names(z)[names(z)=="measured_p"]<-"observed_p";names(z)[names(z)=="measured_FDR"]<-"observed_FDR"
    z$matched_comparison<-TRUE
  }else {
    a<-c1$association_adj2;if(is.null(a))a<-c1$association
    p<-c1$pgs_incident
    conv<-function(d,pgs=FALSE){if(is.null(d)||!nrow(d))return(data.frame(feature=character()))
      pre<-if(pgs)"pgs"else"observed";o<-data.frame(feature=d$term)
      o[[paste0(pre,"_beta")]]<-d$beta;o[[paste0(pre,"_p")]]<-d$p.value;o[[paste0(pre,"_FDR")]]<-d$FDR;o}
    z<-merge(conv(a),conv(p,TRUE),by="feature",all=TRUE);z$matched_comparison<-rep(FALSE,nrow(z))
  }
  features<-unique(c(z$feature,loci$feature))
  if(!length(features)){z<-data.frame(feature=character());attr(z,"loci")<-loci;return(z)}
  z<-merge(data.frame(feature=features),z,by="feature",all.x=TRUE)
  if(!nrow(z)){attr(z,"loci")<-loci;return(z)}
  for(n in c("observed_beta","observed_p","observed_FDR","pgs_beta","pgs_p","pgs_FDR"))if(!n%in%names(z))z[[n]]<-NA_real_
  if(!"matched_comparison"%in%names(z))z$matched_comparison<-FALSE
  z$matched_comparison[is.na(z$matched_comparison)]<-FALSE
  z$observed_supported<-is.finite(z$observed_FDR)&z$observed_FDR<.05
  z$pgs_supported<-is.finite(z$pgs_FDR)&z$pgs_FDR<.05
  z$pgs_observed_sign_match<-is.finite(z$observed_beta)&is.finite(z$pgs_beta)&z$observed_beta*z$pgs_beta>=0
  z$triangulation_class<-pgs_classify(z$observed_beta,z$pgs_beta,z$observed_FDR,z$pgs_FDR)
  z$n_loci_tested<-z$n_robust_loci<-z$n_region_clusters<-0L;z$robust_PP4<-NA_real_;z$coloc_supported<-FALSE
  z$coloc_summary<-"Not tested / unavailable"
  for(i in seq_len(nrow(z))) {
    a<-loci[loci$feature==z$feature[i],,drop=FALSE];if(!nrow(a))next
    z$n_loci_tested[i]<-sum(a$status%in%"ok");z$n_robust_loci[i]<-sum(a$robust_coloc)
    z$n_region_clusters[i]<-length(unique(a$region_cluster[a$robust_coloc]))
    z$robust_PP4[i]<-if(any(is.finite(a$PP.H4_robust_min)))max(a$PP.H4_robust_min,na.rm=TRUE)else NA_real_
    z$coloc_supported[i]<-any(a$robust_coloc)
    z$coloc_summary[i]<-paste(unique(paste(a$locus_class,a$locus_evidence,sep=": ")),collapse="; ")
  }
  z$inherited_locus_pattern<-z$observed_supported&z$pgs_supported&z$pgs_observed_sign_match&z$coloc_supported
  z$reactive_compatible_pattern<-FALSE # Absence of genetic support is not evidence of consequence.
  z$interpretation<-ifelse(z$matched_comparison,
    "Identical-sample associations; opposite directions do not establish antagonistic pleiotropy or causal partition",
    "Legacy different-sample screening only; run pgs_focus before comparing direction or statistical strength")
  if(is.list(focus)&&is.data.frame(focus$calibration)&&nrow(focus$calibration))z<-merge(z,focus$calibration,by="feature",all.x=TRUE)
  z<-z[order(!(z$triangulation_class=="Both supported: opposite"),z$pgs_p,z$observed_p),,drop=FALSE]
  attr(z,"loci")<-loci;z
}
read_c3_pgs_integration <- function(layer,outdir,coloc_summary) {
  get<-function(module,file){p<-file.path(outdir,module,file);if(file.exists(p))tryCatch(readRDS(p),error=function(e)list())else list()}
  x<-make_c3_pgs_integration(get("c1_correlate","c1.res.rds"),get("c2_cause","c2.res.rds"),
    coloc_summary,get0("H4_STRONG",ifnotfound=.70),get("c1_correlate","c1.pgs_focus.rds"))
  rd<-file.path(outdir,"c3_coloc");dir.create(rd,recursive=TRUE,showWarnings=FALSE)
  loci<-attr(x,"loci");clusters<-if(nrow(loci))unique(loci[,c("region_cluster","chr","cluster_features","cluster_robust_features")])else data.frame()
  pgs_focus_export(list(features=x,loci=loci,region_clusters=clusters),rd,"c3.pgs_focus")
  # Refresh the established CSV name as well; downstream readers must not see
  # an old "PGS weak" label after a focused rerun.
  data.table::fwrite(as.data.frame(x),file.path(rd,"c3.pgs_observed_coloc_triangulation.csv"))
  x
}
plot_c3_pgs_integration <- function(x,outdir) {
  if(!nrow(x))return(invisible(NULL))
  rd<-file.path(outdir,"c3_coloc");dir.create(rd,recursive=TRUE,showWarnings=FALSE)
  pp<-pgs_locus_panels(x,attr(x,"loci"))
  if(length(pp))pgs_save_panels(pp,rd,"c3.Fig7.pgs_coloc_triangulation",
    list(features=as.data.frame(x),loci=attr(x,"loci")),"Measured / PGS discordance and locus evidence")
  invisible(x)
}
