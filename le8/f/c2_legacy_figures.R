# Restored verbatim plotting functions from le8_review_backup_20260913.
# Kept alongside the new Figure 5/6 comparisons; no estimator changes.
le8_plot_legacy_mr_overview <- function(mr,assoc,layer){
  local_name<-if(layer=="protein")"cis"else"local";distal_name<-if(layer=="protein")"trans"else"distal"
  wide<-mr|>select(exposure,analysis,b,se,pval)|>pivot_wider(names_from=analysis,values_from=c(b,se,pval))|>
    left_join(assoc|>select(exposure=term,obs_beta=beta,obs_se=std.error,obs_p=p.value),by="exposure")
  bl<-paste0("b_",local_name);sl<-paste0("se_",local_name);pl<-paste0("pval_",local_name);bd<-paste0("b_",distal_name);sd<-paste0("se_",distal_name);pd<-paste0("pval_",distal_name)
  a<-wide|>filter(is.finite(.data[[bl]]),is.finite(obs_beta))|>mutate(label=ifelse(min_rank(pmin(obs_p,.data[[pl]],na.rm=TRUE))<=12,exposure,NA_character_))
  pA<-if(!nrow(a))blank_plot(paste0("a. Observational versus ",local_name," MR"),"No paired estimate") else ggplot(a,aes(obs_beta,.data[[bl]]))+geom_hline(yintercept=0,color="grey70")+geom_vline(xintercept=0,color="grey70")+
    geom_smooth(method="lm",se=TRUE,color="grey35",fill="grey85",linewidth=.7)+geom_point(aes(color=.data[[pl]]<.05),size=2)+ggrepel::geom_text_repel(aes(label=label),size=2.8,fontface="bold",seed=1,max.overlaps=20,na.rm=TRUE)+
    scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="grey70"),guide="none")+labs(title=paste0("a. Observational versus ",local_name,"-MR effects"),x="Observational effect",y=paste0(local_name," MR effect"))+theme_5c(11)
  b<-wide|>filter(is.finite(.data[[bl]]),is.finite(.data[[bd]]))|>mutate(label=ifelse(min_rank(pmin(.data[[pl]],.data[[pd]],na.rm=TRUE))<=12,exposure,NA_character_))
  pB<-if(!nrow(b))blank_plot(paste0("b. ",local_name," versus ",distal_name),"No trait had both instrument classes") else ggplot(b,aes(.data[[bl]],.data[[bd]]))+geom_abline(slope=1,intercept=0,linetype=2,color="grey45")+
    geom_hline(yintercept=0,color="grey70")+geom_vline(xintercept=0,color="grey70")+geom_smooth(method="lm",se=TRUE,color="grey35",fill="grey85",linewidth=.7)+
    geom_point(aes(color=sign(.data[[bl]])==sign(.data[[bd]])),size=2)+ggrepel::geom_text_repel(aes(label=label),size=2.8,fontface="bold",seed=2,max.overlaps=20,na.rm=TRUE)+
    scale_color_manual(values=c(`TRUE`="#1B9E77",`FALSE`="#D95F02"),guide="none")+labs(title=paste0("b. ",local_name," versus ",distal_name," MR"),x=paste0(local_name," effect"),y=paste0(distal_name," effect"))+theme_5c(11)
  ev<-bind_rows(
    assoc|>transmute(exposure=term,evidence="Observational",effect=beta,p=p.value),
    mr|>filter(analysis%in%c(local_name,distal_name))|>
      transmute(exposure,evidence=ifelse(analysis==local_name,paste0("MR: ",local_name),paste0("MR: ",distal_name)),effect=b,p=pval))|>
    filter(is.finite(p),is.finite(effect))|>group_by(evidence)|>mutate(q=p.adjust(p,"BH"))|>ungroup()
  top_ev<-ev|>group_by(exposure)|>summarise(best_p=min(p),.groups="drop")|>slice_min(best_p,n=28,with_ties=FALSE)|>pull(exposure)
  ev<-ev|>filter(exposure%in%top_ev)|>mutate(score=stable_neglog10_p(p),
      signed_score=sign(effect)*pmin(score,12),significant=q<.05,
      evidence=factor(evidence,levels=c("Observational",paste0("MR: ",local_name),paste0("MR: ",distal_name))),
      exposure=factor(exposure,levels=rev(top_ev)))
  pC<-if(!nrow(ev))blank_plot("c. Evidence matrix")else ggplot(ev,aes(evidence,exposure))+
    geom_tile(aes(fill=signed_score),color="white",linewidth=.35)+
    geom_point(data=ev|>filter(significant),shape=8,size=2.1,color="black")+
    scale_fill_gradient2(low="#3F78A8",mid="white",high="#C86B4A",midpoint=0,limits=c(-12,12),
      name="signed\n-log10(P)")+
    labs(title="c. Signed cross-evidence matrix",
      subtitle="Blue = protective; orange = risk-increasing; asterisk = within-analysis FDR < 0.05; scale capped at 12",
      x=NULL,y=NULL)+theme_5c(8)+theme(legend.position="right")
  list(plot=(pA|pB)/pC+plot_layout(heights=c(1,.9)),wide=wide)
}

le8_plot_legacy_qtl_variance <- function(mr,layer){
  cls<-if(layer=="protein")c("cis","trans")else c("local","distal")
  ps<-map(cls,function(cc){
    d<-mr|>filter(analysis==cc,is.finite(r2_median),n_IV>0)|>
      mutate(r2_med=100*r2_median,r2_lo=100*r2_q25,r2_hi=100*r2_q75,r2_p90_plot=100*r2_p90)|>
      slice_max(r2_p90_plot,n=20,with_ties=FALSE)|>arrange(r2_med)|>mutate(exposure=factor(exposure,levels=exposure))
    if(!nrow(d))return(blank_plot(paste0(toupper(substr(cc,1,1)),substr(cc,2,nchar(cc))," instruments"),"No valid instrument set"))
    ggplot(d,aes(r2_med,exposure))+
      geom_segment(aes(x=r2_med,xend=r2_p90_plot,yend=exposure),color="grey72",linewidth=.7)+
      geom_errorbarh(aes(xmin=r2_lo,xmax=r2_hi),height=.16,color=ifelse(cc%in%c("cis","local"),"#D95F02","#12AEB5"),linewidth=.8)+
      geom_point(aes(size=pmin(n_IV,500)),color=ifelse(cc%in%c("cis","local"),"#D95F02","#12AEB5"))+
      geom_text(aes(x=r2_p90_plot,label=sprintf("median %.3f%%; K=%d",r2_med,n_IV)),hjust=-.05,size=2.35,fontface="bold")+
      scale_x_continuous(expand=expansion(mult=c(.02,.34)))+scale_size_continuous(range=c(1.6,4.8),name="IV count")+
      labs(title=paste0(cc," QTL instruments"),subtitle="Per-IV partial R²: point = median; interval = IQR; grey tail = 90th percentile",
           x="Per-instrument partial phenotypic R² (%)",y=NULL)+theme_5c(8)
  })
  (ps[[1]]|ps[[2]])+
    plot_annotation(title="QTL instrument-level variance architecture",
      subtitle="This figure does not sum R² and does not use 1 - product(1 - R²); values cannot saturate merely because K is large",
      caption="Instrument-level R² distribution",
      theme=theme(plot.title=element_text(face="bold",size=15),
                  plot.subtitle=element_text(size=10,color="grey30"),
                  plot.caption=element_text(size=8,color="grey40",hjust=1)))
}

