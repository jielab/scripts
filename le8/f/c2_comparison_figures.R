# Effect/CI comparison inspired by Koprulu et al., Cell 2026, Figure 5.
# Estimates are this run's MR results; the display does not re-estimate MR.
le8_c2_paired_estimates<-function(mr,layer) {
  classes<-if(layer=='protein')c('cis','trans')else c('local','distal')
  d<-mr|>filter(analysis%in%classes)
  if(anyDuplicated(d[,c('exposure','analysis')]))stop('Multiple primary MR rows for one exposure/class; specify the primary estimator before plotting')
  for(nm in c('n_IV','FDR_analysis','LD_status','instrument_snps'))if(!nm%in%names(d))d[[nm]]<-NA
  d<-d|>group_by(analysis)|>mutate(FDR_analysis=coalesce(FDR_analysis,p.adjust(pval,'BH')))|>ungroup()
  wide<-d|>select(exposure,analysis,b,se,pval,n_IV,FDR_analysis,LD_status,instrument_snps)|>
    pivot_wider(names_from=analysis,values_from=-c(exposure,analysis))
  for(field in c('b','se','pval','n_IV','FDR_analysis','LD_status','instrument_snps')) {
    for(i in 1:2) {
      nm<-paste0(field,'_',classes[i]);if(!nm%in%names(wide))wide[[nm]]<-NA
      wide[[paste0(field,if(i==1)'_primary'else'_distal')]]<-wide[[nm]]
    }
  }
  wide<-wide|>mutate(paired=is.finite(b_primary)&is.finite(b_distal)&is.finite(se_primary)&
      is.finite(se_distal)&se_primary>0&se_distal>0,
    effect_difference=b_primary-b_distal,difference_se=sqrt(se_primary^2+se_distal^2),
    Q_difference=ifelse(paired,(effect_difference/difference_se)^2,NA_real_),
    I2_difference=ifelse(paired,ifelse(Q_difference>0,pmax(0,(Q_difference-1)/Q_difference)*100,0),NA_real_),
    p_heterogeneity=ifelse(paired,2*pnorm(-abs(effect_difference/difference_se)),NA_real_),
    FDR_heterogeneity=p.adjust(p_heterogeneity,'BH'),
    heterogeneity=case_when(!paired~'Unavailable',p_heterogeneity>=.01~'No detected difference',
      sign(b_primary)!=sign(b_distal)~'Different effects; opposite signs',TRUE~'Different effects; same sign'),
    primary_supported=FDR_analysis_primary<.05,
    lo_primary=b_primary-1.96*se_primary,hi_primary=b_primary+1.96*se_primary,
    lo_distal=b_distal-1.96*se_distal,hi_distal=b_distal+1.96*se_distal)
  wide
}

le8_plot_cis_trans_comparison<-function(mr,assoc,layer) {
  cls<-if(layer=='protein')c('cis','trans')else c('local','distal')
  wide<-le8_c2_paired_estimates(mr,layer)
  wide<-wide|>left_join(assoc|>select(exposure=term,obs_beta=beta,obs_se=std.error,obs_p=p.value),by='exposure')
  focus<-wide|>filter(paired,primary_supported)
  selection<-paste0(cls[1],' MR FDR < 0.05 with an estimable ',cls[2],' effect')
  if(!nrow(focus)) {focus<-wide|>filter(paired);selection<-'All estimable pairs; no primary-class FDR discovery'}
  if(!nrow(focus))return(list(plot=blank_plot('Paired genetic effect comparison','No biomarker has both instrument classes'),wide=wide))
  anchors<-if(layer=='protein')c('PCSK9','LPA','GDF15','NTPROBNP','MMP12')else
    c('L_VLDL_TG.pct','L_VLDL_TG','Total_TG','ApoB','Glucose')
  labels<-head(unique(c(intersect(anchors,focus$exposure),focus$exposure[order(focus$pval_primary)])),10)
  focus$label<-ifelse(focus$exposure%in%labels,focus$exposure,NA_character_)
  pal<-c('Different effects; opposite signs'='#BB4560','Different effects; same sign'='#DA953F',
    'No detected difference'='#3B7D9B')
  pa<-ggplot(focus,aes(b_primary,b_distal))+
    geom_abline(slope=1,intercept=0,linetype=2,color='grey55',linewidth=.5)+
    geom_hline(yintercept=0,color='grey85',linewidth=.4)+geom_vline(xintercept=0,color='grey85',linewidth=.4)+
    geom_segment(aes(x=lo_primary,xend=hi_primary,yend=b_distal,color=heterogeneity),alpha=.32,linewidth=.4)+
    geom_segment(aes(y=lo_distal,yend=hi_distal,xend=b_primary,color=heterogeneity),alpha=.32,linewidth=.4)+
    geom_point(aes(color=heterogeneity),size=2.8,alpha=.95)+
    ggrepel::geom_text_repel(aes(label=label),size=3.1,seed=SEED,max.overlaps=Inf,na.rm=TRUE,
      box.padding=.5,min.segment.length=0,color='grey20')+
    scale_color_manual(values=pal,drop=FALSE)+labs(title=paste0('A. ',cls[1],' and ',cls[2],' effects on ',Y),
      subtitle=paste0(nrow(focus),' biomarkers | ',selection),
      x=paste0(cls[1],' MR effect'),y=paste0(cls[2],' MR effect'),color=NULL)+theme_5c(11)+
    theme(panel.grid=element_blank())
  # Fixed anchors first, followed by the strongest primary-class signals.
  candidates<-wide|>filter(paired)|>arrange(pval_primary,exposure)
  examples<-head(unique(c(intersect(anchors,candidates$exposure),candidates$exposure)),12)
  forest<-bind_rows(lapply(1:2,function(i) {
    suffix<-if(i==1)'primary'else'distal'
    candidates|>filter(exposure%in%examples)|>transmute(exposure,class=cls[i],
      b=.data[[paste0('b_',suffix)]],se=.data[[paste0('se_',suffix)]],n_IV=.data[[paste0('n_IV_',suffix)]],
      pval=.data[[paste0('pval_',suffix)]],
      y=match(exposure,rev(examples))+if(i==1).16 else -.16)
  }))|>mutate(lo=b-1.96*se,hi=b+1.96*se)
  row_labels<-candidates[match(rev(examples),candidates$exposure),]|>
    transmute(label=paste0(exposure,'  [',n_IV_primary,' / ',n_IV_distal,']'))|>pull(label)
  pb<-ggplot(forest,aes(b,y,color=class))+
    geom_hline(yintercept=seq(.5,length(examples)+.5,by=1),color='grey94',linewidth=.35)+
    geom_vline(xintercept=0,color='grey65',linetype=2,linewidth=.5)+
    geom_segment(aes(x=lo,xend=hi,yend=y),linewidth=.8)+geom_point(size=2.6)+
    scale_color_manual(values=setNames(c('#D88D2D','#337EAB'),cls),breaks=cls)+
    scale_y_continuous(breaks=seq_along(examples),labels=row_labels,expand=expansion(add=.5))+
    labs(title='B. Paired estimates and 95% confidence intervals',
      subtitle=paste0('Fixed anchors + strongest ',cls[1],' signals; [',cls[1],' / ',cls[2],' IV counts]'),
      x='MR effect (disease GWAS scale)',y=NULL,color='Instrument class')+theme_5c(11)+theme(panel.grid=element_blank())
  wide$shown_scatter<-wide$exposure%in%focus$exposure;wide$shown_forest<-wide$exposure%in%examples
  wide$heterogeneity_assumption<-'Approximate independent-estimate comparison; covariance / cross-class LD not modelled'
  caption<-paste0('Points and intervals are the saved primary MR estimates and 95% CIs. Dashed diagonal: equal effects. ',
    'Colour uses exploratory P(diff) < 0.01, with Var(diff) = SE1^2 + SE2^2; cross-class covariance is unavailable. ',
    'Q (1 df), I² and BH-adjusted P values are in c2.cis_trans_comparison.csv. ',
    sum(focus$n_IV_primary==1&focus$n_IV_distal==1),' / ',nrow(focus),' scatter pairs use one IV in both classes. ',
    'The current instrument set has no protein-specific / phenome pleiotropy filter. ',
    'Primary-class MR significance alone does not establish colocalization or causal validity.')
  list(plot=(pa|pb)+plot_layout(widths=c(1.05,1))+plot_annotation(caption=caption),wide=wide,forest=forest)
}

le8_plot_instrument_comparison<-function(mr,layer) {
  cls<-if(layer=='protein')c('cis','trans')else c('local','distal')
  d<-mr|>filter(analysis%in%cls)
  paired<-d|>select(exposure,analysis,r2_median)|>pivot_wider(names_from=analysis,values_from=r2_median)
  if(!all(cls%in%names(paired)))return(blank_plot('Instrument architecture','One instrument class is unavailable'))
  paired<-paired|>filter(is.finite(.data[[cls[1]]]),is.finite(.data[[cls[2]]]),.data[[cls[1]]]>0,.data[[cls[2]]]>0)
  anchors<-if(layer=='protein')c('PCSK9','LPA','GDF15','NTPROBNP','MMP12')else c('ApoB','Glucose','Total_TG','L_VLDL_TG.pct')
  paired$label<-ifelse(paired$exposure%in%anchors,paired$exposure,NA_character_)
  pa<-ggplot(paired,aes(100*.data[[cls[1]]],100*.data[[cls[2]]]))+
    geom_abline(slope=1,intercept=0,linetype=2,color='grey65')+geom_point(color='#397F9D',alpha=.55,size=2)+
    ggrepel::geom_text_repel(aes(label=label),size=3,seed=SEED,na.rm=TRUE,max.overlaps=Inf)+
    scale_x_log10()+scale_y_log10()+labs(title='A. Instrument strength across the same biomarkers',
      subtitle=paste0(nrow(paired),' paired sets; median per-IV variance, without summing across SNPs'),
      x=paste0(cls[1],' median per-IV partial R² (%)'),y=paste0(cls[2],' median per-IV partial R² (%)'))+
    theme_5c(11)+theme(panel.grid=element_blank())
  d<-d|>mutate(status=case_when(!is.finite(b)|!is.finite(pval)~'No estimate',
    grepl('fallback',LD_status,fixed=TRUE)~'Single IV: LD unavailable',n_IV==1~'Single eligible IV',TRUE~'Multiple retained IVs'))
  counts<-d|>count(analysis,status)|>mutate(analysis=factor(analysis,levels=rev(cls)))
  pb<-ggplot(counts,aes(n,analysis,fill=status))+geom_col(width=.55)+
    geom_text(aes(label=ifelse(n>=10,n,'')),position=position_stack(vjust=.5),size=3.5,color='white')+
    scale_fill_manual(values=c('Single IV: LD unavailable'='#BC6575','Single eligible IV'='#6897AD',
      'Multiple retained IVs'='#44836F','No estimate'='#B7BDC4'))+
    labs(title='B. What the saved MR models could use',subtitle='Retained instrument sets after harmonization and LD handling',
      x='Biomarker–instrument-class records',y=NULL,fill=NULL)+theme_5c(11)+theme(panel.grid=element_blank())
  (pa|pb)+plot_annotation(caption='Per-IV partial R² is not total genetic variance explained. LD-unavailable fallback is shown explicitly; it is not evidence that the biological architecture is monogenic.')
}
