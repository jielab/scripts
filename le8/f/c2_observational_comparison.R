# Fig6-style triangulation: independent cis/local MR, incident Cox and prevalent
# logistic results. This adaptation uses within-analysis FDR, not the paper's
# high-confidence MR + coloc gate. No choice of minimum P across cis/trans.
le8_c2_observational_data<-function(mr,incident,prevalent,layer) {
  primary<-if(layer=='protein')'cis'else'local'
  m<-mr|>filter(analysis==primary)
  if(anyDuplicated(m$exposure))stop('Duplicate primary MR estimates')
  if(!'FDR_analysis'%in%names(m))m$FDR_analysis<-p.adjust(m$pval,'BH')
  m<-m|>transmute(exposure,b_mr=b,se_mr=se,p_mr=pval,q_mr=FDR_analysis)
  obs<-function(d,suffix) {
    if(anyDuplicated(d$term))stop('Duplicate observational estimates')
    if(!'FDR'%in%names(d))d$FDR<-p.adjust(d$p.value,'BH')
    z<-d|>transmute(exposure=term,b=beta,se=std.error,p=p.value,q=FDR)
    names(z)[-1]<-paste0(names(z)[-1],'_',suffix);z
  }
  z<-full_join(m,obs(incident,'incident'),by='exposure')|>full_join(obs(prevalent,'prevalent'),by='exposure')
  for(k in c('mr','incident','prevalent')) {
    z[[paste0('available_',k)]]<-is.finite(z[[paste0('b_',k)]])&is.finite(z[[paste0('se_',k)]])&
      z[[paste0('se_',k)]]>0&is.finite(z[[paste0('q_',k)]])
    z[[paste0('supported_',k)]]<-z[[paste0('available_',k)]]&z[[paste0('q_',k)]]<.05
  }
  z
}

le8_c2_support_counts<-function(z) {
  levels<-c('Concordant','Discordant','No FDR support','Unavailable')
  rows<-list()
  for(endpoint in c('incident','prevalent'))for(direction in c('MR to observational','Observational to MR')) {
    source<-if(direction=='MR to observational')'mr'else endpoint
    target<-if(direction=='MR to observational')endpoint else'mr'
    d<-z[z[[paste0('supported_',source)]],,drop=FALSE]
    state<-ifelse(!d[[paste0('available_',target)]],'Unavailable',
      ifelse(!d[[paste0('supported_',target)]],'No FDR support',
        ifelse(sign(d[[paste0('b_',source)]])==sign(d[[paste0('b_',target)]]),'Concordant','Discordant')))
    counts<-table(factor(state,levels=levels))
    rows[[length(rows)+1L]]<-tibble(endpoint=endpoint,direction=direction,status=levels,
      n=as.integer(counts),denominator=nrow(d),fraction=if(nrow(d))as.integer(counts)/nrow(d)else NA_real_)
  }
  bind_rows(rows)
}

le8_plot_mr_incident_prevalent<-function(mr,incident,prevalent,layer) {
  z<-le8_c2_observational_data(mr,incident,prevalent,layer)
  counts<-le8_c2_support_counts(z)
  if(!any(z$available_mr)||!any(z$available_incident)||!any(z$available_prevalent))
    return(list(plot=blank_plot('MR, incident and prevalent comparison','One or more result sources unavailable'),data=z,counts=counts))
  pal<-c('Concordant'='#E37768','Discordant'='#587DAB','No FDR support'='#BDC5C4','Unavailable'='#ECECEC')
  bars<-lapply(c('incident','prevalent'),function(ep) {
    d<-counts|>filter(endpoint==ep)|>mutate(status=factor(status,levels=names(pal)),
      row=paste0(ifelse(direction=='MR to observational','MR → observed','Observed → MR'),'  (n=',denominator,')'))
    ggplot(d,aes(fraction,row,fill=status))+geom_col(width=.56,position=position_stack(reverse=TRUE),color='white',linewidth=.3)+
      geom_text(aes(label=ifelse(is.finite(fraction)&fraction>=.065,n,'')),
        position=position_stack(vjust=.5,reverse=TRUE),size=3.5,color='grey20')+
      scale_fill_manual(values=pal,drop=FALSE)+scale_x_continuous(labels=scales::label_percent(),limits=c(0,1),expand=expansion(mult=c(0,.01)))+
      labs(title=paste0(if(ep=='incident')'A. Incident'else'B. Prevalent',' evidence convergence'),
        subtitle='Source set: FDR < 0.05; opposite evidence tested in the same biomarkers',x='Fraction of source discoveries',y=NULL,fill=NULL)+theme_5c(11)
  })
  anchors<-if(layer=='protein')c('PCSK9','LPA','GDF15','NTPROBNP','MMP12')else
    c('L_VLDL_TG.pct','ApoB','Glucose','GlycA','Total_TG')
  eligible<-z|>filter(available_incident|available_prevalent)|>arrange(q_mr,p_mr,exposure)
  observed<-eligible|>arrange(p_incident)|>slice_head(n=3)|>pull(exposure)
  chosen<-head(unique(c(intersect(anchors,eligible$exposure),observed,eligible$exposure)),16)
  z$shown_forest<-z$exposure%in%chosen
  forests<-lapply(c('mr','incident','prevalent'),function(k) {
    d<-z|>filter(shown_forest)|>transmute(exposure=factor(exposure,levels=rev(chosen)),
      b=.data[[paste0('b_',k)]],se=.data[[paste0('se_',k)]],available=.data[[paste0('available_',k)]],
      supported=.data[[paste0('supported_',k)]],b_mr,supported_mr)
    d<-d|>mutate(status=case_when(!available~'Unavailable',!supported~'No FDR support',
      k=='mr'~'MR FDR support',!supported_mr~'Observed FDR support',
      sign(b)==sign(b_mr)~'Concordant',TRUE~'Discordant'),lo=b-1.96*se,hi=b+1.96*se)
    title<-switch(k,mr=paste0('C. ',if(layer=='protein')'cis'else'local',' Mendelian randomization'),
      incident='D. Incident disease · Cox',prevalent='E. Prevalent disease · logistic')
    ggplot(d,aes(y=exposure))+geom_hline(yintercept=seq_along(chosen),color='grey94',linewidth=.4)+
      geom_vline(xintercept=0,color='grey65',linetype=2,linewidth=.4)+
      geom_segment(data=filter(d,available),aes(x=lo,xend=hi,yend=exposure,color=status),linewidth=.8)+
      geom_point(data=filter(d,available),aes(x=b,color=status,shape=supported),size=2.7)+
      geom_point(data=filter(d,!available),aes(x=0),shape=4,color='grey65',size=2)+
      scale_y_discrete(drop=FALSE)+scale_color_manual(values=c(pal,'MR FDR support'='#866AA3','Observed FDR support'='#B38D43'))+
      scale_shape_manual(values=c('TRUE'=16,'FALSE'=1))+guides(color='none',shape='none')+
      labs(title=title,subtitle=if(k=='mr')'Genetic effect; 95% CI'else'Measured baseline omic level; 95% CI',
        x=if(k=='incident')'Log hazard ratio'else'Log odds ratio',y=NULL)+theme_5c(11)+theme(panel.grid=element_blank())
  })
  caption<-paste('Fixed anchors + top incident signals + primary MR ranking; the same rows are retained in all three forests.',
    'Filled points: within-analysis FDR < 0.05; open points: no FDR support; ×: unavailable.',
    'Purple: MR support; gold: observed support without MR support; coral/blue: concordant/discordant supported effects.',
    'Cox HRs and logistic/MR ORs are different estimands; separate horizontal scales are used.',
    'MR significance alone is not causal validation. Unlike the reference Figure 6, this descriptive adaptation does not impose a colocalization gate.')
  list(plot=wrap_plots(c(bars,forests),design='AAABBB\nCCDDEE',heights=c(.55,1.5))+
    plot_annotation(caption=caption),data=z,counts=counts)
}

le8_emit_restored_c2<-function(mr,assoc,layer,outdir) {
  rd<-le8_job_dir(outdir,'c2_cause')
  old<-le8_plot_legacy_mr_overview(mr,assoc,layer)
  save_plot(old$plot,'c2.Fig20.observational_mr_overview.png',20,15,outdir=outdir)
  save_plot(le8_plot_legacy_qtl_variance(mr,layer),'c2.Fig21.qtl_variance_ranked.png',19,12,outdir=outdir)
  prefix<-if(layer=='protein')'pwas'else'mwas'
  files<-file.path(outdir,'c1_correlate',paste0(prefix,c('_incident_adj2.csv','_prevalent_adj2.csv')))
  if(all(file.exists(files))) {
    z<-le8_plot_mr_incident_prevalent(mr,as_tibble(data.table::fread(files[1])),as_tibble(data.table::fread(files[2])),layer)
    save_plot(z$plot,'c2.Fig22.mr_incident_prevalent.png',22,14,outdir=outdir)
    write_raw_csv(z$data,'c2.mr_incident_prevalent.csv',rd)
    write_raw_csv(z$counts,'c2.mr_incident_prevalent_support.csv',rd)
  } else save_plot(blank_plot('MR, incident and prevalent comparison','C1 incident/prevalent tables unavailable'),
    'c2.Fig22.mr_incident_prevalent.png',22,14,outdir=outdir)
}
