make_plot<-function(d,kind,metric,title='GRID'){
  if(!nrow(d))return(empty_plot())
  if(nrow(d)>1000)return(empty_plot('More than 1,000 rows: refine the table filters before plotting.'))
  theme<-theme_minimal(base_size=13)+theme(legend.position='bottom',plot.title=element_text(size=14,face='bold'),panel.grid.minor=element_blank())
  method<-first_col(d,c('method','model'));if(is.null(method)){d$.method<-'Result';method<-'.method'}
  if(kind=='abm_selector_curve'){
    z<-reshape(d[,c('purpose','selector','actual_coverage','candidate_gain_same_subset','selector_delta_loss_vs_control')],
      varying=c('candidate_gain_same_subset','selector_delta_loss_vs_control'),v.names='value',timevar='question',
      times=c('candidate gain vs fallback','selector effect vs matched control'),direction='long')
    return(ggplot(z,aes(actual_coverage,value,colour=question,linetype=purpose))+geom_hline(yintercept=0,linetype=3)+geom_line()+geom_point()+
      labs(x='Actual coverage',y='Primary-loss difference',title=title)+theme)
  }
  if(kind=='abm_primary_contrasts'){
    z<-reshape(d[,intersect(c('purpose','selector','candidate_minus_fallback_loss','selector_minus_control_loss'),names(d)),drop=FALSE],
      varying=c('candidate_minus_fallback_loss','selector_minus_control_loss'),v.names='value',timevar='contrast',
      times=c('candidate - fallback','selector - matched control'),direction='long')
    return(ggplot(z,aes(value,purpose,shape=contrast))+geom_vline(xintercept=0,linetype=2)+geom_point(size=3)+labs(x='Primary-loss difference',y=NULL,title=title)+theme)
  }
  if(kind=='abm_reference_cleaning'){
    return(ggplot(d,aes(delta_loss,contrast))+geom_vline(xintercept=0,linetype=2)+geom_point(size=3)+
      geom_segment(aes(x=lower,xend=upper,yend=contrast),na.rm=TRUE)+labs(x='Paired loss difference (negative favours first model)',y=NULL,title=title)+theme)
  }
  if(kind%in%c('contrasts','legacy_contrasts')){
    value<-first_col(d,c('delta','difference'));lo<-first_col(d,c('lower','lower95'));hi<-first_col(d,c('upper','upper95'))
    if(is.null(value))return(empty_plot())
    labels<-intersect(c('contrast','method','reference','subset'),names(d));d$.label<-apply(d[,labels,drop=FALSE],1,paste,collapse=' | ')
    p<-ggplot(d,aes(x=.data[[value]],y=.label))+geom_point(size=2)+geom_vline(xintercept=0,linetype=2)
    if(!is.null(lo)&&!is.null(hi))p<-p+geom_segment(aes(x=.data[[lo]],xend=.data[[hi]],yend=.label))
    return(p+labs(x=paste('Paired difference —',if('metric'%in%names(d))paste(unique(d$metric),collapse='; ')else metric),y=NULL,title=title)+theme)
  }
  if(!metric%in%names(d))return(empty_plot('Choose a numeric metric for this table.'))
  if(kind=='confidence')return(ggplot(d,aes(x=mean_confidence,y=.data[[metric]],colour=.data[[method]],group=.data[[method]]))+geom_point()+geom_line()+labs(x='Pre-outcome confidence: reference support (0–100)',y=metric,title=title)+theme)
  if(kind=='legacy_confidence')return(ggplot(d,aes(x=confidence_bin,y=.data[[metric]],group=.data[[method]],colour=.data[[method]]))+geom_point()+labs(x='Legacy confidence bin',y=metric,title=title)+theme)
  if(kind=='legacy_coverage')return(ggplot(d,aes(x=actual_coverage_among_eligible,y=.data[[metric]],colour=.data[[method]],group=.data[[method]]))+geom_line()+geom_point()+labs(x='Legacy coverage among eligible subjects',y=metric,title=title)+theme)
  if(kind=='metrics'&&'subset'%in%names(d)){
    if(all(grepl('^support_quality:',d$subset)))return(ggplot(d,aes(actual_coverage,.data[[metric]],colour=.data[[method]],group=.data[[method]]))+geom_point()+geom_line()+labs(x='Actual coverage of original test roster',y=metric,title=title)+theme)
    shapes<-setNames(rep(c(16,17,15,3,7,8),length.out=length(unique(d$subset))),sort(unique(d$subset)))
    p<-ggplot(d,aes(x=.data[[method]],y=.data[[metric]],shape=subset))+geom_point(size=3)+scale_shape_manual(values=shapes)+coord_flip()+labs(x=NULL,y=metric,title=title)+theme
    return(p)
  }
  if(kind=='legacy_performance'){
    p<-ggplot(d,aes(x=.data[[method]],y=.data[[metric]]))+geom_point(size=3)+coord_flip()
    if(metric=='estimate'&&all(c('lower95','upper95')%in%names(d)))p<-p+geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.15)
    recorded<-if('metric'%in%names(d))paste(unique(d$metric),collapse='; ')else metric
    return(p+labs(x=NULL,y=if(metric=='estimate')recorded else metric,title=title)+theme)
  }
  x<-first_col(d,c('feature','role','trait','name','reason','stratum','method','model','target','parameter'))
  if(is.null(x)){d$.row<-seq_len(nrow(d));x<-'.row'}
  if(nrow(d)>1000)return(empty_plot('More than 1,000 rows: refine the table filters before plotting.'))
  ggplot(d,aes(x=.data[[x]],y=.data[[metric]]))+geom_point()+labs(x=x,y=metric,title=title)+theme
}

# An explicit all-ancestry/all-metric choice produces facets, not pooled effects.
make_scoped_plot<-function(d,kind,metric,title='GRID'){
  p<-make_plot(d,kind,metric,title)
  fc<-intersect(c('ancestry','target','metric'),names(d))
  fc<-fc[vapply(fc,function(k)length(unique(d[[k]][!is.na(d[[k]])]))>1L,logical(1))]
  if(length(fc)&&nrow(d)){
    d$.scope<-apply(d[,fc,drop=FALSE],1,paste,collapse=' | ')
    p<-make_plot(d,kind,metric,title)
    if('.scope'%in%names(p$data))p<-p+facet_wrap(~.scope,scales='free')
  }
  p
}
