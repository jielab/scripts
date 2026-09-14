# Paper_t2dm Fig6f / design slide20: descriptive association ribbons.
# Proxy display selection uses LE8 replication only, never disease P values.
le8_c4_pathway_data<-function(sets,med,n_per_pillar=2L) {
  d<-sets$YS_edges
  needed<-c('feature','component','r_rep','FDR_rep')
  if(!is.data.frame(d)||!all(needed%in%names(d))||!nrow(d))return(tibble())
  d<-d|>filter(is.finite(r_rep),is.finite(FDR_rep),FDR_rep<.05)|>
    arrange(component,desc(abs(r_rep)),feature)|>distinct(component,feature,.keep_all=TRUE)|>
    group_by(component)|>slice_head(n=n_per_pillar)|>ungroup()
  m<-sets$membership|>select(feature,disease_beta,disease_p,disease_FDR)
  if(anyDuplicated(m$feature))stop('Duplicate biomarker membership in pathway plot')
  d<-d|>left_join(m,by='feature')|>mutate(risk=case_when(!is.finite(disease_p)~'Disease association unavailable',
    !is.finite(disease_FDR)|disease_FDR>=.05~'No FDR disease association',
    disease_beta>0~'Higher disease risk',TRUE~'Lower disease risk'),
    left_sign=ifelse(r_rep>0,'Positive association','Inverse association'),
    right_sign=case_when(risk=='Higher disease risk'~'Positive association',risk=='Lower disease risk'~'Inverse association',TRUE~'No supported association'))
  d$mediation_supported<-FALSE
  if(is.data.frame(med)&&all(c('component','feature','FDR_indirect')%in%names(med))) {
    hits<-med|>filter(is.finite(FDR_indirect),FDR_indirect<.05)
    d$mediation_supported<-paste(d$component,d$feature)%in%paste(hits$component,hits$feature)
  }
  d
}

le8_plot_c4_pathways<-function(sets,med,layer,outdir) {
  n<-as.integer(Sys.getenv('C4_FLOW_PER_PILLAR',unset='2'))
  if(is.na(n)||n<1||n>10)stop('C4_FLOW_PER_PILLAR must be 1..10')
  d<-le8_c4_pathway_data(sets,med,n)
  rd<-le8_job_dir(outdir,'c4_connect')
  write_raw_csv(d,'c4.lifestyle_omics_risk_display.csv',rd)
  if(!nrow(d))return(save_plot(blank_plot('LE8–omics–disease paths','No replicated LE8 proxy available'),
    'c4.Fig15.lifestyle_omics_risk.png',18,13,outdir=outdir))
  pillar_names<-c(diet='Diet',pa='Physical activity',smoke='Nicotine exposure',sleep='Sleep',
    bmi='Body mass index',nonhdl='Blood lipids',hba1c='Blood glucose',bp='Blood pressure')
  d<-d|>mutate(pillar=coalesce(unname(pillar_names[sub('[.]pts$','',component)]),component))
  d$.id<-seq_len(nrow(d))
  orders<-list(unique(d$pillar),unique(d$feature),intersect(c('Higher disease risk','Lower disease risk',
    'No FDR disease association','Disease association unavailable'),d$risk))
  cols<-c('pillar','feature','risk');max_height<-nrow(d)+.45*(max(lengths(orders))-1)
  nodes<-list();lanes<-list()
  for(stage in 1:3) {
    cursor<-(max_height-(nrow(d)+.45*(length(orders[[stage]])-1)))/2
    for(name in orders[[stage]]) {
      ids<-d$.id[d[[cols[stage]]]==name];h<-length(ids)
      nodes[[length(nodes)+1L]]<-tibble(stage=stage,name=name,x=stage,ymin=cursor,ymax=cursor+h)
      lanes[[length(lanes)+1L]]<-tibble(.id=ids,stage=stage,lo=cursor+seq_along(ids)-1,hi=cursor+seq_along(ids))
      cursor<-cursor+h+.45
    }
  }
  nodes<-bind_rows(nodes);lanes<-bind_rows(lanes)
  ribbons<-list()
  for(stage in 1:2)for(id in d$.id) {
    a<-lanes|>filter(.id==id,stage==!!stage);b<-lanes|>filter(.id==id,stage==!!(stage+1))
    t<-seq(0,1,length.out=60);s<-3*t^2-2*t^3
    ribbons[[length(ribbons)+1L]]<-tibble(path=paste(stage,id),x=c(stage+.06+t*.88,rev(stage+.06+t*.88)),
      y=c(a$lo+(b$lo-a$lo)*s,rev(a$hi+(b$hi-a$hi)*s)),
      sign=if(stage==1)d$left_sign[id]else d$right_sign[id])
  }
  nodes$fill<-'#C7CDCF'
  for(i in seq_len(nrow(nodes)))if(nodes$stage[i]==1) {
    component<-d$component[match(nodes$name[i],d$pillar)]
    nodes$fill[i]<-if(component%in%names(cols_le8))cols_le8[[component]]else'#5086A1'
  }
  nodes$fill[nodes$name=='Higher disease risk']<-'#DB726D'
  nodes$fill[nodes$name=='Lower disease risk']<-'#6D9ABD'
  nodes$label<-nodes$name
  supported<-unique(d$feature[d$mediation_supported]);nodes$label[nodes$stage==2&nodes$name%in%supported]<-
    paste0(nodes$label[nodes$stage==2&nodes$name%in%supported],' *')
  p<-ggplot(bind_rows(ribbons),aes(x,y,group=path))+
    geom_polygon(aes(fill=sign),alpha=.42,color=NA)+
    scale_fill_manual(values=c('Positive association'='#E99194','Inverse association'='#89B8D8','No supported association'='#D5DADB'),drop=FALSE)+
    geom_rect(data=nodes,aes(xmin=x-.06,xmax=x+.06,ymin=ymin,ymax=ymax),inherit.aes=FALSE,
      fill=nodes$fill,color='grey55',linewidth=.4)+
    geom_label(data=nodes|>filter(stage==2),aes(x=x+.09,y=(ymin+ymax)/2,label=label),inherit.aes=FALSE,
      hjust=0,size=3,linewidth=0,fill='white',alpha=.93)+
    geom_text(data=nodes|>filter(stage==1),aes(x=x-.10,y=(ymin+ymax)/2,label=label),inherit.aes=FALSE,hjust=1,size=3.5)+
    geom_text(data=nodes|>filter(stage==3),aes(x=x+.10,y=(ymin+ymax)/2,label=label),inherit.aes=FALSE,hjust=0,size=3.5)+
    scale_y_reverse(breaks=NULL)+scale_x_continuous(breaks=NULL)+coord_cartesian(xlim=c(.45,3.65),clip='off')+theme_void()+
    labs(title=paste('LE8-supervised',if(layer=='protein')'protein'else'metabolite','paths to',Y),
      subtitle=paste0('Up to ',n,' strongest replicated proxies per pillar; display selection does not use disease results'),fill=NULL,x=NULL,y=NULL,
      caption=paste('LE8 score → omics → disease association. Higher LE8 scores indicate healthier status.',
        'Ribbon colour is the sign of each association; width counts displayed links, not effect size or proportion mediated.',
        '* At least one displayed component–biomarker path has mediation FDR < 0.05. These are associational paths, not proven causal mediation.',
        'All eligible proxies remain in the module tables.'))+
    theme(plot.title=element_text(face='bold',size=15),plot.subtitle=element_text(size=11),legend.position='bottom',
      plot.caption=element_text(size=9,hjust=0),plot.margin=margin(20,40,20,40),
      axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
  save_plot(p,'c4.Fig15.lifestyle_omics_risk.png',18,max(10,min(20,n_distinct(d$feature)*.55)),outdir=outdir)
  invisible(d)
}
