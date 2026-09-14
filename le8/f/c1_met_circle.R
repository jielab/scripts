# Metabolite associations by biochemical class, without genomic coordinates.
# Inspired by the category-wise radial bars in Li et al., Nature Medicine 2026,
# Fig2 (s41591-025-04105-8); all plotted observations come from the current run.
le8_met_circle_data<-function(d) {
  d<-d|>filter(is.finite(p.value),p.value>=0,p.value<=1)
  # Keep the derived columns even when every association is unavailable.
  if(!nrow(d))return(d|>mutate(neglog10_FDR=numeric(),FDR=numeric(),
    circle_group=character(),significant=logical()))
  logp<-log(pmax(d$p.value,.Machine$double.xmin))
  if('statistic'%in%names(d)) {
    recover<-d$p.value==0 & is.finite(d$statistic)
    logp[recover]<-log(2)+pnorm(abs(d$statistic[recover]),lower.tail=FALSE,log.p=TRUE)
  }
  # BH family includes all tested metabolites, before selecting significant bars.
  ix<-order(logp);m<-length(logp)
  logq<-numeric(m);logq[ix]<-pmin(0,rev(cummin(rev(logp[ix]+log(m/seq_len(m))))))
  d$neglog10_FDR<--logq/log(10);d$FDR<-exp(logq)
  if(!'super_group'%in%names(d))d$super_group<-met_super_group(d$group)
  d|>mutate(circle_group=coalesce(super_group,'Other metabolites'),
    significant=is.finite(beta)&neglog10_FDR> -log10(.05))
}

plot_met_circle<-function(d,title,label_n=Inf,fdr_cap=50,beta_limit=.3,show_legend=TRUE) {
  all<-le8_met_circle_data(d)
  if(!nrow(all))return(blank_plot(title,'No valid metabolite P values available'))
  d<-all|>filter(significant)
  if(!nrow(d))return(blank_plot(title,'No metabolites with BH FDR < 0.05'))
  stopifnot(is.finite(fdr_cap),fdr_cap>0,is.finite(beta_limit),beta_limit>0)
  group_order<-c('Lipoprotein lipids','Cholesterol and apolipoproteins','Fatty acids',
    'Amino acids','Energy metabolism','Other metabolites')
  group_names<-c('Lipoprotein lipids'='Lipoprotein lipids',
    'Cholesterol and apolipoproteins'='Cholesterol / Apo','Fatty acids'='Fatty acids',
    'Amino acids'='Amino acids','Energy metabolism'='Energy','Other metabolites'='Other')
  d<-d|>arrange(match(circle_group,group_order),term)
  # Blank angular slots separate categories and leave room for a radial scale.
  parts<-split(d,factor(d$circle_group,levels=unique(d$circle_group)))
  offset<-6L
  for(i in seq_along(parts)) {
    parts[[i]]$id<-offset+seq_len(nrow(parts[[i]]));offset<-max(parts[[i]]$id)+3L
  }
  z<-bind_rows(parts);slots<-offset+3L
  z<-z|>mutate(angle0=90-360*(id-.5)/slots,hjust=ifelse(angle0< -90,1,0),
    angle=ifelse(angle0< -90,angle0+180,angle0),height=pmin(neglog10_FDR,fdr_cap),
    lab=ifelse(rank(-neglog10_FDR,ties.method='first')<=label_n,sub('^met_','',term),NA_character_))
  bands<-z|>group_by(circle_group)|>summarise(x1=min(id)-.45,x2=max(id)+.45,.groups='drop')|>
    mutate(x=(x1+x2)/2,angle=(180-360*(x-.5)/slots+90)%%180-90,
      label=unname(group_names[circle_group]))
  ticks<-data.frame(y=c(0,fdr_cap/2,fdr_cap),label=c('0',format(fdr_cap/2),format(fdr_cap)))
  p<-ggplot(z,aes(id,height))+
    geom_col(aes(fill=beta),width=.88,na.rm=TRUE)+
    geom_segment(data=bands,aes(x=x1,xend=x2,y=-fdr_cap*.02,yend=-fdr_cap*.02),
      inherit.aes=FALSE,linewidth=.5,color='grey35')+
    geom_text(data=bands,aes(x=x,y=-fdr_cap*.15,label=label,angle=angle),
      inherit.aes=FALSE,size=2.3,color='grey20')+
    geom_text(aes(y=height+fdr_cap*.025,label=lab,angle=angle,hjust=hjust),size=1.85,
      color='grey35',na.rm=TRUE)+
    annotate('segment',x=3,xend=3,y=0,yend=fdr_cap,linewidth=.3,color='grey60')+
    geom_segment(data=ticks,aes(x=2,xend=4,y=y,yend=y),inherit.aes=FALSE,linewidth=.3,color='grey55')+
    geom_text(data=ticks,aes(x=1,y=y,label=label),inherit.aes=FALSE,hjust=1,size=2.2,color='grey40')+
    annotate('text',x=3,y=fdr_cap*1.10,label='-log10(FDR)',size=2.8,hjust=.5)+
    scale_x_continuous(limits=c(.5,slots+.5),expand=c(0,0),breaks=NULL)+
    scale_y_continuous(limits=c(-fdr_cap*.5,fdr_cap*1.35),expand=c(0,0),breaks=NULL)+
    scale_fill_gradient2(low='#0878BF',mid='white',high='#FF5A52',midpoint=0,
      limits=c(-beta_limit,beta_limit),oob=scales::squish,name='Log effect per SD',
      breaks=c(-beta_limit,0,beta_limit))+
    coord_polar(clip='off')+labs(title=title,x=NULL,y=NULL,
      subtitle=paste0(nrow(d),' / ',nrow(all),' metabolites: FDR < 0.05; ',sum(z$neglog10_FDR>fdr_cap),
        ' bars capped at ',fdr_cap))+theme_void(11)+
    theme(plot.title=element_text(face='bold',hjust=.5),legend.position=if(show_legend)'bottom'else'none',
      axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),
      plot.margin=margin(12,15,12,15))
  attr(p,'le8_circle_data')<-all
  p
}
