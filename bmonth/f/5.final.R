source(file.path(Sys.getenv('BMONTH_FDIR'),'0.common.R'))
suppressPackageStartupMessages({library(ggplot2);library(ggrepel)})
arfile<-file.path(rawdir,'main_associations.rds')
if(file.exists(arfile)){
  ar<-as_tibble(readRDS(arfile))%>%filter(status=='ok',is.finite(HR),is.finite(conf.low),is.finite(conf.high),conf.low>0)
  if(nrow(ar)){p<-ggplot(ar,aes(HR,interaction(outcome,exposure),xmin=conf.low,xmax=conf.high))+geom_vline(xintercept=1,linetype=2)+geom_errorbar(width=.18,orientation='y')+geom_point(size=2.4)+scale_x_log10()+labs(x='Hazard ratio',y=NULL,title='Birth month/location and incident disease risk',subtitle='Adjusted for age, sex, education, Townsend deprivation index and assessment centre')+theme_bw(base_size=11);ggsave(file.path(figdir,'Fig2_main_associations.png'),p,width=10,height=6,dpi=300)}
}
pf<-file.path(rawdir,'phewas.rds');jf<-file.path(rawdir,'phewas_joint.rds')
if(file.exists(pf)){source(file.path(Sys.getenv('BMONTH_FDIR'),'phewas_plot.R'));plot_phewas_figures(readRDS(pf),if(file.exists(jf))readRDS(jf) else NULL)}
mf<-file.path(rawdir,'mediation.rds')
if(file.exists(mf)){m<-as_tibble(readRDS(mf))%>%filter(is.finite(p_boot),is.finite(indirect_lo),is.finite(indirect_hi));if(nrow(m)){top<-m%>%arrange(p_boot)%>%slice_head(n=30);p<-ggplot(top,aes(indirect_beta,reorder(mediator,indirect_beta),xmin=indirect_lo,xmax=indirect_hi,color=layer))+geom_vline(xintercept=0,linetype=2)+geom_errorbar(width=.15,orientation='y')+geom_point()+facet_grid(outcome~exposure,scales='free_y',space='free_y')+labs(x='Indirect coefficient product (a×b)',y=NULL,title='Exploratory omics mediation',subtitle='Associational decomposition; not causal natural indirect effects')+theme_bw(base_size=9);ggsave(file.path(figdir,'Fig3_mediation.png'),p,width=14,height=10,dpi=300)}}
cat('Final figures: ',figdir,'\n',sep='')
