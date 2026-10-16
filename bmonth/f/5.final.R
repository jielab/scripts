source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
suppressPackageStartupMessages({library(ggplot2);library(ggrepel)})

arfile<-file.path(rawdir,"main_associations.rds")
if(file.exists(arfile)){
  ar<-as_tibble(readRDS(arfile))%>%
    filter(status=="ok",is.finite(HR),is.finite(conf.low),is.finite(conf.high),conf.low>0)
  if(nrow(ar)){
    p<-ggplot(ar,aes(HR,interaction(outcome,exposure),xmin=conf.low,xmax=conf.high))+
      geom_vline(xintercept=1,linetype=2)+geom_errorbar(width=.18,orientation="y")+
      geom_point(size=2.4)+scale_x_log10()+
      labs(x="Hazard ratio",y=NULL,title="Birth month and incident disease risk",
           subtitle="Primary model: age, sex, education, adult TDI, assessment centre, birth north/east and birthplace TDI proxy")+
      theme_bw(base_size=11)
    ggsave(file.path(figdir,"Fig2_main_associations.png"),p,width=10,height=6,dpi=300)
  }
}

pf<-file.path(rawdir,"phewas.rds");jf<-file.path(rawdir,"phewas_joint.rds")
if(file.exists(pf)){
  source(file.path(Sys.getenv("BMONTH_FDIR"),"phewas_plot.R"))
  plot_phewas_figures(readRDS(pf),if(file.exists(jf))readRDS(jf) else NULL)
}

df<-file.path(rawdir,"decomposition_scalar.rds")
if(file.exists(df)){
  d<-as_tibble(readRDS(df))
  if(nrow(d))d<-d%>%filter(is.finite(OR),is.finite(p))
  if(nrow(d)){
    top_names<-d%>%group_by(exposure,phenotype)%>%summarise(pmin=min(p),.groups="drop")%>%
      arrange(pmin)%>%slice_head(n=20)%>%transmute(key=paste(exposure,phenotype,sep="|"))
    d<-d%>%mutate(key=paste(exposure,phenotype,sep="|"))%>%filter(key%in%top_names$key)
    p<-ggplot(d,aes(stage,OR,group=key))+geom_hline(yintercept=1,linetype=2)+
      geom_line(alpha=.45)+geom_point()+scale_y_log10()+
      facet_wrap(~exposure,scales="free_y")+
      labs(x=NULL,y="Odds ratio (common D4 complete-case sample)",
           title="Four-stage decomposition of leading PheWAS signals",
           subtitle="D1 core → D2 adult SES → D3 birth geography/deprivation → D4 migration distance")+
      theme_bw(base_size=9)+theme(axis.text.x=element_text(angle=30,hjust=1))
    ggsave(file.path(figdir,"Fig3_decomposition.png"),p,width=13,height=8,dpi=300)
  }
}

tf<-file.path(rawdir,"temperature_cycle_attenuation.rds")
if(file.exists(tf)){
  t<-as_tibble(readRDS(tf))
  if(nrow(t))t<-t%>%filter(is.finite(cycle_p_before),is.finite(cycle_p_after))
  if(nrow(t)){
    top<-t%>%arrange(cycle_p_before)%>%slice_head(n=30)
    long<-bind_rows(
      top%>%transmute(phenotype,mediator,model="Before temperature",p=cycle_p_before),
      top%>%transmute(phenotype,mediator,model="After temperature",p=cycle_p_after))
    p<-ggplot(long,aes(-log10(pmax(p,.Machine$double.xmin)),reorder(phenotype,-log10(pmax(p,.Machine$double.xmin))),shape=model))+
      geom_point(size=2)+facet_wrap(~mediator,scales="free_y")+
      labs(x="-log10(P) for 2-df annual cycle",y=NULL,
           title="Temperature attenuation of birth-month seasonality")+
      theme_bw(base_size=9)
    ggsave(file.path(figdir,"Fig4_temperature_attenuation.png"),p,width=12,height=9,dpi=300)
  }
}

mf<-file.path(rawdir,"mediation.rds")
if(file.exists(mf)){
  m<-as_tibble(readRDS(mf))%>%filter(is.finite(p_boot),is.finite(indirect_lo),is.finite(indirect_hi))
  if(nrow(m)){
    top<-m%>%arrange(p_boot)%>%slice_head(n=30)
    p<-ggplot(top,aes(indirect_beta,reorder(mediator,indirect_beta),
                      xmin=indirect_lo,xmax=indirect_hi,color=layer))+
      geom_vline(xintercept=0,linetype=2)+geom_errorbar(width=.15,orientation="y")+
      geom_point()+facet_grid(outcome~exposure,scales="free_y",space="free_y")+
      labs(x="Indirect coefficient product (a×b)",y=NULL,
           title="Exploratory omics mediation",
           subtitle="Associational decomposition; not causal natural indirect effects")+
      theme_bw(base_size=9)
    ggsave(file.path(figdir,"Fig5_omics_mediation.png"),p,width=14,height=10,dpi=300)
  }
}
cat("Final figures: ",figdir,"\n",sep="")
