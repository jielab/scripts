source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
suppressPackageStartupMessages({library(ggplot2);library(ggrepel)})

for (name in c("Fig2_main_associations", "FigS_exploratory_associations", "Fig3_decomposition",
  "FigS_exploratory_decomposition", "FigS_geography_alcohol_audit", "Fig4_temperature_attenuation",
  "Fig5_omics_mediation", "FigS_exploratory_omics_mediation")) remove_result_figure(name)

arfile<-file.path(rawdir,"main_associations.rds")
if(file.exists(arfile)){
  associations<-as_tibble(readRDS(arfile))%>%
    filter(status=="ok",is.finite(HR),is.finite(conf.low),is.finite(conf.high),conf.low>0)
  for (exploratory in c(FALSE, TRUE)) {
    ar <- associations[(associations$exposure %in% EXPLORATORY_ALL) == exploratory,,drop=FALSE]
    if(nrow(ar)){
      p<-ggplot(ar,aes(HR,interaction(outcome,exposure),xmin=conf.low,xmax=conf.high))+
        geom_vline(xintercept=1,linetype=2)+geom_errorbar(width=.18,orientation="y")+
        geom_point(size=2.4)+scale_x_log10()+
        labs(x="Hazard ratio",y=NULL,
             title=if(exploratory) "Exploratory birth-order traits and incident disease risk" else "Birth month and incident disease risk",
             subtitle=if(exploratory) "D2: age, sex, centre, education and adult TDI; trait-specific eligibility" else
               "Primary: age, sex, education, adult TDI, centre, birth north/east and birthplace TDI proxy")+
        theme_bw(base_size=11)
      name <- if(exploratory) "FigS_exploratory_associations" else "Fig2_main_associations"
      save_result_figure(p,name,ar,width=10,height=6)
    }
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
  for (exploratory in c(FALSE, TRUE)) {
    subset <- d[(d$exposure %in% EXPLORATORY_ALL) == exploratory,,drop=FALSE]
    if(nrow(subset)){
      top_names<-subset%>%group_by(exposure,phenotype)%>%summarise(pmin=min(p),.groups="drop")%>%
        arrange(pmin)%>%slice_head(n=30)%>%transmute(key=paste(exposure,phenotype,sep="|"))
      d2<-subset%>%mutate(key=paste(exposure,phenotype,sep="|"))%>%filter(key%in%top_names$key)
      p<-ggplot(d2,aes(stage,OR,group=key))+geom_hline(yintercept=1,linetype=2)+
        geom_line(alpha=.45)+geom_point()+scale_y_log10()+
        facet_wrap(~exposure,scales="free_y")+
        labs(x=NULL,y="Odds ratio (common D5 complete-case sample)",
             title=if(exploratory) "Exploratory birth-order signals: five-stage decomposition" else
               "Five-stage decomposition of leading PheWAS signals",
             subtitle="D1 core → D2 adult SES → D3 birth geography → D4 birthplace TDI → D5 migration")+
        theme_bw(base_size=9)+theme(axis.text.x=element_text(angle=30,hjust=1))
      name <- if(exploratory) "FigS_exploratory_decomposition" else "Fig3_decomposition"
      save_result_figure(p,name,d2,width=13,height=8)
    }
  }

  # Dedicated audit of the OLD geography-alcohol finding.
  a<-d%>%filter(exposure%in%c("birth_north_100km","birth_east_100km"),
                phenotype%in%c("317","317.1","317.11"))
  if(nrow(a)){
    a$label <- factor(a$phenotype,levels=c("317","317.1","317.11"),
                      labels=c("Alcohol-related disorders","Alcoholism","Alcoholic liver damage"))
    p<-ggplot(a,aes(stage,OR,group=interaction(exposure,label),shape=exposure))+
      geom_hline(yintercept=1,linetype=2)+geom_line()+
      geom_point(size=2.2)+facet_wrap(~label,scales="free_y")+
      scale_y_log10()+
      labs(x=NULL,y="Odds ratio per 100 km",
           title="Why did the OLD alcohol signal disappear?",
           subtitle="Same participants at every stage: D2 reproduces OLD-style geography PheWAS; D3/D4 isolate geography and birthplace-TDI attenuation")+
      theme_bw(base_size=9)+theme(axis.text.x=element_text(angle=30,hjust=1))
    save_result_figure(p,"FigS_geography_alcohol_audit",a,width=12,height=6.5,dpi=320)
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
    p<-ggplot(long,aes(-log10(pmax(p,.Machine$double.xmin)),
                       reorder(phenotype,-log10(pmax(p,.Machine$double.xmin))),shape=model))+
      geom_point(size=2)+facet_wrap(~mediator,scales="free_y")+
      labs(x="-log10(P) for 2-df annual cycle",y=NULL,
           title="Temperature attenuation of birth-month seasonality")+
      theme_bw(base_size=9)
    save_result_figure(p,"Fig4_temperature_attenuation",long,width=12,height=9)
  }
}

mf<-file.path(rawdir,"mediation.rds")
if(file.exists(mf)){
  mediation<-as_tibble(readRDS(mf))%>%filter(is.finite(p_boot),is.finite(indirect_lo),is.finite(indirect_hi))
  for (exploratory in c(FALSE, TRUE)) {
    m <- mediation[(mediation$exposure %in% EXPLORATORY_ALL) == exploratory,,drop=FALSE]
    if(nrow(m)){
      top<-m%>%arrange(p_boot)%>%slice_head(n=30)
      p<-ggplot(top,aes(indirect_beta,reorder(mediator,indirect_beta),
                        xmin=indirect_lo,xmax=indirect_hi,color=layer))+
        geom_vline(xintercept=0,linetype=2)+geom_errorbar(width=.15,orientation="y")+
        geom_point()+facet_grid(outcome~exposure,scales="free_y",space="free_y")+
        labs(x="Indirect coefficient product (a×b)",y=NULL,
             title=if(exploratory) "Exploratory birth-order traits: omics mediation" else "Exploratory omics mediation",
             subtitle="Associational decomposition; not causal natural indirect effects")+
        theme_bw(base_size=9)
      name <- if(exploratory) "FigS_exploratory_omics_mediation" else "Fig5_omics_mediation"
      save_result_figure(p,name,top,width=14,height=10)
    }
  }
}
cat("Final figures: ",figdir,"\n",sep="")
