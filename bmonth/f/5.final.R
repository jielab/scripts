source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
suppressPackageStartupMessages({library(ggplot2);library(ggrepel)})

for (name in c("Fig2_main_associations","FigS_exploratory_associations","Fig3_decomposition",
  "FigS_context_decomposition","FigS_exploratory_decomposition","FigS_geography_alcohol_audit",
  "FigS_birth_order_family_size_sensitivity","Fig4_temperature_attenuation",
  "Fig5_omics_mediation","FigS_exploratory_omics_mediation")) remove_result_figure(name)

arfile<-file.path(rawdir,"main_associations.rds")
if(length(OUTCOMES) && file.exists(arfile)){
  associations<-read_survival_result(arfile,c("exposure","status","HR","conf.low","conf.high"),empty_associations())%>%
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

# Family-size sensitivity figure for birth-order Cox analyses.
sf<-file.path(rawdir,"birth_order_family_size_sensitivity.rds")
if(length(OUTCOMES) && file.exists(sf) && file.exists(arfile)){
  prim<-read_survival_result(arfile,c("exposure","status","HR","conf.low","conf.high"),empty_associations())%>%
    filter(exposure%in%FAMILY_SIZE_SENS_TRAITS,status=="ok")%>%
    mutate(model="D2 primary")
  fam<-read_survival_result(sf,c("exposure","status","HR","conf.low","conf.high"),empty_associations())%>%
    filter(exposure%in%FAMILY_SIZE_SENS_TRAITS,status=="ok")%>%
    mutate(model="D2 + categorical sibship size")
  z<-bind_rows(prim,fam)%>%
    filter(is.finite(HR),is.finite(conf.low),is.finite(conf.high),conf.low>0)
  if(nrow(z)){
    p<-ggplot(z,aes(HR,interaction(outcome,exposure),xmin=conf.low,xmax=conf.high,shape=model))+
      geom_vline(xintercept=1,linetype=2)+
      geom_errorbar(width=.16,orientation="y",position=position_dodge(width=.4))+
      geom_point(size=2.3,position=position_dodge(width=.4))+scale_x_log10()+
      labs(x="Hazard ratio",y=NULL,title="Birth-order sensitivity to family size",
           subtitle="Sibship size entered categorically as 2, 3, 4, 5, or 6+; border12 omitted because family size is fixed")+
      theme_bw(base_size=10)
    save_result_figure(p,"FigS_birth_order_family_size_sensitivity",z,width=11,height=6.5)
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
  if("sample_set"%in%names(d))d<-d%>%filter(sample_set=="common_D5")
  if(nrow(d))d<-d%>%filter(is.finite(OR),is.finite(p))

  tracks<-list(
    birth_month=list(
      data=d%>%filter(analysis_track=="birth_month"),
      name="Fig3_decomposition",
      title="Birth-month signals: five-stage decomposition"),
    birthplace_context=list(
      data=d%>%filter(analysis_track=="birthplace_context"),
      name="FigS_context_decomposition",
      title="Birthplace-context signals: five-stage decomposition"),
    exploratory_birth_order=list(
      data=d%>%filter(analysis_track=="exploratory_birth_order"),
      name="FigS_exploratory_decomposition",
      title="Exploratory birth-order signals: five-stage decomposition")
  )

  for(obj in tracks){
    subset<-obj$data
    if(!nrow(subset))next
    top_names<-subset%>%group_by(exposure,phenotype)%>%
      summarise(pmin=min(p),.groups="drop")%>%arrange(pmin)%>%
      slice_head(n=30)%>%transmute(key=paste(exposure,phenotype,sep="|"))
    d2<-subset%>%mutate(key=paste(exposure,phenotype,sep="|"))%>%
      filter(key%in%top_names$key)
    p<-ggplot(d2,aes(stage,OR,group=key))+geom_hline(yintercept=1,linetype=2)+
      geom_line(alpha=.45)+geom_point()+scale_y_log10()+
      facet_wrap(~exposure,scales="free_y")+
      labs(x=NULL,y="Odds ratio (common D5 complete-case sample)",
           title=obj$title,
           subtitle="D1 core → D2 adult SES → D3 birth geography → D4 birthplace TDI → D5 migration")+
      theme_bw(base_size=9)+theme(axis.text.x=element_text(angle=30,hjust=1))
    save_result_figure(p,obj$name,d2,width=13,height=8)
  }

  a<-d%>%filter(exposure%in%c("birth_north_100km","birth_east_100km"),
                phenotype%in%c("317","317.1","317.11"))
  if(nrow(a)){
    a$label<-factor(a$phenotype,levels=c("317","317.1","317.11"),
      labels=c("Alcohol-related disorders","Alcoholism","Alcoholic liver damage"))
    p<-ggplot(a,aes(stage,OR,group=interaction(exposure,label),shape=exposure))+
      geom_hline(yintercept=1,linetype=2)+geom_line()+geom_point(size=2.2)+
      facet_wrap(~label,scales="free_y")+scale_y_log10()+
      labs(x=NULL,y="Odds ratio per 100 km",
           title="Attenuation of the birthplace-geography alcohol signal",
           subtitle="D2 uses the OLD covariate specification on the common D5 sample; full-sample OLD-style replication is the separate context PheWAS")+
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
           title="Annual-cycle evidence with temperature adjustment (exploratory)")+
      theme_bw(base_size=9)
    save_result_figure(p,"Fig4_temperature_attenuation",long,width=12,height=9)
  }
}

mf<-file.path(rawdir,"mediation.rds")
if(length(OUTCOMES) && file.exists(mf)){
  source(file.path(Sys.getenv("BMONTH_FDIR"),"mediation.R"))
  mediation<-read_survival_result(mf,c("exposure","mediator","layer","p_boot","indirect_beta","indirect_lo","indirect_hi"),empty_mediation())
  if(nrow(mediation))mediation<-mediation%>%filter(is.finite(p_boot),is.finite(indirect_lo),is.finite(indirect_hi))
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
# Exposure-focused displays from the locked follow-up, including unsuccessful
# rows in the companion data rather than presenting them as negative findings.
for(name in c("Fig_birth_month_cycle","Fig_birth_order_family_size"))remove_result_figure(name)
curve_file <- file.path(rawdir,"followup_cycle_monthly_curves.csv")
if(length(OUTCOMES) && file.exists(curve_file)) {
  curves <- read_survival_result(curve_file,c("analysis","month","HR","conf.low","conf.high"))
  if(nrow(curves)) {
    show <- curves[curves$analysis=="D4_joint_cycle",,drop=FALSE]
    if(nrow(show)) {
      p <- ggplot(show,aes(month,HR))+
        geom_hline(yintercept=1,linetype=2,color="grey60")+
        geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#477caa",alpha=.18)+
        geom_line(color="#245b8a",linewidth=.8)+geom_point(color="#245b8a",size=1.4)+
        facet_wrap(~outcome)+scale_x_continuous(breaks=1:12,labels=month.abb)+
        labs(x="Birth month",y="Hazard ratio relative to January",
          title="Birth-month annual cycle",
          subtitle="Joint negative-cosine + sine model; pointwise 95% confidence intervals")+
        theme_bw(base_size=11)+theme(axis.text.x=element_text(angle=45,hjust=1))
      save_result_figure(p,"Fig_birth_month_cycle",curves,width=12,height=5)
      joint_file <- file.path(rawdir,"followup_cycle_cox.csv")
      if(file.exists(joint_file))
        writexl::write_xlsx(list(monthly_curves=curves,
          joint_tests=read.csv(joint_file,check.names=FALSE)),
          file.path(figdir,"Fig_birth_month_cycle.xlsx"))
    }
  }
}
target_file <- file.path(rawdir,"followup_targeted.csv")
if(file.exists(target_file)) {
  target <- read.csv(target_file,check.names=FALSE,colClasses=c(phenotype="character"))
  if(nrow(target)) {
    paired <- target[target$analysis%in%c("full_D2","full_D2_family_size"),,drop=FALSE]
    paired$model <- ifelse(paired$analysis=="full_D2","D2","D2 + family size")
    show <- paired[grepl("^ok_",paired$status)&is.finite(paired$OR)&
      is.finite(paired$conf.low)&is.finite(paired$conf.high),,drop=FALSE]
    if(nrow(show)) {
      p <- ggplot(show,aes(OR,phenotype,xmin=conf.low,xmax=conf.high,color=model))+
        geom_vline(xintercept=1,linetype=2,color="grey60")+
        geom_errorbar(width=.15,orientation="y",position=position_dodge(width=.5))+
        geom_point(position=position_dodge(width=.5))+
        facet_wrap(~exposure,scales="free_x")+scale_x_log10()+
        labs(x="Odds ratio",y="Locked exploratory phecode",color=NULL,
          title="Birth order and family-size adjustment",
          subtitle="Same participants within each comparison; border12 already fixes family size at two")+
        theme_bw(base_size=10)+theme(legend.position="bottom")
      save_result_figure(p,"Fig_birth_order_family_size",paired,width=12,height=8)
    }
  }
}
cat("Final figures: ",figdir,"\n",sep="")
