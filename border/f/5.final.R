source(file.path(Sys.getenv("BORDER_FDIR"), "0.common.R"))
figure_note <- Sys.getenv("BORDER_FIGURE_NOTE", "")
suppressPackageStartupMessages({library(ggplot2); library(ggrepel)})

for(name in c("Fig1_main_associations.png","Fig2_main_associations.png","Fig3_mediation.png")) {
  path <- file.path(figdir,name); if(file.exists(path)) unlink(path)
}

ar <- if(file.exists(file.path(rawdir,"main_associations.rds"))) readRDS(file.path(rawdir,"main_associations.rds")) else tibble()
if(nrow(ar)) ar <- ar %>% filter(is.finite(HR), is.finite(conf.low), is.finite(conf.high), conf.low > 0)
if(nrow(ar)) {
  ar$label <- sprintf("%.2f (%.2f–%.2f)",ar$HR,ar$conf.low,ar$conf.high)
  p <- ggplot(ar,aes(HR,interaction(outcome,exposure),xmin=conf.low,xmax=conf.high))+
    geom_vline(xintercept=1,linetype=2)+geom_errorbar(width=.18, orientation="y")+geom_point(size=2.5)+
    scale_x_log10()+labs(x="Hazard ratio",y=NULL,title="Birth order and incident disease risk",subtitle="Adjusted for age, sex, education, Townsend deprivation index and assessment centre")+
    theme_bw(base_size=11)
  p <- p + labs(caption=figure_note)
  ggsave(file.path(figdir,"Fig2_main_associations.png"),p,width=9,height=5.5,dpi=300)
}

pf <- file.path(rawdir,"phewas.rds")
if(file.exists(pf)) {
  source(file.path(Sys.getenv("BORDER_FDIR"), "phewas_plot.R"))
  plot_phewas_figures(readRDS(pf))
}

mf<-file.path(rawdir,"mediation.rds")
if(file.exists(mf)) {
  m<-as_tibble(readRDS(mf))
  m <- m %>% filter(is.finite(p_boot),is.finite(indirect_lo),is.finite(indirect_hi))
  if(nrow(m)) {
    top<-m%>%filter(is.finite(p_boot))%>%arrange(p_boot)%>%slice_head(n=30)
    p<-ggplot(top,aes(indirect_beta,reorder(mediator,indirect_beta),xmin=indirect_lo,xmax=indirect_hi,color=layer))+
      geom_vline(xintercept=0,linetype=2)+geom_errorbar(width=.15, orientation="y")+geom_point()+facet_grid(outcome~exposure,scales="free_y",space="free_y")+
      labs(x="Indirect coefficient product (a×b)",y=NULL,title="Exploratory omics mediation",subtitle="Associational decomposition; not causal natural indirect effects")+theme_bw(base_size=9)
    p <- p + labs(caption=figure_note)
    ggsave(file.path(figdir,"Fig3_mediation.png"),p,width=13,height=10,dpi=300)
  }
}
cat("Final figures: ",figdir,"\n",sep="")
