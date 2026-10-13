# Rendering only. Use the shared helper so all callers use the same symbols.
suppressPackageStartupMessages(source(file.path(dir0,'scripts','0f','plotting.R')))
plot_phewas_figures <- function(res, sample_info=NULL) {
  res<-as_tibble(res)
  obsolete<-c('Fig2_phewas.png',paste0('Fig1_phewas_',ALL_EXPOSURES,'.png'))
  if(is.null(sample_info)) {
    cohort<-readRDS(file.path(rawdir,'border.cohort.rds'))
    phenotype_ids<-readRDS(file.path(indir,'Rdata','PheWAS.rds'))$id
    types<-c(border12='binary',border1n='binary',border_rank='ordinal',border_relative='continuous')
    sample_info<-setNames(lapply(unique(res$exposure),function(x)
      phewas_sample_info(cohort,x,covars,phenotype_ids,type=types[[x]])),unique(res$exposure))
  }
  xmap<-make_phewas_xmap(res)
  summaries<-samples<-list()
  for(name in obsolete){path<-file.path(figdir,name);if(file.exists(path))unlink(path)}
  for(x in unique(res$exposure)) {
    d<-res[res$exposure==x,,drop=FALSE];si<-sample_info[[x]]
    p<-plot_phewas_results(d,title=paste0(x,' PheWAS'),sample_info=si,xmap=xmap)
    if(is.null(p)){message('[PheWAS plot] No estimable tests for ',x,'; no figure generated');next}
    d<-d[is.finite(d$p)&d$p>=0&d$p<=1,,drop=FALSE];threshold<-unique(d$bonferroni_threshold)
    note<-Sys.getenv('BORDER_FIGURE_NOTE','')
    # Only isolated synthetic tests set this marker; production has no footer.
    if(nzchar(note))p<-p+ggplot2::labs(subtitle=note)
    path<-file.path(figdir,paste0('Fig1_phewas_',x,'.png'))
    ggplot2::ggsave(path,p,width=13,height=6.5,dpi=320,bg='white')
    message('[PheWAS plot] Saved ',path,' (',si$label,')')
    summaries[[x]]<-d %>% mutate(category=tolower(trimws(category))) %>% group_by(category) %>%
      summarise(n_tested=n(),n_bonf=sum(p<threshold),.groups='drop') %>% mutate(exposure=x,.before=1)
    samples[[x]]<-tibble(exposure=x,type=si$type,N=si$N,cases=si$cases,controls=si$controls,title_counts=si$label)
  }
  empty<-tibble(exposure=character(),category=character(),n_tested=integer(),n_bonf=integer())
  write_csv2(if(length(summaries))bind_rows(summaries) else empty,'phewas_category_summary.csv')
  if(length(samples))write_csv2(bind_rows(samples),'phewas_plot_samples.csv')
  invisible(NULL)
}
