source(file.path(Sys.getenv('BMONTH_FDIR'),'0.common.R'));source(file.path(Sys.getenv('BMONTH_FDIR'),'phewas_plot.R'))
pf<-file.path(rawdir,'phewas.rds');if(!file.exists(pf))stop('Saved PheWAS missing');jf<-file.path(rawdir,'phewas_joint.rds');plot_phewas_figures(readRDS(pf),if(file.exists(jf))readRDS(jf) else NULL)
