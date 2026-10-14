#!/usr/bin/env Rscript
# Aggregate-only browser; no participant RDS, pickle or model objects are loaded.
suppressPackageStartupMessages({library(shiny);library(DT);library(DBI);library(RSQLite);library(ggplot2);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
opt<-function(name,default){i<-match(name,args);if(is.na(i))return(default);if(i==length(args))stop(name,' needs a value');args[i+1L]}
filearg<-grep('^--file=',commandArgs(),value=TRUE)
HERE<-if(length(filearg))dirname(normalizePath(sub('^--file=','',filearg[1])))else normalizePath('.')
DB_PATH<-normalizePath(opt('--db',Sys.getenv('GRID_EVIDENCE_DB','/mnt/d/analysis/grid/evidence/evidence.sqlite')),mustWork=TRUE)
PORT<-as.integer(opt('--port','3839'));HOST<-opt('--host','127.0.0.1')
if(!is.finite(PORT)||PORT<1024||PORT>65535)stop('Invalid port')
if(!HOST%in%c('127.0.0.1','localhost','::1')&&Sys.getenv('GRID_SHINY_ALLOW_NETWORK')!='TRUE')stop('Network binding is disabled without explicit permission')
for(f in c('data.R','plots.R','ui.R','server.R'))source(file.path(HERE,f),local=TRUE)
if('--review'%in%args){
  con<-open_store(DB_PATH)
  catalog<-DBI::dbReadTable(con,'catalogue')
  for(i in seq_len(nrow(catalog))){d<-read_values(con,catalog[i,,drop=FALSE],limit=200);p<-make_plot(d,catalog$plot_kind[i],choose_metric(d),catalog$name[i]);ggplot2::ggplot_build(p)}
  DBI::dbDisconnect(con);cat('All registered table reads and bounded plot builds passed. This is not a browser interaction test.\n');quit(status=0)
}
shiny::runApp(shinyApp(ui,server),host=HOST,port=PORT,launch.browser=FALSE)
