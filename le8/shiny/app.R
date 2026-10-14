#!/usr/bin/env Rscript
# Read-only aggregate explorer. Never loads individual/model RDS or joblib files.
suppressPackageStartupMessages({library(shiny);library(bslib);library(DT);library(DBI);library(RSQLite);library(ggplot2);library(plotly)})
args <- commandArgs(trailingOnly=TRUE)
opt <- function(name, default) {i<-match(name,args);if(is.na(i)) return(default);if(i==length(args))stop(name,' requires a value');args[i+1L]}
script <- grep('^--file=',commandArgs(),value=TRUE)
HERE <- if(length(script)) dirname(normalizePath(sub('^--file=','',script[1]))) else normalizePath('.')
bundled_db <- file.path(HERE,'..','evidence','evidence.sqlite')
default_db <- if(file.exists(bundled_db)) bundled_db else file.path(Sys.getenv('LE8_ANALYSIS_ROOT','/mnt/d/analysis/le8'),'evidence','evidence.sqlite')
DB_PATH <- normalizePath(opt('--db',Sys.getenv('LE8_EVIDENCE_DB',default_db)),mustWork=TRUE)
PORT <- as.integer(opt('--port',Sys.getenv('LE8_SHINY_PORT','3839')));HOST<-opt('--host',Sys.getenv('LE8_SHINY_HOST','127.0.0.1'))
if(!is.finite(PORT)||PORT<1024||PORT>65535)stop('Invalid port')
if(!HOST%in%c('127.0.0.1','localhost','::1') && Sys.getenv('LE8_SHINY_ALLOW_NETWORK')!='TRUE')stop('Network binding requires LE8_SHINY_ALLOW_NETWORK=TRUE')
for(f in c('data.R','plots.R','ui.R','server.R'))source(file.path(HERE,f),local=TRUE)
if('--review'%in%args) {
  con<-open_store();catalogue<-read_catalogue(con)
  if(nrow(catalogue))for(i in seq_len(nrow(catalogue)))query_table(con,catalogue[i,,drop=FALSE],limit=1L)
  DBI::dbDisconnect(con);cat('Shiny database and table access checks passed\n');quit(status=0)
}
shiny::runApp(shinyApp(ui,server),port=PORT,host=HOST,launch.browser=FALSE)
