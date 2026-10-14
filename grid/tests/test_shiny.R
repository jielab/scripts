#!/usr/bin/env Rscript
# Run explicitly in an R environment. This test is NOT claimed as executed in the
# package-generation environment, which had no Rscript.
args<-commandArgs();f<-grep('^--file=',args,value=TRUE)
ROOT<-normalizePath(file.path(dirname(sub('^--file=','',f[1])),'..'))
suppressPackageStartupMessages({library(shiny);library(DT);library(DBI);library(RSQLite);library(ggplot2);library(jsonlite);library(testthat)})
for(name in c('data.R','plots.R','server.R'))source(file.path(ROOT,'shiny',name),local=FALSE)
DB_PATH<-tempfile(fileext='.sqlite');con<-dbConnect(SQLite(),DB_PATH)
a<-strrep('a',20);b<-strrep('b',20)
d<-data.frame(method=c('ABM_clean','ABM_full','ABM_clean','ABM_full'),ancestry='ALL',subset=c('all','all','selected','selected'),actual_coverage=c(1,1,.5,.5),MSE=c(2,3,1,2))
dbWriteTable(con,paste0('t_',a),d);dbWriteTable(con,paste0('t_',b),d[1,,drop=FALSE])
cat<-data.frame(table_id=c(a,b),sql_table=paste0('t_',c(a,b)),trait='height',module='Twin test',name=c('metrics','one'),rows=c(4,1),source='fixture.tsv',source_sha256='fixture',run_id='test',plot_kind='metrics',status='available',legacy=FALSE,note='synthetic')
dbWriteTable(con,'catalogue',cat);dbWriteTable(con,'metadata',data.frame(synthetic_only=TRUE));dbWriteTable(con,'audit',data.frame(source=character(),status=character(),detail=character()));dbDisconnect(con)
shiny::testServer(server,{
  session$setInputs(trait='height',module='Twin test',table=a,ancestry='ALL',subset='ALL',source_metric='ALL',methods='__ALL__',metric='MSE');session$flushReact()
  stopifnot(nrow(source_values())==4L)
  g<-generation();session$setInputs(table_filter_state=list(generation=g,rows=c(2L,4L),initialized=TRUE));session$flushReact()
  stopifnot(nrow(filtered())==2L,all(filtered()$method=='ABM_full'))
  figure<-current_plot();stopifnot(nrow(figure$data)==2L)
  path<-tempfile();write_filtered(path);export<-read.delim(path);stopifnot(nrow(export)==2L,all(export$method=='ABM_full'))
  # rows_all includes ALL filtered pages; rows_current must not limit the export.
  session$setInputs(table_view_rows_current=2L);session$flushReact();stopifnot(nrow(filtered())==2L)
  session$setInputs(table_filter_state=list(generation=g,rows=integer(),initialized=TRUE));session$flushReact();stopifnot(nrow(filtered())==0L)
  session$setInputs(table=b);session$flushReact();stopifnot(nrow(source_values())==1L,nrow(filtered())==1L)
  # A late event from the previous table cannot filter the new table.
  session$setInputs(table_filter_state=list(generation=g,rows=integer(),initialized=TRUE));session$flushReact();stopifnot(nrow(filtered())==1L)
  session$setInputs(methods=character());session$flushReact();stopifnot(nrow(source_values())==0L)
})
unlink(DB_PATH)
d<-data.frame(method='M',ancestry=rep(c('EUR','AFR'),501),MSE=1)
invisible(ggplot2::ggplot_build(make_scoped_plot(d,'table','MSE')))
d<-data.frame(method='M',subset=paste0('coverage',seq_len(12)),MSE=seq_len(12))
points<-ggplot2::ggplot_build(make_plot(d,'metrics','MSE'))$data[[1]]
stopifnot(nrow(points)==12L,all(!is.na(points$shape)))
cat('Shiny filter/zero-row/generation/download checks passed. JavaScript event delivery still merits a real-browser check.\n')
