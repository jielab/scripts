#!/usr/bin/env Rscript
# Run from the le8 directory with RSQLite and the viewer dependencies installed.
suppressPackageStartupMessages({library(shiny);library(DT);library(DBI);library(RSQLite);library(ggplot2);library(plotly);library(testthat)})
for(f in c('data.R','plots.R','server.R'))source(file.path('shiny',f))
DB_PATH<-tempfile(fileext='.sqlite')
con<-dbConnect(SQLite(),DB_PATH)
d<-expand.grid(feature=c('A','B'),model=c('m1','m2'),analysis=c('a1','a2'),metric=c('u','v'),rep=1:3,stringsAsFactors=FALSE)
d$row_key<-seq_len(nrow(d));d$beta<-d$row_key/100;d$se<-.05
d$task<-'known_cad_state'
catalogue<-data.frame(table_id=c('one','two'),sql_table=paste0('data_',c(strrep('1',20),strrep('2',20))),
  trait='cvd_cad',layer='met',module='Molecular deviation',name=c('effects','other effects'),
  rows=nrow(d),status='completed',run_id=c('run1','run2'),source='synthetic',source_sha256='abc',
  note='synthetic',feature_column='feature',plot_kind='effect')
dbWriteTable(con,'catalogue',catalogue)
dbWriteTable(con,'metadata',data.frame(status='synthetic'))
dbWriteTable(con,'audit',data.frame(status='ok'))
dbWriteTable(con,'feature_index',expand.grid(table_id=c('one','two'),feature=c('A','B')))
dbWriteTable(con,catalogue$sql_table[1],d)
d2<-d;d2$row_key<-d2$row_key+1000
dbWriteTable(con,catalogue$sql_table[2],d2)
dbDisconnect(con)

testServer(server, {
  session$setInputs(trait='all',layer='all',module='all',table='one',feature='ALL',nav='explore')
  session$flushReact()
  expect_equal(nrow(values()),48L)
  expect_equal(nrow(panel_values$explore()),48L) # NULL means not initialized.
  old_generation<-generation()
  rows<-which(d$model=='m1')
  session$setInputs(explore_table_rows_all=rows,
    explore_table_filter_state=list(generation=generation(),rows=rows,initialized=TRUE))
  expect_equal(panel_values$explore()$row_key,d$row_key[rows])
  expect_equal(current_plot()$data$row_key,d$row_key[rows])
  expect_equal(nrow(values()),48L) # The filtered subset never feeds its own DT.
  downloaded<-read.delim(output$explore_data,check.names=FALSE)
  expect_equal(downloaded$row_key,d$row_key[rows])
  session$setInputs(explore_table_rows_current=rows[1:3])
  expect_equal(nrow(read.delim(output$explore_data)),24L) # All filtered pages.
  rows2<-which(d$model=='m1'&d$analysis=='a2'&d$metric=='u')
  session$setInputs(explore_table_rows_all=rows2,
    explore_table_filter_state=list(generation=generation(),rows=rows2,initialized=TRUE))
  expect_equal(current_plot()$data$row_key,read.delim(output$explore_data)$row_key)
  expect_equal(nrow(panel_values$explore()),6L)
  session$setInputs(explore_table_rows_all=integer(),
    explore_table_filter_state=list(generation=generation(),rows=integer(),initialized=TRUE))
  expect_equal(nrow(panel_values$explore()),0L)
  expect_equal(nrow(read.delim(output$explore_data)),0L)
  expect_true(is.null(current_plot()$data)||inherits(current_plot()$data,'waiver'))
  session$setInputs(explore_table_rows_all=seq_len(48),
    explore_table_filter_state=list(generation=generation(),rows=seq_len(48),initialized=TRUE))
  expect_equal(nrow(read.delim(output$explore_data)),48L)
  # Feature changes reject stale rows, even when their numbers fit in the new scope.
  session$setInputs(feature='A')
  expect_gt(generation(),old_generation)
  expect_equal(nrow(panel_values$explore()),24L)
  session$setInputs(explore_table_filter_state=list(generation=old_generation,rows=c(1,2),initialized=TRUE))
  expect_equal(nrow(panel_values$explore()),24L)
  session$setInputs(explore_table_filter_state=list(generation=generation(),rows=c(1,3),initialized=TRUE))
  expect_equal(nrow(read.delim(output$explore_data)),2L)
  session$setInputs(table='two')
  expect_equal(nrow(panel_values$explore()),24L)
  expect_true(all(panel_values$explore()$row_key>1000))
  # Each of the three panels has its own synchronized plot/download row set.
  for(id in c('explore','twin','deviation')) {
    args<-setNames(list(list(generation=generation(),rows=c(2,4,6),initialized=TRUE)),paste0(id,'_table_filter_state'))
    do.call(session$setInputs,args)
    session$setInputs(nav=id)
    expect_equal(panel_plot(id)$data$row_key,panel_values[[id]]()$row_key)
    expect_equal(read.delim(output[[paste0(id,'_data')]])$row_key,panel_values[[id]]()$row_key)
  }
})

test_that('all three molecular plots label the declared task; old tasks stay neutral', {
  for(task in c('known_cad_state','objective_plaque','coronary_calcium',NA_character_)) {
    a<-data.frame(feature=c('A','A'),state=c(0,1),mean_deviation=c(0,.6),
      lower=c(0,.5),upper=c(.5,1),N=c(20,15),beta=c(.6,.7),lo=c(.1,.2),hi=c(1,1.1),
      covariates_before=3,covariates_after=2,task=task)
    if(is.na(task))a$task<-NULL
    label<-deviation_task_label(a)
    for(kind in c('deviation_state','deviation_histogram','effect')) {
      p<-make_plot(a,kind,'synthetic')
      expect_equal(if(kind=='effect')p$labels$subtitle else p$labels$colour,label)
      if(is.na(task)||task=='known_cad_state')expect_false(grepl('客观病变|客观冠脉',paste(unlist(p$labels),collapse=' ')))
      expect_s3_class(ggplot_build(p),'ggplot_built')
      expect_s3_class(ggplotly(p),'plotly')
    }
  }
})
test_that('legacy effect display and download metadata stay neutral', {
  con<-dbConnect(SQLite(),DB_PATH)
  legacy<-data.frame(feature='A',beta=.5,effect_measure='adjusted mean Z difference, objective case minus control')
  dbWriteTable(con,catalogue$sql_table[1],legacy,overwrite=TRUE)
  record<-catalogue[1,,drop=FALSE];record$name<-'feature_state_association'
  viewed<-query_table(con,record)
  expect_true(is.na(viewed$task))
  expect_equal(viewed$task_label,'记录状态/结局状态')
  expect_false(grepl('objective',viewed$effect_measure))
  # The historical SQLite source remains byte-for-byte equivalent at table level.
  expect_equal(dbReadTable(con,record$sql_table),legacy)
  dbDisconnect(con)
})
unlink(DB_PATH)
cat('Shiny testServer: filters, pagination, empty result, scope switching, downloads and task plots passed.\n')
