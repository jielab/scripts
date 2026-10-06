suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(data.table);library(survival)})
script <- sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1]);root<-dirname(dirname(normalizePath(script)))
walk <- function(x,only=NULL) {
 if(missing(x)||(!is.call(x)&&!is.expression(x))) return(invisible(NULL))
 if(is.call(x)&&identical(x[[1]],as.name('<-'))&&length(x)==3L&&is.symbol(x[[2]])&&is.call(x[[3]])&&identical(x[[3]][[1]],as.name('function'))) {
  if(is.null(only)||as.character(x[[2]])%in%only) eval(x,.GlobalEnv)
  return(invisible(NULL))
 }
 for(z in as.list(x)) walk(z,only)
}
for(f in c('0.common.R','c1.correlate.R')) walk(parse(file.path(root,'f',f)))
if(exists('split',envir=.GlobalEnv,inherits=FALSE))rm(split,envir=.GlobalEnv)
SEED<-2026L;N_CORES<-8L;Y<-'test';LE8_REPLACE<-FALSE
scratch<-tempfile('le8-final-R-',tmpdir='/tmp');dir.create(scratch);Sys.setenv(LE8_ANALYSIS_ROOT=scratch,LE8_RESOURCE_MEMORY_GIB='24',LE8_PHASE_CORES='8')
set.seed(14);n<-1200L;features<-paste0('F',1:64)
base<-data.frame(eid=paste0('p',1:n),age=runif(n,40,70),sex=rbinom(n,1,.5),family=rep(1:(n/2),each=2),
 test.t2e=sample(1:12,n,TRUE),test.Yt2e=rbinom(n,1,.4),prevalent_status=rbinom(n,1,.3),.omic_overlap=seq_len(n)<=900)
base$.attained_entry<-base$age;base$.attained_exit<-base$age+base$test.t2e
scores<-data.frame(eid=base$eid)
for(j in 1:64) {scores[[features[j]]]<-rnorm(n)+.2*base$test.Yt2e;scores[sample(n,20+j),features[j]]<-NA_real_}
dat<-cbind(base,scores[,-1]);score_map<-setNames(features,features);covars<-c('age','sex');outcome<-Y;tvar<-'test.t2e';evar<-'test.Yt2e'
# Evaluate the exact six-model closure embedded in production run_c1_pgs_scan.
walk(body(run_c1_pgs_scan),'fit_feature')
passed<-0L;results<-list();timing<-list()
for(workers in c(1,2,4,8)) {
 Sys.setenv(LE8_PWAS_WORKERS=workers,LE8_PGS_WORKERS=workers)
 start<-proc.time()[['elapsed']]
 z<-list(incident=cox_scan(dat,features,covars,Y),delayed=cox_scan_delayed_entry(dat,features,'sex',Y),prev=logistic_scan(dat,features,covars,'prevalent_status'))
 if(workers==1)reference<-z else stopifnot(isTRUE(all.equal(z,reference,tolerance=1e-12)))
 stopifnot(length(unique(z$incident$N_total))>1L,all(is.finite(z$incident$std.error)))
 # Robust Cox is exercised through the same feature scheduler; production Cox formulas are left unchanged.
 robust<-le8_dynamic_map(as.list(features),function(f) {
  form<-as.formula(paste('Surv(test.t2e,test.Yt2e)~',f,'+age+sex+cluster(family)'))
  coef(summary(coxph(form,dat,ties='efron')))[1,]
 },workers=workers,task_ids=features)
 if(workers==1)robust_ref<-robust else stopifnot(isTRUE(all.equal(robust,robust_ref,tolerance=1e-12)))
 rows<-c1_pgs_checkpoint_scan(score_map,fit_feature,paste0('six-model-',workers),'prot')
 if(workers==1)pgs_reference<-rows else stopifnot(isTRUE(all.equal(rows,pgs_reference,tolerance=1e-12)))
 timing[[length(timing)+1]]<-data.frame(workers=workers,features=64,elapsed_seconds=proc.time()[['elapsed']]-start)
 passed<-passed+3L;cat('PASS workers=',workers,' PWAS/delayed/logistic and robust dispatch and PGS six models\n',sep='')
}
Sys.setenv(LE8_PGS_WORKERS='4');sig<-'six-model-4';cache<-le8_cache_dir('c1_pgs','prot',sig);files<-list.files(cache,full.names=TRUE)
rows<-c1_pgs_checkpoint_scan(score_map,function(f)stop('Unexpected recomputation'),sig,'prot');stopifnot(identical(rows,pgs_reference))
unlink(file.path(cache,paste0(pgs_hash('F3'),'.rds')))
rows<-c1_pgs_checkpoint_scan(score_map,function(f) {stopifnot(f=='F3');fit_feature(f)},sig,'prot');stopifnot(identical(rows,pgs_reference));passed<-passed+1L
for(workers in c(1,2,4,8)) {
 z<-le8_dynamic_map(as.list(features),function(f) {stopifnot(Sys.getenv('OMP_NUM_THREADS')=='1');c(runif(3),le8_dynamic_map(list(1),function(x)Sys.getpid(),workers=8)[[1]]==Sys.getpid())},workers=workers,task_ids=features)
 if(workers==1)random_ref<-z else stopifnot(identical(z,random_ref))
}
z<-le8_dynamic_map(list(1,2),function(x)NULL,workers=2,task_ids=c('a','b'));stopifnot(length(z)==2,all(vapply(z,is.null,logical(1))));passed<-passed+1L
fail<-try(le8_dynamic_map(list(1,2),function(x)stop('intentional'),workers=2),silent=TRUE);stopifnot(inherits(fail,'try-error'));passed<-passed+1L
fail<-try(le8_parallel_workers(4,4,20,4,24,2),silent=TRUE);stopifnot(inherits(fail,'try-error'));passed<-passed+1L
write.csv(bind_rows(timing),file.path(scratch,'timings.csv'),row.names=FALSE)
cat('PASS count=',passed,'; fixture and timings=',scratch,'\n',sep='')
