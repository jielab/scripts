args<-commandArgs(TRUE);stopifnot(length(args)==1)
script<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1]);root<-dirname(dirname(normalizePath(script)))
load_defs<-function(root){
 e<-new.env(parent=baseenv())
 walk<-function(x){
  if(missing(x)||(!is.call(x)&&!is.expression(x)))return(NULL)
  if(is.call(x)&&identical(x[[1]],as.name('<-'))&&length(x)==3L&&is.symbol(x[[2]])&&is.call(x[[3]])&&identical(x[[3]][[1]],as.name('function'))){eval(x,e);return(NULL)}
  for(v in as.list(x))walk(v)
 }
 for(file in c('0.common.R','c1.correlate.R'))walk(parse(file.path(root,'f',file)))
 e
}
old<-load_defs(args[1]);current<-load_defs(root)
names<-c('read_all','read_prot','read_met','le8_select_phenotypes','le8_rebuild_baseline','filter_analysis_cohort','make_outcome','t2e','add_attained_age_time','make_prevalent_status','cox_scan','cox_scan_delayed_entry','logistic_scan','prevalent_duration_scan','landmark_incident_scan','risk_window_scan','same_sample_attenuation','std_num','run_c1_pgs_scan','map_c1_pgs_columns','c1_association_contract','c1_pgs_signature')
for(nm in names)stopifnot(identical(formals(old[[nm]]),formals(current[[nm]])),identical(body(old[[nm]]),body(current[[nm]])))
cat('PASS ',length(names),' C1 numerical/input contract function bodies unchanged; dispatch functions are deliberately excluded\n',sep='')
