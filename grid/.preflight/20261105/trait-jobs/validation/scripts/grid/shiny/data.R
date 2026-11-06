open_store<-function(path)DBI::dbConnect(RSQLite::SQLite(),path,flags=RSQLite::SQLITE_RO)
read_values<-function(con,record,limit=100000L){
  name<-record$sql_table
  if(length(name)!=1||!grepl('^t_[a-f0-9]{20}$',name))stop('Invalid registered table')
  d<-DBI::dbGetQuery(con,paste0('SELECT * FROM ',DBI::dbQuoteIdentifier(con,name),' LIMIT ?'),params=list(as.integer(limit+1L)))
  truncated<-nrow(d)>limit
  d<-head(d,limit);attr(d,'truncated')<-truncated;d
}
first_col<-function(d,choices){x<-intersect(choices,names(d));if(length(x))x[1]else NULL}
metric_columns<-function(d){
  primary<-intersect(c('MSE','Brier','prediction_R2','AUC','RMSE','MAE','AUPRC','LogLoss','MSE_over_train_variance','estimate','observed_loss','delta','difference'),names(d))
  unique(c(primary,names(d)[vapply(d,is.numeric,logical(1))]))
}
choose_metric<-function(d){x<-metric_columns(d);if(length(x))x[1]else ''}
filtered_rows<-function(d,state,generation){
  if(is.null(state)||!isTRUE(state$initialized)||!identical(as.integer(state$generation),as.integer(generation)))return(d)
  ix<-as.integer(unlist(state$rows,use.names=FALSE));ix<-unique(ix[!is.na(ix)&ix>=1L&ix<=nrow(d)])
  d[ix,,drop=FALSE]
}
empty_plot<-function(message='No finite results in this selection. Missing is not zero.')ggplot()+annotate('text',x=0,y=0,label=message,size=4)+theme_void()

abm_records<-function(catalogue,trait,run){
  if(is.null(trait)||is.null(run))return(catalogue[FALSE,,drop=FALSE])
  catalogue[catalogue$module=='abm' & catalogue$trait==trait & catalogue$run_id==run,,drop=FALSE]
}
abm_read<-function(con,catalogue,trait,run,name){
  r<-abm_records(catalogue,trait,run);r<-r[r$name==name,,drop=FALSE]
  if(!nrow(r))return(data.frame(status='unavailable_missing_table',reason=paste(name,trait,run)))
  if(nrow(r)!=1L)stop('Ambiguous run/table identity')
  d<-read_values(con,r);d$run_id<-rep(run,nrow(d));d$trait<-rep(trait,nrow(d));d
}
abm_plot_data<-function(d,view){
  if(!nrow(d))return(data.frame())
  keys<-switch(view,
    predictability_curve=c('Selected baseline loss'='matched.selected_base_loss','Control baseline loss'='matched.control_base_loss','Selected - control'='matched.predictability_loss_difference'),
    benefit_curve=c('Selected gain'='matched.selected_gain_vs_base','Control gain'='matched.control_gain_vs_base','Gain enrichment'='matched.gain_enrichment'),
    policy_performance=c('Policy gain vs always-base'='policy_gain_vs_base','Policy gain vs always-candidate'='policy_gain_vs_candidate'),
    count_incremental_contrasts=NULL)
  if(view=='count_incremental_contrasts'){
    if(!all(c('mse_reduction','ci_low','ci_high','objective','contrast','status')%in%names(d)))return(data.frame())
    return(data.frame(x=paste(d$objective,d$contrast),effect='Assessment target MSE reduction',
      value=d$mse_reduction,low=d$ci_low,high=d$ci_high,status=d$status))
  }
  if(view=='policy_performance'&&'method'%in%names(d))d<-d[grepl('^ABM\\(',d$method),,drop=FALSE]
  rows<-lapply(names(keys),function(label){
    key<-keys[[label]];v<-paste0(key,'.estimate');if(!v%in%names(d))return(NULL)
    get<-function(s)if(paste0(key,s)%in%names(d))d[[paste0(key,s)]]else rep(NA,nrow(d))
    data.frame(x=if('nominal_coverage'%in%names(d))d$nominal_coverage else rep('Full query policy',nrow(d)),
      effect=label,value=d[[v]],low=get('.ci_low'),high=get('.ci_high'),status=get('.status'))
  });do.call(rbind,rows)
}
