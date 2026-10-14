open_store<-function(path)DBI::dbConnect(RSQLite::SQLite(),path,flags=RSQLite::SQLITE_RO)
private_columns<-function(d){n<-tolower(gsub('[^A-Za-z0-9]','',names(d)));any(n%in%c('eid','iid','fid','id','participantid','individualid','donorid','familyid','group','queryid','targetid'))}
read_values<-function(con,record,limit=100000L){
  name<-record$sql_table
  if(length(name)!=1||!grepl('^t_[a-f0-9]{20}$',name))stop('Invalid registered table')
  d<-DBI::dbGetQuery(con,paste0('SELECT * FROM ',DBI::dbQuoteIdentifier(con,name),' LIMIT ?'),params=list(as.integer(limit+1L)))
  if(private_columns(d))stop('Private columns are not allowed in this viewer')
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
