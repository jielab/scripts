# Only the generated SQLite allowlist is queried. No arbitrary source paths.
open_store <- function() DBI::dbConnect(RSQLite::SQLite(),DB_PATH,flags=RSQLite::SQLITE_RO)
private_names <- function(n) {
  x<-tolower(gsub('[^a-zA-Z0-9]','',n))
  allowed<-if(all(c('runid','modelid','primaryid','referenceid','n','unochorizon')%in%x))'referenceid' else character()
  any(x%in%setdiff(c('eid','iid','fid','id','id1','id2','participantid','individualid','patientid','queryid','donorid','referenceid','sampleid','ukbid','birthdate'),allowed))
}
read_catalogue <- function(con) {
  d<-DBI::dbReadTable(con,'catalogue')
  if(nrow(d)&&any(!grepl('^data_[0-9a-f]{20}$',d$sql_table)))stop('Invalid table allowlist')
  d
}
query_table <- function(con,record,feature='ALL',limit=100000L,trait='all',layer='all') {
  stopifnot(nrow(record)==1L,grepl('^data_[0-9a-f]{20}$',record$sql_table))
  table<-as.character(DBI::dbQuoteIdentifier(con,record$sql_table))
  sql<-paste('SELECT * FROM',table);params<-list();clauses<-character();fields<-DBI::dbListFields(con,record$sql_table)
  if(nzchar(record$feature_column)&&!is.null(feature)&&feature!='ALL') {
    col<-as.character(DBI::dbQuoteIdentifier(con,record$feature_column))
    clauses<-c(clauses,paste0('CAST(',col,' AS TEXT) = ?'));params<-c(params,list(feature))
  }
  tc<-intersect(c('Y','trait'),fields)
  if(length(tc)&&!is.null(trait)&&trait!='all'){clauses<-c(clauses,paste0(DBI::dbQuoteIdentifier(con,tc[1])," IN (?, 'all')"));params<-c(params,list(trait))}
  if('layer'%in%fields&&!is.null(layer)&&layer!='all'){clauses<-c(clauses,paste0(DBI::dbQuoteIdentifier(con,'layer')," IN (?, 'all')"));params<-c(params,list(layer))}
  if(length(clauses))sql<-paste(sql,'WHERE',paste(clauses,collapse=' AND '))
  sql<-paste(sql,'LIMIT',as.integer(limit)+1L)
  d<-if(length(params))DBI::dbGetQuery(con,sql,params=params)else DBI::dbGetQuery(con,sql)
  if(private_names(names(d)))stop('Private columns detected; table withheld')
  attr(d,'truncated')<-nrow(d)>limit
  if(nrow(d)>limit)d<-d[seq_len(limit),,drop=FALSE]
  d
}
features_for <- function(con,id) DBI::dbGetQuery(con,'SELECT feature FROM feature_index WHERE table_id=? ORDER BY feature',params=list(id))$feature
empty_plot <- function(text='No estimate for this selection; missing evidence is not zero.') ggplot()+annotate('text',x=0,y=0,label=text,size=4)+theme_void()
numeric_cols <- function(d) names(d)[vapply(d,is.numeric,logical(1))]
