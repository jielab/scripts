suppressPackageStartupMessages({
  library(shiny); library(bslib); library(plotly); library(DT); library(DBI); library(RSQLite); library(data.table); library(scales)
})

`%||%` <- function(x,y) if (is.null(x) || length(x)==0 || is.na(x) || !nzchar(x)) y else x
ofile <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
app_dir <- normalizePath(dirname(ofile %||% getwd()), mustWork = FALSE)

default_db <- normalizePath(file.path(Sys.getenv("GU_ANALYSIS_ROOT", "/mnt/d/analysis/gu"), "Rshiny", "gu.sqlite"), mustWork = FALSE)
db_path <- Sys.getenv("GU_SQLITE", unset = default_db)
if (!file.exists(db_path)) {
  stop("GU SQLite database not found: ", db_path, "\nRun: ./gu.sh normalize")
}

chr_order <- c(as.character(1:22), "X", "Y", "MT")
fmt_bp <- function(x) ifelse(x >= 1e6, paste0(round(x/1e6,2)," Mb"), ifelse(x>=1e3,paste0(round(x/1e3,1)," kb"),paste0(x," bp")))
q <- function(con, sql, params = NULL) {
  if (is.null(params) || length(params) == 0L) DBI::dbGetQuery(con, sql)
  else DBI::dbGetQuery(con, sql, params = params)
}
scalar_q <- function(con, sql, params=NULL) { z <- q(con,sql,params); if(nrow(z)) z[[1]][1] else NA }

reduce_intervals <- function(d) {
  if (!nrow(d)) return(d[, .(start,end)])
  x <- as.data.table(d)[order(start,end), .(start=as.numeric(start),end=as.numeric(end))]
  out <- vector("list", nrow(x)); k <- 0L; s <- x$start[1]; e <- x$end[1]
  if (nrow(x)>1) for(i in 2:nrow(x)) {
    if (x$start[i] <= e) e <- max(e,x$end[i]) else { k<-k+1L; out[[k]]<-c(s,e); s<-x$start[i]; e<-x$end[i] }
  }
  k<-k+1L; out[[k]]<-c(s,e)
  m <- do.call(rbind,out[seq_len(k)])
  data.table(start=as.numeric(m[,1]), end=as.numeric(m[,2]))
}

interval_bp <- function(d) { r <- reduce_intervals(as.data.table(d)); if(!nrow(r)) 0 else sum(r$end-r$start) }
intersection_bp <- function(a,b) {
  a <- reduce_intervals(as.data.table(a)); b <- reduce_intervals(as.data.table(b)); i<-1L;j<-1L;z<-0
  while(i<=nrow(a) && j<=nrow(b)) { z<-z+max(0,min(a$end[i],b$end[j])-max(a$start[i],b$start[j])); if(a$end[i]<b$end[j]) i<-i+1L else j<-j+1L }
  z
}
